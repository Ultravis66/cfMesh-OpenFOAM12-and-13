/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "refineBoundaryLayers.H"
#include "boundaryLayerConstraintPlanner.H"
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "polyMeshGenAddressing.H"
#include "polyMeshGen2DEngine.H"
#include "VRWGraphList.H"
#include "meshSurfacePartitioner.H"
#include "detectBoundaryLayers.H"

#include "labelledPair.H"
#include "labelledScalar.H"

# ifdef USE_OMP
#include <omp.h>
#include <map>
#include <utility>
# endif

//#define DEBUGLayer

# ifdef DEBUGLayer
#include "OFstream.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool refineBoundaryLayers::analyseLayers()
{
    const meshSurfaceEngine& mse = surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& facePatch = mse.boundaryFacePatches();

    meshSurfacePartitioner mPart(mse);
    //- Partial-patch BL policy:
    //- use the existing effective requested layer count as the scope.
    //- requestedLayers > 1  => partial detected layer is allowed
    //- requestedLayers <= 1 => preserve legacy all-or-nothing behaviour
    const PtrList<boundaryPatch>& partialBoundaries = mesh_.boundaries();
    boolList allowPartialLayerPatch(partialBoundaries.size(), false);

    forAll(partialBoundaries, patchI)
    {
        label requestedLayers = globalNumLayers_;
        const word pName = partialBoundaries[patchI].patchName();

        std::map<word, label>::const_iterator it =
            numLayersForPatch_.find(pName);

        if( it != numLayersForPatch_.end() )
            requestedLayers = it->second;

        allowPartialLayerPatch[patchI] = (requestedLayers > 1);

        if( allowPartialLayerPatch[patchI] )
        {
            Info << "BLPARTIALPATCH eligible patch " << pName
                 << " requestedLayers=" << requestedLayers << endl;
        }
    }

    detectBoundaryLayers dbl
    (
        mPart,
        is2DMesh_,
        allowPartialLayerPatch
    );

    const label nGroups = dbl.nDistinctLayers();
    const labelList& faceInLayer = dbl.faceInLayer();

    //- get the hair edges
    splitEdges_ = dbl.hairEdges();

    //- Transaction diagnostic: distinguish "no layers detected" from
    //- "layers detected but hair-edge generation failed".  This runs for
    //- both the ordinary pass and a restored pass-2 mesh.
    label nDetectedLayerFaces = 0;
    forAll(faceInLayer, bfI)
        if( faceInLayer[bfI] >= 0 )
            ++nDetectedLayerFaces;

    Info << "BLANALYSE"
         << " meshPoints=" << mesh_.points().size()
         << " meshFaces=" << mesh_.faces().size()
         << " meshCells=" << mesh_.cells().size()
         << " boundaryFaces=" << bFaces.size()
         << " nGroups=" << nGroups
         << " layerFaces=" << nDetectedLayerFaces
         << " hairEdges=" << splitEdges_.size()
         << endl;

    //- Zero hair edges cannot support refineBoundaryLayers.  Previously
    //- the validation loop below could vacuously succeed when every
    //- faceInLayer entry was negative, allowing generateNewFaces() to run
    //- against empty split-edge metadata.
    //
    //- Do NOT mark refinementValid_ false here: "nothing detected to refine"
    //- and "metadata became structurally inconsistent" are different states.
    //- The caller separately checks refinementCompleted(), which remains
    //- false because done_ is set only after generateNewCells() completes.
    if( splitEdges_.size() == 0 )
    {
        Info << "BLANALYSE_NO_HAIR_EDGES"
             << " nGroups=" << nGroups
             << " layerFaces=" << nDetectedLayerFaces
             << " boundaryFaces=" << bFaces.size()
             << " -- refinement cannot execute"
             << endl;
        return false;
    }

    # ifdef DEBUGLayer
    OFstream file("hairEdges.vtk");

    //- write the header
    file << "# vtk DataFile Version 3.0\n";
    file << "vtk output\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    //- write points
    file << "POINTS " << 2*splitEdges_.size() << " float\n";
    forAll(splitEdges_, seI)
    {
        const point& p = mse.mesh().points()[splitEdges_[seI].start()];

        file << p.x() << ' ' << p.y() << ' ' << p.z() << nl;

        const point op = mse.mesh().points()[splitEdges_[seI].end()];

        file << op.x() << ' ' << op.y() << ' ' << op.z() << nl;
    }

    //- write lines
    file << "\nLINES " << splitEdges_.size()
         << " " << 3*splitEdges_.size() << nl;
    forAll(splitEdges_, eI)
    {
        file << 2 << " " << 2*eI << " " << (2*eI+1) << nl;
    }

    file << "\n";
    # endif

    //- create point to split edges addressing
    splitEdgesAtPoint_.reverseAddressing(splitEdges_);

    //- check if the layer is valid
    bool validLayer(true);
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40)
    # endif
    forAll(faceInLayer, bfI)
    {
        if( faceInLayer[bfI] < 0 )
            continue;

        const face& bf = bFaces[bfI];

        forAll(bf, pI)
            if( splitEdgesAtPoint_.sizeOfRow(bf[pI]) == 0 )
                validLayer = false;
    }

    # ifdef DEBUGLayer
    Info << "Number of independent layers in the mesh is " << nGroups << endl;
    Info << "Is valid layer " << validLayer << endl;
    # endif

    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

    //- create patch name to index addressing
    std::map<word, label> patchNameToIndex;
    forAll(boundaries, patchI)
        patchNameToIndex[boundaries[patchI].patchName()] = patchI;

    //- check layer labels over a patch
    layerAtPatch_.setSize(boundaries.size());
    forAll(layerAtPatch_, i)
        layerAtPatch_[i].clear();
    List<DynList<label> > groupsAtPatch(boundaries.size());
    forAll(faceInLayer, bfI)
        groupsAtPatch[facePatch[bfI]].appendIfNotIn(faceInLayer[bfI]);

    //- set the information which patches have an extruded layer
    forAll(groupsAtPatch, patchI)
    {
        const DynList<label>& layers = groupsAtPatch[patchI];

        forAll(layers, i)
        {
            if( layers[i] < 0 )
            {
                if
                (
                    patchI >= 0 &&
                    patchI < label(allowPartialLayerPatch.size()) &&
                    allowPartialLayerPatch[patchI]
                )
                {
                    //- Keep the valid detected portion of an explicitly
                    //- requested BL patch. The local unsupported faces
                    //- remain faceInLayer < 0 and are protected below.
                    continue;
                }

                //- Legacy all-or-nothing behaviour for non-BL patches.
                layerAtPatch_[patchI].clear();
                break;
            }
            else
            {
                layerAtPatch_[patchI].append(layers[i]);
            }
        }
    }

    # ifdef DEBUGLayer
    Info << "Layer at patch " << layerAtPatch_ << endl;
    # endif

    //- set the information which patches are a single boundary layer face
    patchesInLayer_.setSize(nGroups);
    forAll(layerAtPatch_, patchI)
    {
        const DynList<label>& layers = layerAtPatch_[patchI];

        forAll(layers, i)
            patchesInLayer_[layers[i]].append
            (
                boundaries[patchI].patchName()
            );
    }

    # ifdef DEBUGLayer
    Info << "Patches in layer " << patchesInLayer_ << endl;
    # endif

    //- set the number of boundary layers for each patch
    labelList nLayersAtPatch(layerAtPatch_.size(), -1);
    boolList protectedValue(layerAtPatch_.size(), false);

    forAll(patchesInLayer_, layerI)
    {
        const DynList<word>& layerPatches = patchesInLayer_[layerI];

        label maxNumLayers(1);
        bool hasLocalValue(false);

        //- find the maximum requested number of layers over the layer
        forAll(layerPatches, lpI)
        {
            const word pName = layerPatches[lpI];

            std::map<word, label>::const_iterator it =
                numLayersForPatch_.find(pName);

            if( it != numLayersForPatch_.end() )
            {
                const label ptchI = patchNameToIndex[pName];

                // CFMitch/shared-modern hard no-BL invariant:
                //
                // setNumberOfLayersForPatch() documents nLayers < 2
                // as disabling boundary layers on the patch.  The original
                // cfMesh layer-group logic allowed such a patch to inherit
                // maxNumLayers from another patch in the same detected
                // layer group, which could silently turn an explicit
                // nLayers=0 surface into a full BL-generating surface.
                //
                // Internally one layer means "retain the unsplit parent",
                // therefore normalize every explicit nLayers < 2 request
                // to one and protect it from group maximum propagation.
                //
                // This is shared by legacyEnhanced and constraintPlanner.
                // classicCfMesh uses its separate frozen implementation.
                if( it->second < 2 )
                {
                    nLayersAtPatch[ptchI] = 1;
                    protectedValue[ptchI] = true;
                    hasLocalValue = true;
                }
                //- check if the layer is interrupted at this patch
                else if
                (
                    discontinuousLayersForPatch_.find(pName) !=
                    discontinuousLayersForPatch_.end()
                )
                {
                    //- set the number of layers and lock this location
                    nLayersAtPatch[ptchI] = it->second;
                    protectedValue[ptchI] = true;
                    hasLocalValue = true;
                }
                else
                {
                    //- take the maximum number of layers
                    maxNumLayers = Foam::max(maxNumLayers, it->second);
                    hasLocalValue = true;
                }
            }
        }

        //- apply the global value if no local values exist
        if( !hasLocalValue )
            maxNumLayers = globalNumLayers_;

        //- apply the maximum number of ayer of all unprotected patches
        forAll(layerPatches, lpI)
        {
            const label ptchI = patchNameToIndex[layerPatches[lpI]];

            if( !protectedValue[ptchI] )
                nLayersAtPatch[ptchI] = maxNumLayers;
        }
    }

    if( is2DMesh_ )
    {
        polyMeshGen2DEngine mesh2DEngine(mesh_);
        const boolList& zMinPoint = mesh2DEngine.zMinPoints();
        const boolList& zMaxPoint = mesh2DEngine.zMaxPoints();

        const faceList::subList& bFaces = mse.boundaryFaces();

        boolList allZMax(mesh_.boundaries().size(), true);
        boolList allZMin(mesh_.boundaries().size(), true);

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(bFaces, bfI)
        {
            const face& bf = bFaces[bfI];

            forAll(bf, pI)
            {
                if( !zMinPoint[bf[pI]] )
                    allZMin[facePatch[bfI]] = false;
                if( !zMaxPoint[bf[pI]] )
                    allZMax[facePatch[bfI]] = false;
            }
        }

        //- mark empty patches as already used
        forAll(allZMin, patchI)
        {
            if( allZMin[patchI] ^ allZMax[patchI] )
            {
                nLayersAtPatch[patchI] = -1;
                layerAtPatch_[patchI].clear();
            }
        }
    }

    //- perform reduction over all processors
    reduce(nLayersAtPatch, maxOp<labelList>());

    // CFMitch/shared-modern hard no-BL audit.
    // An explicitly disabled patch must resolve to <= 1 and may never
    // inherit an active layer count from its detected layer group.
    {
        label nHardNoBLPatches = 0;
        label nHardNoBLSafe = 0;
        label nHardNoBLViolations = 0;

        forAll(nLayersAtPatch, patchI)
        {
            const word& pName = boundaries[patchI].patchName();

            const std::map<word, label>::const_iterator it =
                numLayersForPatch_.find(pName);

            if
            (
                it == numLayersForPatch_.end()
             || it->second >= 2
            )
                continue;

            ++nHardNoBLPatches;

            if( nLayersAtPatch[patchI] <= 1 )
                ++nHardNoBLSafe;
            else
                ++nHardNoBLViolations;
        }

        Info
            << "BL_HARD_NOBL_AUDIT"
            << " requestedPatches=" << nHardNoBLPatches
            << " safe=" << nHardNoBLSafe
            << " violations=" << nHardNoBLViolations
            << endl;
    }

    # ifdef DEBUGLayer
    Pout << "nLayersAtPatch " << nLayersAtPatch << endl;
    # endif

    //- set the number of boundary layers which shall be generated above
    //- each boundary face
    nLayersAtBndFace_.setSize(facePatch.size());
    nLayersAtBndFace_ = globalNumLayers_;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(nLayersAtBndFace_, bfI)
    {
        const label patchI = facePatch[bfI];

        //- A locally unsupported detected-layer face must never be sent
        //- into boundary-layer refinement. This is the critical guard
        //- which lets the rest of a configured BL patch refine safely.
        if( faceInLayer[bfI] < 0 )
        {
            nLayersAtBndFace_[bfI] = 1;
            continue;
        }

        if( nLayersAtPatch[patchI] < 0 )
        {
            nLayersAtBndFace_[bfI] = 1;
        }
        else
        {
            nLayersAtBndFace_[bfI] = nLayersAtPatch[patchI];

            if( specialMode_ )
            {
                ++nLayersAtBndFace_[bfI];
            }
        }
    }

    // BL/BL junction ramp: faces touching sharp junction points get
    // reduced layer count to prevent degenerate cells at blade/hub/shroud.
    // Ring 0 (junction face)    -> 1 layer
    // Ring 1 (neighbor face)    -> max 2 layers
    // Ring 2+ resumes full nLayers
    //
    // Experimental baseline synchronization:
    // the installed test library used for the 15-layer contact-smoother
    // A/B has this historical BL/BL layer-count ramp disabled.  Keep the
    // source behavior identical while testing the contact sweep transaction.
    const bool enableBlblJunctionLayerCountRamp = false;

    if( !enableBlblJunctionLayerCountRamp )
        Info << "BL/BL junction ramp: DISABLED A/B" << endl;

    if
    (
        enableBlblJunctionLayerCountRamp
     && blblJunctionPoints_.size() > 0
    )
    {
        const meshSurfaceEngine& mseLoc = surfaceEngine();
        const VRWGraph& ptFaces = mseLoc.pointFaces();

        // Build mesh-point to boundary-point map
        labelList meshToBnd(mesh_.points().size(), -1);
        const labelList& bPts = mseLoc.boundaryPoints();
        forAll(bPts, bpI)
            meshToBnd[bPts[bpI]] = bpI;

        // Ring 0: faces directly touching junction points -> 1 layer
        // Only reduce faces that already have more than 1 layer
        boolList ring0face(nLayersAtBndFace_.size(), false);
        forAllConstIter(labelHashSet, blblJunctionPoints_, it)
        {
            const label bpI = it.key();
            if( bpI < 0 || bpI >= label(ptFaces.size()) ) continue;
            // Skip finite ramp seeds -- not hard junction caps
            if( rampSeedPoints_.size() > bpI && rampSeedPoints_[bpI] ) continue;
            forAllRow(ptFaces, bpI, pfI)
            {
                const label bfI = ptFaces(bpI, pfI);
                if( bfI < 0 || bfI >= label(nLayersAtBndFace_.size()) ) continue;
                // Virtual topology takes priority -- skip VT-handled faces
                if( bfI < label(vtFaceRing_.size()) && vtFaceRing_[bfI] >= 0 ) continue;
                if( nLayersAtBndFace_[bfI] > 1 )
                {
                    ring0face[bfI] = true;
                    nLayersAtBndFace_[bfI] = 1;
                }
            }
        }

        // Ring 1: faces sharing a point with ring0 faces -> max 2 layers
        // Only reduce faces that already have more than 2 layers
        const faceList::subList& bFacesLoc = mseLoc.boundaryFaces();
        boolList ring1face(nLayersAtBndFace_.size(), false);
        forAll(nLayersAtBndFace_, bfI)
        {
            if( !ring0face[bfI] ) continue;
            const face& f = bFacesLoc[bfI];
            forAll(f, pI)
            {
                const label meshPtI = f[pI];
                if( meshPtI < 0 || meshPtI >= label(meshToBnd.size()) ) continue;
                const label bpI = meshToBnd[meshPtI];
                if( bpI < 0 || bpI >= label(ptFaces.size()) ) continue;
                forAllRow(ptFaces, bpI, pfI)
                {
                    const label nbfI = ptFaces(bpI, pfI);
                    if( nbfI < 0 || nbfI >= label(nLayersAtBndFace_.size()) ) continue;
                    if( ring0face[nbfI] ) continue;
                    // Virtual topology takes priority -- skip VT-handled faces
                    if( nbfI < label(vtFaceRing_.size()) && vtFaceRing_[nbfI] >= 0 ) continue;
                    if( nLayersAtBndFace_[nbfI] > 2 )
                    {
                        ring1face[nbfI] = true;
                        nLayersAtBndFace_[nbfI] = 2;
                    }
                }
            }
        }

        label nRing0 = 0, nRing1 = 0;
        forAll(nLayersAtBndFace_, bfI)
        {
            if( ring0face[bfI] ) ++nRing0;
            else if( ring1face[bfI] ) ++nRing1;
        }
        Info << "BL/BL junction ramp: ring0=" << nRing0
             << " faces forced to 1 layer, ring1=" << nRing1
             << " faces capped at 2 layers" << endl;
    }

    // Acute corner face cap: controlled by acuteCornerCapLayers_ member
    // Set via setAcuteCornerCapLayers() before refineLayers()
    Info << "Acute corner face cap switch: " << acuteCornerCapLayers_
         << " acutePts=" << blblAcuteCornerPoints_.size() << endl;
    if( acuteCornerCapLayers_ && blblAcuteCornerPoints_.size() > 0 )
    {
        const meshSurfaceEngine& mseLoc = surfaceEngine();
        const VRWGraph& ptFaces = mseLoc.pointFaces();
        const faceList::subList& bFacesLoc = mseLoc.boundaryFaces();
        labelList meshToBnd2(mesh_.points().size(), -1);
        const labelList& bPts2 = mseLoc.boundaryPoints();
        forAll(bPts2, bpI)
            meshToBnd2[bPts2[bpI]] = bpI;

        boolList acRing0(nLayersAtBndFace_.size(), false);
        forAllConstIter(labelHashSet, blblJunctionPoints_, it)
        {
            const label bpI = it.key();
            if( bpI < 0 || bpI >= label(ptFaces.size()) ) continue;
            // Skip finite ramp seeds -- not hard junction caps
            if( rampSeedPoints_.size() > bpI && rampSeedPoints_[bpI] ) continue;
            forAllRow(ptFaces, bpI, pfI)
            {
                const label bfI = ptFaces(bpI, pfI);
                if( bfI < 0 || bfI >= label(nLayersAtBndFace_.size()) ) continue;
                if( nLayersAtBndFace_[bfI] > 1 )
                    nLayersAtBndFace_[bfI] = 1;
                acRing0[bfI] = true;
            }
        }

        boolList acRing1(nLayersAtBndFace_.size(), false);
        forAll(nLayersAtBndFace_, bfI)
        {
            if( !acRing0[bfI] ) continue;
            const face& f = bFacesLoc[bfI];
            forAll(f, pI)
            {
                const label meshPtI = f[pI];
                if( meshPtI < 0 || meshPtI >= label(meshToBnd2.size()) ) continue;
                const label bpI = meshToBnd2[meshPtI];
                if( bpI < 0 || bpI >= label(ptFaces.size()) ) continue;
                forAllRow(ptFaces, bpI, pfI)
                {
                    const label nbfI = ptFaces(bpI, pfI);
                    if( nbfI < 0 || nbfI >= label(nLayersAtBndFace_.size()) ) continue;
                    if( acRing0[nbfI] ) continue;
                    if( nLayersAtBndFace_[nbfI] > 1 )
                        nLayersAtBndFace_[nbfI] = 1;
                    acRing1[nbfI] = true;
                }
            }
        }

        boolList acRing2(nLayersAtBndFace_.size(), false);
        forAll(nLayersAtBndFace_, bfI)
        {
            if( !acRing1[bfI] ) continue;
            const face& f = bFacesLoc[bfI];
            forAll(f, pI)
            {
                const label meshPtI = f[pI];
                if( meshPtI < 0 || meshPtI >= label(meshToBnd2.size()) ) continue;
                const label bpI = meshToBnd2[meshPtI];
                if( bpI < 0 || bpI >= label(ptFaces.size()) ) continue;
                forAllRow(ptFaces, bpI, pfI)
                {
                    const label nbfI = ptFaces(bpI, pfI);
                    if( nbfI < 0 || nbfI >= label(nLayersAtBndFace_.size()) ) continue;
                    if( acRing0[nbfI] || acRing1[nbfI] ) continue;
                    if( nLayersAtBndFace_[nbfI] > 2 )
                        nLayersAtBndFace_[nbfI] = 2;
                    acRing2[nbfI] = true;
                }
            }
        }

        label nAC0=0, nAC1=0, nAC2=0;
        forAll(nLayersAtBndFace_, bfI)
        {
            if( acRing0[bfI] ) ++nAC0;
            else if( acRing1[bfI] ) ++nAC1;
            else if( acRing2[bfI] ) ++nAC2;
        }
        Info << "Acute corner face cap: ring0=" << nAC0
             << " ring1=" << nAC1
             << " ring2=" << nAC2
             << " faces capped" << endl;
    }

    // Provenance/junction BL retraction: cap selected boundary faces to
    // a maximum layer count. Generalizes the old binary
    // forceSingleLayerFaces_ logic into a tapered ring cap.
    // Applied after all global/patch layer-count logic so it always wins.
    actuallyCappedFaces_.clear();
    actuallyCappedProviderPlansAtFace_.clear();

    if( forcedMaxLayersAtFace_.size() )
    {
        label nCapped = 0;
        label nAlreadyAtOrBelow = 0;
        label nOutOfRange = 0;
        label nInvalidCap = 0;
        label nCappedNoProvider = 0;

        forAllConstIter(Map<label>, forcedMaxLayersAtFace_, it)
        {
            const label bfI = it.key();
            const label maxLayers = it();

            if( bfI < 0 || bfI >= label(nLayersAtBndFace_.size()) )
            {
                ++nOutOfRange;
                continue;
            }

            if( maxLayers < 0 )
            {
                ++nInvalidCap;
                continue;
            }

            if( nLayersAtBndFace_[bfI] > maxLayers )
            {
                nLayersAtBndFace_[bfI] = maxLayers;
                ++nCapped;

                actuallyCappedFaces_.insert(bfI);

                if( forcedMaxLayerProviderPlansAtFace_.found(bfI) )
                {
                    actuallyCappedProviderPlansAtFace_.insert
                    (
                        bfI,
                        forcedMaxLayerProviderPlansAtFace_[bfI]
                    );
                }
                else
                {
                    ++nCappedNoProvider;
                }
            }
            else
            {
                ++nAlreadyAtOrBelow;
            }
        }

        Info << "refineBoundaryLayers: forced max-layer caps applied: "
             << "capped=" << nCapped
             << " alreadyAtOrBelow=" << nAlreadyAtOrBelow
             << " outOfRange=" << nOutOfRange
             << " invalidCap=" << nInvalidCap
             << " requested=" << forcedMaxLayersAtFace_.size()
             << " effectiveFaces=" << actuallyCappedFaces_.size()
             << " providerFaces="
             << actuallyCappedProviderPlansAtFace_.size()
             << " cappedNoProvider=" << nCappedNoProvider
             << endl;
    }

    # ifdef DEBUGLayer
    forAll(nLayersAtBndFace_, bfI)
    Pout << "Boundary face " << bfI << " in patch "
        << facePatch[bfI] << " num layers " << nLayersAtBndFace_[bfI] << endl;
    //::exit(1);
    # endif

    return validLayer;
}

void refineBoundaryLayers::generateNewVertices()
{
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    pointFieldPMG& points = mesh_.points();

    const meshSurfaceEngine& mse = surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const VRWGraph& pointFaces = mse.pointFaces();
    const labelList& facePatch = mse.boundaryFacePatches();
    const labelList& bp = mse.bp();

    // CFMitch v0 architecture entry point.
    // The planner is behaviour-neutral at this stage.  The
    // constraintPlanner mode intentionally retains legacyEnhanced
    // execution so exact mesh equivalence can be established before
    // responsibilities are migrated into the planner.
    boundaryLayerConstraintPlanner cfmitchPlanner
    (
        mesh_,
        boundaryLayerArchitecture_
    );

    cfmitchPlanner.report();

    //- allocate the data from storing parameters applying to a split edge
    LongList<scalar> firstLayerThickness(splitEdges_.size());
    LongList<scalar> thicknessRatio(splitEdges_.size());
    labelLongList nNodesAtEdge(splitEdges_.size());
    labelLongList nLayersAtEdge(splitEdges_.size());

    // Per-split-edge thickness scale from BLRepairPlan.
    // Min scale wins if multiple capped faces touch the same edge.
    LongList<scalar> forcedThicknessScaleAtEdge(splitEdges_.size(), scalar(1.0));

    if( forcedThicknessScaleAtFace_.size() )
    {
        const faceList::subList& bFaces = mse.boundaryFaces();

        label nScaledFaces = 0;
        label nScaledEdges = 0;
        label nOutOfRangeFaces = 0;

        forAllConstIter(Map<scalar>, forcedThicknessScaleAtFace_, it)
        {
            const label bfI = it.key();
            const scalar scale =
                Foam::max(scalar(0.0), Foam::min(scalar(1.0), it()));

            if( bfI < 0 || bfI >= label(bFaces.size()) )
            {
                ++nOutOfRangeFaces;
                continue;
            }

            ++nScaledFaces;

            const face& bf = bFaces[bfI];

            forAll(bf, fpI)
            {
                const label bpI = bf[fpI];

                if( bpI < 0 || bpI >= label(splitEdgesAtPoint_.size()) )
                    continue;

                forAllRow(splitEdgesAtPoint_, bpI, seRowI)
                {
                    const label seI = splitEdgesAtPoint_(bpI, seRowI);

                    if( seI < 0 || seI >= label(forcedThicknessScaleAtEdge.size()) )
                        continue;

                    const scalar oldScale = forcedThicknessScaleAtEdge[seI];
                    const scalar newScale = Foam::min(oldScale, scale);

                    if( newScale < oldScale - SMALL )
                        ++nScaledEdges;

                    forcedThicknessScaleAtEdge[seI] = newScale;
                }
            }
        }

        //- Do the requested faces actually own split edges? If their
        //- points have no splitEdgesAtPoint_ rows, no cap and no scale can
        //- ever reach an edge, and the face is invisible to refinement.
        //- analyseLayers() only validates faces with faceInLayer >= 0, so
        //- provenance seed faces outside a detected layer are never checked.
        {
            label nFacesNoRows = 0, nPtsNoRows = 0, nPtsWithRows = 0;
            forAllConstIter(Map<scalar>, forcedThicknessScaleAtFace_, dit)
            {
                const label dbfI = dit.key();
                if( dbfI < 0 || dbfI >= label(bFaces.size()) ) continue;
                const face& dbf = bFaces[dbfI];
                bool anyRow = false;
                forAll(dbf, dpI)
                {
                    const label gp = dbf[dpI];
                    if( gp < 0 || gp >= label(splitEdgesAtPoint_.size()) )
                        { ++nPtsNoRows; continue; }
                    if( splitEdgesAtPoint_.sizeOfRow(gp) == 0 ) ++nPtsNoRows;
                    else { ++nPtsWithRows; anyRow = true; }
                }
                if( !anyRow ) ++nFacesNoRows;
            }
            Info << "refineBoundaryLayers: seed split-edge coverage:"
                 << " facesWithNoSplitEdges=" << nFacesNoRows
                 << " pointsWithRows=" << nPtsWithRows
                 << " pointsWithoutRows=" << nPtsNoRows
                 << " splitEdgesAtPointSize=" << splitEdgesAtPoint_.size()
                 << " splitEdges=" << splitEdges_.size() << endl;
        }

        Info << "refineBoundaryLayers: forced thickness scales prepared: "
             << "scaledFaces=" << nScaledFaces
             << " scaledEdgesUpdates=" << nScaledEdges
             << " outOfRangeFaces=" << nOutOfRangeFaces
             << " requestedFaces=" << forcedThicknessScaleAtFace_.size()
             << endl;
    }

    if( cfmitchPlanner.plannerEnabled() )
    {
        const boundaryLayerPlan cfmitchPlan =
            cfmitchPlanner.solveLayerCounts
            (
                mse,
                splitEdges_,
                neutralLayerScaleAtMeshPoint_,
                nLayersAtBndFace_,
                globalThicknessRatio_,
                thicknessRatioForPatch_,
                vtFaceRing_,
                actuallyCappedFaces_,
                constraintPlannerMaxLayerStep_
            );

        nLayersAtBndFace_ =
            cfmitchPlan.faceLayers();

        Info
            << "CFMITCH: resolved layer-count plan applied"
            << " faces=" << nLayersAtBndFace_.size()
            << endl;

        // -------------------------------------------------------------
        // CFMitch V3.8 -- pre-topology edge-hex compatibility planner.
        //
        // refineEdgeHexCell can refine a parent touched by two active
        // BL boundary faces only when each generating face has one
        // well-defined opposite quadrilateral terminal face.
        //
        // A six-face polyhedron is not necessarily a topological hex.
        // In particular, a polygonal opposite face can contain an
        // edge-conforming chain of N-1 intermediate vertices while
        // remaining a single face. Assigning that whole polygon to one
        // of N terminal children leaves exactly N-1 open child edges.
        //
        // Detect that condition using only the unmodified Q0 topology.
        // Both generating boundary faces are locally terminated at one
        // layer, then the normal constraint planner is re-run so the
        // requested layer count recovers smoothly according to
        // constraintPlannerMaxLayerStep_.
        // -------------------------------------------------------------
        {
            const labelList& structuralFaceOwners =
                mse.faceOwners();

            const cellListPMG& structuralCells =
                mesh_.cells();

            const faceListPMG& structuralFaces =
                mesh_.faces();

            const PtrList<boundaryPatch>& structuralBoundaries =
                mesh_.boundaries();

            const label structuralStartBoundary =
                structuralBoundaries.size()
              ? structuralBoundaries[0].patchStart()
              : structuralFaces.size();

            auto structuralCommonEdgeCount =
            [] (const face& a, const face& b) -> label
            {
                label count = 0;

                forAll(a, aeI)
                {
                    const label a0 = a[aeI];
                    const label a1 = a[(aeI+1) % a.size()];
                    const label alo = Foam::min(a0, a1);
                    const label ahi = Foam::max(a0, a1);

                    forAll(b, beI)
                    {
                        const label b0 = b[beI];
                        const label b1 = b[(beI+1) % b.size()];

                        if
                        (
                            alo == Foam::min(b0, b1)
                         && ahi == Foam::max(b0, b1)
                        )
                        {
                            ++count;
                        }
                    }
                }

                return count;
            };

            label structuralMaxLayers = 1;

            forAll(nLayersAtBndFace_, bfI)
            {
                structuralMaxLayers =
                    Foam::max
                    (
                        structuralMaxLayers,
                        nLayersAtBndFace_[bfI]
                    );
            }

            const label structuralMaxPasses =
                structuralMaxLayers + 2;

            labelHashSet structuralAllSeeds;
            label structuralPass = 0;
            bool structuralChanged = true;

            while
            (
                structuralChanged
             && structuralPass < structuralMaxPasses
            )
            {
                ++structuralPass;
                structuralChanged = false;

                labelList activeCount
                (
                    structuralCells.size(),
                    0
                );

                labelList activeBf0
                (
                    structuralCells.size(),
                    -1
                );

                labelList activeBf1
                (
                    structuralCells.size(),
                    -1
                );

                forAll(structuralFaceOwners, bfI)
                {
                    if
                    (
                        bfI < 0
                     || bfI >= label(nLayersAtBndFace_.size())
                     || nLayersAtBndFace_[bfI] <= 1
                    )
                    {
                        continue;
                    }

                    const label cellI =
                        structuralFaceOwners[bfI];

                    if
                    (
                        cellI < 0
                     || cellI >= label(structuralCells.size())
                    )
                    {
                        continue;
                    }

                    if( activeCount[cellI] == 0 )
                    {
                        activeBf0[cellI] = bfI;
                    }
                    else if( activeCount[cellI] == 1 )
                    {
                        activeBf1[cellI] = bfI;
                    }

                    ++activeCount[cellI];
                }

                labelHashSet structuralPassSeeds;

                label nType2Parents = 0;
                label nBadParentFaceCount = 0;
                label nBadActiveFaces = 0;
                label nBadActiveCommonEdges = 0;
                label nBadOppositeCount = 0;
                label nNonQuadOpposite = 0;

                forAll(activeCount, cellI)
                {
                    if( activeCount[cellI] != 2 )
                        continue;

                    ++nType2Parents;

                    const label bf0 = activeBf0[cellI];
                    const label bf1 = activeBf1[cellI];

                    bool compatible = true;

                    const cell& parent =
                        structuralCells[cellI];

                    if( parent.size() != 6 )
                    {
                        compatible = false;
                        ++nBadParentFaceCount;
                    }

                    const label activeFaceI0 =
                        structuralStartBoundary + bf0;

                    const label activeFaceI1 =
                        structuralStartBoundary + bf1;

                    if
                    (
                        compatible
                     &&
                        (
                            activeFaceI0 < 0
                         || activeFaceI1 < 0
                         || activeFaceI0 >= label(structuralFaces.size())
                         || activeFaceI1 >= label(structuralFaces.size())
                        )
                    )
                    {
                        compatible = false;
                        ++nBadActiveFaces;
                    }

                    bool parentHasActive0 = false;
                    bool parentHasActive1 = false;

                    if( compatible )
                    {
                        forAll(parent, localFaceI)
                        {
                            parentHasActive0 =
                                parentHasActive0
                             || parent[localFaceI] == activeFaceI0;

                            parentHasActive1 =
                                parentHasActive1
                             || parent[localFaceI] == activeFaceI1;
                        }

                        if
                        (
                            !parentHasActive0
                         || !parentHasActive1
                        )
                        {
                            compatible = false;
                            ++nBadActiveFaces;
                        }
                    }

                    if
                    (
                        compatible
                     && structuralCommonEdgeCount
                        (
                            structuralFaces[activeFaceI0],
                            structuralFaces[activeFaceI1]
                        ) != 1
                    )
                    {
                        compatible = false;
                        ++nBadActiveCommonEdges;
                    }

                    label oppositeFaceI0 = -1;
                    label oppositeFaceI1 = -1;
                    label nOpposite0 = 0;
                    label nOpposite1 = 0;

                    if( compatible )
                    {
                        const face& active0 =
                            structuralFaces[activeFaceI0];

                        const face& active1 =
                            structuralFaces[activeFaceI1];

                        forAll(parent, localFaceI)
                        {
                            const label candidateFaceI =
                                parent[localFaceI];

                            if
                            (
                                candidateFaceI < 0
                             || candidateFaceI >=
                                label(structuralFaces.size())
                            )
                            {
                                compatible = false;
                                ++nBadActiveFaces;
                                break;
                            }

                            if
                            (
                                candidateFaceI != activeFaceI0
                             && structuralCommonEdgeCount
                                (
                                    active0,
                                    structuralFaces[candidateFaceI]
                                ) == 0
                            )
                            {
                                oppositeFaceI0 = candidateFaceI;
                                ++nOpposite0;
                            }

                            if
                            (
                                candidateFaceI != activeFaceI1
                             && structuralCommonEdgeCount
                                (
                                    active1,
                                    structuralFaces[candidateFaceI]
                                ) == 0
                            )
                            {
                                oppositeFaceI1 = candidateFaceI;
                                ++nOpposite1;
                            }
                        }
                    }

                    if
                    (
                        compatible
                     &&
                        (
                            nOpposite0 != 1
                         || nOpposite1 != 1
                        )
                    )
                    {
                        compatible = false;
                        ++nBadOppositeCount;
                    }

                    if
                    (
                        compatible
                     &&
                        (
                            structuralFaces[oppositeFaceI0].size() != 4
                         || structuralFaces[oppositeFaceI1].size() != 4
                        )
                    )
                    {
                        compatible = false;
                        ++nNonQuadOpposite;
                    }

                    if( !compatible )
                    {
                        structuralPassSeeds.insert(bf0);
                        structuralPassSeeds.insert(bf1);
                    }
                }

                Info
                    << "CFMITCH V3.8 EDGEHEX PLAN GUARD:"
                    << " pass=" << structuralPass
                    << " type2Parents=" << nType2Parents
                    << " seedFaces=" << structuralPassSeeds.size()
                    << " badParentFaceCount="
                    << nBadParentFaceCount
                    << " badActiveFaces="
                    << nBadActiveFaces
                    << " badActiveCommonEdges="
                    << nBadActiveCommonEdges
                    << " badOppositeCount="
                    << nBadOppositeCount
                    << " nonQuadOpposite="
                    << nNonQuadOpposite
                    << endl;

                if( structuralPassSeeds.empty() )
                    break;

                structuralChanged = true;

                forAllConstIter
                (
                    labelHashSet,
                    structuralPassSeeds,
                    seedIt
                )
                {
                    const label bfI = seedIt.key();

                    if
                    (
                        bfI < 0
                     || bfI >= label(nLayersAtBndFace_.size())
                    )
                    {
                        continue;
                    }

                    nLayersAtBndFace_[bfI] = 1;
                    structuralAllSeeds.insert(bfI);
                }

                const boundaryLayerPlan structuralPlan =
                    cfmitchPlanner.solveLayerCounts
                    (
                        mse,
                        splitEdges_,
                        neutralLayerScaleAtMeshPoint_,
                        nLayersAtBndFace_,
                        globalThicknessRatio_,
                        thicknessRatioForPatch_,
                        vtFaceRing_,
                        structuralAllSeeds,
                        constraintPlannerMaxLayerStep_
                    );

                nLayersAtBndFace_ =
                    structuralPlan.faceLayers();

                // Structural termination is a hard upper bound. The
                // planner is expected to preserve it; enforce it again
                // so a future planner implementation cannot reactivate
                // an unsupported edge-hex parent.
                forAllConstIter
                (
                    labelHashSet,
                    structuralAllSeeds,
                    seedIt
                )
                {
                    const label bfI = seedIt.key();

                    if
                    (
                        bfI >= 0
                     && bfI < label(nLayersAtBndFace_.size())
                    )
                    {
                        nLayersAtBndFace_[bfI] = 1;
                    }
                }
            }

            Info
                << "CFMITCH V3.8 EDGEHEX PLAN SUMMARY:"
                << " passes=" << structuralPass
                << " structuralSeedFaces="
                << structuralAllSeeds.size()
                << " maxPasses="
                << structuralMaxPasses
                << endl;
        }
    }
    else
    {
        // Dense stable mesh-point lookup for BL/neutral height scales.
        // Built serially; the OMP split-edge loop below is read-only.
        scalarField neutralScaleByMeshPoint
        (
            mesh_.points().size(),
            scalar(1)
        );

        label nNeutralScaleMapped = 0;

        forAllConstIter
        (
            Map<scalar>,
            neutralLayerScaleAtMeshPoint_,
            it
        )
        {
            const label meshPtI = it.key();

            if
            (
                meshPtI < 0
             || meshPtI >= label(neutralScaleByMeshPoint.size())
            )
                continue;

            neutralScaleByMeshPoint[meshPtI] =
                Foam::max
                (
                    scalar(0),
                    Foam::min(scalar(1), it())
                );

            ++nNeutralScaleMapped;
        }

        Info << "refineBoundaryLayers: neutral scale metadata:"
             << " input=" << neutralLayerScaleAtMeshPoint_.size()
             << " mapped=" << nNeutralScaleMapped
             << endl;

        // BL/neutral height-aware FACE layer-count adaptation.
        //
        // neutralScaleByMeshPoint is keyed by the ORIGINAL stable wall point.
        // After BL extrusion that point is normally the interior endpoint of a
        // split hair; the current boundary faces contain the opposite/top
        // endpoint. Transfer the scale along the hair first, then modify the
        // canonical nLayersAtBndFace_ topology.
        scalarField neutralScaleByCurrentBp
        (
            mse.boundaryPoints().size(),
            scalar(1)
        );

        boolList neutralScalePresentAtCurrentBp
        (
            mse.boundaryPoints().size(),
            false
        );

        label nNeutralScaleTransferEdges = 0;
        label nNeutralScaleTransferBpUpdates = 0;

        forAll(splitEdges_, seI)
        {
            const edge& se = splitEdges_[seI];

            scalar edgeScale = scalar(1);
            bool hasEdgeScale = false;

            const label ep0 = se.start();
            const label ep1 = se.end();

            if
            (
                ep0 >= 0
             && ep0 < label(neutralScaleByMeshPoint.size())
             && neutralScaleByMeshPoint[ep0]
                    < scalar(1) - SMALL
            )
            {
                edgeScale = neutralScaleByMeshPoint[ep0];
                hasEdgeScale = true;
            }

            if
            (
                ep1 >= 0
             && ep1 < label(neutralScaleByMeshPoint.size())
             && neutralScaleByMeshPoint[ep1]
                    < scalar(1) - SMALL
            )
            {
                edgeScale =
                    hasEdgeScale
                  ? Foam::min
                    (
                        edgeScale,
                        neutralScaleByMeshPoint[ep1]
                    )
                  : neutralScaleByMeshPoint[ep1];

                hasEdgeScale = true;
            }

            if( !hasEdgeScale )
                continue;

            ++nNeutralScaleTransferEdges;

            // Transfer to whichever split-edge endpoint is a CURRENT
            // boundary point. Usually this is e.start(), but deliberately
            // handle both endpoints so the logic is orientation-independent.
            const label endpoints[2] = {ep0, ep1};

            for( label ei=0; ei<2; ++ei )
            {
                const label meshPtI = endpoints[ei];

                if
                (
                    meshPtI < 0
                 || meshPtI >= label(bp.size())
                )
                    continue;

                const label currBpI = bp[meshPtI];

                if
                (
                    currBpI < 0
                 || currBpI >= label(neutralScaleByCurrentBp.size())
                )
                    continue;

                if( neutralScalePresentAtCurrentBp[currBpI] )
                {
                    neutralScaleByCurrentBp[currBpI] =
                        Foam::min
                        (
                            neutralScaleByCurrentBp[currBpI],
                            edgeScale
                        );
                }
                else
                {
                    neutralScaleByCurrentBp[currBpI] =
                        edgeScale;

                    neutralScalePresentAtCurrentBp[currBpI] =
                        true;
                }

                ++nNeutralScaleTransferBpUpdates;
            }
        }

        label nNeutralScaleCappedFaces = 0;
        label minNeutralFaceLayers = labelMax;
        label maxNeutralFaceLayers = 0;

        // Seeds for topology compatibility propagation.
        // Only faces DIRECTLY reduced by the BL/neutral height rule may seed
        // this ramp. Existing low-layer gap/termination/BLBL faces must not
        // initiate propagation.
        boolList neutralScaleDirectCappedFace
        (
            nLayersAtBndFace_.size(),
            false
        );

        forAll(nLayersAtBndFace_, bfI)
        {
            const label requestedNLayers =
                nLayersAtBndFace_[bfI];

            if( requestedNLayers <= 1 )
                continue;

            if( bfI < 0 || bfI >= label(bFaces.size()) )
                continue;

            const face& bf = bFaces[bfI];

            scalar faceScale = scalar(1);
            bool hasFaceScale = false;

            forAll(bf, fpI)
            {
                const label meshPtI = bf[fpI];

                if
                (
                    meshPtI < 0
                 || meshPtI >= label(bp.size())
                )
                    continue;

                const label currBpI = bp[meshPtI];

                if
                (
                    currBpI < 0
                 || currBpI >= label(neutralScaleByCurrentBp.size())
                 || !neutralScalePresentAtCurrentBp[currBpI]
                )
                    continue;

                const scalar sPt =
                    neutralScaleByCurrentBp[currBpI];

                faceScale =
                    hasFaceScale
                  ? Foam::min(faceScale, sPt)
                  : sPt;

                hasFaceScale = true;
            }

            if( !hasFaceScale )
                continue;

            scalar ratio = globalThicknessRatio_;

            if
            (
                bfI >= 0
             && bfI < label(facePatch.size())
            )
            {
                const label patchI = facePatch[bfI];

                if
                (
                    patchI >= 0
                 && patchI < label(boundaries.size())
                )
                {
                    const word& patchName =
                        boundaries[patchI].patchName();

                    const std::map<word, scalar>::const_iterator rIt =
                        thicknessRatioForPatch_.find(patchName);

                    if( rIt != thicknessRatioForPatch_.end() )
                        ratio = rIt->second;
                }
            }

            const scalar r =
                Foam::max(ratio, scalar(1e-12));

            scalar totalWeight = scalar(0);
            scalar w = scalar(1);

            for
            (
                label li=0;
                li<requestedNLayers;
                ++li
            )
            {
                totalWeight += w;
                w *= r;
            }

            const scalar availableWeight =
                faceScale * totalWeight;

            label faceCap = 1;
            scalar cumulative = scalar(0);
            w = scalar(1);

            const scalar fitTol =
                scalar(100) * SMALL
              * Foam::max(scalar(1), totalWeight);

            for
            (
                label li=1;
                li<=requestedNLayers;
                ++li
            )
            {
                cumulative += w;

                if
                (
                    cumulative
                 <= availableWeight + fitTol
                )
                {
                    faceCap = li;
                }
                else
                {
                    break;
                }

                w *= r;
            }

            faceCap =
                Foam::max
                (
                    label(1),
                    Foam::min(requestedNLayers, faceCap)
                );

            if( faceCap < nLayersAtBndFace_[bfI] )
            {
                nLayersAtBndFace_[bfI] = faceCap;
                neutralScaleDirectCappedFace[bfI] = true;

                ++nNeutralScaleCappedFaces;

                minNeutralFaceLayers =
                    Foam::min
                    (
                        minNeutralFaceLayers,
                        faceCap
                    );

                maxNeutralFaceLayers =
                    Foam::max
                    (
                        maxNeutralFaceLayers,
                        faceCap
                    );
            }
        }

        // BL/neutral FACE layer-count compatibility ramp.
        //
        // The directly height-capped faces may differ sharply from neighboring
        // full-layer faces (for example 2 -> 15). Propagate ONLY outward from
        // those direct seeds and ONLY across faces on the same boundary patch.
        //
        // Enforce:
        //
        //     neighbourLayers <= sourceLayers + 1
        //
        // one topological ring per pass. This creates a deterministic
        // 2->3->4->...->15 recovery without allowing unrelated existing
        // one-layer faces to seed a global reduction.
        boolList neutralRampFrontier
        (
            neutralScaleDirectCappedFace
        );

        boolList neutralRampAdjustedFace
        (
            nLayersAtBndFace_.size(),
            false
        );

        label nNeutralRampAdjustedFaces = 0;
        label nNeutralRampUpdates = 0;
        label nNeutralRampPasses = 0;

        // The required number of rings cannot exceed the maximum layer-count
        // difference. Add a small guard margin rather than using an arbitrary
        // fixed iteration count.
        label neutralRampMaxLayers = 1;

        forAll(nLayersAtBndFace_, bfI)
        {
            neutralRampMaxLayers =
                Foam::max
                (
                    neutralRampMaxLayers,
                    nLayersAtBndFace_[bfI]
                );
        }

        const label neutralRampMaxPasses =
            neutralRampMaxLayers + 2;

        bool neutralRampChanged = true;

        while
        (
            neutralRampChanged
         && nNeutralRampPasses < neutralRampMaxPasses
        )
        {
            neutralRampChanged = false;
            ++nNeutralRampPasses;

            boolList nextFrontier
            (
                nLayersAtBndFace_.size(),
                false
            );

            // Snapshot the layer-count field so one pass advances exactly
            // one boundary-face ring regardless of face ordering.
            const labelList layersBeforePass
            (
                nLayersAtBndFace_
            );

            forAll(neutralRampFrontier, bfI)
            {
                if( !neutralRampFrontier[bfI] )
                    continue;

                if
                (
                    bfI < 0
                 || bfI >= label(bFaces.size())
                 || bfI >= label(facePatch.size())
                )
                    continue;

                const label sourcePatch =
                    facePatch[bfI];

                const label sourceLayers =
                    layersBeforePass[bfI];

                // Permit at most a two-layer change per topological ring.
                // This keeps the neutral transition smooth while reducing the
                // spatial footprint of the compatibility ramp.
                const label neighbourCap =
                    sourceLayers + 2;

                const face& f =
                    bFaces[bfI];

                forAll(f, fpI)
                {
                    const label meshPtI =
                        f[fpI];

                    if
                    (
                        meshPtI < 0
                     || meshPtI >= label(bp.size())
                    )
                        continue;

                    const label bpI =
                        bp[meshPtI];

                    if
                    (
                        bpI < 0
                     || bpI >= label(pointFaces.size())
                    )
                        continue;

                    forAllRow(pointFaces, bpI, pfI)
                    {
                        const label nbfI =
                            pointFaces(bpI, pfI);

                        if
                        (
                            nbfI < 0
                         || nbfI >= label(nLayersAtBndFace_.size())
                         || nbfI >= label(facePatch.size())
                         || nbfI == bfI
                        )
                            continue;

                        // Never cross patch boundaries.
                        if( facePatch[nbfI] != sourcePatch )
                            continue;

                        // Existing virtual-topology treatment has priority.
                        if
                        (
                            nbfI < label(vtFaceRing_.size())
                         && vtFaceRing_[nbfI] >= 0
                        )
                            continue;

                        if
                        (
                            layersBeforePass[nbfI]
                         <= neighbourCap
                        )
                            continue;

                        if
                        (
                            nLayersAtBndFace_[nbfI]
                         > neighbourCap
                        )
                        {
                            nLayersAtBndFace_[nbfI] =
                                neighbourCap;

                            nextFrontier[nbfI] = true;
                            neutralRampChanged = true;
                            ++nNeutralRampUpdates;

                            if
                            (
                                !neutralRampAdjustedFace[nbfI]
                            )
                            {
                                neutralRampAdjustedFace[nbfI] =
                                    true;

                                ++nNeutralRampAdjustedFaces;
                            }
                        }
                    }
                }
            }

            neutralRampFrontier.transfer
            (
                nextFrontier
            );
        }

        label nNeutralRampSeeds = 0;

        forAll(neutralScaleDirectCappedFace, bfI)
        {
            if( neutralScaleDirectCappedFace[bfI] )
                ++nNeutralRampSeeds;
        }

        Info << "BL/neutral FACE layer compatibility ramp:"
             << " maxStep=2"
             << " seedFaces=" << nNeutralRampSeeds
             << " adjustedFaces=" << nNeutralRampAdjustedFaces
             << " updates=" << nNeutralRampUpdates
             << " passes=" << nNeutralRampPasses
             << " maxPasses=" << neutralRampMaxPasses
             << endl;

        Info << "BL/neutral hair-scale transfer:"
             << " sourceEdges=" << nNeutralScaleTransferEdges
             << " bpUpdates=" << nNeutralScaleTransferBpUpdates
             << endl;

        Info << "BL/neutral height-aware FACE layer caps:"
             << " cappedFaces=" << nNeutralScaleCappedFaces;

        if( nNeutralScaleCappedFaces > 0 )
        {
            Info << " minLayers=" << minNeutralFaceLayers
                 << " maxLayers=" << maxNeutralFaceLayers;
        }

        Info << endl;

    }

    //- count the number of vertices for each split edge
    # ifdef USE_OMP
    const label nThreads = omp_get_max_threads();
    # else
    const label nThreads = 1;
    # endif

    boolList cappedEdge(splitEdges_.size(), false);
    # ifdef USE_OMP
    # pragma omp parallel num_threads(nThreads)
    # endif
    {
        //- start counting vertices at each thread
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(splitEdges_, seI)
        {
            const edge& e = splitEdges_[seI];

            //- get the requested number of boundary layers
            label nLayers(1);
            scalar ratio(globalThicknessRatio_);
            scalar thickness(globalMaxThicknessFirstLayer_);
            bool overridenThickness(false);

            const label bpI = bp[e.start()];

            forAllRow(pointFaces, bpI, pfI)
            {
                const label bfI = pointFaces(bpI, pfI);
                const label pos = help::positionOfEdgeInFace(e, bFaces[bfI]);
                if( pos >= 0 )
                    continue;

                const word& patchName =
                    boundaries[facePatch[bfI]].patchName();

                //- overrride the global value with the maximum number of layers
                //- at this edge
                nLayers = Foam::max(nLayers, nLayersAtBndFace_[bfI]);

                //- override with the maximum ratio
                const std::map<word, scalar>::const_iterator rIt =
                    thicknessRatioForPatch_.find(patchName);
                if( rIt != thicknessRatioForPatch_.end() )
                {
                    ratio = rIt->second;
                }

                //- override with the minimum thickness set for this edge
                const std::map<word, scalar>::const_iterator tIt =
                    maxThicknessForPatch_.find(patchName);
                if( tIt != maxThicknessForPatch_.end() )
                {
                    if( overridenThickness )
                    {
                        thickness = Foam::min(thickness, tIt->second);
                    }
                    else
                    {
                        thickness = tIt->second;
                        overridenThickness = true;
                    }
                }
            }

            // Option B split-edge cap: local gap-zone + loser-side check.
            // Checks BOTH edge endpoints (orientation-independent).
            // Uses mesh point labels (stable) + patch names (stable).
            // Only caps where edge touches gap action zone AND loser patch.
            // Minimum effective cap is 1 (ring0 already topology-suppressed).
            if( gapActionPoints_.size() > 0 && gapLoserPatchNames_.size() > 0 )
            {
                // Gap action points are always boundary surface points (e.start()).
                // Interior BL points (e.end()) are never in gapActionPoints_.
                const bool edgeInGap =
                    gapActionPoints_.found(e.start());
                if( edgeInGap )
                {
                    label edgeCap = -1;
                    forAllRow(pointFaces, bpI, pfI)
                    {
                        const label bfI = pointFaces(bpI, pfI);
                        if( bfI < 0 || bfI >= label(facePatch.size()) ) continue;
                        const word& pName =
                            boundaries[facePatch[bfI]].patchName();
                        bool isLoser = false;
                        forAll(gapLoserPatchNames_, pi)
                            if( gapLoserPatchNames_[pi] == pName )
                                { isLoser = true; break; }
                        if( !isLoser ) continue;
                        const label cap = Foam::max(label(1), gapRing1MaxLayers_);
                        edgeCap = (edgeCap < 0) ?
                            cap : Foam::min(edgeCap, cap);
                    }
                    if( edgeCap > 0 && edgeCap < nLayers )
                    {
                        nLayers = edgeCap;
                        cappedEdge[seI] = true;
                    }
                }
            }

            // BL/termination edge cap: reduce layer count at inlet/outlet junctions.
            // Checks BOTH endpoints. No patch-name check needed -- blTerminationEdgePoints_
            // already contains only BL/no-BL transition points by construction.
            // Default cap=3 means disabled (min(3,3)=3, no effect).
            // Start conservative: blTerminationRing1MaxLayers=2 reduces 3->2 only.
            if( blTerminationEdgePoints_.size() > 0
             && blTerminationRing1MaxLayers_ < 3 )
            {
                const bool edgeAtTermination =
                    blTerminationEdgePoints_.found(e.start())
                 || blTerminationEdgePoints_.found(e.end());
                if( edgeAtTermination )
                {
                    const label cap =
                        Foam::max(label(1), blTerminationRing1MaxLayers_);
                    if( cap < nLayers )
                    {
                        nLayers = cap;
                        cappedEdge[seI] = true;
                    }
                }
            }

            //- store the information
            firstLayerThickness[seI] = thickness;

            if( forcedThicknessScaleAtEdge[seI] < scalar(1.0) - SMALL )
                firstLayerThickness[seI] *= forcedThicknessScaleAtEdge[seI];

            thicknessRatio[seI] = ratio;
            nLayersAtEdge[seI] = nLayers;

            if( !specialMode_ )
            {
                nNodesAtEdge[seI] = nLayers + 1;
            }
            else
            {
                nNodesAtEdge[seI] = 3;
            }
        }
    }

    // Split-edge layer-count compatibility smoothing.
    // Only runs when a cap actually reduced at least one split edge.
    // Prevents 3->1 or 3->2 isolated topology cliffs.
    // Does not propagate on geometries with intentional layer differences.
    {
        label nInitiallyCapped = 0;
        forAll(cappedEdge, seI)
            if( cappedEdge[seI] )
                ++nInitiallyCapped;

        if( nInitiallyCapped > 0 )
        {
            const label nSplitEdges = splitEdges_.size();
            bool changed = true;
            label nAdjusted = 0;
            label nPasses = 0;
            const label maxPasses = 20;

            while( changed && nPasses < maxPasses )
            {
                changed = false;
                ++nPasses;

                forAll(splitEdges_, seI)
                {
                    const edge& e = splitEdges_[seI];
                    const label nL = nLayersAtEdge[seI];

                    forAllRow(splitEdgesAtPoint_, e.start(), i)
                    {
                        const label seJ = splitEdgesAtPoint_(e.start(), i);
                        if( seJ == seI || seJ < 0 || seJ >= nSplitEdges )
                            continue;
                        if( nLayersAtEdge[seJ] > nL + 1 )
                        {
                            nLayersAtEdge[seJ] = nL + 1;
                            nNodesAtEdge[seJ] =
                                specialMode_ ? 3 : nLayersAtEdge[seJ] + 1;
                            changed = true;
                            ++nAdjusted;
                        }
                    }

                    forAllRow(splitEdgesAtPoint_, e.end(), i)
                    {
                        const label seJ = splitEdgesAtPoint_(e.end(), i);
                        if( seJ == seI || seJ < 0 || seJ >= nSplitEdges )
                            continue;
                        if( nLayersAtEdge[seJ] > nL + 1 )
                        {
                            nLayersAtEdge[seJ] = nL + 1;
                            nNodesAtEdge[seJ] =
                                specialMode_ ? 3 : nLayersAtEdge[seJ] + 1;
                            changed = true;
                            ++nAdjusted;
                        }
                    }
                }
            }

            Info << "Split-edge layer-count smoothing: "
                 << "initialCapped=" << nInitiallyCapped
                 << " adjusted=" << nAdjusted
                 << " passes=" << nPasses
                 << endl;
        }
    }

    if( Pstream::parRun() )
    {
        //- transfer the information over all processor for edges
        //- at inter-processor boundaries
        const labelLongList& globalEdgeLabel =
            mesh_.addressingData().globalEdgeLabel();
        const VRWGraph& edgeAtProcs = mesh_.addressingData().edgeAtProcs();
        const Map<label>& globalToLocal =
            mesh_.addressingData().globalToLocalEdgeAddressing();
        const DynList<label>& neiProcs = mesh_.addressingData().edgeNeiProcs();
        const edgeList& edges = mesh_.addressingData().edges();
        const VRWGraph& pointEdges = mesh_.addressingData().pointEdges();

        //- exchange point number of layers
        std::map<label, LongList<labelPair> > exchangeNumLayers;
        std::map<label, LongList<labelPair> > exchangeNumNodesAtEdge;
        std::map<label, LongList<labelledScalar> > exchangeThickness;
        std::map<label, LongList<labelledScalar> > exchangeRatio;
        forAll(neiProcs, i)
        {
            exchangeNumNodesAtEdge.insert
            (
                std::make_pair(neiProcs[i], LongList<labelPair>())
            );
            exchangeNumLayers.insert
            (
                std::make_pair(neiProcs[i], LongList<labelPair>())
            );
            exchangeThickness.insert
            (
                std::make_pair(neiProcs[i], LongList<labelledScalar>())
            );
            exchangeRatio.insert
            (
                std::make_pair(neiProcs[i], LongList<labelledScalar>())
            );
        }

        //- exchange the number of layers
        forAll(splitEdges_, seI)
        {
            const edge& se = splitEdges_[seI];

            const label s = se.start();
            label edgeI(-1);
            forAllRow(pointEdges, s, peI)
            {
                const label eI = pointEdges(s, peI);

                if( edges[eI] == se )
                {
                    edgeI = eI;
                    break;
                }
            }

            const label geI = globalEdgeLabel[edgeI];

            if( globalToLocal.found(geI) )
            {
                forAllRow(edgeAtProcs, edgeI, i)
                {
                    const label neiProc = edgeAtProcs(edgeI, i);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    exchangeNumNodesAtEdge[neiProc].append
                    (
                        labelPair(geI, nNodesAtEdge[seI])
                    );
                    exchangeNumLayers[neiProc].append
                    (
                        labelPair(geI, nLayersAtEdge[seI])
                    );
                    exchangeThickness[neiProc].append
                    (
                        labelledScalar(geI, firstLayerThickness[seI])
                    );
                    exchangeRatio[neiProc].append
                    (
                        labelledScalar(geI, thicknessRatio[seI])
                    );
                }
            }
        }

        //- exchange number of nodes at split edge
        LongList<labelPair> receivedNumLayers;
        help::exchangeMap(exchangeNumNodesAtEdge, receivedNumLayers);

        forAll(receivedNumLayers, i)
        {
            const labelPair& lp = receivedNumLayers[i];
            if( !globalToLocal.found(lp.first()) ) continue;
            const label eI = globalToLocal[lp.first()];
            const edge& e = edges[eI];
            label seI(-1);
            forAllRow(splitEdgesAtPoint_, e.start(), i)
            {
                const label seJ = splitEdgesAtPoint_(e.start(), i);
                if( splitEdges_[seJ] == e )
                {
                    seI = seJ;
                    break;
                }
            }
            nNodesAtEdge[seI] = std::max(nNodesAtEdge[seI], lp.second());
        }

        //- exchange number of layers
        receivedNumLayers.clear();
        help::exchangeMap(exchangeNumLayers, receivedNumLayers);

        forAll(receivedNumLayers, i)
        {
            const labelPair& lp = receivedNumLayers[i];
            if( !globalToLocal.found(lp.first()) ) continue;
            const label eI = globalToLocal[lp.first()];
            const edge& e = edges[eI];
            label seI(-1);
            forAllRow(splitEdgesAtPoint_, e.start(), i)
            {
                const label seJ = splitEdgesAtPoint_(e.start(), i);
                if( splitEdges_[seJ] == e )
                {
                    seI = seJ;
                    break;
                }
            }
            nLayersAtEdge[seI] = std::max(nLayersAtEdge[seI], lp.second());
        }

        //- exchange thickness ratio
        LongList<labelledScalar> receivedScalar;
        help::exchangeMap(exchangeRatio, receivedScalar);

        forAll(receivedScalar, i)
        {
            const labelledScalar& ls = receivedScalar[i];
            if( !globalToLocal.found(ls.scalarLabel()) ) continue;
            const label eI = globalToLocal[ls.scalarLabel()];
            const edge& e = edges[eI];
            label seI(-1);
            forAllRow(splitEdgesAtPoint_, e.start(), i)
            {
                const label seJ = splitEdgesAtPoint_(e.start(), i);
                if( splitEdges_[seJ] == e )
                {
                    seI = seJ;
                    break;
                }
            }
            thicknessRatio[seI] = std::max(thicknessRatio[seI], ls.value());
        }

        //- exchange maximum thickness of the first layer
        receivedScalar.clear();
        help::exchangeMap(exchangeThickness, receivedScalar);

        forAll(receivedScalar, i)
        {
            const labelledScalar& ls = receivedScalar[i];
            if( !globalToLocal.found(ls.scalarLabel()) ) continue;
            const label eI = globalToLocal[ls.scalarLabel()];
            const edge& e = edges[eI];
            label seI(-1);
            forAllRow(splitEdgesAtPoint_, e.start(), i)
            {
                const label seJ = splitEdgesAtPoint_(e.start(), i);
                if( splitEdges_[seJ] == e )
                {
                    seI = seJ;
                    break;
                }
            }
            firstLayerThickness[seI] =
                std::min(firstLayerThickness[seI], ls.value());
        }
    }

    //- calculate the number of additional vertices which will be generated
    //- on edges of the mesh
    DynList<label> numPointsAtThread;
    numPointsAtThread.setSize(nThreads);
    numPointsAtThread = 0;

    # ifdef USE_OMP
    # pragma omp parallel for num_threads(nThreads) schedule(static, 1)
    # endif
    forAll(nNodesAtEdge, seI)
    {
        # ifdef USE_OMP
        const label threadI = omp_get_thread_num();
        # else
        const label threadI(0);
        # endif

        numPointsAtThread[threadI] += nNodesAtEdge[seI] - 2;
    }

    //- allocate the space in a graph storing ids of points on a split edge
    newVerticesForSplitEdge_.setSizeAndRowSize(nNodesAtEdge);

    //- calculate the number of points which will be generated
    //- on split edges
    label numPoints = points.size();
    forAll(numPointsAtThread, threadI)
    {
        const label nPts = numPointsAtThread[threadI];
        numPointsAtThread[threadI] = numPoints;
        numPoints += nPts;
    }

    points.setSize(numPoints);

    # ifdef DEBUGLayer
    Info << "Generating split vertices" << endl;
    # endif

    //- generate vertices on split edges
    # ifdef USE_OMP
    # pragma omp parallel num_threads(nThreads)
    # endif
    {
        # ifdef USE_OMP
        const label threadI = omp_get_thread_num();
        # else
        const label threadI(0);
        # endif

        label& nPoints = numPointsAtThread[threadI];

        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(splitEdges_, seI)
        {
            const edge& e = splitEdges_[seI];

            const vector v = e.vec(points);
            const scalar magv = mag(v);

            const label nLayers = newVerticesForSplitEdge_.sizeOfRow(seI) - 1;

            scalar firstThickness = magv / nLayersAtEdge[seI];
            if( thicknessRatio[seI] > (1. + SMALL) )
            {
                firstThickness =
                    magv /
                    (
                        (1 - Foam::pow(thicknessRatio[seI], nLayersAtEdge[seI]))
                        / (1.0 - thicknessRatio[seI])
                    );

                # ifdef DEBUGLayer
                Pout << "Thread " << threadI << endl;
                Pout << "Generating vertices at split edge "
                     << " start point " << points[e.start()]
                     << " end point " << points[e.end()] << endl;
                Pout << "Edge length " << magv << endl;
                Pout << "Thickness of the first layer "
                     << firstThickness << endl;
                # endif
            }

            firstThickness =
                Foam::min
                (
                    Foam::max(firstLayerThickness[seI], SMALL),
                    firstThickness
                );

            if( specialMode_ )
            {
                scalar t = firstThickness;

                for(label i=1;i<nLayersAtEdge[seI]-1;++i)
                    t += firstThickness * Foam::pow(thicknessRatio[seI], i);

                firstThickness = t;
            }

            //- generate vertices for this edge
            newVerticesForSplitEdge_(seI, 0) = e.start();

            // C3: degenerate hair edge (zero-length, BL/BL junction point)
            // Force all intermediate vertices to surface point to produce
            // clean zero-height wedge topology instead of machine-epsilon cells
            if( magv < 1000.0*SMALL )
            {
                for(label pI=1;pI<nLayers;++pI)
                    newVerticesForSplitEdge_(seI, pI) = e.start();
                newVerticesForSplitEdge_(seI, nLayers) = e.start();
                continue;
            }

            scalar param = firstThickness;
            const vector vec = v / (magv + VSMALL);

            for(label pI=1;pI<nLayers;++pI)
            {
                //- generate the new vertex
                const point newP = points[e.start()] + param * vec;

                # ifdef DEBUGLayer
                Pout << "Split edge " << seI << " edge points " << e
                    << " start point " << points[e.start()]
                    << " end point " << points[e.end()]
                    << " param " << param
                    << " new point " << nPoints
                    << " has coordinates " << newP << endl;
                # endif

                param += firstThickness * Foam::pow(thicknessRatio[seI], pI);

                newVerticesForSplitEdge_(seI, pI) = nPoints;
                points[nPoints++] = newP;
            }

            newVerticesForSplitEdge_(seI, nLayers) = e.end();
        }
    }

    if( specialMode_ )
    {
        //- set the number of layers to 2
        forAll(nLayersAtBndFace_, bfI)
            if( nLayersAtBndFace_[bfI] > 1 )
                nLayersAtBndFace_[bfI] = 2;
    }

    # ifdef DEBUGLayer
    for(label procI=0;procI<Pstream::nProcs();++procI)
    {
        if( procI == Pstream::myProcNo() )
        {
            forAll(splitEdges_, seI)
            {
                Pout << "\nSplit edge " << seI << " nodes " << splitEdges_[seI]
                    << " coordinates " << points[splitEdges_[seI][0]]
                    << " " << points[splitEdges_[seI][1]]
                    << " has new points "
                    << newVerticesForSplitEdge_[seI] << endl;

                forAllRow(newVerticesForSplitEdge_, seI, i)
                    Pout << "Point " << i << " on edge ha coordinates "
                         << points[newVerticesForSplitEdge_(seI, i)] << endl;
            }
        }

        returnReduce(1, sumOp<label>());
    }

    Info << "Finished generating vertices at split edges" << endl;
    //::exit(1);
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
