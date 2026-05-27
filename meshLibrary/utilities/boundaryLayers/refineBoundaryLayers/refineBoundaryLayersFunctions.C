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
    detectBoundaryLayers dbl(mPart, is2DMesh_);

    const label nGroups = dbl.nDistinctLayers();
    const labelList& faceInLayer = dbl.faceInLayer();

    //- get the hair edges
    splitEdges_ = dbl.hairEdges();

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
                //- check if the layer is interrupted at this patch
                if(
                    discontinuousLayersForPatch_.find(pName) !=
                    discontinuousLayersForPatch_.end()
                )
                {
                    //- set the number of layers and lock this location
                    nLayersAtPatch[patchNameToIndex[pName]] = it->second;
                    protectedValue[patchNameToIndex[pName]] = true;
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
    if( blblJunctionPoints_.size() > 0 )
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
            forAllRow(ptFaces, bpI, pfI)
            {
                const label bfI = ptFaces(bpI, pfI);
                if( bfI < 0 || bfI >= label(nLayersAtBndFace_.size()) ) continue;
                // Virtual topology takes priority — skip VT-handled faces
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
                    // Virtual topology takes priority — skip VT-handled faces
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

    //- allocate the data from storing parameters applying to a split edge
    LongList<scalar> firstLayerThickness(splitEdges_.size());
    LongList<scalar> thicknessRatio(splitEdges_.size());
    labelLongList nNodesAtEdge(splitEdges_.size());
    labelLongList nLayersAtEdge(splitEdges_.size());

    //- count the number of vertices for each split edge
    # ifdef USE_OMP
    const label nThreads = omp_get_max_threads();
    # else
    const label nThreads = 1;
    # endif

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

            //- store the information
            firstLayerThickness[seI] = thickness;
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
