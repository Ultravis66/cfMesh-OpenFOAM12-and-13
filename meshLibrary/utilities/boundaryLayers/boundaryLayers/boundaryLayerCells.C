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

#include "boundaryLayers.H"
#include "OFstream.H"
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "helperFunctionsPar.H"
#include "demandDrivenData.H"
#include "VRWGraphList.H"

#include "labelledPair.H"
#include "HashSet.H"

#include <map>

//#define DEBUGLayer

# ifdef DEBUGLayer
#include "polyMeshGenAddressing.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayers::createLayerCells(const labelList& patchLabels)
{
    Info << "Starting creating layer cells" << endl;

    const meshSurfaceEngine& mse = surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const edgeList& edges = mse.edges();
    const VRWGraph& faceEdges = mse.faceEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();
    const labelList& boundaryFacePatches = mse.boundaryFacePatches();
    const labelList& faceOwners = mse.faceOwners();
    const labelList& bp = mse.bp();
    const VRWGraph& pointFaces = mse.pointFaces();
    const pointFieldPMG& points = mesh_.points();

    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const VRWGraph& pointPatches = mPart.pointPatches();

    //- mark patches which will be extruded into layer cells
    boolList treatPatches(mesh_.boundaries().size(), false);
    forAll(patchLabels, patchI)
    {
        const label pLabel = patchLabels[patchI];
        forAll(treatPatchesWithPatch_[pLabel], i)
            treatPatches[treatPatchesWithPatch_[pLabel][i]] = true;
    }

    //- create new faces at parallel boundaries
    const Map<label>* otherProcPatchPtr(NULL);
    const Map<label>* otherFaceProcPtr(NULL);

    if( Pstream::parRun() )
    {
        createNewFacesParallel(treatPatches);

        otherProcPatchPtr = &mse.otherEdgeFacePatch();
        otherFaceProcPtr = &mse.otherEdgeFaceAtProc();
    }

    //- create lists for new boundary faces
    VRWGraph newBoundaryFaces;
    labelLongList newBoundaryOwners;
    labelLongList newBoundaryPatches;

    //- create storage for new cells
    VRWGraphList cellsToAdd;

    //- create layer cells and store boundary faces
    const label nOldCells = mesh_.cells().size();

    //- CAP_FACE_USAGE_AUDIT: CSV for all cap-touched faces (local OFstream)
    const bool writeCapFaceAtlas =
        useHardBLBLCapCells_ && !capSideVrtMap_.empty();
    OFstream capFaceAuditOs("blCapCellFaceUsageAtlas.csv");
    if( writeCapFaceAtlas )
        capFaceAuditOs << "bfI,patchName,pKey,nPoints,nCapSide,nCollapsed,"
                       << "fullyCollapsed" << nl;

    //- REDUCED-CELL DRY-RUN ATLAS: classify collapse patterns
    //- Gated by useHardBLBLReducedCells_ (no topology change yet)
    const bool writeReducedDryRun =
        useHardBLBLReducedCells_ && !capSideVrtMap_.empty();
    OFstream capReducedDryRunOs("blCapReducedCellDryRunAtlas.csv");
    if( writeReducedDryRun )
        capReducedDryRunOs
            << "bfI,patchName,pKey,nPoints,nCapSide,nCollapsed,"
            << "cellClass,collapsedMask,"
            << "nCandidateFaces,nValidFaces,nInvalidFaces" << nl;

    forAll(bFaces, bfI)
    {
        if( treatPatches[boundaryFacePatches[bfI]] )
        {
            const face& f = bFaces[bfI];

            const label pKey = patchKey_[boundaryFacePatches[bfI]];
            const label bfIPatch = boundaryFacePatches[bfI];
            //- Patch-aware top vertex resolver: checks capSideVrtMap_
            //- keyed by (meshPointI, patchI) first, then falls back
            //- to findNewNodeLabel. Required because hub+blade share
            //- pKey=0 in the same treatPatchesWithPatch_ group.
            auto capAwareTopLabel = [&](const label baseLabel) -> label
            {
                if( !capSideVrtMap_.empty() )
                {
                    const std::pair<label,label> capKey(baseLabel, bfIPatch);
                    auto it = capSideVrtMap_.find(capKey);
                    if( it != capSideVrtMap_.end() ) return it->second;
                }
                return findNewNodeLabel(baseLabel, pKey);
            };

            // Guard 0: skip faces with missing extruded vertices.
            // findNewNodeLabel returns -1 for transition-zone points
            // that never got a new vertex created. Building a layer
            // cell from such a face creates degenerate zero-volume cells.
            {
                bool missingVertex = false;
                forAll(f, pI)
                {
                    if( findNewNodeLabel(f[pI], pKey) < 0 )
                    {
                        missingVertex = true;
                        break;
                    }
                }
                if( missingVertex )
                {
                    static label nMissing = 0;
                    if( ++nMissing <= 100 )
                        Info << "Skipping layer face: missing extruded"
                             << " vertex bfI=" << bfI
                             << " pKey=" << pKey << endl;
                    continue;
                }
            }

            // Franjo TODO: build proper wedge/pyramid/reduced cells
            // at sharp BL/BL feature curves instead of collapsed prisms.
            // Implementation: appendValidFace canonicalizer + reduced
            // topology for REDUCED_CAP_CELL pattern faces.
            // Gate: useHardBLBLReducedCells_ (default false).

            // CAP_FACE_USAGE_AUDIT: write to local CSV OFstream
            if( writeCapFaceAtlas )
            {
                label nCapSide = 0;
                label nCollapsed = 0;
                forAll(f, pI)
                {
                    const label topLabel = findNewNodeLabel(f[pI], pKey);
                    //- Key by patchI (hub+blade share pKey=0)
                    const std::pair<label,label> capKey(f[pI], boundaryFacePatches[bfI]);
                    if( capSideVrtMap_.find(capKey) != capSideVrtMap_.end() )
                        ++nCapSide;
                    bool collapsed = false;
                    if( topLabel == f[pI] ) collapsed = true;
                    if( !collapsed
                     && topLabel >= 0
                     && topLabel < label(mesh_.points().size())
                     && f[pI] >= 0
                     && f[pI] < label(mesh_.points().size()) )
                    {
                        if( mag(mesh_.points()[topLabel]-mesh_.points()[f[pI]])
                            < scalar(1e-10) ) collapsed = true;
                    }
                    if( collapsed ) ++nCollapsed;
                }
                if( nCapSide > 0 )
                {
                    const label patchI = boundaryFacePatches[bfI];
                    const word pName =
                        (patchI>=0 && patchI<label(patchNames_.size())) ?
                        patchNames_[patchI] : word("?");
                    capFaceAuditOs
                        << bfI << "," << pName << "," << pKey << ","
                        << f.size() << "," << nCapSide << "," << nCollapsed << ","
                        << (nCollapsed==label(f.size()) ? "1" : "0") << nl;

                    if( writeReducedDryRun )
                    {
                        word cellClass = "NORMAL_PRISM";
                        if( nCollapsed == label(f.size()) )
                            cellClass = "INVALID_FULL_COLLAPSE";
                        else if( nCollapsed > 0 )
                            cellClass = "REDUCED_CAP_CELL";

                        //- Build collapsed mask and top face
                        word collapsedMask = "";
                        DynList<label> topFace;
                        forAll(f, pI)
                        {
                            const std::pair<label,label> dryCapKey(f[pI], boundaryFacePatches[bfI]);
                            auto dryIt = capSideVrtMap_.find(dryCapKey);
                            const label topLabel = (dryIt != capSideVrtMap_.end()) ?
                                dryIt->second : findNewNodeLabel(f[pI], pKey);
                            bool coll = (topLabel == f[pI]);
                            if( !coll && topLabel>=0
                             && topLabel<label(mesh_.points().size())
                             && f[pI]>=0
                             && f[pI]<label(mesh_.points().size()) )
                                if( mag(mesh_.points()[topLabel]-mesh_.points()[f[pI]])
                                    < scalar(1e-10) ) coll = true;
                            collapsedMask += (coll ? "1" : "0");
                            topFace.append(coll ? f[pI] : topLabel);
                        }

                        label nCandidateFaces = 0;
                        label nValidFaces = 0;

                        //- Helper lambda: unique-label + area check
                        auto faceValid = [&](const DynList<label>& sf) -> bool
                        {
                            if( sf.size() < 3 ) return false;
                            //- Unique label check
                            forAll(sf, i)
                                for(label j=0;j<i;++j)
                                    if(sf[i]==sf[j]) return false;
                            //- Area check
                            vector areaVec = vector::zero;
                            const point& p0a = mesh_.points()[sf[0]];
                            for(label ai=1;ai<sf.size()-1;++ai)
                                areaVec += (mesh_.points()[sf[ai]]-p0a)
                                         ^ (mesh_.points()[sf[ai+1]]-p0a);
                            return mag(areaVec) > VSMALL;
                        };

                        //- Base face
                        ++nCandidateFaces;
                        {
                            DynList<label> bf2;
                            forAll(f, pI) bf2.append(f[pI]);
                            if( faceValid(bf2) ) ++nValidFaces;
                        }

                        //- Top face
                        ++nCandidateFaces;
                        {
                            DynList<label> tf;
                            forAll(topFace, pI) tf.append(topFace[pI]);
                            if( faceValid(tf) ) ++nValidFaces;
                        }

                        //- Side faces
                        forAll(f, pI)
                        {
                            ++nCandidateFaces;
                            const label b0 = f[pI];
                            const label b1 = f.nextLabel(pI);
                            const label t1 = topFace[(pI+1)%topFace.size()];
                            const label t0 = topFace[pI];
                            DynList<label> sf;
                            sf.append(b0);
                            if( b1!=b0 ) sf.append(b1);
                            if( t1!=b1 && t1!=b0 ) sf.append(t1);
                            if( t0!=t1 && t0!=b0 && t0!=b1 ) sf.append(t0);
                            if( faceValid(sf) ) ++nValidFaces;
                        }

                        const label nInvalidFaces =
                            nCandidateFaces - nValidFaces;
                        capReducedDryRunOs
                            << bfI << "," << pName << "," << pKey << ","
                            << f.size() << "," << nCapSide << "," << nCollapsed << ","
                            << cellClass << "," << collapsedMask << ","
                            << nCandidateFaces << "," << nValidFaces
                            << "," << nInvalidFaces << nl;
                    }
                }
            }

            DynList<DynList<label> > cellFaces;
            DynList<label> newF;

            //- Reduced cap cell builder (Franjo TODO implementation).
            //- Gate: useHardBLBLReducedCells_ (default false).
            bool builtReducedCell = false;
            if( useHardBLBLReducedCells_ && !capSideVrtMap_.empty() )
            {
                label nCapSideRC = 0;
                label nCollapsedRC = 0;
                DynList<label> topFaceRC;
                //- Use outer capAwareTopLabel lambda (patchI-keyed)
                auto isCollapsedTop = [&](const label base, const label top) -> bool
                {
                    if( top == base ) return true;
                    if( top>=0 && top<label(mesh_.points().size())
                     && base>=0 && base<label(mesh_.points().size()) )
                        return mag(mesh_.points()[top]-mesh_.points()[base])
                               < scalar(1e-10);
                    return true;
                };
                forAll(f, pI)
                {
                    const label baseLabel = f[pI];
                    const label topLabel = capAwareTopLabel(baseLabel);
                    const std::pair<label,label> capKey(baseLabel, bfIPatch);
                    if( capSideVrtMap_.find(capKey) != capSideVrtMap_.end() )
                        ++nCapSideRC;
                    const bool coll = isCollapsedTop(baseLabel, topLabel);
                    if( coll ) ++nCollapsedRC;
                    topFaceRC.append(coll ? baseLabel : topLabel);
                }
                if( nCapSideRC > 0 && nCollapsedRC > 0
                 && nCollapsedRC < label(f.size()) )
                {
                    //- Snapshot boundary face lists for rollback
                    const label nBndFacesSnap = newBoundaryFaces.size();
                    const label nBndOwnSnap = newBoundaryOwners.size();
                    const label nBndPatchSnap = newBoundaryPatches.size();

                    //- makeValidFace: simplify duplicates first, then
                    //- validate unique-label count and nonzero area.
                    //- Returns cleaned face in validFace, true if valid.
                    auto makeValidFace =
                    [&](const DynList<label>& cand,
                        DynList<label>& validFace) -> bool
                    {
                        validFace.clear();
                        //- Remove consecutive duplicates
                        forAll(cand, i)
                        {
                            const label lbl = cand[i];
                            if( validFace.size()==0
                             || validFace[validFace.size()-1] != lbl )
                                validFace.append(lbl);
                        }
                        //- Remove last if same as first (closed loop)
                        if( validFace.size() > 1
                         && validFace[validFace.size()-1] == validFace[0] )
                            validFace.setSize(validFace.size()-1);
                        if( validFace.size() < 3 ) return false;
                        //- Unique label check
                        forAll(validFace, i)
                            for(label j=0;j<i;++j)
                                if(validFace[i]==validFace[j]) return false;
                        //- Nonzero area check
                        vector areaVec = vector::zero;
                        const point& p0a = mesh_.points()[validFace[0]];
                        for(label ai=1;ai<validFace.size()-1;++ai)
                            areaVec += (mesh_.points()[validFace[ai]]-p0a)
                                     ^ (mesh_.points()[validFace[ai+1]]-p0a);
                        return mag(areaVec) > VSMALL;
                    };

                    //- Base face (reversed boundary face)
                    DynList<label> baseCand;
                    baseCand.append(f[0]);
                    for(label pI=f.size()-1;pI>0;--pI)
                        baseCand.append(f[pI]);
                    DynList<label> validBaseF;
                    if( makeValidFace(baseCand, validBaseF) )
                        cellFaces.append(validBaseF);

                    //- Top face (collapsed top_i replaced by base_i)
                    DynList<label> validTopF;
                    const bool topOK = makeValidFace(topFaceRC, validTopF);
                    if( topOK )
                    {
                        cellFaces.append(validTopF);
                        newBoundaryFaces.appendList(validTopF);
                        newBoundaryOwners.append(cellsToAdd.size() + nOldCells);
                        newBoundaryPatches.append(boundaryFacePatches[bfI]);
                    }

                    //- Side faces + edge-boundary patch logic
                    bool internalSideMismatch = false;
                    forAll(f, pI)
                    {
                        if( internalSideMismatch ) break;
                        const label b0 = f[pI];
                        const label b1 = f.nextLabel(pI);
                        const label t1 = topFaceRC[(pI+1)%f.size()];
                        const label t0 = topFaceRC[pI];
                        DynList<label> sideCand;
                        sideCand.append(b0);
                        sideCand.append(b1);
                        sideCand.append(t1);
                        sideCand.append(t0);
                        const label edgeI = faceEdges(bfI, pI);
                        bool treatedInternal = false;
                        if( edgeFaces.sizeOfRow(edgeI) == 2 )
                        {
                            label neiFaceT = edgeFaces(edgeI, 0);
                            if( neiFaceT == bfI ) neiFaceT = edgeFaces(edgeI, 1);
                            if( treatPatches[boundaryFacePatches[neiFaceT]] )
                                treatedInternal = true;
                        }
                        DynList<label> validSideF;
                        const bool sideOK = makeValidFace(sideCand, validSideF);
                        if( treatedInternal && (!sideOK || validSideF.size() != 4) )
                        {
                            internalSideMismatch = true;
                            break;
                        }
                        if( sideOK )
                        {
                            cellFaces.append(validSideF);
                            if( edgeFaces.sizeOfRow(edgeI) == 2 )
                            {
                                label neiFace2 = edgeFaces(edgeI, 0);
                                if( neiFace2 == bfI ) neiFace2 = edgeFaces(edgeI, 1);
                                if( !treatPatches[boundaryFacePatches[neiFace2]] )
                                {
                                    newBoundaryFaces.appendList(validSideF);
                                    newBoundaryOwners.append(cellsToAdd.size() + nOldCells);
                                    newBoundaryPatches.append(boundaryFacePatches[neiFace2]);
                                }
                            }
                            else if( edgeFaces.sizeOfRow(edgeI) == 1 )
                            {
                                const Map<label>& otherProcPatch = *otherProcPatchPtr;
                                if( !treatPatches[otherProcPatch[edgeI]] )
                                {
                                    newBoundaryFaces.appendList(validSideF);
                                    newBoundaryOwners.append(cellsToAdd.size() + nOldCells);
                                    newBoundaryPatches.append(otherProcPatch[edgeI]);
                                }
                            }
                        }
                    }
                    static label nReducedAttempt = 0;
                    static label nReducedBuilt = 0;
                    static label nRollbackTopo = 0;
                    static label nRollbackGeom = 0;
                    static label nRollbackInternal = 0;
                    static label nRollbackFewFaces = 0;
                    ++nReducedAttempt;

                    //- Check 1: topological closure
                    //- Every undirected edge must appear exactly twice.
                    auto cellLooksClosed =
                    [&](const DynList<DynList<label> >& cFaces) -> bool
                    {
                        std::map<std::pair<label,label>,label> edgeUse;
                        forAll(cFaces, fi)
                        {
                            const DynList<label>& cf = cFaces[fi];
                            if( cf.size() < 3 ) return false;
                            forAll(cf, i)
                            {
                                const label a = cf[i];
                                const label b = cf[(i+1)%cf.size()];
                                if( a==b ) return false;
                                ++edgeUse[std::make_pair(Foam::min(a,b),Foam::max(a,b))];
                            }
                        }
                        for( auto it=edgeUse.begin(); it!=edgeUse.end(); ++it )
                            if( it->second != 2 ) return false;
                        return true;
                    };

                    //- Check 2: geometric openness
                    //- Sum of oriented face area vectors must be ~0.
                    auto cellGeomClosed =
                    [&](const DynList<DynList<label> >& cFaces) -> bool
                    {
                        vector areaSum = vector::zero;
                        scalar areaMagSum = scalar(0);
                        forAll(cFaces, fi)
                        {
                            const DynList<label>& cf = cFaces[fi];
                            vector fArea = vector::zero;
                            const point& p0f = mesh_.points()[cf[0]];
                            for(label ai=1;ai<cf.size()-1;++ai)
                                fArea += (mesh_.points()[cf[ai]]-p0f)
                                       ^ (mesh_.points()[cf[ai+1]]-p0f);
                            areaSum += fArea;
                            areaMagSum += mag(fArea);
                        }
                        if( areaMagSum < VSMALL ) return false;
                        return mag(areaSum)/areaMagSum < scalar(1e-4);
                    };

                    const bool fewFacesOK = (cellFaces.size() >= 4);
                    const bool topoOK = fewFacesOK
                        && !internalSideMismatch
                        && cellLooksClosed(cellFaces);
                    const bool geomOK = topoOK && cellGeomClosed(cellFaces);

                    if( !fewFacesOK ) ++nRollbackFewFaces;
                    else if( internalSideMismatch ) ++nRollbackInternal;
                    else if( !topoOK ) ++nRollbackTopo;
                    else if( !geomOK ) ++nRollbackGeom;

                    if( fewFacesOK && !internalSideMismatch && topoOK && geomOK )
                    {
                        builtReducedCell = true;
                        ++nReducedBuilt;
                        if
                        (
                            nReducedAttempt == 1
                         || nReducedAttempt % 100 == 0
                         || nReducedBuilt == 1
                         || nReducedBuilt % 100 == 0
                        )
                        {
                            Info << "Reduced cap cells: attempt=" << nReducedAttempt
                                 << " built=" << nReducedBuilt
                                 << " rollbackTopo=" << nRollbackTopo
                                 << " rollbackGeom=" << nRollbackGeom
                                 << " rollbackInternal=" << nRollbackInternal
                                 << " rollbackFewFaces=" << nRollbackFewFaces
                                 << endl;
                        }
                    }
                    else
                    {
                        //- Roll back all speculative additions
                        cellFaces.clear();
                        newBoundaryFaces.setSize(nBndFacesSnap);
                        newBoundaryOwners.setSize(nBndOwnSnap);
                        newBoundaryPatches.setSize(nBndPatchSnap);
                    }
                }
            }

            if( !builtReducedCell )
            {
            //- Normal prism path
            //- store the current boundary face
            newF.clear();
            newF.append(f[0]);
            for(label pI=f.size()-1;pI>0;--pI)
                newF.append(f[pI]);
            cellFaces.append(newF);
            //- create parallel face (patch-aware cap vertex)
            forAll(f, pI)
                newF[pI] = capAwareTopLabel(f[pI]);
            cellFaces.append(newF);
            newBoundaryFaces.appendList(newF);
            newBoundaryOwners.append(cellsToAdd.size() + nOldCells);
            newBoundaryPatches.append(boundaryFacePatches[bfI]);
            //- create quad faces
            newF.setSize(4);
            forAll(f, pI)
            {
                newF[0] = f[pI];
                newF[1] = f.nextLabel(pI);
                newF[2] = capAwareTopLabel(newF[1]);
                newF[3] = capAwareTopLabel(f[pI]);
                cellFaces.append(newF);
                //- check if the face is at the boundary
                //- of the treated partitions
                const label edgeI = faceEdges(bfI, pI);
                if( edgeFaces.sizeOfRow(edgeI) == 2 )
                {
                    label neiFace = edgeFaces(edgeI, 0);
                    if( neiFace == bfI )
                        neiFace = edgeFaces(edgeI, 1);

                    if( !treatPatches[boundaryFacePatches[neiFace]] )
                    {
                        newBoundaryFaces.appendList(newF);
                        newBoundaryOwners.append(cellsToAdd.size() + nOldCells);
                        newBoundaryPatches.append(boundaryFacePatches[neiFace]);
                    }
                }
                else if( edgeFaces.sizeOfRow(edgeI) == 1 )
                {
                    const Map<label>& otherProcPatch = *otherProcPatchPtr;
                    if( !treatPatches[otherProcPatch[edgeI]] )
                    {
                        //- face is a new boundary face
                        newBoundaryFaces.appendList(newF);
                        newBoundaryOwners.append(cellsToAdd.size() + nOldCells);
                        newBoundaryPatches.append(otherProcPatch[edgeI]);
                    }
                }
            }

            # ifdef DEBUGLayer
            Info << "Adding cell " << cellFaces << endl;
            # endif

            } //- end if( !builtReducedCell )

            cellsToAdd.appendGraph(cellFaces);
        }
        else
        {
            # ifdef DEBUGLayer
            Info << "Storing original boundary face "
                << bfI << " into patch " << boundaryFacePatches[bfI] << endl;
            # endif

            newBoundaryFaces.appendList(bFaces[bfI]);
            newBoundaryOwners.append(faceOwners[bfI]);
            newBoundaryPatches.append(boundaryFacePatches[bfI]);
        }
    }

    //- data for parallel execution
    boolList procPoint;
    LongList<DynList<label, 4> > pointProcFaces;
    LongList<labelPair> faceAtPatches;
    if( Pstream::parRun() )
    {
        procPoint.setSize(nPoints_);
        procPoint = false;

        const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();
        const labelList& bPoints = mse.boundaryPoints();

        for
        (
            Map<label>::const_iterator iter=globalToLocal.begin();
            iter!=globalToLocal.end();
            ++iter
        )
        {
            const label bpI = iter();
            procPoint[bPoints[bpI]] = true;
        }
    }

    //- create cells at edges
    forAll(edgeFaces, edgeI)
    {
        //- do not consider edges with no faces attached to it
        if( edgeFaces.sizeOfRow(edgeI) == 0 )
            continue;

        //- cells are generated at the processor with the lowest label
        if(
            (edgeFaces.sizeOfRow(edgeI) == 1) &&
            (otherFaceProcPtr->operator[](edgeI) < Pstream::myProcNo())
        )
            continue;

        //- check if the edge is a feature edge
        const label patchI = boundaryFacePatches[edgeFaces(edgeI, 0)];

        label patchJ(-1);
        if( Pstream::parRun() && otherProcPatchPtr && otherProcPatchPtr->found(edgeI) )
        {
            patchJ = otherProcPatchPtr->operator[](edgeI);
        }
        else if( edgeFaces.sizeOfRow(edgeI) >= 2 )
        {
            patchJ = boundaryFacePatches[edgeFaces(edgeI, 1)];
        }
        else
        {
            continue;  // guard: single-face edge with no parallel info
        }

        if( patchI == patchJ )
            continue;

        //- check if the faces attached to the edge have different keys
        const label pKeyI = patchKey_[patchI];
        const label pKeyJ = patchKey_[patchJ];

        if( pKeyI < 0 || pKeyJ < 0 )
        {
            continue;
            FatalErrorIn
            (
                "void boundaryLayers::createLayerCells(const labelList&)"
            ) << "Patch key is negative at concave edge" << abort(FatalError);
        }

        if( pKeyI == pKeyJ )
            continue;

        const edge& e = edges[edgeI];
        if( otherVrts_.find(e.start()) == otherVrts_.end() )
            continue;
        if( otherVrts_.find(e.end()) == otherVrts_.end() )
            continue;

        //- generate faces of the bnd layer cell
        FixedList<FixedList<label, 4>, 6> cellFaces;
        createNewCellFromEdge(e, pKeyI, pKeyJ, cellFaces);

        //- store boundary faces
        newBoundaryFaces.appendList(cellFaces[1]);
        newBoundaryOwners.append(nOldCells+cellsToAdd.size());
        newBoundaryPatches.append(patchJ);

        newBoundaryFaces.appendList(cellFaces[3]);
        newBoundaryOwners.append(nOldCells+cellsToAdd.size());
        newBoundaryPatches.append(patchI);

        //- check if face 5 is a boundary face or at an inter-processor boundary
        const label bps = bp[e.start()];
        if( bps < 0 ) continue;  // guard: invalid boundary point
        label unusedPatch(-1);
        forAllRow(pointPatches, bps, i)
        {
            const label ptchI = pointPatches(bps, i);

            if( ptchI == patchI )
                continue;
            if( ptchI == patchJ )
                continue;
            if( unusedPatch != -1 )
            {
                unusedPatch = -1;
                break;
            }

            unusedPatch = ptchI;
        }

        if( unusedPatch != -1 && treatedPatch_[unusedPatch] )
        {
            //- add a face in the empty patch in case of 2D layer generation
            newBoundaryFaces.appendList(cellFaces[5]);
            newBoundaryOwners.append(nOldCells+cellsToAdd.size());
            newBoundaryPatches.append(unusedPatch);
        }
        else if( Pstream::parRun() && procPoint[e.start()] )
        {
            //- add a face at inter-pocessor boundary
            pointProcFaces.append(cellFaces[5]);
            faceAtPatches.append(labelPair(patchI, patchJ));
        }

        //- check if face 4 is a boundary face or at an inter-processor boundary
        const label bpe = bp[e.end()];
        if( bpe < 0 ) continue;  // guard: invalid boundary point
        unusedPatch = -1;
        forAllRow(pointPatches, bpe, i)
        {
            const label ptchI = pointPatches(bpe, i);

            if( ptchI == patchI )
                continue;
            if( ptchI == patchJ )
                continue;
            if( unusedPatch != -1 )
            {
                unusedPatch = -1;
                break;
            }

            unusedPatch = ptchI;
        }

        if( unusedPatch != -1 && treatedPatch_[unusedPatch] )
        {
            //- add a face in the empty patch in case of 2D layer generation
            newBoundaryFaces.appendList(cellFaces[4]);
            newBoundaryOwners.append(nOldCells+cellsToAdd.size());
            newBoundaryPatches.append(unusedPatch);
        }
        else if( Pstream::parRun() && procPoint[e.end()] )
        {
            //- add a face at inter-pocessor boundary
            pointProcFaces.append(cellFaces[4]);
            faceAtPatches.append(labelPair(patchI, patchJ));
        }

        # ifdef DEBUGLayer
        Info << "Adding new cell at edge " << cellFaces << endl;
        # endif

        //- append cell to the queue
        cellsToAdd.appendGraph(cellFaces);
    }

    //- create cells for corner nodes
    typedef std::map<std::pair<label, label>, label> mPairToLabelType;
    typedef std::map<label, mPairToLabelType> mPointsType;
    typedef std::map<label, DynList<label, 3> > ppType;

    ppType nodePatches;
    labelHashSet parPoint;

    if( Pstream::parRun() )
    {
        const labelList& bPoints = mse.boundaryPoints();
        const VRWGraph& pProcs = mse.bpAtProcs();
        const labelList& globalPointLabel = mse.globalBoundaryPointLabel();
        const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();

        std::map<label, labelLongList> facesToSend;

        typedef std::map<label, DynList<DynList<label, 8>, 8> > ppfType;

        ppfType parPointFaces;
        ppType parPointPatches;

        forAllConstIter(mPointsType, otherVrts_, iter)
        {
            //- skip points on feature edges
            if( iter->second.size() == 2 )
                continue;

            const label bpI = bp[iter->first];

            if( pProcs.sizeOfRow(bpI) != 0 )
            {
                parPoint.insert(iter->first);

                //- point is at a parallel boundary
                label pMin = pProcs(bpI, 0);
                forAllRow(pProcs, bpI, i)
                {
                    const label prI = pProcs(bpI, i);

                    if( facesToSend.find(prI) == facesToSend.end() )
                        facesToSend.insert
                        (
                            std::make_pair(prI, labelLongList())
                        );

                    if( prI < pMin )
                        pMin = prI;
                }

                if( Pstream::myProcNo() == pMin )
                {
                    DynList<label, 3>& pPatches = parPointPatches[bpI];
                    pPatches.setSize(pointFaces.sizeOfRow(bpI));

                    DynList<DynList<label, 8>, 8>& pFaces = parPointFaces[bpI];
                    pFaces.setSize(pPatches.size());

                    forAllRow(pointFaces, bpI, pfI)
                    {
                        const label bfI = pointFaces(bpI, pfI);
                        const face& bf = bFaces[bfI];

                        pPatches[pfI] = boundaryFacePatches[bfI];

                        DynList<label, 8>& bfCopy = pFaces[pfI];
                        bfCopy.setSize(bf.size());
                        forAll(bf, pI)
                            bfCopy[pI] = globalPointLabel[bp[bf[pI]]];
                    }

                    continue;
                }

                labelLongList& stp = facesToSend[pMin];

                //- send the data to the processor with the lowest label
                //- data is flatenned as follows
                //- 1. the number of faces and global point label
                //- 2. number of points in the face
                //- 3. patch label
                //- 4. global labels of face points
                stp.append(globalPointLabel[bpI]);
                stp.append(pointFaces.sizeOfRow(bpI));
                forAllRow(pointFaces, bpI, pfI)
                {
                    const label bfI = pointFaces(bpI, pfI);
                    const face& bf = bFaces[bfI];

                    stp.append(bf.size());
                    stp.append(boundaryFacePatches[bfI]);
                    forAll(bf, pI)
                        stp.append(globalPointLabel[bp[bf[pI]]]);
                }
            }
        }

        //- exchange data with other processors
        labelLongList receivedData;
        help::exchangeMap(facesToSend, receivedData);

        label counter(0);
        while( counter < receivedData.size() )
        {
            const label gpI_recv = receivedData[counter++];
            const label nFaces = receivedData[counter++];
            if( !globalToLocal.found(gpI_recv) )
            {
                // Skip all face data for this point.
                // Each face packet: face_size + patch_label + face_size points
                for(label fI=0;fI<nFaces;++fI)
                {
                    const label fSize = receivedData[counter++];
                    ++counter; // patch label
                    counter += fSize; // face points
                }
                continue;
            }
            const label bpI = globalToLocal[gpI_recv];
            for(label fI=0;fI<nFaces;++fI)
            {
                DynList<label, 8> f(receivedData[counter++]);
                parPointPatches[bpI].append(receivedData[counter++]);
                forAll(f, pI)
                    f[pI] = receivedData[counter++];
                parPointFaces[bpI].append(f);
            }
        }

        //- sort faces sharing corners at the parallel boundaries
        forAllIter(ppfType, parPointFaces, iter)
        {
            DynList<DynList<label, 8>, 8>& pFaces = iter->second;
            DynList<label, 3>& fPatches = parPointPatches[iter->first];
            const label gpI = globalPointLabel[iter->first];

            for(label i=0;i<pFaces.size();++i)
            {
                const DynList<label, 8>& bf = pFaces[i];
                const label pos = bf.containsAtPosition(gpI);
                const edge e(bf[pos], bf[bf.fcIndex(pos)]);

                for(label j=i+1;j<pFaces.size();++j)
                {
                    const DynList<label, 8>& obf = pFaces[j];
                    if( obf.contains(e.start()) && obf.contains(e.end()) )
                    {
                        DynList<label, 8> add;
                        add = pFaces[i+1];
                        pFaces[i+1] = pFaces[j];
                        pFaces[j] = add;

                        const label pAdd = fPatches[i+1];
                        fPatches[i+1] = fPatches[j];
                        fPatches[j] = pAdd;
                        break;
                    }
                }
            }

            DynList<label, 3> patchIDs;
            forAll(fPatches, fpI)
                patchIDs.appendIfNotIn(fPatches[fpI]);

            nodePatches.insert(std::make_pair(bPoints[iter->first], patchIDs));
        }
    }

    //- sort out point which are not at inter-processor boundaries
    forAllConstIter(mPointsType, otherVrts_, iter)
    {
        if( iter->second.size() == 2 )
            continue;

        if( parPoint.found(iter->first) )
            continue;

        const label bpI = bp[iter->first];

        //- ensure correct orientation
        DynList<label> pFaces(pointFaces.sizeOfRow(bpI));
        forAll(pFaces, fI)
            pFaces[fI] = pointFaces(bpI, fI);

        for(label i=0;i<pFaces.size();++i)
        {
            const face& bf = bFaces[pFaces[i]];
            const edge e = bf.faceEdge(bf.which(iter->first));

            for(label j=i+1;j<pFaces.size();++j)
            {
                const face& obf = bFaces[pFaces[j]];
                if(
                    (obf.which(e.start()) >= 0) &&
                    (obf.which(e.end()) >= 0)
                )
                {
                    const label add = pFaces[i+1];
                    pFaces[i+1] = pFaces[j];
                    pFaces[j] = add;
                    break;
                }
            }
        }

        DynList<label, 3> patchIDs;
        forAll(pFaces, patchI)
        {
            patchIDs.appendIfNotIn(boundaryFacePatches[pFaces[patchI]]);
        }

        nodePatches.insert(std::make_pair(iter->first, patchIDs));
    }

    //- create layer cells for corner nodes
    forAllIter(ppType, nodePatches, iter)
    {
        const DynList<label, 3>& patchIDs = iter->second;
        DynList<label, 3> pKeys;
        forAll(patchIDs, patchI)
        {
            const label pKey = patchKey_[patchIDs[patchI]];

            if( pKey < 0 )
                continue;

            pKeys.appendIfNotIn(pKey);
        }

        if( pKeys.size() != 3 )
            continue;

        # ifdef DEBUGLayer
        Pout << "Creating corner cell at point " << iter->first << endl;
        # endif

        FixedList<FixedList<label, 4>, 6> cellFaces;
        createNewCellFromNode(iter->first, pKeys, cellFaces);

        //- store boundary faces
        newBoundaryFaces.appendList(cellFaces[1]);
        newBoundaryOwners.append(nOldCells+cellsToAdd.size());
        newBoundaryPatches.append(patchIDs[0]);

        newBoundaryFaces.appendList(cellFaces[3]);
        newBoundaryOwners.append(nOldCells+cellsToAdd.size());
        newBoundaryPatches.append(patchIDs[1]);

        newBoundaryFaces.appendList(cellFaces[5]);
        newBoundaryOwners.append(nOldCells+cellsToAdd.size());
        newBoundaryPatches.append(patchIDs[2]);

        if( Pstream::parRun() )
        {
            if( procPoint[iter->first] )
            {
                pointProcFaces.append(cellFaces[0]);
                faceAtPatches.append(labelPair(patchIDs[1], patchIDs[2]));

                pointProcFaces.append(cellFaces[2]);
                faceAtPatches.append(labelPair(patchIDs[0], patchIDs[2]));

                pointProcFaces.append(cellFaces[4]);
                faceAtPatches.append(labelPair(patchIDs[0], patchIDs[1]));
            }
        }

        # ifdef DEBUGLayer
        Info << "Adding corner cell " << cellFaces << endl;
        # endif

        //- append cell to the queue
        cellsToAdd.appendGraph(cellFaces);
    }

    if( Pstream::parRun() )
    {
        //- create faces at parallel boundaries created from
        //- points at parallel boundaries
        createNewFacesFromPointsParallel
        (
            pointProcFaces,
            faceAtPatches
        );
    }

    //- create mesh modifier
    polyMeshGenModifier meshModifier(mesh_);

    meshModifier.addCells(cellsToAdd);

    cellsToAdd.clear();
    meshModifier.reorderBoundaryFaces();
    meshModifier.replaceBoundary
    (
        patchNames_,
        newBoundaryFaces,
        newBoundaryOwners,
        newBoundaryPatches
    );

    PtrList<boundaryPatch>& boundaries = meshModifier.boundariesAccess();
    forAll(boundaries, patchI)
        boundaries[patchI].patchType() = patchTypes_[patchI];

    //- delete meshSurfaceEngine
    this->clearOut();

    Info << "Finished creating layer cells" << endl;
}

void boundaryLayers::createNewFacesFromPointsParallel
(
    const LongList<DynList<label, 4> >& faceCandidates,
    const LongList<labelPair>& candidatePatches
)
{
    const meshSurfaceEngine& mse = this->surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    const labelList& bp = mse.bp();
    const VRWGraph& bpAtProcs = mse.bpAtProcs();
    const labelList& globalPointLabel = mse.globalBoundaryPointLabel();
    const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();

    labelList otherFaceProc(faceCandidates.size(), -1);
    //- some faces may appear more than once
    //- such faces are ordinary internal faces
    VRWGraph pointFaceCandidates(nPoints_);
    forAll(faceCandidates, fI)
    {
        forAll(faceCandidates[fI], pI)
            pointFaceCandidates.append(faceCandidates[fI][pI], fI);
    }

    boolList duplicateFace(faceCandidates.size(), false);
    List<labelledPair> pointOfOrigin(faceCandidates.size());
    std::map<labelledPair, label> pointOfOriginToFaceLabel;
    forAll(faceCandidates, fI)
    {
        const DynList<label, 4>& f = faceCandidates[fI];

        const label pointI = f[0];

        const labelledPair lp
        (
            globalPointLabel[bp[pointI]],
            Pair<label>
            (
                patchKey_[candidatePatches[fI][0]],
                patchKey_[candidatePatches[fI][1]]
            )
        );

        if(
            pointOfOriginToFaceLabel.find(lp) != pointOfOriginToFaceLabel.end()
        )
        {
            duplicateFace[fI] = true;
            pointOfOrigin[fI] = lp;
            duplicateFace[pointOfOriginToFaceLabel[lp]] = true;
            continue;
        }

        pointOfOrigin[fI] = lp;

        pointOfOriginToFaceLabel.insert(std::make_pair(lp, fI));
    }

    //- find the processor patch for each processor boundary face
    //- the key of the algorithm is the point from which the face was created
    //- by sending the point label and the associated patches, it will be
    //- possible to find the other processor containing that face
    std::map<label, LongList<labelledPair> > exchangeData;
    const DynList<label>& neiProcs = mse.bpNeiProcs();
    forAll(neiProcs, procI)
    {
        const label neiProcI = neiProcs[procI];

        if( neiProcI == Pstream::myProcNo() )
            continue;

        if( exchangeData.find(neiProcI) == exchangeData.end() )
            exchangeData.insert
            (
                std::make_pair(neiProcI, LongList<labelledPair>())
            );
    }

    forAll(faceCandidates, fI)
    {
        if( duplicateFace[fI] )
            continue;

        const label bpI = bp[faceCandidates[fI][0]];

        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProcNo = bpAtProcs(bpI, procI);
            if( neiProcNo == Pstream::myProcNo() )
                continue;

            LongList<labelledPair>& dataToSend = exchangeData[neiProcNo];
            dataToSend.append(pointOfOrigin[fI]);
        }
    }

    //- exchange the data with other processors
    std::map<label, List<labelledPair> > receivedMap;
    help::exchangeMap(exchangeData, receivedMap);
    exchangeData.clear();

    for
    (
        std::map<label, List<labelledPair> >::const_iterator
        iter=receivedMap.begin();
        iter!=receivedMap.end();
        ++iter
    )
    {
        const List<labelledPair>& receivedData = iter->second;

        forAll(receivedData, i)
        {
            const labelledPair& lpp = receivedData[i];
            const label gpI = lpp.pairLabel();
            if( !globalToLocal.found(gpI) ) continue;
            const label pointI = bPoints[globalToLocal[gpI]];
            const labelPair& lp = lpp.pair();

            forAllRow(pointFaceCandidates, pointI, i)
            {
                const label fI = pointFaceCandidates(pointI, i);
                const DynList<label, 4>& f = faceCandidates[fI];

                const labelPair pk
                (
                    patchKey_[candidatePatches[fI][0]],
                    patchKey_[candidatePatches[fI][1]]
                );

                const labelPair rpk
                (
                    patchKey_[candidatePatches[fI][1]],
                    patchKey_[candidatePatches[fI][0]]
                );

                if(
                    (f[0] == pointI) && ((pk == lp) || (rpk == lp))
                )
                {
                    //- found the processor containing other face
                    auto fIt = pointOfOriginToFaceLabel.find(lpp);
                    if( fIt != pointOfOriginToFaceLabel.end() )
                        otherFaceProc[fIt->second] = iter->first;
                }
            }
        }
    }
    receivedMap.clear();

    //- sort the points in ascending order
    //- this ensures the correct order of faces at the processor boundaries
    sort(pointOfOrigin);

    Map<label> otherProcToProcPatch;
    forAll(mesh_.procBoundaries(), patchI)
    {
        const processorBoundaryPatch& wp = mesh_.procBoundaries()[patchI];
        otherProcToProcPatch.insert(wp.neiProcNo(), patchI);
    }

    //- store processor faces
    VRWGraph newProcFaces;
    labelLongList newProc;

    forAll(pointOfOrigin, i)
    {
        const label fI = pointOfOriginToFaceLabel[pointOfOrigin[i]];

        if( duplicateFace[fI] || (otherFaceProc[fI] == -1) )
            continue;

        if( !otherProcToProcPatch.found(otherFaceProc[fI]) )
        {
            otherProcToProcPatch.insert
            (
                otherFaceProc[fI],
                polyMeshGenModifier(mesh_).addProcessorPatch
                (
                    otherFaceProc[fI]
                )
            );
        }

        newProcFaces.appendList(faceCandidates[fI]);
        newProc.append(otherProcToProcPatch[otherFaceProc[fI]]);
    }

    polyMeshGenModifier(mesh_).addProcessorFaces(newProcFaces, newProc);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
