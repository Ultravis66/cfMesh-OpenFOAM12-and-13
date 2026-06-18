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

#include "demandDrivenData.H"
#include "meshSurfaceEngineModifier.H"
#include "meshSurfaceMapper.H"
#include "meshSurfacePartitioner.H"
#include "meshOctree.H"
#include "triSurf.H"
#include "helperFunctionsPar.H"
#include "helperFunctions.H"
#include "OFstream.H"

#include <map>

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGMapping

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Private member functions

void meshSurfaceMapper::selectNodesAtParallelBnd(const labelLongList& selNodes)
{
    if( !Pstream::parRun() )
        return;

    std::map<label, labelLongList> exchangeData;
    const DynList<label>& neiProcs = surfaceEngine_.bpNeiProcs();
    forAll(neiProcs, i)
        exchangeData.insert(std::make_pair(neiProcs[i], labelLongList()));

    const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
    const labelList& globalPointLabel =
        surfaceEngine_.globalBoundaryPointLabel();
    const Map<label>& globalToLocal =
        surfaceEngine_.globalToLocalBndPointAddressing();

    boolList selectedNode(bpAtProcs.size(), false);

    forAll(selNodes, i)
    {
        const label bpI = selNodes[i];

        selectedNode[bpI] = true;

        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            exchangeData[neiProc].append(globalPointLabel[bpI]);
        }
    }

    //- exchange data
    labelLongList receivedData;
    help::exchangeMap(exchangeData, receivedData);

    forAll(receivedData, i)
    {
            if( !globalToLocal.found(receivedData[i]) ) { continue; }
        if( !selectedNode[globalToLocal[receivedData[i]]] )
        {
            selectedNode[globalToLocal[receivedData[i]]] = true;
            const_cast<labelLongList&>(selNodes).append
            (
                globalToLocal[receivedData[i]]
            );
        }
    }
}

void meshSurfaceMapper::mapToSmallestDistance(LongList<parMapperHelper>& parN)
{
    if( !Pstream::parRun() )
        return;

    std::map<label, LongList<parMapperHelper> > exchangeData;
    const DynList<label>& neiProcs = surfaceEngine_.bpNeiProcs();
    forAll(neiProcs, i)
        exchangeData.insert
        (
            std::make_pair(neiProcs[i], LongList<parMapperHelper>())
        );

    const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
    const labelList& globalPointLabel =
        surfaceEngine_.globalBoundaryPointLabel();
    const Map<label>& globalToLocal =
        surfaceEngine_.globalToLocalBndPointAddressing();

    Map<label> bpToList(parN.size());

    forAll(parN, i)
    {
        const label bpI = parN[i].globalLabel();
        bpToList.insert(bpI, i);

        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            exchangeData[neiProc].append
            (
                parMapperHelper
                (
                    parN[i].coordinates(),
                    parN[i].movingDistance(),
                    globalPointLabel[bpI],
                    parN[i].pointPatch()
                )
            );
        }
    }

    //- exchange data
    LongList<parMapperHelper> receivedData;
    help::exchangeMap(exchangeData, receivedData);

    //- select the point with the smallest moving distance
    meshSurfaceEngineModifier surfModifier(surfaceEngine_);
    forAll(receivedData, i)
    {
        const parMapperHelper& ph = receivedData[i];

        if( !globalToLocal.found(ph.globalLabel()) ) continue;
        const label bpI = globalToLocal[ph.globalLabel()];

        parMapperHelper& phOrig = parN[bpToList[bpI]];
        // Select the candidate with the SMALLEST moving distance --
        // was previously reversed (selected largest), causing wrong
        // projection choices at parallel processor boundaries.
        if( ph.movingDistance() < phOrig.movingDistance() )
        {
            surfModifier.moveBoundaryVertexNoUpdate(bpI, ph.coordinates());
            phOrig = ph;
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceMapper::mapNodeToPatch(const label bpI, const label patchI)
{
    label patch, nt;
    point mapPoint;
    scalar dSq;

    const pointFieldPMG& points = surfaceEngine_.points();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const point p = points[bPoints[bpI]];

    if( patchI < 0 )
    {
        meshOctree_.findNearestSurfacePoint(mapPoint, dSq, nt, patch, p);
    }
    else
    {
        meshOctree_.findNearestSurfacePointInRegion
        (
            mapPoint,
            dSq,
            nt,
            patchI,
            p
        );
    }

    meshSurfaceEngineModifier surfModifier(surfaceEngine_);
    surfModifier.moveBoundaryVertex(bpI, mapPoint);
}

void meshSurfaceMapper::mapVerticesOntoSurface()
{
    Info << "Mapping vertices onto surface" << endl;

    labelLongList nodesToMap(surfaceEngine_.boundaryPoints().size());
    forAll(nodesToMap, i)
        nodesToMap[i] = i;

    mapVerticesOntoSurface(nodesToMap);

    Info << "Finished mapping vertices onto surface" << endl;
}

void meshSurfaceMapper::mapVerticesOntoSurface(const labelLongList& nodesToMap)
{
    const labelList& boundaryPoints = surfaceEngine_.boundaryPoints();
    const pointFieldPMG& points = surfaceEngine_.points();

    const VRWGraph* bpAtProcsPtr(NULL);
    if( Pstream::parRun() )
        bpAtProcsPtr = &surfaceEngine_.bpAtProcs();

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);
    LongList<parMapperHelper> parallelBndNodes;

    // Store old positions for validity-check revert
    pointField oldPositions(nodesToMap.size());
    forAll(nodesToMap, i)
        oldPositions[i] = points[boundaryPoints[nodesToMap[i]]];

    // Build filtered node list excluding BL/no-BL protected points.
    // Protected points (boundary-point indices) must not be moved by
    // generic nearest-surface projection -- they are constrained to
    // their feature curve by the BL/no-BL transition system.
    // If protectedPoints_ is empty this is a no-op (default behaviour).
    labelLongList filteredNodes;
    if( protectedPoints_.empty() )
    {
        filteredNodes = nodesToMap;
    }
    else
    {
        filteredNodes.setSize(nodesToMap.size());
        label nFiltered = 0;
        forAll(nodesToMap, i)
        {
            if( !protectedPoints_.found(nodesToMap[i]) )
                filteredNodes[nFiltered++] = nodesToMap[i];
        }
        filteredNodes.setSize(nFiltered);
        Info << "mapVerticesOntoSurface: excluded "
             << (nodesToMap.size() - nFiltered)
             << " protected BL/no-BL interface points" << endl;
    }

    // Patch-constrained projection: multi-patch points (npp > 1) get
    // their unconstrained nearest-surface result validated against their
    // patch membership set. If the returned patch is not in the set,
    // re-project into each valid patch region and keep the closest hit.
    // Single-patch interior points use the fast unconstrained path unchanged.
    const VRWGraph& ppLocal = meshPartitioner().pointPatches();

    # ifdef USE_OMP
    const label size = filteredNodes.size();
    # pragma omp parallel for if( size > 1000 ) shared(parallelBndNodes) \
    schedule(dynamic, Foam::max(1, size / (3 * omp_get_max_threads())))
    # endif
    forAll(filteredNodes, i)
    {
        const label bpI = filteredNodes[i];

        # ifdef DEBUGMapping
        Info << nl << "Mapping vertex " << bpI << " with coordinates "
            << points[boundaryPoints[bpI]] << endl;
        # endif

        label patch, nt;
        point mapPoint;
        scalar dSq;

        meshOctree_.findNearestSurfacePoint
        (
            mapPoint,
            dSq,
            nt,
            patch,
            points[boundaryPoints[bpI]]
        );

        // Patch-constraint validation for multi-patch points only.
        const label npp = ppLocal.sizeOfRow(bpI);
        if( npp > 1 )
        {
            bool patchValid = false;
            for( label ppI = 0; ppI < npp; ++ppI )
            {
                if( ppLocal(bpI, ppI) == patch )
                { patchValid = true; break; }
            }
            if( !patchValid )
            {
                scalar bestDsq   = GREAT;
                point  bestPt    = mapPoint;
                label  bestPatch = patch;
                label  bestNt    = nt;
                for( label ppI = 0; ppI < npp; ++ppI )
                {
                    point  rPt;
                    scalar rDsq;
                    label  rNt;
                    const label region = ppLocal(bpI, ppI);
                    meshOctree_.findNearestSurfacePointInRegion
                    (
                        rPt,
                        rDsq,
                        rNt,
                        region,
                        points[boundaryPoints[bpI]]
                    );
                    if( (rDsq < bestDsq - SMALL) ||
                        (Foam::mag(rDsq - bestDsq) <= SMALL && region < bestPatch) )
                    {
                        bestDsq   = rDsq;
                        bestPt    = rPt;
                        bestPatch = region;
                        bestNt    = rNt;
                    }
                }
                if( bestDsq < GREAT/10 )
                {
                    mapPoint = bestPt;
                    dSq      = bestDsq;
                    patch    = bestPatch;
                    nt       = bestNt;
                }
            }
        }

        surfaceModifier.moveBoundaryVertexNoUpdate(bpI, mapPoint);

        if( bpAtProcsPtr && bpAtProcsPtr->sizeOfRow(bpI) )
        {
            # ifdef USE_OMP
            # pragma omp critical
            # endif
            parallelBndNodes.append
            (
                parMapperHelper
                (
                    mapPoint,
                    dSq,
                    bpI,
                    patch
                )
            );
        }

        # ifdef DEBUGMapping
        Info << "Mapped point " << points[boundaryPoints[bpI]] << endl;
        # endif
    }

    //- make sure that the points are at the nearest location on the surface
    mapToSmallestDistance(parallelBndNodes);

    // Validity check: serial pass after all moves complete
    // surfaceEngine_ data accessed outside OMP - thread safe
    {
        const VRWGraph& pFaces = surfaceEngine_.pointFaces();
        const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();
        const labelList& faceOwners = surfaceEngine_.faceOwners();
        const pointFieldPMG& pts = surfaceEngine_.points();
        const cellListPMG& cells = surfaceEngine_.mesh().cells();
        const faceListPMG& allFaces = surfaceEngine_.mesh().faces();
        label nInvalid = 0;
        rejectedBpI_.clear();
        // Only evaluate points that were actually moved.
        // nodesToMap includes protected/unmoved points which can have
        // pre-existing invalid faces -- flagging them seeds repairRejectedPoints()
        // with a large 2-ring neighborhood near corners/junctions.
        forAll(filteredNodes, i)
        {
            const label bpI = filteredNodes[i];
            bool validMove = true;
            forAllRow(pFaces, bpI, pfI)
            {
                const label bfI = pFaces(bpI, pfI);
                const face& f = bFaces[bfI];
                point fc = point::zero;
                forAll(f, fpI) fc += pts[f[fpI]];
                fc /= scalar(f.size());
                vector fn = vector::zero;
                const point& p0 = pts[f[0]];
                for(label fpI=1; fpI<f.size()-1; ++fpI)
                    fn += (pts[f[fpI]]-p0)^(pts[f[fpI+1]]-p0);
                const label cellI = faceOwners[bfI];
                point cc = point::zero;
                const cell& cll = cells[cellI];
                forAll(cll, cfI)
                {
                    const face& cf = allFaces[cll[cfI]];
                    point cfc = point::zero;
                    forAll(cf, cpI) cfc += pts[cf[cpI]];
                    cc += cfc / scalar(cf.size());
                }
                cc /= scalar(cll.size());
                // Relative threshold: reject only if pyramid height is
                // meaningfully negative relative to local face size.
                // Fixed SMALL threshold rejects valid moves at small cells
                // near patch boundaries. Scale by face normal magnitude
                // (proportional to face area) for cell-size independence.
                const scalar faceScale = Foam::mag(fn) * scalar(1e-6);
                if( (fn & (fc - cc)) <= -faceScale )
                { validMove = false; break; }
            }
            if( !validMove )
            {
                rejectedBpI_.append(bpI);
                // EXPERIMENTAL: binary search backtracking.
                // Gated off until proposedMoveIsValid includes
                // skewness/non-ortho comparison against old position.
                // Enabling this improves visual spike removal but
                // worsens skewness/non-ortho metrics until validator
                // is strengthened. See feature/bisection-quality branch.
                const bool useBisectionBacktracking = false;
                if( useBisectionBacktracking )
                {
                    const point snapPos = pts[boundaryPoints[bpI]];
                    const point& origPos = oldPositions[i];
                    point bestPos = origPos;
                    scalar alpha = 0.5;
                    for( label bisI = 0; bisI < 6; ++bisI )
                    {
                        const point candidate =
                            origPos + alpha*(snapPos - origPos);
                        if( proposedMoveIsValid(bpI, candidate, origPos,
                                scalar(4.0)*magSqr(snapPos - origPos)) )
                        {
                            bestPos = candidate;
                            break;
                        }
                        alpha *= 0.5;
                    }
                    surfaceModifier.moveBoundaryVertexNoUpdate(bpI, bestPos);
                }
                else
                {
                    surfaceModifier.moveBoundaryVertexNoUpdate
                        (bpI, oldPositions[i]);
                }
                ++nInvalid;
            }
        }
        if( nInvalid > 0 )
            Info << "[ValidityCheck] reverted " << nInvalid
                 << " invalid surface moves" << endl;

        // Write validity-rejected points to VTK (targeted diagnostic)
        // Gate: only write when processing full boundary (nodesToMap
        // size matches all boundary points = first global projection).
        if( !rejectedBpI_.empty() &&
            nodesToMap.size() ==
                label(surfaceEngine_.boundaryPoints().size()) )
        {
            static label rejCallCount = 0;
            ++rejCallCount;
            const labelList& bPtsR = surfaceEngine_.boundaryPoints();
            const pointFieldPMG& ptsR = surfaceEngine_.points();
            const VRWGraph& ptPtsR = surfaceEngine_.pointPoints();
            labelHashSet rejSet;
            forAll(rejectedBpI_, k) rejSet.insert(rejectedBpI_[k]);
            labelHashSet ringSet(rejSet);
            forAll(rejectedBpI_, k)
                forAllRow(ptPtsR, rejectedBpI_[k], nI)
                    ringSet.insert(ptPtsR(rejectedBpI_[k], nI));
            auto writeVTKCloud = [&](const word& name,
                const DynamicList<label>& bpList)
            {
                fileName fName(name + word("_call")
                    + Foam::name(rejCallCount) + word(".vtk"));
                OFstream os(fName);
                os<<"# vtk DataFile Version 2.0\n"<<name
                  <<"\nASCII\nDATASET POLYDATA\n";
                os<<"POINTS "<<bpList.size()<<" float\n";
                forAll(bpList, k)
                {
                    const point& p = ptsR[bPtsR[bpList[k]]];
                    os<<p.x()<<" "<<p.y()<<" "<<p.z()<<"\n";
                }
                os<<"VERTICES "<<bpList.size()
                  <<" "<<2*bpList.size()<<"\n";
                for(label k=0;k<bpList.size();++k)
                    os<<"1 "<<k<<"\n";
                Info<<"[Diag] "<<bpList.size()
                    <<" pts -> "<<fName<<endl;
            };
            writeVTKCloud("validityRejectedPoints", rejectedBpI_);
            // Write ring as labelList
            DynamicList<label> ringList;
            forAllConstIter(labelHashSet, ringSet, it)
                ringList.append(it.key());
            writeVTKCloud("rejectedNeighbourhood", ringList);
        }

    }

    //- re-calculate face normals, point normals, etc.
    surfaceModifier.updateGeometry(nodesToMap);
}

void meshSurfaceMapper::mapVerticesOntoSurfacePatches()
{
    Info << "Mapping vertices with respect to surface patches" << endl;

    labelLongList nodesToMap(surfaceEngine_.boundaryPoints().size());
    forAll(nodesToMap, i)
        nodesToMap[i] = i;

    mapVerticesOntoSurfacePatches(nodesToMap);

    Info << "Finished mapping vertices with respect to surface patches" << endl;
}

void meshSurfaceMapper::mapVerticesOntoSurfacePatches
(
    const labelLongList& nodesToMap
)
{
    const meshSurfacePartitioner& mPart = meshPartitioner();
    const labelHashSet& cornerPoints = mPart.corners();
    const labelHashSet& edgePoints = mPart.edgePoints();
    const VRWGraph& pointPatches = mPart.pointPatches();

    boolList treatedPoint(surfaceEngine_.boundaryPoints().size(), false);

    // ---------------------------------------------------------------
    // Phase 1: boundary point topology classification
    // Diagnostics only -- do not change projection behavior here.
    //   0 = single-patch surface point
    //   1 = two-patch feature-edge point
    //   2 = multi-patch corner/junction point
    //   3 = BL/no-BL transition point
    //   4 = non-manifold / unclassified point
    // ---------------------------------------------------------------
    const label CLS_SINGLE   = 0;
    const label CLS_TWOPATCH = 1;
    const label CLS_CORNER   = 2;
    const label CLS_BLNOBL     = 3;
    const label CLS_NONMANIF   = 4;
    const label CLS_BL_NEUTRAL = 5;

    labelList pointClass(surfaceEngine_.boundaryPoints().size(), CLS_SINGLE);

    label nSingle = 0;
    label nTwoPatch = 0;
    label nCorner = 0;
    label nBlNoBl = 0;
    label nNonManif = 0;
    label nBlNeutral = 0;

    forAll(nodesToMap, i)
    {
        const label bpI = nodesToMap[i];
        const label nPatches = pointPatches.sizeOfRow(bpI);

        if( protectedPoints_.found(bpI) )
        {
            pointClass[bpI] = CLS_BLNOBL;
            ++nBlNoBl;
        }
        else if( nPatches == 0 )
        {
            pointClass[bpI] = CLS_NONMANIF;
            ++nNonManif;
        }
        else if( nPatches == 1 )
        {
            pointClass[bpI] = CLS_SINGLE;
            ++nSingle;
        }
        else if( nPatches == 2 )
        {
            if( !blNeutralPoints_.empty() && blNeutralPoints_.found(bpI) )
            {
                pointClass[bpI] = CLS_BL_NEUTRAL;
                ++nBlNeutral;
            }
            else
            {
                pointClass[bpI] = CLS_TWOPATCH;
                ++nTwoPatch;
            }
        }
        else
        {
            pointClass[bpI] = CLS_CORNER;
            ++nCorner;
        }
    }

    Info << "BoundaryPointClassifier: "
         << nSingle    << " single-patch, "
         << nTwoPatch  << " two-patch-edge, "
         << nCorner    << " multi-patch-corner, "
         << nBlNoBl    << " BL/no-BL transition, "
         << nBlNeutral << " BL/neutral-edge, "
         << nNonManif  << " non-manifold" << endl;

    // Phase 2: VTK diagnostic output for topology classes.
    // Diagnostics only -- no movement or projection behavior changes.
    if( nCorner > 0 || nTwoPatch > 0 || nBlNoBl > 0 || nNonManif > 0 || nBlNeutral > 0 )
    {
        static label callCount = 0;
        ++callCount;

        const labelList& bPts = surfaceEngine_.boundaryPoints();
        const pointFieldPMG& allPts = surfaceEngine_.points();

        auto writePointClassVTK =
        [&](const word& name, const label cls, const label count)
        {
            if( count <= 0 ) return;
            fileName fName(name + word("_call") + Foam::name(callCount) + word(".vtk"));
            OFstream os(fName);
            os << "# vtk DataFile Version 2.0\n";
            os << name << "\n";
            os << "ASCII\n";
            os << "DATASET POLYDATA\n";
            os << "POINTS " << count << " float\n";
            forAll(nodesToMap, i)
            {
                const label bpI = nodesToMap[i];
                if( pointClass[bpI] != cls ) continue;
                const point& p = allPts[bPts[bpI]];
                os << p.x() << " " << p.y() << " " << p.z() << "\n";
            }
            os << "VERTICES " << count << " " << 2*count << "\n";
            for(label k = 0; k < count; ++k)
                os << "1 " << k << "\n";
            Info << "Wrote " << count << " " << name
                 << " points to " << fName << endl;
        };

        writePointClassVTK("twoPatchEdgePoints", CLS_TWOPATCH, nTwoPatch);
        writePointClassVTK("multiPatchCornerPoints", CLS_CORNER, nCorner);
        writePointClassVTK("blNoBlTransitionPoints", CLS_BLNOBL, nBlNoBl);
        writePointClassVTK("blNeutralEdgePoints", CLS_BL_NEUTRAL, nBlNeutral);
        writePointClassVTK("nonManifoldPoints", CLS_NONMANIF, nNonManif);
    }

    // ---------------------------------------------------------------
    // Phase 2: Feature curve extraction
    // Build list of feature edge segments -- boundary edges where
    // adjacent faces belong to different patches.
    // ---------------------------------------------------------------
    const edgeList& meshEdges = surfaceEngine_.edges();
    const VRWGraph& edgeFaces = surfaceEngine_.edgeFaces();
    const labelList& facePatch = surfaceEngine_.boundaryFacePatches();
    const pointFieldPMG& allPoints = surfaceEngine_.points();

    // Compute patch diversity immediately after facePatch is available
    label nPatchedFaces2 = 0;
    label nUniquePatchVals = 0;
    {
        labelHashSet seen;
        forAll(facePatch, fI)
        {
            if( facePatch[fI] >= 0 ) ++nPatchedFaces2;
            seen.insert(facePatch[fI]);
        }
        nUniquePatchVals = seen.size();
    }

    struct FeatureSeg
    {
        point p0, p1;
        label patchA, patchB;
    };
    DynList<FeatureSeg> featureSegs;

    if( nUniquePatchVals > 1 )
    {
        // Phase 2A: face-patch data available -- extract from edgeFaces
        forAll(meshEdges, eI)
        {
            if( edgeFaces.sizeOfRow(eI) != 2 ) continue;
            const label fA = edgeFaces(eI, 0);
            const label fB = edgeFaces(eI, 1);
            const label pA = facePatch[fA];
            const label pB = facePatch[fB];
            // Bug #8: guard unassigned boundary faces (facePatch==-1).
            // A -1 patch on either side creates a spurious feature segment
            // that corrupts corner/edge snap targets at unpatched faces.
            if( pA < 0 || pB < 0 ) continue;
            if( pA == pB ) continue;
            FeatureSeg seg;
            seg.p0 = allPoints[meshEdges[eI][0]];
            seg.p1 = allPoints[meshEdges[eI][1]];
            seg.patchA = pA;
            seg.patchB = pB;
            featureSegs.append(seg);
        }
    }
    else
    {
        // Phase 2B: face-patch data unavailable -- extract from pointPatches
        // An edge is a feature edge if both endpoints share exactly the
        // same set of 2 patch IDs in their pointPatches.
        // Build global-point to boundary-point reverse map
        const labelList& bp = surfaceEngine_.boundaryPoints();
        Map<label> globalToBP;
        forAll(bp, bpI)
            globalToBP.insert(bp[bpI], bpI);

        forAll(meshEdges, eI)
        {
            const label gp0 = meshEdges[eI][0];
            const label gp1 = meshEdges[eI][1];
            Map<label>::const_iterator it0 = globalToBP.find(gp0);
            Map<label>::const_iterator it1 = globalToBP.find(gp1);
            if( it0 == globalToBP.end() || it1 == globalToBP.end() ) continue;
            const label bpI0 = it0();
            const label bpI1 = it1();
            if( pointPatches.sizeOfRow(bpI0) != 2 ) continue;
            if( pointPatches.sizeOfRow(bpI1) != 2 ) continue;
            // Check both endpoints share the same two patches
            const label p0a = pointPatches(bpI0, 0);
            const label p0b = pointPatches(bpI0, 1);
            const label p1a = pointPatches(bpI1, 0);
            const label p1b = pointPatches(bpI1, 1);
            bool match =
                (p0a==p1a && p0b==p1b) ||
                (p0a==p1b && p0b==p1a);
            if( !match ) continue;
            FeatureSeg seg;
            seg.p0 = allPoints[gp0];
            seg.p1 = allPoints[gp1];
            seg.patchA = p0a;
            seg.patchB = p0b;
            featureSegs.append(seg);
        }
    }

    Info << "Feature curve extraction: "
         << featureSegs.size() << " feature edge segments"
         << " (patched faces: " << nPatchedFaces2
         << ", unique patch IDs: " << nUniquePatchVals << ")" << endl;

    // Helper: project point p onto nearest point on segment (a,b)
    auto projectOntoSegment = [](const point& p, const point& a, const point& b) -> point
    {
        const vector ab = b - a;
        const scalar len2 = magSqr(ab);
        if( len2 < VSMALL ) return a;
        const scalar tRaw = ((p-a) & ab) / len2;
        const scalar t = Foam::max(scalar(0), Foam::min(scalar(1), tRaw));
        return a + t*ab;
    };

    // ---------------------------------------------------------------
    // Phase 3 dry-run: feature curve projection statistics
    // For every TWO_PATCH_EDGE point, find nearest feature segment
    // matching its patch pair. Report distance stats. No movement.
    // ---------------------------------------------------------------
    if( featureSegs.size() > 0 && nTwoPatch > 0 )
    {
        const labelList& bp2 = surfaceEngine_.boundaryPoints();
        const pointFieldPMG& pts2 = surfaceEngine_.points();
        label nMatched = 0, nMissed = 0;
        scalar sumDist = 0, maxDist = 0;
        label nLargeMoves = 0;

        forAll(nodesToMap, i)
        {
            const label bpI = nodesToMap[i];
            if( pointClass[bpI] != CLS_TWOPATCH ) continue;

            const point& p = pts2[bp2[bpI]];
            const label pA = pointPatches(bpI, 0);
            const label pB = pointPatches(bpI, 1);

            // Find nearest feature segment matching this patch pair
            scalar bestDSq(GREAT);
            point bestProj(p);
            bool found = false;

            forAll(featureSegs, sI)
            {
                const FeatureSeg& seg = featureSegs[sI];
                bool match =
                    (seg.patchA==pA && seg.patchB==pB) ||
                    (seg.patchA==pB && seg.patchB==pA);
                if( !match ) continue;
                const point proj = projectOntoSegment(p, seg.p0, seg.p1);
                const scalar dSq2 = magSqr(proj - p);
                if( dSq2 < bestDSq )
                { bestDSq = dSq2; bestProj = proj; found = true; }
            }

            if( found )
            {
                ++nMatched;
                const scalar d = Foam::sqrt(bestDSq);
                sumDist += d;
                if( d > maxDist ) maxDist = d;
                if( d > 0.001 ) ++nLargeMoves;
            }
            else
                ++nMissed;
        }

        Info << "FeatureProjectionDryRun:"
             << " two-patch checked=" << (nMatched+nMissed)
             << " matched=" << nMatched
             << " missed=" << nMissed
             << " avgDist=" << (nMatched>0 ? sumDist/nMatched : 0)
             << " maxDist=" << maxDist
             << " largeMoves(>1mm)=" << nLargeMoves << endl;

        // Corner dry-run: inspect multi-patch corner relation to
        // incident feature-curve endpoints. No movement.
        if( nCorner > 0 )
        {
            Info << "CornerDryRun:" << endl;
            forAll(nodesToMap, i)
            {
                const label bpI = nodesToMap[i];
                if( pointClass[bpI] != CLS_CORNER ) continue;
                const point& p = pts2[bp2[bpI]];
                scalar bestSegDSq(GREAT);
                point bestSegProj(p);
                scalar bestEndDSq(GREAT);
                point bestEndPoint(p);
                label nIncidentSegs = 0;
                Info << "  corner bpI=" << bpI << " p=" << p << " patches=(";
                forAllRow(pointPatches, bpI, ppI)
                    Info << pointPatches(bpI, ppI)
                         << (ppI+1<pointPatches.sizeOfRow(bpI) ? "," : "");
                Info << ")" << endl;
                forAll(featureSegs, sI)
                {
                    const FeatureSeg& seg = featureSegs[sI];
                    bool hasA=false, hasB=false;
                    forAllRow(pointPatches, bpI, ppI)
                    {
                        const label pI = pointPatches(bpI, ppI);
                        if( pI==seg.patchA ) hasA=true;
                        if( pI==seg.patchB ) hasB=true;
                    }
                    if( !hasA || !hasB ) continue;
                    ++nIncidentSegs;
                    const point proj = projectOntoSegment(p, seg.p0, seg.p1);
                    const scalar dSqSeg = magSqr(proj - p);
                    if( dSqSeg < bestSegDSq ) { bestSegDSq=dSqSeg; bestSegProj=proj; }
                    const scalar dSq0 = magSqr(seg.p0 - p);
                    const scalar dSq1 = magSqr(seg.p1 - p);
                    if( dSq0 < bestEndDSq ) { bestEndDSq=dSq0; bestEndPoint=seg.p0; }
                    if( dSq1 < bestEndDSq ) { bestEndDSq=dSq1; bestEndPoint=seg.p1; }
                }
                Info << "    incidentSegs=" << nIncidentSegs
                     << " nearestSegDist=" << Foam::sqrt(bestSegDSq)
                     << " nearestEndDist=" << Foam::sqrt(bestEndDSq)
                     << " nearestEndpoint=" << bestEndPoint << endl;
            }
        }
    }

    //- find corner and edge points
    labelLongList selectedCorners, selectedEdges;
    forAll(nodesToMap, i)
    {
        if( cornerPoints.found(nodesToMap[i]) )
        {
            treatedPoint[nodesToMap[i]] = true;
            selectedCorners.append(nodesToMap[i]);
        }
        else if( edgePoints.found(nodesToMap[i]) )
        {
            treatedPoint[nodesToMap[i]] = true;
            selectedEdges.append(nodesToMap[i]);
        }
    }

    //- map the remaining selected points
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const pointFieldPMG& points = surfaceEngine_.points();

    const VRWGraph* bpAtProcsPtr(NULL);
    if( Pstream::parRun() )
        bpAtProcsPtr = &surfaceEngine_.bpAtProcs();

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);
    LongList<parMapperHelper> parallelBndNodes;

    # ifdef USE_OMP
    const label size = nodesToMap.size();
    # pragma omp parallel for if( size > 1000 ) shared(parallelBndNodes) \
    schedule(dynamic, Foam::max(1, size / (3 * omp_get_max_threads())))
    # endif
    forAll(nodesToMap, nI)
    {
        const label bpI = nodesToMap[nI];

        if( treatedPoint[bpI] )
            continue;

        const point& p = points[bPoints[bpI]];
        point mapPoint;
        scalar dSq;
        label nt;

        // Multi-patch corner: snap to nearest incident feature endpoint.
        // Constrained topology motion -- not nearest-patch projection.
        if( pointClass[bpI] == CLS_CORNER && featureSegs.size() > 0 )
        {
            const VRWGraph& ppGraph = surfaceEngine_.pointPoints();
            scalar localLen(GREAT);
            forAllRow(ppGraph, bpI, ppI)
            {
                const label nbpI = ppGraph(bpI, ppI);
                localLen = Foam::min(localLen, mag(points[bPoints[nbpI]] - p));
            }
            scalar bestEndDSq(GREAT);
            point bestEndPt(p);

            // Pass 1: strict -- both segment patches in corner patch set
            scalar strictEndDSq(GREAT);
            point strictEndPt(p);
            forAll(featureSegs, sI)
            {
                const FeatureSeg& seg = featureSegs[sI];
                bool hasA=false, hasB=false;
                forAllRow(pointPatches, bpI, ppI)
                {
                    const label pI = pointPatches(bpI, ppI);
                    if( pI==seg.patchA ) hasA=true;
                    if( pI==seg.patchB ) hasB=true;
                }
                if( !hasA || !hasB ) continue;
                const scalar d0 = magSqr(seg.p0 - p);
                const scalar d1 = magSqr(seg.p1 - p);
                if( d0 < strictEndDSq ) { strictEndDSq=d0; strictEndPt=seg.p0; }
                if( d1 < strictEndDSq ) { strictEndDSq=d1; strictEndPt=seg.p1; }
            }
            if( strictEndDSq < GREAT/2.0 )
            { bestEndDSq=strictEndDSq; bestEndPt=strictEndPt; }
            else
            {
                // Pass 2: relaxed -- one patch matches, tight distance cap
                const scalar relaxedCap = 0.1*localLen;
                forAll(featureSegs, sI)
                {
                    const FeatureSeg& seg = featureSegs[sI];
                    bool hasAny=false;
                    forAllRow(pointPatches, bpI, ppI)
                    {
                        const label pI = pointPatches(bpI, ppI);
                        if( pI==seg.patchA || pI==seg.patchB ) hasAny=true;
                    }
                    if( !hasAny ) continue;
                    const scalar d0 = magSqr(seg.p0 - p);
                    const scalar d1 = magSqr(seg.p1 - p);
                    if( d0 < bestEndDSq && Foam::sqrt(d0) <= relaxedCap )
                    { bestEndDSq=d0; bestEndPt=seg.p0; }
                    if( d1 < bestEndDSq && Foam::sqrt(d1) <= relaxedCap )
                    { bestEndDSq=d1; bestEndPt=seg.p1; }
                }
            }
            const scalar snapDist = Foam::sqrt(bestEndDSq);
            // Cap at 2x local edge length -- dry-run showed snap distances
            // are 0.0002-0.0004m, well within typical cell size of 0.001-0.003m
            const scalar maxSnap = (localLen < GREAT/2.0) ? 2.0*localLen : GREAT;
            if( snapDist <= maxSnap )
            {
                surfaceModifier.moveBoundaryVertexNoUpdate(bpI, bestEndPt);
                treatedPoint[bpI] = true;
                continue;
            }
        }

        // For BL/no-BL interface points use patch-constrained projection
        // onto the BL-side patch only, preventing projection across
        // the feature curve onto the wrong patch.
        Map<label>::const_iterator protIt =
            protectedPointPatches_.find(bpI);
        Map<label>::const_iterator neutIt =
            blNeutralPointPatches_.find(bpI);

        if( protIt != protectedPointPatches_.end() && protIt() >= 0 )
        {
            // BL/no-BL protected point -- constrain to BL-side patch only
            meshOctree_.findNearestSurfacePointInRegion
            (mapPoint, dSq, nt, protIt(), p);
        }
        else if( neutIt != blNeutralPointPatches_.end() && neutIt() >= 0 )
        {
            // BL/neutral point (blade/periodic junction) --
            // constrain projection to BL-side patch only.
            meshOctree_.findNearestSurfacePointInRegion
            (mapPoint, dSq, nt, neutIt(), p);
        }
        else
        {
            // Generic single-patch projection only.
            // Multi-patch points must not be projected onto pointPatches(bpI,0)
            // because that pulls feature-edge/corner points off their intersection.
            // Edge/corner handling is done later by mapEdgeNodes()/mapCorners().
            // Also skips non-manifold (sizeOfRow==0) cleanly.
            if( pointPatches.sizeOfRow(bpI) != 1 )
                continue;

            meshOctree_.findNearestSurfacePointInRegion
            (
                mapPoint,
                dSq,
                nt,
                pointPatches(bpI, 0),
                p
            );
        }

        surfaceModifier.moveBoundaryVertexNoUpdate(bpI, mapPoint);

        if( bpAtProcsPtr && bpAtProcsPtr->sizeOfRow(bpI) )
        {
            # ifdef USE_OMP
            # pragma omp critical
            # endif
            {
                parallelBndNodes.append
                (
                    parMapperHelper
                    (
                        mapPoint,
                        dSq,
                        bpI,
                        -1
                    )
                );
            }
        }

        # ifdef DEBUGMapping
        Info << "Mapped point " << points[boundaryPoints[bpI]] << endl;
        # endif
    }

    //- map vertices at inter-processor boundaries to the nearest location
    //- on the surface
    mapToSmallestDistance(parallelBndNodes);

    //- update face normals, point normals, etc.
    surfaceModifier.updateGeometry(nodesToMap);

    //- map edge nodes
    mapEdgeNodes(selectedEdges);

    // Phase 3: corner snap -- project CLS_CORNER points to nearest
    // incident feature curve endpoint before mapCorners runs.
    // Points successfully snapped are removed from selectedCorners.
    if( featureSegs.size() > 0 )
    {
        const VRWGraph& ppGraph2 = surfaceEngine_.pointPoints();
        const bool debugCornerSnap = false;  // set true to diagnose corner snap
        labelLongList remainingCorners;
        label nSnapped = 0;
        label nRejectedSnaps = 0;
        forAll(selectedCorners, cI)
        {
            const label bpI = selectedCorners[cI];
            if( debugCornerSnap )
            {
            Info << "  selCorner bpI=" << bpI
                 << " cls=" << pointClass[bpI]
                 << " CLS_CORNER=" << CLS_CORNER << endl;
            }
            if( pointClass[bpI] != CLS_CORNER )
            { remainingCorners.append(bpI); continue; }
            const point& cp = points[bPoints[bpI]];
            scalar localLen(GREAT);
            forAllRow(ppGraph2, bpI, ppI)
            {
                const label nbpI = ppGraph2(bpI, ppI);
                localLen = Foam::min(localLen, mag(points[bPoints[nbpI]] - cp));
            }
            scalar bestEndDSq(GREAT);
            point bestEndPt(cp);
            forAll(featureSegs, sI)
            {
                const FeatureSeg& seg = featureSegs[sI];
                bool hasA=false, hasB=false;
                forAllRow(pointPatches, bpI, ppI)
                {
                    const label pI = pointPatches(bpI, ppI);
                    if( pI==seg.patchA ) hasA=true;
                    if( pI==seg.patchB ) hasB=true;
                }
                if( !hasA || !hasB ) continue;
                const scalar d0 = magSqr(seg.p0 - cp);
                const scalar d1 = magSqr(seg.p1 - cp);
                if( d0 < bestEndDSq ) { bestEndDSq=d0; bestEndPt=seg.p0; }
                if( d1 < bestEndDSq ) { bestEndDSq=d1; bestEndPt=seg.p1; }
            }
            const scalar snapDist = Foam::sqrt(bestEndDSq);
            const scalar maxSnap = (localLen < GREAT/2.0) ? 2.0*localLen : GREAT;
            if( debugCornerSnap )
            {
            Info << "  CornerSnapActive snapDist=" << snapDist
                 << " maxSnap=" << maxSnap
                 << " localLen=" << localLen
                 << " bestEndDSq=" << bestEndDSq
                 << " GREAT/2=" << GREAT/2.0 << endl;
            }
            const bool snapCondA = (snapDist <= maxSnap);
            const bool snapCondB = (bestEndDSq < GREAT/2.0);
            if( debugCornerSnap )
            {
            Info << "  SnapConditions: A=" << snapCondA << " B=" << snapCondB << endl;
            }
            if( snapCondA && snapCondB )
            {
                // Snapshot, move, validate pyramid heights, revert if invalid.
                // Phase 3 is serial, so pts is stable during this check.
                const label globalPtI = bPoints[bpI];
                const point oldPos = points[globalPtI];

                // Use the same damped corner motion as mapCorners().
                // Full endpoint snaps can shear BL/prism stacks and create
                // visible bulges at blade/periodic or wall/flow junctions.
                const vector deltaV = bestEndPt - cp;
                const point relaxedEndPt = cp + cornerSnapRelaxation_ * deltaV;

                surfaceModifier.moveBoundaryVertexNoUpdate(bpI, relaxedEndPt);

                bool validSnap = true;
                const VRWGraph& pFacesV = surfaceEngine_.pointFaces();
                const faceList::subList& bFacesV = surfaceEngine_.boundaryFaces();
                const labelList& faceOwnersV = surfaceEngine_.faceOwners();
                const pointFieldPMG& ptsV = surfaceEngine_.points();
                const cellListPMG& cellsV = surfaceEngine_.mesh().cells();
                const faceListPMG& allFacesV = surfaceEngine_.mesh().faces();

                forAllRow(pFacesV, bpI, pfI)
                {
                    const label bfI = pFacesV(bpI, pfI);
                    const face& fV = bFacesV[bfI];

                    point fcV = point::zero;
                    forAll(fV, fpI) fcV += ptsV[fV[fpI]];
                    fcV /= scalar(fV.size());

                    vector fnV = vector::zero;
                    const point& p0V = ptsV[fV[0]];
                    for(label fpI=1; fpI<fV.size()-1; ++fpI)
                        fnV += (ptsV[fV[fpI]] - p0V)
                              ^(ptsV[fV[fpI+1]] - p0V);

                    const label cellIV = faceOwnersV[bfI];
                    point ccV = point::zero;
                    const cell& cllV = cellsV[cellIV];

                    forAll(cllV, cfI)
                    {
                        const face& cfV = allFacesV[cllV[cfI]];
                        point cfcV = point::zero;
                        forAll(cfV, cpI) cfcV += ptsV[cfV[cpI]];
                        ccV += cfcV / scalar(cfV.size());
                    }

                    ccV /= scalar(cllV.size());

                    const scalar hV = fnV & (fcV - ccV);

                    // L3 threshold. localLen is already a length.
                    const scalar fsV = Foam::max
                    (
                        Foam::mag(fnV) * localLen * scalar(1e-6),
                        Foam::pow(localLen, 3) * scalar(1e-12)
                    );

                    if( hV <= -fsV )
                    {
                        validSnap = false;
                        break;
                    }
                }

                if( validSnap )
                {
                    ++nSnapped;

                    if( debugCornerSnap )
                    {
                        Info << "  SNAPPED bpI=" << bpI << endl;
                    }
                }
                else
                {
                    surfaceModifier.moveBoundaryVertexNoUpdate(bpI, oldPos);
                    remainingCorners.append(bpI);
                    ++nRejectedSnaps;
                }
            }
            else
            {
                remainingCorners.append(bpI);
            }
        }
        Info << "CornerSnap: snapped " << nSnapped
             << " of " << selectedCorners.size()
             << " corner points to feature endpoints"
             << ", rejected " << nRejectedSnaps
             << " invalid snaps" << endl;

        // Replace with remaining unsnapped corners for mapCorners().
        // Previously this assignment was trapped inside a malformed
        // debug-only if block, so snapped corners could be reprocessed
        // by mapCorners() and the endpoint snap could be undone.
        selectedCorners = remainingCorners;
    }

    //- map corner vertices
    mapCorners(selectedCorners);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceMapper::repairRejectedPoints()
{
    if( rejectedBpI_.empty() ) return;

    Info << "repairRejectedPoints: repairing "
         << rejectedBpI_.size() << " rejected points" << endl;

    const meshSurfacePartitioner& mPart = meshPartitioner();
    const labelHashSet& cornerPts = mPart.corners();
    const labelHashSet& edgePts   = mPart.edgePoints();
    const VRWGraph& pPatches      = mPart.pointPatches();
    const labelList& bPoints      = surfaceEngine_.boundaryPoints();
    const pointFieldPMG& points   = surfaceEngine_.points();
    const VRWGraph& pointPoints   = surfaceEngine_.pointPoints();
    meshSurfaceEngineModifier sMod(surfaceEngine_);

    // Build mapping distance for movement cap
    labelLongList allBndPts(bPoints.size());
    forAll(allBndPts, i) allBndPts[i] = i;
    scalarList mappingDist;
    findMappingDistance(allBndPts, mappingDist);

    // Expand rejected set by 2 rings
    labelHashSet repairSet;
    forAll(rejectedBpI_, k) repairSet.insert(rejectedBpI_[k]);
    labelHashSet ring1(repairSet);
    forAllConstIter(labelHashSet, ring1, it)
        forAllRow(pointPoints, it.key(), nI)
            repairSet.insert(pointPoints(it.key(), nI));
    labelHashSet ring2(repairSet);
    forAllConstIter(labelHashSet, ring2, it)
        forAllRow(pointPoints, it.key(), nI)
            repairSet.insert(pointPoints(it.key(), nI));

    Info << "repairRejectedPoints: repair neighbourhood = "
         << repairSet.size() << " points" << endl;

    label nMoved = 0;
    label nSkipped = 0;

    forAllConstIter(labelHashSet, repairSet, it)
    {
        const label bpI = it.key();

        // Lock corners, BL/no-BL, non-manifold
        if( cornerPts.found(bpI) ) { ++nSkipped; continue; }
        if( !protectedPoints_.empty() &&
             protectedPoints_.found(bpI) ) { ++nSkipped; continue; }

        const point& oldPos = points[bPoints[bpI]];

        if( edgePts.found(bpI) )
        {
            // TWO_PATCH_EDGE: project to nearest feature curve point
            // Use both patches -- try each, keep minimum distance
            if( pPatches.sizeOfRow(bpI) < 2 ) continue;
            point bestPt = oldPos;
            scalar bestDSq = GREAT;
            for(label pi=0; pi<pPatches.sizeOfRow(bpI); ++pi)
            {
                const label patch = pPatches(bpI, pi);
                point candidate; scalar dSq; label nt;
                meshOctree_.findNearestSurfacePointInRegion
                    (candidate, dSq, nt, patch, oldPos);
                if( dSq < bestDSq &&
                    proposedMoveIsValid(bpI, candidate, oldPos,
                        mappingDist[bpI]) )
                { bestDSq = dSq; bestPt = candidate; }
            }
            if( bestPt != oldPos )
            { sMod.moveBoundaryVertexNoUpdate(bpI, bestPt); ++nMoved; }
        }
        else if( pPatches.sizeOfRow(bpI) == 1 )
        {
            // SINGLE_PATCH: Laplacian smooth + re-project to own patch
            const label myPatch = pPatches(bpI, 0);
            point avg = point::zero;
            label nNei = 0;
            forAllRow(pointPoints, bpI, nI)
            {
                const label nbI = pointPoints(bpI, nI);
                bool sameP = false;
                forAllRow(pPatches, nbI, ppI)
                    if( pPatches(nbI, ppI) == myPatch )
                        { sameP = true; break; }
                if( sameP )
                    { avg += points[bPoints[nbI]]; ++nNei; }
            }
            if( nNei < 2 ) continue;
            avg /= scalar(nNei);
            point projected; scalar dSq; label nt;
            meshOctree_.findNearestSurfacePointInRegion
                (projected, dSq, nt, myPatch, avg);
            if( proposedMoveIsValid(bpI, projected, oldPos,
                    mappingDist[bpI]) )
            { sMod.moveBoundaryVertexNoUpdate(bpI, projected); ++nMoved; }
        }
    }

    // Update geometry for all repaired points
    labelLongList repairList;
    forAllConstIter(labelHashSet, repairSet, it)
        repairList.append(it.key());
    sMod.updateGeometry(repairList);

    Info << "repairRejectedPoints: moved=" << nMoved
         << " skipped=" << nSkipped << endl;
}

void meshSurfaceMapper::smoothSinglePatchPoints(const label nIterations)
{
    Info << "smoothSinglePatchPoints: " << nIterations
         << " iterations" << endl;

    const meshSurfacePartitioner& mPart = meshPartitioner();
    const labelHashSet& cornerPts = mPart.corners();
    const labelHashSet& edgePts   = mPart.edgePoints();
    const VRWGraph& pPatches      = mPart.pointPatches();
    const labelList& bPoints      = surfaceEngine_.boundaryPoints();
    const pointFieldPMG& points   = surfaceEngine_.points();
    const VRWGraph& pointPoints   = surfaceEngine_.pointPoints();

    meshSurfaceEngineModifier sMod(surfaceEngine_);

    // Build mapping distance for movement cap
    labelLongList allBndPts(bPoints.size());
    forAll(allBndPts, i) allBndPts[i] = i;
    scalarList mappingDist;
    findMappingDistance(allBndPts, mappingDist);

    label nMoved = 0;
    label nRejected = 0;

    for(label iter = 0; iter < nIterations; ++iter)
    {
        nMoved = 0; nRejected = 0;
        forAll(bPoints, bpI)
        {
            // Lock: corners, edge points, BL/no-BL, non-manifold
            if( cornerPts.found(bpI) ) continue;
            if( edgePts.found(bpI)   ) continue;
            if( !protectedPoints_.empty() &&
                 protectedPoints_.found(bpI) ) continue;
            // Only smooth single-patch points
            if( pPatches.sizeOfRow(bpI) != 1 ) continue;
            const label myPatch = pPatches(bpI, 0);

            // Laplacian: average of same-patch neighbors
            point avg = point::zero;
            label nNei = 0;
            forAllRow(pointPoints, bpI, nI)
            {
                const label nbI = pointPoints(bpI, nI);
                // Only include neighbors on same patch
                bool sameP = false;
                forAllRow(pPatches, nbI, ppI)
                    if( pPatches(nbI, ppI) == myPatch )
                        { sameP = true; break; }
                if( sameP )
                {
                    avg += points[bPoints[nbI]];
                    ++nNei;
                }
            }
            if( nNei < 2 ) continue;
            avg /= scalar(nNei);

            // Project averaged position back onto own patch
            point projected;
            scalar dSq;
            label nt;
            meshOctree_.findNearestSurfacePointInRegion
            (
                projected, dSq, nt, myPatch,
                avg
            );

            // Accept only if quality valid
            const point& oldPos = points[bPoints[bpI]];
            if( proposedMoveIsValid(bpI, projected, oldPos,
                                    mappingDist[bpI]) )
            {
                sMod.moveBoundaryVertexNoUpdate(bpI, projected);
                ++nMoved;
            }
            else
            {
                ++nRejected;
            }
        }
        sMod.updateGeometry(allBndPts);
        Info << "  smoothSinglePatchPoints iter " << iter
             << ": moved=" << nMoved
             << " rejected=" << nRejected << endl;
        if( nMoved == 0 ) break;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
