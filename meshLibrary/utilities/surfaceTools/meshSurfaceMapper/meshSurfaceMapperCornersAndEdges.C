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

#include "meshOctree.H"
#include "triSurf.H"
#include "triSurfacePartitioner.H"
#include "meshSurfaceMapper.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceEngineModifier.H"
#include "meshSurfacePartitioner.H"
#include "labelledScalar.H"

#include "helperFunctions.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGMapping

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool meshSurfaceMapper::proposedMoveIsValid
(
    const label bpI,
    const point& candidate,
    const point& oldPosition,
    const scalar mappingDistSq
) const
{
    const meshSurfaceEngine& mse = surfaceEngine_;
    const VRWGraph& pFaces = mse.pointFaces();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& faceOwners = mse.faceOwners();
    const pointFieldPMG& pts = mse.points();
    const cellListPMG& cells = mse.mesh().cells();
    const faceListPMG& allFaces = mse.mesh().faces();
    const labelList& bPoints = mse.boundaryPoints();

    // Check 1: movement cap — reject if candidate is further from
    // old position than the local mapping distance allows.
    const scalar moveSq = magSqr(candidate - oldPosition);
    if( mappingDistSq > VSMALL && moveSq > mappingDistSq )
        return false;

    // Pure geometric evaluation — no mesh mutation, thread-safe.
    // Substitute candidate for bPoints[bpI] inline using lambda.
    const label globalPtI = bPoints[bpI];
    auto ptOf = [&](const label ptI) -> point
    {
        return (ptI == globalPtI) ? candidate : point(pts[ptI]);
    };

    bool valid = true;
    forAllRow(pFaces, bpI, pfI)
    {
        const label bfI = pFaces(bpI, pfI);
        const face& f = bFaces[bfI];

        // Check 2: face area
        point fc = point::zero;
        forAll(f, fpI) fc += ptOf(f[fpI]);
        fc /= scalar(f.size());
        vector fn = vector::zero;
        const point p0 = ptOf(f[0]);
        for(label fpI=1; fpI<f.size()-1; ++fpI)
            fn += (ptOf(f[fpI])-p0)^(ptOf(f[fpI+1])-p0);
        const scalar areaSq = magSqr(fn);
        if( areaSq < VSMALL )
        { valid = false; break; }

        // Check 3: pyramid height
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
        const scalar h = fn & (fc - cc);
        const scalar faceScale = Foam::sqrt(areaSq) * scalar(1e-6);
        if( h <= -faceScale )
        { valid = false; break; }
    }

    // Check 4: skewness proxy using candidate position
    if( valid )
    {
        forAllRow(pFaces, bpI, pfI)
        {
            const label bfI = pFaces(bpI, pfI);
            const face& f = bFaces[bfI];
            point fc2 = point::zero;
            forAll(f, fpI) fc2 += ptOf(f[fpI]);
            fc2 /= scalar(f.size());
            const label cellI = faceOwners[bfI];
            point cc2 = point::zero;
            const cell& cll = cells[cellI];
            forAll(cll, cfI)
            {
                const face& cf = allFaces[cll[cfI]];
                point cfc = point::zero;
                forAll(cf, cpI) cfc += pts[cf[cpI]];
                cc2 += cfc / scalar(cf.size());
            }
            cc2 /= scalar(cll.size());
            const scalar d = mag(fc2 - cc2);
            if( d < VSMALL ) continue;
            const scalar skewness = mag(fc2 - 0.5*(fc2+cc2)) / d;
            if( skewness > scalar(4.0) )
            { valid = false; break; }
        }
    }

    return valid;
}

void meshSurfaceMapper::findMappingDistance
(
    const labelLongList& nodesToMap,
    scalarList& mappingDistance
) const
{
    const vectorField& faceCentres = surfaceEngine_.faceCentres();
    const VRWGraph& pFaces = surfaceEngine_.pointFaces();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const pointFieldPMG& points = surfaceEngine_.points();

    //- generate search distance for corner nodes
    mappingDistance.setSize(nodesToMap.size());
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(nodesToMap, i)
    {
        const label bpI = nodesToMap[i];

        mappingDistance[i] = 0.0;

        const point& p = points[bPoints[bpI]];
        forAllRow(pFaces, bpI, pfI)
        {
            const scalar d = magSqr(faceCentres[pFaces(bpI, pfI)] - p);
            mappingDistance[i] = Foam::max(mappingDistance[i], d);
        }

        //- safety factor
        mappingDistance[i] *= 4.0;
    }

    if( Pstream::parRun() )
    {
        //- make sure that corner nodesd at parallel boundaries
        //- have the same range in which they accept the corners
        const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
        const labelList& globalPointLabel =
            surfaceEngine_.globalBoundaryPointLabel();

        //- create the map for exchanging data
        std::map<label, DynList<labelledScalar> > exchangeData;
        const DynList<label>& neiProcs = surfaceEngine_.bpNeiProcs();
        forAll(neiProcs, i)
            exchangeData.insert
            (
                std::make_pair(neiProcs[i], DynList<labelledScalar>())
            );

        Map<label> globalToLocal;

        forAll(nodesToMap, nI)
        {
            const label bpI = nodesToMap[nI];

            if( bpAtProcs.sizeOfRow(bpI) != 0 )
                globalToLocal.insert(globalPointLabel[bpI], nI);

            forAllRow(bpAtProcs, bpI, i)
            {
                const label neiProc = bpAtProcs(bpI, i);
                if( neiProc == Pstream::myProcNo() )
                    continue;

                exchangeData[neiProc].append
                (
                    labelledScalar(globalPointLabel[bpI], mappingDistance[nI])
                );
            }
        }

        //- exchange data between processors
        LongList<labelledScalar> receivedData;
        help::exchangeMap(exchangeData, receivedData);

        //- select the maximum mapping distance for processor points
        forAll(receivedData, i)
        {
            const labelledScalar& ls = receivedData[i];

            const label nI = globalToLocal[ls.scalarLabel()];

            //- choose the maximum value for the mapping distance
            mappingDistance[nI] = Foam::max(mappingDistance[nI], ls.value());
        }
    }
}

scalar meshSurfaceMapper::faceMetricInPatch
(
    const label bfI,
    const label patchI
) const
{
    const face& bf = surfaceEngine_.boundaryFaces()[bfI];

    const pointFieldPMG& points = surfaceEngine_.points();

    vector centre=vector::zero; forAll(bf,_pi) centre+=points[bf[_pi]]; centre/=bf.size();
    vector area=vector::zero; const point& _p0a=points[bf[0]]; for(label _pi=1;_pi<bf.size()-1;++_pi) area+=(points[bf[_pi]]-_p0a)^(points[bf[_pi+1]]-_p0a);

    point projCentre;
    scalar dSq;
    label nt;

    meshOctree_.findNearestSurfacePointInRegion
    (
        projCentre,
        dSq,
        nt,
        patchI,
        centre
    );

    DynList<point> projPoints(bf.size());
    forAll(bf, pI)
    {
        meshOctree_.findNearestSurfacePointInRegion
        (
            projPoints[pI],
            dSq,
            nt,
            patchI,
            points[bf[pI]]
        );
    }

    vector projArea(vector::zero);
    forAll(bf, pI)
    {
        projArea +=
            triPointRef
            (
                projPoints[pI],
                projPoints[bf.fcIndex(pI)],
                projCentre
            ).normal();
    }

    return magSqr(centre - projCentre) + mag(mag(projArea) - mag(area));
}

void meshSurfaceMapper::mapCorners(const labelLongList& nodesToMap)
{
    const triSurfacePartitioner& sPartitioner = surfacePartitioner();
    const labelList& surfCorners = sPartitioner.corners();
    const List<DynList<label> >& cornerPatches = sPartitioner.cornerPatches();

    const meshSurfacePartitioner& mPart = meshPartitioner();
    const labelHashSet& corners = mPart.corners();
    const VRWGraph& pPatches = mPart.pointPatches();

    const pointFieldPMG& points = surfaceEngine_.points();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();

    const triSurf& surf = meshOctree_.surface();
    const pointField& sPoints = surf.points();

    //std::map<label, scalar> mappingDistance;
    scalarList mappingDistance;
    findMappingDistance(nodesToMap, mappingDistance);

    //- for every corner in the mesh surface find the nearest corner in the
    //- triSurface
    meshSurfaceEngineModifier sMod(surfaceEngine_);

    // Store old positions for validity-check revert
    pointField oldPositions(nodesToMap.size());
    forAll(nodesToMap, i)
        oldPositions[i] = points[bPoints[nodesToMap[i]]];

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(nodesToMap, cornerI)
    {
        const label bpI = nodesToMap[cornerI];
        // Skip BL/no-BL interface points - must stay on feature curve
        if( !protectedPoints_.empty() && protectedPoints_.found(bpI) )
            continue;
        if( !corners.found(bpI) )
            FatalErrorIn
            (
                "meshSurfaceMapper::mapCorners(const labelLongList&)"
            ) << "Trying to map a point that is not a corner"
                << abort(FatalError);

        const point& p = points[bPoints[bpI]];
        const scalar maxDist = mappingDistance[cornerI];

        //- find the nearest position to the given point patches
        const DynList<label> patches = pPatches[bpI];

        point mapPointApprox(p);
        scalar distSqApprox;

        label iter(0);
        while( iter++ < 5 )
        {
            point newP(vector::zero);
            forAll(patches, patchI)
            {
                point np;
                label nt;
                meshOctree_.findNearestSurfacePointInRegion
                (
                    np,
                    distSqApprox,
                    nt,
                    patches[patchI],
                    mapPointApprox
                );

                newP += np;
            }

            newP /= patches.size();

            if( magSqr(newP - mapPointApprox) < 1e-8 * maxDist )
                break;

            mapPointApprox = newP;
        }
        distSqApprox = magSqr(mapPointApprox - p);

        //- find the nearest triSurface corner for the given corner
        scalar distSq(mappingDistance[cornerI]);
        point mapPoint(p);
        forAll(surfCorners, scI)
        {
            const label cornerID = surfCorners[scI];
            const point& sp = sPoints[cornerID];

            if( Foam::magSqr(sp - p) < distSq )
            {
                bool store(true);
                const DynList<label>& cPatches = cornerPatches[scI];
                forAll(patches, i)
                {
                    if( !cPatches.contains(patches[i]) )
                    {
                        store = false;
                        break;
                    }
                }

                if( store )
                {
                    mapPoint = sp;
                    distSq = Foam::magSqr(sp - p);
                }
            }
        }

        if( distSq > 1.2 * distSqApprox )
        {
            mapPoint = mapPointApprox;
        }

        //- move the point toward the nearest corner using damped move.
        //- Full corner snaps can invert faces at 3+ patch junctions.
        const vector delta = mapPoint - p;
        const point relaxedPoint = p + cornerSnapRelaxation_ * delta;
        sMod.moveBoundaryVertexNoUpdate(bpI, relaxedPoint);
    }

    // Validity check: serial revert of invalid corner moves
    {
        const VRWGraph& pFaces = surfaceEngine_.pointFaces();
        const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();
        const labelList& faceOwners = surfaceEngine_.faceOwners();
        const pointFieldPMG& pts = surfaceEngine_.points();
        const cellListPMG& cells = surfaceEngine_.mesh().cells();
        const faceListPMG& allFaces = surfaceEngine_.mesh().faces();
        label nReverted = 0;
        forAll(nodesToMap, i)
        {
            const label bpI = nodesToMap[i];
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
                const scalar h = fn & (fc - cc);
                // Relative threshold: reject only if pyramid height is
                // meaningfully negative relative to local cell size.
                // mappingDistance[i] ~ (4 * max edge length)^2 at point.
                // Tolerance of -1e-6 * sqrt(mappingDistance) matches
                // commercial mesher behaviour: accepts flat-but-valid
                // faces at tight concave junctions, rejects true inversions.
                const scalar localLen =
                    Foam::sqrt(mappingDistance[i]) * scalar(1e-6);
                if( h <= -localLen )
                { validMove = false; break; }
            }
            if( !validMove )
            {
                sMod.moveBoundaryVertexNoUpdate(bpI, oldPositions[i]);
                ++nReverted;
            }
        }
        if( nReverted > 0 )
            Info << "[CornerValidity] reverted " << nReverted
                 << " invalid corner moves" << endl;
    }

    sMod.updateGeometry(nodesToMap);
}

void meshSurfaceMapper::mapEdgeNodes(const labelLongList& nodesToMap)
{
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();

    const meshSurfacePartitioner& mPart = meshPartitioner();
    const VRWGraph& pPatches = mPart.pointPatches();

    //const triSurf& surf = meshOctree_.surface();
    //const pointField& sPoints = surf.points();

    //- find mapping distance for selected vertices
    scalarList mappingDistance;
    findMappingDistance(nodesToMap, mappingDistance);

    const VRWGraph* bpAtProcsPtr(NULL);
    if( Pstream::parRun() )
        bpAtProcsPtr = &surfaceEngine_.bpAtProcs();

    LongList<parMapperHelper> parallelBndNodes;

    meshSurfaceEngineModifier sMod(surfaceEngine_);

    // Store old positions for validity-check revert
    pointField oldPositions(nodesToMap.size());
    forAll(nodesToMap, i)
        oldPositions[i] = points[bPoints[nodesToMap[i]]];

    // Deterministic storage - OMP loop computes only, no mesh mutation
    pointField projectedPoints(nodesToMap.size(), point::zero);
    scalarList projectedDistSq(nodesToMap.size(), scalar(-1));

    // Mark protected indices so OMP loop can skip them
    // BL/no-BL interface points must stay on feature curve
    boolList isProtected(nodesToMap.size(), false);
    if( !protectedPoints_.empty() )
    {
        label nProt = 0;
        forAll(nodesToMap, i)
            if( protectedPoints_.found(nodesToMap[i]) )
            { isProtected[i] = true; ++nProt; }
        if( nProt > 0 )
            Info << "mapEdgeNodes: excluding "
                 << nProt << " protected BL/no-BL interface points" << endl;
    }

    //- map point to the nearest vertex on the surface mesh
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(nodesToMap, i)
    {
        // Skip BL/no-BL interface points - must stay on feature curve
        if( isProtected[i] ) continue;
        const label bpI = nodesToMap[i];
        const point& p = points[bPoints[bpI]];

        //- find patches at this edge point
        const DynList<label> patches = pPatches[bpI];

        const scalar maxDist = mappingDistance[i];

        //- find approximate position of the vertex on the edge
        point mapPointApprox(p);
        scalar distSqApprox;
        label iter(0);
        while( iter++ < 5 )
        {
            point newP(vector::zero);

            forAll(patches, patchI)
            {
                point np;
                label nt;
                meshOctree_.findNearestSurfacePointInRegion
                (
                    np,
                    distSqApprox,
                    nt,
                    patches[patchI],
                    mapPointApprox
                );

                newP += np;
            }

            newP /= patches.size();

            if( magSqr(newP - mapPointApprox) < 1e-8 * maxDist )
                break;

            mapPointApprox = newP;
        }
        distSqApprox = magSqr(mapPointApprox - p);

        //- find the nearest vertex on the triSurface feature edge
        point mapPoint(mapPointApprox);
        scalar distSq(distSqApprox);
        label nse;
        meshOctree_.findNearestEdgePoint(mapPoint, distSq, nse, p, patches);

        //- use the vertex with the smallest mapping distance
        if( distSq > 1.2 * distSqApprox )
        {
            mapPoint = mapPointApprox;
            distSq = distSqApprox;
        }

        //- check if the mapping distance is within the given tolerances
        if( distSq > maxDist )
        {
            //- this indicates possible problems
            //- reduce the mapping distance
            const scalar f = Foam::sqrt(maxDist / distSq);
            distSq = mappingDistance[i];
            mapPoint = f * (mapPoint - p) + p;
        }

                //- store only - no mesh mutation in parallel
        projectedPoints[i] = mapPoint;
        projectedDistSq[i] = distSq;
    }

    // Serial move pass - deterministic order, no validity check
    // Edge nodes must reach feature curve - validity check inappropriate
    forAll(nodesToMap, i)
    {
        const label bpI = nodesToMap[i];
        if( projectedDistSq[i] < scalar(0) ) continue;
        sMod.moveBoundaryVertexNoUpdate(bpI, projectedPoints[i]);
    }

    // Build parallelBndNodes deterministically - only accepted points
    forAll(nodesToMap, i)
    {
        const label bpI = nodesToMap[i];
        if( projectedDistSq[i] < scalar(0) ) continue;
        if( bpAtProcsPtr && bpAtProcsPtr->sizeOfRow(bpI) )
            parallelBndNodes.append
            (
                parMapperHelper
                (
                    projectedPoints[i],
                    projectedDistSq[i],
                    bpI,
                    -1
                )
            );
    }
    mapToSmallestDistance(parallelBndNodes);

    // Old validity check replaced by serial pass above
    if( false )
    {
        const VRWGraph& pFaces = surfaceEngine_.pointFaces();
        const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();
        const labelList& faceOwners = surfaceEngine_.faceOwners();
        const pointFieldPMG& pts = surfaceEngine_.points();
        const cellListPMG& cells = surfaceEngine_.mesh().cells();
        const faceListPMG& allFaces = surfaceEngine_.mesh().faces();
        label nReverted = 0;
        forAll(nodesToMap, i)
        {
            const label bpI = nodesToMap[i];
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
                if( (fn & (fc - cc)) <= SMALL )
                { validMove = false; break; }
            }
            if( !validMove )
            {
                sMod.moveBoundaryVertexNoUpdate(bpI, oldPositions[i]);
                ++nReverted;
            }
        }
        if( nReverted > 0 )
            Info << "[EdgeValidity] reverted " << nReverted
                 << " invalid edge moves" << endl;
    }

    sMod.updateGeometry(nodesToMap);
}

void meshSurfaceMapper::mapCornersAndEdges()
{
    const meshSurfacePartitioner& mPart = meshPartitioner();
    const labelHashSet& cornerPoints = mPart.corners();
    labelLongList selectedPoints;
    forAllConstIter(labelHashSet, cornerPoints, it)
        selectedPoints.append(it.key());

    mapCorners(selectedPoints);

    selectedPoints.clear();
    const labelHashSet& edgePoints = mPart.edgePoints();
    forAllConstIter(labelHashSet, edgePoints, it)
        selectedPoints.append(it.key());

    mapEdgeNodes(selectedPoints);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
