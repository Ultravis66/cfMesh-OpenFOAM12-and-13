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
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "helperFunctionsPar.H"
#include "demandDrivenData.H"

#include "labelledPoint.H"
#include "labelledScalar.H"

#include <map>

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGLayer

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayers::findPatchVertices
(
    const boolList& treatPatches,
    List<direction>& pVertices
) const
{
    const meshSurfaceEngine& mse = surfaceEngine();
    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const VRWGraph& pPatches = mPart.pointPatches();

    pVertices.setSize(pPatches.size());
    pVertices = NONE;

    # ifdef USE_OMP
    # pragma omp parallel for if( pPatches.size() > 1000 ) \
    schedule(dynamic, Foam::max(10, pPatches.size()/(2*omp_get_num_threads())))
    # endif
    forAll(pPatches, bpI)
    {
        bool hasTreated(false);
        bool hasNotTreated(false);

        forAllRow(pPatches, bpI, patchI)
        {
            const label patch = pPatches(bpI, patchI);
            if( treatPatches[patch] )
            {
                hasTreated = true;
            }
            else
            {
                hasNotTreated = true;
            }
        }

        if( hasTreated )
        {
            pVertices[bpI] |= PATCHNODE;

            if( hasNotTreated )
                pVertices[bpI] |= EDGENODE;
        }
    }

    if( Pstream::parRun() )
    {
        const VRWGraph& bpAtProcs = mse.bpAtProcs();
        forAll(pVertices, bpI)
            if( pVertices[bpI] && (bpAtProcs.sizeOfRow(bpI) != 0) )
                pVertices[bpI] |= PARALLELBOUNDARY;
    }
}

point boundaryLayers::createNewVertex
(
    const label bpI,
    const boolList& treatPatches,
    const List<direction>& patchVertex
) const
{
    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const vectorField& pNormals = mse.pointNormals();
    const VRWGraph& pFaces = mse.pointFaces();
    const labelList& boundaryFacePatches = mse.boundaryFacePatches();
    const VRWGraph& pointPoints = mse.pointPoints();

    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const VRWGraph& pPatches = mPart.pointPatches();

    const pointFieldPMG& points = mesh_.points();

    # ifdef DEBUGLayer
    Info << "Creating new vertex for boundary vertex " << bpI << endl;
    Info << "Global vertex label " << bPoints[bpI] << endl;
    # endif

    vector normal(vector::zero);
    scalar dist(VGREAT);
    const point& p = points[bPoints[bpI]];
    if( patchVertex[bpI] & EDGENODE )
    {
        # ifdef DEBUGLayer
        Info << "Vertex is on the border" << endl;
        # endif

        DynList<label> otherPatches;
        forAllRow(pPatches, bpI, patchI)
            if( !treatPatches[pPatches(bpI, patchI)] )
                otherPatches.appendIfNotIn
                (
                    pPatches(bpI, patchI)
                );

        if( otherPatches.size() == 1 )
        {
            //- vertex is on an edge (or BL+BL+neutral corner)
            # ifdef DEBUGLayer
            Info << "Vertex is on an edge" << endl;
            # endif

            // BL+BL+neutral corner: use per-patch conservative normal
            // Do NOT apply layerScale_ here — common code below handles it once
            if( blblCornerPoints_.found(bpI) )
            {
                Map<vector> patchNormals;
                forAllRow(pFaces, bpI, pfI)
                {
                    const label faceI = pFaces(bpI, pfI);
                    const label patchLabel = boundaryFacePatches[faceI];
                    if( patchLabel < 0 || patchLabel >= label(treatPatches.size()) ) continue;
                    if( !treatPatches[patchLabel] ) continue;
                    const face& f = bFaces[faceI];
                    if( f.size() < 3 ) continue;
                    vector fn = vector::zero;
                    const point& p0 = points[f[0]];
                    for(label pi=1; pi<f.size()-1; ++pi)
                        fn += (points[f[pi]] - p0) ^ (points[f[pi+1]] - p0);
                    if( patchNormals.found(patchLabel) )
                        patchNormals[patchLabel] += fn;
                    else
                        patchNormals.insert(patchLabel, fn);
                }
                if( patchNormals.size() >= 2 )
                {
                    scalar bestDist = GREAT;
                    vector bestNormal = vector::zero;
                    forAllConstIter(Map<vector>, patchNormals, it)
                    {
                        vector n = it();
                        const scalar magN = mag(n);
                        if( magN < VSMALL ) continue;
                        n /= magN;
                        scalar localDist = VGREAT;
                        forAllRow(pointPoints, bpI, ppI)
                        {
                            const label bpJ = pointPoints(bpI, ppI);
                            const vector vec = points[bPoints[bpJ]] - p;
                            const scalar d = 0.5 * mag(vec & n);
                            if( d < localDist ) localDist = d;
                        }
                        if( localDist < bestDist )
                        {
                            bestDist = localDist;
                            bestNormal = n;
                        }
                    }
                    if( mag(bestNormal) > VSMALL )
                    {
                        normal = bestNormal;
                        dist = bestDist;
                        // layerScale_ applied once by common code below
                    }
                    else
                        normal = pNormals[bpI];
                }
                else
                    normal = pNormals[bpI];
            }
            else
            //- zero-dist for BL-transition edge points
            if( terminateLayersAtConcaveEdges_
             && layerScale_.size() > bpI
             && layerScale_[bpI] < 0.01 )
            {
                dist = 0.0;
            }
            else
            {
            vector v(vector::zero);

            forAllRow(pFaces, bpI, pfI)
            {
                const face& f = bFaces[pFaces(bpI, pfI)];
                const label patchLabel =
                    boundaryFacePatches[pFaces(bpI, pfI)];

                if( treatPatches[patchLabel] )
                {
                    { vector _n=vector::zero; const point& _p0=points[f[0]]; for(label _pi=1;_pi<f.size()-1;++_pi) _n+=(points[f[_pi]]-_p0)^(points[f[_pi+1]]-_p0); normal += _n; }
                }
                else
                {
                    { vector _n=vector::zero; const point& _p0=points[f[0]]; for(label _pi=1;_pi<f.size()-1;++_pi) _n+=(points[f[_pi]]-_p0)^(points[f[_pi+1]]-_p0); v += _n; }
                }
            }

            const scalar magV = mag(v) + VSMALL;
            v /= magV;

            // For BL/no-BL transition zone: skip projection, use pure wall normal
            if( terminateLayersAtConcaveEdges_
             && layerScale_.size() > bpI
             && layerScale_[bpI] < 0.99 )
            {
                const scalar magN = mag(normal) + VSMALL;
                normal /= magN;
            }
            else
            {
                normal -= (normal & v) * v;
            }
            const scalar magN = mag(normal) + VSMALL;
            normal /= magN;

            // Geometry preservation: skip neighbor dist clamping for
            // BL ramp zone points. layerScale_ already handles reduction.
            // The neighbor loop pulls dist to near-zero for transition
            // points near flat inlet/outlet/periodic faces causing warts.
            if( !terminateLayersAtConcaveEdges_
             || layerScale_.size() <= bpI
             || layerScale_[bpI] >= 0.99 )
            {
            forAllRow(pointPoints, bpI, ppI)
            {
                if( patchVertex[pointPoints(bpI, ppI)] )
                    continue;

                const vector vec = points[bPoints[pointPoints(bpI, ppI)]] - p;
                const scalar prod = 0.5 * mag(vec & normal);

                if( prod < dist )
                    dist = prod;
            }
            }
            }  // closes zero-dist else block
        }
        else if( otherPatches.size() == 2 )
        {
            # ifdef DEBUGLayer
            Info << "Vertex is a corner" << endl;
            # endif

            label otherVertex(-1);
            forAllRow(pointPoints, bpI, ppI)
            {
                const label bpJ = pointPoints(bpI, ppI);

                bool found(true);
                forAll(otherPatches, opI)
                    if( !pPatches.contains(bpJ, otherPatches[opI]) )
                    {
                        found = false;
                        break;
                    }

                if( found )
                {
                    otherVertex = bpJ;
                    break;
                }
            }

            if( otherVertex == -1 )
            {
                FatalErrorIn
                (
                    "void boundaryLayers::createNewVertices"
                    "("
                        "const boolList& treatPatches,"
                        "labelList& newLabelForVertex"
                    ")"
                ) << "Cannot find moving vertex!" << exit(FatalError);
            }

            if( terminateLayersAtConcaveEdges_
             && layerScale_.size() > bpI
             && layerScale_[bpI] < 0.01 )
            {
                dist = 0.0;
            }
            else
            {
                //- normal vector is co-linear with that edge
                normal = p - points[bPoints[otherVertex]];
                dist = 0.5 * mag(normal) + VSMALL;
                normal /= 2.0 * dist;
            }
        }
        else
        {
            // Multi-patch singularity: 3+ non-treated patches meet here
            // (e.g. blade/root/periodic triple junction). A single prism
            // extrusion direction is not well-defined. Use surface normal
            // with a minimal safe distance to avoid zero-volume collapse
            // while keeping the extrusion within the local cell geometry.
            normal = pNormals[bpI];

            // Adaptive safe extrusion at multi-patch singularity.
            scalar minEdge(GREAT);
            forAllRow(pointPoints, bpI, ppI)
            {
                const label bpJ = pointPoints(bpI, ppI);
                const scalar d = mag(points[bPoints[bpJ]] - p);
                if( d > VSMALL && d < minEdge )
                    minEdge = d;
            }

            scalar candidateDist =
                (minEdge < GREAT)
              ? Foam::max(scalar(0.02) * minEdge, scalar(100) * VSMALL)
              : scalar(0.0);

            bool accepted(false);

            for(label attempt=0; attempt<8 && candidateDist > VSMALL; ++attempt)
            {
                const point candidate = p - candidateDist * normal;
                bool crossesGuardPlane(false);

                forAllRow(pFaces, bpI, pfI)
                {
                    const label faceI = pFaces(bpI, pfI);
                    const label patchI = boundaryFacePatches[faceI];
                    if( patchI < 0 || patchI >= label(treatPatches.size()) )
                        continue;
                    if( treatPatches[patchI] )
                        continue;
                    const face& f = bFaces[faceI];
                    if( f.size() < 3 )
                        continue;
                    vector fn(vector::zero);
                    const point& fp0 = points[f[0]];
                    for(label pi=1; pi<f.size()-1; ++pi)
                        fn += (points[f[pi]] - fp0) ^ (points[f[pi+1]] - fp0);
                    if( mag(fn) < VSMALL )
                        continue;
                    fn /= mag(fn);
                    point fc(point::zero);
                    forAll(f, fi)
                        fc += points[f[fi]];
                    fc /= scalar(f.size());
                    const scalar s0 = (p - fc) & fn;
                    const scalar s1 = (candidate - fc) & fn;
                    if( mag(s0) > SMALL && s0 * s1 < scalar(0) )
                    {
                        crossesGuardPlane = true;
                        break;
                    }
                }

                if( !crossesGuardPlane )
                {
                    dist = candidateDist;
                    accepted = true;
                    break;
                }

                candidateDist *= scalar(0.5);
            }

            # ifdef DEBUGLayer
            if( accepted )
            {
                Info << "Multi-patch singularity at bpI=" << bpI
                     << " p=" << p
                     << ": accepted dist=" << dist
                     << " (attempt " << attempt << ")" << endl;
            }
            else
            {
                Info << "Multi-patch singularity at bpI=" << bpI
                     << " p=" << p
                     << ": adaptive bisection failed, dist=0" << endl;
            }
            # endif

            if( !accepted )
                dist = 0.0;
        }

        //- limit distances
        forAllRow(pFaces, bpI, pfI)
        {
            const label faceLabel = pFaces(bpI, pfI);
            if( otherPatches.contains(boundaryFacePatches[faceLabel]) )
            {
                const face& f = bFaces[faceLabel];
                const label pos = f.which(bPoints[bpI]);

                if( pos != -1 )
                {
                    const point& ep1 = points[f.prevLabel(pos)];
                    const point& ep2 = points[f.nextLabel(pos)];

                    const scalar dst =
                        help::distanceOfPointFromTheEdge(ep1, ep2, p);

                    if( dst < dist )
                        dist = 0.9 * dst;
                }
                else
                {
                    FatalErrorIn
                    (
                        "void boundaryLayers::createNewVertices"
                        "("
                            "const boolList& treatPatches,"
                            "labelList& newLabelForVertex"
                        ") const"
                    ) << "Face does not contains this vertex!"
                        << abort(FatalError);
                }
            }
        }
    }
    else
    {
        // BL/BL junction: per-patch normals, conservative min-dist selection
        if( blblJunctionPoints_.found(bpI) )
        {
            Map<vector> patchNormals;
            forAllRow(pFaces, bpI, pfI)
            {
                const label faceI = pFaces(bpI, pfI);
                const label patchLabel = boundaryFacePatches[faceI];
                if( patchLabel < 0 || patchLabel >= label(treatPatches.size()) ) continue;
                if( !treatPatches[patchLabel] ) continue;
                const face& f = bFaces[faceI];
                if( f.size() < 3 ) continue;
                vector fn = vector::zero;
                const point& p0 = points[f[0]];
                for(label pi=1; pi<f.size()-1; ++pi)
                    fn += (points[f[pi]] - p0) ^ (points[f[pi+1]] - p0);
                if( patchNormals.found(patchLabel) )
                    patchNormals[patchLabel] += fn;
                else
                    patchNormals.insert(patchLabel, fn);
            }
            if( patchNormals.size() >= 2 )
            {
                scalar bestDist = GREAT;
                vector bestNormal = vector::zero;
                forAllConstIter(Map<vector>, patchNormals, it)
                {
                    vector n = it();
                    const scalar magN = mag(n);
                    if( magN < VSMALL ) continue;
                    n /= magN;
                    scalar localDist = VGREAT;
                    forAllRow(pointPoints, bpI, ppI)
                    {
                        const label bpJ = pointPoints(bpI, ppI);
                        const vector vec = points[bPoints[bpJ]] - p;
                        const scalar d = 0.5 * mag(vec & n);
                        if( d < localDist ) localDist = d;
                    }
                    if( localDist < bestDist )
                    {
                        bestDist = localDist;
                        bestNormal = n;
                    }
                }
                if( mag(bestNormal) > VSMALL )
                {
                    normal = bestNormal;
                    dist = Foam::min(dist, bestDist);
                }
                else
                    normal = pNormals[bpI];
            }
            else
                normal = pNormals[bpI];
        }
        else
        {
            // Feature-aware normal: only average faces from treated patches
            vector patchNormal(vector::zero);
            forAllRow(pFaces, bpI, pfI)
            {
                const label patchLabel = boundaryFacePatches[pFaces(bpI, pfI)];
                if( treatPatches[patchLabel] )
                {
                    const face& f = bFaces[pFaces(bpI, pfI)];
                    vector _n=vector::zero;
                    const point& _p0=points[f[0]];
                    for(label _pi=1;_pi<f.size()-1;++_pi)
                        _n+=(points[f[_pi]]-_p0)^(points[f[_pi+1]]-_p0);
                    patchNormal += _n;
                }
            }
            const scalar magPN = mag(patchNormal);
            if( magPN > VSMALL )
                normal = patchNormal / magPN;
            else
                normal = pNormals[bpI];
        }

        forAllRow(pointPoints, bpI, ppI)
        {
            const scalar d =
            0.5 * mag
            (
                points[bPoints[pointPoints(bpI, ppI)]] -
                p
            );

            if( d < dist )
                dist = d;
        }
    }

    //- create new vertex
    # ifdef DEBUGLayer
    Info << "Normal for vertex " << bpI << " is " << normal << endl;
    Info << "Distance is " << dist << endl;
    # endif

    // Apply layerScale_ ramp at BL/no-BL transition zones
    const scalar rawDist = dist;
    if( terminateLayersAtConcaveEdges_ && layerScale_.size() > bpI )
    {
        const scalar oldDist = dist;
        dist *= layerScale_[bpI];
    }
    if( dist > SMALL )
        dist = Foam::max(dist, VSMALL);
    else if( terminateLayersAtConcaveEdges_
          && layerScale_.size() > bpI
          && layerScale_[bpI] > 0.01 )
    {
        // Geometry forced a nonzero-scale point to near-zero dist.
        // Use rawDist-scaled floor to avoid collapsed layer cells.
        dist = Foam::max
        (
            scalar(1e-6) * rawDist,
            scalar(100) * VSMALL
        );
    }
    else
        dist = 0.0;

    point newP = p - dist * normal;
    if( help::isnan(newP) || help::isinf(newP) )
        return p;
    // BL/neutral crossing clamp: prevent extrusion across periodic/symmetry planes.
    // Fires for blNeutralEdgePoints_ regardless of layerScale_.
    if( !blNeutralEdgePoints_.empty() && blNeutralEdgePoints_.found(bpI) )
    {
        forAllRow(pFaces, bpI, pfI)
        {
            const label faceI = pFaces(bpI, pfI);
            const label patchI = boundaryFacePatches[faceI];
            if( patchI < 0 || patchI >= label(patchNames_.size()) ) continue;
            // Only check neutral patches (patchRole_ == 2)
            if( patchRole_.size() <= patchI ) continue;
            if( patchRole_[patchI] != 2 ) continue;
            const face& f = bFaces[faceI];
            vector fn = vector::zero;
            const point& fp0 = points[f[0]];
            for(label pi=1; pi<f.size()-1; ++pi)
                fn += (points[f[pi]]-fp0)^(points[f[pi+1]]-fp0);
            if( mag(fn) < VSMALL ) continue;
            fn /= mag(fn);
            point fc = point::zero;
            forAll(f, fi) fc += points[f[fi]];
            fc /= scalar(f.size());
            const scalar s0 = (p - fc) & fn;
            const scalar s1 = (newP - fc) & fn;
            if( mag(s0) > SMALL && s0*s1 < scalar(0) )
            {
                // Bisect extrusion distance — hard clamp to p creates
                // zero-thickness BL cells and astronomical aspect ratios.
                scalar localDist = dist;
                point candidate = newP;
                for(label attempt=0; attempt<6; ++attempt)
                {
                    localDist *= scalar(0.5);
                    candidate = p - localDist*normal;
                    const scalar sCand = (candidate - fc) & fn;
                    if( s0*sCand >= scalar(0) )
                    {
                        newP = candidate;
                        break;
                    }
                }
                break;
            }
        }
    }
    // Robust candidate-point clamping near BL/no-BL termination patches
    // Check ALL no-BL faces — use most restrictive constraint.
    // A point must not move further from any no-BL face than its
    // original position. This prevents BL extrusion through inlet,
    // outlet, and periodic surfaces at corner junctions.
    if( terminateLayersAtConcaveEdges_
     && layerScale_.size() > bpI
     && layerScale_[bpI] < 0.99 )
    {
        forAllRow(pFaces, bpI, pfI)
        {
            const label faceI = pFaces(bpI, pfI);
            const label patchI = boundaryFacePatches[faceI];
            if( patchI < 0 || patchI >= label(patchRole_.size()) ) continue;
            if( patchRole_.size() <= label(patchI) || patchRole_[patchI] != 1 ) continue;  // explicit termination only
            const face& f = bFaces[faceI];
            vector fn = vector::zero;
            const point& fp0 = points[f[0]];
            for(label pi=1; pi<f.size()-1; ++pi)
                fn += (points[f[pi]]-fp0)^(points[f[pi+1]]-fp0);
            if( mag(fn) < VSMALL ) continue;
            fn /= mag(fn);
            point fc = point::zero;
            forAll(f, fi) fc += points[f[fi]];
            fc /= scalar(f.size());
            // Signed distance from no-BL face plane
            // fn points outward from no-BL patch into the domain
            // p should be on the positive side (inside domain)
            // If newP goes to negative side, it crossed through the surface
            const scalar s0 = (p - fc) & fn;
            const scalar s1 = (newP - fc) & fn;
            // Only clamp if point actually crossed the face plane
            // s0 > 0 means original point is on correct side
            // s1 < 0 means extruded point crossed to wrong side
            if( mag(s0) > SMALL && s0*s1 < scalar(0) )
            {
                // Bisect extrusion distance — hard clamp to p creates
                // zero-thickness BL cells and astronomical aspect ratios.
                scalar localDist = dist;
                point candidate = newP;
                for(label attempt=0; attempt<6; ++attempt)
                {
                    localDist *= scalar(0.5);
                    candidate = p - localDist*normal;
                    const scalar sCand = (candidate - fc) & fn;
                    if( s0*sCand >= scalar(0) )
                    {
                        newP = candidate;
                        break;
                    }
                }
                break;
            }
        }
    }
    return newP;
}

void boundaryLayers::createNewVertices(const boolList& treatPatches)
{
    Info << "Creating vertices for layer cells" << endl;

    List<direction> patchVertex;
    findPatchVertices(treatPatches, patchVertex);

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();

    //- the following is needed for parallel runs
    //- it is ugly, but must stay for now :(
    if( Pstream::parRun() )
    {
        mse.pointNormals();
        mse.pointPoints();
    }

    pointFieldPMG& points = mesh_.points();

    label nExtrudedVertices(0);
    forAll(patchVertex, bpI)
        if( patchVertex[bpI] )
            ++nExtrudedVertices;

    points.setSize(points.size() + nExtrudedVertices);

    labelLongList procPoints;
    forAll(bPoints, bpI)
        if( patchVertex[bpI] )
        {
            if( patchVertex[bpI] & PARALLELBOUNDARY )
            {
                procPoints.append(bpI);
                continue;
            }

            points[nPoints_] = createNewVertex(bpI, treatPatches, patchVertex);
            newLabelForVertex_[bPoints[bpI]] = nPoints_;
            ++nPoints_;
        }

    if( Pstream::parRun() )
    {
        createNewPartitionVerticesParallel
        (
            procPoints,
            patchVertex,
            treatPatches
        );

        createNewEdgeVerticesParallel
        (
            procPoints,
            patchVertex,
            treatPatches
        );
    }

    //- swap coordinates of new and old points
    forAll(bPoints, bpI)
    {
        const label pLabel = newLabelForVertex_[bPoints[bpI]];
        if( pLabel != -1 )
        {
            const point p = points[pLabel];
            points[pLabel] = points[bPoints[bpI]];
            points[bPoints[bpI]] = p;
        }
    }

    if( nPoints_ != points.size() )
        FatalErrorIn
        (
            "void boundaryLayers::createNewVertices("
            "const meshSurfaceEngine& mse,"
            "const boolList& treatPatches,"
            "labelList& newLabelForVertex)"
        ) << "Number of vertices " << nPoints_
            << " does not match the list size "
            << abort(FatalError);

    // Local topology-aware layer rollback.
    // Restricted to topology-sensitive points only:
    // BL/neutral edges, BL/no-BL edges, BL/BL junction points.
    // Uses correct sign logic and area-scaled tolerance.
    {
        const label maxRollbackIter = 5;
        const scalar dampFactor = 0.5;
        const VRWGraph& ptFacesRB = mse.pointFaces();
        const faceList::subList& bFacesRB = mse.boundaryFaces();
        label nRolledBack = 0;

        // Build restricted rollback set: topology-sensitive points only
        labelHashSet rollbackSet;
        forAllConstIter(labelHashSet, blNeutralEdgePoints_, it)
            rollbackSet.insert(it.key());
        forAllConstIter(labelHashSet, blNoBlEdgePoints_, it)
            rollbackSet.insert(it.key());
        forAllConstIter(labelHashSet, blblJunctionPoints_, it)
            rollbackSet.insert(it.key());

        for(label iter=0; iter<maxRollbackIter; ++iter)
        {
            label nBad = 0;
            forAll(bPoints, bpI)  // Fix: was procPoints (serial no-op)
            {
                if( !rollbackSet.found(bpI) ) continue;
                const label meshPtI = bPoints[bpI];
                const label origPtI = newLabelForVertex_[meshPtI];
                if( origPtI < 0 ) continue;

                // Value copies — avoid reference aliasing during mutation
                const point layerPt = points[meshPtI];
                const point basePt  = points[origPtI];

                bool bad = false;
                forAllRow(ptFacesRB, bpI, pfI)
                {
                    const face& f = bFacesRB[ptFacesRB(bpI, pfI)];
                    point fc = point::zero;
                    forAll(f, fi) fc += points[f[fi]];
                    fc /= scalar(f.size());
                    vector fn = vector::zero;
                    const point& fp0 = points[f[0]];
                    for(label pi=1; pi<f.size()-1; ++pi)
                        fn += (points[f[pi]]-fp0)^(points[f[pi+1]]-fp0);
                    const scalar areaMag = mag(fn);
                    if( areaMag < VSMALL ) continue;
                    const scalar tol = 1e-12 * areaMag;
                    const scalar volLayer = (layerPt - fc) & fn;
                    const scalar volBase  = (basePt  - fc) & fn;
                    // Bad: base and layer on opposite sides of face plane
                    if( mag(volBase) > tol && volBase*volLayer < -tol )
                    { bad = true; break; }
                }

                if( bad )
                {
                    points[meshPtI] = dampFactor*layerPt
                                   + (1.0-dampFactor)*basePt;
                    ++nBad;
                    ++nRolledBack;
                }
            }
            if( nBad == 0 ) break;
        }

        if( nRolledBack > 0 )
            Info << "Layer rollback: " << nRolledBack
                 << " topology-sensitive vertices relaxed" << endl;
    }

    Info << "Finished creating layer vertices" << endl;
}

void boundaryLayers::createNewVertices(const labelList& patchLabels)
{
    otherVrts_.clear();

    patchKey_.setSize(mesh_.boundaries().size());
    patchKey_ = -1;

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    //- populate zeroDistPoints_ for BL-transition vertices
    if( terminateLayersAtConcaveEdges_ )
    {
        boolList _skipPoint(bPoints.size(), false);
        markConcaveEdgePoints(_skipPoint);
    }

    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const VRWGraph& pPatches = mPart.pointPatches();

    //- the following is needed for parallel runs
    //- it is ugly, but must stay for now :(
    mse.boundaryFaces();
    mse.pointNormals();
    mse.pointFaces();
    mse.pointPoints();

    pointFieldPMG& points = mesh_.points();
    boolList treatPatches(mesh_.boundaries().size());
    List<direction> patchVertex(bPoints.size());

    //- make sure than the points are never re-allocated during the process
    points.reserve(points.size() + 2 * bPoints.size());

    //- generate new layer vertices for each patch
    forAll(patchLabels, patchI)
    {
        const label pLabel = patchLabels[patchI];
        treatPatches = false;

        bool treat(true);
        forAll(treatPatchesWithPatch_[pLabel], pI)
        {
            const label otherPatch = treatPatchesWithPatch_[pLabel][pI];
            treatPatches[otherPatch] = true;

            if( patchKey_[otherPatch] == -1 )
            {
                patchKey_[otherPatch] = patchI;
            }
            else
            {
                treat = false;
            }
        }

        if( !treat )
            continue;

        const label pKey = patchKey_[pLabel];

        //- classify vertices belonging to this patch
        findPatchVertices(treatPatches, patchVertex);

        //- create indices and allocate maps for new points
        labelLongList procPoints, patchPoints;
        forAll(bPoints, bpI)
        {
            if( !patchVertex[bpI] )
                continue;

            //- skip vertices at parallel boundaries
            if( patchVertex[bpI] & PARALLELBOUNDARY )
            {
                procPoints.append(bpI);

                continue;
            }

            patchPoints.append(bpI);
            const label pointI = bPoints[bpI];

            if( patchVertex[bpI] & EDGENODE )
            {
                if( otherVrts_.find(pointI) == otherVrts_.end() )
                {
                    std::map<std::pair<label, label>, label> m;

                    otherVrts_.insert(std::make_pair(pointI, m));
                }

                std::pair<label, label> pr(pKey, pKey);
                otherVrts_[pointI].insert(std::make_pair(pr, nPoints_++));
            }
            else
            {
                //- this the only new point
                newLabelForVertex_[pointI] = nPoints_++;
            }
        }

        //- set the size of points
        points.setSize(nPoints_);

        //- calculate coordinates of new points
        // Serial: otherVrts_ is shared std::map — concurrent find/insert/[]
        // is not thread-safe. Races corrupt BL edge/corner vertex positions.
        forAll(patchPoints, i)
        {
            const label bpI = patchPoints[i];

            const label pointI = bPoints[bpI];

            //- create new point
            const point p = createNewVertex(bpI, treatPatches, patchVertex);

            if( patchVertex[bpI] & EDGENODE )
            {
                //- set the new point or an edge point
                if( otherVrts_.find(pointI) == otherVrts_.end() )
                {
                    std::map<std::pair<label, label>, label> m;

                    otherVrts_.insert(std::make_pair(pointI, m));
                }

                std::pair<label, label> pr(pKey, pKey);
                std::map<std::pair<label,label>,label>::const_iterator npIter =
                    otherVrts_[pointI].find(pr);
                if( npIter == otherVrts_[pointI].end() )
                {
                    FatalErrorIn("boundaryLayers::createNewVertices")
                        << "Missing otherVrts_ entry for point " << pointI
                        << abort(FatalError);
                }
                const label npI = npIter->second;
                points[npI] = p;
            }
            else
            {
                //- set the new point
                points[newLabelForVertex_[pointI]] = p;
            }
        }

        if( Pstream::parRun() )
        {
            points.setSize(nPoints_+procPoints.size());

            createNewPartitionVerticesParallel
            (
                procPoints,
                patchVertex,
                treatPatches
            );

            createNewEdgeVerticesParallel
            (
                procPoints,
                patchVertex,
                treatPatches
            );
        }
    }

    //- create missing vertices for edge and corner vertices
    //- they should be stored in the otherNodes map
    forAll(bPoints, bpI)
    {
        const label pointI = bPoints[bpI];

        if( otherVrts_.find(pointI) == otherVrts_.end() )
            continue;

        const point& p = points[pointI];

        std::map<std::pair<label, label>, label>& m = otherVrts_[pointI];
        DynList<label> usedPatches;
        DynList<label> newNodeLabel;
        DynList<vector> newPatchPenetrationVector;
        forAllRow(pPatches, bpI, patchI)
        {
            const label pKey = patchKey_[pPatches(bpI, patchI)];
            const std::pair<label, label> pr(pKey, pKey);
            const std::map<std::pair<label, label>, label>::const_iterator it =
                m.find(pr);
            if( (it != m.end()) && !usedPatches.contains(pKey) )
            {
                usedPatches.append(pKey);
                newNodeLabel.append(it->second);
                newPatchPenetrationVector.append(points[it->second] - p);
            }
        }

        if( newNodeLabel.size() == 1 )
        {
            //- only one patch is treated
            newLabelForVertex_[pointI] = newNodeLabel[0];
            otherVrts_.erase(pointI);
        }
        else if( newNodeLabel.size() == 2 )
        {
            //- point is located at an extrusion edge
            //- create the new position for the existing point
            point newP(p);
            newP += newPatchPenetrationVector[0];
            newP += newPatchPenetrationVector[1];

            if( !help::isnan(newP) && !help::isinf(newP) )
            {
                points.append(newP);
            }
            else
            {
                points.append(p);
            }
            newLabelForVertex_[pointI] = nPoints_;
            ++nPoints_;
        }
        else if( newNodeLabel.size() == 3 )
        {
            //- point is located at an extrusion corner
            //- create 3 points and the new position for the existing point
            point newP(p);
            for(label i=0;i<3;++i)
            {
                newP += newPatchPenetrationVector[i];
                for(label j=i+1;j<3;++j)
                {
                    const point np =
                        p + newPatchPenetrationVector[i] +
                        newPatchPenetrationVector[j];

                    if( !help::isnan(np) && !help::isinf(np) )
                    {
                        points.append(np);
                    }
                    else
                    {
                        points.append(p);
                    }

                    m.insert
                    (
                        std::make_pair
                        (
                            std::make_pair(usedPatches[i], usedPatches[j]),
                            nPoints_
                        )
                    );
                    ++nPoints_;
                }
            }

            //- create new position for the existing point
            if( !help::isnan(newP) && !help::isinf(newP) )
            {
                points.append(newP);
            }
            else
            {
                points.append(p);
            }

            newLabelForVertex_[pointI] = nPoints_;
            ++nPoints_;
        }
        else
        {
            FatalErrorIn
            (
                "void boundaryLayers::createNewVertices("
                "const labelList& patchLabels, labelLongList& newLabelForVertex,"
                "std::map<label, std::map<std::pair<label, label>, label> >&)"
            ) << "Boundary node " << bpI << " is not at an edge!"
                << abort(FatalError);
        }
    }

    //- swap coordinates of new and old points
    # ifdef USE_OMP
    # pragma omp parallel for if( bPoints.size() > 1000 ) \
    schedule(dynamic, 100)
    # endif
    forAll(bPoints, bpI)
    {
        const label pLabel = newLabelForVertex_[bPoints[bpI]];

        if( pLabel != -1 )
        {
            const point p = points[pLabel];
            points[pLabel] = points[bPoints[bpI]];
            points[bPoints[bpI]] = p;
        }
    }

    // Local topology-aware layer rollback (second createNewVertices).
    {
        const meshSurfaceEngine& mseRB = surfaceEngine();
        const VRWGraph& ptFacesRB = mseRB.pointFaces();
        const faceList::subList& bFacesRB = mseRB.boundaryFaces();
        const label maxRollbackIter = 5;
        const scalar dampFactor = 0.5;
        label nRolledBack = 0;

        labelHashSet rollbackSet;
        forAllConstIter(labelHashSet, blNeutralEdgePoints_, it)
            rollbackSet.insert(it.key());
        forAllConstIter(labelHashSet, blNoBlEdgePoints_, it)
            rollbackSet.insert(it.key());
        forAllConstIter(labelHashSet, blblJunctionPoints_, it)
            rollbackSet.insert(it.key());

        for(label iter=0; iter<maxRollbackIter; ++iter)
        {
            label nBad = 0;
            forAll(bPoints, bpI)
            {
                if( !rollbackSet.found(bpI) ) continue;
                const label meshPtI = bPoints[bpI];
                const label origPtI = newLabelForVertex_[meshPtI];
                if( origPtI < 0 ) continue;
                const point layerPt = points[meshPtI];
                const point basePt  = points[origPtI];
                bool bad = false;
                forAllRow(ptFacesRB, bpI, pfI)
                {
                    const face& f = bFacesRB[ptFacesRB(bpI, pfI)];
                    point fc = point::zero;
                    forAll(f, fi) fc += points[f[fi]];
                    fc /= scalar(f.size());
                    vector fn = vector::zero;
                    const point& fp0 = points[f[0]];
                    for(label pi=1; pi<f.size()-1; ++pi)
                        fn += (points[f[pi]]-fp0)^(points[f[pi+1]]-fp0);
                    const scalar areaMag = mag(fn);
                    if( areaMag < VSMALL ) continue;
                    const scalar tol = 1e-12 * areaMag;
                    const scalar volLayer = (layerPt - fc) & fn;
                    const scalar volBase  = (basePt  - fc) & fn;
                    if( mag(volBase) > tol && volBase*volLayer < -tol )
                    { bad = true; break; }
                }
                if( bad )
                {
                    points[meshPtI] = dampFactor*layerPt
                                   + (1.0-dampFactor)*basePt;
                    ++nBad;
                    ++nRolledBack;
                }
            }
            if( nBad == 0 ) break;
        }
        if( nRolledBack > 0 )
            Info << "Layer rollback (pass2): " << nRolledBack
                 << " topology-sensitive vertices relaxed" << endl;
    }
}

void boundaryLayers::createNewPartitionVerticesParallel
(
    const labelLongList& procPoints,
    const List<direction>& pVertices,
    const boolList& /*treatPatches*/
)
{
    if( !Pstream::parRun() )
        return;

    if( returnReduce(procPoints.size(), sumOp<label>()) == 0 )
        return;

    const meshSurfaceEngine& mse = surfaceEngine();
    pointFieldPMG& points = mesh_.points();
    const labelList& bPoints = mse.boundaryPoints();
    const VRWGraph& pointPoints = mse.pointPoints();
    const VRWGraph& bpAtProcs = mse.bpAtProcs();
    const labelList& globalPointLabel = mse.globalBoundaryPointLabel();
    const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();

    scalarField penetrationDistances(bPoints.size(), VGREAT);

    std::map<label, LongList<labelledScalar> > exchangeDistances;

    forAll(procPoints, pointI)
    {
        const label bpI = procPoints[pointI];
        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            if( exchangeDistances.find(neiProc) == exchangeDistances.end() )
            {
                exchangeDistances.insert
                (
                    std::make_pair(neiProc, LongList<labelledScalar>())
                );
            }
        }

        if( pVertices[bpI] & EDGENODE )
            continue;

        scalar dist(VGREAT);
        const point& p = points[bPoints[bpI]];
        forAllRow(pointPoints, bpI, ppI)
        {
            const scalar d =
                0.5 * mag(points[bPoints[pointPoints(bpI, ppI)]] - p);

            if( d < dist )
                dist = d;
        }

        penetrationDistances[bpI] = dist;

        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            exchangeDistances[neiProc].append
            (
                labelledScalar(globalPointLabel[bpI], dist)
            );
        }
    }

    //- exchange distances with other processors
    LongList<labelledScalar> receivedData;
    help::exchangeMap(exchangeDistances, receivedData);
    forAll(receivedData, i)
    {
        if( !globalToLocal.found(receivedData[i].scalarLabel()) ) continue;
        const label bpI = globalToLocal[receivedData[i].scalarLabel()];

        if( penetrationDistances[bpI] > receivedData[i].value() )
            penetrationDistances[bpI] = receivedData[i].value();
    }

    //- Finally, create the points
    // Override point normals at BL/no-BL interface points to use
    // BL-side faces only. The global pointNormals() averages all
    // adjacent faces including no-BL patch faces, which tilts the
    // extrusion normal toward the no-BL surface and causes protrusions.
    vectorField pNormals = mse.pointNormals();
    if( !blNoBlEdgePoints_.empty() )
    {
        const VRWGraph& pFaces = mse.pointFaces();
        const labelList& boundaryFacePatches = mse.boundaryFacePatches();
        const faceList::subList& bFaces = mse.boundaryFaces();
        const pointFieldPMG& pts = mesh_.points();
        const labelList& bPoints = mse.boundaryPoints();

        forAllConstIter(labelHashSet, blNoBlEdgePoints_, iter)
        {
            const label bpI = iter.key();

            // Find the BL-side patch for this point
            Map<label>::const_iterator patchIt = blNoBlPointPatch_.find(bpI);
            if( patchIt == blNoBlPointPatch_.end() || patchIt() < 0 )
                continue;
            const label blPatch = patchIt();

            // Recompute normal using only BL-side faces
            vector blNormal(vector::zero);
            forAllRow(pFaces, bpI, pfI)
            {
                if( boundaryFacePatches[pFaces(bpI, pfI)] != blPatch )
                    continue;
                const face& f = bFaces[pFaces(bpI, pfI)];
                vector fn = vector::zero;
                const point& p0 = pts[f[0]];
                for(label pi=1; pi<f.size()-1; ++pi)
                    fn += (pts[f[pi]]-p0)^(pts[f[pi+1]]-p0);
                blNormal += fn;
            }
            const scalar magN = mag(blNormal);
            if( magN > VSMALL )
                pNormals[bpI] = blNormal / magN;
        }
    }

    // Override point normals at BL/neutral interface points (blade/periodic)
    // Same logic as BL/no-BL: use only BL-side faces to compute normal.
    // Prevents extrusion direction tilting toward periodic plane.
    if( !blNeutralEdgePoints_.empty() )
    {
        const VRWGraph& pFaces2 = mse.pointFaces();
        const labelList& bFacePatches2 = mse.boundaryFacePatches();
        const faceList::subList& bFaces2 = mse.boundaryFaces();
        const pointFieldPMG& pts2 = mesh_.points();

        forAllConstIter(labelHashSet, blNeutralEdgePoints_, iter)
        {
            const label bpI = iter.key();

            Map<label>::const_iterator patchIt = blNeutralPointPatch_.find(bpI);
            if( patchIt == blNeutralPointPatch_.end() || patchIt() < 0 )
                continue;
            const label blPatch = patchIt();

            vector blNormal(vector::zero);
            forAllRow(pFaces2, bpI, pfI)
            {
                if( bFacePatches2[pFaces2(bpI, pfI)] != blPatch )
                    continue;
                const face& f = bFaces2[pFaces2(bpI, pfI)];
                vector fn = vector::zero;
                const point& p0 = pts2[f[0]];
                for(label pi=1; pi<f.size()-1; ++pi)
                    fn += (pts2[f[pi]]-p0)^(pts2[f[pi+1]]-p0);
                blNormal += fn;
            }
            const scalar magN = mag(blNormal);
            if( magN > VSMALL )
                pNormals[bpI] = blNormal / magN;
        }
    }

    forAll(procPoints, pointI)
    {
        const label bpI = procPoints[pointI];

        if( pVertices[bpI] & EDGENODE )
            continue;

        const point& p = points[bPoints[bpI]];
        scalar layerDist = penetrationDistances[bpI];
        if( terminateLayersAtConcaveEdges_
         && layerScale_.size() > bpI )
            layerDist *= layerScale_[bpI];
        const point np = p - pNormals[bpI] * layerDist;
        if( !help::isnan(np) && !help::isinf(np) )
        {
            points[nPoints_] = np;
        }
        else
        {
            points[nPoints_] = p;
        }
        newLabelForVertex_[bPoints[bpI]] = nPoints_;
        ++nPoints_;
    }
}

void boundaryLayers::createNewEdgeVerticesParallel
(
    const labelLongList& procPoints,
    const List<direction>& pVertices,
    const boolList& treatPatches
)
{
    if( !Pstream::parRun() )
        return;

    if( returnReduce(procPoints.size(), sumOp<label>()) == 0 )
        return;

    const meshSurfaceEngine& mse = surfaceEngine();
    pointFieldPMG& points = mesh_.points();
    const labelList& bPoints = mse.boundaryPoints();
    const VRWGraph& pointPoints = mse.pointPoints();
    const VRWGraph& bpAtProcs = mse.bpAtProcs();
    const labelList& globalPointLabel = mse.globalBoundaryPointLabel();
    const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();

    DynList<label> neiProcs;
    labelLongList edgePoints;
    Map<label> bpToEdgePoint;
    forAll(procPoints, pointI)
    {
        const label bpI = procPoints[pointI];
        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            neiProcs.appendIfNotIn(neiProc);
        }

        if( pVertices[bpI] & EDGENODE )
        {
            bpToEdgePoint.insert(bpI, edgePoints.size());
            edgePoints.append(bpI);
        }
    }

    if( returnReduce(edgePoints.size(), sumOp<label>()) == 0 )
        return;

    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const VRWGraph& pPatches = mPart.pointPatches();

    const VRWGraph& pFaces = mse.pointFaces();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& boundaryFacePatches = mse.boundaryFacePatches();

    scalarField dist(edgePoints.size(), VGREAT);
    vectorField normal(edgePoints.size(), vector::zero);
    vectorField v(edgePoints.size(), vector::zero);

    label pKey(-1);
    if( patchKey_.size() )
    {
        forAll(treatPatches, patchI)
            if( treatPatches[patchI] )
            {
                pKey = patchKey_[patchI];
                break;
            }
    }

    forAll(edgePoints, epI)
    {
        const label bpI = edgePoints[epI];
        const point& p = points[bPoints[bpI]];

        //- find patches for the given point
        DynList<label> otherPatches;
        forAllRow(pPatches, bpI, patchI)
            if( !treatPatches[pPatches(bpI, patchI)] )
                otherPatches.appendIfNotIn(pPatches(bpI, patchI));

        //- find local values of normals and v
        if( otherPatches.size() == 1 )
        {
            forAllRow(pFaces, bpI, pfI)
            {
                const face& f = bFaces[pFaces(bpI, pfI)];
                const label patchLabel =
                    boundaryFacePatches[pFaces(bpI, pfI)];

                if( treatPatches[patchLabel] )
                {
                    { vector _n=vector::zero; const point& _p0=points[f[0]]; for(label _pi=1;_pi<f.size()-1;++_pi) _n+=(points[f[_pi]]-_p0)^(points[f[_pi+1]]-_p0); normal[epI] += _n; }
                }
                else
                {
                    { vector _n=vector::zero; const point& _p0=points[f[0]]; for(label _pi=1;_pi<f.size()-1;++_pi) _n+=(points[f[_pi]]-_p0)^(points[f[_pi+1]]-_p0); v[epI] += _n; }
                }
            }
        }
        else if( otherPatches.size() == 2 )
        {
            label otherVertex(-1);
            forAllRow(pointPoints, bpI, ppI)
            {
                const label bpJ = pointPoints(bpI, ppI);

                bool found(true);
                forAll(otherPatches, opI)
                    if( !pPatches.contains(bpJ, otherPatches[opI]) )
                    {
                        found = false;
                        break;
                    }

                if( found )
                {
                    otherVertex = bpJ;
                    break;
                }
            }

            if( otherVertex == -1 )
                continue;

            //- normal vector is co-linear with that edge
            normal[epI] = p - points[bPoints[otherVertex]];
            dist[epI] = mag(normal[epI]);
        }
        else
        {
            // Multi-patch singularity in parallel edge vertex creation.
            // Force local termination — zero extrusion at this point.
            normal[epI] = mse.pointNormals()[bpI];
            dist[epI] = 0.0;
        }
    }

    //- prepare normals and v for sending to other procs
    std::map<label, LongList<labelledPoint> > exchangeNormals;
    forAll(neiProcs, procI)
        exchangeNormals.insert
        (
            std::make_pair(neiProcs[procI], LongList<labelledPoint>())
        );

    forAll(edgePoints, epI)
    {
        const label bpI = edgePoints[epI];

        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            //- store values in the list for sending
            LongList<labelledPoint>& dataToSend = exchangeNormals[neiProc];
            dataToSend.append
            (
                labelledPoint(globalPointLabel[bpI], normal[epI])
            );
            dataToSend.append(labelledPoint(globalPointLabel[bpI], v[epI]));
        }
    }

    //- exchange data with other processors
    LongList<labelledPoint> receivedData;
    help::exchangeMap(exchangeNormals, receivedData);
    exchangeNormals.clear();

    label counter(0);
    while( counter < receivedData.size() )
    {
        const labelledPoint& otherNormal = receivedData[counter++];
        const labelledPoint& otherV = receivedData[counter++];

        if( !globalToLocal.found(otherNormal.pointLabel()) ) continue;
        const label bpI = globalToLocal[otherNormal.pointLabel()];
        normal[bpToEdgePoint[bpI]] += otherNormal.coordinates();
        v[bpToEdgePoint[bpI]] += otherV.coordinates();
    }

    //- calculate normals
    forAll(normal, epI)
    {
        const label bpI = edgePoints[epI];

        //- find patches for the given point
        DynList<label> otherPatches;
        forAllRow(pPatches, bpI, patchI)
            if( !treatPatches[pPatches(bpI, patchI)] )
                otherPatches.appendIfNotIn(pPatches(bpI, patchI));

        if( otherPatches.size() == 1 )
        {
            const scalar magV = mag(v[epI]) + VSMALL;
            v[epI] /= magV;
            normal[epI] -= (normal[epI] & v[epI]) * v[epI];
        }

        const scalar magN = mag(normal[epI]) + VSMALL;
        normal[epI] /= magN;
    }

    //- calculate distances
    forAll(edgePoints, epI)
    {
        const label bpI = edgePoints[epI];
        const point& p = points[bPoints[bpI]];

        //- find patches for the given point
        DynList<label> otherPatches;
        forAllRow(pPatches, bpI, patchI)
            if( !treatPatches[pPatches(bpI, patchI)] )
                otherPatches.appendIfNotIn(pPatches(bpI, patchI));

        if( otherPatches.size() == 1 )
        {
            forAllRow(pointPoints, bpI, ppI)
            {
                if( pVertices[pointPoints(bpI, ppI)] )
                    continue;

                const vector vec = points[bPoints[pointPoints(bpI, ppI)]] - p;
                const scalar prod = 0.5 * mag(vec & normal[epI]);

                if( prod < dist[epI] )
                    dist[epI] = prod;
            }
        }

        //- limit distances
        forAllRow(pFaces, bpI, pfI)
        {
            const label faceLabel = pFaces(bpI, pfI);
            if( otherPatches.contains(boundaryFacePatches[faceLabel]) )
            {
                const face& f = bFaces[faceLabel];
                const label pos = f.which(bPoints[bpI]);

                if( pos != -1 )
                {
                    const point& ep1 = points[f.prevLabel(pos)];
                    const point& ep2 = points[f.nextLabel(pos)];

                    const scalar dst =
                        help::distanceOfPointFromTheEdge(ep1, ep2, p);

                    if( dst < dist[epI] )
                        dist[epI] = 0.9 * dst;
                }
                else
                {
                    FatalErrorIn
                    (
                        "void boundaryLayers::createNewEdgeVerticesParallel"
                        "("
                            "const labelLongList& procPoints,"
                            "const List<direction>& pVertices,"
                            "const boolList& treatPatches,"
                            "labelList& newLabelForVertex"
                        ") const"
                    ) << "Face does not contains this vertex!"
                        << abort(FatalError);
                }
            }
        }
    }

    //- exchange distances with other processors
    std::map<label, LongList<labelledScalar> > exchangeDistances;
    forAll(neiProcs, procI)
        exchangeDistances.insert
        (
            std::make_pair(neiProcs[procI], LongList<labelledScalar>())
        );

    forAll(edgePoints, epI)
    {
        const label bpI = edgePoints[epI];
        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            LongList<labelledScalar>& ls = exchangeDistances[neiProc];
            ls.append(labelledScalar(globalPointLabel[bpI], dist[epI]));
        }
    }

    //- exchange distances with other processors
    LongList<labelledScalar> receivedDistances;
    help::exchangeMap(exchangeDistances, receivedDistances);
    exchangeDistances.clear();

    forAll(receivedDistances, i)
    {
        if( !globalToLocal.found(receivedDistances[i].scalarLabel()) ) continue;
        const label bpI = globalToLocal[receivedDistances[i].scalarLabel()];
        const label epI = bpToEdgePoint[bpI];
        if( dist[epI] > receivedDistances[i].value() )
            dist[epI] = receivedDistances[i].value();
    }

    //- Finally, create new points
    forAll(edgePoints, epI)
    {
        const label bpI = edgePoints[epI];

        const point& p = points[bPoints[bpI]];
        const point np = p - normal[epI] * dist[epI];
        if( !help::isnan(np) && !help::isinf(np) )
        {
            points[nPoints_] = np;
        }
        else
        {
            points[nPoints_] = p;
        }

        if( pKey == -1 )
        {
            //- extrusion for one patch in a single go
            newLabelForVertex_[bPoints[bpI]] = nPoints_;
        }
        else
        {
            const label pointI = bPoints[bpI];

            if( otherVrts_.find(pointI) == otherVrts_.end() )
            {
                std::map<std::pair<label, label>, label> m;
                otherVrts_.insert(std::make_pair(pointI, m));
            }

            std::pair<label, label> pr(pKey, pKey);
            otherVrts_[pointI].insert(std::make_pair(pr, nPoints_));
        }
        ++nPoints_;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
