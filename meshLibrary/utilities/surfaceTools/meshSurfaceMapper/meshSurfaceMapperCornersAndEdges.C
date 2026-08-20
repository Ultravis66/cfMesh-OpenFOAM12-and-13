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

    // Check 1: movement cap -- reject if candidate is further from
    // old position than the local mapping distance allows.
    const scalar moveSq = magSqr(candidate - oldPosition);
    if( mappingDistSq > VSMALL && moveSq > mappingDistSq )
        return false;

    // Pure geometric evaluation -- no mesh mutation, thread-safe.
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
            forAll(cf, cpI) cfc += ptOf(cf[cpI]);
            cc += cfc / scalar(cf.size());
        }
        cc /= scalar(cll.size());
        const scalar h = fn & (fc - cc);

        // L3 threshold: h = fn & (fc - cc) has area*length units.
        // Use the local mapping distance as the length scale, matching
        // the post-move validity guard used later in this file.
        const scalar localLen = Foam::sqrt(mappingDistSq + VSMALL);
        const scalar faceScale = Foam::max
        (
            Foam::sqrt(areaSq) * localLen * scalar(1e-6),
            Foam::pow(localLen, 3) * scalar(1e-12)
        );

        if( h <= -faceScale )
        { valid = false; break; }
    }

    // Check 4 removed: skewness proxy fc2-0.5*(fc2+cc2) always equals
    // 0.5*(fc2-cc2), giving ratio always 0.5 -- never rejected anything.
    // Rely on movement cap + face area + pyramid height checks above.

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

            if( !globalToLocal.found(ls.scalarLabel()) ) continue;
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

// v7.1 forensic diagnostic (SOL-specified): brute-force A/B triSurf
// graph probe for two hard-coded control junctions. Serial, read-only,
// no mesh mutation, no OMP. Answers: does a raw ALL_3 (ball touching
// all three required patches) triSurf vertex exist near this query
// point, and if not, what patch-pair vertices exist instead and how
// far apart are they? Distinguishes "findNearestCorner() rejects a
// real ALL_3 vertex" from "no ALL_3 vertex exists at all" -- these
// require different fixes (corner-recognition criteria vs. surface
// welding/topology repair).
// Forward declaration -- probeTriSurfJunction (below) calls
// probeCornerEdgeFan (defined further below it) automatically on every
// ALL_3 vertex it finds.
static void probeCornerEdgeFan
(
    const triSurf& surf,
    const label surfacePointId,
    const DynList<label>& requiredPatches,
    const point& expectedPos,
    const label reqBlade,
    const label reqPeriodic,
    const label reqWall
);

// v7.1 Phase 2a (SOL): isolated, hard-coded prototype repair for
// vertex 3748's blade_3|shroud split contact-line interface ONLY.
// Completely separate from the live pipeline -- builds a genuinely
// mutable COPY of the surface data via triSurfModifier, remaps the
// two known open-edge endpoints (3766, 3767) to a shared midpoint,
// constructs a FRESH triSurf from the modified data (clean addressing
// from birth, per SOL -- avoids any stale-cache question entirely),
// then re-probes the repaired copy to verify the exact false->true
// transition. Never touches meshOctree_.surface() / surfacePtr_ /
// the real pipeline in any way. Diagnostic-only in the sense that it
// changes nothing about mesh generation -- it operates entirely on a
// throwaway local copy for verification purposes.
static void probeTriSurfJunction
(
    const triSurf& surf,
    const point& queryPoint,
    const scalar mappingDistanceSq,
    const label queryBpI,
    const DynList<label>& requiredPatches,
    const labelHashSet& bladeIds,
    const labelHashSet& periodicIds,
    const labelHashSet& wallIds
)
{
    const pointField& surfPts = surf.points();
    const VRWGraph& ptFacets = surf.pointFacets();
    const VRWGraph& ptEdges = surf.pointEdges();
    const LongList<labelledTri>& facets = surf.facets();
    const wordList& patchNames = surf.patchNames();

    const scalar mappingRadius = Foam::sqrt(mappingDistanceSq);
    const scalar inspectRadius = 2.0 * mappingRadius;
    const scalar inspectRadiusSq = inspectRadius * inspectRadius;

    std::string reqNames;
    forAll(requiredPatches, i)
    {
        if( i > 0 ) reqNames += " ";
        const label rp = requiredPatches[i];
        reqNames += (rp >= 0 && rp < label(patchNames.size())) ?
            std::string(patchNames[rp].c_str()) : std::string("?");
    }

    Info << "[TriSurfProbeSummary] queryBpI=" << queryBpI
         << " required=(" << reqNames.c_str() << ")"
         << " queryPos=" << queryPoint
         << " mappingRadius=" << mappingRadius
         << " inspectionRadius=" << inspectRadius
         << endl;

    label nNearby = 0;
    label nAll3 = 0;
    label nBladePeriodic = 0, nBladeWall = 0, nPeriodicWall = 0;

    // For pairwise-separation computation on the three partial classes.
    // requiredPatches is expected size 3: [blade, periodic, wall].
    // Store the NEAREST vertex position seen for each partial class.
    bool haveBladePeriodic = false, haveBladeWall = false, havePeriodicWall = false;
    point posBladePeriodic(vector::zero), posBladeWall(vector::zero), posPeriodicWall(vector::zero);
    scalar bestBladePeriodicDist = GREAT, bestBladeWallDist = GREAT, bestPeriodicWallDist = GREAT;

    // v7.1 REAL FIX (SOL review): requiredPatches is NOT guaranteed to
    // be ordered [blade, periodic, wall] -- it is simply the set of
    // patches touching this corner, in whatever order the mesh
    // partitioner stored them. The original code assumed positional
    // order, causing a real mislabeling bug (a shroud+periodic_1 vertex
    // was printed as BLADE_PERIODIC because requiredPatches[0]
    // happened to be a non-blade patch ID that passed the hasBlade
    // check against the WRONG assumed role). Fixed by deriving each
    // role explicitly via caller-supplied role membership -- this
    // function now takes the actual blade/periodic/wall ID sets rather
    // than inferring them from position.
    label reqBlade = -1, reqPeriodic = -1, reqWall = -1;
    forAll(requiredPatches, ri)
    {
        const label p = requiredPatches[ri];
        if( bladeIds.found(p) && reqBlade < 0 ) reqBlade = p;
        else if( periodicIds.found(p) && reqPeriodic < 0 ) reqPeriodic = p;
        else if( wallIds.found(p) && reqWall < 0 ) reqWall = p;
    }

    forAll(surfPts, spI)
    {
        const scalar dSq = magSqr(surfPts[spI] - queryPoint);
        if( dSq > inspectRadiusSq ) continue;
        ++nNearby;

        const scalar dist = Foam::sqrt(dSq);

        // Raw topology: which patches touch this vertex via incident facets.
        labelHashSet incidentPatchSet;
        forAllRow(ptFacets, spI, fI)
        {
            const label triI = ptFacets(spI, fI);
            if( triI < 0 || triI >= facets.size() ) continue;
            incidentPatchSet.insert(facets[triI].region());
        }

        const label nEdgesHere = (spI < ptEdges.size()) ? ptEdges.sizeOfRow(spI) : 0;

        const bool hasBlade = reqBlade >= 0 && incidentPatchSet.found(reqBlade);
        const bool hasPeriodic = reqPeriodic >= 0 && incidentPatchSet.found(reqPeriodic);
        const bool hasWall = reqWall >= 0 && incidentPatchSet.found(reqWall);

        std::string classStr;
        if( hasBlade && hasPeriodic && hasWall )
        {
            classStr = "ALL_3"; ++nAll3;
            // v7.1 (SOL review): automatically probe the edge fan of
            // EVERY nearby ALL_3 vertex found, not just a hardcoded
            // ID from a prior run. Eliminates cross-run index-space
            // risk entirely and covers every control point (including
            // blade_2) in one run.
            probeCornerEdgeFan
            (
                surf, spI, requiredPatches, surfPts[spI],
                reqBlade, reqPeriodic, reqWall
            );

            // v7.1 Phase 2a prototype (phase2aStitchPrototype, hardcoded
            // to vertex 3748) REMOVED during code review -- confirmed
            // safe (read-only relative to surf) but superseded by the
            // generalized, role-based Phase 2c mechanism
            // (phase2cMultiJunctionStitch in cartesianMeshGenerator.C),
            // which now does this work for real, on any qualifying
            // vertex, not just a hardcoded one. Keeping it running was
            // pure redundant overhead on every single run.
        }
        else if( hasBlade && hasPeriodic )
        {
            classStr = "BLADE_PERIODIC"; ++nBladePeriodic;
            if( dist < bestBladePeriodicDist )
            { bestBladePeriodicDist = dist; posBladePeriodic = surfPts[spI]; haveBladePeriodic = true; }
        }
        else if( hasBlade && hasWall )
        {
            classStr = "BLADE_WALL"; ++nBladeWall;
            if( dist < bestBladeWallDist )
            { bestBladeWallDist = dist; posBladeWall = surfPts[spI]; haveBladeWall = true; }
        }
        else if( hasPeriodic && hasWall )
        {
            classStr = "PERIODIC_WALL"; ++nPeriodicWall;
            if( dist < bestPeriodicWallDist )
            { bestPeriodicWallDist = dist; posPeriodicWall = surfPts[spI]; havePeriodicWall = true; }
        }
        else if( hasBlade ) classStr = "BLADE_ONLY";
        else if( hasPeriodic ) classStr = "PERIODIC_ONLY";
        else if( hasWall ) classStr = "WALL_ONLY";
        else classStr = "OTHER";

        std::string patchListStr;
        label pCount = 0;
        forAllConstIter(labelHashSet, incidentPatchSet, it)
        {
            if( pCount++ > 0 ) patchListStr += " ";
            const label pr = it.key();
            patchListStr += (pr >= 0 && pr < label(patchNames.size())) ?
                std::string(patchNames[pr].c_str()) : std::string("?");
        }

        Info << "[TriSurfProbeVertex] queryBpI=" << queryBpI
             << " surfacePointId=" << spI
             << " distance=" << dist
             << " pos=" << surfPts[spI]
             << " nFacets=" << ptFacets.sizeOfRow(spI)
             << " nEdges=" << nEdgesHere
             << " patches=(" << patchListStr.c_str() << ")"
             << " class=" << classStr.c_str()
             << endl;
    }

    Info << "[TriSurfProbeSummary] bpI=" << queryBpI
         << " nearbyVertices=" << nNearby
         << " all3=" << nAll3
         << " bladePeriodic=" << nBladePeriodic
         << " bladeWall=" << nBladeWall
         << " periodicWall=" << nPeriodicWall
         << endl;

    if( haveBladePeriodic && haveBladeWall )
        Info << "[TriSurfProbePairwise] bpI=" << queryBpI
             << " |bladePeriodic-bladeWall|="
             << Foam::sqrt(magSqr(posBladePeriodic - posBladeWall)) << endl;
    if( haveBladePeriodic && havePeriodicWall )
        Info << "[TriSurfProbePairwise] bpI=" << queryBpI
             << " |bladePeriodic-periodicWall|="
             << Foam::sqrt(magSqr(posBladePeriodic - posPeriodicWall)) << endl;
    if( haveBladeWall && havePeriodicWall )
        Info << "[TriSurfProbePairwise] bpI=" << queryBpI
             << " |bladeWall-periodicWall|="
             << Foam::sqrt(magSqr(posBladeWall - posPeriodicWall)) << endl;
}

// v7.1 forensic diagnostic (SOL-specified): replicates
// meshOctree::findNearestCorner()'s EXACT edge-walking state machine
// for one specific vertex, logging every incident edge in pointEdges()
// order with its facet count, facet patch names, and action
// (NONMANIFOLD_WOULD_BREAK / SAME_REGION / PATCH_BOUNDARY), plus the
// accumulating qualifying-edge count and collected-patch set. Verified
// against actual source at meshOctreeFindNearestSurfacePoint.C:500-548
// before writing this -- same eFacets.sizeOfRow(eI)!=2 -> break logic,
// same region-difference test, same nEdges>2 threshold. Read-only,
// serial, no mesh mutation.
static void probeCornerEdgeFan
(
    const triSurf& surf,
    const label surfacePointId,
    const DynList<label>& requiredPatches,
    const point& expectedPos,
    const label reqBlade,
    const label reqPeriodic,
    const label reqWall
)
{
    const pointField& surfPts = surf.points();
    if( surfacePointId < 0 || surfacePointId >= surfPts.size() )
    {
        Info << "[TriSurfCornerEdgeProbe] ID_OUT_OF_RANGE surfacePointId="
             << surfacePointId << endl;
        return;
    }

    const point& actualPos = surfPts[surfacePointId];
    const scalar idCheckSq = magSqr(actualPos - expectedPos);
    if( idCheckSq > sqr(scalar(1e-10)) )
    {
        Info << "[TriSurfCornerEdgeProbe] ID_MISMATCH"
             << " surfacePointId=" << surfacePointId
             << " actualPos=" << actualPos
             << " expectedPos=" << expectedPos
             << " delta=" << Foam::sqrt(idCheckSq)
             << endl;
        return;
    }

    const VRWGraph& ptEdges = surf.pointEdges();
    const VRWGraph& edgeFacetsG = surf.edgeFacets();
    const LongList<edge>& edges = surf.edges();
    const LongList<labelledTri>& facets = surf.facets();
    const wordList& patchNames = surf.patchNames();

    auto patchName = [&](const label pr) -> std::string
    {
        return (pr >= 0 && pr < label(patchNames.size())) ?
            std::string(patchNames[pr].c_str()) : std::string("?");
    };

    const label totalIncidentEdges = ptEdges.sizeOfRow(surfacePointId);

    DynList<label> realNodePatches;
    label realQualifying(0);
    bool realBrokeEarly = false;
    label realBreakOrdinal = -1;
    label realBreakEdgeId = -1;
    label realBreakEdgeFacetCount = -1;

    DynList<label> skipNodePatches;
    label skipQualifying(0);

    // v7.1 (SOL): open (nFacets==1) edge tracking for pairwise comparison.
    DynList<label> openEdgeIds;
    DynList<label> openEdgeOtherPt;
    DynList<label> openEdgeOwningPatch;

    // v7.1 (SOL): tracks whether each of the 3 expected pairwise
    // interfaces already exists as a proper 2-facet cross-region edge.
    bool validBladePeriodic = false;
    bool validBladeWall = false;
    bool validPeriodicWall = false;

    Info << "[TriSurfCornerEdgeProbe] BEGIN surfacePointId=" << surfacePointId
         << " totalIncidentEdges=" << totalIncidentEdges << endl;

    forAllRow(ptEdges, surfacePointId, ord)
    {
        const label eI = ptEdges(surfacePointId, ord);
        const label nFacetsHere = edgeFacetsG.sizeOfRow(eI);

        std::string facetIdsStr, facetPatchesStr;
        forAllRow(edgeFacetsG, eI, fI)
        {
            if( fI > 0 ) { facetIdsStr += " "; facetPatchesStr += " "; }
            const label triI = edgeFacetsG(eI, fI);
            facetIdsStr += std::to_string(triI);
            if( triI >= 0 && triI < facets.size() )
                facetPatchesStr += patchName(facets[triI].region());
            else
                facetPatchesStr += "?";
        }

        std::string action;
        if( nFacetsHere != 2 )
        {
            action = "NONMANIFOLD_WOULD_BREAK";
        }
        else
        {
            const label r0 = facets[edgeFacetsG(eI, 0)].region();
            const label r1 = facets[edgeFacetsG(eI, 1)].region();
            action = (r0 != r1) ? "PATCH_BOUNDARY" : "SAME_REGION";
        }

        Info << "[TriSurfCornerEdgeProbe] surfacePointId=" << surfacePointId
             << " edgeOrdinal=" << ord
             << " edgeI=" << eI
             << " endpoints=(" << edges[eI].start() << " " << edges[eI].end() << ")"
             << " nFacets=" << nFacetsHere
             << " facetIds=(" << facetIdsStr.c_str() << ")"
             << " facetPatches=(" << facetPatchesStr.c_str() << ")"
             << " action=" << action.c_str()
             << endl;

        bool qualifiesCrossRegion = false;
        label rA = -1, rB = -1;
        if( nFacetsHere == 2 )
        {
            rA = facets[edgeFacetsG(eI, 0)].region();
            rB = facets[edgeFacetsG(eI, 1)].region();
            qualifiesCrossRegion = (rA != rB);

            // v7.1 (SOL): record which expected interface this proper
            // 2-facet edge represents, regardless of realBrokeEarly --
            // we want to know the FULL fan's interface completeness for
            // classification, not just what the real algorithm saw
            // before its break.
            if( qualifiesCrossRegion )
            {
                const bool hasBlade = (rA == reqBlade) || (rB == reqBlade);
                const bool hasPeriodic = (rA == reqPeriodic) || (rB == reqPeriodic);
                const bool hasWall = (rA == reqWall) || (rB == reqWall);
                if( hasBlade && hasPeriodic ) validBladePeriodic = true;
                if( hasBlade && hasWall ) validBladeWall = true;
                if( hasPeriodic && hasWall ) validPeriodicWall = true;
            }
        }

        if( realBrokeEarly == false )
        {
            if( nFacetsHere != 2 )
            {
                realBrokeEarly = true;
                realBreakOrdinal = ord;
                realBreakEdgeId = eI;
                realBreakEdgeFacetCount = nFacetsHere;
            }
            else if( qualifiesCrossRegion )
            {
                ++realQualifying;
                realNodePatches.appendIfNotIn(rA);
                realNodePatches.appendIfNotIn(rB);
            }
        }

        if( nFacetsHere == 2 && qualifiesCrossRegion )
        {
            ++skipQualifying;
            skipNodePatches.appendIfNotIn(rA);
            skipNodePatches.appendIfNotIn(rB);
        }

        // v7.1 (SOL): track open (nFacets==1) edges for the follow-up
        // pairwise-comparison test -- does the missing blade-wall
        // interface actually exist as two coincident-but-separate
        // dangling edge chains?
        if( nFacetsHere == 1 )
        {
            const label otherPt =
                (edges[eI].start() == surfacePointId) ?
                edges[eI].end() : edges[eI].start();
            const label owningPatch =
                (edgeFacetsG.sizeOfRow(eI) > 0) ?
                facets[edgeFacetsG(eI, 0)].region() : -1;
            openEdgeIds.append(eI);
            openEdgeOtherPt.append(otherPt);
            openEdgeOwningPatch.append(owningPatch);
        }
    }

    // v7.1 (SOL): pairwise comparison of every pair of open edges.
    // Tests the "split contact line" hypothesis: does the missing
    // blade-wall interface actually exist as two coincident-but-
    // separate dangling edge chains, rather than one shared edge?
    // No mesh mutation -- pure geometric measurement.
    // v7.1 (SOL): role-aware missing-interface pairing. Determine which
    // of the 3 expected pairwise interfaces (blade-periodic, blade-wall,
    // periodic-wall) already exist as proper 2-facet cross-region edges
    // (tracked via validBladePeriodic/validBladeWall/validPeriodicWall,
    // set during the main edge loop above whenever qualifiesCrossRegion
    // fired with the appropriate region pair). Only pair open edges
    // whose owning patches correspond to a MISSING interface -- this
    // prevents comparing unrelated open-edge pairs (e.g. two edges from
    // two separate, unrelated cracks at the same vertex), which caused
    // the spurious 80-100 degree CLASS_C results in the first sweep.
    label nBladeOpen = 0, nPeriodicOpen = 0, nWallOpen = 0;
    forAll(openEdgeOwningPatch, oi)
    {
        const label p = openEdgeOwningPatch[oi];
        if( p == reqBlade ) ++nBladeOpen;
        else if( p == reqPeriodic ) ++nPeriodicOpen;
        else if( p == reqWall ) ++nWallOpen;
    }

    auto tryPair = [&](const label roleA, const label roleB, const word& ifaceName)
    {
        forAll(openEdgeOwningPatch, ia)
        {
            if( openEdgeOwningPatch[ia] != roleA ) continue;
            forAll(openEdgeOwningPatch, ib)
            {
                if( openEdgeOwningPatch[ib] != roleB ) continue;

                const label ptA = openEdgeOtherPt[ia];
                const label ptB = openEdgeOtherPt[ib];
                if( ptA < 0 || ptA >= surfPts.size() ) continue;
                if( ptB < 0 || ptB >= surfPts.size() ) continue;

                const point& otherPosA = surfPts[ptA];
                const point& otherPosB = surfPts[ptB];
                const point& triplePos = surfPts[surfacePointId];

                const vector dirA = otherPosA - triplePos;
                const vector dirB = otherPosB - triplePos;
                const scalar lenA = Foam::sqrt(magSqr(dirA));
                const scalar lenB = Foam::sqrt(magSqr(dirB));

                const scalar otherEndpointSep =
                    Foam::sqrt(magSqr(otherPosA - otherPosB));

                scalar angleDeg = -1;
                if( lenA > VSMALL && lenB > VSMALL )
                {
                    const scalar cosAngle = (dirA & dirB) / (lenA * lenB);
                    const scalar clamped =
                        Foam::max(scalar(-1), Foam::min(scalar(1), cosAngle));
                    angleDeg = Foam::acos(clamped) * 180.0 / Foam::constant::mathematical::pi;
                }

                const scalar lengthRatio =
                    (lenB > VSMALL) ? (lenA / lenB) : scalar(-1);
                const scalar sepMm = otherEndpointSep * 1000.0;

                std::string prelimClass;
                if( sepMm < 0.5 && angleDeg < 0.1 && lengthRatio > 0.9 && lengthRatio < 1.111 )
                    prelimClass = "CLASS_A_NEAR_EXACT_SPLIT";
                else if( sepMm < 5.0 && angleDeg < 2.0 )
                    prelimClass = "CLASS_B_SEGMENTATION_MISMATCH";
                else
                    prelimClass = "CLASS_C_LIKELY_UNRELATED_OR_DIVERGENT";

                Info << "[TriSurfOpenEdgePairRoleAware] surfacePointId=" << surfacePointId
                     << " missingInterface=" << ifaceName.c_str()
                     << " edgeA=" << openEdgeIds[ia]
                     << " otherPtA=" << ptA
                     << " lengthA=" << lenA
                     << " edgeB=" << openEdgeIds[ib]
                     << " otherPtB=" << ptB
                     << " lengthB=" << lenB
                     << " separationMm=" << sepMm
                     << " angleDeg=" << angleDeg
                     << " lengthRatio=" << lengthRatio
                     << " prelimClass=" << prelimClass.c_str()
                     << endl;
            }
        }
    };

    if( !validBladePeriodic && nBladeOpen > 0 && nPeriodicOpen > 0 )
        tryPair(reqBlade, reqPeriodic, word("blade-periodic"));
    if( !validBladeWall && nBladeOpen > 0 && nWallOpen > 0 )
        tryPair(reqBlade, reqWall, word("blade-wall"));
    if( !validPeriodicWall && nPeriodicOpen > 0 && nWallOpen > 0 )
        tryPair(reqPeriodic, reqWall, word("periodic-wall"));

    // Vertex-level classification (SOL's six categories).
    label nMissingInterfaces = 0;
    if( !validBladePeriodic ) ++nMissingInterfaces;
    if( !validBladeWall ) ++nMissingInterfaces;
    if( !validPeriodicWall ) ++nMissingInterfaces;

    std::string vertexClass;
    if( nMissingInterfaces == 0 )
        vertexClass = "COMPLETE_MANIFOLD";
    else if( nMissingInterfaces == 1 )
        vertexClass = "SPLIT_ONE_INTERFACE";
    else if( nMissingInterfaces == 2 )
        vertexClass = "SPLIT_TWO_INTERFACES";
    else
        vertexClass = "SPLIT_THREE_INTERFACES";

    Info << "[TriSurfVertexClass] surfacePointId=" << surfacePointId
         << " validBladePeriodic=" << (validBladePeriodic ? 1 : 0)
         << " validBladeWall=" << (validBladeWall ? 1 : 0)
         << " validPeriodicWall=" << (validPeriodicWall ? 1 : 0)
         << " openBlade=" << nBladeOpen
         << " openPeriodic=" << nPeriodicOpen
         << " openWall=" << nWallOpen
         << " vertexClass=" << vertexClass.c_str()
         << endl;

    auto collectStr = [&](const DynList<label>& pats) -> std::string
    {
        std::string s;
        forAll(pats, i)
        {
            if( i > 0 ) s += " ";
            s += patchName(pats[i]);
        }
        return s;
    };

    auto countPresent = [&](const DynList<label>& pats) -> label
    {
        label n = 0;
        forAll(requiredPatches, i)
            if( pats.contains(requiredPatches[i]) ) ++n;
        return n;
    };

    const label realPresent = countPresent(realNodePatches);
    const label skipPresent = countPresent(skipNodePatches);
    const bool realEligible =
        (realQualifying > 2) && (realPresent >= requiredPatches.size());
    const bool skipEligible =
        (skipQualifying > 2) && (skipPresent >= requiredPatches.size());

    Info << "[TriSurfCornerEdgeProbe] SUMMARY surfacePointId=" << surfacePointId
         << " totalIncidentEdges=" << totalIncidentEdges
         << " realAlgorithm_brokeEarly=" << (realBrokeEarly ? 1 : 0)
         << " realAlgorithm_breakOrdinal=" << realBreakOrdinal
         << " realAlgorithm_breakEdgeId=" << realBreakEdgeId
         << " realAlgorithm_breakEdgeFacetCount=" << realBreakEdgeFacetCount
         << " realAlgorithm_qualifying=" << realQualifying
         << " realAlgorithm_collectedPatches=(" << collectStr(realNodePatches).c_str() << ")"
         << " realAlgorithm_requiredPresent=" << realPresent << "/" << requiredPatches.size()
         << " realAlgorithm_ACCEPT=" << (realEligible ? 1 : 0)
         << " skipMalformedSim_qualifying=" << skipQualifying
         << " skipMalformedSim_collectedPatches=(" << collectStr(skipNodePatches).c_str() << ")"
         << " skipMalformedSim_requiredPresent=" << skipPresent << "/" << requiredPatches.size()
         << " skipMalformedSim_wouldACCEPT=" << (skipEligible ? 1 : 0)
         << endl;
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

    // Flag for non-corner point detected inside OMP loop.
    // Cannot call FatalError inside parallel region (UB) -- check after.
    label foundNonCorner = 0;

    // v7: census storage. Pre-sized before the parallel loop, each
    // thread writes only to its own unique cornerI index -- same
    // race-free pattern established and confirmed safe in v1-v3.
    struct TJCensusRecord
    {
        bool valid;
        label bpI;
        point pos;
        bool isTarget;
        DynList<word> patchNames;
        TJCensusRecord() : valid(false), bpI(-1), pos(vector::zero), isTarget(false) {}
    };
    List<TJCensusRecord> censusRecords;
    if( tripleJunctionCornerFixEnabled_ )
        censusRecords.setSize(nodesToMap.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(nodesToMap, cornerI)
    {
        const label bpI = nodesToMap[cornerI];
        // Skip BL/no-BL interface points - must stay on feature curve
        if( !protectedPoints_.empty() && protectedPoints_.found(bpI) )
            continue;
        // Skip BL/neutral points - must stay on BL-side patch
        if( !blNeutralPoints_.empty() && blNeutralPoints_.found(bpI) )
            continue;
        if( !corners.found(bpI) )
        {
            // FatalError inside OMP is UB -- set flag, report after loop.
            # ifdef USE_OMP
            # pragma omp atomic write
            # endif
            foundNonCorner = 1;
            continue;
        }

        const point& p = points[bPoints[bpI]];
        const scalar maxDist = mappingDistance[cornerI];

        //- find the nearest position to the given point patches
        const DynList<label> patches = pPatches[bpI];
        if( patches.size() == 0 ) continue;  // guard: no patches = no valid snap

        // v7: narrow blade+periodic+hub/shroud fix. Generic 3+ patch
        // corners are still frozen (unchanged production behavior) --
        // only the EXACT target junction type gets the fix, and only
        // when tripleJunctionCornerFixEnabled_ is explicitly true.
        //
        // CENSUS (v7): runs inside the OMP loop, read-only (patch name
        // lookup only, same safe call used since v1 selectiveWeld code),
        // stores into per-thread-indexed arrays. Answers SOL Q1: does
        // the blade+periodic+hub/shroud junction even exist among the
        // mapper 3+ patch corners at this pipeline stage?
        if( patches.size() >= 3 )
        {
            if( tripleJunctionCornerFixEnabled_ )
            {
                label nBlade = 0;
                label nPeriodic = 0;
                label nHub = 0;
                label nShroud = 0;
                label nFlow = 0;
                forAll(patches, pI)
                {
                    const label patchI = patches[pI];
                    nBlade    += tjBladePatchIds_.found(patchI);
                    nPeriodic += tjPeriodicPatchIds_.found(patchI);
                    nHub      += tjHubPatchIds_.found(patchI);
                    nShroud   += tjShroudPatchIds_.found(patchI);
                    nFlow     += tjInletOutletPatchIds_.found(patchI);
                }
                const bool isTargetTJ =
                    patches.size() == 3
                 && nBlade == 1
                 && nPeriodic == 1
                 && (nHub + nShroud) == 1
                 && nFlow == 0;

                censusRecords[cornerI].valid = true;
                censusRecords[cornerI].bpI = bpI;
                censusRecords[cornerI].pos = p;
                censusRecords[cornerI].isTarget = isTargetTJ;
                const wordList& allNames = meshOctree_.surface().patchNames();
                forAll(patches, pI)
                {
                    const label patchI = patches[pI];
                    if( patchI >= 0 && patchI < label(allNames.size()) )
                        censusRecords[cornerI].patchNames.append(allNames[patchI]);
                    else
                        censusRecords[cornerI].patchNames.append(word("UNKNOWN"));
                }
            }
            continue;
        }

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

    if( foundNonCorner )
    {
        FatalErrorIn
        (
            "meshSurfaceMapper::mapCorners(const labelLongList&)"
        ) << "Trying to map a point that is not a corner"
            << abort(FatalError);
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
                // Tolerance based on face area magnitude -- dimensionally
                // consistent with h = fn & (fc - cc) which has area*length units.
                // Tolerance: volume-scale hybrid.
                // h = fn & (fc - cc) has area*length scale.
                // face area * local length as main tolerance,
                // local-volume floor for very small faces.
                const scalar localLen =
                    Foam::sqrt(mappingDistance[i] + VSMALL);
                const scalar faceScale = Foam::max
                (
                    Foam::mag(fn) * localLen * scalar(1e-6),
                    Foam::pow(localLen, 3) * scalar(1e-12)
                );
                if( h <= -faceScale )
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

    // v7: SERIAL pass. tripleJunctionCornerFixEnabled_ gates everything
    // below -- unset call sites see zero behavior change (SOL v6 lesson).
    if( tripleJunctionCornerFixEnabled_ )
    {
        // Census (SOL request): one line per generic 3+ patch corner,
        // reporting real patch names, BEFORE any filtering/fix logic.
        label nGeneric = 0;
        label nTarget = 0;
        forAll(censusRecords, i)
        {
            const TJCensusRecord& rec = censusRecords[i];
            if( rec.valid == false ) continue;
            ++nGeneric;
            if( rec.isTarget ) ++nTarget;

            std::string nameStr;
            forAll(rec.patchNames, pI)
            {
                if( pI > 0 ) nameStr += " ";
                nameStr += std::string(rec.patchNames[pI].c_str());
            }

            Info << "[TripleJunctionCensus] " << tripleJunctionDiagTag_.c_str()
                 << " bpI=" << rec.bpI
                 << " patches=(" << nameStr.c_str() << ")"
                 << " target=" << (rec.isTarget ? 1 : 0)
                 << endl;
        }
        Info << "[TripleJunctionCensus] SUMMARY " << tripleJunctionDiagTag_.c_str()
             << ": generic3plus=" << nGeneric
             << " targetBladePeriodicHubShroud=" << nTarget
             << endl;

        // v7 FIX PASS: only for target=1 corners. Same mechanism as v6
        // (findNearestCorner + mappingDistance-bounded gate + damped
        // move + pyramid validity check), narrowed to the exact
        // blade+periodic+hub/shroud junction type only. On acceptance,
        // record the corner's immediate mesh neighbors (via
        // pointPoints(), real mesh connectivity, not triSurf) so
        // mapEdgeNodes() can report what production logic does to them.
        tjWatchedBoundaryPointToSourceCorner_.clear();
        tjEdgeWatchArmed_ = false;
        const VRWGraph& pPoints = surfaceEngine_.pointPoints();
        label nWatchOverlaps = 0;
        const VRWGraph& pFaces = surfaceEngine_.pointFaces();
        const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();
        const labelList& faceOwners = surfaceEngine_.faceOwners();
        const pointFieldPMG& pts = surfaceEngine_.points();
        const cellListPMG& cells = surfaceEngine_.mesh().cells();
        const faceListPMG& allFaces = surfaceEngine_.mesh().faces();

        label nTargetMoved = 0;
        forAll(censusRecords, i)
        {
            const TJCensusRecord& rec = censusRecords[i];
            if( rec.valid == false ) continue;
            if( rec.isTarget == false ) continue;

            const label bpI = rec.bpI;
            const point p = rec.pos;
            const DynList<label> cPatches = pPatches[bpI];

            const scalar mappingDistanceSq = mappingDistance[i];

            // Code review cleanup: this forensic sweep (probeTriSurfJunction,
            // a brute-force scan of every surface point per target corner)
            // was built for one-time investigation and is expensive. Now
            // independently gated by tripleJunctionDeepDiagnosticEnabled_
            // so re-enabling corner-snapping does not silently re-enable
            // this scan on every run.
            if( tripleJunctionDeepDiagnosticEnabled_ )
            {
                labelHashSet wallIds = tjHubPatchIds_;
                forAllConstIter(labelHashSet, tjShroudPatchIds_, wIt)
                    wallIds.insert(wIt.key());
                probeTriSurfJunction
                (
                    meshOctree_.surface(),
                    p,
                    mappingDistanceSq,
                    bpI,
                    cPatches,
                    tjBladePatchIds_,
                    tjPeriodicPatchIds_,
                    wallIds
                );
            }

            point cornerP(p);
            scalar cornerDistSq;
            label surfaceCornerId;
            const bool cornerOk = meshOctree_.findNearestCorner
            (
                cornerP, cornerDistSq, surfaceCornerId, p, cPatches
            );

            // v7.1 (SOL): diagnostic-only terminal-outcome classification
            // for EVERY target corner. Zero behavior change -- logs the
            // branch the EXISTING logic already takes.
            std::string v71NameStr;
            forAll(rec.patchNames, pI)
            {
                if( pI > 0 ) v71NameStr += " ";
                v71NameStr += std::string(rec.patchNames[pI].c_str());
            }

            if( cornerOk == false )
            {
                Info << "[TripleJunctionFixV71] bpI=" << bpI
                     << " meshPointI=" << bPoints[bpI]
                     << " names=(" << v71NameStr.c_str() << ")"
                     << " from=" << p
                     << " mappingDistance=" << mappingDistanceSq
                     << " result=TRUE_CORNER_NOT_FOUND"
                     << endl;
                continue;
            }
            if( cornerDistSq > mappingDistanceSq )
            {
                Info << "[TripleJunctionFixV71] bpI=" << bpI
                     << " meshPointI=" << bPoints[bpI]
                     << " names=(" << v71NameStr.c_str() << ")"
                     << " from=" << p
                     << " target=" << cornerP
                     << " targetDist=" << Foam::sqrt(cornerDistSq)
                     << " mappingDistance=" << Foam::sqrt(mappingDistanceSq)
                     << " surfaceCornerId=" << surfaceCornerId
                     << " result=TOO_FAR"
                     << endl;
                continue;
            }

            const vector delta = cornerP - p;
            const point relaxedPoint = p + cornerSnapRelaxation_ * delta;

            sMod.moveBoundaryVertexNoUpdate(bpI, relaxedPoint);

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
                const scalar localLen = Foam::sqrt(mappingDistanceSq + VSMALL);
                const scalar faceScale = Foam::max
                (
                    Foam::mag(fn) * localLen * scalar(1e-6),
                    Foam::pow(localLen, 3) * scalar(1e-12)
                );
                if( h <= -faceScale )
                { validMove = false; break; }
            }

            if( validMove == false )
            {
                sMod.moveBoundaryVertexNoUpdate(bpI, p);
                Info << "[TripleJunctionFixV71] bpI=" << bpI
                     << " meshPointI=" << bPoints[bpI]
                     << " names=(" << v71NameStr.c_str() << ")"
                     << " from=" << p
                     << " target=" << cornerP
                     << " targetDist=" << Foam::sqrt(cornerDistSq)
                     << " appliedDist="
                     << Foam::sqrt(Foam::max(magSqr(relaxedPoint - p), scalar(0)))
                     << " surfaceCornerId=" << surfaceCornerId
                     << " result=PYRAMID_REVERT"
                     << endl;
                continue;
            }

            ++nTargetMoved;
            const vector appliedVec = relaxedPoint - p;
            const scalar appliedDisp = Foam::sqrt(Foam::max(magSqr(appliedVec), scalar(0)));

            Info << "[TripleJunctionFixV71] bpI=" << bpI
                 << " meshPointI=" << bPoints[bpI]
                 << " names=(" << v71NameStr.c_str() << ")"
                 << " from=" << p
                 << " target=" << cornerP
                 << " targetDist=" << Foam::sqrt(cornerDistSq)
                 << " appliedDist=" << appliedDisp
                 << " surfaceCornerId=" << surfaceCornerId
                 << " result=ACCEPT"
                 << endl;

            Info << "[TripleJunctionFixV7] ACCEPT"
                 << " bpI=" << bpI
                 << " meshPointI=" << bPoints[bpI]
                 << " patches=" << cPatches
                 << " from=" << p
                 << " target=" << cornerP
                 << " applied=" << appliedVec
                 << " appliedDist=" << appliedDisp
                 << endl;

            // Record immediate mesh neighbors for mapEdgeNodes() to
            // report on. pointPoints() gives real mesh connectivity.
            forAllRow(pPoints, bpI, nI)
            {
                const label neighborBpI = pPoints(bpI, nI);
                if( tjWatchedBoundaryPointToSourceCorner_.found(neighborBpI) )
                {
                    if( tjWatchedBoundaryPointToSourceCorner_[neighborBpI] != bpI )
                    {
                        Info << "[TripleJunctionEdgeWatch] OVERLAP"
                             << " edgeBpI=" << neighborBpI
                             << " firstCorner="
                             << tjWatchedBoundaryPointToSourceCorner_[neighborBpI]
                             << " secondCorner=" << bpI << endl;
                        ++nWatchOverlaps;
                    }
                }
                else
                {
                    tjWatchedBoundaryPointToSourceCorner_.insert(neighborBpI, bpI);
                }
            }
        }
        if( tjWatchedBoundaryPointToSourceCorner_.size() > 0 )
            tjEdgeWatchArmed_ = true;
        Info << "[TripleJunctionFixV7] SUMMARY " << tripleJunctionDiagTag_.c_str()
             << ": targetCorners=" << nTarget
             << " targetMoved=" << nTargetMoved
             << " watchedNeighbors=" << tjWatchedBoundaryPointToSourceCorner_.size()
             << " watchOverlaps=" << nWatchOverlaps
             << endl;
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

    // v7: watch-list storage for the incident-edge-neighbor
    // instrumentation. tjWatchedBoundaryPointToSourceCorner_ was populated
    // by the fix pass in mapCorners() (same mapper object, same
    // mapCornersAndEdges() call). POD-only storage (no DynList/string/
    // I/O) so this is safe to write from inside the OMP loop below --
    // per SOL: instrument the ACTUAL production calculation, no new
    // octree queries, no additional measurement-perturbation risk.
    boolList tjIsWatched(nodesToMap.size(), false);
    labelList tjSourceCornerBpI(nodesToMap.size(), label(-1));
    pointField tjApproxPos(nodesToMap.size(), point::zero);
    scalarList tjDistSqApprox(nodesToMap.size(), scalar(-1));
    pointField tjEdgePos(nodesToMap.size(), point::zero);
    scalarList tjEdgeDistSq(nodesToMap.size(), scalar(-1));
    // v7 decision codes (per SOL review -- must cover every production
    // exit path, so a leftover TJ_UNSEEN at output time is itself a bug
    // signal meaning we missed a path):
    //   -1 = TJ_UNSEEN (should never appear in final output)
    //    0 = EDGE_ACCEPTED
    //    1 = RATIO_REJECTED  (distSq > 1.2*distSqApprox)
    //    2 = NO_PATCHES
    //    3 = PROTECTED
    //    4 = EDGE_NOT_FOUND  (findNearestEdgePoint returned false --
    //        production still proceeds using mapPointApprox-derived
    //        mapPoint/distSq, but the true feature-edge search failed)
    labelList tjDecision(nodesToMap.size(), label(-1));
    const bool watchTJEdges =
        tripleJunctionCornerFixEnabled_
     && tjEdgeWatchArmed_
     && !tjWatchedBoundaryPointToSourceCorner_.empty();
    if( watchTJEdges )
    {
        forAll(nodesToMap, i)
        {
            const label bpI = nodesToMap[i];
            if( tjWatchedBoundaryPointToSourceCorner_.found(bpI) )
            {
                tjIsWatched[i] = true;
                tjSourceCornerBpI[i] = tjWatchedBoundaryPointToSourceCorner_[bpI];
            }
        }
    }

    //- map point to the nearest vertex on the surface mesh
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(nodesToMap, i)
    {
        // Skip BL/no-BL interface points - must stay on feature curve
        if( isProtected[i] )
        {
            if( tjIsWatched[i] ) tjDecision[i] = 3;
            continue;
        }
        const label bpI = nodesToMap[i];
        const point& p = points[bPoints[bpI]];

        //- find patches at this edge point
        const DynList<label> patches = pPatches[bpI];
        if( patches.size() == 0 )
        {
            if( tjIsWatched[i] ) tjDecision[i] = 2;
            continue;  // guard: no patches = no valid snap
        }

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
        const bool edgeFoundOk =
            meshOctree_.findNearestEdgePoint(mapPoint, distSq, nse, p, patches);

        // v7 tap: EDGE_NOT_FOUND path (SOL review -- previously
        // uninstrumented; production silently proceeds using
        // mapPointApprox-derived values when this search fails, which
        // our earlier decision codes could not distinguish from a
        // genuine acceptance).
        if( tjIsWatched[i] && edgeFoundOk == false )
        {
            tjApproxPos[i] = mapPointApprox;
            tjDistSqApprox[i] = distSqApprox;
            tjEdgePos[i] = point::zero;
            tjEdgeDistSq[i] = scalar(-1);
            tjDecision[i] = 4;
        }

        //- use the vertex with the smallest mapping distance
        // Do NOT fall back to mapPointApprox if findNearestEdgePoint fails.
        // mapPointApprox is an average of patch projections -- not guaranteed
        // to lie on the feature curve. Near blade/hub/periodic corners this
        // creates protrusions. If true edge projection fails, skip this point.
        if( distSq > 1.2 * distSqApprox )
        {
            // v7 tap: capture production's RATIO_REJECTED outcome for
            // watched points, using values already computed above
            // (mapPointApprox, distSqApprox, mapPoint, distSq). No new
            // octree queries -- exactly the values production itself
            // just calculated.
            if( tjIsWatched[i] )
            {
                tjApproxPos[i] = mapPointApprox;
                tjDistSqApprox[i] = distSqApprox;
                tjEdgePos[i] = mapPoint;
                tjEdgeDistSq[i] = distSq;
                tjDecision[i] = 1;
            }
            projectedPoints[i] = point::zero;
            projectedDistSq[i] = scalar(-1);
            continue;
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

        // v7 tap: capture production's EDGE_ACCEPTED outcome for
        // watched points. Same no-new-queries principle.
        if( tjIsWatched[i] )
        {
            tjApproxPos[i] = mapPointApprox;
            tjDistSqApprox[i] = distSqApprox;
            tjEdgePos[i] = mapPoint;
            tjEdgeDistSq[i] = distSq;
            tjDecision[i] = 0;
        }
    }

    // Serial move pass - deterministic order, no validity check
    // Edge nodes must reach feature curve - validity check inappropriate
    forAll(nodesToMap, i)
    {
        const label bpI = nodesToMap[i];
        if( projectedDistSq[i] < scalar(0) ) continue;
        sMod.moveBoundaryVertexNoUpdate(bpI, projectedPoints[i]);
    }

    // v7: serial dump of watched incident-edge-neighbor outcomes.
    // Uses ONLY values already captured from production computation
    // above -- no new octree queries, per SOL requirement.
    if( tripleJunctionCornerFixEnabled_ )
    {
        label nWatched = 0;
        forAll(nodesToMap, i)
            if( tjIsWatched[i] ) ++nWatched;

        if( nWatched > 0 )
        {
            const char* decisionNames[6] =
            {
                "EDGE_ACCEPTED", "RATIO_REJECTED", "NO_PATCHES",
                "PROTECTED", "EDGE_NOT_FOUND", "TJ_UNSEEN"
            };

            label nAccepted = 0, nRatioRejected = 0, nNoPatches = 0;
            label nProtected = 0, nEdgeNotFound = 0, nUnseen = 0;

            forAll(nodesToMap, i)
            {
                if( !tjIsWatched[i] ) continue;
                const label bpI = nodesToMap[i];
                const label rawDecision = tjDecision[i];
                const label decisionIdx =
                    (rawDecision >= 0 && rawDecision <= 4) ? rawDecision : 5;

                switch(decisionIdx)
                {
                    case 0: ++nAccepted; break;
                    case 1: ++nRatioRejected; break;
                    case 2: ++nNoPatches; break;
                    case 3: ++nProtected; break;
                    case 4: ++nEdgeNotFound; break;
                    default: ++nUnseen; break;
                }

                const point& curPos = points[bPoints[bpI]];
                const scalar appliedDispSq =
                    (decisionIdx == 0) ?
                    magSqr(projectedPoints[i] - curPos) : scalar(-1);

                Info << "[TripleJunctionEdgeNeighbor] " << tripleJunctionDiagTag_.c_str()
                     << " cornerBpI=" << tjSourceCornerBpI[i]
                     << " edgeBpI=" << bpI
                     << " edgeMeshPointI=" << bPoints[bpI]
                     << " currentPos=" << curPos
                     << " mapPointApprox=" << tjApproxPos[i]
                     << " distSqApprox=" << tjDistSqApprox[i]
                     << " edgeTarget=" << tjEdgePos[i]
                     << " edgeDistSq=" << tjEdgeDistSq[i]
                     << " decision=" << decisionNames[decisionIdx]
                     << " appliedDist="
                     << (appliedDispSq >= 0 ? Foam::sqrt(appliedDispSq) : scalar(-1))
                     << endl;

                if( decisionIdx == 5 )
                    Info << "[TripleJunctionEdgeNeighbor] WARNING: bpI=" << bpI
                         << " reached output with TJ_UNSEEN -- missed a"
                         << " production exit path in this instrumentation"
                         << endl;
            }

            Info << "[TripleJunctionEdgeNeighbor] SUMMARY " << tripleJunctionDiagTag_.c_str()
                 << ": watchedEdgeNodes=" << nWatched
                 << " accepted=" << nAccepted
                 << " ratioRejected=" << nRatioRejected
                 << " edgeNotFound=" << nEdgeNotFound
                 << " noPatches=" << nNoPatches
                 << " protectedPts=" << nProtected
                 << " unseen=" << nUnseen
                 << endl;
        }
        else
        {
            Info << "[TripleJunctionEdgeNeighbor] SUMMARY " << tripleJunctionDiagTag_.c_str()
                 << ": watchedEdgeNodes=0 (no watched points this call)" << endl;
        }

        // v7: one-shot consumption (SOL review) -- clear so a future
        // mapEdgeNodes() call on this mapper without a preceding armed
        // mapCorners() call sees an empty/unarmed watch set, not stale
        // attribution from this call.
        tjWatchedBoundaryPointToSourceCorner_.clear();
        tjEdgeWatchArmed_ = false;
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
    // Old edge validity check removed -- edge nodes must reach
    // feature curves; hard revert inappropriate here.
    // If bad pyramids persist near feature edges, add bisection.

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
