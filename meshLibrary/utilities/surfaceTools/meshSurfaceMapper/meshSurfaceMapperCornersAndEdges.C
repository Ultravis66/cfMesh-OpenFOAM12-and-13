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
#include "sourceFeatureGraph.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceEngineModifier.H"
#include "polyMeshGenAddressing.H"
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
    boolList mappingAccepted;
    mapEdgeNodes(nodesToMap, mappingAccepted);
}


void meshSurfaceMapper::mapEdgeNodes
(
    const labelLongList& nodesToMap,
    boolList& mappingAccepted
)
{
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();

    pointField desiredPositions
    (
        nodesToMap.size(),
        point::zero
    );

    forAll(nodesToMap, i)
    {
        desiredPositions[i] =
            points[bPoints[nodesToMap[i]]];
    }

    mapEdgeNodes
    (
        nodesToMap,
        desiredPositions,
        mappingAccepted
    );
}


void meshSurfaceMapper::mapEdgeNodes
(
    const labelLongList& nodesToMap,
    const pointField& desiredPositions,
    boolList& mappingAccepted
)
{
    mappingAccepted.setSize(nodesToMap.size());
    mappingAccepted = false;

    if( desiredPositions.size() != nodesToMap.size() )
    {
        WarningIn
        (
            "void meshSurfaceMapper::mapEdgeNodes"
            "(const labelLongList&, const pointField&, boolList&)"
        )
            << "desiredPositions size " << desiredPositions.size()
            << " does not match nodesToMap size "
            << nodesToMap.size()
            << ". No points will be moved."
            << endl;

        return;
    }

    const pointFieldPMG& points = surfaceEngine_.points();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();

    const meshSurfacePartitioner& mPart = meshPartitioner();
    const VRWGraph& pPatches = mPart.pointPatches();

    // Diagnostic only: census the SOURCE triSurface representation of
    // the physical patch pairs entering this mapper call.
    //
    // A native source feature requires exactly two incident facets with
    // different region IDs. Open nFacets==1 edges are deliberately NOT
    // native cross-region features and are invisible to the current
    // findNearestEdgePoint() feature search.
    static bool triSurfPairCensusPrinted(false);

    if( !triSurfPairCensusPrinted )
    {
        triSurfPairCensusPrinted = true;

        const triSurf& censusSurf = meshOctree_.surface();
        const VRWGraph& censusEdgeFacets = censusSurf.edgeFacets();
        const LongList<labelledTri>& censusFacets = censusSurf.facets();
        const LongList<edge>& censusEdges = censusSurf.edges();
        const pointField& censusPoints = censusSurf.points();
        const wordList& censusPatchNames = censusSurf.patchNames();

        const label nCensusPatches = censusPatchNames.size();

        labelList sharedCrossRegion
        (
            nCensusPatches*nCensusPatches,
            0
        );

        labelList openEdgesByPatch(nCensusPatches, 0);

        // Keep the source IDs/owners of all true open edges so we can
        // measure whether separate patch-boundary chains geometrically
        // represent the same physical seam.
        DynList<label> censusOpenEdgeIds;
        DynList<label> censusOpenEdgePatch;

        label nNonManifoldEdges = 0;

        forAll(censusEdgeFacets, edgeI)
        {
            const label nEF = censusEdgeFacets.sizeOfRow(edgeI);

            if( nEF == 1 )
            {
                const label f0 = censusEdgeFacets(edgeI, 0);

                if( f0 >= 0 && f0 < censusFacets.size() )
                {
                    const label r0 = censusFacets[f0].region();

                    if( r0 >= 0 && r0 < nCensusPatches )
                    {
                        ++openEdgesByPatch[r0];
                        censusOpenEdgeIds.append(edgeI);
                        censusOpenEdgePatch.append(r0);
                    }
                }

                continue;
            }

            if( nEF != 2 )
            {
                if( nEF > 2 )
                    ++nNonManifoldEdges;

                continue;
            }

            const label f0 = censusEdgeFacets(edgeI, 0);
            const label f1 = censusEdgeFacets(edgeI, 1);

            if
            (
                f0 < 0 || f1 < 0 ||
                f0 >= censusFacets.size() ||
                f1 >= censusFacets.size()
            )
            {
                continue;
            }

            label r0 = censusFacets[f0].region();
            label r1 = censusFacets[f1].region();

            if( r0 == r1 )
                continue;

            if
            (
                r0 < 0 || r1 < 0 ||
                r0 >= nCensusPatches ||
                r1 >= nCensusPatches
            )
            {
                continue;
            }

            if( r1 < r0 )
            {
                const label tmp = r0;
                r0 = r1;
                r1 = tmp;
            }

            ++sharedCrossRegion[r0*nCensusPatches + r1];
        }

        // Which patch pairs are actually represented by the volume
        // feature points in this mapper call?
        labelList usedVolumePair
        (
            nCensusPatches*nCensusPatches,
            0
        );

        forAll(nodesToMap, i)
        {
            const label bpI = nodesToMap[i];

            if( pPatches.sizeOfRow(bpI) != 2 )
                continue;

            label p0 = pPatches(bpI, 0);
            label p1 = pPatches(bpI, 1);

            if
            (
                p0 < 0 || p1 < 0 ||
                p0 >= nCensusPatches ||
                p1 >= nCensusPatches
            )
            {
                continue;
            }

            if( p1 < p0 )
            {
                const label tmp = p0;
                p0 = p1;
                p1 = tmp;
            }

            ++usedVolumePair[p0*nCensusPatches + p1];
        }

        for(label p0=0; p0<nCensusPatches; ++p0)
        {
            for(label p1=p0+1; p1<nCensusPatches; ++p1)
            {
                const label pairI =
                    p0*nCensusPatches + p1;

                if( usedVolumePair[pairI] == 0 )
                    continue;

                Info
                    << "[TRISURF_PAIR_CENSUS]"
                    << " pair="
                    << censusPatchNames[p0] << "|"
                    << censusPatchNames[p1]
                    << " volumeFeaturePoints="
                    << usedVolumePair[pairI]
                    << " nativeSharedEdges="
                    << sharedCrossRegion[pairI]
                    << endl;
            }
        }

        for(label patchI=0; patchI<nCensusPatches; ++patchI)
        {
            if( openEdgesByPatch[patchI] == 0 )
                continue;

            Info
                << "[TRISURF_OPEN_CENSUS]"
                << " patch=" << censusPatchNames[patchI]
                << " openEdges=" << openEdgesByPatch[patchI]
                << endl;
        }

        // Diagnostic-only geometric pairing census for open source edges.
        //
        // For every volume-used A/B patch pair having open source edges
        // on both sides, query BOTH directions (A->B and B->A). Each open
        // edge selects the opposite-patch open edge having the smallest
        // normalized midpoint-to-segment separation.
        //
        // Report:
        //   sepNorm      : transverse separation / min(local edge lengths)
        //   align        : |unit tangent dot product|
        //   overlapFrac  : projected overlap / min(local edge lengths)
        //
        // This does NOT create virtual features and does NOT move the mesh.
        for(label p0=0; p0<nCensusPatches; ++p0)
        {
            for(label p1=p0+1; p1<nCensusPatches; ++p1)
            {
                const label pairI = p0*nCensusPatches + p1;

                if( usedVolumePair[pairI] == 0 )
                    continue;

                if
                (
                    openEdgesByPatch[p0] == 0 ||
                    openEdgesByPatch[p1] == 0
                )
                {
                    continue;
                }

                label nQueries = 0;
                label nFound = 0;

                scalar sumBestSepNorm = 0;
                scalar minBestSepNorm = GREAT;
                scalar maxBestSepNorm = 0;

                scalar sumBestAlign = 0;
                scalar minBestAlign = GREAT;

                scalar sumBestOverlap = 0;
                scalar minBestOverlap = GREAT;

                label nNearExact = 0;
                label nStrong = 0;

                forAll(censusOpenEdgeIds, oi)
                {
                    const label ownerA = censusOpenEdgePatch[oi];

                    if( ownerA != p0 && ownerA != p1 )
                        continue;

                    const label ownerB =
                        (ownerA == p0) ? p1 : p0;

                    const label eAI = censusOpenEdgeIds[oi];
                    const edge& eA = censusEdges[eAI];

                    const point& a0 = censusPoints[eA.start()];
                    const point& a1 = censusPoints[eA.end()];

                    const vector dA = a1 - a0;
                    const scalar lenA = Foam::sqrt(magSqr(dA));

                    if( lenA < ROOTVSMALL )
                        continue;

                    ++nQueries;

                    scalar bestSepNorm = GREAT;
                    scalar bestAlign = 0;
                    scalar bestOverlap = 0;
                    bool foundBest = false;

                    forAll(censusOpenEdgeIds, oj)
                    {
                        if( censusOpenEdgePatch[oj] != ownerB )
                            continue;

                        const label eBI = censusOpenEdgeIds[oj];
                        const edge& eB = censusEdges[eBI];

                        const point& b0 = censusPoints[eB.start()];
                        const point& b1 = censusPoints[eB.end()];

                        const vector dB = b1 - b0;
                        const scalar lenB = Foam::sqrt(magSqr(dB));

                        if( lenB < ROOTVSMALL )
                            continue;

                        const scalar minLen =
                            Foam::min(lenA, lenB);

                        if( minLen < ROOTVSMALL )
                            continue;

                        const scalar align =
                            Foam::min
                            (
                                scalar(1),
                                Foam::mag(dA & dB) /
                                (lenA*lenB + VSMALL)
                            );

                        const point midA = 0.5*(a0 + a1);
                        const point midB = 0.5*(b0 + b1);

                        const point qAonB =
                            help::nearestPointOnTheEdgeExact
                            (
                                b0, b1, midA
                            );

                        const point qBonA =
                            help::nearestPointOnTheEdgeExact
                            (
                                a0, a1, midB
                            );

                        const scalar sep =
                            Foam::min
                            (
                                Foam::sqrt(magSqr(midA - qAonB)),
                                Foam::sqrt(magSqr(midB - qBonA))
                            );

                        const scalar sepNorm =
                            sep / (minLen + VSMALL);

                        // Project B endpoints onto the infinite A line,
                        // intersect that parameter interval with A=[0,1].
                        const scalar lenASq = magSqr(dA);
                        const scalar tB0A =
                            ((b0-a0) & dA) / (lenASq + VSMALL);
                        const scalar tB1A =
                            ((b1-a0) & dA) / (lenASq + VSMALL);

                        const scalar loA =
                            Foam::max
                            (
                                scalar(0),
                                Foam::min(tB0A, tB1A)
                            );

                        const scalar hiA =
                            Foam::min
                            (
                                scalar(1),
                                Foam::max(tB0A, tB1A)
                            );

                        const scalar overlapA =
                            Foam::max(scalar(0), hiA-loA)*lenA;

                        // Symmetric projection: A endpoints onto B.
                        const scalar lenBSq = magSqr(dB);
                        const scalar tA0B =
                            ((a0-b0) & dB) / (lenBSq + VSMALL);
                        const scalar tA1B =
                            ((a1-b0) & dB) / (lenBSq + VSMALL);

                        const scalar loB =
                            Foam::max
                            (
                                scalar(0),
                                Foam::min(tA0B, tA1B)
                            );

                        const scalar hiB =
                            Foam::min
                            (
                                scalar(1),
                                Foam::max(tA0B, tA1B)
                            );

                        const scalar overlapB =
                            Foam::max(scalar(0), hiB-loB)*lenB;

                        const scalar overlapFrac =
                            Foam::min
                            (
                                scalar(1),
                                0.5*(overlapA + overlapB) /
                                (minLen + VSMALL)
                            );

                        // Primary selection is geometric separation.
                        // For effectively tied separation, prefer better
                        // tangent alignment.
                        if
                        (
                            !foundBest ||
                            sepNorm < bestSepNorm ||
                            (
                                Foam::mag(sepNorm-bestSepNorm) < 1e-12 &&
                                align > bestAlign
                            )
                        )
                        {
                            foundBest = true;
                            bestSepNorm = sepNorm;
                            bestAlign = align;
                            bestOverlap = overlapFrac;
                        }
                    }

                    if( !foundBest )
                        continue;

                    ++nFound;

                    sumBestSepNorm += bestSepNorm;
                    minBestSepNorm =
                        Foam::min(minBestSepNorm, bestSepNorm);
                    maxBestSepNorm =
                        Foam::max(maxBestSepNorm, bestSepNorm);

                    sumBestAlign += bestAlign;
                    minBestAlign =
                        Foam::min(minBestAlign, bestAlign);

                    sumBestOverlap += bestOverlap;
                    minBestOverlap =
                        Foam::min(minBestOverlap, bestOverlap);

                    // Diagnostic tiers only. These are dimensionless and
                    // are NOT production seam-acceptance tolerances.
                    if
                    (
                        bestSepNorm <= 1e-3 &&
                        bestAlign >= 0.999 &&
                        bestOverlap >= 0.5
                    )
                    {
                        ++nNearExact;
                    }

                    if
                    (
                        bestSepNorm <= 2e-2 &&
                        bestAlign >= 0.98 &&
                        bestOverlap >= 0.25
                    )
                    {
                        ++nStrong;
                    }
                }

                if( nFound == 0 )
                    continue;

                Info
                    << "[TRISURF_OPEN_PAIR]"
                    << " pair="
                    << censusPatchNames[p0] << "|"
                    << censusPatchNames[p1]
                    << " volumeFeaturePoints="
                    << usedVolumePair[pairI]
                    << " nativeSharedEdges="
                    << sharedCrossRegion[pairI]
                    << " openA=" << openEdgesByPatch[p0]
                    << " openB=" << openEdgesByPatch[p1]
                    << " queries=" << nQueries
                    << " found=" << nFound
                    << " bestSepNorm(min/avg/max)="
                    << minBestSepNorm << "/"
                    << sumBestSepNorm/scalar(nFound) << "/"
                    << maxBestSepNorm
                    << " bestAlign(min/avg)="
                    << minBestAlign << "/"
                    << sumBestAlign/scalar(nFound)
                    << " bestOverlap(min/avg)="
                    << minBestOverlap << "/"
                    << sumBestOverlap/scalar(nFound)
                    << " nearExact=" << nNearExact
                    << " strong=" << nStrong
                    << endl;
            }
        }

        Info
            << "[TRISURF_CENSUS_SUMMARY]"
            << " nonManifoldEdges=" << nNonManifoldEdges
            << endl;
    }

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

    // Diagnostic-only mapEdgeNodes outcome storage.
    // Each OMP iteration writes only its own index.
    //
    //  -1 unseen
    //   1 BL/no-BL protected
    //   2 no patches
    //   3 true feature edge not found
    //   4 ratio rejected
    //   5 true edge target outside mapping distance
    //   6 accepted
    labelList edgeDiagDecision(nodesToMap.size(), label(-1));
    boolList edgeDiagNeutral(nodesToMap.size(), false);

    pointField edgeDiagOverTarget(nodesToMap.size(), point::zero);
    scalarList edgeDiagOverTrueDist(nodesToMap.size(), scalar(-1));
    scalarList edgeDiagOverAllowedDist(nodesToMap.size(), scalar(-1));
    scalarList edgeDiagOverRatio(nodesToMap.size(), scalar(-1));


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
    // Diagnostic-only BL/neutral membership using this mapper's
    // CURRENT boundary-point addressing.
    if( !blNeutralPoints_.empty() )
    {
        forAll(nodesToMap, i)
        {
            if( blNeutralPoints_.found(nodesToMap[i]) )
                edgeDiagNeutral[i] = true;
        }
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
            edgeDiagDecision[i] = 1;
            if( tjIsWatched[i] ) tjDecision[i] = 3;
            continue;
        }
        const label bpI = nodesToMap[i];

        // liveP is the actual current mesh state.
        // queryP is only the geometric smoothing/search intent.
        const point& liveP = points[bPoints[bpI]];
        const point& queryP = desiredPositions[i];

        //- find patches at this edge point
        const DynList<label> patches = pPatches[bpI];
        if( patches.size() == 0 )
        {
            edgeDiagDecision[i] = 2;
            if( tjIsWatched[i] ) tjDecision[i] = 2;
            continue;  // guard: no patches = no valid snap
        }

        const scalar maxDist = mappingDistance[i];

        //- find approximate position of the vertex on the edge
        point mapPointApprox(queryP);
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
        distSqApprox = magSqr(mapPointApprox - queryP);

        //- find the nearest vertex on the triSurface feature edge
        point mapPoint(mapPointApprox);
        scalar distSq(distSqApprox);
        label nse;
        const bool edgeFoundOk =
            meshOctree_.findNearestEdgePoint
            (
                mapPoint,
                distSq,
                nse,
                queryP,
                patches
            );

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

          // Fail closed if the true feature-edge search failed.
          //
          // mapPointApprox is only an average of independent patch
          // projections. It is not guaranteed to lie on the true
          // patch-intersection curve.
          if( !edgeFoundOk )
          {
            edgeDiagDecision[i] = 3;
              projectedPoints[i] = point::zero;
              projectedDistSq[i] = scalar(-1);
              continue;
          }

          if( distSq > 1.2 * distSqApprox )
        {
            edgeDiagDecision[i] = 4;
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

        //- A successful feature-edge target must remain on the feature.
        //- Do not truncate an over-distance move by interpolating from p
        //- toward mapPoint: unless p already lies on the same straight
        //- feature segment, that interpolation produces an off-feature
        //- point and destroys the boundary-intersection constraint.
        //
        //- If the true feature target is outside the allowed mapping
        //- distance, fail closed and leave this vertex unchanged.
        const scalar actualMoveDistSq =
            magSqr(mapPoint - liveP);

        if( actualMoveDistSq > maxDist )
        {
            edgeDiagDecision[i] = 5;

            const scalar trueDist =
                Foam::sqrt
                (
                    Foam::max
                    (
                        actualMoveDistSq,
                        scalar(0)
                    )
                );

            const scalar allowedDist =
                Foam::sqrt(Foam::max(maxDist, scalar(0)));

            edgeDiagOverTarget[i] = mapPoint;
            edgeDiagOverTrueDist[i] = trueDist;
            edgeDiagOverAllowedDist[i] = allowedDist;
            edgeDiagOverRatio[i] =
                trueDist / Foam::max(allowedDist, VSMALL);

            projectedPoints[i] = point::zero;
            projectedDistSq[i] = scalar(-1);
            continue;
        }

                //- store only - no mesh mutation in parallel
        edgeDiagDecision[i] = 6;
        projectedPoints[i] = mapPoint;

        // Synchronization/commit distance is the real mesh movement,
        // not the query-position residual to the feature.
        projectedDistSq[i] = actualMoveDistSq;

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
    // Diagnostic-only population census. No geometry behavior change.
    {
        label nProtected = 0;
        label nNoPatches = 0;
        label nEdgeNotFound = 0;
        label nRatioRejected = 0;
        label nOverDistance = 0;
        label nAccepted = 0;
        label nUnseen = 0;

        label nNeutralInput = 0;
        label nNeutralOverDistance = 0;
        label nOrdinaryOverDistance = 0;
        label nNeutralAccepted = 0;

        scalar sumOverRatio = 0;
        scalar maxOverRatio = 0;

        label nPrinted = 0;

        forAll(nodesToMap, i)
        {
            if( edgeDiagNeutral[i] )
                ++nNeutralInput;

            switch( edgeDiagDecision[i] )
            {
                case 1:
                    ++nProtected;
                    break;

                case 2:
                    ++nNoPatches;
                    break;

                case 3:
                    ++nEdgeNotFound;
                    break;

                case 4:
                    ++nRatioRejected;
                    break;

                case 5:
                {
                    ++nOverDistance;

                    if( edgeDiagNeutral[i] )
                        ++nNeutralOverDistance;
                    else
                        ++nOrdinaryOverDistance;

                    sumOverRatio += edgeDiagOverRatio[i];

                    maxOverRatio =
                        Foam::max
                        (
                            maxOverRatio,
                            edgeDiagOverRatio[i]
                        );

                    if( nPrinted < 8 )
                    {
                        const label bpI = nodesToMap[i];
                        const DynList<label> patches = pPatches[bpI];

                        Info
                            << "[EDGEMAP_OVERDIST]"
                            << " bpI=" << bpI
                            << " meshPointI=" << bPoints[bpI]
                            << " neutral="
                            << (edgeDiagNeutral[i] ? 1 : 0)
                            << " patches=" << patches
                            << " from=" << oldPositions[i]
                            << " trueTarget="
                            << edgeDiagOverTarget[i]
                            << " trueDist="
                            << edgeDiagOverTrueDist[i]
                            << " allowedDist="
                            << edgeDiagOverAllowedDist[i]
                            << " ratio="
                            << edgeDiagOverRatio[i]
                            << endl;

                        ++nPrinted;
                    }

                    break;
                }

                case 6:
                    ++nAccepted;

                    if( edgeDiagNeutral[i] )
                        ++nNeutralAccepted;

                    break;

                default:
                    ++nUnseen;
                    break;
            }
        }

        // Diagnostic only: classify edge-map outcomes by the number of
        // surface patches incident to the boundary point.
        //
        // Buckets:
        //   0 = no patches
        //   1 = one patch
        //   2 = two patches
        //   3 = three or more patches
        //
        // This does not alter search, projection, acceptance, or motion.
        label edgeNotFoundByPatchCount[4] = {0, 0, 0, 0};
        label ratioRejectedByPatchCount[4] = {0, 0, 0, 0};
        label overDistanceByPatchCount[4] = {0, 0, 0, 0};
        label acceptedByPatchCount[4] = {0, 0, 0, 0};

        forAll(nodesToMap, i)
        {
            const label decision = edgeDiagDecision[i];

            if
            (
                decision != 3 &&
                decision != 4 &&
                decision != 5 &&
                decision != 6
            )
            {
                continue;
            }

            const label bpI = nodesToMap[i];
            const label nPatches = pPatches[bpI].size();

            label bucket = 3;
            if( nPatches <= 0 )
                bucket = 0;
            else if( nPatches == 1 )
                bucket = 1;
            else if( nPatches == 2 )
                bucket = 2;

            if( decision == 3 )
                ++edgeNotFoundByPatchCount[bucket];
            else if( decision == 4 )
                ++ratioRejectedByPatchCount[bucket];
            else if( decision == 5 )
                ++overDistanceByPatchCount[bucket];
            else if( decision == 6 )
                ++acceptedByPatchCount[bucket];
        }

        // ------------------------------------------------------------
        // Diagnostic-only SourceFeatureGraph recovery audit.
        //
        // Production has already completed the OMP mapping calculation.
        // Query only points classified as EDGE_NOT_FOUND (decision 3).
        //
        // Nothing below modifies:
        //   projectedPoints
        //   projectedDistSq
        //   mappingAccepted
        //   mesh point coordinates
        //
        // Virtual seam diagnostic policy:
        //   sourceDist / local edge length <= 0.25
        //   seam gap / min edge length      <= 0.02
        //   |tA.tB|                         >= 0.98
        //
        // SourceFeatureGraph also requires both open-edge targets to lie
        // within the existing production mappingDistance[i].
        // ------------------------------------------------------------
        {
            const sourceFeatureGraph& sfGraph = sourceFeatures();

            const triSurf& sfSurf = meshOctree_.surface();
            const wordList& sfPatchNames = sfSurf.patchNames();
            const edgeLongList& sfEdges = sfSurf.edges();
            const pointField& sfPoints = sfSurf.points();

            const label nSfPatches = sfPatchNames.size();
            const label nSfPairs = nSfPatches*nSfPatches;

            labelList sfFailed(nSfPairs, 0);
            labelList sfNative(nSfPairs, 0);
            labelList sfVirtual(nSfPairs, 0);
            labelList sfUnresolved(nSfPairs, 0);

            // Diagnostic decomposition of unresolved graph queries.
            labelList sfNoCandidate(nSfPairs, 0);
            labelList sfGeometryRejected(nSfPairs, 0);
            labelList sfTargetOverDistance(nSfPairs, 0);
            labelList sfSourceCapOnly(nSfPairs, 0);
            labelList sfNativeOverDistance(nSfPairs, 0);
            labelList sfOtherUnresolved(nSfPairs, 0);

            label sfTotalFailed = 0;
            label sfTotalNative = 0;
            label sfTotalVirtual = 0;
            label sfTotalUnresolved = 0;

            label sfTotalNoCandidate = 0;
            label sfTotalGeometryRejected = 0;
            label sfTotalTargetOverDistance = 0;
            label sfTotalSourceCapOnly = 0;
            label sfTotalNativeOverDistance = 0;
            label sfTotalOtherUnresolved = 0;

            forAll(nodesToMap, i)
            {
                if( edgeDiagDecision[i] != 3 )
                    continue;

                const label bpI = nodesToMap[i];

                if( pPatches.sizeOfRow(bpI) != 2 )
                    continue;

                label patchA = pPatches(bpI, 0);
                label patchB = pPatches(bpI, 1);

                if
                (
                    patchA < 0 ||
                    patchB < 0 ||
                    patchA >= nSfPatches ||
                    patchB >= nSfPatches ||
                    patchA == patchB
                )
                {
                    continue;
                }

                if( patchB < patchA )
                {
                    const label tmp = patchA;
                    patchA = patchB;
                    patchB = tmp;
                }

                const label pairI =
                    patchA*nSfPatches + patchB;

                ++sfFailed[pairI];
                ++sfTotalFailed;

                const point& queryPoint =
                    points[bPoints[bpI]];

                sourceFeatureGraph::featureHit hit;

                const bool recovered =
                    sfGraph.findNearestFeaturePoint
                    (
                        hit,
                        queryPoint,
                        patchA,
                        patchB,
                        mappingDistance[i],
                        scalar(0.02),
                        scalar(0.98)
                    );

                if( !recovered )
                {
                    ++sfUnresolved[pairI];
                    ++sfTotalUnresolved;

                    if( !hit.found )
                    {
                        ++sfNoCandidate[pairI];
                        ++sfTotalNoCandidate;
                    }
                    else if
                    (
                        hit.type ==
                        sourceFeatureGraph::NATIVE_FEATURE
                    )
                    {
                        // Native candidate exists, but the graph rejected
                        // it. In v0 the native admissibility gate is the
                        // production mapping-distance limit.
                        ++sfNativeOverDistance[pairI];
                        ++sfTotalNativeOverDistance;
                    }
                    else if
                    (
                        hit.type ==
                        sourceFeatureGraph::VIRTUAL_OPEN_SEAM &&
                        hit.sourceEdgeA >= 0 &&
                        hit.sourceEdgeB >= 0 &&
                        hit.sourceEdgeA < sfEdges.size() &&
                        hit.sourceEdgeB < sfEdges.size()
                    )
                    {
                        const edge& eA =
                            sfEdges[hit.sourceEdgeA];

                        const edge& eB =
                            sfEdges[hit.sourceEdgeB];

                        const scalar lenA =
                            Foam::sqrt
                            (
                                magSqr
                                (
                                    sfPoints[eA.end()]
                                  - sfPoints[eA.start()]
                                )
                            );

                        const scalar lenB =
                            Foam::sqrt
                            (
                                magSqr
                                (
                                    sfPoints[eB.end()]
                                  - sfPoints[eB.start()]
                                )
                            );

                        const scalar sourceDistANorm =
                            Foam::sqrt(hit.sourceDistASq)
                          / (lenA + VSMALL);

                        const scalar sourceDistBNorm =
                            Foam::sqrt(hit.sourceDistBSq)
                          / (lenB + VSMALL);

                        const bool geometryOk =
                        (
                            Foam::max
                            (
                                sourceDistANorm,
                                sourceDistBNorm
                            ) <= scalar(0.25)
                         && hit.gapNorm <= scalar(0.02)
                         && hit.alignment >= scalar(0.98)
                        );

                        const bool targetWithinDistance =
                        (
                            mappingDistance[i] <= VSMALL ||
                            hit.distSq <= mappingDistance[i]
                        );

                        const bool sourcesWithinDistance =
                        (
                            mappingDistance[i] <= VSMALL ||
                            (
                                hit.sourceDistASq <= mappingDistance[i] &&
                                hit.sourceDistBSq <= mappingDistance[i]
                            )
                        );

                        if( !geometryOk )
                        {
                            ++sfGeometryRejected[pairI];
                            ++sfTotalGeometryRejected;
                        }
                        else if( !targetWithinDistance )
                        {
                            ++sfTargetOverDistance[pairI];
                            ++sfTotalTargetOverDistance;
                        }
                        else if( !sourcesWithinDistance )
                        {
                            // Canonical virtual target is inside the
                            // production mapping radius, but one of the
                            // two support projections is outside it.
                            //
                            // This isolates an overly conservative
                            // SourceFeatureGraph v0 distance policy from
                            // a genuine production target-distance reject.
                            ++sfSourceCapOnly[pairI];
                            ++sfTotalSourceCapOnly;
                        }
                        else
                        {
                            ++sfOtherUnresolved[pairI];
                            ++sfTotalOtherUnresolved;
                        }
                    }
                    else
                    {
                        ++sfOtherUnresolved[pairI];
                        ++sfTotalOtherUnresolved;
                    }

                    continue;
                }

                if
                (
                    hit.type ==
                    sourceFeatureGraph::NATIVE_FEATURE
                )
                {
                    ++sfNative[pairI];
                    ++sfTotalNative;
                }
                else if
                (
                    hit.type ==
                    sourceFeatureGraph::VIRTUAL_OPEN_SEAM
                )
                {
                    ++sfVirtual[pairI];
                    ++sfTotalVirtual;
                }
                else
                {
                    ++sfUnresolved[pairI];
                    ++sfTotalUnresolved;
                }
            }

            for(label patchA=0; patchA<nSfPatches; ++patchA)
            {
                for
                (
                    label patchB=patchA+1;
                    patchB<nSfPatches;
                    ++patchB
                )
                {
                    const label pairI =
                        patchA*nSfPatches + patchB;

                    if( sfFailed[pairI] == 0 )
                        continue;

                    Info
                        << "[SOURCE_FEATURE_RECOVERY]"
                        << " pair="
                        << sfPatchNames[patchA] << "|"
                        << sfPatchNames[patchB]
                        << " edgeNotFound="
                        << sfFailed[pairI]
                        << " native="
                        << sfNative[pairI]
                        << " virtual="
                        << sfVirtual[pairI]
                        << " unresolved="
                        << sfUnresolved[pairI]
                        << " noCandidate="
                        << sfNoCandidate[pairI]
                        << " geometryRejected="
                        << sfGeometryRejected[pairI]
                        << " targetOverDistance="
                        << sfTargetOverDistance[pairI]
                        << " sourceCapOnly="
                        << sfSourceCapOnly[pairI]
                        << " nativeOverDistance="
                        << sfNativeOverDistance[pairI]
                        << " other="
                        << sfOtherUnresolved[pairI]
                        << endl;
                }
            }

            Info
                << "[SOURCE_FEATURE_RECOVERY_SUMMARY]"
                << " edgeNotFound=" << sfTotalFailed
                << " native=" << sfTotalNative
                << " virtual=" << sfTotalVirtual
                << " unresolved=" << sfTotalUnresolved
                << " noCandidate=" << sfTotalNoCandidate
                << " geometryRejected=" << sfTotalGeometryRejected
                << " targetOverDistance=" << sfTotalTargetOverDistance
                << " sourceCapOnly=" << sfTotalSourceCapOnly
                << " nativeOverDistance=" << sfTotalNativeOverDistance
                << " other=" << sfTotalOtherUnresolved
                << endl;
        }

        Info
            << "[EDGEMAP_PATCHCOUNT]"
            << " edgeNotFound(p0/p1/p2/p3p)="
            << edgeNotFoundByPatchCount[0] << "/"
            << edgeNotFoundByPatchCount[1] << "/"
            << edgeNotFoundByPatchCount[2] << "/"
            << edgeNotFoundByPatchCount[3]
            << " ratioRejected(p0/p1/p2/p3p)="
            << ratioRejectedByPatchCount[0] << "/"
            << ratioRejectedByPatchCount[1] << "/"
            << ratioRejectedByPatchCount[2] << "/"
            << ratioRejectedByPatchCount[3]
            << " overDistance(p0/p1/p2/p3p)="
            << overDistanceByPatchCount[0] << "/"
            << overDistanceByPatchCount[1] << "/"
            << overDistanceByPatchCount[2] << "/"
            << overDistanceByPatchCount[3]
            << " accepted(p0/p1/p2/p3p)="
            << acceptedByPatchCount[0] << "/"
            << acceptedByPatchCount[1] << "/"
            << acceptedByPatchCount[2] << "/"
            << acceptedByPatchCount[3]
            << endl;

        // Diagnostic only: classify mapping outcomes by the number of
        // actual cross-patch feature edges attached to the boundary point.
        //
        // Buckets:
        //   0 = no local feature edge
        //   1 = dangling / endpoint-like feature support
        //   2 = ordinary continuous feature curve
        //   3 = three or more feature edges
        label edgeNotFoundByFeatureCount[4] = {0, 0, 0, 0};
        label ratioRejectedByFeatureCount[4] = {0, 0, 0, 0};
        label overDistanceByFeatureCount[4] = {0, 0, 0, 0};
        label acceptedByFeatureCount[4] = {0, 0, 0, 0};

        forAll(nodesToMap, i)
        {
            const label decision = edgeDiagDecision[i];

            if
            (
                decision != 3 &&
                decision != 4 &&
                decision != 5 &&
                decision != 6
            )
            {
                continue;
            }

            const label bpI = nodesToMap[i];
            const label nFeatureEdges =
                mPart.numberOfFeatureEdgesAtPoint(bpI);

            label bucket = 3;
            if( nFeatureEdges <= 0 )
                bucket = 0;
            else if( nFeatureEdges == 1 )
                bucket = 1;
            else if( nFeatureEdges == 2 )
                bucket = 2;

            if( decision == 3 )
                ++edgeNotFoundByFeatureCount[bucket];
            else if( decision == 4 )
                ++ratioRejectedByFeatureCount[bucket];
            else if( decision == 5 )
                ++overDistanceByFeatureCount[bucket];
            else if( decision == 6 )
                ++acceptedByFeatureCount[bucket];
        }

        // Diagnostic only: cross-tab exact mapper outcomes by the
        // physical two-patch pair carried by each boundary point.
        //
        // Outcome slots:
        //   0 = EDGE_NOT_FOUND
        //   1 = RATIO_REJECTED
        //   2 = OVER_DISTANCE
        //   3 = ACCEPTED
        const wordList& edgeMapPatchNames =
            meshOctree_.surface().patchNames();

        const label nEdgeMapPatches = edgeMapPatchNames.size();

        labelList edgeMapPairOutcome
        (
            4*nEdgeMapPatches*nEdgeMapPatches,
            0
        );

        label nInvalidEdgeMapPairs = 0;

        forAll(nodesToMap, i)
        {
            const label decision = edgeDiagDecision[i];

            label outcome = -1;
            if( decision == 3 )
                outcome = 0;
            else if( decision == 4 )
                outcome = 1;
            else if( decision == 5 )
                outcome = 2;
            else if( decision == 6 )
                outcome = 3;
            else
                continue;

            const label bpI = nodesToMap[i];

            if( pPatches.sizeOfRow(bpI) != 2 )
            {
                ++nInvalidEdgeMapPairs;
                continue;
            }

            label p0 = pPatches(bpI, 0);
            label p1 = pPatches(bpI, 1);

            if
            (
                p0 < 0 || p1 < 0 ||
                p0 >= nEdgeMapPatches ||
                p1 >= nEdgeMapPatches
            )
            {
                ++nInvalidEdgeMapPairs;
                continue;
            }

            if( p1 < p0 )
            {
                const label tmp = p0;
                p0 = p1;
                p1 = tmp;
            }

            const label pairI =
                p0*nEdgeMapPatches + p1;

            ++edgeMapPairOutcome[4*pairI + outcome];
        }

        for(label p0=0; p0<nEdgeMapPatches; ++p0)
        {
            for(label p1=p0; p1<nEdgeMapPatches; ++p1)
            {
                const label pairI =
                    p0*nEdgeMapPatches + p1;

                const label nNotFound =
                    edgeMapPairOutcome[4*pairI + 0];
                const label nRatio =
                    edgeMapPairOutcome[4*pairI + 1];
                const label nOver =
                    edgeMapPairOutcome[4*pairI + 2];
                const label nAccept =
                    edgeMapPairOutcome[4*pairI + 3];

                const label total =
                    nNotFound + nRatio + nOver + nAccept;

                if( total == 0 )
                    continue;

                Info
                    << "[EDGEMAP_PAIR]"
                    << " pair="
                    << edgeMapPatchNames[p0] << "|"
                    << edgeMapPatchNames[p1]
                    << " ids=" << p0 << "|" << p1
                    << " input=" << total
                    << " edgeNotFound=" << nNotFound
                    << " ratioRejected=" << nRatio
                    << " overDistance=" << nOver
                    << " accepted=" << nAccept
                    << endl;
            }
        }

        if( nInvalidEdgeMapPairs )
        {
            Info
                << "[EDGEMAP_PAIR_INVALID]"
                << " count=" << nInvalidEdgeMapPairs
                << endl;
        }

        Info
            << "[EDGEMAP_FEATURECOUNT]"
            << " edgeNotFound(e0/e1/e2/e3p)="
            << edgeNotFoundByFeatureCount[0] << "/"
            << edgeNotFoundByFeatureCount[1] << "/"
            << edgeNotFoundByFeatureCount[2] << "/"
            << edgeNotFoundByFeatureCount[3]
            << " ratioRejected(e0/e1/e2/e3p)="
            << ratioRejectedByFeatureCount[0] << "/"
            << ratioRejectedByFeatureCount[1] << "/"
            << ratioRejectedByFeatureCount[2] << "/"
            << ratioRejectedByFeatureCount[3]
            << " overDistance(e0/e1/e2/e3p)="
            << overDistanceByFeatureCount[0] << "/"
            << overDistanceByFeatureCount[1] << "/"
            << overDistanceByFeatureCount[2] << "/"
            << overDistanceByFeatureCount[3]
            << " accepted(e0/e1/e2/e3p)="
            << acceptedByFeatureCount[0] << "/"
            << acceptedByFeatureCount[1] << "/"
            << acceptedByFeatureCount[2] << "/"
            << acceptedByFeatureCount[3]
            << endl;

        // Diagnostic only: point-conditioned correspondence between
        // EDGE_NOT_FOUND volume feature points and open source boundaries.
        //
        // For each failed two-patch point p on A/B:
        //   1. Find nearest OPEN source edge owned by A.
        //   2. Find nearest OPEN source edge owned by B.
        //   3. Compare the two projected source positions and tangents.
        //
        // This asks directly whether the feature which native A/B lookup
        // failed to find is represented instead by two separate open source
        // boundaries at that SAME physical location.
        //
        // No mesh mutation and no additional production octree queries.
        {
            const triSurf& corrSurf = meshOctree_.surface();
            const VRWGraph& corrEdgeFacets = corrSurf.edgeFacets();
            const LongList<labelledTri>& corrFacets = corrSurf.facets();
            const LongList<edge>& corrEdges = corrSurf.edges();
            const pointField& corrPoints = corrSurf.points();
            const wordList& corrPatchNames = corrSurf.patchNames();

            const label nCorrPatches = corrPatchNames.size();
            const label nPairSlots =
                nCorrPatches*nCorrPatches;

            // Collect source OPEN edges once.
            DynList<label> corrOpenEdgeIds;
            DynList<label> corrOpenEdgePatch;
            labelList corrOpenCount(nCorrPatches, 0);

            forAll(corrEdgeFacets, eI)
            {
                if( corrEdgeFacets.sizeOfRow(eI) != 1 )
                    continue;

                const label fI = corrEdgeFacets(eI, 0);

                if( fI < 0 || fI >= corrFacets.size() )
                    continue;

                const label patchI =
                    corrFacets[fI].region();

                if
                (
                    patchI < 0 ||
                    patchI >= nCorrPatches
                )
                {
                    continue;
                }

                corrOpenEdgeIds.append(eI);
                corrOpenEdgePatch.append(patchI);
                ++corrOpenCount[patchI];
            }

            labelList nFailed(nPairSlots, 0);
            labelList nBothOpen(nPairSlots, 0);
            labelList nNearBoth(nPairSlots, 0);
            labelList nGapSmall(nPairSlots, 0);
            labelList nAligned(nPairSlots, 0);
            labelList nStrong(nPairSlots, 0);
            labelList nNearExact(nPairSlots, 0);

            scalarList sumDistANorm(nPairSlots, scalar(0));
            scalarList sumDistBNorm(nPairSlots, scalar(0));
            scalarList sumGapNorm(nPairSlots, scalar(0));
            scalarList sumAlign(nPairSlots, scalar(0));

            scalarList maxDistANorm(nPairSlots, scalar(0));
            scalarList maxDistBNorm(nPairSlots, scalar(0));
            scalarList maxGapNorm(nPairSlots, scalar(0));
            scalarList minAlign(nPairSlots, GREAT);

            forAll(nodesToMap, i)
            {
                // Production decision 3 == EDGE_NOT_FOUND.
                if( edgeDiagDecision[i] != 3 )
                    continue;

                const label bpI = nodesToMap[i];

                if( pPatches.sizeOfRow(bpI) != 2 )
                    continue;

                label pA = pPatches(bpI, 0);
                label pB = pPatches(bpI, 1);

                if
                (
                    pA < 0 || pB < 0 ||
                    pA >= nCorrPatches ||
                    pB >= nCorrPatches ||
                    pA == pB
                )
                {
                    continue;
                }

                if( pB < pA )
                {
                    const label tmp = pA;
                    pA = pB;
                    pB = tmp;
                }

                const label pairI =
                    pA*nCorrPatches + pB;

                ++nFailed[pairI];

                const point& pQuery =
                    points[bPoints[bpI]];

                bool foundA = false;
                bool foundB = false;

                scalar bestASq = GREAT;
                scalar bestBSq = GREAT;

                point qA(pQuery);
                point qB(pQuery);

                scalar lenA = 0;
                scalar lenB = 0;

                vector dirA(vector::zero);
                vector dirB(vector::zero);

                forAll(corrOpenEdgeIds, oi)
                {
                    const label owner =
                        corrOpenEdgePatch[oi];

                    if( owner != pA && owner != pB )
                        continue;

                    const label eI =
                        corrOpenEdgeIds[oi];

                    const edge& e = corrEdges[eI];

                    const point& e0 =
                        corrPoints[e.start()];
                    const point& e1 =
                        corrPoints[e.end()];

                    const vector d = e1-e0;
                    const scalar len =
                        Foam::sqrt(magSqr(d));

                    if( len < ROOTVSMALL )
                        continue;

                    const point q =
                        help::nearestPointOnTheEdgeExact
                        (
                            e0,
                            e1,
                            pQuery
                        );

                    const scalar dSq =
                        magSqr(q-pQuery);

                    if( owner == pA && dSq < bestASq )
                    {
                        foundA = true;
                        bestASq = dSq;
                        qA = q;
                        lenA = len;
                        dirA = d;
                    }

                    if( owner == pB && dSq < bestBSq )
                    {
                        foundB = true;
                        bestBSq = dSq;
                        qB = q;
                        lenB = len;
                        dirB = d;
                    }
                }

                if( !foundA || !foundB )
                    continue;

                ++nBothOpen[pairI];

                const scalar dA =
                    Foam::sqrt(bestASq);
                const scalar dB =
                    Foam::sqrt(bestBSq);

                const scalar dANorm =
                    dA/(lenA + VSMALL);
                const scalar dBNorm =
                    dB/(lenB + VSMALL);

                const scalar minLen =
                    Foam::min(lenA, lenB);

                const scalar gap =
                    Foam::sqrt(magSqr(qA-qB));

                const scalar gapNorm =
                    gap/(minLen + VSMALL);

                const scalar align =
                    Foam::min
                    (
                        scalar(1),
                        Foam::mag(dirA & dirB) /
                        (lenA*lenB + VSMALL)
                    );

                sumDistANorm[pairI] += dANorm;
                sumDistBNorm[pairI] += dBNorm;
                sumGapNorm[pairI] += gapNorm;
                sumAlign[pairI] += align;

                maxDistANorm[pairI] =
                    Foam::max
                    (
                        maxDistANorm[pairI],
                        dANorm
                    );

                maxDistBNorm[pairI] =
                    Foam::max
                    (
                        maxDistBNorm[pairI],
                        dBNorm
                    );

                maxGapNorm[pairI] =
                    Foam::max
                    (
                        maxGapNorm[pairI],
                        gapNorm
                    );

                minAlign[pairI] =
                    Foam::min
                    (
                        minAlign[pairI],
                        align
                    );

                // Diagnostic categories only. These are NOT proposed
                // production tolerances.
                const bool nearBoth =
                    Foam::max(dANorm, dBNorm) <= 0.25;

                const bool gapSmall =
                    gapNorm <= 0.02;

                const bool aligned =
                    align >= 0.98;

                if( nearBoth )
                    ++nNearBoth[pairI];

                if( gapSmall )
                    ++nGapSmall[pairI];

                if( aligned )
                    ++nAligned[pairI];

                if( nearBoth && gapSmall && aligned )
                    ++nStrong[pairI];

                if
                (
                    Foam::max(dANorm, dBNorm) <= 0.05 &&
                    gapNorm <= 1e-3 &&
                    align >= 0.999
                )
                {
                    ++nNearExact[pairI];
                }
            }

            for(label pA=0; pA<nCorrPatches; ++pA)
            {
                for(label pB=pA+1; pB<nCorrPatches; ++pB)
                {
                    const label pairI =
                        pA*nCorrPatches + pB;

                    if( nFailed[pairI] == 0 )
                        continue;

                    Info
                        << "[EDGEMAP_OPEN_CORRESPONDENCE]"
                        << " pair="
                        << corrPatchNames[pA] << "|"
                        << corrPatchNames[pB]
                        << " edgeNotFound="
                        << nFailed[pairI]
                        << " sourceOpenA="
                        << corrOpenCount[pA]
                        << " sourceOpenB="
                        << corrOpenCount[pB]
                        << " bothOpen="
                        << nBothOpen[pairI]
                        << " nearBoth="
                        << nNearBoth[pairI]
                        << " gapSmall="
                        << nGapSmall[pairI]
                        << " aligned="
                        << nAligned[pairI]
                        << " strong="
                        << nStrong[pairI]
                        << " nearExact="
                        << nNearExact[pairI];

                    if( nBothOpen[pairI] > 0 )
                    {
                        const scalar denom =
                            scalar(nBothOpen[pairI]);

                        Info
                            << " distANorm(avg/max)="
                            << sumDistANorm[pairI]/denom
                            << "/"
                            << maxDistANorm[pairI]
                            << " distBNorm(avg/max)="
                            << sumDistBNorm[pairI]/denom
                            << "/"
                            << maxDistBNorm[pairI]
                            << " gapNorm(avg/max)="
                            << sumGapNorm[pairI]/denom
                            << "/"
                            << maxGapNorm[pairI]
                            << " align(min/avg)="
                            << minAlign[pairI]
                            << "/"
                            << sumAlign[pairI]/denom;
                    }

                    Info << endl;
                }
            }
        }

        Info
            << "[EDGEMAPSUMMARY]"
            << " input=" << nodesToMap.size()
            << " protected=" << nProtected
            << " neutralInput=" << nNeutralInput
            << " edgeNotFound=" << nEdgeNotFound
            << " ratioRejected=" << nRatioRejected
            << " overDistanceRejected=" << nOverDistance
            << " neutralOverDistance=" << nNeutralOverDistance
            << " ordinaryOverDistance=" << nOrdinaryOverDistance
            << " accepted=" << nAccepted
            << " neutralAccepted=" << nNeutralAccepted
            << " unseen=" << nUnseen
            << " avgOverDistanceRatio="
            <<
            (
                nOverDistance
              ? sumOverRatio / scalar(nOverDistance)
              : scalar(0)
            )
            << " maxOverDistanceRatio=" << maxOverRatio
            << endl;
    }



    // ------------------------------------------------------------
    // DIAGNOSTIC ONLY:
    // Probe the existing production edge-move path against raw signed
    // incident-cell volumes BEFORE committing each existing move.
    //
    // No movement is rejected or modified here.
    // ------------------------------------------------------------
    {
        const polyMeshGen& probeMesh = surfaceEngine_.mesh();

        const pointFieldPMG& probePoints = probeMesh.points();
        const faceListPMG& probeFaces = probeMesh.faces();
        const cellListPMG& probeCells = probeMesh.cells();
        const labelList& probeOwner = probeMesh.owner();

        const VRWGraph& probePointCells =
            probeMesh.addressingData().pointCells();

        const labelList& probeBoundaryPoints =
            surfaceEngine_.boundaryPoints();

        auto probeFaceGeometry =
        [&]
        (
            const face& f,
            const label movedGlobalPtI,
            const point& candidate,
            const bool substituteCandidate,
            point& fCtr,
            vector& fArea
        )
        {
            auto pt =
            [&]
            (
                const label ptI
            ) -> point
            {
                if
                (
                    substituteCandidate &&
                    ptI == movedGlobalPtI
                )
                {
                    return candidate;
                }

                return point(probePoints[ptI]);
            };

            const label nPoints = f.size();

            if( nPoints == 3 )
            {
                const point p0 = pt(f[0]);
                const point p1 = pt(f[1]);
                const point p2 = pt(f[2]);

                fCtr = (1.0/3.0)*(p0 + p1 + p2);
                fArea = 0.5*((p1 - p0)^(p2 - p0));
                return;
            }

            vector sumN = vector::zero;
            scalar sumA = 0.0;
            vector sumAc = vector::zero;

            point fCentre = pt(f[0]);

            for(label pi=1; pi<nPoints; ++pi)
            {
                fCentre += pt(f[pi]);
            }

            fCentre /= nPoints;

            for(label pi=0; pi<nPoints; ++pi)
            {
                const point curPoint = pt(f[pi]);
                const point nextPoint = pt(f.nextLabel(pi));

                const vector c =
                    curPoint + nextPoint + fCentre;

                const vector n =
                    (nextPoint - curPoint)^
                    (fCentre - curPoint);

                const scalar a = mag(n);

                sumN += n;
                sumA += a;
                sumAc += a*c;
            }

            fCtr =
                (1.0/3.0)*sumAc/(sumA + VSMALL);

            fArea = 0.5*sumN;
        };


        auto probeRawCellVolume =
        [&]
        (
            const label cellI,
            const label movedGlobalPtI,
            const point& candidate,
            const bool substituteCandidate
        ) -> scalar
        {
            const cell& c = probeCells[cellI];

            List<point> localFCentres(c.size());
            List<vector> localFAreas(c.size());

            vector cEst = vector::zero;

            forAll(c, cfI)
            {
                const label faceI = c[cfI];

                point fc;
                vector fa;

                probeFaceGeometry
                (
                    probeFaces[faceI],
                    movedGlobalPtI,
                    candidate,
                    substituteCandidate,
                    fc,
                    fa
                );

                localFCentres[cfI] = fc;
                localFAreas[cfI] = fa;

                cEst += fc;
            }

            cEst /= scalar(c.size());

            scalar cellVol = 0.0;

            forAll(c, cfI)
            {
                const label faceI = c[cfI];

                scalar pyr3Vol =
                    localFAreas[cfI] &
                    (localFCentres[cfI] - cEst);

                if( probeOwner[faceI] != cellI )
                {
                    pyr3Vol *= -1.0;
                }

                cellVol += pyr3Vol;
            }

            return cellVol/3.0;
        };


        label nAppliedMoves = 0;
        label nMovesTouchingExistingBad = 0;
        label nMovesCreatingNewNeg = 0;
        label nMovesRepairingNeg = 0;
        label nRejectedCreatingNewNeg = 0;

        label nNewNegCellEvents = 0;
        label nRepairedCellEvents = 0;

        scalar minPositiveCandidateRatio = GREAT;

        forAll(nodesToMap, i)
        {
            const label bpI = nodesToMap[i];

            if( projectedDistSq[i] < scalar(0) )
            {
                continue;
            }

            ++nAppliedMoves;

            const label globalPtI =
                probeBoundaryPoints[bpI];

            const point candidate =
                projectedPoints[i];

            bool touchesExistingBad = false;
            bool createsNewNeg = false;
            bool repairsNeg = false;

            const label nIncidentCells =
                probePointCells.sizeOfRow(globalPtI);

            forAllRow(probePointCells, globalPtI, pcI)
            {
                const label cellI =
                    probePointCells(globalPtI, pcI);

                const scalar oldVol =
                    probeRawCellVolume
                    (
                        cellI,
                        globalPtI,
                        candidate,
                        false
                    );

                const scalar newVol =
                    probeRawCellVolume
                    (
                        cellI,
                        globalPtI,
                        candidate,
                        true
                    );

                if( oldVol < VSMALL )
                {
                    touchesExistingBad = true;
                }

                if
                (
                    oldVol >= VSMALL &&
                    newVol < VSMALL
                )
                {
                    createsNewNeg = true;
                    ++nNewNegCellEvents;
                }

                if
                (
                    oldVol < VSMALL &&
                    newVol >= VSMALL
                )
                {
                    repairsNeg = true;
                    ++nRepairedCellEvents;
                }

                if
                (
                    oldVol >= VSMALL &&
                    newVol >= VSMALL
                )
                {
                    const scalar ratio =
                        newVol/(oldVol + VSMALL);

                    minPositiveCandidateRatio =
                        Foam::min
                        (
                            minPositiveCandidateRatio,
                            ratio
                        );
                }
            }

            if( touchesExistingBad )
            {
                ++nMovesTouchingExistingBad;
            }

            if( createsNewNeg )
            {
                ++nMovesCreatingNewNeg;
            }

            if( repairsNeg )
            {
                ++nMovesRepairingNeg;
            }

            // Production edge transaction v1:
            // never allow a currently-positive incident cell to become
            // zero/negative. Full target or reject; no bisection yet.
            if( createsNewNeg )
            {
                ++nRejectedCreatingNewNeg;

                // Mark this production proposal as not committed so the
                // existing parallel synchronization path also skips it.
                projectedDistSq[i] = scalar(-1);
                mappingAccepted[i] = false;

                continue;
            }

            sMod.moveBoundaryVertexNoUpdate
            (
                bpI,
                projectedPoints[i]
            );

            mappingAccepted[i] = true;
        }

        Info
            << "[PRODUCTION_EDGE_VOLUME_PROBE]"
            << " appliedMoves="
            << nAppliedMoves
            << " movesTouchingExistingBad="
            << nMovesTouchingExistingBad
            << " movesCreatingNewNeg="
            << nMovesCreatingNewNeg
            << " newNegCellEvents="
            << nNewNegCellEvents
            << " rejectedCreatingNewNeg="
            << nRejectedCreatingNewNeg
            << " movesRepairingNeg="
            << nMovesRepairingNeg
            << " repairedCellEvents="
            << nRepairedCellEvents
            << " minPositiveCandidateRatio="
            <<
            (
                minPositiveCandidateRatio < GREAT
              ? minPositiveCandidateRatio
              : scalar(-1)
            )
            << endl;
    }

    // ------------------------------------------------------------
    // SourceFeatureGraph v2 behavioral experiment:
    //
    // Recover production EDGE_NOT_FOUND points only when:
    //   1. SourceFeatureGraph finds an admissible VIRTUAL_OPEN_SEAM,
    //   2. the existing proposedMoveIsValid() guard passes, AND
    //   3. every cell incident to the moved mesh point retains a
    //      strictly positive raw signed volume using the SAME volume
    //      convention as polyMeshGenChecks::checkCellVolumes().
    //
    // Existing successful/native edge mappings have already been moved
    // by the serial pass above. Virtual recoveries are then tested and
    // committed ONE AT A TIME. The raw volume evaluator uses live point
    // coordinates directly, not cached face/cell geometry, so each
    // subsequent candidate sees all earlier accepted moves.
    //
    // Full-target-or-reject only. No chord bisection.
    // ------------------------------------------------------------
    {
        const sourceFeatureGraph& sfg = sourceFeatures();

        const triSurf& sfgSurf = meshOctree_.surface();
        const wordList& sfgPatchNames = sfgSurf.patchNames();
        const label nSfgPatches = sfgPatchNames.size();

        const polyMeshGen& sfgMesh = surfaceEngine_.mesh();

        const pointFieldPMG& sfgPoints = sfgMesh.points();
        const faceListPMG& sfgFaces = sfgMesh.faces();
        const cellListPMG& sfgCells = sfgMesh.cells();
        const labelList& sfgOwner = sfgMesh.owner();

        // Topological addressing only. Safe to reuse while coordinates move.
        const VRWGraph& sfgPointCells =
            sfgMesh.addressingData().pointCells();

        const labelList& sfgBoundaryPoints =
            surfaceEngine_.boundaryPoints();

        // Exact local face geometry copied from
        // polyMeshGenAddressing::makeFaceCentresAndAreas(), except that
        // one selected mesh point may be replaced by candidate.
        auto sfgFaceGeometry =
        [&]
        (
            const face& f,
            const label movedGlobalPtI,
            const point& candidate,
            const bool substituteCandidate,
            point& fCtr,
            vector& fArea
        )
        {
            auto pt =
            [&]
            (
                const label ptI
            ) -> point
            {
                if
                (
                    substituteCandidate &&
                    ptI == movedGlobalPtI
                )
                {
                    return candidate;
                }

                return point(sfgPoints[ptI]);
            };

            const label nPoints = f.size();

            if( nPoints == 3 )
            {
                const point p0 = pt(f[0]);
                const point p1 = pt(f[1]);
                const point p2 = pt(f[2]);

                fCtr = (1.0/3.0)*(p0 + p1 + p2);
                fArea = 0.5*((p1 - p0)^(p2 - p0));

                return;
            }

            vector sumN = vector::zero;
            scalar sumA = 0.0;
            vector sumAc = vector::zero;

            point fCentre = pt(f[0]);

            for(label pi=1; pi<nPoints; ++pi)
            {
                fCentre += pt(f[pi]);
            }

            fCentre /= nPoints;

            for(label pi=0; pi<nPoints; ++pi)
            {
                const point curPoint = pt(f[pi]);
                const point nextPoint = pt(f.nextLabel(pi));

                const vector c =
                    curPoint + nextPoint + fCentre;

                const vector n =
                    (nextPoint - curPoint)^
                    (fCentre - curPoint);

                const scalar a = mag(n);

                sumN += n;
                sumA += a;
                sumAc += a*c;
            }

            fCtr =
                (1.0/3.0)*sumAc/(sumA + VSMALL);

            fArea = 0.5*sumN;
        };

        // Raw signed volume using the exact orientation convention from
        // polyMeshGenChecks::checkCellVolumes().
        //
        // IMPORTANT: all faces are recalculated from current point
        // coordinates. This deliberately avoids cached face geometry,
        // which is stale until the normal updateGeometry() at the end
        // of mapEdgeNodes().
        auto sfgRawCellVolume =
        [&]
        (
            const label cellI,
            const label movedGlobalPtI,
            const point& candidate,
            const bool substituteCandidate
        ) -> scalar
        {
            const cell& c = sfgCells[cellI];

            List<point> localFCentres(c.size());
            List<vector> localFAreas(c.size());

            vector cEst = vector::zero;

            forAll(c, cfI)
            {
                const label faceI = c[cfI];

                point fc;
                vector fa;

                sfgFaceGeometry
                (
                    sfgFaces[faceI],
                    movedGlobalPtI,
                    candidate,
                    substituteCandidate,
                    fc,
                    fa
                );

                localFCentres[cfI] = fc;
                localFAreas[cfI] = fa;

                cEst += fc;
            }

            cEst /= scalar(c.size());

            scalar cellVol = 0.0;

            forAll(c, cfI)
            {
                const label faceI = c[cfI];

                scalar pyr3Vol =
                    localFAreas[cfI] &
                    (localFCentres[cfI] - cEst);

                if( sfgOwner[faceI] != cellI )
                {
                    pyr3Vol *= -1.0;
                }

                cellVol += pyr3Vol;
            }

            return cellVol/3.0;
        };


        const label nSfgPairs =
            nSfgPatches*nSfgPatches;

        labelList acceptedByPair(nSfgPairs, 0);
        labelList volumeRejectedByPair(nSfgPairs, 0);
        labelList faceRejectedByPair(nSfgPairs, 0);

        label nEdgeNotFoundSeen = 0;
        label nGraphUnresolved = 0;
        label nNativeDeferred = 0;
        label nVirtualCandidate = 0;

        label nAccepted = 0;

        label nRejectedFaceOnly = 0;
        label nRejectedVolumeOnly = 0;
        label nRejectedFaceAndVolume = 0;
        label nRejectedPreExistingBad = 0;
        label nRejectedNoCells = 0;
        label nRejectedActualMoveDistance = 0;

        scalar minAcceptedVolRatio = GREAT;
        scalar minCandidateVolRatio = GREAT;
        scalar maxAcceptedDist = scalar(0);

        forAll(nodesToMap, i)
        {
            // Production EDGE_NOT_FOUND only.
            if( edgeDiagDecision[i] != 3 )
            {
                continue;
            }

            ++nEdgeNotFoundSeen;

            const label bpI = nodesToMap[i];

            if( pPatches.sizeOfRow(bpI) != 2 )
            {
                ++nGraphUnresolved;
                continue;
            }

            label patchA = pPatches(bpI, 0);
            label patchB = pPatches(bpI, 1);

            if
            (
                patchA < 0 ||
                patchB < 0 ||
                patchA >= nSfgPatches ||
                patchB >= nSfgPatches ||
                patchA == patchB
            )
            {
                ++nGraphUnresolved;
                continue;
            }

            if( patchB < patchA )
            {
                const label tmp = patchA;
                patchA = patchB;
                patchB = tmp;
            }

            const label pairI =
                patchA*nSfgPatches + patchB;

            const label globalPtI =
                sfgBoundaryPoints[bpI];

            // COPY live position. Do not retain a reference across commit.
            const point oldPosition =
                sfgPoints[globalPtI];

            // Geometry search follows the smoothing intent. The mesh itself
            // remains at oldPosition until the final feature target passes
            // the transaction.
            const point queryPosition =
                desiredPositions[i];

            sourceFeatureGraph::featureHit hit;

            const bool recovered =
                sfg.findNearestFeaturePoint
                (
                    hit,
                    queryPosition,
                    patchA,
                    patchB,
                    mappingDistance[i],
                    scalar(0.02),
                    scalar(0.98)
                );

            if( !recovered )
            {
                ++nGraphUnresolved;
                continue;
            }

            if
            (
                hit.type ==
                sourceFeatureGraph::NATIVE_FEATURE
            )
            {
                // Separate experiment later.
                ++nNativeDeferred;
                continue;
            }

            if
            (
                hit.type !=
                sourceFeatureGraph::VIRTUAL_OPEN_SEAM
            )
            {
                ++nGraphUnresolved;
                continue;
            }

            // SourceFeatureGraph search distance is measured from
            // queryPosition. For smoothing-intent callers queryPosition
            // differs from the actual live mesh point, so independently
            // enforce the real live-point displacement authority here.
            const scalar actualVirtualMoveDistSq =
                magSqr(hit.position - oldPosition);

            if
            (
                actualVirtualMoveDistSq >
                mappingDistance[i]
            )
            {
                ++nRejectedActualMoveDistance;
                continue;
            }

            ++nVirtualCandidate;

            // Existing mapper quality proxy. This operates on the current
            // live point state and substitutes this candidate locally.
            const bool faceGateOK =
                proposedMoveIsValid
                (
                    bpI,
                    hit.position,
                    oldPosition,
                    mappingDistance[i]
                );

            // Exact raw signed-volume gate.
            bool volumeGateOK = true;
            bool preExistingBad = false;

            scalar candidateMinRatio = GREAT;

            const label nIncidentCells =
                sfgPointCells.sizeOfRow(globalPtI);

            if( nIncidentCells == 0 )
            {
                volumeGateOK = false;
                ++nRejectedNoCells;
            }
            else
            {
                forAllRow(sfgPointCells, globalPtI, pcI)
                {
                    const label cellI =
                        sfgPointCells(globalPtI, pcI);

                    const scalar oldVol =
                        sfgRawCellVolume
                        (
                            cellI,
                            globalPtI,
                            hit.position,
                            false
                        );

                    const scalar newVol =
                        sfgRawCellVolume
                        (
                            cellI,
                            globalPtI,
                            hit.position,
                            true
                        );

                    // For this first transaction experiment, do not touch
                    // a cell which is already zero/negative. This makes
                    // attribution unambiguous and fail-closed.
                    if( oldVol < VSMALL )
                    {
                        preExistingBad = true;
                        volumeGateOK = false;
                        break;
                    }

                    const scalar ratio =
                        newVol/(oldVol + VSMALL);

                    candidateMinRatio =
                        Foam::min
                        (
                            candidateMinRatio,
                            ratio
                        );

                    // Same sign criterion as checkCellVolumes().
                    if( newVol < VSMALL )
                    {
                        volumeGateOK = false;
                        break;
                    }
                }
            }

            if( candidateMinRatio < GREAT )
            {
                minCandidateVolRatio =
                    Foam::min
                    (
                        minCandidateVolRatio,
                        candidateMinRatio
                    );
            }

            if( preExistingBad )
            {
                ++nRejectedPreExistingBad;
            }

            if( !faceGateOK )
            {
                ++faceRejectedByPair[pairI];
            }

            if( !volumeGateOK )
            {
                ++volumeRejectedByPair[pairI];
            }

            if( !faceGateOK || !volumeGateOK )
            {
                if( !faceGateOK && !volumeGateOK )
                {
                    ++nRejectedFaceAndVolume;
                }
                else if( !faceGateOK )
                {
                    ++nRejectedFaceOnly;
                }
                else
                {
                    ++nRejectedVolumeOnly;
                }

                continue;
            }

            // Transaction commit.
            //
            // Commit immediately so every later candidate is evaluated
            // against the state containing all earlier accepted moves.
            sMod.moveBoundaryVertexNoUpdate
            (
                bpI,
                hit.position
            );

            projectedPoints[i] = hit.position;

            // hit.distSq is measured from the geometric query position.
            // Synchronization uses the already-validated actual live-point
            // displacement.
            projectedDistSq[i] =
                actualVirtualMoveDistSq;

            mappingAccepted[i] = true;

            ++acceptedByPair[pairI];
            ++nAccepted;

            if( candidateMinRatio < GREAT )
            {
                minAcceptedVolRatio =
                    Foam::min
                    (
                        minAcceptedVolRatio,
                        candidateMinRatio
                    );
            }

            maxAcceptedDist =
                Foam::max
                (
                    maxAcceptedDist,
                    Foam::sqrt
                    (
                        Foam::max
                        (
                            actualVirtualMoveDistSq,
                            scalar(0)
                        )
                    )
                );
        }


        for(label patchA=0; patchA<nSfgPatches; ++patchA)
        {
            for
            (
                label patchB=patchA+1;
                patchB<nSfgPatches;
                ++patchB
            )
            {
                const label pairI =
                    patchA*nSfgPatches + patchB;

                if
                (
                    acceptedByPair[pairI] == 0 &&
                    faceRejectedByPair[pairI] == 0 &&
                    volumeRejectedByPair[pairI] == 0
                )
                {
                    continue;
                }

                Info
                    << "[SOURCE_FEATURE_VOLUME_GATE]"
                    << " pair="
                    << sfgPatchNames[patchA]
                    << "|"
                    << sfgPatchNames[patchB]
                    << " accepted="
                    << acceptedByPair[pairI]
                    << " faceRejected="
                    << faceRejectedByPair[pairI]
                    << " volumeRejected="
                    << volumeRejectedByPair[pairI]
                    << endl;
            }
        }

        Info
            << "[SOURCE_FEATURE_VOLUME_GATE_SUMMARY]"
            << " edgeNotFoundSeen="
            << nEdgeNotFoundSeen
            << " virtualCandidate="
            << nVirtualCandidate
            << " accepted="
            << nAccepted
            << " graphUnresolved="
            << nGraphUnresolved
            << " nativeDeferred="
            << nNativeDeferred
            << " rejectFaceOnly="
            << nRejectedFaceOnly
            << " rejectVolumeOnly="
            << nRejectedVolumeOnly
            << " rejectFaceAndVolume="
            << nRejectedFaceAndVolume
            << " rejectPreExistingBad="
            << nRejectedPreExistingBad
            << " rejectNoCells="
            << nRejectedNoCells
            << " rejectActualMoveDistance="
            << nRejectedActualMoveDistance
            << " minCandidateVolRatio="
            <<
            (
                minCandidateVolRatio < GREAT
              ? minCandidateVolRatio
              : scalar(-1)
            )
            << " minAcceptedVolRatio="
            <<
            (
                minAcceptedVolRatio < GREAT
              ? minAcceptedVolRatio
              : scalar(-1)
            )
            << " maxAcceptedDist="
            << maxAcceptedDist
            << endl;
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
