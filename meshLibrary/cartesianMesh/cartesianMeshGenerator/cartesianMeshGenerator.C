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

#include "cartesianMeshGenerator.H"
#include "triSurf.H"
#include "triSurfacePatchManipulator.H"
#include "demandDrivenData.H"
#include "meshOctreeCreator.H"
#include "cartesianMeshExtractor.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceEngineModifier.H"
#include "polyMeshGenAddressing.H"
#include "boundaryLayers.H"
#include "meshSurfaceMapper.H"
#include "edgeExtractor.H"
#include "meshSurfaceEdgeExtractorNonTopo.H"
#include "meshOptimizer.H"
#include "meshSurfaceOptimizer.H"
#include "topologicalCleaner.H"
#include "boundaryLayers.H"
#include "refineBoundaryLayers.H"
#include "renameBoundaryPatches.H"
#include "checkMeshDict.H"
#include "checkCellConnectionsOverFaces.H"
#include "checkIrregularSurfaceConnections.H"
#include "checkNonMappableCellConnections.H"
#include "OFstream.H"
#include "checkBoundaryFacesSharingTwoEdges.H"
#include "triSurfaceMetaData.H"
#include "polyMeshGenChecks.H"
#include "triSurfaceCleanupDuplicates.H"
#include "polyMeshGenGeometryModification.H"
#include "surfaceMeshGeometryModification.H"

//#define DEBUG

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private member functions  * * * * * * * * * * * * //

// Non-mutating scan for near-coincident vertices
// Returns count of unique vertex pairs closer than tolerance
static label scanNearCoincidentPoints
(
    const triSurf& surf,
    const meshOctree& octree,
    const scalar tolerance,
    const dictionary& meshDict,
    const bool verbose = false
)
{
    const pointField& pts = surf.points();
    const wordList pNames = surf.patchNames();

    // Build BL patch index set from meshDict patchBoundaryLayers
    labelHashSet blPatchIndices;
    if( meshDict.isDict("boundaryLayers") )
    {
        const dictionary& bndL = meshDict.subDict("boundaryLayers");
        if( bndL.isDict("patchBoundaryLayers") )
        {
            const dictionary& pbl = bndL.subDict("patchBoundaryLayers");
            forAll(pNames, pi)
            {
                if( pbl.isDict(pNames[pi]) )
                {
                    const dictionary& pd = pbl.subDict(pNames[pi]);
                    const label nL = pd.found("nLayers") ?
                        readLabel(pd.lookup("nLayers")) : 0;
                    if( nL > 0 )
                        blPatchIndices.insert(pi);
                }
            }
        }
        // Also check global nLayers
        if( bndL.found("nLayers") )
        {
            const label globalN = readLabel(bndL.lookup("nLayers"));
            if( globalN > 0 )
                forAll(pNames, pi)
                    blPatchIndices.insert(pi);
        }
    }
    std::set<std::pair<label,label>> countedPairs;
    label nPairs = 0;

    // Build point-to-patch membership once
    List<labelHashSet> pointPatches(surf.nPoints());
    forAll(surf, triI)
    {
        const labelledTri& tri = surf[triI];
        const label region = tri.region();
        forAll(tri, vi)
        {
            const label pI = tri[vi];
            if( pI >= 0 && pI < label(pointPatches.size()) )
                pointPatches[pI].insert(region);
        }
    }

    for(label leafI=0; leafI<octree.numberOfLeaves(); ++leafI)
    {
        DynList<label> ct;
        octree.containedTriangles(leafI, ct);

        std::set<label> points;
        forAll(ct, ctI)
        {
            const label triI = ct[ctI];
            if( triI < 0 || triI >= label(surf.size()) ) continue;
            const labelledTri& tri = surf[triI];
            forAll(tri, i)
                points.insert(tri[i]);
        }

        for
        (
            std::set<label>::const_iterator it = points.begin();
            it != points.end();
            ++it
        )
        {
            const label pointI = *it;
            std::set<label>::const_iterator nIt = it;
            ++nIt;
            for(; nIt != points.end(); ++nIt)
            {
                const label pointJ = *nIt;
                if( magSqr(pts[pointI] - pts[pointJ]) < sqr(tolerance) )
                {
                    const label a = Foam::min(pointI, pointJ);
                    const label b = Foam::max(pointI, pointJ);
                    if( countedPairs.insert(std::make_pair(a, b)).second )
                    {
                        ++nPairs;

                        // Collect patch names -- always, not just verbose
                        DynList<word> pNamesA, pNamesB;
                        forAllConstIter(labelHashSet, pointPatches[a], it2)
                        {
                            const label r = it2.key();
                            if( r >= 0 && r < label(pNames.size()) )
                                pNamesA.append(pNames[r]);
                        }
                        forAllConstIter(labelHashSet, pointPatches[b], it2)
                        {
                            const label r = it2.key();
                            if( r >= 0 && r < label(pNames.size()) )
                                pNamesB.append(pNames[r]);
                        }

                        const scalar d = mag(pts[pointI]-pts[pointJ]);
                        const scalar ratio = d / tolerance;

                        // Conservative classification
                        word classification = "UNKNOWN";
                        bool samePatch = false;
                        forAll(pNamesA, pi)
                            forAll(pNamesB, pj)
                                if( pNamesA[pi] == pNamesB[pj] )
                                    samePatch = true;

                        // Check if both points are on BL wall patches
                        bool bothBLWall = false;
                        if( pNamesA.size() && pNamesB.size() )
                        {
                            bool aIsBL = false, bIsBL = false;
                            forAll(pNamesA, pii)
                            {
                                forAll(pNames, ni)
                                    if( pNames[ni] == pNamesA[pii]
                                     && blPatchIndices.found(ni) )
                                        aIsBL = true;
                            }
                            forAll(pNamesB, pii)
                            {
                                forAll(pNames, ni)
                                    if( pNames[ni] == pNamesB[pii]
                                     && blPatchIndices.found(ni) )
                                        bIsBL = true;
                            }
                            bothBLWall = aIsBL && bIsBL;
                        }

                        if( ratio < 0.20 )
                            // Feature <20% of tolerance -- almost certainly defect
                            classification = "WELD_SAFE_RATIO";
                        else if( !pNamesA.size() || !pNamesB.size() )
                            // Missing patch data -- never auto-weld
                            classification = "FLAG_INCOMPLETE";
                        else if( samePatch && ratio < 0.80 && pNamesA.size() && pNamesB.size() )
                            // Same patch, moderate ratio -- surface sliver
                            classification = "WELD_SAME_PATCH";
                        else if( bothBLWall && ratio < 0.50 )
                            // Both BL wall patches, moderate ratio -- TE knife-edge
                            classification = "WELD_BL_WALL_PAIR";
                        else if( ratio > 0.80 )
                            // High ratio -- may be legitimate geometry
                            classification = "FLAG_HIGH_RATIO";
                        else
                            // Cross-patch or unknown -- flag for review
                            classification = "FLAG_CROSS_PATCH_CANDIDATE";

                        if( verbose )
                        {
                            Info << "  Near-coincident pair: "
                                 << a << " " << b
                                 << " d=" << d*1000.0 << "mm"
                                 << " ratio=" << ratio
                                 << " class=" << classification
                                 << " at " << pts[pointI];
                            Info << " pA(" << pNamesA.size() << ")=(";
                            forAll(pNamesA, i) Info << pNamesA[i] << " ";
                            Info << ")";
                            Info << " pB(" << pNamesB.size() << ")=(";
                            forAll(pNamesB, i) Info << pNamesB[i] << " ";
                            Info << ")" << endl;
                        }
                    }
                }
            }
        }
    }
    return nPairs;
}

void cartesianMeshGenerator::createCartesianMesh()
{
    //- create polyMesh from octree boxes
    cartesianMeshExtractor cme(*octreePtr_, meshDict_, mesh_);

    if( meshDict_.found("decomposePolyhedraIntoTetsAndPyrs") )
    {
        if( readBool(meshDict_.lookup("decomposePolyhedraIntoTetsAndPyrs")) )
            cme.decomposeSplitHexes();
    }

    cme.createMesh();
}

void cartesianMeshGenerator::surfacePreparation()
{
    //- removes unnecessary cells and morph the boundary
    //- such that there is only one boundary face per cell
    //- It also checks topology of cells after morphing is performed
    bool changed;

    do
    {
        changed = false;

        checkIrregularSurfaceConnections checkConnections(mesh_);
        if( checkConnections.checkAndFixIrregularConnections() )
            changed = true;

        if( checkNonMappableCellConnections(mesh_).removeCells() )
            changed = true;

        if( checkCellConnectionsOverFaces(mesh_).checkCellGroups() )
            changed = true;
    } while( changed );

    checkBoundaryFacesSharingTwoEdges(mesh_).improveTopology();
}

void cartesianMeshGenerator::mapMeshToSurface()
{
    //- calculate mesh surface
    meshSurfaceEngine mse(mesh_);

    //- pre-map mesh surface
    meshSurfaceMapper mapper(mse, *octreePtr_);
    mapper.preMapVertices(0);

    //- map mesh surface on the geometry surface
    mapper.mapVerticesOntoSurface();
    //- targeted repair of validity-rejected points before corner snap
    mapper.repairRejectedPoints();

    //- snap corner and edge vertices onto feature curves
    //- early pass: stabilizes features before surface optimizer runs
    mapper.mapCornersAndEdges();

    //- constrained surface smoothing: redistribute single-patch
    //- points around snapped features before untangling
    mapper.smoothSinglePatchPoints(3);

    //- untangle surface faces
    meshSurfaceOptimizer(mse, *octreePtr_).untangleSurface();
}

void cartesianMeshGenerator::extractPatches()
{
    edgeExtractor extractor(mesh_, *octreePtr_);

    Info << "Extracting edges" << endl;
    extractor.extractEdges();

    extractor.updateMeshPatches();
}

void cartesianMeshGenerator::mapEdgesAndCorners()
{
    if( !blNoBlEdgePoints_.empty() || !blNeutralEdgePoints_.empty() )
    {
        meshSurfaceEdgeExtractorNonTopo
        (
            mesh_,
            *octreePtr_,
            blNoBlEdgePoints_,
            blNoBlPointPatch_
        );
    }
    else
    {
        meshSurfaceEdgeExtractorNonTopo(mesh_, *octreePtr_);
    }
}

void cartesianMeshGenerator::optimiseMeshSurface()
{
    meshSurfaceEngine mse(mesh_);
    meshSurfaceOptimizer(mse, *octreePtr_).optimizeSurface();
}

void cartesianMeshGenerator::detectGapPoints
(
    boundaryLayers& bl
)
{
    if( !meshDict_.isDict("boundaryLayers") )
        return;
    const dictionary& bndL = meshDict_.subDict("boundaryLayers");
    if( !bndL.found("detectGaps") )
        return;
    if( !Switch(bndL.lookup("detectGaps")) )
        return;
    if( !bndL.found("gapPatchPairs") )
    {
        Info << "Gap detection: no gapPatchPairs defined, skipping" << endl;
        return;
    }

    const scalar gapThreshold =
        bndL.found("gapThreshold") ?
        readScalar(bndL.lookup("gapThreshold")) : 5e-4;

    // Read explicit patch pairs
    const List<Pair<word>> pairList(bndL.lookup("gapPatchPairs"));

    Info << "Gap detection: threshold=" << gapThreshold
         << " m, " << pairList.size() << " patch pairs" << endl;

    if( !octreePtr_ )
    {
        Info << "Gap detection: octreePtr_ is null, skipping" << endl;
        return;
    }

    const triSurf& surf = octreePtr_->surface();
    const wordList pNames = surf.patchNames();

    // Build word->index map
    Map<label> patchNameToIdx;
    forAll(pNames, pi)
        patchNameToIdx.insert(pi, pi);
    // Rebuild as name->index
    HashTable<label> nameToIdx;
    forAll(pNames, pi)
        nameToIdx.insert(pNames[pi], pi);

    // Build set of pairs as index pairs
    List<Pair<label>> idxPairs;
    forAll(pairList, pI)
    {
        const word& nameA = pairList[pI].first();
        const word& nameB = pairList[pI].second();
        if( !nameToIdx.found(nameA) )
        {
            Info << "Gap detection: patch " << nameA << " not found, skipping pair" << endl;
            continue;
        }
        if( !nameToIdx.found(nameB) )
        {
            Info << "Gap detection: patch " << nameB << " not found, skipping pair" << endl;
            continue;
        }
        idxPairs.append(Pair<label>(nameToIdx[nameA], nameToIdx[nameB]));
        Info << "Gap pair: " << nameA << " <-> " << nameB << endl;
    }

    if( idxPairs.empty() )
    {
        Info << "Gap detection: no valid patch pairs, skipping" << endl;
        return;
    }

    meshSurfaceEngine mse(mesh_);
    const labelList& bPoints = mse.boundaryPoints();
    const pointFieldPMG& points = mesh_.points();
    const meshSurfacePartitioner mPart(mse);
    const VRWGraph& pPatches = mPart.pointPatches();

    // Gap conflict arbitration: determine loser side for each pair BEFORE scan.
    // Mode autoThickness: loser = patch with higher requested total BL thickness.
    // Mode manualSuppressPatches: loser = user-provided gapSuppressPatches list.
    // The same loser set controls both gap-point seeding and ring1 suppression.
    const word conflictMode =
        bndL.found("gapConflictMode") ?
        word(bndL.lookup("gapConflictMode")) : word("autoThickness");

    // Helper: compute requested total BL thickness for a named patch from meshDict.
    // total = maxFirstLayerThickness * (1 - ratio^nLayers) / (1 - ratio)
    auto requestedThickness = [&](const word& pName) -> scalar
    {
        scalar firstT = 1e-4;
        scalar ratio  = 1.3;
        label  nL     = 0;
        if( bndL.found("maxFirstLayerThickness") )
            firstT = readScalar(bndL.lookup("maxFirstLayerThickness"));
        if( bndL.found("thicknessRatio") )
            ratio = readScalar(bndL.lookup("thicknessRatio"));
        if( bndL.found("nLayers") )
            nL = readLabel(bndL.lookup("nLayers"));
        if( bndL.isDict("patchBoundaryLayers") )
        {
            const dictionary& pbl = bndL.subDict("patchBoundaryLayers");
            if( pbl.isDict(pName) )
            {
                const dictionary& pd = pbl.subDict(pName);
                if( pd.found("maxFirstLayerThickness") )
                    firstT = readScalar(pd.lookup("maxFirstLayerThickness"));
                if( pd.found("thicknessRatio") )
                    ratio = readScalar(pd.lookup("thicknessRatio"));
                if( pd.found("nLayers") )
                    nL = readLabel(pd.lookup("nLayers"));
            }
        }
        if( nL <= 0 ) return 0.0;
        if( mag(ratio - 1.0) < SMALL )
            return firstT * scalar(nL);
        return firstT * (1.0 - Foam::pow(ratio, scalar(nL))) / (1.0 - ratio);
    };

    labelHashSet loserPatches;
    DynList<word> loserPatchNamesDyn;  // built at detection time from STL names
    if( conflictMode == "manualSuppressPatches" )
    {
        // Legacy manual mode: gapSuppressPatches list is the loser side.
        if( bndL.found("gapSuppressPatches") )
        {
            wordList suppressNames(bndL.lookup("gapSuppressPatches"));
            forAll(suppressNames, si)
                if( nameToIdx.found(suppressNames[si]) )
                {
                    loserPatches.insert(nameToIdx[suppressNames[si]]);
                    loserPatchNamesDyn.appendIfNotIn(suppressNames[si]);
                }
            Info << "Gap conflict: manual suppress-side patches: "
                 << suppressNames << endl;
        }
    }
    else
    {
        // autoThickness: loser = higher requested BL thickness per pair.
        // Uses idxPairs (already validated) -> pNames, not pairList, to avoid
        // index mismatch when pairs are skipped during validation.
        forAll(idxPairs, pI)
        {
            const label idxA = idxPairs[pI].first();
            const label idxB = idxPairs[pI].second();
            if( idxA < 0 || idxA >= label(pNames.size()) ) continue;
            if( idxB < 0 || idxB >= label(pNames.size()) ) continue;
            const word& nameA = pNames[idxA];
            const word& nameB = pNames[idxB];
            const scalar tA = requestedThickness(nameA);
            const scalar tB = requestedThickness(nameB);
            if( tA <= SMALL && tB <= SMALL )
            {
                Info << "Gap conflict arbitration: pair (" << nameA
                     << " <-> " << nameB
                     << ") zero BL thickness both sides; skipping" << endl;
                continue;
            }
            const label loserIdx  = (tA >= tB) ? idxA  : idxB;
            const word& loserName = (tA >= tB) ? nameA : nameB;
            loserPatches.insert(loserIdx);
            loserPatchNamesDyn.appendIfNotIn(loserName);
            Info << "Gap conflict arbitration: pair (" << nameA
                 << " <-> " << nameB
                 << ") thickness=(" << tA << " " << tB << ")"
                 << " loser=" << loserName
                 << " [autoThickness]" << endl;
        }
        if( loserPatches.size() == 0 )
            Info << "Gap conflict arbitration: WARNING no losers found"
                 << " -- gap points will not be loser-side gated" << endl;
    }

    labelHashSet gapPoints;
    label nScanned = 0;
    label nGapHits = 0;

    forAll(bPoints, bpI)
    {
        const point& pt = points[bPoints[bpI]];
        bool isGapPoint = false;

        forAll(idxPairs, pairI)
        {
            const label pIdxA = idxPairs[pairI].first();
            const label pIdxB = idxPairs[pairI].second();

            bool onA = false, onB = false;
            forAllRow(pPatches, bpI, pI)
            {
                if( pPatches(bpI, pI) == pIdxA ) onA = true;
                if( pPatches(bpI, pI) == pIdxB ) onB = true;
            }
            if( !onA && !onB ) continue;

            // Gap detection is symmetric: seed gap points on both sides of
            // the narrow clearance. loserPatches_ is retained as side metadata
            // for later loser-side layer capping/suppression, but it must not
            // shrink the detected geometric danger zone.

            ++nScanned;
            const label searchPatch = onA ? pIdxB : pIdxA;

            point nearest;
            scalar distSq = GREAT;
            label nearestTri = -1;
            octreePtr_->findNearestSurfacePointInRegion
            (
                nearest, distSq, nearestTri, searchPatch, pt
            );

            if( distSq < sqr(gapThreshold) )
            {
                isGapPoint = true;
                ++nGapHits;
                break;
            }
        }

        if( isGapPoint )
            gapPoints.insert(bPoints[bpI]);
    }

    Info << "Gap detection: scanned " << nScanned
         << " candidate points, found " << gapPoints.size()
         << " gap points from " << nGapHits << " hits" << endl;

    // Pass symmetric zone points for classification (both sides).
    if( gapPoints.size() > 0 )
        bl.setGapZonePoints(gapPoints);

    // Pass loser-side action points for BL suppression.
    // gapPoints_ in boundaryLayers means "action/suppression points" --
    // not the full symmetric geometric zone.
    if( gapPoints.size() > 0 && loserPatches.size() > 0 )
    {
        labelHashSet loserGapPoints;

        // Reuse already-built reverse map: mesh point -> boundary point index.
        labelList meshToBnd(mesh_.points().size(), -1);
        forAll(bPoints, bpI)
            meshToBnd[bPoints[bpI]] = bpI;

        forAllConstIter(labelHashSet, gapPoints, it)
        {
            const label meshPtI = it.key();
            if( meshPtI < 0 || meshPtI >= label(meshToBnd.size()) ) continue;
            const label bpI = meshToBnd[meshPtI];
            if( bpI < 0 ) continue;
            forAllRow(pPatches, bpI, pI)
            {
                if( loserPatches.found(pPatches(bpI, pI)) )
                {
                    loserGapPoints.insert(meshPtI);
                    break;
                }
            }
        }

        if( loserGapPoints.size() > 0 )
            bl.setGapPoints(loserGapPoints);

        Info << "Gap action points: " << loserGapPoints.size()
             << " loser-side of " << gapPoints.size()
             << " zone points" << endl;
    }
    else if( gapPoints.size() > 0 )
    {
        // No arbitration data: fall back to legacy symmetric action.
        bl.setGapPoints(gapPoints);
    }

    if( loserPatches.size() > 0 )
    {
        bl.setGapLoserPatches(loserPatches);
        // Pass names built at detection time -- avoids STL/polyMesh index mismatch.
        wordList loserPatchNames(loserPatchNamesDyn.size());
        forAll(loserPatchNamesDyn, i)
            loserPatchNames[i] = loserPatchNamesDyn[i];
        bl.setGapLoserPatchNames(loserPatchNames);
        Info << "Gap loser patch names: " << loserPatchNames << endl;
    }
}

void cartesianMeshGenerator::detectTripleJunctions
(
    boundaryLayers& bl
)
{
    if( !meshDict_.isDict("boundaryLayers") )
        return;
    const dictionary& bndL = meshDict_.subDict("boundaryLayers");
    if( !bndL.found("tripleJunctionSuppressPatches") )
        return;
    if( !bndL.found("tripleJunctionWallPatches") )
        return;
    if( !bndL.found("tripleJunctionNeutralPatches") )
        return;

    const wordList suppressNames(bndL.lookup("tripleJunctionSuppressPatches"));
    const wordList wallNames(bndL.lookup("tripleJunctionWallPatches"));
    const wordList neutralNames(bndL.lookup("tripleJunctionNeutralPatches"));

    if( suppressNames.empty() || wallNames.empty() || neutralNames.empty() )
        return;

    if( !octreePtr_ )
    {
        Info << "Triple-junction detection: octreePtr_ null, skipping" << endl;
        return;
    }

    const triSurf& surf = octreePtr_->surface();
    const wordList pNames = surf.patchNames();

    HashTable<label> nameToIdx;
    forAll(pNames, pi)
        nameToIdx.insert(pNames[pi], pi);

    labelHashSet suppressIdx;
    forAll(suppressNames, si)
        if( nameToIdx.found(suppressNames[si]) )
            suppressIdx.insert(nameToIdx[suppressNames[si]]);

    labelHashSet wallIdx;
    forAll(wallNames, wi)
        if( nameToIdx.found(wallNames[wi]) )
            wallIdx.insert(nameToIdx[wallNames[wi]]);

    labelHashSet neutralIdx;
    forAll(neutralNames, ni)
        if( nameToIdx.found(neutralNames[ni]) )
            neutralIdx.insert(nameToIdx[neutralNames[ni]]);

    if( suppressIdx.empty() || wallIdx.empty() || neutralIdx.empty() )
    {
        Info << "Triple-junction detection: one or more patch sets not found "
             << "in surface -- skipping" << endl;
        return;
    }

    meshSurfaceEngine mse(mesh_);
    const labelList& bPoints = mse.boundaryPoints();
    const meshSurfacePartitioner mPart(mse);
    const VRWGraph& pPatches = mPart.pointPatches();

    labelHashSet triplePoints;

    forAll(bPoints, bpI)
    {
        if( pPatches.sizeOfRow(bpI) < 3 ) continue;

        bool onSuppress = false;
        bool onWall     = false;
        bool onNeutral  = false;

        forAllRow(pPatches, bpI, pI)
        {
            const label pIdx = pPatches(bpI, pI);
            if( suppressIdx.found(pIdx) ) onSuppress = true;
            if( wallIdx.found(pIdx) )     onWall     = true;
            if( neutralIdx.found(pIdx) )  onNeutral  = true;
        }

        if( onSuppress && onWall && onNeutral )
            triplePoints.insert(bPoints[bpI]);
    }

    Info << "Triple-junction detection: found " << triplePoints.size()
         << " triple-junction points"
         << " (suppress=" << suppressNames
         << " wall=" << wallNames
         << " neutral=" << neutralNames << ")" << endl;

    // Store triple-junction seeds in boundaryLayers for planner + gated exclusion
    bl.addTripleJunctionPoints(triplePoints);
}

void cartesianMeshGenerator::generateBoundaryLayers()
{
    //- add boundary layers
    boundaryLayers bl(mesh_, meshDict_);
    bl.terminateLayersAtConcaveEdges();

    // Gap/proximity closure: detect tight BL/BL patch proximity and suppress
    // BL locally before createNewVertices builds the prism graph.
    detectGapPoints(bl);
    detectTripleJunctions(bl);
    bl.reportBLTransitionSeeds();
    bl.buildBLTransitionPlan();

    bl.addLayerForAllPatches();
    // Capture layerScale for post-replaceBoundaries coverage report
    blLayerScale_ = bl.layerScale();

    // Capture junction points for handoff to refineBoundaryLayers
    blblJunctionPoints_ = bl.junctionEdgePoints();
    blblAcuteCornerPoints_ = bl.blblAcuteCornerPoints();
    vtFaceRing_ = bl.vtFaceRing();
    blGapActionPoints_    = bl.gapPoints();
    blGapLoserPatchNames_ = bl.gapLoserPatchNames();
    Info << "Acute BL+BL+neutral corners captured: "
         << blblAcuteCornerPoints_.size() << endl;

    // Capture BL/no-BL transition edge points (boundary-point indices)
    // for handoff to post-BL mapper instances
    blNoBlEdgePoints_ = bl.blNoBlEdgePoints();
    blNeutralEdgePoints_ = bl.blNeutralEdgePoints();
    blNeutralPointPatch_ = bl.blNeutralPointPatch();
    blNoBlPointPatch_ = bl.blNoBlPointPatch();
    Info << "BL/no-BL edge points captured for mapper exclusion: "
         << blNoBlEdgePoints_.size() << endl;
}

void cartesianMeshGenerator::refBoundaryLayers()
{
    if( meshDict_.isDict("boundaryLayers") )
    {
        refineBoundaryLayers refLayers(mesh_);

        refineBoundaryLayers::readSettings(meshDict_, refLayers);

        // Pass BL/BL junction points for wedge topology
        refLayers.setBlblJunctionPoints(blblJunctionPoints_);
        refLayers.setBlblAcuteCornerPoints(blblAcuteCornerPoints_);
        refLayers.setVtFaceRing(vtFaceRing_);
        // Pass gap conflict data for Option B local split-edge layer capping.
        // Uses stable mesh-point labels + patch names -- no face-index mismatch.
        if( blGapActionPoints_.size() > 0 )
        {
            refLayers.setGapActionPoints(blGapActionPoints_);
            refLayers.setGapLoserPatchNames(blGapLoserPatchNames_);
            // Read ring max-layer knobs from meshDict
            label ring1Max = 1;
            label ring2Max = 2;
            if( meshDict_.isDict("boundaryLayers") )
            {
                const dictionary& bndL = meshDict_.subDict("boundaryLayers");
                if( bndL.found("gapLoserRing1MaxLayers") )
                    ring1Max = readLabel(bndL.lookup("gapLoserRing1MaxLayers"));
                if( bndL.found("gapLoserRing2MaxLayers") )
                    ring2Max = readLabel(bndL.lookup("gapLoserRing2MaxLayers"));
            }
            refLayers.setGapRingMaxLayers(ring1Max, ring2Max);
        }

        // Pass BL/termination edge points for inlet/outlet layer-count cap.
        // IMPORTANT: blNoBlEdgePoints_ is stored as boundary-point indices (bpI).
        // refineBoundaryLayers splitEdges use global mesh point labels, so convert
        // bpI -> mesh point label here. Do not change blNoBlEdgePoints_ globally,
        // because mapper/projection paths may expect bpI indexing.
        if( !blNoBlEdgePoints_.empty() )
        {
            const meshSurfaceEngine mseTerm(mesh_);
            const labelList& bPointsTerm = mseTerm.boundaryPoints();
            labelHashSet blTerminationMeshPoints;
            label nBadBp = 0;
            forAllConstIter(labelHashSet, blNoBlEdgePoints_, it)
            {
                const label bpI = it.key();
                if( bpI >= 0 && bpI < label(bPointsTerm.size()) )
                    blTerminationMeshPoints.insert(bPointsTerm[bpI]);
                else
                    ++nBadBp;
            }
            if( blTerminationMeshPoints.size() > 0 )
                refLayers.setBlTerminationEdgePoints(blTerminationMeshPoints);
            label tRing1Max = 3;  // disabled by default
            label tRing2Max = 3;
            if( meshDict_.isDict("boundaryLayers") )
            {
                const dictionary& bndL = meshDict_.subDict("boundaryLayers");
                if( bndL.found("blTerminationRing1MaxLayers") )
                    tRing1Max = readLabel(bndL.lookup("blTerminationRing1MaxLayers"));
                if( bndL.found("blTerminationRing2MaxLayers") )
                    tRing2Max = readLabel(bndL.lookup("blTerminationRing2MaxLayers"));
            }
            refLayers.setBlTerminationRingMaxLayers(tRing1Max, tRing2Max);
            Info << "BL termination edge cap: "
                 << blNoBlEdgePoints_.size() << " bp-index edge points, "
                 << blTerminationMeshPoints.size() << " mesh points"
                 << " ring1max=" << tRing1Max
                 << " ring2max=" << tRing2Max
                 << " badBp=" << nBadBp << endl;
        }
        {
            bool capLayers = false;
            if( meshDict_.isDict("boundaryLayers") )
            {
                const dictionary& bndL = meshDict_.subDict("boundaryLayers");
                if( bndL.found("acuteCornerCapLayers") )
                    capLayers = bool(Switch(bndL.lookup("acuteCornerCapLayers")));
            }
            refLayers.setAcuteCornerCapLayers(capLayers);
            Info << "Acute corner face cap: " << (capLayers ? "enabled" : "disabled") << endl;
        }

        nPointsBeforeBL_ = mesh_.points().size();
        refLayers.refineLayers();

        refLayers.pointsInBndLayer(blPoints_);

        {
            mesh_.clearAddressingData();
            const bool hadUnusedBefore =
                polyMeshGenChecks::checkPoints(mesh_, false);
            if( hadUnusedBefore )
            {
                labelHashSet negBefore;
                polyMeshGenChecks::checkCellVolumes(mesh_, false, &negBefore);
                polyMeshGenModifier(mesh_).removeUnusedVertices();
                mesh_.clearAddressingData();
                labelHashSet negAfter;
                polyMeshGenChecks::checkCellVolumes(mesh_, false, &negAfter);
                const bool hasUnusedAfter =
                    polyMeshGenChecks::checkPoints(mesh_, false);
                Info << "Post-BL cleanup: removeUnusedVertices unusedPoints "
                     << "bad->" << (hasUnusedAfter ? "bad" : "ok")
                     << " negVol " << negBefore.size() << "->" << negAfter.size()
                     << endl;
                if( negAfter.size() > negBefore.size() )
                    Info << "Post-BL cleanup warning: removeUnusedVertices "
                         << "worsened negative-volume count" << endl;
            }
            else
            {
                Info << "Post-BL cleanup: no unused vertices found" << endl;
            }
        }


        // 3-gate post-refinement diagnostic (non-mutating)
        {
            Pout << nl
                 << "### ENTERING 3-GATE POST-REFINEMENT DIAGNOSTIC ###"
                 << nl << endl;
            Pout << "3-gate diagnostic: checking post-refinement BL quality" << endl;

            labelHashSet badCells;
            labelHashSet badPyramidFaces;
            labelHashSet nonOrthoFaces;

            const bool hasNegVol =
                polyMeshGenChecks::checkCellVolumes(mesh_, false, &badCells);

            const bool hasBadPyramids =
                polyMeshGenChecks::checkFacePyramids
                (mesh_, false, -SMALL, &badPyramidFaces);

            // Gate 3: severe non-orthogonality above 85 degrees
            const bool hasNonOrtho =
                polyMeshGenChecks::checkFaceDotProduct
                (mesh_, false, 85.0, &nonOrthoFaces);

            Pout << "3-gate diagnostic results:" << endl;
            Pout << "  Gate 1 (neg vol cells):    " << badCells.size() << endl;
            Pout << "  Gate 2 (bad pyramids):     " << badPyramidFaces.size() << endl;
            Pout << "  Gate 3 (non-ortho >85deg): " << nonOrthoFaces.size() << endl;

            // Final bad-pyramid classifier. Non-mutating diagnostic only.
            // This classifies the remaining post-refinement bad faces by
            // nearby patch contact so we know whether the residual defect is
            // periodic-transition, blade-wall, hub/shroud, or generic volume.
            if( badPyramidFaces.size() > 0 )
            {
                const faceListPMG& faces = mesh_.faces();
                const pointFieldPMG& points = mesh_.points();
                const labelList& owner = mesh_.owner();
                const labelList& neighbour = mesh_.neighbour();
                const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

                boolList cellTouchesPeriodic(mesh_.cells().size(), false);
                boolList cellTouchesBlade(mesh_.cells().size(), false);
                boolList cellTouchesHub(mesh_.cells().size(), false);
                boolList cellTouchesShroud(mesh_.cells().size(), false);
                boolList cellTouchesInletOutlet(mesh_.cells().size(), false);

                forAll(boundaries, patchI)
                {
                    const word& pName = boundaries[patchI].patchName();

                    const bool isPeriodic =
                        (pName == "periodic_1" || pName == "periodic_2");
                    const bool isBlade = pName.find("blade") == 0;
                    const bool isHub = (pName == "hub");
                    const bool isShroud = (pName == "shroud");
                    const bool isInletOutlet =
                        (pName == "inlet" || pName == "outlet");

                    if( !isPeriodic && !isBlade && !isHub
                     && !isShroud && !isInletOutlet )
                        continue;

                    const label start = boundaries[patchI].patchStart();
                    const label size  = boundaries[patchI].patchSize();

                    for(label pfI=0; pfI<size; ++pfI)
                    {
                        const label faceJ = start + pfI;
                        if( faceJ < 0 || faceJ >= label(owner.size()) )
                            continue;

                        const label cellI = owner[faceJ];
                        if( cellI < 0 || cellI >= label(cellTouchesBlade.size()) )
                            continue;

                        if( isPeriodic )    cellTouchesPeriodic[cellI] = true;
                        if( isBlade )       cellTouchesBlade[cellI] = true;
                        if( isHub )         cellTouchesHub[cellI] = true;
                        if( isShroud )      cellTouchesShroud[cellI] = true;
                        if( isInletOutlet ) cellTouchesInletOutlet[cellI] = true;
                    }
                }

                label nPeriodicClass(0);
                label nBladeClass(0);
                label nHubClass(0);
                label nShroudClass(0);
                label nInletOutletClass(0);
                label nGenericClass(0);

                // Candidate subset for future localized Gate 2 repair.
                // These are final bad-pyramid faces adjacent to periodic/
                // neutral patches, separated from the mixed global bad set.
                labelHashSet gate2PeriodicBadFaces;

                label nPrinted(0);

                forAllConstIter(labelHashSet, badPyramidFaces, it)
                {
                    const label faceI = it.key();
                    if( faceI < 0 || faceI >= label(faces.size()) )
                        continue;

                    const label ownCell =
                        faceI < label(owner.size()) ? owner[faceI] : -1;
                    const label neiCell =
                        faceI < label(neighbour.size()) ? neighbour[faceI] : -1;

                    const bool touchPeriodic =
                        (ownCell >= 0 && ownCell < label(cellTouchesPeriodic.size())
                      && cellTouchesPeriodic[ownCell])
                     || (neiCell >= 0 && neiCell < label(cellTouchesPeriodic.size())
                      && cellTouchesPeriodic[neiCell]);

                    const bool touchBlade =
                        (ownCell >= 0 && ownCell < label(cellTouchesBlade.size())
                      && cellTouchesBlade[ownCell])
                     || (neiCell >= 0 && neiCell < label(cellTouchesBlade.size())
                      && cellTouchesBlade[neiCell]);

                    const bool touchHub =
                        (ownCell >= 0 && ownCell < label(cellTouchesHub.size())
                      && cellTouchesHub[ownCell])
                     || (neiCell >= 0 && neiCell < label(cellTouchesHub.size())
                      && cellTouchesHub[neiCell]);

                    const bool touchShroud =
                        (ownCell >= 0 && ownCell < label(cellTouchesShroud.size())
                      && cellTouchesShroud[ownCell])
                     || (neiCell >= 0 && neiCell < label(cellTouchesShroud.size())
                      && cellTouchesShroud[neiCell]);

                    const bool touchInletOutlet =
                        (ownCell >= 0 && ownCell < label(cellTouchesInletOutlet.size())
                      && cellTouchesInletOutlet[ownCell])
                     || (neiCell >= 0 && neiCell < label(cellTouchesInletOutlet.size())
                      && cellTouchesInletOutlet[neiCell]);

                    if( touchPeriodic )
                    {
                        ++nPeriodicClass;
                        gate2PeriodicBadFaces.insert(faceI);
                    }
                    else if( touchBlade ) ++nBladeClass;
                    else if( touchHub ) ++nHubClass;
                    else if( touchShroud ) ++nShroudClass;
                    else if( touchInletOutlet ) ++nInletOutletClass;
                    else ++nGenericClass;

                    if( nPrinted < 25 )
                    {
                        const face& f = faces[faceI];

                        point c(point::zero);
                        forAll(f, fpI)
                        {
                            const label pI = f[fpI];
                            if( pI >= 0 && pI < label(points.size()) )
                                c += points[pI];
                        }

                        if( f.size() > 0 )
                            c /= scalar(f.size());

                        const scalar r = Foam::sqrt(c.x()*c.x() + c.y()*c.y());
                        const scalar theta =
                            Foam::atan2(c.y(), c.x()) * 180.0 / constant::mathematical::pi;

                        Pout << "  Gate2 badFace faceI=" << faceI
                             << " nPts=" << f.size()
                             << " owner=" << ownCell
                             << " neighbour=" << neiCell
                             << " centroid=" << c
                             << " r=" << r
                             << " theta=" << theta
                             << " class="
                             << (touchPeriodic ? "periodic" :
                                 touchBlade ? "blade" :
                                 touchHub ? "hub" :
                                 touchShroud ? "shroud" :
                                 touchInletOutlet ? "inletOutlet" : "generic")
                             << endl;

                        ++nPrinted;
                    }
                }

                Pout << "  Gate2 bad pyramid classes:"
                     << " periodic=" << nPeriodicClass
                     << " blade=" << nBladeClass
                     << " hub=" << nHubClass
                     << " shroud=" << nShroudClass
                     << " inletOutlet=" << nInletOutletClass
                     << " generic=" << nGenericClass
                     << endl;

                Pout << "  Gate2 periodic-local repair candidate faces: "
                     << gate2PeriodicBadFaces.size()
                     << endl;

                // Gate 2 final repair pass. This is intentionally motion-only:
                // snapshot points, attempt a local low-quality optimization,
                // then accept only if bad pyramids improve and no hard quality
                // metric regresses. Do not flip faces here; remaining Gate 2
                // defects include boundary faces and mixed transition classes.
                const bool gate2UnusedBefore =
                    polyMeshGenChecks::checkPoints(mesh_, false);

                labelHashSet gate2NegBefore;
                polyMeshGenChecks::checkCellVolumes(mesh_, false, &gate2NegBefore);

                labelHashSet gate2OpenBefore;
                polyMeshGenChecks::checkClosedCells(mesh_, false, 0.5, &gate2OpenBefore);

                labelHashSet gate2NonOrthoBefore;
                polyMeshGenChecks::checkFaceDotProduct
                (mesh_, false, 85.0, &gate2NonOrthoBefore);

                scalarField gate2SkewBefore;
                polyMeshGenChecks::checkFaceSkewness(mesh_, gate2SkewBefore);
                const scalar gate2MaxSkewBefore =
                    gate2SkewBefore.size() > 0
                  ? max(gate2SkewBefore)
                  : scalar(0.0);

                const pointField gate2PointsBefore(mesh_.points());

                meshOptimizer gate2Optimizer(mesh_);
                if( gate2PeriodicBadFaces.size() > 0 )
                {
                    const faceListPMG& allFaces = mesh_.faces();
                    const labelList& own = mesh_.owner();
                    const labelList& nei = mesh_.neighbour();
                    const cellListPMG& cells = mesh_.cells();
                    labelHashSet freeCells;
                    forAllConstIter(labelHashSet, gate2PeriodicBadFaces, it)
                    {
                        const label faceI = it.key();
                        freeCells.insert(own[faceI]);
                        if( faceI < label(nei.size()) && nei[faceI] >= 0 )
                            freeCells.insert(nei[faceI]);
                    }
                    labelHashSet freePoints;
                    forAllConstIter(labelHashSet, freeCells, cit)
                    {
                        const cell& c = cells[cit.key()];
                        forAll(c, fI)
                        {
                            const face& f = allFaces[c[fI]];
                            forAll(f, pI)
                                freePoints.insert(f[pI]);
                        }
                    }
                    labelLongList lockedPts;
                    const pointFieldPMG& pts = mesh_.points();
                    forAll(pts, pointI)
                        if( !freePoints.found(pointI) )
                            lockedPts.append(pointI);
                    gate2Optimizer.lockPoints(lockedPts);
                    Info << "Gate2 local repair: locking " << lockedPts.size()
                         << " points, freeing " << freePoints.size()
                         << " points in owner/neighbour cells" << endl;
                }

                if( gate2NegBefore.size() > 0 )
                {
                    Info << "Gate2 local repair: skipped -- "
                         << gate2NegBefore.size()
                         << " negVol cells present, optimizer unsafe" << endl;
                }
                else
                {
                    gate2Optimizer.optimizeLowQualityFaces(3);
                }

                const bool gate2UnusedAfter =
                    polyMeshGenChecks::checkPoints(mesh_, false);

                labelHashSet gate2BadAfter;
                polyMeshGenChecks::checkFacePyramids
                (mesh_, false, -SMALL, &gate2BadAfter);

                labelHashSet gate2NegAfter;
                polyMeshGenChecks::checkCellVolumes(mesh_, false, &gate2NegAfter);

                labelHashSet gate2OpenAfter;
                polyMeshGenChecks::checkClosedCells(mesh_, false, 0.5, &gate2OpenAfter);

                labelHashSet gate2NonOrthoAfter;
                polyMeshGenChecks::checkFaceDotProduct
                (mesh_, false, 85.0, &gate2NonOrthoAfter);

                scalarField gate2SkewAfter;
                polyMeshGenChecks::checkFaceSkewness(mesh_, gate2SkewAfter);
                const scalar gate2MaxSkewAfter =
                    gate2SkewAfter.size() > 0
                  ? max(gate2SkewAfter)
                  : scalar(0.0);

                const bool gate2RepairOK =
                    gate2BadAfter.size() < badPyramidFaces.size()
                 && (!gate2UnusedAfter || gate2UnusedBefore)
                 && gate2NegBefore.size() == 0
                 && gate2NegAfter.size() == 0
                 && gate2OpenAfter.size() == 0
                 && gate2NonOrthoAfter.size() <= gate2NonOrthoBefore.size()
                 && gate2MaxSkewAfter <= gate2MaxSkewBefore
                 && gate2MaxSkewAfter <= scalar(20.0);

                if( gate2RepairOK )
                {
                    Pout << "  Gate2 final repair accepted: badPyramids "
                         << badPyramidFaces.size() << "->" << gate2BadAfter.size()
                         << " unusedPoints " << (gate2UnusedBefore ? "bad" : "ok")
                         << "->" << (gate2UnusedAfter ? "bad" : "ok")
                         << " negVol " << gate2NegBefore.size() << "->" << gate2NegAfter.size()
                         << " openCells " << gate2OpenBefore.size() << "->" << gate2OpenAfter.size()
                         << " nonOrtho85 " << gate2NonOrthoBefore.size() << "->" << gate2NonOrthoAfter.size()
                         << " skew " << gate2MaxSkewBefore << "->" << gate2MaxSkewAfter
                         << endl;

                    badPyramidFaces = gate2BadAfter;
                    badCells = gate2NegAfter;
                }
                else
                {
                    Pout << "  Gate2 final repair rejected: badPyramids "
                         << badPyramidFaces.size() << "->" << gate2BadAfter.size()
                         << " unusedPoints " << (gate2UnusedBefore ? "bad" : "ok")
                         << "->" << (gate2UnusedAfter ? "bad" : "ok")
                         << " negVol " << gate2NegBefore.size() << "->" << gate2NegAfter.size()
                         << " openCells " << gate2OpenBefore.size() << "->" << gate2OpenAfter.size()
                         << " nonOrtho85 " << gate2NonOrthoBefore.size() << "->" << gate2NonOrthoAfter.size()
                         << " skew " << gate2MaxSkewBefore << "->" << gate2MaxSkewAfter
                         << " -- rolling back"
                         << endl;

                    polyMeshGenModifier gate2MeshModifier(mesh_);
                    pointFieldPMG& pts = gate2MeshModifier.pointsAccess();
                    pts = gate2PointsBefore;
                    mesh_.clearAddressingData();
                }
            }

            if( hasNegVol || hasBadPyramids || hasNonOrtho )
            {
                const meshSurfaceEngine mse(mesh_);
                const labelList& bFaceOwner = mse.faceOwners();
                const faceList::subList& bFaces = mse.boundaryFaces();
                const label nBndFaces = bFaces.size();
                const label nInternalFaces = mesh_.nInternalFaces();

                labelHashSet rollbackBndFaces;

                forAllConstIter(labelHashSet, badCells, it)
                {
                    const label cellI = it.key();
                    for(label bfI=0; bfI<nBndFaces; ++bfI)
                        if( bFaceOwner[bfI] == cellI )
                            rollbackBndFaces.insert(bfI);
                }

                forAllConstIter(labelHashSet, badPyramidFaces, it)
                {
                    const label faceI = it.key();
                    if( faceI >= nInternalFaces )
                        rollbackBndFaces.insert(faceI - nInternalFaces);
                }

                forAllConstIter(labelHashSet, nonOrthoFaces, it)
                {
                    const label faceI = it.key();
                    if( faceI >= nInternalFaces )
                        rollbackBndFaces.insert(faceI - nInternalFaces);
                }

                Pout << "  Rollback candidates: "
                     << rollbackBndFaces.size()
                     << " boundary faces would be targeted" << endl;
            }
            else
            {
                Pout << "  All gates passed - no rollback needed" << endl;
            }
        }

        meshOptimizer mOpt(mesh_);
        mOpt.lockPoints(blPoints_);
        mOpt.untangleBoundaryLayer();
        // Post-refinement BL optimisation -- optional, off by default.
        // Enable with: postRefineBLOptimisation true; in boundaryLayers dict.
        bool postRefineBLOpt = false;
        if( meshDict_.isDict("boundaryLayers") )
        {
            const dictionary& bndL = meshDict_.subDict("boundaryLayers");
            if( bndL.found("postRefineBLOptimisation") )
                postRefineBLOpt = Switch(bndL.lookup("postRefineBLOptimisation"));
        }
        if( postRefineBLOpt )
        {
            Info << "Running post-refinement boundary-layer optimisation" << endl;
            mOpt.optimizeBoundaryLayer(modSurfacePtr_==NULL);
        }
        Info << "refBoundaryLayers: stored "
             << blPoints_.size()
             << " BL interior points" << endl;
    }
}

void cartesianMeshGenerator::optimiseFinalMesh()
{
    finalUntangleRejected_ = false;

    //- untangle the surface if needed
    bool enforceConstraints(false);
    if( meshDict_.found("enforceGeometryConstraints") )
    {
        enforceConstraints =
            readBool(meshDict_.lookup("enforceGeometryConstraints"));
    }

        bool lockAcuteCorners = false;
        if( meshDict_.isDict("boundaryLayers") )
        {
            const dictionary& bndL =
                meshDict_.subDict("boundaryLayers");
            if( bndL.found("lockAcuteCornerPoints") )
                lockAcuteCorners =
                    bool(Switch(bndL.lookup("lockAcuteCornerPoints")));
        }

    {
        meshSurfaceEngine mse(mesh_);
        meshSurfaceOptimizer surfOpt(mse, *octreePtr_);


        if( enforceConstraints )
            surfOpt.enforceConstraints();


        //- lock acute BL+BL+neutral corners: prevent optimizer
        //- from moving these points across patch boundaries
        //- controlled by meshDict: lockAcuteCornerPoints true/false
        if( lockAcuteCorners && blblAcuteCornerPoints_.size() )
        {
            Info << "Locking " << blblAcuteCornerPoints_.size()
                 << " acute BL+BL+neutral corner points in optimizer" << endl;
            surfOpt.setAcuteCornerPoints(blblAcuteCornerPoints_);
        }
        surfOpt.optimizeSurface();
    }

    //- final optimisation
    meshOptimizer optimizer(mesh_);
    if( enforceConstraints )
        optimizer.enforceConstraints();

    // Compute acute corner global point list once -- reused before both
    // optimizeMeshFV and untangleMeshFV since optimizeBoundaryLayer
    // calls removeUserConstraints() internally, wiping the first lock.
    labelLongList acuteGlobalPts;
    if( lockAcuteCorners && !blblAcuteCornerPoints_.empty() )
    {
        const meshSurfaceEngine mseForLock(mesh_);
        const labelList& bPoints = mseForLock.boundaryPoints();
        forAllConstIter(labelHashSet, blblAcuteCornerPoints_, it)
        {
            const label bpI = it.key();
            if( bpI >= 0 && bpI < label(bPoints.size()) )
                acuteGlobalPts.append(bPoints[bpI]);
        }
        if( acuteGlobalPts.size() > 0 )
        {
            optimizer.lockPoints(acuteGlobalPts);
            Info << "optimiseFinalMesh: locked "
                 << acuteGlobalPts.size()
                 << " acute corner points in volume optimizer" << endl;
        }
    }

    // Surface-constrained optimizer: gate by meshDict switch.
    bool constrainOptimizerBoundary = false;
    if( meshDict_.isDict("boundaryLayers") )
    {
        const dictionary& bndLC = meshDict_.subDict("boundaryLayers");
        if( bndLC.found("constrainOptimizerBoundaryMotion") )
            constrainOptimizerBoundary =
                Switch(bndLC.lookup("constrainOptimizerBoundaryMotion"));
    }

    if( constrainOptimizerBoundary && octreePtr_ )
    {
        Info << "Surface-constrained optimizer: building boundary mapping"
             << endl;
        meshSurfaceEngine mseConstraint(mesh_);
        meshSurfacePartitioner mPartConstraint(mseConstraint);
        labelLongList globalToBp(mesh_.points().size(), -1);
        const labelList& bPtsC = mseConstraint.boundaryPoints();
        forAll(bPtsC, bpI)
            globalToBp[bPtsC[bpI]] = bpI;

        // Build feature curve tangent vectors indexed by boundary point.
        // Zero vector = not a feature curve point.
        // Corner points get zero tangent and are handled by corner locking.
        const label nBp = bPtsC.size();
        vectorField featureTangents(nBp, vector::zero);
        {
            const edgeList& edges = mseConstraint.edges();
            const VRWGraph& bpEdges = mseConstraint.boundaryPointEdges();
            const labelHashSet& featEdges = mPartConstraint.featureEdges();
            const labelHashSet& edgePts = mPartConstraint.edgePoints();
            const pointFieldPMG& pts = mesh_.points();
            const labelList& bp = mseConstraint.bp();

            label nTangents = 0;

            forAllConstIter(labelHashSet, edgePts, it)
            {
                const label bpI = it.key();

                if( bpI < 0 || bpI >= nBp )
                    continue;

                label nbr0 = -1;
                label nbr1 = -1;

                forAllRow(bpEdges, bpI, eI)
                {
                    const label beI = bpEdges(bpI, eI);

                    if( !featEdges.found(beI) )
                        continue;

                    const edge& e = edges[beI];

                    // edges store mesh point indices -- convert to bp indices
                    const label ep0 = e.start();
                    const label ep1 = e.end();

                    if( ep0 < 0 || ep0 >= bp.size() || ep1 < 0 || ep1 >= bp.size() )
                        continue;

                    const label otherBp0 = bp[ep0];
                    const label otherBp1 = bp[ep1];

                    if( otherBp0 < 0 || otherBp1 < 0 )
                        continue;

                    label otherBp = -1;
                    if( otherBp0 == bpI )
                        otherBp = otherBp1;
                    else if( otherBp1 == bpI )
                        otherBp = otherBp0;

                    if( otherBp < 0 || otherBp >= nBp )
                        continue;

                    if( nbr0 == -1 )
                        nbr0 = otherBp;
                    else if( nbr1 == -1 && otherBp != nbr0 )
                        nbr1 = otherBp;
                }

                vector t = vector::zero;

                if( nbr0 != -1 && nbr1 != -1 )
                {
                    t = pts[bPtsC[nbr1]] - pts[bPtsC[nbr0]];
                }
                else if( nbr0 != -1 )
                {
                    t = pts[bPtsC[nbr0]] - pts[bPtsC[bpI]];
                }

                if( magSqr(t) > VSMALL )
                {
                    featureTangents[bpI] = t / mag(t);
                    ++nTangents;
                }
            }

            Info << "Surface-constrained optimizer: built "
                 << nTangents << " feature curve tangents" << endl;
        }

        Info << "Surface-constrained optimizer: active -- "
             << bPtsC.size() << " boundary points mapped, "
             << mPartConstraint.edgePoints().size() << " feature curve pts, "
             << mPartConstraint.corners().size() << " corner pts locked"
             << endl;

        optimizer.setSurfaceConstraint
        (
            octreePtr_,
            &mPartConstraint.pointPatches(),
            &globalToBp,
            &mPartConstraint.corners(),
            &featureTangents
        );
        optimizer.optimizeMeshFV();
        optimizer.optimizeLowQualityFaces();
        optimizer.setSurfaceConstraint(NULL, NULL, NULL, NULL, NULL);
    }
    else
    {
        optimizer.optimizeMeshFV();
        optimizer.optimizeLowQualityFaces();
    }

    // Post-optimizer surface re-projection: re-project drifted single-patch
    // boundary points while octree is still live. Tolerance-filtered and
    // validated -- same pattern as snapSurfaceBeforeBLRefinement.
    if( octreePtr_ )
    {
        scalar reprojTol = 1e-6;
        if( meshDict_.isDict("boundaryLayers") )
        {
            const dictionary& bndLR = meshDict_.subDict("boundaryLayers");
            if( bndLR.found("postBLSnapTolerance") )
                reprojTol = readScalar(bndLR.lookup("postBLSnapTolerance"));
        }
        const scalar reprojTolSq = sqr(reprojTol);
        Info << "Post-optimizer re-projection tolerance = "
             << reprojTol << " m" << endl;

        meshSurfaceEngine mseReproj(mesh_);
        meshSurfaceMapper mapperReproj(mseReproj, *octreePtr_);
        if( !blNoBlEdgePoints_.empty() )
        {
            mapperReproj.setProtectedPoints(blNoBlEdgePoints_);
            mapperReproj.setProtectedPointPatches(blNoBlPointPatch_);
        }
        const labelList& bPtsR = mseReproj.boundaryPoints();
        const meshSurfacePartitioner mPartR(mseReproj);
        const VRWGraph& pPatchesR = mPartR.pointPatches();

        // Per-patch drift diagnostics
        Map<label> drift1e6ByPatch;
        Map<label> drift1e5ByPatch;
        Map<label> drift1e4ByPatch;
        Map<scalar> maxDriftByPatch;

        // Build protected point set: skip points attached to existing
        // bad cells or bad pyramid faces -- those zones are fragile and
        // any move near them risks creating new negVol.
        labelHashSet negCellsBefore, pyrFacesBefore;
        polyMeshGenChecks::checkCellVolumes(mesh_, false, &negCellsBefore);
        polyMeshGenChecks::checkFacePyramids
            (mesh_, false, -SMALL, &pyrFacesBefore);

        labelHashSet protectedPts;
        {
            const faceListPMG& allFaces = mesh_.faces();
            const cellListPMG& allCells = mesh_.cells();
            const VRWGraph& ptPts = mesh_.addressingData().pointPoints();

            // Points on bad pyramid faces + two-ring neighbors
            labelHashSet pyrPts;
            forAllConstIter(labelHashSet, pyrFacesBefore, it)
            {
                const face& f = allFaces[it.key()];
                forAll(f, fpI) pyrPts.insert(f[fpI]);
            }
            // Ring 1
            labelHashSet pyrPtsRing1;
            forAllConstIter(labelHashSet, pyrPts, it)
            {
                protectedPts.insert(it.key());
                forAllRow(ptPts, it.key(), nI)
                {
                    protectedPts.insert(ptPts(it.key(), nI));
                    pyrPtsRing1.insert(ptPts(it.key(), nI));
                }
            }
            // Ring 2
            forAllConstIter(labelHashSet, pyrPtsRing1, it)
            {
                protectedPts.insert(it.key());
                forAllRow(ptPts, it.key(), nI)
                    protectedPts.insert(ptPts(it.key(), nI));
            }

            // Points on negative volume cells + one-ring neighbors
            labelHashSet negPts;
            forAllConstIter(labelHashSet, negCellsBefore, it)
            {
                const cell& c = allCells[it.key()];
                forAll(c, cfI)
                {
                    const face& f = allFaces[c[cfI]];
                    forAll(f, fpI) negPts.insert(f[fpI]);
                }
            }
            forAllConstIter(labelHashSet, negPts, it)
            {
                protectedPts.insert(it.key());
                forAllRow(ptPts, it.key(), nI)
                    protectedPts.insert(ptPts(it.key(), nI));
            }
        }
        Info << "Post-optimizer re-projection: protecting "
             << protectedPts.size()
             << " points near existing bad cells/faces (incl. 1-ring)" << endl;

        labelLongList reprojPoints;
        forAll(bPtsR, bpI)
        {
            const label meshPtI = bPtsR[bpI];
            if( meshPtI < 0 || meshPtI >= label(mesh_.points().size()) )
                continue;
            if( pPatchesR.sizeOfRow(bpI) != 1 )
                continue;
            // Skip points attached to fragile cells/faces
            if( protectedPts.found(meshPtI) )
                continue;
            const label patchI = pPatchesR(bpI, 0);
            point testPt;
            scalar testDsq;
            label testNt;
            octreePtr_->findNearestSurfacePointInRegion
            (
                testPt,
                testDsq,
                testNt,
                patchI,
                mesh_.points()[meshPtI]
            );
            const scalar drift = Foam::sqrt(testDsq);
            if( drift > 1e-6 )
            {
                if( !drift1e6ByPatch.found(patchI) )
                    drift1e6ByPatch.set(patchI, 0);
                ++drift1e6ByPatch[patchI];
            }
            if( drift > 1e-5 )
            {
                if( !drift1e5ByPatch.found(patchI) )
                    drift1e5ByPatch.set(patchI, 0);
                ++drift1e5ByPatch[patchI];
            }
            if( drift > 1e-4 )
            {
                if( !drift1e4ByPatch.found(patchI) )
                    drift1e4ByPatch.set(patchI, 0);
                ++drift1e4ByPatch[patchI];
            }
            if( !maxDriftByPatch.found(patchI) ||
                drift > maxDriftByPatch[patchI] )
                maxDriftByPatch.set(patchI, drift);
            if( testDsq < reprojTolSq )
                continue;
            reprojPoints.append(bpI);
        }

        // Print per-patch drift summary
        const PtrList<boundaryPatch>& bPatches = mesh_.boundaries();
        Info << "Post-optimizer drift by patch:" << endl;
        forAll(bPatches, pI)
        {
            const label n6 = drift1e6ByPatch.found(pI) ?
                drift1e6ByPatch[pI] : 0;
            const label n5 = drift1e5ByPatch.found(pI) ?
                drift1e5ByPatch[pI] : 0;
            const label n4 = drift1e4ByPatch.found(pI) ?
                drift1e4ByPatch[pI] : 0;
            const scalar mxD = maxDriftByPatch.found(pI) ?
                maxDriftByPatch[pI] : scalar(0);
            if( n6 > 0 )
                Info << "  patch " << bPatches[pI].patchName()
                     << ": >1e-6=" << n6
                     << " >1e-5=" << n5
                     << " >1e-4=" << n4
                     << " max=" << mxD
                     << endl;
        }

        Info << "Post-optimizer re-projection candidates: "
             << reprojPoints.size() << endl;

        // Limited-displacement re-projection: move each drifted point
        // a fraction of the way to the STL per pass. Cells gradually
        // reshape toward the correct position rather than being inverted
        // in one full snap. Read step fraction from meshDict.
        scalar reprojStepFraction = 0.02;
        label reprojMaxPasses = 5;
        if( meshDict_.isDict("boundaryLayers") )
        {
            const dictionary& bndLR2 =
                meshDict_.subDict("boundaryLayers");
            if( bndLR2.found("postOptimizerReprojStepFraction") )
                reprojStepFraction =
                    readScalar(bndLR2.lookup("postOptimizerReprojStepFraction"));
            if( bndLR2.found("postOptimizerReprojPasses") )
                reprojMaxPasses =
                    readLabel(bndLR2.lookup("postOptimizerReprojPasses"));
        }

        if( reprojPoints.size() > 0 )
        {
            meshSurfaceEngineModifier surfModR(mseReproj);
            label totalAccepted = 0;
            label totalRolledBack = 0;

            for( label passI = 0; passI < reprojMaxPasses; ++passI )
            {
                // Snapshot for this pass
                const pointField passPtsBefore(mesh_.points());
                labelHashSet negPassBefore, pyrPassBefore;
                polyMeshGenChecks::checkCellVolumes
                    (mesh_, false, &negPassBefore);
                polyMeshGenChecks::checkFacePyramids
                    (mesh_, false, -SMALL, &pyrPassBefore);

                // Apply limited displacement toward STL
                label nMoved = 0;
                forAll(reprojPoints, rpI)
                {
                    const label bpI = reprojPoints[rpI];
                    const label meshPtI = bPtsR[bpI];
                    if( meshPtI < 0 ||
                        meshPtI >= label(mesh_.points().size()) )
                        continue;
                    const label patchI = pPatchesR(bpI, 0);
                    point snapPt;
                    scalar snapDsq;
                    label snapNt;
                    octreePtr_->findNearestSurfacePointInRegion
                    (
                        snapPt,
                        snapDsq,
                        snapNt,
                        patchI,
                        mesh_.points()[meshPtI]
                    );
                    if( snapDsq < reprojTolSq )
                        continue;
                    const vector disp =
                        snapPt - mesh_.points()[meshPtI];
                    const point limitedPt =
                        mesh_.points()[meshPtI]
                        + reprojStepFraction * disp;
                    surfModR.moveBoundaryVertexNoUpdate(bpI, limitedPt);
                    ++nMoved;
                }
                surfModR.updateGeometry(reprojPoints);
                mesh_.clearAddressingData();

                labelHashSet negPassAfter, pyrPassAfter;
                polyMeshGenChecks::checkCellVolumes
                    (mesh_, false, &negPassAfter);
                polyMeshGenChecks::checkFacePyramids
                    (mesh_, false, -SMALL, &pyrPassAfter);

                if( negPassAfter.size() > negPassBefore.size() )
                {
                    Info << "Post-optimizer re-projection pass "
                         << passI << " rejected: negVol "
                         << negPassBefore.size() << "->"
                         << negPassAfter.size()
                         << " badPyramids "
                         << pyrPassBefore.size() << "->"
                         << pyrPassAfter.size()
                         << " -- rolling back" << endl;
                    pointField& pts = mesh_.points();
                    pts = passPtsBefore;
                    mesh_.clearAddressingData();
                    ++totalRolledBack;
                    break;
                }
                else
                {
                    Info << "Post-optimizer re-projection pass "
                         << passI << " accepted: moved " << nMoved
                         << " negVol "
                         << negPassBefore.size() << "->"
                         << negPassAfter.size()
                         << " badPyramids "
                         << pyrPassBefore.size() << "->"
                         << pyrPassAfter.size();
                    if( pyrPassAfter.size() > pyrPassBefore.size() )
                        Info << " WARNING: badPyramids increased"
                             << " during limited correction";
                    Info << endl;
                    ++totalAccepted;
                }
            }
            Info << "Post-optimizer re-projection: "
                 << totalAccepted << " passes accepted, "
                 << totalRolledBack << " rolled back" << endl;
        }
    }
    deleteDemandDrivenData(octreePtr_);

    optimizer.optimizeBoundaryLayer(modSurfacePtr_==NULL);

    // Second low-quality face pass after BL refinement -- targets skew
    // introduced by boundary layer cells that weren't present pre-BL.
    optimizer.optimizeLowQualityFaces();

    // Post-BL validity audit: find incorrectly oriented faces and attempt
    // conservative face-flip repair with full accept/reject validation.
    {
        labelHashSet badFaces;
        polyMeshGenChecks::checkFacePyramids(mesh_, false, -SMALL, &badFaces);

        if( badFaces.size() > 0 )
        {
            Info << "Post-BL audit: " << badFaces.size()
                 << " incorrectly oriented faces -- attempting validated face-flip repair"
                 << endl;

            labelHashSet negBefore;
            polyMeshGenChecks::checkCellVolumes(mesh_, false, &negBefore);

            labelHashSet openBefore;
            polyMeshGenChecks::checkClosedCells(mesh_, false, 0.5, &openBefore);

            // Stage 1: targeted Laplacian smoothing on bad face neighbourhood
            {
                scalarField skewS1Before;
                polyMeshGenChecks::checkFaceSkewness(mesh_, skewS1Before);
                const scalar maxSkewS1Before =
                    skewS1Before.size() > 0 ? max(skewS1Before) : scalar(0.0);

                const pointField pointsBefore(mesh_.points());
                optimizer.optimizeLowQualityFaces(3);

                labelHashSet badStage1;
                polyMeshGenChecks::checkFacePyramids(mesh_, false, -SMALL, &badStage1);
                labelHashSet negStage1;
                polyMeshGenChecks::checkCellVolumes(mesh_, false, &negStage1);
                labelHashSet openStage1;
                polyMeshGenChecks::checkClosedCells(mesh_, false, 0.5, &openStage1);
                scalarField skewS1After;
                polyMeshGenChecks::checkFaceSkewness(mesh_, skewS1After);
                const scalar maxSkewS1After =
                    skewS1After.size() > 0 ? max(skewS1After) : scalar(0.0);

                const bool stage1OK =
                    badStage1.size() <= badFaces.size()
                 && negStage1.size() <= negBefore.size()
                 && openStage1.size() <= openBefore.size()
                 && maxSkewS1After <= scalar(20.0);

                if( stage1OK )
                {
                    Info << "Post-BL stage1 Laplacian accepted: bad "
                         << badFaces.size() << "->" << badStage1.size()
                         << " skew " << maxSkewS1Before << "->" << maxSkewS1After
                         << endl;
                    badFaces   = badStage1;
                    negBefore  = negStage1;
                    openBefore = openStage1;
                }
                else
                {
                    Info << "Post-BL stage1 rejected: bad "
                         << badFaces.size() << "->" << badStage1.size()
                         << " skew " << maxSkewS1Before << "->" << maxSkewS1After
                         << " -- rolling back" << endl;
                    polyMeshGenModifier meshModifier2(mesh_);
                    pointFieldPMG& pts = meshModifier2.pointsAccess();
                    pts = pointsBefore;
                    mesh_.clearAddressingData();
                }
            }

            // Stage 2: face flip on remaining bad internal faces
            polyMeshGenModifier meshMod(mesh_);
            faceListPMG& faces = meshMod.facesAccess();

            label nFixed(0);
            label nRejected(0);
            label nInvalid(0);

            label currentBad  = badFaces.size();
            label currentNeg  = negBefore.size();
            label currentOpen = openBefore.size();

            const labelList& owner = mesh_.owner();
            const labelList& neighbour = mesh_.neighbour();

            // Build cell -> boundary patch contact map. Periodic transition
            // bad faces are usually internal triangular faces adjacent to
            // cells touching neutral/periodic patches. These are not ordinary
            // orientation errors; flipping them tends to open the adjacent
            // transition cells.
            bool periodicTransitionProtection(true);
            bool periodicTransitionSkipFlip(true);

            wordHashSet periodicTransitionPatches;
            periodicTransitionPatches.insert("periodic_1");
            periodicTransitionPatches.insert("periodic_2");

            if( meshDict_.isDict("boundaryLayers") )
            {
                const dictionary& bndL = meshDict_.subDict("boundaryLayers");

                if( bndL.found("periodicTransitionProtection") )
                    periodicTransitionProtection =
                        bool(Switch(bndL.lookup("periodicTransitionProtection")));

                if( bndL.found("periodicTransitionSkipFlip") )
                    periodicTransitionSkipFlip =
                        bool(Switch(bndL.lookup("periodicTransitionSkipFlip")));

                if( bndL.found("periodicTransitionPatches") )
                {
                    periodicTransitionPatches.clear();

                    const wordList pNames(bndL.lookup("periodicTransitionPatches"));
                    forAll(pNames, pI)
                        periodicTransitionPatches.insert(pNames[pI]);
                }
            }

            const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
            boolList cellTouchesPeriodic(mesh_.cells().size(), false);

            if( periodicTransitionProtection )
            {
                forAll(boundaries, patchI)
                {
                    const word& pName = boundaries[patchI].patchName();

                    if( !periodicTransitionPatches.found(pName) )
                        continue;

                    const label start = boundaries[patchI].patchStart();
                    const label size  = boundaries[patchI].patchSize();

                    for(label pfI=0; pfI<size; ++pfI)
                    {
                        const label faceJ = start + pfI;
                        if( faceJ < 0 || faceJ >= label(owner.size()) ) continue;

                        const label c = owner[faceJ];
                        if( c >= 0 && c < label(cellTouchesPeriodic.size()) )
                            cellTouchesPeriodic[c] = true;
                    }
                }
            }

            Info << "Post-BL periodic transition protection: "
                 << (periodicTransitionProtection ? "enabled" : "disabled")
                 << ", skipFlip=" << (periodicTransitionSkipFlip ? "true" : "false")
                 << ", nPatches=" << periodicTransitionPatches.size()
                 << ", patches=(";

            forAllConstIter(wordHashSet, periodicTransitionPatches, pIt)
                Info << " " << pIt.key();

            Info << " )" << endl;

            label nPeriodicTransitionSkipped(0);

            forAllConstIter(labelHashSet, badFaces, it)
            {
                const label faceI = it.key();

                if( faceI < 0 || faceI >= label(faces.size()) )
                { ++nInvalid; continue; }

                // Only attempt flips on internal faces. Boundary faces have
                // only an owner-side pyramid; reversing them usually makes
                // patch orientation worse instead of repairing a cell.
                if( faceI >= label(neighbour.size()) || neighbour[faceI] < 0 )
                { ++nRejected; continue; }

                const label ownCell = owner[faceI];
                const label neiCell = neighbour[faceI];

                const bool periodicTransitionFace =
                    faces[faceI].size() <= 4
                 && (
                        (ownCell >= 0 && ownCell < label(cellTouchesPeriodic.size())
                      && cellTouchesPeriodic[ownCell])
                     || (neiCell >= 0 && neiCell < label(cellTouchesPeriodic.size())
                      && cellTouchesPeriodic[neiCell])
                    );

                if
                (
                    periodicTransitionProtection
                 && periodicTransitionSkipFlip
                 && periodicTransitionFace
                )
                {
                    ++nPeriodicTransitionSkipped;
                    ++nRejected;

                    if( nPeriodicTransitionSkipped <= 10 )
                    {
                        Info << "Post-BL audit: periodic transition bad face faceI="
                             << faceI
                             << " owner=" << ownCell
                             << " neighbour=" << neiCell
                             << " nPts=" << faces[faceI].size()
                             << " -- skipping flip" << endl;
                    }

                    continue;
                }

                faces[faceI] = faces[faceI].reverseFace();
                mesh_.clearAddressingData();

                labelHashSet badAfter;
                polyMeshGenChecks::checkFacePyramids(mesh_, false, -SMALL, &badAfter);

                labelHashSet negAfter;
                polyMeshGenChecks::checkCellVolumes(mesh_, false, &negAfter);

                labelHashSet openAfter;
                polyMeshGenChecks::checkClosedCells(mesh_, false, 0.5, &openAfter);

                const bool acceptFlip =
                    badAfter.size() < currentBad
                 && negAfter.size() <= currentNeg
                 && openAfter.size() <= currentOpen;

                if( acceptFlip )
                {
                    ++nFixed;
                    currentBad  = badAfter.size();
                    currentNeg  = negAfter.size();
                    currentOpen = openAfter.size();
                    Info << "Post-BL audit: accepted flip faceI=" << faceI
                         << " bad=" << currentBad
                         << " neg=" << currentNeg
                         << " open=" << currentOpen << endl;
                }
                else
                {
                    if( nRejected < 5 )
                    {
                        const labelList& own = mesh_.owner();
                        const labelList& nei = mesh_.neighbour();
                        Info << "Post-BL audit: rejected flip faceI=" << faceI
                             << " bad " << currentBad << "->" << badAfter.size()
                             << " neg " << currentNeg << "->" << negAfter.size()
                             << " open " << currentOpen << "->" << openAfter.size()
                             << " owner=" << own[faceI]
                             << " neighbour="
                             << (faceI < label(nei.size()) ? nei[faceI] : -1)
                             << endl;
                    }
                    faces[faceI] = faces[faceI].reverseFace();
                    mesh_.clearAddressingData();
                    ++nRejected;
                }
            }

            Info << "Post-BL audit: fixed=" << nFixed
                 << " rejected=" << nRejected
                 << " invalid=" << nInvalid
                 << " periodicTransitionSkipped=" << nPeriodicTransitionSkipped
                 << " remainingBad=" << currentBad
                 << " negVol=" << currentNeg
                 << " openCells=" << currentOpen << endl;
        }
        else
        {
            Info << "Post-BL audit: all face pyramids OK" << endl;
        }
    }

    // Final untangle intentionally runs after optimizeBoundaryLayer has
    // cleared user constraints via removeUserConstraints(). Acute corner
    // locks are useful during broad optimization phases but hard-locking
    // them during final untangle prevents closure of cells adjacent to
    // junctions. The natural constraint reset is the correct policy here.
    //
    // Protect with accept/reject rollback -- untangle can occasionally
    // make junction cells worse (skew 233, neg vol) on bad OMP paths.
    {
        labelHashSet badBefore;
        polyMeshGenChecks::checkFacePyramids(mesh_, false, -SMALL, &badBefore);

        labelHashSet negBefore;
        polyMeshGenChecks::checkCellVolumes(mesh_, false, &negBefore);

        labelHashSet openBefore;
        polyMeshGenChecks::checkClosedCells(mesh_, false, 0.5, &openBefore);

        scalarField skewBefore;
        polyMeshGenChecks::checkFaceSkewness(mesh_, skewBefore);
        const scalar maxSkewBefore =
            skewBefore.size() > 0 ? max(skewBefore) : scalar(0.0);

        const pointField pointsBefore(mesh_.points());

        // Always attempt untangleMeshFV -- it exists to fix inverted cells.
        // optimizeMeshFV can create negVol cells; skipping the untangler
        // when negVol>0 is self-defeating. The rollback below rejects the
        // result if untangle makes things worse (negVol increases).
        {
            optimizer.untangleMeshFV();

            labelHashSet badAfter;
            polyMeshGenChecks::checkFacePyramids(mesh_, false, -SMALL, &badAfter);
            labelHashSet negAfter;
            polyMeshGenChecks::checkCellVolumes(mesh_, false, &negAfter);
            labelHashSet openAfter;
            polyMeshGenChecks::checkClosedCells(mesh_, false, 0.5, &openAfter);
            scalarField skewAfter;
            polyMeshGenChecks::checkFaceSkewness(mesh_, skewAfter);
            const scalar maxSkewAfter =
                skewAfter.size() > 0 ? max(skewAfter) : scalar(0.0);
            const bool skewOK =
                maxSkewAfter <= scalar(20.0)
             && maxSkewAfter <= scalar(2.0) *
                    Foam::max(maxSkewBefore, scalar(1.0));
            const bool untangleOK =
                badAfter.size() <= badBefore.size()
             && negAfter.size() <= negBefore.size()
             && openAfter.size() <= openBefore.size()
             && skewOK;
            if( !untangleOK )
            {
                Info << "Final untangle rejected: badFaces "
                     << badBefore.size() << "->" << badAfter.size()
                     << " negVol " << negBefore.size() << "->" << negAfter.size()
                     << " openCells "
                     << openBefore.size() << "->" << openAfter.size()
                     << " -- rolling back" << endl;
                polyMeshGenModifier meshModifier(mesh_);
                pointFieldPMG& pts = meshModifier.pointsAccess();
                pts = pointsBefore;
                mesh_.clearAddressingData();
                finalUntangleRejected_ = true;
            }
            else
            {
                Info << "Final untangle accepted: badFaces "
                     << badBefore.size() << "->" << badAfter.size()
                     << " negVol " << negBefore.size() << "->" << negAfter.size()
                     << " openCells "
                     << openBefore.size() << "->" << openAfter.size()
                     << endl;
            }
        }
    }

    mesh_.clearAddressingData();

    if( modSurfacePtr_ )
    {
        polyMeshGenGeometryModification meshMod(mesh_, meshDict_);

        //- revert the mesh into the original space
        meshMod.revertGeometryModification();

        //- delete modified surface mesh
        deleteDemandDrivenData(modSurfacePtr_);
    }
}

void cartesianMeshGenerator::projectSurfaceAfterBackScaling()
{
    if( !meshDict_.found("anisotropicSources") )
        return;

    deleteDemandDrivenData(octreePtr_);
    octreePtr_ = new meshOctree(*surfacePtr_);

    meshOctreeCreator
    (
        *octreePtr_,
        meshDict_
    ).createOctreeWithRefinedBoundary(20, 30);

    //- calculate mesh surface
    meshSurfaceEngine mse(mesh_);

    //- pre-map mesh surface
    meshSurfaceMapper mapper(mse, *octreePtr_);

    //- map mesh surface on the geometry surface
    mapper.mapVerticesOntoSurface();

    optimiseFinalMesh();
}

void cartesianMeshGenerator::replaceBoundaries()
{
    renameBoundaryPatches rbp(mesh_, meshDict_);
}

void cartesianMeshGenerator::renumberMesh()
{
    polyMeshGenModifier(mesh_).renumberMesh();
}

void cartesianMeshGenerator::snapSurfaceBeforeBLRefinement()
{
    if( !meshDict_.isDict("boundaryLayers") )
        return;
    const dictionary& bndL = meshDict_.subDict("boundaryLayers");
    if( !bndL.found("postBLSnap") )
        return;
    if( !Switch(bndL.lookup("postBLSnap")) )
        return;

    Info << "Pre-BL snap: re-projecting surface after volume optimisation" << endl;

    // Rebuild octree -- deleted by optimiseFinalMesh
    meshOctree* snapOctreePtr = new meshOctree(*surfacePtr_);
    meshOctreeCreator
    (
        *snapOctreePtr,
        meshDict_
    ).createOctreeWithRefinedBoundary(20, 30);

    meshSurfaceEngine mse(mesh_);
    meshSurfaceMapper mapper(mse, *snapOctreePtr);

    // Protect BL/no-BL and BL/neutral interface points
    if( !blNoBlEdgePoints_.empty() )
    {
        mapper.setProtectedPoints(blNoBlEdgePoints_);
        mapper.setProtectedPointPatches(blNoBlPointPatch_);
    }
    if( !blNeutralEdgePoints_.empty() )
    {
        mapper.setBLNeutralPoints(blNeutralEdgePoints_);
        mapper.setBLNeutralPointPatches(blNeutralPointPatch_);
    }

    // Only snap true single-patch points -- generic nearest-surface
    // projection is safe only for pure wall-face interior points.
    // Multi-patch edges/corners/junctions need feature-aware mapping.
    const labelList& bPoints = mse.boundaryPoints();
    const meshSurfacePartitioner mPart(mse);
    const VRWGraph& pPatches = mPart.pointPatches();

    // Read snap tolerance from meshDict -- default 1e-6 m
    scalar snapTol = 1e-6;
    if( bndL.found("postBLSnapTolerance") )
        snapTol = readScalar(bndL.lookup("postBLSnapTolerance"));
    const scalar snapTolSq = sqr(snapTol);
    Info << "Pre-BL snap displacement tolerance = " << snapTol << " m" << endl;

    labelLongList snapPoints;
    label nSinglePatch = 0;
    label nAlreadyOnSurface = 0;
    label nSnapCandidates = 0;
    label nDrift1e9 = 0, nDrift1e8 = 0, nDrift1e7 = 0;
    label nDrift1e6 = 0, nDrift1e5 = 0, nDrift1e4 = 0;
    scalar maxDrift = 0.0;
    forAll(bPoints, bpI)
    {
        const label meshPtI = bPoints[bpI];
        if( meshPtI < 0 || meshPtI >= label(mesh_.points().size()) )
            continue;
        if( pPatches.sizeOfRow(bpI) != 1 )
            continue;
        ++nSinglePatch;
        point testPt;
        scalar testDsq;
        label testNt;
        snapOctreePtr->findNearestSurfacePointInRegion
        (
            testPt,
            testDsq,
            testNt,
            pPatches(bpI, 0),
            mesh_.points()[meshPtI]
        );
        const scalar drift = Foam::sqrt(testDsq);
        maxDrift = Foam::max(maxDrift, drift);
        if( drift > 1e-9 ) ++nDrift1e9;
        if( drift > 1e-8 ) ++nDrift1e8;
        if( drift > 1e-7 ) ++nDrift1e7;
        if( drift > 1e-6 ) ++nDrift1e6;
        if( drift > 1e-5 ) ++nDrift1e5;
        if( drift > 1e-4 ) ++nDrift1e4;
        if( testDsq < snapTolSq )
        {
            ++nAlreadyOnSurface;
            continue;
        }
        snapPoints.append(bpI);
        ++nSnapCandidates;
    }

    Info << "Pre-BL snap drift histogram: singlePatch=" << nSinglePatch
         << " >1e-9=" << nDrift1e9
         << " >1e-8=" << nDrift1e8
         << " >1e-7=" << nDrift1e7
         << " >1e-6=" << nDrift1e6
         << " >1e-5=" << nDrift1e5
         << " >1e-4=" << nDrift1e4
         << " max=" << maxDrift
         << endl;
    Info << "Pre-BL snap candidates: singlePatch=" << nSinglePatch
         << " alreadyOnSurface=" << nAlreadyOnSurface
         << " snapCandidates=" << nSnapCandidates
         << endl;

    if( snapPoints.size() == 0 )
    {
        deleteDemandDrivenData(snapOctreePtr);
        return;
    }

    // Snapshot for rollback
    const pointField pointsBeforeSnap(mesh_.points());
    labelHashSet negBefore, pyrBefore;
    polyMeshGenChecks::checkCellVolumes(mesh_, false, &negBefore);
    polyMeshGenChecks::checkFacePyramids(mesh_, false, -SMALL, &pyrBefore);

    mapper.mapVerticesOntoSurface(snapPoints);
    mesh_.clearAddressingData();

    labelHashSet negAfter, pyrAfter;
    polyMeshGenChecks::checkCellVolumes(mesh_, false, &negAfter);
    polyMeshGenChecks::checkFacePyramids(mesh_, false, -SMALL, &pyrAfter);

    const label pyrIncrease =
        label(pyrAfter.size()) - label(pyrBefore.size());
    label allowedPyrMin = 25;
    scalar allowedPyrFrac = 0.001;
    if( bndL.found("postBLSnapAllowedPyramidIncrease") )
        allowedPyrMin =
            readLabel(bndL.lookup("postBLSnapAllowedPyramidIncrease"));
    if( bndL.found("postBLSnapAllowedPyramidIncreaseFraction") )
        allowedPyrFrac =
            readScalar(bndL.lookup("postBLSnapAllowedPyramidIncreaseFraction"));
    const label allowedPyrIncrease =
        Foam::max
        (
            allowedPyrMin,
            label(allowedPyrFrac * Foam::max(label(1), label(pyrBefore.size())))
        );
    if( negAfter.size() > negBefore.size()
     || pyrIncrease > allowedPyrIncrease )
    {
        Info << "Pre-BL snap rejected: negVol "
             << negBefore.size() << "->" << negAfter.size()
             << " badPyramids " << pyrBefore.size() << "->" << pyrAfter.size()
             << " allowedPyrIncrease " << allowedPyrIncrease
             << " -- rolling back" << endl;
        pointField& pts = mesh_.points();
        pts = pointsBeforeSnap;
        mesh_.clearAddressingData();
    }
    else
    {
        Info << "Pre-BL snap accepted: negVol "
             << negBefore.size() << "->" << negAfter.size()
             << " badPyramids " << pyrBefore.size() << "->" << pyrAfter.size()
             << " allowedPyrIncrease " << allowedPyrIncrease
             << endl;
    }

    deleteDemandDrivenData(snapOctreePtr);
}

void cartesianMeshGenerator::generateMesh()
{
    try
    {
        if( controller_.runCurrentStep("templateGeneration") )
        {
            createCartesianMesh();
        }

        if( controller_.runCurrentStep("surfaceTopology") )
        {
            surfacePreparation();
        }

        if( controller_.runCurrentStep("patchAssignment") )
        {
            // Patch assignment moved before surface projection so that
            // mapVerticesOntoSurface has valid patch identity available.
            // edgeExtractor uses only mesh topology + octree -- no
            // dependency on projected surface positions.
            extractPatches();
        }

        if( controller_.runCurrentStep("surfaceProjection") )
        {
            mapMeshToSurface();
            // Re-run patch assignment after projection to correct any
            // misassignments that occurred on the unprojected hex mesh.
            extractPatches();
        }

        if( controller_.runCurrentStep("edgeExtraction") )
        {
            // Detect BL/no-BL transition edge points before any snapping
            // so all mapper instances in this block can exclude them
            // from generic nearest-surface projection.
            // Uses meshDict nLayersForPatch only -- no mesh modification.
            {
                boundaryLayers blDetect(mesh_, meshDict_);
                blDetect.detectBLNoBlTransitionEdges();
                blNoBlEdgePoints_ = blDetect.blNoBlEdgePoints();
                blNoBlPointPatch_ = blDetect.blNoBlPointPatch();
                blNeutralEdgePoints_ = blDetect.blNeutralEdgePoints();
                blNeutralPointPatch_ = blDetect.blNeutralPointPatch();
                Info << "Edge extraction: BL/no-BL interface points protected: "
                     << blNoBlEdgePoints_.size() << endl;
                Info << "Edge extraction: BL/neutral interface points detected: "
                     << blNeutralEdgePoints_.size() << endl;

                // Write BL/neutral edge points to VTK for spatial verification
                // Enable with: writeDiagnosticVTK true; in meshDict
                bool writeDiagVTK = false;
                if( meshDict_.found("writeDiagnosticVTK") )
                    writeDiagVTK = Switch(meshDict_.lookup("writeDiagnosticVTK"));
                if( writeDiagVTK && blNeutralEdgePoints_.size() > 0 )
                {
                    const meshSurfaceEngine mseVtk(mesh_);
                    const labelList& bPts = mseVtk.boundaryPoints();
                    const pointFieldPMG& allPts = mesh_.points();
                    const label nNeutral = blNeutralEdgePoints_.size();
                    OFstream osVtk("blNeutralEdgePoints_predetect.vtk");
                    osVtk << "# vtk DataFile Version 2.0\n";
                    osVtk << "blNeutralEdgePoints\n";
                    osVtk << "ASCII\n";
                    osVtk << "DATASET POLYDATA\n";
                    osVtk << "POINTS " << nNeutral << " float\n";
                    forAll(bPts, bpI)
                    {
                        if( !blNeutralEdgePoints_.found(bpI) ) continue;
                        const point& p = allPts[bPts[bpI]];
                        osVtk << p.x() << " " << p.y() << " " << p.z() << "\n";
                    }
                    osVtk << "VERTICES " << nNeutral << " " << 2*nNeutral << "\n";
                    for(label k=0; k<nNeutral; ++k)
                        osVtk << "1 " << k << "\n";
                    Info << "Wrote " << nNeutral
                         << " blNeutralEdgePoints to blNeutralEdgePoints_predetect.vtk" << endl;
                }
            }

            mapEdgesAndCorners();

            optimiseMeshSurface();


            // Step 1: snap corner points first (damped relaxation)
            // Single pass -- full BL/no-BL and BL/neutral protection.
            {
                scalar cornerSnapRelax = 0.25;
                if( meshDict_.isDict("boundaryLayers") )
                {
                    const dictionary& bndL =
                        meshDict_.subDict("boundaryLayers");
                    if( bndL.found("cornerSnapRelaxation") )
                        cornerSnapRelax = readScalar
                        (
                            bndL.lookup("cornerSnapRelaxation")
                        );
                }
                meshSurfaceEngine mse(mesh_);
                meshSurfacePartitioner mPart(mse);
                meshSurfaceMapper mapper(mse, *octreePtr_);
                mapper.setCornerSnapRelaxation(cornerSnapRelax);
                if( !blNoBlEdgePoints_.empty() )
                {
                    mapper.setProtectedPoints(blNoBlEdgePoints_);
                    mapper.setProtectedPointPatches(blNoBlPointPatch_);
                }
                if( !blNeutralEdgePoints_.empty() )
                {
                    mapper.setBLNeutralPoints(blNeutralEdgePoints_);
                    mapper.setBLNeutralPointPatches(blNeutralPointPatch_);
                }
                const labelHashSet& corners = mPart.corners();
                labelLongList cornerPts;
                forAllConstIter(labelHashSet, corners, it)
                    cornerPts.append(it.key());
                Info << "Ordered snap: snapping "
                     << cornerPts.size()
                     << " corner points (relax=" << cornerSnapRelax
                     << ")" << endl;
                mapper.mapCorners(cornerPts);
            }

            // Step 2: snap non-corner edge points after corners
            {
                meshSurfaceEngine mse(mesh_);
                meshSurfacePartitioner mPart(mse);
                meshSurfaceMapper mapper(mse, *octreePtr_);
                if( !blNoBlEdgePoints_.empty() )
                {
                    mapper.setProtectedPoints(blNoBlEdgePoints_);
                    mapper.setProtectedPointPatches(blNoBlPointPatch_);
                }
                if( !blNeutralEdgePoints_.empty() )
                {
                    mapper.setBLNeutralPoints(blNeutralEdgePoints_);
                    mapper.setBLNeutralPointPatches(blNeutralPointPatch_);
                }
                const labelHashSet& edgePoints = mPart.edgePoints();
                const labelHashSet& corners = mPart.corners();
                labelLongList edgePts;
                forAllConstIter(labelHashSet, edgePoints, it)
                {
                    const label bpI = it.key();
                    if( !corners.found(bpI) )
                        edgePts.append(bpI);
                }
                Info << "Ordered snap: snapping "
                     << edgePts.size()
                     << " non-corner edge points" << endl;
                mapper.mapEdgeNodes(edgePts);
            }
        }

        // Stage-gated bad-cell lineage writer.
        // Writes negVolCellCentres_<stage>.csv at 5 pipeline checkpoints.
        // Comparing stages tells us which step creates the bad cells.
        bool writeLineageDiagnostics = false;
        if( meshDict_.found("writeLineageDiagnostics") )
            writeLineageDiagnostics =
                Switch(meshDict_.lookup("writeLineageDiagnostics"));

        auto writeLineageCSV = [&](const std::string& stageName)
        {
            if( !writeLineageDiagnostics ) return;

            wordList enabledStages;
            if( meshDict_.found("writeLineageStages") )
                enabledStages = wordList(meshDict_.lookup("writeLineageStages"));
            else
                enabledStages = wordList(2);

            if( enabledStages.size() == 2 && enabledStages[0].empty() )
            {
                enabledStages[0] = word("postRefBL");
                enabledStages[1] = word("final");
            }

            bool stageEnabled = false;
            forAll(enabledStages, si)
            {
                if( enabledStages[si] == word(stageName.c_str()) )
                {
                    stageEnabled = true;
                    break;
                }
            }

            if( !stageEnabled ) return;

            mesh_.clearAddressingData();
            labelHashSet stageCells;
            polyMeshGenChecks::checkCellVolumes(mesh_, false, &stageCells);
            if( stageCells.size() == 0 ) return;
            const pointFieldPMG& pts  = mesh_.points();
            const cellListPMG&   cells = mesh_.cells();
            const faceListPMG&   faces = mesh_.faces();
            const std::string fname =
                "negVolCellCentres_" + stageName + ".csv";
            OFstream stageFile(fname);
            stageFile << "cellI,cx,cy,cz,nFaces,nUniquePoints" << nl;
            forAllConstIter(labelHashSet, stageCells, it)
            {
                const label cellI = it.key();
                if( cellI < 0 || cellI >= label(cells.size()) ) continue;
                const cell& c = cells[cellI];
                labelHashSet uniquePts;
                forAll(c, cfI)
                {
                    const label faceI = c[cfI];
                    if( faceI < 0 || faceI >= label(faces.size()) ) continue;
                    const face& f = faces[faceI];
                    forAll(f, fpI) uniquePts.insert(f[fpI]);
                }
                point cc = point::zero;
                label nUnique = 0;
                forAllConstIter(labelHashSet, uniquePts, pit)
                {
                    const label pI = pit.key();
                    if( pI < 0 || pI >= label(pts.size()) ) continue;
                    cc += pts[pI]; ++nUnique;
                }
                if( nUnique > 0 ) cc /= scalar(nUnique);
                stageFile << cellI << ","
                          << cc.x() << "," << cc.y() << "," << cc.z() << ","
                          << c.size() << "," << nUnique << nl;
            }
            Info << "Lineage [" << stageName.c_str() << "]: "
                 << stageCells.size() << " negVol cells -> "
                 << fname.c_str() << endl;
        };

        if( controller_.runCurrentStep("boundaryLayerGeneration") )
        {
            writeLineageCSV("preBL");
            generateBoundaryLayers();
            writeLineageCSV("postBLCreate");
        }

        if( controller_.runCurrentStep("meshOptimisation") )
        {
            optimiseFinalMesh();

            projectSurfaceAfterBackScaling();
            writeLineageCSV("postOptimize");
        }
        snapSurfaceBeforeBLRefinement();
        if( controller_.runCurrentStep("boundaryLayerRefinement") )
        {
            if( finalUntangleRejected_ )
            {
                Info << "refBoundaryLayers: skipped -- final untangle was rejected, mesh state unsafe" << endl;
            }
            else
            {
                refBoundaryLayers();
                writeLineageCSV("postRefBL");
            }
        }

        // Validate immediately before renumbering.
        {
            mesh_.clearAddressingData();
            labelHashSet negBeforeRenumber;
            labelHashSet pyrBeforeRenumber;
            polyMeshGenChecks::checkCellVolumes(mesh_, false, &negBeforeRenumber);
            polyMeshGenChecks::checkFacePyramids(mesh_, false, -SMALL, &pyrBeforeRenumber);
            Info << "Pre-renumber validation: negVol=" << negBeforeRenumber.size()
                 << " badPyramids=" << pyrBeforeRenumber.size() << endl;
        }

        renumberMesh();

        // Compact any points orphaned by renumbering before validation.
        {
            mesh_.clearAddressingData();
            const bool unusedAfterRenumber =
                polyMeshGenChecks::checkPoints(mesh_, false);
            if( unusedAfterRenumber )
            {
                Info << "Post-renumber cleanup: removing unused vertices" << endl;
                polyMeshGenModifier(mesh_).removeUnusedVertices();
                mesh_.clearAddressingData();
                const bool unusedAfterCleanup =
                    polyMeshGenChecks::checkPoints(mesh_, false);
                Info << "Post-renumber cleanup: unusedPoints bad->"
                     << (unusedAfterCleanup ? "bad" : "ok") << endl;
            }
            else
            {
                Info << "Post-renumber cleanup: no unused vertices found" << endl;
            }
        }

        // Validate immediately after renumbering.
        {
            mesh_.clearAddressingData();
            labelHashSet negAfterRenumber;
            labelHashSet pyrAfterRenumber;
            polyMeshGenChecks::checkCellVolumes(mesh_, false, &negAfterRenumber);
            polyMeshGenChecks::checkFacePyramids(mesh_, false, -SMALL, &pyrAfterRenumber);
            Info << "Post-renumber validation: negVol=" << negAfterRenumber.size()
                 << " badPyramids=" << pyrAfterRenumber.size() << endl;
        }

        replaceBoundaries();

        // FINAL DEBUG:
        // Validate mesh immediately after boundary replacement/renaming.
        // This tells us whether final checkMesh failures are already present
        // in-memory before mesh_.write(), or appear only after write/read.
        {
            labelHashSet finalNegCells;
            labelHashSet finalBadPyrFaces;

            polyMeshGenChecks::checkCellVolumes(mesh_, false, &finalNegCells);
            polyMeshGenChecks::checkFacePyramids
            (
                mesh_,
                false,
                -SMALL,
                &finalBadPyrFaces
            );

            Info << "FINAL internal validation after replaceBoundaries: "
                 << "negVol=" << finalNegCells.size()
                 << " badPyramids=" << finalBadPyrFaces.size()
                 << endl;
            writeLineageCSV("final");

            // Write final negative-volume cell centres for spatial diagnostics.
            if( finalNegCells.size() > 0 )
            {
                const pointFieldPMG& pts = mesh_.points();
                const cellListPMG& cells = mesh_.cells();
                const faceListPMG& faces = mesh_.faces();

                // Build attribution data for negVol spatial diagnosis.
                const meshSurfaceEngine mseNV(mesh_);
                const labelList& bPointsNV   = mseNV.boundaryPoints();
                const labelList& facePatchNV = mseNV.boundaryFacePatches();
                const VRWGraph& pointFacesNV = mseNV.pointFaces();
                const PtrList<boundaryPatch>& boundariesNV = mesh_.boundaries();

                // O(1) reverse map: mesh point label -> boundary point index
                labelList meshToBpNV(pts.size(), -1);
                forAll(bPointsNV, bpI)
                {
                    const label mpI = bPointsNV[bpI];
                    if( mpI >= 0 && mpI < label(meshToBpNV.size()) )
                        meshToBpNV[mpI] = bpI;
                }

                // Loser patch index set for O(1) lookup
                labelHashSet loserPatchIdxNV;
                forAll(boundariesNV, pI)
                    forAll(blGapLoserPatchNames_, ni)
                        if( boundariesNV[pI].patchName() == blGapLoserPatchNames_[ni] )
                            loserPatchIdxNV.insert(pI);

                // Gap action point positions for distance computation
                List<point> gapPtPos(blGapActionPoints_.size());
                {
                    label gi = 0;
                    forAllConstIter(labelHashSet, blGapActionPoints_, it)
                    {
                        const label mpI = it.key();
                        if( mpI >= 0 && mpI < label(pts.size()) )
                            gapPtPos[gi++] = pts[mpI];
                    }
                    gapPtPos.setSize(gi);
                }

                // Triple junction point positions
                List<point> tjPtPos(blblJunctionPoints_.size());
                {
                    label ti = 0;
                    forAllConstIter(labelHashSet, blblJunctionPoints_, it)
                    {
                        const label mpI = it.key();
                        if( mpI >= 0 && mpI < label(pts.size()) )
                            tjPtPos[ti++] = pts[mpI];
                    }
                    tjPtPos.setSize(ti);
                }

                OFstream negVolFile("negVolCellCentres.csv");
                negVolFile << "cellI,cx,cy,cz,nFaces,nUniquePoints,"
                           << "nearestPatch,isLoserPatch,"
                           << "distGapAction,distTripleJunction" << nl;

                forAllConstIter(labelHashSet, finalNegCells, it)
                {
                    const label cellI = it.key();
                    if( cellI < 0 || cellI >= label(cells.size()) )
                        continue;
                    const cell& c = cells[cellI];

                    // Collect unique mesh points for this cell
                    labelHashSet uniquePts;
                    forAll(c, cfI)
                    {
                        const label faceI = c[cfI];
                        if( faceI < 0 || faceI >= label(faces.size()) ) continue;
                        const face& f = faces[faceI];
                        forAll(f, fpI)
                            uniquePts.insert(f[fpI]);
                    }

                    // Compute cell centre
                    point cc = point::zero;
                    label nUnique = 0;
                    forAllConstIter(labelHashSet, uniquePts, pit)
                    {
                        const label pointI = pit.key();
                        if( pointI < 0 || pointI >= label(pts.size()) ) continue;
                        cc += pts[pointI];
                        ++nUnique;
                    }
                    if( nUnique > 0 )
                        cc /= scalar(nUnique);

                    // Nearest patch attribution using reverse map + all adjacent faces.
                    // Prefer loser patch if any adjacent face touches one.
                    label nearestPatchI = -1;
                    bool isLoser = false;
                    forAllConstIter(labelHashSet, uniquePts, pit)
                    {
                        if( isLoser ) break;
                        const label pointI = pit.key();
                        if( pointI < 0 || pointI >= label(meshToBpNV.size()) ) continue;
                        const label bpI = meshToBpNV[pointI];
                        if( bpI < 0 ) continue;
                        forAllRow(pointFacesNV, bpI, pfI)
                        {
                            const label bfI = pointFacesNV(bpI, pfI);
                            if( bfI < 0 || bfI >= label(facePatchNV.size()) ) continue;
                            const label pI = facePatchNV[bfI];
                            if( pI < 0 || pI >= label(boundariesNV.size()) ) continue;
                            if( nearestPatchI < 0 )
                                nearestPatchI = pI;
                            if( loserPatchIdxNV.found(pI) )
                            {
                                nearestPatchI = pI;
                                isLoser = true;
                                break;
                            }
                        }
                    }
                    const word nearestPatch =
                        nearestPatchI >= 0 ?
                        boundariesNV[nearestPatchI].patchName() :
                        word("unknown");

                    // Distance to nearest gap action point
                    scalar distGap = GREAT;
                    forAll(gapPtPos, gi)
                    {
                        const scalar d = mag(gapPtPos[gi] - cc);
                        if( d < distGap ) distGap = d;
                    }

                    // Distance to nearest triple junction point
                    scalar distTJ = GREAT;
                    forAll(tjPtPos, ti)
                    {
                        const scalar d = mag(tjPtPos[ti] - cc);
                        if( d < distTJ ) distTJ = d;
                    }

                    negVolFile << cellI << ","
                               << cc.x() << ","
                               << cc.y() << ","
                               << cc.z() << ","
                               << c.size() << ","
                               << nUnique << ","
                               << nearestPatch << ","
                               << (isLoser ? 1 : 0) << ","
                               << distGap << ","
                               << distTJ << nl;
                }
                Info << "negVol cell centres written to negVolCellCentres.csv"
                     << endl;
            }
        }

        // BL effective coverage report
        if( blLayerScale_.size() > 0 )
        {
            const meshSurfaceEngine mse(mesh_);
            const labelList& bPoints = mse.boundaryPoints();
            const VRWGraph& pointFaces = mse.pointFaces();
            const labelList& facePatch = mse.boundaryFacePatches();
            const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

            labelList patchTotal(boundaries.size(), 0);
            labelList patchBL(boundaries.size(), 0);

            forAll(bPoints, bpI)
            {
                const bool hasBL =
                    bpI < label(blLayerScale_.size())
                 && blLayerScale_[bpI] >= 0.01;

                labelHashSet countedPatches;
                forAllRow(pointFaces, bpI, pfI)
                {
                    const label pI = facePatch[pointFaces(bpI, pfI)];
                    if( pI < 0 || pI >= label(boundaries.size()) ) continue;
                    if( countedPatches.found(pI) ) continue;
                    countedPatches.insert(pI);
                    ++patchTotal[pI];
                    if( hasBL ) ++patchBL[pI];
                }
            }

            Info << "BL effective coverage per patch:" << endl;
            forAll(boundaries, patchI)
            {
                const scalar pct =
                    patchTotal[patchI] > 0
                  ? 100.0 * patchBL[patchI] / patchTotal[patchI]
                  : 0.0;
                Info << "  " << boundaries[patchI].patchName()
                     << ": " << patchBL[patchI] << "/" << patchTotal[patchI]
                     << " (" << pct << "%)" << endl;
            }
        }

        controller_.workflowCompleted();
    }
    catch(const std::string& message)
    {
        Info << message << endl;
    }
    catch(...)
    {
        WarningIn
        (
            "void cartesianMeshGenerator::generateMesh()"
        ) << "Meshing process terminated!" << endl;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

cartesianMeshGenerator::cartesianMeshGenerator(const Time& time)
:
    db_(time),
    surfacePtr_(NULL),
    modSurfacePtr_(NULL),
    meshDict_
    (
        IOobject
        (
            "meshDict",
            db_.system(),
            db_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    octreePtr_(NULL),
    mesh_(time),
    controller_(mesh_),
    finalUntangleRejected_(false),
    nPointsBeforeBL_(0)
{
    checkMeshDict cmd(meshDict_);

    fileName surfaceFile = meshDict_.lookup("surfaceFile");
    if( Pstream::parRun() )
        surfaceFile = ".."/surfaceFile;

    surfacePtr_ = new triSurf(db_.path()/surfaceFile);

    //- save meta data with the mesh (surface mesh + its topology info)
    triSurfaceMetaData sMetaData(*surfacePtr_);
    const dictionary& surfMetaDict = sMetaData.metaData();

    mesh_.metaData().add("surfaceFile", surfaceFile, true);
    mesh_.metaData().add("surfaceMeta", surfMetaDict, true);

    if( surfacePtr_->featureEdges().size() != 0 )
    {
        //- create surface patches based on the feature edges
        //- and update the meshDict based on the given data
        triSurfacePatchManipulator manipulator(*surfacePtr_);

        const triSurf* surfaceWithPatches =
            manipulator.surfaceWithPatches(&meshDict_);

        //- delete the old surface and assign the new one
        deleteDemandDrivenData(surfacePtr_);
        surfacePtr_ = surfaceWithPatches;
    }

    // Geometry preprocessing
    // Controls: active, reportOnly, weldNearPoints, weldTolerance,
    //   autoScaleWithCellSize, minFeatureToCellRatio, writeDiagnostics
    if( meshDict_.isDict("geometryPreprocessing") )
    {
        const dictionary& preProcDict =
            meshDict_.subDict("geometryPreprocessing");
        const bool active =
            preProcDict.found("active") ?
            bool(Switch(preProcDict.lookup("active"))) : false;
        if( active )
        {
            const bool reportOnly =
                preProcDict.found("reportOnly") ?
                bool(Switch(preProcDict.lookup("reportOnly"))) : false;
            const bool weldNearPoints =
                preProcDict.found("weldNearPoints") ?
                bool(Switch(preProcDict.lookup("weldNearPoints"))) : true;
            const bool autoScale =
                preProcDict.found("autoScaleWithCellSize") ?
                bool(Switch(preProcDict.lookup("autoScaleWithCellSize"))) : false;
            const scalar minRatio =
                preProcDict.found("minFeatureToCellRatio") ?
                readScalar(preProcDict.lookup("minFeatureToCellRatio")) : 0.20;
            const bool writeDiag =
                preProcDict.found("writeDiagnostics") ?
                bool(Switch(preProcDict.lookup("writeDiagnostics"))) : false;
            scalar weldTol =
                preProcDict.found("weldTolerance") ?
                readScalar(preProcDict.lookup("weldTolerance")) : 1e-4;
            // Auto-scale: weldTol = min(specified, ratio*minCellSize)
            if( autoScale )
            {
                scalar minCellSize = GREAT;
                if( meshDict_.found("minCellSize") )
                    minCellSize = readScalar(meshDict_.lookup("minCellSize"));
                else if( meshDict_.found("maxCellSize") )
                    minCellSize = readScalar(meshDict_.lookup("maxCellSize"));
                const scalar autoTol = minRatio * minCellSize;
                weldTol = Foam::min(weldTol, autoTol);
                Info << "Geometry preprocessing: autoScaleWithCellSize"
                     << " minCellSize=" << minCellSize
                     << " ratio=" << minRatio
                     << " => weldTolerance=" << weldTol << " m" << endl;
            }
            const label nPtsBefore = surfacePtr_->points().size();
            const label nTriBefore = surfacePtr_->size();
            Info << "Geometry preprocessing:" << endl;
            Info << "  Surface: " << nPtsBefore << " points, "
                 << nTriBefore << " facets" << endl;
            Info << "  weldTolerance: " << weldTol << " m" << endl;
            if( reportOnly )
            {
                // Non-mutating scan -- no geometry modification
                meshOctree* scanOctree = new meshOctree(*surfacePtr_);
                meshOctreeCreator
                (
                    *scanOctree,
                    meshDict_
                ).createOctreeBoxes();
                const label nFound = scanNearCoincidentPoints
                (
                    *surfacePtr_,
                    *scanOctree,
                    weldTol,
                    meshDict_,
                    writeDiag
                );
                deleteDemandDrivenData(scanOctree);
                Info << "  Near-coincident pairs found: " << nFound << endl;
                Info << "  No geometry modified (reportOnly true)" << endl;
            }
            else if( weldNearPoints )
            {
                const bool selectiveWeld =
                    preProcDict.found("selectiveWeld") ?
                    bool(Switch(preProcDict.lookup("selectiveWeld"))) : false;

                if( selectiveWeld )
                {
                    // Phase 2C: selective weld -- direct point scan, no octree.
                    // Building a temporary octree here triggers OMP-parallel lazy
                    // cache population on surfacePtr_ mutable members, racing with
                    // the main octree build and corrupting mesh quality.
                    Info << "  selectiveWeld: building approved pair list" << endl;

                    const pointField& pts = surfacePtr_->points();
                    const wordList pNames = surfacePtr_->patchNames();

                    // Build BL patch index set -- honor both per-patch and global nLayers
                    labelHashSet blPatchIdx;
                    if( meshDict_.isDict("boundaryLayers") )
                    {
                        const dictionary& bndL = meshDict_.subDict("boundaryLayers");
                        if( bndL.isDict("patchBoundaryLayers") )
                        {
                            const dictionary& pbl = bndL.subDict("patchBoundaryLayers");
                            forAll(pNames, pi)
                            {
                                if( pbl.isDict(pNames[pi]) )
                                {
                                    const dictionary& pd = pbl.subDict(pNames[pi]);
                                    const label nL = pd.found("nLayers") ?
                                        readLabel(pd.lookup("nLayers")) : 0;
                                    if( nL > 0 ) blPatchIdx.insert(pi);
                                }
                            }
                        }
                        // Global nLayers fallback -- mark all non-explicitly-zero patches as BL
                        if( bndL.found("nLayers") )
                        {
                            const label globalN = readLabel(bndL.lookup("nLayers"));
                            if( globalN > 0 )
                            {
                                forAll(pNames, pi)
                                {
                                    bool explicitlyZero = false;
                                    if( bndL.isDict("patchBoundaryLayers") )
                                    {
                                        const dictionary& pbl =
                                            bndL.subDict("patchBoundaryLayers");
                                        if( pbl.isDict(pNames[pi]) )
                                        {
                                            const dictionary& pd =
                                                pbl.subDict(pNames[pi]);
                                            if( pd.found("nLayers")
                                             && readLabel(pd.lookup("nLayers")) == 0 )
                                                explicitlyZero = true;
                                        }
                                    }
                                    if( !explicitlyZero )
                                        blPatchIdx.insert(pi);
                                }
                            }
                        }
                    }

                    // Build point-to-patch membership
                    List<labelHashSet> ptPatches(surfacePtr_->nPoints());
                    forAll(*surfacePtr_, triI)
                    {
                        const labelledTri& tri = (*surfacePtr_)[triI];
                        forAll(tri, vi)
                            if( tri[vi] >= 0 && tri[vi] < label(ptPatches.size()) )
                                ptPatches[tri[vi]].insert(tri.region());
                    }

                    // Build edge-adjacency set: pairs sharing a triangle edge must never weld
                    std::set<std::pair<label,label>> adjacentPairs;
                    forAll(*surfacePtr_, triI)
                    {
                        const labelledTri& tri = (*surfacePtr_)[triI];
                        for(label i=0; i<3; ++i)
                        {
                            const label v0 = tri[i];
                            const label v1 = tri[(i+1)%3];
                            if( v0 < 0 || v1 < 0 ) continue;
                            const label va = Foam::min(v0, v1);
                            const label vb = Foam::max(v0, v1);
                            if( va != vb )
                                adjacentPairs.insert(std::make_pair(va,vb));
                        }
                    }
                    Info << "  selectiveWeld: adjacency set built, "
                         << label(adjacentPairs.size()) << " edges" << endl;

                    // Build termination patch set once
                    wordHashSet termPatches;
                    if( meshDict_.isDict("boundaryLayers") )
                    {
                        const dictionary& bndL = meshDict_.subDict("boundaryLayers");
                        if( bndL.found("terminationPatches") )
                        {
                            wordList tp(bndL.lookup("terminationPatches"));
                            forAll(tp, i) termPatches.insert(tp[i]);
                        }
                    }

                    labelLongList newPointLabel(surfacePtr_->nPoints());
                    forAll(newPointLabel, pI) newPointLabel[pI] = pI;

                    label nFound = 0;
                    label nAdjacencyRejected = 0;
                    label nNeutralRejected = 0;
                    label nCrossPatchRejected = 0;

                    // Direct O(n^2) scan -- safe, deterministic, no octree needed
                    // At typical surface sizes (~5000 pts) this is <10M comparisons
                    for(label pI=0; pI<label(pts.size()); ++pI)
                    {
                        for(label pJ=pI+1; pJ<label(pts.size()); ++pJ)
                        {
                            if( magSqr(pts[pI]-pts[pJ]) >= sqr(weldTol) ) continue;

                            ++nFound;
                            const label a = pI;
                            const label b = pJ;

                            DynList<word> pNA, pNB;
                            forAllConstIter(labelHashSet, ptPatches[a], it2)
                            {
                                const label r = it2.key();
                                if( r >= 0 && r < label(pNames.size()) )
                                    pNA.append(pNames[r]);
                            }
                            forAllConstIter(labelHashSet, ptPatches[b], it2)
                            {
                                const label r = it2.key();
                                if( r >= 0 && r < label(pNames.size()) )
                                    pNB.append(pNames[r]);
                            }

                            const scalar ratio = mag(pts[pI]-pts[pJ]) / weldTol;

                            // Reject if points share a triangle edge
                            const bool adjacent =
                                adjacentPairs.count(std::make_pair(a,b)) > 0;
                            if( adjacent ) { ++nAdjacencyRejected; continue; }

                            // Reject neutral/termination patch touches
                            bool touchesNeutral = false;
                            forAll(pNA, pi)
                            {
                                bool isBLorTerm = false;
                                forAll(pNames, ni)
                                    if( pNames[ni] == pNA[pi] )
                                        if( blPatchIdx.found(ni) || termPatches.found(pNA[pi]) )
                                            isBLorTerm = true;
                                if( !isBLorTerm ) touchesNeutral = true;
                            }
                            forAll(pNB, pi)
                            {
                                bool isBLorTerm = false;
                                forAll(pNames, ni)
                                    if( pNames[ni] == pNB[pi] )
                                        if( blPatchIdx.found(ni) || termPatches.found(pNB[pi]) )
                                            isBLorTerm = true;
                                if( !isBLorTerm ) touchesNeutral = true;
                            }
                            if( touchesNeutral )
                            {
                                ++nNeutralRejected;
                                continue;
                            }

                            // Require at least one shared BL patch.
                            // Identical-patch-set is too strict for triple-junction
                            // points (blade/hub/periodic) where the two near-coincident
                            // points legitimately have different patch memberships.
                            // Neutral/termination touches already rejected above.
                            bool sharedBLPatch = false;
                            forAllConstIter(labelHashSet, ptPatches[a], iterA)
                            {
                                if( blPatchIdx.found(iterA.key())
                                 && ptPatches[b].found(iterA.key()) )
                                {
                                    sharedBLPatch = true;
                                    break;
                                }
                            }
                            if( !sharedBLPatch )
                            {
                                ++nCrossPatchRejected;
                                continue;
                            }

                            if( !pNA.size() || !pNB.size() ) continue;

                            if( ratio < 0.50 )
                            {
                                newPointLabel[b] = a;
                                Info << "  Approved weld pair: " << a
                                     << " " << b
                                     << " ratio=" << ratio << endl;
                            }
                        }
                    }

                    Info << "  selectiveWeld: " << nFound << " pairs scanned" << endl;
                    Info << "  selectiveWeld: adjacency-rejected = "
                         << nAdjacencyRejected << endl;
                    Info << "  selectiveWeld: neutral-rejected = "
                         << nNeutralRejected << endl;
                    Info << "  selectiveWeld: cross-patch-rejected = "
                         << nCrossPatchRejected << endl;

                    bool anyApproved = false;
                    forAll(newPointLabel, pI)
                        if( newPointLabel[pI] != pI ) { anyApproved = true; break; }

                    if( anyApproved )
                    {
                        triSurfaceCleanupDuplicates cleaner(*surfacePtr_, weldTol);
                        cleaner.mergeApprovedPairs(newPointLabel);
                    }
                    else
                    {
                        Info << "  selectiveWeld: no approved pairs -- surface unchanged" << endl;
                    }
                }
                else
                {
                    // Full weld â proven safe path
                    meshOctree* preprocOctree = new meshOctree(*surfacePtr_);
                    meshOctreeCreator
                    (
                        *preprocOctree,
                        meshDict_
                    ).createOctreeBoxes();
                    triSurfaceCleanupDuplicates cleaner(*preprocOctree, weldTol);
                    cleaner.mergeIdentities();
                    deleteDemandDrivenData(preprocOctree);
                }
                const label nPtsAfter = surfacePtr_->points().size();
                const label nTriAfter = surfacePtr_->size();
                Info << "  After:  " << nPtsAfter << " points, "
                     << nTriAfter << " facets" << endl;
                Info << "  Welded: " << (nPtsBefore - nPtsAfter)
                     << " points, removed "
                     << (nTriBefore - nTriAfter)
                     << " degenerate facets" << endl;
                if( writeDiag )
                    Info << "  writeDiagnostics: VTK output not yet implemented" << endl;
            }
        }
    }

    if( meshDict_.found("anisotropicSources") )
    {
        surfaceMeshGeometryModification surfMod(*surfacePtr_, meshDict_);

        modSurfacePtr_ = surfMod.modifyGeometry();

        octreePtr_ = new meshOctree(*modSurfacePtr_);
    }
    else
    {
        octreePtr_ = new meshOctree(*surfacePtr_);
    }

    meshOctreeCreator(*octreePtr_, meshDict_).createOctreeBoxes();

    generateMesh();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cartesianMeshGenerator::~cartesianMeshGenerator()
{
    deleteDemandDrivenData(surfacePtr_);
    deleteDemandDrivenData(modSurfacePtr_);
    deleteDemandDrivenData(octreePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void cartesianMeshGenerator::writeMesh() const
{
    mesh_.write();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
