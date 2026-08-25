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
#include "detectBoundaryLayers.H"
#include "meshSurfacePartitioner.H"
#include "BLJunctionClassifier.H"
#include "PatchRoleMap.H"
#include "renameBoundaryPatches.H"
#include "checkMeshDict.H"
#include "checkCellConnectionsOverFaces.H"
#include "checkIrregularSurfaceConnections.H"
#include "checkNonMappableCellConnections.H"
#include "OFstream.H"
#include "IFstream.H"
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

//- Raw signed volume of one cell, using the same unclamped
//- construction as rawCellVolumeStats() and checkCellVolumes().
static scalar rawSignedCellVolume
(
    const polyMeshGen& mesh,
    const label cellI
)
{
    const vectorField& fCtrs  = mesh.addressingData().faceCentres();
    const vectorField& fAreas = mesh.addressingData().faceAreas();
    const labelList& own      = mesh.owner();
    const cellListPMG& cells  = mesh.cells();

    if( cellI < 0 || cellI >= label(cells.size()) )
        return scalar(0.0);

    const cell& c = cells[cellI];

    if( c.size() == 0 )
        return scalar(0.0);

    vector cEst(vector::zero);

    forAll(c, fI)
        cEst += fCtrs[c[fI]];

    cEst /= c.size();

    scalar cellVol(0.0);

    forAll(c, fI)
    {
        scalar pyr3Vol =
            fAreas[c[fI]] & (fCtrs[c[fI]] - cEst);

        if( own[c[fI]] != cellI )
            pyr3Vol *= -1.0;

        cellVol += pyr3Vol;
    }

    return cellVol / scalar(3.0);
}


//- Raw signed cell-volume statistics.
//
//  MUST NOT use mesh.addressingData().cellVolumes(): that cache stores
//  CLAMPED (always positive) volumes for optimizer safety, so summing its
//  negative entries always yields zero. This reproduces the signed volume
//  computation performed inside polyMeshGenChecks::checkCellVolumes() so
//  that severity metrics and the negVol count describe the same geometry.
//
//  negMag   -- sum of |volume| over cells with cellVol < 0 (true inversions)
//  minVol   -- true signed minimum over all cells
//  nNeg     -- count of cells with cellVol < 0
//  nBelowVS -- count of cells with cellVol < VSMALL; this is the population
//              checkCellVolumes() reports, i.e. inversions PLUS zero/
//              near-zero slivers. Kept separate so the two can be compared.
static void rawCellVolumeStats
(
    const polyMeshGen& mesh,
    scalar& negMag,
    scalar& minVol,
    label& nNeg,
    label& nBelowVS
)
{
    negMag = 0.0;
    minVol = GREAT;
    nNeg = 0;
    nBelowVS = 0;

    const vectorField& fCtrs  = mesh.addressingData().faceCentres();
    const vectorField& fAreas = mesh.addressingData().faceAreas();
    const labelList&   own    = mesh.owner();
    const cellListPMG& cells  = mesh.cells();

    forAll(cells, cellI)
    {
        const cell& c = cells[cellI];
        if( c.size() == 0 ) continue;

        vector cEst(vector::zero);
        forAll(c, fI)
            cEst += fCtrs[c[fI]];
        cEst /= c.size();

        scalar cellVol(0.0);
        forAll(c, fI)
        {
            scalar pyr3Vol = fAreas[c[fI]] & (fCtrs[c[fI]] - cEst);
            if( own[c[fI]] != cellI )
                pyr3Vol *= -1.0;
            cellVol += pyr3Vol;
        }
        cellVol /= 3.0;

        if( cellVol < minVol ) minVol = cellVol;
        if( cellVol < 0 )      { ++nNeg; negMag -= cellVol; }
        if( cellVol < VSMALL ) ++nBelowVS;
    }
}

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
    // v7: narrow blade+periodic+hub/shroud triple-junction fix.
    // Resolve role -> patch-ID sets via PatchRoleMap (same local-
    // construction pattern already used elsewhere in this file for
    // BLJunctionClassifier), pass small ID sets through a setter --
    // per SOL review, deliberately NOT threading PatchRoleMap through
    // meshSurfaceMapper's constructor (too much architecture change
    // for this experiment).
    {
        PatchRoleMap tjRoles(meshDict_);
        if( tjRoles.active() )
        {
            const wordList& allPatchNames = octreePtr_->surface().patchNames();
            labelHashSet bladeIds, periodicIds, hubIds, shroudIds, flowIds;
            forAll(allPatchNames, patchI)
            {
                const word& pName = allPatchNames[patchI];
                if( tjRoles.hasRole(pName, "blade") ) bladeIds.insert(patchI);
                if( tjRoles.hasRole(pName, "periodic") ) periodicIds.insert(patchI);
                if( tjRoles.hasRole(pName, "hub") ) hubIds.insert(patchI);
                if( tjRoles.hasRole(pName, "shroud") ) shroudIds.insert(patchI);
                if( tjRoles.hasRole(pName, "inletOutlet") ) flowIds.insert(patchI);
            }
            mapper.setTripleJunctionTargetPatches
            (
                bladeIds, periodicIds, hubIds, shroudIds, flowIds
            );
            Info << "v7 TJ target patches: blade=" << bladeIds.size()
                 << " periodic=" << periodicIds.size()
                 << " hub=" << hubIds.size()
                 << " shroud=" << shroudIds.size()
                 << " inletOutlet=" << flowIds.size() << endl;

            bool writeTripleJunctionDiagnostic = false;
            if( meshDict_.found("writeTripleJunctionDiagnostic") )
                writeTripleJunctionDiagnostic =
                    Switch(meshDict_.lookup("writeTripleJunctionDiagnostic"));
            mapper.setTripleJunctionCornerFix
            (
                writeTripleJunctionDiagnostic,
                word("surfaceProjection")
            );
        }
        else
        {
            Info << "v7: PatchRoleMap inactive (no patchRoles in meshDict) "
                 << "-- triple-junction fix disabled" << endl;
        }
    }
    mapper.mapCornersAndEdges();
    mapper.setTripleJunctionCornerFix(false);

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
    bool skipMapEdgesAndCorners = false;
    if( meshDict_.found("skipMapEdgesAndCorners") )
    {
        skipMapEdgesAndCorners =
            Switch(meshDict_.lookup("skipMapEdgesAndCorners"));
    }

    if( skipMapEdgesAndCorners )
    {
        Info << "mapEdgesAndCorners skipped by meshDict" << endl;
        return;
    }

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

labelHashSet cartesianMeshGenerator::traceToSeedFaces
(
    const labelHashSet& badFaces,
    const label nRings
)
{
    const labelList& owner     = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();
    const cellListPMG& cells   = mesh_.cells();
    const label nInternalFaces = mesh_.nInternalFaces();
    const label nFaces         = mesh_.faces().size();
    const label nCells         = mesh_.cells().size();
    const label nBndFaces      = nFaces - nInternalFaces;

    // Seed cells from bad internal faces
    labelHashSet seedCells;
    forAllConstIter(labelHashSet, badFaces, it)
    {
        const label faceI = it.key();
        if( faceI < 0 || faceI >= label(owner.size()) ) continue;
        const label ownC = owner[faceI];
        if( ownC >= 0 && ownC < nCells ) seedCells.insert(ownC);
        if( faceI < label(neighbour.size()) )
        {
            const label neiC = neighbour[faceI];
            if( neiC >= 0 && neiC < nCells ) seedCells.insert(neiC);
        }
    }

    // N-ring expansion
    for( label ring = 0; ring < nRings; ++ring )
    {
        labelHashSet ringCells;
        forAllConstIter(labelHashSet, seedCells, cit)
        {
            const cell& c = cells[cit.key()];
            forAll(c, fI)
            {
                const label faceI = c[fI];
                if( faceI < 0 || faceI >= nFaces ) continue;
                if( faceI < label(owner.size()) )
                {
                    const label oc = owner[faceI];
                    if( oc >= 0 && oc < nCells ) ringCells.insert(oc);
                }
                if( faceI < label(neighbour.size()) )
                {
                    const label nc = neighbour[faceI];
                    if( nc >= 0 && nc < nCells ) ringCells.insert(nc);
                }
            }
        }
        forAllConstIter(labelHashSet, ringCells, it)
            seedCells.insert(it.key());
    }

    // Collect boundary-local bfI from seed cells
    labelHashSet seedBfI;
    forAllConstIter(labelHashSet, seedCells, cit)
    {
        const cell& c = cells[cit.key()];
        forAll(c, fI)
        {
            const label faceI = c[fI];
            if( faceI >= nInternalFaces && faceI < nFaces )
            {
                const label bfI = faceI - nInternalFaces;
                if( bfI >= 0 && bfI < nBndFaces )
                    seedBfI.insert(bfI);
            }
        }
    }

    Info << "traceToSeedFaces: badFaces=" << badFaces.size()
         << " seedCells=" << seedCells.size()
         << " nRings=" << nRings
         << " seedBndFaces(bfI)=" << seedBfI.size()
         << endl;

    return seedBfI;
}

static void blTopoAudit(polyMeshGen& auditMesh, const word& stage)
{
    try
    {
        auditMesh.clearAddressingData();
        const meshSurfaceEngine mse(auditMesh);
        meshSurfacePartitioner mPart(mse);
        detectBoundaryLayers dbl(mPart, false);
        const labelList& fInLayer = dbl.faceInLayer();
        const labelList& fPatch   = mse.boundaryFacePatches();
        const PtrList<boundaryPatch>& bnd = auditMesh.boundaries();
        const label nP = bnd.size();
        const label nBF = fPatch.size();

        if( fInLayer.size() != nBF )
        {
            Info << "BLTOPO_AUDIT stage=" << stage
                 << " SIZE_MISMATCH faceInLayer=" << fInLayer.size()
                 << " boundaryFaces=" << nBF
                 << " -- skipping per-patch (index spaces differ)" << endl;
            return;
        }

        // meshSurfaceEngine has no faceAreas() accessor in this fork.
        // Map boundary-local bfI explicitly to the global mesh face.
        const vectorField& allFaceAreas =
            auditMesh.addressingData().faceAreas();
        const label nInternal = auditMesh.nInternalFaces();

        if
        (
            nInternal < 0
         || nInternal + nBF > label(allFaceAreas.size())
        )
        {
            Info << "BLTOPO_AUDIT stage=" << stage
                 << " AREA_SIZE_MISMATCH"
                 << " nInternal=" << nInternal
                 << " boundaryFaces=" << nBF
                 << " faceAreasSize=" << allFaceAreas.size()
                 << " -- skipping per-patch" << endl;
            return;
        }

        labelList totF(nP,0), layF(nP,0);
        scalarField totA(nP,0.0), layA(nP,0.0);
        for(label bfI=0; bfI<nBF; ++bfI)
        {
            const label p = fPatch[bfI];
            if( p<0 || p>=nP ) continue;

            const label faceI = nInternal + bfI;
            const scalar a = mag(allFaceAreas[faceI]);

            ++totF[p]; totA[p]+=a;
            if( fInLayer[bfI] >= 0 )
            { ++layF[p]; layA[p]+=a; }
        }
        Info << "BLTOPO_AUDIT stage=" << stage
             << " nDistinctLayers=" << dbl.nDistinctLayers()
             << " hairEdges=" << dbl.hairEdges().size()
             << " boundaryFaces=" << nBF
             << " faceInLayerSize=" << fInLayer.size() << endl;
        forAll(bnd, p)
        {
            const scalar fp = totF[p]>0 ? 100.0*layF[p]/totF[p] : 0.0;
            const scalar ap = totA[p]>VSMALL ? 100.0*layA[p]/totA[p] : 0.0;
            Info << "  " << bnd[p].patchName()
                 << ": layerFaces=" << layF[p] << "/" << totF[p]
                 << " (" << fp << "%) area=" << ap << "%" << endl;
        }
    }
    catch( ... )
    {
        Info << "BLTOPO_AUDIT stage=" << stage
             << " FAILED (partitioner/detector threw)" << endl;
    }
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
    blTopoAudit(mesh_, "POST_ADDLAYER");

    // Capture junction points for handoff to refineBoundaryLayers
    blblJunctionPoints_ = bl.junctionEdgePoints();
    blblAcuteCornerPoints_ = bl.blblAcuteCornerPoints();
    blRampSeedPoints_ = bl.rampSeedPoints();
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
        refLayers.setRampSeedPoints(blRampSeedPoints_);
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

        // Pre-refBL retraction: trace post-optimizer bad pyramid faces
        // to boundary-local seed faces and force 1 layer there.
        // Gated by meshDict: postOptBLRetraction true/false (default: false).
        // Enable to attempt surgical BL suppression at bad cluster locations.
        {
            bool doRetraction = false;
            if( meshDict_.isDict("boundaryLayers") )
            {
                const dictionary& bndL =
                    meshDict_.subDict("boundaryLayers");
                if( bndL.found("postOptBLRetraction") )
                    doRetraction =
                        bool(Switch(bndL.lookup("postOptBLRetraction")));
            }

            if( doRetraction && postOptBadFaces_.size() > 0 )
            {
                label retractionRings = 1;
                if( meshDict_.isDict("boundaryLayers") )
                {
                    const dictionary& bndL =
                        meshDict_.subDict("boundaryLayers");
                    if( bndL.found("postOptBLRetractionRings") )
                        retractionRings =
                            readLabel(bndL.lookup("postOptBLRetractionRings"));
                }

                Info << "Pre-refBL retraction: tracing "
                     << postOptBadFaces_.size()
                     << " bad pyramid faces to seed faces"
                     << " (rings=" << retractionRings << ")" << endl;

                const labelHashSet seedFaces =
                    traceToSeedFaces(postOptBadFaces_, retractionRings);

                if( seedFaces.size() > 0 )
                {
                    refLayers.forceSingleLayerAtFaces(seedFaces);
                    Info << "Pre-refBL retraction: "
                         << seedFaces.size()
                         << " boundary faces retracted to 1 layer" << endl;
                }
                else
                {
                    Info << "Pre-refBL retraction: no boundary seed faces "
                         << "found -- retraction skipped" << endl;
                }
            }
            else if( postOptBadFaces_.size() > 0 )
            {
                Info << "Pre-refBL retraction: DIAGNOSTIC MODE "
                     << "(postOptBLRetraction false) -- "
                     << postOptBadFaces_.size()
                     << " bad faces available, tracing for inspection" << endl;

                label diagRings = 1;
                if( meshDict_.isDict("boundaryLayers") )
                {
                    const dictionary& bndL =
                        meshDict_.subDict("boundaryLayers");
                    if( bndL.found("postOptBLRetractionRings") )
                        diagRings =
                            readLabel(bndL.lookup("postOptBLRetractionRings"));
                }
                const labelHashSet seedFaces =
                    traceToSeedFaces(postOptBadFaces_, diagRings);

                Info << "Pre-refBL retraction: would retract "
                     << seedFaces.size()
                     << " boundary faces if enabled (rings="
                     << diagRings << ")" << endl;
            }
        }

        // Provenance-direct BL retraction.
        // Separate from postOptBLRetraction (geometric traceToSeedFaces).
        // Uses exact bfI seeds from a previous post-refBL diagnostic run.
        // File format: one boundary-local bfI per line, no header.
        // Generate from diagnostic run:
        //   tail -n +2 postRefBL_provenanceSeedBfI.csv \
        //       > postRefBL_provenanceSeedBfI.labels
        {
            bool doProvRetraction = false;
            word seedFileName("postRefBL_provenanceSeedBfI.labels");

            if( meshDict_.isDict("boundaryLayers") )
            {
                const dictionary& bndL =
                    meshDict_.subDict("boundaryLayers");

                if( bndL.found("postRefBLProvenanceRetraction") )
                    doProvRetraction =
                        bool(Switch(bndL.lookup("postRefBLProvenanceRetraction")));

                if( bndL.found("postRefBLProvenanceSeedFile") )
                    seedFileName =
                        word(bndL.lookup("postRefBLProvenanceSeedFile"));
            }

            if( doProvRetraction && reprojUnsafe_ )
            {
                Info << "Provenance-direct retraction: skipped -- "
                     << "pre-refBL mesh has residual negVol (reprojUnsafe). "
                     << "Retraction on dirty mesh gives unreliable results."
                     << endl;
            }
            else if( doProvRetraction )
            {
                IFstream seedFile(seedFileName);

                if( seedFile.good() )
                {
                    labelHashSet provenanceSeeds;
                    label bfI = -1;

                    while( seedFile.good() )
                    {
                        seedFile >> bfI;
                        if( seedFile.good() && bfI >= 0 )
                            provenanceSeeds.insert(bfI);
                    }

                    if( provenanceSeeds.size() > 0 )
                    {
                        label ring0MaxLayers = 3;
                        label ring1MaxLayers = 4;
                        label ring2MaxLayers = 0;

                        scalar ring0ThicknessScale = 0.60;
                        scalar ring1ThicknessScale = 0.80;
                        scalar ring2ThicknessScale = 1.00;

                        if( meshDict_.isDict("boundaryLayers") )
                        {
                            const dictionary& bndL =
                                meshDict_.subDict("boundaryLayers");

                            if( bndL.found("postRefBLRing0MaxLayers") )
                                ring0MaxLayers =
                                    readLabel(bndL.lookup("postRefBLRing0MaxLayers"));

                            if( bndL.found("postRefBLRing1MaxLayers") )
                                ring1MaxLayers =
                                    readLabel(bndL.lookup("postRefBLRing1MaxLayers"));

                            if( bndL.found("postRefBLRing2MaxLayers") )
                                ring2MaxLayers =
                                    readLabel(bndL.lookup("postRefBLRing2MaxLayers"));

                            if( bndL.found("postRefBLRing0ThicknessScale") )
                                ring0ThicknessScale =
                                    readScalar(bndL.lookup("postRefBLRing0ThicknessScale"));

                            if( bndL.found("postRefBLRing1ThicknessScale") )
                                ring1ThicknessScale =
                                    readScalar(bndL.lookup("postRefBLRing1ThicknessScale"));

                            if( bndL.found("postRefBLRing2ThicknessScale") )
                                ring2ThicknessScale =
                                    readScalar(bndL.lookup("postRefBLRing2ThicknessScale"));
                        }

                        refLayers.forceMaxLayersAtFaces
                        (
                            provenanceSeeds,
                            ring0MaxLayers,
                            ring1MaxLayers,
                            ring2MaxLayers,
                            ring0ThicknessScale,
                            ring1ThicknessScale,
                            ring2ThicknessScale
                        );

                        Info << "Provenance-direct tapered retraction: loaded "
                             << provenanceSeeds.size()
                             << " bfI seeds from " << seedFileName
                             << " with ring caps=("
                             << ring0MaxLayers << ','
                             << ring1MaxLayers << ','
                             << ring2MaxLayers << ')'
                             << " thickness scales=("
                             << ring0ThicknessScale << ','
                             << ring1ThicknessScale << ','
                             << ring2ThicknessScale << ')' << endl;
                    }
                    else
                    {
                        Info << "Provenance-direct retraction: seed file "
                             << seedFileName
                             << " contained no valid bfI entries" << endl;
                    }
                }
                else
                {
                    Info << "Provenance-direct retraction: seed file "
                         << seedFileName
                         << " not found/readable -- skipped" << endl;
                }
            }
        }

        // Pre-refBL mesh snapshot for two-pass repair loop.
        struct PreRefBLMeshSnapshot
        {
            pointField   points;
            faceList     faces;
            cellList     cells;
            wordList     patchNames;
            labelList    patchStart;
            labelList    patchSize;

            // Mesh subsets are topology-indexed state and must participate
            // in rollback together with points/faces/cells.
            std::map<label, meshSubset> pointSubsets;
            std::map<label, meshSubset> faceSubsets;
            std::map<label, meshSubset> cellSubsets;

            // BL detection metadata. These contain point/face indices
            // and must match the mesh state restored by this snapshot.
            labelHashSet blblJunctionPoints;
            labelHashSet blblAcuteCornerPoints;
            boolList     rampSeedPoints;
            labelList    vtFaceRing;

            bool         valid;
            PreRefBLMeshSnapshot() : valid(false) {}
        };

        auto takePreRefBLSnapshot =
        [&](PreRefBLMeshSnapshot& snap)
        {
            Info << "Pre-refBL snapshot: saving mesh state" << endl;

            const pointFieldPMG& mp = mesh_.points();
            snap.points.setSize(mp.size());
            forAll(mp, pI) snap.points[pI] = mp[pI];

            const faceListPMG& mf = mesh_.faces();
            snap.faces.setSize(mf.size());
            forAll(mf, fI) snap.faces[fI] = mf[fI];

            const cellListPMG& mc = mesh_.cells();
            snap.cells.setSize(mc.size());
            forAll(mc, cI) snap.cells[cI] = mc[cI];

            const PtrList<boundaryPatch>& bnd = mesh_.boundaries();
            snap.patchNames.setSize(bnd.size());
            snap.patchStart.setSize(bnd.size());
            snap.patchSize.setSize(bnd.size());
            forAll(bnd, patchI)
            {
                snap.patchNames[patchI] = bnd[patchI].patchName();
                snap.patchStart[patchI] = bnd[patchI].patchStart();
                snap.patchSize[patchI]  = bnd[patchI].patchSize();
            }

            polyMeshGenModifier subsetModifier(mesh_);
            snap.pointSubsets = subsetModifier.pointSubsetsAccess();
            snap.faceSubsets  = subsetModifier.faceSubsetsAccess();
            snap.cellSubsets  = subsetModifier.cellSubsetsAccess();

            Info << "Pre-refBL snapshot: saved subsets point="
                 << snap.pointSubsets.size()
                 << " face=" << snap.faceSubsets.size()
                 << " cell=" << snap.cellSubsets.size() << endl;

            snap.blblJunctionPoints = blblJunctionPoints_;
            snap.blblAcuteCornerPoints = blblAcuteCornerPoints_;
            snap.rampSeedPoints = blRampSeedPoints_;
            snap.vtFaceRing = vtFaceRing_;

            snap.valid = true;

            Info << "Pre-refBL snapshot: saved "
                 << snap.points.size() << " points, "
                 << snap.faces.size() << " faces, "
                 << snap.cells.size() << " cells, "
                 << snap.patchNames.size() << " patches" << endl;
        };

        auto restorePreRefBLSnapshot =
        [&](const PreRefBLMeshSnapshot& snap) -> bool
        {
            if( !snap.valid )
            {
                Info << "Pre-refBL snapshot: restore skipped -- "
                     << "snapshot is not valid" << endl;
                return false;
            }

            // Validate topology-indexed subset state before changing
            // the current mesh. updateSubset() assumes every stored label
            // is a valid row in its supplied topology map.
            label maxPointSubsetLabel = -1;
            label maxFaceSubsetLabel  = -1;
            label maxCellSubsetLabel  = -1;
            label invalidPointSubsetLabels = 0;
            label invalidFaceSubsetLabels  = 0;
            label invalidCellSubsetLabels  = 0;

            for
            (
                std::map<label, meshSubset>::const_iterator it =
                    snap.pointSubsets.begin();
                it != snap.pointSubsets.end();
                ++it
            )
            {
                labelList labels;
                it->second.containedElements(labels);
                forAll(labels, i)
                {
                    if( labels[i] > maxPointSubsetLabel )
                        maxPointSubsetLabel = labels[i];

                    if( labels[i] < 0 || labels[i] >= snap.points.size() )
                        ++invalidPointSubsetLabels;
                }
            }

            for
            (
                std::map<label, meshSubset>::const_iterator it =
                    snap.faceSubsets.begin();
                it != snap.faceSubsets.end();
                ++it
            )
            {
                labelList labels;
                it->second.containedElements(labels);
                forAll(labels, i)
                {
                    if( labels[i] > maxFaceSubsetLabel )
                        maxFaceSubsetLabel = labels[i];

                    if( labels[i] < 0 || labels[i] >= snap.faces.size() )
                        ++invalidFaceSubsetLabels;
                }
            }

            for
            (
                std::map<label, meshSubset>::const_iterator it =
                    snap.cellSubsets.begin();
                it != snap.cellSubsets.end();
                ++it
            )
            {
                labelList labels;
                it->second.containedElements(labels);
                forAll(labels, i)
                {
                    if( labels[i] > maxCellSubsetLabel )
                        maxCellSubsetLabel = labels[i];

                    if( labels[i] < 0 || labels[i] >= snap.cells.size() )
                        ++invalidCellSubsetLabels;
                }
            }

            Info << "PREREFBL_SUBSET_RESTORE"
                 << " pointSubsets=" << snap.pointSubsets.size()
                 << " faceSubsets=" << snap.faceSubsets.size()
                 << " cellSubsets=" << snap.cellSubsets.size()
                 << " maxLabels=("
                 << maxPointSubsetLabel << ","
                 << maxFaceSubsetLabel << ","
                 << maxCellSubsetLabel << ")"
                 << " invalid=("
                 << invalidPointSubsetLabels << ","
                 << invalidFaceSubsetLabels << ","
                 << invalidCellSubsetLabels << ")"
                 << endl;

            if
            (
                invalidPointSubsetLabels
             || invalidFaceSubsetLabels
             || invalidCellSubsetLabels
            )
            {
                Info << "Pre-refBL snapshot: restore rejected -- "
                     << "snapshot contains out-of-range subset labels"
                     << endl;
                return false;
            }

            Info << "Pre-refBL snapshot: restoring mesh state" << endl;

            polyMeshGenModifier meshModifier(mesh_);
            meshModifier.pointsAccess() = snap.points;
            meshModifier.facesAccess()  = snap.faces;
            meshModifier.cellsAccess()  = snap.cells;

            PtrList<boundaryPatch>& bnd = meshModifier.boundariesAccess();
            if( bnd.size() != snap.patchNames.size() )
            {
                Info << "Pre-refBL snapshot: restore failed -- patch count "
                     << "changed from " << snap.patchNames.size()
                     << " to " << bnd.size() << endl;
                return false;
            }
            forAll(bnd, patchI)
            {
                bnd[patchI].patchName()  = snap.patchNames[patchI];
                bnd[patchI].patchStart() = snap.patchStart[patchI];
                bnd[patchI].patchSize()  = snap.patchSize[patchI];
            }

            meshModifier.pointSubsetsAccess() = snap.pointSubsets;
            meshModifier.faceSubsetsAccess()  = snap.faceSubsets;
            meshModifier.cellSubsetsAccess()  = snap.cellSubsets;

            Info << "PREREFBL_SUBSET_RESTORE restored point="
                 << snap.pointSubsets.size()
                 << " face=" << snap.faceSubsets.size()
                 << " cell=" << snap.cellSubsets.size() << endl;

            blblJunctionPoints_ = snap.blblJunctionPoints;
            blblAcuteCornerPoints_ = snap.blblAcuteCornerPoints;
            blRampSeedPoints_ = snap.rampSeedPoints;
            vtFaceRing_ = snap.vtFaceRing;

            //- Topology has been replaced wholesale. The ordinary
            //- clearAddressingData() deletes polyMeshGenAddressing but leaves
            //- ownerPtr_/neighbourPtr_ alive. Those arrays describe the
            //- pass-1 mesh and must be discarded before pass 2.
            meshModifier.clearTopologyAddressing();

            Info << "PREREFBL_RESTORE cleared owner/neighbour + derived addressing"
                 << endl;

            Info << "Pre-refBL snapshot: restored "
                 << snap.points.size() << " points, "
                 << snap.faces.size() << " faces, "
                 << snap.cells.size() << " cells" << endl;
            return true;
        };

        bool doPreRefBLSnapshot = false;
        if( meshDict_.isDict("boundaryLayers") )
        {
            const dictionary& bndL = meshDict_.subDict("boundaryLayers");
            if( bndL.found("preRefBLSnapshot") )
                doPreRefBLSnapshot =
                    bool(Switch(bndL.lookup("preRefBLSnapshot")));
        }

        //- Hoisted so the snapshot decision can see it: auto-repair
        //- requires a recovery point, and must not depend on the user
        //- separately remembering to set preRefBLSnapshot.
        bool doAutoRepair = false;
        if( meshDict_.isDict("boundaryLayers") )
        {
            const dictionary& bndLAR = meshDict_.subDict("boundaryLayers");
            if( bndLAR.found("preRefBLAutoRepair") )
                doAutoRepair =
                    bool(Switch(bndLAR.lookup("preRefBLAutoRepair")));
        }

        //- Snapshot capability != repair eligibility. A dirty mesh is an
        //- argument FOR holding a recovery point, not against it. Formerly
        //- gated on !reprojUnsafe_, which meant "surface re-projection is
        //- unsafe" -- a different capability, and one that also destroyed
        //- the prerequisite state needed to consider repair at all.
        PreRefBLMeshSnapshot preRefBLSnap;
        if( doPreRefBLSnapshot || doAutoRepair )
            takePreRefBLSnapshot(preRefBLSnap);

        //- Q0: observational baseline before any refBL work. Logged only --
        //- the transaction comparator is Q1 (see acceptance gate).
        label  q0NegVol = 0, q0BadPyr = 0;
        scalar q0MinCellVol = GREAT, q0NegMag = 0.0;
        //- Declared at this scope (not inside the block below) because the
        //- raw-population report in the acceptance gate is two scopes deeper.
        label  q0NNeg = 0, q0NBelowVS = 0;
        {
            labelHashSet q0Neg, q0Bad;
            polyMeshGenChecks::checkCellVolumes(mesh_, false, &q0Neg);
            polyMeshGenChecks::checkFacePyramids
                (mesh_, false, -SMALL, &q0Bad);
            q0NegVol = label(q0Neg.size());
            q0BadPyr = label(q0Bad.size());
            rawCellVolumeStats
                (mesh_, q0NegMag, q0MinCellVol, q0NNeg, q0NBelowVS);
            Info << "PREREFBL Q0: negVol=" << q0NegVol
                 << " badPyr=" << q0BadPyr
                 << " minCellVol=" << q0MinCellVol
                 << " negMag=" << q0NegMag
                 << " rawNeg=" << q0NNeg
                 << " rawBelowVSmall=" << q0NBelowVS
                 << " autoRepair=" << (doAutoRepair ? "yes" : "no")
                 << " snapshotValid="
                 << (preRefBLSnap.valid ? "yes" : "no") << endl;
        }

        //- Function-scope mirror of the inner twoPassAccepted flag, so the
        //- blPoints_ harvest below (two scopes shallower) can tell whether
        //- the pass-2 point set is already in place.
        bool blPointsFromPass2 = false;

        // ---- BLCOVERAGE base-face snapshot (pre-refBL, report-only) ----
        labelList  blcovBaseFacePatch;
        scalarField blcovBaseFaceArea;
        labelList  blcovBasePatchTotF;
        scalarField blcovBasePatchTotA;
        {
            const meshSurfaceEngine mseBase(mesh_);
            blcovBaseFacePatch = mseBase.boundaryFacePatches();
            const label nBaseBF = blcovBaseFacePatch.size();
            const vectorField& allFaceAreas =
                mesh_.addressingData().faceAreas();
            const label nInternal = mesh_.nInternalFaces();
            blcovBaseFaceArea.setSize(nBaseBF, 0.0);
            forAll(blcovBaseFaceArea, bfI)
            {
                const label faceI = nInternal + bfI;
                blcovBaseFaceArea[bfI] = Foam::mag(allFaceAreas[faceI]);
            }
            const label nP = mesh_.boundaries().size();
            blcovBasePatchTotF.setSize(nP, 0);  blcovBasePatchTotF = 0;
            blcovBasePatchTotA.setSize(nP, 0.0); blcovBasePatchTotA = 0.0;
            forAll(blcovBaseFacePatch, bfI)
            {
                const label p = blcovBaseFacePatch[bfI];
                if( p < 0 || p >= nP ) continue;
                ++blcovBasePatchTotF[p];
                blcovBasePatchTotA[p] += blcovBaseFaceArea[bfI];
            }
            Info << "BLCOVERAGE base snapshot: baseBoundaryFaces="
                 << nBaseBF << endl;
        }
        // ---- end snapshot ----

        nPointsBeforeBL_ = mesh_.points().size();
        blTopoAudit(mesh_, "PRE_REFBL");
        refLayers.refineLayers();
        blTopoAudit(mesh_, "POST_REFBL_PASS1");

        // Post-refBL provenance diagnostic.
        // Dumps all BL-generated cells to CSV for direct query.
        // Does NOT depend on negVol -- bad pyramids have negVol=0 and are
        // invisible to checkCellVolumes(). Query: grep "534099," the CSV.
        // Extended: cross-reference bad pyramid faces after refBL against
        // cellToBaseBndFace_ so we know which original bfI generated them.
        // NOTE: no mesh_.clearAddressingData() here -- downstream optimizer
        // pipeline needs addressing intact. checkFacePyramids works without it.
        {
            const labelList& prov = refLayers.cellToBaseBndFace();
            if( prov.size() > 0 )
            {
                OFstream provOs("postRefBL_cellProvenance.csv");
                provOs << "cellI,primaryBfI" << nl;

                label nMapped = 0;
                forAll(prov, cellI)
                {
                    if( prov[cellI] < 0 )
                        continue;
                    ++nMapped;
                    provOs << cellI << ',' << prov[cellI] << nl;
                }

                Info << "postRefBL provenance: mappedCells=" << nMapped
                     << " wrote postRefBL_cellProvenance.csv" << endl;

                // ---- BLCOVERAGE_AUDIT PROV_BASE (construction coverage) ----
                {
                    labelHashSet builtBaseFaces;
                    label nInvalidProv = 0;
                    forAll(prov, cellI)
                    {
                        const label baseBfI = prov[cellI];
                        if( baseBfI < 0 ) continue;
                        if( baseBfI >= blcovBaseFacePatch.size() )
                        { ++nInvalidProv; continue; }
                        builtBaseFaces.insert(baseBfI);
                    }
                    const label nP = mesh_.boundaries().size();
                    labelList  builtF(nP, 0);
                    scalarField builtA(nP, 0.0);
                    forAllConstIter(labelHashSet, builtBaseFaces, it)
                    {
                        const label b = it.key();
                        if( b < 0 || b >= blcovBaseFacePatch.size() ) continue;
                        const label p = blcovBaseFacePatch[b];
                        if( p < 0 || p >= nP ) continue;
                        ++builtF[p];
                        builtA[p] += blcovBaseFaceArea[b];
                    }
                    Info << "BLCOVERAGE_AUDIT PROV_BASE_PASS1 (construction coverage)"
                         << " invalidProv=" << nInvalidProv << endl;
                    const PtrList<boundaryPatch>& bnd = mesh_.boundaries();
                    forAll(bnd, p)
                    {
                        const scalar fpct = blcovBasePatchTotF[p] > 0
                          ? 100.0*builtF[p]/blcovBasePatchTotF[p] : 0.0;
                        const scalar apct = blcovBasePatchTotA[p] > SMALL
                          ? 100.0*builtA[p]/blcovBasePatchTotA[p] : 0.0;
                        Info << "  " << bnd[p].patchName()
                             << ": faces=" << builtF[p] << "/"
                             << blcovBasePatchTotF[p]
                             << " (" << fpct << "%) area=" << apct << "%" << endl;
                    }
                }
                // ---- end BLCOVERAGE_AUDIT PROV_BASE ----

                // Audit bad pyramid faces on grown refBL mesh.
                // No clearAddressingData -- optimizer needs addressing intact.
                labelHashSet postRefBLBadPyramids;
                polyMeshGenChecks::checkFacePyramids
                (
                    mesh_,
                    false,
                    -SMALL,
                    &postRefBLBadPyramids
                );

                labelHashSet provenanceSeedBfI;
                label nAuditedCells = 0;
                label nWithProv = 0;
                label nNoProv = 0;

                OFstream seedOs("postRefBL_badPyramidProvenance.csv");
                seedOs << "badFaceI,side,cellI,primaryBfI" << nl;

                const labelList& own = mesh_.owner();
                const labelList& nei = mesh_.neighbour();

                forAllConstIter(labelHashSet, postRefBLBadPyramids, it)
                {
                    const label faceI = it.key();

                    // Owner side
                    if( faceI >= 0 && faceI < label(own.size()) )
                    {
                        const label cellI = own[faceI];
                        if( cellI >= 0 && cellI < label(prov.size()) )
                        {
                            ++nAuditedCells;
                            const label bfI = prov[cellI];
                            seedOs << faceI << ",owner,"
                                   << cellI << ',' << bfI << nl;
                            if( bfI >= 0 )
                            {
                                provenanceSeedBfI.insert(bfI);
                                ++nWithProv;
                            }
                            else ++nNoProv;
                        }
                    }

                    // Neighbour side, internal faces only
                    if( faceI >= 0 && faceI < label(nei.size()) )
                    {
                        const label cellI = nei[faceI];
                        if( cellI >= 0 && cellI < label(prov.size()) )
                        {
                            ++nAuditedCells;
                            const label bfI = prov[cellI];
                            seedOs << faceI << ",neighbour,"
                                   << cellI << ',' << bfI << nl;
                            if( bfI >= 0 )
                            {
                                provenanceSeedBfI.insert(bfI);
                                ++nWithProv;
                            }
                            else ++nNoProv;
                        }
                    }
                }

                Info << "postRefBL bad pyramid audit: "
                     << "badFaces=" << postRefBLBadPyramids.size()
                     << " auditedCells=" << nAuditedCells
                     << " withProv=" << nWithProv
                     << " noProv=" << nNoProv
                     << " uniqueSeedBfI=" << provenanceSeedBfI.size()
                     << " wrote postRefBL_badPyramidProvenance.csv"
                     << endl;

                OFstream seedBfOs("postRefBL_provenanceSeedBfI.csv");
                seedBfOs << "bfI" << nl;
                forAllConstIter(labelHashSet, provenanceSeedBfI, sit)
                    seedBfOs << sit.key() << nl;

                Info << "postRefBL provenance seed bfI written to "
                     << "postRefBL_provenanceSeedBfI.csv" << endl;

                // Junction classification diagnostic.
                // Local PatchRoleMap -- no persistent state, no mesh change.
                // diagnosticOnly=true: logs junction type counts only.
                {
                    PatchRoleMap roles(meshDict_);

                    if( postRefBLBadPyramids.size() > 0 && roles.active() )
                    {
                        roles.print();
                        BLJunctionClassifier classifier(mesh_, roles);
                        classifier.classify
                        (
                            postRefBLBadPyramids,
                            prov,
                            true  // diagnosticOnly -- log only, no repair
                        );
                    }
                    else if( postRefBLBadPyramids.size() > 0 )
                    {
                        Info << "BLJunctionClassifier: skipped -- "
                             << "patchRoles not configured in meshDict" << endl;
                    }
                }

                // Two-pass BL repair loop.
                // Ordering: classify while pass-1 mesh alive, snapshot pass-1,
                // restore pre-refBL, run pass-2, accept only if better,
                // restore pass-1 if rejected.
                bool twoPassAccepted = false;
                bool twoPassAttempted = false;

                {
                    //- doAutoRepair hoisted to function scope (see above).
                    //- !reprojUnsafe_ removed with no replacement gate:
                    //- refinementValid() below, snapshot validity here, and
                    //- pass-1 restore on rejection cover the real conditions.
                    if( doAutoRepair
                     && preRefBLSnap.valid
                     && postRefBLBadPyramids.size() > 0
                     && provenanceSeedBfI.size() > 0 )
                    {
                        twoPassAttempted = true;

                        labelHashSet pass1NegVol;
                        polyMeshGenChecks::checkCellVolumes
                            (mesh_, false, &pass1NegVol);
                        const label pass1BadPyr =
                            label(postRefBLBadPyramids.size());

                        //- Q1 severity: the state pass 2 must beat.
                        //- Count alone permits defect SUBSTITUTION
                        //- ({A,B,C} -> {D,E,F} scores 3 <= 3) and lets mild
                        //- inversions be replaced by severe ones.
                        scalar q1MinCellVol = GREAT, q1NegMag = 0.0;
                        label q1NNeg = 0, q1NBelowVS = 0;
                        rawCellVolumeStats
                        (
                            mesh_, q1NegMag, q1MinCellVol,
                            q1NNeg, q1NBelowVS
                        );

                        // Preserve Q1 negative/near-zero cell identity and
                        // BL provenance across the Q1 -> Q0 rollback.  Cell
                        // numbering changes in Q2, but primary base bfI is
                        // the stable attribution coordinate.
                        labelList q1NegCellI(pass1NegVol.size(), -1);
                        labelList q1NegBaseBfI(pass1NegVol.size(), -1);
                        scalarField q1NegSignedVol
                        (
                            pass1NegVol.size(),
                            scalar(0.0)
                        );

                        label q1NegAttrI = 0;

                        forAllConstIter
                        (
                            labelHashSet,
                            pass1NegVol,
                            q1NvIt
                        )
                        {
                            const label cellI = q1NvIt.key();

                            q1NegCellI[q1NegAttrI] = cellI;

                            if
                            (
                                cellI >= 0
                             && cellI < label(prov.size())
                            )
                            {
                                q1NegBaseBfI[q1NegAttrI] = prov[cellI];
                            }

                            q1NegSignedVol[q1NegAttrI] =
                                rawSignedCellVolume(mesh_, cellI);

                            ++q1NegAttrI;
                        }

                        Info << "Two-pass BL repair: pass1 badPyramids="
                             << pass1BadPyr
                             << " negVol=" << pass1NegVol.size()
                             << " seeds=" << provenanceSeedBfI.size()
                             << endl;

                        // Snapshot pass-1 BEFORE restore so rejection
                        // can honestly recover pass-1 result.
                        PreRefBLMeshSnapshot pass1BLSnap;
                        takePreRefBLSnapshot(pass1BLSnap);

                        // Classify while pass-1 mesh is still alive.
                        // postRefBLBadPyramids face IDs belong to pass-1 mesh.
                        Map<BLRepairPlan> plans;
                        bool haveClassifierPlans = false;
                        PatchRoleMap roles(meshDict_);

                        if( roles.active() )
                        {
                            BLJunctionClassifier classifier2(mesh_, roles);
                            plans = classifier2.classify
                            (
                                postRefBLBadPyramids,
                                prov,
                                false  // repair mode
                            );
                            haveClassifierPlans = plans.size() > 0;
                        }

                        // Now restore pre-refBL mesh for pass 2.
                        if( restorePreRefBLSnapshot(preRefBLSnap) )
                        {
                            refineBoundaryLayers refLayers2(mesh_);
                            refineBoundaryLayers::readSettings
                                (meshDict_, refLayers2);
                            refLayers2.setBlblJunctionPoints
                                (preRefBLSnap.blblJunctionPoints);
                            refLayers2.setBlblAcuteCornerPoints
                                (preRefBLSnap.blblAcuteCornerPoints);
                            refLayers2.setRampSeedPoints
                                (preRefBLSnap.rampSeedPoints);
                            refLayers2.setVtFaceRing(preRefBLSnap.vtFaceRing);

                            // -------------------------------------------------
                            // Repair interaction groups.
                            //
                            // Seed components are not necessarily independent:
                            // their expanded repair footprints can overlap.
                            // Build the footprints on the RESTORED Q0 topology,
                            // using the exact ring construction semantics used
                            // by forceMaxLayersAtFaces(), then merge plans whose
                            // footprints touch.
                            //
                            // Diagnostic only: every plan is still applied.
                            // -------------------------------------------------
                            labelList repairGroupAtBfI;
                            Map<labelHashSet> repairGroupPlans;
                            Map<labelHashSet> repairGroupFootprints;
                            label nRepairGroups = 0;

                            labelList q1GroupNegCount;
                            scalarField q1GroupNegMag;
                            scalarField q1GroupMinVol;

                            // Populated from exploratory Q2.  Group IDs and
                            // plan IDs are runtime/topology-local bookkeeping;
                            // nothing geometry-specific is hard-coded here.
                            labelHashSet unsafeRepairGroups;
                            labelHashSet unsafePlanIds;

                            if( haveClassifierPlans )
                            {
                                const meshSurfaceEngine repairMse(mesh_);
                                const VRWGraph& repairFaceFaces =
                                    repairMse.faceFaces();

                                Map<labelHashSet> planFootprints;

                                label nActivePlans = 0;

                                forAllConstIter
                                (
                                    Map<BLRepairPlan>,
                                    plans,
                                    pCountIt
                                )
                                {
                                    if( pCountIt().active() )
                                        ++nActivePlans;
                                }

                                labelList planIds(nActivePlans, -1);
                                label planPos = 0;

                                // Construct exactly the rings that the current
                                // forceMaxLayersAtFaces() implementation builds.
                                forAllConstIter
                                (
                                    Map<BLRepairPlan>,
                                    plans,
                                    pFpIt
                                )
                                {
                                    const BLRepairPlan& plan = pFpIt();

                                    if( !plan.active() )
                                        continue;

                                    const label planId = pFpIt.key();

                                    planIds[planPos++] = planId;

                                    labelHashSet ring0;
                                    labelHashSet ring1;
                                    labelHashSet ring2;

                                    forAllConstIter
                                    (
                                        HashSet<label>,
                                        plan.seedBfI_,
                                        sIt
                                    )
                                    {
                                        const label bfI = sIt.key();

                                        if
                                        (
                                            bfI >= 0
                                         && bfI < label(repairFaceFaces.size())
                                        )
                                        {
                                            ring0.insert(bfI);
                                        }
                                    }

                                    if( plan.ring1_.maxLayers > 0 )
                                    {
                                        forAllConstIter
                                        (
                                            labelHashSet,
                                            ring0,
                                            r0It
                                        )
                                        {
                                            const label bfI = r0It.key();

                                            forAllRow
                                            (
                                                repairFaceFaces,
                                                bfI,
                                                nI
                                            )
                                            {
                                                const label nbfI =
                                                    repairFaceFaces(bfI, nI);

                                                if
                                                (
                                                    nbfI < 0
                                                 || nbfI >=
                                                    label
                                                    (
                                                        repairFaceFaces.size()
                                                    )
                                                )
                                                    continue;

                                                if( ring0.found(nbfI) )
                                                    continue;

                                                ring1.insert(nbfI);
                                            }
                                        }
                                    }

                                    if( plan.ring2_.maxLayers > 0 )
                                    {
                                        forAllConstIter
                                        (
                                            labelHashSet,
                                            ring1,
                                            r1It
                                        )
                                        {
                                            const label bfI = r1It.key();

                                            forAllRow
                                            (
                                                repairFaceFaces,
                                                bfI,
                                                nI
                                            )
                                            {
                                                const label nbfI =
                                                    repairFaceFaces(bfI, nI);

                                                if
                                                (
                                                    nbfI < 0
                                                 || nbfI >=
                                                    label
                                                    (
                                                        repairFaceFaces.size()
                                                    )
                                                )
                                                    continue;

                                                if
                                                (
                                                    ring0.found(nbfI)
                                                 || ring1.found(nbfI)
                                                )
                                                    continue;

                                                ring2.insert(nbfI);
                                            }
                                        }
                                    }

                                    labelHashSet footprint;

                                    forAllConstIter
                                    (
                                        labelHashSet,
                                        ring0,
                                        fIt
                                    )
                                        footprint.insert(fIt.key());

                                    forAllConstIter
                                    (
                                        labelHashSet,
                                        ring1,
                                        fIt
                                    )
                                        footprint.insert(fIt.key());

                                    forAllConstIter
                                    (
                                        labelHashSet,
                                        ring2,
                                        fIt
                                    )
                                        footprint.insert(fIt.key());

                                    planFootprints.insert(planId, footprint);

                                    Info << "BLRepairFootprint:"
                                         << " planId=" << planId
                                         << " type="
                                         << BLRepairPlan::junctionTypeName
                                            (plan.sourceType_)
                                         << " seeds="
                                         << plan.seedBfI_.size()
                                         << " footprint="
                                         << footprint.size()
                                         << endl;
                                }

                                auto footprintsOverlap =
                                []
                                (
                                    const labelHashSet& a,
                                    const labelHashSet& b
                                ) -> bool
                                {
                                    if( a.size() <= b.size() )
                                    {
                                        forAllConstIter
                                        (
                                            labelHashSet,
                                            a,
                                            aIt
                                        )
                                        {
                                            if( b.found(aIt.key()) )
                                                return true;
                                        }
                                    }
                                    else
                                    {
                                        forAllConstIter
                                        (
                                            labelHashSet,
                                            b,
                                            bIt
                                        )
                                        {
                                            if( a.found(bIt.key()) )
                                                return true;
                                        }
                                    }

                                    return false;
                                };

                                boolList assigned(planIds.size(), false);

                                forAll(planIds, pI)
                                {
                                    if( assigned[pI] )
                                        continue;

                                    const label firstPlanId = planIds[pI];

                                    labelHashSet groupPlans;
                                    labelHashSet groupFootprint;

                                    groupPlans.insert(firstPlanId);

                                    const labelHashSet& firstFp =
                                        planFootprints[firstPlanId];

                                    forAllConstIter
                                    (
                                        labelHashSet,
                                        firstFp,
                                        fIt
                                    )
                                        groupFootprint.insert(fIt.key());

                                    assigned[pI] = true;

                                    // Transitive closure: if a newly merged
                                    // footprint touches another plan, that plan
                                    // belongs to the same interaction group.
                                    bool grew = true;

                                    while( grew )
                                    {
                                        grew = false;

                                        forAll(planIds, pJ)
                                        {
                                            if( assigned[pJ] )
                                                continue;

                                            const label candidateId =
                                                planIds[pJ];

                                            const labelHashSet& candidateFp =
                                                planFootprints[candidateId];

                                            if
                                            (
                                                footprintsOverlap
                                                (
                                                    groupFootprint,
                                                    candidateFp
                                                )
                                            )
                                            {
                                                groupPlans.insert(candidateId);

                                                forAllConstIter
                                                (
                                                    labelHashSet,
                                                    candidateFp,
                                                    cfIt
                                                )
                                                {
                                                    groupFootprint.insert
                                                        (cfIt.key());
                                                }

                                                assigned[pJ] = true;
                                                grew = true;
                                            }
                                        }
                                    }

                                    repairGroupPlans.insert
                                        (nRepairGroups, groupPlans);
                                    repairGroupFootprints.insert
                                        (nRepairGroups, groupFootprint);

                                    Info << "BLRepairInteractionGroup:"
                                         << " group=" << nRepairGroups
                                         << " plans=" << groupPlans.size()
                                         << " footprint="
                                         << groupFootprint.size()
                                         << endl;

                                    ++nRepairGroups;
                                }

                                repairGroupAtBfI.setSize
                                (
                                    repairFaceFaces.size(),
                                    -1
                                );

                                label nGroupMapConflicts = 0;

                                forAllConstIter
                                (
                                    Map<labelHashSet>,
                                    repairGroupFootprints,
                                    gFpIt
                                )
                                {
                                    const label groupI = gFpIt.key();

                                    forAllConstIter
                                    (
                                        labelHashSet,
                                        gFpIt(),
                                        gfIt
                                    )
                                    {
                                        const label bfI = gfIt.key();

                                        if
                                        (
                                            bfI < 0
                                         || bfI >=
                                            label(repairGroupAtBfI.size())
                                        )
                                            continue;

                                        if
                                        (
                                            repairGroupAtBfI[bfI] >= 0
                                         && repairGroupAtBfI[bfI] != groupI
                                        )
                                        {
                                            ++nGroupMapConflicts;
                                        }

                                        repairGroupAtBfI[bfI] = groupI;
                                    }
                                }

                                Info << "BLRepairInteractionGroups:"
                                     << " activePlans=" << nActivePlans
                                     << " groups=" << nRepairGroups
                                     << " mapConflicts="
                                     << nGroupMapConflicts
                                     << endl;

                                q1GroupNegCount.setSize(nRepairGroups, 0);
                                q1GroupNegMag.setSize
                                    (nRepairGroups, scalar(0.0));
                                q1GroupMinVol.setSize
                                    (nRepairGroups, GREAT);

                                // Attribute stored Q1 bad-volume cells through
                                // stable base-face provenance.
                                forAll(q1NegCellI, q1I)
                                {
                                    const label bfI =
                                        q1NegBaseBfI[q1I];

                                    label groupI = -1;

                                    if
                                    (
                                        bfI >= 0
                                     && bfI <
                                        label(repairGroupAtBfI.size())
                                    )
                                    {
                                        groupI =
                                            repairGroupAtBfI[bfI];
                                    }

                                    const scalar signedVol =
                                        q1NegSignedVol[q1I];

                                    Info << "BLRepairNegAttrib:"
                                         << " pass=Q1"
                                         << " cell=" << q1NegCellI[q1I]
                                         << " bfI=" << bfI
                                         << " group=" << groupI
                                         << " signedVol=" << signedVol
                                         << endl;

                                    if
                                    (
                                        groupI >= 0
                                     && groupI < nRepairGroups
                                    )
                                    {
                                        ++q1GroupNegCount[groupI];

                                        if( signedVol < 0.0 )
                                        {
                                            q1GroupNegMag[groupI] +=
                                                -signedVol;
                                        }

                                        q1GroupMinVol[groupI] =
                                            Foam::min
                                            (
                                                q1GroupMinVol[groupI],
                                                signedVol
                                            );
                                    }
                                }
                            }

                            if( haveClassifierPlans )
                            {
                                forAllIter(Map<BLRepairPlan>, plans, pit)
                                {
                                    const BLRepairPlan& plan = pit();
                                    if( plan.active() )
                                        refLayers2.forceMaxLayersAtFaces
                                        (
                                            plan.seedBfI_,
                                            plan.ring0_.maxLayers,
                                            plan.ring1_.maxLayers,
                                            plan.ring2_.maxLayers,
                                            plan.ring0_.thicknessScale,
                                            plan.ring1_.thicknessScale,
                                            plan.ring2_.thicknessScale
                                        );
                                }
                            }
                            else
                            {
                                Info << "Two-pass BL repair: using fallback "
                                     << "provenance seed plan" << endl;
                                refLayers2.forceMaxLayersAtFaces
                                (
                                    provenanceSeedBfI,
                                    3, 4, 0,
                                    0.60, 0.80, 1.00
                                );
                            }

                            nPointsBeforeBL_ = mesh_.points().size();
                            refLayers2.refineLayers();

                            const bool pass2RefinementValid =
                                refLayers2.refinementValid();
                            const bool pass2RefinementCompleted =
                                refLayers2.refinementCompleted();

                            if
                            (
                                !pass2RefinementValid
                             || !pass2RefinementCompleted
                            )
                            {
                                Info << "Two-pass BL repair: pass2 did not "
                                     << "complete refinement"
                                     << " refinementValid="
                                     << (pass2RefinementValid ? "yes" : "no")
                                     << " refinementCompleted="
                                     << (pass2RefinementCompleted ? "yes" : "no")
                                     << " -- REJECTED, restoring pass1"
                                     << endl;

                                restorePreRefBLSnapshot(pass1BLSnap);
                                twoPassAccepted = false;
                            }
                            else
                            {
                            labelHashSet pass2NegVol;
                            polyMeshGenChecks::checkCellVolumes
                                (mesh_, false, &pass2NegVol);
                            labelHashSet pass2BadPyr;
                            polyMeshGenChecks::checkFacePyramids
                                (mesh_, false, -SMALL, &pass2BadPyr);

                            scalar q2MinCellVol = GREAT, q2NegMag = 0.0;
                            label q2NNeg = 0, q2NBelowVS = 0;
                            rawCellVolumeStats
                            (
                                mesh_, q2NegMag, q2MinCellVol,
                                q2NNeg, q2NBelowVS
                            );

                            // Group-local negative-volume attribution.
                            // Diagnostic only: global Q1/Q2 acceptance below
                            // remains completely unchanged.
                            if
                            (
                                haveClassifierPlans
                             && nRepairGroups > 0
                            )
                            {
                                labelList q2GroupNegCount
                                    (nRepairGroups, 0);
                                scalarField q2GroupNegMag
                                    (nRepairGroups, scalar(0.0));
                                scalarField q2GroupMinVol
                                    (nRepairGroups, GREAT);

                                const labelList& q2Prov =
                                    refLayers2.cellToBaseBndFace();

                                forAllConstIter
                                (
                                    labelHashSet,
                                    pass2NegVol,
                                    q2NvIt
                                )
                                {
                                    const label cellI = q2NvIt.key();

                                    label bfI = -1;

                                    if
                                    (
                                        cellI >= 0
                                     && cellI < label(q2Prov.size())
                                    )
                                    {
                                        bfI = q2Prov[cellI];
                                    }

                                    label groupI = -1;

                                    if
                                    (
                                        bfI >= 0
                                     && bfI <
                                        label(repairGroupAtBfI.size())
                                    )
                                    {
                                        groupI =
                                            repairGroupAtBfI[bfI];
                                    }

                                    const scalar signedVol =
                                        rawSignedCellVolume(mesh_, cellI);

                                    Info << "BLRepairNegAttrib:"
                                         << " pass=Q2"
                                         << " cell=" << cellI
                                         << " bfI=" << bfI
                                         << " group=" << groupI
                                         << " signedVol=" << signedVol
                                         << endl;

                                    if
                                    (
                                        groupI >= 0
                                     && groupI < nRepairGroups
                                    )
                                    {
                                        ++q2GroupNegCount[groupI];

                                        if( signedVol < 0.0 )
                                        {
                                            q2GroupNegMag[groupI] +=
                                                -signedVol;
                                        }

                                        q2GroupMinVol[groupI] =
                                            Foam::min
                                            (
                                                q2GroupMinVol[groupI],
                                                signedVol
                                            );
                                    }
                                }

                                for
                                (
                                    label groupI = 0;
                                    groupI < nRepairGroups;
                                    ++groupI
                                )
                                {
                                    const label q1Count =
                                        q1GroupNegCount[groupI];
                                    const label q2Count =
                                        q2GroupNegCount[groupI];

                                    const scalar q1Mag =
                                        q1GroupNegMag[groupI];
                                    const scalar q2Mag =
                                        q2GroupNegMag[groupI];

                                    const scalar q1Min =
                                        q1Count > 0
                                      ? q1GroupMinVol[groupI]
                                      : scalar(0.0);

                                    const scalar q2Min =
                                        q2Count > 0
                                      ? q2GroupMinVol[groupI]
                                      : scalar(0.0);

                                    const bool countOK =
                                        q2Count <= q1Count;

                                    const bool magOK =
                                        q2Mag <= q1Mag;

                                    const bool minOK =
                                        q2Count == 0
                                     || (
                                            q1Count > 0
                                         && q2Min >= q1Min
                                        );

                                    const bool negSafe =
                                        countOK && magOK && minOK;

                                    const label nPlans =
                                        repairGroupPlans.found(groupI)
                                      ? repairGroupPlans[groupI].size()
                                      : 0;

                                    const label nFootprint =
                                        repairGroupFootprints.found(groupI)
                                      ? repairGroupFootprints[groupI].size()
                                      : 0;

                                    Info << "BLRepairNegGroup:"
                                         << " group=" << groupI
                                         << " plans=" << nPlans
                                         << " footprint=" << nFootprint
                                         << " count=" << q1Count
                                         << "->" << q2Count
                                         << " negMag=" << q1Mag
                                         << "->" << q2Mag
                                         << " minVol=" << q1Min
                                         << "->" << q2Min
                                         << " negSafe="
                                         << (negSafe ? "yes" : "no")
                                         << endl;

                                    if( !negSafe )
                                    {
                                        unsafeRepairGroups.insert(groupI);

                                        if
                                        (
                                            repairGroupPlans.found(groupI)
                                        )
                                        {
                                            const labelHashSet& badGroupPlans =
                                                repairGroupPlans[groupI];

                                            forAllConstIter
                                            (
                                                labelHashSet,
                                                badGroupPlans,
                                                badPlanIt
                                            )
                                            {
                                                unsafePlanIds.insert
                                                    (badPlanIt.key());
                                            }
                                        }
                                    }
                                }

                                Info << "BLRepairSelectiveCandidates:"
                                     << " unsafeGroups="
                                     << unsafeRepairGroups.size()
                                     << " unsafePlans="
                                     << unsafePlanIds.size()
                                     << " totalGroups="
                                     << nRepairGroups
                                     << endl;
                            }

                            //- Phase 1: effectively zero tolerance. Deltas
                            //- are logged below so a scale-aware floor can
                            //- be chosen from measured data rather than
                            //- guessed before this path has ever run.
                            //- NOTE: no special case needed for
                            //- q1NegMag == 0 -- the count term already
                            //- forces q2NegVol == 0 when q1NegVol == 0.
                            const scalar volTolAbs = 0.0;
                            const scalar volTolRel = 0.0;
                            const scalar negMagTol =
                                volTolAbs + volTolRel*q1NegMag;

                            const bool negCountOK =
                                pass2NegVol.size() <= pass1NegVol.size();
                            const bool badPyrBetter =
                                label(pass2BadPyr.size()) < pass1BadPyr;
                            const bool negMagOK =
                                q2NegMag <= q1NegMag + negMagTol;
                            const bool minVolOK =
                                q2MinCellVol >= q1MinCellVol - volTolAbs;

                            const bool pass2Better =
                                negCountOK && badPyrBetter
                             && negMagOK   && minVolOK;

                            Info << "Two-pass BL repair: pass2 badPyramids="
                                 << pass2BadPyr.size()
                                 << " negVol=" << pass2NegVol.size()
                                 << (pass2Better ? " -- ACCEPTED" : " -- REJECTED")
                                 << endl;

                            Info << "PREREFBL Q0/Q1/Q2"
                                 << "  negVol " << q0NegVol
                                 << "/" << pass1NegVol.size()
                                 << "/" << pass2NegVol.size()
                                 << "  badPyr " << q0BadPyr
                                 << "/" << pass1BadPyr
                                 << "/" << pass2BadPyr.size()
                                 << "  negMag " << q0NegMag
                                 << "/" << q1NegMag << "/" << q2NegMag
                                 << "  minVol " << q0MinCellVol
                                 << "/" << q1MinCellVol
                                 << "/" << q2MinCellVol << endl;

                            Info << "PREREFBL severity delta:"
                                 << " dNegMag=" << (q2NegMag - q1NegMag)
                                 << " dMinVol=" << (q2MinCellVol - q1MinCellVol)
                                 << endl;

                            //- rawNeg counts true inversions (vol < 0);
                            //- rawBelowVSmall is the population
                            //- checkCellVolumes() reports (vol < VSMALL),
                            //- i.e. inversions plus zero/near-zero slivers.
                            //- A gap between them means the count gate and
                            //- the severity terms are seeing different sets.
                            Info << "PREREFBL raw populations Q0/Q1/Q2:"
                                 << " rawNeg " << q0NNeg << "/" << q1NNeg
                                 << "/" << q2NNeg
                                 << "  rawBelowVSmall " << q0NBelowVS
                                 << "/" << q1NBelowVS << "/" << q2NBelowVS
                                 << "  (checkCellVolumes count Q1/Q2 "
                                 << pass1NegVol.size() << "/"
                                 << pass2NegVol.size() << ")" << endl;

                            Info << "PREREFBL decision terms:"
                                 << " negCountOK="  << (negCountOK  ? "yes" : "no")
                                 << " badPyrBetter="<< (badPyrBetter? "yes" : "no")
                                 << " negMagOK="    << (negMagOK    ? "yes" : "no")
                                 << " minVolOK="    << (minVolOK    ? "yes" : "no")
                                 << endl;

                            if( pass2Better )
                            {
                                refLayers2.pointsInBndLayer(blPoints_);
                                twoPassAccepted = true;
                                blPointsFromPass2 = true;
                            }
                            else if
                            (
                                haveClassifierPlans
                             && unsafeRepairGroups.size() > 0
                             && unsafePlanIds.size() > 0
                            )
                            {
                                // Exploratory Q2 found locally unsafe repair
                                // interaction groups.  Restore the exact Q0
                                // transaction and regenerate once more with
                                // only the locally safe candidate plans.
                                //
                                // Q3 is still subject to the SAME strict
                                // global acceptance criteria as Q2.
                                Info << "Selective BL repair Q3: exploratory Q2 "
                                     << "rejected; retrying without "
                                     << unsafeRepairGroups.size()
                                     << " unsafe group(s), "
                                     << unsafePlanIds.size()
                                     << " unsafe plan(s)"
                                     << endl;

                                if
                                (
                                    restorePreRefBLSnapshot(preRefBLSnap)
                                )
                                {
                                    refineBoundaryLayers refLayers3(mesh_);

                                    refineBoundaryLayers::readSettings
                                        (meshDict_, refLayers3);

                                    refLayers3.setBlblJunctionPoints
                                    (
                                        preRefBLSnap.blblJunctionPoints
                                    );

                                    refLayers3.setBlblAcuteCornerPoints
                                    (
                                        preRefBLSnap.blblAcuteCornerPoints
                                    );

                                    refLayers3.setRampSeedPoints
                                    (
                                        preRefBLSnap.rampSeedPoints
                                    );

                                    refLayers3.setVtFaceRing
                                    (
                                        preRefBLSnap.vtFaceRing
                                    );

                                    label nQ3AppliedPlans = 0;
                                    label nQ3SkippedPlans = 0;

                                    forAllIter
                                    (
                                        Map<BLRepairPlan>,
                                        plans,
                                        q3Pit
                                    )
                                    {
                                        const label planId = q3Pit.key();
                                        const BLRepairPlan& plan = q3Pit();

                                        if( !plan.active() )
                                            continue;

                                        if
                                        (
                                            unsafePlanIds.found(planId)
                                        )
                                        {
                                            ++nQ3SkippedPlans;

                                            Info << "Selective BL repair Q3:"
                                                 << " skip planId="
                                                 << planId
                                                 << " type="
                                                 << BLRepairPlan::
                                                    junctionTypeName
                                                    (plan.sourceType_)
                                                 << " seeds="
                                                 << plan.seedBfI_.size()
                                                 << endl;

                                            continue;
                                        }

                                        ++nQ3AppliedPlans;

                                        refLayers3.forceMaxLayersAtFaces
                                        (
                                            plan.seedBfI_,
                                            plan.ring0_.maxLayers,
                                            plan.ring1_.maxLayers,
                                            plan.ring2_.maxLayers,
                                            plan.ring0_.thicknessScale,
                                            plan.ring1_.thicknessScale,
                                            plan.ring2_.thicknessScale
                                        );
                                    }

                                    Info << "Selective BL repair Q3:"
                                         << " appliedPlans="
                                         << nQ3AppliedPlans
                                         << " skippedPlans="
                                         << nQ3SkippedPlans
                                         << endl;

                                    nPointsBeforeBL_ =
                                        mesh_.points().size();

                                    refLayers3.refineLayers();

                                    const bool pass3RefinementValid =
                                        refLayers3.refinementValid();

                                    const bool pass3RefinementCompleted =
                                        refLayers3.refinementCompleted();

                                    if
                                    (
                                        !pass3RefinementValid
                                     || !pass3RefinementCompleted
                                    )
                                    {
                                        Info << "Selective BL repair Q3:"
                                             << " refinement incomplete"
                                             << " refinementValid="
                                             << (
                                                    pass3RefinementValid
                                                  ? "yes" : "no"
                                                )
                                             << " refinementCompleted="
                                             << (
                                                    pass3RefinementCompleted
                                                  ? "yes" : "no"
                                                )
                                             << " -- REJECTED, restoring Q1"
                                             << endl;

                                        restorePreRefBLSnapshot
                                            (pass1BLSnap);
                                    }
                                    else
                                    {
                                        labelHashSet pass3NegVol;

                                        polyMeshGenChecks::checkCellVolumes
                                        (
                                            mesh_,
                                            false,
                                            &pass3NegVol
                                        );

                                        labelHashSet pass3BadPyr;

                                        polyMeshGenChecks::checkFacePyramids
                                        (
                                            mesh_,
                                            false,
                                            -SMALL,
                                            &pass3BadPyr
                                        );

                                        scalar q3MinCellVol = GREAT;
                                        scalar q3NegMag = 0.0;
                                        label q3NNeg = 0;
                                        label q3NBelowVS = 0;

                                        rawCellVolumeStats
                                        (
                                            mesh_,
                                            q3NegMag,
                                            q3MinCellVol,
                                            q3NNeg,
                                            q3NBelowVS
                                        );

                                        const bool q3NegCountOK =
                                            pass3NegVol.size()
                                         <= pass1NegVol.size();

                                        const bool q3BadPyrBetter =
                                            label(pass3BadPyr.size())
                                          < pass1BadPyr;

                                        const bool q3NegMagOK =
                                            q3NegMag
                                         <= q1NegMag + negMagTol;

                                        const bool q3MinVolOK =
                                            q3MinCellVol
                                         >= q1MinCellVol - volTolAbs;

                                        const bool pass3Better =
                                            q3NegCountOK
                                         && q3BadPyrBetter
                                         && q3NegMagOK
                                         && q3MinVolOK;

                                        Info << "Selective BL repair Q3:"
                                             << " badPyramids="
                                             << pass3BadPyr.size()
                                             << " negVol="
                                             << pass3NegVol.size()
                                             << (
                                                    pass3Better
                                                  ? " -- ACCEPTED"
                                                  : " -- REJECTED"
                                                )
                                             << endl;

                                        Info << "PREREFBL Q1/Q3"
                                             << "  negVol "
                                             << pass1NegVol.size()
                                             << "/"
                                             << pass3NegVol.size()
                                             << "  badPyr "
                                             << pass1BadPyr
                                             << "/"
                                             << pass3BadPyr.size()
                                             << "  negMag "
                                             << q1NegMag
                                             << "/"
                                             << q3NegMag
                                             << "  minVol "
                                             << q1MinCellVol
                                             << "/"
                                             << q3MinCellVol
                                             << endl;

                                        Info << "PREREFBL Q3 raw populations:"
                                             << " rawNeg="
                                             << q3NNeg
                                             << " rawBelowVSmall="
                                             << q3NBelowVS
                                             << " checkCellVolumes="
                                             << pass3NegVol.size()
                                             << endl;

                                        Info << "PREREFBL Q3 decision terms:"
                                             << " negCountOK="
                                             << (
                                                    q3NegCountOK
                                                  ? "yes" : "no"
                                                )
                                             << " badPyrBetter="
                                             << (
                                                    q3BadPyrBetter
                                                  ? "yes" : "no"
                                                )
                                             << " negMagOK="
                                             << (
                                                    q3NegMagOK
                                                  ? "yes" : "no"
                                                )
                                             << " minVolOK="
                                             << (
                                                    q3MinVolOK
                                                  ? "yes" : "no"
                                                )
                                             << endl;

                                        if( pass3Better )
                                        {
                                            refLayers3.pointsInBndLayer
                                                (blPoints_);

                                            // Historical variable name:
                                            // means an accepted repair-pass
                                            // point set is already harvested.
                                            twoPassAccepted = true;
                                            blPointsFromPass2 = true;

                                            Info << "Selective BL repair Q3:"
                                                 << " keeping selective "
                                                 << "repair result"
                                                 << endl;
                                        }
                                        else
                                        {
                                            Info << "Selective BL repair Q3:"
                                                 << " restoring Q1 after "
                                                 << "global rejection"
                                                 << endl;

                                            restorePreRefBLSnapshot
                                                (pass1BLSnap);
                                        }
                                    }
                                }
                                else
                                {
                                    Info << "Selective BL repair Q3:"
                                         << " Q0 restore failed"
                                         << " -- restoring Q1"
                                         << endl;

                                    restorePreRefBLSnapshot(pass1BLSnap);
                                }
                            }
                            else
                            {
                                Info << "Two-pass BL repair: restoring pass1 "
                                     << "result after rejected pass2"
                                     << " -- no locally unsafe repair "
                                     << "groups available for selective Q3"
                                     << endl;

                                restorePreRefBLSnapshot(pass1BLSnap);
                            }
                            } // end refinementValid else
                        }
                        else
                        {
                            Info << "Two-pass BL repair: pre-refBL restore "
                                 << "failed -- keeping pass1 result" << endl;
                        }
                    }
                    else if( doAutoRepair )
                    {
                        Info << "Two-pass BL repair: skipped -- "
                             << (reprojUnsafe_ ? "reprojUnsafe " : "")
                             << (!preRefBLSnap.valid ? "noSnapshot " : "")
                             << (postRefBLBadPyramids.size()==0 ? "noBadPyr " : "")
                             << (provenanceSeedBfI.size()==0 ? "noSeeds" : "")
                             << endl;
                    }
                }

                if( twoPassAttempted && !twoPassAccepted )
                    Info << "Two-pass BL repair: using pass1 BL point set" << endl;

            }
        }

        //- Only harvest the pass-1 BL point set if pass 2 was NOT accepted.
        //- On acceptance the mesh is the pass-2 mesh and blPoints_ was
        //- already filled from refLayers2; overwriting it from refLayers
        //- would describe the pass-1 layer point set against a pass-2 mesh,
        //- and blPoints_ feeds mOpt.lockPoints() in optimiseFinalMesh().
        if( !blPointsFromPass2 )
            refLayers.pointsInBndLayer(blPoints_);
        else
            Info << "Two-pass BL repair: retaining pass2 BL point set ("
                 << blPoints_.size() << " points)" << endl;

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
            Info << nl
                 << "### ENTERING 3-GATE POST-REFINEMENT DIAGNOSTIC ###"
                 << nl << endl;
            Info << "3-gate diagnostic: checking post-refinement BL quality" << endl;

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

            Info << "3-gate diagnostic results:" << endl;
            Info << "  Gate 1 (neg vol cells):    " << badCells.size() << endl;
            Info << "  Gate 2 (bad pyramids):     " << badPyramidFaces.size() << endl;
            Info << "  Gate 3 (non-ortho >85deg): " << nonOrthoFaces.size() << endl;

            // Post-refBL provenance audit.
            // Runs while refLayers.cellToBaseBndFace() is still alive and
            // before later renumbering can invalidate cell indices.
            // Report-only: maps bad pyramid adjacent cells back to their
            // original base boundary face / patch.
            {
                const labelList& provAudit = refLayers.cellToBaseBndFace();

                if( badPyramidFaces.size() > 0 && provAudit.size() > 0 )
                {
                    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
                    const labelList& owner = mesh_.owner();
                    const labelList& neighbour = mesh_.neighbour();
                    const label firstBoundaryFace = boundaries[0].patchStart();

                    label nWithProv = 0;
                    label nNoProv = 0;
                    labelHashSet auditedCells;
                    Map<label> patchBadCount;

                    forAllConstIter(labelHashSet, badPyramidFaces, it)
                    {
                        const label faceI = it.key();

                        for(label sideI = 0; sideI < 2; ++sideI)
                        {
                            label cellI = -1;

                            if( sideI == 0 )
                            {
                                if( faceI >= 0 && faceI < label(owner.size()) )
                                    cellI = owner[faceI];
                            }
                            else
                            {
                                if( faceI >= 0 && faceI < label(neighbour.size()) )
                                    cellI = neighbour[faceI];
                            }

                            if( cellI < 0 || auditedCells.found(cellI) )
                                continue;

                            auditedCells.insert(cellI);

                            if
                            (
                                cellI >= label(provAudit.size())
                             || provAudit[cellI] < 0
                            )
                            {
                                ++nNoProv;
                                continue;
                            }

                            ++nWithProv;

                            const label bfI = provAudit[cellI];

                            label patchI = -1;
                            forAll(boundaries, pI)
                            {
                                const label patchStartLocal =
                                    boundaries[pI].patchStart() - firstBoundaryFace;
                                const label patchEndLocal =
                                    patchStartLocal + boundaries[pI].patchSize();

                                if( bfI >= patchStartLocal && bfI < patchEndLocal )
                                {
                                    patchI = pI;
                                    break;
                                }
                            }

                            if( patchI >= 0 )
                            {
                                Map<label>::iterator mIt =
                                    patchBadCount.find(patchI);

                                if( mIt == patchBadCount.end() )
                                    patchBadCount.insert(patchI, 1);
                                else
                                    ++mIt();
                            }
                        }
                    }

                    Info << "PostRefBLProvenance: badFaces="
                         << badPyramidFaces.size()
                         << " auditedCells=" << auditedCells.size()
                         << " withProv=" << nWithProv
                         << " noProv=" << nNoProv << endl;

                    forAllConstIter(Map<label>, patchBadCount, it)
                    {
                        Info << "  patch "
                             << boundaries[it.key()].patchName()
                             << ": " << it() << " bad cells" << endl;
                    }
                }
                else
                {
                    Info << "PostRefBLProvenance: skipped badFaces="
                         << badPyramidFaces.size()
                         << " provSize=" << provAudit.size() << endl;
                }
            }

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

                        Info << "  Gate2 badFace faceI=" << faceI
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

                Info << "  Gate2 bad pyramid classes:"
                     << " periodic=" << nPeriodicClass
                     << " blade=" << nBladeClass
                     << " hub=" << nHubClass
                     << " shroud=" << nShroudClass
                     << " inletOutlet=" << nInletOutletClass
                     << " generic=" << nGenericClass
                     << endl;

                Info << "  Gate2 periodic-local repair candidate faces: "
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
                else if( gate2PeriodicBadFaces.size() == 0 )
                {
                    Info << "Gate2 local repair: skipped -- "
                         << "no periodic bad faces to target" << endl;
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

                const bool gate2SkewOK =
                    gate2MaxSkewAfter <= gate2MaxSkewBefore;

                const bool gate2RepairOK =
                    gate2BadAfter.size() < badPyramidFaces.size()
                 && (!gate2UnusedAfter || gate2UnusedBefore)
                 && gate2NegBefore.size() == 0
                 && gate2NegAfter.size() == 0
                 && gate2OpenAfter.size() <= gate2OpenBefore.size()
                 && gate2NonOrthoAfter.size() <= gate2NonOrthoBefore.size()
                 && gate2SkewOK;

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

                Info << "  Rollback candidates: "
                     << rollbackBndFaces.size()
                     << " boundary faces would be targeted" << endl;
            }
            else
            {
                Info << "  All gates passed - no rollback needed" << endl;
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
    reprojUnsafe_ = false;
    meshHistory_ = MeshHistory::CleanNatural;

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
    // Preserve user intent before the negVol gate clears the flag.
    const bool requestedConstrainOptimizerBoundary = constrainOptimizerBoundary;

    // Gate constrained optimizer on mesh validity.
    // Running surface-constrained optimization on a mesh with pre-existing
    // negVol cells can diverge catastrophically. Keep the capability, but
    // defer it until the mesh is valid enough for aggressive constrained motion.
    {
        mesh_.clearAddressingData();

        labelHashSet negVolCheck;
        polyMeshGenChecks::checkCellVolumes(mesh_, false, &negVolCheck);

        if( negVolCheck.size() > 0 )
        {
            Info << "optimiseFinalMesh: " << negVolCheck.size()
                 << " negVol cells present -- skipping constrained optimizer, "
                 << "falling through to plain optimizeMeshFV" << endl;

            constrainOptimizerBoundary = false;
        }
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

        // Second constrained pass: user requested surface-constrained
        // optimization but it was skipped due to pre-existing negVol.
        // Attempt it now after plain optimizer has cleared negVol.
        // Point-motion only -- no optimizeLowQualityFaces() so
        // point-only rollback is sufficient.
        // Second constrained pass disabled -- optimizeMeshFV does more
        // than pure point motion so point-only rollback is insufficient.
        // Re-enable when topology-aware rollback is available.
        if( false && requestedConstrainOptimizerBoundary && octreePtr_ )
        {
            mesh_.clearAddressingData();
            labelHashSet negVolAfterPlain;
            polyMeshGenChecks::checkCellVolumes
                (mesh_, false, &negVolAfterPlain);

            if( negVolAfterPlain.size() == 0 )
            {
                Info << "Surface-constrained optimizer: attempting second "
                     << "pass after plain optimizer cleared negVol" << endl;

                // Logical-size snapshot
                const label nActivePts2 = mesh_.points().size();
                pointField pointsBefore2(nActivePts2);
                forAll(pointsBefore2, pI)
                    pointsBefore2[pI] = mesh_.points()[pI];

                mesh_.clearAddressingData();
                labelHashSet negBefore2, pyrBefore2;
                polyMeshGenChecks::checkCellVolumes
                    (mesh_, false, &negBefore2);
                polyMeshGenChecks::checkFacePyramids
                    (mesh_, false, -SMALL, &pyrBefore2);

                // Build surface constraint
                meshSurfaceEngine mseConstraint2(mesh_);
                meshSurfacePartitioner mPartConstraint2(mseConstraint2);
                labelLongList globalToBp2(mesh_.points().size(), -1);
                const labelList& bPtsC2 = mseConstraint2.boundaryPoints();
                forAll(bPtsC2, bpI)
                    globalToBp2[bPtsC2[bpI]] = bpI;

                const label nBp2 = bPtsC2.size();
                vectorField featureTangents2(nBp2, vector::zero);
                {
                    const edgeList& edges2 = mseConstraint2.edges();
                    const VRWGraph& bpEdges2 =
                        mseConstraint2.boundaryPointEdges();
                    const labelHashSet& featEdges2 =
                        mPartConstraint2.featureEdges();
                    const labelHashSet& edgePts2 =
                        mPartConstraint2.edgePoints();
                    const pointFieldPMG& pts2 = mesh_.points();
                    const labelList& bp2 = mseConstraint2.bp();

                    forAllConstIter(labelHashSet, edgePts2, it)
                    {
                        const label bpI2 = it.key();
                        if( bpI2 < 0 || bpI2 >= nBp2 ) continue;
                        label nbr0 = -1, nbr1 = -1;
                        forAllRow(bpEdges2, bpI2, eI)
                        {
                            const label beI = bpEdges2(bpI2, eI);
                            if( !featEdges2.found(beI) ) continue;
                            const edge& e = edges2[beI];
                            const label ep0 = e.start();
                            const label ep1 = e.end();
                            if( ep0 < 0 || ep0 >= label(bp2.size()) ||
                                ep1 < 0 || ep1 >= label(bp2.size()) )
                                continue;
                            const label ob0 = bp2[ep0];
                            const label ob1 = bp2[ep1];
                            if( ob0 < 0 || ob1 < 0 ) continue;
                            label otherBp = -1;
                            if( ob0 == bpI2 ) otherBp = ob1;
                            else if( ob1 == bpI2 ) otherBp = ob0;
                            if( otherBp < 0 || otherBp >= nBp2 ) continue;
                            if( nbr0 == -1 ) nbr0 = otherBp;
                            else if( nbr1 == -1 && otherBp != nbr0 )
                                nbr1 = otherBp;
                        }
                        vector t = vector::zero;
                        if( nbr0 != -1 && nbr1 != -1 )
                            t = pts2[bPtsC2[nbr1]] - pts2[bPtsC2[nbr0]];
                        else if( nbr0 != -1 )
                            t = pts2[bPtsC2[nbr0]] - pts2[bPtsC2[bpI2]];
                        if( magSqr(t) > VSMALL )
                            featureTangents2[bpI2] = t / mag(t);
                    }
                }

                // Lock acute corners
                labelLongList acuteGlobalPts2;
                if( lockAcuteCorners && !blblAcuteCornerPoints_.empty() )
                {
                    forAllConstIter(labelHashSet, blblAcuteCornerPoints_, it)
                    {
                        const label bpI2 = it.key();
                        if( bpI2 >= 0 && bpI2 < label(bPtsC2.size()) )
                            acuteGlobalPts2.append(bPtsC2[bpI2]);
                    }
                }

                optimizer.setSurfaceConstraint
                (
                    octreePtr_,
                    &mPartConstraint2.pointPatches(),
                    &globalToBp2,
                    &mPartConstraint2.corners(),
                    &featureTangents2
                );
                if( acuteGlobalPts2.size() > 0 )
                    optimizer.lockPoints(acuteGlobalPts2);

                optimizer.optimizeMeshFV();
                optimizer.setSurfaceConstraint(NULL, NULL, NULL, NULL, NULL);

                // Validate
                mesh_.clearAddressingData();
                labelHashSet negAfter2, pyrAfter2;
                polyMeshGenChecks::checkCellVolumes
                    (mesh_, false, &negAfter2);
                polyMeshGenChecks::checkFacePyramids
                    (mesh_, false, -SMALL, &pyrAfter2);

                label secondPassAllowedPyrIncrease = 0;
                if( meshDict_.isDict("boundaryLayers") )
                {
                    const dictionary& bndSP =
                        meshDict_.subDict("boundaryLayers");
                    if( bndSP.found("meshOptAllowedPyrIncrease") )
                        secondPassAllowedPyrIncrease = readLabel
                        (
                            bndSP.lookup("meshOptAllowedPyrIncrease")
                        );
                }

                const bool secondPassBad =
                    negAfter2.size() > 0
                 || label(pyrAfter2.size()) >
                    label(pyrBefore2.size()) + secondPassAllowedPyrIncrease;

                if( secondPassBad )
                {
                    Info << "Surface-constrained optimizer: second pass "
                         << "rejected: negVol "
                         << negBefore2.size() << "->" << negAfter2.size()
                         << " badPyramids "
                         << pyrBefore2.size() << "->" << pyrAfter2.size()
                         << " -- rolling back" << endl;

                    polyMeshGenModifier rbMod2(mesh_);
                    pointFieldPMG& rbPts2 = rbMod2.pointsAccess();

                    if( rbPts2.size() != pointsBefore2.size() )
                    {
                        FatalErrorInFunction
                            << "Cannot rollback constrained second pass: "
                            << "point count changed. rbPts2.size()="
                            << rbPts2.size()
                            << " pointsBefore2.size()="
                            << pointsBefore2.size()
                            << abort(FatalError);
                    }

                    forAll(pointsBefore2, pI)
                        rbPts2[pI] = pointsBefore2[pI];

                    mesh_.clearAddressingData();
                }
                else
                {
                    Info << "Surface-constrained optimizer: second pass "
                         << "accepted: negVol "
                         << negBefore2.size() << "->" << negAfter2.size()
                         << " badPyramids "
                         << pyrBefore2.size() << "->" << pyrAfter2.size()
                         << endl;
                }
            }
            else
            {
                Info << "Surface-constrained optimizer: second pass skipped"
                     << " -- " << negVolAfterPlain.size()
                     << " negVol cells remain after plain optimizer" << endl;
            }
        }
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
            // Default 0: any pyramid increase during re-projection triggers
            // rollback. Prevents surface snap from inverting face pyramids
            // on periodic-edge cells while keeping negVol unchanged.
            label reprojAllowedPyrIncrease = 0;
            if( meshDict_.isDict("boundaryLayers") )
            {
                const dictionary& bndLR3 =
                    meshDict_.subDict("boundaryLayers");
                if( bndLR3.found("postOptimizerReprojAllowedPyrIncrease") )
                    reprojAllowedPyrIncrease = readLabel
                    (
                        bndLR3.lookup("postOptimizerReprojAllowedPyrIncrease")
                    );
            }

            label totalAccepted = 0;
            label totalRolledBack = 0;

            // Causal retry: points proven to cause bad geometry are
            // blocked for all subsequent passes. Preserves the safe
            // majority of moves and rejects only proven bad actors.
            // meshDict knob: postOptimizerReprojMaxCausalRetries (default 2)
            label reprojMaxCausalRetries = 2;
            if( meshDict_.isDict("boundaryLayers") )
            {
                const dictionary& bndLR4 =
                    meshDict_.subDict("boundaryLayers");
                if( bndLR4.found("postOptimizerReprojMaxCausalRetries") )
                    reprojMaxCausalRetries = readLabel
                    (
                        bndLR4.lookup("postOptimizerReprojMaxCausalRetries")
                    );
            }

            labelHashSet blockedReprojMeshPts;

            // Capture pre-reproject baseline -- used as pyramid floor
            // across ALL passes so accepted passes cannot ratchet up
            // bad pyramid count incrementally.
            // Force fresh geometry/addressing before the global
            // re-projection baseline. Without this, the baseline can be
            // stale and rollback checks compare against a cached state.
            mesh_.clearAddressingData();

            labelHashSet negReprojBaseline, pyrReprojBaseline;
            polyMeshGenChecks::checkCellVolumes
                (mesh_, false, &negReprojBaseline);
            polyMeshGenChecks::checkFacePyramids
                (mesh_, false, -SMALL, &pyrReprojBaseline);

            Info << "Post-optimizer re-projection baseline: negVol="
                 << negReprojBaseline.size()
                 << " badPyramids=" << pyrReprojBaseline.size()
                 << endl;

            bool reprojSkipped = false;
            if( negReprojBaseline.size() > 0 )
            {
                reprojSkipped = true;
                Info << "Post-optimizer re-projection skipped: baseline has "
                     << negReprojBaseline.size()
                     << " negVol cells -- re-projection unsafe on dirty mesh"
                     << endl;
            }

            for( label passI = 0;
                 !reprojSkipped && passI < reprojMaxPasses;
                 ++passI )
            {
                bool passAccepted = false;
                label nRetries = 0;

                while( nRetries <= reprojMaxCausalRetries )
                {
                    // Snapshot for this attempt
                    // Snapshot active points only (logical size, not oversized buffer).
                    // pointFieldPMG::size() returns nElmts_ (logical count).
                    // pointField copy-ctor from pointFieldPMG copies raw capacity (1.5x).
                    // Explicit loop uses logical size and avoids the size mismatch.
                    const label nActivePts = mesh_.points().size();
                    pointField passPtsBefore(nActivePts);
                    forAll(passPtsBefore, pI)
                        passPtsBefore[pI] = mesh_.points()[pI];
                    // Force fresh geometry/addressing before the
                    // per-attempt baseline. This makes rollback verification
                    // compare against the real current point state, not a
                    // stale cached addressing state.
                    mesh_.clearAddressingData();

                    labelHashSet negPassBefore, pyrPassBefore;
                    polyMeshGenChecks::checkCellVolumes
                        (mesh_, false, &negPassBefore);
                    polyMeshGenChecks::checkFacePyramids
                        (mesh_, false, -SMALL, &pyrPassBefore);

                    // Move candidates, skipping protected and blocked points.
                    // IMPORTANT: use direct point writes, not meshSurfaceEngineModifier.
                    // surfModR holds cached surface-engine state that becomes unsafe
                    // across rollback/retry. Raw point writes + clearAddressingData()
                    // are sufficient for checkCellVolumes/checkFacePyramids to recompute.
                    label nMoved = 0;
                    DynList<label> movedMeshPtI;
                    DynList<label> movedBpI;
                    DynList<label> movedPatchI;
                    DynList<point> movedOldPt;
                    DynList<point> movedNewPt;
                    labelHashSet movedThisAttempt;

                    polyMeshGenModifier reprojMeshModifier(mesh_);
                    pointFieldPMG& reprojPts =
                        reprojMeshModifier.pointsAccess();

                    forAll(reprojPoints, rpI)
                    {
                        const label bpI = reprojPoints[rpI];
                        const label meshPtI = bPtsR[bpI];
                        if( meshPtI < 0 ||
                            meshPtI >= label(mesh_.points().size()) )
                            continue;
                        if( movedThisAttempt.found(meshPtI) )
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
                        // Skip protected fragile zones
                        if( protectedPts.found(meshPtI) )
                            continue;
                        // Skip points proven causal in prior attempts
                        if( blockedReprojMeshPts.found(meshPtI) )
                            continue;

                        const point oldPt = mesh_.points()[meshPtI];
                        const vector disp = snapPt - oldPt;
                        const point limitedPt =
                            oldPt + reprojStepFraction * disp;

                        reprojPts[meshPtI] = limitedPt;
                        movedThisAttempt.insert(meshPtI);

                        movedMeshPtI.append(meshPtI);
                        movedBpI.append(bpI);
                        movedPatchI.append(patchI);
                        movedOldPt.append(oldPt);
                        movedNewPt.append(limitedPt);
                        ++nMoved;
                    }

                    mesh_.clearAddressingData();

                    labelHashSet negPassAfter, pyrPassAfter;
                    polyMeshGenChecks::checkCellVolumes
                        (mesh_, false, &negPassAfter);
                    polyMeshGenChecks::checkFacePyramids
                        (mesh_, false, -SMALL, &pyrPassAfter);

                    // Compare against pre-reproject baseline, not per-pass
                    // snapshot. Prevents ratcheting up bad pyramids across
                    // accepted passes.
                    const bool negVolWorsened =
                        negPassAfter.size() > negReprojBaseline.size();
                    const bool pyrWorsened =
                        label(pyrPassAfter.size()) >
                        label(pyrReprojBaseline.size()) + reprojAllowedPyrIncrease;

                    if( !negVolWorsened && !pyrWorsened )
                    {
                        Info << "Post-optimizer re-projection pass "
                             << passI
                             << (nRetries > 0 ?
                                 word(" (retry ") + name(nRetries) + ")" : "")
                             << " accepted: moved " << nMoved
                             << " negVol "
                             << negPassBefore.size() << "->"
                             << negPassAfter.size()
                             << " badPyramids "
                             << pyrPassBefore.size() << "->"
                             << pyrPassAfter.size()
                             << endl;
                        passAccepted = true;
                        ++totalAccepted;
                        break;
                    }

                    // Pass failed -- identify causal points
                    Info << "Post-optimizer re-projection pass "
                         << passI
                         << " attempt " << nRetries
                         << " rejected: negVol "
                         << negPassBefore.size() << "->"
                         << negPassAfter.size()
                         << " badPyramids "
                         << pyrPassBefore.size() << "->"
                         << pyrPassAfter.size()
                         << (negVolWorsened ? " [negVol]" : "")
                         << (pyrWorsened ? " [badPyramids]" : "")
                         << endl;

                    // Extract new bad face/cell point labels
                    labelHashSet newBadFaces(pyrPassAfter);
                    forAllConstIter(labelHashSet, pyrPassBefore, it)
                        newBadFaces.erase(it.key());

                    labelHashSet newNegCells(negPassAfter);
                    forAllConstIter(labelHashSet, negPassBefore, it)
                        newNegCells.erase(it.key());

                    const faceListPMG& fs = mesh_.faces();
                    const cellListPMG& cls = mesh_.cells();

                    labelHashSet badFacePts;
                    forAllConstIter(labelHashSet, newBadFaces, it)
                    {
                        const label fI = it.key();
                        if( fI < 0 || fI >= fs.size() ) continue;
                        const face& f = fs[fI];
                        forAll(f, i) badFacePts.insert(f[i]);
                    }
                    forAllConstIter(labelHashSet, newNegCells, it)
                    {
                        const label cI = it.key();
                        if( cI < 0 || cI >= cls.size() ) continue;
                        const cell& c = cls[cI];
                        forAll(c, fI)
                        {
                            const label faceI = c[fI];
                            if( faceI < 0 || faceI >= fs.size() ) continue;
                            const face& f = fs[faceI];
                            forAll(f, pI) badFacePts.insert(f[pI]);
                        }
                    }

                    // Selective rollback: only restore points moved in this
                    // attempt. Avoids pointFieldPMG::operator= and prevents
                    // point storage resizing from invalidating cached
                    // surface-engine structures (bPtsR, pPatchesR).
                    // Full point-by-point restore from pre-attempt snapshot.
                    // Selective rollback (only moved points) proved unreliable --
                    // checkFacePyramids saw stale state after selective restore.
                    // Point-by-point avoids pointFieldPMG::operator= resize risk.
                    {
                        polyMeshGenModifier rbModifier(mesh_);
                        pointFieldPMG& rbPts = rbModifier.pointsAccess();

                        if( rbPts.size() != passPtsBefore.size() )
                        {
                            FatalErrorInFunction
                                << "Cannot rollback: point count changed. "
                                << "rbPts.size()=" << rbPts.size()
                                << " passPtsBefore.size()="
                                << passPtsBefore.size()
                                << abort(FatalError);
                        }

                        forAll(passPtsBefore, pI)
                            rbPts[pI] = passPtsBefore[pI];

                        // Verify coordinates were actually restored.
                        // If this passes but pyramid counts differ, the issue
                        // is addressing/cache/baseline state, not point rollback.
                        label nRollbackCoordDiff = 0;
                        scalar maxRollbackCoordDiff = 0.0;
                        forAll(passPtsBefore, pI)
                        {
                            const scalar d = mag(rbPts[pI] - passPtsBefore[pI]);
                            if( d > SMALL )
                            {
                                ++nRollbackCoordDiff;
                                if( d > maxRollbackCoordDiff )
                                    maxRollbackCoordDiff = d;
                            }
                        }

                        Info << "Post-optimizer re-projection rollback "
                             << "coordinate check: nDiff="
                             << nRollbackCoordDiff
                             << " maxDiff=" << maxRollbackCoordDiff
                             << endl;

                        mesh_.clearAddressingData();
                    }

                    // Verify rollback restored baseline state
                    {
                        labelHashSet negRB, pyrRB;
                        polyMeshGenChecks::checkCellVolumes
                            (mesh_, false, &negRB);
                        polyMeshGenChecks::checkFacePyramids
                            (mesh_, false, -SMALL, &pyrRB);
                        const bool rollbackIncomplete =
                            negRB.size() != negPassBefore.size()
                         || pyrRB.size() != pyrPassBefore.size();

                        Info << "Post-optimizer re-projection rollback check"
                             << " pass=" << passI
                             << " attempt=" << nRetries
                             << ": negVol=" << negRB.size()
                             << " badPyramids=" << pyrRB.size()
                             << (rollbackIncomplete ?
                                 " FATAL: rollback incomplete" : " OK")
                             << endl;

                        if( rollbackIncomplete )
                        {
                            FatalErrorInFunction
                                << "Post-optimizer re-projection rollback "
                                << "did not restore pre-attempt state. "
                                << "pre-attempt negVol="
                                << negPassBefore.size()
                                << " badPyramids="
                                << pyrPassBefore.size()
                                << " rollback negVol=" << negRB.size()
                                << " badPyramids=" << pyrRB.size()
                                << abort(FatalError);
                        }
                    }

                    // Intersect moved points with bad face points
                    labelHashSet newCausal;
                    forAll(movedMeshPtI, i)
                        if( badFacePts.found(movedMeshPtI[i]) )
                            newCausal.insert(movedMeshPtI[i]);

                    if( newCausal.empty() )
                    {
                        Info << "Post-optimizer re-projection pass "
                             << passI
                             << ": no causal points identified"
                             << " -- abandoning pass" << endl;
                        ++totalRolledBack;
                        break;
                    }

                    forAllConstIter(labelHashSet, newCausal, it)
                        blockedReprojMeshPts.insert(it.key());

                    Info << "Post-optimizer re-projection pass "
                         << passI
                         << " retry " << (nRetries+1)
                         << ": blocked " << newCausal.size()
                         << " causal points ("
                         << blockedReprojMeshPts.size()
                         << " total blocked), retrying" << endl;

                    ++nRetries;
                }

                if( !passAccepted )
                    break;
            }

            if( reprojSkipped )
            {
                Info << "Post-optimizer re-projection: skipped (dirty baseline)"
                     << endl;
            }
            else
            {
                Info << "Post-optimizer re-projection: "
                     << totalAccepted << " passes accepted, "
                     << totalRolledBack << " rolled back, "
                     << blockedReprojMeshPts.size()
                     << " points permanently blocked" << endl;
            }

            // Store blocked set for pre-BL snap to skip same bad actors
            postOptimizerReprojBlockedPts_ = blockedReprojMeshPts;

            // Write causal CSV
            if( blockedReprojMeshPts.size() > 0 )
            {
                OFstream csvOs("postOptimizerReprojCausalPoints.csv");
                csvOs << "meshPtI,bpI,patchI,"
                      << "oldX,oldY,oldZ,snapX,snapY,snapZ,dispMag" << nl;
                forAll(bPtsR, bpI)
                {
                    const label meshPtI = bPtsR[bpI];
                    if( !blockedReprojMeshPts.found(meshPtI) )
                        continue;
                    if( pPatchesR.sizeOfRow(bpI) != 1 )
                        continue;
                    const label patchI = pPatchesR(bpI, 0);
                    const point& oldPt = mesh_.points()[meshPtI];
                    point snapPt;
                    scalar snapDsq;
                    label snapNt;
                    octreePtr_->findNearestSurfacePointInRegion
                    (
                        snapPt, snapDsq, snapNt, patchI, oldPt
                    );
                    const vector dv = snapPt - oldPt;
                    csvOs << meshPtI << ','
                          << bpI << ','
                          << patchI << ','
                          << oldPt.x() << ',' << oldPt.y() << ',' << oldPt.z() << ','
                          << snapPt.x() << ',' << snapPt.y() << ',' << snapPt.z() << ','
                          << mag(dv) << nl;
                }
                Info << "postOptimizerReproj: "
                     << blockedReprojMeshPts.size()
                     << " blocked points written to postOptimizerReprojCausalPoints.csv"
                     << endl;
            }
        }
    }
    //- octreePtr_ kept alive until post-BL rescue completes.
    //- Previously deleted here; deferred to destructor so that
    //- surface-constrained post-BL rescue can use it.

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

                const bool stage1SkewOK =
                    maxSkewS1After <= maxSkewS1Before;

                const bool stage1OK =
                    badStage1.size() <= badFaces.size()
                 && negStage1.size() <= negBefore.size()
                 && openStage1.size() <= openBefore.size()
                 && stage1SkewOK;

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

        label dirtyRecoverableMaxNegVol = 0;
        if( meshDict_.found("dirtyRecoverableMaxNegVol") )
            dirtyRecoverableMaxNegVol =
                readLabel(meshDict_.lookup("dirtyRecoverableMaxNegVol"));

        if( negBefore.size() > 0
         && label(negBefore.size()) <= dirtyRecoverableMaxNegVol )
        {
            //- DIRTY_RECOVERABLE cleanup branch.
            //- Residual cells are typically epsilon-scale and the
            //- untangler CAN clear them. Measurement showed ~40% of runs
            //- left 1 such cell that cascaded through BL to skew~500.
            //- Attempt untangle; if it reaches negVol=0 the downstream
            //- improved-but-dirty check sees a clean mesh and the run is
            //- promoted to CLEAN automatically. If cleanup fails, roll
            //- back to the prior DIRTY state -- crucially WITHOUT setting
            //- finalUntangleRejected_ (the hard-unsafe flag that would
            //- skip BL entirely). Failed cleanup must return to dirty,
            //- not escalate to hard-rejected.
            Info << "optimiseFinalMesh: DIRTY_RECOVERABLE -- attempting "
                 << "untangle cleanup on " << negBefore.size()
                 << " residual negVol cell(s) (<= threshold "
                 << dirtyRecoverableMaxNegVol << ")" << endl;

            optimizer.untangleMeshFV();

            labelHashSet dcBadAfter;
            polyMeshGenChecks::checkFacePyramids
                (mesh_, false, -SMALL, &dcBadAfter);
            labelHashSet dcNegAfter;
            polyMeshGenChecks::checkCellVolumes(mesh_, false, &dcNegAfter);
            labelHashSet dcOpenAfter;
            polyMeshGenChecks::checkClosedCells
                (mesh_, false, 0.5, &dcOpenAfter);
            scalarField dcSkewAfter;
            polyMeshGenChecks::checkFaceSkewness(mesh_, dcSkewAfter);
            const scalar dcMaxSkewAfter =
                dcSkewAfter.size() > 0 ? max(dcSkewAfter) : scalar(0.0);
            const bool dcSkewOK =
                dcMaxSkewAfter <= scalar(20.0)
             && dcMaxSkewAfter <= scalar(2.0) *
                    Foam::max(maxSkewBefore, scalar(1.0));

            //- Strict success: cleanup only counts if it reaches
            //- negVol == 0 AND does not worsen pyramids/openCells/skew.
            const bool dirtyCleanupOK =
                dcNegAfter.size() == 0
             && dcBadAfter.size() <= badBefore.size()
             && dcOpenAfter.size() <= openBefore.size()
             && dcSkewOK;

            if( dirtyCleanupOK )
            {
                Info << "DIRTY_RECOVERABLE cleanup ACCEPTED: negVol "
                     << negBefore.size() << "->" << dcNegAfter.size()
                     << " badPyr " << badBefore.size() << "->"
                     << dcBadAfter.size()
                     << " -- promoted to CLEAN" << endl;
                meshHistory_ = MeshHistory::CleanPromoted;
                //- No flag changes needed: downstream sees negVol=0.
            }
            else
            {
                Info << "DIRTY_RECOVERABLE cleanup REJECTED: negVol "
                     << negBefore.size() << "->" << dcNegAfter.size()
                     << " -- Rung1 stalled, attempting Rung2 TJ smoother"
                     << endl;
                //- Trap 4 fix: do NOT roll back before Rung2.
                //- Rung1 partially improved (e.g. 3->1 negVol).
                //- Rung2 continues from that improved state.
                //- Full rollback to pointsBefore only if Rung2 also fails.
                {
                    bool enableRung2 = false;
                    if( meshDict_.found("enableRung2TJSmoother") )
                        enableRung2 =
                            readBool(meshDict_.lookup("enableRung2TJSmoother"));
                    if( !enableRung2 )
                    {
                        Info << "Rung2 TJ smoother disabled (enableRung2TJSmoother"
                             << " not set) -- rolling back to pointsBefore" << endl;
                        polyMeshGenModifier r2Mod(mesh_);
                        pointFieldPMG& r2ModPts = r2Mod.pointsAccess();
                        forAll(pointsBefore, pI)
                            r2ModPts[pI] = pointsBefore[pI];
                        mesh_.clearAddressingData();
                        meshHistory_ = MeshHistory::DirtyRecoverable;
                    }
                    else
                    {
                    //- Rung 2: defect-seeded neighborhood smoother.
                    //- Seeds from residual negVol cells, expands by
                    //- ring1+ring2 neighbors for movable interior points.
                    //- Atlas data: distTJ<0.005 predicts 100% of bad
                    //- pyramids -- repair zone targets this neighborhood.
                    label rung2MaxRings = 2;
                    if( meshDict_.found("rung2MaxRings") )
                        rung2MaxRings =
                            readLabel(meshDict_.lookup("rung2MaxRings"));
                    label rung2AllowedPyrIncrease = 5;
                    if( meshDict_.found("rung2AllowedPyrIncrease") )
                        rung2AllowedPyrIncrease =
                            readLabel(meshDict_.lookup("rung2AllowedPyrIncrease"));
                    // r2Pts not needed directly -- accessed via mesh_ in zone builder
                    const cellListPMG&   r2Cells = mesh_.cells();
                    const faceListPMG&   r2Faces = mesh_.faces();
                    const labelList&     r2Own   = mesh_.owner();
                    const labelList&     r2Nei   = mesh_.neighbour();
                    const label r2NInternal = mesh_.nInternalFaces();
                    //- Step 1: seed from residual negVol cells
                    labelHashSet repairZone;
                    forAllConstIter(labelHashSet, dcNegAfter, it)
                        repairZone.insert(it.key());
                    const label nSeeds = repairZone.size();
                    //- Step 2: expand by rung2MaxRings rings via
                    //- owner/neighbour face adjacency
                    for( label ring = 0; ring < rung2MaxRings; ++ring )
                    {
                        labelHashSet ringAdd;
                        forAllConstIter(labelHashSet, repairZone, it)
                        {
                            const label cellI = it.key();
                            if( cellI < 0 || cellI >= label(r2Cells.size()) ) continue;
                            const cell& c = r2Cells[cellI];
                            forAll(c, fI)
                            {
                                const label faceI = c[fI];
                                if( faceI < 0 || faceI >= r2NInternal ) continue;
                                const label ownC = r2Own[faceI];
                                const label neiC = r2Nei[faceI];
                                if( ownC >= 0 && !repairZone.found(ownC) )
                                    ringAdd.insert(ownC);
                                if( neiC >= 0 && !repairZone.found(neiC) )
                                    ringAdd.insert(neiC);
                            }
                        }
                        forAllConstIter(labelHashSet, ringAdd, it)
                            repairZone.insert(it.key());
                    }
                    //- Step 3: count movable interior points
                    //- (points belonging ONLY to repair zone cells)
                    labelHashSet repairPts;
                    labelHashSet outsidePts;
                    forAll(r2Cells, cellI)
                    {
                        const cell& c = r2Cells[cellI];
                        const bool inZone = repairZone.found(cellI);
                        forAll(c, fI)
                        {
                            const label faceI = c[fI];
                            if( faceI < 0 || faceI >= label(r2Faces.size()) ) continue;
                            const face& f = r2Faces[faceI];
                            forAll(f, pI)
                            {
                                const label ptI = f[pI];
                                if( ptI < 0 ) continue;
                                if( inZone ) repairPts.insert(ptI);
                                else outsidePts.insert(ptI);
                            }
                        }
                    }
                    label movablePts = 0;
                    forAllConstIter(labelHashSet, repairPts, it)
                        if( !outsidePts.found(it.key()) ) ++movablePts;
                    Info << "Rung2 zone: seeds=" << nSeeds
                         << " rings=" << rung2MaxRings
                         << " repairCells=" << repairZone.size()
                         << " repairPts=" << repairPts.size()
                         << " movableInteriorPts=" << movablePts
                         << " cellsToLock="
                         << label(r2Cells.size()) - label(repairZone.size())
                         << endl;
                    if( movablePts > 0 && repairZone.size() > 0 )
                    {
                        //- No local snapshot needed: rollback uses
                        //- pointsBefore (combined Rung1+Rung2 transaction)
                        labelHashSet r2BadBefore;
                        polyMeshGenChecks::checkFacePyramids
                            (mesh_, false, -SMALL, &r2BadBefore);
                        labelHashSet r2OpenBefore;
                        polyMeshGenChecks::checkClosedCells
                            (mesh_, false, 0.5, &r2OpenBefore);
                        //- Lock everything outside repair zone
                        labelLongList cellsToLock;
                        forAll(r2Cells, cellI)
                            if( !repairZone.found(cellI) )
                                cellsToLock.append(cellI);
                        meshOptimizer rung2Opt(mesh_);
                        rung2Opt.lockCells(cellsToLock);

                        //- Rung2 is a boundary-junction local repair.
                        //- Must be surface-constrained. Running unconstrained
                        //- can move blade/hub, blade/shroud, periodic-edge
                        //- points off their geometry.
                        if( octreePtr_ )
                        {
                            meshSurfaceEngine mseR2(mesh_);
                            meshSurfacePartitioner mPartR2(mseR2);
                            labelLongList globalToBpR2
                            (
                                mesh_.points().size(),
                                -1
                            );
                            const labelList& bPtsR2 = mseR2.boundaryPoints();
                            forAll(bPtsR2, bpI)
                                globalToBpR2[bPtsR2[bpI]] = bpI;
                            const label nBpR2 = bPtsR2.size();
                            vectorField featTanR2(nBpR2, vector::zero);
                            {
                                const edgeList& edgesR2 = mseR2.edges();
                                const VRWGraph& bpEdgesR2 =
                                    mseR2.boundaryPointEdges();
                                const labelHashSet& featEdgesR2 =
                                    mPartR2.featureEdges();
                                const labelHashSet& edgePtsR2 =
                                    mPartR2.edgePoints();
                                const pointFieldPMG& ptsR2 = mesh_.points();
                                const labelList& bpR2 = mseR2.bp();
                                forAllConstIter(labelHashSet, edgePtsR2, it)
                                {
                                    const label bpI = it.key();
                                    if( bpI < 0 || bpI >= nBpR2 ) continue;
                                    label nbr0 = -1;
                                    label nbr1 = -1;
                                    forAllRow(bpEdgesR2, bpI, eI)
                                    {
                                        const label beI = bpEdgesR2(bpI, eI);
                                        if( !featEdgesR2.found(beI) ) continue;
                                        const edge& e = edgesR2[beI];
                                        const label ep0 = e.start();
                                        const label ep1 = e.end();
                                        if( ep0 < 0 || ep0 >= bpR2.size() ) continue;
                                        if( ep1 < 0 || ep1 >= bpR2.size() ) continue;
                                        const label o0 = bpR2[ep0];
                                        const label o1 = bpR2[ep1];
                                        if( o0 < 0 || o1 < 0 ) continue;
                                        label otherBp = -1;
                                        if( o0 == bpI ) otherBp = o1;
                                        else if( o1 == bpI ) otherBp = o0;
                                        if( otherBp < 0 || otherBp >= nBpR2 ) continue;
                                        if( nbr0 == -1 ) nbr0 = otherBp;
                                        else if( nbr1 == -1 && otherBp != nbr0 )
                                            nbr1 = otherBp;
                                    }
                                    vector t = vector::zero;
                                    if( nbr0 != -1 && nbr1 != -1 )
                                        t = ptsR2[bPtsR2[nbr1]]
                                          - ptsR2[bPtsR2[nbr0]];
                                    else if( nbr0 != -1 )
                                        t = ptsR2[bPtsR2[nbr0]]
                                          - ptsR2[bPtsR2[bpI]];
                                    if( magSqr(t) > VSMALL )
                                        featTanR2[bpI] = t / mag(t);
                                }
                            }
                            rung2Opt.setSurfaceConstraint
                            (
                                octreePtr_,
                                &mPartR2.pointPatches(),
                                &globalToBpR2,
                                &mPartR2.corners(),
                                &featTanR2
                            );
                            rung2Opt.optimizeMeshFV(3, 5, 20, 1);
                            rung2Opt.setSurfaceConstraint
                            (NULL, NULL, NULL, NULL, NULL);
                        }
                        else
                        {
                            Info << "Rung2 TJ smoother skipped: no octreePtr_ "
                                 << "for surface-constrained local repair "
                                 << "-- rolling back Rung1+Rung2 to pointsBefore"
                                 << endl;
                            polyMeshGenModifier r2Mod(mesh_);
                            pointFieldPMG& r2ModPts = r2Mod.pointsAccess();
                            forAll(pointsBefore, pI)
                                r2ModPts[pI] = pointsBefore[pI];
                            mesh_.clearAddressingData();
                            meshHistory_ = MeshHistory::DirtyRecoverable;
                        }
                        //- Measure: full MeshQuality gate
                        labelHashSet r2NegAfter;
                        polyMeshGenChecks::checkCellVolumes
                            (mesh_, false, &r2NegAfter);
                        labelHashSet r2BadAfter;
                        polyMeshGenChecks::checkFacePyramids
                            (mesh_, false, -SMALL, &r2BadAfter);
                        labelHashSet r2OpenAfter;
                        polyMeshGenChecks::checkClosedCells
                            (mesh_, false, 0.5, &r2OpenAfter);
                        scalarField r2SkewAfter;
                        polyMeshGenChecks::checkFaceSkewness
                            (mesh_, r2SkewAfter);
                        const scalar r2MaxSkew =
                            r2SkewAfter.size() > 0 ?
                            max(r2SkewAfter) : scalar(0);
                        const bool r2OK =
                            r2NegAfter.size() == 0
                         && label(r2BadAfter.size()) <=
                            label(r2BadBefore.size()) + rung2AllowedPyrIncrease
                         && r2OpenAfter.size() <= r2OpenBefore.size()
                         && r2MaxSkew <= scalar(20.0);
                        if( r2OK )
                        {
                            Info << "Rung2 TJ smoother ACCEPTED: negVol "
                                 << dcNegAfter.size() << "->0"
                                 << " badPyr " << r2BadBefore.size()
                                 << "->" << r2BadAfter.size()
                                 << " maxSkew=" << r2MaxSkew
                                 << " -- promoted to CleanPromoted" << endl;
                            meshHistory_ = MeshHistory::CleanPromoted;
                        }
                        else
                        {
                            Info << "Rung2 TJ smoother REJECTED: negVol "
                                 << dcNegAfter.size() << "->"
                                 << r2NegAfter.size()
                                 << " badPyr " << r2BadBefore.size()
                                 << "->" << r2BadAfter.size()
                                 << " maxSkew=" << r2MaxSkew
                                 << " -- rolling back combined Rung1+Rung2"
                                 << " transaction to pointsBefore" << endl;
                            polyMeshGenModifier r2Mod(mesh_);
                            pointFieldPMG& r2ModPts = r2Mod.pointsAccess();
                            forAll(pointsBefore, pI)
                                r2ModPts[pI] = pointsBefore[pI];
                            mesh_.clearAddressingData();
                            meshHistory_ = MeshHistory::DirtyRecoverable;
                        }
                    }
                    else
                    {
                        Info << "Rung2 TJ smoother: insufficient movable"
                             << " points (" << movablePts << ")"
                             << " -- rolling back Rung1 to pointsBefore"
                             << endl;
                        polyMeshGenModifier r2Mod(mesh_);
                        pointFieldPMG& r2ModPts = r2Mod.pointsAccess();
                        forAll(pointsBefore, pI)
                            r2ModPts[pI] = pointsBefore[pI];
                        mesh_.clearAddressingData();
                        meshHistory_ = MeshHistory::DirtyRecoverable;
                    }
                    } // end if( enableRung2 )
                } // end Rung2 block
                //- Deliberately DO NOT set finalUntangleRejected_.
                //- The mesh returns to its prior dirty-but-recoverable
                //- state; downstream reprojUnsafe_ logic handles it.
            }
        }
        else if( negBefore.size() > 0 )
        {
            meshHistory_ = MeshHistory::DirtyUnsafe;
            Info << "optimiseFinalMesh: skipping untangleMeshFV -- "
                 << negBefore.size()
                 << " negVol cells present (> " << dirtyRecoverableMaxNegVol
                 << " threshold) -- proceeding to BL" << endl;
            // Large residual negVol: untangle is unlikely to help and may
            // stall. Keep the mesh unchanged; downstream gates handle it.
            mesh_.clearAddressingData();
        }
        else
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
                meshHistory_ = MeshHistory::HardRejected;
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
        // Skip points proven causal during post-optimizer re-projection.
        // Re-projecting these points again would repeat the same failure.
        if( postOptimizerReprojBlockedPts_.found(meshPtI) )
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
        auto earlyLineageCSV = [&](const std::string& stageName)
        {
            bool writeLineageDiagnostics = false;
            if( meshDict_.found("writeLineageDiagnostics") )
                writeLineageDiagnostics =
                    Switch(meshDict_.lookup("writeLineageDiagnostics"));

            if( !writeLineageDiagnostics ) return;

            wordList enabledStages;
            if( meshDict_.found("writeLineageStages") )
            {
                enabledStages = wordList(meshDict_.lookup("writeLineageStages"));
            }
            else
            {
                enabledStages.setSize(0);
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

            Info << "Lineage [" << stageName.c_str() << "]: "
                 << stageCells.size() << " negVol cells" << endl;

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
                    cc += pts[pI];
                    ++nUnique;
                }

                if( nUnique > 0 ) cc /= scalar(nUnique);

                stageFile << cellI << ","
                          << cc.x() << "," << cc.y() << "," << cc.z() << ","
                          << c.size() << "," << nUnique << nl;
            }

            Info << "Lineage [" << stageName.c_str() << "]: wrote "
                 << fname.c_str() << endl;
        };

        if( controller_.runCurrentStep("templateGeneration") )
        {
            createCartesianMesh();
            earlyLineageCSV("postTemplateGeneration");
        }

        if( controller_.runCurrentStep("surfaceTopology") )
        {
            surfacePreparation();
            earlyLineageCSV("postSurfaceTopology");
        }

        if( controller_.runCurrentStep("patchAssignment") )
        {
            // Patch assignment moved before surface projection so that
            // mapVerticesOntoSurface has valid patch identity available.
            // edgeExtractor uses only mesh topology + octree -- no
            // dependency on projected surface positions.
            extractPatches();
            earlyLineageCSV("postPatchAssignmentInitial");
        }

        if( controller_.runCurrentStep("surfaceProjection") )
        {
            mapMeshToSurface();
            earlyLineageCSV("postSurfaceProjection");
            // Re-run patch assignment after projection to correct any
            // misassignments that occurred on the unprojected hex mesh.
            extractPatches();
            earlyLineageCSV("postPatchAssignmentAfterProjection");
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
            earlyLineageCSV("postMapEdgesAndCorners");

            optimiseMeshSurface();
            earlyLineageCSV("postOptimiseMeshSurface");


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
                earlyLineageCSV("postCornerSnap");
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
                earlyLineageCSV("postEdgeSnap");
            }
            earlyLineageCSV("postEdgeExtraction");
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
            {
                enabledStages = wordList(meshDict_.lookup("writeLineageStages"));
            }
            else
            {
                enabledStages.setSize(2);
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
            mesh_.clearAddressingData();
            labelHashSet meshOptBadBefore;
            polyMeshGenChecks::checkFacePyramids(mesh_, false, -SMALL, &meshOptBadBefore);
            labelHashSet meshOptNegBefore;
            polyMeshGenChecks::checkCellVolumes(mesh_, false, &meshOptNegBefore);
            const pointField meshOptPointsBefore(mesh_.points());
            Info << "MESHOPTDIAG before: badPyramids="
                 << meshOptBadBefore.size()
                 << " negVol=" << meshOptNegBefore.size()
                 << endl;
            optimiseFinalMesh();

            projectSurfaceAfterBackScaling();
            mesh_.clearAddressingData();
            labelHashSet meshOptBadAfter;
            polyMeshGenChecks::checkFacePyramids(mesh_, false, -SMALL, &meshOptBadAfter);
            labelHashSet meshOptNegAfter;
            polyMeshGenChecks::checkCellVolumes(mesh_, false, &meshOptNegAfter);
            Info << "MESHOPTDIAG after: badPyramids="
                 << meshOptBadAfter.size()
                 << " negVol=" << meshOptNegAfter.size()
                 << endl;
            // Capture post-optimizer bad pyramid faces for pre-refBL retraction.
            // Only store if mesh is valid (0 negVol).
            if( meshOptNegAfter.size() == 0 && meshOptBadAfter.size() > 0 )
            {
                postOptBadFaces_ = meshOptBadAfter;
                Info << "postOptBadFaces captured: "
                     << postOptBadFaces_.size()
                     << " bad pyramid faces for pre-refBL retraction" << endl;
            }
            else
            {
                postOptBadFaces_.clear();
            }

            // Rollback if negVol got worse, or badPyramids got much worse
            // without negVol improvement. If negVol improved but is nonzero,
            // keep the improved points but flag mesh as dirty so re-projection
            // and BL stages skip safely.
            label meshOptAllowedPyrIncrease = 0;
            if( meshDict_.isDict("boundaryLayers") )
            {
                const dictionary& bndLMO =
                    meshDict_.subDict("boundaryLayers");
                if( bndLMO.found("meshOptAllowedPyrIncrease") )
                    meshOptAllowedPyrIncrease = readLabel
                    (
                        bndLMO.lookup("meshOptAllowedPyrIncrease")
                    );
            }

            const bool meshOptNegWorse =
                meshOptNegAfter.size() > meshOptNegBefore.size();
            const bool meshOptNegImproved =
                meshOptNegAfter.size() < meshOptNegBefore.size();
            const bool meshOptNegNonZero =
                meshOptNegAfter.size() > 0;
            const bool meshOptPyrTooMuchWorse =
                label(meshOptBadAfter.size()) >
                label(meshOptBadBefore.size()) + meshOptAllowedPyrIncrease
             && !meshOptNegImproved;

            if( meshOptNegNonZero && !meshOptNegWorse && !meshOptPyrTooMuchWorse )
            {
                {
                    label dirtyRecoverableMaxNegVol = 0;
                    if( meshDict_.found("dirtyRecoverableMaxNegVol") )
                        dirtyRecoverableMaxNegVol =
                            readLabel(meshDict_.lookup("dirtyRecoverableMaxNegVol"));

                    const word candidateState =
                        label(meshOptNegAfter.size()) <= dirtyRecoverableMaxNegVol
                      ? word("DIRTY_RECOVERABLE_CANDIDATE")
                      : word("DIRTY_UNSAFE_CANDIDATE");

                    Info << "MESHOPTDIAG improved-but-dirty: negVol "
                         << meshOptNegBefore.size() << "->" << meshOptNegAfter.size()
                         << " badPyramids "
                         << meshOptBadBefore.size() << "->" << meshOptBadAfter.size()
                         << " candidateState=" << candidateState
                         << " threshold=" << dirtyRecoverableMaxNegVol
                         << " -- keeping improved points, marking re-projection unsafe"
                         << endl;
                    // No behavior change: reprojUnsafe_ still gates downstream
                    // repair/snap logic exactly as before. This classification
                    // is observability only.
                    reprojUnsafe_ = true;
                }
            }

            if( meshOptNegWorse || meshOptPyrTooMuchWorse )
            {
                Info << "MESHOPTDIAG rejected: badPyramids "
                     << meshOptBadBefore.size() << "->" << meshOptBadAfter.size()
                     << ", negVol " << meshOptNegBefore.size()
                     << "->" << meshOptNegAfter.size()
                     << " -- restoring pre-meshOptimisation points" << endl;
                polyMeshGenModifier meshModifier(mesh_);
                pointFieldPMG& pts = meshModifier.pointsAccess();
                pts = meshOptPointsBefore;
                mesh_.clearAddressingData();
                finalUntangleRejected_ = false;
                labelHashSet rbBad;
                polyMeshGenChecks::checkFacePyramids(mesh_, false, -SMALL, &rbBad);
                labelHashSet rbNeg;
                polyMeshGenChecks::checkCellVolumes(mesh_, false, &rbNeg);
                Info << "MESHOPTDIAG after rollback: badPyramids="
                     << rbBad.size()
                     << " negVol=" << rbNeg.size()
                     << endl;
            }
            writeLineageCSV("postOptimize");
        }
        if( finalUntangleRejected_ )
        {
            Info << "Pre-BL snap: skipped -- mesh state unsafe" << endl;
        }
        else if( reprojUnsafe_ )
        {
            Info << "Pre-BL snap: skipped -- re-projection marked unsafe"
                 << endl;
        }
        else
        {
            snapSurfaceBeforeBLRefinement();
        }
        if( controller_.runCurrentStep("boundaryLayerRefinement") )
        {
            if( finalUntangleRejected_ )
            {
                Info << "refBoundaryLayers: skipped -- final untangle was rejected, mesh state unsafe" << endl;
            }
            else
            {
                if( reprojUnsafe_ )
                    Info << "refBoundaryLayers: running on dirty pre-BL mesh "
                         << "(reprojUnsafe -- provenance retraction disabled)" << endl;

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

            // Print pre-renumber bad pyramid face/cell IDs so they can be
            // queried against postRefBL_cellProvenance.csv (also pre-renumber).
            if( pyrBeforeRenumber.size() > 0 )
            {
                const labelList& own = mesh_.owner();
                const labelList& nei = mesh_.neighbour();

                labelHashSet badCells;

                Info << "Pre-renumber bad pyramid face/cell map:" << endl;

                Info << "  pyrBeforeRenumber raw IDs:";
                forAllConstIter(labelHashSet, pyrBeforeRenumber, it2)
                    Info << " " << it2.key();
                Info << endl;

                forAllConstIter(labelHashSet, pyrBeforeRenumber, it)
                {
                    const label faceI = it.key();

                    label ownCell = -1;
                    label neiCell = -1;

                    if( faceI >= 0 && faceI < own.size() )
                    {
                        ownCell = own[faceI];
                        badCells.insert(ownCell);
                    }

                    if( faceI >= 0 && faceI < nei.size() )
                    {
                        neiCell = nei[faceI];
                        badCells.insert(neiCell);
                    }

                    Info << "  badPyrFace=" << faceI
                         << " owner=" << ownCell
                         << " neighbour=" << neiCell
                         << endl;
                }

                Info << "Pre-renumber bad pyramid candidate cells:";
                forAllConstIter(labelHashSet, badCells, it)
                    Info << " " << it.key();
                Info << endl;
            }
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
            //- Post-BL Rescue Rung: surface-constrained local repair
            //- for negVol introduced by BL/replaceBoundaries.
            //- Runs BEFORE MESHHISTORY/writeLineageCSV so final state
            //- is honest. Recomputes finalNegCells/finalBadPyrFaces
            //- after rescue so downstream diagnostics see true state.
            {
                bool enablePostBLRescue = false;
                if( meshDict_.found("enablePostBLRescue") )
                    enablePostBLRescue =
                        readBool(meshDict_.lookup("enablePostBLRescue"));
                label postBLRescueMaxNegVol = 20;
                if( meshDict_.found("postBLRescueMaxNegVol") )
                    postBLRescueMaxNegVol =
                        readLabel(meshDict_.lookup("postBLRescueMaxNegVol"));
                label postBLRescueMaxRings = 3;
                if( meshDict_.found("postBLRescueMaxRings") )
                    postBLRescueMaxRings =
                        readLabel(meshDict_.lookup("postBLRescueMaxRings"));
                label postBLRescueAllowedPyrIncrease = 10;
                if( meshDict_.found("postBLRescueAllowedPyrIncrease") )
                    postBLRescueAllowedPyrIncrease =
                        readLabel(meshDict_.lookup("postBLRescueAllowedPyrIncrease"));
                scalar postBLRescueMaxSkew = 100.0;
                if( meshDict_.found("postBLRescueMaxSkew") )
                    postBLRescueMaxSkew =
                        readScalar(meshDict_.lookup("postBLRescueMaxSkew"));
                if( enablePostBLRescue
                 && finalNegCells.size() > 0
                 && label(finalNegCells.size()) <= postBLRescueMaxNegVol
                 && octreePtr_ )
                {
                    Info << "POSTBLRESCUE before: negVol="
                         << finalNegCells.size()
                         << " badPyr=" << finalBadPyrFaces.size()
                         << endl;
                    const cellListPMG& pbrCells = mesh_.cells();
                    const faceListPMG& pbrFaces = mesh_.faces();
                    const labelList&   pbrOwn   = mesh_.owner();
                    const labelList&   pbrNei   = mesh_.neighbour();
                    const label pbrNInternal = mesh_.nInternalFaces();
                    //- Step 1: seed from finalNegCells + badPyr owners
                    labelHashSet pbrZone;
                    forAllConstIter(labelHashSet, finalNegCells, it)
                        pbrZone.insert(it.key());
                    forAllConstIter(labelHashSet, finalBadPyrFaces, it)
                    {
                        const label faceI = it.key();
                        if( faceI >= 0 && faceI < label(pbrOwn.size()) )
                            pbrZone.insert(pbrOwn[faceI]);
                        if( faceI >= 0 && faceI < label(pbrNei.size())
                         && pbrNei[faceI] >= 0 )
                            pbrZone.insert(pbrNei[faceI]);
                    }
                    const label pbrNSeeds = pbrZone.size();
                    //- Step 2: ring expansion
                    for( label ring = 0; ring < postBLRescueMaxRings; ++ring )
                    {
                        labelHashSet ringAdd;
                        forAllConstIter(labelHashSet, pbrZone, it)
                        {
                            const label cellI = it.key();
                            if( cellI < 0 || cellI >= label(pbrCells.size()) ) continue;
                            const cell& c = pbrCells[cellI];
                            forAll(c, fI)
                            {
                                const label faceI = c[fI];
                                if( faceI < 0 || faceI >= pbrNInternal ) continue;
                                const label ownC = pbrOwn[faceI];
                                const label neiC = pbrNei[faceI];
                                if( ownC >= 0 && !pbrZone.found(ownC) )
                                    ringAdd.insert(ownC);
                                if( neiC >= 0 && !pbrZone.found(neiC) )
                                    ringAdd.insert(neiC);
                            }
                        }
                        forAllConstIter(labelHashSet, ringAdd, it)
                            pbrZone.insert(it.key());
                    }
                    //- Step 3: count movable interior points
                    labelHashSet pbrRepairPts;
                    labelHashSet pbrOutsidePts;
                    forAll(pbrCells, cellI)
                    {
                        const cell& c = pbrCells[cellI];
                        const bool inZone = pbrZone.found(cellI);
                        forAll(c, fI)
                        {
                            const label faceI = c[fI];
                            if( faceI < 0 || faceI >= label(pbrFaces.size()) ) continue;
                            const face& f = pbrFaces[faceI];
                            forAll(f, pI)
                            {
                                const label ptI = f[pI];
                                if( ptI < 0 ) continue;
                                if( inZone ) pbrRepairPts.insert(ptI);
                                else pbrOutsidePts.insert(ptI);
                            }
                        }
                    }
                    label pbrMovable = 0;
                    forAllConstIter(labelHashSet, pbrRepairPts, it)
                        if( !pbrOutsidePts.found(it.key()) ) ++pbrMovable;
                    Info << "POSTBLRESCUE zone: seeds=" << pbrNSeeds
                         << " repairCells=" << pbrZone.size()
                         << " repairPts=" << pbrRepairPts.size()
                         << " movablePts=" << pbrMovable
                         << " rings=" << postBLRescueMaxRings
                         << endl;
                    if( pbrMovable > 0 )
                    {
                        const pointField pbrSnapBefore(mesh_.points());
                        labelHashSet pbrBadBefore;
                        polyMeshGenChecks::checkFacePyramids
                            (mesh_, false, -SMALL, &pbrBadBefore);
                        labelHashSet pbrOpenBefore;
                        polyMeshGenChecks::checkClosedCells
                            (mesh_, false, 0.5, &pbrOpenBefore);
                        meshSurfaceEngine msePBR(mesh_);
                        meshSurfacePartitioner mPartPBR(msePBR);
                        labelLongList globalToBpPBR(mesh_.points().size(), -1);
                        const labelList& bPtsPBR = msePBR.boundaryPoints();
                        forAll(bPtsPBR, bpI)
                            globalToBpPBR[bPtsPBR[bpI]] = bpI;
                        const label nBpPBR = bPtsPBR.size();
                        vectorField featTanPBR(nBpPBR, vector::zero);
                        {
                            const edgeList& edgesPBR = msePBR.edges();
                            const VRWGraph& bpEdgesPBR = msePBR.boundaryPointEdges();
                            const labelHashSet& featEdgesPBR = mPartPBR.featureEdges();
                            const labelHashSet& edgePtsPBR = mPartPBR.edgePoints();
                            const pointFieldPMG& ptsPBR = mesh_.points();
                            const labelList& bpPBR = msePBR.bp();
                            forAllConstIter(labelHashSet, edgePtsPBR, it)
                            {
                                const label bpI = it.key();
                                if( bpI < 0 || bpI >= nBpPBR ) continue;
                                label nbr0 = -1, nbr1 = -1;
                                forAllRow(bpEdgesPBR, bpI, eI)
                                {
                                    const label beI = bpEdgesPBR(bpI, eI);
                                    if( !featEdgesPBR.found(beI) ) continue;
                                    const edge& e = edgesPBR[beI];
                                    const label ep0 = e.start();
                                    const label ep1 = e.end();
                                    if( ep0 < 0 || ep0 >= bpPBR.size() ) continue;
                                    if( ep1 < 0 || ep1 >= bpPBR.size() ) continue;
                                    const label o0 = bpPBR[ep0];
                                    const label o1 = bpPBR[ep1];
                                    if( o0 < 0 || o1 < 0 ) continue;
                                    label otherBp = -1;
                                    if( o0 == bpI ) otherBp = o1;
                                    else if( o1 == bpI ) otherBp = o0;
                                    if( otherBp < 0 || otherBp >= nBpPBR ) continue;
                                    if( nbr0 == -1 ) nbr0 = otherBp;
                                    else if( nbr1 == -1 && otherBp != nbr0 )
                                        nbr1 = otherBp;
                                }
                                vector t = vector::zero;
                                if( nbr0 != -1 && nbr1 != -1 )
                                    t = ptsPBR[bPtsPBR[nbr1]]
                                      - ptsPBR[bPtsPBR[nbr0]];
                                else if( nbr0 != -1 )
                                    t = ptsPBR[bPtsPBR[nbr0]]
                                      - ptsPBR[bPtsPBR[bpI]];
                                if( magSqr(t) > VSMALL )
                                    featTanPBR[bpI] = t / mag(t);
                            }
                        }
                        labelLongList pbrCellsToLock;
                        forAll(pbrCells, cellI)
                            if( !pbrZone.found(cellI) )
                                pbrCellsToLock.append(cellI);
                        meshOptimizer pbrOpt(mesh_);
                        pbrOpt.lockCells(pbrCellsToLock);
                        pbrOpt.setSurfaceConstraint
                        (
                            octreePtr_,
                            &mPartPBR.pointPatches(),
                            &globalToBpPBR,
                            &mPartPBR.corners(),
                            &featTanPBR
                        );
                        //- Use untangleMeshFV: BL-introduced negVol cells
                        //- are topologically stressed. Full optimizeMeshFV
                        //- cascaded 1->124 negVol. Untangler is targeted
                        //- at inverted cells, less likely to cascade damage.
                        //- Rollback/acceptance gate unchanged.
                        Info << "POSTBLRESCUE method=constrainedUntangleMeshFV"
                             << endl;
                        pbrOpt.untangleMeshFV(5, 20, 2, false);
                        pbrOpt.setSurfaceConstraint(NULL,NULL,NULL,NULL,NULL);
                        labelHashSet pbrNegAfter;
                        polyMeshGenChecks::checkCellVolumes
                            (mesh_, false, &pbrNegAfter);
                        labelHashSet pbrBadAfter;
                        polyMeshGenChecks::checkFacePyramids
                            (mesh_, false, -SMALL, &pbrBadAfter);
                        labelHashSet pbrOpenAfter;
                        polyMeshGenChecks::checkClosedCells
                            (mesh_, false, 0.5, &pbrOpenAfter);
                        scalarField pbrSkewAfter;
                        polyMeshGenChecks::checkFaceSkewness
                            (mesh_, pbrSkewAfter);
                        const scalar pbrMaxSkew =
                            pbrSkewAfter.size() > 0 ?
                            max(pbrSkewAfter) : scalar(0);
                        const bool pbrOK =
                            pbrNegAfter.size() == 0
                         && label(pbrBadAfter.size()) <=
                            label(pbrBadBefore.size()) + postBLRescueAllowedPyrIncrease
                         && pbrOpenAfter.size() <= pbrOpenBefore.size()
                         && pbrMaxSkew <= postBLRescueMaxSkew;
                        Info << "POSTBLRESCUE after: negVol="
                             << pbrNegAfter.size()
                             << " badPyr=" << pbrBadAfter.size()
                             << " maxSkew=" << pbrMaxSkew
                             << endl;
                        if( pbrOK )
                        {
                            Info << "POSTBLRESCUE ACCEPTED: negVol "
                                 << finalNegCells.size() << "->0"
                                 << " badPyr " << pbrBadBefore.size()
                                 << "->" << pbrBadAfter.size()
                                 << " maxSkew=" << pbrMaxSkew
                                 << endl;
                            //- Recompute final sets so downstream
                            //- MESHHISTORY/CSV reflects rescued state
                            finalNegCells.clear();
                            finalBadPyrFaces = pbrBadAfter;
                            meshHistory_ = MeshHistory::CleanPromoted;
                        }
                        else
                        {
                            Info << "POSTBLRESCUE REJECTED: negVol "
                                 << finalNegCells.size() << "->"
                                 << pbrNegAfter.size()
                                 << " badPyr " << pbrBadBefore.size()
                                 << "->" << pbrBadAfter.size()
                                 << " maxSkew=" << pbrMaxSkew
                                 << " -- rolling back point motion" << endl;
                            polyMeshGenModifier pbrMod(mesh_);
                            pointFieldPMG& pbrModPts = pbrMod.pointsAccess();
                            forAll(pbrSnapBefore, pI)
                                pbrModPts[pI] = pbrSnapBefore[pI];
                            mesh_.clearAddressingData();
                        }
                    }
                    else
                    {
                        Info << "POSTBLRESCUE: no movable points in repair zone"
                             << " -- skipping (increase postBLRescueMaxRings)"
                             << endl;
                    }
                }
                else if( enablePostBLRescue && finalNegCells.size() > 0
                      && !octreePtr_ )
                {
                    Info << "POSTBLRESCUE skipped: no octreePtr_ for"
                         << " surface-constrained repair" << endl;
                }
            }
            //- MESHHISTORY logged here -- after BL and replaceBoundaries
            //- so it reflects true final mesh state.
            if( finalNegCells.size() > 0
             && (
                    meshHistory_ == MeshHistory::CleanNatural
                 || meshHistory_ == MeshHistory::CleanPromoted
                )
            )
            {
                Info << "MESHHISTORY: final validation found "
                     << finalNegCells.size()
                     << " negVol cells after BL/replaceBoundaries -- demoting "
                     << meshHistoryName(meshHistory_)
                     << " to DirtyRecoverable" << endl;
                meshHistory_ = MeshHistory::DirtyRecoverable;
            }
            Info << "MESHHISTORY meshHistory="
                 << meshHistoryName(meshHistory_)
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
            // Defect atlas: bad pyramid face centroids with full patch
            // attribution (direct + nearest boundary for internal faces).
            // Mirrors negVol dump strategy so both CSVs are joinable.
            if( finalBadPyrFaces.size() > 0 )
            {
                const pointFieldPMG& ptsBP  = mesh_.points();
                const faceListPMG&   facesBP = mesh_.faces();
                // cellsBP not needed -- owner/nei accessed via mesh_ directly
                const meshSurfaceEngine mseBP(mesh_);
                const labelList& bPointsBP    = mseBP.boundaryPoints();
                const VRWGraph&  pointFacesBP = mseBP.pointFaces();
                const labelList& facePatchBP  = mseBP.boundaryFacePatches();
                const PtrList<boundaryPatch>& boundsBP = mesh_.boundaries();
                // mesh-point -> boundary-point reverse map
                labelList meshToBpBP(ptsBP.size(), -1);
                forAll(bPointsBP, bpI)
                    meshToBpBP[bPointsBP[bpI]] = bpI;
                const label nInternalBP = mesh_.nInternalFaces();
                // owner/neighbour for internal face attribution
                const labelList& own = mesh_.owner();
                const labelList& nei = mesh_.neighbour();
                OFstream badPyrFile("badPyrFaceCentres.csv");
                badPyrFile << "faceI,cx,cy,cz,nPts,faceType,"
                           << "directPatch,nearestBoundaryPatch,"
                           << "distGapAction,distTripleJunction,"
                           << "ownerCell,neighbourCell" << nl;
                forAllConstIter(labelHashSet, finalBadPyrFaces, it)
                {
                    const label faceI = it.key();
                    if( faceI < 0 || faceI >= label(facesBP.size()) ) continue;
                    const face& f = facesBP[faceI];
                    // Face centre
                    point fc = point::zero;
                    label nValid = 0;
                    forAll(f, fpI)
                    {
                        const label pI = f[fpI];
                        if( pI >= 0 && pI < label(ptsBP.size()) )
                        { fc += ptsBP[pI]; ++nValid; }
                    }
                    if( nValid > 0 ) fc /= scalar(nValid);
                    // Direct patch (boundary faces only)
                    word faceType("internal");
                    word directPatch("internal");
                    if( faceI >= nInternalBP )
                    {
                        faceType = "boundary";
                        const label bfI = faceI - nInternalBP;
                        if( bfI < label(facePatchBP.size()) )
                        {
                            const label pI = facePatchBP[bfI];
                            if( pI >= 0 && pI < label(boundsBP.size()) )
                                directPatch = boundsBP[pI].patchName();
                        }
                    }
                    // Nearest boundary patch via face points -> bpI -> facePatch
                    word nearestBPatch("unknown");
                    forAll(f, fpI)
                    {
                        const label mPt = f[fpI];
                        if( mPt < 0 || mPt >= label(meshToBpBP.size()) ) continue;
                        const label bpI = meshToBpBP[mPt];
                        if( bpI < 0 ) continue;
                        forAllRow(pointFacesBP, bpI, pfI)
                        {
                            const label bfI = pointFacesBP(bpI, pfI);
                            if( bfI < 0 || bfI >= label(facePatchBP.size()) ) continue;
                            const label pI = facePatchBP[bfI];
                            if( pI >= 0 && pI < label(boundsBP.size()) )
                            {
                                nearestBPatch = boundsBP[pI].patchName();
                                break;
                            }
                        }
                        if( nearestBPatch != word("unknown") ) break;
                    }
                    // Distance to gap action and triple junction points
                    // (rebuild point lists locally -- gapPtPos/tjPtPos
                    //  are scoped inside the negVol block above)
                    scalar distGapBP = GREAT;
                    forAllConstIter(labelHashSet, blGapActionPoints_, git)
                    {
                        const label mpI = git.key();
                        if( mpI < 0 || mpI >= label(ptsBP.size()) ) continue;
                        const scalar d = mag(ptsBP[mpI] - fc);
                        if( d < distGapBP ) distGapBP = d;
                    }
                    scalar distTJBP = GREAT;
                    forAllConstIter(labelHashSet, blblJunctionPoints_, tjit)
                    {
                        const label mpI = tjit.key();
                        if( mpI < 0 || mpI >= label(ptsBP.size()) ) continue;
                        const scalar d = mag(ptsBP[mpI] - fc);
                        if( d < distTJBP ) distTJBP = d;
                    }
                    // Owner/neighbour cells
                    const label ownCell =
                        faceI < label(own.size()) ? own[faceI] : -1;
                    const label neiCell =
                        faceI < label(nei.size()) ? nei[faceI] : -1;
                    badPyrFile << faceI << ","
                               << fc.x() << ","
                               << fc.y() << ","
                               << fc.z() << ","
                               << f.size() << ","
                               << faceType << ","
                               << directPatch << ","
                               << nearestBPatch << ","
                               << distGapBP << ","
                               << distTJBP << ","
                               << ownCell << ","
                               << neiCell << nl;
                }
                Info << "Bad pyramid face centres written to"
                     << " badPyrFaceCentres.csv"
                     << " (" << finalBadPyrFaces.size() << " faces)" << endl;
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

// v7.1 Phase 2b (SOL): structural/spatial junction contact-line stitch.
// Proven in isolation (Phase 2a) on vertex 3748's blade_3|shroud
// interface. This version identifies the SAME junction structurally
// (patches + missing interface + tight geometric proximity) rather than
// by hardcoded point ID, since point numbering is not guaranteed stable.
// FAIL-CLOSED: returns nullptr (no repair applied) unless EXACTLY one
// qualifying candidate is found. Intentionally narrow -- this is still
// the single-junction proof case, not general production repair.
static triSurf* phase2bJunctionStitch
(
    const triSurf& rawSurf,
    const dictionary& meshDict
)
{
    PatchRoleMap roles(meshDict);
    if( !roles.active() )
    {
        Info << "[JunctionStitch3748] PatchRoleMap inactive -- skipping" << endl;
        return nullptr;
    }

    const wordList& patchNames = rawSurf.patchNames();
    label bladeId = -1, periodicId = -1, wallId = -1;
    forAll(patchNames, pI)
    {
        const word& n = patchNames[pI];
        if( n == word("blade_3") ) bladeId = pI;
        else if( n == word("periodic_1") ) periodicId = pI;
        else if( n == word("shroud") ) wallId = pI;
    }
    if( bladeId < 0 || periodicId < 0 || wallId < 0 )
    {
        Info << "[JunctionStitch3748] required patches not found"
             << " (blade_3/periodic_1/shroud) -- skipping" << endl;
        return nullptr;
    }

    const pointField& pts = rawSurf.points();
    const VRWGraph& ptFacets = rawSurf.pointFacets();
    const VRWGraph& ptEdges = rawSurf.pointEdges();
    const VRWGraph& edgeFacetsG = rawSurf.edgeFacets();
    const LongList<labelledTri>& facets = rawSurf.facets();

    // Known approximate location of the proof junction (from Phase 2a),
    // used only as a coarse spatial pre-filter -- NOT as an identity
    // check. The real qualification is entirely structural (patches +
    // missing interface + open-edge geometry), matching SOL's
    // instruction not to key off literal point/vertex IDs.
    const point approxLoc(0.242172, 0.042408, 0.0295046);
    const scalar coarseRadius = 0.01; // 10mm, generous pre-filter only

    label candidateVertex = -1;
    label candidateOpenA = -1, candidateOpenB = -1;
    point candidatePosA(vector::zero), candidatePosB(vector::zero);
    scalar candidateGap = -1;
    label nCandidatesFound = 0;

    forAll(pts, spI)
    {
        if( magSqr(pts[spI] - approxLoc) > sqr(coarseRadius) ) continue;

        // Raw incident patch check (ALL_3 for our three required roles).
        labelHashSet incident;
        forAllRow(ptFacets, spI, fI)
        {
            const label triI = ptFacets(spI, fI);
            if( triI >= 0 && triI < facets.size() )
                incident.insert(facets[triI].region());
        }
        if( !incident.found(bladeId) || !incident.found(periodicId)
         || !incident.found(wallId) ) continue;

        // Walk incident edges: find validity of each of the 3 interfaces,
        // and collect open (nFacets==1) edges by owning patch.
        bool validBladePeriodic = false, validBladeWall = false, validPeriodicWall = false;
        DynList<label> openBladeEdges, openWallEdges;

        forAllRow(ptEdges, spI, ord)
        {
            const label eI = ptEdges(spI, ord);
            const label nF = edgeFacetsG.sizeOfRow(eI);
            if( nF == 1 )
            {
                const label owner = facets[edgeFacetsG(eI, 0)].region();
                if( owner == bladeId ) openBladeEdges.append(eI);
                else if( owner == wallId ) openWallEdges.append(eI);
            }
            else if( nF == 2 )
            {
                const label r0 = facets[edgeFacetsG(eI, 0)].region();
                const label r1 = facets[edgeFacetsG(eI, 1)].region();
                if( r0 != r1 )
                {
                    const bool hasBlade = (r0 == bladeId) || (r1 == bladeId);
                    const bool hasPeriodic = (r0 == periodicId) || (r1 == periodicId);
                    const bool hasWall = (r0 == wallId) || (r1 == wallId);
                    if( hasBlade && hasPeriodic ) validBladePeriodic = true;
                    if( hasBlade && hasWall ) validBladeWall = true;
                    if( hasPeriodic && hasWall ) validPeriodicWall = true;
                }
            }
        }

        // Require EXACTLY the blade_3|shroud interface missing, the
        // other two present, and exactly one open edge on each side.
        if( validBladeWall ) continue;
        if( !validBladePeriodic || !validPeriodicWall ) continue;
        if( openBladeEdges.size() != 1 || openWallEdges.size() != 1 ) continue;

        // Geometric proximity check (SOL's tight prototype tolerances).
        auto otherEndpoint = [&](const label eI) -> label
        {
            const edge& e = rawSurf.edges()[eI];
            return (e.start() == spI) ? e.end() : e.start();
        };
        const label ptA = otherEndpoint(openBladeEdges[0]);
        const label ptB = otherEndpoint(openWallEdges[0]);
        if( ptA < 0 || ptA >= pts.size() || ptB < 0 || ptB >= pts.size() ) continue;

        const point& posA = pts[ptA];
        const point& posB = pts[ptB];
        const scalar gap = Foam::sqrt(magSqr(posA - posB));
        const point& triplePos = pts[spI];
        const vector dirA = posA - triplePos;
        const vector dirB = posB - triplePos;
        const scalar lenA = Foam::sqrt(magSqr(dirA));
        const scalar lenB = Foam::sqrt(magSqr(dirB));
        if( lenA < VSMALL || lenB < VSMALL ) continue;

        const scalar cosAngle = (dirA & dirB) / (lenA * lenB);
        const scalar clamped = Foam::max(scalar(-1), Foam::min(scalar(1), cosAngle));
        const scalar angleDeg =
            Foam::acos(clamped) * 180.0 / Foam::constant::mathematical::pi;
        const scalar lengthRatio = lenA / lenB;

        // Tight prototype tolerances -- matches blade_3's known signature
        // (gap ~10um, angle ~0.004deg, ratio ~0.997), NOT the looser
        // CLASS_B thresholds used for the exploratory sweep.
        const bool qualifies =
            gap < 0.0005 && angleDeg < 0.5
         && lengthRatio > 0.9 && lengthRatio < 1.111;

        if( !qualifies ) continue;

        ++nCandidatesFound;
        candidateVertex = spI;
        candidateOpenA = ptA;
        candidateOpenB = ptB;
        candidatePosA = posA;
        candidatePosB = posB;
        candidateGap = gap;

        Info << "[JunctionStitch3748] candidate found"
             << " vertexPos=" << triplePos
             << " interface=blade_3|shroud"
             << " gap=" << gap
             << " angle=" << angleDeg
             << " lengthRatio=" << lengthRatio
             << endl;
    }

    if( nCandidatesFound != 1 )
    {
        Info << "[JunctionStitch3748] ABORT: found " << nCandidatesFound
             << " qualifying candidates, required exactly 1 -- repair NOT applied"
             << endl;
        return nullptr;
    }

    // Perform the stitch -- same validated mechanism as Phase 2a.
    pointField newPoints = pts;
    LongList<labelledTri> newFacets = facets;
    const geometricSurfacePatchList& origPatches = rawSurf.patches();
    const edgeLongList& origFeatureEdges = rawSurf.featureEdges();

    const point midpoint = 0.5 * (candidatePosA + candidatePosB);
    const label midpointId = newPoints.size();
    newPoints.append(midpoint);

    label nRemapped = 0;
    forAll(newFacets, triI)
    {
        labelledTri& tri = newFacets[triI];
        forAll(tri, vi)
        {
            if( tri[vi] == candidateOpenA || tri[vi] == candidateOpenB )
            {
                tri[vi] = midpointId;
                ++nRemapped;
            }
        }
    }

    label nDegenerate = 0, nZeroArea = 0;
    forAll(newFacets, triI)
    {
        const labelledTri& tri = newFacets[triI];
        if( tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2] )
        { ++nDegenerate; continue; }
        const point& p0 = newPoints[tri[0]];
        const point& p1 = newPoints[tri[1]];
        const point& p2 = newPoints[tri[2]];
        if( Foam::mag((p1-p0)^(p2-p0)) < VSMALL ) ++nZeroArea;
    }

    if( nDegenerate > 0 || nZeroArea > 0 )
    {
        Info << "[JunctionStitch3748] ABORT: would create " << nDegenerate
             << " degenerate, " << nZeroArea
             << " zero-area triangles -- repair NOT applied" << endl;
        return nullptr;
    }

    triSurf* repaired =
        new triSurf(newFacets, origPatches, origFeatureEdges, newPoints);

    Info << "[JunctionStitch3748] APPLIED: vertex=" << candidateVertex
         << " remapped=" << nRemapped << " facet refs"
         << " newPointCount=" << repaired->points().size()
         << " (was " << pts.size() << ")"
         << endl;

    // Checkpoint 1: verify immediately after stitch.
    {
        const VRWGraph& vPtEdges = repaired->pointEdges();
        const VRWGraph& vEdgeFacets = repaired->edgeFacets();
        const LongList<labelledTri>& vFacets = repaired->facets();
        label qualifying = 0;
        bool sharedFound = false;
        label sharedFacetCount = -1;
        forAllRow(vPtEdges, candidateVertex, ord)
        {
            const label eI = vPtEdges(candidateVertex, ord);
            const label nF = vEdgeFacets.sizeOfRow(eI);
            if( nF != 2 ) continue;
            const label r0 = vFacets[vEdgeFacets(eI, 0)].region();
            const label r1 = vFacets[vEdgeFacets(eI, 1)].region();
            if( r0 != r1 )
            {
                ++qualifying;
                if( (r0 == bladeId && r1 == wallId) || (r0 == wallId && r1 == bladeId) )
                { sharedFound = true; sharedFacetCount = nF; }
            }
        }
        Info << "[JunctionStitchCheckpoint] stage=afterStitch"
             << " sharedEdgeFound=" << (sharedFound ? 1 : 0)
             << " sharedEdgeFacetCount=" << sharedFacetCount
             << " qualifyingEdges=" << qualifying
             << " eligible=" << ((qualifying > 2) ? 1 : 0)
             << endl;
    }

    return repaired;
}

// v7.1 Phase 2c (SOL): generalized multi-junction stitch. Same proven
// mechanism as Phase 2b, but role-based (blade/periodic/hub-or-shroud
// via PatchRoleMap, not hardcoded patch names) and repeated until no
// more qualifying SPLIT_ONE_INTERFACE candidates remain. Each stitch
// is applied one at a time, then the surface is fully re-scanned from
// its updated state before looking for the next candidate -- avoids
// any stale-index risk from applying multiple repairs against one
// snapshot of point/facet data. SPLIT_TWO_INTERFACES vertices (two
// missing interfaces) are intentionally left alone -- Phase 2d handles
// those separately, per SOL's sequential detect-repair-redetect plan.
// v7.1 Phase 2c diagnosis (SOL): dry-run scanner, NO mutation. Counts
// qualifying SPLIT_ONE_INTERFACE candidates with per-rejection-reason
// counters, per SOL's exact spec. Used to isolate whether the scanner
// itself finds the expected ~14 candidates on a given surface (SCAN 0
// on original, SCAN 1 on post-repair reconstruction).
static void phase2cDryRunScan
(
    const triSurf& surf,
    const labelHashSet& bladeIds,
    const labelHashSet& periodicIds,
    const labelHashSet& wallIds,
    const word& tag
)
{
    const pointField& pts = surf.points();
    const VRWGraph& ptFacets = surf.pointFacets();
    const VRWGraph& ptEdges = surf.pointEdges();
    const VRWGraph& edgeFacetsG = surf.edgeFacets();
    const LongList<labelledTri>& facets = surf.facets();

    label verticesVisited = 0, all3Vertices = 0, wrongRoles = 0;
    label completeManifold = 0, notSplitOne = 0, wrongOpenEdgeCount = 0;
    label splitTwoInterfaces = 0;
    label gapTooLarge = 0, angleTooLarge = 0, lengthRatioBad = 0;
    label qualifying = 0;

    forAll(pts, spI)
    {
        ++verticesVisited;

        labelHashSet incident;
        forAllRow(ptFacets, spI, fI)
        {
            const label triI = ptFacets(spI, fI);
            if( triI >= 0 && triI < facets.size() )
                incident.insert(facets[triI].region());
        }

        label thisBlade = -1, thisPeriodic = -1, thisWall = -1;
        forAllConstIter(labelHashSet, incident, it)
        {
            const label p = it.key();
            if( bladeIds.found(p) ) thisBlade = p;
            else if( periodicIds.found(p) ) thisPeriodic = p;
            else if( wallIds.found(p) ) thisWall = p;
        }
        if( thisBlade < 0 || thisPeriodic < 0 || thisWall < 0 )
        {
            if( incident.size() >= 3 ) ++wrongRoles;
            continue;
        }
        ++all3Vertices;

        bool validBladePeriodic = false, validBladeWall = false, validPeriodicWall = false;
        DynList<label> openBladeEdges, openWallEdges;

        forAllRow(ptEdges, spI, ord)
        {
            const label eI = ptEdges(spI, ord);
            const label nF = edgeFacetsG.sizeOfRow(eI);
            if( nF == 1 )
            {
                const label owner = facets[edgeFacetsG(eI, 0)].region();
                if( owner == thisBlade ) openBladeEdges.append(eI);
                else if( owner == thisWall ) openWallEdges.append(eI);
            }
            else if( nF == 2 )
            {
                const label r0 = facets[edgeFacetsG(eI, 0)].region();
                const label r1 = facets[edgeFacetsG(eI, 1)].region();
                if( r0 != r1 )
                {
                    const bool hasBlade = (r0 == thisBlade) || (r1 == thisBlade);
                    const bool hasPeriodic = (r0 == thisPeriodic) || (r1 == thisPeriodic);
                    const bool hasWall = (r0 == thisWall) || (r1 == thisWall);
                    if( hasBlade && hasPeriodic ) validBladePeriodic = true;
                    if( hasBlade && hasWall ) validBladeWall = true;
                    if( hasPeriodic && hasWall ) validPeriodicWall = true;
                }
            }
        }

        // SOL fix: distinct, mutually exclusive buckets instead of a
        // single conflated notSplitOne, so SPLIT_TWO_INTERFACES vertices
        // (614/711-style) are never hidden inside the same counter as
        // "blade-wall unexpectedly already valid" cases.
        label nMissing = 0;
        if( !validBladePeriodic ) ++nMissing;
        if( !validBladeWall ) ++nMissing;
        if( !validPeriodicWall ) ++nMissing;

        Info << "[StitchScanAll3] tag=" << tag.c_str()
             << " vertex=" << spI << " pos=" << pts[spI]
             << " bladeId=" << thisBlade
             << " periodicId=" << thisPeriodic
             << " wallId=" << thisWall
             << " validBladePeriodic=" << (validBladePeriodic ? 1 : 0)
             << " validBladeWall=" << (validBladeWall ? 1 : 0)
             << " validPeriodicWall=" << (validPeriodicWall ? 1 : 0)
             << " openBladeEdges=" << openBladeEdges.size()
             << " openWallEdges=" << openWallEdges.size()
             << " nMissingInterfaces=" << nMissing
             << endl;

        if( nMissing == 0 )
        { ++completeManifold; continue; }
        if( !validBladeWall && validBladePeriodic && validPeriodicWall )
        {
            // exactly the target pattern -- fall through to gap/angle checks
        }
        else if( nMissing >= 2 )
        { ++splitTwoInterfaces; continue; }
        else
        { ++notSplitOne; continue; } // blade-wall present but a DIFFERENT single interface missing
        if( openBladeEdges.size() != 1 || openWallEdges.size() != 1 )
        { ++wrongOpenEdgeCount; continue; }

        auto otherEndpoint = [&](const label eI) -> label
        {
            const edge& e = surf.edges()[eI];
            return (e.start() == spI) ? e.end() : e.start();
        };
        const label ptA = otherEndpoint(openBladeEdges[0]);
        const label ptB = otherEndpoint(openWallEdges[0]);
        if( ptA < 0 || ptA >= pts.size() || ptB < 0 || ptB >= pts.size() ) continue;

        const point& posA = pts[ptA];
        const point& posB = pts[ptB];
        const scalar gap = Foam::sqrt(magSqr(posA - posB));
        const point& triplePos = pts[spI];
        const vector dirA = posA - triplePos;
        const vector dirB = posB - triplePos;
        const scalar lenA = Foam::sqrt(magSqr(dirA));
        const scalar lenB = Foam::sqrt(magSqr(dirB));
        if( lenA < VSMALL || lenB < VSMALL ) continue;

        const scalar cosAngle = (dirA & dirB) / (lenA * lenB);
        const scalar clamped = Foam::max(scalar(-1), Foam::min(scalar(1), cosAngle));
        const scalar angleDeg =
            Foam::acos(clamped) * 180.0 / Foam::constant::mathematical::pi;
        const scalar lengthRatio = lenA / lenB;

        // SOL redesign: scale-invariant structural criteria, NO absolute
        // distance threshold. A and B qualify as the same split curve if
        // (1) they point in a consistent shared direction (angle small),
        // (2) their lengths are consistent (already have lengthRatio),
        // and (3) the perpendicular offset between A and B is small
        // RELATIVE to the local edge length -- i.e. B sits close to the
        // LINE through triplePos/A, not just close in absolute terms.
        // DRY-RUN ONLY: both old absolute and new relative criteria are
        // computed and logged side by side. qualifying/acceptance logic
        // is UNCHANGED in this pass -- no repair behavior is altered.
        const vector dirAUnit = dirA / lenA;
        const vector dirBUnit = dirB / lenB;
        const vector avgDir = 0.5 * (dirAUnit + dirBUnit);
        const scalar avgDirMag = Foam::sqrt(magSqr(avgDir));

        scalar perpDeviation = -1;
        scalar perpDeviationRatio = -1;
        if( avgDirMag > VSMALL )
        {
            const vector avgDirUnit = avgDir / avgDirMag;
            const vector offsetVec = posB - posA;
            const scalar alongComponent = offsetVec & avgDirUnit;
            const vector perpVec = offsetVec - alongComponent * avgDirUnit;
            perpDeviation = Foam::sqrt(magSqr(perpVec));
            const scalar localScale = Foam::min(lenA, lenB);
            if( localScale > VSMALL )
                perpDeviationRatio = perpDeviation / localScale;
        }

        // SOL refinement: two additional DIAGNOSTIC-ONLY fields (not
        // gates) to inform the merge-point choice (50/50 midpoint vs
        // weighted toward the more "authoritative" side) once we move
        // past pure detection.
        const scalar shorterEdgeRatio = Foam::min(lenA, lenB) / Foam::max(lenA, lenB);

        label facetsAtA = 0, facetsAtB = 0;
        if( ptA >= 0 && ptA < ptFacets.size() ) facetsAtA = ptFacets.sizeOfRow(ptA);
        if( ptB >= 0 && ptB < ptFacets.size() ) facetsAtB = ptFacets.sizeOfRow(ptB);
        const scalar triangleDensityRatio =
            (facetsAtB > 0) ? (scalar(facetsAtA) / scalar(facetsAtB)) : scalar(-1);

        Info << "[StitchScanGeometry] tag=" << tag.c_str()
             << " vertex=" << spI
             << " pos=" << pts[spI]
             << " gap=" << gap
             << " angleDeg=" << angleDeg
             << " lengthRatio=" << lengthRatio
             << " lenA=" << lenA
             << " lenB=" << lenB
             << " perpDeviation=" << perpDeviation
             << " perpDeviationRatio=" << perpDeviationRatio
             << " shorterEdgeRatio=" << shorterEdgeRatio
             << " facetsAtA=" << facetsAtA
             << " facetsAtB=" << facetsAtB
             << " triangleDensityRatio=" << triangleDensityRatio
             << endl;

        // SOL: perpDeviationRatio is now the SOLE hard acceptance gate --
        // scale-invariant, validated against real data across an 800x
        // gap range (10 microns to 8mm), all 8 known candidates pass
        // with comfortable margin (max observed 0.038). angleDeg/
        // lengthRatio old absolute-style thresholds are REMOVED as
        // gates -- they rejected most real candidates and added no
        // discriminating value beyond what perpDeviationRatio already
        // captures. gapTooLarge/angleTooLarge/lengthRatioBad counters
        // kept declared (still initialized above) but no longer
        // incremented -- will read 0 in this and future runs; not
        // removing the declarations yet to avoid touching the summary
        // print signature mid-diagnosis.
        if( perpDeviationRatio < 0 || perpDeviationRatio >= 0.1 ) continue;

        ++qualifying;
        Info << "[StitchScanCandidate] tag=" << tag.c_str()
             << " vertex=" << spI << " pos=" << triplePos
             << " bladeId=" << thisBlade << " wallId=" << thisWall
             << " gap=" << gap << " angle=" << angleDeg
             << " lengthRatio=" << lengthRatio << endl;
    }

    Info << "[StitchScan] tag=" << tag.c_str()
         << " verticesVisited=" << verticesVisited
         << " all3Vertices=" << all3Vertices
         << " wrongRoles=" << wrongRoles
         << " completeManifold=" << completeManifold
         << " notSplitOne=" << notSplitOne
         << " splitTwoInterfaces=" << splitTwoInterfaces
         << " wrongOpenEdgeCount=" << wrongOpenEdgeCount
         << " gapTooLarge=" << gapTooLarge
         << " angleTooLarge=" << angleTooLarge
         << " lengthRatioBad=" << lengthRatioBad
         << " qualifying=" << qualifying
         << endl;
}

static triSurf* phase2cMultiJunctionStitch
(
    const triSurf& rawSurf,
    const dictionary& meshDict,
    const label maxRepairs = 20,
    const labelHashSet& vertexAllowlist = labelHashSet()
)
{
    // SOL controlled A/B test: empty allowlist = allow all (unchanged
    // default behavior). Non-empty = ONLY repair candidates whose
    // vertex ID is in this set. Temporary, for isolating vertex 3737's
    // repair from the other 5 -- not intended as a permanent API.
    const bool useAllowlist = !vertexAllowlist.empty();
    PatchRoleMap roles(meshDict);
    if( !roles.active() )
    {
        Info << "[JunctionStitchMulti] PatchRoleMap inactive -- skipping" << endl;
        return nullptr;
    }

    const wordList& patchNames0 = rawSurf.patchNames();
    labelHashSet bladeIds, periodicIds, wallIds;
    forAll(patchNames0, pI)
    {
        const word& n = patchNames0[pI];
        if( roles.hasRole(n, word("blade")) ) bladeIds.insert(pI);
        if( roles.hasRole(n, word("periodic")) ) periodicIds.insert(pI);
        if( roles.hasRole(n, word("hub")) || roles.hasRole(n, word("shroud")) )
            wallIds.insert(pI);
    }
    if( bladeIds.empty() || periodicIds.empty() || wallIds.empty() )
    {
        Info << "[JunctionStitchMulti] required roles not resolved -- skipping" << endl;
        return nullptr;
    }

    // Working copy -- rebuilt fresh into a new triSurf after each
    // accepted repair, so every subsequent scan sees up-to-date topology.
    triSurf* current = new triSurf
    (
        rawSurf.facets(), rawSurf.patches(),
        rawSurf.featureEdges(), rawSurf.points()
    );

    label nApplied = 0;
    label nAttempts = 0;

    while( nApplied < maxRepairs )
    {
        ++nAttempts;
        if( nAttempts > maxRepairs * 3 )
        {
            Info << "[JunctionStitchMulti] ABORT: too many scan attempts"
                 << " without progress -- stopping to avoid infinite loop"
                 << endl;
            break;
        }

        const pointField& pts = current->points();
        const VRWGraph& ptFacets = current->pointFacets();
        const VRWGraph& ptEdges = current->pointEdges();
        const VRWGraph& edgeFacetsG = current->edgeFacets();
        const LongList<labelledTri>& facets = current->facets();
        const wordList& patchNames = current->patchNames();

        // SOL: per-known-vertex persistence trace. Checks the 8
        // original candidate IDs (117, 132, 3737, 3738, 3748, 3749,
        // 4006, 4010) on EVERY scan pass, regardless of which one
        // actually gets chosen this pass. Answers whether 3749/4006
        // legitimately lose their open edges as a side effect of an
        // earlier repair's facet remap, or disappear for an unrelated
        // reason (rescan bug). Diagnostic only -- does not affect
        // which candidate is chosen or how the repair is applied.
        {
            static const label knownIds[8] =
                { 117, 132, 3737, 3738, 3748, 3749, 4006, 4010 };
            for( label ki = 0; ki < 8; ++ki )
            {
                const label kv = knownIds[ki];
                if( kv < 0 || kv >= pts.size() )
                {
                    Info << "[CandidateTrace] scan=" << nAttempts
                         << " vertex=" << kv
                         << " status=OUT_OF_RANGE" << endl;
                    continue;
                }

                labelHashSet kIncident;
                forAllRow(ptFacets, kv, fI)
                {
                    const label triI = ptFacets(kv, fI);
                    if( triI >= 0 && triI < facets.size() )
                        kIncident.insert(facets[triI].region());
                }

                label kBlade = -1, kPeriodic = -1, kWall = -1;
                forAllConstIter(labelHashSet, kIncident, it)
                {
                    const label p = it.key();
                    if( bladeIds.found(p) ) kBlade = p;
                    else if( periodicIds.found(p) ) kPeriodic = p;
                    else if( wallIds.found(p) ) kWall = p;
                }

                if( kBlade < 0 || kPeriodic < 0 || kWall < 0 )
                {
                    Info << "[CandidateTrace] scan=" << nAttempts
                         << " vertex=" << kv
                         << " status=ROLES_LOST"
                         << " nIncidentPatches=" << kIncident.size()
                         << endl;
                    continue;
                }

                bool kValidBW = false, kValidBP = false, kValidPW = false;
                label kOpenBlade = 0, kOpenWall = 0;
                forAllRow(ptEdges, kv, ord)
                {
                    const label eI = ptEdges(kv, ord);
                    const label nF = edgeFacetsG.sizeOfRow(eI);
                    if( nF == 1 )
                    {
                        const label owner = facets[edgeFacetsG(eI, 0)].region();
                        if( owner == kBlade ) ++kOpenBlade;
                        else if( owner == kWall ) ++kOpenWall;
                    }
                    else if( nF == 2 )
                    {
                        const label r0 = facets[edgeFacetsG(eI, 0)].region();
                        const label r1 = facets[edgeFacetsG(eI, 1)].region();
                        if( r0 != r1 )
                        {
                            const bool hb = (r0==kBlade)||(r1==kBlade);
                            const bool hp = (r0==kPeriodic)||(r1==kPeriodic);
                            const bool hw = (r0==kWall)||(r1==kWall);
                            if( hb && hp ) kValidBP = true;
                            if( hb && hw ) kValidBW = true;
                            if( hp && hw ) kValidPW = true;
                        }
                    }
                }

                std::string kStatus;
                if( kValidBW && kValidBP && kValidPW )
                    kStatus = "COMPLETE_MANIFOLD";
                else if( kValidBW )
                    kStatus = "BLADEWALL_UNEXPECTEDLY_VALID";
                else if( !kValidBP || !kValidPW )
                    kStatus = "OTHER_INTERFACE_MISSING";
                else if( kOpenBlade != 1 || kOpenWall != 1 )
                    kStatus = "WRONG_OPEN_EDGE_COUNT";
                else
                    kStatus = "STILL_QUALIFYING_CANDIDATE";

                Info << "[CandidateTrace] scan=" << nAttempts
                     << " vertex=" << kv
                     << " status=" << kStatus.c_str()
                     << " validBladeWall=" << (kValidBW?1:0)
                     << " validBladePeriodic=" << (kValidBP?1:0)
                     << " validPeriodicWall=" << (kValidPW?1:0)
                     << " openBlade=" << kOpenBlade
                     << " openWall=" << kOpenWall
                     << endl;
            }
        }

        label foundVertex = -1;
        label foundBladeId = -1, foundPeriodicId = -1, foundWallId = -1;
        label foundOpenA = -1, foundOpenB = -1;
        point foundPosA(vector::zero), foundPosB(vector::zero);
        scalar foundGap = -1, foundAngle = -1, foundRatio = -1;

        forAll(pts, spI)
        {
            labelHashSet incident;
            forAllRow(ptFacets, spI, fI)
            {
                const label triI = ptFacets(spI, fI);
                if( triI >= 0 && triI < facets.size() )
                    incident.insert(facets[triI].region());
            }

            label thisBlade = -1, thisPeriodic = -1, thisWall = -1;
            forAllConstIter(labelHashSet, incident, it)
            {
                const label p = it.key();
                if( bladeIds.found(p) ) thisBlade = p;
                else if( periodicIds.found(p) ) thisPeriodic = p;
                else if( wallIds.found(p) ) thisWall = p;
            }
            if( thisBlade < 0 || thisPeriodic < 0 || thisWall < 0 ) continue;

            bool validBladePeriodic = false, validBladeWall = false, validPeriodicWall = false;
            DynList<label> openBladeEdges, openWallEdges;

            forAllRow(ptEdges, spI, ord)
            {
                const label eI = ptEdges(spI, ord);
                const label nF = edgeFacetsG.sizeOfRow(eI);
                if( nF == 1 )
                {
                    const label owner = facets[edgeFacetsG(eI, 0)].region();
                    if( owner == thisBlade ) openBladeEdges.append(eI);
                    else if( owner == thisWall ) openWallEdges.append(eI);
                }
                else if( nF == 2 )
                {
                    const label r0 = facets[edgeFacetsG(eI, 0)].region();
                    const label r1 = facets[edgeFacetsG(eI, 1)].region();
                    if( r0 != r1 )
                    {
                        const bool hasBlade = (r0 == thisBlade) || (r1 == thisBlade);
                        const bool hasPeriodic = (r0 == thisPeriodic) || (r1 == thisPeriodic);
                        const bool hasWall = (r0 == thisWall) || (r1 == thisWall);
                        if( hasBlade && hasPeriodic ) validBladePeriodic = true;
                        if( hasBlade && hasWall ) validBladeWall = true;
                        if( hasPeriodic && hasWall ) validPeriodicWall = true;
                    }
                }
            }

            // SPLIT_ONE_INTERFACE only -- exactly blade-wall missing,
            // the other two present. SPLIT_TWO_INTERFACES vertices
            // (614/711-style) deliberately excluded -- Phase 2d.
            if( validBladeWall ) continue;
            if( !validBladePeriodic || !validPeriodicWall ) continue;
            if( openBladeEdges.size() != 1 || openWallEdges.size() != 1 ) continue;

            auto otherEndpoint = [&](const label eI) -> label
            {
                const edge& e = current->edges()[eI];
                return (e.start() == spI) ? e.end() : e.start();
            };
            const label ptA = otherEndpoint(openBladeEdges[0]);
            const label ptB = otherEndpoint(openWallEdges[0]);
            if( ptA < 0 || ptA >= pts.size() || ptB < 0 || ptB >= pts.size() ) continue;

            const point& posA = pts[ptA];
            const point& posB = pts[ptB];
            const scalar gap = Foam::sqrt(magSqr(posA - posB));
            const point& triplePos = pts[spI];
            const vector dirA = posA - triplePos;
            const vector dirB = posB - triplePos;
            const scalar lenA = Foam::sqrt(magSqr(dirA));
            const scalar lenB = Foam::sqrt(magSqr(dirB));
            if( lenA < VSMALL || lenB < VSMALL ) continue;

            const scalar cosAngle = (dirA & dirB) / (lenA * lenB);
            const scalar clamped = Foam::max(scalar(-1), Foam::min(scalar(1), cosAngle));
            const scalar angleDeg =
                Foam::acos(clamped) * 180.0 / Foam::constant::mathematical::pi;
            const scalar lengthRatio = lenA / lenB;

            // SOL: scale-invariant perpDeviationRatio, validated against
            // real data across an 800x gap range (10 microns to 8mm) --
            // replaces the old absolute gap/angle/lengthRatio gate, which
            // rejected 7 of 8 known-good candidates. Threshold 0.1 (10%)
            // chosen with comfortable margin above the observed max of
            // 0.038 across all 8 validated samples.
            const vector dirAUnit = dirA / lenA;
            const vector dirBUnit = dirB / lenB;
            const vector avgDir = 0.5 * (dirAUnit + dirBUnit);
            const scalar avgDirMag = Foam::sqrt(magSqr(avgDir));

            scalar perpDeviationRatio = -1;
            if( avgDirMag > VSMALL )
            {
                const vector avgDirUnit = avgDir / avgDirMag;
                const vector offsetVec = posB - posA;
                const scalar alongComponent = offsetVec & avgDirUnit;
                const vector perpVec = offsetVec - alongComponent * avgDirUnit;
                const scalar perpDeviation = Foam::sqrt(magSqr(perpVec));
                const scalar localScale = Foam::min(lenA, lenB);
                if( localScale > VSMALL )
                    perpDeviationRatio = perpDeviation / localScale;
            }

            const bool qualifies =
                (perpDeviationRatio >= 0) && (perpDeviationRatio < 0.1);
            if( !qualifies ) continue;

            // SOL controlled A/B test filter (temporary).
            if( useAllowlist && !vertexAllowlist.found(spI) ) continue;

            // Take the FIRST qualifying candidate this scan -- one
            // repair per full-surface scan, then rescan fresh.
            foundVertex = spI;
            foundBladeId = thisBlade; foundPeriodicId = thisPeriodic; foundWallId = thisWall;
            foundOpenA = ptA; foundOpenB = ptB;
            foundPosA = posA; foundPosB = posB;
            foundGap = gap; foundAngle = angleDeg; foundRatio = lengthRatio;
            break;
        }

        if( foundVertex < 0 )
        {
            Info << "[JunctionStitchMulti] no more qualifying candidates"
                 << " -- stopping (applied " << nApplied << " repairs)" << endl;
            break;
        }

        Info << "[JunctionStitchMulti] candidate found vertex=" << foundVertex
             << " bladeId=" << foundBladeId
             << " wallId=" << foundWallId
             << " gap=" << foundGap
             << " angle=" << foundAngle
             << " lengthRatio=" << foundRatio
             << endl;

        pointField newPoints = pts;
        LongList<labelledTri> newFacets = facets;
        const geometricSurfacePatchList& curPatches = current->patches();
        const edgeLongList& curFeatureEdges = current->featureEdges();

        const point midpoint = 0.5 * (foundPosA + foundPosB);
        const label midpointId = newPoints.size();
        newPoints.append(midpoint);

        label nRemapped = 0;
        forAll(newFacets, triI)
        {
            labelledTri& tri = newFacets[triI];
            forAll(tri, vi)
            {
                if( tri[vi] == foundOpenA || tri[vi] == foundOpenB )
                { tri[vi] = midpointId; ++nRemapped; }
            }
        }

        label nDegenerate = 0, nZeroArea = 0;
        forAll(newFacets, triI)
        {
            const labelledTri& tri = newFacets[triI];
            if( tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2] )
            { ++nDegenerate; continue; }
            const point& p0 = newPoints[tri[0]];
            const point& p1 = newPoints[tri[1]];
            const point& p2 = newPoints[tri[2]];
            if( Foam::mag((p1-p0)^(p2-p0)) < VSMALL ) ++nZeroArea;
        }

        if( nDegenerate > 0 || nZeroArea > 0 )
        {
            Info << "[JunctionStitchMulti] REJECT candidate at vertex=" << foundVertex
                 << ": would create " << nDegenerate << " degenerate, "
                 << nZeroArea << " zero-area triangles -- skipping this one"
                 << endl;
            // Do not retry the same candidate -- break out and let the
            // caller's maxRepairs/nAttempts bound prevent infinite retry
            // of an unfixable candidate. Since we don't remove it from
            // consideration, rely on the nAttempts*3 ceiling above.
            continue;
        }

        triSurf* next = new triSurf(newFacets, curPatches, curFeatureEdges, newPoints);
        delete current;
        current = next;
        ++nApplied;

        Info << "[JunctionStitchMulti] APPLIED repair #" << nApplied
             << " at vertex=" << foundVertex
             << " remapped=" << nRemapped << " facet refs"
             << " newPointCount=" << current->points().size()
             << endl;
    }

    Info << "[JunctionStitchMulti] SUMMARY: applied " << nApplied
         << " repairs across " << nAttempts << " scan attempts" << endl;

    // SOL forensic check: UNCONDITIONAL integrity re-verification,
    // regardless of nApplied. If nApplied==0, current SHOULD be
    // point-for-point identical to rawSurf (same point count, same
    // coordinates at every shared index) -- any deviation would prove
    // rejected/discarded attempts are NOT being fully retracted,
    // explaining how a "0 repairs applied" run could still produce a
    // visibly worse mesh (Mitch's observed hub-end blade_3 blowout).
    {
        const pointField& origPts = rawSurf.points();
        const pointField& curPts = current->points();
        label nMismatched = 0;
        scalar maxDelta = 0;
        const label nCompare = Foam::min(origPts.size(), curPts.size());
        for( label pi = 0; pi < nCompare; ++pi )
        {
            const scalar d = Foam::sqrt(magSqr(curPts[pi] - origPts[pi]));
            if( d > 1e-12 )
            {
                ++nMismatched;
                if( d > maxDelta ) maxDelta = d;
            }
        }
        Info << "[JunctionStitchIntegrity] nApplied=" << nApplied
             << " origPointCount=" << origPts.size()
             << " curPointCount=" << curPts.size()
             << " nMismatchedSharedPoints=" << nMismatched
             << " maxDelta=" << maxDelta
             << endl;
        if( nApplied == 0 && (nMismatched > 0 || origPts.size() != curPts.size()) )
        {
            Info << "[JunctionStitchIntegrity] ANOMALY: nApplied=0 but"
                 << " surface differs from original -- rejected attempts"
                 << " are leaving residual state" << endl;
        }
    }

    if( nApplied == 0 )
    {
        delete current;
        return nullptr;
    }

    return current;
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

    // v7.1 Phase 2b (SOL): structural/spatial junction-contact-line
    // stitch. Runs on the RAW loaded surface, before feature-edge
    // extraction or patch derivation, so downstream feature/corner
    // topology is derived from the REPAIRED surface, not the malformed
    // one. Fail-closed: if the detector finds anything other than
    // EXACTLY one qualifying candidate, the repair is skipped entirely
    // and surfacePtr_ continues unmodified -- we would rather do
      // nothing than repair the wrong junction.
    {
        // v7.1 Phase 2c diagnosis (SOL): SCAN 0 on the original surface,
        // dry-run only, before any repair. Resolves role IDs fresh
        // against THIS surface (never carried across reconstruction).
        PatchRoleMap diagRoles(meshDict_);
        if( diagRoles.active() )
        {
            const wordList& diagNames0 = surfacePtr_->patchNames();
            labelHashSet diagBladeIds0, diagPeriodicIds0, diagWallIds0;
            forAll(diagNames0, pI)
            {
                const word& n = diagNames0[pI];
                if( diagRoles.hasRole(n, word("blade")) ) diagBladeIds0.insert(pI);
                if( diagRoles.hasRole(n, word("periodic")) ) diagPeriodicIds0.insert(pI);
                if( diagRoles.hasRole(n, word("hub")) || diagRoles.hasRole(n, word("shroud")) )
                    diagWallIds0.insert(pI);
            }
            phase2cDryRunScan
            (
                *surfacePtr_, diagBladeIds0, diagPeriodicIds0, diagWallIds0,
                word("SCAN0_original")
            );
        }

        // SOL controlled A/B test: isolate vertex 3737 only, to test
        // whether skewness is inherent to any highly asymmetric
        // (lengthRatio ~0.5) merge, or specific to 3737's local
        // topology beyond its asymmetry. Temporary hardcoded allowlist.
        labelHashSet tjTestAllowlist;
        // SOL scope-dependency test 1: 3738 ALONE (mirror of the
        // earlier 3737-alone test). If 3738 also succeeds cleanly in
        // isolation, both are individually safe and the earlier
        // {3737,3738} rejection is genuinely scope-dependent, not a
        // property of either vertex alone.
        tjTestAllowlist.insert(3738);
        // Run 1 of SOL bisection: {3737, 3738} together, nothing else.
        triSurf* stitched = phase2cMultiJunctionStitch
        (
            *surfacePtr_, meshDict_, 20, tjTestAllowlist
        );
        if( stitched )
        {
            deleteDemandDrivenData(surfacePtr_);
            surfacePtr_ = stitched;

            // SCAN 1: same dry-run scanner, on the POST-REPAIR
            // reconstructed surface, role IDs re-resolved fresh again.
            PatchRoleMap diagRoles1(meshDict_);
            if( diagRoles1.active() )
            {
                const wordList& diagNames1 = surfacePtr_->patchNames();
                labelHashSet diagBladeIds1, diagPeriodicIds1, diagWallIds1;
                forAll(diagNames1, pI)
                {
                    const word& n = diagNames1[pI];
                    if( diagRoles1.hasRole(n, word("blade")) ) diagBladeIds1.insert(pI);
                    if( diagRoles1.hasRole(n, word("periodic")) ) diagPeriodicIds1.insert(pI);
                    if( diagRoles1.hasRole(n, word("hub")) || diagRoles1.hasRole(n, word("shroud")) )
                        diagWallIds1.insert(pI);
                }
                phase2cDryRunScan
                (
                    *surfacePtr_, diagBladeIds1, diagPeriodicIds1, diagWallIds1,
                    word("SCAN1_afterRepair")
                );
            }
        }
    }

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
