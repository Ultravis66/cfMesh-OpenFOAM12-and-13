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
    const bool verbose = false
)
{
    const pointField& pts = surf.points();
    const wordList pNames = surf.patchNames();
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

                        // Collect patch names — always, not just verbose
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

                        if( ratio < 0.20 )
                            classification = "WELD_SAFE_RATIO";
                        else if( samePatch && ratio < 0.80 )
                            classification = "WELD_SAME_PATCH";
                        else if( !pNamesA.size() || !pNamesB.size() )
                            classification = "FLAG_INCOMPLETE";
                        else if( ratio > 0.80 )
                            classification = "FLAG_HIGH_RATIO";
                        else
                            classification = "FLAG_CROSS_PATCH_CANDIDATE";

                        if( verbose )
                        {
                            Info << "  Near-coincident pair: "
                                 << a << " " << b
                                 << " d=" << d*1000.0 << "mm"
                                 << " ratio=" << ratio
                                 << " class=" << classification
                                 << " at " << pts[pointI]
                                 << " pA=" << pNamesA
                                 << " pB=" << pNamesB
                                 << endl;
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

void cartesianMeshGenerator::generateBoundaryLayers()
{
    //- add boundary layers
    boundaryLayers bl(mesh_, meshDict_);
    bl.terminateLayersAtConcaveEdges();
    bl.addLayerForAllPatches();
    // Capture junction points for handoff to refineBoundaryLayers
    blblJunctionPoints_ = bl.junctionEdgePoints();
    blblAcuteCornerPoints_ = bl.blblAcuteCornerPoints();
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

        refLayers.refineLayers();

             refLayers.pointsInBndLayer(blPoints_);


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
        // Post-refinement BL optimisation — optional, off by default.
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
    //- untangle the surface if needed
    bool enforceConstraints(false);
    if( meshDict_.found("enforceGeometryConstraints") )
    {
        enforceConstraints =
            readBool(meshDict_.lookup("enforceGeometryConstraints"));
    }

    {
        meshSurfaceEngine mse(mesh_);
        meshSurfaceOptimizer surfOpt(mse, *octreePtr_);


        if( enforceConstraints )
            surfOpt.enforceConstraints();


        //- lock acute BL+BL+neutral corners: prevent optimizer
        //- from moving these points across patch boundaries
        //- controlled by meshDict: lockAcuteCornerPoints true/false
        bool lockAcuteCorners = false;
        if( meshDict_.isDict("boundaryLayers") )
        {
            const dictionary& bndL =
                meshDict_.subDict("boundaryLayers");
            if( bndL.found("lockAcuteCornerPoints") )
                lockAcuteCorners =
                    bool(Switch(bndL.lookup("lockAcuteCornerPoints")));
        }
        if( lockAcuteCorners && blblAcuteCornerPoints_.size() )
        {
            Info << "Locking " << blblAcuteCornerPoints_.size()
                 << " acute BL+BL+neutral corner points in optimizer" << endl;
            surfOpt.setAcuteCornerPoints(blblAcuteCornerPoints_);
        }
        surfOpt.optimizeSurface();
    }

    deleteDemandDrivenData(octreePtr_);

    //- final optimisation
    meshOptimizer optimizer(mesh_);
    if( enforceConstraints )
        optimizer.enforceConstraints();

    optimizer.optimizeMeshFV();
    optimizer.optimizeLowQualityFaces();
    optimizer.optimizeBoundaryLayer(modSurfacePtr_==NULL);
    optimizer.untangleMeshFV();

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
            // edgeExtractor uses only mesh topology + octree — no
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
            // Uses meshDict nLayersForPatch only — no mesh modification.
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
            // Single pass — full BL/no-BL and BL/neutral protection.
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

        if( controller_.runCurrentStep("boundaryLayerGeneration") )
        {
            generateBoundaryLayers();
        }

        if( controller_.runCurrentStep("meshOptimisation") )
        {
            optimiseFinalMesh();

            projectSurfaceAfterBackScaling();
        }

        if( controller_.runCurrentStep("boundaryLayerRefinement") )
        {
            refBoundaryLayers();
        }

        // Post-BL geometry snap - optional, controlled by meshDict
        {
            bool doSnap(false);
            if( meshDict_.isDict("boundaryLayers") )
            {
                const dictionary& bndL =
                    meshDict_.subDict("boundaryLayers");
                if( bndL.found("postBLSnap") )
                    doSnap = Switch(bndL.lookup("postBLSnap"));
            }
            if( doSnap )
            {
                Info << "Post-BL snap: excluded "
                     << blPoints_.size()
                     << " BL interior points" << endl;
                // Rebuild local octree - global octreePtr_ deleted by optimiseFinalMesh
                meshOctree* snapOctreePtr = new meshOctree(*surfacePtr_);
                meshOctreeCreator
                (
                    *snapOctreePtr,
                    meshDict_
                ).createOctreeWithRefinedBoundary(20, 30);
                meshSurfaceEngine mse(mesh_);
                meshSurfaceMapper mapper(mse, *snapOctreePtr);
                // Exclude BL/no-BL interface points from generic snapping
                // These points are constrained to their feature curve
                // and must not be moved by nearest-surface projection
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
                const labelList& bPoints = mse.boundaryPoints();
                boolList isBLPoint(mesh_.points().size(), false);
                forAll(blPoints_, i)
                    isBLPoint[blPoints_[i]] = true;
                labelLongList outerBndPoints;
                forAll(bPoints, bpI)
                    if( !isBLPoint[bPoints[bpI]] )
                        outerBndPoints.append(bpI);
                Info << "Post-BL snap: projecting "
                     << outerBndPoints.size()
                     << " outer boundary points onto STL" << endl;
                if( outerBndPoints.size() > 0 )
                    mapper.mapVerticesOntoSurface(outerBndPoints);
                deleteDemandDrivenData(snapOctreePtr);
            }
        }

        renumberMesh();

        replaceBoundaries();

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
    controller_(mesh_)
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
                // Non-mutating scan — no geometry modification
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
                    writeDiag
                );
                deleteDemandDrivenData(scanOctree);
                Info << "  Near-coincident pairs found: " << nFound << endl;
                Info << "  No geometry modified (reportOnly true)" << endl;
            }
            else if( weldNearPoints )
            {
                meshOctree* preprocOctree = new meshOctree(*surfacePtr_);
                meshOctreeCreator
                (
                    *preprocOctree,
                    meshDict_
                ).createOctreeBoxes();
                triSurfaceCleanupDuplicates cleaner(*preprocOctree, weldTol);
                cleaner.mergeIdentities();
                deleteDemandDrivenData(preprocOctree);
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
