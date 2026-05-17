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
#include "checkBoundaryFacesSharingTwoEdges.H"
#include "triSurfaceMetaData.H"
#include "polyMeshGenGeometryModification.H"
#include "surfaceMeshGeometryModification.H"

//#define DEBUG

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private member functions  * * * * * * * * * * * * //

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

    //- snap corner and edge vertices onto feature curves
    //- prevents mesh protrusion at patch intersection edges
    mapper.mapCornersAndEdges();

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
    meshSurfaceEdgeExtractorNonTopo(mesh_, *octreePtr_);
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
    Info << "DEBUG: terminateLayersAtConcaveEdges called" << endl;
    bl.addLayerForAllPatches();
    // Capture junction points for handoff to refineBoundaryLayers
    blblJunctionPoints_ = bl.junctionEdgePoints();
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

        meshOptimizer mOpt(mesh_);
        mOpt.lockPoints(blPoints_);
        mOpt.untangleBoundaryLayer();
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

    if( true )
    {
        meshSurfaceEngine mse(mesh_);
        meshSurfaceOptimizer surfOpt(mse, *octreePtr_);

        if( enforceConstraints )
            surfOpt.enforceConstraints();

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

        if( controller_.runCurrentStep("surfaceProjection") )
        {
            mapMeshToSurface();
        }

        if( controller_.runCurrentStep("patchAssignment") )
        {
            extractPatches();
        }

        if( controller_.runCurrentStep("edgeExtraction") )
        {
            mapEdgesAndCorners();

            optimiseMeshSurface();

            // Step 1: snap corner points first
            {
                meshSurfaceEngine mse(mesh_);
                meshSurfacePartitioner mPart(mse);
                meshSurfaceMapper mapper(mse, *octreePtr_);
                const labelHashSet& corners = mPart.corners();
                labelLongList cornerPts;
                forAllConstIter(labelHashSet, corners, it)
                    cornerPts.append(it.key());
                Info << "Ordered snap: snapping "
                     << cornerPts.size()
                     << " corner points" << endl;
                mapper.mapCorners(cornerPts);
            }

            // Step 1: snap corner points first (damped relaxation)
            {
                // Read relaxation factor from meshDict, default 0.25
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
    if( true )
    {
        checkMeshDict cmd(meshDict_);
    }

    fileName surfaceFile = meshDict_.lookup("surfaceFile");
    if( Pstream::parRun() )
        surfaceFile = ".."/surfaceFile;

    surfacePtr_ = new triSurf(db_.path()/surfaceFile);

    if( true )
    {
        //- save meta data with the mesh (surface mesh + its topology info)
        triSurfaceMetaData sMetaData(*surfacePtr_);
        const dictionary& surfMetaDict = sMetaData.metaData();

        mesh_.metaData().add("surfaceFile", surfaceFile, true);
        mesh_.metaData().add("surfaceMeta", surfMetaDict, true);
    }

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
