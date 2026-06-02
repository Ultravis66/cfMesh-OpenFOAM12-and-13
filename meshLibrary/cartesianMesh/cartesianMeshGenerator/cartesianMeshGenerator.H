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

Class
    cartesianMeshGenerator

Description
    Creates cartesian mesh from the octree

SourceFiles
    cartesianMeshGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef cartesianMeshGenerator_H
#define cartesianMeshGenerator_H

#include "polyMeshGen.H"
#include "Map.H"
#include "HashSet.H"
#include "IOdictionary.H"
#include "workflowControls.H"
//#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class triSurf;
class boundaryLayers;
class meshOctree;
class Time;

/*---------------------------------------------------------------------------*\
                           Class cartesianMeshGenerator Declaration
\*---------------------------------------------------------------------------*/

class cartesianMeshGenerator
{
    // Private data
        //- reference to Time
        const Time& db_;

        //- pointer to the surface
        const triSurf* surfacePtr_;

        //- pointer to the modified surface
        const triSurf* modSurfacePtr_;

        //- IOdictionary containing information about cell sizes, etc..
        IOdictionary meshDict_;


        //- pointer to the octree
        meshOctree* octreePtr_;

        //- mesh
        polyMeshGen mesh_;

        //- workflow controller
        workflowControls controller_;

        //- BL/BL junction points passed from boundaryLayers to refineBoundaryLayers
        labelHashSet blblJunctionPoints_;
        labelHashSet blblAcuteCornerPoints_;
        labelList vtFaceRing_;

        //- BL/no-BL transition edge points (boundary-point indices)
        //- passed from boundaryLayers to post-BL mapper instances
        //- to exclude from generic nearest-surface projection
        labelHashSet blNoBlEdgePoints_;

        //- BL-side patch index per interface point
        //- used for patch-constrained projection instead of global search
        Map<label> blNoBlPointPatch_;

        //- BL/neutral edge points (blade/periodic junctions)
        labelHashSet blNeutralEdgePoints_;

        //- BL-side patch per neutral edge point
        Map<label> blNeutralPointPatch_;

        //- mesh points inside the boundary layer (from refineBoundaryLayers)
        labelLongList blPoints_;
        label nPointsBeforeBL_;
        scalarField blLayerScale_;

        //- true if optimiseFinalMesh() rejected and rolled back final untangle
        bool finalUntangleRejected_;

    // Private member functions
        //- create cartesian mesh
        void createCartesianMesh();

        //- prepare mesh surface
        void surfacePreparation();

        //- map mesh to the surface and untangle surface
        void mapMeshToSurface();

        //- capture edges and corners
        void extractPatches();

        //- map points at feature edges and corners
        void mapEdgesAndCorners();

        //- optimise surface mesh
        void optimiseMeshSurface();

        //- add boundary layers
        //- Detect tight proximity/gap regions between BL patches and
        //- suppress BL locally through boundaryLayers::setGapPoints().
        void detectGapPoints(boundaryLayers& bl);
        void detectTripleJunctions(boundaryLayers& bl);

        void generateBoundaryLayers();

        //- refine boundary layers
        void refBoundaryLayers();

        //- mesh optimisation
        void optimiseFinalMesh();

        //- re-project points back on the surface mesh after back-scaling
        void projectSurfaceAfterBackScaling();

        //- replace boundaries
        void replaceBoundaries();

        //- renumber the mesh
        void renumberMesh();

        //- generate mesh
        void generateMesh();

        //- Disallow default bitwise copy construct
        cartesianMeshGenerator(const cartesianMeshGenerator&);

        //- Disallow default bitwise assignment
        void operator=(const cartesianMeshGenerator&);

public:

    // Constructors

        //- Construct from time
        cartesianMeshGenerator(const Time&);

    // Destructor

        ~cartesianMeshGenerator();


    // Member Functions

        //- write the mesh
        void writeMesh() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
