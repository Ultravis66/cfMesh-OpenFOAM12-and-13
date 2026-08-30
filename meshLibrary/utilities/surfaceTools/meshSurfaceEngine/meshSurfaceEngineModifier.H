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
    meshSurfaceEngineModifier

Description
    Modifier for the meshSurfaceEngine

SourceFiles
    meshSurfaceEngine.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceEngineModifier_H
#define meshSurfaceEngineModifier_H

#include "meshSurfaceEngine.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class meshSurfaceEngineModifier Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceEngineModifier
{
    // Private data
        //- reference to the meshSurfaceEngine
        meshSurfaceEngine& surfaceEngine_;

        //- Disallow default bitwise copy construct
        meshSurfaceEngineModifier(const meshSurfaceEngineModifier&);

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceEngineModifier&);

        //- Raw signed volume of one cell using live mesh-point
        //- coordinates, optionally substituting one point by a candidate.
        //- Uses the same face-geometry/orientation convention as
        //- polyMeshGenChecks::checkCellVolumes().
        scalar rawCellVolumeWithPoint
        (
            const label cellI,
            const label movedGlobalPtI,
            const point& candidate,
            const bool substituteCandidate
        ) const;

public:
    
    // Constructors

        //- Construct from reference to meshSurfaceEngine
        meshSurfaceEngineModifier(meshSurfaceEngine& surfaceEngine);

        //- Construct from const reference to meshSurfaceEngine
        meshSurfaceEngineModifier(const meshSurfaceEngine& surfaceEngine);

    // Destructor

        ~meshSurfaceEngineModifier();

    // Member Functions
        
        //- Check one proposed boundary-point move against the live
        //- volume mesh. Returns false if any currently-positive incident
        //- cell would become non-positive. Existing bad cells may be
        //- repaired; no relative-volume quality floor is imposed yet.
        bool candidatePreservesPositiveCellVolumes
        (
            const label bpI,
            const point& candidate,
            bool& touchesExistingBad,
            scalar& minPositiveRatio
        ) const;

        //- relocate the selected boundary vertex
        void moveBoundaryVertexNoUpdate(const label bpI, const point& newP);
        
        //- relocate the selected boundary vertex and update geometry data
        void moveBoundaryVertex(const label bpI, const point& newP);
        
        //- update normals of boundary vertices at processor boundaries
        inline void updateVertexNormals()
        {
            if( surfaceEngine_.pointNormalsPtr_ )
                surfaceEngine_.updatePointNormalsAtProcBoundaries();
        }
        
        //- updates faceCentres, face normal, and point normals
        //- of the boundary faces attached to the selected boundary points
        void updateGeometry(const labelLongList&);
        void updateGeometry();
        
        //- makes sure that all surface vertices at parallel boundaries
        //- have the same coordinates
        void syncVerticesAtParallelBoundaries();
        void syncVerticesAtParallelBoundaries(const labelLongList&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
