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

#include "demandDrivenData.H"
#include "meshSurfaceEdgeExtractorNonTopo.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceMapper.H"

#include "correctEdgesBetweenPatches.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceEdgeExtractorNonTopo::decomposeBoundaryFaces()
{
    correctEdgesBetweenPatches featureEdges(mesh_);
    // GeomFix: store new centroid points for targeted projection
    patchCorrectionPoints_ = featureEdges.patchCorrectionPoints();
    Info << "[GeomFix] patchCorrection created "
         << patchCorrectionPoints_.size()
         << " new centroid points" << endl;
}

void meshSurfaceEdgeExtractorNonTopo::remapBoundaryPoints()
{
    meshSurfaceEngine mse(mesh_);
    meshSurfaceMapper mapper(mse, meshOctree_);

    // Pass BL/no-BL protected points if provided
    if( !protectedPoints_.empty() )
    {
        mapper.setProtectedPoints(protectedPoints_);
        mapper.setProtectedPointPatches(protectedPointPatches_);
        Info << "remapBoundaryPoints: protecting "
             << protectedPoints_.size()
             << " BL/no-BL interface points" << endl;
    }

    mapper.mapVerticesOntoSurfacePatches();

    // GeomFix: project new centroid points onto STL
    if( patchCorrectionPoints_.size() > 0 )
    {
        const labelList& bPoints = mse.boundaryPoints();
        labelList meshToBnd(mesh_.points().size(), -1);
        forAll(bPoints, bpI)
            meshToBnd[bPoints[bpI]] = bpI;
        labelLongList bpToProject;
        forAll(patchCorrectionPoints_, i)
        {
            const label ptI = patchCorrectionPoints_[i];
            if( ptI < meshToBnd.size() && meshToBnd[ptI] >= 0 )
                bpToProject.append(meshToBnd[ptI]);
        }
        Info << "[GeomFix] projecting "
             << bpToProject.size()
             << " centroid points onto STL" << endl;
        if( bpToProject.size() > 0 )
            mapper.mapVerticesOntoSurface(bpToProject);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
