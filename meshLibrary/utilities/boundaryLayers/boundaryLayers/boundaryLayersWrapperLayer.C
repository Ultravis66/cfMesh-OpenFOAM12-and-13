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
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayers::addWrapperLayer()
{
    createOTopologyLayers();

    if( treatedPatch_.size() && treatedPatch_[0] ) return;

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    boolList treatPatches(mesh_.boundaries().size(), true);

    pointFieldPMG& points = mesh_.points();

    // Fix: was populating local newLabelForVertex instead of member
    // newLabelForVertex_ -- createNewFacesAndCells uses the member
    const label oldNPoints = points.size();
    nPoints_ = oldNPoints;
    newLabelForVertex_.setSize(oldNPoints);
    newLabelForVertex_ = -1;
    points.setSize(oldNPoints + bPoints.size());

    forAll(bPoints, bpI)
    {
        const label oldPointI = bPoints[bpI];
        points[nPoints_] = points[oldPointI];
        newLabelForVertex_[oldPointI] = nPoints_;
        ++nPoints_;
    }

    createNewFacesAndCells(treatPatches);

    forAll(treatPatches, patchI)
        if( treatPatches[patchI] )
            treatedPatch_[patchI] = true;

    //- delete surface engine
    clearOut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
