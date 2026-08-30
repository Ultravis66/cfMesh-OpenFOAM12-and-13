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

#include "correctEdgesBetweenPatches.H"
#include "demandDrivenData.H"
#include "meshSurfaceEngine.H"
#include "decomposeCells.H"
#include "polyMeshGenChecks.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const meshSurfaceEngine& correctEdgesBetweenPatches::meshSurface() const
{
    if( !msePtr_ )
        msePtr_ = new meshSurfaceEngine(mesh_);

    return *msePtr_;
}

//- delete mesh surface
void correctEdgesBetweenPatches::clearMeshSurface()
{
    deleteDemandDrivenData(msePtr_);
}

void correctEdgesBetweenPatches::replaceBoundary()
{
    clearMeshSurface();

    polyMeshGenModifier(mesh_).replaceBoundary
    (
        patchNames_,
        newBoundaryFaces_,
        newBoundaryOwners_,
        newBoundaryPatches_
    );
}

void correctEdgesBetweenPatches::decomposeCorrectedCells()
{
    if( decompose_ )
    {
        clearMeshSurface();

        decomposeCells dc(mesh_);
        dc.decomposeMesh(decomposeCell_);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
correctEdgesBetweenPatches::correctEdgesBetweenPatches(polyMeshGen& mesh)
:
    mesh_(mesh),
    msePtr_(NULL),
    patchNames_(mesh.boundaries().size()),
    patchTypes_(mesh.boundaries().size()),
    newBoundaryFaces_(),
    newBoundaryOwners_(),
    newBoundaryPatches_(),
    decomposeCell_(mesh_.cells().size(), false),
    decompose_(false)
{
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    forAll(boundaries, patchI)
    {
        patchNames_[patchI] = boundaries[patchI].patchName();
        patchTypes_[patchI] = boundaries[patchI].patchType();
    }

    // Diagnostic only: isolate where the original cfMesh
    // correctEdgesBetweenPatches topology pipeline first creates
    // non-positive cells.
    auto topologyLineage =
    [&](const word& stageName)
    {
        labelHashSet negVolCells;

        polyMeshGenChecks::checkCellVolumes
        (
            mesh_,
            false,
            &negVolCells
        );

        // Diagnostic only:
        // reproduce OpenFOAM's repeated-neighbour semantics at each
        // correctEdgesBetweenPatches topology checkpoint.
        const cellListPMG& lineageCells =
            mesh_.cells();

        const labelList& lineageOwner =
            mesh_.owner();

        const labelList& lineageNeighbour =
            mesh_.neighbour();

        const label lineageNInternal =
            mesh_.nInternalFaces();

        label lineageDuplicatePairs = 0;
        label lineageMasterCells = 0;


        forAll(lineageCells, cellI)
        {
            const cell& c =
                lineageCells[cellI];

            labelHashSet seenNeighbours;
            labelHashSet duplicateNeighbours;

            bool masterHasDuplicate = false;


            forAll(c, cfI)
            {
                const label faceI =
                    c[cfI];

                if( faceI >= lineageNInternal )
                    continue;


                label otherCell = -1;

                if
                (
                    lineageOwner[faceI]
                 == cellI
                )
                {
                    otherCell =
                        lineageNeighbour[faceI];
                }
                else if
                (
                    lineageNeighbour[faceI]
                 == cellI
                )
                {
                    otherCell =
                        lineageOwner[faceI];
                }
                else
                {
                    FatalErrorIn
                    (
                        "correctEdgesBetweenPatches "
                        "topologyLineage"
                    )
                        << "Internal face "
                        << faceI
                        << " listed in cell "
                        << cellI
                        << " but owner/neighbour are "
                        << lineageOwner[faceI]
                        << " and "
                        << lineageNeighbour[faceI]
                        << abort(FatalError);
                }


                // Match OpenFOAM's lower-labelled master-cell
                // convention so each unordered pair is counted once.
                if( cellI >= otherCell )
                    continue;


                if
                (
                    seenNeighbours.found
                    (
                        otherCell
                    )
                )
                {
                    if
                    (
                        !duplicateNeighbours.found
                        (
                            otherCell
                        )
                    )
                    {
                        duplicateNeighbours.insert
                        (
                            otherCell
                        );

                        ++lineageDuplicatePairs;
                        masterHasDuplicate = true;
                    }
                }
                else
                {
                    seenNeighbours.insert
                    (
                        otherCell
                    );
                }
            }


            if( masterHasDuplicate )
                ++lineageMasterCells;
        }


        Info
            << "[CORRECT_EDGES_TOPOLOGY_LINEAGE]"
            << " stage=" << stageName
            << " negVol=" << negVolCells.size()
            << " duplicatePairs="
            << lineageDuplicatePairs
            << " duplicateMasterCells="
            << lineageMasterCells
            << endl;
    };

    topologyLineage("entry");

    //decomposeProblematicFaces();

    // EXPERIMENT:
    // Preserve the original internal polygon interfaces instead of
    // fan-decomposing them in decomposeConcaveFaces().
    //
    // Diagnostic lineage has shown that this historical operation
    // creates 308 duplicate-neighbour interfaces from an entry state
    // containing zero duplicate-neighbour pairs.
    Info
        << "[EXPERIMENT_SKIP_CONCAVE_FACE_DECOMP]"
        << " preserving original internal polygon interfaces"
        << endl;

    //decomposeConcaveFaces();

    topologyLineage("afterDecomposeConcaveFaces");

    patchCorrection();

    topologyLineage("afterPatchCorrection");

    decomposeCorrectedCells();

    topologyLineage("afterDecomposeCorrectedCells");
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

correctEdgesBetweenPatches::~correctEdgesBetweenPatches()
{
    deleteDemandDrivenData(msePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
