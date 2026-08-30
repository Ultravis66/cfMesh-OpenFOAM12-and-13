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

#include "decomposeCells.H"
#include "demandDrivenData.H"
#include "polyMeshGenAddressing.H"
#include "meshSurfaceEngine.H"
#include "decomposeFaces.H"
#include "labelLongList.H"
#include "polyMeshGenChecks.H"

//#define DEBUGDecompose

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void decomposeCells::decomposeMesh(const boolList& decomposeCell)
{
    if( decomposeCell.size() != mesh_.cells().size() )
    {
        FatalErrorIn
        (
            "void decomposeCells::decomposeMesh(const boolList&)"
        )   << "Incorrect decomposeCell size " << decomposeCell.size()
            << ", mesh has " << mesh_.cells().size() << " cells"
            << abort(FatalError);
    }

    // Diagnostic only. Check volumes only at complete mesh states.
    // Do not inspect between parent removal and child addition.
    auto decomposeVolumeLineage =
    [&](const word& stageName)
    {
        labelHashSet negVolCells;

        polyMeshGenChecks::checkCellVolumes
        (
            mesh_,
            false,
            &negVolCells
        );

        Info
            << "[DECOMPOSE_VOLUME_LINEAGE]"
            << " stage=" << stageName
            << " negVol=" << negVolCells.size()
            << endl;
    };

    decomposeVolumeLineage("entry");

    label nRequested = 0;

    forAll(decomposeCell, cellI)
    {
        if( decomposeCell[cellI] )
            ++nRequested;
    }

    if( nRequested == 0 )
    {
        Info
            << "DECOMPTOPPAIR requested=0"
            << " duplicatePairs=0"
            << " selectedCells=0"
            << endl;

        return;
    }


    // Diagnostic only:
    // Count repeated-neighbour topology using the same lower-labelled
    // master-cell convention as OpenFOAM checkUpperTriangular().
    //
    // References are acquired fresh on every call because
    // checkFaceConnections() rebuilds the face/cell topology.
    auto duplicatePairLineage =
    [&]
    (
        const word& stageName
    )
    {
        const cellListPMG& dCells =
            mesh_.cells();

        const labelList& dOwner =
            mesh_.owner();

        const labelList& dNeighbour =
            mesh_.neighbour();

        const label dNInternal =
            mesh_.nInternalFaces();


        label nDuplicatePairs = 0;
        label nMasterCells = 0;


        forAll(dCells, cellI)
        {
            const cell& c = dCells[cellI];

            labelHashSet seenNeighbours;
            labelHashSet duplicateNeighbours;

            bool hasDuplicate = false;


            forAll(c, cfI)
            {
                const label faceI = c[cfI];

                if( faceI >= dNInternal )
                    continue;


                label otherCell = -1;

                if( dOwner[faceI] == cellI )
                {
                    otherCell =
                        dNeighbour[faceI];
                }
                else if
                (
                    dNeighbour[faceI] == cellI
                )
                {
                    otherCell =
                        dOwner[faceI];
                }
                else
                {
                    FatalErrorIn
                    (
                        "decomposeCells duplicatePairLineage"
                    )
                        << "Internal face " << faceI
                        << " listed in cell " << cellI
                        << " but owner/neighbour are "
                        << dOwner[faceI] << " and "
                        << dNeighbour[faceI]
                        << abort(FatalError);
                }


                // Match OpenFOAM: lower-labelled cell owns
                // the unordered pair check.
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

                        ++nDuplicatePairs;
                        hasDuplicate = true;
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


            if( hasDuplicate )
                ++nMasterCells;
        }


        Info
            << "[DECOMPOSE_DUPLICATE_LINEAGE]"
            << " stage=" << stageName
            << " uniquePairs=" << nDuplicatePairs
            << " masterCells=" << nMasterCells
            << endl;
    };


    duplicatePairLineage("entry");


    // -----------------------------------------------------------------
    // Phase 1:
    // Keep the useful historical face-connectivity repair.
    //
    // This can split faces which would otherwise produce invalid local
    // face connectivity, but it does not delete or renumber cells.
    // -----------------------------------------------------------------

    checkFaceConnections(decomposeCell);

    decomposeVolumeLineage("afterCheckFaceConnections");

    duplicatePairLineage("afterCheckFaceConnections");


    if( decomposeCell.size() != mesh_.cells().size() )
    {
        FatalErrorIn
        (
            "void decomposeCells::decomposeMesh(const boolList&)"
        )   << "Cell count changed unexpectedly during "
            << "checkFaceConnections()"
            << abort(FatalError);
    }


    const cellListPMG& cells = mesh_.cells();
    const faceListPMG& faces = mesh_.faces();

    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();

    const label nInternalFaces = mesh_.nInternalFaces();


    // -----------------------------------------------------------------
    // Phase 2:
    // Preflight the historical candidate population geometrically.
    //
    // This does NOT mean all safe candidates will be decomposed.
    // It merely tells the topology selector which endpoint of an
    // offending neighbour pair is eligible for one-apex decomposition.
    // -----------------------------------------------------------------

    boolList apexSafe(decomposeCell.size(), false);
    List<scalar> apexMargin(decomposeCell.size(), -GREAT);

    label nSafeRequested = 0;
    label nUnsafeRequested = 0;

    forAll(decomposeCell, cellI)
    {
        if( !decomposeCell[cellI] )
            continue;

        point apex(point::zero);
        scalar relativeMargin = -GREAT;

        if
        (
            findValidPyramidApex
            (
                cellI,
                apex,
                relativeMargin
            )
        )
        {
            scalar minChildVolume = GREAT;

            if
            (
                exactPyramidChildrenPositive
                (
                    cellI,
                    apex,
                    minChildVolume
                )
            )
            {
                apexSafe[cellI] = true;
                apexMargin[cellI] = relativeMargin;
                ++nSafeRequested;
            }
            else
            {
                apexMargin[cellI] = relativeMargin;
                ++nUnsafeRequested;

                Info
                    << "[DECOMPOSE_EXACT_CHILD_REJECT]"
                    << " cell=" << cellI
                    << " apexMargin=" << relativeMargin
                    << " minChildVolume=" << minChildVolume
                    << endl;
            }
        }
        else
        {
            apexMargin[cellI] = relativeMargin;
            ++nUnsafeRequested;
        }
    }


    // Diagnostic only: provide a topology-independent spatial
    // fingerprint for cells across later renumbering/reconstruction.
    // Use the arithmetic mean of UNIQUE cell vertices.
    auto diagnosticCellVertexCentre =
    [&]
    (
        const label cellI
    ) -> point
    {
        const cell& dc = mesh_.cells()[cellI];
        const faceListPMG& df = mesh_.faces();
        const pointFieldPMG& dp = mesh_.points();

        labelHashSet cellPoints;

        forAll(dc, cfI)
        {
            const face& f = df[dc[cfI]];

            forAll(f, fpI)
                cellPoints.insert(f[fpI]);
        }

        if( cellPoints.empty() )
            return point::zero;

        point c(point::zero);

        forAllConstIter(labelHashSet, cellPoints, it)
            c += dp[it.key()];

        c /= scalar(cellPoints.size());

        return c;
    };


    // -----------------------------------------------------------------
    // Phase 3:
    // Find the ACTUAL topology defect that requires cell decomposition:
    //
    //     two neighbouring cells connected through >1 internal face
    //
    // Process every offending pair once.  Instead of decomposing the
    // complete broad decomposeCell mask, select ONE geometrically safe
    // endpoint of each offending pair.
    //
    // If both endpoints are eligible, choose the endpoint with the
    // larger star-kernel margin.
    //
    // A single selected cell can resolve several offending pairs.
    // -----------------------------------------------------------------

    boolList topologyDecomposeCell(decomposeCell.size(), false);

    label nDuplicatePairs = 0;
    label nInScopePairs = 0;
    label nResolvedPairs = 0;
    label nUnresolvedUnsafePairs = 0;
    label nOutOfScopePairs = 0;

    label nSelectionPrint = 0;
    label nUnresolvedPrint = 0;


    forAll(cells, cellI)
    {
        const cell& c = cells[cellI];

        labelHashSet seenNeighbours;
        labelHashSet duplicateSet;

        DynList<label, 16> duplicateNeighbours;


        forAll(c, cfI)
        {
            const label faceI = c[cfI];

            if( faceI >= nInternalFaces )
                continue;

            label otherCell = -1;

            if( owner[faceI] == cellI )
            {
                otherCell = neighbour[faceI];
            }
            else if( neighbour[faceI] == cellI )
            {
                otherCell = owner[faceI];
            }
            else
            {
                FatalErrorIn
                (
                    "void decomposeCells::decomposeMesh"
                    "(const boolList&)"
                )   << "Internal face " << faceI
                    << " listed in cell " << cellI
                    << " but owner/neighbour are "
                    << owner[faceI] << " and "
                    << neighbour[faceI]
                    << abort(FatalError);
            }


            if( seenNeighbours.found(otherCell) )
            {
                if( !duplicateSet.found(otherCell) )
                {
                    duplicateSet.insert(otherCell);
                    duplicateNeighbours.append(otherCell);
                }
            }
            else
            {
                seenNeighbours.insert(otherCell);
            }
        }


        forAll(duplicateNeighbours, dnI)
        {
            const label otherCell =
                duplicateNeighbours[dnI];

            // Process this unordered pair only once.
            if( cellI > otherCell )
                continue;

            ++nDuplicatePairs;


            const bool thisRequested =
                decomposeCell[cellI];

            const bool otherRequested =
                decomposeCell[otherCell];


            if( !thisRequested && !otherRequested )
            {
                ++nOutOfScopePairs;
                continue;
            }

            ++nInScopePairs;


            const bool thisSafe =
                thisRequested && apexSafe[cellI];

            const bool otherSafe =
                otherRequested && apexSafe[otherCell];


            label chosenCell = -1;


            if( thisSafe && otherSafe )
            {
                // Pick the more robust one-apex decomposition.
                if
                (
                    apexMargin[cellI]
                  >= apexMargin[otherCell]
                )
                {
                    chosenCell = cellI;
                }
                else
                {
                    chosenCell = otherCell;
                }
            }
            else if( thisSafe )
            {
                chosenCell = cellI;
            }
            else if( otherSafe )
            {
                chosenCell = otherCell;
            }


            if( chosenCell >= 0 )
            {
                topologyDecomposeCell[chosenCell] = true;
                ++nResolvedPairs;

                if( nSelectionPrint < 20 )
                {
                    Info
                        << "DECOMPTOPPAIR_SELECT"
                        << " pair=("
                        << cellI << ' ' << otherCell << ')'
                        << " chosen=" << chosenCell
                        << " margin="
                        << apexMargin[chosenCell]
                        << endl;

                    ++nSelectionPrint;
                }
            }
            else
            {
                ++nUnresolvedUnsafePairs;

                // Diagnostic only: print EVERY unresolved unordered pair,
                // not merely the first 20, and attach a spatial fingerprint
                // which can be correlated after later cell renumbering.
                {
                    label nSharedFaces = 0;

                    forAll(c, cfI)
                    {
                        const label faceI = c[cfI];

                        if( faceI >= nInternalFaces )
                            continue;

                        label nbrCell = -1;

                        if( owner[faceI] == cellI )
                            nbrCell = neighbour[faceI];
                        else if( neighbour[faceI] == cellI )
                            nbrCell = owner[faceI];

                        if( nbrCell == otherCell )
                            ++nSharedFaces;
                    }

                    const point centreA =
                        diagnosticCellVertexCentre(cellI);

                    const point centreB =
                        diagnosticCellVertexCentre(otherCell);

                    const point mid =
                        0.5*(centreA + centreB);

                    Info
                        << "[DECOMPTOPPAIR_UNRESOLVED_ALL]"
                        << " pair=("
                        << cellI << ' ' << otherCell << ')'
                        << " nSharedFaces=" << nSharedFaces
                        << " requested=("
                        << thisRequested << ' '
                        << otherRequested << ')'
                        << " safe=("
                        << thisSafe << ' '
                        << otherSafe << ')'
                        << " margin=("
                        << apexMargin[cellI] << ' '
                        << apexMargin[otherCell] << ')'
                        << " centreA=" << centreA
                        << " centreB=" << centreB
                        << " mid=" << mid
                        << endl;
                }

                if( nUnresolvedPrint < 20 )
                {
                    Info
                        << "DECOMPTOPPAIR_UNRESOLVED"
                        << " pair=("
                        << cellI << ' ' << otherCell << ')'
                        << " requested=("
                        << thisRequested << ' '
                        << otherRequested << ')'
                        << " safe=("
                        << thisSafe << ' '
                        << otherSafe << ')'
                        << " margin=("
                        << apexMargin[cellI] << ' '
                        << apexMargin[otherCell] << ')'
                        << endl;

                    ++nUnresolvedPrint;
                }
            }
        }
    }


    label nSelectedCells = 0;

    forAll(topologyDecomposeCell, cellI)
    {
        if( topologyDecomposeCell[cellI] )
            ++nSelectedCells;
    }


    label nSelectedPrint = 0;

    forAll(topologyDecomposeCell, cellI)
    {
        if
        (
            topologyDecomposeCell[cellI]
         && nSelectedPrint < 50
        )
        {
            Info
                << "[DECOMPOSE_SELECTED_PARENT]"
                << " cell=" << cellI
                << " apexMargin=" << apexMargin[cellI]
                << endl;

            ++nSelectedPrint;
        }
    }

    Info
        << "DECOMPTOPPAIR"
        << " requested=" << nRequested
        << " safeRequested=" << nSafeRequested
        << " unsafeRequested=" << nUnsafeRequested
        << " duplicatePairs=" << nDuplicatePairs
        << " inScopePairs=" << nInScopePairs
        << " resolvedPairs=" << nResolvedPairs
        << " unresolvedUnsafePairs="
        << nUnresolvedUnsafePairs
        << " outOfScopePairs=" << nOutOfScopePairs
        << " selectedCells=" << nSelectedCells
        << endl;


    // Nothing requires safe volume decomposition.
    if( nSelectedCells == 0 )
    {
        Info
            << "DECOMPTOPPAIR: retaining all corrected "
            << "parent polyhedra"
            << endl;

        return;
    }


    // -----------------------------------------------------------------
    // Phase 4:
    // Pyramid-decompose ONLY the minimal topology-driven cell subset.
    //
    // findTopVertex() independently re-runs the robust apex validation,
    // so even this selected mask remains fail-closed.
    // -----------------------------------------------------------------

    createPointsAndCellFaces(topologyDecomposeCell);

    storeBoundaryFaces(topologyDecomposeCell);

    removeDecomposedCells(topologyDecomposeCell);

    // Diagnostic only: addCells() appends facesOfNewCells_ in order.
    const label firstNewCell = mesh_.cells().size();
    const label nChildRecords = facesOfNewCells_.size();

    Info
        << "[DECOMPOSE_CHILD_RANGE]"
        << " firstNewCell=" << firstNewCell
        << " nChildRecords=" << nChildRecords
        << endl;

    addNewCells();

    // Diagnostic only: identify the exact final negative children and
    // map their final cell labels back to facesOfNewCells_ record IDs.
    labelHashSet postAddNegVolCells;

    polyMeshGenChecks::checkCellVolumes
    (
        mesh_,
        false,
        &postAddNegVolCells
    );

    for
    (
        label cellI=0;
        cellI<mesh_.cells().size();
        ++cellI
    )
    {
        if( !postAddNegVolCells.found(cellI) )
            continue;

        const label record =
            cellI - firstNewCell;

        const bool isNewChild =
        (
            record >= 0
         && record < nChildRecords
        );

        Info
            << "[DECOMPOSE_BAD_CHILD]"
            << " cell=" << cellI
            << " isNewChild=" << isNewChild
            << " record=" << record
            << endl;
    }

    decomposeVolumeLineage("afterAddNewCells");


    # ifdef DEBUGDecompose
    mesh_.addressingData().checkMesh();
    # endif
}


void decomposeCells::checkFaceConnections(const boolList& decomposeCell)
{
    const faceListPMG& faces = mesh_.faces();
    const cellListPMG& cells = mesh_.cells();

    boolList decomposeFace(faces.size(), false);
    forAll(cells, cellI)
    {
        if( decomposeCell[cellI] )
        {
            DynList<label, 32> vrt;
            DynList<edge, 64> edges;
            DynList<DynList<label, 8> > faceEdges;
            DynList<DynList<label, 2>, 64> edgeFaces;

            findAddressingForCell(cellI, vrt, edges, faceEdges, edgeFaces);

            forAll(faceEdges, fI)
            {
                const DynList<label, 8>& fEdges = faceEdges[fI];

                labelHashSet neiFaces;
                forAll(fEdges, feI)
                {
                    label neiFace = edgeFaces[fEdges[feI]][0];
                    if( neiFace == fI )
                        neiFace = edgeFaces[fEdges[feI]][1];

                    if( neiFaces.found(neiFace) )
                    {
                        decomposeFace[cells[cellI][fI]] = true;
                    }
                    else
                    {
                        neiFaces.insert(neiFace);
                    }
                }
            }
        }
    }

    if( Pstream::parRun() )
    {
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh_.procBoundaries();

        //- send information to the neighbour processor
        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            boolList decFace(procBoundaries[patchI].patchSize(), false);
            const label size = decFace.size();

            for(label i=0;i<size;++i)
            {
                if( decomposeFace[start+i] )
                    decFace[i] = true;
            }

            OPstream toOtherProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo(),
                decFace.byteSize()
            );

            toOtherProc << decFace;
        }

        //- receive information from the neighbour processor
        forAll(procBoundaries, patchI)
        {
            boolList decFace;

            IPstream fromOtherProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo()
            );

            fromOtherProc >> decFace;

            const label start = procBoundaries[patchI].patchStart();
            forAll(decFace, i)
            {
                if( decFace[i] )
                    decomposeFace[start+i] = true;
            }
        }
    }

    // Diagnostic only:
    // classify exactly which original faces are about to be fan-decomposed.
    {
        const label nInternal =
            mesh_.nInternalFaces();

        label nMarkedInternal = 0;
        label nMarkedNonInternal = 0;

        label nInternalTri = 0;
        label nInternalQuad = 0;
        label nInternalPent = 0;
        label nInternalHex = 0;
        label nInternalHept = 0;
        label nInternalGt7 = 0;

        label nProducedInternalTriangles = 0;


        forAll(decomposeFace, faceI)
        {
            if( !decomposeFace[faceI] )
                continue;


            if( faceI < nInternal )
            {
                ++nMarkedInternal;

                const label nVerts =
                    faces[faceI].size();

                nProducedInternalTriangles +=
                    nVerts;


                if( nVerts == 3 )
                    ++nInternalTri;
                else if( nVerts == 4 )
                    ++nInternalQuad;
                else if( nVerts == 5 )
                    ++nInternalPent;
                else if( nVerts == 6 )
                    ++nInternalHex;
                else if( nVerts == 7 )
                    ++nInternalHept;
                else
                    ++nInternalGt7;
            }
            else
            {
                ++nMarkedNonInternal;
            }
        }


        Info
            << "[CHECKFACECONNECTIONS_MARKED]"
            << " internalFaces="
            << nMarkedInternal
            << " nonInternalFaces="
            << nMarkedNonInternal
            << " tri="
            << nInternalTri
            << " quad="
            << nInternalQuad
            << " pent="
            << nInternalPent
            << " hex="
            << nInternalHex
            << " hept="
            << nInternalHept
            << " gt7="
            << nInternalGt7
            << " producedInternalTriangles="
            << nProducedInternalTriangles
            << endl;
    }


    // EXPERIMENT:
    // Keep the connectivity diagnostic above, but preserve the original
    // polygon interfaces rather than fan-decomposing each marked face.
    //
    // The historical fan operation converts one owner/neighbour interface
    // into N triangular interfaces between the same two cells, thereby
    // creating duplicate-neighbour topology which is subsequently repaired
    // by whole-cell pyramid decomposition.
    Info
        << "[EXPERIMENT_SKIP_CHECKFACECONNECTIONS_FAN]"
        << " preserving "
        << "marked polygon interfaces"
        << endl;

    //decomposeFaces(mesh_).decomposeMeshFaces(decomposeFace);
}

void decomposeCells::createPointsAndCellFaces(const boolList& decomposeCell)
{
    facesOfNewCells_.clear();

    forAll(decomposeCell, cI)
        if( decomposeCell[cI] )
        {
            decomposeCellIntoPyramids(cI);
        }
}

void decomposeCells::storeBoundaryFaces(const boolList& /*decomposeCell*/)
{
    meshSurfaceEngine mse(mesh_);
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& facePatch = mse.boundaryFacePatches();

    forAll(bFaces, bfI)
    {
        newBoundaryFaces_.appendList(bFaces[bfI]);
        newBoundaryPatches_.append(facePatch[bfI]);
    }
}

void decomposeCells::removeDecomposedCells(const boolList& decomposeCell)
{
    # ifdef DEBUGDecompose
    Info << "Number of cells before removal " << mesh_.cells().size() << endl;
    # endif

    polyMeshGenModifier meshModifier(mesh_);
    meshModifier.removeCells(decomposeCell, false);

    # ifdef DEBUGDecompose
    Info << "Number of cells after removal " << mesh_.cells().size() << endl;
    # endif
}

void decomposeCells::addNewCells()
{
    Info << "Adding new cells " << endl;
    polyMeshGenModifier(mesh_).addCells(facesOfNewCells_);
    facesOfNewCells_.clear();
    Info << "Reordering bnd faces" << endl;
    polyMeshGenModifier(mesh_).reorderBoundaryFaces();

    Info << "Finding bnd faces" << endl;
    const faceListPMG& faces = mesh_.faces();
    const labelList& owner = mesh_.owner();
    const VRWGraph& pointFaces = mesh_.addressingData().pointFaces();

    labelLongList newBoundaryOwners;

    forAll(newBoundaryFaces_, faceI)
    {
        face bf(newBoundaryFaces_.sizeOfRow(faceI));
        forAllRow(newBoundaryFaces_, faceI, pI)
            bf[pI] = newBoundaryFaces_(faceI, pI);

        # ifdef DEBUGDecompose
        Info << "Finding cell for boundary face " << bf << endl;
        bool found(false);
        forAllRow(pointFaces, bf[0], pfI)
            if( bf == faces[pointFaces(bf[0], pfI)] )
                found = true;
        if( !found )
            FatalErrorIn
            (
                "void decomposeCells::addNewCells()"
            ) << "Face " << bf << " does not exist in the mesh"
                << abort(FatalError);
        #endif

        forAllRow(pointFaces, bf[0], pfI)
        {
            const label fLabel = pointFaces(bf[0], pfI);
            if( (mesh_.faceIsInPatch(fLabel) != -1) && (bf == faces[fLabel]) )
            {
                # ifdef DEBUGDecompose
                Info << "Boundary face " << bf << " is in cell "
                << owner[fLabel]] << endl;
                # endif

                newBoundaryOwners.append(owner[fLabel]);
            }
        }
    }

    polyMeshGenModifier(mesh_).replaceBoundary
    (
        patchNames_,
        newBoundaryFaces_,
        newBoundaryOwners,
        newBoundaryPatches_
    );

    polyMeshGenModifier(mesh_).removeUnusedVertices();
    polyMeshGenModifier(mesh_).clearAll();

    PtrList<boundaryPatch>& boundaries =
        polyMeshGenModifier(mesh_).boundariesAccess();
    forAll(boundaries, patchI)
        boundaries[patchI].patchType() = patchTypes_[patchI];
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// ************************************************************************* //
