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

#include "refineBoundaryLayers.H"
#include "polyMeshGenAddressing.H"
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "demandDrivenData.H"
#include "pyramidPointFaceRef.H"
#include <map>
#include <set>
#include <utility>
#include <cstring>

//#define DEBUGLayer

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool refineBoundaryLayers::generateNewCellsPrism
(
    const label cellI,
    DynList<DynList<DynList<label, 8>, 10>, 64>& cellsFromCell
) const
{
    cellsFromCell.clear();

    // --------------------------------------------------------------
    // CFMitch V3.5 -- fail-closed prism structural preflight.
    //
    // A topology-retreat candidate can convert a multidirectional
    // BL-intersection parent into refType 1 while leaving face/split-edge
    // metadata that does not satisfy the classical prism contract.
    //
    // Historically generateNewCellsPrism() indexed this metadata without
    // validation.  An unresolved opposite face or split edge therefore
    // became an index of -1 and caused SIGSEGV.
    //
    // Valid parents follow the original construction path unchanged.
    // --------------------------------------------------------------
    auto prismPreflightFail =
    [&]
    (
        const char* reason,
        const label d0,
        const label d1,
        const label d2,
        const label d3
    ) -> bool
    {
        Info
            << "CFMITCH V3.5 PRISM PREFLIGHT FAIL:"
            << " cell=" << cellI
            << " reason=" << reason
            << " d0=" << d0
            << " d1=" << d1
            << " d2=" << d2
            << " d3=" << d3
            << endl;

        cellsFromCell.clear();
        return false;
    };

    if
    (
        cellI < 0
     || cellI >= label(mesh_.cells().size())
    )
    {
        return prismPreflightFail
        (
            "cellOutOfRange",
            cellI,
            mesh_.cells().size(),
            -1,
            -1
        );
    }

    if( mesh_.boundaries().size() == 0 )
    {
        return prismPreflightFail
        (
            "noBoundaryPatches",
            -1,
            -1,
            -1,
            -1
        );
    }

    const cell& c = mesh_.cells()[cellI];
    const labelList& owner = mesh_.owner();

    # ifdef DEBUGLayer
    Pout << "New cells from cell " << cellI << endl;
    # endif

    const label startBoundary = mesh_.boundaries()[0].patchStart();

    //- find the number of lyers for this cell
    label nLayers(1), baseFace(-1);
    forAll(c, fI)
    {
        const label bfI = c[fI] - startBoundary;

        if( (bfI < 0) || (bfI >= nLayersAtBndFace_.size()) )
            continue;

        if( nLayersAtBndFace_[bfI] < 2 )
            continue;

        # ifdef DEBUGLayer
        Pout << "Boundary face " << bfI << endl;
        # endif

        nLayers = nLayersAtBndFace_[bfI];
        baseFace = fI;
    }

    if
    (
        baseFace < 0
     || baseFace >= label(c.size())
     || nLayers < 2
    )
    {
        return prismPreflightFail
        (
            "activeBaseFaceMissing",
            baseFace,
            nLayers,
            c.size(),
            -1
        );
    }

    const faceListPMG& preflightFaces =
        mesh_.faces();

    label preflightOtherBaseFace = -1;
    label nOtherBaseCandidates = 0;

    forAll(c, preflightLocalFaceI)
    {
        const label sourceFaceI =
            c[preflightLocalFaceI];

        if
        (
            sourceFaceI < 0
         || sourceFaceI >= label(preflightFaces.size())
         || sourceFaceI >= label(owner.size())
         || sourceFaceI >= label(facesFromFace_.size())
        )
        {
            return prismPreflightFail
            (
                "sourceFaceOutOfRange",
                preflightLocalFaceI,
                sourceFaceI,
                preflightFaces.size(),
                facesFromFace_.size()
            );
        }

        const label nDerived =
            facesFromFace_.sizeOfRow(sourceFaceI);

        if( nDerived <= 0 )
        {
            return prismPreflightFail
            (
                "sourceFaceHasNoDerivedFaces",
                preflightLocalFaceI,
                sourceFaceI,
                nDerived,
                -1
            );
        }

        forAllRow
        (
            facesFromFace_,
            sourceFaceI,
            preflightDerivedI
        )
        {
            const label derivedFaceI =
                facesFromFace_
                (
                    sourceFaceI,
                    preflightDerivedI
                );

            if
            (
                derivedFaceI < 0
             || derivedFaceI >= label(newFaces_.size())
            )
            {
                return prismPreflightFail
                (
                    "derivedFaceOutOfRange",
                    sourceFaceI,
                    preflightDerivedI,
                    derivedFaceI,
                    newFaces_.size()
                );
            }
        }

        if
        (
            preflightLocalFaceI != baseFace
         && nDerived == 1
        )
        {
            preflightOtherBaseFace =
                preflightLocalFaceI;

            ++nOtherBaseCandidates;
        }
    }

    if( nOtherBaseCandidates != 1 )
    {
        return prismPreflightFail
        (
            "oppositeBaseFaceCount",
            nOtherBaseCandidates,
            preflightOtherBaseFace,
            baseFace,
            nLayers
        );
    }

    # ifdef DEBUGLayer
    Pout << "Number of layers " << nLayers << endl;
    Pout << "Base face " << baseFace << " has points "
         << mesh_.faces()[c[baseFace]] << endl;
    forAll(c, fI)
    {
        Pout << "Faces from face " << fI << " are "
             << facesFromFace_[c[fI]] << endl;

        forAllRow(facesFromFace_, c[fI], i)
            Pout << "Face " << facesFromFace_(c[fI], i)
                 << " is " << newFaces_[facesFromFace_(c[fI], i)] << endl;
    }
    # endif

    //- set the number of layers
    cellsFromCell.setSize(nLayers);

    //- distribute existing faces into new cells
    label otherBaseFace(preflightOtherBaseFace);
    forAll(c, fI)
    {
        if( fI == baseFace )
        {
            const label faceI = facesFromFace_(c[fI], 0);
            DynList<label, 8> f;
            f = newFaces_[faceI];
            cellsFromCell[nLayers-1].append(f);
        }
        else if( facesFromFace_.sizeOfRow(c[fI]) == 1 )
        {
            const label faceI = facesFromFace_(c[fI], 0);
            otherBaseFace = fI;
            DynList<label, 8> f;
            f = newFaces_[faceI];
            cellsFromCell[0].append(f);
        }
        else
        {
            forAllRow(facesFromFace_, c[fI], cfI)
            {
                const label nfI = facesFromFace_(c[fI], cfI);

                DynList<label, 8> cf;
                cf = newFaces_[nfI];

                if( owner[c[fI]] != cellI )
                    cf = help::reverseFace(cf);

                cellsFromCell[Foam::max(nLayers-1-cfI, 0)].append(cf);
            }
        }
    }

    //- generate missing faces
    const faceListPMG& faces = mesh_.faces();
    const face& bf = faces[c[baseFace]];
    const face& obf = faces[c[otherBaseFace]];

    if( bf.size() < 3 || obf.size() < 3 )
    {
        return prismPreflightFail
        (
            "baseFaceTooSmall",
            bf.size(),
            obf.size(),
            baseFace,
            otherBaseFace
        );
    }

    for(label layerI=1;layerI<nLayers;++layerI)
    {
        //- create new face from points at the same height
        DynList<label, 8> cf;
        forAll(bf, pI)
        {
            const label pointI = bf[pI];

            # ifdef DEBUGLayer
            Pout << "Split edges at point " << pointI << " are "
                 << splitEdgesAtPoint_[pointI] << endl;
            # endif

            if
            (
                pointI < 0
             || pointI >= label(splitEdgesAtPoint_.size())
            )
            {
                return prismPreflightFail
                (
                    "basePointOutOfRange",
                    pointI,
                    splitEdgesAtPoint_.size(),
                    layerI,
                    pI
                );
            }

            label seI(-1);
            if( splitEdgesAtPoint_.sizeOfRow(pointI) == 1 )
            {
                seI = splitEdgesAtPoint_(pointI, 0);
            }
            else
            {
                forAllRow(splitEdgesAtPoint_, pointI, sepI)
                {
                    const label seJ =
                        splitEdgesAtPoint_(pointI, sepI);

                    if
                    (
                        seJ < 0
                     || seJ >= label(splitEdges_.size())
                    )
                        continue;

                    const edge& se = splitEdges_[seJ];

                    if( obf.which(se.end()) >= 0 || obf.which(se.start()) >= 0 )
                    {
                        seI = seJ;
                        break;
                    }
                }
            }

            if
            (
                seI < 0
             || seI >= label(splitEdges_.size())
             || seI >= label(newVerticesForSplitEdge_.size())
            )
            {
                return prismPreflightFail
                (
                    "splitEdgeUnresolved",
                    pointI,
                    layerI,
                    seI,
                    splitEdgesAtPoint_.sizeOfRow(pointI)
                );
            }

            if
            (
                newVerticesForSplitEdge_.sizeOfRow(seI)
             <= layerI
            )
            {
                return prismPreflightFail
                (
                    "splitEdgeRowTooShort",
                    seI,
                    layerI,
                    newVerticesForSplitEdge_.sizeOfRow(seI),
                    nLayers
                );
            }

            const label generatedPointI =
                newVerticesForSplitEdge_
                (
                    seI,
                    layerI
                );

            if
            (
                generatedPointI < 0
             || generatedPointI >= label(mesh_.points().size())
            )
            {
                return prismPreflightFail
                (
                    "generatedPointOutOfRange",
                    seI,
                    layerI,
                    generatedPointI,
                    mesh_.points().size()
                );
            }

            cf.append(generatedPointI);
        }

        // C3: skip degenerate intermediate face where all vertices are
        // the same point (zero-length hair edge at BL/BL junction).
        // This produces correct wedge topology instead of zero-volume cells.
        bool degenerateFace = true;
        for(label pI=1;pI<cf.size();++pI)
            if( cf[pI] != cf[0] )
            { degenerateFace = false; break; }
        if( degenerateFace )
        {
            //- DEGFACE_AUDIT: diagnostic only, no behavior change
            static label nDegFaceSkips = 0;
            if( nDegFaceSkips < 10 )
            {
                const label baseGlobalFaceI = c[baseFace];
                const label startBnd = mesh_.boundaries()[0].patchStart();
                const label bfI_audit = baseGlobalFaceI - startBnd;
                Info << "DEGFACE_AUDIT:"
                     << " cellI=" << cellI
                     << " layerI=" << layerI
                     << " nLayers=" << nLayers
                     << " baseFace=" << baseFace
                     << " bfI=" << bfI_audit
                     << " cf0=" << cf[0]
                     << " cfSize=" << cf.size()
                     << endl;
            }
            ++nDegFaceSkips;
            if( nDegFaceSkips == 10 )
                Info << "DEGFACE_AUDIT: (further suppressed)" << endl;
            continue;
        }

        //- add faces to cells
        cellsFromCell[nLayers-layerI].append(cf);
        cellsFromCell[nLayers-1-layerI].append(cf);
    }

    # ifdef DEBUGLayer
    Pout << "New cells from cell " << cellI << " are " << cellsFromCell << endl;
    //::exit(1);

    Pout << "1. Newly generated cells " << cellsFromCell << endl;

    //- check if all generated cells are topologically closed
    forAll(cellsFromCell, cI)
    {
        const DynList<DynList<label, 8>, 10>& cellFaces = cellsFromCell[cI];

        DynList<edge, 12> edges;
        DynList<label, 12> nAppearances;

        forAll(cellFaces, fI)
        {
            const DynList<label, 8>& f = cellFaces[fI];

            forAll(f, eI)
            {
                const edge e(f[eI], f.fcElement(eI));

                const label pos = edges.containsAtPosition(e);

                if( pos < 0 )
                {
                    edges.append(e);
                    nAppearances.append(1);
                }
                else
                {
                    ++nAppearances[pos];
                }
            }
        }

        forAll(nAppearances, eI)
            if( nAppearances[eI] != 2 )
            {
                Pout << "Prism cell " << cI << " edge " << edges[eI]
                    << " is present " << nAppearances[eI] << " times!" << endl;
                abort(FatalError);
            }
    }
    # endif

    return true;
}

void refineBoundaryLayers::storeFacesIntoCells
(
    const label faceI,
    const bool reverseOrientation,
    const label normalDirection,
    const bool maxCoordinate,
    const label nLayersI,
    const label nLayersJ,
    const label nLayersK,
    DynList<DynList<DynList<label, 4>, 6>, 256>& cellsFromCell
) const
{
    DynList<DynList<label> > faceFaces;
    sortFaceFaces(faceI, faceFaces, reverseOrientation);

    const label maxI = nLayersI - 1;
    const label maxJ = nLayersJ - 1;
    const label maxK = nLayersK - 1;

    # ifdef DEBUGLayer
    Pout << "Storing new faces from face " << faceI
         << " reverseOrientation = " << reverseOrientation
         << " normal direction " << normalDirection
         << " maxCoordinate " << maxCoordinate << endl;
    Pout << "faceFaces " << faceFaces << endl;
    # endif

    label i(-1), j(-1), k(-1);

    forAll(faceFaces, nI)
    {
        forAll(faceFaces[nI], nJ)
        {
            const label nfI = faceFaces[nI][nJ];

            # ifdef DEBUGLayer
            Pout << "nI = " << nI << " nJ = " << nJ << endl;
            # endif

            if( normalDirection == 0 )
            {
                //- k is const
                i = Foam::min(nI, maxI);
                j = Foam::min(nJ, maxJ);
                k = maxCoordinate?maxK:0;
            }
            else if( normalDirection == 1 )
            {
                //- j is const
                i = Foam::min(nJ, maxI);
                j = maxCoordinate?maxJ:0;
                k = Foam::min(nI, maxK);
            }
            else if( normalDirection == 2 )
            {
                //- i is const
                i = maxCoordinate?maxI:0;
                j = Foam::min(nI, maxJ);
                k = Foam::min(nJ, maxK);
            }

            //- store the face into a new cell
            const label cI
            (
                j * nLayersI +
                k * nLayersI * nLayersJ +
                i
            );

            # ifdef DEBUGLayer
            Pout << "Storing face " << newFaces_[nfI]
                 << " i = " << i << " j = " << j << " k = " << k
                 << "\n cell label " << cI << endl;
            # endif

            cellsFromCell[cI].append(newFaces_[nfI]);
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void refineBoundaryLayers::refineEdgeHexCell::determineFacesInDirections()
{
    const labelList& nLayersAtBndFace = bndLayers_.nLayersAtBndFace_;
    const polyMeshGen& mesh = bndLayers_.mesh_;
    const faceListPMG& faces = mesh.faces();
    const cell& c = mesh.cells()[cellI_];

    # ifdef DEBUGLayer
    Pout << "Generating new cells from edge cell " << cellI_ << endl;
    # endif

    const PtrList<boundaryPatch>& bnd = mesh.boundaries();
    const label startBoundary = bnd[0].patchStart();

    //- find the number of layers for this cell
    FixedList<label, 2> layersInDirection(-1), dirFace;
    label currDir(0);

    FixedList<bool, 6> determinedFace(false);

    forAll(c, fI)
    {
        const label bfI = c[fI] - startBoundary;

        if( (bfI < 0) || (bfI >= nLayersAtBndFace.size()) )
            continue;

        # ifdef DEBUGLayer
        Pout << "Boundary face " << bfI << endl;
        # endif

        if( nLayersAtBndFace[bfI] < 2 )
            continue;

        layersInDirection[currDir] = nLayersAtBndFace[bfI];
        dirFace[currDir] = fI;
        ++currDir;
    }

    //- set the number of newly create cells
    nLayersI_ = layersInDirection[0];
    nLayersJ_ = layersInDirection[1];
    cellsFromCell_.setSize(nLayersI_ * nLayersJ_);

    //- find the shared edge between the boundary faces
    const edge commonEdge =
        help::sharedEdge(faces[c[dirFace[0]]], faces[c[dirFace[1]]]);

    //- faces at i = const in the local coordinate system
    faceInDirection_[4] = dirFace[0];
    determinedFace[dirFace[0]] = true;
    forAll(c, fI)
    {
        if( determinedFace[fI] )
            continue;

        if( !help::shareAnEdge(faces[c[dirFace[0]]], faces[c[fI]]) )
        {
            faceInDirection_[5] = fI;
            determinedFace[fI] = true;
            break;
        }
    }

    //- faces k = const in the local coordinate system
    faceInDirection_[2] = dirFace[1];
    determinedFace[dirFace[1]] = true;
    forAll(c, fI)
    {
        if( determinedFace[fI] )
            continue;

        if( !help::shareAnEdge(faces[c[dirFace[1]]], faces[c[fI]]) )
        {
            faceInDirection_[3] = fI;
            determinedFace[fI] = true;
            break;
        }
    }

    # ifdef DEBUGLayer
    Pout << "Common edge " << commonEdge << endl;
    Pout << "Donor face " << dirFace[0] << endl;
    Pout << "Donor face points " << faces[c[dirFace[0]]] << endl;
    # endif

    //- find the face attached to the starting point of the edge and
    //- the face attached to the end point of the edge
    forAll(c, fI)
    {
        if( determinedFace[fI] )
            continue;

        if(
            (faces[c[fI]].which(commonEdge.start()) >= 0) &&
            (help::positionOfEdgeInFace(commonEdge, faces[c[fI]]) < 0)
        )
            faceInDirection_[0] = fI;

        if(
            (faces[c[fI]].which(commonEdge.end()) >= 0) &&
            (help::positionOfEdgeInFace(commonEdge, faces[c[fI]]) < 0)
        )
            faceInDirection_[1] = fI;
    }

    //- check the orientation of faces
    const labelList& owner = mesh.owner();

    //- checking face at direction k = 0
    faceOrientation_[0] = owner[c[faceInDirection_[0]]] == cellI_?true:false;

    //- checking face in direction k = 1
    faceOrientation_[1] = owner[c[faceInDirection_[1]]] == cellI_?false:true;

    //- set orientation flag for face in direction j = 0
    faceOrientation_[2] = true;

    //- checking face in direction j = nLayersJ_
    faceOrientation_[3] = owner[c[faceInDirection_[3]]] == cellI_?false:true;

    //- set orientation flag for face in direction i = 0
    faceOrientation_[4] = true;

    //- checking face in direction i = nLayersI_
    faceOrientation_[5] = owner[c[faceInDirection_[5]]] == cellI_?false:true;

    # ifdef DEBUGLayer
    Pout << "Face at start " << faces[c[faceInDirection_[0]]] << endl;
    Pout << "Face at end " << faces[c[faceInDirection_[1]]] << endl;
    forAll(faceInDirection_, i)
        Pout << "Face in direction " << i << " is "
             << faces[c[faceInDirection_[i]]]
             << " orientation " << faceOrientation_[i] << endl;
    # endif
}

void refineBoundaryLayers::refineEdgeHexCell::populateExistingFaces()
{
    const cell& c = bndLayers_.mesh_.cells()[cellI_];
    const VRWGraph& facesFromFace = bndLayers_.facesFromFace_;
    const VRWGraph& newFaces = bndLayers_.newFaces_;

    cellsFromCell_.setSize(nLayersI_ * nLayersJ_);
    forAll(cellsFromCell_, cI)
        cellsFromCell_[cI].clear();

    //- store new faces at k = 0
    bndLayers_.storeFacesIntoCells
    (
        c[faceInDirection_[0]], faceOrientation_[0],
        0, 0,
        nLayersI_, nLayersJ_, 1,
        cellsFromCell_
    );

    //- store new faces at k = 1
    bndLayers_.storeFacesIntoCells
    (
        c[faceInDirection_[1]], faceOrientation_[1],
        0, 1,
        nLayersI_, nLayersJ_, 1,
        cellsFromCell_
    );

    //- store new faces at j = 0
    forAllRow(facesFromFace, c[faceInDirection_[2]], i)
    {
        const label faceI = facesFromFace(c[faceInDirection_[2]], i);
        cellsFromCell_[i].append(newFaces[faceI]);
    }

    //- store faces at j = nLayersJ
    const label maxJ = nLayersJ_ - 1;
    forAllRow(facesFromFace, c[faceInDirection_[3]], i)
    {
        const label faceI = facesFromFace(c[faceInDirection_[3]], i);
        cellsFromCell_[i + maxJ * nLayersI_].append(newFaces[faceI]);
    }

    //- store new faces at i = 0
    forAllRow(facesFromFace, c[faceInDirection_[4]], j)
    {
        const label faceI = facesFromFace(c[faceInDirection_[4]], j);
        cellsFromCell_[j * nLayersI_].append(newFaces[faceI]);
    }

    //- store new faces at i = nLayersI
    const label maxI = nLayersI_ - 1;
    forAllRow(facesFromFace, c[faceInDirection_[5]], j)
    {
        const label faceI = facesFromFace(c[faceInDirection_[5]], j);
        cellsFromCell_[j * nLayersI_ + maxI].append(newFaces[faceI]);
    }

    # ifdef DEBUGLayer
    Pout << "New cells after populating existing faces "
         << cellsFromCell_ << endl;
    # endif
}

void refineBoundaryLayers::refineEdgeHexCell::generateMissingFaces()
{
    const cell& c = bndLayers_.mesh_.cells()[cellI_];

    //- fill up the matrix of points for this cell
    //- the matrix is used for generation of new cells
    FixedList<DynList<DynList<label> >, 2> cellPoints;

    //- fill in the data for a cross-split faces
    bndLayers_.sortFacePoints
    (
        c[faceInDirection_[0]],
        cellPoints[0],
        faceOrientation_[0]
    );
    bndLayers_.sortFacePoints
    (
        c[faceInDirection_[1]],
        cellPoints[1],
        faceOrientation_[1]
    );

    //- generate new internal faces for this cell
    //- generate faces with normal in the i direction
    const label maxI = nLayersI_ - 1;
    const label maxJ = nLayersJ_ - 1;

    for(label i=1;i<nLayersI_;++i)
    {
        for(label j=0;j<nLayersJ_;++j)
        {
            const label own = j * nLayersI_ + i - 1;
            const label nei = own + 1;

            if( j < maxJ )
            {
                //- generate a quad face
                FixedList<label, 4> mf;

                //- populate the points form cellPoints
                mf[0] = cellPoints[0][i][j];
                mf[1] = cellPoints[0][i][j+1];
                mf[2] = cellPoints[1][i][j+1];
                mf[3] = cellPoints[1][i][j];

                # ifdef DEBUGLayer
                Pout << "1. Adding missing face " << mf
                     << " to cells " << own << " and " << nei << endl;
                # endif

                cellsFromCell_[own].append(mf);
                cellsFromCell_[nei].append(help::reverseFace(mf));
            }
            else
            {
                DynList<label> mf;
                for(label index=j;index<cellPoints[0][i].size();++index)
                    mf.append(cellPoints[0][i][index]);
                for(label index=cellPoints[1][i].size()-1;index>=j;--index)
                    mf.append(cellPoints[1][i][index]);

                # ifdef DEBUGLayer
                Pout << "2. Adding missing face " << mf
                     << " to cells " << own << " and " << nei << endl;
                # endif

                cellsFromCell_[own].append(mf);
                cellsFromCell_[nei].append(help::reverseFace(mf));
            };
        }
    }

    //- generate faces with the normal in j direction
    for(label i=0;i<nLayersI_;++i)
    {
        for(label j=1;j<nLayersJ_;++j)
        {
            const label nei = j * nLayersI_ + i;
            const label own = (j - 1) * nLayersI_ + i;

            if( i < maxI )
            {
                //- generate a quad face
                FixedList<label, 4> mf;

                //- populate the points form cellPoints
                mf[0] = cellPoints[0][i][j];
                mf[1] = cellPoints[1][i][j];
                mf[2] = cellPoints[1][i+1][j];
                mf[3] = cellPoints[0][i+1][j];

                # ifdef DEBUGLayer
                Pout << "3. Adding missing face " << mf
                     << " to cells " << own << " and " << nei << endl;
                # endif

                cellsFromCell_[own].append(mf);
                cellsFromCell_[nei].append(help::reverseFace(mf));
            }
            else
            {
                DynList<label> mf;
                for(label index=i;index<cellPoints[1].size();++index)
                    mf.append(cellPoints[1][index][j]);
                for(label index=cellPoints[0].size()-1;index>=i;--index)
                    mf.append(cellPoints[0][index][j]);

                # ifdef DEBUGLayer
                Pout << "4. Adding missing face " << mf
                     << " to cells " << own << " and " << nei << endl;
                # endif

                cellsFromCell_[own].append(mf);
                cellsFromCell_[nei].append(help::reverseFace(mf));
            };
        }
    }

    # ifdef DEBUGLayer
    Pout << "Cell " << cellI_ << " new cells are " << cellsFromCell_ << endl;
    //::exit(1);
    # endif
}

refineBoundaryLayers::refineEdgeHexCell::refineEdgeHexCell
(
    const label cellI,
    const refineBoundaryLayers& ref
)
:
    cellI_(cellI),
    nLayersI_(),
    nLayersJ_(),
    cellsFromCell_(),
    bndLayers_(ref),
    faceInDirection_(),
    faceOrientation_(),
    cellPoints_()
{
    determineFacesInDirections();

    populateExistingFaces();

    generateMissingFaces();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void refineBoundaryLayers::refineCornerHexCell::determineFacesInDirections()
{
    const polyMeshGen& mesh = bndLayers_.mesh_;
    const cell& c = mesh.cells()[cellI_];
    const faceListPMG& faces = mesh.faces();
    const labelList& nLayersAtBndFace = bndLayers_.nLayersAtBndFace_;

    # ifdef DEBUGLayer
    Pout << "Generating new cells from corner hex cell " << cellI_ << endl;
    Pout << "Cell faces " << c << endl;
    # endif

    const label startBoundary = mesh.boundaries()[0].patchStart();

    //- find the number of layers for this cell
    FixedList<label, 3> layersInDirection(-1), dirFace;
    FixedList<bool, 6> usedDirection(false);
    label currDir(0);

    forAll(c, fI)
    {
        const label bfI = c[fI] - startBoundary;

        if( (bfI < 0) || (bfI >= nLayersAtBndFace.size()) )
            continue;

        # ifdef DEBUGLayer
        Pout << "Boundary face " << bfI << endl;
        # endif

        if( nLayersAtBndFace[bfI] < 2 )
            continue;

        usedDirection[fI] = true;
        layersInDirection[currDir] = nLayersAtBndFace[bfI];
        dirFace[currDir] = fI;
        ++currDir;
    }

    //- find a common point for all three boundary faces
    FixedList<DynList<label, 4>, 3> bndFaces;
    forAll(dirFace, i)
    {
        bndFaces[i] = faces[c[dirFace[i]]];
    }

    const label commonPoint = help::sharedVertex(bndFaces);

    # ifdef DEBUGLayer
    Pout << "Used directions " << usedDirection << endl;
    Pout << "Layers in direction " << layersInDirection << endl;
    Pout << "dirFace " << dirFace << endl;
    Pout << "Common point " << commonPoint << endl;

    forAll(dirFace, i)
        Pout << "bnd face " << i << " is " << faces[c[dirFace[i]]] << endl;
    # endif

    //- find the position of the common point in each boundary face
    const edgeLongList& splitEdges = bndLayers_.splitEdges_;
    const VRWGraph& splitEdgesAtPoint = bndLayers_.splitEdgesAtPoint_;

    const face& baseFace = faces[c[dirFace[0]]];
    const label posInBndFace = baseFace.which(commonPoint);

    //- find split edges starting at the commonPoints
    forAllRow(splitEdgesAtPoint, commonPoint, i)
    {
        const edge& se = splitEdges[splitEdgesAtPoint(commonPoint, i)];

        if( se == baseFace.faceEdge(posInBndFace) )
        {
            //- this edge is in j direction
            splitEdgeInDirection_[1] = splitEdgesAtPoint(commonPoint, i);
        }
        else if( se == baseFace.faceEdge(baseFace.rcIndex(posInBndFace)) )
        {
            //- this edge is in i diretion
            splitEdgeInDirection_[0] = splitEdgesAtPoint(commonPoint, i);
        }
        else if( splitEdgesAtPoint.sizeOfRow(commonPoint) == 3 )
        {
            //- this point is in k direction
            splitEdgeInDirection_[2] = splitEdgesAtPoint(commonPoint, i);
        }
        else
        {
            //- this situation is not allowed
            FatalErrorIn
            (
                "void refineBoundaryLayers::refineCornerHexCell::"
                "determineFacesInDirections()"
            ) << "Cannot refine layer for cell " << cellI_ << abort(FatalError);
        }
    }

    # ifdef DEBUGLayer
    const VRWGraph& newVerticesForSplitEdge =
        bndLayers_.newVerticesForSplitEdge_;
    forAll(splitEdgeInDirection_, i)
        Pout << "Split edge in direction " << i << " has nodes "
             << splitEdges[splitEdgeInDirection_[i]]
             << " number of points on split edge "
             << newVerticesForSplitEdge.sizeOfRow(splitEdgeInDirection_[i])
             << endl;
    # endif

    //- find the direction od other boundary faces
    //- in the local coordinate system
    FixedList<label, 3> permutation;
    permutation[0] = 0;

    label helper = help::positionOfEdgeInFace
    (
        baseFace.faceEdge(baseFace.rcIndex(posInBndFace)),
        faces[c[dirFace[1]]]
    );

    if( helper >= 0 )
    {
        permutation[1] = 1;
        permutation[2] = 2;
    }
    else
    {
        permutation[1] = 2;
        permutation[2] = 1;
    }

    //- find the number of layers and a split in each direction
    nLayersI_ = layersInDirection[permutation[2]];
    nLayersJ_ = layersInDirection[permutation[1]];
    nLayersK_ = layersInDirection[permutation[0]];

    //- determine the directions of cell faces
    //- store boundary faces first. Their normals point in the wrong direction
    //- face at k = 0
    faceInDirection_[0] = dirFace[permutation[0]];
    faceOrientation_[0] = true;
    //- face at j = 0
    faceInDirection_[2] = dirFace[permutation[1]];
    faceOrientation_[2] = true;
    //- face at i = 0
    faceInDirection_[4] = dirFace[permutation[2]];
    faceOrientation_[4] = true;

    //- find directions of other faces and thrie orientation
    const labelList& owner = mesh.owner();
    forAll(c, fI)
    {
        if( usedDirection[fI] )
            continue;

        const bool orientation = owner[c[fI]]==cellI_?false:true;

        if( !help::shareAnEdge(faces[c[fI]], faces[c[faceInDirection_[0]]]) )
        {
            //- face at k = nLayersK_
            faceInDirection_[1] = fI;
            faceOrientation_[1] = orientation;
        }
        else if
        (
            !help::shareAnEdge(faces[c[fI]], faces[c[faceInDirection_[2]]])
        )
        {
            //- face at j = nLayersJ_
            faceInDirection_[3] = fI;
            faceOrientation_[3] = orientation;
        }
        else if
        (
            !help::shareAnEdge(faces[c[fI]], faces[c[faceInDirection_[4]]])
        )
        {
            //- face at i = nLayersI_
            faceInDirection_[5] = fI;
            faceOrientation_[5] = orientation;
        }
    }

    # ifdef DEBUGLayer
    forAll(faceInDirection_, i)
        Pout << "Face in direction " << i
             << " is " << faces[c[faceInDirection_[i]]]
             << " orientation " << faceOrientation_[i] << endl;
    Pout << "nLayersI = " << nLayersI_
         << " nLayersJ = " << nLayersJ_
         << " nLayersK = " << nLayersK_ << endl;
    # endif
}

void refineBoundaryLayers::refineCornerHexCell::populateExistingFaces()
{
    const cell& c = bndLayers_.mesh_.cells()[cellI_];

    //- set the number of cells
    cellsFromCell_.setSize(nLayersI_ * nLayersJ_ * nLayersK_);
    forAll(cellsFromCell_, i)
        cellsFromCell_[i].clear();

    //- add new faces from existing faces into new cells
    forAll(faceInDirection_, dirI)
    {
        bndLayers_.storeFacesIntoCells
        (
            c[faceInDirection_[dirI]], faceOrientation_[dirI],
            dirI / 2, dirI % 2,
            nLayersI_, nLayersJ_, nLayersK_,
            cellsFromCell_
        );
    }

    # ifdef DEBUGLayer
    Pout << "cellsFromCell_ before new faces " << cellsFromCell_ << endl;
    //::exit(1);
    # endif
}

void refineBoundaryLayers::refineCornerHexCell::generateNewPoints()
{
    const cell& c = bndLayers_.mesh_.cells()[cellI_];

    //- allocate space for points generated inside the cell
    cellPoints_.setSize(nLayersI_+1);
    forAll(cellPoints_, i)
    {
        cellPoints_[i].setSize(nLayersJ_+1);

        forAll(cellPoints_[i], j)
        {
            cellPoints_[i][j].setSize(nLayersK_+1);
            cellPoints_[i][j] = -1;
        }
    }

    //- collect information about points generated on faces of the cell
    forAll(faceInDirection_, dirI)
    {
        bndLayers_.sortFacePoints
        (
            c[faceInDirection_[dirI]],
            facePoints_[dirI],
            faceOrientation_[dirI]
        );
    }

    # ifdef DEBUGLayer
    Pout << "Face points " << facePoints_ << endl;
    # endif

    //- fill in cellPoints at the boundary
    forAll(cellPoints_, i)
    {
        forAll(cellPoints_[i], j)
        {
            cellPoints_[i][j][0] = facePoints_[0][i][j];
            cellPoints_[i][j][nLayersK_] = facePoints_[1][i][j];
        }
    }

    forAll(cellPoints_, i)
    {
        forAll(cellPoints_[i][0], k)
        {
            cellPoints_[i][0][k] = facePoints_[2][k][i];
            cellPoints_[i][nLayersJ_][k] = facePoints_[3][k][i];
        }
    }

    forAll(cellPoints_[0], j)
    {
        forAll(cellPoints_[0][j], k)
        {
            cellPoints_[0][j][k] = facePoints_[4][j][k];
            cellPoints_[nLayersI_][j][k] = facePoints_[5][j][k];
        }
    }

    //- useful data for generating missing points
    const edgeLongList& splitEdges = bndLayers_.splitEdges_;
    const edge& seDirI = splitEdges[splitEdgeInDirection_[0]];
    const edge& seDirJ = splitEdges[splitEdgeInDirection_[1]];
    const edge& seDirK = splitEdges[splitEdgeInDirection_[2]];
    const VRWGraph& ptsAtEdge = bndLayers_.newVerticesForSplitEdge_;

    //- const references to vertices of the cell ordered in a local
    //- i, j, k coordinate system
    pointFieldPMG& points = bndLayers_.mesh_.points();
    const point v000 = points[seDirI.start()];
    const point v100 = points[seDirI.end()];
    const point v110 = points[facePoints_[0].lastElement().lastElement()];
    const point v010 = points[seDirJ.end()];
    const point v001 = points[seDirK.end()];
    const point v101 = points[facePoints_[1].lastElement()[0]];
    const point v111 = points[facePoints_[1].lastElement().lastElement()];
    const point v011 = points[facePoints_[1][0].lastElement()];

    for(label i=1;i<nLayersI_;++i)
    {
        const scalar u
        (
            Foam::mag
            (
                points[ptsAtEdge(splitEdgeInDirection_[0], i)] -
                points[seDirI.start()]
            ) /
            seDirI.mag(points)
        );

        for(label j=1;j<nLayersJ_;++j)
        {
            const scalar v
            (
                Foam::mag
                (
                    points[ptsAtEdge(splitEdgeInDirection_[1], j)] -
                    points[seDirJ.start()]
                ) /
                seDirJ.mag(points)
            );

            for(label k=1;k<nLayersK_;++k)
            {
                const scalar w
                (
                    Foam::mag
                    (
                        points[ptsAtEdge(splitEdgeInDirection_[2], k)] -
                        points[seDirK.start()]
                    ) /
                    seDirK.mag(points)
                );

                # ifdef DEBUGLayer
                Pout << "Generating point in corner cell local coordinates "
                     << "u = " << u << " v = " << v << " w = " << w << endl;
                # endif

                //- calculate coordinates of the new vertex
                const point newP =
                    (1.0 - u) * (1.0 - v) * (1.0 - w) * v000 +
                    u * (1.0 - v) * (1.0 - w) * v100 +
                    u * v * (1.0 - w) * v110 +
                    (1.0 - u) * v * (1.0 - w) * v010 +
                    (1.0 - u) * (1.0 - v) * w * v001 +
                    u * (1.0 - v) * w * v101 +
                    u * v * w * v111 +
                    (1.0 - u) * v * w * v011;

                # ifdef DEBUGLayer
                Pout << "New point " << points.size() << " in corner hex "
                    << "has coordinates " << newP << endl;
                # endif

                //- add the point to the mesh
                cellPoints_[i][j][k] = points.size();
                points.append(newP);
            }
        }
    }

    # ifdef DEBUGLayer
    Pout << "New cell points " << cellPoints_ << endl;
    //::exit(1);
    # endif
}

void refineBoundaryLayers::refineCornerHexCell::generateMissingFaces()
{
    //- generate face in direction i
    for(label i=1;i<nLayersI_;++i)
    {
        //- generate quad faces
        for(label j=0;j<nLayersJ_;++j)
        {
            for(label k=0;k<nLayersK_;++k)
            {
                //- skip generating last face because it might not be a quad
                if( (j == (nLayersJ_-1)) && (k == (nLayersK_-1)) )
                    continue;

                const label own
                (
                    k * nLayersI_ * nLayersJ_ +
                    j * nLayersI_ +
                    i - 1
                );
                const label nei = own + 1;

                FixedList<label, 4> mf;

                mf[0] = cellPoints_[i][j][k];
                mf[1] = cellPoints_[i][j+1][k];
                mf[2] = cellPoints_[i][j+1][k+1];
                mf[3] = cellPoints_[i][j][k+1];

                cellsFromCell_[own].append(mf);
                cellsFromCell_[nei].append(help::reverseFace(mf));
            }
        }

        //- generate faces which might not be a quads
        DynList<label> mf;

        mf.append(cellPoints_[i][nLayersJ_-1][nLayersK_-1]);

        //- this face might not be a quad
        //- add points fom the last face in direction j
        const DynList<DynList<label> >& f3 = facePoints_[3];
        for(label index=nLayersK_-1;index<f3.size()-1;++index)
            mf.append(f3[index][i]);

        //- add points from the last face in direction k
        const DynList<DynList<label> >& f1 = facePoints_[1];
        for(label index=f1[i].size()-1;index>=nLayersJ_-1;--index)
            mf.append(f1[i][index]);

        const label own
        (
            (nLayersK_-1) * nLayersI_ * nLayersJ_ +
            (nLayersJ_-1) * nLayersI_ +
            i - 1
        );

        const label nei = own + 1;

        # ifdef DEBUGLayer
        Pout << "Additional face in direction i = " << i
             << " j = " << (nLayersJ_-1)
             << " has owner " << own
             << " neighbour " << nei << " with nodes " << mf << endl;
        # endif

        cellsFromCell_[own].append(mf);
        cellsFromCell_[nei].append(help::reverseFace(mf));
    }

    //- generate faces in direction j
    for(label j=1;j<nLayersJ_;++j)
    {
        //- generate quad faces
        for(label i=0;i<nLayersI_;++i)
        {
            for(label k=0;k<nLayersK_;++k)
            {
                //- skip generating late face because it might not be a quad
                if( (i == (nLayersI_-1)) && (k == (nLayersK_-1)) )
                    continue;

                const label own
                (
                    k * nLayersI_ * nLayersJ_ +
                    (j-1) * nLayersI_ +
                    i
                );

                const label nei
                (
                    k * nLayersI_ * nLayersJ_ +
                    j * nLayersI_ +
                    i
                );

                FixedList<label, 4> mf;

                mf[0] = cellPoints_[i][j][k];
                mf[1] = cellPoints_[i][j][k+1];
                mf[2] = cellPoints_[i+1][j][k+1];
                mf[3] = cellPoints_[i+1][j][k];

                cellsFromCell_[own].append(mf);
                cellsFromCell_[nei].append(help::reverseFace(mf));
            }
        }

        //- generate a face which might not be a quad
        DynList<label> mf;

        mf.append(cellPoints_[nLayersI_-1][j][nLayersK_-1]);

        //- add points from the last face in direction k
        const DynList<DynList<label> >& fp1 = facePoints_[1];
        for(label index=nLayersI_-1;index<fp1.size()-1;++index)
            mf.append(fp1[index][j]);

        //- add points from the last face in direction i
        const DynList<DynList<label> >& fp5 = facePoints_[5];
        for(label index=fp5[j].size()-1;index>=nLayersK_-1;--index)
            mf.append(fp5[j][index]);

        const label own
        (
            (nLayersK_-1) * nLayersI_ * nLayersJ_ +
            (j-1) * nLayersI_ +
            (nLayersI_ - 1)
        );

        const label nei
        (
            (nLayersK_-1) * nLayersI_ * nLayersJ_ +
            j * nLayersI_ +
            (nLayersI_ - 1)
        );

        # ifdef DEBUGLayer
        Pout << "Additional face at i = " << (nLayersI_-1)
             << " j = " << j << " k = " << (nLayersK_-1)
             << " has owner " << own
             << " neighbour " << nei << " with nodes " << mf << endl;
        # endif

        cellsFromCell_[own].append(mf);
        cellsFromCell_[nei].append(help::reverseFace(mf));
    }

    //- generate faces in direction k
    for(label k=1;k<nLayersK_;++k)
    {
        //- generate quad faces
        for(label i=0;i<nLayersI_;++i)
        {
            for(label j=0;j<nLayersJ_;++j)
            {
                //- skip the last face because it might not be a quad
                if( (i == (nLayersI_-1)) && (j == (nLayersJ_-1)) )
                    continue;

                const label own
                (
                    (k-1) * nLayersI_ * nLayersJ_ +
                    j * nLayersI_ +
                    i
                );

                const label nei
                (
                    k * nLayersI_ * nLayersJ_ +
                    j * nLayersI_ +
                    i
                );

                FixedList<label, 4> mf;

                mf[0] = cellPoints_[i][j][k];
                mf[1] = cellPoints_[i+1][j][k];
                mf[2] = cellPoints_[i+1][j+1][k];
                mf[3] = cellPoints_[i][j+1][k];

                cellsFromCell_[own].append(mf);
                cellsFromCell_[nei].append(help::reverseFace(mf));
            }
        }

        //- generate a face which might not be a quad
        DynList<label> mf;

        mf.append(cellPoints_[nLayersI_-1][nLayersJ_-1][k]);

        //- this face might not be a quad
        //- add points from the last face in direction i
        const DynList<DynList<label> >& fp5 = facePoints_[5];
        for(label index=nLayersJ_-1;index<fp5.size()-1;++index)
            mf.append(fp5[index][k]);

        //- add points from the last face in direction j
        const DynList<DynList<label> >& fp3 = facePoints_[3];
        for(label index=fp3[k].size()-1;index>=nLayersI_-1;--index)
            mf.append(fp3[k][index]);

        const label own
        (
            (k-1) * nLayersI_ * nLayersJ_ +
            (nLayersJ_-1) * nLayersI_ +
            (nLayersI_ - 1)
        );

        const label nei
        (
            k * nLayersI_ * nLayersJ_ +
            (nLayersJ_-1) * nLayersI_ +
            (nLayersI_ - 1)
        );

        # ifdef DEBUGLayer
        Pout << "Additional face at position i = " << (nLayersI_-1)
             << " j = " << (nLayersJ_-1) << " k = " << k
             << " has owner " << own
             << " neighbour " << nei << " with nodes " << mf << endl;
        # endif

        cellsFromCell_[own].append(mf);
        cellsFromCell_[nei].append(help::reverseFace(mf));
    }

    # ifdef DEBUGLayer
    Pout << "Generated cells " << cellsFromCell_ << endl;

    forAll(cellsFromCell_, cI)
    {
        const DynList<DynList<label, 4>, 6>& cellFaces = cellsFromCell_[cI];

        DynList<edge, 12> edges;
        DynList<label, 12> nAppearances;

        forAll(cellFaces, fI)
        {
            const DynList<label, 4>& f = cellFaces[fI];

            forAll(f, eI)
            {
                const edge e(f[eI], f.fcElement(eI));

                const label pos = edges.containsAtPosition(e);

                if( pos < 0 )
                {
                    edges.append(e);
                    nAppearances.append(1);
                }
                else
                {
                    ++nAppearances[pos];
                }
            }
        }

        forAll(nAppearances, eI)
            if( nAppearances[eI] != 2 )
            {
                Pout << "Edge hex cell " << cI << " edge " << edges[eI]
                    << " is present " << nAppearances[eI] << " times!" << endl;
                abort(FatalError);
            }
    }

    //::exit(1);
    # endif
}

refineBoundaryLayers::refineCornerHexCell::refineCornerHexCell
(
    const label cellI,
    const refineBoundaryLayers& ref
)
:
    cellI_(cellI),
    nLayersI_(),
    nLayersJ_(),
    nLayersK_(),
    splitEdgeInDirection_(),
    cellsFromCell_(),
    bndLayers_(ref),
    faceInDirection_(),
    faceOrientation_(),
    facePoints_(),
    cellPoints_()
{
    determineFacesInDirections();

    populateExistingFaces();

    generateNewPoints();

    generateMissingFaces();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void refineBoundaryLayers::generateNewCells()
{
    labelList nCellsFromCell(mesh_.cells().size(), 1);
    labelList refType(mesh_.cells().size(), 0);
    labelList cellToBfI(mesh_.cells().size(), -1);

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& faceOwners = mse.faceOwners();
    const labelList& childSweepFacePatch =
        mse.boundaryFacePatches();

    // Owning copy survives deleteDemandDrivenData(msePtr_) and lets the
    // exact-volume birth audit recover the source BL patch afterwards.
    const labelList exactVolumeFacePatch(childSweepFacePatch);

    const PtrList<boundaryPatch>& childSweepBoundaries =
        mesh_.boundaries();

    //- calculate the number new cells generated from a cell
    forAll(faceOwners, bfI)
    {
        const label cellI = faceOwners[bfI];

        nCellsFromCell[cellI] *= nLayersAtBndFace_[bfI];

        if( nLayersAtBndFace_[bfI] > 1 )
        {
            ++refType[cellI];
            if( cellToBfI[cellI] < 0 )
                cellToBfI[cellI] = bfI;
        }
    }

    //- add cells which shall be refined in a subset
    if( cellSubsetName_ != "" )
    {
        label subsetI = mesh_.cellSubsetIndex(cellSubsetName_);
        if( subsetI >= 0 )
            Warning << "The subset with name " << cellSubsetName_
                    << " already exists. Skipping generation of a new subset"
                    << endl;

        subsetI = mesh_.addCellSubset(cellSubsetName_);

        forAll(nCellsFromCell, cI)
            if( nCellsFromCell[cI] > 1 )
                mesh_.addCellToSubset(subsetI, cI);
    }

    //- check the number of cells which will be generated
    label nNewCells(0);
    forAll(nCellsFromCell, cellI)
        nNewCells += (nCellsFromCell[cellI] - 1);

    # ifdef DEBUGLayer
    forAll(nCellsFromCell, cellI)
    {
        Pout << "\nCell " << cellI << endl;
        Pout << "nCellsFromCell " << nCellsFromCell[cellI] << endl;
        Pout << "Ref type " << refType[cellI] << endl;
    }
    #  endif

    const label totalNumNewCells = returnReduce(nNewCells, sumOp<label>());
    Info << "Number of newly generated cells " << totalNumNewCells << endl;

    //- create mesh modifier
    polyMeshGenModifier meshModifier(mesh_);
    faceListPMG& faces = meshModifier.facesAccess();

    const label numFacesBefore = newFaces_.size();

    //- set the number of cells to the new value
    cellListPMG& cells = meshModifier.cellsAccess();
    label nCells = cells.size();
    cells.setSize(nCells+nNewCells);

    // Exact-volume provenance.
    //
    // A refined parent's local child zero may reuse the original cell label,
    // while later children receive new labels.  Preserve the mapping now so
    // it remains available after newCellsFromCell is discarded and faces
    // are relabelled.
    labelList exactVolumeParent(cells.size(), -1);
    labelList exactVolumeLocalChild(cells.size(), -1);
    labelList exactVolumeRefType(cells.size(), -1);

    // BL_CHILD_SWEEP_AUDIT
    //
    // Diagnostic only.
    //
    // For triangular type-1 prism parents, inspect the exact six points
    // used by every discrete child interval.  Unlike the earlier coarse
    // contact-sweep test, this uses the independently generated
    // newVerticesForSplitEdge_ rows consumed by generateNewCellsPrism().
    //
    // A child is bad if its complete linear-wedge Jacobian is singular
    // somewhere, or if its orientation differs from the full parent sweep.
    label nChildSweepParentsChecked = 0;
    label nChildSweepChildrenChecked = 0;
    label nChildSweepBadChildren = 0;
    label nChildSweepParentUnsafe = 0;
    label nChildSweepSkippedParents = 0;
    label nChildSweepBadPrinted = 0;

    auto childSweepState =
    [&]
    (
        const point& B0,
        const point& B1,
        const point& B2,
        const point& T0,
        const point& T1,
        const point& T2,
        scalar& normalizedMargin
    ) -> label
    {
        const vector H0 = T0 - B0;
        const vector H1 = T1 - B1;
        const vector H2 = T2 - B2;

        const vector A0 = B1 - B0;
        const vector Bv0 = B2 - B0;

        const vector dA = H1 - H0;
        const vector dB = H2 - H0;

        scalar L = scalar(0);

        L = Foam::max(L, mag(B1-B0));
        L = Foam::max(L, mag(B2-B0));
        L = Foam::max(L, mag(B2-B1));

        L = Foam::max(L, mag(T1-T0));
        L = Foam::max(L, mag(T2-T0));
        L = Foam::max(L, mag(T2-T1));

        L = Foam::max(L, mag(H0));
        L = Foam::max(L, mag(H1));
        L = Foam::max(L, mag(H2));

        if( L <= VSMALL )
        {
            normalizedMargin = scalar(0);
            return 0;
        }

        const scalar L3 = L*L*L;

        // Numerical-zero tolerance only, not a quality criterion.
        const scalar tol =
            scalar(1e-12)*L3 + VSMALL;

        scalar globalMin = GREAT;
        scalar globalMax = -GREAT;

        vector H[3];
        H[0] = H0;
        H[1] = H1;
        H[2] = H2;

        for(label k=0; k<3; ++k)
        {
            const scalar c0 =
                (A0 ^ Bv0) & H[k];

            const scalar c1 =
                ((dA ^ Bv0) + (A0 ^ dB)) & H[k];

            const scalar c2 =
                (dA ^ dB) & H[k];

            auto sample =
            [&](const scalar t)
            {
                const scalar q =
                    c0 + c1*t + c2*t*t;

                globalMin =
                    Foam::min(globalMin, q);

                globalMax =
                    Foam::max(globalMax, q);
            };

            sample(scalar(0));
            sample(scalar(1));

            if( mag(c2) > VSMALL )
            {
                const scalar tStar =
                    -c1/(scalar(2)*c2);

                if
                (
                    tStar > scalar(0)
                 && tStar < scalar(1)
                )
                    sample(tStar);
            }
        }

        if( globalMin > tol )
        {
            normalizedMargin = globalMin/L3;
            return 1;
        }

        if( globalMax < -tol )
        {
            normalizedMargin = -globalMax/L3;
            return -1;
        }

        normalizedMargin = scalar(0);
        return 0;
    };


    auto auditExactPrismChildren =
    [&]
    (
        const label parentCellI
    )
    {
        if
        (
            parentCellI < 0
         || parentCellI >= label(cells.size())
        )
        {
            ++nChildSweepSkippedParents;
            return;
        }

        const cell& ac = cells[parentCellI];

        const label startBoundary =
            mesh_.boundaries()[0].patchStart();

        label nLayers = 1;
        label baseFace = -1;

        forAll(ac, fI)
        {
            const label bfI =
                ac[fI] - startBoundary;

            if
            (
                bfI < 0
             || bfI >= label(nLayersAtBndFace_.size())
            )
                continue;

            if( nLayersAtBndFace_[bfI] < 2 )
                continue;

            nLayers =
                nLayersAtBndFace_[bfI];

            baseFace = fI;
        }

        if
        (
            baseFace < 0
         || nLayers < 2
        )
        {
            ++nChildSweepSkippedParents;
            return;
        }

        label otherBaseFace = -1;

        forAll(ac, fI)
        {
            if( fI == baseFace )
                continue;

            if
            (
                facesFromFace_.sizeOfRow(ac[fI]) == 1
            )
                otherBaseFace = fI;
        }

        if( otherBaseFace < 0 )
        {
            ++nChildSweepSkippedParents;
            return;
        }

        const face& bf =
            faces[ac[baseFace]];

        const face& obf =
            faces[ac[otherBaseFace]];

        // This first audit targets the triangular-prism/wedge population,
        // including the known surviving deep-BL failure.  Polygonal type-1
        // parents remain untouched and are counted as skipped.
        if( bf.size() != 3 )
        {
            ++nChildSweepSkippedParents;
            return;
        }

        label seI[3];
        seI[0] = -1;
        seI[1] = -1;
        seI[2] = -1;

        for(label pI=0; pI<3; ++pI)
        {
            const label pointI = bf[pI];

            if
            (
                pointI < 0
             || pointI >= label(splitEdgesAtPoint_.size())
            )
            {
                ++nChildSweepSkippedParents;
                return;
            }

            if
            (
                splitEdgesAtPoint_.sizeOfRow(pointI) == 1
            )
            {
                seI[pI] =
                    splitEdgesAtPoint_(pointI, 0);
            }
            else
            {
                forAllRow
                (
                    splitEdgesAtPoint_,
                    pointI,
                    sepI
                )
                {
                    const label seJ =
                        splitEdgesAtPoint_(pointI, sepI);

                    if
                    (
                        seJ < 0
                     || seJ >= label(splitEdges_.size())
                    )
                        continue;

                    const edge& se =
                        splitEdges_[seJ];

                    if
                    (
                        obf.which(se.end()) >= 0
                     || obf.which(se.start()) >= 0
                    )
                    {
                        seI[pI] = seJ;
                        break;
                    }
                }
            }

            if
            (
                seI[pI] < 0
             || seI[pI] >= label(splitEdges_.size())
             || newVerticesForSplitEdge_.sizeOfRow
                (
                    seI[pI]
                ) <= nLayers
            )
            {
                ++nChildSweepSkippedParents;
                return;
            }
        }

        const pointFieldPMG& auditPoints =
            mesh_.points();

        // Local provenance for the targeted forensic dump.
        // Keep this independent of the later child-sweep provenance block.
        const label badParentDumpBfI =
            ac[baseFace] - startBoundary;

        label badParentDumpPatchI = -1;

        if
        (
            badParentDumpBfI >= 0
         && badParentDumpBfI < label(childSweepFacePatch.size())
        )
            badParentDumpPatchI =
                childSweepFacePatch[badParentDumpBfI];

        word badParentDumpPatchName("?");

        if
        (
            badParentDumpPatchI >= 0
         && badParentDumpPatchI < label(childSweepBoundaries.size())
        )
            badParentDumpPatchName =
                childSweepBoundaries[badParentDumpPatchI].patchName();

        // Temporary forensic dump for the exact Q1 birth parent.
        // Diagnostic only; remove after the mechanism is established.
        if( parentCellI == 1218619 )
        {
            Info
                << "BL_BAD_PARENT_GEOM"
                << " parent=" << parentCellI
                << " bfI=" << badParentDumpBfI
                << " patch=" << badParentDumpPatchName
                << " nLayers=" << nLayers
                << " splitEdges=("
                << seI[0] << " "
                << seI[1] << " "
                << seI[2] << ")"
                << endl;

            Info
                << "BL_BAD_PARENT_CELL_FACES"
                << " parent=" << parentCellI
                << " cellFaces=" << ac
                << endl;

            forAll(ac, srcLocalFI)
            {
                const label srcFaceI =
                    ac[srcLocalFI];

                const label srcBfI =
                    srcFaceI - startBoundary;

                // We only need physical boundary lineage here.
                if
                (
                    srcBfI < 0
                 || srcBfI >= label(childSweepFacePatch.size())
                )
                    continue;

                const label srcPatchI =
                    childSweepFacePatch[srcBfI];

                word srcPatchName("?");

                if
                (
                    srcPatchI >= 0
                 && srcPatchI < label(childSweepBoundaries.size())
                )
                    srcPatchName =
                        childSweepBoundaries[srcPatchI].patchName();

                Info
                    << "BL_BAD_PARENT_SOURCE_FACE"
                    << " parent=" << parentCellI
                    << " localFace=" << srcLocalFI
                    << " faceI=" << srcFaceI
                    << " bfI=" << srcBfI
                    << " patch=" << srcPatchName
                    << " oldPts=" << faces[srcFaceI]
                    << " nDerived="
                    << facesFromFace_.sizeOfRow(srcFaceI)
                    << endl;

                forAllRow
                (
                    facesFromFace_,
                    srcFaceI,
                    srcDerivedI
                )
                {
                    const label nfI =
                        facesFromFace_
                        (
                            srcFaceI,
                            srcDerivedI
                        );

                    if
                    (
                        nfI < 0
                     || nfI >= label(newFaces_.size())
                    )
                        continue;

                    Info
                        << "BL_BAD_PARENT_DERIVED_FACE"
                        << " parent=" << parentCellI
                        << " sourceFaceI=" << srcFaceI
                        << " patch=" << srcPatchName
                        << " derivedLocal=" << srcDerivedI
                        << " newFaceI=" << nfI
                        << " pts=" << newFaces_[nfI]
                        << endl;
                }
            }


            for(label edgeSlot=0; edgeSlot<3; ++edgeSlot)
            {
                const label splitEdgeI = seI[edgeSlot];
                const edge& se = splitEdges_[splitEdgeI];

                Info
                    << "BL_BAD_PARENT_EDGE"
                    << " parent=" << parentCellI
                    << " slot=" << edgeSlot
                    << " seI=" << splitEdgeI
                    << " edge=("
                    << se.start() << " "
                    << se.end() << ")"
                    << " start="
                    << auditPoints[se.start()]
                    << " row1="
                    << auditPoints
                       [
                           newVerticesForSplitEdge_
                           (
                               splitEdgeI,
                               1
                           )
                       ]
                    << " end="
                    << auditPoints[se.end()]
                    << " length="
                    << mag
                       (
                           auditPoints[se.end()]
                         - auditPoints[se.start()]
                       )
                    << endl;
            }
        }

        auto splitPoint =
        [&]
        (
            const label edgeSlot,
            const label rowI
        ) -> const point&
        {
            const label pI =
                newVerticesForSplitEdge_
                (
                    seI[edgeSlot],
                    rowI
                );

            return auditPoints[pI];
        };


        // Full parent sweep using the same three split edges.
        scalar parentMargin = scalar(0);

        const label parentState =
            childSweepState
            (
                splitPoint(0, 0),
                splitPoint(1, 0),
                splitPoint(2, 0),

                splitPoint(0, nLayers),
                splitPoint(1, nLayers),
                splitPoint(2, nLayers),

                parentMargin
            );

        ++nChildSweepParentsChecked;

        const label baseBfI =
            ac[baseFace] - startBoundary;

        label patchI = -1;

        if
        (
            baseBfI >= 0
         && baseBfI < label(childSweepFacePatch.size())
        )
            patchI =
                childSweepFacePatch[baseBfI];

        word patchName("?");

        if
        (
            patchI >= 0
         && patchI < label(childSweepBoundaries.size())
        )
            patchName =
                childSweepBoundaries[patchI].patchName();

        if( parentState == 0 )
        {
            ++nChildSweepParentUnsafe;

            if( nChildSweepBadPrinted < 20 )
            {
                ++nChildSweepBadPrinted;

                Info
                    << "BL_CHILD_SWEEP_PARENT_UNSAFE"
                    << " parentCell=" << parentCellI
                    << " bfI=" << baseBfI
                    << " patch=" << patchName
                    << " nLayers=" << nLayers
                    << " parentMargin="
                    << parentMargin
                    << " splitEdges=("
                    << seI[0] << " "
                    << seI[1] << " "
                    << seI[2] << ")"
                    << endl;
            }

            return;
        }


        auto edgeFraction =
        [&]
        (
            const label edgeSlot,
            const label rowI
        ) -> scalar
        {
            const label splitEdgeI =
                seI[edgeSlot];

            const edge& se =
                splitEdges_[splitEdgeI];

            const vector ev =
                auditPoints[se.end()]
              - auditPoints[se.start()];

            const scalar ev2 = ev & ev;

            if( ev2 <= VSMALL )
                return scalar(0);

            const label pI =
                newVerticesForSplitEdge_
                (
                    splitEdgeI,
                    rowI
                );

            return
                (
                    (
                        auditPoints[pI]
                      - auditPoints[se.start()]
                    )
                  & ev
                ) / ev2;
        };


        for
        (
            label layerI=0;
            layerI<nLayers;
            ++layerI
        )
        {
            scalar childMargin = scalar(0);

            const label childState =
                childSweepState
                (
                    splitPoint(0, layerI),
                    splitPoint(1, layerI),
                    splitPoint(2, layerI),

                    splitPoint(0, layerI+1),
                    splitPoint(1, layerI+1),
                    splitPoint(2, layerI+1),

                    childMargin
                );

            ++nChildSweepChildrenChecked;

            if( childState == parentState )
                continue;

            ++nChildSweepBadChildren;

            const scalar l0 =
                edgeFraction(0, layerI);
            const scalar l1 =
                edgeFraction(1, layerI);
            const scalar l2 =
                edgeFraction(2, layerI);

            const scalar u0 =
                edgeFraction(0, layerI+1);
            const scalar u1 =
                edgeFraction(1, layerI+1);
            const scalar u2 =
                edgeFraction(2, layerI+1);

            const scalar lowerMin =
                Foam::min(l0, Foam::min(l1, l2));

            const scalar lowerMax =
                Foam::max(l0, Foam::max(l1, l2));

            const scalar upperMin =
                Foam::min(u0, Foam::min(u1, u2));

            const scalar upperMax =
                Foam::max(u0, Foam::max(u1, u2));

            if( nChildSweepBadPrinted < 20 )
            {
                ++nChildSweepBadPrinted;

                Info
                    << "BL_CHILD_SWEEP_BAD"
                    << " parentCell=" << parentCellI
                    << " bfI=" << baseBfI
                    << " patch=" << patchName
                    << " nLayers=" << nLayers
                    << " layerInterval="
                    << layerI << "->" << layerI+1
                    << " localChild="
                    << nLayers-1-layerI
                    << " parentState="
                    << parentState
                    << " childState="
                    << childState
                    << " parentMargin="
                    << parentMargin
                    << " childMargin="
                    << childMargin
                    << " splitEdges=("
                    << seI[0] << " "
                    << seI[1] << " "
                    << seI[2] << ")"
                    << " lowerT=("
                    << l0 << " "
                    << l1 << " "
                    << l2 << ")"
                    << " upperT=("
                    << u0 << " "
                    << u1 << " "
                    << u2 << ")"
                    << " lowerSpread="
                    << lowerMax-lowerMin
                    << " upperSpread="
                    << upperMax-upperMin
                    << endl;
            }
        }
    };


    // REFINE_CHILD_CLOSURE_AUDIT
    // Diagnostic only.  Check generated child shells before face-label
    // consolidation/reconstruction can alter their connectivity.
    label nGeneratedChildrenChecked = 0;
    label nGeneratedChildrenBad = 0;
    label nGeneratedBadType1 = 0;
    label nGeneratedBadType2 = 0;
    label nGeneratedBadType3 = 0;
    label nGeneratedBadEdges = 0;
    label nGeneratedMalformedFaces = 0;

    auto auditGeneratedChild =
    [&]
    (
        const auto& childFaces,
        const label parentCellI,
        const label localChildI,
        const label childRefType
    ) -> bool
    {
        ++nGeneratedChildrenChecked;

        std::map<std::pair<label,label>, label> edgeUse;

        bool bad = false;
        label badEdgesThisChild = 0;
        label malformedThisChild = 0;

        forAll(childFaces, fI)
        {
            const auto& f = childFaces[fI];

            if( f.size() < 3 )
            {
                bad = true;
                ++malformedThisChild;
                continue;
            }

            forAll(f, pI)
            {
                const label a = f[pI];
                const label b = f[(pI+1)%f.size()];

                if( a == b )
                    bad = true;

                ++edgeUse
                [
                    std::make_pair
                    (
                        Foam::min(a,b),
                        Foam::max(a,b)
                    )
                ];
            }
        }

        for
        (
            std::map<std::pair<label,label>, label>::const_iterator
                iter=edgeUse.begin();
            iter!=edgeUse.end();
            ++iter
        )
        {
            if( iter->second != 2 )
            {
                bad = true;
                ++badEdgesThisChild;
            }
        }

        if( bad )
        {
            ++nGeneratedChildrenBad;
            nGeneratedBadEdges += badEdgesThisChild;
            nGeneratedMalformedFaces += malformedThisChild;

            if( childRefType == 1 )
                ++nGeneratedBadType1;
            else if( childRefType == 2 )
                ++nGeneratedBadType2;
            else if( childRefType == 3 )
                ++nGeneratedBadType3;

            if( nGeneratedChildrenBad <= 50 )
            {
                Info << "REFINE_CHILD_CLOSURE_BAD"
                     << " parent=" << parentCellI
                     << " childLocal=" << localChildI
                     << " refType=" << childRefType
                     << " nFaces=" << childFaces.size()
                     << " badEdges=" << badEdgesThisChild
                     << " malformedFaces=" << malformedThisChild
                     << endl;
            }
        }

        return !bad;
    };

    // ==================================================================
    // CFMITCH V2.9a PROSPECTIVE FRONT OPTIMIZER
    //
    // Birth-time quality control for type-1 BL refinement.
    //
    // IMPORTANT:
    //
    //   * splitEdges_[seI].start() and .end() remain FIXED.
    //   * Only newly-created interior split-row points move.
    //   * This runs after generateNewFaces(), so prospective child topology
    //     is available, but before generateNewCells() commits those children.
    //
    // The objective differs fundamentally from the post-construction V1C/
    // V1D repair:
    //
    //   - zero negative prospective cells is mandatory;
    //   - fewer bad prospective pyramids wins;
    //   - at equal bad-pyramid count, a better worst pyramid margin wins.
    //
    // A seed transaction is rejected if its moved generated points touch
    // a refType-2/3 parent.  Those junction/intersection topologies require
    // their own prospective evaluator and are deliberately fail-closed here.
    // ==================================================================
    {
        pointFieldPMG& v29Points =
            mesh_.points();

        const labelList& v29Owner =
            mesh_.owner();

        const labelList& v29Neighbour =
            mesh_.neighbour();

        const label v29StartBoundary =
            mesh_.boundaries()[0].patchStart();


        struct V29Score
        {
            label invalid;
            label negative;
            label badPyr;

            scalar negMag;
            scalar minPyr;
            scalar minPositiveVol;

            V29Score()
            :
                invalid(0),
                negative(0),
                badPyr(0),
                negMag(0),
                minPyr(GREAT),
                minPositiveVol(GREAT)
            {}
        };


        // --------------------------------------------------------------
        // OpenFOAM-parity polygon face centre and area, but operating on
        // temporary DynList faces rather than committed mesh faces.
        // --------------------------------------------------------------

        auto v29FaceCentreArea =
        [&]
        (
            const auto& f,
            vector& fCtr,
            vector& fArea
        ) -> bool
        {
            const label nFp =
                f.size();

            if( nFp < 3 )
                return false;

            if( nFp == 3 )
            {
                const point& p0 =
                    v29Points[f[0]];

                const point& p1 =
                    v29Points[f[1]];

                const point& p2 =
                    v29Points[f[2]];

                fArea =
                    scalar(0.5)
                   *((p1-p0)^(p2-p0));

                fCtr =
                    (scalar(1)/scalar(3))
                   *(p0+p1+p2);

                return
                    (
                        !help::isnan(fCtr)
                     && !help::isinf(fCtr)
                     && !help::isnan(fArea)
                     && !help::isinf(fArea)
                    );
            }

            point pAvg(vector::zero);

            for(label fpI=0; fpI<nFp; ++fpI)
                pAvg += v29Points[f[fpI]];

            pAvg /= scalar(nFp);

            vector sumA(vector::zero);

            for(label fpI=0; fpI<nFp; ++fpI)
            {
                const label nextI =
                    (fpI+1)%nFp;

                const point& fp =
                    v29Points[f[fpI]];

                const point& fpNext =
                    v29Points[f[nextI]];

                sumA +=
                    (fpNext-fp)^(pAvg-fp);
            }

            const scalar magSumA =
                mag(sumA);

            if( magSumA <= VSMALL )
                return false;

            const vector sumAHat =
                sumA/(magSumA + VSMALL);

            scalar sumAn =
                scalar(0);

            vector sumAnc(vector::zero);

            for(label fpI=0; fpI<nFp; ++fpI)
            {
                const label nextI =
                    (fpI+1)%nFp;

                const point& fp =
                    v29Points[f[fpI]];

                const point& fpNext =
                    v29Points[f[nextI]];

                const vector a =
                    (fpNext-fp)^(pAvg-fp);

                const vector c =
                    fp + fpNext + pAvg;

                const scalar an =
                    a & sumAHat;

                sumAn += an;
                sumAnc += an*c;
            }

            fArea =
                scalar(0.5)*sumA;

            if( sumAn > VSMALL )
            {
                fCtr =
                    (scalar(1)/scalar(3))
                   *sumAnc/sumAn;
            }
            else
            {
                fCtr = pAvg;
            }

            return
                (
                    !help::isnan(fCtr)
                 && !help::isinf(fCtr)
                 && !help::isnan(fArea)
                 && !help::isinf(fArea)
                );
        };


        // --------------------------------------------------------------
        // Prospective cell quality.
        //
        // childFaces are face POINT LABELS, already oriented relative to
        // the prospective cell by generateNewCellsPrism() or by the type-0
        // builder below.
        //
        // The centre calculation mirrors primitiveMesh:
        //
        //       cEst = average(face centres)
        //       pyr3 = Sf & (Cf-cEst)
        //       C = sum(pyr3*(.75Cf+.25cEst))/sum(pyr3)
        //
        // A positive outward-oriented face pyramid has:
        //
        //       Sf & (Cf-C) > 0
        //
        // so this directly provides the checkMesh-facing pyramid margin.
        // --------------------------------------------------------------

        auto v29EvaluateChild =
        [&]
        (
            const auto& childFaces,
            V29Score& score
        ) -> bool
        {
            const label nCf =
                childFaces.size();

            if( nCf < 4 )
            {
                ++score.invalid;
                return false;
            }

            List<vector> faceCtr(nCf);
            List<vector> faceArea(nCf);

            point cEst(vector::zero);

            for(label cfI=0; cfI<nCf; ++cfI)
            {
                if
                (
                    !v29FaceCentreArea
                    (
                        childFaces[cfI],
                        faceCtr[cfI],
                        faceArea[cfI]
                    )
                )
                {
                    ++score.invalid;
                    return false;
                }

                cEst += faceCtr[cfI];
            }

            cEst /= scalar(nCf);

            // ------------------------------------------------------
            // TOPOLOGICALLY ORIENTED prospective cell evaluation.
            //
            // childFaces MUST already be oriented outward relative to THIS
            // prospective cell.  Orientation is never inferred from the
            // candidate geometry here.
            //
            // For type-1 cells v29BuildOrientedType1Children() below fixes
            // the two cross-layer faces using layer topology while retaining
            // generateNewCellsPrism()'s exact lateral polygons.
            //
            // For type-0 cells v29BuildType0Cell() already reverses every
            // derived source face according to the original owner relation.
            //
            // With cell-local outward faces, this is directly equivalent to
            // primitiveMesh::makeCellCentresAndVols() after applying the
            // owner/neighbour sign for the cell.
            // ------------------------------------------------------

            vector weightedCentre(vector::zero);
            scalar vol3 = scalar(0);

            for(label cfI=0; cfI<nCf; ++cfI)
            {
                const scalar pyr3 =
                    faceArea[cfI]
                  & (
                        faceCtr[cfI]
                       -cEst
                    );

                const vector pc =
                    scalar(0.75)*faceCtr[cfI]
                  + scalar(0.25)*cEst;

                weightedCentre +=
                    pyr3*pc;

                vol3 += pyr3;
            }

            point cellCtr(cEst);

            if( Foam::mag(vol3) > VSMALL )
                cellCtr = weightedCentre/vol3;

            if
            (
                help::isnan(cellCtr)
             || help::isinf(cellCtr)
            )
            {
                ++score.invalid;
                return false;
            }

            const scalar cellVol =
                vol3/scalar(3);

            if( cellVol <= scalar(0) )
            {
                ++score.negative;

                score.negMag +=
                    Foam::mag(cellVol);
            }
            else
            {
                score.minPositiveVol =
                    Foam::min
                    (
                        score.minPositiveVol,
                        cellVol
                    );
            }


            // With outward cell-local orientation:
            //
            //     Sf & (Cf-C) / 3
            //
            // is the positive pyramid margin.  This is equivalent to
            // -pyramidPointFaceRef(...).mag() for an OpenFOAM owner cell.
            for(label cfI=0; cfI<nCf; ++cfI)
            {
                const scalar pyrMargin =
                    (
                        faceArea[cfI]
                      & (
                            faceCtr[cfI]
                           -cellCtr
                        )
                    )
                   /scalar(3);

                score.minPyr =
                    Foam::min
                    (
                        score.minPyr,
                        pyrMargin
                    );

                if( pyrMargin < -SMALL )
                    ++score.badPyr;
            }

            return true;
        };


        // --------------------------------------------------------------
        // Exact split-edge set used by a type-1 parent.
        //
        // This deliberately mirrors generateNewCellsPrism().
        // --------------------------------------------------------------

        auto v29Type1Edges =
        [&]
        (
            const label parentCellI,
            DynList<label, 16>& parentEdges
        ) -> bool
        {
            parentEdges.clear();

            if
            (
                parentCellI < 0
             || parentCellI >= nCells
             || refType[parentCellI] != 1
            )
                return false;

            const cell& ac =
                cells[parentCellI];

            label baseFace = -1;
            label nLayers = 1;

            forAll(ac, fI)
            {
                const label bfI =
                    ac[fI] - v29StartBoundary;

                if
                (
                    bfI < 0
                 || bfI >= label(nLayersAtBndFace_.size())
                )
                    continue;

                if( nLayersAtBndFace_[bfI] < 2 )
                    continue;

                nLayers =
                    nLayersAtBndFace_[bfI];

                baseFace = fI;
            }

            if
            (
                baseFace < 0
             || nLayers < 2
            )
                return false;

            label otherBaseFace = -1;

            forAll(ac, fI)
            {
                if( fI == baseFace )
                    continue;

                if
                (
                    facesFromFace_.sizeOfRow
                    (
                        ac[fI]
                    ) == 1
                )
                    otherBaseFace = fI;
            }

            if( otherBaseFace < 0 )
                return false;

            const face& bf =
                faces[ac[baseFace]];

            const face& obf =
                faces[ac[otherBaseFace]];

            if( bf.size() < 3 )
                return false;

            forAll(bf, pI)
            {
                const label pointI =
                    bf[pI];

                if
                (
                    pointI < 0
                 || pointI >=
                    label(splitEdgesAtPoint_.size())
                )
                    return false;

                label seI = -1;

                if
                (
                    splitEdgesAtPoint_.
                        sizeOfRow(pointI) == 1
                )
                {
                    seI =
                        splitEdgesAtPoint_
                        (
                            pointI,
                            0
                        );
                }
                else
                {
                    forAllRow
                    (
                        splitEdgesAtPoint_,
                        pointI,
                        sepI
                    )
                    {
                        const label seJ =
                            splitEdgesAtPoint_
                            (
                                pointI,
                                sepI
                            );

                        if
                        (
                            seJ < 0
                         || seJ >= label(splitEdges_.size())
                        )
                            continue;

                        const edge& se =
                            splitEdges_[seJ];

                        if
                        (
                            obf.which(se.end()) >= 0
                         || obf.which(se.start()) >= 0
                        )
                        {
                            seI = seJ;
                            break;
                        }
                    }
                }

                if
                (
                    seI < 0
                 || seI >= label(splitEdges_.size())
                 || newVerticesForSplitEdge_.
                    sizeOfRow(seI) < 3
                )
                    return false;

                parentEdges.appendIfNotIn(seI);
            }

            return
                (
                    parentEdges.size()
                 == bf.size()
                );
        };


        // --------------------------------------------------------------
        // Exact type-1 prospective child topology with CELL-LOCAL
        // topological orientation.
        //
        // generateNewCellsPrism() is still used to obtain the exact lateral
        // polygons, including NONQUAD_EDGE_CONFORM / inserted vertices.
        //
        // Its two cross-layer face orderings cannot be trusted directly:
        // the same ordering is deliberately appended to both adjacent
        // children.  Recover those two orientations from the known
        // CORE<->WALL layer topology instead of candidate geometry.
        //
        // Also impose a hard triangle-fan sweep/Jacobian guard using the
        // existing childSweepState().  Any sign change or singular triangle
        // rejects the prospective state before pyramid scoring.
        // --------------------------------------------------------------

        auto v29BuildOrientedType1Children =
        [&]
        (
            const label parentCellI,
            DynList
            <
                DynList
                <
                    DynList<label,8>,
                    10
                >,
                64
            >& children
        ) -> bool
        {
            children.clear();

            DynList<label,16> parentEdges;

            if
            (
                !v29Type1Edges
                (
                    parentCellI,
                    parentEdges
                )
            )
                return false;

            const label nHairs =
                parentEdges.size();

            if( nHairs < 3 )
                return false;

            const label firstEdgeI =
                parentEdges[0];

            if
            (
                firstEdgeI < 0
             || firstEdgeI >= label(splitEdges_.size())
            )
                return false;

            const label nLayers =
                newVerticesForSplitEdge_.
                    sizeOfRow(firstEdgeI) - 1;

            if( nLayers < 2 )
                return false;


            // All hairs of one type-1 parent must use the same number of
            // rows.  Otherwise the ordinary prism constructor itself is
            // not a legal model for this prospective parent.
            forAll(parentEdges, peI)
            {
                const label seI =
                    parentEdges[peI];

                if
                (
                    seI < 0
                 || seI >= label(splitEdges_.size())
                 || newVerticesForSplitEdge_.
                    sizeOfRow(seI) != nLayers + 1
                )
                    return false;
            }


            auto rowPointLabel =
            [&]
            (
                const label hairI,
                const label rowI
            ) -> label
            {
                return
                    newVerticesForSplitEdge_
                    (
                        parentEdges[hairI],
                        rowI
                    );
            };


            auto rowPoint =
            [&]
            (
                const label hairI,
                const label rowI
            ) -> const point&
            {
                return
                    v29Points
                    [
                        rowPointLabel
                        (
                            hairI,
                            rowI
                        )
                    ];
            };


            // ----------------------------------------------------------
            // HARD SWEEP/JACOBIAN GUARD.
            //
            // Fan triangulation is tied to the base-face/hair ordering:
            //
            //     (0,1,2), (0,2,3), (0,3,4), ...
            //
            // For each fan triangle, every discrete child interval must
            // retain the orientation of the complete original sweep.
            // ----------------------------------------------------------

            for
            (
                label fanI=1;
                fanI<nHairs-1;
                ++fanI
            )
            {
                scalar parentMargin =
                    scalar(0);

                const label parentState =
                    childSweepState
                    (
                        rowPoint(0,    0),
                        rowPoint(fanI, 0),
                        rowPoint(fanI+1, 0),

                        rowPoint(0,    nLayers),
                        rowPoint(fanI, nLayers),
                        rowPoint(fanI+1, nLayers),

                        parentMargin
                    );

                if( parentState == 0 )
                    return false;


                for
                (
                    label layerI=0;
                    layerI<nLayers;
                    ++layerI
                )
                {
                    scalar childMargin =
                        scalar(0);

                    const label childState =
                        childSweepState
                        (
                            rowPoint(0,    layerI),
                            rowPoint(fanI, layerI),
                            rowPoint(fanI+1, layerI),

                            rowPoint(0,    layerI+1),
                            rowPoint(fanI, layerI+1),
                            rowPoint(fanI+1, layerI+1),

                            childMargin
                        );

                    if
                    (
                        childState == 0
                     || childState != parentState
                    )
                        return false;
                }
            }


            // Obtain the EXACT prospective polygons cfMesh will consume.
            if
            (
                !generateNewCellsPrism
                (
                    parentCellI,
                    children
                )
            )
                return false;

            if
            (
                children.size() != nLayers
            )
                return false;


            // ----------------------------------------------------------
            // Fix ONLY the two cross-layer face orientations.
            //
            // Type-1 child numbering:
            //
            //     child 0          = core-side child
            //     child N-1        = wall-side child
            //
            // Therefore:
            //
            //     wallRow = N-1-child
            //     coreRow = wallRow+1
            //
            // Base-face/hair ordering is outward at the wall side.
            // Thus:
            //
            //     wall cross-face  = base ordering
            //     core cross-face  = reversed base ordering
            //
            // Lateral polygons are left untouched because
            // generateNewCellsPrism() already reverses them according to
            // the original parent owner relation.
            // ----------------------------------------------------------

            forAll(children, childI)
            {
                const label wallRow =
                    nLayers - 1 - childI;

                const label coreRow =
                    wallRow + 1;

                if
                (
                    wallRow < 0
                 || coreRow > nLayers
                )
                    return false;


                DynList<label,8> wallFace;
                DynList<label,8> coreFace;

                for
                (
                    label hairI=0;
                    hairI<nHairs;
                    ++hairI
                )
                {
                    wallFace.append
                    (
                        rowPointLabel
                        (
                            hairI,
                            wallRow
                        )
                    );
                }

                for
                (
                    label hairI=nHairs-1;
                    hairI>=0;
                    --hairI
                )
                {
                    coreFace.append
                    (
                        rowPointLabel
                        (
                            hairI,
                            coreRow
                        )
                    );
                }


                label wallFaceI = -1;
                label coreFaceI = -1;

                forAll(children[childI], cfI)
                {
                    const DynList<label,8>& f =
                        children[childI][cfI];

                    if
                    (
                        help::areFacesEqual
                        (
                            f,
                            wallFace
                        )
                    )
                    {
                        if( wallFaceI >= 0 )
                            return false;

                        wallFaceI = cfI;
                    }

                    if
                    (
                        help::areFacesEqual
                        (
                            f,
                            coreFace
                        )
                    )
                    {
                        if( coreFaceI >= 0 )
                            return false;

                        coreFaceI = cfI;
                    }
                }


                if
                (
                    wallFaceI < 0
                 || coreFaceI < 0
                 || wallFaceI == coreFaceI
                )
                    return false;


                children[childI][wallFaceI] =
                    wallFace;

                children[childI][coreFaceI] =
                    coreFace;
            }


            return true;
        };


        // --------------------------------------------------------------
        // Prospective type-0 cell.
        //
        // facesFromFace_ has already been generated.  Preserve the
        // original source-face orientation relative to this cell.
        // --------------------------------------------------------------

        auto v29BuildType0Cell =
        [&]
        (
            const label cellI,
            DynList<DynList<label,8>, 10>& childFaces
        ) -> bool
        {
            childFaces.clear();

            if
            (
                cellI < 0
             || cellI >= nCells
             || refType[cellI] != 0
            )
                return false;

            const cell& c =
                cells[cellI];

            forAll(c, cfI)
            {
                const label sourceFaceI =
                    c[cfI];

                if
                (
                    sourceFaceI < 0
                 || sourceFaceI >=
                    label(facesFromFace_.size())
                 || sourceFaceI >=
                    label(v29Owner.size())
                )
                    return false;

                if
                (
                    facesFromFace_.
                        sizeOfRow(sourceFaceI) == 0
                )
                    return false;

                forAllRow
                (
                    facesFromFace_,
                    sourceFaceI,
                    dfI
                )
                {
                    const label newFaceI =
                        facesFromFace_
                        (
                            sourceFaceI,
                            dfI
                        );

                    if
                    (
                        newFaceI < 0
                     || newFaceI >=
                        label(newFaces_.size())
                    )
                        return false;

                    DynList<label,8> f;

                    f =
                        newFaces_[newFaceI];

                    if
                    (
                        v29Owner[sourceFaceI]
                     != cellI
                    )
                    {
                        f =
                            help::reverseFace(f);
                    }

                    childFaces.append(f);
                }
            }

            return
                childFaces.size() >= 4;
        };


        // --------------------------------------------------------------
        // Prospective quality of one ORIGINAL coarse parent.
        //
        // refType 0 => one replacement cell
        // refType 1 => exact future BL children from generateNewCellsPrism
        //
        // refType 2/3 intentionally fail closed.
        // --------------------------------------------------------------

        auto v29EvaluateParent =
        [&]
        (
            const label parentCellI,
            V29Score& score
        ) -> bool
        {
            if
            (
                parentCellI < 0
             || parentCellI >= nCells
            )
            {
                ++score.invalid;
                return false;
            }

            if( refType[parentCellI] == 0 )
            {
                DynList
                <
                    DynList<label,8>,
                    10
                > childFaces;

                if
                (
                    !v29BuildType0Cell
                    (
                        parentCellI,
                        childFaces
                    )
                )
                {
                    ++score.invalid;
                    return false;
                }

                return
                    v29EvaluateChild
                    (
                        childFaces,
                        score
                    );
            }

            if( refType[parentCellI] == 1 )
            {
                DynList
                <
                    DynList
                    <
                        DynList<label,8>,
                        10
                    >,
                    64
                > children;

                if
                (
                    !v29BuildOrientedType1Children
                    (
                        parentCellI,
                        children
                    )
                )
                {
                    ++score.invalid;
                    return false;
                }

                forAll(children, childI)
                {
                    if
                    (
                        !v29EvaluateChild
                        (
                            children[childI],
                            score
                        )
                    )
                        return false;
                }

                return true;
            }

            ++score.invalid;
            return false;
        };


        // --------------------------------------------------------------
        // Wall-adjacent prospective child seed predicate.
        //
        // Type-1 local child numbering is CORE -> WALL, therefore the last
        // prospective child is the wall-adjacent child.
        // --------------------------------------------------------------

        auto v29WallChildScore =
        [&]
        (
            const label parentCellI,
            V29Score& score
        ) -> bool
        {
            if
            (
                parentCellI < 0
             || parentCellI >= nCells
             || refType[parentCellI] != 1
            )
                return false;

            DynList
            <
                DynList
                <
                    DynList<label,8>,
                    10
                >,
                64
            > children;

            if
            (
                !v29BuildOrientedType1Children
                (
                    parentCellI,
                    children
                )
            )
                return false;

            return
                v29EvaluateChild
                (
                    children
                    [
                        children.size()-1
                    ],
                    score
                );
        };


        // --------------------------------------------------------------
        // newFace -> source original face.
        // --------------------------------------------------------------

        labelList v29NewFaceSource
        (
            newFaces_.size(),
            -1
        );

        for
        (
            label sourceFaceI=0;
            sourceFaceI<label(facesFromFace_.size());
            ++sourceFaceI
        )
        {
            forAllRow
            (
                facesFromFace_,
                sourceFaceI,
                dfI
            )
            {
                const label newFaceI =
                    facesFromFace_
                    (
                        sourceFaceI,
                        dfI
                    );

                if
                (
                    newFaceI >= 0
                 && newFaceI <
                    label(v29NewFaceSource.size())
                 && v29NewFaceSource[newFaceI] < 0
                )
                {
                    v29NewFaceSource[newFaceI] =
                        sourceFaceI;
                }
            }
        }


        // Point -> generated face reverse addressing.
        VRWGraph v29PointNewFaces;

        v29PointNewFaces.reverseAddressing
        (
            newFaces_
        );


        auto v29AffectedParents =
        [&]
        (
            const DynList<label,16>& seedEdges,
            std::set<label>& affected
        ) -> bool
        {
            affected.clear();

            forAll(seedEdges, eeI)
            {
                const label seI =
                    seedEdges[eeI];

                if
                (
                    seI < 0
                 || seI >= label(splitEdges_.size())
                )
                    return false;

                const label rowSize =
                    newVerticesForSplitEdge_.
                        sizeOfRow(seI);

                if( rowSize < 3 )
                    return false;

                for
                (
                    label rowI=1;
                    rowI<rowSize-1;
                    ++rowI
                )
                {
                    const label pointI =
                        newVerticesForSplitEdge_
                        (
                            seI,
                            rowI
                        );

                    if
                    (
                        pointI < 0
                     || pointI >=
                        label(v29PointNewFaces.size())
                    )
                        return false;

                    forAllRow
                    (
                        v29PointNewFaces,
                        pointI,
                        pfI
                    )
                    {
                        const label newFaceI =
                            v29PointNewFaces
                            (
                                pointI,
                                pfI
                            );

                        if
                        (
                            newFaceI < 0
                         || newFaceI >=
                            label(v29NewFaceSource.size())
                        )
                            continue;

                        const label sourceFaceI =
                            v29NewFaceSource[newFaceI];

                        if
                        (
                            sourceFaceI < 0
                         || sourceFaceI >=
                            label(v29Owner.size())
                        )
                            continue;

                        const label own =
                            v29Owner[sourceFaceI];

                        if
                        (
                            own >= 0
                         && own < nCells
                        )
                            affected.insert(own);

                        if
                        (
                            sourceFaceI <
                            label(v29Neighbour.size())
                        )
                        {
                            const label nei =
                                v29Neighbour[sourceFaceI];

                            if
                            (
                                nei >= 0
                             && nei < nCells
                            )
                                affected.insert(nei);
                        }
                    }
                }
            }

            return !affected.empty();
        };


        auto v29EvaluateSet =
        [&]
        (
            const std::set<label>& affected,
            V29Score& score
        ) -> bool
        {
            for
            (
                std::set<label>::const_iterator
                    cIt=affected.begin();
                cIt!=affected.end();
                ++cIt
            )
            {
                const label cellI =
                    *cIt;

                if
                (
                    refType[cellI] != 0
                 && refType[cellI] != 1
                )
                {
                    ++score.invalid;
                    return false;
                }

                if
                (
                    !v29EvaluateParent
                    (
                        cellI,
                        score
                    )
                )
                    return false;
            }

            return true;
        };


        // --------------------------------------------------------------
        // Quality comparison.
        //
        // Negative prospective cells are never introduced.
        //
        // A lower bad-pyramid count wins, provided:
        //
        //   * minimum positive cell volume retains >= 50% baseline;
        //   * worst pyramid margin is not made >2x worse.
        //
        // At equal bad-pyramid count the worst margin must improve.
        // --------------------------------------------------------------

        auto v29CandidateBetter =
        [&]
        (
            const V29Score& trial,
            const V29Score& base
        ) -> bool
        {
            // v2.9a is intentionally strict:
            //
            //   * never introduce a prospective negative cell;
            //   * retain at least half the baseline minimum raw volume;
            //   * MUST remove at least one prospective bad pyramid.
            //
            // Equal bad-pyramid count is NOT an improvement in v2.9a.
            // This prevents large row deformations being accepted for
            // numerically insignificant margin changes.

            if
            (
                trial.invalid != 0
             || base.invalid != 0
            )
                return false;

            if
            (
                base.negative != 0
             || trial.negative != 0
            )
                return false;

            if
            (
                base.minPositiveVol < GREAT/scalar(2)
             && trial.minPositiveVol <
                scalar(0.50)*base.minPositiveVol
            )
                return false;

            if( trial.badPyr >= base.badPyr )
                return false;

            // Permit a lower total bad-pyramid count to win, but prevent
            // trading it for an extreme worsening of the worst survivor.
            if
            (
                base.minPyr < scalar(0)
             && trial.minPyr <
                scalar(2.0)*base.minPyr
            )
                return false;

            return true;
        };


        label v29Seeds = 0;
        label v29Attempted = 0;
        label v29Accepted = 0;
        label v29SkippedMixed = 0;
        label v29SkippedInvalid = 0;
        label v29Trials = 0;

        label v29Reported = 0;


        const scalar v29Amplitudes[] =
        {
            scalar(-0.01), scalar( 0.01),
            scalar(-0.02), scalar( 0.02),
            scalar(-0.03), scalar( 0.03),
            scalar(-0.05), scalar( 0.05),
            scalar(-0.08), scalar( 0.08),
            scalar(-0.12), scalar( 0.12),
            scalar(-0.16), scalar( 0.16),
            scalar(-0.20), scalar( 0.20)
        };

        const label nV29Amplitudes =
            sizeof(v29Amplitudes)
           /sizeof(v29Amplitudes[0]);


        for
        (
            label parentCellI=0;
            parentCellI<nCells;
            ++parentCellI
        )
        {
            if( refType[parentCellI] != 1 )
                continue;


            // Seed only on a positive-volume prospective WALL child whose
            // baseline geometry contains a bad pyramid.
            V29Score wallScore;

            if
            (
                !v29WallChildScore
                (
                    parentCellI,
                    wallScore
                )
            )
                continue;

            if
            (
                wallScore.invalid != 0
             || wallScore.negative != 0
             || wallScore.badPyr == 0
            )
                continue;

            ++v29Seeds;


            DynList<label,16> seedEdges;

            if
            (
                !v29Type1Edges
                (
                    parentCellI,
                    seedEdges
                )
            )
            {
                ++v29SkippedInvalid;
                continue;
            }


            std::set<label> affected;

            if
            (
                !v29AffectedParents
                (
                    seedEdges,
                    affected
                )
            )
            {
                ++v29SkippedInvalid;
                continue;
            }

            affected.insert(parentCellI);


            bool mixedUnsupported = false;

            for
            (
                std::set<label>::const_iterator
                    cIt=affected.begin();
                cIt!=affected.end();
                ++cIt
            )
            {
                if
                (
                    refType[*cIt] != 0
                 && refType[*cIt] != 1
                )
                {
                    mixedUnsupported = true;
                    break;
                }
            }

            if( mixedUnsupported )
            {
                ++v29SkippedMixed;
                continue;
            }


            V29Score baseline;

            if
            (
                !v29EvaluateSet
                (
                    affected,
                    baseline
                )
             || baseline.invalid != 0
             || baseline.negative != 0
            )
            {
                ++v29SkippedInvalid;
                continue;
            }


            // Snapshot every unique interior split-row point.
            std::map<label,point> originalPositions;

            bool snapshotValid = true;

            forAll(seedEdges, eeI)
            {
                const label seI =
                    seedEdges[eeI];

                const edge& se =
                    splitEdges_[seI];

                const label rowSize =
                    newVerticesForSplitEdge_.
                        sizeOfRow(seI);

                if( rowSize < 3 )
                {
                    snapshotValid = false;
                    break;
                }

                for
                (
                    label rowI=1;
                    rowI<rowSize-1;
                    ++rowI
                )
                {
                    const label pointI =
                        newVerticesForSplitEdge_
                        (
                            seI,
                            rowI
                        );

                    if
                    (
                        pointI == se.start()
                     || pointI == se.end()
                    )
                    {
                        snapshotValid = false;
                        break;
                    }

                    std::map<label,point>::const_iterator
                        oldIt =
                            originalPositions.find(pointI);

                    if
                    (
                        oldIt != originalPositions.end()
                    )
                    {
                        // A generated point participating in more than one
                        // seed hair is not a legal v2.9a scalar variable.
                        snapshotValid = false;
                        break;
                    }

                    originalPositions.insert
                    (
                        std::make_pair
                        (
                            pointI,
                            v29Points[pointI]
                        )
                    );
                }

                if( !snapshotValid )
                    break;
            }

            if
            (
                !snapshotValid
             || originalPositions.empty()
            )
            {
                ++v29SkippedInvalid;
                continue;
            }


            auto restoreOriginal =
            [&]()
            {
                for
                (
                    std::map<label,point>::const_iterator
                        pIt=originalPositions.begin();
                    pIt!=originalPositions.end();
                    ++pIt
                )
                {
                    v29Points[pIt->first] =
                        pIt->second;
                }
            };


            auto applyAmplitude =
            [&]
            (
                const scalar amplitude
            ) -> bool
            {
                restoreOriginal();

                forAll(seedEdges, eeI)
                {
                    const label seI =
                        seedEdges[eeI];

                    const edge& se =
                        splitEdges_[seI];

                    const vector edgeVec =
                        v29Points[se.end()]
                       -v29Points[se.start()];

                    const scalar edgeMagSqr =
                        edgeVec & edgeVec;

                    if( edgeMagSqr <= VSMALL )
                    {
                        restoreOriginal();
                        return false;
                    }

                    const label rowSize =
                        newVerticesForSplitEdge_.
                            sizeOfRow(seI);

                    const label firstPointI =
                        newVerticesForSplitEdge_
                        (
                            seI,
                            1
                        );

                    std::map<label,point>::const_iterator
                        firstIt =
                            originalPositions.find
                            (
                                firstPointI
                            );

                    if
                    (
                        firstIt ==
                        originalPositions.end()
                    )
                    {
                        restoreOriginal();
                        return false;
                    }

                    const scalar firstT =
                        (
                            (
                                firstIt->second
                               -v29Points[se.start()]
                            )
                          & edgeVec
                        )
                       /(edgeMagSqr + VSMALL);

                    if
                    (
                        firstT <= scalar(0)
                     || firstT >= scalar(1)
                    )
                    {
                        restoreOriginal();
                        return false;
                    }

                    scalar prevOriginalT =
                        scalar(0);

                    scalar prevWarpedT =
                        scalar(0);

                    for
                    (
                        label rowI=1;
                        rowI<rowSize-1;
                        ++rowI
                    )
                    {
                        const label pointI =
                            newVerticesForSplitEdge_
                            (
                                seI,
                                rowI
                            );

                        std::map<label,point>::const_iterator
                            pIt =
                                originalPositions.find
                                (
                                    pointI
                                );

                        if
                        (
                            pIt ==
                            originalPositions.end()
                        )
                        {
                            restoreOriginal();
                            return false;
                        }

                        const scalar t =
                            (
                                (
                                    pIt->second
                                   -v29Points[se.start()]
                                )
                              & edgeVec
                            )
                           /(edgeMagSqr + VSMALL);

                        if
                        (
                            !(t > prevOriginalT)
                         || !(t < scalar(1))
                        )
                        {
                            restoreOriginal();
                            return false;
                        }

                        const scalar decayBase =
                            Foam::max
                            (
                                scalar(0),
                                (
                                    scalar(1)-t
                                )
                               /(
                                    scalar(1)
                                   -firstT
                                   +VSMALL
                                )
                            );

                        const scalar deltaT =
                            amplitude
                           *firstT
                           *Foam::pow
                            (
                                decayBase,
                                scalar(4)
                            );

                        const scalar warpedT =
                            t + deltaT;

                        if
                        (
                            !(warpedT > prevWarpedT)
                         || !(warpedT < scalar(1))
                        )
                        {
                            restoreOriginal();
                            return false;
                        }

                        v29Points[pointI] =
                            v29Points[se.start()]
                          + warpedT*edgeVec;

                        prevOriginalT = t;
                        prevWarpedT = warpedT;
                    }
                }

                return true;
            };


            ++v29Attempted;

            bool foundBetter = false;

            V29Score bestScore =
                baseline;

            scalar bestAmplitude =
                scalar(0);

            std::map<label,point> bestPositions;


            for
            (
                label ampI=0;
                ampI<nV29Amplitudes;
                ++ampI
            )
            {
                const scalar amplitude =
                    v29Amplitudes[ampI];

                if
                (
                    !applyAmplitude
                    (
                        amplitude
                    )
                )
                    continue;

                ++v29Trials;

                V29Score trial;

                if
                (
                    !v29EvaluateSet
                    (
                        affected,
                        trial
                    )
                )
                    continue;

                if
                (
                    !v29CandidateBetter
                    (
                        trial,
                        bestScore
                    )
                )
                    continue;

                foundBetter = true;

                bestScore =
                    trial;

                bestAmplitude =
                    amplitude;

                bestPositions.clear();

                for
                (
                    std::map<label,point>::const_iterator
                        pIt=originalPositions.begin();
                    pIt!=originalPositions.end();
                    ++pIt
                )
                {
                    bestPositions.insert
                    (
                        std::make_pair
                        (
                            pIt->first,
                            v29Points[pIt->first]
                        )
                    );
                }
            }


            restoreOriginal();


            if( foundBetter )
            {
                for
                (
                    std::map<label,point>::const_iterator
                        pIt=bestPositions.begin();
                    pIt!=bestPositions.end();
                    ++pIt
                )
                {
                    v29Points[pIt->first] =
                        pIt->second;
                }

                ++v29Accepted;

                if( v29Reported < 50 )
                {
                    ++v29Reported;

                    Info
                        << "CFMITCH V2.9a.1 FRONT ACCEPT"
                        << " parent=" << parentCellI
                        << " hairs=" << seedEdges.size()
                        << " affectedParents="
                        << affected.size()
                        << " amplitude="
                        << bestAmplitude
                        << " badPyr="
                        << baseline.badPyr
                        << "->"
                        << bestScore.badPyr
                        << " minPyr="
                        << baseline.minPyr
                        << "->"
                        << bestScore.minPyr
                        << " minPositiveVol="
                        << baseline.minPositiveVol
                        << "->"
                        << bestScore.minPositiveVol
                        << endl;
                }
            }
        }


        Info
            << "CFMITCH V2.9a.1 PROSPECTIVE FRONT SUMMARY:"
            << " seeds=" << v29Seeds
            << " attempted=" << v29Attempted
            << " accepted=" << v29Accepted
            << " trials=" << v29Trials
            << " skippedMixed=" << v29SkippedMixed
            << " skippedInvalid=" << v29SkippedInvalid
            << endl;
    }


    //- provenance map: newCellI -> bfI that generated it (-1 if not a BL cell)
    cellToBaseBndFace_.setSize(nCells+nNewCells, -1);
    forAll(cellToBfI, cI)
        if( cellToBfI[cI] >= 0 )
            cellToBaseBndFace_[cI] = cellToBfI[cI];

    //- start creating new cells
    //- store the information which new cells were generated from
    //- an existing cell
    VRWGraph newCellsFromCell(refType.size());

    VRWGraph pointNewFaces;
    pointNewFaces.reverseAddressing(newFaces_);

    forAll(nCellsFromCell, cellI)
    {
        if( refType[cellI] == 0 )
        {
            //- this cell is not refined
            //- update face labels
            newCellsFromCell.append(cellI, cellI);

            cell& c = cells[cellI];

            //- copy the new faces of this cell
            DynList<label, 64> newC;
            forAll(c, fI)
            {
                forAllRow(facesFromFace_, c[fI], cfI)
                    newC.append(facesFromFace_(c[fI], cfI));
            }

            // TYPE0_PROVENANCE_AUDIT
            //
            // Target the first simple bad cell and one larger bad cell
            // previously identified by REFINE_PRE_RELABEL_CLOSURE.
            // At this point:
            //
            //   c     = original cell face labels
            //   newC  = replacement newFaces_ labels
            //
            // Nothing has yet been written back to c.
            if( cellI == 951557 || cellI == 951730 )
            {
                std::map<std::pair<label,label>, label> edgeUse;

                forAll(newC, nfLocalI)
                {
                    const label nfI = newC[nfLocalI];

                    if
                    (
                        nfI < 0
                     || nfI >= label(newFaces_.size())
                    )
                    {
                        Info << "TYPE0_PROVENANCE_BAD_REF"
                             << " cell=" << cellI
                             << " newFace=" << nfI
                             << endl;
                        continue;
                    }

                    const label nPts = newFaces_.sizeOfRow(nfI);

                    for(label pI=0; pI<nPts; ++pI)
                    {
                        const label a = newFaces_(nfI, pI);
                        const label b =
                            newFaces_(nfI, (pI+1)%nPts);

                        ++edgeUse
                        [
                            std::make_pair
                            (
                                Foam::min(a,b),
                                Foam::max(a,b)
                            )
                        ];
                    }
                }

                label nBadEdges = 0;
                for
                (
                    std::map<std::pair<label,label>, label>::const_iterator
                        iter=edgeUse.begin();
                    iter!=edgeUse.end();
                    ++iter
                )
                {
                    if( iter->second != 2 )
                        ++nBadEdges;
                }

                Info << "TYPE0_PROVENANCE_BEGIN"
                     << " cell=" << cellI
                     << " originalFaces=" << c.size()
                     << " replacementFaces=" << newC.size()
                     << " badEdges=" << nBadEdges
                     << endl;

                // Print every original face and every replacement face
                // originating from it.  These cells are tiny, so this is
                // deliberately exhaustive.
                forAll(c, oldLocalI)
                {
                    const label oldFaceI = c[oldLocalI];

                    Info << "TYPE0_PROVENANCE_OLD"
                         << " cell=" << cellI
                         << " local=" << oldLocalI
                         << " oldFace=" << oldFaceI;

                    if
                    (
                        oldFaceI >= 0
                     && oldFaceI < label(faces.size())
                    )
                    {
                        Info << " oldPts=" << faces[oldFaceI];
                    }
                    else
                    {
                        Info << " oldPts=OUT_OF_RANGE";
                    }

                    Info << " replacements=(";

                    if
                    (
                        oldFaceI >= 0
                     && oldFaceI < label(facesFromFace_.size())
                    )
                    {
                        forAllRow
                        (
                            facesFromFace_,
                            oldFaceI,
                            repI
                        )
                        {
                            Info << ' '
                                 << facesFromFace_(oldFaceI, repI);
                        }
                    }

                    Info << " )" << endl;

                    if
                    (
                        oldFaceI < 0
                     || oldFaceI >= label(facesFromFace_.size())
                    )
                        continue;

                    forAllRow(facesFromFace_, oldFaceI, repI)
                    {
                        const label nfI =
                            facesFromFace_(oldFaceI, repI);

                        Info << "TYPE0_PROVENANCE_REP"
                             << " cell=" << cellI
                             << " oldFace=" << oldFaceI
                             << " repLocal=" << repI
                             << " newFace=" << nfI;

                        if
                        (
                            nfI >= 0
                         && nfI < label(newFaces_.size())
                        )
                        {
                            Info << " newPts="
                                 << newFaces_[nfI];
                        }
                        else
                        {
                            Info << " newPts=OUT_OF_RANGE";
                        }

                        Info << endl;
                    }
                }

                // Print all unmatched edges and identify exactly which
                // original/replacement face contributes each edge.
                for
                (
                    std::map<std::pair<label,label>, label>::const_iterator
                        iter=edgeUse.begin();
                    iter!=edgeUse.end();
                    ++iter
                )
                {
                    if( iter->second == 2 )
                        continue;

                    const label ea = iter->first.first;
                    const label eb = iter->first.second;

                    Info << "TYPE0_PROVENANCE_BAD_EDGE"
                         << " cell=" << cellI
                         << " edge=(" << ea << ' ' << eb << ')'
                         << " use=" << iter->second
                         << endl;

                    forAll(c, oldLocalI)
                    {
                        const label oldFaceI = c[oldLocalI];

                        if
                        (
                            oldFaceI < 0
                         || oldFaceI >= label(facesFromFace_.size())
                        )
                            continue;

                        forAllRow
                        (
                            facesFromFace_,
                            oldFaceI,
                            repI
                        )
                        {
                            const label nfI =
                                facesFromFace_(oldFaceI, repI);

                            if
                            (
                                nfI < 0
                             || nfI >= label(newFaces_.size())
                            )
                                continue;

                            const label nPts =
                                newFaces_.sizeOfRow(nfI);

                            bool containsBadEdge = false;

                            for(label pI=0; pI<nPts; ++pI)
                            {
                                const label a =
                                    newFaces_(nfI, pI);
                                const label b =
                                    newFaces_
                                    (
                                        nfI,
                                        (pI+1)%nPts
                                    );

                                if
                                (
                                    Foam::min(a,b) == ea
                                 && Foam::max(a,b) == eb
                                )
                                {
                                    containsBadEdge = true;
                                    break;
                                }
                            }

                            if( containsBadEdge )
                            {
                                Info
                                    << "TYPE0_PROVENANCE_EDGE_SOURCE"
                                    << " cell=" << cellI
                                    << " edge=("
                                    << ea << ' ' << eb << ')'
                                    << " oldFace=" << oldFaceI
                                    << " replacementFace=" << nfI
                                    << " oldPts=";

                                if
                                (
                                    oldFaceI >= 0
                                 && oldFaceI < label(faces.size())
                                )
                                    Info << faces[oldFaceI];
                                else
                                    Info << "OUT_OF_RANGE";

                                Info << " newPts="
                                     << newFaces_[nfI]
                                     << endl;
                            }
                        }
                    }
                }

                Info << "TYPE0_PROVENANCE_END"
                     << " cell=" << cellI
                     << endl;
            }

            //- update the cell
            c.setSize(newC.size());
            forAll(c, fI)
                c[fI] = newC[fI];
        }
        else if( refType[cellI] == 1 )
        {
            //- generate new cells from this prism refined in one direction
            DynList<DynList<DynList<label, 8>, 10>, 64> cellsFromCell;

            // Report-only audit of the exact discrete split geometry.
            auditExactPrismChildren(cellI);

            if
            (
                !generateNewCellsPrism
                (
                    cellI,
                    cellsFromCell
                )
            )
            {
                refinementValid_ = false;

                WarningIn
                (
                    "void refineBoundaryLayers::generateNewCells()"
                )
                    << "CFMitch V3.5 rejected structurally invalid "
                    << "type-1 prism parent " << cellI
                    << " -- aborting this refinement transaction"
                    << endl;

                return;
            }

            forAll(cellsFromCell, cI)
            {
                const DynList<DynList<label, 8>, 10>& nc = cellsFromCell[cI];

                if( cellI == 1218619 && cI == 14 )
                {
                    Info
                        << "BL_BAD_CHILD_SHELL"
                        << " parent=" << cellI
                        << " localChild=" << cI
                        << " nFaces=" << nc.size()
                        << endl;

                    forAll(nc, badChildFI)
                    {
                        Info
                            << "BL_BAD_CHILD_FACE"
                            << " parent=" << cellI
                            << " localChild=" << cI
                            << " localFace=" << badChildFI
                            << " pts=" << nc[badChildFI]
                            << endl;
                    }
                }

                auditGeneratedChild(nc, cellI, cI, refType[cellI]);

                const label newCellI = cI==0?cellI:nCells++;
                cellToBaseBndFace_[newCellI] = cellToBfI[cellI];

                exactVolumeParent[newCellI] = cellI;
                exactVolumeLocalChild[newCellI] = cI;
                exactVolumeRefType[newCellI] = 1;

                newCellsFromCell.append(cellI, newCellI);

                cell& c = cells[newCellI];
                c.setSize(nc.size());

                //- find face labels for this cell
                forAll(nc, fI)
                {
                    const DynList<label, 8>& nf = nc[fI];

                    label faceLabel(-1);
                    forAllRow(pointNewFaces, nf[0], pfI)
                    {
                        const label nfI = pointNewFaces(nf[0], pfI);

                        if( help::areFacesEqual(nf, newFaces_[nfI]) )
                        {
                            c[fI] = nfI;
                            faceLabel = nfI;
                            break;
                        }
                    }

                    if( faceLabel < 0 )
                    {
                        forAll(nf, pI)
                            pointNewFaces.append(nf[pI], newFaces_.size());
                        c[fI] = newFaces_.size();
                        newFaces_.appendList(nf);
                    }
                }
            }
        }
        else if( refType[cellI] == 2 )
        {
            // ---------------------------------------------------------
            // CFMitch V3.6 -- type-2 edge-hex structural preflight.
            //
            // refType==2 only proves that two boundary faces request
            // refinement.  It does NOT prove that the parent is the
            // six-faced topological hex required by refineEdgeHexCell.
            //
            // Validate that assumption before the constructor indexes
            // its FixedList<...,6> directional storage.
            // ---------------------------------------------------------

            // CFMitch V3.8: retain the original parent while generated
            // child zero reuses and overwrites cells[cellI].
            const cell type2Parent(cells[cellI]);

            FixedList<label, 2> type2BfI(-1);
            FixedList<label, 2> type2BoundaryLocalFace(-1);

            label type2ActiveBoundaryFaces = 0;
            bool type2ParentValid = true;
            word type2FailureReason("none");

            if( type2Parent.size() != 6 )
            {
                type2ParentValid = false;
                type2FailureReason = "parentNotSixFaced";
            }

            const label type2StartBoundary =
                childSweepBoundaries.size()
              ? childSweepBoundaries[0].patchStart()
              : -1;

            std::map<std::pair<label,label>, label>
                type2ParentEdgeUse;

            if( type2ParentValid )
            {
                forAll(type2Parent, parentLocalFI)
                {
                    const label parentFaceI =
                        type2Parent[parentLocalFI];

                    if
                    (
                        parentFaceI < 0
                     || parentFaceI >= label(faces.size())
                    )
                    {
                        type2ParentValid = false;
                        type2FailureReason =
                            "parentFaceOutOfRange";
                        break;
                    }

                    const face& parentFace =
                        faces[parentFaceI];

                    if( parentFace.size() < 3 )
                    {
                        type2ParentValid = false;
                        type2FailureReason =
                            "parentFaceTooSmall";
                        break;
                    }

                    forAll(parentFace, parentPointI)
                    {
                        const label a =
                            parentFace[parentPointI];

                        const label b =
                            parentFace
                            [
                                (parentPointI+1)
                              % parentFace.size()
                            ];

                        if
                        (
                            a < 0
                         || b < 0
                         || a >= label(mesh_.points().size())
                         || b >= label(mesh_.points().size())
                         || a == b
                        )
                        {
                            type2ParentValid = false;
                            type2FailureReason =
                                "parentFaceBadEdge";
                            break;
                        }

                        ++type2ParentEdgeUse
                        [
                            std::make_pair
                            (
                                Foam::min(a,b),
                                Foam::max(a,b)
                            )
                        ];
                    }

                    if( !type2ParentValid )
                        break;

                    const label bfI =
                        parentFaceI - type2StartBoundary;

                    if
                    (
                        bfI >= 0
                     && bfI < label(nLayersAtBndFace_.size())
                     && nLayersAtBndFace_[bfI] > 1
                    )
                    {
                        if( type2ActiveBoundaryFaces < 2 )
                        {
                            type2BfI
                            [
                                type2ActiveBoundaryFaces
                            ] = bfI;

                            type2BoundaryLocalFace
                            [
                                type2ActiveBoundaryFaces
                            ] = parentLocalFI;
                        }

                        ++type2ActiveBoundaryFaces;
                    }
                }
            }

            if
            (
                type2ParentValid
             && type2ActiveBoundaryFaces != 2
            )
            {
                type2ParentValid = false;
                type2FailureReason =
                    "activeBoundaryFaceCount";
            }

            if( type2ParentValid )
            {
                for
                (
                    std::map
                    <
                        std::pair<label,label>,
                        label
                    >::const_iterator edgeIt =
                        type2ParentEdgeUse.begin();

                    edgeIt != type2ParentEdgeUse.end();
                    ++edgeIt
                )
                {
                    if( edgeIt->second != 2 )
                    {
                        type2ParentValid = false;
                        type2FailureReason =
                            "parentShellOpen";
                        break;
                    }
                }
            }

            label type2CommonEdgeCount = 0;

            if( type2ParentValid )
            {
                const face& activeFace0 =
                    faces
                    [
                        type2Parent
                        [
                            type2BoundaryLocalFace[0]
                        ]
                    ];

                const face& activeFace1 =
                    faces
                    [
                        type2Parent
                        [
                            type2BoundaryLocalFace[1]
                        ]
                    ];

                forAll(activeFace0, edge0I)
                {
                    const label a0 =
                        activeFace0[edge0I];

                    const label b0 =
                        activeFace0
                        [
                            (edge0I+1)
                          % activeFace0.size()
                        ];

                    const label lo0 = Foam::min(a0,b0);
                    const label hi0 = Foam::max(a0,b0);

                    forAll(activeFace1, edge1I)
                    {
                        const label a1 =
                            activeFace1[edge1I];

                        const label b1 =
                            activeFace1
                            [
                                (edge1I+1)
                              % activeFace1.size()
                            ];

                        if
                        (
                            lo0 == Foam::min(a1,b1)
                         && hi0 == Foam::max(a1,b1)
                        )
                        {
                            ++type2CommonEdgeCount;
                        }
                    }
                }

                if( type2CommonEdgeCount != 1 )
                {
                    type2ParentValid = false;
                    type2FailureReason =
                        "refinedFacesDoNotShareOneEdge";
                }
            }

            if( !type2ParentValid )
            {
                refinementValid_ = false;

                WarningIn
                (
                    "void refineBoundaryLayers::generateNewCells()"
                )
                    << "CFMITCH V3.6 EDGEHEX PREFLIGHT FAIL:"
                    << " parent=" << cellI
                    << " nFaces=" << type2Parent.size()
                    << " activeBoundaryFaces="
                    << type2ActiveBoundaryFaces
                    << " bfI=("
                    << type2BfI[0] << " "
                    << type2BfI[1] << ")"
                    << " commonEdges="
                    << type2CommonEdgeCount
                    << " reason="
                    << type2FailureReason
                    << " -- rejecting refinement transaction"
                    << endl;

                return;
            }

            //- generate new cells from a topologically validated hex
            //- where two boundary-layer directions intersect.
            refineEdgeHexCell refEdgeHex(cellI, *this);

            const DynList
            <
                DynList<DynList<label, 4>, 6>,
                256
            >& cellsFromCell =
                refEdgeHex.newCells();

            forAll(cellsFromCell, cI)
            {
                const DynList<DynList<label, 4>, 6>& nc =
                    cellsFromCell[cI];

                if
                (
                    !auditGeneratedChild
                    (
                        nc,
                        cellI,
                        cI,
                        refType[cellI]
                    )
                )
                {
                    // -------------------------------------------------
                    // CFMitch V3.7 -- first-failure edge-hex forensic.
                    //
                    // Diagnostic only.  V3.6 still rejects the child
                    // before commitment and restores Q0.
                    // -------------------------------------------------

                    const label forensicLayers0 =
                        (
                            type2BfI[0] >= 0
                         && type2BfI[0] <
                            label(nLayersAtBndFace_.size())
                        )
                      ? nLayersAtBndFace_[type2BfI[0]]
                      : -1;

                    const label forensicLayers1 =
                        (
                            type2BfI[1] >= 0
                         && type2BfI[1] <
                            label(nLayersAtBndFace_.size())
                        )
                      ? nLayersAtBndFace_[type2BfI[1]]
                      : -1;

                    Info
                        << "CFMITCH V3.7 EDGEHEX FORENSIC BEGIN:"
                        << " parent=" << cellI
                        << " childLocal=" << cI
                        << " parentFaces=" << type2Parent.size()
                        << " generatedChildren="
                        << cellsFromCell.size()
                        << " bfI=("
                        << type2BfI[0] << " "
                        << type2BfI[1] << ")"
                        << " localBoundaryFaces=("
                        << type2BoundaryLocalFace[0] << " "
                        << type2BoundaryLocalFace[1] << ")"
                        << " layers=("
                        << forensicLayers0 << " "
                        << forensicLayers1 << ")"
                        << endl;

                    forAll(type2Parent, forensicParentLocalFI)
                    {
                        const label forensicSourceFaceI =
                            type2Parent[forensicParentLocalFI];

                        const label forensicBfI =
                            forensicSourceFaceI
                          - type2StartBoundary;

                        label forensicRequestedLayers = -1;

                        if
                        (
                            forensicBfI >= 0
                         && forensicBfI <
                            label(nLayersAtBndFace_.size())
                        )
                        {
                            forensicRequestedLayers =
                                nLayersAtBndFace_
                                [
                                    forensicBfI
                                ];
                        }

                        Info
                            << "CFMITCH V3.7 EDGEHEX PARENT FACE:"
                            << " parent=" << cellI
                            << " localFace="
                            << forensicParentLocalFI
                            << " sourceFace="
                            << forensicSourceFaceI
                            << " bfI=" << forensicBfI
                            << " requestedLayers="
                            << forensicRequestedLayers;

                        if
                        (
                            forensicSourceFaceI >= 0
                         && forensicSourceFaceI <
                            label(faces.size())
                        )
                        {
                            Info
                                << " oldPoints="
                                << faces[forensicSourceFaceI];
                        }
                        else
                        {
                            Info << " oldPoints=OUT_OF_RANGE";
                        }

                        label forensicDerivedCount = 0;

                        if
                        (
                            forensicSourceFaceI >= 0
                         && forensicSourceFaceI <
                            label(facesFromFace_.size())
                        )
                        {
                            forensicDerivedCount =
                                facesFromFace_.sizeOfRow
                                (
                                    forensicSourceFaceI
                                );
                        }

                        Info
                            << " derivedCount="
                            << forensicDerivedCount
                            << endl;

                        if
                        (
                            forensicSourceFaceI < 0
                         || forensicSourceFaceI >=
                            label(facesFromFace_.size())
                        )
                        {
                            continue;
                        }

                        forAllRow
                        (
                            facesFromFace_,
                            forensicSourceFaceI,
                            forensicDerivedI
                        )
                        {
                            const label forensicDerivedFaceI =
                                facesFromFace_
                                (
                                    forensicSourceFaceI,
                                    forensicDerivedI
                                );

                            Info
                                << "CFMITCH V3.7 EDGEHEX DERIVED FACE:"
                                << " parent=" << cellI
                                << " sourceLocal="
                                << forensicParentLocalFI
                                << " sourceFace="
                                << forensicSourceFaceI
                                << " derivedLocal="
                                << forensicDerivedI
                                << " derivedFace="
                                << forensicDerivedFaceI;

                            if
                            (
                                forensicDerivedFaceI >= 0
                             && forensicDerivedFaceI <
                                label(newFaces_.size())
                            )
                            {
                                Info
                                    << " points="
                                    << newFaces_
                                       [
                                           forensicDerivedFaceI
                                       ];
                            }
                            else
                            {
                                Info << " points=OUT_OF_RANGE";
                            }

                            Info << endl;
                        }
                    }

                    forAll(nc, forensicChildLocalFI)
                    {
                        label forensicMatchedSource = -1;
                        label forensicMatchedDerived = -1;

                        forAll(type2Parent, forensicParentLocalFI)
                        {
                            const label forensicSourceFaceI =
                                type2Parent[forensicParentLocalFI];

                            if
                            (
                                forensicSourceFaceI < 0
                             || forensicSourceFaceI >=
                                label(facesFromFace_.size())
                            )
                            {
                                continue;
                            }

                            forAllRow
                            (
                                facesFromFace_,
                                forensicSourceFaceI,
                                forensicDerivedI
                            )
                            {
                                const label forensicDerivedFaceI =
                                    facesFromFace_
                                    (
                                        forensicSourceFaceI,
                                        forensicDerivedI
                                    );

                                if
                                (
                                    forensicDerivedFaceI >= 0
                                 && forensicDerivedFaceI <
                                    label(newFaces_.size())
                                 && help::areFacesEqual
                                    (
                                        nc[forensicChildLocalFI],
                                        newFaces_
                                        [
                                            forensicDerivedFaceI
                                        ]
                                    )
                                )
                                {
                                    forensicMatchedSource =
                                        forensicSourceFaceI;

                                    forensicMatchedDerived =
                                        forensicDerivedFaceI;
                                }
                            }
                        }

                        Info
                            << "CFMITCH V3.7 EDGEHEX CHILD FACE:"
                            << " parent=" << cellI
                            << " childLocal=" << cI
                            << " faceLocal="
                            << forensicChildLocalFI
                            << " matchedSource="
                            << forensicMatchedSource
                            << " matchedDerived="
                            << forensicMatchedDerived
                            << " points="
                            << nc[forensicChildLocalFI]
                            << endl;
                    }

                    std::map
                    <
                        std::pair<label,label>,
                        label
                    > forensicEdgeUse;

                    forAll(nc, forensicChildLocalFI)
                    {
                        const DynList<label, 4>& forensicFace =
                            nc[forensicChildLocalFI];

                        forAll(forensicFace, forensicPointI)
                        {
                            const label a =
                                forensicFace[forensicPointI];

                            const label b =
                                forensicFace
                                [
                                    (forensicPointI+1)
                                  % forensicFace.size()
                                ];

                            ++forensicEdgeUse
                            [
                                std::make_pair
                                (
                                    Foam::min(a,b),
                                    Foam::max(a,b)
                                )
                            ];
                        }
                    }

                    for
                    (
                        std::map
                        <
                            std::pair<label,label>,
                            label
                        >::const_iterator forensicEdgeIt =
                            forensicEdgeUse.begin();

                        forensicEdgeIt !=
                            forensicEdgeUse.end();
                        ++forensicEdgeIt
                    )
                    {
                        if( forensicEdgeIt->second == 2 )
                            continue;

                        const label forensicA =
                            forensicEdgeIt->first.first;

                        const label forensicB =
                            forensicEdgeIt->first.second;

                        Info
                            << "CFMITCH V3.7 EDGEHEX BAD EDGE:"
                            << " parent=" << cellI
                            << " childLocal=" << cI
                            << " edge=("
                            << forensicA << " "
                            << forensicB << ")"
                            << " use="
                            << forensicEdgeIt->second;

                        if
                        (
                            forensicA >= 0
                         && forensicB >= 0
                         && forensicA <
                            label(mesh_.points().size())
                         && forensicB <
                            label(mesh_.points().size())
                        )
                        {
                            Info
                                << " p0="
                                << mesh_.points()[forensicA]
                                << " p1="
                                << mesh_.points()[forensicB];
                        }

                        Info << " childFaces=(";

                        forAll(nc, forensicChildLocalFI)
                        {
                            const DynList<label, 4>& forensicFace =
                                nc[forensicChildLocalFI];

                            bool forensicContainsEdge = false;

                            forAll(forensicFace, forensicPointI)
                            {
                                const label a =
                                    forensicFace[forensicPointI];

                                const label b =
                                    forensicFace
                                    [
                                        (forensicPointI+1)
                                      % forensicFace.size()
                                    ];

                                if
                                (
                                    Foam::min(a,b) == forensicA
                                 && Foam::max(a,b) == forensicB
                                )
                                {
                                    forensicContainsEdge = true;
                                    break;
                                }
                            }

                            if( forensicContainsEdge )
                            {
                                Info << " "
                                     << forensicChildLocalFI;
                            }
                        }

                        Info << " )" << endl;
                    }

                    Info
                        << "CFMITCH V3.7 EDGEHEX FORENSIC END:"
                        << " parent=" << cellI
                        << " childLocal=" << cI
                        << endl;

                    refinementValid_ = false;

                    WarningIn
                    (
                        "void refineBoundaryLayers::generateNewCells()"
                    )
                        << "CFMITCH V3.6 EDGEHEX CHILD FAIL:"
                        << " parent=" << cellI
                        << " childLocal=" << cI
                        << " bfI=("
                        << type2BfI[0] << " "
                        << type2BfI[1] << ")"
                        << " -- rejecting before child commitment"
                        << endl;

                    return;
                }

                # ifdef DEBUGLayer
                Pout << "Adding cell " << (cI==0?cellI:nCells)
                     << " originating from cell " << cellI << endl;
                # endif

                const label newCellI = cI==0?cellI:nCells++;
                cellToBaseBndFace_[newCellI] = cellToBfI[cellI];

                newCellsFromCell.append(cellI, newCellI);

                cell& c = cells[newCellI];
                c.setSize(nc.size());

                //- find face labels for this cell
                forAll(nc, fI)
                {
                    const DynList<label, 4>& nf = nc[fI];

                    label faceLabel(-1);
                    forAllRow(pointNewFaces, nf[0], pfI)
                    {
                        const label nfI = pointNewFaces(nf[0], pfI);

                        if( help::areFacesEqual(nf, newFaces_[nfI]) )
                        {
                            c[fI] = nfI;
                            faceLabel = nfI;
                            break;
                        }
                    }

                    if( faceLabel < 0 )
                    {
                        forAll(nf, pI)
                            pointNewFaces.append(nf[pI], newFaces_.size());
                        c[fI] = newFaces_.size();
                        newFaces_.appendList(nf);
                    }
                }
            }
        }
        else if( refType[cellI] == 3 )
        {
            //- generate new cells from a hex at a corner where three
            //- layers intersect
            //- generate mostly hex cells
            refineCornerHexCell refCell(cellI, *this);
            const DynList<DynList<DynList<label, 4>, 6>, 256>& cellsFromCell =
                refCell.newCells();

            //- new points have been generated
            pointNewFaces.setSize(mesh_.points().size());

            //- recognise face cells in the graph of new faces
            forAll(cellsFromCell, cI)
            {
                const DynList<DynList<label, 4>, 6>& nc = cellsFromCell[cI];

                auditGeneratedChild(nc, cellI, cI, refType[cellI]);

                const label newCellI = cI==0?cellI:nCells++;
                cellToBaseBndFace_[newCellI] = cellToBfI[cellI];

                newCellsFromCell.append(cellI, newCellI);

                cell& c = cells[newCellI];
                c.setSize(nc.size());

                //- find face labels for this cell
                forAll(nc, fI)
                {
                    const DynList<label, 4>& nf = nc[fI];

                    label faceLabel(-1);
                    forAllRow(pointNewFaces, nf[0], pfI)
                    {
                        const label nfI = pointNewFaces(nf[0], pfI);

                        if( help::areFacesEqual(nf, newFaces_[nfI]) )
                        {
                            c[fI] = nfI;
                            faceLabel = nfI;
                            break;
                        }
                    }

                    if( faceLabel < 0 )
                    {
                        forAll(nf, pI)
                            pointNewFaces.append(nf[pI], newFaces_.size());
                        c[fI] = newFaces_.size();
                        newFaces_.appendList(nf);
                    }
                }
            }
        }
        else
        {
            FatalErrorIn
            (
                "void refineBoundaryLayers::generateNewCells()"
            ) << "Cannot refine boundary layer for cell "
              << cellI << abort(FatalError);
        }
    }

    Info << "REFINE_CHILD_CLOSURE"
         << " checked=" << nGeneratedChildrenChecked
         << " bad=" << nGeneratedChildrenBad
         << " badType1=" << nGeneratedBadType1
         << " badType2=" << nGeneratedBadType2
         << " badType3=" << nGeneratedBadType3
         << " badEdges=" << nGeneratedBadEdges
         << " malformedFaces=" << nGeneratedMalformedFaces
         << endl;

    if( nGeneratedChildrenBad > 0 )
    {
        refinementValid_ = false;

        WarningIn
        (
            "void refineBoundaryLayers::generateNewCells()"
        )
            << "CFMITCH V4.4 GENERATED CHILD CLOSURE REJECT:"
            << " badChildren=" << nGeneratedChildrenBad
            << " badType1=" << nGeneratedBadType1
            << " badType2=" << nGeneratedBadType2
            << " badType3=" << nGeneratedBadType3
            << " badEdges=" << nGeneratedBadEdges
            << " malformedFaces=" << nGeneratedMalformedFaces
            << " -- rejecting before face relabel"
            << endl;

        return;
    }

    // REFINE_PRE_RELABEL_CLOSURE_AUDIT
    //
    // Diagnostic only. At this point every cell already references
    // newFaces_, but the final face-list reconstruction/newFaceLabel
    // renumbering has NOT happened yet.
    //
    // This catches both:
    //   - refined cells after face consolidation, and
    //   - refType==0 cells reconstructed from facesFromFace_.
    {
        labelLongList badCellIds;

        label nBadCellsType0 = 0;
        label nBadCellsType1 = 0;
        label nBadCellsType2 = 0;
        label nBadCellsType3 = 0;
        label nBadCellsAppended = 0;

        label nBadEdgesTotal = 0;
        label nBadFaceRefs = 0;
        label nDegenerateFaces = 0;

        forAll(cells, cellI)
        {
            const cell& c = cells[cellI];

            std::map<std::pair<label,label>, label> edgeUse;

            bool bad = false;
            label badEdgesThisCell = 0;

            forAll(c, cfI)
            {
                const label faceI = c[cfI];

                if
                (
                    faceI < 0
                 || faceI >= label(newFaces_.size())
                )
                {
                    bad = true;
                    ++nBadFaceRefs;
                    continue;
                }

                const label nPts = newFaces_.sizeOfRow(faceI);

                if( nPts < 3 )
                {
                    bad = true;
                    ++nDegenerateFaces;
                    continue;
                }

                for(label pI=0; pI<nPts; ++pI)
                {
                    const label a = newFaces_(faceI, pI);
                    const label b =
                        newFaces_(faceI, (pI+1)%nPts);

                    if( a == b )
                        bad = true;

                    ++edgeUse
                    [
                        std::make_pair
                        (
                            Foam::min(a,b),
                            Foam::max(a,b)
                        )
                    ];
                }
            }

            for
            (
                std::map<std::pair<label,label>, label>::const_iterator
                    iter=edgeUse.begin();
                iter!=edgeUse.end();
                ++iter
            )
            {
                if( iter->second != 2 )
                {
                    bad = true;
                    ++badEdgesThisCell;
                }
            }

            if( bad )
            {
                badCellIds.append(cellI);
                nBadEdgesTotal += badEdgesThisCell;

                label rt = -1;

                if( cellI < label(refType.size()) )
                {
                    rt = refType[cellI];

                    if( rt == 0 )
                        ++nBadCellsType0;
                    else if( rt == 1 )
                        ++nBadCellsType1;
                    else if( rt == 2 )
                        ++nBadCellsType2;
                    else if( rt == 3 )
                        ++nBadCellsType3;
                }
                else
                {
                    ++nBadCellsAppended;
                }

                if( badCellIds.size() <= 100 )
                {
                    Info << "REFINE_PRE_RELABEL_BAD"
                         << " cell=" << cellI
                         << " refType=" << rt
                         << " nFaces=" << c.size()
                         << " badEdges=" << badEdgesThisCell
                         << endl;
                }
            }
        }

        Info << "REFINE_PRE_RELABEL_CLOSURE"
             << " cells=" << cells.size()
             << " badCells=" << badCellIds.size()
             << " badType0=" << nBadCellsType0
             << " badType1=" << nBadCellsType1
             << " badType2=" << nBadCellsType2
             << " badType3=" << nBadCellsType3
             << " badAppended=" << nBadCellsAppended
             << " badEdges=" << nBadEdgesTotal
             << " badFaceRefs=" << nBadFaceRefs
             << " degenerateFaces=" << nDegenerateFaces
             << endl;

        if( badCellIds.size() > 0 )
        {
            refinementValid_ = false;

            WarningIn
            (
                "void refineBoundaryLayers::generateNewCells()"
            )
                << "CFMITCH V4.4 PRE-RELABEL CLOSURE REJECT:"
                << " badCells=" << badCellIds.size()
                << " badType0=" << nBadCellsType0
                << " badType1=" << nBadCellsType1
                << " badType2=" << nBadCellsType2
                << " badType3=" << nBadCellsType3
                << " badAppended=" << nBadCellsAppended
                << " badEdges=" << nBadEdgesTotal
                << " badFaceRefs=" << nBadFaceRefs
                << " degenerateFaces=" << nDegenerateFaces
                << " -- rejecting before owner/neighbour relabel"
                << endl;

            return;
        }
    }

    //- check the orientation of new faces, because owner and neighbour cells
    //- may require a face to be flipped
    const label nOrigInternalFaces = mesh_.nInternalFaces();

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        const labelList& owner = mesh_.owner();
        const labelList& neighbour = mesh_.neighbour();

        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 40)
        # endif
        for(label fI=0;fI<nOrigInternalFaces;++fI)
        {
            const label own = owner[fI];
            const label nei = neighbour[fI];

            if( facesFromFace_.sizeOfRow(fI) == 1 )
                continue;

            forAllRow(facesFromFace_, fI, cfI)
            {
                const label nfI = facesFromFace_(fI, cfI);

                //- find the new owner and neighbour cells of the new face
                label newOwner(-1), newNeighbour(-1);
                forAllRow(newCellsFromCell, own, cI)
                {
                    const cell& cOwn = cells[newCellsFromCell(own, cI)];

                    const label pos = help::positionInList(nfI, cOwn);

                    if( pos >= 0 )
                    {
                        newOwner = newCellsFromCell(own, cI);
                        break;
                    }
                }

                forAllRow(newCellsFromCell, nei, cI)
                {
                    const cell& cNei = cells[newCellsFromCell(nei, cI)];

                    const label pos = help::positionInList(nfI, cNei);

                    if( pos >= 0 )
                    {
                        newNeighbour = newCellsFromCell(nei, cI);
                        break;
                    }
                }

                if( newOwner > newNeighbour )
                {
                    DynList<label> rf;
                    rf.setSize(newFaces_.sizeOfRow(nfI));

                    forAll(rf, i)
                        rf[i] = newFaces_(nfI, i);

                    rf = help::reverseFace(rf);

                    newFaces_.setRow(nfI, rf);
                }
            }
        }
    }

    //- update cell sets
    mesh_.updateCellSubsets(newCellsFromCell);
    newCellsFromCell.setSize(0);

    //- point-faces addressing is not needed any more
    pointNewFaces.setSize(0);

    //- copy newFaces to the mesh
    # ifdef DEBUGLayer
    Pout << "Copying internal faces " << endl;
    Pout << "Original number of internal faces " << nOrigInternalFaces << endl;
    # endif

    //- store internal faces originating from existing faces
    labelLongList newFaceLabel(newFaces_.size());
    faces.setSize(newFaces_.size());

    label currFace = 0;
    label nInternalRelabelMismatch = 0;

    for(label faceI=0;faceI<nOrigInternalFaces;++faceI)
    {
        forAllRow(facesFromFace_, faceI, ffI)
        {
            const label newFaceI = facesFromFace_(faceI, ffI);

            face& f = faces[currFace];

            // newFaceLabel maps an index in newFaces_ to its final
            // face index in the reconstructed mesh.  Using currFace
            // as both key and value is only valid accidentally when
            // newFaceI == currFace.
            if( newFaceI != currFace )
                ++nInternalRelabelMismatch;

            newFaceLabel[newFaceI] = currFace;
            ++currFace;

            f.setSize(newFaces_.sizeOfRow(newFaceI));

            forAll(f, pI)
                f[pI] = newFaces_(newFaceI, pI);
        }
    }

    Info << "REFINE_RELABEL internalDerivedMismatch="
         << nInternalRelabelMismatch
         << " internalDerivedFaces=" << currFace
         << endl;

    //- store newly-generated internal faces
    # ifdef DEBUGLayer
    Pout << "Copying newly generated internal faces" << endl;
    Pout << "nNewInternalFaces " << currFace << endl;
    Pout << "numFacesBefore " << numFacesBefore << endl;
    Pout << "Total number of faces " << newFaces_.size() << endl;
    # endif

    for(label faceI=numFacesBefore;faceI<newFaces_.size();++faceI)
    {
        newFaceLabel[faceI] = currFace;
        face& f = faces[currFace];
        ++currFace;

        f.setSize(newFaces_.sizeOfRow(faceI));

        forAll(f, pI)
            f[pI] = newFaces_(faceI, pI);
    }

    //- store new boundary faces
    # ifdef DEBUGLayer
    Pout << "Copying boundary faces " << endl;
    Pout << "currFace " << currFace << endl;
    Pout << "Faces size " << faces.size() << endl;
    Pout << "Initial number of faces " << facesFromFace_.size() << endl;
    # endif

    PtrList<boundaryPatch>& boundaries = meshModifier.boundariesAccess();
    forAll(boundaries, patchI)
    {
        const label start = boundaries[patchI].patchStart();
        const label size = boundaries[patchI].patchSize();

        const label newStart = currFace;
        label nNewFacesInPatch(0);
        for(label fI=0;fI<size;++fI)
        {
            const label faceI = start + fI;

            forAllRow(facesFromFace_, faceI, nfI)
            {
                face& f = faces[currFace];

                //- update the new label
                const label origFaceI = facesFromFace_(faceI, nfI);
                newFaceLabel[origFaceI] = currFace;
                facesFromFace_(faceI, nfI) = currFace;
                ++currFace;

                //- copy the face into the mesh
                f.setSize(newFaces_.sizeOfRow(origFaceI));
                forAll(f, pI)
                    f[pI] = newFaces_(origFaceI, pI);

                ++nNewFacesInPatch;
            }
        }

        //- update patch
        boundaries[patchI].patchStart() = newStart;
        boundaries[patchI].patchSize() = nNewFacesInPatch;
    }

    if( Pstream::parRun() )
    {
        # ifdef DEBUGLayer
        Pout << "Copying processor faces" << endl;
        # endif

        //- copy faces at inter-processor boundaries
        PtrList<processorBoundaryPatch>& procBoundaries =
            meshModifier.procBoundariesAccess();

        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const label size = procBoundaries[patchI].patchSize();

            const label newStart = currFace;
            label nNewFacesInPatch(0);
            for(label fI=0;fI<size;++fI)
            {
                const label faceI = start + fI;
                forAllRow(facesFromFace_, faceI, nfI)
                {
                    face& f = faces[currFace];

                    //- update the new label
                    const label origFaceI = facesFromFace_(faceI, nfI);
                    newFaceLabel[origFaceI] = currFace;
                    facesFromFace_(faceI, nfI) = currFace;
                    ++currFace;

                    //- copy the face into the mesh
                    f.setSize(newFaces_.sizeOfRow(origFaceI));
                    forAll(f, pI)
                        f[pI] = newFaces_(origFaceI, pI);

                    ++nNewFacesInPatch;
                }
            }

            //- update patch
            procBoundaries[patchI].patchStart() = newStart;
            procBoundaries[patchI].patchSize() = nNewFacesInPatch;
        }
    }

    # ifdef DEBUGLayer
    Pout << "Faces after refinement " << faces << endl;
    Pout << "newFaceLabel " << newFaceLabel << endl;
    # endif

    //- update face subsets
    mesh_.updateFaceSubsets(facesFromFace_);
    facesFromFace_.setSize(0);
    newFaces_.setSize(0);

    //- update cells to match the faces
    # ifdef DEBUGLayer
    Pout << "Updating cells to match new faces" << endl;
    # endif

    forAll(cells, cellI)
    {
        cell& c = cells[cellI];

        forAll(c, fI)
            c[fI] = newFaceLabel[c[fI]];
    }

    // REFINE_POST_RELABEL_CLOSURE_AUDIT
    // Diagnostic only. At this point cells reference the reconstructed
    // mesh face list, so test the actual committed cell shells.
    {
        labelLongList badCellIds;
        label nBadEdgesTotal = 0;
        label nBadFaceRefs = 0;
        label nDegenerateFaces = 0;

        forAll(cells, cellI)
        {
            const cell& c = cells[cellI];

            std::map<std::pair<label,label>, label> edgeUse;

            bool bad = false;
            label badEdgesThisCell = 0;

            forAll(c, cfI)
            {
                const label faceI = c[cfI];

                if( faceI < 0 || faceI >= faces.size() )
                {
                    bad = true;
                    ++nBadFaceRefs;
                    continue;
                }

                const face& f = faces[faceI];

                if( f.size() < 3 )
                {
                    bad = true;
                    ++nDegenerateFaces;
                    continue;
                }

                forAll(f, pI)
                {
                    const label a = f[pI];
                    const label b = f[(pI+1)%f.size()];

                    if( a == b )
                        bad = true;

                    ++edgeUse
                    [
                        std::make_pair
                        (
                            Foam::min(a,b),
                            Foam::max(a,b)
                        )
                    ];
                }
            }

            for
            (
                std::map<std::pair<label,label>, label>::const_iterator
                    iter=edgeUse.begin();
                iter!=edgeUse.end();
                ++iter
            )
            {
                if( iter->second != 2 )
                {
                    bad = true;
                    ++badEdgesThisCell;
                }
            }

            if( bad )
            {
                badCellIds.append(cellI);
                nBadEdgesTotal += badEdgesThisCell;
            }
        }

        Info << "REFINE_POST_RELABEL_CLOSURE"
             << " cells=" << cells.size()
             << " badCells=" << badCellIds.size()
             << " badEdges=" << nBadEdgesTotal
             << " badFaceRefs=" << nBadFaceRefs
             << " degenerateFaces=" << nDegenerateFaces
             << endl;

        if( badCellIds.size() )
        {
            Info << "REFINE_POST_RELABEL_BAD_IDS ids=(";

            const label nPrint =
                Foam::min(label(badCellIds.size()), label(100));

            for(label i=0; i<nPrint; ++i)
            {
                if( i ) Info << ',';
                Info << badCellIds[i];
            }

            if( badCellIds.size() > nPrint )
                Info << ",...";

            Info << ')' << endl;
        }
    }

    # ifdef DEBUGLayer
    Pout << "Cleaning mesh " << endl;
    # endif

    //- delete all adressing which is no longer up-to-date

    // ==================================================================
    // CFMITCH V4.5 FINAL FACE-INCIDENCE GATE
    //
    // Cell-shell closure alone does not prove that the final global
    // topology is valid.  Every final internal face must appear in
    // exactly two cells and every boundary face in exactly one.
    //
    // This check runs after final face relabeling and cell-face remapping
    // but before addressing caches are rebuilt from the committed mesh.
    // ==================================================================
    {
        labelLongList finalFaceUse(faces.size());
        labelLongList finalFirstCell(faces.size());

        forAll(finalFaceUse, faceI)
        {
            finalFaceUse[faceI] = 0;
            finalFirstCell[faceI] = -1;
        }

        label nBadCellFaceRefs = 0;

        forAll(cells, cellI)
        {
            const cell& c = cells[cellI];

            forAll(c, cfI)
            {
                const label faceI = c[cfI];

                if
                (
                    faceI < 0
                 || faceI >= label(faces.size())
                )
                {
                    ++nBadCellFaceRefs;
                    continue;
                }

                if( finalFaceUse[faceI] == 0 )
                    finalFirstCell[faceI] = cellI;

                ++finalFaceUse[faceI];
            }
        }

        const label finalBoundaryStart =
            boundaries.size()
          ? boundaries[0].patchStart()
          : label(faces.size());

        label nInternalUseBad = 0;
        label nBoundaryUseBad = 0;

        label nUse0 = 0;
        label nUse1 = 0;
        label nUse2 = 0;
        label nUse3Plus = 0;

        label nBadPrinted = 0;
        label firstBadInternalFace = -1;
        label firstBadBoundaryFace = -1;

        forAll(finalFaceUse, faceI)
        {
            const label nUse = finalFaceUse[faceI];

            if( nUse == 0 )
                ++nUse0;
            else if( nUse == 1 )
                ++nUse1;
            else if( nUse == 2 )
                ++nUse2;
            else
                ++nUse3Plus;

            const label expectedUse =
                faceI < finalBoundaryStart ? 2 : 1;

            if( nUse == expectedUse )
                continue;

            if( faceI < finalBoundaryStart )
            {
                ++nInternalUseBad;

                if( firstBadInternalFace < 0 )
                    firstBadInternalFace = faceI;
            }
            else
            {
                ++nBoundaryUseBad;

                if( firstBadBoundaryFace < 0 )
                    firstBadBoundaryFace = faceI;
            }

            if( nBadPrinted < 40 )
            {
                Info
                    << "CFMITCH V4.5 FACE INCIDENCE BAD:"
                    << " face=" << faceI
                    << " class="
                    << (
                           faceI < finalBoundaryStart
                         ? "internal"
                         : "boundary"
                       )
                    << " expectedUse=" << expectedUse
                    << " actualUse=" << nUse
                    << " firstCell=" << finalFirstCell[faceI]
                    << " points=" << faces[faceI]
                    << endl;

                ++nBadPrinted;
            }
        }

        const bool finalFaceIncidenceValid =
            nBadCellFaceRefs == 0
         && nInternalUseBad == 0
         && nBoundaryUseBad == 0
         && currFace == label(faces.size());

        Info
            << "CFMITCH V4.5 FINAL FACE INCIDENCE:"
            << " valid="
            << (finalFaceIncidenceValid ? "yes" : "no")
            << " faces=" << faces.size()
            << " currFace=" << currFace
            << " boundaryStart=" << finalBoundaryStart
            << " use0=" << nUse0
            << " use1=" << nUse1
            << " use2=" << nUse2
            << " use3plus=" << nUse3Plus
            << " badInternal=" << nInternalUseBad
            << " badBoundary=" << nBoundaryUseBad
            << " badCellFaceRefs=" << nBadCellFaceRefs
            << " firstBadInternal=" << firstBadInternalFace
            << " firstBadBoundary=" << firstBadBoundaryFace
            << endl;

        if( !finalFaceIncidenceValid )
        {
            refinementValid_ = false;

            WarningIn
            (
                "void refineBoundaryLayers::generateNewCells()"
            )
                << "CFMITCH V4.5 FINAL FACE INCIDENCE REJECT:"
                << " refusing to commit globally inconsistent "
                << "owner/neighbour topology"
                << endl;

            return;
        }
    }

    meshModifier.clearAll();
    deleteDemandDrivenData(msePtr_);

    // REFINE_EXACT_VOLUME_BIRTH_AUDIT
    //
    // Diagnostic only.
    //
    // Reproduce polyMeshGenChecks::checkCellVolumes() exactly, immediately
    // after generateNewCells() has finalized face labels, winding, and
    // owner/neighbour topology.
    //
    // This is deliberately NOT based on the positive/clamped addressing
    // cellVolumes() cache.
    {
        const vectorField& exactFCtrs =
            mesh_.addressingData().faceCentres();

        const vectorField& exactFAreas =
            mesh_.addressingData().faceAreas();

        const labelList& exactOwner =
            mesh_.owner();

        const cellListPMG& exactCells =
            mesh_.cells();

        label nBelowVSmall = 0;
        label nTrueNegative = 0;

        label nType1Checked = 0;
        label nType1BelowVSmall = 0;
        label nType1TrueNegative = 0;

        label nDetailedPrinted = 0;

        // BL_VALIDITY_REPAIR_V1C
        // True-negative generated type-1 children discovered by the
        // existing exact-volume parity scan.
        DynList<label, 32> blV1NegativeType1Cells;

        scalar minVolume = GREAT;
        scalar minType1Volume = GREAT;

        forAll(exactCells, exactCellI)
        {
            const cell& ec = exactCells[exactCellI];

            if( ec.size() == 0 )
                continue;

            vector cEst(vector::zero);

            forAll(ec, fI)
                cEst += exactFCtrs[ec[fI]];

            cEst /= ec.size();

            scalar cellVol = scalar(0);

            forAll(ec, fI)
            {
                scalar pyr3Vol =
                    exactFAreas[ec[fI]]
                  & (exactFCtrs[ec[fI]] - cEst);

                if( exactOwner[ec[fI]] != exactCellI )
                    pyr3Vol *= scalar(-1);

                cellVol += pyr3Vol;
            }

            cellVol /= scalar(3);

            minVolume =
                Foam::min(minVolume, cellVol);

            if( cellVol < VSMALL )
                ++nBelowVSmall;

            if( cellVol < scalar(0) )
                ++nTrueNegative;


            const bool isGeneratedType1 =
            (
                exactCellI >= 0
             && exactCellI < label(exactVolumeRefType.size())
             && exactVolumeRefType[exactCellI] == 1
            );

            if( !isGeneratedType1 )
                continue;

            ++nType1Checked;

            minType1Volume =
                Foam::min(minType1Volume, cellVol);

            if( cellVol < VSMALL )
                ++nType1BelowVSmall;

            if( cellVol < scalar(0) )
            {
                ++nType1TrueNegative;
                blV1NegativeType1Cells.append(exactCellI);
            }


            if
            (
                cellVol < VSMALL
             && nDetailedPrinted < 20
            )
            {
                ++nDetailedPrinted;

                const label parentCell =
                    exactVolumeParent[exactCellI];

                const label localChild =
                    exactVolumeLocalChild[exactCellI];

                label bfI = -1;

                if
                (
                    exactCellI >= 0
                 && exactCellI
                    < label(cellToBaseBndFace_.size())
                )
                    bfI =
                        cellToBaseBndFace_[exactCellI];

                label patchI = -1;

                if
                (
                    bfI >= 0
                 && bfI < label(exactVolumeFacePatch.size())
                )
                    patchI =
                        exactVolumeFacePatch[bfI];

                word patchName("?");

                if
                (
                    patchI >= 0
                 && patchI
                    < label(mesh_.boundaries().size())
                )
                    patchName =
                        mesh_.boundaries()[patchI].patchName();


                label nLayers = -1;

                if
                (
                    bfI >= 0
                 && bfI < label(nLayersAtBndFace_.size())
                )
                    nLayers =
                        nLayersAtBndFace_[bfI];


                Info
                    << "REFINE_EXACT_VOLUME_BAD"
                    << " cell=" << exactCellI
                    << " parent=" << parentCell
                    << " localChild=" << localChild
                    << " refType=1"
                    << " bfI=" << bfI
                    << " patch=" << patchName
                    << " nLayers=" << nLayers
                    << " volume=" << cellVol
                    << " nFaces=" << ec.size()
                    << " cEst=" << cEst
                    << endl;


                // Print the exact signed 3*pyramid contributions used in
                // checkCellVolumes().  Their sum / 3 is cellVol.
                Info
                    << "REFINE_EXACT_VOLUME_PYR3"
                    << " cell=" << exactCellI
                    << " contributions=(";

                forAll(ec, fI)
                {
                    const label faceI = ec[fI];

                    scalar pyr3Vol =
                        exactFAreas[faceI]
                      & (exactFCtrs[faceI] - cEst);

                    if( exactOwner[faceI] != exactCellI )
                        pyr3Vol *= scalar(-1);

                    if( fI )
                        Info << ' ';

                    Info
                        << faceI
                        << ':'
                        << pyr3Vol;
                }

                Info << ')' << endl;
            }
        }


        // Match checkCellVolumes() parallel semantics for the summary.
        reduce(minVolume, minOp<scalar>());
        reduce(nBelowVSmall, sumOp<label>());
        reduce(nTrueNegative, sumOp<label>());

        reduce(minType1Volume, minOp<scalar>());
        reduce(nType1Checked, sumOp<label>());
        reduce(nType1BelowVSmall, sumOp<label>());
        reduce(nType1TrueNegative, sumOp<label>());


        Info
            << "REFINE_EXACT_VOLUME_SUMMARY"
            << " allCells=" << exactCells.size()
            << " belowVSmall=" << nBelowVSmall
            << " trueNegative=" << nTrueNegative
            << " minVolume=" << minVolume
            << " type1Checked=" << nType1Checked
            << " type1BelowVSmall=" << nType1BelowVSmall
            << " type1TrueNegative=" << nType1TrueNegative
            << " minType1Volume=" << minType1Volume
            << endl;


        // ==============================================================
        // BL_VALIDITY_REPAIR_V1C
        //
        // Coherent split-edge chain repair.
        //
        // v0 moved individual generated row points independently.  That
        // can create a kink between adjacent BL rows and visually jagged
        // contact lines.
        //
        // v1c retains the v1b first-row-focused coherent hair repair:
        //
        //     deltaT(row) =
        //         amplitude * firstT
        //       * ((1-t)/(1-firstT))^4
        //
        // Therefore:
        //     wall endpoint t=0 is fixed
        //     first generated row receives maximum displacement
        //     deeper-row displacement monotonically decays
        //     BL/core endpoint t=1 is fixed
        //
        // v1c additionally evaluates every candidate using OpenFOAM-parity
        // face centres, face areas, signed cell centres, pyramids,
        // orthogonality and skewness.  All admissible candidate hairs and
        // amplitudes are searched before the best safe candidate is chosen.
        //
        // A complete split-edge chain is committed atomically.
        // ==============================================================

        label blV1InitialNegative =
            blV1NegativeType1Cells.size();

        label blV1Fixed = 0;
        label blV1Unresolved = 0;
        label blV1CommittedChains = 0;

        // CFMitch v2.5:
        //
        // V1C was historically gated by true-negative type-1 children.
        // That prevents positive-volume BL children with invalid oriented
        // face pyramids from ever reaching the quality-aware hair solver.
        //
        // Enter unconditionally so a repair population can also be built
        // from CURRENT bad-pyramid geometry once split-edge lineage exists.
        label blV25PyramidAdditionalSeeds = 0;
        label blV25RepairSeeds = 0;
        label blV25TargetPyrRejects = 0;

        // CFMitch v2.7 diagnostic:
        // classify the bad OpenFOAM pyramid faces belonging to unresolved
        // wall-adjacent type-1 children.
        //
        // No behaviour depends on these counters.
        label blV27WallChildCells = 0;
        label blV27WallChildBadFaces = 0;
        label blV27WallChildWallBaseBad = 0;
        label blV27WallChildInternalOuterBad = 0;
        label blV27WallChildLateralBad = 0;
        label blV27WallChildUnknownBad = 0;
        label blV27WallChildMultiBad = 0;

        // CFMitch v2.7.1:
        // bounded coherent wall-face repair for positive-volume,
        // type-1, wall-adjacent, WALL_BASE-only pyramid defects.
        label blV271FaceBreathEligible = 0;
        label blV271FaceBreathAttempted = 0;
        label blV271FaceBreathSkipped = 0;
        label blV271FaceBreathTrials = 0;
        label blV271FaceBreathVolumePass = 0;
        label blV271FaceBreathQualityPass = 0;
        label blV271FaceBreathTargetReject = 0;
        label blV271FaceBreathQualityReject = 0;
        label blV271FaceBreathFixed = 0;

        // Diagnostic-only quality-rejection classifier.
        label blV271RejectNewPyramid = 0;
        label blV271RejectWorsePyramid = 0;
        label blV271RejectNewSkew = 0;
        label blV271RejectWorseSkew = 0;
        label blV271RejectNewNonOrtho90 = 0;
        label blV271RejectWorseNonOrtho = 0;
        label blV271RejectOFGeometry = 0;
        label blV271RejectOther = 0;

        label blV271RejectOnTargetCellFace = 0;
        label blV271RejectOnOtherStarFace = 0;

        // Diagnostic-only:
        // one record for each unique
        //
        //     (repair target, moved split edge, worse pyramid face)
        //
        // No repair behavior depends on this state.
        std::set
        <
            std::pair
            <
                label,
                std::pair<label,label>
            >
        > blV1WorsePyrDiagSeen;

        {
            pointFieldPMG& v1Points =
                mesh_.points();

            const faceListPMG& v1Faces =
                mesh_.faces();

            const cellListPMG& v1Cells =
                mesh_.cells();

            const labelList& v1Neighbour =
                mesh_.neighbour();


            // ----------------------------------------------------------
            // Exact local face-centre and area-vector calculation.
            // Mirrors polyMeshGenAddressing::makeFaceCentresAndAreas().
            // ----------------------------------------------------------

            auto v1FaceCentreArea =
            [&]
            (
                const face& f,
                vector& fCtr,
                vector& fArea
            )
            {
                const label nPoints = f.size();

                if( nPoints == 3 )
                {
                    const point& p0 = v1Points[f[0]];
                    const point& p1 = v1Points[f[1]];
                    const point& p2 = v1Points[f[2]];

                    fCtr =
                        (scalar(1)/scalar(3))
                       *(p0+p1+p2);

                    fArea =
                        scalar(0.5)
                       *((p1-p0)^(p2-p0));

                    return;
                }

                vector sumN(vector::zero);
                scalar sumA = scalar(0);
                vector sumAc(vector::zero);

                point fCentre =
                    v1Points[f[0]];

                for(label pi=1; pi<nPoints; ++pi)
                    fCentre += v1Points[f[pi]];

                fCentre /= scalar(nPoints);

                for(label pi=0; pi<nPoints; ++pi)
                {
                    const point& curr =
                        v1Points[f[pi]];

                    const point& next =
                        v1Points[f.nextLabel(pi)];

                    const vector c =
                        curr + next + fCentre;

                    const vector n =
                        (next-curr)^(fCentre-curr);

                    const scalar a =
                        mag(n);

                    sumN += n;
                    sumA += a;
                    sumAc += a*c;
                }

                fCtr =
                    (scalar(1)/scalar(3))
                   *sumAc
                   /(sumA + VSMALL);

                fArea =
                    scalar(0.5)*sumN;
            };


            // ----------------------------------------------------------
            // Exact local raw signed cell volume.
            // Mirrors polyMeshGenChecks::checkCellVolumes().
            // ----------------------------------------------------------

            auto v1CellVolume =
            [&]
            (
                const label cellI
            ) -> scalar
            {
                if
                (
                    cellI < 0
                 || cellI >= label(v1Cells.size())
                )
                    return -GREAT;

                const cell& c =
                    v1Cells[cellI];

                if( c.size() == 0 )
                    return -GREAT;

                vector cEst(vector::zero);

                forAll(c, cfI)
                {
                    const label faceI =
                        c[cfI];

                    if
                    (
                        faceI < 0
                     || faceI >= label(v1Faces.size())
                    )
                        return -GREAT;

                    vector fc(vector::zero);
                    vector fa(vector::zero);

                    v1FaceCentreArea
                    (
                        v1Faces[faceI],
                        fc,
                        fa
                    );

                    cEst += fc;
                }

                cEst /= scalar(c.size());

                scalar cellVol =
                    scalar(0);

                forAll(c, cfI)
                {
                    const label faceI =
                        c[cfI];

                    vector fc(vector::zero);
                    vector fa(vector::zero);

                    v1FaceCentreArea
                    (
                        v1Faces[faceI],
                        fc,
                        fa
                    );

                    scalar pyr3Vol =
                        fa & (fc-cEst);

                    if( exactOwner[faceI] != cellI )
                        pyr3Vol *= scalar(-1);

                    cellVol += pyr3Vol;
                }

                return cellVol/scalar(3);
            };


            // ----------------------------------------------------------
            // OpenFOAM-parity local face centre and area vector.
            //
            // Mirrors face::areaAndCentre() used by primitiveMesh.
            //
            // This deliberately differs from v1FaceCentreArea() for
            // general polygon faces: OpenFOAM weights triangle centres
            // using signed area projected onto the resultant face normal,
            // not triangle-area magnitude.
            // ----------------------------------------------------------

            auto v1OFFaceCentreArea =
            [&]
            (
                const face& f,
                vector& fCtr,
                vector& fArea
            ) -> bool
            {
                const label nPoints = f.size();

                if( nPoints < 3 )
                    return false;

                if( nPoints == 3 )
                {
                    const point& p0 = v1Points[f[0]];
                    const point& p1 = v1Points[f[1]];
                    const point& p2 = v1Points[f[2]];

                    fArea =
                        scalar(0.5)
                       *((p1-p0)^(p2-p0));

                    fCtr =
                        (scalar(1)/scalar(3))
                       *(p0+p1+p2);

                    return true;
                }

                point pAvg(vector::zero);

                forAll(f, pi)
                    pAvg += v1Points[f[pi]];

                pAvg /= scalar(nPoints);

                vector sumA(vector::zero);

                forAll(f, pi)
                {
                    const point& fp =
                        v1Points[f[pi]];

                    const point& fpNext =
                        v1Points[f.nextLabel(pi)];

                    const vector a =
                        (fpNext-fp)^(pAvg-fp);

                    sumA += a;
                }

                const vector sumAHat =
                    normalised(sumA);

                scalar sumAn = scalar(0);
                vector sumAnc(vector::zero);

                forAll(f, pi)
                {
                    const point& fp =
                        v1Points[f[pi]];

                    const point& fpNext =
                        v1Points[f.nextLabel(pi)];

                    const vector a =
                        (fpNext-fp)^(pAvg-fp);

                    const vector c =
                        fp + fpNext + pAvg;

                    const scalar an =
                        a & sumAHat;

                    sumAn += an;
                    sumAnc += an*c;
                }

                fArea =
                    scalar(0.5)*sumA;

                if( sumAn > vSmall )
                {
                    fCtr =
                        (scalar(1)/scalar(3))
                       *sumAnc/sumAn;
                }
                else
                {
                    fCtr = pAvg;
                }

                return true;
            };


            // ----------------------------------------------------------
            // OpenFOAM-parity signed cell centre and volume.
            //
            // Mirrors primitiveMesh::makeCellCentresAndVols().
            //
            // IMPORTANT:
            // This is ONLY a local quality evaluator.  It must NOT replace
            // cfMesh's globally cached bounded surrogate cell centre,
            // which intentionally clamps negative pyramid weights to keep
            // defective intermediate meshes numerically finite.
            // ----------------------------------------------------------

            auto v1OFCellCentreVolume =
            [&]
            (
                const label cellI,
                point& cellCtr,
                scalar& cellVol
            ) -> bool
            {
                if
                (
                    cellI < 0
                 || cellI >= label(v1Cells.size())
                )
                    return false;

                const cell& c =
                    v1Cells[cellI];

                if( c.empty() )
                    return false;

                point cEst(vector::zero);

                forAll(c, cfI)
                {
                    const label faceI =
                        c[cfI];

                    if
                    (
                        faceI < 0
                     || faceI >= label(v1Faces.size())
                    )
                        return false;

                    vector fc(vector::zero);
                    vector fa(vector::zero);

                    if
                    (
                        !v1OFFaceCentreArea
                        (
                            v1Faces[faceI],
                            fc,
                            fa
                        )
                    )
                        return false;

                    cEst += fc;
                }

                cEst /= scalar(c.size());

                vector weightedCentre(vector::zero);
                scalar vol3 = scalar(0);

                forAll(c, cfI)
                {
                    const label faceI =
                        c[cfI];

                    vector fc(vector::zero);
                    vector fa(vector::zero);

                    if
                    (
                        !v1OFFaceCentreArea
                        (
                            v1Faces[faceI],
                            fc,
                            fa
                        )
                    )
                        return false;

                    scalar pyr3Vol =
                        fa & (fc-cEst);

                    if( exactOwner[faceI] != cellI )
                        pyr3Vol *= scalar(-1);

                    const vector pc =
                        scalar(0.75)*fc
                      + scalar(0.25)*cEst;

                    weightedCentre +=
                        pyr3Vol*pc;

                    vol3 += pyr3Vol;
                }

                if( Foam::mag(vol3) > vSmall )
                    cellCtr = weightedCentre/vol3;
                else
                    cellCtr = cEst;

                cellVol =
                    vol3/scalar(3);

                return true;
            };


            // ----------------------------------------------------------
            // OpenFOAM-parity single-face finite-volume quality.
            //
            // Internal-face skew mirrors
            // meshCheck::faceSkewness().
            //
            // Boundary-face skew mirrors
            // meshCheck::boundaryFaceSkewness().
            //
            // pyrMargin:
            //     owner    = -ownerPyramidVolume
            //     neighbour= +neighbourPyramidVolume
            //
            // Larger pyrMargin is better.  checkMesh considers the face
            // bad when either oriented pyramid margin is below -SMALL.
            // ----------------------------------------------------------

            auto v1OFFaceQuality =
            [&]
            (
                const label faceI,
                scalar& skewness,
                scalar& orthogonality,
                scalar& pyrMargin,
                bool& badPyramid,
                scalar& ownerVol,
                scalar& neighbourVol
            ) -> bool
            {
                if
                (
                    faceI < 0
                 || faceI >= label(v1Faces.size())
                )
                    return false;

                vector fc(vector::zero);
                vector fa(vector::zero);

                if
                (
                    !v1OFFaceCentreArea
                    (
                        v1Faces[faceI],
                        fc,
                        fa
                    )
                )
                    return false;

                const label ownCellI =
                    exactOwner[faceI];

                point ownCc(vector::zero);

                if
                (
                    !v1OFCellCentreVolume
                    (
                        ownCellI,
                        ownCc,
                        ownerVol
                    )
                )
                    return false;

                const scalar ownerPyrVol =
                    pyramidPointFaceRef
                    (
                        v1Faces[faceI],
                        ownCc
                    ).mag(v1Points);

                pyrMargin =
                    -ownerPyrVol;

                badPyramid =
                    ownerPyrVol > SMALL;

                const label neiCellI =
                    v1Neighbour[faceI];

                if( neiCellI >= 0 )
                {
                    point neiCc(vector::zero);

                    if
                    (
                        !v1OFCellCentreVolume
                        (
                            neiCellI,
                            neiCc,
                            neighbourVol
                        )
                    )
                        return false;

                    const scalar neighbourPyrVol =
                        pyramidPointFaceRef
                        (
                            v1Faces[faceI],
                            neiCc
                        ).mag(v1Points);

                    pyrMargin =
                        Foam::min
                        (
                            pyrMargin,
                            neighbourPyrVol
                        );

                    if( neighbourPyrVol < -SMALL )
                        badPyramid = true;

                    const vector Cpf =
                        fc-ownCc;

                    const vector d =
                        neiCc-ownCc;

                    const vector sv =
                        Cpf
                      - (
                            (fa & Cpf)
                           /((fa & d) + rootVSmall)
                        )*d;

                    const vector svHat =
                        sv/(mag(sv) + rootVSmall);

                    scalar fd =
                        scalar(0.2)*mag(d)
                      + rootVSmall;

                    const face& f =
                        v1Faces[faceI];

                    forAll(f, pi)
                    {
                        fd =
                            Foam::max
                            (
                                fd,
                                Foam::mag
                                (
                                    svHat
                                  & (
                                        v1Points[f[pi]]
                                       -fc
                                    )
                                )
                            );
                    }

                    skewness =
                        mag(sv)/fd;

                    orthogonality =
                        (d & fa)
                       /(mag(d)*mag(fa) + rootVSmall);
                }
                else
                {
                    neighbourVol = GREAT;

                    const vector Cpf =
                        fc-ownCc;

                    vector normal = fa;
                    normal /=
                        mag(normal) + rootVSmall;

                    const vector d =
                        normal*(normal & Cpf);

                    const vector sv =
                        Cpf
                      - (
                            (fa & Cpf)
                           /((fa & d) + rootVSmall)
                        )*d;

                    const vector svHat =
                        sv/(mag(sv) + rootVSmall);

                    scalar fd =
                        scalar(0.4)*mag(d)
                      + rootVSmall;

                    const face& f =
                        v1Faces[faceI];

                    forAll(f, pi)
                    {
                        fd =
                            Foam::max
                            (
                                fd,
                                Foam::mag
                                (
                                    svHat
                                  & (
                                        v1Points[f[pi]]
                                       -fc
                                    )
                                )
                            );
                    }

                    skewness =
                        mag(sv)/fd;

                    // Non-orthogonality is an owner-neighbour metric.
                    orthogonality = GREAT;
                }

                return true;
            };


            // ----------------------------------------------------------
            // point -> split-edge map for generated interior vertices.
            // Ambiguous points are not legal repair variables.
            // ----------------------------------------------------------

            std::map<label,label> v1PointEdge;
            std::set<label> v1AmbiguousPoints;

            for
            (
                label seI=0;
                seI<label(splitEdges_.size());
                ++seI
            )
            {
                const label rowSize =
                    newVerticesForSplitEdge_.
                        sizeOfRow(seI);

                for
                (
                    label rowI=1;
                    rowI<rowSize-1;
                    ++rowI
                )
                {
                    const label pointI =
                        newVerticesForSplitEdge_
                        (
                            seI,
                            rowI
                        );

                    std::map<label,label>::iterator it =
                        v1PointEdge.find(pointI);

                    if( it == v1PointEdge.end() )
                    {
                        v1PointEdge[pointI] = seI;
                    }
                    else if( it->second != seI )
                    {
                        v1AmbiguousPoints.insert(pointI);
                    }
                }
            }

            for
            (
                std::set<label>::const_iterator
                    it=v1AmbiguousPoints.begin();
                it!=v1AmbiguousPoints.end();
                ++it
            )
            {
                v1PointEdge.erase(*it);
            }


            // ----------------------------------------------------------
            // CFMitch v2.5 cell-side pyramid predicate.
            //
            // This deliberately tests the pyramid orientation relative to
            // THIS cell, rather than assigning both sides of a globally bad
            // face to the repair solver.
            //
            // owner:
            //     pyramidPointFaceRef(...).mag() > SMALL  => bad
            //
            // neighbour:
            //     pyramidPointFaceRef(...).mag() < -SMALL => bad
            //
            // This is the same orientation convention used by
            // v1OFFaceQuality() above.
            // ----------------------------------------------------------

            auto v1CellBadPyramidCount =
            [&]
            (
                const label cellI
            ) -> label
            {
                if
                (
                    cellI < 0
                 || cellI >= label(v1Cells.size())
                )
                    return labelMax;

                point cellCtr(vector::zero);
                scalar cellVol = -GREAT;

                if
                (
                    !v1OFCellCentreVolume
                    (
                        cellI,
                        cellCtr,
                        cellVol
                    )
                )
                    return labelMax;

                const cell& c =
                    v1Cells[cellI];

                label nBad = 0;

                forAll(c, cfI)
                {
                    const label faceI =
                        c[cfI];

                    if
                    (
                        faceI < 0
                     || faceI >= label(v1Faces.size())
                    )
                        return labelMax;

                    const scalar pyrVol =
                        pyramidPointFaceRef
                        (
                            v1Faces[faceI],
                            cellCtr
                        ).mag(v1Points);

                    if
                    (
                        faceI < label(exactOwner.size())
                     && exactOwner[faceI] == cellI
                    )
                    {
                        if( pyrVol > SMALL )
                            ++nBad;
                    }
                    else if
                    (
                        faceI < label(v1Neighbour.size())
                     && v1Neighbour[faceI] == cellI
                    )
                    {
                        if( pyrVol < -SMALL )
                            ++nBad;
                    }
                    else
                    {
                        // Cell/face addressing disagreement means this is
                        // not a safe V1C quality target.
                        return labelMax;
                    }
                }

                return nBad;
            };


            // ----------------------------------------------------------
            // CFMitch v2.7 diagnostic:
            // face-role classification for an unresolved wall-adjacent
            // type-1 child.
            //
            // Hair orientation has already been established by
            // detectBoundaryLayers:
            //
            //     splitEdges_[seI].start() = boundary/root point
            //     splitEdges_[seI].end()   = interior/core-side point
            //
            // This diagnostic deliberately uses the same cell-side
            // OpenFOAM pyramid predicate as v1CellBadPyramidCount().
            // ----------------------------------------------------------

            labelHashSet v1HairStartPoints;

            forAll(splitEdges_, seI)
            {
                v1HairStartPoints.insert
                (
                    splitEdges_[seI].start()
                );
            }


            auto v1ReportWallChildBadFaceRoles =
            [&]
            (
                const label cellI
            )
            {
                if
                (
                    cellI < 0
                 || cellI >= label(v1Cells.size())
                 || cellI >= label(exactVolumeParent.size())
                 || cellI >= label(exactVolumeLocalChild.size())
                 || cellI >= label(exactVolumeRefType.size())
                )
                {
                    return;
                }

                if( exactVolumeRefType[cellI] != 1 )
                    return;

                const label parentI =
                    exactVolumeParent[cellI];

                const label localChildI =
                    exactVolumeLocalChild[cellI];

                if
                (
                    parentI < 0
                 || localChildI < 0
                )
                    return;

                point cellCtr(vector::zero);
                scalar cellVol = -GREAT;

                if
                (
                    !v1OFCellCentreVolume
                    (
                        cellI,
                        cellCtr,
                        cellVol
                    )
                )
                    return;

                const cell& c =
                    v1Cells[cellI];

                label nBadThisCell = 0;

                forAll(c, cfI)
                {
                    const label faceI =
                        c[cfI];

                    if
                    (
                        faceI < 0
                     || faceI >= label(v1Faces.size())
                     || faceI >= label(exactOwner.size())
                     || faceI >= label(v1Neighbour.size())
                    )
                        continue;

                    const scalar pyrVol =
                        pyramidPointFaceRef
                        (
                            v1Faces[faceI],
                            cellCtr
                        ).mag(v1Points);

                    bool badForThisCell = false;
                    bool thisCellIsOwner = false;

                    if( exactOwner[faceI] == cellI )
                    {
                        thisCellIsOwner = true;

                        if( pyrVol > SMALL )
                            badForThisCell = true;
                    }
                    else if( v1Neighbour[faceI] == cellI )
                    {
                        if( pyrVol < -SMALL )
                            badForThisCell = true;
                    }
                    else
                    {
                        continue;
                    }

                    if( !badForThisCell )
                        continue;

                    ++nBadThisCell;
                    ++blV27WallChildBadFaces;

                    word role("UNKNOWN");

                    // First identify the actual wall/base face from
                    // topology: every point on it is a known hair root.
                    bool allHairStarts = true;

                    const face& f =
                        v1Faces[faceI];

                    forAll(f, fpI)
                    {
                        if
                        (
                            !v1HairStartPoints.found
                            (
                                f[fpI]
                            )
                        )
                        {
                            allHairStarts = false;
                            break;
                        }
                    }

                    if( allHairStarts )
                    {
                        role = "WALL_BASE";
                        ++blV27WallChildWallBaseBad;
                    }
                    else
                    {
                        label otherCellI = -1;

                        if( thisCellIsOwner )
                        {
                            otherCellI =
                                v1Neighbour[faceI];
                        }
                        else
                        {
                            otherCellI =
                                exactOwner[faceI];
                        }

                        bool classifiedInternalOuter =
                            false;

                        if
                        (
                            otherCellI >= 0
                         && otherCellI <
                            label(exactVolumeParent.size())
                         && otherCellI <
                            label(exactVolumeLocalChild.size())
                         && otherCellI <
                            label(exactVolumeRefType.size())
                         && exactVolumeRefType[otherCellI] == 1
                         && exactVolumeParent[otherCellI] ==
                            parentI
                         && exactVolumeLocalChild[otherCellI] ==
                            localChildI - 1
                        )
                        {
                            role = "INTERNAL_OUTER";
                            classifiedInternalOuter = true;
                            ++blV27WallChildInternalOuterBad;
                        }

                        if( !classifiedInternalOuter )
                        {
                            // For the first wall-adjacent child, every
                            // non-wall face other than its coreward sibling
                            // interface is a lateral face unless addressing
                            // says that another same-parent relationship
                            // exists unexpectedly.
                            if
                            (
                                otherCellI >= 0
                             && otherCellI <
                                label(exactVolumeParent.size())
                             && exactVolumeParent[otherCellI] ==
                                parentI
                            )
                            {
                                role = "UNKNOWN";
                                ++blV27WallChildUnknownBad;
                            }
                            else
                            {
                                role = "LATERAL";
                                ++blV27WallChildLateralBad;
                            }
                        }
                    }

                    Info
                        << "CFMITCH V2.7 WALLCHILD BADFACE:"
                        << " cell=" << cellI
                        << " parent=" << parentI
                        << " bfI="
                        << cellToBaseBndFace_[cellI]
                        << " localChild="
                        << localChildI
                        << " face=" << faceI
                        << " role=" << role
                        << " side="
                        << (
                            thisCellIsOwner
                          ? "owner"
                          : "neighbour"
                           )
                        << " pyrVol=" << pyrVol
                        << " nPts=" << f.size()
                        << endl;

                    if( role == "WALL_BASE" )
                    {
                        Info
                            << "CFMITCH V2.7.2 WALLBASE ROOTS:"
                            << " cell=" << cellI
                            << " parent=" << parentI
                            << " bfI="
                            << cellToBaseBndFace_[cellI]
                            << " face=" << faceI
                            << " roots=";

                        forAll(f, fpI)
                        {
                            if( fpI )
                                Info << ',';

                            Info << f[fpI];
                        }

                        Info << endl;
                    }
                }

                if( nBadThisCell > 0 )
                {
                    ++blV27WallChildCells;

                    if( nBadThisCell > 1 )
                        ++blV27WallChildMultiBad;

                    Info
                        << "CFMITCH V2.7 WALLCHILD CELL:"
                        << " cell=" << cellI
                        << " parent=" << parentI
                        << " bfI="
                        << cellToBaseBndFace_[cellI]
                        << " localChild="
                        << localChildI
                        << " badFaces="
                        << nBadThisCell
                        << endl;
                }
            };


            // ----------------------------------------------------------
            // Generalized repair seed population.
            //
            // Start with the historical negative-volume seeds.
            // ----------------------------------------------------------

            boolList v1RepairCellMask
            (
                v1Cells.size(),
                false
            );

            forAll(blV1NegativeType1Cells, negI)
            {
                const label cellI =
                    blV1NegativeType1Cells[negI];

                if
                (
                    cellI >= 0
                 && cellI < label(v1RepairCellMask.size())
                )
                {
                    v1RepairCellMask[cellI] = true;
                }
            }


            // ----------------------------------------------------------
            // Add positive/otherwise non-negative pyramid-defect cells,
            // but ONLY when the cell actually contains an unambiguous
            // generated split-edge point.  This establishes that V1C has
            // a legitimate construction degree of freedom for the cell.
            // ----------------------------------------------------------

            label nV25GeneratedEdgeCells = 0;
            label nV25PyramidEvalFailed = 0;

            forAll(v1Cells, cellI)
            {
                const cell& c =
                    v1Cells[cellI];

                bool hasGeneratedSplitPoint = false;

                forAll(c, cfI)
                {
                    const label faceI =
                        c[cfI];

                    if
                    (
                        faceI < 0
                     || faceI >= label(v1Faces.size())
                    )
                        continue;

                    const face& f =
                        v1Faces[faceI];

                    forAll(f, fpI)
                    {
                        if
                        (
                            v1PointEdge.find(f[fpI])
                         != v1PointEdge.end()
                        )
                        {
                            hasGeneratedSplitPoint = true;
                            break;
                        }
                    }

                    if( hasGeneratedSplitPoint )
                        break;
                }

                if( !hasGeneratedSplitPoint )
                    continue;

                ++nV25GeneratedEdgeCells;

                const label nBadPyr =
                    v1CellBadPyramidCount(cellI);

                if( nBadPyr == labelMax )
                {
                    ++nV25PyramidEvalFailed;
                    continue;
                }

                if
                (
                    nBadPyr > 0
                 && !v1RepairCellMask[cellI]
                )
                {
                    v1RepairCellMask[cellI] = true;
                    ++blV25PyramidAdditionalSeeds;
                }
            }


            // Deterministic ascending-cell repair order.
            DynList<label, 256> blV1RepairCells;

            forAll(v1RepairCellMask, cellI)
            {
                if( v1RepairCellMask[cellI] )
                    blV1RepairCells.append(cellI);
            }

            blV25RepairSeeds =
                blV1RepairCells.size();


            Info
                << "CFMITCH V2.5 BIRTH QUALITY SEEDS:"
                << " negative="
                << blV1InitialNegative
                << " pyramidAdditional="
                << blV25PyramidAdditionalSeeds
                << " totalRepairCells="
                << blV25RepairSeeds
                << " generatedEdgeCells="
                << nV25GeneratedEdgeCells
                << " pyramidEvalFailed="
                << nV25PyramidEvalFailed
                << endl;


            // ----------------------------------------------------------
            // Candidate split edges are those containing an interior
            // generated point belonging to a V1C repair cell.
            // ----------------------------------------------------------

            std::set<label> v1CandidateEdges;

            forAll(blV1RepairCells, badI)
            {
                const label badCellI =
                    blV1RepairCells[badI];

                if
                (
                    badCellI < 0
                 || badCellI >= label(v1Cells.size())
                )
                    continue;

                const cell& c =
                    v1Cells[badCellI];

                forAll(c, cfI)
                {
                    const face& f =
                        v1Faces[c[cfI]];

                    forAll(f, fpI)
                    {
                        std::map<label,label>::const_iterator
                            it =
                                v1PointEdge.find(f[fpI]);

                        if( it != v1PointEdge.end() )
                            v1CandidateEdges.insert(it->second);
                    }
                }
            }


            // ----------------------------------------------------------
            // Build complete incident-cell stars for all candidate
            // split-edge chains with ONE global cell scan.
            // ----------------------------------------------------------

            std::map<label,label> v1PointCandidateEdge;

            std::map<label, std::set<label> >
                v1EdgeCells;

            for
            (
                std::set<label>::const_iterator
                    eIt=v1CandidateEdges.begin();
                eIt!=v1CandidateEdges.end();
                ++eIt
            )
            {
                const label seI = *eIt;

                v1EdgeCells[seI];

                const label rowSize =
                    newVerticesForSplitEdge_.
                        sizeOfRow(seI);

                for
                (
                    label rowI=1;
                    rowI<rowSize-1;
                    ++rowI
                )
                {
                    const label pointI =
                        newVerticesForSplitEdge_
                        (
                            seI,
                            rowI
                        );

                    v1PointCandidateEdge[pointI] =
                        seI;
                }
            }


            forAll(v1Cells, cellI)
            {
                std::set<label> touchedEdges;

                const cell& c =
                    v1Cells[cellI];

                forAll(c, cfI)
                {
                    const face& f =
                        v1Faces[c[cfI]];

                    forAll(f, fpI)
                    {
                        std::map<label,label>::const_iterator
                            pIt =
                                v1PointCandidateEdge.find
                                (
                                    f[fpI]
                                );

                        if
                        (
                            pIt
                         != v1PointCandidateEdge.end()
                        )
                            touchedEdges.insert
                            (
                                pIt->second
                            );
                    }
                }

                for
                (
                    std::set<label>::const_iterator
                        eIt=touchedEdges.begin();
                    eIt!=touchedEdges.end();
                    ++eIt
                )
                {
                    v1EdgeCells[*eIt].
                        insert(cellI);
                }
            }


            // Search smallest amplitudes first.
            static const scalar v1AmplitudeMag[] =
            {
                scalar(0.01),
                scalar(0.02),
                scalar(0.03),
                scalar(0.04),
                scalar(0.05),
                scalar(0.06),
                scalar(0.08),
                scalar(0.10),
                scalar(0.12),
                scalar(0.16),
                scalar(0.18),
                scalar(0.20),
                scalar(0.25),
                scalar(0.30),
                scalar(0.35),
                scalar(0.40),
                scalar(0.425),
                scalar(0.45),
                scalar(0.475),
                scalar(0.50)
            };

            static const label nV1AmplitudeMag =
                sizeof(v1AmplitudeMag)
               /sizeof(v1AmplitudeMag[0]);


            forAll(blV1RepairCells, badI)
            {
                const label badCellI =
                    blV1RepairCells[badI];

                scalar badVolBefore =
                    v1CellVolume(badCellI);

                label badPyrBefore = 0;

                if( badVolBefore > scalar(0) )
                {
                    badPyrBefore =
                        v1CellBadPyramidCount(badCellI);

                    if( badPyrBefore == labelMax )
                    {
                        ++blV1Unresolved;

                        Info
                            << "CFMITCH V2.5 BIRTH TARGET SKIP:"
                            << " cell=" << badCellI
                            << " reason=pyramidEvaluationFailed"
                            << endl;

                        continue;
                    }
                }

                // A previously committed shared-chain transaction may have
                // fixed this positive-volume pyramid seed already.
                if
                (
                    badVolBefore > scalar(0)
                 && badPyrBefore == 0
                )
                {
                    ++blV1Fixed;
                    continue;
                }

                const bool requirePyramidRepair =
                    (
                        badVolBefore > scalar(0)
                     && badPyrBefore > 0
                    );


                // Candidate split edges touching this bad child.
                std::set<label> badEdges;

                const cell& badCell =
                    v1Cells[badCellI];

                forAll(badCell, cfI)
                {
                    const face& f =
                        v1Faces[badCell[cfI]];

                    forAll(f, fpI)
                    {
                        std::map<label,label>::const_iterator it =
                            v1PointEdge.find(f[fpI]);

                        if( it != v1PointEdge.end() )
                            badEdges.insert(it->second);
                    }
                }


                bool found = false;

                label bestEdge = -1;
                scalar bestAmplitude = scalar(0);
                scalar bestTargetVol = -GREAT;
                scalar bestMinRatio = -GREAT;

                List<point> bestPositions;

                scalar bestTargetOFVol = -GREAT;
                scalar bestMaxOFSkew = GREAT;
                scalar bestMinOFOrtho = -GREAT;
                scalar bestMinOFRatio = -GREAT;
                scalar bestMaxMove = GREAT;
                label bestBadPyrCount = labelMax;

                // v1d fallback bookkeeping.  Stage 1 remains the existing
                // exhaustive single-hair v1c search.  These fields are used
                // only if Stage 1 finds no safe candidate and a face-derived
                // two-hair transaction succeeds.
                bool bestIsPair = false;

                label bestEdgeA = -1;
                label bestEdgeB = -1;

                scalar bestAmplitudeA = scalar(0);
                scalar bestAmplitudeB = scalar(0);

                List<point> bestPositionsA;
                List<point> bestPositionsB;


                for
                (
                    std::set<label>::const_iterator
                        eIt=badEdges.begin();
                    eIt!=badEdges.end();
                    ++eIt
                )
                {
                    const label seI = *eIt;

                    const edge& se =
                        splitEdges_[seI];

                    const point edgeStart =
                        v1Points[se.start()];

                    const point edgeEnd =
                        v1Points[se.end()];

                    const vector edgeVec =
                        edgeEnd-edgeStart;

                    const scalar edgeMagSqr =
                        magSqr(edgeVec);

                    if( edgeMagSqr <= VSMALL )
                        continue;


                    const label rowSize =
                        newVerticesForSplitEdge_.
                            sizeOfRow(seI);

                    if( rowSize < 3 )
                        continue;


                    List<point> originalPositions
                    (
                        rowSize-2
                    );

                    List<scalar> originalT
                    (
                        rowSize
                    );

                    originalT[0] = scalar(0);
                    originalT[rowSize-1] = scalar(1);

                    bool chainValid = true;

                    for
                    (
                        label rowI=1;
                        rowI<rowSize-1;
                        ++rowI
                    )
                    {
                        const label pointI =
                            newVerticesForSplitEdge_
                            (
                                seI,
                                rowI
                            );

                        originalPositions[rowI-1] =
                            v1Points[pointI];

                        originalT[rowI] =
                            (
                                (v1Points[pointI]-edgeStart)
                              & edgeVec
                            )
                           /(edgeMagSqr + VSMALL);

                        if
                        (
                            !(originalT[rowI] > originalT[rowI-1])
                         || !(originalT[rowI] < scalar(1))
                        )
                        {
                            chainValid = false;
                            break;
                        }
                    }

                    if( !chainValid )
                        continue;


                    std::map
                    <
                        label,
                        std::set<label>
                    >::const_iterator starIt =
                        v1EdgeCells.find(seI);

                    if
                    (
                        starIt == v1EdgeCells.end()
                     || starIt->second.empty()
                    )
                        continue;


                    // Baseline star volumes.
                    std::map<label,scalar> baselineVolumes;

                    for
                    (
                        std::set<label>::const_iterator
                            cIt=starIt->second.begin();
                        cIt!=starIt->second.end();
                        ++cIt
                    )
                    {
                        baselineVolumes[*cIt] =
                            v1CellVolume(*cIt);
                    }


                    // OpenFOAM-parity baseline cell volumes.
                    std::map<label,scalar> baselineOFVolumes;

                    for
                    (
                        std::set<label>::const_iterator
                            cIt=starIt->second.begin();
                        cIt!=starIt->second.end();
                        ++cIt
                    )
                    {
                        point cc(vector::zero);
                        scalar cv = -GREAT;

                        if
                        (
                            !v1OFCellCentreVolume
                            (
                                *cIt,
                                cc,
                                cv
                            )
                        )
                        {
                            chainValid = false;
                            break;
                        }

                        baselineOFVolumes[*cIt] = cv;
                    }

                    if( !chainValid )
                        continue;


                    // Every face whose quality can change because this
                    // candidate hair moves either a face vertex or an
                    // incident cell centre.
                    std::set<label> affectedFaces;

                    for
                    (
                        std::set<label>::const_iterator
                            cIt=starIt->second.begin();
                        cIt!=starIt->second.end();
                        ++cIt
                    )
                    {
                        const cell& starCell =
                            v1Cells[*cIt];

                        forAll(starCell, cfI)
                            affectedFaces.insert(starCell[cfI]);
                    }


                    std::map<label,scalar> baselineOFSkew;
                    std::map<label,scalar> baselineOFOrtho;
                    std::map<label,scalar> baselineOFPyrMargin;
                    std::map<label,bool> baselineOFPyrBad;
                    std::map<label,bool> baselineOFEligible;

                    bool baselineQualityValid = true;

                    for
                    (
                        std::set<label>::const_iterator
                            fIt=affectedFaces.begin();
                        fIt!=affectedFaces.end();
                        ++fIt
                    )
                    {
                        const label faceI = *fIt;

                        scalar skew = GREAT;
                        scalar ortho = -GREAT;
                        scalar pyrMargin = -GREAT;
                        scalar ownVol = -GREAT;
                        scalar neiVol = -GREAT;
                        bool pyrBad = true;

                        if
                        (
                            !v1OFFaceQuality
                            (
                                faceI,
                                skew,
                                ortho,
                                pyrMargin,
                                pyrBad,
                                ownVol,
                                neiVol
                            )
                        )
                        {
                            baselineQualityValid = false;
                            break;
                        }

                        baselineOFSkew[faceI] =
                            skew;

                        baselineOFOrtho[faceI] =
                            ortho;

                        baselineOFPyrMargin[faceI] =
                            pyrMargin;

                        baselineOFPyrBad[faceI] =
                            pyrBad;

                        baselineOFEligible[faceI] =
                            (
                                ownVol > scalar(0)
                             && (
                                    v1Neighbour[faceI] < 0
                                 || neiVol > scalar(0)
                                )
                            );
                    }

                    if( !baselineQualityValid )
                        continue;


                    auto restoreChain =
                    [&]()
                    {
                        for
                        (
                            label rowI=1;
                            rowI<rowSize-1;
                            ++rowI
                        )
                        {
                            const label pointI =
                                newVerticesForSplitEdge_
                                (
                                    seI,
                                    rowI
                                );

                            v1Points[pointI] =
                                originalPositions[rowI-1];
                        }
                    };


                    auto applyAmplitude =
                    [&]
                    (
                        const scalar amplitude
                    ) -> bool
                    {
                        // BLValidityRepair v1b:
                        //
                        // amplitude is the RELATIVE correction applied to
                        // the first generated row.  That absolute shift is
                        // then monotonically tapered to zero towards the
                        // BL/core endpoint.
                        //
                        // Unlike v1a, no deeper point can move farther than
                        // the first generated row.
                        const scalar firstT =
                            originalT[1];

                        if
                        (
                            firstT <= scalar(0)
                         || firstT >= scalar(1)
                        )
                        {
                            restoreChain();
                            return false;
                        }

                        scalar prevT = scalar(0);

                        for
                        (
                            label rowI=1;
                            rowI<rowSize-1;
                            ++rowI
                        )
                        {
                            const scalar t =
                                originalT[rowI];

                            const scalar decayBase =
                                Foam::max
                                (
                                    scalar(0),
                                    (scalar(1)-t)
                                   /(scalar(1)-firstT + VSMALL)
                                );

                            const scalar deltaT =
                                amplitude
                               *firstT
                               *Foam::pow
                                (
                                    decayBase,
                                    scalar(4)
                                );

                            const scalar warpedT =
                                t + deltaT;

                            if
                            (
                                !(warpedT > prevT)
                             || !(warpedT < originalT[rowI+1])
                            )
                            {
                                restoreChain();
                                return false;
                            }

                            const label pointI =
                                newVerticesForSplitEdge_
                                (
                                    seI,
                                    rowI
                                );

                            v1Points[pointI] =
                                edgeStart
                              + warpedT*edgeVec;

                            prevT = warpedT;
                        }

                        return true;
                    };


                    for
                    (
                        label ampI=0;
                        ampI<nV1AmplitudeMag;
                        ++ampI
                    )
                    {
                        for(label sign=-1; sign<=1; sign+=2)
                        {
                            restoreChain();

                            const scalar amplitude =
                                scalar(sign)
                               *v1AmplitudeMag[ampI];

                            if
                            (
                                !applyAmplitude(amplitude)
                            )
                                continue;


                            const scalar targetVol =
                                v1CellVolume(badCellI);

                            // Require a meaningful positive margin, not
                            // merely a numerical sign crossing.
                            const scalar targetFloor =
                                scalar(0.20)
                               *Foam::mag(badVolBefore);

                            if( targetVol <= targetFloor )
                                continue;


                            bool safe = true;
                            scalar minPositiveRatio = GREAT;

                            for
                            (
                                std::set<label>::const_iterator
                                    cIt=starIt->second.begin();
                                cIt!=starIt->second.end();
                                ++cIt
                            )
                            {
                                const label starCellI =
                                    *cIt;

                                const scalar oldV =
                                    baselineVolumes[starCellI];

                                const scalar newV =
                                    v1CellVolume(starCellI);


                                if( oldV > scalar(0) )
                                {
                                    // No neighbour is allowed to collapse
                                    // while fixing the target.
                                    if( newV <= scalar(0) )
                                    {
                                        safe = false;
                                        break;
                                    }

                                    const scalar ratio =
                                        newV
                                       /(oldV + VSMALL);

                                    minPositiveRatio =
                                        Foam::min
                                        (
                                            minPositiveRatio,
                                            ratio
                                        );

                                    if( ratio < scalar(0.75) )
                                    {
                                        safe = false;
                                        break;
                                    }
                                }
                                else
                                {
                                    // Existing negative neighbours may
                                    // remain for their own later repair,
                                    // but this chain may not worsen them.
                                    if( newV < oldV )
                                    {
                                        safe = false;
                                        break;
                                    }
                                }
                            }


                            if( !safe )
                                continue;


                            // ------------------------------------------
                            // OpenFOAM-parity signed-volume star.
                            // ------------------------------------------

                            point targetOFCentre(vector::zero);
                            scalar targetOFVol = -GREAT;

                            if
                            (
                                !v1OFCellCentreVolume
                                (
                                    badCellI,
                                    targetOFCentre,
                                    targetOFVol
                                )
                            )
                                continue;

                            if( targetOFVol <= targetFloor )
                                continue;


                            scalar minOFPositiveRatio = GREAT;

                            for
                            (
                                std::set<label>::const_iterator
                                    cIt=starIt->second.begin();
                                cIt!=starIt->second.end();
                                ++cIt
                            )
                            {
                                const label starCellI =
                                    *cIt;

                                const scalar oldOFV =
                                    baselineOFVolumes[starCellI];

                                point newOFCentre(vector::zero);
                                scalar newOFV = -GREAT;

                                if
                                (
                                    !v1OFCellCentreVolume
                                    (
                                        starCellI,
                                        newOFCentre,
                                        newOFV
                                    )
                                )
                                {
                                    safe = false;
                                    break;
                                }

                                if( oldOFV > scalar(0) )
                                {
                                    if( newOFV <= scalar(0) )
                                    {
                                        safe = false;
                                        break;
                                    }

                                    const scalar ratio =
                                        newOFV
                                       /(oldOFV + VSMALL);

                                    minOFPositiveRatio =
                                        Foam::min
                                        (
                                            minOFPositiveRatio,
                                            ratio
                                        );

                                    if( ratio < scalar(0.75) )
                                    {
                                        safe = false;
                                        break;
                                    }
                                }
                                else if( newOFV < oldOFV )
                                {
                                    safe = false;
                                    break;
                                }
                            }

                            if( !safe )
                                continue;


                            // ------------------------------------------
                            // Solver-facing local quality.
                            //
                            // A face whose baseline owner/neighbour cell
                            // was already non-positive does not have a
                            // trustworthy baseline FV metric.  Once the
                            // trial makes that face fully positive, require
                            // absolute safe quality.
                            //
                            // For an already-valid baseline face, permit
                            // existing bad quality to remain only when it
                            // does not regress.
                            // ------------------------------------------

                            static const scalar v1CSkewLimit =
                                scalar(4);

                            // Hard FV admissibility boundary:
                            // orthogonality < 0 corresponds to >90 deg.
                            // Faces in the 70-90 deg range remain eligible
                            // but are penalised by best-candidate ranking.
                            static const scalar v1COrtho90 =
                                scalar(0);

                            bool ofSafe = true;
                            const char* rejectReason = "none";

                            scalar trialMaxSkew = scalar(0);
                            scalar trialMinOrtho = GREAT;
                            label trialBadPyrCount = 0;

                            for
                            (
                                std::set<label>::const_iterator
                                    fIt=affectedFaces.begin();
                                fIt!=affectedFaces.end();
                                ++fIt
                            )
                            {
                                const label faceI = *fIt;

                                scalar trialSkew = GREAT;
                                scalar trialOrtho = -GREAT;
                                scalar trialPyrMargin = -GREAT;
                                scalar trialOwnVol = -GREAT;
                                scalar trialNeiVol = -GREAT;
                                bool trialPyrBad = true;

                                if
                                (
                                    !v1OFFaceQuality
                                    (
                                        faceI,
                                        trialSkew,
                                        trialOrtho,
                                        trialPyrMargin,
                                        trialPyrBad,
                                        trialOwnVol,
                                        trialNeiVol
                                    )
                                )
                                {
                                    ofSafe = false;
                                    rejectReason = "ofGeometry";
                                    break;
                                }

                                const bool trialEligible =
                                    (
                                        trialOwnVol > scalar(0)
                                     && (
                                            v1Neighbour[faceI] < 0
                                         || trialNeiVol > scalar(0)
                                        )
                                    );

                                // Faces still attached to an existing
                                // negative cell are handled when that cell
                                // receives its own repair.  Their skew and
                                // orthogonality are not meaningful yet.
                                if( !trialEligible )
                                    continue;

                                if( trialPyrBad )
                                    ++trialBadPyrCount;

                                trialMaxSkew =
                                    Foam::max
                                    (
                                        trialMaxSkew,
                                        trialSkew
                                    );

                                if( v1Neighbour[faceI] >= 0 )
                                {
                                    trialMinOrtho =
                                        Foam::min
                                        (
                                            trialMinOrtho,
                                            trialOrtho
                                        );
                                }

                                const bool baselineEligible =
                                    baselineOFEligible[faceI];

                                if( !baselineEligible )
                                {
                                    // The candidate has converted a face
                                    // next to an invalid cell into fully
                                    // positive FV geometry.  It must now
                                    // satisfy absolute production safety.
                                    if( trialPyrBad )
                                    {
                                        Info
                                            << "BL_VALIDITY_REPAIR_V1C_PYRFAIL"
                                            << " cell=" << badCellI
                                            << " parent="
                                            << exactVolumeParent[badCellI]
                                            << " localChild="
                                            << exactVolumeLocalChild[badCellI]
                                            << " splitEdge=" << seI
                                            << " amplitude=" << amplitude
                                            << " face=" << faceI
                                            << " owner="
                                            << exactOwner[faceI]
                                            << " neighbour="
                                            << v1Neighbour[faceI]
                                            << " baselineEligible="
                                            << baselineEligible
                                            << " baselinePyrBad="
                                            << baselineOFPyrBad[faceI]
                                            << " baselinePyrMargin="
                                            << baselineOFPyrMargin[faceI]
                                            << " trialPyrMargin="
                                            << trialPyrMargin
                                            << " trialSkew="
                                            << trialSkew
                                            << " trialOrtho="
                                            << trialOrtho
                                            << " facePoints="
                                            << v1Faces[faceI]
                                            << " pointEdges=(";

                                        const face& diagFace =
                                            v1Faces[faceI];

                                        forAll(diagFace, dpi)
                                        {
                                            const label pointI =
                                                diagFace[dpi];

                                            std::map<label,label>::
                                                const_iterator peIt =
                                                    v1PointEdge.find
                                                    (
                                                        pointI
                                                    );

                                            Info
                                                << pointI
                                                << ":";

                                            if
                                            (
                                                peIt
                                             != v1PointEdge.end()
                                            )
                                            {
                                                Info << peIt->second;
                                            }
                                            else
                                            {
                                                Info << -1;
                                            }

                                            if
                                            (
                                                dpi
                                              < diagFace.size()-1
                                            )
                                                Info << ",";
                                        }

                                        Info << ")" << endl;

                                        ofSafe = false;
                                        rejectReason = "newPyramid";
                                        break;
                                    }

                                    if( trialSkew > v1CSkewLimit )
                                    {
                                        ofSafe = false;
                                        rejectReason = "newSkew";
                                        break;
                                    }

                                    if
                                    (
                                        v1Neighbour[faceI] >= 0
                                     && trialOrtho < v1COrtho90
                                    )
                                    {
                                        ofSafe = false;
                                        rejectReason = "newNonOrtho90";
                                        break;
                                    }

                                    continue;
                                }


                                const bool basePyrBad =
                                    baselineOFPyrBad[faceI];

                                const scalar basePyrMargin =
                                    baselineOFPyrMargin[faceI];

                                if
                                (
                                    !basePyrBad
                                 && trialPyrBad
                                )
                                {
                                    Info
                                        << "BL_VALIDITY_REPAIR_V1C_PYRFAIL"
                                        << " cell=" << badCellI
                                        << " parent="
                                        << exactVolumeParent[badCellI]
                                        << " localChild="
                                        << exactVolumeLocalChild[badCellI]
                                        << " splitEdge=" << seI
                                        << " amplitude=" << amplitude
                                        << " face=" << faceI
                                        << " owner="
                                        << exactOwner[faceI]
                                        << " neighbour="
                                        << v1Neighbour[faceI]
                                        << " baselineEligible="
                                        << baselineEligible
                                        << " baselinePyrBad="
                                        << basePyrBad
                                        << " baselinePyrMargin="
                                        << basePyrMargin
                                        << " trialPyrMargin="
                                        << trialPyrMargin
                                        << " trialSkew="
                                        << trialSkew
                                        << " trialOrtho="
                                        << trialOrtho
                                        << " facePoints="
                                        << v1Faces[faceI]
                                        << " pointEdges=(";

                                    const face& diagFace =
                                        v1Faces[faceI];

                                    forAll(diagFace, dpi)
                                    {
                                        const label pointI =
                                            diagFace[dpi];

                                        std::map<label,label>::
                                            const_iterator peIt =
                                                v1PointEdge.find
                                                (
                                                    pointI
                                                );

                                        Info
                                            << pointI
                                            << ":";

                                        if
                                        (
                                            peIt
                                         != v1PointEdge.end()
                                        )
                                        {
                                            Info << peIt->second;
                                        }
                                        else
                                        {
                                            Info << -1;
                                        }

                                        if
                                        (
                                            dpi
                                          < diagFace.size()-1
                                        )
                                            Info << ",";
                                    }

                                    Info << ")" << endl;

                                    ofSafe = false;
                                    rejectReason = "newPyramid";
                                    break;
                                }

                                if
                                (
                                    basePyrBad
                                 && trialPyrBad
                                )
                                {
                                    const scalar pyrTol =
                                        scalar(1e-12)
                                       *(
                                            Foam::mag(basePyrMargin)
                                          + SMALL
                                        );

                                    if
                                    (
                                        trialPyrMargin
                                      < basePyrMargin-pyrTol
                                    )
                                    {
                                        // Diagnostic-only:
                                        // expose the EXISTING bad pyramid
                                        // which blocks this candidate.
                                        //
                                        // Deduplicate by:
                                        //
                                        //     target cell
                                        //     moved split edge
                                        //     blocking face
                                        //
                                        // so the amplitude sweep does not
                                        // emit thousands of equivalent lines.
                                        const std::pair
                                        <
                                            label,
                                            std::pair<label,label>
                                        > diagKey =
                                            std::make_pair
                                            (
                                                badCellI,
                                                std::make_pair
                                                (
                                                    seI,
                                                    faceI
                                                )
                                            );

                                        if
                                        (
                                            blV1WorsePyrDiagSeen.insert
                                            (
                                                diagKey
                                            ).second
                                        )
                                        {
                                            const face& diagFace =
                                                v1Faces[faceI];

                                            bool movedEdgeOnFace =
                                                false;

                                            std::set<label>
                                                generatedFaceEdges;

                                            forAll(diagFace, dpi)
                                            {
                                                const label pointI =
                                                    diagFace[dpi];

                                                std::map<label,label>::
                                                    const_iterator peIt =
                                                        v1PointEdge.find
                                                        (
                                                            pointI
                                                        );

                                                if
                                                (
                                                    peIt
                                                 != v1PointEdge.end()
                                                )
                                                {
                                                    generatedFaceEdges.insert
                                                    (
                                                        peIt->second
                                                    );

                                                    if
                                                    (
                                                        peIt->second
                                                     == seI
                                                    )
                                                    {
                                                        movedEdgeOnFace =
                                                            true;
                                                    }
                                                }
                                            }

                                            const bool touchesTarget =
                                                (
                                                    exactOwner[faceI]
                                                 == badCellI
                                                 || v1Neighbour[faceI]
                                                 == badCellI
                                                );

                                            // CFMitch production-quiet:
                                            // per-trial PYRWORSE forensic
                                            // output suppressed.  Quality
                                            // rejection behavior unchanged.
                                            (void)touchesTarget;
                                            (void)movedEdgeOnFace;
                                            (void)generatedFaceEdges;
                                        }

                                        ofSafe = false;
                                        rejectReason = "worsePyramid";
                                        break;
                                    }
                                }


                                const scalar baseSkew =
                                    baselineOFSkew[faceI];

                                const scalar skewTol =
                                    scalar(1e-10)
                                   *(
                                        scalar(1)
                                      + Foam::mag(baseSkew)
                                    );

                                if( baseSkew <= v1CSkewLimit )
                                {
                                    if( trialSkew > v1CSkewLimit )
                                    {
                                        ofSafe = false;
                                        rejectReason = "newSkew";
                                        break;
                                    }
                                }
                                else if
                                (
                                    trialSkew
                                  > baseSkew+skewTol
                                )
                                {
                                    ofSafe = false;
                                    rejectReason = "worseSkew";
                                    break;
                                }


                                if( v1Neighbour[faceI] >= 0 )
                                {
                                    const scalar baseOrtho =
                                        baselineOFOrtho[faceI];

                                    const scalar orthoTol =
                                        scalar(1e-12);

                                    if( baseOrtho >= v1COrtho90 )
                                    {
                                        if
                                        (
                                            trialOrtho < v1COrtho90
                                        )
                                        {
                                            ofSafe = false;
                                            rejectReason =
                                                "newNonOrtho90";
                                            break;
                                        }
                                    }
                                    else if
                                    (
                                        trialOrtho
                                      < baseOrtho-orthoTol
                                    )
                                    {
                                        ofSafe = false;
                                        rejectReason =
                                            "worseNonOrtho";
                                        break;
                                    }
                                }
                            }


                            if( !ofSafe )
                            {
                                // CFMitch production-quiet:
                                // per-trial reject output suppressed.
                                // rejectReason and all acceptance logic
                                // remain unchanged.

                                continue;
                            }


                            // CFMitch V2.5 target-pyramid admissibility.
                            //
                            // Existing V1C permits a pre-existing bad
                            // pyramid to remain when it does not worsen.
                            // That is appropriate while rescuing a negative
                            // cell, but it is insufficient when the pyramid
                            // itself is the reason this positive cell was
                            // seeded.
                            //
                            // For a quality-only target, require this trial
                            // to eliminate every bad pyramid side belonging
                            // to the target cell before it may enter the
                            // candidate ranking.
                            if( requirePyramidRepair )
                            {
                                const label trialTargetBadPyr =
                                    v1CellBadPyramidCount(badCellI);

                                if( trialTargetBadPyr != 0 )
                                {
                                    ++blV25TargetPyrRejects;
                                    continue;
                                }
                            }


                            scalar trialMaxMove = scalar(0);

                            for
                            (
                                label rowI=1;
                                rowI<rowSize-1;
                                ++rowI
                            )
                            {
                                const label pointI =
                                    newVerticesForSplitEdge_
                                    (
                                        seI,
                                        rowI
                                    );

                                trialMaxMove =
                                    Foam::max
                                    (
                                        trialMaxMove,
                                        mag
                                        (
                                            v1Points[pointI]
                                          - originalPositions[rowI-1]
                                        )
                                    );
                            }


                            Info
                                << "BL_VALIDITY_REPAIR_V1C_CANDIDATE"
                                << " cell=" << badCellI
                                << " splitEdge=" << seI
                                << " amplitude=" << amplitude
                                << " targetVol=" << targetVol
                                << " targetOFVol=" << targetOFVol
                                << " maxOFSkew=" << trialMaxSkew
                                << " minOFOrtho=" << trialMinOrtho
                                << " badOFPyr="
                                << trialBadPyrCount
                                << " maxPhysicalMove="
                                << trialMaxMove
                                << " minPositiveStarRatio="
                                << minPositiveRatio
                                << " minOFPositiveStarRatio="
                                << minOFPositiveRatio
                                << endl;


                            // Deterministic lexicographic quality ranking:
                            //
                            //   1) fewer bad OF pyramids
                            //   2) lower maximum OF skew
                            //   3) higher minimum OF orthogonality
                            //   4) smaller physical movement
                            //   5) better raw-volume star preservation
                            bool better = !found;

                            if( found )
                            {
                                if
                                (
                                    trialBadPyrCount
                                  < bestBadPyrCount
                                )
                                {
                                    better = true;
                                }
                                else if
                                (
                                    trialBadPyrCount
                                 == bestBadPyrCount
                                )
                                {
                                    if
                                    (
                                        trialMaxSkew
                                      < bestMaxOFSkew
                                       - scalar(1e-12)
                                    )
                                    {
                                        better = true;
                                    }
                                    else if
                                    (
                                        Foam::mag
                                        (
                                            trialMaxSkew
                                          - bestMaxOFSkew
                                        )
                                      <= scalar(1e-12)
                                    )
                                    {
                                        if
                                        (
                                            trialMinOrtho
                                          > bestMinOFOrtho
                                           + scalar(1e-12)
                                        )
                                        {
                                            better = true;
                                        }
                                        else if
                                        (
                                            Foam::mag
                                            (
                                                trialMinOrtho
                                              - bestMinOFOrtho
                                            )
                                          <= scalar(1e-12)
                                        )
                                        {
                                            if
                                            (
                                                trialMaxMove
                                              < bestMaxMove
                                               - scalar(1e-15)
                                            )
                                            {
                                                better = true;
                                            }
                                            else if
                                            (
                                                Foam::mag
                                                (
                                                    trialMaxMove
                                                  - bestMaxMove
                                                )
                                              <= scalar(1e-15)
                                             && minPositiveRatio
                                              > bestMinRatio
                                            )
                                            {
                                                better = true;
                                            }
                                        }
                                    }
                                }
                            }


                            if( better )
                            {
                                found = true;
                                bestEdge = seI;
                                bestAmplitude = amplitude;
                                bestTargetVol = targetVol;
                                bestTargetOFVol = targetOFVol;
                                bestMinRatio = minPositiveRatio;
                                bestMinOFRatio = minOFPositiveRatio;
                                bestMaxOFSkew = trialMaxSkew;
                                bestMinOFOrtho = trialMinOrtho;
                                bestBadPyrCount =
                                    trialBadPyrCount;
                                bestMaxMove =
                                    trialMaxMove;

                                bestPositions.setSize(rowSize-2);

                                for
                                (
                                    label rowI=1;
                                    rowI<rowSize-1;
                                    ++rowI
                                )
                                {
                                    const label pointI =
                                        newVerticesForSplitEdge_
                                        (
                                            seI,
                                            rowI
                                        );

                                    bestPositions[rowI-1] =
                                        v1Points[pointI];
                                }
                            }
                        }
                    }


                    restoreChain();
                }


                // ======================================================
                // CFMITCH V2.7.1 WALL-FACE COHERENT REPAIR
                //
                // The v2.7 forensic classifier established that the
                // dominant unresolved type-1 quality population is:
                //
                //   - positive volume
                //   - wallLayer == 0
                //   - exactly one bad pyramid face
                //   - that face is the original WALL_BASE
                //
                // Existing V1C moves one complete hair at a time.  For a
                // WALL_BASE-only failure that can unnecessarily twist an
                // otherwise-valid first cross-section.
                //
                // This bounded stage moves every hair rooted on the actual
                // wall/base face with the SAME normalized amplitude.
                //
                // There is therefore only one scalar search variable:
                //
                //     alpha
                //
                // rather than the combinatorial alphaA/alphaB search used
                // by V1D_PAIR.
                //
                // Every candidate is evaluated over the union of all
                // participating incident-cell stars with the same
                // production safety policy as V1C/V1D:
                //
                //   * target raw/OF volume remains meaningful positive
                //   * positive neighbour volumes remain positive
                //   * no positive neighbour loses more than 25% volume
                //   * no new/worse bad pyramid
                //   * no new/worse skew
                //   * no new/worse >90-degree nonorthogonality
                //   * quality target reaches zero bad pyramid faces
                //
                // No topology or layer-count change occurs here.
                // ======================================================

                bool v271FaceBreathCommitted = false;

                if
                (
                    !found
                 && requirePyramidRepair
                 && badCellI >= 0
                 && badCellI < label(exactVolumeParent.size())
                 && badCellI < label(exactVolumeLocalChild.size())
                 && badCellI < label(exactVolumeRefType.size())
                 && exactVolumeRefType[badCellI] == 1
                )
                {
                    const label v271BfI =
                        cellToBaseBndFace_[badCellI];

                    const label v271LocalChild =
                        exactVolumeLocalChild[badCellI];

                    bool v271Eligible =
                        (
                            v271BfI >= 0
                         && v271BfI <
                            label(nLayersAtBndFace_.size())
                         && v271LocalChild >= 0
                        );

                    label v271N = -1;
                    label v271WallLayer = -1;

                    if( v271Eligible )
                    {
                        v271N =
                            nLayersAtBndFace_[v271BfI];

                        v271WallLayer =
                            v271N
                          - 1
                          - v271LocalChild;

                        if
                        (
                            v271N < 2
                         || v271WallLayer != 0
                        )
                            v271Eligible = false;
                    }


                    // --------------------------------------------------
                    // Find every pyramid-bad face belonging specifically
                    // to THIS target cell, using exactly the same owner /
                    // neighbour sign convention as
                    // v1CellBadPyramidCount().
                    //
                    // For this first experiment require exactly one bad
                    // face, and require every vertex of that face to be a
                    // known hair root.  That is our topology-derived
                    // WALL_BASE classification.
                    // --------------------------------------------------

                    label v271BadFaceI = -1;
                    label v271BadFaceCount = 0;
                    label v271WallBaseBadCount = 0;

                    point v271TargetCentre(vector::zero);
                    scalar v271TargetOFVolBefore = -GREAT;

                    if
                    (
                        v271Eligible
                     && !v1OFCellCentreVolume
                        (
                            badCellI,
                            v271TargetCentre,
                            v271TargetOFVolBefore
                        )
                    )
                    {
                        v271Eligible = false;
                    }

                    if( v271Eligible )
                    {
                        const cell& v271Cell =
                            v1Cells[badCellI];

                        forAll(v271Cell, cfI)
                        {
                            const label faceI =
                                v271Cell[cfI];

                            if
                            (
                                faceI < 0
                             || faceI >= label(v1Faces.size())
                             || faceI >= label(exactOwner.size())
                             || faceI >= label(v1Neighbour.size())
                            )
                            {
                                v271Eligible = false;
                                break;
                            }

                            const scalar pyrVol =
                                pyramidPointFaceRef
                                (
                                    v1Faces[faceI],
                                    v271TargetCentre
                                ).mag(v1Points);

                            bool badForTarget = false;

                            if( exactOwner[faceI] == badCellI )
                            {
                                if( pyrVol > SMALL )
                                    badForTarget = true;
                            }
                            else if
                            (
                                v1Neighbour[faceI] == badCellI
                            )
                            {
                                if( pyrVol < -SMALL )
                                    badForTarget = true;
                            }
                            else
                            {
                                v271Eligible = false;
                                break;
                            }

                            if( !badForTarget )
                                continue;

                            ++v271BadFaceCount;

                            bool allHairStarts = true;

                            const face& f =
                                v1Faces[faceI];

                            forAll(f, fpI)
                            {
                                if
                                (
                                    !v1HairStartPoints.found
                                    (
                                        f[fpI]
                                    )
                                )
                                {
                                    allHairStarts = false;
                                    break;
                                }
                            }

                            if( allHairStarts )
                            {
                                ++v271WallBaseBadCount;
                                v271BadFaceI = faceI;
                            }
                        }
                    }

                    if
                    (
                        v271Eligible
                     && (
                            v271BadFaceCount != 1
                         || v271WallBaseBadCount != 1
                         || v271BadFaceI < 0
                        )
                    )
                    {
                        v271Eligible = false;
                    }


                    if( v271Eligible )
                    {
                        ++blV271FaceBreathEligible;

                        const face& v271WallFace =
                            v1Faces[v271BadFaceI];

                        labelHashSet v271BadCellPoints;

                        const cell& v271BadCell =
                            v1Cells[badCellI];

                        forAll(v271BadCell, cfI)
                        {
                            const face& f =
                                v1Faces[v271BadCell[cfI]];

                            forAll(f, fpI)
                                v271BadCellPoints.insert(f[fpI]);
                        }


                        // ----------------------------------------------
                        // Derive exactly one participating hair from each
                        // wall/base vertex.
                        //
                        // The hair must:
                        //   - start at this wall vertex
                        //   - contain generated rows
                        //   - have its first generated point in the target
                        //     first child
                        //
                        // Ambiguous mappings fail closed.
                        // ----------------------------------------------

                        std::set<label> v271Edges;
                        bool v271ParticipantsValid = true;

                        forAll(v271WallFace, fpI)
                        {
                            const label wallPointI =
                                v271WallFace[fpI];

                            if
                            (
                                wallPointI < 0
                             || wallPointI >=
                                label(splitEdgesAtPoint_.size())
                            )
                            {
                                v271ParticipantsValid = false;
                                break;
                            }

                            label matchedEdge = -1;
                            label nMatches = 0;

                            forAllRow
                            (
                                splitEdgesAtPoint_,
                                wallPointI,
                                sepI
                            )
                            {
                                const label seI =
                                    splitEdgesAtPoint_
                                    (
                                        wallPointI,
                                        sepI
                                    );

                                if
                                (
                                    seI < 0
                                 || seI >= label(splitEdges_.size())
                                )
                                    continue;

                                const edge& se =
                                    splitEdges_[seI];

                                if( se.start() != wallPointI )
                                    continue;

                                const label rowSize =
                                    newVerticesForSplitEdge_.
                                        sizeOfRow(seI);

                                if( rowSize < 3 )
                                    continue;

                                const label firstPointI =
                                    newVerticesForSplitEdge_
                                    (
                                        seI,
                                        1
                                    );

                                if
                                (
                                    !v271BadCellPoints.found
                                    (
                                        firstPointI
                                    )
                                )
                                    continue;

                                matchedEdge = seI;
                                ++nMatches;
                            }

                            if( nMatches != 1 )
                            {
                                v271ParticipantsValid = false;
                                break;
                            }

                            v271Edges.insert(matchedEdge);
                        }

                        if
                        (
                            v271ParticipantsValid
                         && label(v271Edges.size()) !=
                            label(v271WallFace.size())
                        )
                        {
                            v271ParticipantsValid = false;
                        }


                        if( !v271ParticipantsValid )
                        {
                            ++blV271FaceBreathSkipped;

                            Info
                                << "CFMITCH V2.7.1 WALLFACE SKIP:"
                                << " cell=" << badCellI
                                << " parent="
                                << exactVolumeParent[badCellI]
                                << " bfI=" << v271BfI
                                << " localChild="
                                << v271LocalChild
                                << " reason=participantMapping"
                                << " wallFacePts="
                                << v271WallFace.size()
                                << " mappedEdges="
                                << v271Edges.size()
                                << endl;
                        }
                        else
                        {
                            ++blV271FaceBreathAttempted;


                            // ------------------------------------------
                            // Capture every complete participating chain.
                            // ------------------------------------------

                            std::map<label,point>
                                v271EdgeStart;

                            std::map<label,vector>
                                v271EdgeVec;

                            std::map<label,label>
                                v271RowSize;

                            std::map<label,List<point> >
                                v271OriginalPositions;

                            std::map<label,List<scalar> >
                                v271OriginalT;

                            // Generated point labels must be unique across
                            // participating hairs for this first prototype.
                            // Shared/degenerate chains fail closed.
                            std::map<label,label>
                                v271GeneratedPointEdge;

                            bool v271CaptureValid = true;


                            for
                            (
                                std::set<label>::const_iterator
                                    eIt=v271Edges.begin();
                                eIt!=v271Edges.end();
                                ++eIt
                            )
                            {
                                const label edgeI = *eIt;

                                const edge& se =
                                    splitEdges_[edgeI];

                                const point edgeStart =
                                    v1Points[se.start()];

                                const point edgeEnd =
                                    v1Points[se.end()];

                                const vector edgeVec =
                                    edgeEnd-edgeStart;

                                const scalar edgeMagSqr =
                                    magSqr(edgeVec);

                                if( edgeMagSqr <= VSMALL )
                                {
                                    v271CaptureValid = false;
                                    break;
                                }

                                const label rowSize =
                                    newVerticesForSplitEdge_.
                                        sizeOfRow(edgeI);

                                if( rowSize < 3 )
                                {
                                    v271CaptureValid = false;
                                    break;
                                }

                                List<point> originalPositions
                                (
                                    rowSize-2
                                );

                                List<scalar> originalT
                                (
                                    rowSize
                                );

                                originalT[0] =
                                    scalar(0);

                                originalT[rowSize-1] =
                                    scalar(1);

                                for
                                (
                                    label rowI=1;
                                    rowI<rowSize-1;
                                    ++rowI
                                )
                                {
                                    const label pointI =
                                        newVerticesForSplitEdge_
                                        (
                                            edgeI,
                                            rowI
                                        );

                                    std::map<label,label>::
                                        iterator gpIt =
                                            v271GeneratedPointEdge.find
                                            (
                                                pointI
                                            );

                                    if
                                    (
                                        gpIt !=
                                        v271GeneratedPointEdge.end()
                                     && gpIt->second != edgeI
                                    )
                                    {
                                        v271CaptureValid = false;
                                        break;
                                    }

                                    v271GeneratedPointEdge[pointI] =
                                        edgeI;

                                    originalPositions[rowI-1] =
                                        v1Points[pointI];

                                    originalT[rowI] =
                                        (
                                            (v1Points[pointI]-edgeStart)
                                          & edgeVec
                                        )
                                       /(edgeMagSqr + VSMALL);

                                    if
                                    (
                                        !(originalT[rowI] >
                                            originalT[rowI-1])
                                     || !(originalT[rowI] <
                                            scalar(1))
                                    )
                                    {
                                        v271CaptureValid = false;
                                        break;
                                    }
                                }

                                if( !v271CaptureValid )
                                    break;

                                v271EdgeStart[edgeI] =
                                    edgeStart;

                                v271EdgeVec[edgeI] =
                                    edgeVec;

                                v271RowSize[edgeI] =
                                    rowSize;

                                v271OriginalPositions[edgeI] =
                                    originalPositions;

                                v271OriginalT[edgeI] =
                                    originalT;
                            }


                            if( !v271CaptureValid )
                            {
                                ++blV271FaceBreathSkipped;

                                Info
                                    << "CFMITCH V2.7.1 WALLFACE SKIP:"
                                    << " cell=" << badCellI
                                    << " parent="
                                    << exactVolumeParent[badCellI]
                                    << " bfI=" << v271BfI
                                    << " localChild="
                                    << v271LocalChild
                                    << " reason=chainCapture"
                                    << " edges="
                                    << v271Edges.size()
                                    << endl;
                            }
                            else
                            {
                                auto v271Restore =
                                [&]()
                                {
                                    for
                                    (
                                        std::set<label>::const_iterator
                                            eIt=v271Edges.begin();
                                        eIt!=v271Edges.end();
                                        ++eIt
                                    )
                                    {
                                        const label edgeI =
                                            *eIt;

                                        const label rowSize =
                                            v271RowSize[edgeI];

                                        const List<point>&
                                            originalPositions =
                                                v271OriginalPositions
                                                [
                                                    edgeI
                                                ];

                                        for
                                        (
                                            label rowI=1;
                                            rowI<rowSize-1;
                                            ++rowI
                                        )
                                        {
                                            const label pointI =
                                                newVerticesForSplitEdge_
                                                (
                                                    edgeI,
                                                    rowI
                                                );

                                            v1Points[pointI] =
                                                originalPositions
                                                [
                                                    rowI-1
                                                ];
                                        }
                                    }
                                };


                                auto v271Apply =
                                [&]
                                (
                                    const scalar amplitude
                                ) -> bool
                                {
                                    for
                                    (
                                        std::set<label>::const_iterator
                                            eIt=v271Edges.begin();
                                        eIt!=v271Edges.end();
                                        ++eIt
                                    )
                                    {
                                        const label edgeI =
                                            *eIt;

                                        const point& edgeStart =
                                            v271EdgeStart[edgeI];

                                        const vector& edgeVec =
                                            v271EdgeVec[edgeI];

                                        const label rowSize =
                                            v271RowSize[edgeI];

                                        const List<scalar>&
                                            originalT =
                                                v271OriginalT[edgeI];

                                        const scalar firstT =
                                            originalT[1];

                                        if
                                        (
                                            firstT <= scalar(0)
                                         || firstT >= scalar(1)
                                        )
                                        {
                                            v271Restore();
                                            return false;
                                        }

                                        scalar prevT =
                                            scalar(0);

                                        for
                                        (
                                            label rowI=1;
                                            rowI<rowSize-1;
                                            ++rowI
                                        )
                                        {
                                            const scalar t =
                                                originalT[rowI];

                                            const scalar decayBase =
                                                Foam::max
                                                (
                                                    scalar(0),
                                                    (scalar(1)-t)
                                                   /(
                                                        scalar(1)-firstT
                                                      + VSMALL
                                                    )
                                                );

                                            const scalar deltaT =
                                                amplitude
                                               *firstT
                                               *Foam::pow
                                                (
                                                    decayBase,
                                                    scalar(4)
                                                );

                                            const scalar warpedT =
                                                t + deltaT;

                                            if
                                            (
                                                !(warpedT > prevT)
                                             || !(warpedT <
                                                    originalT[rowI+1])
                                            )
                                            {
                                                v271Restore();
                                                return false;
                                            }

                                            const label pointI =
                                                newVerticesForSplitEdge_
                                                (
                                                    edgeI,
                                                    rowI
                                                );

                                            v1Points[pointI] =
                                                edgeStart
                                              + warpedT*edgeVec;

                                            prevT =
                                                warpedT;
                                        }
                                    }

                                    return true;
                                };


                                // --------------------------------------
                                // Union of all participating hair stars.
                                // --------------------------------------

                                std::set<label> v271UnionStar;
                                bool v271BaselineValid = true;

                                for
                                (
                                    std::set<label>::const_iterator
                                        eIt=v271Edges.begin();
                                    eIt!=v271Edges.end();
                                    ++eIt
                                )
                                {
                                    std::map
                                    <
                                        label,
                                        std::set<label>
                                    >::const_iterator starIt =
                                        v1EdgeCells.find(*eIt);

                                    if
                                    (
                                        starIt == v1EdgeCells.end()
                                     || starIt->second.empty()
                                    )
                                    {
                                        v271BaselineValid = false;
                                        break;
                                    }

                                    v271UnionStar.insert
                                    (
                                        starIt->second.begin(),
                                        starIt->second.end()
                                    );
                                }


                                std::map<label,scalar>
                                    v271BaselineVolumes;

                                std::map<label,scalar>
                                    v271BaselineOFVolumes;


                                if( v271BaselineValid )
                                {
                                    for
                                    (
                                        std::set<label>::const_iterator
                                            cIt=v271UnionStar.begin();
                                        cIt!=v271UnionStar.end();
                                        ++cIt
                                    )
                                    {
                                        v271BaselineVolumes[*cIt] =
                                            v1CellVolume(*cIt);

                                        point cc(vector::zero);
                                        scalar cv = -GREAT;

                                        if
                                        (
                                            !v1OFCellCentreVolume
                                            (
                                                *cIt,
                                                cc,
                                                cv
                                            )
                                        )
                                        {
                                            v271BaselineValid = false;
                                            break;
                                        }

                                        v271BaselineOFVolumes[*cIt] =
                                            cv;
                                    }
                                }


                                std::set<label>
                                    v271AffectedFaces;

                                if( v271BaselineValid )
                                {
                                    for
                                    (
                                        std::set<label>::const_iterator
                                            cIt=v271UnionStar.begin();
                                        cIt!=v271UnionStar.end();
                                        ++cIt
                                    )
                                    {
                                        const cell& starCell =
                                            v1Cells[*cIt];

                                        forAll(starCell, sfI)
                                        {
                                            v271AffectedFaces.insert
                                            (
                                                starCell[sfI]
                                            );
                                        }
                                    }
                                }


                                std::map<label,scalar>
                                    v271BaselineOFSkew;

                                std::map<label,scalar>
                                    v271BaselineOFOrtho;

                                std::map<label,scalar>
                                    v271BaselineOFPyrMargin;

                                std::map<label,bool>
                                    v271BaselineOFPyrBad;

                                std::map<label,bool>
                                    v271BaselineOFEligible;


                                if( v271BaselineValid )
                                {
                                    for
                                    (
                                        std::set<label>::const_iterator
                                            fIt=
                                                v271AffectedFaces.begin();
                                        fIt!=
                                            v271AffectedFaces.end();
                                        ++fIt
                                    )
                                    {
                                        const label faceI =
                                            *fIt;

                                        scalar skew = GREAT;
                                        scalar ortho = -GREAT;
                                        scalar pyrMargin = -GREAT;
                                        scalar ownVol = -GREAT;
                                        scalar neiVol = -GREAT;
                                        bool pyrBad = true;

                                        if
                                        (
                                            !v1OFFaceQuality
                                            (
                                                faceI,
                                                skew,
                                                ortho,
                                                pyrMargin,
                                                pyrBad,
                                                ownVol,
                                                neiVol
                                            )
                                        )
                                        {
                                            v271BaselineValid = false;
                                            break;
                                        }

                                        v271BaselineOFSkew[faceI] =
                                            skew;

                                        v271BaselineOFOrtho[faceI] =
                                            ortho;

                                        v271BaselineOFPyrMargin[faceI] =
                                            pyrMargin;

                                        v271BaselineOFPyrBad[faceI] =
                                            pyrBad;

                                        v271BaselineOFEligible[faceI] =
                                            (
                                                ownVol > scalar(0)
                                             && (
                                                    v1Neighbour[faceI] < 0
                                                 || neiVol > scalar(0)
                                                )
                                            );
                                    }
                                }


                                if( !v271BaselineValid )
                                {
                                    ++blV271FaceBreathSkipped;

                                    Info
                                        << "CFMITCH V2.7.1 WALLFACE SKIP:"
                                        << " cell=" << badCellI
                                        << " parent="
                                        << exactVolumeParent[badCellI]
                                        << " bfI=" << v271BfI
                                        << " localChild="
                                        << v271LocalChild
                                        << " reason=baseline"
                                        << " edges="
                                        << v271Edges.size()
                                        << endl;
                                }
                                else
                                {
                                    const scalar v271TargetFloor =
                                        scalar(0.20)
                                       *Foam::mag(badVolBefore);

                                    static const scalar
                                        v271SkewLimit =
                                            scalar(4);

                                    static const scalar
                                        v271Ortho90 =
                                            scalar(0);


                                    bool v271Found = false;

                                    scalar v271BestAmplitude =
                                        scalar(0);

                                    scalar v271BestTargetVol =
                                        -GREAT;

                                    scalar v271BestTargetOFVol =
                                        -GREAT;

                                    scalar v271BestMinRatio =
                                        -GREAT;

                                    scalar v271BestMinOFRatio =
                                        -GREAT;

                                    scalar v271BestMaxSkew =
                                        GREAT;

                                    scalar v271BestMinOrtho =
                                        -GREAT;

                                    scalar v271BestMaxMove =
                                        GREAT;

                                    label v271BestBadPyrCount =
                                        labelMax;

                                    std::map
                                    <
                                        label,
                                        List<point>
                                    >
                                        v271BestPositions;


                                    label v271Trials = 0;
                                    label v271VolumePass = 0;
                                    label v271QualityPass = 0;
                                    label v271TargetReject = 0;
                                    label v271QualityReject = 0;


                                    for
                                    (
                                        label ampI=0;
                                        ampI<nV1AmplitudeMag;
                                        ++ampI
                                    )
                                    {
                                        for
                                        (
                                            label sign=-1;
                                            sign<=1;
                                            sign+=2
                                        )
                                        {
                                            v271Restore();

                                            const scalar amplitude =
                                                scalar(sign)
                                               *v1AmplitudeMag[ampI];

                                            if
                                            (
                                                !v271Apply(amplitude)
                                            )
                                                continue;

                                            ++v271Trials;
                                            ++blV271FaceBreathTrials;


                                            // --------------------------
                                            // Raw cfMesh volume gate.
                                            // --------------------------

                                            const scalar targetVol =
                                                v1CellVolume
                                                (
                                                    badCellI
                                                );

                                            if
                                            (
                                                targetVol
                                             <= v271TargetFloor
                                            )
                                                continue;

                                            bool safe = true;

                                            scalar minPositiveRatio =
                                                GREAT;

                                            for
                                            (
                                                std::set<label>::
                                                    const_iterator
                                                    cIt=
                                                        v271UnionStar.begin();
                                                cIt!=
                                                    v271UnionStar.end();
                                                ++cIt
                                            )
                                            {
                                                const label starCellI =
                                                    *cIt;

                                                const scalar oldV =
                                                    v271BaselineVolumes
                                                    [
                                                        starCellI
                                                    ];

                                                const scalar newV =
                                                    v1CellVolume
                                                    (
                                                        starCellI
                                                    );

                                                if
                                                (
                                                    oldV > scalar(0)
                                                )
                                                {
                                                    if
                                                    (
                                                        newV <= scalar(0)
                                                    )
                                                    {
                                                        safe = false;
                                                        break;
                                                    }

                                                    const scalar ratio =
                                                        newV
                                                       /(oldV + VSMALL);

                                                    minPositiveRatio =
                                                        Foam::min
                                                        (
                                                            minPositiveRatio,
                                                            ratio
                                                        );

                                                    if
                                                    (
                                                        ratio
                                                      < scalar(0.75)
                                                    )
                                                    {
                                                        safe = false;
                                                        break;
                                                    }
                                                }
                                                else if
                                                (
                                                    newV < oldV
                                                )
                                                {
                                                    safe = false;
                                                    break;
                                                }
                                            }

                                            if( !safe )
                                                continue;


                                            // --------------------------
                                            // OpenFOAM signed-volume gate.
                                            // --------------------------

                                            point targetOFCentre
                                            (
                                                vector::zero
                                            );

                                            scalar targetOFVol =
                                                -GREAT;

                                            if
                                            (
                                                !v1OFCellCentreVolume
                                                (
                                                    badCellI,
                                                    targetOFCentre,
                                                    targetOFVol
                                                )
                                            )
                                                continue;

                                            if
                                            (
                                                targetOFVol
                                             <= v271TargetFloor
                                            )
                                                continue;

                                            scalar minOFPositiveRatio =
                                                GREAT;

                                            for
                                            (
                                                std::set<label>::
                                                    const_iterator
                                                    cIt=
                                                        v271UnionStar.begin();
                                                cIt!=
                                                    v271UnionStar.end();
                                                ++cIt
                                            )
                                            {
                                                const label starCellI =
                                                    *cIt;

                                                const scalar oldOFV =
                                                    v271BaselineOFVolumes
                                                    [
                                                        starCellI
                                                    ];

                                                point newOFCentre
                                                (
                                                    vector::zero
                                                );

                                                scalar newOFV =
                                                    -GREAT;

                                                if
                                                (
                                                    !v1OFCellCentreVolume
                                                    (
                                                        starCellI,
                                                        newOFCentre,
                                                        newOFV
                                                    )
                                                )
                                                {
                                                    safe = false;
                                                    break;
                                                }

                                                if
                                                (
                                                    oldOFV > scalar(0)
                                                )
                                                {
                                                    if
                                                    (
                                                        newOFV <= scalar(0)
                                                    )
                                                    {
                                                        safe = false;
                                                        break;
                                                    }

                                                    const scalar ratio =
                                                        newOFV
                                                       /(oldOFV + VSMALL);

                                                    minOFPositiveRatio =
                                                        Foam::min
                                                        (
                                                            minOFPositiveRatio,
                                                            ratio
                                                        );

                                                    if
                                                    (
                                                        ratio
                                                      < scalar(0.75)
                                                    )
                                                    {
                                                        safe = false;
                                                        break;
                                                    }
                                                }
                                                else if
                                                (
                                                    newOFV < oldOFV
                                                )
                                                {
                                                    safe = false;
                                                    break;
                                                }
                                            }

                                            if( !safe )
                                                continue;


                                            ++v271VolumePass;
                                            ++blV271FaceBreathVolumePass;


                                            // --------------------------
                                            // OpenFOAM-parity quality gate.
                                            // --------------------------

                                            bool ofSafe = true;

                                            const char* v271RejectReason =
                                                "none";

                                            label v271RejectFace =
                                                -1;

                                            scalar trialMaxSkew =
                                                scalar(0);

                                            scalar trialMinOrtho =
                                                GREAT;

                                            label trialBadPyrCount =
                                                0;


                                            for
                                            (
                                                std::set<label>::
                                                    const_iterator
                                                    fIt=
                                                        v271AffectedFaces.begin();
                                                fIt!=
                                                    v271AffectedFaces.end();
                                                ++fIt
                                            )
                                            {
                                                const label faceI =
                                                    *fIt;

                                                scalar trialSkew =
                                                    GREAT;

                                                scalar trialOrtho =
                                                    -GREAT;

                                                scalar trialPyrMargin =
                                                    -GREAT;

                                                scalar trialOwnVol =
                                                    -GREAT;

                                                scalar trialNeiVol =
                                                    -GREAT;

                                                bool trialPyrBad =
                                                    true;


                                                if
                                                (
                                                    !v1OFFaceQuality
                                                    (
                                                        faceI,
                                                        trialSkew,
                                                        trialOrtho,
                                                        trialPyrMargin,
                                                        trialPyrBad,
                                                        trialOwnVol,
                                                        trialNeiVol
                                                    )
                                                )
                                                {
                                                    ofSafe = false;
                                                    v271RejectReason =
                                                        "ofGeometry";
                                                    v271RejectFace =
                                                        faceI;
                                                    break;
                                                }


                                                const bool trialEligible =
                                                    (
                                                        trialOwnVol
                                                      > scalar(0)
                                                     && (
                                                            v1Neighbour
                                                            [
                                                                faceI
                                                            ] < 0
                                                         || trialNeiVol
                                                          > scalar(0)
                                                        )
                                                    );

                                                if( !trialEligible )
                                                    continue;


                                                if( trialPyrBad )
                                                {
                                                    ++trialBadPyrCount;
                                                }


                                                trialMaxSkew =
                                                    Foam::max
                                                    (
                                                        trialMaxSkew,
                                                        trialSkew
                                                    );


                                                if
                                                (
                                                    v1Neighbour[faceI]
                                                  >= 0
                                                )
                                                {
                                                    trialMinOrtho =
                                                        Foam::min
                                                        (
                                                            trialMinOrtho,
                                                            trialOrtho
                                                        );
                                                }


                                                const bool
                                                    baselineEligible =
                                                        v271BaselineOFEligible
                                                        [
                                                            faceI
                                                        ];


                                                if( !baselineEligible )
                                                {
                                                    if( trialPyrBad )
                                                    {
                                                        ofSafe = false;
                                                        v271RejectReason =
                                                            "newPyramid";
                                                        v271RejectFace =
                                                            faceI;
                                                        break;
                                                    }

                                                    if
                                                    (
                                                        trialSkew
                                                      > v271SkewLimit
                                                    )
                                                    {
                                                        ofSafe = false;
                                                        v271RejectReason =
                                                            "newSkew";
                                                        v271RejectFace =
                                                            faceI;
                                                        break;
                                                    }

                                                    if
                                                    (
                                                        v1Neighbour[faceI]
                                                      >= 0
                                                     && trialOrtho
                                                      < v271Ortho90
                                                    )
                                                    {
                                                        ofSafe = false;
                                                        v271RejectReason =
                                                            "newNonOrtho90";
                                                        v271RejectFace =
                                                            faceI;
                                                        break;
                                                    }

                                                    continue;
                                                }


                                                const bool basePyrBad =
                                                    v271BaselineOFPyrBad
                                                    [
                                                        faceI
                                                    ];

                                                const scalar
                                                    basePyrMargin =
                                                        v271BaselineOFPyrMargin
                                                        [
                                                            faceI
                                                        ];


                                                if
                                                (
                                                    !basePyrBad
                                                 && trialPyrBad
                                                )
                                                {
                                                    ofSafe = false;
                                                    v271RejectReason =
                                                        "newPyramid";
                                                    v271RejectFace =
                                                        faceI;
                                                    break;
                                                }


                                                if
                                                (
                                                    basePyrBad
                                                 && trialPyrBad
                                                )
                                                {
                                                    const scalar pyrTol =
                                                        scalar(1e-12)
                                                       *(
                                                            Foam::mag
                                                            (
                                                                basePyrMargin
                                                            )
                                                          + SMALL
                                                        );

                                                    if
                                                    (
                                                        trialPyrMargin
                                                      <
                                                        basePyrMargin
                                                       -pyrTol
                                                    )
                                                    {
                                                        ofSafe = false;
                                                        v271RejectReason =
                                                            "worsePyramid";
                                                        v271RejectFace =
                                                            faceI;
                                                        break;
                                                    }
                                                }


                                                const scalar baseSkew =
                                                    v271BaselineOFSkew
                                                    [
                                                        faceI
                                                    ];

                                                const scalar skewTol =
                                                    scalar(1e-10)
                                                   *(
                                                        scalar(1)
                                                      + Foam::mag
                                                        (
                                                            baseSkew
                                                        )
                                                    );


                                                if
                                                (
                                                    baseSkew
                                                  <= v271SkewLimit
                                                )
                                                {
                                                    if
                                                    (
                                                        trialSkew
                                                      > v271SkewLimit
                                                    )
                                                    {
                                                        ofSafe = false;
                                                        v271RejectReason =
                                                            "newSkew";
                                                        v271RejectFace =
                                                            faceI;
                                                        break;
                                                    }
                                                }
                                                else if
                                                (
                                                    trialSkew
                                                  > baseSkew+skewTol
                                                )
                                                {
                                                    ofSafe = false;
                                                    v271RejectReason =
                                                        "worseSkew";
                                                    v271RejectFace =
                                                        faceI;
                                                    break;
                                                }


                                                if
                                                (
                                                    v1Neighbour[faceI]
                                                  >= 0
                                                )
                                                {
                                                    const scalar
                                                        baseOrtho =
                                                            v271BaselineOFOrtho
                                                            [
                                                                faceI
                                                            ];

                                                    const scalar orthoTol =
                                                        scalar(1e-12);


                                                    if
                                                    (
                                                        baseOrtho
                                                      >= v271Ortho90
                                                    )
                                                    {
                                                        if
                                                        (
                                                            trialOrtho
                                                          < v271Ortho90
                                                        )
                                                        {
                                                            ofSafe = false;
                                                            v271RejectReason =
                                                                "newNonOrtho90";
                                                            v271RejectFace =
                                                                faceI;
                                                            break;
                                                        }
                                                    }
                                                    else if
                                                    (
                                                        trialOrtho
                                                      < baseOrtho
                                                       -orthoTol
                                                    )
                                                    {
                                                        ofSafe = false;
                                                        v271RejectReason =
                                                            "worseNonOrtho";
                                                        v271RejectFace =
                                                            faceI;
                                                        break;
                                                    }
                                                }
                                            }


                                            if( !ofSafe )
                                            {
                                                ++v271QualityReject;
                                                ++blV271FaceBreathQualityReject;

                                                if
                                                (
                                                    std::strcmp
                                                    (
                                                        v271RejectReason,
                                                        "newPyramid"
                                                    ) == 0
                                                )
                                                {
                                                    ++blV271RejectNewPyramid;
                                                }
                                                else if
                                                (
                                                    std::strcmp
                                                    (
                                                        v271RejectReason,
                                                        "worsePyramid"
                                                    ) == 0
                                                )
                                                {
                                                    ++blV271RejectWorsePyramid;
                                                }
                                                else if
                                                (
                                                    std::strcmp
                                                    (
                                                        v271RejectReason,
                                                        "newSkew"
                                                    ) == 0
                                                )
                                                {
                                                    ++blV271RejectNewSkew;
                                                }
                                                else if
                                                (
                                                    std::strcmp
                                                    (
                                                        v271RejectReason,
                                                        "worseSkew"
                                                    ) == 0
                                                )
                                                {
                                                    ++blV271RejectWorseSkew;
                                                }
                                                else if
                                                (
                                                    std::strcmp
                                                    (
                                                        v271RejectReason,
                                                        "newNonOrtho90"
                                                    ) == 0
                                                )
                                                {
                                                    ++blV271RejectNewNonOrtho90;
                                                }
                                                else if
                                                (
                                                    std::strcmp
                                                    (
                                                        v271RejectReason,
                                                        "worseNonOrtho"
                                                    ) == 0
                                                )
                                                {
                                                    ++blV271RejectWorseNonOrtho;
                                                }
                                                else if
                                                (
                                                    std::strcmp
                                                    (
                                                        v271RejectReason,
                                                        "ofGeometry"
                                                    ) == 0
                                                )
                                                {
                                                    ++blV271RejectOFGeometry;
                                                }
                                                else
                                                {
                                                    ++blV271RejectOther;
                                                }


                                                bool rejectOnTargetFace =
                                                    false;

                                                if
                                                (
                                                    v271RejectFace >= 0
                                                )
                                                {
                                                    const cell&
                                                        targetCell =
                                                            v1Cells
                                                            [
                                                                badCellI
                                                            ];

                                                    forAll
                                                    (
                                                        targetCell,
                                                        tcfI
                                                    )
                                                    {
                                                        if
                                                        (
                                                            targetCell[tcfI]
                                                         ==
                                                            v271RejectFace
                                                        )
                                                        {
                                                            rejectOnTargetFace =
                                                                true;
                                                            break;
                                                        }
                                                    }
                                                }

                                                if( rejectOnTargetFace )
                                                {
                                                    ++blV271RejectOnTargetCellFace;
                                                }
                                                else
                                                {
                                                    ++blV271RejectOnOtherStarFace;
                                                }

                                                continue;
                                            }


                                            // The quality-only target must
                                            // actually be repaired, not just
                                            // made "no worse".
                                            const label
                                                trialTargetBadPyr =
                                                    v1CellBadPyramidCount
                                                    (
                                                        badCellI
                                                    );

                                            if
                                            (
                                                trialTargetBadPyr != 0
                                            )
                                            {
                                                ++v271TargetReject;
                                                ++blV271FaceBreathTargetReject;
                                                continue;
                                            }


                                            ++v271QualityPass;
                                            ++blV271FaceBreathQualityPass;


                                            // --------------------------
                                            // Maximum physical movement.
                                            // --------------------------

                                            scalar trialMaxMove =
                                                scalar(0);

                                            for
                                            (
                                                std::set<label>::
                                                    const_iterator
                                                    eIt=
                                                        v271Edges.begin();
                                                eIt!=
                                                    v271Edges.end();
                                                ++eIt
                                            )
                                            {
                                                const label edgeI =
                                                    *eIt;

                                                const label rowSize =
                                                    v271RowSize[edgeI];

                                                const List<point>&
                                                    originalPositions =
                                                        v271OriginalPositions
                                                        [
                                                            edgeI
                                                        ];

                                                for
                                                (
                                                    label rowI=1;
                                                    rowI<rowSize-1;
                                                    ++rowI
                                                )
                                                {
                                                    const label pointI =
                                                        newVerticesForSplitEdge_
                                                        (
                                                            edgeI,
                                                            rowI
                                                        );

                                                    trialMaxMove =
                                                        Foam::max
                                                        (
                                                            trialMaxMove,
                                                            mag
                                                            (
                                                                v1Points
                                                                [
                                                                    pointI
                                                                ]
                                                              -
                                                                originalPositions
                                                                [
                                                                    rowI-1
                                                                ]
                                                            )
                                                        );
                                                }
                                            }


                                            // Same deterministic ranking
                                            // used by V1C/V1D.
                                            bool better =
                                                !v271Found;

                                            if( v271Found )
                                            {
                                                if
                                                (
                                                    trialBadPyrCount
                                                  <
                                                    v271BestBadPyrCount
                                                )
                                                {
                                                    better = true;
                                                }
                                                else if
                                                (
                                                    trialBadPyrCount
                                                 ==
                                                    v271BestBadPyrCount
                                                )
                                                {
                                                    if
                                                    (
                                                        trialMaxSkew
                                                      <
                                                        v271BestMaxSkew
                                                       -scalar(1e-12)
                                                    )
                                                    {
                                                        better = true;
                                                    }
                                                    else if
                                                    (
                                                        Foam::mag
                                                        (
                                                            trialMaxSkew
                                                          -
                                                            v271BestMaxSkew
                                                        )
                                                      <= scalar(1e-12)
                                                    )
                                                    {
                                                        if
                                                        (
                                                            trialMinOrtho
                                                          >
                                                            v271BestMinOrtho
                                                           +scalar(1e-12)
                                                        )
                                                        {
                                                            better = true;
                                                        }
                                                        else if
                                                        (
                                                            Foam::mag
                                                            (
                                                                trialMinOrtho
                                                              -
                                                                v271BestMinOrtho
                                                            )
                                                          <= scalar(1e-12)
                                                        )
                                                        {
                                                            if
                                                            (
                                                                trialMaxMove
                                                              <
                                                                v271BestMaxMove
                                                               -scalar(1e-15)
                                                            )
                                                            {
                                                                better = true;
                                                            }
                                                            else if
                                                            (
                                                                Foam::mag
                                                                (
                                                                    trialMaxMove
                                                                  -
                                                                    v271BestMaxMove
                                                                )
                                                              <= scalar(1e-15)
                                                             &&
                                                                minPositiveRatio
                                                              >
                                                                v271BestMinRatio
                                                            )
                                                            {
                                                                better = true;
                                                            }
                                                        }
                                                    }
                                                }
                                            }


                                            if( better )
                                            {
                                                v271Found = true;

                                                v271BestAmplitude =
                                                    amplitude;

                                                v271BestTargetVol =
                                                    targetVol;

                                                v271BestTargetOFVol =
                                                    targetOFVol;

                                                v271BestMinRatio =
                                                    minPositiveRatio;

                                                v271BestMinOFRatio =
                                                    minOFPositiveRatio;

                                                v271BestMaxSkew =
                                                    trialMaxSkew;

                                                v271BestMinOrtho =
                                                    trialMinOrtho;

                                                v271BestBadPyrCount =
                                                    trialBadPyrCount;

                                                v271BestMaxMove =
                                                    trialMaxMove;

                                                v271BestPositions.clear();

                                                for
                                                (
                                                    std::set<label>::
                                                        const_iterator
                                                        eIt=
                                                            v271Edges.begin();
                                                    eIt!=
                                                        v271Edges.end();
                                                    ++eIt
                                                )
                                                {
                                                    const label edgeI =
                                                        *eIt;

                                                    const label rowSize =
                                                        v271RowSize
                                                        [
                                                            edgeI
                                                        ];

                                                    List<point>
                                                        positions
                                                        (
                                                            rowSize-2
                                                        );

                                                    for
                                                    (
                                                        label rowI=1;
                                                        rowI<rowSize-1;
                                                        ++rowI
                                                    )
                                                    {
                                                        const label pointI =
                                                            newVerticesForSplitEdge_
                                                            (
                                                                edgeI,
                                                                rowI
                                                            );

                                                        positions[rowI-1] =
                                                            v1Points[pointI];
                                                    }

                                                    v271BestPositions
                                                    [
                                                        edgeI
                                                    ] = positions;
                                                }
                                            }
                                        }
                                    }


                                    // All trial geometry must be gone before
                                    // committing the selected transaction.
                                    v271Restore();


                                    Info
                                        << "CFMITCH V2.7.1 WALLFACE SEARCH:"
                                        << " cell=" << badCellI
                                        << " parent="
                                        << exactVolumeParent[badCellI]
                                        << " bfI=" << v271BfI
                                        << " localChild="
                                        << v271LocalChild
                                        << " wallLayer="
                                        << v271WallLayer
                                        << " edges="
                                        << v271Edges.size()
                                        << " trials="
                                        << v271Trials
                                        << " volumePass="
                                        << v271VolumePass
                                        << " qualityPass="
                                        << v271QualityPass
                                        << " rejectQuality="
                                        << v271QualityReject
                                        << " rejectTargetPyramid="
                                        << v271TargetReject
                                        << " found="
                                        << (
                                            v271Found
                                          ? "true"
                                          : "false"
                                           )
                                        << endl;


                                    if( v271Found )
                                    {
                                        scalar maxPhysicalMove =
                                            scalar(0);

                                        for
                                        (
                                            std::set<label>::
                                                const_iterator
                                                eIt=
                                                    v271Edges.begin();
                                            eIt!=
                                                v271Edges.end();
                                            ++eIt
                                        )
                                        {
                                            const label edgeI =
                                                *eIt;

                                            const label rowSize =
                                                v271RowSize[edgeI];

                                            const List<point>&
                                                bestPositions =
                                                    v271BestPositions
                                                    [
                                                        edgeI
                                                    ];

                                            for
                                            (
                                                label rowI=1;
                                                rowI<rowSize-1;
                                                ++rowI
                                            )
                                            {
                                                const label pointI =
                                                    newVerticesForSplitEdge_
                                                    (
                                                        edgeI,
                                                        rowI
                                                    );

                                                maxPhysicalMove =
                                                    Foam::max
                                                    (
                                                        maxPhysicalMove,
                                                        mag
                                                        (
                                                            bestPositions
                                                            [
                                                                rowI-1
                                                            ]
                                                          -
                                                            v1Points[pointI]
                                                        )
                                                    );

                                                v1Points[pointI] =
                                                    bestPositions[rowI-1];
                                            }
                                        }


                                        const scalar committedVol =
                                            v1CellVolume(badCellI);

                                        const label committedBadPyr =
                                            v1CellBadPyramidCount
                                            (
                                                badCellI
                                            );


                                        if
                                        (
                                            committedVol > scalar(0)
                                         && committedBadPyr == 0
                                        )
                                        {
                                            ++blV1Fixed;

                                            blV1CommittedChains +=
                                                v271Edges.size();

                                            ++blV271FaceBreathFixed;

                                            v271FaceBreathCommitted =
                                                true;

                                            Info
                                                << "CFMITCH V2.7.1 WALLFACE MOVE:"
                                                << " cell="
                                                << badCellI
                                                << " parent="
                                                << exactVolumeParent
                                                   [
                                                       badCellI
                                                   ]
                                                << " bfI="
                                                << v271BfI
                                                << " localChild="
                                                << v271LocalChild
                                                << " edges="
                                                << v271Edges.size()
                                                << " amplitude="
                                                << v271BestAmplitude
                                                << " maxPhysicalMove="
                                                << maxPhysicalMove
                                                << " oldCellVol="
                                                << badVolBefore
                                                << " targetVol="
                                                << v271BestTargetVol
                                                << " committedVol="
                                                << committedVol
                                                << " targetOFVol="
                                                << v271BestTargetOFVol
                                                << " maxOFSkew="
                                                << v271BestMaxSkew
                                                << " minOFOrtho="
                                                << v271BestMinOrtho
                                                << " badOFPyr="
                                                << v271BestBadPyrCount
                                                << " minPositiveStarRatio="
                                                << v271BestMinRatio
                                                << " minOFPositiveStarRatio="
                                                << v271BestMinOFRatio
                                                << endl;
                                        }
                                        else
                                        {
                                            // Defensive fail-closed restore.
                                            v271Restore();

                                            Info
                                                << "CFMITCH V2.7.1 WALLFACE COMMIT_REJECT:"
                                                << " cell="
                                                << badCellI
                                                << " committedVol="
                                                << committedVol
                                                << " committedBadPyr="
                                                << committedBadPyr
                                                << endl;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }


                // A successful coherent face transaction has already been
                // committed atomically and accounted for.  Move directly
                // to the next repair seed.
                if( v271FaceBreathCommitted )
                    continue;


                // ======================================================
                // BL_VALIDITY_REPAIR_V1D_PAIR
                //
                // Stage 2 fallback: if no coherent single hair can make
                // the negative child solver-safe, search pairs of hairs
                // which actually occur together on a face of the bad
                // child.
                //
                // No Rotor37 edge ids are encoded here.  Candidate pairs
                // are derived entirely from topology + v1PointEdge.
                //
                // Each trial:
                //   - restores BOTH complete chains
                //   - moves BOTH complete chains coherently
                //   - evaluates the UNION of their incident-cell stars
                //   - evaluates every face of that union with the same
                //     OpenFOAM-parity quality rules used by v1c
                //   - survives only as an atomic two-chain transaction
                // ======================================================

                // CFMitch v2.6:
                //
                // The paired-hair fallback is now available to both:
                //   - historical negative-volume repair seeds
                //   - positive-volume bad-pyramid quality seeds
                //
                // Quality-only seeds receive an additional target-cell
                // pyramid-elimination gate below before a pair candidate
                // is permitted into ranking.
                // CFMitch v2.7.1:
                //
                // Positive-volume pyramid-quality seeds do not enter the
                // combinatorial paired-hair fallback.  Eligible type-1
                // wall-adjacent WALL_BASE-only targets have already received
                // the bounded coherent wall-face search above.
                //
                // Preserve V1D_PAIR only for its original negative-volume
                // rescue role.
                if( !found && !requirePyramidRepair )
                {
                    std::set<std::pair<label,label> > pairEdges;


                    // --------------------------------------------------
                    // Discover hair pairs which share an actual face of
                    // the negative child.
                    // --------------------------------------------------

                    forAll(badCell, cfI)
                    {
                        const face& f =
                            v1Faces[badCell[cfI]];

                        std::set<label> faceEdges;

                        forAll(f, fpI)
                        {
                            std::map<label,label>::const_iterator peIt =
                                v1PointEdge.find(f[fpI]);

                            if( peIt != v1PointEdge.end() )
                                faceEdges.insert(peIt->second);
                        }

                        for
                        (
                            std::set<label>::const_iterator
                                aIt=faceEdges.begin();
                            aIt!=faceEdges.end();
                            ++aIt
                        )
                        {
                            std::set<label>::const_iterator bIt = aIt;
                            ++bIt;

                            for
                            (
                                ;
                                bIt!=faceEdges.end();
                                ++bIt
                            )
                            {
                                pairEdges.insert
                                (
                                    std::make_pair
                                    (
                                        *aIt,
                                        *bIt
                                    )
                                );
                            }
                        }
                    }


                    // Keep the first paired search deliberately bounded.
                    // Both signs are tested for both hairs.
                    static const scalar v1DPairAmplitudeMag[] =
                    {
                        scalar(0.01),
                        scalar(0.02),
                        scalar(0.03),
                        scalar(0.04),
                        scalar(0.05),
                        scalar(0.06),
                        scalar(0.08),
                        scalar(0.10),
                        scalar(0.12),
                        scalar(0.16),
                        scalar(0.20),
                        scalar(0.25),
                        scalar(0.30),
                        scalar(0.35),
                        scalar(0.40),
                        scalar(0.425),
                        scalar(0.45),
                        scalar(0.475),
                        scalar(0.50)
                    };

                    static const label nV1DPairAmplitudeMag =
                        sizeof(v1DPairAmplitudeMag)
                       /sizeof(v1DPairAmplitudeMag[0]);

                    // At least one member of a pair must remain a
                    // compensating correction.  Either edge may be the
                    // large mover; no topology-specific ordering is used.
                    static const scalar v1DPairSecondaryMax =
                        scalar(0.20);


                    label v1DPairTrials = 0;
                    label v1DPairVolumePass = 0;
                    label v1DPairQualityPass = 0;

                    label v1DRejectNewPyramid = 0;
                    label v1DRejectWorsePyramid = 0;
                    label v1DRejectNewSkew = 0;
                    label v1DRejectWorseSkew = 0;
                    label v1DRejectNewNonOrtho90 = 0;
                    label v1DRejectWorseNonOrtho = 0;
                    label v1DRejectOFGeometry = 0;

                    // CFMitch v2.6:
                    // Pair candidate passed the ordinary union-star
                    // admissibility tests but failed to eliminate every
                    // bad pyramid side belonging to the quality target.
                    label v1DRejectTargetPyramid = 0;

                    label v1DRejectOther = 0;


                    // --------------------------------------------------
                    // Generic chain capture helper.
                    // --------------------------------------------------

                    auto capturePairChain =
                    [&]
                    (
                        const label edgeI,
                        point& edgeStart,
                        vector& edgeVec,
                        label& rowSize,
                        List<point>& originalPositions,
                        List<scalar>& originalT
                    ) -> bool
                    {
                        const edge& se =
                            splitEdges_[edgeI];

                        edgeStart =
                            v1Points[se.start()];

                        const point edgeEnd =
                            v1Points[se.end()];

                        edgeVec =
                            edgeEnd-edgeStart;

                        const scalar edgeMagSqr =
                            magSqr(edgeVec);

                        if( edgeMagSqr <= VSMALL )
                            return false;

                        rowSize =
                            newVerticesForSplitEdge_.
                                sizeOfRow(edgeI);

                        if( rowSize < 3 )
                            return false;

                        originalPositions.setSize
                        (
                            rowSize-2
                        );

                        originalT.setSize
                        (
                            rowSize
                        );

                        originalT[0] =
                            scalar(0);

                        originalT[rowSize-1] =
                            scalar(1);

                        for
                        (
                            label rowI=1;
                            rowI<rowSize-1;
                            ++rowI
                        )
                        {
                            const label pointI =
                                newVerticesForSplitEdge_
                                (
                                    edgeI,
                                    rowI
                                );

                            originalPositions[rowI-1] =
                                v1Points[pointI];

                            originalT[rowI] =
                                (
                                    (v1Points[pointI]-edgeStart)
                                  & edgeVec
                                )
                               /(edgeMagSqr + VSMALL);

                            if
                            (
                                !(originalT[rowI] > originalT[rowI-1])
                             || !(originalT[rowI] < scalar(1))
                            )
                                return false;
                        }

                        return true;
                    };


                    // --------------------------------------------------
                    // Generic coherent-chain amplitude helper.
                    //
                    // Caller owns restoration if this returns false.
                    // --------------------------------------------------

                    auto applyPairChain =
                    [&]
                    (
                        const label edgeI,
                        const point& edgeStart,
                        const vector& edgeVec,
                        const label rowSize,
                        const List<scalar>& originalT,
                        const scalar amplitude
                    ) -> bool
                    {
                        const scalar firstT =
                            originalT[1];

                        if
                        (
                            firstT <= scalar(0)
                         || firstT >= scalar(1)
                        )
                            return false;

                        scalar prevT =
                            scalar(0);

                        for
                        (
                            label rowI=1;
                            rowI<rowSize-1;
                            ++rowI
                        )
                        {
                            const scalar t =
                                originalT[rowI];

                            const scalar decayBase =
                                Foam::max
                                (
                                    scalar(0),
                                    (scalar(1)-t)
                                   /(scalar(1)-firstT + VSMALL)
                                );

                            const scalar deltaT =
                                amplitude
                               *firstT
                               *Foam::pow
                                (
                                    decayBase,
                                    scalar(4)
                                );

                            const scalar warpedT =
                                t + deltaT;

                            if
                            (
                                !(warpedT > prevT)
                             || !(warpedT < originalT[rowI+1])
                            )
                                return false;

                            const label pointI =
                                newVerticesForSplitEdge_
                                (
                                    edgeI,
                                    rowI
                                );

                            v1Points[pointI] =
                                edgeStart
                              + warpedT*edgeVec;

                            prevT =
                                warpedT;
                        }

                        return true;
                    };


                    for
                    (
                        std::set<std::pair<label,label> >::const_iterator
                            pairIt=pairEdges.begin();
                        pairIt!=pairEdges.end();
                        ++pairIt
                    )
                    {
                        const label edgeA =
                            pairIt->first;

                        const label edgeB =
                            pairIt->second;


                        std::map
                        <
                            label,
                            std::set<label>
                        >::const_iterator starAIt =
                            v1EdgeCells.find(edgeA);

                        std::map
                        <
                            label,
                            std::set<label>
                        >::const_iterator starBIt =
                            v1EdgeCells.find(edgeB);

                        if
                        (
                            starAIt == v1EdgeCells.end()
                         || starBIt == v1EdgeCells.end()
                         || starAIt->second.empty()
                         || starBIt->second.empty()
                        )
                            continue;


                        point edgeStartA(vector::zero);
                        point edgeStartB(vector::zero);

                        vector edgeVecA(vector::zero);
                        vector edgeVecB(vector::zero);

                        label rowSizeA = 0;
                        label rowSizeB = 0;

                        List<point> originalPositionsA;
                        List<point> originalPositionsB;

                        List<scalar> originalTA;
                        List<scalar> originalTB;


                        if
                        (
                            !capturePairChain
                            (
                                edgeA,
                                edgeStartA,
                                edgeVecA,
                                rowSizeA,
                                originalPositionsA,
                                originalTA
                            )
                         || !capturePairChain
                            (
                                edgeB,
                                edgeStartB,
                                edgeVecB,
                                rowSizeB,
                                originalPositionsB,
                                originalTB
                            )
                        )
                            continue;


                        auto restorePair =
                        [&]()
                        {
                            for
                            (
                                label rowI=1;
                                rowI<rowSizeA-1;
                                ++rowI
                            )
                            {
                                const label pointI =
                                    newVerticesForSplitEdge_
                                    (
                                        edgeA,
                                        rowI
                                    );

                                v1Points[pointI] =
                                    originalPositionsA[rowI-1];
                            }

                            for
                            (
                                label rowI=1;
                                rowI<rowSizeB-1;
                                ++rowI
                            )
                            {
                                const label pointI =
                                    newVerticesForSplitEdge_
                                    (
                                        edgeB,
                                        rowI
                                    );

                                v1Points[pointI] =
                                    originalPositionsB[rowI-1];
                            }
                        };


                        // ----------------------------------------------
                        // Union of both complete incident-cell stars.
                        // ----------------------------------------------

                        std::set<label> unionStar =
                            starAIt->second;

                        unionStar.insert
                        (
                            starBIt->second.begin(),
                            starBIt->second.end()
                        );


                        std::map<label,scalar>
                            pairBaselineVolumes;

                        std::map<label,scalar>
                            pairBaselineOFVolumes;


                        bool pairBaselineValid = true;

                        for
                        (
                            std::set<label>::const_iterator
                                cIt=unionStar.begin();
                            cIt!=unionStar.end();
                            ++cIt
                        )
                        {
                            pairBaselineVolumes[*cIt] =
                                v1CellVolume(*cIt);

                            point cc(vector::zero);
                            scalar cv = -GREAT;

                            if
                            (
                                !v1OFCellCentreVolume
                                (
                                    *cIt,
                                    cc,
                                    cv
                                )
                            )
                            {
                                pairBaselineValid = false;
                                break;
                            }

                            pairBaselineOFVolumes[*cIt] =
                                cv;
                        }

                        if( !pairBaselineValid )
                            continue;


                        // ----------------------------------------------
                        // All faces whose direct geometry OR owner/
                        // neighbour centre may change.
                        // ----------------------------------------------

                        std::set<label> pairAffectedFaces;

                        for
                        (
                            std::set<label>::const_iterator
                                cIt=unionStar.begin();
                            cIt!=unionStar.end();
                            ++cIt
                        )
                        {
                            const cell& starCell =
                                v1Cells[*cIt];

                            forAll(starCell, sfI)
                                pairAffectedFaces.insert(starCell[sfI]);
                        }


                        std::map<label,scalar>
                            pairBaselineOFSkew;

                        std::map<label,scalar>
                            pairBaselineOFOrtho;

                        std::map<label,scalar>
                            pairBaselineOFPyrMargin;

                        std::map<label,bool>
                            pairBaselineOFPyrBad;

                        std::map<label,bool>
                            pairBaselineOFEligible;


                        for
                        (
                            std::set<label>::const_iterator
                                fIt=pairAffectedFaces.begin();
                            fIt!=pairAffectedFaces.end();
                            ++fIt
                        )
                        {
                            const label faceI =
                                *fIt;

                            scalar skew = GREAT;
                            scalar ortho = -GREAT;
                            scalar pyrMargin = -GREAT;
                            scalar ownVol = -GREAT;
                            scalar neiVol = -GREAT;
                            bool pyrBad = true;

                            if
                            (
                                !v1OFFaceQuality
                                (
                                    faceI,
                                    skew,
                                    ortho,
                                    pyrMargin,
                                    pyrBad,
                                    ownVol,
                                    neiVol
                                )
                            )
                            {
                                pairBaselineValid = false;
                                break;
                            }

                            pairBaselineOFSkew[faceI] =
                                skew;

                            pairBaselineOFOrtho[faceI] =
                                ortho;

                            pairBaselineOFPyrMargin[faceI] =
                                pyrMargin;

                            pairBaselineOFPyrBad[faceI] =
                                pyrBad;

                            pairBaselineOFEligible[faceI] =
                                (
                                    ownVol > scalar(0)
                                 && (
                                        v1Neighbour[faceI] < 0
                                     || neiVol > scalar(0)
                                    )
                                );
                        }

                        if( !pairBaselineValid )
                            continue;


                        const scalar pairTargetFloor =
                            scalar(0.20)
                           *Foam::mag(badVolBefore);

                        static const scalar v1DPairSkewLimit =
                            scalar(4);

                        static const scalar v1DPairOrtho90 =
                            scalar(0);


                        for
                        (
                            label ampAI=0;
                            ampAI<nV1DPairAmplitudeMag;
                            ++ampAI
                        )
                        {
                            for
                            (
                                label signA=-1;
                                signA<=1;
                                signA+=2
                            )
                            {
                                const scalar amplitudeA =
                                    scalar(signA)
                                   *v1DPairAmplitudeMag[ampAI];


                                for
                                (
                                    label ampBI=0;
                                    ampBI<nV1DPairAmplitudeMag;
                                    ++ampBI
                                )
                                {
                                    const scalar pairMagA =
                                        v1DPairAmplitudeMag[ampAI];

                                    const scalar pairMagB =
                                        v1DPairAmplitudeMag[ampBI];

                                    if
                                    (
                                        pairMagA > v1DPairSecondaryMax
                                     && pairMagB > v1DPairSecondaryMax
                                    )
                                    {
                                        continue;
                                    }

                                    for
                                    (
                                        label signB=-1;
                                        signB<=1;
                                        signB+=2
                                    )
                                    {
                                        const scalar amplitudeB =
                                            scalar(signB)
                                           *v1DPairAmplitudeMag[ampBI];

                                        restorePair();


                                        if
                                        (
                                            !applyPairChain
                                            (
                                                edgeA,
                                                edgeStartA,
                                                edgeVecA,
                                                rowSizeA,
                                                originalTA,
                                                amplitudeA
                                            )
                                        )
                                        {
                                            restorePair();
                                            continue;
                                        }

                                        if
                                        (
                                            !applyPairChain
                                            (
                                                edgeB,
                                                edgeStartB,
                                                edgeVecB,
                                                rowSizeB,
                                                originalTB,
                                                amplitudeB
                                            )
                                        )
                                        {
                                            restorePair();
                                            continue;
                                        }


                                        ++v1DPairTrials;


                                        // ----------------------------------
                                        // Raw cfMesh signed-volume gate.
                                        // ----------------------------------

                                        const scalar targetVol =
                                            v1CellVolume(badCellI);

                                        if
                                        (
                                            targetVol
                                         <= pairTargetFloor
                                        )
                                            continue;


                                        bool safe = true;

                                        scalar minPositiveRatio =
                                            GREAT;


                                        for
                                        (
                                            std::set<label>::const_iterator
                                                cIt=unionStar.begin();
                                            cIt!=unionStar.end();
                                            ++cIt
                                        )
                                        {
                                            const label starCellI =
                                                *cIt;

                                            const scalar oldV =
                                                pairBaselineVolumes
                                                [
                                                    starCellI
                                                ];

                                            const scalar newV =
                                                v1CellVolume
                                                (
                                                    starCellI
                                                );

                                            if( oldV > scalar(0) )
                                            {
                                                if( newV <= scalar(0) )
                                                {
                                                    safe = false;
                                                    break;
                                                }

                                                const scalar ratio =
                                                    newV
                                                   /(oldV + VSMALL);

                                                minPositiveRatio =
                                                    Foam::min
                                                    (
                                                        minPositiveRatio,
                                                        ratio
                                                    );

                                                if
                                                (
                                                    ratio
                                                  < scalar(0.75)
                                                )
                                                {
                                                    safe = false;
                                                    break;
                                                }
                                            }
                                            else if( newV < oldV )
                                            {
                                                safe = false;
                                                break;
                                            }
                                        }

                                        if( !safe )
                                            continue;


                                        // ----------------------------------
                                        // OpenFOAM signed-volume gate.
                                        // ----------------------------------

                                        point targetOFCentre
                                        (
                                            vector::zero
                                        );

                                        scalar targetOFVol =
                                            -GREAT;

                                        if
                                        (
                                            !v1OFCellCentreVolume
                                            (
                                                badCellI,
                                                targetOFCentre,
                                                targetOFVol
                                            )
                                        )
                                            continue;

                                        if
                                        (
                                            targetOFVol
                                         <= pairTargetFloor
                                        )
                                            continue;


                                        scalar minOFPositiveRatio =
                                            GREAT;


                                        for
                                        (
                                            std::set<label>::const_iterator
                                                cIt=unionStar.begin();
                                            cIt!=unionStar.end();
                                            ++cIt
                                        )
                                        {
                                            const label starCellI =
                                                *cIt;

                                            const scalar oldOFV =
                                                pairBaselineOFVolumes
                                                [
                                                    starCellI
                                                ];

                                            point newOFCentre
                                            (
                                                vector::zero
                                            );

                                            scalar newOFV =
                                                -GREAT;

                                            if
                                            (
                                                !v1OFCellCentreVolume
                                                (
                                                    starCellI,
                                                    newOFCentre,
                                                    newOFV
                                                )
                                            )
                                            {
                                                safe = false;
                                                break;
                                            }

                                            if
                                            (
                                                oldOFV
                                              > scalar(0)
                                            )
                                            {
                                                if
                                                (
                                                    newOFV
                                                 <= scalar(0)
                                                )
                                                {
                                                    safe = false;
                                                    break;
                                                }

                                                const scalar ratio =
                                                    newOFV
                                                   /(oldOFV + VSMALL);

                                                minOFPositiveRatio =
                                                    Foam::min
                                                    (
                                                        minOFPositiveRatio,
                                                        ratio
                                                    );

                                                if
                                                (
                                                    ratio
                                                  < scalar(0.75)
                                                )
                                                {
                                                    safe = false;
                                                    break;
                                                }
                                            }
                                            else if
                                            (
                                                newOFV < oldOFV
                                            )
                                            {
                                                safe = false;
                                                break;
                                            }
                                        }

                                        if( !safe )
                                            continue;


                                        ++v1DPairVolumePass;


                                        // ----------------------------------
                                        // OpenFOAM-parity quality gate over
                                        // the ENTIRE union-star face set.
                                        // ----------------------------------

                                        bool ofSafe = true;

                                        scalar trialMaxSkew =
                                            scalar(0);

                                        scalar trialMinOrtho =
                                            GREAT;

                                        label trialBadPyrCount =
                                            0;

                                        label rejectFace =
                                            -1;

                                        const char* rejectReason =
                                            "none";


                                        for
                                        (
                                            std::set<label>::const_iterator
                                                fIt=
                                                    pairAffectedFaces.begin();
                                            fIt!=
                                                pairAffectedFaces.end();
                                            ++fIt
                                        )
                                        {
                                            const label faceI =
                                                *fIt;

                                            scalar trialSkew =
                                                GREAT;

                                            scalar trialOrtho =
                                                -GREAT;

                                            scalar trialPyrMargin =
                                                -GREAT;

                                            scalar trialOwnVol =
                                                -GREAT;

                                            scalar trialNeiVol =
                                                -GREAT;

                                            bool trialPyrBad =
                                                true;


                                            if
                                            (
                                                !v1OFFaceQuality
                                                (
                                                    faceI,
                                                    trialSkew,
                                                    trialOrtho,
                                                    trialPyrMargin,
                                                    trialPyrBad,
                                                    trialOwnVol,
                                                    trialNeiVol
                                                )
                                            )
                                            {
                                                ofSafe = false;
                                                rejectFace = faceI;
                                                rejectReason =
                                                    "ofGeometry";
                                                break;
                                            }


                                            const bool trialEligible =
                                                (
                                                    trialOwnVol
                                                  > scalar(0)
                                                 && (
                                                        v1Neighbour
                                                        [
                                                            faceI
                                                        ] < 0
                                                     || trialNeiVol
                                                      > scalar(0)
                                                    )
                                                );

                                            if( !trialEligible )
                                                continue;


                                            if( trialPyrBad )
                                                ++trialBadPyrCount;


                                            trialMaxSkew =
                                                Foam::max
                                                (
                                                    trialMaxSkew,
                                                    trialSkew
                                                );


                                            if
                                            (
                                                v1Neighbour[faceI]
                                              >= 0
                                            )
                                            {
                                                trialMinOrtho =
                                                    Foam::min
                                                    (
                                                        trialMinOrtho,
                                                        trialOrtho
                                                    );
                                            }


                                            const bool baselineEligible =
                                                pairBaselineOFEligible
                                                [
                                                    faceI
                                                ];


                                            if( !baselineEligible )
                                            {
                                                if( trialPyrBad )
                                                {
                                                    ofSafe = false;
                                                    rejectFace = faceI;
                                                    rejectReason =
                                                        "newPyramid";
                                                    break;
                                                }

                                                if
                                                (
                                                    trialSkew
                                                  > v1DPairSkewLimit
                                                )
                                                {
                                                    ofSafe = false;
                                                    rejectFace = faceI;
                                                    rejectReason =
                                                        "newSkew";
                                                    break;
                                                }

                                                if
                                                (
                                                    v1Neighbour[faceI]
                                                  >= 0
                                                 && trialOrtho
                                                  < v1DPairOrtho90
                                                )
                                                {
                                                    ofSafe = false;
                                                    rejectFace = faceI;
                                                    rejectReason =
                                                        "newNonOrtho90";
                                                    break;
                                                }

                                                continue;
                                            }


                                            const bool basePyrBad =
                                                pairBaselineOFPyrBad
                                                [
                                                    faceI
                                                ];

                                            const scalar basePyrMargin =
                                                pairBaselineOFPyrMargin
                                                [
                                                    faceI
                                                ];


                                            if
                                            (
                                                !basePyrBad
                                             && trialPyrBad
                                            )
                                            {
                                                ofSafe = false;
                                                rejectFace = faceI;
                                                rejectReason =
                                                    "newPyramid";
                                                break;
                                            }


                                            if
                                            (
                                                basePyrBad
                                             && trialPyrBad
                                            )
                                            {
                                                const scalar pyrTol =
                                                    scalar(1e-12)
                                                   *(
                                                        Foam::mag
                                                        (
                                                            basePyrMargin
                                                        )
                                                      + SMALL
                                                    );

                                                if
                                                (
                                                    trialPyrMargin
                                                  < basePyrMargin-pyrTol
                                                )
                                                {
                                                    ofSafe = false;
                                                    rejectFace = faceI;
                                                    rejectReason =
                                                        "worsePyramid";
                                                    break;
                                                }
                                            }


                                            const scalar baseSkew =
                                                pairBaselineOFSkew
                                                [
                                                    faceI
                                                ];

                                            const scalar skewTol =
                                                scalar(1e-10)
                                               *(
                                                    scalar(1)
                                                  + Foam::mag(baseSkew)
                                                );


                                            if
                                            (
                                                baseSkew
                                              <= v1DPairSkewLimit
                                            )
                                            {
                                                if
                                                (
                                                    trialSkew
                                                  > v1DPairSkewLimit
                                                )
                                                {
                                                    ofSafe = false;
                                                    rejectFace = faceI;
                                                    rejectReason =
                                                        "newSkew";
                                                    break;
                                                }
                                            }
                                            else if
                                            (
                                                trialSkew
                                              > baseSkew+skewTol
                                            )
                                            {
                                                ofSafe = false;
                                                rejectFace = faceI;
                                                rejectReason =
                                                    "worseSkew";
                                                break;
                                            }


                                            if
                                            (
                                                v1Neighbour[faceI]
                                              >= 0
                                            )
                                            {
                                                const scalar baseOrtho =
                                                    pairBaselineOFOrtho
                                                    [
                                                        faceI
                                                    ];

                                                const scalar orthoTol =
                                                    scalar(1e-12);


                                                if
                                                (
                                                    baseOrtho
                                                  >= v1DPairOrtho90
                                                )
                                                {
                                                    if
                                                    (
                                                        trialOrtho
                                                      < v1DPairOrtho90
                                                    )
                                                    {
                                                        ofSafe = false;
                                                        rejectFace = faceI;
                                                        rejectReason =
                                                            "newNonOrtho90";
                                                        break;
                                                    }
                                                }
                                                else if
                                                (
                                                    trialOrtho
                                                  < baseOrtho-orthoTol
                                                )
                                                {
                                                    ofSafe = false;
                                                    rejectFace = faceI;
                                                    rejectReason =
                                                        "worseNonOrtho";
                                                    break;
                                                }
                                            }
                                        }


                                        if( !ofSafe )
                                        {
                                            if
                                            (
                                                std::strcmp
                                                (
                                                    rejectReason,
                                                    "newPyramid"
                                                ) == 0
                                            )
                                            {
                                                ++v1DRejectNewPyramid;
                                            }
                                            else if
                                            (
                                                std::strcmp
                                                (
                                                    rejectReason,
                                                    "worsePyramid"
                                                ) == 0
                                            )
                                            {
                                                ++v1DRejectWorsePyramid;
                                            }
                                            else if
                                            (
                                                std::strcmp
                                                (
                                                    rejectReason,
                                                    "newSkew"
                                                ) == 0
                                            )
                                            {
                                                ++v1DRejectNewSkew;
                                            }
                                            else if
                                            (
                                                std::strcmp
                                                (
                                                    rejectReason,
                                                    "worseSkew"
                                                ) == 0
                                            )
                                            {
                                                ++v1DRejectWorseSkew;
                                            }
                                            else if
                                            (
                                                std::strcmp
                                                (
                                                    rejectReason,
                                                    "newNonOrtho90"
                                                ) == 0
                                            )
                                            {
                                                ++v1DRejectNewNonOrtho90;
                                            }
                                            else if
                                            (
                                                std::strcmp
                                                (
                                                    rejectReason,
                                                    "worseNonOrtho"
                                                ) == 0
                                            )
                                            {
                                                ++v1DRejectWorseNonOrtho;
                                            }
                                            else if
                                            (
                                                std::strcmp
                                                (
                                                    rejectReason,
                                                    "ofGeometry"
                                                ) == 0
                                            )
                                            {
                                                ++v1DRejectOFGeometry;
                                            }
                                            else
                                            {
                                                ++v1DRejectOther;
                                            }

                                            continue;
                                        }


                                        // CFMitch v2.6:
                                        // For a positive-volume quality
                                        // seed, the target pyramid defect
                                        // itself must be completely removed.
                                        //
                                        // The ordinary V1D face gate above
                                        // protects the full union star from
                                        // NEW or WORSE quality.  This extra
                                        // predicate makes the repair
                                        // objective explicit for the target
                                        // cell rather than merely allowing
                                        // its existing bad pyramid to remain.
                                        if( requirePyramidRepair )
                                        {
                                            const label
                                                trialTargetBadPyr =
                                                    v1CellBadPyramidCount
                                                    (
                                                        badCellI
                                                    );

                                            if( trialTargetBadPyr != 0 )
                                            {
                                                ++v1DRejectTargetPyramid;

                                                // Aggregate with the v2.5
                                                // single-hair target reject
                                                // count for the final V1C
                                                // summary as well.
                                                ++blV25TargetPyrRejects;

                                                continue;
                                            }
                                        }


                                        ++v1DPairQualityPass;


                                        // ----------------------------------
                                        // Candidate motion magnitude.
                                        // ----------------------------------

                                        scalar trialMaxMove =
                                            scalar(0);


                                        for
                                        (
                                            label rowI=1;
                                            rowI<rowSizeA-1;
                                            ++rowI
                                        )
                                        {
                                            const label pointI =
                                                newVerticesForSplitEdge_
                                                (
                                                    edgeA,
                                                    rowI
                                                );

                                            trialMaxMove =
                                                Foam::max
                                                (
                                                    trialMaxMove,
                                                    mag
                                                    (
                                                        v1Points[pointI]
                                                      - originalPositionsA
                                                        [
                                                            rowI-1
                                                        ]
                                                    )
                                                );
                                        }


                                        for
                                        (
                                            label rowI=1;
                                            rowI<rowSizeB-1;
                                            ++rowI
                                        )
                                        {
                                            const label pointI =
                                                newVerticesForSplitEdge_
                                                (
                                                    edgeB,
                                                    rowI
                                                );

                                            trialMaxMove =
                                                Foam::max
                                                (
                                                    trialMaxMove,
                                                    mag
                                                    (
                                                        v1Points[pointI]
                                                      - originalPositionsB
                                                        [
                                                            rowI-1
                                                        ]
                                                    )
                                                );
                                        }


                                        // Same deterministic ranking as
                                        // v1c, now over paired candidates.
                                        bool better =
                                            !found;


                                        if( found )
                                        {
                                            if
                                            (
                                                trialBadPyrCount
                                              < bestBadPyrCount
                                            )
                                            {
                                                better = true;
                                            }
                                            else if
                                            (
                                                trialBadPyrCount
                                             == bestBadPyrCount
                                            )
                                            {
                                                if
                                                (
                                                    trialMaxSkew
                                                  < bestMaxOFSkew
                                                   - scalar(1e-12)
                                                )
                                                {
                                                    better = true;
                                                }
                                                else if
                                                (
                                                    Foam::mag
                                                    (
                                                        trialMaxSkew
                                                      - bestMaxOFSkew
                                                    )
                                                  <= scalar(1e-12)
                                                )
                                                {
                                                    if
                                                    (
                                                        trialMinOrtho
                                                      > bestMinOFOrtho
                                                       + scalar(1e-12)
                                                    )
                                                    {
                                                        better = true;
                                                    }
                                                    else if
                                                    (
                                                        Foam::mag
                                                        (
                                                            trialMinOrtho
                                                          - bestMinOFOrtho
                                                        )
                                                      <= scalar(1e-12)
                                                    )
                                                    {
                                                        if
                                                        (
                                                            trialMaxMove
                                                          < bestMaxMove
                                                           - scalar(1e-15)
                                                        )
                                                        {
                                                            better = true;
                                                        }
                                                        else if
                                                        (
                                                            Foam::mag
                                                            (
                                                                trialMaxMove
                                                              - bestMaxMove
                                                            )
                                                          <= scalar(1e-15)
                                                         && minPositiveRatio
                                                          > bestMinRatio
                                                        )
                                                        {
                                                            better = true;
                                                        }
                                                    }
                                                }
                                            }
                                        }


                                        if( better )
                                        {
                                            found = true;
                                            bestIsPair = true;

                                            bestEdgeA =
                                                edgeA;

                                            bestEdgeB =
                                                edgeB;

                                            bestAmplitudeA =
                                                amplitudeA;

                                            bestAmplitudeB =
                                                amplitudeB;

                                            bestTargetVol =
                                                targetVol;

                                            bestTargetOFVol =
                                                targetOFVol;

                                            bestMinRatio =
                                                minPositiveRatio;

                                            bestMinOFRatio =
                                                minOFPositiveRatio;

                                            bestMaxOFSkew =
                                                trialMaxSkew;

                                            bestMinOFOrtho =
                                                trialMinOrtho;

                                            bestBadPyrCount =
                                                trialBadPyrCount;

                                            bestMaxMove =
                                                trialMaxMove;


                                            bestPositionsA.setSize
                                            (
                                                rowSizeA-2
                                            );

                                            for
                                            (
                                                label rowI=1;
                                                rowI<rowSizeA-1;
                                                ++rowI
                                            )
                                            {
                                                const label pointI =
                                                    newVerticesForSplitEdge_
                                                    (
                                                        edgeA,
                                                        rowI
                                                    );

                                                bestPositionsA[rowI-1] =
                                                    v1Points[pointI];
                                            }


                                            bestPositionsB.setSize
                                            (
                                                rowSizeB-2
                                            );

                                            for
                                            (
                                                label rowI=1;
                                                rowI<rowSizeB-1;
                                                ++rowI
                                            )
                                            {
                                                const label pointI =
                                                    newVerticesForSplitEdge_
                                                    (
                                                        edgeB,
                                                        rowI
                                                    );

                                                bestPositionsB[rowI-1] =
                                                    v1Points[pointI];
                                            }


                                            Info
                                                << "BL_VALIDITY_REPAIR_V1D_PAIR_CANDIDATE"
                                                << " cell="
                                                << badCellI
                                                << " parent="
                                                << exactVolumeParent
                                                   [
                                                       badCellI
                                                   ]
                                                << " localChild="
                                                << exactVolumeLocalChild
                                                   [
                                                       badCellI
                                                   ]
                                                << " edgeA="
                                                << edgeA
                                                << " amplitudeA="
                                                << amplitudeA
                                                << " edgeB="
                                                << edgeB
                                                << " amplitudeB="
                                                << amplitudeB
                                                << " targetVol="
                                                << targetVol
                                                << " targetOFVol="
                                                << targetOFVol
                                                << " maxOFSkew="
                                                << trialMaxSkew
                                                << " minOFOrtho="
                                                << trialMinOrtho
                                                << " badOFPyr="
                                                << trialBadPyrCount
                                                << " maxPhysicalMove="
                                                << trialMaxMove
                                                << " minPositiveStarRatio="
                                                << minPositiveRatio
                                                << " minOFPositiveStarRatio="
                                                << minOFPositiveRatio
                                                << endl;
                                        }
                                    }
                                }
                            }
                        }


                        restorePair();
                    }


                    Info
                        << "BL_VALIDITY_REPAIR_V1D_PAIR_SUMMARY"
                        << " cell=" << badCellI
                        << " parent="
                        << exactVolumeParent[badCellI]
                        << " localChild="
                        << exactVolumeLocalChild[badCellI]
                        << " candidatePairs="
                        << pairEdges.size()
                        << " trials="
                        << v1DPairTrials
                        << " volumePass="
                        << v1DPairVolumePass
                        << " qualityPass="
                        << v1DPairQualityPass
                        << " rejectNewPyramid="
                        << v1DRejectNewPyramid
                        << " rejectWorsePyramid="
                        << v1DRejectWorsePyramid
                        << " rejectNewSkew="
                        << v1DRejectNewSkew
                        << " rejectWorseSkew="
                        << v1DRejectWorseSkew
                        << " rejectNewNonOrtho90="
                        << v1DRejectNewNonOrtho90
                        << " rejectWorseNonOrtho="
                        << v1DRejectWorseNonOrtho
                        << " rejectOFGeometry="
                        << v1DRejectOFGeometry
                        << " rejectTargetPyramid="
                        << v1DRejectTargetPyramid
                        << " rejectOther="
                        << v1DRejectOther
                        << " found="
                        << (bestIsPair ? "true" : "false")
                        << endl;
                }


                if( found && bestIsPair )
                {
                    const label rowSizeA =
                        newVerticesForSplitEdge_.
                            sizeOfRow(bestEdgeA);

                    const label rowSizeB =
                        newVerticesForSplitEdge_.
                            sizeOfRow(bestEdgeB);


                    scalar maxPhysicalMove =
                        scalar(0);


                    for
                    (
                        label rowI=1;
                        rowI<rowSizeA-1;
                        ++rowI
                    )
                    {
                        const label pointI =
                            newVerticesForSplitEdge_
                            (
                                bestEdgeA,
                                rowI
                            );

                        maxPhysicalMove =
                            Foam::max
                            (
                                maxPhysicalMove,
                                mag
                                (
                                    bestPositionsA[rowI-1]
                                  - v1Points[pointI]
                                )
                            );

                        v1Points[pointI] =
                            bestPositionsA[rowI-1];
                    }


                    for
                    (
                        label rowI=1;
                        rowI<rowSizeB-1;
                        ++rowI
                    )
                    {
                        const label pointI =
                            newVerticesForSplitEdge_
                            (
                                bestEdgeB,
                                rowI
                            );

                        maxPhysicalMove =
                            Foam::max
                            (
                                maxPhysicalMove,
                                mag
                                (
                                    bestPositionsB[rowI-1]
                                  - v1Points[pointI]
                                )
                            );

                        v1Points[pointI] =
                            bestPositionsB[rowI-1];
                    }


                    const scalar committedVol =
                        v1CellVolume(badCellI);


                    Info
                        << "BL_VALIDITY_REPAIR_V1D_PAIR_MOVE"
                        << " cell=" << badCellI
                        << " parent="
                        << exactVolumeParent[badCellI]
                        << " localChild="
                        << exactVolumeLocalChild[badCellI]
                        << " bfI="
                        << cellToBaseBndFace_[badCellI]
                        << " edgeA="
                        << bestEdgeA
                        << " amplitudeA="
                        << bestAmplitudeA
                        << " edgeB="
                        << bestEdgeB
                        << " amplitudeB="
                        << bestAmplitudeB
                        << " maxPhysicalMove="
                        << maxPhysicalMove
                        << " oldCellVol="
                        << badVolBefore
                        << " targetVol="
                        << bestTargetVol
                        << " committedVol="
                        << committedVol
                        << " minPositiveStarRatio="
                        << bestMinRatio
                        << " targetOFVol="
                        << bestTargetOFVol
                        << " maxOFSkew="
                        << bestMaxOFSkew
                        << " minOFOrtho="
                        << bestMinOFOrtho
                        << " badOFPyr="
                        << bestBadPyrCount
                        << " minOFPositiveStarRatio="
                        << bestMinOFRatio
                        << endl;


                    if( committedVol > scalar(0) )
                    {
                        ++blV1Fixed;

                        // This transaction commits two coherent chains.
                        blV1CommittedChains += 2;
                    }
                    else
                    {
                        ++blV1Unresolved;
                    }
                }
                else if( found )
                {
                    const label rowSize =
                        newVerticesForSplitEdge_.
                            sizeOfRow(bestEdge);

                    scalar maxPhysicalMove = scalar(0);

                    for
                    (
                        label rowI=1;
                        rowI<rowSize-1;
                        ++rowI
                    )
                    {
                        const label pointI =
                            newVerticesForSplitEdge_
                            (
                                bestEdge,
                                rowI
                            );

                        maxPhysicalMove =
                            Foam::max
                            (
                                maxPhysicalMove,
                                mag
                                (
                                    bestPositions[rowI-1]
                                  - v1Points[pointI]
                                )
                            );

                        v1Points[pointI] =
                            bestPositions[rowI-1];
                    }


                    const scalar committedVol =
                        v1CellVolume(badCellI);


                    Info
                        << "BL_VALIDITY_REPAIR_V1C_MOVE"
                        << " cell=" << badCellI
                        << " parent="
                        << exactVolumeParent[badCellI]
                        << " localChild="
                        << exactVolumeLocalChild[badCellI]
                        << " bfI="
                        << cellToBaseBndFace_[badCellI]
                        << " splitEdge="
                        << bestEdge
                        << " amplitude="
                        << bestAmplitude
                        << " maxPhysicalMove="
                        << maxPhysicalMove
                        << " oldCellVol="
                        << badVolBefore
                        << " targetVol="
                        << bestTargetVol
                        << " committedVol="
                        << committedVol
                        << " minPositiveStarRatio="
                        << bestMinRatio
                        << " targetOFVol="
                        << bestTargetOFVol
                        << " maxOFSkew="
                        << bestMaxOFSkew
                        << " minOFOrtho="
                        << bestMinOFOrtho
                        << " badOFPyr="
                        << bestBadPyrCount
                        << " minOFPositiveStarRatio="
                        << bestMinOFRatio
                        << endl;


                    if( committedVol > scalar(0) )
                    {
                        ++blV1Fixed;
                        ++blV1CommittedChains;
                    }
                    else
                    {
                        ++blV1Unresolved;
                    }
                }
                else
                {
                    ++blV1Unresolved;

                    // CFMitch v2.7 diagnostic:
                    // classify the actual OpenFOAM-bad face(s) only when
                    // this is an unresolved positive-pyramid type-1 child
                    // at wallLayer zero.
                    if
                    (
                        requirePyramidRepair
                     && badCellI >= 0
                     && badCellI <
                        label(exactVolumeRefType.size())
                     && exactVolumeRefType[badCellI] == 1
                     && badCellI <
                        label(exactVolumeLocalChild.size())
                    )
                    {
                        const label diagnosticBfI =
                            cellToBaseBndFace_[badCellI];

                        if
                        (
                            diagnosticBfI >= 0
                         && diagnosticBfI <
                            label(nLayersAtBndFace_.size())
                        )
                        {
                            const label diagnosticN =
                                nLayersAtBndFace_
                                [
                                    diagnosticBfI
                                ];

                            const label diagnosticLocalChild =
                                exactVolumeLocalChild
                                [
                                    badCellI
                                ];

                            const label diagnosticWallLayer =
                                diagnosticN
                              - 1
                              - diagnosticLocalChild;

                            if( diagnosticWallLayer == 0 )
                            {
                                v1ReportWallChildBadFaceRoles
                                (
                                    badCellI
                                );
                            }
                        }
                    }

                    // CFMitch v2.7.1:
                    //
                    // v2.7 quality-derived N retreat is intentionally
                    // disabled.
                    //
                    // Source forensics established that type-1 localChild
                    // numbering runs from the core side toward the wall:
                    //
                    //     localChild = 0      -> outer/core-side child
                    //     localChild = N - 1  -> wall-adjacent child
                    //
                    // Therefore the former
                    //
                    //     proposedCap = localChild
                    //
                    // did not represent a safe wall-side prefix and must
                    // not generate planner constraints.
                    //
                    // qualityMaxLayersAtFace_ remains available for a future
                    // correctly-defined front/termination policy.

                    Info
                        << "BL_VALIDITY_REPAIR_V1C_UNRESOLVED"
                        << " cell=" << badCellI
                        << " parent="
                        << exactVolumeParent[badCellI]
                        << " localChild="
                        << exactVolumeLocalChild[badCellI]
                        << " bfI="
                        << cellToBaseBndFace_[badCellI]
                        << " volume="
                        << badVolBefore
                        << " candidateEdges="
                        << badEdges.size()
                        << endl;
                }
            }


            Info
                << "CFMITCH V2.7.1 WALLFACE REPAIR SUMMARY:"
                << " eligible="
                << blV271FaceBreathEligible
                << " attempted="
                << blV271FaceBreathAttempted
                << " skipped="
                << blV271FaceBreathSkipped
                << " trials="
                << blV271FaceBreathTrials
                << " volumePass="
                << blV271FaceBreathVolumePass
                << " qualityPass="
                << blV271FaceBreathQualityPass
                << " rejectQuality="
                << blV271FaceBreathQualityReject
                << " rejectTargetPyramid="
                << blV271FaceBreathTargetReject
                << " fixed="
                << blV271FaceBreathFixed
                << " rejectNewPyramid="
                << blV271RejectNewPyramid
                << " rejectWorsePyramid="
                << blV271RejectWorsePyramid
                << " rejectNewSkew="
                << blV271RejectNewSkew
                << " rejectWorseSkew="
                << blV271RejectWorseSkew
                << " rejectNewNonOrtho90="
                << blV271RejectNewNonOrtho90
                << " rejectWorseNonOrtho="
                << blV271RejectWorseNonOrtho
                << " rejectOFGeometry="
                << blV271RejectOFGeometry
                << " rejectOther="
                << blV271RejectOther
                << " rejectOnTargetFace="
                << blV271RejectOnTargetCellFace
                << " rejectOnOtherStarFace="
                << blV271RejectOnOtherStarFace
                << endl;


            if( blV1CommittedChains )
            {
                // Coordinates changed; force all normal downstream checks
                // to rebuild exact geometry from the committed positions.
                mesh_.clearAddressingData();
            }
        }


        Info
            << "BL_VALIDITY_REPAIR_V1C_SUMMARY"
            << " initialNegative="
            << blV1InitialNegative
            << " pyramidAdditional="
            << blV25PyramidAdditionalSeeds
            << " repairSeeds="
            << blV25RepairSeeds
            << " fixed="
            << blV1Fixed
            << " unresolved="
            << blV1Unresolved
            << " committedChains="
            << blV1CommittedChains
            << " targetPyrRejects="
            << blV25TargetPyrRejects
            << endl;


        label blV27ZeroLayerFaces = 0;
        label blV27MinCap = labelMax;
        label blV27MaxCap = -1;

        forAllConstIter
        (
            Map<label>,
            qualityMaxLayersAtFace_,
            qIt
        )
        {
            const label cap = qIt();

            blV27MinCap =
                Foam::min(blV27MinCap, cap);

            blV27MaxCap =
                Foam::max(blV27MaxCap, cap);

            if( cap == 0 )
                ++blV27ZeroLayerFaces;
        }

        if( qualityMaxLayersAtFace_.empty() )
            blV27MinCap = -1;

        Info
            << "CFMITCH V2.7 QUALITY LAYER CAPS:"
            << " faces="
            << qualityMaxLayersAtFace_.size()
            << " zeroLayerFaces="
            << blV27ZeroLayerFaces
            << " minCap="
            << blV27MinCap
            << " maxCap="
            << blV27MaxCap
            << endl;

        Info
            << "CFMITCH V2.7 WALLCHILD FACE SUMMARY:"
            << " cells="
            << blV27WallChildCells
            << " badFaces="
            << blV27WallChildBadFaces
            << " wallBase="
            << blV27WallChildWallBaseBad
            << " internalOuter="
            << blV27WallChildInternalOuterBad
            << " lateral="
            << blV27WallChildLateralBad
            << " unknown="
            << blV27WallChildUnknownBad
            << " multiBadCells="
            << blV27WallChildMultiBad
            << endl;

    }

    # ifdef DEBUGLayer
    for(label procI=0;procI<Pstream::nProcs();++procI)
    {
        if( procI == Pstream::myProcNo() )
        {
            forAll(cells, cellI)
            {
                const cell& c = cells[cellI];

                DynList<edge> edges;
                DynList<label> nAppearances;
                forAll(c, fI)
                {
                    const face& f = faces[c[fI]];

                    forAll(f, eI)
                    {
                        const edge e = f.faceEdge(eI);

                        const label pos = edges.containsAtPosition(e);

                        if( pos < 0 )
                        {
                            edges.append(e);
                            nAppearances.append(1);
                        }
                        else
                        {
                            ++nAppearances[pos];
                        }
                    }
                }

                bool badCell(false);
                forAll(nAppearances, i)
                    if( nAppearances[i] != 2 )
                    {
                        badCell = true;
                        break;

                    }

                if( badCell )
                {
                    Pout << "Cell " << cellI
                         << " is not topologically closed" << endl;

                    forAll(c, fI)
                        Pout << "Face " << c[fI] << " with points "
                             << faces[c[fI]] << endl;
                    Pout << "Cell edges " << edges << endl;
                    Pout << "nAppearances " << nAppearances << endl;
                    ::exit(1);
                }
            }
        }

        returnReduce(1, sumOp<label>());
    }

    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();
    const label nInternalFaces = mesh_.nInternalFaces();

    for(label procI=0;procI<Pstream::nProcs();++procI)
    {
        if( procI == Pstream::myProcNo() )
        {
            forAll(faces, faceI)
            {
                if( faceI < nInternalFaces && neighbour[faceI] < 0 )
                {
                    Pout << "Num interface faces " << nInternalFaces
                         << " current face " << faceI
                         << " face points " << faces[faceI] << endl;
                    ::exit(1);
                }
                Pout << "Face " << faceI << " owner " << owner[faceI]
                     << " neighbour " << neighbour[faceI]
                     << " face points " << faces[faceI] << endl;
            }

            forAll(cells, cellI)
                Pout << "Cell " << cellI << " has faces "
                     << cells[cellI] << endl;
        }

        returnReduce(procI, maxOp<label>());
    }
    # endif

    Info << "Finished generating new cells " << endl;

    Info
        << "BL_CHILD_SWEEP_SUMMARY"
        << " triParentsChecked="
        << nChildSweepParentsChecked
        << " childrenChecked="
        << nChildSweepChildrenChecked
        << " badChildren="
        << nChildSweepBadChildren
        << " parentUnsafe="
        << nChildSweepParentUnsafe
        << " skippedParents="
        << nChildSweepSkippedParents
        << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
