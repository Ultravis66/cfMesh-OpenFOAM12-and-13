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

\*----------------------p-----------------------------------------------------*/

#include "decomposeCells.H"
#include "helperFunctions.H"
#include "triFace.H"
#include "polyMeshGenAddressing.H"

//#define DEBUGDecompose

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void decomposeCells::findAddressingForCell
(
    const label cellI,
    DynList<label, 32>& vrt,
    DynList<edge, 64>& edges,
    DynList<DynList<label, 8> >& faceEdges,
    DynList<DynList<label, 2>, 64>& edgeFaces
) const
{
    const cell& c = mesh_.cells()[cellI];

    vrt.clear();
    edges.clear();
    edgeFaces.clear();
    faceEdges.setSize(c.size());

    const faceListPMG& faces = mesh_.faces();
    forAll(faceEdges, feI)
    {
        faceEdges[feI].setSize(faces[c[feI]].size());
        faceEdges[feI] = -1;
    }

    forAll(c, fI)
    {
        const face& f = faces[c[fI]];

        forAll(f, eI)
        {
            const edge e = f.faceEdge(eI);

            bool store(true);
            forAll(vrt, vI)
                if( vrt[vI] == f[eI] )
                {
                    store = false;
                    break;
                }
            if( store )
            {
                vrt.append(f[eI]);
            }

            //- check if the edge alreready exists
            store = true;

            forAll(edges, eJ)
                if( e == edges[eJ] )
                {
                    store = false;
                    faceEdges[fI][eI] = eJ;
                    edgeFaces[eJ].append(fI);
                    break;
                }

            if( store )
            {
                faceEdges[fI][eI] = edges.size();
                DynList<label, 2> ef;
                ef.append(fI);
                edgeFaces.append(ef);
                edges.append(e);
            }
        }
    }

    //- check if the cell is topologically closed
    forAll(edgeFaces, efI)
        if( edgeFaces[efI].size() != 2 )
        {
            forAll(c, fI)
                Info << "Face " << c[fI] << " is " << faces[c[fI]] << endl;

            Info << "Edges " << edges << endl;
            Info << "faceEdges " << faceEdges << endl;
            Info << "edgeFaces " << edgeFaces << endl;
            mesh_.write();
            FatalErrorIn
            (
                "void decomposeCells::findAddressingForCell"
                "(const label, DynList<label, 32>&, DynList<edge, 64>&"
                ", DynList<DynList<label, 8> >&"
                ", DynList<DynList<label, 2>, 64>&) const"
            ) << "Cell " << cellI << " is not topologically closed!"
                << abort(FatalError);
        }
}

bool decomposeCells::findValidPyramidApex
(
    const label cellI,
    point& apex,
    scalar& relativeMargin
) const
{
    const cell& c = mesh_.cells()[cellI];
    const faceListPMG& faces = mesh_.faces();
    const labelList& owner = mesh_.owner();
    const pointFieldPMG& points = mesh_.points();

    const labelList cp = c.labels(faces);

    relativeMargin = -GREAT;
    apex = point::zero;

    if( cp.empty() )
        return false;

    // Historical cfMesh starting point.
    point seed(point::zero);

    forAll(cp, cpI)
        seed += points[cp[cpI]];

    seed /= scalar(cp.size());

    // Cell-relative length scale.
    scalar localScale = VSMALL;

    forAll(cp, cpI)
    {
        localScale = Foam::max
        (
            localScale,
            mag(points[cp[cpI]] - seed)
        );
    }

    if( localScale <= VSMALL )
        return false;

    const scalar feasibilityTol = Foam::max
    (
        scalar(1e-12)*localScale,
        scalar(100)*VSMALL
    );

    // Returns true when x can be projected into the intersection of
    // all normalized face-triangle half-spaces with the requested
    // positive interior margin.
    //
    // For owner faces the stored face orientation is outward from
    // this cell. For neighbour faces it is reversed.
    auto projectForMargin =
    [&]
    (
        const scalar targetMargin,
        point& x
    ) -> bool
    {
        const label maxProjectionPasses = 160;

        for
        (
            label pass=0;
            pass<maxProjectionPasses;
            ++pass
        )
        {
            scalar maxDeficit = 0.0;

            forAll(c, cfI)
            {
                const label faceI = c[cfI];
                const face& f = faces[faceI];

                if( f.size() < 3 )
                    return false;

                const bool cellIsOwner =
                    (owner[faceI] == cellI);

                const point& p0 = points[f[0]];

                // Fan-triangle constraints intentionally make this
                // conservative for warped/non-planar polygon faces.
                for
                (
                    label fpI=1;
                    fpI<f.size()-1;
                    ++fpI
                )
                {
                    const point& p1 = points[f[fpI]];
                    const point& p2 = points[f[fpI+1]];

                    vector outward =
                        (p1 - p0) ^ (p2 - p0);

                    if( !cellIsOwner )
                        outward = -outward;

                    const scalar nMag = mag(outward);

                    if( nMag <= VSMALL )
                        return false;

                    outward /= nMag;

                    const scalar insideDistance =
                        -((x - p0) & outward);

                    const scalar deficit =
                        targetMargin - insideDistance;

                    if( deficit > 0.0 )
                    {
                        // Orthogonal projection into this half-space.
                        x -= deficit*outward;

                        maxDeficit =
                            Foam::max(maxDeficit, deficit);
                    }
                }
            }

            if( maxDeficit <= feasibilityTol )
                break;
        }

        // Independent final verification. The projection loop itself
        // is never trusted as proof of feasibility.
        scalar minInsideDistance = GREAT;
        label nConstraints = 0;

        forAll(c, cfI)
        {
            const label faceI = c[cfI];
            const face& f = faces[faceI];

            if( f.size() < 3 )
                return false;

            const bool cellIsOwner =
                (owner[faceI] == cellI);

            const point& p0 = points[f[0]];

            for
            (
                label fpI=1;
                fpI<f.size()-1;
                ++fpI
            )
            {
                const point& p1 = points[f[fpI]];
                const point& p2 = points[f[fpI+1]];

                vector outward =
                    (p1 - p0) ^ (p2 - p0);

                if( !cellIsOwner )
                    outward = -outward;

                const scalar nMag = mag(outward);

                if( nMag <= VSMALL )
                    return false;

                outward /= nMag;

                const scalar insideDistance =
                    -((x - p0) & outward);

                minInsideDistance =
                    Foam::min
                    (
                        minInsideDistance,
                        insideDistance
                    );

                ++nConstraints;
            }
        }

        if( nConstraints == 0 )
            return false;

        return
            minInsideDistance
         >= targetMargin - feasibilityTol;
    };


    // First prove that a zero-margin kernel exists at all.
    point feasiblePoint(seed);

    if( !projectForMargin(0.0, feasiblePoint) )
    {
        apex = feasiblePoint;
        relativeMargin = -1.0;
        return false;
    }


    // Approximate the Chebyshev-centre margin by binary searching the
    // largest common positive offset for which all face half-spaces
    // remain feasible.
    //
    // localScale is a conservative upper bound for useful clearance.
    scalar lowerMargin = 0.0;
    scalar upperMargin = localScale;

    point bestPoint(feasiblePoint);

    const label nMarginIterations = 28;

    for
    (
        label iter=0;
        iter<nMarginIterations;
        ++iter
    )
    {
        const scalar trialMargin =
            0.5*(lowerMargin + upperMargin);

        point trialPoint(bestPoint);

        if
        (
            projectForMargin
            (
                trialMargin,
                trialPoint
            )
        )
        {
            lowerMargin = trialMargin;
            bestPoint = trialPoint;
        }
        else
        {
            upperMargin = trialMargin;
        }
    }

    apex = bestPoint;

    relativeMargin =
        lowerMargin / localScale;


    // Fail closed on extremely narrow kernels.
    //
    // 1e-6 is dimensionless: the closest supporting face must be at
    // least one-millionth of the local cell scale from the common apex.
    // This is deliberately conservative; unsafe cells remain legal
    // polyhedra instead of being forced into near-zero-height pyramids.
    const scalar minimumRelativeMargin = scalar(1e-6);

    if( relativeMargin < minimumRelativeMargin )
        return false;


    // Final independent proof at the accepted margin. Use a slightly
    // smaller target to avoid classifying pure last-bit binary-search
    // noise as a geometry failure.
    point finalPoint(apex);

    if
    (
        !projectForMargin
        (
            scalar(0.999)*lowerMargin,
            finalPoint
        )
    )
    {
        return false;
    }

    apex = finalPoint;

    return true;
}



bool decomposeCells::exactPyramidChildrenPositive
(
    const label cellI,
    const point& apex,
    scalar& minChildVolume
) const
{
    const cell& c = mesh_.cells()[cellI];
    const faceListPMG& faces = mesh_.faces();
    const labelList& owner = mesh_.owner();
    const pointFieldPMG& points = mesh_.points();

    // Use the exact existing polygon representation which the final
    // children will retain as their base faces.
    const vectorField& faceCentres =
        mesh_.addressingData().faceCentres();

    const vectorField& faceAreas =
        mesh_.addressingData().faceAreas();

    minChildVolume = GREAT;

    forAll(c, cfI)
    {
        const label faceI = c[cfI];
        const face& f = faces[faceI];

        if( f.size() < 3 )
        {
            minChildVolume = -GREAT;
            return false;
        }

        const bool cellIsOwner =
            (owner[faceI] == cellI);

        // ----------------------------------------------------
        // Prospective child corresponding to this parent face:
        //
        //     original polygon base
        //     + one triangular side for every base edge
        //
        // All geometry below is oriented OUTWARD from this
        // prospective child.  This is equivalent to the signed
        // owner/neighbour treatment in checkCellVolumes().
        // ----------------------------------------------------

        const point baseCentre =
            faceCentres[faceI];

        vector baseArea =
            faceAreas[faceI];

        if( !cellIsOwner )
            baseArea = -baseArea;


        // checkCellVolumes() first estimates the cell centre as
        // the arithmetic mean of all face centres.
        point cEst(baseCentre);

        forAll(f, pI)
        {
            const point& pCurrent =
                points[f[pI]];

            const point& pNext =
                points[f.nextLabel(pI)];

            const point sideCentre =
                (1.0/3.0)
               *(pCurrent + pNext + apex);

            cEst += sideCentre;
        }

        cEst /= scalar(f.size() + 1);


        // Base contribution.
        scalar childVol3 =
            baseArea & (baseCentre - cEst);


        // Side-triangle contributions.  Match
        // decomposeCellIntoPyramids() orientation exactly.
        forAll(f, pI)
        {
            const point& pCurrent =
                points[f[pI]];

            const point& pNext =
                points[f.nextLabel(pI)];

            point p0;
            point p1;
            point p2;

            if( cellIsOwner )
            {
                // triFaces:
                //     next, current, apex
                p0 = pNext;
                p1 = pCurrent;
                p2 = apex;
            }
            else
            {
                // reverse triFace used for neighbour-side base:
                //     next, apex, current
                p0 = pNext;
                p1 = apex;
                p2 = pCurrent;
            }

            const point sideCentre =
                (1.0/3.0)*(p0 + p1 + p2);

            const vector sideArea =
                0.5*((p1 - p0) ^ (p2 - p0));

            childVol3 +=
                sideArea & (sideCentre - cEst);
        }


        const scalar childVolume =
            childVol3/3.0;

        minChildVolume =
            Foam::min
            (
                minChildVolume,
                childVolume
            );

        // Exact parity with checkCellVolumes().
        if( childVolume < VSMALL )
            return false;
    }

    return true;
}



label decomposeCells::findTopVertex
(
    const label cellI,
    const DynList<label, 32>& /*vrt*/,
    const DynList<edge, 64>& /*edges*/,
    const DynList<DynList<label, 2>, 64>& /*edgeFaces*/
)
{
    point apex(point::zero);
    scalar relativeMargin = -GREAT;

    // This call occurs only for cells which already passed the exact
    // same non-mutating preflight in decomposeMesh(). Re-run it here
    // fail-closed so no unchecked apex can ever enter the mesh.
    scalar minChildVolume = GREAT;

    if
    (
        !findValidPyramidApex
        (
            cellI,
            apex,
            relativeMargin
        )
     ||
        !exactPyramidChildrenPositive
        (
            cellI,
            apex,
            minChildVolume
        )
    )
    {
        FatalErrorIn
        (
            "label decomposeCells::findTopVertex"
            "(const label, const DynList<label,32>&,"
            " const DynList<edge,64>&,"
            " const DynList<DynList<label,2>,64>&)"
        )   << "Robust pyramid preflight mismatch for cell "
            << cellI
            << ", relativeMargin=" << relativeMargin
            << ", minChildVolume=" << minChildVolume
            << ". Refusing unchecked decomposition."
            << abort(FatalError);
    }

    pointFieldPMG& pointsAccess = mesh_.points();

    const label topVertex = pointsAccess.size();
    pointsAccess.append(apex);

    # ifdef DEBUGDecompose
    Info
        << "Robust top vertex " << topVertex
        << " cell=" << cellI
        << " relativeMargin=" << relativeMargin
        << endl;
    # endif

    return topVertex;
}


void decomposeCells::decomposeCellIntoPyramids(const label cellI)
{
    const cellListPMG& cells = mesh_.cells();
    const faceListPMG& faces = mesh_.faces();
    const labelList& owner = mesh_.owner();

    const cell& c = cells[cellI];

    # ifdef DEBUGDecompose
    Info << "Starting decomposing cell " << cellI << endl;
    Info << "Cell consists of faces " << c << endl;
    forAll(c, fI)
        Info << "Face " << c[fI] << " is " << faces[c[fI]] << endl;
    # endif

    //- calculate edges, faceEdges and edgeFaces addressing
    DynList<label, 32> vrt;
    DynList<edge, 64> edges;
    DynList<DynList<label, 8> > faceEdges;
    faceEdges.setSize(c.size());
    DynList<DynList<label, 2>, 64> edgeFaces;
    findAddressingForCell(cellI, vrt, edges, faceEdges, edgeFaces);

    // find a vertex which will be the top of the pyramids
    //- if there exist a corner vertex which is in 3 or more patches then
    //- it is selected as the top vertex
    label topVertex = findTopVertex(cellI, vrt, edges, edgeFaces);

    //- start generating pyramids
    forAll(c, fI)
    {
        # ifdef DEBUGDecompose
        Info << "Face " << faces[c[fI]] << " is a base face" << endl;
        #endif
        const face& f = faces[c[fI]];
        DynList<DynList<label, 8> > cellFaces;
        cellFaces.setSize(f.size() + 1);

        DynList<triFace> triFaces;
        triFaces.setSize(f.size());
        forAll(triFaces, pI)
        {
            triFaces[pI][0] = f.nextLabel(pI);
            triFaces[pI][1] = f[pI];
            triFaces[pI][2] = topVertex;
        }

        label cfI(0);
        if( owner[c[fI]] == cellI )
        {
            cellFaces[cfI++] = faces[c[fI]];

            forAll(triFaces, tfI)
            {
                cellFaces[cfI++] = triFaces[tfI];
            }
        }
        else
        {
            cellFaces[cfI++] = faces[c[fI]].reverseFace();

            forAll(triFaces, tfI)
            {
                triFace rf;
                rf[0] = triFaces[tfI][0];
                rf[1] = triFaces[tfI][2];
                rf[2] = triFaces[tfI][1];
                cellFaces[cfI++] = rf;
            }
        }

        # ifdef DEBUGDecompose
        Info << "Cell for face is " << cellFaces << endl;

        DynList<edge, 64> cEdges;
        DynList<DynList<label, 2>, 64> eFaces;
        forAll(cellFaces, fI)
        {
            const DynList<label, 8>& f = cellFaces[fI];
            forAll(f, eI)
            {
                const edge e(f[eI], f[(eI+1)%f.size()]);

                const label pos = cEdges.contains(e);

                if( pos < 0 )
                {
                    cEdges.append(e);
                    DynList<label, 2> ef;
                    ef.append(fI);
                    eFaces.append(ef);
                }
                else
                {
                    eFaces[pos].append(fI);
                }
            }
        }

        forAll(eFaces, eI)
            if( eFaces[eI].size() != 2 )
                Pout << "This pyrmid is not topologically closed" << endl;
        # endif

        const label childRecord = facesOfNewCells_.size();

        // Diagnostic only: evaluate EVERY prospective quad-base
        // pyramid using the exact final polygon representation and both
        // explicit diagonal representations.  Print only geometrically
        // suspicious/disagreeing candidates.
        if( cellFaces[0].size() == 4 )
        {
            const pointFieldPMG& pts = mesh_.points();

            auto calcVirtualFaceGeometry =
            [&]
            (
                const face& vf,
                point& fc,
                vector& fa
            )
            {
                const label nPoints = vf.size();

                if( nPoints == 3 )
                {
                    fc =
                        (1.0/3.0)
                       *(
                            pts[vf[0]]
                          + pts[vf[1]]
                          + pts[vf[2]]
                        );

                    fa =
                        0.5
                       *(
                            (pts[vf[1]] - pts[vf[0]])
                          ^ (pts[vf[2]] - pts[vf[0]])
                        );

                    return;
                }

                vector sumN = vector::zero;
                scalar sumA = 0.0;
                vector sumAc = vector::zero;

                point fCentre = pts[vf[0]];

                for
                (
                    label pi=1;
                    pi<nPoints;
                    ++pi
                )
                {
                    fCentre += pts[vf[pi]];
                }

                fCentre /= scalar(nPoints);

                for
                (
                    label pi=0;
                    pi<nPoints;
                    ++pi
                )
                {
                    const label nextI = (pi+1)%nPoints;

                    const point& nextPoint =
                        pts[vf[nextI]];

                    const vector c =
                        pts[vf[pi]]
                      + nextPoint
                      + fCentre;

                    const vector n =
                        (nextPoint - pts[vf[pi]])
                      ^ (fCentre - pts[vf[pi]]);

                    const scalar a = mag(n);

                    sumN += n;
                    sumA += a;
                    sumAc += a*c;
                }

                fc =
                    (1.0/3.0)
                   *sumAc/(sumA + VSMALL);

                fa = 0.5*sumN;
            };


            auto copyCellFace =
            [&]
            (
                const DynList<label, 8>& src
            ) -> face
            {
                face dst(src.size());

                forAll(src, pI)
                    dst[pI] = src[pI];

                return dst;
            };


            // cellFaces are already oriented outward for this prospective
            // child, so no owner/neighbour sign lookup is required here.
            auto virtualCellVolume =
            [&]
            (
                const faceList& virtualFaces
            ) -> scalar
            {
                point cEst(point::zero);

                forAll(virtualFaces, vfI)
                {
                    point fc(point::zero);
                    vector fa(vector::zero);

                    calcVirtualFaceGeometry
                    (
                        virtualFaces[vfI],
                        fc,
                        fa
                    );

                    cEst += fc;
                }

                cEst /= scalar(virtualFaces.size());

                scalar volume = 0.0;

                forAll(virtualFaces, vfI)
                {
                    point fc(point::zero);
                    vector fa(vector::zero);

                    calcVirtualFaceGeometry
                    (
                        virtualFaces[vfI],
                        fc,
                        fa
                    );

                    volume += fa & (fc - cEst);
                }

                return volume/3.0;
            };


            // Actual representation:
            // intact quad base + triangular side faces.
            faceList polyFaces(cellFaces.size());

            forAll(cellFaces, vfI)
            {
                polyFaces[vfI] =
                    copyCellFace(cellFaces[vfI]);
            }

            const scalar polyVol =
                virtualCellVolume(polyFaces);


            const face base =
                copyCellFace(cellFaces[0]);

            faceList fan02Faces
            (
                cellFaces.size() + 1
            );

            faceList fan13Faces
            (
                cellFaces.size() + 1
            );


            // Base diagonal 0--2.
            fan02Faces[0].setSize(3);
            fan02Faces[0][0] = base[0];
            fan02Faces[0][1] = base[1];
            fan02Faces[0][2] = base[2];

            fan02Faces[1].setSize(3);
            fan02Faces[1][0] = base[0];
            fan02Faces[1][1] = base[2];
            fan02Faces[1][2] = base[3];


            // Base diagonal 1--3.
            fan13Faces[0].setSize(3);
            fan13Faces[0][0] = base[1];
            fan13Faces[0][1] = base[2];
            fan13Faces[0][2] = base[3];

            fan13Faces[1].setSize(3);
            fan13Faces[1][0] = base[1];
            fan13Faces[1][1] = base[3];
            fan13Faces[1][2] = base[0];


            for
            (
                label vfI=1;
                vfI<cellFaces.size();
                ++vfI
            )
            {
                fan02Faces[vfI+1] =
                    copyCellFace(cellFaces[vfI]);

                fan13Faces[vfI+1] =
                    copyCellFace(cellFaces[vfI]);
            }


            const scalar fan02Vol =
                virtualCellVolume(fan02Faces);

            const scalar fan13Vol =
                virtualCellVolume(fan13Faces);


            const bool polyBad =
                polyVol < VSMALL;

            const bool fan02Bad =
                fan02Vol < VSMALL;

            const bool fan13Bad =
                fan13Vol < VSMALL;


            // Print only candidates which are already invalid in the
            // exact polygon representation OR where representation choice
            // changes the validity classification.
            if
            (
                polyBad
             || polyBad != fan02Bad
             || polyBad != fan13Bad
            )
            {
                Info
                    << "[DECOMPOSE_QUAD_CHILD_GEOM]"
                    << " record=" << childRecord
                    << " parent=" << cellI
                    << " parentFaceI=" << fI
                    << " baseFace=" << c[fI]
                    << " ownerSide="
                    << (owner[c[fI]] == cellI)
                    << " polyVol=" << polyVol
                    << " fan02Vol=" << fan02Vol
                    << " fan13Vol=" << fan13Vol
                    << " polyBad=" << polyBad
                    << " fan02Bad=" << fan02Bad
                    << " fan13Bad=" << fan13Bad
                    << " apex=" << pts[topVertex]
                    << " basePoints=" << f
                    << endl;
            }
        }

        Info
            << "[DECOMPOSE_CHILD_RECORD]"
            << " record=" << childRecord
            << " parent=" << cellI
            << " parentFaceI=" << fI
            << " baseFace=" << c[fI]
            << " baseNVerts=" << f.size()
            << " ownerSide="
            << (owner[c[fI]] == cellI)
            << " topVertex=" << topVertex
            << " apex=" << mesh_.points()[topVertex]
            << " basePoints=" << f
            << endl;

        facesOfNewCells_.appendGraph(cellFaces);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// ************************************************************************* //
