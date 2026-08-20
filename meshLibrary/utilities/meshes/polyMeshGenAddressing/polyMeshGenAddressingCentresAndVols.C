/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
     \\/     M anipulation  |
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
    polyMeshGenAddressing

Description
    Efficient cell-centre calculation using face-addressing, face-centres and
    face-areas.

\*---------------------------------------------------------------------------*/

#include "polyMeshGenAddressing.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void polyMeshGenAddressing::calcCellCentresAndVols() const
{
    if( cellCentresPtr_ || cellVolumesPtr_ )
    {
        FatalErrorIn("polyMeshGenAddressing::calcCellCentresAndVols() const")
            << "Cell centres or cell volumes already calculated"
            << abort(FatalError);
    }

    const cellListPMG& cells = mesh_.cells();

    // set the accumulated cell centre to zero vector
    cellCentresPtr_ = new vectorField(cells.size());
    vectorField& cellCtrs = *cellCentresPtr_;

    // Initialise cell volumes to 0
    cellVolumesPtr_ = new scalarField(cells.size());
    scalarField& cellVols = *cellVolumesPtr_;

    // Make centres and volumes
    makeCellCentresAndVols(faceCentres(), faceAreas(), cellCtrs, cellVols);
}

void polyMeshGenAddressing::makeCellCentresAndVols
(
    const vectorField& fCtrs,
    const vectorField& fAreas,
    vectorField& cellCtrs,
    scalarField& cellVols
) const
{
    const labelList& own = mesh_.owner();
    const cellListPMG& cells = mesh_.cells();
    const label nCells = cells.size();

    # ifdef USE_OMP
    # pragma omp parallel for if( nCells > 1000 )
    # endif
    for(label cellI=0;cellI<nCells;++cellI)
    {
        const cell& c = cells[cellI];

        //- approximate the centre first
        vector cEst(vector::zero);
        forAll(c, fI)
            cEst += fCtrs[c[fI]];

        cEst /= c.size();

        //- start evaluating the volume and the cell centre
        vector cellCentre(vector::zero);
        scalar cellVol(0.0);

        forAll(c, fI)
        {
            // Calculate 3*face-pyramid volume
            scalar pyr3Vol = (fAreas[c[fI]] & (fCtrs[c[fI]] - cEst));

            if( own[c[fI]] != cellI )
                pyr3Vol *= -1.0;

            // Addressing deliberately uses the SAME positive surrogate
            // weight for BOTH the centre numerator and its denominator.
            //
            // Splitting them (positive weights in the numerator, signed in
            // the denominator) lets a mixed cell whose pyramids nearly
            // cancel produce numerator ~2e-10 over denominator VSMALL
            // (2.225e-308) -> cell centre ~1e+298. Downstream mag() then
            // squares that and raises SIGFPE.
            //
            // For a mixed/inverted cell there is no useful physical
            // centroid; this bounded surrogate is what stock cfMesh used so
            // optimizer/addressing consumers receive a finite local centre.
            //
            // IMPORTANT: honest signed-volume validation is NOT performed
            // through this cache. polyMeshGenChecks::checkCellVolumes()
            // independently recomputes raw signed pyramid volumes.
            pyr3Vol = Foam::max(pyr3Vol, VSMALL);

            // Calculate face-pyramid centre
            const vector pc = (3.0/4.0)*fCtrs[c[fI]] + (1.0/4.0)*cEst;

            // Numerator and denominator use identical weights.
            cellCentre += pyr3Vol*pc;
            cellVol += pyr3Vol;
        }

        cellCtrs[cellI] = cellCentre / cellVol;
        cellVols[cellI] = cellVol / 3.0;
    }

    // Temporary diagnostic for the 4167a0a centroid-overflow investigation.
    //
    // Do NOT use mag(cellCtrs[cellI]) here: vector mag() squares components
    // and is exactly the arithmetic which SIGFPEs on an astronomical centre.
    // Foam::mag(scalar) is fabs-like and does not square.
    scalar maxCentreCmpt = 0.0;
    label maxCentreCell = -1;
    label nAbsurdCentreComponents = 0;

    forAll(cellCtrs, cellI)
    {
        for(direction cmpt=0; cmpt<vector::nComponents; ++cmpt)
        {
            const scalar value = cellCtrs[cellI][cmpt];

            if( value != value || Foam::mag(value) > VGREAT )
            {
                ++nAbsurdCentreComponents;
                continue;
            }

            const scalar absValue = Foam::mag(value);
            if( absValue > maxCentreCmpt )
            {
                maxCentreCmpt = absValue;
                maxCentreCell = cellI;
            }
        }
    }

    Info << "CELLCTRDIAG maxComponent=" << maxCentreCmpt
         << " cell=" << maxCentreCell
         << " nAbsurdComponents=" << nAbsurdCentreComponents;

    if( maxCentreCell >= 0 )
        Info << " centre=" << cellCtrs[maxCentreCell];

    Info << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const vectorField& polyMeshGenAddressing::cellCentres() const
{
    if( !cellCentresPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const vectorField& polyMeshGenAddressing::cellCentres() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcCellCentresAndVols();
    }

    return *cellCentresPtr_;
}

const scalarField& polyMeshGenAddressing::cellVolumes() const
{
    if( !cellVolumesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const scalarField& polyMeshGenAddressing::cellVolumes() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcCellCentresAndVols();
    }

    return *cellVolumesPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
