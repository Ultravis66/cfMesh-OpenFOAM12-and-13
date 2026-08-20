/*---------------------------------------------------------------------------*\
    skewRank -- diagnostic utility (project-local, not part of OpenFOAM)

    Reads a polyMesh and computes per-boundary-face skewness using the
    EXACT same meshCheck::boundaryFaceSkewness / meshCheck::faceSkewness
    functions checkMesh itself calls (linked from the stock meshCheck
    library, formula not reimplemented). Ranks all faces by skewness,
    prints the top N with face ID, centroid, owner cell ID, and distance
    to a small hardcoded list of known triple-junction repair positions
    from tonight's session, to directly test whether the worst-skew
    face is spatially associated with any specific repair.

    Read-only. Never modifies the mesh.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "primitiveMeshCheck.H"
#include "SortableList.H"

using namespace Foam;

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createPolyMesh.H"

    const label nTop = 20;

    const pointField& pts = mesh.points();
    const vectorField& fCtrs = mesh.faceCentres();
    const vectorField& fAreas = mesh.faceAreas();
    const vectorField& cellCtrs = mesh.cellCentres();
    const labelList& owner = mesh.faceOwner();

    const label nInternal = mesh.nInternalFaces();
    const label nTotal = mesh.nFaces();

    scalarField skew(nTotal, 0.0);

    // Internal faces: meshCheck::faceSkewness (owner + neighbour).
    const labelList& neighbour = mesh.faceNeighbour();
    for (label facei = 0; facei < nInternal; ++facei)
    {
        skew[facei] = meshCheck::faceSkewness
        (
            mesh, pts, fCtrs, fAreas, facei,
            cellCtrs[owner[facei]], cellCtrs[neighbour[facei]]
        );
    }

    // Boundary faces: meshCheck::boundaryFaceSkewness (owner only).
    for (label facei = nInternal; facei < nTotal; ++facei)
    {
        skew[facei] = meshCheck::boundaryFaceSkewness
        (
            mesh, pts, fCtrs, fAreas, facei,
            cellCtrs[owner[facei]]
        );
    }

    // Sort descending by skewness, keep face IDs.
    labelList faceIds(nTotal);
    forAll(faceIds, i) faceIds[i] = i;

    SortableList<scalar> sortedSkew(skew);
    const labelList& sortOrder = sortedSkew.indices();

    Info << "[SkewRank] total faces=" << nTotal
         << " internal=" << nInternal
         << " boundary=" << (nTotal - nInternal) << endl;

    // Known repair / reference positions from tonight's session.
    struct RefPoint { word name; point pos; };
    List<RefPoint> refPoints(7);
    refPoints[0] = {word("repair_117"),  point(0.175602, -0.0309635, 0.00222714)};
    refPoints[1] = {word("repair_132"),  point(0.175998, -0.0310332, 0.00398567)};
    refPoints[2] = {word("repair_3737"), point(0.241605, -0.0426015, 0.0314922)};
    refPoints[3] = {word("repair_3738"), point(0.241605,  0.0426015, 0.0314922)};
    refPoints[4] = {word("repair_3748"), point(0.241979,  0.0426674, 0.0300462)};
    refPoints[5] = {word("repair_4010"), point(0.248608, -0.0438363, 0.00592228)};
    // Approximate blade_4 sliver region (from earlier session notes,
    // shroud-side, near vertex 4005/4009 neighbourhood).
    refPoints[6] = {word("blade4_region"), point(0.2486, 0.0437, 0.0060)};

    // Print top N worst faces, from the END of the ascending sort.
    for (label rank = 0; rank < nTop && rank < nTotal; ++rank)
    {
        const label idx = sortOrder[nTotal - 1 - rank];
        const scalar sVal = sortedSkew[nTotal - 1 - rank];
        const point& centre = fCtrs[idx];
        const label ownCell = owner[idx];
        const bool isBoundary = (idx >= nInternal);
        const label neiCell = isBoundary ? -1 : neighbour[idx];

        Info << "[SkewRank] rank=" << rank
             << " faceI=" << idx
             << " skew=" << sVal
             << " centre=" << centre
             << " owner=" << ownCell
             << " neighbour=" << neiCell
             << " boundary=" << (isBoundary ? 1 : 0);

        // Distance to each known reference point.
        scalar bestDist = GREAT;
        word bestName("none");
        forAll(refPoints, ri)
        {
            const scalar d = Foam::sqrt(magSqr(centre - refPoints[ri].pos));
            if (d < bestDist)
            {
                bestDist = d;
                bestName = refPoints[ri].name;
            }
        }
        Info << " nearestRef=" << bestName.c_str()
             << " nearestRefDist=" << bestDist
             << endl;
    }

    Info << "End" << endl;
    return 0;
}

// ************************************************************************* //
