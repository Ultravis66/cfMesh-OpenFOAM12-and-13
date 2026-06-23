/*---------------------------------------------------------------------------*\
    BL contact line atlas -- edge-based boundary contact classification.
    Diagnostic only: no layerScale_ or policy changes.
    Output: blContactLineAtlas.csv

    v1: edge-pair classification + first-pass TripleJunction detection.
    TODO v2: GapClearance, PatchRoleMap roles, policy columns.
\*---------------------------------------------------------------------------*/

#include "boundaryLayers.H"
#include "meshSurfaceEngine.H"
#include "meshSurfacePartitioner.H"
#include "OFstream.H"

namespace Foam
{

void boundaryLayers::writeBLContactLineAtlas() const
{
    const meshSurfaceEngine& mse = surfaceEngine();
    const edgeList& edges = mse.edges();
    const VRWGraph& edgeFaces = mse.edgeFaces();
    const labelList& boundaryFacePatches = mse.boundaryFacePatches();
    const vectorField& faceNormals = mse.faceNormals();
    const labelList& bPoints = mse.boundaryPoints();
    const pointFieldPMG& points = mesh_.points();
    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const VRWGraph& pPatches = mPart.pointPatches();

    // Build global-to-boundary-point map
    Map<label> globalToBP;
    forAll(bPoints, bpI)
        globalToBP.insert(bPoints[bpI], bpI);

    // Patch role helpers from patchRole_
    const label nPatches = patchNames_.size();
    boolList isBL(nPatches, false);
    boolList isTerm(nPatches, false);
    boolList isNeutral(nPatches, false);
    forAll(patchNames_, pI)
    {
        if( pI >= label(patchRole_.size()) ) continue;
        if( patchRole_[pI] == 0 ) isBL[pI] = true;
        else if( patchRole_[pI] == 1 ) isTerm[pI] = true;
        else if( patchRole_[pI] == 2 ) isNeutral[pI] = true;
    }

    OFstream os("blContactLineAtlas.csv");
    os << "edgeI,bp0,bp1,xMid,yMid,zMid,"
       << "patchA,patchB,"
       << "nBL,nTerm,nNeutral,"
       << "minDot,angleDeg,"
       << "contactType,"
       << "tripleJunction" << nl;

    label nWritten = 0;

    forAll(edges, eI)
    {
        if( edgeFaces.sizeOfRow(eI) != 2 ) continue;

        const label fA = edgeFaces(eI, 0);
        const label fB = edgeFaces(eI, 1);

        if( fA < 0 || fB < 0 ) continue;
        if( fA >= label(boundaryFacePatches.size()) ) continue;
        if( fB >= label(boundaryFacePatches.size()) ) continue;

        const label pA = boundaryFacePatches[fA];
        const label pB = boundaryFacePatches[fB];

        if( pA < 0 || pB < 0 ) continue;
        if( pA >= nPatches || pB >= nPatches ) continue;
        if( pA == pB ) continue;

        // Count roles
        label nBL = 0, nTerm = 0, nNeutral = 0;
        if( isBL[pA] ) ++nBL; else if( isTerm[pA] ) ++nTerm; else ++nNeutral;
        if( isBL[pB] ) ++nBL; else if( isTerm[pB] ) ++nTerm; else ++nNeutral;

        if( nBL == 0 ) continue;

        // Dihedral angle
        vector nA = faceNormals[fA];
        vector nB = faceNormals[fB];
        nA /= mag(nA) + VSMALL;
        nB /= mag(nB) + VSMALL;
        const scalar minDot = nA & nB;
        const scalar angleDeg =
            Foam::acos(Foam::max(scalar(-1), Foam::min(scalar(1), minDot)))
            * 180.0 / M_PI;

        // Edge endpoints
        const edge& e = edges[eI];
        Map<label>::const_iterator it0 = globalToBP.find(e[0]);
        Map<label>::const_iterator it1 = globalToBP.find(e[1]);
        const label bp0 = (it0 != globalToBP.end()) ? it0() : -1;
        const label bp1 = (it1 != globalToBP.end()) ? it1() : -1;

        // Triple junction: endpoint touches 3+ BL patches
        bool tripleJunction = false;
        for( label ei = 0; ei < 2; ++ei )
        {
            const label bpI = (ei == 0) ? bp0 : bp1;
            if( bpI < 0 ) continue;
            label nBLPatches = 0;
            forAllRow(pPatches, bpI, pI)
            {
                const label patchI = pPatches(bpI, pI);
                if( patchI >= 0 && patchI < nPatches && isBL[patchI] )
                    ++nBLPatches;
            }
            if( nBLPatches >= 3 ) tripleJunction = true;
        }

        // Classify contact type
        word contactType = "Unknown";
        if( tripleJunction )
        {
            contactType = "TripleJunction";
        }
        else if( nBL == 2 )
        {
            if( angleDeg < 15.0 )       contactType = "SmoothBLBL";
            else if( angleDeg < 40.0 )  contactType = "MildBLBL";
            else if( angleDeg < 75.0 )  contactType = "ModerateBLBL";
            else                        contactType = "HardBLBL";
        }
        else if( nBL == 1 && nTerm == 1 )
        {
            const label termPatch = isTerm[pA] ? pA : pB;
            const word tName =
                (termPatch < nPatches) ?
                patchNames_[termPatch] : word("unknown");
            if( tName == "inlet" || tName == "outlet" )
                contactType = "FlowBoundaryTermination";
            else
                contactType = "LayerTermination";
        }
        else if( nBL == 1 && nNeutral == 1 )
        {
            contactType = "PeriodicSeam";
        }

        // Edge midpoint
        const point mid = 0.5*(points[e[0]] + points[e[1]]);

        const word nameA =
            (pA < nPatches) ? patchNames_[pA] : word("unknown");
        const word nameB =
            (pB < nPatches) ? patchNames_[pB] : word("unknown");

        os << eI << ","
           << bp0 << ","
           << bp1 << ","
           << mid.x() << ","
           << mid.y() << ","
           << mid.z() << ","
           << nameA << ","
           << nameB << ","
           << nBL << ","
           << nTerm << ","
           << nNeutral << ","
           << minDot << ","
           << angleDeg << ","
           << contactType << ","
           << (tripleJunction ? "1" : "0") << nl;

        ++nWritten;
    }

    Info << "BL contact line atlas: wrote " << nWritten
         << " contact edges to blContactLineAtlas.csv" << endl;
}

} // End namespace Foam
