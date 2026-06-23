/*---------------------------------------------------------------------------*\
    BL contact line atlas -- edge-based boundary contact classification.
    Output: blContactLineAtlas.csv (diagnostic only, no policy changes)

    buildBLContactPointPolicy(): builds per-point contact class map,
    consumed by blblSharp suppression block when useContactLinePolicy=true.

    Contact class enum:
    0=Unclassified  1=TripleJunction  2=HardBLBL  3=ModerateBLBL
    4=MildBLBL  5=SmoothBLBL  6=FlowBoundaryTermination  7=PeriodicSeam
    Min-class wins (lower=more restrictive). Ignore class 0 in min logic.
\*---------------------------------------------------------------------------*/

#include "boundaryLayers.H"
#include "meshSurfaceEngine.H"
#include "meshSurfacePartitioner.H"
#include "OFstream.H"

namespace Foam
{

namespace
{
    label contactTypeToClass(const word& contactType)
    {
        if( contactType == "TripleJunction" )          return 1;
        if( contactType == "HardBLBL" )                return 2;
        if( contactType == "ModerateBLBL" )            return 3;
        if( contactType == "MildBLBL" )                return 4;
        if( contactType == "SmoothBLBL" )              return 5;
        if( contactType == "FlowBoundaryTermination" ) return 6;
        if( contactType == "PeriodicSeam" )            return 7;
        return 0;
    }
}

// ---------------------------------------------------------------------------
// Helper: classify a boundary edge given role counts and angle
// ---------------------------------------------------------------------------
static word classifyContactEdge
(
    const label nBL,
    const label nTerm,
    const label nNeutral,
    const scalar angleDeg,
    const bool tripleJunction,
    const word& termPatchName
)
{
    if( tripleJunction )             return "TripleJunction";
    if( nBL == 2 )
    {
        if( angleDeg < 15.0 )        return "SmoothBLBL";
        if( angleDeg < 40.0 )        return "MildBLBL";
        if( angleDeg < 75.0 )        return "ModerateBLBL";
        return "HardBLBL";
    }
    if( nBL == 1 && nTerm == 1 )
    {
        if( termPatchName == "inlet" || termPatchName == "outlet" )
            return "FlowBoundaryTermination";
        return "LayerTermination";
    }
    if( nBL == 1 && nNeutral == 1 )  return "PeriodicSeam";
    return "Unknown";
}

// ---------------------------------------------------------------------------
// Core edge iteration -- shared by both atlas writer and policy builder
// ---------------------------------------------------------------------------
struct ContactEdgeResult
{
    label  bp0, bp1;
    label  pA, pB;
    scalar minDot, angleDeg;
    bool   tripleJunction;
    word   contactType;
};

static bool classifyEdge
(
    const label eI,
    const edgeList& edges,
    const VRWGraph& edgeFaces,
    const labelList& boundaryFacePatches,
    const vectorField& faceNormals,
    const Map<label>& globalToBP,
    const VRWGraph& pPatches,
    const boolList& isBL,
    const boolList& isTerm,
    const boolList& isNeutral,
    const wordList& patchNames,
    const label nPatches,
    ContactEdgeResult& result
)
{
    if( edgeFaces.sizeOfRow(eI) != 2 ) return false;

    const label fA = edgeFaces(eI, 0);
    const label fB = edgeFaces(eI, 1);
    if( fA < 0 || fB < 0 ) return false;
    if( fA >= label(boundaryFacePatches.size()) ) return false;
    if( fB >= label(boundaryFacePatches.size()) ) return false;

    result.pA = boundaryFacePatches[fA];
    result.pB = boundaryFacePatches[fB];
    if( result.pA < 0 || result.pB < 0 ) return false;
    if( result.pA >= nPatches || result.pB >= nPatches ) return false;
    if( result.pA == result.pB ) return false;

    label nBL = 0, nTerm = 0, nNeutral = 0;
    if( isBL[result.pA] ) ++nBL;
    else if( isTerm[result.pA] ) ++nTerm;
    else ++nNeutral;
    if( isBL[result.pB] ) ++nBL;
    else if( isTerm[result.pB] ) ++nTerm;
    else ++nNeutral;
    if( nBL == 0 ) return false;

    vector nA = faceNormals[fA]; nA /= mag(nA) + VSMALL;
    vector nB = faceNormals[fB]; nB /= mag(nB) + VSMALL;
    result.minDot = nA & nB;
    result.angleDeg =
        Foam::acos(Foam::max(scalar(-1), Foam::min(scalar(1), result.minDot)))
        * 180.0 / M_PI;

    const edge& e = edges[eI];
    Map<label>::const_iterator it0 = globalToBP.find(e[0]);
    Map<label>::const_iterator it1 = globalToBP.find(e[1]);
    result.bp0 = (it0 != globalToBP.end()) ? it0() : -1;
    result.bp1 = (it1 != globalToBP.end()) ? it1() : -1;

    result.tripleJunction = false;
    for( label ei = 0; ei < 2; ++ei )
    {
        const label bpI = (ei == 0) ? result.bp0 : result.bp1;
        if( bpI < 0 ) continue;
        label nBLPatches = 0;
        forAllRow(pPatches, bpI, pI)
        {
            const label patchI = pPatches(bpI, pI);
            if( patchI >= 0 && patchI < nPatches && isBL[patchI] )
                ++nBLPatches;
        }
        if( nBLPatches >= 3 ) result.tripleJunction = true;
    }

    const label termPatch =
        isTerm[result.pA] ? result.pA :
        (isTerm[result.pB] ? result.pB : -1);
    const word termPatchName =
        (termPatch >= 0 && termPatch < nPatches) ?
        patchNames[termPatch] : word("unknown");

    result.contactType = classifyContactEdge
    (
        nBL, nTerm, nNeutral,
        result.angleDeg,
        result.tripleJunction,
        termPatchName
    );

    return true;
}

// ---------------------------------------------------------------------------
// Build per-point policy map (consumed by blblSharp suppression)
// ---------------------------------------------------------------------------
void boundaryLayers::buildBLContactPointPolicy() const
{
    blContactPointClass_.clear();

    const meshSurfaceEngine& mse = surfaceEngine();
    const edgeList& edges = mse.edges();
    const VRWGraph& edgeFaces = mse.edgeFaces();
    const labelList& boundaryFacePatches = mse.boundaryFacePatches();
    const vectorField& faceNormals = mse.faceNormals();
    const labelList& bPoints = mse.boundaryPoints();
    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const VRWGraph& pPatches = mPart.pointPatches();

    Map<label> globalToBP;
    forAll(bPoints, bpI)
        globalToBP.insert(bPoints[bpI], bpI);

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

    label nClassified = 0;

    forAll(edges, eI)
    {
        ContactEdgeResult r;
        if( !classifyEdge(eI, edges, edgeFaces, boundaryFacePatches,
                faceNormals, globalToBP, pPatches,
                isBL, isTerm, isNeutral, patchNames_, nPatches, r) )
            continue;

        const label cls = contactTypeToClass(r.contactType);
        if( cls == 0 ) continue;

        for( label ei = 0; ei < 2; ++ei )
        {
            const label bpI = (ei == 0) ? r.bp0 : r.bp1;
            if( bpI < 0 ) continue;
            if( !blContactPointClass_.found(bpI) )
            {
                blContactPointClass_.insert(bpI, cls);
                ++nClassified;
            }
            else if( cls < blContactPointClass_[bpI] )
            {
                blContactPointClass_[bpI] = cls;
            }
        }
    }

    Info << "BL contact-line policy: classified " << nClassified
         << " boundary points from edge atlas" << endl;
}

// ---------------------------------------------------------------------------
// Write diagnostic CSV atlas
// ---------------------------------------------------------------------------
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

    Map<label> globalToBP;
    forAll(bPoints, bpI)
        globalToBP.insert(bPoints[bpI], bpI);

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
        ContactEdgeResult r;
        if( !classifyEdge(eI, edges, edgeFaces, boundaryFacePatches,
                faceNormals, globalToBP, pPatches,
                isBL, isTerm, isNeutral, patchNames_, nPatches, r) )
            continue;

        const edge& e = edges[eI];
        const point mid = 0.5*(points[e[0]] + points[e[1]]);
        const word nameA = (r.pA < nPatches) ? patchNames_[r.pA] : word("unknown");
        const word nameB = (r.pB < nPatches) ? patchNames_[r.pB] : word("unknown");

        label nBL = 0, nTerm = 0, nNeutral = 0;
        if( isBL[r.pA] ) ++nBL; else if( isTerm[r.pA] ) ++nTerm; else ++nNeutral;
        if( isBL[r.pB] ) ++nBL; else if( isTerm[r.pB] ) ++nTerm; else ++nNeutral;

        os << eI << ","
           << r.bp0 << ","
           << r.bp1 << ","
           << mid.x() << ","
           << mid.y() << ","
           << mid.z() << ","
           << nameA << ","
           << nameB << ","
           << nBL << ","
           << nTerm << ","
           << nNeutral << ","
           << r.minDot << ","
           << r.angleDeg << ","
           << r.contactType << ","
           << (r.tripleJunction ? "1" : "0") << nl;

        ++nWritten;
    }

    Info << "BL contact line atlas: wrote " << nWritten
         << " contact edges to blContactLineAtlas.csv" << endl;
}

} // End namespace Foam
