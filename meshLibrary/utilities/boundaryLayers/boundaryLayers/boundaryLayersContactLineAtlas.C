/*---------------------------------------------------------------------------*\
    BL contact line atlas -- edge-based boundary contact classification.
    Output: blContactLineAtlas.csv (diagnostic only, no policy changes)

    buildBLContactPointPolicy(): builds per-point contact class map,
    consumed by blblSharp suppression block when useContactLinePolicy=true.

    buildBLCapCellCandidates(): identifies HardBLBL blade/endwall contacts
    for cap cell topology (diagnostic only, gated by useHardBLBLCapCells).

    writeCapCellRoutingAtlas(): proves pKey/otherVrts_ routing before
    any asymmetric vertex assignment.
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
// Build per-point policy map
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

// ---------------------------------------------------------------------------
// Build HardBLBL cap cell candidate map (no geometry change)
// ---------------------------------------------------------------------------
void boundaryLayers::buildBLCapCellCandidates() const
{
    blCapCellEndwallPatch_.clear();
    blCapCellBladePatch_.clear();
    blCapCellEndwallPKey_.clear();
    blCapCellBladePKey_.clear();

    if( !useHardBLBLCapCells_ ) return;

    const meshSurfaceEngine& mse = surfaceEngine();
    const edgeList& edges = mse.edges();
    const VRWGraph& edgeFaces = mse.edgeFaces();
    const labelList& boundaryFacePatches = mse.boundaryFacePatches();
    const vectorField& faceNormals = mse.faceNormals();
    const labelList& bPoints = mse.boundaryPoints();
    Map<label> globalToBP;
    forAll(bPoints, bpI)
        globalToBP.insert(bPoints[bpI], bpI);

    const label nPatches = patchNames_.size();
    boolList isBL(nPatches, false);
    boolList isEndwall(nPatches, false);
    boolList isBlade(nPatches, false);

    forAll(patchNames_, pI)
    {
        if( pI >= label(patchRole_.size()) ) continue;
        if( patchRole_[pI] == 0 )
        {
            isBL[pI] = true;
            const word& pn = patchNames_[pI];
            if( pn == "hub" || pn == "shroud" )
                isEndwall[pI] = true;
            else if( pn.size() >= 5 && pn(0,5) == "blade" )
                isBlade[pI] = true;
        }
    }

    label nCandidates = 0;
    label nDuplicateSame = 0;
    label nConflicts = 0;

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
        if( pA < 0 || pB < 0 || pA >= nPatches || pB >= nPatches ) continue;
        if( pA == pB || !isBL[pA] || !isBL[pB] ) continue;

        // HardBLBL only (angle > 75deg)
        vector nA = faceNormals[fA]; nA /= mag(nA) + VSMALL;
        vector nB = faceNormals[fB]; nB /= mag(nB) + VSMALL;
        const scalar minDot = nA & nB;
        const scalar angleDeg =
            Foam::acos(Foam::max(scalar(-1), Foam::min(scalar(1), minDot)))
            * 180.0 / M_PI;
        if( angleDeg < 75.0 ) continue;

        // Must be blade/endwall pair
        if( !((isEndwall[pA] && isBlade[pB]) || (isEndwall[pB] && isBlade[pA])) )
            continue;

        const label endwallPatch = isEndwall[pA] ? pA : pB;
        const label bladePatch   = isBlade[pA]   ? pA : pB;

        const edge& e = edges[eI];
        for( label ei = 0; ei < 2; ++ei )
        {
            const label gp = (ei == 0) ? e[0] : e[1];
            Map<label>::const_iterator it = globalToBP.find(gp);
            if( it == globalToBP.end() ) continue;
            const label bpI = it();

            if( !blCapCellEndwallPatch_.found(bpI) )
            {
                blCapCellEndwallPatch_.insert(bpI, endwallPatch);
                blCapCellBladePatch_.insert(bpI, bladePatch);
                blCapCellEndwallPKey_.insert(bpI, -1);
                blCapCellBladePKey_.insert(bpI, -1);
                ++nCandidates;
            }
            else
            {
                //- Point already assigned -- check for conflicts
                const bool sameEndwall =
                    (blCapCellEndwallPatch_[bpI] == endwallPatch);
                const bool sameBlade =
                    (blCapCellBladePatch_[bpI] == bladePatch);
                if( sameEndwall && sameBlade )
                    ++nDuplicateSame;
                else
                {
                    ++nConflicts;
                    if( nConflicts <= 3 )
                        Info << "CAP_CONFLICT: bpI=" << bpI
                             << " existing=(" << patchNames_[blCapCellEndwallPatch_[bpI]]
                             << "+" << patchNames_[blCapCellBladePatch_[bpI]] << ")"
                             << " new=(" << patchNames_[endwallPatch]
                             << "+" << patchNames_[bladePatch] << ")" << endl;
                }
            }
        }
    }

    Info << "BL cap cell candidates: " << nCandidates
         << " new, " << nDuplicateSame << " duplicate-same,"
         << " " << nConflicts << " conflicts" << endl;
}

// ---------------------------------------------------------------------------
// Write cap cell pKey routing diagnostic CSV
// Called after patchKey_ is populated inside createNewVertices
// ---------------------------------------------------------------------------
void boundaryLayers::writeCapCellRoutingAtlas() const
{
    if( blCapCellEndwallPatch_.empty() ) return;

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    const pointFieldPMG& points = mesh_.points();

    OFstream os("blCapCellRoutingAtlas.csv");
    os << "bpI,meshPointI,x,y,z,"
       << "endwallPatch,bladePatch,"
       << "endwallPKey,bladePKey,"
       << "otherVrtsSize,"
       << "findForEndwallPKey,findForBladePKey,"
       << "endwallVertLabel,bladeVertLabel,"
       << "distEndwall,distBlade,distFindEndwall,distFindBlade,"
       << "findEWisEWentry,findEWisBLentry,findBLisEWentry,findBLisBLentry"
       << nl;

    forAllConstIter(Map<label>, blCapCellEndwallPatch_, it)
    {
        const label bpI = it.key();
        const label endwallPatch = it();
        const label bladePatch =
            blCapCellBladePatch_.found(bpI) ? blCapCellBladePatch_[bpI] : -1;
        const label endwallPKey =
            blCapCellEndwallPKey_.found(bpI) ? blCapCellEndwallPKey_[bpI] : -1;
        const label bladePKey =
            blCapCellBladePKey_.found(bpI) ? blCapCellBladePKey_[bpI] : -1;

        if( bpI < 0 || bpI >= label(bPoints.size()) ) continue;
        const label meshPtI = bPoints[bpI];
        const point& pt = (meshPtI >= 0 && meshPtI < label(points.size())) ?
            points[meshPtI] : point(Zero);

        label otherVrtsSize = 0;
        label findForEndwall = -1;
        label findForBlade = -1;
        label endwallVertLabel = -1;
        label bladeVertLabel = -1;

        const std::map<label,std::map<std::pair<label,label>,label>>::const_iterator
            oit = otherVrts_.find(meshPtI);
        if( oit != otherVrts_.end() )
        {
            otherVrtsSize = label(oit->second.size());
            if( endwallPKey >= 0 )
                findForEndwall = findNewNodeLabel(meshPtI, endwallPKey);
            if( bladePKey >= 0 )
                findForBlade = findNewNodeLabel(meshPtI, bladePKey);
            if( endwallPKey >= 0 )
            {
                const std::pair<label,label> pr(endwallPKey, endwallPKey);
                auto eit2 = oit->second.find(pr);
                if( eit2 != oit->second.end() )
                    endwallVertLabel = eit2->second;
            }
            if( bladePKey >= 0 )
            {
                const std::pair<label,label> pr(bladePKey, bladePKey);
                auto bit2 = oit->second.find(pr);
                if( bit2 != oit->second.end() )
                    bladeVertLabel = bit2->second;
            }
        }

        const word ewName = (endwallPatch >= 0 && endwallPatch < label(patchNames_.size())) ?
            patchNames_[endwallPatch] : word("?");
        const word blName = (bladePatch >= 0 && bladePatch < label(patchNames_.size())) ?
            patchNames_[bladePatch] : word("?");

        // Compute distances from original point
        const scalar distEW = (endwallVertLabel >= 0 && endwallVertLabel < label(points.size())) ?
            mag(points[endwallVertLabel] - pt) : scalar(-1);
        const scalar distBL = (bladeVertLabel >= 0 && bladeVertLabel < label(points.size())) ?
            mag(points[bladeVertLabel] - pt) : scalar(-1);
        const scalar distFindEW = (findForEndwall >= 0 && findForEndwall < label(points.size())) ?
            mag(points[findForEndwall] - pt) : scalar(-1);
        const scalar distFindBL = (findForBlade >= 0 && findForBlade < label(points.size())) ?
            mag(points[findForBlade] - pt) : scalar(-1);
        // Routing equality: does findNewNodeLabel return EW or BL entry?
        const label findEWisEW = (findForEndwall >= 0 && findForEndwall == endwallVertLabel) ? 1:0;
        const label findEWisBL = (findForEndwall >= 0 && findForEndwall == bladeVertLabel) ? 1:0;
        const label findBLisEW = (findForBlade >= 0 && findForBlade == endwallVertLabel) ? 1:0;
        const label findBLisBL = (findForBlade >= 0 && findForBlade == bladeVertLabel) ? 1:0;
        os << bpI << "," << meshPtI << ","
           << pt.x() << "," << pt.y() << "," << pt.z() << ","
           << ewName << "," << blName << ","
           << endwallPKey << "," << bladePKey << ","
           << otherVrtsSize << ","
           << findForEndwall << "," << findForBlade << ","
           << endwallVertLabel << "," << bladeVertLabel << ","
           << distEW << "," << distBL << ","
           << distFindEW << "," << distFindBL << ","
           << findEWisEW << "," << findEWisBL << ","
           << findBLisEW << "," << findBLisBL << nl;
    }

    Info << "BL cap cell routing atlas written to blCapCellRoutingAtlas.csv"
         << endl;
}

} // End namespace Foam
