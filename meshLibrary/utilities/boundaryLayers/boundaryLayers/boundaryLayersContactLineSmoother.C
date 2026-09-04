/*---------------------------------------------------------------------------*\
  Contact-line height limiter for boundary layer vertices.

  Runs after createNewVertices() coordinate swap, before rollback/repair.
  Detects local height spikes and collapses along contact-line points and
  optionally corrects them.

  Fixes vs v1:
  - contactMask prevents averaging across unrelated contact families
  - zeroDistPoints_ anchor prevents lifting truly collapsed singular points
  - Method is non-const (mutates mesh_.points() directly)
  - Atlas write gated by writeContactLineHeightSmootherAtlas_
\*---------------------------------------------------------------------------*/

#include "boundaryLayers.H"
#include "meshSurfaceEngine.H"
#include "OFstream.H"
#include "helperFunctions.H"

#include <algorithm>

namespace Foam
{

void boundaryLayers::smoothContactLineHeights
(
    const word& passName,
    const labelList& bPoints,
    const VRWGraph& pointPoints,
    const boolList& treatPatches,
    const faceList::subList& bFaces,
    const VRWGraph& pointFaces,
    const labelList& boundaryFacePatches
)
{
    // Do NOT call surfaceEngine() here.  Coordinates have already been
    // swapped.  All surface topology used here was captured pre-swap and
    // passed by the caller.
    pointFieldPMG& points = mesh_.points();

    // Contact-family bit mask:
    // bit 1 = BL/no-BL termination
    // bit 2 = BL/neutral seam
    // bit 4 = BL/BL junction
    // bit 8 = BL/BL corner
    labelHashSet candidateSet;
    Map<label> contactMask;

    auto insertMasked = [&](const labelHashSet& src, const label bit)
    {
        forAllConstIter(labelHashSet, src, it)
        {
            const label bpI = it.key();
            candidateSet.insert(bpI);

            if( contactMask.found(bpI) )
                contactMask[bpI] |= bit;
            else
                contactMask.insert(bpI, bit);
        }
    };

    insertMasked(blNoBlEdgePoints_,    1);
    insertMasked(blNeutralEdgePoints_, 2);
    insertMasked(blblJunctionPoints_,  4);
    insertMasked(blblCornerPoints_,    8);

    const label nCandidates = candidateSet.size();

    Info << "Contact-line height smoother (" << passName << "): "
         << nCandidates << " candidate points" << endl;

    if( nCandidates == 0 )
        return;

    // Deterministic processing order.  Neighbor targets remain frozen to
    // the iteration-entry heightMap, while geometric transactions observe
    // already-accepted earlier moves.
    labelList orderedCandidates(nCandidates);
    label orderedI = 0;

    forAllConstIter(labelHashSet, candidateSet, it)
        orderedCandidates[orderedI++] = it.key();

    std::sort(orderedCandidates.begin(), orderedCandidates.end());

    Map<scalar> heightMap;
    Map<vector> dirMap;
    Map<label> basePtMap;

    forAll(orderedCandidates, cI)
    {
        const label bpI = orderedCandidates[cI];

        if( bpI < 0 || bpI >= label(bPoints.size()) )
            continue;

        const label meshPtI = bPoints[bpI];

        if( meshPtI < 0 || meshPtI >= label(points.size()) )
            continue;

        const label basePtI = newLabelForVertex_[meshPtI];

        if( basePtI < 0 || basePtI >= label(points.size()) )
            continue;

        const point base = points[basePtI];
        const point layer = points[meshPtI];
        const vector hVec = layer - base;
        const scalar h = mag(hVec);

        heightMap.insert(bpI, h);
        dirMap.insert
        (
            bpI,
            h > VSMALL ? hVec/h : vector::zero
        );
        basePtMap.insert(bpI, basePtI);
    }

    // ------------------------------------------------------------------
    // Swept triangular-prism Jacobian state.
    //
    // For a triangular wall face B0,B1,B2 and its layer face T0,T1,T2,
    //
    //   X(u,v,t) =
    //       B(u,v) + t [T(u,v)-B(u,v)]
    //
    // is the standard six-node wedge extrusion map.
    //
    // At fixed t the Jacobian is linear over (u,v), therefore its extrema
    // over the reference triangle occur at its three vertices.  At each
    // vertex the Jacobian is quadratic in t.  Evaluating t=0, t=1 and the
    // quadratic stationary point therefore checks the complete swept
    // triangular prism, not just one sampled layer position.
    //
    // Return:
    //   +1 : strictly positive orientation over the full sweep
    //   -1 : strictly negative orientation over the full sweep
    //    0 : zero/crossing/degenerate sweep
    //
    // Orientation itself is not prescribed here; a candidate must preserve
    // the old non-zero orientation.  This makes the guard independent of
    // boundary-face winding convention.
    // ------------------------------------------------------------------
    auto triangleSweepState =
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

        // Scale-aware numerical zero only.  This is deliberately tiny;
        // quality thresholds belong elsewhere.  This gate prevents actual
        // orientation loss, not merely low-quality but valid cells.
        const scalar tol = scalar(1e-12)*L3 + VSMALL;

        scalar globalMin = GREAT;
        scalar globalMax = -GREAT;

        FixedList<vector,3> H;
        H[0] = H0;
        H[1] = H1;
        H[2] = H2;

        for(label k=0; k<3; ++k)
        {
            const vector& Hk = H[k];

            const scalar c0 = (A0 ^ Bv0) & Hk;
            const scalar c1 =
                ((dA ^ Bv0) + (A0 ^ dB)) & Hk;
            const scalar c2 = (dA ^ dB) & Hk;

            auto sampleQuadratic =
            [&](const scalar t)
            {
                const scalar q = c0 + c1*t + c2*t*t;
                globalMin = Foam::min(globalMin, q);
                globalMax = Foam::max(globalMax, q);
            };

            sampleQuadratic(scalar(0));
            sampleQuadratic(scalar(1));

            if( mag(c2) > VSMALL )
            {
                const scalar tStar = -c1/(scalar(2)*c2);

                if( tStar > scalar(0) && tStar < scalar(1) )
                    sampleQuadratic(tStar);
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


    // Evaluate a proposed move against every REGULAR treated BL face incident
    // to bpI.  Transition/reduced faces with collapsed or missing base/top
    // mappings are deliberately left to their specialized topology code.
    //
    // For arbitrary polygonal faces we use a deterministic fan
    // decomposition.  Every fan triangle must preserve its previous sweep
    // orientation.  This is exact for triangular faces and convex polygons
    // represented by the fan; it is conservative for more complicated
    // polygon geometry.
    auto evaluateContactMove =
    [&]
    (
        const label bpI,
        const label meshPtI,
        const point& oldLayer,
        const point& trialLayer,
        bool& oldAllSafe,
        label& badFace,
        label& badPatch,
        scalar& oldWorstMargin,
        scalar& trialWorstMargin,
        label& nRegularFaces
    ) -> bool
    {
        oldAllSafe = true;
        badFace = -1;
        badPatch = -1;
        oldWorstMargin = GREAT;
        trialWorstMargin = GREAT;
        nRegularFaces = 0;

        bool trialAcceptable = true;

        if( bpI < 0 || bpI >= label(pointFaces.size()) )
            return true;

        forAllRow(pointFaces, bpI, pfI)
        {
            const label bfI = pointFaces(bpI, pfI);

            if
            (
                bfI < 0
             || bfI >= label(bFaces.size())
             || bfI >= label(boundaryFacePatches.size())
            )
                continue;

            const label patchI = boundaryFacePatches[bfI];

            if
            (
                patchI < 0
             || patchI >= label(treatPatches.size())
             || !treatPatches[patchI]
            )
                continue;

            if( patchI >= label(patchKey_.size()) )
                continue;

            const face& f = bFaces[bfI];

            if( f.size() < 3 )
                continue;

            const label pKey = patchKey_[patchI];

            List<point> basePts(f.size());
            List<point> oldTopPts(f.size());
            List<point> trialTopPts(f.size());

            bool regularFace = true;

            forAll(f, fpI)
            {
                const label topLabel = f[fpI];

                if( topLabel < 0 || topLabel >= label(points.size()) )
                {
                    regularFace = false;
                    break;
                }

                // Patch-aware lookup falls back to the ordinary resolver and
                // also remains valid for asymmetric-cap bookkeeping.
                const label baseLabel =
                    findNewNodeLabelForPatch(topLabel, patchI, pKey);

                if
                (
                    baseLabel < 0
                 || baseLabel >= label(points.size())
                 || baseLabel == topLabel
                )
                {
                    regularFace = false;
                    break;
                }

                const point baseP = points[baseLabel];
                const point oldTopP =
                    topLabel == meshPtI ? oldLayer : points[topLabel];

                if( mag(oldTopP-baseP) <= scalar(1e-12) )
                {
                    regularFace = false;
                    break;
                }

                basePts[fpI] = baseP;
                oldTopPts[fpI] = oldTopP;
                trialTopPts[fpI] =
                    topLabel == meshPtI ? trialLayer : points[topLabel];
            }

            // Specialized collapsed/reduced topology is not modeled by an
            // ordinary swept prism and must not be falsely frozen by this
            // regular-prism guard.
            if( !regularFace )
                continue;

            ++nRegularFaces;

            for(label ti=1; ti<label(f.size())-1; ++ti)
            {
                scalar oldMargin = scalar(0);
                scalar trialMargin = scalar(0);

                const label oldState =
                    triangleSweepState
                    (
                        basePts[0],
                        basePts[ti],
                        basePts[ti+1],
                        oldTopPts[0],
                        oldTopPts[ti],
                        oldTopPts[ti+1],
                        oldMargin
                    );

                const label trialState =
                    triangleSweepState
                    (
                        basePts[0],
                        basePts[ti],
                        basePts[ti+1],
                        trialTopPts[0],
                        trialTopPts[ti],
                        trialTopPts[ti+1],
                        trialMargin
                    );

                oldWorstMargin =
                    Foam::min(oldWorstMargin, oldMargin);
                trialWorstMargin =
                    Foam::min(trialWorstMargin, trialMargin);

                if( oldState == 0 )
                    oldAllSafe = false;

                bool triAccept = false;

                if( oldState == 0 )
                {
                    // Existing local sweep is already singular/folded:
                    // only accept a move that makes this triangle strictly
                    // admissible.
                    triAccept = (trialState != 0);
                }
                else
                {
                    // Existing sweep is admissible: orientation must be
                    // preserved and may never cross zero.
                    triAccept = (trialState == oldState);
                }

                if( !triAccept )
                {
                    trialAcceptable = false;

                    if( badFace < 0 )
                    {
                        badFace = bfI;
                        badPatch = patchI;
                    }
                }
            }
        }

        // No regular prism faces were available to test.  Preserve legacy
        // behavior and leave such geometry to the specialized transition
        // machinery rather than pretending this validator modeled it.
        if( nRegularFaces == 0 )
        {
            oldAllSafe = true;
            return true;
        }

        return trialAcceptable;
    };


    autoPtr<OFstream> osPtr;

    if( writeContactLineHeightSmootherAtlas_ )
    {
        osPtr.reset
        (
            new OFstream
            (
                "blContactLineHeightSmootherAtlas_"
              + passName
              + ".csv"
            )
        );

        osPtr()
            << "bpI,meshPtI,basePtI,contactMask,cls,"
            << "x,y,z,h,nNei,hMin,hAvg,hMax,ratioToAvg,"
            << "isZeroSeed,wouldSpikeClamp,wouldPinchClamp,hTarget"
            << nl;
    }

    label nSpike = 0;
    label nPinch = 0;
    label nMoved = 0;
    label nAnchored = 0;

    label nTxnTested = 0;
    label nFullAccepted = 0;
    label nBacktracked = 0;
    label nRejected = 0;
    label nPreventedUnsafe = 0;
    label nOldUnsafe = 0;
    label nOldUnsafeLegacyAccepted = 0;
    label nNoRegularFace = 0;
    label nGuardPrint = 0;

    scalar minAcceptedFraction = scalar(1);

    forAll(orderedCandidates, cI)
    {
        const label bpI = orderedCandidates[cI];

        if( !heightMap.found(bpI) )
            continue;

        const label meshPtI = bPoints[bpI];
        const label basePtI = basePtMap[bpI];
        const scalar h = heightMap[bpI];
        const point& base = points[basePtI];

        const label myMask =
            contactMask.found(bpI) ? contactMask[bpI] : 0;

        scalar hSum = scalar(0);
        scalar hMin = GREAT;
        scalar hMax = -GREAT;
        label nNei = 0;

        forAllRow(pointPoints, bpI, ppI)
        {
            const label nbpI = pointPoints(bpI, ppI);

            if( !heightMap.found(nbpI) )
                continue;

            if( !contactMask.found(nbpI) )
                continue;

            if( (myMask & contactMask[nbpI]) == 0 )
                continue;

            const scalar hn = heightMap[nbpI];

            hSum += hn;
            hMin = Foam::min(hMin, hn);
            hMax = Foam::max(hMax, hn);
            ++nNei;
        }

        if( nNei < 2 )
            continue;

        const scalar hAvg = hSum/scalar(nNei);

        const scalar ratioToAvg =
            hAvg > VSMALL ? h/hAvg : scalar(1);

        label cls = 0;

        Map<label>::const_iterator clsIt =
            blContactPointClass_.find(bpI);

        if( clsIt != blContactPointClass_.end() )
            cls = clsIt();

        const bool isRampSeed =
        (
            bpI >= 0
         && bpI < label(rampSeedPoints_.size())
         && rampSeedPoints_[bpI]
        );

        const bool isZeroSeed =
        (
            bpI < label(zeroDistPoints_.size())
         && zeroDistPoints_[bpI]
        )
        && !isRampSeed;

        scalar hTarget = h;
        bool wouldSpike = false;
        bool wouldPinch = false;

        if( hAvg > VSMALL )
        {
            if( h > contactLineSmootherMaxGrowth_*hAvg )
            {
                hTarget =
                    contactLineSmootherMaxGrowth_*hAvg;

                wouldSpike = true;
                ++nSpike;
            }
            else if
            (
                h
              < contactLineSmootherMinCollapse_*hAvg
            )
            {
                hTarget =
                    contactLineSmootherMinCollapse_*hAvg;

                wouldPinch = true;
                ++nPinch;
            }
        }

        hTarget =
            Foam::max
            (
                hTarget,
                contactLineSmootherMinHeight_
            );

        if
        (
            writeContactLineHeightSmootherAtlas_
         && osPtr.valid()
        )
        {
            osPtr()
                << bpI << ","
                << meshPtI << ","
                << basePtI << ","
                << myMask << ","
                << cls << ","
                << base.x() << ","
                << base.y() << ","
                << base.z() << ","
                << h << ","
                << nNei << ","
                << hMin << ","
                << hAvg << ","
                << hMax << ","
                << ratioToAvg << ","
                << (isZeroSeed ? "1" : "0") << ","
                << (wouldSpike ? "1" : "0") << ","
                << (wouldPinch ? "1" : "0") << ","
                << hTarget
                << nl;
        }

        if( !useContactLineHeightSmoother_ )
            continue;

        // Periodic seams are coupled multi-hair structures.
        // Independent scalar height smoothing can preserve the local
        // coarse swept-prism orientation while degrading the later
        // subdivided seam-cell geometry.  Leave these points unchanged.
        if( cls == 1 || cls == 7 || isZeroSeed )
        {
            ++nAnchored;
            continue;
        }

        if( h <= VSMALL )
            continue;

        if( !wouldSpike && !wouldPinch )
            continue;

        const scalar limitedH =
            h
          + contactLineSmootherMaxMoveFraction_
           *(hTarget-h);

        const vector& dir = dirMap[bpI];

        const point oldLayer = points[meshPtI];
        const point fullCandidate =
            base + limitedH*dir;

        if
        (
            help::isnan(fullCandidate)
         || help::isinf(fullCandidate)
        )
        {
            ++nRejected;
            continue;
        }

        ++nTxnTested;

        bool oldAllSafe = true;
        label badFace = -1;
        label badPatch = -1;
        scalar oldWorstMargin = GREAT;
        scalar fullWorstMargin = GREAT;
        label nRegularFaces = 0;

        const bool fullSafe =
            evaluateContactMove
            (
                bpI,
                meshPtI,
                oldLayer,
                fullCandidate,
                oldAllSafe,
                badFace,
                badPatch,
                oldWorstMargin,
                fullWorstMargin,
                nRegularFaces
            );

        if( nRegularFaces == 0 )
            ++nNoRegularFace;

        if( fullSafe )
        {
            points[meshPtI] = fullCandidate;
            ++nFullAccepted;
            ++nMoved;
            continue;
        }

        // The coarse sweep was already non-admissible before this
        // smoother move.  This validator is not authoritative for such
        // geometry, so preserve the legacy smoother result and report it.
        // Only an old-safe -> candidate-unsafe transition may veto/backtrack.
        if( !oldAllSafe )
        {
            ++nOldUnsafe;
            ++nOldUnsafeLegacyAccepted;

            points[meshPtI] = fullCandidate;
            ++nMoved;

            if( nGuardPrint < 20 )
            {
                ++nGuardPrint;

                const word patchName =
                (
                    badPatch >= 0
                 && badPatch < label(patchNames_.size())
                )
                ? patchNames_[badPatch]
                : word("?");

                Info
                    << "CONTACT_SWEEP_DIAG"
                    << " pass=" << passName
                    << " bpI=" << bpI
                    << " meshPtI=" << meshPtI
                    << " mask=" << myMask
                    << " bfI=" << badFace
                    << " patch=" << patchName
                    << " oldSafe=0"
                    << " fullSafe=0"
                    << " action=legacyAccept"
                    << " oldMargin=" << oldWorstMargin
                    << " fullMargin=" << fullWorstMargin
                    << endl;
            }

            continue;
        }

        ++nPreventedUnsafe;

        scalar acceptedFraction = scalar(0);

        // Only a genuinely admissible old state provides a guaranteed safe
        // lower bracket for bisection.  If the old sweep is already singular
        // and the full proposal did not repair it, leave the point unchanged.
        if( oldAllSafe )
        {
            scalar lo = scalar(0);
            scalar hi = scalar(1);

            for(label lsI=0; lsI<12; ++lsI)
            {
                const scalar mid =
                    scalar(0.5)*(lo+hi);

                const point trial =
                    oldLayer
                  + mid*(fullCandidate-oldLayer);

                bool midOldSafe = true;
                label midBadFace = -1;
                label midBadPatch = -1;
                scalar midOldMargin = GREAT;
                scalar midTrialMargin = GREAT;
                label midRegularFaces = 0;

                const bool midSafe =
                    evaluateContactMove
                    (
                        bpI,
                        meshPtI,
                        oldLayer,
                        trial,
                        midOldSafe,
                        midBadFace,
                        midBadPatch,
                        midOldMargin,
                        midTrialMargin,
                        midRegularFaces
                    );

                if( midSafe )
                    lo = mid;
                else
                    hi = mid;
            }

            acceptedFraction = lo;
        }

        if( acceptedFraction > scalar(0) )
        {
            const point accepted =
                oldLayer
              + acceptedFraction
               *(fullCandidate-oldLayer);

            if
            (
                !help::isnan(accepted)
             && !help::isinf(accepted)
            )
            {
                points[meshPtI] = accepted;
                ++nBacktracked;
                ++nMoved;

                minAcceptedFraction =
                    Foam::min
                    (
                        minAcceptedFraction,
                        acceptedFraction
                    );
            }
            else
            {
                acceptedFraction = scalar(0);
                ++nRejected;
            }
        }
        else
        {
            ++nRejected;
        }

        if( nGuardPrint < 20 )
        {
            ++nGuardPrint;

            const word patchName =
            (
                badPatch >= 0
             && badPatch < label(patchNames_.size())
            )
            ? patchNames_[badPatch]
            : word("?");

            Info
                << "CONTACT_SWEEP_GUARD"
                << " pass=" << passName
                << " bpI=" << bpI
                << " meshPtI=" << meshPtI
                << " mask=" << myMask
                << " bfI=" << badFace
                << " patch=" << patchName
                << " oldSafe=" << (oldAllSafe ? 1 : 0)
                << " fullSafe=0"
                << " oldMargin=" << oldWorstMargin
                << " fullMargin=" << fullWorstMargin
                << " acceptedFraction="
                << acceptedFraction
                << endl;
        }
    }

    Info << "Contact-line height smoother (" << passName << "): "
         << "spikes=" << nSpike
         << " pinches=" << nPinch
         << " anchored=" << nAnchored
         << " moved=" << nMoved
         << " (smoother "
         << (useContactLineHeightSmoother_ ? "ON" : "OFF")
         << ")" << endl;

    if( useContactLineHeightSmoother_ )
    {
        Info
            << "CONTACT_SWEEP_TRANSACTION"
            << " pass=" << passName
            << " tested=" << nTxnTested
            << " fullAccepted=" << nFullAccepted
            << " backtracked=" << nBacktracked
            << " rejected=" << nRejected
            << " preventedUnsafe=" << nPreventedUnsafe
            << " oldUnsafe=" << nOldUnsafe
            << " oldUnsafeLegacyAccepted="
            << nOldUnsafeLegacyAccepted
            << " noRegularFace=" << nNoRegularFace
            << " minAcceptedFraction="
            << (nBacktracked > 0
                ? minAcceptedFraction
                : scalar(1))
            << endl;
    }
}

} // End namespace Foam
