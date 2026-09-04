/*---------------------------------------------------------------------------*\
  CFMitch - cfMesh Boundary-Layer Constraint Planner

  V2 height-field atlas.

  Diagnostic only.

  Samples the actual as-born extrusion height immediately after the
  extrusion coordinate swap and before contact-line smoothing or rollback.
\*---------------------------------------------------------------------------*/

#include "boundaryLayers.H"
#include "OFstream.H"

namespace Foam
{

void boundaryLayers::writeCFMitchHeightAtlas
(
    const word& passName,
    const labelList& bPoints,
    const VRWGraph& pointPoints,
    const boolList& treatPatches,
    const labelList& boundaryFacePatches,
    const VRWGraph& pointFaces
) const
{
    if( !cfmitchHeightDiagnostics_ )
        return;

    const pointFieldPMG& points = mesh_.points();

    const label nBP = bPoints.size();

    // ------------------------------------------------------------
    // Phase 1:
    // Resolve current as-born hair height and unique active BL patch
    // for every boundary point.
    // ------------------------------------------------------------

    scalarField height(nBP, scalar(-1));
    labelList primaryPatch(nBP, -1);
    labelList nActivePatches(nBP, 0);

    forAll(bPoints, bpI)
    {
        const label meshPtI = bPoints[bpI];

        if
        (
            meshPtI < 0
         || meshPtI >= label(newLabelForVertex_.size())
        )
            continue;

        const label basePtI =
            newLabelForVertex_[meshPtI];

        if
        (
            basePtI < 0
         || basePtI >= label(points.size())
         || meshPtI >= label(points.size())
        )
            continue;

        height[bpI] =
            mag(points[meshPtI] - points[basePtI]);

        DynList<label> active;

        forAllRow(pointFaces, bpI, pfI)
        {
            const label bfI =
                pointFaces(bpI, pfI);

            if
            (
                bfI < 0
             || bfI >= label(boundaryFacePatches.size())
            )
                continue;

            const label patchI =
                boundaryFacePatches[bfI];

            if
            (
                patchI < 0
             || patchI >= label(treatPatches.size())
             || !treatPatches[patchI]
            )
                continue;

            active.appendIfNotIn(patchI);
        }

        nActivePatches[bpI] = active.size();

        if( active.size() == 1 )
            primaryPatch[bpI] = active[0];
    }


    // ------------------------------------------------------------
    // Resolve smooth-wall eligibility before neighbour statistics.
    //
    // IMPORTANT:
    // A clean interior point must not be judged against a neighbouring
    // point which is itself constrained by a termination, periodic seam,
    // BL/BL junction, corner, ramp or layerScale taper.
    //
    // This mask is deliberately conservative.  CFMitch v2 smoothing will
    // initially be allowed only across edges whose BOTH endpoints satisfy
    // this mask.
    // ------------------------------------------------------------

    boolList baseSmoothCandidate(nBP, false);

    label nBaseSmoothCandidates = 0;
    label nBaseHeightLe1e12 = 0;
    label nBaseHeightLe1e10 = 0;
    label nBaseHeightLe1e8  = 0;
    label nBaseHeightLe1e6  = 0;

    scalar minBaseSmoothHeight = GREAT;
    scalar maxBaseSmoothHeight = scalar(0);

    forAll(bPoints, bpI)
    {
        if( height[bpI] < scalar(0) )
            continue;

        const label patchI =
            primaryPatch[bpI];

        const bool isNoBl =
            blNoBlEdgePoints_.found(bpI);

        const bool isNeutral =
            blNeutralEdgePoints_.found(bpI);

        const bool isBlbl =
            blblJunctionPoints_.found(bpI);

        const bool isCorner =
            blblCornerPoints_.found(bpI);

        const bool isAcute =
            blblAcuteCornerPoints_.found(bpI);

        const bool isRamp =
        (
            bpI >= 0
         && bpI < label(rampSeedPoints_.size())
         && rampSeedPoints_[bpI]
        );

        const label suppressReason =
        (
            bpI >= 0
         && bpI < label(blSuppressReason_.size())
        )
        ? blSuppressReason_[bpI]
        : -1;

        const scalar ls =
        (
            bpI >= 0
         && bpI < label(layerScale_.size())
        )
        ? layerScale_[bpI]
        : scalar(1);

        const bool candidate =
        (
            patchI >= 0
         && nActivePatches[bpI] == 1
         && !isNoBl
         && !isNeutral
         && !isBlbl
         && !isCorner
         && !isAcute
         && !isRamp
         && suppressReason == 0
         && ls >= scalar(1) - SMALL
        );

        if( !candidate )
            continue;

        baseSmoothCandidate[bpI] = true;

        ++nBaseSmoothCandidates;

        minBaseSmoothHeight =
            Foam::min
            (
                minBaseSmoothHeight,
                height[bpI]
            );

        maxBaseSmoothHeight =
            Foam::max
            (
                maxBaseSmoothHeight,
                height[bpI]
            );

        if( height[bpI] <= scalar(1e-12) )
            ++nBaseHeightLe1e12;

        if( height[bpI] <= scalar(1e-10) )
            ++nBaseHeightLe1e10;

        if( height[bpI] <= scalar(1e-8) )
            ++nBaseHeightLe1e8;

        if( height[bpI] <= scalar(1e-6) )
            ++nBaseHeightLe1e6;
    }


    // ------------------------------------------------------------
    // Phase 2:
    // Same-wall neighbour statistics.
    //
    // Only compare points that each belong to exactly one active BL
    // patch and that active patch is the same patch.
    //
    // This prevents wall/periodic, BL/BL and other contacts from
    // contaminating the smooth-wall height-gradient measurement.
    // ------------------------------------------------------------

    fileName atlasName
    (
        "cfmitchHeightAtlas_" + passName + ".csv"
    );

    OFstream os(atlasName);

    os
        << "bpI,meshPoint,patch,nActivePatches,"
        << "baseX,baseY,baseZ,"
        << "topX,topY,topZ,"
        << "height,layerScale,"
        << "samePatchNeighbours,"
        << "neighMinH,neighMeanH,neighMaxH,"
        << "localHeightRatio,maxRelativeJump,"
        << "blNoBl,blNeutral,blblJunction,blblCorner,"
        << "acuteCorner,rampSeed,suppressReason,"
        << "smoothInteriorCandidate,"
        << "birthBaseMeshPoint,birthComponents,"
        << "birthMaxComponentH,birthSumComponentH,"
        << "birthActualOverMax,birthActualOverSum,"
        << "birthMinActualComponentCos,birthMinPairCos,"
        << "birthSumResidualRel"
        << nl;

    label nRows = 0;
    label nSinglePatch = 0;
    label nSmoothCandidates = 0;

    label nRatio125 = 0;
    label nRatio150 = 0;
    label nRatio200 = 0;
    label nRatio400 = 0;

    scalar worstRatio = scalar(1);
    label worstBpI = -1;

    forAll(bPoints, bpI)
    {
        if( height[bpI] < scalar(0) )
            continue;

        const label meshPtI = bPoints[bpI];

        if
        (
            meshPtI < 0
         || meshPtI >= label(newLabelForVertex_.size())
        )
            continue;

        const label basePtI =
            newLabelForVertex_[meshPtI];

        if
        (
            basePtI < 0
         || basePtI >= label(points.size())
         || meshPtI >= label(points.size())
        )
            continue;

        const point& topP =
            points[meshPtI];

        const point& baseP =
            points[basePtI];

        // --------------------------------------------------------
        // CFMitch v2.7.2 diagnostic:
        //
        // At this pre-smoother observation point the coordinate
        // swap has already occurred:
        //
        //   meshPtI = central/as-born extruded point
        //   basePtI = original wall/base coordinate
        //
        // For multi-treatment edge/corner points, diagonal
        // otherVrts_ entries (pKey,pKey) are the independently
        // generated treatment-group vertices.  Their vectors from
        // baseP are therefore the component penetration vectors
        // which were algebraically combined during cleanup.
        //
        // Pair-key entries (pKeyI,pKeyJ) are derived combination
        // vertices and are intentionally NOT counted as independent
        // components.
        // --------------------------------------------------------

        const vector birthActualVector =
            topP - baseP;

        const scalar birthActualH =
            mag(birthActualVector);

        DynList<vector, 4> birthComponentVectors;

        const auto ovIt =
            otherVrts_.find(meshPtI);

        if( ovIt != otherVrts_.end() )
        {
            const auto& ovMap =
                ovIt->second;

            for
            (
                auto it = ovMap.begin();
                it != ovMap.end();
                ++it
            )
            {
                if( it->first.first != it->first.second )
                    continue;

                const label componentPtI =
                    it->second;

                if
                (
                    componentPtI < 0
                 || componentPtI >= label(points.size())
                )
                    continue;

                birthComponentVectors.append
                (
                    points[componentPtI] - baseP
                );
            }
        }

        // Ordinary/single-treatment points have no surviving
        // otherVrts_ diagonal map.  Treat the actual as-born hair
        // itself as their single component.
        if( birthComponentVectors.size() == 0 )
        {
            birthComponentVectors.append
            (
                birthActualVector
            );
        }

        const label birthComponents =
            birthComponentVectors.size();

        scalar birthMaxComponentH =
            scalar(0);

        scalar birthSumComponentH =
            scalar(0);

        vector birthComponentSum =
            vector::zero;

        scalar birthMinActualComponentCos =
            scalar(1);

        forAll(birthComponentVectors, bcI)
        {
            const vector& cv =
                birthComponentVectors[bcI];

            const scalar ch =
                mag(cv);

            birthMaxComponentH =
                Foam::max
                (
                    birthMaxComponentH,
                    ch
                );

            birthSumComponentH +=
                ch;

            birthComponentSum +=
                cv;

            if
            (
                birthActualH > VSMALL
             && ch > VSMALL
            )
            {
                const scalar c =
                    (birthActualVector & cv)
                  / (birthActualH * ch);

                birthMinActualComponentCos =
                    Foam::min
                    (
                        birthMinActualComponentCos,
                        c
                    );
            }
        }

        scalar birthMinPairCos =
            scalar(1);

        if( birthComponents >= 2 )
        {
            for
            (
                label bcI = 0;
                bcI < birthComponents;
                ++bcI
            )
            {
                const scalar hi =
                    mag(birthComponentVectors[bcI]);

                if( hi <= VSMALL )
                    continue;

                for
                (
                    label bcJ = bcI + 1;
                    bcJ < birthComponents;
                    ++bcJ
                )
                {
                    const scalar hj =
                        mag(birthComponentVectors[bcJ]);

                    if( hj <= VSMALL )
                        continue;

                    const scalar c =
                        (
                            birthComponentVectors[bcI]
                          & birthComponentVectors[bcJ]
                        )
                      / (hi * hj);

                    birthMinPairCos =
                        Foam::min
                        (
                            birthMinPairCos,
                            c
                        );
                }
            }
        }

        const scalar birthActualOverMax =
            birthActualH
          / Foam::max
            (
                birthMaxComponentH,
                VSMALL
            );

        const scalar birthActualOverSum =
            birthActualH
          / Foam::max
            (
                birthSumComponentH,
                VSMALL
            );

        // This should be approximately zero wherever our
        // interpretation of the diagonal component vectors is exact.
        // Normalize by sum of component magnitudes so cancellation
        // does not turn the diagnostic into a huge ratio.
        const scalar birthSumResidualRel =
            mag
            (
                birthActualVector
              - birthComponentSum
            )
          / Foam::max
            (
                birthSumComponentH,
                VSMALL
            );

        const label patchI =
            primaryPatch[bpI];

        scalar neighMin = GREAT;
        scalar neighMax = scalar(-1);
        scalar neighSum = scalar(0);
        scalar maxRelativeJump = scalar(0);
        label nNeigh = 0;

        if( patchI >= 0 )
        {
            ++nSinglePatch;

            if( baseSmoothCandidate[bpI] )
            forAllRow(pointPoints, bpI, ppI)
            {
                const label nbpI =
                    pointPoints(bpI, ppI);

                if
                (
                    nbpI < 0
                 || nbpI >= nBP
                 || primaryPatch[nbpI] != patchI
                 || height[nbpI] < scalar(0)
                 || !baseSmoothCandidate[nbpI]
                )
                    continue;

                const scalar hn =
                    height[nbpI];

                neighMin =
                    Foam::min(neighMin, hn);

                neighMax =
                    Foam::max(neighMax, hn);

                neighSum += hn;
                ++nNeigh;

                const scalar denom =
                    Foam::max
                    (
                        Foam::max(height[bpI], hn),
                        VSMALL
                    );

                maxRelativeJump =
                    Foam::max
                    (
                        maxRelativeJump,
                        mag(height[bpI] - hn) / denom
                    );
            }
        }

        scalar neighMean = scalar(-1);
        scalar localRatio = scalar(1);

        if( nNeigh > 0 )
        {
            neighMean =
                neighSum / scalar(nNeigh);

            const scalar localMin =
                Foam::min(height[bpI], neighMin);

            const scalar localMax =
                Foam::max(height[bpI], neighMax);

            // Anything below 1e-12 m is already classified separately
            // as a collapsed/tiny candidate height.  Do not let it turn
            // this diagnostic ratio into an IEEE/VSMALL artefact.
            const scalar ratioFloor = scalar(1e-12);

            localRatio =
                localMax / Foam::max(localMin, ratioFloor);
        }
        else
        {
            neighMin = scalar(-1);
            neighMax = scalar(-1);
        }

        const bool isNoBl =
            blNoBlEdgePoints_.found(bpI);

        const bool isNeutral =
            blNeutralEdgePoints_.found(bpI);

        const bool isBlbl =
            blblJunctionPoints_.found(bpI);

        const bool isCorner =
            blblCornerPoints_.found(bpI);

        const bool isAcute =
            blblAcuteCornerPoints_.found(bpI);

        const bool isRamp =
        (
            bpI >= 0
         && bpI < label(rampSeedPoints_.size())
         && rampSeedPoints_[bpI]
        );

        const label suppressReason =
        (
            bpI >= 0
         && bpI < label(blSuppressReason_.size())
        )
        ? blSuppressReason_[bpI]
        : -1;

        const scalar ls =
        (
            bpI >= 0
         && bpI < label(layerScale_.size())
        )
        ? layerScale_[bpI]
        : scalar(1);

        const bool smoothInterior =
        (
            baseSmoothCandidate[bpI]
         && nNeigh > 0
        );

        if( smoothInterior )
        {
            ++nSmoothCandidates;

            if( localRatio > scalar(1.25) )
                ++nRatio125;

            if( localRatio > scalar(1.50) )
                ++nRatio150;

            if( localRatio > scalar(2.00) )
                ++nRatio200;

            if( localRatio > scalar(4.00) )
                ++nRatio400;

            if( localRatio > worstRatio )
            {
                worstRatio = localRatio;
                worstBpI = bpI;
            }
        }

        word patchName("MULTI_OR_NONE");

        if
        (
            patchI >= 0
         && patchI < label(patchNames_.size())
        )
        {
            patchName =
                patchNames_[patchI];
        }

        os
            << bpI << ','
            << meshPtI << ','
            << patchName << ','
            << nActivePatches[bpI] << ','

            << baseP.x() << ','
            << baseP.y() << ','
            << baseP.z() << ','

            << topP.x() << ','
            << topP.y() << ','
            << topP.z() << ','

            << height[bpI] << ','
            << ls << ','

            << nNeigh << ','
            << neighMin << ','
            << neighMean << ','
            << neighMax << ','
            << localRatio << ','
            << maxRelativeJump << ','

            << (isNoBl ? 1 : 0) << ','
            << (isNeutral ? 1 : 0) << ','
            << (isBlbl ? 1 : 0) << ','
            << (isCorner ? 1 : 0) << ','
            << (isAcute ? 1 : 0) << ','
            << (isRamp ? 1 : 0) << ','
            << suppressReason << ','
            << (smoothInterior ? 1 : 0) << ','

            << basePtI << ','
            << birthComponents << ','
            << birthMaxComponentH << ','
            << birthSumComponentH << ','
            << birthActualOverMax << ','
            << birthActualOverSum << ','
            << birthMinActualComponentCos << ','
            << birthMinPairCos << ','
            << birthSumResidualRel
            << nl;

        ++nRows;
    }

    Info
        << "CFMITCH V2 HEIGHT ATLAS:"
        << " file=" << atlasName
        << " rows=" << nRows
        << " singlePatch=" << nSinglePatch
        << " baseSmoothCandidates=" << nBaseSmoothCandidates
        << " smoothCandidates=" << nSmoothCandidates
        << " baseMinH=" << minBaseSmoothHeight
        << " baseMaxH=" << maxBaseSmoothHeight
        << " baseH<=1e-12=" << nBaseHeightLe1e12
        << " baseH<=1e-10=" << nBaseHeightLe1e10
        << " baseH<=1e-8=" << nBaseHeightLe1e8
        << " baseH<=1e-6=" << nBaseHeightLe1e6
        << " ratio>1.25=" << nRatio125
        << " ratio>1.50=" << nRatio150
        << " ratio>2=" << nRatio200
        << " ratio>4=" << nRatio400
        << " worstRatio=" << worstRatio
        << " worstBpI=" << worstBpI
        << endl;
}


void boundaryLayers::applyCFMitchHeightCompatibility
(
    const word& passName,
    const labelList& bPoints,
    const VRWGraph& pointPoints,
    const boolList& treatPatches,
    const labelList& boundaryFacePatches,
    const VRWGraph& pointFaces
)
{
    if( !cfmitchHeightSmoothing_ )
        return;

    pointFieldPMG& points = mesh_.points();

    const label nBP =
        bPoints.size();

    const scalar maxEdgeRatio =
        Foam::max
        (
            scalar(1.01),
            cfmitchHeightMaxEdgeRatio_
        );

    const scalar logMaxEdgeRatio =
        Foam::log(maxEdgeRatio);

    const scalar maxMoveFraction =
        Foam::max
        (
            scalar(0),
            Foam::min
            (
                scalar(0.20),
                cfmitchHeightMaxMoveFraction_
            )
        );

    const label nIterations =
        Foam::max
        (
            label(0),
            cfmitchHeightSmoothingIterations_
        );


    Info
        << "CFMITCH V2.3 GRAPH HEIGHT SOLVE:"
        << " pass=" << passName
        << " enabled=1"
        << " maxEdgeRatio=" << maxEdgeRatio
        << " maxMoveFraction=" << maxMoveFraction
        << " iterations=" << nIterations
        << endl;


    if
    (
        nIterations == 0
     || maxMoveFraction <= SMALL
    )
        return;


    // ------------------------------------------------------------
    // Active treated-patch classification.
    // ------------------------------------------------------------

    labelList primaryPatch(nBP, -1);
    labelList nActivePatches(nBP, 0);

    forAll(bPoints, bpI)
    {
        DynList<label> active;

        forAllRow(pointFaces, bpI, pfI)
        {
            const label bfI =
                pointFaces(bpI, pfI);

            if
            (
                bfI < 0
             || bfI >= label(boundaryFacePatches.size())
            )
                continue;

            const label patchI =
                boundaryFacePatches[bfI];

            if
            (
                patchI < 0
             || patchI >= label(treatPatches.size())
             || !treatPatches[patchI]
            )
                continue;

            active.appendIfNotIn(patchI);
        }

        nActivePatches[bpI] =
            active.size();

        if( active.size() == 1 )
            primaryPatch[bpI] = active[0];
    }


    // ------------------------------------------------------------
    // Ordinary smooth-wall candidate classification.
    // ------------------------------------------------------------

    boolList smoothCandidate(nBP, false);

    label nSmoothCandidates = 0;

    forAll(bPoints, bpI)
    {
        const label meshPtI =
            bPoints[bpI];

        if
        (
            meshPtI < 0
         || meshPtI >= label(newLabelForVertex_.size())
        )
            continue;

        const label basePtI =
            newLabelForVertex_[meshPtI];

        if
        (
            basePtI < 0
         || basePtI >= label(points.size())
         || meshPtI >= label(points.size())
        )
            continue;

        const bool isNoBl =
            blNoBlEdgePoints_.found(bpI);

        const bool isNeutral =
            blNeutralEdgePoints_.found(bpI);

        const bool isBlbl =
            blblJunctionPoints_.found(bpI);

        const bool isCorner =
            blblCornerPoints_.found(bpI);

        const bool isAcute =
            blblAcuteCornerPoints_.found(bpI);

        const bool isRamp =
        (
            bpI >= 0
         && bpI < label(rampSeedPoints_.size())
         && rampSeedPoints_[bpI]
        );

        const label suppressReason =
        (
            bpI >= 0
         && bpI < label(blSuppressReason_.size())
        )
        ? blSuppressReason_[bpI]
        : -1;

        const scalar ls =
        (
            bpI >= 0
         && bpI < label(layerScale_.size())
        )
        ? layerScale_[bpI]
        : scalar(1);

        if
        (
            primaryPatch[bpI] >= 0
         && nActivePatches[bpI] == 1
         && !isNoBl
         && !isNeutral
         && !isBlbl
         && !isCorner
         && !isAcute
         && !isRamp
         && suppressReason == 0
         && ls >= scalar(1) - SMALL
        )
        {
            smoothCandidate[bpI] = true;
            ++nSmoothCandidates;
        }
    }


    // ------------------------------------------------------------
    // Historical-raw-height verification.
    //
    // Hraw = 0.5 * minimum connected boundary-edge length.
    //
    // A point is allowed to move only when its actual as-born hair
    // still matches this historical sizing heuristic.  Anything
    // shortened by another geometry/crossing constraint remains
    // immutable.
    // ------------------------------------------------------------

    boolList graphEligible(nBP, false);

    label nGraphEligible = 0;
    label nRawMismatch = 0;
    label nNoRawNeighbour = 0;

    scalar worstRelativeRawMismatch =
        scalar(0);


    forAll(bPoints, bpI)
    {
        if( !smoothCandidate[bpI] )
            continue;

        const label meshPtI =
            bPoints[bpI];

        const label basePtI =
            newLabelForVertex_[meshPtI];

        const point& baseP =
            points[basePtI];

        scalar rawHalfEdge = GREAT;
        label nRawNeighbours = 0;


        forAllRow(pointPoints, bpI, ppI)
        {
            const label nbpI =
                pointPoints(bpI, ppI);

            if
            (
                nbpI < 0
             || nbpI >= nBP
            )
                continue;

            const label nbMeshPtI =
                bPoints[nbpI];

            if
            (
                nbMeshPtI < 0
             || nbMeshPtI >= label(newLabelForVertex_.size())
            )
                continue;

            const label nbBasePtI =
                newLabelForVertex_[nbMeshPtI];

            if
            (
                nbBasePtI < 0
             || nbBasePtI >= label(points.size())
            )
                continue;

            const scalar d =
                scalar(0.5)
              * mag
                (
                    points[nbBasePtI]
                  - baseP
                );

            rawHalfEdge =
                Foam::min
                (
                    rawHalfEdge,
                    d
                );

            ++nRawNeighbours;
        }


        if
        (
            nRawNeighbours == 0
         || rawHalfEdge >= GREAT
         || rawHalfEdge <= VSMALL
        )
        {
            ++nNoRawNeighbour;
            continue;
        }


        const scalar bornH =
            mag
            (
                points[meshPtI]
              - points[basePtI]
            );

        const scalar absMismatch =
            mag(bornH - rawHalfEdge);

        const scalar relativeMismatch =
            absMismatch
          / Foam::max
            (
                rawHalfEdge,
                VSMALL
            );

        worstRelativeRawMismatch =
            Foam::max
            (
                worstRelativeRawMismatch,
                relativeMismatch
            );


        const scalar matchTolerance =
            Foam::max
            (
                scalar(1e-12),
                scalar(1e-6)
              * rawHalfEdge
            );


        if( absMismatch <= matchTolerance )
        {
            graphEligible[bpI] = true;
            ++nGraphEligible;
        }
        else
        {
            ++nRawMismatch;
        }
    }


    Info
        << "CFMITCH V2.3 GRAPH HEIGHT ELIGIBILITY:"
        << " pass=" << passName
        << " smoothCandidates=" << nSmoothCandidates
        << " graphEligible=" << nGraphEligible
        << " rawMismatch=" << nRawMismatch
        << " noRawNeighbour=" << nNoRawNeighbour
        << " worstRelativeRawMismatch="
        << worstRelativeRawMismatch
        << endl;


    // ------------------------------------------------------------
    // Symmetric log-height graph projection.
    //
    // For every eligible SAME-PATCH edge:
    //
    //      |log(Hi)-log(Hj)| <= log(maxEdgeRatio)
    //
    // A violating edge contributes equal-and-opposite log-space
    // corrections to its endpoints.
    //
    // Multiple incident edge corrections are averaged before the
    // Jacobi commit.
    //
    // Final point movement is bounded by maxMoveFraction.
    // ------------------------------------------------------------

    for(label iter=0; iter<nIterations; ++iter)
    {
        scalarField oldHeight
        (
            nBP,
            scalar(-1)
        );

        vectorField oldDirection
        (
            nBP,
            vector::zero
        );


        forAll(bPoints, bpI)
        {
            if( !graphEligible[bpI] )
                continue;

            const label meshPtI =
                bPoints[bpI];

            const label basePtI =
                newLabelForVertex_[meshPtI];

            const vector hVec =
                points[meshPtI]
              - points[basePtI];

            const scalar h =
                mag(hVec);

            if( h <= VSMALL )
                continue;

            oldHeight[bpI] =
                h;

            oldDirection[bpI] =
                hVec / h;
        }


        scalarField correctionSum
        (
            nBP,
            scalar(0)
        );

        labelList correctionCount
        (
            nBP,
            label(0)
        );


        label nEdges = 0;
        label nViolatingEdges = 0;

        scalar worstEdgeRatioBefore =
            scalar(1);

        scalar sumViolationExcess =
            scalar(0);


        forAll(bPoints, bpI)
        {
            if
            (
                !graphEligible[bpI]
             || oldHeight[bpI] <= VSMALL
            )
                continue;

            const label patchI =
                primaryPatch[bpI];


            forAllRow(pointPoints, bpI, ppI)
            {
                const label nbpI =
                    pointPoints(bpI, ppI);

                // Unique undirected edge.
                if( nbpI <= bpI )
                    continue;

                if
                (
                    nbpI < 0
                 || nbpI >= nBP
                 || !graphEligible[nbpI]
                 || primaryPatch[nbpI] != patchI
                 || oldHeight[nbpI] <= VSMALL
                )
                    continue;


                ++nEdges;


                const scalar logHi =
                    Foam::log
                    (
                        oldHeight[bpI]
                    );

                const scalar logHj =
                    Foam::log
                    (
                        oldHeight[nbpI]
                    );

                const scalar delta =
                    logHj - logHi;

                const scalar absDelta =
                    mag(delta);

                const scalar edgeRatio =
                    Foam::exp(absDelta);


                worstEdgeRatioBefore =
                    Foam::max
                    (
                        worstEdgeRatioBefore,
                        edgeRatio
                    );


                if
                (
                    absDelta
                 <= logMaxEdgeRatio
                )
                    continue;


                ++nViolatingEdges;


                const scalar excess =
                    absDelta
                  - logMaxEdgeRatio;

                sumViolationExcess +=
                    excess;


                // Exact symmetric pair projection would move each
                // endpoint by excess/2 in opposite directions.
                const scalar halfCorrection =
                    scalar(0.5)
                  * excess;


                if( delta > scalar(0) )
                {
                    // j is taller than i.
                    correctionSum[bpI] +=
                        halfCorrection;

                    correctionSum[nbpI] -=
                        halfCorrection;
                }
                else
                {
                    // i is taller than j.
                    correctionSum[bpI] -=
                        halfCorrection;

                    correctionSum[nbpI] +=
                        halfCorrection;
                }


                ++correctionCount[bpI];
                ++correctionCount[nbpI];
            }
        }


        boolList movePoint
        (
            nBP,
            false
        );

        pointField proposedTop
        (
            nBP,
            point::zero
        );


        label nPointsCorrected = 0;
        label nRaised = 0;
        label nLowered = 0;

        scalar maxRelativeMove =
            scalar(0);

        scalar maxAbsoluteMove =
            scalar(0);

        scalar sumAbsoluteMove =
            scalar(0);

        scalar maxAbsLogCorrection =
            scalar(0);


        forAll(bPoints, bpI)
        {
            if
            (
                !graphEligible[bpI]
             || oldHeight[bpI] <= VSMALL
             || correctionCount[bpI] == 0
            )
                continue;


            scalar logCorrection =
                correctionSum[bpI]
              / scalar(correctionCount[bpI]);


            maxAbsLogCorrection =
                Foam::max
                (
                    maxAbsLogCorrection,
                    mag(logCorrection)
                );


            const scalar oldH =
                oldHeight[bpI];

            scalar requestedH =
                oldH
              * Foam::exp(logCorrection);


            // Bound physical movement around current hair height.
            const scalar minH =
                oldH
              * (
                    scalar(1)
                  - maxMoveFraction
                );

            const scalar maxH =
                oldH
              * (
                    scalar(1)
                  + maxMoveFraction
                );


            const scalar newH =
                Foam::max
                (
                    minH,
                    Foam::min
                    (
                        maxH,
                        requestedH
                    )
                );


            const scalar absMove =
                mag(newH - oldH);


            if
            (
                absMove
             <= scalar(1e-12) * oldH
            )
                continue;


            const label meshPtI =
                bPoints[bpI];

            const label basePtI =
                newLabelForVertex_[meshPtI];


            proposedTop[bpI] =
                points[basePtI]
              + newH
              * oldDirection[bpI];


            movePoint[bpI] = true;

            ++nPointsCorrected;


            if( newH > oldH )
                ++nRaised;
            else
                ++nLowered;


            const scalar relMove =
                absMove / oldH;


            maxRelativeMove =
                Foam::max
                (
                    maxRelativeMove,
                    relMove
                );

            maxAbsoluteMove =
                Foam::max
                (
                    maxAbsoluteMove,
                    absMove
                );

            sumAbsoluteMove +=
                absMove;
        }


        // Jacobi commit.
        forAll(bPoints, bpI)
        {
            if( !movePoint[bpI] )
                continue;

            const label meshPtI =
                bPoints[bpI];

            if
            (
                meshPtI >= 0
             && meshPtI < label(points.size())
            )
            {
                points[meshPtI] =
                    proposedTop[bpI];
            }
        }


        const scalar avgAbsoluteMove =
        (
            nPointsCorrected > 0
          ? sumAbsoluteMove
          / scalar(nPointsCorrected)
          : scalar(0)
        );


        const scalar avgViolationExcess =
        (
            nViolatingEdges > 0
          ? sumViolationExcess
          / scalar(nViolatingEdges)
          : scalar(0)
        );


        Info
            << "CFMITCH V2.3 GRAPH HEIGHT ITER:"
            << " pass=" << passName
            << " iter=" << iter
            << " edges=" << nEdges
            << " violatingEdges=" << nViolatingEdges
            << " correctedPoints=" << nPointsCorrected
            << " raised=" << nRaised
            << " lowered=" << nLowered
            << " worstEdgeRatioBefore="
            << worstEdgeRatioBefore
            << " avgViolationLogExcess="
            << avgViolationExcess
            << " maxAbsLogCorrection="
            << maxAbsLogCorrection
            << " maxRelativeMove="
            << maxRelativeMove
            << " maxAbsoluteMove="
            << maxAbsoluteMove
            << " avgAbsoluteMove="
            << avgAbsoluteMove
            << endl;


        if( nPointsCorrected == 0 )
            break;
    }
}




} // End namespace Foam
