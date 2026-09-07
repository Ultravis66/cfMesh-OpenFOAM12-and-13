/*---------------------------------------------------------------------------*\
  CFMitch - cfMesh Boundary-Layer Constraint Planner
\*---------------------------------------------------------------------------*/

#include "boundaryLayerConstraintPlanner.H"
#include "boundaryPatch.H"
#include "IOstreams.H"
#include "error.H"

namespace Foam
{

boundaryLayerConstraintPlanner::architectureType
boundaryLayerConstraintPlanner::parseArchitecture
(
    const word& architectureNameValue
)
{
    if( architectureNameValue == "classicCfMesh" )
        return CLASSIC_CFMESH;

    if( architectureNameValue == "legacyEnhanced" )
        return LEGACY_ENHANCED;

    if( architectureNameValue == "constraintPlanner" )
        return CONSTRAINT_PLANNER;

    FatalErrorIn
    (
        "boundaryLayerConstraintPlanner::parseArchitecture(const word&)"
    )
        << "Unknown boundaryLayerArchitecture '"
        << architectureNameValue << "'" << nl
        << "Valid values are:" << nl
        << "    classicCfMesh" << nl
        << "    legacyEnhanced" << nl
        << "    constraintPlanner"
        << exit(FatalError);

    return LEGACY_ENHANCED;
}


word boundaryLayerConstraintPlanner::architectureName
(
    const architectureType architecture
)
{
    switch( architecture )
    {
        case CLASSIC_CFMESH:
            return word("classicCfMesh");

        case LEGACY_ENHANCED:
            return word("legacyEnhanced");

        case CONSTRAINT_PLANNER:
            return word("constraintPlanner");
    }

    return word("unknown");
}


boundaryLayerConstraintPlanner::boundaryLayerConstraintPlanner
(
    polyMeshGen& mesh,
    const word& architectureNameValue
)
:
    mesh_(mesh),
    architecture_(parseArchitecture(architectureNameValue))
{}


void boundaryLayerConstraintPlanner::report() const
{
    Info
        << "CFMITCH: Boundary Layer Constraint Planner v1"
        << " architecture="
        << architectureName(architecture_)
        << " meshPoints=" << mesh_.points().size()
        << endl;

    if( architecture_ == CLASSIC_CFMESH )
    {
        FatalErrorIn
        (
            "boundaryLayerConstraintPlanner::report()"
        )
            << "boundaryLayerArchitecture classicCfMesh was selected, "
            << "but the verified classic cfMesh execution path has not "
            << "yet been wired." << nl
            << "Use legacyEnhanced or constraintPlanner."
            << exit(FatalError);
    }

    if( architecture_ == LEGACY_ENHANCED )
    {
        Info
            << "CFMITCH: legacyEnhanced path selected; "
            << "constraint planner inactive"
            << endl;
    }
    else
    {
        Info
            << "CFMITCH: constraintPlanner selected; "
            << "resolved boundary-layer plan active"
            << endl;
    }
}


boundaryLayerPlan
boundaryLayerConstraintPlanner::solveLayerCounts
(
    const meshSurfaceEngine& mse,
    const LongList<edge>& splitEdges,
    const Map<scalar>& stableHeightScaleAtMeshPoint,
    const labelList& requestedFaceLayers,
    const scalar globalThicknessRatio,
    const std::map<word, scalar>& thicknessRatioForPatch,
    const labelList& vtFaceRing,
    const labelHashSet& directLayerCapFaces,
    const label maxLayerStep
) const
{
    if( maxLayerStep < 1 )
    {
        FatalErrorIn
        (
            "boundaryLayerConstraintPlanner::solveLayerCounts(...)"
        )
            << "maxLayerStep must be >= 1, received "
            << maxLayerStep
            << exit(FatalError);
    }

    boundaryLayerPlan plan(requestedFaceLayers);
    labelList& faceLayers = plan.faceLayers();

    const PtrList<boundaryPatch>& boundaries =
        mesh_.boundaries();

    const faceList::subList& bFaces =
        mse.boundaryFaces();

    const VRWGraph& pointFaces =
        mse.pointFaces();

    const labelList& facePatch =
        mse.boundaryFacePatches();

    // mesh-point -> current boundary-point addressing
    const labelList& bp =
        mse.bp();


    // -----------------------------------------------------------------
    // Phase 1: stable point constraints
    // -----------------------------------------------------------------

    scalarField heightScaleByMeshPoint
    (
        mesh_.points().size(),
        scalar(1)
    );

    label nStableMapped = 0;

    forAllConstIter
    (
        Map<scalar>,
        stableHeightScaleAtMeshPoint,
        it
    )
    {
        const label meshPtI = it.key();

        if
        (
            meshPtI < 0
         || meshPtI >= label(heightScaleByMeshPoint.size())
        )
            continue;

        heightScaleByMeshPoint[meshPtI] =
            Foam::max
            (
                scalar(0),
                Foam::min(scalar(1), it())
            );

        ++nStableMapped;
    }


    // -----------------------------------------------------------------
    // Phase 2: map stable point constraints through current hair edges
    //          onto current boundary points.
    // -----------------------------------------------------------------

    scalarField heightScaleByCurrentBp
    (
        mse.boundaryPoints().size(),
        scalar(1)
    );

    boolList heightScalePresentAtCurrentBp
    (
        mse.boundaryPoints().size(),
        false
    );

    label nSourceEdges = 0;
    label nBpUpdates = 0;

    forAll(splitEdges, seI)
    {
        const edge& se = splitEdges[seI];

        const label ep0 = se.start();
        const label ep1 = se.end();

        scalar edgeScale = scalar(1);
        bool hasEdgeScale = false;

        if
        (
            ep0 >= 0
         && ep0 < label(heightScaleByMeshPoint.size())
         && heightScaleByMeshPoint[ep0] < scalar(1) - SMALL
        )
        {
            edgeScale = heightScaleByMeshPoint[ep0];
            hasEdgeScale = true;
        }

        if
        (
            ep1 >= 0
         && ep1 < label(heightScaleByMeshPoint.size())
         && heightScaleByMeshPoint[ep1] < scalar(1) - SMALL
        )
        {
            edgeScale =
                hasEdgeScale
              ? Foam::min
                (
                    edgeScale,
                    heightScaleByMeshPoint[ep1]
                )
              : heightScaleByMeshPoint[ep1];

            hasEdgeScale = true;
        }

        if( !hasEdgeScale )
            continue;

        ++nSourceEdges;

        const label endpoints[2] = {ep0, ep1};

        for( label ei=0; ei<2; ++ei )
        {
            const label meshPtI = endpoints[ei];

            if
            (
                meshPtI < 0
             || meshPtI >= label(bp.size())
            )
                continue;

            const label currBpI = bp[meshPtI];

            if
            (
                currBpI < 0
             || currBpI >= label(heightScaleByCurrentBp.size())
            )
                continue;

            if( heightScalePresentAtCurrentBp[currBpI] )
            {
                heightScaleByCurrentBp[currBpI] =
                    Foam::min
                    (
                        heightScaleByCurrentBp[currBpI],
                        edgeScale
                    );
            }
            else
            {
                heightScaleByCurrentBp[currBpI] =
                    edgeScale;

                heightScalePresentAtCurrentBp[currBpI] =
                    true;
            }

            ++nBpUpdates;
        }
    }


    // -----------------------------------------------------------------
    // Phase 3: convert available height into a hard face-layer upper
    //          bound while preserving near-wall geometric progression.
    // -----------------------------------------------------------------

    boolList directCappedFace
    (
        faceLayers.size(),
        false
    );

    label nDirectCappedFaces = 0;
    label nExplicitLayerCapSeeds = 0;
    label minDirectLayers = labelMax;
    label maxDirectLayers = 0;

    forAllConstIter
    (
        labelHashSet,
        directLayerCapFaces,
        seedIt
    )
    {
        const label bfI = seedIt.key();

        if( bfI < 0 || bfI >= label(faceLayers.size()) )
            continue;

        directCappedFace[bfI] = true;
        ++nDirectCappedFaces;
        ++nExplicitLayerCapSeeds;

        minDirectLayers =
            Foam::min(minDirectLayers, faceLayers[bfI]);

        maxDirectLayers =
            Foam::max(maxDirectLayers, faceLayers[bfI]);
    }

    forAll(faceLayers, bfI)
    {
        const label requestedNLayers =
            faceLayers[bfI];

        if( requestedNLayers <= 1 )
            continue;

        if( bfI < 0 || bfI >= label(bFaces.size()) )
            continue;

        const face& bf = bFaces[bfI];

        scalar faceScale = scalar(1);
        bool hasFaceScale = false;

        forAll(bf, fpI)
        {
            const label meshPtI = bf[fpI];

            if
            (
                meshPtI < 0
             || meshPtI >= label(bp.size())
            )
                continue;

            const label currBpI = bp[meshPtI];

            if
            (
                currBpI < 0
             || currBpI >= label(heightScaleByCurrentBp.size())
             || !heightScalePresentAtCurrentBp[currBpI]
            )
                continue;

            const scalar sPt =
                heightScaleByCurrentBp[currBpI];

            faceScale =
                hasFaceScale
              ? Foam::min(faceScale, sPt)
              : sPt;

            hasFaceScale = true;
        }

        if( !hasFaceScale )
            continue;

        scalar ratio = globalThicknessRatio;

        if( bfI < label(facePatch.size()) )
        {
            const label patchI = facePatch[bfI];

            if
            (
                patchI >= 0
             && patchI < label(boundaries.size())
            )
            {
                const word& patchName =
                    boundaries[patchI].patchName();

                const std::map<word, scalar>::const_iterator rIt =
                    thicknessRatioForPatch.find(patchName);

                if( rIt != thicknessRatioForPatch.end() )
                    ratio = rIt->second;
            }
        }

        const scalar r =
            Foam::max(ratio, scalar(1e-12));

        scalar totalWeight = scalar(0);
        scalar w = scalar(1);

        for
        (
            label li=0;
            li<requestedNLayers;
            ++li
        )
        {
            totalWeight += w;
            w *= r;
        }

        const scalar availableWeight =
            faceScale * totalWeight;

        label faceCap = 1;
        scalar cumulative = scalar(0);
        w = scalar(1);

        const scalar fitTol =
            scalar(100) * SMALL
          * Foam::max(scalar(1), totalWeight);

        for
        (
            label li=1;
            li<=requestedNLayers;
            ++li
        )
        {
            cumulative += w;

            if
            (
                cumulative
             <= availableWeight + fitTol
            )
            {
                faceCap = li;
            }
            else
            {
                break;
            }

            w *= r;
        }

        faceCap =
            Foam::max
            (
                label(1),
                Foam::min(requestedNLayers, faceCap)
            );

        if( faceCap < faceLayers[bfI] )
        {
            faceLayers[bfI] = faceCap;
            if( !directCappedFace[bfI] )
            {
                directCappedFace[bfI] = true;
                ++nDirectCappedFaces;
            }

            minDirectLayers =
                Foam::min
                (
                    minDirectLayers,
                    faceCap
                );

            maxDirectLayers =
                Foam::max
                (
                    maxDirectLayers,
                    faceCap
                );
        }
    }


    // -----------------------------------------------------------------
    // Phase 4: solve same-patch layer-count compatibility.
    //
    // Constraints:
    //
    // final(face) <= local hard upper bound
    //
    // final(neighbour) <= final(face) + maxLayerStep
    //
    // The relaxation only propagates from direct constrained seeds.
    // Existing unrelated low-layer faces therefore cannot seed the solve.
    // -----------------------------------------------------------------

    boolList frontier(directCappedFace);

    boolList compatibilityAdjustedFace
    (
        faceLayers.size(),
        false
    );

    label nCompatibilityAdjustedFaces = 0;
    label nCompatibilityUpdates = 0;
    label nCompatibilityPasses = 0;

    label maxRequestedLayers = 1;

    forAll(faceLayers, bfI)
    {
        maxRequestedLayers =
            Foam::max
            (
                maxRequestedLayers,
                faceLayers[bfI]
            );
    }

    const label maxPasses =
        maxRequestedLayers + 2;

    bool changed = true;

    while
    (
        changed
     && nCompatibilityPasses < maxPasses
    )
    {
        changed = false;
        ++nCompatibilityPasses;

        boolList nextFrontier
        (
            faceLayers.size(),
            false
        );

        // One topological ring per pass, independent of face ordering.
        const labelList layersBeforePass(faceLayers);

        forAll(frontier, bfI)
        {
            if( !frontier[bfI] )
                continue;

            if
            (
                bfI < 0
             || bfI >= label(bFaces.size())
             || bfI >= label(facePatch.size())
            )
                continue;

            const label sourcePatch =
                facePatch[bfI];

            const label neighbourCap =
                layersBeforePass[bfI] + maxLayerStep;

            const face& f = bFaces[bfI];

            forAll(f, fpI)
            {
                const label meshPtI = f[fpI];

                if
                (
                    meshPtI < 0
                 || meshPtI >= label(bp.size())
                )
                    continue;

                const label bpI = bp[meshPtI];

                if
                (
                    bpI < 0
                 || bpI >= label(pointFaces.size())
                )
                    continue;

                forAllRow(pointFaces, bpI, pfI)
                {
                    const label nbfI =
                        pointFaces(bpI, pfI);

                    if
                    (
                        nbfI < 0
                     || nbfI >= label(faceLayers.size())
                     || nbfI >= label(facePatch.size())
                     || nbfI == bfI
                    )
                        continue;

                    // Compatibility is currently same-patch only.
                    if( facePatch[nbfI] != sourcePatch )
                        continue;

                    // Existing virtual-topology decisions have priority.
                    if
                    (
                        nbfI < label(vtFaceRing.size())
                     && vtFaceRing[nbfI] >= 0
                    )
                        continue;

                    if
                    (
                        layersBeforePass[nbfI]
                     <= neighbourCap
                    )
                        continue;

                    if( faceLayers[nbfI] > neighbourCap )
                    {
                        faceLayers[nbfI] =
                            neighbourCap;

                        nextFrontier[nbfI] = true;
                        changed = true;
                        ++nCompatibilityUpdates;

                        if
                        (
                            !compatibilityAdjustedFace[nbfI]
                        )
                        {
                            compatibilityAdjustedFace[nbfI] =
                                true;

                            ++nCompatibilityAdjustedFaces;
                        }
                    }
                }
            }
        }

        frontier.transfer(nextFrontier);
    }


    Info
        << "CFMITCH PLAN layerCounts:"
        << " maxStep=" << maxLayerStep
        << " stableInput="
        << stableHeightScaleAtMeshPoint.size()
        << " stableMapped=" << nStableMapped
        << " sourceEdges=" << nSourceEdges
        << " bpUpdates=" << nBpUpdates
        << " explicitLayerCapSeeds="
        << nExplicitLayerCapSeeds
        << " directCappedFaces=" << nDirectCappedFaces;

    if( nDirectCappedFaces > 0 )
    {
        Info
            << " minDirectLayers=" << minDirectLayers
            << " maxDirectLayers=" << maxDirectLayers;
    }

    Info
        << " compatibilityAdjustedFaces="
        << nCompatibilityAdjustedFaces
        << " compatibilityUpdates="
        << nCompatibilityUpdates
        << " passes=" << nCompatibilityPasses
        << " maxPasses=" << maxPasses
        << endl;


    return plan;
}

} // End namespace Foam
