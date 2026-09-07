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
#include "meshSurfaceEngine.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const meshSurfaceEngine& refineBoundaryLayers::surfaceEngine() const
{
    if( !msePtr_ )
        msePtr_ = new meshSurfaceEngine(mesh_);

    return *msePtr_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

refineBoundaryLayers::refineBoundaryLayers(polyMeshGen& mesh)
:
    mesh_(mesh),
    msePtr_(NULL),
    globalNumLayers_(1),
    globalThicknessRatio_(1.0),
    globalMaxThicknessFirstLayer_(VGREAT),
    numLayersForPatch_(),
    thicknessRatioForPatch_(),
    maxThicknessForPatch_(),
    discontinuousLayersForPatch_(),
    cellSubsetName_(),
    done_(false),
    is2DMesh_(false),
    specialMode_(false),
    refinementValid_(true),
    boundaryLayerArchitecture_("legacyEnhanced"),
    constraintPlannerMaxLayerStep_(2),
    cfmitchV52aFrontCensus_(false),
    cfmitchV52bApplyFrontCaps_(false),
    nLayersAtBndFace_(),
    cellToBaseBndFace_(),
    qualityMaxLayersAtFace_(),
    splitEdges_(),
    splitEdgesAtPoint_(),
    newVerticesForSplitEdge_(),
    facesFromFace_(),
    newFaces_(),
    acuteCornerCapLayers_(false),
    gapActionPoints_(),
    gapLoserPatchNames_(),
    gapRing1MaxLayers_(1),
    gapRing2MaxLayers_(2),
    blTerminationEdgePoints_(),
    blTerminationRing1MaxLayers_(3),
    blTerminationRing2MaxLayers_(3)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

refineBoundaryLayers::~refineBoundaryLayers()
{
    deleteDemandDrivenData(msePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void refineBoundaryLayers::avoidRefinement()
{
    if( done_ )
    {
        FatalErrorIn
        (
            "void refineBoundaryLayers::avoidRefinement()"
        ) << "refineLayers is already executed" << exit(FatalError);
    }

    globalNumLayers_ = 1;
    numLayersForPatch_.clear();
}

void refineBoundaryLayers::activate2DMode()
{
    if( done_ )
    {
        FatalErrorIn
        (
            "void refineBoundaryLayers::activate2DMode()"
        ) << "refineLayers is already executed" << exit(FatalError);
    }

    is2DMesh_ = true;
}

void refineBoundaryLayers::setGlobalNumberOfLayers(const label nLayers)
{
    if( done_ )
    {
        FatalErrorIn
        (
            "void refineBoundaryLayers::setGlobalNumberOfLayers(const label)"
        ) << "refineLayers is already executed" << exit(FatalError);
    }

    if( nLayers < 2 )
    {
        WarningIn
        (
            "void refineBoundaryLayers::setGlobalNumberOfLayers(const label)"
        ) << "The specified global number of boundary layers is less than 2"
          << endl;

        return;
    }

    globalNumLayers_ = nLayers;
}

void refineBoundaryLayers::setGlobalThicknessRatio(const scalar thicknessRatio)
{
    if( done_ )
    {
        FatalErrorIn
        (
            "void refineBoundaryLayers::setGlobalThicknessRatio(const scalar)"
        ) << "refineLayers is already executed" << exit(FatalError);
    }

    if( thicknessRatio < 1.0 )
    {
        WarningIn
        (
            "void refineBoundaryLayers::setGlobalThicknessRatio(const scalar)"
        ) << "The specified global thickness ratio is less than 1.0" << endl;

        return;
    }

    globalThicknessRatio_ = thicknessRatio;
}

void refineBoundaryLayers::setGlobalMaxThicknessOfFirstLayer
(
    const scalar maxThickness
)
{
    if( done_ )
    {
        FatalErrorIn
        (
            "void refineBoundaryLayers::setGlobalMaxThicknessOfFirstLayer"
            "(const scalar)"
        ) << "refineLayers is already executed" << exit(FatalError);
    }

    if( maxThickness <= 0.0 )
    {
        WarningIn
        (
            "void refineBoundaryLayers::setGlobalMaxThicknessOfFirstLayer"
            "(const scalar)"
        ) << "The specified global maximum thickness of the first"
          << " boundary layer is negative!!" << endl;

        return;
    }

    globalMaxThicknessFirstLayer_ = maxThickness;
}

void refineBoundaryLayers::setNumberOfLayersForPatch
(
    const word& patchName,
    const label nLayers
)
{
    if( done_ )
    {
        FatalErrorIn
        (
            "void refineBoundaryLayers::setNumberOfLayersForPatch"
            "(const word&, const label)"
        ) << "refineLayers is already executed" << exit(FatalError);
    }

    if( nLayers < 2 )
    {
        WarningIn
        (
            "void refineBoundaryLayers::setNumberOfLayersForPatch"
            "(const word&, const label)"
        ) << "The specified number of boundary layers for patch " << patchName
          << " is less than 2... boundary layers disabled for this patch!" << endl;
    }

    const labelList matchedIDs = mesh_.findPatches(patchName);

    forAll(matchedIDs, matchI)
    {
        numLayersForPatch_[mesh_.getPatchName(matchedIDs[matchI])] = nLayers;
    }
}

void refineBoundaryLayers::setThicknessRatioForPatch
(
    const word& patchName,
    const scalar thicknessRatio
)
{
    if( done_ )
    {
        FatalErrorIn
        (
            "void refineBoundaryLayers::setThicknessRatioForPatch"
            "(const word&, const scalar)"
        ) << "refineLayers is already executed" << exit(FatalError);
    }

    if( thicknessRatio < 1.0 )
    {
        WarningIn
        (
            "void refineBoundaryLayers::setThicknessRatioForPatch"
            "(const word&, const scalar)"
        ) << "The specified thickness ratio for patch " << patchName
          << " is less than 1.0" << endl;

        return;
    }

    const labelList matchedIDs = mesh_.findPatches(patchName);

    forAll(matchedIDs, matchI)
    {
        const word pName = mesh_.getPatchName(matchedIDs[matchI]);
        thicknessRatioForPatch_[pName] = thicknessRatio;
    }
}

void refineBoundaryLayers::setMaxThicknessOfFirstLayerForPatch
(
    const word& patchName,
    const scalar maxThickness
)
{
    if( done_ )
    {
        FatalErrorIn
        (
            "void refineBoundaryLayers::setMaxThicknessOfFirstLayerForPatch"
            "(const word&, const scalar)"
        ) << "refineLayers is already executed" << exit(FatalError);
    }

    if( maxThickness <= 0.0 )
    {
        WarningIn
        (
            "void refineBoundaryLayers::setGlobalMaxThicknessOfFirstLayer"
            "(const word&, const scalar)"
        ) << "The specified maximum thickness of the first boundary layer "
          << "for patch " << patchName << " is negative!!" << endl;

        return;
    }

    const labelList matchedIDs = mesh_.findPatches(patchName);

    forAll(matchedIDs, matchI)
    {
        const word pName = mesh_.getPatchName(matchedIDs[matchI]);
        maxThicknessForPatch_[pName] = maxThickness;
    }
}

void refineBoundaryLayers::setInteruptForPatch(const word& patchName)
{
    if( done_ )
    {
        FatalErrorIn
        (
            "void refineBoundaryLayers::setInteruptForPatch(const word&)"
        ) << "refineLayers is already executed" << exit(FatalError);
    }

    const labelList matchedIDs = mesh_.findPatches(patchName);

    forAll(matchedIDs, matchI)
    {
        const word pName = mesh_.getPatchName(matchedIDs[matchI]);
        discontinuousLayersForPatch_.insert(pName);
    }
}

void refineBoundaryLayers::setCellSubset(const word subsetName)
{
    if( done_ )
    {
        FatalErrorIn
        (
            "void refineBoundaryLayers::setCellSubset(const word)"
        ) << "refineLayers is already executed" << exit(FatalError);
    }

    cellSubsetName_ = subsetName;
}

void refineBoundaryLayers::activateSpecialMode()
{
    specialMode_ = true;
}

void refineBoundaryLayers::refineLayers()
{
    bool refinePatch(false);
    for
    (
        std::map<word, label>::const_iterator it=numLayersForPatch_.begin();
        it!=numLayersForPatch_.end();
        ++it
    )
        if( it->second > 1 )
            refinePatch = true;

    if( (globalNumLayers_ < 2) && !refinePatch )
        return;

    Info << "Starting refining boundary layers" << endl;

    if( done_ )
    {
        WarningIn
        (
            "void refineBoundaryLayers::refineLayers()"
        ) << "Boundary layers are already refined! "
          << "Stopping refinement" << endl;

        return;
    }

    if( !analyseLayers() )
    {
        WarningIn
        (
            "void refineBoundaryLayers::refineLayers()"
        ) << "Boundary layers do not exist in the mesh! Cannot refine" << endl;

        return;
    }

    generateNewVertices();

    generateNewFaces();

    //- Fail closed: if face refinement hit inconsistent split-edge
    //- metadata, abort before generateNewCells() mutates more topology.
    //- The caller (two-pass loop) checks refinementValid() and rejects
    //- this pass, restoring the previous valid mesh.
    if( !refinementValid_ )
    {
        WarningIn("void refineBoundaryLayers::refineLayers()")
            << "Boundary-layer refinement metadata invalid -- "
            << "aborting refinement attempt before cell generation"
            << endl;
        return;
    }

    generateNewCells();

    // CFMitch V3.5: generateNewCells() may reject a structurally invalid
    // prism parent.  Stop immediately so the owning mesh transaction can
    // restore its previous snapshot.
    if( !refinementValid_ )
    {
        WarningIn("void refineBoundaryLayers::refineLayers()")
            << "Boundary-layer cell generation rejected -- "
            << "leaving refinementCompleted=false for rollback"
            << endl;

        return;
    }

    done_ = true;

    Info << "Finished refining boundary layers" << endl;
}

void refineBoundaryLayers::pointsInBndLayer(labelLongList& layerPoints)
{
    layerPoints.clear();

    boolList pointInLayer(mesh_.points().size(), false);

    forAll(newVerticesForSplitEdge_, seI)
    {
        forAllRow(newVerticesForSplitEdge_, seI, i)
            pointInLayer[newVerticesForSplitEdge_(seI, i)] = true;
    }

    forAll(pointInLayer, pointI)
        if( pointInLayer[pointI] )
            layerPoints.append(pointI);
}

void refineBoundaryLayers::pointsInBndLayer(const word subsetName)
{
    label sId = mesh_.pointSubsetIndex(subsetName);
    if( sId < 0 )
        sId = mesh_.addPointSubset(subsetName);

    forAll(newVerticesForSplitEdge_, seI)
    {
        forAllRow(newVerticesForSplitEdge_, seI, i)
            mesh_.addPointToSubset(sId, newVerticesForSplitEdge_(seI, i));
    }
}

void refineBoundaryLayers::readSettings
(
    const dictionary& meshDict,
    refineBoundaryLayers& refLayers
)
{
    if( meshDict.isDict("boundaryLayers") )
    {
        const dictionary& bndLayers = meshDict.subDict("boundaryLayers");

        // CFMitch boundary-layer architecture selector.
        // Default remains legacyEnhanced for backward compatibility.
        if( bndLayers.found("boundaryLayerArchitecture") )
        {
            word architectureName;
            bndLayers.lookup("boundaryLayerArchitecture")
                >> architectureName;

            refLayers.boundaryLayerArchitecture_ =
                architectureName;
        }

        // CFMitch planner-specific settings.
        if( bndLayers.isDict("constraintPlanner") )
        {
            const dictionary& plannerDict =
                bndLayers.subDict("constraintPlanner");

            if( plannerDict.found("maxLayerStep") )
            {
                const label maxLayerStep =
                    readLabel
                    (
                        plannerDict.lookup("maxLayerStep")
                    );

                if( maxLayerStep < 1 )
                {
                    FatalErrorIn
                    (
                        "refineBoundaryLayers::readSettings"
                    )
                        << "boundaryLayers.constraintPlanner."
                        << "maxLayerStep must be >= 1"
                        << exit(FatalError);
                }

                refLayers.constraintPlannerMaxLayerStep_ =
                    maxLayerStep;
            }
        }

        // CFMitch V5.2a: read-only birth-stage marching contract census.
        // Default OFF.  It is deliberately a boundaryLayers setting since
        // the owning refiner, rather than cartesianMeshGenerator, performs
        // the prospective evaluation.
        if( bndLayers.found("cfmitchV52aFrontCensus") )
        {
            refLayers.cfmitchV52aFrontCensus_ =
                bool
                (
                    Switch
                    (
                        bndLayers.lookup("cfmitchV52aFrontCensus")
                    )
                );
        }

        if( bndLayers.found("cfmitchV52bApplyFrontCaps") )
        {
            refLayers.cfmitchV52bApplyFrontCaps_ =
                bool
                (
                    Switch
                    (
                        bndLayers.lookup
                        (
                            "cfmitchV52bApplyFrontCaps"
                        )
                    )
                );

            if( refLayers.cfmitchV52bApplyFrontCaps_ )
                refLayers.cfmitchV52aFrontCensus_ = true;
        }

        //- read global properties
        if( bndLayers.found("nLayers") )
        {
            const label nLayers = readLabel(bndLayers.lookup("nLayers"));
            refLayers.setGlobalNumberOfLayers(nLayers);
        }
        if( bndLayers.found("thicknessRatio") )
        {
            const scalar ratio =
                readScalar(bndLayers.lookup("thicknessRatio"));
            refLayers.setGlobalThicknessRatio(ratio);
        }
        if( bndLayers.found("maxFirstLayerThickness") )
        {
            const scalar maxFirstThickness =
                readScalar(bndLayers.lookup("maxFirstLayerThickness"));
            refLayers.setGlobalMaxThicknessOfFirstLayer(maxFirstThickness);
        }

        //- consider specified patches for exclusion from boundary layer creation
        //- implemented by setting the number of layers to 1
        if( bndLayers.found("excludedPatches") )
        {
            const wordList patchNames(bndLayers.lookup("excludedPatches"));

            forAll(patchNames, patchI)
            {
                const word pName = patchNames[patchI];

                refLayers.setNumberOfLayersForPatch(pName, 1);
            }
        }

        //- patch-based properties
        if( bndLayers.isDict("patchBoundaryLayers") )
        {
            const dictionary& patchBndLayers =
                bndLayers.subDict("patchBoundaryLayers");
            const wordList patchNames = patchBndLayers.toc();

            forAll(patchNames, patchI)
            {
                const word pName = patchNames[patchI];

                if( patchBndLayers.isDict(pName) )
                {
                    const dictionary& patchDict =
                        patchBndLayers.subDict(pName);

                    if( patchDict.found("nLayers") )
                    {
                        const label nLayers =
                            readLabel(patchDict.lookup("nLayers"));

                        refLayers.setNumberOfLayersForPatch(pName, nLayers);
                    }
                    if( patchDict.found("thicknessRatio") )
                    {
                        const scalar ratio =
                            readScalar(patchDict.lookup("thicknessRatio"));
                        refLayers.setThicknessRatioForPatch(pName, ratio);
                    }
                    if( patchDict.found("maxFirstLayerThickness") )
                    {
                        const scalar maxFirstThickness =
                            readScalar
                            (
                                patchDict.lookup("maxFirstLayerThickness")
                            );
                        refLayers.setMaxThicknessOfFirstLayerForPatch
                        (
                            pName,
                            maxFirstThickness
                        );
                    }
                    if( patchDict.found("allowDiscontinuity") )
                    {
                        const bool allowDiscontinuity =
                            readBool(patchDict.lookup("allowDiscontinuity"));

                        if( allowDiscontinuity )
                            refLayers.setInteruptForPatch(pName);
                    }
                }
                else
                {
                    Warning << "Cannot refine layer for patch "
                        << patchNames[patchI] << endl;
                }
            }
        }
    }
    else
    {
        //- the layer will not be refined
        refLayers.avoidRefinement();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void refineBoundaryLayers::setQualityMaxLayersAtFaces
(
    const Map<label>& caps
)
{
    label nImported = 0;
    label nAdded = 0;
    label nLowered = 0;
    label nZero = 0;
    label nInvalid = 0;

    forAllConstIter(Map<label>, caps, it)
    {
        const label bfI = it.key();
        const label requestedCap = it();

        if( bfI < 0 || requestedCap < 0 )
        {
            ++nInvalid;
            continue;
        }

        ++nImported;

        label effectiveCap = requestedCap;

        if( qualityMaxLayersAtFace_.found(bfI) )
        {
            const label oldCap =
                qualityMaxLayersAtFace_[bfI];

            effectiveCap =
                Foam::min(oldCap, requestedCap);

            if( effectiveCap < oldCap )
                ++nLowered;

            qualityMaxLayersAtFace_[bfI] =
                effectiveCap;
        }
        else
        {
            qualityMaxLayersAtFace_.insert
            (
                bfI,
                effectiveCap
            );

            ++nAdded;
        }

        // Merge directly into the existing final cap field.
        //
        // analyseLayers() already accepts zero as a valid cap and
        // applies forcedMaxLayersAtFace_ before CFMitch resolves
        // layer-count compatibility.
        if( forcedMaxLayersAtFace_.found(bfI) )
        {
            forcedMaxLayersAtFace_[bfI] =
                Foam::min
                (
                    forcedMaxLayersAtFace_[bfI],
                    effectiveCap
                );
        }
        else
        {
            forcedMaxLayersAtFace_.insert
            (
                bfI,
                effectiveCap
            );
        }

        if( effectiveCap == 0 )
            ++nZero;
    }

    Info
        << "CFMITCH V2.7 QUALITY CAP IMPORT:"
        << " requested=" << caps.size()
        << " imported=" << nImported
        << " added=" << nAdded
        << " lowered=" << nLowered
        << " zero=" << nZero
        << " invalid=" << nInvalid
        << " effectiveQualityFaces="
        << qualityMaxLayersAtFace_.size()
        << endl;
}


void refineBoundaryLayers::forceSingleLayerAtFaces
(
    const labelHashSet& faces
)
{
    // Backwards-compatible wrapper: old behavior means seed faces capped
    // to exactly one layer and no neighbour-ring expansion.
    forceMaxLayersAtFaces(faces, 1, 0, 0);
}

void refineBoundaryLayers::forceMaxLayersAtFaces
(
    const labelHashSet& faces,
    const label ring0MaxLayers,
    const label ring1MaxLayers,
    const label ring2MaxLayers,
    const scalar ring0ThicknessScale,
    const scalar ring1ThicknessScale,
    const scalar ring2ThicknessScale,
    const label sourcePlanId
)
{
    const meshSurfaceEngine& mse = surfaceEngine();
    const VRWGraph& faceFaces = mse.faceFaces();
    const label nBndFaces = faceFaces.size();

    labelHashSet ring0;
    labelHashSet ring1;
    labelHashSet ring2;

    forAllConstIter(labelHashSet, faces, it)
    {
        const label bfI = it.key();
        if( bfI >= 0 && bfI < nBndFaces )
            ring0.insert(bfI);
    }

    // CFMitch V4.1 -- thickness-only repair requests must also
    // build their smoothing rings.  Ring 1 is additionally required
    // whenever ring 2 is requested.
    const bool ring2Requested =
        ring2MaxLayers > 0
     || ring2ThicknessScale < scalar(1.0) - SMALL;

    const bool ring1Requested =
        ring1MaxLayers > 0
     || ring1ThicknessScale < scalar(1.0) - SMALL
     || ring2Requested;

    if( ring1Requested )
    {
        forAllConstIter(labelHashSet, ring0, it)
        {
            const label bfI = it.key();
            forAllRow(faceFaces, bfI, nI)
            {
                const label nbfI = faceFaces(bfI, nI);
                if( nbfI < 0 || nbfI >= nBndFaces ) continue;
                if( ring0.found(nbfI) ) continue;
                ring1.insert(nbfI);
            }
        }
    }

    if( ring2Requested )
    {
        forAllConstIter(labelHashSet, ring1, it)
        {
            const label bfI = it.key();
            forAllRow(faceFaces, bfI, nI)
            {
                const label nbfI = faceFaces(bfI, nI);
                if( nbfI < 0 || nbfI >= nBndFaces ) continue;
                if( ring0.found(nbfI) || ring1.found(nbfI) ) continue;
                ring2.insert(nbfI);
            }
        }
    }

    label nAdded = 0;
    label nLowered = 0;
    label nThicknessAdded = 0;
    label nThicknessLowered = 0;

    // Track which source plan supplied the final minimum cap.
    // sourcePlanId==-1 identifies a non-BLRepairPlan caller.
    auto updateCapProviders =
    [&](const label bfI, const bool existed,
        const label oldCap, const label requestedCap)
    {
        // Provider provenance follows the REQUEST that supplied the
        // winning minimum cap, not min(oldCap, requestedCap).
        //
        // requested < old: stricter request becomes sole winner
        // requested == old: equal winning request shares provenance
        // requested > old: weaker request contributes nothing
        if( !existed || requestedCap < oldCap )
        {
            if( forcedMaxLayerProviderPlansAtFace_.found(bfI) )
            {
                forcedMaxLayerProviderPlansAtFace_[bfI].clear();
                forcedMaxLayerProviderPlansAtFace_[bfI].insert(sourcePlanId);
            }
            else
            {
                labelHashSet providers;
                providers.insert(sourcePlanId);
                forcedMaxLayerProviderPlansAtFace_.insert(bfI, providers);
            }
        }
        else if( requestedCap == oldCap )
        {
            if( !forcedMaxLayerProviderPlansAtFace_.found(bfI) )
            {
                labelHashSet providers;
                providers.insert(sourcePlanId);
                forcedMaxLayerProviderPlansAtFace_.insert(bfI, providers);
            }
            else
            {
                forcedMaxLayerProviderPlansAtFace_[bfI].insert(sourcePlanId);
            }
        }
    };

    auto applyThicknessScale =
    [&](const labelHashSet& ringFaces, const scalar scale)
    {
        if( scale >= scalar(1.0) - SMALL )
            return;

        const scalar clippedScale =
            Foam::max(scalar(0.0), Foam::min(scalar(1.0), scale));

        forAllConstIter(labelHashSet, ringFaces, it)
        {
            const label bfI = it.key();

            if( forcedThicknessScaleAtFace_.found(bfI) )
            {
                const scalar oldScale = forcedThicknessScaleAtFace_[bfI];
                const scalar newScale = Foam::min(oldScale, clippedScale);

                if( newScale < oldScale - SMALL )
                    ++nThicknessLowered;

                forcedThicknessScaleAtFace_[bfI] = newScale;
            }
            else
            {
                forcedThicknessScaleAtFace_.insert(bfI, clippedScale);
                ++nThicknessAdded;
            }
        }
    };

    forAllConstIter(labelHashSet, ring0, it)
    {
        const label bfI = it.key();
        if( ring0MaxLayers <= 0 ) continue;
        if( forcedMaxLayersAtFace_.found(bfI) )
        {
            const label oldCap = forcedMaxLayersAtFace_[bfI];
            const label newCap = Foam::min(oldCap, ring0MaxLayers);
            if( newCap < oldCap ) ++nLowered;
            updateCapProviders(bfI, true, oldCap, ring0MaxLayers);
            forcedMaxLayersAtFace_[bfI] = newCap;
        }
        else
        {
            forcedMaxLayersAtFace_.insert(bfI, ring0MaxLayers);
            updateCapProviders
            (
                bfI,
                false,
                ring0MaxLayers,
                ring0MaxLayers
            );
            ++nAdded;
        }
    }

    forAllConstIter(labelHashSet, ring1, it)
    {
        const label bfI = it.key();
        if( ring1MaxLayers <= 0 ) continue;
        if( forcedMaxLayersAtFace_.found(bfI) )
        {
            const label oldCap = forcedMaxLayersAtFace_[bfI];
            const label newCap = Foam::min(oldCap, ring1MaxLayers);
            if( newCap < oldCap ) ++nLowered;
            updateCapProviders(bfI, true, oldCap, ring1MaxLayers);
            forcedMaxLayersAtFace_[bfI] = newCap;
        }
        else
        {
            forcedMaxLayersAtFace_.insert(bfI, ring1MaxLayers);
            updateCapProviders
            (
                bfI,
                false,
                ring1MaxLayers,
                ring1MaxLayers
            );
            ++nAdded;
        }
    }

    forAllConstIter(labelHashSet, ring2, it)
    {
        const label bfI = it.key();
        if( ring2MaxLayers <= 0 ) continue;
        if( forcedMaxLayersAtFace_.found(bfI) )
        {
            const label oldCap = forcedMaxLayersAtFace_[bfI];
            const label newCap = Foam::min(oldCap, ring2MaxLayers);
            if( newCap < oldCap ) ++nLowered;
            updateCapProviders(bfI, true, oldCap, ring2MaxLayers);
            forcedMaxLayersAtFace_[bfI] = newCap;
        }
        else
        {
            forcedMaxLayersAtFace_.insert(bfI, ring2MaxLayers);
            updateCapProviders
            (
                bfI,
                false,
                ring2MaxLayers,
                ring2MaxLayers
            );
            ++nAdded;
        }
    }

    applyThicknessScale(ring0, ring0ThicknessScale);
    applyThicknessScale(ring1, ring1ThicknessScale);
    applyThicknessScale(ring2, ring2ThicknessScale);

    Info << "refineBoundaryLayers: forceMaxLayersAtFaces: "
         << "seedFaces=" << faces.size()
         << " validRing0=" << ring0.size()
         << " ring1=" << ring1.size()
         << " ring2=" << ring2.size()
         << " caps=(" << ring0MaxLayers << ','
                     << ring1MaxLayers << ','
                     << ring2MaxLayers << ')'
         << " totalCappedFaces=" << forcedMaxLayersAtFace_.size()
         << " added=" << nAdded
         << " lowered=" << nLowered
         << " thicknessScales=(" << ring0ThicknessScale << ','
                                 << ring1ThicknessScale << ','
                                 << ring2ThicknessScale << ')'
         << " totalScaledFaces=" << forcedThicknessScaleAtFace_.size()
         << " thicknessAdded=" << nThicknessAdded
         << " thicknessLowered=" << nThicknessLowered
         << endl;
}

void refineBoundaryLayers::setBlblJunctionPoints
(
    const labelHashSet& pts
)
{
    blblJunctionPoints_ = pts;
    Info << "refineBoundaryLayers: received "
         << blblJunctionPoints_.size()
         << " BL/BL junction points" << endl;
}
void refineBoundaryLayers::setBlblAcuteCornerPoints
(
    const labelHashSet& pts
)
{
    blblAcuteCornerPoints_ = pts;
    Info << "refineBoundaryLayers: received "
         << blblAcuteCornerPoints_.size()
         << " acute BL+BL+neutral corner points" << endl;
}

void refineBoundaryLayers::setRampSeedPoints
(
    const boolList& pts
)
{
    rampSeedPoints_ = pts;
    label nRamp = 0;
    forAll(rampSeedPoints_, i)
        if( rampSeedPoints_[i] ) ++nRamp;
    Info << "refineBoundaryLayers: received "
         << nRamp
         << " finite ramp seed points (excluded from junction caps)" << endl;
}
void refineBoundaryLayers::setVtFaceRing
(
    const labelList& ring
)
{
    vtFaceRing_ = ring;
    label nHandled = 0;
    forAll(vtFaceRing_, fI)
        if( vtFaceRing_[fI] >= 0 ) ++nHandled;
    Info << "refineBoundaryLayers: received vtFaceRing, "
         << nHandled << " faces under virtual topology control" << endl;
}

} // End namespace Foam

// ************************************************************************* //
