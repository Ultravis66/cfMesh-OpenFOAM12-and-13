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
    nLayersAtBndFace_(),
    cellToBaseBndFace_(),
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
    const scalar ring2ThicknessScale
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

    if( ring1MaxLayers > 0 )
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

    if( ring2MaxLayers > 0 )
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
            forcedMaxLayersAtFace_[bfI] = newCap;
        }
        else
        {
            forcedMaxLayersAtFace_.insert(bfI, ring0MaxLayers);
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
            forcedMaxLayersAtFace_[bfI] = newCap;
        }
        else
        {
            forcedMaxLayersAtFace_.insert(bfI, ring1MaxLayers);
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
            forcedMaxLayersAtFace_[bfI] = newCap;
        }
        else
        {
            forcedMaxLayersAtFace_.insert(bfI, ring2MaxLayers);
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
