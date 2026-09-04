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

#include "classicRefineBoundaryLayers.H"
#include "meshSurfaceEngine.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const meshSurfaceEngine& classicRefineBoundaryLayers::surfaceEngine() const
{
    if( !msePtr_ )
        msePtr_ = new meshSurfaceEngine(mesh_);

    return *msePtr_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

classicRefineBoundaryLayers::classicRefineBoundaryLayers(polyMeshGen& mesh)
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
    nLayersAtBndFace_(),
    splitEdges_(),
    splitEdgesAtPoint_(),
    newVerticesForSplitEdge_(),
    facesFromFace_(),
    newFaces_()
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

classicRefineBoundaryLayers::~classicRefineBoundaryLayers()
{
    deleteDemandDrivenData(msePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void classicRefineBoundaryLayers::avoidRefinement()
{
    if( done_ )
    {
        FatalErrorIn
        (
            "void classicRefineBoundaryLayers::avoidRefinement()"
        ) << "refineLayers is already executed" << exit(FatalError);
    }

    globalNumLayers_ = 1;
    numLayersForPatch_.clear();
}

void classicRefineBoundaryLayers::activate2DMode()
{
    if( done_ )
    {
        FatalErrorIn
        (
            "void classicRefineBoundaryLayers::activate2DMode()"
        ) << "refineLayers is already executed" << exit(FatalError);
    }

    is2DMesh_ = true;
}

void classicRefineBoundaryLayers::setGlobalNumberOfLayers(const label nLayers)
{
    if( done_ )
    {
        FatalErrorIn
        (
            "void classicRefineBoundaryLayers::setGlobalNumberOfLayers(const label)"
        ) << "refineLayers is already executed" << exit(FatalError);
    }

    if( nLayers < 2 )
    {
        WarningIn
        (
            "void classicRefineBoundaryLayers::setGlobalNumberOfLayers(const label)"
        ) << "The specified global number of boundary layers is less than 2"
          << endl;

        return;
    }

    globalNumLayers_ = nLayers;
}

void classicRefineBoundaryLayers::setGlobalThicknessRatio(const scalar thicknessRatio)
{
    if( done_ )
    {
        FatalErrorIn
        (
            "void classicRefineBoundaryLayers::setGlobalThicknessRatio(const scalar)"
        ) << "refineLayers is already executed" << exit(FatalError);
    }

    if( thicknessRatio < 1.0 )
    {
        WarningIn
        (
            "void classicRefineBoundaryLayers::setGlobalThicknessRatio(const scalar)"
        ) << "The specified global thickness ratio is less than 1.0" << endl;

        return;
    }

    globalThicknessRatio_ = thicknessRatio;
}

void classicRefineBoundaryLayers::setGlobalMaxThicknessOfFirstLayer
(
    const scalar maxThickness
)
{
    if( done_ )
    {
        FatalErrorIn
        (
            "void classicRefineBoundaryLayers::setGlobalMaxThicknessOfFirstLayer"
            "(const scalar)"
        ) << "refineLayers is already executed" << exit(FatalError);
    }

    if( maxThickness <= 0.0 )
    {
        WarningIn
        (
            "void classicRefineBoundaryLayers::setGlobalMaxThicknessOfFirstLayer"
            "(const scalar)"
        ) << "The specified global maximum thickness of the first"
          << " boundary layer is negative!!" << endl;

        return;
    }

    globalMaxThicknessFirstLayer_ = maxThickness;
}

void classicRefineBoundaryLayers::setNumberOfLayersForPatch
(
    const word& patchName,
    const label nLayers
)
{
    if( done_ )
    {
        FatalErrorIn
        (
            "void classicRefineBoundaryLayers::setNumberOfLayersForPatch"
            "(const word&, const label)"
        ) << "refineLayers is already executed" << exit(FatalError);
    }

    if( nLayers < 2 )
    {
        WarningIn
        (
            "void classicRefineBoundaryLayers::setNumberOfLayersForPatch"
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

void classicRefineBoundaryLayers::setThicknessRatioForPatch
(
    const word& patchName,
    const scalar thicknessRatio
)
{
    if( done_ )
    {
        FatalErrorIn
        (
            "void classicRefineBoundaryLayers::setThicknessRatioForPatch"
            "(const word&, const scalar)"
        ) << "refineLayers is already executed" << exit(FatalError);
    }

    if( thicknessRatio < 1.0 )
    {
        WarningIn
        (
            "void classicRefineBoundaryLayers::setThicknessRatioForPatch"
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

void classicRefineBoundaryLayers::setMaxThicknessOfFirstLayerForPatch
(
    const word& patchName,
    const scalar maxThickness
)
{
    if( done_ )
    {
        FatalErrorIn
        (
            "void classicRefineBoundaryLayers::setMaxThicknessOfFirstLayerForPatch"
            "(const word&, const scalar)"
        ) << "refineLayers is already executed" << exit(FatalError);
    }

    if( maxThickness <= 0.0 )
    {
        WarningIn
        (
            "void classicRefineBoundaryLayers::setGlobalMaxThicknessOfFirstLayer"
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

void classicRefineBoundaryLayers::setInteruptForPatch(const word& patchName)
{
    if( done_ )
    {
        FatalErrorIn
        (
            "void classicRefineBoundaryLayers::setInteruptForPatch(const word&)"
        ) << "refineLayers is already executed" << exit(FatalError);
    }

    const labelList matchedIDs = mesh_.findPatches(patchName);

    forAll(matchedIDs, matchI)
    {
        const word pName = mesh_.getPatchName(matchedIDs[matchI]);
        discontinuousLayersForPatch_.insert(pName);
    }
}

void classicRefineBoundaryLayers::setCellSubset(const word subsetName)
{
    if( done_ )
    {
        FatalErrorIn
        (
            "void classicRefineBoundaryLayers::setCellSubset(const word)"
        ) << "refineLayers is already executed" << exit(FatalError);
    }

    cellSubsetName_ = subsetName;
}

void classicRefineBoundaryLayers::activateSpecialMode()
{
    specialMode_ = true;
}

void classicRefineBoundaryLayers::refineLayers()
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
            "void classicRefineBoundaryLayers::refineLayers()"
        ) << "Boundary layers are already refined! "
          << "Stopping refinement" << endl;

        return;
    }

    if( !analyseLayers() )
    {
        WarningIn
        (
            "void classicRefineBoundaryLayers::refineLayers()"
        ) << "Boundary layers do not exist in the mesh! Cannot refine" << endl;

        return;
    }

    generateNewVertices();

    generateNewFaces();

    generateNewCells();

    done_ = true;

    Info << "Finished refining boundary layers" << endl;
}

void classicRefineBoundaryLayers::pointsInBndLayer(labelLongList& layerPoints)
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

void classicRefineBoundaryLayers::pointsInBndLayer(const word subsetName)
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

void classicRefineBoundaryLayers::readSettings
(
    const dictionary& meshDict,
    classicRefineBoundaryLayers& refLayers
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

} // End namespace Foam

// ************************************************************************* //
