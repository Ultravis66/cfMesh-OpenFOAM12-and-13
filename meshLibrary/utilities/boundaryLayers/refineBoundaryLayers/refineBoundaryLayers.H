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

Class
    refineBoundaryLayers

Description
    Refine existing boundary layers

SourceFiles
    refineBoundaryLayers.C
    refineBoundaryLayersFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef refineBoundaryLayers_H
#define refineBoundaryLayers_H

#include "polyMeshGenModifier.H"
#include "meshSurfaceEngine.H"
#include "DynList.H"
#include "labelLongList.H"
#include "labelPair.H"

#include <map>
#include <set>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
                Class refineBoundaryLayers Declaration
\*---------------------------------------------------------------------------*/

class refineBoundaryLayers
{
        //- Reference to the mesh
        polyMeshGen& mesh_;

        //- pointer to mesh surface engine
        mutable meshSurfaceEngine* msePtr_;

        //- global number of boundary layers
        label globalNumLayers_;

        //- global thickness ratio
        scalar globalThicknessRatio_;

        //- global maximum thickness of the first layer
        scalar globalMaxThicknessFirstLayer_;

        //- number of boundary layers for user-selected patches
        std::map<word, label> numLayersForPatch_;

        //- local thickness ratio for selected patches
        std::map<word, scalar> thicknessRatioForPatch_;

        //- local maximum layer thickness for selected patches
        std::map<word, scalar> maxThicknessForPatch_;

        //- allow discontinuous layers for patch
        std::set<word> discontinuousLayersForPatch_;

        //- the name of a subset which shall contain the indices
        //- of the cells in the boundary layer
        word cellSubsetName_;

        //- check whether the refinement is already executed
        bool done_;

        //- a flag whether a 2D mesh generation is active or not
        bool is2DMesh_;

        //- shall the layer be refined into two layers with the first one
        //- matching the thickness of the n-1 ones
        bool specialMode_;

        //- information about existing boundary layers at patches
        //- only available layers
        List<DynList<label> > layerAtPatch_;

        //- which patches are part of a single layer
        List<DynList<word> > patchesInLayer_;

        //- a containing the number of layers which shall be generated above
        //- a boundary face
        labelList nLayersAtBndFace_;

        //- boundary faces forced to single layer treatment
        labelHashSet forceSingleLayerFaces_;
        //- gap conflict arbitration data for local split-edge layer capping.
        //- gapActionPoints_: mesh point labels of loser-side gap zone points.
        //- gapLoserPatchNames_: loser-side patch names.
        //- Caps are applied only where BOTH conditions are met (local + side).
        labelHashSet gapActionPoints_;
        wordList     gapLoserPatchNames_;
        label        gapRing1MaxLayers_;   // default 1
        label        gapRing2MaxLayers_;   // default 2

        //- BL/termination edge cap: reduce layer count at inlet/outlet junctions.
        //- blTerminationEdgePoints_: mesh point labels on BL/no-BL edges.
        //- No patch-name check needed -- locality is sufficient.
        labelHashSet blTerminationEdgePoints_;
        label        blTerminationRing1MaxLayers_;  // default 3 (disabled)
        label        blTerminationRing2MaxLayers_;  // default 3 (disabled)
        labelHashSet blblAcuteCornerPoints_;
        bool acuteCornerCapLayers_;

        //- BL/BL junction points from boundaryLayers - used for wedge topology
        labelHashSet blblJunctionPoints_;

        //- Virtual topology face ring from boundaryLayers
        //- -1=untouched, 0/1/2=ring. BL/BL ramp skips handled faces.
        labelList vtFaceRing_;

        //- a list of edges which shall be refined
        LongList<edge> splitEdges_;

        //- split edges at point
        VRWGraph splitEdgesAtPoint_;

        //- new vertices for on edges which shall be refined
        VRWGraph newVerticesForSplitEdge_;

        //- a graph containing information which new faces were generated
        //- from an existing face
        VRWGraph facesFromFace_;

        //- a graph containing faces after layer refinement
        VRWGraph newFaces_;

    // Private member functions
        //- Return reference to meshSurfaceEngine
        const meshSurfaceEngine& surfaceEngine() const;

        //- analyse layers to check their topology
        bool analyseLayers();

        //- generate new points on edges, faces and in cells
        void generateNewVertices();

        //- refine a given face and return the new faces
        //- generates new points at cross-split faces
        void refineFace
        (
            const face& f,
            const FixedList<label, 2>& nLayersInDirection,
            DynList<DynList<label, 4>, 128>& newFaces
        );

        //- generate a matrix of points generated by splitting a face
        //- and return them in the local i, j system of the face
        void sortFacePoints
        (
            const label faceI,
            DynList<DynList<label> >& facePoints,
            const label transpose = false
        ) const;

        //- generate a matrix of faces generated by splitting a face
        //- and return them in the local i, j, system of the face
        void sortFaceFaces
        (
            const label faceI,
            DynList<DynList<label> >& faceFaces,
            const label transpose = false
        ) const;

        //- map split edges onto a cell
        void generateNewFaces();

        //- generate new cells for a prism with one boundary face
        void generateNewCellsPrism
        (
            const label cellI,
            DynList<DynList<DynList<label, 8>, 10>, 64>& cellsFromCell
        ) const;

        //- a helper function which stores faces generated from
        //- an existing face into new cells
        void storeFacesIntoCells
        (
            const label faceI,
            const bool reverseOrientation,
            const label normalDirection,
            const bool maxCoordinate,
            const label nLayersI,
            const label nLayersJ,
            const label nLayersK,
            DynList<DynList<DynList<label, 4>, 6>, 256>& cellsFromCell
        ) const;

        //- generate new cells and add them to the mesh
        void generateNewCells();

    // Nested classes
        class refineEdgeHexCell
        {
            // Private data
                //- label of cell
                const label cellI_;

                //- number of cells in local direction i
                label nLayersI_;

                //- number of cells in locatiol direction j
                label nLayersJ_;

                //- container for new cells
                DynList<DynList<DynList<label, 4>, 6>, 256> cellsFromCell_;

                //- const reference to the boundary layer class
                const refineBoundaryLayers& bndLayers_;

                //- faces sorted into directions of a hex shape
                FixedList<label, 6> faceInDirection_;

                //- information about orientation of faces
                //- false means the orientation as expected
                //- true means wrong orientation
                FixedList<bool, 6> faceOrientation_;

                //- points on cross-split faces
                FixedList<DynList<DynList<label> >, 2> cellPoints_;

            // Private member functions
                //- populate faceInDirection_nad wrongFaceOrientation_
                void determineFacesInDirections();

                //- populate new cells with new faces generated from already
                //- existing faces
                void populateExistingFaces();

                //- generate new internal faces and tore them to new cells
                void generateMissingFaces();

        public:

            // Constructor
                //- construct from cell label and the refineBoundaryLayers
                refineEdgeHexCell
                (
                    const label cellI,
                    const refineBoundaryLayers& ref
                );

            // Public member functions
                inline const DynList<DynList<DynList<label, 4>, 6>, 256>&
                newCells() const
                {
                    return cellsFromCell_;
                }
        };

        class refineCornerHexCell
        {
            // Private data
                //- label of cell
                const label cellI_;

                //- number of cells in local direction i
                label nLayersI_;

                //- number of cells in local direction j
                label nLayersJ_;

                //- number of cells in local direction k
                label nLayersK_;

                //- split edge in directions
                FixedList<label, 3> splitEdgeInDirection_;

                //- container for new cells
                DynList<DynList<DynList<label, 4>, 6>, 256> cellsFromCell_;

                //- reference to the boundary layer class
                const refineBoundaryLayers& bndLayers_;

                //- faces sorted into directions of a hex shape
                FixedList<label, 6> faceInDirection_;

                //- information about orientation of faces
                //- false means the orientation as expected
                //- true means wrong orientation
                FixedList<bool, 6> faceOrientation_;

                //- points on cross-split faces
                FixedList<DynList<DynList<label> >, 6> facePoints_;

                //- points inside the cell
                DynList<DynList<DynList<label> > > cellPoints_;

            // Private member functions
                //- populate faceInDirection_nad wrongFaceOrientation_
                void determineFacesInDirections();

                //- populate new cells with new faces generated from already
                //- existing faces
                void populateExistingFaces();

                //- generate missing points inside the cell
                void generateNewPoints();

                //- generate new internal faces and tore them to new cells
                void generateMissingFaces();

        public:

            // Constructor
                //- construct from cell label and the refineBoundaryLayers
                refineCornerHexCell
                (
                    const label cellI,
                    const refineBoundaryLayers& ref
                );

            // Public member functions
                inline const DynList<DynList<DynList<label, 4>, 6>, 256>&
                newCells() const
                {
                    return cellsFromCell_;
                }
        };

    // Private member functions

        //- Disallow bitwise copy construct
        refineBoundaryLayers(const refineBoundaryLayers&);

        //- Disallow bitwise assignment
        void operator=(const refineBoundaryLayers&);

public:

    // Constructors

        //- Construct from mesh reference
        refineBoundaryLayers(polyMeshGen& mesh);

    // Destructor
        ~refineBoundaryLayers();

    // Public member functions

        // Settings
        //- set no refinement flag
        void avoidRefinement();

        //- activate 2D layer refinement
        void activate2DMode();

        //- set the global number of boundary layers
        void setGlobalNumberOfLayers(const label nLayers);

        //- set the global thickness ratio (default is 1)
        void setGlobalThicknessRatio(const scalar thicknessRatio);

        //- set the maximum thickness of the first boundary layer
        void setGlobalMaxThicknessOfFirstLayer(const scalar maxThickness);

        //- set the number of layers for a patch
        //- the settings override the global settings
        void setNumberOfLayersForPatch
        (
            const word& patchName,
            const label nLayers
        );

        //- set the thickness ratio for a patch
        //- it overrides the global settings
        void setThicknessRatioForPatch
        (
            const word& patchName,
            const scalar thicknessRatio
        );

        //- set the maximum thickness of the first layer for a patch
        void setMaxThicknessOfFirstLayerForPatch
        (
            const word& patchName,
            const scalar maxThickness
        );

        //- set whether the settings for a given patch are valid for the
        //- patch only, or whether they extend over th whole sheet
        //- the selected patch belongs to
        //- the default behaviour is to apply the patch settings to the whole
        //- sheet
        void setInteruptForPatch(const word& patchName);

        //- set the flag to store the boundary layer cells in a subset
        void setCellSubset(const word subsetName = "layerCells");

        //- this flag refinement of a layer into a single one matching
        //- the thickness of the n layers. This functionality is intended to
        //- create additional layer that shall be improved by the boundary
        //- layer smoother in order to leave the inside points intact
        void activateSpecialMode();

    // Worker
        //- performs refinement based on the given settings
        void refineLayers();

    // Post information
        //- provide the list of points in the boundary layer
        void pointsInBndLayer(labelLongList&);

        //- create a subset containing points in the boundary layer
        void pointsInBndLayer(const word subsetName="pointsInBndLayer");

    // Static member functions
        //- read the settings from dictionary
        static void readSettings(const dictionary&, refineBoundaryLayers&);

        //- force single layer on mixed-displacement boundary faces
        void forceSingleLayerAtFaces(const labelHashSet& faces);

        //- set BL/BL junction points for wedge topology emission
        void setBlblJunctionPoints(const labelHashSet& pts);
        void setBlblAcuteCornerPoints(const labelHashSet& pts);
        void setVtFaceRing(const labelList& ring);
        void setAcuteCornerCapLayers(const bool cap) { acuteCornerCapLayers_ = cap; }
        //- set gap conflict arbitration data for local split-edge layer capping
        void setGapActionPoints(const labelHashSet& pts)
        {
            gapActionPoints_ = pts;
        }
        void setGapLoserPatchNames(const wordList& names)
        {
            gapLoserPatchNames_ = names;
        }
        void setGapRingMaxLayers(const label ring1Max, const label ring2Max)
        {
            gapRing1MaxLayers_ = ring1Max;
            gapRing2MaxLayers_ = ring2Max;
        }
        //- set BL/termination edge points for inlet/outlet layer-count cap
        void setBlTerminationEdgePoints(const labelHashSet& pts)
        {
            blTerminationEdgePoints_ = pts;
        }
        void setBlTerminationRingMaxLayers(const label ring1Max, const label ring2Max)
        {
            blTerminationRing1MaxLayers_ = ring1Max;
            blTerminationRing2MaxLayers_ = ring2Max;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
