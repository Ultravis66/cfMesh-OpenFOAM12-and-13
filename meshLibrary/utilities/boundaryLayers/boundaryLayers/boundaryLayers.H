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
    boundaryLayers

Description
    Creates boundary layers

SourceFiles
    boundaryLayers.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryLayers_H
#define boundaryLayers_H

#include "polyMeshGenModifier.H"
#include "meshSurfaceEngine.H"
#include "meshSurfacePartitioner.H"
#include "DynList.H"
#include "labelLongList.H"
#include "Map.H"
#include "labelPair.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshSurfaceEngine;
class meshSurfacePartitioner;

/*---------------------------------------------------------------------------*\
                           Class boundaryLayers Declaration
\*---------------------------------------------------------------------------*/

class boundaryLayers
{
        //- Reference to the mesh
        polyMeshGen& mesh_;

        //- pointer to mesh surface engine
        mutable meshSurfaceEngine* msePtr_;

        //- poiner to meshSurfacePartitioner
        mutable meshSurfacePartitioner* meshPartitionerPtr_;

        //- shall I create patch-wise layers (true)
        //- O-topology layer (false)
        bool patchWiseLayers_;

        //- shall the layers be terminated at concave edges (true)
        bool terminateLayersAtConcaveEdges_;
        mutable boolList zeroDistPoints_;
        mutable scalarField layerScale_;

        //- meshDict-configurable BL transition parameters
        scalar blblFeatureAngleDeg_;   // degrees, default 40
        scalar blblCornerAcuteThreshold_; // dot product threshold, default 0.3
        scalar layerScaleRing1_;       // default 0.25
        scalar layerScaleRing2_;       // default 0.50
        scalar layerScaleRing3_;       // default 0.75
        scalar layerScaleRing4_;       // default 1.00
        scalar layerScaleRing5_;       // default 0.80
        scalar layerScaleRing6_;       // default 1.00
        // Acute corner local taper — stronger suppression seeded from blblAcuteCornerPoints_
        scalar acuteCornerRing0_;      // default 0.0
        scalar acuteCornerRing1_;      // default 0.0
        scalar acuteCornerRing2_;      // default 0.05
        scalar acuteCornerRing3_;      // default 0.15
        scalar acuteCornerRing4_;      // default 0.35
        scalar acuteCornerRing5_;      // default 0.60
        scalar acuteCornerRing6_;      // default 1.00
        // Virtual topology BL exclusion/taper around pathological
        // acute triple-junctions. Internal mesher view only:
        // no STL mutation, no patch renaming, no BL graph surgery.
        bool virtualTopologyExclusion_; // default false
        scalar virtualTopoRing0_;       // default 0.0  full exclusion
        scalar virtualTopoRing1_;       // default 0.0  full exclusion
        scalar virtualTopoRing2_;       // default 0.05 taper onset

        // Gap face-ring exclusion controls
        bool   gapFaceRingExclusion_;  // default true
        bool   tripleJunctionFaceRingExclusion_;  // default false
        scalar gapFaceRing0Scale_;     // default 0.02
        scalar gapFaceRing1Scale_;     // default 0.05
        scalar gapFaceRing2Scale_;     // default 0.20

        //- is it a 2D mesh
        bool is2DMesh_;

        //- patch names
        wordList patchNames_;

        //- patch types
        wordList patchTypes_;

        //- helper which contains information if a boundary layer
        //- has already been extruded for a given patch
        boolList treatedPatch_;

        //- extrude patches with patch
        List<DynList<label> > treatPatchesWithPatch_;

        //- label of a new node (helper)
        labelLongList newLabelForVertex_;

        //- map storing labels of new vertices created at node
        //- and corner vertices
        std::map<label, std::map<std::pair<label, label>, label> > otherVrts_;

        //- a key assigned to each patch. It is needed to search in otherVrts_
        labelList patchKey_;

        //- number of BL layers per patch (0 = no BL = termination patch)
        labelList nLayersForPatch_;

        //- patch role classification
        //  0 = BL_PATCH, 1 = TERMINATION_PATCH, 2 = NEUTRAL_PATCH
        labelList patchRole_;

        //- boundary faces with mixed zero/nonzero layer displacement
        //- passed to refineBoundaryLayers for single-layer treatment
        mutable labelHashSet mixedDisplacementFaces_;

        //- boundary points on sharp BL/BL junction feature curves
        //- populated by markConcaveEdgePoints, consumed by C2/C3
        mutable labelHashSet blblJunctionPoints_;
        mutable labelHashSet blblCornerPoints_;
        mutable labelHashSet blblAcuteCornerPoints_;

        //- Per-face virtual topology ring assignment (-1=untouched, 0/1/2=ring)
        //- Populated by applyVirtualTopologyExclusion, consumed by refineBoundaryLayers
        mutable labelList vtFaceRing_;

        //- gap points for BL suppression in thin clearance regions
        labelHashSet gapPoints_;
        //- triple-junction seed points (blade/hub-shroud/periodic singularities)
        labelHashSet tripleJunctionPoints_;

        //- face-level BL suppression mask — set by gap/VT detection
        //- createNewFacesAndCells skips faces where this is true
        mutable boolList suppressLayerAtBndFace_;

        //- edges where at least one adjacent patch requests BL and
        //- at least one adjacent patch requests zero layers.
        //- populated by markConcaveEdgePoints, consumed by surface mapper.
        mutable labelHashSet blNoBlEdges_;

        //- boundary points on BL/no-BL transition edges
        mutable labelHashSet blNoBlEdgePoints_;

        //- patch index for each BL/no-BL interface point
        //- stores the BL-side patch index for patch-constrained projection
        mutable Map<label> blNoBlPointPatch_;

        //- BL/neutral edge points (BL patch meets periodic/symmetry)
        mutable labelHashSet blNeutralEdgePoints_;

        //- BL-side patch index per neutral edge point
        mutable Map<label> blNeutralPointPatch_;

        //- number of vertices in the mesh
        label nPoints_;

        //- has the geometry been analysed
        bool geometryAnalysed_;

    // Private member functions
        //- Return const reference to meshSurfaceEngine
        const meshSurfaceEngine& surfaceEngine() const;

        //- return const reference to meshSurfacePartitioner
        const meshSurfacePartitioner& surfacePartitioner() const;

        //- find if any other patches need to be treated together
        //- with the given one
        void findPatchesToBeTreatedTogether();

        //- find vertices of the selected patches
        void findPatchVertices
        (
            const boolList& treatPatches,
            List<direction>& patchVertex
        ) const;

        //- check and correct the topology of boundary faces where
        //- the layers terminate
        void checkTopologyOfBoundaryFaces(const labelList& patchLabels);

        //- create new vertex
        point createNewVertex
        (
            const label bpI,
            const boolList& treatPatches,
            const List<direction>& pVertices
        ) const;

        //- create new vertices for the selected patches
        void createNewVertices
        (
            const boolList& treatPatches
        );

        //- create new vertices such that layers for selected patches
        //- are generated in a single run
        void createNewVertices
        (
            const labelList& patchLabels
        );

        //- creates new vertices for vertices at parallel boundaries
        void createNewPartitionVerticesParallel
        (
            const labelLongList& procPoints,
            const List<direction>& pVertices,
            const boolList& treatPatches
        );

        //- creates new vertices from vertices at parallel boundaries
        //- which are also at the border of the treated partitions
        void createNewEdgeVerticesParallel
        (
            const labelLongList& procPoints,
            const List<direction>& pVertices,
            const boolList& treatPatches
        );

        //- create a layer of cells
        void createNewFacesAndCells
        (
            const boolList& treatPatches
        );

        //- create new faces at parallel boundaries
        //- faces are extruded from edges
        void createNewFacesParallel
        (
            const boolList& treatPatches
        );

        //- create new faces at parallel boundaries
        //- faces are created from points at parallel boundaries
        //- the function takes a reference to the faces which are the candidates
        //- to create faces at parallel boundaries
        void createNewFacesFromPointsParallel
        (
            const LongList<DynList<label, 4> >& faceCandidates,
            const LongList<labelPair>& candidatePatches
        );

        //- create layer cells in one go
        //- this is much faster than layer by layer
        void createLayerCells(const labelList& patchLabels);

        //- helper function finding a new face label for multiply extruded nodes
        inline label findNewNodeLabel
        (
            const label pointI,
            const label pKey
        ) const;

        //- creating hex cells near feature edges
        inline void createNewCellFromEdge
        (
            const edge& e,
            const label pKeyI,
            const label pKeyJ,
            FixedList<FixedList<label, 4>, 6>& cellFaces
        ) const;

        //- creating hex cells near corners
        inline void createNewCellFromNode
        (
            const label pointI,
            const DynList<label, 3>& pKeys,
            FixedList<FixedList<label, 4>, 6>& cellFaces
        ) const;

        //- create a bnd layer for a given patch
        void addLayerForPatch(const label patchLabel);

        //- delete meshSurfaceEngine
        inline void clearOut()
        {
            deleteDemandDrivenData(msePtr_);
            deleteDemandDrivenData(meshPartitionerPtr_);
        }

    // Enumerators
        enum vertexTypes
        {
            NONE = 0,
            PATCHNODE = 1,
            EDGENODE = 2,
            CORNERNODE = 4,
            PARALLELBOUNDARY = 8
        };

        //- Disallow bitwise copy construct
        boundaryLayers(const boundaryLayers&);

        //- Disallow bitwise assignment
        void operator=(const boundaryLayers&);

public:

    // Constructors

        //- Construct from mesh reference
        boundaryLayers(polyMeshGen& mesh, const dictionary& meshDict);

    // Destructor
        ~boundaryLayers();

        //- return faces with mixed zero/nonzero layer displacement
        const labelHashSet& mixedDisplacementFaces() const
        {
            return mixedDisplacementFaces_;
        }

        //- return boundary points on sharp BL/BL junction curves
        const labelHashSet& junctionEdgePoints() const
        {
            return blblJunctionPoints_;
        }

        const labelHashSet& blblCornerPoints() const
        {
            return blblCornerPoints_;
        }

        const labelHashSet& blblAcuteCornerPoints() const
        {
            return blblAcuteCornerPoints_;
        }

        //- Return per-boundary-point intended layer scale
        //- (0=suppressed, 1=full BL). Use layerScale < 0.01 as suppressed
        //- threshold. Reports effective BL coverage, not vertex existence.
        const scalarField& layerScale() const
        {
            return layerScale_;
        }

        //- Return virtual topology face ring assignments
        const labelList& vtFaceRing() const
        {
            return vtFaceRing_;
        }

        //- return edges where BL and no-BL patches meet
        const labelHashSet& blNoBlEdges() const
        {
            return blNoBlEdges_;
        }

        //- detect BL/no-BL transition edges and populate blNoBlEdgePoints_
        //- lightweight: reads meshDict + mesh surface only, no mesh modification
        //- safe to call before addLayerForAllPatches()
        void detectBLNoBlTransitionEdges() const;

        //- return boundary points on BL/no-BL transition edges
        const labelHashSet& blNoBlEdgePoints() const
        {
            return blNoBlEdgePoints_;
        }

        //- return patch index map for BL/no-BL interface points
        const Map<label>& blNoBlPointPatch() const
        {
            return blNoBlPointPatch_;
        }

        //- return BL/neutral edge points (blade/periodic junctions)
        const labelHashSet& blNeutralEdgePoints() const
        {
            return blNeutralEdgePoints_;
        }

        //- return BL-side patch per neutral edge point
        const Map<label>& blNeutralPointPatch() const
        {
            return blNeutralPointPatch_;
        }

        //- set gap points for BL suppression in thin clearance regions
        void setGapPoints(const labelHashSet& gpts)
        {
            gapPoints_ = gpts;
        }
        void addGapPoints(const labelHashSet& gpts)
        {
            forAllConstIter(labelHashSet, gpts, iter)
                gapPoints_.insert(iter.key());
        }
        void setTripleJunctionPoints(const labelHashSet& pts)
        {
            tripleJunctionPoints_ = pts;
        }
        void addTripleJunctionPoints(const labelHashSet& pts)
        {
            forAllConstIter(labelHashSet, pts, iter)
                tripleJunctionPoints_.insert(iter.key());
        }

    // Public member functions
        //- adds layer for a given patch
        void addLayerForPatch(const word& patchName);

        //- create O-topology layers (used as flag)
        void createOTopologyLayers();

        //- terminate boundary layers at concave edges (used as a flag)
        void terminateLayersAtConcaveEdges();
        void markConcaveEdgePoints(boolList& skipPoint) const;
        //- Apply virtual topology BL exclusion zones around acute
        //- BL triple-junctions before any layer graph is constructed.
        void applyVirtualTopologyExclusion() const;
        //- Apply face-ring BL suppression/taper around gap points
        void applyGapFaceRingExclusion() const;
        void reportBLTransitionSeeds() const;
        void buildBLTransitionPlan() const;

        //- avoid generating layers for empty patches in case of a 2D mesh
        //- used as a flag prior to addLayerForAllPatches
        void activate2DMode();

        //- add layers for all patches
        void addLayerForAllPatches();

        //- add wrapper layer
        //- this function is intended for usage before surface recovery
        void addWrapperLayer();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "boundaryLayersI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
