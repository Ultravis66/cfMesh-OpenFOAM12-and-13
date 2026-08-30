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
    meshSurfaceMapper

Description
    Maps vertices to the nearest point on the geometry surface

SourceFiles
    meshSurfaceMapper.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceMapper_H
#define meshSurfaceMapper_H

#include "boolList.H"

#include "labelList.H"
#include "pointField.H"
#include "labelLongList.H"
#include "HashSet.H"
#include "Map.H"
#include "parMapperHelper.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class meshSurfaceEngine;
class meshSurfacePartitioner;
class triSurfacePartitioner;
class sourceFeatureGraph;

/*---------------------------------------------------------------------------*\
                    Class meshSurfaceMapper Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceMapper
{
    // Private data
        //- mesh surface
        const meshSurfaceEngine& surfaceEngine_;

        //- reference to the octree
        const meshOctree& meshOctree_;

        //- relaxation factor for corner snapping (0=no move, 1=full)
        scalar cornerSnapRelaxation_ = 0.25;

        // ------------------------------------------------------------
        // v7 (2026-08-09): narrow blade+periodic+hub/shroud triple-
        // junction fix, built after v3-v6 established that (a) a real
        // findNearestCorner()-based fix mechanism works locally but
        // (b) applying it to ANY 3+ patch corner (v5/v6) destabilizes
        // the mesh, with damage traced to the surface-mapping stage
        // BEFORE BL creation or the optimizer even run. SOL hypothesis:
        // a triple junction is not a free point but the anchor of a
        // constrained local feature-edge network; moving the corner
        // without its incident edges creates a topological kink that
        // mapEdgeNodes()'s independent accept/reject logic amplifies.
        // v7 narrows the fix to EXACTLY blade+periodic+(hub OR shroud)
        // corners (no 4+ patch, no inlet/outlet) and instruments what
        // mapEdgeNodes() actually does to the corner's incident edge
        // points afterward, using only production-computed values (no
        // new octree queries, avoiding the v3 measurement-perturbation
        // class of bug).
        // ------------------------------------------------------------

        //- EXPLICIT fix-enable gate, separate from any diagnostic tag.
        //- Defaults false (v6 lesson: an unenabled call site must run
        //- ORIGINAL frozen behavior, never the fix, by default).
        bool tripleJunctionCornerFixEnabled_ = false;

        //- tag for log-line labeling only, does not gate behavior.
        word tripleJunctionDiagTag_;

        //- patch-ID sets defining the narrow target junction type.
        //- Populated via setTripleJunctionTargetPatches(). Empty by
        //- default -- an empty set means "no patches of this role",
        //- so an unpopulated mapper safely classifies zero corners as
        //- target junctions (falls through to original frozen behavior
        //- for everything, same fail-safe principle as the main gate).
        labelHashSet tjBladePatchIds_;
        labelHashSet tjPeriodicPatchIds_;
        labelHashSet tjHubPatchIds_;
        labelHashSet tjShroudPatchIds_;
        labelHashSet tjInletOutletPatchIds_;

        //- key/value label space made explicit after an
        //- earlier ambiguous comment ("bPoints[bpI]-space") was flagged
        //- as a correctness risk -- confirmed against actual call sites
        //- in meshSurfaceMapperCornersAndEdges.C, both key and value are
        //- BOUNDARY-POINT indices (bpI, into surfaceEngine_.boundaryPoints()),
        //- NOT mesh-point labels (bPoints[bpI]).
        //-   key:   boundary-point index (bpI) of an immediate neighbor
        //-          of an accepted target-corner move
        //-   value: boundary-point index (bpI) of the source triple
        //-          corner that move originated from
        //- Populated by the fix pass in THIS mapCorners() call, consumed
        //- by mapEdgeNodes() in the SAME mapCornersAndEdges() sequence
        //- to instrument what production edge logic does to them --
        //- observe only, no behavior change. Cleared at the start of
        //- each mapCorners() call, and again after mapEdgeNodes()
        //- consumes it (one-shot, per SOL review) to avoid stale
        //- diagnostic attribution if these functions are ever called
        //- out of the expected corners-then-edges sequence.
        Map<label> tjWatchedBoundaryPointToSourceCorner_;

        //- v7: makes the watcher state one-shot. Only true between a
        //- mapCorners() call that accepted >=1 target-corner move and
        //- the immediately following mapEdgeNodes() call consuming it.
        bool tjEdgeWatchArmed_ = false;

        //- Code review finding: probeTriSurfJunction() (a heavy,
        //- brute-force forensic scan over every surface point, built
        //- for one-time investigation) was unconditionally embedded in
        //- the v7 fix pass, meaning re-enabling corner-snapping via
        //- tripleJunctionCornerFixEnabled_ would ALSO silently re-run
        //- the full diagnostic scan on every target corner, every
        //- run. Decoupled: this flag now independently gates the
        //- expensive diagnostic probe. Defaults false -- corner
        //- snapping can be enabled without the diagnostic overhead.
        bool tripleJunctionDeepDiagnosticEnabled_ = false;

        //- mesh surface partitioner
        mutable const meshSurfacePartitioner* surfaceEnginePartitionerPtr_;
        const bool deletePartitioner_;

        //- triSurface partitioner
        mutable triSurfacePartitioner* surfPartitionerPtr_;

        //- Unified source-geometry feature representation.
        //- Lazily constructed from meshOctree_.surface().
        //- Infrastructure only; no mapping behaviour change.
        mutable sourceFeatureGraph* sourceFeatureGraphPtr_;

        //- minimum pyramid height for validity check (default 1e-10)
        scalar snapMinPyramidHeight_ = 1e-10;

        //- boundary-point indices excluded from generic surface projection
        //- set via setProtectedPoints(); empty by default (no exclusion)
        labelHashSet protectedPoints_;
        labelHashSet blNeutralPoints_;
        Map<label> blNeutralPointPatches_;

        //- Boundary points rejected by validity check in last
        //- mapVerticesOntoSurface call -- used by repairRejectedPoints
        mutable DynamicList<label> rejectedBpI_;

        //- BL-side patch index per protected boundary point
        //- used for patch-constrained projection instead of global search
        //- set via setProtectedPointPatches(); empty by default
        Map<label> protectedPointPatches_;

    // Private member functions
        //- create and return mesh surface partitioner
        void createMeshSurfacePartitioner() const;
        inline const meshSurfacePartitioner& meshPartitioner() const
        {
            if( !surfaceEnginePartitionerPtr_ )
                createMeshSurfacePartitioner();

            return *surfaceEnginePartitionerPtr_;
        }

        //- create and return triSurfacePartitioner
        void createTriSurfacePartitioner() const;
        inline const triSurfacePartitioner& surfacePartitioner() const
        {
            if( !surfPartitionerPtr_ )
                createTriSurfacePartitioner();

            return *surfPartitionerPtr_;
        }

        //- create source feature graph
        void createSourceFeatureGraph() const;

        //- return lazily constructed source feature graph
        const sourceFeatureGraph& sourceFeatures() const;

        //- delete demand-driven mapper addressing
        void clearOut();

    // Private member functions
        //- map corner nodes to the boundary


        //- find mapping distance for selected points
        void findMappingDistance
        (
            const labelLongList& nodesToMap,
            scalarList & mappingDistance
        ) const;

        //- calculate face metric
        scalar faceMetricInPatch(const label bfI, const label patchI) const;

        //- Quality evaluator: returns true if moving boundary point bpI
        //- to candidate leaves all incident faces valid.
        //- Checks: pyramid height, face area, movement cap.
        //- Used by smoothing, bisection, and corner snap.
        bool proposedMoveIsValid
        (
            const label bpI,
            const point& candidate,
            const point& oldPosition,
            const scalar mappingDistSq
        ) const;

    // Private member functions needed for parallel execution

        //- check if nodes at parallel boundaries are selected at all processors
        void selectNodesAtParallelBnd(const labelLongList&);

        //- map to the smallest distance
        //- makes sense for parallel calculations, only
        void mapToSmallestDistance(LongList<parMapperHelper>&);

        //- Disallow default bitwise copy construct
        meshSurfaceMapper(const meshSurfaceMapper&);

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceMapper&);

public:

    // Constructors

        //- Construct from meshSurfaceEngine and octree
        meshSurfaceMapper(const meshSurfaceEngine&, const meshOctree&);

        //- Construct from meshSurfacePartitioner and octree
        meshSurfaceMapper(const meshSurfacePartitioner&, const meshOctree&);

    // Destructor

        ~meshSurfaceMapper();

    // Member Functions
        //- map the given surface node on the selected patch
        //- this does not make sense for vertices at parallel boundaries
        void mapNodeToPatch(const label bpI, const label patchI = -1);

        //- projects surface vertices onto their nearest location
        //- on the surface mesh
        void mapVerticesOntoSurface();

        //- projects selected surface vertices to their nearest location
        //- on the surface mesh
        void mapVerticesOntoSurface(const labelLongList& nodesToMap);

        //- projects corner and edge vertices onto their nearest location
        //- on the surface mesh
        void mapCornersAndEdges();

        //- set corner snap relaxation factor (0=no move, 1=full snap)
        void setCornerSnapRelaxation(const scalar relax)
        {
            cornerSnapRelaxation_ = relax;
        }

        //- v7: EXPLICIT enable/disable for the narrow triple-junction
        //- fix. Defaults false. Must be true only immediately around
        //- the one intended call site.
        void setTripleJunctionCornerFix
        (
            const bool enabled,
            const word& tag = word("default")
        )
        {
            tripleJunctionCornerFixEnabled_ = enabled;
            tripleJunctionDiagTag_ = tag;
        }

        //- Code review addition: independently gates the expensive
        //- probeTriSurfJunction() forensic scan. Separate from the
        //- corner-fix-enable flag above -- see member comment.
        void setTripleJunctionDeepDiagnostic(const bool enabled)
        {
            tripleJunctionDeepDiagnosticEnabled_ = enabled;
        }

        //- v7: define the narrow target junction type via explicit
        //- patch-ID sets (resolved by the CALLER, e.g. via PatchRoleMap
        //- in cartesianMeshGenerator -- NOT by threading PatchRoleMap
        //- through this class, per SOL review: too much architecture
        //- for this experiment). A corner is a "target" triple junction
        //- iff it touches EXACTLY 3 patches: one from bladeIds, one from
        //- periodicIds, one from (hubIds union shroudIds), and none
        //- from inletOutletIds. Any set may be left empty if that role
        //- is irrelevant to the current geometry.
        void setTripleJunctionTargetPatches
        (
            const labelHashSet& bladeIds,
            const labelHashSet& periodicIds,
            const labelHashSet& hubIds,
            const labelHashSet& shroudIds,
            const labelHashSet& inletOutletIds
        )
        {
            // Fail-closed disjointness check. A role-
            // classification bug placing one patch ID in two role sets
            // would otherwise silently misclassify junctions downstream
            // -- this project has already hit at least one real role-
            // ordering/classification bug this session, so this is not
            // a hypothetical risk.
            const List<const labelHashSet*> roleSets
            {
                &bladeIds, &periodicIds, &hubIds, &shroudIds, &inletOutletIds
            };
            for( label a = 0; a < roleSets.size(); ++a )
            {
                for( label b = a+1; b < roleSets.size(); ++b )
                {
                    forAllConstIter(labelHashSet, *roleSets[a], it)
                    {
                        if( roleSets[b]->found(it.key()) )
                        {
                            FatalErrorInFunction
                                << "setTripleJunctionTargetPatches: patch ID "
                                << it.key() << " appears in more than one"
                                << " role set -- role sets must be disjoint."
                                << exit(FatalError);
                        }
                    }
                }
            }

            tjBladePatchIds_ = bladeIds;
            tjPeriodicPatchIds_ = periodicIds;
            tjHubPatchIds_ = hubIds;
            tjShroudPatchIds_ = shroudIds;
            tjInletOutletPatchIds_ = inletOutletIds;
        }

        //- set boundary-point indices that must not be moved
        //- by generic nearest-surface projection.
        //- label space: boundary-point indices (into surfaceEngine_.boundaryPoints())
        //- defaults to empty -- no effect unless explicitly set
        void setProtectedPoints(const labelHashSet& pts)
        {
            protectedPoints_ = pts;
        }

        //- set BL/neutral edge points (blade/periodic junctions)
        void setBLNeutralPoints(const labelHashSet& pts)
        {
            blNeutralPoints_ = pts;
        }

        //- set BL-side patch per neutral edge point
        void setBLNeutralPointPatches(const Map<label>& patches)
        {
            blNeutralPointPatches_ = patches;
        }

        //- set patch index map for patch-constrained projection
        //- of BL/no-BL interface points
        void setProtectedPointPatches(const Map<label>& patches)
        {
            protectedPointPatches_ = patches;
        }

        //- projects corner vertices onto the nearest corner on the surface
        void mapCorners(const labelLongList& nodesToMap);

        //- projects selected edge vertices onto their nearest
        //- locations on the surface mesh
        void mapEdgeNodes(const labelLongList& nodesToMap);

        //- Same feature-edge mapping, additionally reporting whether each
        //- requested point completed a valid constrained mapping.
        //- mappingAccepted is indexed like nodesToMap.
        void mapEdgeNodes
        (
            const labelLongList& nodesToMap,
            boolList& mappingAccepted
        );

        //- Map edge nodes using desiredPositions only as the geometric
        //- smoothing/search intent. The live mesh point is not moved before
        //- feature resolution. Final feature targets are committed against
        //- the current live volume state.
        //- desiredPositions is indexed like nodesToMap.
        void mapEdgeNodes
        (
            const labelLongList& nodesToMap,
            const pointField& desiredPositions,
            boolList& mappingAccepted
        );

        //- projects surface vertices onto the surface with respect
        //- to the surface patch they belong to. Edges and corner are respected
        void mapVerticesOntoSurfacePatches();

        //- projects selected surface vertices onto the surface with respect
        //- to the surface patch they belong to. Edges and corner are respected
        void mapVerticesOntoSurfacePatches(const labelLongList& nodesToMap);

        //- a combination of mapping and smoothing intended for better
        //- feature capturing
        void preMapVertices(const label nIterations = 3);

        //- Constrained surface smoothing: redistributes single-patch
        //- boundary points using patch-aware Laplacian smoothing.
        //- Feature-edge, corner, BL/no-BL, and non-manifold points
        //- are locked. Moves accepted only if proposedMoveIsValid.
        void smoothSinglePatchPoints(const label nIterations = 3);

        //- Targeted repair of points rejected by mapVerticesOntoSurface.
        //- Smooths rejected points + 2-ring neighbours, constrained by
        //- proposedMoveIsValid. Single-patch: Laplacian + re-project.
        //- Two-patch-edge: project to nearest feature curve point.
        //- Corners/BL/nonmanifold: locked.
        void repairRejectedPoints();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
