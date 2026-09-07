/*---------------------------------------------------------------------------*\
  CFMitch - cfMesh Boundary-Layer Constraint Planner
\*---------------------------------------------------------------------------*/

#ifndef boundaryLayerConstraintPlanner_H
#define boundaryLayerConstraintPlanner_H

#include "polyMeshGen.H"
#include "meshSurfaceEngine.H"
#include "boundaryLayerPlan.H"
#include "LongList.H"
#include "Map.H"
#include "HashSet.H"
#include "edge.H"
#include "word.H"

#include <map>

namespace Foam
{

class boundaryLayerConstraintPlanner
{
public:

    enum architectureType
    {
        CLASSIC_CFMESH,
        LEGACY_ENHANCED,
        CONSTRAINT_PLANNER
    };


private:

    polyMeshGen& mesh_;

    architectureType architecture_;


    static architectureType parseArchitecture
    (
        const word& architectureName
    );


public:

    boundaryLayerConstraintPlanner
    (
        polyMeshGen& mesh,
        const word& architectureName
    );

    ~boundaryLayerConstraintPlanner() = default;


    architectureType architecture() const
    {
        return architecture_;
    }


    bool plannerEnabled() const
    {
        return architecture_ == CONSTRAINT_PLANNER;
    }


    static word architectureName
    (
        const architectureType architecture
    );


    void report() const;


    //- Resolve height-compatible boundary-face layer counts and then
    //- enforce same-patch layer-count compatibility.
    //
    //- Stable constraints are keyed by persistent mesh-point labels.
    //- Current face indices are deliberately rebuilt from mse.
    boundaryLayerPlan solveLayerCounts
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
    ) const;
};

} // End namespace Foam

#endif
