/*---------------------------------------------------------------------------*\
Description
    BLRepairPlan: transactional boundary-layer junction repair descriptor.

    Encapsulates the complete repair strategy for a set of defective
    boundary faces: provenance seed faces, per-ring layer caps, per-ring
    thickness scales, junction type classification, and corner locking.

    Junction types are role-based (from PatchRoleMap), never patch-name
    specific. The same C++ logic works for any geometry.
\*---------------------------------------------------------------------------*/

#ifndef BLRepairPlan_H
#define BLRepairPlan_H

#include "HashSet.H"

#include "scalar.H"
#include "word.H"
#include "Map.H"
#include "IOstreams.H"

namespace Foam
{

class BLRepairPlan
{
public:

    // Junction type -- role-based, never patch-name specific
    enum junctionType
    {
        WALL_PERIODIC,     // wall + periodic
        WALL_WALL,         // two wall patches meeting
        HUB_FLOW,          // hub + inlet or outlet
        SHROUD_FLOW,       // shroud + inlet or outlet
        GENERIC_TRIPLE,    // unclassified triple junction
        BASE_MESH_DEFECT,  // no BL provenance
        UNKNOWN
    };

    // Per-ring repair parameters
    struct ringParams
    {
        label  maxLayers;       //- 0 = uncapped
        scalar thicknessScale;  //- 1.0 = full thickness

        ringParams() : maxLayers(0), thicknessScale(1.0) {}

        ringParams(const label ml, const scalar ts)
        : maxLayers(ml), thicknessScale(ts) {}
    };

    // Data
    HashSet<label> seedBfI_;
    ringParams   ring0_;
    ringParams   ring1_;
    ringParams   ring2_;
    junctionType sourceType_;
    word         sourceReason_;
    bool         lockCornerPoints_;

    // Constructors
    BLRepairPlan()
    :
        seedBfI_(),
        ring0_(), ring1_(), ring2_(),
        sourceType_(UNKNOWN),
        sourceReason_("unknown"),
        lockCornerPoints_(false)
    {}

    BLRepairPlan
    (
        const HashSet<label>& seeds,
        const ringParams& r0,
        const ringParams& r1,
        const ringParams& r2,
        const junctionType jt = UNKNOWN,
        const word& reason = "manual",
        const bool lockCorners = false
    )
    :
        seedBfI_(seeds),
        ring0_(r0), ring1_(r1), ring2_(r2),
        sourceType_(jt),
        sourceReason_(reason),
        lockCornerPoints_(lockCorners)
    {}

    // Member functions
    bool active() const
    {
        return seedBfI_.size() > 0
            && (ring0_.maxLayers > 0 || ring1_.maxLayers > 0);
    }

    bool hasThicknessScaling() const
    {
        return ring0_.thicknessScale < 1.0 - SMALL
            || ring1_.thicknessScale < 1.0 - SMALL
            || ring2_.thicknessScale < 1.0 - SMALL;
    }

    static word junctionTypeName(const junctionType jt)
    {
        switch(jt)
        {
            case WALL_PERIODIC:    return "wall-periodic";
            case WALL_WALL:        return "wall-wall";
            case HUB_FLOW:         return "hub-flow";
            case SHROUD_FLOW:      return "shroud-flow";
            case GENERIC_TRIPLE:   return "generic-triple";
            case BASE_MESH_DEFECT: return "base-mesh-defect";
            default:               return "unknown";
        }
    }

    void print() const
    {
        Info << "BLRepairPlan:"
             << " seeds=" << seedBfI_.size()
             << " type=" << junctionTypeName(sourceType_)
             << " reason=" << sourceReason_
             << " ring0(maxL=" << ring0_.maxLayers
             << " thk=" << ring0_.thicknessScale << ")"
             << " ring1(maxL=" << ring1_.maxLayers
             << " thk=" << ring1_.thicknessScale << ")"
             << " ring2(maxL=" << ring2_.maxLayers
             << " thk=" << ring2_.thicknessScale << ")"
             << " lockCorners=" << lockCornerPoints_
             << endl;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
