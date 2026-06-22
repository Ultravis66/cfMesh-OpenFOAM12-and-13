/*---------------------------------------------------------------------------*\
    cfMesh extension: mesh quality / stage-outcome types
    Provides MeshQuality snapshot, StageOutcome enum, and MeshHistory enum.
    These replace ad-hoc booleans and loose strings for repair-stage state.
\*---------------------------------------------------------------------------*/
#ifndef meshQualityTypes_H
#define meshQualityTypes_H

#include "label.H"
#include "scalar.H"

namespace Foam
{

// ---------------------------------------------------------------------------
// Scalar quality snapshot at any pipeline stage
// ---------------------------------------------------------------------------
struct MeshQuality
{
    label negVol;       // cells with negative volume
    label badPyramids;  // faces failing pyramid orientation test
    label openCells;    // non-closed cells
    scalar maxSkew;     // maximum face skewness

    MeshQuality()
    :
        negVol(0),
        badPyramids(0),
        openCells(0),
        maxSkew(0.0)
    {}

    MeshQuality(label nv, label bp, label oc, scalar sk)
    :
        negVol(nv),
        badPyramids(bp),
        openCells(oc),
        maxSkew(sk)
    {}

    // True if this quality is no worse than 'ref' on the critical metrics.
    // negVol must not increase; badPyramids may increase by allowedPyrIncrease.
    bool noWorseThan
    (
        const MeshQuality& ref,
        const label allowedPyrIncrease = 0
    ) const
    {
        return negVol <= ref.negVol
            && badPyramids <= ref.badPyramids + allowedPyrIncrease
            && openCells  <= ref.openCells;
    }

    bool isClean() const { return negVol == 0; }
};


// ---------------------------------------------------------------------------
// Outcome of a single repair or optimization stage
// ---------------------------------------------------------------------------
enum class StageOutcome
{
    NotRun,                  // stage was gated off / skipped
    Accepted,                // ran and accepted (quality improved or held)
    RejectedWorse,           // ran, quality worsened, rolled back
    RejectedInvalidMetadata, // ran, metadata inconsistent, rolled back
    SkippedUnsafe            // skipped because upstream state was unsafe
};

inline const char* stageOutcomeName(StageOutcome o)
{
    switch(o)
    {
        case StageOutcome::NotRun:                  return "NotRun";
        case StageOutcome::Accepted:                return "Accepted";
        case StageOutcome::RejectedWorse:           return "RejectedWorse";
        case StageOutcome::RejectedInvalidMetadata: return "RejectedInvalidMetadata";
        case StageOutcome::SkippedUnsafe:           return "SkippedUnsafe";
        default:                                    return "Unknown";
    }
}


// ---------------------------------------------------------------------------
// Coarse history label for the whole mesh at pipeline exit
// ---------------------------------------------------------------------------
enum class MeshHistory
{
    CleanNatural,     // optimizer reached negVol=0 without dirty cleanup
    CleanPromoted,    // dirty cleanup cleared negVol; topologically clean
                      //   but junction neighborhood may be wounded
    DirtyRecoverable, // residual negVol <= threshold; cleanup failed/skipped
    DirtyUnsafe,      // residual negVol > threshold; hard unsafe
    HardRejected      // finalUntangleRejected_: BL skipped entirely
};

inline const char* meshHistoryName(MeshHistory h)
{
    switch(h)
    {
        case MeshHistory::CleanNatural:     return "CleanNatural";
        case MeshHistory::CleanPromoted:    return "CleanPromoted";
        case MeshHistory::DirtyRecoverable: return "DirtyRecoverable";
        case MeshHistory::DirtyUnsafe:      return "DirtyUnsafe";
        case MeshHistory::HardRejected:     return "HardRejected";
        default:                            return "Unknown";
    }
}

} // End namespace Foam

#endif
