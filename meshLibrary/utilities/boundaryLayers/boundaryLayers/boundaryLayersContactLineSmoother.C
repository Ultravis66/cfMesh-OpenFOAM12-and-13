/*---------------------------------------------------------------------------*\
  Contact-line height limiter for boundary layer vertices.

  Runs after createNewVertices() coordinate swap, before rollback/repair.
  Detects local height spikes and collapses along contact-line points and
  optionally corrects them.

  Fixes vs v1:
  - contactMask prevents averaging across unrelated contact families
  - zeroDistPoints_ anchor prevents lifting truly collapsed singular points
  - Method is non-const (mutates mesh_.points() directly)
  - Atlas write gated by writeContactLineHeightSmootherAtlas_
\*---------------------------------------------------------------------------*/

#include "boundaryLayers.H"
#include "meshSurfaceEngine.H"
#include "OFstream.H"
#include "helperFunctions.H"

namespace Foam
{

void boundaryLayers::smoothContactLineHeights
(
    const word& passName,
    const labelList& bPoints,
    const VRWGraph& pointPoints
)
{
    // Do NOT call surfaceEngine() here -- it would rebuild against
    // swapped coordinates and corrupt mesh state.
    // bPoints and pointPoints are passed from the call site.
    pointFieldPMG& points = mesh_.points();

    // Build candidate set with contact-family bitmask
    // bit 1 = BL/no-BL (inlet/outlet termination)
    // bit 2 = BL/neutral (periodic seam)
    // bit 4 = BL/BL junction (hub/blade, shroud/blade)
    // bit 8 = BL/BL corner
    labelHashSet candidateSet;
    Map<label> contactMask;

    auto insertMasked = [&](const labelHashSet& src, label bit)
    {
        forAllConstIter(labelHashSet, src, it)
        {
            const label bpI = it.key();
            candidateSet.insert(bpI);
            if( contactMask.found(bpI) )
                contactMask[bpI] |= bit;
            else
                contactMask.insert(bpI, bit);
        }
    };

    insertMasked(blNoBlEdgePoints_,    1);
    insertMasked(blNeutralEdgePoints_, 2);
    insertMasked(blblJunctionPoints_,  4);
    insertMasked(blblCornerPoints_,    8);

    const label nCandidates = candidateSet.size();
    Info << "Contact-line height smoother (" << passName << "): "
         << nCandidates << " candidate points" << endl;

    if( nCandidates == 0 ) return;

    // Collect heights from current point positions
    Map<scalar> heightMap;
    Map<vector> dirMap;
    Map<label>  basePtMap;

    forAllConstIter(labelHashSet, candidateSet, it)
    {
        const label bpI = it.key();
        if( bpI < 0 || bpI >= label(bPoints.size()) ) continue;
        const label meshPtI = bPoints[bpI];
        const label basePtI = newLabelForVertex_[meshPtI];
        if( basePtI < 0 || basePtI >= label(points.size()) ) continue;
        if( meshPtI < 0 || meshPtI >= label(points.size()) ) continue;

        const point base  = points[basePtI];
        const point layer = points[meshPtI];
        const vector hVec = layer - base;
        const scalar h    = mag(hVec);

        heightMap.insert(bpI, h);
        dirMap.insert(bpI, h > VSMALL ? hVec/h : vector::zero);
        basePtMap.insert(bpI, basePtI);
    }

    // Open atlas CSV if requested
    autoPtr<OFstream> osPtr;
    if( writeContactLineHeightSmootherAtlas_ )
    {
        osPtr.reset(new OFstream(
            "blContactLineHeightSmootherAtlas_" + passName + ".csv"));
        osPtr() << "bpI,meshPtI,basePtI,contactMask,cls,"
                << "x,y,z,h,nNei,hMin,hAvg,hMax,ratioToAvg,"
                << "isZeroSeed,wouldSpikeClamp,wouldPinchClamp,hTarget" << nl;
    }

    label nSpike = 0, nPinch = 0, nMoved = 0, nAnchored = 0;

    forAllConstIter(labelHashSet, candidateSet, it)
    {
        const label bpI = it.key();
        if( !heightMap.found(bpI) ) continue;

        const label meshPtI = bPoints[bpI];
        const label basePtI = basePtMap[bpI];
        const scalar h = heightMap[bpI];
        const point& base = points[basePtI];

        const label myMask = contactMask.found(bpI) ? contactMask[bpI] : 0;

        // Collect neighbor heights -- same contact family only
        scalar hSum = 0; scalar hMin = GREAT; scalar hMax = -GREAT;
        label nNei = 0;
        forAllRow(pointPoints, bpI, ppI)
        {
            const label nbpI = pointPoints(bpI, ppI);
            if( !heightMap.found(nbpI) ) continue;
            if( !contactMask.found(nbpI) ) continue;
            // Must share at least one contact family bit
            if( (myMask & contactMask[nbpI]) == 0 ) continue;
            const scalar hn = heightMap[nbpI];
            hSum += hn;
            hMin = Foam::min(hMin, hn);
            hMax = Foam::max(hMax, hn);
            ++nNei;
        }

        if( nNei < 2 ) continue;

        const scalar hAvg = hSum / scalar(nNei);
        const scalar ratioToAvg = (hAvg > VSMALL) ? h/hAvg : scalar(1);

        // Determine atlas class
        label cls = 0;
        Map<label>::const_iterator clsIt = blContactPointClass_.find(bpI);
        if( clsIt != blContactPointClass_.end() ) cls = clsIt();

        // Zero-seed anchor: don't lift truly collapsed singular points
        // rampSeedPoints_ have zeroDistPoints_=true for topology contract
        // but are finite-height -- they CAN be lifted by pinch correction.
        const bool isRampSeed =
            (bpI >= 0
          && bpI < label(rampSeedPoints_.size())
          && rampSeedPoints_[bpI]);
        const bool isZeroSeed =
            (bpI < label(zeroDistPoints_.size()) && zeroDistPoints_[bpI])
            && !isRampSeed;

        // Compute target
        scalar hTarget = h;
        bool wouldSpike = false;
        bool wouldPinch = false;

        if( hAvg > VSMALL )
        {
            if( h > contactLineSmootherMaxGrowth_ * hAvg )
            {
                hTarget = contactLineSmootherMaxGrowth_ * hAvg;
                wouldSpike = true;
                ++nSpike;
            }
            else if( h < contactLineSmootherMinCollapse_ * hAvg )
            {
                hTarget = contactLineSmootherMinCollapse_ * hAvg;
                wouldPinch = true;
                ++nPinch;
            }
        }
        hTarget = Foam::max(hTarget, contactLineSmootherMinHeight_);

        // Write CSV row
        if( writeContactLineHeightSmootherAtlas_ && osPtr.valid() )
        {
            osPtr() << bpI << "," << meshPtI << "," << basePtI << ","
                    << myMask << "," << cls << ","
                    << base.x() << "," << base.y() << "," << base.z() << ","
                    << h << "," << nNei << ","
                    << hMin << "," << hAvg << "," << hMax << ","
                    << ratioToAvg << ","
                    << (isZeroSeed ? "1" : "0") << ","
                    << (wouldSpike ? "1" : "0") << ","
                    << (wouldPinch ? "1" : "0") << ","
                    << hTarget << nl;
        }

        // Apply if behavior mode enabled
        if( !useContactLineHeightSmoother_ ) continue;

        // Anchor: never move TripleJunction or zero-seed points
        if( cls == 1 || isZeroSeed )
        {
            ++nAnchored;
            continue;
        }

        if( h <= VSMALL ) continue;
        if( !wouldSpike && !wouldPinch ) continue;

        // Damped move toward target along born extrusion direction
        const scalar limitedH =
            h + contactLineSmootherMaxMoveFraction_ * (hTarget - h);

        const vector& dir = dirMap[bpI];
        const point newLayer = base + limitedH * dir;

        if( !help::isnan(newLayer) && !help::isinf(newLayer) )
        {
            points[meshPtI] = newLayer;
            ++nMoved;
        }
    }

    Info << "Contact-line height smoother (" << passName << "): "
         << "spikes=" << nSpike
         << " pinches=" << nPinch
         << " anchored=" << nAnchored
         << " moved=" << nMoved
         << " (smoother " << (useContactLineHeightSmoother_ ? "ON" : "OFF")
         << ")" << endl;
}

} // End namespace Foam
