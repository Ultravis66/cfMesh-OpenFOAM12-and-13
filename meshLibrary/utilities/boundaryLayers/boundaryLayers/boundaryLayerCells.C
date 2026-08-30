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

#include "boundaryLayers.H"
#include "OFstream.H"
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "helperFunctionsPar.H"
#include "demandDrivenData.H"
#include "VRWGraphList.H"

#include "labelledPair.H"
#include "HashSet.H"

#include <map>
#include <unordered_map>
#include <set>

//#define DEBUGLayer

# ifdef DEBUGLayer
#include "polyMeshGenAddressing.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayers::createLayerCells(const labelList& patchLabels)
{
    Info << "Starting creating layer cells" << endl;

    const meshSurfaceEngine& mse = surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const edgeList& edges = mse.edges();
    const VRWGraph& faceEdges = mse.faceEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();
    const labelList& boundaryFacePatches = mse.boundaryFacePatches();
    const labelList& faceOwners = mse.faceOwners();
    const labelList& bp = mse.bp();
    const VRWGraph& pointFaces = mse.pointFaces();
    const pointFieldPMG& points = mesh_.points();

    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const VRWGraph& pointPatches = mPart.pointPatches();

    //- mark patches which will be extruded into layer cells
    boolList treatPatches(mesh_.boundaries().size(), false);
    forAll(patchLabels, patchI)
    {
        const label pLabel = patchLabels[patchI];
        forAll(treatPatchesWithPatch_[pLabel], i)
            treatPatches[treatPatchesWithPatch_[pLabel][i]] = true;
    }

    //- create new faces at parallel boundaries
    const Map<label>* otherProcPatchPtr(NULL);
    const Map<label>* otherFaceProcPtr(NULL);

    if( Pstream::parRun() )
    {
        createNewFacesParallel(treatPatches);

        otherProcPatchPtr = &mse.otherEdgeFacePatch();
        otherFaceProcPtr = &mse.otherEdgeFaceAtProc();
    }

    //- create lists for new boundary faces
    VRWGraph newBoundaryFaces;

    // Global top-vertex collapse rules discovered while building
    // canonical triangle transition cells. If one normal cell
    // collapses top T -> base B, every queued cell/boundary face
    // touching T must use B before addCells().
    std::map<label,label> globalNormalTopSubst;
    // Demand collector: for non-base (shared) substitution candidates,
    // track EVERY distinct target requested for a given raw vertex
    // across the whole construction pass. A raw vertex is only safe
    // to substitute globally if ALL demands agree on one target.
    // Conflicting demands must fall back to base-collapse instead.
    std::map<label, std::set<label>> normalTopSubstDemands;
    labelLongList newBoundaryOwners;
    labelLongList newBoundaryPatches;

    //- create storage for new cells
    VRWGraphList cellsToAdd;

    //- create layer cells and store boundary faces
    const label nOldCells = mesh_.cells().size();

    //- CAP_FACE_USAGE_AUDIT: CSV for all cap-touched faces (local OFstream)
    const bool writeCapFaceAtlas =
        useHardBLBLCapCells_ && !capSideVrtMap_.empty();
    OFstream capFaceAuditOs("blCapCellFaceUsageAtlas.csv");
    if( writeCapFaceAtlas )
        capFaceAuditOs << "bfI,patchName,pKey,nPoints,nCapSide,nCollapsed,"
                       << "fullyCollapsed" << nl;

    //- REDUCED-CELL DRY-RUN ATLAS: classify collapse patterns
    //- Gated by useHardBLBLReducedCells_ (no topology change yet)
    const bool writeReducedDryRun =
        useHardBLBLReducedCells_ && !capSideVrtMap_.empty();
    OFstream capReducedDryRunOs("blCapReducedCellDryRunAtlas.csv");
    if( writeReducedDryRun )
        capReducedDryRunOs
            << "bfI,patchName,pKey,nPoints,nCapSide,nCollapsed,"
            << "cellClass,collapsedMask,"
            << "nCandidateFaces,nValidFaces,nInvalidFaces" << nl;

    //- TWO-PASS CANONICAL FACE SYSTEM
    //- Pass 1: classify all treated-treated internal edges
    //- Pass 2: decide per-face whether to reduce or fallback
    //- Pass 3 (cell loop): build using pre-agreed canonical faces
    enum EdgeState { EDGE_QUAD=0, EDGE_DROP=1, EDGE_TRIANGLE=2, EDGE_UNSAFE=3, EDGE_ASYM_STITCH=4 };
    Map<label> edgeStateMap;
    Map<DynList<label>> edgeCanonicalFace;
    // Per-edge shared top-vertex substitution, computed ONCE during
    // classification (not independently per-cell during construction).
    // Maps: for a given EDGE_TRIANGLE edgeIA, (rawTopLabel -> canonTopLabel)
    // for whichever side's raw top vertex was dropped from the canonical
    // triangle. Both sides look up the SAME entry, guaranteeing identical
    // substitution regardless of which side computes first.
    // Keyed by (edgeIA, rawTopLabel) -- NOT raw label alone, since the
    // same vertex label can appear across multiple different edges,
    // and raw-label-only keying lets one edge's substitution silently
    // overwrite/alias another edge's intended substitution.
    std::map<std::pair<label,label>, label> edgeTopSubst;
    std::map<std::pair<label,label>, std::set<label>> edgeTopSubstDemands;
    // Diagnostic/staging storage for EDGE_ASYM_STITCH candidates.
    // Not used by downstream topology until the activation patch.
    Map<DynList<label>> edgeStitchTriA;
    Map<DynList<label>> edgeStitchTriB;
    boolList faceReducible(bFaces.size(), false);
    boolList triangleBypassFace(bFaces.size(), false);
    label nEdgeQuad=0,nEdgeDrop=0,nEdgeTriangle=0,nEdgeUnsafe=0;
    label nFaceReduce=0,nFaceFallbackTri=0,nFaceFallbackUnsafe=0;
    label nAsymStitchCandidate=0; // diagnostic only -- not yet acted upon

    const bool enableReducedCellTopology =
        useHardBLBLReducedCells_ && !capSideVrtMap_.empty();
    Info << "PAIRKEY_BUILD_MARKER edgeTopSubst_pairkey_v2 active"
         << " enableReducedCellTopology=" << enableReducedCellTopology
         << " useHardBLBLReducedCells=" << useHardBLBLReducedCells_
         << " capSideVrtMapSize=" << capSideVrtMap_.size()
         << endl;

    auto buildTopFaceTP =
    [&](const face& bf, const label bfPatch, DynList<label>& topOut)
    {
        topOut.clear();
        const label pKeyBF = patchKey_[bfPatch];
        forAll(bf, pI2)
        {
            label tl = findNewNodeLabelForPatch(bf[pI2], bfPatch, pKeyBF);
            bool coll = (tl==bf[pI2]);
            if( !coll && tl>=0 && tl<label(mesh_.points().size())
             && bf[pI2]>=0 && bf[pI2]<label(mesh_.points().size()) )
                if( mag(mesh_.points()[tl]-mesh_.points()[bf[pI2]])
                    < scalar(1e-10) ) coll = true;
            topOut.append(coll ? bf[pI2] : tl);
        }
    };

    auto buildSideFaceTP =
    [](const label b0, const label b1,
       const label t0, const label t1,
       DynList<label>& sf)
    {
        sf.clear();
        sf.append(b0);
        if( b1!=b0 ) sf.append(b1);
        if( t1!=b1 && t1!=b0 ) sf.append(t1);
        if( t0!=t1 && t0!=b0 && t0!=b1 ) sf.append(t0);
    };

    auto sameLabelSetTP =
    [](const DynList<label>& a, const DynList<label>& b) -> bool
    {
        if( a.size()!=b.size() ) return false;
        forAll(a,i)
        {
            bool found=false;
            forAll(b,j) if(a[i]==b[j]){found=true;break;}
            if(!found) return false;
        }
        return true;
    };

    // EDGE_ASYM_STITCH: build two triangles that tile a quad-vs-degenerate
    // side-face pair. b0,b1 is the shared base edge (identical on both
    // sides by construction). t0,t1 are the QUAD side's own two cap
    // vertices (t1 corresponds to b1, t0 corresponds to b0). Always
    // computed from the quad side's own labels, in a fixed deterministic
    // order, so both adjacent cells produce byte-identical triangles
    // regardless of which one is doing the computing.
    auto buildStitchTriangles =
    [](const label b0, const label b1, const label t0, const label t1,
       DynList<label>& triA, DynList<label>& triB)
    {
        triA.clear();
        triA.append(b0); triA.append(b1); triA.append(t1);
        triB.clear();
        triB.append(b0); triB.append(t1); triB.append(t0);
    };

    // Geometric sanity check for a stitch triangle: unique labels,
    // nonzero area. Used before ever accepting a stitch triangle into
    // a cell's face list.
    auto validStitchTri =
    [&](const DynList<label>& tri) -> bool
    {
        if( tri.size() != 3 ) return false;
        if( tri[0]==tri[1] || tri[1]==tri[2] || tri[0]==tri[2] ) return false;
        if( tri[0] < 0 || tri[0] >= label(mesh_.points().size()) ) return false;
        if( tri[1] < 0 || tri[1] >= label(mesh_.points().size()) ) return false;
        if( tri[2] < 0 || tri[2] >= label(mesh_.points().size()) ) return false;
        const point& p0 = mesh_.points()[tri[0]];
        const point& p1 = mesh_.points()[tri[1]];
        const point& p2 = mesh_.points()[tri[2]];
        return mag((p1-p0)^(p2-p0)) > VSMALL;
    };

    //- ASYMMETRIC EDGE ATLAS: diagnostic for EDGE_UNSAFE cases
    const bool writeAsymAtlas =
        enableReducedCellTopology && useHardBLBLCapVertexInsertion_;
    OFstream capAsymOs("blCapAsymmetricEdgeAtlas.csv");
    if( writeAsymAtlas )
        capAsymOs
            << "edgeI,bfI,bfPatch,neiFace,neiPatch,"
            << "thisSize,neiSize,unionSize,"
            << "this0,this1,this2,this3,"
            << "nei0,nei1,nei2,nei3,"
            << "union0,union1,union2,union3,union4,"
            << "hasCapThis,hasCapNei,stateGuess" << nl;

    if( enableReducedCellTopology )
    {
        forAll(bFaces, bfI)
        {
            const label bfIPatch = boundaryFacePatches[bfI];
            if( !treatPatches[bfIPatch] ) continue;
            const face& f = bFaces[bfI];
            bool hasCapSide = false;
            forAll(f, pI2)
            {
                const std::pair<label,label> ck(f[pI2], bfIPatch);
                if( capSideVrtMap_.find(ck)!=capSideVrtMap_.end() )
                { hasCapSide=true; break; }
            }
            if( !hasCapSide ) continue;
            DynList<label> topThis;
            buildTopFaceTP(f, bfIPatch, topThis);
            forAll(f, pIA)
            {
                const label edgeIA = faceEdges(bfI, pIA);
                if( edgeFaces.sizeOfRow(edgeIA)!=2 ) continue;
                label neiFaceIA = edgeFaces(edgeIA,0);
                if( neiFaceIA==bfI ) neiFaceIA=edgeFaces(edgeIA,1);
                const label neiPatchIA = boundaryFacePatches[neiFaceIA];
                if( !treatPatches[neiPatchIA] ) continue;
                if( edgeStateMap.found(edgeIA) ) continue;
                const label b0=f[pIA], b1=f.nextLabel(pIA);
                const label t0=topThis[pIA], t1=topThis[(pIA+1)%f.size()];
                DynList<label> sfThis;
                buildSideFaceTP(b0,b1,t0,t1,sfThis);
                const face& fNei = bFaces[neiFaceIA];
                DynList<label> topNei;
                buildTopFaceTP(fNei, neiPatchIA, topNei);
                DynList<label> sfNei;
                label nt0Outer = -1, nt1Outer = -1;
                forAll(fNei, pN)
                {
                    const label nb0=fNei[pN], nb1=fNei.nextLabel(pN);
                    if(!((nb0==b0&&nb1==b1)||(nb0==b1&&nb1==b0))) continue;
                    const label nt0=topNei[pN], nt1=topNei[(pN+1)%fNei.size()];
                    nt0Outer = nt0; nt1Outer = nt1;
                    buildSideFaceTP(nb0,nb1,nt0,nt1,sfNei);
                    break;
                }
                label estate = EDGE_UNSAFE;
                if( sfThis.size()<3 && sfNei.size()<3 )
                { estate=EDGE_DROP; ++nEdgeDrop; }
                else if( sameLabelSetTP(sfThis,sfNei) )
                {
                    if( sfThis.size()==4 )
                    { estate=EDGE_QUAD; ++nEdgeQuad; }
                    else if( sfThis.size()==3 )
                    {
                        estate=EDGE_TRIANGLE; ++nEdgeTriangle;
                        edgeCanonicalFace.insert(edgeIA, sfThis);
                        // Precompute the shared top-vertex substitution
                        // ONCE here, using this side's real t0/t1 and the
                        // neighbor's real nt0/nt1. Both sides will look up
                        // the SAME edgeTopSubst entries later during
                        // construction, instead of each independently
                        // computing findCanonTop from their own local view.
                        auto inCanon = [&](const label v) -> bool
                        {
                            forAll(sfThis, ci)
                            {
                                if( sfThis[ci] == v ) return true;
                            }
                            return false;
                        };
                        label realTop = -1;
                        forAll(sfThis, ci)
                        {
                            if( sfThis[ci] != b0 && sfThis[ci] != b1 )
                            {
                                realTop = sfThis[ci];
                                break;
                            }
                        }
                        if( realTop >= 0 )
                        {
                            const label pKeyCls = patchKey_[bfIPatch];
                            const label freshT0 = findNewNodeLabelForPatch(b0, bfIPatch, pKeyCls);
                            const label freshT1 = findNewNodeLabelForPatch(b1, bfIPatch, pKeyCls);
                            const label pKeyNei = patchKey_[neiPatchIA];
                            const label freshNt0 = findNewNodeLabelForPatch(b0, neiPatchIA, pKeyNei);
                            const label freshNt1 = findNewNodeLabelForPatch(b1, neiPatchIA, pKeyNei);

                            if( !inCanon(freshT0) )
                            {
                                edgeTopSubstDemands[std::make_pair(edgeIA, freshT0)].insert(realTop);
                                Info << "EDGETOPSUBST_DEMAND edgeIA=" << edgeIA
                                     << " raw=" << freshT0 << " realTop=" << realTop
                                     << " side=this0" << endl;
                            }
                            if( !inCanon(freshT1) )
                            {
                                edgeTopSubstDemands[std::make_pair(edgeIA, freshT1)].insert(realTop);
                                Info << "EDGETOPSUBST_DEMAND edgeIA=" << edgeIA
                                     << " raw=" << freshT1 << " realTop=" << realTop
                                     << " side=this1" << endl;
                            }
                            if( !inCanon(freshNt0) )
                            {
                                edgeTopSubstDemands[std::make_pair(edgeIA, freshNt0)].insert(realTop);
                                Info << "EDGETOPSUBST_DEMAND edgeIA=" << edgeIA
                                     << " raw=" << freshNt0 << " realTop=" << realTop
                                     << " side=nei0" << endl;
                            }
                            if( !inCanon(freshNt1) )
                            {
                                edgeTopSubstDemands[std::make_pair(edgeIA, freshNt1)].insert(realTop);
                                Info << "EDGETOPSUBST_DEMAND edgeIA=" << edgeIA
                                     << " raw=" << freshNt1 << " realTop=" << realTop
                                     << " side=nei1" << endl;
                            }
                        }
                    }
                    else { estate=EDGE_UNSAFE; ++nEdgeUnsafe; }
                }
                else
                {
                    estate=EDGE_UNSAFE; ++nEdgeUnsafe;
                    // Quad-vs-degenerate size mismatch: our target stitch
                    // pattern. estate stays EDGE_UNSAFE for now (mutual
                    // reducibility / cell builder behavior unchanged).
                    // We only compute and store the candidate triangles
                    // here so a later, separate activation step can use
                    // them without needing to redo this classification.
                    const bool isQuadVsDegenerate =
                        (sfThis.size()==4 && sfNei.size()==2)
                     || (sfThis.size()==2 && sfNei.size()==4);
                    if( isQuadVsDegenerate )
                    {
                        ++nAsymStitchCandidate;
                        // sfThis/sfNei built via buildSideFaceTP(b0,b1,t0,t1,sf)
                        // which appends in order [b0,b1,t1,t0]. The quad
                        // side is whichever has size 4.
                        const DynList<label>& quadSf =
                            (sfThis.size()==4) ? sfThis : sfNei;
                        if( quadSf.size() == 4 )
                        {
                            const label qb0 = quadSf[0];
                            const label qb1 = quadSf[1];
                            const label qt1 = quadSf[2];
                            const label qt0 = quadSf[3];
                            DynList<label> triA, triB;
                            buildStitchTriangles(qb0,qb1,qt0,qt1,triA,triB);
                            if( validStitchTri(triA) && validStitchTri(triB) )
                            {
                                edgeStitchTriA.insert(edgeIA, triA);
                                edgeStitchTriB.insert(edgeIA, triB);
                            }
                        }
                    }
                }
                edgeStateMap.insert(edgeIA, estate);

                //- Write asymmetric atlas for UNSAFE edges
                if( writeAsymAtlas && estate==EDGE_UNSAFE
                 && sfThis.size()>=3 && sfNei.size()>=3 )
                {
                    //- Build union face: unique labels from both sides
                    DynList<label> unionF;
                    forAll(sfThis,k)
                    {
                        bool dup=false;
                        forAll(unionF,j) if(unionF[j]==sfThis[k]){dup=true;break;}
                        if(!dup) unionF.append(sfThis[k]);
                    }
                    forAll(sfNei,k)
                    {
                        bool dup=false;
                        forAll(unionF,j) if(unionF[j]==sfNei[k]){dup=true;break;}
                        if(!dup) unionF.append(sfNei[k]);
                    }
                    //- Check if either face has a cap-side vertex
                    bool hasCapThis=false, hasCapNei=false;
                    forAll(sfThis,k)
                    {
                        const std::pair<label,label> ck(sfThis[k],bfIPatch);
                        if(capSideVrtMap_.find(ck)!=capSideVrtMap_.end())
                        {hasCapThis=true;break;}
                    }
                    forAll(sfNei,k)
                    {
                        const std::pair<label,label> ck(sfNei[k],neiPatchIA);
                        if(capSideVrtMap_.find(ck)!=capSideVrtMap_.end())
                        {hasCapNei=true;break;}
                    }
                    //- Guess state
                    word stateGuess = "UNSAFE_OTHER";
                    if( sfThis.size()!=sfNei.size()
                     && (hasCapThis||hasCapNei)
                     && unionF.size()<=6 )
                        stateGuess = "ASYM_CAP_CANDIDATE";
                    else if( sfThis.size()!=sfNei.size() )
                        stateGuess = "UNSAFE_SIZE_MISMATCH";
                    else
                        stateGuess = "UNSAFE_LABEL_MISMATCH";
                    capAsymOs
                        << edgeIA << "," << bfI << "," << bfIPatch
                        << "," << neiFaceIA << "," << neiPatchIA
                        << "," << sfThis.size() << "," << sfNei.size()
                        << "," << unionF.size()
                        << "," << (sfThis.size()>0?sfThis[0]:-1)
                        << "," << (sfThis.size()>1?sfThis[1]:-1)
                        << "," << (sfThis.size()>2?sfThis[2]:-1)
                        << "," << (sfThis.size()>3?sfThis[3]:-1)
                        << "," << (sfNei.size()>0?sfNei[0]:-1)
                        << "," << (sfNei.size()>1?sfNei[1]:-1)
                        << "," << (sfNei.size()>2?sfNei[2]:-1)
                        << "," << (sfNei.size()>3?sfNei[3]:-1)
                        << "," << (unionF.size()>0?unionF[0]:-1)
                        << "," << (unionF.size()>1?unionF[1]:-1)
                        << "," << (unionF.size()>2?unionF[2]:-1)
                        << "," << (unionF.size()>3?unionF[3]:-1)
                        << "," << (unionF.size()>4?unionF[4]:-1)
                        << "," << (hasCapThis?"1":"0")
                        << "," << (hasCapNei?"1":"0")
                        << "," << stateGuess << nl;
                }
            }
        }
        Info << "Two-pass edge states: QUAD=" << nEdgeQuad
             << " DROP=" << nEdgeDrop
             << " TRIANGLE=" << nEdgeTriangle
             << " UNSAFE=" << nEdgeUnsafe
             << " ASYMSTITCHCANDIDATE_TOTAL=" << nAsymStitchCandidate
             << " ASYMSTITCH_STORED=" << edgeStitchTriA.size()
             << endl;

        forAll(bFaces, bfI)
        {
            const label bfIPatch = boundaryFacePatches[bfI];
            if( !treatPatches[bfIPatch] ) continue;
            const face& f = bFaces[bfI];
            bool hasCapSide = false;
            forAll(f, pI2)
            {
                const std::pair<label,label> ck(f[pI2], bfIPatch);
                if( capSideVrtMap_.find(ck)!=capSideVrtMap_.end() )
                { hasCapSide=true; break; }
            }
            if( !hasCapSide ) continue;
            //- Check partial collapse: only reduce if some top vertices
            //- collapse (REDUCED_CAP_CELL pattern, not NORMAL_PRISM)
            label nCollapsedP2 = 0;
            const label pKeyP2 = patchKey_[bfIPatch];
            forAll(f, pI2c)
            {
                label tl = findNewNodeLabelForPatch(f[pI2c], bfIPatch, pKeyP2);
                bool coll = (tl==f[pI2c]);
                if( !coll && tl>=0 && tl<label(mesh_.points().size())
                 && f[pI2c]>=0 && f[pI2c]<label(mesh_.points().size()) )
                    if( mag(mesh_.points()[tl]-mesh_.points()[f[pI2c]])
                        < scalar(1e-10) ) coll = true;
                if( coll ) ++nCollapsedP2;
            }
            //- Only REDUCED_CAP_CELL pattern qualifies
            if( nCollapsedP2==0 || nCollapsedP2==label(f.size()) ) continue;

            bool reducible=true;
            forAll(f, pIA)
            {
                const label edgeIA = faceEdges(bfI, pIA);
                if( edgeFaces.sizeOfRow(edgeIA)!=2 ) continue;
                label neiFaceIA = edgeFaces(edgeIA,0);
                if( neiFaceIA==bfI ) neiFaceIA=edgeFaces(edgeIA,1);
                if( !treatPatches[boundaryFacePatches[neiFaceIA]] ) continue;
                if( !edgeStateMap.found(edgeIA) )
                { reducible=false; ++nFaceFallbackUnsafe; break; }
                const label es = edgeStateMap[edgeIA];
                if( es==EDGE_UNSAFE )
                { reducible=false; ++nFaceFallbackUnsafe; break; }
            }
            if( reducible ) { faceReducible[bfI]=true; ++nFaceReduce; }
        }
        Info << "Two-pass face decisions: REDUCE=" << nFaceReduce
             << " FALLBACK_UNSAFE=" << nFaceFallbackUnsafe << endl;

        //- MUTUAL REDUCIBILITY PASS: a face can only reduce if
        //- ALL treated neighbors sharing TRIANGLE/DROP edges are
        //- also reducible. Iterate until stable.
        bool anyChange = true;
        label nMutualFallback = 0;
        while( anyChange )
        {
            anyChange = false;
            forAll(bFaces, bfI)
            {
                if( !faceReducible[bfI] ) continue;
                const face& f = bFaces[bfI];
                forAll(f, pIA)
                {
                    const label edgeIA = faceEdges(bfI, pIA);
                    if( edgeFaces.sizeOfRow(edgeIA)!=2 ) continue;
                    label neiFaceIA = edgeFaces(edgeIA,0);
                    if( neiFaceIA==bfI ) neiFaceIA=edgeFaces(edgeIA,1);
                    if( !treatPatches[boundaryFacePatches[neiFaceIA]] ) continue;
                    if( !edgeStateMap.found(edgeIA) ) continue;
                    const label es = edgeStateMap[edgeIA];

                    if( es==EDGE_TRIANGLE && !faceReducible[neiFaceIA] )
                    {
                        static label nTriangleBypassFaces = 0;
                        // Controlled transition-stencil scale test.
                        // v1 proof used 1.  Try 3 before replacing this
                        // with a true connected-component selector.
                        const label triangleBypassFaceLimit = 3;

                        if( !triangleBypassFace[bfI]
                         && nTriangleBypassFaces < triangleBypassFaceLimit )
                        {
                            triangleBypassFace[bfI] = true;
                            ++nTriangleBypassFaces;
                            Info << "TRIANGLE_BYPASS_FACE"
                                 << " bfI=" << bfI
                                 << " bfPatch=" << boundaryFacePatches[bfI]
                                 << " firstEdgeIA=" << edgeIA
                                 << " firstNeiFace=" << neiFaceIA
                                 << " firstNeiPatch=" << boundaryFacePatches[neiFaceIA]
                                 << " count=" << nTriangleBypassFaces
                                 << endl;
                        }

                        if( triangleBypassFace[bfI] )
                        {
                            continue;
                        }
                    }

                    //- Only check TRIANGLE and DROP edges
                    //- (QUAD edges don't need mutual reducibility)
                    if( es!=EDGE_TRIANGLE && es!=EDGE_DROP ) continue;
                    //- Neighbor must also be reducible
                    if( !faceReducible[neiFaceIA] )
                    {
                        faceReducible[bfI] = false;
                        --nFaceReduce;
                        ++nMutualFallback;
                        anyChange = true;
                        break;
                    }
                }
            }
        }
        Info << "Mutual reducibility pass: "
             << nMutualFallback << " faces demoted, "
             << nFaceReduce << " faces remain reducible" << endl;
    }

    //- Legacy interface atlas disabled
    const bool writeInterfaceAtlas = false;
    OFstream capInterfaceOs("blCapReducedInterfaceAtlas.csv");
    if( writeInterfaceAtlas )
        capInterfaceOs
            << "bfI,bfIPatch,edgeI,pI,neiFace,neiPatch,"
            << "thisSideSize,neiSideSize,compatible,"
            << "this0,this1,this2,this3,"
            << "nei0,nei1,nei2,nei3" << nl;


    // Resolve per-edge substitution demands before construction.
    {
        label nEdgeDemandUnanimous = 0;
        label nEdgeDemandConflict = 0;

        for
        (
            std::map<std::pair<label,label>, std::set<label>>::const_iterator dIt =
                edgeTopSubstDemands.begin();
            dIt != edgeTopSubstDemands.end();
            ++dIt
        )
        {
            const std::pair<label,label>& key = dIt->first;
            const std::set<label>& targets = dIt->second;

            if( targets.size() == 1 )
            {
                edgeTopSubst[key] = *targets.begin();
                ++nEdgeDemandUnanimous;
                Info << "EDGETOPSUBST_RESOLVED"
                     << " edgeIA=" << key.first
                     << " raw=" << key.second
                     << " subst=" << *targets.begin()
                     << endl;
            }
            else
            {
                ++nEdgeDemandConflict;
                Info << "EDGETOPSUBST_CONFLICT"
                     << " edgeIA=" << key.first
                     << " raw=" << key.second
                     << " nTargets=" << label(targets.size())
                     << endl;
            }
        }

        Info << "EDGETOPSUBST_RESOLUTION"
             << " unanimous=" << nEdgeDemandUnanimous
             << " conflict=" << nEdgeDemandConflict
             << endl;

        // CLUSTER-LOCAL ACTIVATION: group resolved edge demands by raw
        // vertex label. A raw vertex is only activated for substitution
        // if EVERY edge referencing it (within edgeTopSubst) resolved to
        // the SAME target -- i.e. it behaves as an isolated 2-sided
        // relationship, not a true multi-way junction. Raw vertices
        // touched by more than one DISTINCT target across their edges
        // are true multi-way junctions and get suppressed entirely
        // (fall back to base-collapse), rather than trying a partial
        // pairwise fix that creates cross-face inconsistency.
        std::map<label, std::set<label>> rawVertexTargets;
        for
        (
            std::map<std::pair<label,label>, label>::const_iterator eIt =
                edgeTopSubst.begin();
            eIt != edgeTopSubst.end();
            ++eIt
        )
        {
            rawVertexTargets[eIt->first.second].insert(eIt->second);
        }
        std::set<label> rawVertexClusterSuppressed;
        label nClusterSafe = 0;
        label nClusterSuppressed = 0;
        for
        (
            std::map<label, std::set<label>>::const_iterator cIt =
                rawVertexTargets.begin();
            cIt != rawVertexTargets.end();
            ++cIt
        )
        {
            if( cIt->second.size() == 1 )
            {
                ++nClusterSafe;
            }
            else
            {
                ++nClusterSuppressed;
                rawVertexClusterSuppressed.insert(cIt->first);
                Info << "EDGETOPSUBST_CLUSTER_SUPPRESSED raw=" << cIt->first
                     << " nTargets=" << label(cIt->second.size())
                     << endl;
            }
        }

        // Diagnostic export: suppressed multi-target raw vertices.
        // These are the Q3 cases deliberately falling back to base-collapse.
        // Load this CSV in ParaView and overlay it on the visible
        // protrusions/pinches to test whether the remaining artifact is
        // located at the suppressed multi-way clusters.
        {
            OFstream suppressedOs("suppressedClusterVertices.csv");
            suppressedOs << "vertexLabel,x,y,z,nTargets,targetLabels" << nl;

            for
            (
                std::set<label>::const_iterator sIt3 =
                    rawVertexClusterSuppressed.begin();
                sIt3 != rawVertexClusterSuppressed.end();
                ++sIt3
            )
            {
                const label v = *sIt3;
                if( v >= 0 && v < label(mesh_.points().size()) )
                {
                    const point& p = mesh_.points()[v];
                    const std::set<label>& tgts = rawVertexTargets[v];

                    suppressedOs << v << ',' << p.x() << ',' << p.y()
                                 << ',' << p.z() << ','
                                 << label(tgts.size()) << ',' << '"';

                    label ti = 0;
                    for
                    (
                        std::set<label>::const_iterator tIt = tgts.begin();
                        tIt != tgts.end();
                        ++tIt
                    )
                    {
                        if( ti++ ) suppressedOs << ';';
                        suppressedOs << *tIt;
                    }

                    suppressedOs << '"' << nl;
                }
            }

            Info << "Suppressed cluster vertex coordinates written to"
                 << " suppressedClusterVertices.csv ("
                 << rawVertexClusterSuppressed.size() << " vertices)"
                 << endl;
        }
        // Remove suppressed raw-vertex entries from edgeTopSubst so
        // construction's lookup naturally falls back to base-collapse
        // for them, while safe (isolated, unanimous) entries remain
        // active and get applied.
        for
        (
            std::map<std::pair<label,label>, label>::iterator eIt2 =
                edgeTopSubst.begin();
            eIt2 != edgeTopSubst.end();
        )
        {
            if( rawVertexClusterSuppressed.count(eIt2->first.second) > 0 )
                eIt2 = edgeTopSubst.erase(eIt2);
            else
                ++eIt2;
        }
        Info << "EDGETOPSUBST_CLUSTER_RESOLUTION"
             << " safe=" << nClusterSafe
             << " suppressed=" << nClusterSuppressed
             << " remainingActive=" << edgeTopSubst.size()
             << endl;
    }

    // ---- BLCELLGEOM_AUDIT counters (report-only) ----
    const label clqNP = mesh_.boundaries().size();
    labelList clqTot(clqNP,0);
    labelList clqPlainInvalid(clqNP,0), clqPatchInvalid(clqNP,0);
    labelList clqPlainSame(clqNP,0),    clqPatchSame(clqNP,0);
    labelList clqPlainNear(clqNP,0),    clqPatchNear(clqNP,0);
    labelList clqDifferent(clqNP,0),    clqCapTouched(clqNP,0);
    labelList clqPlainFullColl(clqNP,0), clqPatchFullColl(clqNP,0);
    labelList clqPlainPartial(clqNP,0),  clqPatchPartial(clqNP,0);
    labelList clqPlainHealthy(clqNP,0),  clqPatchHealthy(clqNP,0);
    scalarField clqPlainMin(clqNP,GREAT), clqPatchMin(clqNP,GREAT);
    scalarField clqPlainSum(clqNP,0.0),   clqPatchSum(clqNP,0.0);
    labelList clqPlainNvalid(clqNP,0),    clqPatchNvalid(clqNP,0);
    label clqHubExamples = 0;
    // ---- end BLCELLGEOM_AUDIT counters ----

    forAll(bFaces, bfI)
    {
        if( treatPatches[boundaryFacePatches[bfI]] )
        {
            const face& f = bFaces[bfI];

            const label pKey = patchKey_[boundaryFacePatches[bfI]];
            const label bfIPatch = boundaryFacePatches[bfI];

            // ---- BLCELLGEOM_AUDIT per-face measurement ----
            {
                const label clqNPts = mesh_.points().size();

                if( bfIPatch < 0 || bfIPatch >= clqNP )
                {
                    Info << "BLCELLGEOM_AUDIT INVALID_PATCH"
                         << " bfI=" << bfI
                         << " bfIPatch=" << bfIPatch << endl;
                    continue;
                }

                const label pp = bfIPatch;
                label pInv=0,pSame=0,pNear=0,qInv=0,qSame=0,qNear=0,nDiff=0,nCap=0;
                forAll(f, pI)
                {
                    const label base = f[pI];
                    const label plain = findNewNodeLabel(base, pKey);
                    const label patchL = findNewNodeLabelForPatch(base, bfIPatch, pKey);
                    if( plain != patchL ) ++nDiff;
                    const std::pair<label,label> capKey(base, bfIPatch);
                    if( capSideVrtMap_.find(capKey) != capSideVrtMap_.end() ) ++nCap;
                    if( plain < 0 || plain >= clqNPts )
                    {
                        ++pInv;
                    }
                    else
                    {
                        const scalar d = (base>=0 && base<clqNPts)
                            ? mag(mesh_.points()[plain]-mesh_.points()[base])
                            : scalar(0);
                        if( plain == base ) ++pSame;
                        else if( d < scalar(1e-10) ) ++pNear;
                        clqPlainMin[pp] = Foam::min(clqPlainMin[pp], d);
                        clqPlainSum[pp] += d;
                        ++clqPlainNvalid[pp];
                    }
                    if( patchL < 0 || patchL >= clqNPts )
                    {
                        ++qInv;
                    }
                    else
                    {
                        const scalar d = (base>=0 && base<clqNPts)
                            ? mag(mesh_.points()[patchL]-mesh_.points()[base])
                            : scalar(0);
                        if( patchL == base ) ++qSame;
                        else if( d < scalar(1e-10) ) ++qNear;
                        clqPatchMin[pp] = Foam::min(clqPatchMin[pp], d);
                        clqPatchSum[pp] += d;
                        ++clqPatchNvalid[pp];
                    }
                }
                const label nv = f.size();
                ++clqTot[pp];
                clqPlainInvalid[pp]+=pInv; clqPatchInvalid[pp]+=qInv;
                clqPlainSame[pp]+=pSame;   clqPatchSame[pp]+=qSame;
                clqPlainNear[pp]+=pNear;   clqPatchNear[pp]+=qNear;
                clqDifferent[pp]+=nDiff;   clqCapTouched[pp]+=(nCap>0?1:0);
                const label pBad = pInv + pSame + pNear;
                const label qBad = qInv + qSame + qNear;

                if( pBad == nv ) ++clqPlainFullColl[pp];
                else if( pBad == 0 ) ++clqPlainHealthy[pp];
                else ++clqPlainPartial[pp];

                if( qBad == nv ) ++clqPatchFullColl[pp];
                else if( qBad == 0 ) ++clqPatchHealthy[pp];
                else ++clqPatchPartial[pp];
                const word clqPnm =
                    (bfIPatch>=0 && bfIPatch<label(patchNames_.size()))
                    ? patchNames_[bfIPatch] : word("?");
                if( clqPnm == word("hub") && clqHubExamples < 5 )
                {
                    ++clqHubExamples;
                    Info << "BLCELLGEOM_HUBEX face=" << bfI << " pKey=" << pKey;
                    forAll(f, pI)
                    {
                        const label base=f[pI];
                        const label plain=findNewNodeLabel(base,pKey);
                        const label patchL=findNewNodeLabelForPatch(base,bfIPatch,pKey);
                        scalar dp =
                            (plain>=0 && plain<clqNPts && base>=0 && base<clqNPts)
                            ? mag(mesh_.points()[plain]-mesh_.points()[base])
                            : scalar(-1);

                        scalar dq =
                            (patchL>=0 && patchL<clqNPts && base>=0 && base<clqNPts)
                            ? mag(mesh_.points()[patchL]-mesh_.points()[base])
                            : scalar(-1);

                        Info << " | base=" << base
                             << " plain=" << plain << "(d=" << dp << ")"
                             << " patch=" << patchL << "(d=" << dq << ")";
                    }
                    Info << endl;
                }
            }
            // ---- end BLCELLGEOM_AUDIT per-face measurement ----
            //- Patch-aware top vertex resolver: checks capSideVrtMap_
            //- keyed by (meshPointI, patchI) first, then falls back
            //- to findNewNodeLabel. Required because hub+blade share
            //- pKey=0 in the same treatPatchesWithPatch_ group.
            auto capAwareTopLabel = [&](const label baseLabel) -> label
            {
                return findNewNodeLabelForPatch(baseLabel, bfIPatch, pKey);
            };

            // Guard 0: skip faces with missing extruded vertices.
            // findNewNodeLabel returns -1 for transition-zone points
            // that never got a new vertex created. Building a layer
            // cell from such a face creates degenerate zero-volume cells.
            {
                bool missingVertex = false;
                forAll(f, pI)
                {
                    if( findNewNodeLabel(f[pI], pKey) < 0 )
                    {
                        missingVertex = true;
                        break;
                    }
                }
                if( missingVertex )
                {
                    static label nMissing = 0;
                    if( ++nMissing <= 100 )
                        Info << "Skipping layer face: missing extruded"
                             << " vertex bfI=" << bfI
                             << " pKey=" << pKey << endl;
                    continue;
                }
            }

            // Franjo TODO: build proper wedge/pyramid/reduced cells
            // at sharp BL/BL feature curves instead of collapsed prisms.
            // Implementation: appendValidFace canonicalizer + reduced
            // topology for REDUCED_CAP_CELL pattern faces.
            // Gate: useHardBLBLReducedCells_ (default false).

            // CAP_FACE_USAGE_AUDIT: write to local CSV OFstream
            if( writeCapFaceAtlas )
            {
                label nCapSide = 0;
                label nCollapsed = 0;
                forAll(f, pI)
                {
                    const label topLabel = findNewNodeLabelForPatch(
                        f[pI], boundaryFacePatches[bfI], pKey);
                    //- Cap-side detection: resolved label differs from base
                    {
                        const std::pair<label,label> capKey(f[pI], boundaryFacePatches[bfI]);
                        if( capSideVrtMap_.find(capKey) != capSideVrtMap_.end() )
                            ++nCapSide;
                    }
                    bool collapsed = false;
                    if( topLabel == f[pI] ) collapsed = true;
                    if( !collapsed
                     && topLabel >= 0
                     && topLabel < label(mesh_.points().size())
                     && f[pI] >= 0
                     && f[pI] < label(mesh_.points().size()) )
                    {
                        if( mag(mesh_.points()[topLabel]-mesh_.points()[f[pI]])
                            < scalar(1e-10) ) collapsed = true;
                    }
                    if( collapsed ) ++nCollapsed;
                }
                if( nCapSide > 0 )
                {
                    const label patchI = boundaryFacePatches[bfI];
                    const word pName =
                        (patchI>=0 && patchI<label(patchNames_.size())) ?
                        patchNames_[patchI] : word("?");
                    capFaceAuditOs
                        << bfI << "," << pName << "," << pKey << ","
                        << f.size() << "," << nCapSide << "," << nCollapsed << ","
                        << (nCollapsed==label(f.size()) ? "1" : "0") << nl;

                    if( writeReducedDryRun )
                    {
                        word cellClass = "NORMAL_PRISM";
                        if( nCollapsed == label(f.size()) )
                            cellClass = "INVALID_FULL_COLLAPSE";
                        else if( nCollapsed > 0 )
                            cellClass = "REDUCED_CAP_CELL";

                        //- Build collapsed mask and top face
                        word collapsedMask = "";
                        DynList<label> topFace;
                        forAll(f, pI)
                        {
                            const label topLabel = findNewNodeLabelForPatch(
                                f[pI], boundaryFacePatches[bfI], pKey);
                            bool coll = (topLabel == f[pI]);
                            if( !coll && topLabel>=0
                             && topLabel<label(mesh_.points().size())
                             && f[pI]>=0
                             && f[pI]<label(mesh_.points().size()) )
                                if( mag(mesh_.points()[topLabel]-mesh_.points()[f[pI]])
                                    < scalar(1e-10) ) coll = true;
                            collapsedMask += (coll ? "1" : "0");
                            topFace.append(coll ? f[pI] : topLabel);
                        }

                        label nCandidateFaces = 0;
                        label nValidFaces = 0;

                        //- Helper lambda: unique-label + area check
                        auto faceValid = [&](const DynList<label>& sf) -> bool
                        {
                            if( sf.size() < 3 ) return false;
                            //- Unique label check
                            forAll(sf, i)
                                for(label j=0;j<i;++j)
                                    if(sf[i]==sf[j]) return false;
                            //- Area check
                            vector areaVec = vector::zero;
                            const point& p0a = mesh_.points()[sf[0]];
                            for(label ai=1;ai<sf.size()-1;++ai)
                                areaVec += (mesh_.points()[sf[ai]]-p0a)
                                         ^ (mesh_.points()[sf[ai+1]]-p0a);
                            return mag(areaVec) > VSMALL;
                        };

                        //- Base face
                        ++nCandidateFaces;
                        {
                            DynList<label> bf2;
                            forAll(f, pI) bf2.append(f[pI]);
                            if( faceValid(bf2) ) ++nValidFaces;
                        }

                        //- Top face
                        ++nCandidateFaces;
                        {
                            DynList<label> tf;
                            forAll(topFace, pI) tf.append(topFace[pI]);
                            if( faceValid(tf) ) ++nValidFaces;
                        }

                        //- Side faces
                        forAll(f, pI)
                        {
                            ++nCandidateFaces;
                            const label b0 = f[pI];
                            const label b1 = f.nextLabel(pI);
                            const label t1 = topFace[(pI+1)%topFace.size()];
                            const label t0 = topFace[pI];
                            DynList<label> sf;
                            sf.append(b0);
                            if( b1!=b0 ) sf.append(b1);
                            if( t1!=b1 && t1!=b0 ) sf.append(t1);
                            if( t0!=t1 && t0!=b0 && t0!=b1 ) sf.append(t0);
                            if( faceValid(sf) ) ++nValidFaces;
                        }

                        const label nInvalidFaces =
                            nCandidateFaces - nValidFaces;
                        capReducedDryRunOs
                            << bfI << "," << pName << "," << pKey << ","
                            << f.size() << "," << nCapSide << "," << nCollapsed << ","
                            << cellClass << "," << collapsedMask << ","
                            << nCandidateFaces << "," << nValidFaces
                            << "," << nInvalidFaces << nl;
                    }
                }
            }

            //- INTERFACE ATLAS: check treated-treated internal edges
            if( writeInterfaceAtlas )
            {
                const label bfIPatchIA = boundaryFacePatches[bfI];
                //- Only process cap-touched faces
                bool hasCapSide = false;
                forAll(f, pI2)
                {
                    const std::pair<label,label> ck(f[pI2], bfIPatchIA);
                    if( capSideVrtMap_.find(ck) != capSideVrtMap_.end() )
                    { hasCapSide = true; break; }
                }
                if( hasCapSide )
                {
                    const label pKeyIA = patchKey_[bfIPatchIA];
                    //- Build top face for this bfI using patch-aware lookup
                    DynList<label> topIA;
                    forAll(f, pI2)
                    {
                        label tl = findNewNodeLabelForPatch(
                            f[pI2], bfIPatchIA, pKeyIA);
                        bool coll = (tl == f[pI2]);
                        if( !coll && tl>=0 && tl<label(mesh_.points().size())
                         && f[pI2]>=0 && f[pI2]<label(mesh_.points().size()) )
                            if( mag(mesh_.points()[tl]-mesh_.points()[f[pI2]])
                                < scalar(1e-10) ) coll = true;
                        topIA.append(coll ? f[pI2] : tl);
                    }
                    //- Helper: simplify side face (same rule as builder)
                    auto simpleSide =
                    [](const label b0, const label b1,
                       const label t1, const label t0,
                       DynList<label>& sf)
                    {
                        sf.clear();
                        sf.append(b0);
                        if( b1!=b0 ) sf.append(b1);
                        if( t1!=b1 && t1!=b0 ) sf.append(t1);
                        if( t0!=t1 && t0!=b0 && t0!=b1 ) sf.append(t0);
                    };
                    //- Helper: same label set (no operator== on HashSet)
                    auto sameLabelSet =
                    [](const DynList<label>& a,
                       const DynList<label>& b) -> bool
                    {
                        if( a.size() != b.size() ) return false;
                        forAll(a, i)
                        {
                            bool found = false;
                            forAll(b, j) if(a[i]==b[j]){found=true;break;}
                            if(!found) return false;
                        }
                        return true;
                    };
                    forAll(f, pIA)
                    {
                        const label edgeIA = faceEdges(bfI, pIA);
                        if( edgeFaces.sizeOfRow(edgeIA) != 2 ) continue;
                        label neiFaceIA = edgeFaces(edgeIA, 0);
                        if( neiFaceIA == bfI )
                            neiFaceIA = edgeFaces(edgeIA, 1);
                        const label neiPatchIA = boundaryFacePatches[neiFaceIA];
                        if( !treatPatches[neiPatchIA] ) continue;
                        //- This side face
                        const label b0 = f[pIA];
                        const label b1 = f.nextLabel(pIA);
                        const label t1 = topIA[(pIA+1)%f.size()];
                        const label t0 = topIA[pIA];
                        DynList<label> sfThis;
                        simpleSide(b0, b1, t1, t0, sfThis);
                        //- Neighbor top face
                        const face& fNei = bFaces[neiFaceIA];
                        const label neiPKeyIA = patchKey_[neiPatchIA];
                        DynList<label> topNei;
                        forAll(fNei, pN)
                        {
                            label tlN = findNewNodeLabelForPatch(
                                fNei[pN], neiPatchIA, neiPKeyIA);
                            bool cN = (tlN == fNei[pN]);
                            if( !cN && tlN>=0 && tlN<label(mesh_.points().size())
                             && fNei[pN]>=0 && fNei[pN]<label(mesh_.points().size()) )
                                if( mag(mesh_.points()[tlN]-mesh_.points()[fNei[pN]])
                                    < scalar(1e-10) ) cN = true;
                            topNei.append(cN ? fNei[pN] : tlN);
                        }
                        //- Find matching edge in neighbor and build its side
                        DynList<label> sfNei;
                        forAll(fNei, pN)
                        {
                            const label nb0 = fNei[pN];
                            const label nb1 = fNei.nextLabel(pN);
                            if( !((nb0==b0&&nb1==b1)||(nb0==b1&&nb1==b0)) )
                                continue;
                            const label nt1 = topNei[(pN+1)%fNei.size()];
                            const label nt0 = topNei[pN];
                            simpleSide(nb0, nb1, nt1, nt0, sfNei);
                            break;
                        }
                        const bool compat =
                            sfThis.size() >= 3
                         && sfNei.size() >= 3
                         && sameLabelSet(sfThis, sfNei);
                        capInterfaceOs
                            << bfI << "," << bfIPatchIA
                            << "," << edgeIA << "," << pIA
                            << "," << neiFaceIA << "," << neiPatchIA
                            << "," << sfThis.size() << "," << sfNei.size()
                            << "," << (compat?"1":"0")
                            << "," << (sfThis.size()>0?sfThis[0]:-1)
                            << "," << (sfThis.size()>1?sfThis[1]:-1)
                            << "," << (sfThis.size()>2?sfThis[2]:-1)
                            << "," << (sfThis.size()>3?sfThis[3]:-1)
                            << "," << (sfNei.size()>0?sfNei[0]:-1)
                            << "," << (sfNei.size()>1?sfNei[1]:-1)
                            << "," << (sfNei.size()>2?sfNei[2]:-1)
                            << "," << (sfNei.size()>3?sfNei[3]:-1)
                            << nl;
                    }
                }
            }

            DynList<DynList<label> > cellFaces;
            DynList<label> newF;

            //- Reduced cap cell builder (Franjo TODO implementation).
            //- Gate: useHardBLBLReducedCells_ (default false).
            bool builtReducedCell = false;
            if( enableReducedCellTopology && useHardBLBLReducedCells_
             && !capSideVrtMap_.empty() && faceReducible[bfI] )
            {
                label nCapSideRC = 0;
                label nCollapsedRC = 0;
                DynList<label> topFaceRC;
                //- Use outer capAwareTopLabel lambda (patchI-keyed)
                auto isCollapsedTop = [&](const label base, const label top) -> bool
                {
                    if( top == base ) return true;
                    if( top>=0 && top<label(mesh_.points().size())
                     && base>=0 && base<label(mesh_.points().size()) )
                        return mag(mesh_.points()[top]-mesh_.points()[base])
                               < scalar(1e-10);
                    return true;
                };
                forAll(f, pI)
                {
                    const label baseLabel = f[pI];
                    const label topLabel = capAwareTopLabel(baseLabel);
                    //- Cap-side detection via map lookup
                    {
                        const std::pair<label,label> capKey(baseLabel, bfIPatch);
                        if( capSideVrtMap_.find(capKey) != capSideVrtMap_.end() )
                            ++nCapSideRC;
                    }
                    const bool coll = isCollapsedTop(baseLabel, topLabel);
                    if( coll ) ++nCollapsedRC;
                    topFaceRC.append(coll ? baseLabel : topLabel);
                }
                if( nCapSideRC > 0 && nCollapsedRC > 0
                 && nCollapsedRC < label(f.size()) )
                {
                    //- Snapshot boundary face lists for rollback
                    const label nBndFacesSnap = newBoundaryFaces.size();
                    const label nBndOwnSnap = newBoundaryOwners.size();
                    const label nBndPatchSnap = newBoundaryPatches.size();

                    //- makeValidFace: simplify duplicates first, then
                    //- validate unique-label count and nonzero area.
                    //- Returns cleaned face in validFace, true if valid.
                    auto makeValidFace =
                    [&](const DynList<label>& cand,
                        DynList<label>& validFace) -> bool
                    {
                        validFace.clear();
                        //- Remove consecutive duplicates
                        forAll(cand, i)
                        {
                            const label lbl = cand[i];
                            if( validFace.size()==0
                             || validFace[validFace.size()-1] != lbl )
                                validFace.append(lbl);
                        }
                        //- Remove last if same as first (closed loop)
                        if( validFace.size() > 1
                         && validFace[validFace.size()-1] == validFace[0] )
                            validFace.setSize(validFace.size()-1);
                        if( validFace.size() < 3 ) return false;
                        //- Unique label check
                        forAll(validFace, i)
                            for(label j=0;j<i;++j)
                                if(validFace[i]==validFace[j]) return false;
                        //- Nonzero area check
                        vector areaVec = vector::zero;
                        const point& p0a = mesh_.points()[validFace[0]];
                        for(label ai=1;ai<validFace.size()-1;++ai)
                            areaVec += (mesh_.points()[validFace[ai]]-p0a)
                                     ^ (mesh_.points()[validFace[ai+1]]-p0a);
                        return mag(areaVec) > VSMALL;
                    };

                    //- Base face (reversed boundary face)
                    DynList<label> baseCand;
                    baseCand.append(f[0]);
                    for(label pI=f.size()-1;pI>0;--pI)
                        baseCand.append(f[pI]);
                    DynList<label> validBaseF;
                    if( makeValidFace(baseCand, validBaseF) )
                        cellFaces.append(validBaseF);

                    //- Top face (collapsed top_i replaced by base_i)
                    DynList<label> validTopF;
                    const bool topOK = makeValidFace(topFaceRC, validTopF);
                    if( topOK )
                    {
                        cellFaces.append(validTopF);
                        newBoundaryFaces.appendList(validTopF);
                        newBoundaryOwners.append(cellsToAdd.size() + nOldCells);
                        newBoundaryPatches.append(boundaryFacePatches[bfI]);
                    }

                    //- Side faces + edge-boundary patch logic
                    bool internalSideMismatch = false;
                    forAll(f, pI)
                    {
                        if( internalSideMismatch ) break;
                        const label b0 = f[pI];
                        const label b1 = f.nextLabel(pI);
                        const label t1 = topFaceRC[(pI+1)%f.size()];
                        const label t0 = topFaceRC[pI];
                        DynList<label> sideCand;
                        sideCand.append(b0);
                        sideCand.append(b1);
                        sideCand.append(t1);
                        sideCand.append(t0);
                        const label edgeI = faceEdges(bfI, pI);
                        label neiFaceTI = -1;
                        bool treatedInternal = false;
                        if( edgeFaces.sizeOfRow(edgeI) == 2 )
                        {
                            neiFaceTI = edgeFaces(edgeI, 0);
                            if( neiFaceTI == bfI ) neiFaceTI = edgeFaces(edgeI, 1);
                            if( treatPatches[boundaryFacePatches[neiFaceTI]] )
                                treatedInternal = true;
                        }
                        DynList<label> validSideF;
                        const bool sideOK = makeValidFace(sideCand, validSideF);
                        if( treatedInternal && neiFaceTI >= 0 )
                        {
                            //- Use pre-computed edge state from two-pass pre-pass
                            if( !edgeStateMap.found(edgeI) )
                            { internalSideMismatch=true; break; }
                            const label preEs = edgeStateMap[edgeI];
                            if( preEs==EDGE_UNSAFE )
                            { internalSideMismatch=true; break; }
                            if( preEs==EDGE_DROP )
                            { continue; } //- drop degenerate side face
                            if( preEs==EDGE_TRIANGLE )
                            {
                                //- Two-pass pre-pass confirmed compatible.
                                //- faceReducible[bfI] ensures neighbor also reduces.
                                //- Use local validSideF for correct per-cell orientation.
                                if( !sideOK || validSideF.size() != 3 )
                                { internalSideMismatch=true; break; }
                                cellFaces.append(validSideF);
                                continue;
                            }
                            if( preEs==EDGE_ASYM_STITCH )
                            {
                                //- Emit the two pre-computed, pre-validated stitch
                                //- triangles (same stored labels used by whichever
                                //- side -- this face or neighbor -- looks them up
                                //- for this edgeI, guaranteeing matching internal
                                //- boundary labels on both cells).
                                if( !edgeStitchTriA.found(edgeI)
                                 || !edgeStitchTriB.found(edgeI) )
                                { internalSideMismatch=true; break; }
                                cellFaces.append(edgeStitchTriA[edgeI]);
                                cellFaces.append(edgeStitchTriB[edgeI]);
                                continue;
                            }
                            //- EDGE_QUAD: fall through to sideOK append
                        }
                        if( sideOK )
                        {
                            cellFaces.append(validSideF);
                            if( edgeFaces.sizeOfRow(edgeI) == 2 )
                            {
                                label neiFace2 = edgeFaces(edgeI, 0);
                                if( neiFace2 == bfI ) neiFace2 = edgeFaces(edgeI, 1);
                                if( !treatPatches[boundaryFacePatches[neiFace2]] )
                                {
                                    newBoundaryFaces.appendList(validSideF);
                                    newBoundaryOwners.append(cellsToAdd.size() + nOldCells);
                                    newBoundaryPatches.append(boundaryFacePatches[neiFace2]);
                                }
                            }
                            else if( edgeFaces.sizeOfRow(edgeI) == 1 )
                            {
                                const Map<label>& otherProcPatch = *otherProcPatchPtr;
                                if( !treatPatches[otherProcPatch[edgeI]] )
                                {
                                    newBoundaryFaces.appendList(validSideF);
                                    newBoundaryOwners.append(cellsToAdd.size() + nOldCells);
                                    newBoundaryPatches.append(otherProcPatch[edgeI]);
                                }
                            }
                        }
                    }
                    static label nReducedAttempt = 0;
                    static label nReducedBuilt = 0;
                    static label nRollbackTopo = 0;
                    static label nRollbackGeom = 0;
                    static label nRollbackInternal = 0;
                    static label nRollbackFewFaces = 0;
                    ++nReducedAttempt;

                    //- Check 1: topological closure
                    //- Every undirected edge must appear exactly twice.
                    auto cellLooksClosed =
                    [&](const DynList<DynList<label> >& cFaces) -> bool
                    {
                        std::map<std::pair<label,label>,label> edgeUse;
                        forAll(cFaces, fi)
                        {
                            const DynList<label>& cf = cFaces[fi];
                            if( cf.size() < 3 ) return false;
                            forAll(cf, i)
                            {
                                const label a = cf[i];
                                const label b = cf[(i+1)%cf.size()];
                                if( a==b ) return false;
                                ++edgeUse[std::make_pair(Foam::min(a,b),Foam::max(a,b))];
                            }
                        }
                        for( auto it=edgeUse.begin(); it!=edgeUse.end(); ++it )
                            if( it->second != 2 ) return false;
                        return true;
                    };

                    //- Check 2: geometric openness
                    //- Sum of oriented face area vectors must be ~0.
                    auto cellGeomClosed =
                    [&](const DynList<DynList<label> >& cFaces) -> bool
                    {
                        vector areaSum = vector::zero;
                        scalar areaMagSum = scalar(0);
                        forAll(cFaces, fi)
                        {
                            const DynList<label>& cf = cFaces[fi];
                            vector fArea = vector::zero;
                            const point& p0f = mesh_.points()[cf[0]];
                            for(label ai=1;ai<cf.size()-1;++ai)
                                fArea += (mesh_.points()[cf[ai]]-p0f)
                                       ^ (mesh_.points()[cf[ai+1]]-p0f);
                            areaSum += fArea;
                            areaMagSum += mag(fArea);
                        }
                        if( areaMagSum < VSMALL ) return false;
                        return mag(areaSum)/areaMagSum < scalar(1e-4);
                    };

                    const bool fewFacesOK = (cellFaces.size() >= 4);
                    const bool topoOK = fewFacesOK
                        && !internalSideMismatch
                        && cellLooksClosed(cellFaces);
                    const bool geomOK = topoOK && cellGeomClosed(cellFaces);

                    if( !fewFacesOK ) ++nRollbackFewFaces;
                    else if( internalSideMismatch ) ++nRollbackInternal;
                    else if( !topoOK ) ++nRollbackTopo;
                    else if( !geomOK ) ++nRollbackGeom;

                    if( fewFacesOK && !internalSideMismatch && topoOK && geomOK )
                    {
                        builtReducedCell = true;
                        ++nReducedBuilt;
                        if
                        (
                            nReducedAttempt == 1
                         || nReducedAttempt % 100 == 0
                         || nReducedBuilt == 1
                         || nReducedBuilt % 100 == 0
                        )
                        {
                            Info << "Reduced cap cells: attempt=" << nReducedAttempt
                                 << " built=" << nReducedBuilt
                                 << " rollbackTopo=" << nRollbackTopo
                                 << " rollbackGeom=" << nRollbackGeom
                                 << " rollbackInternal=" << nRollbackInternal
                                 << " rollbackFewFaces=" << nRollbackFewFaces
                                 << endl;
                        }
                    }
                    else
                    {
                        //- Roll back all speculative additions
                        cellFaces.clear();
                        newBoundaryFaces.setSize(nBndFacesSnap);
                        newBoundaryOwners.setSize(nBndOwnSnap);
                        newBoundaryPatches.setSize(nBndPatchSnap);
                    }
                }
            }

            bool normalCellUsedCanonicalTri = false;
            bool normalCanonTriTopoOK = true;
            bool normalCanonTriGeomOK = true;
            std::map<label,label> normalTopSubst;
            label normalTopCellFaceI = -1;
            label normalTopBoundaryRow = -1;

            if( !builtReducedCell )
            {
            //- Normal prism path
            //- store the current boundary face
            newF.clear();
            newF.append(f[0]);
            for(label pI=f.size()-1;pI>0;--pI)
                newF.append(f[pI]);
            cellFaces.append(newF);
            //- create parallel face
            forAll(f, pI)
                newF[pI] = findNewNodeLabel(f[pI], pKey);
            normalTopCellFaceI = cellFaces.size();
            normalTopBoundaryRow = newBoundaryFaces.size();
            cellFaces.append(newF);
            newBoundaryFaces.appendList(newF);
            newBoundaryOwners.append(cellsToAdd.size() + nOldCells);
            newBoundaryPatches.append(boundaryFacePatches[bfI]);

            //- create quad faces
            newF.setSize(4);
            forAll(f, pI)
            {
                newF[0] = f[pI];
                newF[1] = f.nextLabel(pI);
                newF[2] = findNewNodeLabel(newF[1], pKey);
                newF[3] = findNewNodeLabel(f[pI], pKey);

                // Pairing fix experiment: if this normal-prism side
                // edge is a treated-treated EDGE_TRIANGLE, emit the same
                // canonical 3-label side face used by the reduced-cell
                // path. Otherwise one side can emit a reduced triangle
                // while the other emits a raw quad, creating addCells
                // orphans and open cells.
                const label edgeIA = faceEdges(bfI, pI);
                bool emittedCanonicalTri = false;

                if
                (
                    edgeFaces.sizeOfRow(edgeIA) == 2
                 && edgeStateMap.found(edgeIA)
                 && edgeStateMap[edgeIA] == EDGE_TRIANGLE
                 && edgeCanonicalFace.found(edgeIA)
                )
                {
                    label neiFaceIA = edgeFaces(edgeIA, 0);
                    if( neiFaceIA == bfI ) neiFaceIA = edgeFaces(edgeIA, 1);

                    if( treatPatches[boundaryFacePatches[neiFaceIA]]
                     && faceReducible[neiFaceIA] )
                    {
                        // edgeTopSubst is now pre-resolved before construction.
                        // Do not populate it here; construction only reads it.
                        static label nNormalSideCanonTri = 0;
                        ++nNormalSideCanonTri;
                        if( nNormalSideCanonTri <= 20
                         || nNormalSideCanonTri % 100 == 0 )
                        {
                            Info << "NORMAL_SIDE_CANON_TRI"
                                 << " bfI=" << bfI
                                 << " bfPatch=" << boundaryFacePatches[bfI]
                                 << " edgeIA=" << edgeIA
                                 << " neiFace=" << neiFaceIA
                                 << " neiPatch=" << boundaryFacePatches[neiFaceIA]
                                 << " count=" << nNormalSideCanonTri
                                 << endl;
                        }

                        const DynList<label>& canonF = edgeCanonicalFace[edgeIA];

                        auto canonContains = [&](const label v) -> bool
                        {
                            forAll(canonF, ci)
                            {
                                if( canonF[ci] == v ) return true;
                            }
                            return false;
                        };

                        // Raw side quad convention here is [b0,b1,t1,t0].
                        // canonF is a verified 3-vertex triangle. If our
                        // raw top vertex isn't in canonF, substitute to
                        // the ACTUAL correct top vertex in canonF (not
                        // base), guarded by verifying canonF genuinely
                        // contains our two base labels plus one top.
                        auto findCanonTop = [&](const label b0v, const label b1v) -> label
                        {
                            forAll(canonF, ci)
                            {
                                if( canonF[ci] != b0v && canonF[ci] != b1v )
                                    return canonF[ci];
                            }
                            return label(-1);
                        };
                        {
                            label nBaseInCanon = 0;
                            forAll(canonF, ci)
                            {
                                if( canonF[ci] == newF[0] || canonF[ci] == newF[1] )
                                    ++nBaseInCanon;
                            }
                            const label probeCanonTop =
                                findCanonTop(newF[0], newF[1]);
                            if( nBaseInCanon != 2 || probeCanonTop < 0 )
                            {
                                Info << "NORMAL_CANON_TRI_SUBST_NO_CANONTOP"
                                     << " bfI=" << bfI
                                     << " edgeIA=" << edgeIA
                                     << " b0=" << newF[0]
                                     << " b1=" << newF[1]
                                     << " t1=" << newF[2]
                                     << " t0=" << newF[3]
                                     << " canonSize=" << canonF.size()
                                     << " nBaseInCanon=" << nBaseInCanon
                                     << endl;
                            }
                        }
                        // Use SHARED, once-computed edgeTopSubst (built once at
                        // classification, keyed by (edgeIA, rawTopLabel))
                        // instead of independently recomputing per-cell,
                        // guaranteeing both sides substitute identically.
                        // SAFETY: only base-collapse fallback (OLD, proven-safe
                        // behavior) propagates GLOBALLY. Non-base (shared)
                        // substitution stays LOCAL to this cell only --
                        // CLUSTER_TRACE confirmed global non-base
                        // substitution creates NEW orphaned faces elsewhere.
                        if( !canonContains(newF[2]) )
                        {
                            std::map<std::pair<label,label>,label>::const_iterator sIt =
                                edgeTopSubst.find(std::make_pair(edgeIA, newF[2]));
                            const bool haveShared = (sIt != edgeTopSubst.end());
                            const label subst =
                                haveShared ? sIt->second : newF[1];
                            normalTopSubst[newF[2]] = subst;
                            if( haveShared )
                                normalTopSubstDemands[newF[2]].insert(subst);
                            else
                                globalNormalTopSubst[newF[2]] = subst;
                            Info << "NORMAL_CANON_TRI_TOP_SUBST"
                                 << " bfI=" << bfI
                                 << " edgeIA=" << edgeIA
                                 << " top=" << newF[2]
                                 << " subst=" << subst
                                 << " usedShared=" << (haveShared ? 1 : 0)
                                 << endl;
                        }

                        if( !canonContains(newF[3]) )
                        {
                            std::map<std::pair<label,label>,label>::const_iterator sIt =
                                edgeTopSubst.find(std::make_pair(edgeIA, newF[3]));
                            const bool haveShared = (sIt != edgeTopSubst.end());
                            const label subst =
                                haveShared ? sIt->second : newF[0];
                            normalTopSubst[newF[3]] = subst;
                            if( haveShared )
                                normalTopSubstDemands[newF[3]].insert(subst);
                            else
                                globalNormalTopSubst[newF[3]] = subst;
                            Info << "NORMAL_CANON_TRI_TOP_SUBST"
                                 << " bfI=" << bfI
                                 << " top=" << newF[3]
                                 << " subst=" << subst
                                 << " usedShared=" << (haveShared ? 1 : 0)
                                 << endl;
                        }

                        cellFaces.append(canonF);
                        emittedCanonicalTri = true;
                        normalCellUsedCanonicalTri = true;
                    }
                }

                if( emittedCanonicalTri ) continue;

                cellFaces.append(newF);
                //- check if the face is at the boundary
                //- of the treated partitions
                const label edgeI = faceEdges(bfI, pI);
                if( edgeFaces.sizeOfRow(edgeI) == 2 )
                {
                    label neiFace = edgeFaces(edgeI, 0);
                    if( neiFace == bfI )
                        neiFace = edgeFaces(edgeI, 1);

                    if( !treatPatches[boundaryFacePatches[neiFace]] )
                    {
                        newBoundaryFaces.appendList(newF);
                        newBoundaryOwners.append(cellsToAdd.size() + nOldCells);
                        newBoundaryPatches.append(boundaryFacePatches[neiFace]);
                    }
                }
                else if( edgeFaces.sizeOfRow(edgeI) == 1 )
                {
                    const Map<label>& otherProcPatch = *otherProcPatchPtr;
                    if( !treatPatches[otherProcPatch[edgeI]] )
                    {
                        //- face is a new boundary face
                        newBoundaryFaces.appendList(newF);
                        newBoundaryOwners.append(cellsToAdd.size() + nOldCells);
                        newBoundaryPatches.append(otherProcPatch[edgeI]);
                    }
                }
            }

            # ifdef DEBUGLayer
            Info << "Adding cell " << cellFaces << endl;
            # endif

            } //- end if( !builtReducedCell )

            if( normalCellUsedCanonicalTri && !normalTopSubst.empty() )
            {
                auto substituteAndCleanFace =
                [&](const DynList<label>& inF, DynList<label>& outF) -> bool
                {
                    outF.clear();

                    forAll(inF, ii)
                    {
                        label v = inF[ii];
                        std::map<label,label>::const_iterator sit = normalTopSubst.find(v);
                        if( sit != normalTopSubst.end() ) v = sit->second;

                        bool already = false;
                        forAll(outF, oi)
                        {
                            if( outF[oi] == v )
                            {
                                already = true;
                                break;
                            }
                        }

                        if( !already ) outF.append(v);
                    }

                    return outF.size() >= 3;
                };

                DynList<DynList<label> > cleanedCellFaces;

                forAll(cellFaces, fi)
                {
                    DynList<label> cleanedF;
                    if( substituteAndCleanFace(cellFaces[fi], cleanedF) )
                    {
                        cleanedCellFaces.append(cleanedF);
                    }
                    else
                    {
                        Info << "NORMAL_CANON_TRI_DROPPED_FACE"
                             << " bfI=" << bfI
                             << " fi=" << fi
                             << endl;
                    }
                }

                cellFaces.clear();
                forAll(cleanedCellFaces, fi)
                {
                    cellFaces.append(cleanedCellFaces[fi]);
                }

                if( normalTopBoundaryRow >= 0 && normalTopCellFaceI >= 0 )
                {
                    DynList<label> cleanedTopF;
                    if( substituteAndCleanFace(cleanedCellFaces[normalTopCellFaceI], cleanedTopF) )
                    {
                        newBoundaryFaces.setRow(normalTopBoundaryRow, cleanedTopF);
                        Info << "NORMAL_CANON_TRI_TOPFACE_REWRITE"
                             << " bfI=" << bfI
                             << " row=" << normalTopBoundaryRow
                             << " size=" << cleanedTopF.size()
                             << endl;
                    }
                }
            }

            if( normalCellUsedCanonicalTri )
            {
                // Orient faces of the transition cell consistently.
                // Topological closure only requires every undirected
                // edge to appear twice. Geometric closure additionally
                // requires each shared edge to be traversed in opposite
                // directions by the two incident faces. Canonical side
                // triangles copied from the reduced neighbor can have
                // the correct labels but the wrong orientation for this
                // normal-side transition cell.
                auto reverseDynFace = [](DynList<label>& cf)
                {
                    for(label lo=0, hi=cf.size()-1; lo<hi; ++lo, --hi)
                    {
                        const label tmp = cf[lo];
                        cf[lo] = cf[hi];
                        cf[hi] = tmp;
                    }
                };

                auto sharedEdgeSameDirection =
                [&](const DynList<label>& aF, const DynList<label>& bF,
                    bool& sameDir) -> bool
                {
                    forAll(aF, ai)
                    {
                        const label a0 = aF[ai];
                        const label a1 = aF[(ai+1)%aF.size()];

                        forAll(bF, bi)
                        {
                            const label b0 = bF[bi];
                            const label b1 = bF[(bi+1)%bF.size()];

                            if
                            (
                                (a0 == b0 && a1 == b1)
                             || (a0 == b1 && a1 == b0)
                            )
                            {
                                sameDir = (a0 == b0 && a1 == b1);
                                return true;
                            }
                        }
                    }

                    return false;
                };

                boolList faceOriented(cellFaces.size(), false);
                if( cellFaces.size() ) faceOriented[0] = true;

                bool changedOrient = true;
                label nOrientFlips = 0;
                while( changedOrient )
                {
                    changedOrient = false;

                    forAll(cellFaces, fi)
                    {
                        if( !faceOriented[fi] ) continue;

                        forAll(cellFaces, fj)
                        {
                            if( faceOriented[fj] ) continue;

                            bool sameDir = false;
                            if( sharedEdgeSameDirection(cellFaces[fi], cellFaces[fj], sameDir) )
                            {
                                if( sameDir )
                                {
                                    reverseDynFace(cellFaces[fj]);
                                    ++nOrientFlips;
                                    Info << "NORMAL_CANON_TRI_ORIENT_FLIP"
                                         << " bfI=" << bfI
                                         << " fi=" << fj
                                         << endl;
                                }

                                faceOriented[fj] = true;
                                changedOrient = true;
                            }
                        }
                    }
                }

                label nUnorientedFaces = 0;
                forAll(faceOriented, fi)
                    if( !faceOriented[fi] ) ++nUnorientedFaces;

                if( nOrientFlips || nUnorientedFaces )
                {
                    Info << "NORMAL_CANON_TRI_ORIENT_SUMMARY"
                         << " bfI=" << bfI
                         << " flips=" << nOrientFlips
                         << " unoriented=" << nUnorientedFaces
                         << endl;
                }

                normalCanonTriTopoOK = true;
                std::map<std::pair<label,label>, label> edgeUse;

                forAll(cellFaces, fi)
                {
                    const DynList<label>& cf = cellFaces[fi];
                    if( cf.size() < 3 ) normalCanonTriTopoOK = false;

                    forAll(cf, i)
                    {
                        const label a = cf[i];
                        const label b = cf[(i+1)%cf.size()];
                        if( a == b ) normalCanonTriTopoOK = false;
                        ++edgeUse[std::make_pair(Foam::min(a,b), Foam::max(a,b))];
                    }
                }

                label nBadEdgeUse = 0;
                for
                (
                    std::map<std::pair<label,label>, label>::const_iterator
                        iter = edgeUse.begin();
                    iter != edgeUse.end();
                    ++iter
                )
                {
                    if( iter->second != 2 )
                    {
                        normalCanonTriTopoOK = false;
                        ++nBadEdgeUse;
                    }
                }

                vector areaSum = vector::zero;
                scalar areaMagSum = scalar(0);

                forAll(cellFaces, fi)
                {
                    const DynList<label>& cf = cellFaces[fi];
                    vector fArea = vector::zero;
                    const point& p0f = mesh_.points()[cf[0]];

                    for(label ai=1; ai<cf.size()-1; ++ai)
                    {
                        fArea += (mesh_.points()[cf[ai]] - p0f)
                               ^ (mesh_.points()[cf[ai+1]] - p0f);
                    }

                    areaSum += fArea;
                    areaMagSum += mag(fArea);
                }

                normalCanonTriGeomOK =
                    areaMagSum >= VSMALL
                 && mag(areaSum)/areaMagSum < scalar(1e-4);

                static label nNormalCanonTriCellClosureDiag = 0;
                ++nNormalCanonTriCellClosureDiag;

                Info << "NORMAL_CANON_TRI_CELL_CLOSURE"
                     << " bfI=" << bfI
                     << " bfPatch=" << boundaryFacePatches[bfI]
                     << " nFaces=" << cellFaces.size()
                     << " normalCanonTriTopoOK=" << (normalCanonTriTopoOK ? 1 : 0)
                     << " normalCanonTriGeomOK=" << (normalCanonTriGeomOK ? 1 : 0)
                     << " badEdgeUse=" << nBadEdgeUse
                     << " openness="
                     << (areaMagSum > VSMALL ? mag(areaSum)/areaMagSum : GREAT)
                     << " count=" << nNormalCanonTriCellClosureDiag
                     << endl;

                if( (!normalCanonTriTopoOK || !normalCanonTriGeomOK) && nNormalCanonTriCellClosureDiag <= 20 )
                {
                    forAll(cellFaces, fi)
                    {
                        const DynList<label>& cf = cellFaces[fi];
                        Info << "NORMAL_CANON_TRI_CELL_FACE"
                             << " bfI=" << bfI
                             << " fi=" << fi
                             << " size=" << cf.size()
                             << " labels=(";
                        forAll(cf, pi)
                        {
                            if( pi ) Info << ',';
                            Info << cf[pi];
                        }
                        Info << ')' << endl;
                    }

                    label nPrintedBadEdges = 0;
                    for
                    (
                        std::map<std::pair<label,label>, label>::const_iterator
                            iter = edgeUse.begin();
                        iter != edgeUse.end();
                        ++iter
                    )
                    {
                        if( iter->second != 2 && nPrintedBadEdges < 20 )
                        {
                            ++nPrintedBadEdges;
                            Info << "NORMAL_CANON_TRI_BAD_EDGE"
                                 << " bfI=" << bfI
                                 << " a=" << iter->first.first
                                 << " b=" << iter->first.second
                                 << " use=" << iter->second
                                 << endl;
                        }
                    }
                }
            }

            // ACCEPT/REJECT GATE: topoOK/geomOK were computed above but
            // never actually gated cellsToAdd.appendGraph() -- meaning
            // internally non-manifold cells (badEdgeUse>0, e.g. an edge
            // used 3 times instead of exactly 2) were queued anyway,
            // creating orphaned faces downstream. Reject here instead.
            if( normalCanonTriTopoOK && normalCanonTriGeomOK )
            {
                cellsToAdd.appendGraph(cellFaces);
            }
            else
            {
                Info << "NORMAL_CANON_TRI_CELL_REJECTED bfI=" << bfI
                     << " -- falling back to raw boundary face"
                     << endl;
                newBoundaryFaces.appendList(bFaces[bfI]);
                newBoundaryOwners.append(faceOwners[bfI]);
                newBoundaryPatches.append(boundaryFacePatches[bfI]);
            }
        }
        else
        {
            # ifdef DEBUGLayer
            Info << "Storing original boundary face "
                << bfI << " into patch " << boundaryFacePatches[bfI] << endl;
            # endif

            newBoundaryFaces.appendList(bFaces[bfI]);
            newBoundaryOwners.append(faceOwners[bfI]);
            newBoundaryPatches.append(boundaryFacePatches[bfI]);
        }
    }

    // ---- BLCELLGEOM_AUDIT summary (report-only) ----
    Info << "BLCELLGEOM_AUDIT (base->top distance, both resolvers, at cell build)"
         << endl;

    forAll(clqTot, p)
    {
        if( clqTot[p] == 0 )
            continue;

        const word pnm =
            (p >= 0 && p < label(patchNames_.size()))
            ? patchNames_[p]
            : word("?");

        const scalar pMean =
            clqPlainNvalid[p] > 0
            ? clqPlainSum[p] / clqPlainNvalid[p]
            : scalar(0);

        const scalar qMean =
            clqPatchNvalid[p] > 0
            ? clqPatchSum[p] / clqPatchNvalid[p]
            : scalar(0);

        const scalar pMin =
            clqPlainNvalid[p] > 0
            ? clqPlainMin[p]
            : scalar(-1);

        const scalar qMin =
            clqPatchNvalid[p] > 0
            ? clqPatchMin[p]
            : scalar(-1);

        Info << "  " << pnm
             << ": faces=" << clqTot[p]

             << " | PLAIN"
             << " inv=" << clqPlainInvalid[p]
             << " same=" << clqPlainSame[p]
             << " near=" << clqPlainNear[p]
             << " fullColl=" << clqPlainFullColl[p]
             << " partial=" << clqPlainPartial[p]
             << " healthy=" << clqPlainHealthy[p]
             << " minD=" << pMin
             << " meanD=" << pMean

             << " | PATCH"
             << " inv=" << clqPatchInvalid[p]
             << " same=" << clqPatchSame[p]
             << " near=" << clqPatchNear[p]
             << " fullColl=" << clqPatchFullColl[p]
             << " partial=" << clqPatchPartial[p]
             << " healthy=" << clqPatchHealthy[p]
             << " minD=" << qMin
             << " meanD=" << qMean

             << " | resolverDiff=" << clqDifferent[p]
             << " capTouched=" << clqCapTouched[p]
             << endl;
    }
    // ---- end BLCELLGEOM_AUDIT summary ----

    //- data for parallel execution
    boolList procPoint;
    LongList<DynList<label, 4> > pointProcFaces;
    LongList<labelPair> faceAtPatches;
    if( Pstream::parRun() )
    {
        procPoint.setSize(nPoints_);
        procPoint = false;

        const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();
        const labelList& bPoints = mse.boundaryPoints();

        for
        (
            Map<label>::const_iterator iter=globalToLocal.begin();
            iter!=globalToLocal.end();
            ++iter
        )
        {
            const label bpI = iter();
            procPoint[bPoints[bpI]] = true;
        }
    }

    //- create cells at edges
    forAll(edgeFaces, edgeI)
    {
        //- do not consider edges with no faces attached to it
        if( edgeFaces.sizeOfRow(edgeI) == 0 )
            continue;

        //- cells are generated at the processor with the lowest label
        if(
            (edgeFaces.sizeOfRow(edgeI) == 1) &&
            (otherFaceProcPtr->operator[](edgeI) < Pstream::myProcNo())
        )
            continue;

        //- check if the edge is a feature edge
        const label patchI = boundaryFacePatches[edgeFaces(edgeI, 0)];

        label patchJ(-1);
        if( Pstream::parRun() && otherProcPatchPtr && otherProcPatchPtr->found(edgeI) )
        {
            patchJ = otherProcPatchPtr->operator[](edgeI);
        }
        else if( edgeFaces.sizeOfRow(edgeI) >= 2 )
        {
            patchJ = boundaryFacePatches[edgeFaces(edgeI, 1)];
        }
        else
        {
            continue;  // guard: single-face edge with no parallel info
        }

        if( patchI == patchJ )
            continue;

        //- check if the faces attached to the edge have different keys
        const label pKeyI = patchKey_[patchI];
        const label pKeyJ = patchKey_[patchJ];

        if( pKeyI < 0 || pKeyJ < 0 )
        {
            continue;
            FatalErrorIn
            (
                "void boundaryLayers::createLayerCells(const labelList&)"
            ) << "Patch key is negative at concave edge" << abort(FatalError);
        }

        if( pKeyI == pKeyJ )
            continue;

        const edge& e = edges[edgeI];
        if( otherVrts_.find(e.start()) == otherVrts_.end() )
            continue;
        if( otherVrts_.find(e.end()) == otherVrts_.end() )
            continue;

        //- generate faces of the bnd layer cell
        FixedList<FixedList<label, 4>, 6> cellFaces;
        createNewCellFromEdge(e, pKeyI, pKeyJ, cellFaces);

        //- store boundary faces
        newBoundaryFaces.appendList(cellFaces[1]);
        newBoundaryOwners.append(nOldCells+cellsToAdd.size());
        newBoundaryPatches.append(patchJ);

        newBoundaryFaces.appendList(cellFaces[3]);
        newBoundaryOwners.append(nOldCells+cellsToAdd.size());
        newBoundaryPatches.append(patchI);

        //- check if face 5 is a boundary face or at an inter-processor boundary
        const label bps = bp[e.start()];
        if( bps < 0 ) continue;  // guard: invalid boundary point
        label unusedPatch(-1);
        forAllRow(pointPatches, bps, i)
        {
            const label ptchI = pointPatches(bps, i);

            if( ptchI == patchI )
                continue;
            if( ptchI == patchJ )
                continue;
            if( unusedPatch != -1 )
            {
                unusedPatch = -1;
                break;
            }

            unusedPatch = ptchI;
        }

        if( unusedPatch != -1 && treatedPatch_[unusedPatch] )
        {
            //- add a face in the empty patch in case of 2D layer generation
            newBoundaryFaces.appendList(cellFaces[5]);
            newBoundaryOwners.append(nOldCells+cellsToAdd.size());
            newBoundaryPatches.append(unusedPatch);
        }
        else if( Pstream::parRun() && procPoint[e.start()] )
        {
            //- add a face at inter-pocessor boundary
            pointProcFaces.append(cellFaces[5]);
            faceAtPatches.append(labelPair(patchI, patchJ));
        }

        //- check if face 4 is a boundary face or at an inter-processor boundary
        const label bpe = bp[e.end()];
        if( bpe < 0 ) continue;  // guard: invalid boundary point
        unusedPatch = -1;
        forAllRow(pointPatches, bpe, i)
        {
            const label ptchI = pointPatches(bpe, i);

            if( ptchI == patchI )
                continue;
            if( ptchI == patchJ )
                continue;
            if( unusedPatch != -1 )
            {
                unusedPatch = -1;
                break;
            }

            unusedPatch = ptchI;
        }

        if( unusedPatch != -1 && treatedPatch_[unusedPatch] )
        {
            //- add a face in the empty patch in case of 2D layer generation
            newBoundaryFaces.appendList(cellFaces[4]);
            newBoundaryOwners.append(nOldCells+cellsToAdd.size());
            newBoundaryPatches.append(unusedPatch);
        }
        else if( Pstream::parRun() && procPoint[e.end()] )
        {
            //- add a face at inter-pocessor boundary
            pointProcFaces.append(cellFaces[4]);
            faceAtPatches.append(labelPair(patchI, patchJ));
        }

        # ifdef DEBUGLayer
        Info << "Adding new cell at edge " << cellFaces << endl;
        # endif

        //- append cell to the queue
        cellsToAdd.appendGraph(cellFaces);
    }

    //- create cells for corner nodes
    typedef std::map<std::pair<label, label>, label> mPairToLabelType;
    typedef std::map<label, mPairToLabelType> mPointsType;
    typedef std::map<label, DynList<label, 3> > ppType;

    ppType nodePatches;
    labelHashSet parPoint;

    if( Pstream::parRun() )
    {
        const labelList& bPoints = mse.boundaryPoints();
        const VRWGraph& pProcs = mse.bpAtProcs();
        const labelList& globalPointLabel = mse.globalBoundaryPointLabel();
        const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();

        std::map<label, labelLongList> facesToSend;

        typedef std::map<label, DynList<DynList<label, 8>, 8> > ppfType;

        ppfType parPointFaces;
        ppType parPointPatches;

        forAllConstIter(mPointsType, otherVrts_, iter)
        {
            //- skip points on feature edges
            if( iter->second.size() == 2 )
                continue;

            const label bpI = bp[iter->first];

            if( pProcs.sizeOfRow(bpI) != 0 )
            {
                parPoint.insert(iter->first);

                //- point is at a parallel boundary
                label pMin = pProcs(bpI, 0);
                forAllRow(pProcs, bpI, i)
                {
                    const label prI = pProcs(bpI, i);

                    if( facesToSend.find(prI) == facesToSend.end() )
                        facesToSend.insert
                        (
                            std::make_pair(prI, labelLongList())
                        );

                    if( prI < pMin )
                        pMin = prI;
                }

                if( Pstream::myProcNo() == pMin )
                {
                    DynList<label, 3>& pPatches = parPointPatches[bpI];
                    pPatches.setSize(pointFaces.sizeOfRow(bpI));

                    DynList<DynList<label, 8>, 8>& pFaces = parPointFaces[bpI];
                    pFaces.setSize(pPatches.size());

                    forAllRow(pointFaces, bpI, pfI)
                    {
                        const label bfI = pointFaces(bpI, pfI);
                        const face& bf = bFaces[bfI];

                        pPatches[pfI] = boundaryFacePatches[bfI];

                        DynList<label, 8>& bfCopy = pFaces[pfI];
                        bfCopy.setSize(bf.size());
                        forAll(bf, pI)
                            bfCopy[pI] = globalPointLabel[bp[bf[pI]]];
                    }

                    continue;
                }

                labelLongList& stp = facesToSend[pMin];

                //- send the data to the processor with the lowest label
                //- data is flatenned as follows
                //- 1. the number of faces and global point label
                //- 2. number of points in the face
                //- 3. patch label
                //- 4. global labels of face points
                stp.append(globalPointLabel[bpI]);
                stp.append(pointFaces.sizeOfRow(bpI));
                forAllRow(pointFaces, bpI, pfI)
                {
                    const label bfI = pointFaces(bpI, pfI);
                    const face& bf = bFaces[bfI];

                    stp.append(bf.size());
                    stp.append(boundaryFacePatches[bfI]);
                    forAll(bf, pI)
                        stp.append(globalPointLabel[bp[bf[pI]]]);
                }
            }
        }

        //- exchange data with other processors
        labelLongList receivedData;
        help::exchangeMap(facesToSend, receivedData);

        label counter(0);
        while( counter < receivedData.size() )
        {
            const label gpI_recv = receivedData[counter++];
            const label nFaces = receivedData[counter++];
            if( !globalToLocal.found(gpI_recv) )
            {
                // Skip all face data for this point.
                // Each face packet: face_size + patch_label + face_size points
                for(label fI=0;fI<nFaces;++fI)
                {
                    const label fSize = receivedData[counter++];
                    ++counter; // patch label
                    counter += fSize; // face points
                }
                continue;
            }
            const label bpI = globalToLocal[gpI_recv];
            for(label fI=0;fI<nFaces;++fI)
            {
                DynList<label, 8> f(receivedData[counter++]);
                parPointPatches[bpI].append(receivedData[counter++]);
                forAll(f, pI)
                    f[pI] = receivedData[counter++];
                parPointFaces[bpI].append(f);
            }
        }

        //- sort faces sharing corners at the parallel boundaries
        forAllIter(ppfType, parPointFaces, iter)
        {
            DynList<DynList<label, 8>, 8>& pFaces = iter->second;
            DynList<label, 3>& fPatches = parPointPatches[iter->first];
            const label gpI = globalPointLabel[iter->first];

            for(label i=0;i<pFaces.size();++i)
            {
                const DynList<label, 8>& bf = pFaces[i];
                const label pos = bf.containsAtPosition(gpI);
                const edge e(bf[pos], bf[bf.fcIndex(pos)]);

                for(label j=i+1;j<pFaces.size();++j)
                {
                    const DynList<label, 8>& obf = pFaces[j];
                    if( obf.contains(e.start()) && obf.contains(e.end()) )
                    {
                        DynList<label, 8> add;
                        add = pFaces[i+1];
                        pFaces[i+1] = pFaces[j];
                        pFaces[j] = add;

                        const label pAdd = fPatches[i+1];
                        fPatches[i+1] = fPatches[j];
                        fPatches[j] = pAdd;
                        break;
                    }
                }
            }

            DynList<label, 3> patchIDs;
            forAll(fPatches, fpI)
                patchIDs.appendIfNotIn(fPatches[fpI]);

            nodePatches.insert(std::make_pair(bPoints[iter->first], patchIDs));
        }
    }

    //- sort out point which are not at inter-processor boundaries
    forAllConstIter(mPointsType, otherVrts_, iter)
    {
        if( iter->second.size() == 2 )
            continue;

        if( parPoint.found(iter->first) )
            continue;

        const label bpI = bp[iter->first];

        //- ensure correct orientation
        DynList<label> pFaces(pointFaces.sizeOfRow(bpI));
        forAll(pFaces, fI)
            pFaces[fI] = pointFaces(bpI, fI);

        for(label i=0;i<pFaces.size();++i)
        {
            const face& bf = bFaces[pFaces[i]];
            const edge e = bf.faceEdge(bf.which(iter->first));

            for(label j=i+1;j<pFaces.size();++j)
            {
                const face& obf = bFaces[pFaces[j]];
                if(
                    (obf.which(e.start()) >= 0) &&
                    (obf.which(e.end()) >= 0)
                )
                {
                    const label add = pFaces[i+1];
                    pFaces[i+1] = pFaces[j];
                    pFaces[j] = add;
                    break;
                }
            }
        }

        DynList<label, 3> patchIDs;
        forAll(pFaces, patchI)
        {
            patchIDs.appendIfNotIn(boundaryFacePatches[pFaces[patchI]]);
        }

        nodePatches.insert(std::make_pair(iter->first, patchIDs));
    }

    //- create layer cells for corner nodes
    forAllIter(ppType, nodePatches, iter)
    {
        const DynList<label, 3>& patchIDs = iter->second;
        DynList<label, 3> pKeys;
        forAll(patchIDs, patchI)
        {
            const label pKey = patchKey_[patchIDs[patchI]];

            if( pKey < 0 )
                continue;

            pKeys.appendIfNotIn(pKey);
        }

        if( pKeys.size() != 3 )
            continue;

        # ifdef DEBUGLayer
        Pout << "Creating corner cell at point " << iter->first << endl;
        # endif

        FixedList<FixedList<label, 4>, 6> cellFaces;
        createNewCellFromNode(iter->first, pKeys, cellFaces);

        //- store boundary faces
        newBoundaryFaces.appendList(cellFaces[1]);
        newBoundaryOwners.append(nOldCells+cellsToAdd.size());
        newBoundaryPatches.append(patchIDs[0]);

        newBoundaryFaces.appendList(cellFaces[3]);
        newBoundaryOwners.append(nOldCells+cellsToAdd.size());
        newBoundaryPatches.append(patchIDs[1]);

        newBoundaryFaces.appendList(cellFaces[5]);
        newBoundaryOwners.append(nOldCells+cellsToAdd.size());
        newBoundaryPatches.append(patchIDs[2]);

        if( Pstream::parRun() )
        {
            if( procPoint[iter->first] )
            {
                pointProcFaces.append(cellFaces[0]);
                faceAtPatches.append(labelPair(patchIDs[1], patchIDs[2]));

                pointProcFaces.append(cellFaces[2]);
                faceAtPatches.append(labelPair(patchIDs[0], patchIDs[2]));

                pointProcFaces.append(cellFaces[4]);
                faceAtPatches.append(labelPair(patchIDs[0], patchIDs[1]));
            }
        }

        # ifdef DEBUGLayer
        Info << "Adding corner cell " << cellFaces << endl;
        # endif

        //- append cell to the queue
        cellsToAdd.appendGraph(cellFaces);
    }

    if( Pstream::parRun() )
    {
        //- create faces at parallel boundaries created from
        //- points at parallel boundaries
        createNewFacesFromPointsParallel
        (
            pointProcFaces,
            faceAtPatches
        );
    }

    // CLUSTER_TRACE: forensic tracer for hardcoded label set,
    // temporary diagnostic only, not a fix.
    auto debugHitLabel = [&](const label v) -> bool
    {
        return v == 1706262 || v == 1706257 || v == 1706265
            || v == 1376401 || v == 1376402 || v == 1706263;
    };

    auto debugContainsLabelDyn = [&](const DynList<label>& ff) -> bool
    {
        forAll(ff, i)
        {
            if( debugHitLabel(ff[i]) ) return true;
        }
        return false;
    };

    auto debugPrintFaceDyn =
    [&](const char* tag, const label cellI, const label faceI, const DynList<label>& ff)
    {
        if( debugContainsLabelDyn(ff) )
        {
            Info << "CLUSTER_TRACE " << tag
                 << " cell=" << cellI
                 << " face=" << faceI
                 << " size=" << ff.size()
                 << " labels=(";
            forAll(ff, i)
            {
                Info << ff[i];
                if( i+1 < ff.size() ) Info << ",";
            }
            Info << ")" << endl;
        }
    };

    auto debugPrintFaceKey =
    [&](const char* tag, const label q, const label b, const label e, const std::vector<label>& k)
    {
        bool hit = false;
        for(size_t i=0; i<k.size(); ++i)
        {
            if( debugHitLabel(k[i]) ) { hit = true; break; }
        }
        if( hit )
        {
            Info << "CLUSTER_TRACE " << tag
                 << " q=" << q
                 << " b=" << b
                 << " e=" << e
                 << " size=" << label(k.size())
                 << " labels=(";
            for(size_t i=0; i<k.size(); ++i)
            {
                Info << k[i];
                if( i+1 < k.size() ) Info << ",";
            }
            Info << ")" << endl;
        }
    };

    // Resolve substitution demands: a raw vertex is only safe to
    // substitute GLOBALLY if every demand for it agreed on the SAME
    // target. Conflicting demands are dropped for global purposes;
    // the per-cell normalTopSubst still applies locally regardless.
    {
        label nDemandUnanimous = 0;
        label nDemandConflict = 0;
        for
        (
            std::map<label, std::set<label>>::const_iterator dIt =
                normalTopSubstDemands.begin();
            dIt != normalTopSubstDemands.end();
            ++dIt
        )
        {
            const label rawV = dIt->first;
            const std::set<label>& targets = dIt->second;
            if( targets.size() == 1 )
            {
                ++nDemandUnanimous;
                // Do NOT promote shared/non-base canonical-triangle
                // substitutions into the global row-rewrite map.
                // They were already applied locally to the canonical
                // transition face. Promoting them globally can rewrite
                // unrelated queued/boundary rows and create orphan
                // interface topology under different construction order.
            }
            else
            {
                ++nDemandConflict;
                Info << "NORMALTOPSUBST_DEMAND_CONFLICT rawV=" << rawV
                     << " nTargets=" << label(targets.size())
                     << " targets=(";
                label ti = 0;
                for
                (
                    std::set<label>::const_iterator tIt = targets.begin();
                    tIt != targets.end();
                    ++tIt, ++ti
                )
                {
                    if( ti ) Info << ",";
                    Info << *tIt;
                }
                Info << ") -- rejected, no global substitution applied" << endl;
            }
        }
        Info << "NORMALTOPSUBST_DEMAND_RESOLUTION"
             << " unanimous=" << nDemandUnanimous
             << " conflict=" << nDemandConflict
             << endl;
    }

    // FINALCELL_CLOSURE_AUDIT
    // Diagnostic only. Audit every queued cell shell before and after
    // GLOBAL_NORMAL_TOP_SUBST. INTERFACEAUDIT validates face pairing,
    // but cannot detect a face that is absent from a cell altogether.
    auto auditQueuedCellClosure =
    [&](const char* stage)
    {
        label nBadCells = 0;
        label nBadEdgesTotal = 0;
        label nMalformedFacesTotal = 0;

        for(label cellI=0; cellI<cellsToAdd.size(); ++cellI)
        {
            std::map<std::pair<label,label>, label> edgeUse;

            bool malformedFace = false;
            label nMalformedFaces = 0;

            for
            (
                label faceI=0;
                faceI<cellsToAdd.sizeOfGraph(cellI);
                ++faceI
            )
            {
                const label nVerts =
                    cellsToAdd.sizeOfRow(cellI, faceI);

                if( nVerts < 3 )
                {
                    malformedFace = true;
                    ++nMalformedFaces;
                    continue;
                }

                for(label pI=0; pI<nVerts; ++pI)
                {
                    const label a =
                        cellsToAdd(cellI, faceI, pI);

                    const label b =
                        cellsToAdd
                        (
                            cellI,
                            faceI,
                            (pI+1)%nVerts
                        );

                    if( a == b )
                    {
                        malformedFace = true;
                    }

                    ++edgeUse
                    [
                        std::make_pair
                        (
                            Foam::min(a,b),
                            Foam::max(a,b)
                        )
                    ];
                }
            }

            label nBadEdges = 0;

            for
            (
                std::map<std::pair<label,label>, label>
                    ::const_iterator iter=edgeUse.begin();
                iter!=edgeUse.end();
                ++iter
            )
            {
                if( iter->second != 2 )
                {
                    ++nBadEdges;
                }
            }

            if( malformedFace || nBadEdges )
            {
                ++nBadCells;
                nBadEdgesTotal += nBadEdges;
                nMalformedFacesTotal += nMalformedFaces;

                if( nBadCells <= 100 )
                {
                    Info
                        << "FINALCELL_CLOSURE_BAD"
                        << " stage=" << stage
                        << " queuedCell=" << cellI
                        << " predictedCell="
                        << (nOldCells + cellI)
                        << " nFaces="
                        << cellsToAdd.sizeOfGraph(cellI)
                        << " malformedFaces="
                        << nMalformedFaces
                        << " badEdges="
                        << nBadEdges
                        << endl;

                    label nPrintedEdges = 0;

                    for
                    (
                        std::map<std::pair<label,label>, label>
                            ::const_iterator iter=edgeUse.begin();
                        iter!=edgeUse.end();
                        ++iter
                    )
                    {
                        if
                        (
                            iter->second != 2
                         && nPrintedEdges < 20
                        )
                        {
                            ++nPrintedEdges;

                            Info
                                << "FINALCELL_CLOSURE_EDGE"
                                << " stage=" << stage
                                << " queuedCell=" << cellI
                                << " predictedCell="
                                << (nOldCells + cellI)
                                << " a=" << iter->first.first
                                << " b=" << iter->first.second
                                << " use=" << iter->second
                                << endl;
                        }
                    }
                }
            }
        }

        Info
            << "FINALCELL_CLOSURE_AUDIT"
            << " stage=" << stage
            << " queuedCells=" << cellsToAdd.size()
            << " badCells=" << nBadCells
            << " badEdgesTotal=" << nBadEdgesTotal
            << " malformedFacesTotal="
            << nMalformedFacesTotal
            << endl;
    };

    auditQueuedCellClosure("preSubst");

    //- create mesh modifier
    if( !globalNormalTopSubst.empty() )
    {
        Info << "GLOBAL_NORMAL_TOP_SUBST_BEGIN"
             << " nRules=" << globalNormalTopSubst.size()
             << endl;

        auto substituteAndCleanQueuedFace =
        [&](const VRWGraph& g, const label rowI, DynList<label>& outF) -> bool
        {
            outF.clear();

            for(label pI=0; pI<g.sizeOfRow(rowI); ++pI)
            {
                label v = g(rowI, pI);
                std::map<label,label>::const_iterator sit = globalNormalTopSubst.find(v);
                if( sit != globalNormalTopSubst.end() ) v = sit->second;

                bool already = false;
                forAll(outF, oi)
                {
                    if( outF[oi] == v )
                    {
                        already = true;
                        break;
                    }
                }

                if( !already ) outF.append(v);
            }

            return outF.size() >= 3;
        };

        label nCellRowsChanged = 0;
        label nCellRowsDropped = 0;

        VRWGraphList cleanedCellsToAdd;

        for(label cellI=0; cellI<cellsToAdd.size(); ++cellI)
        {
            DynList<DynList<label> > cleanedCellFaces;

            for(label faceI=0; faceI<cellsToAdd.sizeOfGraph(cellI); ++faceI)
            {
                {
                    DynList<label> rawF;
                    for(label pI2=0; pI2<cellsToAdd.sizeOfRow(cellI, faceI); ++pI2)
                        rawF.append(cellsToAdd(cellI, faceI, pI2));
                    debugPrintFaceDyn("queued_before", cellI, faceI, rawF);
                }
                DynList<label> cleanedF;

                for(label pI=0; pI<cellsToAdd.sizeOfRow(cellI, faceI); ++pI)
                {
                    label v = cellsToAdd(cellI, faceI, pI);
                    std::map<label,label>::const_iterator sit = globalNormalTopSubst.find(v);
                    if( sit != globalNormalTopSubst.end() ) v = sit->second;

                    bool already = false;
                    forAll(cleanedF, oi)
                    {
                        if( cleanedF[oi] == v )
                        {
                            already = true;
                            break;
                        }
                    }

                    if( !already ) cleanedF.append(v);
                }
                debugPrintFaceDyn("queued_after", cellI, faceI, cleanedF);

                if( cleanedF.size() >= 3 )
                {
                    bool changed = cleanedF.size() != cellsToAdd.sizeOfRow(cellI, faceI);
                    if( !changed )
                    {
                        for(label pI=0; pI<cleanedF.size(); ++pI)
                        {
                            if( cleanedF[pI] != cellsToAdd(cellI, faceI, pI) )
                            {
                                changed = true;
                                break;
                            }
                        }
                    }

                    if( changed ) ++nCellRowsChanged;
                    cleanedCellFaces.append(cleanedF);
                }
                else
                {
                    ++nCellRowsDropped;
                    Info << "GLOBAL_NORMAL_TOP_SUBST_BAD_CELL_ROW"
                         << " cell=" << cellI
                         << " face=" << faceI
                         << endl;
                }
            }

            cleanedCellsToAdd.appendGraph(cleanedCellFaces);
        }

        cellsToAdd = cleanedCellsToAdd;

        label nBndRowsChanged = 0;
        label nBndRowsDropped = 0;

        for(label rowI=0; rowI<newBoundaryFaces.size(); ++rowI)
        {
            {
                DynList<label> rawBndF;
                for(label pI2=0; pI2<newBoundaryFaces.sizeOfRow(rowI); ++pI2)
                    rawBndF.append(newBoundaryFaces(rowI, pI2));
                debugPrintFaceDyn("boundary_before", -1, rowI, rawBndF);
            }
            DynList<label> cleanedF;
            if( substituteAndCleanQueuedFace(newBoundaryFaces, rowI, cleanedF) )
            {
                debugPrintFaceDyn("boundary_after", -1, rowI, cleanedF);
                bool changed = cleanedF.size() != newBoundaryFaces.sizeOfRow(rowI);
                if( !changed )
                {
                    for(label pI=0; pI<cleanedF.size(); ++pI)
                    {
                        if( cleanedF[pI] != newBoundaryFaces(rowI, pI) )
                        {
                            changed = true;
                            break;
                        }
                    }
                }

                if( changed ) ++nBndRowsChanged;
                newBoundaryFaces.setRow(rowI, cleanedF);
            }
            else
            {
                ++nBndRowsDropped;
                Info << "GLOBAL_NORMAL_TOP_SUBST_BAD_BND_ROW"
                     << " row=" << rowI
                     << endl;
            }
        }

        Info << "GLOBAL_NORMAL_TOP_SUBST_DONE"
             << " nRules=" << globalNormalTopSubst.size()
             << " cellRowsChanged=" << nCellRowsChanged
             << " cellRowsDropped=" << nCellRowsDropped
             << " bndRowsChanged=" << nBndRowsChanged
             << " bndRowsDropped=" << nBndRowsDropped
             << endl;
    }

    auditQueuedCellClosure("postSubst");

    polyMeshGenModifier meshModifier(mesh_);

    //- Diagnostic only: audit the addCells() face-pairing contract.
    //  A queued face is valid if:
    //   (a) it appears twice across cellsToAdd: internal face between new cells;
    //   (b) it appears once in cellsToAdd and once in newBoundaryFaces: new boundary;
    //   (c) it appears once in cellsToAdd and already exists in mesh_.faces():
    //       old surface/base face being converted into an internal face.
    //  Anything else is a true suspicious orphan/mismatch candidate.
    {
        typedef std::vector<label> FaceKey;

        struct FaceKeyHash
        {
            std::size_t operator()(const FaceKey& k) const
            {
                std::size_t h = 1469598103934665603ULL;
                for(std::size_t i=0; i<k.size(); ++i)
                {
                    h ^= std::size_t(k[i]);
                    h *= 1099511628211ULL;
                }
                return h;
            }
        };

        // Canonical polygon key preserving cyclic edge connectivity.
        //
        // Equivalent:
        //   (a b c d)
        //   (b c d a)
        //   (d c b a)
        //
        // Not equivalent:
        //   (a c b d)
        //
        // The previous sorted-label key discarded edge adjacency and
        // could therefore declare topologically different polygons equal.
        auto canonicalFaceKey =
        [](FaceKey k) -> FaceKey
        {
            if( k.size() < 2 )
                return k;

            const label n = label(k.size());

            FaceKey best;
            bool haveBest = false;

            // Test all cyclic rotations in both directions.
            for(label reverse=0; reverse<2; ++reverse)
            {
                for(label startI=0; startI<n; ++startI)
                {
                    FaceKey candidate;
                    candidate.reserve(k.size());

                    for(label offset=0; offset<n; ++offset)
                    {
                        label idx;

                        if( !reverse )
                        {
                            idx = (startI + offset) % n;
                        }
                        else
                        {
                            idx = startI - offset;

                            while( idx < 0 )
                                idx += n;

                            idx %= n;
                        }

                        candidate.push_back(k[idx]);
                    }

                    if( !haveBest || candidate < best )
                    {
                        best = candidate;
                        haveBest = true;
                    }
                }
            }

            return best;
        };

        std::unordered_map<FaceKey, label, FaceKeyHash> queuedUse;
        std::unordered_map<FaceKey, label, FaceKeyHash> boundaryUse;
        std::unordered_map<FaceKey, label, FaceKeyHash> existingUse;

        forAll(cellsToAdd, cI)
        {
            for(label fI=0; fI<cellsToAdd.sizeOfGraph(cI); ++fI)
            {
                FaceKey k;
                k.reserve(cellsToAdd.sizeOfRow(cI, fI));

                for(label pI=0; pI<cellsToAdd.sizeOfRow(cI, fI); ++pI)
                {
                    k.push_back(cellsToAdd(cI, fI, pI));
                }

                k = canonicalFaceKey(k);
                ++queuedUse[k];
            }
        }

        for(label fI=0; fI<newBoundaryFaces.size(); ++fI)
        {
            FaceKey k;
            k.reserve(newBoundaryFaces.sizeOfRow(fI));

            for(label pI=0; pI<newBoundaryFaces.sizeOfRow(fI); ++pI)
            {
                k.push_back(newBoundaryFaces(fI, pI));
            }

            k = canonicalFaceKey(k);
            ++boundaryUse[k];
        }

        const faceListPMG& oldFaces = mesh_.faces();

        forAll(oldFaces, fI)
        {
            const face& f = oldFaces[fI];

            FaceKey k;
            k.reserve(f.size());

            forAll(f, pI)
            {
                k.push_back(f[pI]);
            }

            k = canonicalFaceKey(k);
            ++existingUse[k];
        }

        label nAuditBad = 0;
        label nAuditInternalOK = 0;
        label nAuditBoundaryOK = 0;
        label nAuditExistingOK = 0;

        for
        (
            std::unordered_map<FaceKey, label, FaceKeyHash>::const_iterator
                iter = queuedUse.begin();
            iter != queuedUse.end();
            ++iter
        )
        {
            const FaceKey& k = iter->first;
            const label q = iter->second;

            std::unordered_map<FaceKey, label, FaceKeyHash>::const_iterator
                bIter = boundaryUse.find(k);

            std::unordered_map<FaceKey, label, FaceKeyHash>::const_iterator
                eIter = existingUse.find(k);

            const label b =
                (bIter == boundaryUse.end()) ? label(0) : bIter->second;

            const label e =
                (eIter == existingUse.end()) ? label(0) : eIter->second;

            const bool okInternal = (q == 2 && b == 0 && e == 0);
            const bool okBoundary = (q == 1 && b == 1 && e == 0);
            const bool okExisting = (q == 1 && b == 0 && e > 0);

            if( okInternal )
            {
                ++nAuditInternalOK;
            }
            else if( okBoundary )
            {
                ++nAuditBoundaryOK;
            }
            else if( okExisting )
            {
                ++nAuditExistingOK;
            }
            else
            {
                ++nAuditBad;

                if( nAuditBad <= 100 )
                {
                    debugPrintFaceKey("audit_orphan_key", q, b, e, k);
                    Info << "INTERFACEAUDIT_ORPHAN"
                         << " q=" << q
                         << " b=" << b
                         << " e=" << e
                         << " size=" << label(k.size())
                         << " labels=(";

                    for(label i=0; i<label(k.size()); ++i)
                    {
                        if( i ) Info << ',';
                        Info << k[i];
                    }

                    Info << ')' << endl;
                }
            }
        }

        Info << "INTERFACEAUDIT summary"
             << " queuedKeys=" << label(queuedUse.size())
             << " boundaryKeys=" << label(boundaryUse.size())
             << " existingKeys=" << label(existingUse.size())
             << " internalOK=" << nAuditInternalOK
             << " boundaryOK=" << nAuditBoundaryOK
             << " existingOK=" << nAuditExistingOK
             << " bad=" << nAuditBad
             << endl;

        if( nAuditBad != 0 )
        {
            FatalErrorInFunction
                << "INTERFACEAUDIT_FATAL_BAD_TOPOLOGY"
                << " bad=" << nAuditBad
                << " -- refusing meshModifier.addCells(cellsToAdd) because queued transition topology is not pair-consistent."
                << exit(FatalError);
        }

    }

    // POSTCOMMIT_ZIPUP_AUDIT
    //
    // Diagnostic only. Whole-mesh topological cell-closure audit.
    // A valid closed polyhedral shell requires every undirected edge
    // of every cell to occur exactly twice across that cell's faces.
    //
    // This intentionally duplicates the essential checkCellsZipUp()
    // invariant locally so this forensic diagnostic has no dependency
    // on polyMeshGenChecks linkage/namespace details.
    auto auditMeshZipUp =
    [&](const char* stage)
    {
        const cellListPMG& auditCells = mesh_.cells();
        const faceListPMG& auditFaces = mesh_.faces();

        labelLongList badCellIds;

        label nBadEdgesTotal = 0;
        label nBadFaceRefs = 0;
        label nDegenerateFaces = 0;

        forAll(auditCells, cellI)
        {
            const cell& c = auditCells[cellI];

            std::map<std::pair<label,label>, label> edgeUse;

            bool cellBad = false;
            label nBadEdgesCell = 0;

            forAll(c, cfI)
            {
                const label faceI = c[cfI];

                if( faceI < 0 || faceI >= auditFaces.size() )
                {
                    cellBad = true;
                    ++nBadFaceRefs;
                    continue;
                }

                const face& f = auditFaces[faceI];

                if( f.size() < 3 )
                {
                    cellBad = true;
                    ++nDegenerateFaces;
                    continue;
                }

                forAll(f, pI)
                {
                    const label a = f[pI];
                    const label b = f[(pI+1)%f.size()];

                    if( a == b )
                    {
                        cellBad = true;
                    }

                    ++edgeUse
                    [
                        std::make_pair
                        (
                            Foam::min(a,b),
                            Foam::max(a,b)
                        )
                    ];
                }
            }

            for
            (
                std::map<std::pair<label,label>, label>
                    ::const_iterator iter=edgeUse.begin();
                iter!=edgeUse.end();
                ++iter
            )
            {
                if( iter->second != 2 )
                {
                    cellBad = true;
                    ++nBadEdgesCell;
                }
            }

            if( cellBad )
            {
                badCellIds.append(cellI);
                nBadEdgesTotal += nBadEdgesCell;
            }
        }

        Info << "POSTCOMMIT_ZIPUP"
             << " stage=" << stage
             << " cells=" << auditCells.size()
             << " badCells=" << badCellIds.size()
             << " badEdgesTotal=" << nBadEdgesTotal
             << " badFaceRefs=" << nBadFaceRefs
             << " degenerateFaces=" << nDegenerateFaces
             << endl;

        if( badCellIds.size() )
        {
            Info << "POSTCOMMIT_ZIPUP_IDS"
                 << " stage=" << stage
                 << " ids=(";

            const label nPrint =
                Foam::min(label(badCellIds.size()), label(100));

            for(label i=0; i<nPrint; ++i)
            {
                if( i ) Info << ',';
                Info << badCellIds[i];
            }

            if( badCellIds.size() > nPrint )
                Info << ",...";

            Info << ')' << endl;
        }
    };

    auditMeshZipUp("preAdd");

    meshModifier.addCells(cellsToAdd);

    auditMeshZipUp("postAdd");

    cellsToAdd.clear();

    meshModifier.reorderBoundaryFaces();

    auditMeshZipUp("postReorder");

    meshModifier.replaceBoundary
    (
        patchNames_,
        newBoundaryFaces,
        newBoundaryOwners,
        newBoundaryPatches
    );

    auditMeshZipUp("postReplace");

    PtrList<boundaryPatch>& boundaries = meshModifier.boundariesAccess();
    forAll(boundaries, patchI)
        boundaries[patchI].patchType() = patchTypes_[patchI];

    //- delete meshSurfaceEngine
    this->clearOut();

    Info << "Finished creating layer cells" << endl;
}

void boundaryLayers::createNewFacesFromPointsParallel
(
    const LongList<DynList<label, 4> >& faceCandidates,
    const LongList<labelPair>& candidatePatches
)
{
    const meshSurfaceEngine& mse = this->surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    const labelList& bp = mse.bp();
    const VRWGraph& bpAtProcs = mse.bpAtProcs();
    const labelList& globalPointLabel = mse.globalBoundaryPointLabel();
    const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();

    labelList otherFaceProc(faceCandidates.size(), -1);
    //- some faces may appear more than once
    //- such faces are ordinary internal faces
    VRWGraph pointFaceCandidates(nPoints_);
    forAll(faceCandidates, fI)
    {
        forAll(faceCandidates[fI], pI)
            pointFaceCandidates.append(faceCandidates[fI][pI], fI);
    }

    boolList duplicateFace(faceCandidates.size(), false);
    List<labelledPair> pointOfOrigin(faceCandidates.size());
    std::map<labelledPair, label> pointOfOriginToFaceLabel;
    forAll(faceCandidates, fI)
    {
        const DynList<label, 4>& f = faceCandidates[fI];

        const label pointI = f[0];

        const labelledPair lp
        (
            globalPointLabel[bp[pointI]],
            Pair<label>
            (
                patchKey_[candidatePatches[fI][0]],
                patchKey_[candidatePatches[fI][1]]
            )
        );

        if(
            pointOfOriginToFaceLabel.find(lp) != pointOfOriginToFaceLabel.end()
        )
        {
            duplicateFace[fI] = true;
            pointOfOrigin[fI] = lp;
            duplicateFace[pointOfOriginToFaceLabel[lp]] = true;
            continue;
        }

        pointOfOrigin[fI] = lp;

        pointOfOriginToFaceLabel.insert(std::make_pair(lp, fI));
    }

    //- find the processor patch for each processor boundary face
    //- the key of the algorithm is the point from which the face was created
    //- by sending the point label and the associated patches, it will be
    //- possible to find the other processor containing that face
    std::map<label, LongList<labelledPair> > exchangeData;
    const DynList<label>& neiProcs = mse.bpNeiProcs();
    forAll(neiProcs, procI)
    {
        const label neiProcI = neiProcs[procI];

        if( neiProcI == Pstream::myProcNo() )
            continue;

        if( exchangeData.find(neiProcI) == exchangeData.end() )
            exchangeData.insert
            (
                std::make_pair(neiProcI, LongList<labelledPair>())
            );
    }

    forAll(faceCandidates, fI)
    {
        if( duplicateFace[fI] )
            continue;

        const label bpI = bp[faceCandidates[fI][0]];

        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProcNo = bpAtProcs(bpI, procI);
            if( neiProcNo == Pstream::myProcNo() )
                continue;

            LongList<labelledPair>& dataToSend = exchangeData[neiProcNo];
            dataToSend.append(pointOfOrigin[fI]);
        }
    }

    //- exchange the data with other processors
    std::map<label, List<labelledPair> > receivedMap;
    help::exchangeMap(exchangeData, receivedMap);
    exchangeData.clear();

    for
    (
        std::map<label, List<labelledPair> >::const_iterator
        iter=receivedMap.begin();
        iter!=receivedMap.end();
        ++iter
    )
    {
        const List<labelledPair>& receivedData = iter->second;

        forAll(receivedData, i)
        {
            const labelledPair& lpp = receivedData[i];
            const label gpI = lpp.pairLabel();
            if( !globalToLocal.found(gpI) ) continue;
            const label pointI = bPoints[globalToLocal[gpI]];
            const labelPair& lp = lpp.pair();

            forAllRow(pointFaceCandidates, pointI, i)
            {
                const label fI = pointFaceCandidates(pointI, i);
                const DynList<label, 4>& f = faceCandidates[fI];

                const labelPair pk
                (
                    patchKey_[candidatePatches[fI][0]],
                    patchKey_[candidatePatches[fI][1]]
                );

                const labelPair rpk
                (
                    patchKey_[candidatePatches[fI][1]],
                    patchKey_[candidatePatches[fI][0]]
                );

                if(
                    (f[0] == pointI) && ((pk == lp) || (rpk == lp))
                )
                {
                    //- found the processor containing other face
                    auto fIt = pointOfOriginToFaceLabel.find(lpp);
                    if( fIt != pointOfOriginToFaceLabel.end() )
                        otherFaceProc[fIt->second] = iter->first;
                }
            }
        }
    }
    receivedMap.clear();

    //- sort the points in ascending order
    //- this ensures the correct order of faces at the processor boundaries
    sort(pointOfOrigin);

    Map<label> otherProcToProcPatch;
    forAll(mesh_.procBoundaries(), patchI)
    {
        const processorBoundaryPatch& wp = mesh_.procBoundaries()[patchI];
        otherProcToProcPatch.insert(wp.neiProcNo(), patchI);
    }

    //- store processor faces
    VRWGraph newProcFaces;
    labelLongList newProc;

    forAll(pointOfOrigin, i)
    {
        const label fI = pointOfOriginToFaceLabel[pointOfOrigin[i]];

        if( duplicateFace[fI] || (otherFaceProc[fI] == -1) )
            continue;

        if( !otherProcToProcPatch.found(otherFaceProc[fI]) )
        {
            otherProcToProcPatch.insert
            (
                otherFaceProc[fI],
                polyMeshGenModifier(mesh_).addProcessorPatch
                (
                    otherFaceProc[fI]
                )
            );
        }

        newProcFaces.appendList(faceCandidates[fI]);
        newProc.append(otherProcToProcPatch[otherFaceProc[fI]]);
    }

    polyMeshGenModifier(mesh_).addProcessorFaces(newProcFaces, newProc);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
