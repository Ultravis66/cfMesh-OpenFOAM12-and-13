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
#include "meshSurfaceEngine.H"
#include "demandDrivenData.H"
#include "helperFunctions.H"
#include "helperFunctionsPar.H"
#include "meshSurfaceCheckInvertedVertices.H"
#include "meshSurfaceCheckEdgeTypes.H"
#include "meshSurfacePartitioner.H"
#include "polyMeshGen2DEngine.H"

#include "labelledPoint.H"
#include <map>
#include <string>
#include <set>
#include "OFstream.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGLayer

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace
{
    inline bool boundaryPointsCanShareBLRamp
    (
        const VRWGraph& pPatches,
        const boolList& isBLPatch,
        const label bpI,
        const label nbpI
    )
    {
        // Same patch -- always safe
        forAllRow(pPatches, bpI, pi)
        {
            const label patchI = pPatches(bpI, pi);
            forAllRow(pPatches, nbpI, qi)
                if( patchI == pPatches(nbpI, qi) )
                    return true;
        }

        // Different patches: allow only if both touch BL wall patches
        bool bpHasBL(false), nbpHasBL(false);
        forAllRow(pPatches, bpI, pi)
        {
            const label pI = pPatches(bpI, pi);
            if( pI >= 0 && pI < label(isBLPatch.size()) && isBLPatch[pI] )
            { bpHasBL = true; break; }
        }
        forAllRow(pPatches, nbpI, qi)
        {
            const label pI = pPatches(nbpI, qi);
            if( pI >= 0 && pI < label(isBLPatch.size()) && isBLPatch[pI] )
            { nbpHasBL = true; break; }
        }
        return bpHasBL && nbpHasBL;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const meshSurfaceEngine& boundaryLayers::surfaceEngine() const
{
    if( !msePtr_ )
        msePtr_ = new meshSurfaceEngine(mesh_);

    return *msePtr_;
}

const meshSurfacePartitioner& boundaryLayers::surfacePartitioner() const
{
    if( !meshPartitionerPtr_ )
        meshPartitionerPtr_ = new meshSurfacePartitioner(surfaceEngine());

    return *meshPartitionerPtr_;
}

void boundaryLayers::findPatchesToBeTreatedTogether()
{
    if( geometryAnalysed_ )
        return;

    // Pre-mark zero-layer patches as treated so they are excluded
    // from concave-edge grouping with BL patches
    forAll(treatedPatch_, patchI)
        if( patchI < nLayersForPatch_.size()
         && nLayersForPatch_[patchI] == 0 )
            treatedPatch_[patchI] = true;

    forAll(treatPatchesWithPatch_, patchI)
        treatPatchesWithPatch_[patchI].append(patchI);

    const meshSurfaceEngine& mse = surfaceEngine();

    const pointFieldPMG& points = mesh_.points();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const edgeList& edges = mse.edges();
    const VRWGraph& eFaces = mse.edgeFaces();
    const labelList& boundaryFacePatches = mse.boundaryFacePatches();

    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const VRWGraph& pPatches = mPart.pointPatches();

    //- patches must be treated together if there exist a corner where
    //- more than three patches meet
    const labelHashSet& corners = mPart.corners();
    forAllConstIter(labelHashSet, corners, it)
    {
        const label bpI = it.key();

        // Diagnostic: log exact triple-junction corners (nFeat==3)
        // These are blade/hub/periodic junctions missed by > 3 condition.
        // Evaluate before changing to >= 3 to avoid over-grouping BL patches.
        if( mPart.numberOfFeatureEdgesAtPoint(bpI) == 3 )
        {
            Info << "BL triple-corner bpI=" << bpI
                 << " nPatches=" << pPatches.sizeOfRow(bpI)
                 << " patches=(";
            forAllRow(pPatches, bpI, ppI)
                Info << pPatches(bpI, ppI)
                     << (ppI+1<pPatches.sizeOfRow(bpI) ? "," : "");
            Info << ")" << endl;
        }
        if( mPart.numberOfFeatureEdgesAtPoint(bpI) > 3 )
        {
            labelHashSet commonPatches;
            DynList<label> allPatches;

            forAllRow(pPatches, bpI, patchI)
            {
                const DynList<label>& tpwp =
                    treatPatchesWithPatch_[pPatches(bpI, patchI)];

                forAll(tpwp, pJ)
                {
                    if( commonPatches.found(tpwp[pJ]) )
                        continue;

                    commonPatches.insert(tpwp[pJ]);
                    allPatches.append(tpwp[pJ]);
                }
            }

            forAllRow(pPatches, bpI, patchI)
                treatPatchesWithPatch_[pPatches(bpI, patchI)] = allPatches;

            # ifdef DEBUGLayer
            Info << "Corner " << bpI << " is shared by patches "
                << pPatches[bpI] << endl;
            Info << "All patches " << allPatches << endl;
            # endif
        }
    }

    //- patches must be treated together for concave geometries
    //- edgeClassification map counts the number of convex and concave edges
    //- for a given patch. The first counts convex edges and the second counts
    //- concave ones. If the number of concave edges is of the considerable
    //- percentage, it is treated as O-topology
    meshSurfaceCheckInvertedVertices vertexCheck(mse);
    const labelHashSet& invertedVertices = vertexCheck.invertedVertices();

    std::map<std::pair<label, label>, Pair<label> > edgeClassification;
    forAll(eFaces, eI)
    {
        if( eFaces.sizeOfRow(eI) != 2 )
            continue;

        //- check if the any of the face vertices is tangled
        const edge& e = edges[eI];
        if
        (
            !is2DMesh_ &&
            (invertedVertices.found(e[0]) || invertedVertices.found(e[1]))
        )
            continue;

        const label patch0 = boundaryFacePatches[eFaces(eI, 0)];
        const label patch1 = boundaryFacePatches[eFaces(eI, 1)];
        if( patch0 != patch1 )
        {
            std::pair<label, label> pp
            (
                Foam::min(patch0, patch1),
                Foam::max(patch0, patch1)
            );
            if( edgeClassification.find(pp) == edgeClassification.end() )
                edgeClassification.insert
                (
                    std::make_pair(pp, Pair<label>(0, 0))
                );

            const face& f1 = bFaces[eFaces(eI, 0)];
            const face& f2 = bFaces[eFaces(eI, 1)];

            if
            (
                !help::isSharedEdgeConvex(points, f1, f2) ||
                (help::angleBetweenFaces(points, f1, f2) > 0.75 * M_PI)
            )
            {
                ++edgeClassification[pp].second();
            }
            else
            {
                ++edgeClassification[pp].first();
            }
        }
    }

    if( Pstream::parRun() )
    {
        const labelList& bPoints = mse.boundaryPoints();

        //- check faces over processor edges
        const labelList& globalEdgeLabel = mse.globalBoundaryEdgeLabel();
        const Map<label>& globalToLocal = mse.globalToLocalBndEdgeAddressing();

        const DynList<label>& neiProcs = mse.beNeiProcs();
        const Map<label>& otherProcPatches = mse.otherEdgeFacePatch();
        const Map<label>& otherFaceProc = mse.otherEdgeFaceAtProc();

        //- send faces sharing processor edges to other processors
        //- faces are flattened into a single contiguous array
        const labelList& bp = mse.bp();
        const labelList& globalPointLabel = mse.globalBoundaryPointLabel();
        const Map<label>& globalPointToLocal =
            mse.globalToLocalBndPointAddressing();

        std::map<label, LongList<labelledPoint> > exchangePoints;
        forAll(neiProcs, procI)
        {
            exchangePoints.insert
            (
                std::make_pair(neiProcs[procI], LongList<labelledPoint>())
            );
        }

        //- store faces for sending
        forAllConstIter(Map<label>, otherFaceProc, it)
        {
            const label beI = it.key();

            if( eFaces.sizeOfRow(beI) == 0 )
                continue;

            const edge& e = edges[beI];

            if
            (
                !is2DMesh_ &&
                (invertedVertices.found(e[0]) || invertedVertices.found(e[1]))
            )
                continue;

            //- do not send data if the face on other processor
            //- is in the same patch
            if( otherProcPatches[beI] == boundaryFacePatches[eFaces(beI, 0)] )
                continue;

            const face& f = bFaces[eFaces(beI, 0)];

            const label neiProc = it();

            //- each face is sent as follows
            //- 1. global edge label
            //- 2. number of face nodes
            //- 3. faces nodes and vertex coordinates
            LongList<labelledPoint>& dps = exchangePoints[neiProc];
            dps.append(labelledPoint(globalEdgeLabel[beI], point()));
            dps.append(labelledPoint(f.size(), point()));
            forAll(f, pI)
            {
                dps.append
                (
                    labelledPoint
                    (
                        globalPointLabel[bp[f[pI]]],
                        points[f[pI]]
                    )
                );
            }
        }

        LongList<labelledPoint> receivedData;
        help::exchangeMap(exchangePoints, receivedData);

        //- receive faces from other processors
        Map<label> transferredPointToLocal;

        label counter(0);
        while( counter < receivedData.size() )
        {
            const label geI = receivedData[counter++].pointLabel();
            if( !globalToLocal.found(geI) )
            {
                // Skip -- unknown edge label on this processor.
                // Consume face-size + face-point entries to keep counter in sync.
                const label fSize = receivedData[counter++].pointLabel();
                counter += fSize;
                continue;
            }
            const label beI = globalToLocal[geI];

            DynList<label> f(receivedData[counter++].pointLabel());
            forAll(f, pI)
            {
                const labelledPoint& lp = receivedData[counter++];

                if( globalPointToLocal.found(lp.pointLabel()) )
                {
                    //- this point already exist on this processor
                    f[pI] = bPoints[globalPointToLocal[lp.pointLabel()]];
                }
                else
                {
                    //- this point does not exist on this processor
                    //- add it to the local list of points
                    //- it will be deleted when this procedure is finished
                    if( !transferredPointToLocal.found(lp.pointLabel()) )
                    {
                        //- this point has not yet been received
                        transferredPointToLocal.insert
                        (
                            lp.pointLabel(),
                            points.size()
                        );
                        mesh_.points().append(lp.coordinates());
                    }

                    f[pI] = transferredPointToLocal[lp.pointLabel()];
                }
            }

            const face& bf = bFaces[eFaces(beI, 0)];

            const label patch0 = boundaryFacePatches[eFaces(beI, 0)];
            const label patch1 = otherProcPatches[beI];

            std::pair<label, label> pp
            (
                Foam::min(patch0, patch1),
                Foam::max(patch0, patch1)
            );
            if( edgeClassification.find(pp) == edgeClassification.end() )
                edgeClassification.insert
                (
                    std::make_pair(pp, Pair<label>(0, 0))
                );

            if(
                (otherFaceProc[beI] > Pstream::myProcNo()) &&
                (
                    !help::isSharedEdgeConvex(points, bf, f) ||
                    (help::angleBetweenFaces(points, bf, f) > 0.75 * M_PI)
                )
            )
            {
                ++edgeClassification[pp].second();
            }
            else if( otherFaceProc[beI] > Pstream::myProcNo() )
            {
                ++edgeClassification[pp].first();
            }
        }

        //- set the size of points back to their original number
        mesh_.points().setSize(nPoints_);
    }

    std::map<std::pair<label, label>, Pair<label> >::const_iterator it;
    for(it=edgeClassification.begin();it!=edgeClassification.end();++it)
    {
        const std::pair<label, label>& edgePair = it->first;
        const Pair<label>& nConvexAndConcave = it->second;

        if( nConvexAndConcave.second() != 0 )
        {
            //- number of concave edges is greater than the number
            //- of the convex ones. Treat patches together.
            const label patch0 = edgePair.first;
            const label patch1 = edgePair.second;

            //- avoid adding unused patches in case of 2D meshing
            if( treatedPatch_[patch0] || treatedPatch_[patch1] )
                continue;

            treatPatchesWithPatch_[patch0].append(patch1);
            treatPatchesWithPatch_[patch1].append(patch0);
        }
    }

    if( Pstream::parRun() )
    {
        //- make sure that all processors have the same graph
        labelLongList flattenedPatches;
        forAll(treatPatchesWithPatch_, patchI)
        {
            if( treatPatchesWithPatch_[patchI].size() <= 1 )
                continue;

            flattenedPatches.append(patchI);
            flattenedPatches.append(treatPatchesWithPatch_[patchI].size());
            forAll(treatPatchesWithPatch_[patchI], itemI)
                flattenedPatches.append(treatPatchesWithPatch_[patchI][itemI]);
        }

        labelListList procPatches(Pstream::nProcs());
        procPatches[Pstream::myProcNo()].setSize(flattenedPatches.size());
        forAll(flattenedPatches, i)
            procPatches[Pstream::myProcNo()][i] = flattenedPatches[i];
        Pstream::gatherList(procPatches);
        Pstream::scatterList(procPatches);

        forAll(procPatches, procI)
        {
            if( procI == Pstream::myProcNo() )
            continue;

            const labelList& cPatches = procPatches[procI];
            label counter(0);

            while( counter < cPatches.size() )
            {
                const label patchI = cPatches[counter++];
                const label size = cPatches[counter++];
                for(label i=0;i<size;++i)
                    treatPatchesWithPatch_[patchI].appendIfNotIn
                    (
                        cPatches[counter++]
                    );
            }
        }
    }

    //- final adjusting of patches which shall be treated together
    boolList confirmed(treatPatchesWithPatch_.size(), false);
    forAll(treatPatchesWithPatch_, patchI)
    {
        if( treatPatchesWithPatch_[patchI].size() <= 1 )
        {
            confirmed[patchI] = true;
            continue;
        }

        if( confirmed[patchI] )
            continue;

        std::set<label> commonPatches;
        commonPatches.insert(patchI);

        DynList<label> front;
        front.append(patchI);
        confirmed[patchI] = true;

        while( front.size() )
        {
            const label fPatch = front.removeLastElement();

            forAll(treatPatchesWithPatch_[fPatch], i)
            {
                const label patchJ = treatPatchesWithPatch_[fPatch][i];

                if( confirmed[patchJ] )
                    continue;

                front.append(patchJ);
                confirmed[patchJ] = true;
                commonPatches.insert(patchJ);
                forAll(treatPatchesWithPatch_[patchJ], j)
                    commonPatches.insert(treatPatchesWithPatch_[patchJ][j]);
            }
        }

        forAllConstIter(std::set<label>, commonPatches, it)
        {
            const label patchJ = *it;

            treatPatchesWithPatch_[patchJ].clear();
            forAllConstIter(std::set<label>, commonPatches, iter)
                treatPatchesWithPatch_[patchJ].append(*iter);
        }
    }

    # ifdef DEBUGLayer
    for(it=edgeClassification.begin();it!=edgeClassification.end();++it)
    {
        const std::pair<label, label>& edgePair = it->first;
        const Pair<label>& nConvexAndConcave = it->second;
        Info << "Pair of patches " << edgePair.first << " "
            << edgePair.second << " is " << nConvexAndConcave << endl;
    }

    Info << "Patch names " << patchNames_ << endl;
    Info << "Treat patches with patch " << treatPatchesWithPatch_ << endl;

    label layerI(0), subsetId;
    boolList usedPatch(treatPatchesWithPatch_.size(), false);
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

    forAll(treatPatchesWithPatch_, patchI)
    {
        if( usedPatch[patchI] || (boundaries[patchI].patchSize() == 0) )
            continue;

        Info << "Adding layer subset " << layerI
             << " for patch " << patchI << endl;
        usedPatch[patchI] = true;
        subsetId = mesh_.addFaceSubset("layer_"+help::scalarToText(layerI));
        ++layerI;

        forAll(treatPatchesWithPatch_[patchI], i)
        {
            const label cPatch = treatPatchesWithPatch_[patchI][i];
            usedPatch[cPatch] = true;

            label start = boundaries[cPatch].patchStart();
            const label size = boundaries[cPatch].patchSize();
            for(label i=0;i<size;++i)
                mesh_.addFaceToSubset(subsetId, start++);
        }
    }

    mesh_.write();
    # endif

    geometryAnalysed_ = true;
}

void boundaryLayers::addLayerForPatch(const label patchLabel)
{
    if( treatedPatch_[patchLabel] )
        return;

    // Skip patches explicitly configured with nLayers==0
    if( patchLabel < nLayersForPatch_.size()
     && nLayersForPatch_[patchLabel] == 0 )
    {
        treatedPatch_[patchLabel] = true;
        return;
    }

    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

    if( returnReduce(boundaries[patchLabel].patchSize(), sumOp<label>()) == 0 )
        return;

    boolList treatPatches(boundaries.size(), false);
    if( patchWiseLayers_ )
    {
        forAll(treatPatchesWithPatch_[patchLabel], pI)
            treatPatches[treatPatchesWithPatch_[patchLabel][pI]] = true;
    }
    else
    {
        forAll(treatedPatch_, patchI)
            if( !treatedPatch_[patchI] )
                treatPatches[patchI] = true;
    }

    newLabelForVertex_.setSize(nPoints_);
    newLabelForVertex_ = -1;
    otherVrts_.clear();
    patchKey_.clear();

    createNewVertices(treatPatches);

    createNewFacesAndCells(treatPatches);

    forAll(treatPatches, patchI)
        if( treatPatches[patchI] )
            treatedPatch_[patchI] = true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh reference
boundaryLayers::boundaryLayers
(
    polyMeshGen& mesh,
    const dictionary& meshDict
)
:
    mesh_(mesh),
    msePtr_(NULL),
    meshPartitionerPtr_(NULL),
    patchWiseLayers_(true),
    terminateLayersAtConcaveEdges_(false),
    is2DMesh_(false),
    patchNames_(),
    patchTypes_(),
    treatedPatch_(),
    treatPatchesWithPatch_(),
    newLabelForVertex_(),
    otherVrts_(),
    patchKey_(),
    nPoints_(mesh.points().size()),
    geometryAnalysed_(false),
    blblFeatureAngleDeg_(40.0),
    blSharpEdgeAngleDeg_(75.0),
    blblCornerAcuteThreshold_(0.3),
    layerScaleRing1_(0.25),
    layerScaleRing2_(0.50),
    layerScaleRing3_(0.75),
    layerScaleRing4_(0.60),
    layerScaleRing5_(0.80),
    layerScaleRing6_(1.00),
    acuteCornerRing0_(0.0),
    acuteCornerRing1_(0.0),
    acuteCornerRing2_(0.05),
    acuteCornerRing3_(0.15),
    acuteCornerRing4_(0.35),
    acuteCornerRing5_(0.60),
    acuteCornerRing6_(1.00),
    virtualTopologyExclusion_(false),
    virtualTopoRing0_(0.0),
    virtualTopoRing1_(0.0),
    virtualTopoRing2_(0.05),
    gapFaceRingExclusion_(true),
    tripleJunctionFaceRingExclusion_(false),
    gapFaceRing0Scale_(0.02),
    gapFaceRing1Scale_(0.05),
    gapFaceRing2Scale_(0.20),
    gapFaceRing3Scale_(0.50),
    tripleJunctionProtectedRing0Scale_(1.0),
    gapLoserPatches_(),
    gapLoserPatchNames_(),
    gapLoserRing1Suppress_(true),
    gapLoserRing1MaxLayers_(1),
    gapLoserRing2MaxLayers_(2)
{
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    patchNames_.setSize(boundaries.size());
    patchTypes_.setSize(boundaries.size());
    forAll(boundaries, patchI)
    {
        patchNames_[patchI] = boundaries[patchI].patchName();
        patchTypes_[patchI] = boundaries[patchI].patchType();
    }

    treatedPatch_.setSize(boundaries.size());
    treatedPatch_ = false;

    treatPatchesWithPatch_.setSize(boundaries.size());

    // Per-patch nLayers: 0 means no BL (termination patch)
    nLayersForPatch_.setSize(boundaries.size(), 0);
    // Patch role: 0=BL, 1=TERMINATION, 2=NEUTRAL
    // Default: all nLayers==0 patches start as TERMINATION
    // User can override with terminationPatches / neutralPatches in meshDict
    patchRole_.setSize(boundaries.size(), 2); // default NEUTRAL
    if( meshDict.isDict("boundaryLayers") )
    {
        const dictionary& bndLayers = meshDict.subDict("boundaryLayers");
        label globalNLayers(0);
        if( bndLayers.found("nLayers") )
            globalNLayers = readLabel(bndLayers.lookup("nLayers"));
        nLayersForPatch_ = globalNLayers;

        // Ramp parameters - optional, defaults match hardcoded values
        if( bndLayers.found("blblFeatureAngleDeg") )
            blblFeatureAngleDeg_ =
                readScalar(bndLayers.lookup("blblFeatureAngleDeg"));
        if( bndLayers.found("blSharpEdgeAngleDeg") )
            blSharpEdgeAngleDeg_ =
                readScalar(bndLayers.lookup("blSharpEdgeAngleDeg"));
        if( bndLayers.found("blblCornerAcuteThreshold") )
            blblCornerAcuteThreshold_ =
                readScalar(bndLayers.lookup("blblCornerAcuteThreshold"));
        if( bndLayers.found("layerScaleRing1") )
            layerScaleRing1_ =
                readScalar(bndLayers.lookup("layerScaleRing1"));
        if( bndLayers.found("layerScaleRing2") )
            layerScaleRing2_ =
                readScalar(bndLayers.lookup("layerScaleRing2"));
        if( bndLayers.found("layerScaleRing3") )
            layerScaleRing3_ =
                readScalar(bndLayers.lookup("layerScaleRing3"));
        if( bndLayers.found("layerScaleRing4") )
            layerScaleRing4_ =
                readScalar(bndLayers.lookup("layerScaleRing4"));
        if( bndLayers.found("layerScaleRing5") )
            layerScaleRing5_ =
                readScalar(bndLayers.lookup("layerScaleRing5"));
        if( bndLayers.found("layerScaleRing6") )
            layerScaleRing6_ =
                readScalar(bndLayers.lookup("layerScaleRing6"));
        if( bndLayers.found("acuteCornerRing0") )
            acuteCornerRing0_ = readScalar(bndLayers.lookup("acuteCornerRing0"));
        if( bndLayers.found("acuteCornerRing1") )
            acuteCornerRing1_ = readScalar(bndLayers.lookup("acuteCornerRing1"));
        if( bndLayers.found("acuteCornerRing2") )
            acuteCornerRing2_ = readScalar(bndLayers.lookup("acuteCornerRing2"));
        if( bndLayers.found("acuteCornerRing3") )
            acuteCornerRing3_ = readScalar(bndLayers.lookup("acuteCornerRing3"));
        if( bndLayers.found("acuteCornerRing4") )
            acuteCornerRing4_ = readScalar(bndLayers.lookup("acuteCornerRing4"));
        if( bndLayers.found("acuteCornerRing5") )
            acuteCornerRing5_ = readScalar(bndLayers.lookup("acuteCornerRing5"));
        if( bndLayers.found("acuteCornerRing6") )
            acuteCornerRing6_ = readScalar(bndLayers.lookup("acuteCornerRing6"));

        if( bndLayers.found("virtualTopologyExclusion") )
            virtualTopologyExclusion_ =
                Switch(bndLayers.lookup("virtualTopologyExclusion"));
        if( bndLayers.found("virtualTopoRing0") )
            virtualTopoRing0_ = readScalar(bndLayers.lookup("virtualTopoRing0"));
        if( bndLayers.found("virtualTopoRing1") )
            virtualTopoRing1_ = readScalar(bndLayers.lookup("virtualTopoRing1"));
        if( bndLayers.found("virtualTopoRing2") )
            virtualTopoRing2_ = readScalar(bndLayers.lookup("virtualTopoRing2"));
        if( bndLayers.found("gapFaceRingExclusion") )
            gapFaceRingExclusion_ = Switch(bndLayers.lookup("gapFaceRingExclusion"));
        if( bndLayers.found("tripleJunctionFaceRingExclusion") )
            tripleJunctionFaceRingExclusion_ =
                Switch(bndLayers.lookup("tripleJunctionFaceRingExclusion"));
        if( bndLayers.found("tripleJunctionSuppressPatches") )
        {
            const wordList suppressNames
            (
                bndLayers.lookup("tripleJunctionSuppressPatches")
            );
            const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
            forAll(boundaries, pI)
                forAll(suppressNames, sI)
                    if( boundaries[pI].patchName() == suppressNames[sI] )
                        tripleJunctionSuppressPatches_.insert(pI);

            Info << "Triple-junction suppress patches: found "
                 << tripleJunctionSuppressPatches_.size()
                 << " patch indices from " << suppressNames << endl;
        }
        if( bndLayers.found("gapFaceRing0Scale") )
            gapFaceRing0Scale_ = readScalar(bndLayers.lookup("gapFaceRing0Scale"));
        if( bndLayers.found("gapFaceRing1Scale") )
            gapFaceRing1Scale_ = readScalar(bndLayers.lookup("gapFaceRing1Scale"));
        if( bndLayers.found("gapFaceRing2Scale") )
            gapFaceRing2Scale_ = readScalar(bndLayers.lookup("gapFaceRing2Scale"));
        if( bndLayers.found("gapFaceRing3Scale") )
            gapFaceRing3Scale_ = readScalar(bndLayers.lookup("gapFaceRing3Scale"));
        if( bndLayers.found("gapLoserRing1Suppress") )
            gapLoserRing1Suppress_ = Switch(bndLayers.lookup("gapLoserRing1Suppress"));
        if( bndLayers.found("gapLoserRing1MaxLayers") )
            gapLoserRing1MaxLayers_ = readLabel(bndLayers.lookup("gapLoserRing1MaxLayers"));
        if( bndLayers.found("gapLoserRing2MaxLayers") )
            gapLoserRing2MaxLayers_ = readLabel(bndLayers.lookup("gapLoserRing2MaxLayers"));
        if( bndLayers.found("tripleJunctionProtectedRing0Scale") )
            tripleJunctionProtectedRing0Scale_ =
                readScalar(bndLayers.lookup("tripleJunctionProtectedRing0Scale"));

        if( bndLayers.isDict("patchBoundaryLayers") )
        {
            const dictionary& patchBndLayers =
                bndLayers.subDict("patchBoundaryLayers");
            forAll(boundaries, patchI)
            {
                const word& nm = boundaries[patchI].patchName();
                if( patchBndLayers.isDict(nm) )
                {
                    const dictionary& pd = patchBndLayers.subDict(nm);
                    if( pd.found("nLayers") )
                        nLayersForPatch_[patchI] =
                            readLabel(pd.lookup("nLayers"));
                }
            }
        }

        // Assign patch roles after all nLayers are set
        // BL patches: nLayers > 0
        // Termination patches: explicitly listed in terminationPatches
        // Neutral patches: everything else (periodic, symmetry, etc)
        wordList terminationPatchNames;
        if( bndLayers.found("terminationPatches") )
            terminationPatchNames = wordList(bndLayers.lookup("terminationPatches"));

        forAll(boundaries, patchI)
        {
            if( nLayersForPatch_[patchI] > 0 )
            {
                patchRole_[patchI] = 0; // BL_PATCH
            }
            else
            {
                // Check if explicitly listed as termination
                bool isTermination = false;
                forAll(terminationPatchNames, tI)
                    if( terminationPatchNames[tI] == boundaries[patchI].patchName() )
                    { isTermination = true; break; }
                patchRole_[patchI] = isTermination ? 1 : 2; // TERMINATION or NEUTRAL
            }
        }

        // Report patch role assignment
        label nBLR=0, nTermR=0, nNeutR=0;
        forAll(boundaries, patchI)
        {
            if( patchRole_[patchI] == 0 ) ++nBLR;
            else if( patchRole_[patchI] == 1 ) ++nTermR;
            else ++nNeutR;
        }
        Info << "Patch role assignment: "
             << nBLR << " BL, "
             << nTermR << " termination, "
             << nNeutR << " neutral" << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

boundaryLayers::~boundaryLayers()
{
    clearOut();

    if( Pstream::parRun() )
        polyMeshGenModifier(mesh_).removeUnusedVertices();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayers::addLayerForPatch(const word& patchName)
{
    if( !geometryAnalysed_ )
        findPatchesToBeTreatedTogether();

    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

    forAll(boundaries, patchI)
        if( boundaries[patchI].patchName() == patchName )
            addLayerForPatch(patchI);
}

void boundaryLayers::createOTopologyLayers()
{
    patchWiseLayers_ = false;
}

void boundaryLayers::terminateLayersAtConcaveEdges()
{
    terminateLayersAtConcaveEdges_ = true;
}

void boundaryLayers::detectBLNoBlTransitionEdges() const
{
    const meshSurfaceEngine& mse = surfaceEngine();
    const VRWGraph& edgeFaces = mse.edgeFaces();
    const labelList& boundaryFacePatches = mse.boundaryFacePatches();
    const edgeList& edges = mse.edges();
    const labelList& bPoints = mse.boundaryPoints();

    // Classify patches from patchRole_ -- single source of truth
    boolList isBLPatch(patchNames_.size(), false);
    forAll(patchNames_, patchI)
        if( patchI < label(patchRole_.size()) && patchRole_[patchI] == 0 )
            isBLPatch[patchI] = true;

    // Build global-to-boundary-point reverse map once, O(N)
    Map<label> globalToBP;
    forAll(bPoints, bpI)
        globalToBP.insert(bPoints[bpI], bpI);

    blNoBlEdges_.clear();
    blNoBlEdgePoints_.clear();
    blNoBlPointPatch_.clear();
    blNeutralEdgePoints_.clear();
    blNeutralPointPatch_.clear();

    forAll(edges, eI)
    {
        bool hasBL   = false;
        bool hasNoBL = false;
        label blPatchI = -1;

        forAllRow(edgeFaces, eI, efI)
        {
            const label faceI  = edgeFaces(eI, efI);
            const label patchI = boundaryFacePatches[faceI];
            if( patchI < 0 || patchI >= label(isBLPatch.size()) )
                continue;
            if( patchRole_.size() > 0 && patchI < label(patchRole_.size()) )
            {
                if( patchRole_[patchI] == 0 )
                {
                    hasBL = true;
                    blPatchI = patchI;
                }
                else if( patchRole_[patchI] == 1 )
                {
                    hasNoBL = true;
                }
                // patchRole_==2 (neutral) ignored here
            }
        }

        if( !hasBL || !hasNoBL ) continue;

        blNoBlEdges_.insert(eI);

        const edge& e = edges[eI];
        for( label ei = 0; ei < 2; ++ei )
        {
            const label gp = (ei == 0) ? e[0] : e[1];
            Map<label>::const_iterator it = globalToBP.find(gp);
            if( it == globalToBP.end() ) continue;
            const label bpI = it();
            blNoBlEdgePoints_.insert(bpI);
            // Store BL-side patch for patch-constrained projection
            // Only insert if not already stored (first BL patch wins)
            if( !blNoBlPointPatch_.found(bpI) )
                blNoBlPointPatch_.insert(bpI, blPatchI);
        }
    }

    // Detect BL/neutral edge points: two-patch points where
    // one patch is BL and the other is neutral (periodic/symmetry).
    // These are NOT termination points but need special handling
    // to prevent BL extrusion across the neutral boundary.
    forAll(edges, eI)
    {
        if( edgeFaces.sizeOfRow(eI) != 2 ) continue;
        const label fA = edgeFaces(eI, 0);
        const label fB = edgeFaces(eI, 1);
        if( fA < 0 || fA >= boundaryFacePatches.size() ) continue;
        if( fB < 0 || fB >= boundaryFacePatches.size() ) continue;
        const label pA = boundaryFacePatches[fA];
        const label pB = boundaryFacePatches[fB];
        if( pA < 0 || pA >= label(isBLPatch.size()) ) continue;
        if( pB < 0 || pB >= label(isBLPatch.size()) ) continue;
        if( pA == pB ) continue;
        const bool blA = isBLPatch[pA];
        const bool blB = isBLPatch[pB];
        // One BL, one no-BL
        if( !((blA && !blB) || (!blA && blB)) ) continue;
        // Check role: neutral = patchRole 2
        const label noBlPatch = blA ? pB : pA;
        if( patchRole_.size() > 0
         && noBlPatch < label(patchRole_.size())
         && patchRole_[noBlPatch] != 2 ) continue; // not neutral
        const edge& e = edges[eI];
        for( label ei = 0; ei < 2; ++ei )
        {
            const label gp = (ei==0) ? e[0] : e[1];
            Map<label>::const_iterator it = globalToBP.find(gp);
            if( it == globalToBP.end() ) continue;
            const label bpI = it();
            blNeutralEdgePoints_.insert(bpI);
            if( !blNeutralPointPatch_.found(bpI) )
                blNeutralPointPatch_.insert(bpI, blA ? pA : pB);
        }
    }

    Info << "BL/no-BL transition edge pre-detection: "
         << blNoBlEdges_.size() << " transition edges, "
         << blNoBlEdgePoints_.size() << " interface points" << endl;
    Info << "BL/neutral edge points detected: "
         << blNeutralEdgePoints_.size()
         << " (blade/periodic-style junctions)" << endl;
}

void boundaryLayers::markConcaveEdgePoints(boolList& skipPoint) const
{
    const meshSurfaceEngine& mse = surfaceEngine();
    const VRWGraph& edgeFaces = mse.edgeFaces();
    const labelList& boundaryFacePatches = mse.boundaryFacePatches();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const pointFieldPMG& points = mesh_.points();
    const edgeList& edges = mse.edges();
    const labelList& bPoints = mse.boundaryPoints();
    const VRWGraph& pointPoints = mse.pointPoints();
    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const VRWGraph& pPatches = mPart.pointPatches();

    // Classify patches using patchRole_ -- the single source of truth.
    // patchRole_: 0=BL, 1=TERMINATION, 2=NEUTRAL
    // Previously derived from nLayersForPatch_==0 which treated neutral
    // patches (periodic/symmetry) the same as termination (inlet/outlet).
    // That caused over-suppression at blade/periodic junctions.
    boolList isBLPatch(patchNames_.size(), false);
    boolList isTerminationPatch(patchNames_.size(), false);
    forAll(patchNames_, patchI)
    {
        if( patchI < patchRole_.size() )
        {
            isBLPatch[patchI]          = (patchRole_[patchI] == 0);
            isTerminationPatch[patchI] = (patchRole_[patchI] == 1);
            // patchRole_==2 (neutral) is neither BL nor termination
            // -- no suppression triggered at periodic/symmetry junctions
        }
    }

    // Mark which boundary points belong to at least one BL patch
    boolList boundaryPointIsBL(bPoints.size(), false);
    forAll(bPoints, bpI)
        forAllRow(pPatches, bpI, pI)
        {
            const label patchI = pPatches(bpI, pI);
            if( patchI >= 0 && patchI < isBLPatch.size()
             && isBLPatch[patchI] )
            {
                boundaryPointIsBL[bpI] = true;
                break;
            }
        }

    labelList meshToBnd(mesh_.points().size(), -1);
    forAll(bPoints, bpI)
        meshToBnd[bPoints[bpI]] = bpI;

    // Initialize scale fields
    layerScale_.setSize(bPoints.size(), 1.0);
    blSuppressReason_.setSize(bPoints.size(), 0);
    blRampRing_.setSize(bPoints.size(), 0);
    blRampSeedReason_.setSize(bPoints.size(), 0);
    zeroDistPoints_.setSize(bPoints.size(), false);
    boolList zeroPts(bPoints.size(), false);

    // Inject externally detected gap points (mesh point labels) into zeroPts.
    if( gapPoints_.size() > 0 )
    {
        labelList meshToBnd(mesh_.points().size(), -1);
        forAll(bPoints, bpI)
            meshToBnd[bPoints[bpI]] = bpI;
        label nGapSuppressed = 0;
        forAllConstIter(labelHashSet, gapPoints_, it)
        {
            const label meshPtI = it.key();
            if( meshPtI < 0 || meshPtI >= label(meshToBnd.size()) ) continue;
            const label bpI = meshToBnd[meshPtI];
            if( bpI >= 0 && bpI < label(bPoints.size()) )
            {
                zeroDistPoints_[bpI] = true;
                zeroPts[bpI] = true;
                layerScale_[bpI] = 0.0;
                blSuppressReason_[bpI] = 1;
                ++nGapSuppressed;
            }
        }
        Info << "Gap detection: suppressed "
             << nGapSuppressed
             << " BL points in thin clearance regions" << endl;

        // Gap taper rings: BFS outward from gap points to create smooth
        // BL restart instead of hard zero cliff
        const VRWGraph& ptPts = pointPoints;
        const label nBP = bPoints.size();

        // Ring 1: immediate neighbours -- strong suppression
        boolList gapRing0(nBP, false);
        forAllConstIter(labelHashSet, gapPoints_, it)
        {
            const label meshPtI = it.key();
            if( meshPtI < 0 || meshPtI >= label(meshToBnd.size()) ) continue;
            const label bpI = meshToBnd[meshPtI];
            if( bpI >= 0 && bpI < nBP ) gapRing0[bpI] = true;
        }

        boolList gapRing1(nBP, false);
        forAll(gapRing0, bpI)
        {
            if( !gapRing0[bpI] ) continue;
            forAllRow(ptPts, bpI, ppI)
            {
                const label nbpI = ptPts(bpI, ppI);
                if( nbpI < 0 || nbpI >= nBP ) continue;
                if( gapRing0[nbpI] ) continue;
                gapRing1[nbpI] = true;
                layerScale_[nbpI] = Foam::min(layerScale_[nbpI], 0.05);
            }
        }

        // Ring 2: taper onset
        boolList gapRing2(nBP, false);
        forAll(gapRing1, bpI)
        {
            if( !gapRing1[bpI] ) continue;
            forAllRow(ptPts, bpI, ppI)
            {
                const label nbpI = ptPts(bpI, ppI);
                if( nbpI < 0 || nbpI >= nBP ) continue;
                if( gapRing0[nbpI] || gapRing1[nbpI] ) continue;
                gapRing2[nbpI] = true;
                layerScale_[nbpI] = Foam::min(layerScale_[nbpI], 0.20);
            }
        }

        // Ring 3: gentle restart (scale via meshDict gapFaceRing3Scale)
        boolList gapRing3(nBP, false);
        forAll(gapRing2, bpI)
        {
            if( !gapRing2[bpI] ) continue;
            forAllRow(ptPts, bpI, ppI)
            {
                const label nbpI = ptPts(bpI, ppI);
                if( nbpI < 0 || nbpI >= nBP ) continue;
                if( gapRing0[nbpI] || gapRing1[nbpI] || gapRing2[nbpI] ) continue;
                gapRing3[nbpI] = true;
                layerScale_[nbpI] = Foam::min(layerScale_[nbpI], gapFaceRing3Scale_);
            }
        }
    }

    // Mark exact transition edge points
    label nTransitionEdges = 0;
    forAll(edges, edgeI)
    {
        if( edgeFaces.sizeOfRow(edgeI) != 2 ) continue;
        const label f0 = edgeFaces(edgeI, 0);
        const label f1 = edgeFaces(edgeI, 1);
        if( f0 < 0 || f0 >= boundaryFacePatches.size() ) continue;
        if( f1 < 0 || f1 >= boundaryFacePatches.size() ) continue;
        const label patch0 = boundaryFacePatches[f0];
        const label patch1 = boundaryFacePatches[f1];
        if( patch0 < 0 || patch0 >= patchNames_.size() ) continue;
        if( patch1 < 0 || patch1 >= patchNames_.size() ) continue;
        if( patch0 == patch1 ) continue;
        const bool bl0 = isBLPatch[patch0];
        const bool bl1 = isBLPatch[patch1];
        const bool term0 = isTerminationPatch[patch0];
        const bool term1 = isTerminationPatch[patch1];
        if( (bl0 && term1) || (bl1 && term0) )
        {
            const edge& e = edges[edgeI];
            const label bp0 = meshToBnd[e[0]];
            const label bp1 = meshToBnd[e[1]];
            if( bp0 >= 0 )
            {
                zeroDistPoints_[bp0] = true;
                layerScale_[bp0] = 0.02;
                if( bp0 >= 0 && bp0 < blSuppressReason_.size() )
                    blSuppressReason_[bp0] = 2;
                zeroPts[bp0] = true;
            }
            if( bp1 >= 0 )
            {
                zeroDistPoints_[bp1] = true;
                layerScale_[bp1] = 0.02;
                if( bp1 >= 0 && bp1 < blSuppressReason_.size() )
                    blSuppressReason_[bp1] = 2;
                zeroPts[bp1] = true;
            }
            ++nTransitionEdges;
        }
    }

    // Triple-junction suppression:
    // Points on 2+ BL patches + 1+ termination patch are geometrically
    // overconstrained - suppress BL and let ramp handle transition
    label nTriple = 0;
    label nPts2 = 0, nPts3 = 0, nPts4plus = 0;
    forAll(bPoints, bpI)
    {
        if( !boundaryPointIsBL[bpI] ) continue;
        if( zeroPts[bpI] ) continue;
        label nPatches = 0, nBLPatches = 0, nTermPatches = 0, nNeutPatches = 0;
        DynList<label> seenPatches;
        forAllRow(pPatches, bpI, pI)
        {
            const label patchI = pPatches(bpI, pI);
            if( patchI < 0 || patchI >= label(patchNames_.size()) ) continue;
            if( seenPatches.contains(patchI) ) continue;
            seenPatches.append(patchI);
            ++nPatches;
            if( patchI < label(nLayersForPatch_.size())
             && nLayersForPatch_[patchI] > 0 )
                ++nBLPatches;
            else if( patchRole_.size() > patchI && patchRole_[patchI] == 2 )
                ++nNeutPatches;
            else
                ++nTermPatches;
        }
        if( nPatches == 2 ) ++nPts2;
        else if( nPatches == 3 ) ++nPts3;
        else if( nPatches > 3 ) ++nPts4plus;
        // BL + termination corners
        if( nPatches >= 3 && nBLPatches >= 2 && nTermPatches >= 1 )
        {
            zeroDistPoints_[bpI] = true;
            layerScale_[bpI] = 0.02;
            blSuppressReason_[bpI] = 3;
            zeroPts[bpI] = true;
            ++nTriple;
        }
        // BL+BL+neutral corners (blade/shroud/periodic, blade/hub/periodic)
        // Measure patch-normal spread to classify acute vs mild corners
        else if( nPatches >= 3 && nBLPatches >= 2 && nNeutPatches >= 1 )
        {
            // Compute average normal per BL patch at this corner
            const VRWGraph& ptFaces = mse.pointFaces();
            const faceList::subList& bFaces = mse.boundaryFaces();
            const labelList& bFacePatches = mse.boundaryFacePatches();

            Map<vector> blPatchNormals;
            forAllRow(ptFaces, bpI, pfI)
            {
                const label faceI = ptFaces(bpI, pfI);
                const label patchI = bFacePatches[faceI];
                if( patchI < 0 || patchI >= label(patchRole_.size()) ) continue;
                if( patchRole_[patchI] != 0 ) continue; // BL patches only
                const face& f = bFaces[faceI];
                if( f.size() < 3 ) continue;
                vector fn = vector::zero;
                const point& fp0 = points[f[0]];
                for(label fi=1; fi<f.size()-1; ++fi)
                    fn += (points[f[fi]]-fp0)^(points[f[fi+1]]-fp0);
                if( blPatchNormals.found(patchI) )
                    blPatchNormals[patchI] += fn;
                else
                    blPatchNormals.insert(patchI, fn);
            }

            // Find minimum dot product between any two BL patch normals
            scalar minDot = GREAT;
            DynList<vector> blNorms;
            forAllConstIter(Map<vector>, blPatchNormals, it)
            {
                const scalar magN = mag(it());
                if( magN > VSMALL )
                    blNorms.append(it() / magN);
            }
            for(label i=0; i<blNorms.size(); ++i)
                for(label j=i+1; j<blNorms.size(); ++j)
                    minDot = Foam::min(minDot, blNorms[i] & blNorms[j]);

            // Acute corner: BL patch normals diverge strongly
            // Threshold read from meshDict via blblCornerAcuteThreshold_
            // Default 0.3 (~73 degrees). Lower = more conservative.
            const bool isAcute = (minDot < blblCornerAcuteThreshold_ && minDot < GREAT);

            zeroDistPoints_[bpI] = true;
            layerScale_[bpI] = isAcute ? 0.0 : 0.15;
            zeroPts[bpI] = true;
            blblCornerPoints_.insert(bpI);
            if( isAcute )
                blblAcuteCornerPoints_.insert(bpI);
            ++nTriple;

            # ifdef DEBUGLayer
            Info << "BL+BL+neutral corner bpI=" << bpI
                 << " minDot=" << minDot
                 << " acute=" << isAcute
                 << " layerScale=" << layerScale_[bpI] << endl;
            # endif
        }
    }
    Info << "BL triple-junction stats: "
         << "nPts2=" << nPts2
         << " nPts3=" << nPts3
         << " nPts4plus=" << nPts4plus
         << " suppressed=" << nTriple
         << " blblCornerPts=" << blblCornerPoints_.size()
         << " acuteCornerPts=" << blblAcuteCornerPoints_.size() << endl;

    // Topology-aware scale assignment:
    // Upgrade transition point suppression based on point topology class.
    // Corner points (3+ patches) at BL/no-BL transitions are geometrically
    // overconstrained -- full suppress regardless of angle.
    // Two-patch edge points get full suppress only if patch angle is sharp.
    {
        const scalar cosSharp =
            Foam::cos(blSharpEdgeAngleDeg_ * M_PI / 180.0);
        Info << "BL sharp-edge suppression: angle=" << blSharpEdgeAngleDeg_
             << " deg (cos=" << cosSharp << ")" << endl;
        label nCornerSuppressed = 0;
        label nEdgeSuppressed = 0;
        // Diagnostic histogram of patch normal dot products
        label nDotNeg = 0, nDot0to05 = 0, nDot05to07 = 0, nDot07plus = 0;
        forAll(bPoints, bpI)
        {
            if( !zeroPts[bpI] ) continue;
            if( layerScale_[bpI] <= 0.0 ) continue;
            // Count unique patches at this point
            DynList<label> ptPatchList;
            label nBLPt = 0, nTermPt = 0, nNeutPt = 0;
            forAllRow(pPatches, bpI, pI)
            {
                const label patchI = pPatches(bpI, pI);
                if( patchI < 0 || patchI >= label(patchNames_.size()) ) continue;
                if( ptPatchList.contains(patchI) ) continue;
                ptPatchList.append(patchI);
                const label role = patchRole_[patchI];
                if( role == 0 ) ++nBLPt;
                else if( role == 1 ) ++nTermPt;
                else ++nNeutPt;
            }
            const label nPt = ptPatchList.size();
            // Corner: 3+ patches with BL + explicit termination -- full suppress
            // Neutral patches (periodic etc) do not trigger suppression
            if( nPt >= 3 && nBLPt >= 1 && nTermPt >= 1 )
            {
                layerScale_[bpI] = 0.0;
                blSuppressReason_[bpI] = 4;
                ++nCornerSuppressed;

                if( nCornerSuppressed <= 200 )
                {
                    const point& p = points[bPoints[bpI]];
                    Info << "TERMDIAG corner bpI=" << bpI
                         << " nPt=" << nPt
                         << " nBL=" << nBLPt
                         << " nTerm=" << nTermPt
                         << " x=" << p.x()
                         << " y=" << p.y()
                         << " z=" << p.z()
                         << endl;
                }

                continue;
            }

            // Two-patch edge: suppress only if BL meets explicit termination
            // patch at sharp angle. Neutral patches never trigger suppression.
            if( nPt == 2 && nBLPt >= 1 && nTermPt >= 1
             && blNoBlEdgePoints_.found(bpI) )
            {
                // Compute normals for the two patches
                const VRWGraph& ptFaces2 = mse.pointFaces();
                DynList<vector> patchNormals;
                forAll(ptPatchList, pi)
                {
                    vector n = vector::zero;
                    label nf = 0;
                    forAllRow(ptFaces2, bpI, pfI)
                    {
                        const label faceI = ptFaces2(bpI, pfI);
                        if( boundaryFacePatches[faceI] != ptPatchList[pi] ) continue;
                        const face& f = bFaces[faceI];
                        vector fn = vector::zero;
                        const point& fp0 = points[f[0]];
                        for(label fi=1; fi<f.size()-1; ++fi)
                            fn += (points[f[fi]]-fp0)^(points[f[fi+1]]-fp0);
                        if( mag(fn) > VSMALL ) { n += fn/mag(fn); ++nf; }
                    }
                    if( nf > 0 ) n /= scalar(nf);
                    if( mag(n) > VSMALL ) n /= mag(n);
                    patchNormals.append(n);
                }
                if( patchNormals.size() == 2 )
                {
                    const scalar dotProd = patchNormals[0] & patchNormals[1];
                    if( dotProd < 0 ) ++nDotNeg;
                    else if( dotProd < 0.5 ) ++nDot0to05;
                    else if( dotProd < 0.707 ) ++nDot05to07;
                    else ++nDot07plus;
                    if( dotProd < cosSharp )
                    {
                        layerScale_[bpI] = 0.0;
                        blSuppressReason_[bpI] = 5;
                        ++nEdgeSuppressed;
                    }
                }
            }
        }
        Info << "Topology-aware BL suppression: "
             << nCornerSuppressed << " corner points fully suppressed, "
             << nEdgeSuppressed << " sharp edge points fully suppressed" << endl;
        Info << "Normal dot-product histogram: "
             << "dot<0: " << nDotNeg
             << " 0-0.5: " << nDot0to05
             << " 0.5-0.707: " << nDot05to07
             << " >0.707: " << nDot07plus << endl;
    }

    // BL/BL sharp-junction suppression:
    // Points touching 2+ BL patches where normals diverge sharply
    // (blade+hub, blade+shroud) create degenerate layer cells.
    // Threshold: 40 degrees between patch normals.
    {
        blblJunctionPoints_.clear();
        blblJunctionClass_.clear();
        const scalar cosThresh = Foam::cos(blblFeatureAngleDeg_ * M_PI / 180.0);
        label nBLBL = 0;
        label nBLBLHard = 0;
        label nBLBLModerate = 0;
        label nBLBLMild = 0;
        const VRWGraph& ptFaces = mse.pointFaces();
        forAll(bPoints, bpI)
        {
            if( zeroPts[bpI] ) continue;
            if( !boundaryPointIsBL[bpI] ) continue;

            // Collect unique BL patches at this point
            DynList<label> blPatches;
            forAllRow(pPatches, bpI, pI)
            {
                const label patchI = pPatches(bpI, pI);
                if( patchI >= 0
                 && patchI < label(nLayersForPatch_.size())
                 && nLayersForPatch_[patchI] > 0 )
                    blPatches.appendIfNotIn(patchI);
            }
            if( blPatches.size() < 2 ) continue;

            // Compute average face normal per BL patch at this point
            DynList<vector> avgNormals;
            forAll(blPatches, i)
            {
                vector n = vector::zero;
                label nf = 0;
                forAllRow(ptFaces, bpI, pfI)
                {
                    const label faceI = ptFaces(bpI, pfI);
                    if( boundaryFacePatches[faceI] != blPatches[i] )
                        continue;
                    const face& f = bFaces[faceI];
                    vector fn = vector::zero;
                    const point& fp0 = points[f[0]];
                    for(label pi=1; pi<f.size()-1; ++pi)
                        fn += (points[f[pi]]-fp0)^(points[f[pi+1]]-fp0);
                    if( mag(fn) > VSMALL )
                    {
                        n += fn / mag(fn);
                        ++nf;
                    }
                }
                if( nf > 0 ) n /= scalar(nf);
                if( mag(n) > VSMALL ) n /= mag(n);
                avgNormals.append(n);
            }

            // Suppress if any pair of patch normals diverges sharply
            bool sharpJunction = false;
            for(label i=0; i<avgNormals.size()-1; ++i)
                for(label j=i+1; j<avgNormals.size(); ++j)
                    if( (avgNormals[i] & avgNormals[j]) < cosThresh )
                        sharpJunction = true;

            if( sharpJunction )
            {
                //- Commercial-grade blblSharp policy v1:
                //- Graduated taper instead of hard zero for all junctions.
                //- Atlas: blblSharp=2599 dominant dropout + 1640 ramp victims.
                scalar minDot = scalar(1.0);
                for(label i=0; i<avgNormals.size()-1; ++i)
                    for(label j=i+1; j<avgNormals.size(); ++j)
                    {
                        const scalar d = avgNormals[i] & avgNormals[j];
                        if( d < minDot ) minDot = d;
                    }
                const scalar cosHard = Foam::cos(scalar(75.0)*M_PI/180.0);
                const scalar cosMild = Foam::cos(scalar(60.0)*M_PI/180.0);
                scalar taperFloor = scalar(0.0);
                label junctionClass = 0; // 0=hard 1=moderate 2=mild
                if( minDot < cosHard )
                {
                    taperFloor = scalar(0.0);   // hard corner -- full suppress
                    junctionClass = 0;
                    ++nBLBLHard;
                }
                else if( minDot < cosMild )
                {
                    taperFloor = scalar(0.25);  // moderate -- quarter layer
                    junctionClass = 1;
                    ++nBLBLModerate;
                }
                else
                {
                    taperFloor = scalar(0.50);  // mild -- half layer
                    junctionClass = 2;
                    ++nBLBLMild;
                }
                layerScale_[bpI] = taperFloor;
                blSuppressReason_[bpI] = 6;
                if( taperFloor < 0.01 )
                {
                    zeroDistPoints_[bpI] = true;
                    zeroPts[bpI] = true;
                }
                else
                {
                    zeroDistPoints_[bpI] = false;
                    zeroPts[bpI] = false;
                }
                ++nBLBL;
                blblJunctionPoints_.insert(bpI);
                if( blblJunctionClass_.found(bpI) )
                    blblJunctionClass_[bpI] = junctionClass;
                else
                    blblJunctionClass_.insert(bpI, junctionClass);
            }
        }
        Info << "BL/BL sharp-junction suppression: "
             << nBLBL << " points total"
             << " (hard=" << nBLBLHard
             << " moderate=" << nBLBLModerate
             << " mild=" << nBLBLMild << ")"
             << " " << blblJunctionPoints_.size() << " junction pts captured"
             << endl;
    }

    // Detect generic layer-termination edges where at least one adjacent
    // patch requests boundary layers and at least one adjacent patch requests
    // zero layers. These edges are protected because unconstrained extrusion
    // across mixed layer/no-layer patch junctions creates invalid cells.
    // Fully general: no patch names, no face-count assumptions.
    // Works for manifold, non-manifold, open, and multi-region edges.
    {
        blNoBlEdges_.clear();
        blNoBlEdgePoints_.clear();
        blNoBlPointPatch_.clear();

        // Build global-to-boundary-point reverse map once, O(N)
        Map<label> globalToBP;
        forAll(bPoints, bpI)
            globalToBP.insert(bPoints[bpI], bpI);

        forAll(edges, eI)
        {
            // Collect all patch IDs adjacent to this edge.
            // Only BL + explicit TERMINATION edges are suppression edges.
            // Neutral patches (periodic/symmetry) do not trigger zero BL here.
            bool hasBL   = false;
            bool hasNoBL = false;
            label blPatchI = -1;
            label noBLPatchI = -1;

            forAllRow(edgeFaces, eI, efI)
            {
                const label faceI  = edgeFaces(eI, efI);
                const label patchI = boundaryFacePatches[faceI];
                if( patchI < 0 || patchI >= label(patchRole_.size()) )
                    continue;
                if( patchRole_[patchI] == 0 )
                {
                    hasBL = true;
                    blPatchI = patchI;
                }
                else if( patchRole_[patchI] == 1 )
                {
                    hasNoBL = true;
                    noBLPatchI = patchI;
                }
            }
            if( !hasBL || !hasNoBL ) continue;

            blNoBlEdges_.insert(eI);

            const edge& e = edges[eI];
            for( label ei = 0; ei < 2; ++ei )
            {
                const label gp = (ei == 0) ? e[0] : e[1];
                Map<label>::const_iterator it = globalToBP.find(gp);
                if( it == globalToBP.end() ) continue;
                const label bpI = it();
                if( blNoBlEdgePoints_.insert(bpI) )
                {
                    bool flowTermination = false;
                    if( noBLPatchI >= 0 && noBLPatchI < label(patchNames_.size()) )
                    {
                        const word& noBLName = patchNames_[noBLPatchI];
                        flowTermination =
                            (noBLName == "inlet" || noBLName == "outlet");
                    }
                    const scalar blNoBLTaper =
                        flowTermination ? scalar(0.25) : scalar(0.0);
                    zeroDistPoints_[bpI] = (blNoBLTaper < 0.01);
                    layerScale_[bpI] = blNoBLTaper;
                    blSuppressReason_[bpI] = 7;
                    zeroPts[bpI] = (blNoBLTaper < 0.01);
                }
                // Store BL-side patch for patch-constrained projection
                if( !blNoBlPointPatch_.found(bpI) )
                    blNoBlPointPatch_.insert(bpI, blPatchI);
            }
        }

        Info << "BL/no-BL transition edge detection: "
             << blNoBlEdges_.size() << " transition edges, "
             << blNoBlEdgePoints_.size() << " interface points locked" << endl;
    }

    // BL dropout reason atlas CSV dump.
    // Emits one row per suppressed/thin boundary point.
    // reasonCode: 0=none 1=gap 2=transEdge 3=tripleJunction 4=corner
    //             5=sharpEdge 6=blblSharp 7=blNoBLEdge
    {
        const meshSurfaceEngine& mseAtlas = surfaceEngine();
        const labelList& bPtsAtlas = mseAtlas.boundaryPoints();
        const pointField& ptsAtlas = mseAtlas.mesh().points();
        meshSurfacePartitioner mPartAtlas(mseAtlas);
        const VRWGraph& pPatchesAtlas = mPartAtlas.pointPatches();
        const PtrList<boundaryPatch>& bndAtlas = mesh_.boundaries();
        OFstream atlasOs("blDropoutReasonAtlas.csv");
        atlasOs << "bpI,meshPointI,x,y,z,layerScale,"
                << "reasonCode,reasonName,"
                << "rampRing,rampSeedReason,rampSeedName,"
                << "patches" << nl;
        const char* reasonNames[] = {
            "none","gap","transEdge","tripleJunction",
            "corner","sharpEdge","blblSharp","blNoBLEdge"
        };
        label nDumped = 0;
        forAll(bPtsAtlas, bpI)
        {
            const scalar ls =
                (layerScale_.size() > bpI) ? layerScale_[bpI] : scalar(1.0);
            const label rc =
                (blSuppressReason_.size() > bpI) ? blSuppressReason_[bpI] : 0;
            if( rc == 0 && ls >= 0.99 ) continue;
            const label meshPtI = bPtsAtlas[bpI];
            const point& pt =
                (meshPtI >= 0 && meshPtI < label(ptsAtlas.size())) ?
                ptsAtlas[meshPtI] : point(Zero);
            word patchStr("unknown");
            if( pPatchesAtlas.sizeOfRow(bpI) > 0 )
            {
                patchStr = "";
                labelHashSet seenP;
                forAllRow(pPatchesAtlas, bpI, pI)
                {
                    const label patchI = pPatchesAtlas(bpI, pI);
                    if( patchI < 0 || patchI >= label(bndAtlas.size()) ) continue;
                    if( seenP.found(patchI) ) continue;
                    seenP.insert(patchI);
                    if( patchStr.size() > 0 ) patchStr += "+";
                    patchStr += bndAtlas[patchI].patchName();
                }
                if( patchStr.size() == 0 ) patchStr = "unknown";
            }
            const label rcClamped = (rc >= 0 && rc <= 7) ? rc : 0;
            const label rampRing =
                (blRampRing_.size() > bpI) ? blRampRing_[bpI] : 0;
            const label rampSeedRc =
                (blRampSeedReason_.size() > bpI) ? blRampSeedReason_[bpI] : 0;
            const label rampSeedClamped =
                (rampSeedRc >= 0 && rampSeedRc <= 7) ? rampSeedRc : 0;
            atlasOs << bpI << ","
                    << meshPtI << ","
                    << pt.x() << ","
                    << pt.y() << ","
                    << pt.z() << ","
                    << ls << ","
                    << rcClamped << ","
                    << reasonNames[rcClamped] << ","
                    << rampRing << ","
                    << rampSeedClamped << ","
                    << reasonNames[rampSeedClamped] << ","
                    << patchStr << nl;
            ++nDumped;
        }
        Info << "BL dropout reason atlas: wrote " << nDumped
             << " suppressed/thin points to blDropoutReasonAtlas.csv" << endl;
    }


    // Topology-role VTK diagnostic dump
    // Roles: 0=SINGLE_PATCH 1=TWO_PATCH_EDGE 2=MULTI_PATCH_CORNER
    //        3=BL_NOBL_TRANSITION 4=BLBL_JUNCTION 5=ACUTE_CORNER
    {
        const meshSurfaceEngine& mseDiag = surfaceEngine();
        const labelList& bPtsDiag = mseDiag.boundaryPoints();
        const pointField& ptsDiag = mseDiag.mesh().points();
        meshSurfacePartitioner mPartDiag(mseDiag);
        const labelHashSet& cornersDiag = mPartDiag.corners();
        const labelHashSet& edgePtsDiag = mPartDiag.edgePoints();
        OFstream vtkFile("boundaryPointRoles.vtk");
        vtkFile << "# vtk DataFile Version 3.0\n";
        vtkFile << "Boundary point topology roles\n";
        vtkFile << "ASCII\n";
        vtkFile << "DATASET POLYDATA\n";
        vtkFile << "POINTS " << bPtsDiag.size() << " float\n";
        forAll(bPtsDiag, bpI)
        {
            const point& p = ptsDiag[bPtsDiag[bpI]];
            vtkFile << p.x() << ' ' << p.y() << ' ' << p.z() << '\n';
        }
        vtkFile << "\nVERTICES " << bPtsDiag.size()
                << " " << 2*bPtsDiag.size() << '\n';
        forAll(bPtsDiag, bpI)
            vtkFile << "1 " << bpI << '\n';
        vtkFile << "\nPOINT_DATA " << bPtsDiag.size() << '\n';
        vtkFile << "SCALARS role int 1\n";
        vtkFile << "LOOKUP_TABLE default\n";
        label nSingle=0, nEdge=0, nCorner=0, nBlNoBl=0, nBlBl=0, nAcute=0;
        forAll(bPtsDiag, bpI)
        {
            int role = 0;
            if( blblAcuteCornerPoints_.found(bpI) )      role = 5;
            else if( blblJunctionPoints_.found(bpI) )    role = 4;
            else if( blNoBlEdgePoints_.found(bpI) )      role = 3;
            else if( cornersDiag.found(bpI) )            role = 2;
            else if( edgePtsDiag.found(bpI) )            role = 1;
            vtkFile << role << '\n';
            switch(role)
            {
                case 0: ++nSingle; break;
                case 1: ++nEdge;   break;
                case 2: ++nCorner; break;
                case 3: ++nBlNoBl; break;
                case 4: ++nBlBl;   break;
                case 5: ++nAcute;  break;
            }
        }
        Info << "[Diag] boundaryPointRoles.vtk written" << endl;
        Info << "[Diag] Role counts:"
             << " SINGLE_PATCH=" << nSingle
             << " TWO_PATCH_EDGE=" << nEdge
             << " MULTI_PATCH_CORNER=" << nCorner
             << " BL_NOBL=" << nBlNoBl
             << " BLBL_JUNCTION=" << nBlBl
             << " ACUTE_CORNER=" << nAcute
             << endl;
    }

    // Ring 1: neighbors of zero points on BL patches -> 0.25
    boolList ring1(bPoints.size(), false);
    forAll(bPoints, bpI)
    {
        if( !zeroPts[bpI] ) continue;
        forAllRow(pointPoints, bpI, ppI)
        {
            const label nbpI = pointPoints(bpI, ppI);
            if( nbpI < 0 || nbpI >= label(bPoints.size()) ) continue;
            if( zeroPts[nbpI] ) continue;
            if( !boundaryPointIsBL[nbpI] ) continue;
            // Cross-patch: allow propagation but cap scale at 5% to prevent
            // bleed while avoiding abrupt BL termination / high skew.
            if( !boundaryPointsCanShareBLRamp(pPatches, isBLPatch, bpI, nbpI) ) continue;
            ring1[nbpI] = true;
            layerScale_[nbpI] = Foam::min(layerScale_[nbpI], layerScaleRing1_);
            if( blRampRing_[nbpI] == 0 )
            {
                blRampRing_[nbpI] = 1;
                blRampSeedReason_[nbpI] =
                    (blSuppressReason_[bpI] > 0) ?
                    blSuppressReason_[bpI] : blRampSeedReason_[bpI];
            }
        }
    }

    // Ring 2: neighbors of ring1 on BL patches -> 0.5
    boolList ring2(bPoints.size(), false);
    forAll(bPoints, bpI)
    {
        if( !ring1[bpI] ) continue;
        forAllRow(pointPoints, bpI, ppI)
        {
            const label nbpI = pointPoints(bpI, ppI);
            if( nbpI < 0 || nbpI >= label(bPoints.size()) ) continue;
            if( zeroPts[nbpI] || ring1[nbpI] ) continue;
            if( !boundaryPointIsBL[nbpI] ) continue;
            // Cross-patch: allow propagation but cap scale at 5% to prevent
            // bleed while avoiding abrupt BL termination / high skew.
            if( !boundaryPointsCanShareBLRamp(pPatches, isBLPatch, bpI, nbpI) ) continue;
            ring2[nbpI] = true;
            layerScale_[nbpI] = Foam::min(layerScale_[nbpI], layerScaleRing2_);
            if( blRampRing_[nbpI] == 0 )
            {
                blRampRing_[nbpI] = 2;
                blRampSeedReason_[nbpI] =
                    (blSuppressReason_[bpI] > 0) ?
                    blSuppressReason_[bpI] : blRampSeedReason_[bpI];
            }
        }
    }

    // Ring 3: neighbors of ring2 on BL patches -> 0.75
    boolList ring3(bPoints.size(), false);
    forAll(bPoints, bpI)
    {
        if( !ring2[bpI] ) continue;
        forAllRow(pointPoints, bpI, ppI)
        {
            const label nbpI = pointPoints(bpI, ppI);
            if( nbpI < 0 || nbpI >= label(bPoints.size()) ) continue;
            if( zeroPts[nbpI] || ring1[nbpI] || ring2[nbpI] ) continue;
            if( !boundaryPointIsBL[nbpI] ) continue;
            // Cross-patch: allow propagation but cap scale at 5% to prevent
            // bleed while avoiding abrupt BL termination / high skew.
            if( !boundaryPointsCanShareBLRamp(pPatches, isBLPatch, bpI, nbpI) ) continue;
            ring3[nbpI] = true;
            layerScale_[nbpI] = Foam::min(layerScale_[nbpI], layerScaleRing3_);
            if( blRampRing_[nbpI] == 0 )
            {
                blRampRing_[nbpI] = 3;
                blRampSeedReason_[nbpI] =
                    (blSuppressReason_[bpI] > 0) ?
                    blSuppressReason_[bpI] : blRampSeedReason_[bpI];
            }
        }
    }

    // Ring 4: neighbors of ring3 on BL patches -> 0.90
    boolList ring4(bPoints.size(), false);
    forAll(bPoints, bpI)
    {
        if( !ring3[bpI] ) continue;
        forAllRow(pointPoints, bpI, ppI)
        {
            const label nbpI = pointPoints(bpI, ppI);
            if( nbpI < 0 || nbpI >= label(bPoints.size()) ) continue;
            if( zeroPts[nbpI] || ring1[nbpI] || ring2[nbpI] || ring3[nbpI] ) continue;
            if( !boundaryPointIsBL[nbpI] ) continue;
            // Cross-patch: allow propagation but cap scale at 5% to prevent
            // bleed while avoiding abrupt BL termination / high skew.
            if( !boundaryPointsCanShareBLRamp(pPatches, isBLPatch, bpI, nbpI) ) continue;
            ring4[nbpI] = true;
            layerScale_[nbpI] = Foam::min(layerScale_[nbpI], layerScaleRing4_);
            if( blRampRing_[nbpI] == 0 )
            {
                blRampRing_[nbpI] = 4;
                blRampSeedReason_[nbpI] =
                    (blSuppressReason_[bpI] > 0) ?
                    blSuppressReason_[bpI] : blRampSeedReason_[bpI];
            }
        }
    }

    boolList ring5(bPoints.size(), false);
    forAll(bPoints, bpI)
    {
        if( !ring4[bpI] ) continue;
        forAllRow(pointPoints, bpI, ppI)
        {
            const label nbpI = pointPoints(bpI, ppI);
            if( nbpI < 0 || nbpI >= label(bPoints.size()) ) continue;
            if( zeroPts[nbpI] || ring1[nbpI] || ring2[nbpI] || ring3[nbpI] || ring4[nbpI] ) continue;
            if( !boundaryPointIsBL[nbpI] ) continue;
            // Cross-patch: allow propagation but cap scale at 5% to prevent
            // bleed while avoiding abrupt BL termination / high skew.
            if( !boundaryPointsCanShareBLRamp(pPatches, isBLPatch, bpI, nbpI) ) continue;
            ring5[nbpI] = true;
            layerScale_[nbpI] = Foam::min(layerScale_[nbpI], layerScaleRing5_);
            if( blRampRing_[nbpI] == 0 )
            {
                blRampRing_[nbpI] = 5;
                blRampSeedReason_[nbpI] =
                    (blSuppressReason_[bpI] > 0) ?
                    blSuppressReason_[bpI] : blRampSeedReason_[bpI];
            }
        }
    }
    boolList ring6(bPoints.size(), false);
    forAll(bPoints, bpI)
    {
        if( !ring5[bpI] ) continue;
        forAllRow(pointPoints, bpI, ppI)
        {
            const label nbpI = pointPoints(bpI, ppI);
            if( nbpI < 0 || nbpI >= label(bPoints.size()) ) continue;
            if( zeroPts[nbpI] || ring1[nbpI] || ring2[nbpI] || ring3[nbpI] || ring4[nbpI] || ring5[nbpI] ) continue;
            if( !boundaryPointIsBL[nbpI] ) continue;
            // Cross-patch: allow propagation but cap scale at 5% to prevent
            // bleed while avoiding abrupt BL termination / high skew.
            if( !boundaryPointsCanShareBLRamp(pPatches, isBLPatch, bpI, nbpI) ) continue;
            ring6[nbpI] = true;
            layerScale_[nbpI] = Foam::min(layerScale_[nbpI], layerScaleRing6_);
            if( blRampRing_[nbpI] == 0 )
            {
                blRampRing_[nbpI] = 6;
                blRampSeedReason_[nbpI] =
                    (blSuppressReason_[bpI] > 0) ?
                    blSuppressReason_[bpI] : blRampSeedReason_[bpI];
            }
        }
    }
    label nZero = 0, nRing1 = 0, nRing2 = 0, nRing3 = 0, nRing4 = 0, nRing5 = 0, nRing6 = 0;
    forAll(bPoints, bpI)
    {
        if( zeroPts[bpI] ) ++nZero;
        else if( ring1[bpI] ) ++nRing1;
        else if( ring2[bpI] ) ++nRing2;
        else if( ring3[bpI] ) ++nRing3;
        else if( ring4[bpI] ) ++nRing4;
        else if( ring5[bpI] ) ++nRing5;
        else if( ring6[bpI] ) ++nRing6;
    }
    Info << "BL layerScale ramp: zero=" << nZero
         << " ring1=" << nRing1
         << " ring2=" << nRing2
         << " ring3=" << nRing3
         << " ring4=" << nRing4
         << " ring5=" << nRing5
         << " ring6=" << nRing6
         << endl;
    // Acute corner local taper
    // Seeded from blblAcuteCornerPoints_ -- stronger suppression than general ramp
    // Propagates along BL wall patches only
    if( blblAcuteCornerPoints_.size() > 0 )
    {
        // Apply ring0 scale to exact acute corner points
        forAllConstIter(labelHashSet, blblAcuteCornerPoints_, it)
        {
            const label bpI = it.key();
            if( bpI >= 0 && bpI < label(bPoints.size()) )
                layerScale_[bpI] = Foam::min(layerScale_[bpI], acuteCornerRing0_);
        }

        // Build acute corner seed set for ring propagation
        boolList acuteZero(bPoints.size(), false);
        forAllConstIter(labelHashSet, blblAcuteCornerPoints_, it)
        {
            const label bpI = it.key();
            if( bpI >= 0 && bpI < label(bPoints.size()) )
                acuteZero[bpI] = true;
        }

        // Ring 1 from acute corners
        boolList acuteRing1(bPoints.size(), false);
        forAll(bPoints, bpI)
        {
            if( !acuteZero[bpI] ) continue;
            forAllRow(pointPoints, bpI, ppI)
            {
                const label nbpI = pointPoints(bpI, ppI);
                if( nbpI < 0 || nbpI >= label(bPoints.size()) ) continue;
                if( acuteZero[nbpI] ) continue;
                if( !boundaryPointIsBL[nbpI] ) continue;
            // Cross-patch: allow propagation but cap scale at 5% to prevent
            // bleed while avoiding abrupt BL termination / high skew.
            const bool crossPatch =
                !boundaryPointsCanShareBLRamp(pPatches, isBLPatch, bpI, nbpI);
                if( crossPatch ) continue;
                acuteRing1[nbpI] = true;
                layerScale_[nbpI] = Foam::min(layerScale_[nbpI], acuteCornerRing1_);
            }
        }

        // Ring 2
        boolList acuteRing2(bPoints.size(), false);
        forAll(bPoints, bpI)
        {
            if( !acuteRing1[bpI] ) continue;
            forAllRow(pointPoints, bpI, ppI)
            {
                const label nbpI = pointPoints(bpI, ppI);
                if( nbpI < 0 || nbpI >= label(bPoints.size()) ) continue;
                if( acuteZero[nbpI] || acuteRing1[nbpI] ) continue;
                if( !boundaryPointIsBL[nbpI] ) continue;
            // Cross-patch: allow propagation but cap scale at 5% to prevent
            // bleed while avoiding abrupt BL termination / high skew.
            const bool crossPatch =
                !boundaryPointsCanShareBLRamp(pPatches, isBLPatch, bpI, nbpI);
                if( crossPatch ) continue;
                acuteRing2[nbpI] = true;
                layerScale_[nbpI] = Foam::min(layerScale_[nbpI], acuteCornerRing2_);
            }
        }

        // Ring 3
        boolList acuteRing3(bPoints.size(), false);
        forAll(bPoints, bpI)
        {
            if( !acuteRing2[bpI] ) continue;
            forAllRow(pointPoints, bpI, ppI)
            {
                const label nbpI = pointPoints(bpI, ppI);
                if( nbpI < 0 || nbpI >= label(bPoints.size()) ) continue;
                if( acuteZero[nbpI] || acuteRing1[nbpI] || acuteRing2[nbpI] ) continue;
                if( !boundaryPointIsBL[nbpI] ) continue;
            // Cross-patch: allow propagation but cap scale at 5% to prevent
            // bleed while avoiding abrupt BL termination / high skew.
            const bool crossPatch =
                !boundaryPointsCanShareBLRamp(pPatches, isBLPatch, bpI, nbpI);
                if( crossPatch ) continue;
                acuteRing3[nbpI] = true;
                layerScale_[nbpI] = Foam::min(layerScale_[nbpI], acuteCornerRing3_);
            }
        }

        // Ring 4
        boolList acuteRing4(bPoints.size(), false);
        forAll(bPoints, bpI)
        {
            if( !acuteRing3[bpI] ) continue;
            forAllRow(pointPoints, bpI, ppI)
            {
                const label nbpI = pointPoints(bpI, ppI);
                if( nbpI < 0 || nbpI >= label(bPoints.size()) ) continue;
                if( acuteZero[nbpI] || acuteRing1[nbpI] || acuteRing2[nbpI] || acuteRing3[nbpI] ) continue;
                if( !boundaryPointIsBL[nbpI] ) continue;
            // Cross-patch: allow propagation but cap scale at 5% to prevent
            // bleed while avoiding abrupt BL termination / high skew.
            const bool crossPatch =
                !boundaryPointsCanShareBLRamp(pPatches, isBLPatch, bpI, nbpI);
                if( crossPatch ) continue;
                acuteRing4[nbpI] = true;
                layerScale_[nbpI] = Foam::min(layerScale_[nbpI], acuteCornerRing4_);
            }
        }

        // Ring 5
        boolList acuteRing5(bPoints.size(), false);
        forAll(bPoints, bpI)
        {
            if( !acuteRing4[bpI] ) continue;
            forAllRow(pointPoints, bpI, ppI)
            {
                const label nbpI = pointPoints(bpI, ppI);
                if( nbpI < 0 || nbpI >= label(bPoints.size()) ) continue;
                if( acuteZero[nbpI] || acuteRing1[nbpI] || acuteRing2[nbpI] || acuteRing3[nbpI] || acuteRing4[nbpI] ) continue;
                if( !boundaryPointIsBL[nbpI] ) continue;
            // Cross-patch: allow propagation but cap scale at 5% to prevent
            // bleed while avoiding abrupt BL termination / high skew.
            const bool crossPatch =
                !boundaryPointsCanShareBLRamp(pPatches, isBLPatch, bpI, nbpI);
                if( crossPatch ) continue;
                acuteRing5[nbpI] = true;
                layerScale_[nbpI] = Foam::min(layerScale_[nbpI], acuteCornerRing5_);
            }
        }

        // Ring 6
        boolList acuteRing6(bPoints.size(), false);
        forAll(bPoints, bpI)
        {
            if( !acuteRing5[bpI] ) continue;
            forAllRow(pointPoints, bpI, ppI)
            {
                const label nbpI = pointPoints(bpI, ppI);
                if( nbpI < 0 || nbpI >= label(bPoints.size()) ) continue;
                if( acuteZero[nbpI] || acuteRing1[nbpI] || acuteRing2[nbpI] || acuteRing3[nbpI] || acuteRing4[nbpI] || acuteRing5[nbpI] ) continue;
                if( !boundaryPointIsBL[nbpI] ) continue;
            // Cross-patch: allow propagation but cap scale at 5% to prevent
            // bleed while avoiding abrupt BL termination / high skew.
            const bool crossPatch =
                !boundaryPointsCanShareBLRamp(pPatches, isBLPatch, bpI, nbpI);
                if( crossPatch ) continue;
                acuteRing6[nbpI] = true;
                layerScale_[nbpI] = Foam::min(layerScale_[nbpI], acuteCornerRing6_);
            }
        }

        label nAR0=0, nAR1=0, nAR2=0, nAR3=0, nAR4=0, nAR5=0, nAR6=0;
        forAll(bPoints, bpI)
        {
            if( acuteZero[bpI] ) ++nAR0;
            else if( acuteRing1[bpI] ) ++nAR1;
            else if( acuteRing2[bpI] ) ++nAR2;
            else if( acuteRing3[bpI] ) ++nAR3;
            else if( acuteRing4[bpI] ) ++nAR4;
            else if( acuteRing5[bpI] ) ++nAR5;
            else if( acuteRing6[bpI] ) ++nAR6;
        }
        Info << "Acute corner taper: ring0=" << nAR0
             << " ring1=" << nAR1
             << " ring2=" << nAR2
             << " ring3=" << nAR3
             << " ring4=" << nAR4
             << " ring5=" << nAR5
             << " ring6=" << nAR6
             << endl;
    }

    // Virtual topology: zero/taper layerScale_ at acute triple-junctions.
    // blblAcuteCornerPoints_ and layerScale_ fully built above.
    // Must run before return -- vertex creation reads layerScale_.
    applyVirtualTopologyExclusion();

    // Gap face-ring: taper layerScale_ around gap contact points
    applyGapFaceRingExclusion();

    Info << "terminateLayersAtConcaveEdges: marked "
         << nTransitionEdges << " BL-transition edges." << endl;
}

void boundaryLayers::applyVirtualTopologyExclusion() const
{
    if( !virtualTopologyExclusion_ )
        return;

    if( blblAcuteCornerPoints_.size() == 0 )
        return;

    const meshSurfaceEngine& mse = surfaceEngine();

    const labelList& bPoints = mse.boundaryPoints();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const VRWGraph& pointFaces = mse.pointFaces();

    const label nBP = bPoints.size();
    const label nBF = bFaces.size();

    if( layerScale_.size() != nBP )
        return;

    // Reverse map: mesh point label -> boundary point label.
    labelList meshToBnd(mesh_.points().size(), -1);
    forAll(bPoints, bpI)
        meshToBnd[bPoints[bpI]] = bpI;

    // faceRing: -1=untouched, 0=seed, 1=ring1, 2=ring2
    labelList faceRing(nBF, -1);

    // Ring 0: all boundary faces touching any acute corner seed point.
    forAllConstIter(labelHashSet, blblAcuteCornerPoints_, it)
    {
        const label bpI = it.key();
        if( bpI < 0 || bpI >= nBP ) continue;
        forAllRow(pointFaces, bpI, pfI)
        {
            const label bfI = pointFaces(bpI, pfI);
            if( bfI < 0 || bfI >= nBF ) continue;
            faceRing[bfI] = 0;
        }
    }

    // Build ring0 point set.
    boolList ring0pt(nBP, false);
    forAll(faceRing, bfI)
    {
        if( faceRing[bfI] != 0 ) continue;
        const face& f = bFaces[bfI];
        forAll(f, pI)
        {
            const label meshPtI = f[pI];
            if( meshPtI < 0 || meshPtI >= label(meshToBnd.size()) ) continue;
            const label bpI = meshToBnd[meshPtI];
            if( bpI >= 0 && bpI < nBP ) ring0pt[bpI] = true;
        }
    }

    // Ring 1: all faces touching ring0 points.
    forAll(ring0pt, bpI)
    {
        if( !ring0pt[bpI] ) continue;
        forAllRow(pointFaces, bpI, pfI)
        {
            const label bfI = pointFaces(bpI, pfI);
            if( bfI < 0 || bfI >= nBF ) continue;
            if( faceRing[bfI] >= 0 ) continue;
            faceRing[bfI] = 1;
        }
    }

    // Build ring1 point set.
    boolList ring1pt(nBP, false);
    forAll(faceRing, bfI)
    {
        if( faceRing[bfI] != 1 ) continue;
        const face& f = bFaces[bfI];
        forAll(f, pI)
        {
            const label meshPtI = f[pI];
            if( meshPtI < 0 || meshPtI >= label(meshToBnd.size()) ) continue;
            const label bpI = meshToBnd[meshPtI];
            if( bpI >= 0 && bpI < nBP ) ring1pt[bpI] = true;
        }
    }

    // Ring 2: all faces touching ring1 points.
    forAll(ring1pt, bpI)
    {
        if( !ring1pt[bpI] ) continue;
        forAllRow(pointFaces, bpI, pfI)
        {
            const label bfI = pointFaces(bpI, pfI);
            if( bfI < 0 || bfI >= nBF ) continue;
            if( faceRing[bfI] >= 0 ) continue;
            faceRing[bfI] = 2;
        }
    }

    // Apply layerScale_ -- smallest ring wins for shared points.
    forAll(faceRing, bfI)
    {
        const label ring = faceRing[bfI];
        if( ring < 0 ) continue;
        scalar scale = 1.0;
        if(      ring == 0 ) scale = virtualTopoRing0_;
        else if( ring == 1 ) scale = virtualTopoRing1_;
        else if( ring == 2 ) scale = virtualTopoRing2_;
        const face& f = bFaces[bfI];
        forAll(f, pI)
        {
            const label meshPtI = f[pI];
            if( meshPtI < 0 || meshPtI >= label(meshToBnd.size()) ) continue;
            const label bpI = meshToBnd[meshPtI];
            if( bpI < 0 || bpI >= nBP ) continue;
            layerScale_[bpI] = Foam::min(layerScale_[bpI], scale);
        }
    }

    // Export face ring for refineBoundaryLayers -- BL/BL ramp skips VT-handled faces
    vtFaceRing_ = faceRing;

    // Diagnostics
    label nF0=0, nF1=0, nF2=0;
    labelList ptRing(nBP, -1);
    forAll(faceRing, bfI)
    {
        const label ring = faceRing[bfI];
        if(      ring == 0 ) ++nF0;
        else if( ring == 1 ) ++nF1;
        else if( ring == 2 ) ++nF2;
        if( ring < 0 ) continue;
        const face& f = bFaces[bfI];
        forAll(f, pI)
        {
            const label meshPtI = f[pI];
            if( meshPtI < 0 || meshPtI >= label(meshToBnd.size()) ) continue;
            const label bpI = meshToBnd[meshPtI];
            if( bpI < 0 || bpI >= nBP ) continue;
            if( ptRing[bpI] < 0 || ring < ptRing[bpI] ) ptRing[bpI] = ring;
        }
    }
    label nP0=0, nP1=0, nP2=0;
    forAll(ptRing, bpI)
    {
        if(      ptRing[bpI] == 0 ) ++nP0;
        else if( ptRing[bpI] == 1 ) ++nP1;
        else if( ptRing[bpI] == 2 ) ++nP2;
    }

    Info << "VirtualTopology face-ring exclusion:"
         << " faces(r0=" << nF0 << " r1=" << nF1 << " r2=" << nF2 << ")"
         << " pts(r0=" << nP0 << " r1=" << nP1 << " r2=" << nP2 << ")"
         << " scale=(" << virtualTopoRing0_
         << " " << virtualTopoRing1_
         << " " << virtualTopoRing2_ << ")"
         << endl;
}

void boundaryLayers::applyGapFaceRingExclusion() const
{
    const bool useTriple =
        tripleJunctionFaceRingExclusion_ && tripleJunctionPoints_.size() > 0;
    if( !gapFaceRingExclusion_ && !useTriple )
        return;
    if( gapPoints_.size() == 0 && !useTriple )
        return;

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& bPoints     = mse.boundaryPoints();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const VRWGraph& pointFaces   = mse.pointFaces();
    const labelList& facePatch   = mse.boundaryFacePatches();

    const label nBP = bPoints.size();
    const label nBF = bFaces.size();

    if( layerScale_.size() != nBP ) return;

    // Reverse map: mesh point -> boundary point
    labelList meshToBnd(mesh_.points().size(), -1);
    forAll(bPoints, bpI)
        meshToBnd[bPoints[bpI]] = bpI;

    // faceRing: -1=untouched, 0=seed, 1=ring1, 2=ring2
    labelList faceRing(nBF, -1);

    // Ring 0: boundary faces touching gap points and/or triple-junction points
    auto seedRing0 = [&](const labelHashSet& seedPts)
    {
        forAllConstIter(labelHashSet, seedPts, it)
        {
            const label meshPtI = it.key();
            if( meshPtI < 0 || meshPtI >= label(meshToBnd.size()) ) continue;
            const label bpI = meshToBnd[meshPtI];
            if( bpI < 0 || bpI >= nBP ) continue;
            forAllRow(pointFaces, bpI, pfI)
            {
                const label bfI = pointFaces(bpI, pfI);
                if( bfI < 0 || bfI >= nBF ) continue;
                faceRing[bfI] = 0;
            }
        }
    };
    if( gapFaceRingExclusion_ )
        seedRing0(gapPoints_);
    if( useTriple )
    {
        // Patch-restricted seeding: only suppress faces on the designated
        // suppress-side patches (hub/shroud). Blade and periodic faces
        // adjacent to the triple junction must NOT be suppressed --
        // they need BL for correct turbomachinery flow resolution.
        if( tripleJunctionSuppressPatches_.size() > 0 )
        {
            forAllConstIter(labelHashSet, tripleJunctionPoints_, it)
            {
                const label meshPtI = it.key();
                if( meshPtI < 0 || meshPtI >= label(meshToBnd.size()) ) continue;
                const label bpI = meshToBnd[meshPtI];
                if( bpI < 0 || bpI >= nBP ) continue;
                forAllRow(pointFaces, bpI, pfI)
                {
                    const label bfI = pointFaces(bpI, pfI);
                    if( bfI < 0 || bfI >= nBF ) continue;
                    // Only suppress faces on the suppress-side patches
                    if( !tripleJunctionSuppressPatches_.found(facePatch[bfI]) )
                        continue;
                    faceRing[bfI] = 0;
                }
            }
        }
        else
        {
            Info << "Triple-junction face-ring exclusion: "
                 << "no suppress patch filter defined -- skipping" << endl;
        }
    }

    // Build protected-side taper mask: blade/periodic faces at triple
    // junctions NOT on suppress-side patches. Taper only -- no topology suppression.
    boolList protectedTripleRing0Face(nBF, false);
    if( useTriple && tripleJunctionSuppressPatches_.size() > 0 )
    {
        forAllConstIter(labelHashSet, tripleJunctionPoints_, it)
        {
            const label meshPtI = it.key();
            if( meshPtI < 0 || meshPtI >= label(meshToBnd.size()) ) continue;
            const label bpI = meshToBnd[meshPtI];
            if( bpI < 0 || bpI >= nBP ) continue;
            forAllRow(pointFaces, bpI, pfI)
            {
                const label bfI = pointFaces(bpI, pfI);
                if( bfI < 0 || bfI >= nBF ) continue;
                // Only taper faces NOT on suppress-side patches
                if( tripleJunctionSuppressPatches_.found(facePatch[bfI]) )
                    continue;
                // Don't override an existing hard suppress
                if( faceRing[bfI] == 0 )
                    continue;
                protectedTripleRing0Face[bfI] = true;
            }
        }
    }

    // Build ring0 point set
    boolList ring0pt(nBP, false);
    forAll(faceRing, bfI)
    {
        if( faceRing[bfI] != 0 ) continue;
        const face& f = bFaces[bfI];
        forAll(f, pI)
        {
            const label meshPtI = f[pI];
            if( meshPtI < 0 || meshPtI >= label(meshToBnd.size()) ) continue;
            const label bpI = meshToBnd[meshPtI];
            if( bpI >= 0 && bpI < nBP ) ring0pt[bpI] = true;
        }
    }

    // Ring 1
    forAll(ring0pt, bpI)
    {
        if( !ring0pt[bpI] ) continue;
        forAllRow(pointFaces, bpI, pfI)
        {
            const label bfI = pointFaces(bpI, pfI);
            if( bfI < 0 || bfI >= nBF ) continue;
            if( faceRing[bfI] >= 0 ) continue;
            faceRing[bfI] = 1;
        }
    }

    // Build ring1 point set
    boolList ring1pt(nBP, false);
    forAll(faceRing, bfI)
    {
        if( faceRing[bfI] != 1 ) continue;
        const face& f = bFaces[bfI];
        forAll(f, pI)
        {
            const label meshPtI = f[pI];
            if( meshPtI < 0 || meshPtI >= label(meshToBnd.size()) ) continue;
            const label bpI = meshToBnd[meshPtI];
            if( bpI >= 0 && bpI < nBP ) ring1pt[bpI] = true;
        }
    }

    // Ring 2
    forAll(ring1pt, bpI)
    {
        if( !ring1pt[bpI] ) continue;
        forAllRow(pointFaces, bpI, pfI)
        {
            const label bfI = pointFaces(bpI, pfI);
            if( bfI < 0 || bfI >= nBF ) continue;
            if( faceRing[bfI] >= 0 ) continue;
            faceRing[bfI] = 2;
        }
    }

    // Ring 3 face tracking
    boolList ring2pt(nBP, false);
    forAll(faceRing, bfI)
    {
        if( faceRing[bfI] != 2 ) continue;
        const face& f = bFaces[bfI];
        forAll(f, pI)
        {
            const label meshPtI = f[pI];
            if( meshPtI < 0 || meshPtI >= label(meshToBnd.size()) ) continue;
            const label bpI = meshToBnd[meshPtI];
            if( bpI >= 0 && bpI < nBP ) ring2pt[bpI] = true;
        }
    }

    forAll(ring2pt, bpI)
    {
        if( !ring2pt[bpI] ) continue;
        forAllRow(pointFaces, bpI, pfI)
        {
            const label bfI = pointFaces(bpI, pfI);
            if( bfI < 0 || bfI >= nBF ) continue;
            if( faceRing[bfI] >= 0 ) continue;
            faceRing[bfI] = 3;
        }
    }

    // Populate face-level suppression mask.
    // Ring0: always suppress topology (gap too tight for any BL).
    // Ring1 on loser-side patches: suppress topology to give winner-side BL room.
    // createNewFacesAndCells will skip suppressed faces entirely.
    suppressLayerAtBndFace_.setSize(nBF, false);
    const bool applyLoserRing1 =
        gapLoserRing1Suppress_ && gapLoserPatches_.size() > 0;
    label nLoserRing1Suppressed = 0;
    forAll(faceRing, bfI)
    {
        if( faceRing[bfI] == 0 )
        {
            suppressLayerAtBndFace_[bfI] = true;
        }
        else if( faceRing[bfI] == 1 && applyLoserRing1 )
        {
            // Only suppress ring1 on the loser side of the gap conflict.
            // Winner side keeps BL topology through the transition.
            if( gapLoserPatches_.found(facePatch[bfI]) )
            {
                suppressLayerAtBndFace_[bfI] = true;
                ++nLoserRing1Suppressed;
            }
        }
    }
    if( nLoserRing1Suppressed > 0 )
        Info << "Gap conflict arbitration: ring1 topology suppressed on loser side: "
             << nLoserRing1Suppressed << " faces" << endl;

    // Gap layer-count cap metadata is passed to refineBoundaryLayers via
    // gapPoints_ (mesh point labels) + gapLoserPatchNames() + ring max knobs.
    // refineBoundaryLayers applies the cap locally at split-edge generation time
    // using stable mesh-point-label + patch-name identification (Option B).
    if( gapLoserPatches_.size() > 0 )
        Info << "Gap layer-count caps: loser-side patches="
             << gapLoserPatches_.size()
             << " gap action points=" << gapPoints_.size()
             << " (ring1max=" << gapLoserRing1MaxLayers_
             << " ring2max=" << gapLoserRing2MaxLayers_ << ")" << endl;

    // Apply layerScale_ -- scales configurable via meshDict
    forAll(faceRing, bfI)
    {
        const label ring = faceRing[bfI];
        if( ring < 0 ) continue;
        scalar scale = 1.0;
        if(      ring == 0 ) scale = gapFaceRing0Scale_;
        else if( ring == 1 ) scale = gapFaceRing1Scale_;
        else if( ring == 2 ) scale = gapFaceRing2Scale_;
        else if( ring == 3 ) scale = gapFaceRing3Scale_;
        const face& f = bFaces[bfI];
        forAll(f, pI)
        {
            const label meshPtI = f[pI];
            if( meshPtI < 0 || meshPtI >= label(meshToBnd.size()) ) continue;
            const label bpI = meshToBnd[meshPtI];
            if( bpI < 0 || bpI >= nBP ) continue;
            layerScale_[bpI] = Foam::min(layerScale_[bpI], scale);
        }
    }

    // Apply protected-side triple-junction taper.
    // These faces are not topology-suppressed; only local BL height is reduced.
    label nProtectedTripleR0 = 0;
    forAll(protectedTripleRing0Face, bfI)
    {
        if( !protectedTripleRing0Face[bfI] )
            continue;

        ++nProtectedTripleR0;

        const face& f = bFaces[bfI];
        forAll(f, pI)
        {
            const label meshPtI = f[pI];
            if( meshPtI < 0 || meshPtI >= label(meshToBnd.size()) ) continue;

            const label bpI = meshToBnd[meshPtI];
            if( bpI < 0 || bpI >= nBP ) continue;

            layerScale_[bpI] =
                Foam::min(layerScale_[bpI], tripleJunctionProtectedRing0Scale_);
        }
    }

    // Diagnostics
    label nF0=0, nF1=0, nF2=0, nF3=0;
    forAll(faceRing, bfI)
    {
        if(      faceRing[bfI] == 0 ) ++nF0;
        else if( faceRing[bfI] == 1 ) ++nF1;
        else if( faceRing[bfI] == 2 ) ++nF2;
        else if( faceRing[bfI] == 3 ) ++nF3;
    }
    Info << "Gap face-ring exclusion:"
         << " faces(r0=" << nF0 << " r1=" << nF1
         << " r2=" << nF2 << " r3=" << nF3 << ")"
         << " scale=(" << gapFaceRing0Scale_
         << " " << gapFaceRing1Scale_
         << " " << gapFaceRing2Scale_
         << " " << gapFaceRing3Scale_ << ")" << endl;

    if( nProtectedTripleR0 > 0 )
        Info << "Triple-junction protected-side taper:"
             << " faces=" << nProtectedTripleR0
             << " scale=" << tripleJunctionProtectedRing0Scale_
             << endl;
}

void boundaryLayers::reportBLTransitionSeeds() const
{
    Info << "BL transition seed summary:" << nl
         << "  gap points:                    " << gapPoints_.size() << nl
         << "  triple-junction points:        " << tripleJunctionPoints_.size() << nl
         << "  tripleJunctionFaceRingExclusion: "
         << (tripleJunctionFaceRingExclusion_ ? "true" : "false")
         << endl;
}

void boundaryLayers::buildBLTransitionPlan() const
{
    // Diagnostic only -- computes and reports ring face counts by patch.
    // No topology mutation.
    if( gapPoints_.size() == 0 && tripleJunctionPoints_.size() == 0 )
    {
        Info << "BL transition planner: no seeds, skipping" << endl;
        return;
    }

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& bPoints        = mse.boundaryPoints();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const VRWGraph& pointFaces      = mse.pointFaces();
    const labelList& facePatch      = mse.boundaryFacePatches();

    const label nBP = bPoints.size();
    const label nBF = bFaces.size();

    labelList meshToBnd(mesh_.points().size(), -1);
    forAll(bPoints, bpI)
        meshToBnd[bPoints[bpI]] = bpI;

    labelList faceRing(nBF, -1);

    auto seedRing0 = [&](const labelHashSet& seedPts)
    {
        forAllConstIter(labelHashSet, seedPts, it)
        {
            const label meshPtI = it.key();
            if( meshPtI < 0 || meshPtI >= label(meshToBnd.size()) ) continue;
            const label bpI = meshToBnd[meshPtI];
            if( bpI < 0 || bpI >= nBP ) continue;
            forAllRow(pointFaces, bpI, pfI)
            {
                const label bfI = pointFaces(bpI, pfI);
                if( bfI < 0 || bfI >= nBF ) continue;
                faceRing[bfI] = 0;
            }
        }
    };
    seedRing0(gapPoints_);
    seedRing0(tripleJunctionPoints_);

    // Ring 1
    boolList ring0pt(nBP, false);
    forAll(faceRing, bfI)
    {
        if( faceRing[bfI] != 0 ) continue;
        const face& f = bFaces[bfI];
        forAll(f, pI)
        {
            const label meshPtI = f[pI];
            if( meshPtI < 0 || meshPtI >= label(meshToBnd.size()) ) continue;
            const label bpI = meshToBnd[meshPtI];
            if( bpI >= 0 && bpI < nBP ) ring0pt[bpI] = true;
        }
    }
    forAll(ring0pt, bpI)
    {
        if( !ring0pt[bpI] ) continue;
        forAllRow(pointFaces, bpI, pfI)
        {
            const label bfI = pointFaces(bpI, pfI);
            if( bfI < 0 || bfI >= nBF ) continue;
            if( faceRing[bfI] >= 0 ) continue;
            faceRing[bfI] = 1;
        }
    }

    // Ring 2
    boolList ring1pt(nBP, false);
    forAll(faceRing, bfI)
    {
        if( faceRing[bfI] != 1 ) continue;
        const face& f = bFaces[bfI];
        forAll(f, pI)
        {
            const label meshPtI = f[pI];
            if( meshPtI < 0 || meshPtI >= label(meshToBnd.size()) ) continue;
            const label bpI = meshToBnd[meshPtI];
            if( bpI >= 0 && bpI < nBP ) ring1pt[bpI] = true;
        }
    }
    forAll(ring1pt, bpI)
    {
        if( !ring1pt[bpI] ) continue;
        forAllRow(pointFaces, bpI, pfI)
        {
            const label bfI = pointFaces(bpI, pfI);
            if( bfI < 0 || bfI >= nBF ) continue;
            if( faceRing[bfI] >= 0 ) continue;
            faceRing[bfI] = 2;
        }
    }

    // Count by ring and patch
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    const label nPatches = boundaries.size();
    List<label> patchR0(nPatches, 0);
    List<label> patchR1(nPatches, 0);
    List<label> patchR2(nPatches, 0);
    label nF0=0, nF1=0, nF2=0;

    forAll(faceRing, bfI)
    {
        const label ring = faceRing[bfI];
        const label pI   = facePatch[bfI];
        if( ring == 0 )
        {
            ++nF0;
            if( pI >= 0 && pI < nPatches ) ++patchR0[pI];
        }
        else if( ring == 1 )
        {
            ++nF1;
            if( pI >= 0 && pI < nPatches ) ++patchR1[pI];
        }
        else if( ring == 2 )
        {
            ++nF2;
            if( pI >= 0 && pI < nPatches ) ++patchR2[pI];
        }
    }

    Info << "BL transition planner (diagnostic -- no mutation):" << nl
         << "  ring0 suppress faces: " << nF0 << nl
         << "  ring1 cap-to-1 faces: " << nF1 << nl
         << "  ring2 cap-to-2 faces: " << nF2 << nl
         << "  per-patch breakdown:" << endl;
    forAll(boundaries, pI)
    {
        if( patchR0[pI] + patchR1[pI] + patchR2[pI] == 0 ) continue;
        Info << "    " << boundaries[pI].patchName()
             << ": r0=" << patchR0[pI]
             << " r1=" << patchR1[pI]
             << " r2=" << patchR2[pI] << nl;
    }
    Info << "  (tripleJunctionFaceRingExclusion="
         << (tripleJunctionFaceRingExclusion_ ? "true" : "false")
         << " -- topology unchanged)" << endl;
}

void boundaryLayers::reportBLPlanningPerPatch() const
{
    const meshSurfaceEngine& mse = surfaceEngine();

    const labelList& bPoints        = mse.boundaryPoints();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& facePatch      = mse.boundaryFacePatches();

    const label nBP = bPoints.size();

    // Reverse map mesh point -> boundary point
    labelList meshToBnd(mesh_.points().size(), -1);
    forAll(bPoints, bpI)
    {
        const label meshPtI = bPoints[bpI];
        if( meshPtI >= 0 && meshPtI < label(meshToBnd.size()) )
            meshToBnd[meshPtI] = bpI;
    }

    // Build gap and triple-junction boundary-point sets.
    // Use gapZonePoints_ (symmetric, both sides) for transition classification
    // so the full geometric danger zone is visible to the transition detector.
    // Fall back to gapPoints_ if no zone points available (legacy behavior).
    boolList isGapPoint(nBP, false);
    const labelHashSet& zoneSet =
        gapZonePoints_.size() > 0 ? gapZonePoints_ : gapPoints_;
    forAllConstIter(labelHashSet, zoneSet, it)
    {
        const label meshPtI = it.key();
        if( meshPtI < 0 || meshPtI >= label(meshToBnd.size()) ) continue;

        const label bpI = meshToBnd[meshPtI];
        if( bpI >= 0 && bpI < nBP )
            isGapPoint[bpI] = true;
    }

    boolList isTriplePoint(nBP, false);
    forAllConstIter(labelHashSet, tripleJunctionPoints_, it)
    {
        const label meshPtI = it.key();
        if( meshPtI < 0 || meshPtI >= label(meshToBnd.size()) ) continue;

        const label bpI = meshToBnd[meshPtI];
        if( bpI >= 0 && bpI < nBP )
            isTriplePoint[bpI] = true;
    }

    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    const label nPatches = boundaries.size();

    List<label> nPtTotal(nPatches, 0);
    List<label> nMultiPatchPts(nPatches, 0);
    List<label> nScale0(nPatches, 0);
    List<label> nScale002(nPatches, 0);
    List<label> nScale005(nPatches, 0);
    List<label> nScaleLt1(nPatches, 0);
    List<label> nGapSup(nPatches, 0);
    List<label> nTripleSup(nPatches, 0);

    List<label> nFaceTotal(nPatches, 0);
    List<label> nFaceSup(nPatches, 0);

    // Multi-patch ownership is exactly what we are debugging.
    // Count each boundary point against every patch it belongs to, not
    // only pPatches(bpI,0), otherwise blade/hub/periodic junctions are hidden.
    const meshSurfacePartitioner mPart(mse);
    const VRWGraph& pPatches = mPart.pointPatches();

    forAll(bPoints, bpI)
    {
        const scalar sc =
            layerScale_.size() > bpI ? layerScale_[bpI] : scalar(1.0);

        const bool isMultiPatch = pPatches.sizeOfRow(bpI) > 1;

        forAllRow(pPatches, bpI, ppi)
        {
            const label pI = pPatches(bpI, ppi);
            if( pI < 0 || pI >= nPatches ) continue;

            ++nPtTotal[pI];

            if( isMultiPatch )
                ++nMultiPatchPts[pI];

            if( sc <= 0.0 )       ++nScale0[pI];
            if( sc <= scalar(0.02) ) ++nScale002[pI];
            if( sc <= scalar(0.05) ) ++nScale005[pI];
            if( sc <  scalar(1.0) )  ++nScaleLt1[pI];

            if( isGapPoint[bpI] )
                ++nGapSup[pI];

            if( isTriplePoint[bpI] )
                ++nTripleSup[pI];
        }
    }

    forAll(bFaces, bfI)
    {
        if( bfI < 0 || bfI >= label(facePatch.size()) ) continue;

        const label pI = facePatch[bfI];
        if( pI < 0 || pI >= nPatches ) continue;

        ++nFaceTotal[pI];

        if( suppressLayerAtBndFace_.size() > bfI
         && suppressLayerAtBndFace_[bfI] )
        {
            ++nFaceSup[pI];
        }
    }

    Info << "BL planning per-patch audit:" << nl;

    forAll(boundaries, pI)
    {
        if( nPtTotal[pI] == 0 && nFaceTotal[pI] == 0 )
            continue;

        const label nLayers =
            pI < label(nLayersForPatch_.size()) ? nLayersForPatch_[pI] : -1;

        Info << "  " << boundaries[pI].patchName() << ":" << nl
             << "    nLayers=" << nLayers
             << " points total=" << nPtTotal[pI]
             << " multiPatchPts=" << nMultiPatchPts[pI]
             << " scale0=" << nScale0[pI]
             << " scale<=0.02=" << nScale002[pI]
             << " scale<=0.05=" << nScale005[pI]
             << " scale<1=" << nScaleLt1[pI]
             << " gapPts=" << nGapSup[pI]
             << " triplePts=" << nTripleSup[pI] << nl
             << "    faces total=" << nFaceTotal[pI]
             << " faceSup=" << nFaceSup[pI] << nl;
    }
}

void boundaryLayers::activate2DMode()
{
    polyMeshGen2DEngine mesh2DEngine(mesh_);
    const boolList& zMinPoint = mesh2DEngine.zMinPoints();
    const boolList& zMaxPoint = mesh2DEngine.zMaxPoints();

    const faceList::subList& bFaces = surfaceEngine().boundaryFaces();
    const labelList& facePatch = surfaceEngine().boundaryFacePatches();

    boolList allZMax(mesh_.boundaries().size(), true);
    boolList allZMin(mesh_.boundaries().size(), true);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];

        forAll(bf, pI)
        {
            if( !zMinPoint[bf[pI]] )
                allZMin[facePatch[bfI]] = false;
            if( !zMaxPoint[bf[pI]] )
                allZMax[facePatch[bfI]] = false;
        }
    }

    //- mark empty patches as already used
    forAll(allZMin, patchI)
    {
        if( allZMin[patchI] ^ allZMax[patchI] )
        {
            treatedPatch_[patchI] = true;
        }
    }

    forAll(treatPatchesWithPatch_, patchI)
    {
        DynList<label>& patches = treatPatchesWithPatch_[patchI];

        for(label i=patches.size()-1;i>=0;--i)
            if( treatedPatch_[patches[i]] )
                patches.removeElement(i);
    }

    is2DMesh_ = true;
}

void boundaryLayers::addLayerForAllPatches()
{
    if( !geometryAnalysed_ )
        findPatchesToBeTreatedTogether();

    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

    if( !patchWiseLayers_ )
    {
        forAll(boundaries, patchI)
        {
            if( patchI < nLayersForPatch_.size()
             && nLayersForPatch_[patchI] == 0 )
                continue;
            addLayerForPatch(patchI);
        }
    }
    else
    {
        newLabelForVertex_.setSize(nPoints_);
        newLabelForVertex_ = -1;
        otherVrts_.clear();
        patchKey_.clear();

        //- avoid generating bnd layer at empty patches in case of 2D meshing
        //- also skip patches with nLayers==0 (inlet, outlet, periodic, etc)
        forAll(treatedPatch_, patchI)
            if( patchI < nLayersForPatch_.size()
             && nLayersForPatch_[patchI] == 0 )
                treatedPatch_[patchI] = true;
        label counter(0);
        forAll(treatedPatch_, patchI)
            if( !treatedPatch_[patchI] )
                ++counter;

        labelList treatedPatches(counter);
        counter = 0;
        forAll(treatedPatch_, i)
            if( !treatedPatch_[i] )
                treatedPatches[counter++] = i;

        //- per-patch BL planning audit before vertex creation
        reportBLPlanningPerPatch();

        //- suppress faces whose multi-patch singular points would produce
        //- zero-thickness extrusion before constructing the BL vertex graph
        suppressFailedSingularityExtrusions(treatedPatches);

        //- audit again after singularity suppression
        reportBLPlanningPerPatch();

        //- create bnd layer vertices
        createNewVertices(treatedPatches);

        //- create bnd layer cells
        createLayerCells(treatedPatches);

        //- OF12 improvement: quality-based BL rollback disabled (O(n^2) hang)
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
