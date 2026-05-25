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

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGLayer

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

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
            const label beI =
                globalToLocal[receivedData[counter++].pointLabel()];

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
    blblCornerAcuteThreshold_(0.3),
    layerScaleRing1_(0.25),
    layerScaleRing2_(0.50),
    layerScaleRing3_(0.75),
    layerScaleRing4_(0.60),
    layerScaleRing5_(0.80),
    layerScaleRing6_(1.00)
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

    // Classify patches from patchRole_ — single source of truth
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

    // Classify patches using patchRole_ — the single source of truth.
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
            // — no suppression triggered at periodic/symmetry junctions
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
                ++nGapSuppressed;
            }
        }
        Info << "Gap detection: suppressed "
             << nGapSuppressed
             << " BL points in thin clearance regions" << endl;
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
                zeroPts[bp0] = true;
            }
            if( bp1 >= 0 )
            {
                zeroDistPoints_[bp1] = true;
                layerScale_[bp1] = 0.02;
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
            layerScale_[bpI] = isAcute ? 0.02 : 0.15;
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
         << " suppressed=" << nTriple << endl;

    // Topology-aware scale assignment:
    // Upgrade transition point suppression based on point topology class.
    // Corner points (3+ patches) at BL/no-BL transitions are geometrically
    // overconstrained — full suppress regardless of angle.
    // Two-patch edge points get full suppress only if patch angle is sharp.
    {
        const scalar cosSharp = Foam::cos(scalar(75.0) * M_PI / 180.0);
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
            // Corner: 3+ patches with BL + explicit termination — full suppress
            // Neutral patches (periodic etc) do not trigger suppression
            if( nPt >= 3 && nBLPt >= 1 && nTermPt >= 1 )
            {
                layerScale_[bpI] = 0.0;
                ++nCornerSuppressed;
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
        const scalar cosThresh = Foam::cos(blblFeatureAngleDeg_ * M_PI / 180.0);
        label nBLBL = 0;
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
                zeroDistPoints_[bpI] = true;
                layerScale_[bpI] = 0.0;
                zeroPts[bpI] = true;
                ++nBLBL;
                // C1: persist junction point for topology fix in C2/C3
                blblJunctionPoints_.insert(bpI);
            }
        }
        Info << "BL/BL sharp-junction suppression: "
             << nBLBL << " points suppressed, "
             << blblJunctionPoints_.size() << " junction points captured" << endl;
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
                    zeroDistPoints_[bpI] = true;
                    layerScale_[bpI] = 0.0;
                    zeroPts[bpI] = true;
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
            ring1[nbpI] = true;
            layerScale_[nbpI] = Foam::min(layerScale_[nbpI], layerScaleRing1_);
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
            ring2[nbpI] = true;
            layerScale_[nbpI] = Foam::min(layerScale_[nbpI], layerScaleRing2_);
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
            ring3[nbpI] = true;
            layerScale_[nbpI] = Foam::min(layerScale_[nbpI], layerScaleRing3_);
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
            ring4[nbpI] = true;
            layerScale_[nbpI] = Foam::min(layerScale_[nbpI], layerScaleRing4_);
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
            ring5[nbpI] = true;
            layerScale_[nbpI] = Foam::min(layerScale_[nbpI], layerScaleRing5_);
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
            ring6[nbpI] = true;
            layerScale_[nbpI] = Foam::min(layerScale_[nbpI], layerScaleRing6_);
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
    Info << "terminateLayersAtConcaveEdges: marked "
         << nTransitionEdges << " BL-transition edges." << endl;
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
