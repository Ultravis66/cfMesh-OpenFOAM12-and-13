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

#include "demandDrivenData.H"
#include "meshSurfaceOptimizer.H"
#include "meshSurfaceEngineModifier.H"
#include "meshSurfaceCheckInvertedVertices.H"
#include "meshOctree.H"
#include "triangle.H"
#include "helperFunctionsPar.H"
#include "meshSurfaceMapper.H"
#include "meshSurfaceMapper2D.H"
#include "polyMeshGen2DEngine.H"
#include "polyMeshGenAddressing.H"
#include "polyMeshGenChecks.H"
#include "labelledPoint.H"
#include "FIFOStack.H"

#include <map>
#include <stdexcept>

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label meshSurfaceOptimizer::findInvertedVertices
(
    boolList& smoothVertex,
    const label nAdditionalLayers
) const
{
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const VRWGraph& pPoints = surfaceEngine_.pointPoints();

    if( smoothVertex.size() != bPoints.size() )
    {
        smoothVertex.setSize(bPoints.size());
        smoothVertex = true;
    }

    label nInvertedTria(0);

    //- check the vertices at the surface
    //- mark the ones where the mesh is tangled
    meshSurfaceCheckInvertedVertices vrtCheck(*partitionerPtr_, smoothVertex);
    const labelHashSet& inverted = vrtCheck.invertedVertices();

    smoothVertex = false;
    forAll(bPoints, bpI)
    {
        if( inverted.found(bPoints[bpI]) )
        {
            ++nInvertedTria;
            smoothVertex[bpI] = true;
        }
    }

    if( Pstream::parRun() )
        reduce(nInvertedTria, sumOp<label>());
    Info << "Number of inverted boundary faces is " << nInvertedTria << endl;

    if( nInvertedTria == 0 )
        return 0;

    //- add additional layers around inverted points
    for(label i=0;i<nAdditionalLayers;++i)
    {
        boolList originallySelected = smoothVertex;
        forAll(smoothVertex, bpI)
            if( originallySelected[bpI] )
                forAllRow(pPoints, bpI, ppI)
                    smoothVertex[pPoints(bpI, ppI)] = true;

        if( Pstream::parRun() )
        {
            //- exchange global labels of inverted points
            const labelList& globalPointLabel =
                surfaceEngine_.globalBoundaryPointLabel();
            const Map<label>& globalToLocal =
                surfaceEngine_.globalToLocalBndPointAddressing();
            const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
            const DynList<label>& neiProcs = surfaceEngine_.bpNeiProcs();

            std::map<label, labelLongList> shareData;
            forAll(neiProcs, procI)
                shareData.insert
                (
                    std::make_pair(neiProcs[procI], labelLongList())
                );

            forAllConstIter(Map<label>, globalToLocal, iter)
            {
                const label bpI = iter();

                if( !smoothVertex[bpI] )
                    continue;

                forAllRow(bpAtProcs, bpI, procI)
                {
                    const label neiProc = bpAtProcs(bpI, procI);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    shareData[neiProc].append(globalPointLabel[bpI]);
                }
            }

            //- exchange data with other processors
            labelLongList receivedData;
            help::exchangeMap(shareData, receivedData);

            forAll(receivedData, j)
            {
                if( !globalToLocal.found(receivedData[j]) ) continue;
                const label bpI = globalToLocal[receivedData[j]];

                smoothVertex[bpI] = true;
            }
        }
    }

    return nInvertedTria;
}

void meshSurfaceOptimizer::smoothEdgePoints
(
    const labelLongList& edgePoints,
    const labelLongList& procEdgePoints
)
{
    List<LongList<labelledPoint> > newPositions(1);
    # ifdef USE_OMP
    newPositions.setSize(omp_get_max_threads());
    # endif

    //- smooth edge vertices
    # ifdef USE_OMP
    # pragma omp parallel num_threads(newPositions.size())
    # endif
    {
        # ifdef USE_OMP
        LongList<labelledPoint>& newPos =
            newPositions[omp_get_thread_num()];
        # else
        LongList<labelledPoint>& newPos = newPositions[0];
        # endif

        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 40)
        # endif
        forAll(edgePoints, i)
        {
            const label bpI = edgePoints[i];

            if( vertexType_[bpI] & (PROCBND | LOCKED) )
                continue;

            newPos.append(labelledPoint(bpI, newEdgePositionLaplacian(bpI)));
        }
    }

    if( Pstream::parRun() )
        edgeNodeDisplacementParallel(procEdgePoints);

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);

    if( !Pstream::parRun() )
    {
        label nCandidates = 0;
        label nAccepted = 0;
        label nRejectedVolumeGate = 0;
        label nTouchingExistingBad = 0;

        scalar minAcceptedPositiveRatio = GREAT;

        // Deterministic serial transaction for feature-edge smoothing.
        // The edge proposal itself is unchanged here. This pass only
        // prevents a currently-positive incident volume cell from being
        // made non-positive by the proposed edge-point motion.
        forAll(newPositions, threadI)
        {
            const LongList<labelledPoint>& newPos =
                newPositions[threadI];

            forAll(newPos, i)
            {
                ++nCandidates;

                const label bpI =
                    newPos[i].pointLabel();

                const point& candidate =
                    newPos[i].coordinates();

                bool touchesExistingBad = false;
                scalar minPositiveRatio = GREAT;

                const bool preservesVolume =
                    surfaceModifier.candidatePreservesPositiveCellVolumes
                    (
                        bpI,
                        candidate,
                        touchesExistingBad,
                        minPositiveRatio
                    );

                if( touchesExistingBad )
                {
                    ++nTouchingExistingBad;
                }

                if( !preservesVolume )
                {
                    ++nRejectedVolumeGate;
                    continue;
                }

                surfaceModifier.moveBoundaryVertexNoUpdate
                (
                    bpI,
                    candidate
                );

                ++nAccepted;

                if( minPositiveRatio < GREAT )
                {
                    minAcceptedPositiveRatio =
                        Foam::min
                        (
                            minAcceptedPositiveRatio,
                            minPositiveRatio
                        );
                }
            }
        }

        Info
            << "[EDGE_SMOOTH_VOLUME_TRANSACTION]"
            << " candidates=" << nCandidates
            << " accepted=" << nAccepted
            << " rejectedVolumeGate=" << nRejectedVolumeGate
            << " touchingExistingBad=" << nTouchingExistingBad
            << " minAcceptedPositiveRatio="
            <<
            (
                minAcceptedPositiveRatio < GREAT
              ? minAcceptedPositiveRatio
              : scalar(-1)
            )
            << endl;
    }
    else
    {
        // Preserve legacy MPI behaviour until the parallel transaction
        // and shared-point synchronization path is explicitly audited.
        forAll(newPositions, threadI)
        {
            const LongList<labelledPoint>& newPos =
                newPositions[threadI];

            forAll(newPos, i)
            {
                surfaceModifier.moveBoundaryVertexNoUpdate
                (
                    newPos[i].pointLabel(),
                    newPos[i].coordinates()
                );
            }
        }
    }

    surfaceModifier.updateGeometry(edgePoints);
}

void meshSurfaceOptimizer::smoothLaplacianFC
(
    const labelLongList& selectedPoints,
    const labelLongList& selectedProcPoints,
    const bool transform
)
{
    List<LongList<labelledPoint> > newPositions(1);
    # ifdef USE_OMP
    newPositions.setSize(omp_get_max_threads());
    # endif

    # ifdef USE_OMP
    # pragma omp parallel num_threads(newPositions.size())
    # endif
    {
        # ifdef USE_OMP
        LongList<labelledPoint>& newPos =
            newPositions[omp_get_thread_num()];
        # else
        LongList<labelledPoint>& newPos = newPositions[0];
        # endif

        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 40)
        # endif
        forAll(selectedPoints, i)
        {
            const label bpI = selectedPoints[i];

            if( vertexType_[bpI] & (PROCBND | LOCKED) )
                continue;

            newPos.append
            (
                labelledPoint(bpI, newPositionLaplacianFC(bpI, transform))
            );
        }
    }

    if( Pstream::parRun() )
        nodeDisplacementLaplacianFCParallel(selectedProcPoints, transform);

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);

    if( !Pstream::parRun() )
    {
        label nCandidates = 0;
        label nAccepted = 0;
        label nRejectedVolumeGate = 0;
        label nTouchingExistingBad = 0;

        scalar minAcceptedPositiveRatio = GREAT;

        // Deterministic serial transaction.
        // Candidate positions were calculated above in parallel, but
        // each commit is checked against the current live volume mesh,
        // including all earlier accepted moves in this pass.
        forAll(newPositions, threadI)
        {
            const LongList<labelledPoint>& newPos =
                newPositions[threadI];

            forAll(newPos, i)
            {
                ++nCandidates;

                const label bpI =
                    newPos[i].pointLabel();

                const point& candidate =
                    newPos[i].coordinates();

                bool touchesExistingBad = false;
                scalar minPositiveRatio = GREAT;

                const bool preservesVolume =
                    surfaceModifier.candidatePreservesPositiveCellVolumes
                    (
                        bpI,
                        candidate,
                        touchesExistingBad,
                        minPositiveRatio
                    );

                if( touchesExistingBad )
                {
                    ++nTouchingExistingBad;
                }

                if( !preservesVolume )
                {
                    ++nRejectedVolumeGate;
                    continue;
                }

                surfaceModifier.moveBoundaryVertexNoUpdate
                (
                    bpI,
                    candidate
                );

                ++nAccepted;

                if( minPositiveRatio < GREAT )
                {
                    minAcceptedPositiveRatio =
                        Foam::min
                        (
                            minAcceptedPositiveRatio,
                            minPositiveRatio
                        );
                }
            }
        }

        Info
            << "[LAPLACIAN_FC_VOLUME_TRANSACTION]"
            << " transform=" << transform
            << " candidates=" << nCandidates
            << " accepted=" << nAccepted
            << " rejectedVolumeGate=" << nRejectedVolumeGate
            << " touchingExistingBad=" << nTouchingExistingBad
            << " minAcceptedPositiveRatio="
            <<
            (
                minAcceptedPositiveRatio < GREAT
              ? minAcceptedPositiveRatio
              : scalar(-1)
            )
            << endl;
    }
    else
    {
        // Preserve legacy MPI behaviour until parallel transactional
        // synchronization is explicitly audited.
        forAll(newPositions, threadI)
        {
            const LongList<labelledPoint>& newPos =
                newPositions[threadI];

            forAll(newPos, i)
            {
                surfaceModifier.moveBoundaryVertexNoUpdate
                (
                    newPos[i].pointLabel(),
                    newPos[i].coordinates()
                );
            }
        }
    }

    surfaceModifier.updateGeometry(selectedPoints);
}

void meshSurfaceOptimizer::smoothSurfaceOptimizer
(
    const labelLongList& selectedPoints,
    const labelLongList& selectedProcPoints
)
{
    //- create partTriMesh is it is not yet present
    this->triMesh();

    //- update coordinates of the triangulation
    updateTriMesh(selectedPoints);

    pointField newPositions(selectedPoints.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 20)
    # endif
    forAll(selectedPoints, i)
    {
        const label bpI = selectedPoints[i];

        newPositions[i] = newPositionSurfaceOptimizer(bpI);
    }

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);

    if( !Pstream::parRun() )
    {
        label nCandidates = 0;
        label nAccepted = 0;
        label nRejectedVolumeGate = 0;
        label nTouchingExistingBad = 0;

        scalar minAcceptedPositiveRatio = GREAT;

        // Deterministic serial transaction for surface-optimizer proposals.
        // The optimization target itself is unchanged. Each candidate is
        // evaluated against the current live volume state, including all
        // earlier accepted moves in this pass.
        forAll(newPositions, i)
        {
            ++nCandidates;

            const label bpI = selectedPoints[i];
            const point& candidate = newPositions[i];

            bool touchesExistingBad = false;
            scalar minPositiveRatio = GREAT;

            const bool preservesVolume =
                surfaceModifier.candidatePreservesPositiveCellVolumes
                (
                    bpI,
                    candidate,
                    touchesExistingBad,
                    minPositiveRatio
                );

            if( touchesExistingBad )
            {
                ++nTouchingExistingBad;
            }

            if( !preservesVolume )
            {
                ++nRejectedVolumeGate;
                continue;
            }

            surfaceModifier.moveBoundaryVertexNoUpdate
            (
                bpI,
                candidate
            );

            ++nAccepted;

            if( minPositiveRatio < GREAT )
            {
                minAcceptedPositiveRatio =
                    Foam::min
                    (
                        minAcceptedPositiveRatio,
                        minPositiveRatio
                    );
            }
        }

        Info
            << "[SURFACE_OPT_VOLUME_TRANSACTION]"
            << " candidates=" << nCandidates
            << " accepted=" << nAccepted
            << " rejectedVolumeGate=" << nRejectedVolumeGate
            << " touchingExistingBad=" << nTouchingExistingBad
            << " minAcceptedPositiveRatio="
            <<
            (
                minAcceptedPositiveRatio < GREAT
              ? minAcceptedPositiveRatio
              : scalar(-1)
            )
            << endl;
    }
    else
    {
        // Preserve legacy parallel behaviour until shared-point
        // transactional synchronization is explicitly audited.
        forAll(newPositions, i)
        {
            const label bpI = selectedPoints[i];

            surfaceModifier.moveBoundaryVertexNoUpdate
            (
                bpI,
                newPositions[i]
            );
        }
    }

    //- ensure that vertices at inter-processor boundaries are at the same
    //- location at all processors
    surfaceModifier.syncVerticesAtParallelBoundaries(selectedProcPoints);

    //- update geometry addressing for moved points
    surfaceModifier.updateGeometry(selectedPoints);
}

bool meshSurfaceOptimizer::untangleSurface
(
    const labelLongList& selectedBoundaryPoints,
    const label nAdditionalLayers
)
{
    Info << "Starting untangling the surface of the volume mesh" << endl;

    bool changed(false);

    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const pointFieldPMG& points = surfaceEngine_.points();
    surfaceEngine_.pointFaces();
    surfaceEngine_.faceCentres();
    surfaceEngine_.pointPoints();
    surfaceEngine_.boundaryFacePatches();
    surfaceEngine_.pointNormals();
    surfaceEngine_.boundaryPointEdges();

    if( Pstream::parRun() )
    {
        surfaceEngine_.bpAtProcs();
        surfaceEngine_.globalToLocalBndPointAddressing();
        surfaceEngine_.globalBoundaryPointLabel();
        surfaceEngine_.bpNeiProcs();
    }

    boolList smoothVertex(bPoints.size(), false);
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(selectedBoundaryPoints, i)
    {
        if( vertexType_[selectedBoundaryPoints[i]] & LOCKED )
            continue;

        smoothVertex[selectedBoundaryPoints[i]] = true;
    }

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);

    meshSurfaceMapper* mapperPtr = NULL;
    if( octreePtr_ )
        mapperPtr = new meshSurfaceMapper(*partitionerPtr_, *octreePtr_);

    bool remapVertex(true);
    label nInvertedTria;
    label nGlobalIter(0);

    // Distinguish the many untangleSurface() calls made during one mesh run.
    static label untangleVolumeCallCounter = 0;
    const label untangleVolumeCallId =
        ++untangleVolumeCallCounter;

    // Diagnostic-only raw volume lineage for untangleSurface.
    // Cache coherency is now maintained by meshSurfaceEngineModifier.
    auto untangleVolumeLineage =
    [&](const word& stageName, const label outerI, const label innerI)
    {
        labelHashSet negVolCells;

        polyMeshGenChecks::checkCellVolumes
        (
            surfaceEngine_.mesh(),
            false,
            &negVolCells
        );

        Info << "[UNTANGLE_VOLUME_LINEAGE]"
             << " call=" << untangleVolumeCallId
             << " outer=" << outerI
             << " inner=" << innerI
             << " stage=" << stageName
             << " negVol=" << negVolCells.size()
             << endl;
    };

    untangleVolumeLineage("entry", -1, -1);

    labelLongList procBndPoints, movedPoints;
    labelLongList procEdgePoints, movedEdgePoints;

    label minNumInverted(bPoints.size());
    FIFOStack<label> nInvertedHistory;
    pointField minInvertedPoints(bPoints.size());

    do
    {
        label nIter(0), nAfterRefresh(0);

        do
        {
            nInvertedTria =
                findInvertedVertices(smoothVertex, nAdditionalLayers);

            untangleVolumeLineage
            (
                "iterEntry",
                nGlobalIter,
                nIter
            );

            if( nInvertedTria == 0 )
            {
                break;
            }
            else if( enforceConstraints_ && !remapVertex )
            {
                polyMeshGen& mesh =
                    const_cast<polyMeshGen&>(surfaceEngine_.mesh());

                const label subsetId =
                    mesh.addPointSubset(badPointsSubsetName_);

                forAll(smoothVertex, bpI)
                    if( smoothVertex[bpI] )
                        mesh.addPointToSubset(subsetId, bPoints[bpI]);

                WarningIn
                (
                    "bool meshSurfaceOptimizer::untangleSurface"
                    "(const labelLongList&, const label)"
                ) << "Writing mesh with " << badPointsSubsetName_
                  << " subset. These points cannot be untangled"
                  << " without sacrificing geometry constraints. Exitting.."
                  << endl;

                returnReduce(1, sumOp<label>());

                throw std::logic_error
                (
                    "bool meshSurfaceOptimizer::untangleSurface"
                    "(const labelLongList&, const label)"
                    "Cannot untangle mesh!!"
                );
            }

            //- find the min number of inverted points and
            //- add the last number to the stack
            if( nInvertedTria < minNumInverted )
            {
                minNumInverted = nInvertedTria;
                nAfterRefresh = 0;

                # ifdef USE_OMP
                # pragma omp parallel for schedule(dynamic, 100)
                # endif
                forAll(bPoints, bpI)
                    minInvertedPoints[bpI] = points[bPoints[bpI]];
            }

            //- count the number of iteration after the last minimum occurence
            ++nAfterRefresh;

            //- update the stack
            nInvertedHistory.push(nInvertedTria);
            if( nInvertedHistory.size() > 2 )
                nInvertedHistory.pop();

            //- check if the number of inverted points reduces
            bool minimumInStack(false);
            forAllConstIter(FIFOStack<label>, nInvertedHistory, it)
                if( it() == minNumInverted )
                    minimumInStack = true;

            //- stop if the procedure does not minimise
            //- the number of inverted points
            if( !minimumInStack || (nAfterRefresh > 2) )
                break;

            //- find points which will be handled by the smoothers
            changed = true;

            procBndPoints.clear();
            movedPoints.clear();
            procEdgePoints.clear();
            movedEdgePoints.clear();

            forAll(bPoints, bpI)
            {
                if( !smoothVertex[bpI] )
                    continue;

                if( vertexType_[bpI] & PARTITION )
                {
                    movedPoints.append(bpI);

                    if( vertexType_[bpI] & PROCBND )
                        procBndPoints.append(bpI);
                }
                else if( vertexType_[bpI] & EDGE )
                {
                    movedEdgePoints.append(bpI);

                    if( vertexType_[bpI] & PROCBND )
                        procEdgePoints.append(bpI);
                }
            }

            //- Feature-edge motion: resolve smoothing intent onto the
            //- feature BEFORE modifying the live mesh.
            //
            //- The old path committed an off-feature Laplacian midpoint,
            //- then remapped it, then blindly restored a stale pre-smoothing
            //- position when remapping failed. Neighbouring points may have
            //- moved meanwhile, so that rollback was not a valid transaction.
            //
            //- In the serial transactional path below, the midpoint is only
            //- a geometric intent. mapEdgeNodes() resolves the final exact
            //- native/virtual feature target while the point remains at its
            //- current live position, validates the final move, and performs
            //- at most one commit.
            const bool directFeatureTransaction =
                remapVertex
             && mapperPtr
             && transactionalFeatureOptimization_
             && !Pstream::parRun();

            if( directFeatureTransaction )
            {
                pointField desiredEdgePositions
                (
                    movedEdgePoints.size(),
                    point::zero
                );

                forAll(movedEdgePoints, epI)
                {
                    desiredEdgePositions[epI] =
                        newEdgePositionLaplacian
                        (
                            movedEdgePoints[epI]
                        );
                }

                // No mesh mutation has occurred yet. This checkpoint must
                // therefore match iterEntry apart from unrelated prior state.
                untangleVolumeLineage
                (
                    "afterEdgeIntentBeforeFeatureTransaction",
                    nGlobalIter,
                    nIter
                );

                boolList mappingAccepted;

                mapperPtr->mapEdgeNodes
                (
                    movedEdgePoints,
                    desiredEdgePositions,
                    mappingAccepted
                );

                label nAcceptedFeature = 0;

                forAll(mappingAccepted, epI)
                {
                    if( mappingAccepted[epI] )
                    {
                        ++nAcceptedFeature;
                    }
                }

                Info
                    << "[EDGE_1DOF_TRANSACTION]"
                    << " phase=untangle"
                    << " input=" << movedEdgePoints.size()
                    << " accepted=" << nAcceptedFeature
                    << " rejected="
                    << (movedEdgePoints.size() - nAcceptedFeature)
                    << endl;
            }
            else
            {
                // Legacy path retained for MPI, disabled transactions, and
                // remapVertex=false. The serial remap=true transactional path
                // above never creates an off-feature intermediate state.
                pointField oldEdgePositions(movedEdgePoints.size());

                forAll(movedEdgePoints, epI)
                {
                    oldEdgePositions[epI] =
                        points[bPoints[movedEdgePoints[epI]]];
                }

                smoothEdgePoints
                (
                    movedEdgePoints,
                    procEdgePoints
                );

                untangleVolumeLineage
                (
                    "afterEdgeSmoothingBeforeRemap",
                    nGlobalIter,
                    nIter
                );

                if( remapVertex && mapperPtr )
                {
                    if
                    (
                        transactionalFeatureOptimization_
                     && !Pstream::parRun()
                    )
                    {
                        boolList mappingAccepted;

                        mapperPtr->mapEdgeNodes
                        (
                            movedEdgePoints,
                            mappingAccepted
                        );

                        label nRolledBack = 0;

                        forAll(movedEdgePoints, epI)
                        {
                            if( mappingAccepted[epI] )
                                continue;

                            surfaceModifier.moveBoundaryVertexNoUpdate
                            (
                                movedEdgePoints[epI],
                                oldEdgePositions[epI]
                            );

                            ++nRolledBack;
                        }

                        Info
                            << "[EDGEOPT_TRANSACTION]"
                            << " phase=untangle"
                            << " input=" << movedEdgePoints.size()
                            << " accepted="
                            << (movedEdgePoints.size() - nRolledBack)
                            << " rolledBack=" << nRolledBack
                            << endl;
                    }
                    else
                    {
                        mapperPtr->mapEdgeNodes
                        (
                            movedEdgePoints
                        );
                    }
                }
            }

            surfaceModifier.updateGeometry(movedEdgePoints);

            untangleVolumeLineage
            (
                "afterEdgeTransaction",
                nGlobalIter,
                nIter
            );

            //- use laplacian smoothing
            smoothLaplacianFC(movedPoints, procBndPoints);
            surfaceModifier.updateGeometry(movedPoints);

            untangleVolumeLineage
            (
                "afterPartitionLaplacian",
                nGlobalIter,
                nIter
            );

            //- use surface optimizer
            smoothSurfaceOptimizer(movedPoints, procBndPoints);

            untangleVolumeLineage
            (
                "afterPartitionSurfaceOptimizer",
                nGlobalIter,
                nIter
            );

            if( remapVertex && mapperPtr )
                mapperPtr->mapVerticesOntoSurface(movedPoints);

            //- update normals and other geometric data
            surfaceModifier.updateGeometry(movedPoints);

            untangleVolumeLineage
            (
                "afterPartitionRemap",
                nGlobalIter,
                nIter
            );

        } while( nInvertedTria && (++nIter < 20) );

        if( nInvertedTria > 0 )
        {
            //- use the combination with the minimum number of inverted points
            meshSurfaceEngineModifier sMod(surfaceEngine_);
            forAll(minInvertedPoints, bpI)
                sMod.moveBoundaryVertexNoUpdate(bpI, minInvertedPoints[bpI]);

            sMod.updateGeometry();

            untangleVolumeLineage
            (
                "afterMinInvertedRestore",
                nGlobalIter,
                nIter
            );
        }

        if( nInvertedTria )
        {
            Info << "Smoothing remaining inverted vertices " << endl;

            movedPoints.clear();
            procBndPoints.clear();
            forAll(smoothVertex, bpI)
                if( smoothVertex[bpI] )
                {
                    // Fallback: only smooth PARTITION points here.
                    // EDGE points must not be Laplacian-smoothed as surface
                    // points -- that moves them off feature curves.
                    // mapVerticesOntoSurface is also unsafe for edge points
                    // (global projection); use patch-constrained version.
                    if( !(vertexType_[bpI] & PARTITION) ) continue;

                    movedPoints.append(bpI);

                    if( vertexType_[bpI] & PROCBND )
                        procBndPoints.append(bpI);
                }

            smoothLaplacianFC(movedPoints, procBndPoints, false);

            untangleVolumeLineage
            (
                "afterFallbackPartitionLaplacian",
                nGlobalIter,
                nIter
            );

            // Use patch-constrained projection -- mapVerticesOntoSurface
            // uses global nearest-surface which is unsafe near junctions.
            if( remapVertex && mapperPtr )
                mapperPtr->mapVerticesOntoSurfacePatches(movedPoints);

            //- update normals and other geometric data
            surfaceModifier.updateGeometry(movedPoints);

            untangleVolumeLineage
            (
                "afterFallbackPartitionRemap",
                nGlobalIter,
                nIter
            );

            if( nGlobalIter > 5 )
                remapVertex = false;
        }

    } while( nInvertedTria && (++nGlobalIter < 10) );

    untangleVolumeLineage("exit", nGlobalIter, -1);

    deleteDemandDrivenData(mapperPtr);

    if( nInvertedTria != 0 )
    {
        //- the procedure has given up without success
        //- there exist some remaining inverted faces in the mesh
        polyMeshGen& mesh =
            const_cast<polyMeshGen&>(surfaceEngine_.mesh());

        label subsetId = mesh.pointSubsetIndex(badPointsSubsetName_);
        if( subsetId >= 0 )
            mesh.removePointSubset(subsetId);
        subsetId = mesh.addPointSubset(badPointsSubsetName_);

        forAll(smoothVertex, bpI)
            if( smoothVertex[bpI] )
                mesh.addPointToSubset(subsetId, bPoints[bpI]);
    }

    Info << "Finished untangling the surface of the volume mesh" << endl;

    return changed;
}

bool meshSurfaceOptimizer::untangleSurface(const label nAdditionalLayers)
{
    labelLongList selectedPts(surfaceEngine_.boundaryPoints().size());
    forAll(selectedPts, i)
        selectedPts[i] = i;

    return untangleSurface(selectedPts, nAdditionalLayers);
}

void meshSurfaceOptimizer::optimizeSurface(const label nIterations)
{
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const pointFieldPMG& points = surfaceEngine_.points();

    //- needed for parallel execution
    surfaceEngine_.pointFaces();
    surfaceEngine_.faceCentres();
    surfaceEngine_.pointPoints();
    surfaceEngine_.boundaryPointEdges();
    surfaceEngine_.boundaryFacePatches();
    surfaceEngine_.pointNormals();
    surfaceEngine_.boundaryPointEdges();

    meshSurfaceMapper* mapperPtr = NULL;
    if( octreePtr_ )
        mapperPtr = new meshSurfaceMapper(*partitionerPtr_, *octreePtr_);

    labelLongList procBndPoints, edgePoints, partitionPoints, procPoints;
    forAll(bPoints, bpI)
    {
        if( vertexType_[bpI] & LOCKED )
            continue;

        if( vertexType_[bpI] & EDGE )
        {
            edgePoints.append(bpI);

            if( vertexType_[bpI] & PROCBND )
                procBndPoints.append(bpI);
        }
        else if( vertexType_[bpI] & PARTITION )
        {
            partitionPoints.append(bpI);

            if( vertexType_[bpI] & PROCBND )
                procPoints.append(bpI);
        }
    }

    //- optimize edge vertices
    Info << "Optimizing edges. Iteration:" << flush;
    for(label i=0;i<nIterations;++i)
    {
        Info << "." << flush;

        meshSurfaceEngineModifier bMod(surfaceEngine_);

        // Snapshot the last valid feature positions before the unconstrained
        // Laplacian edge proposal.
        pointField oldEdgePositions(edgePoints.size());
        forAll(edgePoints, epI)
        {
            oldEdgePositions[epI] =
                points[bPoints[edgePoints[epI]]];
        }

        smoothEdgePoints(edgePoints, procBndPoints);

        //- Project vertices back onto the true feature boundary.
        //
        // In serial, treat smoothing + constrained projection as one
        // transaction.  A rejected projection must not leave the preceding
        // off-feature Laplacian proposal in the mesh.
        if( mapperPtr )
        {
            if
            (
                transactionalFeatureOptimization_
             && !Pstream::parRun()
            )
            {
                boolList mappingAccepted;
                mapperPtr->mapEdgeNodes(edgePoints, mappingAccepted);

                label nRolledBack = 0;

                forAll(edgePoints, epI)
                {
                    if( mappingAccepted[epI] )
                        continue;

                    bMod.moveBoundaryVertexNoUpdate
                    (
                        edgePoints[epI],
                        oldEdgePositions[epI]
                    );

                    ++nRolledBack;
                }

                Info
                    << "[EDGEOPT_TRANSACTION]"
                    << " phase=optimize"
                    << " iteration=" << i
                    << " input=" << edgePoints.size()
                    << " accepted="
                    << (edgePoints.size() - nRolledBack)
                    << " rolledBack=" << nRolledBack
                    << endl;
            }
            else
            {
                // Keep existing parallel behaviour pending an explicit audit
                // of mapToSmallestDistance() cross-rank acceptance semantics.
                mapperPtr->mapEdgeNodes(edgePoints);
            }
        }

        //- update the geometry information
        bMod.updateGeometry(edgePoints);
    }
    Info << endl;

    //- delete the mapper
    deleteDemandDrivenData(mapperPtr);

    //- optimize positions of surface vertices which are not on surface edges
    Info << "Optimizing surface vertices. Iteration:";
    for(label i=0;i<nIterations;++i)
    {
        smoothLaplacianFC(partitionPoints, procPoints, true);

        smoothSurfaceOptimizer(partitionPoints, procPoints);

        Info << "." << flush;
    }

    Info << endl;

    untangleSurface(0);
}

void meshSurfaceOptimizer::optimizeSurface2D(const label nIterations)
{
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const edgeList& edges = surfaceEngine_.edges();
    const labelList& bp = surfaceEngine_.bp();

    polyMeshGen2DEngine mesh2DEngine
    (
        const_cast<polyMeshGen&>(surfaceEngine_.mesh())
    );
    const boolList& zMinPoint = mesh2DEngine.zMinPoints();

    //- needed for parallel execution
    surfaceEngine_.pointFaces();
    surfaceEngine_.faceCentres();
    surfaceEngine_.pointPoints();
    surfaceEngine_.boundaryPointEdges();
    surfaceEngine_.boundaryFacePatches();
    surfaceEngine_.pointNormals();

    labelLongList procBndPoints, movedPoints, activeEdges, updatePoints;
    forAll(edges, beI)
    {
        const edge& e = edges[beI];

        if( zMinPoint[e.start()] ^ zMinPoint[e.end()] )
        {
            label bpI = bp[e.start()];
            if( !zMinPoint[e.start()] )
                bpI = bp[e.end()];

            if( vertexType_[bpI] & EDGE )
            {
                activeEdges.append(beI);

                updatePoints.append(bp[e.start()]);
                updatePoints.append(bp[e.end()]);

                movedPoints.append(bpI);

                if( vertexType_[bpI] & PROCBND )
                    procBndPoints.append(bpI);
            }
        }
    }

    meshSurfaceMapper2D* mapperPtr = NULL;
    if( octreePtr_ )
        mapperPtr = new meshSurfaceMapper2D(surfaceEngine_, *octreePtr_);

    //- optimize edge vertices
    meshSurfaceEngineModifier bMod(surfaceEngine_);

    Info << "Optimizing edges. Iteration:" << flush;
    for(label i=0;i<nIterations;++i)
    {
        Info << "." << flush;

        smoothEdgePoints(movedPoints, procBndPoints);

        //- move points with maximum z coordinate
        mesh2DEngine.correctPoints();

        //- map boundary edges to the surface
        if( mapperPtr )
            mapperPtr->mapVerticesOntoSurfacePatches(activeEdges);

        //- update normal, centres, etc, after the surface has been modified
        bMod.updateGeometry(updatePoints);
    }
    Info << endl;

    //- optimize Pts of surface vertices which are not on surface edges
    procBndPoints.clear();
    movedPoints.clear();
    forAll(bPoints, bpI)
        if( zMinPoint[bPoints[bpI]] && (vertexType_[bpI] & PARTITION) )
        {
            movedPoints.append(bpI);

            if( vertexType_[bpI] & PROCBND )
                procBndPoints.append(bpI);
        }
    Info << "Optimizing surface vertices. Iteration:";
    for(label i=0;i<nIterations;++i)
    {
        Info << "." << flush;

        smoothSurfaceOptimizer(movedPoints, procBndPoints);

        //- move the points which are not at minimum z coordinate
        mesh2DEngine.correctPoints();

        //- update geometrical data due to movement of vertices
        bMod.updateGeometry(movedPoints);
    }

    Info << endl;

    deleteDemandDrivenData(mapperPtr);
}

void meshSurfaceOptimizer::untangleSurface2D()
{
    const polyMeshGen& mesh = surfaceEngine_.mesh();
    const faceListPMG& faces = mesh.faces();
    const VRWGraph& pointFaces = mesh.addressingData().pointFaces();

    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const labelList& bp = surfaceEngine_.bp();

    polyMeshGen2DEngine mesh2DEngine(const_cast<polyMeshGen&>(mesh));
    const boolList& zMinPoint = mesh2DEngine.zMinPoints();
    const boolList& activeFace = mesh2DEngine.activeFace();

    //- needed for parallel execution
    surfaceEngine_.pointFaces();
    surfaceEngine_.faceCentres();
    surfaceEngine_.pointPoints();
    surfaceEngine_.boundaryPointEdges();
    surfaceEngine_.boundaryFacePatches();
    surfaceEngine_.pointNormals();

    boolList activeBoundaryPoint(bPoints.size());
    boolList changedFace(activeFace.size(), true);

    label iterationI(0);
    do
    {
        labelHashSet badFaces;
        const label nBadFaces =
            polyMeshGenChecks::findBadFaces
            (
                mesh,
                badFaces,
                false,
                &changedFace
            );

        Info << "Iteration " << iterationI
             << ". Number of bad faces " << nBadFaces << endl;

        if( nBadFaces == 0 )
            break;

        //- update active points and faces affected by the movement
        //- of active points
        activeBoundaryPoint = false;
        changedFace = false;
        forAllConstIter(labelHashSet, badFaces, it)
        {
            const face& f = faces[it.key()];

            forAll(f, pI)
            {
                if( zMinPoint[f[pI]] )
                {
                    activeBoundaryPoint[bp[f[pI]]] = true;

                    forAllRow(pointFaces, f[pI], pfI)
                        changedFace[pointFaces(f[pI], pfI)] = true;
                }
            }
        }

        if( Pstream::parRun() )
        {
            const Map<label>& globalToLocal =
                surfaceEngine_.globalToLocalBndPointAddressing();
            const DynList<label>& neiProcs = surfaceEngine_.bpNeiProcs();
            const VRWGraph& bpNeiProcs = surfaceEngine_.bpAtProcs();

            std::map<label, labelLongList> exchangeData;
            forAll(neiProcs, i)
                exchangeData[neiProcs[i]].clear();

            //- collect active points at inter-processor boundaries
            forAllConstIter(Map<label>, globalToLocal, it)
            {
                const label bpI = it();

                if( activeBoundaryPoint[bpI] )
                {
                    forAllRow(bpNeiProcs, bpI, i)
                    {
                        const label neiProc = bpNeiProcs(bpI, i);

                        if( neiProc == Pstream::myProcNo() )
                            continue;

                        exchangeData[neiProc].append(it.key());
                    }
                }
            }

            //- exchange active points among the processors
            labelLongList receivedData;
            help::exchangeMap(exchangeData, receivedData);

            //- ensure that all processors have the same Pts active
            forAll(receivedData, i)
            {
                if( !globalToLocal.found(receivedData[i]) ) continue;
                const label bpI = globalToLocal[receivedData[i]];

                //- activate this boundary point
                activeBoundaryPoint[bpI] = true;

                //- set the changeFaces for the faces attached to this point
                forAllRow(pointFaces, bPoints[bpI], pfI)
                    changedFace[pointFaces(bPoints[bpI], pfI)] = true;
            }
        }

        //- apply smoothing to the activated points
        meshSurfaceEngineModifier bMod(surfaceEngine_);

        labelLongList movedPts, procBndPts, edgePts, procEdgePts;
        forAll(bPoints, bpI)
        {
            if( !activeBoundaryPoint[bpI] )
                continue;

            if( vertexType_[bpI] & EDGE )
            {
                edgePts.append(bpI);

                if( vertexType_[bpI] & PROCBND )
                    procEdgePts.append(bpI);
            }
            else if( vertexType_[bpI] & PARTITION )
            {
                movedPts.append(bpI);

                if( vertexType_[bpI] & PROCBND )
                    procBndPts.append(bpI);
            }
        }

        for(label i=0;i<5;++i)
        {
            smoothEdgePoints(edgePts, procEdgePts);

            bMod.updateGeometry(edgePts);

            smoothSurfaceOptimizer(movedPts, procBndPts);

            bMod.updateGeometry(movedPts);
        }

        //- move the points which are not at minimum z coordinate
        mesh2DEngine.correctPoints();

        //- update geometrical data due to movement of vertices
        bMod.updateGeometry();

        //- update cell centres and face centres
        const_cast<polyMeshGenAddressing&>
        (
            mesh.addressingData()
        ).updateGeometry(changedFace);

    } while( ++iterationI < 20 );

    //- delete invalid data
    mesh.clearAddressingData();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
