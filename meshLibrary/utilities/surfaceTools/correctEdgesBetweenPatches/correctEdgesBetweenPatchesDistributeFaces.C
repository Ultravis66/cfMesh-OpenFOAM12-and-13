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
#include "correctEdgesBetweenPatches.H"
#include "decomposeFaces.H"
#include "triSurface.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceCheckEdgeTypes.H"
#include "meshSurfacePartitioner.H"
#include "helperFunctions.H"
#include "helperFunctionsPar.H"

#include <map>

# ifdef USE_OMP
#include <omp.h>
#include "polyMeshGenChecks.H"
# endif

//#define DEBUGMapping

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void correctEdgesBetweenPatches::decomposeProblematicFaces()
{
    Info << "Decomposing problematic faces" << endl;
    const meshSurfaceEngine& mse = meshSurface();
    const labelList& bp = mse.bp();
    const edgeList& edges = mse.edges();
    const VRWGraph& pointEdges = mse.boundaryPointEdges();
    const VRWGraph& faceEdges = mse.faceEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();
    const labelList& facePatches = mse.boundaryFacePatches();

    //- mark feature edges
    boolList featureBndEdge(edgeFaces.size(), false);

    forAll(edgeFaces, beI)
    {
        if( edgeFaces.sizeOfRow(beI) != 2 )
            continue;

        if( facePatches[edgeFaces(beI, 0)] != facePatches[edgeFaces(beI, 1)] )
            featureBndEdge[beI] = true;
    }

    if( Pstream::parRun() )
    {
        //- find feature edges at parallel boundaries and propagate the
        //- information to all processors
        const Map<label>& globalToLocalEdge =
            mse.globalToLocalBndEdgeAddressing();
        const VRWGraph& beAtProcs = mse.beAtProcs();
        const Map<label>& otherProcPatch = mse.otherEdgeFacePatch();
        forAllConstIter(Map<label>, globalToLocalEdge, it)
        {
            const label beI = it();

            if( edgeFaces.sizeOfRow(beI) != 1 )
                continue;
            if( facePatches[edgeFaces(beI, 0)] != otherProcPatch[beI] )
                featureBndEdge[beI] = true;
        }

        //- propagate information to all processors that need this information
        std::map<label, labelLongList> exchangeData;
        forAll(mse.beNeiProcs(), i)
            exchangeData.insert
            (
                std::make_pair(mse.beNeiProcs()[i], labelLongList())
            );

        //- append labels of feature edges that need to be sent to other
        //- processors sharing that edge
        forAllConstIter(Map<label>, globalToLocalEdge, it)
        {
            const label beI = it();

            if( featureBndEdge[beI] )
            {
                forAllRow(beAtProcs, beI, i)
                {
                    const label procI = beAtProcs(beI, i);
                    if( procI == Pstream::myProcNo() )
                        continue;

                    exchangeData[procI].append(it.key());
                }
            }
        }

        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        label counter(0);
        while( counter < receivedData.size() )
        {
            const label geI = receivedData[counter++];
            if( !globalToLocalEdge.found(geI) ) continue;
            if( !globalToLocalEdge.found(geI) ) continue;
            featureBndEdge[globalToLocalEdge[geI]] = true;
        }
    }

    const polyMeshGen& mesh = mse.mesh();
    const faceListPMG& faces = mesh.faces();
    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();

    boolList decomposeFace(faces.size(), false);
    label nDecomposedFaces(0);

    //- decompose internal faces with more than one feature edge
    const label nIntFaces = mesh.nInternalFaces();
    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided) reduction(+ : nDecomposedFaces)
    # endif
    for(label faceI=0;faceI<nIntFaces;++faceI)
    {
        const face& f = faces[faceI];

        label nFeatureEdges(0);

        forAll(f, eI)
        {
            const edge e = f.faceEdge(eI);

            const label bs = bp[e[0]];
            const label be = bp[e[1]];
            if( (bs != -1) && (be != -1) )
            {
                //- check if this edge is a boundary edges and a feature edge
                forAllRow(pointEdges, bs, i)
                {
                    const label beI = pointEdges(bs, i);

                    if( (edges[beI] == e) && featureBndEdge[beI] )
                        ++nFeatureEdges;
                }
            }
        }

        if( nFeatureEdges > 1 )
        {
            ++nDecomposedFaces;
            decomposeFace[faceI] = true;
            decomposeCell_[owner[faceI]] = true;
            decomposeCell_[neighbour[faceI]] = true;
        }
    }

    //- decompose boundary faces in case the feature edges are not connected
    //- into a single open chain of edges
    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided) reduction(+ : nDecomposedFaces)
    # endif
    forAll(faceEdges, bfI)
    {
        boolList featureEdge(faceEdges.sizeOfRow(bfI), false);

        forAllRow(faceEdges, bfI, feI)
            if( featureBndEdge[faceEdges(bfI, feI)] )
                featureEdge[feI] = true;

        if( !help::areElementsInChain(featureEdge) )
        {
            ++nDecomposedFaces;
            decomposeFace[nIntFaces+bfI] = true;
            decomposeCell_[owner[nIntFaces+bfI]] = true;
        }
    }

    if( Pstream::parRun() )
    {
        //- decompose processor faces having more than one feature edge
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh.procBoundaries();

        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const label end = start + procBoundaries[patchI].patchSize();

            # ifdef USE_OMP
            # pragma omp parallel for schedule(guided) \
            reduction(+ : nDecomposedFaces)
            # endif
            for(label faceI=start;faceI<end;++faceI)
            {
                const face& f = faces[faceI];

                label nFeatureEdges(0);

                forAll(f, eI)
                {
                    const edge e = f.faceEdge(eI);

                    const label bs = bp[e[0]];
                    const label be = bp[e[1]];
                    if( (bs != -1) && (be != -1) )
                    {
                        //- check if this edge is a boundary edge
                        //- and a feature edge
                        forAllRow(pointEdges, bs, i)
                        {
                            const label beI = pointEdges(bs, i);

                            if( (edges[beI] == e) && featureBndEdge[beI] )
                                ++nFeatureEdges;
                        }
                    }
                }

                if( nFeatureEdges > 1 )
                {
                    ++nDecomposedFaces;
                    decomposeFace[faceI] = true;
                    decomposeCell_[owner[faceI]] = true;
                }
            }
        }
    }

    reduce(nDecomposedFaces, sumOp<label>());

    if( nDecomposedFaces != 0 )
    {
        Info << nDecomposedFaces << " faces decomposed into triangles" << endl;

        decompose_ = true;
        decomposeFaces df(mesh_);
        df.decomposeMeshFaces(decomposeFace);

        clearMeshSurface();
        mesh_.clearAddressingData();
    }

    Info << "Finished decomposing problematic faces" << endl;
}

void correctEdgesBetweenPatches::decomposeConcaveFaces()
{
    const meshSurfaceEngine& mse = meshSurface();
    const labelList& bPoints = mse.boundaryPoints();
    const labelList& bp = mse.bp();
    const edgeList& edges = mse.edges();
    const VRWGraph& edgeFaces = mse.edgeFaces();
    const VRWGraph& bpEdges = mse.boundaryPointEdges();
    const labelList& facePatch = mse.boundaryFacePatches();

    //- classify edges at the surface
    meshSurfaceCheckEdgeTypes edgeChecker(mse);
    const List<direction>& edgeType = edgeChecker.edgeTypes();

    //- find concave points
    boolList concavePoint(bPoints.size(), false);

    labelList edgeInPatch(edges.size(), -1);

    direction problematicTypes = 0;
    problematicTypes |= meshSurfaceCheckEdgeTypes::CONCAVEEDGE;
    problematicTypes |= meshSurfaceCheckEdgeTypes::UNDETERMINED;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(edgeType, eI)
    {
        if( edgeType[eI] & meshSurfaceCheckEdgeTypes::PATCHEDGE )
        {
            if( edgeFaces.sizeOfRow(eI) )
                edgeInPatch[eI] = facePatch[edgeFaces(eI, 0)];

            continue;
        }

        if( edgeType[eI] & problematicTypes)
        {
            const edge& e = edges[eI];

            concavePoint[bp[e.start()]] = true;
            concavePoint[bp[e.end()]] = true;
        }
    }

    if( Pstream::parRun() )
    {
        const Map<label>& globalToLocal =
            mse.globalToLocalBndEdgeAddressing();
        const VRWGraph& beAtProcs = mse.beAtProcs();

        std::map<label, labelLongList> exchangeData;
        forAll(mse.beNeiProcs(), i)
            exchangeData[mse.beNeiProcs()[i]].clear();

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label beI = it();

            if( edgeInPatch[beI] < 0 )
                continue;

            forAllRow(beAtProcs, beI, i)
            {
                const label neiProc = beAtProcs(beI, i);

                if( neiProc == Pstream::myProcNo() )
                    continue;

                labelLongList& dts = exchangeData[neiProc];

                dts.append(it.key());
                dts.append(edgeInPatch[beI]);
            }
        }

        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        for(label i=0;i<receivedData.size();)
        {
            const label geI = receivedData[i++];
            const label patchI = receivedData[i++];
            if( !globalToLocal.found(geI) ) continue;
            const label beI = globalToLocal[geI];
            if( edgeInPatch[beI] == -1 )
            {
                edgeInPatch[beI] = patchI;
            }
            else if( edgeInPatch[beI] != patchI )
            {
                FatalErrorIn
                (
                    "void correctEdgesBetweenPatches::decomposeConcaveFaces()"
                ) << "Invalid patch!" << abort(FatalError);
            }
        }
    }

    //- decompose internal faces attached to concave vertices which have two
    //- or more edges at the boundary
    const faceListPMG& faces = mesh_.faces();
    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();

    boolList decomposeFace(faces.size(), false);
    label nDecomposed(0);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100) \
    reduction(+ : nDecomposed)
    # endif
    for(label faceI=0;faceI<mesh_.nInternalFaces();++faceI)
    {
        const face& f = faces[faceI];

        bool hasConcave(false);
        label nBndEdges(0);
        DynList<label> bndEdgePatches;

        forAll(f, pI)
        {
            const label bpI = bp[f[pI]];

            if( bpI < 0 )
                continue;

            if( concavePoint[bpI] )
                hasConcave = true;

            //- points is at a concave edge
            //- count the number of boundary edge
            const edge e = f.faceEdge(pI);

            forAllRow(bpEdges, bpI, bpeI)
            {
                const label beI = bpEdges(bpI, bpeI);
                const edge& ee = edges[beI];

                if( e == ee )
                {
                    ++nBndEdges;
                    bndEdgePatches.appendIfNotIn(edgeInPatch[beI]);
                    break;
                }
            }
        }

        if( hasConcave && (nBndEdges > 1) && (bndEdgePatches.size() > 1) )
        {
            //- the face has two or more edges at the boundary
            //- Hence, it is marked for decomposition
            decomposeFace[faceI] = true;

            decomposeCell_[owner[faceI]] = true;
            decomposeCell_[neighbour[faceI]] = true;

            ++nDecomposed;
        }
    }

    //- finally, perform decomposition of marked faces
    if( returnReduce(nDecomposed, sumOp<label>()) != 0 )
    {
        Info << "Decomposing " << nDecomposed << " internal faces" << endl;
        decomposeFaces(mesh_).decomposeMeshFaces(decomposeFace);

        decompose_ = true;

        clearMeshSurface();
        mesh_.clearAddressingData();
    }
}

void correctEdgesBetweenPatches::patchCorrection()
{
    Info << "Performing patch correction" << endl;
    newPatchCorrectionPoints_.clear();

    const meshSurfaceEngine& mse = meshSurface();

    meshSurfacePartitioner surfacePartitioner(mse);

    const labelList& bPoints = mse.boundaryPoints();
    const VRWGraph& faceEdges = mse.faceEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& bp = mse.bp();
    const labelList& boundaryFaceOwners = mse.faceOwners();
    const labelList& facePatches = mse.boundaryFacePatches();

    //- set flag 1 to corner vertices, flag 2 to edge vertices
    List<direction> nodeType(bPoints.size(), direction(0));

    //- set corner flags
    const labelHashSet& corners = surfacePartitioner.corners();
    forAllConstIter(labelHashSet, corners, it)
        nodeType[it.key()] |= 1;

    //- set flgs to edge vertices
    const labelHashSet& edgePoints = surfacePartitioner.edgePoints();
    forAllConstIter(labelHashSet, edgePoints, it)
        nodeType[it.key()] |= 2;

    //- set flags for feature edges
    boolList featureEdge(edgeFaces.size(), false);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided)
    # endif
    forAll(edgeFaces, eI)
    {
        if( edgeFaces.sizeOfRow(eI) != 2 )
            continue;

        if( facePatches[edgeFaces(eI, 0)] != facePatches[edgeFaces(eI, 1)] )
            featureEdge[eI] = true;
    }

    if( Pstream::parRun() )
    {
        //- set flags for edges at parallel boundaries
        const Map<label>& otherProcPatch = mse.otherEdgeFacePatch();

        forAllConstIter(Map<label>, otherProcPatch, it)
            if( facePatches[edgeFaces(it.key(), 0)] != it() )
                featureEdge[it.key()] = true;
    }

    //- decompose bad faces into triangles
    newBoundaryFaces_.clear();
    newBoundaryOwners_.clear();
    newBoundaryPatches_.clear();
    face triF(3);

    label nDecomposedFaces(0);

    // Diagnostic only: quantify whether face::centre() is a
    // geometrically admissible fan apex for >4-sided patch faces.
    label nFanApex = 0;
    label nFanApexOutsideAabb = 0;
    label worstFanApexBfI = -1;
    scalar maxFanApexOffsetRatio = -1.0;

    // Diagnostic only: retain provenance of every fan decomposition
    // so any post-patchCorrection negative owner can be tied directly
    // back to the boundary face/apex which created it.
    DynamicList<label> fanDiagBfI;
    DynamicList<label> fanDiagOwner;
    DynamicList<label> fanDiagPatch;
    DynamicList<label> fanDiagNVerts;
    DynamicList<point> fanDiagApex;

    // Diagnostic only: provenance for legacy 4-sided split.
    DynamicList<label> quadDiagBfI;
    DynamicList<label> quadDiagOwner;
    DynamicList<label> quadDiagPatch;
    DynamicList<label> quadDiagCorner;
    DynamicList<label> quadDiagA;
    DynamicList<label> quadDiagB;

    // General patch-correction transaction state.
    //
    // 0 = unchanged
    // 1 = >4-sided arithmetic-apex fan
    // 2 = quad split
    labelList patchCorrectionType(bFaces.size(), 0);
    labelList quadCorrectionCorner(bFaces.size(), -1);
    labelList fanCorrectionApexPoint(bFaces.size(), -1);

    // 0 = legacy i--i+2 diagonal
    // 1 = alternate (i+1)--(i-1) diagonal
    // -1 = preserve original quad
    labelList quadCorrectionChoice(bFaces.size(), 0);

    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];

        bool store(true);

        forAll(bf, i)
        {
            if(
                (nodeType[bp[bf[i]]] == direction(2)) &&
                featureEdge[faceEdges(bfI, i)] &&
                featureEdge[faceEdges(bfI, bf.rcIndex(i))]
            )
            {
                if( bf.size() > 4 )
                {
                    store = false;
                    ++nDecomposedFaces;
                    decomposeCell_[boundaryFaceOwners[bfI]] = true;
                    decompose_ = true;

                    //- decompose into triangles
                    const point p = bf.centre(mesh_.points());

                    // Diagnostic only: compare the OpenFOAM
                    // area-weighted face centre with the raw vertex
                    // cloud of this polygon.
                    point pAvg = point::zero;

                    scalar minX = GREAT;
                    scalar minY = GREAT;
                    scalar minZ = GREAT;
                    scalar maxX = -GREAT;
                    scalar maxY = -GREAT;
                    scalar maxZ = -GREAT;
                    scalar maxEdgeLen = 0.0;

                    forAll(bf, fpI)
                    {
                        const point& fp =
                            mesh_.points()[bf[fpI]];

                        pAvg += fp;

                        minX = Foam::min(minX, fp.x());
                        minY = Foam::min(minY, fp.y());
                        minZ = Foam::min(minZ, fp.z());

                        maxX = Foam::max(maxX, fp.x());
                        maxY = Foam::max(maxY, fp.y());
                        maxZ = Foam::max(maxZ, fp.z());

                        const point& fpNext =
                            mesh_.points()
                            [
                                bf[bf.fcIndex(fpI)]
                            ];

                        maxEdgeLen =
                            Foam::max
                            (
                                maxEdgeLen,
                                mag(fpNext - fp)
                            );
                    }

                    pAvg /= scalar(bf.size());

                    const point aabbMin(minX, minY, minZ);
                    const point aabbMax(maxX, maxY, maxZ);

                    const scalar aabbDiag =
                        mag(aabbMax - aabbMin);

                    const scalar localScale =
                        Foam::max
                        (
                            Foam::max(aabbDiag, maxEdgeLen),
                            VSMALL
                        );

                    const scalar centreOffset =
                        mag(p - pAvg);

                    const scalar centreOffsetRatio =
                        centreOffset/localScale;

                    const bool outsideAabb =
                    (
                        p.x() < minX - SMALL ||
                        p.x() > maxX + SMALL ||
                        p.y() < minY - SMALL ||
                        p.y() > maxY + SMALL ||
                        p.z() < minZ - SMALL ||
                        p.z() > maxZ + SMALL
                    );

                    ++nFanApex;

                    if( outsideAabb )
                    {
                        ++nFanApexOutsideAabb;
                    }

                    if
                    (
                        centreOffsetRatio >
                        maxFanApexOffsetRatio
                    )
                    {
                        maxFanApexOffsetRatio =
                            centreOffsetRatio;

                        worstFanApexBfI = bfI;
                    }

                    if
                    (
                        outsideAabb ||
                        centreOffsetRatio > 0.25
                    )
                    {
                        Info
                            << "[PATCH_APEX_DIAG]"
                            << " bfI=" << bfI
                            << " owner="
                            << boundaryFaceOwners[bfI]
                            << " patch="
                            << facePatches[bfI]
                            << " nVerts=" << bf.size()
                            << " pAvg=" << pAvg
                            << " pCtr=" << p
                            << " aabbMin=" << aabbMin
                            << " aabbMax=" << aabbMax
                            << " aabbDiag=" << aabbDiag
                            << " maxEdgeLen=" << maxEdgeLen
                            << " centreOffset="
                            << centreOffset
                            << " offsetRatio="
                            << centreOffsetRatio
                            << " outsideAabb="
                            << outsideAabb
                            << endl;
                    }

                    fanDiagBfI.append(bfI);
                    fanDiagOwner.append(boundaryFaceOwners[bfI]);
                    fanDiagPatch.append(facePatches[bfI]);
                    fanDiagNVerts.append(bf.size());
                    fanDiagApex.append(pAvg);

                    triF[2] = mesh_.points().size();

                    // Use the arithmetic vertex mean as the actual
                    // fan apex. Unlike face::centre(), this remains
                    // bounded by the local vertex cloud and cannot
                    // fly hundreds of local face lengths away when
                    // projected triangle-area contributions cancel.
                    mesh_.points().append(pAvg);

                    patchCorrectionType[bfI] = 1;
                    fanCorrectionApexPoint[bfI] = triF[2];

                    // Record the new apex for the existing targeted
                    // post-topology surface projection.
                    newPatchCorrectionPoints_.append(triF[2]);

                    static label nCentroid = 0;
                    if( ++nCentroid <= 10 )
                    {
                        Info
                            << "[GeomFix] new centroid point "
                            << triF[2]
                            << " at " << pAvg
                            << " legacyFaceCentre=" << p
                            << " patch="
                            << facePatches[bfI]
                            << endl;
                    }

                    forAll(bf, j)
                    {
                        triF[0] = bf[j];
                        triF[1] = bf.nextLabel(j);

                        newBoundaryFaces_.appendList(triF);
                        newBoundaryOwners_.append(boundaryFaceOwners[bfI]);
                        newBoundaryPatches_.append(facePatches[bfI]);
                    }

                    break;
                }
                else if( bf.size() == 4 )
                {
                    store = false;

                    // Defer the actual diagonal choice until all corrected
                    // quads belonging to the same owner cell are known.
                    // Choices belonging to one owner are evaluated together
                    // against the exact prospective owner-cell volume.
                    patchCorrectionType[bfI] = 2;
                    quadCorrectionCorner[bfI] = i;

                    // Retain provenance while this algorithm is being
                    // validated.
                    quadDiagBfI.append(bfI);
                    quadDiagOwner.append(boundaryFaceOwners[bfI]);
                    quadDiagPatch.append(facePatches[bfI]);
                    quadDiagCorner.append(i);
                    quadDiagA.append(bf[i]);
                    quadDiagB.append(bf[(i+2)%4]);

                    break;
                }
            }
        }

        if( store )
        {
            //- face has not been altered
            newBoundaryFaces_.appendList(bf);
            newBoundaryOwners_.append(boundaryFaceOwners[bfI]);
            newBoundaryPatches_.append(facePatches[bfI]);
        }
    }


    // -----------------------------------------------------------------
    // General quad-diagonal transaction.
    //
    // Boundary replacement reconstructs the physical boundary faces in
    // every affected owner cell.  Therefore choose all quad diagonals
    // belonging to an owner together, before replaceBoundary() mutates
    // anything.
    //
    // The prospective volume below duplicates the face-centre/face-area
    // representation and raw signed-volume arithmetic used by
    // polyMeshGenChecks::checkCellVolumes().
    // -----------------------------------------------------------------

    const faceListPMG& allFaces = mesh_.faces();
    const cellListPMG& allCells = mesh_.cells();
    const labelList& allOwners = mesh_.owner();
    const pointFieldPMG& allPoints = mesh_.points();

    const label boundaryStart =
        mesh_.boundaries()[0].patchStart();

    const label boundaryEnd =
        boundaryStart + bFaces.size();


    auto calcVirtualFaceGeometry =
    [&]
    (
        const face& vf,
        point& fc,
        vector& fa
    )
    {
        const label nPoints = vf.size();

        if( nPoints == 3 )
        {
            fc =
                (1.0/3.0)
               *(
                    allPoints[vf[0]]
                  + allPoints[vf[1]]
                  + allPoints[vf[2]]
                );

            fa =
                0.5
               *(
                    (allPoints[vf[1]] - allPoints[vf[0]])
                  ^ (allPoints[vf[2]] - allPoints[vf[0]])
                );

            return;
        }

        vector sumN = vector::zero;
        scalar sumA = 0.0;
        vector sumAc = vector::zero;

        point fCentre = allPoints[vf[0]];

        for(label pi=1; pi<nPoints; ++pi)
            fCentre += allPoints[vf[pi]];

        fCentre /= scalar(nPoints);

        for(label pi=0; pi<nPoints; ++pi)
        {
            const point& nextPoint =
                allPoints[vf[(pi+1)%nPoints]];

            const vector c =
                allPoints[vf[pi]]
              + nextPoint
              + fCentre;

            const vector n =
                (nextPoint - allPoints[vf[pi]])
              ^ (fCentre - allPoints[vf[pi]]);

            const scalar a = mag(n);

            sumN += n;
            sumA += a;
            sumAc += a*c;
        }

        fc =
            (1.0/3.0)
           *sumAc/(sumA + VSMALL);

        fa = 0.5*sumN;
    };


    auto appendFaceGeometry =
    [&]
    (
        const face& vf,
        const bool reverseOrientation,
        DynamicList<point>& centres,
        DynamicList<vector>& areas
    )
    {
        point fc(point::zero);
        vector fa(vector::zero);

        calcVirtualFaceGeometry(vf, fc, fa);

        if( reverseOrientation )
            fa = -fa;

        centres.append(fc);
        areas.append(fa);
    };


    auto prospectiveOwnerVolume =
    [&]
    (
        const label ownerCell,
        const DynamicList<label>& ownerQuadBfI,
        const label choiceMask,
        const bool preserveAllQuads
    ) -> scalar
    {
        DynamicList<point> centres;
        DynamicList<vector> areas;

        const cell& c = allCells[ownerCell];

        forAll(c, cfI)
        {
            const label faceI = c[cfI];

            const bool physicalBoundary =
            (
                faceI >= boundaryStart
             && faceI < boundaryEnd
            );

            if( !physicalBoundary )
            {
                appendFaceGeometry
                (
                    allFaces[faceI],
                    allOwners[faceI] != ownerCell,
                    centres,
                    areas
                );

                continue;
            }


            const label bfI =
                faceI - boundaryStart;

            const face& bf = bFaces[bfI];

            if( patchCorrectionType[bfI] == 0 )
            {
                appendFaceGeometry
                (
                    bf,
                    false,
                    centres,
                    areas
                );

                continue;
            }


            if( patchCorrectionType[bfI] == 1 )
            {
                const label apexPoint =
                    fanCorrectionApexPoint[bfI];

                if( apexPoint < 0 )
                    return -GREAT;

                face tf(3);

                forAll(bf, j)
                {
                    tf[0] = bf[j];
                    tf[1] = bf.nextLabel(j);
                    tf[2] = apexPoint;

                    appendFaceGeometry
                    (
                        tf,
                        false,
                        centres,
                        areas
                    );
                }

                continue;
            }


            // Quad correction.
            if( preserveAllQuads )
            {
                appendFaceGeometry
                (
                    bf,
                    false,
                    centres,
                    areas
                );

                continue;
            }


            label localQuadI = -1;

            forAll(ownerQuadBfI, qI)
            {
                if( ownerQuadBfI[qI] == bfI )
                {
                    localQuadI = qI;
                    break;
                }
            }

            if( localQuadI < 0 )
                return -GREAT;


            const label i =
                quadCorrectionCorner[bfI];

            if( i < 0 )
                return -GREAT;


            const bool alternate =
                (choiceMask & (label(1) << localQuadI));

            face t0(3);
            face t1(3);

            if( !alternate )
            {
                // Legacy diagonal i -- i+2.
                t0[0] = bf[i];
                t0[1] = bf.nextLabel(i);
                t0[2] = bf[(i+2)%4];

                t1[0] = bf[i];
                t1[1] = bf[(i+2)%4];
                t1[2] = bf.prevLabel(i);
            }
            else
            {
                // Alternate diagonal (i+1) -- (i-1).
                t0[0] = bf[i];
                t0[1] = bf.nextLabel(i);
                t0[2] = bf.prevLabel(i);

                t1[0] = bf.nextLabel(i);
                t1[1] = bf[(i+2)%4];
                t1[2] = bf.prevLabel(i);
            }

            appendFaceGeometry
            (
                t0,
                false,
                centres,
                areas
            );

            appendFaceGeometry
            (
                t1,
                false,
                centres,
                areas
            );
        }


        if( centres.size() == 0 )
            return -GREAT;


        point cEst(point::zero);

        forAll(centres, fI)
            cEst += centres[fI];

        cEst /= scalar(centres.size());


        scalar cellVol = 0.0;

        forAll(centres, fI)
        {
            cellVol +=
                areas[fI] & (centres[fI] - cEst);
        }

        return cellVol/3.0;
    };


    labelHashSet quadOwnerSet;

    forAll(quadDiagOwner, qI)
        quadOwnerSet.insert(quadDiagOwner[qI]);


    label nQuadOwnerGroups = 0;
    label nQuadLegacySafeGroups = 0;
    label nQuadAlternateGroups = 0;
    label nQuadAlternateFaces = 0;
    label nQuadPreservedGroups = 0;

    const label maxEnumeratedQuadsPerOwner = 8;


    forAllConstIter(labelHashSet, quadOwnerSet, ownIt)
    {
        const label ownerCell = ownIt.key();

        DynamicList<label> ownerQuadBfI;

        forAll(quadDiagOwner, qI)
        {
            if( quadDiagOwner[qI] == ownerCell )
                ownerQuadBfI.append(quadDiagBfI[qI]);
        }

        const label nQuads =
            ownerQuadBfI.size();

        if( nQuads == 0 )
            continue;

        ++nQuadOwnerGroups;


        if( nQuads > maxEnumeratedQuadsPerOwner )
        {
            const scalar preserveVol =
                prospectiveOwnerVolume
                (
                    ownerCell,
                    ownerQuadBfI,
                    0,
                    true
                );

            if( preserveVol < VSMALL )
            {
                FatalErrorIn
                (
                    "void correctEdgesBetweenPatches::patchCorrection()"
                )
                    << "Cannot preserve positive owner cell "
                    << ownerCell
                    << " when quad transaction exceeds complexity guard."
                    << " nQuads=" << nQuads
                    << " preserveVol=" << preserveVol
                    << abort(FatalError);
            }

            forAll(ownerQuadBfI, qI)
            {
                const label bfI = ownerQuadBfI[qI];

                patchCorrectionType[bfI] = 0;
                quadCorrectionChoice[bfI] = -1;
            }

            ++nQuadPreservedGroups;

            Info
                << "[PATCH_QUAD_OWNER_TRANSACTION]"
                << " owner=" << ownerCell
                << " nQuads=" << nQuads
                << " action=preserve"
                << " reason=complexityGuard"
                << " preserveVol=" << preserveVol
                << endl;

            continue;
        }


        const scalar legacyVolume =
            prospectiveOwnerVolume
            (
                ownerCell,
                ownerQuadBfI,
                0,
                false
            );


        label selectedMask = -1;
        scalar selectedVolume = -GREAT;


        // Preserve historical behaviour whenever it is already valid.
        if( legacyVolume >= VSMALL )
        {
            selectedMask = 0;
            selectedVolume = legacyVolume;
            ++nQuadLegacySafeGroups;
        }
        else
        {
            const label nCombinations =
                label(1) << nQuads;

            for
            (
                label mask=1;
                mask<nCombinations;
                ++mask
            )
            {
                const scalar candidateVolume =
                    prospectiveOwnerVolume
                    (
                        ownerCell,
                        ownerQuadBfI,
                        mask,
                        false
                    );

                if
                (
                    candidateVolume >= VSMALL
                 && candidateVolume > selectedVolume
                )
                {
                    selectedMask = mask;
                    selectedVolume = candidateVolume;
                }
            }
        }


        if( selectedMask < 0 )
        {
            // Fail closed: retain the original quads rather than creating
            // a known non-positive owner cell.
            const scalar preserveVol =
                prospectiveOwnerVolume
                (
                    ownerCell,
                    ownerQuadBfI,
                    0,
                    true
                );

            if( preserveVol < VSMALL )
            {
                FatalErrorIn
                (
                    "void correctEdgesBetweenPatches::patchCorrection()"
                )
                    << "No positive quad decomposition or preservation "
                    << "exists for owner " << ownerCell
                    << ". legacyVol=" << legacyVolume
                    << " preserveVol=" << preserveVol
                    << abort(FatalError);
            }

            forAll(ownerQuadBfI, qI)
            {
                const label bfI = ownerQuadBfI[qI];

                patchCorrectionType[bfI] = 0;
                quadCorrectionChoice[bfI] = -1;
            }

            ++nQuadPreservedGroups;

            Info
                << "[PATCH_QUAD_OWNER_TRANSACTION]"
                << " owner=" << ownerCell
                << " nQuads=" << nQuads
                << " action=preserve"
                << " legacyVol=" << legacyVolume
                << " preserveVol=" << preserveVol
                << endl;

            continue;
        }


        label nAlternateThisOwner = 0;

        forAll(ownerQuadBfI, qI)
        {
            const label bfI = ownerQuadBfI[qI];

            const label choice =
                (
                    selectedMask
                  & (label(1) << qI)
                )
              ? 1
              : 0;

            quadCorrectionChoice[bfI] = choice;

            if( choice == 1 )
            {
                ++nAlternateThisOwner;
                ++nQuadAlternateFaces;
            }

            ++nDecomposedFaces;
        }

        decomposeCell_[ownerCell] = true;
        decompose_ = true;


        if( selectedMask != 0 )
        {
            ++nQuadAlternateGroups;

            Info
                << "[PATCH_QUAD_OWNER_TRANSACTION]"
                << " owner=" << ownerCell
                << " nQuads=" << nQuads
                << " action=alternate"
                << " selectedMask=" << selectedMask
                << " alternateFaces=" << nAlternateThisOwner
                << " legacyVol=" << legacyVolume
                << " selectedVol=" << selectedVolume
                << endl;
        }
    }


    Info
        << "[PATCH_QUAD_TRANSACTION_SUMMARY]"
        << " ownerGroups=" << nQuadOwnerGroups
        << " legacySafeGroups=" << nQuadLegacySafeGroups
        << " alternateGroups=" << nQuadAlternateGroups
        << " alternateFaces=" << nQuadAlternateFaces
        << " preservedGroups=" << nQuadPreservedGroups
        << endl;


    // -----------------------------------------------------------------
    // Rebuild the complete proposed boundary exactly once using the
    // owner-cell-selected quad choices.
    // -----------------------------------------------------------------

    newBoundaryFaces_.clear();
    newBoundaryOwners_.clear();
    newBoundaryPatches_.clear();


    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];

        if( patchCorrectionType[bfI] == 0 )
        {
            newBoundaryFaces_.appendList(bf);
            newBoundaryOwners_.append
            (
                boundaryFaceOwners[bfI]
            );
            newBoundaryPatches_.append
            (
                facePatches[bfI]
            );

            continue;
        }


        if( patchCorrectionType[bfI] == 1 )
        {
            const label apexPoint =
                fanCorrectionApexPoint[bfI];

            if( apexPoint < 0 )
            {
                FatalErrorIn
                (
                    "void correctEdgesBetweenPatches::patchCorrection()"
                )
                    << "Missing fan apex for boundary face "
                    << bfI
                    << abort(FatalError);
            }

            face tf(3);
            tf[2] = apexPoint;

            forAll(bf, j)
            {
                tf[0] = bf[j];
                tf[1] = bf.nextLabel(j);

                newBoundaryFaces_.appendList(tf);
                newBoundaryOwners_.append
                (
                    boundaryFaceOwners[bfI]
                );
                newBoundaryPatches_.append
                (
                    facePatches[bfI]
                );
            }

            continue;
        }


        const label i =
            quadCorrectionCorner[bfI];

        const label choice =
            quadCorrectionChoice[bfI];

        if( i < 0 || (choice != 0 && choice != 1) )
        {
            FatalErrorIn
            (
                "void correctEdgesBetweenPatches::patchCorrection()"
            )
                << "Invalid committed quad state for bfI="
                << bfI
                << " corner=" << i
                << " choice=" << choice
                << abort(FatalError);
        }


        face t0(3);
        face t1(3);

        if( choice == 0 )
        {
            t0[0] = bf[i];
            t0[1] = bf.nextLabel(i);
            t0[2] = bf[(i+2)%4];

            t1[0] = bf[i];
            t1[1] = bf[(i+2)%4];
            t1[2] = bf.prevLabel(i);
        }
        else
        {
            t0[0] = bf[i];
            t0[1] = bf.nextLabel(i);
            t0[2] = bf.prevLabel(i);

            t1[0] = bf.nextLabel(i);
            t1[1] = bf[(i+2)%4];
            t1[2] = bf.prevLabel(i);
        }

        newBoundaryFaces_.appendList(t0);
        newBoundaryOwners_.append
        (
            boundaryFaceOwners[bfI]
        );
        newBoundaryPatches_.append
        (
            facePatches[bfI]
        );

        newBoundaryFaces_.appendList(t1);
        newBoundaryOwners_.append
        (
            boundaryFaceOwners[bfI]
        );
        newBoundaryPatches_.append
        (
            facePatches[bfI]
        );
    }


    Info
        << "[PATCH_APEX_SUMMARY]"
        << " fanApexCount=" << nFanApex
        << " outsideAabb=" << nFanApexOutsideAabb
        << " maxOffsetRatio=" << maxFanApexOffsetRatio
        << " worstBfI=" << worstFanApexBfI
        << endl;

    reduce(decompose_, maxOp<bool>());

    if( returnReduce(nDecomposedFaces, sumOp<label>()) != 0 )
    {
        replaceBoundary();
        clearMeshSurface();
        mesh_.clearAddressingData();

        labelHashSet patchNegVolCells;

        polyMeshGenChecks::checkCellVolumes
        (
            mesh_,
            false,
            &patchNegVolCells
        );

        label matchedFanOwners = 0;
        label matchedFanRecords = 0;
        label matchedQuadOwners = 0;
        label matchedQuadRecords = 0;

        labelHashSet matchedAnyOwner;

        forAll(fanDiagOwner, i)
        {
            if( patchNegVolCells.found(fanDiagOwner[i]) )
            {
                ++matchedFanRecords;
                matchedAnyOwner.insert(fanDiagOwner[i]);

                bool firstForOwner = true;

                for(label j=0; j<i; ++j)
                {
                    if
                    (
                        fanDiagOwner[j] == fanDiagOwner[i]
                     && patchNegVolCells.found(fanDiagOwner[j])
                    )
                    {
                        firstForOwner = false;
                        break;
                    }
                }

                if( firstForOwner )
                {
                    ++matchedFanOwners;
                }

                Info
                    << "[PATCH_RESIDUAL_FAN_OWNER]"
                    << " cell=" << fanDiagOwner[i]
                    << " bfI=" << fanDiagBfI[i]
                    << " patch=" << fanDiagPatch[i]
                    << " nVerts=" << fanDiagNVerts[i]
                    << " apex=" << fanDiagApex[i]
                    << endl;
            }
        }

        labelHashSet matchedQuadOwnerSet;

        forAll(quadDiagOwner, i)
        {
            if( patchNegVolCells.found(quadDiagOwner[i]) )
            {
                ++matchedQuadRecords;
                matchedQuadOwnerSet.insert(quadDiagOwner[i]);
                matchedAnyOwner.insert(quadDiagOwner[i]);

                Info
                    << "[PATCH_RESIDUAL_QUAD_OWNER]"
                    << " cell=" << quadDiagOwner[i]
                    << " bfI=" << quadDiagBfI[i]
                    << " patch=" << quadDiagPatch[i]
                    << " cornerI=" << quadDiagCorner[i]
                    << " diagonalPoints="
                    << quadDiagA[i] << "|" << quadDiagB[i]
                    << endl;
            }
        }

        matchedQuadOwners = matchedQuadOwnerSet.size();

        Info
            << "[PATCH_RESIDUAL_OWNER_SUMMARY]"
            << " negVol=" << patchNegVolCells.size()
            << " fanRecords=" << fanDiagOwner.size()
            << " matchedFanOwners=" << matchedFanOwners
            << " matchedFanRecords=" << matchedFanRecords
            << " quadRecords=" << quadDiagOwner.size()
            << " matchedQuadOwners=" << matchedQuadOwners
            << " matchedQuadRecords=" << matchedQuadRecords
            << " matchedAnyOwners=" << matchedAnyOwner.size()
            << " unmatchedNegOwners="
            << patchNegVolCells.size() - matchedAnyOwner.size()
            << endl;
    }

    Info << "Finished with patch correction" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
