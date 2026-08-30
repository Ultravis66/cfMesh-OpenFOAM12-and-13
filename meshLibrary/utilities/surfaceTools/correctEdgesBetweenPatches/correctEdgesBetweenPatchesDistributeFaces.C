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
                    ++nDecomposedFaces;
                    decomposeCell_[boundaryFaceOwners[bfI]] = true;
                    decompose_ = true;

                    // Diagnostic only: the legacy quad split chooses
                    // the diagonal from the detected feature-corner
                    // vertex i to the opposite vertex i+2.
                    quadDiagBfI.append(bfI);
                    quadDiagOwner.append(boundaryFaceOwners[bfI]);
                    quadDiagPatch.append(facePatches[bfI]);
                    quadDiagCorner.append(i);
                    quadDiagA.append(bf[i]);
                    quadDiagB.append(bf[(i+2)%4]);

                    // Forensic experiment: the residual negative
                    // owners were traced exactly to these three quad
                    // decompositions. Test the alternate diagonal only
                    // on those faces.
                    const bool useAlternateDiagonal =
                    (
                        bfI == 124493
                     || bfI == 219116
                     || bfI == 511093
                    );

                    if( useAlternateDiagonal )
                    {
                        Info
                            << "[QUAD_ALT_FORENSIC]"
                            << " bfI=" << bfI
                            << " owner=" << boundaryFaceOwners[bfI]
                            << " patch=" << facePatches[bfI]
                            << " legacyDiagonal="
                            << bf[i] << "|" << bf[(i+2)%4]
                            << " alternateDiagonal="
                            << bf.nextLabel(i)
                            << "|" << bf.prevLabel(i)
                            << endl;

                        // Alternate diagonal:
                        // (i+1) -------- (i-1)
                        //
                        // Preserve original polygon orientation.
                        triF[0] = bf[i];
                        triF[1] = bf.nextLabel(i);
                        triF[2] = bf.prevLabel(i);

                        newBoundaryFaces_.appendList(triF);
                        newBoundaryOwners_.append
                        (
                            boundaryFaceOwners[bfI]
                        );
                        newBoundaryPatches_.append
                        (
                            facePatches[bfI]
                        );

                        triF[0] = bf.nextLabel(i);
                        triF[1] = bf[(i+2)%4];
                        triF[2] = bf.prevLabel(i);

                        newBoundaryFaces_.appendList(triF);
                        newBoundaryOwners_.append
                        (
                            boundaryFaceOwners[bfI]
                        );
                        newBoundaryPatches_.append
                        (
                            facePatches[bfI]
                        );
                    }
                    else
                    {
                        // Legacy diagonal i -- i+2.
                        triF[0] = bf[i];

                        triF[1] = bf.nextLabel(i);
                        triF[2] = bf[(i+2)%4];
                        newBoundaryFaces_.appendList(triF);
                        newBoundaryOwners_.append
                        (
                            boundaryFaceOwners[bfI]
                        );
                        newBoundaryPatches_.append
                        (
                            facePatches[bfI]
                        );

                        triF[1] = bf[(i+2)%4];
                        triF[2] = bf.prevLabel(i);
                        newBoundaryFaces_.appendList(triF);
                        newBoundaryOwners_.append
                        (
                            boundaryFaceOwners[bfI]
                        );
                        newBoundaryPatches_.append
                        (
                            facePatches[bfI]
                        );
                    }

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
