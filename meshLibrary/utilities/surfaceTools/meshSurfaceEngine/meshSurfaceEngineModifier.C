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

#include "meshSurfaceEngineModifier.H"
#include "polyMeshGenModifier.H"
#include "polyMeshGenAddressing.H"
#include "demandDrivenData.H"

#include "labelledPoint.H"
#include "helperFunctionsPar.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshSurfaceEngineModifier::meshSurfaceEngineModifier
(
    meshSurfaceEngine& surfaceEngine
)
:
    surfaceEngine_(surfaceEngine)
{}

meshSurfaceEngineModifier::meshSurfaceEngineModifier
(
    const meshSurfaceEngine& surfaceEngine
)
:
    surfaceEngine_(const_cast<meshSurfaceEngine&>(surfaceEngine))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshSurfaceEngineModifier::~meshSurfaceEngineModifier()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar meshSurfaceEngineModifier::rawCellVolumeWithPoint
(
    const label cellI,
    const label movedGlobalPtI,
    const point& candidate,
    const bool substituteCandidate
) const
{
    const polyMeshGen& mesh = surfaceEngine_.mesh_;

    const pointFieldPMG& points = mesh.points();
    const faceListPMG& faces = mesh.faces();
    const cellListPMG& cells = mesh.cells();
    const labelList& owner = mesh.owner();

    const cell& c = cells[cellI];

    List<point> localFCentres(c.size());
    List<vector> localFAreas(c.size());

    vector cEst = vector::zero;

    forAll(c, cfI)
    {
        const label faceI = c[cfI];
        const face& f = faces[faceI];

        auto pt =
        [&]
        (
            const label ptI
        ) -> point
        {
            if
            (
                substituteCandidate &&
                ptI == movedGlobalPtI
            )
            {
                return candidate;
            }

            return point(points[ptI]);
        };

        point fCtr;
        vector fArea;

        const label nPoints = f.size();

        if( nPoints == 3 )
        {
            const point p0 = pt(f[0]);
            const point p1 = pt(f[1]);
            const point p2 = pt(f[2]);

            fCtr = (1.0/3.0)*(p0 + p1 + p2);
            fArea = 0.5*((p1 - p0)^(p2 - p0));
        }
        else
        {
            vector sumN = vector::zero;
            scalar sumA = 0.0;
            vector sumAc = vector::zero;

            point fCentre = pt(f[0]);

            for(label pi=1; pi<nPoints; ++pi)
            {
                fCentre += pt(f[pi]);
            }

            fCentre /= nPoints;

            for(label pi=0; pi<nPoints; ++pi)
            {
                const point curPoint = pt(f[pi]);
                const point nextPoint = pt(f.nextLabel(pi));

                const vector fc3 =
                    curPoint + nextPoint + fCentre;

                const vector n =
                    (nextPoint - curPoint)^
                    (fCentre - curPoint);

                const scalar a = mag(n);

                sumN += n;
                sumA += a;
                sumAc += a*fc3;
            }

            fCtr =
                (1.0/3.0)*sumAc/(sumA + VSMALL);

            fArea = 0.5*sumN;
        }

        localFCentres[cfI] = fCtr;
        localFAreas[cfI] = fArea;

        cEst += fCtr;
    }

    cEst /= scalar(c.size());

    scalar cellVol = 0.0;

    forAll(c, cfI)
    {
        const label faceI = c[cfI];

        scalar pyr3Vol =
            localFAreas[cfI] &
            (localFCentres[cfI] - cEst);

        if( owner[faceI] != cellI )
        {
            pyr3Vol *= -1.0;
        }

        cellVol += pyr3Vol;
    }

    return cellVol/3.0;
}


bool meshSurfaceEngineModifier::candidatePreservesPositiveCellVolumes
(
    const label bpI,
    const point& candidate,
    bool& touchesExistingBad,
    scalar& minPositiveRatio
) const
{
    touchesExistingBad = false;
    minPositiveRatio = GREAT;

    const polyMeshGen& mesh = surfaceEngine_.mesh_;
    const labelList& boundaryPoints =
        surfaceEngine_.boundaryPoints();

    if( bpI < 0 || bpI >= boundaryPoints.size() )
    {
        return false;
    }

    const label globalPtI = boundaryPoints[bpI];

    const VRWGraph& pointCells =
        mesh.addressingData().pointCells();

    if
    (
        globalPtI < 0 ||
        globalPtI >= pointCells.size() ||
        pointCells.sizeOfRow(globalPtI) == 0
    )
    {
        return false;
    }

    forAllRow(pointCells, globalPtI, pcI)
    {
        const label cellI =
            pointCells(globalPtI, pcI);

        const scalar oldVol =
            rawCellVolumeWithPoint
            (
                cellI,
                globalPtI,
                candidate,
                false
            );

        const scalar newVol =
            rawCellVolumeWithPoint
            (
                cellI,
                globalPtI,
                candidate,
                true
            );

        if( oldVol < VSMALL )
        {
            touchesExistingBad = true;
        }

        if
        (
            oldVol >= VSMALL &&
            newVol < VSMALL
        )
        {
            return false;
        }

        if
        (
            oldVol >= VSMALL &&
            newVol >= VSMALL
        )
        {
            minPositiveRatio =
                Foam::min
                (
                    minPositiveRatio,
                    newVol/(oldVol + VSMALL)
                );
        }
    }

    return true;
}


void meshSurfaceEngineModifier::moveBoundaryVertexNoUpdate
(
    const label bpI,
    const point& newP
)
{
    surfaceEngine_.mesh_.points()[surfaceEngine_.boundaryPoints()[bpI]] = newP;
}

void meshSurfaceEngineModifier::moveBoundaryVertex
(
    const label bpI,
    const point& newP
)
{
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    pointFieldPMG& points = surfaceEngine_.mesh_.points();
    points[bPoints[bpI]] = newP;

    if( surfaceEngine_.faceCentresPtr_ )
    {
        vectorField& faceCentres = *surfaceEngine_.faceCentresPtr_;
        const VRWGraph& pFaces = surfaceEngine_.pointFaces();
        const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();

        forAllRow(pFaces, bpI, pfI)
        {
            const label bfI = pFaces(bpI, pfI);

            {
                const face& bf = bFaces[bfI];
                if( bf.empty() )
                {
                    faceCentres[bfI] = vector::zero;
                    continue;
                }
                vector c = vector::zero;
                forAll(bf, pI) c += points[bf[pI]];
                faceCentres[bfI] = c / bf.size();
            }
        }
    }

    if( surfaceEngine_.faceNormalsPtr_ )
    {
        vectorField& faceNormals = *surfaceEngine_.faceNormalsPtr_;
        const VRWGraph& pFaces = surfaceEngine_.pointFaces();
        const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();

        forAllRow(pFaces, bpI, pfI)
        {
            const label bfI = pFaces(bpI, pfI);

            {
                const face& bf = bFaces[bfI];
                if( bf.size() < 3 )
                {
                    faceNormals[bfI] = vector::zero;
                    continue;
                }
                vector n = vector::zero;
                const point& p0 = points[bf[0]];
                for(label pI=1; pI<bf.size()-1; ++pI)
                    n += (points[bf[pI]]-p0)^(points[bf[pI+1]]-p0);
                faceNormals[bfI] = n;
            }
        }
    }

    if( surfaceEngine_.pointNormalsPtr_ && surfaceEngine_.faceNormalsPtr_ )
    {
        const vectorField& faceNormals = *surfaceEngine_.faceNormalsPtr_;
        const VRWGraph& pFaces = surfaceEngine_.pointFaces();
        const VRWGraph& pPoints = surfaceEngine_.pointPoints();

        vectorField& pn = *surfaceEngine_.pointNormalsPtr_;
        vector n(vector::zero);
        forAllRow(pFaces, bpI, pfI)
            n += faceNormals[pFaces(bpI, pfI)];

        const scalar l = mag(n);
        if( l > VSMALL )
        {
            n /= l;
        }
        else
        {
            n = vector::zero;
        }

        pn[bpI] = n;

        //- change normal of vertices connected to bpI
        forAllRow(pPoints, bpI, ppI)
        {
            const label bpJ = pPoints(bpI, ppI);
            n = vector::zero;
            forAllRow(pFaces, bpJ, pfI)
                n += faceNormals[pFaces(bpJ, pfI)];

            const scalar d = mag(n);
            if( d > VSMALL )
            {
                n /= d;
            }
            else
            {
                n = vector::zero;
            }

            pn[bpJ] = n;
        }
    }
}

void meshSurfaceEngineModifier::syncVerticesAtParallelBoundaries()
{
    if( !Pstream::parRun() )
        return;

    const Map<label>& globalToLocal =
        surfaceEngine_.globalToLocalBndPointAddressing();
    labelLongList syncNodes;
    forAllConstIter(Map<label>, globalToLocal, it)
        syncNodes.append(it());

    syncVerticesAtParallelBoundaries(syncNodes);
}

void meshSurfaceEngineModifier::syncVerticesAtParallelBoundaries
(
    const labelLongList& syncNodes
)
{
    if( !Pstream::parRun() )
        return;

    const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
    const labelList& globalLabel =
        surfaceEngine_.globalBoundaryPointLabel();
    const Map<label>& globalToLocal =
        surfaceEngine_.globalToLocalBndPointAddressing();
    const DynList<label>& neiProcs = surfaceEngine_.bpNeiProcs();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const pointFieldPMG& points = surfaceEngine_.mesh().points();

    std::map<label, LongList<labelledPoint> > exchangeData;
    forAll(neiProcs, i)
        exchangeData.insert
        (
            std::make_pair(neiProcs[i], LongList<labelledPoint>())
        );

    //- construct the map
    forAll(syncNodes, snI)
    {
        const label bpI = syncNodes[snI];

        if( bpAtProcs.sizeOfRow(bpI) == 0 )
            continue;

        point p = points[bPoints[bpI]] / bpAtProcs.sizeOfRow(bpI);
        moveBoundaryVertexNoUpdate(bpI, p);

        forAllRow(bpAtProcs, bpI, i)
        {
            const label neiProc = bpAtProcs(bpI, i);
            if( neiProc == Pstream::myProcNo() )
                continue;

            exchangeData[neiProc].append(labelledPoint(globalLabel[bpI], p));
        }
    }

    //- exchange the data with other processors
    LongList<labelledPoint> receivedData;
    help::exchangeMap(exchangeData, receivedData);

    //- adjust the coordinates
    forAll(receivedData, i)
    {
        const labelledPoint& lp = receivedData[i];
        if( !globalToLocal.found(lp.pointLabel()) )
            continue;
        if( !globalToLocal.found(lp.pointLabel()) ) continue;
        const label bpI = globalToLocal[lp.pointLabel()];
        const point newP = points[bPoints[bpI]] + lp.coordinates();
        moveBoundaryVertexNoUpdate(bpI, newP);
    }
}

void meshSurfaceEngineModifier::updateGeometry
(
    const labelLongList& updateBndNodes
)
{
    const pointFieldPMG& points = surfaceEngine_.points();
    const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();
    const VRWGraph& pFaces = surfaceEngine_.pointFaces();
    const labelList& bp = surfaceEngine_.bp();

    boolList updateFaces(bFaces.size(), false);
    // Serial: multiple bpI share faces -- updateFaces[...]=true races
    forAll(updateBndNodes, i)
    {
        const label bpI = updateBndNodes[i];
        forAllRow(pFaces, bpI, j)
            updateFaces[pFaces(bpI, j)] = true;
    }

    if( surfaceEngine_.faceCentresPtr_ )
    {
        vectorField& faceCentres = *surfaceEngine_.faceCentresPtr_;

        # ifdef USE_OMP
        # pragma omp parallel for if( updateFaces.size() > 1000 ) \
        schedule(dynamic, 100)
        # endif
        forAll(updateFaces, bfI)
        {
            if( updateFaces[bfI] )
            {
                const face& bf = bFaces[bfI];
                if( bf.empty() )
                {
                    faceCentres[bfI] = vector::zero;
                    continue;
                }
                vector c = vector::zero;
                forAll(bf, pI) c += points[bf[pI]];
                faceCentres[bfI] = c / bf.size();
            }
        }
    }

    if( surfaceEngine_.faceNormalsPtr_ )
    {
        vectorField& faceNormals = *surfaceEngine_.faceNormalsPtr_;

        # ifdef USE_OMP
        # pragma omp parallel for if( updateFaces.size() > 1000 ) \
        schedule(dynamic, 100)
        # endif
        forAll(updateFaces, bfI)
        {
            if( updateFaces[bfI] )
            {
                const face& bf = bFaces[bfI];
                if( bf.size() < 3 )
                {
                    faceNormals[bfI] = vector::zero;
                    continue;
                }
                vector n = vector::zero;
                const point& p0 = points[bf[0]];
                for(label pI=1; pI<bf.size()-1; ++pI)
                    n += (points[bf[pI]]-p0)^(points[bf[pI+1]]-p0);
                faceNormals[bfI] = n;
            }
        }
    }

    if( surfaceEngine_.pointNormalsPtr_ )
    {
        const vectorField& faceNormals = surfaceEngine_.faceNormals();

        boolList updateBndPoint(pFaces.size(), false);
        // Serial: multiple bpI share face points -- updateBndPoint[...]=true races
        forAll(updateBndNodes, i)
        {
            const label bpI = updateBndNodes[i];

            forAllRow(pFaces, bpI, pfI)
            {
                const face& bf = bFaces[pFaces(bpI, pfI)];

                forAll(bf, pI)
                {
                    const label bpJ = bp[bf[pI]];
                    if( bpJ >= 0 )
                        updateBndPoint[bpJ] = true;
                }
            }
        }

        vectorField& pn = *surfaceEngine_.pointNormalsPtr_;
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 100)
        # endif
        forAll(updateBndPoint, bpI)
        {
            if( !updateBndPoint[bpI] )
                continue;

            vector n(vector::zero);
            forAllRow(pFaces, bpI, pfI)
                n += faceNormals[pFaces(bpI, pfI)];

            const scalar l = mag(n);
            if( l > VSMALL )
            {
                n /= l;
            }
            else
            {
                n = vector::zero;
            }

            pn[bpI] = n;
        }

        if( Pstream::parRun() )
        {
            //- update point normals at inter-processor boundaries
            const Map<label>& globalToLocal =
                surfaceEngine_.globalToLocalBndPointAddressing();
            const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
            const DynList<label>& neiProcs = surfaceEngine_.bpNeiProcs();

            //- make sure that the points ar updated on all processors
            std::map<label, labelLongList> exchangeNodeLabels;
            forAll(neiProcs, i)
                exchangeNodeLabels[neiProcs[i]].clear();

            forAllConstIter(Map<label>, globalToLocal, it)
            {
                const label bpI = it();

                if( updateBndPoint[bpI] )
                {
                    forAllRow(bpAtProcs, bpI, i)
                    {
                        const label neiProc = bpAtProcs(bpI, i);

                        if( neiProc == Pstream::myProcNo() )
                            continue;

                        exchangeNodeLabels[neiProc].append(it.key());
                    }
                }
            }

            labelLongList receivedNodes;
            help::exchangeMap(exchangeNodeLabels, receivedNodes);

            forAll(receivedNodes, i)
            {
                if( !globalToLocal.found(receivedNodes[i]) )
                    continue;
                if( !globalToLocal.found(receivedNodes[i]) ) continue;
                updateBndPoint[globalToLocal[receivedNodes[i]]] = true;
            }


            //- start updating point normals
            std::map<label, LongList<labelledPoint> > exchangeData;
            forAll(neiProcs, i)
                exchangeData[neiProcs[i]].clear();

            //- prepare data for sending
            forAllConstIter(Map<label>, globalToLocal, iter)
            {
                const label bpI = iter();

                if( !updateBndPoint[bpI] )
                    continue;

                vector& n = pn[bpI];
                n = vector::zero;

                forAllRow(pFaces, bpI, pfI)
                    n += faceNormals[pFaces(bpI, pfI)];

                forAllRow(bpAtProcs, bpI, procI)
                {
                    const label neiProc = bpAtProcs(bpI, procI);
                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    exchangeData[neiProc].append(labelledPoint(iter.key(), n));
                }
            }

            //- exchange data with other procs
            LongList<labelledPoint> receivedData;
            help::exchangeMap(exchangeData, receivedData);

            forAll(receivedData, i)
            {
                if( !globalToLocal.found(receivedData[i].pointLabel()) )
                    continue;
                if( !globalToLocal.found(receivedData[i].pointLabel()) ) continue;
                const label bpI = globalToLocal[receivedData[i].pointLabel()];
                pn[bpI] += receivedData[i].coordinates();
            }

            //- normalize vectors
            forAllConstIter(Map<label>, globalToLocal, it)
            {
                const label bpI = it();

                if( !updateBndPoint[bpI] )
                    continue;

                vector normal = pn[bpI];
                const scalar d = mag(normal);
                if( d > VSMALL )
                {
                    normal /= d;
                }
                else
                {
                    normal = vector::zero;
                }

                pn[bpI] = normal;
            }
        }
    }
    // Keep polyMeshGenAddressing geometry coherent with point motion.
    //
    // meshSurfaceEngine owns a separate set of geometry caches, updated
    // above.  polyMeshGenAddressing may also already contain cached
    // face centres/areas and cell centres/volumes.  If that demand-driven
    // object exists, update every FULL-MESH face incident to each moved
    // global mesh point.  Do not instantiate addressing solely for this
    // update.
    polyMeshGen& mesh = surfaceEngine_.mesh_;

    if( mesh.hasAddressingData() )
    {
        const polyMeshGenAddressing& addressing =
            mesh.addressingData();

        const VRWGraph& meshPointFaces =
            addressing.pointFaces();

        const labelList& bPoints =
            surfaceEngine_.boundaryPoints();

        boolList changedFace(mesh.faces().size(), false);

        forAll(updateBndNodes, i)
        {
            const label bpI = updateBndNodes[i];

            if( bpI < 0 || bpI >= bPoints.size() )
                continue;

            const label globalPtI = bPoints[bpI];

            if
            (
                globalPtI < 0 ||
                globalPtI >= meshPointFaces.size()
            )
            {
                continue;
            }

            forAllRow(meshPointFaces, globalPtI, pfI)
            {
                const label faceI =
                    meshPointFaces(globalPtI, pfI);

                if
                (
                    faceI >= 0 &&
                    faceI < changedFace.size()
                )
                {
                    changedFace[faceI] = true;
                }
            }
        }

        const_cast<polyMeshGenAddressing&>
        (
            addressing
        ).updateGeometry(changedFace);
    }

}

void meshSurfaceEngineModifier::updateGeometry()
{
    labelLongList updateBndNodes(surfaceEngine_.boundaryPoints().size());

    # ifdef USE_OMP
    # pragma omp parallel for if( updateBndNodes.size() > 10000 )
    # endif
    forAll(updateBndNodes, bpI)
        updateBndNodes[bpI] = bpI;

    updateGeometry(updateBndNodes);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
