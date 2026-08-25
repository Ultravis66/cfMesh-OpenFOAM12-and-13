/*---------------------------------------------------------------------------*\
Description
    BLJunctionClassifier: classifies post-refBL bad pyramid faces by
    junction type using PatchRoleMap, then generates a BLRepairPlan.

    Classification is neighborhood-based: collects roles from the source
    bfI patch AND its immediate boundary face neighbors via faceFaces().
    This correctly identifies blade-periodic junctions even when the
    source bfI belongs to a blade patch.

    Pipeline:
      1. bad face -> owner/neighbour cell -> source bfI (provenance)
      2. source bfI -> source patch roles
      3. faceFaces(bfI) -> neighbor patch roles
      4. combined role set -> junction type
      5. group source bfI by type and connected component
         -> one BLRepairPlan per connected defect component

    First integration: diagnostic only. Log classification counts.
    Do not apply repair until classification is verified correct.
\*---------------------------------------------------------------------------*/

#ifndef BLJunctionClassifier_H
#define BLJunctionClassifier_H

#include "HashSet.H"

#include "BLRepairPlan.H"
#include "PatchRoleMap.H"
#include "polyMeshGen.H"
#include "meshSurfaceEngine.H"
#include "Map.H"
#include "PtrList.H"
#include "boundaryPatch.H"

namespace Foam
{

class BLJunctionClassifier
{
    polyMeshGen&  mesh_;
    const PatchRoleMap& roles_;

    //- bfI -> patch index
    labelList bfIToPatch_;

    void buildBfIToPatch()
    {
        const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
        const label nInternalFaces = mesh_.nInternalFaces();

        label nBndFaces = 0;
        forAll(boundaries, pI)
            nBndFaces += boundaries[pI].patchSize();

        bfIToPatch_.setSize(nBndFaces, -1);

        forAll(boundaries, pI)
        {
            const label start = boundaries[pI].patchStart();
            const label size  = boundaries[pI].patchSize();

            for(label pfI=0; pfI<size; ++pfI)
            {
                const label globalFaceI = start + pfI;
                const label bfI = globalFaceI - nInternalFaces;

                if( bfI >= 0 && bfI < label(bfIToPatch_.size()) )
                    bfIToPatch_[bfI] = pI;
            }
        }
    }

    //- Collect all roles around a source bfI from source patch
    //- and immediate face neighbors via faceFaces()
    HashSet<word> neighborhoodRoles
    (
        const label bfI,
        const VRWGraph& faceFaces,
        const PtrList<boundaryPatch>& boundaries
    ) const
    {
        HashSet<word> localRoles;

        auto addPatchRoles = [&](const label bfJ)
        {
            if( bfJ < 0 || bfJ >= label(bfIToPatch_.size()) ) return;
            const label pI = bfIToPatch_[bfJ];
            if( pI < 0 || pI >= label(boundaries.size()) ) return;
            const word& pName = boundaries[pI].patchName();
            HashSet<word> pRoles = roles_.rolesOf(pName);
            forAllConstIter(HashSet<word>, pRoles, rIt)
                localRoles.insert(rIt.key());
        };

        // Source face
        addPatchRoles(bfI);

        // Immediate face neighbors
        if( bfI >= 0 && bfI < label(faceFaces.size()) )
            forAllRow(faceFaces, bfI, nI)
                addPatchRoles(faceFaces(bfI, nI));

        return localRoles;
    }

    //- Classify combined role set to junction type
    BLRepairPlan::junctionType classifyRoles
    (
        const HashSet<word>& localRoles
    ) const
    {
        if( localRoles.found("wall") && localRoles.found("periodic") )
            return BLRepairPlan::WALL_PERIODIC;

        if( localRoles.found("hub") && localRoles.found("inletOutlet") )
            return BLRepairPlan::HUB_FLOW;

        if( localRoles.found("shroud") && localRoles.found("inletOutlet") )
            return BLRepairPlan::SHROUD_FLOW;

        if( localRoles.found("wall") )
            return BLRepairPlan::WALL_WALL;

        return BLRepairPlan::GENERIC_TRIPLE;
    }

    //- Default ring params by junction type
    void defaultRingParams
    (
        const BLRepairPlan::junctionType jType,
        BLRepairPlan::ringParams& r0,
        BLRepairPlan::ringParams& r1,
        BLRepairPlan::ringParams& r2
    ) const
    {
        switch(jType)
        {
            case BLRepairPlan::WALL_PERIODIC:
                // Most aggressive: blade/periodic corners
                r0 = BLRepairPlan::ringParams(2, 0.5);
                r1 = BLRepairPlan::ringParams(3, 0.75);
                r2 = BLRepairPlan::ringParams(4, 1.0);
                break;

            case BLRepairPlan::HUB_FLOW:
            case BLRepairPlan::SHROUD_FLOW:
                // Moderate: annular inlet/outlet rings
                r0 = BLRepairPlan::ringParams(2, 0.6);
                r1 = BLRepairPlan::ringParams(4, 0.8);
                r2 = BLRepairPlan::ringParams(0, 1.0);
                break;

            case BLRepairPlan::WALL_WALL:
                // Conservative: blade/hub or blade/shroud
                r0 = BLRepairPlan::ringParams(2, 1.0);
                r1 = BLRepairPlan::ringParams(4, 1.0);
                r2 = BLRepairPlan::ringParams(0, 1.0);
                break;

            default:
                // Generic-triple disabled for mechanism isolation.
                // Preserve the normal 3-layer BL and apply no thickness scaling.
                r0 = BLRepairPlan::ringParams(3, 1.0);
                r1 = BLRepairPlan::ringParams(0, 1.0);
                r2 = BLRepairPlan::ringParams(0, 1.0);
                break;
        }
    }

public:

    BLJunctionClassifier
    (
        polyMeshGen& mesh,
        const PatchRoleMap& roles
    )
    :
        mesh_(mesh),
        roles_(roles)
    {
        buildBfIToPatch();
    }

    //- Classify bad faces and generate repair plans by junction type.
    //- diagnosticOnly=true: log only, do not populate plan seed sets.
    Map<BLRepairPlan> classify
    (
        const HashSet<label>& badFaces,
        const labelList& cellToBaseBndFace,
        const bool diagnosticOnly = true
    ) const
    {
        const labelList& own = mesh_.owner();
        const labelList& nei = mesh_.neighbour();
        const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

        // Build surface engine for faceFaces connectivity
        const meshSurfaceEngine mse(mesh_);
        const VRWGraph& faceFaces = mse.faceFaces();

        // Count by junction type
        Map<label> typeCounts;
        Map<HashSet<label>> typeToSeeds;

        label nWithProv = 0;
        label nNoProv = 0;

        forAllConstIter(HashSet<label>, badFaces, it)
        {
            const label faceI = it.key();

            for( label side = 0; side < 2; ++side )
            {
                label cellI = -1;
                if( side == 0 )
                {
                    if( faceI < 0 || faceI >= label(own.size()) ) continue;
                    cellI = own[faceI];
                }
                else
                {
                    if( faceI < 0 || faceI >= label(nei.size()) ) continue;
                    cellI = nei[faceI];
                }

                if( cellI < 0 || cellI >= label(cellToBaseBndFace.size()) )
                    continue;

                const label bfI = cellToBaseBndFace[cellI];

                if( bfI < 0 )
                {
                    ++nNoProv;
                    continue;
                }

                ++nWithProv;

                // Neighborhood-based classification
                const HashSet<word> localRoles =
                    neighborhoodRoles(bfI, faceFaces, boundaries);

                const BLRepairPlan::junctionType jType =
                    classifyRoles(localRoles);

                const label jKey = label(jType);

                if( !typeCounts.found(jKey) )
                    typeCounts.insert(jKey, 0);
                typeCounts[jKey]++;

                if( !diagnosticOnly )
                {
                    if( !typeToSeeds.found(jKey) )
                        typeToSeeds.insert(jKey, HashSet<label>());
                    typeToSeeds[jKey].insert(bfI);
                }
            }
        }

        // Log classification summary
        Info << "BLJunctionClassifier: badFaces=" << badFaces.size()
             << " withProv=" << nWithProv
             << " noProv=" << nNoProv
             << endl;

        forAllConstIter(Map<label>, typeCounts, it)
        {
            const BLRepairPlan::junctionType jType =
                BLRepairPlan::junctionType(it.key());
            Info << "  " << BLRepairPlan::junctionTypeName(jType)
                 << ": " << it() << " faces" << endl;
        }

        // Build repair plans (empty if diagnosticOnly).
        //
        // Each junction type may contain multiple spatially disconnected
        // defect regions.  Do not collapse those regions into one repair
        // transaction: split the provenance seed boundary faces into
        // connected components using boundary-face adjacency.
        //
        // The returned Map key is only a unique plan/component ID.  Callers
        // consume the BLRepairPlan value and do not attach semantics to the
        // key.
        Map<BLRepairPlan> plans;

        if( !diagnosticOnly )
        {
            label nextPlanId = 0;

            forAllIter(Map<HashSet<label>>, typeToSeeds, it)
            {
                const BLRepairPlan::junctionType jType =
                    BLRepairPlan::junctionType(it.key());

                const HashSet<label>& seeds = it();

                BLRepairPlan::ringParams r0, r1, r2;
                defaultRingParams(jType, r0, r1, r2);

                // Flood only through seed faces of this junction type.
                // This deliberately defines connectivity on the defect seeds
                // themselves; repair-footprint overlap is a separate policy
                // decision and is not introduced in this decomposition step.
                HashSet<label> visitedSeeds;
                label nComponents = 0;

                forAllConstIter(HashSet<label>, seeds, seedIt)
                {
                    const label startBfI = seedIt.key();

                    if( visitedSeeds.found(startBfI) )
                        continue;

                    HashSet<label> component;

                    // A component can contain at most all seeds of this type.
                    // Use a simple bounded FIFO rather than depending on
                    // DynamicList stack semantics.
                    labelList queue(seeds.size(), -1);
                    label qHead = 0;
                    label qTail = 0;

                    visitedSeeds.insert(startBfI);
                    component.insert(startBfI);
                    queue[qTail++] = startBfI;

                    while( qHead < qTail )
                    {
                        const label bfI = queue[qHead++];

                        if( bfI < 0 || bfI >= label(faceFaces.size()) )
                            continue;

                        forAllRow(faceFaces, bfI, nI)
                        {
                            const label nbBfI = faceFaces(bfI, nI);

                            if
                            (
                                seeds.found(nbBfI)
                             && !visitedSeeds.found(nbBfI)
                            )
                            {
                                visitedSeeds.insert(nbBfI);
                                component.insert(nbBfI);

                                // qTail cannot exceed seeds.size() because
                                // every queued face is inserted into
                                // visitedSeeds exactly once.
                                queue[qTail++] = nbBfI;
                            }
                        }
                    }

                    const label componentI = nComponents++;
                    const label planId = nextPlanId++;

                    Info << "BLRepairComponent:"
                         << " planId=" << planId
                         << " type="
                         << BLRepairPlan::junctionTypeName(jType)
                         << " component=" << componentI
                         << " seeds=" << component.size()
                         << endl;

                    BLRepairPlan plan
                    (
                        component,
                        r0, r1, r2,
                        jType,
                        BLRepairPlan::junctionTypeName(jType)
                    );

                    plan.print();
                    plans.insert(planId, plan);
                }

                Info << "BLRepairComponents:"
                     << " type="
                     << BLRepairPlan::junctionTypeName(jType)
                     << " components=" << nComponents
                     << " seeds=" << seeds.size()
                     << " visited=" << visitedSeeds.size()
                     << endl;
            }
        }

        return plans;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
