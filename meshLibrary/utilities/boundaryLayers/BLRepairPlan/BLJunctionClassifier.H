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
      5. group source bfI by type -> one BLRepairPlan per type

    First integration: diagnostic only. Log classification counts.
    Do not apply repair until classification is verified correct.
\*---------------------------------------------------------------------------*/

#ifndef BLJunctionClassifier_H
#define BLJunctionClassifier_H

#include "BLRepairPlan.H"
#include "PatchRoleMap.H"
#include "polyMeshGen.H"
#include "meshSurfaceEngine.H"
#include "labelHashSet.H"
#include "wordHashSet.H"
#include "Map.H"
#include "PtrList.H"
#include "boundaryPatch.H"

namespace Foam
{

class BLJunctionClassifier
{
    const polyMeshGen&  mesh_;
    const PatchRoleMap& roles_;

    //- bfI -> patch index
    labelList bfIToPatch_;

    void buildBfIToPatch()
    {
        const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
        label nBndFaces = 0;
        forAll(boundaries, pI)
            nBndFaces += boundaries[pI].patchSize();

        bfIToPatch_.setSize(nBndFaces, -1);
        label bfI = 0;
        forAll(boundaries, pI)
            for(label f=0; f<boundaries[pI].patchSize(); ++f)
                bfIToPatch_[bfI++] = pI;
    }

    //- Collect all roles around a source bfI from source patch
    //- and immediate face neighbors via faceFaces()
    wordHashSet neighborhoodRoles
    (
        const label bfI,
        const VRWGraph& faceFaces,
        const PtrList<boundaryPatch>& boundaries
    ) const
    {
        wordHashSet localRoles;

        auto addPatchRoles = [&](const label bfJ)
        {
            if( bfJ < 0 || bfJ >= label(bfIToPatch_.size()) ) return;
            const label pI = bfIToPatch_[bfJ];
            if( pI < 0 || pI >= label(boundaries.size()) ) return;
            const word& pName = boundaries[pI].patchName();
            wordHashSet pRoles = roles_.rolesOf(pName);
            forAllConstIter(wordHashSet, pRoles, rIt)
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
        const wordHashSet& localRoles
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
                r0 = BLRepairPlan::ringParams(3, 0.6);
                r1 = BLRepairPlan::ringParams(4, 0.8);
                r2 = BLRepairPlan::ringParams(0, 1.0);
                break;

            case BLRepairPlan::WALL_WALL:
                // Conservative: blade/hub or blade/shroud
                r0 = BLRepairPlan::ringParams(3, 0.7);
                r1 = BLRepairPlan::ringParams(4, 0.9);
                r2 = BLRepairPlan::ringParams(0, 1.0);
                break;

            default:
                // Generic: minimal intervention
                r0 = BLRepairPlan::ringParams(3, 0.8);
                r1 = BLRepairPlan::ringParams(0, 1.0);
                r2 = BLRepairPlan::ringParams(0, 1.0);
                break;
        }
    }

public:

    BLJunctionClassifier
    (
        const polyMeshGen& mesh,
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
        const labelHashSet& badFaces,
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
        Map<labelHashSet> typeToSeeds;

        label nWithProv = 0;
        label nNoProv = 0;

        forAllConstIter(labelHashSet, badFaces, it)
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
                const wordHashSet localRoles =
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
                        typeToSeeds.insert(jKey, labelHashSet());
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

        // Build repair plans (empty if diagnosticOnly)
        Map<BLRepairPlan> plans;

        if( !diagnosticOnly )
        {
            forAllIter(Map<labelHashSet>, typeToSeeds, it)
            {
                const BLRepairPlan::junctionType jType =
                    BLRepairPlan::junctionType(it.key());

                BLRepairPlan::ringParams r0, r1, r2;
                defaultRingParams(jType, r0, r1, r2);

                BLRepairPlan plan
                (
                    it(),
                    r0, r1, r2,
                    jType,
                    BLRepairPlan::junctionTypeName(jType)
                );

                plan.print();
                plans.insert(it.key(), plan);
            }
        }

        return plans;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
