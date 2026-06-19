/*---------------------------------------------------------------------------*\
Description
    PatchRoleMap: maps patch names to semantic roles for junction
    classification. Roles are user-defined in meshDict patchRoles block.
    C++ code uses roles only -- never hardcoded patch names.

    A patch may have multiple roles.

    meshDict example:
        boundaryLayers
        {
            patchRoles
            {
                wall        (blade_1 blade_2 blade_3 blade_4 hub shroud);
                blade       (blade_1 blade_2 blade_3 blade_4);
                hub         (hub);
                shroud      (shroud);
                periodic    (periodic_1 periodic_2);
                inletOutlet (inlet outlet);
            }
        }
\*---------------------------------------------------------------------------*/

#ifndef PatchRoleMap_H
#define PatchRoleMap_H

#include "dictionary.H"
#include "wordList.H"
#include "HashTable.H"
#include "wordHashSet.H"

namespace Foam
{

class PatchRoleMap
{
    //- role name -> set of patch names
    HashTable<wordHashSet> roleToPatches_;

    //- patch name -> set of role names (multi-role support)
    HashTable<wordHashSet> patchToRoles_;

public:

    PatchRoleMap()
    {}

    //- Construct from meshDict boundaryLayers.patchRoles block
    PatchRoleMap(const dictionary& meshDict)
    {
        if( !meshDict.isDict("boundaryLayers") )
            return;

        const dictionary& bndL = meshDict.subDict("boundaryLayers");

        if( !bndL.isDict("patchRoles") )
            return;

        const dictionary& roles = bndL.subDict("patchRoles");
        const wordList roleNames = roles.toc();

        forAll(roleNames, rI)
        {
            const word& role = roleNames[rI];
            const wordList patches(roles.lookup(role));

            if( !roleToPatches_.found(role) )
                roleToPatches_.insert(role, wordHashSet());

            forAll(patches, pI)
            {
                const word& patchName = patches[pI];

                roleToPatches_[role].insert(patchName);

                if( !patchToRoles_.found(patchName) )
                    patchToRoles_.insert(patchName, wordHashSet());

                patchToRoles_[patchName].insert(role);
            }
        }

        Info << "PatchRoleMap: loaded "
             << roleToPatches_.size() << " roles, "
             << patchToRoles_.size() << " patches" << endl;
    }

    //- Return true if any roles are mapped
    bool active() const
    {
        return patchToRoles_.size() > 0;
    }

    //- Return true if patch has the given role
    bool hasRole(const word& patchName, const word& role) const
    {
        if( !patchToRoles_.found(patchName) )
            return false;

        return patchToRoles_[patchName].found(role);
    }

    //- Backwards-friendly alias
    bool isRole(const word& patchName, const word& role) const
    {
        return hasRole(patchName, role);
    }

    //- Return all roles for a patch, empty if unknown
    wordHashSet rolesOf(const word& patchName) const
    {
        if( patchToRoles_.found(patchName) )
            return patchToRoles_[patchName];

        return wordHashSet();
    }

    //- Return true if two patches form the requested role pair,
    //- independent of order.
    bool hasRolePair
    (
        const word& patchA,
        const word& patchB,
        const word& roleA,
        const word& roleB
    ) const
    {
        return
            (hasRole(patchA, roleA) && hasRole(patchB, roleB))
         || (hasRole(patchA, roleB) && hasRole(patchB, roleA));
    }

    //- Print summary
    void print() const
    {
        Info << "PatchRoleMap:" << endl;

        forAllConstIter(HashTable<wordHashSet>, roleToPatches_, it)
        {
            Info << "  role=" << it.key() << " patches=(";

            forAllConstIter(wordHashSet, it(), pIt)
                Info << " " << pIt.key();

            Info << " )" << endl;
        }
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
