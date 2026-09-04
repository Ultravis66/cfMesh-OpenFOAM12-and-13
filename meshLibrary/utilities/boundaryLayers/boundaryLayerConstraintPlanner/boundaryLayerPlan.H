/*---------------------------------------------------------------------------*\
  CFMitch - resolved boundary-layer plan
\*---------------------------------------------------------------------------*/

#ifndef boundaryLayerPlan_H
#define boundaryLayerPlan_H

#include "labelList.H"

namespace Foam
{

class boundaryLayerPlan
{
    //- Final number of layers assigned to each current boundary face
    labelList faceLayers_;


public:

    explicit boundaryLayerPlan
    (
        const labelList& requestedFaceLayers
    )
    :
        faceLayers_(requestedFaceLayers)
    {}


    labelList& faceLayers()
    {
        return faceLayers_;
    }


    const labelList& faceLayers() const
    {
        return faceLayers_;
    }
};

} // End namespace Foam

#endif
