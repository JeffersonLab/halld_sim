
#include "RegisterGlueXModels.h"
#include <iostream>

#include "EvtGenBase/EvtModel.hh"

#include "EvtEtaDalitz_GlueX.h"
#include "EvtEtaPiPiGamma.h"

namespace GlueX_EvtGen 
{

  void RegisterGlueXModels() 
  {
    EvtModel &modelist = EvtModel::instance();

    // All new models must be added here, or EvtGen won't know them
    modelist.registerModel(new EvtEtaDalitz_GlueX);
    modelist.registerModel(new EvtEtaPiPiGamma);
  }


}
