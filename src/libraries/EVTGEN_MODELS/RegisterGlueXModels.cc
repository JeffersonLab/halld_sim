
#include "RegisterGlueXModels.h"
#include <iostream>

#include "EvtGenBase/EvtModel.hh"

#include "EvtEtaDalitz_GlueX.h"
#include "EvtEtaDalitz_3Pi0_GlueX.h"
#include "EvtEtaPiPiGamma.h"
#include "EvtEtaPiPiGamma_GlueX.h"

namespace GlueX_EvtGen 
{

  void RegisterGlueXModels() 
  {
    EvtModel &modelist = EvtModel::instance();

    // All new models must be added here, or EvtGen won't know them
    modelist.registerModel(new EvtEtaDalitz_GlueX);
    modelist.registerModel(new EvtEtaDalitz_3Pi0_GlueX);
    modelist.registerModel(new EvtEtaPiPiGamma);
  }


}
