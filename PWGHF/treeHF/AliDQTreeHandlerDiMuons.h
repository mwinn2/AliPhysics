#ifndef ALIHFTREEHANDLERD0TOKPI_H
#define ALIHFTREEHANDLERD0TOKPI_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerDiMuons
// \brief helper class to handle a tree for Dimuon cut optimisation and MVA analyses
// \authors:
// M. Winn, mwinn@cern.ch
//heavily based on AliHFTreeHandlerD0toKpi
/////////////////////////////////////////////////////////////

#include "AliHFTreeHandler.h"

class AliHFTreeHandlerDiMuons : public AliHFTreeHandler
{
  public:
    
    enum isDzeroDzeroBar {
      kDzeroTopo         = BIT(11),//Bit meaning to be understood
        kDzeroBarTopo      = BIT(12),
        kDzeroPID          = BIT(13),
        kDzeroBarPID       = BIT(14),
        kDzeroComb         = BIT(15),
        kDzeroBarComb      = BIT(16),
        kDzeroTopoFilt     = BIT(17),
        kDzeroBarTopoFilt  = BIT(18),
        kDzeroPIDFilt      = BIT(19),
        kDzeroBarPIDFilt   = BIT(20),
        kDzeroCombFilt     = BIT(21),
        kDzeroBarCombFilt  = BIT(22)
    };
    
    AliHFTreeHandlerD0toKpi();
    AliHFTreeHandlerD0toKpi(int PIDopt);

    virtual ~AliHFTreeHandlerD0toKpi();

    virtual TTree* BuildTree(TString name="tree", TString title="tree");
    virtual bool SetVariables(int runnumber, unsigned int eventID, float ptgen, AliAODRecoDecayHF* cand, float bfiled, int masshypo=0, AliPIDResponse *pidrespo=nullptr);
    
    void SetIsDzeroDzeroBar(int isSel, int isSelTopo, int isSelPID, int isSelFilt, int isSelTopoFilt, int isSelPIDFilt);

  private:

    float fImpParProng[knMaxProngs]; ///prong impact parameter
    float fCosThetaStar; /// candidate costhetastar
    float fImpParProd; /// daughter impact-parameter product
    float fNormd0MeasMinusExp; ///candidate topomatic variable
    float fImpParErrProng[knMaxProngs]; ///error on prongs rphi impact param [cm]
    float fNormDecayLength; ///candidate normalised decay length

    /// \cond CLASSIMP
    ClassDef(AliHFTreeHandlerD0toKpi,3); /// 
    /// \endcond
};
#endif
