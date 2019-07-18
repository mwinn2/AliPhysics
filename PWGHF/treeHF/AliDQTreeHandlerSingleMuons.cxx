/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliDQTreeHandlerDiMuons
// \brief helper class to handle a tree for D0 cut optimisation and MVA analyses
// \authors:
// M. Winn, mwinn@cern.ch
// heavily based on AliHFTreeHandler
/////////////////////////////////////////////////////////////

#include <cmath>
#include <limits>
#include <TString.h>
#include "AliDQTreeHandlerSingleMuons.h"
//#include "AliAODRecoDecayHF2Prong.h"//TBC

/// \cond CLASSIMP
ClassImp(AliDQTreeHandlerSingleMuons);
/// \endcond

//________________________________________________________________
AliDQTreeHandlerSingleMuons::AliDQTreeHandlerSingleMuons():
  TObject(),
  fTreeVar(nullptr),
  fNProngs(-1),
  fCandType(0),
  fInvMass(-9999.),
  fPt(-9999.),
  fPtGen(-9999.),
  fY(-9999.),
  fEta(-9999.),
  fPhi(-9999.),
  fDecayLength(-9999.),
  fDecayLengthNorm(-9999.),
  fDecayLengthZ(-9999.),
  fNormDecayLengthZ(-9999.),
  fMuonDCA(-9999.),
  fMuonP(-9999.),
  fMuonPt(-9999.),
  fMuonEta(-9999.),
  fMuonPhi(-9999.),
  fMuonChi2perNDF(-9999.),
  fMUONMCHcls(-1),
  fMIDPID(-9999.),
  fPidOpt(-1),
  fSingleTrackOpt(-1),
  fFillOnlySignal(false),
  fIsMCGenTree(false),
  fTrackInAcceptance(false),
  fEvID(9999),
  fRunNumber(9999),
  fRunNumberPrevCand(9999)
{
  //
  // Default constructor
  //

  fNProngs=1; // --> cannot be changed
  //keep for the moment only the single muon
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
    fMFTChi2perNDFProng[iProng] = -9999.;
    fMFTclsProng[iProng] = -9999.;
    fImpParProng[iProng] = -9999.;
    fImpParErrProng[iProng] = -9999.;
  }
}

//________________________________________________________________
AliDQTreeHandlerSingleMuons::~AliHFTreeHandlerSingleMuons()
{
  //
  // Default Destructor
  //
  
  if(fTreeVar) delete fTreeVar;
}

//________________________________________________________________
TTree* AliDQTreeHandlerSingleMuons::BuildTree(TString name, TString title) 
{
  fIsMCGenTree=false;
  
  if(fTreeVar) {
    delete fTreeVar;
    fTreeVar=nullptr;
  }

  //Set Common variables
  fTreeVar = new TTree(name.Data(),title.Data());
  fTreeVar->Branch("run_number",&fRunNumber);
  fTreeVar->Branch("ev_id",&fEvID);
  fTreeVar->Branch("cand_type",&fCandType);
  fTreeVar->Branch("pt_cand",&fPt);
  fTreeVar->Branch("y_cand",&fY);
  fTreeVar->Branch("eta_cand",&fEta);
  fTreeVar->Branch("phi_cand",&fPhi);
  fTreeVar->Branch("d_len",&fDecayLength);
  fTreeVar->Branch("norm_dl",&fDecayLengthNorm);
  fTreeVar->Branch("d_len_z",&fDecayLengthZ);
  fTreeVar->Branch("norm_dl_z",&fNormDecayLengthZ);
  fTreeVar->Branch("muon_dca",&fMuonDCA);
  fTreeVar->Branch("muon_p",&fMuonP);
  fTreeVar->Branch("muon_pt",&fMuonPt);
  fTreeVar->Branch("muon_eta",&fMuonEta);
  fTreeVar->Branch("muon_phi",&fMuonPhi);
  fTreeVar->Branch("muon_chi2perndf",&fMuonChi2perNDF);
  fTreeVar->Branch("muon_muonmchcls",&fMuonMCHcls);
  fTreeVar->Branch("muon_midpid",&fMIDPID);
  
  //set common variables TODO
  //  AddSingleTrackBranches();


  //set PID variables
  //to be replaced
  //  if(fPidOpt!=kNoPID) AddPidBranches(true,true,false,true,true); 

  return fTreeVar;
}
//________________________________________________________________
TTree* AliDQTreeHandlerSingleMuons::BuildTreeMCGen(TString name, TString title) {

  fIsMCGenTree = true;

  if(fTreeVar) {
    delete fTreeVar;
    fTreeVar=nullptr;
  }
  fTreeVar = new TTree(name.Data(),title.Data());
  fTreeVar->Branch("run_number",&fRunNumber);
  fTreeVar->Branch("ev_id",&fEvID);
  fTreeVar->Branch("cand_type",&fCandType);
  fTreeVar->Branch("pt_cand",&fPt);
  fTreeVar->Branch("y_cand",&fY);
  fTreeVar->Branch("eta_cand",&fEta);
  fTreeVar->Branch("phi_cand",&fPhi);
  fTreeVar->Branch("track_in_acc",&fInAcceptance);

  return fTreeVar;

}
//________________________________________________________________
bool AliHFTreeHandlerDiMuons::SetVariables(int runnumber, unsigned int eventID, float ptgen, AliAODTrack* cand, float bfield, int masshypo) 
//needs to be changed to standard Dimuons, is there already an object like that?
{
  fIsMCGenTree=false;

  if(!cand) return false;
  if(fFillOnlySignal) { //if fill only signal and not signal candidate, do not store
    if(!fCandType&kSignal) return true;
  }
  fNCandidates++;
  fRunNumber=runnumber;
  fEvID=eventID;
  fPtGen=ptgen;
  
  //topological variables
  //common
  fPt=cand->Pt();// TODO all to be changed from a dimuon
  fY=cand->Y(1);//TBC: 1 should be muon
  fEta=cand->Eta();
  fPhi=cand->Phi();
  fDecayLength= 0.0; //to be seen how to do best to get it from AOD in principle, just vector length between PV and seconary vertex
  fDecayLengthXY= 0.0; //cand->DecayLengthXY(); idem
  fNormDecayLengthXY= 0.0; //cand->NormalizedDecayLengthXY();
  fDCA=cand->DCA();
  fNormDecayLength=cand->NormalizedDecayLength();
  
  //D0 -> Kpi variables

  fInvMass= 0;//to be done with two prongs

  //single prong variables TBD
  //PID to be done
  //  bool setpid = SetPidVars(prongtracks,pidrespo,true,true,false,true,true);
  // if(!setpid) return false;

  return true;
}

