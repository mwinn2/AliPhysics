/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliDQTreeHandlerSingleMuons
// \brief helper class to handle a tree for D0 cut optimisation and MVA analyses
// \authors:
// M. Winn, mwinn@cern.ch
// heavily based on AliHFTreeHandler
/////////////////////////////////////////////////////////////

#include <cmath>
#include <limits>
#include <TString.h>
#include "AliDQTreeHandlerSingleMuons.h"

/// \cond CLASSIMP
ClassImp(AliDQTreeHandlerSingleMuons);
/// \endcond

//________________________________________________________________
AliDQTreeHandlerSingleMuons::AliDQTreeHandlerSingleMuons():
  TObject(),
  fTreeVar(nullptr),
  fNCandidates(0),
  fCandType(0),
  fInvMass(-9999.),
  fPt(-9999.),
  fP(-9999.),
  fPtGen(-9999.),
  fY(-9999.),
  fEta(-9999.),
  fPhi(-9999.),
  fCharge(-9999.),
  fDCA(-9999),
  fDCAxy(-9999.),
  fDCAz(-9999.),
  fMuonChi2perNDF(-9999.),
  fRatAbsorberEnd(-9999.),
  fHasMFT(false),
  fMuonMatchChi2perNDF(-9999.),
  fMuonMCHcls(-1),
  fMIDPID(-1),
  fPidOpt(-1),
  fSingleTrackOpt(-1),
  fFillOnlySignal(false),
  fIsMCGenTree(false),
  fTrackInAcceptance(false),
  fabsPDGmother(9999),
  fabsPDGgdmother(9999),
  fEvID(9999),
  fRunNumber(9999),
  fRunNumberPrevCand(9999),
  fMFTChi2perNDF(-999.),
  fMFTcls(-1)
{
  //
  // Default constructor
  //
}

//________________________________________________________________
AliDQTreeHandlerSingleMuons::~AliDQTreeHandlerSingleMuons()
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

  //Set Common variables for muons
  fTreeVar = new TTree(name.Data(),title.Data());
  fTreeVar->Branch("run_number",&fRunNumber);
  fTreeVar->Branch("ev_id",&fEvID);
  fTreeVar->Branch("n_cand",&fNCandidates);
  fTreeVar->Branch("cand_type",&fCandType);
  fTreeVar->Branch("pt_cand",&fPt);
  fTreeVar->Branch("p_cand",&fP);
  fTreeVar->Branch("y_cand",&fY);
  fTreeVar->Branch("eta_cand",&fEta);
  fTreeVar->Branch("phi_cand",&fPhi);
  fTreeVar->Branch("charge_cand",&fCharge);
  fTreeVar->Branch("muon_dca",&fDCA);
  fTreeVar->Branch("muon_dcaxy",&fDCAxy);
  fTreeVar->Branch("muon_dcaz",&fDCAz);
  fTreeVar->Branch("muon_chi2perndf",&fMuonChi2perNDF);
  fTreeVar->Branch("muon_muonmchcls",&fMuonMCHcls);
  fTreeVar->Branch("muon_midpid",&fMIDPID);
  
  //set single track variables specific to muon arm
  AddSingleTrackBranches();

  //set PID variables
  if(fPidOpt!=kNoMID) AddPIDBranches(); 

  return fTreeVar;
}
//_______________________________________________________________
void AliDQTreeHandlerSingleMuons::AddSingleTrackBranches(){

  if(fSingleTrackOpt==kNoSingleTrackVars) return;
  
  fTreeVar->Branch("muon_chi2perndf", &fMuonChi2perNDF);
  fTreeVar->Branch("muon_ratabsorberend", &fRatAbsorberEnd);
  fTreeVar->Branch("muon_fMFTChi2perNDF", &fMuonChi2perNDF);
  fTreeVar->Branch("muon_fMFTcls", &fMFTcls);
  
}
//_______________________________________________________________
void AliDQTreeHandlerSingleMuons::AddPIDBranches(){

  fTreeVar->Branch("muon_midpid", &fMIDPID);
  fTreeVar->Branch("muon_midchi2", &fMuonMatchChi2perNDF);
  fTreeVar->Branch("muon_hasmft", &fHasMFT);
  
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
  fTreeVar->Branch("charge_cand",fCharge);
  fTreeVar->Branch("track_in_acc",&fTrackInAcceptance);
  fTreeVar->Branch("absPDGmother",&fabsPDGmother);
  fTreeVar->Branch("absPDGgdmother",&fabsPDGmother);
  return fTreeVar;

}
//________________________________________________________________
bool AliDQTreeHandlerSingleMuons::SetVariables(int runnumber, unsigned int eventID, float ptgen, AliAODTrack* cand, float bfield) 
{
  fIsMCGenTree=false;

  if(!cand) return false;
  if(fFillOnlySignal) { //if fill only signal and not signal candidate, do not store
    if(!(fCandType&kSignal)) return true;
  }
  fNCandidates++;
  fRunNumber=runnumber;
  fEvID=eventID;
  fPtGen=ptgen;
  
  //topological variables
  //common
  fPt=cand->Pt();
  fP=cand->P();
  fY=cand->Y(1);//TBC: 1 should be muon
  fEta=cand->Eta();
  fPhi=cand->Phi();
  fCharge=cand->Charge();
  fDCA=cand->DCA();
  fDCAxy=TMath::Sqrt(cand->XAtDCA()*cand->XAtDCA()+cand->YAtDCA()*cand->YAtDCA());//in AOD the same as DCA!
  fDCAz=cand->ZAtDCA();
  //MFT related variables to be filled
  
  fInvMass= 0; //to be done with two prongs

  SetSingleTrackVars(cand);
  bool setpid = SetPidVars(cand);
  if(!setpid) return false;

  return true;
}
//________________________________________________________________
bool AliDQTreeHandlerSingleMuons::SetMCGenVariables(int runnumber, unsigned int eventID, AliAODMCParticle* mcpart, unsigned int pdgmother, unsigned int pdggdmother ){
  

  if(!mcpart) return false;
  if(!(fCandType&kSignal)) return true;

  fRunNumber = runnumber;
  fEvID = eventID;
  fPt = mcpart->Pt();
  fY = mcpart->Y();
  fEta = mcpart->Eta();
  fPhi = mcpart->Phi();
  fabsPDGmother = pdgmother;
  fabsPDGgdmother = pdggdmother;

  return true;
}
//________________________________________________________________
void AliDQTreeHandlerSingleMuons::SetCandidateType(bool issignal, bool isbkg, bool ischarm, bool isbeauty, bool isew, bool isquarkonia)
{ 
  if(issignal) fCandType |= kSignal;
  else fCandType &= ~kSignal;
  if(isbkg && !fIsMCGenTree) fCandType |= kBkg;
  else fCandType &= ~kBkg;
  if(ischarm) fCandType |= kCharm;
  else fCandType &= ~kCharm;
  if(isbeauty) fCandType |= kBeauty;
  else fCandType &= ~kBeauty;
  if(isew) fCandType |= kEW;
  else fCandType &= ~kEW;
  if(isquarkonia && !fIsMCGenTree) fCandType |= kQuarkonia;
  else fCandType &= ~kQuarkonia;
}
//________________________________________________________________
bool AliDQTreeHandlerSingleMuons::SetSingleTrackVars(AliAODTrack* muon){
  //muon arm specific single track variables

  if(fSingleTrackOpt==kNoSingleTrackVars) return true;

  if(fSingleTrackOpt==kSingleTrackVars) {
    //chi2perNDF
    fMuonChi2perNDF = muon->Chi2perNDF();//check if correct for muons
    //R at entrance of absorber
    fRatAbsorberEnd = muon->GetRAtAbsorberEnd();
    //add number of muon clusters TODO
    //add number of MFT clusters TODO
  }


  return true;
}
//________________________________________________________________
bool AliDQTreeHandlerSingleMuons::SetPidVars(AliAODTrack* track)
{
  fMIDPID = track->GetMatchTrigger();
  fMuonMatchChi2perNDF= (float) track->GetChi2MatchTrigger();
  fHasMFT = track->IsMuonGlobalTrack();
  return kTRUE;
  
}

