/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliDQTreeHandlerDiMuons
// \brief helper class to handle a tree for Dimuons (quarkonia, low-masses, tobe promoted for 3 muons)
// \authors:
// M. Winn, mwinn@cern.ch
// heavily based on AliHFTreeHandlerDKpi
/////////////////////////////////////////////////////////////

#include <cmath>
#include <limits>
#include <TString.h>
#include "AliDQTreeHandlerDiMuons.h"

/// \cond CLASSIMP
ClassImp(AliDQTreeHandlerDiMuons);
/// \endcond

//________________________________________________________________
AliDQTreeHandlerDiMuons::AliDQTreeHandlerDiMuons():
  TObject(),
  fTreeVar(nullptr),
  fNProngs(-1),
  fNCandidates(0),
  fCandType(0),
  fInvMass(-9999.),
  fPt(-9999.),
  fP(-9999.),
  fPtGen(-9999.),
  fY(-9999.),
  fEta(-9999.),
  fPhi(-9999.),
  fDecayLength(-9999.),
  fNormDecayLength(-9999.),
  fDecayLengthZ(-9999.),
  fNormDecayLengthZ(-9999.),
  fPseudopropertimeZ(-9999.),
  fPseudopropertimeZres(-9999.),
  fPseudopropertimeXY(-9999.),
  fPseudopropertimeXYres(-9999.),
  fCosP(-9999.),
  fCosPz(-9999.),
  fImpParZ(-9999.),
  fXF(-9999.),
  fCostCS(-9999.),
  fCostHe(-9999.),
  fPhiCS(-9999.),
  fPhiHe(-9999.),
  fPidOpt(-1),
  fSingleTrackOpt(true),
  fFillOnlySignal(false),
  fIsMCGenTree(false),
  fDauInAcceptance(true),
  fEvID(0),
  fRunNumber(-1),
  fRunNumberPrevCand(-1),
  fCosThetaStar(-9999.),
  fImpParProd(-9999.),
  fNormjpsiMeasMinusExp(-9999.)
{
  //
  // Default constructor
  //
  fNProngs=2; // --> cannot be changed
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
    fCharge[iProng] = -999.;
    fDCA[iProng] = -9999.;
    fDCAz[iProng] = -999.;
    fPProng[iProng] = 9999.;
    fMCHPProng[iProng] = 9999.;
    fPtProng[iProng] = 9999.;
    fMCHPtProng[iProng] = 9999.;
    fEtaProng[iProng] = 9999.;
    fMCHEtaProng[iProng] = 9999.;
    fPhiProng[iProng] = 9999.;
    fMCHPhiProng[iProng] = 9999.;
    fNMFTclsProng[iProng] = -1;
    fMFTclsMapProng[iProng] = -1;
    fNMCHclsProng[iProng] = -1;
    fMIDPID[iProng] = true;
    fMIDchi2perNDF[iProng] = -9999.;
    fImpParProng[iProng] = -9999.;
    fImpParErrProng[iProng] = -9999.;
    fMFTChi2perNDF[iProng] = -9999.;
    fMuonChi2perNDF[iProng] = -999.;
    fMuonMatchChi2perNDF[iProng] = -999.;
    fRatAbsorberEnd[iProng] = -999.;
    fHasMFTProng[iProng] = true;
    fabsPDGProng[iProng] = 999;
    fabsPDGProngmother[iProng] = 999;
    fabsPDGPronggdmother[iProng] = 999;
  }
}
//________________________________________________________________
AliDQTreeHandlerDiMuons::~AliDQTreeHandlerDiMuons()
{
  //
  // Default Destructor
  //
}

//________________________________________________________________
TTree* AliDQTreeHandlerDiMuons::BuildTree(TString name, TString title) 
{
  fIsMCGenTree=false;
  
  if(fTreeVar) {
    delete fTreeVar;
    fTreeVar=nullptr;
  }
  fTreeVar = new TTree(name.Data(),title.Data());

  //set common variables
  // need to be replaced, simple things momentum eta etc. pp
  //  AddCommonDmesonVarBranches();
  //CONTINUE HERE!
  //set Dimuon variables
  fTreeVar->Branch("cos_t_star",&fCosThetaStar);
  fTreeVar->Branch("imp_par_prod",&fImpParProd);
  fTreeVar->Branch("max_norm_Jpsiexp",&fNormjpsiMeasMinusExp);//TODO change naming
  fTreeVar->Branch("norm_dl",&fNormDecayLength);
  for(unsigned int iProng=0; iProng<fNProngs; iProng++){
    fTreeVar->Branch(Form("imp_par_prong%d",iProng),&fImpParProng[iProng]);
    fTreeVar->Branch(Form("imp_par_err_prong%d",iProng),&fImpParErrProng[iProng]);
  }

  //set single-track variables
  //  AddSingleTrackBranches();//need different variables, TPCcluster, cross rows chi2 ITS, ITs clusters
  //TPC -> MCH, ITS -> MFT, PID -> MID

  //set PID variables
  //to be replaced
  //  if(fPidOpt!=kNoPID) AddPidBranches(true,true,false,true,true); 

  return fTreeVar;
}

//________________________________________________________________
bool AliDQTreeHandlerDiMuons::SetVariables(int runnumber, unsigned int eventID, float ptgen, AliAODDimuon* cand, float bfield) 
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
  fInvMass=cand->M();
  fPt=cand->Pt();
  fP=cand->P();
  fY=cand->Y();//assumes rec. pair mass
  fEta=cand->Eta();
  fPhi=cand->Phi();
  fDecayLength=0.0; //cand->DecayLength(); to be implemented in AliAODMuon class
  fNormDecayLengthZ=0.0;//idem
  fDecayLengthZ=0.0;//idem
  fNormDecayLengthZ=0.0;//idem
  fCosP=0.0;//idem
  fCosPz=0.0;//idem
  fXF=cand->XF();
  fCostCS=cand->CostCS();
  fCostHe=cand->CostHe();
  fPhiCS=cand->PhiCS();
  fPhiHe=cand->PhiHe();
  fNormjpsiMeasMinusExp=0.0;///to be checked based on HFcode
  fPseudopropertimeZ=0.0;///to be implemented
  fPseudopropertimeZres=0.0;///to be implemented
  fPseudopropertimeXY=0.0;///to be implemented
  fPseudopropertimeXYres=0.0;///to be implemented
  
  //single track variables
  AliAODTrack* prongtracks[2];
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) prongtracks[iProng] = (AliAODTrack*)cand->GetMu(iProng);//what defines order?
 
  bool setsingletrack = SetSingleTrackVars(prongtracks);  
  if(!setsingletrack) return false;

  //pid variables
  if(fPidOpt==kNoMID) return true;

  bool setpid = SetPiDVars(prongtracks);
  if(!setpid) return false;

  return true;
}
//________________________________________________________________
TTree* AliDQTreeHandlerDiMuons::BuildTreeMCGen(TString name, TString title) {

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
  fTreeVar->Branch("charge_cand",&fCharge);
  fTreeVar->Branch("daughters_in_acc",&fDauInAcceptance);
  //TODO  fTreeVar->Branch("absPDGmother",&fabsPDGProngmother);
  //fTreeVar->Branch("absPDGgdmother",&fabsPDGgdmother);
  return fTreeVar;
  
}
//________________________________________________________________
bool AliDQTreeHandlerDiMuons::SetMCGenVariables(int runnumber, unsigned int eventID, AliAODMCParticle* mcpart, unsigned int pdgmother, unsigned int pdggdmother){
  if(!mcpart) return false;
  if(!(fCandType&kSignal)) return true;

  //TODO: pass also daugther array
  fRunNumber = runnumber;
  fEvID = eventID;
  fPt = mcpart->Pt();
  fY = mcpart->Y();
  fEta = mcpart->Eta();
  fPhi = mcpart->Phi();
  //  fabsPDGmother = pdgmother;
  //  fabsPDGgdmother = pdggdmother;
  
  return true;
}
//________________________________________________________________
void AliDQTreeHandlerDiMuons::SetCandidateType(bool issignal, bool isbkg, bool isprompt, bool isbeauty, bool ischarmonia, bool isbottomonia, bool isew){
  if(issignal) fCandType |= kSignal;
  else fCandType &= ~kSignal;
  if(isbkg && ! fIsMCGenTree) fCandType |= kBkg;
  else fCandType &= ~kBkg;
  if(isprompt) fCandType |= kPrompt;
  else fCandType &= ~kPrompt;
  if(isbeauty) fCandType |= kBeauty;
  else fCandType &= ~kBeauty;
  if(ischarmonia) fCandType |= kCharmonium;
  else fCandType &= ~kCharmonium;
  if(isbottomonia) fCandType |= kBottomonium;
  else fCandType &= ~kBottomonium;
  if(isew && !fIsMCGenTree) fCandType |= kEW;
  else fCandType &= ~kEW;
  
}
//________________________________________________________________
void AliDQTreeHandlerDiMuons::AddSingleTrackBranches(){

  if(fSingleTrackOpt==kNoSingleTrackVars) return;

  for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
    //question how is ordering done? why not with charge?
    if(fSingleTrackOpt==kRedSingleTrackVars) {
      fTreeVar->Branch(Form("charge_prong%d",iProng), &fCharge[iProng]);
      fTreeVar->Branch(Form("pt_prong%d",iProng), &fPtProng[iProng]);
      fTreeVar->Branch(Form("eta_prong%d",iProng), &fEtaProng[iProng]);
      fTreeVar->Branch(Form("phi_prong%d",iProng), &fPhiProng[iProng]);
      fTreeVar->Branch(Form("p_prong%d",iProng), fPProng[iProng]);
      fTreeVar->Branch(Form("dca_prong%d",iProng), fDCA[iProng]);
    }
    else if(fSingleTrackOpt==kAllSingleTrackVars) {
      fTreeVar->Branch(Form("charge_prong%d",iProng), &fCharge[iProng]);
      fTreeVar->Branch(Form("pt_prong%d",iProng),&fPtProng[iProng]);
      fTreeVar->Branch(Form("pt_mch_prong%d",iProng),&fMCHPtProng[iProng]);
      fTreeVar->Branch(Form("eta_prong%d",iProng), &fEtaProng[iProng]);
      fTreeVar->Branch(Form("eta_mch_prong%d",iProng), &fMCHEtaProng[iProng]);
      fTreeVar->Branch(Form("phi_prong%d",iProng), &fPhiProng[iProng]);
      fTreeVar->Branch(Form("phi_mch_prong%d",iProng), &fMCHPhiProng[iProng]);
      fTreeVar->Branch(Form("p_prong%d",iProng), &fPProng[iProng]);
      fTreeVar->Branch(Form("p_mch_prong%d",iProng), &fMCHPProng[iProng]);
      fTreeVar->Branch(Form("dca_prong%d",iProng), fDCA[iProng]);
      fTreeVar->Branch(Form("dcaz_prong%d",iProng), fDCAz[iProng]);
      fTreeVar->Branch(Form("MFTcls_prong%d",iProng), &fNMFTclsProng[iProng]);
      fTreeVar->Branch(Form("MFTclsmap_prong%d",iProng), &fMFTclsMapProng[iProng]);
      fTreeVar->Branch(Form("MCHcls_prong%d",iProng), &fNMCHclsProng[iProng]);
      fTreeVar->Branch(Form("imPar_prong%d",iProng), &fImpParProng[iProng]);
      fTreeVar->Branch(Form("imParErr_prong%d",iProng), &fImpParErrProng[iProng]);
      fTreeVar->Branch(Form("chi2perndf_mft_prong%d",iProng), &fMFTChi2perNDF[iProng]);
      fTreeVar->Branch(Form("chi2perndf_mch_prong%d",iProng), &fMuonChi2perNDF[iProng]);
      fTreeVar->Branch(Form("ratabsorber_prong%d",iProng), &fRatAbsorberEnd[iProng]);
      fTreeVar->Branch(Form("hasMFT_prong%d",iProng), &fHasMFTProng[iProng]);
      fTreeVar->Branch(Form("absPDG_prong%d",iProng), &fabsPDGProngmother[iProng]);
      fTreeVar->Branch(Form("absPDG_mother_prong%d",iProng), &fabsPDGProng[iProng]);
      fTreeVar->Branch(Form("absPDG_gdmother_prong%d",iProng), &fabsPDGPronggdmother[iProng]);
    }
  }
}
//________________________________________________________________
void AliDQTreeHandlerDiMuons::AddPiDBranches(){
  if(fPidOpt==kNoMID) return;
  //TODO add further variables  to do BDT

  for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
    fTreeVar->Branch(Form("muon_midpid_prong%d",iProng),&fMIDPID[iProng]);
    fTreeVar->Branch(Form("muon_midchi2",iProng), &fMuonMatchChi2perNDF[iProng]);
  }
}
//________________________________________________________________
bool AliDQTreeHandlerDiMuons::SetSingleTrackVars(AliAODTrack* tracks[]){

  if(fSingleTrackOpt==kNoSingleTrackVars) return true;
  
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
    if(!tracks[iProng]) {
      AliWarning("Prong track not found!");
      return false;
    }
  }

  for(unsigned int iProng=0; iProng<fNProngs; iProng++){
  
  if(fSingleTrackOpt==kRedSingleTrackVars) {
    fCharge[iProng]=tracks[iProng]->Charge();
    fPtProng[iProng]=tracks[iProng]->Pt();
    fEtaProng[iProng]=tracks[iProng]->Eta();
    fPhiProng[iProng]=tracks[iProng]->Phi();
    fPProng[iProng]=tracks[iProng]->P();
    fDCA[iProng]=tracks[iProng]->DCA();
    
  }
  else if(fSingleTrackOpt==kAllSingleTrackVars) {
    fCharge[iProng]=tracks[iProng]->Charge();
    fPtProng[iProng]=tracks[iProng]->Pt();
    fMCHPtProng[iProng]=0.0;//TODO to be seen
    fEtaProng[iProng]=tracks[iProng]->Eta();
    //TODO MCH eta?
    fPhiProng[iProng]=tracks[iProng]->Phi();
    fMCHPhiProng[iProng]=tracks[iProng]->Eta();
    fPProng[iProng]=tracks[iProng]->P();
    fMCHPProng[iProng]=0.0;//TODO to be seen
    fDCA[iProng]=tracks[iProng]->DCA();
    fDCAz[iProng]=tracks[iProng]->ZAtDCA();//TODO is this the coordinate or the distance to PV?
    fNMFTclsProng[iProng]=0.0;//TODO to be seen
    fMFTclsMapProng[iProng]=0.0;//TODO to be seen
    fNMCHclsProng[iProng]=0.0;//TODO to be seen
    fImpParProng[iProng]=0.0;//TODO D0
    fImpParErrProng[iProng]=0.0;//TODO D0
    fMFTChi2perNDF[iProng]=0.0;//TODO D0
    fMuonChi2perNDF[iProng]=tracks[iProng]->Chi2perNDF();//TODO to be seen
    fRatAbsorberEnd[iProng]=tracks[iProng]->GetRAtAbsorberEnd();
    fHasMFTProng[iProng]=tracks[iProng]->IsMuonGlobalTrack();

  }

  }
  return true;
  
}
//________________________________________________________________
bool AliDQTreeHandlerDiMuons::SetPiDVars(AliAODTrack* tracks[]){

 for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
   if(!tracks[iProng]) {
     AliWarning("Prong track not found!");
     return false;
   }
 }
  for(unsigned int iProng=0; iProng<fNProngs; iProng++){
    fMIDPID[iProng] = tracks[iProng]->GetMatchTrigger();
    fMuonMatchChi2perNDF[iProng] = (float) tracks[iProng]->GetChi2MatchTrigger();
 
  }  
  return true;
}


