/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliDQTreeHandlerDiMuons
// \brief helper class to handle a tree for D0 cut optimisation and MVA analyses
// \authors:
// M. Winn, mwinn@cern.ch
// heavily based on AliHFTreeHandlerD0toKpi
/////////////////////////////////////////////////////////////

#include <cmath>
#include <limits>
#include <TString.h>
#include "AliDQTreeHandlerDiMuons.h"
//#include "AliAODRecoDecayHF2Prong.h": use different framework

/// \cond CLASSIMP
ClassImp(AliDQTreeHandlerDiMuons);
/// \endcond

//________________________________________________________________
AliDQTreeHandlerDiMuons::AliDQTreeHandlerDiMuons():
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
  fCosP(-9999.),
  fCosPz(-9999.),
  fPidOpt(-1),
  fMIDchi2perNDF(-9999.),
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
    fDCA[iProng] = -9999.;
    fPProng[iProng] = 9999.;
    fMCHPProng[iProng] = 9999.;
    fPtProng[iProng] = 9999.;
    fEtaProng[iProng] = 9999.;
    fPhiProng[iProng] = 9999.;
    fNMFTclsProng[iProng] = -1;
    fMFTclsMapProng[iProng] = -1;
    fNMCHclsProng[iProng] = -1;
    fImpParProng[iProng] = -9999.;
    fImpParErrProng[iProng] = -9999.;
  }
}
//________________________________________________________________
AliHFTreeHandlerDiMuons::~AliHFTreeHandlerDiMuons()
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
  fTreeVar->Branch("max_norm_d0d0exp",&fNormd0MeasMinusExp);//TODO change naming
  fTreeVar->Branch("norm_dl",&fNormDecayLength);
  for(unsigned int iProng=0; iProng<fNProngs; iProng++){
    fTreeVar->Branch(Form("imp_par_prong%d",iProng),&fImpParProng[iProng]);
    fTreeVar->Branch(Form("imp_par_err_prong%d",iProng),&fImpParErrProng[iProng]);
  }

  //set single-track variables
  AddSingleTrackBranches();//need different variables, TPCcluster, cross rows chi2 ITS, ITs clusters
  //TPC -> MCH, ITS -> MFT, PID -> MID

  //set PID variables
  //to be replaced
  //  if(fPidOpt!=kNoPID) AddPidBranches(true,true,false,true,true); 

  return fTreeVar;
}

//________________________________________________________________
bool AliHFTreeHandlerDiMuons::SetVariables(int runnumber, unsigned int eventID, float ptgen, AliAODRecoDecayHF* cand, float bfield, int masshypo, AliPIDResponse *pidrespo) //what is AliAODRecoDecayHF? to be replaced by AliAODRecoDimuon?
//needs to be changed to standard Dimuons, is there already an object like that?
{
  fIsMCGenTree=false;

  if(!cand) return false;
  if(fFillOnlySignal) { //if fill only signal and not signal candidate, do not store
    if(!(fCandType&kSignal || fCandType&kRefl)) return true;
  }
  fNCandidates++;
  fRunNumber=runnumber;
  fEvID=eventID;
  fPtGen=ptgen;
  
  //topological variables
  //common
  fPt=cand->Pt();// TODO all to be changed from a dimuon
  fY=cand->Y(421);//not clear
  fEta=cand->Eta();
  fPhi=cand->Phi();
  fDecayLength=cand->DecayLength();
  fDecayLengthXY=cand->DecayLengthXY();
  fNormDecayLengthXY=cand->NormalizedDecayLengthXY();
  fCosP=cand->CosPointingAngle();
  fCosPXY=cand->CosPointingAngleXY();
  fImpParXY=cand->ImpParXY();
  fDCA=cand->GetDCA();
  fNormd0MeasMinusExp=ComputeMaxd0MeasMinusExp(cand,bfield);
  fNormDecayLength=cand->NormalizedDecayLength();
  
  //D0 -> Kpi variables
  fImpParProd=((AliAODRecoDecayHF2Prong*)cand)->Prodd0d0();
  if(masshypo==0) { //D0 -> Kpi
    fInvMass=((AliAODRecoDecayHF2Prong*)cand)->InvMassD0();
    fCosThetaStar=((AliAODRecoDecayHF2Prong*)cand)->CosThetaStarD0();
  }
  else if(masshypo==1) { //D0 -> piK
    fInvMass=((AliAODRecoDecayHF2Prong*)cand)->InvMassD0bar();
    fCosThetaStar=((AliAODRecoDecayHF2Prong*)cand)->CosThetaStarD0bar();
  }
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
    fImpParProng[iProng]=cand->Getd0Prong(iProng);
    fImpParErrProng[iProng]=cand->Getd0errProng(iProng);
  }
    
  //single track variables
  AliAODTrack* prongtracks[2];
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) prongtracks[iProng] = (AliAODTrack*)cand->GetDaughter(iProng);
  bool setsingletrack = SetSingleTrackVars(prongtracks);  
  if(!setsingletrack) return false;

  //pid variables
  if(fPidOpt==kNoPID) return true;

  bool setpid = SetPidVars(prongtracks,pidrespo,true,true,false,true,true);
  if(!setpid) return false;

  return true;
}

//________________________________________________________________
void AliHFTreeHandlerDiMuons::SetIsDzeroDzeroBar(int isSel, int isSelTopo, int isSelPID, int isSelFilt, int isSelTopoFilt, int isSelPIDFilt) {
    
    //analysis cuts
    //combined selection PID topo
    if(isSel==0){
        fCandType &= ~kDzeroComb;
        fCandType &= ~kDzeroBarComb;
    }
    else if(isSel==1){
        fCandType |= kDzeroComb;
        fCandType &= ~kDzeroBarComb;
    }
    else if(isSel==2){
        fCandType |= kDzeroBarComb;
        fCandType &= ~kDzeroComb;
    }
    else if(isSel==3){
        fCandType |= kDzeroComb;
        fCandType |= kDzeroBarComb;
    }
    //topol selection
    if(isSelTopo==0){
        fCandType &= ~kDzeroTopo;
        fCandType &= ~kDzeroBarTopo;
    }
    else if(isSelTopo==1){
        fCandType |= kDzeroTopo;
        fCandType &= ~kDzeroBarTopo;
    }
    else if(isSelTopo==2){
        fCandType |= kDzeroBarTopo;
        fCandType &= ~kDzeroTopo;
    }
    else if(isSelTopo==3){
        fCandType |= kDzeroTopo;
        fCandType |= kDzeroBarTopo;
    }
    //PID selection
    if(isSelPID==0){
        fCandType &= ~kDzeroPID;
        fCandType &= ~kDzeroBarPID;
    }
    else if(isSelPID==1){
        fCandType |= kDzeroPID;
        fCandType &= ~kDzeroBarPID;
    }
    else if(isSelPID==2){
        fCandType |= kDzeroBarPID;
        fCandType &= ~kDzeroPID;
    }
    else if(isSelPID==3){
        fCandType |= kDzeroPID;
        fCandType |= kDzeroBarPID;
    }
    
    //filtering cuts
    //combined selection PID topo
    if(isSelFilt==0){
        fCandType &= ~kDzeroCombFilt;
        fCandType &= ~kDzeroBarCombFilt;
    }
    else if(isSelFilt==1){
        fCandType |= kDzeroCombFilt;
        fCandType &= ~kDzeroBarCombFilt;
    }
    else if(isSelFilt==2){
        fCandType |= kDzeroBarCombFilt;
        fCandType &= ~kDzeroCombFilt;
    }
    else if(isSelFilt==3){
        fCandType |= kDzeroCombFilt;
        fCandType |= kDzeroBarCombFilt;
    }
    //topol selection
    if(isSelTopoFilt==0){
        fCandType &= ~kDzeroTopoFilt;
        fCandType &= ~kDzeroBarTopoFilt;
    }
    else if(isSelTopoFilt==1){
        fCandType |= kDzeroTopoFilt;
        fCandType &= ~kDzeroBarTopoFilt;
    }
    else if(isSelTopoFilt==2){
        fCandType |= kDzeroBarTopoFilt;
        fCandType &= ~kDzeroTopoFilt;
    }
    else if(isSelTopoFilt==3){
        fCandType |= kDzeroTopoFilt;
        fCandType |= kDzeroBarTopoFilt;
    }
    //PID selection
    if(isSelPIDFilt==0){
        fCandType &= ~kDzeroPIDFilt;
        fCandType &= ~kDzeroBarPIDFilt;
    }
    else if(isSelPIDFilt==1){
        fCandType |= kDzeroPIDFilt;
        fCandType &= ~kDzeroBarPIDFilt;
    }
    else if(isSelPIDFilt==2){
        fCandType |= kDzeroBarPIDFilt;
        fCandType &= ~kDzeroPIDFilt;
    }
    else if(isSelPIDFilt==3){
        fCandType |= kDzeroPIDFilt;
        fCandType |= kDzeroBarPIDFilt;
    }    
}
