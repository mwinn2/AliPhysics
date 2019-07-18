#ifndef ALIDQTREEHANDLERDIMUONS_H
#define ALIDQTREEHANDLERDIMUONS_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliDQTreeHandlerDiMuons
// \brief helper class to handle a tree for Dimuon cut optimisation and MVA analyses
// \authors:
// M. Winn, mwinn@cern.ch
//heavily based on AliHFTreeHandler
/////////////////////////////////////////////////////////////

#include <TTree.h>
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"

class AliDQTreeHandlerDiMuons : public TObject
{
  public:
    
  enum candtype{//check if name is
    kMFTmatchnone,
    kMFTmatchsingle,//
    kMFTmatchboth//
  };

  enum optpid {//shall this mean, not filled or not required?
    kNoMID,
    kMIDsingle,
    kMIDboth
  };

  enum optstrack {
    kNoSingleTrackVars, //single-track vars off
    kRedSingleTrackVars, // only pt, p, eta, phi
    kAllSingleTrackVars // all single-track vars
  }
    
    AliDQTreeHandlerDiMuons();

    virtual ~AliDQTreeHandlerDiMuons();

    TTree* BuildTree(TString name, TString title);
    bool SetVariables(int runnumber, unsigned int eventID, float ptgen, AliAODRecoDecayHF* cand, float bfiled, int masshypo=0);
    //ToDo check mass hypothesis ting
    
    TTree* BuildTreeMCGen(TString name, TString title);
    bool SetMCGenVariables(int runnumber, unsigned int eventID, AliAODMCParticle* mcpart);

    void FillTree() {
      if(fFillOnlySignal && !(fCandType&kSignal)) {
	fCandType=0;
      }
      else {
	fTreeVar->Fill();
	fCandType=0;
	fRunNumberPrevCand = fRunNumber;
      }
    }

    void SetOptPID(int PIDopt) {fPidOpt=PIDopt;}
    void SetOptSingleTrackVars(int opt) {fSingleTrackOpt=opt;}
    void SetFillOnlySignal(bool fillopt=true) {fFillOnlySignal=fillopt;}

    void SetCandidateType(bool issignal, bool isbkg, bool isprompt, bool isFD, bool isreflected);
    void SetIsSelectedStd(bool isselected, bool isselectedTopo, bool isselectedPID, bool isselectedTracks) {
      if(isselected) fCandType |= kSelected;
      else fCandType &= ~kSelected;
      if(isselectedTopo) fCandType |= fSelectedTopo;
      else fCandType & =~kSelectedPID;
      if(isselectedTracks) fCandType |= kSelectedTracks;
      else fCandType &= ~kSelectedTracks;
    }
      

    void SetInAcceptance(bool inacc = true ) {fInAcceptance=inacc;}

    static bool IsSelectedStd(int candtype) {
      if(candtype&1) return true;
      return false;
    }
    static bool IsSignal(int candtype) {
      if(candtype>>1&1) return true;
      return false;
    }
    static bool IsBkg(int candtype){
      if(candtype>>2&1) return true;
      return false;
    }
    static bool IsPrompt(int candtype){
      if(candtype>>3&1) return true;
      return false;
    }
    static bool IsFD(int candtype){
      if(candtype>>4&1) return true;
      return false;
    }

    static bool IsSelectedStdTopo(int candtype){
      if(candtype>>6&1) return true;
      return false;
    }

    static bool IsSelectedStdPID(int candtype){
      if(candtype>>7&1) return true;
    }
    static bool IsSelectedStdTracks(int candtype) {
      if(candtype>>8&1) return true;
      return false;
    }
    
    
  protected:

    //constant variables
    statc const unsigned int knMaxProngs = 5;
    //Pront: 2x Muons + MFT standalone tracklet

    const float kCSPEED = 2.99792457999999984e-02;

    void AddSingleTrackBranches();
    void AddPIDBranches();
    bool SetSingleTrackVars(AliAODTrack* tracks[]);
    bool SetPIDVars(AliAODTrack* tracks[]);

    TTree* fTreeVar; //tree with variables
    unsigned int fNProngs; ///number of prongs: muon + MFT tracklets
    unsigned int fNCandidates; /// number of candidates in one fill (event)
    unsigned int fCandType; /// flag for candidate type (bit map above)
    float fInvMass; ///candidate invariant mass: by default dimuon mass
    float fPt; //candidate pt
    float fPtGen; ///generated candidate pt
    float fY; ///candidate rapidity
    float fEta; ///candidate pseudorapidity
    float fPhi; //candidate azimuthal angle
    float fDecayLength; ///candidate decay length
    float fDecayLengthNorm; ///candidate decay length normalised, only with at least two tracks
    float fDecayLengthZ;////candidate decay length in the longitudinal
    float fNormDecayLengthZ; ///candidate normalised decay length in the beam direction
    float fCosP; //candidate cosine of pointing angle
    float fCosPz; ///candidate cosine of pointing angle in the longitudinal direction
    //CONTINUE HERE
    float fImpParZ; ///candidate impact parameter in the longitudinal direction
    float fDCA[knMaxProngs]; /// DCA of candidates prongs
    float fPProng[knMaxProngs]; ///prong momentum
    float fMCHPProng[knMaxProngs]; ///prong MCH momentum
    float fPtProng[knMaxProngs]; ///prong pt
    float fEtaProng[knMaxProngs];///prong pseudorapidity
    float fPhiProng[knMaxProngs];///prong azimuthal angle
    int fNMFTclsProng[knMaxProngs];///prong track number of MFT clusters
    int fMFTclsMapProng[knMaxProngs];///prong track MFT cluster map
    float fNMCHclsProng[knMaxProngs];///prong track number of MCH clusters
    int fPidOpt; ///option for PID variables
    float fMIDchi2perNDF; ///PID variable for matching
    ///anything else for PID?
    int fSingleTrackOpt; ///option for single-track variables
    bool fFillOnlySignal; ///flag to enable only signal filling
    bool fIsMCGenTree; ////flag to know if is a tree for MC generated particles
    bool fDauInAcceptance; ///flag to know if the daughter are in acceptance in case of MC gen
    unsigned int fEvID; ///event ID correspodning to the one set in FTreeEvChar
    int fRunNumber; ////run number
    int fRunNumberPrevCand; ///run number of previous candidate
    
    float fImpParProng[knMaxProngs]; ///prong impact parameter
    float fCosThetaStar; /// candidate costhetastar
    float fImpParProd; /// daughter impact-parameter product
    float fNormjpsiMeasMinusExp; ///candidate topomatic variable
    float fImpParErrProng[knMaxProngs]; ///error on prongs z impact param [cm]

    /// \cond CLASSIMP
    ClassDef(AliHFTreeHandlerD0toKpi,3); /// 
    /// \endcond
};
#endif
