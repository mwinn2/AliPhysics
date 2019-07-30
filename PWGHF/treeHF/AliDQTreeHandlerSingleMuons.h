#ifndef ALIDQTREEHANDLERSINGLEMUONS_H
#define ALIDQTREEHANDLERSINGLEMUONS_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliDQTreeHandlerSingleMuons
// \brief helper class to handle a tree for Dimuon cut optimisation and MVA analyses
// \authors:
// M. Winn, mwinn@cern.ch
//heavily based on AliHFTreeHandlerD0toKpi
/////////////////////////////////////////////////////////////

#include <TTree.h>
#include "AliAODTrack.h"
#include "AliAODMCParticle.h" 

class AliDQTreeHandlerSingleMuons : public TObject
{
  public:

  enum candtype{//TODO introduce things that make sense for single muons
    kSelected         = BIT(0),
    kSignal           = BIT(1),
    kBkg              = BIT(2),
    kCharm            = BIT(3),
    kBeauty           = BIT(4),
    kEW               = BIT(5),
    kQuarkonia        = BIT(6),
    kSelectedTopo     = BIT(7),
    kSelectedPID      = BIT(8),
    kSelectedTracks   = BIT(9),
  };
  
  
  //optPID something about MID
  //for Run 2 just trigger match
  enum optpid {// check if 
    kNoMID, //trigger
    kMID
  };
  
  enum optsingletrack {//TODO: does it make sense?
    kNoSingleTrackVars, // single-track vars off
    kSingleTrackVars // only pT, p, eta, phi
  };
    
    AliDQTreeHandlerSingleMuons();

    virtual ~AliDQTreeHandlerSingleMuons();

    TTree* BuildTree(TString name, TString title);
    bool SetVariables(int runnumber, unsigned int eventID, float ptgen, AliAODTrack* cand, float bfield);

    TTree* BuildTreeMCGen(TString name, TString title);
    bool SetMCGenVariables(int runnumber, unsigned int eventID, AliAODMCParticle* mcpart, unsigned int pdgmother, unsigned int pdggdmother);

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

    void SetCandidateType(bool issignal, bool isbkg, bool ischarm, bool isbeauty, bool isew, bool isquarkonia);
    void SetIsSelectedStd(bool isselected, bool isselectedTopo, bool isselectedPID, bool isselectedTracks) {
      if(isselected) fCandType |= kSelected;
      else fCandType &= ~kSelected;
      if(isselectedTopo) fCandType |= kSelectedTopo;
      else fCandType &= ~kSelectedTopo;
      if(isselectedPID) fCandType |= kSelectedPID;
      else fCandType &= ~kSelectedPID;
      if(isselectedTracks) fCandType |= kSelectedTracks;
      else fCandType &= ~kSelectedTracks;
    }

    void SetInAcceptance(bool inacc = true) {fTrackInAcceptance=inacc;}
    
    static bool IsSelectedStd(int candtype) {
      if(candtype&1) return true;
      return false;
    }
    static bool IsSignal(int candtype) {
      if(candtype>>1&1) return true;
      return false;
    }
    static bool IsBkg(int candtype) {
      if(candtype>>2&1) return true;
      return false;
    }
    static bool IsCharm(int candtype) {
      if(candtype>>3&1) return true;
      return false;
    }
    static bool IsBeauty(int candtype) {
      if(candtype>>4&1) return true;
      return false;
    }
    static bool IsQuarkonia(int candtype) {
      if(candtype>>5&1) return true;
      return false;
    }
    static bool IsSelectedStdTopo(int candtype) {
      if(candtype>>6&1) return true;
      return false;
    }
    static bool IsSelectedStdPID(int candtype) {
      if(candtype>>7&1) return true;
      return false;
    }
    static bool IsSelectedStdTracks(int candtype) {
      if(candtype>>8&1) return true;
      return false;
    }

 protected:

    
    const float kCSPEED = 2.99792457999999984e-02; // cm / ps

    void AddSingleTrackBranches();
    void AddPIDBranches();
    bool SetSingleTrackVars(AliAODTrack* muon);
    bool SetPidVars(AliAODTrack* track);

    
    TTree* fTreeVar; /// tree with variables
    unsigned int fNCandidates; /// number of candidates in one fill (event)
    int fCandType; /// flag for candidate type (bit map above)
    float fInvMass; ///candidate invariant mass: for MFT tracklet + muon: could estimate momentum from transverse plane balance of both tracks w.r.t. line of flight -> get momentum tracklet, get mass of 2-track system, similar for following variables...
    float fPt; ///candidate pt
    float fP; //candidate momentum
    float fPtGen; ///generated candidate pt
    float fY; ///candidate rapidity
    float fEta; ///candidate pseudorapidity
    float fPhi; //candidate azimuthal angle
    Short_t fCharge;//candidate charge
    float fDCA; //DCA of muon is dcaxy at the moment in aod
    float fDCAxy;//
    float fDCAz;
    float fMuonChi2perNDF;///muon track chi2/ndf
    float fRatAbsorberEnd;///R at end of absorber
    bool fHasMFT;/// muon track has a MFT tracklet
    float fMuonMatchChi2perNDF;//muon chi2 for matching with trigger stations
    float fMFTChi2perNDF;///MFT chi2perNDF for all prongs
    int fMFTcls;///number of prong MFT clusters
    int fMuonMCHcls;///number of Muon chamber cluster
    int fMIDPID;//just one variable for Muon ID to be expanded! should be a neural network, which optimises muon-id
    int fPidOpt;///option for PID variables
    int fSingleTrackOpt; ///option for single-track variables
    bool fFillOnlySignal; ///flag to enable only signal filling
    bool fIsMCGenTree; ///flag to know if is a tree for MC generated particles
    bool fTrackInAcceptance; ///flag to know if the track is in acceptance in case of MC gen
    unsigned int fabsPDGmother; ///abs PDG mother particle number for MC
    unsigned int fabsPDGgdmother;///abs PDG grandmother particle number for MC
    unsigned int fEvID; ///event ID corresponding to  the one set in fTreeEvChar
    int fRunNumber; ///run number
    int fRunNumberPrevCand; ///run number of previous candidate
    

    float fImpPar; ///prong impact parameter
    float fImpParErr; ///error on prongs rphi impact param [cm]

    /// \cond CLASSIMP
    ClassDef(AliDQTreeHandlerSingleMuons,3); ///  what does the number mean afterwards?
    /// \endcond
};
#endif
