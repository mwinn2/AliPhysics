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

#include <TTree.h>
#include "AliAODTrack.h"
#include "AliAODMCParticle.h" 

class AliDQTreeHandlerSingleMuons : public TObject
{
  public:

  enum candtype{
    kMFTmatch //empty for now for Run 2
  };
  
  
  //optPID something about MID
  enum{
    kNoMID, //trigger
    kMID
  };
  
     enum optsingletrack {
       kNoSingleTrackVars, // single-track vars off
       kRedSingleTrackVars, // only pT, p, eta, phi
       kAllSingleTrackVars // all single-track vars
     };
    
    AliDQTreeHandlerSingleMuons();

    virtual ~AliDQTreeHandlerSingleMuons();

    TTree* BuildTree(TString name, TString title);
    bool SetVariables(int runnumber, unsigned int eventID, float ptgen, AliAODTrack* cand, float bfiled, int masshypo=0);

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
      if(isselectedTopo) fCandType |= kSelectedTopo;
      else fCandType & = ~kSelectedTopo;
      if(isselectedPID) fCandType |= kSelectedPID;
      else fCandType & = ~kSelectedPID;
      if(isselectedTracks) fCandType |= kSelectedTracks;
      else fCandType &= ~kSelectedTracks;
    }

    void SetInAcceptance(bool inacc = true) {fInAcceptance=inacc;}
    
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
    static bool IsPrompt(int candtype) {
      if(candtype>>3&1) return true;
      return false;
    }
    static bool IsFD(int candtype) {
      if(candtype>>4&1) return true;
      return false;
    }
    //    static bool IsRefl(int cand) makes no sense for single
    static bool IsSelectedStdTopo(int candtype) {
      if(candtype>>6&1) return true;
      return false;
    }
    static bool IsSelectedStdPID(int candtype) {
      if(candtype>>7&1) return true;
      return false;
    }
    static bool IsSelectedStdTracks(int candtype) {
      of(candtype>>8&1) return true;
      return false;
    }

 protected:

    //constant variables
    static const unsgined int knMaxProngs = 3;
    //Prong: Muon + MFT standalone tracklet

    
    const float kCSPEED = 2.99792457999999984e-02; // cm / ps

    //TODO
    //    void AddSingleTrackBranches();
    //  void AddPIDBranches();
    //TODO
    //    bool SetSingleTrackVars(AliAODTrack* tracks[]);
    // bool SetPidVars(AliAODTrack* tracks[]);

    
    TTree* fTreeVar; /// tree with variables
    unsigned int fNProngs;///number of prongs: muon + MFT trackles
    unsigned int fNCandidates; /// number of candidates in one fill (event)
    inf fCandType; /// flag for candidate type (bit map above)
    float fInvMass; ///candidate invariant mass: for MFT tracklet + muon: could estimate momentum from transverse plane balance of both tracks w.r.t. line of flight -> get momentum tracklet, get mass of 2-track system, similar for following variables...
    float fPt; ///candidate pt
    float fPtGen; ///generated candidate pt
    float fY; ///candidate rapidity
    float fEta; ///candidate pseudorapidity
    float fPhi; //candidate azimuthal angle
    float fDecayLength; ///candidate decay length, only with at least 2 prongs useful
    float fDecayLengthNorm; ///candidate decay length normalised, only with at least 2 prongs
    float fDecayLengthZ; ///candidate decay length in the longitudinal direction
    float fNormDecayLengthZ; ///candidate normalised decay length in the beam direction
    //    float fCosP; ///candidate cosine of pointing angle (makes only sense if some dummy hypothesis about momenta of MFT tracks, since semileptonic, not too much sense.), same for impact parameter on candidate level, rather muons
    float fMuonDCA; //DCA of muon
    float fMuonP; ///muon momentum
    float fMuonPt;//muon pt
    float fMuonEta;///Muon pseudorapidity
    float fMuonPhi;///Muon phi
    float fMuonChi2perNDF;///muon track chi2/ndf
    float fMFTChi2perNDFProng[knMaxProngs];///MFT chi2perNDF for all prongs
    int fMFTclsProng[knMaxProngs];///number of prong MFT clusters
    int fMuonMCHcls;///number of Muon chamber cluster
    float fMIDPID;//just one variable for Muon ID to be expanded! should be a neural network, which optimises muon-id
    int fPidOpt;///option for PID variables
    int fSingleTrackOpt; ///option for single-track variables
    bool fFillOnlySignal; ///flag to enable only signal filling
    bool fIsMCGenTree; ///flag to know if is a tree for MC generated particles
    bool fTrackInAcceptance; ///flag to know if the track is in acceptance in case of MC gen
    unsigned int fEvID; ///event ID corresponding to  the one set in fTreeEvChar
    int fRunNumber; ///run number
    int fRunNumberPrevCand; ///run number of previous candidate
    

    float fImpParProng[knMaxProngs]; ///prong impact parameter
    float fImpParErrProng[knMaxProngs]; ///error on prongs rphi impact param [cm]

    /// \cond CLASSIMP
    ClassDef(AliDQTreeHandlerSingleMuons,3); ///  what does the number mean afterwards?
    /// \endcond
};
#endif
