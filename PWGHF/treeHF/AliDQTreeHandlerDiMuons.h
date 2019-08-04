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
#include "AliAODDimuon.h"

class AliDQTreeHandlerDiMuons : public TObject
{
  public:
  //TODO add B_c
  enum candtype{//check if name is
   kSelected         = BIT(0),
   kSignal           = BIT(1),
   kBkg              = BIT(2),
   kPrompt           = BIT(3),
   kBeauty           = BIT(4),
   kCharmonium       = BIT(5),
   kBottomonium      = BIT(6),
   kEW               = BIT(7),
   kSelectedTopo     = BIT(8),
   kSelectedPID      = BIT(9),
   kSelectedTracks   = BIT(10),    
  };

  enum optpid {//shall this mean, not filled or not required?
    kNoMID,
    kMID
  };

  enum optstrack {
    kNoSingleTrackVars, //single-track vars off
    kRedSingleTrackVars, // only pt, p, eta, phi
    kAllSingleTrackVars // all single-track vars
  };
    
    AliDQTreeHandlerDiMuons();
    AliDQTreeHandlerDiMuons(int PIDopt);
    
    virtual ~AliDQTreeHandlerDiMuons();

    TTree* BuildTree(TString name, TString title);
    bool SetVariables(int runnumber, unsigned int eventID, float ptgen, AliAODDimuon* cand, float bfield);
    //ToDo check mass hypothesis ting
    
    TTree* BuildTreeMCGen(TString name, TString title);
    bool SetMCGenVariables(int runnumber, unsigned int eventID, AliAODMCParticle* mcpart, unsigned int pdgmother, unsigned int pdggdmother);//todo pass also daugther tracks?

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

    void SetCandidateType(bool issignal, bool isbkg, bool isprompt, bool isbeauty, bool ischarmonia, bool isbottomonia, bool isew);
    void SetIsSelectedStd(bool isselected, bool isselectedTopo, bool isselectedPID, bool isselectedTracks) {
      if(isselected) fCandType |= kSelected;
      else fCandType &= ~kSelected;
      if(isselectedTopo) fCandType |= kSelectedTopo;
      else fCandType &= ~kSelectedPID;
      if(isselectedTracks) fCandType |= kSelectedTracks;
      else fCandType &= ~kSelectedTracks;
    }
      

    void SetInAcceptance(bool inacc = true ) {fDauInAcceptance=inacc;}

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
    static bool IsBeauty(int candtype){
      if(candtype>>4&1) return true;
      return false;
    }

    static bool IsCharmonium(int candtype){
      if(candtype>>5&1) return true;
      return false;
    }

    static bool IsBottomonium(int candtype){
      if(candtype>>6&1) return true;
      return false;
    }

    static bool IsEW(int candtype){
      if(candtype>>7&1) return true;
      return false;
    }

    static bool IsSelectedStdTopo(int candtype){
      if(candtype>>8&1) return true;
      return false;
    }

    static bool IsSelectedStdPID(int candtype){
      if(candtype>>9&1) return true;
    }
    static bool IsSelectedStdTracks(int candtype) {
      if(candtype>>10&1) return true;
      return false;
    }
    
    
  protected:

    
    //constant variables
    static const unsigned int knMaxProngs = 5;
    //Pront: 2x Muons + MFT standalone tracklet

    const float kCSPEED = 2.99792457999999984e-02;

    void AddSingleTrackBranches();
    void AddPiDBranches();
    bool SetSingleTrackVars(AliAODTrack* tracks[]);
    bool SetPiDVars(AliAODTrack* tracks[]);

    TTree* fTreeVar; //tree with variables
    unsigned int fNProngs; ///number of prongs: muon + MFT tracklets
    unsigned int fNCandidates; /// number of candidates in one fill (event)
    unsigned int fCandType; /// flag for candidate type (bit map above)
    float fInvMass; ///candidate invariant mass: by default dimuon mass
    float fPt; //candidate pt
    float fP; // candidate momentum
    float fPtGen; ///generated candidate pt
    float fY; ///candidate rapidity
    float fEta; ///candidate pseudorapidity
    float fPhi; //candidate azimuthal angle
    float fDecayLength; ///candidate decay length
    float fNormDecayLength; ///candidate decay length normalised, only with at least two tracks
    float fDecayLengthZ;////candidate decay length in the longitudinal
    float fNormDecayLengthZ; ///candidate normalised decay length in the beam direction
    float fPseudopropertimeZ;///pseudopropertime along z
    float fPseudopropertimeZres;///pseudopropertime along z, resolution
    float fPseudopropertimeXY;///pseudopropertime along transvers. plane, XY
    float fPseudopropertimeXYres;///pseudopropertime along transvers. plane XY res.
    float fCosP; //candidate cosine of pointing angle
    float fCosPz; ///candidate cosine of pointing angle in the longitudinal direction
    float fImpParZ; ///candidate impact parameter in the longitudinal direction
    float fXF;/// Feynman x
    float fCostCS;/// Cosinus of the Collins-Sope polar decay angle
    float fCostHe; ///Cosinus of the Helicity polar decay angle
    float fPhiCS;////Azimuthal angle in the Hellicity frame
    float fPhiHe;////Azimuthal angle in the Hellicity frame
    ///what Jackson etc. name
    int fCharge[knMaxProngs];
    float fDCA[knMaxProngs]; /// DCA of candidates prongs
    float fDCAz[knMaxProngs];/// DCAz of candidate prongs
    float fPProng[knMaxProngs]; ///prong momentum
    float fMCHPProng[knMaxProngs]; ///prong MCH momentum
    float fPtProng[knMaxProngs]; ///prong pt
    float fMCHPtProng[knMaxProngs];//prong MCH pt
    float fEtaProng[knMaxProngs];///prong pseudorapidity
    float fMCHEtaProng[knMaxProngs];///prong MCH eta
    float fPhiProng[knMaxProngs];///prong azimuthal angle
    float fMCHPhiProng[knMaxProngs];///prong MCH phi angle
    int fNMFTclsProng[knMaxProngs];///prong track number of MFT clusters
    int fMFTclsMapProng[knMaxProngs];///prong track MFT cluster map
    int fNMCHclsProng[knMaxProngs];///prong track MFT cluster map
    int fPidOpt; ///option for PID variables
    bool fMIDPID[knMaxProngs];
    float fMIDchi2perNDF[knMaxProngs]; ///PID variable for matching
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
    float fMFTChi2perNDF[knMaxProngs];///number of prong MFT clusters
    float fMuonChi2perNDF[knMaxProngs];///numberof muon chi2 clusters
    float fMuonMatchChi2perNDF[knMaxProngs];///numberof muon chi2 clusters
    float fRatAbsorberEnd[knMaxProngs];///R at end of absorber
    bool fHasMFTProng[knMaxProngs];///number of prong MFT clusters
    unsigned int fabsPDGProng[knMaxProngs];///abs. pdg number of prong
    unsigned int fabsPDGProngmother[knMaxProngs];/// abs. pdg number of mother
    unsigned int fabsPDGPronggdmother[knMaxProngs];/// abs. pdg gd mother
    /// \cond CLASSIMP
    ClassDef(AliDQTreeHandlerDiMuons,3); /// 
    /// \endcond
};
#endif
