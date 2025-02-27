#ifndef ALI_ANALYSIS_TASK_SE_LBTOLBPI4_H
#define ALI_ANALYSIS_TASK_SE_LBTOLBPI4_H

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TH2F.h>
#include <TArrayD.h>
#include "TClonesArray.h"
#include "AliAnalysisTaskSE.h"
#include "AliRDHFCutsLctopKpi.h"
#include "AliRDHFCuts.h"
#include "AliAODMCHeader.h"
#include "AliVertexingHFUtils.h"

#include "AliAODRecoDecayHF3Prong.h"
#include <TH1F.h>

//#include "TMVAClassification_BDT_pt4to7.class.C"
//#include "TMVAClassification_BDT_pt7to10.class.C"
//#include "TMVAClassification_BDT_pt10to14.class.C"
//#include "TMVAClassification_BDT_pt14to9999.class.C"
//#include "IClassifierReader.C" 


class TNtuple;
class TGraph;
class TList;
class AliAODTrack;
class TClonesArray;
class TObjArray;
class AliESDVertex;
class AliVVertex;

class AliAnalysisTaskSELbtoLcpi4:public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskSELbtoLcpi4();
  AliAnalysisTaskSELbtoLcpi4(const char *name,
    Bool_t fillntuple,
    AliRDHFCutsLctopKpi *lccutsanal,
    AliRDHFCutsLctopKpi *lccutsprod,
    Int_t ndebug);
  virtual ~AliAnalysisTaskSELbtoLcpi4();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  // virtual void UserCreateOutputHistos();
  virtual void Init();
  //  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  // options to fill Ntuple
  void SetFillNtupleSignal(Bool_t val = kTRUE) {fFillNtupleSignal = val;}
  void SetFillNtupleBackgroundRotated(Bool_t val = kTRUE) {fFillNtupleBackgroundRotated = val;}
  void SetFillNtupleBackgroundNonRotated(Bool_t val = kTRUE) {fFillNtupleBackgroundNonRotated = val;}

 private:
  AliAnalysisTaskSELbtoLcpi4(const AliAnalysisTaskSELbtoLcpi4&);
  AliAnalysisTaskSELbtoLcpi4& operator=(const AliAnalysisTaskSELbtoLcpi4&); 

  // Helper functions
  //AliESDVertex* RecalculateVertex(const AliVVertex *old,TObjArray *tracks,Double_t bField);
  //aggiunto
  void AddDaughterRefs(AliAODVertex *v,const AliVEvent *event, const TObjArray *trkArray) const;
  
  Int_t CheckMCLc(AliAODRecoDecayHF3Prong *d, TClonesArray* arrayMC);
  Int_t CheckMCpartPIONaf(AliAODTrack *p, TClonesArray* arrayMC);
  void  FillLbHists(AliAODRecoDecayHF2Prong *part,Int_t lb,AliAODMCHeader *mcHeader,TClonesArray* arrayMC,AliAODTrack *p,AliAODRecoDecayHF3Prong *d);
  void  FillLbHistsnr(AliAODRecoDecayHF2Prong *part,Int_t lb,AliAODMCHeader *mcHeader,TClonesArray* arrayMC,AliAODTrack *p,AliAODRecoDecayHF3Prong *d);
  void FillHistos(AliAODRecoDecayHF3Prong* d,TClonesArray* arrayMC,AliAODEvent *ev,AliAODMCHeader *mcHeader);
  Bool_t CheckGenerator(AliAODTrack *p,AliAODRecoDecayHF3Prong *d,AliAODMCHeader *mcHeader,TClonesArray* arrayMC);
  Int_t IsSelectedLbMY(TObject* obj,Int_t selectionLevel,Int_t lb,Int_t isRot, Bool_t isHijing) const;
  void CheckMCKine(TClonesArray *mcs); 
  Bool_t CountLc(AliAODRecoDecayHF3Prong* Lc, AliAODTrack* pion, TClonesArray* arrayMC,Int_t motherLabelLc,Int_t motherLabelpione);
  Bool_t IsCandidateInjected(AliAODRecoDecayHF *part, AliAODMCHeader *header,TClonesArray *arrayMC);
  Int_t IsTrackInjected(AliAODTrack *part,AliAODMCHeader *header,TClonesArray *arrayMC);
  TObjArray* GetArrayCandRotated(AliAODEvent* ev,AliAODRecoDecayHF2Prong *decay,TClonesArray* arrayMC,Int_t nRot);
  AliAODVertex* RecalculateVertex(const AliVVertex *primary,TObjArray *tracks,Double_t bField);


  AliPIDResponse *fPIDResponse;
  TList   *fOutput; //! list send on output slot 3
  TH1F    *fHistNEvents; //!hist. for No. of events 
  TH1F    *fHistNEventsCuts; //!hist. for No. of events cuts
  TH1F    *fHistNEventsCutsLb; //!hist. for No. of events Lb
  AliRDHFCutsLctopKpi *fRDCutsAnalysisLc; //Cuts for Analysis
  AliRDHFCutsLctopKpi *fRDCutsProductionLb; //Production Cuts
  TList *fListCuts; //list of cuts


  Double_t fBzkG;                     // z component of magnetic field
  AliAODVertex *fvtx1;                // primary vertex
  TList    *fOutputList;            //! output slot 8
  TH1F     *fInvMassLbSign0;     //!histogram with invariant mass of Lb signal its upgrade
  TH1F     *fInvMassLbSign1;     //!histogram with invariant mass of Lb signal its upgrade
  TH1F     *fInvMassLbSign2;     //!histogram with invariant mass of Lb signal its upgrade
  TH1F     *fInvMassLbSign3;     //!histogram with invariant mass of Lb signal its upgrade
  TH1F     *fInvMassLbSign4;     //!histogram with invariant mass of Lb signal its upgrade
  TH1F     *fInvMassLbSign5;     //!histogram with invariant mass of Lb signal its upgrade
  TH1I      *fSelMC;
  TH1I      *fCountLc;
   TNtuple *fNtupleLambdacUPG; //! output ntuple  
   //TNtuple *fNtupleDiffD0rot; //! output ntuple 
  //TH2F *fMassHistBDT[8*kMaxPtBins]; //!hist. for inv mass vs BDT response


  //IClassifierReader *fBDTReader[4]; //!BDT reader for standalone class
  AliAODVertex* ReconstructSecondaryVertex(TObjArray *trkArray,Double_t &dispersion,Bool_t useTRefArray=kTRUE) const; //Reconstruct vertex of B

  Bool_t fFillNtupleSignal; /// flag to fill ntuple with signal candidates
  Bool_t fFillNtupleBackgroundRotated; /// flag to fill ntuple with background rotated candidates
  Bool_t fFillNtupleBackgroundNonRotated; /// flag to fill ntuple with background non rotated candidates

  ClassDef(AliAnalysisTaskSELbtoLcpi4,2);
};

#endif

