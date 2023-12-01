#ifndef CalculateQC_H
#define CalculateQC_H

#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <complex>
#include <cmath>
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TDirectoryFile.h"
#include "TMatrixD.h"
#include "TTree.h"
#include "TMath.h"
#include "Event.h"


class CalculateFlow
{
public:
  
  CalculateFlow();
  CalculateFlow(const char* name);
  ~CalculateFlow() = default;
  
  virtual void UserCreateOutputObjects();
  virtual void UserExec();
  virtual void Terminate(Int_t Nevents);
  
  virtual void Make(Event* fEvent);
  void InitializeArraysForFlowQC();
  void InitializeArraysForFlowEP();
  void InitializeArraysForFlowRP();
  void InitializeArraysForFlowGF();
  void InitializeArraysForQA();
  virtual void ResetEventByEventQuantities();
  
  virtual void CalculateFlowQC();
  virtual void CalculateFlowEP();
  virtual void CalculateFlowEP1();
  virtual void CalculateFlowRP();
  virtual void CalculateFlowGF();
  virtual void FinalizeFlowQC();
  virtual void FinalizeFlowEP();
  virtual void FinalizeFlowRP();
  virtual void FinalizeFlowGF();
  virtual void FinalizeQA();
  
  Int_t GetCRCCenBin(Double_t Centrality);
  Double_t GetSumPro(TProfile *pro, Int_t bin);
  Double_t GetWeightedCorrelations(TProfile *profile, Int_t pt);
  Double_t GetWeightedCorrelationsError(TProfile *profile, Int_t pt);
  
  void SetCentralityEBE(Double_t const c) {this->fCentralityEBE = c;};
  Double_t GetCentralityEBE() const {return this->fCentralityEBE;};
  void SetEvent(Event* e) {this->fEvent = e;};
  
  TList* GetFlowQCList() {return this->fFlowQCList;};
  TList* GetFlowEPList() {return this->fFlowEPList;};
  TList* GetFlowRPList() {return this->fFlowRPList;};
  TList* GetFlowGFList() {return this->fFlowGFList;};
  
  TList* GetQAList() {return this->fQAList;};
  TList* GetSpectraList() {return this->fSpectraList;} // Particle spectra
  
  Event* GetEvent() {return this->fEvent;};
  virtual std::complex<double> ucN(const Int_t n, const TArrayI& h, Int_t ptb);
  virtual std::complex<double> ucN2(const Int_t n, TArrayI& h, TArrayI& cnt, Int_t ptb);
  
  void SetmaxPtCut(Double_t maxPt) {this->maxPtCut = maxPt;};
  void SetminPtCut(Double_t minPt) {this->minPtCut = minPt;};
  void SetminNtrackCut(Double_t minNtrk) {this->minNtracks = minNtrk;};
  void SetmaxEtaCut(Double_t maxEta) {this->maxEtaCut = maxEta;};
  void SetEtaDiff(Bool_t etaflag) {this->EtaDiff = etaflag;};
  void SetdoQA(Bool_t bflag) {this->doQA = bflag;};
  
private:
  Event* fEvent;
  
  Int_t fCenBin = -1;
  Int_t fCRCnCen = 10;
  Double_t fCentralityEBE;
  Double_t fCenWeightEbE = 1; // In MC, set to 1 for now.
  Double_t wPhiEta = 1; // In MC, set to 1 for now.
  Double_t wPhi = 1; // In MC, set to 1 for now.
  Double_t wPt = 1; // In MC, set to 1 for now.
  Double_t wEta = 1; // In MC, set to 1 for now.
  Double_t wTrack = 1; // In MC, set to 1 for now.
  Double_t *fCRCPtBins;
  Double_t *fBins;
  Bool_t doQA = kFALSE;
  Double_t trkWgt = 1;
  const static Int_t fNParticles = 4; // Incl, Pions, Protons, Kaons
  
  // QA Histograms
  TList *fQAList;
  TH1D *fMultChargedParticlesDistribution;
  TH1D *fNumberOfParticipantsDistribution;
  TH1D *fQ2Distribution;
  TH1D *fQ2TPCDistribution;
  TH1D *fQ2V0ADistribution;
  TH1D *fQ2V0CDistribution;
  TH1D *fPtChargedParticlesDistribution;
  TH1D *fEtaChargedParticlesDistribution;
  TH1D *fPhiChargedParticlesDistribution;
  TH2D *fEtaPhiChargedParticlesDistribution;
  TH1D *fPionsPtSpectra;
  TH1D *fPionsEtaSpectra;
  TH1D *fPionsPhiSpectra;
  TH1D *fPosPionsPtSpectra;
  TH1D *fPosPionsEtaSpectra;
  TH1D *fPosPionsPhiSpectra;
  TH1D *fAntiPionsPtSpectra;
  TH1D *fAntiPionsEtaSpectra;
  TH1D *fAntiPionsPhiSpectra;
  TH1D *fKaonsPtSpectra;
  TH1D *fKaonsEtaSpectra;
  TH1D *fKaonsPhiSpectra;
  TH1D *fPosKaonsPtSpectra;
  TH1D *fPosKaonsEtaSpectra;
  TH1D *fPosKaonsPhiSpectra;
  TH1D *fAntiKaonsPtSpectra;
  TH1D *fAntiKaonsEtaSpectra;
  TH1D *fAntiKaonsPhiSpectra;
  TH1D *fProtonsPtSpectra;
  TH1D *fProtonsEtaSpectra;
  TH1D *fProtonsPhiSpectra;
  TH1D *fPosProtonsPtSpectra;
  TH1D *fPosProtonsEtaSpectra;
  TH1D *fPosProtonsPhiSpectra;
  TH1D *fAntiProtonsPtSpectra;
  TH1D *fAntiProtonsEtaSpectra;
  TH1D *fAntiProtonsPhiSpectra;
  
  // Particle Spectra
  virtual void FinalizeSpectra(Int_t Nevents);
  TList *fSpectraList;
  
  TProfile *fChargedParticleSpectra;
  TH1D *fPosPionsSpectra;
  TH1D *fPosKaonsSpectra;
  TH1D *fPosProtonsSpectra;
  TH1D *fNegPionsSpectra;
  TH1D *fNegKaonsSpectra;
  TH1D *fNegProtonsSpectra;
  
  // Flow GF part
  TList *fFlowGFList;
  const static Int_t fkFlowGFNHarm = 4;
  const static Int_t fkFlowGFNOrde = 4;
  const static Int_t fFlowGFCenBin = 10;
  const static Int_t fkGFPtB = 8;
  TMatrixD *fReQGF; // fReQ[m][k] = sum_{i=1}^{M} w_{i}^{k} cos(m*phi_{i})
  TMatrixD *fImQGF; // fImQ[m][k] = sum_{i=1}^{M} w_{i}^{k} sin(m*phi_{i})
  TMatrixD *fReQGFPt[fkGFPtB]; // fReQ[m][k] = sum_{i=1}^{M} w_{i}^{k} cos(m*phi_{i})
  TMatrixD *fImQGFPt[fkGFPtB]; // fImQ[m][k] = sum_{i=1}^{M} w_{i}^{k} sin(m*phi_{i})
  
  TProfile *fFlowGFIntCorPro[fkFlowGFNHarm][fkFlowGFNOrde]; //
  TH1D *fFlowGFIntCorHist[fkFlowGFNHarm][fkFlowGFNOrde]; //
  TH1D *fFlowGFIntCumHist[fkFlowGFNHarm][fkFlowGFNOrde]; //
  TH1D *fFlowGFIntFinalHist[fkFlowGFNHarm][fkFlowGFNOrde]; //
  
  TProfile *fFlowGFIntCovPro[fkFlowGFNHarm][fkFlowGFNOrde][fkFlowGFNOrde]; //
  TH1D *fFlowGFIntCovHist[fkFlowGFNHarm][fkFlowGFNOrde][fkFlowGFNOrde]; //
  
  TProfile *fFlowGFMixedCorPro[fkFlowGFNHarm][fkFlowGFNHarm]; //
  TH1D *fFlowGFMixedCorHist[fkFlowGFNHarm][fkFlowGFNHarm]; //
  TH1D *fFlowGFMixedFinalHist[fkFlowGFNHarm][fkFlowGFNHarm]; //
  
  TProfile *fFlowGFIntCorProPtB[fkGFPtB][fkFlowGFNHarm][fkFlowGFNOrde]; //
  TH1D *fFlowGFIntCorHistPtB[fkGFPtB][fkFlowGFNHarm][fkFlowGFNOrde]; //
  TProfile *fFlowGFIntCovProPtB[fkGFPtB][fkFlowGFNHarm][fkFlowGFNOrde][fkFlowGFNOrde]; //
  TH1D *fFlowGFIntCovHistPtB[fkGFPtB][fkFlowGFNHarm][fkFlowGFNOrde][fkFlowGFNOrde]; //
  
  
  // Flow QC part
  TList *fFlowQCList;
  const static Int_t fFlowNHarm = 4;
  const static Int_t fFlowNHarmMax = 10; // WARNING: MIN (2*fFlowNHarm+2)
  const static Int_t fQVecPower = 5; //wordt dit alleen gebruikt voor de weights??
  const static Int_t fCharge = 3;
  Int_t fPtDiffNBins; //
  Int_t fEtaDiffNBins;
  Int_t fNBins;
  
  
  TH1D *fPOIPtDiffQRe[fQVecPower][fFlowNHarmMax][fCharge]; // real part
  TH1D *fPOIPtDiffQIm[fQVecPower][fFlowNHarmMax][fCharge]; // imaginary part
  TH1D *fPOIPtDiffMul[fQVecPower][fFlowNHarmMax][fCharge]; // imaginary part
  
  const static Int_t fkFlowQCnIntCorPro = 5;
  TProfile *fFlowQCIntCorPro[fFlowNHarm][fkFlowQCnIntCorPro][fCharge]; //
  TH1D *fFlowQCIntCorHist[fFlowNHarm][fkFlowQCnIntCorPro][fCharge]; //
  TH1D *fFlowQCIntFlow2Hist[fFlowNHarm][fkFlowQCnIntCorPro][fCharge];
  TH1D *fFlowQCIntFlow4Hist[fFlowNHarm][fkFlowQCnIntCorPro][fCharge];
  TH1D *fFlowQCIntCumHist[fFlowNHarm][fkFlowQCnIntCorPro][fCharge];
  
  const static Int_t fFlowQCNPro = 4;
  const static Int_t fCRCMaxnCen = 10;
  const static Int_t fFlowQCNCov = 8;
  TProfile *fFlowQCCorPro[fCRCMaxnCen][fFlowNHarm][fFlowQCNPro][fCharge];
  TProfile *fFlowQCCorCovPro[fCRCMaxnCen][fFlowNHarm][fFlowQCNCov][fCharge];
  TH1D *fFlowQCCorHist[fCRCMaxnCen][fFlowNHarm][fFlowQCNPro][fCharge]; // <<2'>>, [CRCBin][eg]
  TH1D *fFlowQCCorCovHist[fCRCMaxnCen][fFlowNHarm][fFlowQCNCov][fCharge]; // histo for covariances
  TH1D *fFlowQCFinalPtDifHist[fCRCMaxnCen][fFlowNHarm][fFlowQCNCov][fCharge]; //
  
  Int_t fFlowQCCenBin;
  
  const static Int_t fFlowQCNRef = 14;
  TProfile *fFlowQCRefCorPro[fFlowNHarm][fFlowQCNRef][fCharge]; //
  TH1D *fFlowQCRefCorHist[fFlowNHarm][fFlowQCNRef][fCharge]; //
  TH1D *fFlowQCRefCorFinal[fFlowNHarm][4][fCharge]; //
  
  //Flow EP & RP Part
  
  TH1D *fPOIDiffQRe[fFlowNHarmMax][fNParticles][fCharge]; // real part
  TH1D *fPOIDiffQIm[fFlowNHarmMax][fNParticles][fCharge]; // imaginary part
  TH1D *fPOIDiffMul[fFlowNHarmMax][fNParticles][fCharge];
  
  TH1D *fPOIDiffQRe_neg[fFlowNHarmMax][fNParticles][fCharge]; // real part
  TH1D *fPOIDiffQIm_neg[fFlowNHarmMax][fNParticles][fCharge]; // imaginary part
  TH1D *fPOIDiffMul_neg[fFlowNHarmMax][fNParticles][fCharge];
  
  TH1D *fPOIDiffQRe_pos[fFlowNHarmMax][fNParticles][fCharge]; // real part
  TH1D *fPOIDiffQIm_pos[fFlowNHarmMax][fNParticles][fCharge]; // imaginary part
  TH1D *fPOIDiffMul_pos[fFlowNHarmMax][fNParticles][fCharge];
  
  TList *fFlowEPList;

  TH1D *fEPEPresolutionPro[fFlowNHarmMax];
  TProfile *fFlowEPIntPro[fFlowNHarmMax][fNParticles][fCharge];
  TH1D *fFlowEPIntCorHist[fFlowNHarmMax][fNParticles][fCharge];
  TH1D *fFlowEPIntFlow2Hist[fFlowNHarmMax][fNParticles][fCharge];
  
  TProfile *fFlowEPCorPro[fFlowNHarmMax][fNParticles][fCharge];
  TH1D *fFlowEPDiffFlow2Hist[fFlowNHarmMax][fNParticles][fCharge];
  
  TList *fFlowRPList;

  TProfile *fFlowRPIntPro_pos[fFlowNHarmMax][fNParticles][fCharge];
  TProfile *fFlowRPIntPro_neg[fFlowNHarmMax][fNParticles][fCharge];
  TProfile *fFlowRPIntPro[fFlowNHarmMax][fNParticles][fCharge];
  
  TH1D *fFlowRPIntFlow2Hist_pos[fFlowNHarmMax][fNParticles][fCharge];
  TH1D *fFlowRPIntFlow2Hist_neg[fFlowNHarmMax][fNParticles][fCharge];
  TH1D *fFlowRPIntFlow2Hist[fFlowNHarmMax][fNParticles][fCharge];
  
  TProfile *fFlowRPCorPro[fFlowNHarmMax][fNParticles][fCharge];
  TH1D *fFlowRPDiffFlow2Hist[fFlowNHarmMax][fNParticles][fCharge];
  
  
  // Cuts:
  Double_t maxPtCut = 9999; // the data is in GeV? Double check.
  Double_t minPtCut = 0;
  Double_t maxEtaCut = 99; //|eta|<maxEtaCut
  Int_t minNtracks = 0;
  Bool_t EtaDiff= kTRUE;
  
  
  // added
  Int_t fNumberOfParticipants=0;
  Double_t fImpactParameter=-999;
  
  Double_t QRe_EP[fFlowNHarmMax];
  Double_t QIm_EP[fFlowNHarmMax];
  Double_t Mul_EP[fFlowNHarmMax];
  
};
#endif
