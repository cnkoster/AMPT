#include <fstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TMath.h"
#include "Event.h"
#include "Particle.h"
#include "CalculateFlow.h"
#include "Stopwatch_header.h"

using namespace std;
Stopwatch sw;

/// Remove Qvecpowers for SPM!! And add finalize flow!

CalculateFlow::CalculateFlow ()
{
  std::cout<<"default CalculateFlow constructor"<<'\n';
  
}

//=====================================================================================================

CalculateFlow::CalculateFlow(const char* name):fQAList(NULL),fSpectraList(NULL),fFlowQCList(NULL), fFlowSPMList(NULL), fFlowGFList(NULL), fFlowQCCenBin(10), fReQGF(NULL),fImQGF(NULL)
{
  std::cout<<"CalculateFlow constructor"<<'\n';
  
  fFlowQCList = new TList();
  fFlowQCList->SetName("fFlowQCList");
  fFlowQCList->SetOwner(kTRUE);
  
  fFlowSPMList = new TList();
  fFlowSPMList->SetName("fFlowEPMList"); //Scalar Product Method
  fFlowSPMList->SetOwner(kTRUE);
  
  fFlowGFList = new TList();
  fFlowGFList->SetName("fFlowGFList");
  fFlowGFList->SetOwner(kTRUE);
  
  fQAList = new TList();
  fQAList->SetName("fQAList");
  fQAList->SetOwner(kTRUE);
  
  // Particle Spectra
  fSpectraList = new TList();
  fSpectraList->SetName("fSpectraList");
  fSpectraList->SetOwner(kTRUE);
  
//  InitializeArraysForFlowQC();
  InitializeArraysForFlowSPM();
 // InitializeArraysForFlowGF();
  InitializeArraysForQA();
  
  for(Int_t i=0; i<fkGFPtB; i++) {
      fReQGFPt[i] = NULL;
      fImQGFPt[i] = NULL;
  }
  
}

void CalculateFlow::InitializeArraysForQA()
{
  fMultChargedParticlesDistribution = NULL;
  fNumberOfParticipantsDistribution = NULL;
  fPtChargedParticlesDistribution = NULL;
  fEtaChargedParticlesDistribution = NULL;
  fPhiChargedParticlesDistribution = NULL;
  
  fEtaPhiChargedParticlesDistribution = NULL;
  
  fPionsPtSpectra = NULL;
  fPionsEtaSpectra = NULL;
  fPionsPhiSpectra = NULL;
  fPosPionsPtSpectra = NULL;
  fPosPionsEtaSpectra = NULL;
  fPosPionsPhiSpectra = NULL;
  fAntiPionsPtSpectra = NULL;
  fAntiPionsEtaSpectra = NULL;
  fAntiPionsPhiSpectra = NULL;
  fKaonsPtSpectra = NULL;
  fKaonsEtaSpectra = NULL;
  fKaonsPhiSpectra = NULL;
  fPosKaonsPtSpectra = NULL;
  fPosKaonsEtaSpectra = NULL;
  fPosKaonsPhiSpectra = NULL;
  fAntiKaonsPtSpectra = NULL;
  fAntiKaonsEtaSpectra = NULL;
  fAntiKaonsPhiSpectra = NULL;
  fProtonsPtSpectra = NULL;
  fProtonsEtaSpectra = NULL;
  fProtonsPhiSpectra = NULL;
  fPosProtonsPtSpectra = NULL;
  fPosProtonsEtaSpectra = NULL;
  fPosProtonsPhiSpectra = NULL;
  fAntiProtonsPtSpectra = NULL;
  fAntiProtonsEtaSpectra = NULL;
  fAntiProtonsPhiSpectra = NULL;
  
  fChargedParticleSpectra = NULL;
  fPosPionsSpectra = NULL;
  fPosKaonsSpectra = NULL;
  fPosProtonsSpectra = NULL;
  fNegPionsSpectra = NULL;
  fNegKaonsSpectra = NULL;
  fNegProtonsSpectra = NULL;
}

void CalculateFlow::InitializeArraysForFlowQC()
{
  for(Int_t charge=0; charge<fCharge;charge++){
    for (Int_t c=0;c<fQVecPower;c++) {                  //fQVecPower
      for (Int_t h=0;h<fFlowNHarmMax;h++) {           //fFlowNHarmMax
        fPOIPtDiffQRe[c][h][charge] = NULL;         //POI Pt Diff Q Re [fQVecPower][fFlowHarmonic]
        fPOIPtDiffQIm[c][h][charge] = NULL;
        fPOIPtDiffMul[c][h][charge] = NULL;
      }
    }
    
    for(Int_t i=0; i<fFlowNHarm; i++) {
      for(Int_t j=0; j<fkFlowQCnIntCorPro; j++) {
        // charge=0 for all part. charge=1 for Charge+ charge=2 for Charge-
        fFlowQCIntCorPro[i][j][charge]= NULL;
        fFlowQCIntCorHist[i][j][charge] = NULL;
        fFlowQCIntFlow2Hist[i][j][charge] = NULL;
        fFlowQCIntFlow4Hist[i][j][charge] = NULL;
        fFlowQCIntCumHist[i][j][charge] = NULL;
        
      }
    }
    
    // reference flow
    for(Int_t i=0; i<fFlowNHarm; i++) {
      for(Int_t j=0; j<fFlowQCNRef; j++) {
        fFlowQCRefCorPro[i][j][charge] = NULL;
        fFlowQCRefCorHist[i][j][charge] = NULL;
      }
      for(Int_t j=0; j<4; j++) {
        fFlowQCRefCorFinal[i][j][charge] = NULL;
      }
    }
    
    // differential flow
    for (Int_t h=0; h<fCRCMaxnCen; h++) {
      for(Int_t i=0; i<fFlowNHarm; i++) {
        for(Int_t j=0; j<fFlowQCNPro; j++) {
          fFlowQCCorPro[h][i][j][charge] = NULL;
          fFlowQCCorHist[h][i][j][charge] = NULL;
        }
        for(Int_t k=0; k<fFlowQCNCov; k++) {
          fFlowQCCorCovPro[h][i][k][charge] = NULL;
          fFlowQCCorCovHist[h][i][k][charge] = NULL;
          fFlowQCFinalPtDifHist[h][i][k][charge] = NULL;
        }
      }
    }
  }
}
//=====================================================================================================

void CalculateFlow::InitializeArraysForFlowGF()
{
    for (Int_t h=0; h<fkFlowGFNHarm; h++) {
        for(Int_t i=0; i<fkFlowGFNOrde; i++) {
            fFlowGFIntCorPro[h][i] = NULL;
            fFlowGFIntCorHist[h][i] = NULL;
            fFlowGFIntCumHist[h][i] = NULL;
            fFlowGFIntFinalHist[h][i] = NULL;
            for(Int_t k=0; k<fkFlowGFNOrde; k++) {
                fFlowGFIntCovPro[h][i][k] = NULL;
                fFlowGFIntCovHist[h][i][k] = NULL;
            }
            for(Int_t s=0; s<fkGFPtB; s++) {
                fFlowGFIntCorProPtB[s][h][i] = NULL;
                fFlowGFIntCorHistPtB[s][h][i] = NULL;
                for(Int_t k=0; k<fkFlowGFNOrde; k++) {
                    fFlowGFIntCovProPtB[s][h][i][k] = NULL;
                    fFlowGFIntCovHistPtB[s][h][i][k] = NULL;
                }
            }
        }
    }
    
    for (Int_t h=0; h<fkFlowGFNHarm; h++) {
        for(Int_t i=0; i<fkFlowGFNHarm; i++) {
            fFlowGFMixedCorPro[h][i] = NULL;
            fFlowGFMixedCorHist[h][i] = NULL;
            fFlowGFMixedFinalHist[h][i] = NULL;
        }
    }
} // end of CalculateFlow::InitializeArraysForFlowGF()


//=====================================================================================================

void CalculateFlow::InitializeArraysForFlowSPM()
{
  
    for (Int_t h=0;h<fFlowNHarmMax;h++) {
      
      QRe_EP[h] = 0;                  //POI Pt Diff Q Re [fQVecPower][fFlowHarmonic]
      QIm_EP[h] = 0;
      Mul_EP[h] = 0;
      
      for(Int_t charge=0; charge<fCharge;charge++){
      //fFlowNHarmMax
    
      fFlowSPMIntPro[h][charge] = NULL;
      fFlowSPMIntFlow2Hist[h][charge] = NULL;
      
      fFlowSPM1IntPro[h][charge] = NULL;
      fFlowSPM1IntFlow2Hist[h][charge] = NULL;
      
      fPOISPMPtDiffQRe_pos[h][charge] = NULL;                  //POI Pt Diff Q Re [fQVecPower][fFlowHarmonic]
      fPOISPMPtDiffQIm_pos[h][charge] = NULL;
      fPOISPMPtDiffMul_pos[h][charge] = NULL;
      
      fPOISPMPtDiffQRe_neg[h][charge] = NULL;                  //POI Pt Diff Q Re [fQVecPower][fFlowHarmonic]
      fPOISPMPtDiffQIm_neg[h][charge] = NULL;
      fPOISPMPtDiffMul_neg[h][charge] = NULL;
      
    }
    
  }
}
//=====================================================================================================
void CalculateFlow::InitializeArraysForFlowEPM()
{
  for(Int_t charge=0; charge<fCharge;charge++){
    for (Int_t h=0;h<fFlowNHarmMax;h++) {
      
      fFlowEPMIntPro_pos[h][charge] = NULL;
      fFlowEPMIntPro_neg[h][charge] = NULL;
      
      fFlowEPMIntFlow2Hist_pos[h][charge] = NULL;
      fFlowEPMIntFlow2Hist_neg[h][charge] = NULL;
      fFlowEPMPtFlow2Hist[h][charge] = NULL;
    }
  }
}
//=====================================================================================================
void CalculateFlow::UserCreateOutputObjects() {
  
  fReQGF = new TMatrixD(21,9);
  fImQGF = new TMatrixD(21,9);
  for(Int_t i=0; i<fkGFPtB; i++) {
      fReQGFPt[i] = new TMatrixD(21,9);
      fImQGFPt[i] = new TMatrixD(21,9);
  }
  
  Double_t pTbinEdge[59] = {0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.5,5,5.5,6,6.5,7,8,9,10,11,12,13,14,15,16,18,20};
  
  
  Double_t pTPionbinEdge[42] = {0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3};
  
  Double_t pTKaonbinEdge[37] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3};
  
  Double_t pTProtonbinEdge[43] = {0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.2,3.4,3.6,3.8,4,4.2,4.4,4.6};
  
  Double_t xbins[12] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.};
  
  // QA Hists
  fMultChargedParticlesDistribution = new TH1D("fMultChargedParticlesDistribution", "fMultChargedParticlesDistribution", 200, 0, 4000);
  fNumberOfParticipantsDistribution = new TH1D("fNumberOfParticipantsDistribution", "fNumberOfParticipantsDistribution", 421, -0.5, 420.5);
  
  //fPtChargedParticlesDistribution = new TH1D("fPtChargedParticlesDistribution", "fPtChargedParticlesDistribution", 500, 0, 20);
  fPtChargedParticlesDistribution = new TH1D("fPtChargedParticlesDistribution", "fPtChargedParticlesDistribution", 58, pTbinEdge);
  fEtaChargedParticlesDistribution = new TH1D("fEtaChargedParticlesDistribution", "fEtaChargedParticlesDistribution", 200, -5, 5);
  fPhiChargedParticlesDistribution = new TH1D("fPhiChargedParticlesDistribution", "fPhiChargedParticlesDistribution", 200, -0.1, 2*TMath::Pi()+0.1);
  
  
//  fEtaPhiChargedParticlesDistribution = new TH2D("fEtaPhiChargedParticlesDistribution", "fEtaPhiChargedParticlesDistribution",200, -0.1, 2*TMath::Pi()+0.1, 200, -0.8, 0.8);
  
  // std::cout<<"Spectra is defined"<<std::endl;
  
  //fPionsPtSpectra = new TH1D("fPionsPtSpectra", "fPionsPtSpectra", 200, 0, 2);
  fPionsPtSpectra = new TH1D("fPionsPtSpectra", "fPionsPtSpectra", 41, pTPionbinEdge);
  fPionsPtSpectra->Sumw2();
  fPionsEtaSpectra = new TH1D("fPionsEtaSpectra", "fPionsEtaSpectra", 200, -5, 5);
  fPionsEtaSpectra->Sumw2();
  fPionsPhiSpectra = new TH1D("fPionsPhiSpectra", "fPionsPhiSpectra", 200, -0.1, 2*TMath::Pi()+0.1);
  fPionsPhiSpectra->Sumw2();
  
  fPosPionsPtSpectra = new TH1D("fPosPionsPtSpectra", "fPosPionsPtSpectra", 41, pTPionbinEdge);
  fPosPionsPtSpectra->Sumw2();
  fPosPionsEtaSpectra = new TH1D("fPosPionsEtaSpectra", "fPosPionsEtaSpectra", 200, -5, 5);
  fPosPionsEtaSpectra->Sumw2();
  fPosPionsPhiSpectra = new TH1D("fPosPionsPhiSpectra", "fPosPionsPhiSpectra", 200, -0.1, 2*TMath::Pi()+0.1);
  fPosPionsPhiSpectra->Sumw2();
  
  fAntiPionsPtSpectra = new TH1D("fAntiPionsPtSpectra", "fAntiPionsPtSpectra", 41, pTPionbinEdge);
  fAntiPionsPtSpectra->Sumw2();
  fAntiPionsEtaSpectra = new TH1D("fAntiPionsEtaSpectra", "fAntiPionsEtaSpectra", 200, -5, 5);
  fAntiPionsEtaSpectra->Sumw2();
  fAntiPionsPhiSpectra = new TH1D("fAntiPionsPhiSpectra", "fAntiPionsPhiSpectra", 200, -0.1, 2*TMath::Pi()+0.1);
  fAntiPionsPhiSpectra->Sumw2();
  
  //fKaonsPtSpectra = new TH1D("fKaonsPtSpectra", "fKaonsPtSpectra", 300, 0, 3);
  fKaonsPtSpectra = new TH1D("fKaonsPtSpectra", "fKaonsPtSpectra", 36, pTKaonbinEdge);
  fKaonsPtSpectra->Sumw2();
  fKaonsEtaSpectra = new TH1D("fKaonsEtaSpectra", "fKaonsEtaSpectra", 200, -5, 5);
  fKaonsEtaSpectra->Sumw2();
  fKaonsPhiSpectra = new TH1D("fKaonsPhiSpectra", "fKaonsPhiSpectra", 200, -0.1, 2*TMath::Pi()+0.1);
  fKaonsPhiSpectra->Sumw2();
  
  fPosKaonsPtSpectra = new TH1D("fPosKaonsPtSpectra", "fPosKaonsPtSpectra", 36, pTKaonbinEdge);
  fPosKaonsPtSpectra->Sumw2();
  fPosKaonsEtaSpectra = new TH1D("fPosKaonsEtaSpectra", "fPosKaonsEtaSpectra", 200, -5, 5);
  fPosKaonsEtaSpectra->Sumw2();
  fPosKaonsPhiSpectra = new TH1D("fPosKaonsPhiSpectra", "fPosKaonsPhiSpectra", 200, -0.1, 2*TMath::Pi()+0.1);
  fPosKaonsPhiSpectra->Sumw2();
  
  fAntiKaonsPtSpectra = new TH1D("fAntiKaonsPtSpectra", "fAntiKaonsPtSpectra", 36, pTKaonbinEdge);
  fAntiKaonsPtSpectra->Sumw2();
  fAntiKaonsEtaSpectra = new TH1D("fAntiKaonsEtaSpectra", "fAntiKaonsEtaSpectra", 200, -5, 5);
  fAntiKaonsEtaSpectra->Sumw2();
  fAntiKaonsPhiSpectra = new TH1D("fAntiKaonsPhiSpectra", "fAntiKaonsPhiSpectra", 200, -0.1, 2*TMath::Pi()+0.1);
  fAntiKaonsPhiSpectra->Sumw2();
  
  //fProtonsPtSpectra = new TH1D("fProtonsPtSpectra", "fProtonsPtSpectra", 400, 0, 4);
  fProtonsPtSpectra = new TH1D("fProtonsPtSpectra", "fProtonsPtSpectra", 42, pTProtonbinEdge);
  fProtonsPtSpectra->Sumw2();
  fProtonsEtaSpectra = new TH1D("fProtonsEtaSpectra", "fProtonsEtaSpectra", 200, -5, 5);
  fProtonsEtaSpectra->Sumw2();
  fProtonsPhiSpectra = new TH1D("fProtonsPhiSpectra", "fProtonsPhiSpectra", 200, -0.1, 2*TMath::Pi()+0.1);
  fProtonsPhiSpectra->Sumw2();
  
  fPosProtonsPtSpectra = new TH1D("fPosProtonsPtSpectra", "fPosProtonsPtSpectra", 42, pTProtonbinEdge);
  fPosProtonsPtSpectra->Sumw2();
  fPosProtonsEtaSpectra = new TH1D("fPosProtonsEtaSpectra", "fPosProtonsEtaSpectra", 200, -5, 5);
  fPosProtonsEtaSpectra->Sumw2();
  fPosProtonsPhiSpectra = new TH1D("fPosProtonsPhiSpectra", "fPosProtonsPhiSpectra", 200, -0.1, 2*TMath::Pi()+0.1);
  fPosProtonsPhiSpectra->Sumw2();
  
  fAntiProtonsPtSpectra = new TH1D("fAntiProtonsPtSpectra", "fAntiProtonsPtSpectra", 42, pTProtonbinEdge);
  fAntiProtonsPtSpectra->Sumw2();
  fAntiProtonsEtaSpectra = new TH1D("fAntiProtonsEtaSpectra", "fAntiProtonsEtaSpectra", 200, -5, 5);
  fAntiProtonsEtaSpectra->Sumw2();
  fAntiProtonsPhiSpectra = new TH1D("fAntiProtonsPhiSpectra", "fAntiProtonsPhiSpectra", 200, -0.1, 2*TMath::Pi()+0.1);
  fAntiProtonsPhiSpectra->Sumw2();
  
  fQAList->Add(fMultChargedParticlesDistribution);
  fQAList->Add(fNumberOfParticipantsDistribution);
  
  fQAList->Add(fPtChargedParticlesDistribution);
  fQAList->Add(fEtaChargedParticlesDistribution);
  fQAList->Add(fPhiChargedParticlesDistribution);
  //fQAList->Add(fEtaPhiChargedParticlesDistribution);
  //std::cout<<"Spectra is added"<<std::endl;
  fQAList->Add(fPionsPtSpectra);
  fQAList->Add(fPionsEtaSpectra);
  fQAList->Add(fPionsPhiSpectra);
  fQAList->Add(fPosPionsPtSpectra);
  fQAList->Add(fPosPionsEtaSpectra);
  fQAList->Add(fPosPionsPhiSpectra);
  fQAList->Add(fAntiPionsPtSpectra);
  fQAList->Add(fAntiPionsEtaSpectra);
  fQAList->Add(fAntiPionsPhiSpectra);
  fQAList->Add(fKaonsPtSpectra);
  fQAList->Add(fKaonsEtaSpectra);
  fQAList->Add(fKaonsPhiSpectra);
  fQAList->Add(fPosKaonsPtSpectra);
  fQAList->Add(fPosKaonsEtaSpectra);
  fQAList->Add(fPosKaonsPhiSpectra);
  fQAList->Add(fAntiKaonsPtSpectra);
  fQAList->Add(fAntiKaonsEtaSpectra);
  fQAList->Add(fAntiKaonsPhiSpectra);
  fQAList->Add(fProtonsPtSpectra);
  fQAList->Add(fProtonsEtaSpectra);
  fQAList->Add(fProtonsPhiSpectra);
  fQAList->Add(fPosProtonsPtSpectra);
  fQAList->Add(fPosProtonsEtaSpectra);
  fQAList->Add(fPosProtonsPhiSpectra);
  fQAList->Add(fAntiProtonsPtSpectra);
  fQAList->Add(fAntiProtonsEtaSpectra);
  fQAList->Add(fAntiProtonsPhiSpectra);
  
  // Particle spectra
  fChargedParticleSpectra = new TProfile("fChargedParticleSpectra","fChargedParticleSpectra",11, xbins,"s");
  fPosPionsSpectra = new TH1D("fPosPionsSpectra", "fPosPionsSpectra", 58, pTbinEdge);
  fPosKaonsSpectra = new TH1D("fPosKaonsSpectra", "fPosKaonsSpectra", 58, pTbinEdge);
  fPosProtonsSpectra = new TH1D("fPosProtonsSpectra", "fPosProtonsSpectra", 58, pTbinEdge);
  fNegPionsSpectra = new TH1D("fNegPionsSpectra", "fNegPionsSpectra", 58, pTbinEdge);
  fNegKaonsSpectra = new TH1D("fNegKaonsSpectra", "fNegKaonsSpectra", 58, pTbinEdge);
  fNegProtonsSpectra = new TH1D("fNegProtonsSpectra", "fNegProtonsSpectra", 58, pTbinEdge);
  
  fSpectraList->Add(fChargedParticleSpectra);
  fSpectraList->Add(fPosPionsSpectra);
  fSpectraList->Add(fPosKaonsSpectra);
  fSpectraList->Add(fPosProtonsSpectra);
  fSpectraList->Add(fNegPionsSpectra);
  fSpectraList->Add(fNegKaonsSpectra);
  fSpectraList->Add(fNegProtonsSpectra);
  
  // Calculate FlowQC Hists
  // choose for eta diff or pt diff:
  
  fPtDiffNBins = 36-8; //for pt 20 for eta
  fEtaDiffNBins = 5;
  
  fCRCPtBins = new Double_t[37-8];
  Double_t PtBins[] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.25,1.5,1.75,2.,2.25,2.5,2.75,3.,3.25,3.5,3.75,4.,4.5,5.,5.5,6.,7.,8.,9.,10.};//,12.,14.,17.,20.,25.,30.,40.,50.};
  Double_t EtaBins[] = {-0.8,-0.48,-0.16,0.16,0.48,0.8};
  Double_t ImPaBins[] = {3.72, 5.23, 7.31, 8.88, 10.20, 11.38, 12.47, 13.50, 14.51, 15.0};
  
  
  if(EtaDiff){
    std::cout<< "We have eta bins" <<std::endl;
    fNBins = fEtaDiffNBins;
    fBins = new Double_t[fNBins+1];
    for(Int_t r=0; r<9; r++) {
      fBins[r] = EtaBins[r];
    }
  }
  
  if(!EtaDiff){
    std::cout<< "We have pt bins" <<std::endl;
    fNBins = fPtDiffNBins;
    fBins = new Double_t[fNBins+1];
    for(Int_t r=0; r<37-8; r++) {
      fBins[r] = PtBins[r];
    }
  }
  
  // so for each power of the Q-vector and each flowHarmonic make a Th1D
  // In output vinden we: [0][0] tot [13][13] ?
  
  
  for (Int_t c=0;c<fQVecPower;c++) {
    for (Int_t h=0;h<fFlowNHarmMax;h++) {
      for(Int_t charge=0; charge<fCharge; charge++){
        fPOIPtDiffQRe[c][h][charge] = new TH1D(Form("fPOIPtDiffQRe[%d][%d][%d]",c,h,charge),Form("fPOIPtDiffQRe[%d][%d][%d]",c,h,charge), fNBins, fBins);
        //   fFlowQCList->Add(fPOIPtDiffQRe[c][h]);
        fPOIPtDiffQIm[c][h][charge] = new TH1D(Form("fPOIPtDiffQIm[%d][%d][%d]",c,h,charge),Form("fPOIPtDiffQIm[%d][%d][%d]",c,h,charge), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
        //  fFlowQCList->Add(fPOIPtDiffQIm[c][h]);
        fPOIPtDiffMul[c][h][charge] = new TH1D(Form("fPOIPtDiffMul[%d][%d][%d]",c,h,charge),Form("fPOIPtDiffMul[%d][%d][%d]",c,h,charge), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
        //   fFlowQCList->Add(fPOIPtDiffMul[c][h]);
      }
    }
  }
  
  
  for(Int_t i=0; i<fFlowNHarm; i++) {
    for(Int_t j=0; j<fkFlowQCnIntCorPro; j++) {
      for(Int_t charge=0; charge<fCharge; charge++){
        fFlowQCIntCorPro[i][j][charge] = new TProfile(Form("fFlowQCIntCorPro[%d][%d][%d]",i,j,charge),Form("fFlowQCIntCorPro[%d][%d][%d]",i,j,charge),9,ImPaBins,"s"); //here we changed the bins to nParticipants
        fFlowQCIntCorPro[i][j][charge]->Sumw2();
        //fFlowQCList->Add(fFlowQCIntCorPro[i][j]);
        fFlowQCIntCorHist[i][j][charge] = new TH1D(Form("fFlowQCIntCorHist[%d][%d][%d]",i,j,charge),Form("fFlowQCIntCorHist[%d][%d][%d]",i,j,charge),9,ImPaBins);
        fFlowQCIntCorHist[i][j][charge]->Sumw2();
        //   fFlowQCList->Add(fFlowQCIntCorHist[i][j][charge]);
        
        fFlowQCIntFlow2Hist[i][j][charge] = new TH1D(Form("fFlowQCIntFlow2Hist[%d][%d][%d]",i,j,charge),Form("fFlowQCIntFlow2Hist[%d][%d][%d]",i,j,charge),9,ImPaBins);
        fFlowQCIntFlow2Hist[i][j][charge]->Sumw2();
        
        
        fFlowQCIntFlow4Hist[i][j][charge] = new
        TH1D(Form("fFlowQCIntFlow4Hist[%d][%d][%d]",i,j,charge),Form("fFlowQCIntFlow4Hist[%d][%d][%d]",i,j,charge),9,ImPaBins);
        fFlowQCIntFlow4Hist[i][j][charge]->Sumw2();
        
        if(j==0){
          fFlowQCList->Add(fFlowQCIntFlow2Hist[i][j][charge]);}
          //        fFlowQCList->Add(fFlowQCIntFlow4Hist[i][j][charge]);}
          
          fFlowQCIntCumHist[i][j][charge] = new TH1D(Form("fFlowQCIntCumHist[%d][%d][%d]",i,j,charge),Form("fFlowQCIntCumHist[%d][%d][%d]",i,j,charge),9,ImPaBins);
          fFlowQCIntCumHist[i][j][charge]->Sumw2();
          //   fFlowQCList->Add(fFlowQCIntCumHist[i][j][charge]);
        }
      }
    }
  
  
  // reference flow
  for(Int_t i=0; i<fFlowNHarm; i++) {
    for(Int_t j=0; j<fFlowQCNRef; j++) {
      for(Int_t charge=0; charge<fCharge; charge++){
        fFlowQCRefCorPro[i][j][charge] = new TProfile(Form("fFlowQCRefCorPro[%d][%d][%d]",i,j,charge),Form("fFlowQCRefCorPro[%d][%d][%d]",i,j,charge),9,ImPaBins,"s");
        fFlowQCRefCorPro[i][j][charge]->Sumw2();
        //fFlowQCList->Add(fFlowQCRefCorPro[i][j][charge]);
        fFlowQCRefCorHist[i][j][charge] = new TH1D(Form("fFlowQCRefCorHist[%d][%d][%d]",i,j,charge),Form("fFlowQCRefCorHist[%d][%d][%d]",i,j,charge),9,ImPaBins);
        fFlowQCRefCorHist[i][j][charge]->Sumw2();
        // fFlowQCList->Add(fFlowQCRefCorHist[i][j][charge]);
      }
    }
    for(Int_t j=0; j<4; j++) {
      for(Int_t charge=0; charge<fCharge; charge++){
        fFlowQCRefCorFinal[i][j][charge] = new TH1D(Form("fFlowQCRefCorFinal[%d][%d][%d]",i,j,charge),Form("fFlowQCRefCorFinal[%d][%d][%d]",i,j,charge),9,ImPaBins);
        fFlowQCRefCorFinal[i][j][charge]->Sumw2();
        //  fFlowQCList->Add(fFlowQCRefCorFinal[i][j][charge]);
      }
    }
  }
  
  // differential flow
  for (Int_t h=0; h<fCRCnCen; h++) {
    for(Int_t i=0; i<fFlowNHarm; i++) {
      for(Int_t j=0; j<fFlowQCNPro; j++) {
        for(Int_t charge=0; charge<fCharge; charge++){
          fFlowQCCorPro[h][i][j][charge] = new TProfile(Form("fFlowQCCorPro[%d][%d][%d][%d]",h,i,j,charge),Form("fFlowQCCorPro[%d][%d][%d][%d]",h,i,j,charge), fNBins, fBins,"s");
          fFlowQCCorPro[h][i][j][charge]->Sumw2();
          //      fFlowQCList->Add(fFlowQCCorPro[h][i][j]);
          fFlowQCCorHist[h][i][j][charge] = new TH1D(Form("fFlowQCCorHist[%d][%d][%d][%d]",h,i,j,charge),Form("fFlowQCCorHist[%d][%d][%d][%d]",h,i,j,charge), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
          fFlowQCCorHist[h][i][j][charge]->Sumw2();
          //  fFlowQCList->Add(fFlowQCCorHist[h][i][j]);
        }
      }
      
      for(Int_t k=0; k<fFlowQCNCov; k++) {
        for(Int_t charge=0; charge<fCharge; charge++){
          fFlowQCCorCovPro[h][i][k][charge] = new TProfile(Form("fFlowQCCorCovPro[%d][%d][%d][%d]",h,i,k,charge),Form("fFlowQCCorCovPro[%d][%d][%d][%d]",h,i,k,charge), fNBins, fBins,"s");
          fFlowQCCorCovPro[h][i][k][charge]->Sumw2();
          //  fFlowQCList->Add(fFlowQCCorCovPro[h][i][k]);
          fFlowQCCorCovHist[h][i][k][charge] = new TH1D(Form("fFlowQCCorCovHist[%d][%d][%d][%d]",h,i,k,charge),Form("fFlowQCCorCovHist[%d][%d][%d][%d]",h,i,k,charge), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
          fFlowQCCorCovHist[h][i][k][charge]->Sumw2();
          //    fFlowQCList->Add(fFlowQCCorCovHist[h][i][k]);
          fFlowQCFinalPtDifHist[h][i][k][charge] = new TH1D(Form("fFlowQCFinalPtDifHist[%d][%d][%d][%d]",h,i,k,charge),Form("fFlowQCFinalPtDifHist[%d][%d][%d][%d]",h,i,k,charge), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
          fFlowQCFinalPtDifHist[h][i][k][charge]->Sumw2();
          
          if(h==GetCRCCenBin(fCentralityEBE)){
            fFlowQCList->Add(fFlowQCFinalPtDifHist[GetCRCCenBin(fCentralityEBE)][i][k][charge]);} //we say cen=1 so the histo is filled only for h=6
        }
      }
    }
  }
  
  // for CalculateFlowGF
  for (Int_t h=0; h<fkFlowGFNHarm; h++) {
      for(Int_t i=0; i<fkFlowGFNOrde; i++) {
          fFlowGFIntCorPro[h][i] = new TProfile(Form("fFlowGFIntCorPro[%d][%d]",h,i),Form("fFlowGFIntCorPro[%d][%d]",h,i),9,ImPaBins,"s");//fFlowGFCenBin,0.,100.,"s");
          fFlowGFIntCorPro[h][i]->Sumw2();
          //      fFlowGFList->Add(fFlowGFIntCorPro[h][i]);
          fFlowGFIntCorHist[h][i] = new TH1D(Form("fFlowGFIntCorHist[%d][%d]",h,i),Form("fFlowGFIntCorHist[%d][%d]",h,i),9,ImPaBins);//fFlowGFCenBin,0.,100.);
          fFlowGFIntCorHist[h][i]->Sumw2();
          //     fFlowGFList->Add(fFlowGFIntCorHist[h][i]);
          fFlowGFIntCumHist[h][i] = new TH1D(Form("fFlowGFIntCumHist[%d][%d]",h,i),Form("fFlowGFIntCumHist[%d][%d]",h,i),9,ImPaBins);//fFlowGFCenBin,0.,100.);
          fFlowGFIntCumHist[h][i]->Sumw2();
          //         fFlowGFList->Add(fFlowGFIntCumHist[h][i]);
          fFlowGFIntFinalHist[h][i] = new TH1D(Form("fFlowGFIntFinalHist[%d][%d]",h,i),Form("fFlowGFIntFinalHist[%d][%d]",h,i),9,ImPaBins);//fFlowGFCenBin,0.,100.);
          fFlowGFIntFinalHist[h][i]->Sumw2();
          fFlowGFList->Add(fFlowGFIntFinalHist[h][i]);
          for(Int_t k=0; k<fkFlowGFNOrde; k++) {
              fFlowGFIntCovPro[h][i][k] = new TProfile(Form("fFlowGFIntCovPro[%d][%d][%d]",h,i,k),Form("fFlowGFIntCovPro[%d][%d][%d]",h,i,k),9,ImPaBins,"s");//fFlowGFCenBin,0.,100.,"s");
              fFlowGFIntCovPro[h][i][k]->Sumw2();
              //           fFlowGFList->Add(fFlowGFIntCovPro[h][i][k]);
              fFlowGFIntCovHist[h][i][k] = new TH1D(Form("fFlowGFIntCovHist[%d][%d][%d]",h,i,k),Form("fFlowGFIntCovHist[%d][%d][%d]",h,i,k),9,ImPaBins);//fFlowGFCenBin,0.,100.);
              fFlowGFIntCovHist[h][i][k]->Sumw2();
              //          fFlowGFList->Add(fFlowGFIntCovHist[h][i][k]);
          }
          for(Int_t s=0; s<fkGFPtB; s++) {
              fFlowGFIntCorProPtB[s][h][i] = new TProfile(Form("fFlowGFIntCorProPtB[%d][%d][%d]",s,h,i),Form("fFlowGFIntCorProPtB[%d][%d][%d]",s,h,i),9,ImPaBins,"s");//fFlowGFCenBin,0.,100.,"s");
              fFlowGFIntCorProPtB[s][h][i]->Sumw2();
              //            fFlowGFList->Add(fFlowGFIntCorProPtB[s][h][i]);
              fFlowGFIntCorHistPtB[s][h][i] = new TH1D(Form("fFlowGFIntCorHistPtB[%d][%d][%d]",s,h,i),Form("fFlowGFIntCorHistPtB[%d][%d][%d]",s,h,i),9,ImPaBins);//fFlowGFCenBin,0.,100.);
              fFlowGFIntCorHistPtB[s][h][i]->Sumw2();
              //            fFlowGFList->Add(fFlowGFIntCorHistPtB[s][h][i]);
              for(Int_t k=0; k<fkFlowGFNOrde; k++) {
                  fFlowGFIntCovProPtB[s][h][i][k] = new TProfile(Form("fFlowGFIntCovProPtB[%d][%d][%d][%d]",s,h,i,k),Form("fFlowGFIntCovProPtB[%d][%d][%d][%d]",s,h,i,k),9,ImPaBins,"s");//fFlowGFCenBin,0.,100.,"s");
                  fFlowGFIntCovProPtB[s][h][i][k]->Sumw2();
                  //  fFlowGFList->Add(fFlowGFIntCovProPtB[s][h][i][k]);
                  fFlowGFIntCovHistPtB[s][h][i][k] = new TH1D(Form("fFlowGFIntCovHistPtB[%d][%d][%d][%d]",s,h,i,k),Form("fFlowGFIntCovHistPtB[%d][%d][%d][%d]",s,h,i,k),9,ImPaBins);//fFlowGFCenBin,0.,100.);
                  fFlowGFIntCovHistPtB[s][h][i][k]->Sumw2();
                  //fFlowGFList->Add(fFlowGFIntCovHistPtB[s][h][i][k]);
              }
          }
      }
  }
  
  for (Int_t h=0; h<fkFlowGFNHarm; h++) {
      for(Int_t i=0; i<fkFlowGFNHarm; i++) {
          fFlowGFMixedCorPro[h][i] = new TProfile(Form("fFlowGFMixedCorPro[%d][%d]",h,i),Form("fFlowGFMixedCorPro[%d][%d]",h,i),9,ImPaBins);//fFlowGFCenBin,0.,100.,"s");
          fFlowGFMixedCorPro[h][i]->Sumw2();
          //  fFlowGFList->Add(fFlowGFMixedCorPro[h][i]);
          fFlowGFMixedCorHist[h][i] = new TH1D(Form("fFlowGFMixedCorHist[%d][%d]",h,i),Form("fFlowGFMixedCorHist[%d][%d]",h,i),9,ImPaBins);//fFlowGFCenBin,0.,100.);
          fFlowGFMixedCorHist[h][i]->Sumw2();
          //   fFlowGFList->Add(fFlowGFMixedCorHist[h][i]);
          fFlowGFMixedFinalHist[h][i] = new TH1D(Form("fFlowGFMixedFinalHist[%d][%d]",h,i),Form("fFlowGFMixedFinalHist[%d][%d]",h,i),9,ImPaBins);//fFlowGFCenBin,0.,100.);
          fFlowGFMixedFinalHist[h][i]->Sumw2();
          //   fFlowGFList->Add(fFlowGFMixedFinalHist[h][i]);
      }
  }
  
  // Calculate FlowSPM Hists --------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // so for each power of the Q-vector and each flowHarmonic make a Th1D
  // In output vinden we: [0][0] tot [13][13] ?
  
  
  
  for (Int_t h=0;h<fFlowNHarm;h++) {
    fSPMEPresolutionPro[h] = new TH1D(Form("fSPMEPresolutionPro[%d]",h),Form("fSPMEPresolutionPro[%d]",h),1000,-TMath::Pi()-0.2,TMath::Pi()+0.2);
    fSPMEPresolutionPro[h]->Sumw2();
//    fFlowSPMList->Add(fSPMEPresolutionPro[h]);
    
    for(Int_t charge=0; charge<fCharge; charge++){
      
      fPOISPMPtDiffQRe_pos[h][charge] = new TH1D(Form("fPOISPMPtDiffQRe_pos[%d][%d]",h,charge),Form("fPOISPMPtDiffQRe_pos[%d][%d]",h,charge), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
      fPOISPMPtDiffQRe_pos[h][charge]->Sumw2();
      fPOISPMPtDiffQIm_pos[h][charge] = new TH1D(Form("fPOISPMPtDiffQIm_pos[%d][%d]",h,charge),Form("fPOISPMPtDiffQIm_pos[%d][%d]",h,charge), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
      fPOISPMPtDiffQIm_pos[h][charge]->Sumw2();
      fPOISPMPtDiffMul_pos[h][charge] = new TH1D(Form("fPOISPMPtDiffMul_pos[%d][%d]",h,charge),Form("fPOISPMPtDiffMul_pos[%d][%d]",h,charge), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
      fPOISPMPtDiffMul_pos[h][charge]->Sumw2();
      
      
      fPOISPMPtDiffQRe_neg[h][charge] = new TH1D(Form("fPOISPMPtDiffQRe_neg[%d][%d]",h,charge),Form("fPOISPMPtDiffQRe_neg[%d][%d]",h,charge), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
      fPOISPMPtDiffQRe_neg[h][charge]->Sumw2();
      fPOISPMPtDiffQIm_neg[h][charge] = new TH1D(Form("fPOISPMPtDiffQIm_neg[%d][%d]",h,charge),Form("fPOISPMPtDiffQIm_neg[%d][%d]",h,charge), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
      fPOISPMPtDiffQIm_neg[h][charge]->Sumw2();
      fPOISPMPtDiffMul_neg[h][charge] = new TH1D(Form("fPOISPMPtDiffMul_neg[%d][%d]",h,charge),Form("fPOISPMPtDiffMul_neg[%d][%d]",h,charge), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
      fPOISPMPtDiffMul_neg[h][charge]->Sumw2();
      
      
      
      fRFPSPMPtDiffQRe_V0A[h][charge] = new TH1D(Form("fRFPSPMPtDiffQRe_V0A[%d][%d]",h,charge),Form("fRFPSPMPtDiffQRe_V0A[%d][%d]",h,charge),fNBins,-10,10);//,fPtDiffNBins,fCRCPtBins);
      fRFPSPMPtDiffQRe_V0A[h][charge]->Sumw2();
      fRFPSPMPtDiffQIm_V0A[h][charge] = new TH1D(Form("fRFPSPMPtDiffQIm_V0A[%d][%d]",h,charge),Form("fRFPSPMPtDiffQIm_V0A[%d][%d]",h,charge),fNBins,-10,10);//,fPtDiffNBins,fCRCPtBins);
      fRFPSPMPtDiffQIm_V0A[h][charge]->Sumw2();
      fRFPSPMPtDiffMul_V0A[h][charge] = new TH1D(Form("fRFPSPMPtDiffMul_V0A[%d][%d]",h,charge),Form("fRFPSPMPtDiffMul_V0A[%d][%d]",h,charge),fNBins,-10,10);//,fPtDiffNBins,fCRCPtBins);
      fRFPSPMPtDiffMul_V0A[h][charge]->Sumw2();
      fRFPSPMPtDiffQRe_V0C[h][charge] = new TH1D(Form("fRFPSPMPtDiffQRe_V0C[%d][%d]",h,charge),Form("fRFPSPMPtDiffQRe_V0C[%d][%d]",h,charge),fNBins,-10,10);//,fPtDiffNBins,fCRCPtBins);
      fRFPSPMPtDiffQRe_V0C[h][charge]->Sumw2();
      fRFPSPMPtDiffQIm_V0C[h][charge] = new TH1D(Form("fRFPSPMPtDiffQIm_V0C[%d][%d]",h,charge),Form("fRFPSPMPtDiffQIm_V0C[%d][%d]",h,charge),fNBins,-10,10);//,fPtDiffNBins,fCRCPtBins);
      fRFPSPMPtDiffQIm_V0C[h][charge]->Sumw2();
      fRFPSPMPtDiffMul_V0C[h][charge] = new TH1D(Form("fRFPSPMPtDiffMul_V0C[%d][%d]",h,charge),Form("fRFPSPMPtDiffMul_V0C[%d][%d]",h,charge),fNBins,-10,10);//,fPtDiffNBins,fCRCPtBins);
      fRFPSPMPtDiffMul_V0C[h][charge]->Sumw2();
      
      fFlowSPMIntPro[h][charge]= new TProfile(Form("fFlowSPMIntPro[%d][%d]",h,charge),Form("fFlowSPMIntPro[%d][%d]",h,charge),9,ImPaBins,"s");
      fFlowSPMIntPro[h][charge]->Sumw2();

      
      fFlowSPMIntFlow2Hist[h][charge] = new TH1D(Form("fFlowSPMIntFlow2Hist[%d][%d]",h,charge),Form("fFlowSPMIntFlow2Hist[%d][%d]",h,charge),9,ImPaBins);
      fFlowSPMIntFlow2Hist[h][charge]->Sumw2();
      fFlowSPMList->Add(fFlowSPMIntFlow2Hist[h][charge]);
      
      fFlowSPM1IntPro[h][charge]= new TProfile(Form("fSPM1EPresolutionPro[%d][%d]",h,charge),Form("fSPM1EPresolutionPro%d][%d]",h,charge),9,ImPaBins,"s");
      fFlowSPM1IntPro[h][charge]->Sumw2();
      
      fFlowSPM1IntFlow2Hist[h][charge] = new TH1D(Form("fFlowSPM1IntFlow2Hist[%d][%d]",h,charge),Form("fFlowSPM1IntFlow2Hist[%d][%d]",h,charge),9,ImPaBins);
      fFlowSPM1IntFlow2Hist[h][charge]->Sumw2();
      //            fFlowSPMList->Add(fFlowSPM1IntFlow2Hist[h][charge]);
      
    }
  }
  
  
  for (Int_t h=0;h<fFlowNHarm;h++) {
    for(Int_t charge=0; charge<fCharge; charge++){
      // Event Plane Method (SPECTATORS +EPangle=0)
      
      fFlowEPMIntPro_pos[h][charge]= new TProfile(Form("fFlowEPMIntPro_pos[%d][%d]",h,charge),Form("fFlowEPMIntPro_pos[%d][%d]",h,charge),9,ImPaBins,"s");
      fFlowEPMIntPro_pos[h][charge]->Sumw2();
      
      fFlowEPMIntPro_neg[h][charge]= new TProfile(Form("fFlowEPMIntPro_neg[%d][%d]",h,charge),Form("fFlowEPMIntPro_neg[%d][%d]",h,charge),9,ImPaBins,"s");
      fFlowEPMIntPro_neg[h][charge]->Sumw2();
      
      fFlowEPMIntPro[h][charge]= new TProfile(Form("fFlowEPMIntPro[%d][%d]",h,charge),Form("fFlowEPMIntPro[%d][%d]",h,charge),9,ImPaBins,"s");
      fFlowEPMIntPro[h][charge]->Sumw2();
      
      fFlowEPMIntFlow2Hist_pos[h][charge] = new TH1D(Form("fFlowEPMIntFlow2Hist_pos[%d][%d]",h,charge),Form("fFlowEPMIntFlow2Hist_pos[%d][%d]",h,charge),9,ImPaBins); //,fPtDiffNBins,fCRCPtBins);
      fFlowEPMIntFlow2Hist_pos[h][charge]->Sumw2();
//      fFlowSPMList->Add(fFlowEPMIntFlow2Hist_pos[h][charge]);
      
      fFlowEPMIntFlow2Hist_neg[h][charge] = new TH1D(Form("fFlowEPMIntFlow2Hist_neg[%d][%d]",h,charge),Form("fFlowEPMIntFlow2Hist_neg[%d][%d]",h,charge),9,ImPaBins); //,fPtDiffNBins,fCRCPtBins);
      fFlowEPMIntFlow2Hist_neg[h][charge]->Sumw2();
//      fFlowSPMList->Add(fFlowEPMIntFlow2Hist_neg[h][charge]);
      
      fFlowEPMIntFlow2Hist[h][charge] = new TH1D(Form("fFlowEPMIntFlow2Hist[%d][%d]",h,charge),Form("fFlowEPMIntFlow2Hist[%d][%d]",h,charge),9,ImPaBins); //,fPtDiffNBins,fCRCPtBins);
      fFlowEPMIntFlow2Hist[h][charge]->Sumw2();
//      fFlowSPMList->Add(fFlowEPMIntFlow2Hist[h][charge]);
      
      fFlowEPMCorPro[h][charge]= new TProfile(Form("fFlowEPMCorPro[%d][%d]",h,charge),Form("fFlowEPMCorPro[%d][%d]",h,charge), fNBins, fBins, "s");//, fNBins, fBins,"s");
      fFlowEPMCorPro[h][charge]->Sumw2();
      
      fFlowEPMPtFlow2Hist[h][charge] = new TH1D(Form("fFlowEPMPtFlow2Hist[%d][%d]",h,charge),Form("fFlowEPMPtFlow2Hist[%d][%d]",h,charge), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
      fFlowEPMPtFlow2Hist[h][charge]->Sumw2();
      fFlowSPMList->Add(fFlowEPMPtFlow2Hist[h][charge]);
      
      
    }
  }
  
  
  
  
}


// Double_t centRange[10] = {0,10,20,30,40,50,60,70,80,90};


//=====================================================================================================

void CalculateFlow::UserExec()
{
  Make(fEvent);
}

//=====================================================================================================

void CalculateFlow::Make(Event* anEvent) {
  
  sw.tick(); //start timer
  Int_t nPrim = 0;
  Int_t Pid = -99999;
  Double_t dPhi = 0.; // azimuthal angle in the laboratory frame
  Double_t dPt  = 0.; // transverse momentum
  Double_t dEta = 0.; // pseudorapidity
  Int_t dCharge = -2; // charge
  Double_t totaltime=0;
  fNumberOfParticipants = 0;
  Stopwatch sw1;
  Double_t fQ2TPC = 0, fQ2V0A = 0., fQ2V0C = 0., fQ2 = 0.;
  
  fCenBin = GetCRCCenBin(fCentralityEBE);
  nPrim = anEvent->getNtrack();
  //std::cout<<"nPrim is: " <<nPrim<< "   "<< "minTraks is: "<< minNtracks<<std::endl;
  fNumberOfParticipants = anEvent->getnPart();
  fImpactParameter = anEvent->getb();
  // std::cout<<fImpactParameter<<std::endl;
  Int_t cw = 0;
  Int_t nChargedParticles = 0;
  // multiplicity for charged particles
  Double_t fNumOfPos = 0;
  Double_t fNumOfNeg = 0;
  Double_t xval;
  

  
  if (nPrim < minNtracks) return;
  sw.tick();
  for(Int_t i=0;i<nPrim;i++) {
    if(anEvent->getParticle(i).getSpecflag()<0) continue;
    Pid = anEvent->getParticle(i).getPid();
    dPhi = anEvent->getParticle(i).getPhi();
    dPt = anEvent->getParticle(i).getPt();
    dEta = anEvent->getParticle(i).getRapidity();  // in MC, the rapidity is stored rather than pesudorapidity (eta)
    dCharge = anEvent->getParticle(i).getCharge();
    
    // std::cout<<"before cuts "<<anEvent->getParticle(i).getSpecflag()<<std::endl
    
    
    if (dPt > maxPtCut) continue;
    if (dPt < minPtCut) continue;
    
    for (Int_t h=0;h<3;h++) {
      QRe_EP[h] += pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi);
      QIm_EP[h] += pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi);
      Mul_EP[h] += pow(wPhiEta,1);
    }
    


    
    if (dCharge == 0) continue;
    cw = (dCharge > 0. ? 0 : 1);
  
    
    if(EtaDiff) xval = dEta;
    if(!EtaDiff) xval = dPt;
    
    if (TMath::Abs(dEta) > maxEtaCut) continue;
    
    nChargedParticles++; //inside |eta|<0.8
    
    // ====== Plot some QA ============
    // spectra for pions, kaons and protons
    if (doQA) {
      fPtChargedParticlesDistribution->Fill(dPt);
      fEtaChargedParticlesDistribution->Fill(dEta);
      fPhiChargedParticlesDistribution->Fill(dPhi);
      //      fEtaPhiChargedParticlesDistribution->Fill(dPhi,dEta);
      
      // Pions
      if (Pid == 211 || Pid == -211) {
        fPionsPtSpectra->Fill(dPt);
        fPionsEtaSpectra->Fill(dEta);
        fPionsPhiSpectra->Fill(dPhi);
        if (Pid == 211) {
          fPosPionsPtSpectra->Fill(dPt);
          fPosPionsEtaSpectra->Fill(dEta);
          fPosPionsPhiSpectra->Fill(dPhi);
          fPosPionsSpectra->Fill(dPt);
        }
        if (Pid == -211) {
          fAntiPionsPtSpectra->Fill(dPt);
          fAntiPionsEtaSpectra->Fill(dEta);
          fAntiPionsPhiSpectra->Fill(dPhi);
          fNegPionsSpectra->Fill(dPt);
        }
      }
      // Kaons
      if (Pid == 321 || Pid == -321) {
        fKaonsPtSpectra->Fill(dPt);
        fKaonsEtaSpectra->Fill(dEta);
        fKaonsPhiSpectra->Fill(dPhi);
        if (Pid == 321) {
          fPosKaonsPtSpectra->Fill(dPt);
          fPosKaonsEtaSpectra->Fill(dEta);
          fPosKaonsPhiSpectra->Fill(dPhi);
          fPosKaonsSpectra->Fill(dPt);
        }
        if (Pid == -321) {
          fAntiKaonsPtSpectra->Fill(dPt);
          fAntiKaonsEtaSpectra->Fill(dEta);
          fAntiKaonsPhiSpectra->Fill(dPhi);
          fNegKaonsSpectra->Fill(dPt);
        }
      }
      // Protons
      if (Pid == 2212 || Pid == -2212) {
        fProtonsPtSpectra->Fill(dPt);
        fProtonsEtaSpectra->Fill(dEta);
        fProtonsPhiSpectra->Fill(dPhi);
        if (Pid == 2212) {
          fPosProtonsPtSpectra->Fill(dPt);
          fPosProtonsEtaSpectra->Fill(dEta);
          fPosProtonsPhiSpectra->Fill(dPhi);
          fPosProtonsSpectra->Fill(dPt);
        }
        if (Pid == -2212) {
          fAntiProtonsPtSpectra->Fill(dPt);
          fAntiProtonsEtaSpectra->Fill(dEta);
          fAntiProtonsPhiSpectra->Fill(dPhi);
          fNegProtonsSpectra->Fill(dPt);
        }
      }
    }
    
    
    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // if (dEta > 0) continue; //Determine v1 at negative eta only!
    // if (dEta < 0) continue; //Determine v1 at positive eta only!
    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
    // ====== for calculateFlowGF (generic framework) =====
    // Generic Framework: Calculate Re[Q_{m*n,k}] and Im[Q_{m*n,k}] for this event (m = 1,2,...,12, k = 0,1,...,8):
    //    for(Int_t m=0;m<21;m++) {
    //      for(Int_t k=0;k<9;k++) {
    //  (*fReQGF)(m,k) += pow(wPhiEta*wPhi*wPt*wEta*wTrack,k)*TMath::Cos(m*dPhi);
    //  (*fImQGF)(m,k) += pow(wPhiEta*wPhi*wPt*wEta*wTrack,k)*TMath::Sin(m*dPhi);
    //      }
    //    }
    //
    //    for(Int_t ptb=0; ptb<fkGFPtB; ptb++) {
    //      if(ptb==0 && dPt>0.5) continue;
    //      if(ptb==1 && (dPt<0.5 || dPt>1.)) continue;
    //      if(ptb==2 && (dPt<1. || dPt>2.)) continue;
    //      if(ptb==3 && dPt<2.) continue;
    //      if(ptb==4 && (dPt<1. || dPt>2.5)) continue;
    //      if(ptb==5 && dPt<2.5) continue;
    //      if(ptb==6 && (dPt<1. || dPt>3.)) continue;
    //      if(ptb==7 && dPt<3.) continue;
    //      for(Int_t m=0;m<21;m++) {
    //  for(Int_t k=0;k<9;k++) {
    //    (*fReQGFPt[ptb])(m,k) += pow(wPhiEta*wPhi*wPt*wEta*wTrack,k)*TMath::Cos(m*dPhi);
    //    (*fImQGFPt[ptb])(m,k) += pow(wPhiEta*wPhi*wPt*wEta*wTrack,k)*TMath::Sin(m*dPhi);
    //  }
    //      }
    //    }
    //
    // ====== for calculateFlowQC =========
    // Also POI SPM!
    
    
        for (Int_t k=0; k<fQVecPower; k++) {
          for (Int_t h=0;h<fFlowNHarmMax;h++) {
            fPOIPtDiffQRe[k][h][0]->Fill(xval,pow(wPhiEta,k)*TMath::Cos((h+1.)*dPhi));
            fPOIPtDiffQIm[k][h][0]->Fill(xval,pow(wPhiEta,k)*TMath::Sin((h+1.)*dPhi));
            fPOIPtDiffMul[k][h][0]->Fill(xval,pow(wPhiEta,k));
    
    
    
            if(dCharge>0){
              //if(dPhi>TMath::Pi)
              fPOIPtDiffQRe[k][h][1]->Fill(xval,pow(wPhiEta,k)*TMath::Cos((h+1.)*dPhi));
              fPOIPtDiffQIm[k][h][1]->Fill(xval,pow(wPhiEta,k)*TMath::Sin((h+1.)*dPhi));
              fPOIPtDiffMul[k][h][1]->Fill(xval,pow(wPhiEta,k));
    
            }
    
    
            if(dCharge<0){
              fPOIPtDiffQRe[k][h][2]->Fill(xval,pow(wPhiEta,k)*TMath::Cos((h+1.)*dPhi));
              fPOIPtDiffQIm[k][h][2]->Fill(xval,pow(wPhiEta,k)*TMath::Sin((h+1.)*dPhi));
              fPOIPtDiffMul[k][h][2]->Fill(xval,pow(wPhiEta,k));
            }
    
          }
        }
    
    if(dEta>=0){
      for (Int_t h=0;h<3;h++) {
        fPOISPMPtDiffQRe_pos[h][0]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
        fPOISPMPtDiffQIm_pos[h][0]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
        fPOISPMPtDiffMul_pos[h][0]->Fill(xval,pow(wPhiEta,1));
        
        
        
        if(dCharge>0){
          //if(dPhi>TMath::Pi)
          fPOISPMPtDiffQRe_pos[h][1]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
          fPOISPMPtDiffQIm_pos[h][1]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
          fPOISPMPtDiffMul_pos[h][1]->Fill(xval,pow(wPhiEta,1));
          
        }
        
        
        if(dCharge<0){
          fPOISPMPtDiffQRe_pos[h][2]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
          fPOISPMPtDiffQIm_pos[h][2]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
          fPOISPMPtDiffMul_pos[h][2]->Fill(xval,pow(wPhiEta,1));}
        
      }
    }
    
        if(dEta<0){
          for (Int_t h=0;h<3;h++) {
            fPOISPMPtDiffQRe_neg[h][0]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
            fPOISPMPtDiffQIm_neg[h][0]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
            fPOISPMPtDiffMul_neg[h][0]->Fill(xval,pow(wPhiEta,1));
    
    
    
            if(dCharge>0){
              //if(dPhi>TMath::Pi)
              fPOISPMPtDiffQRe_neg[h][1]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
              fPOISPMPtDiffQIm_neg[h][1]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
              fPOISPMPtDiffMul_neg[h][1]->Fill(xval,pow(wPhiEta,1));
    
            }
    
    
            if(dCharge<0){
              fPOISPMPtDiffQRe_neg[h][2]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
              fPOISPMPtDiffQIm_neg[h][2]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
              fPOISPMPtDiffMul_neg[h][2]->Fill(xval,pow(wPhiEta,1));
    
            }
    
          }
        }
    
  }
  
  if (doQA) {
    fMultChargedParticlesDistribution->Fill(nChargedParticles);
    fNumberOfParticipantsDistribution->Fill(fNumberOfParticipants);
    fChargedParticleSpectra->Fill(fCentralityEBE, nChargedParticles/(2*maxEtaCut));
  }
  sw.tock(); //stop timer
  
  // cout << "Time for looping over event: "<< sw.takeTime() << endl;
  // cout <<"Generic Framework matrix loop: " << totaltime << endl;
  
  //sw.tick();
  //CalculateFlowQC();
 // CalculateFlowGF();
//  CalculateFlowSPM();
  //    CalculateFlowSPM1();
  CalculateFlowEPM();
  
  // sw.tock();
  // cout << "Time for CalculateFlowCQ() "<< sw.takeTime() << endl;
  ResetEventByEventQuantities();
}


//=====================================================================================================

void CalculateFlow::ResetEventByEventQuantities()
{
  //FlowGF
  fReQGF->Zero();
  fImQGF->Zero();
  for(Int_t i=0; i<fkGFPtB; i++) {
    fReQGFPt[i]->Zero();
    fImQGFPt[i]->Zero();
  }
  
  //FlowQC
  for(Int_t c=0;c<fQVecPower;c++) {
    for (Int_t h=0;h<fFlowNHarmMax;h++) {
      for (Int_t charge=0; charge<fCharge; charge++){
        if(fPOIPtDiffQRe[c][h][charge]) fPOIPtDiffQRe[c][h][charge]->Reset();
        if(fPOIPtDiffQIm[c][h][charge]) fPOIPtDiffQIm[c][h][charge]->Reset();
        if(fPOIPtDiffMul[c][h][charge]) fPOIPtDiffMul[c][h][charge]->Reset();
      }
    }
  }
  
  //FlowSPM & EPM
  for (Int_t h=0;h<fFlowNHarm;h++) {
    
    if(QRe_EP[h]!=0) QRe_EP[h]=0;
    if(QIm_EP[h]!=0) QIm_EP[h]=0;
    if(Mul_EP[h]!=0) Mul_EP[h]=0;

    
    for (Int_t charge=0; charge<fCharge; charge++){
    
      if(fPOISPMPtDiffQRe_pos[h][charge]) fPOISPMPtDiffQRe_pos[h][charge]->Reset();
      if(fPOISPMPtDiffQIm_pos[h][charge]) fPOISPMPtDiffQIm_pos[h][charge]->Reset();
      if(fPOISPMPtDiffMul_pos[h][charge]) fPOISPMPtDiffMul_pos[h][charge]->Reset();
      
      if(fPOISPMPtDiffQRe_neg[h][charge]) fPOISPMPtDiffQRe_neg[h][charge]->Reset();
      if(fPOISPMPtDiffQIm_neg[h][charge]) fPOISPMPtDiffQIm_neg[h][charge]->Reset();
      if(fPOISPMPtDiffMul_neg[h][charge]) fPOISPMPtDiffMul_neg[h][charge]->Reset();
      
      
      
    }
  }
  
  
}

//=====================================================================================================

void CalculateFlow::Terminate(Int_t Nevents)
{
  FinalizeSpectra(Nevents);
 // FinalizeFlowQC();
//  FinalizeFlowGF();
//  FinalizeFlowSPM();
  FinalizeFlowEPM();
  FinalizeQA();
}

void CalculateFlow::CalculateFlowGF()
{
    Int_t order[fkFlowGFNOrde] = {0};
    for(Int_t k=0; k<fkFlowGFNOrde; k++) {
        order[k] = 2*(k+1);
    }
    
    Double_t dMult = (*fReQGF)(0,0);
    
    for(Int_t hr=0; hr<fkFlowGFNHarm; hr++) {
        
        Double_t CorrOrd[fkFlowGFNOrde] = {0.};
        Double_t WeigOrd[fkFlowGFNOrde] = {0.};
        
        for(Int_t no=0; no<fkFlowGFNOrde; no++) {
            
            if(dMult<order[no]) continue;
            
            //here change h+2 to h+1
            TArrayI harmonics(order[no]);
            Int_t halforder = (Int_t)order[no]/2;
            for(Int_t k=0; k<order[no]; k++) {
                if(k<halforder) harmonics[k] = -(hr+1);
                else            harmonics[k] = hr+1;
            }
            TArrayI emptyness(order[no]);
            for(Int_t k=0; k<order[no]; k++) emptyness[k] = 0;
            
            std::complex<double> N = this->ucN(order[no], harmonics, -1);
            std::complex<double> D = this->ucN(order[no], emptyness, -1);
            
            if(D.real()>0.) {
                fFlowGFIntCorPro[hr][no]->Fill(fImpactParameter,N.real()/D.real(),D.real()*fCenWeightEbE);
                CorrOrd[no] = N.real()/D.real();
                WeigOrd[no] = D.real();
            }
            
        } // end of for(Int_t no=0; no<fkFlowGFNOrde; no++)
        
        for(Int_t no=0; no<fkFlowGFNOrde; no++) {
            for(Int_t no2=0; no2<fkFlowGFNOrde; no2++) {
                if(WeigOrd[no]>0. && WeigOrd[no2]>0.) {
                    fFlowGFIntCovPro[hr][no][no2]->Fill(fImpactParameter,CorrOrd[no]*CorrOrd[no2],WeigOrd[no]*WeigOrd[no2]*fCenWeightEbE*fCenWeightEbE);
                }
            }
        }
        
    } // end of for(Int_t hr=0; hr<fkFlowGFNHarm; hr++)
    
    for(Int_t hr=0; hr<fkFlowGFNHarm; hr++) {
        for(Int_t hr2=0; hr2<fkFlowGFNHarm; hr2++) {
            
            if(dMult<4) continue;
            //here change h+2 to h+1
            TArrayI harmonics(4);
            harmonics[0] = hr+1;
            harmonics[1] = hr2+1;
            harmonics[2] = -(hr+1);
            harmonics[3] = -(hr2+1);
            
            TArrayI emptyness(4);
            for(Int_t k=0; k<4; k++) emptyness[k] = 0;
            
            std::complex<double> N = this->ucN(4, harmonics, -1);
            std::complex<double> D = this->ucN(4, emptyness, -1);
            
            if(D.real()>0.) {
                fFlowGFMixedCorPro[hr][hr2]->Fill(fImpactParameter,N.real()/D.real(),D.real()*fCenWeightEbE);
            }
        }
    }
    
    // in wide pt bins
    for(Int_t i=0; i<fkGFPtB; i++) {
        Double_t dMult = (*fReQGFPt[i])(0,0);
        
        for(Int_t hr=0; hr<fkFlowGFNHarm; hr++) {
            
            Double_t CorrOrd[fkFlowGFNOrde] = {0.};
            Double_t WeigOrd[fkFlowGFNOrde] = {0.};
            
            for(Int_t no=0; no<fkFlowGFNOrde; no++) {
                
                if(dMult<order[no]) continue;
                //here change h+2 to h+1
                TArrayI harmonics(order[no]);
                Int_t halforder = (Int_t)order[no]/2;
                for(Int_t k=0; k<order[no]; k++) {
                    if(k<halforder) harmonics[k] = -(hr+1);
                    else            harmonics[k] = hr+1;
                }
                TArrayI emptyness(order[no]);
                for(Int_t k=0; k<order[no]; k++) emptyness[k] = 0;
                
                std::complex<double> N = this->ucN(order[no], harmonics, i);
                std::complex<double> D = this->ucN(order[no], emptyness, i);
                
                if(D.real()>0.) {
                    fFlowGFIntCorProPtB[i][hr][no]->Fill(fImpactParameter,N.real()/D.real(),D.real()*fCenWeightEbE);
                    CorrOrd[no] = N.real()/D.real();
                    WeigOrd[no] = D.real();
                }
                
            } // end of for(Int_t no=0; no<fkFlowGFNOrde; no++)
            
            for(Int_t no=0; no<fkFlowGFNOrde; no++) {
                for(Int_t no2=0; no2<fkFlowGFNOrde; no2++) {
                    if(WeigOrd[no]>0. && WeigOrd[no2]>0.) {
                        fFlowGFIntCovProPtB[i][hr][no][no2]->Fill(fImpactParameter,CorrOrd[no]*CorrOrd[no2],WeigOrd[no]*WeigOrd[no2]*fCenWeightEbE*fCenWeightEbE);
                    }
                }
            }
            
        } // end of for(Int_t hr=0; hr<fkFlowGFNHarm; hr++)
    } // end of for(Int_t i=0; i<fkGFPtB; i++)
    
} // end of CalculateFlow::CalculateFlowGF()



//=======================================================================================================================


void CalculateFlow::CalculateFlowSPM()
{
  

  Float_t QRe, QIm, Mu;
  Float_t Denom_pty;
  Float_t x_QQ, y_QQ, x_uQ, y_uQ;
  
  Float_t EventPlane;
  
  Float_t cosEvPl, sinEvPl;
  Float_t cosPhi, sinPhi;
  Float_t EP_res;
  
  //********************************************************************
  // Calculate Q vectors for the spectators (P&T)
  // Then determine correlations with the POI (u)
  // *******************************************************************
  
    
  for (Int_t hr=0;hr<fFlowNHarm;hr++) {
    
    for(Int_t charge=0;charge<fCharge;charge++){
      // ********************************************************************
      // pT-integrated: {2} *************************************************
      // ********************************************************************

      QRe=0.; QIm=0.;
      Mu = 0.;
      
      
      Denom_pty=0.;
      
      EventPlane = TMath::ATan2(QIm_EP[hr], QRe_EP[hr]);

      // std::cout<< " EP_res " << EP_res << std::endl;
      
      cosEvPl = TMath::Cos(EventPlane);
      sinEvPl = TMath::Sin(EventPlane);
      
      for(Int_t pt=0; pt<fNBins; pt++) {
        //std::cout<<fPOIPtDiffQRe[1][hr][charge]->GetBinContent(pt+1)<<std::endl;
        QRe += fPOISPMPtDiffQRe_pos[hr][charge]->GetBinContent(pt+1);
        QIm += fPOISPMPtDiffQIm_pos[hr][charge]->GetBinContent(pt+1);
        Mu += fPOISPMPtDiffMul_pos[0][charge]->GetBinContent(pt+1);
      }
      
      
      x_QQ = QRe*cosEvPl + QIm*sinEvPl;
      //            std::cout<< " ------------------ " << std::endl;
      //
      //            std::cout<< QRe << "  " << QIm << "  " << Mu << std::endl;
      //            std::cout<< QRe_RFP_V0A << "  " << QRe_RFP_V0C << "  " << QIm_RFP_V0A << std::endl;
      
      //            x_QQ = ( (QRe_RFP_V0A*QRe) * (QRe_RFP_V0C*QRe_RFP_V0A) )/(QRe*QRe_RFP_V0C);
      //            y_QQ = -1*( (QIm_RFP_V0A*QIm) * (QIm_RFP_V0C*QIm_RFP_V0A) )/(QIm*QIm_RFP_V0C);
      //            x_uQ = (QRe*QRe_RFP_V0A);
      //            y_uQ = -1*(QIm*QIm_RFP_V0A);
      
//
//      x_QQ = ( (QRe_RFP_V0A*QRe)/Mu * (QRe_RFP_V0C*QRe)/Mu )/(QRe_RFP_V0A*QRe_RFP_V0C);
//      y_QQ = ( (QIm_RFP_V0A*QIm)/Mu * (QIm_RFP_V0C*QIm)/Mu )/(QIm_RFP_V0A*QIm_RFP_V0C);
//      x_uQ = (QRe*QRe)/Mu;
//      y_uQ = (QIm*QIm)/Mu;
      
      
      if(TMath::Abs(Mu)>0.){
        Denom_pty = (x_QQ) / Mu; }
      
      //            std::cout<<Denom_pty<<std::endl;
      
      fFlowSPMIntPro[hr][charge]->Fill(fImpactParameter, Denom_pty,1.); //1 for weights
      
      
      
    } //end of for(Int_t charge=0;charge<fCharge;charge++)
    fSPMEPresolutionPro[hr]->Fill(EventPlane);
  }//end of nHarm
}



//=====================================================================================================
void CalculateFlow::CalculateFlowSPM1()
{
  // verwacht dat deze niet werkt? Use q-vectors in V0A and V0C acceptance
  // From the directed flow paper.
  
  
  Float_t QRe_RFP_V0A, QIm_RFP_V0A, Mu_RFP_V0A;
  Float_t QRe_RFP_V0C, QIm_RFP_V0C, Mu_RFP_V0C;
  Float_t QRe, QIm, Mu;
  Float_t Denom_pty, Denom_ptx;
  Float_t x_QQ, y_QQ, x_uQ_V0A, y_uQ_V0A,x_uQ_V0C, y_uQ_V0C;
  Float_t v_V0A=0., v_V0C=0.;
  
  Float_t EvevntPlane, EventPlaneV0A, EventPlaneV0C;
  
  Float_t cosEvPlV0A, sinEvPlV0A;
  Float_t cosEvPlV0C, sinEvPlV0C;
  Float_t cosPhi, sinPhi;
  Float_t EP_res;
  
  //********************************************************************
  // Calculate Q vectors for the spectators (P&T)
  // Then determine correlations with the POI (u)
  // *******************************************************************
  for (Int_t hr=0;hr<fFlowNHarm;hr++) {
    
    for(Int_t charge=0;charge<fCharge;charge++){
      // ********************************************************************
      // pT-integrated: {2} *************************************************
      // ********************************************************************
      
      
      QRe_RFP_V0A=0.; QIm_RFP_V0A=0.;
      Mu_RFP_V0A = 0.;
      
      QRe_RFP_V0C=0.; QIm_RFP_V0C=0.;
      Mu_RFP_V0C = 0.;
      
      QRe=0.; QIm=0.;
      Mu = 0.;
      
      
      Denom_pty=0.; Denom_ptx=0.;
      
      
      
      for(Int_t pt=0; pt<fNBins; pt++) {
        
        QRe_RFP_V0A += fRFPSPMPtDiffQRe_V0A[hr][charge]->GetBinContent(pt+1); // Cos((hr+1.)*dPhi) --> Re(u)
        QIm_RFP_V0A += fRFPSPMPtDiffQIm_V0A[hr][charge]->GetBinContent(pt+1); // Sin((hr+1.)*dPhi) --> Im(u)
        Mu_RFP_V0A +=  fRFPSPMPtDiffMul_V0A[0][charge]->GetBinContent(pt+1);
        
        QRe_RFP_V0C += fRFPSPMPtDiffQRe_V0C[hr][charge]->GetBinContent(pt+1); // Cos((hr+1.)*dPhi) --> Re(u)
        QIm_RFP_V0C += fRFPSPMPtDiffQIm_V0C[hr][charge]->GetBinContent(pt+1); // Sin((hr+1.)*dPhi) --> Im(u)
        Mu_RFP_V0C +=  fRFPSPMPtDiffMul_V0C[0][charge]->GetBinContent(pt+1);
        //  std::cout<<fPOIPtDiffQRe[1][hr][charge]->GetBinCenter(pt+1)<<std::endl;
        
        QRe += fPOIPtDiffQRe[1][hr][charge]->GetBinContent(pt+1);
        QIm += fPOIPtDiffQIm[1][hr][charge]->GetBinContent(pt+1);
        Mu += fPOIPtDiffMul[1][0][charge]->GetBinContent(pt+1);}
      
      x_QQ = TMath::Sqrt(TMath::Abs(QIm_RFP_V0A * QIm_RFP_V0C));
      y_QQ = TMath::Sqrt(TMath::Abs(QIm_RFP_V0A * QIm_RFP_V0C));
      x_uQ_V0A = (QRe*QRe_RFP_V0A)/Mu;
      x_uQ_V0C = (QRe*QRe_RFP_V0C)/Mu;
      y_uQ_V0A = (QIm*QIm_RFP_V0A)/Mu;
      y_uQ_V0C = (QIm*QIm_RFP_V0C)/Mu;
      
      
      
      if(TMath::Abs(Mu)>0. && TMath::Abs(x_QQ) > 0 && TMath::Abs(y_QQ) > 0){
        v_V0A = (1/TMath::Sqrt(2) * (x_uQ_V0A/x_QQ + y_uQ_V0A/y_QQ));
        v_V0C = (-1/TMath::Sqrt(2) * (x_uQ_V0C/x_QQ + y_uQ_V0C/y_QQ));
        
        fFlowSPM1IntPro[hr][charge]->Fill(fImpactParameter, (v_V0A + v_V0C)/2., 1.); //1 for weights
      }
    } //end of for(Int_t charge=0;charge<fCharge;charge++)
  }//end of nHarm
}

//=======================================================================================================================
void CalculateFlow::CalculateFlowEPM()
{
  
  Float_t QRe, QIm, Mu;
  Float_t QRe_pos, QIm_pos, Mu_pos, QRe_neg, QIm_neg, Mu_neg;
  Float_t qpRe, qpIm, qpM, meanPtdiff;
  Double_t FillPtBin = 0;
  Float_t v_pos, v_neg, v_tot;
  
  //!!!!!!!!!!!!!!!!!! NOTE: What we do here only makes sense for odd v2!!
  //!
  //********************************************************************
  // Calculate Q vectors for the spectators (P&T)
  // Then determine correlations with the POI (u)
  // *******************************************************************
  for (Int_t h=0;h<fFlowNHarm;h++) {
    
    for(Int_t charge=0;charge<fCharge;charge++){
      // ********************************************************************
      // pT-integrated ******************************************************
      // ********************************************************************
      
      QRe=0.; QIm=0.;
      Mu = 0.;
      
      QRe_pos =0; QIm_pos=0; Mu_pos=0;
      QRe_neg =0; QIm_neg=0; Mu_neg=0;
      
      for(Int_t pt=0; pt<fNBins; pt++) {
        QRe_pos += fPOISPMPtDiffQRe_pos[h][charge]->GetBinContent(pt+1);
        QIm_pos += fPOISPMPtDiffQIm_pos[h][charge]->GetBinContent(pt+1);
        Mu_pos += fPOISPMPtDiffMul_pos[h][charge]->GetBinContent(pt+1);
        
        
        //Only for eta<0
        QRe_neg += fPOISPMPtDiffQRe_neg[h][charge]->GetBinContent(pt+1);
        QIm_neg += fPOISPMPtDiffQIm_neg[h][charge]->GetBinContent(pt+1);
        Mu_neg += fPOISPMPtDiffMul_neg[h][charge]->GetBinContent(pt+1);
        
        
      }
      
      
      if(TMath::Abs(Mu_neg)>0.){
        v_neg = QRe_neg/Mu_neg;
        fFlowEPMIntPro_neg[h][charge]->Fill(fImpactParameter, v_neg , 1.);
        fFlowEPMIntPro[h][charge]->Fill(fImpactParameter, (QRe_neg+QRe_pos)/(Mu_pos+Mu_neg),1.);
      }
      
      if(TMath::Abs(Mu_pos)>0.) {
        v_pos = QRe_pos/Mu_pos;
        fFlowEPMIntPro_pos[h][charge]->Fill(fImpactParameter, v_pos , 1.); //1 for weights
      }
      
      
      // store pt-differential flow ******************************************************************************************************************************************
      for(Int_t pt=0; pt<fNBins; pt++) {
        
        FillPtBin = fPOIPtDiffQRe[1][h][charge]->GetBinCenter(pt+1);
        qpRe=0.; qpIm=0.; qpM=0.;
        qpRe = fPOIPtDiffQRe[1][h][charge]->GetBinContent(pt+1);
        qpIm = fPOIPtDiffQIm[1][h][charge]->GetBinContent(pt+1);
        qpM = fPOIPtDiffMul[1][h][charge]->GetBinContent(pt+1);
        
        
        
        if(qpM>0) {
          //          if(h==0){
          //            if(FillPtBin < 0 ) {meanPtdiff = ( (qpRe*QReInt)/qpM ) / (TMath::Abs(QReInt));}
          //            if(FillPtBin >= 0 ) {meanPtdiff = ( (qpRe*QRe)/qpM ) / (TMath::Abs(QRe));}
          //          }
          
          //if(h>=0){
          meanPtdiff = qpRe/qpM;//(( (qpRe*(QRe_pos+QRe_neg))/qpM ) / (TMath::Abs(QRe_pos+QRe_neg)) ) + (( -1*(qpIm*(QIm_pos+QIm_neg))/qpM ) / (TMath::Abs(QIm_pos+QIm_neg)) ) ;//}
          //  std::cout<<meanPtdiff<< " " << FillPtBin <<std::endl;
          fFlowEPMCorPro[h][charge]->Fill(FillPtBin, meanPtdiff, 1.);            // ADD: fPOIEPMPtDiffQRe[h]
        }
        
      }//end pt
      
    } //end charge loop
  }//end harm loop
  
  
}

void CalculateFlow::CalculateFlowQC()
{
  Double_t FillPtBin = 0.;
  Double_t IQC2[fFlowNHarm] = {0.};
  Double_t IQC4[fFlowNHarm] = {0.};
  Double_t IQM2=0., IQM4=0.;
  Double_t WQM2=0., WQM4=0.;
  Double_t dQC2=0., dQC4=0.;
  Double_t dQM2=0., dQM4=0.;
  Double_t WdQM2=0., WdQM4=0., WdQM2EG=0., WdQM2EGB=0.;
  Double_t QRe=0., QIm=0., Q2Re2=0., Q2Im2=0., QRe3=0., QIm3=0., QM0=0., QM=0., QM2=0., QM3=0., QM4=0.;
  Double_t qpRe0=0., qpIm0=0., qpRe2=0., qpIm2=0., qp2Re=0., qp2Im=0., qpM0=0., qpM=0., qpM2=0., qpM3=0.;
  Double_t WqpM0=0., WqpAM=0.;
  Bool_t Q2f=kFALSE, Q4f=kFALSE, dQ2f=kFALSE, dQ4f=kFALSE;
  Bool_t WeigMul = kTRUE; //(fCorrWeightTPC==kMultiplicity ? kTRUE : kFALSE);

  for(Int_t charge=0;charge<fCharge;charge++){
    //********************************************************************
    // 0= inclusive 1=positive 2= negative *******************************
    // *******************************************************************


    for(Int_t hr=0; hr<fFlowNHarm; hr++) {
      // ********************************************************************
      // pT-integrated: {2}, {4} ********************************************
      // ********************************************************************

      // store reference flow (2 and 4p) ***********************************
      QRe=0.; QIm=0.; Q2Re2=0.; Q2Im2=0.; QRe3=0.; QIm3=0.;
      QM0=0.; QM=0.; QM2=0.; QM3=0.; QM4=0.;
      Q2f=kFALSE; Q4f=kFALSE;

      //here change h+2 to h+1 (hr+1->hr, 2hr+3-> 2hr+1)
      for(Int_t pt=0; pt<fNBins; pt++) {
        QRe += fPOIPtDiffQRe[1][hr][charge]->GetBinContent(pt+1); // Cos((hr+1.)*dPhi)
        QIm += fPOIPtDiffQIm[1][hr][charge]->GetBinContent(pt+1); // Sin((hr+1.)*dPhi)
        Q2Re2 += fPOIPtDiffQRe[2][2*hr+1][charge]->GetBinContent(pt+1); // w^2*Cos((hr+1.)*dPhi)
        Q2Im2 += fPOIPtDiffQIm[2][2*hr+1][charge]->GetBinContent(pt+1); // w^2*Sin((hr+1.)*dPhi)
        QRe3 += fPOIPtDiffQRe[3][hr][charge]->GetBinContent(pt+1); // w^3*Cos((hr+1.)*dPhi)
        QIm3 += fPOIPtDiffQIm[3][hr][charge]->GetBinContent(pt+1); // w^3*Sin((hr+1.)*dPhi)

        QM0 += fPOIPtDiffMul[0][0][charge]->GetBinContent(pt+1); // w^0
        QM  += fPOIPtDiffMul[1][0][charge]->GetBinContent(pt+1); // w^1
        QM2 += fPOIPtDiffMul[2][0][charge]->GetBinContent(pt+1); // w^2
        QM3 += fPOIPtDiffMul[3][0][charge]->GetBinContent(pt+1); // w^3
        QM4 += fPOIPtDiffMul[4][0][charge]->GetBinContent(pt+1); // w^4
      }

      IQM2 = QM*QM-QM2;
      WQM2 = (WeigMul? IQM2 : 1.);
      if(QM0>1) {
        IQC2[hr] = (QRe*QRe+QIm*QIm-QM2)/IQM2; // <2> = |Qn,1|^2-S1,2/(S2,1-S1,2)
        fFlowQCIntCorPro[hr][0][charge]->Fill(fImpactParameter,IQC2[hr],WQM2*fCenWeightEbE); // nPart vs. <2>
        fFlowQCRefCorPro[hr][0][charge]->Fill(fImpactParameter,IQC2[hr],WQM2*fCenWeightEbE); // nPart vs. <2>
        Q2f = kTRUE;
      }

      IQM4 = QM*QM*QM*QM - 6.*QM2*QM*QM + 8.*QM3*QM + 3.*QM2*QM2 - 6.*QM4;
      WQM4 = (WeigMul? IQM4 : 1.);
      if(QM0>3) {
        IQC4[hr] = ((QRe*QRe+QIm*QIm)*(QRe*QRe+QIm*QIm)                     // |Q_n,1|^4
                    - 2.*(QRe*QRe*Q2Re2+2.*QRe*QIm*Q2Im2-QIm*QIm*Q2Re2)     //
                    + 8.*(QRe3*QRe+QIm3*QIm)
                    + (Q2Re2*Q2Re2+Q2Im2*Q2Im2)
                    - 4.*QM2*(QRe*QRe+QIm*QIm)
                    - 6.*QM4 + 2.*QM2*QM2) / IQM4;

        fFlowQCIntCorPro[hr][1][charge]->Fill(fImpactParameter,IQC4[hr],WQM4*fCenWeightEbE);
        fFlowQCRefCorPro[hr][1][charge]->Fill(fImpactParameter,IQC4[hr],WQM4*fCenWeightEbE);
        Q4f = kTRUE;
      }

      // product of correlations or covariances
      if(Q2f && Q4f) {
        fFlowQCIntCorPro[hr][2][charge]->Fill(fImpactParameter,IQC2[hr]*IQC4[hr],WQM2*WQM4*fCenWeightEbE);
        fFlowQCRefCorPro[hr][13][charge]->Fill(fImpactParameter,IQC2[hr]*IQC4[hr],WQM2*WQM4*fCenWeightEbE);
      }


      // ********************************************************************
      // pT-differential: {2}, {4} ******************************************
      // ********************************************************************

      // store pt-differential flow ****************************************
      for(Int_t pt=0; pt<fNBins; pt++) {

        FillPtBin = fPOIPtDiffQRe[1][hr][charge]->GetBinCenter(pt+1);
        
        qpRe0=0.; qpIm0=0.; qpRe2=0.; qpIm2=0.; qp2Re=0.; qp2Im=0.; qpM0=0.; qpM=0.; qpM2=0.; qpM3=0.;
        
        //here change h+2 to h+1 (hr+1->hr, 2hr+3-> 2hr+1)
        //std::cout<<fPOIPtDiffQRe[0][hr][charge]->GetNbinsX()<<std::endl;
        
        qpRe0 = fPOIPtDiffQRe[0][hr][charge]->GetBinContent(pt+1);
        qpIm0 = fPOIPtDiffQIm[0][hr][charge]->GetBinContent(pt+1);
        qpRe2 = fPOIPtDiffQRe[2][hr][charge]->GetBinContent(pt+1);
        qpIm2 = fPOIPtDiffQIm[2][hr][charge]->GetBinContent(pt+1);
        qp2Re = fPOIPtDiffQRe[1][2*hr+1][charge]->GetBinContent(pt+1);
        qp2Im = fPOIPtDiffQIm[1][2*hr+1][charge]->GetBinContent(pt+1);

        qpM0 = fPOIPtDiffMul[0][0][charge]->GetBinContent(pt+1);
        qpM  = fPOIPtDiffMul[1][0][charge]->GetBinContent(pt+1);
        qpM2 = fPOIPtDiffMul[2][0][charge]->GetBinContent(pt+1);
        qpM3 = fPOIPtDiffMul[3][0][charge]->GetBinContent(pt+1);

        //if(hr==0) {//////////////////////////////////////////////////////////////////////////////////
        // fFlowFlowectra->Fill(fImpactParameter,FillPtBin,qpM*fCenWeightEbE);
        //}

        dQM2 = qpM0*QM-qpM;
        WdQM2 = (WeigMul? dQM2 : 1.);

        if(qpM0>0 && QM0>0) {
          dQC2 = (qpRe0*QRe+qpIm0*QIm-qpM)/dQM2;
         // std::cout<<"pt bin: " << FillPtBin<< " filled with : " << dQC2 << std::endl;
          fFlowQCCorPro[fCenBin][hr][1][charge]->Fill(FillPtBin,dQC2,WdQM2*fCenWeightEbE);
          dQ2f = kTRUE;
        }

        dQM4 = qpM0*(QM*QM*QM-3.*QM*QM2+2.*QM3)-3.*(qpM*(QM*QM-QM2)+2.*(qpM3-qpM2*QM));
        WdQM4 = (WeigMul? dQM4 : 1.);
        //@Shi dQM4 has to be nonzero
        if(qpM0>0 && QM0>3 && dQM4!=0) {
          dQC4 = ((pow(QRe,2.)+pow(QIm,2.))*(qpRe0*QRe+qpIm0*QIm)
                  - qp2Re*(pow(QRe,2.)-pow(QIm,2.))
                  - 2.*qp2Im*QRe*QIm
                  - qpRe0*(QRe*Q2Re2+QIm*Q2Im2)
                  + qpIm0*(QIm*Q2Re2-QRe*Q2Im2)
                  - 2.*QM2*(qpRe0*QRe+qpIm0*QIm)
                  - 2.*(pow(QRe,2.)+pow(QIm,2.))*qpM
                  + 6.*(qpRe2*QRe+qpIm2*QIm)
                  + 1.*(qp2Re*Q2Re2+qp2Im*Q2Im2)
                  + 2.*(qpRe0*QRe3+qpIm0*QIm3)
                  + 2.*qpM*QM2
                  - 6.*qpM3) / dQM4;
          fFlowQCCorPro[fCenBin][hr][2][charge]->Fill(FillPtBin,dQC4,WdQM4*fCenWeightEbE);
          dQ4f = kTRUE;
        }

        // product of correlations or covariances
        if(Q2f && dQ2f) fFlowQCCorCovPro[fCenBin][hr][0][charge]->Fill(FillPtBin,IQC2[hr]*dQC2,WQM2*WdQM2*fCenWeightEbE);
        if(Q4f && dQ2f) fFlowQCCorCovPro[fCenBin][hr][1][charge]->Fill(FillPtBin,IQC4[hr]*dQC2,WQM4*WdQM2*fCenWeightEbE);
        if(Q2f && dQ4f) fFlowQCCorCovPro[fCenBin][hr][2][charge]->Fill(FillPtBin,IQC2[hr]*dQC4,WQM2*WdQM4*fCenWeightEbE);
        if(dQ2f && dQ4f) fFlowQCCorCovPro[fCenBin][hr][3][charge]->Fill(FillPtBin,dQC2*dQC4,WdQM2*WdQM4*fCenWeightEbE);
        if(Q4f && dQ4f) fFlowQCCorCovPro[fCenBin][hr][4][charge]->Fill(FillPtBin,IQC4[hr]*dQC4,WQM4*WdQM4*fCenWeightEbE);

      } // end of for(Int_t pt=0; pt<fCRCnPtBin; pt++)

    } // end of for(Int_t hr=0; hr<fFlowNHarm; hr++)



    // pt diff get fFlowQCCorHist ==============================


    // Pt-DIFFERENTIAL
    for(Int_t hr=0; hr<fFlowNHarm; hr++) {
      for(Int_t j=0; j<fFlowQCNRef; j++) {
        for(Int_t pt=1;pt<=fNBins;pt++) { // fNBins was: fFlowQCRefCorPro[hr][j][charge]->GetNbinsX()
          Double_t stats[6]={0.};
          fFlowQCRefCorPro[hr][j][charge]->GetXaxis()->SetRange(pt,pt);
          fFlowQCRefCorPro[hr][j][charge]->GetStats(stats);
          Double_t SumWeig   = stats[0];
          Double_t SumWeigSq  = stats[1];
          Double_t SumTwo  = stats[4];
          Double_t SumTwoSq = stats[5];

          if(SumWeig>0.) {
            Double_t Corr = SumTwo/SumWeig;
            Double_t SqCorr = SumTwoSq/SumWeig;
            Double_t Weig = SumWeig;
            Double_t SqWeig = SumWeigSq;
            Double_t spread=0., termA=0., termB=0.;
            if(SqCorr-pow(Corr,2.)>=0.) { spread = pow(SqCorr-pow(Corr,2.),0.5); }
            if(TMath::Abs(Weig)>0.) { termA = (pow(SqWeig,0.5)/Weig); }
            if(1.-pow(termA,2.)>0.) { termB = 1./pow(1.-pow(termA,2.),0.5); }
            Double_t CorrErr = termA*spread*termB; // final error (unbiased estimator for standard deviation)
            if(CorrErr) {
              fFlowQCRefCorHist[hr][j][charge]->SetBinContent(pt,Corr);
              fFlowQCRefCorHist[hr][j][charge]->SetBinError(pt,CorrErr);
            }
          }
        } // end of for(Int_t pt=1;pt<=100;pt++)
        fFlowQCRefCorPro[hr][j][charge]->GetXaxis()->SetRange(1,fFlowQCRefCorPro[hr][j][charge]->GetNbinsX());
      } // end of for(Int_t j=0; j<5; j++)


      for (Int_t h=0; h<fCRCnCen; h++) {

        // STORE IN HISTOGRAMS

        for(Int_t j=0; j<fFlowQCNPro; j++) {
          for(Int_t pt=1;pt<=fNBins;pt++) { //was 13

            Double_t stats[6]={0.};
            fFlowQCCorPro[h][hr][j][charge]->GetXaxis()->SetRange(pt,pt);
            fFlowQCCorPro[h][hr][j][charge]->GetStats(stats);
            Double_t SumWeig   = stats[0]; //SUM(w)
            Double_t SumWeigSq  = stats[1]; //SUM(w^2)
            Double_t SumTwo  = stats[4]; //SUM(w*y)
            Double_t SumTwoSq = stats[5]; //SUM(w*y^2)

            if(SumWeig>0.) {
              Double_t Corr = SumTwo/SumWeig;
              Double_t SqCorr = SumTwoSq/SumWeig;
              Double_t Weig = SumWeig;
              Double_t SqWeig = SumWeigSq;
              Double_t spread=0., termA=0., termB=0.;
              if(SqCorr-pow(Corr,2.)>=0.) { spread = pow(SqCorr-pow(Corr,2.),0.5); } //sigma_y
              if(TMath::Abs(Weig)>0.) { termA = (pow(SqWeig,0.5)/Weig); }
              if(1.-pow(termA,2.)>0.) { termB = 1./pow(1.-pow(termA,2.),0.5); }
              Double_t CorrErr = termA*spread*termB; // final error (unbiased estimator for standard deviation)
              if(CorrErr) {
                fFlowQCCorHist[h][hr][j][charge]->SetBinContent(pt,Corr);
                //std::cout<<"for h="<<h<<" hr="<<hr<<" j="<<j<<" | Corr="<<Corr<< " for pt="<<pt<<"\n";
                fFlowQCCorHist[h][hr][j][charge]->SetBinError(pt,CorrErr);
              }
            }

          } // end of for(Int_t pt=1;pt<=fNBins;pt++)
          fFlowQCCorPro[h][hr][j][charge]->GetXaxis()->SetRange(1,416);
        }

        // reference flow
        // 2- and 4-particle cumulants
        Double_t QC2    = fFlowQCRefCorHist[hr][0][charge]->GetBinContent(h+1);
        Double_t QC2E = fFlowQCRefCorHist[hr][0][charge]->GetBinError(h+1);
        Double_t QC4    = fFlowQCRefCorHist[hr][1][charge]->GetBinContent(h+1);
        Double_t QC4E = fFlowQCRefCorHist[hr][1][charge]->GetBinError(h+1);
        Double_t Cn2 = QC2;
        Double_t Cn2E = QC2E;
        Double_t wCov24 = fFlowQCRefCorHist[hr][13][charge]->GetBinContent(h+1);
        Double_t Cn4 = QC4-2.*QC2*QC2;
        Double_t Cn4Esq = 16.*pow(QC2,2.)*pow(QC2E,2) + pow(QC4E,2.) - 8.*QC2*wCov24;

        Double_t Cos1 = fFlowQCRefCorHist[hr][3][charge]->GetBinContent(h+1); // <<cos(n*phi1)>>
        Double_t Sin1 = fFlowQCRefCorHist[hr][4][charge]->GetBinContent(h+1); // <<sin(n*phi1)>>
        Double_t Sin1P2 = fFlowQCRefCorHist[hr][5][charge]->GetBinContent(h+1);
        Double_t Cos1P2 = fFlowQCRefCorHist[hr][6][charge]->GetBinContent(h+1);
        Double_t Sin1M2M3 = fFlowQCRefCorHist[hr][7][charge]->GetBinContent(h+1);
        Double_t Cos1M2M3 = fFlowQCRefCorHist[hr][8][charge]->GetBinContent(h+1);
        // change vocabulary, to be changed
        Double_t cosP1nPhi = fFlowQCRefCorHist[hr][3][charge]->GetBinContent(h+1); // <<cos(n*phi1)>>
        Double_t sinP1nPhi = fFlowQCRefCorHist[hr][4][charge]->GetBinContent(h+1); // <<sin(n*phi1)>>
        Double_t sinP1nPhi1P1nPhi2 = fFlowQCRefCorHist[hr][5][charge]->GetBinContent(h+1); //sin(n*(phi1+phi2))
        Double_t cosP1nPhi1P1nPhi2 = fFlowQCRefCorHist[hr][6][charge]->GetBinContent(h+1);  //cos(n*(phi1+phi2))
        Double_t sinP1nPhi1M1nPhi2M1nPhi3 = fFlowQCRefCorHist[hr][7][charge]->GetBinContent(h+1);  //sin(n*(phi1-phi2-phi3))
        Double_t cosP1nPhi1M1nPhi2M1nPhi3 = fFlowQCRefCorHist[hr][8][charge]->GetBinContent(h+1); //cos(n*(phi1-phi2-phi3))

        fFlowQCRefCorFinal[hr][0][charge]->SetBinContent(h+1,Cn2);
        fFlowQCRefCorFinal[hr][0][charge]->SetBinError(h+1,Cn2E);

        if(Cn4Esq>0.) {
          Double_t Cn4E = pow(Cn4Esq,0.5);
          fFlowQCRefCorFinal[hr][3][charge]->SetBinContent(h+1,Cn4);
          fFlowQCRefCorFinal[hr][3][charge]->SetBinError(h+1,Cn4E);
          if(Cn4<0.) {
            Double_t Flow4 = pow(fabs(Cn4),0.25);
            Double_t Flow4E = fabs(Flow4/(4.*Cn4))*Cn4E;
            fFlowQCRefCorFinal[hr][1][charge]->SetBinContent(h+1,Flow4);
            fFlowQCRefCorFinal[hr][1][charge]->SetBinError(h+1,Flow4E);
          }
        }

        // pt-differential
        for(Int_t pt=1; pt<=fNBins; pt++) {
          Double_t qp2    = fFlowQCCorHist[h][hr][1][charge]->GetBinContent(pt);
          //std::cout<<"for h="<<h<<" hr="<<hr<<" j="<<1<<" | pt: "<<pt<<" qp2="<<qp2<<std::endl;
          Double_t qp2E = fFlowQCCorHist[h][hr][1][charge]->GetBinError(pt);
          Double_t qp4    = fFlowQCCorHist[h][hr][2][charge]->GetBinContent(pt);
          Double_t qp4E = fFlowQCCorHist[h][hr][2][charge]->GetBinError(pt);
          Double_t Dn2 = qp2;
          Double_t Dn2E = qp2E;
          Double_t Dn4 = qp4-2.*qp2*QC2;
          Double_t wCovTwoFourReduced = fFlowQCCorCovHist[h][hr][1][charge]->GetBinContent(pt);
          Double_t wCovTwoReducedFourReduced = fFlowQCCorCovHist[h][hr][4][charge]->GetBinContent(pt);
          Double_t Dn4Esq = 4.*pow(QC2,2.)*pow(qp2E,2) + 4.*pow(qp2,2.)*pow(QC2E,2) + pow(qp4E,2.) - 4.*qp2*wCovTwoFourReduced - 4.*QC2*wCovTwoReducedFourReduced;


          fFlowQCFinalPtDifHist[h][hr][5][charge]->SetBinContent(pt,Dn2);
          fFlowQCFinalPtDifHist[h][hr][5][charge]->SetBinError(pt,Dn2E);

          if(Cn2) {
            Double_t Flow2 = Dn2/sqrt(fabs(Cn2));
            //std::cout<<"for h="<<h<<" hr="<<hr<<" j="<<1<<" | pt: "<<pt<<" Flow2="<<Flow2<<std::endl;
            Double_t Flow2E = 0.;
            // change vocabulary, to be changed
            Double_t two = QC2; //std::cout<<two<<std::endl;
            Double_t twoError = QC2E; //std::cout<<twoError<<std::endl;
            Double_t twoReduced = qp2; //std::cout<<twoReduced<<std::endl;
            Double_t twoReducedError = qp2E; //std::cout<<twoReducedError<<std::endl;
            Double_t wCovTwoTwoReduced = fFlowQCCorCovHist[h][hr][0][charge]->GetBinContent(pt);
            Double_t v2PrimeErrorSquared = (1./4.)*pow(two,-3.)*(pow(twoReduced,2.)*pow(twoError,2.)
                                                                 + 4.*pow(two,2.)*pow(twoReducedError,2.)
                                                                 - 4.*two*twoReduced*wCovTwoTwoReduced);
            if(v2PrimeErrorSquared>0.){Flow2E = pow(v2PrimeErrorSquared,0.5);}

            
            if(Flow2E>0.) {
//              std::cout<<"1766: komen we hier -> Ptdiff flow 2? zo ja Flow2 = "<<Flow2<<std::endl;
              
              fFlowQCFinalPtDifHist[h][hr][0][charge]->SetBinContent(pt,Flow2);
              fFlowQCFinalPtDifHist[h][hr][0][charge]->SetBinError(pt,Flow2E);
            }
          }

          if(Dn4Esq>0.) {
            Double_t Dn4E = pow(Dn4Esq,0.5);
            fFlowQCFinalPtDifHist[h][hr][6][charge]->SetBinContent(pt,Dn4);
            fFlowQCFinalPtDifHist[h][hr][6][charge]->SetBinError(pt,Dn4E);
          }

          if(Cn4Esq>0.) {
            Double_t Flow4 = - Dn4/pow(fabs(Cn4),0.75);
            Double_t Flow4E = 0.;
            // change vocabulary, to be changed
            Double_t two = QC2;
            Double_t twoError = QC2E;
            Double_t twoReduced = qp2;
            Double_t twoReducedError = qp2E;
            Double_t four = QC4;
            Double_t fourError = QC4E;
            Double_t fourReduced = qp4;
            Double_t fourReducedError = qp4E;
            Double_t wCovTwoTwoReduced = fFlowQCCorCovHist[h][hr][0][charge]->GetBinContent(pt);
            Double_t wCovTwoFourReduced = fFlowQCCorCovHist[h][hr][1][charge]->GetBinContent(pt);
            Double_t wCovFourTwoReduced = fFlowQCCorCovHist[h][hr][2][charge]->GetBinContent(pt);
            Double_t wCovFourFourReduced = fFlowQCCorCovHist[h][hr][3][charge]->GetBinContent(pt);
            Double_t wCovTwoReducedFourReduced = fFlowQCCorCovHist[h][hr][4][charge]->GetBinContent(pt);
            Double_t wCovTwoFour = fFlowQCRefCorHist[hr][13][charge]->GetBinContent(h+1);

            Double_t v4PrimeErrorSquared = 0.;
            if(2.*pow(two,2.)-four>0.) {
              v4PrimeErrorSquared = pow(2.*pow(two,2.)-four,-7./2.)
              * (pow(2.*pow(two,2.)*twoReduced-3.*two*fourReduced+2.*four*twoReduced,2.)*pow(twoError,2.)
                 + (9./16.)*pow(2.*two*twoReduced-fourReduced,2.)*pow(fourError,2.)
                 + 4.*pow(two,2.)*pow(2.*pow(two,2.)-four,2.)*pow(twoReducedError,2.)
                 + pow(2.*pow(two,2.)-four,2.)*pow(fourReducedError,2.)
                 - (3./2.)*(2.*two*twoReduced-fourReduced)
                 * (2.*pow(two,2.)*twoReduced-3.*two*fourReduced+2.*four*twoReduced)*wCovTwoFour
                 - 4.*two*(2.*pow(two,2.)-four)
                 * (2.*pow(two,2.)*twoReduced-3.*two*fourReduced+2.*four*twoReduced)*wCovTwoTwoReduced
                 + 2.*(2.*pow(two,2.)-four)
                 * (2.*pow(two,2.)*twoReduced-3.*two*fourReduced+2.*four*twoReduced)*wCovTwoFourReduced
                 + 3.*two*(2.*pow(two,2.)-four)*(2.*two*twoReduced-fourReduced)*wCovFourTwoReduced
                 - (3./2.)*(2.*pow(two,2.)-four)*(2.*two*twoReduced-fourReduced)*wCovFourFourReduced
                 - 4.*two*pow(2.*pow(two,2.)-four,2.)*wCovTwoReducedFourReduced);
            }
            if(v4PrimeErrorSquared>0.){Flow4E = pow(v4PrimeErrorSquared,0.5);}

            if(Flow4E>0.) {
              fFlowQCFinalPtDifHist[h][hr][1][charge]->SetBinContent(pt,Flow4);
              fFlowQCFinalPtDifHist[h][hr][1][charge]->SetBinError(pt,Flow4E);
            }
          }
        } // end of for(Int_t pt=1; pt<=fNBins; pt++) {
      } // end of for (Int_t h=0; h<fCRCnCen; h++) {
    } // end of for(Int_t hr=0; hr<fFlowNHarm; hr++)
  } //end of for(Int_t charge=0;charge<fCharge;charge++)
}

//=======================================================================================================================

//=======================================================================================================================
void CalculateFlow::FinalizeFlowQC()
{
  std::cout << "Finalizing Flow QC"<<'\n'<<endl;
  std::cout << "Print bins where v_n{2} or v_n{4} cannot be determined"<<endl;
  cout << "*************************************" << "\n";

  for(Int_t charge=0; charge<fCharge; charge++){
    for(Int_t hr=0; hr<fFlowNHarm; hr++) {
      // Pt-INTEGRATED
      // STORE IN HISTOGRAMS
      // 2- and 4-particle cumulants
      for(Int_t j=0; j<fkFlowQCnIntCorPro; j++) {
        for(Int_t pt=1;pt<=fFlowQCIntCorPro[hr][j][charge]->GetNbinsX();pt++) { //pt is hier helamaal geen pt
          Double_t stats[6]={0.};
          fFlowQCIntCorPro[hr][j][charge]->GetXaxis()->SetRange(pt,pt);
          fFlowQCIntCorPro[hr][j][charge]->GetStats(stats);
          LongDouble_t SumWeig   = stats[0];
          LongDouble_t SumWeigSq  = stats[1];
          LongDouble_t SumTwo  = stats[4];
          LongDouble_t SumTwoSq = stats[5];

          if(SumWeig>0.) {
            LongDouble_t Corr = SumTwo/SumWeig;
            LongDouble_t SqCorr = SumTwoSq/SumWeig;
            LongDouble_t Weig = SumWeig;
            LongDouble_t SqWeig = SumWeigSq;
            LongDouble_t spread=0., termA=0., termB=0.;
            if(SqCorr-pow(Corr,2.)>=0.) { spread = pow(SqCorr-pow(Corr,2.),0.5); }
            if(TMath::Abs(Weig)>0.) { termA = (pow(SqWeig,0.5)/Weig); }
            if(1.-pow(termA,2.)>0.) { termB = 1./pow(1.-pow(termA,2.),0.5); }
            LongDouble_t CorrErr = termA*spread*termB; // final error (unbiased estimator for standard deviation)

            if(CorrErr) {
              fFlowQCIntCorHist[hr][j][charge]->SetBinContent(pt,Corr); //= <<2>> = cn{2} = vn{2}^2
              fFlowQCIntCorHist[hr][j][charge]->SetBinError(pt,CorrErr);

              if(Corr>0.){
                Double_t Flow2 = pow(fabs(Corr),0.5);
                Double_t Flow2E = fabs(Flow2/(4.*Corr))*CorrErr;

                if(j==0){
                  fFlowQCIntFlow2Hist[hr][0][charge]->SetBinContent(pt,Flow2); //neem sqrt(Corr) -> v_n{2}
                  fFlowQCIntFlow2Hist[hr][0][charge]->SetBinError(pt,Flow2E);}
              }

            }
          }
        } // end of for(Int_t pt=1;pt<=100;pt++)
        fFlowQCIntCorPro[hr][j][charge]->GetXaxis()->SetRange(1,fFlowQCIntCorPro[hr][j][charge]->GetNbinsX());
      } // end of for(Int_t j=0; j<5; j++)

    } // end of for(Int_t hr=0; hr<fFlowNHarm; hr++)
  }// end of for(Int_t charge=0; charge<fCharge; charge++)



  // FINALISE (calculate flow)
  for(Int_t charge=0; charge<fCharge; charge++){
    for(Int_t hr=0; hr<fFlowNHarm; hr++) {
      // calculate covariance
      for(Int_t pt=1; pt<=fFlowQCIntCorHist[hr][0][charge]->GetNbinsX(); pt++) {
        // average reduced correlations:
        Double_t two = fFlowQCIntCorHist[hr][0][charge]->GetBinContent(pt); // <<2>>
        Double_t four = fFlowQCIntCorHist[hr][1][charge]->GetBinContent(pt); // <<4>>
        // sum of weights for reduced correlation:
        Double_t sumOfWeightsForTwo = GetSumPro(fFlowQCIntCorPro[hr][0][charge],pt); // sum_{i=1}^{N} w_{<2>}
        Double_t sumOfWeightsForFour = GetSumPro(fFlowQCIntCorPro[hr][1][charge],pt); // sum_{i=1}^{N} w_{<4>}
        // product of weights for reduced correlation:
        Double_t productOfWeightsForTwoFour = GetSumPro(fFlowQCIntCorPro[hr][2][charge],pt); // sum_{i=1}^{N} w_{<2>}w_{<4>}
        // products for differential flow:
        Double_t twoFour = fFlowQCIntCorHist[hr][2][charge]->GetBinContent(pt); // <<2><4>>

        // <2>,<4>:
        Double_t term1 = productOfWeightsForTwoFour;
        Double_t term2 = sumOfWeightsForTwo;
        Double_t term3 = sumOfWeightsForFour;
        if(term2*term3>0.)
        {
          Double_t denominator = 1.-term1/(term2*term3);
          Double_t prefactor = term1/(term2*term3);
          if(TMath::Abs(denominator)>1.e-6)
          {
            Double_t covTwoFour = (twoFour-two*four)/denominator;
            Double_t wCovTwoFour = covTwoFour*prefactor;
            fFlowQCIntCorHist[hr][2][charge]->SetBinContent(pt,wCovTwoFour);
          }
        }
      } // end of for(Int_t pt=1;pt<=fNBins;pt++)

      // 2- and 4-particle cumulants
      for(Int_t pt=1; pt<=fFlowQCIntCorHist[hr][0][charge]->GetNbinsX(); pt++) {
        Double_t QC2    = fFlowQCIntCorHist[hr][0][charge]->GetBinContent(pt); //<<2>>
        Double_t QC2E   = fFlowQCIntCorHist[hr][0][charge]->GetBinError(pt);
        Double_t QC4    = fFlowQCIntCorHist[hr][1][charge]->GetBinContent(pt); //<<4>>
        Double_t QC4E   = fFlowQCIntCorHist[hr][1][charge]->GetBinError(pt);
        Double_t wCov24 = fFlowQCIntCorHist[hr][2][charge]->GetBinContent(pt); //<<2>><<4>>
        Double_t Cn2 = QC2;
        Double_t Cn2E = QC2E;
        Double_t Cn4 = QC4-2.*QC2*QC2;
        Double_t Cn4Esq = 16.*pow(QC2,2.)*pow(QC2E,2) + pow(QC4E,2.) - 8.*QC2*wCov24;




        if(Cn2<0.){
          // std::cout<<"QC2: "<<QC2<<" QC2E: "<<QC2E<<" QC4E: "<<QC4E<<" wCov42: "<<wCov24<<endl;
          // std::cout<<" Cn4Esq = 16.*pow(QC2,2.)*pow(QC2E,2) + pow(QC4E,2.) - 8.*QC2*wCov24 "<<endl;
          std::cout<<"v_"<<hr+1<<"{2}   "<<" bin is: "<<pt<< " charge is: " <<charge<<" Cn2: "<<Cn2<<" Cn2E: "<<QC2E<<"\n";}
        if(Cn4>0. || Cn4Esq<0){
          //  std::cout<<"QC2: "<<QC2<<" QC2E: "<<QC2E<<" QC4E: "<<QC4E<<" wCov42: "<<wCov24<<endl;
          //  std::cout<<" Cn4Esq = 16.*pow(QC2,2.)*pow(QC2E,2) + pow(QC4E,2.) - 8.*QC2*wCov24 "<<endl;
          std::cout<<"v_"<<hr+1<<"{4}   "<<" bin is: "<<pt<< " charge is: " <<charge<<" Cn4: "<<Cn4<<" Cn4Esq: "<<Cn4Esq<<"\n";}


        fFlowQCIntCumHist[hr][0][charge]->SetBinContent(pt,Cn2);
        fFlowQCIntCumHist[hr][0][charge]->SetBinError(pt,Cn2E);

        if(Cn4Esq>0.) {
          Double_t Cn4E = pow(Cn4Esq,0.5);
          fFlowQCIntCumHist[hr][1][charge]->SetBinContent(pt,Cn4);
          fFlowQCIntCumHist[hr][1][charge]->SetBinError(pt,Cn4E);

          if (Cn4<0.) {
            Double_t Flow4 = pow(fabs(Cn4),0.25);
            Double_t Flow4E = fabs(Flow4/(4.*Cn4))*Cn4E;

            fFlowQCIntCorHist[hr][2][charge]->SetBinContent(pt,Flow4);
            fFlowQCIntCorHist[hr][2][charge]->SetBinError(pt,Flow4E);

            fFlowQCIntFlow4Hist[hr][0][charge]->SetBinContent(pt,Flow4); //neem sqrt(Corr) -> v_n{2}
            fFlowQCIntFlow4Hist[hr][0][charge]->SetBinError(pt,Flow4E);
          }

          //                    } else {
          //                        fFlowQCIntCorHist[hr][2][charge]->SetBinContent(pt,0.);
          //                        fFlowQCIntCorHist[hr][2][charge]->SetBinError(pt,0.);
          //
          //                       // fFlowQCIntFlow4Hist[hr][0][charge]->SetBinContent(pt,0.); //neem sqrt(Corr) -> v_n{2}
          //                       // fFlowQCIntFlow4Hist[hr][0][charge]->SetBinError(pt,0.); //
          //
          //                    }
        }
      }
    } // end of for(Int_t hr=0; hr<fFlowNHarm; hr++)
  } // end of for(Int_t charge=0; charge<fCharge; charge++)
}

//=======================================================================================================================
void CalculateFlow::FinalizeFlowGF()
{
    cout << "*************************************" << "\n";
    cout << "\n";
    cout << "Finalizing Flow with Generic Framework";
    cout << "\n";
    cout << "\n";
    
    for (Int_t h=0; h<fkFlowGFNHarm; h++) {
        for(Int_t i=0; i<fkFlowGFNOrde; i++) {
            for(Int_t pt=1; pt<=fFlowGFIntCorPro[h][i]->GetNbinsX(); pt++) {
                Double_t stats[6]={0.};
                fFlowGFIntCorPro[h][i]->GetXaxis()->SetRange(pt,pt);
                fFlowGFIntCorPro[h][i]->GetStats(stats);
                Double_t SumWeig   = stats[0];
                Double_t SumWeigSq  = stats[1];
                Double_t SumTwo  = stats[4];
                Double_t SumTwoSq = stats[5];
                if(SumWeig>0.) {
                    Double_t Corr = SumTwo/SumWeig;
                    Double_t SqCorr = SumTwoSq/SumWeig;
                    Double_t Weig = SumWeig;
                    Double_t SqWeig = SumWeigSq;
                    Double_t spread=0., termA=0., termB=0.;
                    if(SqCorr-pow(Corr,2.)>=0.) { spread = pow(SqCorr-pow(Corr,2.),0.5); }
                    if(TMath::Abs(Weig)>0.) { termA = (pow(SqWeig,0.5)/Weig); }
                    if(1.-pow(termA,2.)>0.) { termB = 1./pow(1.-pow(termA,2.),0.5); }
                    Double_t CorrErr = termA*spread*termB; // final error (unbiased estimator for standard deviation)
                    if(CorrErr) {
                        fFlowGFIntCorHist[h][i]->SetBinContent(pt,Corr);
                        fFlowGFIntCorHist[h][i]->SetBinError(pt,CorrErr);
                    }
                }
            } // end of for(Int_t pt=1;pt<=fNBins;pt++)
            fFlowGFIntCorPro[h][i]->GetXaxis()->SetRange(1,fFlowGFIntCorPro[h][i]->GetNbinsX());
        } // end of for(Int_t i=0; i<fkFlowGFNOrde; i++)
    }
    
    for (Int_t h=0; h<fkFlowGFNHarm; h++) {
        for(Int_t i=0; i<fkFlowGFNOrde; i++) {
            for(Int_t k=0; k<fkFlowGFNOrde; k++) {
                for(Int_t pt=1; pt<=fFlowGFIntCovPro[h][i][k]->GetNbinsX(); pt++) {
                    // correlations:
                    Double_t A = fFlowGFIntCorHist[h][i]->GetBinContent(pt); // <<A>>
                    Double_t B = fFlowGFIntCorHist[h][k]->GetBinContent(pt); // <<B>>
                    // sum of weights for correlation:
                    Double_t sumOfWeightsForA = GetSumPro(fFlowGFIntCorPro[h][i],pt); // sum_{i=1}^{N} w_{<A>}
                    Double_t sumOfWeightsForB = GetSumPro(fFlowGFIntCorPro[h][k],pt); // sum_{i=1}^{N} w_{<B>}
                    // products for correlations:
                    Double_t AB = fFlowGFIntCovPro[h][i][k]->GetBinContent(pt); // <<A><B>>
                    // sum of weights for product of correlation:
                    Double_t productOfWeightsForAB = GetSumPro(fFlowGFIntCovPro[h][i][k],pt); // sum_{i=1}^{N} w_{<A>}w_{<B>}
                    // <A>,<B>:
                    Double_t term1 = productOfWeightsForAB;
                    Double_t term2 = sumOfWeightsForA;
                    Double_t term3 = sumOfWeightsForB;
                    if(term2*term3>0.)
                    {
                        Double_t denominator = 1.-term1/(term2*term3);
                        Double_t prefactor = term1/(term2*term3);
                        if(TMath::Abs(denominator)>1.e-6)
                        {
                            Double_t covAB = (AB-A*B)/denominator;
                            Double_t wCovAB = covAB*prefactor;
                            fFlowGFIntCovHist[h][i][k]->SetBinContent(pt,wCovAB);
                        }
                    }
                } // end of for(Int_t pt=1; pt<=fFlowGFIntCovPro[h][i][k]->GetNbinsX(); pt++)
                fFlowGFIntCovPro[h][i][k]->GetXaxis()->SetRange(1,fFlowGFIntCovPro[h][i][k]->GetNbinsX());
                fFlowGFIntCorPro[h][k]->GetXaxis()->SetRange(1,fFlowGFIntCorPro[h][k]->GetNbinsX());
            } // end of for(Int_t k=0; k<fkFlowGFNOrde; k++)
            fFlowGFIntCorPro[h][i]->GetXaxis()->SetRange(1,fFlowGFIntCorPro[h][i]->GetNbinsX());
        }
    }
    
    for (Int_t h=0; h<fkFlowGFNHarm; h++) {
        for(Int_t pt=1; pt<=fFlowGFIntCorHist[h][0]->GetNbinsX(); pt++) {
            
            // Correlations:
            Double_t two = fFlowGFIntCorHist[h][0]->GetBinContent(pt); // <<2>>
            Double_t four = fFlowGFIntCorHist[h][1]->GetBinContent(pt); // <<4>>
            Double_t six = fFlowGFIntCorHist[h][2]->GetBinContent(pt); // <<6>>
            Double_t eight = fFlowGFIntCorHist[h][3]->GetBinContent(pt); // <<8>>
            // Statistical errors of average 2-, 4-, 6- and 8-particle azimuthal correlations:
            Double_t twoError = fFlowGFIntCorHist[h][0]->GetBinError(pt); // statistical error of <2>
            Double_t fourError = fFlowGFIntCorHist[h][1]->GetBinError(pt); // statistical error of <4>
            Double_t sixError = fFlowGFIntCorHist[h][2]->GetBinError(pt); // statistical error of <6>
            Double_t eightError = fFlowGFIntCorHist[h][3]->GetBinError(pt); // statistical error of <8>
            
            // Q-cumulants:
            Double_t qc2 = 0.; // QC{2}
            Double_t qc4 = 0.; // QC{4}
            Double_t qc6 = 0.; // QC{6}
            Double_t qc8 = 0.; // QC{8}
            if(TMath::Abs(two) > 0.){qc2 = two;}
            if(TMath::Abs(four) > 0.){qc4 = four-2.*pow(two,2.);}
            if(TMath::Abs(six) > 0.){qc6 = six-9.*two*four+12.*pow(two,3.);}
            if(TMath::Abs(eight) > 0.){qc8 = eight-16.*two*six-18.*pow(four,2.)+144.*pow(two,2.)*four-144.*pow(two,4.);}
            // Statistical errors of Q-cumulants:
            Double_t qc2Error = 0.;
            Double_t qc4Error = 0.;
            Double_t qc6Error = 0.;
            Double_t qc8Error = 0.;
            // Squared statistical errors of Q-cumulants:
            //Double_t qc2ErrorSquared = 0.;
            Double_t qc4ErrorSquared = 0.;
            Double_t qc6ErrorSquared = 0.;
            Double_t qc8ErrorSquared = 0.;
            // covariances:
            Double_t wCov24 = fFlowGFIntCovHist[h][0][1]->GetBinContent(pt);
            Double_t wCov26 = fFlowGFIntCovHist[h][0][2]->GetBinContent(pt);
            Double_t wCov28 = fFlowGFIntCovHist[h][0][3]->GetBinContent(pt);
            Double_t wCov46 = fFlowGFIntCovHist[h][1][2]->GetBinContent(pt);
            Double_t wCov48 = fFlowGFIntCovHist[h][1][3]->GetBinContent(pt);
            Double_t wCov68 = fFlowGFIntCovHist[h][2][3]->GetBinContent(pt);
            
            // Statistical error of QC{2}:
            qc2Error = twoError;
            // Statistical error of QC{4}:
            qc4ErrorSquared = 16.*pow(two,2.)*pow(twoError,2.)+pow(fourError,2.)
            - 8.*two*wCov24;
            if(qc4ErrorSquared>0.) {
                qc4Error = pow(qc4ErrorSquared,0.5);
            }
            // Statistical error of QC{6}:
            qc6ErrorSquared = 81.*pow(4.*pow(two,2.)-four,2.)*pow(twoError,2.)
            + 81.*pow(two,2.)*pow(fourError,2.)
            + pow(sixError,2.)
            - 162.*two*(4.*pow(two,2.)-four)*wCov24
            + 18.*(4.*pow(two,2.)-four)*wCov26
            - 18.*two*wCov46;
            if(qc6ErrorSquared>0.) {
                qc6Error = pow(qc6ErrorSquared,0.5);
            }
            // Statistical error of QC{8}:
            qc8ErrorSquared = 256.*pow(36.*pow(two,3.)-18.*four*two+six,2.)*pow(twoError,2.)
            + 1296.*pow(4.*pow(two,2.)-four,2.)*pow(fourError,2.)
            + 256.*pow(two,2.)*pow(sixError,2.)
            + pow(eightError,2.)
            - 1152.*(36.*pow(two,3.)-18.*four*two+six)*(4.*pow(two,2.)-four)*wCov24
            + 512.*two*(36.*pow(two,3.)-18.*four*two+six)*wCov26
            - 32.*(36.*pow(two,3.)-18.*four*two+six)*wCov28
            - 1152.*two*(4.*pow(two,2.)-four)*wCov46
            + 72.*(4.*pow(two,2.)-four)*wCov48
            - 32.*two*wCov68;
            if(qc8ErrorSquared>0.) {
                qc8Error = pow(qc8ErrorSquared,0.5);
            }
            // Store the cumulants:
            fFlowGFIntCumHist[h][0]->SetBinContent(pt,qc2);
            fFlowGFIntCumHist[h][0]->SetBinError(pt,qc2Error);
            fFlowGFIntCumHist[h][1]->SetBinContent(pt,qc4);
            fFlowGFIntCumHist[h][1]->SetBinError(pt,qc4Error);
            fFlowGFIntCumHist[h][2]->SetBinContent(pt,qc6);
            fFlowGFIntCumHist[h][2]->SetBinError(pt,qc6Error);
            fFlowGFIntCumHist[h][3]->SetBinContent(pt,qc8);
            fFlowGFIntCumHist[h][3]->SetBinError(pt,qc8Error);
            
            // Reference flow estimates:
            Double_t v2 = 0.; // v{2,QC}
            Double_t v4 = 0.; // v{4,QC}
            Double_t v6 = 0.; // v{6,QC}
            Double_t v8 = 0.; // v{8,QC}
            // Reference flow statistical errors:
            Double_t v2Error = 0.; // v{2,QC} stat. error
            Double_t v4Error = 0.; // v{4,QC} stat. error
            Double_t v6Error = 0.; // v{6,QC} stat. error
            Double_t v8Error = 0.; // v{8,QC} stat. error
            // calculate flow
            if(qc2>=0.){v2 = pow(qc2,0.5);}
            if(qc4<=0.){v4 = pow(-1.*qc4,1./4.);}
            if(qc6>=0.){v6 = pow((1./4.)*qc6,1./6.);}
            if(qc8<=0.){v8 = pow((-1./33.)*qc8,1./8.);}
            // Calculate stat. error for reference flow estimates from stat. error of Q-cumulants:
            if(qc2>0.){v2Error = (1./2.)*pow(qc2,-0.5)*qc2Error;}
            if(qc4<0.){v4Error = (1./4.)*pow(-qc4,-3./4.)*qc4Error;}
            if(qc6>0.){v6Error = (1./6.)*pow(2.,-1./3.)*pow(qc6,-5./6.)*qc6Error;}
            if(qc8<0.){v8Error = (1./8.)*pow(33.,-1./8.)*pow(-qc8,-7./8.)*qc8Error;}
            // Store the results:
            if(h==0){
                std::cout<<"qc2: "<<qc2<<" qc4: "<<qc4<<" qc6: "<<qc6<<" qc8: "<<qc8<<"\n";
            }
            if(qc2>0.) {
              
                fFlowGFIntFinalHist[h][0]->SetBinContent(pt,v2);        //v_n{2}
                fFlowGFIntFinalHist[h][0]->SetBinError(pt,v2Error);
            }
            if(qc4<0.) {
                fFlowGFIntFinalHist[h][1]->SetBinContent(pt,v4);
                fFlowGFIntFinalHist[h][1]->SetBinError(pt,v4Error);
            }
            if(qc6>0.) {
                fFlowGFIntFinalHist[h][2]->SetBinContent(pt,v6);
                fFlowGFIntFinalHist[h][2]->SetBinError(pt,v6Error);
            }
            if(qc8<0.) {
                fFlowGFIntFinalHist[h][3]->SetBinContent(pt,v8);
                fFlowGFIntFinalHist[h][3]->SetBinError(pt,v8Error);
            }
            
        }
    }
    
    // MIXED HARMONICS ***********************************************************
    
    for (Int_t h=0; h<fkFlowGFNHarm; h++) {
        for(Int_t i=0; i<fkFlowGFNHarm; i++) {
            
            for(Int_t pt=1; pt<=fFlowGFMixedCorPro[h][i]->GetNbinsX(); pt++) {
                
                Double_t stats[6]={0.};
                fFlowGFMixedCorPro[h][i]->GetXaxis()->SetRange(pt,pt);
                fFlowGFMixedCorPro[h][i]->GetStats(stats);
                Double_t SumWeig   = stats[0];
                Double_t SumWeigSq  = stats[1];
                Double_t SumTwo  = stats[4];
                Double_t SumTwoSq = stats[5];
                
                if(SumWeig>0.) {
                    Double_t Corr = SumTwo/SumWeig;
                    Double_t SqCorr = SumTwoSq/SumWeig;
                    Double_t Weig = SumWeig;
                    Double_t SqWeig = SumWeigSq;
                    Double_t spread=0., termA=0., termB=0.;
                    if(SqCorr-pow(Corr,2.)>=0.) { spread = pow(SqCorr-pow(Corr,2.),0.5); }
                    if(TMath::Abs(Weig)>0.) { termA = (pow(SqWeig,0.5)/Weig); }
                    if(1.-pow(termA,2.)>0.) { termB = 1./pow(1.-pow(termA,2.),0.5); }
                    Double_t CorrErr = termA*spread*termB; // final error (unbiased estimator for standard deviation)
                    if(CorrErr) {
                        fFlowGFMixedCorHist[h][i]->SetBinContent(pt,Corr);
                        fFlowGFMixedCorHist[h][i]->SetBinError(pt,CorrErr);
                    }
                }
                
            } // end of for(Int_t pt=1;pt<=fNBins;pt++)
            fFlowGFMixedCorPro[h][i]->GetXaxis()->SetRange(1,fFlowGFMixedCorPro[h][i]->GetNbinsX());
        }
    }
    
    for (Int_t h=0; h<fkFlowGFNHarm; h++) {
        for(Int_t i=0; i<fkFlowGFNHarm; i++) {
            
            for(Int_t pt=1; pt<=fFlowGFIntCorHist[h][0]->GetNbinsX(); pt++) {
                // Correlations:
                Double_t twoA = fFlowGFIntCumHist[h][0]->GetBinContent(pt); // <<2A>>
                Double_t twoB = fFlowGFIntCumHist[i][0]->GetBinContent(pt); // <<2B>>
                Double_t four = fFlowGFMixedCorHist[h][i]->GetBinContent(pt); // <<4>>
                // Statistical errors:
                Double_t twoAError = fFlowGFIntCumHist[h][0]->GetBinError(pt); // statistical error of <2A>
                Double_t twoBError = fFlowGFIntCumHist[i][0]->GetBinError(pt); // statistical error of <2B>
                Double_t fourError = fFlowGFMixedCorHist[h][i]->GetBinError(pt); // statistical error of <4>
                // Symmetric Cumulants:
                Double_t SC = four - twoA*twoB;
                Double_t SCErrorSquared = pow(fourError,2.) + pow(twoA*twoBError,2.) + pow(twoB*twoAError,2.); // TBI
                // Store the results:
                if(SCErrorSquared>0.) {
                    fFlowGFMixedFinalHist[h][i]->SetBinContent(pt,SC);
                    fFlowGFMixedFinalHist[h][i]->SetBinError(pt,pow(SCErrorSquared,0.5));
                }
            }  // end of for(Int_t pt=1;pt<=fNBins;pt++)
            
        }
    }
    
    // in wide pt bins
    for(Int_t s=0; s<fkGFPtB; s++) {
        
        if(!fFlowGFIntCorProPtB[0][0][0]) continue;
        
        for (Int_t h=0; h<fkFlowGFNHarm; h++) {
            for(Int_t i=0; i<fkFlowGFNOrde; i++) {
                for(Int_t pt=1; pt<=fFlowGFIntCorProPtB[s][h][i]->GetNbinsX(); pt++) {
                    Double_t stats[6]={0.};
                    fFlowGFIntCorProPtB[s][h][i]->GetXaxis()->SetRange(pt,pt);
                    fFlowGFIntCorProPtB[s][h][i]->GetStats(stats);
                    Double_t SumWeig   = stats[0];
                    Double_t SumWeigSq  = stats[1];
                    Double_t SumTwo  = stats[4];
                    Double_t SumTwoSq = stats[5];
                    if(SumWeig>0.) {
                        Double_t Corr = SumTwo/SumWeig;
                        Double_t SqCorr = SumTwoSq/SumWeig;
                        Double_t Weig = SumWeig;
                        Double_t SqWeig = SumWeigSq;
                        Double_t spread=0., termA=0., termB=0.;
                        if(SqCorr-pow(Corr,2.)>=0.) { spread = pow(SqCorr-pow(Corr,2.),0.5); }
                        if(TMath::Abs(Weig)>0.) { termA = (pow(SqWeig,0.5)/Weig); }
                        if(1.-pow(termA,2.)>0.) { termB = 1./pow(1.-pow(termA,2.),0.5); }
                        Double_t CorrErr = termA*spread*termB; // final error (unbiased estimator for standard deviation)
                        if(CorrErr) {
                            fFlowGFIntCorHistPtB[s][h][i]->SetBinContent(pt,Corr);
                            fFlowGFIntCorHistPtB[s][h][i]->SetBinError(pt,CorrErr);
                        }
                    }
                } // end of for(Int_t pt=1;pt<=fNBins;pt++)
                fFlowGFIntCorProPtB[s][h][i]->GetXaxis()->SetRange(1,fFlowGFIntCorProPtB[s][h][i]->GetNbinsX());
            } // end of for(Int_t i=0; i<fkFlowGFNOrde; i++)
        }
        
        for (Int_t h=0; h<fkFlowGFNHarm; h++) {
            for(Int_t i=0; i<fkFlowGFNOrde; i++) {
                for(Int_t k=0; k<fkFlowGFNOrde; k++) {
                    for(Int_t pt=1; pt<=fFlowGFIntCovProPtB[s][h][i][k]->GetNbinsX(); pt++) {
                        // correlations:
                        Double_t A = fFlowGFIntCorHistPtB[s][h][i]->GetBinContent(pt); // <<A>>
                        Double_t B = fFlowGFIntCorHistPtB[s][h][k]->GetBinContent(pt); // <<B>>
                        // sum of weights for correlation:
                        Double_t sumOfWeightsForA = GetSumPro(fFlowGFIntCorProPtB[s][h][i],pt); // sum_{i=1}^{N} w_{<A>}
                        Double_t sumOfWeightsForB = GetSumPro(fFlowGFIntCorProPtB[s][h][k],pt); // sum_{i=1}^{N} w_{<B>}
                        // products for correlations:
                        Double_t AB = fFlowGFIntCovProPtB[s][h][i][k]->GetBinContent(pt); // <<A><B>>
                        // sum of weights for product of correlation:
                        Double_t productOfWeightsForAB = GetSumPro(fFlowGFIntCovProPtB[s][h][i][k],pt); // sum_{i=1}^{N} w_{<A>}w_{<B>}
                        // <A>,<B>:
                        Double_t term1 = productOfWeightsForAB;
                        Double_t term2 = sumOfWeightsForA;
                        Double_t term3 = sumOfWeightsForB;
                        if(term2*term3>0.)
                        {
                            Double_t denominator = 1.-term1/(term2*term3);
                            Double_t prefactor = term1/(term2*term3);
                            if(TMath::Abs(denominator)>1.e-6)
                            {
                                Double_t covAB = (AB-A*B)/denominator;
                                Double_t wCovAB = covAB*prefactor;
                                fFlowGFIntCovHistPtB[s][h][i][k]->SetBinContent(pt,wCovAB);
                            }
                        }
                    } // end of for(Int_t pt=1; pt<=fFlowGFIntCovProPtB[s][h][i][k]->GetNbinsX(); pt++)
                    fFlowGFIntCovProPtB[s][h][i][k]->GetXaxis()->SetRange(1,fFlowGFIntCovProPtB[s][h][i][k]->GetNbinsX());
                    fFlowGFIntCorProPtB[s][h][k]->GetXaxis()->SetRange(1,fFlowGFIntCorProPtB[s][h][k]->GetNbinsX());
                } // end of for(Int_t k=0; k<fkFlowGFNOrde; k++)
                fFlowGFIntCorProPtB[s][h][i]->GetXaxis()->SetRange(1,fFlowGFIntCorProPtB[s][h][i]->GetNbinsX());
            }
        }
    }
    
    cout << "*************************************" << "\n";
    cout << "\n";
    
} // end of void CalculateFlow::FinalizeFlowGF()

//=====================================================================================================
//
//void CalculateFlow::FinalizeFlowQC()
//{
//    std::cout << "Finalizing Flow QC"<<'\n';
//    for(Int_t charge=0; charge<fCharge; charge++){
//        for(Int_t hr=0; hr<fFlowNHarm; hr++) {
//            // Pt-INTEGRATED
//            // STORE IN HISTOGRAMS
//            // 2- and 4-particle cumulants
//            for(Int_t j=0; j<fkFlowQCnIntCorPro; j++) {
//                for(Int_t pt=1;pt<=fFlowQCIntCorPro[hr][j][charge]->GetNbinsX();pt++) { //pt is hier helamaal geen pt
//                    Double_t stats[6]={0.};
//                    fFlowQCIntCorPro[hr][j][charge]->GetXaxis()->SetRange(pt,pt);
//                    fFlowQCIntCorPro[hr][j][charge]->GetStats(stats);
//                    LongDouble_t SumWeig   = stats[0];
//                    LongDouble_t SumWeigSq  = stats[1];
//                    LongDouble_t SumTwo  = stats[4];
//                    LongDouble_t SumTwoSq = stats[5];
//
//                    if(SumWeig>0.) {
//                        LongDouble_t Corr = SumTwo/SumWeig;
//                        LongDouble_t SqCorr = SumTwoSq/SumWeig;
//                        LongDouble_t Weig = SumWeig;
//                        LongDouble_t SqWeig = SumWeigSq;
//                        LongDouble_t spread=0., termA=0., termB=0.;
//                        if(SqCorr-pow(Corr,2.)>=0.) { spread = pow(SqCorr-pow(Corr,2.),0.5); }
//                        if(TMath::Abs(Weig)>0.) { termA = (pow(SqWeig,0.5)/Weig); }
//                        if(1.-pow(termA,2.)>0.) { termB = 1./pow(1.-pow(termA,2.),0.5); }
//                        LongDouble_t CorrErr = termA*spread*termB; // final error (unbiased estimator for standard deviation)
//
//                        if(CorrErr) {
//                            fFlowQCIntCorHist[hr][j][charge]->SetBinContent(pt,Corr); //neem sqrt(Corr) -> v_n{2}
//                            fFlowQCIntCorHist[hr][j][charge]->SetBinError(pt,CorrErr);
//
//                            if(Corr>0.){
//                                Double_t Flow2 = pow(fabs(Corr),0.5);
//                                Double_t Flow2E = fabs(Flow2/(4.*Corr))*CorrErr;
//
//                                if(j==0){
//                                    fFlowQCIntFlow2Hist[hr][0][charge]->SetBinContent(pt,Flow2); //neem sqrt(Corr) -> v_n{2}
//                                    fFlowQCIntFlow2Hist[hr][0][charge]->SetBinError(pt,Flow2E);}
//                            }
//
//                        }
//                    }
//                } // end of for(Int_t pt=1;pt<=100;pt++)
//                fFlowQCIntCorPro[hr][j][charge]->GetXaxis()->SetRange(1,fFlowQCIntCorPro[hr][j][charge]->GetNbinsX());
//            } // end of for(Int_t j=0; j<5; j++)
//
//        } // end of for(Int_t hr=0; hr<fFlowNHarm; hr++)
//    }// end of for(Int_t charge=0; charge<fCharge; charge++)
//
//
//
//    // FINALISE (calculate flow)
//    for(Int_t charge=0; charge<fCharge; charge++){
//        for(Int_t hr=0; hr<fFlowNHarm; hr++) {
//            // calculate covariance
//            for(Int_t pt=1; pt<=fFlowQCIntCorHist[hr][0][charge]->GetNbinsX(); pt++) {
//                // average reduced correlations:
//                Double_t two = fFlowQCIntCorHist[hr][0][charge]->GetBinContent(pt); // <<2>>
//                Double_t four = fFlowQCIntCorHist[hr][1][charge]->GetBinContent(pt); // <<4>>
//                // sum of weights for reduced correlation:
//                Double_t sumOfWeightsForTwo = GetSumPro(fFlowQCIntCorPro[hr][0][charge],pt); // sum_{i=1}^{N} w_{<2>}
//                Double_t sumOfWeightsForFour = GetSumPro(fFlowQCIntCorPro[hr][1][charge],pt); // sum_{i=1}^{N} w_{<4>}
//                // product of weights for reduced correlation:
//                Double_t productOfWeightsForTwoFour = GetSumPro(fFlowQCIntCorPro[hr][2][charge],pt); // sum_{i=1}^{N} w_{<2>}w_{<4>}
//                // products for differential flow:
//                Double_t twoFour = fFlowQCIntCorHist[hr][2][charge]->GetBinContent(pt); // <<2><4>>
//
//                // <2>,<4>:
//                Double_t term1 = productOfWeightsForTwoFour;
//                Double_t term2 = sumOfWeightsForTwo;
//                Double_t term3 = sumOfWeightsForFour;
//                if(term2*term3>0.)
//                {
//                    Double_t denominator = 1.-term1/(term2*term3);
//                    Double_t prefactor = term1/(term2*term3);
//                    if(TMath::Abs(denominator)>1.e-6)
//                    {
//                        Double_t covTwoFour = (twoFour-two*four)/denominator;
//                        Double_t wCovTwoFour = covTwoFour*prefactor;
//                        fFlowQCIntCorHist[hr][2][charge]->SetBinContent(pt,wCovTwoFour);
//                    }
//                }
//            } // end of for(Int_t pt=1;pt<=fNBins;pt++)
//
//            // 2- and 4-particle cumulants
//            for(Int_t pt=1; pt<=fFlowQCIntCorHist[hr][0][charge]->GetNbinsX(); pt++) {
//                Double_t QC2    = fFlowQCIntCorHist[hr][0][charge]->GetBinContent(pt);
//                Double_t QC2E   = fFlowQCIntCorHist[hr][0][charge]->GetBinError(pt);
//                Double_t QC4    = fFlowQCIntCorHist[hr][1][charge]->GetBinContent(pt);
//                Double_t QC4E   = fFlowQCIntCorHist[hr][1][charge]->GetBinError(pt);
//                Double_t wCov24 = fFlowQCIntCorHist[hr][2][charge]->GetBinContent(pt);
//                Double_t Cn2 = QC2;
//                Double_t Cn2E = QC2E;
//                Double_t Cn4 = QC4-2.*QC2*QC2;
//                Double_t Cn4Esq = 16.*pow(QC2,2.)*pow(QC2E,2) + pow(QC4E,2.) - 8.*QC2*wCov24;
//
//                std::cout<<" harmonic is: "<< hr+1 << " bin is: "<<pt<< "Cn4: "<<Cn4<<" Cn4Esq: "<<Cn4Esq<<"\n";
//                std::cout<<"-2QC2^2 is: "<<-2*QC2*QC2<<" QC4: "<<QC4<<endl;
//
//                fFlowQCIntCumHist[hr][0][charge]->SetBinContent(pt,Cn2);
//                fFlowQCIntCumHist[hr][0][charge]->SetBinError(pt,Cn2E);
//
//                if(Cn4Esq>0.) {
//                    Double_t Cn4E = pow(Cn4Esq,0.5);
//                    fFlowQCIntCumHist[hr][1][charge]->SetBinContent(pt,Cn4);
//                    fFlowQCIntCumHist[hr][1][charge]->SetBinError(pt,Cn4E);
//                    if (Cn4<0.) {
//                        Double_t Flow4 = pow(fabs(Cn4),0.25);
//                        Double_t Flow4E = fabs(Flow4/(4.*Cn4))*Cn4E;
//
//                        fFlowQCIntCorHist[hr][2][charge]->SetBinContent(pt,Flow4);
//                        fFlowQCIntCorHist[hr][2][charge]->SetBinError(pt,Flow4E);
//
//                        fFlowQCIntFlow4Hist[hr][0][charge]->SetBinContent(pt,Flow4); //neem sqrt(Corr) -> v_n{2}
//                        fFlowQCIntFlow4Hist[hr][0][charge]->SetBinError(pt,Flow4E);
//
//                    } else {
//                        fFlowQCIntCorHist[hr][2][charge]->SetBinContent(pt,0.);
//                        fFlowQCIntCorHist[hr][2][charge]->SetBinError(pt,0.);
//
//                        fFlowQCIntFlow4Hist[hr][0][charge]->SetBinContent(pt,0.); //neem sqrt(Corr) -> v_n{2}
//                        fFlowQCIntFlow4Hist[hr][0][charge]->SetBinError(pt,0.);
//
//                    }
//                }
//            }
//        } // end of for(Int_t hr=0; hr<fFlowNHarm; hr++)
//    } // end of for(Int_t charge=0; charge<fCharge; charge++)
//}

//=====================================================================================================

void CalculateFlow::FinalizeFlowSPM()
{
  std::cout << "Finalizing Flow SPM"<<'\n'<<endl;
  std::cout << "Average correlations over all events"<<endl;
  cout << "*************************************" << "\n";
  
  Float_t v_p =0.;
  Float_t v_t=0.;
  
  Float_t EP_res;
  
  //    for(Int_t charge=0; charge<fCharge; charge++){
  //        for(Int_t pt=1;pt<=fSPMEPresolutionPro[0][charge]->GetNbinsX();pt++) {
  //
  //        }
  //    }
  
  
  for(Int_t charge=0; charge<fCharge; charge++){
    for (Int_t h=0;h<fFlowNHarm;h++) {
      for(Int_t pt=1;pt<=fFlowSPMIntPro[h][charge]->GetNbinsX();pt++) {
        
        
        Float_t Corr_QQ_y = 0; Double_t CorrErr_QQ_y = 0;
        
        
        Corr_QQ_y = GetWeightedCorrelations(fFlowSPMIntPro[h][charge], pt);
        CorrErr_QQ_y = GetWeightedCorrelationsError(fFlowSPMIntPro[h][charge], pt);
        
      //  EP_res = GetWeightedCorrelations(fSPMEPresolutionPro[h][charge],pt);
        
        fFlowSPMIntFlow2Hist[h][charge]->SetBinContent(pt, Corr_QQ_y);
        fFlowSPMIntFlow2Hist[h][charge]->SetBinError(pt, CorrErr_QQ_y);
      }
    }
  }// end of for(Int_t charge=0; charge<fCharge; charge++)
//
//  for(Int_t charge=0; charge<fCharge; charge++){
//    for (Int_t h=0;h<fFlowNHarm;h++) {
//      for(Int_t pt=1;pt<=fFlowSPMIntPro[h][charge]->GetNbinsX();pt++) {
//
//
//        Float_t Corr_QQ_y = 0; Double_t CorrErr_QQ_y = 0;
//
//        Corr_QQ_y = GetWeightedCorrelations(fFlowSPM1IntPro[h][charge], pt);
//        CorrErr_QQ_y = GetWeightedCorrelationsError(fFlowSPM1IntPro[h][charge], pt);
//
//        //EP_res = GetWeightedCorrelations(fSPMEPresolutionPro[h][charge],pt);
//
//        fFlowSPM1IntFlow2Hist[h][charge]->SetBinContent(pt, Corr_QQ_y);
//        fFlowSPM1IntFlow2Hist[h][charge]->SetBinError(pt, CorrErr_QQ_y);
//      }
//    }
//  }//
}

//=====================================================================================================



void CalculateFlow::FinalizeFlowEPM()
{
  std::cout << "Finalizing Flow EPM"<<'\n'<<endl;
  std::cout << "Average correlations over all events"<<endl;
  cout << "*************************************" << "\n";
  
  Float_t v_p =0.;
  Float_t v_t=0.;
  
  for (Int_t h=0;h<fFlowNHarm;h++) {
    for(Int_t charge=0; charge<fCharge; charge++){
      
      
      for(Int_t pt=1;pt<=fFlowEPMIntPro_pos[h][charge]->GetNbinsX();pt++) {
        Float_t Corr = 0; Double_t CorrErr = 0;
        Corr = GetWeightedCorrelations(fFlowEPMIntPro_pos[h][charge], pt);
        CorrErr = GetWeightedCorrelationsError(fFlowEPMIntPro_pos[h][charge], pt);
        fFlowEPMIntFlow2Hist_pos[h][charge]->SetBinContent(pt, Corr); //Let op!! Hier moet min als eta<0 in make function
        fFlowEPMIntFlow2Hist_pos[h][charge]->SetBinError(pt, CorrErr);
      }
      
      for(Int_t pt=1;pt<=fFlowEPMIntPro_neg[h][charge]->GetNbinsX();pt++) {
        Float_t Corr = 0; Double_t CorrErr = 0;
        Corr = GetWeightedCorrelations(fFlowEPMIntPro_neg[h][charge], pt);
        CorrErr = GetWeightedCorrelationsError(fFlowEPMIntPro_neg[h][charge], pt);
        fFlowEPMIntFlow2Hist_neg[h][charge]->SetBinContent(pt, Corr); //Let op!! Hier moet min als eta<0 in make function
        fFlowEPMIntFlow2Hist_neg[h][charge]->SetBinError(pt, CorrErr);
      }
      
      for(Int_t pt=1;pt<=fFlowEPMIntPro[h][charge]->GetNbinsX();pt++) {
        Float_t Corr = 0; Double_t CorrErr = 0;
        Corr = GetWeightedCorrelations(fFlowEPMIntPro[h][charge], pt);
        CorrErr = GetWeightedCorrelationsError(fFlowEPMIntPro[h][charge], pt);
        fFlowEPMIntFlow2Hist[h][charge]->SetBinContent(pt, Corr); //Let op!! Hier moet min als eta<0 in make function
        fFlowEPMIntFlow2Hist[h][charge]->SetBinError(pt, CorrErr);
      }
      
      
      for(Int_t pt=1;pt<=fFlowEPMCorPro[h][charge]->GetNbinsX();pt++) {
        Float_t Corr = 0; Double_t CorrErr = 0;
        Corr= GetWeightedCorrelations(fFlowEPMCorPro[h][charge], pt);
        CorrErr = GetWeightedCorrelationsError(fFlowEPMCorPro[h][charge], pt);
        if(Corr && CorrErr){
          fFlowEPMPtFlow2Hist[h][charge]->SetBinContent(pt, Corr); //Let op!! Hier moet min als eta<0 in make function
          fFlowEPMPtFlow2Hist[h][charge]->SetBinError(pt, CorrErr);}
        else {
          fFlowEPMPtFlow2Hist[h][charge]->SetBinContent(pt, 0); //Let op!! Hier moet min als eta<0 in make function
        }
      }
      
    }  // end of for(Int_t charge=0; charge<fCharge; charge++)
  }// end of pt
}


void CalculateFlow::FinalizeSpectra(Int_t Nevents) {
  for(Int_t pt=1; pt<=fPosPionsSpectra->GetNbinsX(); pt++) {
    Double_t Normalization = 1/(2*Nevents*maxEtaCut*(fPosPionsSpectra->GetBinWidth(pt)));
    
    Double_t PosPionContent = fPosPionsSpectra->GetBinContent(pt);
    Double_t PosPionError = fPosPionsSpectra->GetBinError(pt);
    Double_t NegPionContent = fNegPionsSpectra->GetBinContent(pt);
    Double_t NegPionError = fNegPionsSpectra->GetBinError(pt);
    
    Double_t PosKaonContent = fPosKaonsSpectra->GetBinContent(pt);
    Double_t PosKaonError = fPosKaonsSpectra->GetBinError(pt);
    Double_t NegKaonContent = fNegKaonsSpectra->GetBinContent(pt);
    Double_t NegKaonError = fNegKaonsSpectra->GetBinError(pt);
    
    Double_t PosProtonContent = fPosProtonsSpectra->GetBinContent(pt);
    Double_t PosProtonError = fPosProtonsSpectra->GetBinError(pt);
    Double_t NegProtonContent = fNegProtonsSpectra->GetBinContent(pt);
    Double_t NegProtonError = fNegProtonsSpectra->GetBinError(pt);
    
    fPosPionsSpectra->SetBinContent(pt, PosPionContent*Normalization);
    fPosPionsSpectra->SetBinError(pt, PosPionError*Normalization);
    fNegPionsSpectra->SetBinContent(pt, NegPionContent*Normalization);
    fNegPionsSpectra->SetBinError(pt, NegPionError*Normalization);
    
    fPosKaonsSpectra->SetBinContent(pt, PosKaonContent*Normalization);
    fPosKaonsSpectra->SetBinError(pt, PosKaonError*Normalization);
    fNegKaonsSpectra->SetBinContent(pt, NegKaonContent*Normalization);
    fNegKaonsSpectra->SetBinError(pt, NegKaonError*Normalization);
    
    fPosProtonsSpectra->SetBinContent(pt, PosProtonContent*Normalization);
    fPosProtonsSpectra->SetBinError(pt, PosProtonError*Normalization);
    fNegProtonsSpectra->SetBinContent(pt, NegProtonContent*Normalization);
    fNegProtonsSpectra->SetBinError(pt, NegProtonError*Normalization);
    
  }
}

//=======================================================================================================================
void CalculateFlow::FinalizeQA()
{
  cout << "*************************************" << "\n";
  cout << "\n";
  cout << "Finalizing the QA Pt spectra" << "\n";
  cout << "Correcting spectra for binwidth"<< "\n";
  cout << "Denk ook nog aan Nevents en bv 2pi (nu nog niet voor gecorrigeerd)";
  cout << "\n";
  cout << "\n";
  
  fPtChargedParticlesDistribution->Scale(1,"width");
  fPionsPtSpectra->Scale(1,"width");
  fPosPionsPtSpectra->Scale(1,"width");
  fAntiPionsPtSpectra->Scale(1,"width");
  fKaonsPtSpectra->Scale(1,"width");
  fPosKaonsPtSpectra->Scale(1,"width");
  fAntiKaonsPtSpectra->Scale(1,"width");
  fProtonsPtSpectra->Scale(1,"width");
  fPosProtonsPtSpectra->Scale(1,"width");
  fAntiProtonsPtSpectra->Scale(1,"width");
  
  
  
}




Int_t CalculateFlow::GetCRCCenBin(Double_t Centrality)
{
  Int_t CenBin=-1;
  if (Centrality>0. && Centrality<5.) CenBin=0;
  if (Centrality>=5. && Centrality<10.) CenBin=1;
  if (Centrality>10. && Centrality<20.) CenBin=2;
  if (Centrality>20. && Centrality<30.) CenBin=3;
  if (Centrality>30. && Centrality<40.) CenBin=4;
  if (Centrality>40. && Centrality<50.) CenBin=5;
  if (Centrality>50. && Centrality<60.) CenBin=6;
  if (Centrality>60. && Centrality<70.) CenBin=7;
  if (Centrality>70. && Centrality<80.) CenBin=8;
  if (Centrality>80. && Centrality<90.) CenBin=9;
  if (Centrality>90. && Centrality<100.) CenBin=10;
  if (CenBin>=fCRCnCen) CenBin=-1;
  if (fCRCnCen==1) CenBin=0;
  return CenBin;
} // end of CalculateFlow::GetCRCCenBin(Double_t Centrality)

//=====================================================================================================

Double_t CalculateFlow::GetSumPro(TProfile *pro, Int_t bin) //get sum of weigts in bin bin from TProfile pro
{
  Double_t stats[6]={0.};
  pro->GetXaxis()->SetRange(bin,bin);
  pro->GetStats(stats);
  pro->GetXaxis()->SetRange(1,pro->GetNbinsX());
  return stats[0];
}

//=====================================================================================================

Double_t CalculateFlow::GetWeightedCorrelations(TProfile *profile, Int_t pt)
{
  LongDouble_t CorrErr=0., Corr=0.;
  
  
  Double_t stats[6]={0.};
  profile->GetXaxis()->SetRange(pt,pt);
  profile->GetStats(stats);
  LongDouble_t SumWeig   = stats[0];      //Sum(w)
  LongDouble_t SumTwo  = stats[4];        //Sum(w*y)
  
  //std::cout<< SumWeig << "  " << SumTwo << std::endl;
  
  if(SumWeig>0.) {
    Corr = SumTwo/SumWeig;
    
  }
  // std::cout<<"corr: "<<Corr<<" CorrErr: "<<CorrErr<<std::endl;
  
  return Corr;
}

Double_t CalculateFlow::GetWeightedCorrelationsError(TProfile *profile, Int_t pt)
{
  LongDouble_t CorrErr=0., Corr=0.;
  
  
  Double_t stats[6]={0.};
  profile->GetXaxis()->SetRange(pt,pt);
  profile->GetStats(stats);
  LongDouble_t SumWeig   = stats[0];      //Sum(w)
  LongDouble_t SumWeigSq  = stats[1];     //Sum(w^2)
  LongDouble_t SumTwo  = stats[4];        //Sum(w*y)
  LongDouble_t SumTwoSq = stats[5];       //Sum(wy^2)
  
  //std::cout<< SumWeig << "  " << SumTwo << std::endl;
  
  if(SumWeig>0.) {
    Corr = SumTwo/SumWeig;
    LongDouble_t SqCorr = SumTwoSq/SumWeig;
    LongDouble_t Weig = SumWeig;
    LongDouble_t SqWeig = SumWeigSq;
    LongDouble_t spread=0., termA=0., termB=0.;
    if(SqCorr-pow(Corr,2.)>=0.) { spread = pow(SqCorr-pow(Corr,2.),0.5); }
    if(TMath::Abs(Weig)>0.) { termA = (pow(SqWeig,0.5)/Weig); }
    if(1.-pow(termA,2.)>0.) { termB = 1./pow(1.-pow(termA,2.),0.5); }
    CorrErr = termA*spread*termB; // final error (unbiased estimator for standard deviation)
    
    
  }
  
  
  return CorrErr;
}



//=====================================================================================================

std::complex<double> CalculateFlow::ucN(const Int_t n, const TArrayI& h, Int_t ptb=-1)
{
  TArrayI cnt(n);
  for (Int_t i = 0; i < n; i++) {
    cnt[i] = 1;
  }
  TArrayI hh(h);

  return ucN2(n, hh, cnt, ptb);
}

std::complex<double> CalculateFlow::ucN2(const Int_t n, TArrayI& h, TArrayI& cnt, Int_t ptb=-1)
{
  Int_t j = n-1;
  std::complex<double> c;
  if(ptb<0) {
    if(h[j] >= 0) {
      c = std::complex<double>((*fReQGF)(h[j],cnt[j]),(*fImQGF)(h[j],cnt[j]));
    } else {
      c = std::complex<double>((*fReQGF)(-h[j],cnt[j]),-(*fImQGF)(-h[j],cnt[j]));
    }
  } else {
    if(h[j] >= 0) {
      c = std::complex<double>((*fReQGFPt[ptb])(h[j],cnt[j]),(*fImQGFPt[ptb])(h[j],cnt[j]));
    } else {
      c = std::complex<double>((*fReQGFPt[ptb])(-h[j],cnt[j]),-(*fImQGFPt[ptb])(-h[j],cnt[j]));
    }
  }

  if (n == 1) return c;

  c *= ucN2(j, h, cnt, ptb);

  if (cnt[j] > 1) return c;

  for (Int_t i = 0; i < (n-1); i++) {
    h[i] += h[j];
    cnt[i] = cnt[i] + 1;
    Double_t factor = 1.*(cnt[i]-1);
    c -= factor * ucN2(j, h, cnt, ptb);
    cnt[i]--;
    h[i] -= h[j];
  }
  return c;
}
