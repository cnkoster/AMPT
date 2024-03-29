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

/// Remove Qvecpowers for EP!! And add finalize flow!

CalculateFlow::CalculateFlow ()
{
  std::cout<<"default CalculateFlow constructor"<<'\n';
  
}

//=====================================================================================================

CalculateFlow::CalculateFlow(const char* name):fQAList(NULL),fSpectraList(NULL),fFlowQCList(NULL), fFlowEPList(NULL), fFlowRPList(NULL), fFlowGFList(NULL), fFlowQCCenBin(10), fReQGF(NULL), fImQGF(NULL)
{
  std::cout<<"CalculateFlow constructor"<<'\n';
  
  fFlowQCList = new TList();
  fFlowQCList->SetName("fFlowQCList");
  fFlowQCList->SetOwner(kTRUE);
  
  fFlowEPList = new TList();
  fFlowEPList->SetName("fFlowEPList"); // Event Plane
  fFlowEPList->SetOwner(kTRUE);
  
  fFlowRPList = new TList();
  fFlowRPList->SetName("fFlowRPList"); // Reference Plane ( =0 in AMPT)
  fFlowRPList->SetOwner(kTRUE);
  
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
  // InitializeArraysForFlowGF();
//  InitializeArraysForFlowEPRP();
  //  InitializeArraysForFlowRP();
  
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
      for(Int_t p=0;p<fNParticles;p++) {
          for(Int_t charge=0; charge<fCharge;charge++){//fFlowNHarmMax
            
            for (Int_t c=0;c<fQVecPower;c++) {                  //fQVecPower
            for (Int_t h=0;h<fFlowNHarmMax;h++) {
              fPOIPtDiffQRe[c][h][p][charge] = NULL;         //POI Pt Diff Q Re [fQVecPower][fFlowHarmonic]
              fPOIPtDiffQIm[c][h][p][charge] = NULL;
              fPOIPtDiffMul[c][h][p][charge] = NULL;
            }
          }
        
        for(Int_t i=0; i<fFlowNHarm; i++) {
          for(Int_t j=0; j<fkFlowQCnIntCorPro; j++) {
            // charge=0 for all part. charge=1 for Charge+ charge=2 for Charge-
            fFlowQCIntCorPro[i][j][p][charge]= NULL;
            fFlowQCIntCorHist[i][j][p][charge] = NULL;
            fFlowQCIntFlow2Hist[i][j][p][charge] = NULL;
            fFlowQCIntFlow4Hist[i][j][p][charge] = NULL;
            fFlowQCIntCumHist[i][j][p][charge] = NULL;
          }
        }
        
        // reference flow
        for(Int_t i=0; i<fFlowNHarm; i++) {
          for(Int_t j=0; j<fFlowQCNRef; j++) {
            fFlowQCRefCorPro[i][j][p][charge] = NULL;
            fFlowQCRefCorHist[i][j][p][charge] = NULL;
          }
          for(Int_t j=0; j<4; j++) {
            fFlowQCRefCorFinal[i][j][p][charge] = NULL;
          }
        }
        
        // differential flow
        for (Int_t h=0; h<fCRCMaxnCen; h++) {
          for(Int_t i=0; i<fFlowNHarm; i++) {
            for(Int_t j=0; j<fFlowQCNPro; j++) {
              fFlowQCCorPro[h][i][j][p][charge] = NULL;
              fFlowQCCorHist[h][i][j][p][charge] = NULL;
            }
            for(Int_t k=0; k<fFlowQCNCov; k++) {
              fFlowQCCorCovPro[h][i][k][p][charge] = NULL;
              fFlowQCCorCovHist[h][i][k][p][charge] = NULL;
              fFlowQCFinalPtDifHist[h][i][k][p][charge] = NULL;
            }
          }
        }
      } // end of c
  }// end of p
  
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

void CalculateFlow::InitializeArraysForFlowEPRP()
{
  
  for (Int_t h=0;h<fFlowNHarmMax;h++) {
    
    QRe_EP[h] = 0;                  //POI Pt Diff Q Re [fQVecPower][fFlowHarmonic]
    QIm_EP[h] = 0;
    Mul_EP[h] = 0;
    
    for(Int_t p=0; p<fNParticles;p++){
      for(Int_t c=0; c<fCharge;c++){
        
        fFlowEPIntPro[h][p][c] = NULL;
        fFlowEPIntFlow2Hist[h][p][c] = NULL;
        
        fFlowEPIntPro_pos[h][p][c] = NULL;
        fFlowEPIntPro_neg[h][p][c] = NULL;
        
        fFlowRPIntPro_pos[h][p][c] = NULL;
        fFlowRPIntPro_neg[h][p][c] = NULL;
        
        fFlowRPIntFlow2Hist_pos[h][p][c] = NULL;
        fFlowRPIntFlow2Hist_neg[h][p][c] = NULL;
        fFlowRPDiffFlow2Hist[h][p][c] = NULL;
        
        fPOIDiffQRe[h][p][c] = NULL;
        fPOIDiffQIm[h][p][c] = NULL;
        fPOIDiffMul[h][p][c] = NULL;
        
        fPOIDiffQRe_pos[h][p][c] = NULL;
        fPOIDiffQIm_pos[h][p][c] = NULL;
        fPOIDiffMul_pos[h][p][c] = NULL;
        
        fPOIDiffQRe_neg[h][p][c] = NULL;
        fPOIDiffQIm_neg[h][p][c] = NULL;
        fPOIDiffMul_neg[h][p][c] = NULL;
        
      }
      
    }
  }
}
//=====================================================================================================
void CalculateFlow::InitializeArraysForFlowRP()
{
  for(Int_t c=0; c<fCharge;c++){
    for(Int_t p=0; p<fNParticles; p++){
      for (Int_t h=0;h<fFlowNHarmMax;h++) {
        
        fFlowRPIntPro_pos[h][p][c] = NULL;
        fFlowRPIntPro_neg[h][p][c] = NULL;
        
        fFlowRPIntFlow2Hist_pos[h][p][c] = NULL;
        fFlowRPIntFlow2Hist_neg[h][p][c] = NULL;
        fFlowRPDiffFlow2Hist[h][p][c] = NULL;
        
        fPOIDiffQRe[h][p][c] = NULL;
        fPOIDiffQIm[h][p][c] = NULL;
        fPOIDiffMul[h][p][c] = NULL;
        
        fPOIDiffQRe_pos[h][p][c] = NULL;
        fPOIDiffQIm_pos[h][p][c] = NULL;
        fPOIDiffMul_pos[h][p][c] = NULL;
        
        fPOIDiffQRe_neg[h][p][c] = NULL;
        fPOIDiffQIm_neg[h][p][c] = NULL;
        fPOIDiffMul_neg[h][p][c] = NULL;
      }
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
  fPDGSpectra = new TH1D("fPDGSpectra", "fPDGSpectra", 8000, -4000, 4000);

  
  fSpectraList->Add(fChargedParticleSpectra);
  fSpectraList->Add(fPosPionsSpectra);
  fSpectraList->Add(fPosKaonsSpectra);
  fSpectraList->Add(fPosProtonsSpectra);
  fSpectraList->Add(fNegPionsSpectra);
  fSpectraList->Add(fNegKaonsSpectra);
  fSpectraList->Add(fNegProtonsSpectra);
  fSpectraList->Add(fPDGSpectra);
  
  // Calculate FlowQC Hists
  // choose for eta diff or pt diff:
  
  fPtDiffNBins = 36-8; //for pt 20 for eta
  fEtaDiffNBins = 5;
  
  fCRCPtBins = new Double_t[37-8];
  Double_t PtBins[] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.25,1.5,1.75,2.,2.25,2.5,2.75,3.,3.25,3.5,3.75,4.,4.5,5.,5.5,6.,7.,8.,9.,10.};//,12.,14.,17.,20.,25.,30.,40.,50.};
  Double_t EtaBins[] = {-0.8,-0.48,-0.16,0.16,0.48,0.8};
  Double_t ImPaBins[] = {3.72, 5.23, 7.31, 8.88, 10.20, 11.38,12.47, 13.50, 14.51, 15.0}; // 
  
  
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
  
  
  
  for (Int_t p=0;p<fNParticles;p++) {
    for(Int_t charge=0; charge<fCharge; charge++){
      for (Int_t c=0;c<fQVecPower;c++) {
        for (Int_t h=0;h<fFlowNHarmMax;h++) {
          fPOIPtDiffQRe[c][h][p][charge] = new TH1D(Form("fPOIPtDiffQRe[%d][%d][%d][%d]",c,h,p,charge),Form("fPOIPtDiffQRe[%d][%d][%d][%d]",c,h,p,charge), fNBins, fBins);
          //   fFlowQCList->Add(fPOIPtDiffQRe[c][h]);
          fPOIPtDiffQIm[c][h][p][charge] = new TH1D(Form("fPOIPtDiffQIm[%d][%d][%d][%d]",c,h,p,charge),Form("fPOIPtDiffQIm[%d][%d][%d][%d]",c,h,p,charge), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
          //  fFlowQCList->Add(fPOIPtDiffQIm[c][h]);
          fPOIPtDiffMul[c][h][p][charge] = new TH1D(Form("fPOIPtDiffMul[%d][%d][%d][%d]",c,h,p,charge),Form("fPOIPtDiffMul[%d][%d][%d][%d]",c,h,p,charge), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
          //   fFlowQCList->Add(fPOIPtDiffMul[c][h]);
        }
      }
      
      
      for(Int_t i=0; i<fFlowNHarm; i++) {
        for(Int_t j=0; j<fkFlowQCnIntCorPro; j++) {
          fFlowQCIntCorPro[i][j][p][charge] = new TProfile(Form("fFlowQCIntCorPro[%d][%d][%d][%d]",i,j,p,charge),Form("fFlowQCIntCorPro[%d][%d][%d][%d]",i,j,p,charge),9-fCbins,ImPaBins,"s"); //here we changed the bins to nParticipants
          fFlowQCIntCorPro[i][j][p][charge]->Sumw2();
          //fFlowQCList->Add(fFlowQCIntCorPro[i][j]);
          fFlowQCIntCorHist[i][j][p][charge] = new TH1D(Form("fFlowQCIntCorHist[%d][%d][%d][%d]",i,j,p,charge),Form("fFlowQCIntCorHist[%d][%d][%d][%d]",i,j,p,charge),9-fCbins,ImPaBins);
          fFlowQCIntCorHist[i][j][p][charge]->Sumw2();
          //   fFlowQCList->Add(fFlowQCIntCorHist[i][j][p][charge]);
          
          fFlowQCIntFlow2Hist[i][j][p][charge] = new TH1D(Form("fFlowQCIntFlow2Hist[%d][%d][%d][%d]",i,j,p,charge),Form("fFlowQCIntFlow2Hist[%d][%d][%d][%d]",i,j,p,charge),9-fCbins,ImPaBins);
          fFlowQCIntFlow2Hist[i][j][p][charge]->Sumw2();
          
          
          fFlowQCIntFlow4Hist[i][j][p][charge] = new
          TH1D(Form("fFlowQCIntFlow4Hist[%d][%d][%d][%d]",i,j,p,charge),Form("fFlowQCIntFlow4Hist[%d][%d][%d][%d]",i,j,p,charge),9-fCbins,ImPaBins);
          fFlowQCIntFlow4Hist[i][j][p][charge]->Sumw2();
          
          if(j==0){
            fFlowQCList->Add(fFlowQCIntFlow2Hist[i][j][p][charge]);}
          //        fFlowQCList->Add(fFlowQCIntFlow4Hist[i][j][p][charge]);}
          
          fFlowQCIntCumHist[i][j][p][charge] = new TH1D(Form("fFlowQCIntCumHist[%d][%d][%d][%d]",i,j,p,charge),Form("fFlowQCIntCumHist[%d][%d][%d][%d]",i,j,p,charge),9-fCbins,ImPaBins);
          fFlowQCIntCumHist[i][j][p][charge]->Sumw2();
          //   fFlowQCList->Add(fFlowQCIntCumHist[i][j][p][charge]);
        }
      }
      
      
      // reference flow //Error kan hier komen
      for(Int_t i=0; i<fFlowNHarm; i++) {
        for(Int_t j=0; j<fFlowQCNRef; j++) {
          fFlowQCRefCorPro[i][j][p][charge] = new TProfile(Form("fFlowQCRefCorPro[%d][%d][%d][%d]",i,j,p,charge),Form("fFlowQCRefCorPro[%d][%d][%d][%d]",i,j,p,charge),9-fCbins,ImPaBins,"s");
          fFlowQCRefCorPro[i][j][p][charge]->Sumw2();
          //fFlowQCList->Add(fFlowQCRefCorPro[i][j][p][charge]);
          fFlowQCRefCorHist[i][j][p][charge] = new TH1D(Form("fFlowQCRefCorHist[%d][%d][%d][%d]",i,j,p,charge),Form("fFlowQCRefCorHist[%d][%d][%d][%d]",i,j,p,charge),9-fCbins,ImPaBins);
          fFlowQCRefCorHist[i][j][p][charge]->Sumw2();
          // fFlowQCList->Add(fFlowQCRefCorHist[i][j][p][charge]);
  
          fFlowQCRefCorFinal[i][j][p][charge] = new TH1D(Form("fFlowQCRefCorFinal[%d][%d][%d][%d]",i,j,p,charge),Form("fFlowQCRefCorFinal[%d][%d][%d][%d]",i,j,p,charge),9-fCbins,ImPaBins);
          fFlowQCRefCorFinal[i][j][p][charge]->Sumw2();
          //  fFlowQCList->Add(fFlowQCRefCorFinal[i][j][p][charge]);
        }
      }
      
      
      // differential flow
      for (Int_t h=0; h<fCRCnCen; h++) {
        for(Int_t i=0; i<fFlowNHarm; i++) {
          for(Int_t j=0; j<fFlowQCNPro; j++) {
            fFlowQCCorPro[h][i][j][p][charge] = new TProfile(Form("fFlowQCCorPro[%d][%d][%d][%d][%d]",h,i,j,p,charge),Form("fFlowQCCorPro[%d][%d][%d][%d][%d]",h,i,j,p,charge), fNBins, fBins,"s");
            fFlowQCCorPro[h][i][j][p][charge]->Sumw2();
            //      fFlowQCList->Add(fFlowQCCorPro[h][i][j]);
            fFlowQCCorHist[h][i][j][p][charge] = new TH1D(Form("fFlowQCCorHist[%d][%d][%d][%d][%d]",h,i,j,p,charge),Form("fFlowQCCorHist[%d][%d][%d][%d][%d]",h,i,j,p,charge), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
            fFlowQCCorHist[h][i][j][p][charge]->Sumw2();
            //  fFlowQCList->Add(fFlowQCCorHist[h][i][j]);
          }
          
          
          for(Int_t k=0; k<fFlowQCNCov; k++) {
            fFlowQCCorCovPro[h][i][k][p][charge] = new TProfile(Form("fFlowQCCorCovPro[%d][%d][%d][%d][%d]",h,i,k,p,charge),Form("fFlowQCCorCovPro[%d][%d][%d][%d][%d]",h,i,k,p,charge), fNBins, fBins,"s");
            fFlowQCCorCovPro[h][i][k][p][charge]->Sumw2();
            //  fFlowQCList->Add(fFlowQCCorCovPro[h][i][k]);
            fFlowQCCorCovHist[h][i][k][p][charge] = new TH1D(Form("fFlowQCCorCovHist[%d][%d][%d][%d][%d]",h,i,k,p,charge),Form("fFlowQCCorCovHist[%d][%d][%d][%d][%d]",h,i,k,p,charge), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
            fFlowQCCorCovHist[h][i][k][p][charge]->Sumw2();
            //    fFlowQCList->Add(fFlowQCCorCovHist[h][i][k]);
            fFlowQCFinalPtDifHist[h][i][k][p][charge] = new TH1D(Form("fFlowQCFinalPtDifHist[%d][%d][%d][%d][%d]",h,i,k,p,charge),Form("fFlowQCFinalPtDifHist[%d][%d][%d][%d][%d]",h,i,k,p,charge), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
            fFlowQCFinalPtDifHist[h][i][k][p][charge]->Sumw2();
            
            if(h==GetCRCCenBin(fCentralityEBE)){
              if(k==0 || k==1){fFlowQCList->Add(fFlowQCFinalPtDifHist[h][i][k][p][charge]);} //Only print for v{2} and v{4} [k] k=5,6 give Dn2, Dn4
            }
          }
        }
      }
      
    } //end of charge
  } //end of particles
    
  
  // for CalculateFlowGF
  for (Int_t h=0; h<fkFlowGFNHarm; h++) {
    for(Int_t i=0; i<fkFlowGFNOrde; i++) {
      fFlowGFIntCorPro[h][i] = new TProfile(Form("fFlowGFIntCorPro[%d][%d]",h,i),Form("fFlowGFIntCorPro[%d][%d]",h,i),9-fCbins,ImPaBins,"s");//fFlowGFCenBin,0.,100.,"s");
      fFlowGFIntCorPro[h][i]->Sumw2();
      //      fFlowGFList->Add(fFlowGFIntCorPro[h][i]);
      fFlowGFIntCorHist[h][i] = new TH1D(Form("fFlowGFIntCorHist[%d][%d]",h,i),Form("fFlowGFIntCorHist[%d][%d]",h,i),9-fCbins,ImPaBins);//fFlowGFCenBin,0.,100.);
      fFlowGFIntCorHist[h][i]->Sumw2();
      //     fFlowGFList->Add(fFlowGFIntCorHist[h][i]);
      fFlowGFIntCumHist[h][i] = new TH1D(Form("fFlowGFIntCumHist[%d][%d]",h,i),Form("fFlowGFIntCumHist[%d][%d]",h,i),9-fCbins,ImPaBins);//fFlowGFCenBin,0.,100.);
      fFlowGFIntCumHist[h][i]->Sumw2();
      //         fFlowGFList->Add(fFlowGFIntCumHist[h][i]);
      fFlowGFIntFinalHist[h][i] = new TH1D(Form("fFlowGFIntFinalHist[%d][%d]",h,i),Form("fFlowGFIntFinalHist[%d][%d]",h,i),9-fCbins,ImPaBins);//fFlowGFCenBin,0.,100.);
      fFlowGFIntFinalHist[h][i]->Sumw2();
      fFlowGFList->Add(fFlowGFIntFinalHist[h][i]);
      for(Int_t k=0; k<fkFlowGFNOrde; k++) {
        fFlowGFIntCovPro[h][i][k] = new TProfile(Form("fFlowGFIntCovPro[%d][%d][%d]",h,i,k),Form("fFlowGFIntCovPro[%d][%d][%d]",h,i,k),9-fCbins,ImPaBins,"s");//fFlowGFCenBin,0.,100.,"s");
        fFlowGFIntCovPro[h][i][k]->Sumw2();
        //           fFlowGFList->Add(fFlowGFIntCovPro[h][i][k]);
        fFlowGFIntCovHist[h][i][k] = new TH1D(Form("fFlowGFIntCovHist[%d][%d][%d]",h,i,k),Form("fFlowGFIntCovHist[%d][%d][%d]",h,i,k),9-fCbins,ImPaBins);//fFlowGFCenBin,0.,100.);
        fFlowGFIntCovHist[h][i][k]->Sumw2();
        //          fFlowGFList->Add(fFlowGFIntCovHist[h][i][k]);
      }
      for(Int_t s=0; s<fkGFPtB; s++) {
        fFlowGFIntCorProPtB[s][h][i] = new TProfile(Form("fFlowGFIntCorProPtB[%d][%d][%d]",s,h,i),Form("fFlowGFIntCorProPtB[%d][%d][%d]",s,h,i),9-fCbins,ImPaBins,"s");//fFlowGFCenBin,0.,100.,"s");
        fFlowGFIntCorProPtB[s][h][i]->Sumw2();
        //            fFlowGFList->Add(fFlowGFIntCorProPtB[s][h][i]);
        fFlowGFIntCorHistPtB[s][h][i] = new TH1D(Form("fFlowGFIntCorHistPtB[%d][%d][%d]",s,h,i),Form("fFlowGFIntCorHistPtB[%d][%d][%d]",s,h,i),9-fCbins,ImPaBins);//fFlowGFCenBin,0.,100.);
        fFlowGFIntCorHistPtB[s][h][i]->Sumw2();
        //            fFlowGFList->Add(fFlowGFIntCorHistPtB[s][h][i]);
        for(Int_t k=0; k<fkFlowGFNOrde; k++) {
          fFlowGFIntCovProPtB[s][h][i][k] = new TProfile(Form("fFlowGFIntCovProPtB[%d][%d][%d][%d]",s,h,i,k),Form("fFlowGFIntCovProPtB[%d][%d][%d][%d]",s,h,i,k),9-fCbins,ImPaBins,"s");//fFlowGFCenBin,0.,100.,"s");
          fFlowGFIntCovProPtB[s][h][i][k]->Sumw2();
          //  fFlowGFList->Add(fFlowGFIntCovProPtB[s][h][i][k]);
          fFlowGFIntCovHistPtB[s][h][i][k] = new TH1D(Form("fFlowGFIntCovHistPtB[%d][%d][%d][%d]",s,h,i,k),Form("fFlowGFIntCovHistPtB[%d][%d][%d][%d]",s,h,i,k),9-fCbins,ImPaBins);//fFlowGFCenBin,0.,100.);
          fFlowGFIntCovHistPtB[s][h][i][k]->Sumw2();
          //fFlowGFList->Add(fFlowGFIntCovHistPtB[s][h][i][k]);
        }
      }
    }
  }
  
  for (Int_t h=0; h<fkFlowGFNHarm; h++) {
    for(Int_t i=0; i<fkFlowGFNHarm; i++) {
      fFlowGFMixedCorPro[h][i] = new TProfile(Form("fFlowGFMixedCorPro[%d][%d]",h,i),Form("fFlowGFMixedCorPro[%d][%d]",h,i),9-fCbins,ImPaBins);//fFlowGFCenBin,0.,100.,"s");
      fFlowGFMixedCorPro[h][i]->Sumw2();
      //  fFlowGFList->Add(fFlowGFMixedCorPro[h][i]);
      fFlowGFMixedCorHist[h][i] = new TH1D(Form("fFlowGFMixedCorHist[%d][%d]",h,i),Form("fFlowGFMixedCorHist[%d][%d]",h,i),9-fCbins,ImPaBins);//fFlowGFCenBin,0.,100.);
      fFlowGFMixedCorHist[h][i]->Sumw2();
      //   fFlowGFList->Add(fFlowGFMixedCorHist[h][i]);
      fFlowGFMixedFinalHist[h][i] = new TH1D(Form("fFlowGFMixedFinalHist[%d][%d]",h,i),Form("fFlowGFMixedFinalHist[%d][%d]",h,i),9-fCbins,ImPaBins);//fFlowGFCenBin,0.,100.);
      fFlowGFMixedFinalHist[h][i]->Sumw2();
      //   fFlowGFList->Add(fFlowGFMixedFinalHist[h][i]);
    }
  }
  
  // Calculate FlowEP Hists --------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // so for each power of the Q-vector and each flowHarmonic make a Th1D
  // In output vinden we: [0][0] tot [13][13] ?
  
  
  
  for (Int_t h=0;h<fFlowNHarm;h++) {
    fEPEPresolutionPro[h] = new TH1D(Form("fEPEPresolutionPro[%d]",h),Form("fEPEPresolutionPro[%d]",h),1000,-TMath::Pi()-0.2,TMath::Pi()+0.2);
    //fEPEPresolutionPro[h]->Sumw2();
    fFlowEPList->Add(fEPEPresolutionPro[h]);
    
    for(Int_t p=0; p<fNParticles;p++){
      for(Int_t c=0; c<fCharge; c++){
        
        fPOIDiffQRe[h][p][c] = new TH1D(Form("fPOIDiffQRe[%d][%d][%d]",h,p,c),Form("fPOIDiffQRe[%d][%d][%d]",h,p,c), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
        fPOIDiffQRe[h][p][c]->Sumw2();
        fPOIDiffQIm[h][p][c] = new TH1D(Form("fPOIDiffQIm[%d][%d][%d]",h,p,c),Form("fPOIDiffQIm[%d][%d][%d]",h,p,c), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
        fPOIDiffQIm[h][p][c]->Sumw2();
        fPOIDiffMul[h][p][c] = new TH1D(Form("fPOIDiffMul[%d][%d][%d]",h,p,c),Form("fPOIDiffMul[%d][%d][%d]",h,p,c), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
        fPOIDiffMul[h][p][c]->Sumw2();
        
        fPOIDiffQRe_pos[h][p][c] = new TH1D(Form("fPOIDiffQRe_pos[%d][%d][%d]",h,p,c),Form("fPOIDiffQRe_pos[%d][%d][%d]",h,p,c), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
        fPOIDiffQRe_pos[h][p][c]->Sumw2();
        fPOIDiffQIm_pos[h][p][c] = new TH1D(Form("fPOIDiffQIm_pos[%d][%d][%d]",h,p,c),Form("fPOIDiffQIm_pos[%d][%d][%d]",h,p,c), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
        fPOIDiffQIm_pos[h][p][c]->Sumw2();
        fPOIDiffMul_pos[h][p][c] = new TH1D(Form("fPOIDiffMul_pos[%d][%d][%d]",h,p,c),Form("fPOIDiffMul_pos[%d][%d][%d]",h,p,c), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
        fPOIDiffMul_pos[h][p][c]->Sumw2();
        
        
        fPOIDiffQRe_neg[h][p][c] = new TH1D(Form("fPOIDiffQRe_neg[%d][%d][%d]",h,p,c),Form("fPOIDiffQRe_neg[%d][%d][%d]",h,p,c), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
        fPOIDiffQRe_neg[h][p][c]->Sumw2();
        fPOIDiffQIm_neg[h][p][c] = new TH1D(Form("fPOIDiffQIm_neg[%d][%d][%d]",h,p,c),Form("fPOIDiffQIm_neg[%d][%d][%d]",h,p,c), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
        fPOIDiffQIm_neg[h][p][c]->Sumw2();
        fPOIDiffMul_neg[h][p][c] = new TH1D(Form("fPOIDiffMul_neg[%d][%d][%d]",h,p,c),Form("fPOIDiffMul_neg[%d][%d][%d]",h,p,c), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
        fPOIDiffMul_neg[h][p][c]->Sumw2();
        
        
        fFlowEPIntPro[h][p][c]= new TProfile(Form("fFlowEPIntPro[%d][%d][%d]",h,p,c),Form("fFlowEPIntPro[%d][%d][%d]",h,p,c),9-fCbins,ImPaBins,"s");
        fFlowEPIntPro[h][p][c]->Sumw2();
        
        fFlowEPIntPro_neg[h][p][c]= new TProfile(Form("fFlowEPIntPro_neg[%d][%d][%d]",h,p,c),Form("fFlowEPIntPro_neg[%d][%d][%d]",h,p,c),9-fCbins,ImPaBins,"s");
        fFlowEPIntPro_neg[h][p][c]->Sumw2();
        
        fFlowEPIntPro_pos[h][p][c]= new TProfile(Form("fFlowEPIntPro_pos[%d][%d][%d]",h,p,c),Form("fFlowEPIntPro_pos[%d][%d][%d]",h,p,c),9-fCbins,ImPaBins,"s");
        fFlowEPIntPro_pos[h][p][c]->Sumw2();
        
        
        fFlowEPIntFlow2Hist[h][p][c] = new TH1D(Form("fFlowEPIntFlow2Hist[%d][%d][%d]",h,p,c),Form("fFlowEPIntFlow2Hist[%d][%d][%d]",h,p,c),9-fCbins,ImPaBins);
        fFlowEPIntFlow2Hist[h][p][c]->Sumw2();
        //      fFlowEPList->Add(fFlowEPIntFlow2Hist[h][p][c]);
        
        fFlowEPIntFlow2Hist_pos[h][p][c] = new TH1D(Form("fFlowEPIntFlow2Hist_pos[%d][%d][%d]",h,p,c),Form("fFlowEPIntFlow2Hist_pos[%d][%d][%d]",h,p,c),9-fCbins,ImPaBins);
        fFlowEPIntFlow2Hist_pos[h][p][c]->Sumw2();
        //      fFlowEPList->Add(fFlowEPIntFlow2Hist[h][p][c]);
        
        fFlowEPIntFlow2Hist_neg[h][p][c] = new TH1D(Form("fFlowEPIntFlow2Hist_neg[%d][%d][%d]",h,p,c),Form("fFlowEPIntFlow2Hist_neg[%d][%d][%d]",h,p,c),9-fCbins,ImPaBins);
        fFlowEPIntFlow2Hist_neg[h][p][c]->Sumw2();
        //      fFlowEPList->Add(fFlowEPIntFlow2Hist[h][p][c]);
        
        fFlowEPCorPro[h][p][c]= new TProfile(Form("fFlowEPCorPro[%d][%d][%d]",h,p,c),Form("fFlowEPCorPro[%d][%d][%d]",h,p,c), fNBins, fBins, "s");//, fNBins, fBins,"s");
        fFlowEPCorPro[h][p][c]->Sumw2();
        
        fFlowEPDiffFlow2Hist[h][p][c] = new TH1D(Form("fFlowEPDiffFlow2Hist[%d][%d][%d]",h,p,c),Form("fFlowEPDiffFlow2Hist[%d][%d][%d]",h,p,c), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
        fFlowEPDiffFlow2Hist[h][p][c]->Sumw2();
        fFlowEPList->Add(fFlowEPDiffFlow2Hist[h][p][c]);
        
      }//end of c
    }//end of p
  }//end of h
  
  // COS(NPHI)
  for (Int_t h=0;h<fFlowNHarm;h++) {
    for(Int_t p=0; p<fNParticles; p++){
      for(Int_t c=0; c<fCharge; c++){
        // Event Plane Method (SPECTATORS +EPangle=0)
        
        fFlowRPIntPro_pos[h][p][c]= new TProfile(Form("fFlowRPIntPro_pos[%d][%d][%d]",h,p,c),Form("fFlowRPIntPro_pos[%d][%d][%d]",h,p,c),9-fCbins,ImPaBins,"s");
        fFlowRPIntPro_pos[h][p][c]->Sumw2();
        
        fFlowRPIntPro_neg[h][p][c]= new TProfile(Form("fFlowRPIntPro_neg[%d][%d][%d]",h,p,c),Form("fFlowRPIntPro_neg[%d][%d][%d]",h,p,c),9-fCbins,ImPaBins,"s");
        fFlowRPIntPro_neg[h][p][c]->Sumw2();
        
        fFlowRPIntPro[h][p][c]= new TProfile(Form("fFlowRPIntPro[%d][%d][%d]",h,p,c),Form("fFlowRPIntPro[%d][%d][%d]",h,p,c),9-fCbins,ImPaBins,"s");
        fFlowRPIntPro[h][p][c]->Sumw2();
        
        fFlowRPIntFlow2Hist_pos[h][p][c] = new TH1D(Form("fFlowRPIntFlow2Hist_pos[%d][%d][%d]",h,p,c),Form("fFlowRPIntFlow2Hist_pos[%d][%d][%d]",h,p,c),9-fCbins,ImPaBins); //,fPtDiffNBins,fCRCPtBins);
        fFlowRPIntFlow2Hist_pos[h][p][c]->Sumw2();
        //      fFlowEPList->Add(fFlowRPIntFlow2Hist_pos[h][particle][charge]);
        
        fFlowRPIntFlow2Hist_neg[h][p][c] = new TH1D(Form("fFlowRPIntFlow2Hist_neg[%d][%d][%d]",h,p,c),Form("fFlowRPIntFlow2Hist_neg[%d][%d][%d]",h,p,c),9-fCbins,ImPaBins); //,fPtDiffNBins,fCRCPtBins);
        fFlowRPIntFlow2Hist_neg[h][p][c]->Sumw2();
        //      fFlowEPList->Add(fFlowRPIntFlow2Hist_neg[h][p][c]);
        
        fFlowRPIntFlow2Hist[h][p][c] = new TH1D(Form("fFlowRPIntFlow2Hist[%d][%d][%d]",h,p,c),Form("fFlowRPIntFlow2Hist[%d][%d][%d]",h,p,c),9-fCbins,ImPaBins); //,fPtDiffNBins,fCRCPtBins);
        fFlowRPIntFlow2Hist[h][p][c]->Sumw2();
        //      fFlowEPList->Add(fFlowRPIntFlow2Hist[h][p][c]);
        
        fFlowRPCorPro[h][p][c]= new TProfile(Form("fFlowRPCorPro[%d][%d][%d]",h,p,c),Form("fFlowRPCorPro[%d][%d][%d]",h,p,c), fNBins, fBins, "s");//, fNBins, fBins,"s");
        fFlowRPCorPro[h][p][c]->Sumw2();
        
        fFlowRPDiffFlow2Hist[h][p][c] = new TH1D(Form("fFlowRPDiffFlow2Hist[%d][%d][%d]",h,p,c),Form("fFlowRPDiffFlow2Hist[%d][%d][%d]",h,p,c), fNBins, fBins);//,fPtDiffNBins,fCRCPtBins);
        fFlowRPDiffFlow2Hist[h][p][c]->Sumw2();
        fFlowRPList->Add(fFlowRPDiffFlow2Hist[h][p][c]);
        
      }//end of charge
    }//end of particle
  }//end of h
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
//  std::cout<< " Centrality bin is: "<< fCenBin <<std::endl;
  nPrim = anEvent->getNtrack();
  fNumberOfParticipants = anEvent->getnPart();
  fImpactParameter = anEvent->getb();
  Int_t cw = 0;
  Int_t nChargedParticles = 0;
  // multiplicity for charged particles
  Double_t fNumOfPos = 0;
  Double_t fNumOfNeg = 0;
  Double_t xval; // used for pt or eta diff flow
  
  
  if (nPrim < minNtracks) return;
  sw.tick();
  for(Int_t i=0;i<nPrim;i++) {
    if(anEvent->getParticle(i).getSpecflag()<0) continue;
    Pid = anEvent->getParticle(i).getPid();
    dPhi = anEvent->getParticle(i).getPhi();
    dPt = anEvent->getParticle(i).getPt();
    dEta = anEvent->getParticle(i).getRapidity();  // in MC, the rapidity is stored rather than pesudorapidity (eta)
    dCharge = anEvent->getParticle(i).getCharge();
  
    if(TMath::Abs(Pid)<23 && Pid>-23)continue;
    
    for (Int_t h=0;h<4;h++) { // let op: vorige run was dit NA pt cut (4/12/23)
      QRe_EP[h] += pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi);
      QIm_EP[h] += pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi);
      Mul_EP[h] += pow(wPhiEta,1);
    }
    
    
    if (dPt > maxPtCut) continue;
    if (dPt < minPtCut) continue;
    
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
      fPDGSpectra->Fill(Pid); //check what particles you include !!
      
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
    
//
//    for (Int_t k=0; k<fQVecPower; k++) {
//      for (Int_t h=0;h<fFlowNHarmMax;h++) {
//        fPOIPtDiffQRe[k][h][0][0]->Fill(xval,pow(wPhiEta,k)*TMath::Cos((h+1.)*dPhi));
//        fPOIPtDiffQIm[k][h][0][0]->Fill(xval,pow(wPhiEta,k)*TMath::Sin((h+1.)*dPhi));
//        fPOIPtDiffMul[k][h][0][0]->Fill(xval,pow(wPhiEta,k));
//
//
//        if(dCharge>0){
//          //if(dPhi>TMath::Pi)
//          fPOIPtDiffQRe[k][h][0][1]->Fill(xval,pow(wPhiEta,k)*TMath::Cos((h+1.)*dPhi));
//          fPOIPtDiffQIm[k][h][0][1]->Fill(xval,pow(wPhiEta,k)*TMath::Sin((h+1.)*dPhi));
//          fPOIPtDiffMul[k][h][0][1]->Fill(xval,pow(wPhiEta,k));
//
//        }
//
//
//        if(dCharge<0){
//          fPOIPtDiffQRe[k][h][0][2]->Fill(xval,pow(wPhiEta,k)*TMath::Cos((h+1.)*dPhi));
//          fPOIPtDiffQIm[k][h][0][2]->Fill(xval,pow(wPhiEta,k)*TMath::Sin((h+1.)*dPhi));
//          fPOIPtDiffMul[k][h][0][2]->Fill(xval,pow(wPhiEta,k));
//        }
//
//        //PIONS
//        if (Pid == 211 || Pid == -211) {
//          fPOIPtDiffQRe[k][h][1][0]->Fill(xval,pow(wPhiEta,k)*TMath::Cos((h+1.)*dPhi));
//          fPOIPtDiffQIm[k][h][1][0]->Fill(xval,pow(wPhiEta,k)*TMath::Sin((h+1.)*dPhi));
//          fPOIPtDiffMul[k][h][1][0]->Fill(xval,pow(wPhiEta,k));
//
//          if(dCharge>0){
//            //if(dPhi>TMath::Pi)
//            fPOIPtDiffQRe[k][h][1][1]->Fill(xval,pow(wPhiEta,k)*TMath::Cos((h+1.)*dPhi));
//            fPOIPtDiffQIm[k][h][1][1]->Fill(xval,pow(wPhiEta,k)*TMath::Sin((h+1.)*dPhi));
//            fPOIPtDiffMul[k][h][1][1]->Fill(xval,pow(wPhiEta,k));
//
//          }
//          if(dCharge<0){
//            fPOIPtDiffQRe[k][h][1][2]->Fill(xval,pow(wPhiEta,k)*TMath::Cos((h+1.)*dPhi));
//            fPOIPtDiffQIm[k][h][1][2]->Fill(xval,pow(wPhiEta,k)*TMath::Sin((h+1.)*dPhi));
//            fPOIPtDiffMul[k][h][1][2]->Fill(xval,pow(wPhiEta,k));
//          }
//        }
//
//        //KAONS
//        if (Pid == 321 || Pid == -321) {
//          fPOIPtDiffQRe[k][h][2][0]->Fill(xval,pow(wPhiEta,k)*TMath::Cos((h+1.)*dPhi));
//          fPOIPtDiffQIm[k][h][2][0]->Fill(xval,pow(wPhiEta,k)*TMath::Sin((h+1.)*dPhi));
//          fPOIPtDiffMul[k][h][2][0]->Fill(xval,pow(wPhiEta,k));
//
//          if(dCharge>0){
//            //if(dPhi>TMath::Pi)
//            fPOIPtDiffQRe[k][h][2][1]->Fill(xval,pow(wPhiEta,k)*TMath::Cos((h+1.)*dPhi));
//            fPOIPtDiffQIm[k][h][2][1]->Fill(xval,pow(wPhiEta,k)*TMath::Sin((h+1.)*dPhi));
//            fPOIPtDiffMul[k][h][2][1]->Fill(xval,pow(wPhiEta,k));
//
//          }
//          if(dCharge<0){
//            fPOIPtDiffQRe[k][h][2][2]->Fill(xval,pow(wPhiEta,k)*TMath::Cos((h+1.)*dPhi));
//            fPOIPtDiffQIm[k][h][2][2]->Fill(xval,pow(wPhiEta,k)*TMath::Sin((h+1.)*dPhi));
//            fPOIPtDiffMul[k][h][2][2]->Fill(xval,pow(wPhiEta,k));
//          }
//        }
//
//        //PROTONS
//        if (Pid == 2212 || Pid == -2212) {
//          fPOIPtDiffQRe[k][h][3][0]->Fill(xval,pow(wPhiEta,k)*TMath::Cos((h+1.)*dPhi));
//          fPOIPtDiffQIm[k][h][3][0]->Fill(xval,pow(wPhiEta,k)*TMath::Sin((h+1.)*dPhi));
//          fPOIPtDiffMul[k][h][3][0]->Fill(xval,pow(wPhiEta,k));
//
//          if(dCharge>0){
//            //if(dPhi>TMath::Pi)
//            fPOIPtDiffQRe[k][h][3][1]->Fill(xval,pow(wPhiEta,k)*TMath::Cos((h+1.)*dPhi));
//            fPOIPtDiffQIm[k][h][3][1]->Fill(xval,pow(wPhiEta,k)*TMath::Sin((h+1.)*dPhi));
//            fPOIPtDiffMul[k][h][3][1]->Fill(xval,pow(wPhiEta,k));
//
//          }
//          if(dCharge<0){
//            fPOIPtDiffQRe[k][h][3][2]->Fill(xval,pow(wPhiEta,k)*TMath::Cos((h+1.)*dPhi));
//            fPOIPtDiffQIm[k][h][3][2]->Fill(xval,pow(wPhiEta,k)*TMath::Sin((h+1.)*dPhi));
//            fPOIPtDiffMul[k][h][3][2]->Fill(xval,pow(wPhiEta,k));
//          }
//        }
//
//      } // end of h loop
//    } // end of k loop
//
//
//    // ====== for calculateFlowEP/RP =========
//
//    // All hadrons +/-
//    for (Int_t h=0;h<4;h++) {
//
//      fPOIDiffQRe[h][0][0]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//      fPOIDiffQIm[h][0][0]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//      fPOIDiffMul[h][0][0]->Fill(xval,pow(wPhiEta,1));
//
//      if(dCharge>0){
//        //if(dPhi>TMath::Pi)
//        fPOIDiffQRe[h][0][1]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//        fPOIDiffQIm[h][0][1]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//        fPOIDiffMul[h][0][1]->Fill(xval,pow(wPhiEta,1));
//
//      }
//      if(dCharge<0){
//        fPOIDiffQRe[h][0][2]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//        fPOIDiffQIm[h][0][2]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//        fPOIDiffMul[h][0][2]->Fill(xval,pow(wPhiEta,1));
//      }
//
//
//
//      if(dEta>=0){
//        for (Int_t h=0;h<4;h++) {
//          fPOIDiffQRe_pos[h][0][0]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//          fPOIDiffQIm_pos[h][0][0]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//          fPOIDiffMul_pos[h][0][0]->Fill(xval,pow(wPhiEta,1));
//
//
//
//          if(dCharge>0){
//            //if(dPhi>TMath::Pi)
//            fPOIDiffQRe_pos[h][0][1]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//            fPOIDiffQIm_pos[h][0][1]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//            fPOIDiffMul_pos[h][0][1]->Fill(xval,pow(wPhiEta,1));
//
//          }
//
//
//          if(dCharge<0){
//            fPOIDiffQRe_pos[h][0][2]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//            fPOIDiffQIm_pos[h][0][2]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//            fPOIDiffMul_pos[h][0][2]->Fill(xval,pow(wPhiEta,1));}
//
//        }
//      }
//
//      if(dEta<0){ //check even: moet die min hier?
//        for (Int_t h=0;h<4;h++) {
//          fPOIDiffQRe_neg[h][0][0]->Fill(xval,pow(-wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//          fPOIDiffQIm_neg[h][0][0]->Fill(xval,pow(-wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//          fPOIDiffMul_neg[h][0][0]->Fill(xval,pow(-wPhiEta,1));
//
//
//
//          if(dCharge>0){
//            //if(dPhi>TMath::Pi)
//            fPOIDiffQRe_neg[h][0][1]->Fill(xval,pow(-wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//            fPOIDiffQIm_neg[h][0][1]->Fill(xval,pow(-wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//            fPOIDiffMul_neg[h][0][1]->Fill(xval,pow(-wPhiEta,1));
//
//          }
//
//
//          if(dCharge<0){
//            fPOIDiffQRe_neg[h][0][2]->Fill(xval,pow(-wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//            fPOIDiffQIm_neg[h][0][2]->Fill(xval,pow(-wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//            fPOIDiffMul_neg[h][0][2]->Fill(xval,pow(-wPhiEta,1));
//
//          }
//
//        }
//      }
//
//
//      //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
//      //PID
//
//      // Pions
//      if (Pid == 211 || Pid == -211) {
//        fPOIDiffQRe[h][1][0]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//        fPOIDiffQIm[h][1][0]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//        fPOIDiffMul[h][1][0]->Fill(xval,pow(wPhiEta,1));
//        if(dCharge>0){
//          //if(dPhi>TMath::Pi)
//          fPOIDiffQRe[h][1][1]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//          fPOIDiffQIm[h][1][1]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//          fPOIDiffMul[h][1][1]->Fill(xval,pow(wPhiEta,1));
//
//        }
//        if(dCharge<0){
//          fPOIDiffQRe[h][1][2]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//          fPOIDiffQIm[h][1][2]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//          fPOIDiffMul[h][1][2]->Fill(xval,pow(wPhiEta,1));
//        }
//        if(dEta>=0){
//          if(dCharge>0){
//            fPOIDiffQRe_pos[h][1][1]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//            fPOIDiffQIm_pos[h][1][1]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//            fPOIDiffMul_pos[h][1][1]->Fill(xval,pow(wPhiEta,1));}
//          if(dCharge<0){
//            fPOIDiffQRe_pos[h][1][2]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//            fPOIDiffQIm_pos[h][1][2]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//            fPOIDiffMul_pos[h][1][2]->Fill(xval,pow(wPhiEta,1));}
//        }
//
//        if(dEta<0){
//          if(dCharge>0){
//            //if(dPhi>TMath::Pi)
//            fPOIDiffQRe_neg[h][1][1]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//            fPOIDiffQIm_neg[h][1][1]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//            fPOIDiffMul_neg[h][1][1]->Fill(xval,pow(wPhiEta,1));}
//          if(dCharge<0){
//            fPOIDiffQRe_neg[h][1][2]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//            fPOIDiffQIm_neg[h][1][2]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//            fPOIDiffMul_neg[h][1][2]->Fill(xval,pow(wPhiEta,1));}
//        }
//      }//end of Pid ±221
//
//      // Kaons
//      if (Pid == 321 || Pid == -321) {
//        fPOIDiffQRe[h][2][0]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//        fPOIDiffQIm[h][2][0]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//        fPOIDiffMul[h][2][0]->Fill(xval,pow(wPhiEta,1));
//        if(dCharge>0){
//          //if(dPhi>TMath::Pi)
//          fPOIDiffQRe[h][2][1]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//          fPOIDiffQIm[h][2][1]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//          fPOIDiffMul[h][2][1]->Fill(xval,pow(wPhiEta,1));
//        }
//        if(dCharge<0){
//          fPOIDiffQRe[h][2][2]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//          fPOIDiffQIm[h][2][2]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//          fPOIDiffMul[h][2][2]->Fill(xval,pow(wPhiEta,1));
//        }
//        if(dEta>=0){
//          if(dCharge>0){
//            fPOIDiffQRe_pos[h][2][1]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//            fPOIDiffQIm_pos[h][2][1]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//            fPOIDiffMul_pos[h][2][1]->Fill(xval,pow(wPhiEta,1));}
//          if(dCharge<0){
//            fPOIDiffQRe_pos[h][2][2]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//            fPOIDiffQIm_pos[h][2][2]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//            fPOIDiffMul_pos[h][2][2]->Fill(xval,pow(wPhiEta,1));}
//        }
//
//        if(dEta<0){
//          if(dCharge>0){
//            fPOIDiffQRe_neg[h][2][1]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//            fPOIDiffQIm_neg[h][2][1]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//            fPOIDiffMul_neg[h][2][1]->Fill(xval,pow(wPhiEta,1));}
//          if(dCharge<0){
//            fPOIDiffQRe_neg[h][2][2]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//            fPOIDiffQIm_neg[h][2][2]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//            fPOIDiffMul_neg[h][2][2]->Fill(xval,pow(wPhiEta,1));}
//        }
//      } //end of Pid ±321
//
//      // Protons
//      if (Pid == 2212 || Pid == -2212) {
//        fPOIDiffQRe[h][3][0]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//        fPOIDiffQIm[h][3][0]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//        fPOIDiffMul[h][3][0]->Fill(xval,pow(wPhiEta,1));
//
//        if(dCharge>0){
//          //if(dPhi>TMath::Pi)
//          fPOIDiffQRe[h][3][1]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//          fPOIDiffQIm[h][3][1]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//          fPOIDiffMul[h][3][1]->Fill(xval,pow(wPhiEta,1));
//
//        }
//        if(dCharge<0){
//          fPOIDiffQRe[h][3][2]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//          fPOIDiffQIm[h][3][2]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//          fPOIDiffMul[h][3][2]->Fill(xval,pow(wPhiEta,1));
//        }
//        if(dEta>=0){
//          if(dCharge>0){
//            //if(dPhi>TMath::Pi)
//            fPOIDiffQRe_pos[h][3][1]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//            fPOIDiffQIm_pos[h][3][1]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//            fPOIDiffMul_pos[h][3][1]->Fill(xval,pow(wPhiEta,1));}
//          if(dCharge<0){
//            fPOIDiffQRe_pos[h][3][2]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//            fPOIDiffQIm_pos[h][3][2]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//            fPOIDiffMul_pos[h][3][2]->Fill(xval,pow(wPhiEta,1));}
//        }
//        if(dEta<0){
//          if(dCharge>0){
//            //if(dPhi>TMath::Pi)
//            fPOIDiffQRe_neg[h][3][1]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//            fPOIDiffQIm_neg[h][3][1]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//            fPOIDiffMul_neg[h][3][1]->Fill(xval,pow(wPhiEta,1));}
//          if(dCharge<0){
//            fPOIDiffQRe_neg[h][3][2]->Fill(xval,pow(wPhiEta,1)*TMath::Cos((h+1.)*dPhi));
//            fPOIDiffQIm_neg[h][3][2]->Fill(xval,pow(wPhiEta,1)*TMath::Sin((h+1.)*dPhi));
//            fPOIDiffMul_neg[h][3][2]->Fill(xval,pow(wPhiEta,1));}
//        }
//      }//end of Pid ±2212
//
//    }// end of h loop
//    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
//
  } // end of event  loop

  if (doQA) {
    fMultChargedParticlesDistribution->Fill(nChargedParticles);
    fNumberOfParticipantsDistribution->Fill(fNumberOfParticipants);
    fChargedParticleSpectra->Fill(fCentralityEBE, nChargedParticles/(2*maxEtaCut));
  }
  sw.tock(); //stop timer

  // cout << "Time for looping over event: "<< sw.takeTime() << endl;
  // cout <<"Generic Framework matrix loop: " << totaltime << endl;

  //sw.tick();
//  CalculateFlowQC();
  // CalculateFlowGF();
//  CalculateFlowEP();
  //    CalculateFlowEP1();
//  CalculateFlowRP();

  // sw.tock();
  // cout << "Time for CalculateFlowCQ() "<< sw.takeTime() << endl;
  ResetEventByEventQuantities();
}
//
//
////=====================================================================================================
//
void CalculateFlow::ResetEventByEventQuantities()
{
//  //FlowGF
//  fReQGF->Zero();
//  fImQGF->Zero();
//  for(Int_t i=0; i<fkGFPtB; i++) {
//    fReQGFPt[i]->Zero();
//    fImQGFPt[i]->Zero();
//  }
//
//  //FlowQC
//  for(Int_t c=0;c<fQVecPower;c++) {
//    for (Int_t h=0;h<fFlowNHarmMax;h++) {
//      for (Int_t p=0; p<fNParticles; p++){
//        for (Int_t charge=0; charge<fCharge; charge++){
//          if(fPOIPtDiffQRe[c][h][p][charge]) fPOIPtDiffQRe[c][h][p][charge]->Reset();
//          if(fPOIPtDiffQIm[c][h][p][charge]) fPOIPtDiffQIm[c][h][p][charge]->Reset();
//          if(fPOIPtDiffMul[c][h][p][charge]) fPOIPtDiffMul[c][h][p][charge]->Reset();
//        }
//      }
//    }
//  }
//
//  //FlowEP & RP
//  for (Int_t h=0;h<fFlowNHarmMax;h++) {
//
//    if(QRe_EP[h]!=0) QRe_EP[h]=0;
//    if(QIm_EP[h]!=0) QIm_EP[h]=0;
//    if(Mul_EP[h]!=0) Mul_EP[h]=0;
//
//    for(Int_t p=0;p<fNParticles;p++){
//      for (Int_t charge=0; charge<fCharge; charge++){
//
//        if(fPOIDiffQRe_pos[h][p][charge]) fPOIDiffQRe_pos[h][p][charge]->Reset();
//        if(fPOIDiffQIm_pos[h][p][charge]) fPOIDiffQIm_pos[h][p][charge]->Reset();
//        if(fPOIDiffMul_pos[h][p][charge]) fPOIDiffMul_pos[h][p][charge]->Reset();
//
//        if(fPOIDiffQRe_neg[h][p][charge]) fPOIDiffQRe_neg[h][p][charge]->Reset();
//        if(fPOIDiffQIm_neg[h][p][charge]) fPOIDiffQIm_neg[h][p][charge]->Reset();
//        if(fPOIDiffMul_neg[h][p][charge]) fPOIDiffMul_neg[h][p][charge]->Reset();
//
//        if(fPOIDiffQRe[h][p][charge]) fPOIDiffQRe[h][p][charge]->Reset();
//        if(fPOIDiffQIm[h][p][charge]) fPOIDiffQIm[h][p][charge]->Reset();
//        if(fPOIDiffMul[h][p][charge]) fPOIDiffMul[h][p][charge]->Reset();
//      }
//    }
//  }
  
  
}

//=====================================================================================================

void CalculateFlow::Terminate(Int_t Nevents)
{
  FinalizeSpectra(Nevents);
//   FinalizeFlowQC();
//  //  FinalizeFlowGF();
//  FinalizeFlowEP();
//  FinalizeFlowRP();
//  FinalizeQA();
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


void CalculateFlow::CalculateFlowEP()
{
  
  Double_t FillPtBin = 0;
  Float_t QRe, QIm, Mu;
  Float_t QRe_pos, QIm_pos, Mu_pos, QRe_neg, QIm_neg, Mu_neg;
  Float_t vx, vx_pos, vx_neg;
  Float_t x, x_pos, x_neg;
  
  Float_t EventPlane;
  
  Float_t cosEvPl, sinEvPl;
  Float_t cosPhi, sinPhi;
  Float_t EP_res;
  
  Float_t qpRe, qpIm, qpM, meanPtdiff;
  
  //********************************************************************
  // Calculate Q vectors for the spectators (P&T)
  // Then determine correlations with the POI (u)
  // *******************************************************************
  
  
  for (Int_t hr=0;hr<4;hr++) {
    for(Int_t p=0; p<fNParticles;p++){
      for(Int_t charge=0;charge<fCharge;charge++){
        // ********************************************************************
        // pT-integrated: {2} *************************************************
        // ********************************************************************
        
        QRe=0.; QIm=0.;Mu = 0.;
        QRe_pos=0.; QIm_pos=0.;Mu_pos = 0.;
        QRe_neg=0.; QIm_neg=0.;Mu_neg = 0.;
        
        
        vx=0.;vx_pos=0.;vx_neg=0.;
        
        EventPlane = TMath::ATan2(QIm_EP[hr], QRe_EP[hr])/(hr+1.);
        
        cosEvPl = TMath::Cos(EventPlane);
        sinEvPl = TMath::Sin(EventPlane);

        for(Int_t pt=0; pt<fNBins; pt++) {
          QRe += fPOIDiffQRe[hr][p][charge]->GetBinContent(pt+1);
          QIm += fPOIDiffQIm[hr][p][charge]->GetBinContent(pt+1);
          Mu += fPOIDiffMul[0][p][charge]->GetBinContent(pt+1);
          
          QRe_pos += fPOIDiffQRe_pos[hr][p][charge]->GetBinContent(pt+1);
          QIm_pos += fPOIDiffQIm_pos[hr][p][charge]->GetBinContent(pt+1);
          Mu_pos += fPOIDiffMul_pos[hr][p][charge]->GetBinContent(pt+1);

          //Only for eta<0
          QRe_neg += fPOIDiffQRe_neg[hr][p][charge]->GetBinContent(pt+1);
          QIm_neg += fPOIDiffQIm_neg[hr][p][charge]->GetBinContent(pt+1);
          Mu_neg += fPOIDiffMul_neg[hr][p][charge]->GetBinContent(pt+1);

        }

        x = QRe*cosEvPl + QIm*sinEvPl;
        x_pos = QRe_pos*cosEvPl + QIm_pos*sinEvPl;
        x_neg = QRe_neg*cosEvPl + QIm_neg*sinEvPl;

        if(TMath::Abs(Mu)>0.){
          vx = (x) / Mu; }
        
        if(TMath::Abs(Mu_pos)>0.){
          vx_pos = (x_pos) / Mu_pos; }
        
        if(TMath::Abs(Mu_neg)>0.){
          vx_neg = (x_neg) / Mu_neg; }

        //            std::cout<<Denom_pty<<std::endl;

        fFlowEPIntPro[hr][p][charge]->Fill(fImpactParameter, vx,1.); //1 for weights
        fFlowEPIntPro_pos[hr][p][charge]->Fill(fImpactParameter, vx_pos,1.); //1 for weights
        fFlowEPIntPro_neg[hr][p][charge]->Fill(fImpactParameter, vx_neg,1.); //1 for weights
        
        
        // store pt-differential flow ******************************************************************************************************************************************
        for(Int_t pt=0; pt<fNBins; pt++) {
          
          FillPtBin = fPOIDiffQRe[hr][p][charge]->GetBinCenter(pt+1);
          qpRe=0.; qpIm=0.; qpM=0.;
          qpRe = fPOIDiffQRe[hr][p][charge]->GetBinContent(pt+1);
          qpIm = fPOIDiffQIm[hr][p][charge]->GetBinContent(pt+1);
          qpM = fPOIDiffMul[hr][p][charge]->GetBinContent(pt+1);
          
          
          if(qpM>0) {
            meanPtdiff = qpRe*cosEvPl + qpIm*sinEvPl;
            fFlowEPCorPro[hr][p][charge]->Fill(FillPtBin, meanPtdiff/qpM, 1.);            // ADD: fPOIRPDiffDiffQRe[h]
            
          }
          
        }//end pt
        
        
      } //end of for(Int_t charge=0;charge<fCharge;charge++)
    } //end of p
    fEPEPresolutionPro[hr]->Fill(EventPlane);
  }//end of nHarm
}


void CalculateFlow::CalculateFlowRP()
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
  for (Int_t h=0;h<4;h++) {
    for(Int_t p=0; p<fNParticles;p++){
      //      std::cout<<" p is : "<<p<<std::endl;
      for(Int_t charge=0;charge<fCharge;charge++){
        // ********************************************************************
        // pT-integrated ******************************************************
        // ********************************************************************
        //        std::cout<<" p in charge is : "<<p<<std::endl;
        QRe=0.; QIm=0.;
        Mu = 0.;

        QRe_pos =0; QIm_pos=0; Mu_pos=0;
        QRe_neg =0; QIm_neg=0; Mu_neg=0;

        for(Int_t pt=0; pt<fNBins; pt++) {
          QRe_pos += fPOIDiffQRe_pos[h][p][charge]->GetBinContent(pt+1);
          QIm_pos += fPOIDiffQIm_pos[h][p][charge]->GetBinContent(pt+1);
          Mu_pos += fPOIDiffMul_pos[h][p][charge]->GetBinContent(pt+1);

          QRe_neg += fPOIDiffQRe_neg[h][p][charge]->GetBinContent(pt+1);
          QIm_neg += fPOIDiffQIm_neg[h][p][charge]->GetBinContent(pt+1);
          Mu_neg += fPOIDiffMul_neg[h][p][charge]->GetBinContent(pt+1);


        }


        if(TMath::Abs(Mu_neg)>0.){
          v_neg = QRe_neg/Mu_neg;
          fFlowRPIntPro_neg[h][p][charge]->Fill(fImpactParameter, v_neg , 1.);
          fFlowRPIntPro[h][p][charge]->Fill(fImpactParameter, (QRe_neg+QRe_pos)/(Mu_pos+Mu_neg),1.);
        }

        if(TMath::Abs(Mu_pos)>0.) {
          v_pos = QRe_pos/Mu_pos;
          fFlowRPIntPro_pos[h][p][charge]->Fill(fImpactParameter, v_pos , 1.); //1 for weights
        }
        
        
        // store pt-differential flow ******************************************************************************************************************************************
        
        for(Int_t pt=0; pt<fNBins; pt++) {
          
          FillPtBin = fPOIDiffQRe[h][p][charge]->GetBinCenter(pt+1);
          qpRe=0.; qpIm=0.; qpM=0.;
          qpRe = fPOIDiffQRe[h][p][charge]->GetBinContent(pt+1);
          qpIm = fPOIDiffQIm[h][p][charge]->GetBinContent(pt+1);
          qpM = fPOIDiffMul[h][p][charge]->GetBinContent(pt+1);
          
//          std::cout<<" p in ptdiff is : "<<p<< " and qpM is: "<< qpM<<std::endl;
          if(qpM>0) {
            meanPtdiff = qpRe/qpM;
            fFlowRPCorPro[h][p][charge]->Fill(FillPtBin, meanPtdiff, 1.);            // ADD: fPOIRPDiffDiffQRe[h]
            //            std::cout<< "particle: " << p << " h="<<h<< " charge="<<charge<<" v="<<meanPtdiff<<std::endl;
          }
          
        }//end pt
        
      } //end charge loop
    }//end p loop
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
  
  for (Int_t p =0 ; p<fNParticles;p++){
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
          QRe += fPOIPtDiffQRe[1][hr][p][charge]->GetBinContent(pt+1); // Cos((hr+1.)*dPhi)
          QIm += fPOIPtDiffQIm[1][hr][p][charge]->GetBinContent(pt+1); // Sin((hr+1.)*dPhi)
          Q2Re2 += fPOIPtDiffQRe[2][2*hr+1][p][charge]->GetBinContent(pt+1); // w^2*Cos((hr+1.)*dPhi)
          Q2Im2 += fPOIPtDiffQIm[2][2*hr+1][p][charge]->GetBinContent(pt+1); // w^2*Sin((hr+1.)*dPhi)
          QRe3 += fPOIPtDiffQRe[3][hr][p][charge]->GetBinContent(pt+1); // w^3*Cos((hr+1.)*dPhi)
          QIm3 += fPOIPtDiffQIm[3][hr][p][charge]->GetBinContent(pt+1); // w^3*Sin((hr+1.)*dPhi)
          
          QM0 += fPOIPtDiffMul[0][0][p][charge]->GetBinContent(pt+1); // w^0
          QM  += fPOIPtDiffMul[1][0][p][charge]->GetBinContent(pt+1); // w^1
          QM2 += fPOIPtDiffMul[2][0][p][charge]->GetBinContent(pt+1); // w^2
          QM3 += fPOIPtDiffMul[3][0][p][charge]->GetBinContent(pt+1); // w^3
          QM4 += fPOIPtDiffMul[4][0][p][charge]->GetBinContent(pt+1); // w^4
        }
        
        IQM2 = QM*QM-QM2;
        WQM2 = (WeigMul? IQM2 : 1.);
        if(QM0>1) {
          IQC2[hr] = (QRe*QRe+QIm*QIm-QM2)/IQM2; // <2> = |Qn,1|^2-S1,2/(S2,1-S1,2)
          fFlowQCIntCorPro[hr][0][p][charge]->Fill(fImpactParameter,IQC2[hr],WQM2*fCenWeightEbE); // nPart vs. <2>
          fFlowQCRefCorPro[hr][0][p][charge]->Fill(fImpactParameter,IQC2[hr],WQM2*fCenWeightEbE); // nPart vs. <2>
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
          
          fFlowQCIntCorPro[hr][1][p][charge]->Fill(fImpactParameter,IQC4[hr],WQM4*fCenWeightEbE);
          fFlowQCRefCorPro[hr][1][p][charge]->Fill(fImpactParameter,IQC4[hr],WQM4*fCenWeightEbE);
          Q4f = kTRUE;
        }
        
        // product of correlations or covariances
        if(Q2f && Q4f) {
          fFlowQCIntCorPro[hr][2][p][charge]->Fill(fImpactParameter,IQC2[hr]*IQC4[hr],WQM2*WQM4*fCenWeightEbE);
          fFlowQCRefCorPro[hr][13][p][charge]->Fill(fImpactParameter,IQC2[hr]*IQC4[hr],WQM2*WQM4*fCenWeightEbE);
        }
        
        
        // ********************************************************************
        // pT-differential: {2}, {4} ******************************************
        // ********************************************************************
        
        // store pt-differential flow ****************************************
        for(Int_t pt=0; pt<fNBins; pt++) {
          
          FillPtBin = fPOIPtDiffQRe[1][hr][p][charge]->GetBinCenter(pt+1);
          
          qpRe0=0.; qpIm0=0.; qpRe2=0.; qpIm2=0.; qp2Re=0.; qp2Im=0.; qpM0=0.; qpM=0.; qpM2=0.; qpM3=0.;
          
          //here change h+2 to h+1 (hr+1->hr, 2hr+3-> 2hr+1)
          //std::cout<<fPOIPtDiffQRe[0][hr][p][charge]->GetNbinsX()<<std::endl;
          
          qpRe0 = fPOIPtDiffQRe[0][hr][p][charge]->GetBinContent(pt+1);
          qpIm0 = fPOIPtDiffQIm[0][hr][p][charge]->GetBinContent(pt+1);
          qpRe2 = fPOIPtDiffQRe[2][hr][p][charge]->GetBinContent(pt+1);
          qpIm2 = fPOIPtDiffQIm[2][hr][p][charge]->GetBinContent(pt+1);
          qp2Re = fPOIPtDiffQRe[1][2*hr+1][p][charge]->GetBinContent(pt+1);
          qp2Im = fPOIPtDiffQIm[1][2*hr+1][p][charge]->GetBinContent(pt+1);
          
          qpM0 = fPOIPtDiffMul[0][0][p][charge]->GetBinContent(pt+1);
          qpM  = fPOIPtDiffMul[1][0][p][charge]->GetBinContent(pt+1);
          qpM2 = fPOIPtDiffMul[2][0][p][charge]->GetBinContent(pt+1);
          qpM3 = fPOIPtDiffMul[3][0][p][charge]->GetBinContent(pt+1);
          
          //if(hr==0) {//////////////////////////////////////////////////////////////////////////////////
          // fFlowFlowectra->Fill(fImpactParameter,FillPtBin,qpM*fCenWeightEbE);
          //}
          
          dQM2 = qpM0*QM-qpM;
          WdQM2 = (WeigMul? dQM2 : 1.);
          
          if(qpM0>0 && QM0>0) {
            dQC2 = (qpRe0*QRe+qpIm0*QIm-qpM)/dQM2;
//            std::cout<<"pt bin: " << FillPtBin<< " filled with : " << dQC2 << " In cent bin: "<<fCenBin<<std::endl;
            fFlowQCCorPro[fCenBin][hr][1][p][charge]->Fill(FillPtBin,dQC2,WdQM2*fCenWeightEbE);
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
            fFlowQCCorPro[fCenBin][hr][2][p][charge]->Fill(FillPtBin,dQC4,WdQM4*fCenWeightEbE);
            dQ4f = kTRUE;
          }
          
          // product of correlations or covariances
          if(Q2f && dQ2f) fFlowQCCorCovPro[fCenBin][hr][0][p][charge]->Fill(FillPtBin,IQC2[hr]*dQC2,WQM2*WdQM2*fCenWeightEbE);
          if(Q4f && dQ2f) fFlowQCCorCovPro[fCenBin][hr][1][p][charge]->Fill(FillPtBin,IQC4[hr]*dQC2,WQM4*WdQM2*fCenWeightEbE);
          if(Q2f && dQ4f) fFlowQCCorCovPro[fCenBin][hr][2][p][charge]->Fill(FillPtBin,IQC2[hr]*dQC4,WQM2*WdQM4*fCenWeightEbE);
          if(dQ2f && dQ4f) fFlowQCCorCovPro[fCenBin][hr][3][p][charge]->Fill(FillPtBin,dQC2*dQC4,WdQM2*WdQM4*fCenWeightEbE);
          if(Q4f && dQ4f) fFlowQCCorCovPro[fCenBin][hr][4][p][charge]->Fill(FillPtBin,IQC4[hr]*dQC4,WQM4*WdQM4*fCenWeightEbE);
          
        } // end of for(Int_t pt=0; pt<fCRCnPtBin; pt++)
        
      } // end of for(Int_t hr=0; hr<fFlowNHarm; hr++)
      
      
      
      // pt diff get fFlowQCCorHist ==============================
      
      
      // Pt-DIFFERENTIAL
      for(Int_t hr=0; hr<fFlowNHarm; hr++) {
        for(Int_t j=0; j<fFlowQCNRef; j++) {
          for(Int_t pt=1;pt<=fNBins;pt++) { // fNBins was: fFlowQCRefCorPro[hr][j][p][charge]->GetNbinsX()
            Double_t stats[6]={0.};
            fFlowQCRefCorPro[hr][j][p][charge]->GetXaxis()->SetRange(pt,pt);
            fFlowQCRefCorPro[hr][j][p][charge]->GetStats(stats);
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
                fFlowQCRefCorHist[hr][j][p][charge]->SetBinContent(pt,Corr);
                fFlowQCRefCorHist[hr][j][p][charge]->SetBinError(pt,CorrErr);
              }
            }
          } // end of for(Int_t pt=1;pt<=100;pt++)
          fFlowQCRefCorPro[hr][j][p][charge]->GetXaxis()->SetRange(1,fFlowQCRefCorPro[hr][j][p][charge]->GetNbinsX());
        } // end of for(Int_t j=0; j<5; j++)
        
        
        for (Int_t h=0; h<fCRCnCen; h++) {
          
          // STORE IN HISTOGRAMS
          
          for(Int_t j=0; j<fFlowQCNPro; j++) {
            for(Int_t pt=1;pt<=fNBins;pt++) { //was 13
              
              Double_t stats[6]={0.};
              fFlowQCCorPro[h][hr][j][p][charge]->GetXaxis()->SetRange(pt,pt);
              fFlowQCCorPro[h][hr][j][p][charge]->GetStats(stats);
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
                  fFlowQCCorHist[h][hr][j][p][charge]->SetBinContent(pt,Corr);
//                  std::cout<<"for h="<<h<<" hr="<<hr<<" j="<<j<<" | Corr="<<Corr<< " for pt="<<pt<<"\n";
                  fFlowQCCorHist[h][hr][j][p][charge]->SetBinError(pt,CorrErr);
                }
              }
              
            } // end of for(Int_t pt=1;pt<=fNBins;pt++)
            fFlowQCCorPro[h][hr][j][p][charge]->GetXaxis()->SetRange(1,416);
          }
          
          // reference flow
          // Look for flow in bin h -> Depends on centrality
          // 2- and 4-particle cumulants
          
          Double_t QC2    = fFlowQCRefCorHist[hr][0][p][charge]->GetBinContent(h+1);
//          std::cout<<"Cn2: "<<QC2<<" for h = "<<h<<std::endl;
          Double_t QC2E = fFlowQCRefCorHist[hr][0][p][charge]->GetBinError(h+1);
          Double_t QC4    = fFlowQCRefCorHist[hr][1][p][charge]->GetBinContent(h+1);
          Double_t QC4E = fFlowQCRefCorHist[hr][1][p][charge]->GetBinError(h+1);
          Double_t Cn2 = QC2;
          Double_t Cn2E = QC2E;
          Double_t wCov24 = fFlowQCRefCorHist[hr][13][p][charge]->GetBinContent(h+1);
          Double_t Cn4 = QC4-2.*QC2*QC2;
          Double_t Cn4Esq = 16.*pow(QC2,2.)*pow(QC2E,2) + pow(QC4E,2.) - 8.*QC2*wCov24;
          
          Double_t Cos1 = fFlowQCRefCorHist[hr][3][p][charge]->GetBinContent(h+1); // <<cos(n*phi1)>>
          Double_t Sin1 = fFlowQCRefCorHist[hr][4][p][charge]->GetBinContent(h+1); // <<sin(n*phi1)>>
          Double_t Sin1P2 = fFlowQCRefCorHist[hr][5][p][charge]->GetBinContent(h+1);
          Double_t Cos1P2 = fFlowQCRefCorHist[hr][6][p][charge]->GetBinContent(h+1);
          Double_t Sin1M2M3 = fFlowQCRefCorHist[hr][7][p][charge]->GetBinContent(h+1);
          Double_t Cos1M2M3 = fFlowQCRefCorHist[hr][8][p][charge]->GetBinContent(h+1);
          // change vocabulary, to be changed
          Double_t cosP1nPhi = fFlowQCRefCorHist[hr][3][p][charge]->GetBinContent(h+1); // <<cos(n*phi1)>>
          Double_t sinP1nPhi = fFlowQCRefCorHist[hr][4][p][charge]->GetBinContent(h+1); // <<sin(n*phi1)>>
          Double_t sinP1nPhi1P1nPhi2 = fFlowQCRefCorHist[hr][5][p][charge]->GetBinContent(h+1); //sin(n*(phi1+phi2))
          Double_t cosP1nPhi1P1nPhi2 = fFlowQCRefCorHist[hr][6][p][charge]->GetBinContent(h+1);  //cos(n*(phi1+phi2))
          Double_t sinP1nPhi1M1nPhi2M1nPhi3 = fFlowQCRefCorHist[hr][7][p][charge]->GetBinContent(h+1);  //sin(n*(phi1-phi2-phi3))
          Double_t cosP1nPhi1M1nPhi2M1nPhi3 = fFlowQCRefCorHist[hr][8][p][charge]->GetBinContent(h+1); //cos(n*(phi1-phi2-phi3))
          
          fFlowQCRefCorFinal[hr][0][p][charge]->SetBinContent(h+1,Cn2);
          fFlowQCRefCorFinal[hr][0][p][charge]->SetBinError(h+1,Cn2E);
          
          if(Cn4Esq>0.) {
            Double_t Cn4E = pow(Cn4Esq,0.5);
            fFlowQCRefCorFinal[hr][3][p][charge]->SetBinContent(h+1,Cn4);
            fFlowQCRefCorFinal[hr][3][p][charge]->SetBinError(h+1,Cn4E);
            if(Cn4<0.) {
              Double_t Flow4 = pow(fabs(Cn4),0.25);
              Double_t Flow4E = fabs(Flow4/(4.*Cn4))*Cn4E;
              fFlowQCRefCorFinal[hr][1][p][charge]->SetBinContent(h+1,Flow4);
              fFlowQCRefCorFinal[hr][1][p][charge]->SetBinError(h+1,Flow4E);
            }
          }
          
          // pt-differential
          for(Int_t pt=1; pt<=fNBins; pt++) {
            Double_t qp2    = fFlowQCCorHist[h][hr][1][p][charge]->GetBinContent(pt);
//            if(qp2){std::cout<<"for h="<<h<<" hr="<<hr<<" j="<<1<<" | pt: "<<pt<<" qp2="<<qp2<<" Cn2: "<<Cn2<<std::endl;}
            Double_t qp2E = fFlowQCCorHist[h][hr][1][p][charge]->GetBinError(pt);
            Double_t qp4    = fFlowQCCorHist[h][hr][2][p][charge]->GetBinContent(pt);
            Double_t qp4E = fFlowQCCorHist[h][hr][2][p][charge]->GetBinError(pt);
            Double_t Dn2 = qp2;
            Double_t Dn2E = qp2E;
            Double_t Dn4 = qp4-2.*qp2*QC2;
            Double_t wCovTwoFourReduced = fFlowQCCorCovHist[h][hr][1][p][charge]->GetBinContent(pt);
            Double_t wCovTwoReducedFourReduced = fFlowQCCorCovHist[h][hr][4][p][charge]->GetBinContent(pt);
            Double_t Dn4Esq = 4.*pow(QC2,2.)*pow(qp2E,2) + 4.*pow(qp2,2.)*pow(QC2E,2) + pow(qp4E,2.) - 4.*qp2*wCovTwoFourReduced - 4.*QC2*wCovTwoReducedFourReduced;
            
            
            fFlowQCFinalPtDifHist[h][hr][5][p][charge]->SetBinContent(pt,Dn2);
            fFlowQCFinalPtDifHist[h][hr][5][p][charge]->SetBinError(pt,Dn2E);
            
            if(Cn2) {
              Double_t Flow2 = Dn2/sqrt(fabs(Cn2));
//              std::cout<<"for h="<<h<<" hr="<<hr<<" j="<<1<<" | pt: "<<pt<<" Flow2="<<Flow2<<std::endl;
              Double_t Flow2E = 0.;
              // change vocabulary, to be changed
              Double_t two = QC2; //std::cout<<two<<std::endl;
              Double_t twoError = QC2E; //std::cout<<twoError<<std::endl;
              Double_t twoReduced = qp2; //std::cout<<twoReduced<<std::endl;
              Double_t twoReducedError = qp2E; //std::cout<<twoReducedError<<std::endl;
              Double_t wCovTwoTwoReduced = fFlowQCCorCovHist[h][hr][0][p][charge]->GetBinContent(pt);
              Double_t v2PrimeErrorSquared = (1./4.)*pow(two,-3.)*(pow(twoReduced,2.)*pow(twoError,2.)
                                                                   + 4.*pow(two,2.)*pow(twoReducedError,2.)
                                                                   - 4.*two*twoReduced*wCovTwoTwoReduced);
              if(v2PrimeErrorSquared>0.){Flow2E = pow(v2PrimeErrorSquared,0.5);}
              
              
              if(Flow2E>0.) {
//                std::cout<< "Flow2: h="<<h<<" hr = " <<hr<< " p="<<p<<" c="<<charge<<" voor v="<<Flow2<<std::endl;
                fFlowQCFinalPtDifHist[h][hr][0][p][charge]->SetBinContent(pt,Flow2);
                fFlowQCFinalPtDifHist[h][hr][0][p][charge]->SetBinError(pt,Flow2E);
              }
            }
            
            if(Dn4Esq>0.) {
              Double_t Dn4E = pow(Dn4Esq,0.5);
              fFlowQCFinalPtDifHist[h][hr][6][p][charge]->SetBinContent(pt,Dn4);
              fFlowQCFinalPtDifHist[h][hr][6][p][charge]->SetBinError(pt,Dn4E);
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
              Double_t wCovTwoTwoReduced = fFlowQCCorCovHist[h][hr][0][p][charge]->GetBinContent(pt);
              Double_t wCovTwoFourReduced = fFlowQCCorCovHist[h][hr][1][p][charge]->GetBinContent(pt);
              Double_t wCovFourTwoReduced = fFlowQCCorCovHist[h][hr][2][p][charge]->GetBinContent(pt);
              Double_t wCovFourFourReduced = fFlowQCCorCovHist[h][hr][3][p][charge]->GetBinContent(pt);
              Double_t wCovTwoReducedFourReduced = fFlowQCCorCovHist[h][hr][4][p][charge]->GetBinContent(pt);
              Double_t wCovTwoFour = fFlowQCRefCorHist[hr][13][p][charge]->GetBinContent(h+1);
              
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
                fFlowQCFinalPtDifHist[h][hr][1][p][charge]->SetBinContent(pt,Flow4);
                fFlowQCFinalPtDifHist[h][hr][1][p][charge]->SetBinError(pt,Flow4E);
              }
            }
          } // end of for(Int_t pt=1; pt<=fNBins; pt++) {
        } // end of for (Int_t h=0; h<fCRCnCen; h++) {
      } // end of for(Int_t hr=0; hr<fFlowNHarm; hr++)
    } //end of for(Int_t charge=0;charge<fCharge;charge++)
  } // end of for (Int_t p =0 ; p<fNParticles;p++)
}

//=======================================================================================================================

//=======================================================================================================================
void CalculateFlow::FinalizeFlowQC()
{
  std::cout << "Finalizing Flow QC"<<'\n'<<endl;
  std::cout << "Print bins where v_n{2} or v_n{4} cannot be determined"<<endl;
  cout << "*************************************" << "\n";
  
  for (Int_t p=0; p<fNParticles;p++){
    for(Int_t charge=0; charge<fCharge; charge++){
      for(Int_t hr=0; hr<fFlowNHarm; hr++) {
        // Pt-INTEGRATED
        // STORE IN HISTOGRAMS
        // 2- and 4-particle cumulants
        for(Int_t j=0; j<fkFlowQCnIntCorPro; j++) {
          for(Int_t pt=1;pt<=fFlowQCIntCorPro[hr][j][p][charge]->GetNbinsX();pt++) { //pt is hier helamaal geen pt
            Double_t stats[6]={0.};
            fFlowQCIntCorPro[hr][j][p][charge]->GetXaxis()->SetRange(pt,pt);
            fFlowQCIntCorPro[hr][j][p][charge]->GetStats(stats);
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
                fFlowQCIntCorHist[hr][j][p][charge]->SetBinContent(pt,Corr); //= <<2>> = cn{2} = vn{2}^2
                fFlowQCIntCorHist[hr][j][p][charge]->SetBinError(pt,CorrErr);
                
                if(Corr>0.){
                  Double_t Flow2 = pow(fabs(Corr),0.5);
                  Double_t Flow2E = fabs(Flow2/(4.*Corr))*CorrErr;
                  
                  if(j==0){
                    fFlowQCIntFlow2Hist[hr][0][p][charge]->SetBinContent(pt,Flow2); //neem sqrt(Corr) -> v_n{2}
                    fFlowQCIntFlow2Hist[hr][0][p][charge]->SetBinError(pt,Flow2E);}
                }
                
              }
            }
          } // end of for(Int_t pt=1;pt<=100;pt++)
          fFlowQCIntCorPro[hr][j][p][charge]->GetXaxis()->SetRange(1,fFlowQCIntCorPro[hr][j][p][charge]->GetNbinsX());
        } // end of for(Int_t j=0; j<5; j++)
        
      } // end of for(Int_t hr=0; hr<fFlowNHarm; hr++)
    }// end of for(Int_t charge=0; charge<fCharge; charge++)
  }//end of p
  
  
  
  // FINALISE (calculate flow)
  for (Int_t p=0; p<fNParticles;p++){
    for(Int_t charge=0; charge<fCharge; charge++){
      for(Int_t hr=0; hr<fFlowNHarm; hr++) {
        // calculate covariance
        for(Int_t pt=1; pt<=fFlowQCIntCorHist[hr][0][p][charge]->GetNbinsX(); pt++) {
          // average reduced correlations:
          Double_t two = fFlowQCIntCorHist[hr][0][p][charge]->GetBinContent(pt); // <<2>>
          Double_t four = fFlowQCIntCorHist[hr][1][p][charge]->GetBinContent(pt); // <<4>>
          // sum of weights for reduced correlation:
          Double_t sumOfWeightsForTwo = GetSumPro(fFlowQCIntCorPro[hr][0][p][charge],pt); // sum_{i=1}^{N} w_{<2>}
          Double_t sumOfWeightsForFour = GetSumPro(fFlowQCIntCorPro[hr][1][p][charge],pt); // sum_{i=1}^{N} w_{<4>}
          // product of weights for reduced correlation:
          Double_t productOfWeightsForTwoFour = GetSumPro(fFlowQCIntCorPro[hr][2][p][charge],pt); // sum_{i=1}^{N} w_{<2>}w_{<4>}
          // products for differential flow:
          Double_t twoFour = fFlowQCIntCorHist[hr][2][p][charge]->GetBinContent(pt); // <<2><4>>
          
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
              fFlowQCIntCorHist[hr][2][p][charge]->SetBinContent(pt,wCovTwoFour);
            }
          }
        } // end of for(Int_t pt=1;pt<=fNBins;pt++)
        
        // 2- and 4-particle cumulants
        for(Int_t pt=1; pt<=fFlowQCIntCorHist[hr][0][p][charge]->GetNbinsX(); pt++) {
          Double_t QC2    = fFlowQCIntCorHist[hr][0][p][charge]->GetBinContent(pt); //<<2>>
          Double_t QC2E   = fFlowQCIntCorHist[hr][0][p][charge]->GetBinError(pt);
          Double_t QC4    = fFlowQCIntCorHist[hr][1][p][charge]->GetBinContent(pt); //<<4>>
          Double_t QC4E   = fFlowQCIntCorHist[hr][1][p][charge]->GetBinError(pt);
          Double_t wCov24 = fFlowQCIntCorHist[hr][2][p][charge]->GetBinContent(pt); //<<2>><<4>>
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
          
          
          fFlowQCIntCumHist[hr][0][p][charge]->SetBinContent(pt,Cn2);
          fFlowQCIntCumHist[hr][0][p][charge]->SetBinError(pt,Cn2E);
          
          if(Cn4Esq>0.) {
            Double_t Cn4E = pow(Cn4Esq,0.5);
            fFlowQCIntCumHist[hr][1][p][charge]->SetBinContent(pt,Cn4);
            fFlowQCIntCumHist[hr][1][p][charge]->SetBinError(pt,Cn4E);
            
            if (Cn4<0.) {
              Double_t Flow4 = pow(fabs(Cn4),0.25);
              Double_t Flow4E = fabs(Flow4/(4.*Cn4))*Cn4E;
              
              fFlowQCIntCorHist[hr][2][p][charge]->SetBinContent(pt,Flow4);
              fFlowQCIntCorHist[hr][2][p][charge]->SetBinError(pt,Flow4E);
              
              fFlowQCIntFlow4Hist[hr][0][p][charge]->SetBinContent(pt,Flow4); //neem sqrt(Corr) -> v_n{2}
              fFlowQCIntFlow4Hist[hr][0][p][charge]->SetBinError(pt,Flow4E);
            }
            
            //                    } else {
            //                        fFlowQCIntCorHist[hr][2][p][charge]->SetBinContent(pt,0.);
            //                        fFlowQCIntCorHist[hr][2][p][charge]->SetBinError(pt,0.);
            //
            //                       // fFlowQCIntFlow4Hist[hr][0][p][charge]->SetBinContent(pt,0.); //neem sqrt(Corr) -> v_n{2}
            //                       // fFlowQCIntFlow4Hist[hr][0][p][charge]->SetBinError(pt,0.); //
            //
            //                    }
          }
        }
      } // end of for(Int_t hr=0; hr<fFlowNHarm; hr++)
    } // end of for(Int_t charge=0; charge<fCharge; charge++)
  } // end of p
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

void CalculateFlow::FinalizeFlowEP()
{
  std::cout << "Finalizing Flow EP"<<'\n'<<endl;
  std::cout << "Average correlations over all events"<<endl;
  cout << "*************************************" << "\n";
  
  Float_t v_p =0.;
  Float_t v_t=0.;
  
  Float_t EP_res;
  
  
  for(Int_t charge=0; charge<fCharge; charge++){
    for (Int_t p=0;p<fNParticles; p++){
      for (Int_t h=0;h<4;h++) {
        for(Int_t pt=1;pt<=fFlowEPIntPro[h][p][charge]->GetNbinsX();pt++) {
          Float_t Corr_QQ_y = 0; Double_t CorrErr_QQ_y = 0;
          Corr_QQ_y = GetWeightedCorrelations(fFlowEPIntPro[h][p][charge], pt);
          CorrErr_QQ_y = GetWeightedCorrelationsError(fFlowEPIntPro[h][p][charge], pt);
          fFlowEPIntFlow2Hist[h][p][charge]->SetBinContent(pt, Corr_QQ_y);
          fFlowEPIntFlow2Hist[h][p][charge]->SetBinError(pt, CorrErr_QQ_y);
        }
        
        for(Int_t pt=1;pt<=fFlowEPIntPro_pos[h][p][charge]->GetNbinsX();pt++) {
          Float_t Corr_QQ_y = 0; Double_t CorrErr_QQ_y = 0;
          Corr_QQ_y = GetWeightedCorrelations(fFlowEPIntPro_pos[h][p][charge], pt);
          CorrErr_QQ_y = GetWeightedCorrelationsError(fFlowEPIntPro_pos[h][p][charge], pt);
          fFlowEPIntFlow2Hist_pos[h][p][charge]->SetBinContent(pt, Corr_QQ_y);
          fFlowEPIntFlow2Hist_pos[h][p][charge]->SetBinError(pt, CorrErr_QQ_y);
        }
        
        for(Int_t pt=1;pt<=fFlowEPIntPro_neg[h][p][charge]->GetNbinsX();pt++) {
          Float_t Corr_QQ_y = 0; Double_t CorrErr_QQ_y = 0;
          Corr_QQ_y = GetWeightedCorrelations(fFlowEPIntPro_neg[h][p][charge], pt);
          CorrErr_QQ_y = GetWeightedCorrelationsError(fFlowEPIntPro_neg[h][p][charge], pt);
          fFlowEPIntFlow2Hist_neg[h][p][charge]->SetBinContent(pt, Corr_QQ_y);
          fFlowEPIntFlow2Hist_neg[h][p][charge]->SetBinError(pt, CorrErr_QQ_y);
        }
        
        for(Int_t pt=1;pt<=fFlowEPCorPro[h][p][charge]->GetNbinsX();pt++) {
          Float_t Corr = 0; Double_t CorrErr = 0;
          Corr= GetWeightedCorrelations(fFlowEPCorPro[h][p][charge], pt);
          CorrErr = GetWeightedCorrelationsError(fFlowEPCorPro[h][p][charge], pt);
          if(Corr && CorrErr){
            fFlowEPDiffFlow2Hist[h][p][charge]->SetBinContent(pt, Corr);
            fFlowEPDiffFlow2Hist[h][p][charge]->SetBinError(pt, CorrErr);}
          else {
            fFlowEPDiffFlow2Hist[h][p][charge]->SetBinContent(pt, 0);
          }
        }
        
      } // end of harm
    }// end of p
  }// end of charge
}

//=====================================================================================================



void CalculateFlow::FinalizeFlowRP()
{
  std::cout << "Finalizing Flow RP"<<'\n'<<endl;
  std::cout << "Average correlations over all events"<<endl;
  cout << "*************************************" << "\n";
  
  Float_t v_p =0.;
  Float_t v_t=0.;
  
  for (Int_t h=0;h<4;h++) {
    for (Int_t p=0;p<fNParticles; p++){
      for(Int_t charge=0; charge<fCharge; charge++){
        
        
        for(Int_t pt=1;pt<=fFlowRPIntPro_pos[h][p][charge]->GetNbinsX();pt++) {
          Float_t Corr = 0; Double_t CorrErr = 0;
          Corr = GetWeightedCorrelations(fFlowRPIntPro_pos[h][p][charge], pt);
          CorrErr = GetWeightedCorrelationsError(fFlowRPIntPro_pos[h][p][charge], pt);
          fFlowRPIntFlow2Hist_pos[h][p][charge]->SetBinContent(pt, Corr); //Let op!! Hier moet min als eta<0 in make function
          fFlowRPIntFlow2Hist_pos[h][p][charge]->SetBinError(pt, CorrErr);
        }

        for(Int_t pt=1;pt<=fFlowRPIntPro_neg[h][p][charge]->GetNbinsX();pt++) {
          Float_t Corr = 0; Double_t CorrErr = 0;
          Corr = GetWeightedCorrelations(fFlowRPIntPro_neg[h][p][charge], pt);
          CorrErr = GetWeightedCorrelationsError(fFlowRPIntPro_neg[h][p][charge], pt);
          fFlowRPIntFlow2Hist_neg[h][p][charge]->SetBinContent(pt, Corr); //Let op!! Hier moet min als eta<0 in make function
          fFlowRPIntFlow2Hist_neg[h][p][charge]->SetBinError(pt, CorrErr);
        }

        for(Int_t pt=1;pt<=fFlowRPIntPro[h][p][charge]->GetNbinsX();pt++) {
          Float_t Corr = 0; Double_t CorrErr = 0;
          Corr = GetWeightedCorrelations(fFlowRPIntPro[h][p][charge], pt);
          CorrErr = GetWeightedCorrelationsError(fFlowRPIntPro[h][p][charge], pt);
          fFlowRPIntFlow2Hist[h][p][charge]->SetBinContent(pt, Corr); //Let op!! Hier moet min als eta<0 in make function
          fFlowRPIntFlow2Hist[h][p][charge]->SetBinError(pt, CorrErr);
        }

        
        for(Int_t pt=1;pt<=fFlowRPCorPro[h][p][charge]->GetNbinsX();pt++) {
          Float_t Corr = 0; Double_t CorrErr = 0;
          Corr= GetWeightedCorrelations(fFlowRPCorPro[h][p][charge], pt);
          CorrErr = GetWeightedCorrelationsError(fFlowRPCorPro[h][p][charge], pt);
          if(Corr && CorrErr){
            fFlowRPDiffFlow2Hist[h][p][charge]->SetBinContent(pt, Corr); //Let op!! Hier moet min als eta<0 in make function
            fFlowRPDiffFlow2Hist[h][p][charge]->SetBinError(pt, CorrErr);}
          else {
            fFlowRPDiffFlow2Hist[h][p][charge]->SetBinContent(pt, 0); //Let op!! Hier moet min als eta<0 in make function
          }
        }
        
      }  // end of for(Int_t charge=0; charge<fCharge; charge++)
    } // end of for (Int_t p=0;p<fNParticles; p++)
  }// end of for (Int_t h=0;h<fFlowNHarm;h++)
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
  /// Let op: We hebben maar 9 bins for impact parameter en we starten bij b van 5% centrality. So we skip the most central events.
  Int_t CenBin=-1;
//  if (Centrality>0. && Centrality<5.) CenBin=0;
  if (Centrality>=5. && Centrality<10.) CenBin=1-1;
  if (Centrality>=10. && Centrality<20.) CenBin=2-1;
  if (Centrality>=20. && Centrality<30.) CenBin=3-1;
  if (Centrality>=30. && Centrality<60.) CenBin=4-1;
  if (Centrality>=40. && Centrality<50.) CenBin=5-1;
  if (Centrality>=50. && Centrality<60.) CenBin=6-1;
  if (Centrality>=60. && Centrality<70.) CenBin=7-1;
  if (Centrality>=70. && Centrality<80.) CenBin=8-1;
  if (Centrality>=80. && Centrality<90.) CenBin=9-1;
  if (Centrality>=90. && Centrality<100.) CenBin=10-1;
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
