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
//test
//====================Spectators===================//
//60-70: 25% @ 376, 50% @ 382, 75% @ 388, 100% @ 499
//50-60: 25% @ 348, 50% @ 356, 75% @ 365, 100% @ 499
//40-50: 25% @ 311, 50% @ 323, 75% @ 332, 100% @ 499
//30-40: 25% @ 266, 50% @ 280, 75% @ 291, 100% @ 499
//20-30: 25% @ 204, 50% @ 223, 75% @ 239, 100% @ 499
//10-20: 25% @ 124, 50% @ 144, 75% @ 167, 100% @ 499
//=================================================//

//========================q2=======================//
//60-70: 25% @ 1.17, 50% @ 1.79, 75% @ 2.46, 100% @ 20
//50-60: 25% @ 1.55, 50% @ 2.31, 75% @ 3.13, 100% @ 20
//40-50: 25% @ 2.02, 50% @ 2.95, 75% @ 3.88, 100% @ 20
//30-40: 25% @ 2.46, 50% @ 3.42, 75% @ 4.39, 100% @ 20
//20-30: 25% @ 2.47, 50% @ 3.56, 75% @ 4.64, 100% @ 20
//10-20: 25% @ 2.19, 50% @ 3.20, 75% @ 4.34, 100% @ 20
//========================q2=======================//

void runFlow(Bool_t etaFlag = kTRUE, TString centrality="", Double_t gCentrality=1., Int_t iGroupMin=0, Int_t iGroupMax=2000, Double_t gSpectatorMin=0, Double_t gSpectatorMax=500) {
  
  std::cout<< "here??" <<std::endl;
  
  TFile *f= new TFile(Form("/data/alice/nkoster/TreeOutput_Group0-6000_Cent30_50.root"));
  //"/dcache/alice/nkoster/PhD/AMPT_out/Run2_Energy_PbPb/nEvents100/TreeOutput/TreeOutput_Group%i-%i.root",iGroupMin, iGroupMax));
  //"/dcache/alice/nkoster/PhD/AMPT_out/Run2_Energy_PbPb/Cent%s/TreeOutput_Group%i-%i_Cent%s.root",centrality.Data(),iGroupMin, iGroupMax+1000, centrality.Data()));
  //"/dcache/alice/nkoster/PhD/AMPT_out/Run2_Energy_PbPb/nEvents100/TreeOutput/TreeOutput_Group%i-%i.root",iGroupMin, iGroupMax));
  //"/dcache/alice/nkoster/PhD/AMPT_out/Run2_Energy_PbPb/Cent%s/TreeOutput_Group%i-%i_Cent%s.root",centrality.Data(),iGroupMin, iGroupMax, centrality.Data()));
  
  if (!f) {
    std::cout<<"Input file does not exist"<<"\n";
    return;
  } else {
    std::cout<<"Input file exists"<<"\n";
  }
  TTree *tree = (TTree*)f->Get("EventTree");
  
  // tree->AddFile(Form("5.02TeV/Centrality%s/tree_PaperPreProduction_Group%d.root",centrality.Data(),i));
  
  // create a pointer to an event object. This will be used
  // to read the branch values.
  Event *event = new Event();
  
  // get the branch and set the branch address
  TBranch *bnevent = tree->GetBranch("event");
  bnevent->SetAddress(&event);
  
  Long64_t nevent = 100;//tree->GetEntries();
  Int_t nselected = 0;
  Int_t nb = 0;
  Double_t nSpectators = 0;
  Double_t nPart = 0;
  //  Double_t nPartTotal=416;
  
  
  CalculateFlow *fQC = new CalculateFlow("CalculateFlow");
  
  TString diff;
  if(etaFlag) diff="eta";
  else diff="pt";
  
  fQC->SetEtaDiff(etaFlag); // pt differential flow -> kTRUE for eta diff flow.
  fQC->SetCentralityEBE(gCentrality);
  fQC->UserCreateOutputObjects();
  fQC->SetmaxPtCut(10);
  fQC->SetminPtCut(0.2);
  fQC->SetminNtrackCut(500);
  fQC->SetmaxEtaCut(0.8);
  fQC->SetdoQA(kTRUE);
  
  
  
  for (Long64_t i=0;i<nevent;i++) {
    if (i%100 == 0) std::cout<<"Event: "<<i<<"\n";
    
    bnevent->GetEntry(i); //this is the branch event
    
    std::cout<<"1"<<std::endl;
    fQC->SetEvent(event);
    std::cout<<"2"<<std::endl; 
    fQC->UserExec();
    
    
  }
  fQC->Terminate(nevent);
  
  // Save list holding histogram with weights:
  TFile *fResultsFile = new TFile(Form("AnalysisResults_Group0-6000_etaDiff_PID_3050.root"),"RECREATE");
  //"AnalysisResults_Group%i-%i_%sDiff_All.root",iGroupMin, iGroupMax, diff.Data()),"RECREATE");
  //"AnalysisResults_Group%i-%i_%sDiff_Full_Cent%s.root",iGroupMin, iGroupMax, diff.Data(), centrality.Data()),"RECREATE");
  fResultsFile->WriteObject(fQC->GetQAList(),"QAList","SingleKey");
  fResultsFile->WriteObject(fQC->GetSpectraList(),"SpectraList","SingleKey");
  //fResultsFile->WriteObject(fQC->GetFlowQCList(),"FLowQCList","SingleKey");
  //fResultsFile->WriteObject(fQC->GetFlowGFList(),"FLowGFList","SingleKey");
  fResultsFile->WriteObject(fQC->GetFlowEPList(),"FlowEPList","SingleKey");
  fResultsFile->WriteObject(fQC->GetFlowRPList(),"FlowRPList","SingleKey");
  fResultsFile->Close();
  
  
}
