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


void runFlow(Bool_t etaFlag = kTRUE, TString centrality="", Double_t gCentrality=35, Int_t iGroupMin=0, Int_t iGroupMax=3000, Double_t gSpectatorMin=0, Double_t gSpectatorMax=500) {
  
  TFile *f= new TFile(Form("/dcache/alice/nkoster/PhD/AMPT_out/Run2_Energy_PbPb/Cent%s/TreeOutput_Group%i-%i_Cent%s.root",centrality.Data(),iGroupMin, iGroupMax, centrality.Data()));
                           //"/data/alice/nkoster/TreeOutput_Group0-6000_Cent30_60.root"));
  //"/dcache/alice/nkoster/PhD/AMPT_out/Run2_Energy_PbPb/nEvents100/TreeOutput/TreeOutput_Group%i-%i.root",iGroupMin, iGroupMax));
  
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
    

    fQC->SetEvent(event);
   
    fQC->UserExec();
    
    
  }
  fQC->Terminate(nevent);
  
  // Save list holding histogram with weights:
  TFile *fResultsFile = new TFile(Form("AnalysisResults_Group%i-%i_%sDiff_Full_Cent%s.root",iGroupMin, iGroupMax, diff.Data(), centrality.Data()),"RECREATE");
                                       //"AnalysisResults_Group0-6000_%sDiff_PID_3060_QC.root", diff.Data()),"RECREATE");
  //"AnalysisResults_Group%i-%i_%sDiff_All.root",iGroupMin, iGroupMax, diff.Data()),"RECREATE");
  fResultsFile->WriteObject(fQC->GetQAList(),"QAList","SingleKey");
  fResultsFile->WriteObject(fQC->GetSpectraList(),"SpectraList","SingleKey");
  fResultsFile->WriteObject(fQC->GetFlowQCList(),"FLowQCList","SingleKey");
  //fResultsFile->WriteObject(fQC->GetFlowGFList(),"FLowGFList","SingleKey");
  fResultsFile->WriteObject(fQC->GetFlowEPList(),"FlowEPList","SingleKey");
  fResultsFile->WriteObject(fQC->GetFlowRPList(),"FlowRPList","SingleKey");
  fResultsFile->Close();
  
  
}
