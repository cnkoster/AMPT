#include "TFile.h"
#include "TList.h"
#include "TH1.h"

Int_t GetCentBin(Float_t ImpactParameter)
{
    Int_t CenBin=-1;
//    std::cout<<"b is "<<ImpactParameter<<std::endl;
    
    if (ImpactParameter>0. && ImpactParameter<3.72) CenBin=0;       //0-5%
    if (ImpactParameter>=3.72 && ImpactParameter<5.23) CenBin=1;    //5-10%
    if (ImpactParameter>5.23 && ImpactParameter<7.31) CenBin=2;     //10-20%
    if (ImpactParameter>7.31 && ImpactParameter<8.88) CenBin=3;      //20-30%
    if (ImpactParameter>8.88 && ImpactParameter<10.20) CenBin=4;     //30-40%
    if (ImpactParameter>10.20 && ImpactParameter<11.38) CenBin=5;     //40-50%
    if (ImpactParameter>11.38 && ImpactParameter<12.47) CenBin=6;     //50-60%
    if (ImpactParameter>12.47 && ImpactParameter<13.50) CenBin=7;     //60-70%
    if (ImpactParameter>13.50 && ImpactParameter<14.51) CenBin=8;     //70-80%
    if (ImpactParameter>14.51) CenBin=9;                            //80=100%
//    std::cout<<"Centrality bin: " << CenBin <<std::endl;
    
    
    return CenBin;
}

void SetHistogramSettings(TH1* Hist, Int_t color, Int_t MStyle, Int_t FStyle){
    
    Hist->SetLineColor(color);
    Hist->SetMarkerColor(color);
    Hist->SetMarkerStyle(MStyle);
    Hist->SetMarkerSize(1);
    Hist->SetFillStyle(FStyle);
    Hist->SetFillColor(color);
}

TH1D* Switch_b_to_Centrality(TH1* hist, Int_t n){
    Double_t centbin[] = {5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 100.0};
    Int_t nCentBins = 9;
    TH1D* centHist = new TH1D(Form("v%i_vs_Centrality",n), Form("v%i_vs_Centrality",n), nCentBins, centbin);
    
    Double_t b=0.;
    Int_t Centrality_Bin=0;

    for(Int_t i=0;i<=hist->GetNbinsX();i++){
        b = hist->GetBinCenter(i+1);
        Centrality_Bin = GetCentBin(b);
        centHist->SetBinContent(Centrality_Bin,hist->GetBinContent(i+1));
        centHist->SetBinError(Centrality_Bin, hist->GetBinError(i+1));
    }
        
        return centHist;
}


void comparator(Int_t iMin=0, Int_t iMax=1000) {
    std::vector<int> colors = {kPink+7, 38, kGreen+3, kOrange+7, kBlack, kViolet+1, kCyan+1};
    
  TFile *f = TFile::Open("/data/alice/nkoster/AMPT_out/Run2_Energy/nEvent100/EPM/AnalysisResults_Group0-2000_fullEta_EPMSPM.root");
                           //"/data/alice/nkoster/AMPT_out/Run2_Energy/nEvents100/AnalysisResultsQC_Group0-1000_v1-4_ImpactParam.root");
  if (!f || f->IsZombie()) {
            cout << "Error: Failed to open ROOT file." << endl;
            f->Close();
            return;
        }
    // AnalysisResults_Group0-1000_testEta.root");
    //Form("/data/alice/nkoster/AMPT_out/Run2_Energy/nEvents100/AnalysisResults_Group%i-%i_%s.root",iMin,iMax,addon));
    TList *list = (TList*)f->Get("FLowSPMList");
    TH1D* Hist1 = (TH1D*)list->FindObject("fFlowEPMIntFlow2Hist[0][0]");
  if (!Hist1) {
    cout << "Error: Failed to open histo" << endl;
    //           f->Close();
    //           return;
  }
    
    TH1D* Hist2 = (TH1D*)list->FindObject("fFlowSPMIntFlow2Hist[0][0]");
    TH1D* Hist3 = (TH1D*)list->FindObject("fFlowSPM1IntFlow2Hist[0][0]");
   // TH1D* Hist4 = (TH1D*)list->FindObject("fFlowQCIntFlow2Hist[3][0][0]");
    
    
    Double_t centbin[] = {5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 100.0};
    Int_t nCentBins = 8;
    TH1D* v2Hist = Switch_b_to_Centrality(Hist1,1); //new TH1D("v2_vs_Centrality", "v2_vs_Centrality", nCentBins, centbin);
    TH1D* v3Hist = Switch_b_to_Centrality(Hist2,2); //new TH1D("v3_vs_Centrality", "v3_vs_Centrality", nCentBins, centbin);
    TH1D* v4Hist = Switch_b_to_Centrality(Hist3,3); //new TH1D("v4_vs_Centrality", "v4_vs_Centrality", nCentBins, centbin);
    
    
//    //Rootfile data article https://arxiv.org/pdf/1306.4145.pdf
//    TFile *data = new TFile(Form("HEPData-ins1666817-v1-Table_1.root")); //193 for v2{2}
//    if (!data || data->IsZombie()) {
//           cout << "Error: Failed to open ROOT file." << endl;
//           f->Close();
//           return;
//       }
//
//    TFile *datav3 = new TFile(Form("HEPData-ins1666817-v1-Table_3.root"));
//    if (!datav3 || datav3->IsZombie()) {
//           cout << "Error: Failed to open ROOT file." << endl;
//           f->Close();
//           return;
//       }
//
//    TFile *datav4 = new TFile(Form("HEPData-ins1666817-v1-Table_4.root"));
//    if (!datav4|| datav4->IsZombie()) {
//           cout << "Error: Failed to open ROOT file." << endl;
//           f->Close();
//           return;
//       }
//
//    TH1F* dataHist_v2 = (TH1F*)data->Get("Table 1/Hist1D_y1"); //193 for v2{2}
//    TH1F* statHist_v2 = (TH1F*)data->Get("Table 1/Hist1D_y1_e1"); //193 for v2{2} (e2=sys)
//    TH1F* sysHist_v2 = (TH1F*)data->Get("Table 1/Hist1D_y1_e2");
//
//
//    for(Int_t i=1;i<dataHist_v2->GetNbinsX()+1;i++){
//        statHist_v2->SetBinContent(i,dataHist_v2->GetBinContent(i));
//        sysHist_v2->SetBinContent(i,dataHist_v2->GetBinContent(i));
//    }
//
//    if (!dataHist_v2) {
//          cout << "Error: Histogram v2 not found.. " << endl;
//          f->Close();
//          data->Close();
//          return;
//      }
//
//
//    TH1F* dataHist_v3 = (TH1F*)datav3->Get("Table 3/Hist1D_y1");
//    TH1F* statHist_v3 = (TH1F*)datav3->Get("Table 3/Hist1D_y1_e1");
//    TH1F* sysHist_v3 = (TH1F*)datav3->Get("Table 3/Hist1D_y1_e2");
//
//    for(Int_t i=1;i<dataHist_v3->GetNbinsX()+1;i++){
//        statHist_v3->SetBinContent(i,dataHist_v3->GetBinContent(i));
//        sysHist_v3->SetBinContent(i,dataHist_v3->GetBinContent(i));
//    }
//    if (!dataHist_v3) {
//          cout << "Error: Histogram v3 not found.. " << endl;
//          f->Close();
//          data->Close();
//          return;
//      }
//
//
//    TH1F* dataHist_v4 = (TH1F*)datav4->Get("Table 4/Hist1D_y1");
//    TH1F* statHist_v4 = (TH1F*)datav4->Get("Table 4/Hist1D_y1_e1");
//    TH1F* sysHist_v4 = (TH1F*)datav4->Get("Table 4/Hist1D_y1_e2");
//
//    for(Int_t i=1;i<dataHist_v4->GetNbinsX()+1;i++){
//        statHist_v4->SetBinContent(i,dataHist_v4->GetBinContent(i));
//        sysHist_v4->SetBinContent(i,dataHist_v4->GetBinContent(i));
//    }
//    if (!dataHist_v4) {
//          cout << "Error: Histogram v4 not found.. " << endl;
//          f->Close();
//          data->Close();
//          return;
//      }
    
    //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    //-=-=-=-=-=-=-= Draw v2, v3, v4 -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
    auto c = c11("v2_v3_v4");

    auto p = (TPad*)gROOT->FindObject(TString::Format("pv2_v3_v4"));
    p->cd();
   // p->SetLogy();

    auto leg = new TLegend(0.65, 0.6, 0.9, 0.8, "#splitline{Pb-Pb #sqrt{s_{NN}} = 5.02 TeV }{|#eta| < 0.8, 0.2 < p_{T} < 3 GeV}");
    leg->SetFillStyle(0);
    leg->SetTextSize(18);
    leg->SetNColumns(2);
    leg->SetColumnSeparation(0.001);
    p->cd();



    SetHistogramSettings(v2Hist, colors[0], 8, 0);
//    SetHistogramSettings(sysHist_v2, colors[1], 33, 0);
//    SetHistogramSettings(statHist_v2, colors[1], 33, 0);
    
    SetHistogramSettings(v3Hist, colors[2], 8, 0);
//    SetHistogramSettings(sysHist_v3, colors[3], 33, 0);
//    SetHistogramSettings(statHist_v3, colors[3], 33, 0);
    
    SetHistogramSettings(v4Hist, colors[4], 8, 0);
//    SetHistogramSettings(sysHist_v4, colors[5], 33, 0);
//    SetHistogramSettings(statHist_v4, colors[5], 33, 0);


    v2Hist->GetYaxis()->SetTitleOffset(1.2);
    v2Hist->GetXaxis()->SetNdivisions(506);
    v2Hist->GetYaxis()->SetNdivisions(504);

    v2Hist->SetTitle("v_{1}{EPM}");
    v3Hist->SetTitle("v_{1}{SPM}");
    v4Hist->SetTitle("v_{1}{SPM}");
    v2Hist->GetXaxis()->SetTitle("Centrality (%)");
    v2Hist->GetYaxis()->SetRangeUser(-0.004,0.005);
    v2Hist->GetYaxis()->SetTitle("v_{n}");

    
//    leg->AddEntry(v2Hist, "AMPT", "");
//    leg->AddEntry(v2Hist, "Data", "");
    
    leg->AddEntry(v2Hist, v2Hist->GetTitle(), "p");
   // leg->AddEntry(statHist_v2, "v_{2}{2,|#Delta#eta|>1}", "p");
    
    
    leg->AddEntry(v3Hist, v3Hist->GetTitle(), "p");
    //leg->AddEntry(statHist_v3, "v_{3}{2,|#Delta#eta|>1}", "p");
    
    leg->AddEntry(v4Hist, v4Hist->GetTitle(), "p");
   // leg->AddEntry(statHist_v4, "v_{4}{2,|#Delta#eta|>1}", "p");

    v2Hist->Draw("");
    //sysHist_v2->Draw("E2 SAME");
    //statHist_v2->Draw("E1 SAME");
    
    v3Hist->Draw("SAME");
   // sysHist_v3->Draw("E2 SAME");
   // statHist_v3->Draw("E1 SAME");
    
    v4Hist->Draw("SAME");
    //sysHist_v4->Draw("E2 SAME");
   // statHist_v4->Draw("E1 SAME");
    
    leg->Draw("same");
    
    //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    //-=-=-=-=-=-=-= Draw v1 -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
//    TH1D* v1Hist = Switch_b_to_Centrality(Hist_v1,1);
//
//    auto c1 = c11("v1");
//
//    auto p1 = (TPad*)gROOT->FindObject(TString::Format("pv1"));
//    p1->cd();
//    p1->SetLogy();
//
//    auto leg1= new TLegend(0.65, 0.6, 0.9, 0.8, "#splitline{Pb-Pb #sqrt{s_{NN}} = 5.02 TeV }{|#eta| < 0.8, 0.2 < p_{T} GeV}");
//    leg1->SetFillStyle(0);
//    leg1->SetTextSize(18);
//   // leg->SetNColumns(2);
//   // leg->SetColumnSeparation(0.001);
//    p1->cd();
//
//
//
//    SetHistogramSettings(v1Hist, colors[0], 8, 0);
////    SetHistogramSettings(sysHist_v2, colors[1], 33, 0);
////    SetHistogramSettings(statHist_v2, colors[1], 33, 0);
//
//    v1Hist->GetYaxis()->SetTitleOffset(1.2);
//    v1Hist->GetXaxis()->SetNdivisions(506);
//    v1Hist->GetYaxis()->SetNdivisions(504);
//
//    v1Hist->SetTitle("v_{1}{2}");
//
//    v1Hist->GetXaxis()->SetTitle("Centrality (%)");
//    v1Hist->GetYaxis()->SetRangeUser(5e-3,1);
//    v1Hist->GetYaxis()->SetTitle("v_{n}");
//
//
//    leg1->AddEntry(v1Hist, "AMPT", "");
//    //leg->AddEntry(v2Hist, "Data", "");
//
//    leg1->AddEntry(v1Hist, v2Hist->GetTitle(), "p");
//    //leg->AddEntry(statHist_v2, "v_{2}{2,|#Delta#eta|>1}", "p");
//
////
////    leg->AddEntry(v3Hist, v3Hist->GetTitle(), "p");
////    leg->AddEntry(statHist_v3, "v_{3}{2,|#Delta#eta|>1}", "p");
////
////    leg->AddEntry(v4Hist, v4Hist->GetTitle(), "p");
////    leg->AddEntry(statHist_v4, "v_{4}{2,|#Delta#eta|>1}", "p");
//
//    v1Hist->Draw("");
////    sysHist_v2->Draw("E2 SAME");
////    statHist_v2->Draw("E1 SAME");
////
////    v3Hist->Draw("SAME");
////    sysHist_v3->Draw("E2 SAME");
////    statHist_v3->Draw("E1 SAME");
////
////    v4Hist->Draw("SAME");
////    sysHist_v4->Draw("E2 SAME");
////    statHist_v4->Draw("E1 SAME");
////
////    leg->Draw("same");
////
//    //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//    //-=-=-=-=-=-=-= Write to outputfile -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//    //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
    TFile* output = TFile::Open(TString::Format("outComp_AMPT_Data.root"), "recreate");
    c->Write("v2_v3_v4");
//    c1->Write("v1");
    
    output->Close();
    delete output;
    
    
}



