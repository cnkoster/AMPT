#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <stdio.h>
#include "Particle.h"
#include "Event.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
using namespace std;

Float_t ParticleCharge(TDatabasePDG *pdg, int pid){
    if (pid == 42 || pid == -42){
        Float_t charge = (pid < 0) ? -1.0 : 1.0;
        return charge;
    } else {
        TParticlePDG *part = pdg->GetParticle(pid);
        Float_t charge = part->Charge()/3.0; //waarom gedeeld door 3?
        return charge;
    }
}
//Run macro in folder that contains the Group{i} folders!! 


void amptTTree_v6spec(Int_t Njobs = 2001, Int_t Nsplit = 1, Int_t nEventsperJob = 100){
    // Global variables
    TDatabasePDG *pdg = new TDatabasePDG();
    //Int_t Njobs = 5;// number of groups to read the data from:q
    Int_t j=0;// int for the loop over the particles in an event
    Int_t k=0; // int for the loop over zpc particles
    Int_t eventEnd = 0;// int to specify the line in the file at which the event ends
    Int_t Nevent = 0;
    // Int_t Nsplit = 1;
    Int_t Ngroups = Njobs/Nsplit;
    Int_t beginJob, endJob;
    Float_t epsilon=0.0001f;
    Int_t Anumber = 136;
    string Centrality[1] = {" "}; //{"30_40", "40_50", "50_60"};
    
    Bool_t fileFilled = false;
    TH1F *hEventCounter = new TH1F("hEventCounter","Number of events",10,0,1);
    TH1F *hPtDist = new TH1F("hPtDist", "Pt distribution (no spectators)",100,0,10);
    TH1F *hPtDist_specIncluded = new TH1F("hPtDist_specIncluded", "Pt distribution (including spectators)",100,0,10);
    
    // Loop over subjobs defined by Nsplit
    // for {a : b} loop over elements in b where a is element
    
    for (string cent : Centrality){
        for (Int_t batch = 0; batch<Nsplit; batch++) {
            Nevent=0;
            beginJob = batch*Ngroups+1; //= 0,n+1, 2n+1
            endJob = (batch+1)*Ngroups; //=n n, 2n, 3n
         //   cout<< "we starten de for loop over batches"<<endl;
            // Create a .root file to fill with the tree
           // string treePath = Form("/dcache/alice/nkoster/PhD/AMPT_out/Run2_Energy_PbPb/Cent%s/TreeOutput_Group%i-%i_Cent%s.root", cent.c_str(),  beginJob-1,endJob-1, cent.c_str());
            string treePath = Form("/dcache/alice/nkoster/PhD/AMPT_out/Run2_Energy_PbPb/nEvents100/TreeOutput/TreeOutput_Group%i-%i.root", beginJob-1,endJob-1);
            cout << " path to rootfile: " << treePath << endl;
            
            TFile *f = new TFile(treePath.c_str(),"recreate");
            
            if (f->IsZombie()){
                cout<<"Can't open the file"<< endl;
                return;
            }
            
            // Create the tree and it's branches
            TTree *eventTree = new TTree("EventTree","AMPT simulation");
            
            // create the event: see Event.h for class reference
            Event event;
            
            // add event branch to TTree
            eventTree->Branch("event",&event);
            
            
            
            // Loop over the different groups
            for (Int_t i=beginJob; i<=endJob; i++){
                fileFilled = false;
                j = 0;
                eventEnd = 0;
                
               // if (i==54) {
                //    cout << " iGroup "<<i-1<<"is empty -> SKIP " <<endl;
                //    continue;
               // }
                // Open the .dat file
                
                string Path = Form("/dcache/alice/nkoster/PhD/AMPT_out/Run2_Energy_PbPb/nEvents100/Group"); //, cent.c_str());
                
                string FilePath = Path + to_string(i-1) + "/ampt.dat"; //add $i to get folder Group{i}/ampt.dat file
                
                cout << " Opening: " << FilePath.c_str() << endl;
                
                FILE *fp = fopen(FilePath.c_str(),"r"); //open ampt.dat
                
                if (fp == NULL){
                    cout<<"Can't load in AMPT data"<<endl;
                    cout << " iGroup "<<i-1<<"is empty -> SKIP " <<endl;
                    continue;
                }
                
                cout << " Opened file: " << FilePath.c_str() << endl;
                
                vector<Particle> particle; //create particle see Particle.h
                
                
                Int_t pid{0}, nTrack{0}, n_Part{0}, nPartP{0}, nPartT{0}, nElP{0}, nInelP{0}, nElT{0}, nInelT{0}, seqN{0}, status{-1}, pidPart{-1},  pidOrig{-1};
                Int_t specFlag{-999};
                Float_t px{-999.}, py{-999.}, pz{-999.}, mass{-999.}, charge{-999.}, x{-999.}, y{-999.}, z{-999.}, t{-999.}, b{-999.}, xPart{-999.}, yPart{-999.}, zPart{-999.};
        
                
                // Loop over the lines ampt.dat
                char line[81];
                while (fgets(line,81,fp)){
                    fileFilled = true;
                    // If its the start of an event, get the global event parameters
                    if (j==eventEnd){
                        if (j!=0){
                            if(nTrack!=particle.size()){
                            cout<<"PAS OP! nTracks is niet gelijk aan aantal particles in vector!! nEvent: "<<Nevent<<endl;
                                break;
                            }
                            Nevent++;
                            //
                            
                            event = Event(nTrack,b,nPartP,nPartT,nElP,nInelP,nElT,nInelT,particle); //Define event here when particle vector is full
                            // kan ook hier toevoegen welke centralilty bin je wil!!!!! 
                            eventTree->Fill();
                            particle.clear();
                        }
                        // read in event line of ampt.dat file
                        sscanf(&line[0],"%*i %*i %i %f %i %i %i %i %i %i %*i",&(nTrack),&(b),&(nPartP),&(nPartT),&(nElP),&(nInelP),&(nElT),&(nInelT));
                        // determine end of event by number of tracks in event
                        eventEnd = j + nTrack+1;
                        //cout<< "nTrack is: " << nTrack << "\n";
                        hEventCounter->Fill(0.5); // fill histo to keep track of number of events
                       // particle.clear(); //moet hier (na event = Event......) anders zit je nog in je vorige event en verwijder je daar al je particles...
                        j++;
                        continue;
                    }
                    
                    // Get the particle variables
                    sscanf(&line[0],"%i %f %f %f %f %f %f %f %f",&pid,&px,&py,&pz,&mass,&x,&y,&z,&t);
                    specFlag = ((abs(px) < epsilon) && (abs(py) < epsilon) && (pid==2212 || pid==2112)) ? -1 : 1; // set -1 for spectator particles (included in event data) -1 means particle = spectator
                  //  if(specFlag<0){std::cout<<"We've got one! SPECTATOR pz = "<< pz << " | event is: "<<Nevent<<std::endl;}
                    charge = ParticleCharge(pdg, pid);
                    auto part = Particle(pid,specFlag,px,py,pz,mass,charge,x,y,z,t); //add all particles including spectator. Use flag later to eliminate spectators.
                    particle.push_back(part);
                    hPtDist_specIncluded->Fill(part.getPt());
                    if(!specFlag){hPtDist->Fill(part.getPt());} //fill pt distribution WITHOUT spectators
                    j++;
                    
                }
                
             //   cout<<"get event properties: " << event.getParticle(4).getPt()<<endl;
               // cout<<"number of particles in vector: "<<particle.size()<<"\n";
                if(nTrack!=particle.size()){
                cout<<"PAS OP! nTracks is niet gelijk aan aantal particles in vector!! nEvent: "<<Nevent<<endl;
                    break;
                }
                event = Event(nTrack,b,nPartP,nPartT,nElP,nInelP,nElT,nInelT,particle);
                eventTree->Fill();
                particle.clear();
                Nevent++;
                cout << "Nevent is " << Nevent << " group number is " << i-1 << endl;
                
                if (!fileFilled) continue;
               
                fclose(fp);
                
                //
            }
            
            cout << "Now we are filling the TTree" << endl;
            eventTree->Print();
            eventTree->Write();
            hEventCounter->Write();
            hPtDist->Write();
            hPtDist_specIncluded->Write();
            f->Close();
            delete f;
        }
    }
}
