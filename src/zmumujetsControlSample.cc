#define zmumujetsControlSample_cxx
#include "EmanTreeAnalysis.h"
//C or C++ header files
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip> //for input/output manipulators
//ROOT header files
#include <TAxis.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TVirtualFitter.h>
//my headers
#include "functionsForAnalysis.h"
#include "myClasses.h"
#include "thresholds.h"
//#include "edimarcoTreeFriend.h"

using namespace std;
using namespace myAnalyzerTEman;

#ifdef zmumujetsControlSample_cxx

zmumujetsControlSample::zmumujetsControlSample(TTree *tree, const char* inputSuffix) : edimarcoTree(tree) {
  //cout <<"check in constructor "<<endl;
  suffix = inputSuffix;
  Init(tree);

}

#endif

#define LUMI 5                         // integrated luminosity in fb^-1 to which yields are normalized

void zmumujetsControlSample::loop(vector< Double_t > &yRow, vector< Double_t > &eRow)
{

   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);  
   // warning: in Emanuele's trees non integer values are float

   fChain->SetBranchStatus("weight",1);   // includes k-factor

   fChain->SetBranchStatus("nMu10V",1);  // # of muons passing loose selection
   fChain->SetBranchStatus("nEle10V",1);  // # of electrons passing loose selection for electron veto
   fChain->SetBranchStatus("nGamma15V",1);  // # of photons passing loose selection for photon veto
   fChain->SetBranchStatus("nMu20T",1);  // # of muons passing tight selection (pt > 20 + everything else)

   fChain->SetBranchStatus("dphijj",1);          // dphi between 1st and 2nd jet, 999 if second jet doesn't exist
   fChain->SetBranchStatus("jetclean",1);      // 1 if jet is cleaned, 0 otherwise
   fChain->SetBranchStatus("nJet",1);         // # of jets with pt > 25 && |eta| < 2.5
   fChain->SetBranchStatus("nJet30",1);         // # of jets with pt > 30 && |eta| < 2.4
   fChain->SetBranchStatus("nJet30a",1);       // # of jets with pt > 30 && |eta| < 4.7 
   fChain->SetBranchStatus("Jet_pt",1);  
   fChain->SetBranchStatus("Jet_eta",1);  
 
   fChain->SetBranchStatus("LepGood_pdgId",1);  // must be 13 for muons ( -13 for mu+)
   fChain->SetBranchStatus("LepGood_pt",1);
   fChain->SetBranchStatus("LepGood_eta",1);
   //fChain->SetBranchStatus("LepGood_charge",1);
   fChain->SetBranchStatus("LepGood_tightId",1);
   fChain->SetBranchStatus("LepGood_relIso04",1);
   fChain->SetBranchStatus("genLep_pdgId",1);
   //fChain->SetBranchStatus("m2l",1);  // m(ll)  (I can compute it myself, maybe it's better)
   fChain->SetBranchStatus("mZ1",1);  // best m(ll) SF/OS

   fChain->SetBranchStatus("GenPart_pdgId",1);
   fChain->SetBranchStatus("GenPart_motherId",1);
   //fChain->SetBranchStatus("GenPart_motherIndex",1);

   fChain->SetBranchStatus("metNoMu_pt",1);   

   Double_t metBinEdges[] = {200., 250., 300., 350., 400., 500., 650., 1000.};
   Int_t nMetBins = (sizeof(metBinEdges)/sizeof(Double_t)) - 1;

   vector<Double_t> metCut;
   metCut.push_back(250);
   metCut.push_back(300);
   metCut.push_back(350);
   metCut.push_back(400);
   metCut.push_back(500);
   selection mumetC[metCut.size()];
   for (Int_t i = 0; i < metCut.size(); i++) {
     mumetC[i].set(Form("mumetC[%i]",i),Form("mumet > %3.0lf",metCut.at(i)));
   }
   selection njetsC("njetsC",Form("njets <= %i",NJETS),"pt > 30; |eta| < 4.7");   // using nJet30a
   selection njetsEmanC("njetsEmanC","njets","1 or 2 jets, cleaning, pt > 30, |eta| < 2.4");       // using nJet30
   selection jet1ptC("jet1ptC",Form("jet1pt > %4.0lf",(Double_t)J1PT));
   selection jet1etaC("jet1etaC",Form("|jet1eta| < %2.1lf",J1ETA));
   selection jetCleaningC("jetCleaningC","jet cleaning");               // jetclean is 1 if cleaning is passed, 0 otherwise. It's applied to first jet and , if any, to the second
   selection jet2etaC("jet2etaC",Form("|jet2eta| < %2.1lf",J2ETA),"only if njets = 2");
   selection jet1jet2dphiC("jet1jet2dphiC",Form("|jet1jet2dphi| < %1.1lf",J1J2DPHI),"only if njets = 2");
   selection jjdphiEmanC("jjdphiEmanC","jjdphi < 2.5","only if var \"nJet\" = 2 (not nJet30)");
   selection eLooseVetoC("eLooseVetoC","electrons veto");
   //selection ntausC("ntausC","ntaus","=",NTAUS);
   selection gammaLooseVetoC("gammaLooseVetoC","photons veto");
   //selection mumet200C("mumet200C","mumet",">",200);  it's built-in in our sample
   // additional selections for control sample
   selection oppChargeMuonsC("oppChargeMuonsC","OS muons");
   selection oppChargeLeptonsC("oppChargeLeptonsC","OS/SF leptons");
   selection twomuonsC("twomuonsC","muons");
   selection muLooseVetoC("muLooseVetoC","muons veto");
   selection twomuLooseC("twomuLooseC","2 loose muons");
   // selection nmuTightIdC("nmuTightIdC","nTightMuons",">",0,"tight ID");
   // selection nmuTightIdForEfficiencyC("nmuTightIdForEfficiencyC","nTightMuons","=",2,"tight ID");
   selection mu1tightC("mu1tightC","leading muon tight");   // with Emanuele's tree, it means both tight id and pt > 20 and |eta| < 2.4 (see friend tree)
   selection mu1tightIdC("mu1tightIdC","leading muon tight","tight ID + relIso04 (as Emanuele)");
   selection twomuTightC("twomuTightC","2 tight muons");
   selection mu1ptC("mu1ptC","mu1pt > 20","leading pt muon");
   selection mu2ptC("mu2ptC","mu2Pt > 10");
   selection mu1etaC("mu1etaC","|mu1eta| < 2.4","leading pt muon");  
   selection mu2etaC("mu2etaC","|mu2eta| < 2.4");
   selection mumuInvMassC("mumuInvMassC","muons mass in [60,120]");
   selection genTausC("genTausC","taus generated");               // 11, 13, 15 for e, mu, tau  ( - sign for antiparticles)              
   selection genMuonsC("genMuonsC","muons generated");            // 11, 13, 15 for e, mu, tau  ( - sign for antiparticles)  
   //selection recoMuIdC("recoMuIdC","|recoMuons|","=",1); 

   selection::checkMaskLength();
   //selection::printActiveSelections(cout);

   mask zmumujetsControlSample("Z->mumu control sample selection flow as Emanuele's");
   zmumujetsControlSample.append(twomuLooseC.get2ToId() + oppChargeLeptonsC.get2ToId());
   zmumujetsControlSample.append(twomuonsC.get2ToId());
   zmumujetsControlSample.append(mu1tightIdC.get2ToId());
   zmumujetsControlSample.append(mumuInvMassC.get2ToId());
   zmumujetsControlSample.append(njetsEmanC.get2ToId());
   zmumujetsControlSample.append(jet1ptC.get2ToId());
   zmumujetsControlSample.append(jjdphiEmanC.get2ToId());
   zmumujetsControlSample.append(eLooseVetoC.get2ToId());
   zmumujetsControlSample.append(gammaLooseVetoC.get2ToId());

   Double_t nTotalWeightedEvents = 0.0;    // total events (including weights)

   TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 
   TH1::StatOverflows();                 //enable use of underflows and overflows for statistics computation 
   //TVirtualFitter::SetDefaultFitter("Minuit");

   Int_t Hcolor[] = {1,2,3,4,5,6,7,8,9,12,18,30,38,41,42,46,47,49};    

   TH1D *HmumuInvMass[nMetBins];
   for (Int_t i = 0; i < nMetBins; i++) {
     HmumuInvMass[i] = new TH1D(Form("HmumuInvMass[%i]",i),"",30,60.,120.);
   } 

   TH1D *HzmumujetsYieldsMetBin = new TH1D("HzmumujetsYieldsMetBin","yields of Zmumu control sample in bins of met;MET;# of events",
                                            nMetBins,metBinEdges);

   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"zmumujetsControlSample::loop()"<<endl;
   cout<<"nentries = "<<nentries<<endl;   

   Long64_t nbytes = 0, nb = 0;
   for (Int_t jentry=0; jentry<nentries; jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;   

     UInt_t eventMask = 0; 
     Double_t newwgt = weight * LUMI;

     nTotalWeightedEvents += newwgt;  // counting events with weights

     eventMask += njetsC.addToMask(nJet30a <= NJETS);
     eventMask += njetsEmanC.addToMask( (nJet30 == 1 || nJet30 == 2) && jetclean > 0.5);
     eventMask += jjdphiEmanC.addToMask( nJet30 == 1 || (nJet == 2 && abs(dphijj) < 2.5));
     eventMask += jet1ptC.addToMask(Jet_pt[0] > J1PT);               
     //in Emanuele's tree we have vectors: [0] is the first jet, [1] is the second and so on (ordered in pt)
     eventMask += jet1etaC.addToMask(fabs(Jet_eta[0]) < J1ETA);
     eventMask += jetCleaningC.addToMask(jetclean > 0.5);     
     // jetclean is 1 if cleaning is passed, 0 otherwise. It's applied to first jet and , if any, to the second
     if (nJet30a == 2) {
       eventMask += jet2etaC.addToMask(fabs(Jet_eta[1]) < J2ETA);
       eventMask += jet1jet2dphiC.addToMask(fabs(dphijj) < J1J2DPHI);
     } else if (nJet30a == 1) {
       eventMask += jet2etaC.get2ToId();
       eventMask += jet1jet2dphiC.get2ToId();
     }
     eventMask += eLooseVetoC.addToMask(nEle10V == 0);
     //eventMask += ntausC.addToMask(ntaus);
     eventMask += gammaLooseVetoC.addToMask(nGamma15V == 0);
     //eventMask += mumet200C.addToMask(mumet);
     for (Int_t i = 0; i <  metCut.size(); i++) {
       eventMask += mumetC[i].addToMask(metNoMu_pt > metCut[i]);
     }
     eventMask += oppChargeMuonsC.addToMask( (LepGood_pdgId[0] * LepGood_pdgId[1]) == -169);  // |pdgID| = 13 for muons
     eventMask += oppChargeLeptonsC.addToMask( (LepGood_pdgId[0] + LepGood_pdgId[1]) == 0);
     eventMask += twomuonsC.addToMask((fabs(LepGood_pdgId[0]) == 13) && (fabs(LepGood_pdgId[1]) == 13));
     eventMask += twomuLooseC.addToMask(nMu10V == 2);
     eventMask += muLooseVetoC.addToMask(nMu10V == 0);
     eventMask += mu1tightC.addToMask((LepGood_tightId[0] == 1) && (fabs(LepGood_pdgId[0]) == 13) && (LepGood_relIso04[0] < 0.12 ) && (fabs(LepGood_eta[0]) < 2.4) && (LepGood_pt[0] > 20));
     eventMask += mu1tightIdC.addToMask((LepGood_tightId[0] == 1) && (LepGood_relIso04[0] < 0.12 ) && (fabs(LepGood_pdgId[0]) == 13));
     eventMask += mu1ptC.addToMask((LepGood_pt[0] > 20) && (fabs(LepGood_pdgId[0]) == 13)); 
     eventMask += mu1etaC.addToMask( (fabs(LepGood_eta[0]) < 2.4) && (fabs(LepGood_pdgId[0]) == 13) );
     eventMask += twomuTightC.addToMask(nMu20T == 2);
     eventMask += mu2ptC.addToMask((LepGood_pt[1] > 10) && (fabs(LepGood_pdgId[1]) == 13));
     eventMask += mu2etaC.addToMask((fabs(LepGood_eta[1]) < 2.4) && (fabs(LepGood_pdgId[1]) == 13));
     eventMask += mumuInvMassC.addToMask((mZ1 > 60) && (mZ1 < 120));
     eventMask += genMuonsC.addToMask( myPartGenAlgo(nGenPart, GenPart_pdgId, GenPart_motherId, 13, 23) );
     eventMask += genTausC.addToMask( myPartGenAlgo(nGenPart, GenPart_pdgId, GenPart_motherId, 15, 23) );  

     zmumujetsControlSample.countEvents(eventMask,newwgt);

     // filling histogram with yields at the end of the selection in bins of met
     if ( ((eventMask & zmumujetsControlSample.globalMask.back()) == zmumujetsControlSample.globalMask.back()) ) {  
       // this histogram holds the final yields in bins of MET
       HzmumujetsYieldsMetBin->Fill(metNoMu_pt,newwgt);     
     }

     if ((metNoMu_pt > metBinEdges[0]) && (metNoMu_pt < metBinEdges[nMetBins])) {

       Int_t bin = myGetBin(metNoMu_pt,metBinEdges,nMetBins);
       //find step where cut on invariant mass is added the first time
       //Int_t index = zmumujetsControlSample.whichStepHas(mumuInvMassC.get2ToId()); 
       // if (((eventMask & zmumujetsControlSample.globalMask[index]) == zmumujetsControlSample.globalMask[index]) && (index < zmumujetsControlSample.getMaskSize())) {
       // 	 // this histogram holds the invariant mass distribution (one for each met bin)
       // 	 HmumuInvMass[bin]->Fill(mZ1,newwgt);       
       // }
       if (((eventMask & zmumujetsControlSample.globalMask.back()) == zmumujetsControlSample.globalMask.back())) {
	 // this histogram holds the invariant mass distribution (one for each met bin)
   	 HmumuInvMass[bin]->Fill(mZ1,newwgt);       
       }
              
     } 

   }

   cout<<endl;   
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &zmumujetsControlSample);
   cout<<endl;   

   for (Int_t i = 0; i < nMetBins; i++) {
     cout<<" mumet in ["<<metBinEdges[i]<<" , "<<metBinEdges[i+1]<<"] :     HmumuInvMass["<<i<<"]->GetSumOfWeights() = ";
     cout<<HmumuInvMass[i]->GetSumOfWeights()<<endl;
   }


   selection::saveYieldsAndEfficiency(nTotalWeightedEvents, &zmumujetsControlSample, yRow, eRow);

   string fname = "zmumujetsCSYields_" + suffix + ".txt";
   const char* filename = fname.c_str();   //this file is like a logbook
   //char answer = '\0';
   //answer = myAskToSaveFile(filename);
   char answer = 'y';

   if (answer == 'y') {

     ofstream myfile(filename,ios::out);

     if ( !myfile.is_open() ) {
       cout<<"Error: unable to open file "<<filename<<" !"<<endl;

     } else {
       //when writing on file, the date is printed as well unless an error occurs
       string command = "date>>" + fname;      
       if ( system(command.c_str()) != 0) cout<<"Error during \"system("<<command<<")\" call"<<endl;  

       myfile<<endl;      
       selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &zmumujetsControlSample);
       myfile<<endl;
       myfile.close();

     }
 
   }

   fname = "zmumujetsCSYields_" + suffix + ".tex";
   const char *texfname = fname.c_str();
   
   //answer = myAskToSaveFile(texfname);

   if (answer == 'y') {

     //creating a .tex file to build tables with data
     FILE *fp;	 
     //cout<<"name given to file: "<<texfname<<endl;
     if ( (fp=fopen(texfname,"w")) == NULL) {
       cout<<"Error: '"<<texfname<<"' not opened"<<endl;
     } else {
       cout<<"creating file '"<<texfname<<"' ..."<<endl;
       myAddDefaultPackages(fp,texfname);
       fprintf(fp,"\\begin{document}\n");
       fprintf(fp,"\n");
       string commentInTable;             
       commentInTable = "Note that cuts on second jet are applied only if a second jet exists with $p_t$ > 30\\,GeV.";
       makeTableTex(fp, LUMI, nTotalWeightedEvents, &zmumujetsControlSample);
       fprintf(fp,"\\end{document}\n");      
       fclose(fp);
     }

   }
   fname = "histZmumujetsInvMassCS_" + suffix + ".root";
   const char *histFileName = fname.c_str();
   cout<<"Saving histograms in file \""<<histFileName<<"\" ..."<<endl;
   TFile *histFile = new TFile(histFileName,"RECREATE");

   if (!histFile->IsOpen()) {

     cout<<"Error: file \""<<histFileName<<"\" was not opened."<<endl;

   } else {

     for (Int_t i = 0; i < nMetBins; i++) {
       HmumuInvMass[i]->Write();
     }

     histFile->Close();
     
   }

   delete histFile;

   fname = "histZmumujetsYieldsMetBinCS_" + suffix + ".root";
   const char *YieldsFileName = fname.c_str();
   cout<<"Saving histogram \""<<HzmumujetsYieldsMetBin->GetName()<<"\" in file \""<<YieldsFileName<<"\" ..."<<endl;
   TFile *YieldsFile = new TFile(YieldsFileName,"RECREATE");

   if (!YieldsFile->IsOpen()) cout<<"Error: file \""<<YieldsFileName<<"\" was not opened."<<endl;
   else HzmumujetsYieldsMetBin->Write();

   cout<<"MET < "<<metBinEdges[0]<<" : yield = ";
   cout<<HzmumujetsYieldsMetBin->GetBinContent(0)<<" +/- "<<HzmumujetsYieldsMetBin->GetBinError(0)<<endl;
   for (Int_t i = 0; i < nMetBins; i++ ) {  //from underflow to overflow bin 
     cout<<"MET in ["<<metBinEdges[i]<<","<<metBinEdges[i+1]<<"] : yield = ";
     cout<<HzmumujetsYieldsMetBin->GetBinContent(i+1)<<" +/- "<<HzmumujetsYieldsMetBin->GetBinError(i+1)<<endl;
   }
   cout<<"MET > "<<metBinEdges[nMetBins]<<" : yield = ";
   cout<<HzmumujetsYieldsMetBin->GetBinContent(nMetBins + 1)<<" +/- "<<HzmumujetsYieldsMetBin->GetBinError(nMetBins + 1)<<endl;

   YieldsFile->Close();
     
   delete YieldsFile;

   delete HzmumujetsYieldsMetBin;
   for (Int_t i = 0; i < nMetBins; i++) {
     delete HmumuInvMass[i];
   } 
   

}
