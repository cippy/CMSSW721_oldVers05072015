#define znunujetsAna_cxx
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
#include <TVector2.h>
#include <TVirtualFitter.h>
//my headers
#include "functionsForAnalysis.h"
#include "myClasses.h"

using namespace std;
using namespace myAnalyzerTEman;

#ifdef znunujetsAna_cxx

znunujetsAna::znunujetsAna(TTree *tree) : edimarcoTree(tree) {
  //cout <<"check in constructor "<<endl;
  Init(tree);

}

#endif

void znunujetsAna::loop(const char* configFileName)
{

   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);  
   // warning: in Emanuele's trees non integer values are float

   fChain->SetBranchStatus("weight",1);   // includes k-factor

   fChain->SetBranchStatus("nMu10V",1);  // # of muons passing loose selection for muon veto
   fChain->SetBranchStatus("nEle10V",1);  // # of electrons passing loose selection for electron veto
   fChain->SetBranchStatus("nTau15V",1);
   fChain->SetBranchStatus("nGamma15V",1);  // # of photons passing loose selection for photon veto

   fChain->SetBranchStatus("dphijj",1);          // dphi between 1st and 2nd jet, 999 if second jet doesn't exist
   fChain->SetBranchStatus("jetclean",1);      // 1 if jet is cleaned, 0 otherwise
   fChain->SetBranchStatus("nJet30",1);         // # of jets with pt > 30 && |eta| < 2.4
   fChain->SetBranchStatus("nJet",1);         // # of jets with pt > 25 && |eta| < 2.5
   fChain->SetBranchStatus("Jet_pt",1);  
   fChain->SetBranchStatus("Jet_eta",1); 

   // fChain->SetBranchStatus("met_pt",1);
   // fChain->SetBranchStatus("met_phi",1);
   fChain->SetBranchStatus("metNoMu_pt",1);   // likely this will coincide with the pt of the Z(nunu)
   
   char ROOT_FNAME[50];
   char TXT_FNAME[50];
   char TEX_FNAME[50]; 

   strcpy(ROOT_FNAME,"znunujetsAna.root");
   strcpy(TXT_FNAME,"znunujetsAna.txt");
   strcpy(TEX_FNAME,"znunujetsAna.tex");
   
   Double_t LUMI;
   Int_t NJETS;
   Double_t J1PT;
   Double_t J1ETA;
   Double_t J2PT;
   Double_t J2ETA;
   Double_t J1J2DPHI;
   Int_t LEP_PDG_ID;

   ifstream inputFile(configFileName);

   if (inputFile.is_open()) {

     Double_t value;
     string parameterName;
     vector<Double_t> parameterValue;

     mySpaces(cout,2);
     cout << "Printing content of " << configFileName << " file" << endl;
     mySpaces(cout,1);

     while (inputFile >> parameterName >> value) {

       parameterValue.push_back(value);
       cout << setw(20) << parameterName << setw(7) << value << endl;

     }
     
     // following variables are initialized with values in the file configFileName
     LUMI = parameterValue[0];
     NJETS = (Int_t) parameterValue[1];
     J1PT = parameterValue[2];
     J1ETA = parameterValue[3];
     J2PT = parameterValue[4];
     J2ETA = parameterValue[5];
     J1J2DPHI = parameterValue[6];
     LEP_PDG_ID = (Int_t) parameterValue[7];
     mySpaces(cout,2);

     inputFile.close();
                                                                                                                         
   } else {

     cout << "Error: could not open file " << configFileName << endl;
     exit(EXIT_FAILURE);

   }

   //Double_t metBinEdges[] = {200., 250., 300., 350., 400., 500., 650., 1000.};
   Double_t metBinEdges[] = {200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 750., 850., 1000.};
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
   selection njetsEmanC("njetsEmanC","njets","1 or 2 jets, cleaning, pt > 30, |eta| < 2.4");       // using nJet30
   selection jet1ptC("jet1ptC",Form("jet1pt > %4.0lf",(Double_t)J1PT));
   selection jet1etaC("jet1etaC",Form("|jet1eta| < %2.1lf",J1ETA));
   selection jet2etaC("jet2etaC",Form("|jet2eta| < %2.1lf",J2ETA),Form("only if njets = %i",NJETS));
   selection jjdphiEmanC("jjdphiEmanC",Form("jjdphi < %1.1lf",J1J2DPHI),Form("only if njets = %i",NJETS));
   selection eLooseVetoC("eLooseVetoC","electrons veto");
   selection muLooseVetoC("muLooseVeto","muons veto");
   selection gammaLooseVetoC("gammaLooseVetoC","photons veto");
   selection tauLooseVetoC("tauLooseVetoC","taus veto");
   selection mumet200C("mumet200C","mumet > 200"); 
   
   selection::checkMaskLength();
   selection::printActiveSelections(cout);

   mask znunujetsMonojetSelection("Z->mumu control sample with selection flow as Emanuele's");
   znunujetsMonojetSelection.append(mumet200C.get2ToId());
   znunujetsMonojetSelection.append(njetsEmanC.get2ToId());
   znunujetsMonojetSelection.append(jet1ptC.get2ToId());
   znunujetsMonojetSelection.append(jjdphiEmanC.get2ToId());
   znunujetsMonojetSelection.append(eLooseVetoC.get2ToId());
   znunujetsMonojetSelection.append(muLooseVetoC.get2ToId());
   znunujetsMonojetSelection.append(tauLooseVetoC.get2ToId());
   znunujetsMonojetSelection.append(gammaLooseVetoC.get2ToId());

   mask znunujetsMonojetSelectionNoMuVeto("Z->mumu control sample with selection flow as Emanuele's, no Mu veto");
   znunujetsMonojetSelectionNoMuVeto.append(mumet200C.get2ToId());
   znunujetsMonojetSelectionNoMuVeto.append(njetsEmanC.get2ToId());
   znunujetsMonojetSelectionNoMuVeto.append(jet1ptC.get2ToId());
   znunujetsMonojetSelectionNoMuVeto.append(jjdphiEmanC.get2ToId());
   znunujetsMonojetSelectionNoMuVeto.append(eLooseVetoC.get2ToId());
   znunujetsMonojetSelectionNoMuVeto.append(tauLooseVetoC.get2ToId());
   znunujetsMonojetSelectionNoMuVeto.append(gammaLooseVetoC.get2ToId());

   mask znunujetsMonojetSelectionNoEleVeto("Z->mumu control sample with selection flow as Emanuele's, no Ele veto");
   znunujetsMonojetSelectionNoEleVeto.append(mumet200C.get2ToId());
   znunujetsMonojetSelectionNoEleVeto.append(njetsEmanC.get2ToId());
   znunujetsMonojetSelectionNoEleVeto.append(jet1ptC.get2ToId());
   znunujetsMonojetSelectionNoEleVeto.append(jjdphiEmanC.get2ToId());
   znunujetsMonojetSelectionNoEleVeto.append(muLooseVetoC.get2ToId());
   znunujetsMonojetSelectionNoEleVeto.append(tauLooseVetoC.get2ToId());
   znunujetsMonojetSelectionNoEleVeto.append(gammaLooseVetoC.get2ToId());

   cout << "Opening file " <<ROOT_FNAME<< endl;

   TFile *rootFile = new TFile(ROOT_FNAME,"RECREATE");
   if (!rootFile || !rootFile->IsOpen()) {
     cout<<"Error: file \""<<ROOT_FNAME<<"\" was not opened."<<endl;
     exit(EXIT_FAILURE);
   }

   Double_t nTotalWeightedEvents = 0.0;    // total events (including weights)

   TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 
   TH1::StatOverflows();                 //enable use of underflows and overflows for statistics computation 
   TVirtualFitter::SetDefaultFitter("Minuit");

   TH1D *HzlljetsYieldsMetBin = new TH1D("HzlljetsYieldsMetBin","yields of Z#nu#nu control sample in bins of met;#slash{E}_{T};# of events",nMetBins,metBinEdges);
   TH1D *HzlljetsYieldsMetBinNoMuVeto = new TH1D("HzlljetsYieldsMetBinNoMuVeto","yields of Z#nu#nu control sample in bins of met, no mu veto;#slash{E}_{T};# of events",nMetBins,metBinEdges);
   TH1D *HzlljetsYieldsMetBinNoEleVeto = new TH1D("HzlljetsYieldsMetBinNoEleVeto","yields of Z#nu#nu control sample in bins of met, no ele veto;#slash{E}_{T};# of events",nMetBins,metBinEdges);

   // saving histograms with bin edges of other histograms used (e.g. content of metBinEdges array ...)
   TH1D *HmetBinEdges = new TH1D("HmetBinEdges","bin edges for met distributions",nMetBins+1,0.0,nMetBins+1);
   for (Int_t i = 0; i <= nMetBins; i++) {
     HmetBinEdges->SetBinContent(i+1,metBinEdges[i]);
   }

   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"znunujetsAna::loop()"<<endl;
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

     eventMask += njetsEmanC.addToMask( (nJet30 == 1 || nJet30 == 2) && jetclean > 0.5);
     eventMask += jjdphiEmanC.addToMask( nJet30 == 1 || (nJet == 2 && abs(dphijj) < J1J2DPHI));
     eventMask += jet1ptC.addToMask(Jet_pt[0] > J1PT);               
     eventMask += eLooseVetoC.addToMask(nEle10V == 0);
     eventMask += muLooseVetoC.addToMask(nMu10V == 0);
     eventMask += tauLooseVetoC.addToMask(nTau15V == 0);
     eventMask += gammaLooseVetoC.addToMask(nGamma15V == 0);
     eventMask += mumet200C.addToMask(metNoMu_pt > 200);
     for (Int_t i = 0; i <  metCut.size(); i++) {
       eventMask += mumetC[i].addToMask(metNoMu_pt > metCut[i]);
     }
     
     znunujetsMonojetSelection.countEvents(eventMask, newwgt);
     znunujetsMonojetSelectionNoMuVeto.countEvents(eventMask, newwgt);
     znunujetsMonojetSelectionNoEleVeto.countEvents(eventMask, newwgt);

     // filling histogram with yields at the end of the selection in bins of met
     if ( ((eventMask & znunujetsMonojetSelection.globalMask.back()) == znunujetsMonojetSelection.globalMask.back()) ) {  
       // this histogram holds the final yields in bins of MET
       HzlljetsYieldsMetBin->Fill(metNoMu_pt,newwgt);     
     }
     if ( ((eventMask & znunujetsMonojetSelectionNoMuVeto.globalMask.back()) == znunujetsMonojetSelectionNoMuVeto.globalMask.back()) ) {  
       // this histogram holds the final yields in bins of MET
       HzlljetsYieldsMetBinNoMuVeto->Fill(metNoMu_pt,newwgt);     
     }
     if ( ((eventMask & znunujetsMonojetSelectionNoEleVeto.globalMask.back()) == znunujetsMonojetSelectionNoEleVeto.globalMask.back()) ) {  
       // this histogram holds the final yields in bins of MET
       HzlljetsYieldsMetBinNoEleVeto->Fill(metNoMu_pt,newwgt);     
     }
     
   }

   mySpaces(cout,2);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &znunujetsMonojetSelection);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &znunujetsMonojetSelectionNoMuVeto);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &znunujetsMonojetSelectionNoEleVeto);
   cout<<endl;   

   mySpaces(cout,2);
   myPrintYieldsMetBinInStream(cout, HzlljetsYieldsMetBin, metBinEdges, nMetBins);

   cout<<"creating file '"<<TXT_FNAME<<"' ..."<<endl;
   ofstream myfile(TXT_FNAME,ios::out);

   if ( !myfile.is_open() ) {

     cout<<"Error: unable to open file "<<TXT_FNAME<<" !"<<endl;
     exit(EXIT_FAILURE);
     
   }
      
   //selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &lep_Acc_Eff);      
   selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &znunujetsMonojetSelection);
   selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &znunujetsMonojetSelectionNoMuVeto);
   selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &znunujetsMonojetSelectionNoEleVeto);
   mySpaces(myfile,2);
   myPrintYieldsMetBinInStream(myfile, HzlljetsYieldsMetBin, metBinEdges, nMetBins);

   rootFile->Write();

   rootFile->Close();
   delete rootFile;

   //creating a .tex file to build tables with data
   FILE *fp;
   fp = fopen(TEX_FNAME,"w");

   if ( fp == NULL)  cout<<"Error: '"<<TEX_FNAME<<"' not opened"<<endl;
   else {

     cout<<"creating file '"<<TEX_FNAME<<"' ..."<<endl;
     myAddDefaultPackages(fp,TEX_FNAME);
     fprintf(fp,"\\begin{document}\n");
     fprintf(fp,"\n");
     string commentInTable;       
     //makeTableTex(fp, LUMI, nTotalWeightedEvents, &mu_Acc_Eff, commentInTable);      
     commentInTable = "Note that cuts on second jet are applied only if a second jet exists with $p_t$ > 30\\,GeV.";
     makeTableTex(fp, LUMI, nTotalWeightedEvents, &znunujetsMonojetSelection,commentInTable);
     makeTableTex(fp, LUMI, nTotalWeightedEvents, &znunujetsMonojetSelectionNoMuVeto,commentInTable);
     makeTableTex(fp, LUMI, nTotalWeightedEvents, &znunujetsMonojetSelectionNoEleVeto,commentInTable);
     fprintf(fp,"\\end{document}\n");      
     fclose(fp);

   }

   // end of tex file


}
