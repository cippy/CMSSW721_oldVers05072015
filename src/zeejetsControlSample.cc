#define zeejetsControlSample_cxx
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

#ifdef zeejetsControlSample_cxx

zeejetsControlSample::zeejetsControlSample(TTree *tree, const char* inputSuffix) : edimarcoTree(tree) {
  //cout <<"check in constructor "<<endl;
  suffix = inputSuffix;
  Init(tree);

}

#endif

#define LUMI 5                         // integrated luminosity in fb^-1 to which yields are normalized
#define ELE_PDGID 11
#define HLT_ELECTRONS 1   // flag to include HLT in selection: if 1, trigger is applied

void zeejetsControlSample::loop(vector< Double_t > &yRow, vector< Double_t > &eRow)
{

   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);  
   // warning: in Emanuele's trees non integer values are float

   fChain->SetBranchStatus("weight",1);   // includes k-factor

   fChain->SetBranchStatus("nMu10V",1);  // # of muons passing loose selection
   fChain->SetBranchStatus("nEle10V",1);  // # of electrons passing loose selection for electron veto
   fChain->SetBranchStatus("nGamma15V",1);  // # of photons passing loose selection for photon veto
   fChain->SetBranchStatus("nMu20T",1);  // # of muons passing tight selection (pt > 20 + everything else)
   fChain->SetBranchStatus("nEle20T",1);  // # of electrons passing tight selection (pt > 20 + everything else)

   fChain->SetBranchStatus("dphijj",1);          // dphi between 1st and 2nd jet, 999 if second jet doesn't exist
   fChain->SetBranchStatus("jetclean",1);      // 1 if jet is cleaned, 0 otherwise
   fChain->SetBranchStatus("nJet",1);         // # of jets with pt > 25 && |eta| < 2.5
   fChain->SetBranchStatus("nJet30",1);         // # of jets with pt > 30 && |eta| < 2.4
   fChain->SetBranchStatus("nJet30a",1);       // # of jets with pt > 30 && |eta| < 4.7 
   fChain->SetBranchStatus("Jet_pt",1);  
   fChain->SetBranchStatus("Jet_eta",1);  
 
   fChain->SetBranchStatus("nLepGood",1);
   fChain->SetBranchStatus("LepGood_pdgId",1);  // must be 13 for muons ( -13 for mu+), 11 for electrons and 15 for taus
   fChain->SetBranchStatus("LepGood_pt",1);
   fChain->SetBranchStatus("LepGood_phi",1);
   fChain->SetBranchStatus("LepGood_eta",1);
   //fChain->SetBranchStatus("LepGood_charge",1);
   fChain->SetBranchStatus("LepGood_tightId",1);
   fChain->SetBranchStatus("LepGood_relIso04",1);
   fChain->SetBranchStatus("ngenLep",1);
   fChain->SetBranchStatus("genLep_pdgId",1);
   //fChain->SetBranchStatus("m2l",1);  // m(ll)  (I can compute it myself, maybe it's better)
   fChain->SetBranchStatus("mZ1",1);  // best m(ll) SF/OS

   fChain->SetBranchStatus("LepGood_pdgId",1);  // must be 13 for muons ( -13 for mu+), 11 for electrons and 15 for taus
   fChain->SetBranchStatus("LepGood_pt",1);
   fChain->SetBranchStatus("LepGood_phi",1);
   fChain->SetBranchStatus("LepGood_eta",1);
   //fChain->SetBranchStatus("LepGood_charge",1);
   fChain->SetBranchStatus("LepGood_tightId",1);
   fChain->SetBranchStatus("LepGood_relIso04",1);
   fChain->SetBranchStatus("genLep_pdgId",1);
   //fChain->SetBranchStatus("m2l",1);  // m(ll)  (I can compute it myself, maybe it's better)
   fChain->SetBranchStatus("mZ1",1);  // best m(ll) SF/OS

   fChain->SetBranchStatus("met_pt",1);
   fChain->SetBranchStatus("met_phi",1); 

   Double_t metBinEdges[] = {200., 250., 300., 350., 400., 500., 650., 1000.};
   Int_t nMetBins = (sizeof(metBinEdges)/sizeof(Double_t)) - 1;

   vector<Double_t> metCut;
   metCut.push_back(250);
   metCut.push_back(300);
   metCut.push_back(350);
   metCut.push_back(400);
   metCut.push_back(500);
   selection elemetC[metCut.size()];
   for (Int_t i = 0; i < metCut.size(); i++) {
     elemetC[i].set(Form("elemetC[%i]",i),Form("elemet > %3.0lf",metCut.at(i)));
   }
   selection njetsC("njetsC",Form("njets <= %i",NJETS),"pt > 30; |eta| < 4.7");   // using nJet30a
   selection njetsEmanC("njetsEmanC","njets","1 or 2 jets, cleaning, pt > 30, |eta| < 2.4");       // using nJet30
   selection jet1ptC("jet1ptC",Form("jet1pt > %4.0lf",(Double_t)J1PT));
   selection jet1etaC("jet1etaC",Form("|jet1eta| < %2.1lf",J1ETA));
   // jetclean is 1 if cleaning is passed, 0 otherwise. It's applied to first jet and , if any, to the second
   selection jetCleaningC("jetCleaningC","jet cleaning");           
   selection jet2etaC("jet2etaC",Form("|jet2eta| < %2.1lf",J2ETA),"only if njets = 2");
   selection jet1jet2dphiC("jet1jet2dphiC",Form("|jet1jet2dphi| < %1.1lf",J1J2DPHI),"only if njets = 2");
   selection jjdphiEmanC("jjdphiEmanC","jjdphi < 2.5","only if var \"nJet\" = 2 (not nJet30)");
   selection muLooseVetoC("muLooseVeto","muons veto");
   selection gammaLooseVetoC("gammaLooseVetoC","photons veto");
   selection elemet200C("elemet200C","elemet > 200"); 
   // additional selections for control sample
   selection oppChargeElectronsC("oppChargeElectronsC","OS electrons");
   selection oppChargeLeptonsC("oppChargeLeptonsC","OS/SF leptons");
   selection twoelectronsC("twoelectronsC","2 electrons");
   selection eLooseVetoC("eLooseVetoC","electrons veto");
   selection twoeleLooseC("twoeleLooseC","2 loose electrons");
   // selection neleTightIdC("neleTightIdC","nTightElectrons",">",0,"tight ID");
   // selection neleTightIdForEfficiencyC("neleTightIdForEfficiencyC","nTightElectrons","=",2,"tight ID");
   //selection ele1tightC("ele1tightC","leading electron tight");   // with Emanuele's tree, it means both tight id and pt > 20 and |eta| < 2.4 (see friend tree)
   selection ele1tightIdC("ele1tightIdC","leading electron tight","tight ID + relIso04 (as Emanuele)");
   selection ele2tightIdC("ele2TightIdC","trailing electron tight","tight ID + relIso04 (as Emanuele)");
   selection twoeleTightC("twoeleTightC","2 tight electrons");
   selection ele1ptC("ele1ptC","ele1pt > 32","leading pt electron");
   selection ele2ptC("ele2ptC","ele2Pt > 20");
   selection ele1etaC("ele1etaC","|ele1eta| < 2.4","leading pt electron");  
   selection ele2etaC("ele2etaC","|ele2eta| < 2.4");
   selection eeInvMassC("eeInvMassC","ee mass in [60,120]");
   //selection eeHLTselectionC("eeHLTselectionC","ee HLT selection","2 tightID electrons, pt1 > 30, pt2 > 18, |eta| < 2.5");
   // selection genTausC("genTausC","taus generated");               // 11, 13, 15 for e, mu, tau  ( - sign for antiparticles)              
   // selection genElectronsC("genElectronsC","electrons generated");            // 11, 13, 15 for e, mu, tau  ( - sign for antiparticles)  

   if (HLT_ELECTRONS) cout << " ==== HLT applied ==== " << endl;

   selection::checkMaskLength();
   //selection::printActiveSelections(cout);

   mask zeejetsControlSample("Z->ee control sample with selection flow as Emanuele's");
   zeejetsControlSample.append(twoeleLooseC.get2ToId() + oppChargeLeptonsC.get2ToId());
   zeejetsControlSample.append(twoelectronsC.get2ToId());
   zeejetsControlSample.append(ele1tightIdC.get2ToId() + ele2tightIdC.get2ToId());
   zeejetsControlSample.append(eeInvMassC.get2ToId());
   zeejetsControlSample.append(elemet200C.get2ToId());
   zeejetsControlSample.append(njetsEmanC.get2ToId());
   zeejetsControlSample.append(jet1ptC.get2ToId());
   zeejetsControlSample.append(jjdphiEmanC.get2ToId());
   zeejetsControlSample.append(muLooseVetoC.get2ToId());
   zeejetsControlSample.append(gammaLooseVetoC.get2ToId());

   Double_t nTotalWeightedEvents = 0.0;    // total events (including weights)

   TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 
   TH1::StatOverflows();                 //enable use of underflows and overflows for statistics computation 
   //TVirtualFitter::SetDefaultFitter("Minuit");

   Int_t Hcolor[] = {1,2,3,4,5,6,7,8,9,12,18,30,38,41,42,46,47,49};    

   TH1D *HeeInvMass[nMetBins];
   for (Int_t i = 0; i < nMetBins; i++) {
     HeeInvMass[i] = new TH1D(Form("HeeInvMass[%i]",i),"",30,60.,120.);
   } 

   TH1D *HzeejetsYieldsMetBin = new TH1D("HzeejetsYieldsMetBin","yields of Zee control sample in bins of met;MET;# of events",
                                            nMetBins,metBinEdges);

   TVector2 e1, e2, met, ele, eleVectorSum;  // ele is any electron (mu1 is the first and mu2 is the second valid in the list sorted by pt)
   // following indices refer to the leading pair of OS/SF in the list of LepGood. They are initialized with 0 and 1 by default, but can be set with function
   // myGetPairIndexInArray (see functionsForAnalysis.cc for reference). 
   // When myGetPairIndexInArray() is called, the index of "correct" particles will be used. If they are not found (e.g. a pair of OS/SF is mismeasured as 2 e+), 
   // indices are set as 0 and 1 (and following selection asking lep[0] and lep[1] to be OS or whatever will fail).
   Int_t firstIndex = 0;
   Int_t secondIndex = 1;

   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"zeejetsControlSample::loop()"<<endl;
   cout<<"nentries = "<<nentries<<endl;   

   Long64_t nbytes = 0, nb = 0;
   for (Int_t jentry=0; jentry<nentries; jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;   

     // I find the indices corresponding to the 2 leading lepton
     //cout<<"entry : "<<jentry<<endl;
     myGetPairIndexInArray(ELE_PDGID, nLepGood, LepGood_pdgId, firstIndex, secondIndex);

     //analysis is performed on events that passed the HLT trigger for electrons (the flag HLT_ELECTRONS must also be set to 1)
     if ( !( HLT_ELECTRONS && (LepGood_tightId[firstIndex] == 1) && (fabs(LepGood_pdgId[firstIndex]) == ELE_PDGID) && 
	  (LepGood_tightId[secondIndex] == 1) && (fabs(LepGood_pdgId[secondIndex]) == ELE_PDGID) && 
	  ( (LepGood_pdgId[firstIndex] + LepGood_pdgId[secondIndex]) == 0) &&
	  (fabs(LepGood_eta[firstIndex]) < 2.5) && (fabs(LepGood_eta[secondIndex]) < 2.5) && 
	     (LepGood_pt[secondIndex] > 30) && (LepGood_pt[secondIndex] > 18) ) ) continue;

     UInt_t eventMask = 0; 
     Double_t newwgt = weight * LUMI;

     nTotalWeightedEvents += newwgt;  // counting events with weights

     e1.SetMagPhi(LepGood_pt[firstIndex],LepGood_phi[firstIndex]);
     e2.SetMagPhi(LepGood_pt[secondIndex],LepGood_phi[secondIndex]);
     met.SetMagPhi(met_pt,met_phi);
     eleVectorSum.SetMagPhi(0.0,0.0);   // for each event it must be initialized to 0

     for (Int_t i = 0; i < nLepGood; i++) {
       if (fabs(LepGood_pdgId[i]) == ELE_PDGID) {
	 ele.SetMagPhi(LepGood_pt[i],LepGood_phi[i]);
	 eleVectorSum += ele;
       }
     }

     Double_t metNoEle = (e1 + e2 + met).Mod();

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
     eventMask += muLooseVetoC.addToMask(nMu10V == 0);
     eventMask += gammaLooseVetoC.addToMask(nGamma15V == 0);
     eventMask += elemet200C.addToMask(metNoEle > 200);
     for (Int_t i = 0; i <  metCut.size(); i++) {
       eventMask += elemetC[i].addToMask(metNoEle > metCut[i]);
     }
     eventMask += oppChargeElectronsC.addToMask( (LepGood_pdgId[firstIndex] * LepGood_pdgId[secondIndex]) == -121);  // |pdgID| = 11 for electrons
     eventMask += oppChargeLeptonsC.addToMask( (LepGood_pdgId[firstIndex] + LepGood_pdgId[secondIndex]) == 0);
     eventMask += twoelectronsC.addToMask((fabs(LepGood_pdgId[firstIndex]) == ELE_PDGID) && (fabs(LepGood_pdgId[secondIndex]) == ELE_PDGID));
     eventMask += twoeleLooseC.addToMask(nEle10V == 2);    
     //eventMask += ele1tightC.addToMask((LepGood_tightId[firstIndex] == 1) && (fabs(LepGood_pdgId[firstIndex]) == ELE_PDGID) && 
     //(LepGood_relIso04[firstIndex] < 0.12 ) && (fabs(LepGood_eta[firstIndex]) < 2.4) && (LepGood_pt[firstIndex] > 32));
     eventMask += ele1tightIdC.addToMask((LepGood_tightId[firstIndex] == 1) && (LepGood_relIso04[firstIndex] < 0.12 ) && 
					 (fabs(LepGood_pdgId[firstIndex]) == ELE_PDGID));
     eventMask += ele2tightIdC.addToMask((LepGood_tightId[secondIndex] == 1) && (LepGood_relIso04[secondIndex] < 0.12 ) && 
					 (fabs(LepGood_pdgId[secondIndex]) == ELE_PDGID));
     eventMask += ele1ptC.addToMask((LepGood_pt[firstIndex] > 32) && (fabs(LepGood_pdgId[firstIndex]) == ELE_PDGID)); 
     eventMask += ele1etaC.addToMask( (fabs(LepGood_eta[firstIndex]) < 2.4) && (fabs(LepGood_pdgId[firstIndex]) == ELE_PDGID) );
     eventMask += twoeleTightC.addToMask(nEle20T == 2);
     eventMask += ele2ptC.addToMask((LepGood_pt[secondIndex] > 20) && (fabs(LepGood_pdgId[secondIndex]) == ELE_PDGID));
     eventMask += ele2etaC.addToMask((fabs(LepGood_eta[secondIndex]) < 2.4) && (fabs(LepGood_pdgId[secondIndex]) == ELE_PDGID));
     eventMask += eeInvMassC.addToMask((mZ1 > 60) && (mZ1 < 120));
     // eventMask += genElectronsC.addToMask( myPartGenAlgo(nGenPart, GenPart_pdgId, GenPart_motherId, ELE_PDGID, 23) );
     // eventMask += genTausC.addToMask( myPartGenAlgo(nGenPart, GenPart_pdgId, GenPart_motherId, 15, 23) );  

     zeejetsControlSample.countEvents(eventMask,newwgt);

     // filling histogram with yields at the end of the selection in bins of met
     if ( ((eventMask & zeejetsControlSample.globalMask.back()) == zeejetsControlSample.globalMask.back()) ) {  
       // this histogram holds the final yields in bins of MET
       HzeejetsYieldsMetBin->Fill(metNoEle,newwgt);     
     }

     if ((metNoEle > metBinEdges[0]) && (metNoEle < metBinEdges[nMetBins])) {

       Int_t bin = myGetBin(metNoEle,metBinEdges,nMetBins);
       if (((eventMask & zeejetsControlSample.globalMask.back()) == zeejetsControlSample.globalMask.back()) ) {
	 // this histogram holds the invariant mass distribution (one for each met bin)
   	 HeeInvMass[bin]->Fill(mZ1,newwgt);       
       }
              
     } 

   }

   cout<<endl;   
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &zeejetsControlSample);
   cout<<endl;   

   for (Int_t i = 0; i < nMetBins; i++) {
     cout<<" elemet in ["<<metBinEdges[i]<<" , "<<metBinEdges[i+1]<<"] :     HeeInvMass["<<i<<"]->GetSumOfWeights() = ";
     cout<<HeeInvMass[i]->GetSumOfWeights()<<endl;
   }


   selection::saveYieldsAndEfficiency(nTotalWeightedEvents, &zeejetsControlSample, yRow, eRow);

   string fname = "zeejetsCSYields_" + suffix + ".txt";
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
       selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &zeejetsControlSample);
       myfile<<endl;
       myfile.close();

     }
 
   }

   fname = "zeejetsCSYields_" + suffix + ".tex";
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
       makeTableTex(fp, LUMI, nTotalWeightedEvents, &zeejetsControlSample);
       fprintf(fp,"\\end{document}\n");      
       fclose(fp);
     }

   }
   fname = "histZeejetsInvMassCS_" + suffix + ".root";
   const char *histFileName = fname.c_str();
   cout<<"Saving histograms in file \""<<histFileName<<"\" ..."<<endl;
   TFile *histFile = new TFile(histFileName,"RECREATE");

   if (!histFile->IsOpen()) {

     cout<<"Error: file \""<<histFileName<<"\" was not opened."<<endl;

   } else {

     for (Int_t i = 0; i < nMetBins; i++) {
       HeeInvMass[i]->Write();
     }

     histFile->Close();
     
   }

   delete histFile;

   fname = "histZeejetsYieldsMetBinCS_" + suffix + ".root";
   const char *YieldsFileName = fname.c_str();
   cout<<"Saving histogram \""<<HzeejetsYieldsMetBin->GetName()<<"\" in file \""<<YieldsFileName<<"\" ..."<<endl;
   TFile *YieldsFile = new TFile(YieldsFileName,"RECREATE");

   if (!YieldsFile->IsOpen()) cout<<"Error: file \""<<YieldsFileName<<"\" was not opened."<<endl;
   else HzeejetsYieldsMetBin->Write();

   cout<<"MET < "<<metBinEdges[0]<<" : yield = ";
   cout<<HzeejetsYieldsMetBin->GetBinContent(0)<<" +/- "<<HzeejetsYieldsMetBin->GetBinError(0)<<endl;
   for (Int_t i = 0; i < nMetBins; i++ ) {  //from underflow to overflow bin 
     cout<<"MET in ["<<metBinEdges[i]<<","<<metBinEdges[i+1]<<"] : yield = ";
     cout<<HzeejetsYieldsMetBin->GetBinContent(i+1)<<" +/- "<<HzeejetsYieldsMetBin->GetBinError(i+1)<<endl;
   }
   cout<<"MET > "<<metBinEdges[nMetBins]<<" : yield = ";
   cout<<HzeejetsYieldsMetBin->GetBinContent(nMetBins + 1)<<" +/- "<<HzeejetsYieldsMetBin->GetBinError(nMetBins + 1)<<endl;

   YieldsFile->Close();
     
   delete YieldsFile;

   delete HzeejetsYieldsMetBin;
   for (Int_t i = 0; i < nMetBins; i++) {
     delete HeeInvMass[i];
   } 
   

}
