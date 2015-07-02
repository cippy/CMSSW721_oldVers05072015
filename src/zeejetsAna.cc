#define zeejetsAna_cxx
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
#include "thresholds.h"
//#include "edimarcoTreeFriend.h"

using namespace std;
using namespace myAnalyzerTEman;

#ifdef zeejetsAna_cxx

zeejetsAna::zeejetsAna(TTree *tree) : edimarcoTree(tree) {
  //cout <<"check in constructor "<<endl;
  Init(tree);

}

#endif

#define LUMI 5                         // integrated luminosity in fb^-1 to which yields are normalized
#define ELE_PDGID 11
#define HLT_ELECTRONS 1   // flag to include HLT in selection: if 1, trigger is applied

void zeejetsAna::loop()
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

   fChain->SetBranchStatus("nGenPart",1);
   fChain->SetBranchStatus("GenPart_pdgId",1);
   fChain->SetBranchStatus("GenPart_motherId",1);
   fChain->SetBranchStatus("GenPart_pt",1);
   fChain->SetBranchStatus("GenPart_phi",1);
   //fChain->SetBranchStatus("GenPart_motherIndex",1);

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
   selection ele1tightIdC("ele1tightIdC","leading electron tight","tight ID + relIso04 (as Emanuele)");
   selection ele2tightIdC("ele2TightIdC","trailing electron tight","tight ID + relIso04 (as Emanuele)");
   selection twoeleTightC("twoeleTightC","2 tight electrons");
   selection ele1ptC("ele1ptC","ele1pt > 32","leading pt electron");
   selection ele2ptC("ele2ptC","ele2Pt > 20");
   selection ele1etaC("ele1etaC","|ele1eta| < 2.4","leading pt electron");  
   selection ele2etaC("ele2etaC","|ele2eta| < 2.4");
   selection eeInvMassC("eeInvMassC","ee mass in [60,120]");
   //selection eeHLTselectionC("eeHLTselectionC","ee HLT selection","2 tightID electrons, pt1 > 30, pt2 > 18, |eta| < 2.5");
   selection genTausC("genTausC","taus generated");               // 11, 13, 15 for e, mu, tau  ( - sign for antiparticles)              
   selection genElectronsC("genElectronsC","electrons generated");            // 11, 13, 15 for e, mu, tau  ( - sign for antiparticles)  

   if (HLT_ELECTRONS) cout << " ==== HLT applied ==== " << endl;

   selection::checkMaskLength();
   selection::printActiveSelections(cout);
   
   UInt_t maskMonoJetSelection = njetsEmanC.get2ToId() + jet1ptC.get2ToId() + jjdphiEmanC.get2ToId() +
                                   muLooseVetoC.get2ToId() + gammaLooseVetoC.get2ToId() + elemet200C.get2ToId();

   UInt_t maskElectronAcceptance = ele1ptC.get2ToId() + ele2ptC.get2ToId() + ele1etaC.get2ToId() + ele2etaC.get2ToId() +
                                    eeInvMassC.get2ToId(); 

   UInt_t maskElectronEfficiency = ele1tightIdC.get2ToId() + ele2tightIdC.get2ToId(); 

   UInt_t maskEleAxe = ele1ptC.get2ToId() + ele2ptC.get2ToId() + ele1etaC.get2ToId() + ele2etaC.get2ToId() +
                                            eeInvMassC.get2ToId() + twoeleLooseC.get2ToId() + ele1tightIdC.get2ToId(); 

   UInt_t maskEleGen = genElectronsC.get2ToId();
   UInt_t maskTauGen = genTausC.get2ToId();

   mask ele_Acc_Eff("ee acceptance and efficiency"); 
   ele_Acc_Eff.append(maskEleGen);
   ele_Acc_Eff.append(maskMonoJetSelection);
   ele_Acc_Eff.append(maskElectronAcceptance);
   ele_Acc_Eff.append(maskElectronEfficiency);

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

   mask zeejetsControlSampleGenEle("Z->ee control sample (e generated) with selection flow as Emanuele's");
   zeejetsControlSampleGenEle.append(maskEleGen);
   zeejetsControlSampleGenEle.append(twoeleLooseC.get2ToId() + oppChargeLeptonsC.get2ToId());
   zeejetsControlSampleGenEle.append(twoelectronsC.get2ToId());
   zeejetsControlSampleGenEle.append(ele1tightIdC.get2ToId() + ele2tightIdC.get2ToId());
   zeejetsControlSampleGenEle.append(eeInvMassC.get2ToId());
   zeejetsControlSampleGenEle.append(elemet200C.get2ToId());
   zeejetsControlSampleGenEle.append(njetsEmanC.get2ToId());
   zeejetsControlSampleGenEle.append(jet1ptC.get2ToId());
   zeejetsControlSampleGenEle.append(jjdphiEmanC.get2ToId());
   zeejetsControlSampleGenEle.append(muLooseVetoC.get2ToId());
   zeejetsControlSampleGenEle.append(gammaLooseVetoC.get2ToId());

   mask tautaubkgInZee("tau tau background in Z->ee control sample");
   tautaubkgInZee.append(maskTauGen);
   tautaubkgInZee.append(twoeleLooseC.get2ToId() + oppChargeLeptonsC.get2ToId());
   tautaubkgInZee.append(twoelectronsC.get2ToId());
   tautaubkgInZee.append(ele1tightIdC.get2ToId() + ele2tightIdC.get2ToId());
   tautaubkgInZee.append(eeInvMassC.get2ToId());
   tautaubkgInZee.append(elemet200C.get2ToId());
   tautaubkgInZee.append(njetsEmanC.get2ToId());
   tautaubkgInZee.append(jet1ptC.get2ToId());
   tautaubkgInZee.append(jjdphiEmanC.get2ToId());
   tautaubkgInZee.append(muLooseVetoC.get2ToId());
   tautaubkgInZee.append(gammaLooseVetoC.get2ToId());

   mask *ele_acc_eff[nMetBins];
   for ( Int_t i = 0; i < nMetBins; i++) {
     ele_acc_eff[i] = new mask;
     ele_acc_eff[i]->setName(Form("ele_acc_eff:  %3.0lf < met < %3.0lf",metBinEdges[i], metBinEdges[i+1]));
     ele_acc_eff[i]->append(maskEleGen);
     ele_acc_eff[i]->append(maskMonoJetSelection);
     ele_acc_eff[i]->append(maskElectronAcceptance);
     ele_acc_eff[i]->append(maskElectronEfficiency);
   }

   Double_t nTotalWeightedEvents = 0.0;    // total events (including weights)

   TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 
   TH1::StatOverflows();                 //enable use of underflows and overflows for statistics computation 
   TVirtualFitter::SetDefaultFitter("Minuit");

   Int_t Hcolor[] = {1,2,3,4,5,6,7,8,9,12,18,30,38,41,42,46,47,49};    

   TH1D *HeeInvMass[nMetBins];
   TH1D *HzeejetsInvMassMetBinGenEle[nMetBins];
   TH1D *HzeejetsInvMassMetBinGenTau[nMetBins];
   for (Int_t i = 0; i < nMetBins; i++) {
     HeeInvMass[i] = new TH1D(Form("HeeInvMass[%i]",i),"",30,60.,120.);
     HzeejetsInvMassMetBinGenEle[i] = new TH1D(Form("HzeejetsInvMassMetBinGenEle[%i]",i),"",30,60.,120.);
     HzeejetsInvMassMetBinGenTau[i] = new TH1D(Form("HzeejetsInvMassMetBinGenTau[%i]",i),"",30,60.,120.);
   } 

   TH1D *HzeejetsYieldsMetBin = new TH1D("HzeejetsYieldsMetBin","yields of Zee control sample in bins of met;#slash{E}_{T};# of events",
					   nMetBins,metBinEdges);
   TH1D *HzeejetsYieldsMetBinGenEle = new TH1D("HzeejetsYieldsMetBinGenEle","yields of Zee control sample (Z->ee) in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdges);
   TH1D *HzeejetsYieldsMetBinGenTau = new TH1D("HzeejetsYieldsMetBinGenTau","yields of Zee control sample (Z->#tau#tau) in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdges);

   TH1D *HZtoEleEleRecoPt = new TH1D("HZtoEleEleRecoPt","",101,0.,1010);
   TH1D *HZtoEleEleGenPt = new TH1D("HZtoEleEleGenPt","",101,0.,1010);
   // this is the histogram with reco/gen
   TH1D *HZtoEleElePt_RecoGenRatio = new TH1D("HZtoEleElePt_RecoGenRatio","",101,0.,1010.);
   // histogram of reco/gen distribution function
   TH1D *HZtoEleElePt_RecoGenRatio_pdf = new TH1D("HZtoEleElePt_RecoGenRatio_pdf","",100,0.,2.);
   TH1D *HZtoEleEleRecoPt_MetBin[nMetBins];
   TH1D *HZtoEleEleGenPt_MetBin[nMetBins];
   TH1D *HZtoEleElePt_RecoGenRatio_MetBin[nMetBins];
   TH1D *HZtoEleElePt_RecoGenRatio_pdf_MetBin[nMetBins];
   for (Int_t i = 0; i < nMetBins; i++) {
     HZtoEleEleRecoPt_MetBin[i] = new TH1D(Form("HZtoEleEleRecoPt_MetBin[%i]",i),"",101,0.,1010.);
     HZtoEleEleGenPt_MetBin[i] = new TH1D(Form("HZtoEleEleGenPt_MetBin[%i]",i),"",101,0.,1010.);
     HZtoEleElePt_RecoGenRatio_MetBin[i] = new TH1D(Form("HZtoEleElePt_RecoGenRatio_MetBin[%i]",i),"",101,0.,1010.);
     HZtoEleElePt_RecoGenRatio_pdf_MetBin[i] = new TH1D(Form("HZtoEleElePt_RecoGenRatio_pdf_MetBin[%i]",i),"",100,0.,2.);
   } 

   TVector2 e1, e2, met, ele, eleVectorSum;  // ele is any electron (mu1 is the first and mu2 is the second valid in the list sorted by pt)
   // following indices refer to the leading pair of OS/SF in the list of LepGood. They are initialized with 0 and 1 by default, but can be set with function
   // myGetPairIndexInArray (see functionsForAnalysis.cc for reference). 
   // When myGetPairIndexInArray() is called, the index of "correct" particles will be used. If they are not found (e.g. a pair of OS/SF is mismeasured as 2 e+), 
   // indices are set as 0 and 1 (and following selection asking lep[0] and lep[1] to be OS or whatever will fail).
   Int_t firstIndex = 0;
   Int_t secondIndex = 1;

   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"zeejetsAna::loop()"<<endl;
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
	     (LepGood_pt[secondIndex] > 30) && (LepGood_pt[secondIndex] > 18) ) )  continue;

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

     Double_t metNoEle = (eleVectorSum + met).Mod();

     eventMask += njetsC.addToMask(nJet30a <= NJETS);
     eventMask += njetsEmanC.addToMask( (nJet30 == 1 || nJet30 == 2) && jetclean > 0.5);
     eventMask += jjdphiEmanC.addToMask( nJet30 == 1 || (nJet == 2 && abs(dphijj) < J1J2DPHI));
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
     eventMask += genElectronsC.addToMask( myPartGenAlgo(nGenPart, GenPart_pdgId, GenPart_motherId, ELE_PDGID, 23) );
     eventMask += genTausC.addToMask( myPartGenAlgo(nGenPart, GenPart_pdgId, GenPart_motherId, 15, 23) );  

     ele_Acc_Eff.countEvents(eventMask,newwgt);
     zeejetsControlSample.countEvents(eventMask,newwgt);
     zeejetsControlSampleGenEle.countEvents(eventMask,newwgt);
     tautaubkgInZee.countEvents(eventMask, newwgt);

     Double_t ZtoEleEleRecoPt = (e1 + e2).Mod();
     // since when I use the following index I will ask for 2 ele from Z, it's enough to use directly the genZ instead of 2 genEle
     // thus I look for a Z whose pdgId is 23
     Int_t Z_index = myGetPartIndex(23,nGenPart,GenPart_pdgId);
     Double_t ZtoEleEleGenPt = GenPart_pt[Z_index];

     // filling histogram with yields at the end of the selection in bins of met
     if ( ((eventMask & zeejetsControlSample.globalMask.back()) == zeejetsControlSample.globalMask.back()) ) {  
       // this histogram holds the final yields in bins of MET
       HzeejetsYieldsMetBin->Fill(metNoEle,newwgt);     
     }
     if ( ((eventMask & zeejetsControlSampleGenEle.globalMask.back()) == zeejetsControlSampleGenEle.globalMask.back()) ) {  
       // this histogram holds the final yields in bins of MET
       HzeejetsYieldsMetBinGenEle->Fill(metNoEle,newwgt);
       HZtoEleEleRecoPt->Fill(ZtoEleEleRecoPt,newwgt);
       HZtoEleEleGenPt->Fill(ZtoEleEleGenPt,newwgt);
       if (ZtoEleEleGenPt != 0) HZtoEleElePt_RecoGenRatio_pdf->Fill(ZtoEleEleRecoPt/ZtoEleEleGenPt,newwgt);
     }
     if ( ((eventMask & tautaubkgInZee.globalMask.back()) == tautaubkgInZee.globalMask.back()) ) {  
       // this histogram holds the final yields in bins of MET
       HzeejetsYieldsMetBinGenTau->Fill(metNoEle,newwgt);  
     }

     if ((metNoEle> metBinEdges[0]) && (metNoEle < metBinEdges[nMetBins])) {

       Int_t bin = myGetBin(metNoEle,metBinEdges,nMetBins);
       ele_acc_eff[bin]->countEvents(eventMask,newwgt);
       if (((eventMask & zeejetsControlSample.globalMask.back()) == zeejetsControlSample.globalMask.back())) {
	 // this histogram holds the invariant mass distribution (one for each met bin)
   	 HeeInvMass[bin]->Fill(mZ1,newwgt);	
       }
       if ( ((eventMask & zeejetsControlSampleGenEle.globalMask.back()) == zeejetsControlSampleGenEle.globalMask.back()) ) {  
	 HzeejetsInvMassMetBinGenEle[bin]->Fill(mZ1,newwgt); 
	 HZtoEleEleRecoPt_MetBin[bin]->Fill(ZtoEleEleRecoPt,newwgt);
	 HZtoEleEleGenPt_MetBin[bin]->Fill(ZtoEleEleGenPt,newwgt);
	 if (ZtoEleEleGenPt != 0) HZtoEleElePt_RecoGenRatio_pdf_MetBin[bin]->Fill(ZtoEleEleRecoPt/ZtoEleEleGenPt,newwgt);
       }
       if ( ((eventMask & tautaubkgInZee.globalMask.back()) == tautaubkgInZee.globalMask.back()) ) {  
	 HzeejetsInvMassMetBinGenTau[bin]->Fill(mZ1,newwgt);   
       }
              
     } 

   }

   cout<<endl;   
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &ele_Acc_Eff);
   for (Int_t i = 0; i < nMetBins; i++ ) {
     selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, ele_acc_eff[i] );
   }
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &zeejetsControlSample);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &zeejetsControlSampleGenEle);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &tautaubkgInZee);
   cout<<endl;   

   for (Int_t i = 0; i < nMetBins; i++) {
     cout<<" elemet in ["<<metBinEdges[i]<<" , "<<metBinEdges[i+1]<<"] :     HeeInvMass["<<i<<"]->GetSumOfWeights() = ";
     cout<<HeeInvMass[i]->GetSumOfWeights()<<endl;
   }

   cout << "Printing acceptance and efficiency." << endl;
   TH1D* Hacc = new TH1D("Hacc","",nMetBins,metBinEdges);
   TH1D* Heff = new TH1D("Heff","",nMetBins,metBinEdges);
   for (Int_t i = 0; i < nMetBins; i++) {
     // do not merge with previous loop: I want to print them after the previous loop because I will copy and paste this output to make acc & eff table
     Double_t acc = ele_acc_eff[i]->nEvents[2]/ele_acc_eff[i]->nEvents[1];
     Double_t eff = ele_acc_eff[i]->nEvents[3]/ele_acc_eff[i]->nEvents[2];
     Double_t accStatErr = sqrt(acc * (1 - acc) / ele_acc_eff[i]->nEvents[1]);
     Double_t effStatErr = sqrt(eff * (1 - eff) / ele_acc_eff[i]->nEvents[2]);
     Hacc->SetBinContent(i+1,acc);
     Hacc->SetBinError(i+1,accStatErr);
     Heff->SetBinContent(i+1,eff);
     Heff->SetBinError(i+1,effStatErr);
     cout<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" ";
     cout<<acc<<" "<<accStatErr<<" ";
     cout<<eff<<" "<<effStatErr<<endl;
   }
   Hacc->GetXaxis()->SetTitle("#slash{E}_{T} [GeV]");
   Hacc->GetYaxis()->SetTitle("Acceptance");
   Hacc->GetYaxis()->CenterTitle();
   Heff->GetXaxis()->SetTitle("#slash{E}_{T} [GeV]");
   Heff->GetYaxis()->SetTitle("efficiency");
   Heff->GetYaxis()->CenterTitle();

   TH1D* Hacceff = new TH1D("Hacceff","",nMetBins,metBinEdges);
   cout << "Printing acceptance * efficiency" << endl;
   for (Int_t i = 0; i < nMetBins; i++) {
     // do not merge with previous loop: I want to print them after the previous loop because I will copy and paste this output to make acc * eff table
     Double_t acceff = ele_acc_eff[i]->nEvents[3]/ele_acc_eff[i]->nEvents[1];
     Double_t acceffStatErr = sqrt(acceff * (1 - acceff) / ele_acc_eff[i]->nEvents[1]);
     Hacceff->SetBinContent(i+1,acceff);
     Hacceff->SetBinError(i+1,acceffStatErr);
     cout<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" ";
     cout<<acceff<<" "<<acceffStatErr<<endl;;
   }
   Hacceff->GetXaxis()->SetTitle("#slash{E}_{T} [GeV]");
   Hacceff->GetYaxis()->SetTitle("A #times #epsilon");
   Hacceff->GetYaxis()->CenterTitle();
   Hacceff->SaveAs("hist_ee_AccEff_metBin.root");

   const char* accAndEff_filename = "hist_ee_accEff_all_metBin.root";
   TFile *accAndEff_file = new TFile(accAndEff_filename,"RECREATE");
   if (!accAndEff_file->IsOpen()) {
     cout<<"Error: file \""<<accAndEff_filename<<"\" was not opened."<<endl;
   } else {
     for (Int_t i = 0; i < nMetBins; i++) {
       Hacc->Write();
       Heff->Write();
       Hacceff->Write();
     }
     accAndEff_file->Close();
   }
   delete accAndEff_file;

   char answer = '\0';
   string fName = "zeejetsAnaYields.txt";
   const char* filename = fName.c_str();   //this file is like a logbook

   answer = myAskToSaveFile(filename);

   if (answer == 'y') {

     ofstream myfile(filename,ios::app);

     if ( !myfile.is_open() ) {
       cout<<"Error: unable to open file "<<filename<<" !"<<endl;

     } else {
       //when writing on file, the date is printed as well unless an error occurs
       string command = "date>>" + fName;
       //when writing on file, the date is printed as well unless an error occurs
       if ( system(command.c_str()) != 0) cout<<"Error during \"system("<<command<<")\" call"<<endl; 
       myfile<<endl;
       myfile<<endl;
       selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &ele_Acc_Eff);      
       selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &zeejetsControlSample);
       selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &zeejetsControlSampleGenEle);
       selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &tautaubkgInZee);
       myfile<<endl;
       myfile<<endl;
       myfile.close();

     }
 
   }

   const char *texfname = "zeejetsAnaYields.tex";
   
   answer = myAskToSaveFile(texfname);

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
       makeTableTex(fp, LUMI, nTotalWeightedEvents, &ele_Acc_Eff, commentInTable);      
       commentInTable = "Note that cuts on second jet are applied only if a second jet exists with $p_t$ > 30\\,GeV.";
       makeTableTex(fp, LUMI, nTotalWeightedEvents, &zeejetsControlSample,commentInTable);
       makeTableTex(fp, LUMI, nTotalWeightedEvents, &zeejetsControlSampleGenEle,commentInTable);
       makeTableTex(fp, LUMI, nTotalWeightedEvents, &tautaubkgInZee,commentInTable);
       fprintf(fp,"\\end{document}\n");      
       fclose(fp);
     }

   }

   const char *histFileName = "histZeejetsAnaInvMass.root";
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

   const char *histFileName2 = "histZeejetsInvMassCS_ZJetsToee.root";
   cout<<"Saving histograms in file \""<<histFileName2<<"\" ..."<<endl;
   TFile *histFile2 = new TFile(histFileName2,"RECREATE");

   if (!histFile2->IsOpen()) {

     cout<<"Error: file \""<<histFileName2<<"\" was not opened."<<endl;

   } else {

     for (Int_t i = 0; i < nMetBins; i++) {
       HzeejetsInvMassMetBinGenEle[i]->Write();
     }

     histFile2->Close();
     
   }

   delete histFile2;

   const char *histFileName3 = "histZeejetsInvMassCS_ZJetsToTauTau.root";
   cout<<"Saving histograms in file \""<<histFileName3<<"\" ..."<<endl;
   TFile *histFile3 = new TFile(histFileName3,"RECREATE");

   if (!histFile3->IsOpen()) {

     cout<<"Error: file \""<<histFileName3<<"\" was not opened."<<endl;

   } else {

     for (Int_t i = 0; i < nMetBins; i++) {
       HzeejetsInvMassMetBinGenTau[i]->Write();
     }

     histFile3->Close();
     
   }

   delete histFile3;

   const char *YieldsFileName = "histZeejetsAnaYieldsMetBin.root";
   cout<<"Saving histogram \""<<HzeejetsYieldsMetBin->GetName()<<"\" in file \""<<YieldsFileName<<"\" ..."<<endl;
   TFile *YieldsFile = new TFile(YieldsFileName,"RECREATE");

   if (!YieldsFile->IsOpen()) cout<<"Error: file \""<<YieldsFileName<<"\" was not opened."<<endl;
   else HzeejetsYieldsMetBin->Write();

   myPrintYieldsMetBin(HzeejetsYieldsMetBin,metBinEdges,nMetBins);
   
   YieldsFile->Close();
     
   delete YieldsFile;

   HzeejetsYieldsMetBinGenEle->SaveAs("histZeejetsAnaYieldsMetBinGenEle.root");
   myPrintYieldsMetBin(HzeejetsYieldsMetBinGenEle,metBinEdges,nMetBins);
   myPrintYieldsMetBin(HzeejetsYieldsMetBinGenEle,metBinEdges,nMetBins,"histZeejetsAnaYieldsMetBinGenEle.dat");

   HzeejetsYieldsMetBinGenTau->SaveAs("histZeejetsAnaYieldsMetBinGenTau.root");

   // I add overflow bin's content in the last bin for all histograms where that is needed
   // for those histogram filled with Divide() method, it's not done as long as it was already done on the histograms given as
   // argument to the Divide() method
   myAddOverflowInLastBin(HZtoEleEleRecoPt);
   myAddOverflowInLastBin(HZtoEleEleGenPt);

   // saving results about PtZReco/PtZGen
   HZtoEleEleRecoPt->SaveAs("histZeeRecoPt.root");
   HZtoEleEleGenPt->SaveAs("histZeeGenPt.root");
   HZtoEleElePt_RecoGenRatio->Divide(HZtoEleEleRecoPt,HZtoEleEleGenPt);
   HZtoEleElePt_RecoGenRatio->SaveAs("histZeePt_RecoGenRatio.root");
   HZtoEleElePt_RecoGenRatio_pdf->SaveAs("histZeePt_RecoGenRatio_pdf.root");

   // save binned histograms in another file
   const char *histFileName4 = "histZeePt_RecoGenRatio_MetBin.root";
   cout<<"Saving histograms in file \""<<histFileName4<<"\" ..."<<endl;
   TFile *histFile4 = new TFile(histFileName4,"RECREATE");

   if (!histFile4->IsOpen()) {

     cout<<"Error: file \""<<histFileName4<<"\" was not opened."<<endl;

   } else {

     for (Int_t i = 0; i < nMetBins; i++) {
       myAddOverflowInLastBin(HZtoEleEleRecoPt_MetBin[i]);
       myAddOverflowInLastBin(HZtoEleEleGenPt_MetBin[i]);
       HZtoEleElePt_RecoGenRatio_MetBin[i]->Divide(HZtoEleEleRecoPt_MetBin[i],HZtoEleEleGenPt_MetBin[i]);
       HZtoEleElePt_RecoGenRatio_MetBin[i]->Write();
     }

     histFile4->Close();
     
   }

   delete histFile4;

   const char *histFileName5 = "histZeePt_RecoGenRatio_pdf_MetBin.root";
   cout<<"Saving histograms in file \""<<histFileName5<<"\" ..."<<endl;
   TFile *histFile5 = new TFile(histFileName5,"RECREATE");

   if (!histFile5->IsOpen()) {

     cout<<"Error: file \""<<histFileName5<<"\" was not opened."<<endl;

   } else {

     for (Int_t i = 0; i < nMetBins; i++) {
       HZtoEleElePt_RecoGenRatio_pdf_MetBin[i]->Write();
     }

     histFile5->Close();
     
   }

   delete histFile5;

   const char *histFileName6 = "histZeeRecoPt_MetBin.root";
   cout<<"Saving histograms in file \""<<histFileName6<<"\" ..."<<endl;
   TFile *histFile6 = new TFile(histFileName6,"RECREATE");

   if (!histFile6->IsOpen()) {

     cout<<"Error: file \""<<histFileName6<<"\" was not opened."<<endl;

   } else {

     for (Int_t i = 0; i < nMetBins; i++) {
       HZtoEleEleRecoPt_MetBin[i]->Write();
     }

     histFile6->Close();
     
   }

   delete histFile6;

   const char *histFileName7 = "histZeeGenPt_MetBin.root";
   cout<<"Saving histograms in file \""<<histFileName7<<"\" ..."<<endl;
   TFile *histFile7 = new TFile(histFileName7,"RECREATE");

   if (!histFile7->IsOpen()) {

     cout<<"Error: file \""<<histFileName7<<"\" was not opened."<<endl;

   } else {

     for (Int_t i = 0; i < nMetBins; i++) {
       HZtoEleEleGenPt_MetBin[i]->Write();
     }

     histFile7->Close();
     
   }

   delete histFile7;

}
