#define zmumujetsAna_cxx
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

#ifdef zmumujetsAna_cxx

zmumujetsAna::zmumujetsAna(TTree *tree) : edimarcoTree(tree) {
  //cout <<"check in constructor "<<endl;
  Init(tree);

}

#endif

//#define LUMI 5                         // integrated luminosity in fb^-1 to which yields are normalized
#define MU_PDGID 13

void zmumujetsAna::loop()
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
   fChain->SetBranchStatus("genLep_pt",1);
   fChain->SetBranchStatus("genLep_eta",1);
   //fChain->SetBranchStatus("m2l",1);  // m(ll)  (I can compute it myself, maybe it's better)
   fChain->SetBranchStatus("mZ1",1);  // best m(ll) SF/OS

   fChain->SetBranchStatus("nGenPart",1);
   fChain->SetBranchStatus("GenPart_pdgId",1);
   fChain->SetBranchStatus("GenPart_motherId",1);
   fChain->SetBranchStatus("GenPart_pt",1);
   fChain->SetBranchStatus("GenPart_phi",1);
   fChain->SetBranchStatus("GenPart_mass",1);
   fChain->SetBranchStatus("GenPart_motherIndex",1);

   fChain->SetBranchStatus("met_pt",1);
   fChain->SetBranchStatus("met_phi",1);

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
   selection gammaLooseVetoC("gammaLooseVetoC","photons veto");
   // additional selections for control sample
   selection oppChargeLeptonsC("oppChargeLeptonsC","OS/SF leptons");
   selection twomuonsC("twomuonsC","muons");
   selection muLooseVetoC("muLooseVetoC","muons veto");
   selection twomuLooseC("twomuLooseC","2 loose muons");
   selection mu1tightIdC("mu1tightIdC","leading muon tight","tight ID + relIso04 (as Emanuele)");
   selection twomuTightC("twomuTightC","2 tight muons");
   selection mu1ptC("mu1ptC","mu1pt > 20","leading pt muon");
   selection mu2ptC("mu2ptC","mu2pt > 10");
   selection mu1etaC("mu1etaC","|mu1eta| < 2.4","leading pt muon");  
   selection mu2etaC("mu2etaC","|mu2eta| < 2.4");
   selection mumuInvMassC("mumuInvMassC","muons mass in [60,120]");
   selection genTausC("genTausC","taus generated");               // 11, 13, 15 for e, mu, tau  ( - sign for antiparticles)              
   selection genMuonsC("genMuonsC","muons generated");            // 11, 13, 15 for e, mu, tau  ( - sign for antiparticles)  
   selection acceptanceC("acceptanceC","acceptance cuts");
   selection efficiencyC("efficiencyC","efficiency cuts");

   selection::checkMaskLength();
   selection::printActiveSelections(cout);

   UInt_t maskMonoJetSelection = njetsEmanC.get2ToId() + jet1ptC.get2ToId() + jjdphiEmanC.get2ToId() +
                                   eLooseVetoC.get2ToId() + gammaLooseVetoC.get2ToId();

   UInt_t maskMuonAcceptance = mu1ptC.get2ToId() + mu2ptC.get2ToId() + mu1etaC.get2ToId() + mu2etaC.get2ToId() +
                                                    mumuInvMassC.get2ToId(); 

   UInt_t maskMuonEfficiency = twomuLooseC.get2ToId() + mu1tightIdC.get2ToId(); 

   UInt_t maskMuAxe = mu1ptC.get2ToId() + mu2ptC.get2ToId() + mu1etaC.get2ToId() + mu2etaC.get2ToId() +
                                            mumuInvMassC.get2ToId() + twomuLooseC.get2ToId() + mu1tightIdC.get2ToId();

   UInt_t maskMuGen = genMuonsC.get2ToId();
   UInt_t maskTauGen = genTausC.get2ToId();    

   mask mu_Acc_Eff("mumu acceptance and efficiency"); 
   mu_Acc_Eff.append(maskMuGen);
   mu_Acc_Eff.append(maskMonoJetSelection);
   mu_Acc_Eff.append(acceptanceC.get2ToId());
   mu_Acc_Eff.append(efficiencyC.get2ToId());

   mask zmumujetsControlSample("Z->mumu control sample with selection flow as Emanuele's");
   zmumujetsControlSample.append(twomuLooseC.get2ToId() + oppChargeLeptonsC.get2ToId());
   zmumujetsControlSample.append(twomuonsC.get2ToId());
   zmumujetsControlSample.append(mu1tightIdC.get2ToId());
   zmumujetsControlSample.append(mumuInvMassC.get2ToId());
   zmumujetsControlSample.append(njetsEmanC.get2ToId());
   zmumujetsControlSample.append(jet1ptC.get2ToId());
   zmumujetsControlSample.append(jjdphiEmanC.get2ToId());
   zmumujetsControlSample.append(eLooseVetoC.get2ToId());
   zmumujetsControlSample.append(gammaLooseVetoC.get2ToId());

   mask zmumujetsControlSampleGenMu("Z->mumu control sample (mu generated) with selection flow as Emanuele's");
   zmumujetsControlSampleGenMu.append(maskMuGen);
   zmumujetsControlSampleGenMu.append(twomuLooseC.get2ToId() + oppChargeLeptonsC.get2ToId());
   zmumujetsControlSampleGenMu.append(twomuonsC.get2ToId());
   zmumujetsControlSampleGenMu.append(mu1tightIdC.get2ToId());
   zmumujetsControlSampleGenMu.append(mumuInvMassC.get2ToId());
   zmumujetsControlSampleGenMu.append(njetsEmanC.get2ToId());
   zmumujetsControlSampleGenMu.append(jet1ptC.get2ToId());
   zmumujetsControlSampleGenMu.append(jjdphiEmanC.get2ToId());
   zmumujetsControlSampleGenMu.append(eLooseVetoC.get2ToId());
   zmumujetsControlSampleGenMu.append(gammaLooseVetoC.get2ToId());

   mask tautaubkgInZmumu("tau tau background in Z->mumu control sample");
   tautaubkgInZmumu.append(maskTauGen);
   tautaubkgInZmumu.append(twomuLooseC.get2ToId() + oppChargeLeptonsC.get2ToId());
   tautaubkgInZmumu.append(twomuonsC.get2ToId());
   tautaubkgInZmumu.append(mu1tightIdC.get2ToId());
   tautaubkgInZmumu.append(mumuInvMassC.get2ToId());
   tautaubkgInZmumu.append(njetsEmanC.get2ToId());
   tautaubkgInZmumu.append(jet1ptC.get2ToId());
   tautaubkgInZmumu.append(jjdphiEmanC.get2ToId());
   tautaubkgInZmumu.append(eLooseVetoC.get2ToId());
   tautaubkgInZmumu.append(gammaLooseVetoC.get2ToId());

   mask *mu_acc_eff[nMetBins];
   for ( Int_t i = 0; i < nMetBins; i++) {
     mu_acc_eff[i] = new mask;
     mu_acc_eff[i]->setName(Form("mu_acc_eff:  %3.0lf < met < %3.0lf",metBinEdges[i], metBinEdges[i+1]));
     mu_acc_eff[i]->append(maskMuGen);
     mu_acc_eff[i]->append(maskMonoJetSelection);
     mu_acc_eff[i]->append(maskMuonAcceptance);
     mu_acc_eff[i]->append(maskMuonEfficiency);
   }


   Double_t nTotalWeightedEvents = 0.0;    // total events (including weights)

   TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 
   TH1::StatOverflows();                 //enable use of underflows and overflows for statistics computation 
   TVirtualFitter::SetDefaultFitter("Minuit");

   Int_t Hcolor[] = {1,2,3,4,5,6,7,8,9,12,18,30,38,41,42,46,47,49};    

   TH1D *HmumuInvMass[nMetBins];
   TH1D *HzmumujetsInvMassMetBinGenMu[nMetBins];
   TH1D *HzmumujetsInvMassMetBinGenTau[nMetBins];
   for (Int_t i = 0; i < nMetBins; i++) {
     HmumuInvMass[i] = new TH1D(Form("HmumuInvMass[%i]",i),"",30,60.,120.);
     HzmumujetsInvMassMetBinGenMu[i] = new TH1D(Form("HzmumujetsInvMassMetBinGenMu[%i]",i),"",30,60.,120.);
     HzmumujetsInvMassMetBinGenTau[i] = new TH1D(Form("HzmumujetsInvMassMetBinGenTau[%i]",i),"",30,60.,120.);
   } 

   TH1D *HzmumujetsYieldsMetBin = new TH1D("HzmumujetsYieldsMetBin","yields of Zmumu control sample in bins of met;#slash{E}_{T};# of events",
					   nMetBins,metBinEdges);
   TH1D *HzmumujetsYieldsMetBinGenMu = new TH1D("HzmumujetsYieldsMetBinGenMu","yields of Zmumu control sample (Z->#mu#mu) in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdges);
   TH1D *HzmumujetsYieldsMetBinGenTau = new TH1D("HzmumujetsYieldsMetBinGenTau","yields of Zmumu control sample (Z->#tau#tau) in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdges);

   TH1D *HZtoMuMuRecoPt = new TH1D("HZtoMuMuRecoPt","",101,0.,1010);
   TH1D *HZtoMuMuGenPt = new TH1D("HZtoMuMuGenPt","",101,0.,1010);
   // this is the histogram with reco/gen
   TH1D *HZtoMuMuPt_RecoGenRatio = new TH1D("HZtoMuMuPt_RecoGenRatio","",101,0.,1010.);
   // histogram of reco/gen distribution function
   TH1D *HZtoMuMuPt_RecoGenRatio_pdf = new TH1D("HZtoMuMuPt_RecoGenRatio_pdf","",100,0.,2.);
   TH1D *HZtoMuMuRecoPt_MetBin[nMetBins];
   TH1D *HZtoMuMuGenPt_MetBin[nMetBins];
   TH1D *HZtoMuMuPt_RecoGenRatio_MetBin[nMetBins];
   TH1D *HZtoMuMuPt_RecoGenRatio_pdf_MetBin[nMetBins];
   for (Int_t i = 0; i < nMetBins; i++) {
     HZtoMuMuRecoPt_MetBin[i] = new TH1D(Form("HZtoMuMuRecoPt_MetBin[%i]",i),"",101,0.,1010.);
     HZtoMuMuGenPt_MetBin[i] = new TH1D(Form("HZtoMuMuGenPt_MetBin[%i]",i),"",101,0.,1010.);
     HZtoMuMuPt_RecoGenRatio_MetBin[i] = new TH1D(Form("HZtoMuMuPt_RecoGenRatio_MetBin[%i]",i),"",101,0.,1010.);
     HZtoMuMuPt_RecoGenRatio_pdf_MetBin[i] = new TH1D(Form("HZtoMuMuPt_RecoGenRatio_pdf_MetBin[%i]",i),"",100,0.,2.);
   } 

   TVector2 mu1, mu2;

   // following indices refer to the leading pair of OS/SF in the list of LepGood. They are initialized with 0 and 1 
   //by default, but can be set with function
   // myGetPairIndexInArray (see functionsForAnalysis.cc for reference). 
   // When myGetPairIndexInArray() is called, the index of "correct" particles will be used. If they are not found 
   //(e.g. a pair of OS/SF is mismeasured as 2 mu+), 
   // indices are set as 0 and 1 (and following selection asking lep[0] and lep[1] to be OS or whatever will fail).
   Int_t firstIndex = 0;
   Int_t secondIndex = 1;

   Int_t firstIndexGen = 0;
   Int_t secondIndexGen = 1;
   Int_t genMuFound = 0;
   Int_t Z_index = 0; 

   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"zmumujetsAna::loop()"<<endl;
   cout<<"nentries = "<<nentries<<endl;   

   Long64_t nbytes = 0, nb = 0;
   for (Int_t jentry=0; jentry<nentries; jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;   

     UInt_t eventMask = 0; 
     Double_t newwgt = weight * LUMI;

     genMuFound = myPartGenAlgo(nGenPart, GenPart_pdgId, GenPart_motherId, MU_PDGID, 23, firstIndexGen, secondIndexGen);
     Z_index = GenPart_motherIndex[firstIndexGen];
     // I find the indices corresponding to the 2 leading lepton
     //cout<<"entry : "<<jentry<<endl;
     myGetPairIndexInArray(MU_PDGID, nLepGood, LepGood_pdgId, firstIndex, secondIndex);

     mu1.SetMagPhi(LepGood_pt[firstIndex],LepGood_phi[firstIndex]);
     mu2.SetMagPhi(LepGood_pt[secondIndex],LepGood_phi[secondIndex]);

     nTotalWeightedEvents += newwgt;  // counting events with weights

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
     //eventMask += ntausC.addToMask(ntaus);
     eventMask += gammaLooseVetoC.addToMask(nGamma15V == 0);
     //eventMask += mumet200C.addToMask(mumet);
     for (Int_t i = 0; i <  metCut.size(); i++) {
       eventMask += mumetC[i].addToMask(metNoMu_pt > metCut[i]);
     }
     eventMask += oppChargeLeptonsC.addToMask( (LepGood_pdgId[firstIndex] + LepGood_pdgId[secondIndex]) == 0);
     eventMask += twomuonsC.addToMask((fabs(LepGood_pdgId[firstIndex]) == MU_PDGID) && (fabs(LepGood_pdgId[secondIndex]) == MU_PDGID));
     eventMask += twomuLooseC.addToMask(nMu10V == 2);
     eventMask += muLooseVetoC.addToMask(nMu10V == 0);
     eventMask += mu1tightIdC.addToMask((LepGood_tightId[firstIndex] == 1) && (LepGood_relIso04[firstIndex] < 0.12 ) && 
					(fabs(LepGood_pdgId[firstIndex]) == MU_PDGID));
     eventMask += mu1ptC.addToMask((LepGood_pt[firstIndex] > 20) && (fabs(LepGood_pdgId[firstIndex]) == MU_PDGID)); 
     eventMask += mu1etaC.addToMask( (fabs(LepGood_eta[firstIndex]) < 2.4) && (fabs(LepGood_pdgId[firstIndex]) == MU_PDGID) );
     eventMask += twomuTightC.addToMask(nMu20T == 2);
     eventMask += mu2ptC.addToMask((LepGood_pt[secondIndex] > 10) && (fabs(LepGood_pdgId[secondIndex]) == MU_PDGID));
     eventMask += mu2etaC.addToMask((fabs(LepGood_eta[secondIndex]) < 2.4) && (fabs(LepGood_pdgId[secondIndex]) == MU_PDGID));
     eventMask += mumuInvMassC.addToMask((mZ1 > 60) && (mZ1 < 120));
     eventMask += genMuonsC.addToMask(genMuFound );
     eventMask += genTausC.addToMask( myPartGenAlgo(nGenPart, GenPart_pdgId, GenPart_motherId, 15, 23) );  
     eventMask += acceptanceC.addToMask( genMuFound && (GenPart_pt[firstIndexGen] > 10) && (GenPart_pt[secondIndexGen] > 10) &&
					 (fabs(GenPart_eta[firstIndexGen]) < 2.4) && (fabs(GenPart_eta[secondIndexGen]) < 2.4) &&
					(GenPart_mass[Z_index] > 60) && (GenPart_mass[Z_index] < 120) );
     eventMask += efficiencyC.addToMask( (nMu10V == 2) && (LepGood_tightId[firstIndex] == 1) && 
							   (LepGood_relIso04[firstIndex] < 0.12 ) &&
							   (fabs(LepGood_pdgId[firstIndex]) == MU_PDGID) );

     Double_t ZtoMuMuRecoPt = (mu1 + mu2).Mod();
     //Z_index = myGetPartIndex(23,nGenPart,GenPart_pdgId);
     // since when I use the following index I will ask for 2 mu from Z, it's enough to use directly the genZ instead of 2 genMu
     // thus I look for a Z whose pdgId is 23
     Double_t ZtoMuMuGenPt = GenPart_pt[Z_index];

     mu_Acc_Eff.countEvents(eventMask,newwgt);
     zmumujetsControlSample.countEvents(eventMask,newwgt);
     zmumujetsControlSampleGenMu.countEvents(eventMask,newwgt);
     tautaubkgInZmumu.countEvents(eventMask, newwgt);

     // filling histogram with yields and invariant mass at the end of the selection in bins of met
     if ( ((eventMask & zmumujetsControlSample.globalMask.back()) == zmumujetsControlSample.globalMask.back()) ) {  
       // this histogram holds the final yields in bins of MET
       HzmumujetsYieldsMetBin->Fill(metNoMu_pt,newwgt);     
     }
     if ( ((eventMask & zmumujetsControlSampleGenMu.globalMask.back()) == zmumujetsControlSampleGenMu.globalMask.back()) ) {  
       // this histogram holds the final yields in bins of MET
       HzmumujetsYieldsMetBinGenMu->Fill(metNoMu_pt,newwgt);
       HZtoMuMuRecoPt->Fill(ZtoMuMuRecoPt,newwgt);
       HZtoMuMuGenPt->Fill(ZtoMuMuGenPt,newwgt);
       if (ZtoMuMuGenPt != 0) HZtoMuMuPt_RecoGenRatio_pdf->Fill(ZtoMuMuRecoPt/ZtoMuMuGenPt,newwgt);
     }
     if ( ((eventMask & tautaubkgInZmumu.globalMask.back()) == tautaubkgInZmumu.globalMask.back()) ) {  
       // this histogram holds the final yields in bins of MET
       HzmumujetsYieldsMetBinGenTau->Fill(metNoMu_pt,newwgt);  
     }

     if ((metNoMu_pt > metBinEdges[0]) && (metNoMu_pt < metBinEdges[nMetBins])) {

       Int_t bin = myGetBin(metNoMu_pt,metBinEdges,nMetBins);
       mu_acc_eff[bin]->countEvents(eventMask,newwgt);
       if ((eventMask & zmumujetsControlSample.globalMask.back()) == zmumujetsControlSample.globalMask.back()) {
	 // this histogram holds the invariant mass distribution (one for each met bin)
   	 HmumuInvMass[bin]->Fill(mZ1,newwgt);   
       }
       if ( ((eventMask & zmumujetsControlSampleGenMu.globalMask.back()) == zmumujetsControlSampleGenMu.globalMask.back()) ) {  
	 HzmumujetsInvMassMetBinGenMu[bin]->Fill(mZ1,newwgt); 
	 HZtoMuMuRecoPt_MetBin[bin]->Fill(ZtoMuMuRecoPt,newwgt);
	 HZtoMuMuGenPt_MetBin[bin]->Fill(ZtoMuMuGenPt,newwgt);
	 if (ZtoMuMuGenPt != 0) HZtoMuMuPt_RecoGenRatio_pdf_MetBin[bin]->Fill(ZtoMuMuRecoPt/ZtoMuMuGenPt,newwgt);
       }
       if ( ((eventMask & tautaubkgInZmumu.globalMask.back()) == tautaubkgInZmumu.globalMask.back()) ) {  
       	 HzmumujetsInvMassMetBinGenTau[bin]->Fill(mZ1,newwgt);   
       }
     
     } 

   }

   cout<<endl;   
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &mu_Acc_Eff);
   for (Int_t i = 0; i < nMetBins; i++ ) {
     selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, mu_acc_eff[i] );
   }
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &zmumujetsControlSample);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &zmumujetsControlSampleGenMu);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &tautaubkgInZmumu);
   cout<<endl;   

   for (Int_t i = 0; i < nMetBins; i++) {
     cout<<" mumet in ["<<metBinEdges[i]<<" , "<<metBinEdges[i+1]<<"] :     HmumuInvMass["<<i<<"]->GetSumOfWeights() = ";
     cout<<HmumuInvMass[i]->GetSumOfWeights()<<endl;
   }

   cout << "Printing acceptance and efficiency." << endl;
   TH1D* Hacc = new TH1D("Hacc","",nMetBins,metBinEdges);
   TH1D* Heff = new TH1D("Heff","",nMetBins,metBinEdges);
   for (Int_t i = 0; i < nMetBins; i++) {
     // do not merge with previous loop: I want to print them after the previous loop because I will copy and paste this output to make acc & eff table
     Double_t acc = mu_acc_eff[i]->nEvents[2]/mu_acc_eff[i]->nEvents[1];
     Double_t eff = mu_acc_eff[i]->nEvents[3]/mu_acc_eff[i]->nEvents[2];
     Double_t accStatErr = sqrt(acc * (1 - acc) / mu_acc_eff[i]->nEvents[1]);
     Double_t effStatErr = sqrt(eff * (1 - eff) / mu_acc_eff[i]->nEvents[2]);
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
     Double_t acceff = mu_acc_eff[i]->nEvents[3]/mu_acc_eff[i]->nEvents[1];
     Double_t acceffStatErr = sqrt(acceff * (1 - acceff) / mu_acc_eff[i]->nEvents[1]);
     Hacceff->SetBinContent(i+1,acceff);
     Hacceff->SetBinError(i+1,acceffStatErr);
     cout<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" ";
     cout<<acceff<<" "<<acceffStatErr<<endl;;
   }
   Hacceff->GetXaxis()->SetTitle("#slash{E}_{T} [GeV]");
   Hacceff->GetYaxis()->SetTitle("A #times #epsilon");
   Hacceff->GetYaxis()->CenterTitle();
   Hacceff->SaveAs("hist_mumu_AccEff_metBin.root");

   TH1D *H_BR_ratio = new TH1D("H_BR_ratio","BR(Z->#nu#nu/BR(Z->#mu#mu)",nMetBins,metBinEdges);

   for(Int_t i = 0; i <= nMetBins; i++) {
     H_BR_ratio->SetBinContent(i,RATIO_BR_ZINV_ZLL);
     H_BR_ratio->SetBinError(i,UNC_RATIO_BR_ZINV_ZLL);
   }
   TH1D *HzvvEstimate = new TH1D("HzvvEstimate","yields of Z->#nu#nu estimated as N(Z->#mu#mu) * BR_ratio / (A*#varepsilon)",nMetBins,metBinEdges);
   HzvvEstimate->Multiply(HzmumujetsYieldsMetBinGenMu,H_BR_ratio);
   HzvvEstimate->Divide(Hacceff);
   HzvvEstimate->SaveAs("hist_zvvestimate _ZmumuSample.root");
   myPrintYieldsMetBin(HzvvEstimate,metBinEdges,nMetBins,"histZmumujetsAna_ZnunuYieldsEstimate.dat");
   delete H_BR_ratio; //no need to save it

   const char* accAndEff_filename = "hist_mumu_accEff_all_metBin.root";
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
   string fName = "zmumujetsAnaYields.txt";
   const char* filename = fName.c_str();   //this file is like a logbook

   answer = myAskToSaveFile(filename);

   if (answer == 'y') {
     
     ofstream myfile(filename,ios::app);

     if ( !myfile.is_open() ) {
       cout<<"Error: unable to open file "<<filename<<" !"<<endl;

     } else {
       string command = "date>>" + fName;
       //when writing on file, the date is printed as well unless an error occurs
       if ( system(command.c_str()) != 0) cout<<"Error during \"system("<<command<<")\" call"<<endl;  
       myfile<<endl;
       myfile<<endl;
       selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &mu_Acc_Eff);      
       selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &zmumujetsControlSample);
       selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &zmumujetsControlSampleGenMu);
       selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &tautaubkgInZmumu);
       myfile<<endl;
       myfile<<endl;
       myfile.close();

     }
 
   }

   const char *texfname = "zmumujetsAnaYields.tex";
   
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
       makeTableTex(fp, LUMI, nTotalWeightedEvents, &mu_Acc_Eff, commentInTable);      
       commentInTable = "Note that cuts on second jet are applied only if a second jet exists with $p_t$ > 30\\,GeV.";
       makeTableTex(fp, LUMI, nTotalWeightedEvents, &zmumujetsControlSample,commentInTable);
       makeTableTex(fp, LUMI, nTotalWeightedEvents, &zmumujetsControlSampleGenMu,commentInTable);
       makeTableTex(fp, LUMI, nTotalWeightedEvents, &tautaubkgInZmumu,commentInTable);
       fprintf(fp,"\\end{document}\n");      
       fclose(fp);
     }

   }

   const char *histFileName = "histZmumujetsAnaInvMass.root";
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

   const char *histFileName2 = "histZmumujetsInvMassCS_ZJetsToMuMu.root";
   cout<<"Saving histograms in file \""<<histFileName2<<"\" ..."<<endl;
   TFile *histFile2 = new TFile(histFileName2,"RECREATE");

   if (!histFile2->IsOpen()) {

     cout<<"Error: file \""<<histFileName2<<"\" was not opened."<<endl;

   } else {

     for (Int_t i = 0; i < nMetBins; i++) {
       HzmumujetsInvMassMetBinGenMu[i]->Write();
     }

     histFile2->Close();
     
   }

   delete histFile2;

   const char *histFileName3 = "histZmumujetsInvMassCS_ZJetsToTauTau.root";
   cout<<"Saving histograms in file \""<<histFileName3<<"\" ..."<<endl;
   TFile *histFile3 = new TFile(histFileName3,"RECREATE");

   if (!histFile3->IsOpen()) {

     cout<<"Error: file \""<<histFileName3<<"\" was not opened."<<endl;

   } else {

     for (Int_t i = 0; i < nMetBins; i++) {
       HzmumujetsInvMassMetBinGenTau[i]->Write();
     }

     histFile3->Close();
     
   }

   delete histFile3;				       
				  
   const char *YieldsFileName = "histZmumujetsAnaYieldsMetBin.root";
   cout<<"Saving histogram \""<<HzmumujetsYieldsMetBin->GetName()<<"\" in file \""<<YieldsFileName<<"\" ..."<<endl;
   TFile *YieldsFile = new TFile(YieldsFileName,"RECREATE");

   if (!YieldsFile->IsOpen()) cout<<"Error: file \""<<YieldsFileName<<"\" was not opened."<<endl;
   else HzmumujetsYieldsMetBin->Write();

   myPrintYieldsMetBin(HzmumujetsYieldsMetBin,metBinEdges,nMetBins);

   YieldsFile->Close();
     
   delete YieldsFile;

   HzmumujetsYieldsMetBinGenMu->SaveAs("histZmumujetsAnaYieldsMetBinGenMu.root");
   myPrintYieldsMetBin(HzmumujetsYieldsMetBinGenMu,metBinEdges,nMetBins);
   myPrintYieldsMetBin(HzmumujetsYieldsMetBinGenMu,metBinEdges,nMetBins,"histZmumujetsAnaYieldsMetBinGenMu.dat");

   HzmumujetsYieldsMetBinGenTau->SaveAs("histZmumujetsAnaYieldsMetBinGenTau.root");

   // I add overflow bin's content in the last bin for all histograms where that is needed
   // for those histogram filled with Divide() method, it's not done as long as it was already done on the histograms given as
   // argument to the Divide() method
   myAddOverflowInLastBin(HZtoMuMuRecoPt);
   myAddOverflowInLastBin(HZtoMuMuGenPt);

   // saving results about PtZReco/PtZGen
   HZtoMuMuRecoPt->SaveAs("histZmumuRecoPt.root");
   HZtoMuMuGenPt->SaveAs("histZmumuGenPt.root");
   HZtoMuMuPt_RecoGenRatio->Divide(HZtoMuMuRecoPt,HZtoMuMuGenPt);
   HZtoMuMuPt_RecoGenRatio->SaveAs("histZmumuPt_RecoGenRatio.root");
   HZtoMuMuPt_RecoGenRatio_pdf->SaveAs("histZmumuPt_RecoGenRatio_pdf.root");

   // save binned histograms in another file
   const char *histFileName4 = "histZmumuPt_RecoGenRatio_MetBin.root";
   cout<<"Saving histograms in file \""<<histFileName4<<"\" ..."<<endl;
   TFile *histFile4 = new TFile(histFileName4,"RECREATE");

   if (!histFile4->IsOpen()) {

     cout<<"Error: file \""<<histFileName4<<"\" was not opened."<<endl;

   } else {

     for (Int_t i = 0; i < nMetBins; i++) {
       myAddOverflowInLastBin(HZtoMuMuRecoPt_MetBin[i]);
       myAddOverflowInLastBin(HZtoMuMuGenPt_MetBin[i]);
       HZtoMuMuPt_RecoGenRatio_MetBin[i]->Divide(HZtoMuMuRecoPt_MetBin[i],HZtoMuMuGenPt_MetBin[i]);
       HZtoMuMuPt_RecoGenRatio_MetBin[i]->Write();
     }

     histFile4->Close();
     
   }

   delete histFile4;

   const char *histFileName5 = "histZmumuPt_RecoGenRatio_pdf_MetBin.root";
   cout<<"Saving histograms in file \""<<histFileName5<<"\" ..."<<endl;
   TFile *histFile5 = new TFile(histFileName5,"RECREATE");

   if (!histFile5->IsOpen()) {

     cout<<"Error: file \""<<histFileName5<<"\" was not opened."<<endl;

   } else {

     for (Int_t i = 0; i < nMetBins; i++) {
       HZtoMuMuPt_RecoGenRatio_pdf_MetBin[i]->Write();
     }

     histFile5->Close();
     
   }

   delete histFile5;

   const char *histFileName6 = "histZmumuRecoPt_MetBin.root";
   cout<<"Saving histograms in file \""<<histFileName6<<"\" ..."<<endl;
   TFile *histFile6 = new TFile(histFileName6,"RECREATE");

   if (!histFile6->IsOpen()) {

     cout<<"Error: file \""<<histFileName6<<"\" was not opened."<<endl;

   } else {

     for (Int_t i = 0; i < nMetBins; i++) {
       HZtoMuMuRecoPt_MetBin[i]->Write();
     }

     histFile6->Close();
     
   }

   delete histFile6;

   const char *histFileName7 = "histZmumuGenPt_MetBin.root";
   cout<<"Saving histograms in file \""<<histFileName7<<"\" ..."<<endl;
   TFile *histFile7 = new TFile(histFileName7,"RECREATE");

   if (!histFile7->IsOpen()) {

     cout<<"Error: file \""<<histFileName7<<"\" was not opened."<<endl;

   } else {

     for (Int_t i = 0; i < nMetBins; i++) {
       HZtoMuMuGenPt_MetBin[i]->Write();
     }

     histFile7->Close();
     
   }

   delete histFile7;

}
