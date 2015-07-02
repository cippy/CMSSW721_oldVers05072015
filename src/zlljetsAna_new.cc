#define zlljetsAna_new_cxx
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

#define NVTXS  21                           // # of points for study of u_par and u_perp vs # of reconstructed vertices nvtx
#define FIRST_NVTX  7                    // starting number of vertices for met study   

#ifdef zlljetsAna_new_cxx

zlljetsAna_new::zlljetsAna_new(TTree *tree) : edimarcoTree(tree) {
  //cout <<"check in constructor "<<endl;
  Init(tree);

}

#endif

void zlljetsAna_new::loop(const char* configFileName)
{

   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);  
   // warning: in Emanuele's trees non integer values are float

   fChain->SetBranchStatus("weight",1);   // includes k-factor

   fChain->SetBranchStatus("nMu10V",1);  // # of muons passing loose selection
   fChain->SetBranchStatus("nEle10V",1);  // # of electrons passing loose selection for electron veto
   fChain->SetBranchStatus("nGamma15V",1);  // # of photons passing loose selection for photon veto
   //fChain->SetBranchStatus("nMu20T",1);  // # of muons passing tight selection (pt > 20 + everything else)
   fChain->SetBranchStatus("nTau15V",1);

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
   fChain->SetBranchStatus("GenPart_eta",1);
   fChain->SetBranchStatus("GenPart_phi",1);
   fChain->SetBranchStatus("GenPart_mass",1);
   fChain->SetBranchStatus("GenPart_motherIndex",1);

   fChain->SetBranchStatus("met_pt",1);
   fChain->SetBranchStatus("met_phi",1);

   fChain->SetBranchStatus("metNoMu_pt",1);
   fChain->SetBranchStatus("metNoMu_phi",1);

   fChain->SetBranchStatus("nVert",1);  // number of good vertices

   char ROOT_FNAME[50];
   char TXT_FNAME[50];
   char TEX_FNAME[50];
   char FLAVOUR[10];                   // e.g. "ele", "mu"
   char LL_FLAVOUR[10];             // e.g. "ee", "mumu"
   char CONTROL_SAMPLE[10];   // e.g. "Z-->ee"

   Double_t LUMI;
   Int_t NJETS;
   Double_t J1PT;
   Double_t J1ETA;
   Double_t J2PT;
   Double_t J2ETA;
   Double_t J1J2DPHI;
   Double_t RATIO_BR_ZINV_ZLL;
   Double_t UNC_RATIO_BR_ZINV_ZLL;
   Int_t LEP_PDG_ID;
   Int_t NEG_LEP_PDG_ID2;
   Double_t LEP1PT;
   Double_t LEP2PT;
   Double_t LEP1ETA;
   Double_t LEP2ETA;
   Double_t DILEPMASS_LOW;
   Double_t DILEPMASS_UP;
   Double_t LEP_ISO_04;
   Double_t GENLEP1PT;
   Double_t GENLEP2PT;
   Double_t GENLEP1ETA;
   Double_t GENLEP2ETA;
   Double_t GEN_ZMASS_LOW;
   Double_t GEN_ZMASS_UP;
   Int_t TAU_VETO_FLAG;
   Int_t HLT_FLAG;
   Double_t HLT_LEP1PT;
   Double_t HLT_LEP2PT;
   Double_t HLT_LEP1ETA;
   Double_t HLT_LEP2ETA;

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
     RATIO_BR_ZINV_ZLL = parameterValue[7];
     UNC_RATIO_BR_ZINV_ZLL = parameterValue[8];
     LEP_PDG_ID = (Int_t) parameterValue[9];
     NEG_LEP_PDG_ID2 = (Int_t) parameterValue[10];
     LEP1PT = parameterValue[11];
     LEP2PT = parameterValue[12];
     LEP1ETA = parameterValue[13];
     LEP2ETA = parameterValue[14];
     DILEPMASS_LOW = parameterValue[15];
     DILEPMASS_UP = parameterValue[16];
     LEP_ISO_04 = parameterValue[17];
     GENLEP1PT = parameterValue[18];
     GENLEP2PT = parameterValue[19];
     GENLEP1ETA = parameterValue[20];
     GENLEP2ETA = parameterValue[21];
     GEN_ZMASS_LOW = parameterValue[22];
     GEN_ZMASS_UP = parameterValue[23];
     TAU_VETO_FLAG = (Int_t) parameterValue[24];
     HLT_FLAG = (Int_t) parameterValue[25];
     HLT_LEP1PT = parameterValue[26];
     HLT_LEP2PT = parameterValue[27];
     HLT_LEP1ETA = parameterValue[28];
     HLT_LEP2ETA = parameterValue[29];
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

   selection njetsC("njetsC",Form("njets <= %i",NJETS),"pt > 30; |eta| < 4.7");   // using nJet30a
   selection njetsEmanC("njetsEmanC","njets","1 or 2 jets, cleaning, pt > 30, |eta| < 2.4");       // using nJet30
   selection jet1ptC("jet1ptC",Form("jet1pt > %4.0lf",(Double_t)J1PT));
   selection jet1etaC("jet1etaC",Form("|jet1eta| < %2.1lf",J1ETA));
   selection jet2etaC("jet2etaC",Form("|jet2eta| < %2.1lf",J2ETA),Form("only if njets = %i",NJETS));
   selection jjdphiEmanC("jjdphiEmanC",Form("jjdphi < %1.1lf",J1J2DPHI),Form("only if njets = %i",NJETS));
   selection gammaLooseVetoC("gammaLooseVetoC","photons veto");
   selection tauLooseVetoC("tauLooseVetoC","tau veto");
   // additional selections for control sample
   selection oppChargeLeptonsC("oppChargeLeptonsC","OS/SF leptons");
   selection invMassC("invMassC",Form("mass in [%3.0lf,%3.0lf]",DILEPMASS_LOW,DILEPMASS_UP));
   selection metNoLepC[metCut.size()];
   // following selections are set differently in the next "if" statements depending on the lepton flavour 
   selection lepLooseVetoC;
   selection twoLeptonsC;;
   selection twoLepLooseC;;
   selection lep1tightIdIso04C;;
   selection twoLepTightC;
   selection lep1ptC;
   selection lep2ptC;
   selection lep1etaC;  
   selection lep2etaC;
   selection genLepC;  
   // the following are only for electrons
   selection metNoEle200C;
   selection lep2tightIdIso04C;

   TVector2 l1, l2;
   TVector2 recoZTV2, metNoLepTV2;
   TVector2 met, ele, eleVectorSum;    // ele is any electron to compute MetNoEle, for muons it's not needed because it's already in the tree

   TLorentzVector l1gen, l2gen, Zgen;

   // following indices refer to the leading pair of OS/SF in the list of LepGood. They are initialized with 0 and 1 by default, but can be set with function
   // myGetPairIndexInArray (see functionsForAnalysis.cc for reference). 
   // When myGetPairIndexInArray() is called, the index of "correct" particles will be used. If they are not found (e.g. a pair of OS/SF is mismeasured as 2 mu+), 
   // indices are set as 0 and 1 (and following selection asking lep[0] and lep[1] to be OS or whatever will fail).
   Int_t firstIndex = 0;
   Int_t secondIndex = 1;

   Int_t firstIndexGen = 0;
   Int_t secondIndexGen = 1;
   Int_t recoLepFound = 0;
   Int_t genLepFound = 0;
   Int_t Z_index = 0; 

   // following 2 variable are used for acceptance and efficiency selection, define below in the loop: if selection is passed they are set to 1, otherwise they are set to 0
   Int_t acceptanceSelectionDef = 0;
   Int_t efficiencySelectionDef = 0;

   Float_t *ptr_nLepLoose = NULL;    // depending on lepton flavour in Z-->ll, it will point to different branches
   Float_t *ptr_nLep10V = NULL;   

   Float_t *ptr_metNoLep = NULL;       // only needed for muons, it will point to the branches with the metNoMu_pt, then metNoLep = *ptr_metNoLep (metNoLep defined below)
   Float_t *ptr_metNoLepPhi = NULL;  

   Float_t nLepLoose = 0.0;               // this variable and the following should be an integer, but in Emanuele's trees they are float, so I keep them as such
   Float_t nLep10V = 0.0;
   Double_t metNoLep = 0.0;        // this variable will be assigned with *ptr_metNoLep, where the pointer will point to the branch metNoMu_pt for mu, and with a hand-defined variable for e
   Double_t metNoLepPhi = 0.0;   // same story as above

   if (fabs(LEP_PDG_ID) == 13) {  // if we have Z -> mumu do stuff...
     
     strcpy(ROOT_FNAME,"zmumujetsAna.root");
     strcpy(TXT_FNAME,"zmumujetsAna.txt");
     strcpy(TEX_FNAME,"zmumujetsAna.tex");
     strcpy(FLAVOUR,"mu");
     strcpy(LL_FLAVOUR,"mumu");
     strcpy(CONTROL_SAMPLE,"Z-->mumu");
         
     ptr_nLepLoose = &nMu10V;                      // ask 2 muons
     ptr_nLep10V = &nEle10V;                         // veto on electrons
     ptr_metNoLep = &metNoMu_pt;               // for muons  get this variable from the tree 
     ptr_metNoLepPhi = &metNoMu_phi;         // for muons  get this variable from the tree

     for (Int_t i = 0; i < metCut.size(); i++) {
       metNoLepC[i].set(Form("metNoMuC[%i]",i),Form("metNoMu > %3.0lf",metCut.at(i)));
     }

     lepLooseVetoC.set("eLooseVetoC","electrons veto");
     twoLeptonsC.set("twomuonsC","muons");
     twoLepLooseC.set("twomuLooseC","2 loose muons");
     lep1tightIdIso04C.set("mu1tightIdIso04C","leading muon tight","tight ID + relIso04 (as Emanuele)");
     twoLepTightC.set("twomuTightC","2 tight muons");
     lep1ptC.set("mu1ptC",Form("mu1pt > %3.0lf",LEP1PT),"leading muon pt");
     lep2ptC.set("mu2ptC",Form("mu2pt > %3.0lf",LEP2PT),"trailing muon pt");
     lep1etaC.set("mu1etaC",Form("|mu1eta| < %1.1lf",LEP1ETA),"leading muon eta");  
     lep2etaC.set("mu2etaC",Form("|mu2eta| < %1.1lf",LEP2ETA),"trailing muon eta");
     genLepC.set("genMuonsC","muons generated");     

   } else if (fabs(LEP_PDG_ID) == 11) {   // if we have Z -> ee do different stuff...

     strcpy(ROOT_FNAME,"zeejetsAna.root");
     strcpy(TXT_FNAME,"zeejetsAna.txt");
     strcpy(TEX_FNAME,"zeejetsAna.tex");
     strcpy(FLAVOUR,"ele");
     strcpy(LL_FLAVOUR,"ee");
     strcpy(CONTROL_SAMPLE,"Z-->ee");

     ptr_nLepLoose = &nEle10V;                      // ask 2 electrons
     ptr_nLep10V = &nMu10V;                         // veto on muons   

     for (Int_t i = 0; i < metCut.size(); i++) {
       metNoLepC[i].set(Form("metNoEleC[%i]",i),Form("metNoEle > %3.0lf",metCut.at(i)));
     }
     
     lepLooseVetoC.set("muLooseVetoC","muons veto");
     twoLeptonsC.set("twoelectronsC","electrons");
     twoLepLooseC.set("twoeleLooseC","2 loose electrons");
     lep1tightIdIso04C.set("ele1tightIdIso04C","leading electron tight","tight ID + relIso04 (as Emanuele)");
     twoLepTightC.set("twoeleTightC","2 tight electrons");
     lep1ptC.set("ele1ptC",Form("ele1pt > %3.0lf",LEP1PT),"leading electron pt");
     lep2ptC.set("ele2ptC",Form("ele2pt > %3.0lf",LEP2PT),"trailing electron pt");
     lep1etaC.set("ele1etaC",Form("|ele1eta| < %1.1lf",LEP1ETA),"leading electron eta");  
     lep2etaC.set("ele2etaC",Form("|ele2eta| < %1.1lf",LEP2ETA),"trailing electron eta");
     genLepC.set("genElectronsC","electrons generated");     

     // the following are only for electrons
     metNoEle200C.set("metNoEle200C","metNoEle > 200");
     lep2tightIdIso04C.set("ele2tightIdIso04C","trailing electron tight","tight ID + relIso04 (as Emanuele)");

   }

   selection genTausC("genTausC","taus generated");                       
   selection acceptanceC("acceptanceC","acceptance cuts");
   selection efficiencyC("efficiencyC","efficiency cuts");

   selection::checkMaskLength();
   selection::printActiveSelections(cout);

   UInt_t maskMonoJetSelection = njetsEmanC.get2ToId() + jet1ptC.get2ToId() + jjdphiEmanC.get2ToId() +
                                   lepLooseVetoC.get2ToId() + gammaLooseVetoC.get2ToId();

   if ( TAU_VETO_FLAG ) maskMonoJetSelection += tauLooseVetoC.get2ToId();

   UInt_t maskTightTag;

   mask zlljetsControlSample(Form("%s control sample with selection flow as Emanuele's",CONTROL_SAMPLE));

   mask zlljetsControlSampleGenLep(Form("%s control sample (%s gen ) with selection flow as Emanuele's",CONTROL_SAMPLE,FLAVOUR));
   zlljetsControlSampleGenLep.append(genLepC.get2ToId());

   mask tautaubkgInZll(Form("tau tau background in %s control sample",CONTROL_SAMPLE));
   tautaubkgInZll.append(genTausC.get2ToId());

   if (fabs(LEP_PDG_ID) == 13) {  

     maskTightTag = lep1tightIdIso04C.get2ToId();  // for now tight requirements on pt and eta are already included in the loose condition because they coincide (not true for electrons)
   
     zlljetsControlSample.append(twoLepLooseC.get2ToId() + oppChargeLeptonsC.get2ToId());
     zlljetsControlSample.append(twoLeptonsC.get2ToId());
     zlljetsControlSample.append(maskTightTag);
     zlljetsControlSample.append(invMassC.get2ToId());
     
     zlljetsControlSampleGenLep.append(twoLepLooseC.get2ToId() + oppChargeLeptonsC.get2ToId());
     zlljetsControlSampleGenLep.append(twoLeptonsC.get2ToId());
     zlljetsControlSampleGenLep.append(maskTightTag); 
     zlljetsControlSampleGenLep.append(invMassC.get2ToId());
   
     tautaubkgInZll.append(twoLepLooseC.get2ToId() + oppChargeLeptonsC.get2ToId());
     tautaubkgInZll.append(twoLeptonsC.get2ToId());
     tautaubkgInZll.append(maskTightTag);
     tautaubkgInZll.append(invMassC.get2ToId());

   } else if (fabs(LEP_PDG_ID) == 11) {  

     maskTightTag = lep1tightIdIso04C.get2ToId() + lep2tightIdIso04C.get2ToId() + lep1ptC.get2ToId() + lep2ptC.get2ToId() + lep1etaC.get2ToId() + lep2etaC.get2ToId();

     zlljetsControlSample.append(metNoEle200C.get2ToId());
     zlljetsControlSample.append(oppChargeLeptonsC.get2ToId()); // skip loose requirement because I wil ask the tight one for both
     zlljetsControlSample.append(twoLeptonsC.get2ToId());
     zlljetsControlSample.append(maskTightTag);
     zlljetsControlSample.append(invMassC.get2ToId());
       
     zlljetsControlSampleGenLep.append(metNoEle200C.get2ToId());
     zlljetsControlSampleGenLep.append(oppChargeLeptonsC.get2ToId());
     zlljetsControlSampleGenLep.append(twoLeptonsC.get2ToId());
     zlljetsControlSampleGenLep.append(maskTightTag);
     zlljetsControlSampleGenLep.append(invMassC.get2ToId());
     
     tautaubkgInZll.append(metNoEle200C.get2ToId());
     tautaubkgInZll.append(oppChargeLeptonsC.get2ToId());
     tautaubkgInZll.append(twoLeptonsC.get2ToId());
     tautaubkgInZll.append(maskTightTag);
     tautaubkgInZll.append(invMassC.get2ToId());

   }
 
   zlljetsControlSample.append(njetsEmanC.get2ToId());
   zlljetsControlSample.append(jet1ptC.get2ToId());
   zlljetsControlSample.append(jjdphiEmanC.get2ToId());
   zlljetsControlSample.append(lepLooseVetoC.get2ToId());
   zlljetsControlSample.append(gammaLooseVetoC.get2ToId());
   
   zlljetsControlSampleGenLep.append(njetsEmanC.get2ToId());
   zlljetsControlSampleGenLep.append(jet1ptC.get2ToId());
   zlljetsControlSampleGenLep.append(jjdphiEmanC.get2ToId());
   zlljetsControlSampleGenLep.append(lepLooseVetoC.get2ToId());
   zlljetsControlSampleGenLep.append(gammaLooseVetoC.get2ToId());

   if (TAU_VETO_FLAG) {
     zlljetsControlSample.append(tauLooseVetoC.get2ToId());
     zlljetsControlSampleGenLep.append(tauLooseVetoC.get2ToId());
   }

   tautaubkgInZll.append(njetsEmanC.get2ToId());
   tautaubkgInZll.append(jet1ptC.get2ToId());
   tautaubkgInZll.append(jjdphiEmanC.get2ToId());
   tautaubkgInZll.append(lepLooseVetoC.get2ToId());
   tautaubkgInZll.append(gammaLooseVetoC.get2ToId());

   mask *lep_acc_eff[nMetBins];
   for ( Int_t i = 0; i < nMetBins; i++) {
     lep_acc_eff[i] = new mask;
     lep_acc_eff[i]->setName(Form("%s_acc_eff:  %3.0lf < met < %3.0lf",FLAVOUR,metBinEdges[i], metBinEdges[i+1]));
     lep_acc_eff[i]->append(genLepC.get2ToId());
     lep_acc_eff[i]->append(maskMonoJetSelection);
     lep_acc_eff[i]->append(acceptanceC.get2ToId());
     lep_acc_eff[i]->append(efficiencyC.get2ToId());
   }

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

   //Int_t Hcolor[] = {1,2,3,4,5,6,7,8,9,12,18,30,38,41,42,46,47,49};       

   TH1D *HzlljetsYieldsMetBin = new TH1D("HzlljetsYieldsMetBin",Form("yields of %s control sample in bins of met;#slash{E}_{T};# of events",
								     CONTROL_SAMPLE),nMetBins,metBinEdges);
   TH1D *HzlljetsYieldsMetBinGenLep = new TH1D("HzlljetsYieldsMetBinGenLep",Form("yields of %s control sample (%s gen) in bins of met; #slash{E}_{T};# of events",CONTROL_SAMPLE,CONTROL_SAMPLE),nMetBins,metBinEdges);
   TH1D *HzlljetsYieldsMetBinGenTau = new TH1D("HzlljetsYieldsMetBinGenTau",Form("yields of %s control sample (Z->#tau#tau gen) in bins of met; #slash{E}_{T};# of events",CONTROL_SAMPLE),nMetBins,metBinEdges);
   TH1D *HzvvEstimate = new TH1D("HzvvEstimate",Form("yields of Z->#nu#nu estimated as N(%s) * BR_ratio / (A*#varepsilon)",CONTROL_SAMPLE),nMetBins,metBinEdges);

   TH1D *HZtoLLRecoPt = new TH1D("HZtoLLRecoPt","",101,0.,1010);
   TH1D *HZtoLLGenPt = new TH1D("HZtoLLGenPt","",101,0.,1010);
   // this is the histogram with reco/gen
   TH1D *HZtoLLPt_RecoGenRatio = new TH1D("HZtoLLPt_RecoGenRatio","",101,0.,1010.);
   // histogram of reco/gen distribution function
   TH1D *HZtoLLPt_RecoGenRatio_pdf = new TH1D("HZtoLLPt_RecoGenRatio_pdf","",100,0.,2.);

   TH1D* Hacc = new TH1D("Hacc","",nMetBins,metBinEdges);
   TH1D* Heff = new TH1D("Heff","",nMetBins,metBinEdges);
   TH1D* Hacceff = new TH1D("Hacceff","",nMetBins,metBinEdges);

   TH1D *HinvMass[nMetBins];
   TH1D *HzlljetsInvMassMetBinGenLep[nMetBins];
   TH1D *HzlljetsInvMassMetBinGenTau[nMetBins];
   TH1D *HZtoLLRecoPt_MetBin[nMetBins];
   TH1D *HZtoLLGenPt_MetBin[nMetBins];
   TH1D *HZtoLLPt_RecoGenRatio_MetBin[nMetBins];
   TH1D *HZtoLLPt_RecoGenRatio_pdf_MetBin[nMetBins];

   for (Int_t i = 0; i < nMetBins; i++) {

     HinvMass[i] = new TH1D(Form("HinvMass_met%2.0lfTo%2.0lf",metBinEdges[i],metBinEdges[i+1]),"",30,DILEPMASS_LOW,DILEPMASS_UP);
     HzlljetsInvMassMetBinGenLep[i] = new TH1D(Form("HzlljetsInvMassMetBinGenLep_met%2.0lfTo%2.0lf",metBinEdges[i],metBinEdges[i+1]),"",30,DILEPMASS_LOW,DILEPMASS_UP);
     HzlljetsInvMassMetBinGenTau[i] = new TH1D(Form("HzlljetsInvMassMetBinGenTau_met%2.0lfTo%2.0lf",metBinEdges[i],metBinEdges[i+1]),"",30,DILEPMASS_LOW,DILEPMASS_UP);
     HZtoLLRecoPt_MetBin[i] = new TH1D(Form("HZtoLLRecoPt_MetBin_met%2.0lfTo%2.0lf",metBinEdges[i],metBinEdges[i+1]),"",101,0.,1010.);
     HZtoLLGenPt_MetBin[i] = new TH1D(Form("HZtoLLGenPt_MetBin_met%2.0lfTo%2.0lf",metBinEdges[i],metBinEdges[i+1]),"",101,0.,1010.);
     HZtoLLPt_RecoGenRatio_MetBin[i] = new TH1D(Form("HZtoLLPt_RecoGenRatio_MetBin_met%2.0lfTo%2.0lf",metBinEdges[i],metBinEdges[i+1]),"",101,0.,1010.);
     HZtoLLPt_RecoGenRatio_pdf_MetBin[i] = new TH1D(Form("HZtoLLPt_RecoGenRatio_pdf_MetBin_met%2.0lfTo%2.0lf",metBinEdges[i],metBinEdges[i+1]),"",100,0.,2.);

   } 

   TH1D *H_uPerp_VS_Nvtx[NVTXS];
   TH1D *H_uPar_VS_Nvtx[NVTXS]; 
 
   for (Int_t i = 0; i < NVTXS; i++) {

     H_uPerp_VS_Nvtx[i] = new TH1D(Form("H_uPerp_VS_Nvtx_nvtx%i",FIRST_NVTX+i),"",80,-200,200);  // 5 GeV bins 
     H_uPar_VS_Nvtx[i] = new TH1D(Form("H_uPar_VS_Nvtx_nvtx%i",FIRST_NVTX+i),"",80,-200,200);  // 5 GeV bins

   }

   //Double_t ZptBinEdges[] = {250., 260., 270., 280., 290., 310., 330., 350., 370., 390., 410., 430., 450., 470., 500., 530., 560, 600., 640., 700., 800.};
   Double_t ZptBinEdges[] = {250., 260., 270., 280., 290., 310., 330., 350., 370., 400., 430., 460., 490., 530., 570, 610., 650., 700., 800.};
   Int_t nBinsForResponse = sizeof(ZptBinEdges)/sizeof(Double_t) - 1;  //number of bins is n-1 where n is the number of ZptBinEdges's elements

   TH1D *H_uPerp_VS_ZpT[nBinsForResponse];  
   TH1D *H_uPar_VS_ZpT[nBinsForResponse]; 
   TH1D *H_uPar_ZpT_ratio[nBinsForResponse];  // for the response curve
   TH1D *HZptBinned[nBinsForResponse];
   //the following histograms will give the distribution of met|| / wzpt. The mean value will be used to create the response curve, that is (<met|| / wzpt>) vs wzpt
   // for each point, wzpt will be taken as the average wzpt in the range considered
 
   for (Int_t i = 0; i < nBinsForResponse; i++) {   

     //HZptBinned[i] are histograms with 5 bins in the range given by ZptBinEdges[i] and ZptBinEdges[i+1]
     // the mean wzpt in each bin will be computed as the histogram's mean
     HZptBinned[i] = new TH1D(Form("HZptBinned_ZpT%2.0lfTo%2.0lf",ZptBinEdges[i],ZptBinEdges[i+1]),"",5,ZptBinEdges[i],ZptBinEdges[i+1]); 
     H_uPar_ZpT_ratio[i] = new TH1D(Form("H_uPar_ZpT_ratio_ZpT%2.0lfTo%2.0lf",ZptBinEdges[i],ZptBinEdges[i+1]),"",50,0.0,2.0); 
     H_uPerp_VS_ZpT[i] = new TH1D(Form("H_uPerp_VS_ZpT_ZpT%2.0lfTo%2.0lf",ZptBinEdges[i],ZptBinEdges[i+1]),"",40,-200,200); 
     H_uPar_VS_ZpT[i] = new TH1D(Form("H_uPar_VS_ZpT_ZpT%2.0lfTo%2.0lf",ZptBinEdges[i],ZptBinEdges[i+1]),"",40,-200,200); 

   }


   // saving histograms with bin edges of other histograms used (e.g. content of metBinEdges array ...)
   TH1D *HmetBinEdges = new TH1D("HmetBinEdges","bin edges for met distributions",nMetBins+1,0.0,nMetBins+1);
   for (Int_t i = 0; i <= nMetBins; i++) {
     HmetBinEdges->SetBinContent(i+1,metBinEdges[i]);
   }

   TH1D *HZptBinEdges = new TH1D("HZptBinEdges","bin edges for ZpT distributions",nBinsForResponse+1,0.0,nBinsForResponse+1);
   for (Int_t i = 0; i <= nBinsForResponse; i++) {
     HZptBinEdges->SetBinContent(i+1,ZptBinEdges[i]);
   }

   TH1D *Hnvtx = new TH1D("Hnvtx","# of vertices for studies of variables as a function of nvtx",NVTXS,0.0,NVTXS);
   for (Int_t i = 0; i < NVTXS; i++) {              // watch out: differently from above, i < NVTXS, not <=, because if NVTXS = 3 I need 3 points, not 4
     Hnvtx->SetBinContent(i+1,FIRST_NVTX+i);
   }

   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"zlljetsAna_new::loop()"<<endl;
   cout<<"nentries = "<<nentries<<endl;   

   Long64_t nbytes = 0, nb = 0;

   for (Int_t jentry=0; jentry<nentries; jentry++) {

     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;   

     UInt_t eventMask = 0; 
     Double_t newwgt = weight * LUMI;

     nLepLoose = *ptr_nLepLoose;          
     nLep10V = *ptr_nLep10V;

     genLepFound = myPartGenAlgo(nGenPart, GenPart_pdgId, GenPart_motherId, LEP_PDG_ID, 23, firstIndexGen, secondIndexGen, Z_index, GenPart_motherIndex); 
     l1gen.SetPtEtaPhiM(GenPart_pt[firstIndexGen],GenPart_eta[firstIndexGen],GenPart_phi[firstIndexGen],GenPart_mass[firstIndexGen]);
     l2gen.SetPtEtaPhiM(GenPart_pt[secondIndexGen],GenPart_eta[secondIndexGen],GenPart_phi[secondIndexGen],GenPart_mass[secondIndexGen]);
     Zgen = l1gen + l2gen;
     Double_t ZgenMass = Zgen.Mag();
     recoLepFound = myGetPairIndexInArray(LEP_PDG_ID, nLepGood, LepGood_pdgId, firstIndex, secondIndex);
     // Z_PDGID = 23
     Z_index = GenPart_motherIndex[firstIndexGen];
     //Z_index = myGetPartIndex(23, nGenPart, GenPart_pdgId);
     // I find the indices corresponding to the 2 leading lepton
     //cout<<"entry : "<<jentry<<endl;
     l1.SetMagPhi(LepGood_pt[firstIndex],LepGood_phi[firstIndex]);
     l2.SetMagPhi(LepGood_pt[secondIndex],LepGood_phi[secondIndex]);
     recoZTV2 = l1 + l2;
     Double_t ZtoLLRecoPt = recoZTV2.Mod();
     //Double_t ZtoLLGenPt = GenPart_pt[Z_index];
     Double_t ZtoLLGenPt = Zgen.Pt();

     if (fabs(LEP_PDG_ID) == 13) { 

       if ( (GenPart_pt[firstIndexGen] > GENLEP1PT) && (GenPart_pt[secondIndexGen] > GENLEP2PT) && ( fabs(GenPart_eta[firstIndexGen]) < GENLEP1ETA) && ( fabs(GenPart_eta[secondIndexGen]) < GENLEP2ETA) && (ZgenMass > GEN_ZMASS_LOW) && (ZgenMass < GEN_ZMASS_UP) )  acceptanceSelectionDef = 1;
       else acceptanceSelectionDef = 0;

       if ( (nLepLoose == 2) && (LepGood_tightId[firstIndex] == 1) && (LepGood_relIso04[firstIndex] < LEP_ISO_04 ) && (fabs(LepGood_pdgId[firstIndex]) == LEP_PDG_ID) ) efficiencySelectionDef = 1;
       else efficiencySelectionDef = 0;

       metNoLep = *ptr_metNoLep;         
       metNoLepPhi = *ptr_metNoLepPhi;    

     } else if (fabs(LEP_PDG_ID) == 11) { 

       if ( HLT_FLAG ) {

	 if ( !( (LepGood_tightId[firstIndex] == 1) && (LepGood_tightId[secondIndex] == 1) && 
		 ( (LepGood_pdgId[firstIndex] * LepGood_pdgId[secondIndex] ) == NEG_LEP_PDG_ID2) && 
		 (fabs(LepGood_eta[firstIndex]) < HLT_LEP1ETA) && (fabs(LepGood_eta[secondIndex]) < HLT_LEP2ETA) && 
		 (LepGood_pt[secondIndex] > HLT_LEP1PT) && (LepGood_pt[secondIndex] > HLT_LEP2PT) ) )  continue;  //if HLT is not passed, skip everything in loop 

       }

       met.SetMagPhi(met_pt,met_phi);
       eleVectorSum.SetMagPhi(0.0,0.0);   // for each event it must be initialized to 0

       for (Int_t i = 0; i < nLepGood; i++) {
	 if (fabs(LepGood_pdgId[i]) == LEP_PDG_ID) {
	   ele.SetMagPhi(LepGood_pt[i],LepGood_phi[i]);
	   eleVectorSum += ele;
	 }
       }

       met += eleVectorSum;
       metNoLep = met.Mod();  // for electrons we define it by hand, for muons we use the variable in the tree
       metNoLepPhi = met.Phi();

       if ( (GenPart_pt[firstIndexGen] > GENLEP1PT) && (GenPart_pt[secondIndexGen] > GENLEP2PT) &&
       	    ( fabs(GenPart_eta[firstIndexGen]) < GENLEP1ETA) && ( fabs(GenPart_eta[secondIndexGen]) < GENLEP2ETA) &&
       	    (ZgenMass > GEN_ZMASS_LOW) && (ZgenMass < GEN_ZMASS_UP) ) acceptanceSelectionDef = 1;
       else acceptanceSelectionDef = 0;

       if ( (LepGood_tightId[firstIndex] == 1) && (LepGood_tightId[secondIndex] == 1) &&
       	    (LepGood_relIso04[firstIndex] < LEP_ISO_04 ) && (LepGood_relIso04[secondIndex] < LEP_ISO_04 ) &&
       	    ((LepGood_pdgId[firstIndex] * LepGood_pdgId[secondIndex]) == NEG_LEP_PDG_ID2) ) efficiencySelectionDef = 1;
       else efficiencySelectionDef = 0;

     }

     if (recoLepFound) {  // following is done if two OS leptons are found (otherwise there would be no Z)
		
       //recoZTV2.SetMagPhi(ZtoLLRecoPt, (l1 + l2).Phi());
       metNoLepTV2.SetMagPhi((Double_t)metNoLep,(Double_t)metNoLepPhi);
       Double_t dphiMetNoLepZ = metNoLepTV2.DeltaPhi(recoZTV2);

       Double_t u_par = metNoLep * TMath::Cos(dphiMetNoLepZ);
       Double_t u_perp = metNoLep * TMath::Sin(dphiMetNoLepZ);

       if (ZtoLLRecoPt > 250) {   //this corresponds to trigger efficiency plateaux (actually it would be mumet, but mumet recoils against wzpt so their pt should be of the same order of magnitude). In any case, note that in our tree is always mumet > 200 (for electrons there is no such a cut)

	 Int_t nvtxBin = nVert - FIRST_NVTX;
	 Int_t lastnvtx = NVTXS + FIRST_NVTX;

	 if ((nvtxBin >= 0) && (nVert < lastnvtx)) {

	   if (ZtoLLRecoPt < 500) {                       // (met||-wzpt) distribution's width depends on pt, thus I use this range

	     H_uPar_VS_Nvtx[nvtxBin]->Fill(u_par - ZtoLLRecoPt,newwgt);
	 
	   }

	   H_uPerp_VS_Nvtx[nvtxBin]->Fill(u_perp,newwgt);       
    
	 }

       }            // end of if (ZtoLLRecoPt > 250)
	
       /**************************************************/
       // computing met responses
       /**************************************************/

       // first of all I make sure that wzpt is in the appropriate range
       if ((ZtoLLRecoPt > ZptBinEdges[0]) && (ZtoLLRecoPt < ZptBinEdges[nBinsForResponse])) {

	 Int_t respBin = myGetBin(ZtoLLRecoPt,ZptBinEdges,nBinsForResponse);
	 //cout<<"bin = "<<bin<<endl;
	 HZptBinned[respBin]->Fill(ZtoLLRecoPt,newwgt);        
	 H_uPar_ZpT_ratio[respBin]->Fill(u_par/ZtoLLRecoPt,newwgt);     //the mean value of this histogram is the response
	 H_uPerp_VS_ZpT[respBin]->Fill(u_perp,newwgt);
	 H_uPar_VS_ZpT[respBin]->Fill(u_par - ZtoLLRecoPt,newwgt);

       }

     }		// end of if (recoLepFound)
		 
     nTotalWeightedEvents += newwgt;  // counting events with weights

     eventMask += njetsC.addToMask(nJet30a <= NJETS);
     // jetclean is 1 if cleaning is passed, 0 otherwise. It's applied to first jet and , if any, to the second
     eventMask += njetsEmanC.addToMask( (nJet30 == 1 || nJet30 == 2) && jetclean > 0.5);
     //eventMask += njetsEmanC.addToMask( ((nJet30a == 1 ) || (nJet30a == 2 && Jet_eta[1] < JET2ETA)) && Jet_eta[0] < JET1ETA && jetclean > 0.5);
     eventMask += jjdphiEmanC.addToMask( nJet30 == 1 || (nJet == 2 && fabs(dphijj) < J1J2DPHI));
     //eventMask += jjdphiEmanC.addToMask( (nJet30a == 1 && Jet_eta[0] < 2.5) || (nJet30a == 2 && abs(dphijj) < J1J2DPHI));
     eventMask += jet1ptC.addToMask(Jet_pt[0] > J1PT);               
     //in Emanuele's tree we have vectors: [0] is the first jet, [1] is the second and so on (ordered in pt)
     eventMask += jet1etaC.addToMask(fabs(Jet_eta[0]) < J1ETA);

     eventMask += lepLooseVetoC.addToMask(nLep10V == 0);
     eventMask += gammaLooseVetoC.addToMask(nGamma15V == 0);
     eventMask += tauLooseVetoC.addToMask(nTau15V == 0);
     for (Int_t i = 0; i <  metCut.size(); i++) {
       eventMask += metNoLepC[i].addToMask(metNoLep > metCut[i]);
     }
     eventMask += oppChargeLeptonsC.addToMask( (LepGood_pdgId[firstIndex] * LepGood_pdgId[secondIndex]) == NEG_LEP_PDG_ID2);
     eventMask += twoLeptonsC.addToMask((fabs(LepGood_pdgId[firstIndex]) == LEP_PDG_ID) && (fabs(LepGood_pdgId[secondIndex]) == LEP_PDG_ID));
     eventMask += twoLepLooseC.addToMask(nLepLoose == 2);
     eventMask += lep1tightIdIso04C.addToMask((LepGood_tightId[firstIndex] == 1) && (LepGood_relIso04[firstIndex] < LEP_ISO_04 ) && (fabs(LepGood_pdgId[firstIndex]) == LEP_PDG_ID));
 
     if (fabs(LEP_PDG_ID) == 11) { 

       eventMask += lep2tightIdIso04C.addToMask((LepGood_tightId[secondIndex] == 1) && (LepGood_relIso04[secondIndex] < LEP_ISO_04 ) && (fabs(LepGood_pdgId[secondIndex]) == LEP_PDG_ID));
       eventMask += metNoEle200C.addToMask(metNoLep > 200);
       //if ((jentry%50000) == 0) cout << "jentry = " <<jentry<<endl;
       
     }

     eventMask += lep1ptC.addToMask((LepGood_pt[firstIndex] > LEP1PT) && (fabs(LepGood_pdgId[firstIndex]) == LEP_PDG_ID)); 
     eventMask += lep1etaC.addToMask( (fabs(LepGood_eta[firstIndex]) < LEP1ETA) && (fabs(LepGood_pdgId[firstIndex]) == LEP_PDG_ID) );
     eventMask += lep2ptC.addToMask((LepGood_pt[secondIndex] > LEP2PT) && (fabs(LepGood_pdgId[secondIndex]) == LEP_PDG_ID));
     eventMask += lep2etaC.addToMask((fabs(LepGood_eta[secondIndex]) < LEP2ETA) && (fabs(LepGood_pdgId[secondIndex]) == LEP_PDG_ID));
     eventMask += invMassC.addToMask((mZ1 > DILEPMASS_LOW) && (mZ1 < DILEPMASS_UP));
     eventMask += genLepC.addToMask( genLepFound );
     eventMask += genTausC.addToMask( myPartGenAlgo(nGenPart, GenPart_pdgId, GenPart_motherId, 15, 23) );  // tau pdg id = 15, Z pdg id = 23 
     eventMask += acceptanceC.addToMask( acceptanceSelectionDef );
     eventMask += efficiencyC.addToMask( efficiencySelectionDef );

     zlljetsControlSample.countEvents(eventMask,newwgt);
     zlljetsControlSampleGenLep.countEvents(eventMask,newwgt);
     tautaubkgInZll.countEvents(eventMask, newwgt);

     // filling histogram with yields and invariant mass at the end of the selection in bins of met
     if ( ((eventMask & zlljetsControlSample.globalMask.back()) == zlljetsControlSample.globalMask.back()) ) {  
       // this histogram holds the final yields in bins of MET
       HzlljetsYieldsMetBin->Fill(metNoLep,newwgt);    
     }

     if ( ((eventMask & zlljetsControlSampleGenLep.globalMask.back()) == zlljetsControlSampleGenLep.globalMask.back()) ) {
       // this histogram holds the final yields in bins of MET
       HzlljetsYieldsMetBinGenLep->Fill(metNoLep,newwgt);
       HZtoLLRecoPt->Fill(ZtoLLRecoPt,newwgt);
       HZtoLLGenPt->Fill(ZtoLLGenPt,newwgt);
       if (ZtoLLGenPt != 0) HZtoLLPt_RecoGenRatio_pdf->Fill(ZtoLLRecoPt/ZtoLLGenPt,newwgt);
     }

     if ( ((eventMask & tautaubkgInZll.globalMask.back()) == tautaubkgInZll.globalMask.back()) ) {  
       // this histogram holds the final yields in bins of MET
       HzlljetsYieldsMetBinGenTau->Fill(metNoLep,newwgt);  
     }

     if ((metNoLep > metBinEdges[0]) && (metNoLep < metBinEdges[nMetBins])) {

       Int_t bin = myGetBin(metNoLep,metBinEdges,nMetBins);
       lep_acc_eff[bin]->countEvents(eventMask,newwgt);
       
       if ((eventMask & zlljetsControlSample.globalMask.back()) == zlljetsControlSample.globalMask.back()) {
	 // this histogram holds the invariant mass distribution (one for each met bin)
	 HinvMass[bin]->Fill(mZ1,newwgt);   
       }

       if ( ((eventMask & zlljetsControlSampleGenLep.globalMask.back()) == zlljetsControlSampleGenLep.globalMask.back()) ) {  
	 HzlljetsInvMassMetBinGenLep[bin]->Fill(mZ1,newwgt); 
	 HZtoLLRecoPt_MetBin[bin]->Fill(ZtoLLRecoPt,newwgt);
	 HZtoLLGenPt_MetBin[bin]->Fill(ZtoLLGenPt,newwgt);
	 if (ZtoLLGenPt != 0) HZtoLLPt_RecoGenRatio_pdf_MetBin[bin]->Fill(ZtoLLRecoPt/ZtoLLGenPt,newwgt);
       }

       if ( ((eventMask & tautaubkgInZll.globalMask.back()) == tautaubkgInZll.globalMask.back()) ) {  
	 HzlljetsInvMassMetBinGenTau[bin]->Fill(mZ1,newwgt);   
       }
     
     } 

   }     // end of loop on entries

   
   /************************************/
   //                    MET|| & MET_|_ VS NVTX & ZpT
   /************************************/

   //resolution vs nvtx

   Double_t xValues[NVTXS];
   Double_t yValues[NVTXS];
   Double_t yValuesErr[NVTXS];
   for (Int_t i = 0; i < NVTXS; i++) {
     xValues[i] = i + FIRST_NVTX;
     yValues[i] = H_uPar_VS_Nvtx[i]->GetRMS();
     yValuesErr[i] = H_uPar_VS_Nvtx[i]->GetRMSError();
   }

   TGraphErrors *GresolutionMetNoLepParZvsNvtx = new TGraphErrors(NVTXS,xValues,yValues,0,yValuesErr);
   GresolutionMetNoLepParZvsNvtx->SetTitle("resolution || from histogram's RMS");
   GresolutionMetNoLepParZvsNvtx->Draw("AP");
   GresolutionMetNoLepParZvsNvtx->SetMarkerStyle(7);  // 7 is a medium dot
   GresolutionMetNoLepParZvsNvtx->GetXaxis()->SetTitle("nvtx");
   GresolutionMetNoLepParZvsNvtx->GetYaxis()->SetTitle("#sigma (u_{||}) [GeV]");
   GresolutionMetNoLepParZvsNvtx->GetYaxis()->SetTitleOffset(1.4); 
   GresolutionMetNoLepParZvsNvtx->SetName("gr_resolution_uPar_vs_Nvtx");
   GresolutionMetNoLepParZvsNvtx->Write();

   for (Int_t i = 0; i < NVTXS; i++) {
     yValues[i] = H_uPerp_VS_Nvtx[i]->GetRMS();
     yValuesErr[i] = H_uPerp_VS_Nvtx[i]->GetRMSError();
   }

   TGraphErrors *GresolutionMetNoLepOrtZvsNvtx = new TGraphErrors(NVTXS,xValues,yValues,0,yValuesErr);
   GresolutionMetNoLepOrtZvsNvtx->SetTitle("resolution _|_ from histogram's RMS");
   GresolutionMetNoLepOrtZvsNvtx->Draw("AP");
   GresolutionMetNoLepOrtZvsNvtx->SetMarkerStyle(7);
   GresolutionMetNoLepOrtZvsNvtx->GetXaxis()->SetTitle("nvtx");
   GresolutionMetNoLepOrtZvsNvtx->GetYaxis()->SetTitle("#sigma (u_#perp ) [GeV]");
   GresolutionMetNoLepOrtZvsNvtx->GetYaxis()->SetTitleOffset(1.4); 
   GresolutionMetNoLepOrtZvsNvtx->SetName("gr_resolution_uPerp_vs_Nvtx");
   GresolutionMetNoLepOrtZvsNvtx->Write();

   // response curve

   Double_t response[nBinsForResponse];
   Double_t responseErr[nBinsForResponse];
   Double_t meanZpt[nBinsForResponse];
   Double_t meanZptErr[nBinsForResponse];

   for (Int_t i = 0; i < nBinsForResponse; i++) {
     meanZpt[i] = HZptBinned[i]->GetMean();
     meanZptErr[i] = HZptBinned[i]->GetMeanError();
     response[i] = H_uPar_ZpT_ratio[i]->GetMean();
     responseErr[i] = H_uPar_ZpT_ratio[i]->GetMeanError();
     //cout<<i<<" meanZpt = "<<meanZpt[i]<<" +/- "<<meanZptErr[i]<<"    response = "<<response[i]<<" +/- "<<responseErr[i]<<endl;
   }

   TGraphErrors *GresponseCurve = new TGraphErrors(nBinsForResponse,meanZpt,response,0,responseErr);
   GresponseCurve->SetTitle("response curve");
   GresponseCurve->Draw("AP");
   GresponseCurve->SetMarkerStyle(7);    // 7 is a medium dot
   GresponseCurve->GetXaxis()->SetTitle("ZpT [GeV]");
   GresponseCurve->GetYaxis()->SetTitle(" < u_{||} / ZpT >");
   GresponseCurve->GetYaxis()->SetRangeUser(0.6, 1.1);
   GresponseCurve->GetYaxis()->SetTitleOffset(1.4); 
   GresponseCurve->SetName("gr_responseCurve");
   GresponseCurve->Write();

   // resolution vs ZpT

   Double_t resoMetNoLepParZvsZpt[nBinsForResponse];
   Double_t resoMetNoLepParZvsZptErr[nBinsForResponse];
   Double_t resoMetNoLepOrtZvsZpt[nBinsForResponse];
   Double_t resoMetNoLepOrtZvsZptErr[nBinsForResponse];

   for (Int_t i = 0; i < nBinsForResponse; i++) {
     resoMetNoLepParZvsZpt[i] = H_uPar_VS_ZpT[i]->GetRMS();
     resoMetNoLepParZvsZptErr[i] = H_uPar_VS_ZpT[i]->GetRMSError();
     resoMetNoLepOrtZvsZpt[i] = H_uPerp_VS_ZpT[i]->GetRMS();
     resoMetNoLepOrtZvsZptErr[i] = H_uPerp_VS_ZpT[i]->GetRMSError();
   }

   TGraphErrors *GresolutionMetNoLepParZvsZpt = new TGraphErrors(nBinsForResponse,meanZpt,resoMetNoLepParZvsZpt,0,resoMetNoLepParZvsZptErr);
   GresolutionMetNoLepParZvsZpt->SetTitle("resolution || from histogram's RMS");
   GresolutionMetNoLepParZvsZpt->Draw("AP");
   GresolutionMetNoLepParZvsZpt->SetMarkerStyle(7);  // 7 is a medium dot
   GresolutionMetNoLepParZvsZpt->GetXaxis()->SetTitle("Zpt [GeV]");
   GresolutionMetNoLepParZvsZpt->GetYaxis()->SetTitle("#sigma (u_{||}) [GeV]");
   GresolutionMetNoLepParZvsZpt->GetYaxis()->SetTitleOffset(1.2); 
   GresolutionMetNoLepParZvsZpt->SetName("gr_resolution_uPar_vs_ZpT");
   GresolutionMetNoLepParZvsZpt->Write();

   TGraphErrors *GresolutionMetNoLepOrtZvsZpt = new TGraphErrors(nBinsForResponse,meanZpt,resoMetNoLepOrtZvsZpt,0,resoMetNoLepOrtZvsZptErr);
   GresolutionMetNoLepOrtZvsZpt->SetTitle("resolution _|_ from histogram's RMS");
   GresolutionMetNoLepOrtZvsZpt->Draw("AP");
   GresolutionMetNoLepOrtZvsZpt->SetMarkerStyle(7);
   GresolutionMetNoLepOrtZvsZpt->GetXaxis()->SetTitle("Zpt [GeV]");
   GresolutionMetNoLepOrtZvsZpt->GetYaxis()->SetTitle("#sigma (u_#perp ) [GeV]");
   GresolutionMetNoLepOrtZvsZpt->GetYaxis()->SetTitleOffset(1.2); 
   GresolutionMetNoLepOrtZvsZpt->SetName("gr_resolution_uPerp_vs_ZpT");
   GresolutionMetNoLepOrtZvsZpt->Write();

   // end of TGraphs

   mySpaces(cout,2);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &zlljetsControlSample);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &zlljetsControlSampleGenLep);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &tautaubkgInZll);
   // for (Int_t i = 0; i < nMetBins; i++) {
   //   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, lep_acc_eff[i] );
   // }

   mySpaces(cout,2);
   myPrintYieldsMetBinInStream(cout, HzlljetsYieldsMetBinGenLep, metBinEdges, nMetBins);
 
   cout<<"creating file '"<<TXT_FNAME<<"' ..."<<endl;
   ofstream myfile(TXT_FNAME,ios::out);

   if ( !myfile.is_open() ) {

     cout<<"Error: unable to open file "<<TXT_FNAME<<" !"<<endl;
     exit(EXIT_FAILURE);
     
   }
          
   selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &zlljetsControlSample);
   selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &zlljetsControlSampleGenLep);
   selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &tautaubkgInZll);
   mySpaces(myfile,2);
   myPrintYieldsMetBinInStream(myfile, HzlljetsYieldsMetBinGenLep, metBinEdges, nMetBins);

   Int_t stepMonojetSelection_In_lepAccEff = lep_acc_eff[0]->whichStepHas(maskMonoJetSelection);
   Int_t stepAcceptance_In_lepAccEff = lep_acc_eff[0]->whichStepHas(acceptanceC.get2ToId());
   Int_t stepEfficiency_In_lepAccEff = lep_acc_eff[0]->whichStepHas(efficiencyC.get2ToId());
   // cout<<"step: MJ     acc     eff"<<endl;
   // cout<<stepMonojetSelection_In_lepAccEff<<stepAcceptance_In_lepAccEff<<stepEfficiency_In_lepAccEff<<endl;
   Double_t acc, eff, accStatErr, effStatErr, acceff, acceffStatErr;

   mySpaces(cout,2);
   mySpaces(myfile,2);
   cout << "Printing acceptance and efficiency." << endl;
   cout << "MET [GeV]     acc     acc_unc     eff     eff_unc" <<endl;
   myfile << "MET [GeV]     acc     acc_unc     eff     eff_unc" <<endl;

   for (Int_t i = 0; i < nMetBins; i++) {
  
     acc = lep_acc_eff[i]->nEvents[stepAcceptance_In_lepAccEff]/lep_acc_eff[i]->nEvents[stepMonojetSelection_In_lepAccEff];
     eff = lep_acc_eff[i]->nEvents[stepEfficiency_In_lepAccEff]/lep_acc_eff[i]->nEvents[stepAcceptance_In_lepAccEff];
     accStatErr = sqrt(acc * (1 - acc) / lep_acc_eff[i]->nEvents[stepMonojetSelection_In_lepAccEff]);
     effStatErr = sqrt(eff * (1 - eff) / lep_acc_eff[i]->nEvents[stepAcceptance_In_lepAccEff]);
     Hacc->SetBinContent(i+1,acc);
     Hacc->SetBinError(i+1,accStatErr);
     Heff->SetBinContent(i+1,eff);
     Heff->SetBinError(i+1,effStatErr);

     cout<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" "<<acc<<" "<<accStatErr<<" "<<eff<<" "<<effStatErr<<endl;
     myfile<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" "<<acc<<" "<<accStatErr<<" "<<eff<<" "<<effStatErr<<endl;

   }

   mySpaces(cout,2);
   mySpaces(myfile,2);
   cout << "Printing acceptance * efficiency" << endl;
   cout << "MET [GeV]     acc*eff     acc*eff_unc" <<endl;
   myfile << "MET [GeV]     acc*eff     acc*eff_unc" <<endl;

   for (Int_t i = 0; i < nMetBins; i++) {
     // do not merge with previous loop: I want to print them after the previous loop because I might copy and paste this output to make acc * eff table
     
     acceff = lep_acc_eff[i]->nEvents[stepEfficiency_In_lepAccEff]/lep_acc_eff[i]->nEvents[stepMonojetSelection_In_lepAccEff];
     acceffStatErr = sqrt(acceff * (1 - acceff) / lep_acc_eff[i]->nEvents[stepMonojetSelection_In_lepAccEff]);
     Hacceff->SetBinContent(i+1,acceff);
     Hacceff->SetBinError(i+1,acceffStatErr);

     cout<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" "<<acceff<<" "<<acceffStatErr<<endl;
     myfile<<(Int_t)metBinEdges[i]<<"-"<<(Int_t)metBinEdges[i+1]<<" "<<acceff<<" "<<acceffStatErr<<endl;
  
   }
   
   // now get Z(inv) estimate as N_Zvv = N_Zll * R / (A*e), R being BR(Zvv)/BR(Zll) where l is either mu or e (R ~ 6)
   TH1D *H_BR_ratio = new TH1D("H_BR_ratio",Form("BR(Z#nu#nu/BR(%s)",CONTROL_SAMPLE),nMetBins,metBinEdges);

   for(Int_t i = 0; i <= nMetBins; i++) {
     H_BR_ratio->SetBinContent(i,RATIO_BR_ZINV_ZLL);
     H_BR_ratio->SetBinError(i,UNC_RATIO_BR_ZINV_ZLL);
   }

   HzvvEstimate->Multiply(HzlljetsYieldsMetBinGenLep,H_BR_ratio);
   HzvvEstimate->Divide(Hacceff);
   delete H_BR_ratio; //no need to save it

   // I add overflow bin's content in the last bin for all histograms where that is needed
   // for those histogram filled with Divide() method, it's not done as long as it was already done on the histograms given as
   // argument to the Divide() method
   myAddOverflowInLastBin(HZtoLLRecoPt);
   myAddOverflowInLastBin(HZtoLLGenPt);  
   HZtoLLPt_RecoGenRatio->Divide(HZtoLLRecoPt,HZtoLLGenPt);

   for (Int_t i = 0; i < nMetBins; i++) {
       myAddOverflowInLastBin(HZtoLLRecoPt_MetBin[i]);
       myAddOverflowInLastBin(HZtoLLGenPt_MetBin[i]);
       HZtoLLPt_RecoGenRatio_MetBin[i]->Divide(HZtoLLRecoPt_MetBin[i],HZtoLLGenPt_MetBin[i]);  
   }

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
     makeTableTex(fp, LUMI, nTotalWeightedEvents, &zlljetsControlSample,commentInTable);
     makeTableTex(fp, LUMI, nTotalWeightedEvents, &zlljetsControlSampleGenLep,commentInTable);
     makeTableTex(fp, LUMI, nTotalWeightedEvents, &tautaubkgInZll,commentInTable);
     fprintf(fp,"\\end{document}\n");      
     fclose(fp);

   }

   // end of tex file

}



