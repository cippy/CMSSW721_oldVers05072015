#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <cmath>
#include <math.h>

#include <TTree.h>
#include <TChain.h>
#include <TROOT.h>
#include <TFile.h>
#include "EmanTreeAnalysis.h"
#include "whichApplication.h"

using namespace myAnalyzerTEman;

using namespace std;


int main(int argc, char* argv[]) {


//================ Creating chain 
 
#if Application == 1 

  if (argc < 2) {
    std::cout << "Not enough arguments: launch as -> "; 
    std::cout << argv[0] << " inputfile " <<std::endl;
    exit(EXIT_FAILURE);
  }
  std::cout<<"Using Emanuele's trees with skim on mumet > 200"<<std::endl;

  char inputFileName[200];
  char buffer1[200];
  char buffer2[200];
  char rootFileToChain[500];
  char rootFriendFileToChain[500];

  //vector< vector<Double_t> > matrix;
  vector< Double_t > yieldsRow;
  vector< Double_t > efficiencyRow;
  Int_t nSample = 0;
  std::vector<std::string> sampleName;
  sampleName.push_back("DYJetsToLL");
  sampleName.push_back("Top");
  sampleName.push_back("QCD");
  sampleName.push_back("GJets");
  sampleName.push_back("WJetsToLNu");
  sampleName.push_back("ZJetsToNuNu");
  std::vector<std::string> selectionDefinition;
  selectionDefinition.push_back("entry point");
  selectionDefinition.push_back("2 lept SF/OS");
  selectionDefinition.push_back("muons");
  selectionDefinition.push_back("tight Tag");
  selectionDefinition.push_back("mll");
  selectionDefinition.push_back("njets");
  selectionDefinition.push_back("jet1pt > 110");
  selectionDefinition.push_back("jetjetdphi");
  selectionDefinition.push_back("electron veto");
  selectionDefinition.push_back("photon veto");

  strcpy(inputFileName,argv[1]);

  ifstream *inputFile = new ifstream(inputFileName);
  while( !(inputFile->eof()) ){
    inputFile->getline(buffer1,200);
    inputFile->getline(buffer2,200);
    if ((!strstr(buffer1,"#") && !(strspn(buffer1," ") == strlen(buffer1))) && 
          (!strstr(buffer2,"#") && !(strspn(buffer2," ") == strlen(buffer2))))
      {
        sscanf(buffer1,"%s",rootFileToChain);
	sscanf(buffer2,"%s",rootFriendFileToChain);
        
	std::cout << "Creating chain ..." << std::endl;
	TChain* chain = new TChain("tree");
	chain->Add(TString(rootFileToChain));

        std::cout << "Adding friend to chain ..." << std::endl;
	TChain* chFriend = new TChain("mjvars/t");
	chain->AddFriend("mjvars/t",TString(rootFriendFileToChain));
	
	if(!chain) {
	  std::cout << "Error: chain not created. End of programme" << std::endl;
	  exit(EXIT_FAILURE);
	}
	nSample++;
	std::cout << "analysing sample n " << nSample << " : " << sampleName[nSample-1]  << std::endl; 
	std::cout<<chain->GetEntries()<<std::endl;      
	//================ Run Analysis
	//zmumujetsAna tree( chain );
	zmumujetsControlSample tree( chain , sampleName[nSample-1].c_str());
	tree.loop(yieldsRow, efficiencyRow); 
	//matrix.push_back(yieldsRow);
	//matrix.push_back(efficiencyRow);
	delete chain;
	delete chFriend;

      }
    	
  }

  if ( yieldsRow.size() != efficiencyRow.size() ) {
    std::cout << "Warning:  different number of steps for yields and efficiencies." << std::endl; 
    std::cout << "Will use the bigger one." << std::endl; 
  }
  Int_t selectionSize = (yieldsRow.size() >= efficiencyRow.size()) ? (yieldsRow.size()/nSample) : (efficiencyRow.size()/nSample);
  FILE* fp;
  const char* finalFileName = "zmumujetsCSandBkg_yieldsEff.dat";
  if ( (fp=fopen(finalFileName,"w")) == NULL) {
    cout<<"Error: '"<<finalFileName<<"' not opened"<<endl;
  } else {
    cout<<"creating file '"<<finalFileName<<"' ..."<<endl;
    fprintf(fp,"#    step         ");
    for(Int_t i = 0; i < nSample; i++) {
      fprintf(fp,"%-16s ",sampleName[i].c_str());
    }
    fprintf(fp,"\n");
    for (Int_t i = 0; i < selectionSize; i++) {
      fprintf(fp,"%-16s",selectionDefinition[i].c_str());
      for(Int_t j = 0; j < nSample; j++) {
	if (yieldsRow.at( i + j * selectionSize) < 1000) fprintf(fp,"%7.2lf ",yieldsRow.at( i + j * selectionSize));	   
	else fprintf(fp,"%7.0lf ",yieldsRow.at( i + j * selectionSize));
	fprintf(fp,"%5.1lf%%   ",(100 * efficiencyRow.at( i + j * selectionSize)));
      } 
      fprintf(fp,"\n");
    }
    fclose(fp);
  }

  inputFile->close();
  delete inputFile;

#endif

#if Application == 2

  std::cout<<"Using Emanuele's trees with skim on mumet > 200"<<std::endl;
  std::cout << "Creating chain ..." << std::endl;
  TChain* chain = new TChain("tree");
  chain->Add("/cmshome/ciprianim/edimarcoTree/tree_metNoMuSkim200/DYJetsToLL_M50/tree.root");
  TChain* chFriend = new TChain("mjvars/t");
  chain->AddFriend("mjvars/t","/cmshome/ciprianim/edimarcoTree/tree_metNoMuSkim200/DYJetsToLL_M50/evVarFriend_DYJetsToLL_M50.root");

  if(!chain) {
    std::cout << "Error: chain not created. End of programme" << std::endl;
    exit(EXIT_FAILURE);
  } 

  std::cout<<chain->GetEntries()<<std::endl;      
  //================ Run Analysis

  zmumujetsAna tree( chain );
  tree.loop();
  
  delete chain;
  delete chFriend;

#endif

#if Application == 3

  std::cout<<"Using Emanuele's trees with skim on 2 leptons"<<std::endl;
  std::cout << "Creating chain ..." << std::endl;
  TChain* chain = new TChain("tree");
  chain->Add("/cmshome/ciprianim/edimarcoTree/tree_2LepGoodSkimVeto/DYJetsToLL_M50/tree.root");
  TChain* chFriend = new TChain("mjvars/t");
  chain->AddFriend("mjvars/t","/cmshome/ciprianim/edimarcoTree/tree_2LepGoodSkimVeto/DYJetsToLL_M50/evVarFriend_DYJetsToLL_M50.root");

  if(!chain) {
    std::cout << "Error: chain not created. End of programme" << std::endl;
    exit(EXIT_FAILURE);
  } 

  std::cout<<chain->GetEntries()<<std::endl;      
  //================ Run Analysis

  zeejetsAna tree( chain );
  tree.loop();
  
  delete chain;
  delete chFriend;

#endif

#if Application == 4 

  if (argc < 2) {
    std::cout << "Not enough arguments: launch as -> "; 
    std::cout << argv[0] << " inputfile " <<std::endl;
    exit(EXIT_FAILURE);
  }
  std::cout<<"Using Emanuele's trees with skim on 2 leptons"<<std::endl;

  char inputFileName[200];
  char buffer1[200];
  char buffer2[200];
  char rootFileToChain[500];
  char rootFriendFileToChain[500];

  //vector< vector<Double_t> > matrix;
  vector< Double_t > yieldsRow;
  vector< Double_t > efficiencyRow;
  Int_t nSample = 0;
  std::vector<std::string> sampleName;
  sampleName.push_back("DYJetsToLL");
  sampleName.push_back("Top");
  sampleName.push_back("QCD");
  //sampleName.push_back("GJets");
  sampleName.push_back("WJetsToLNu");
  sampleName.push_back("ZJetsToNuNu");
  std::vector<std::string> selectionDefinition;
  selectionDefinition.push_back("entry point");
  selectionDefinition.push_back("2 lept SF/OS");
  selectionDefinition.push_back("electrons");
  selectionDefinition.push_back("tight Tag");
  selectionDefinition.push_back("mll");
  selectionDefinition.push_back("elemet > 200");
  selectionDefinition.push_back("njets");
  selectionDefinition.push_back("jet1pt > 110");
  selectionDefinition.push_back("jetjetdphi");
  selectionDefinition.push_back("muon veto");
  selectionDefinition.push_back("photon veto");

  strcpy(inputFileName,argv[1]);

  ifstream *inputFile = new ifstream(inputFileName);
  while( !(inputFile->eof()) ){
    inputFile->getline(buffer1,200);
    inputFile->getline(buffer2,200);
    if ((!strstr(buffer1,"#") && !(strspn(buffer1," ") == strlen(buffer1))) && 
          (!strstr(buffer2,"#") && !(strspn(buffer2," ") == strlen(buffer2))))
      {
        sscanf(buffer1,"%s",rootFileToChain);
	sscanf(buffer2,"%s",rootFriendFileToChain);
        
	std::cout << "Creating chain ..." << std::endl;
	TChain* chain = new TChain("tree");
	chain->Add(TString(rootFileToChain));

        std::cout << "Adding friend to chain ..." << std::endl;
	TChain* chFriend = new TChain("mjvars/t");
	chain->AddFriend("mjvars/t",TString(rootFriendFileToChain));
	
	if(!chain) {
	  std::cout << "Error: chain not created. End of programme" << std::endl;
	  exit(EXIT_FAILURE);
	}
	nSample++;
	std::cout << "analysing sample n " << nSample << " : " << sampleName[nSample-1]  << std::endl; 
	std::cout<<chain->GetEntries()<<std::endl;      
	//================ Run Analysis
	//zeejetsAna tree( chain );
	zeejetsControlSample tree( chain , sampleName[nSample-1].c_str());
	tree.loop(yieldsRow, efficiencyRow); 
	//matrix.push_back(yieldsRow);
	//matrix.push_back(efficiencyRow);
	delete chain;
	delete chFriend;

      }
    	
  }

  if ( yieldsRow.size() != efficiencyRow.size() ) {
    std::cout << "Warning:  different number of steps for yields and efficiencies." << std::endl; 
    std::cout << "Will use the bigger one." << std::endl; 
  }
  Int_t selectionSize = (yieldsRow.size() >= efficiencyRow.size()) ? (yieldsRow.size()/nSample) : (efficiencyRow.size()/nSample);
  FILE* fp;
  const char* finalFileName = "zeejetsCSandBkg_yieldsEff.dat";
  if ( (fp=fopen(finalFileName,"w")) == NULL) {
    cout<<"Error: '"<<finalFileName<<"' not opened"<<endl;
  } else {
    cout<<"creating file '"<<finalFileName<<"' ..."<<endl;
    fprintf(fp,"#    step         ");
    for(Int_t i = 0; i < nSample; i++) {
      fprintf(fp,"%-16s ",sampleName[i].c_str());
    }
    fprintf(fp,"\n");
    for (Int_t i = 0; i < selectionSize; i++) {
      fprintf(fp,"%-16s",selectionDefinition[i].c_str());
      for(Int_t j = 0; j < nSample; j++) {
	if (yieldsRow.at( i + j * selectionSize) < 1000) fprintf(fp,"%7.2lf ",yieldsRow.at( i + j * selectionSize));	   
	else fprintf(fp,"%7.0lf ",yieldsRow.at( i + j * selectionSize));
	fprintf(fp,"%5.1lf%%   ",(100 * efficiencyRow.at( i + j * selectionSize)));
      } 
      fprintf(fp,"\n");
    }
    fclose(fp);
  }

  inputFile->close();
  delete inputFile;

#endif

#if Application == 5

  std::cout<<"Using Emanuele's trees with skim on 2 leptons"<<std::endl;
  std::cout << "Creating chain ..." << std::endl;
  TChain* chain = new TChain("tree");
  chain->Add("/cmshome/ciprianim/edimarcoTree/tree_2LepGoodSkimVeto/DYJetsToLL_M50/tree.root");
  TChain* chFriend = new TChain("mjvars/t");
  chain->AddFriend("mjvars/t","/cmshome/ciprianim/edimarcoTree/tree_2LepGoodSkimVeto/DYJetsToLL_M50/evVarFriend_DYJetsToLL_M50.root");

  if(!chain) {
    std::cout << "Error: chain not created. End of programme" << std::endl;
    exit(EXIT_FAILURE);
  } 

  std::cout<<chain->GetEntries()<<std::endl;      
  //================ Run Analysis

  zmumujetsAna_LepSk tree( chain );
  tree.loop();
  
  delete chain;
  delete chFriend;

#endif

#if Application == 6

  if (argc < 2) {
    std::cout << "Not enough arguments: launch as -> "; 
    std::cout << argv[0] << " inputfile " <<std::endl;
    std::cout << "inputfile is a configuration file containing thresholds and other parameters" << std::endl;
    exit(EXIT_FAILURE);
  }

  char configFileName[200];
  std::strcpy(configFileName,argv[1]);

  std::cout<<"Using Emanuele's trees with skim on mumet > 200"<<std::endl;
  std::cout << "Creating chain ..." << std::endl;
  TChain* chain = new TChain("tree");
  chain->Add("/cmshome/ciprianim/edimarcoTree/tree_metNoMuSkim200/ZJetsToNuNu/tree.root");
  TChain* chFriend = new TChain("mjvars/t");
  chain->AddFriend("mjvars/t","/cmshome/ciprianim/edimarcoTree/tree_metNoMuSkim200/ZJetsToNuNu/evVarFriend_ZJetsToNuNu.root");

  if(!chain) {
    std::cout << "Error: chain not created. End of programme" << std::endl;
    exit(EXIT_FAILURE);
  } 

  std::cout<<chain->GetEntries()<<std::endl;      
  //================ Run Analysis

  znunujetsAna tree( chain );
  tree.loop(configFileName);
  
  delete chain;
  delete chFriend;

#endif

#if Application == 7

#if defined MUON

  std::cout<<"Using Emanuele's trees with skim on mumet > 200"<<std::endl;
  std::cout << "Creating chain ..." << std::endl;
  TChain* chain = new TChain("tree");
  chain->Add("/cmshome/ciprianim/edimarcoTree/tree_metNoMuSkim200/DYJetsToLL_M50/tree.root");
  TChain* chFriend = new TChain("mjvars/t");
  chain->AddFriend("mjvars/t","/cmshome/ciprianim/edimarcoTree/tree_metNoMuSkim200/DYJetsToLL_M50/evVarFriend_DYJetsToLL_M50.root");

#elif defined ELECTRON

  std::cout<<"Using Emanuele's trees with skim on 2 leptons"<<std::endl;
  std::cout << "Creating chain ..." << std::endl;
  TChain* chain = new TChain("tree");
  chain->Add("/cmshome/ciprianim/edimarcoTree/tree_2LepGoodSkimVeto/DYJetsToLL_M50/tree.root");
  TChain* chFriend = new TChain("mjvars/t");
  chain->AddFriend("mjvars/t","/cmshome/ciprianim/edimarcoTree/tree_2LepGoodSkimVeto/DYJetsToLL_M50/evVarFriend_DYJetsToLL_M50.root");

#endif

  if(!chain) {
    std::cout << "Error: chain not created. End of programme" << std::endl;
    exit(EXIT_FAILURE);
  } 

  std::cout<<chain->GetEntries()<<std::endl;      
  //================ Run Analysis

  zlljetsAna tree( chain );
  tree.loop();
  
  delete chain;
  delete chFriend;

#endif


#if Application == 8

  if (argc < 2) {
    std::cout << "Not enough arguments: launch as -> "; 
    std::cout << argv[0] << " inputfile " <<std::endl;
    std::cout << "inputfile is a configuration file containing thresholds and other parameters" << std::endl;
    exit(EXIT_FAILURE);
  }

  char configFileName[200];
  std::strcpy(configFileName,argv[1]);

  ifstream inputFile(configFileName);
  Double_t muonOrElectronOrNeutrino_PDGID;

   if (inputFile.is_open()) {

     Double_t value;
     string parameterName;
     while (inputFile >> parameterName >> value) {

       if (parameterName == "LEP_PDG_ID") {

	 muonOrElectronOrNeutrino_PDGID = (Int_t) value;
	 std::cout << "lepton_pdgID = " << value <<std::endl;

	 if (fabs(muonOrElectronOrNeutrino_PDGID) == 13) {
	   std::cout << "Analysis of Z->mumu" << std::endl;
	 } else if (fabs(muonOrElectronOrNeutrino_PDGID) == 11) {
	   std::cout << "Analysis of Z->ee" << std::endl;
	 } else if (fabs(muonOrElectronOrNeutrino_PDGID) == 12 || fabs(muonOrElectronOrNeutrino_PDGID) == 14 || fabs(muonOrElectronOrNeutrino_PDGID) == 16) {
	   std::cout << "Analysis of Z->nunu (any flavour)" << std::endl;
	 }

       }

     }

     inputFile.close();
                                                                                                                         
   } else {

     cout << "Error: could not open file " << configFileName << endl;
     exit(EXIT_FAILURE);

   }

   TChain* chain = NULL;
   TChain* chFriend = NULL;

   if (fabs(muonOrElectronOrNeutrino_PDGID) == 13) {

     std::cout<<"Using Emanuele's trees with skim on mumet > 200"<<std::endl;
     std::cout << "Creating chain ..." << std::endl;
     chain = new TChain("tree");
     chain->Add("/cmshome/ciprianim/edimarcoTree/tree_metNoMuSkim200/DYJetsToLL_M50/tree.root");
     chFriend = new TChain("mjvars/t");
     chain->AddFriend("mjvars/t","/cmshome/ciprianim/edimarcoTree/tree_metNoMuSkim200/DYJetsToLL_M50/evVarFriend_DYJetsToLL_M50.root");

   } else if (fabs(muonOrElectronOrNeutrino_PDGID) == 11) {

     std::cout<<"Using Emanuele's trees with skim on 2 leptons (biased Axe)"<<std::endl;
     std::cout << "Creating chain ..." << std::endl;
     chain = new TChain("tree");
     chain->Add("/cmshome/ciprianim/edimarcoTree/tree_2LepGoodSkimVeto/DYJetsToLL_M50/tree.root");
     chFriend = new TChain("mjvars/t");
     chain->AddFriend("mjvars/t","/cmshome/ciprianim/edimarcoTree/tree_2LepGoodSkimVeto/DYJetsToLL_M50/evVarFriend_DYJetsToLL_M50.root");

   } else if (fabs(muonOrElectronOrNeutrino_PDGID) == 12 || fabs(muonOrElectronOrNeutrino_PDGID) == 14 || fabs(muonOrElectronOrNeutrino_PDGID) == 16) {

     std::cout<<"Using Emanuele's trees with skim on mumet > 200"<<std::endl;
     std::cout << "Creating chain ..." << std::endl;
     chain = new TChain("tree");
     chain->Add("/cmshome/ciprianim/edimarcoTree/tree_metNoMuSkim200/ZJetsToNuNu/tree.root");
     chFriend = new TChain("mjvars/t");
     chain->AddFriend("mjvars/t","/cmshome/ciprianim/edimarcoTree/tree_metNoMuSkim200/ZJetsToNuNu/evVarFriend_ZJetsToNuNu.root");

   }

  if(!chain) {
    std::cout << "Error: chain not created. End of programme" << std::endl;
    exit(EXIT_FAILURE);
  } 

  std::cout<<chain->GetEntries()<<std::endl;      
  //================ Run Analysis
  if (fabs(muonOrElectronOrNeutrino_PDGID) == 13 || fabs(muonOrElectronOrNeutrino_PDGID) == 11) {

    zlljetsAna_new tree( chain );
    tree.loop(configFileName);

  } else if (fabs(muonOrElectronOrNeutrino_PDGID) == 12 || fabs(muonOrElectronOrNeutrino_PDGID) == 14 || fabs(muonOrElectronOrNeutrino_PDGID) == 16) {

    znunujetsAna tree( chain );
    tree.loop(configFileName);

  }

  delete chain;
  delete chFriend;

#endif

#if Application == 9

  if (argc < 2) {
    std::cout << "Not enough arguments: launch as -> "; 
    std::cout << argv[0] << " inputfile " <<std::endl;
    std::cout << "inputfile is a configuration file containing thresholds and other parameters" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (argc == 3 && std::strcmp("l",argv[2])) {

    std::cout << "not valid option " << argv[2] <<std::endl;
    exit(EXIT_FAILURE);

  }

  char configFileName[200];
  std::strcpy(configFileName,argv[1]);

  ifstream inputFile(configFileName);
  Double_t muonOrElectron_PDGID;

   if (inputFile.is_open()) {

     Double_t value;
     string parameterName;
     while (inputFile >> parameterName >> value) {

       if (parameterName == "LEP_PDG_ID") {

	 muonOrElectron_PDGID = (Int_t) value;
	 std::cout << "lepton_pdgID = " << value <<std::endl;

	 if (fabs(muonOrElectron_PDGID) == 13) {
	   std::cout << "Analysis of Z->mumu" << std::endl;
	 } else if (fabs(muonOrElectron_PDGID) == 11) {
	   std::cout << "Analysis of Z->ee" << std::endl;
	 }

       }

     }

     inputFile.close();
                                                                                                                         
   } else {

     cout << "Error: could not open file " << configFileName << endl;
     exit(EXIT_FAILURE);

   }

   std::cout<<"Using Emanuele's trees with no skim"<<std::endl;
   std::cout << "Creating chain ..." << std::endl;
   TChain* chain = new TChain("tree");
   chain->Add("/cmshome/ciprianim/edimarcoTree/tree_noSkim/DYJetsToLL_M50/tree.root");
   TChain* chFriend = new TChain("mjvars/t");
   chain->AddFriend("mjvars/t","/cmshome/ciprianim/edimarcoTree/tree_noSkim/DYJetsToLL_M50/evVarFriend_DYJetsToLL_M50.root");

  if(!chain) {
    std::cout << "Error: chain not created. End of programme" << std::endl;
    exit(EXIT_FAILURE);
  } 

  std::cout<<chain->GetEntries()<<std::endl;      
  //================ Run Analysis

  if (!(argc == 3 && std::strcmp("l",argv[2]))) {

    zlljets_Axe_noSkim_light tree( chain );
    tree.loop(configFileName);

  } else {

    zlljets_Axe_noSkim tree( chain );
    tree.loop(configFileName);

  }
  
  delete chain;
  delete chFriend;

#endif


  return 0;
}

