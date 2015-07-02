#ifndef WHICH_APPLICATION_H
#define WHICH_APPLICATION_H
// choose application in main.cc
//
// 1 -> zmumujetsControlSample
// 2 -> zmumujetsAna
// 3 -> zeejetsAna                       using tree_2LepGoodSkimVeto
// 4 -> zeejetsControlSample                using tree_2LepGoodSkimVeto
// 5 -> zmumujetsAna_LepSk                using tree_2LepGoodSkimVeto
// 6 -> znunujetsAna  
// 7 -> zlljetsAna                   2 and 3, depending on LEPTON directives below    
// 8 -> zlljetsAna_new           code modified    
// 9 -> zlljets_Axe_noSkim          trees without skim are different from previous ones          

#define Application 9

// define leptons flavour for zlljetsAna
// 1 -> muons
// 2-> electrons
#define LEPTON 1

#if LEPTON == 1
#define MUON
#elif LEPTON == 2
#define ELECTRON
#endif


#endif  // WHICH_APPLICATION_H
