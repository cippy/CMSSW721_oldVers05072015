//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat May 16 19:50:26 2015 by ROOT version 5.34/23
// from TTree t/t
// found on file: evVarFriend_DYJetsToLL_M50_HT200to400.root
//////////////////////////////////////////////////////////

#ifndef edimarcoTreeFriend_h
#define edimarcoTreeFriend_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class edimarcoTreeFriend {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         nMu10V;
   Float_t         nMu20T;
   Float_t         nEle10V;
   Float_t         nEle20T;
   Float_t         nTau15V;
   Float_t         nGamma15V;
   Float_t         nGamma175T;
   Float_t         dphijj;
   Float_t         jetclean;
   Float_t         weight;

   // List of branches
   TBranch        *b_nMu10V;   //!
   TBranch        *b_nMu20T;   //!
   TBranch        *b_nEle10V;   //!
   TBranch        *b_nEle20T;   //!
   TBranch        *b_nTau15V;   //!
   TBranch        *b_nGamma15V;   //!
   TBranch        *b_nGamma175T;   //!
   TBranch        *b_dphijj;   //!
   TBranch        *b_jetclean;   //!
   TBranch        *b_weight;   //!

   edimarcoTreeFriend(TTree *tree=0);
   virtual ~edimarcoTreeFriend();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef edimarcoTreeFriend_cxx
edimarcoTreeFriend::edimarcoTreeFriend(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("evVarFriend_DYJetsToLL_M50_HT200to400.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("evVarFriend_DYJetsToLL_M50_HT200to400.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("evVarFriend_DYJetsToLL_M50_HT200to400.root:/mjvars");
      dir->GetObject("t",tree);

   }
   Init(tree);
}

edimarcoTreeFriend::~edimarcoTreeFriend()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t edimarcoTreeFriend::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t edimarcoTreeFriend::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void edimarcoTreeFriend::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nMu10V", &nMu10V, &b_nMu10V);
   fChain->SetBranchAddress("nMu20T", &nMu20T, &b_nMu20T);
   fChain->SetBranchAddress("nEle10V", &nEle10V, &b_nEle10V);
   fChain->SetBranchAddress("nEle20T", &nEle20T, &b_nEle20T);
   fChain->SetBranchAddress("nTau15V", &nTau15V, &b_nTau15V);
   fChain->SetBranchAddress("nGamma15V", &nGamma15V, &b_nGamma15V);
   fChain->SetBranchAddress("nGamma175T", &nGamma175T, &b_nGamma175T);
   fChain->SetBranchAddress("dphijj", &dphijj, &b_dphijj);
   fChain->SetBranchAddress("jetclean", &jetclean, &b_jetclean);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   Notify();
}

Bool_t edimarcoTreeFriend::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void edimarcoTreeFriend::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t edimarcoTreeFriend::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef edimarcoTreeFriend_cxx
