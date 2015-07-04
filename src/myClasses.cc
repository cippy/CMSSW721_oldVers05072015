#include "myClasses.h"

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

// #include <TAxis.h>
// #include <TChain.h>
// #include <TCanvas.h>
// #include <TChain.h>
// #include <TFile.h>
// #include <TH1F.h>
// #include <TH2.h>
#include <TMath.h>
#include <TROOT.h>
// #include <TStyle.h>
// #include <TTree.h>
// #include <TVector2.h>
// #include <TVector3.h>

using namespace std;


mask::mask() {

  name_ = "";

}

mask::mask(string name) {

  name_ = name;

}

mask::~mask() {

  //cout<<"~mask() called for "<<name_<<endl;

}

void mask::append(const UInt_t x) {

  this->singleMask.push_back(x);
  if (this->singleMask.size() != 1) this->globalMask.push_back(this->singleMask.back() + this->globalMask.back());
  else  this->globalMask.push_back(this->singleMask.back());
  this->nEvents.push_back(0.0);
  this->nEventsErr2.push_back(0.0);

}


// void mask::countEvents(const UInt_t &eventMask, const Double_t &weight) {
   
//   for (Int_t i = 0; i < this->getMaskSize(); i++) {
//     if ( (eventMask & this->globalMask[i]) == this->globalMask[i]) {                //  exploiting bit-bit &
//       this->nEvents[i] += weight;  
//       this->nEventsErr2[i] += weight*weight;      
//     }
//   }

// }


void mask::countEvents(const UInt_t &eventMask) {
   
  //globalMask includes the previous step: once the eventMask doesn't match the globalMask of i-th step, it's useless to continue
  Int_t i = 0;
  Int_t size = this->getMaskSize();
  
  while ( ((eventMask & this->globalMask[i]) == this->globalMask[i]) && i < size) {                //  exploiting bit-bit &
    this->nEvents[i] += 1;  
    this->nEventsErr2[i] += 1;
    i++;      
  }

}

void mask::countEvents(const UInt_t &eventMask, const Double_t &weight) {
   
  //globalMask includes the previous step: once the eventMask doesn't match the globalMask of i-th step, it's useless to continue
  Int_t i = 0;
  Int_t size = this->getMaskSize();
  
  while ( ((eventMask & this->globalMask[i]) == this->globalMask[i]) && i < size) {                //  exploiting bit-bit 
    this->nEvents[i] += weight;  
    this->nEventsErr2[i] += weight*weight;
    i++;      
  }

}

Int_t mask::whichStepHas(const UInt_t &a) const {

  Int_t size = this->getMaskSize();
  Int_t index = size;
  Int_t i = 0;

  while( (index == size) && (i < size)) {
    if( (this->singleMask[i] & a) == a ) index = i;
    i++; 
  }

  if (index != size) return index;   
  else {
    cout << " Error: step not found in the mask " << endl;
    exit(EXIT_FAILURE);
  }

}

 //----------------------------------------------------------------------------------

// implementation of methods of class cut 

 //static data members
 Int_t cut::nCuts_ = 0;
 vector<cut*> cut::listOfCuts; 

 void cut::printCutFlow(ostream & myOutStream, const Int_t cutSteps, const UInt_t *singleCutMask) {

   //this function prints the cut flow on the appropriate ofstream (can be a file or cout). Since it's a member function of class cut, it only needs to get the number of cut steps. For now it needs the array of masks which is not a class data member yet

   myOutStream<<"**************************"<<endl;
   myOutStream<<"*          CUTS FLOW          *"<<endl;
   myOutStream<<"**************************"<<endl;
   myOutStream<<"-----------------------------------------------------------------------------------"<<endl;  
   myOutStream<<"Printing list of cuts applied at each step"<<endl;
   for (Int_t i = 0; i < cutSteps; i++) {
     myOutStream<<"-----------------------------------"<<endl;
     myOutStream<<setw(2)<<(i+1)<<endl;
     for (Int_t j = 0; j < cut::getNCuts(); j++) {
       if ((singleCutMask[i] >> j) & 1) {
 	myOutStream<<"      ";
 	cut::listOfCuts[j]->print(myOutStream,1); 
       }
     }
   }
   myOutStream<<"-----------------------------------"<<endl;  

 }

 void cut::printCutFlow(ostream & myOutStream, const vector<UInt_t> singleCutMask) {

   // this function prints the cut flow on the appropriate ofstream (can be a file or cout). Since it's a member function of class cut, it only needs to get the number of cut steps. For now it needs the array of masks which is not a class data member yet

  myOutStream<<"**************************"<<endl;
  myOutStream<<"*          CUTS FLOW          *"<<endl;
  myOutStream<<"**************************"<<endl;
  myOutStream<<"-----------------------------------------------------------------------------------"<<endl;  
  myOutStream<<"Printing list of cuts applied at each step"<<endl;
  for (Int_t i = 0; i < singleCutMask.size(); i++) {
    myOutStream<<"-----------------------------------"<<endl;
    myOutStream<<setw(2)<<(i+1)<<endl;
    for (Int_t j = 0; j < cut::getNCuts(); j++) {
      if ((singleCutMask[i] >> j) & 1) {
	myOutStream<<"      ";
	cut::listOfCuts[j]->print(myOutStream,1); 
      }
    }
  }
  myOutStream<<"-----------------------------------"<<endl;  

}

void cut::printCutFlow(ostream & myOutStream, const mask *m) {

  //this function prints the cut flow on the appropriate ofstream (can be a file or cout). Since it's a member function of class cut, it only needs to get the number of cut steps. For now it needs the array of masks which is not a class data member yet

  myOutStream<<"**************************"<<endl;
  myOutStream<<"*          CUTS FLOW          *"<<endl;
  myOutStream<<"**************************"<<endl;
  myOutStream<<"Mask's name --> "<<m->getName()<<endl;
  myOutStream<<"-----------------------------------------------------------------------------------"<<endl;  
  myOutStream<<"Printing list of cuts applied at each step"<<endl;
  for (Int_t i = 0; i < m->getMaskSize(); i++) {
    myOutStream<<"-----------------------------------"<<endl;
    myOutStream<<setw(2)<<(i+1)<<endl;
    for (Int_t j = 0; j < cut::getNCuts(); j++) {
      if ((m->singleMask[i] >> j) & 1) {
	myOutStream<<"      ";
	cut::listOfCuts[j]->print(myOutStream,1); 
      }
    }
  }
  myOutStream<<"-----------------------------------"<<endl;  

}

void cut::printCutFlowAndYields(ostream & myOutStream, const Double_t lumi, const Double_t nwentries, const vector<Double_t> eventsInStep, const vector<UInt_t> singleCutMask) {

  myOutStream<<"**************************"<<endl;
  myOutStream<<"*        CUTS FLOW & YIELDS        *"<<endl;
  myOutStream<<"**************************"<<endl;
  myOutStream<<"-----------------------------------------------------------------------------------"<<endl;  
  //nwentries is the number of events taking all weights into account (cuts due to triggers were not taken into account)
  myOutStream<<"nwentries:  weighted total number of entries = "<<nwentries<<endl;
  myOutStream<<"n:               number of events after i-th cut"<<endl;
  myOutStream<<"aR:            absolute ratio = n(i)/nwentries"<<endl;
  myOutStream<<"rR:             relative ratio = n(i)/n(i-1)"<<endl;
  myOutStream<<"**************************"<<endl;
  myOutStream<<"data normalised to "<<lumi<<" fb^-1"<<endl;
  myOutStream<<"==========================================================="<<endl;  
  myOutStream<<setw(5)<<right<<"cut"<<setw(25)<<"definition"<<setw(12)<<"n"<<setw(12)<<"aR"<<setw(8)<<"rR"<<endl;    
  myOutStream<<"==========================================================="<<endl;  
  for (Int_t i = 0; i < eventsInStep.size(); i++) {
    if (i == 0) { 
      myOutStream<<setw(5)<<right<<i+1<<setw(25)<<""<<setw(12)<<eventsInStep[i]<<fixed<<setprecision(4)<<setw(12)<<eventsInStep[i]/nwentries<<
      setw(8)<<eventsInStep[i]/nwentries<<endl;  
    } else {
      myOutStream<<setw(5)<<right<<i+1<<setw(25)<<""<<setw(12)<<eventsInStep[i]<<fixed<<setprecision(4)<<setw(12)<<eventsInStep[i]/nwentries<<
      setw(8)<<eventsInStep[i]/eventsInStep[i-1]<<endl;
    }
    for (Int_t j = 0; j < cut::getNCuts(); j++) {
      if ((singleCutMask[i] >> j) & 1) {
	// flag for comment is 0 (not enough room to write comments during this phase)
	myOutStream<<setw(5)<<" "<<setw(25)<<right<<cut::listOfCuts[j]->getCutDefinition()<<endl;;  
      }
    }
    myOutStream<<"-----------------------------------------------------------------------------------"<<endl;  
  }
  myOutStream<<"-----------------------------------------------------------------------------------"<<endl;  
  myOutStream<<endl;

}

void cut::printCutFlowAndYields(ostream & myOutStream, const Double_t lumi, const Double_t nwentries, const Double_t eventsInStep, const UInt_t singleCutMask) {

  myOutStream<<"**************************"<<endl;
  myOutStream<<"*        CUTS FLOW & YIELDS        *"<<endl;
  myOutStream<<"**************************"<<endl;
  myOutStream<<"-----------------------------------------------------------------------------------"<<endl;  
  //nwentries is the number of events taking all weights into account (cuts due to triggers were not taken into account)
  myOutStream<<"nwentries:  weighted total number of entries = "<<nwentries<<endl;
  myOutStream<<"n:               number of events after cut"<<endl;
  myOutStream<<"aR:            absolute ratio = n/nwentries"<<endl;
  myOutStream<<"**************************"<<endl;
  myOutStream<<"data normalised to "<<lumi<<" fb^-1"<<endl;
  myOutStream<<"==========================================================="<<endl;  
  myOutStream<<setw(5)<<right<<"cut"<<setw(25)<<"definition"<<setw(12)<<"n"<<setw(12)<<"aR"<<endl;    
  myOutStream<<"==========================================================="<<endl;   
  myOutStream<<setw(5)<<right<<1<<setw(25)<<""<<setw(12)<<eventsInStep<<fixed<<setprecision(4)<<setw(12)<<eventsInStep/nwentries<<endl;  
  for (Int_t j = 0; j < cut::getNCuts(); j++) {
    if ((singleCutMask >> j) & 1) {
      myOutStream<<setw(5)<<" "<<setw(25)<<cut::listOfCuts[j]->getCutDefinition()<<endl;
    } 
  }
  myOutStream<<"-----------------------------------------------------------------------------------"<<endl;  
  myOutStream<<endl;

}

void cut::printCutFlowAndYields(ostream & myOutStream, const Double_t lumi, const Double_t nwentries, const mask *m) {

  myOutStream<<"**************************"<<endl;
  myOutStream<<"*        CUTS FLOW & YIELDS        *"<<endl;
  myOutStream<<"**************************"<<endl;
  myOutStream<<"Mask's name --> "<<m->getName()<<endl;
  myOutStream<<"-----------------------------------------------------------------------------------"<<endl;  
  //nwentries is the number of events taking all weights into account (cuts due to triggers were not taken into account)
  myOutStream<<"nwentries:  weighted total number of entries = "<<nwentries<<endl;
  myOutStream<<"n:               number of events after i-th cut"<<endl;
  myOutStream<<"aR:            absolute ratio = n(i)/nwentries"<<endl;
  myOutStream<<"rR:             relative ratio = n(i)/n(i-1)"<<endl;
  myOutStream<<"**************************"<<endl;
  myOutStream<<"data normalised to "<<lumi<<" fb^-1"<<endl;
  myOutStream<<"==========================================================="<<endl;  
  myOutStream<<setw(5)<<right<<"cut"<<setw(25)<<"definition"<<setw(12)<<"n"<<setw(12)<<"aR"<<setw(8)<<"rR"<<endl;    
  myOutStream<<"==========================================================="<<endl;  
  for (Int_t i = 0; i < m->getMaskSize(); i++) {
    if (i == 0) { 
      myOutStream<<setw(5)<<right<<i+1<<setw(25)<<""<<setw(12)<<m->nEvents[i]<<fixed<<setprecision(4)<<setw(12)<<m->nEvents[i]/nwentries<<
      setw(8)<<m->nEvents[i]/nwentries<<endl;  
    } else {
      myOutStream<<setw(5)<<right<<i+1<<setw(25)<<""<<setw(12)<<m->nEvents[i]<<fixed<<setprecision(4)<<setw(12)<<m->nEvents[i]/nwentries<<
      setw(8)<<m->nEvents[i]/m->nEvents[i-1]<<endl;
    }
    for (Int_t j = 0; j < cut::getNCuts(); j++) {
      if ((m->singleMask[i] >> j) & 1) {
	// flag for comment is 0 (not enough room to write comments during this phase)
	myOutStream<<setw(5)<<" "<<setw(25)<<right<<cut::listOfCuts[j]->getCutDefinition()<<endl;;  
      }
    }
    myOutStream<<"-----------------------------------------------------------------------------------"<<endl;  
  }
  myOutStream<<"-----------------------------------------------------------------------------------"<<endl;  
  myOutStream<<endl;

}

void cut::printActiveCuts(ostream & myOutStream) {

  myOutStream<<"------------------------------------------------------------------------------------------"<<endl;
  myOutStream<<"Printing list of activated cuts"<<endl;
  myOutStream<<"------------------------------------------------------------------------------------------"<<endl;
  for (Int_t i = 0; i < cut::getNCuts(); i++ ) {
    if ( cut::listOfCuts[i]->isActive() ) {
      cut::listOfCuts[i]->printAllInfo(myOutStream);
    }
  }
  myOutStream<<"------------------------------------------------------------------------------------------"<<endl;

}

void cut::checkMaskLength() {

  if (cut::getNCuts() > 8*sizeof(UInt_t)) {
     cout<<"Warning: not enough bits in the mask to accomodate all "<<cut::getNCuts()<<" cuts (max is "<<8*sizeof(UInt_t)<<").\n End of programme."<<endl;
     exit(EXIT_FAILURE);
   }

}

cut::cut(const char *cut_name, const char *var_name, const char *condition, const Double_t threshold, const string comment) {
  flag_ = true;
  cut_ =cut_name;
  var_ = var_name;
  cond_ = condition;
  thr_ = threshold;
  thr2_ = 0.;
  comment_ = comment;
  id_ = nCuts_;
  twoToId_ = (UInt_t) TMath::Power(2.,id_);
  nCuts_++;
  listOfCuts.push_back(this);
}

cut::cut(const char *cut_name, const char *var_name, const char *condition, const Double_t threshold) {
  flag_ = true;
  cut_ =cut_name;
  var_ = var_name;
  cond_ = condition;
  thr_ = threshold;
  thr2_ =  0.;
  comment_ = "";
  id_ = nCuts_;
  twoToId_ = (UInt_t) TMath::Power(2.,id_);
  nCuts_++;
  listOfCuts.push_back(this);
}

cut::cut(const char *cut_name, const char *var_name, const char *condition) {
  flag_ =false;
  cut_ =cut_name;
  var_ = var_name;
  cond_ = condition;
  thr_ = 0.;
  thr2_ = 0.;
  comment_ = "";
  id_ = nCuts_;
  twoToId_ = (UInt_t) TMath::Power(2.,id_);
  nCuts_++;
  listOfCuts.push_back(this);
}

cut::cut(const char *cut_name, const char *var_name) {
  flag_ =false;
  cut_ =cut_name;
  var_ = var_name;
  cond_ = " ";
  thr_ = 0.;
  thr2_ =  0.;
  comment_ = "";
  id_ = nCuts_;
  twoToId_ = (UInt_t) TMath::Power(2.,id_);
  nCuts_++;
  listOfCuts.push_back(this);
}

cut::cut(const char *cut_name, const char *var_name, const char *condition, const Double_t threshold, const Double_t thr2, const string comment) {
  flag_ = true;
  cut_ =cut_name;
  var_ = var_name;
  cond_ = condition;
  thr_ = (threshold <= thr2) ? threshold : thr2;
  thr2_ = (threshold >= thr2) ? threshold : thr2;
  comment_ = comment;
  id_ = nCuts_;
  twoToId_ = (UInt_t) TMath::Power(2.,id_);
  nCuts_++;
  listOfCuts.push_back(this);
}

cut::cut(const char *cut_name, const char *var_name, const char *condition, const Double_t threshold, const Double_t thr2) {
  flag_ = true;
  cut_ =cut_name;
  var_ = var_name;
  cond_ = condition;
  thr_ = (threshold <= thr2) ? threshold : thr2;
  thr2_ = (threshold >= thr2) ? threshold : thr2;
  comment_ = "";
  id_ = nCuts_;
  twoToId_ = (UInt_t) TMath::Power(2.,id_);
  nCuts_++;
  listOfCuts.push_back(this);
}

cut::~cut() {
  //cout<<"~cut() called for "<<cut_<<endl;
  nCuts_--;
}

string cut::getCutDefinition() { 

  // following if is to have, as an example, 20.0 written as 20
  if (cond_ != "€") {
    if (fmod(thr_,1.) == 0) return (var_ + " " + cond_ + " " + Form("%4i",(Int_t)thr_)); 
    else return (var_ + " " + cond_ + " " + Form("%4.1lf",thr_)); 
  } else {
    if ((fmod(thr_,1.) == 0) && (fmod(thr2_,1.) == 0)) return (var_ + " " + cond_ + " " + Form("[%i, %i]",(Int_t)thr_,(Int_t)thr2_)); 
    else return (var_ + " " + cond_ + " " + Form("[%lf,%lf]",thr_,thr2_)); 
  }

}

void cut::printAllInfo(ostream & myOutStream) const {

  if (cond_ != "€") {
    myOutStream<<setw(20)<<left<<cut_<<": id = "<<setw(2)<<right<<id_<<" | "<<setw(18)<<left<<var_<<" "<<setw(2)<<cond_<<" ";
    myOutStream<<setw(4)<<right<<thr_<<"   "<<left<<comment_<<endl; 
  } else {
    myOutStream<<setw(20)<<left<<cut_<<": id = "<<setw(2)<<right<<id_<<" | "<<setw(18)<<left<<var_<<" "<<setw(2)<<">"<<" ";
    myOutStream<<setw(4)<<right<<thr_<<"   "<<left<<comment_<<endl;   
    myOutStream<<setw(32)<<" "<<setw(18)<<left<<var_<<" "<<setw(2)<<"<"<<" ";
    myOutStream<<setw(4)<<right<<thr2_<<"   "<<left<<comment_<<endl;
  }

 }

void cut::print(ostream & myOutStream = cout, Bool_t addComment = 0) const {

  //print cut definition on the right ofstream (cout is default), if addComment = 1, cut  comment is printed if any
  if (addComment) {
    if (cond_ != "€") {
      myOutStream<<setw(18)<<left<<var_<<" "<<setw(2)<<cond_<<" "<<setw(4)<<right<<thr_<<"   "<<left<<comment_<<endl; 
    } else {
      myOutStream<<setw(18)<<left<<var_<<" "<<setw(2)<<cond_<<" "<<setw(4)<<right<<Form("[%4.1lf, %4.1lf]",thr_,thr2_)<<"   "<<left<<comment_<<endl;
    }
  } else {
    if (cond_ != "in") {
      myOutStream<<setw(18)<<left<<var_<<" "<<setw(2)<<cond_<<" "<<setw(4)<<right<<thr_<<endl; 
    } else {
      myOutStream<<setw(18)<<left<<var_<<" "<<setw(2)<<cond_<<" "<<setw(4)<<right<<Form("[%4.1lf, %4.1lf]",thr_,thr2_)<<endl;
    }
  }
  
}

Bool_t cut::isPassed(Double_t input) {

  if (cond_ == "<") return (input < thr_) ? true : false;
  else if (cond_== ">") return (input > thr_) ? true : false;
  else if (cond_== "=") return (input == thr_) ? true : false;
  else if (cond_ == "<=") return (input <= thr_) ? true : false;
  else if (cond_ == ">=") return (input >= thr_) ? true : false;
  else if (cond_ == "€") return ((input > thr_) && (input < thr2_)) ? true : false;  
  else {
    cout<<"Error in cut::isPassed()! No condition fulfilled."<<endl;
    cout<<"End of programme"<<endl;
    exit(EXIT_FAILURE);
  }

}



//----------------------------------------------------------------------------------

// implementation of methods of class selection


//static data members
 Int_t selection::nSelections_ = 0;
 vector<selection*> selection::listOfSelections; 

void selection::printSelectionFlow(ostream & myOutStream, const mask *m) {

  //this function prints the selection flow on the appropriate ofstream (can be a file or cout). Since it's a member function of class selection, it only needs to get the number of selection steps. For now it needs the array of masks which is not a class data member yet

  myOutStream<<"**************************"<<endl;
  myOutStream<<"*          SELECTIONS FLOW          *"<<endl;
  myOutStream<<"**************************"<<endl;
  myOutStream<<"Mask's name --> "<<m->getName()<<endl;
  myOutStream<<"-----------------------------------------------------------------------------------"<<endl;  
  myOutStream<<"Printing list of selections applied at each step"<<endl;
  for (Int_t i = 0; i < m->getMaskSize(); i++) {
    myOutStream<<"-----------------------------------"<<endl;
    myOutStream<<setw(2)<<(i+1)<<endl;
    for (Int_t j = 0; j < selection::getNSelections(); j++) {
      if ((m->singleMask[i] >> j) & 1) {
	myOutStream<<"      ";
	selection::listOfSelections[j]->print(myOutStream,1); 
      }
    }
  }
  myOutStream<<"-----------------------------------"<<endl;  

}

void selection::printSelectionFlowAndYields(ostream & myOutStream, const Double_t lumi, const Double_t nwentries, const Double_t eventsInStep, const UInt_t singleCutMask) {

  myOutStream<<"**************************"<<endl;
  myOutStream<<"*        SELECTIONS FLOW & YIELDS        *"<<endl;
  myOutStream<<"**************************"<<endl;
  myOutStream<<"-----------------------------------------------------------------------------------"<<endl;  
  //nwentries is the number of events taking all weights into account (cuts due to triggers were not taken into account)
  myOutStream<<"nwentries:  weighted total number of entries = "<<nwentries<<endl;
  myOutStream<<"n:               number of events after selection"<<endl;
  myOutStream<<"aR:            absolute ratio = n/nwentries"<<endl;
  myOutStream<<"**************************"<<endl;
  myOutStream<<"data normalised to "<<lumi<<" fb^-1"<<endl;
  myOutStream<<"==========================================================="<<endl;  
  myOutStream<<setw(5)<<right<<"step"<<setw(25)<<"definition"<<setw(12)<<"n"<<setw(12)<<"aR"<<endl;    
  myOutStream<<"==========================================================="<<endl;   
  myOutStream<<setw(5)<<right<<1<<setw(25)<<""<<setw(12)<<eventsInStep<<fixed<<setprecision(4)<<setw(12)<<eventsInStep/nwentries<<endl;  
  for (Int_t j = 0; j < selection::getNSelections(); j++) {
    if ((singleCutMask >> j) & 1) {
      myOutStream<<setw(5)<<" "<<setw(25)<<selection::listOfSelections[j]->getDefinition()<<endl;
    } 
  }
  myOutStream<<"-----------------------------------------------------------------------------------"<<endl;  
  myOutStream<<endl;

}

void selection::printSelectionFlowAndYields(ostream & myOutStream, const Double_t lumi, const Double_t nwentries, const mask *m) {

  myOutStream<<"**************************"<<endl;
  myOutStream<<"*        SELECTIONS FLOW & YIELDS        *"<<endl;
  myOutStream<<"**************************"<<endl;
  myOutStream<<"Mask's name --> "<<m->getName()<<endl;
  myOutStream<<"-----------------------------------------------------------------------------------"<<endl;  
  //nwentries is the number of events taking all weights into account (selections due to triggers were not taken into account)
  myOutStream<<"nwentries:  weighted total number of entries = "<<nwentries<<endl;
  myOutStream<<"n:               number of events after i-th selection"<<endl;
  myOutStream<<"aR:            absolute ratio = n(i)/nwentries"<<endl;
  myOutStream<<"rR:             relative ratio = n(i)/n(i-1)"<<endl;
  myOutStream<<"**************************"<<endl;
  myOutStream<<"data normalised to "<<lumi<<" fb^-1"<<endl;
  myOutStream<<"==========================================================="<<endl;  
  myOutStream<<setw(5)<<right<<"step"<<setw(25)<<"definition"<<setw(12)<<"n"<<setw(12)<<"aR"<<setw(8)<<"rR"<<endl;    
  myOutStream<<"==========================================================="<<endl;  
  for (Int_t i = 0; i < m->getMaskSize(); i++) {
    if (i == 0) { 
      myOutStream<<setw(5)<<right<<i+1<<setw(25)<<""<<setw(12)<<m->nEvents[i]<<fixed<<setprecision(4)<<setw(12)<<m->nEvents[i]/nwentries<<
      setw(8)<<m->nEvents[i]/nwentries<<endl;  
    } else {
      // avoid division by 0
      if( m->nEvents[i-1] == 0 ) {
	myOutStream<<setw(5)<<right<<i+1<<setw(25)<<""<<setw(12)<<m->nEvents[i]<<fixed<<setprecision(4)<<setw(12)<<m->nEvents[i]/nwentries<<
	  setw(8)<<1.000<<endl;
      } else {
	myOutStream<<setw(5)<<right<<i+1<<setw(25)<<""<<setw(12)<<m->nEvents[i]<<fixed<<setprecision(4)<<setw(12)<<m->nEvents[i]/nwentries<<
	  setw(8)<<m->nEvents[i]/m->nEvents[i-1]<<endl;
      }
    }
    for (Int_t j = 0; j < selection::getNSelections(); j++) {
      if ((m->singleMask[i] >> j) & 1) {
	myOutStream<<setw(5)<<" "<<setw(25)<<right<<selection::listOfSelections[j]->getDefinition()<<endl;;  
      }
    }
    myOutStream<<"-----------------------------------------------------------------------------------"<<endl;  
  }
  myOutStream<<"-----------------------------------------------------------------------------------"<<endl;  
  myOutStream<<endl;

}

void selection::saveYieldsAndEfficiency(const Double_t nTotal, const mask *m, vector<Double_t> &yRow, vector<Double_t> &eRow) {

  yRow.push_back(nTotal);
  eRow.push_back(1.0000);
  for(Int_t i = 0; i < m->getMaskSize(); i++) {
    yRow.push_back(m->nEvents[i]);
    if (i == 0) eRow.push_back(m->nEvents[i]/nTotal);
    else if( (i != 0) && (m->nEvents[i-1] == 0) ) eRow.push_back(1.0000);
    else eRow.push_back(m->nEvents[i]/m->nEvents[i-1]);
  }
  
}

void selection::printActiveSelections(ostream & myOutStream) {

  myOutStream<<"------------------------------------------------------------------------------------------"<<endl;
  myOutStream<<"Printing list of activated selections"<<endl;
  myOutStream<<"------------------------------------------------------------------------------------------"<<endl;
  for (Int_t i = 0; i < selection::getNSelections(); i++ ) {
    if ( selection::listOfSelections[i]->isActive() ) {
      selection::listOfSelections[i]->printAllInfo(myOutStream);
    }
  }
  myOutStream<<"------------------------------------------------------------------------------------------"<<endl;

}

void selection::checkMaskLength() {

  if (selection::getNSelections() > 8*sizeof(UInt_t)) {
    cout<<"Warning: not enough bits in the mask to accomodate all "<<selection::getNSelections()<<" selections (max is "<<8*sizeof(UInt_t)<<")"<<endl;
    cout<<"End of programme."<<endl;
    exit(EXIT_FAILURE);
  }

}

selection::selection() {

  flag_ = false;
  name_ = "";
  definition_ = "";
  comment_ = "";
  id_ = nSelections_;
  twoToId_ = (UInt_t) TMath::Power(2.,id_);
  nSelections_++;
  listOfSelections.push_back(this);

}

selection::selection(const char *selection_name, const char *definition, const string comment) {

  flag_ = true;
  name_ = selection_name;
  definition_ = definition;
  comment_ = comment;
  id_ = nSelections_;
  twoToId_ = (UInt_t) TMath::Power(2.,id_);
  nSelections_++;
  listOfSelections.push_back(this);

  checkMaskLength();

}

selection::selection(const char *selection_name, const char *definition) {

  flag_ = true;
  name_ = selection_name;
  definition_ = definition;
  comment_ = "";
  id_ = nSelections_;
  twoToId_ = (UInt_t) TMath::Power(2.,id_);
  nSelections_++;
  listOfSelections.push_back(this);

  checkMaskLength();

}

selection::~selection() {
  //cout<<"~selection() called for "<<name_<<endl;
  nSelections_--;
}

void selection::set(const char* name, const char* definition, const char* comment) {

  flag_ = true;
  name_ = name;
  definition_ = definition;
  comment_ = comment;
  // id_ = nSelections_;
  // twoToId_ = (UInt_t) TMath::Power(2.,id_);
  // nSelections_++;
  // listOfSelections.push_back(this);
 
  checkMaskLength();

}

void selection::set(const char* name, const char* definition) {

  flag_ = true;
  name_ = name;
  definition_ = definition;
  // id_ = nSelections_;
  // twoToId_ = (UInt_t) TMath::Power(2.,id_);
  // nSelections_++;
  // listOfSelections.push_back(this);

  checkMaskLength();

}

void selection::printAllInfo(ostream & myOutStream) const {

  myOutStream<<setw(20)<<left<<name_<<": id = "<<setw(2)<<right<<id_<<" | "<<setw(25)<<definition_<<"   "<<left<<comment_<<endl;   
 
 }

void selection::print(ostream & myOutStream = cout, Bool_t addComment = 0) const {

  //print selection definition on the right ofstream (cout is default), if addComment = 1, selection  comment is printed if any
  if (addComment) {
    myOutStream<<setw(25)<<left<<definition_<<left<<"   "<<comment_<<endl; 
  } else {
    myOutStream<<setw(25)<<left<<definition_<<endl;
  }
  
}

Bool_t selection::isPassed(Bool_t input) {

  return input;

}

