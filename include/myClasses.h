#ifndef MyClasses_h
#define MyClasses_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>

class mask {
 
 public:
  mask();
  mask(std::string maskName);
  ~mask();
  std::vector<UInt_t> singleMask;
  std::vector<UInt_t> globalMask;
  std::vector<Double_t> nEvents;
  std::vector<Double_t> nEventsErr2;
  std::string getName() const { return name_; }
  Int_t getMaskSize() const { return singleMask.size(); }
  Double_t getEvents2(const Int_t i) const { return nEvents.at(i) * nEvents.at(i); }
  Double_t getEventsErr(const Int_t i) const { return sqrt(nEventsErr2.at(i)); }
  Double_t getEventsErr2(const Int_t i) const { return nEventsErr2.at(i); }
  void setName(std::string name) { name_ = name; } 
  void append(const UInt_t);
  void countEvents(const UInt_t &, const Double_t &);
  Int_t whichStepHas(const UInt_t &) const;

 private:
  std::string name_;
  

};


class cut {

 public:
  cut(const char *cut_name, const char *var_name, const char *condition, const Double_t threshold, const std::string comment);
  cut(const char *cut_name, const char *var_name, const char *condition, const Double_t threshold);
  cut(const char *cut_name, const char *var_name, const char *condition);
  cut(const char *cut_name, const char *var_name);
  cut(const char *cut_name, const char *var_name, const char *condition, const Double_t threshold, const Double_t thr2, const std::string comment);
  cut(const char *cut_name, const char *var_name, const char *condition, const Double_t threshold, const Double_t thr2);
   ~cut();
   // setters
   //void setFlagTrue() { flag_ = true; }
   //void setFlagFalse() { flag_ = false; }
   void setCondition(std::string condition) { cond_ = condition; }
   void setThreshold(Double_t thr) { thr_ = thr; }
   void setThreshold2(Double_t thr) { thr2_ = thr; }
   // getters
   std::string getCutName() const { return cut_; }
   std::string getVarName() const { return var_; }
   std::string getCondition() const { return cond_; }
   std::string getCutDefinition();
   Double_t getThreshold() const { return thr_; }
   Double_t getThreshold2() const { return thr2_; }
   Int_t getId() const { return id_; }
   UInt_t get2ToId() const { return twoToId_; }
   // other member functions
   void print(ostream &, Bool_t) const;
   void printAllInfo(ostream &) const;
   Bool_t isPassed(Double_t input);
   UInt_t addToMask(Double_t input) {return (this->isPassed(input)) ? (this->get2ToId()) : 0; }
   Bool_t isActive() const { return flag_; }

   static Int_t getNCuts() { return nCuts_; }
   static std::vector<cut*> listOfCuts;
   static void printCutFlow(ostream &, const Int_t, const UInt_t *);
   static void printCutFlow(ostream &, const std::vector<UInt_t>);
   static void printCutFlow(ostream &, const mask *);
   static void printCutFlowAndYields(ostream &, const Double_t, const Double_t, const std::vector<Double_t>, const std::vector<UInt_t>);
   static void printCutFlowAndYields(ostream &, const Double_t, const Double_t, const Double_t, const UInt_t);
   static void printCutFlowAndYields(ostream & myOutStream, const Double_t, const Double_t, const mask *);
   static void printActiveCuts(ostream &);
   static void checkMaskLength();
  
 private:
   Bool_t flag_;
   std::string cut_;
   std::string var_;
   std::string cond_;
   //std::string thr_;
   Double_t thr_;
   Double_t thr2_;
   std::string comment_;
   Int_t id_;          //identifies cut in the mask (is 0 for the first cut, 1 for the second, 2 for the third ecc...)
   UInt_t twoToId_;
   static Int_t nCuts_; // total number of variables on which a cut is applied (or now, if I need different thresholds for the same variable, each one counts as a different cut

};



// work in progress

/*
class cutInterval : public cut {

 public:
  cutInterval(Bool_t flag, const char *cut_name, const char *var_name, const char *condition, const Double_t lowThr, const Double_t upThr, const std::string comment);
  cutInterval(Bool_t flag, const char *cut_name, const char *var_name, const char *condition, const Double_t lowThr, const Double_t upThr);

 private:
  Double_t upThr_;

};

*/

class selection {

 public:
  selection();
  selection(const char *selection_name, const char *definition, const std::string comment);
  selection(const char *selection_name, const char *definition);
   ~selection();
   /* void setName(const char* name) { name_ = name; } */
   /* void setDefinition(const char* definition) { definition_ = definition; } */
   void set(const char* name, const char* definition, const char* comment);
   void set(const char* name, const char* definition);
   void setComment(const char* comment) { comment_ = comment; }
   void setFlagTrue() { flag_ = true; }
   void setFlagFalse() { flag_ = false; }
   std::string getName() const { return name_; }
   std::string getDefinition() const {return definition_; }
   Int_t getId() const { return id_; }
   UInt_t get2ToId() const { return twoToId_; }
   // other member functions
   void print(ostream &, Bool_t) const;
   void printAllInfo(ostream &) const;
   Bool_t isPassed(Bool_t input);
   UInt_t addToMask(Bool_t input) {return input ? (this->get2ToId()) : 0; }
   Bool_t isActive() const { return flag_; }

   static Int_t getNSelections() { return nSelections_; }
   static std::vector<selection*> listOfSelections;
   static void printSelectionFlow(ostream &, const mask *);
   static void printSelectionFlowAndYields(ostream & myOutStream, const Double_t, const Double_t, const Double_t, const UInt_t);
   static void printSelectionFlowAndYields(ostream & myOutStream, const Double_t, const Double_t, const mask *);
   static void saveYieldsAndEfficiency(const Double_t, const mask *, std::vector<Double_t>&, std::vector<Double_t>&);
   static void printActiveSelections(ostream &);
   static void checkMaskLength();
  
 private:
   Bool_t flag_;
   std::string name_;
   std::string definition_;
   std::string comment_;
   Int_t id_;          //identifies selection in the mask (is 0 for the first selection, 1 for the second, 2 for the third ecc...)
   UInt_t twoToId_;
   static Int_t nSelections_; // total number of variables on which a selection is applied (or now, if I need different thresholds for the same variable, each one counts as a different selection

};

#endif 
