#ifndef _FUNCTIONS_FOR_ANALYSIS
#define _FUNCTIONS_FOR_ANALYSIS

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "myClasses.h"

void myDashes(ostream &);
void mySpaces(ostream &, Int_t );

Double_t myGetPhi(Double_t, Double_t);

void myFillOverflowBin(TH1D*, const Double_t);

void myFillOverflowBinW(TH1D*, const Double_t, const Double_t);

void myDrawOverflow(const TH1D *, const char*); 
void myDrawOverflow(const TH1D *, const char*, const Int_t); 
TH1D * myOverflowInLastBin(const TH1D *); 
void myOverflowInLastBin2(TH1D *, const TH1D *); 
void myAddOverflowInLastBin(TH1D *);

Double_t myRejectionFactor(Int_t, Int_t);

Double_t myCrystalBall(double* , double* );

Double_t my2sideCrystalBall(double* , double* );

Double_t myResolutionFunction(double* , double* );

Double_t myResolutionFunctionNoB(double* , double* );

void myAddDefaultPackages(FILE *, const char* );
void makeTableTex(FILE *, const Double_t, const Double_t, const mask*, const std::string);
void makeTableTex(FILE *, const Double_t, const Double_t, const mask*);
void makeTableTexNoEff(FILE *, const Double_t, const Double_t, const mask*, const std::string);
void makeTableTexNoEff(FILE *, const Double_t, const Double_t, const mask*);

char myAskToSaveFile(const char*);

Int_t myGetBin(Double_t, Double_t*, Int_t);

void myPrintEventYields(ostream &, const Double_t, const Int_t, const Double_t *);

void myPrintEventYields(ostream &, const Double_t, const std::vector<Double_t>);

Int_t myPartGenAlgo(const Int_t, const Int_t*, const Int_t*, const Int_t, const Int_t);
Int_t myPartGenAlgo(const Int_t, const Int_t*, const Int_t*, const Int_t, const Int_t, Int_t &, Int_t &);
Int_t myPartGenAlgo(const Int_t, const Int_t*, const Int_t*, const Int_t, const Int_t, Int_t &, Int_t &, Int_t &, const Int_t*);

Int_t myGetPairIndexInArray(const Int_t, const Int_t, const Int_t *, Int_t &, Int_t &);
Int_t myGetPartIndex(const Int_t, const Int_t, const Int_t*);

void myPrintYieldsMetBin(const TH1D* histo, const Double_t *metBinEdges, const Int_t nMetBins, const char* YieldsFileName); 
void myPrintYieldsMetBin(const TH1D* histo, const Double_t *metBinEdges, const Int_t nMetBins); 
void myPrintYieldsMetBinInStream(ostream &, const TH1D* histo, const Double_t *metBinEdges, const Int_t nMetBins); 

#endif
