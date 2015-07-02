#ifndef CUTS_THRESHOLDS
#define CUTS_THRESHOLDS

#include "whichApplication.h"

#define LUMI 5         // integrated luminosity in fb^-1 to which yields are normalized
#define RATIO_BR_ZINV_ZLL 5.942
#define UNC_RATIO_BR_ZINV_ZLL 0.019
#define TAU_VETO_FLAG false      // decide whether to include the tau veto in the selection: true is yes, false is no

//thesholds for cuts and other parameters
#define NJETS 2
#define J1PT 110
#define J1ETA 2.5
#define J2PT 30
#define J2ETA 2.5
#define J1J2DPHI 2.5

#if defined MUON

#define ROOT_FNAME "zmumujetsAna.root"
#define TXT_FNAME "zmumujetsAna.txt"
#define TEX_FNAME "zmumujetsAna.tex"

#define FLAVOUR "mu"
#define LL_FLAVOUR "mumu"
#define CONTROL_SAMPLE "Z-->mumu"
#define LEP_PDG_ID 13
#define NEG_LEP_PDG_ID2 -169   // -13*13, used to see if I have 2 OS muons
#define LEP1PT 10
#define LEP2PT 10
#define LEP1ETA 2.4
#define LEP2ETA 2.4
#define DILEPMASS_LOW 60
#define DILEPMASS_UP 120
#define LEP_ISO_04 0.12
/* #define GENLEP1PT LEP1PT */
/* #define GENLEP2PT LEP2PT */
/* #define GENLEP1ETA LEP1ETA */
/* #define GENLEP2ETA LEP2ETA */
/* #define GEN_ZMASS_LOW DILEPMASS_LOW */
/* #define GEN_ZMASS_UP DILEPMASS_UP */
#define GENLEP1PT 10
#define GENLEP2PT 10
#define GENLEP1ETA 2.4
#define GENLEP2ETA 2.4
#define GEN_ZMASS_LOW 60
#define GEN_ZMASS_UP 120

#elif defined ELECTRON

#define ROOT_FNAME "zeejetsAna.root"
#define TXT_FNAME "zeejetsAna.txt"
#define TEX_FNAME "zeejetsAna.tex"

#define HLT_ELECTRON true
//HLT parameter
#define HLT_LEP1PT 30
#define HLT_LEP2PT 18
#define HLT_LEP1ETA 2.5
#define HLT_LEP2ETA 2.5

#define FLAVOUR "ele"
#define LL_FLAVOUR "ee"
#define CONTROL_SAMPLE "Z-->ee"
#define LEP_PDG_ID 11
#define NEG_LEP_PDG_ID2 -121   // -11*11, used to see if I have 2 OS electrons
#define LEP1PT 32
#define LEP2PT 20
#define LEP1ETA 2.4
#define LEP2ETA 2.4
#define DILEPMASS_LOW 60
#define DILEPMASS_UP 120
#define LEP_ISO_04 0.12
#define GENLEP1PT LEP1PT
#define GENLEP2PT LEP2PT
#define GENLEP1ETA LEP1ETA
#define GENLEP2ETA LEP2ETA
#define GEN_ZMASS_LOW DILEPMASS_LOW
#define GEN_ZMASS_UP DILEPMASS_UP

#endif

#endif


