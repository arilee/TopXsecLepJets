#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TObject.h"
#include "TStopwatch.h"
#include "TLegend.h"
#include <vector>
#include <string>
#include <fstream>
#include <vector>
#include <sys/stat.h>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <set>


/***********************
 Files and Directories
************************/
TString dirnameIn= "/home/brochero/ttbar/TopCodettbarSingleLepton/TopResults/";

//TString fl  = "hSF-FixCSV_Tree_LepJets_PHYS14-13TeV-PU20bx25_1000pb-1";  
TString fl  = "hSF-QCD-NewSel_Tree_LepJets_PHYS14-13TeV-PU20bx25_1000pb-1";  

/******************
     Samples
******************/
typedef struct {
  
  //[histo][channel];
  
  TH1F     *hist[15][5];
  THStack  *mc[15][5];
  
} Histograms;


Histograms DoubleMuon;
Histograms DoubleElectron;
Histograms MuEG;  
Histograms Data;

Histograms DYJets;
Histograms DY;

Histograms QCD;

Histograms TbarW;
Histograms TW;
Histograms SingleT;

Histograms TTbar;
Histograms TTJets;
Histograms MCatNLO;
Histograms TTJetsTauola;
Histograms Powheg;


Histograms TTbarBkg;
Histograms WJets;
Histograms NONWZ;

Histograms WW;
Histograms WZ;
Histograms ZZ;
Histograms VV;


Histograms ZJets;
Histograms Z;
Histograms ZDY;

Histograms ratio;     //Ratio

Histograms MCStack ; //MC Stack
Histograms MCStackBkg ; 
Histograms MCStackTTJets ; 
Histograms MCStackMCatNLO ; 
Histograms MCStackTTJetsTauola ;
Histograms MCStackPowheg ;

/******************
  Additional SF
*******************/
Float_t SFmu = 1.0;
Float_t SFe  = 1.0;

/******************
  Efficiencies
*****************/
// lepton selection model (ISO, ID, trigger)->09/12/2013
double lept_IDISO  =  0.0; 
double lept_Trigger=  0.0; 
double lept_uncer  =  0.0; 

//float SF_BR_uncer = 1.7/100.;
float SF_BR_uncer = 0.0;

//float XsecTTbar_uncer = 15.0/100.; 
float XsecTTbar_uncer = 0.0; // Measured Xsec 
float XsecTWVV_uncer  = 0.0;

double SF_uncer[5];//btag[bin]

double lumi_error = 0.0;

double theoryQ2_error      =0.0;
double theoryMatching_error=0.0;

double dd_dd_uncer     = 0.0; //DY Data Driven
double dd_tautau_uncer = 0.0;
double fake_uncer      = 0.0;

//Uncertainties->09/12/2013
double LES_uncer= 0.0;
double JES_uncer= 0.0;
double JER_uncer= 0.0;
double PU_uncer = 0.0;
  
