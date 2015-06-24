#ifndef __CINT__

#include<string>
#include<iostream>
#include<sstream>
#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include<set>
#include<vector>

#endif

// Root
#include "TString.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TTree.h"
//#include "TKey.h"
//#include "TPrint.h"
//#include <exception>
#include <sys/stat.h>

#ifndef __CINT__

void display_usage()
{
  std::cout << "\033[1;37musage:\033[1;m skimfile cutindex [options]" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "    -i inputfile  Input file without .root" << std::endl;
  std::cout << "    -o name in the output file \"h_\"" << std::endl;
  std::cout << "    -d Input file directory. Default directory: InputTrees" << std::endl;
  std::cout << "    -h                 displays this help message and exits " << std::endl;
  std::cout << "" << std::endl;
}


// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const TString currentDateTime() {
  time_t     now = time(0);
  struct tm  tstruct;
  char       buf[80];
  tstruct = *localtime(&now);
  // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
  // for more information about date/time format
  strftime(buf, sizeof(buf), "%Y-%m-%d at %X", &tstruct);

  return buf;
}

int main(int argc, const char* argv[]){

  gSystem->Load("libTree");

  bool   _syst      = false;
  bool	 _tr_unc    = false;
  bool	 _idiso_unc = false;
  const char * _output   = 0;
  const char * _input    = 0;
  const char * _dir      = "/home/brochero/ttbar/TopTrees_CATuples/";
  const char * _tr       = 0;
  const char * _idiso    = 0;

  // Arguments used
  //std::set<int> usedargs;
  //Parsing input options
  if(argc == 1){
    display_usage();
    return -1;
  }

  else{
      //Argumet 1 must be a valid input fileName
      for (int i = 1; i < argc; i++){
	if( strcmp(argv[i],"-i") == 0 ){
	  _input = argv[i+1];
	  i++;
	}
	if( strcmp(argv[i],"-d") == 0 ){
	  _dir = argv[i+1];
	  i++;
	}
	if( strcmp(argv[i],"-o") == 0 ){
	  _output= argv[i+1];
	  i++;
	}
	if( strcmp(argv[i],"-h") == 0 ||
	    strcmp(argv[i],"--help") == 0 ){
	  display_usage();
	  return 0;
	}
      }
  }//else
  if( _input ==0 ){
    std::cerr << "\033[1;31mskimfile ERROR:\033[1;m The '-i' option is mandatory!"
	      << std::endl;
    display_usage();
    return -1;
  }
  
  // reassigning
  TString fname(_input);
  TString hname(_output);
  TString fdir(_dir);

  
  TChain theTree("ttbarSingleLepton/AnalysisTree"); 
  
  std::cout << "-----------------------------------------" << std::endl;
  std::cout << "-----------------------------------------" << std::endl;
  std::cout << "Signal: ";
  std::cout << fname + ".root" << std::endl;
  
  theTree.Add(fdir + fname + ".root");

  
  int Event,Channel,PV,GoodPV;
  float PUWeight; // Temporal
  float MET,MET_Phi;

  float Lep_px, Lep_py, Lep_pz, Lep_E;
  std::vector<float> *Jet_px=0, *Jet_py=0, *Jet_pz=0, *Jet_E=0;
  std::vector<bool>  *Jet_LooseID=0;
  std::vector<float> *Jet_CSV=0;
  
  theTree.SetBranchAddress( "event",    &Event );
  theTree.SetBranchAddress( "PUWeight", &PUWeight );
  theTree.SetBranchAddress( "channel",  &Channel );

  theTree.SetBranchAddress( "PV",     &PV );
  theTree.SetBranchAddress( "GoodPV",  &GoodPV );

  theTree.SetBranchAddress( "MET",     &MET );
  theTree.SetBranchAddress( "MET_phi", &MET_Phi );

  theTree.SetBranchAddress( "lepton_px", &Lep_px );
  theTree.SetBranchAddress( "lepton_py", &Lep_py );
  theTree.SetBranchAddress( "lepton_pz", &Lep_pz );
  theTree.SetBranchAddress( "lepton_E",  &Lep_E );

  theTree.SetBranchAddress( "jet_px", &Jet_px );
  theTree.SetBranchAddress( "jet_py", &Jet_py );
  theTree.SetBranchAddress( "jet_pz", &Jet_pz );
  theTree.SetBranchAddress( "jet_E",  &Jet_E );
  theTree.SetBranchAddress( "jet_LooseID",  &Jet_LooseID );

  theTree.SetBranchAddress( "jet_CSV",  &Jet_CSV );

  
  /*********************************
             Histograms
  **********************************/
  
  TH1F *hPV[4][2];
  TH1F *hMET[4][2],*hMET_Phi[4][2];
  TH1F *hLepPt[4][2],*hLepEta[4][2],*hLepPhi[4][2];
  TH1F *hNJets[4][2],*hHT[4][2],*hNBtagJets[4][2];
  TH1F *CSV[4][4][2];

  TString namech[2];
  namech[0]="mujets";
  namech[1]="ejets";

  TString namecut[4];
  namecut[0]="lepton";
  namecut[1]="4Jets";
  namecut[2]="MET";
  namecut[3]="2btag";

  TString titlenamech[2];
  titlenamech[0]="#mu+Jets";
  titlenamech[1]="e+Jets";
  
  for(int j=0; j<4; j++){   // Cut
    for(int i=0; i<2; i++){ // Channel
      hPV[j][i]         = new TH1F("hPV_"+namech[i]+"_"+namecut[j],"PV Distribution  " + titlenamech[i],30,0,30);
      
      hMET[j][i]        = new TH1F("hMET_"+namech[i]+"_"+namecut[j],"#slash{E}_{T} " + titlenamech[i],40,0,200);
      hMET_Phi[j][i]    = new TH1F("hMET_Phi_"+namech[i]+"_"+namecut[j],"#Phi_{#slash{E}_{T}} " + titlenamech[i],160,-4,4);
      
      hLepPt [j][i]    = new TH1F("hLepPt_"  +namech[i] + "_" + namecut[j], "Lepton p_{T} " + titlenamech[i],50,0.0,250.0);
      hLepEta[j][i]    = new TH1F("hLepEta_" +namech[i] + "_" + namecut[j], "#eta_{Lep} " + titlenamech[i],50,-2.5,2.5);
      hLepPhi[j][i]    = new TH1F("hLepPhi_" +namech[i] + "_" + namecut[j], "#phi_{Lep} " + titlenamech[i],100,-5,5);
      
      hNJets[j][i]      = new TH1F("hNJets_"+namech[i]+"_"+namecut[j],"Jet multiplicity " + titlenamech[i],9,-0.5,8.5);

      hNBtagJets[j][i]  = new TH1F("hNBtagJets_"+namech[i]+"_"+namecut[j],"b-tag jet multiplicity " + titlenamech[i],9,-0.5,8.5);

      hHT[j][i]         = new TH1F("hHT_"+namech[i]+"_"+namecut[j],"H_{T} " + titlenamech[i],300,0,600);

      TString jetn[4];
      jetn[0]= "Jet0"; 
      jetn[1]= "Jet1"; 
      jetn[2]= "Jet2"; 
      jetn[3]= "Jet3"; 

      for(int ij=0; ij<4; ij++){
	CSV[ij][j][i]         = new TH1F("hCSV_" + jetn[ij] + "_" + namech[i] + "_" + namecut[j],"CSV " + jetn[ij] + " " + titlenamech[i],80,0,1);
      }
    }//for(i)
  }//for(j)


  TStopwatch sw;
  sw.Start(kTRUE);

  // Number de events for <pT Reweight>
  //          [Cut][Channel]
  float SF_pTweight[4][2]={0,0,0,0,
			   0,0,0,0};

  // Number de events for acceptance
  //          [Cut][Channel]
  int AccEvent[4][2]={0,0,0,0,
		      0,0,0,0};
  // Number de events for acceptance
  //          [Cut][Channel]
  float EffEvent[4][2]={0.0,0.0,0.0,0.0,
			0.0,0.0,0.0,0.0};
  


  std::cout << "--- Processing: " << theTree.GetEntries() << " events" << std::endl;
  for (Long64_t ievt=0; ievt<theTree.GetEntries();ievt++) {
     
    theTree.GetEntry( ievt ); 
    //////////////////////////////////////////////////////
    int step = theTree.GetEntries()/50;
    if (ievt%(step) == 0){
      float progress=(ievt)/(theTree.GetEntries()*1.0);
      int barWidth = 50;
      
      std::cout << "[";
      int pos = barWidth * progress;
    
      for (int i = 0; i < barWidth; ++i) {
      	if (i < pos) std::cout << "=";
      	else if (i == pos) std::cout << ">";
      	else std::cout << " ";
      }
      
      std::cout << "] " << int(progress * 100.0) << " %\r";
      std::cout.flush();
    }
    ////////////////////////////////////////////////////////

    int NJets,NBtagJets;
    
    TLorentzVector Lep;
    std::vector<TLorentzVector> Jet;
    std::vector<TLorentzVector> bJet;
         
    Lep.SetPxPyPzE(Lep_px,Lep_py,Lep_pz,Lep_E);
    if(Lep.Pt()<30) continue; // Lep pT >30GeV

    // Jets 
    NJets     = 0;
    NBtagJets = 0;

    for(int ijet=0; ijet < Jet_px->size(); ijet++){

      TLorentzVector jet;
      jet.SetPxPyPzE((*Jet_px)[ijet],(*Jet_py)[ijet],(*Jet_pz)[ijet],(*Jet_E)[ijet]);
      
      if(jet.Pt()>25 &&         // Jet pT > 25GeV
	 !(*Jet_LooseID)[ijet]){ // Loose ID
	
	Jet.push_back(jet);
	NJets++; // Number of jets
	
	if((*Jet_CSV)[ijet] > 0.814){ // CSVM. Luca b-tagging code to apply SF?
	  bJet.push_back(jet);
	  NBtagJets++; // Number of b-tagged jets
	} // if(b-tag)
      } // if(Jet_pT)
    }// for(jets)

          
    /***************************
            Selection
    ***************************/
    int                                    cut=0; // Single Lepton
    if(NJets>3)                            cut=1; // + 4Jets 
    if(NJets>3 && MET>30.0)                cut=2; // + MET
    if(NJets>3 && MET>30.0 && NBtagJets>1) cut=3; // + 2 Btag


   // /*******************
   //    Fill Histograms
   //  ******************/
     
    float PUWeight_event=PUWeight;

    for(int icut=0; icut<cut+1; icut++){

      // PUWeight reset for each cut
      PUWeight=PUWeight_event;
  
      /************************************************************
       Scale Factors (Luminosity)
      *************************************************************/
      // Number of events from https://cmsweb.cern.ch/das
      // Cross Sections at 13 TeV:
      // - ttbar from https://twiki.cern.ch/twiki/bin/view/Sandbox/FullNNLOcrossSections 
      //              https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
      // - backgrounds from https://twiki.cern.ch/twiki/bin/view/CMS/StandardModelCrossSectionsat13TeV

      int Lumi=1000; //pb
      // PUWeight = PUWeight*Lumi*(1.0/N_Gen_events)*(Xsec)*(Br)
      if(fname.Contains("QCD"))   PUWeight = PUWeight * Lumi * (1.0/4777926.0)  * (866600000.0) * (0.00044);     // [pb] (cross section) * (Filter Eff)
      if(fname.Contains("ZJets")) PUWeight = PUWeight * Lumi * (1.0/2829164.0)  * (3591.6) * (0.03*3);           // [pb]
      if(fname.Contains("WJets")) PUWeight = PUWeight * Lumi * (1.0/10017930)   * (61526.7) * (0.1*3);           // [pb]
      if(fname.Contains("tW"))    PUWeight = PUWeight * Lumi * (1.0/986100.0)   * (35.6);                        // [pb]
      if(fname.Contains("tbarW")) PUWeight = PUWeight * Lumi * (1.0/971800.0)   * (35.6);                        // [pb]
      if(fname.Contains("ttbar")) PUWeight = PUWeight * Lumi * (1.0/25446993.0) * (827.1) * (0.108*3) * (0.67);  // [pb]

      /*************************************************************/

      /******************
          Acceptace
      ******************/
      AccEvent[icut][Channel]++;
      EffEvent[icut][Channel] = EffEvent[icut][Channel] + PUWeight;


      hPV[icut][Channel]->Fill(PV,PUWeight);
      
      hMET[icut][Channel]->Fill(MET,PUWeight);
      hMET_Phi[icut][Channel]->Fill(MET_Phi,PUWeight);
      
      hLepPt[icut][Channel]->Fill(Lep.Pt(),PUWeight);
      hLepEta[icut][Channel]->Fill(Lep.Eta(),PUWeight);
      hLepPhi[icut][Channel]->Fill(Lep.Phi(),PUWeight);
      
      hNJets[icut][Channel]->Fill(NJets,PUWeight); 
      hNBtagJets[icut][Channel]->Fill(NBtagJets,PUWeight);

    for(int ijet=0; ijet < Jet_px->size(); ijet++){
      if (ijet<4) CSV[ijet][icut][Channel]->Fill((*Jet_CSV)[ijet],PUWeight);
    }

    }//for(icuts) 
    
  Jet.clear();    
  bJet.clear();    

  }//for(events)


  delete Jet_px;
  delete Jet_py; 
  delete Jet_pz;
  delete Jet_E;
  delete Jet_CSV;

  // Get elapsed time
  sw.Stop();
  std::cout << "==================================================] 100% " << std::endl;
  std::cout << "--- End of event loop: "; sw.Print();
  

  //Acceptance-Efficiency
  std::cout << "--------  Acceptace  --------" << std::endl;
  std::cout << "Number of RAW-mu+Jets events:" << std::endl;
  std::cout << "lepton: "   << AccEvent[0][0] << std::endl;
  std::cout << "4 Jets: "   << AccEvent[1][0] << std::endl;
  std::cout << "MET: "      << AccEvent[2][0] << std::endl;
  std::cout << "2 btag: "   << AccEvent[3][0] << std::endl;

  std::cout << "--------  Efficiency  --------" << std::endl;
  std::cout << "Number of Weigthed-mu+Jets events:" << std::endl;
  std::cout << "lepton: "   << EffEvent[0][0] << " +/- " << sqrt(AccEvent[0][0])*EffEvent[0][0]/AccEvent[0][0] << std::endl;
  std::cout << "4 Jets: "   << EffEvent[1][0] << " +/- " << sqrt(AccEvent[1][0])*EffEvent[1][0]/AccEvent[1][0] << std::endl;
  std::cout << "MET: "      << EffEvent[2][0] << " +/- " << sqrt(AccEvent[2][0])*EffEvent[2][0]/AccEvent[2][0] << std::endl;
  std::cout << "2 btag: "   << EffEvent[3][0] << " +/- " << sqrt(AccEvent[3][0])*EffEvent[3][0]/AccEvent[3][0] << std::endl;
  
  std::cout << "--------  Acceptace  --------" << std::endl;
  std::cout << "Number of RAW-e+Jets events:" << std::endl;
  std::cout << "lepton: "   << AccEvent[0][1] << std::endl;
  std::cout << "4 Jets: "   << AccEvent[1][1] << std::endl;
  std::cout << "MET: "      << AccEvent[2][1] << std::endl;
  std::cout << "2 btag: "   << AccEvent[3][1] << std::endl;

  std::cout << "--------  Efficiency  --------" << std::endl;
  std::cout << "Number of Weigthed-e+Jets events: " << std::endl;
  std::cout << "lepton: "   << EffEvent[0][1] << " +/- " << sqrt(AccEvent[0][1])*EffEvent[0][1]/AccEvent[0][1] << std::endl;
  std::cout << "4 Jets: "   << EffEvent[1][1] << " +/- " << sqrt(AccEvent[1][1])*EffEvent[1][1]/AccEvent[1][1] << std::endl;
  std::cout << "MET: "      << EffEvent[2][1] << " +/- " << sqrt(AccEvent[2][1])*EffEvent[2][1]/AccEvent[2][1] << std::endl;
  std::cout << "2 btag: "   << EffEvent[3][1] << " +/- " << sqrt(AccEvent[3][1])*EffEvent[3][1]/AccEvent[3][1] << std::endl;
  
  //Output Dir
  TString dirname="TopResults";   
  // make a dir if it does not exist!!
  struct stat st;
  if(stat(dirname,&st) != 0) system("mkdir " + dirname);
  
  TString systname="";
  bool matchname=false;
  
  TString temp_fname=fname + ".root";
  for(int i=0;i<temp_fname.Sizeof();i++){
    if (i>2){
      if (temp_fname[i-3]=='-' && 
  	  temp_fname[i-2]=='1' && 
  	  temp_fname[i-1]=='_') matchname=true;
    }
    if (temp_fname[i]=='.') matchname=false;
    if (matchname) systname.Append(temp_fname[i]);
  }
  
  // Yields
  TString Yieldfile=dirname + "/Yields_" + hname + ".h";
  FILE* fyields = fopen(Yieldfile, "a");
  
  fprintf(fyields,"// %s Sample on %s \n", (fname + ".root").Data() , currentDateTime().Data());
  fprintf(fyields,"// %s version \n", hname.Data());
  fprintf(fyields,"float  %s[2][4];//[channel][Cut] \n", systname.Data());
  fprintf(fyields,"float  err_%s[2][4]; //[channel][Cut] \n", systname.Data());
  fprintf(fyields,"// Channel: [0]=mu+Jets [1]=e+Jets \n");
  fprintf(fyields,"// Cut [0]=lepton [1]= Jets [2]=MET [3]=btag \n");
  for(int ch=0;ch<2;ch++){
    for(int cut=0;cut<4;cut++){
      fprintf(fyields,"%s[%i][%i] = %.3f ; \n", systname.Data(), ch, cut, EffEvent[cut][ch]);
      if(AccEvent[cut][ch]!=0.0) fprintf(fyields,"err_%s[%i][%i] = %.3f ; \n", systname.Data(), ch, cut, sqrt(AccEvent[cut][ch])*EffEvent[cut][ch]/AccEvent[cut][ch]);
      else fprintf(fyields,"err_%s[%i][%i] = 0.0 ; \n", systname.Data(), ch, cut);
    }
  }
  fclose(fyields);
  
  std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  std::cout << "Yields saved into " << Yieldfile << " file" << std::endl;
  std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;    
  
  // --- Write histograms
  
  TString outfname=dirname + "/hSF-" + hname + "_" + fname + ".root";
  TFile *target  = new TFile(outfname,"RECREATE" );  

  for(int j=0; j<4; j++){
    for(int i=0; i<2; i++){
      
      hPV[j][i]->Write();
    
      hMET[j][i]->Write();
      hMET_Phi[j][i]->Write();
      
      hLepPt[j][i]->Write();
      hLepEta[j][i]->Write();
      hLepPhi[j][i]->Write();

      hNJets[j][i]->Write();
      hNBtagJets[j][i]->Write();            

      for(int ij=0; ij<4; ij++){
        CSV[ij][j][i]->Write();
      }

    }//for(i)

  }//for(j)
  
  std::cout << "File saved as " << outfname << std::endl;

}


#endif

