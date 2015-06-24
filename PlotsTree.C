#include "General.h"


void PlotsTree(TString plots="1btag", int ch=0, TString sys="", bool EvtNorm=false){

  
  TString files  = dirnameIn + fl;  
  
  /****************
      Channel
   ***************/
  
  TString channel;
  if(ch==0) channel="mujets";
  if(ch==1) channel="ejets";  
  if(ch==2) channel="lepjet"; 
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1000);
  gStyle->SetOptStat("emruo");
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);


  /****************
      tdr Style
   ***************/
  gROOT->ProcessLine(".L /home/brochero/ttbar/TopCodettbarSingleLepton/tdrStyle.C");
  setTDRStyle();

  Int_t chatch = 1756;
  TColor *color = new TColor(chatch, 0.3, 0.5, 0.5, "", 0.45); // alpha = 0.5
  

  if(plots=="dilepton") dd_dd_uncer=0.0;

  int histos=8;

  /****************
       DY - DD
   ***************/  

  ZDY=loadhistograms(plots,files + "_ZJets.root");
   
  scalehistograms(ZDY, 0, SFmu); // Scale Factors
  scalehistograms(ZDY, 1, SFe);   // Scale Factors

  ZDY=addhistograms(ZDY); // All Channels
  Z=ZDY;

  setupdraw(Z, kAzure-2);



  /****************
     MC Sample
  ***************/  

  TW=loadhistograms(plots,files + "_tW" + sys + ".root");

  scalehistograms(TW, 0, SFmu); // Scale Factors
  scalehistograms(TW, 1, SFe);   // Scale Factors

  TW=addhistograms(TW); // All Channels


  // TbarW=loadhistograms(plots,files + "_tbarW" + sys + ".root");

  // scalehistograms(TbarW, 0, SFmu); // Scale Factors
  // scalehistograms(TbarW, 1, SFe);   // Scale Factors

  // TbarW=addhistograms(TbarW); // All Channels

  for(int j=0; j<histos; j++){
    for(int k=0; k<3; k++){
      
      SingleT.hist[j][k]=(TH1F*)TW.hist[j][k]->Clone();
      //SingleT.hist[j][k]->Add(SingleT.hist[j][k],TbarW.hist[j][k]);
      
    }
  }
  

  setupdraw(SingleT, kPink-3); //Select Histogram Color

  ///////////////////////////////////////////////////////////////

  TTbar=loadhistograms(plots,files + "_ttbar" + sys + ".root"); // 7TeV
  //TTbar=loadhistograms(plots,files + "_TTJets_MadSpin" + sys + ".root"); // 8TeV

  scalehistograms(TTbar, 0, SFmu);
  scalehistograms(TTbar, 1, SFe);
  
  TTbar=addhistograms(TTbar);
  setupdraw(TTbar, kRed+1);

  //////////////////////////////////////////////////////////////
  WJets=loadhistograms(plots,files + "_WJets.root");

  scalehistograms(WJets, 0, SFmu); // Scale Factors
  scalehistograms(WJets, 1, SFe);   // Scale Factors

  WJets=addhistograms(WJets); // All Channels

  for(int j=0; j<histos; j++){
    for(int k=0; k<3; k++){

      NONWZ.hist[j][k]=(TH1F*)WJets.hist[j][k]->Clone();// Only W+Jets

    }
  }  
  setupdraw(NONWZ, kGreen-3); //Select Histogram Color
  
  //////////////////////////////////////////////////////////////
  QCD=loadhistograms(plots,files + "_QCD.root");

  scalehistograms(QCD, 0, SFmu);  // Scale Factors
  scalehistograms(QCD, 1, SFe);   // Scale Factors

  QCD=addhistograms(QCD); // All Channels

  setupdraw(QCD, kOrange-3); //Select Histogram Color
  

  ////////////////////////////////////////////////////////////////////
  /////////////////////     THStack     //////////////////////////////
  ////////////////////////////////////////////////////////////////////

  for(int j=0; j<histos; j++){
    
    TString variable;

    for(int k=0; k<3; k++){
      MCStack.mc[j][k]=new THStack(variable, "");
      MCStack.mc[j][k]->SetHistogram(TTbar.hist[j][k]);
    }
  }

  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  MCStack=addstack(MCStack,TTbar);

  MCStack=addstack(MCStack,Z);

  MCStack=addstack(MCStack,SingleT);

  MCStack=addstack(MCStack,NONWZ);

  MCStack=addstack(MCStack,QCD);

  ///////////////////////////////////////////////////////////////////

  TCanvas *cPlots[8];//histos
  
  for(int i=0;i<histos;i++){
    char can[200];
    sprintf(can,"canvas%i",i);
    cPlots[i]=new TCanvas(can,"Preselection Plots");
    cPlots[i]->Divide(1,2);  
  }

  
  TLegend *leg[8];
  TPad    *pad[8][2]; 
  
  for(int i=0; i<histos; i++){

    // Temporal
    Data.hist[i][ch] = (TH1F*) MCStack.mc[i][ch]->GetStack()->Last()->Clone();


    leg[i]=createlegend(leg[i],plots);//each Plot
    
    //Plot Pad
    pad[i][0] = (TPad*)cPlots[i]->GetPad(1); 
    pad[i][0]->SetPad(0.01, 0.23, 0.99, 0.99);
    pad[i][0]->SetTopMargin(0.1);
    pad[i][0]->SetRightMargin(0.04);
    
    //Ratio Pad
    pad[i][1] = (TPad*)cPlots[i]->GetPad(2);
    pad[i][1]->SetPad(0.01, 0.02, 0.99, 0.3);
    pad[i][1]->SetGridx();
    pad[i][1]->SetGridy();
    pad[i][1]->SetTopMargin(0.05);
    pad[i][1]->SetBottomMargin(0.4);
    pad[i][1]->SetRightMargin(0.04);

    pad[i][0]->cd();
 
    MCStack.mc[i][ch]->Draw("hist");
    MCStack.mc[i][ch]->GetYaxis()->SetTitle("Events");
    MCStack.mc[i][ch]->GetYaxis()->SetTitleOffset(1.2);
    MCStack.mc[i][ch]->GetYaxis()->SetTitleSize(0.07);
    MCStack.mc[i][ch]->GetYaxis()->SetLabelSize(0.055);
    MCStack.mc[i][ch]->GetYaxis()->SetNdivisions(607);
    //MCStack.mc[i][ch]->GetYaxis()->SetLabelSize(0.05);
    TGaxis *hYaxis = (TGaxis*)MCStack.mc[i][ch]->GetYaxis();
    //hYaxis->SetMaxDigits(3);
    MCStack.mc[i][ch]->GetXaxis()->SetLabelSize(0.0);
    MCStack.mc[i][ch]->GetXaxis()->SetTitle("");

    
    // Number of events in the Legend
    // leg[i]=addlegend(leg[i],"Data      ", Data.hist[i][ch]);
    // leg[i]=addlegend(leg[i],"VV        ", VV.hist[i][ch]);
    // leg[i]=addlegend(leg[i],"Non W/Z   ", NONWZ.hist[i][ch]);
    // leg[i]=addlegend(leg[i],"Single Top", SingleT.hist[i][ch]);
    // leg[i]=addlegend(leg[i],"Z/#gamma* l^{+}l^{-}  ", Z.hist[i][ch]);
    // leg[i]=addlegend(leg[i],"t#bar{t}        ", TTbar.hist[i][ch]);
    // leg[i]->Draw("SAME");


    
    float maxh=Data.hist[i][ch]->GetMaximum();
    if(maxh<MCStack.mc[i][ch]->GetMaximum()) maxh=MCStack.mc[i][ch]->GetMaximum();
    MCStack.mc[i][ch]->SetMaximum(1.7*maxh);
    
    /******************************************
     ******************************************
                  Uncertainties
    *******************************************
    *******************************************/

    TH1F *totalmc = (TH1F*) MCStack.mc[i][ch]->GetStack()->Last()->Clone();

    for (int nbin=0; nbin<totalmc->GetNbinsX(); nbin++){


      /***************************
       Theoretical  Uncertainties
      ***************************/

      float ThTTbarUnc=sqrt((TTbar.hist[i][ch]->GetBinContent(nbin+1)*theoryQ2_error)*
			    (TTbar.hist[i][ch]->GetBinContent(nbin+1)*theoryQ2_error)+
			    (TTbar.hist[i][ch]->GetBinContent(nbin+1)*theoryMatching_error)*
			    (TTbar.hist[i][ch]->GetBinContent(nbin+1)*theoryMatching_error)			   
			    );
      
      float ThSingleTUnc=sqrt((SingleT.hist[i][ch]->GetBinContent(nbin+1)*theoryQ2_error)*
			      (SingleT.hist[i][ch]->GetBinContent(nbin+1)*theoryQ2_error)+
			      (SingleT.hist[i][ch]->GetBinContent(nbin+1)*theoryMatching_error)*
			      (SingleT.hist[i][ch]->GetBinContent(nbin+1)*theoryMatching_error)			   
			      );
      
      float ThVVUnc=0.0;

      float ThTotalUnc=sqrt(ThTTbarUnc*ThTTbarUnc + ThSingleTUnc*ThSingleTUnc + ThVVUnc*ThVVUnc);
      
      /***************************
           Other Uncertainties
      ***************************/
    
      float LumUnc=totalmc->GetBinContent(nbin+1)*lumi_error;

      float XsecTTbarUnc=sqrt((TTbar.hist[i][ch]->GetBinContent(nbin+1)*SF_BR_uncer)*(TTbar.hist[i][ch]->GetBinContent(nbin+1)*SF_BR_uncer) + (TTbar.hist[i][ch]->GetBinContent(nbin+1)*XsecTTbar_uncer )*(TTbar.hist[i][ch]->GetBinContent(nbin+1)*XsecTTbar_uncer));
      float XsecSingleTUnc=SingleT.hist[i][ch]->GetBinContent(nbin+1)*XsecTWVV_uncer;
      float XsecVVUnc=0.0*XsecTWVV_uncer;
      float XsecTotalUnc=sqrt(XsecTTbarUnc*XsecTTbarUnc + XsecSingleTUnc*XsecSingleTUnc + XsecVVUnc*XsecVVUnc);
		      
      float PUUnc=totalmc->GetBinContent(nbin+1)*PU_uncer;

      float SFIDISOUnc=totalmc->GetBinContent(nbin+1)*lept_uncer;

      float DYDDUnc= Z.hist[i][ch]->GetBinContent(nbin+1)*dd_dd_uncer;
		        
      float btagUnc;
      
      btagUnc=0.0;      

      float EnergyUnc=sqrt((totalmc->GetBinContent(nbin+1)*LES_uncer)*
			   (totalmc->GetBinContent(nbin+1)*LES_uncer) +
			   (totalmc->GetBinContent(nbin+1)*JES_uncer)*
			   (totalmc->GetBinContent(nbin+1)*JES_uncer) +
			   (totalmc->GetBinContent(nbin+1)*JER_uncer)*
			   (totalmc->GetBinContent(nbin+1)*JER_uncer)
			   );

      float PlotUncer=0.0;
      // if (plots=="dilepton") PlotUncer = sqrt(totalmc->GetBinError(nbin+1)*totalmc->GetBinError(nbin+1) + 
      // 					      LumUnc*LumUnc + 
      // 					      PUUnc*PUUnc + 
      // 					      SFIDISOUnc*SFIDISOUnc + 
      // 					      XsecTotalUnc*XsecTotalUnc );
      // if (plots=="2Jets")    PlotUncer = sqrt(totalmc->GetBinError(nbin+1)*totalmc->GetBinError(nbin+1) + 
      // 					      LumUnc*LumUnc + 
      // 					      PUUnc*PUUnc + 
      // 					      SFIDISOUnc*SFIDISOUnc + 
      // 					      XsecTotalUnc*XsecTotalUnc + 
      // 					      DYDDUnc*DYDDUnc);
      // if (plots=="MET")      PlotUncer = sqrt(totalmc->GetBinError(nbin+1)*totalmc->GetBinError(nbin+1) + 
      // 					      LumUnc*LumUnc + 
      // 					      PUUnc*PUUnc + 
      // 					      SFIDISOUnc*SFIDISOUnc + 
      // 					      XsecTotalUnc*XsecTotalUnc + 
      // 					      DYDDUnc*DYDDUnc);
      // if (plots=="1btag" || (plots=="MET" && i==6))    PlotUncer = sqrt(totalmc->GetBinError(nbin+1)*totalmc->GetBinError(nbin+1) + 
      // 									LumUnc*LumUnc + 
      // 									PUUnc*PUUnc + 
      // 									SFIDISOUnc*SFIDISOUnc + 
      // 									XsecTotalUnc*XsecTotalUnc + 
      // 									DYDDUnc*DYDDUnc + 
      // 									EnergyUnc*EnergyUnc + 
      // 									btagUnc*btagUnc);

      if (plots=="lepton")   PlotUncer=totalmc->GetBinError(nbin+1);
      if (plots=="4Jets")    PlotUncer=totalmc->GetBinError(nbin+1);
      if (plots=="MET")      PlotUncer=totalmc->GetBinError(nbin+1);
      if (plots=="2btag")    PlotUncer=totalmc->GetBinError(nbin+1);
      
      totalmc->SetBinError(nbin+1,PlotUncer);
    }

    TGraphErrors *thegraph = new TGraphErrors(totalmc);
    thegraph->SetName("thegraph");
    // thegraph->SetFillStyle(3004);
    // thegraph->SetFillColor(1);
    thegraph->SetFillStyle(1001);
    thegraph->SetFillColor(chatch);
    thegraph->SetLineColor(chatch);

    thegraph->Draw("e2SAME");


    Data.hist[i][ch]->Sumw2();
    Data.hist[i][ch]->SetMarkerStyle(20);
    float MarkerDataSize=0.7;
    Data.hist[i][ch]->SetMarkerSize(MarkerDataSize);
    Data.hist[i][ch]->Draw("SAME");

    /***********************
             Ratio
     **********************/    

    //Graph Ratio Clone
    TH1F *Ra;
    Ra=(TH1F*)Data.hist[i][ch]->Clone();
    Ra->Divide(totalmc);
    ratio.hist[i][ch]=Ra;
    
    ratio.hist[i][ch]->SetMarkerStyle(20);
    ratio.hist[i][ch]->SetMarkerSize(MarkerDataSize);
    ratio.hist[i][ch]->SetMarkerColor(1);
    ratio.hist[i][ch]->SetLineColor(1);
    ratio.hist[i][ch]->SetLineWidth(1);
    ratio.hist[i][ch]->SetMaximum(2);
    ratio.hist[i][ch]->SetMinimum(0);
    ratio.hist[i][ch]->SetTitle("");

    ratio.hist[i][ch]->GetYaxis()->SetTitle("Obs/Exp");
    ratio.hist[i][ch]->GetYaxis()->CenterTitle();
    ratio.hist[i][ch]->GetYaxis()->SetTitleOffset(0.45);
    ratio.hist[i][ch]->GetYaxis()->SetTitleSize(0.16);
    ratio.hist[i][ch]->GetYaxis()->SetLabelSize(0.15);
    ratio.hist[i][ch]->GetYaxis()->SetNdivisions(402);
    ratio.hist[i][ch]->GetXaxis()->SetNdivisions(509);
    ratio.hist[i][ch]->GetXaxis()->SetTitleOffset(1.1);
    ratio.hist[i][ch]->GetXaxis()->SetLabelSize(0.20);
    ratio.hist[i][ch]->GetXaxis()->SetTitleSize(0.16);

    ratio.hist[i][ch]->SetMinimum(0.6);
    ratio.hist[i][ch]->SetMaximum(1.4);


    TGraphErrors *thegraphRatio = new TGraphErrors(ratio.hist[i][ch]);
    thegraphRatio->SetFillStyle(1001);
    thegraphRatio->SetFillColor(chatch);
    thegraphRatio->SetName("thegraphRatio");

    //thegraphRatio->Draw("e2SAME");
    for (int nbin=0; nbin<ratio.hist[i][ch]->GetNbinsX(); nbin++) ratio.hist[i][ch]->SetBinError(nbin+1,0.001); 


    leg[i]=addlegend(leg[i],"Data", Data.hist[i][ch]);
    //leg[i]=addlegend(leg[i],"VV", VV.hist[i][ch]);
    leg[i]=addlegend(leg[i],"W+Jets", NONWZ.hist[i][ch]);
    leg[i]=addlegend(leg[i],"Single t", SingleT.hist[i][ch]);
    leg[i]=addlegend(leg[i],"Z+Jets",  Z.hist[i][ch]);
    leg[i]=addlegend(leg[i],"t#bar{t}", TTbar.hist[i][ch]);
    leg[i]=addlegend(leg[i],"QCD(#mu)", QCD.hist[i][ch]);
    leg[i]->AddEntry("thegraph","Uncertainty","f");
    leg[i]->Draw("SAME");


    pad[i][1]->cd();

    ratio.hist[i][ch]->Draw();
    thegraphRatio->Draw("e2");
    ratio.hist[i][ch]->Draw("SAME");

    /***********************
           CMS Legend
     **********************/
    cPlots[i]->cd();
    pad[i][0]->cd();
    if(plots=="lepton") pad[i][0]->SetLogy();

    TString htitleCMSChannel;
    if(ch==0) htitleCMSChannel="#mu^{#pm}+jets channel";
    if(ch==1) htitleCMSChannel="e^{#pm}+jets channel";
    if(ch==2) htitleCMSChannel="l^{#pm}+jets channel";

    titlePr  = new TLatex(-20.,50.,"Preliminary");
    titlePr->SetNDC();
    titlePr->SetTextAlign(12);
    titlePr->SetX(0.25);
    titlePr->SetY(0.93);
    titlePr->SetTextColor(2);
    titlePr->SetTextFont(42);
    titlePr->SetTextSize(0.05);
    titlePr->SetTextSizePixels(24);
    titlePr->Draw("SAME");

    title  = new TLatex(-20.,50.,"CMS #sqrt{s} = 13TeV, L = 1000 pb^{-1}");
    //title  = new TLatex(-20.,50.,"CMS #sqrt{s} = 8TeV, L = 19.5 fb^{-1}");
    title->SetNDC();
    title->SetTextAlign(12);
    title->SetX(0.20);
    title->SetY(0.83);
    title->SetTextFont(42);
    title->SetTextSize(0.05);
    title->SetTextSizePixels(24);
    title->Draw("SAME");

    chtitle  = new TLatex(-20.,50.,htitleCMSChannel);
    chtitle->SetNDC();
    chtitle->SetTextAlign(12);
    chtitle->SetX(0.20);
    chtitle->SetY(0.75);
    chtitle->SetTextFont(42);
    chtitle->SetTextSize(0.05);
    chtitle->SetTextSizePixels(24);
    chtitle->Draw("SAME");

    TString dirfigname_pdf=dirnameIn + "figures/pdf/";
    // make a dir if it does not exist!!
    gSystem->mkdir(dirfigname_pdf,       kTRUE);

    TString dirfigname_png=dirnameIn + "figures/png/";
    // make a dir if it does not exist!!
    gSystem->mkdir(dirfigname_png,       kTRUE);

    char fig[200];
    sprintf(fig,"f_%s",TTbar.hist[i][ch]->GetName());

    if(ch==2){ // Change the name!
      for(int c=0;c<199;c++){
	if(fig[c]=='m' && fig[c+1]=='u'){ 
	  fig[c]='l';	
	  fig[c+1]='e';	
	  fig[c+2]='p';	
	  fig[c+3]='j';	
	  fig[c+4]='e';	
	  fig[c+5]='t';	
	}
      }
    }//if(channel=all)

    // PDF
    dirfigname_pdf=dirfigname_pdf + fig;
    cPlots[i]->SaveAs(dirfigname_pdf + ".pdf");

    // PNG
    dirfigname_png=dirfigname_png + fig;
    cPlots[i]->SaveAs(dirfigname_png + ".png");

  }//for(histograms)

}

/////////////////  Functions  //////////////////////

TLegend* createlegend(TLegend *leg, TString plots=""){

  float legx1=0.68; 
  float legy1=0.58;
  float legx2=0.90;
  float legy2=0.89;

  //if(plots=="dilepton") legx1=0.75;

  leg = new TLegend(legx1,legy1,legx2,legy2);
  //leg->SetTextAlign(32);
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  //leg->SetLineWidth(4);
  //leg->SetTextFont(62); 
  leg->SetTextFont(102); // Events in the leg!
  leg->SetTextSize(0.04);

  return leg;
}

TLegend* addlegend(TLegend *leg, TString name, TH1* histo){

  //Number of events
  /////////////////////////////////////////////
  char com[200]; 
  TString test;
  sprintf(com,"%i",(histo->Integral()));

  int cn=strlen(com);
  int hn=name.Sizeof();
  if (name.Contains("t#"))  hn=3;  
  if (name.Contains("#mu")) hn=7;
  cn=14-(hn+cn); // 14 comes from the max number of character.
  for(int i=0;i<cn;i++){
    name=name + " ";
  }
  
  name=name+com;
  ////////////////////////////////////////////
  if(name.Contains("Data")) leg->AddEntry(histo,name,"PL");
  else leg->AddEntry(histo,name,"F");

  return leg;

}

Histograms addhistogram(Histograms histo, Histograms histoIn){
  
  int histos = 8;

    for(int j=0; j<histos; j++){
      for(int k=0; k<3; k++){
	histo.hist[j][k]->Add(histo.hist[j][k],histoIn.hist[j][k]);	
      }
    }

  return histo;

}

Histograms addstack(Histograms stack, Histograms histoIn){

  int histos = 8;
    for(int j=0; j<histos; j++){
      for(int k=0; k<3; k++){
	stack.mc[j][k]->Add(histoIn.hist[j][k]);	
      }
    }
    
    return stack;
}

void setupdraw(Histograms h, int color) {

  int histos = 8;
  int bin[8];

  for(int j=0; j<histos; j++){
    bin[j]=1;//histos
  }
  
  for(int j=0; j<histos; j++){
    for(int k=0; k<3; k++){
      
      h.hist[j][k]->Rebin(bin[j]);
      
      TString htitle=h.hist[j][k]->GetTitle();

      if(k==0) htitle.Resize(htitle.Sizeof()-10);
      if(k==1) htitle.Resize(htitle.Sizeof()-8);
      if(k==2) htitle.Resize(htitle.Sizeof()-6);
      

      if(htitle.Contains("multiplicity") || htitle.Contains("PV")) h.hist[j][k]->GetXaxis()->SetTitle(htitle);
      else if(htitle.Contains("#Phi") || htitle.Contains("#phi")) h.hist[j][k]->GetXaxis()->SetTitle("|" + htitle + "| [rad]"); 
      else h.hist[j][k]->GetXaxis()->SetTitle(htitle + " [GeV]");

      h.hist[j][k]->SetLineColor(1);
      h.hist[j][k]->SetFillColor(color);
      h.hist[j][k]->SetFillStyle(1001);

      h.hist[j][k]->SetBinContent(h.hist[j][k]->GetNbinsX(),
				     (h.hist[j][k]->GetBinContent(h.hist[j][k]->GetNbinsX()+1)+h.hist[j][k]->GetBinContent(h.hist[j][k]->GetNbinsX())));
      h.hist[j][k]->SetBinContent(h.hist[j][k]->GetNbinsX()+1,0);

      if(j==6) h.hist[6][k]->GetXaxis()->SetTitle("b-jet multiplicity");

    }
  }
  
  
}

Histograms addhistograms(Histograms histoIn){
  
  int histos = 8;
  
  for(int j=0; j<histos; j++){
    histoIn.hist[j][2]=(TH1F*)histoIn.hist[j][0]->Clone();
    
    TString htitle=histoIn.hist[j][0]->GetTitle();
    htitle.Resize(htitle.Sizeof()-8);

    histoIn.hist[j][2]->SetTitle(htitle);
    histoIn.hist[j][2]->Add(histoIn.hist[j][2],histoIn.hist[j][1]);

  }

  return histoIn;
}

void scalehistograms(Histograms histoIn, int channel, float SF){
  
  int histos = 8;
  
    for(int j=0; j<histos; j++){
      histoIn.hist[j][channel]->Scale(SF);
    }
    
}


Histograms loadhistograms(TString plots,TString namefile){
  
  Histograms histofile;
  
  TFile *file=NULL;//new TFile(namefile);
  
  file = TFile::Open(namefile);
  cout << "loading " << plots << " " << namefile << endl; 
  
  TString channel[2];
  channel[0]="mujets";
  channel[1]="ejets";
  
  for(int ch=0;ch<2;ch++){

    histofile.hist[0][ch]=(TH1F*)file->Get("hPV_"+channel[ch]+"_"+plots);

    histofile.hist[1][ch]=(TH1F*)file->Get("hMET_"+channel[ch]+"_"+plots);
    histofile.hist[2][ch]=(TH1F*)file->Get("hCSV_Jet0_"+channel[ch]+"_"+plots);
    //histofile.hist[2][ch]=(TH1F*)file->Get("hMET_Phi_"+channel[ch]+"_"+plots);

    histofile.hist[3][ch]=(TH1F*)file->Get("hNJets_"+channel[ch]+"_"+plots);
    histofile.hist[4][ch]=(TH1F*)file->Get("hNBtagJets_"+channel[ch]+"_"+plots);    

    histofile.hist[5][ch]=(TH1F*)file->Get("hLepPt_"+channel[ch]+"_"+plots);
    histofile.hist[6][ch]=(TH1F*)file->Get("hLepPhi_"+channel[ch]+"_"+plots);
    histofile.hist[7][ch]=(TH1F*)file->Get("hLepEta_"+channel[ch]+"_"+plots);

  }

  cout << "All the histograms have been loading successfully!!!"  << endl; 

  return histofile;

  file.Close();

}

