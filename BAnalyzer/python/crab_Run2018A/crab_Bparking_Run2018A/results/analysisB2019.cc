//Program to extract data (variables and vector format) of the root file and create histograms. This root file is made of others root files, therefore, has multiple entries which we access with a diferent method than a normal root file.

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <stdio.h>      /* printf, fgets */
#include <stdlib.h>     /* atof */
#include <math.h>       /* sin */
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include <vector>
#ifndef ROOT_TLatex
#define ROOT_TLatex
#ifndef ROOTiosfwd
#include "Riosfwd.h"
#endif
#ifndef ROOT_TText
#include "TText.h"
#endif
#ifndef ROOT_TAttLine
#include "TAttLine.h"
#endif
//#include "RooCBShape.h"  //Crystal Ball
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
//#include "RooRealVar.h"
//#include "RooDataSet.h"
//#include "RooGaussian.h"
#include "TCanvas.h"
//#include "RooPlot.h"
#include "TAxis.h"
//using namespace RooFit ;


//gSystem->AddIncludePath("/home/raphael/Desktop/analysis_2018");

using namespace std;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int analysisB2019()
{
	//call a file for a histogram style (Optional)
	//gROOT->LoadMacro("styleTDR.C"); 
	//setTDRStyle();

	//counters for file f1 (pythia)
	int count_Total_Events_pythia = 0;
	
	//Lorentz Vector
	TLorentzVector mu_1;
	TLorentzVector mu_2;	
	
	//Variaveis	and Vectors
	//int Total_Events = 0;

	std::vector<double>* D0mass = 0.;
	std::vector<double>* Dsmass = 0.;
	std::vector<double>* D0_VtxProb = 0.;
	std::vector<double>* D0_VtxPosx = 0.;
	std::vector<double>* D0_VtxPosy = 0.;
	std::vector<double>* D0_VtxPosz = 0.;
	std::vector<double>* D0_Vtxerrx = 0.;
	std::vector<double>* D0_Vtxerry = 0.;
	std::vector<double>* D0_Vtxerrz = 0.;
	std::vector<double>* D0eta = 0.;
	std::vector<double>* D0phi = 0.;
	std::vector<double>* Dseta = 0.;
	std::vector<double>* Dsphi = 0.;
	std::vector<double>* TrkKpt = 0.;
	std::vector<double>* Trkpipt = 0.;
	std::vector<double>* TrkSpt = 0.;
	std::vector<double>* D0pt = 0.;
	std::vector<double>* Dspt = 0.;
	std::vector<double>* DSDeltaR = 0.;
	std::vector<double>* TrkKnhits = 0.;
	std::vector<double>* Trkpinhits = 0.;
	std::vector<double>* TrkSnhits = 0.;
	std::vector<double>* TrkKchi2 = 0.;
	std::vector<double>* Trkpichi2 = 0.;
	std::vector<double>* TrkSchi2 = 0.;
	std::vector<double>* TrkKdxy = 0.;
	std::vector<double>* Trkpidxy = 0.;
	std::vector<double>* TrkSdxy = 0.;
	std::vector<double>* TrkKdz = 0.;
	std::vector<double>* Trkpidz = 0.;
	std::vector<double>* TrkSdz = 0.;
	std::vector<double>* TrkKeta = 0.;
	std::vector<double>* Trkpieta = 0.;
	std::vector<double>* TrkSeta = 0.;
	std::vector<double>* TrkKphi = 0.;
	std::vector<double>* Trkpiphi = 0.;
	std::vector<double>* TrkSphi = 0.;
	std::vector<double>* TrkScharge = 0.;
	std::vector<double>* D0fromDSsXY_vec = 0.;
	std::vector<double>* D0mass = 0.;
	std::vector<double>* D0mass = 0.;
	std::vector<double>* D0mass = 0.;


	//double M = 0.;
	//double Pt = 0.;
	//double Eta = 0.;
	//double Rapidity = 0.;
	
	//*****Creating Histgrams******************************************************	

	//Histogramas cinematics quantities of the D0
	TH1F *D0mass_Histo = new TH1F("D0mass_Histo","D0mass_Histo",100,1.5,2.3);
	D0mass_Histo->SetTitle("Mass distribuition of the D0 ; Mass [GeV] ; Events ");
	D0mass_Histo->SetName("D0mass_Histo");

	TH1F *D0eta_Histo = new TH1F("D0eta_Histo","D0eta_Histo",100,-3,3);
	D0eta_Histo->SetTitle("Pseudo-rapidity distribuition of the D0 ; #eta ; Events ");
	D0eta_Histo->SetName("D0eta_Histo");

	TH1F *D0phi_Histo = new TH1F("D0phi_Histo","D0phi_Histo",100,-5,5);
	D0phi_Histo->SetTitle("#Phi distribuition of the D0 ; #Phi [º] ; Events ");
	D0phi_Histo->SetName("D0phi_Histo");

	//Histogramas cinematics quantities of the D*
	TH1F *Dsmass_Histo = new TH1F("Dsmass_Histo","Dsmass_Histo",100,1.6,2.4);
	Dsmass_Histo->SetTitle("Mass distribuition of the D* ; Mass [GeV] ; Events ");
	Dsmass_Histo->SetName("Dsmass_Histo");

	TH1F *Dseta_Histo = new TH1F("Dseta_Histo","Dseta_Histo",100,-3,3);
	Dseta_Histo->SetTitle("Pseudo-rapidity distribuition of the D* ; #eta ; Events ");
	Dseta_Histo->SetName("Dseta_Histo");

	TH1F *Dsphi_Histo = new TH1F("Dsphi_Histo","Dsphi_Histo",100,-5,5);
	Dsphi_Histo->SetTitle("#Phi distribuition of the D* ; #Phi [º] ; Events ");
	Dsphi_Histo->SetName("Dsphi_Histo");

	//Histogramas cinematics quantities of the SlowPion
	TH1F *TrkSpt_Histo = new TH1F("TrkSpt_Histo","TrkSpt_Histo",100,0,2.5);
	TrkSpt_Histo->SetTitle("pt distribuition of the SlowPion ; Mass [GeV] ; Events ");
	TrkSpt_Histo->SetName("TrkSpt_Histo");

	TH1F *TrkSnhits_Histo = new TH1F("TrkSnhits_Histo","TrkSnhits_Histo",100,0,45);
	TrkSnhits_Histo->SetTitle("nhits distribuition of the SlowPion ; Number of hits ; Events ");
	TrkSnhits_Histo->SetName("TrkSnhits_Histo");

	TH1F *TrkSdxy_Histo = new TH1F("TrkSdxy_Histo","TrkSdxy_Histo",100,-2,2);
	TrkSdxy_Histo->SetTitle("dxy distribuition of the SlowPion ; dx [cm] ; Events ");
	TrkSdxy_Histo->SetName("TrkSdxy_Histo");

	TH1F *TrkSdz_Histo = new TH1F("TrkSdz_Histo","TrkSdz_Histo",100,-4,4);
	TrkSdz_Histo->SetTitle("dz distribuition of the SlowPion ; dz [cm] ; Events ");
	TrkSdz_Histo->SetName("TrkSdz_Histo");

	TH1F *TrkSeta_Histo = new TH1F("TrkSeta_Histo","TrkSeta_Histo",100,-2.5,2.5);
	TrkSeta_Histo->SetTitle("eta distribuition of the SlowPion ; #eta ; Events ");
	TrkSeta_Histo->SetName("TrkSeta_Histo");

	TH1F *TrkSphi_Histo = new TH1F("TrkSphi_Histo","TrkSphi_Histo",100,-4,4);
	TrkSphi_Histo->SetTitle("Phi distribuition of the SlowPion ; #Phi [º] ; Events ");
	TrkSphi_Histo->SetName("TrkSphi_Histo");

	TH1F *TrkSchi2_Histo = new TH1F("TrkSchi2_Histo","TrkSchi2_Histo",10,0,3);
	TrkSchi2_Histo->SetTitle("Chi2 distribuition of the SlowPion ; #Chi^{2} ; Events ");
	TrkSchi2_Histo->SetName("TrkSchi2_Histo");

	//Histogramas cinematics quantities of the Pion
	TH1F *Trkpipt_Histo = new TH1F("Trkpipt_Histo","Trkpipt_Histo",100,0,10);
	Trkpipt_Histo->SetTitle("pt distribuition of the Pion ; Mass [GeV] ; Events ");
	Trkpipt_Histo->SetName("Trkpipt_Histo");

	TH1F *Trkpinhits_Histo = new TH1F("Trkpinhits_Histo","Trkpinhits_Histo",100,0,40);
	Trkpinhits_Histo->SetTitle("nhits distribuition of the Pion ; Number of hits ; Events ");
	Trkpinhits_Histo->SetName("Trkpinhits_Histo");

	TH1F *Trkpidxy_Histo = new TH1F("Trkpidxy_Histo","Trkpidxy_Histo",100,-0.8,0.8);
	Trkpidxy_Histo->SetTitle("dxy distribuition of the Pion ; dx [cm] ; Events ");
	Trkpidxy_Histo->SetName("Trkpidxy_Histo");

	TH1F *Trkpidz_Histo = new TH1F("Trkpidz_Histo","Trkpidz_Histo",100,-1.5,1.5);
	Trkpidz_Histo->SetTitle("dz distribuition of the Pion ; dz [cm] ; Events ");
	Trkpidz_Histo->SetName("Trkpidz_Histo");

	TH1F *Trkpieta_Histo = new TH1F("Trkpieta_Histo","Trkpieta_Histo",100,-2.5,2.5);
	Trkpieta_Histo->SetTitle("eta distribuition of the Pion ; #eta ; Events ");
	Trkpieta_Histo->SetName("Trkpieta_Histo");

	TH1F *Trkpiphi_Histo = new TH1F("Trkpiphi_Histo","Trkpiphi_Histo",100,-4,4);
	Trkpiphi_Histo->SetTitle("phi distribuition of the SlowPion ; #Phi [º] ; Events ");
	Trkpiphi_Histo->SetName("Trkpiphi_Histo");

	TH1F *Trkpichi2_Histo = new TH1F("Trkpichi2_Histo","Trkpichi2_Histo",10,0,3);
	Trkpichi2_Histo->SetTitle("Chi2 distribuition of the Pion ; #Chi^{2} ; Events ");
	Trkpichi2_Histo->SetName("Trkpichi2_Histo");

	//Histogramas cinematics quantities of the Kaon
	TH1F *TrkKpt_Histo = new TH1F("TrkKpt_Histo","TrkKpt_Histo",100,0,10);
	TrkKpt_Histo->SetTitle("pt distribuition of the Kaon ; Mass [GeV] ; Events ");
	TrkKpt_Histo->SetName("TrkKpt_Histo");

	TH1F *TrkKnhits_Histo = new TH1F("TrkKnhits_Histo","TrkKnhits_Histo",100,0,40);
	TrkKnhits_Histo->SetTitle("nhits distribuition of the Kaon ; Number of hits ; Events ");
	TrkKnhits_Histo->SetName("TrkKnhits_Histo");

	TH1F *TrkKdxy_Histo = new TH1F("TrkKdxy_Histo","TrkKdxy_Histo",100,-0.8,0.8);
	TrkKdxy_Histo->SetTitle("dxy distribuition of the Kaon ; dx [cm] ; Events ");
	TrkKdxy_Histo->SetName("TrkKdxy_Histo");

	TH1F *TrkKdz_Histo = new TH1F("TrkKdz_Histo","TrkKdz_Histo",100,-1.5,1.5);
	TrkKdz_Histo->SetTitle("dz distribuition of the Kaon ; dz [cm] ; Events ");
	TrkKdz_Histo->SetName("TrkKdz_Histo");

	TH1F *TrkKeta_Histo = new TH1F("TrkKeta_Histo","TrkKeta_Histo",100,-2.5,2.5);
	TrkKeta_Histo->SetTitle("eta distribuition of the Kaon ; #eta ; Events ");
	TrkKeta_Histo->SetName("TrkKeta_Histo");

	TH1F *TrkKphi_Histo = new TH1F("TrkKphi_Histo","TrkKphi_Histo",100,-4,4);
	TrkKphi_Histo->SetTitle("phi distribuition of the Kaon ; #Phi [º] ; Events ");
	TrkKphi_Histo->SetName("TrkKphi_Histo");

	TH1F *TrkKchi2_Histo = new TH1F("TrkKchi2_Histo","TrkKchi2_Histo",10,0,3);
	TrkKchi2_Histo->SetTitle("Chi2 distribuition of the Kaon ; #Chi^{2} ; Events ");
	TrkKchi2_Histo->SetName("TrkKchi2_Histo");

	//End Histograms
	//-----------------------------------------------------------------------------

	//-------Reading the root file and the tree-------------------------------------	
	TFile *f1 = new TFile("D0DstarData_test.root");
	TTree *t1 = (TTree*)f1->Get("analysis/data");

	//---------------------------------------------------------------------------------
	// addressing the memory to vector and variables for file

	//For Variables
	//TBranch *b_Total_Events = t1->GetBranch("Total_Events");
	//b_Total_Events->SetAddress(&Total_Events);
	

	//For Vectors	
	TBranch *b_D0mass = t1->GetBranch("D0mass");
	b_D0mass->SetAddress(&D0mass);
	TBranch *b_D0eta = t1->GetBranch("D0eta");
	b_D0eta->SetAddress(&D0eta);
	TBranch *b_D0phi = t1->GetBranch("D0phi");
	b_D0phi->SetAddress(&D0phi);

	TBranch *b_Dsmass = t1->GetBranch("Dsmass");
	b_Dsmass->SetAddress(&Dsmass);
	TBranch *b_Dseta = t1->GetBranch("Dseta");
	b_Dseta->SetAddress(&Dseta);
	TBranch *b_Dsphi = t1->GetBranch("Dsphi");
	b_Dsphi->SetAddress(&Dsphi);

	TBranch *b_TrkSpt = t1->GetBranch("TrkSpt");
	b_TrkSpt->SetAddress(&TrkSpt);
	TBranch *b_TrkSnhits = t1->GetBranch("TrkSnhits");
	b_TrkSnhits->SetAddress(&TrkSnhits);
	TBranch *b_TrkSdxy = t1->GetBranch("TrkSdxy");
	b_TrkSdxy->SetAddress(&TrkSdxy);
	TBranch *b_TrkSdz = t1->GetBranch("TrkSdz");
	b_TrkSdz->SetAddress(&TrkSdz);
	TBranch *b_TrkSeta = t1->GetBranch("TrkSeta");
	b_TrkSeta->SetAddress(&TrkSeta);
	TBranch *b_TrkSphi = t1->GetBranch("TrkSphi");
	b_TrkSphi->SetAddress(&TrkSphi);
	TBranch *b_TrkSchi2 = t1->GetBranch("TrkSchi2");
	b_TrkSchi2->SetAddress(&TrkSchi2);

	TBranch *b_Trkpipt = t1->GetBranch("Trkpipt");
	b_Trkpipt->SetAddress(&Trkpipt);
	TBranch *b_Trkpinhits = t1->GetBranch("Trkpinhits");
	b_Trkpinhits->SetAddress(&Trkpinhits);
	TBranch *b_Trkpidxy = t1->GetBranch("Trkpidxy");
	b_Trkpidxy->SetAddress(&Trkpidxy);
	TBranch *b_Trkpidz = t1->GetBranch("Trkpidz");
	b_Trkpidz->SetAddress(&Trkpidz);
	TBranch *b_Trkpieta = t1->GetBranch("Trkpieta");
	b_Trkpieta->SetAddress(&Trkpieta);
	TBranch *b_Trkpiphi = t1->GetBranch("Trkpiphi");
	b_Trkpiphi->SetAddress(&Trkpiphi);
	TBranch *b_Trkpichi2 = t1->GetBranch("Trkpichi2");
	b_Trkpichi2->SetAddress(&Trkpichi2);

	TBranch *b_TrkKpt = t1->GetBranch("TrkKpt");
	b_TrkKpt->SetAddress(&TrkKpt);
	TBranch *b_TrkKnhits = t1->GetBranch("TrkKnhits");
	b_TrkKnhits->SetAddress(&TrkKnhits);
	TBranch *b_TrkKdxy = t1->GetBranch("TrkKdxy");
	b_TrkKdxy->SetAddress(&TrkKdxy);
	TBranch *b_TrkKdz = t1->GetBranch("TrkKdz");
	b_TrkKdz->SetAddress(&TrkKdz);
	TBranch *b_TrkKeta = t1->GetBranch("TrkKeta");
	b_TrkKeta->SetAddress(&TrkKeta);
	TBranch *b_TrkKphi = t1->GetBranch("TrkKphi");
	b_TrkKphi->SetAddress(&TrkKphi);
	TBranch *b_TrkKchi2 = t1->GetBranch("TrkKchi2");
	b_TrkKchi2->SetAddress(&TrkKchi2);


	


	//**********************************************************		
	//Reading Number of tree entries for file f1
	Long64_t nentries = t1->GetEntries();
	cout<< "Number of tree entries: "<< nentries <<std::endl;

	Long64_t GetEntriesFast = t1->GetEntriesFast();
	//cout<< "GetEntriesFast: "<< GetEntriesFast <<std::endl;

	Long64_t nbytes = 0, nb = 0, i=0;
	for (Long64_t jentry=0; jentry < nentries; jentry++) //loop tree entries for file f1
	{
      	Long64_t ientry = t1->LoadTree(jentry);
      	//std::cout << "nentries " << nentries << " ientry " << ientry << " jentry " << jentry <<std::endl;
   
      	if (ientry < 0) break;
		//b_Total_Events->GetEntry(ientry);
		
		//counter_Muon_pythia += Muons;
		//cout << "Total_Events: "<< Total_Events << endl;

       	b_D0mass->GetEntry(ientry);
		b_Dsmass->GetEntry(ientry);
		b_D0eta->GetEntry(ientry);
		b_Dseta->GetEntry(ientry);
		b_D0phi->GetEntry(ientry);
		b_Dsphi->GetEntry(ientry);

		b_TrkSpt->GetEntry(ientry);
		b_TrkSnhits->GetEntry(ientry);
		b_TrkSdxy->GetEntry(ientry);
		b_TrkSdz->GetEntry(ientry);
		b_TrkSeta->GetEntry(ientry);
		b_TrkSphi->GetEntry(ientry);
		b_TrkSchi2->GetEntry(ientry);

		b_Trkpipt->GetEntry(ientry);
		b_Trkpinhits->GetEntry(ientry);
		b_Trkpidxy->GetEntry(ientry);
		b_Trkpidz->GetEntry(ientry);
		b_Trkpieta->GetEntry(ientry);
		b_Trkpiphi->GetEntry(ientry);
		b_Trkpichi2->GetEntry(ientry);

		b_TrkKpt->GetEntry(ientry);
		b_TrkKnhits->GetEntry(ientry);
		b_TrkKdxy->GetEntry(ientry);
		b_TrkKdz->GetEntry(ientry);
		b_TrkKeta->GetEntry(ientry);
		b_TrkKphi->GetEntry(ientry);
		b_TrkKchi2->GetEntry(ientry);
	

		//For D0
		for(Long64_t i=0; i < D0mass->size(); i++)
		{  
			D0mass_Histo->Fill(D0mass->at(i));
  		}

		for(Long64_t i=0; i < D0eta->size(); i++)
		{  
			D0eta_Histo->Fill(D0eta->at(i));
  		}

		for(Long64_t i=0; i < D0phi->size(); i++)
		{  
			D0phi_Histo->Fill(D0phi->at(i));
  		}

		//For D*
		for(Long64_t i=0; i < Dsmass->size(); i++)
		{  
			Dsmass_Histo->Fill(Dsmass->at(i));
  		}

		for(Long64_t i=0; i < Dseta->size(); i++)
		{  
			Dseta_Histo->Fill(Dseta->at(i));
  		}

		for(Long64_t i=0; i < Dsphi->size(); i++)
		{  
			Dsphi_Histo->Fill(Dsphi->at(i));
  		}

		//For SlowPion
		for(Long64_t i=0; i < TrkSpt->size(); i++)
		{  
			TrkSpt_Histo->Fill(TrkSpt->at(i));
  		}

		for(Long64_t i=0; i < TrkSnhits->size(); i++)
		{  
			TrkSnhits_Histo->Fill(TrkSnhits->at(i));
  		}

		for(Long64_t i=0; i < TrkSdxy->size(); i++)
		{  
			TrkSdxy_Histo->Fill(TrkSdxy->at(i));
			//std::cout << "TrkSdxy->at(i): " << TrkSdxy->at(i) << endl;
  		}

		for(Long64_t i=0; i < TrkSdz->size(); i++)
		{  
			TrkSdz_Histo->Fill(TrkSdz->at(i));
  		}

		for(Long64_t i=0; i < TrkSeta->size(); i++)
		{  
			TrkSeta_Histo->Fill(TrkSeta->at(i));
  		}

		for(Long64_t i=0; i < TrkSphi->size(); i++)
		{  
			TrkSphi_Histo->Fill(Dsphi->at(i));
  		}
	
		for(Long64_t i=0; i < TrkSchi2->size(); i++)
		{  
			TrkSchi2_Histo->Fill(TrkSchi2->at(i));
  		}

		//For Pion
		for(Long64_t i=0; i < Trkpipt->size(); i++)
		{  
			Trkpipt_Histo->Fill(Trkpipt->at(i));
  		}

		for(Long64_t i=0; i < Trkpinhits->size(); i++)
		{  
			Trkpinhits_Histo->Fill(Trkpinhits->at(i));
  		}

		for(Long64_t i=0; i < Trkpidxy->size(); i++)
		{  
			Trkpidxy_Histo->Fill(Trkpidxy->at(i));
  		}

		for(Long64_t i=0; i < Trkpidz->size(); i++)
		{  
			Trkpidz_Histo->Fill(Trkpidz->at(i));
  		}

		for(Long64_t i=0; i < Trkpieta->size(); i++)
		{  
			Trkpieta_Histo->Fill(Trkpieta->at(i));
  		}

		for(Long64_t i=0; i < Trkpiphi->size(); i++)
		{  
			Trkpiphi_Histo->Fill(Trkpiphi->at(i));
  		}

		for(Long64_t i=0; i < Trkpichi2->size(); i++)
		{  
			Trkpichi2_Histo->Fill(Trkpichi2->at(i));
  		}

		//For Kon
		for(Long64_t i=0; i < TrkKpt->size(); i++)
		{  
			TrkKpt_Histo->Fill(TrkKpt->at(i));
  		}

		for(Long64_t i=0; i < TrkKnhits->size(); i++)
		{  
			TrkKnhits_Histo->Fill(TrkKnhits->at(i));
  		}

		for(Long64_t i=0; i < TrkKdxy->size(); i++)
		{  
			TrkKdxy_Histo->Fill(TrkKdxy->at(i));
  		}

		for(Long64_t i=0; i < TrkKdz->size(); i++)
		{  
			TrkKdz_Histo->Fill(TrkKdz->at(i));
  		}

		for(Long64_t i=0; i < TrkKeta->size(); i++)
		{  
			TrkKeta_Histo->Fill(TrkKeta->at(i));
  		}

		for(Long64_t i=0; i < TrkKphi->size(); i++)
		{  
			TrkKphi_Histo->Fill(TrkKphi->at(i));
  		}

		for(Long64_t i=0; i < TrkKchi2->size(); i++)
		{  
			TrkKchi2_Histo->Fill(TrkKchi2->at(i));
  		}

			
	}//End loop tree entries for file f1


	//TFile f_analysis("Jpsi.root","recreate"); //Creates root file
	//TTree t_analysis("t_analysis","analise_Tree"); //Creates a Tree
	//t_analysis.Branch("vectorInvariantMassJpsi",&vectorInvariantMassJpsi); //Creates a branch
	//t_analysis.Fill();
	//h_DimuonsOppositeChargeEtaJpsi_M2_Double.Write(); //Write() if a file is open, this function writes a root objectics on it.
	//t_analysis.Write();  //Write in the root file

	//cout << "        " << endl;	
	//cout << "======================================================== " << endl;
	//cout << "count_Total_Events_pythia: "<< count_Total_Events_pythia << endl;	
	//cout << "======================================================== " << endl;
	//cout << "        " << endl;

    //=========================================================================	
	//Creating Canvas
	TCanvas* c1 = new TCanvas("c1","Canvas 1 - behavior of the D",1200,600);
	c1->Divide(2);
	c1->cd(1);
	D0mass_Histo->SetLineColor(kBlue);
	D0mass_Histo->SetMarkerStyle(7);
	D0mass_Histo->SetMarkerStyle(21);
	//D0mass_Histo->SetStats(0);
	D0mass_Histo->SetMarkerColor(kBlue);
	D0mass_Histo->SetFillColor(kBlue);
	D0mass_Histo->Draw("HIST");
	//D0mass_Histo->Draw("e1pSAME");
	TLegend* leg_D0mass = new TLegend(0.82,0.5,0.95,0.65);
   	leg_D0mass->SetFillColor(kWhite);
	leg_D0mass->SetFillStyle(1001);
	leg_D0mass->SetBorderSize(0);
	//leg_D0mass->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_D0mass->AddEntry(D0mass_Histo,"D0mass","L");
	leg_D0mass->Draw();
	//-------------------------------------------------------------------------
	c1->cd(2);
	Dsmass_Histo->SetLineColor(kBlue);
	Dsmass_Histo->SetMarkerStyle(7);
	Dsmass_Histo->SetMarkerStyle(21);
	//Dsmass_Histo->SetStats(0);
	Dsmass_Histo->SetMarkerColor(kBlue);
	Dsmass_Histo->SetFillColor(kBlue);
	Dsmass_Histo->Draw("HIST");
	//Dsmass_Histo->Draw("e1pSAME");
	TLegend* leg_Dsmass = new TLegend(0.82,0.5,0.95,0.65);
   	leg_Dsmass->SetFillColor(kWhite);
	leg_Dsmass->SetFillStyle(1001);
	leg_Dsmass->SetBorderSize(0);
	//leg_Dsmass->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_Dsmass->AddEntry(Dsmass_Histo,"Dsmass","L");
	leg_Dsmass->Draw();

	c1->SaveAs("/eos/user/r/ragomesd/analysisB2019/Canvas1.png");
	//=========================================================================	
	//Creating Canvas
	TCanvas* c2 = new TCanvas("c2","Canvas 2 - behavior of the D",1200,600);
	c2->Divide(2,2);
	c2->cd(1);
	D0eta_Histo->SetLineColor(kBlue);
	D0eta_Histo->SetMarkerStyle(7);
	D0eta_Histo->SetMarkerStyle(21);
	//D0eta_Histo->SetStats(0);
	D0eta_Histo->SetMarkerColor(kBlue);
	D0eta_Histo->SetFillColor(kBlue);
	D0eta_Histo->Draw("HIST");
	//D0eta_Histo->Draw("e1pSAME");
	TLegend* leg_D0eta = new TLegend(0.82,0.5,0.95,0.65);
   	leg_D0eta->SetFillColor(kWhite);
	leg_D0eta->SetFillStyle(1001);
	leg_D0eta->SetBorderSize(0);
	//leg_D0eta->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_D0eta->AddEntry(D0eta_Histo,"D0eta","L");
	leg_D0eta->Draw();
	//-------------------------------------------------------------------------
	c2->cd(2);
	Dseta_Histo->SetLineColor(kBlue);
	Dseta_Histo->SetMarkerStyle(7);
	Dseta_Histo->SetMarkerStyle(21);
	//Dseta_Histo->SetStats(0);
	Dseta_Histo->SetMarkerColor(kBlue);
	Dseta_Histo->SetFillColor(kBlue);
	Dseta_Histo->Draw("HIST");
	//Dsmass_Histo->Draw("e1pSAME");	
	TLegend* leg_Dseta = new TLegend(0.82,0.5,0.95,0.65);
   	leg_Dseta->SetFillColor(kWhite);
	leg_Dseta->SetFillStyle(1001);
	leg_Dseta->SetBorderSize(0);
	//leg_Dsmass->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_Dseta->AddEntry(Dseta_Histo,"Dseta","L");
	leg_Dseta->Draw();
	//-------------------------------------------------------------------------
	c2->cd(3);
	D0phi_Histo->SetLineColor(kBlue);
	D0phi_Histo->SetMarkerStyle(7);
	D0phi_Histo->SetMarkerStyle(21);
	//D0phi_Histo->SetStats(0);
	D0phi_Histo->SetMarkerColor(kBlue);
	D0phi_Histo->SetFillColor(kBlue);
	D0phi_Histo->Draw("HIST");
	//D0eta_Histo->Draw("e1pSAME");
	TLegend* leg_D0phi = new TLegend(0.82,0.5,0.95,0.65);
   	leg_D0phi->SetFillColor(kWhite);
	leg_D0phi->SetFillStyle(1001);
	leg_D0phi->SetBorderSize(0);
	//leg_D0phi->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_D0phi->AddEntry(D0phi_Histo,"D0phi","L");
	leg_D0phi->Draw();
	//-------------------------------------------------------------------------
	c2->cd(4);	
	Dsphi_Histo->SetLineColor(kBlue);
	Dsphi_Histo->SetMarkerStyle(7);
	Dsphi_Histo->SetMarkerStyle(21);
	//Dseta_Histo->SetStats(0);
	Dsphi_Histo->SetMarkerColor(kBlue);
	Dsphi_Histo->SetFillColor(kBlue);
	Dsphi_Histo->Draw("HIST");
	//Dsphi_Histo->Draw("e1pSAME");	
	TLegend* leg_Dsphi = new TLegend(0.82,0.5,0.95,0.65);
   	leg_Dsphi->SetFillColor(kWhite);
	leg_Dsphi->SetFillStyle(1001);
	leg_Dsphi->SetBorderSize(0);
	//leg_Dsmass->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_Dsphi->AddEntry(Dsphi_Histo,"Dsphi","L");
	leg_Dsphi->Draw();
	c2->SaveAs("/eos/user/r/ragomesd/analysisB2019/Canvas2.png");
	//=========================================================================	
	//Creating Canvas
	TCanvas* c3 = new TCanvas("c3","Canvas 3 - behavior of the SlowPion",1200,600);
	c3->Divide(2,2);
	c3->cd(1);
	TrkSpt_Histo->SetLineColor(kBlue);
	TrkSpt_Histo->SetMarkerStyle(7);
	TrkSpt_Histo->SetMarkerStyle(21);
	//TrkSpt_Histo->SetStats(0);
	TrkSpt_Histo->SetMarkerColor(kBlue);
	TrkSpt_Histo->SetFillColor(kBlue);
	TrkSpt_Histo->Draw("HIST");
	//TrkSpt_Histo->Draw("e1pSAME");
	TLegend* leg_TrkSpt = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkSpt->SetFillColor(kWhite);
	leg_TrkSpt->SetFillStyle(1001);
	leg_TrkSpt->SetBorderSize(0);
	//leg_TrkSpt->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkSpt->AddEntry(TrkSpt_Histo,"TrkSpt","L");
	leg_TrkSpt->Draw();
	//-------------------------------------------------------------------------
	c3->cd(2);
	TrkSnhits_Histo->SetLineColor(kBlue);
	TrkSnhits_Histo->SetMarkerStyle(7);
	TrkSnhits_Histo->SetMarkerStyle(21);
	//TrkSnhits_Histo->SetStats(0);
	TrkSnhits_Histo->SetMarkerColor(kBlue);
	TrkSnhits_Histo->SetFillColor(kBlue);
	TrkSnhits_Histo->Draw("HIST");
	//Dsphi_Histo->Draw("e1pSAME");
	TLegend* leg_TrkSnhits = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkSnhits->SetFillColor(kWhite);
	leg_TrkSnhits->SetFillStyle(1001);
	leg_TrkSnhits->SetBorderSize(0);
	//leg_TrkSnhits->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkSnhits->AddEntry(TrkSnhits_Histo,"TrkSnhits","L");
	leg_TrkSnhits->Draw();
	//-------------------------------------------------------------------------
	c3->cd(3);
	TrkSdxy_Histo->SetLineColor(kBlue);
	TrkSdxy_Histo->SetMarkerStyle(7);
	TrkSdxy_Histo->SetMarkerStyle(21);
	//TrkSnhits_Histo->SetStats(0);
	TrkSdxy_Histo->SetMarkerColor(kBlue);
	TrkSdxy_Histo->SetFillColor(kBlue);
	TrkSdxy_Histo->Draw("HIST");
	//TrkSdxy_Histo->Draw("e1pSAME");
	TLegend* leg_TrkSdxy = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkSdxy->SetFillColor(kWhite);
	leg_TrkSdxy->SetFillStyle(1001);
	leg_TrkSdxy->SetBorderSize(0);
	//leg_TrkSdxy->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkSdxy->AddEntry(TrkSdxy_Histo,"TrkSdxy","L");
	leg_TrkSdxy->Draw();
	//-------------------------------------------------------------------------
	c3->cd(4);
	TrkSdz_Histo->SetLineColor(kBlue);
	TrkSdz_Histo->SetMarkerStyle(7);
	TrkSdz_Histo->SetMarkerStyle(21);
	//TrkSdz_Histo->SetStats(0);
	TrkSdz_Histo->SetMarkerColor(kBlue);
	TrkSdz_Histo->SetFillColor(kBlue);
	TrkSdz_Histo->Draw("HIST");
	//TrkSdz_Histo->Draw("e1pSAME");
	TLegend* leg_TrkSdz = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkSdz->SetFillColor(kWhite);
	leg_TrkSdz->SetFillStyle(1001);
	leg_TrkSdz->SetBorderSize(0);
	//leg_TrkSdz->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkSdz->AddEntry(TrkSdz_Histo,"TrkSdz","L");
	leg_TrkSdz->Draw();
	c3->SaveAs("/eos/user/r/ragomesd/analysisB2019/Canvas3.png");
	//=========================================================================	
	//Creating Canvas
	TCanvas* c4 = new TCanvas("c4","Canvas 4 - behavior of the SlowPion",1200,600);
	c4->Divide(2,2);
	c4->cd(1);
	TrkSeta_Histo->SetLineColor(kBlue);
	TrkSeta_Histo->SetMarkerStyle(7);
	TrkSeta_Histo->SetMarkerStyle(21);
	//TrkSeta_Histo->SetStats(0);
	TrkSeta_Histo->SetMarkerColor(kBlue);
	TrkSeta_Histo->SetFillColor(kBlue);
	TrkSeta_Histo->Draw("HIST");
	//TrkSeta_Histo->Draw("e1pSAME");
	TLegend* leg_TrkSeta = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkSeta->SetFillColor(kWhite);
	leg_TrkSeta->SetFillStyle(1001);
	leg_TrkSeta->SetBorderSize(0);
	//leg_TrkSeta->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkSeta->AddEntry(TrkSeta_Histo,"TrkSeta","L");
	leg_TrkSeta->Draw();
	//-------------------------------------------------------------------------
	c4->cd(2);
	TrkSphi_Histo->SetLineColor(kBlue);
	TrkSphi_Histo->SetMarkerStyle(7);
	TrkSphi_Histo->SetMarkerStyle(21);
	//TrkSphi_Histo->SetStats(0);
	TrkSphi_Histo->SetMarkerColor(kBlue);
	TrkSphi_Histo->SetFillColor(kBlue);
	TrkSphi_Histo->Draw("HIST");
	//TrkSphi_Histo->Draw("e1pSAME");
	TLegend* leg_TrkSphi = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkSphi->SetFillColor(kWhite);
	leg_TrkSphi->SetFillStyle(1001);
	leg_TrkSphi->SetBorderSize(0);
	//leg_TrkSphi->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkSphi->AddEntry(TrkSphi_Histo,"TrkSphi","L");
	leg_TrkSphi->Draw();
	//-------------------------------------------------------------------------
	c4->cd(3);
	TrkSchi2_Histo->SetLineColor(kBlue);
	TrkSchi2_Histo->SetMarkerStyle(7);
	TrkSchi2_Histo->SetMarkerStyle(21);
	//TrkSchi2_Histo->SetStats(0);
	TrkSchi2_Histo->SetMarkerColor(kBlue);
	TrkSchi2_Histo->SetFillColor(kBlue);
	TrkSchi2_Histo->Draw("HIST");
	//TrkSchi2_Histo->Draw("e1pSAME");
	TLegend* leg_TrkSchi2 = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkSchi2->SetFillColor(kWhite);
	leg_TrkSchi2->SetFillStyle(1001);
	leg_TrkSchi2->SetBorderSize(0);
	//leg_TrkSchi2->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkSchi2->AddEntry(TrkSchi2_Histo,"TrkSchi2","L");
	leg_TrkSchi2->Draw();
	c4->SaveAs("/eos/user/r/ragomesd/analysisB2019/Canvas4.png");
	//=========================================================================	
	//Creating Canvas
	TCanvas* c5 = new TCanvas("c5","Canvas 5 - behavior of the Pion",1200,600);
	c5->Divide(2,2);
	c5->cd(1);
	Trkpipt_Histo->SetLineColor(kBlue);
	Trkpipt_Histo->SetMarkerStyle(7);
	Trkpipt_Histo->SetMarkerStyle(21);
	//Trkpipt_Histo->SetStats(0);
	Trkpipt_Histo->SetMarkerColor(kBlue);
	Trkpipt_Histo->SetFillColor(kBlue);
	Trkpipt_Histo->Draw("HIST");
	//TrkSpt_Histo->Draw("e1pSAME");
	TLegend* leg_Trkpipt = new TLegend(0.82,0.5,0.95,0.65);
   	leg_Trkpipt->SetFillColor(kWhite);
	leg_Trkpipt->SetFillStyle(1001);
	leg_Trkpipt->SetBorderSize(0);
	//leg_Trkpipt->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_Trkpipt->AddEntry(Trkpipt_Histo,"Trkpipt","L");
	leg_Trkpipt->Draw();
	//-------------------------------------------------------------------------
	c5->cd(2);
	Trkpinhits_Histo->SetLineColor(kBlue);
	Trkpinhits_Histo->SetMarkerStyle(7);
	Trkpinhits_Histo->SetMarkerStyle(21);
	//TrkSnhits_Histo->SetStats(0);
	Trkpinhits_Histo->SetMarkerColor(kBlue);
	Trkpinhits_Histo->SetFillColor(kBlue);
	Trkpinhits_Histo->Draw("HIST");
	//Trkpinhits_Histo->Draw("e1pSAME");
	TLegend* leg_Trkpinhits = new TLegend(0.82,0.5,0.95,0.65);
   	leg_Trkpinhits->SetFillColor(kWhite);
	leg_Trkpinhits->SetFillStyle(1001);
	leg_Trkpinhits->SetBorderSize(0);
	//leg_TrkSnhits->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_Trkpinhits->AddEntry(Trkpinhits_Histo,"Trkpinhits","L");
	leg_Trkpinhits->Draw();
	//-------------------------------------------------------------------------
	c5->cd(3);
	Trkpidxy_Histo->SetLineColor(kBlue);
	Trkpidxy_Histo->SetMarkerStyle(7);
	Trkpidxy_Histo->SetMarkerStyle(21);
	//Trkpidxy_Histo->SetStats(0);
	Trkpidxy_Histo->SetMarkerColor(kBlue);
	Trkpidxy_Histo->SetFillColor(kBlue);
	Trkpidxy_Histo->Draw("HIST");
	//Trkpidxy_Histo->Draw("e1pSAME");
	TLegend* leg_Trkpidxy = new TLegend(0.82,0.5,0.95,0.65);
   	leg_Trkpidxy->SetFillColor(kWhite);
	leg_Trkpidxy->SetFillStyle(1001);
	leg_Trkpidxy->SetBorderSize(0);
	//leg_Trkpidxy->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_Trkpidxy->AddEntry(Trkpidxy_Histo,"Trkpidxy","L");
	leg_Trkpidxy->Draw();
	//-------------------------------------------------------------------------
	c5->cd(4);
	Trkpidz_Histo->SetLineColor(kBlue);
	Trkpidz_Histo->SetMarkerStyle(7);
	Trkpidz_Histo->SetMarkerStyle(21);
	//Trkpidz_Histo->SetStats(0);
	Trkpidz_Histo->SetMarkerColor(kBlue);
	Trkpidz_Histo->SetFillColor(kBlue);
	Trkpidz_Histo->Draw("HIST");
	//Trkpidz_Histo->Draw("e1pSAME");
	TLegend* leg_Trkpidz = new TLegend(0.82,0.5,0.95,0.65);
   	leg_Trkpidz->SetFillColor(kWhite);
	leg_Trkpidz->SetFillStyle(1001);
	leg_Trkpidz->SetBorderSize(0);
	//leg_Trkpidz->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_Trkpidz->AddEntry(Trkpidz_Histo,"Trkpidz","L");
	leg_Trkpidz->Draw();
	c5->SaveAs("/eos/user/r/ragomesd/analysisB2019/Canvas5.png");
	//=========================================================================	
	//Creating Canvas
	TCanvas* c6 = new TCanvas("c6","Canvas 6 - behavior of the Pion",1200,600);
	c6->Divide(2,2);
	c6->cd(1);
	Trkpieta_Histo->SetLineColor(kBlue);
	Trkpieta_Histo->SetMarkerStyle(7);
	Trkpieta_Histo->SetMarkerStyle(21);
	//Trkpieta_Histo->SetStats(0);
	Trkpieta_Histo->SetMarkerColor(kBlue);
	Trkpieta_Histo->SetFillColor(kBlue);
	Trkpieta_Histo->Draw("HIST");
	//Trkpieta_Histo->Draw("e1pSAME");
	TLegend* leg_Trkpieta = new TLegend(0.82,0.5,0.95,0.65);
   	leg_Trkpieta->SetFillColor(kWhite);
	leg_Trkpieta->SetFillStyle(1001);
	leg_Trkpieta->SetBorderSize(0);
	//leg_Trkpieta->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_Trkpieta->AddEntry(Trkpieta_Histo,"Trkpieta","L");
	leg_Trkpieta->Draw();
	//-------------------------------------------------------------------------
	c6->cd(2);
	Trkpiphi_Histo->SetLineColor(kBlue);
	Trkpiphi_Histo->SetMarkerStyle(7);
	Trkpiphi_Histo->SetMarkerStyle(21);
	//TrkTrkpiphiSphi_Histo->SetStats(0);
	Trkpiphi_Histo->SetMarkerColor(kBlue);
	Trkpiphi_Histo->SetFillColor(kBlue);
	Trkpiphi_Histo->Draw("HIST");
	//Trkpiphi_Histo->Draw("e1pSAME");
	TLegend* leg_Trkpiphi = new TLegend(0.82,0.5,0.95,0.65);
   	leg_Trkpiphi->SetFillColor(kWhite);
	leg_Trkpiphi->SetFillStyle(1001);
	leg_Trkpiphi->SetBorderSize(0);
	//leg_Trkpiphi->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_Trkpiphi->AddEntry(Trkpiphi_Histo,"Trkpiphi","L");
	leg_Trkpiphi->Draw();
	//-------------------------------------------------------------------------
	c6->cd(3);
	Trkpichi2_Histo->SetLineColor(kBlue);
	Trkpichi2_Histo->SetMarkerStyle(7);
	Trkpichi2_Histo->SetMarkerStyle(21);
	//TrkTrkpichi2Sphi_Histo->SetStats(0);
	Trkpichi2_Histo->SetMarkerColor(kBlue);
	Trkpichi2_Histo->SetFillColor(kBlue);
	Trkpichi2_Histo->Draw("HIST");
	//Trkpichi2_Histo->Draw("e1pSAME");
	TLegend* leg_Trkpichi2 = new TLegend(0.82,0.5,0.95,0.65);
   	leg_Trkpichi2->SetFillColor(kWhite);
	leg_Trkpichi2->SetFillStyle(1001);
	leg_Trkpichi2->SetBorderSize(0);
	//leg_Trkpichi2->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_Trkpichi2->AddEntry(Trkpichi2_Histo,"Trkpichi2","L");
	leg_Trkpichi2->Draw();
	c6->SaveAs("/eos/user/r/ragomesd/analysisB2019/Canvas6.png");
	//=========================================================================	
	//Creating Canvas
	TCanvas* c7 = new TCanvas("c7","Canvas 7 - behavior of the Kaon",1200,600);
	c7->Divide(2,2);
	c7->cd(1);
	TrkKpt_Histo->SetLineColor(kBlue);
	TrkKpt_Histo->SetMarkerStyle(7);
	TrkKpt_Histo->SetMarkerStyle(21);
	//TrkKpt_Histo->SetStats(0);
	TrkKpt_Histo->SetMarkerColor(kBlue);
	TrkKpt_Histo->SetFillColor(kBlue);
	TrkKpt_Histo->Draw("HIST");
	//TrkSpt_Histo->Draw("e1pSAME");
	TLegend* leg_TrkKpt = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkKpt->SetFillColor(kWhite);
	leg_TrkKpt->SetFillStyle(1001);
	leg_TrkKpt->SetBorderSize(0);
	//leg_TrkKpt->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkKpt->AddEntry(TrkKpt_Histo,"TrkKpt","L");
	leg_TrkKpt->Draw();
	//-------------------------------------------------------------------------
	c7->cd(2);
	TrkKnhits_Histo->SetLineColor(kBlue);
	TrkKnhits_Histo->SetMarkerStyle(7);
	TrkKnhits_Histo->SetMarkerStyle(21);
	//TrkSnhits_Histo->SetStats(0);
	TrkKnhits_Histo->SetMarkerColor(kBlue);
	TrkKnhits_Histo->SetFillColor(kBlue);
	TrkKnhits_Histo->Draw("HIST");
	//TrkKnhits_Histo->Draw("e1pSAME");
	TLegend* leg_TrkKnhits = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkKnhits->SetFillColor(kWhite);
	leg_TrkKnhits->SetFillStyle(1001);
	leg_TrkKnhits->SetBorderSize(0);
	//leg_TrkSnhits->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkKnhits->AddEntry(TrkKnhits_Histo,"TrkKnhits","L");
	leg_TrkKnhits->Draw();
	//-------------------------------------------------------------------------
	c7->cd(3);
	TrkKdxy_Histo->SetLineColor(kBlue);
	TrkKdxy_Histo->SetMarkerStyle(7);
	TrkKdxy_Histo->SetMarkerStyle(21);
	//TrkKdxy_Histo->SetStats(0);
	TrkKdxy_Histo->SetMarkerColor(kBlue);
	TrkKdxy_Histo->SetFillColor(kBlue);
	TrkKdxy_Histo->Draw("HIST");
	//TrkKdxy_Histo->Draw("e1pSAME");
	TLegend* leg_TrkKdxy = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkKdxy->SetFillColor(kWhite);
	leg_TrkKdxy->SetFillStyle(1001);
	leg_TrkKdxy->SetBorderSize(0);
	//leg_TrkKdxy->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkKdxy->AddEntry(TrkKdxy_Histo,"TrkKdxy","L");
	leg_TrkKdxy->Draw();
	//-------------------------------------------------------------------------
	c7->cd(4);
	TrkKdz_Histo->SetLineColor(kBlue);
	TrkKdz_Histo->SetMarkerStyle(7);
	TrkKdz_Histo->SetMarkerStyle(21);
	//TrkKdz_Histo->SetStats(0);
	TrkKdz_Histo->SetMarkerColor(kBlue);
	TrkKdz_Histo->SetFillColor(kBlue);
	TrkKdz_Histo->Draw("HIST");
	//TrkKdz_Histo->Draw("e1pSAME");
	TLegend* leg_TrkKdz = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkKdz->SetFillColor(kWhite);
	leg_TrkKdz->SetFillStyle(1001);
	leg_TrkKdz->SetBorderSize(0);
	//leg_TrkKdz->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkKdz->AddEntry(TrkKdz_Histo,"TrkKdz","L");
	leg_TrkKdz->Draw();
	c7->SaveAs("/eos/user/r/ragomesd/analysisB2019/Canvas7.png");
	//=========================================================================	
	//Creating Canvas
	TCanvas* c8 = new TCanvas("c8","Canvas 8 - behavior of the Kaon",1200,600);
	c8->Divide(2,2);
	c8->cd(1);
	TrkKeta_Histo->SetLineColor(kBlue);
	TrkKeta_Histo->SetMarkerStyle(7);
	TrkKeta_Histo->SetMarkerStyle(21);
	//TrkKeta_Histo->SetStats(0);
	TrkKeta_Histo->SetMarkerColor(kBlue);
	TrkKeta_Histo->SetFillColor(kBlue);
	TrkKeta_Histo->Draw("HIST");
	//TrkKeta_Histo->Draw("e1pSAME");
	TLegend* leg_TrkKeta = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkKeta->SetFillColor(kWhite);
	leg_TrkKeta->SetFillStyle(1001);
	leg_TrkKeta->SetBorderSize(0);
	//leg_TrkKeta->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkKeta->AddEntry(TrkKeta_Histo,"TrkKeta","L");
	leg_TrkKeta->Draw();
	//-------------------------------------------------------------------------
	c8->cd(2);
	TrkKphi_Histo->SetLineColor(kBlue);
	TrkKphi_Histo->SetMarkerStyle(7);
	TrkKphi_Histo->SetMarkerStyle(21);
	//TrkTrkKphiSphi_Histo->SetStats(0);
	TrkKphi_Histo->SetMarkerColor(kBlue);
	TrkKphi_Histo->SetFillColor(kBlue);
	TrkKphi_Histo->Draw("HIST");
	//TrkKphi_Histo->Draw("e1pSAME");
	TLegend* leg_TrkKphi = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkKphi->SetFillColor(kWhite);
	leg_TrkKphi->SetFillStyle(1001);
	leg_TrkKphi->SetBorderSize(0);
	//leg_TrkKphi->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkKphi->AddEntry(TrkKphi_Histo,"TrkKphi","L");
	leg_TrkKphi->Draw();
	//-------------------------------------------------------------------------
	c8->cd(3);
	TrkKchi2_Histo->SetLineColor(kBlue);
	TrkKchi2_Histo->SetMarkerStyle(7);
	TrkKchi2_Histo->SetMarkerStyle(21);
	//TrkTrkKchi2Sphi_Histo->SetStats(0);
	TrkKchi2_Histo->SetMarkerColor(kBlue);
	TrkKchi2_Histo->SetFillColor(kBlue);
	TrkKchi2_Histo->Draw("HIST");
	//TrkKchi2_Histo->Draw("e1pSAME");
	TLegend* leg_TrkKchi2 = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkKchi2->SetFillColor(kWhite);
	leg_TrkKchi2->SetFillStyle(1001);
	leg_TrkKchi2->SetBorderSize(0);
	//leg_TrkKchi2->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkKchi2->AddEntry(TrkKchi2_Histo,"TrkKchi2","L");
	leg_TrkKchi2->Draw();
	c8->SaveAs("/eos/user/r/ragomesd/analysisB2019/Canvas8.png");
}//end program

