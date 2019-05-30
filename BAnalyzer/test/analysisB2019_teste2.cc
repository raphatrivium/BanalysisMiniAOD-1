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
//#ifndef ROOT_TLatex
#define ROOT_TLatex
//#ifndef ROOTiosfwd
#include "Riosfwd.h"
//#endif
//#ifndef ROOT_TText
#include "TText.h"
//#endif
//#ifndef ROOT_TAttLine
#include "TAttLine.h"
//#endif
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

int analysisB2019_teste2()
{
				
	//Variables	and Vectors
	int Total_Events = 0;

	int TotalTracks = 0;
	int TracksAfterTrigger = 0;
	int TracksHasTrackDetails = 0;
	int TracksChargeZero = 0;
	int TracksHighPurity = 0;
	int TracksPDG211 = 0;
	int TracksPtZeroFive = 0;
	int TracksNumberOfHits2 = 0;
	int TracksDxyThree = 0;
	int TracksDzThree = 0;
	int TrackSlowPionCandidates = 0;
	int	Observation = 0;
	int TracksPtZeroSix = 0;
	int TracksChiTwoFive = 0;
	int TracksNumberOfHits5 = 0;
	int TracksNumberOfPixelHits2 = 0;
	int TracksDzOne = 0;
	int TracksDxyZeroOne = 0;
	int TrackKaonPionCandidates = 0;

	int D0AfterLorentzVector = 0;
	int DsMinusD0Zerothree = 0;
	int TransientTrackOfpiK = 0;
	int PointingcosPhi = 0;
	int Significance = 0;
	int D0pTThree = 0;
	int D0Candidates = 0;
	int DsAfterLorentzVector = 0;
	int DsMinusD0 = 0;
	int DsCandidates = 0;
	
	//counters for file f1 (pythia)
	int countTotalEvents = 0;

	int countTotalTracks = 0;
	int countTracksAfterTrigger = 0;
	int countTracksHasTrackDetails = 0;
	int countTracksChargeZero = 0;
	int countTracksHighPurity = 0;
	int countTracksPDG211 = 0;
	int countTracksPtZeroFive = 0;
	int countTracksNumberOfHits2 = 0;
	int countTracksDxyThree = 0;
	int countTracksDzThree = 0;
	int countTrackSlowPionCandidates = 0;
	int countObservation = 0;
	int countTracksPtZeroSix = 0;
	int countTracksChiTwoFive = 0;
	int countTracksNumberOfHits5 = 0;
	int countTracksNumberOfPixelHits2 = 0;
	int countTracksDzOne = 0;
	int countTracksDxyZeroOne = 0;
	int countTrackKaonPionCandidates = 0;

	int countD0AfterLorentzVector = 0;
	int countDsMinusD0Zerothree = 0;
	int countTransientTrackOfpiK = 0;
	int countPointingcosPhi = 0;
	int countSignificance = 0;
	int countD0pTThree = 0;
	int countD0Candidates = 0;
	int countDsAfterLorentzVector = 0;
	int countDsMinusD0 = 0;
	int countDsCandidates = 0;

	//double M = 0.;
	//double Pt = 0.;
	//double Eta = 0.;
	//double Rapidity = 0.;

	//-------Reading the root file and the tree-------------------------------------	
	TFile *f1 = new TFile("D0DstarData_allfiles6.root");
	TTree *t1 = (TTree*)f1->Get("analysis/data");

	
	//=====================================================================		
	//counters
	TBranch *b_Total_Events = t1->GetBranch("Total_Events");
	b_Total_Events->SetAddress(&Total_Events);
	TBranch *b_TotalTracks = t1->GetBranch("TotalTracks");
	b_TotalTracks->SetAddress(&TotalTracks);
	TBranch *b_TracksAfterTrigger = t1->GetBranch("TracksAfterTrigger");
	b_TracksAfterTrigger->SetAddress(&TracksAfterTrigger);
	TBranch *b_TracksHasTrackDetails = t1->GetBranch("TracksHasTrackDetails");
	b_TracksHasTrackDetails->SetAddress(&TracksHasTrackDetails);
	TBranch *b_TracksChargeZero = t1->GetBranch("TracksChargeZero");
	b_TracksChargeZero->SetAddress(&TracksChargeZero);
	TBranch *b_TracksHighPurity = t1->GetBranch("TracksHighPurity");
	b_TracksHighPurity->SetAddress(&TracksHighPurity);
	TBranch *b_TracksPDG211 = t1->GetBranch("TracksPDG211");
	b_TracksPDG211->SetAddress(&TracksPDG211);
	TBranch *b_TracksPtZeroFive = t1->GetBranch("TracksPtZeroFive");
	b_TracksPtZeroFive->SetAddress(&TracksPtZeroFive);
	TBranch *b_TracksNumberOfHits2 = t1->GetBranch("TracksNumberOfHits2");
	b_TracksNumberOfHits2->SetAddress(&TracksNumberOfHits2);
	TBranch *b_TracksDxyThree = t1->GetBranch("TracksDxyThree");
	b_TracksDxyThree->SetAddress(&TracksDxyThree);
	TBranch *b_TracksDzThree = t1->GetBranch("TracksDzThree");
	b_TracksDzThree->SetAddress(&TracksDzThree);
	TBranch *b_TrackSlowPionCandidates = t1->GetBranch("TrackSlowPionCandidates");
	b_TrackSlowPionCandidates->SetAddress(&TrackSlowPionCandidates);
	TBranch *b_Observation = t1->GetBranch("Observation");
	b_Observation->SetAddress(&Observation);
	TBranch *b_TracksPtZeroSix = t1->GetBranch("TracksPtZeroSix");
	b_TracksPtZeroSix->SetAddress(&TracksPtZeroSix);
	TBranch *b_TracksChiTwoFive = t1->GetBranch("TracksChiTwoFive");
	b_TracksChiTwoFive->SetAddress(&TracksChiTwoFive);
	TBranch *b_TracksNumberOfHits5 = t1->GetBranch("TracksNumberOfHits5");
	b_TracksNumberOfHits5->SetAddress(&TracksNumberOfHits5);
	TBranch *b_TracksNumberOfPixelHits2 = t1->GetBranch("TracksNumberOfPixelHits2");
	b_TracksNumberOfPixelHits2->SetAddress(&TracksNumberOfPixelHits2);
	TBranch *b_TracksDzOne = t1->GetBranch("TracksDzOne");
	b_TracksDzOne->SetAddress(&TracksDzOne);
	TBranch *b_TracksDxyZeroOne = t1->GetBranch("TracksDxyZeroOne");
	b_TracksDxyZeroOne->SetAddress(&TracksDxyZeroOne);
	TBranch *b_TrackKaonPionCandidates = t1->GetBranch("TrackKaonPionCandidates");
	b_TrackKaonPionCandidates->SetAddress(&TrackKaonPionCandidates);

	TBranch *b_D0AfterLorentzVector = t1->GetBranch("D0AfterLorentzVector");
	b_D0AfterLorentzVector->SetAddress(&D0AfterLorentzVector);
	TBranch *b_DsMinusD0Zerothree = t1->GetBranch("DsMinusD0Zerothree");
	b_DsMinusD0Zerothree->SetAddress(&DsMinusD0Zerothree); 
	TBranch *b_TransientTrackOfpiK = t1->GetBranch("TransientTrackOfpiK");
	b_TransientTrackOfpiK->SetAddress(&TransientTrackOfpiK);
	TBranch *b_PointingcosPhi = t1->GetBranch("PointingcosPhi");
	b_PointingcosPhi->SetAddress(&PointingcosPhi);
	TBranch *b_Significance = t1->GetBranch("Significance");
	b_Significance->SetAddress(&Significance);
	TBranch *b_D0pTThree = t1->GetBranch("D0pTThree");
	b_D0pTThree->SetAddress(&D0pTThree);
	TBranch *b_D0Candidates = t1->GetBranch("D0Candidates");
	b_D0Candidates->SetAddress(&D0Candidates);
	TBranch *b_DsAfterLorentzVector = t1->GetBranch("DsAfterLorentzVector");
	b_DsAfterLorentzVector->SetAddress(&DsAfterLorentzVector);
	TBranch *b_DsMinusD0 = t1->GetBranch("DsMinusD0");
	b_DsMinusD0->SetAddress(&DsMinusD0);
	TBranch *b_DsCandidates = t1->GetBranch("DsCandidates");
	b_DsCandidates->SetAddress(&DsCandidates);
	

	//**********************************************************		
	//Reading Number of tree entries for file f1
	Long64_t nentries = t1->GetEntries();
	cout<< "Number of tree entries: "<< nentries <<std::endl;

	//Long64_t GetEntriesFast = t1->GetEntriesFast();
	//cout<< "GetEntriesFast: "<< GetEntriesFast <<std::endl;

	//Long64_t nbytes = 0, nb = 0, i=0;
	for (Long64_t jentry=0; jentry < nentries; jentry++) //loop tree entries for file f1
	{
      	Long64_t ientry = t1->LoadTree(jentry);
      	//std::cout << "nentries " << nentries << " ientry " << ientry << " jentry " << jentry <<std::endl;
      	if (ientry < 0) break;

		//if ((jentry/nentries) % (nentries*0.01) == 0) cout <<"*****"<< (jentry*100)/nentries << "per cent done***" << endl;
		//{			
		//}
		double percent = (jentry*100)/nentries;
		if (jentry % 50000 == 0) cout <<"*****"<< percent << "per cent done***" << endl;
		//cout <<"*****entrada"<< jentry << "***" << endl;

		//counters
		b_Total_Events->GetEntry(ientry);				
		b_TotalTracks->GetEntry(ientry);
		b_TracksAfterTrigger->GetEntry(ientry);
		b_TracksHasTrackDetails->GetEntry(ientry);
		b_TracksChargeZero->GetEntry(ientry);
		b_TracksHighPurity->GetEntry(ientry);
		b_TracksPDG211->GetEntry(ientry);
		b_TracksPtZeroFive->GetEntry(ientry);
		b_TracksNumberOfHits2->GetEntry(ientry);
		b_TracksDxyThree->GetEntry(ientry);
		b_TracksDzThree->GetEntry(ientry);
		b_TrackSlowPionCandidates->GetEntry(ientry);
		b_TracksPtZeroSix->GetEntry(ientry); 
		b_Observation->GetEntry(ientry);
		b_TracksChiTwoFive->GetEntry(ientry);
		b_TracksNumberOfHits5->GetEntry(ientry);
		b_TracksNumberOfPixelHits2->GetEntry(ientry);
		b_TracksDxyZeroOne->GetEntry(ientry);
		b_TracksDzOne->GetEntry(ientry);
		b_TrackKaonPionCandidates->GetEntry(ientry);

		b_D0AfterLorentzVector->GetEntry(ientry);
		b_DsMinusD0Zerothree->GetEntry(ientry);
		b_TransientTrackOfpiK->GetEntry(ientry);
		b_PointingcosPhi->GetEntry(ientry);
		b_Significance->GetEntry(ientry);
		b_D0pTThree->GetEntry(ientry);
		b_D0Candidates->GetEntry(ientry);
		b_DsMinusD0->GetEntry(ientry);
		b_DsCandidates->GetEntry(ientry);
		b_DsAfterLorentzVector->GetEntry(ientry);

		countTotalEvents += Total_Events;
		countTotalTracks += TotalTracks;
		countTracksAfterTrigger += TracksAfterTrigger;
		countTracksHasTrackDetails += TracksHasTrackDetails;
		countTracksChargeZero += TracksChargeZero;
		countTracksHighPurity += TracksHighPurity;
		countTracksPDG211 += TracksPDG211;
		countTracksPtZeroFive += TracksPtZeroFive;
		countTracksNumberOfHits2 += TracksNumberOfHits2;
		countTracksDxyThree += TracksDxyThree;
		countTracksDzThree += TracksDzThree;
		countTrackSlowPionCandidates += TrackSlowPionCandidates;
		countTracksPtZeroSix += TracksPtZeroSix;
		countObservation += Observation;
		countTracksChiTwoFive += TracksChiTwoFive;
		countTracksNumberOfHits5 += TracksNumberOfHits5;
		countTracksNumberOfPixelHits2 += TracksNumberOfPixelHits2;
		countTracksDxyZeroOne += TracksDxyZeroOne;
		countTracksDzOne += TracksDzOne;
		countTrackKaonPionCandidates += TrackKaonPionCandidates;

		countD0AfterLorentzVector += D0AfterLorentzVector;
		countDsMinusD0Zerothree += DsMinusD0Zerothree;
		countTransientTrackOfpiK += TransientTrackOfpiK;
		countPointingcosPhi += PointingcosPhi;
		countSignificance += Significance;
		countD0pTThree += D0pTThree;
		countD0Candidates += D0Candidates;
		countDsAfterLorentzVector += DsAfterLorentzVector;
		countDsMinusD0 += DsMinusD0;
		countDsCandidates += DsCandidates;
			
	}//End loop tree entries for file f1
	
	//TFile f_analysis("Jpsi.root","recreate"); //Creates root file
	//TTree t_analysis("t_analysis","analise_Tree"); //Creates a Tree
	//t_analysis.Branch("vectorInvariantMassJpsi",&vectorInvariantMassJpsi); //Creates a branch
	//t_analysis.Fill();
	//h_DimuonsOppositeChargeEtaJpsi_M2_Double.Write(); //Write() if a file is open, this function writes a root objectics on it.
	//t_analysis.Write();  //Write in the root file


	cout << "countTotalEvents: "<< countTotalEvents << endl;
	cout << "        " << endl;	
	cout << "=============Tracks=========================================== " << endl;
	cout << "countTotalTracks: "<< countTotalTracks << endl;
	cout << "countTracksAfterTrigger: "<< countTracksAfterTrigger << endl;
	cout << "countTracksHasTrackDetails: "<< countTracksHasTrackDetails << endl;
	cout << "countTracksChargeZero: "<< countTracksChargeZero << endl;
	cout << "countTracksHighPurity: "<< countTracksHighPurity << endl;
	cout << "countTracksPDG211: "<< countTracksPDG211 << endl;
	cout << "=============Criteria to SlowPion Tracks=======================" << endl;
	cout << "countTracksPtZeroFive: "<< countTracksPtZeroFive << endl;
	cout << "countTracksNumberOfHits2: "<< countTracksNumberOfHits2 << endl;
	cout << "countTracksDxyThree: "<< countTracksDxyThree << endl;
	cout << "countTracksDzThree: "<< countTracksDzThree << endl;
	cout << "countTrackSlowPionCandidates: "<< countTrackSlowPionCandidates << endl;
	cout << "Observation: "<< Observation << endl;
	cout << "=============Criteria to Pions and Kaons Tracks=======================" << endl;
	cout << "countTracksPtZeroSix: "<< countTracksPtZeroSix << endl;
	cout << "countTracksChiTwoFive: "<< countTracksChiTwoFive << endl;
	cout << "countTracksNumberOfHits5: "<< countTracksNumberOfHits5 << endl;
	cout << "countTracksNumberOfPixelHits2: "<< countTracksNumberOfPixelHits2 << endl;
	cout << "countTracksDxyZeroOne: "<< countTracksDxyZeroOne << endl;
	cout << "countTracksDzOne: "<< countTracksDzOne << endl;
	cout << "countTrackKaonPionCandidates: "<< countTrackKaonPionCandidates << endl;
	cout << "=============Criteria to D0 and D* =======================" << endl;
	cout << "countD0AfterLorentzVector: "<< countD0AfterLorentzVector << endl;
	cout << "countDsMinusD0Zerothree: "<< countDsMinusD0Zerothree << endl;
	cout << "countTransientTrackOfpiK: "<< countTransientTrackOfpiK << endl;
	cout << "countPointingcosPhi: "<< countPointingcosPhi << endl;
	cout << "countSignificance: "<< countSignificance << endl;
	cout << "countD0pTThree: "<< countD0pTThree << endl;
	cout << "countD0Candidates: "<< countD0Candidates << endl;
	cout << "countDsAfterLorentzVector: "<< countDsAfterLorentzVector << endl;
	cout << "countDsMinusD0: "<< countDsMinusD0 << endl;
	cout << "countDsCandidates: "<< countDsCandidates << endl;
	
	return 0;
}//end program

