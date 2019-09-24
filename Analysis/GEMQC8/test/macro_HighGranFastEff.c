#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TTree.h>
#include "TGraphAsymmErrors.h"
#include <TBranch.h>
#include <TCanvas.h>
#include <TSpectrum.h>
#include <TF1.h>
#include "TLine.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <TLatex.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

void macro_HighGranFastEff(int run, string configDir)
{
	// Setting variables for min and max displayed efficiency (to be tuned in the analysis if things go wrong...)
	const float min_eff = 0.00;
	const float max_eff = 1.00;

  // Getting the root file

  string filename = "HighGranFastEff_out_run_";
  for (unsigned int i=0; i<(6-to_string(run).size()); i++)
  {
    filename += "0";
  }
  filename += to_string(run) + ".root";
  const char *file = filename.c_str();

  TFile *infile = new TFile(file,"UPDATE");
  if (infile->IsOpen()) cout << "File opened successfully" << endl;

  // Getting the tree

  TTree *tree = (TTree*)infile->Get("HighGranFastEffQC8/tree");
	int nTotEvents = tree->GetEntries();

	// Declaration of name variables

  char *name = new char[40];
  string namename = "";

	// Open stand configuration file for present run & get names + positions of chambers
	string configName = configDir + "StandGeometryConfiguration_run" + to_string(run) + ".csv";
	ifstream standConfigFile (configName);

	string line, split, comma = ",", slash = "/";
	vector<string> chamberName;
	int ChPos = 0;
	vector<int> chamberPos;
	size_t pos = 0;

	if (standConfigFile.is_open())
	{
		while (getline(standConfigFile, line))
		{
			pos = line.find(comma);
			split = line.substr(0, pos);
			if (split == "CH_SERIAL_NUMBER") continue;
			chamberName.push_back(split);
			line.erase(0, pos + comma.length());

			pos = line.find(comma);
			line.erase(0, pos + comma.length());

			pos = line.find(slash);
			split = line.substr(0, pos);
			ChPos = (stoi(split)-1)*2; // (Row-1)*2
			line.erase(0, pos + slash.length());

			pos = line.find(slash);
			split = line.substr(0, pos);
			ChPos += (stoi(split)-1)*10; // (Row-1)*2 + (Col-1)*10
			line.erase(0, pos + slash.length());

			pos = line.find(comma);
			split = line.substr(0, pos);
			if (split == "B") ChPos += 0; // (Row-1)*2 + (Col-1)*10 + 0
			if (split == "T") ChPos += 1; // (Row-1)*2 + (Col-1)*10 + 1
			line.erase(0, pos + comma.length());

			chamberPos.push_back(ChPos);
		}
	}
	else cout << "Error opening file: " << configName << endl;

  //************************ Efficiency by blocks of strips (HGFE) **********************//
	//Get the 3D block histograms
	TH3D *numblocks3D = (TH3D*)infile->Get("HighGranFastEffQC8/num_block");
	TH3D *denomblocks3D = (TH3D*)infile->Get("HighGranFastEffQC8/denom_block");

	//Get the 3D x position block histograms 
	TH3D *numblocks3Dx = (TH3D*)infile->Get("HighGranFastEffQC8/num_block_x");
	TH3D *denomblocks3Dx = (TH3D*)infile->Get("HighGranFastEffQC8/denom_block_x");

	//1D histogram: efficiency per chamber
	TH1D* numblocks1D = new TH1D("numblocks1D","",30,0,30);
	TH1D* denomblocks1D = new TH1D("denomblocks1D","",30,0,30);

	//2D histogram: efficiency per blocks per chamber
	TH2D *Eff2DperBlockperCh[30];
	TH2D *Eff2DxperBlockperCh[30];
	TH2D *Eff2DperVFATperCh[30];
	TH2D *Eff2DxperVFATperCh[30];

	int nstrips = 384;	
	int stripsperblock = 8;
	int nblocks = 0;
	
	nblocks = nstrips / stripsperblock;

	//Overall efficiency per chamber
	for (int ch=0; ch<30; ch++)
	{
	    double total_denom = denomblocks3D->Integral(1,48,1,9,ch+1,ch+1);
	    double total_num   = numblocks3D->Integral(1,48,1,9,ch+1,ch+1);
	    denomblocks1D->SetBinContent(ch+1, total_denom);
	    numblocks1D->SetBinContent(ch+1, total_num);
	}  

	TGraphAsymmErrors *effblocks1Dch = new TGraphAsymmErrors;
	effblocks1Dch->Divide(numblocks1D,denomblocks1D);

	//efficiency per block
	for (int ch=0; ch<30; ch++)
	{
		sprintf(name,"Eff2DperBlockfor_Ch_%u",ch);
		Eff2DperBlockperCh[ch] = new TH2D(name,"",nblocks,0,nblocks,10,0,10);

		for (int b=0; b<nblocks; b++)
		{
			for (int eta=0; eta<10; eta++)
			{

			 double num   = numblocks3D->GetBinContent(b+1,eta+1,ch+1);
 			 double denom = denomblocks3D->GetBinContent(b+1,eta+1,ch+1);
			 double eff = 100. * num / denom;
			 Eff2DperBlockperCh[ch]->SetBinContent(b+1,eta+1,eff);
			}
		}
	}

	//efficiency per block / Trapezoid
	for (int ch=0; ch<30; ch++)
	{
		sprintf(name,"Eff2DxperBlockfor_Ch_%u",ch);
		Eff2DxperBlockperCh[ch] = new TH2D(name,"",500,-50,50,10,0,10);

		for (int b=0; b<500; b++)
		{
			for (int eta=0; eta<10; eta++)
			{
			 	double num   = numblocks3Dx->GetBinContent(b+1,eta+1,ch+1);
 			 	double denom = denomblocks3Dx->GetBinContent(b+1,eta+1,ch+1);
			 	double eff = 100. * num / denom;
				Eff2DxperBlockperCh[ch]->SetBinContent(b+1,eta+1,eff);
			}
		}
	}

	//efficiency per VFAT vs eta
	double vfat_num[3][10];
	double vfat_denom[3][10];

	for (int ch=0; ch<30; ch++)
	{
		sprintf(name,"Eff2DperVFATfor_Ch_%u",ch);
		Eff2DperVFATperCh[ch] = new TH2D(name,"",3,0,3,10,0,10);

		for (int v=0; v<3; v++)
		{
			for (int eta=0; eta<10; eta++)
			{
				vfat_denom[v][eta] = denomblocks3D->Integral(1+16*v,16*(v+1),eta+1,eta+1,ch+1,ch+1);
	   			vfat_num[v][eta]   = numblocks3D->Integral(1+16*v,16*(v+1),eta+1,eta+1,ch+1,ch+1);
			 	double eff = 100. * vfat_num[v][eta] / vfat_denom[v][eta];
				Eff2DperVFATperCh[ch]->SetBinContent(v+1,eta+1,eff);
			}
		}
	}

  // Results for the 30 chambers
  TCanvas *Canvas = new TCanvas("Canvas","Canvas",0,0,800,600);
  TF1 *targetLine = new TF1("targetLine","0.90",0,30);
  targetLine->SetLineColor(kRed);

  ofstream outfile;

	// Plot overall efficiency 1D by integrating blocks
	namename = "eff_per_block_1D_Integral_all_Chambers__HGFE_" + to_string(run);
	effblocks1Dch->SetTitle(namename.c_str());
	effblocks1Dch->GetXaxis()->SetTitle("Chamber position");
	effblocks1Dch->GetYaxis()->SetTitle("Efficiency");
	effblocks1Dch->GetYaxis()->SetRangeUser(min_eff,max_eff);
	effblocks1Dch->SetMarkerStyle(20);
	effblocks1Dch->Draw();
	effblocks1Dch->Write(namename.c_str());
	targetLine->Draw("SAME");
	namename = "efficiency_per_block_1D_Integral_all_Chambers_HGFE_run_" + to_string(run) + ".png";
	Canvas->SaveAs(namename.c_str());
	Canvas->Clear();

	// Plot Numerator vs blocks vs chamber
	namename = "Numerator_Xpos_blocks_3D_run_" + to_string(run);
	numblocks3Dx->SetTitle(namename.c_str());
	numblocks3Dx->GetXaxis()->SetTitle("X position [cm]");
	numblocks3Dx->GetYaxis()->SetTitle("Eta");
	numblocks3Dx->GetZaxis()->SetTitle("Chamber");
	numblocks3Dx->Write(namename.c_str());
	numblocks3Dx->SetLineColor(kRed);
	numblocks3Dx->Draw();
	numblocks3Dx->SetTitle(namename.c_str());
	namename = "Numerator_Xpos_blocks_3D_HGFE_run_" + to_string(run) + ".png";
	Canvas->SaveAs(namename.c_str());
	Canvas->Clear();

	// Plot x pos denominator vs blocks vs chamber
	namename = "Denominator_Xpos_blocks_3D_run_" + to_string(run);
	denomblocks3Dx->SetTitle(namename.c_str());
	denomblocks3Dx->GetXaxis()->SetTitle("X position [cm]");
	denomblocks3Dx->GetYaxis()->SetTitle("Eta");
	denomblocks3Dx->GetZaxis()->SetTitle("Chamber");
	denomblocks3Dx->Write(namename.c_str());
	denomblocks3Dx->SetLineColor(kRed);
	denomblocks3Dx->Draw();
	denomblocks3Dx->SetTitle(namename.c_str());
	namename = "Denominator_blocks_3D_HGFE_run_" + to_string(run) + ".png";
	Canvas->SaveAs(namename.c_str());
	Canvas->Clear();

	// Plot efficiency vs x pos of 4th strip of blocks of 8 strips for each chamber
	for (unsigned int i=0; i<chamberPos.size(); i++)
	 {
	int c = chamberPos[i];

	namename = "Efficiency_blocks_Trapezoid_HGFE_" + chamberName[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
	Eff2DxperBlockperCh[c]->SetTitle(namename.c_str());
	Eff2DxperBlockperCh[c]->GetXaxis()->SetTitle("X position [cm]");
	Eff2DxperBlockperCh[c]->GetYaxis()->SetTitle("Eta");
	Eff2DxperBlockperCh[c]->SetStats(false);
	Eff2DxperBlockperCh[c]->Draw("colz");
	Eff2DxperBlockperCh[c]->Write(namename.c_str());
	namename = "outPlots_Chamber_Pos_" + to_string(chamberPos[i]) + "/Efficiency_blocks_Trapezoid_HGFE_" + to_string(chamberPos[i]) + "_run_" + to_string(run) + ".png";

	Canvas->SaveAs(namename.c_str());
	Canvas->Clear();
	}

	// Plot efficiency per block of 8 strips for each chamber
	for (unsigned int i=0; i<chamberPos.size(); i++)
	 {
	int c = chamberPos[i];

	namename = "Efficiency_blocks_HGFE_" + chamberName[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
	Eff2DperBlockperCh[c]->SetTitle(namename.c_str());
	Eff2DperBlockperCh[c]->GetXaxis()->SetTitle("Block");
	Eff2DperBlockperCh[c]->GetYaxis()->SetTitle("Eta");
	Eff2DperBlockperCh[c]->SetStats(false);
	Eff2DperBlockperCh[c]->Draw("colz");
	Eff2DperBlockperCh[c]->Write(namename.c_str());
	namename = "outPlots_Chamber_Pos_" + to_string(chamberPos[i]) + "/Efficiency_blocks_HGFE_" + to_string(chamberPos[i]) + "_run_" + to_string(run) + ".png";

	Canvas->SaveAs(namename.c_str());
	Canvas->Clear();
	}

	// Plot efficiency per VFAT for each chamber
 	for (unsigned int i=0; i<chamberPos.size(); i++)
 	{
	int c = chamberPos[i];

	namename = "Efficiency_VFAT_HGFE_" + chamberName[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
	Eff2DperVFATperCh[c]->SetTitle(namename.c_str());
	Eff2DperVFATperCh[c]->GetXaxis()->SetTitle("VFAT");
	Eff2DperVFATperCh[c]->GetYaxis()->SetTitle("Eta");
	Eff2DperVFATperCh[c]->SetStats(false);
	Eff2DperVFATperCh[c]->Draw("colz,TEXT0");
	Eff2DperVFATperCh[c]->Write(namename.c_str());
	namename = "outPlots_Chamber_Pos_" + to_string(chamberPos[i]) + "/Efficiency_VFAT_HGFE_" + to_string(chamberPos[i]) + "_run_" + to_string(run) + ".png";

	Canvas->SaveAs(namename.c_str());
	Canvas->Clear();
	}

  infile->Close();

}
