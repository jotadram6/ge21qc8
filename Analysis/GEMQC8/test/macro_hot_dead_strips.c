#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include "TGraphAsymmErrors.h"
#include <TBranch.h>
#include <TCanvas.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <TLatex.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

void macro_hot_dead_strips(int run, string configDir)
{
  // Getting the root file

  string filename = "hot_dead_strips_out_run_";
  for (unsigned int i=0; i<(6-to_string(run).size()); i++)
  {
    filename += "0";
  }
  filename += to_string(run) + ".root";
  const char *file = filename.c_str();

  TFile *infile = new TFile(file,"UPDATE");
  if (infile->IsOpen()) cout << "File opened successfully" << endl;

  // Getting the tree

  TTree *tree = (TTree*)infile->Get("HotDeadStripsQC8/tree");

  // Getting 3D digi histogram

  TH3D *digi3D = (TH3D*)infile->Get("HotDeadStripsQC8/digiStrips");

  // Digi plots per chamber

  char *name = new char[40];
  string namename = "";

  long unsigned int max_digi_occupancy[30];
  long unsigned int digi_occupancy[30];

  TH2D *digi2D[30];

  for (int ch=0; ch<30; ch++)
  {
  	max_digi_occupancy[ch] = 0;
  	digi_occupancy[ch] = 0;

    sprintf(name,"Digi_ch_%u",ch);
    digi2D[ch] = new TH2D(name,"",384,0,384,8,-0.5,7.5);

    for (int eta=0; eta<8; eta++)
    {
      for (int strip=0; strip<384; strip++)
      {
        digi2D[ch]->SetBinContent((strip+1),(eta+1),digi3D->GetBinContent(strip+1,eta+1,ch+1));
        digi_occupancy[ch] = digi3D->GetBinContent(strip+1,eta+1,ch+1);
        if (digi_occupancy[ch] > max_digi_occupancy[ch]) max_digi_occupancy[ch] = digi_occupancy[ch];
      }
    }

    max_digi_occupancy[ch]++; // Just not to have histogram with 0 bins giving errors... Since it is only an upper limit on the histos, +1 doesn't really matter...
  }

  // Plot of digis per strip for each of the chambers

	TH1D *digisPerStripPerCh[30];

  for (int ch=0; ch<30; ch++)
  {
    sprintf(name,"DigisPerStrip_ch_%u",ch);
    digisPerStripPerCh[ch] = new TH1D(name,"",int(max_digi_occupancy[ch]/7),0,max_digi_occupancy[ch]);

    for (int eta=0; eta<8; eta++)
    {
      for (int strip=0; strip<384; strip++)
      {
        digisPerStripPerCh[ch]->Fill(digi3D->GetBinContent(strip+1,eta+1,ch+1));
      }
    }
  }

  // Open stand configuration file for present run & get names + positions of chambers

	string configName = configDir + "StandGeometryConfiguration_run" + to_string(run) + ".csv";
	ifstream standConfigFile (configName);

  string line, split, comma = ",", slash = "/";
	vector<string> chamberName, chamberNamePlot;
	int ChPos = 0;
	vector<int> chamberPos;
	size_t pos = 0, pos_slash = 0;

	if (standConfigFile.is_open())
	{
		while (getline(standConfigFile, line))
		{
			pos = line.find(comma);
			split = line.substr(0, pos);
			if (split == "CH_SERIAL_NUMBER") continue;
			chamberName.push_back(split);
			pos_slash = split.find(slash);
			chamberNamePlot.push_back(split.substr(0,pos_slash)+split.substr(pos_slash+slash.length(),pos));
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

  // Identification of value for being a dead or hot strip per chamber

  long int DeadStripLimitValue[30];
  long int HotStripLimitValue[30];
  const int nsigmas = 5;

  for (int ch=0; ch<30; ch++)
	{
		DeadStripLimitValue[ch] = 0;
		HotStripLimitValue[ch] = max_digi_occupancy[ch];
	}

  TCanvas *Canvas = new TCanvas("Canvas","Canvas",0,0,1000,800);

  for (unsigned int i=0; i<chamberPos.size(); i++)
  {
    int c = chamberPos[i];

    namename = "Digi_PerStrip_PerCh_" + chamberName[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
    digisPerStripPerCh[c]->SetTitle(namename.c_str());
    digisPerStripPerCh[c]->GetXaxis()->SetTitle("digisPerStripPerCh");
    digisPerStripPerCh[c]->GetYaxis()->SetTitle("Counts");
    digisPerStripPerCh[c]->Draw();

    digisPerStripPerCh[c]->GetXaxis()->SetRangeUser(0,20000);

    TF1 *GaussFit = new TF1("GaussFit","gausn",digisPerStripPerCh[c]->GetMean()-3*digisPerStripPerCh[c]->GetStdDev(),digisPerStripPerCh[c]->GetMean()+3*digisPerStripPerCh[c]->GetStdDev());
    digisPerStripPerCh[c]->Fit(GaussFit,"RQ");
    digisPerStripPerCh[c]->Fit(GaussFit,"RQ");
    digisPerStripPerCh[c]->Fit(GaussFit,"RQ");
    GaussFit->Draw("SAME");

    digisPerStripPerCh[c]->GetXaxis()->SetRangeUser(0,max_digi_occupancy[c]);

    if ( (GaussFit->GetParameter(1) + nsigmas*GaussFit->GetParameter(2)) > 0 )
    	HotStripLimitValue[c] = int(GaussFit->GetParameter(1) + nsigmas*GaussFit->GetParameter(2)); // Centroid of the gaussian + n sigmas

    namename = "Pos_" + to_string(chamberPos[i]) + "_Chamber_" + chamberNamePlot[i];
    infile->mkdir(namename.c_str());
    infile->cd(namename.c_str());

    namename = "Digi_PerStrip_PerCh_" + chamberNamePlot[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
    digisPerStripPerCh[c]->Write(namename.c_str());
    namename = "outPlots_Chamber_Pos_" + to_string(chamberPos[i]) + "_" + chamberNamePlot[i] + "/Digi_PerStrip_PerCh_" + to_string(chamberPos[i]) + ".png";
    Canvas->SaveAs(namename.c_str());

    infile->cd();
    delete GaussFit;
    Canvas->Clear();
  }

  // Dead / Hot strips in 2D digi-like and per VFAT plots

  TH2D *deadStrips2D[30];
  TH2D *hotStrips2D[30];
  TH2D *deadStripsVFAT[30];
  TH2D *hotStripsVFAT[30];
  TH2D *allMasksVFAT[30];

  for (int ch=0; ch<30; ch++)
  {
    sprintf(name,"DeadStrips2D_ch_%u",ch);
    deadStrips2D[ch] = new TH2D(name,"",384,0,384,8,-0.5,7.5);
    sprintf(name,"HotStrips2D_ch_%u",ch);
    hotStrips2D[ch] = new TH2D(name,"",384,0,384,8,-0.5,7.5);
    sprintf(name,"DeadStripsVFAT_ch_%u",ch);
    deadStripsVFAT[ch] = new TH2D(name,"",3,0,3,8,-0.5,7.5);
    sprintf(name,"HotStripsVFAT_ch_%u",ch);
    hotStripsVFAT[ch] = new TH2D(name,"",3,0,3,8,-0.5,7.5);
    sprintf(name,"AllMasksVFAT_ch_%u",ch);
    allMasksVFAT[ch] = new TH2D(name,"",3,0,3,8,-0.5,7.5);
  }

  // Dead / Hot strips results in csv files

  ofstream deadfile;
  string outFileName = "DeadStrips_run" + to_string(run) + "_ToDB.csv";
  deadfile.open(outFileName);

  ofstream hotfile;
  outFileName = "HotStrips_run" + to_string(run) + "_ToDB.csv";
  hotfile.open(outFileName);

  string entry = "";

  entry = "RUN_NUMBER," + to_string(run) + "\n";
  deadfile << entry;
  hotfile << entry;
  entry = "CH_SERIAL_NUMBER,GEM_NUM,POSITION,VFAT,CHANNEL,STRIP\n";
  deadfile << entry;
  hotfile << entry;

  for (unsigned int i=0; i<chamberPos.size(); i++)
  {
    int c = chamberPos[i];

    // Compose GEM_NUM

    string GEM_NUM = "11";
    if (c%2 == 0) GEM_NUM += "2"; // Even = Bottom
    if (c%2 == 1) GEM_NUM += "1"; // Odd = Top
    if (c < 10 && c%2 == 0) GEM_NUM += "0" + to_string(c+1);
    if (c < 10 && c%2 == 1) GEM_NUM += "0" + to_string(c);
    if (c >= 10 && c%2 == 0) GEM_NUM += to_string(c+1);
    if (c >= 10 && c%2 == 1) GEM_NUM += to_string(c);

    // Compose POSITION

    int ROW = int((c%10)/2)+1;
    int COLUMN = int(c/10)+1;
    string POSITION = "";
    if (c%2 == 0) POSITION = to_string(ROW) + "/" + to_string(COLUMN) + "/B"; // Even = Bottom
    if (c%2 == 1) POSITION = to_string(ROW) + "/" + to_string(COLUMN) + "/T"; // Odd = Top

    for (int eta=0; eta<8; eta++)
		{
		  for (int strip=0; strip<384; strip++)
		  {
		  	// Calculate VFAT

		  	int phi = int(strip/128);
		  	int VFAT = 8*phi+(7-eta);

		  	// Check if strip is dead or hot

		    if (digi2D[c]->GetBinContent(strip+1,eta+1) == DeadStripLimitValue[c])
		    {
		    	entry = chamberName[i] + "," + GEM_NUM + "," + POSITION + "," + to_string(VFAT) + "," + to_string(-1) + "," + to_string(strip) + "\n";
		    	deadfile << entry;
          deadStrips2D[c]->Fill(strip,eta);
          deadStripsVFAT[c]->Fill(phi,eta);
          allMasksVFAT[c]->Fill(phi,eta);
		    	continue;
		    }
		    if (digi2D[c]->GetBinContent(strip+1,eta+1) > HotStripLimitValue[c])
		    {
		    	entry = chamberName[i] + "," + GEM_NUM + "," + POSITION + "," + to_string(VFAT) + "," + to_string(-1) + "," + to_string(strip) + "\n";
		    	hotfile << entry;
          hotStrips2D[c]->Fill(strip,eta);
          hotStripsVFAT[c]->Fill(phi,eta);
          allMasksVFAT[c]->Fill(phi,eta);
		    }
		  }
		}
    deadStripsVFAT[c]->Scale(100.0/128.0);
    hotStripsVFAT[c]->Scale(100.0/128.0);
    allMasksVFAT[c]->Scale(100.0/128.0);
  }

  deadfile.close();
  hotfile.close();

  // Plotting 2D digi occupancy

  for (unsigned int i=0; i<chamberPos.size(); i++)
  {
    int c = chamberPos[i];

    namename = "Pos_" + to_string(chamberPos[i]) + "_Chamber_" + chamberNamePlot[i];
    infile->cd(namename.c_str());

    namename = "Digi_" + chamberName[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
    digi2D[c]->SetTitle(namename.c_str());
    digi2D[c]->GetXaxis()->SetTitle("Strip Number");
    digi2D[c]->GetYaxis()->SetTitle("ieta");
    digi2D[c]->Draw("colz");
    namename = "Digi_" + chamberNamePlot[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
    digi2D[c]->Write(namename.c_str());
    namename = "outPlots_Chamber_Pos_" + to_string(chamberPos[i]) + "_" + chamberNamePlot[i] + "/Digi_Ch_Pos_" + to_string(chamberPos[i]) + ".png";
    Canvas->SaveAs(namename.c_str());
    Canvas->Clear();

    namename = "DeadStrips_2D_" + chamberName[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
    deadStrips2D[c]->SetTitle(namename.c_str());
    deadStrips2D[c]->GetXaxis()->SetTitle("Strip Number");
    deadStrips2D[c]->GetYaxis()->SetTitle("ieta");
    deadStrips2D[c]->Draw("colz");
    namename = "DeadStrips_2D_" + chamberNamePlot[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
    deadStrips2D[c]->Write(namename.c_str());
    namename = "outPlots_Chamber_Pos_" + to_string(chamberPos[i]) + "_" + chamberNamePlot[i] + "/DeadStrips_2D_Ch_Pos_" + to_string(chamberPos[i]) + ".png";
    Canvas->SaveAs(namename.c_str());
    Canvas->Clear();

    namename = "HotStrips_2D_" + chamberName[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
    hotStrips2D[c]->SetTitle(namename.c_str());
    hotStrips2D[c]->GetXaxis()->SetTitle("Strip Number");
    hotStrips2D[c]->GetYaxis()->SetTitle("ieta");
    hotStrips2D[c]->Draw("colz");
    namename = "HotStrips_2D_" + chamberNamePlot[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
    hotStrips2D[c]->Write(namename.c_str());
    namename = "outPlots_Chamber_Pos_" + to_string(chamberPos[i]) + "_" + chamberNamePlot[i] + "/HotStrips_2D_Ch_Pos_" + to_string(chamberPos[i]) + ".png";
    Canvas->SaveAs(namename.c_str());
    Canvas->Clear();

    namename = "DeadStripsPerVFAT_" + chamberName[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
    deadStripsVFAT[c]->SetTitle(namename.c_str());
    deadStripsVFAT[c]->GetXaxis()->SetTitle("Strip Number");
    deadStripsVFAT[c]->GetYaxis()->SetTitle("ieta");
    deadStripsVFAT[c]->Draw("colz,text0");
    namename = "DeadStripsPerVFAT_" + chamberNamePlot[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
    deadStripsVFAT[c]->Write(namename.c_str());
    namename = "outPlots_Chamber_Pos_" + to_string(chamberPos[i]) + "_" + chamberNamePlot[i] + "/DeadStripsPerVFAT_Ch_Pos_" + to_string(chamberPos[i]) + ".png";
    Canvas->SaveAs(namename.c_str());
    Canvas->Clear();

    namename = "HotStripsPerVFAT_" + chamberName[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
    hotStripsVFAT[c]->SetTitle(namename.c_str());
    hotStripsVFAT[c]->GetXaxis()->SetTitle("Strip Number");
    hotStripsVFAT[c]->GetYaxis()->SetTitle("ieta");
    hotStripsVFAT[c]->Draw("colz,text0");
    namename = "HotStripsPerVFAT_" + chamberNamePlot[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
    hotStripsVFAT[c]->Write(namename.c_str());
    namename = "outPlots_Chamber_Pos_" + to_string(chamberPos[i]) + "_" + chamberNamePlot[i] + "/HotStripsPerVFAT_Ch_Pos_" + to_string(chamberPos[i]) + ".png";
    Canvas->SaveAs(namename.c_str());
    Canvas->Clear();

    namename = "AllMasksVFAT_" + chamberName[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
    allMasksVFAT[c]->SetTitle(namename.c_str());
    allMasksVFAT[c]->GetXaxis()->SetTitle("Strip Number");
    allMasksVFAT[c]->GetYaxis()->SetTitle("ieta");
    allMasksVFAT[c]->Draw("colz,text0");
    namename = "AllMasksVFAT_" + chamberNamePlot[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
    allMasksVFAT[c]->Write(namename.c_str());
    namename = "outPlots_Chamber_Pos_" + to_string(chamberPos[i]) + "_" + chamberNamePlot[i] + "/AllMasksVFAT_Ch_Pos_" + to_string(chamberPos[i]) + ".png";
    Canvas->SaveAs(namename.c_str());
    Canvas->Clear();

    infile->cd();
  }

  standConfigFile.close();
  infile->Close();
}
