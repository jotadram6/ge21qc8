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

void macro_certify_events(int run, string configDir)
{
  // Getting the root file

  string filename = "certify_events_out_run_";
  for (unsigned int i=0; i<(6-to_string(run).size()); i++)
  {
    filename += "0";
  }
  filename += to_string(run) + ".root";
  const char *file = filename.c_str();

  TFile *infile = new TFile(file,"UPDATE");
  if (infile->IsOpen()) cout << "File opened successfully" << endl;

  // Getting the tree

  TTree *tree = (TTree*)infile->Get("CertifyEventsQC8/tree");
	int nTotEvents = tree->GetEntries();

	// Declaration of name variables

  char *name = new char[40];
  string namename = "";

  // Getting digis per ch vs event number histrogram

	TH2D *digisVsEvt = (TH2D*)infile->Get("CertifyEventsQC8/nDigisPerEvtPerCh");

	TH1D *NdigisPerChVsEvt[30];
	TH1D *NdigisPerChDistribution[30];

	for (unsigned int ch=0; ch<30; ch++)
	{
		sprintf(name,"NdigisVsEvent_ch_%u",ch);
		NdigisPerChVsEvt[ch] = new TH1D(name,"",24000,0,24000000);
		sprintf(name,"NdigisDistribution_ch_%u",ch);
		NdigisPerChDistribution[ch] = new TH1D(name,"",1000,1,10001);

		for (int evt=0; evt<24000; evt++)
		{
			NdigisPerChVsEvt[ch]->SetBinContent((evt+1),digisVsEvt->GetBinContent(evt+1,ch+1));
			NdigisPerChDistribution[ch]->Fill(digisVsEvt->GetBinContent(evt+1,ch+1));
		}
	}

	// Getting rechHits per ch vs event number histrogram

	TH2D *recHitsVsEvt = (TH2D*)infile->Get("CertifyEventsQC8/nRecHitsPerEvtPerCh");

	TH1D *NrecHitsPerChVsEvt[30];
	TH1D *NrecHitsPerChDistribution[30];

	for (unsigned int ch=0; ch<30; ch++)
	{
		sprintf(name,"NrecHitsVsEvent_ch_%u",ch);
		NrecHitsPerChVsEvt[ch] = new TH1D(name,"",24000,0,24000000);
		sprintf(name,"NrecHitsDistribution_ch_%u",ch);
		NrecHitsPerChDistribution[ch] = new TH1D(name,"",1000,1,10001);

		for (int evt=0; evt<24000; evt++)
		{
			NrecHitsPerChVsEvt[ch]->SetBinContent((evt+1),recHitsVsEvt->GetBinContent(evt+1,ch+1));
			NrecHitsPerChDistribution[ch]->Fill(recHitsVsEvt->GetBinContent(evt+1,ch+1));
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

  // Results for the 30 chambers

  vector<int> beginBadEvt[chamberPos.size()];
  vector<int> endBadEvt[chamberPos.size()];

  TCanvas *Canvas = new TCanvas("Canvas","Canvas",0,0,1000,800);

  ofstream outfile;

  for (unsigned int i=0; i<chamberPos.size(); i++)
  {
    int c = chamberPos[i];

    namename = "Pos_" + to_string(chamberPos[i]) + "_Chamber_" + chamberNamePlot[i];
    infile->mkdir(namename.c_str());
    infile->cd(namename.c_str());

    // Plotting number of digis per chamber vs evt

		namename = "NdigisVsEvt_" + chamberName[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
		NdigisPerChVsEvt[c]->SetTitle(namename.c_str());
		NdigisPerChVsEvt[c]->GetXaxis()->SetTitle("Evt Number");
		NdigisPerChVsEvt[c]->GetYaxis()->SetTitle("nDigis");
		NdigisPerChVsEvt[c]->GetXaxis()->SetRangeUser(0,nTotEvents);
		NdigisPerChVsEvt[c]->Draw();
    namename = "NdigisVsEvt_" + chamberNamePlot[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
		NdigisPerChVsEvt[c]->Write(namename.c_str());
		namename = "outPlots_Chamber_Pos_" + to_string(chamberPos[i]) + "_" + chamberNamePlot[i] + "/NdigisVsEvt_Ch_Pos_" + to_string(chamberPos[i]) + ".png";
		Canvas->SaveAs(namename.c_str());

		// Plotting number of digis per chamber vs evt

		namename = "NdigisDistribution_" + chamberName[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
		NdigisPerChDistribution[c]->SetTitle(namename.c_str());
		NdigisPerChDistribution[c]->GetXaxis()->SetTitle("nDigis");
		NdigisPerChDistribution[c]->GetYaxis()->SetTitle("Counts");
		NdigisPerChDistribution[c]->Draw();
    namename = "NdigisDistribution_" + chamberNamePlot[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
		NdigisPerChDistribution[c]->Write(namename.c_str());
		namename = "outPlots_Chamber_Pos_" + to_string(chamberPos[i]) + "_" + chamberNamePlot[i] + "/NdigisDistribution_Ch_Pos_" + to_string(chamberPos[i]) + ".png";
		Canvas->SaveAs(namename.c_str());

		// Plotting number of recHits per chamber vs evt

		namename = "NrecHitsVsEvt_" + chamberName[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
		NrecHitsPerChVsEvt[c]->SetTitle(namename.c_str());
		NrecHitsPerChVsEvt[c]->GetXaxis()->SetTitle("Evt Number");
		NrecHitsPerChVsEvt[c]->GetYaxis()->SetTitle("nRecHits");
		NrecHitsPerChVsEvt[c]->GetXaxis()->SetRangeUser(0,nTotEvents);
		NrecHitsPerChVsEvt[c]->Draw();
    namename = "NrecHitsVsEvt_" + chamberNamePlot[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
		NrecHitsPerChVsEvt[c]->Write(namename.c_str());
		namename = "outPlots_Chamber_Pos_" + to_string(chamberPos[i]) + "_" + chamberNamePlot[i] + "/NrecHitsVsEvt_Ch_Pos_" + to_string(chamberPos[i]) + ".png";
		Canvas->SaveAs(namename.c_str());

		// Plotting number of recHits per chamber vs evt

		namename = "NrecHitsDistribution_" + chamberName[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
		NrecHitsPerChDistribution[c]->SetTitle(namename.c_str());
		NrecHitsPerChDistribution[c]->GetXaxis()->SetTitle("nRecHits");
		NrecHitsPerChDistribution[c]->GetYaxis()->SetTitle("Counts");
		NrecHitsPerChDistribution[c]->Draw();
    namename = "NrecHitsDistribution_" + chamberNamePlot[i] + "_in_position_" + to_string(chamberPos[i]) + "_run_" + to_string(run);
		NrecHitsPerChDistribution[c]->Write(namename.c_str());
		namename = "outPlots_Chamber_Pos_" + to_string(chamberPos[i]) + "_" + chamberNamePlot[i] + "/NrecHitsDistribution_Ch_Pos_" + to_string(chamberPos[i]) + ".png";
		Canvas->SaveAs(namename.c_str());

    // Here looking for events with tripped chamber
    int limitRecHitHVtrip = 30;
    int binToEvt = 1000;
    int binBegin = 0, binEnd = 0;

    for (binBegin = 0; binBegin < 24000; binBegin++)
    {
      if (NrecHitsPerChVsEvt[c]->GetBinContent(binBegin+1) < limitRecHitHVtrip)
      {
        cout << "Bad event: binBegin!" << endl;
        for (binEnd = (binBegin+1); binEnd < 24000; binEnd++)
        {
          if (NrecHitsPerChVsEvt[c]->GetBinContent(binEnd+1) > limitRecHitHVtrip)
          {
            cout << "Bad event: binEnd!" << endl;
            break;
          }
        }
        if ((binEnd-binBegin) >= 2)
        {
          if (binBegin < 2) beginBadEvt[i].push_back(0 * binToEvt);
          if (binBegin >= 2) beginBadEvt[i].push_back((binBegin-2) * binToEvt);
          if (binEnd <= (24000-2)) endBadEvt[i].push_back((binEnd+2) * binToEvt);
          if (binEnd > (24000-2)) endBadEvt[i].push_back(24000 * binToEvt);
        }
        binBegin = binEnd;
      }
    }

    Canvas->Clear();
    infile->cd();
  }

	// csv file with certified events per chamber

	string outFileName = "CertifiedEvents_run" + to_string(run) + ".csv";
	outfile.open(outFileName);

	string entry = "";

	entry = "RunNumber," + to_string(run) + "\n";
	outfile << entry;

	for (unsigned int i=0; i<chamberPos.size(); i++)
  {
    int c = chamberPos[i];

		entry = to_string(c) + ",";

    for (unsigned int badIntervals = 0; badIntervals < beginBadEvt[i].size(); badIntervals++)
    {
      entry += to_string(beginBadEvt[i].at(badIntervals)) + "-" + to_string(endBadEvt[i].at(badIntervals));
      if (badIntervals != (beginBadEvt[i].size()-1)) entry += ",";
    }

		entry += "\n";

		outfile << entry;
	}

	outfile.close();

  standConfigFile.close();

  infile->Close();
}
