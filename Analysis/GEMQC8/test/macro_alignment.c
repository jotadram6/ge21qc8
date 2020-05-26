#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include "TGraphErrors.h"
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

void macro_alignment(int run, string dataPath, int step)
{
  // Getting the root file

	string filename = "alignment_out_run_";
	for (unsigned int i=0; i<(6-to_string(run).size()); i++)
	{
		filename += "0";
	}
	filename += to_string(run) + ".root";
	const char *file = filename.c_str();

	TFile *infile = new TFile(file,"UPDATE");
	if (infile->IsOpen()) cout << "File opened successfully" << endl;

  //USeful variable declaration
  int chType[15] = {9,9,9,9,9,9,9,9,9,9,9,9,9,9,9};
  double Ypos[2][8] = {
    {72.747497,53.292499,35.486499,19.329500,4.4980001,-9.008000,-21.43950,-32.79650}, // zero is for long chambers
    {74.7225,58.0175,42.6115,28.5045,15.4475,3.4405,-7.6910,-17.9470} // one is for short chambers
  };

  // Getting the information on chamber type
  string infilename = dataPath + "/StandConfigurationTables/StandGeometryConfiguration_run" + to_string(run) + ".csv";
  ifstream standConfigFile(infilename);

  string line, split, comma = ",";
  size_t pos = 0;
  int ChPos = 0;

  if(standConfigFile.is_open())
  {
    while (getline(standConfigFile, line))
    {
      pos = line.find(comma);
			split = line.substr(0, pos);
			if (split == "CH_SERIAL_NUMBER") continue;
      line.erase(0, pos + comma.length());

      pos = line.find(comma);
			split = line.substr(0, pos);
      ChPos = int(stoi(string(1,split[3])+string(1,split[4]))/2.0);
      line.erase(0, pos + comma.length());

      pos = line.find(comma);
			line.erase(0, pos + comma.length());

      pos = line.find(comma);
			split = line.substr(0, pos);
      if (split == "L") chType[ChPos] = 0;
      if (split == "S") chType[ChPos] = 1;
      line.erase(0, pos + comma.length());
    }
    standConfigFile.close();
  }
  else cout << "Error opening file: " << infilename << endl;

  // Getting the number of events
  TH1D *hevt = (TH1D*)infile->Get("AlignmentQC8/goodVStriggeredEvts");
  int evt_tot_number = hevt->GetBinContent(1);
  cout << "Number of events = " << evt_tot_number << endl;

  // Histogram declaration
  char *histname = new char[20];
  char *histoname = new char[50];
  TH1D *resXperSC[15][8];
  for (int i_SC=0; i_SC<15; i_SC++)
  {
    for (int i_eta=0; i_eta<8; i_eta++)
    {
      sprintf(histname,"resX_SC_%u_%u_eta_%u",(i_SC%5)+1,(i_SC/5)+1,i_eta+1);
      resXperSC[i_SC][i_eta] = new TH1D(histname,"",300,-3,3);
    }
  }

  double resEtaY[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double resEtaYError[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  char *cnvname = new char[20];
  TCanvas *cnvResX[15];
  TCanvas *cnvResCorrPlot[15];

  double dx[15], dx_prev[15];
  double rz[15], rz_prev[15];

  for (int i_SC=0; i_SC<15; i_SC++)
  {
    dx[i_SC] = dx_prev[i_SC] = 0.0;
    rz[i_SC] = rz_prev[i_SC] = 0.0;
  }

  // Getting the dx and rz at the previous step
  string prevoutfilename = dataPath+"StandAligmentTables/StandAlignmentValues_run" +to_string(run) + "_ToDB.csv";
  string word;
  if(step>0){
    ifstream prevFile(prevoutfilename, ios::in);
    if(prevFile.is_open())
    {
      cout << "Reading the alignment factor at the previous step from " << prevoutfilename << endl;
      int i_SC = 0;
      while(getline(prevFile, line))
      {
        pos = line.find(",");
        split = line.substr(0, pos);
        if(split == "RunNumber" || pos == string::npos) continue;
        else if(split == "Position" || pos == string::npos) continue;
        else
        {
          stringstream linestream(line);
          int i = 0;
          while(getline(linestream,word, ','))
					{
            if(i==1) dx_prev[i_SC] = stod(word);
            if(i==6) rz_prev[i_SC] = stod(word);
            i += 1;
          }
          i_SC += 1;
        }
      }
			prevFile.close();
    }
    else cout << "Error opening file: " << prevoutfilename << endl;
  }

  for (int i_SC=0; i_SC<15; i_SC++)
  {
    cout << "prev dx " <<  dx_prev[i_SC];
    cout << " prev rz " <<  rz_prev[i_SC] << endl;
  }

  // Residuals per SC
  cout << "Starting the analysis of run " << run << endl;
  for(int i_SC=0; i_SC<15; i_SC++)
  {
		if (chType[i_SC] == 9) continue;

    sprintf(cnvname,"cnv_SC_%u_%u",(i_SC%5)+1,(i_SC/5)+1);
    cnvResX[i_SC] = new TCanvas(cnvname,cnvname,0,0,1000,600);
    cnvResX[i_SC]->Divide(4,2);
    for(int i_eta=0; i_eta<8; i_eta++)
    {
      cnvResX[i_SC]->cd(i_eta+1);
      sprintf(histoname,"AlignmentQC8/h_resX_eta_%d_%d", i_SC+1, i_eta+1);
      resXperSC[i_SC][i_eta]=(TH1D*)infile->Get(histoname);
      resXperSC[i_SC][i_eta]->Draw();
      TF1 *GaussFit = new TF1("GaussFit","gaus",-3,3);
      resXperSC[i_SC][i_eta]->Fit(GaussFit,"Q");
      GaussFit->Draw("SAME");
      resEtaY[i_eta] = GaussFit->GetParameter(1);
      resEtaYError[i_eta] = GaussFit->GetParError(1);
      delete GaussFit;
    }

    cnvResCorrPlot[i_SC] = new TCanvas(cnvname,cnvname,0,0,1000,600);
    sprintf(histname,"resCorrPlot_SC_%u_%u",(i_SC%5)+1,(i_SC/5)+1);
    TGraphErrors *resCorrPlotSC = new TGraphErrors(8,Ypos[chType[i_SC]],resEtaY,0,resEtaYError);
    resCorrPlotSC->SetTitle(histname);
    resCorrPlotSC->SetMarkerSize(1.5);
    resCorrPlotSC->SetMarkerStyle(21);
    resCorrPlotSC->Draw("ap");
    TF1 *LinFit = new TF1("LinFit","pol1",Ypos[chType[i_SC]][7]-2,Ypos[chType[i_SC]][0]+2);
    resCorrPlotSC->Fit(LinFit,"Q");
    resCorrPlotSC->Write(histname);
    dx[i_SC] = -LinFit->GetParameter(0);
    rz[i_SC] = atan(LinFit->GetParameter(1))*180.0/3.1415926536; // rotation in deg!
    delete LinFit;
  }

  cout << "shiftX = [";
  for (int i_SC=0; i_SC<15; i_SC++)
  {
    cout << dx[i_SC] << ",";
    if (i_SC == 14) cout << dx[i_SC];
    if (i_SC == 4 || i_SC == 9) cout << "\\" << endl;
    if (i_SC == 14) cout << "]\n" << endl;
  }

  cout << "rotationZ = [";
  for (int i_SC=0; i_SC<15; i_SC++)
  {
    cout << rz[i_SC] << ",";
    if (i_SC == 14) cout << rz[i_SC];
    if (i_SC == 4 || i_SC == 9) cout << "\\" << endl;
    if (i_SC == 14) cout << "]\n" << endl;
  }
  // Writing the output in csv format
  string partialoutfilename = "StandAlignmentValues_run" + to_string(run) + "_step" + to_string(step) + ".csv";
  ofstream partialFile(partialoutfilename, std::ios_base::out | std::ios_base::trunc);
  if(partialFile.is_open())
  {
    cout << "Writing the partial alignment factors in " << partialoutfilename << endl;
    partialFile << "RunNumber,"<< run <<",,,,,\n";
    partialFile << "Position,dx(cm),dy(cm),dz(cm),rx(deg),ry(deg),rz(deg)\n";
    for (int i_SC=0; i_SC<15; i_SC++)
    {
      partialFile << (i_SC%5)+1 << "/" << (i_SC/5)+1 << "," << dx[i_SC] << ",0,0,0,0," << rz[i_SC] << "\n";
    }
    partialFile.close();
  }
  else
  {
    cerr << partialoutfilename << " is not correctly opened" << endl;
    exit(EXIT_FAILURE);
  }

  // Writing the output in csv format
  string totaloutfilename = "StandAlignmentValues_run" + to_string(run) + "_ToDB.csv";
  ofstream totalFile(totaloutfilename, std::ios_base::out | std::ios_base::trunc);
  if(totalFile.is_open())
  {
    cout << "Writing the total alignment factors in " << totaloutfilename << endl;
    totalFile << "RunNumber,"<< run <<",,,,,\n";
    totalFile << "Position,dx(cm),dy(cm),dz(cm),rx(deg),ry(deg),rz(deg)\n";
    for (int i_SC=0; i_SC<15; i_SC++)
    {
      totalFile << (i_SC%5)+1 << "/" << (i_SC/5)+1 << "," << dx[i_SC]+dx_prev[i_SC] << ",0,0,0,0," << rz[i_SC]+rz_prev[i_SC] << "\n";
    }
    totalFile.close();
  }
  else
  {
    cerr << totaloutfilename << " is not correctly opened" << endl;
    exit(EXIT_FAILURE);
  }
}
