#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include "TGraphErrors.h"
#include <TBranch.h>
#include <TCanvas.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

void macro_tilt_twist(int run, string runPath, string AlignTablePath){
  // Getting the file
  char *filename = new char[70];
  char run_number[6] = "";
  for(int i=6; i>(int)(log10(run)+1); i--){
    strcat(run_number,"0");
  }
  sprintf(filename, "alignment_out_run_%s%d.root", run_number, run);

  if(gSystem->AccessPathName(filename))
    {
      cout << "ERROR: ROOT file " << filename << " does not exist!" << endl;
      exit(EXIT_FAILURE);
    }
  else
    {
      cout << "---ROOT file " << filename << " successfully opened...Loading" << endl;
    }
  TFile *infile = TFile::Open(filename, "UPDATE");

  //USeful variable declaration
  int i_Eta = 8, i_Col = 3, i_maxSC = 15;
  double dx_tilt[i_maxSC], rz_twist[i_maxSC];
  const double Zpos[10] = {22.66,26.13,45.66,49.13,68.66,72.13,91.66,95.13,114.66,118.13}; // z position of recHits in each column (refers to the single chamber)
  int chidx[i_maxSC];
  double Ypos[3][8] = {
    {72.747497,53.292499,35.486499,19.329500,4.4980001,-9.008000,-21.43950,-32.79650}, // zero is for long chambers
    {74.7225,58.0175,42.6115,28.5045,15.4475,3.4405,-7.6910,-17.9470}, // one is for short chambers
    {-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.}
  };

  // Getting the information on chamber type
  ifstream standConfigFile;
  string configName = runPath+"/configureRun_cfi.py";
  standConfigFile.open(configName);
  string line, split, equal = " = ";
  vector<string> chType;
  char comma = ',';
  size_t pos;
  while (getline(standConfigFile, line))
    {
      pos = line.find(equal);
      split = line.substr(0, pos);
      if (split != "StandConfiguration" || pos == string::npos) continue;
      else if (split == "StandConfiguration")
	{
	  while (getline(standConfigFile, line))
	    {
	      stringstream linestream(line);
	      for (string each; getline(linestream, each, comma); chType.push_back(each));
	      if(chType.at(chType.size()-1)== "\\") chType.pop_back();
	      if(chType.size()==15) continue;
	    }
	}
    }

  for(int i_SC=0; i_SC<i_maxSC; i_SC++)
    {
      if(chType.at(i_SC).substr(0,1)=='0') chidx[i_SC]=2;
      else if(chType.at(i_SC).substr(0,1)=='S') chidx[i_SC]=1;
      else  chidx[i_SC]=0;
  }
  // Getting the TTree
  TString direc = "AlignmentQC8/";
  TTree *datatree = (TTree*)infile->Get(direc+"tree");
  double dx_prev[i_maxSC], rz_prev[i_maxSC];
  for (int i_SC=0; i_SC<i_maxSC; i_SC++)
    {
      dx_prev[i_SC]=90000.;
      rz_prev[i_SC]=90000.;
    }
  // Getting the dx and rz at the previous step
  string prevoutfilename = AlignTablePath+"StandAlignmentValues_run" +to_string(run) + ".csv";
  string word;
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
	      while(getline(linestream,word, ',')){
		if(i==1) dx_prev[i_SC] = stod(word);
		if(i==6) rz_prev[i_SC] = stod(word);
		i += 1;
	      }
	      i_SC += 1;
	    }
	}
    }
  else
    {
      cerr << prevoutfilename << " is not correctly opened" << endl;
      exit(EXIT_FAILURE);
    }

  for (int i_SC=0; i_SC<i_maxSC; i_SC++)
    {
      cout << "prev dx " <<  dx_prev[i_SC];
      cout << " prev rz " <<  rz_prev[i_SC] << endl;
    }

  // Histogram declaration
  char *histname = new char[20];
  char *histoname = new char[20];
  TH1D *hColEtaPxPZ[i_Col][i_Eta];
  for (int i_col=0; i_col<i_Col; i_col++)
    {
      for (int i_eta=0; i_eta<i_Eta; i_eta++)
        {
          sprintf(histname,"hColEtaPxPz_%d_%d",i_col+1,i_eta+1);
	  hColEtaPxPZ[i_col][i_eta] = new TH1D(histname,"",300,-3,3);
        }
    }

  double AngDistrEtaY[i_Eta];
  double AngDistrEtaYError[i_Eta];
  double intercept[i_Col], slope[i_Col];
  char *cnvname = new char[20];
  TCanvas *cnvAngDistrCorrPlot[i_Col];
  for (int i_col=0; i_col<i_Col; i_col++)
    {
      for (int i_eta=0; i_eta<i_Eta; i_eta++)
        {
	  sprintf(cnvname,"cnv_ColEtaPxPz_%d_%d",(i_col)+1,(i_eta)+1);
          sprintf(histoname,"h_PxPz_col_eta_%d_%d",i_col+1,i_eta+1);
	  hColEtaPxPZ[i_col][i_eta]=(TH1D*)infile->Get(direc+histoname);
	  AngDistrEtaY[i_eta]=hColEtaPxPZ[i_col][i_eta]->GetMean();
	  AngDistrEtaYError[i_eta]=hColEtaPxPZ[i_col][i_eta]->GetMeanError();
	}
      cnvAngDistrCorrPlot[i_col] = new TCanvas(cnvname,cnvname,0,0,1000,600);
      sprintf(histname,"angdistrCorrPlot_SC_%u",i_col+1);
      TGraphErrors *angdistrCorrPlotSC = new TGraphErrors(i_Eta,Ypos[chidx[i_col*5]],AngDistrEtaY,0,AngDistrEtaYError);
      angdistrCorrPlotSC->SetTitle(histname);
      angdistrCorrPlotSC->SetMarkerSize(1.5);
      angdistrCorrPlotSC->SetMarkerStyle(21);
      angdistrCorrPlotSC->Draw("ap");
      TF1 *LinFit = new TF1("LinFit","pol1",Ypos[chidx[i_col*5]][i_Eta-1]-2,Ypos[chidx[i_col*5]][0]+2);
      angdistrCorrPlotSC->Fit(LinFit,"Q");
      angdistrCorrPlotSC->Write(histname);
      intercept[i_col] = LinFit->GetParameter(0);
      slope[i_col] = LinFit->GetParameter(1);
      delete LinFit;
    }
  for(int i_SC=0; i_SC<i_maxSC; i_SC++)
    {
      double z0 = (Ypos[chidx[i_SC]][5]-Ypos[chidx[i_SC]][4])/2;
      double z = Ypos[chidx[i_SC]][i_SC%5+1]-Ypos[chidx[i_SC]][i_SC%5];
      dx_tilt[i_SC] = intercept[i_SC/5]*(z-z0);
      rz_twist[i_SC] = asin(slope[i_SC/5]*(z-z0));
    }

  for (int i_SC=0; i_SC<i_maxSC; i_SC++)
    {
      dx_tilt[i_SC] += dx_prev[i_SC];
      rz_twist[i_SC] += rz_prev[i_SC];
    }

  // Writing the output in csv format
  string ofilename = AlignTablePath+"StandAlignmentValues_run" + to_string(run) + "_ToDB.csv";
  ofstream oFile(ofilename, std::ios_base::out | std::ios_base::trunc);
  if (oFile.is_open())
    {
      cout << "Writing the new alignment factors in " << ofilename << endl;
      oFile << "RunNumber,"<< run <<",,,,,\n";
      oFile << "Position,dx(cm),dy(cm),dz(cm),rx(deg),ry(deg),rz(deg)\n";
      for (int i_SC=0; i_SC<i_maxSC; i_SC++)
	{
	  oFile << (i_SC%5)+1 << "/" << (i_SC/5)+1 << "," << dx_tilt[i_SC] << ",0" << ",0," << rz_twist[i_SC] << ",0" << ",0\n";
	}
      oFile.close();
    }
}
