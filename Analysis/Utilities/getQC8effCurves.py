import cx_Oracle
import os, sys, io
import datetime
import time
import csv
import string
import pandas as pd
import numpy as np
import ROOT
import array

def getSCname(det_name):

    if det_name.startswith("GE11-"):
        det_name = det_name.replace("GE11-","GE1/1-")
        
    db_cond = os.environ["GEM_PRODUCTION_DB_COND"]
    db_name = os.environ["GEM_PRODUCTION_DB_NAME"]
    
    db = cx_Oracle.connect(db_cond+db_name) # production DB
    cur = db.cursor()
    query = "select * from CMS_GEM_MUON_VIEW.GEM_SPRCHMBR_PARTS_VIEW_RH"
    cur.execute(query)
    
    for result in cur:
        if (result[3] == det_name and result[11] == "GEM1"):
            SC_name = str(result[0]).replace('GE1/1-','GE11-')
            return SC_name

def getLayName(SC_name,lay):

    if SC_name.startswith("GE11-"):
        SC_name = SC_name.replace("GE11-","GE1/1-")
        
    db_cond = os.environ["GEM_PRODUCTION_DB_COND"]
    db_name = os.environ["GEM_PRODUCTION_DB_NAME"]
    
    db = cx_Oracle.connect(db_cond+db_name) # production DB
    cur = db.cursor()
    query = "select * from CMS_GEM_MUON_VIEW.GEM_SPRCHMBR_PARTS_VIEW_RH"
    cur.execute(query)
    
    for result in cur:
        if (result[0] == SC_name and result[5] == str(lay) and result[11] == "GEM1"):
            lay_name = str(result[3]).replace('GE1/1-','GE11-')
            return lay_name

if __name__ == '__main__':

    ROOT.gROOT.SetBatch(True)
    
    # Define the parser
    import argparse
    from argparse import RawTextHelpFormatter
    parser = argparse.ArgumentParser(description="This script creates \"Efficiency vs gain\" or \"Efficiency vs EDcurrent\" or \"Efficiency vs Threshold\" for the chambers in the specified runs.\n\nHow to use examples:\n\n python getQC8effCurves.py /fullPath/runList.txt gain\n\n python getQC8effCurves.py /fullPath/runList.txt EDcurrent\n\n python getQC8effCurves.py /fullPath/runList.txt Threshold", formatter_class=RawTextHelpFormatter)
    # Positional arguments
    parser.add_argument("runListFile", type=str, help="Please provide the full path to the desired file with list of runs to be used")
    parser.add_argument("xAxisVar", type=str, default="gain", choices=["gain","EDcurrent","Threshold"], help="Please specify the variable to be used for the x axis. Choices: \"gain\" or \"EDcurrent\" or \"Threshold\"")
    args = parser.parse_args()
    
    try:
        infile = open(args.runListFile, 'rb')
    except IOError:
        sys.exit("Could not read file: " + args.runListFile)
    
    maskCol = [None]*3
    runList,maskCol[0],maskCol[1],maskCol[2] = np.loadtxt(infile, delimiter="\t", skiprows=1, unpack=True)
    
    infile.close()
    
    SCList = []
    vfatList = map(str, range(0,24))
 
    for run in runList:
        RunLogsFolder =  "/data/bigdisk/GEM-Data-Taking/GE11_QC8/QC8_RunLogs/Run_{:06d}/".format(int(run))
        for logs in os.listdir(RunLogsFolder):
            chID = logs.split('.')[0].split('_')[2]
            if getSCname(chID) not in SCList:
                SCList.append(getSCname(chID))
                
    for SC in SCList:
        ROOT.gStyle.SetPalette(ROOT.kRainBow)
        c = ROOT.TCanvas("c","",100,100,1200,1200)
        c.SetGrid()
        #legend = ROOT.TLegend(0.15,0.15,0.85,0.4)
        legend = ROOT.TLegend(0.15,0.15,0.85,0.33)
        legend.SetHeader("#splitline{#splitline{}{Gas mixture: Ar/CO_{2} (70/30%)}}{#splitline{}{Super Chamber ID (Layer) == Chambers ID}}")
        x1, y1, xerr1, yerr1 = array.array( 'd' ), array.array( 'd' ), array.array( 'd' ), array.array( 'd' )
        x2, y2, xerr2, yerr2 = array.array( 'd' ), array.array( 'd' ), array.array( 'd' ), array.array( 'd' )
        
        for layer in [1,2]:
            eff = 0.0
            eff_err = 0.0
            chName = getLayName(SC,layer)
            for run in runList:
                chLogFile = "/data/bigdisk/GEM-Data-Taking/GE11_QC8/QC8_RunLogs/Run_{:06d}/".format(int(run)) + "QC8_Run{:06d}_".format(int(run)) + str(chName) +".csv"
                try:
                    chLog = open(chLogFile, 'rb')
                except IOError:
                    print("Could not read file: " + chLogFile)
                    continue
                    #sys.exit("Could not read file: " + chLogFile)
                
                quantity,value = np.genfromtxt(chLog, dtype='str', delimiter="\t", unpack=True)
            
                row = int((str(value[np.where(quantity == "Position")]).split('/')[0])[-1])
                column = int(str(value[np.where(quantity == "Position")]).split('/')[1])
                
                if (str(bin(int(maskCol[column-1][np.where(runList == run)]))[2:].zfill(12))[-(3+2*(row-1)+(layer-1))] == "1"): # if mask value is 1 for that chamber
                    
                    if (args.xAxisVar == "gain" and layer == 1):
                        x1.append(float(value[np.where(quantity == "EffAvgGain")]))
                        xerr1.append(0.0)
                    elif (args.xAxisVar == "EDcurrent" and layer == 1):
                        x1.append(float(value[np.where(quantity == "I0(uA)")]))
                        xerr1.append(0.0)
                    elif (args.xAxisVar == "Threshold" and layer == 1):
                        boolVfat = np.in1d(quantity, vfatList) #returns True if quantity == 0, 1, .. 23
                        thrList = value[boolVfat].astype(np.float) #list of threshold values for the 24 VFATs
                        x1.append(float(np.mean(thrList)))
                        xerr1.append(0.0)
                        
                    if (args.xAxisVar == "gain" and layer == 2):
                        x2.append(float(value[np.where(quantity == "EffAvgGain")]))
                        xerr2.append(0.0)
                    elif (args.xAxisVar == "EDcurrent" and layer == 2):
                        x2.append(float(value[np.where(quantity == "I0(uA)")]))
                        xerr2.append(0.0)
                    elif (args.xAxisVar == "Threshold" and layer == 2):
                        boolVfat = np.in1d(quantity, vfatList) #returns True if quantity == 0, 1, .. 23
                        thrList = value[boolVfat].astype(np.float) #list of threshold values for the 24 VFATs
                        x2.append(float(np.mean(thrList)))
                        xerr2.append(0.0)
                        
                    resultsFile = "/home/gemuser/AnalysisQC8/Results_QC8_validation_run_" + str(int(run)) + "_yesMasks/Average_Efficiency_Per_Chamber.csv"
                    
                    try:
                        results = open(resultsFile, 'rb')
                    except IOError:
                        print("Could not read file: " + resultsFile)
                        continue
                        #sys.exit("Could not read file: " + resultsFile)
                    
                    PositionCMSSW,Position,Chamber,Efficiency,ErrorEfficiency = np.genfromtxt(results, dtype='str', delimiter=",", skip_header=1, unpack=True)
                    
                    eff = float(Efficiency[np.where(Chamber == chName.replace('GE11-','GE1/1-'))])
                    eff_err = float(ErrorEfficiency[np.where(Chamber == chName.replace('GE11-','GE1/1-'))])
                    
                    results.close()
                    
                    """
                    
                    resultsFile = "/home/gemuser/AnalysisQC8/Results_QC8_validation_run_" + str(int(run)) + "_noMasks/Average_Efficiency_Per_Chamber.csv"
                    
                    try:
                        results = open(resultsFile, 'rb')
                    except IOError:
                        sys.exit("Could not read file: " + resultsFile)
                    
                    PositionCMSSW,Position,Chamber,Efficiency,ErrorEfficiency = np.genfromtxt(results, dtype='str', delimiter=",", skip_header=1, unpack=True)
                    
                    if (float(Efficiency[np.where(Chamber == chName.replace('GE11-','GE1/1-'))]) > eff):
                        eff = float(Efficiency[np.where(Chamber == chName.replace('GE11-','GE1/1-'))])
                        eff_err = float(ErrorEfficiency[np.where(Chamber == chName.replace('GE11-','GE1/1-'))])
                    
                    results.close()
                    
                    """
                    
                    if (layer == 1):
                        y1.append(eff)
                        yerr1.append(eff_err)
                    elif (layer == 2):
                        y2.append(eff)
                        yerr2.append(eff_err)
                    
                chLog.close()
                
                
        if (len(x1)==0 or len(x2)==0):
            continue
        
        print SC, runList, getLayName(SC,1), x1, y1, yerr1
        print SC, runList, getLayName(SC,2), x2, y2, yerr2
                
        graph1 = ROOT.TGraphErrors(int(len(x1)), x1, y1, xerr1, yerr1)
        graph2 = ROOT.TGraphErrors(int(len(x2)), x2, y2, xerr2, yerr2)
        
        graph1.SetTitle("Efficiency vs gain")
        if (args.xAxisVar == "gain"):
        	graph1.GetXaxis().SetRangeUser(6000, 23000)
        	graph1.GetXaxis().SetTitle("Avg. Effective Gain")
        if (args.xAxisVar == "EDcurrent"):
        	#graph1.GetXaxis().SetRangeUser(650, 720)
        	graph1.GetXaxis().SetRangeUser(580, 720)
        	graph1.GetXaxis().SetTitle("Equivalent Divider Current (uA)")
        if (args.xAxisVar == "Threshold"):
        	graph1.GetXaxis().SetRangeUser(1, 21)
        	graph1.GetXaxis().SetTitle("Average Threshold (fC)")
        graph1.GetXaxis().SetTitleOffset(1.1)
        graph1.GetXaxis().SetLabelSize(0.03)
        #graph1.GetYaxis().SetRangeUser(0.5, 1.0)
        graph1.GetYaxis().SetRangeUser(-0.3, 1.0)
        graph1.GetYaxis().SetDecimals()
        graph1.GetYaxis().SetLabelSize(0.03)
        graph1.GetYaxis().SetTitleOffset(1.3)
        graph1.GetYaxis().SetTitle("MIP Efficiency")
        graph1.SetLineWidth(1)
        graph1.Sort()
        graph1.Draw("ALP PLC PFC")
        c.Update()
        namename = str(SC) + " (Layer1) == " + str(getLayName(SC,1))
        legend.AddEntry(graph1,namename,"l")
        
        graph2.SetLineWidth(1)
        graph2.Sort()
        graph2.Draw("LP PLC PFC SAME")
        c.Update()
        namename = str(SC) + " (Layer2) == " + str(getLayName(SC,2))
        legend.AddEntry(graph2,namename,"l")
                
        legend.Draw("SAME")
        
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextAngle(0)
        latex.SetTextColor(ROOT.kBlack)
        latex.SetTextFont(42)
        latex.DrawLatex(0.104,0.91,"#scale[0.5]{#bf{CMS} #it{Preliminary}}")
        
        c.SaveAs(str(SC) + "_Eff" + args.xAxisVar + "Curve.png")
            
            
            
            
            
            
            
            
                    
        
        
        




