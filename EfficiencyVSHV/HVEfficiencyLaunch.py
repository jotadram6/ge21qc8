#This script takes in input n file with the names of QC8 chambers
#and calls the function for the monitoring of efficiency

import HVEfficiencyMonitor
import ROOT
import os, sys, io
from datetime import datetime

#ChambersList
#write chambers and run numbers in this format, sepataing chambers
#from runNumbers with a \t
#GE1/1-VII-L-CERN-0001	11,12
#GE1/1-VII-L-CERN-0002	11,12

fileChamber = sys.argv[1]

#fileChamber = "ChambersList.txt"

print ("Efficiency monitor: reading file "+ fileChamber + "\n")

#create the Root file with a name without points or spaces
today = str(datetime.today())
today = today.replace(" ", "_")
today = today.replace(":", "-")
today = today[:-7]

fileRoot = "QC8_Efficiency_monitor_"+ today +".root"
fr = ROOT.TFile( fileRoot,"RECREATE")

#open file ChambersList.txt
fc = open(fileChamber, "r")

#extract from each line of ChambersList.txt the name of the chamber and the run numbers asked for each chamber
for lineNameChamber in fc:
	#print "\n"
	
	lineNameChamber = lineNameChamber[:-1] #remove \n at the end of the line
	
	spaceIdx = lineNameChamber.find('\t')
	chamberName = lineNameChamber[:spaceIdx]
	chamberNameDir = chamberName.replace("/", "-")
	chamberNameDir = chamberNameDir.replace("-", "_")

	stringNumbers = lineNameChamber[(spaceIdx+1):]
	stringNumbers = stringNumbers.replace(","," ")
	runNumberList = stringNumbers.split()

	print "Chamber: "+chamberName
	print "RunNumbers: "+str(runNumberList)

	#transorm this list of strings in a list of numbers
	for x in runNumberList:
		idxElement = runNumberList.index(x)	
		x = int(x)
		runNumberList[idxElement] = x

	 
	print "\n-----------------------------------------------------------------------------------\n"
	#print "Chamber: "+chamberName+" | RunNumbers: "+str(runNumberList)

	#create directories                               	
        firstDir = fr.mkdir( chamberNameDir )
	firstDir.cd()

	#call the function for the Efficiency monitoring
	HVEfficiencyMonitor.HVEfficiencyMonitor( chamberName, runNumberList )
	
fc.close()
fr.Close()

print('\n-------------------------Output--------------------------------')
print(fileRoot+ " has been created.")
print("It is organised in directories: to change directory use DIRNAME->cd()")
print('To draw a TGraph: OBJNAME->Draw()')


