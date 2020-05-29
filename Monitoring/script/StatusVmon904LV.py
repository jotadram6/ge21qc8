#!/usr/bin/env python

import cx_Oracle
import ROOT
import os
import time
from operator import itemgetter
from datetime import datetime
from array import array

#in the DB the DeltaV between pins are saved, not the V from ground

#time interval considered to extract the data
#sta_period = "'2018-07-17 08:00:01'"
#end_period = "'2018-07-17 14:30:01'"

sta_period = raw_input("Insert UTC start time in format YYYY-MM-DD HH:mm:ss\n")
type(sta_period)
end_period = raw_input("Insert UTC end time in format YYYY-MM-DD HH:mm:ss\n")
type(end_period)

#to remove ' and space
start=sta_period.replace(" ", "_")
end=end_period.replace(" ", "_")
start=start.replace(":", "-")
end=end.replace(":", "-")

#add ' at beginning and end to have the date in the format for the query
sta_period = "'" + sta_period + "'"
end_period = "'" + end_period + "'"

# I also include some root histogram in case some plot are needed..
fileName = "QC8_LV_monitor_UTC_start_"+start+"_end_"+end+".root"
f1=ROOT.TFile( fileName,"RECREATE")

#divide the monitoring period in three periods
#the limits between period 1 and 2 is defined by the first  date in which we changed the mapping of LV: 03-Apr-2019 13.10
#the limits between period 2 and 3 is defined by the second date in which we changed the mapping of LV: 15-Apr-2019 08:15

#in the first date only the HV mapping is changed, but the datapoints change. So we define three periods even if we have
#only two mappings

mappingChangeDate = []
firstMappingChange 	= datetime( 2019, 04, 03, 13, 10, 00 )
secondMappingChange 	= datetime( 2019, 04, 15,  8, 15, 00 )
mappingChangeDate.append( firstMappingChange )
mappingChangeDate.append( secondMappingChange )

startDate = datetime(int(start[:4]), int(start[5:7]), int(start[8:10]), int(start[11:13]), int(start[14:16]), int(start[17:]) )
endDate	  = datetime(int(end[:4]), int(end[5:7]), int(end[8:10]), int(end[11:13]), int(end[14:16]), int(end[17:]) )

periodBool = []
numberOfMaps = 3 #date 10-May-2019

#create a vector of numbers full of zero: numbero of zero is equal to the number of maps
for mapIdx in range(numberOfMaps):
	periodBool.append(0)

startIdx = -1
endIdx = -1
#find the index for mapping start and end

if startDate < mappingChangeDate[0]:
	startIdx = 0
elif (startDate > mappingChangeDate[0] and startDate < mappingChangeDate[1]):
	startIdx = 1
elif startDate > mappingChangeDate[1]:
	startIdx = 2
#print "startIdx:", startIdx

if endDate < mappingChangeDate[0]:
	endIdx = 0
elif endDate > mappingChangeDate[0] and endDate < mappingChangeDate[1]:
	endIdx = 1
elif endDate > mappingChangeDate[1]:
	endIdx = 2
#print "endIdx:", endIdx

#decide which periods are active from start and stop index
periodIdx = startIdx
while periodIdx <= endIdx:
	periodBool[ periodIdx ] = 1
	periodIdx = periodIdx + 1

#print periodBool

#validity limits needed for the query of ELEMENT_ID
#find the first one and the last one in the bool vector
firstOne = periodBool.index(1)
lastOne  = max(loc for loc, val in enumerate(periodBool) if val == 1)
#print "firstOne", firstOne
#print "lastOne", lastOne

sinceDelimiter = []
tillDelimiter = []

#three cases:
#only one period section used
#first period used
if firstOne == 0 and lastOne == 0:
	#we don't need the since delimiter
	tillDelimiter.append( mappingChangeDate[0] )
#last period used
elif firstOne == (len(periodBool) - 1) and lastOne == (len(periodBool)-1):
	#we don't need the end delimiter
	sinceDelimiter.append( mappingChangeDate[ len(periodBool)-2 ] )
#more than one section used starting from the first time region
#first and others but not the last used
elif firstOne == 0 and lastOne > 0 and lastOne < (len(periodBool)-1):
	#we need N till and N-1 since
	for idxSince in range( periodBool.count(1)-1):
		sinceDelimiter.append( mappingChangeDate[ idxSince ] )
	for idxTill in range( periodBool.count(1) ):
		tillDelimiter.append( mappingChangeDate[ idxTill ] )
#last and others used but not the first
elif firstOne > 0 and firstOne < (len(periodBool)-1) and lastOne == (len(periodBool)-1):
	for idxSince in range( periodBool.count(1) ):
		sinceDelimiter.append( mappingChangeDate[ idxSince + firstOne-1 ]  )
	for idxTill in range( periodBool.count(1)-1 ):
		tillDelimiter.append( mappingChangeDate[ idxTill + firstOne ] )
#all used ( first, mid and last )
elif firstOne == 0 and lastOne == ( len(periodBool)-1 ):
	sinceDelimiter = mappingChangeDate
	tillDelimiter  = mappingChangeDate
#only mid used
elif firstOne > 0 and firstOne < (len(periodBool)-1) and lastOne > 0 and lastOne < (len(periodBool)-1):
	for idxSince in range( periodBool.count(1) ):
		sinceDelimiter.append( mappingChangeDate[ firstOne + idxSince -1 ] )
	for idxTill in range( periodBool.count(1) ):
		tillDelimiter.append( mappingChangeDate[ firstOne + idxTill ] )


#print "sinceDelimiter", sinceDelimiter
#print "tillDelimiter", tillDelimiter




#boards must be 15 (from 0 to 14) because in the DB the board 15 is not declared
#in the LV we have only board 00 04 08 12 16
insertedBoards = 5
insertedChannels = 6

#insert the list of chambers in DCS convention
howManyChambers = raw_input("Tell me how many chambers you are using\n")
type(howManyChambers)

chamberList = []
print ("Tell me the chambers' names in DCS convention (for example 1-2-Top)")

for cont in range(int(howManyChambers)):
	chamberName = raw_input("Chamber: ")
	type(chamberName)
	chamberName = chamberName.replace("-","_")
	chamberList.append(chamberName)

#print(chamberList)

indexB=0
indexC=0
counter=0
mismatch=0
mismatch2=0
mismatch3=0
mismatch4=0
mismatch5=0
mismatch6=0
mismatch7=0
mismatch8=0
mismatch9=0
mismatch10=0
mismatch11=0
#bits of the error status
nBit = 16

ImonTh1List = []
VmonTh1List = []
SmonTh1List = []
SmonMeaningListList = [] #a list of lists (one list for each channel of each board considered)
			 # each of the inner lists has the meaning of the status

ImonTgraph1List = []
VmonTgraph1List = []
SmonTgraph1List = []

ChamberMapList = [ "1_1_Top", "1_1_Bot", "1_2_Top", "1_2_Bot", "1_3_Top", "1_3_Bot", "2_1_Top", "2_1_Bot", "2_2_Top", "2_2_Bot", "2_3_Top", "2_3_Bot", "3_1_Top", "3_1_Bot", "3_2_Top", "3_2_Bot", "3_3_Top", "3_3_Bot", "4_1_Top", "4_1_Bot", "4_2_Top", "4_2_Bot", "4_3_Top", "4_3_Bot", "5_1_Top", "5_1_Bot", "5_2_Top", "5_2_Bot", "5_3_Top", "5_3_Bot" ]

MainframeMapList1 = [ "904_Shared_mainframe", "904_Shared_mainframe", "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe" ]

BranchControllerList1 = [ "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00" ]

BoardMapList1 = [ "00", "00", "00", "00", "00", "00", "04", "04", "04", "04", "04", "04", "08", "08", "08", "08", "08", "08", "12", "12", "12", "12", "12", "12", "16", "16", "16", "16", "16", "16" ]


MainframeMapList2 = [ "904_Shared_mainframe", "904_Shared_mainframe", "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe" ]

BranchControllerList2 = [ "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00", "00" ]

BoardMapList2 = [ "00", "00", "00", "00", "00", "00", "04", "04", "04", "04", "04", "04", "08", "08", "08", "08", "08", "08", "12", "12", "12", "12", "12", "12", "16", "16", "16", "16", "16", "16" ]


MainframeMapList3 = [ "904_HV_mainframe", "904_HV_mainframe", "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe" ]

BranchControllerList3 = [ "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14" ]

BoardMapList3 = [ "00", "00", "00", "00", "00", "00", "04", "04", "04", "04", "04", "04", "08", "08", "08", "08", "08", "08", "12", "12", "12", "12", "12", "12", "16", "16", "16", "16", "16", "16" ]


ChannelMapList = [ "000", "001", "002", "003", "004", "005", "000", "001", "002", "003", "004", "005", "000", "001", "002", "003", "004", "005", "000", "001", "002", "003", "004", "005", "000", "001", "002", "003", "004", "005" ]

listAllMainframeMappings = []
listAllMainframeMappings.append( MainframeMapList1 )
listAllMainframeMappings.append( MainframeMapList2 )
listAllMainframeMappings.append( MainframeMapList3 )

listAllBranchControllerMappings = []
listAllBranchControllerMappings.append( BranchControllerList1 )
listAllBranchControllerMappings.append( BranchControllerList2 )
listAllBranchControllerMappings.append( BranchControllerList3 )

listAllBoardMappings = []
listAllBoardMappings.append( BoardMapList1 )
listAllBoardMappings.append( BoardMapList2 )
listAllBoardMappings.append( BoardMapList3 )



for indexB in range(len(chamberList)): #loop on the selected boards
	#create directories
        firstDir = f1.mkdir("LV_Chamber"+chamberList[indexB])
        firstDir.cd()

	#print("indexB="+str(indexB))
	#th1List.append(NameB)

	#set bin size
	IMinLV = -10		#A
	IMaxLV = 10		#A
	IResolutionLV = 0.01	#A
	INbinLV = int((IMaxLV-IMinLV)/IResolutionLV)

        Imonh1 = ROOT.TH1F("LV_ImonChamber"+chamberList[indexB]+"_TH1","LV_ImonChamber"+chamberList[indexB]+"_TH1",INbinLV,IMinLV,IMaxLV)
        Imonh1.GetXaxis().SetTitle("I [A]")
        Imonh1.GetYaxis().SetTitle("counts")


	#set bin size
        VMinLV = -50		#V
        VMaxLV = 200		#V
        VResolutionLV = 0.005	#V
        VNbinLV = int((VMaxLV-VMinLV)/VResolutionLV)

        Vmonh1 = ROOT.TH1F("LV_VmonChamber"+chamberList[indexB]+"_TH1","LV_VmonChamber"+chamberList[indexB]+"_TH1",VNbinLV,VMinLV,VMaxLV)
        Vmonh1.GetXaxis().SetTitle("V [V]")
        Vmonh1.GetYaxis().SetTitle("counts")

        ImonTh1List.append(Imonh1)
        VmonTh1List.append(Vmonh1)

	#print(len(ChamberMapList))
	#print(len(BoardMapList))

	#I define 5 status categories
        #type 0: off (bit 0)

        #Smonh1 = ROOT.TH1F("LV_StatusChamber"+chamberList[indexB]+"_TH1","LV_StatusChamber"+chamberList[indexB]+"_TH1",20,-1, 19)
        Smonh1 = ROOT.TH1F("LV_StatusChamber"+chamberList[indexB]+"_TH1","LV_StatusChamber"+chamberList[indexB]+"_TH1", 65536, 0, 65536)
        #Smonh1.GetXaxis().SetTitle("Status cathegory")
        Smonh1.GetXaxis().SetTitle("Status code")
        Smonh1.GetYaxis().SetTitle("counts")

        SmonTh1List.append(Smonh1)

        SmonMeaningList = []


	#find the number of the board from the name
	#create a query string for the maps that must me used
        #usedMaps tells the number of used maps for the period between start and end, not the number of maps charged in the code (said instead by numberOfMaps)
        usedMaps = periodBool.count(1)
        imonNameList	= []
        vmonNameList 	= []
        statusNameList 	= []
	contMaps = 0
        indexTB = 0
        for form1 in ChamberMapList:
        	form2 = chamberList[indexB] #remove -Top or -Bot from the name to make a match
        	#print("form1 e form2:", form1, form2)
        	isMatch = False
        	if form2 == form1:
        		isMatch = True
        		board = "Board"+listAllBoardMappings[startIdx + contMaps][indexTB]#board in the format board04
			channel = "channel"+ChannelMapList[indexTB]
        	if isMatch == True:
        		break
        	indexTB = indexTB +1
        #print(indexTB)
        #print("Board", board, "channel", channel)

	#print ("counter:", counter)

        # this 3 parameters can be input to the script
	#imon_name="'cms_gem_dcs_1:CAEN/904_HV_mainframe/branchController14/easyCrate0/easy"+board+"/"+channel+".actual.iMon'"
	#vmon_name="'cms_gem_dcs_1:CAEN/904_HV_mainframe/branchController14/easyCrate0/easy"+board+"/"+channel+".actual.vMon'"
	#status_name="'cms_gem_dcs_1:CAEN/904_HV_mainframe/branchController14/easyCrate0/easy"+board+"/"+channel+".actual.status'"


        #Three maps are inserted so we have to choose the right one
        for contMaps in range( usedMaps ):
        	imon_name="'cms_gem_dcs_1:CAEN/"+listAllMainframeMappings[startIdx + contMaps][indexTB]+"/branchController"+listAllBranchControllerMappings[startIdx + contMaps][indexTB]+"/easyCrate0/easyBoard"+listAllBoardMappings[startIdx + contMaps][indexTB]+"/"+channel+".actual.iMon'"
        	#this string has to be saved only one time
        	if imonNameList.count( imon_name ) == 1:
        		continue
        	imonNameList.append(imon_name)

        for contMaps in range( usedMaps ):
        	vmon_name="'cms_gem_dcs_1:CAEN/"+listAllMainframeMappings[startIdx + contMaps][indexTB]+"/branchController"+listAllBranchControllerMappings[startIdx + contMaps][indexTB]+"/easyCrate0/easyBoard"+listAllBoardMappings[startIdx + contMaps][indexTB]+"/"+channel+".actual.vMon'"
        	#this string has to be saved only one time
        	if vmonNameList.count( vmon_name ) == 1:
        		continue
        	vmonNameList.append(vmon_name)

        for contMaps in range( usedMaps ):
        	status_name="'cms_gem_dcs_1:CAEN/"+listAllMainframeMappings[startIdx + contMaps][indexTB]+"/branchController"+listAllBranchControllerMappings[startIdx + contMaps][indexTB]+"/easyCrate0/easyBoard"+listAllBoardMappings[startIdx + contMaps][indexTB]+"/"+channel+".actual.status'"
        	#this string has to be saved only one time
        	if statusNameList.count( status_name ) == 1:
        		continue
        	statusNameList.append(status_name)
        	#print (listAllMainframeMappings[startIdx + contMaps][indexTB], " ", listAllBoardMappings[startIdx + contMaps][indexTB], " ", channel)

        print "\nimonNameList", imonNameList
        print "vomnNameList", vmonNameList
        print "statusNameList", statusNameList



        #this is the connection to DB, and the contact point to DB should be input to the script
        #to avoid the have explicitely the pw in the script (lo schema e' il primo campo, la parte a destra e' il server)
        db_cond = os.environ["GEM_DEVELOPMENT_DB_COND"]
	    db_name = os.environ["GEM_DEVELOPMENT_DB_NAME"]

	    db = cx_Oracle.connect(db_cond+db_name) # production DB
	    cur = db.cursor()

        #ELEMENTS is the table in the DB
        #ELEMENT_ID and also ELEMENT_NAME are fields in the table
        #where ELEMENT_NAME is the row of the column ELEMENT_ID

        # this a tytpical query to extract the ELEMENT_ID from the  ELEMENT_NAME.
        # iMon element
		imonIdList = []
        vmonIdList = []
        statusIdList = []
        sinceImonList = []
        sinceVmonList = []
        sinceStatusList = []
        tillImonList = []
        tillVmonList = []
        tillStatusList = []
        imonIDsort = []
        vmonIDsort = []
        statusIDsort = []
        for elementIdIdxImon in range( len(imonNameList) ):
        	# iMon element
        	query = "select ELEMENT_ID, VALID_SINCE, VALID_TILL from ELEMENTS_ALL where ELEMENT_NAME="+imonNameList[ elementIdIdxImon ]
        	cur.execute(query)
        	curImonID = cur
        	for result in curImonID:
        		#check if the lenght of the since and till array are the same
        		#VALID_SINCE is always found in the DB
        		#VALID_TILL is not found for the most recent ELEMENT_ID
        		#in this case I put the end date as its VALID_TILL
        		imonID = result[0]
        		imonIDsince = result[1]
        		imonIDtill = result[2]

        		if imonIDtill is None:
        			imonIDtill = datetime( 2050, 04, 15,  8, 15, 00 )

        		imonIdList.append( imonID )
        		sinceImonList.append( imonIDsince )
        		tillImonList.append( imonIDtill )

        		internalSortList = ( imonID, imonIDsince, imonIDtill )
        		imonIDsort.append( internalSortList )

        imonIDsort = sorted(imonIDsort, key=lambda elementIDimon: elementIDimon[1])

        print "imonIDsort ", imonIDsort

        #I have to reject the ID that are not necessary for our time period
        #write 1 in the period observed in the monitoring
        firstOneImon = -1
        lastOneImon = -1
        imonIDBool = []
        for idxIdImon in range( len(imonIDsort) ):
        	imonIDBool.append(0)
        	#firstOneImon
        	if startDate > imonIDsort[idxIdImon][1] and startDate < imonIDsort[idxIdImon][2]:
        		firstOneImon = idxIdImon
        	#lastOneImon
        	if endDate > imonIDsort[idxIdImon][1] and endDate < imonIDsort[idxIdImon][2]:
        		lastOneImon = idxIdImon

        #print "firstOneImon", firstOneImon
        #print "lastOneImon", lastOneImon

        #I write one in all the active periods
        for timesOneFill in range( lastOneImon - firstOneImon + 1 ):
        	imonIDBool[ timesOneFill + firstOneImon ] = 1

        #print "imonIDsort ", imonIDsort

        imonIDToUse = []
        for idxIdImon in range( len(imonIDsort) ):
        	if imonIDBool[ idxIdImon ] == 1:
        		imonIDToUse.append( imonIDsort[ idxIdImon ][0] )

        #delete duplicate IDs
        imonIDToUse = list(dict.fromkeys( imonIDToUse ))
        print "imonIDToUse", imonIDToUse



        # vMon element
	#find elementID to use for vmon
        for elementIdIdxVmon in range( len(vmonNameList) ):
        	# vMon element
        	query = "select ELEMENT_ID, VALID_SINCE, VALID_TILL from ELEMENTS_ALL where ELEMENT_NAME="+vmonNameList[ elementIdIdxVmon ]
        	cur.execute(query)
        	curVmonID = cur
        	for result in curVmonID:
        		#check if the lenght of the since and till array are the same
        		#VALID_SINCE is always found in the DB
        		#VALID_TILL is not found for the most recent ELEMENT_ID
        		#in this case I put the end date as its VALID_TILL
        		vmonID = result[0]
        		vmonIDsince = result[1]
        		vmonIDtill = result[2]

        		if vmonIDtill is None:
        			vmonIDtill = datetime( 2050, 04, 15,  8, 15, 00 )

        		vmonIdList.append( vmonID )
        		sinceVmonList.append( vmonIDsince )
        		tillVmonList.append( vmonIDtill )

        		internalSortList = ( vmonID, vmonIDsince, vmonIDtill )
        		vmonIDsort.append( internalSortList )

        vmonIDsort = sorted(vmonIDsort, key=lambda elementIDvmon: elementIDvmon[1])

        print "vmonIDsort ", vmonIDsort

        #I have to reject the ID that are not necessary for our time period
        #write 1 in the period observed in the monitoring
        firstOneVmon = -1
        lastOneVmon = -1
        vmonIDBool = []
        for idxIdVmon in range( len(vmonIDsort) ):
        	vmonIDBool.append(0)
        	#firstOneVmon
        	if startDate > vmonIDsort[idxIdVmon][1] and startDate < vmonIDsort[idxIdVmon][2]:
        		firstOneVmon = idxIdVmon
        	#lastOneVmon
        	if endDate > vmonIDsort[idxIdVmon][1] and endDate < vmonIDsort[idxIdVmon][2]:
        		lastOneVmon = idxIdVmon

        #print "firstOneVmon", firstOneVmon
        #print "lastOneVmon", lastOneVmon

        #I write one in all the active periods
        for timesOneFill in range( lastOneVmon - firstOneVmon + 1 ):
        	vmonIDBool[ timesOneFill + firstOneVmon ] = 1

        #print "vmonIDsort ", vmonIDsort

        vmonIDToUse = []
        for idxIdVmon in range( len(vmonIDsort) ):
        	if vmonIDBool[ idxIdVmon ] == 1:
        		vmonIDToUse.append( vmonIDsort[ idxIdVmon ][0] )

        #delete duplicate IDs
        vmonIDToUse = list(dict.fromkeys( vmonIDToUse ))
        print "vmonIDToUse", vmonIDToUse



        # status element
	#find elementID to use for status
        for elementIdIdxStatus in range( len(statusNameList) ):
        	# status element
        	query = "select ELEMENT_ID, VALID_SINCE, VALID_TILL from ELEMENTS_ALL where ELEMENT_NAME="+statusNameList[ elementIdIdxStatus ]
        	cur.execute(query)
        	curStatusID = cur
        	for result in curStatusID:
        		#check if the lenght of the since and till array are the same
        		#VALID_SINCE is always found in the DB
        		#VALID_TILL is not found for the most recent ELEMENT_ID
        		#in this case I put the end date as its VALID_TILL
        		statusID = result[0]
        		statusIDsince = result[1]
        		statusIDtill = result[2]

        		if statusIDtill is None:
        			statusIDtill = datetime( 2050, 04, 15,  8, 15, 00 )

        		statusIdList.append( statusID )
        		sinceStatusList.append( statusIDsince )
        		tillStatusList.append( statusIDtill )

        		internalSortList = ( statusID, statusIDsince, statusIDtill )
        		statusIDsort.append( internalSortList )

        statusIDsort = sorted(statusIDsort, key=lambda elementIDstatus: elementIDstatus[1])

        #print "statusIDsort ", statusIDsort

        #I have to reject the ID that are not necessary for our time period
        #write 1 in the period observed in the monitoring
        firstOneStatus = -1
        lastOneStatus = -1
        statusIDBool = []
        for idxIdStatus in range( len(statusIDsort) ):
        	statusIDBool.append(0)
        	#firstOneStatus
        	if startDate > statusIDsort[idxIdStatus][1] and startDate < statusIDsort[idxIdStatus][2]:
        		firstOneStatus = idxIdStatus
        	#lastOneStatus
        	if endDate > statusIDsort[idxIdStatus][1] and endDate < statusIDsort[idxIdStatus][2]:
        		lastOneStatus = idxIdStatus

        #print "firstOneStatus", firstOneStatus
        #print "lastOneStatus", lastOneStatus

        #I write one in all the active periods
        for timesOneFill in range( lastOneStatus - firstOneStatus + 1 ):
        	statusIDBool[ timesOneFill + firstOneStatus ] = 1

        #print "statusIDsort ", statusIDsort

        statusIDToUse = []
        for idxIdStatus in range( len(statusIDsort) ):
        	if statusIDBool[ idxIdStatus ] == 1:
        		statusIDToUse.append( statusIDsort[ idxIdStatus ][0] )

        #delete duplicate IDs
        statusIDToUse = list(dict.fromkeys( statusIDToUse ))
        print "statusIDToUse", statusIDToUse


	#USEFUL for DEBUG of IDs
	#do the query to fill the Histos with I
        #fill current
	#print "imonNameList ", imonIdList
        #print "sinceImonList ", sinceImonList
        #print "tillImonList ", tillImonList
        #print "vmonNameList ", vmonIdList
        #print "sinceVmonList ", sinceVmonList
        #print "tillVmonList ", tillVmonList
        #print "statusNameList ", statusIdList
        #print "sinceStatusList ", sinceStatusList
        #print "tillStatusList ", tillStatusList
        #print "IMON_ID ", imon_id," VMON_ID ", vmon_id," STATUS_ID ", status_id

        #do the query to fill the Histos with I
        if len( imonIDToUse ) == 0:
        	print "ERROR: len(imonIDToUse)=0"
        	file = open("HVErr.log", "w")
                file.write("ERROR: len(imonIDToUse)=0")
                file.close()
        elif len( imonIDToUse ) == 1:
        	query = "select TS,VALUE_NUMBER from EVENTHISTORY where ELEMENT_ID = "+str(imonIDToUse[0])+" and TS > to_date ("+sta_period+",'YYYY-MM-DD HH24:MI:SS') and TS < to_date ("+end_period+",'YYYY-MM-DD HH24:MI:SS')"
        elif len( imonIDToUse ) > 1:
        	query = "select TS,VALUE_NUMBER from EVENTHISTORY where ("
        	for strIdx in range ( len( imonIDToUse ) ):
        		query = query + "ELEMENT_ID = "+str( imonIDToUse[ strIdx ] )
        		if ( len( imonIDToUse ) > 1 and strIdx < ( len( imonIDToUse ) - 1 ) ):
        			query = query + " or "
        	query = query + ") and TS > to_date ("+sta_period+",'YYYY-MM-DD HH24:MI:SS') and TS < to_date ("+end_period+",'YYYY-MM-DD HH24:MI:SS')"

        print query
        cur.execute(query)
        curimon=cur
        imonRec=[]
        imonOnlyT = array ( 'd' )
        imonOnlyI = array ( 'd' )
        imonOnlyTDate = array ( 'd' )
        contatoreI = 0
        isNotEmptyImon = False
        for result in curimon:
           #start ts: the first read Ts
           if contatoreI == 0:
                #print("dentroI", contatoreI)
                startTsImon = result[0]
           isNotEmptyImon = True
           ImonTh1List[counter].Fill(result[1])
           imonRec+=[(result[0],result[1],"Imon")]
           #print("startTsImon",startTsImon, "present ts", result[0])
           currentTsImon = result[0] #put the present TS value in currentTsImon
           #print("TSImon", currentTsImon-startTsImon)
           #micro=(currentTsImon-startTsImon).microseconds
           tot_secondsImon = (currentTsImon-startTsImon).total_seconds()
           #print("micro=", micro)
           #print("totsecImon=", tot_secondsImon)
           tot_milliImon=tot_secondsImon*1000.
           #print("totmilliImon=", tot_milliImon)
           #print(imonRec[contatoreI])

           #lists with only I and t
           imonOnlyT.append(tot_secondsImon)
           imonOnlyI.append(result[1])

	   imonOnlyTDate.append(4.) #an any float number

           #print currentTsImon
           year = str(currentTsImon)[0:4]
           #print year
           month = str(currentTsImon)[5:7]
           #print month
           day = str(currentTsImon)[8:10]
           #print day
           hour = str(currentTsImon)[11:13]
           #print hour
           minute = str(currentTsImon)[14:16]
           #print minute
           second = str(currentTsImon)[17:19]
           #print second
           micro = str(currentTsImon)[20:26]
           #print micro

           #longString = ROOT.TString( year+"-"+month+"-"+day+" "+hour+":"+minute+":"+second )
           longList =str( year+"-"+month+"-"+day+" "+hour+":"+minute+":"+second )

           da1 = ROOT.TDatime( longList )

           imonOnlyTDate[-1] = da1.Convert()
           floatMicro = "0."+micro
           #print floatMicro
           imonOnlyTDate[contatoreI] = imonOnlyTDate[contatoreI] + float(floatMicro) #add microseconds to times (all ending with .0 because of Convert)
           contatoreI=contatoreI+1

        ImonTh1List[counter].Write()

	#create a TGraph for Imon (I vs time)
        if len(imonOnlyT) != len(imonOnlyI):
        	print ("filling vector different lenght")
		file = open("LVErr.log", "w")
                file.write("ERROR: imonOnlyT and imonOnlyI have different lenght")
                file.close()

	if len(imonOnlyTDate) != len(imonOnlyI):
        	print ("filling vector different lenght")
        	file = open("LVErr.log", "w")
                file.write("ERROR: imonOnlyTDate and imonOnlyI have different lenght")
                file.close()

        #sort the array of imonOnlyT and the imonOnlyV
        #in the case the query is not executed in order (negative times)
        #pair the time with status and the meaning list
        SortListI = []
        for sortCountI in range(len(imonOnlyT)):
        	internalListI = (imonOnlyT[sortCountI], imonOnlyI[sortCountI], imonOnlyTDate[sortCountI])
        	SortListI.append(internalListI)

        #print(SortListI)
        SortListI = sorted(SortListI, key=lambda elementI: elementI[0])
        #print(SortListI)

        for refillI in range(len(imonOnlyT)):
        	imonOnlyT[refillI]=SortListI[refillI][0]
        	imonOnlyI[refillI]=SortListI[refillI][1]
        	imonOnlyTDate[refillI]=SortListI[refillI][2]

        #print(imonOnlyT)
        #print(imonOnlyI)

        #rescale the negative times

	if (len(imonOnlyT))==0:
        	print("imonOnlyT lenght", len(imonOnlyT))
        	print("------------------------------------------------------------------")
        	print("ERROR: there are no I current data for chamber "+ chamberList[indexB])
        	print("------------------------------------------------------------------")

        	file = open("LVErr.log", "w")
        	file.write("ERROR: there are no I current data for chamber "+ chamberList[indexB])
        	file.close()

		#counter = counter + 1
        	#I put an anomalous current value if there is no data
	        imonOnlyT.append(0.)
                imonOnlyTDate.append(0.)
	        imonOnlyI.append(-1000000000)

        	#continue


        negativeStartI = imonOnlyT[0]
        if imonOnlyT[0] < 0:
        	for iterTimeI in range(len(imonOnlyT)):
        		imonOnlyT[iterTimeI] = imonOnlyT[iterTimeI] + negativeStartI*(-1)
        #print(imonOnlyT)


        #Imontg1 = ROOT.TGraph(n,x,y); x and y is the name of arrays with numbers of time and I
        Imontg1 = ROOT.TGraph(len(imonOnlyT),imonOnlyTDate,imonOnlyI)
        Imontg1.SetLineColor(2)
        Imontg1.SetLineWidth(4)
        Imontg1.SetMarkerColor(4)
        Imontg1.SetMarkerStyle(21)
        Imontg1.SetMarkerSize(1)
        Imontg1.SetName("LV_ImonChamber"+chamberList[indexB]+"_UTC_time")
        Imontg1.SetTitle("LV_ImonChamber"+chamberList[indexB]+"_UTC_time")
        #Imontg1.GetXaxis().SetTitle("time [s]")
        Imontg1.GetYaxis().SetTitle("Imon "+chamberList[indexB]+" [A]")
        #Imontg1.Draw("ACP")
	Imontg1.GetXaxis().SetTimeDisplay(1)
        Imontg1.GetXaxis().SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00")
        Imontg1.GetXaxis().SetLabelOffset(0.025)

        #ImonTgraph1List += [Imontg1]
        #ImonTgraph1List.append(Imontg1)
        #ImonTgraph1List[counter].Write()
	Imontg1.Write()

	#fill Voltage
	#do the query to fill the Histos with V
        if len( vmonIDToUse ) == 0:
        	print "ERROR: len(vmonIDToUse)=0"
        	file = open("HVErr.log", "w")
                file.write("ERROR: len(vmonIDToUse)=0")
                file.close()
        elif len( vmonIDToUse ) == 1:
        	query = "select TS,VALUE_NUMBER from EVENTHISTORY where ELEMENT_ID = "+str(vmonIDToUse[0])+" and TS > to_date ("+sta_period+",'YYYY-MM-DD HH24:MI:SS') and TS < to_date ("+end_period+",'YYYY-MM-DD HH24:MI:SS')"
        elif len( vmonIDToUse ) > 1:
        	query = "select TS,VALUE_NUMBER from EVENTHISTORY where ("
        	for strIdx in range ( len( vmonIDToUse ) ):
        		query = query + "ELEMENT_ID = "+str( vmonIDToUse[ strIdx ] )
        		if ( len( vmonIDToUse ) > 1 and strIdx < ( len( vmonIDToUse ) - 1 ) ):
        			query = query + " or "
        	query = query + ") and TS > to_date ("+sta_period+",'YYYY-MM-DD HH24:MI:SS') and TS < to_date ("+end_period+",'YYYY-MM-DD HH24:MI:SS')"

        print query
        cur.execute(query)
        #result[0] is the present TS
        #result[1] is the voltage
        curvmon=cur
        vmonRec=[]
        vmonOnlyT = array ( 'd' )
        vmonOnlyV = array ( 'd' )
        vmonOnlyTDate = array ( 'd' )
        contatoreV = 0
        isNotEmptyVmon = False
        for result in curvmon:
           #start ts: the first read Ts
           if contatoreV == 0:
           	#print("dentroV", contatoreV)
           	startTsVmon = result[0]
           isNotEmptyVmon = True
           VmonTh1List[counter].Fill(result[1])
           vmonRec+=[(result[0],result[1],"Vmon")]
           #print("startTsVmon",startTsVmon, "present ts", result[0])
           currentTsVmon = result[0] #put the present TS value in currentTsVmon
           #print("TSVmon", currentTsVmon-startTsVmon)
           #micro=(currentTsVmon-startTsVmon).microseconds
           tot_secondsVmon = (currentTsVmon-startTsVmon).total_seconds()
           #print("micro=", micro)
           #print("totsecVmon=", tot_secondsVmon)
           tot_milliVmon=tot_secondsVmon*1000.
           #print("totmilliVmon=", tot_milliVmon)
           #print(vmonRec[contatoreV])

           #lists with only V and t
	   vmonOnlyT.append(tot_secondsVmon)
           vmonOnlyV.append(result[1])

	   vmonOnlyTDate.append(4.) #an any float number

           #print currentTsVmon
           year = str(currentTsVmon)[0:4]
           #print year
           month = str(currentTsVmon)[5:7]
           #print month
           day = str(currentTsVmon)[8:10]
           #print day
           hour = str(currentTsVmon)[11:13]
           #print hour
           minute = str(currentTsVmon)[14:16]
           #print minute
           second = str(currentTsVmon)[17:19]
           #print second
           micro = str(currentTsVmon)[20:26]
           #print micro

           #longString = ROOT.TString( year+"-"+month+"-"+day+" "+hour+":"+minute+":"+second )
           longList =str( year+"-"+month+"-"+day+" "+hour+":"+minute+":"+second )

           da1 = ROOT.TDatime( longList )

           vmonOnlyTDate[-1] = da1.Convert()
           floatMicro = "0."+micro
           #print floatMicro
           vmonOnlyTDate[contatoreV] = vmonOnlyTDate[contatoreV] + float(floatMicro) #add microseconds to times (all ending with .0 because of Convert)

	   contatoreV = contatoreV+1
        VmonTh1List[counter].Write()

	#create a TGraph for Imon (V vs time)
        if len(vmonOnlyT) != len(vmonOnlyV):
        	print ("filling vector different lenght")
		file = open("LVErr.log", "w")
                file.write("ERROR: vmonOnlyT and vmonOnlyV have different lenght")
                file.close()

	if len(vmonOnlyTDate) != len(vmonOnlyV):
        	print ("filling vector different lenght")
        	file = open("LVErr.log", "w")
                file.write("ERROR: vmonOnlyTDate and vmonOnlyV have different lenght")
                file.close()


        #sort the array of vmonOnlyT and the vmonOnlyV
        #in the case the query is not executed in order (negative times)
        #pair the time with status and the meaning list
        SortListV = []
        for sortCountV in range(len(vmonOnlyT)):
        	internalListV = (vmonOnlyT[sortCountV], vmonOnlyV[sortCountV], vmonOnlyTDate[sortCountV])
        	SortListV.append(internalListV)

        #print(SortListV)
        SortListV = sorted(SortListV, key=lambda elementV: elementV[0])
        #print(SortListV)

        for refillV in range(len(vmonOnlyT)):
        	vmonOnlyT[refillV]=SortListV[refillV][0]
        	vmonOnlyV[refillV]=SortListV[refillV][1]
        	vmonOnlyTDate[refillV]=SortListV[refillV][2]

        #print(vmonOnlyT)
        #print(vmonOnlyV)

        #rescale the negative times

	if (len(vmonOnlyT))==0:
        	print("vmonOnlyT lenght", len(vmonOnlyT))
        	print("------------------------------------------------------------------")
        	print("ERROR: there are no LV voltage data for chamber "+ chamberList[indexB])
        	print("------------------------------------------------------------------")

        	file = open("LVErr.log", "w")
        	file.write("ERROR: there are no LV voltage data for chamber "+ chamberList[indexB])
        	file.close()

		#counter = counter + 1
		#I put an anomalous voltage value if there is no data
                vmonOnlyT.append(0.)
                vmonOnlyTDate.append(0.)
                vmonOnlyV.append(-1000000000)

        	#continue


        negativeStartV = vmonOnlyT[0]
        if vmonOnlyT[0] < 0:
        	for iterTimeV in range(len(vmonOnlyT)):
        		vmonOnlyT[iterTimeV] = vmonOnlyT[iterTimeV] + negativeStartV*(-1)
        #print(vmonOnlyT)


        #Vmontg1 = ROOT.TGraph(n,x,y); x and y is the name of arrays with numbers of time and V
        Vmontg1 = ROOT.TGraph(len(vmonOnlyT),vmonOnlyTDate,vmonOnlyV)
        Vmontg1.SetLineColor(2)
        Vmontg1.SetLineWidth(4)
        Vmontg1.SetMarkerColor(4)
        Vmontg1.SetMarkerStyle(21)
        Vmontg1.SetMarkerSize(1)
        Vmontg1.SetName("LV_VmonChamber"+chamberList[indexB]+"_UTC_time")
        Vmontg1.SetTitle("LV_VmonChamber"+chamberList[indexB]+"_UTC_time")
        #Vmontg1.GetXaxis().SetTitle("time [s]")
        Vmontg1.GetYaxis().SetTitle("Vmon "+chamberList[indexB]+" [V]")
        #Vmontg1.Draw("ACP")
	Vmontg1.GetXaxis().SetTimeDisplay(1)
        Vmontg1.GetXaxis().SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00")
        Vmontg1.GetXaxis().SetLabelOffset(0.025)

        #VmonTgraph1List += [Vmontg1]
        #VmonTgraph1List.append(Vmontg1)
        #VmonTgraph1List[counter].Write()
	Vmontg1.Write()

	#LV boards (CAEN A3016 o A3016 HP)
        #Bit 0: ON/OFF
        #Bit 1: dont care
        #Bit 2: dont care
        #Bit 3: OverCurrent
        #Bit 4: OverVoltage
        #Bit 5: UnderVoltage
        #Bit 6: dont care
        #Bit 7: Over HVmax
        #Bit 8: dont care
        #Bit 9: Internal Trip
        #Bit 10: Calibration Error
        #Bit 11: Unplugged
        #Bit 12: dont care
        #Bit 13: OverVoltage Protection
        #Bit 14: Power Fail
        #Bit 15: Temperature Error

	#status query
	#do the query to fill the Histos with status
        if len( statusIDToUse ) == 0:
        	print "ERROR: len(statusIDToUse)=0"
        	file = open("HVErr.log", "w")
                file.write("ERROR: len(statusIDToUse)=0")
                file.close()
        elif len( statusIDToUse ) == 1:
        	query = "select TS,VALUE_NUMBER from EVENTHISTORY where ELEMENT_ID = "+str(statusIDToUse[0])+" and TS > to_date ("+sta_period+",'YYYY-MM-DD HH24:MI:SS') and TS < to_date ("+end_period+",'YYYY-MM-DD HH24:MI:SS')"
        elif len( statusIDToUse ) > 1:
        	query = "select TS,VALUE_NUMBER from EVENTHISTORY where ("
        	for strIdx in range ( len( statusIDToUse ) ):
        		query = query + "ELEMENT_ID = "+str( statusIDToUse[ strIdx ] )
        		if ( len( statusIDToUse ) > 1 and strIdx < ( len( statusIDToUse ) - 1 ) ):
        			query = query + " or "
        	query = query + ") and TS > to_date ("+sta_period+",'YYYY-MM-DD HH24:MI:SS') and TS < to_date ("+end_period+",'YYYY-MM-DD HH24:MI:SS')"

        print query
        cur.execute(query)
        curstat=cur
        statRec=[]
        smonOnlyT = array ( 'd' )
        smonOnlyS = array ( 'd' )
        smonDecimalStatus = array ( 'd' )
        smonOnlyTDate = array ( 'd' )
        smonOnlyTDateString = []
        smonOnlyBinStat = []
        smonOnlyMeaningStat = []
        contatoreS = 0
        isNotEmptySmon = False
        for result in curstat:
           #print ("nel for")
           if contatoreS == 0:
           	#print("dentroS", contatoreS)
           	startTsSmon = result[0]
           isNotEmptySmon = True
           statRec+=[(result[0],result[1],"Status")]
           #print ( "time", result[0], "status", result[1]  )
           #print ("statusbin", bin(int(result[1])))
           currentTsSmon = result[0] #put the present TS value in currentTsSmon
           decimalStatus = int(result[1])
           tot_secondsSmon = (currentTsSmon-startTsSmon).total_seconds()
           #use milli seconds to sort the time vector (at the end)
           tot_milliSmon=int(tot_secondsSmon*1000.)

           #print a 16 bit binary status number
           binStat = bin(int(result[1]))[2:] #to take away the 0b in front of the binary number
           lenStat = len(binStat)
           binStat = str(0) * (nBit - lenStat) + binStat
           binStat = "0b"+binStat
           smonOnlyBinStat.append(binStat)
           #print("time", result[0], "status bin", binStat )
           #print("status totseconds", tot_secondsSmon)

           extensibleStat = ""
           if len(binStat) != (nBit + 2) :
           	mismatch = mismatch + 1
           #masks for the bins
           if binStat == "0b0000000000000000": #these are binary numbers
           	channelStat = 0 #zero type of status (OFF)
           	#print("channel status OFF")
           	StatusMeaning = "OFF"
           	#extensibleStat = extensibleStat + StatusMeaning + " "
           	#print(StatusMeaning)

           if binStat == "0b0000000000000001": #these are binary numbers
           	channelStat = 1 #first type of status (ON)
           	#print("channel status ON")
           	StatusMeaning = "ON"
           	#extensibleStat = extensibleStat + StatusMeaning + " "
	   	#print(StatusMeaning)

	   cutBinStr = binStat[-1:]
           if cutBinStr == "0": #if I have OFF
           	extensibleStat = extensibleStat + "OFF" + " "
           elif cutBinStr == "1": #if I have OFF
           	extensibleStat = extensibleStat + "ON" + " "

           #bin produces a string (so the operation >> can be only made only on int)
           #I observe the bin number with bin(shift2)
           #I shift of one bit to delete the bit 0 from the string
	   removedBits = 0 - 1 #negative number
           shift2 = binStat[:removedBits]

           #print("binStat:", binStat, "shift2:", shift2 )
           if len(shift2) != (nBit + 2) + removedBits:
           	mismatch2 = mismatch2 + 1

	   #I have to remove bit 1 and 2 because they are not interesting
	   #len(shift2)-2    -2 because I want the last two bits
	   #print ( "shift2", shift2, "bin 1 and 2", shift2[len(shift2)-2:])

	   #remove bit 1 and 2 : second status cathegory even if it is written mismatch 3: I removed the bits
	   removedBits = removedBits - 2 #negative number
           shift3 = binStat[:removedBits]

           if len(shift3) != (nBit + 2) + removedBits:
           	mismatch3 = mismatch3 + 1

           #for the second status cathegory I need the last two bins of shift3
	   #print ( "shift3", shift3, "bit 3 4 5", shift3[len(shift3)-3:])
           if int(shift3[len(shift3)-3:]) > 0:
           	#print (shift3[len(shift3)-3:])
           	channelStat = 2 #second type of status
           	cutBinStr = shift3[len(shift3)-3:]
           	if cutBinStr[2] == "1": #if I have OVC
           		StatusMeaning = "OVC"
           		extensibleStat = extensibleStat + StatusMeaning + " "
           		#print(StatusMeaning)
           	if cutBinStr[1] == "1": #if I have OVV
           		StatusMeaning = "OVV"
           		extensibleStat = extensibleStat + StatusMeaning + " "
           		#print(StatusMeaning)
		if cutBinStr[0] == "1": #if I have UVV
                	StatusMeaning = "UVV"
                	extensibleStat = extensibleStat + StatusMeaning + " "
                	#print(StatusMeaning)

	   #remove bit 3 4 5
           removedBits = removedBits - 3 #negative number
           shift4 = binStat[:removedBits]

           if len(shift4) != (nBit + 2) + removedBits:
           	mismatch4 = mismatch4 + 1

	   #print ( "shift4", shift4, "bit 6", shift4[len(shift4)-1:])

	   #remove bit 6
           removedBits = removedBits - 1 #negative number
           shift5 = binStat[:removedBits]

           if len(shift5) != (nBit + 2) + removedBits:
           	mismatch5 = mismatch5 + 1

           #for the third status cathegory I need the last four bins of shift5
	   #I dont register the bit 8 beacuse not interesting
           #print ( "shift5", shift5, "bit 7, 8, 9", shift5[len(shift5)-3:])
           if int(shift5[len(shift5)-3:]) > 0:
           	#print (shift5[len(shift5)-3:])
           	channelStat = 3 #third type of status
           	cutBinStr = shift5[len(shift5)-3:]
           	if cutBinStr[2] == "1": #if I have OHVMax
           		StatusMeaning = "OHVMax"
           		extensibleStat = extensibleStat + StatusMeaning + " "
           		#print(StatusMeaning)
           	if cutBinStr[0] == "1": #if I have INTTRIP
           		StatusMeaning = "InTrip"
           		extensibleStat = extensibleStat + StatusMeaning + " "
           		#print(StatusMeaning)


	   #remove bit 7 8 9 to do the fourth status cathegory
           removedBits = removedBits - 3 #negative number
           shift6 = binStat[:removedBits]

           if len(shift6) != (nBit + 2) + removedBits:
           	mismatch6 = mismatch6 + 1

	   #for the fourth status cathegory I need the last bit of shift6
           #print ( "shift6", shift6, "bit 10", shift6[len(shift6)-1:])
           if int(shift6[len(shift6)-1:]) > 0:
           	#print (shift6[len(shift6)-1:])
           	channelStat = 4 #fourth type of status
           	cutBinStr = shift6[len(shift6)-1:]
                if cutBinStr[0] == "1": #if I have Calib Error
           		StatusMeaning = "CalibERR"
           		extensibleStat = extensibleStat + StatusMeaning + " "
           		#print(StatusMeaning)

           #remove bit 10
           removedBits = removedBits - 1 #negative number
           shift7 = binStat[:removedBits]

           if len(shift7) != (nBit + 2) + removedBits:
           	mismatch7 = mismatch7 + 1


           #for the fifth status cathegory I need the last bit of shift7
           #print ( "shift7", shift7, "bit 11", shift7[len(shift7)-1:])
           if int(shift7[len(shift7)-1:]) > 0:
           	#print (shift7[len(shift7)-1:])
           	channelStat = 5 #fifth type of status
           	cutBinStr = shift7[len(shift7)-1:]
          	if cutBinStr[0] == "1": #if I have Unplugged
           		StatusMeaning = "Unplugged"
           		extensibleStat = extensibleStat + StatusMeaning + " "
           		#print(StatusMeaning)

	   #remove bit 11
           removedBits = removedBits - 1 #negative number
           shift8 = binStat[:removedBits]

           if len(shift8) != (nBit + 2) + removedBits:
           	mismatch8 = mismatch8 + 1

           #print ( "shift8", shift8, "bit 12", shift8[len(shift8)-1:])   #bit 12 not interesting

	   #remove bit 12 to do the sixth status cathegory
           removedBits = removedBits - 1 #negative number
           shift9 = binStat[:removedBits]

           if len(shift9) != (nBit + 2) + removedBits:
           	mismatch9 = mismatch9 + 1

           #for the sixth status cathegory I need the last bit of shift9
           #print ( "shift9", shift9, "bit 13", shift9[len(shift9)-1:])
           if int(shift9[len(shift9)-1:]) > 0:
           	#print (shift9[len(shift9)-1:])
           	channelStat = 6 #sixth type of status
           	cutBinStr = shift9[len(shift9)-1:]
           	if cutBinStr[0] == "1": #if I have OVVPROT
           		StatusMeaning = "OVVPROT"
           		extensibleStat = extensibleStat + StatusMeaning + " "
           		#print(StatusMeaning)

	   #remove bit 13 to do the seventh status cathegory
           removedBits = removedBits - 1 #negative number
           shift10 = binStat[:removedBits]

           if len(shift10) != (nBit + 2) + removedBits:
           	mismatch10 = mismatch10 + 1

           #for the seventh status cathegory I need the last bit of shift10
           #print ( "shift10", shift10, "bit 14", shift10[len(shift10)-1:])
           if int(shift10[len(shift10)-1:]) > 0:
           	#print (shift10[len(shift10)-1:])
           	channelStat = 7 #seventh type of status
           	cutBinStr = shift10[len(shift10)-1:]
           	if cutBinStr[0] == "1": #if I have POWFAIL
           		StatusMeaning = "POWFAIL"
           		extensibleStat = extensibleStat + StatusMeaning + " "
           		#print(StatusMeaning)

	   #remove bit 14 to do the eight status cathegory
           removedBits = removedBits - 1 #negative number
           shift11 = binStat[:removedBits]

           if len(shift11) != (nBit + 2) + removedBits:
           	mismatch11 = mismatch11 + 1

           #for the eight status cathegory I need the last bit of shift11
           #print ( "shift11", shift11, "bit 15", shift11[len(shift11)-1:])
           if int(shift11[len(shift11)-1:]) > 0:
           	#print (shift11[len(shift11)-1:])
           	channelStat = 8 #eighth type of status
           	cutBinStr = shift11[len(shift11)-1:]
           	if cutBinStr[0] == "1": #if I have TEMPERR
           		StatusMeaning = "TEMPERR"
           		extensibleStat = extensibleStat + StatusMeaning + " "
           		#print(StatusMeaning)

	   #list S and T for the tgraph status vs time (millisecond int)
           smonOnlyT.append(tot_secondsSmon)
           smonOnlyS.append(channelStat)
           smonDecimalStatus.append( decimalStatus )
           smonOnlyTDateString.append(str(currentTsSmon))
           smonOnlyMeaningStat.append(extensibleStat)

           #th1 for status
           #SmonTh1List[counter].Fill(channelStat)
           SmonTh1List[counter].Fill(decimalStatus)

           SmonMeaning = ("StatChamber"+chamberList[indexB], "time:"+str(result[0]), "Status:"+StatusMeaning )
           SmonMeaningList.append(SmonMeaning)

           smonOnlyTDate.append(4.) #an any float number

           #print currentTsSmon
           year = str(currentTsSmon)[0:4]
           #print year
           month = str(currentTsSmon)[5:7]
           #print month
           day = str(currentTsSmon)[8:10]
           #print day
           hour = str(currentTsSmon)[11:13]
           #print hour
           minute = str(currentTsSmon)[14:16]
           #print minute
           second = str(currentTsSmon)[17:19]
           #print second
           micro = str(currentTsSmon)[20:26]
           #print micro

           #longString = ROOT.TString( year+"-"+month+"-"+day+" "+hour+":"+minute+":"+second )
           longList =str( year+"-"+month+"-"+day+" "+hour+":"+minute+":"+second )

           da1 = ROOT.TDatime( longList )

           smonOnlyTDate[-1] = da1.Convert()
           floatMicro = "0."+micro
           #print floatMicro
           smonOnlyTDate[contatoreS] = smonOnlyTDate[contatoreS] + float(floatMicro) #add microseconds to times (all ending with .0 because of Convert)


           contatoreS=contatoreS+1

        #th1 status
        SmonTh1List[counter].Write()

        #tgraph status
        if len(smonOnlyT) != len( smonOnlyS ):
        	print("!!!!!error tgraph status: filling with lists of different lenght!!!!!")
		file = open("LVErr.log", "w")
                file.write("ERROR: smonOnlyT and smonOnlyS have different lenght")
                file.close()

	if len(smonOnlyTDate) != len( smonOnlyS ):
        	print("!!!!!error tgraph status: filling with lists of different lenght!!!!!")
        	file = open("LVErr.log", "w")
                file.write("ERROR: smonOnlyTDate and smonOnlyS have different lenght")
                file.close()


        #sort the array of smonOnly and the smonOnly
        #in the case the query is not executed in order (negative times)
        #pair the time with status and the meaning list
        SortList = []
        for sortCount in range(len(smonOnlyT)):
        	internalList = (smonOnlyT[sortCount], smonOnlyS[sortCount], SmonMeaningList[sortCount], smonOnlyTDate[sortCount], smonOnlyTDateString[sortCount])
        	SortList.append(internalList)

        #print(SortList)
        SortList = sorted(SortList, key=lambda element: element[0])
        #print(SortList)

        for refill in range(len(smonOnlyT)):
        	smonOnlyT[refill]=SortList[refill][0]
        	smonOnlyS[refill]=SortList[refill][1]
        	SmonMeaningList[refill] = SortList[refill][2]
        	smonOnlyTDate[refill]=SortList[refill][3]
        	smonOnlyTDateString[refill]=SortList[refill][4]

        #print(smonOnlyT)
        #print(smonOnlyS)
        #print(SmonMeaningList)

        #rescale the negative times

	if (len(smonOnlyT))==0:
        	print("smonOnlyT lenght", len(smonOnlyT))
        	print("------------------------------------------------------------------")
        	print("ERROR: there are no status data for chamber "+ chamberList[indexB])
        	print("------------------------------------------------------------------")

        	file = open("LVErr.log", "w")
        	file.write("ERROR: there are no status data for chamber "+ chamberList[indexB])
        	file.close()

		#counter = counter + 1
		#I put a anomalous status value if there is no data
                smonOnlyT.append(0.)
                smonOnlyTDate.append(0.)
                smonOnlyS.append(-1000000000)
                smonDecimalStatus.append(-1000000000)
                smonOnlyTDateString.append("NO TS")
                smonOnlyMeaningStat.append("NO STATUS")
                smonOnlyBinStat.append("NOTHING")

        	#continue


        negativeStart = smonOnlyT[0]
        if smonOnlyT[0] < 0:
        	for iterTime in range(len(smonOnlyT)):
        		smonOnlyT[iterTime] = smonOnlyT[iterTime] + negativeStart*(-1)
        #print(smonOnlyT)

        #Smontg1 = ROOT.TGraph(n,x,y); x and y is the name of arrays with numbers of time and V
        #Smontg1 = ROOT.TGraph(len(smonOnlyT),smonOnlyTDate,smonOnlyS)
        Smontg1 = ROOT.TGraph(len(smonOnlyT),smonOnlyTDate,smonDecimalStatus)
        Smontg1.SetLineColor(2)
        Smontg1.SetLineWidth(4)
        Smontg1.SetMarkerColor(4)
        Smontg1.SetMarkerStyle(21)
        Smontg1.SetMarkerSize(1)
        Smontg1.SetName("LV_StatusChamber"+chamberList[indexB]+"_UTC_time")
        Smontg1.SetTitle("LV_StatusChamber"+chamberList[indexB]+"_UTC_time")
        #Smontg1.GetXaxis().SetTitle("time [s]")
        #Smontg1.GetYaxis().SetTitle("status cathegory "+chamberList[indexB])
        Smontg1.GetYaxis().SetTitle("status code "+chamberList[indexB])
        #Smontg1.Draw("ACP")
	Smontg1.GetXaxis().SetTimeDisplay(1)
        Smontg1.GetXaxis().SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00")
        Smontg1.GetXaxis().SetLabelOffset(0.025)

        #SmonTgraph1List.append(Smontg1)
        #SmonTgraph1List[counter].Write()
	Smontg1.Write()

	#List with all status
        SmonMeaningListList.append(SmonMeaningList)
        #print(SmonMeaningListList[counter])

        #tree for the status
        StatusTree = ROOT.TTree("LV_StatusTree"+chamberList[indexB], "LV_StatusTree"+chamberList[indexB])
        smonRootTimes = ROOT.vector('float')()
        smonRootTimesDate = ROOT.vector('string')()
        smonRootBinStat	= ROOT.vector('string')()
        smonRootDecimalStat = ROOT.vector('string')()
        smonRootMeaningStat = ROOT.vector('string')()

        #StatusTree.Branch( 'TS', smonRootTimes )
        StatusTree.Branch( 'TS', smonRootTimesDate )
        StatusTree.Branch( 'DecimalStat', smonRootDecimalStat )
        StatusTree.Branch( 'BinaryStat', smonRootBinStat )
        StatusTree.Branch( 'MeaningStat', smonRootMeaningStat )

        for lungh in range(len( smonOnlyT )):
        	smonRootTimes.push_back( smonOnlyT[lungh] )
        	smonRootBinStat.push_back( smonOnlyBinStat[lungh] )
        	smonRootDecimalStat.push_back( str(smonDecimalStatus[lungh]) )
        	smonRootMeaningStat.push_back( smonOnlyMeaningStat[lungh] )
        	smonRootTimesDate.push_back( smonOnlyTDateString[lungh] )

        StatusTree.Fill()

        StatusTree.Write()



	#FINAL COUNTER
	counter = counter + 1

f1.Close()

print('\n-------------------------Output--------------------------------')
print( fileName + " has been created.")
print("It is organised in directories: to change directory use DIRNAME->cd()")
print('To draw a TH1 or a TGraph: OBJNAME->Draw()')
print('To scan the root file use for example:\nLV_StatusTree2_2_Top->Scan("","","colsize=26")')
print("ALL MONITOR TIMES ARE IN UTC, DCS TIMES ARE IN CET")
