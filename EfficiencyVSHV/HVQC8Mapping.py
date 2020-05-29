#call the function HVEffPosQC8Stand and convert the position in the stand in the mapped HV channel
#this script gives in output the fisrt part of the string needed for the query of an element in the table ELEMENTS_ALL (it contains the HV values for each )
# the returned string is of the form cms_gem_dcs_1:CAEN/904_HV_mainframe/board00/channel007
import cx_Oracle
from datetime import datetime

import HVEffPosQC8Stand

def HVQC8Mapping ( chamberName, position, runNumber ):
	#print ( chamberName, position, runNumber )

	#take the run number and get the start and end date from the table CMS_GEM_MUON_VIEW.QC8_GEM_CH_VFAT_EFF_VIEW_RH
	db_cond = os.environ["GEM_PRODUCTION_DB_COND"]
    db_name = os.environ["GEM_PRODUCTION_DB_NAME"]

    db = cx_Oracle.connect(db_cond+db_name) # production DB
    cur = db.cursor()

	query = "select RUN_BEGIN, RUN_END from CMS_GEM_MUON_VIEW.QC8_GEM_CH_VFAT_EFF_VIEW_RH where CHAMBER_NAME ='"+str(chamberName)+"' and RUN_NUMBER='"+str(runNumber)+"'"

	#print query
	cur.execute(query)
	curStart = cur
	for result in curStart:
		startDate = result[0]
		endDate = result[1]

	print "runNumber: "+ str( runNumber ) +" | startDate: "+ str( startDate ) +" | endDate"+ str( endDate )



	#make the code understand which mapping use
	mappingChangeDate = []
        firstMappingChange 	= datetime( 2019, 04, 03, 13, 10, 00 )
        secondMappingChange = datetime( 2019, 04, 15,  8, 15, 00 )
        mappingChangeDate.append( firstMappingChange )
        mappingChangeDate.append( secondMappingChange )

        periodBool = []
        numberOfMaps = len(mappingChangeDate)+1 #date 10-May-2019

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

        print "periodBool", periodBool

	#validity limits needed for the query of ELEMENT_ID
        #find the first one and the last one in the bool vector
        firstOne = periodBool.index(1)
        lastOne  = max(loc for loc, val in enumerate(periodBool) if val == 1)
        #print "firstOne", firstOne
        #print "lastOne", lastOne


	position = position.replace("/","_")


	#all the mapping structure
	ChannelMapList = [ "Top_G3Bot", "Top_G3Top", "Top_G2Bot", "Top_G2Top", "Top_G1Bot", "Top_G1Top", "Top_Drift", "Bot_G3Bot", "Bot_G3Top", "Bot_G2Bot", "Bot_G2Top", "Bot_G1Bot", "Bot_G1Top", "Bot_Drift"]

        #MainframeMapList and BoardMapList lists follow the order of ChamberMapList
        ChamberMapList = [ "1_1", "1_2", "1_3", "2_1", "2_2", "2_3", "3_1", "3_2", "3_3", "4_1", "4_2", "4_3", "5_1", "5_2", "5_3" ]

        MainframeMapList1 = [ "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe", "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe" ]
        BoardMapList1 = [ "board00", "board01", "board02", "board03", "board04", "board05", "board06", "board07", "board08", "board09", "board10", "board11", "board12", "board13", "board14" ]

        MainframeMapList2 = [ "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe", "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe" ]
        BoardMapList2 = [ "board00", "board01", "board02", "board03", "board04", "board05", "board06", "board07", "board08", "board10", "board14", "board15", "board11", "board12", "board13" ]

        MainframeMapList3 = [ "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe",  "904_HV_mainframe", "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe",  "904_Shared_mainframe" ]
        BoardMapList3 = [ "board00", "board01", "board02", "board03", "board04", "board05", "board06", "board07", "board08", "board10", "board14", "board15", "board11", "board12", "board13" ]

        channelList = [ "G3Bot", "G3Top", "G2Bot", "G2Top", "G1Bot", "G1Top", "Drift"]

        listAllMainframeMappings = []
        listAllMainframeMappings.append( MainframeMapList1 )
        listAllMainframeMappings.append( MainframeMapList2 )
        listAllMainframeMappings.append( MainframeMapList3 )

        listAllBoardMappings = []
        listAllBoardMappings.append( BoardMapList1 )
        listAllBoardMappings.append( BoardMapList2 )
        listAllBoardMappings.append( BoardMapList3 )


	position = position.replace("/","_")

	#print position[:3]
	chamberIdx = ChamberMapList.index(str(position[:3]))
	#print chamberIdx
	#print ChamberMapList[ chamberIdx ]

	#prepare the string for the query in table GEM_904_COND.ELEMENTS_ALL
	#create a query string for the maps that must me used
        #usedMaps tells the number of used maps for the period between start and end, not the number of maps charged in the code (said instead by numberOfMaps)
        usedMaps = periodBool.count(1)
        vmonNameList 	= []
        contMaps = 0
	for contMaps in range( usedMaps ):
        	vmon_name="'cms_gem_dcs_1:CAEN/"+listAllMainframeMappings[startIdx + contMaps][chamberIdx]+"/"+listAllBoardMappings[startIdx + contMaps][chamberIdx]+"/"
		if position[4] == "T":
        		vmon_name = vmon_name + "channel000"
        	elif position[4] == "B":
        		vmon_name = vmon_name + "channel007"
        	#this string has to be saved only one time
        	if vmonNameList.count( vmon_name ) == 1:
        		continue
        	vmonNameList.append(vmon_name)

	#print vmonNameList


	return vmonNameList



#TEST CODE

#chNAME = "GE1/1-VII-L-CERN-0001"
#runNUMBER = 12
#print HVEffPosQC8Stand.HVEffPosQC8Stand ( chNAME, runNUMBER )
#pos = HVEffPosQC8Stand.HVEffPosQC8Stand ( chNAME, runNUMBER )


#vmonNAMELIST=HVQC8Mapping( chNAME, pos, runNUMBER)
#print "vmonNAMELIST", vmonNAMELIST

#print HVQC8Mapping(chNAME, pos, 1)
