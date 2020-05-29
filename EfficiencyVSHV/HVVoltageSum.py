#this script extracts HV for each voltage channel of the chamber,
#makes a medium of voltage for each channel on the registered time stamps
#and finally sums the seven channels
import HVQC8Mapping
import HVEffPosQC8Stand
import cx_Oracle
import math
from datetime import datetime, timedelta

def HVVoltageSum ( firstPartQueryStringList, runNumber, chamberName ):
	#print "How many strings for query in input:", len( firstPartQueryStringList )

	if len( firstPartQueryStringList ) != 1:
		print "ERROR: len( firstPartQueryStringList )="+str(len( firstPartQueryStringList ))
		print "It should be 1"

	#take only the first string in the list: I shuold have only one string per each runNumber
	firstPartQueryString = firstPartQueryStringList[0]

	#print "first channel called in the HV query: " + firstPartQueryString

	#retrieve the begin and end time of a run_number from the table
	#QC8_GEM_CH_VFAT_EFF_VIEW_RH
	db_cond = os.environ["GEM_PRODUCTION_DB_COND"]
    db_name = os.environ["GEM_PRODUCTION_DB_NAME"]

    db = cx_Oracle.connect(db_cond+db_name) # production DB
    cur = db.cursor()

	query = "select RUN_NUMBER, RUN_BEGIN, RUN_END from CMS_GEM_MUON_VIEW.QC8_GEM_CH_VFAT_EFF_VIEW_RH where RUN_NUMBER="+str(runNumber)+"and CHAMBER_NAME='"+str(chamberName)+"'"
	cur.execute(query)
	curRunDates = cur
	for result in curRunDates:
		run_number = result[0]
		run_begin  = result[1]
		run_end    = result[2]

	#print query
	#print ("run_number:", run_number, "run_begin:", run_begin, "run_end:", run_end )
	#these dates are in the format YYYY-MM-DD HH24:mm:ss, good for the query

	#print run_begin, run_end



	#retrieve HV data on each HV channel (7 channels in a chamber)
	voltageMediumList = []
	for channelIdx in range( 7 ):
		#right format 00 for the HV channel number in the query string
		#print firstPartQueryString, "lastTwoChar", firstPartQueryString[-2:]
		lastTwoChar = int(firstPartQueryString[-2:]) + channelIdx
		if lastTwoChar < 10:
			lastTwoChar = "0"+str(lastTwoChar)
		#print lastTwoChar

		#retrieve HV data for a single HV channel ( HV data are in table ELEMENTS_ALL )
		vmon_name = firstPartQueryString[:-2]+ str(lastTwoChar) + ".actual.vMon'"

		#check if VALID_TILL is NULL
		query = "select ELEMENT_ID, VALID_SINCE, VALID_TILL from ELEMENTS_ALL where ELEMENT_NAME="+vmon_name+ " and VALID_SINCE < to_date ('"+str(run_begin)+"','YYYY-MM-DD HH24:MI:SS') and ( (VALID_TILL > to_date ('"+str(run_end)+"','YYYY-MM-DD HH24:MI:SS') and VALID_TILL IS NOT NULL) or ( VALID_TILL IS NULL ) )"
		#print query
		cur.execute(query)
		curSince = cur
		curSinceLen = 0
		for result in curSince:
			vmon_id = result[0]
			curSinceLen = curSinceLen + 1

		if curSinceLen >1:
			print "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
			print "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
			print "ERROR: more than one HV ELEMENT_ID"
			print "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
			print "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n"

        	#vmon_id = cur.fetchone()[0];

		#print "VMON_ID ", vmon_id

        	query = "select TS,VALUE_NUMBER from EVENTHISTORY where ELEMENT_ID = "+str(vmon_id) +" and TS > to_date ('"+str(run_begin)+"','YYYY-MM-DD HH24:MI:SS') and TS < to_date ('"+str(run_end)+"','YYYY-MM-DD HH24:MI:SS')"
		cur.execute(query)

		voltageSum = 0
		voltageCounter = 0

		curvmon = cur
		for result in curvmon:
			ts = result[0]
			voltage =  result[1]
			print voltage
			voltageSum = voltageSum + voltage
			voltageCounter = voltageCounter + 1

		#print "voltageCounter", voltageCounter
		#print "voltageSum", voltageSum

		#print type(run_begin)

		newBegin = run_begin
		lastBool = False

		#go before the start to look for the last HV refistered (it is the present HV value)
		while voltageCounter == 0:
			lastBool = True
			newEnd   = newBegin
			newBegin = newBegin - timedelta( days = 1 )

			#check that the newBegin does not become  less than the initial validity of the ID for voltage
			query = "select ELEMENT_ID, VALID_SINCE from ELEMENTS_ALL where ELEMENT_NAME="+vmon_name + "and ELEMENT_ID="+str(vmon_id)
			cur.execute(query)

			line = cur.fetchone()
			validSince = line[1]

			if newBegin < validSince:
				newBegin = validSince
				#it would be forced to find something registered at the beginning of the ID
				#this should be sufficient to force it exit the while loop

			query = "select TS,VALUE_NUMBER from EVENTHISTORY where ELEMENT_ID = "+str(vmon_id) +" and TS > to_date ('"+str(newBegin)+"','YYYY-MM-DD HH24:MI:SS') and TS < to_date ('"+str(newEnd)+"','YYYY-MM-DD HH24:MI:SS')"
			cur.execute(query)

			curvmon = cur
        	        for result in curvmon:
                		ts = result[0]
                		voltage =  result[1]
	                	voltageCounter = voltageCounter + 1

		#I register only the last HV scanned
		voltageLast = voltage
		#print "voltageSum", voltageSum, "voltageLast", voltageLast, "voltageCounter", voltageCounter #the only thing important is that it is different from 0

		#average voltage on this channel
		if lastBool == False:
			voltageMedium = voltageSum/voltageCounter
			#calculate standard deviation for the voltage error on the single channel
			squares = 0
			for result in curvmon:
				voltage = result[1]
				squares = squares + ( voltage - voltageMedium )**2

			devStdData = math.sqrt( squares/voltageCounter )
			#error on media of data is devStdData/sqrt(nData)
			devStdMediaHVChannel = devStdData/(math.sqrt( voltageCounter ))
			errHVChannel = devStdMediaHVChannel

		if lastBool == True:
			voltageMedium = voltageLast

			#error on HV from voltage monitor vs Output accuracy (CAEN A1515 datasheet)
			#0.2% +/- 0.2 V +/- 50 ppm/C
			errHVChannel = 0.002*voltageMedium + 0.2

		#print "lastBool", lastBool, "voltageMedium", voltageMedium
		voltageMediumList.append( voltageMedium )

	#print "voltageMediumList", voltageMediumList
	#check number of channels read
	if len( voltageMediumList ) != 7:
		print "ERROR: expected 7 channels in voltageMediumList, "+str(len( voltageMediumList )) + "provided"

	voltageTotChamber = 0
	voltageTotChamberErr = 0
	for channelIdx in range(len(voltageMediumList)):
		voltageTotChamber = voltageTotChamber + voltageMediumList[ channelIdx ]
		#totalHVerror
		voltageTotChamberErr = voltageTotChamberErr + ( errHVChannel )**2

	#sqrt of the sum of squares
	voltageTotChamberErr = math.sqrt( voltageTotChamberErr )

	#print voltageTotChamber
	#print voltageTotChamberErr
	voltageTotChamberList = [ voltageTotChamber, voltageTotChamberErr ]

	return voltageTotChamberList

#TEST CODE
#chamberNameTest = "GE1/1-VII-L-CERN-0001"
#runNumberTest = 12
#chamberNameTest = "GE1/1-VII-L-CERN-0002"
#chamberNameTest = "GE1/1-VII-L-CERN-0002"
#runNumberTest = 12
#print HVEffPosQC8Stand.HVEffPosQC8Stand ( "GE1/1-X-L-CERN-0001",1 )
#pos = HVEffPosQC8Stand.HVEffPosQC8Stand ( chamberNameTest, runNumberTest )
#queryImportStrList = HVQC8Mapping.HVQC8Mapping( chamberNameTest, pos, runNumberTest )


#print queryImportStr

#totVOLTCHAMBER = HVVoltageSum ( queryImportStrList, runNumberTest, chamberNameTest )
#print "totVOLTCHAMBER", totVOLTCHAMBER
