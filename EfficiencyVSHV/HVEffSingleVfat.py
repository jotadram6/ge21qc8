#!/usr/bin/python
import cx_Oracle

#input: chamber name (for example GE1/1-VII-L-CERN-0002), run number and vfat number (0-23)
def HVEffSingleVfat( chamberName, runNumber, vfatNumber ):
	#print ( chamberName, runNumber )

	#connect to the DB to find position of the chamber, and run number  in the table
	db_cond = os.environ["GEM_PRODUCTION_DB_COND"]
    db_name = os.environ["GEM_PRODUCTION_DB_NAME"]

    db = cx_Oracle.connect(db_cond+db_name) # production DB
    cur = db.cursor()

	#print "chamberName", chamberName, "runNumber", runNumber, "vfatNumber", vfatNumber

	#vfatNumber format: 0 1 2 ... 23
	query = "select EFFICIENCY, EFFICIENCY_ERROR from CMS_GEM_MUON_VIEW.QC8_GEM_CH_VFAT_EFF_VIEW_RH where CHAMBER_NAME='"+chamberName+"' and RUN_NUMBER="+str(runNumber)+ " and VFAT_POSN="+str(vfatNumber)
        cur.execute(query)
	curEff = cur
	for result in curEff:
		efficiency 	= result[0]
		effError	= result[1]

	effList = [efficiency, effError]
        #print "CHAMBER_NAME: ", chamber_name, "POSITION: ", position, "RUN_NUMBER: ", run_number

	return effList


#TEST CODE
#chNAME = "GE1/1-VII-L-CERN-0001"
#runNUMBER = 12
#vfatNUMBER = 4
#howMANYVFATS = 24
#HVEffSingleVfat( chNAME, runNUMBER, vfatNUMBER )
#for vfatNUMBER in range( howMANYVFATS ):
#	effLIST = HVEffSingleVfat( chNAME, runNUMBER, vfatNUMBER )
#	print effLIST
