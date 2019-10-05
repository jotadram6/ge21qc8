#!/usr/bin/python
import cx_Oracle
#import ROOT

#input: chamber name (for example GE1/1-VII-L-CERN-0002), run number and vfat number (0-23)
def HVEffPosQC8Stand( chamberName, runNumber ):
	#print ( chamberName, runNumber )

	#connect to the DB to find position of the chamber, and run number  in the table
	db_cond = os.environ["GEM_PRODUCTION_DB_COND"]
    db_name = os.environ["GEM_PRODUCTION_DB_NAME"]

    db = cx_Oracle.connect(db_cond+db_name) # production DB
    cur = db.cursor()

	chamberName = "'"+chamberName+"'"
	query = "select CH_SERIAL_NUMBER, POSITION, RUN_NUMBER from CMS_GEM_MUON_VIEW.QC8_GEM_STAND_GEOMETRY_VIEW_RH where CH_SERIAL_NUMBER="+chamberName+" and RUN_NUMBER="+str(runNumber)
        cur.execute(query)
	curGeom = cur
	for result in curGeom:
		chamber_name 	= result[0]
		position	= result[1]
		run_number	= result[2]

        #print "CHAMBER_NAME: ", chamber_name, "POSITION: ", position, "RUN_NUMBER: ", run_number

	return position


#call the function
#print HVEffPosQC8Stand("GE1/1-X-L-CERN-0001",1)
