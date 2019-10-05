import cx_Oracle
import os, sys, io
import datetime
import time

def getConfigurationTable(run_num,runDateTime):
    print "\nDownloading ConfigurationTable for run {0}".format(run_num)

    db_cond = os.environ["GEM_PRODUCTION_DB_COND"]
    db_name = os.environ["GEM_PRODUCTION_DB_NAME"]

    db = cx_Oracle.connect(db_cond+db_name) # production DB
    cur = db.cursor()

    delta_min = 94608000 # minimum value of delta (this value is the # of seconds in 3 years)
    goodDateTimeForQuery = ""

    dateRun = runDateTime.split('_')[0]
    yearRun = dateRun.split('-')[0]
    monthRun = dateRun.split('-')[1]
    dayRun = dateRun.split('-')[2]
    timeRun = runDateTime.split('_')[1]
    hourRun = timeRun.split('-')[0]
    minutesRun = timeRun.split('-')[1]
    secondsRun = timeRun.split('-')[2]

    query = "select TIME from CMS_GEM_MUON_VIEW.QC8_GEM_STAND_GEOMETRY_VIEW_RH"
    cur.execute(query)
    for result in cur:
        deltaT = (datetime.datetime(int(yearRun),int(monthRun),int(dayRun),int(hourRun),int(minutesRun),int(secondsRun)) - result[0]).total_seconds()
        if (deltaT >= 0 and deltaT < delta_min):
            delta_min = deltaT
            goodDateTimeForQuery = result[0]

    query = "select * from CMS_GEM_MUON_VIEW.QC8_GEM_STAND_GEOMETRY_VIEW_RH"
    cur.execute(query)

    allChambers = []

    for result in cur:
        chamber_name = result[0]
        gem_num      = result[1]
        position     = result[2]
        ch_type      = result[3]
        flip         = result[4]
        amc          = result[5]
        oh           = result[6]
        flow_meter   = result[7]
        time         = result[8]
        run_number   = run_num
        if (time == goodDateTimeForQuery):
            line = str(chamber_name) + "," + str(gem_num) + "," + str(position) + "," + str(ch_type) + "," + str(flip) + "," + str(amc) + "," + str(oh) + "," + str(flow_meter) + "," + str(run_number) + "\n"
            allChambers.append(line)

    orderedPositions = ["1/1/B","1/1/T","2/1/B","2/1/T","3/1/B","3/1/T","4/1/B","4/1/T","5/1/B","5/1/T",\
                        "1/2/B","1/2/T","2/2/B","2/2/T","3/2/B","3/2/T","4/2/B","4/2/T","5/2/B","5/2/T",\
                        "1/3/B","1/3/T","2/3/B","2/3/T","3/3/B","3/3/T","4/3/B","4/3/T","5/3/B","5/3/T"]

    configTablesPath = os.path.abspath("dumpDBtables.py").split('QC8Test')[0] + 'QC8Test/src/Analysis/GEMQC8/data/StandConfigurationTables/'
    outfile_name = configTablesPath + "StandGeometryConfiguration_run{0}.csv".format(run_num)

    with open(outfile_name,"w+") as outfile:
        line = "CH_SERIAL_NUMBER,GEM_NUM,POSITION,CH_TYPE,FLIP,AMC,OH,FLOW_METER,RUN_NUMBER\n"
        print line
        outfile.write(line)
        for pos in range(30):
            for ch in range(len(allChambers)):
                if (allChambers[ch].split(',')[2] == orderedPositions[pos]):
                    print(allChambers[ch])
                    outfile.write(allChambers[ch])

    print "\nSuccesfully done!\n"

def getAlignmentTable(run_num):
    print "\nDownloading AlignmentTable for run {0}".format(run_num)

    db_cond = os.environ["GEM_PRODUCTION_DB_COND"]
    db_name = os.environ["GEM_PRODUCTION_DB_NAME"]

    db = cx_Oracle.connect(db_cond+db_name) # production DB
    cur = db.cursor()

    query = "select * from CMS_GEM_MUON_VIEW.QC8_GEM_ALIGNMENT_VIEW_RH where RUN_NUMBER="+str(run_num)
    cur.execute(query)

    allChambers = []

    for result in cur:
        position   = result[0]
        dx         = result[1]
        dy         = result[2]
        dz         = result[3]
        rx         = result[4]
        ry         = result[5]
        rz         = result[6]
        run_number = result[7]
        line = str(position) + "," + str(dx) + "," + str(dy) + "," + str(dz) + "," + str(rx) + "," + str(ry) + "," + str(rz) + "," + str(run_number) + "\n"
        allChambers.append(line)

    orderedPositions = ["1/1/B","1/1/T","2/1/B","2/1/T","3/1/B","3/1/T","4/1/B","4/1/T","5/1/B","5/1/T",\
                        "1/2/B","1/2/T","2/2/B","2/2/T","3/2/B","3/2/T","4/2/B","4/2/T","5/2/B","5/2/T",\
                        "1/3/B","1/3/T","2/3/B","2/3/T","3/3/B","3/3/T","4/3/B","4/3/T","5/3/B","5/3/T"]

    alignTablesPath = os.path.abspath("dumpDBtables.py").split('QC8Test')[0] + 'QC8Test/src/Analysis/GEMQC8/data/StandAligmentTables/'
    outfile_name = alignTablesPath + "StandAlignmentValues_run{0}.csv".format(run_num)

    with open(outfile_name,"w+") as outfile:
        line = "POSITION,DX_cm,DY_cm,DZ_cm,RX_deg,RY_deg,RZ_deg,RUN_NUMBER\n"
        outfile.write(line)
        for pos in range(30):
            for ch in range(len(allChambers)):
                if (allChambers[ch].split(',')[0] == orderedPositions[pos]):
                    print(allChambers[ch])
                    outfile.write(allChambers[ch])

    print "\nSuccesfully done!\n"

def getHotStripsTable(run_num):
    print "\nDownloading HotStripsTable for run {0}".format(run_num)

    db_cond = os.environ["GEM_PRODUCTION_DB_COND"]
    db_name = os.environ["GEM_PRODUCTION_DB_NAME"]

    db = cx_Oracle.connect(db_cond+db_name) # production DB
    cur = db.cursor()

    query = "select * from CMS_GEM_MUON_VIEW.QC8_GEM_MASKED_STRIPS_HOT_V_RH where RUN_NUMBER="+str(run_num)
    cur.execute(query)

    allChambers = []

    for result in cur:
        chamber_name = result[0]
        gem_num      = result[1]
        position     = result[2]
        vfat         = result[3]
        channel      = result[4]
        strip        = result[5]
        run_number   = result[6]
        line = str(chamber_name) + "," + str(gem_num) + "," + str(position) + "," + str(vfat) + "," + str(channel) + "," + str(strip) + "," + str(run_number) + "\n"
        allChambers.append(line)

    orderedPositions = ["1/1/B","1/1/T","2/1/B","2/1/T","3/1/B","3/1/T","4/1/B","4/1/T","5/1/B","5/1/T",\
                        "1/2/B","1/2/T","2/2/B","2/2/T","3/2/B","3/2/T","4/2/B","4/2/T","5/2/B","5/2/T",\
                        "1/3/B","1/3/T","2/3/B","2/3/T","3/3/B","3/3/T","4/3/B","4/3/T","5/3/B","5/3/T"]

    hotStripsTablesPath = os.path.abspath("dumpDBtables.py").split('QC8Test')[0] + 'QC8Test/src/Analysis/GEMQC8/data/HotStripsTables/'
    outfile_name = hotStripsTablesPath + "HotStrips_run{0}.csv".format(run_num)

    with open(outfile_name,"w+") as outfile:
        line = "CH_SERIAL_NUMBER,GEM_NUMBER,POSITION,VFAT,CHANNEL,STRIP,RUN_NUMBER\n"
        outfile.write(line)
        for pos in range(30):
            for ch in range(len(allChambers)):
                if (allChambers[ch].split(',')[2] == orderedPositions[pos]):
                    print(allChambers[ch])
                    outfile.write(allChambers[ch])

    print "\nSuccesfully done!\n"

def getDeadStripsTable(run_num):
    print "\nDownloading DeadStripsTable for run {0}".format(run_num)

    db_cond = os.environ["GEM_PRODUCTION_DB_COND"]
    db_name = os.environ["GEM_PRODUCTION_DB_NAME"]

    db = cx_Oracle.connect(db_cond+db_name) # production DB
    cur = db.cursor()

    query = "select * from CMS_GEM_MUON_VIEW.QC8_GEM_MASKED_STRIPS_DEAD_RH where RUN_NUMBER="+str(run_num)
    cur.execute(query)

    deadStripsTablesPath = os.path.abspath("dumpDBtables.py").split('QC8Test')[0] + 'QC8Test/src/Analysis/GEMQC8/data/DeadStripsTables/'
    outfile_name = deadStripsTablesPath + "DeadStrips_run{0}.csv".format(run_num)

    allChambers = []

    for result in cur:
        chamber_name = result[0]
        gem_num      = result[1]
        position     = result[2]
        vfat         = result[3]
        channel      = result[4]
        strip        = result[5]
        run_number   = result[6]
        line = str(chamber_name) + "," + str(gem_num) + "," + str(position) + "," + str(vfat) + "," + str(channel) + "," + str(strip) + "," + str(run_number) + "\n"
        allChambers.append(line)

    orderedPositions = ["1/1/B","1/1/T","2/1/B","2/1/T","3/1/B","3/1/T","4/1/B","4/1/T","5/1/B","5/1/T",\
                        "1/2/B","1/2/T","2/2/B","2/2/T","3/2/B","3/2/T","4/2/B","4/2/T","5/2/B","5/2/T",\
                        "1/3/B","1/3/T","2/3/B","2/3/T","3/3/B","3/3/T","4/3/B","4/3/T","5/3/B","5/3/T"]

    with open(outfile_name,"w+") as outfile:
        line = "CH_SERIAL_NUMBER,GEM_NUMBER,POSITION,VFAT,CHANNEL,STRIP,RUN_NUMBER\n"
        outfile.write(line)
        for pos in range(30):
            for ch in range(len(allChambers)):
                if (allChambers[ch].split(',')[2] == orderedPositions[pos]):
                    print(allChambers[ch])
                    outfile.write(allChambers[ch])

    print "\nSuccesfully done!\n"

if __name__ == '__main__':

    # Define the parser
    import argparse
    parser = argparse.ArgumentParser(description="This script dumps the tables from the GEM production DB")
    # Positional arguments
    parser.add_argument("tableType", type=str, choices=["ConfigurationTable","AlignmentTable","HotStripsTable","DeadStripsTable"], help="Specify the table type")
    parser.add_argument("runNumber", type=int, help="Specify the run number")
    parser.add_argument("dateTimeOfRun", type=str, help="Specify the date and time of the run as YYYY-MM-DD_hh-mm-ss. Example: 2019-03-14_09-04-00")
    args = parser.parse_args()

    if args.tableType == "ConfigurationTable":
        getConfigurationTable(args.runNumber,args.dateTimeOfRun)
    elif args.tableType == "AlignmentTable":
        getAlignmentTable(args.runNumber)
    elif args.tableType == "HotStripsTable":
        getHotStripsTable(args.runNumber)
    elif args.tableType == "DeadStripsTable":
        getDeadStripsTable(args.runNumber)
