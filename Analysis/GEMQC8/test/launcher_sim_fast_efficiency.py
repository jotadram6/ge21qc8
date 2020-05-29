import csv
import os
import sys
import io
import subprocess
import time
import datetime

if __name__ == '__main__':

    # Define the parser
    import argparse
    parser = argparse.ArgumentParser(description="QC8 fast efficiency simulation. For any doubt: https://twiki.cern.ch/twiki/bin/view/CMS/GEMCosmicRayAnalysis")
    # Positional arguments
    parser.add_argument("run_number", type=int, help="Specify the run number")
    args = parser.parse_args()

    # Different paths definition
    srcPath = os.path.abspath("launcher_sim_fast_efficiency.py").split('QC8Test')[0]+'QC8Test/src/'
    pyhtonModulesPath = os.path.abspath("launcher_sim_fast_efficiency.py").split('QC8Test')[0]+'QC8Test/src/Analysis/GEMQC8/python/'
    runPath = os.path.abspath("launcher_sim_fast_efficiency.py").split('QC8Test')[0] + 'QC8Test/src/Analysis/GEMQC8/test/'
    configTablesPath = os.path.abspath("launcher_sim_fast_efficiency.py").split('QC8Test')[0] + 'QC8Test/src/Analysis/GEMQC8/data/StandConfigurationTables/'
    outDirPath = os.path.abspath("launcher_sim_fast_efficiency.py").split('QC8Test')[0] + "Results_QC8_sim_fast_efficiency_run_"+str(args.run_number)

    sys.path.insert(0,pyhtonModulesPath)

    import config_creator
    import geometry_files_creator

    # Generate configuration file
    config_creator.configMaker(args.run_number,"noMasks")
    time.sleep(1)

    # Generate geometry files
    geometry_files_creator.geomMaker(args.run_number,"noAlignment")
    time.sleep(1)

    # Compiling after the generation of the geometry files
    scramCommand = "scram build -j 4"
    scramming = subprocess.Popen(scramCommand.split(),stdout=subprocess.PIPE,universal_newlines=True,cwd=srcPath)
    while scramming.poll() is None:
        line = scramming.stdout.readline()
        print(line)
    print scramming.stdout.read()
    scramming.communicate()
    time.sleep(1)

    # Running the CMSSW code
    runCommand = "cmsRun -n 8 runGEMCosmicStand_sim_fast_efficiency.py"
    running = subprocess.Popen(runCommand.split(),stdout=subprocess.PIPE,universal_newlines=True,cwd=runPath)
    while running.poll() is None:
        line = running.stdout.readline()
        print(line)
    print running.stdout.read()
    running.communicate()
    time.sleep(1)

    # Creating folder outside the CMMSW release to put the output files and plots
    #---# Remove old version if want to recreate
    if (os.path.exists(outDirPath)):
        rmDirCommand = "rm -rf "+outDirPath
        rmDir = subprocess.Popen(rmDirCommand.split(),stdout=subprocess.PIPE,universal_newlines=True)
        rmDir.communicate()
    #---# Create the new empty folder
    resDirCommand = "mkdir "+outDirPath
    resDir = subprocess.Popen(resDirCommand.split(),stdout=subprocess.PIPE,universal_newlines=True)
    resDir.communicate()
    time.sleep(1)

    # Create folders for ouput plots per chamber
    import configureRun_cfi as runConfig
    SuperChType = runConfig.StandConfiguration
    ChID = runConfig.ChamberIDs
    for i in range (0,30):
        if (SuperChType[int(i/2)] != '0'):
            plotsDirCommand = "mkdir outPlots_Chamber_Pos_" + str(i) + "_" + ChID[i]
            plotsDirChamber = subprocess.Popen(plotsDirCommand.split(),stdout=subprocess.PIPE,universal_newlines=True,cwd=outDirPath)
            plotsDirChamber.communicate()
    time.sleep(1)

    # Selecting the correct output file, changing the name and moving to the output folder
    out_name = 'out_run_'
    for i in range(6-len(str(args.run_number))):
        out_name = out_name + '0'
    out_name = out_name + str(args.run_number) + '.root'

    mvToDirCommand = "mv fast_efficiency_" + out_name + " " + outDirPath + "/fast_efficiency_" + out_name
    movingToDir = subprocess.Popen(mvToDirCommand.split(),stdout=subprocess.PIPE,universal_newlines=True,cwd=runPath)
    movingToDir.communicate()
    time.sleep(1)

    # Efficiency computation & output
    rootMacroCommand = "root -l -q -b " + runPath + "macro_fast_efficiency.c(" + str(args.run_number) + ",\"" + configTablesPath + "\")"
    rootMacroProcess = subprocess.Popen(rootMacroCommand.split(),stdout=subprocess.PIPE,universal_newlines=True,cwd=outDirPath)
    while rootMacroProcess.poll() is None:
        line = rootMacroProcess.stdout.readline()
        print(line)
    print rootMacroProcess.stdout.read()
    rootMacroProcess.communicate()
    time.sleep(1)
