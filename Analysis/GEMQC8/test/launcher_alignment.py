import csv
import os
import sys
import io
import subprocess
import time
import datetime
import socket

def cmsRunner(split):
  runCommand = "cmsRun -n "+str(split)+" runGEMCosmicStand_alignment.py"
  running = subprocess.Popen(runCommand.split(),stdout=subprocess.PIPE,universal_newlines=True,cwd=runPath)
  while running.poll() is None:
    line = running.stdout.readline()
    print(line)
  print running.stdout.read()
  running.communicate()
  time.sleep(1)

def align_stopper(run_number, step):
  runPath = os.path.abspath("geometry_files_creator.py").split('QC8Test')[0] + 'QC8Test/src/Analysis/GEMQC8/test/'
  sys.path.insert(0,runPath)
  dx = []
  rz = []
  alignmentTablesPath = os.path.abspath("align_stopper.py").split('QC8Test')[0] + 'QC8Test/src/Analysis/GEMQC8/data/StandAligmentTables/'
  infileName = alignmentTablesPath + "StandAlignmentValues_run" + str(run_number) + "_step" + str(step) + ".csv"
  stop = 0
  stop_align = False
  if (os.path.exists(infileName)):
    with open(infileName) as infile:
      for line in infile:
        if ('\n' in line):
          line = line.split('\n')[0]
        if ('\r' in line):
          line = line.split('\r')[0]
        if (line.split(',')[0]=='RunNumber' and line.split(',')[1]!=run_number):
          sys.exit('StandAlignmentValues file has something wrong: run rumber not matching...')
        if (line.split(',')[0]!='RunNumber' and line.split(',')[0]!='Position'):
          dx.append(float(line.split(',')[1]))
          rz.append(float(line.split(',')[6]))

  for dx_ in dx:
    if(dx_ > 0.03):
      stop +=1
  for rz_ in rz:
    if(rz_ > 0.03):
      stop +=1
  if(stop < 1):
    stop_align = True
  return stop_align


if __name__ == '__main__':

  # Define the parser
  import argparse
  parser = argparse.ArgumentParser(description="QC8 data analysis step 4. Software alignment of the chambers in the stand. For any doubt: https://twiki.cern.ch/twiki/bin/view/CMS/GEMCosmicRayAnalysis")
  # Positional arguments
  parser.add_argument("run_number", type=int, help="Specify the run number")
  args = parser.parse_args()

  # Different paths definition
  srcPath = os.path.abspath("launcher_alignment.py").split('QC8Test')[0]+'QC8Test/src/'
  pyhtonModulesPath = os.path.abspath("launcher_alignment.py").split('QC8Test')[0]+'QC8Test/src/Analysis/GEMQC8/python/'
  runPath = os.path.abspath("launcher_alignment.py").split('QC8Test')[0] + 'QC8Test/src/Analysis/GEMQC8/test/'
  configTablesPath = os.path.abspath("launcher_alignment.py").split('QC8Test')[0] + 'QC8Test/src/Analysis/GEMQC8/data/StandConfigurationTables/'
  alignmentTablesPath = os.path.abspath("launcher_alignment.py").split('QC8Test')[0] + 'QC8Test/src/Analysis/GEMQC8/data/StandAligmentTables/'
  resDirPath = os.path.abspath("launcher_alignment.py").split('QC8Test')[0]

  sys.path.insert(0,pyhtonModulesPath)

  import dumpDBtables
  import config_creator
  import geometry_files_creator
  import convertAlignmentTables

  # Retrieve start date and time of the run

  if (socket.gethostname()=="gem904qc8dqm"):

        fpath =  "/data/bigdisk/GEM-Data-Taking/GE11_QC8/Cosmics/run{:06d}/".format(int(args.run_number))
        for x in os.listdir(fpath):
            if x.endswith("ls0001_index000000.raw"):
                file0name = x
                break
            elif x.endswith("chunk_000000.dat"):
                file0name = x
                break
        if (not file0name.endswith("ls0001_index000000.raw") and not file0name.endswith("chunk_000000.dat")):
            print "Check the data files... First file (at least) is missing!"
        t = datetime.datetime.fromtimestamp(os.path.getmtime(fpath+file0name))
        startDateTime = str("{:04d}".format(int(t.year))) + "-" + str("{:02d}".format(int(t.month))) + "-" + str("{:02d}".format(int(t.day))) + "_" + str("{:02d}".format(int(t.hour))) + "-" + str("{:02d}".format(int(t.minute))) + "-" + str("{:02d}".format(int(t.second)))
        time.sleep(1)

  else:

      fpath =  "/eos/cms/store/group/dpg_gem/comm_gem/QC8_Commissioning/run{:06d}/".format(int(args.run_number))
      for x in os.listdir(fpath):
          if x.endswith("ls0001_allindex.raw"):
              file0name = x
              break
          elif x.endswith("chunk_000000.dat"):
              file0name = x
              break
      if (not file0name.endswith("ls0001_index000000.raw") and not file0name.endswith("chunk_000000.dat")):
          print "Check the data files... First file (at least) is missing!"
      startDateTime = file0name.split('_')[3] + "_" + file0name.split('_')[4]
      time.sleep(1)

  # Get stand configuration table from the DB
  if int(args.run_number) > 237:
      dumpDBtables.getConfigurationTable(args.run_number,startDateTime)

  # Generate configuration file
  config_creator.configMaker(args.run_number)
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

  stop_align = False
  step = 0
  docheck = False
  import configureRun_cfi as runConfig
  cores = 8

  # Generate geometry files
  geometry_files_creator.geomMaker(args.run_number, "noAlignment")
  time.sleep(1)

  while not(stop_align or step>5):
    cmsRunner(cores)
    #  # Creating folder outside the CMMSW release to put the output files and plots
    outDirName = "Results_QC8_alignment_run_"+str(args.run_number)
    #---# Remove old version if want to recreate
    if (os.path.exists(resDirPath+outDirName)):
      rmDirCommand = "rm -rf "+outDirName
      rmDir = subprocess.Popen(rmDirCommand.split(),stdout=subprocess.PIPE,universal_newlines=True,cwd=resDirPath)
      rmDir.communicate()
    #---# Create the new empty folder
    resDirCommand = "mkdir "+outDirName
    resDir = subprocess.Popen(resDirCommand.split(),stdout=subprocess.PIPE,universal_newlines=True,cwd=resDirPath)
    resDir.communicate()
    time.sleep(1)

    # Create folders for ouput plots per chamber
    SuperChType = runConfig.StandConfiguration
    alignoutDir = os.path.abspath("launcher_alignment.py").split('QC8Test')[0] + outDirName
    for i in range (0,30):
      if (SuperChType[int(i/2)] != '0'):
        plotsDirCommand = "mkdir outPlots_Chamber_Pos_" + str(i)
        plotsDirChamber = subprocess.Popen(plotsDirCommand.split(),stdout=subprocess.PIPE,universal_newlines=True,cwd=alignoutDir)
        plotsDirChamber.communicate()
    time.sleep(1)

    # Selecting the correct output file, changing the name and moving to the output folder
    out_name = 'out_run_'
    for i in range(6-len(str(args.run_number))):
      out_name = out_name + '0'
    out_name = out_name + str(args.run_number) + '.root'

    mvToDirCommand = "mv alignment_" + out_name + " " + resDirPath+outDirName + "/alignment_" + out_name
    movingToDir = subprocess.Popen(mvToDirCommand.split(),stdout=subprocess.PIPE,universal_newlines=True,cwd=runPath)
    movingToDir.communicate()
    time.sleep(1)

    # Alignment computation & output2
    alignCommand = "root -l -q " + runPath + "macro_alignment.c(" + str(args.run_number) + ",\"" + runPath + "\",\"" + alignmentTablesPath + "\"," + str(step) + ")"
    alignment = subprocess.Popen(alignCommand.split(),stdout=subprocess.PIPE,universal_newlines=True,cwd=alignoutDir)
    print line
    while alignment.poll() is None:
      line = alignment.stdout.readline()
      print(line)
    print alignment.stdout.read()
    alignment.communicate()
    time.sleep(1)

    if(docheck):
      checkCommand = "python " + pyhtonModulesPath + "check.py " + str(args.run_number) + " " + str(step)
      check = subprocess.Popen(checkCommand.split(),stdout=subprocess.PIPE,universal_newlines=True,cwd=alignoutDir)
      print line
      while check.poll() is None:
        line = check.stdout.readline()
        print(line)
      print check.stdout.read()
      check.communicate()
      time.sleep(1)

    stop_align = align_stopper(args.run_number, step)
    print stop_align
    step += 1

    # Generate geometry files
    geometry_files_creator.geomMaker(args.run_number, "forAlignment")
    time.sleep(1)

  # Running the CMSSW code for the last step of alignment
  cmsRunner(cores)

  #  # Creating folder outside the CMMSW release to put the output files and plots
  outDirName = "Results_QC8_alignment_run_"+str(args.run_number)
  #---# Remove old version if want to recreate
  if (os.path.exists(resDirPath+outDirName)):
    rmDirCommand = "rm -rf "+outDirName
    rmDir = subprocess.Popen(rmDirCommand.split(),stdout=subprocess.PIPE,universal_newlines=True,cwd=resDirPath)
    rmDir.communicate()
  #---# Create the new empty folder
  resDirCommand = "mkdir "+outDirName
  resDir = subprocess.Popen(resDirCommand.split(),stdout=subprocess.PIPE,universal_newlines=True,cwd=resDirPath)
  resDir.communicate()
  time.sleep(1)

  # Create folders for ouput plots per chamber
  import configureRun_cfi as runConfig
  SuperChType = runConfig.StandConfiguration
  tilttwistoutDir = os.path.abspath("launcher_alignment.py").split('QC8Test')[0] + outDirName
  for i in range (0,30):
    if (SuperChType[int(i/2)] != '0'):
      plotsDirCommand = "mkdir outPlots_Chamber_Pos_" + str(i)
      plotsDirChamber = subprocess.Popen(plotsDirCommand.split(),stdout=subprocess.PIPE,universal_newlines=True,cwd=tilttwistoutDir)
      plotsDirChamber.communicate()
  time.sleep(1)

  # Selecting the correct output file, changing the name and moving to the output folder
  out_name = 'out_run_'
  for i in range(6-len(str(args.run_number))):
    out_name = out_name + '0'
  out_name = out_name + str(args.run_number) + '.root'

  mvToDirCommand = "mv alignment_" + out_name + " " + resDirPath+outDirName + "/alignment_" + out_name
  movingToDir = subprocess.Popen(mvToDirCommand.split(),stdout=subprocess.PIPE,universal_newlines=True,cwd=runPath)
  movingToDir.communicate()
  time.sleep(1)

  # Alignment computation & output
  tilttwistCommand = "root -l -q " + runPath + "macro_tilt_twist.c(" + str(args.run_number) + ",\"" + runPath + "\",\"" + alignmentTablesPath + "\")"
  tilttwist = subprocess.Popen(tilttwistCommand.split(),stdout=subprocess.PIPE,universal_newlines=True,cwd=tilttwistoutDir)
  while tilttwist.poll() is None:
    line = tilttwist.stdout.readline()
    print(line)
  print tilttwist.stdout.read()
  tilttwist.communicate()
  time.sleep(1)

  # Converting tables ToDB-like into FromDB-like
  convertAlignmentTables.convertAlignment(args.run_number,"alignment")
