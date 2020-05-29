import csv
import os
import sys
import io
import subprocess
import time
import datetime
import socket
import shutil

if __name__ == '__main__':

  # Define the parser
  import argparse
  parser = argparse.ArgumentParser(description="QC8 data analysis step 5. Software alignment of the chambers in the stand. For any doubt: https://twiki.cern.ch/twiki/bin/view/CMS/GEMCosmicRayAnalysis")
  # Positional arguments
  parser.add_argument("input_run", type=int, help="Specify the run number of the input file to be replicated")
  parser.add_argument("output_run", type=str, help="Specify the run number(s) for which replicate the alignment file. If you want to replicate the same alignment file for more than one run, you can write the list separated by commas \"1,2,3,4\", as an interval \"1-4\", or a combination of the two")
  args = parser.parse_args()

  tablePath = os.path.abspath("replicate_alignment_table.py").split('QC8Test')[0] + 'QC8Test/src/Analysis/GEMQC8/data/StandAligmentTables/'

  input_name = tablePath + "StandAlignmentValues_run{}_ToDB.csv".format(int(args.input_run))

  out_runs = []
  commas = args.output_run.count(',')
  for comma in range(commas+1):
    split = args.output_run.split(',')[comma]
    if ('-' in split):
      start = int(split.split('-')[0])
      stop = int(split.split('-')[1])
      for run in range(start,stop+1):
        out_runs.append(int(run))
    else:
      out_runs.append(int(split))

  for run in out_runs:
    # Create a copy with new file name
    output_name = tablePath + "StandAlignmentValues_run{}_ToDB.csv".format(int(run))
    shutil.copy(input_name,output_name)

    # Replace run number in the file
    f = open(output_name)
    contents = f.read()
    newcontents = contents.replace("RunNumber,{},,,,,".format(int(args.input_run)),"RunNumber,{},,,,,".format(int(run)))
    f.close()
    f = open(output_name,"w")
    f.write(newcontents)
    f.close()

    # Converting tables ToDB-like into FromDB-like
    import convertAlignmentTables
    convertAlignmentTables.convertAlignment(run,"alignment")
