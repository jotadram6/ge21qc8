#!/usr/bin/env python
from xml.etree.ElementTree import Element, SubElement, Comment, tostring
import xlrd
import sys
from xmlConversion import generateXMLHeader, generateDataSet, writeToFile, readFile
from xmlConversion import generateXMLDataQuickEfficiencyQC8

def xml_from_csv(csv_file):
  lines = readFile(csv_file)
  for i in range(len(lines)):
    lines[i] = lines[i].split('\n')[0]
  Run = int(lines[0].split(',')[1])
  chamber = str(lines[1].split(',')[1])
  overall_efficiency = float(lines[2].split(',')[1])
  error_efficiency = float(lines[3].split(',')[1])
  Start = '2020-03-14 12:00:00'
  Stop = '2020-03-14 12:00:00'
  comment = ''
  location = ''
  user = ''
  root = generateXMLHeader("QC8_GEM_QUICK_EFFICIENCY","GEM QUICK EFFICIENCY QC8","GEM QUICK EFFICIENCY",str(Run),str(Start),str(Stop),str(comment),str(location),str(user))
  dataSet = generateDataSet(root,comment,"1","GEM Chamber",str(chamber))
  generateXMLDataQuickEfficiencyQC8(dataSet,str(overall_efficiency),str(error_efficiency))
  fileName = csv_file[:-4] + "_Run_" + str(Run) + ".xml"
  writeToFile(fileName, tostring(root))

if __name__ == "__main__":
	fname = sys.argv[1]
	xml_from_csv(fname)
