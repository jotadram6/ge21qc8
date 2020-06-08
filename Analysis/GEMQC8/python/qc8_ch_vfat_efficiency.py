#!/usr/bin/env python
from xml.etree.ElementTree import Element, SubElement, Comment, tostring
import xlrd
import sys
from xmlConversion import generateXMLHeader, generateDataSet, writeToFile, readFile
from xmlConversion import generateXMLDataChVfatEfficiency

def xml_from_csv(csv_file):
  lines = readFile(csv_file)
  for i in range(len(lines)):
    lines[i] = lines[i].split('\n')[0]
  Run = int(lines[0].split(',')[1])
  chamber = str(lines[1].split(',')[1])
  Start = str(lines[2].split(',')[1])
  Stop = str(lines[3].split(',')[1])
  user = ''
  location = ''
  comment = ''
  root = generateXMLHeader("QC8_GEM_CH_VFAT_EFFICIENCY","GEM CH VFAT EFFICIENCY QC8","GEM CH VFAT EFFICIENCY",str(Run),str(Start),str(Stop),str(comment),str(location),str(user))
  dataSet = generateDataSet(root,comment,"1","GEM Chamber",str(chamber))
  for line in range(5,len(lines)):
    vfat_posn = int(lines[line].split(',')[0])
    efficiency = float(lines[line].split(',')[1])
    efficiency_error = float(lines[line].split(',')[2])
    cluster_size_avg = float(lines[line].split(',')[3])
    cluster_size_sigma = float(lines[line].split(',')[4])
    percent_masked = float(lines[line].split(',')[5])
    generateXMLDataChVfatEfficiency(dataSet,str(vfat_posn),str(efficiency), str(efficiency_error),str(cluster_size_avg),str(cluster_size_sigma), str(percent_masked))
    print(str(line-4)+" / "+str(len(lines)-5))
    fileName = csv_file[:-4] + "_Run_" + str(Run) + ".xml"
    writeToFile(fileName, tostring(root))

if __name__ == "__main__":
	fname = sys.argv[1]
	xml_from_csv(fname)
