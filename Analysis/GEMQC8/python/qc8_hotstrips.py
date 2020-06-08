#!/usr/bin/env python
from xml.etree.ElementTree import Element, SubElement, Comment, tostring
import xlrd
import sys
from xmlConversion import generateXMLHeader, generateDataSetMultipleParts, writeToFile, readFile
from xmlConversion import generateXMLDataQC8DeadStrips

def xml_from_csv(csv_file):
  lines = readFile(csv_file)
  for i in range(len(lines)):
    lines[i] = lines[i].split('\n')[0]
  Run = int(lines[0].split(',')[1])
  Start = '2020-03-14 12:00:00'
  Stop = '2020-03-14 12:00:00'
  comment = ''
  location = ''
  user = ''
  root = generateXMLHeader("QC8_GEM_MASKED_STRIPS_HOT","GEM QC8 MASKED STRIPS HOT", "GEM HOT STRIPS QC8",str(Run),str(Start),str(Stop),str(comment),str(location),str(user))
  dataSet = generateDataSetMultipleParts(root,comment,"1")
  for line in range(2,len(lines)):
    ch_serial_number = str(lines[line].split(',')[0])
    gem_number = str(lines[line].split(',')[1])
    position = str(lines[line].split(',')[2])
    vfat = int(lines[line].split(',')[3])
    channel = int(lines[line].split(',')[4])
    strip = int(lines[line].split(',')[5])
    generateXMLDataQC8DeadStrips(dataSet,str(ch_serial_number),str(gem_number),str(position), str(vfat), str(channel),str(strip))
    print(str(line-1)+" / "+str(len(lines)-2))
  fileName = csv_file[:-4] + ".xml"
  writeToFile(fileName,tostring(root))

if __name__ == "__main__":
	fname = sys.argv[1]
	xml_from_csv(fname)
