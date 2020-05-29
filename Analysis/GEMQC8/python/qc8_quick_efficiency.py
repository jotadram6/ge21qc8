#!/usr/bin/env python
from datetime import datetime,date,time
from time import sleep
from xml.etree.ElementTree import Element, SubElement, Comment, tostring
import time
import xlrd
from xlrd import xldate
import re
import sys
from xmlConversion import generateXMLHeader, generateDataSet, writeToFile,writeToFile1
from xmlConversion import generateXMLDatafastamb,generateXMLDatafast,generateXMLDatalongamb,generateXMLDatalong, generateXMLData3,generateXMLData3a,generateXMLData4,generateXMLData4a,generateXMLData5a,generateXMLData5,generateXMLData4s, generateXMLDataChVfatEfficiency, generateXMLDataQuickEfficiencyQC8

def xml_from_excel3(excel_file):
    wb = xlrd.open_workbook(excel_file)
    sh = wb.sheet_by_index(0)
    user = ''
    location = ''
    Start = ''
    Stop = ''
    comment = ''
    chamber = sh.cell(1,1).value
    overall_efficiency = sh.cell(2,1).value
    error_efficiency = sh.cell(3,1).value
    Run = sh.cell(0,1).value
    root = generateXMLHeader("QC8_GEM_QUICK_EFFICIENCY","GEM QUICK EFFICIENCY QC8","GEM QUICK EFFICIENCY",str(Run),str(Start),str(Stop),str(comment),str(location),str(user))
    dataSet = generateDataSet(root,comment,"1","GEM Chamber",str(chamber))
    generateXMLDataQuickEfficiencyQC8(dataSet,str(overall_efficiency),str(error_efficiency))
    fileName=fname[:-5]+"_Run_"+str(Run)+".xml"
    writeToFile(fileName, tostring(root))
if __name__ =="__main__":
	fname = sys.argv[1]
	xml_from_excel3(fname)
