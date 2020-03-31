#!/usr/bin/env python
from datetime import datetime,date,time
from xml.etree.ElementTree import Element, SubElement, Comment, tostring
import time
import xlrd
from xlrd import xldate
import re
import datetime
import sys
from xmlConversion import generateXMLHeader, generateDataSet, writeToFile,writeToFile1
from xmlConversion import generateXMLDatafastamb,generateXMLDatafast,generateXMLDatalongamb,generateXMLDatalong, generateXMLData3,generateXMLData3a,generateXMLData4,generateXMLData4a,generateXMLData5a,generateXMLData5,generateXMLData4s, generateXMLDataChVfatEfficiency

def xml_from_excel3(excel_file):
    wb = xlrd.open_workbook(excel_file)
    sh = wb.sheet_by_index(0)
    user = ''
    location = ''
    Start = sh.cell(2,1).value
    Stop = str(sh.cell(3,1).value)
    comment=''
    chamber = sh.cell(1,1).value
    Run = sh.cell(0,1).value
    root = generateXMLHeader("QC8_GEM_CH_VFAT_EFFICIENCY","GEM CH VFAT EFFICIENCY QC8","GEM CH VFAT EFFICIENCY",str(Run),str(Start),str(Stop),str(comment),str(location),str(user))
    dataSet = generateDataSet(root,comment,"1","GEM Chamber",str(chamber))
    for row in range(5,sh.nrows):
        vfat_posn= sh.row_values(row)[0]
        efficiency= sh.row_values(row)[1]
        efficiency_error =sh.row_values(row)[2]
        cluster_size_avg =sh.row_values(row)[3]
        cluster_size_sigma =sh.row_values(row)[4]
        percent_masked =sh.row_values(row)[5]
        generateXMLDataChVfatEfficiency(dataSet,str(vfat_posn),str(efficiency), str(efficiency_error),str(cluster_size_avg),str(cluster_size_sigma), str(percent_masked))
        print(str(row-4)+" / "+str(sh.nrows-5))
    fileName=fname[:-5]+"_Run_"+str(Run)+".xml"
    writeToFile(fileName, tostring(root))
if __name__ =="__main__":
	fname = sys.argv[1]
	xml_from_excel3(fname)
