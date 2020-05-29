#!/usr/bin/env python
from datetime import datetime,date,time
from time import sleep
from xml.etree.ElementTree import Element, SubElement, Comment, tostring
import time
import xlrd
from xlrd import xldate
import re
import sys
from xmlConversion import generateXMLHeader, generateDataSetMultipleParts, writeToFile,writeToFile1
from xmlConversion import generateXMLDatafastamb,generateXMLDatafast,generateXMLDatalongamb,generateXMLDatalong, generateXMLData3,generateXMLData3a,generateXMLData4,generateXMLData4a,generateXMLData5a,generateXMLData5,generateXMLData4s, generateXMLData5s, generateXMLDataStrips, generateXMLDataAlignment, generateXMLDataQC8DeadStrips

def xml_from_excel5(excel_file):
    wb = xlrd.open_workbook(excel_file)
    sh = wb.sheet_by_index(0)
    user = ''
    location=''
    Start=''
    Stop=''
    comment=''
    Run = sh.cell(0,1).value
    root = generateXMLHeader("QC8_GEM_MASKED_STRIPS_DEAD","GEM QC8 MASKED STRIPS DEAD", "GEM DEAD STRIPS QC8",str(Run),str(Start),str(Stop),str(comment),str(location),str(user))
    dataSet = generateDataSetMultipleParts(root,comment,"1")
    for row in range(2,sh.nrows):
        ch_serial_number= sh.row_values(row)[0]
        gem_number= sh.row_values(row)[1]
        position= sh.row_values(row)[2]
        vfat= sh.row_values(row)[3]
        channel= sh.row_values(row)[4]
        strip= sh.row_values(row)[5]
        generateXMLDataQC8DeadStrips(dataSet,str(ch_serial_number),str(gem_number),str(position), str(vfat), str(channel),str(strip))
        print(str(row-1)+" / "+str(sh.nrows-2))
    writeToFile(fileName,tostring(root))
if __name__ =="__main__":
	fname = sys.argv[1]
	fileName=fname[:-5]+".xml"
	xml_from_excel5(fname)
