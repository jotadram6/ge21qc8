from __future__ import division
import argparse
from xml.etree.ElementTree import Element, SubElement, Comment, tostring
from xml.dom import minidom
from xml.etree import ElementTree
from datetime import datetime
import sys

def writeToFile(path, lines):
    with open(path, mode='w+') as myfile:
        myfile.write(lines)

def readFile(path):
    with open(path) as f:
	    lines=f.readlines()
    return lines

def generateXMLHeader(extensionTableNameText, nameText, runTypeText, runNumberText, runBeginText, runEndText, commentText, locationText, userText):
    root = Element('ROOT')
    root.set('xmlns:xsi', 'http://www.w3.org/2001/XMLSchema-instance')
    header = SubElement(root, 'HEADER')
    type = SubElement(header, 'TYPE')
    extensionTableName = SubElement(type, 'EXTENSION_TABLE_NAME')
    extensionTableName.text = extensionTableNameText
    name = SubElement(type, 'NAME')
    name.text = nameText
    run = SubElement(header, 'RUN')
    runType= SubElement(run, 'RUN_TYPE')
    runType.text= runTypeText
    runNumber = SubElement(run, 'RUN_NUMBER')
    runNumber.text= runNumberText
    runBegin = SubElement(run, 'RUN_BEGIN_TIMESTAMP')
    runBegin.text = runBeginText
    runEnd = SubElement(run,'RUN_END_TIMESTAMP')
    runEnd.text= runEndText
    comment = SubElement(run, 'COMMENT_DESCRIPTION')
    comment.text= commentText
    location = SubElement(run, 'LOCATION')
    location.text = locationText
    user = SubElement(run,'INITIATED_BY_USER')
    user.text = userText
    return root

def generateDataSet(Set,descriptionText,versionText,kindText,serialText):
    dataSet = SubElement(Set, 'DATA_SET')
    description = SubElement(dataSet, 'COMMENT_DESCRIPTION')
    description.text = descriptionText
    version = SubElement(dataSet,'VERSION')
    version.text = versionText
    part = SubElement(dataSet,"PART")
    kind = SubElement(part,"KIND_OF_PART")
    kind.text = kindText
    serial = SubElement(part,"SERIAL_NUMBER")
    serial.text = serialText
    return dataSet

def generateDataSetMultipleParts(Set,descriptionText,versionText):
    dataSet = SubElement(Set, 'DATA_SET')
    description = SubElement(dataSet, 'COMMENT_DESCRIPTION')
    description.text = descriptionText
    version = SubElement(dataSet,'VERSION')
    version.text = versionText
    return dataSet

def generateXMLDataAlignment(dataSetTag,positionText,dxText,dyText,dzText,rxText,ryText,rzText):
    data = SubElement(dataSetTag, 'DATA')
    position = SubElement(data, 'POSITION')
    position.text = positionText
    dx = SubElement(data, 'DX')
    dx.text=dxText
    dy = SubElement(data, 'DY')
    dy.text = dyText
    dz = SubElement(data, 'DZ')
    dz.text=dzText
    rx=SubElement(data,'RX')
    rx.text=rxText
    ry=SubElement(data, 'RY')
    ry.text=ryText
    rz=SubElement(data,'RZ')
    rz.text=rzText

def generateXMLDataQC8DeadStrips(dataSetTag,ch_serial_numberText,gem_numberText,positionText,vfatText,channelText,stripText):
    data = SubElement(dataSetTag, 'DATA')
    ch_serial_number = SubElement(data, 'CH_SERIAL_NUMBER')
    ch_serial_number.text = ch_serial_numberText
    gem_number = SubElement(data, 'GEM_NUMBER')
    gem_number.text=gem_numberText
    position = SubElement(data, 'POSITION')
    position.text = positionText
    vfat = SubElement(data, 'VFAT')
    vfat.text=vfatText
    channel=SubElement(data,'CHANNEL')
    channel.text=channelText
    strip=SubElement(data, 'STRIP')
    strip.text=stripText

def generateXMLDataChVfatEfficiency(dataSetTag,vfat_posnText,efficiencyText,efficiency_errorText, cluster_size_avgText,cluster_size_sigmaText, percent_maskedText):
    data = SubElement(dataSetTag, 'DATA')
    vfat_posn = SubElement(data, 'VFAT_POSN')
    vfat_posn.text = vfat_posnText
    efficiency = SubElement(data, 'EFFICIENCY')
    efficiency.text=efficiencyText
    efficiency_error = SubElement(data, 'EFFICIENCY_ERROR')
    efficiency_error.text = efficiency_errorText
    cluster_size_avg = SubElement(data, 'CLUSTER_SIZE_AVG')
    cluster_size_avg.text = cluster_size_avgText
    cluster_size_sigma = SubElement(data, 'CLUSTER_SIZE_SIGMA')
    cluster_size_sigma.text = cluster_size_sigmaText
    percent_masked = SubElement(data, 'PERCENT_MASKED')
    percent_masked.text = percent_maskedText

def generateXMLDataStandGeoConf(dataSetTag,ch_serialText,positionText,flowmeterText):
    data = SubElement(dataSetTag, 'DATA')
    ch_serial_number = SubElement(data, 'CH_SERIAL_NUMBER')
    ch_serial_number.text = ch_serialText
    position = SubElement(data, 'POSITION')
    position.text=positionText
    flowmeter = SubElement(data, 'FLOW_METER')
    flowmeter.text = flowmeterText

def generateXMLDataQuickEfficiencyQC8(dataSetTag,overall_efficiencyText,error_efficiencyText):
    data = SubElement(dataSetTag, 'DATA')
    overall_efficiency = SubElement(data, 'OVERALLEFFICIENCY')
    overall_efficiency.text = overall_efficiencyText
    error_efficiency = SubElement(data, 'ERROREFFICIENCY')
    error_efficiency.text=error_efficiencyText
