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

def writeToFile1(path, lines):
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
def generateDataSetBlob(Set,descriptionText,file_nameText,versionText,kindText,serialText):
    dataSet = SubElement(Set, 'DATA_SET')
    description = SubElement(dataSet, 'COMMENT_DESCRIPTION')
    description.text = descriptionText
    file_name = SubElement(dataSet, 'DATA_FILE_NAME')
    file_name.text = file_nameText
    version = SubElement(dataSet,'VERSION')
    version.text = versionText
    part = SubElement(dataSet,"PART")
    kind = SubElement(part,"KIND_OF_PART")
    kind.text = kindText
    serial = SubElement(part,"SERIAL_NUMBER")
    serial.text = serialText
    return dataSet
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
#QC2
def generateXMLDatafastamb(dataSetTag,HumidityText):
    data = SubElement(dataSetTag, 'DATA')
    Humidity = SubElement(data, 'HUMIDITY_PERCENT')
    Humidity.text = HumidityText

def generateXMLDatafast(dataSetTag,timeMinutesText,appliedVoltageText,impedanceText,leakageText,sparksText,totsparksText):
    data = SubElement(dataSetTag, 'DATA')
    timeMinutes = SubElement(data, 'TIME_MINUTES')
    timeMinutes.text = timeMinutesText
    appliedVoltage = SubElement(data, 'APPLIED_VOLTAGE')
    appliedVoltage.text=appliedVoltageText
    impedance=SubElement(data,'IMPEDANCE_GOHMS')
    impedance.text=impedanceText
    leakage=SubElement(data, 'LEAKAGE_CURRENT_NA')
    leakage.text=leakageText
    sparks=SubElement(data,'NUM_SPARKS')
    sparks.text=sparksText
    totspark=SubElement(data, 'TOT_NUM_SPARKS')
    totspark.text=totsparksText

def generateXMLDatalongamb(dataSetTag,batchText,RelhumiText,TempdegText):
    data = SubElement(dataSetTag, 'DATA')
    batch = SubElement(data, 'BATCH')
    batch.text = batchText
    Relhumi = SubElement(data, 'REL_HUM_PRCNT')
    Relhumi.text = RelhumiText
    Tempdeg = SubElement(data, 'TEMP_DEG_C')
    Tempdeg.text = TempdegText
   
def generateXMLDatalong(dataSetTag,timesecext,appliedVoltageText,errvoltText,currentText,ErrText):
    data = SubElement(dataSetTag, 'DATA')
    timesec = SubElement(data, 'TIME_SEC')
    timesec.text = timesecText
    appliedVoltage = SubElement(data, 'HV_VOLTS')
    appliedVoltage.text=appliedVoltageText
    errvolt=SubElement(data,'ERR_HV_VOLTS')
    errvolt.text=errvoltText
    current=SubElement(data, 'CURRENT_MICRO_AMP')
    current.text=currentText
    Err=SubElement(data,'ERR_CURNT_MICRO_AMP')
    Err.text=ErrText
#QC3    
def generateXMLData3(dataSetTag,testperformText,timeincrementText,mainpressureText,ambpressureText,temperatureText, increment_hrText, temp_kText):
    data = SubElement(dataSetTag, 'DATA')
    testperform = SubElement(data, 'TEST_TIME')
    testperform.text = testperformText
    timeincrement = SubElement(data, 'INCRMNT_SEC')
    timeincrement.text=timeincrementText
    mainpressure=SubElement(data,'MANF_PRSR_MBAR')
    mainpressure.text=mainpressureText
    ambpressure=SubElement(data, 'AMB_PRSR_MBAR')
    ambpressure.text=ambpressureText
    temperature=SubElement(data,'TEMP_DEGC')
    temperature.text=temperatureText
    increment_hr = SubElement(data, 'INCRMNT_HR')
    increment_hr.text = increment_hrText
    temp_k = SubElement(data, 'TEMP_K')
    temp_k.text = temp_kText
    
def generateXMLData3a(dataSetTag,testperformText,AvgambText,StdambText,AvgpreText,stdpreText,initpreText,finalpreText,durationText,leakText,expoText,expoleakText,elogText,filenameText,commentText,time_constantText):
    data = SubElement(dataSetTag, 'DATA')
    testperform = SubElement(data, 'TEST_DATE')
    testperform.text = testperformText
    Avgamb = SubElement(data, 'AVG_AMBTEMP_DEGC')
    Avgamb.text=AvgambText
    Stdamb = SubElement(data, 'STDDEV_AMBTEMP_DEGC')
    Stdamb.text = StdambText
    Avgpre = SubElement(data, 'AVG_AMBPRSR_MBAR')
    Avgpre.text=AvgpreText
    stdpre=SubElement(data,'STDDEV_AMBPRSR_MBAR')
    stdpre.text=stdpreText
    initpre=SubElement(data, 'INIT_PRSR_MBAR')
    initpre.text=initpreText
    finalpre=SubElement(data,'FINAL_PRSR_MBAR')
    finalpre.text=finalpreText
    duration=SubElement(data,'DURATION_HR')
    duration.text=durationText
    leak=SubElement(data, 'LEAK_RATE_MBAR_HR')
    leak.text=leakText
    expo=SubElement(data,'EXPO_FIT_P0_MBAR')
    expo.text=expoText
    expoleak=SubElement(data,'EXPO_FIT_LEAK_PARAM')
    expoleak.text=expoleakText
    elog=SubElement(data,'ELOG_LINK')
    elog.text=elogText
    filename=SubElement(data,'FILE_NAME')
    filename.text=filenameText
    comment=SubElement(data,'COMMENTS')
    comment.text=commentText
    time_constant= SubElement(data, 'TIME_CONSTANT')
    time_constant.text = time_constantText
#QC4
def generateXMLData4(dataSetTag,timeMinutesText,appliedVoltageText,impedanceText,leakageText,sparksText,rnormText,totsparksText,rateText,errorText):
    data = SubElement(dataSetTag, 'DATA')
    timeMinutes = SubElement(data, 'VSET_VLT')
    timeMinutes.text = timeMinutesText
    appliedVoltage = SubElement(data, 'VMON_VLT')
    appliedVoltage.text=appliedVoltageText
    impedance=SubElement(data,'ISET_UA')
    impedance.text=impedanceText
    leakage=SubElement(data, 'IMON_UA')
    leakage.text=leakageText
    sparks=SubElement(data,'RCALC_MOHM')
    sparks.text=sparksText
    rnorm=SubElement(data, 'RNORM_MOHM')
    rnorm.text=rnormText
    totsparks=SubElement(data,'COUNTS')
    totsparks.text=totsparksText
    rate = SubElement(data, 'RATE_HZ')
    rate.text =rateText 
    error = SubElement(data, 'ERR_RATE_HZ')
    error.text=errorText
#QC4 Configuration
def generateXMLData4a(dataSetTag,timeMinutesText,appliedVoltageText,impedanceText,leakageText,sparksText,rnormText,totsparksText,rateText,errorText,scalText,daqText):
    data = SubElement(dataSetTag, 'DATA')
    timeMinutes = SubElement(data, 'REQUIV_MOHM_MSRD')
    timeMinutes.text = timeMinutesText
    appliedVoltage = SubElement(data, 'PREAMPLIFIER')
    appliedVoltage.text=appliedVoltageText
    impedance=SubElement(data,'AMPLIFIER')
    impedance.text=impedanceText
    leakage=SubElement(data, 'COARSE_GAIN')
    leakage.text=leakageText
    sparks=SubElement(data,'FINE_GAIN')
    sparks.text=sparksText
    rnorm=SubElement(data, 'INT_TIME_NS')
    rnorm.text=rnormText
    totsparks=SubElement(data,'DIFF_TIME_NS')
    totsparks.text=totsparksText
    rate = SubElement(data, 'DISCRIMINATOR')
    rate.text =rateText 
    error = SubElement(data, 'THRSHLD_MV')
    error.text=errorText
    scal = SubElement(data, 'SCALAR')
    scal.text=scalText
    daq = SubElement(data, 'DAQ_TIME_SEC')
    daq.text = daqText
def generateXMLData4s(dataSetTag,testperformText,VmaxText,CurrVmaxText,VdrftText,ReqmsrdText,ReqslpText,diffText,SprsglText,filenameText,elogText,commentText,SprsglErrorText):
    data = SubElement(dataSetTag, 'DATA')
    testperform = SubElement(data, 'TEST_DATE')
    testperform.text = testperformText
    Vmax = SubElement(data, 'VMAX_VLT')
    Vmax.text=VmaxText
    CurrVmax = SubElement(data, 'CURNT_AT_VMAX_UA')
    CurrVmax.text = CurrVmaxText
    Vdrft = SubElement(data, 'VDRFT_VLT')
    Vdrft.text=VdrftText
    Reqmsrd=SubElement(data,'REQUIV_MOHM_MSRD')
    Reqmsrd.text=ReqmsrdText
    Reqslp=SubElement(data, 'REQUIV_MOHM_SLP')
    Reqslp.text=ReqslpText
    diff=SubElement(data,'DIFF_REQUIV_PRCNT')
    diff.text=diffText
    Sprsgl=SubElement(data,'SPR_SGNL_AT_VMAX_HZ')
    Sprsgl.text=SprsglText
    filename=SubElement(data,'FILE_NAME')
    filename.text=filenameText
    elog=SubElement(data,'ELOG_LINK')
    elog.text=elogText
    comment=SubElement(data,'COMMENTS')
    comment.text=commentText
    SprsglError=SubElement(data,'SPR_SGNL_ERR')
    SprsglError.text=SprsglErrorText	
#QC5
def generateXMLData5(dataSetTag,timeMinutesText,tempText,pressText,imonText,vmonText,vdrifText,rateText,errorText,currentText,cerrorText,gainText,gerrorText):
    data = SubElement(dataSetTag, 'DATA')
    timeMinutes = SubElement(data, 'TIME')
    timeMinutes.text = timeMinutesText
    temp = SubElement(data, 'TEMP_DEGC')
    temp.text=tempText
    press=SubElement(data,'PRESSURE_MBAR')
    press.text=pressText
    imon=SubElement(data,'IMON_UA')
    imon.text=imonText
    vmon=SubElement(data, 'VMON_VLT')
    vmon.text=vmonText
    vdrif=SubElement(data,'VDRIFT_VLT')
    vdrif.text=vdrifText
    rate = SubElement(data, 'RATE_HZ')
    rate.text =rateText 
    error = SubElement(data, 'RATE_ERROR_HZ')
    error.text=errorText
    current = SubElement(data, 'CURRENT_AMP')
    current.text=currentText
    cerror = SubElement(data, 'CURRENT_ERROR_AMP')
    cerror.text=cerrorText
    gain = SubElement(data, 'GAIN')
    gain.text=gainText
    gerror = SubElement(data, 'GAIN_ERROR')
    gerror.text=gerrorText
#QC5 Configuration
def generateXMLData5a(dataSetTag,userText,ampText,cgainText,fgainText,itimeText,dtimeText,discText,thrsText,wadjsText,widthText, scalText,daqText,picoText,souText,hvltText,currentText,nbpriText,etaText,gasText,gfracText,flowText,reqText,diviText,activityText,filterstatusText,collimatorstatusText):
    data = SubElement(dataSetTag, 'DATA')
    user = SubElement(data, 'USER_NAME')
    user.text = userText
    amp=SubElement(data,'AMPLIFIER')
    amp.text=ampText
    cgain=SubElement(data, 'COARSE_GAIN')
    cgain.text=cgainText
    fgain=SubElement(data,'FINE_GAIN')
    fgain.text=fgainText
    itime=SubElement(data, 'INT_TIME_NS')
    itime.text=itimeText
    dtime=SubElement(data,'DIFF_TIME_NS')
    dtime.text=dtimeText
    disc = SubElement(data, 'DISCRIMINATOR')
    disc.text =discText 
    thrs = SubElement(data, 'THRSHLD_MV')
    thrs.text=thrsText
    wadjs = SubElement(data, 'WALK_ADJUST_MV')
    wadjs.text = wadjsText
    width = SubElement(data, 'WIDTH_NS')
    width.text = widthText 
    scal = SubElement(data, 'SCALAR')
    scal.text=scalText
    daq = SubElement(data, 'DAQ_TIME_SEC')
    daq.text = daqText
    pico = SubElement(data, 'PICOAMMETER')
    pico.text = picoText
    sou = SubElement(data, 'SOURCE')
    sou.text = souText
    hvlt = SubElement(data, 'HV_VLT')
    hvlt.text = hvltText
    current = SubElement(data, 'CURRENT_UA')
    current.text = currentText
    nbpri = SubElement(data, 'NB_PRIMARIES')
    nbpri.text = nbpriText
    eta = SubElement(data, 'ETA_SECTOR')
    eta.text = etaText
    gas = SubElement(data, 'GAS')
    gas.text = gasText
    gfrac = SubElement(data, 'GAS_FRAC')
    gfrac.text = gfracText
    flow = SubElement(data, 'FLOW_RATE_LHR')
    flow.text = flowText
    req = SubElement(data,'REQUIV_MOHM_MSRD')
    req.text = reqText
    divi = SubElement(data, 'REDIVIDER_MOHM_MSRD')
    divi.text = diviText
    activity = SubElement(data, 'ACTIVITY_BQ')
    activity.text = activityText
    filterstatus = SubElement(data, 'FILTER_STATUS')
    filterstatus.text = filterstatusText
    collimatorstatus = SubElement(data, 'COLLIMATOR_STATUS')
    collimatorstatus.text = collimatorstatusText
def generateXMLData5s(dataSetTag,timeText, AvgambText,AvgpreText,expofitText,expofit2Text,filesText,elogText,commentText, errorTempText, errorPressureText,p5p0Text,p5t0Text):
    data = SubElement(dataSetTag, 'DATA')
    time = SubElement(data, 'TEST_DATE')
    time.text=timeText
    Avgamb = SubElement(data, 'AVG_TEMP_DEGC')
    Avgamb.text=AvgambText
    Avgpre=SubElement(data,'AVG_PRESSURE_MBAR')
    Avgpre.text=AvgpreText
    expofit=SubElement(data,'EXPO_FIT1_PARA1')
    expofit.text=expofitText
    expofit2=SubElement(data,'EXPO_FIT1_PARA2')
    expofit2.text=expofit2Text
    files=SubElement(data,'FILE_NAME')
    files.text=filesText
    elog=SubElement(data,'ELOG_LINK')
    elog.text=elogText
    comment=SubElement(data,'COMMENTS')
    comment.text=commentText
    errorTempt=SubElement(data,'ERR_TEMP_DEGC')
    errorTempt.text=errorTempText
    errorPressure=SubElement(data,'ERR_PRESSURE_MBAR')
    errorPressure.text=errorPressureText
    p5p0=SubElement(data,'P5_P0')
    p5p0.text=p5p0Text
    p5t0=SubElement(data,'P5_T0')
    p5t0.text=p5t0Text
def generateXMLData5unif(dataSetTag,ROOT_DATA_FILE_Text,filesText,elogText,commentText):
    data = SubElement(dataSetTag, 'DATA')
    files=SubElement(data,'FILE_NAME')
    files.text=filesText
    elog=SubElement(data,'ELOG_LINK')
    elog.text=elogText
    comment=SubElement(data,'COMMENTS')
    comment.text=commentText
    root_file = SubElement(data, 'ROOT_DATA_FILE')
    root_file.text = ROOT_DATA_FILE_Text
def generateXMLData5unif_data(dataSetTag,PLOT_DATA_FILE_Text):
    data = SubElement(dataSetTag, 'DATA')
    plot_file = SubElement(data, 'PLOT_DATA_FILE')
    plot_file.text = PLOT_DATA_FILE_Text
def generateXMLData6s(dataSetTag,testperformText,VmaxText,CurrVmaxText,VdrftText,ReqmsrdText,ReqslpText,diffText,sprsglText,tripText,filenameText,elogText,commentText):
    data = SubElement(dataSetTag, 'DATA')
    testperform = SubElement(data, 'STARTING_TIME')
    testperform.text = testperformText
    Vmax = SubElement(data, 'ENDING_TIME')
    Vmax.text=VmaxText
    CurrVmax = SubElement(data, 'STABILITY_TIME_HR')
    CurrVmax.text = CurrVmaxText
    Vdrft = SubElement(data, 'STABILITY_VSET_VLT')
    Vdrft.text=VdrftText
    Reqmsrd=SubElement(data,'STABILITY_VMON_VLT')
    Reqmsrd.text=ReqmsrdText
    Reqslp=SubElement(data, 'STABILITY_ISET_UA')
    Reqslp.text=ReqslpText
    diff=SubElement(data,'STABILITY_IMON_UA')
    diff.text=diffText
    Sprsgl=SubElement(data,'STABILITY_VDRFT_VLT')
    Sprsgl.text=SprsglText
    trip= SubElement(data,'TRIP_TIME_SEC')
    trip.text=tripText
    filename=SubElement(data,'FILE_NAME')
    filename.text=filenameText
    elog=SubElement(data,'ELOG_LINK')
    elog.text=elogText
    comment=SubElement(data,'COMMENTS')
    comment.text=commentText
def generateXMLDataStrips(dataSetTag,sectorText,depthText,positionText,chText,hot_chText,fit_failedText,dead_chText,high_noiseText,high_effText):
    data = SubElement(dataSetTag, 'DATA')
    sector = SubElement(data, 'SECTOR')
    sector.text = sectorText
    depth = SubElement(data, 'DEPTH')
    depth.text=depthText
    position = SubElement(data, 'VFAT_POSN')
    position.text = positionText
    ch = SubElement(data, 'VFAT_CHAN')
    ch.text=chText
    hot_ch=SubElement(data,'HOT_CHAN')
    hot_ch.text=hot_chText
    fit_failed=SubElement(data, 'FIT_FAILED')
    fit_failed.text=fit_failedText
    dead_ch=SubElement(data,'DEAD_CHAN')
    dead_ch.text=dead_chText
    high_noise=SubElement(data,'HIGH_NOISE')
    high_noise.text=high_noiseText
    high_eff= SubElement(data,'HIGH_EFF_PED')
    high_eff.text=high_effText
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
def generateXMLData6a(dataSetTag,reqText):
    data = SubElement(dataSetTag, 'DATA')
    req = SubElement(data, 'REQUIV_MOHM_MSRD')
    req.text=reqText
def generateXMLData6b(dataSetTag,vmon_equ_vltText, imon_equ_uaText,vmon_g3b_vltText,imon_g3b_uaText, vmon_g3t_vltText, imon_g3t_uaText,vmon_g2b_vltText,imon_g2b_uaText,vmon_g2t_vltText,imon_g2t_uaText,vmon_g1b_vltText,imon_g1b_uaText,vmon_g1t_vltText,imon_g1t_uaText,vmon_drift_vltText,imon_drift_uaText):
    data = SubElement(dataSetTag, 'DATA')
    VMON_EQU_VLT = SubElement(data, 'VMON_EQU_VLT')
    VMON_EQU_VLT.text = vmon_equ_vltText
    IMON_EQU_UA = SubElement(data, 'IMON_EQU_UA')
    IMON_EQU_UA.text=imon_equ_uaText
    VMON_G3B_VLT = SubElement(data, 'VMON_G3B_VLT')
    VMON_G3B_VLT.text = vmon_g3b_vltText
    IMON_G3B_UA = SubElement(data, 'IMON_G3B_UA')
    IMON_G3B_UA.text = imon_g3b_uaText
    VMON_G3T_VLT=SubElement(data,'VMON_G3T_VLT')
    VMON_G3T_VLT.text= vmon_g3t_vltText
    IMON_G3T_UA=SubElement(data, 'IMON_G3T_UA')
    IMON_G3T_UA.text= imon_g3t_uaText
    VMON_G2B_VLT=SubElement(data,'VMON_G2B_VLT')
    VMON_G2B_VLT.text= vmon_g2b_vltText
    IMON_G2B_UA=SubElement(data,'IMON_G2B_UA')
    IMON_G2B_UA.text=imon_g2b_uaText
    VMON_G2T_VLT= SubElement(data,'VMON_G2T_VLT')
    VMON_G2T_VLT.text=vmon_g2t_vltText
    IMON_G2T_UA=SubElement(data,'IMON_G2T_UA')
    IMON_G2T_UA.text=imon_g2t_uaText
    VMON_G1B_VLT=SubElement(data,'VMON_G1B_VLT')
    VMON_G1B_VLT.text=vmon_g1b_vltText
    IMON_G1B_UA=SubElement(data,'IMON_G1B_UA')
    IMON_G1B_UA.text=imon_g1b_uaText
    VMON_G1T_VLT= SubElement(data,'VMON_G1T_VLT')
    VMON_G1T_VLT.text=vmon_g1t_vltText
    IMON_G1T_UA=SubElement(data,'IMON_G1T_UA')
    IMON_G1T_UA.text=imon_g1t_uaText
    VMON_DRIFT_VLT=SubElement(data,'VMON_DRIFT_VLT')
    VMON_DRIFT_VLT.text=vmon_drift_vltText
    IMON_DRIFT_UA=SubElement(data,'IMON_DRIFT_UA')
    IMON_DRIFT_UA.text=imon_drift_uaText
def generateXMLData6c(dataSetTag,test_dateText,filenameText,elogText,commentText):
    data = SubElement(dataSetTag, 'DATA')
    test_date = SubElement(data, 'TEST_DATE')
    test_date.text = test_dateText
    filename=SubElement(data,'FILE_NAME')
    filename.text=filenameText
    elog=SubElement(data,'ELOG_LINK')
    elog.text=elogText
    comment=SubElement(data,'COMMENTS')
    comment.text=commentText
def prettify(element):
    rough_string = ElementTree.tostring(element, encoding="UTF-8")#, method="xml")
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")
