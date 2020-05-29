# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: SingleMuPt100_cfi -s GEN,SIM,DIGI,L1,DIGI2RAW,RAW2DIGI,L1Reco,RECO --conditions auto:run2_mc --magField 38T_PostLS1 --datatier GEN-SIM --geometry GEMCosmicStand --eventcontent FEVTDEBUGHLT --era phase2_muon -n 100 --fileout out_reco.root

import sys
import os
import datetime
print datetime.datetime.now()
import socket
import FWCore.ParameterSet.Config as cms
import configureRun_cfi as runConfig

# options
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing('analysis')

run_number = runConfig.RunNumber

options.register("runNum",run_number,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Run number")

options.register("eventsPerJob",-1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "The number of events (in each file)")

options.register('mps',
                 '',
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.int,
                 "List of MPs to process")

options.parseArguments()

# The superchambers in the 15 slots
SuperChType = runConfig.StandConfiguration

print(SuperChType)

# Define and find column type. Default is L. If it is found an S in a column, that column type becomes S.
colType = ['S','S','S']
for col in range(0,3):
	for row in range(0,5):
		if (SuperChType[col*5+row]=='L'):
			colType[col] = 'L'

print(colType)

# Calculation of SuperChSeedingLayers from SuperChType
SuperChSeedingLayers = []

for i in range (0,30):
    SuperChSeedingLayers.append(0)

for j in range (0,3):
    for i in range (5*j,5*(j+1)):
        if (SuperChType[i]!='0'):
            SuperChSeedingLayers[i*2]=1
            SuperChSeedingLayers[i*2+1]=3
            break
    for i in range (5*(j+1)-1,5*j-1,-1):
        if (SuperChType[i]!='0'):
            SuperChSeedingLayers[i*2]=4
            SuperChSeedingLayers[i*2+1]=2
            break

print(SuperChSeedingLayers)

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RECO',eras.phase2_muon)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Analysis.GEMQC8.GeometryGEMCosmicStand_cff')
process.load('Configuration.StandardSequences.MagneticField_0T_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# DEFINITION OF THE SUPERCHAMBERS INSIDE THE STAND
for i in xrange(len(SuperChType)):
    column_row = '_c%d_r%d' % ((i/5)+1, i%5+1)
    if SuperChType[i]=='L' : size = 'L'
    if SuperChType[i]=='S' : size = 'S'
    if SuperChType[i]!='0' :
        geomFile = 'Analysis/GEMQC8/data/GeometryFiles/gem11'+size+column_row+'.xml'
        print(geomFile)
        process.XMLIdealGeometryESSource.geomXMLFiles.append(geomFile)
        print('-> Appended')

# Config importation & settings
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.eventsPerJob))

# Data input folder and file type

if (socket.gethostname()=="gem904qc8dqm"):

    fpath =  "/data/bigdisk/GEM-Data-Taking/GE11_QC8/Cosmics/run{:06d}/".format(int(run_number))

    for x in os.listdir(fpath):
        if x.endswith("ls0001_index000000.raw"):
            dataFileExtension = ".raw"
            uFEDKit = True
            break
        elif x.endswith("chunk_000000.dat"):
            dataFileExtension = ".dat"
            uFEDKit = False
            break
    if (not dataFileExtension.endswith(".raw") and not dataFileExtension.endswith(".dat")):
        print "Check the data files... First file (at least) is missing!"

else:

    fpath =  "/eos/cms/store/group/dpg_gem/comm_gem/QC8_Commissioning/run{:06d}/".format(int(run_number))

    for x in os.listdir(fpath):
        if x.endswith("ls0001_allindex.raw"):
            dataFileExtension = ".raw"
            uFEDKit = True
            break
        elif x.endswith("chunk_000000.dat"):
            dataFileExtension = ".dat"
            uFEDKit = False
            break
    if (not dataFileExtension.endswith(".raw") and not dataFileExtension.endswith(".dat")):
        print "Check the data files... First file (at least) is missing!"

# Input source
process.source = cms.Source("GEMLocalModeDataSource",
                            fileNames = cms.untracked.vstring ([fpath+x for x in sorted(os.listdir(fpath)) if x.endswith(dataFileExtension)]),
                            skipEvents=cms.untracked.uint32(0),
                            fedId = cms.untracked.int32(888),  # which fedID to assign
                            hasFerolHeader = cms.untracked.bool(uFEDKit),
                            runNumber = cms.untracked.int32(run_number),
                            )

process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))

############## DB file #################
from CondCore.CondDB.CondDB_cfi import *
CondDB.DBParameters.authenticationPath = cms.untracked.string('/afs/cern.ch/cms/DB/conddb')

eMapFile = 'GEMeMap_'+colType[0]+colType[1]+colType[2]+'.db'

CondDB.connect = cms.string('sqlite_fip:Analysis/GEMQC8/data/EMapFiles/'+eMapFile)

process.GEMCabling = cms.ESSource("PoolDBESSource",
                                  CondDB,
                                  toGet = cms.VPSet(cms.PSet(record = cms.string('GEMeMapRcd'),
                                                             tag = cms.string('GEMeMap_QC8')
                                                             )
                                                    )
                                  )
####################################

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
# validation event filter
process.load('EventFilter.L1TRawToDigi.validationEventFilter_cfi')

process.load('EventFilter.L1TRawToDigi.tmtFilter_cfi')
process.tmtFilter.mpList = cms.untracked.vint32(options.mps)

# Output definition
strOutput = runConfig.OutputFileName

# Additional output definition
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

# raw to digi
process.load('EventFilter.GEMRawToDigi.muonGEMDigis_cfi')
process.muonGEMDigis.InputLabel = cms.InputTag("source","gemLocalModeDataSource")
process.muonGEMDigis.useDBEMap = True
process.muonGEMDigis.unPackStatusDigis = True

# Getting hot and dead strips files
hotStripsFile = "Analysis/GEMQC8/data/HotStripsTables/Mask_HotStrips_run" + str(run_number) + ".dat"
deadStripsFile = "Analysis/GEMQC8/data/DeadStripsTables/Mask_DeadStrips_run" + str(run_number) + ".dat"

# digi to reco
process.load('RecoLocalMuon.GEMRecHit.gemRecHits_cfi')

process.gemRecHits = cms.EDProducer("GEMRecHitProducer",
                                    recAlgoConfig = cms.PSet(),
                                    recAlgo = cms.string('GEMRecHitStandardAlgo'),
                                    gemDigiLabel = cms.InputTag("muonGEMDigis"),
                                    maskFile = cms.FileInPath(hotStripsFile),
                                    deadFile = cms.FileInPath(deadStripsFile),
                                    applyMasking = cms.bool(True)
                                    )

process.load('RecoMuon.TrackingTools.MuonServiceProxy_cff')

# Certify Events
process.CertifyEventsQC8 = cms.EDProducer('CertifyEventsQC8',
                                         process.MuonServiceProxy,
                                         verboseSimHit = cms.untracked.int32(1),
                                         recHitsInputLabel = cms.InputTag('gemRecHits'),
                                         gemDigiLabel = cms.InputTag("muonGEMDigis"),
                                         maxClusterSize = cms.double(runConfig.maxClusterSize),
                                         minClusterSize = cms.double(runConfig.minClusterSize),
                                         nBinGlobalZR = cms.untracked.vdouble(200,200,200,150,180,250),
                                         RangeGlobalZR = cms.untracked.vdouble(564,572,786,794,786,802,110,260,170,350,100,350)
                                         )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('certify_events_'+strOutput)
                                   )

# Path and EndPath definitions
process.rawTOhits_step = cms.Path(process.muonGEMDigis+process.gemRecHits)
process.certify_events_step = cms.Path(process.CertifyEventsQC8)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
process.schedule = cms.Schedule(process.rawTOhits_step,
                                process.certify_events_step,
                                process.endjob_step
                                )

# enable validation event filtering
process.rawTOhits_step.remove(process.validationEventFilter)
process.rawTOhits_step.remove(process.tmtFilter)
