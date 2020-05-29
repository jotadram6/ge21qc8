import sys
import os
import datetime
print datetime.datetime.now()
import socket
import FWCore.ParameterSet.Config as cms

sys.path.append(os.getenv('CMSSW_BASE')+'/src/Analysis/GEMQC8/test')

import configureRun_cfi as runConfig

run_number = runConfig.RunNumber

maxEvts = 25000

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

process = cms.Process("ISPY")

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

# Get certified events from file
pyhtonModulesPath = os.path.abspath("runGEMCosmicStand_validation.py").split('QC8Test')[0]+'QC8Test/src/Analysis/GEMQC8/python/'
sys.path.insert(1,pyhtonModulesPath)
from readCertEvtsFromFile import GetCertifiedEvents
certEvts = GetCertifiedEvents(run_number)

# Reconstruction of muon track
process.load('RecoMuon.TrackingTools.MuonServiceProxy_cff')
process.MuonServiceProxy.ServiceParameters.Propagators.append('StraightLinePropagator')

process.GEMCosmicMuonForQC8 = cms.EDProducer("GEMCosmicMuonForQC8",
                                             process.MuonServiceProxy,
                                             gemRecHitLabel = cms.InputTag("gemRecHits"),
                                             maxClusterSize = cms.double(runConfig.maxClusterSize),
                                             minClusterSize = cms.double(runConfig.minClusterSize),
                                             trackChi2 = cms.double(runConfig.trackChi2),
                                             trackResX = cms.double(runConfig.trackResX),
                                             trackResY = cms.double(runConfig.trackResY),
                                             MulSigmaOnWindow = cms.double(runConfig.MulSigmaOnWindow),
                                             minNumberOfRecHitsPerTrack = cms.uint32(runConfig.minRecHitsPerTrack),
                                             SuperChamberType = cms.vstring(SuperChType),
                                             SuperChamberSeedingLayers = cms.vdouble(SuperChSeedingLayers),
                                             tripEvents = cms.vstring(certEvts),
                                             MuonSmootherParameters = cms.PSet(PropagatorAlong = cms.string('SteppingHelixPropagatorAny'),
                                                                               PropagatorOpposite = cms.string('SteppingHelixPropagatorAny'),
                                                                               RescalingFactor = cms.double(5.0)
                                                                               ),
                                             )
process.GEMCosmicMuonForQC8.ServiceParameters.GEMLayers = cms.untracked.bool(True)
process.GEMCosmicMuonForQC8.ServiceParameters.CSCLayers = cms.untracked.bool(False)
process.GEMCosmicMuonForQC8.ServiceParameters.RPCLayers = cms.bool(False)

process.rawTOhits_step = cms.Path(process.muonGEMDigis+process.gemRecHits)
process.reconstruction_step = cms.Path(process.GEMCosmicMuonForQC8)

from FWCore.MessageLogger.MessageLogger_cfi import *

process.add_(
    cms.Service("ISpyService",
        outputFileName = cms.untracked.string('qc8.ig'),
        outputIg = cms.untracked.bool(True),
        outputMaxEvents = cms.untracked.int32(maxEvts), # These are the number of events per ig file
        debug = cms.untracked.bool(False)
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(maxEvts) # These are the number of events to cycle through in the input root file
)

process.load("ISpy.Analyzers.ISpyEvent_cfi")
#process.load('ISpy.Analyzers.ISpyGEMRecHit_cfi')
process.load('ISpy.Analyzers.ISpyTrack_cfi')
process.load('ISpy.Analyzers.ISpyTrackingRecHit_cfi')

#process.ISpyGEMRecHit.iSpyGEMRecHitTag = cms.InputTag("gemRecHits")
process.ISpyTrack.iSpyTrackTag = cms.InputTag("GEMCosmicMuonForQC8")
process.ISpyTrackingRecHit.iSpyTrackingRecHitTag = cms.InputTag("GEMCosmicMuonForQC8")

process.iSpy = cms.Path(process.ISpyEvent*
                        #process.ISpyGEMRecHit*
                        process.ISpyTrack*
                        process.ISpyTrackingRecHit
)

process.schedule = cms.Schedule(
    process.rawTOhits_step,
    process.reconstruction_step,
    process.iSpy
)
