import sys
import os
import datetime
print datetime.datetime.now()
import socket
import FWCore.ParameterSet.Config as cms

sys.path.append(os.getenv('CMSSW_BASE')+'/src/Analysis/GEMQC8/test')

import configureRun_cfi as runConfig

run_number = runConfig.RunNumber

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

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.source = cms.Source("EmptySource")

from FWCore.MessageLogger.MessageLogger_cfi import *

process.add_(
    cms.Service("ISpyService",
    outputESFileName = cms.untracked.string('qc8-geometry.ig')
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.load("ISpy.Analyzers.ISpyMuonGeometry_cfi")

process.iSpy = cms.Path(process.ISpyMuonGeometry)
process.schedule = cms.Schedule(process.iSpy)
