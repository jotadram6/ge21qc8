import FWCore.ParameterSet.Config as cms

XMLIdealGeometryESSource = cms.ESSource("XMLIdealGeometryESSource",
    geomXMLFiles = cms.vstring(
    	  'Geometry/CMSCommonData/data/materials.xml',
        'Geometry/CMSCommonData/data/rotations.xml',
        'Geometry/MuonCommonData/data/GEMQC8/cms.xml',
        'Geometry/MuonCommonData/data/gemf/TDR_BaseLine/gemf.xml',
        'Geometry/MuonCommonData/data/GEMQC8/gem11.xml',
        'Geometry/MuonCommonData/data/GEMQC8/muonNumbering.xml',
        'Geometry/MuonCommonData/data/GEMQC8/muonSens.xml',
        'Geometry/MuonCommonData/data/GEMQC8/muonProdCuts.xml',
        'Geometry/MuonCommonData/data/GEMQC8/GEMSpecsFilter.xml',   # Phase-2 Muon
        'Geometry/MuonCommonData/data/GEMQC8/GEMSpecs.xml',
        'Geometry/MuonCommonData/data/GEMQC8/gembox.xml'
        ),
    rootNodeName = cms.string('cms:OCMS')
)
