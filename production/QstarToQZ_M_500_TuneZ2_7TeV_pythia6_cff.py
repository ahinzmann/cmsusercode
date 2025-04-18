import FWCore.ParameterSet.Config as cms

from Configuration.Generator.PythiaUEZ2Settings_cfi import *

generator = cms.EDFilter("Pythia6GeneratorFilter",
	pythiaHepMCVerbosity = cms.untracked.bool(False),
	maxEventsToPrint = cms.untracked.int32(0),
	pythiaPylistVerbosity = cms.untracked.int32(0),
	filterEfficiency = cms.untracked.double(1),
	comEnergy = cms.double(7000.0),
	crossSection = cms.untracked.double(1.870e+2),
	
	PythiaParameters = cms.PSet(
	        pythiaUESettingsBlock,
                processParameters = cms.vstring(
        'MSEL=0 ! (D=1) 0 to select full user control',
        'MSTP(6)=1 ! excited quarks',
        'MSUB(147)=1 ! qg->d*',
        'MSUB(148)=1 ! qg->u*',
        'PMAS(343,1)=500 ! mass of d*',
        'PMAS(344,1)=500 ! mass of u*',
        'RTCM(41)=500 ! Lambda = mass',
        'RTCM(43)=1.0 ! f',
        'RTCM(44)=1.0 ! fp',
        'RTCM(45)=1.0 ! fs',
        '4000001:ALLOFF',
        '4000001:ONIFMATCH 1 23',
        '4000002:ALLOFF',
        '4000002:ONIFMATCH 2 23',
        ),
		parameterSets = cms.vstring(
		        'pythiaUESettings',
			'processParameters')
	)
)

configurationMetadata = cms.untracked.PSet(
	version = cms.untracked.string('\$Revision: 1.2 $'),
	name = cms.untracked.string('\$Source: /cvs/CMSSW/UserCode/hinzmann/production/QstarToQZ_M_500_TuneZ2_7TeV_pythia6_cff.py,v $'),
	annotation = cms.untracked.string('Fall2011 sample with PYTHIA6: Qstar -> qW, TuneZ2')
)

ProductionFilterSequence = cms.Sequence(generator)
