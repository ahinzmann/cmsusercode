import os

minMasses=[1500,1900,2400,2800,3300,3800,4300,5200,6000] # for mass bins 1.9, 2.4, 3.0, 3.6, 4.2, 4.8, 5.4, 6.0, 7.0
maxMasses=[1900,2400,2800,3300,3800,4300,5200,6000,13000] # for mass bins 1.9, 2.4, 3.0, 3.6, 4.2, 4.8, 5.4, 6.0, 7.0
lambdaTes=["",9000,10000,11000,12000,13000,14000,15000,16000,17000,18000,19000,20000,21000,22000,25000,30000]
couplings=[(0,0,0,1),]
samples=[]

for minMass in minMasses:
    for lambdaT in lambdaTes:
        for MD,nED,negInt,opMode in couplings:
            samples+=[('pythia8_add',minMass,maxMasses[minMasses.index(minMass)],lambdaT,MD,nED,negInt,opMode),]
            if lambdaT=="": break

#print samples

version="13TeV_CUETP8M1_Nov2022"

for sample,minMass,maxMass,lambdaT,MD,nED,negInt,opMode in samples:
  
  numjobs=30

  for jobnum in range(numjobs):

    samplename=sample+"_m"+str(minMass)+"_"+str(maxMass)+"_"+str(lambdaT)+"_"+str(MD)+"_"+str(nED)+"_"+str(negInt)+"_"+str(opMode)+"_"+version
    cfg=open(samplename+str(jobnum)+".py","w")
    cfg.writelines("""
import FWCore.ParameterSet.Config as cms

process = cms.Process("GEN")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(300000) )

process.load("Configuration.EventContent.EventContent_cff")
process.out = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    eventAutoFlushCompressedSize = cms.untracked.int32(20971520),
    fileName = cms.untracked.string('"""+samplename+str(jobnum)+""".root'),
    outputCommands = process.AODSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

#process.load("CommonTools.ParticleFlow.PF2PAT_EventContent_cff")
#process.out.outputCommands.extend( process.prunedAODForPF2PATEventContent.outputCommands )

# additional stuff for Maxime: 
process.out.outputCommands.extend(
    [
      'drop *',
      'keep GenEventInfoProduct_*_*_*',
      'keep *_ak4GenJets_*_*',
      #'keep *_ak4CaloJets_*_*',
      #'keep *_ak4JetID_*_*',
      #'keep *_ak4JetExtender_*_*',
      #'keep *_ak5GenJets_*_*',
      #'keep *_ak5CaloJets_*_*',
      #'keep *_ak5JetID_*_*',
      #'keep *_ak5JetExtender_*_*',
      #'keep *_ak7GenJets_*_*',
      #'keep *_ak7CaloJets_*_*',
      #'keep *_ak7JetID_*_*',
      #'keep *_ak7JetExtender_*_*',
      #------- PFJet collections --------
      #'keep *_kt6PFJets_rho_*',
      #'keep *_kt6PFJets_sigma_*',
      #'keep *_ak5PFJets_*_*',
      #'keep *_ak7PFJets_*_*',
      #------- Trigger collections ------
      #'keep edmTriggerResults_TriggerResults_*_*',
      #'keep *_hltTriggerSummaryAOD_*_*',
      #'keep L1GlobalTriggerObjectMapRecord_*_*_*',
      #'keep L1GlobalTriggerReadoutRecord_*_*_*',
      #------- Various collections ------
      #'keep *_EventAuxilary_*_*',
      #'keep *_offlinePrimaryVertices_*_*',
      #'keep *_offlinePrimaryVerticesWithBS_*_*',
      #------- MET collections ----------
      #'keep *_met_*_*',
      #'keep *_pfMet_*_*',
    ])

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Generator_cff')

process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic25ns13TeV2016Collision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('PRIVAT'),
    annotation = cms.untracked.string('PRIVAT'),
    name = cms.untracked.string('PRIVAT')
)

# Input source
process.source = cms.Source("EmptySource")
# Other statements

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mcRun2_asymptotic_v3', '')

process.MessageLogger=cms.Service("MessageLogger",
    destinations = cms.untracked.vstring('logs'),
    logs = cms.untracked.PSet(threshold=cms.untracked.string('WARNING'))
)

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *

process.generator = cms.EDFilter("Pythia8GeneratorFilter",
	comEnergy = cms.double(13000.0),
	crossSection = cms.untracked.double(1e10),
	filterEfficiency = cms.untracked.double(1),
	maxEventsToPrint = cms.untracked.int32(0),
	pythiaHepMCVerbosity = cms.untracked.bool(False),
	pythiaPylistVerbosity = cms.untracked.int32(0),
	
	PythiaParameters = cms.PSet(
                pythia8CommonSettingsBlock,
                pythia8CUEP8M1SettingsBlock,

		processParameters = cms.vstring(
			'PhaseSpace:mHatMin = """+str(minMass)+""" ',
			'PhaseSpace:mHatMax = """+str(maxMass)+""" ',
			'PhaseSpace:pTHatMin = """+str(minMass/10)+""" ',
""")
    if lambdaT=="" and "NonPert" in sample:
        cfg.writelines("""			'HardQCD:all = on ',
			'PartonLevel:MPI = off',
			'HadronLevel:Hadronize = off',
""")
    elif lambdaT=="":
        cfg.writelines("""			'HardQCD:all = on ',
""")
    else:
        cfg.writelines("""			'HardQCD:all = off ',
	                'ExtraDimensionsLED:dijets = on',
			'ExtraDimensionsLED:CutOffmode = 0',
			'ExtraDimensionsLED:LambdaT = """+str(lambdaT)+"""',
			'ExtraDimensionsLED:MD = """+str(MD)+"""',
			'ExtraDimensionsLED:n = """+str(nED)+"""', # Number of extra dimensions
			'ExtraDimensionsLED:nQuarkNew = 5', # outgoing mass-less quark flavours
			'ExtraDimensionsLED:negInt = """+str(negInt)+"""', # Change sign of interferecen term if ==1
			'ExtraDimensionsLED:opMode = """+str(opMode)+"""', # 0=Franceshini paper 1=Giudice paper
""")
    cfg.writelines("""
		),
		parameterSets = cms.vstring('pythia8CommonSettings',
                                            'pythia8CUEP8M1Settings',
                                            'processParameters')
	)
)

process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.Configuration.RecoGenJets_cff")
process.ak4GenJets.jetPtMin="""+str(minMass/10)+"""

process.RandomNumberGeneratorService.generator.initialSeed="""+str(1+jobnum)+"""

process.generation_step = cms.Path(process.generator*process.pgen)
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
process.p = cms.Path(process.genParticles*process.genJetParticles*process.ak4GenJets)#*process.ca08PrunedGenJets
process.endpath = cms.EndPath(process.out)
process.schedule = cms.Schedule(process.generation_step,process.p,process.endpath)
#process.out.outputCommands=cms.untracked.vstring('keep *','drop edmHepMCProduct_generator_*_*','drop *_genParticles*_*_*','drop *_genParticlesForJets*_*_*')
""")
    cfg.close()

    with open(samplename+str(jobnum)+".sh",'w+') as wrapper_script:
            wrapper_script.write("""#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd ~/rivet/CMSSW_10_6_16/src/SubstructureProfessor/Julian
cmsenv
cd /nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis
cmsRun """+samplename+str(jobnum)+""".py
#gfal-copy """+samplename+str(jobnum)+""".root "srm://dcache-se-cms.desy.de:8443/srm/managerv2?SFN=/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/dijet_angular/"
#rm """+samplename+str(jobnum)+""".root
""")
    with open(samplename+str(jobnum)+".submit",'w+') as htc_config:
            htc_config.write("""
#HTC Submission File for GEN sample production
#requirements      =  OpSysAndVer == "SL7"
universe          = vanilla
notification      = Error
notify_user       = andreas.hinzmann@desy.de
initialdir        = /nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/
output            = """+samplename+str(jobnum)+""".o
error             = """+samplename+str(jobnum)+""".e
log               = """+samplename+str(jobnum)+""".log
#Requesting CPU and DISK Memory - default +RequestRuntime of 3h stays unaltered
+RequestRuntime   = 40000
#RequestMemory     = 10G
JobBatchName      = """+samplename+"""
#RequestDisk       = 10G
getenv            = True
executable        = /usr/bin/sh
arguments         = " """+samplename+str(jobnum)+""".sh"
queue 1
""")
    string="condor_submit "+samplename+str(jobnum)+".submit"
    if jobnum%5!=4:
      string+=" &"
    print string
