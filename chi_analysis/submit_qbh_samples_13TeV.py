import os

minMass=1900
samples=[]

version="Feb2025CP5"

for signalMass in [4500,5000,5500,6000,6500,7000,7500,8000]:
  samples+=[('qbh_RS1_'+str(signalMass)+"_"+version,'LHEfile_'+str(signalMass)+'_1_unweighted'),]
for signalMass in [7500,8000,8500,9000,9500,10000,11000,12000]:
  samples+=[('qbh_ADD6_'+str(+signalMass)+"_"+version,'LHEfile_'+str(signalMass)+'_0_unweighted'),]

#print(samples)

print("voms-proxy-init --voms cms --out /afs/desy.de/user/h/hinzmann/run2023/myproxy.pem")

for samplename, filename in samples:
  
  offset=0
  numjobs=1

  for jobnum in range(offset,offset+numjobs):

    cfg=open("submit/"+samplename+str(jobnum)+".py","w")
    cfg.writelines("""
import os
    
import FWCore.ParameterSet.Config as cms

process = cms.Process("GEN")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(300000) )

process.load("Configuration.EventContent.EventContent_cff")
process.out = cms.OutputModule(
    "PoolOutputModule",
    process.AODSIMEventContent,
    fileName = cms.untracked.string('data/"""+samplename+str(jobnum)+""".root'),
    )

process.load("CommonTools.ParticleFlow.PF2PAT_EventContent_cff")
process.out.outputCommands.extend( process.prunedAODForPF2PATEventContent.outputCommands )

# additional stuff for Maxime: 
process.out.outputCommands.extend(
    [
      'keep GenEventInfoProduct_*_*_*',
      'keep LHEEventProduct_*_*_*',
      'keep *_ak4GenJets_*_*',
      'keep *_ak4CaloJets_*_*',
      'keep *_ak4JetID_*_*',
      'keep *_ak4JetExtender_*_*',
      #------- Trigger collections ------
      'keep edmTriggerResults_TriggerResults_*_*',
      'keep *_hltTriggerSummaryAOD_*_*',
      'keep L1GlobalTriggerObjectMapRecord_*_*_*',
      'keep L1GlobalTriggerReadoutRecord_*_*_*',
      #------- Various collections ------
      'keep *_EventAuxilary_*_*',
      'keep *_offlinePrimaryVertices_*_*',
      'keep *_offlinePrimaryVerticesWithBS_*_*',
    ])

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Generator_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('PRIVAT'),
    annotation = cms.untracked.string('PRIVAT'),
    name = cms.untracked.string('PRIVAT')
)

# Input source
process.source = cms.Source("LHESource",
	fileNames = cms.untracked.vstring('file:qbh/"""+filename+""".lhe'),
        skipEvents = cms.untracked.uint32("""+str(jobnum*300000)+"""),      # events to skip
)
# Other statements
#process.GlobalTag.globaltag = 'START42_V12::All'

process.MessageLogger=cms.Service("MessageLogger",
    destinations = cms.untracked.vstring('logs'),
    logs = cms.untracked.PSet(threshold=cms.untracked.string('INFO'))
)

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *

process.generator = cms.EDFilter("Pythia8HadronizerFilter",
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(13000.),
    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CP5Settings',
                                    )
    )
)

process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.Configuration.RecoGenJets_cff")
process.ak4GenJets.jetPtMin="""+str(minMass/10)+"""

#process.RandomNumberGeneratorService.generator.initialSeed="""+str(jobnum)+"""

process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic25ns13TeV2016Collision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mcRun2_asymptotic_v3', '')
process.generation_step = cms.Path(process.generator*process.pgen)
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
process.p = cms.Path(process.genParticles*process.genJetParticles*process.ak4GenJets)#*process.ca08PrunedGenJets
process.endpath = cms.EndPath(process.out)
process.schedule = cms.Schedule(process.generation_step,process.p,process.endpath)

""")
    cfg.close()
    
    with open("submit/"+samplename+str(jobnum)+".sh",'w+') as wrapper_script:
            wrapper_script.write("""#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
source /cvmfs/grid.desy.de/etc/profile.d/grid-ui-env.sh
#cd ~/rivet/CMSSW_10_6_16/src/SubstructureProfessor/Julian
cd /data/dust/user/hinzmann/jetmass/CMSSW_14_1_0_pre4/src/
cmsenv
cd /data/dust/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis
export X509_USER_PROXY=/afs/desy.de/user/h/hinzmann/run2023/myproxy.pem
cmsRun submit/"""+samplename+str(jobnum)+""".py
gfal-copy -f data/"""+samplename+str(jobnum)+""".root davs://dcache-cms-webdav.desy.de:2880/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/dijet_angular/Feb2025CP5/
rm data/"""+samplename+str(jobnum)+""".root
""")
    with open("submit/"+samplename+str(jobnum)+".submit",'w+') as htc_config:
            htc_config.write("""
#HTC Submission File for GEN sample production
#requirements      =  OpSysAndVer == "SL7"
universe          = vanilla
notification      = Error
notify_user       = andreas.hinzmann@desy.de
initialdir        = /data/dust/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/
output            = submit/"""+samplename+str(jobnum)+""".o
error             = submit/"""+samplename+str(jobnum)+""".e
#log               = submit/"""+samplename+str(jobnum)+""".log
#Requesting CPU and DISK Memory - default +RequestRuntime of 3h stays unaltered
+RequestRuntime   = 40000
RequestMemory     = 10G
JobBatchName      = """+samplename+"""
#RequestDisk       = 10G
getenv            = True
executable        = /usr/bin/sh
arguments         = " submit/"""+samplename+str(jobnum)+""".sh"
queue 1
""")
    
    string="condor_submit submit/"+samplename+str(jobnum)+".submit"
    if jobnum%5!=4:
      string+=" &"
    print(string)
