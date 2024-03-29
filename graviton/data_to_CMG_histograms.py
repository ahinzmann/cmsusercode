import os

## import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *

### MASTER FLAGS  ######################################################################

from switches import *

# turn this on if you want to pick a relval in input (see below)
pickRelVal = False

# AK5 sequence with no cleaning is the default
# the other sequences can be turned off with the following flags.
#JOSE: no need to run these guys for what you are up to
runAK5LC = False
runAK7 = False
runCA8 = True

#COLIN: will need to include the event filters in tagging mode

#COLIN : reactivate HPS when bugs corrected
hpsTaus = False

#COLIN: the following leads to rare segmentation faults
doJetPileUpCorrection = True

#patTaus can now be saved even when running the CMG sequence.
doEmbedPFCandidatesInTaus = False

runCMG = True

#add the L2L3Residual corrections only for data
if runOnMC:#MC
    jetCorrections=['L1FastJet','L2Relative','L3Absolute']
else:#Data
    jetCorrections=['L1FastJet','L2Relative','L3Absolute','L2L3Residual']

# process.load("CommonTools.ParticleFlow.Sources.source_ZtoMus_DBS_cfi")
#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

sep_line = "-" * 50
print sep_line
print 'running the following PF2PAT+PAT sequences:'
print '\tAK5'
if runAK5LC: print '\tAK5LC'
if runAK7: print '\tAK7'
if runCA8: print '\tCA8'
print 'embedding in taus: ', doEmbedPFCandidatesInTaus
print 'HPS taus         : ', hpsTaus
print 'produce CMG tuple: ', runCMG
print sep_line


### SOURCE DEFINITION  ################################################################


# process.source.fileNames = cms.untracked.vstring(['/store/relval/CMSSW_4_2_5/RelValTTbar/GEN-SIM-RECO/START42_V12-v1/0113/1C538A2F-799E-E011-8A7E-0026189438BD.root'])

#process.source.fileNames = cms.untracked.vstring(['/store/cmst3/user/hinzmann/fastsim/RSGravitonToWW_kMpl_01_M_1000_Tune23_7TeV_herwiggpp_PFAOD.root'])
process.source.fileNames = cms.untracked.vstring(['/store/cmst3/user/mpierini/DiJetFASTSIM/RSGravitonWW/RSGravitonToWW_kMpl_01_M_1000_Tune23_7TeV_herwiggpp_cff_py_GEN_FASTSIM_HLT_PU.root'])
#process.source.fileNames = cms.untracked.vstring(['file:pythia6_gravitonWW_1000_PFAOD.root'])

# process.load("CMGTools.Common.sources.SingleMu.Run2011A_May10ReReco_v1.AOD.source_cff")
#process.load("CMGTools.Common.sources.HT.Run2011A_May10ReReco_v1.AOD.V2.source_cff")
# process.load("CMGTools.Common.sources.VBF_HToTauTau_M_115_7TeV_powheg_pythia6_tauola.Summer11_PU_S4_START42_V11_v1.AODSIM.V2.source_cff")

if pickRelVal:
    process.source = cms.Source(
        "PoolSource",
        fileNames = cms.untracked.vstring(
        pickRelValInputFiles( cmsswVersion  = 'CMSSW_4_3_0_pre2'
                              , relVal        = 'RelValZmumuJets_Pt_20_300PU1'
                              , globalTag     = 'MC_42_V9_PU_E7TeV_AVE_2_BX2808'
                              , numberOfFiles = -1
                              )
        )
        )


# print "WARNING!!!!!!!!!!!!!!!! remove the following line (see .cfg) before running on the batch!"
#process.source.fileNames = process.source.fileNames[:12]

print 'PF2PAT+PAT+CMG for files:'
print process.source.fileNames

### DEFINITION OF THE PF2PAT+PAT SEQUENCES #############################################

from CMGTools.Common.Tools.getGlobalTag import getGlobalTag
process.GlobalTag.globaltag = cms.string(getGlobalTag(runOnMC))

# load the PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.out.fileName = cms.untracked.string('/tmp/hinzmann/RSGravitonToWW_kMpl_01_M_1000_Tune23_7TeV_herwiggpp_patTuple_PF2PAT.root')

# Configure PAT to use PF2PAT instead of AOD sources
# this function will modify the PAT sequences. It is currently 
# not possible to run PF2PAT+PAT and standart PAT at the same time
from PhysicsTools.PatAlgos.tools.pfTools import *

# ---------------- Sequence AK5 ----------------------


process.eIdSequence = cms.Sequence()

# PF2PAT+PAT sequence 1:
# no lepton cleaning, AK5PFJets

postfixAK5 = "AK5"
jetAlgoAK5="AK5"

#COLIN : we will need to add the L2L3Residual when they become available! also check the other calls to the usePF2PAT function. 
usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgoAK5, runOnMC=runOnMC, postfix=postfixAK5,
          jetCorrections=('AK5PFchs', jetCorrections))


if doJetPileUpCorrection:
    from CommonTools.ParticleFlow.Tools.enablePileUpCorrection import enablePileUpCorrection
    enablePileUpCorrection( process, postfix=postfixAK5)

#configure the taus
from CMGTools.Common.PAT.tauTools import *
if doEmbedPFCandidatesInTaus:
    embedPFCandidatesInTaus( process, postfix=postfixAK5, enable=True )
if hpsTaus:
    adaptPFTaus(process,"hpsPFTau",postfix=postfixAK5)
    #  note that the following disables the tau cleaning in patJets
    adaptSelectedPFJetForHPSTau(process,jetSelection="pt()>15.0",postfix=postfixAK5)
    # currently (Sept 27,2011) there are three sets of tau isolation discriminators better to choose in CMG tuples.
    removeHPSTauIsolation(process,postfix=postfixAK5)

   
# curing a weird bug in PAT..
from CMGTools.Common.PAT.removePhotonMatching import removePhotonMatching
removePhotonMatching( process, postfixAK5 )

getattr(process,"pfNoMuon"+postfixAK5).enable = False 
getattr(process,"pfNoElectron"+postfixAK5).enable = False 
getattr(process,"pfNoTau"+postfixAK5).enable = False 
getattr(process,"pfNoJet"+postfixAK5).enable = True
getattr(process,"pfIsolatedMuons"+postfixAK5).isolationCut = 999999
getattr(process,"pfIsolatedElectrons"+postfixAK5).isolationCut = 999999

process.patJetCorrFactorsAK5.rho = cms.InputTag("kt6PFJets","rho")

# adding vbtf and cic electron IDs
from CMGTools.Common.PAT.addPATElectronID_cff import addPATElectronID
#addPATElectronID( process, postfixAK5 , runOnMC )

# insert the PFMET sifnificance calculation
from CMGTools.Common.PAT.addMETSignificance_cff import addMETSig
#addMETSig( process, postfixAK5 )

# remove non-jet objects
from PhysicsTools.PatAlgos.tools.coreTools import *
removeSpecificPATObjects(process, names=["Electrons","Muons","Photons","Taus"], outputModules=['out'], postfix=postfixAK5)

process.patDefaultSequenceAK5.remove(process.patPFTauIsolationAK5)
process.patDefaultSequenceAK5.remove(process.patShrinkingConePFTauDiscriminationAK5)
process.patDefaultSequenceAK5.remove(process.btaggingAODAK5)
process.patJetsAK5.addBTagInfo=False
process.patJetsAK5.addDiscriminators=False

# ---------------- Sequence AK5LC, lepton x-cleaning ---------------

# PF2PAT+PAT sequence 2:
# lepton cleaning, AK5PFJets. This sequence is a clone of the AK5 sequence defined previously.
# just modifying the x-cleaning parameters, and the isolation cut for x-cleaning

if runAK5LC:
  print 'cloning AK5 sequence to prepare AK5LC sequence...'

  from PhysicsTools.PatAlgos.tools.helpers import cloneProcessingSnippet
  postfixLC = 'LC'
  # just cloning the first sequence, and enabling lepton cleaning 
  cloneProcessingSnippet(process, getattr(process, 'patPF2PATSequence'+postfixAK5), postfixLC)

  postfixAK5LC = postfixAK5+postfixLC
  getattr(process,"pfNoMuon"+postfixAK5LC).enable = True
  getattr(process,"pfNoElectron"+postfixAK5LC).enable = True 
  getattr(process,"pfIsolatedMuons"+postfixAK5LC).isolationCut = 0.2
  getattr(process,"pfIsolatedElectrons"+postfixAK5LC).isolationCut = 0.2

  #COLIN : need to add the VBTF e and mu id

  # configure MET significance
  getattr(process,"PFMETSignificance"+postfixAK5LC).inputPATElectrons = cms.InputTag('patElectrons'+postfixAK5LC)   
  getattr(process,"PFMETSignificance"+postfixAK5LC).inputPATMuons = cms.InputTag('patMuons'+postfixAK5LC)


  print 'cloning AK5 sequence to prepare AK5LC sequence...Done'

# ---------------- Sequence AK7, no lepton x-cleaning ---------------

# PF2PAT+PAT sequence 3
# no lepton cleaning, AK7PFJets

if runAK7:
  postfixAK7 = "AK7"
  jetAlgoAK7="AK7"

  #COLIN : argh! AK7PFchs does not seem to exist yet...
  # Maxime should maybe contact the JEC group if he wants them 
  usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgoAK7, runOnMC=runOnMC, postfix=postfixAK7,
          jetCorrections=('AK7PF', jetCorrections))

  # if doJetPileUpCorrection:
  #    enablePileUpCorrection( process, postfix=postfixAK7)

  # no need for taus in AK7 sequence. could remove the whole tau sequence to gain time?
  # if hpsTaus:
  #    adaptPFTaus(process,"hpsPFTau",postfix=postfixAK7)

  # no top projection: 
  getattr(process,"pfNoMuon"+postfixAK7).enable = False 
  getattr(process,"pfNoElectron"+postfixAK7).enable = False 
  getattr(process,"pfNoTau"+postfixAK7).enable = False 
  getattr(process,"pfNoJet"+postfixAK7).enable = True

  removePhotonMatching( process, postfixAK7 )

  # addPATElectronID( process, postfixAK7 , runOnMC )

  # addMETSig(process,postfixAK7)

# ---------------- Sequence CA8, no lepton x-cleaning ---------------


# PF2PAT+PAT sequence 4
# no lepton cleaning, CA8PFJets

if runCA8:
  postfixCA8 = "CA8"
  jetAlgoCA8="CA8"

  # gen jet reconstruction modules
  #process.load("RecoJets.Configuration.GenJetParticles_cff")
  #process.load("RecoJets.JetProducers.ca4GenJets_cfi")
  #process.ca08GenJets = process.ca4GenJets.clone( rParam = 0.8)

  # pf jet reconstruction module
  process.load("RecoJets.JetProducers.ca4PFJets_cfi")
  process.ca08PFJets = process.ca4PFJets.clone( rParam = cms.double(0.8),
                                              doAreaFastjet = cms.bool(True),
                                              doRhoFastjet = cms.bool(True),
                                              Rho_EtaMax = cms.double(6.0),
                                              Ghost_EtaMax = cms.double(7.0),
					      src = cms.InputTag("pfNoPileUpAK5"),
                                              )

  addJetCollection(process,cms.InputTag('ca08PFJets'),
                  'CA8','PF',
                  doJTA = False,
                  doBTagging = False,
                  jetCorrLabel = ('AK5PFchs', cms.vstring(jetCorrections)),
                  doType1MET = False,
                  doL1Cleaning = False,
                  doL1Counters = False,
                  genJetCollection=cms.InputTag("ca08GenJets"),
                  doJetID = False,
                  jetIdLabel = "ca8"
                  )
  process.patJetsCA8PF.addDiscriminators = cms.bool(False)
  process.patJetsCA8PF.embedCaloTowers = cms.bool(False)
  process.patJetsCA8PF.embedGenJetMatch = cms.bool(False)
  process.patJetsCA8PF.embedGenPartonMatch = cms.bool(False)
  process.patJetsCA8PF.addGenPartonMatch = cms.bool(False)
  process.patJetsCA8PF.getJetMCFlavour = cms.bool(False)
  process.patJetsCA8PF.embedPFCandidates = cms.bool(False)
  process.patJetCorrFactorsCA8PF.rho = cms.InputTag("kt6PFJets","rho")

  process.patJetsCA8PF.embedGenJetMatch = cms.bool(False)
  process.patJetsCA8PF.addGenJetMatch = cms.bool(False)

  #process.recoMyCAJets = cms.Sequence(process.genJetParticles*process.ca08GenJets*process.ca08PFJets*process.patJetGenJetMatchCA8PF*process.kt6PFJets*process.patJetCorrFactorsCA8PF*process.patJetsCA8PF*process.selectedPatJetsCA8PF)
  process.recoMyCAJets = cms.Sequence(process.ca08PFJets*process.patJetCorrFactorsCA8PF*process.patJetsCA8PF*process.selectedPatJetsCA8PF)

  ######
  # Jet pruning
  #####

  from RecoJets.JetProducers.GenJetParameters_cfi import *
  from RecoJets.JetProducers.PFJetParameters_cfi import *
  from RecoJets.JetProducers.CaloJetParameters_cfi import *
  from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
  from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters

  #prune CA08 jets
  process.ca08PrunedPFJets = cms.EDProducer(
                              "SubJetProducer",
                              PFJetParameters.clone(
                                  src = cms.InputTag("pfNoPileUpAK5"),
                                  #src = cms.InputTag(''), #default input, as for the fat CA jet
                                  doAreaFastjet = cms.bool(True),
                                  doRhoFastjet = cms.bool(False)
                                  ),
                              AnomalousCellParameters,
                              SubJetParameters,
                              jetAlgorithm = cms.string("CambridgeAachen"),
                              rParam = cms.double(0.8),
                              jetCollInstanceName=cms.string("subjets")
                              )
  process.ca08PrunedPFJets.nSubjets = cms.int32(2)

  process.ca08PrunedGenJets =  cms.EDProducer(
                             "SubJetProducer",
                             GenJetParameters.clone(
                                 #src = cms.InputTag("genParticlesForJetsNoNu"), #default input
                                 doAreaFastjet = cms.bool(False),
                                 doRhoFastjet = cms.bool(False)
                                 ),
                             AnomalousCellParameters,
                             SubJetParameters,
                             jetAlgorithm = cms.string("CambridgeAachen"),
                             rParam = cms.double(0.8),
                             jetCollInstanceName=cms.string("subjets")
                             )

  process.genjetFilter = cms.EDFilter("CandViewCountFilter",
                                  src = cms.InputTag("ca08PrunedGenJets"),
                                  minNumber = cms.uint32(1)
                                  )

  #add pruned CA08 jets
  addJetCollection(process, 
                 cms.InputTag('ca08PrunedPFJets'),         
                 'CA8Pruned', 'PF',
                 doJTA=False,            # Run Jet-Track association & JetCharge
                 doBTagging=False,       # Run b-tagging
                 jetCorrLabel= ('AK5PFchs', cms.vstring(jetCorrections)),
                 doType1MET=False,
                 doL1Cleaning=False,
                 doL1Counters=False,
          #       genJetCollection = cms.InputTag("ca8GenJetsNoNu"), #why not the corresponding gen level ?
                 genJetCollection = cms.InputTag("ca08PrunedGenJets:subjets"),
                 #    genJetCollection = cms.InputTag("ca08GenJets"),
                 doJetID = False
                 )
  process.patJetsCA8PrunedPF.addDiscriminators = cms.bool(False)
  process.patJetsCA8PrunedPF.embedCaloTowers = cms.bool(False)
  process.patJetsCA8PrunedPF.embedGenJetMatch = cms.bool(False)
  process.patJetsCA8PrunedPF.embedGenPartonMatch = cms.bool(False)
  process.patJetsCA8PrunedPF.addGenPartonMatch = cms.bool(False)
  process.patJetsCA8PrunedPF.getJetMCFlavour = cms.bool(False)
  process.patJetsCA8PrunedPF.embedPFCandidates = cms.bool(False)
  process.patJetCorrFactorsCA8PrunedPF.rho = cms.InputTag("kt6PFJets","rho")

  process.patJetsCA8PrunedPF.embedGenJetMatch = cms.bool(False)
  process.patJetsCA8PrunedPF.addGenJetMatch = cms.bool(False)

  #process.recoMyPrunedCAJets = cms.Sequence(process.ca08PrunedGenJets*process.ca08PrunedPFJets*process.patJetGenJetMatchCA8PrunedPF*process.patJetCorrFactorsCA8PrunedPF*process.patJetsCA8PrunedPF*process.selectedPatJetsCA8PrunedPF)
  process.recoMyPrunedCAJets = cms.Sequence(process.ca08PrunedPFJets*process.patJetCorrFactorsCA8PrunedPF*process.patJetsCA8PrunedPF*process.selectedPatJetsCA8PrunedPF)

  setattr(process,"patPATSequence"+postfixCA8, cms.Sequence( process.recoMyCAJets + process.recoMyPrunedCAJets ))

  ### END JET PRUNING SECTION
  ###############

# ---------------- Common stuff ---------------

process.load('CMGTools.Common.gen_cff')


process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
process.patTrigger.processName = cms.string('*')

### PATH DEFINITION #############################################


# trigger information (no selection)

#process.p = cms.Path( process.patTriggerDefaultSequence )
process.p = cms.Path( )

# filters

process.noscraping = cms.EDFilter("FilterOutScraping",
                                applyfilter = cms.untracked.bool(True),
                                debugOn = cms.untracked.bool(False),
                                numtrack = cms.untracked.uint32(10),
                                thresh = cms.untracked.double(0.25)
                                )

process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')

process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(15), 
                                           maxd0 = cms.double(2) 
                                           )


if not runOnMC:
    from HLTrigger.HLTfilters.hltHighLevel_cfi import *
    trigger_process_name="HLT"
    process.triggerSelectionReference = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::"+trigger_process_name, HLTPaths = ["HLT_HT450_v*","HLT_HT500_v*","HLT_HT550_v*"], throw=False)
    process.triggerSelectionAnalysis1 = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::"+trigger_process_name, HLTPaths = ["HLT_HT600_v*","HLT_FatJetMass850_DR1p1_Deta2p0_v*"], throw=False)
    process.triggerSelectionAnalysis2 = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::"+trigger_process_name, HLTPaths = ["HLT_HT650_v*","HLT_FatJetMass850_DR1p1_Deta2p0_v*"], throw=False)
    process.triggerSelectionAnalysis3 = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::"+trigger_process_name, HLTPaths = ["HLT_HT700_v*","HLT_FatJetMass850_DR1p1_Deta2p0_v*"], throw=False)
    process.triggerSelectionAnalysis4 = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::"+trigger_process_name, HLTPaths = ["HLT_FatJetMass850_DR1p1_Deta2p0_v*"], throw=False)
    process.triggerSelectionAll = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::"+trigger_process_name, HLTPaths = ["HLT_HT450_v*","HLT_HT500_v*","HLT_HT550_v*","HLT_HT600_v*","HLT_HT650*","HLT_HT700_v*","HLT_FatJetMass850_DR1p1_Deta2p0_v*"], throw=False)
    process.p += process.triggerSelectionAll
    process.p += process.noscraping
    process.p += process.HBHENoiseFilter

# gen ---- 

if runOnMC:
    process.p += process.genSequence 

process.load('RecoJets.Configuration.RecoPFJets_cff')
process.kt6PFJets.doRhoFastjet=True
process.p += process.kt6PFJets

# PF2PAT+PAT ---

process.p += getattr(process,"patPF2PATSequence"+postfixAK5)

if runAK5LC:
    process.p += getattr(process,"patPF2PATSequence"+postfixAK5LC) 

if runAK7:
    process.p += getattr(process,"patPF2PATSequence"+postfixAK7) 

if runCA8:
    process.selectedPatJetsAK5.cut="pt()>30"
    process.ak5PFconstituents = cms.EDProducer("JetSlimmer",
                                     JetInputCollection=cms.untracked.InputTag("selectedPatJetsAK5"),
                                     isPFJet=cms.untracked.bool(True),
                                     Verbosity=cms.untracked.bool(False)
                                     )
    process.p+=process.ak5PFconstituents
    process.p += getattr(process,"patPATSequence"+postfixCA8) 

# event cleaning (in tagging mode, no event rejected)

process.load('CMGTools.Common.eventCleaning.eventCleaning_cff')

#process.p += process.eventCleaningSequence

 
# CMG ---

# Analysis Sequence definition
process.load('CMGTools.DiJetHighMass.pfJets_cff')
process.load("CMGTools.Common.skims.cmgPFJetSel_cfi")
process.load("CMGTools.Common.histograms.baseMETHistograms_cfi")

process.baseMETHistograms.histograms.metSumEt = cms.untracked.VPSet(
            cms.untracked.PSet( 
               var = cms.untracked.string('pt()/sumEt()'),
               nbins = cms.untracked.int32(100),
               low = cms.untracked.double(0),
               high = cms.untracked.double(1)
               )
            )

if runCMG:
    
    process.load('CMGTools.Common.analysis_cff')
    # running on PFAOD -> calo objects are not available.
    # we'll need to reactivate caloMET, though
    # process.p += process.analysisSequence

    from CMGTools.Common.Tools.visitorUtils import replacePostfix
    
    if runAK5LC:
        cloneProcessingSnippet(process, getattr(process, 'analysisSequence'), 'AK5LCCMG')
        replacePostfix(getattr(process,"analysisSequenceAK5LCCMG"),'AK5','AK5LC') 
        process.p += process.analysisSequenceAK5LCCMG
        
    if runAK7:
        cloneProcessingSnippet(process, getattr(process, 'jetSequence'), 'AK7CMG')
        replacePostfix(getattr(process,"jetSequenceAK7CMG"),'AK5','AK7') 
        process.p += process.jetSequenceAK7CMG
        process.cmgPFJetSelAK7CMG.cut="pt()>30"

    if runCA8:
        cloneProcessingSnippet(process, getattr(process, 'jetSequence'), 'CA8CMG')
        replacePostfix(getattr(process,"jetSequenceCA8CMG"),'AK5','CA8PF') 
        process.analysisSequence += process.jetSequenceCA8CMG
        cloneProcessingSnippet(process, getattr(process, 'jetSequence'), 'CA8PrunedCMG')
        replacePostfix(getattr(process,"jetSequenceCA8PrunedCMG"),'AK5','CA8PrunedPF') 
	process.cmgPFJetCA8PrunedCMG.cfg.addConstituents=False
	process.cmgPFJetCA8PrunedCMG.cfg.baseJetFactory.addSubjets=True
	process.cmgPFBaseJetCA8PrunedCMG.cfg.addSubjets=True
        process.analysisSequence += process.jetSequenceCA8PrunedCMG
        #cloneProcessingSnippet(process, getattr(process, 'jetCutSummarySequence'), 'CA8CMG')
        #replacePostfix(getattr(process,"jetCutSummarySequenceCA8CMG"),'AK5','CA8') 
        #process.analysisSequence += process.jetCutSummarySequenceCA8CMG
        #cloneProcessingSnippet(process, getattr(process, 'jetHistogramSequence'), 'CA8CMG')
        #replacePostfix(getattr(process,"jetHistogramSequenceCA8CMG"),'AK5','CA8') 
        #process.analysisSequence += process.jetHistogramSequenceCA8CMG

    from CMGTools.Common.Tools.tuneCMGSequences import * 
    tuneCMGSequences(process, postpostfix='CMG')

    #process.p += process.analysisSequence

# Skimming
process.load('CMGTools.DiJetHighMass.skims.selHighMass_cff')
#process.load("CMGTools.Common.runInfoAccounting.runInfoAccounting_cfi")

process.baseMETHistograms.inputCollection = cms.InputTag("cmgPFMET")

process.load("CMGTools.Common.factories.cmgFatJet_cfi")

process.cmgPFJetSel.cut="pt()>30 && abs(eta())<2.5"

process.cmgPFLeadJet.inputCollection = cms.InputTag("cmgPFJetSel")

process.cmgFatJet.cfg.inputCollection = cms.InputTag("cmgPFJetSel")

process.cmgPFDiJet.cuts.dijetKinematics.jetsPhaseSpace.jet1=cms.string('abs(leg1().eta())<2.5 && leg1.getSelection("cuts_looseJetId")') #
process.cmgPFDiJet.cuts.dijetKinematics.jetsPhaseSpace.jet2=cms.string('abs(leg2().eta())<2.5 && leg2.getSelection("cuts_looseJetId")') #
#process.cmgPFDiJet.cuts.dijetKinematics.tightDeltaEta=cms.string('abs(leg1.eta()-leg2.eta()) < 0.7')

process.cmgDiFatJet.cuts.fatDijetKinematics.jetsPhaseSpace.jet1=cms.string('abs(leg1().eta())<2.5')
process.cmgDiFatJet.cuts.fatDijetKinematics.jetsPhaseSpace.jet2=cms.string('abs(leg2().eta())<2.5')

dijetSelection='(mass>890) && getSelection("cuts_dijetKinematics_jetsPhaseSpace") && getSelection("cuts_dijetKinematics_tightDeltaEta")'
process.cmgPFTightDiJet.cut = dijetSelection
process.cmgPFTightDiJetFilter = process.filterHighMass.copy()
process.cmgPFTightDiJetFilter.src=cms.InputTag('cmgPFTightDiJet')

process.cmgPFFatTightDiJet.cut=dijetSelection.replace("dijetKinematics","fatDijetKinematics")
process.cmgPFFatTightDiJet.src=cms.InputTag("cmgDiFatJet")

process.cmgPFMediumDiJet = process.cmgPFTightDiJet.copy()
process.cmgPFMediumDiJet.cut = 'getSelection("cuts_dijetKinematics_mediumMass") && getSelection("cuts_dijetKinematics_jetsPhaseSpace") && getSelection("cuts_dijetKinematics_tightDeltaEta")'
process.cmgPFMediumDiJetFilter = process.cmgPFTightDiJetFilter.copy()
process.cmgPFMediumDiJetFilter.src=cms.InputTag('cmgPFMediumDiJet')

process.cmgPFFatMediumDiJet = process.cmgPFFatTightDiJet.copy()
process.cmgPFFatMediumDiJet.cut = 'getSelection("cuts_fatDijetKinematics_mediumMass") && getSelection("cuts_fatDijetKinematics_jetsPhaseSpace") && getSelection("cuts_fatDijetKinematics_tightDeltaEta")'

process.cmgPFDiJetHistograms.histograms.mass[0].nbins = 7001
process.cmgPFDiJetHistograms.histograms.mass[0].low = 0
process.cmgPFDiJetHistograms.histograms.mass[0].high = 7000
process.cmgPFDiJetHistograms.histograms.mass1[0].nbins = 2001
process.cmgPFDiJetHistograms.histograms.mass1[0].low = 0
process.cmgPFDiJetHistograms.histograms.mass1[0].high = 2000
process.cmgPFDiJetHistograms.histograms.mass2[0].nbins = 2001
process.cmgPFDiJetHistograms.histograms.mass2[0].low = 0
process.cmgPFDiJetHistograms.histograms.mass2[0].high = 2000
process.cmgPFDiJetHistograms.inputCollection=cms.InputTag("cmgPFTightDiJet")
process.cmgPFDiJetHistograms.histograms.deltaPhi[0].nbins = 50
process.cmgPFDiJetHistograms.histograms.deltaEta[0].nbins = 50
process.cmgPFDiJetHistograms.histograms.deltaEta[0].low = -5
process.cmgPFDiJetHistograms.histograms.deltaEta[0].high = 5
process.cmgPFDiJetHistograms.histograms.pt[0].high = 3500
process.cmgPFDiJetHistograms.histograms.pt1[0].high = 3500
process.cmgPFDiJetHistograms.histograms.pt2[0].high = 3500
process.cmgPFDiJetHistograms.histograms.Pt1_vs_Pt2[0].high = 3500
process.cmgPFDiJetHistograms.histograms.Pt1_vs_Pt2[1].high = 3500
process.cmgPFDiJetHistograms.histograms.Eta1_vs_Eta2[0].nbins = 50
process.cmgPFDiJetHistograms.histograms.Eta1_vs_Eta2[0].low = -5
process.cmgPFDiJetHistograms.histograms.Eta1_vs_Eta2[0].high = 5
process.cmgPFDiJetHistograms.histograms.Eta1_vs_Eta2[1].nbins = 50
process.cmgPFDiJetHistograms.histograms.Eta1_vs_Eta2[1].low = -5
process.cmgPFDiJetHistograms.histograms.Eta1_vs_Eta2[1].high = 5
process.cmgPFDiJetHistograms.histograms.Eta1_vs_Phi1[0].nbins = 50
process.cmgPFDiJetHistograms.histograms.Eta1_vs_Phi1[1].nbins = 50
process.cmgPFDiJetHistograms.histograms.Eta1_vs_Phi1[0].low = -5
process.cmgPFDiJetHistograms.histograms.Eta1_vs_Phi1[0].high = 5
process.cmgPFDiJetHistograms.histograms.Eta2_vs_Phi2[0].nbins = 50
process.cmgPFDiJetHistograms.histograms.Eta2_vs_Phi2[1].nbins = 50
process.cmgPFDiJetHistograms.histograms.Eta2_vs_Phi2[0].low = -5
process.cmgPFDiJetHistograms.histograms.Eta2_vs_Phi2[0].high = 5
process.cmgPFDiJetHistograms.histograms.Phi1_vs_Phi2[0].nbins = 50
process.cmgPFDiJetHistograms.histograms.Phi1_vs_Phi2[1].nbins = 50
process.cmgPFDiJetHistograms.histograms.metSumPt = cms.untracked.VPSet(
        cms.untracked.PSet(
            var = cms.untracked.string('sqrt(pow(leg1().px() + leg2().px(),2)+pow(leg1().py() + leg2().py(),2)) / (leg1().pt() + leg2().pt())'),
            nbins = cms.untracked.int32(200),
            low = cms.untracked.double(0),
            high = cms.untracked.double(1),
            title = cms.untracked.string("|#sum #vec{p}_{T}| / #sum p_{T}")
            )
        )

process.dijetPrinter = cms.EDAnalyzer(
    "PhysicsObjectPrinter",
    inputCollection = cms.untracked.InputTag("cmgPFTightDiJet"),
    printSelections = cms.untracked.bool(False),
    reportEvery = cms.untracked.uint32(1)
)
process.dijetPrinter1tag=process.dijetPrinter.clone()
process.dijetPrinter2tag=process.dijetPrinter.clone()

process.dijetPrinterCA8Pruned=process.dijetPrinter.clone()
process.dijetPrinterCA8Pruned.inputCollection=cms.untracked.InputTag("cmgPFDiJetCA8Pruned")
process.dijetPrinterCA8Pruned1tag=process.dijetPrinterCA8Pruned.clone()
process.dijetPrinterCA8Pruned2tag=process.dijetPrinterCA8Pruned.clone()

process.cmgPFFatDiJetHistograms.histograms=process.cmgPFDiJetHistograms.histograms.clone()
process.cmgPFFatDiJetHistograms.inputCollection=cms.InputTag("cmgPFFatTightDiJet")
process.cmgPFFatDiJetHistograms.histograms.dptrel1 = cms.untracked.VPSet(
        cms.untracked.PSet(
            var = cms.untracked.string('(leg1().leadPtr().pt() - leg1().pt()) / leg1().pt()'),
            nbins = cms.untracked.int32(200),
            low = cms.untracked.double(-1),
            high = cms.untracked.double(1),
            title = cms.untracked.string("#Delta p_{T,1} / p_{T,1}")
            )
        )
process.cmgPFFatDiJetHistograms.histograms.dptrel2 = cms.untracked.VPSet(
        cms.untracked.PSet(
            var = cms.untracked.string('(leg2().leadPtr().pt() - leg2().pt()) / leg2().pt()'),
            nbins = cms.untracked.int32(200),
            low = cms.untracked.double(-1),
            high = cms.untracked.double(1),
            title = cms.untracked.string("#Delta p_{T,2} / p_{T,2}")
            )
        )
process.cmgPFFatDiJetHistograms.histograms.dpt1 = cms.untracked.VPSet(
        cms.untracked.PSet(
            var = cms.untracked.string('leg1().leadPtr().pt() - leg1().pt()'),
            nbins = cms.untracked.int32(200),
            low = cms.untracked.double(-500),
            high = cms.untracked.double(500),
            title = cms.untracked.string("#Delta p_{T,1}")
            )
        )
process.cmgPFFatDiJetHistograms.histograms.dpt2 = cms.untracked.VPSet(
        cms.untracked.PSet(
            var = cms.untracked.string('leg2().leadPtr().pt() - leg2().pt()'),
            nbins = cms.untracked.int32(200),
            low = cms.untracked.double(-500),
            high = cms.untracked.double(500),
            title = cms.untracked.string("#Delta p_{T,2}")
            )
        )
process.cmgPFFatDiJetHistograms.histograms.deta1 = cms.untracked.VPSet(
        cms.untracked.PSet(
            var = cms.untracked.string('leg1().leadPtr().eta() - leg1().eta()'),
            nbins = cms.untracked.int32(200),
            low = cms.untracked.double(-2),
            high = cms.untracked.double(2),
            title = cms.untracked.string("#Delta #eta_{1}")
            )
        )
process.cmgPFFatDiJetHistograms.histograms.deta2 = cms.untracked.VPSet(
        cms.untracked.PSet(
            var = cms.untracked.string('leg2().leadPtr().eta() - leg2().eta()'),
            nbins = cms.untracked.int32(200),
            low = cms.untracked.double(-2),
            high = cms.untracked.double(2),
            title = cms.untracked.string("#Delta #eta_{2}")
            )
        )
process.cmgPFFatDiJetHistograms.histograms.dphi1 = cms.untracked.VPSet(
        cms.untracked.PSet(
            var = cms.untracked.string('deltaPhi(leg1().leadPtr().phi() , leg1().phi())'),
            nbins = cms.untracked.int32(200),
            low = cms.untracked.double(0),
            high = cms.untracked.double(3.15),
            title = cms.untracked.string("#Delta #phi_{1}")
            )
        )
process.cmgPFFatDiJetHistograms.histograms.dphi2 = cms.untracked.VPSet(
        cms.untracked.PSet(
            var = cms.untracked.string('deltaPhi(leg2().leadPtr().phi() , leg2().phi())'),
            nbins = cms.untracked.int32(200),
            low = cms.untracked.double(0),
            high = cms.untracked.double(3.15),
            title = cms.untracked.string("#Delta #phi_{2}")
            )
        )
process.cmgPFFatDiJetHistograms.histograms.njets1 = cms.untracked.VPSet(
        cms.untracked.PSet(
            var = cms.untracked.string('leg1().numConstituents()'),
            nbins = cms.untracked.int32(20),
            low = cms.untracked.double(0),
            high = cms.untracked.double(20),
            title = cms.untracked.string("n_{1}")
            )
        )
process.cmgPFFatDiJetHistograms.histograms.njets2 = cms.untracked.VPSet(
        cms.untracked.PSet(
            var = cms.untracked.string('leg2().numConstituents()'),
            nbins = cms.untracked.int32(20),
            low = cms.untracked.double(0),
            high = cms.untracked.double(20),
            title = cms.untracked.string("n_{2}")
            )
        )

process.cmgPFJetSelCA8CMG.cut="pt()>30"
process.cmgPFJetSelCA8PrunedCMG.cut="pt()>30"

process.cmgPFLeadJetCA8Pruned=process.cmgPFLeadJet.copy()
process.cmgPFLeadJetCA8Pruned.inputCollection = cms.InputTag("cmgPFJetSelCA8PrunedCMG")

process.cmgPFDiJetCA8Pruned=process.cmgPFDiJet.copy()
process.cmgPFDiJetCA8Pruned.cfg.leg1Collection = cms.InputTag("cmgPFLeadJetCA8Pruned")
process.cmgPFDiJetCA8Pruned.cfg.leg2Collection = cms.InputTag("cmgPFLeadJetCA8Pruned")

#process.cmgPFDiJetCA8Pruned.cuts.w1tag=cms.string("leg1.mass()>70.")
#process.cmgPFDiJetCA8Pruned.cuts.w2tag=cms.string("leg2.mass()>70.")
#process.cmgPFDiJetCA8Pruned.cuts.w1tag=cms.string("abs(leg1.mass()-85.)<10.")
#process.cmgPFDiJetCA8Pruned.cuts.w2tag=cms.string("abs(leg2.mass()-85.)<10.")
#process.cmgPFDiJetCA8Pruned.cuts.w1tag=cms.string("abs(leg1.mass()-85.)<25.")
#process.cmgPFDiJetCA8Pruned.cuts.w2tag=cms.string("abs(leg2.mass()-85.)<25.")
process.cmgPFDiJetCA8Pruned.cuts.w1mtag=cms.string("abs(leg1.mass()-85.)<15.")
process.cmgPFDiJetCA8Pruned.cuts.w2mtag=cms.string("abs(leg2.mass()-85.)<15.")
process.cmgPFDiJetCA8Pruned.cuts.w1tag=cms.string("abs(leg1.mass()-85.)<15. && max(leg1.subjets()[0].mass(),leg1.subjets()[1].mass())/leg1.mass()<0.25")
process.cmgPFDiJetCA8Pruned.cuts.w2tag=cms.string("abs(leg2.mass()-85.)<15. && max(leg2.subjets()[0].mass(),leg2.subjets()[1].mass())/leg2.mass()<0.25")
process.cmgPFDiJetCA8Pruned.cuts.w1ltag=cms.string("abs(leg1.mass()-85.)<15. && max(leg1.subjets()[0].mass(),leg1.subjets()[1].mass())/leg1.mass()<0.30")
process.cmgPFDiJetCA8Pruned.cuts.w2ltag=cms.string("abs(leg2.mass()-85.)<15. && max(leg2.subjets()[0].mass(),leg2.subjets()[1].mass())/leg2.mass()<0.30")
process.cmgPFDiJetCA8Pruned.cuts.w1ttag=cms.string("abs(leg1.mass()-85.)<15. && max(leg1.subjets()[0].mass(),leg1.subjets()[1].mass())/leg1.mass()<0.20")
process.cmgPFDiJetCA8Pruned.cuts.w2ttag=cms.string("abs(leg2.mass()-85.)<15. && max(leg2.subjets()[0].mass(),leg2.subjets()[1].mass())/leg2.mass()<0.20")

process.cmgPFDiJetCA8PrunedHistograms=process.cmgPFDiJetHistograms.copy()
process.cmgPFDiJetCA8PrunedHistograms.inputCollection=cms.InputTag("cmgPFDiJetCA8Pruned")
process.cmgPFDiJetCA8PrunedHistograms.histograms.massdrop1 = cms.untracked.VPSet(
        cms.untracked.PSet(
            var = cms.untracked.string('max(leg1().subjets()[0].mass(),leg1().subjets()[1].mass())/leg1().mass()'),
            nbins = cms.untracked.int32(100),
            low = cms.untracked.double(0),
            high = cms.untracked.double(1),
            title = cms.untracked.string("massdrop_{1}")
            )
        )
process.cmgPFDiJetCA8PrunedHistograms.histograms.massdrop2 = cms.untracked.VPSet(
        cms.untracked.PSet(
            var = cms.untracked.string('max(leg2().subjets()[0].mass(),leg2().subjets()[1].mass())/leg2().mass()'),
            nbins = cms.untracked.int32(100),
            low = cms.untracked.double(0),
            high = cms.untracked.double(1),
            title = cms.untracked.string("massdrop_{2}")
            )
        )
process.cmgPFDiJetCA8PrunedHistograms.histograms.dR1 = cms.untracked.VPSet(
        cms.untracked.PSet(
            var = cms.untracked.string('sqrt(pow(deltaPhi(leg1().subjets()[0].phi() , leg1().subjets()[1].phi()),2)+pow(leg1().subjets()[0].eta() - leg1().subjets()[1].eta(),2))'),
            nbins = cms.untracked.int32(100),
            low = cms.untracked.double(0),
            high = cms.untracked.double(2),
            title = cms.untracked.string("#Delta R_{1}")
            )
        )
process.cmgPFDiJetCA8PrunedHistograms.histograms.dR2 = cms.untracked.VPSet(
        cms.untracked.PSet(
            var = cms.untracked.string('sqrt(pow(deltaPhi(leg2().subjets()[0].phi() , leg2().subjets()[1].phi()),2)+pow(leg2().subjets()[0].eta() - leg2().subjets()[1].eta(),2))'),
            nbins = cms.untracked.int32(100),
            low = cms.untracked.double(0),
            high = cms.untracked.double(2),
            title = cms.untracked.string("#Delta R_{2}")
            )
        )

process.cmgPFDiJetAK5CA8Pruned=process.cmgPFDiJet.copy()
process.cmgPFDiJetAK5CA8Pruned.cfg.leg1Collection = cms.InputTag("cmgPFLeadJet")
process.cmgPFDiJetAK5CA8Pruned.cfg.leg2Collection = cms.InputTag("cmgPFLeadJetCA8Pruned")

process.cmgPFDiJetAK5CA8PrunedHistograms=process.cmgPFDiJetHistograms.clone()
process.cmgPFDiJetAK5CA8PrunedHistograms.inputCollection=cms.InputTag("cmgPFDiJetAK5CA8Pruned")
process.cmgPFDiJetAK5CA8PrunedHistograms.histograms = cms.untracked.PSet()
process.cmgPFDiJetAK5CA8PrunedHistograms.histograms.dptrel = cms.untracked.VPSet(
        cms.untracked.PSet(
            var = cms.untracked.string('(leg2().pt() - leg1().pt()) / leg1().pt()'),
            nbins = cms.untracked.int32(200),
            low = cms.untracked.double(-1),
            high = cms.untracked.double(1),
            title = cms.untracked.string("#Delta p_{T} / p_{T}")
            )
        )
process.cmgPFDiJetAK5CA8PrunedHistograms.histograms.dpt = cms.untracked.VPSet(
        cms.untracked.PSet(
            var = cms.untracked.string('leg2().pt() - leg1().pt()'),
            nbins = cms.untracked.int32(200),
            low = cms.untracked.double(-500),
            high = cms.untracked.double(500),
            title = cms.untracked.string("#Delta p_{T}")
            )
        )
process.cmgPFDiJetAK5CA8PrunedHistograms.histograms.deta = cms.untracked.VPSet(
        cms.untracked.PSet(
            var = cms.untracked.string('leg2().eta() - leg1().eta()'),
            nbins = cms.untracked.int32(200),
            low = cms.untracked.double(-2),
            high = cms.untracked.double(2),
            title = cms.untracked.string("#Delta #eta")
            )
        )
process.cmgPFDiJetAK5CA8PrunedHistograms.histograms.dphi = cms.untracked.VPSet(
        cms.untracked.PSet(
            var = cms.untracked.string('deltaPhi(leg2().phi() , leg1().phi())'),
            nbins = cms.untracked.int32(200),
            low = cms.untracked.double(0),
            high = cms.untracked.double(3.15),
            title = cms.untracked.string("#Delta #phi")
            )
        )
process.cmgPFDiJetAK5CA8PrunedHistograms.histograms.dR = cms.untracked.VPSet(
        cms.untracked.PSet(
            var = cms.untracked.string('sqrt(pow(deltaPhi(leg2().phi() , leg1().phi()),2)+pow(leg2().eta() - leg1().eta(),2))'),
            nbins = cms.untracked.int32(200),
            low = cms.untracked.double(0),
            high = cms.untracked.double(4),
            title = cms.untracked.string("#Delta R")
            )
        )


process.w1mtag=process.highMass.copy()
process.w1mtag.src=cms.InputTag("cmgPFDiJetCA8Pruned")
process.w1mtag.cut='getSelection("cuts_w1mtag") || getSelection("cuts_w2mtag")'
process.w1mtagFilter = process.filterHighMass.copy()
process.w1mtagFilter.src=cms.InputTag('w1mtag')

process.cmgPFDiJetHistograms1mtag=process.cmgPFDiJetHistograms.copy()
process.cmgPFFatDiJetHistograms1mtag=process.cmgPFFatDiJetHistograms.copy()
process.cmgPFDiJetCA8PrunedHistograms1mtag=process.cmgPFDiJetCA8PrunedHistograms.copy()
process.cmgPFDiJetAK5CA8PrunedHistograms1tag=process.cmgPFDiJetAK5CA8PrunedHistograms.clone()

process.w2mtag=process.highMass.copy()
process.w2mtag.src=cms.InputTag("cmgPFDiJetCA8Pruned")
process.w2mtag.cut='getSelection("cuts_w1mtag") && getSelection("cuts_w2mtag")'
process.w2mtagFilter = process.filterHighMass.copy()
process.w2mtagFilter.src=cms.InputTag('w2mtag')

process.cmgPFDiJetHistograms2mtag=process.cmgPFDiJetHistograms.copy()
process.cmgPFFatDiJetHistograms2mtag=process.cmgPFFatDiJetHistograms.copy()
process.cmgPFDiJetCA8PrunedHistograms2mtag=process.cmgPFDiJetCA8PrunedHistograms.copy()
process.cmgPFDiJetAK5CA8PrunedHistograms2tag=process.cmgPFDiJetAK5CA8PrunedHistograms.clone()

process.w1tag=process.highMass.copy()
process.w1tag.src=cms.InputTag("cmgPFDiJetCA8Pruned")
process.w1tag.cut='getSelection("cuts_w1tag") || getSelection("cuts_w2tag")'
process.w1tagFilter = process.filterHighMass.copy()
process.w1tagFilter.src=cms.InputTag('w1tag')

process.cmgPFDiJetHistograms1tag=process.cmgPFDiJetHistograms.copy()
process.cmgPFFatDiJetHistograms1tag=process.cmgPFFatDiJetHistograms.copy()
process.cmgPFDiJetCA8PrunedHistograms1tag=process.cmgPFDiJetCA8PrunedHistograms.copy()

process.w2tag=process.highMass.copy()
process.w2tag.src=cms.InputTag("cmgPFDiJetCA8Pruned")
process.w2tag.cut='getSelection("cuts_w1tag") && getSelection("cuts_w2tag")'
process.w2tagFilter = process.filterHighMass.copy()
process.w2tagFilter.src=cms.InputTag('w2tag')

process.cmgPFDiJetHistograms2tag=process.cmgPFDiJetHistograms.copy()
process.cmgPFFatDiJetHistograms2tag=process.cmgPFFatDiJetHistograms.copy()
process.cmgPFDiJetCA8PrunedHistograms2tag=process.cmgPFDiJetCA8PrunedHistograms.copy()

process.w1ltag=process.highMass.copy()
process.w1ltag.src=cms.InputTag("cmgPFDiJetCA8Pruned")
process.w1ltag.cut='getSelection("cuts_w1ltag") || getSelection("cuts_w2ltag")'
process.w1ltagFilter = process.filterHighMass.copy()
process.w1ltagFilter.src=cms.InputTag('w1ltag')

process.cmgPFDiJetHistograms1ltag=process.cmgPFDiJetHistograms.copy()
process.cmgPFFatDiJetHistograms1ltag=process.cmgPFFatDiJetHistograms.copy()
process.cmgPFDiJetCA8PrunedHistograms1ltag=process.cmgPFDiJetCA8PrunedHistograms.copy()

process.w2ltag=process.highMass.copy()
process.w2ltag.src=cms.InputTag("cmgPFDiJetCA8Pruned")
process.w2ltag.cut='getSelection("cuts_w1ltag") && getSelection("cuts_w2ltag")'
process.w2ltagFilter = process.filterHighMass.copy()
process.w2ltagFilter.src=cms.InputTag('w2ltag')

process.cmgPFDiJetHistograms2ltag=process.cmgPFDiJetHistograms.copy()
process.cmgPFFatDiJetHistograms2ltag=process.cmgPFFatDiJetHistograms.copy()
process.cmgPFDiJetCA8PrunedHistograms2ltag=process.cmgPFDiJetCA8PrunedHistograms.copy()

process.w1ttag=process.highMass.copy()
process.w1ttag.src=cms.InputTag("cmgPFDiJetCA8Pruned")
process.w1ttag.cut='getSelection("cuts_w1ttag") || getSelection("cuts_w2ttag")'
process.w1ttagFilter = process.filterHighMass.copy()
process.w1ttagFilter.src=cms.InputTag('w1ttag')

process.cmgPFDiJetHistograms1ttag=process.cmgPFDiJetHistograms.copy()
process.cmgPFFatDiJetHistograms1ttag=process.cmgPFFatDiJetHistograms.copy()
process.cmgPFDiJetCA8PrunedHistograms1ttag=process.cmgPFDiJetCA8PrunedHistograms.copy()

process.w2ttag=process.highMass.copy()
process.w2ttag.src=cms.InputTag("cmgPFDiJetCA8Pruned")
process.w2ttag.cut='getSelection("cuts_w1ttag") && getSelection("cuts_w2ttag")'
process.w2ttagFilter = process.filterHighMass.copy()
process.w2ttagFilter.src=cms.InputTag('w2ttag')

process.cmgPFDiJetHistograms2ttag=process.cmgPFDiJetHistograms.copy()
process.cmgPFFatDiJetHistograms2ttag=process.cmgPFFatDiJetHistograms.copy()
process.cmgPFDiJetCA8PrunedHistograms2ttag=process.cmgPFDiJetCA8PrunedHistograms.copy()

process.cmgPFDiJetHistogramsReference=process.cmgPFDiJetHistograms.copy()
process.cmgPFDiJetHistogramsReference.inputCollection=cms.InputTag("cmgPFMediumDiJet")
process.cmgPFFatDiJetHistogramsReference=process.cmgPFFatDiJetHistograms.copy()
process.cmgPFFatDiJetHistogramsReference.inputCollection=cms.InputTag("cmgPFFatMediumDiJet")
process.cmgPFDiJetCA8PrunedHistogramsReference=process.cmgPFDiJetCA8PrunedHistograms.copy()
process.cmgPFDiJetCA8PrunedHistogramsReference.inputCollection=cms.InputTag("cmgPFDiJetCA8Pruned")
process.cmgPFDiJetHistogramsReference1tag=process.cmgPFDiJetHistogramsReference.copy()
process.cmgPFFatDiJetHistogramsReference1tag=process.cmgPFFatDiJetHistogramsReference.copy()
process.cmgPFDiJetCA8PrunedHistogramsReference1tag=process.cmgPFDiJetCA8PrunedHistogramsReference.copy()
process.cmgPFDiJetHistogramsReference2tag=process.cmgPFDiJetHistogramsReference.copy()
process.cmgPFFatDiJetHistogramsReference2tag=process.cmgPFFatDiJetHistogramsReference.copy()
process.cmgPFDiJetCA8PrunedHistogramsReference2tag=process.cmgPFDiJetCA8PrunedHistogramsReference.copy()
process.cmgPFDiJetHistogramsReference1mtag=process.cmgPFDiJetHistogramsReference.copy()
process.cmgPFFatDiJetHistogramsReference1mtag=process.cmgPFFatDiJetHistogramsReference.copy()
process.cmgPFDiJetCA8PrunedHistogramsReference1mtag=process.cmgPFDiJetCA8PrunedHistogramsReference.copy()
process.cmgPFDiJetHistogramsReference2mtag=process.cmgPFDiJetHistogramsReference.copy()
process.cmgPFFatDiJetHistogramsReference2mtag=process.cmgPFFatDiJetHistogramsReference.copy()
process.cmgPFDiJetCA8PrunedHistogramsReference2mtag=process.cmgPFDiJetCA8PrunedHistogramsReference.copy()
process.cmgPFDiJetHistogramsReference1ltag=process.cmgPFDiJetHistogramsReference.copy()
process.cmgPFFatDiJetHistogramsReference1ltag=process.cmgPFFatDiJetHistogramsReference.copy()
process.cmgPFDiJetCA8PrunedHistogramsReference1ltag=process.cmgPFDiJetCA8PrunedHistogramsReference.copy()
process.cmgPFDiJetHistogramsReference2ltag=process.cmgPFDiJetHistogramsReference.copy()
process.cmgPFFatDiJetHistogramsReference2ltag=process.cmgPFFatDiJetHistogramsReference.copy()
process.cmgPFDiJetCA8PrunedHistogramsReference2ltag=process.cmgPFDiJetCA8PrunedHistogramsReference.copy()
process.cmgPFDiJetHistogramsReference1ttag=process.cmgPFDiJetHistogramsReference.copy()
process.cmgPFFatDiJetHistogramsReference1ttag=process.cmgPFFatDiJetHistogramsReference.copy()
process.cmgPFDiJetCA8PrunedHistogramsReference1ttag=process.cmgPFDiJetCA8PrunedHistogramsReference.copy()
process.cmgPFDiJetHistogramsReference2ttag=process.cmgPFDiJetHistogramsReference.copy()
process.cmgPFFatDiJetHistogramsReference2ttag=process.cmgPFFatDiJetHistogramsReference.copy()
process.cmgPFDiJetCA8PrunedHistogramsReference2ttag=process.cmgPFDiJetCA8PrunedHistogramsReference.copy()

process.cmgPFDiJetHistogramsAnalysis1=process.cmgPFDiJetHistogramsReference.copy()
process.cmgPFFatDiJetHistogramsAnalysis1=process.cmgPFFatDiJetHistogramsReference.copy()
process.cmgPFDiJetCA8PrunedHistogramsAnalysis1=process.cmgPFDiJetCA8PrunedHistogramsReference.copy()
process.cmgPFDiJetHistogramsAnalysis11tag=process.cmgPFDiJetHistogramsReference.copy()
process.cmgPFFatDiJetHistogramsAnalysis11tag=process.cmgPFFatDiJetHistogramsReference.copy()
process.cmgPFDiJetCA8PrunedHistogramsAnalysis11tag=process.cmgPFDiJetCA8PrunedHistogramsReference.copy()
process.cmgPFDiJetHistogramsAnalysis12tag=process.cmgPFDiJetHistogramsReference.copy()
process.cmgPFFatDiJetHistogramsAnalysis12tag=process.cmgPFFatDiJetHistogramsReference.copy()
process.cmgPFDiJetCA8PrunedHistogramsAnalysis12tag=process.cmgPFDiJetCA8PrunedHistogramsReference.copy()

process.cmgPFDiJetHistogramsAnalysis2=process.cmgPFDiJetHistogramsReference.copy()
process.cmgPFFatDiJetHistogramsAnalysis2=process.cmgPFFatDiJetHistogramsReference.copy()
process.cmgPFDiJetCA8PrunedHistogramsAnalysis2=process.cmgPFDiJetCA8PrunedHistogramsReference.copy()
process.cmgPFDiJetHistogramsAnalysis21tag=process.cmgPFDiJetHistogramsReference.copy()
process.cmgPFFatDiJetHistogramsAnalysis21tag=process.cmgPFFatDiJetHistogramsReference.copy()
process.cmgPFDiJetCA8PrunedHistogramsAnalysis21tag=process.cmgPFDiJetCA8PrunedHistogramsReference.copy()
process.cmgPFDiJetHistogramsAnalysis22tag=process.cmgPFDiJetHistogramsReference.copy()
process.cmgPFFatDiJetHistogramsAnalysis22tag=process.cmgPFFatDiJetHistogramsReference.copy()
process.cmgPFDiJetCA8PrunedHistogramsAnalysis22tag=process.cmgPFDiJetCA8PrunedHistogramsReference.copy()

process.cmgPFDiJetHistogramsAnalysis3=process.cmgPFDiJetHistogramsReference.copy()
process.cmgPFFatDiJetHistogramsAnalysis3=process.cmgPFFatDiJetHistogramsReference.copy()
process.cmgPFDiJetCA8PrunedHistogramsAnalysis3=process.cmgPFDiJetCA8PrunedHistogramsReference.copy()
process.cmgPFDiJetHistogramsAnalysis31tag=process.cmgPFDiJetHistogramsReference.copy()
process.cmgPFFatDiJetHistogramsAnalysis31tag=process.cmgPFFatDiJetHistogramsReference.copy()
process.cmgPFDiJetCA8PrunedHistogramsAnalysis31tag=process.cmgPFDiJetCA8PrunedHistogramsReference.copy()
process.cmgPFDiJetHistogramsAnalysis32tag=process.cmgPFDiJetHistogramsReference.copy()
process.cmgPFFatDiJetHistogramsAnalysis32tag=process.cmgPFFatDiJetHistogramsReference.copy()
process.cmgPFDiJetCA8PrunedHistogramsAnalysis32tag=process.cmgPFDiJetCA8PrunedHistogramsReference.copy()

process.cmgPFDiJetHistogramsAnalysis4=process.cmgPFDiJetHistogramsReference.copy()
process.cmgPFFatDiJetHistogramsAnalysis4=process.cmgPFFatDiJetHistogramsReference.copy()
process.cmgPFDiJetCA8PrunedHistogramsAnalysis4=process.cmgPFDiJetCA8PrunedHistogramsReference.copy()
process.cmgPFDiJetHistogramsAnalysis41tag=process.cmgPFDiJetHistogramsReference.copy()
process.cmgPFFatDiJetHistogramsAnalysis41tag=process.cmgPFFatDiJetHistogramsReference.copy()
process.cmgPFDiJetCA8PrunedHistogramsAnalysis41tag=process.cmgPFDiJetCA8PrunedHistogramsReference.copy()
process.cmgPFDiJetHistogramsAnalysis42tag=process.cmgPFDiJetHistogramsReference.copy()
process.cmgPFFatDiJetHistogramsAnalysis42tag=process.cmgPFFatDiJetHistogramsReference.copy()
process.cmgPFDiJetCA8PrunedHistogramsAnalysis42tag=process.cmgPFDiJetCA8PrunedHistogramsReference.copy()

# Total process
process.dijetanalysis_fromPAT = cms.Sequence(
    process.cmgPFMET + 
    process.cmgPFMETSel + 
    process.baseMETHistograms +

    process.cmgPFJet + 
    process.cmgPFJetSel +
    process.cmgPFLeadJet +
    process.cmgPFDiJet +
    process.cmgPFMediumDiJet +
    process.cmgPFTightDiJet +
    
    process.cmgPFMediumDiJetFilter +

    #process.cmgFatJet +
    #process.cmgDiFatJet +
    #process.cmgPFFatMediumDiJet +
    #process.cmgPFFatTightDiJet +

    process.cmgPFJetCA8CMG + 
    process.cmgPFJetSelCA8CMG + 

    process.cmgPFJetCA8PrunedCMG + 
    process.cmgPFJetSelCA8PrunedCMG + 
    process.cmgPFLeadJetCA8Pruned +
    process.cmgPFDiJetCA8Pruned
    #process.cmgPFDiJetAK5CA8Pruned
    )
    
process.dijetanalysis_fromCMG = cms.Sequence(
    process.baseMETHistograms +

    process.cmgPFLeadJet +
    process.cmgPFDiJet +
    process.cmgPFMediumDiJet +
    process.cmgPFTightDiJet +

    process.cmgPFMediumDiJetFilter +

    #process.cmgFatJet +
    #process.cmgDiFatJet +
    #process.cmgPFFatMediumDiJet +
    #process.cmgPFFatTightDiJet +

    process.cmgPFLeadJetCA8Pruned +
    process.cmgPFDiJetCA8Pruned
    #process.cmgPFDiJetAK5CA8Pruned
    )

process.dijetanalysis = cms.Sequence(
    process.primaryVertexFilter + 

    process.cmgPFTightDiJetFilter +

    process.cmgPFDiJetHistograms +
    #process.cmgPFFatDiJetHistograms +
    process.cmgPFDiJetCA8PrunedHistograms +
    #process.cmgPFDiJetAK5CA8PrunedHistograms +

    #process.w1mtag +
    #process.w1mtagFilter +
    #process.cmgPFDiJetHistograms1mtag +
    #process.cmgPFFatDiJetHistograms1mtag +
    #process.cmgPFDiJetCA8PrunedHistograms1mtag +

    #process.w1ltag +
    #process.w1ltagFilter +
    #process.cmgPFDiJetHistograms1ltag +
    #process.cmgPFFatDiJetHistograms1ltag +
    #process.cmgPFDiJetCA8PrunedHistograms1ltag +

    process.w1tag +
    process.w1tagFilter +
    process.cmgPFDiJetHistograms1tag +
    #process.cmgPFFatDiJetHistograms1tag +
    process.cmgPFDiJetCA8PrunedHistograms1tag +
    #process.cmgPFDiJetAK5CA8PrunedHistograms1tag +

    #process.w1ttag +
    #process.w1ttagFilter +
    #process.cmgPFDiJetHistograms1ttag +
    #process.cmgPFFatDiJetHistograms1ttag +
    #process.cmgPFDiJetCA8PrunedHistograms1ttag +

    #process.w2mtag +
    #process.w2mtagFilter +
    #process.cmgPFDiJetHistograms2mtag +
    #process.cmgPFFatDiJetHistograms2mtag +
    #process.cmgPFDiJetCA8PrunedHistograms2mtag +

    #process.w2ltag +
    #process.w2ltagFilter +
    #process.cmgPFDiJetHistograms2ltag +
    #process.cmgPFFatDiJetHistograms2ltag +
    #process.cmgPFDiJetCA8PrunedHistograms2ltag +

    process.w2tag +
    process.w2tagFilter +
    process.cmgPFDiJetHistograms2tag +
    #process.cmgPFFatDiJetHistograms2tag +
    process.cmgPFDiJetCA8PrunedHistograms2tag +
    #process.cmgPFDiJetAK5CA8PrunedHistograms2tag

    process.dijetPrinter2tag +
    process.dijetPrinterCA8Pruned2tag

    #process.w2ttag +
    #process.w2ttagFilter +
    #process.cmgPFDiJetHistograms2ttag +
    #process.cmgPFFatDiJetHistograms2ttag +
    #process.cmgPFDiJetCA8PrunedHistograms2ttag
    )

process.trigger_analysis_reference = cms.Sequence(
    process.primaryVertexFilter + 
    
    process.cmgPFDiJetHistogramsReference +
    #process.cmgPFFatDiJetHistogramsReference +
    process.cmgPFDiJetCA8PrunedHistogramsReference +

    #process.w1mtag +
    #process.w1mtagFilter +
    #process.cmgPFDiJetHistogramsReference1mtag +
    #process.cmgPFFatDiJetHistogramsReference1mtag +
    #process.cmgPFDiJetCA8PrunedHistogramsReference1mtag +

    #process.w1ltag +
    #process.w1ltagFilter +
    #process.cmgPFDiJetHistogramsReference1ltag +
    #process.cmgPFFatDiJetHistogramsReference1ltag +
    #process.cmgPFDiJetCA8PrunedHistogramsReference1ltag +

    process.w1tag +
    process.w1tagFilter +
    process.cmgPFDiJetHistogramsReference1tag +
    #process.cmgPFFatDiJetHistogramsReference1tag +
    process.cmgPFDiJetCA8PrunedHistogramsReference1tag +

    #process.w1ttag +
    #process.w1ttagFilter +
    #process.cmgPFDiJetHistogramsReference1ttag +
    #process.cmgPFFatDiJetHistogramsReference1ttag +
    #process.cmgPFDiJetCA8PrunedHistogramsReference1ttag +

    #process.w2mtag +
    #process.w2mtagFilter +
    #process.cmgPFDiJetHistogramsReference2mtag +
    #process.cmgPFFatDiJetHistogramsReference2mtag +
    #process.cmgPFDiJetCA8PrunedHistogramsReference2mtag +

    #process.w2ltag +
    #process.w2ltagFilter +
    #process.cmgPFDiJetHistogramsReference2ltag +
    #process.cmgPFFatDiJetHistogramsReference2ltag +
    #process.cmgPFDiJetCA8PrunedHistogramsReference2ltag +

    process.w2tag +
    process.w2tagFilter +
    process.cmgPFDiJetHistogramsReference2tag +
    #process.cmgPFFatDiJetHistogramsReference2tag +
    process.cmgPFDiJetCA8PrunedHistogramsReference2tag

    #process.w2ttag +
    #process.w2ttagFilter +
    #process.cmgPFDiJetHistogramsReference2ttag +
    #process.cmgPFFatDiJetHistogramsReference2ttag +
    #process.cmgPFDiJetCA8PrunedHistogramsReference2ttag
    )

process.trigger_analysis1 = cms.Sequence(
    process.primaryVertexFilter + 

    process.cmgPFDiJetHistogramsAnalysis1 +
    #process.cmgPFFatDiJetHistogramsAnalysis1 +
    process.cmgPFDiJetCA8PrunedHistogramsAnalysis1 +

    process.w1tag +
    process.w1tagFilter +
    process.cmgPFDiJetHistogramsAnalysis11tag +
    #process.cmgPFFatDiJetHistogramsAnalysis11tag +
    process.cmgPFDiJetCA8PrunedHistogramsAnalysis11tag +

    process.w2tag +
    process.w2tagFilter +
    process.cmgPFDiJetHistogramsAnalysis12tag +
    #process.cmgPFFatDiJetHistogramsAnalysis12tag +
    process.cmgPFDiJetCA8PrunedHistogramsAnalysis12tag
    )

process.trigger_analysis2 = cms.Sequence(
    process.primaryVertexFilter + 

    process.cmgPFDiJetHistogramsAnalysis2 +
    #process.cmgPFFatDiJetHistogramsAnalysis2 +
    process.cmgPFDiJetCA8PrunedHistogramsAnalysis2 +

    process.w1tag +
    process.w1tagFilter +
    process.cmgPFDiJetHistogramsAnalysis21tag +
    #process.cmgPFFatDiJetHistogramsAnalysis21tag +
    process.cmgPFDiJetCA8PrunedHistogramsAnalysis21tag +

    process.w2tag +
    process.w2tagFilter +
    process.cmgPFDiJetHistogramsAnalysis22tag +
    #process.cmgPFFatDiJetHistogramsAnalysis22tag +
    process.cmgPFDiJetCA8PrunedHistogramsAnalysis22tag
    )

process.trigger_analysis3 = cms.Sequence(
    process.primaryVertexFilter + 

    process.cmgPFDiJetHistogramsAnalysis3 +
    #process.cmgPFFatDiJetHistogramsAnalysis3 +
    process.cmgPFDiJetCA8PrunedHistogramsAnalysis3 +

    process.w1tag +
    process.w1tagFilter +
    process.cmgPFDiJetHistogramsAnalysis31tag +
    #process.cmgPFFatDiJetHistogramsAnalysis31tag +
    process.cmgPFDiJetCA8PrunedHistogramsAnalysis31tag +

    process.w2tag +
    process.w2tagFilter +
    process.cmgPFDiJetHistogramsAnalysis32tag +
    #process.cmgPFFatDiJetHistogramsAnalysis32tag +
    process.cmgPFDiJetCA8PrunedHistogramsAnalysis32tag
    )

process.trigger_analysis4 = cms.Sequence(
    process.primaryVertexFilter + 

    process.cmgPFDiJetHistogramsAnalysis4 +
    #process.cmgPFFatDiJetHistogramsAnalysis4 +
    process.cmgPFDiJetCA8PrunedHistogramsAnalysis4 +

    process.w1tag +
    process.w1tagFilter +
    process.cmgPFDiJetHistogramsAnalysis41tag +
    #process.cmgPFFatDiJetHistogramsAnalysis41tag +
    process.cmgPFDiJetCA8PrunedHistogramsAnalysis41tag +

    process.w2tag +
    process.w2tagFilter +
    process.cmgPFDiJetHistogramsAnalysis42tag +
    #process.cmgPFFatDiJetHistogramsAnalysis42tag +
    process.cmgPFDiJetCA8PrunedHistogramsAnalysis42tag
    )

process.p += process.dijetanalysis_fromPAT

### OUTPUT DEFINITION #############################################

# PF2PAT+PAT ---

# Add PF2PAT output to the created file
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning, patTriggerEventContent, patTriggerStandAloneEventContent
process.out.outputCommands = cms.untracked.vstring('drop *',
                                                   *patEventContentNoCleaning
                                                   )
# add trigger information to the pat-tuple
process.out.outputCommands += patTriggerEventContent
process.out.outputCommands += patTriggerStandAloneEventContent

# add gen event content to the pat-tuple (e.g. status 3 GenParticles)
from CMGTools.Common.eventContent.gen_cff import gen 
process.out.outputCommands.extend( gen )

# tuning the PAT event content to our needs
from CMGTools.Common.eventContent.patEventContentCMG_cff import patEventContentCMG
process.out.outputCommands.extend( patEventContentCMG )

# event cleaning results
from CMGTools.Common.eventContent.eventCleaning_cff import eventCleaning
process.out.outputCommands.extend( eventCleaning )

from CMGTools.Common.eventContent.runInfoAccounting_cff import runInfoAccounting
process.out.outputCommands.extend( runInfoAccounting )

if runCA8:
  process.out.outputCommands.extend(cms.untracked.vstring('keep *_ca08PrunedPFJets_*_*', 'keep *_ca08PrunedGenJets_*_*'))

# CMG ---

from CMGTools.Common.eventContent.everything_cff import everything 

process.outcmg = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('CMG_tree.root'),
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = everything,
    dropMetaData = cms.untracked.string('PRIOR')
    )


if runCMG:
    process.outpath += process.outcmg

if runCA8:
  process.outcmg.outputCommands.extend(cms.untracked.vstring('keep *_ak5PFconstituents_*_*'))
  process.outcmg.outputCommands.extend(cms.untracked.vstring('drop *_kt6PFJets_*_*'))
  process.outcmg.outputCommands.extend(cms.untracked.vstring('drop *_kt6PFJetsAK5_*_*'))
  process.outcmg.outputCommands.extend(cms.untracked.vstring('drop *_kt6PFJetsAK7_*_*'))
  process.outcmg.outputCommands.extend(cms.untracked.vstring('drop *_cmgPFBaseJetSel*_*_*'))
  process.outcmg.outputCommands.extend(cms.untracked.vstring('drop *_selectedPatElectronsAK5_*_*'))
  process.outcmg.outputCommands.extend(cms.untracked.vstring('drop *_selectedPatElectronsAK7_*_*'))

process.outcmg.outputCommands.extend(cms.untracked.vstring('keep *_TriggerResults_*_HLT'))
process.outcmg.outputCommands.extend(cms.untracked.vstring('keep *_HBHENoiseFilter_*_*'))

process.load("CMGTools.Common.runInfoAccounting_cff")
process.ria = cms.Sequence(
    process.runInfoAccountingDataSequence
    )
if runOnMC:
    process.ria = cms.Sequence(
        process.runInfoAccountingSequence
    )

process.outpath += process.ria
    

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("CMG_histograms.root"))

# process.Timing = cms.Service("Timing")

# print process.dumpPython()

if runOnMC:
    process.genParticlesStatus3.select=cms.vstring('keep status()==3','keep status()==2 && pdgId()==23','keep status()==2 && numberOfMothers()>0 && mother(0).pdgId()==23','keep status()==2 && pdgId()==24','keep status()==2 && numberOfMothers()>0 && mother(0).pdgId()==24')
    #process.outcmg.outputCommands.extend(cms.untracked.vstring('keep *_genParticles_*_*'))

process.outpath.remove(process.out)
del process.out

if runOnCMG:
    process.p = cms.Path(process.dijetanalysis_fromCMG)
    process.setName_("CMG")
    process.outpath.remove(process.outcmg)
    del process.outcmg
    process.outpath.remove(process.ria)
    del process.ria

    #samplename="428_HT_Run2011A-May10ReReco-v1_vv12"
    #samplename="428_HT_Run2011A-PromptReco-v4_vv12"
    #samplename="428_HT_Run2011A-05Aug2011-v1_vv12"
    #samplename="428_HT_Run2011A-PromptReco-v6_vv12"
    samplename="428_HT_Run2011B-PromptReco-v1_vv12"
    #samplename="QCD_vv5"
    #samplename="428_QCD_Pt-15to3000_Tune23_Flat_7TeV_herwigpp_Fall11-PU_S6_START42_V14B-v2_vv9_3"

    process.TFileService.fileName = cms.string("/tmp/hinzmann/"+samplename+"_metsumpt_histograms.root")

    files='CMG_tree.*root'
    from CMGTools.Production.datasetToSource import *
    process.source = datasetToSource('hinzmann', '/'+samplename, files)
    
    #process.source.fileNames=cms.untracked.vstring([])
    #filenames=os.listdir("/tmp/hinzmann/"+samplename)
    #process.source.fileNames+=["file:///tmp/hinzmann/"+samplename+"/"+filename for filename in filenames]
    
    process.source.eventsToProcess = cms.untracked.VEventRange('176309:1022524-176309:1022524','177074:1171400298-177074:1171400298')
    
    #print process.source.fileNames


if not runOnMC:
    process.p2 = process.p.copy()
    process.p2 += process.triggerSelectionAll
    process.p2 += process.dijetanalysis
    process.p3 = process.p.copy()
    process.p3 += process.triggerSelectionReference
    process.p3 += process.trigger_analysis_reference
    process.p4 = process.p.copy()
    process.p4 += process.triggerSelectionReference
    process.p4 += process.triggerSelectionAnalysis1
    process.p4 += process.trigger_analysis1
    #process.p5 = process.p.copy()
    #process.p5 += process.triggerSelectionReference
    #process.p5 += process.triggerSelectionAnalysis2
    #process.p5 += process.trigger_analysis2
    #process.p6 = process.p.copy()
    #process.p6 += process.triggerSelectionReference
    #process.p6 += process.triggerSelectionAnalysis3
    #process.p6 += process.trigger_analysis3
    process.p7 = process.p.copy()
    process.p7 += process.triggerSelectionReference
    process.p7 += process.triggerSelectionAnalysis4
    process.p7 += process.trigger_analysis4
else:
    process.p2 = process.p.copy()
    process.p2 += process.dijetanalysis
    process.p3 = process.p.copy()
    process.p3 += process.trigger_analysis_reference
