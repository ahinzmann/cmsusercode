#import FWCore.ParameterSet.Config as cms

import copy

import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.ConfigToolBase import *

class MakeAODTrackCandidates(ConfigToolBase):
    """
       ------------------------------------------------------------------
       create selected tracks and a candidate hypothesis on AOD:
       
       process       : process
       label         : output collection will be <'patAOD'+label>
       tracks        : input tracks
       particleType  : particle type (for mass) 
       candSelection : preselection cut on the candidates
       ------------------------------------------------------------------    
       """
    _label='MakeAODTrackCandidates'
    def dumpPython(self):
        
        dumpPython = "\nfrom PhysicsTools.PatAlgos.tools.trackTools import *\n\nmakeAODTrackCandidates(process, "
        dumpPython += "'"+str(self.getvalue('label'))+"'"+", "
        dumpPython += str(self.getvalue('tracks'))+", "
        dumpPython += "'"+str(self.getvalue('particleType'))+"'"+", "
        dumpPython += "'"+str(self.getvalue('candSelection'))+"'"+')\n'
        return dumpPython
    
    def __call__(self,process,
                               label         = 'TrackCands',                 
                               tracks        = cms.InputTag('generalTracks'),
                               particleType  = "pi+",                        
                               candSelection = 'pt > 10'                     
                               ):
    
        self.addParameter('process',process, 'the process')
        self.addParameter('label',label, "output collection will be <'patAOD'+label>")
        self.addParameter('tracks',tracks, 'input tracks')
        self.addParameter('particleType',particleType, 'particle type (for mass)')
        self.addParameter('candSelection',candSelection, 'preselection cut on the candidates')
        
        process=self._parameters['process'].value
        tracks=self._parameters['tracks'].value
        particleType=self._parameters['particleType'].value
        candSelection=self._parameters['candSelection'].value
        process.disableRecording()
        process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi");
    ## add ChargedCandidateProducer from track
        setattr(process, 'patAOD' + label + 'Unfiltered', cms.EDProducer("ConcreteChargedCandidateProducer",
                                                                         src  = tracks,
                                                                         particleType = cms.string(particleType)
                                                                         )
                )
    ## add CandViewSelector with preselection string
        setattr(process, 'patAOD' + label, cms.EDFilter("CandViewSelector",
                                                        src = cms.InputTag('patAOD' + label + 'Unfiltered'),
                                                        cut = cms.string(candSelection)
                                                        )
                )
    ## run production of TrackCandidates at the very beginning of the sequence
        process.patDefaultSequence.replace(process.allLayer1Objects, getattr(process, 'patAOD' + label + 'Unfiltered') * getattr(process, 'patAOD' + label) * process.allLayer1Objects)
        process.enableRecording()
        process.addAction(copy.copy(self))

makeAODTrackCandidates=MakeAODTrackCandidates()

class MakePATTrackCandidates(ConfigToolBase):
    """
       ------------------------------------------------------------------
       create pat track candidates from AOD track collections:
       
       process       : process
       label         : output will be 'all/selectedLayer1'+label
       input         : name of the input collection
       selection     : selection on PAT Layer 1 objects
       isolation     : isolation to use (as 'source': value of dR)
       isoDeposits   : iso deposits
       mcAs          : replicate mc match as the one used by PAT
       on this AOD collection (None=no mc match);
       chosse 'photon', 'electron', 'muon', 'tau',
       'jet', 'met' as input string
       ------------------------------------------------------------------    
       """
    _label='MakePATTrackCandidates'
    def dumpPython(self):
        
        dumpPython = "\nfrom PhysicsTools.PatAlgos.tools.trackTools import *\n\nmakePATTrackCandidates(process, "
        dumpPython += "'"+str(self.getvalue('label'))+"'"+", "
        dumpPython += str(self.getvalue('input'))+", "
        dumpPython += "'"+str(self.getvalue('selection'))+"'"+", "
        dumpPython += str(self.getvalue('isolation'))+", "
        dumpPython += str(self.getvalue('isoDeposits'))+", "
        dumpPython += "'"+str(self.getvalue('mcAs'))+"'"+')\n'
        return dumpPython
    
    def __call__(self,process, 
                 label       = 'TrackCands',                    
                 input       = cms.InputTag('patAODTrackCands'),
                 selection   = 'pt > 10',                       
                 isolation   = {'tracker':0.3, 'ecalTowers':0.3, 'hcalTowers':0.3},  
                 isoDeposits = ['tracker','ecalTowers','hcalTowers'],   
                 mcAs        = 'muons'            
                 ):
       
        self.addParameter('process',process, '')
        self.addParameter('label',label, "output will be 'all/selectedLayer1'+label")
        self.addParameter('input',input, ' name of the input collection')
        self.addParameter('selection',selection, 'selection on PAT Layer 1 objects')
        self.addParameter('isolation',isolation, "isolation to use (as 'source': value of dR)")
        self.addParameter('isoDeposits',isoDeposits, 'iso deposits')
        self.addParameter('mcAs',mcAs, "replicate mc match as the one used by PAT on this AOD collection (None=no mc match); chosse 'photon', 'electron', 'muon', 'tau','jet', 'met' as input string")
        
        process=self._parameters['process'].value
        label=self._parameters['label'].value
        input=self._parameters['input'].value
        selection=self._parameters['selection'].value
        isolation=self._parameters['isolation'].value
        isoDeposits=self._parameters['isoDeposits'].value
        mcAs=self._parameters['mcAs'].value
        process.disableRecording()
        
    ## add allLayer1Tracks to the process
        from PhysicsTools.PatAlgos.producersLayer1.genericParticleProducer_cfi import allLayer1GenericParticles
        setattr(process, 'allLayer1' + label, allLayer1GenericParticles.clone(src = input))
    ## add selectedLayer1Tracks to the process
        setattr(process, 'selectedLayer1' + label, cms.EDFilter("PATGenericParticleSelector",
                                                                src = cms.InputTag("allLayer1"+label),
                                                                cut = cms.string(selection) 
                                                                ) 
                )
    ## add cleanLayer1Tracks to the process
        from PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi import cleanLayer1Tracks
        setattr(process, 'cleanLayer1' + label, cleanLayer1Tracks.clone(src = cms.InputTag('selectedLayer1' + label)))
        
    ## get them as variables, so we can put them in the sequences and/or configure them
        l1cands         = getattr(process, 'allLayer1' + label)
        selectedL1cands = getattr(process, 'selectedLayer1' + label)
        cleanL1cands    = getattr(process, 'cleanLayer1' + label)
        
    ## insert them in sequence, after the electrons
        process.allLayer1Objects.replace(process.allLayer1Electrons, l1cands + process.allLayer1Electrons)
        process.selectedLayer1Objects.replace(process.selectedLayer1Electrons, process.selectedLayer1Electrons + selectedL1cands)
        process.cleanLayer1Objects.replace(process.cleanLayer1Electrons, process.cleanLayer1Electrons + cleanL1cands)
        
    ## add them to the Summary Tables
        process.allLayer1Summary.candidates += [ cms.InputTag("allLayer1"+label) ]
        process.selectedLayer1Summary.candidates += [ cms.InputTag("selectedLayer1"+label) ]
        process.cleanLayer1Summary.candidates += [ cms.InputTag("cleanLayer1"+label) ]
        
    ## isolation: start with empty config
        if(isolation or isoDeposits):
            process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")
            process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")
            process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
        isoModules = []
        runIsoDeps = {'tracker':False, 'caloTowers':False}
            
        for source,deltaR in isolation.items():
        ## loop items in isolation
            if(source == 'tracker'):
                runIsoDeps['tracker'] = True
                l1cands.isolation.tracker = cms.PSet(
                    src    = cms.InputTag('pat'+label+'IsoDepositTracks'),
                    deltaR = cms.double(deltaR),
                    )
            elif(source == 'ecalTowers'):
                runIsoDeps['caloTowers'] = True
                l1cands.isolation.ecal = cms.PSet(
                src    = cms.InputTag('pat'+label+'IsoDepositCaloTowers', 'ecal'),
                deltaR = cms.double(deltaR),
                )
            elif(source == 'hcalTowers'):
                runIsoDeps['caloTowers'] = True
                l1cands.isolation.hcal = cms.PSet(
                    src    = cms.InputTag('pat'+label+'IsoDepositCaloTowers', 'hcal'),
                    deltaR = cms.double(deltaR),
                    )
            
        for source in isoDeposits:
        ## loop items in isoDeposits
            if(source == 'tracker'):
                runIsoDeps['tracker'] = True
                l1cands.isoDeposits.tracker = cms.InputTag('pat'+label+'IsoDepositTracks') 
            elif(source == 'ecalTowers'):
                runIsoDeps['caloTowers'] = True
                l1cands.isoDeposits.ecal = cms.InputTag('pat'+label+'IsoDepositCaloTowers', 'ecal') 
            elif(source == 'hcalTowers'):
                runIsoDeps['caloTowers'] = True
                l1cands.isoDeposits.hcal = cms.InputTag('pat'+label+'IsoDepositCaloTowers', 'hcal')
            
        for dep in [ dep for dep,runme in runIsoDeps.items() if runme == True ]:
            if(dep == 'tracker'):
                from RecoMuon.MuonIsolationProducers.trackExtractorBlocks_cff import MIsoTrackExtractorCtfBlock
                setattr(process, 'pat'+label+'IsoDepositTracks',
                        cms.EDProducer("CandIsoDepositProducer",
                                       src                  = input,
                                       trackType            = cms.string('best'),
                                       MultipleDepositsFlag = cms.bool(False),
                                       ExtractorPSet        = cms.PSet( MIsoTrackExtractorCtfBlock )
                                       )
                        )
                isoModules.append( getattr(process, 'pat'+label+'IsoDepositTracks') )
            elif(dep == 'caloTowers'):
                from RecoMuon.MuonIsolationProducers.caloExtractorByAssociatorBlocks_cff import MIsoCaloExtractorByAssociatorTowersBlock
                setattr(process, 'pat'+label+'IsoDepositCaloTowers',
                        cms.EDProducer("CandIsoDepositProducer",
                                       src                  = input,
                                       trackType            = cms.string('best'),
                                       MultipleDepositsFlag = cms.bool(True),
                                       ExtractorPSet        = cms.PSet( MIsoCaloExtractorByAssociatorTowersBlock )
                                       )
                        )
                isoModules.append( getattr(process, 'pat'+label+'IsoDepositCaloTowers') )
        for m in isoModules:
            process.patDefaultSequence.replace(l1cands, m * l1cands)
        
         # MC
        from PhysicsTools.PatAlgos.tools.helpers import MassSearchParamVisitor
        if(type(mcAs) != type(None)):
            findMatch= []
            findMatch.append(getattr(process, mcAs+'Match'))

        ## clone mc matchiong module of object mcAs and add it to the path
            setattr(process, 'pat'+label+'MCMatch', findMatch[0].clone(src = input))
            process.patDefaultSequence.replace( l1cands, getattr(process, 'pat'+label+'MCMatch') * l1cands)
            l1cands.addGenMatch = True
            l1cands.genParticleMatch = cms.InputTag('pat'+label+'MCMatch')
            process.enableRecording()
            process.addAction(copy.copy(self))
makePATTrackCandidates = MakePATTrackCandidates()

class MakeTrackCandidates(ConfigToolBase):
    """
       ------------------------------------------------------------------
       create selected tracks and a candidate hypothesis on AOD:
       
       process       : process
       label         : output collection will be <'patAOD'+label>
       tracks        : input tracks
       particleType  : particle type (for mass) 
       preselection  : preselection cut on the AOD candidates
       selection     : selection cut on the PAT candidates (for the
       selectedLayer1Candidate collection)
       isolation     : isolation to use (as 'source': value of dR)
       isoDeposits   : iso deposits
       mcAs          : replicate mc match as the one used by PAT
        on this AOD collection (None=no mc match);
        chosse 'photon', 'electron', 'muon', 'tau',
        'jet', 'met' as input string
        ------------------------------------------------------------------    
        """    
    
    _label='MakeTrackCandidates'
    def dumpPython(self):
        
        dumpPython = "\nfrom PhysicsTools.PatAlgos.tools.trackTools import *\n\nmakeTrackCandidates(process, "
        dumpPython += "'"+str(self.getvalue('label'))+"'"+", "
        dumpPython += str(self.getvalue('tracks'))+", "
        dumpPython += "'"+str(self.getvalue('particleType'))+"'"+", "
        dumpPython += "'"+str(self.getvalue('preselection'))+"'"+", "
        dumpPython += "'"+str(self.getvalue('selection'))+"'"+", "
        dumpPython += str(self.getvalue('isolation'))+", "
        dumpPython += str(self.getvalue('isoDeposits'))+", "
        dumpPython += "'"+str(self.getvalue('mcAs'))+"'"+')\n'
        return dumpPython
    
    def __call__(self,process, 
                 label        = 'TrackCands',                 
                 tracks       = cms.InputTag('generalTracks'),
                 particleType = 'pi+',                        
                 preselection = 'pt > 10',                    
                 selection    = 'pt > 10',                     
                 isolation    = {'tracker':   0.3,             
                                 'ecalTowers':0.3,             
                                 'hcalTowers':0.3              
                                 },
                 isoDeposits  = ['tracker','ecalTowers','hcalTowers'],
                 mcAs         = 'muons'
                 ) :

        self.addParameter('process',process, 'the process')
        self.addParameter('label',label, "output collection will be <'patAOD'+label>")
        self.addParameter('tracks',tracks, "input tracks")
        self.addParameter('particleType',particleType, "particle type (for mass)")
        self.addParameter('preselection',preselection, "preselection cut on the AOD candidates")
        self.addParameter('selection',selection, "selection cut on the PAT candidates (for the selectedLayer1Candidate collection)")
        self.addParameter('isolation',isolation, "isolation to use (as 'source': value of dR)")
        self.addParameter('isoDeposits',isoDeposits, " iso deposits")
        self.addParameter('mcAs',mcAs, "replicate mc match as the one used by PAT on this AOD collection (None=no mc match); chosse 'photon', 'electron', 'muon', 'tau', 'jet', 'met' as input string")
        
        process=self._parameters['process'].value
        label=self._parameters['label'].value
        tracks=self._parameters['tracks'].value
        particleType=self._parameters['particleType'].value
        preselection=self._parameters['preselection'].value
        selection=self._parameters['selection'].value
        isolation=self._parameters['isolation'].value
        isoDeposits=self._parameters['isoDeposits'].value
        mcAs=self._parameters['mcAs'].value
        process.disableRecording()
        """
        ------------------------------------------------------------------
        create selected tracks and a candidate hypothesis on AOD:
        
        process       : process
        label         : output collection will be <'patAOD'+label>
        tracks        : input tracks
        particleType  : particle type (for mass) 
        preselection  : preselection cut on the AOD candidates
        selection     : selection cut on the PAT candidates (for the
        selectedLayer1Candidate collection)
        isolation     : isolation to use (as 'source': value of dR)
        tracker     : as muon iso from tracks
        ecalTowers  : as muon iso from calo towers.
        hcalTowers  : as muon iso from calo towers.
        isoDeposits   : iso deposits
        mcAs          : replicate mc match as the one used by PAT
        on this AOD collection (None=no mc match);
        chosse 'photon', 'electron', 'muon', 'tau',
        'jet', 'met' as input string
        ------------------------------------------------------------------    
        """    
    
        makeAODTrackCandidates(process,
                               tracks        = tracks,
                               particleType  = particleType,
                               candSelection = preselection,
                               label         = label
                               ) 
        makePATTrackCandidates(process,
                               label         = label,
                               input         = cms.InputTag('patAOD' + label), 
                               isolation     = isolation,
                               isoDeposits   = isoDeposits,
                               mcAs          = mcAs,
                               selection     = selection
                               )
        process.enableRecording()
        process.addAction(copy.copy(self))
makeTrackCandidates=MakeTrackCandidates()
