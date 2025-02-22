#import FWCore.ParameterSet.Config as cms
import copy
from PhysicsTools.PatAlgos.tools.ConfigToolBase import *
from PhysicsTools.PatAlgos.patEventContent_cff import *

class SwitchOnTrigger(ConfigToolBase):
    """ Tool to switch on pat::Trigger information and to add the needed modules to the patDefaultSequence.
        As trigger information in PAT is quite disc space consuming trigger information is not added per default but left to the user.


    """
    _label='SwitchOnTrigger'
    def dumpPython(self):
        
        dumpPython = "\nfrom PhysicsTools.PatAlgos.tools.trigTools import *\n\nswitchOnTrigger(process) \n "
        return dumpPython
    
    
    def __call__(self, process ):
        """ Enables trigger information in PAT  """

        self.addParameter('process',process, 'the process')
            
        process=self._parameters['process'].value
        process.disableRecording()
    ## add trigger modules to path
        process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
        process.patDefaultSequence += process.patTriggerSequence

    ## configure pat trigger
        process.patTrigger.onlyStandAlone = False
        
    ## add trigger specific event content to PAT event content
        process.out.outputCommands += patTriggerEventContent
        for matchLabel in process.patTriggerEvent.patTriggerMatches:
            process.out.outputCommands += [ 'keep patTriggerObjectsedmAssociation_patTriggerEvent_' + matchLabel + '_*' ]
        process.enableRecording()
        process.addAction(copy.copy(self))

switchOnTrigger=SwitchOnTrigger()

class SwitchOnTriggerStandAlone(ConfigToolBase):
    """ Tool to switch on pat::Trigger information in standalone mode and to add the needed modules to the patDefaultSequence
    """
    
    _label='SwitchOnTriggerStandAlone'
    def dumpPython(self):
        
        dumpPython = "\nfrom PhysicsTools.PatAlgos.tools.trigTools import *\n\nswitchOnTriggerStandAlone(process) \n "
        return dumpPython

    def __call__(self, process ):
    ## add trigger modules to path

        self.addParameter('process',process, 'the process')
        
        process=self._parameters['process'].value
        process.disableRecording()
        process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
        process.patDefaultSequence += process.patTriggerSequence
        
    ## configure pat trigger
        process.patTrigger.onlyStandAlone = True
        process.patTriggerSequence.remove( process.patTriggerEvent )
        process.out.outputCommands += patTriggerStandAloneEventContent
        process.enableRecording()
        process.addAction(copy.copy(self))
               
switchOnTriggerStandAlone=SwitchOnTriggerStandAlone()

class SwitchOnTriggerAll(ConfigToolBase):
    """ Tool to switch on all pat::Trigger information (including standalone mode) and to add the needed modules to the patDefaultSequence.
    """
    _label='SwitchOnTriggerAll'
    def dumpPython(self):
        
        dumpPython = "\nfrom PhysicsTools.PatAlgos.tools.trigTools import *\n\nswitchOnTriggerAll(process) \n "
        return dumpPython

    def __call__( process ):

        self.addParameter('process',process, 'the  process')
        
        process=self._parameters['process'].value
        process.disableRecording()          
        switchOnTrigger( process )
        process.out.outputCommands += patTriggerStandAloneEventContent
        process.enableRecording()
        process.addAction(copy.copy(self))

switchOnTriggerAll=SwitchOnTriggerAll()

class SwitchOnTriggerMatchEmbedding(ConfigToolBase):
    """ Tool to switch on trigger information embedding and to add the needed modules to the patDefaultSequence
    """
    _label='SwitchOnTriggerMatchEmbedding'
    def dumpPython(self):
        
        dumpPython = "\nfrom PhysicsTools.PatAlgos.tools.trigTools import *\n\nswitchOnTriggerMatchEmbedding(process) \n "
        return dumpPython
    
    def __call__(self,process ):
        self.addParameter('process',process, 'the process')
        process=self._parameters['process'].value
        process.disableRecording()         
        process.patTriggerSequence += process.patTriggerMatchEmbedder
        process.out.outputCommands += patEventContentTriggerMatch
        process.enableRecording()
        process.addAction(copy.copy(self))
switchOnTriggerMatchEmbedding=SwitchOnTriggerMatchEmbedding()
