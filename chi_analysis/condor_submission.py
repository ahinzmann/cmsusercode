from glob import glob
import os,subprocess
from numpy import linspace,diff
from numpy.random import randint,seed
class GEN:
    __cfg_py = '/nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/MG_Pythia8_cfg.py'
    __gridpack_locations =['/nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis']

    def __init__(self,name,workdir,**kwargs):
        text_width = len(name)+8
        separator='#'
        print()
        print(''.center(text_width,separator))
        print(''.center(text_width,separator))
        print(' preparing environemt for: '.center(text_width,separator))
        print((' %s '%name).center(text_width,separator))
        print(''.center(text_width,separator))
        print(''.center(text_width,separator))
        print()

        self.__number_jobs = 1 if 'number_jobs' not in kwargs else kwargs['number_jobs']
        self.__number_events = 100 if 'events' not in kwargs else kwargs['events']
        self.__number_events_per_job =  (diff(linspace(0,self.__number_events,self.__number_jobs+1,dtype=int))).tolist()
        print( 'events per job:',self.__number_events_per_job)
        seed(27051993)
        self.__randomseeds = randint(10000000, size=self.__number_jobs).tolist()
        print('using the following random seeds:',self.__randomseeds)
        self.__already_exists = False
        self.workdir = workdir
        self.name = name
        self.gridpack = name if 'gridpack' not in kwargs else kwargs['gridpack']
        print('gridpack:',self.gridpack)
        if(self.__already_exists):
            print('submission env already exists.')
        else:
            self.write_wrapper_script()
            self.write_htc_config()
            self.write_cms_configs()
            with open(self.workdir+'/repr.info', 'w') as outfile:
                outfile.write(repr(self)+'\n')
                outfile.close()

            print('set up submission env for %i job(s)'%self.__number_jobs)
            
    def __repr__(self):
        return("GEN('"+self.name+"','"+self.workdir+"',gridpack='"+self.gridpack+"',events="+str(self.__number_events)+",number_jobs="+str(self.__number_jobs)+",recover=True)")
        
    @property
    def name(self):
        return self.__name

    @name.setter
    def name(self,name):
        if(name in glob(self.__workdir+'*')):
            raise ValueError("GEN process with name '"+name+"' already exists!")
        self.__name = name
        
    @property
    def workdir(self):
        return self.__workdir

    @workdir.setter
    def workdir(self,workdir):
        workdir = os.path.abspath(workdir)
        print('setting up workdir at:',workdir)
        if(os.path.exists(workdir)):
            self.__already_exists = True
        else:
            os.makedirs(workdir)
            for i in range(self.__number_jobs):
                os.makedirs(workdir+'/'+str(i))
        self.__workdir = workdir

    @property
    def gridpack(self):
        return self.__gridpack

    @gridpack.setter
    def gridpack(self,gridpack):
        if self.__already_exists:
            with open(self.workdir+'/repr.info','r') as repr_file:
                for l in repr_file:
                    self.__gridpack = l.split(',')[2].split('=')[-1].replace("'","")
                repr_file.close()
            return
        if '.tar.xz' in gridpack:
            self.__gridpack = gridpack
        else:
            self.__gridpack = gridpack
            print('gridpack location and name was not provided!\n Searching for gridpack containing *'+self.name+'* in the following directories:')
            for iloc,location in enumerate(self.__gridpack_locations):
                print("In "+location)
                print('-> found the following gridpacks:')
                counter = 0
                gridpacks = glob(location+'/*'+self.name+'*.tar.xz')
                for g in gridpacks:                    
                    print('[%i] '%counter+g)
                    counter += 1
                question_text = '([n] continue to next location with)' if iloc < len(self.__gridpack_locations) else '(last location. you have to choose!)'
                choice = ''
                while choice == '':
                    choice = input('Do you want to continue with one of these? %s\n'%(question_text))
                if choice.lower() != 'n':
                    self.__gridpack = gridpacks[int(choice)]
                    break
            if self.__gridpack == self.name:
                raise ValueError('You have to provide a gridpack!')

    def write_cms_configs(self):
        for i in range(self.__number_jobs):
            #os.system("mkdir "+self.workdir+'/'+str(i))
            with open(self.workdir+'/'+str(i)+'/GEN_cfg.py','w+') as GEN_cfg:
                with open(self.__cfg_py,'r') as GEN_cfg_template:
                    for l in GEN_cfg_template:
                        if('OUTPUTNAME' in l):
                            l=l.replace('OUTPUTNAME',self.workdir+'/'+self.name+'_'+str(i))
                        if('GRIDPACKPATH' in l):
                            l=l.replace('GRIDPACKPATH',self.gridpack)
                        if('NEVENTS' in l):
                            l=l.replace('NEVENTS',str(self.__number_events_per_job[i]))
                        if('RANDOMSEED' in l):
                            l=l.replace('RANDOMSEED',str(self.__randomseeds[i]))
                        GEN_cfg.write(l)
    
    def write_wrapper_script(self,fileindex=0):
        with open(self.workdir+'/wrapper_script.sh','w+') as wrapper_script:
            wrapper_script.write("""#!/bin/bash
export SCRIPTDIR="""+'/'.join(self.__cfg_py.split('/')[:-1])+'/'+"""
export WORKDIR="""+self.workdir+"""/${1}

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700
if [ -r CMSSW_10_6_19/src ] ; then
 echo CMSSW release already exists
else
 scram p CMSSW CMSSW_10_6_19
fi
cd CMSSW_10_6_19/src
eval `scram runtime -sh`
cd $WORKDIR/

cmsRun GEN_cfg.py
            """)
        os.system("chmod +x %s/wrapper_script.sh"%self.workdir)

    def write_htc_config(self):
        with open(self.workdir+"/htc_config.submit",'w+') as htc_config:
            htc_config.write("""#HTC Submission File for GEN sample production
fileindex         = $(Process)
#requirements      =  OpSysAndVer == "SL7"
universe          = vanilla
notification      = Error
notify_user       = andreas.hinzmann@desy.de
initialdir        = """+self.workdir+"""/$(fileindex)
output            = GEN_$(fileindex).o$(ClusterId).$(Process)
error             = GEN_$(fileindex).e$(ClusterId).$(Process)
log               = GEN.$(Cluster).log
#Requesting CPU and DISK Memory - default +RequestRuntime of 3h stays unaltered
+RequestRuntime   = 170000
#RequestMemory     = 10G
JobBatchName      = """+self.name+"""
RequestDisk       = 10G
getenv            = True
executable        = """+self.workdir+"""/wrapper_script.sh
arguments         = "$(fileindex)"
queue """+str(self.__number_jobs)+"""
""")

    def submit_jobs(self):
        cmd_list = ['condor_submit',self.workdir+'/htc_config.submit']
        print(cmd_list)
        try:
            process = subprocess.Popen(cmd_list,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            print('submitting..')
            proc_Out = process.communicate()
        except Exception:
            raise
        print('done!')
