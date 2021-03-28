#!/usr/bin/env python3
import os
from condor_submission import GEN
from optparse import OptionParser

parser = OptionParser()
parser.add_option('--submit', action='store_true',default=True, dest='submit',help='Submit jobs')
parser.add_option('--clean', action='store_true',default=False, dest='clean',help='Delete workdirs')
(options, args) = parser.parse_args()

print("preparing scripts and workdirs for submission of GEN sample production...")
gridpack_path = '/nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/'
gridpack_suffix = '_slc7_amd64_gcc700_CMSSW_10_6_19_tarball.tar.xz'

names=[
    'tripleG_QCD_HT2000toInf',
    'tripleG_QCD_HT4000toInf',
    'alp_QCD_HT2000toInf',
    'alp_QCD_HT2000toInf'
]

workdir='/nfs/dust/cms/user/hinzmann/dijetangular/madgraph3/'

if(options.clean):
    print('deleting workdirs..')
    os.system("rm -rf "+workdir)
    exit(0)

for name in names:
    g = GEN(name,workdir+name,gridpack=gridpack_path+name+gridpack_suffix,number_jobs=10,events=100000)
    if(options.submit):
        g.submit_jobs()
