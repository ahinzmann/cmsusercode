#!/bin/sh
#
# source this file.
#
# Configure script to set up enviornment.
#
# Doug Gingrich 30 April-2022
#
export PYTHIAPATH=/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/pythia8/309-315ecb590794cc881e4b029bda922a92
export LHAPDFPATH=/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/lhapdf/6.1.6-giojec7
#
# Setup Pythia stuff.
#
export PYTHIA8DATA=${PYTHIAPATH}/share/Pythia8/xmldoc
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PYTHIAPATH}/lib
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:${PYTHIAPATH}/lib
#
# Setup up LHA PDF stuff.
#
export LHAPATH=${LHAPDFPATH}/share/lhapdf/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${LHAPDFPATH}/lib
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:${LHAPDFPATH}/lib
#
echo "PYTHIAPATH " $PYTHIAPATH
echo "LHAPDFPATH " $LHAPDFPATH
echo "PYTHIA8DATA" $PYTHIA8DATA
echo "LAHPATH    " $LHAPATH
echo "LD_LIBRARY_PATH" $LD_LIBRARY_PATH
echo "DYLD_LIBRARY_PATH" $DYLD_LIBRARY_PATH
#
#EOF

source /cvmfs/grid.desy.de/etc/profile.d/grid-ui-env.sh
source /cvmfs/cms.cern.ch/cmsset_default.sh
#export SCRAM_ARCH=slc7_amd64_gcc700
#cd ~/run2023/CMSSW_10_6_30/src/
cd ~/run2023/CMSSW_14_0_18/src/
cmsenv
cd /data/dust/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/qbh
