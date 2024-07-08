#!/bin/bash

# This script is used to calculate the cross section values for the blackmax generator, it should take in 5 arguments

# 1. planck mass (MD)
# 2. black hole mass (MBH)
# 3. number of extra dimensions (NDIM)
# 4. scenario (MODEL)
# 5. center of mass energy (COMENERGY)

WORKDIR=$PWD

if [[ $# -eq 5 ]]; then
    echo "Calculating cross section for blackmax generator"
    MD=$1
    MBH=$2
    NDIM=$3
    MODEL=$4
    COMENERGY=$5
    echo "Planck mass: $MD"
    echo "Black hole mass: $MBH"
    echo "Number of extra dimensions: $NDIM"
    echo "Scenario: $MODEL"
    echo "Center of mass energy: $COMENERGY"

else
    echo "Invalid number of arguments"
fi

# Set CMSSW version
VERCMSSW=CMSSW_12_4_11_patch3


# Set up the environment
source /cvmfs/cms.cern.ch/cmsset_default.sh
if [ -r $VERCMSSW/src ] ; then
    echo release $VERCMSSW already exists
else
    scram p CMSSW $VERCMSSW
fi

cd $VERCMSSW/src
eval `scram runtime -sh`
scram tool info lhapdf | grep LIBDIR | sed 's/LIBDIR=//'
LHAPDFLIBPATH=`scram tool info lhapdf | grep LIBDIR | sed 's/LIBDIR=//'`
echo "LHAPDF lib path = $LHAPDFLIBPATH"
cd -
mkdir gridpack_workdir
cd gridpack_workdir
# Get LHAPDF lib path


wget https://blackmax.hepforge.org/downloads/?f=BlackMax-2.02.0.tar.gz -O BlackMax-2.02.0.tar.gz
tar -xzf BlackMax-2.02.0.tar.gz
rm BlackMax-2.02.0.tar.gz

cd BlackMax
eval "sed -i '3s/.*/#COMP    = -fno-globals -fno-automatic -finit-local-zero/' Makefile"
eval "sed -i '11s:.*:PDFLIB = ${LHAPDFLIBPATH}:' Makefile"

gmake BlackMax

eval "sed -i '26s/.*/325300/' parameter.txt" # PDF ID
eval "sed -i '6s/.*/${COMENERGY}/' parameter.txt" # COM energy
eval "sed -i '32s/.*/${COMENERGY}/' parameter.txt" # COM energy

cd ../


cd BlackMax
eval "sed -i '8s/.*/${MD}/' parameter.txt" # M_D
eval "sed -i '30s/.*/${MBH}/' parameter.txt" # M_BH
eval "sed -i '16s/.*/${NDIM}/' parameter.txt" # n
eval "sed -i '14s/.*/${MODEL}/' parameter.txt" # hypothesis

eval "sed -i '2s/.*/10000/' parameter.txt" # number of events set to 10000

RANDOMSEED=$RANDOM
echo "Random seed: $RANDOMSEED"
eval "sed -i '60s/.*/${RANDOMSEED}/' parameter.txt" # random seed


./BlackMax 

cp BlackMaxLHArecord.lhe $WORKDIR/BlackMaxLHArecord_${MD}_${MBH}_${NDIM}_${MODEL}_${COMENERGY}.lhe

