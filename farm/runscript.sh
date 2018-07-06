#!/bin/zsh
myDir=$2
orgDir=$PWD
cd $TMP
echo $TMP
. /cvmfs/sft.cern.ch/lcg/contrib/gcc/4.8/x86_64-slc6/setup.sh

#rm -f *.root
#cp $orgDir/../mcEvol .
cp $myDir/mcEvol .
./mcEvol
hadd hTot.root histo*.root
rm histo*.root
#mkdir -p $orgDir/alljobs/histos/$2
cp hTot.root $myDir/histos/mytmd_$1.root
