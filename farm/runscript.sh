#!/bin/zsh
orgDir=$PWD
cd $TMP
echo $TMP
. /cvmfs/sft.cern.ch/lcg/contrib/gcc/4.8/x86_64-slc6/setup.sh

rm -f *.root
$orgDir/../mcEvol
hadd hTot.root histo*.root
rm histo*.root
cp hTot.root $orgDir/histos/mytmd_$1.root
