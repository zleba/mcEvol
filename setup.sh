. /cvmfs/sft.cern.ch/lcg/contrib/gcc/4.8/x86_64-slc6/setup.sh
#cd /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.06.08/x86_64-slc6-gcc48-opt/root/
#. ./bin/thisroot.sh
#cd -

#PATH=/cvmfs/sft.cern.ch/lcg/external/Python/2.7.4/x86_64-slc6-gcc48-opt/bin:$PATH
#PYTHONPATH=/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/py2-numpy/1.12.1-mlhled2/lib/python2.7/site-packages:$PYTHONPATH
#PYTHONPATH=/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/py2-matplotlib/1.5.2-njopjo/lib/python2.7/site-packages:$PYTHONPATH
PYTHONPATH=$PWD/YODA/install/lib64/python2.6/site-packages:$PYTHONPATH
LD_LIBRARY_PATH=$PWD/YODA/install/lib/:$LD_LIBRARY_PATH
