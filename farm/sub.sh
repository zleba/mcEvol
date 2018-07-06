export MYDIR=/nfs/dust/cms/user/zlebcr/tmdOut/zmax1e_6
export NJOBS=300
mkdir -p $MYDIR/logs
mkdir -p $MYDIR/histos
cp ../mcEvol $MYDIR
condor_submit jobs.submit
