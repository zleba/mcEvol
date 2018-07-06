
export MYDIR=/nfs/dust/cms/user/zlebcr/tmdOut

#tags="zmax1e_2 zmax1e_3 zmax1e_4 zmax1e_5"
tags="zmax1e_6"

for i in $tags
do
    hadd histos/${i}.root $MYDIR/$i/histos/*.root

done
