ver=1.6.7

mkdir YODA
cd YODA
wget http://www.hepforge.org/archive/yoda/YODA-${ver}.tar.gz
tar -xvzf YODA-${ver}.tar.gz

mkdir install
instalPath=${PWD}/install
cd YODA-${ver}

sed -i.bak 's/return map(int, s\.strip()\.split("\."))/import re; return map(lambda n: int(re.sub("[^0-9]", "", n)), s.strip().split("."))/' configure  

./configure  --prefix=$instalPath    --without-zlib
make -j2 && make -j2 install
