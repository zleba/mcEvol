
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)



CC=g++47
CFLAGS=-g -std=c++11 -O3 -MMD -MP -I./inc -I./usrInc $(ROOTCFLAGS)   \
                      -pedantic -W -Wall -Wshadow -Wno-long-long -fPIC 



DEPS = inc/Spline.h inc/SplineGen.h inc/alphaSpline.h inc/sudakovSpline.h \
        inc/tmd.h inc/Vec2.h inc/integration.h \
        inc/Sudakov.h inc/SplittingsIntegralSpline.h

SRCS = src/Spline.cpp src/SplineGen.cpp src/alphaSpline.cpp src/sudakovSpline.cpp \
       src/SplittingsLO.cpp src/SplittingsNLO.cpp src/SplittingsNNLO.cpp src/integration.cpp \
       src/Sudakov.cpp    src/tmd.cpp src/main.cpp


OBJS = obj/Spline.o obj/SplineGen.o obj/alphaSpline.o obj/sudakovSpline.o \
       obj/SplittingsLO.o obj/SplittingsNLO.o  obj/SplittingsNNLO.o  obj/integration.o \
       obj/Sudakov.o   obj/tmd.o  obj/main.o



DEP=$(OBJS:.o=.d)


LINKLIBS =     -ldl  $(ROOTLIBS)



obj/%.o: src/%.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

usrObj/%.o: usrSrc/%.cpp
	$(CC) -c -o $@ $< $(CFLAGS)


mcEvol: $(OBJS) 
	$(CC) -g -O3 $^  qcdnum/pij_nlo.f qcdnum/xpij2p.f qcdnum/xpns2p.f  qcdnum/ome.f qcdnum/wgplg.f -lgfortran  $(LINKLIBS)    -o $@ 

ver=1.7.0
yoda: 
	wget http://www.hepforge.org/archive/yoda/YODA-$(ver).tar.gz
	tar -xvzf YODA-$(ver).tar.gz
	mv YODA-$(ver) YODA
	cd YODA; ./configure CC=g++47 --without-zlib



-include $(DEP)

