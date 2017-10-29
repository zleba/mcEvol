
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)

LHALIBDIR = /home/radek/Dropbox/Lund/mcnet/lhapdf/lhadir/lib
LHAINCDIR = /home/radek/Dropbox/Lund/mcnet/lhapdf/lhadir/include


CC=g++47
CFLAGS=-g -std=c++11 -O3 -MMD -MP -I./inc -I./usrInc $(ROOTCFLAGS)   \
                     -I$(LHAINCDIR) -pedantic -W -Wall -Wshadow -Wno-long-long -fPIC 



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

LHApdf = -L/home/radek/Dropbox/Lund/mcnet/lhapdf/lhadir/lib -lLHAPDF 

LINKLIBS =     -ldl  $(ROOTLIBS)



obj/%.o: src/%.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

usrObj/%.o: usrSrc/%.cpp
	$(CC) -c -o $@ $< $(CFLAGS)


myTMD: $(OBJS) 
	$(CC) -g -O3 $^  qcdnum/pij_nlo.f qcdnum/xpij2p.f qcdnum/xpns2p.f  qcdnum/ome.f qcdnum/wgplg.f -lgfortran  $(LINKLIBS)    -o myTMD


-include $(DEP)

disc/a.out: disc/disc.cpp src/integration.cpp
	$(CC) -g -O2 $(CFLAGS) $^  qcdnum/ome.o qcdnum/wgplg.o qcdnum/xa2hgp.o $(LINKLIBS)  -lgfortran  -Wl,-R$(LHALIBDIR)   -o disc/a.out
