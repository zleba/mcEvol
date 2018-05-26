
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)



CC=g++
CFLAGS=-g -std=c++11 -O3 -MMD -MP -I./inc  $(ROOTCFLAGS) -isystem./YODA/install/include/  \
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
	$(CC) -c -o $@ $< $(CFLAGS) -DisROOT=0


mcEvol: $(OBJS) 
	$(CC) -g -O3 $^  qcdnum/pij_nlo.f qcdnum/xpij2p.f qcdnum/xpns2p.f  qcdnum/ome.f qcdnum/wgplg.f -lgfortran  $(LINKLIBS)   ./YODA/install/lib/libYODA.so  -Wl,-rpath=$(PWD)/YODA/install/lib/   -o $@ 

clean:
	rm -f obj/*.o mcEvol


-include $(DEP)

