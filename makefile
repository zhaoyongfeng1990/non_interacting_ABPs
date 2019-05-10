CXX = mpic++
CXXFLAGS = -march=native -ffast-math -std=c++11 -Wno-literal-suffix -fabi-version=0 -O3# -ggdb
LFLAGS = -march=native -ffast-math -lm -std=c++11 -Wno-literal-suffix -fabi-version=0 -O3 #-ggdb

Objects = main.o# snapshot.o aveRadical.o
RandomGen = ./vectorclass/ranvec1.o
RandomObjs = #iterate.o initialization.o growth.o

solid_gas : $(Objects) $(RandomGen) $(RandomObjs)
	$(CXX) -o solid_gas $(Objects) $(RandomGen) $(RandomObjs) $(LFLAGS)

$(Objects) : config.h particles.h #mt64.h
$(RandomGen) : #mt64.h
$(RandomObjs) : #lattice2d.h mt64.hr

.PHONY : clean
clean :
	-rm *.o ./vectorclass/*.o
