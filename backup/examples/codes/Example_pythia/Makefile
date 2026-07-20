EXE=main
GCC=g++


PYTHIAINC=${shell pythia8-config --cxxflags}
PYTHIALIBS=${shell pythia8-config --libs}

ROOTINC=${shell root-config --cflags}
ROOTLIBS=${shell root-config --libs}

FASTJETINC=${shell fastjet-config --cxxflags}
FASTJETLIBS=${shell fastjet-config --libs}




all: $(EXE)

main01.o: main01.cc
	$(GCC) $(PYTHIAINC) $(ROOTINC) $(FASTJETINC) -c main01.cc

runglauber_v2.0_Dict.cxx: runglauber_v2.0.C LinkDef.h
	rootcint -f runglauber_v2.0_Dict.cxx -c runglauber_v2.0.C LinkDef.h

#runglauber.so: runglauber_v2.0_Dict.cxx
#	$(GCC) -fPIC --shared -o runglauber.so runglauber_v2.0_Dict.cxx $(ROOTINC)


#$(EXE): main01.o runglauber.so
#	$(GCC) $(PYTHIAINC) $(ROOTINC) $(FASTJETINC) -o $(EXE) main01.o runglauber.so $(PYTHIALIBS) $(ROOTLIBS) $(FASTJETLIBS) 

$(EXE): main01.o
	$(GCC) $(PYTHIAINC) $(ROOTINC) $(FASTJETINC) -o $(EXE) main01.o $(PYTHIALIBS) $(ROOTLIBS) $(FASTJETLIBS) 


clean:
	rm *.o $(EXE) runglauber.so runglauber_v2.0_Dict.cxx runglauber_v2.0_Dict.h
