EXE=dNdpT
GCC=g++

ROOTINC=${shell root-config --cflags}
ROOTLIBS=${shell root-config --libs}

FASTJETINC=${shell fastjet-config --cxxflags}
FASTJETLIBS=${shell fastjet-config --libs}

PYTHIAINC=${shell pythia8-config --cxxflags}
PYTHIALIBS=${shell pythia8-config --libs}

all: $(EXE)

dNdpT.o: dNdpT.cc
	$(GCC) $(ROOTINC) $(FASTJETINC) -c dNdpT.cc

$(EXE): dNdpT.o
	$(GCC) $(ROOTINC) $(FASTJETINC) -o $(EXE) dNdpT.o $(ROOTLIBS) $(FASTJETLIBS)

.PHONY:clean
clean:
	-rm *.o $(EXE)
