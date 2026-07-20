EXE=lt
GCC=g++ -std=c++11

PYTHIAINC=${shell pythia8-config --cxxflags}
PYTHIALIBS=${shell pythia8-config --libs}

ROOTINC=${shell root-config --cflags}
ROOTLIBS=${shell root-config --libs}

FASTJETINC=${shell fastjet-config --cxxflags}
FASTJETLIBS=${shell fastjet-config --libs}

all: $(EXE)

main100.o: main100.cc
	$(GCC) $(PYTHIAINC) $(ROOTINC) $(FASTJETINC) -c main100.cc

$(EXE): main100.o
	$(GCC) -O3 $(PYTHIAINC) $(ROOTINC) $(FASTJETINC) -o $(EXE)  main100.o  $(PYTHIALIBS) $(ROOTLIBS) $(FASTJETLIBS) 
	#$(GCC) -O3 $(PYTHIAINC) $(ROOTINC) $(FASTJETINC) -o $(EXE) LBT.o maintest.o SoftDrop.o $(PYTHIALIBS) $(ROOTLIBS) $(FASTJETLIBS)
	
.PHONY:clean
clean:
	-rm *.o $(EXE)
