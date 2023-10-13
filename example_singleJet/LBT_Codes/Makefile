EXE=lbt
GCC=g++

PYTHIAINC=${shell pythia8-config --includedir}
PYTHIALIBS=${shell pythia8-config --libdir}

#ROOTINC=${shell root-config --cflags}
#ROOTLIBS=${shell root-config --libs}

FASTJETINC=${shell fastjet-config --cxxflags}
FASTJETLIBS=${shell fastjet-config --libs}

all: $(EXE)

LBT.o: LBT.cpp
#	$(GCC) -I$(PYTHIAINC) $(ROOTINC) $(FASTJETINC) -c LBT.cpp
	$(GCC) -I$(PYTHIAINC) $(FASTJETINC) -c LBT.cpp

main.o: main.cpp
#	$(GCC) -I$(PYTHIAINC) $(ROOTINC) $(FASTJETINC) -c main.cpp
	$(GCC) -I$(PYTHIAINC) $(FASTJETINC) -c main.cpp

$(EXE): LBT.o main.o
#	$(GCC) -I$(PYTHIAINC) $(ROOTINC) $(FASTJETINC) -o $(EXE) LBT.o main.o $(PYTHIALIBS)/libpythia8.a $(ROOTLIBS) $(FASTJETLIBS) 
	$(GCC) -I$(PYTHIAINC) $(FASTJETINC) -o $(EXE) LBT.o main.o $(PYTHIALIBS)/libpythia8.a $(FASTJETLIBS) 

.PHONY:clean
clean:
	-rm *.o $(EXE)
