EXE=lbt
GCC=g++ -std=c++11

PYTHIAINC=${shell pythia8-config --cxxflags}
PYTHIALIBS=${shell pythia8-config --libs}

ROOTINC=${shell root-config --cflags}
ROOTLIBS=${shell root-config --libs}

FASTJETINC=${shell fastjet-config --cxxflags}
FASTJETLIBS=${shell fastjet-config --libs}

all: $(EXE)

LBT.o: LBT.cpp
	$(GCC) $(PYTHIAINC) $(ROOTINC) $(FASTJETINC) -c LBT.cpp

maintest.o: maintest.cpp
	$(GCC) $(PYTHIAINC) $(ROOTINC) $(FASTJETINC) -c maintest.cpp

SoftDrop.o: SoftDrop.cc
	$(GCC) $(PYTHIAINC) $(ROOTINC) $(FASTJETINC) -c SoftDrop.cc

$(EXE): LBT.o maintest.o SoftDrop.o
	$(GCC) -O3 $(PYTHIAINC) $(ROOTINC) $(FASTJETINC) -o $(EXE) LBT.o maintest.o SoftDrop.o $(PYTHIALIBS) $(ROOTLIBS) $(FASTJETLIBS)

.PHONY:clean
clean:
	-rm *.o $(EXE)

#g++ main01.cc -o main01 -I/public/home/heyy/Softwares/PYTHIA_8307/include -O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread  -L/public/home/heyy/Softwares/PYTHIA_8307/lib -Wl,-rpath,/public/home/heyy/Softwares/PYTHIA_8307/lib -lpythia8 -ldl
