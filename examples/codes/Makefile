EXE=lbt
GCC=g++

all: $(EXE)

LBT.o: LBT.cpp
	$(GCC) -c LBT.cpp

main.o: main.cpp
	$(GCC) -c main.cpp

$(EXE): LBT.o main.o
	$(GCC) -o $(EXE) LBT.o main.o 

.PHONY:clean
clean:
	-rm *.o $(EXE)
