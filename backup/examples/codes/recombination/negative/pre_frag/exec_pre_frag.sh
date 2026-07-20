#!/bin/sh                                                                            
rm -rf mm *.o *.dat
echo "***********************"
echo "Let's start compiling. "
echo "***********************"
g++ -o mm main-to-urqmd-all.cc
echo "=========================="
echo "The transformation started"
echo "=========================="
./mm 2000
