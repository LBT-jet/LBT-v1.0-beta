#!/bin/sh                                                                            
make clean
rm -rf output
mkdir output
echo "***********************"
echo "Let's start compiling. "
echo "***********************"
make
echo "=========================="
echo "The fragmentation started"
echo "=========================="
./lt 2000
