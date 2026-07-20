#!/bin/bash
mkdir -p ../Data
make clean
make
cp ./lt ../Data
