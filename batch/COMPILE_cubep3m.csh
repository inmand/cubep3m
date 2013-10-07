##!/bin/csh

module list

cd ../source_threads

make clean
make -f Makefile

cd ../batch

echo "Sourced COMPILE_cubep3m.csh"

