#!/bin/bash
#
g++ -c -g -I/.  matrix_exponential.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling matrix_exponential.cpp"
  exit
fi
rm compiler.txt
#
echo "Library installed as ~/libcpp/$ARCH/matrix_exponential.o"
