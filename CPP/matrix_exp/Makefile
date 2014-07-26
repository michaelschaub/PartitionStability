#TODO SIMPLIFY COMPILING/AUTOMATE WITH VARIABLES ETC...

#Define Compiler -- adapt if you are not using gcc/g++
CC=g++
CF=gfortran

#Compilation flags
CFLAGS=-O3 -std=c++0x
FFLAGS=-O3
#-Wall

all: matrix_exp_input matrix_exp_input_test matrix_exp_test
	
matrix_exp_input: expokit.o
	$(CC) $(CFLAGS) matrix_exp_input.cpp expokit.o -lgfortran -lblas -llapack -o matrix_exp_input.sh

matrix_exp_input_test: expokit.o
	$(CC) $(CFLAGS) matrix_exp_input_test.cpp expokit.o -lgfortran -lblas -llapack -o matrix_exp_input_test.sh

matrix_exp_test:
	$(CC) $(CFLAGS) matrix_exp_test.cpp expokit.o -lgfortran -lblas -llapack -o matrix_exp_test.sh 

expokit.o: 
	$(CF) $(FFLAGS) -c ../../FORTRAN/expokit.f

clean:
	rm *.sh *.o
