#Define Compiler -- adapt if you are not using gcc/g++
CC=g++

#Compilation flags
CFLAGS=-O3 -std=c++0x
#-Wall

all: matrix_exp_test_input matrix_exp_test
	
matrix_exp_test_input:
	#$(CC) $(CFLAGS) matrix_exp_input.cpp -I../cliques -lemon -o matrix_exp_input.sh
	$(CC) $(CFLAGS) matrix_exp_input.cpp r8lib.cpp c8lib.cpp matrix_exponential.cpp -I../cliques -lemon -o matrix_exp_input.sh

matrix_exp_test:
	$(CC) $(CFLAGS) matrix_exp_test.cpp r8lib.cpp c8lib.cpp matrix_exponential.cpp -o matrix_exp_test.sh 

clean:
	rm *.sh 
