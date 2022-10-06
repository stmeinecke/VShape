all:
	g++ -std=c++11 VShape.cpp -o VShape -Ofast -L $$PWD/incl/ -lfftw3 -lm -Wall
