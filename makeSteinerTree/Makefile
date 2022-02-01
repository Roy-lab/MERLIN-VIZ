#compiler: gcc for C programs, g++ for C++ programs
XX = g++
CC = gcc

#compiler flags
CFLAGS = -g

#subset of files

all: clean makeSteinerTree

makeSteinerTree:
	g++ -std=c++11 Framework.cpp Graph.cpp -o makeSteinerTree

clean:
	rm makeSteinerTree

