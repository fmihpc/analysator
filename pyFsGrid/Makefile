DLL=libFsGrid.so
CC=g++
LNK=g++ 
OPT=-O3
CXXFLAGS=-std=c++17  -shared -fPIC


LIBRARY_PREFIX=/home/kostis/
INC_VLSV =-I${LIBRARY_PREFIX}vlsv
LIB_VLSV =-L${LIBRARY_PREFIX}vlsv/ -lvlsv


defaul:libFsGrid

clean: 
	rm  ${DLL} 

libFsGrid: libFsGrid.cpp
	${CC} ${CXXFLAGS} ${OPT}  libFsGrid.cpp   ${INC_VLSV}  ${LIB_VLSV} -o ${DLL}
