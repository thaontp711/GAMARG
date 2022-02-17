

CFLAGS= -O3 -c -Wall
gzFLAGS= -lstdc++ -lgzstream -lz

CXX      = g++   # for Linux RedHat 6.1, g++ version 2.95.2
AR       = ar cr
CPPFLAGS = -I. -O

all: gamarg

gamarg: main.o utilities.o CData.o COneNode.o CGAMARG.o CARG.o Match.o #CVCF.o main.o VCF_utilities.o utilities.o CCDif.o
	g++ ${CPPFLAGS} -O3 -Wall gzstream/gzstream.o main.o CData.o utilities.o COneNode.o Match.o CGAMARG.o CARG.o  -lstdc++ -lpthread -lz  -o gamarg

main.o: main.cpp
	g++ ${CPPFLAGS} -c -Wall main.cpp


utilities.o: utilities.h
	g++ -c -Wall utilities.cpp


CData.o:
	g++ -c -Wall CData.cpp

CGAMARG.o:
	g++ -c -Wall CGAMARG.cpp

CARG.o:
	g++ -c -Wall CARG.cpp

Match.o:
	g++ -c -Wall Match.cpp


COneNode.o:
	g++ -c -Wall COneNode.cpp

gzstream.o : gzstream/gzstream.C
	${CXX} ${CPPFLAGS} -c -o gzstream/gzstream.o gzstream/gzstream.C
libgzstream.a : gzstream.o
	${AR} gzstream/libgzstream.a gzstream/gzstream.o

clean:
	rm -rf *.o gamarg
