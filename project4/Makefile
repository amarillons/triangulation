LIB += -framework Cocoa -framework OpenGL -framework GLUT -framework IOKit -framework Carbon -framework CoreVideo -L/usr/local/lib/ -lglfw3
# LIB+= -framework IOKit -framework Carbon -framework Cocoa -framework CoreVideo -lglfw3

CC = g++
OPT += -O3 -Wall
# OPT += -Wno-logical-op-parentheses
OBJS2 = randomsphere.o

all: prog2

clean:
	rm -vf *.o

prog1:
	$(CC) $(OPT) $(PREDEF) -c delaunay.cpp
	$(CC) delaunay.o $(LIB) -o delaunay.out
	$(CC) $(OPT) $(PREDEF) -c try2.cpp
	$(CC) try2.o $(LIB) -o try2.out
	$(CC) $(OPT) $(PREDEF) -c try3.cpp
	$(CC) try3.o $(LIB) -o try3.out
	$(CC) $(OPT) $(PREDEF) -c try4.cpp
	$(CC) try4.o $(LIB) -o try4.out
	$(CC) $(OPT) $(PREDEF) -c more.cpp
	$(CC) more.o $(LIB) -o more.out
prog2:
	$(CC) $(OPT) $(PREDEF) -c flip.cpp
	$(CC) flip.o $(LIB) -o flip.out
	$(CC) $(OPT) $(PREDEF) -c index.cpp
	$(CC) index.o $(LIB) -o index.out





