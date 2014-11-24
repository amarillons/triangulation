LIB += -framework Cocoa -framework OpenGL -framework GLUT -framework IOKit -framework Carbon -framework CoreVideo -L/usr/local/lib/ -lglfw3
# LIB+= -framework IOKit -framework Carbon -framework Cocoa -framework CoreVideo -lglfw3

CC = g++
OPT += -O3 -Wall
# OPT += -Wno-logical-op-parentheses
OPT += -ferror-limit=1
OBJS2 = randomsphere.o

all: prog1

clean:
	rm -vf *.o

prog1:
	$(CC) $(OPT) $(PREDEF) -c spider.cpp
	$(CC) spider.o $(LIB) -o spider.out





