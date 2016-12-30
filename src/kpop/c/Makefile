OPT = -O3
#OPT = -O3 -m64 -static
#OPT = -Wall -g 

CFLAGS = -Wall -pedantic
CC = gcc
LIBS = -lm

all: target

target: structure 

#valgrind: OPT = -g -O1
#valgrind: CFLAGS += -g
#valgrind: clean
#valgrind: target

#debug: OPT = -g -O1
#debug: CFLAGS += -g
#debug: clean
#debug: target

structure: structure.o params.o datain.o output.o ran.o mymath.o
	$(CC) -o structure structure.o params.o datain.o output.o ran.o mymath.o $(OPT) $(LIBS)

#STRAT: STRAT.o params.o datain.o ran.o mymath.o
#	$(CC) -o STRAT STRAT.o params.o datain.o ran.o mymath.o $(OPT) $(LIBS)

#STRAT.o: STRAT.c
#	$(CC) -c STRAT.c  $(OPT) $(CFLAGS)

structure.o: structure.c
	$(CC) -c structure.c $(OPT) $(CFLAGS)

output.o: output.c
	$(CC) -c output.c $(OPT) $(CFLAGS)

datain.o: datain.c
	$(CC) -c datain.c $(OPT) $(CFLAGS)

params.o: params.c
	$(CC) -c params.c $(OPT) $(CFLAGS)

ran.o: ran.c
	$(CC) -c ran.c $(OPT) $(CFLAGS)

mymath.o: mymath.c
	$(CC) -c mymath.c $(OPT) $(CFLAGS)
clean:
	@rm -f *.o structure
