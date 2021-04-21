<<<<<<< HEAD
FLAGS=  -O3 -mavx -march=native -lm -ftree-vectorize -fopt-info-vec-optimized -fopenmp -pthread -DLIKWID_PERFMON -I/home/soft/likwid/include -L/home/soft/likwid/lib -llikwid
DEBUG= -Wall -g -DDEBUG
=======
FLAGS=  -O3 -mavx -march=native -lm
DEBUG= -Wall -g -DDEBUG 
>>>>>>> 7fdb1cab9e273cdab77817b9688380b6d232c565
CXX=gcc
RM=rm -rf

all: main.c pdelib.o  Makefile
	$(CXX) -o pdeSolver main.c pdelib.o $(FLAGS)

pdelib.o: pdelib.c  pdelib.h Makefile
	$(CXX) -c pdelib.c pdelib.h $(FLAGS)

clean:
	$(RM) pdeSolver *.o

doc: $(OBJ)
	doxygen ./Doxyfile
#doc: *.c trabalho-1.doxy *.h
<<<<<<< HEAD
#	@doxygen trabalho-1.doxy
=======
#	@doxygen trabalho-1.doxy
>>>>>>> 7fdb1cab9e273cdab77817b9688380b6d232c565
