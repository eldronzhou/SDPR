
CC = g++

CFLAGS = -g3 -O3 -DHAVE_INLINE -march=native -Igsl/include -std=c++11 -Wall -Wextra -pedantic

all: SDPR

SDPR: SDPR_io.o LD.o parse_gen.o mcmc.o function_pool.o main.o
	${CC} ${CFLAGS} SDPR_io.o LD.o parse_gen.o mcmc.o function_pool.o main.o -Lgsl/lib/ -Wl,-rpath gsl/lib -lgsl -LMKL/lib/ -Wl,--no-as-needed,-rpath MKL/lib/ -lmkl_rt -lm -lpthread -ldl -o SDPR

SDPR_io.o: SDPR_io.cpp SDPR_io.h
	${CC} ${CFLAGS} -c SDPR_io.cpp

LD.o: LD.cpp LD.h SDPR_io.h
	${CC} ${CFLAGS} -pthread -c LD.cpp

parse_gen.o: parse_gen.cpp parse_gen.h
	${CC} ${CFLAGS} -c parse_gen.cpp

function_pool.o: function_pool.cpp function_pool.h
	${CC} ${CFLAGS} -c function_pool.cpp

mcmc.o: mcmc.cpp mcmc.h parse_gen.h function_pool.h sse_mathfun.h
	${CC} ${CFLAGS} -pthread -c mcmc.cpp

main.o: main.cpp parse_gen.h mcmc.h
	${CC} ${CFLAGS} -c main.cpp

clean:
	rm -f *.o
