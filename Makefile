CPP = g++
CFLAGS = -O3 -march=native -mtune=native -Wall -ggdb -I NFLlib/include/ -I NFLlib/include/nfl -I NFLlib/include/nfl/prng -I include -DNFL_OPTIMIZED=ON -DNTT_AVX2
INCLUDES = include/bench.h include/cpucycles.h
BENCH = src/bench.c src/cpucycles.c
TEST = src/test.c
BLAKE3 = src/blake3/blake3.c src/blake3/blake3_dispatch.c src/blake3/blake3_portable.c src/blake3/blake3_sse2_x86-64_unix.S src/blake3/blake3_sse41_x86-64_unix.S src/blake3/blake3_avx2_x86-64_unix.S src/blake3/blake3_avx512_x86-64_unix.S
LIBS = deps/libnfllib_static.a -lgmp -lmpfr -L deps/ -lflint -lquadmath

all: bdlop bgv ghl pior pibdn linear

bdlop: src/bdlop.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -DUSERS=10000 -c src/bgv.cpp -o bgv1.o
	${CPP} ${CFLAGS} -DMAIN src/bdlop.cpp bgv1.o ${TEST} ${BENCH} -o bdlop ${LIBS}

bgv: src/bgv.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -DUSERS=10000 -DMAIN src/bgv.cpp ${TEST} ${BENCH} -o bgv ${LIBS}

ghl: src/ghl.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -DUSERS=10000 -DMAIN src/ghl.cpp ${TEST} ${BENCH} -o ghl ${LIBS}

pior: src/pior.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -c src/sample_z_small.c -o sample_z_small.o
	${CPP} ${CFLAGS} -DUSERS=10000 -DMAIN src/pior.cpp ${BLAKE3} sample_z_small.o ${TEST} ${BENCH} -o pior ${LIBS}

pibdn: src/pibdn.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -c src/sample_z_small.c -o sample_z_small.o
	${CPP} ${CFLAGS} -c src/sample_z_large.c -o sample_z_large.o
	${CPP} ${CFLAGS} -DUSERS=10000 -DMAIN src/pibdn.cpp ${BLAKE3} sample_z_small.o sample_z_large.o ${TEST} ${BENCH} -o pibdn ${LIBS}

linear: src/linear.cpp src/linear.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -c src/sample_z_small.c -o sample_z_small.o
	${CPP} ${CFLAGS} -c src/bdlop.cpp -o bdlop.o
	${CPP} ${CFLAGS} -DUSERS=10000 -DMAIN src/linear.cpp ${BLAKE3} sample_z_small.o bdlop.o ${TEST} ${BENCH} -o linear ${LIBS}

clean:
	rm *.o bdlop bgv ghl pior pibdn linear
