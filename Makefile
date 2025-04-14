
#CC = gcc
CC = clang
CFLAGS = -Wall  -g  -mtune=native
CVECFLAGS := $(CFLAGS) -mavx2 -ftree-vectorize #-fopt-info-vec-optimized
CFLAGS += -fno-tree-vectorize
NFL_AVX = -DNFL_OPTIMIZED=ON -DNTT_AVX2

EXEC = signature_tests sampling_tests timing main_signature mainprotocol

OBJ = arithmetic.o random.o crt_trees.o preprocessing.o protocolparties.o

HDR = common.h

all: sampling_tests signature_tests

mainprotocol: mainprotocol.o preprocessing.o protocolparties.o crt_trees.o sampling.o random.o arithmetic.o cpucycles.o
	$(CC) $(CFLAGS) -o $@ $^ -lm


arithmetic.o: arithmetic.c random.o $(HDR)
	$(CC) $(CVECFLAGS) -c -o $@ $<

sampling.o: sampling.c random.o $(HDR)
	$(CC) $(CVECFLAGS) -c -o $@ $<

random.o: random.c
	$(CC) $(CVECFLAGS) -maes -c -o $@ $<

%.o: %.c $(HDR)
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(EXEC) *.o *.s
