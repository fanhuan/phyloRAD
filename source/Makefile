GCC         = gcc
CFLAGS_k    = -g -W -Wall -O2 -D_FILE_OFFSET_BITS=64
LIBS        = -lm
CC          = cc
CFLAGS_p    = -O3 -fomit-frame-pointer
GENERIC_SRC = string.h vector.h file_reader.h hashset.h

all:            fitch_kmerX consense treedist fitch_kmerX_long kmer_count kmer_countx kmer_merge ReadsSelector

phylip.o:       phylip.h
dist.o:         phylip.h dist.h
dist_long.o:         phylip.h dist.h
fitch_kmerX.o:  fitch_kmerX.c dist.h phylip.h
fitch_kmerX_long.o:  fitch_kmerX_long.c dist.h phylip.h

fitch_kmerX:    fitch_kmerX.o dist.o phylip.o
	$(CC) $(CFLAGS_p) fitch_kmerX.o dist.o phylip.o $(LIBS) -o fitch_kmerX
fitch_kmerX_long:    fitch_kmerX_long.o dist_long.o phylip.o
	$(CC) $(CFLAGS_p) fitch_kmerX_long.o dist_long.o phylip.o $(LIBS) -o fitch_kmerX_long
cons.o:         cons.c cons.h phylip.h 
consense.o:     consense.c cons.h phylip.h

consense:       consense.o phylip.o cons.o
	$(CC) $(CFLAGS_p) consense.o phylip.o cons.o $(LIBS) -o consense
treedist.o:     treedist.c cons.h phylip.h

treedist:     treedist.o phylip.o cons.o
	$(CC) $(CFLAGS_p) treedist.o cons.o phylip.o $(LIBS) -o treedist

kmer_count: $(GENERIC_SRC) file_reader.c kmer_count.c

	$(GCC) $(CFLAGS_k) -o kmer_count file_reader.c kmer_count.c $(LIBS)

kmer_countx: $(GENERIC_SRC) file_reader.c kmer_countx.c

	$(GCC) $(CFLAGS_k) -o kmer_countx file_reader.c kmer_countx.c $(LIBS)

kmer_merge: $(GENERIC_SRC) file_reader.c filter.c

	$(GCC) $(CFLAGS_k) -o kmer_merge file_reader.c filter.c $(LIBS)

ReadsSelector:
	g++ ReadsSelectorV2.cpp -o ReadsSelector

clean:
	rm -rf *.o *.dSYM





