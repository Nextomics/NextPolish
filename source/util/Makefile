BWADIR = ./bwa
SAMTOOLSDIR = ./samtools
MINIMAP2DIR = ./minimap2
HTSDIR = ../lib/htslib/

CC = gcc
CFLAGS = -Wall -O3 -s
PROGS = seq_split seq_count bwa_ samtools_ minimap2_

all: $(PROGS)

seq_split: kseq.h thpool.h seq_split.c thpool.c
	$(CC) $(CFLAGS) -pthread -o $@ seq_split.c thpool.c -lz 

seq_count: seq_count.c kseq.h
	$(CC) $(CFLAGS) -std=c99 -o $@ seq_count.c -lz

bwa_:
	@+make -C $(BWADIR)

samtools_:
	@+make -C $(SAMTOOLSDIR) HTSDIR=../$(HTSDIR)

minimap2_:
	@+make -C $(MINIMAP2DIR)

clean:
	rm -rf *.o $(PROGS)
	@+make clean -C $(BWADIR)
	@+make clean -C $(MINIMAP2DIR)
	@+make clean -C $(SAMTOOLSDIR) HTSDIR=../$(HTSDIR)
