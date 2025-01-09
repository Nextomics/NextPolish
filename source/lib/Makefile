HTSDIR = ./htslib
HTSLIB = $(HTSDIR)/libhts.a
HTSLIB_CPPFLAGS = -I$(HTSDIR)
LIBS= -lm -lbz2 -llzma -lz

CC       = gcc
AR       = ar
CFLAGS = -Wall -O3 -std=c99 -fPIC -s 
ALL_CFLAGS = -I. $(HTSLIB_CPPFLAGS) $(CFLAGS)
ALL_LIBS = -lpthread $(LIBS)

DEFAULT = clean nextpolish1.so calgs.so nextpolish2.so
OBJECTS = main.o snpvalid.o snpphase.o kmercount.o lgspolish.o scorechain.o contig.o base.o seqlist.o config.o

all: $(DEFAULT)

htslib_:
	@+make -C $(HTSDIR) CPPFLAGS+=-fPIC

nextpolish1.so: htslib_ $(OBJECTS)
	 $(CC) $(CFLAGS) -shared -o $@ $(OBJECTS) $(HTSLIB) $(ALL_LIBS)

nextpolish1: htslib_ $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS) $(HTSLIB) $(ALL_LIBS)

nextpolish2.so: bseq.c align.c bsort.c dag.c ctg_cns.c $(HTSLIB)
	$(CC) $(CFLAGS) -shared -Wno-unused-function -DLGS_CORRECT -D_SVID_SOURCE -o $@ $^ $(ALL_CFLAGS) $(ALL_LIBS)

calgs.so:
	$(CC) $(CFLAGS) -shared -o $@ calgs.c -lz

calgs:
	$(CC) $(CFLAGS) -o $@ calgs.c -lz

%.o : %.c
	$(CC) $(ALL_CFLAGS) -o $@ -c $<

clean_:
	rm -rf *.o *.pyc
	@+make clean -C $(HTSDIR)

clean:
	rm -rf *.o *.so *.pyc $(PROGS) calgs nextpolish1 nextpolish2 __pycache__
	@+make clean -C $(HTSDIR)
