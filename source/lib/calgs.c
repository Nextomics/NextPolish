#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <zlib.h>
#include "../util/kseq.h"
KSEQ_INIT(gzFile, gzread)

uint64_t calgs(const char *file){
	gzFile fp; 
	fp = gzopen(file, "r");
	if (fp == NULL) {
		fprintf(stderr, "Error! %s does not exist!", file); 
		exit(1);
	}
	uint64_t gs = 0;
	kseq_t *seq;
	seq = kseq_init(fp);
	while(kseq_read(seq) >= 0) {
		gs += seq->seq.l;
	}
	kseq_destroy(seq);
	gzclose(fp);
	return gs;
}

int main(int argc, char *argv[])
{
	uint64_t gs = calgs(argv[1]);
	printf("genome size: %lu bp\n", gs);
	return 0;
}
