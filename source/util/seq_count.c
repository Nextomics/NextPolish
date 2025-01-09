#define _GNU_SOURCE
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <getopt.h>
#include <zlib.h>
#include <libgen.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef struct {
	uint64_t depth;
	uint64_t filter_length;
	uint64_t genome_size;
} opt_t;

typedef struct {
	uint32_t *length;
	uint64_t max_size;
	uint64_t len;
	uint64_t total_bases;
	uint64_t filter_len;
	uint64_t filter_bases;
} reads;

typedef struct {
	uint64_t count[11];
	uint32_t length[11];
} nstat;

int cmpfunc (const void * a, const void * b)
{
   return ( *(int*)b - *(int*)a );
}

static void out_stat(opt_t opt, reads *read){
	qsort(read->length, read->len, sizeof(uint32_t), cmpfunc);

	uint64_t bin = read->len/1000;
	if (! bin) bin = 10;
	int bin_step = 1000;
	printf("[Read length histogram ('*' =~ %lu reads)]", bin);
	uint64_t bin_count = 0;
	int bin_index = 1;

	nstat nstat_;
	int nstat_index = 1;
	uint64_t nstat_bases = 0;
	memset(&nstat_, 0, sizeof(nstat));

	uint32_t seed_cutoff = 0;
	uint64_t seed_cov_len = (uint64_t) opt.depth * opt.genome_size;
	for (int i = 0; i < read->len; i ++){
		if (seed_cutoff == 0){
			seed_cov_len -= read->length[i];
			if (seed_cov_len <= read->length[i]) seed_cutoff = read->length[i];
		}

		if (read->length[read->len - 1 - i] < bin_step * bin_index) {
			bin_count ++;
		}else{
			do {
				printf("\n%7d %7d %10lu  ", bin_step * (bin_index - 1), bin_step * bin_index - 1, bin_count);
				for (int j = 0; j < bin_count/bin; j++){
					printf("*");
				}
				bin_count = 0;
				bin_index ++;
			}
			while (read->length[read->len - 1 - i] > bin_step * bin_index);
			bin_count = 1;
		}
		nstat_bases += read->length[i];
		nstat_.count[nstat_index - 1] ++;
		if (nstat_bases >= nstat_index * 0.1 * read->total_bases){
			nstat_.length[nstat_index - 1] = read->length[i];
			nstat_.count[nstat_index] += nstat_.count[nstat_index - 1];
			nstat_index ++;
		}
	}

	printf("\n%7d %7d %10lu  ", bin_step * (bin_index - 1), bin_step * bin_index - 1, bin_count);
	for (int j = 0; j < bin_count/bin; j++){
		printf("*");
	}


	printf("\n\n[Read length stat]\n");
	printf("%5s %20s %10s\n","Types", "Count (#)", "Length (bp)");
	for (int i = 0 ; i < 9; i ++){
		printf("N%-4d %20lu %7u\n",(i + 1) * 10, nstat_.count[i], nstat_.length[i]);
	}

	printf("\n%-8s %20s %20s %10s\n","Types", "Count (#)", "Bases (bp)", "Depth (X)");
	printf("%-8s %20lu %20lu %10.2f\n","Raw", read->len + read->filter_len, \
		read->total_bases + read->filter_bases, (read->total_bases + read->filter_bases)/(float) opt.genome_size);
	printf("%-8s %20lu %20lu %10.2f\n","Filtered", read->filter_len, read->filter_bases, read->filter_bases/(float) opt.genome_size);
	printf("%-8s %20lu %20lu %10.2f\n","Clean", read->len, read->total_bases, read->total_bases/(float) opt.genome_size);
	printf("\n*Suggested length cutoff of reads (genome size: %lu, expected seed depth: %lu) to be corrected: %u bp\n",\
		opt.genome_size, opt.depth, seed_cutoff);
}

static void count_reads(char *line, opt_t opt, reads *read)
{
	gzFile fp;
	fp = gzopen(line, "r");
	if (fp == NULL) {
		fprintf(stderr, "Error! %s does not exist!", line);
		exit(1);
	}
	kseq_t *seq;
	seq = kseq_init(fp);

	while(kseq_read(seq) >= 0) {
		if (seq->seq.l < opt.filter_length ) {
			read->filter_len ++;
			read->filter_bases += seq->seq.l;
		}else{
			read->length[read->len++] = seq->seq.l;
			read->total_bases += seq->seq.l;
			if (read->len >= read->max_size) {
				read->max_size += 1000000;
				read->length = realloc(read->length, read->max_size * sizeof (uint32_t));
			}
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
}

static inline uint64_t mm_parse_num(const char *str)
{
	double x;
	char *p;
	x = strtod(str, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9;
	else if (*p == 'M' || *p == 'm')  x *= 1e6;
	else if (*p == 'K' || *p == 'k')  x *= 1e3;
	return (uint64_t)(x + .499);
}

static int usage()
{
	fprintf(stderr, "Usage: seq_count [options] input.fofn\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, " -f filter length, default 1kb\n");
	fprintf(stderr, " -g genome size, default: 5Mb\n");
	fprintf(stderr, " -d expected seed depth, used to be corrected, default: 35\n");
	return 1;
}

int main(int argc, char *argv[])
{
	int c;
	opt_t opt;
	opt.depth = 35;
	opt.filter_length = 1000;
	opt.genome_size = 5000000;
	while((c = getopt(argc, argv, "f:g:d:")) != -1) {
		switch (c) {
			case 'f':
				opt.filter_length = mm_parse_num(optarg);
				break;
			case 'g':
				opt.genome_size = mm_parse_num(optarg);
				break;
			case 'd':
				opt.depth = mm_parse_num(optarg);
				break;
			default:
				return usage();
		}
	}
	if (optind + 1 > argc) return usage();

	char *fofn = strdup(argv[optind]);
	char *path = dirname(fofn);
	FILE *stream = fopen(argv[optind], "r");
	if (stream == NULL) {
		fprintf(stderr, "Failed open input file list!\n");
		return usage();
	}

	reads read;
	memset(&read, 0, sizeof(reads));
	read.max_size = 50000000;
	read.length = malloc(sizeof(uint32_t) * read.max_size);

	char *line = NULL;
	size_t len = 0;
	ssize_t nread;
	while((nread = getline(&line, &len, stream)) != -1) {
		if (nread != 1 && line[0] != '#') {
			line[strlen(line) - 1] = '\0';
			if (line[0] == '/')
				count_reads(line, opt, &read);
			else {
				char file_path[1024];
				sprintf(file_path, "%s/%s", path, line);
				count_reads(file_path, opt, &read);
			}
		}
	}
	free(line);
	fclose(stream);
	if (fofn) free (fofn);
	if (read.len) out_stat(opt, &read);
	free (read.length);
	return 0;
}
