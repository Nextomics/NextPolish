#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <zlib.h>
#include <libgen.h>
#include <sys/stat.h>
#include <pthread.h>
#include "thpool.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define CHAR_NUM 0x400
#define MAX_READS 0x4c4b40
#define MAX_SGS_READ_LEN 300

typedef struct {
	char (*file)[CHAR_NUM];
	int index;
	int len;
	int max_size;
} f;

typedef struct {
	f f1, f2;
	int thread;
	uint64_t max_mem;
	char out[CHAR_NUM], outdir[CHAR_NUM], outpre[CHAR_NUM];
	uint32_t min_length;
	uint32_t max_length;
	uint64_t max_base;
	int subfile;
	int rmdup;
	int rmn;
} opt_t;

typedef struct {
	gzFile fp;
	uint32_t count;
	int validy;
} fout;

typedef struct {
	char *buf;
	int validy;
	uint32_t buf_index;
	fout *fp;
} buffer;

static pthread_mutex_t buflock;
static pthread_mutex_t fplock;

static int get_buf_index(buffer *buf, int thread){
	pthread_mutex_lock(&buflock);
	int i, j;
	i = 1;
	while (i){
		for (j = 0; j < thread; j ++){
			if (buf[j].validy == 0) {
				i = 0;
				buf[j].validy = 1;
				break;
			}
		}
		if (i) sleep (1);
	}
	pthread_mutex_unlock(&buflock);
	return j;
}

static int get_fp_index(fout *fp, int subfile){
	pthread_mutex_lock(&fplock);
	int i, j, k;
	uint32_t count;
	i = 1;
	k = count = -1;
	while (i){
		for (j = 0; j < subfile; j ++){
			if (fp[j].validy == 0 && (count == -1 || count > fp[j].count)) {
				i = 0;
				k = j;
				count = fp[j].count;
			}
		}
		if (i) sleep (1);
	}
	fp[k].count ++;
	fp[k].validy = 1;
	pthread_mutex_unlock(&fplock);
	return k;
}

static int usage()
{
	fprintf(stderr, "Usage: seq_split [options] input.fofn\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "  -t number of threads to use [8]\n");
	fprintf(stderr, "  -m maximum memory per thread; suffix K/M/G [500M]\n");
	fprintf(stderr, "  -i interleaved paired-end files [1]\n");
	fprintf(stderr, "  -f minimum reads length [50]\n");
	fprintf(stderr, "  -l maximum reads length [inf]\n");
	fprintf(stderr, "  -s total base number to output [inf]\n");
	fprintf(stderr, "  -n subfile number [10]\n");
	fprintf(stderr, "  -N don't discard a read/pair if the read contains N base\n");
	// fprintf(stderr, "  -r remove duplicate reads [0]\n");
	fprintf(stderr, "  -p output directory [input.part]\n");
	fprintf(stderr, "  -d prefix of subfiles [$CWD]\n");
	return 1;
}

static void init_opt(opt_t* opt){
	opt->thread = 8;
	opt->max_mem = 500<<20;
	opt->rmdup = 0;
	opt->min_length = 50;
	opt->max_length = 0;
	opt->max_base = 0;
	opt->subfile = 10;
	opt->rmn = 1;
	getcwd(opt->outdir, sizeof(opt->outdir));
	strcpy(opt->outpre, "input.part");

	opt->f1.max_size = opt->f2.max_size = 20;
	opt->f1.len = opt->f2.len = opt->f1.index = opt->f2.index = 0;
	opt->f1.file = malloc(opt->f1.max_size * sizeof(char [CHAR_NUM]));
	opt->f2.file = malloc(opt->f2.max_size * sizeof(char [CHAR_NUM]));
}

static void realloc_f(f* f_){
	f_->max_size += 20;
	f_->file = realloc(f_->file, f_->max_size * sizeof(char [CHAR_NUM]));
}

static void destroy_opt(opt_t* opt){
	free (opt->f1.file);
	free (opt->f2.file);
}

static buffer *init_buf(opt_t *opt){
	int i;
	buffer *buf = calloc(opt->thread, sizeof(buffer));
	for (i = 0; i < opt->thread; i++){
		buf[i].buf = malloc(opt->max_mem);
	}
	return buf;
}

static void destroy_buf(opt_t *opt, buffer *buf){
	int i;
	for (i = 0; i < opt->thread; i++){
		free (buf[i].buf);
	}
	free (buf);
}

static fout *init_fp(opt_t *opt){
	int i;
	char fn[CHAR_NUM] = {0};
	char *suffix = opt->f2.len ? "fastq.gz" : "fasta.gz";
	fout *fp = calloc(opt->subfile, sizeof(fout));
	for (i = 0; i < opt->subfile; i++){
		sprintf(fn, "%s.%03d.%s", opt->out, i, suffix);
		fp[i].fp = gzopen(fn, "wb");
		if (fp[i].fp == NULL){
			fprintf(stderr, "Failed create subfile: %s\n", fn);
			exit(1);
		}
	}
	return fp;
}

static void destroy_fp(opt_t *opt, fout *fp){
	int i;
	for (i = 0; i < opt->subfile; i++){
		gzclose(fp[i].fp);
	}
	free (fp);
}

static gzFile get_file_h(f *f_){
	gzFile fp;
	if (f_->index >= f_->len){
		fp = NULL;
	}else{
		fp = gzopen(f_->file[f_->index++], "r");
		if (fp == NULL) {
			fprintf(stderr, "Error! %s does not exist!", f_->file[f_->index-1]);
			exit(1);
		}
	}
	return fp;
}


static int init_opt_f(char *fofn, int i, opt_t *opt){
	FILE *stream = fopen(fofn, "r");
	if (stream == NULL) {
		fprintf(stderr, "Failed open input file list: %s\n", fofn);
		exit (1);
	}

	char *fofn_ = strdup(fofn);
	char *dir = dirname(fofn_);
	int t = 0;
	char *line = NULL;
	size_t len = 0;
	ssize_t nread;
	char file_path[CHAR_NUM];
	while((nread = getline(&line, &len, stream)) != -1) {
		if (nread != 1 && line[0] != '#') {
			line[strlen(line) - 1] = '\0';
			if (line[0] == '/'){
				strcpy(file_path, line);
			}else{
				int k = sprintf(file_path, "%s/%s", dir, line);
				file_path[k] = '\0';
			}

			if (i && t%2) {
				strcpy(opt->f2.file[opt->f2.len++], file_path);
				if (opt->f2.len >= opt->f2.max_size) realloc_f(&opt->f2);
 			}else{
				strcpy(opt->f1.file[opt->f1.len++], file_path);
				if (opt->f1.len >= opt->f1.max_size) realloc_f(&opt->f1);
			}
			t ++;
		}
	}
	if (line) free(line);
	free (fofn_);
	fclose(stream);
	return 0;
}


static void save_buf (buffer *buf){
	gzwrite(buf->fp->fp, buf->buf, buf->buf_index);
	buf->fp->validy = 0;
	buf->validy = buf->buf_index = 0;
}

static inline int check_contain_N (kseq_t *seq1, kseq_t *seq2){
	int i;
	for (i = 0; i < seq1->seq.l; i ++){
		if (seq1->seq.s[i] == 'N') return 1;
	}

	if (seq2){
		for (i = 0; i < seq2->seq.l; i ++){
			if (seq2->seq.s[i] == 'N') return 1;
		}
	}
	return 0;
}

static void split_data(opt_t *opt, buffer *buf, fout *fp){

	int i, flag = 0;
	uint64_t total_base, total_read, total_nread;
	total_base = total_read = total_nread = 0;

	gzFile fp1 = get_file_h(&opt->f1);
	gzFile fp2 = get_file_h(&opt->f2);
	
	kseq_t *seq1, *seq2;
	seq1 = seq2 = NULL;
	seq1 = kseq_init(fp1);
	if (fp2) seq2 = kseq_init(fp2);
	uint64_t max_buf_index = opt->max_mem - MAX_READS;

	threadpool thpool_out = thpool_init(opt->thread);
	assert (!pthread_mutex_init(&buflock, NULL));
	assert (!pthread_mutex_init(&fplock, NULL));
	int index = get_buf_index(buf, opt->thread);

	while (fp1){
		while (1){
			flag = kseq_read(seq1);
			total_read ++;
			// total_base += seq1->seq.l;
			if (fp2) {
				kseq_read(seq2);
				// total_base += seq2->seq.l;
			}

			if (opt->max_base && total_base >= opt->max_base) break;

			if (flag < 0 ) {
					kseq_destroy(seq1);
					gzclose(fp1);
				if (fp2) {
					kseq_destroy(seq2);
					gzclose(fp2);
				}
				break;
			}

			if (opt->rmn && seq1->seq.l < MAX_SGS_READ_LEN && ((!fp2) || seq2->seq.l < MAX_SGS_READ_LEN) 
				&& check_contain_N(seq1, seq2)) {
				total_nread ++;
				continue;
			}

			if (seq1->seq.l >= opt->min_length && (seq1->seq.l <= opt->max_length || !opt->max_length)){
				total_base += seq1->seq.l;
				if (fp2 && seq2->seq.l >= opt->min_length && (seq2->seq.l <= opt->max_length || !opt->max_length)){
					buf[index].buf_index += sprintf(buf[index].buf + buf[index].buf_index, "@%s\n%s\n+\n%s\n@%s\n%s\n+\n%s\n", \
						seq1->name.s, seq1->seq.s, seq1->qual.s,
						seq2->name.s, seq2->seq.s, seq2->qual.s);
					total_base += seq2->seq.l;
				}else if (!fp2){
					buf[index].buf_index += sprintf(buf[index].buf + buf[index].buf_index, ">%s\n%s\n", seq1->name.s, seq1->seq.s);
				}
			}
			if (buf[index].buf_index > max_buf_index){
				i = get_fp_index(fp, opt->subfile);
				buf[index].fp = &fp[i];
				thpool_add_work(thpool_out, (void*)save_buf, &buf[index]);
				index = get_buf_index(buf, opt->thread);
			}
		}

		if (opt->max_base && total_base >= opt->max_base) break;
		fp1 = get_file_h(&opt->f1);
		fp2 = get_file_h(&opt->f2);
		seq1 = kseq_init(fp1);
		if (fp2) seq2 = kseq_init(fp2);
	}
	
	i = get_fp_index(fp, opt->subfile);
	buf[index].fp = &fp[i];
	thpool_add_work(thpool_out, (void*)save_buf, &buf[index]);

	thpool_wait_all(thpool_out);
	thpool_destroy(thpool_out);
	pthread_mutex_destroy(&buflock);
	pthread_mutex_destroy(&fplock);
	kseq_destroy(seq1);
	gzclose(fp1);
	if (fp2) {
		kseq_destroy(seq2);
		gzclose(fp2);
	}

	fprintf(stderr, "[STAT] Total used bases:%lu reads count(with N Base):%lu(%lu)\n",total_base, total_read, total_nread);
	if (total_nread > total_read/10) {
		fprintf(stderr, "Too many[%f] reads contains N base, please do QC first.\n", (double)total_nread/total_read);
		exit (1);
	}
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

int main(int argc, char *argv[])
{
	int c;
	int i = 1;
	opt_t opt;
	init_opt(&opt);
	while((c = getopt(argc, argv, "t:m:i:f:l:s:n:d:p:r:N")) != -1) {
	switch (c) {
		case 'i': i = atoi(optarg); break;
		case 'f': opt.min_length = mm_parse_num(optarg); break;
		case 'l': opt.max_length = mm_parse_num(optarg); break;
		case 's': opt.max_base = mm_parse_num(optarg); break;
		case 'n': opt.subfile = atoi(optarg); break;
		case 'r': opt.rmdup = atoi(optarg); break;
		case 'N': opt.rmn = 0; break;
		case 'p': strcpy(opt.outpre, optarg); break;
		case 't': opt.thread = atoi(optarg); assert (opt.thread >= 1); break;
		case 'm': {
			char *q;
			opt.max_mem = strtoll(optarg, &q, 0);
			if (*q == 'k' || *q == 'K') opt.max_mem <<= 10;
			else if (*q == 'm' || *q == 'M') opt.max_mem <<= 20;
			else if (*q == 'g' || *q == 'G') opt.max_mem <<= 30;
			opt.max_mem = opt.max_mem > 2 * MAX_READS ? opt.max_mem : MAX_READS * 2;
			break;
		}
		case 'd':
			if (access(optarg, 0) != 0) {
				if (mkdir(optarg, 0755) == -1) {
				fprintf(stderr, "Failed create output directory!\n");
				exit(1);
				}
			}
			strcpy(opt.outdir, optarg);
			break;
		default: return usage();
		}
	}
	if (optind + 1 > argc) return usage();

	init_opt_f(argv[optind], i, &opt);
	if (opt.f2.len) assert(opt.f1.len == opt.f2.len);

	sprintf(opt.out, "%s/%s", opt.outdir, opt.outpre);

	buffer *buf = init_buf(&opt);
	fout *fp = init_fp(&opt);
	split_data(&opt, buf, fp);
	destroy_opt(&opt);
	destroy_buf(&opt, buf);
	destroy_fp(&opt, fp);
	return 0;
}
