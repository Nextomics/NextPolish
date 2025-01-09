#ifndef CONFIG_H
#define CONFIG_H

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <stdarg.h>
#include <ctype.h>
#include <unistd.h>
#include "htslib/faidx.h"
#include "htslib/bgzf.h"
#include "htslib/sam.h"

#define BASE_SHIFT 4
#define BASE_DEL 3
#define SNP_NUM 2
#define BASE_QUAL 41
#define READ_MAPQ 60
#define MAX_MAPQ 60

typedef struct {
	uint8_t trim_len_edge;
	uint8_t ext_len_edge;
	uint8_t min_map_quality;
	double indel_balance_factor_sgs;
	double min_count_ratio_skip;

	uint8_t min_len_ldr;
	uint8_t min_len_inter_kmer;
	uint8_t max_len_kmer;
	uint8_t max_count_kmer;

	uint8_t min_depth_snp;
	uint8_t min_count_snp;
	int8_t min_count_snp_link;
	double ploidy;
	double indel_balance_factor_lgs;
	double max_indel_factor_lgs;
	double max_snp_factor_lgs;
	double min_snp_factor_sgs;

	int32_t region_count;
	uint32_t count_read_ins_sgs;
	uint32_t max_ins_len_sgs;
	int32_t max_ins_fold_sgs;
	int32_t max_variant_count_lgs;

	double max_clip_ratio_sgs;
	double max_clip_ratio_lgs;
	
	int32_t trace_polish_open;
	int32_t read_tlen;
	int32_t read_len;
	
	char* fastafn;
	char* bamfn;
	char* thirdbamfn;

	/*faidx_t* faidx;
	hts_idx_t* idx;
	hts_idx_t* tidx;
	bam_hdr_t* bamhdr;*/
} Configure;

Configure* config_init(const char * fastafn, const char * bamfn, const char * thirdbamfn);
void config_destory(Configure* config);

hts_idx_t* bam_load_idx(const char* bamfn);
uint32_t bam_tlen(Configure* configure, const char* bamfn);

#ifndef DEBUG_OPEN
#define DEBUG_OPEN

#endif // !DEBUG_H

#define TRACE_STEP 0
#define TRACE_TIME 1
#define TRACE_POLISH 2
#define TRACE_ERROR 3

void print_log(FILE *__restrict __stream, char *format, va_list *ap);

void trace_log(uint32_t trace, char *format, ...);

#endif
