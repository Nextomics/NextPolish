#ifndef BASE_H
#define BASE_H

#include<stdlib.h>
#include "config.h"
#include "seqlist.h"

#define BASE_SHIFT 4

#define FLAG_ZERO 1
#define FLAG_COVERAGE 2
#define FLAG_DEPTH 4
#define FLAG_SNP  8
#define FLAG_THIRD 16
#define FLAG_INSERT 32
#define FLAG_LEFT 64
#define FLAG_RIGHT 128

#define FLAG_ZERO_N 254
#define FLAG_COVERAGE_N 253
#define FLAG_DEPTH_N 251
#define FLAG_SNP_N 247
#define FLAG_THIRD_N 239
#define FLAG_INSERT_N 223
#define FLAG_LEFT_N 191
#define FLAG_RIGHT_N 127

typedef struct {
	uint16_t kmer;
	uint16_t count;
} Kmer;

typedef struct {
	uint8_t base;
	uint16_t kmer;
	double score;
} Score;


typedef struct {
	uint8_t base;
	uint8_t flag;
	uint16_t refkmer;
	uint16_t count;
	SeqList* insert;
	SeqList* data;
	SeqList* score;
} Base;

Base* base_init(uint32_t size);

void base_destory(Base* base);

int32_t comparekmer(void* first, void* second);
void base_add_data(Base* base, uint16_t kmer);
void base_clean_data(Base* base);
double base_get_coverage(Base* base, uint16_t kmer);
int32_t base_get_nlargest(Base * base, Kmer ** maxn, int32_t n);
void base_merge_kmer(Base* base);

int32_t comparescore(void* first, void* second);
void base_add_score(Base* base, uint16_t kmer, double score);
Score* base_get_score(Base* base, uint16_t kmer);
void base_clean_score(Base* base);
Score* base_max_score(Base* base);

void base_set_flag(Base* base, uint8_t flag);

#endif


