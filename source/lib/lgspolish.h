#ifndef TD_SCORE_CHAIN_H
#define TD_SCORE_CHAIN_H

#include "contig.h"

#define KMER_K 3

typedef struct {
	uint8_t base;
	int16_t index;
	int32_t pos;
} TdKmerItem;

typedef struct {
	TdKmerItem kmer[KMER_K];
	uint16_t count;
} TdKmer;

void tdkmer_add_item(TdKmer *tdkmer, uint8_t base, int32_t pos, int16_t index);



PolishResult* lgspolish(const char* tigname, Configure* configure);

void td_region_score_chain(Contig * contig, int32_t start, int32_t end, double rate);

int32_t comparetdkmer(void * first, void * second);
void td_base_add_data(Base * base, TdKmer* kmer);
void td_contig_as_read(Contig * contig, int32_t start, int32_t end);

void td_parse_read(Contig * contig, bam1_t * read, int32_t start, int32_t end);
void td_parse_region(Contig * contig, htsFile* fp, int32_t start, int32_t end, uint8_t filterlevel);

void td_base_add_score(Base * base, TdKmer* tdkmer, double score, uint16_t index);
void td_region_score(Contig * contig, int32_t start, int32_t end, double rate);

void td_region_correct(Contig * contig, int32_t start, int32_t end);

#endif
