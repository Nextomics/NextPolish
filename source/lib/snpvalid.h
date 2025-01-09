#ifndef SNP_VALID_H
#define SNP_VALID_H

#include "contig.h"
#include "kmercount.h"

PolishResult* snp_valid(const char* tigname, Configure* configure);

void fts_spilt_region(Contig* contig, int32_t start, int32_t end, uint8_t flag, SeqList* result);

#endif

