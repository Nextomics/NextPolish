#ifndef KMER_COUNT_H
#define KMER_COUNT_H

#include "contig.h"

typedef struct {
	uint8_t* region;
	int32_t length;
	int32_t qual;
	int32_t mapqual;
	int32_t num;
} KmerScore;

KmerScore* ks_init(int32_t length);
void ks_destory(KmerScore* ks);
void ks_clean(KmerScore* ks, int32_t length);
void ks_append(KmerScore* ks, uint8_t base);
int32_t ks_compare_region(void* firstks, void* secondks);
int32_t ks_compare_length(void* firstks, void* secondks);
int32_t ks_compare_mapqual(void* firstks, void* secondks);
int32_t ks_compare(void* firstks, void* secondks);


PolishResult* kmer_count(const char* tigname, Configure* configure);

SeqList* ss_spilt_region(Contig* contig, SeqList* seqlist, uint8_t flag, uint8_t max);

void ss_kmer_correct(Contig* contig, SeqList* region, SeqList* result, int32_t flagzero);

/*
*Summary: 解析read,检查kmer并存储于regiondata
*Parameters:
*	contig:
*	read:
*	start:
*	end:
*	length: kmer长度
*	regiondata: 存储kmer
*	ks: 存储获取到的kmer的相关信息
*	count: 统计kmer数量，不需要则为NULL
*	left: 需要确认的左边界值，不需要则为-1
*	right: 需要确认的右边界值，不需要则为-1
*	flag: 默认为0
*	flagzero: 1表示当有覆盖时取消base的FLAG_ZERO标记，0则不取消
*Return ks用于循环使用
*/
KmerScore* ss_kmer_get_region(Contig* contig, bam1_t* read, int32_t start, int32_t end, int32_t length, SeqList* regiondata, KmerScore* ks, int32_t* count, int32_t left, int32_t right, int32_t flag, int32_t flagzero);

/*
*Summary: 解析read,获取[start,end]的kmer
*Parameters:
*	contig:
*	read:
*	start:
*	end:
*	ks: 存储获取到的kmer的相关信息
*	left: 需要确认的左边界值，不需要则为-1
*	right: 需要确认的右边界值，不需要则为-1
*	flagzero: 1表示当有覆盖时取消base的FLAG_ZERO标记，0则不取消
*Return 被确认的left和right,0表示read与contig在left和right两个位置的碱基都不相同，2表示在两个位置都相同，1表示有一个位置相同
*/
int32_t ss_parse_read_kmer(Contig* contig, bam1_t* read, int32_t start, int32_t end, KmerScore* ks, int32_t left, int32_t right, int32_t flagzero);



#endif
