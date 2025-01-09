#ifndef SNP_PHASE_H
#define SNP_PHASE_H

#include <math.h>
#include "contig.h"
#include "kmercount.h"

typedef struct {
	int32_t pos;
	int32_t left;
	int32_t right;
	int16_t length;
	int16_t total;
	uint8_t* region[SNP_NUM];
	SeqList* link;
} Snps;

Snps* snps_init(int32_t right, int32_t length);
void snps_destory(Snps* snps);
//void snps_clean(Snps* snps, int32_t length);
void snps_resize(Snps* snps, int32_t length);
void snps_add_region(Snps* snps, uint8_t* region, int32_t i);
int32_t snps_get_index(Snps* snps, uint8_t* region);

typedef struct {
	Snps** data;
	int32_t length;
} SnpsList;

SnpsList* snpslist_init(SeqList* snps);
void snpslist_destory(SnpsList* snpslist);
int32_t snpslist_find(SnpsList * snpslist, int32_t pos);

PolishResult* snp_phase(const char* tigname, Configure* configure);

/*
*Summary: 查找[start,end]上的snp
*Parameters:
*	contig:
*	start:
*	end:
*Return SeqList<Snps>
*/
SeqList* ts_find_snps(Contig* contig, int32_t start, int32_t end);
int32_t ts_check_snps(Contig* contig, int32_t count, double rate, int32_t flag);

/*
*Summary: 过滤SeqList<Snps>中的snp
*Parameters:
*	contig:
*	snps: 通过ts_find_snps找到的snp
*Return snps即为过滤后的结果
*/
void ts_fliter_snps(Contig* contig, SeqList* snps);

/*
*Summary: 寻找两个snp之间的连接，结果通过ts_tranfer_link存储于base->data
*Parameters:
*	contig:
*	snps: 过滤后的snp
*Return
*/
void ts_find_snps_link(Contig* contig, SeqList* snps);
void ts_tranfer_link(Contig* contig, Snps** snp, KmerScore** ks, int16_t* total, int32_t flag);

/*
*Summary: 根据snp之间的连接计算各连接的得分
*Parameters:
*	contig:
*	snps: 过滤后的snp
*Return
*/
void ts_snps_score(Contig* contig, SeqList* snps);

/*
*Summary: 根据得分修正各snp
*Parameters:
*	contig:
*	snps: 过滤后的snp
*Return
*/
void ts_snps_correct(Contig* contig, SeqList* snps);


SeqList* ts_find_snp_region(Contig * contig, SeqList* snps, int32_t gap, uint32_t flag);
void ts_snps_parse_read(Contig * contig, bam1_t * read, int32_t start, int32_t end, uint32_t flag, SeqList* linkdata, KmerScore* ks);
void ts_snps_deal_linkdata(Contig* contig, SeqList* linkdata, SnpsList* snpslist, uint32_t flag);

/*
*Summary: 通过score chain修正低覆盖度区域
*Parameters:
*	contig:
*	nodepth: SeqList<int32_t>,每两个int32_t为一个低覆盖度区间
*Return
*/
void ts_correct_lower_depth(Contig* contig, SeqList* nodepth);
void ts_region_correct(Contig * contig, int32_t start, int32_t end);

int32_t ts_get_nlargest(SeqList * regiondata, KmerScore ** maxn, int32_t n);


#endif
