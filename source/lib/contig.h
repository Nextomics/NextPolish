#ifndef CONTIG_H
#define CONTIG_H

#include "config.h"
#include "base.h"

extern const uint8_t strtobase[];
extern const char basetostr[];

typedef struct {
	int32_t pos;
	int16_t index;
	char curbase;
	char base;
} PolishPoint;

typedef struct {
	char* contig;
	PolishPoint* data;
	int32_t length;
	int32_t datalength;
} PolishResult;

PolishResult* polishresult_init();
void polishresult_destory(PolishResult* polishresult);

typedef struct Contig1 {
	Base** data;
	int32_t length;
	const char* fastafn;
	const char* bamfn;
	const char* thirdbamfn;
	char* tigname;
	int32_t tid;
	int32_t inslength;
	uint8_t(*read_fliter)(struct Contig1*, bam1_t *);
	Configure* configure;

	faidx_t* faidx;
	htsFile* fp;
	htsFile* tfp;
	hts_idx_t* idx;
	hts_idx_t* tidx;
	bam1_t* read;
} Contig;

Contig* contig_init(const char* tigname, Configure* configure, uint32_t size);

Contig* contig_init_data(char* contig, uint32_t length, uint32_t size);

void contig_destory(Contig* contig);

void contig_create_insert(Contig* contig, htsFile* fp, hts_idx_t* idx, int32_t start, int32_t end, uint8_t flag);
void contig_create_insert_region(Contig* contig, SeqList* kmerregion, uint8_t insertflag);
void contig_parse_read_insert(Contig* contig, bam1_t* read, int32_t start, int32_t end, uint8_t flag);

/*
*Summary: 解析read,获取3-kmer数据，存储于base->data
*Parameters:
*	contig:
*	read: 读取的read
*	start: contig 开始位置
*	end: contig 结束位置
*	shift: 一般为4, 左移4位，3-kmer数据总共使用12bit
*Return
*/
void contig_parse_read(Contig* contig, bam1_t* read, int32_t start, int32_t end, uint8_t shift);

/*
*Summary: 根据配置裁剪read,返回裁剪后的开始和结束位置
*Parameters:
*	contig:
*	read:
*	qstart: 裁剪后的开始位置
*	qend: 裁剪后的结束位置
*Return
*/
void contig_cut_read(Contig* contig, bam1_t* read, int32_t* qstart, int32_t* qend);
uint16_t contig_left_kmer(uint16_t kmer, uint8_t base, uint8_t shift);
Base* contig_index(Contig* contig, int32_t i, uint16_t j);

/*
*Summary: contig作为一条read参与score chain计算
*Parameters:
*	contig:
*	start:
*	end:
*Return
*/
void contig_as_read(Contig* contig, int32_t start, int32_t end);

Base* contig_data_next(Contig* contig, int32_t* i, int32_t* j);
Base* contig_data_pre(Contig* contig, int32_t* i, int32_t* j);

/*
*Summary: 根据之前base得分计算当前base的得分，
*Parameters:
*	contig:
*	curbase: 当前base
*	lastbase: 前继base
*Return
*/
void contig_calculate_score(Contig* contig, Base* curbase, Base* lastbase, double rate);

/*
*Summary: 计算[start,end]上的所有base的得分，存储于base->score
*Parameters:
*	contig:
*	start:
*	end:
*Return
*/
void contig_region_score(Contig* contig, int32_t start, int32_t end, double rate);

/*
*Summary: 根据得分修正[start,end]上的碱基
*Parameters:
*	contig:
*	start:
*	end:
*Return
*/
void contig_region_correct(Contig* contig, int32_t start, int32_t end);

typedef void(*findbrim)(Contig* contig, uint8_t flag, int32_t bstart, int32_t bend, int32_t* start, int32_t* end);
void contig_brim_no_extension(Contig* contig, uint8_t flag, int32_t bstart, int32_t bend, int32_t* start, int32_t* end);
void contig_brim_with_extension(Contig* contig, uint8_t flag, int32_t bstart, int32_t bend, int32_t* start, int32_t* end);

/*
*Summary: 获取包含所有base->flag==flag的区间
*Parameters:
*	contig:
*	start:
*	end:
*	gap: 允许的最大间隔
*	con: 允许的最小连续数
*	flag:
*	brimfunc: 定界函数(contig_brim_no_extension/contig_brim_with_extension)
*Return SeqList<int32_t>,每两个int32_t为一个区间
*/
SeqList* contig_get_region(Contig* contig, int32_t start, int32_t end, uint16_t gap, uint16_t con, uint8_t flag, findbrim brimfunc);
void contig_merge_region(Contig* contig, SeqList* seqlist);

void contig_clean_region(Contig* contig, int32_t start, int32_t end);

/*
*Summary: 获取read的裁剪率，(H+S)/read length
*Parameters:
*	read:
*Return 裁剪率
*/
double contig_read_cliprate(bam1_t* read);

/*
*Summary: 过滤read，kmercount,snpphase,snpvalid中二代数据使用此过滤条件
*Parameters:
*	read:
*Return 0,1,2为过滤级别
*/
uint8_t contig_read_fliter(Contig* contig, bam1_t* read);

/*
*Summary: 过滤read，scorechain中二代数据使用此过滤条件
*Parameters:
*	read:
*Return 0,1为过滤级别
*/
uint8_t contig_read_fliter1(Contig* contig, bam1_t* read);

/*
*Summary: 过滤read，三代数据使用此过滤条件
*Parameters:
*	read:
*Return 0,1为过滤级别
*/
uint8_t contig_read_fliter2(Contig* contig, bam1_t* read);

void contig_parse_region(Contig* contig, htsFile* fp, int32_t start, int32_t end, uint8_t filterlevel, uint8_t shift);
void contig_score_correct(Contig* contig, int32_t start, int32_t end, int32_t flag, double rate);

/*
*Summary: 获取contig上[start,end]的值，base->flag==flag的base变小写
*Parameters:
*	contig:
*	start:
*	end:
*	flag: 需要变小写的flag
*Return 返回修正后的contig以及改变的base
*/
PolishResult* contig_get_contig(Contig* contig, int32_t start, int32_t end, uint8_t flag);

int32_t contig_get_length(Contig* contig, int32_t start, int32_t end);
void contig_update_contig(Contig* contig, int32_t start, int32_t end, uint8_t* region, uint16_t index);
void contig_update_flag(Contig* contig, int32_t start, int32_t end, uint16_t flag);
void contig_clean_flag(Contig* contig, int32_t start, int32_t end, uint16_t flag);

hts_itr_t* contig_update_iter(Contig* contig, BGZF* fp, hts_idx_t* idx, hts_itr_t* iter, int32_t start, int32_t end, int32_t* iterend);
int contig_next_iter(Contig * contig, htsFile* fp, hts_idx_t* idx, hts_itr_t ** iter, int32_t start, int32_t end, int32_t * iterend, uint64_t* curr_off, bam1_t** read, int32_t* curr_end, int32_t* nextposend, int32_t flag);

void contig_write_to_file(const char* writefasta, const char* tigname, char * seq, int32_t length, int32_t tag);

typedef PolishResult*(*contig_run)(const char * tigname, Configure* configure);
void contig_total(contig_run func, int32_t step, Configure* configure, int32_t tag);



char* fasta_fetch(faidx_t* faidx, const char* tigname, int* length);

void contig_swap_iter(hts_itr_t* iter);

#endif
