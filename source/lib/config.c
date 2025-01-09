#define _GNU_SOURCE
#include "config.h"

const uint32_t opentrace[] = {
	2,2,2,2 //TRACE_STEP,TRACE_TIME,TRACE_POLISH,TRACE_ERROR
};

Configure * config_init(const char * fastafn, const char * bamfn, const char * thirdbamfn)
{
	Configure *result = calloc(sizeof(Configure), 1);
	result->trim_len_edge = 2;//read左右裁剪的碱基数
	result->ext_len_edge = 2;//区间边界向外延伸的大小
	result->min_map_quality = 0;
	result->indel_balance_factor_sgs = 0.5;//第一步和第二步使用二代数据计算得分时，coverage*indel_balance_factor_sgs
	result->min_count_ratio_skip = 0.8;//第一步中，选择的碱基占比<COVERAGE的变小写

	result->min_len_ldr = 3;//kmercount中区分无覆盖度区域与其他小写区域使用，若连续的小写碱基大于min_len_ldr,则认为是无覆盖度区域
	result->min_len_inter_kmer = 5;//kmercount中寻找kmer区间时，若两个小写碱基之间小于min_len_inter_kmer,则合并
	result->max_len_kmer = 50;//kmercount中kmer区间的最大长度
	result->max_count_kmer = 50;//kmercount中若mapq=60已经有SECOND_KEMR_COUNT个,则不再继续统计

	result->min_depth_snp = 3;//第三步中，对覆盖度<min_depth_snp的区间，使用三代数据进行得分链校正
	result->min_count_snp = 5;//第三步中，若snp的深度<min_count_snp，则添加三代数据进行判断是否是SNP
	result->min_count_snp_link = 5;//第三步中，若两snp之间的连接数<min_count_snp_link，则添加三代数据寻找连接
	result->ploidy = 2;//第三步中，snp coverage*ploidy
	result->indel_balance_factor_lgs = 0.33;//第三步中，使用三代数据进行得分链校正时，coverage*indel_balance_factor_lgs
	result->max_indel_factor_lgs = 0.21;//第三步中，使用三代数据进行得分链校正时，对于插入以及删除,若maxn[1]/maxn[0]>0.21则变小写
	result->max_snp_factor_lgs = 0.53;//第三步中，使用三代数据进行得分链校正时，对于非插入删除，若maxn[1]/maxn[0]>0.53则变小写
	result->min_snp_factor_sgs = 0.34;//第三步中，snp前两种类型的比率，maxn[1]/maxn[0]<min_snp_factor_sgs的判断为非snp点

	result->region_count = 10000;//默认申请内存区间
	result->count_read_ins_sgs = 10000;//统计插入片段大小时，读取的read数
	result->max_ins_len_sgs = 10000;//tlen的最大大小，绝对值>max_ins_len_sgs的不参与统计
	result->max_ins_fold_sgs = 5;//tlen*max_ins_fold_sgs用于过滤read
	result->max_variant_count_lgs = 150000;//三代数据read的最大长度

	result->max_clip_ratio_sgs = 0.15;//用于过滤二代read
	result->max_clip_ratio_lgs = 0.4;//用于过滤三代read

	result->trace_polish_open = 0;//是否打印改变点

	result->fastafn = strdup(fastafn);
	result->bamfn = access(bamfn, 0) ? NULL : strdup(bamfn);
	//result->faidx = fai_load(fastafn);
	if (result->bamfn) {
		result->read_tlen = bam_tlen(result, bamfn)*result->max_ins_fold_sgs;
	}
	else {
		result->read_tlen = 0;
	}
	result->thirdbamfn = access(thirdbamfn, 0) ? NULL : strdup(thirdbamfn);

	//result->idx = bam_load_idx(bamfn);
	//if (thirdbamfn) result->tidx = bam_load_idx(thirdbamfn);
	return result;
}

void config_destory(Configure * config)
{	
	if (config->fastafn) free(config->fastafn);
	if (config->bamfn) free(config->bamfn);
	if (config->thirdbamfn) free(config->thirdbamfn);
	//if (config->bamhdr) bam_hdr_destroy(config->bamhdr);
	//if(config->idx) hts_idx_destroy(config->idx);
	//if (config->tidx) hts_idx_destroy(config->tidx);
	//fai_destroy(config->faidx);
	free(config);
}

inline hts_idx_t * bam_load_idx(const char * bamfn)
{
	int len = strlen(bamfn);
	char* fnidx = calloc(sizeof(char), len + 5);
	sprintf(fnidx, "%s.bai", bamfn);
	hts_idx_t* idx = bam_index_load(fnidx);
	free(fnidx);
	return idx;
}

inline uint32_t bam_tlen(Configure* configure, const char* bamfn)
{
	BGZF* fp = bgzf_open(bamfn, "r");
	bam_hdr_t* bamhdr = bam_hdr_read(fp);
	bam1_t* read = bam_init1();
	uint32_t sum = 0, count = 1;
	configure->read_len = 0;
	while (bam_read1(fp, read) >= 0 && count < configure->count_read_ins_sgs)
	{
		if (read->core.isize > 0 && read->core.isize < configure->max_ins_len_sgs) {
			sum += read->core.isize;
			if (configure->read_len == 0) {
				configure->read_len = read->core.l_qseq;
			}
			count++;
		}
	}
	bam_destroy1(read);
	bam_hdr_destroy(bamhdr);
	bgzf_close(fp);
	return sum / count;
}


void print_log(FILE * __restrict __stream, char * format, va_list* ap)
{
	time_t t;
	struct tm *lt;
	time(&t);
	lt = localtime(&t);

	fprintf(__stream, "[ %02d-%02d-%02d %02d:%02d:%02d ] ", \
		lt->tm_year + 1900, lt->tm_mon + 1, lt->tm_mday, lt->tm_hour, \
		lt->tm_min, lt->tm_sec);
	vfprintf(__stream, format, *ap);
}

void trace_log(uint32_t trace, char * format, ...)
{
	if (opentrace[trace]) {
		va_list ap;
		va_start(ap, format);
		if (opentrace[trace] & 0x1) {
			print_log(stdout, format, &ap);
		}
		else {
			print_log(stderr, format, &ap);
		}
		va_end(ap);
	}
}
