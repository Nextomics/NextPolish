#include "snpvalid.h"

PolishResult * snp_valid(const char * tigname, Configure* configure)
{
	Contig* contig = contig_init(tigname, configure, sizeof(Kmer));
	
	contig->read_fliter = contig_read_fliter;
	SeqList* kmerregion = contig_get_region(contig, 0, contig->length - 1, contig->configure->min_len_inter_kmer, 0, FLAG_ZERO, contig_brim_with_extension);
	if (kmerregion->length > 0) {
		contig_merge_region(contig, kmerregion);
		contig_create_insert_region(contig, kmerregion, 0);
	}
	
	if (kmerregion->length > 0) {
		SeqList* nodepth = ss_spilt_region(contig, kmerregion, FLAG_ZERO, contig->configure->max_len_kmer);
		kmerregion->length = 0;
		ss_kmer_correct(contig, nodepth, kmerregion, 1);
		nodepth->length = 0;
		if (kmerregion->length > 0) {
			int32_t	*p = kmerregion->data;
			for (int i = 0; i < kmerregion->length; i += 2, p += 2) {
				fts_spilt_region(contig, *p, *(p + 1), FLAG_ZERO, nodepth);
			}
			ss_kmer_correct(contig, nodepth, NULL, 0);
		}
		seqlist_destory(nodepth);
	}
	seqlist_destory(kmerregion);
	
	PolishResult* result = contig_get_contig(contig, 0, contig->length - 1, 0);

	contig_destory(contig);

	return result;
}

inline void fts_spilt_region(Contig * contig, int32_t start, int32_t end, uint8_t flag, SeqList * result)
{
	int32_t i = start, j = 0, k, qstart = -1, qend = -1, count, mid;
	Base* p = contig->data[start];
	while (i < end || (i == end && j == 0)) {
		if ((p->flag&flag) == 0) {
			if (qstart == -1) {
				qstart = i;
			}
			qend = i;
		}
		else if (qstart != -1) {
			count = 2;
			if (qstart == start) {
				qend = start;
				count--;
			}
			mid = (qstart + qend) / 2;
			for (k = 0; k < count; k++) {
				seqlist_append(result, &mid);
				if (qstart != qend) {
					mid++;
				}
			}
			qstart = qend = -1;
		}
		p = contig_data_next(contig, &i, &j);
	}
	seqlist_append(result, &end);
}
