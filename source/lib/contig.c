#include "contig.h"


/*#include <strings.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <zlib.h>
#include <assert.h>
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "cram/cram.h"
#include "hts_internal.h"
#include "htslib/hfile.h"
#include "htslib/hts_endian.h"

#include "htslib/khash.h"*/

PolishResult * polishresult_init()
{
	return calloc(sizeof(PolishResult), 1);
}

void polishresult_destory(PolishResult * polishresult)
{
	if (polishresult->contig) free(polishresult->contig);
	if (polishresult->data) free(polishresult->data);
	free(polishresult);
}

Contig * contig_init(const char* tigname, Configure* configure, uint32_t size)
{
	int32_t length = 0;
	faidx_t* faidx = fai_load(configure->fastafn);
	char* seqbase = fasta_fetch(faidx, tigname, &length);
	Contig* contig = contig_init_data(seqbase, length, size);
	contig->faidx = faidx;

	contig->configure = configure;
	contig->idx = NULL;
	contig->fp = NULL;
	contig->tidx = NULL;
	contig->tfp = NULL;
	if (configure->bamfn) {
		contig->idx = bam_load_idx(configure->bamfn);
		contig->fp = hts_open(configure->bamfn, "rb");
	}

	if (configure->thirdbamfn) {
		contig->tfp = hts_open(configure->thirdbamfn, "rb");
		contig->tidx = bam_load_idx(configure->thirdbamfn);
	}
	
	contig->read = bam_init1();

	length = strlen(tigname);
	contig->tigname = calloc(sizeof(char*), length + 1);
	strncpy(contig->tigname, tigname, length);

	contig->fastafn = configure->fastafn;
	contig->bamfn = configure->bamfn;
	contig->thirdbamfn = configure->thirdbamfn;
	
	bam_hdr_t* bamhdr = NULL;
	if (contig->fp) {
		bamhdr = bam_hdr_read(contig->fp->fp.bgzf);
	}
	else {
		bamhdr = bam_hdr_read(contig->tfp->fp.bgzf);
	}
	contig->tid = bam_name2id(bamhdr, tigname);
	contig->inslength = 0;
	//contig->read_fliter = contig_read_fliter;

	bam_hdr_destroy(bamhdr);
	free(seqbase);
	return contig;
}

inline Contig * contig_init_data(char* contig, uint32_t length, uint32_t size)
{
	Contig* result = (Contig*)calloc(sizeof(Contig),1);
	result->data = (Base**)malloc(sizeof(Base*)*length);
	//memset(result->data, 0, sizeof(Base*)*length);
	if (result->data == NULL) {
		fprintf(stderr, "contig data malloc failed,the memory is full!\n");
		exit(1);
	}
	Base** p = result->data;
	char* q = contig;
	for (uint32_t i = 0; i < length; i++, p++, q++) {
		*p = base_init(size);
		if (*q >= 97 && *q <= 122) {
			*q -= 32;
			(*p)->flag |= FLAG_ZERO;
		}
		(*p)->base = strtobase[(uint8_t)*q];
	}
	result->length = length;
	return result;
}

void contig_destory(Contig * contig)
{
	if (contig != NULL) {
		if (contig->data != NULL) {
			Base** p = contig->data;
			for (uint32_t i = 0; i < contig->length; i++) {
				base_destory(*p++);
			}
			free(contig->data);
		}
		if (contig->fp) hts_close(contig->fp);
		if (contig->tfp) hts_close(contig->tfp);
		if (contig->read) bam_destroy1(contig->read);
		if (contig->idx) hts_idx_destroy(contig->idx);
		if (contig->tidx) hts_idx_destroy(contig->tidx);
		fai_destroy(contig->faidx);
		free(contig->tigname);
		free(contig);
	}
}

/*static inline int hts_itr_next1(BGZF *fp, hts_itr_t *iter, void *r, void *data)
{
	int ret, tid, beg, end;
	if (iter == NULL || iter->finished) return -1;
	if (iter->read_rest) {
		if (iter->curr_off) { // seek to the start
			if (bgzf_seek(fp, iter->curr_off, SEEK_SET) < 0) return -1;
			iter->curr_off = 0; // only seek once
		}
		ret = iter->readrec(fp, data, r, &tid, &beg, &end);
		if (ret < 0) iter->finished = 1;
		iter->curr_tid = tid;
		iter->curr_beg = beg;
		iter->curr_end = end;
		return ret;
	}
	// A NULL iter->off should always be accompanied by iter->finished.
	assert(iter->off != NULL);
	for (;;) {
		if (iter->curr_off == 0 || iter->curr_off >= iter->off[iter->i].v) { // then jump to the next chunk
			if (iter->i == iter->n_off - 1) { ret = -1; break; } // no more chunks
			if (iter->i < 0 || iter->off[iter->i].v != iter->off[iter->i + 1].u) { // not adjacent chunks; then seek
				if (bgzf_seek(fp, iter->off[iter->i + 1].u, SEEK_SET) < 0) return -1;
				iter->curr_off = bgzf_tell(fp);
			}
			++iter->i;
		}
		if ((ret = iter->readrec(fp, data, r, &tid, &beg, &end)) >= 0) {
			iter->curr_off = bgzf_tell(fp);
			if (tid != iter->tid || beg >= iter->end) { // no need to proceed
				ret = -1; break;
			}
			else if (end > iter->beg && iter->end > beg) {
				iter->curr_tid = tid;
				iter->curr_beg = beg;
				iter->curr_end = end;
				return ret;
			}
		}
		else break; // end of file or error
	}
	iter->finished = 1;
	return ret;
}*/

void contig_create_insert(Contig * contig, htsFile* fp, hts_idx_t* idx, int32_t start, int32_t end, uint8_t flag)
{
	hts_itr_t* iter = bam_itr_queryi(idx, contig->tid, start, end + 1);
	bam1_t* read = contig->read;
	while (sam_itr_next(fp, iter, read) >= 0) {
		if (contig->read_fliter(contig, read) >= 1) {
			contig_parse_read_insert(contig, read, start, end, flag);
		}
	}
	hts_itr_destroy(iter);
}

void contig_create_insert_region(Contig * contig, SeqList * kmerregion, uint8_t insertflag)
{
	int32_t *p = kmerregion->data, iterend = 0, nextposend, curr_end = 0, flag;
	hts_itr_t *iter = NULL;
	uint64_t curr_off = 0;
	bam1_t *read = contig->read;
	contig->read_fliter = contig_read_fliter;
	for (int i = 0; i < kmerregion->length; i += 2, p += 2) {
		nextposend = i + 2 < kmerregion->length ? *(p + 2) : -1;
		flag = 2;
		while (contig_next_iter(contig, contig->fp, contig->idx, &iter, *p, *(p + 1), &iterend, &curr_off, &read, &curr_end, &nextposend, flag) >= 0) {
			if (contig->read_fliter(contig, read) >= 1) {
				contig_parse_read_insert(contig, read, *p, *(p + 1), insertflag);
			}
			flag = 0;
		}
	}
	if (iter) hts_itr_destroy(iter);
}

inline void contig_parse_read_insert(Contig * contig, bam1_t * read, int32_t start, int32_t end, uint8_t flag)
{
	if (read->core.n_cigar) {
		int32_t pos = read->core.pos, i = 0, j = 0, len;
		uint8_t flag1;
		SeqList* seqlist = NULL;
		Base* base;
		uint32_t *cigar = bam_get_cigar(read);
		for (i = 0; i < read->core.n_cigar; ++i) {
			switch (bam_cigar_op(cigar[i])) {
			case BAM_CMATCH:
			case BAM_CDEL:
				pos += bam_cigar_oplen(cigar[i]);
				break;
			case BAM_CINS:
				if (pos > start && pos <= end) {
					if (flag == 0 || contig->data[pos - 1]->flag&flag) {
						len = bam_cigar_oplen(cigar[i]);
						if (contig->data[pos - 1]->insert == NULL) {
							contig->data[pos - 1]->insert = seqlist_init(sizeof(Base), len);
						}
						seqlist = contig->data[pos - 1]->insert;
						if (seqlist->maxsize < len) {
							seqlist->data = realloc(seqlist->data, len*seqlist->stepsize);
							if (seqlist->data == NULL) {
								fprintf(stderr, "seqlist malloc failed,the memory is full!\n");
								exit(1);
							}
							seqlist->maxsize = len;
						}
						flag1 = contig->data[pos - 1]->flag;
						for (j = seqlist->length; j < len; j++, contig->inslength++) {
							base = base_init(sizeof(Kmer));
							base->flag = flag1;
							seqlist_append(seqlist, base);
							free(base);
						};
					}
				}
				break;
			}
		}
	}
}

void contig_parse_read(Contig * contig, bam1_t * read, int32_t start, int32_t end, uint8_t shift)
{
	if (read->core.n_cigar) {
		/*char* name = bam_get_qname(read);
		char* name1 = "ST-E00600:249:H7NKLALXX:1:1374:8983:11757";
		if (strcmp(name, name1) == 0) {
			int a = 1;
		}*/
		uint16_t kmer = 0x0;
		int32_t pos = read->core.pos, qpos = 0, qstart, qend, i, j, k, len;
		uint8_t curcigar, lastcigar = BAM_CINS;
		uint32_t *cigar = bam_get_cigar(read);
		uint8_t* seq = bam_get_seq(read);
		SeqList* p;
		contig_cut_read(contig, read, &qstart, &qend);
		for (i = 0; i < read->core.n_cigar; ++i) {
			len = bam_cigar_oplen(cigar[i]);
			curcigar = bam_cigar_op(cigar[i]);
			switch (curcigar) {
			case BAM_CMATCH:
			case BAM_CDEL:
				for (j = 0; j < len; j++, pos++) {
					if (pos >= start && pos <= end && qpos >= qstart && qpos <= qend) {
						/*if (pos == 8) {
							int a = 1;
						}*/
						if (lastcigar != BAM_CINS && pos > start && (qpos > qstart || (qpos == qstart && lastcigar == BAM_CDEL))) {
							p = contig->data[pos - 1]->insert;
							if (p != NULL) {
								for (k = 0; k < p->length; k++) {
									kmer = contig_left_kmer(kmer, BASE_DEL, shift);
									/*if (pos == 9648 && kmer == 2083) {
										int a = 1;
									}*/
									base_add_data(seqlist_index(p, k), kmer);
								}
							}
						}
						if (curcigar == BAM_CDEL) {
							kmer = contig_left_kmer(kmer, BASE_DEL, shift);
						}
						else {
							kmer = contig_left_kmer(kmer, bam_seqi(seq, qpos), shift);
						}
						base_add_data(contig->data[pos], kmer);
					}
					if (curcigar != BAM_CDEL) {
						qpos++;
					}
					lastcigar = curcigar;
				}
				break;
			case BAM_CINS:
				if (pos){
					p = contig->data[pos - 1]->insert;
					for (j = 0; j < len; j++, qpos++) {
						if (pos > start && pos <= end && qpos >= qstart && qpos <= qend) {
							kmer = contig_left_kmer(kmer, bam_seqi(seq, qpos), shift);
							base_add_data(seqlist_index(p, j), kmer);
						}
					}
					if (pos > start && pos <= end && qpos > qstart && qpos <= qend + 1) {
						for (; j < p->length; j++) {
							kmer = contig_left_kmer(kmer, BASE_DEL, shift);
							base_add_data(seqlist_index(p, j), kmer);
						}
					}
					lastcigar = curcigar;
				}else{
					qpos += len;
					qstart += len;
					lastcigar = curcigar;
				}
				break;
			case BAM_CHARD_CLIP:
			case BAM_CSOFT_CLIP:
				qpos += len;
				break;
			}
			if (pos > end) {
				break;
			}
		}
	}
}

inline void contig_cut_read(Contig* contig, bam1_t * read, int32_t * qstart, int32_t * qend)
{
	uint32_t *cigar = bam_get_cigar(read);
	int32_t len = bam_cigar_oplen(cigar[0]), addlen = 0;
	uint8_t curcigar = bam_cigar_op(cigar[0]);
	if (curcigar == BAM_CSOFT_CLIP) {
		addlen = len;
	}
	*qstart = contig->configure->trim_len_edge + addlen;
	len = bam_cigar_oplen(cigar[read->core.n_cigar - 1]);
	curcigar = bam_cigar_op(cigar[read->core.n_cigar - 1]);
	addlen = 0;
	if (curcigar == BAM_CSOFT_CLIP) {
		addlen = len;
	}
	*qend = read->core.l_qseq - contig->configure->trim_len_edge - addlen - 1;
	uint8_t* seq = bam_get_seq(read);
	if (contig->configure->trim_len_edge > 0) {
		while (bam_seqi(seq, *qstart) == bam_seqi(seq, (*qstart) - 1)) {
			(*qstart)++;
		}
		while (bam_seqi(seq, *qend) == bam_seqi(seq, (*qend) + 1)) {
			(*qend)--;
		}
	}
}

inline uint16_t contig_left_kmer(uint16_t kmer, uint8_t base, uint8_t shift)
{
	return (kmer & 0xff) << shift | base;
}

inline Base * contig_index(Contig * contig, int32_t i, uint16_t j)
{
	if (j == 0) {
		return contig->data[i];
	}
	return seqlist_index(contig->data[i]->insert, j - 1);
}

void contig_as_read(Contig * contig, int32_t start, int32_t end)
{
	uint16_t kmer = 0x0;
	int32_t i = start, j = 0;
	Base* p = contig->data[start];
	while (i < end || (i == end && j == 0)) {
		p->refkmer = kmer = contig_left_kmer(kmer, p->base, BASE_SHIFT);
		base_add_data(p, kmer);
		p = contig_data_next(contig, &i, &j);
	}
}

inline Base * contig_data_next(Contig * contig, int32_t* i, int32_t* j)
{
	if (*i + 1 >= contig->length) {
		*i = contig->length;
		return contig->data[*i - 1];
	}
	Base* p = contig->data[*i];
	if (p->insert == NULL || p->insert->length == *j) {
		(*i)++;
		*j = 0;
		return contig->data[*i];
	}
	(*j)++;
	return seqlist_index(p->insert, *j - 1);
}

inline Base * contig_data_pre(Contig * contig, int32_t* i, int32_t* j)
{
	if (*i - 1 < 0) {
		*i = -1;
		return contig->data[0];
	}
	Base* p = contig->data[*i];
	if (*j == 0) {
		(*i)--;
		p = contig->data[*i];
		if (p->insert != NULL) {
			*j = p->insert->length;
		}
	}
	else {
		(*j)--;
	}
	if (*j == 0) {
		return p;
	}
	return seqlist_index(p->insert, *j - 1);
}

inline void contig_calculate_score(Contig* contig, Base * curbase, Base * lastbase, double rate)
{
	base_clean_score(curbase);
	Kmer* p = curbase->data->data;
	double score = 0;
	Score* q = NULL;
	uint16_t temp, count, total = curbase->count;
	if (total > 1) {
		total--;
	}
	for (int i = 0; i < curbase->data->length; i++, p++) {
		//score = base_get_score(lastbase, p->kmer >> BASE_SHIFT)->score;
		temp = p->kmer >> BASE_SHIFT;
		if ((temp & 0xf) == 0) {
			score = base_max_score(lastbase)->score;
		}
		else {
			score = base_get_score(lastbase, temp)->score;
		}
		count = p->count;

		if (p->kmer == curbase->refkmer && curbase->count > 1) {
			count--;
		}
		score += count - total * rate;
		q = base_get_score(curbase, p->kmer);
		if (q == NULL || q->score < score) {
			base_add_score(curbase, p->kmer, score);
		}
	}
}

void contig_region_score(Contig * contig, int32_t start, int32_t end, double rate)
{
	int32_t i = start, j = 0;
	Base *temp = base_init(sizeof(Kmer));
	Base *p = temp, *q = contig->data[start];
	Kmer* t = q->data->data;
	for (int i = 0; i < q->data->length; i++, t++) {
		base_add_score(temp, t->kmer >> BASE_SHIFT, 0);
	}
	while (i < end || (i == end && j == 0)) {
		contig_calculate_score(contig, q, p, rate);
		p = q;
		q = contig_data_next(contig, &i, &j);
	}
	base_destory(temp);
}

void contig_region_correct(Contig * contig, int32_t start, int32_t end)
{
	Base* base = contig->data[end];
	Score* score = base_max_score(base);
	int32_t i = end, j = 0;
	while (i > start || (i == start && j == 0))
	{
		base->base = score->base;
		if (base->count == 1) {
			base->flag |= FLAG_ZERO;
		}
		else {
			base->flag &= FLAG_ZERO_N;
		}
		if (base_get_coverage(base, base->base) < contig->configure->min_count_ratio_skip) {
			base->flag |= FLAG_COVERAGE;
		}
		else {
			base->flag &= FLAG_COVERAGE_N;
		}
		base = contig_data_pre(contig, &i, &j);
		score = base_get_score(base, score->kmer >> BASE_SHIFT);
	}
}

void contig_brim_no_extension(Contig* contig, uint8_t flag, int32_t bstart, int32_t bend, int32_t * start, int32_t * end)
{
	*start = *start >= bstart + contig->configure->ext_len_edge ? *start - contig->configure->ext_len_edge : bstart;
	*end = *end <= bend - contig->configure->ext_len_edge ? *end + contig->configure->ext_len_edge : bend;
}

void contig_brim_with_extension(Contig* contig, uint8_t flag, int32_t bstart, int32_t bend, int32_t * start, int32_t * end)
{
	contig_brim_no_extension(contig, flag, bstart, bend, start, end);
	Base **p = &contig->data[*start + 1];
	while (*start > bstart && ((*p)->base == (*(p - 1))->base || ((*(p - 1))->flag&flag) != 0)) {
		(*start)--;
		p--;
	}
	p = &contig->data[*end - 1];
	while (*end < bend && ((*p)->base == (*(p + 1))->base || ((*(p + 1))->flag&flag) != 0)) {
		(*end)++;
		p++;
	}
}

SeqList * contig_get_region(Contig * contig, int32_t start, int32_t end, uint16_t gap, uint16_t con, uint8_t flag, findbrim brimfunc)
{
	SeqList* result = seqlist_init(sizeof(int32_t), contig->configure->region_count);
	int32_t i = start, j = 0, qstart = -1, qend = -1;
	uint16_t pgap = 0, pcon = 0;
	Base* p = contig->data[start];
	while (i < end || (i == end && j == 0)) {
		if ((p->flag&flag) != 0) {
			if (qstart == -1) {
				qstart = i;
				pcon = 1;
			}
			else if (pgap == 0) {
				pcon++;
			}
			else {
				pcon = 1;
			}
			pgap = 0;
			qend = i;
		}
		else if (qstart != -1) {
			pgap++;
			if (pgap > gap) {
				if (pcon > con) {
					brimfunc(contig, flag, start, end, &qstart, &qend);
					seqlist_append(result, &qstart);
					seqlist_append(result, &qend);
					if (qend > i) {
						i = qend;
						j = 0;
					}
				}
				qstart = qend = -1;
			}
		}
		p = contig_data_next(contig, &i, &j);
	}
	if (qstart != -1) {
		brimfunc(contig, flag, start, end, &qstart, &qend);
		seqlist_append(result, &qstart);
		seqlist_append(result, &qend);
	}
	return result;
}

/*SeqList * contig_get_region(Contig * contig, int32_t start, int32_t end,uint16_t gap, uint8_t flag,findbrim brimfunc)
{
	SeqList* result = seqlist_init(sizeof(int32_t), contig->regioncount);
	int32_t i = start, j = 0, qstart = -1, qend = -1;
	Base* p = contig->data[start];
	while(i < end || i == end && j == 0) {
		if ((p->flag&flag)!=0) {
			if (qstart == -1) {
				qstart = i;
			}
			qend = i;
		}
		else if (qstart != -1) {
			brimfunc(contig, flag, start, end, &qstart, &qend);
			seqlist_append(result, &qstart);
			seqlist_append(result, &qend);
			i = qend;
			j = 0;
			qstart = qend = -1;
		}
		p = contig_data_next(contig, &i, &j);
	}
	if (qstart != -1) {
		brimfunc(contig, flag, start, end, &qstart, &qend);
		seqlist_append(result, &qstart);
		seqlist_append(result, &qend);
	}
	return result;
}*/

void contig_merge_region(Contig * contig, SeqList * seqlist)
{
	if (seqlist->length == 0) {
		return;
	}
	int32_t *pstart = seqlist->data, *pend = pstart + 1, *qstart = pstart, *qend = pend, length = 2;
	for (int i = 0; i < seqlist->length; i += 2) {
		if (*pstart >= *qend) {
			qstart += 2;
			qend = qstart + 1;
			if (qstart != pstart) memcpy(qstart, pstart, seqlist->stepsize);
			if (qend != pend) memcpy(qend, pend, seqlist->stepsize);
			length += 2;
		}
		else {
			while (*pstart < *qstart) {
				qstart -= 2;
			}
			qend = qstart + 1;
			*qend = *pend;
		}
		pstart += 2;
		pend = pstart + 1;
	}
	seqlist->length = length;
}

void contig_clean_region(Contig * contig, int32_t start, int32_t end)
{
	int32_t i = start, j = 0;
	Base *p = contig->data[start];
	while (i < end || (i == end && j == 0)) {
		base_clean_data(p);
		p = contig_data_next(contig, &i, &j);
	}
}

inline double contig_read_cliprate(bam1_t * read)
{
	uint32_t *cigar = bam_get_cigar(read);
	int32_t len = bam_cigar_oplen(cigar[0]), addlen = 0;
	uint8_t curcigar = bam_cigar_op(cigar[0]);
	if (curcigar == BAM_CSOFT_CLIP) {
		addlen += len;
	}
	len = bam_cigar_oplen(cigar[read->core.n_cigar - 1]);
	curcigar = bam_cigar_op(cigar[read->core.n_cigar - 1]);
	if (curcigar == BAM_CSOFT_CLIP) {
		addlen += len;
	}
	return read->core.l_qseq > 0 ? addlen / (double)read->core.l_qseq : 0;
}

uint8_t contig_read_fliter(Contig* contig, bam1_t * read)
{
	uint8_t result = 0;
	if ((read->core.flag & 0xC04) == 0) {
		int32_t length = read->core.isize >= 0 ? read->core.isize : -read->core.isize;
		double cliprate = contig_read_cliprate(read);
		if ((length > 0 && length < contig->configure->read_tlen) || cliprate < contig->configure->max_clip_ratio_sgs) {
			result = 1;
			/*if (cliprate < contig->configure->max_clip_ratio_sgs+0.05) {
				result = 2;
			}*/
			if (read->core.qual >= contig->configure->min_map_quality && (cliprate < contig->configure->max_clip_ratio_sgs+0.05)) {
				result = 2;
			}
		}
	}
	return result;
}

uint8_t contig_read_fliter1(Contig * contig, bam1_t * read)
{
	uint8_t result = 0;
	if ((read->core.flag & 0xC04) == 0) {
		result = 1;
		/*if (read->core.qual >= 10) {
			result = 2;
		}*/
	}
	return result;
}

uint8_t contig_read_fliter2(Contig * contig, bam1_t * read)
{
	uint8_t result = 0;
	if ((read->core.flag & 0xD04) == 0 && contig_read_cliprate(read) <= contig->configure->max_clip_ratio_lgs) {
		result = 1;
	}
	return result;
}

void contig_parse_region(Contig * contig, htsFile* fp, int32_t start, int32_t end, uint8_t filterlevel, uint8_t shift)
{
	//BGZF* fp = bgzf_open(bamfn, "r");
	//hts_idx_t* idx = bam_load_idx(bamfn);
	hts_itr_t* iter = bam_itr_queryi(contig->idx, contig->tid, start, end + 1);
	bam1_t* read = contig->read;
	while (sam_itr_next(fp, iter, read) >= 0)
	{
		if (contig->read_fliter(contig, read) == filterlevel) {
			contig_parse_read(contig, read, start, end, shift);
		}
	}
	//bam_destroy1(read);
	hts_itr_destroy(iter);
	//hts_idx_destroy(idx);
	//bgzf_close(fp);
}

void contig_score_correct(Contig * contig, int32_t start, int32_t end, int32_t flag, double rate)
{
	int32_t filterlevel = flag & 0xf, insert = (flag >> 4) & 0xf;

	if ((insert & 0x1) == 0) {
		contig_create_insert(contig, contig->fp, contig->idx, start, end, 0);
	}

	contig_as_read(contig, start, end);

	contig_parse_region(contig, contig->fp, start, end, filterlevel, BASE_SHIFT);

	contig_region_score(contig, start, end, rate);
	contig_region_correct(contig, start, end);

	if (filterlevel == 2) {
		SeqList* nodepth = contig_get_region(contig, start, end, 0, 0, 1, contig_brim_no_extension);
		if (nodepth->length != 0) {
			contig_merge_region(contig, nodepth);
			int32_t* p = nodepth->data;
			for (int i = 0; i < nodepth->length; i += 2, p += 2) {
				contig_parse_region(contig, contig->fp, *p, *(p + 1), 1, BASE_SHIFT);
				contig_region_score(contig, *p, *(p + 1),contig->configure->indel_balance_factor_sgs);
				contig_region_correct(contig, *p, *(p + 1));
			}
		}
		seqlist_destory(nodepth);
	}
}

PolishResult * contig_get_contig(Contig * contig, int32_t start, int32_t end, uint8_t flag)
{

	int32_t i = start, j = 0, length = 0, length1 = 0;
	char* result = calloc(sizeof(char), contig->length + contig->inslength + 1), *q = result, *seqbase = NULL;
	SeqList* data = NULL;
	PolishPoint temp;
	if (contig->configure->trace_polish_open) {
		data = seqlist_init(sizeof(PolishPoint), contig->length / 100);
		seqbase = fasta_fetch(contig->faidx, contig->tigname, &length1);
	}
	if (result != NULL) {
		Base* p = contig->data[start];
		uint8_t sign = 0;
		while (i < end || (i == end && j == 0)) {
			if (p->base == 3) {
				if ((p->flag & flag) != 0) {
					sign = 1;
				}
				if (contig->configure->trace_polish_open&&j == 0) {
					temp.pos = i;
					temp.index = j;
					temp.curbase = '.';
					temp.base = toupper(seqbase[i]);
					seqlist_append(data, &temp);
				}
			}
			else {
				*q = basetostr[p->base];
				if (contig->configure->trace_polish_open) {
					temp.pos = i;
					temp.index = j;
					temp.curbase = *q;
					if (j != 0) {
						temp.base = '.';
						seqlist_append(data, &temp);
					}
					else if (*q != toupper(seqbase[i])) {
						temp.base = toupper(seqbase[i]);
						seqlist_append(data, &temp);
					}
				}
				if (sign || (p->flag & flag) != 0) {
					*q += 32;
					sign = 0;
				}
				q++;
				length++;
			}
			p = contig_data_next(contig, &i, &j);
		}
	}
	result[length] = '\0';
	PolishResult* polishresult = polishresult_init();
	polishresult->contig = result;
	polishresult->length = length;
	if (contig->configure->trace_polish_open) {
		polishresult->data = data->data;
		polishresult->datalength = data->length;
		free(seqbase);
		free(data);
	}
	return polishresult;
}

int32_t contig_get_length(Contig * contig, int32_t start, int32_t end)
{
	int32_t i = start, j = 0, length = 0;
	while (i < end || (i == end && j == 0)) {
		length++;
		contig_data_next(contig, &i, &j);
	}
	return length;
}

void contig_update_contig(Contig * contig, int32_t start, int32_t end, uint8_t * region, uint16_t index)
{
	int32_t i = start, j = 0;
	Base* p = contig->data[start];
	uint8_t* q = region;
	while (i < end || (i == end && j == index)) {
		p->base = *q;
		p = contig_data_next(contig, &i, &j);
		q++;
	}
}

void contig_clean_flag(Contig * contig, int32_t start, int32_t end, uint16_t flag)
{
	int32_t i = start, j = 0;
	Base* p = contig->data[start];
	while (i < end || (i == end && j == 0)) {
		p->flag &= flag;
		p = contig_data_next(contig, &i, &j);
	}
}

void contig_update_flag(Contig * contig, int32_t start, int32_t end, uint16_t flag)
{
	int32_t i = start, j = 0;
	Base* p = contig->data[start];
	while (i < end || (i == end && j == 0)) {
		p->flag |= flag;
		p = contig_data_next(contig, &i, &j);
	}
}

/*
static inline int possibly_expand_bam_data(bam1_t *b, size_t bytes) {
	uint32_t new_len = b->l_data + bytes;

	if (new_len > INT32_MAX || new_len < b->l_data) {
		errno = ENOMEM;
		return -1;
	}
	if (new_len <= b->m_data) return 0;
	return do_realloc_bam_data(b, new_len);
}

int do_realloc_bam_data(bam1_t *b, size_t desired)
{
	uint32_t new_m_data;
	uint8_t *new_data;
	new_m_data = desired;
	kroundup32(new_m_data);
	if (new_m_data < desired) {
		errno = ENOMEM; // Not strictly true but we can't store the size
		return -1;
	}
	new_data = realloc(b->data, new_m_data);
	if (!new_data) return -1;
	b->data = new_data;
	b->m_data = new_m_data;
	return 0;
}

static inline int realloc_bam_data(bam1_t *b, size_t desired)
{
	if (desired <= b->m_data) return 0;
	return do_realloc_bam_data(b, desired);
}

static inline void swap_data(const bam1_core_t *c, int l_data, uint8_t *data, int is_host)
{
	uint32_t *cigar = (uint32_t*)(data + c->l_qname);
	uint32_t i;
	for (i = 0; i < c->n_cigar; ++i) ed_swap_4p(&cigar[i]);
}

static inline int bam_tag2cigar(bam1_t *b, int recal_bin, int give_warning) // return 0 if CIGAR is untouched; 1 if CIGAR is updated with CG
{
	bam1_core_t *c = &b->core;
	uint32_t cigar_st, n_cigar4, CG_st, CG_en, ori_len = b->l_data, *cigar0, CG_len, fake_bytes;
	uint8_t *CG;

	// test where there is a real CIGAR in the CG tag to move
	if (c->n_cigar == 0 || c->tid < 0 || c->pos < 0) return 0;
	cigar0 = bam_get_cigar(b);
	if (bam_cigar_op(cigar0[0]) != BAM_CSOFT_CLIP || bam_cigar_oplen(cigar0[0]) != c->l_qseq) return 0;
	fake_bytes = c->n_cigar * 4;
	if ((CG = bam_aux_get(b, "CG")) == 0) return 0; // no CG tag
	if (CG[0] != 'B' || CG[1] != 'I') return 0; // not of type B,I
	CG_len = le_to_u32(CG + 2);
	if (CG_len < c->n_cigar || CG_len >= 1U << 29) return 0; // don't move if the real CIGAR length is shorter than the fake cigar length

	// move from the CG tag to the right position
	cigar_st = (uint8_t*)cigar0 - b->data;
	c->n_cigar = CG_len;
	n_cigar4 = c->n_cigar * 4;
	CG_st = CG - b->data - 2;
	CG_en = CG_st + 8 + n_cigar4;
	if (possibly_expand_bam_data(b, n_cigar4 - fake_bytes) < 0) return -1;
	b->l_data = b->l_data - fake_bytes + n_cigar4; // we need c->n_cigar-fake_bytes bytes to swap CIGAR to the right place
	memmove(b->data + cigar_st + n_cigar4, b->data + cigar_st + fake_bytes, ori_len - (cigar_st + fake_bytes)); // insert c->n_cigar-fake_bytes empty space to make room
	memcpy(b->data + cigar_st, b->data + (n_cigar4 - fake_bytes) + CG_st + 8, n_cigar4); // copy the real CIGAR to the right place; -fake_bytes for the fake CIGAR
	if (ori_len > CG_en) // move data after the CG tag
		memmove(b->data + CG_st + n_cigar4 - fake_bytes, b->data + CG_en + n_cigar4 - fake_bytes, ori_len - CG_en);
	b->l_data -= n_cigar4 + 8; // 8: CGBI (4 bytes) and CGBI length (4)
	if (recal_bin)
		b->core.bin = hts_reg2bin(b->core.pos, b->core.pos + bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b)), 14, 5);
	if (give_warning)
		hts_log_error("%s encodes a CIGAR with %d operators at the CG tag", bam_get_qname(b), c->n_cigar);
	return 1;
}



int bam_read2(BGZF *fp, bam1_t *b)
{
	bam1_core_t *c = &b->core;
	int32_t block_len, ret, i;
	uint32_t x[8], new_l_data;
	if ((ret = bgzf_read(fp, &block_len, 4)) != 4) {
		if (ret == 0) return -1; // normal end-of-file
		else return -2; // truncated
	}
	if (fp->is_be)
		ed_swap_4p(&block_len);
	if (block_len < 32) return -4;  // block_len includes core data
	if (bgzf_read(fp, x, 32) != 32) return -3;
	if (fp->is_be) {
		for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
	}
	c->tid = x[0]; c->pos = x[1];
	c->bin = x[2] >> 16; c->qual = x[2] >> 8 & 0xff; c->l_qname = x[2] & 0xff;
	c->l_extranul = (c->l_qname % 4 != 0) ? (4 - c->l_qname % 4) : 0;
	if ((uint32_t)c->l_qname + c->l_extranul > 255) // l_qname would overflow
		return -4;
	c->flag = x[3] >> 16; c->n_cigar = x[3] & 0xffff;
	c->l_qseq = x[4];
	c->mtid = x[5]; c->mpos = x[6]; c->isize = x[7];

	new_l_data = block_len - 32 + c->l_extranul;
	if (new_l_data > INT_MAX || c->l_qseq < 0 || c->l_qname < 1) return -4;
	if (((uint64_t)c->n_cigar << 2) + c->l_qname + c->l_extranul
		+ (((uint64_t)c->l_qseq + 1) >> 1) + c->l_qseq > (uint64_t) new_l_data)
		return -4;
	if (realloc_bam_data(b, new_l_data) < 0) return -4;
	b->l_data = new_l_data;

	if (bgzf_read(fp, b->data, c->l_qname) != c->l_qname) return -4;
	for (i = 0; i < c->l_extranul; ++i) b->data[c->l_qname + i] = '\0';
	c->l_qname += c->l_extranul;
	if (b->l_data < c->l_qname ||
		bgzf_read(fp, b->data + c->l_qname, b->l_data - c->l_qname) != b->l_data - c->l_qname)
		return -4;
	if (fp->is_be) swap_data(c, b->l_data, b->data, 0);
	if (bam_tag2cigar(b, 0, 0) < 0)
		return -4;

	if (c->n_cigar > 0) { // recompute "bin" and check CIGAR-qlen consistency
		int rlen, qlen;
		bam_cigar2rqlens(c->n_cigar, bam_get_cigar(b), &rlen, &qlen);
		if ((b->core.flag & BAM_FUNMAP)) rlen = 1;
		b->core.bin = hts_reg2bin(b->core.pos, b->core.pos + rlen, 14, 5);
		// Sanity check for broken CIGAR alignments
		if (c->l_qseq > 0 && !(c->flag & BAM_FUNMAP) && qlen != c->l_qseq) {
			hts_log_error("CIGAR and query sequence lengths differ for %s",
				bam_get_qname(b));
			return -4;
		}
	}

	return 4 + block_len;
}*/

hts_itr_t * contig_update_iter(Contig * contig, BGZF* fp, hts_idx_t* idx, hts_itr_t* iter, int32_t start, int32_t end, int32_t * iterend)
{
	if (iter != NULL) {
		if (end < *iterend) {
			iter->beg = start;
			iter->end = end + 1;
			iter->finished = 0;
		}
		else {
			hts_itr_destroy(iter);
			iter = NULL;
		}
	}
	if (iter == NULL) {
		iter = bam_itr_queryi(idx, contig->tid, start, end + 1);
		bam1_t* read = contig->read;
		if (iter->off && bgzf_seek(fp, iter->off[0].v, SEEK_SET) >= 0) {
			if (bam_read1(fp, read) >= 0 && read->core.tid == contig->tid) {
				*iterend = read->core.pos;
			}
			else {
				*iterend = contig->length;
			}
		}
	}
	return iter;
}

int contig_next_iter(Contig * contig, htsFile* fp, hts_idx_t* idx, hts_itr_t ** iter, int32_t start, int32_t end, int32_t * iterend, uint64_t * curr_off, bam1_t ** read, int32_t * curr_end, int32_t * nextposend, int32_t flag)
{
	if (flag >= 1) {
		*iter = contig_update_iter(contig, fp->fp.bgzf, idx, *iter, start, end, iterend);
		if ((*iter)->curr_off) {
			(*iter)->curr_off = *curr_off;
			(*iter)->curr_end = *curr_end;
			if (bgzf_seek(fp->fp.bgzf, (*iter)->curr_off, SEEK_SET) < 0) {
				(*iter)->curr_off = 0;
			}
			if ((*iter)->curr_off == 0) {
				(*iter)->i = -1;
			}
		}
		else {
			*curr_off = 0;
		}
		if (flag == 1) contig_swap_iter(*iter);
	}
	if ((*iter)->off) {
		flag = sam_itr_next(fp, *iter, *read);
	}
	else {
		flag = -1;
	}
	if (flag >= 0 && (*iter)->curr_end <= *nextposend) {
		*curr_off = (*iter)->curr_off;
		*curr_end = (*iter)->curr_end;
	}
	else {
		*nextposend = -1;
	}
	return flag;
}

void contig_write_to_file(const char* writefasta, const char* tigname, char * seq, int32_t length, int32_t tag)
{
	//char* filename = calloc(sizeof(char), strlen(writefasta) + 10);
	//sprintf(filename, "%s", writefasta);
	//FILE *fp = fopen(filename, "a");
	printf(">%s_%d\n%s\n", tigname, tag, seq);
	/*char* p = seq;
	int32_t len = 0;
	for (int i = 0; i < length; i += 100, p += 100) {
		if (i + 100 > length) {
			len = length - i;
		}
		else {
			len = 100;
		}
		fwrite(p, sizeof(char), len, fp);
		printf(fp, "\n");
	}
	free(filename);
	fclose(fp);*/
}

void contig_total(contig_run func, int32_t step, Configure* configure, int32_t tag)
{
	time_t start, end;

	start = time(0);

	faidx_t* faidx = fai_load(configure->fastafn);
	int len = faidx_nseq(faidx);
	const char* p = strrchr(configure->fastafn, '/');
	if (p == NULL) {
		p = configure->fastafn;
	}
	else {
		p++;
	}
	if (tag == -1) {
		//int a = 1;
		for (int i = 0; i < len; i++) {
			const char* name = faidx_iseq(faidx, i);
			//a = 0;
			/*if (a && strcmp(name, "ctg757_1_2") != 0) {
				continue;
			}
			a = 0;*/
			PolishResult* result = func(name, configure);
			contig_write_to_file(p, name, result->contig, strlen(result->contig), step);
			if (configure->trace_polish_open) {
				PolishPoint* p = result->data;
				for (int i = 0; i < result->datalength; i++, p++) {
					trace_log(TRACE_POLISH, "%s %d %d %c %c\n", name, p->pos, p->index, p->curbase, p->base);
				}
			}
			polishresult_destory(result);
		}
	}
	else {
		const char* name = faidx_iseq(faidx, tag);
		PolishResult* result = func(name, configure);
		contig_write_to_file(p, name, result->contig, strlen(result->contig), step);
		PolishPoint* p = result->data;
		for (int i = 0; i < result->length; i++, p++) {
			trace_log(TRACE_POLISH, "%s\t%d\t%d\t%c\t%c\n", name, p->pos, p->index, p->curbase, p->base);
		}
		polishresult_destory(result);
	}
	config_destory(configure);
	fai_destroy(faidx);

	end = time(0);
	trace_log(TRACE_TIME, "total time:%ds\n", end - start);
}

char * fasta_fetch(faidx_t * faidx, const char * tigname, int* length)
{
	int length1 = faidx_seq_len(faidx, tigname);
	//printf("name:%s,length:%d\n", name, length);
	char buffer[300];
	sprintf(buffer, "%s:%d-%d", tigname, 0, length1);
	//printf("%s\n", buffer);
	//char* buffer1 = "tig00000001:100-200";
	return fai_fetch(faidx, buffer, length);
}

inline void contig_swap_iter(hts_itr_t * iter)
{
	int32_t temp = iter->beg;
	iter->beg = iter->end;
	iter->end = temp;
}



/*

void Testbam::testbam()
{
	hash_map<int, int> a;
	clock_t start, end;
	printf("start\n", BAM_CMATCH);
	//bam_hdr_t* temp = bam_hdr_init();
	//bam_hdr_destroy(temp);
	char* fn = "/export/personal1/daijl/data/NX032/bm.sort.bam";
	BGZF* fp=bgzf_open(fn, "r");
	bam_hdr_t* bamhdr=bam_hdr_read(fp);
	bam1_t* read=bam_init1();
	int result=bam_read1(fp, read);
	if (result < 0) {
		printf("error:%d\n", result);
	}
	char* name = bam_get_qname(read);
	printf("%s,%d\n", name,sizeof(*name));
	kstring_t str;
	str.l = str.m = 0; str.s = NULL;
	if (sam_format1(bamhdr, read, &str) < 0) {
		free(str.s);
		str.s = NULL;
		printf("error samformat\n");
	}
	printf("%s\n", str.s);
	printf("%d\n", bamhdr->target_len[read->core.tid]);
	uint32_t length = bamhdr->target_len[read->core.tid];
	start = clock();
	char* fnidx = "/export/personal1/daijl/data/NX032/bm.sort.bam.bai";
	hts_idx_t* idx=bam_index_load(fnidx);
	hts_itr_t* iter=bam_itr_queryi(idx, read->core.tid, 0, length+1);
	int count = 0;
	while (hts_itr_next(fp, iter, read, 0) !=-1)
	{
		str.l = str.m = 0; str.s = NULL;
		if (sam_format1(bamhdr, read, &str) < 0) {
			free(str.s);
			str.s = NULL;
			printf("error samformat\n");
		}
		//printf("%s\n", str.s);
		free(str.s);
		count += 1;
	}
	end = clock();
	double endtime = (double)(end - start)*1000 / CLOCKS_PER_SEC;
	printf("totalcount:%d\n,needtime:%fms", count,endtime);
	hts_itr_destroy(iter);
	bam_destroy1(read);
	bam_hdr_destroy(bamhdr);
	bgzf_close(fp);
	printf("finish\n");
}
*/
