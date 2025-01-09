#include "lgspolish.h"



void tdkmer_add_item(TdKmer *tdkmer, uint8_t base, int32_t pos, int16_t index)
{
	memmove(&tdkmer->kmer, &tdkmer->kmer[1], sizeof(TdKmerItem)*(KMER_K - 1));
	tdkmer->kmer[KMER_K - 1].base = base;
	tdkmer->kmer[KMER_K - 1].pos = pos;
	tdkmer->kmer[KMER_K - 1].index = index;
}

PolishResult * lgspolish(const char * tigname, Configure * configure)
{
	Contig* contig = contig_init(tigname, configure, sizeof(TdKmer));

	contig->read_fliter = contig_read_fliter2;
	td_region_score_chain(contig, 0, contig->length - 1, contig->configure->indel_balance_factor_lgs);

	PolishResult* result = contig_get_contig(contig, 0, contig->length - 1, 0);

	contig_destory(contig);

	return result;
}

void td_region_score_chain(Contig * contig, int32_t start, int32_t end, double rate)
{
	td_contig_as_read(contig, start, end);

	td_parse_region(contig, contig->tfp, start, end, 1);
	
	td_region_score(contig, start, end, rate);
	td_region_correct(contig, start, end);
}

int32_t comparetdkmer(void * first, void * second)
{
	return memcmp(first, second, sizeof(TdKmerItem)*KMER_K);
}

void td_base_add_data(Base * base, TdKmer * kmer)
{
	TdKmer* result = seqlist_find(base->data, kmer, comparetdkmer);
	if (result == NULL) {
		seqlist_append(base->data, kmer);
	}
	else {
		result->count++;
	}
	base->count++;
}

void td_contig_as_read(Contig * contig, int32_t start, int32_t end)
{
	uint16_t kmer = 0x0;
	TdKmer tdkmer;
	memset(&tdkmer, 0, sizeof(TdKmer));
	tdkmer.count = 1;
	int32_t i = start, j = 0;
	Base* p = contig->data[start];
	while (i < end || (i == end && j == 0)) {
		p->refkmer = kmer = contig_left_kmer(kmer, p->base, BASE_SHIFT);
		tdkmer_add_item(&tdkmer, p->base, i, j);
		td_base_add_data(p, &tdkmer);
		p = contig_data_next(contig, &i, &j);
	}
}

void td_parse_read(Contig * contig, bam1_t * read, int32_t start, int32_t end)
{
	if (read->core.n_cigar) {
		TdKmer tdkmer;
		memset(&tdkmer, 0, sizeof(TdKmer));
		tdkmer.count = 1;
		int32_t pos = read->core.pos, qpos = 0, qstart, qend, i, j, len;
		uint8_t curcigar, flag;
		uint32_t *cigar = bam_get_cigar(read);
		uint8_t* seq = bam_get_seq(read);
		SeqList* p;
		Base* base;
		contig_cut_read(contig, read, &qstart, &qend);
		for (i = 0; i < read->core.n_cigar; ++i) {
			len = bam_cigar_oplen(cigar[i]);
			curcigar = bam_cigar_op(cigar[i]);
			switch (curcigar) {
			case BAM_CMATCH:
			case BAM_CDEL:
				for (j = 0; j < len; j++, pos++) {
					if (pos >= start && pos <= end && qpos >= qstart && qpos <= qend) {
						if (curcigar == BAM_CDEL) {
							tdkmer_add_item(&tdkmer, BASE_DEL, pos, 0);
						}
						else {
							tdkmer_add_item(&tdkmer, bam_seqi(seq, qpos), pos, 0);
						}
						td_base_add_data(contig->data[pos], &tdkmer);
					}
					if (curcigar != BAM_CDEL) {
						qpos++;
					}
				}
				break;
			case BAM_CINS:
				if (contig->data[pos - 1]->insert == NULL) {
					contig->data[pos - 1]->insert = seqlist_init(sizeof(Base), len);
				}
				p = contig->data[pos - 1]->insert;
				if (p->maxsize < len) {
					p->data = realloc(p->data, len*p->stepsize);
					if (p->data == NULL) {
						fprintf(stderr, "seqlist malloc failed,the memory is full!\n");
						exit(1);
					}
					p->maxsize = len;
				}
				flag = contig->data[pos - 1]->flag;
				for (j = p->length; j < len; j++, contig->inslength++) {
					base = base_init(sizeof(TdKmer));
					base->flag = flag;
					seqlist_append(p, base);
					free(base);
				};
				for (j = 0; j < len; j++, qpos++) {
					if (pos > start && pos <= end && qpos >= qstart && qpos <= qend) {
						tdkmer_add_item(&tdkmer, bam_seqi(seq, qpos), pos - 1, j + 1);
						td_base_add_data(seqlist_index(p, j), &tdkmer);
					}
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

void td_parse_region(Contig * contig, htsFile * fp, int32_t start, int32_t end, uint8_t filterlevel)
{
	hts_itr_t* iter = bam_itr_queryi(contig->tidx, contig->tid, start, end + 1);
	bam1_t* read = contig->read;
	while (sam_itr_next(fp, iter, read) >= 0)
	{
		if (contig->read_fliter(contig, read) == filterlevel) {
			td_parse_read(contig, read, start, end);
		}
	}
	hts_itr_destroy(iter);
}

void td_base_add_score(Base * base, TdKmer * tdkmer, double score, uint16_t index)
{
	Score temp = { .base = tdkmer->kmer[KMER_K-1].base,.kmer = index,.score = score };
	Score* result = seqlist_find(base->score, &temp, comparescore);
	if (result == NULL) {
		seqlist_append(base->score, &temp);
	}
	else {
		seqlist_update(base->score, result, &temp);
	}
}

void td_region_score(Contig * contig, int32_t start, int32_t end, double rate)
{
	int32_t i = start, j = 0, k, n;
	Base *p = NULL, *q = contig->data[start];
	TdKmer* ptdkmer;
	TdKmerItem* lastitem, *curitem;
	double score = 0;
	Score* pscore = NULL;
	uint16_t count, total, kmer;
	while (i < end || (i == end && j == 0)) {
		base_clean_score(q);
		ptdkmer = q->data->data;
		total = q->count;
		if (j != 0) {
			p = contig->data[i];
			if (p->insert->length <= 4 || total / (double)p->count < 0.2) {
				total = contig->data[i]->count;
			}
			else {
				total = 1;
				//trace_log(TRACE_ERROR, "%s %d\n", contig->tigname, i);
			}
		}
		if (total > 1) {
			total--;
		}
		for (k = 0; k < q->data->length; k++, ptdkmer++) {
			//score = base_get_score(lastbase, p->kmer >> BASE_SHIFT)->score;
			lastitem = &ptdkmer->kmer[KMER_K - 2];
			curitem = &ptdkmer->kmer[KMER_K - 1];
			if (lastitem->base == 0) {
				if (i > 0) {
					score = base_max_score(contig->data[i - 1])->score;
				}
				else {
					score = 0;
				}
			}
			else {
				score = base_get_score(contig_index(contig,lastitem->pos,lastitem->index), lastitem->base)->score;
			}
			count = ptdkmer->count;

			kmer = 0;
			for (n = 0; n < KMER_K; n++) {
				kmer = contig_left_kmer(kmer, ptdkmer->kmer[n].base, BASE_SHIFT);
			}

			if (kmer == q->refkmer && q->count > 1) {
				count--;
			}
			score += count - total * rate;
			pscore = base_get_score(q, curitem->base);
			if (pscore == NULL || pscore->score < score) {
				td_base_add_score(q, ptdkmer, score, k);
			}
		}
		q = contig_data_next(contig, &i, &j);
	}
}

void td_region_correct(Contig * contig, int32_t start, int32_t end)
{
	Base* base = contig->data[end];
	Score* score = base_max_score(base);
	TdKmer* ptdkmer;
	int32_t i = end, j = 0;
	while (i > start || (i == start && j == 0))
	{
		base->base = score->base;
		ptdkmer = seqlist_index(base->data, score->kmer);
		if (ptdkmer->kmer[KMER_K - 2].base) {
			i = ptdkmer->kmer[KMER_K - 2].pos;
			j = ptdkmer->kmer[KMER_K - 2].index;
			base = contig_index(contig, i, j);
		}
		else{
			i--;
			j = 0;
			if (i >= 0) {
				base = contig->data[i];
			}
		}
		score = base_get_score(base, ptdkmer->kmer[KMER_K - 2].base);
	}
}
