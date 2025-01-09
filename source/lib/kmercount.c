#include "kmercount.h"

inline KmerScore * ks_init(int32_t length)
{
	KmerScore* ks = calloc(sizeof(KmerScore), 1);
	if (length) {
		ks->region = calloc(sizeof(uint8_t), length);
	}
	else {
		ks->region = NULL;
	}
	return ks;
}

inline void ks_destory(KmerScore * ks)
{
	if (ks->region != NULL) {
		free(ks->region);
	}
	free(ks);
}

inline void ks_clean(KmerScore * ks, int32_t length)
{
	if (ks->region == NULL) {
		ks->region = calloc(sizeof(uint8_t), length);
	}
	ks->length = 0;
	ks->qual = 0;
	ks->mapqual = 0;
	ks->num = 0;
}

inline void ks_append(KmerScore * ks, uint8_t base)
{
	ks->region[ks->length++] = base;
}

inline int32_t ks_compare_region(void * firstks, void * secondks)
{
	KmerScore *first = firstks, *second = secondks;
	return memcmp(first->region, second->region, first->length);
}

int32_t ks_compare_length(void * firstks, void * secondks)
{
	KmerScore *first = firstks, *second = secondks;
	if (first->length == second->length) {
		return 0;
	}
	return 1;
}

int32_t ks_compare_mapqual(void * firstks, void * secondks)
{
	KmerScore *first = firstks, *second = secondks;
	if (first->mapqual == second->mapqual) {
		return 0;
	}
	return 1;
}

inline int32_t ks_compare(void * firstks, void * secondks)
{
	if (firstks == secondks) {
		return 0;
	}
	KmerScore *first = firstks, *second = secondks;
	if (first->num > second->num) {
		return 1;
	}
	else if (first->num < second->num) {
		return -1;
	}
	if (first->mapqual > second->mapqual) {
		return 1;
	}
	else if (first->mapqual < second->mapqual) {
		return -1;
	}
	if (first->qual > second->qual) {
		return 1;
	}
	else if (first->qual < second->qual) {
		return -1;
	}
	return 0;
}




PolishResult * kmer_count(const char * tigname, Configure* configure)
{
	Contig* contig = contig_init(tigname, configure, sizeof(Kmer));

	contig->read_fliter = contig_read_fliter;
	SeqList* nodepth = contig_get_region(contig, 0, contig->length - 1, 0, contig->configure->min_len_ldr, 0x1, contig_brim_no_extension);
	SeqList* kmerregion = contig_get_region(contig, 0, contig->length - 1, contig->configure->min_len_inter_kmer, 0, 0x1, contig_brim_with_extension);
	if (kmerregion->length > 0) {
		contig_merge_region(contig, kmerregion);
		contig_create_insert_region(contig, kmerregion, 0);
	}
	if (nodepth->length > 0) {
		contig_merge_region(contig, nodepth);
		contig_create_insert_region(contig, nodepth, 0);
		int32_t* p = nodepth->data;
		for (int i = 0; i < nodepth->length; i += 2, p += 2) {
			contig_score_correct(contig, *p, *(p + 1), 0x12, contig->configure->indel_balance_factor_sgs);
		}
	}
	seqlist_destory(nodepth);

	if (kmerregion->length > 0) {
		nodepth = ss_spilt_region(contig, kmerregion, 0x1, contig->configure->max_len_kmer);
		ss_kmer_correct(contig, nodepth, NULL, 0);
		seqlist_destory(nodepth);
	}
	seqlist_destory(kmerregion);

	PolishResult* result = contig_get_contig(contig, 0, contig->length - 1, FLAG_ZERO);

	contig_destory(contig);

	return result;
}

SeqList * ss_spilt_region(Contig * contig, SeqList * seqlist, uint8_t flag, uint8_t max)
{
	SeqList* result = seqlist_init(seqlist->stepsize, seqlist->length * 2);
	SeqList* temp = seqlist_init(seqlist->stepsize, seqlist->length);
	int32_t	*p = seqlist->data, j, k, qstart, qend, *ptemp;
	Base* q = NULL;
	for (int i = 0; i < seqlist->length; i += 2, p += 2) {
		seqlist_append(result, p);
		if (*(p + 1) - *p > max) {
			j = *p;
			k = 0;
			qstart = qend = -1;
			q = contig->data[j];
			temp->length = 0;
			while (j < *(p + 1) || (j == *(p + 1) && k == 0)) {
				if ((q->flag&flag) != 0) {
					break;
				}
				q = contig_data_next(contig, &j, &k);
			}
			while (j < *(p + 1) || (j == *(p + 1) && k == 0)) {
				if ((q->flag&flag) == 0) {
					if (qstart == -1) {
						qstart = j;
					}
					qend = j;
				}
				else if (qstart != -1) {
					seqlist_append(temp, &qstart);
					seqlist_append(temp, &qend);
					qstart = qend = -1;
				}
				q = contig_data_next(contig, &j, &k);
			}
			ptemp = temp->data;
			for (j = 0; j < temp->length; j += 2, ptemp += 2) {
				k = (*ptemp + *(ptemp + 1)) >> 1;
				seqlist_append(result, &k);
				seqlist_append(result, &k);
			}
		}
		seqlist_append(result, p + 1);
	}
	seqlist_destory(temp);
	return result;
}

void ss_kmer_correct(Contig * contig, SeqList * region, SeqList * nodepth, int32_t flagzero)
{
	int32_t *p = region->data, count = 0;
	SeqList* regiondata = seqlist_init(sizeof(KmerScore), 10);
	KmerScore* ks = NULL;

	int32_t i = 0, j, end1 = 0, end2 = 0, nextposend, nextposend1, curr_end1 = 0, curr_end2 = 0, flag, length;
	int32_t start, end;
	hts_itr_t *iter1 = NULL, *iter2 = NULL;
	uint64_t curroff1 = 0, curroff2 = 0;
	bam1_t* read = bam_init1(), *tdread = bam_init1();

	contig->read_fliter = contig_read_fliter;
	for (i = 0; i < region->length; i += 2, p += 2) {
		start = *p;
		end = *(p + 1);
		length = contig_get_length(contig, start, end);
		count = 0;
		nextposend = i + 2 < region->length ? *(p + 3) : -1;
		flag = 1;
		nextposend1 = nextposend;
		while (contig_next_iter(contig, contig->fp, contig->idx, &iter1, start, end, &end1, &curroff1, &read, &curr_end1, &nextposend1, flag) >= 0) {
			if (contig->read_fliter(contig, read) == 2) {
				ks = ss_kmer_get_region(contig, read, start, end, length, regiondata, ks, 0, -1, -1, 0, flagzero);
				if (ks->mapqual == MAX_MAPQ) {
					count++;
					if (count >= contig->configure->max_count_kmer) {
						break;
					}
				}
				ks_clean(ks, length);
			}
			flag = 0;
		}
		if (regiondata->length == 0) {
			flag = 1;
			nextposend1 = nextposend;
			while (contig_next_iter(contig, contig->fp, contig->idx, &iter2, start, end, &end2, &curroff2, &tdread, &curr_end2, &nextposend1, flag) >= 0) {
				if (contig->read_fliter(contig, read) == 1) {
					ks = ss_kmer_get_region(contig, read, start, end, length, regiondata, ks, 0, -1, -1, 0, flagzero);
					ks_clean(ks, length);
				}
				flag = 0;
			}
		}

		if (regiondata->length > 0) {
			if (flagzero) {
				contig_clean_flag(contig, start, end, FLAG_ZERO_N);
			}
			KmerScore* bestregion = NULL;
			if (count == contig->configure->max_count_kmer) {
				ks->mapqual = MAX_MAPQ * count;
				bestregion = seqlist_find(regiondata, ks, ks_compare_mapqual);
			}
			ks_destory(ks);
			if (bestregion == NULL) {
				bestregion = ks = regiondata->data;
				for (j = 0; j < regiondata->length; j++, ks++) {
					if (ks_compare(bestregion, ks) < 0) {
						bestregion = ks;
					}
				}
			}
			contig_update_contig(contig, start, end, bestregion->region, 0);
			ks = NULL;
		}
		else if (nodepth) {
			seqlist_append(nodepth, &start);
			seqlist_append(nodepth, &end);
		}
		if (ks) ks_destory(ks);
		ks = regiondata->data;
		for (j = 0; j < regiondata->length; j++, ks++) {
			free(ks->region);
		}
		regiondata->length = 0;
		ks = NULL;
	}

	bam_destroy1(read);
	bam_destroy1(tdread);
	if (iter1) hts_itr_destroy(iter1);
	if (iter2) hts_itr_destroy(iter2);

	seqlist_destory(regiondata);
}

//
//SeqList* regiondata = seqlist_init(sizeof(KmerScore), 10);
//BGZF* fp = bgzf_open(contig->bamfn, "r");
//hts_idx_t* idx = bam_load_idx(contig->bamfn);
/*
char* temp = bam_get_qname(read);
if (strcmp(temp, "ST-E00600:249:H7NKLALXX:1:1166:32127:32706") == 0 || strcmp(temp,"ST-E00600:249:H7NKLALXX:1:1619:1877:35399")==0) {
	int a = 1;
}
*/
/*
void ss_kmer_correct(Contig * contig, int32_t start, int32_t end, SeqList* regiondata, hts_itr_t* itert) {
	int32_t length = contig_get_length(contig, start, end);
	KmerScore* ks = NULL;

	hts_itr_t* iter = itert ? itert : bam_itr_queryi(contig->idx, contig->tid, start, end);
	ss_swap_iter(iter);
	bam1_t* read = contig->read;
	int32_t count = 0;
	while (hts_itr_next(contig->fp, iter, read, 0) >= 0)
	{
		if (contig->read_fliter(contig, read) == 2) {
			ks = ss_kmer_get_region(contig, read, start, end, length, regiondata, ks, 0, -1, -1);
			if (ks->mapqual == 60) {
				count++;
				if (count >= 50) {
					break;
				}
			}
		}
	}
	if (regiondata->length == 0) {
		hts_itr_destroy(iter);
		iter = bam_itr_queryi(contig->idx, contig->tid, start, end);
		ss_swap_iter(iter);
		while (hts_itr_next(contig->fp, iter, read, 0) >= 0)
		{
			if (contig->read_fliter(contig, read) == 1) {
				ks = ss_kmer_get_region(contig, read, start, end, length, regiondata, ks, 0, -1, -1);
			}
		}
	}
	if (ks != NULL) {
		ks_destory(ks);
	}
	int32_t i = 0;
	if (regiondata->length > 0) {
		KmerScore* bestregion = ks = regiondata->data;
		for (; i < regiondata->length; i++, ks++) {
			if (ks_compare(bestregion, ks) < 0) {
				bestregion = ks;
			}
		}
		contig_update_contig(contig, start, end, bestregion->region, 0);
	}
	ks = regiondata->data;
	for (i = 0; i < regiondata->length; i++, ks++) {
		free(ks->region);
	}
	regiondata->length = 0;
	if (itert == NULL) {
		hts_itr_destroy(iter);
	}
}
*/
//bam_destroy1(read);
//hts_idx_destroy(idx);
//bgzf_close(fp);

inline KmerScore * ss_kmer_get_region(Contig * contig, bam1_t* read, int32_t start, int32_t end, int32_t length, SeqList * regiondata, KmerScore * ks, int32_t* count, int32_t left, int32_t right, int32_t flag, int32_t flagzero)
{
	if (ks == NULL) {
		ks = ks_init(length);
	}
	int32_t check = 0;
	if (left != -1) {
		check = 2;
	}
	int32_t result = ss_parse_read_kmer(contig, read, start, end, ks, left, right, flagzero);
	if (ks->length == length && result >= check) {
		ks->length += flag;
		KmerScore* p = seqlist_find(regiondata, ks, ks_compare_region);
		if (p == NULL) {
			ks->num = 1;
			seqlist_append(regiondata, ks);
			ks->region = NULL;
		}
		else {
			p->num++;
			p->mapqual += ks->mapqual;
			p->qual += ks->qual;
		}
		if (count) {
			(*count)++;
		}
	}
	else {
		ks->mapqual = 0;
	}
	return ks;
}

inline int32_t ss_parse_read_kmer(Contig * contig, bam1_t* read, int32_t start, int32_t end, KmerScore* ks, int32_t left, int32_t right, int32_t flagzero)
{
	int32_t result = 0;
	if (read->core.n_cigar) {
		/*char* name = bam_get_qname(read);
		char* name1 = "ST-E00600:249:H7NKLALXX:1:1374:8983:11757";
		if (strcmp(name, name1) == 0) {
			int a = 1;
		}*/
		int32_t pos = read->core.pos, qpos = 0, qstart, qend, i, j, k, len, del = 0;
		uint8_t curcigar, lastcigar = BAM_CINS;
		uint32_t *cigar = bam_get_cigar(read);
		uint8_t* seq = bam_get_seq(read);
		uint8_t* qual = bam_get_qual(read);
		SeqList* p;
		contig_cut_read(contig, read, &qstart, &qend);
		ks->mapqual = read->core.qual;
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
									ks_append(ks, BASE_DEL);
									if (flagzero == 0) ((Base*)seqlist_index(p, k))->flag &= FLAG_ZERO_N;
									del++;
								}
							}
						}
						if (curcigar == BAM_CDEL) {
							ks_append(ks, BASE_DEL);
							//del++;
						}
						else {
							ks_append(ks, bam_seqi(seq, qpos));
							ks->qual += qual[qpos];
						}
						if (flagzero == 0) contig->data[pos]->flag &= FLAG_ZERO_N;
					}
					if (left == pos || right == pos) {
						if (bam_seqi(seq, qpos) == contig->data[pos]->base) {
							result++;
						}
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
							ks_append(ks, bam_seqi(seq, qpos));
							ks->qual += qual[qpos];
							if (flagzero == 0) ((Base*)seqlist_index(p, j))->flag &= FLAG_ZERO_N;
						}
					}
					if (pos > start && pos <= end && qpos > qstart && qpos <= qend + 1) {
						for (; j < p->length; j++) {
							ks_append(ks, BASE_DEL);
							if (flagzero == 0) ((Base*)seqlist_index(p, j))->flag &= FLAG_ZERO_N;
							del++;
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
		if (ks->length > 0 && ks->length != del) {
			ks->qual /= ks->length - del;
		}
		else {
			ks->qual = 0;
		}
	}
	return result;
}


