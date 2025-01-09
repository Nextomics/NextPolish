#include "snpphase.h"

inline Snps * snps_init(int32_t right, int32_t length)
{
	Snps* snps = calloc(sizeof(Snps), 1);
	for (int i = 0; i < SNP_NUM; i++) {
		snps->region[i] = calloc(sizeof(uint8_t), length);
	}
	snps->link = seqlist_init(sizeof(KmerScore), SNP_NUM);
	snps->right = right;
	snps->length = length;
	return snps;
}

inline void snps_destory(Snps * snps)
{
	if (snps != NULL) {
		for (int i = 0; i < SNP_NUM; i++) {
			free(snps->region[i]);
		}
		seqlist_destory(snps->link);
		free(snps);
	}
}

inline void snps_resize(Snps * snps, int32_t length)
{
	if (length != snps->length) {
		for (int i = 0; i < SNP_NUM; i++) {
			snps->region[i] = realloc(snps->region[i], sizeof(uint8_t)*length);
		}
	}
}

inline void snps_add_region(Snps * snps, uint8_t * region, int32_t i)
{
	memcpy(snps->region[i], region, snps->length);
}

int32_t snps_get_index(Snps * snps, uint8_t * region)
{
	for (int i = 0; i < SNP_NUM; i++) {
		if (memcmp(snps->region[i], region, snps->length) == 0) {
			return i;
		}
	}
	return -1;
}

inline SnpsList * snpslist_init(SeqList * snps)
{
	SnpsList* result = malloc(sizeof(SnpsList));
	result->data = calloc(sizeof(Snps*), snps->length);
	result->length = snps->length;
	Snps** p = snps->data;
	for (int i = 0; i < snps->length; i++, p++) {
		result->data[i] = *p;
	}
	return result;
}

inline void snpslist_destory(SnpsList * snpslist)
{
	if (snpslist->data) free(snpslist->data);
	free(snpslist);
}

int32_t snpslist_find(SnpsList * snpslist, int32_t pos)
{
	int32_t i = 0, j = snpslist->length - 1, mid, qpos;
	while (i <= j) {
		mid = (i + j) / 2;
		qpos = snpslist->data[mid]->pos;
		if (qpos == pos) {
			return mid;
		}
		if (qpos < pos) {
			i = mid + 1;
		}
		else {
			j = mid - 1;
		}
	}
	return -1;
}

PolishResult* snp_phase(const char * tigname, Configure* configure)
{
	Contig* contig = contig_init(tigname, configure, sizeof(Kmer));
	contig->read_fliter = contig_read_fliter;

	contig_create_insert(contig, contig->fp, contig->idx, 0, contig->length - 1, 0);

	contig_parse_region(contig, contig->fp, 0, contig->length - 1, 2, 16);

	SeqList* snps = ts_find_snps(contig, 0, contig->length - 1);
	
	SeqList* nodepth = contig_get_region(contig, 0, contig->length - 1, contig->configure->ext_len_edge, 0, FLAG_DEPTH, contig_brim_no_extension);
	if (nodepth->length > 0) {
		contig_merge_region(contig, nodepth);
		int32_t	*p = nodepth->data;
		for (int i = 0; i < nodepth->length; i += 2, p += 2) {
			contig_update_flag(contig, *p, *(p + 1), FLAG_INSERT);
		}
	}
	contig->read_fliter = contig_read_fliter2;
	contig_create_insert(contig, contig->tfp, contig->tidx, 0, contig->length - 1, FLAG_INSERT | FLAG_SNP);

	ts_fliter_snps(contig, snps);
	
	if (nodepth->length > 0) {
		//contig_merge_region(contig, nodepth);
		//trace_log(TRACE_STEP, "total %d region need lower depth correct\n", nodepth->length / 2);
		ts_correct_lower_depth(contig, nodepth);
	}
	seqlist_destory(nodepth);
	
	if (snps->length > 1) {
		ts_find_snps_link(contig, snps);
		ts_snps_score(contig, snps);
		ts_snps_correct(contig, snps);
	}
	Snps** p = snps->data;
	for (int i = 0; i < snps->length; i++, p++) {
		snps_destory(*p);
	}
	seqlist_destory(snps);
	
	PolishResult* result = contig_get_contig(contig, 0, contig->length - 1, FLAG_THIRD);

	contig_destory(contig);

	return result;
}

inline SeqList* ts_find_snps(Contig * contig, int32_t start, int32_t end)
{
	SeqList* result = seqlist_init(sizeof(Snps*), contig->configure->region_count);
	int32_t i = start, j = 0, k, count, flag = 0, flag1 = 0, lasti = start, lastj = 0;
	uint8_t base;
	Base* p = contig->data[start];
	Kmer* maxn[SNP_NUM];
	double rate;
	Snps* snps = snps_init(contig->length - 1, 1), **q = NULL;
	while (i < end || (i == end && j == 0)) {
		if (p->count == 0) {
			p->flag |= FLAG_ZERO;
		}
		else {
			p->flag &= FLAG_ZERO_N;
		}
		if (p->count <= contig->configure->min_depth_snp) {
			p->flag |= FLAG_DEPTH;
		}
		else {
			p->flag &= FLAG_DEPTH_N;
		}
		flag = 0;
		if (p->count > 0) {
			count = base_get_nlargest(p, maxn, SNP_NUM);
			if (count == 1) {
				rate = 0;
			}
			else {
				rate = maxn[1]->count / (double)maxn[0]->count;
			}
			flag = ts_check_snps(contig, p->count, rate, maxn[0]->kmer == p->base);
			if (flag == 2) {
				p->base = maxn[0]->kmer;
			}
			else if (flag == 1) {
				if (j == 0 || (contig->data[i]->flag&FLAG_SNP) == 0) {
					contig->data[i]->flag |= FLAG_SNP;
					snps->left = lasti;
					snps->pos = i;
					flag1 = 1;
					for (k = 0; k < count; k++) {
						base = maxn[k]->kmer;
						snps_add_region(snps, &base, k);
					}
					if (count < SNP_NUM) {
						snps_add_region(snps, &p->base, count);
					}
					seqlist_append(result, &snps);
					snps = snps_init(contig->length - 1, 1);
				}
			}
		}
		if (flag != 1 && ((contig->data[i]->flag&FLAG_SNP) == 0) && (contig->data[i]->insert == NULL || j == contig->data[i]->insert->length)) {
			lasti = i;
			if (flag1) {
				q = seqlist_index(result, lastj);
				for (; lastj < result->length; lastj++, q++) {
					(*q)->right = lasti;
				}
				flag1 = 0;
			}
		}
		p = contig_data_next(contig, &i, &j);
	}
	snps_destory(snps);
	return result;
}

inline int32_t ts_check_snps(Contig * contig, int32_t count, double rate, int32_t flag)
{
	if (rate < contig->configure->min_snp_factor_sgs && flag) {
		return 0;
	}
	if (rate == 0 || (count >= contig->configure->min_count_snp && flag == 0 && rate < contig->configure->min_snp_factor_sgs)) {
		return 2;
	}
	return 1;
}

inline void ts_fliter_snps(Contig * contig, SeqList * snps)
{
	if (snps->length > 0) {
		Snps** p = snps->data, *q = NULL, **psnps = p;
		Base* base = NULL;
		int32_t i = 0, j, end1 = 0, end2 = 0, nextend, curr_end1 = 0, curr_end2 = 0, flag, flag1, length, total = 0, count = 0;
		int32_t start, end;
		hts_itr_t *iter1 = NULL, *iter2 = NULL;
		uint64_t curroff1 = 0, curroff2 = 0;
		SeqList* regiondata = seqlist_init(sizeof(KmerScore), 10);
		KmerScore* ks = NULL, *tempks;
		bam1_t* read = bam_init1(), *tdread = bam_init1();
		KmerScore* maxn[SNP_NUM];
		double rate;
		for (; i < snps->length; i++, p++) {
			q = *p;
			base = contig->data[q->pos];
			nextend = i + 1 < snps->length ? (*(p + 1))->pos : -1;
			length = 1;
			start = end = q->pos;
			total = 0;
			flag = 0;
			if (base->insert && base->insert->length > 0) {
				end = q->pos + 1;
				length += base->insert->length + 1;
				contig->read_fliter = contig_read_fliter;
				flag = 1;
				while (contig_next_iter(contig, contig->fp, contig->idx, &iter1, start, end, &end1, &curroff1, &read, &curr_end1, &nextend, flag) >= 0) {
					if (contig->read_fliter(contig, read) == 2) {
						ks = ss_kmer_get_region(contig, read, start, end, length, regiondata, ks, &total, -1, -1, -1, 0);
						ks_clean(ks, length);
					}
					flag = 0;
				}
				flag = 1;
			}
			else {
				ks = ks_init(length);
				total = base->count;
			}
			if (total <= contig->configure->min_count_snp) {
				if (length == 1) {
					Kmer* temp = base->data->data;
					for (j = 0; j < base->data->length; j++, temp++) {
						ks->region[0] = temp->kmer;
						ks->num = temp->count;
						ks->mapqual = READ_MAPQ * ks->num;
						ks->qual = BASE_QUAL * ks->num;
						seqlist_append(regiondata, ks);
						ks->region = calloc(sizeof(uint8_t), length);
					}
				}
				memset(ks->region, BASE_DEL, length);
				flag1 = seqlist_get_index(regiondata, ks, ks_compare_region);
				nextend = i + 1 < snps->length ? (*(p + 1))->right : -1;
				contig->read_fliter = contig_read_fliter2;
				flag = 1;
				while (contig_next_iter(contig, contig->tfp, contig->tidx, &iter2, start, end, &end2, &curroff2, &tdread, &curr_end2, &nextend, flag) >= 0) {
					if (contig->read_fliter(contig, tdread) == 1) {
						ks = ss_kmer_get_region(contig, tdread, start, end, length, regiondata, ks, &total, q->left, q->right, -1, 1);
						ks_clean(ks, length);
					}
					flag = 0;
				}
				flag = 1;
				if (flag1 == -1) {
					memset(ks->region, BASE_DEL, length);
					flag1 = seqlist_get_index(regiondata, ks, ks_compare_region);
					if (flag1 != -1) {
						tempks = seqlist_index(regiondata, flag1);
						free(tempks->region);
						seqlist_remove(regiondata, flag1);
					}
				}
			}
			if (flag) {
				flag1 = ts_get_nlargest(regiondata, maxn, SNP_NUM);
				if (flag1 == 1) {
					rate = 0;
				}
				else {
					rate = maxn[1]->num / (double)maxn[0]->num;
				}
				memset(ks->region, BASE_DEL, length);
				ks->region[0] = base->base;
				flag = ts_check_snps(contig, total, rate, ks_compare_region(maxn[0], ks) == 0);
				if (flag == 1) {
					snps_resize(q, length);
					q->length = maxn[0]->length;
					for (j = 0; j < flag1; j++) {
						snps_add_region(q, maxn[j]->region, j);
					}
					if (flag1 < SNP_NUM) {
						snps_add_region(q, ks->region, flag1);
					}
					*psnps = *p;
					psnps++;
					count++;

				}
				else {
					if (flag == 2) {
						base->base = maxn[0]->region[0];
						contig_update_contig(contig, start, end, maxn[0]->region, -1);
					}
					contig->data[start]->flag &= 0xf7;
					snps_destory(*p);
					*p = NULL;
				}
			}
			else {
				*psnps = *p;
				psnps++;
				count++;
			}

			if (ks) {
				ks_destory(ks);
			}
			ks = regiondata->data;
			for (j = 0; j < regiondata->length; j++, ks++) {
				free(ks->region);
			}
			regiondata->length = 0;
			ks = NULL;
		}
		snps->length = count;
		bam_destroy1(read);
		bam_destroy1(tdread);
		seqlist_destory(regiondata);
		if (iter1) hts_itr_destroy(iter1);
		if (iter2) hts_itr_destroy(iter2);
	}
}

void ts_find_snps_link(Contig * contig, SeqList * snps)
{
	if (snps->length > 1) {
		SeqList* snpsregion = ts_find_snp_region(contig, snps, contig->configure->read_len, FLAG_SNP);
		int32_t *p = snpsregion->data, iterend = 0, nextposend, curr_end = 0, flag;
		hts_itr_t *iter = NULL;
		uint64_t curr_off = 0;
		flag = snps->length / 10 > 1000 ? snps->length / 10 : 1000;
		SeqList* snpslinkdata = seqlist_init(sizeof(KmerScore), flag);
		KmerScore* ks = ks_init(contig->configure->max_variant_count_lgs);
		uint8_t* pks = ks->region;
		SnpsList* snpslist = snpslist_init(snps);
		bam1_t *read = contig->read;
		Snps** psnps = snps->data;
		contig->read_fliter = contig_read_fliter;
		for (int i = 0; i < snpsregion->length; i += 2, p += 2) {
			nextposend = i + 2 < snpsregion->length ? *(p + 2) : -1;
			flag = 2;
			while (contig_next_iter(contig, contig->fp, contig->idx, &iter, *p, *(p + 1), &iterend, &curr_off, &read, &curr_end, &nextposend, flag) >= 0) {
				if (contig->read_fliter(contig, read) == 2) {
					ts_snps_parse_read(contig, read, *p, *(p + 1), 0, snpslinkdata, ks);
					ts_snps_deal_linkdata(contig, snpslinkdata, snpslist, 0);
					snpslinkdata->length = 0;
					ks->length = 0;
					ks->region = pks;
				}
				flag = 0;
			}
		}

		if (iter) hts_itr_destroy(iter);
		iter = NULL;
		psnps++;
		for (int i = 1; i < snps->length; i++, psnps++) {
			if ((*psnps)->total <= contig->configure->min_count_snp_link) {
				contig->data[(*(psnps - 1))->left]->flag |= FLAG_LEFT;
				contig->data[(*(psnps - 1))->pos]->flag |= FLAG_LEFT;
				contig->data[(*(psnps - 1))->right]->flag |= FLAG_RIGHT;
				contig->data[(*psnps)->left]->flag |= FLAG_LEFT;
				contig->data[(*psnps)->pos]->flag |= FLAG_RIGHT;
				contig->data[(*psnps)->right]->flag |= FLAG_RIGHT;
			}
		}
		seqlist_destory(snpsregion);
		snpsregion = ts_find_snp_region(contig, snps, contig->configure->max_variant_count_lgs, 0);

		contig->read_fliter = contig_read_fliter2;
		p = snpsregion->data;
		curr_end = 0;
		for (int i = 0; i < snpsregion->length; i += 2, p += 2) {
			nextposend = i + 2 < snpsregion->length ? *(p + 2) : -1;
			flag = 2;
			while (contig_next_iter(contig, contig->tfp, contig->tidx, &iter, *p, *(p + 1), &iterend, &curr_off, &read, &curr_end, &nextposend, flag) >= 0) {
				if (contig->read_fliter(contig, read) == 1) {
					ts_snps_parse_read(contig, read, *p, *(p + 1), 1, snpslinkdata, ks);
					ts_snps_deal_linkdata(contig, snpslinkdata, snpslist, 1);
					snpslinkdata->length = 0;
					ks->length = 0;
					ks->region = pks;
				}
				flag = 0;
			}
		}

		seqlist_destory(snpsregion);
		snpslist_destory(snpslist);
		ks_destory(ks);
		seqlist_destory(snpslinkdata);
		if (iter) hts_itr_destroy(iter);
	}
}

inline void ts_tranfer_link(Contig * contig, Snps ** snp, KmerScore ** ks, int16_t* total, int32_t flag)
{
	if (ks[0]->length == snp[0]->length&&ks[1]->length == snp[1]->length&&flag == 4) {
		int32_t index = snps_get_index(snp[0], ks[0]->region);
		if (index != -1) {
			index++;
			ks[1]->length = index << BASE_SHIFT;
			index = snps_get_index(snp[1], ks[1]->region);
			if (index != -1) {
				index++;
				ks[1]->length += index;
				KmerScore* p = seqlist_find(snp[1]->link, ks[1], ks_compare_length);
				if (p == NULL) {
					ks[1]->num = 1;
					seqlist_append(snp[1]->link, ks[1]);
				}
				else {
					p->num++;
					p->mapqual += ks[1]->mapqual;
					p->qual += ks[1]->qual;
				}
				(*total)++;
			}
		}
	}
}

void ts_snps_score(Contig * contig, SeqList * snps)
{
	if (snps->length > 1) {
		Snps** p = snps->data, *q = *p;
		KmerScore* ks;
		Score* pscore;
		double score;
		Base* base[SNP_NUM];
		uint16_t link[2][SNP_NUM + 1], temp;
		int32_t i = 1, j, k;
		base[0] = contig->data[(*p)->pos];
		base_clean_score(base[0]);
		for (; i <= SNP_NUM; i++) {
			base_add_score(base[0], i, 0);
		}
		for (i = 1, p++; i < snps->length; i++, p++) {
			q = *p;
			base[0] = contig->data[(*(p - 1))->pos];
			base[1] = contig->data[q->pos];
			base_clean_score(base[1]);
			if (q->link->length) {
				ks = q->link->data;
				memset(link, 0, sizeof(uint16_t) * 2 * (SNP_NUM + 1));
				for (j = 0; j < q->link->length; j++, ks++) {
					temp = ks->length >> BASE_SHIFT;
					score = base_get_score(base[0], temp)->score;
					score += ks->num*log10((ks->mapqual + ks->qual) / (double)ks->num + 2) - q->total / contig->configure->ploidy;
					pscore = base_get_score(base[1], ks->length);
					if (pscore == NULL || pscore->score < score) {
						if (link[0][temp]) {
							if (base_get_score(base[1], link[0][temp])->score >= score) {
								continue;
							}
							link[1][link[0][temp]] = 0;
						}
						if (pscore != NULL) {
							link[0][pscore->kmer >> BASE_SHIFT] = 0;
						}
						base_add_score(base[1], ks->length, score);
						link[0][temp] = ks->length & 0xf;
						link[1][ks->length & 0xf] = temp;
					}
				}

				k = 1;
				for (j = 1; j <= SNP_NUM; j++) {
					if (link[1][j] == 0) {
						for (; k <= SNP_NUM; k++) {
							if (link[0][k] == 0) {
								temp = (k << BASE_SHIFT) + j;
								score = base_get_score(base[0], k)->score;
								score -= q->total / contig->configure->ploidy;
								base_add_score(base[1], temp, score);
								break;
							}
						}
					}
				}
			}
			else {
				for (j = 1; j <= SNP_NUM; j++) {
					base_add_score(base[1], j, 0);
				}
			}
		}
	}
}

void ts_snps_correct(Contig * contig, SeqList * snps)
{
	if (snps->length > 1) {
		int32_t i = snps->length - 1, index;
		Snps** p = seqlist_index(snps, i), *q = NULL;
		Base* base = NULL;
		Score* score = NULL;
		for (; i > 0; i--, p--) {
			q = *p;
			if (q->link->length > 0) {
				if (score == NULL) {
					base = contig->data[q->pos];
					score = base_max_score(base);
					if (q->length == 1) {
						base->base = q->region[score->base - 1][0];
					}
					else {
						contig_update_contig(contig, q->pos, q->pos + 1, q->region[score->base - 1], -1);
					}
				}
				index = (score->kmer >> BASE_SHIFT) - 1;
				q = *(p - 1);
				base = contig->data[q->pos];
				if (q->length == 1) {
					base->base = q->region[index][0];
				}
				else {
					contig_update_contig(contig, q->pos, q->pos + 1, q->region[index], -1);
				}
				if (q->link->length > 0) {
					score = base_get_score(base, index + 1);
				}
				else {
					score = NULL;
				}
			}
		}
	}
}

SeqList * ts_find_snp_region(Contig * contig, SeqList* snps, int32_t gap, uint32_t flag)
{
	SeqList* result = seqlist_init(sizeof(int32_t), snps->length);
	int32_t i = 0, temp;
	uint32_t flag1 = FLAG_LEFT | FLAG_RIGHT, flag2;
	Snps** p = snps->data, *qstart = NULL, *qend = NULL;
	for (; i < snps->length; i++, p++) {
		flag2 = contig->data[(*p)->pos]->flag;
		if ((flag2&flag) || (flag2&flag1)) {
			if (qstart == NULL) {
				qend = qstart = *p;
			}
			else if (flag || (flag2&FLAG_RIGHT)) {
				if (flag) {
					temp = (*p)->pos - qend->pos;
				}
				else {
					temp = (*p)->right - qend->left;
				}
				if (temp < gap) {
					qend = *p;
				}
				else {
					if (qstart != qend) {
						if (flag) {
							temp = qend->pos + 1;
							seqlist_append(result, &qstart->pos);
							seqlist_append(result, &temp);
						}
						else {
							seqlist_append(result, &qstart->left);
							seqlist_append(result, &qend->right);
						}
					}
					if (flag || (flag2&FLAG_LEFT)) {
						qend = qstart = *p;
					}
					else {
						qend = qstart = NULL;
					}
				}
			}
		}
	}
	if (qstart && qstart != qend) {
		if (flag) {
			seqlist_append(result, &qstart->pos);
			seqlist_append(result, &qend->pos);
		}
		else {
			seqlist_append(result, &qstart->left);
			seqlist_append(result, &qend->right);
		}
	}
	return result;
}

void ts_snps_parse_read(Contig * contig, bam1_t * read, int32_t start, int32_t end, uint32_t flagbrim, SeqList* snpslinkdata, KmerScore* ks)
{
	if (read->core.n_cigar) {
		int32_t pos = read->core.pos, qpos = 0, qstart, qend, i, j, k, len, del = 0, curpos = 0, sign = 0, comfirmindex = 0, totallength = contig->configure->max_variant_count_lgs;
		uint8_t curcigar, lastcigar = BAM_CINS, *q = ks->region;
		uint32_t *cigar = bam_get_cigar(read);
		uint16_t flag = FLAG_LEFT | FLAG_RIGHT;
		uint8_t* seq = bam_get_seq(read);
		uint8_t* qual = bam_get_qual(read);
		KmerScore* comfirm = NULL;
		SeqList* p;
		Base* base;
		contig_cut_read(contig, read, &qstart, &qend);
		i = 0;
		ks->mapqual = read->core.qual;
		for (; i < read->core.n_cigar; ++i) {
			len = bam_cigar_oplen(cigar[i]);
			curcigar = bam_cigar_op(cigar[i]);
			if (totallength - len < 0) {
				trace_log(TRACE_ERROR, "[E:over max_variant_count_lgs] read: %s\n",bam_get_qname(read));
				break;
			}
			switch (curcigar) {
			case BAM_CMATCH:
			case BAM_CDEL:
				for (j = 0; j < len; j++, pos++) {
					if (pos >= start && pos <= end && qpos >= qstart && qpos <= qend) {
						base = contig->data[pos];
						if (lastcigar != BAM_CINS && pos > start && (qpos > qstart || (qpos == qstart && lastcigar == BAM_CDEL))) {
							p = contig->data[pos - 1]->insert;
							if (p != NULL && curpos) {
								for (k = 0; k < p->length; k++) {
									ks_append(ks, BASE_DEL);
									totallength--;
									del++;
								}
							}
						}
						if (flagbrim == 0 || (base->flag&flag)) {
							if (base->flag&FLAG_SNP) {
								if (curpos == 0) {
									ks->region = q;
									ks->length = 0;
									ks->num = pos;
									ks->qual = 0;
									del = 0;
									curpos = 1;
									if (flagbrim == 0) sign = 1;
								}
								else {
									sign++;
								}
							}
							else if (flagbrim) {
								if (base->base == bam_seqi(seq, qpos)) {
									sign++;
								}
							}
							else {
								sign++;
							}
							if (curpos) {
								if (curcigar == BAM_CDEL) {
									ks_append(ks, BASE_DEL);
								}
								else {
									ks_append(ks, bam_seqi(seq, qpos));
									ks->qual += qual[qpos];
								}
								totallength--;
								if (ks->num != pos || base->insert == NULL) {
									if (ks->num != pos) {
										if (ks->length != del) {
											ks->qual /= (double)(ks->length - del);
										}
										else {
											ks->qual = 0;
										}
										ks->length--;
									}
									seqlist_append(snpslinkdata, ks);
									q += ks->length;
									curpos = 0;
								}
							}
							if (ks->num != pos) {
								if (base->flag&FLAG_SNP) {
									ks->region = q;
									ks->length = 1;
									ks->num = pos;
									ks->qual = qual[qpos];
									del = 0;
									curpos = 1;
									if (flagbrim == 0) {
										comfirmindex++;
										sign = 1;
									}
								}
								else if (base->flag&FLAG_RIGHT) {
									if (sign == 2) {
										comfirmindex = snpslinkdata->length;
									}
									else {
										comfirm = seqlist_index(snpslinkdata, comfirmindex);
										if (comfirm != NULL) {
											for (k = comfirmindex; k < snpslinkdata->length; k++, comfirm++, comfirmindex++) {
												comfirm->length = 0;
											}
										}
									}
									curpos = 0;
									sign = 0;
									if (base->flag&FLAG_LEFT) sign++;
								}

							}
						}

					}
					if (curcigar != BAM_CDEL) {
						qpos++;
					}
					lastcigar = curcigar;
				}
				break;
			case BAM_CINS:
				if (curpos) {
					if (pos){
						p = contig->data[pos - 1]->insert;
						for (j = 0; j < len; j++, qpos++) {
							if (pos > start && pos <= end && qpos >= qstart && qpos <= qend) {
								ks_append(ks, bam_seqi(seq, qpos));
								ks->qual += qual[qpos];
								totallength--;
							}
						}
						if (pos > start && pos <= end && qpos > qstart && qpos <= qend + 1) {
							for (; j < p->length; j++) {
								ks_append(ks, BASE_DEL);
								totallength--;
								del++;
							}
						}
					}else{
						qpos += len;
						qstart += len;
						lastcigar = curcigar;
					}
				}
				else {
					qpos += len;
				}
				lastcigar = curcigar;
				break;
			case BAM_CHARD_CLIP:
			case BAM_CSOFT_CLIP:
				qpos += len;
				break;
			}
		}
	}
}

void ts_snps_deal_linkdata(Contig * contig, SeqList * linkdata, SnpsList * snpslist, uint32_t flag)
{
	if (linkdata->length > 1) {
		int32_t i = 1, index;
		KmerScore* ks[2], *p = linkdata->data;
		Snps* snps[2];
		for (p++; i < linkdata->length; i++, p++) {
			if (p->length && (p - 1)->length && (flag == 0 || ((contig->data[p->num]->flag&FLAG_RIGHT) && (contig->data[(p - 1)->num]->flag&FLAG_LEFT)))) {
				index = snpslist_find(snpslist, p->num);
				snps[0] = snpslist->data[index - 1];
				snps[1] = snpslist->data[index];
				ks[0] = p - 1;
				ks[1] = p;
				ts_tranfer_link(contig, snps, ks, &snps[1]->total, 4);
			}
		}
	}
}

inline void ts_correct_lower_depth(Contig * contig, SeqList* nodepth)
{
	int32_t *p = nodepth->data, iterend = 0, nextposend, curr_end = 0, flag;
	hts_itr_t *iter = NULL;
	uint64_t curr_off = 0;
	bam1_t *read = contig->read;
	contig->read_fliter = contig_read_fliter;
	for (int i = 0; i < nodepth->length; i += 2, p += 2) {
		nextposend = i + 2 < nodepth->length ? *(p + 2) : -1;
		flag = 2;
		contig_clean_region(contig, *p, *(p + 1));
		contig_as_read(contig, *p, *(p + 1));
		while (contig_next_iter(contig, contig->fp, contig->idx, &iter, *p, *(p + 1), &iterend, &curr_off, &read, &curr_end, &nextposend, flag) >= 0) {
			if (contig->read_fliter(contig, read) == 2) {
				contig_parse_read(contig, read, *p, *(p + 1), BASE_SHIFT);
			}
			flag = 0;
		}
	}

	curr_off = 0;
	curr_end = 0;
	if (iter) hts_itr_destroy(iter);
	iter = NULL;
	contig->read_fliter = contig_read_fliter2;
	p = nodepth->data;
	for (int i = 0; i < nodepth->length; i += 2, p += 2) {
		nextposend = i + 2 < nodepth->length ? *(p + 2) : -1;
		flag = 2;
		while (contig_next_iter(contig, contig->tfp, contig->tidx, &iter, *p, *(p + 1), &iterend, &curr_off, &read, &curr_end, &nextposend, flag) >= 0) {
			if (contig->read_fliter(contig, read) == 1) {
				contig_parse_read(contig, read, *p, *(p + 1), BASE_SHIFT);
			}
			flag = 0;
		}
	}

	p = nodepth->data;
	for (int i = 0; i < nodepth->length; i += 2, p += 2) {
		contig_region_score(contig, *p, *(p + 1),contig->configure->indel_balance_factor_lgs);
		ts_region_correct(contig, *p, *(p + 1));
	}

	if (iter) hts_itr_destroy(iter);
}

void ts_region_correct(Contig * contig, int32_t start, int32_t end)
{
	Base* base = contig->data[end];
	Score* score = base_max_score(base);
	int32_t i = end, j = 0;
	Kmer* maxn[2];
	double rate;
	while (i > start || (i == start && j == 0))
	{
		if (base->flag&FLAG_ZERO || (j == 0 && score->base != BASE_DEL)) {
			base->base = score->base;
		}
		base_merge_kmer(base);
		if (base->data->length >= 2) {
			base_get_nlargest(base, maxn, 2);
			rate = maxn[1]->count / (double)maxn[0]->count;
			if (maxn[0]->kmer != base->base || rate > contig->configure->max_indel_factor_lgs) {
				if (base->base == BASE_DEL || j != 0 || maxn[0]->kmer != base->base || rate > contig->configure->max_snp_factor_lgs) {
					base->flag |= FLAG_THIRD;
				}
				else {
					base->flag &= FLAG_THIRD_N;
				}
			}
		}
		base = contig_data_pre(contig, &i, &j);
		score = base_get_score(base, score->kmer >> BASE_SHIFT);
	}
}

inline int32_t ts_get_nlargest(SeqList * regiondata, KmerScore ** maxn, int32_t n)
{
	int32_t count = 0;
	KmerScore* p = regiondata->data;
	maxn[0] = p;
	if (regiondata->length > 0) {
		p++;
		count++;
		int i = 1, j = 0;
		for (; i < regiondata->length; i++, p++) {
			for (j = count - 1; j >= 0; j--) {
				if (ks_compare(p, maxn[j]) > 0) {
					if (j < n - 1) {
						maxn[j + 1] = maxn[j];
					}
					maxn[j] = p;
				}
				else {
					if (j < n - 1) {
						maxn[j + 1] = p;
					}
					break;
				}
			}
			if (count < n) {
				count++;
			}
		}
	}
	return count;
}
