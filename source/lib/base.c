#include "base.h"

extern const Configure configure;

const char basetostr[] = "=ACMGRSVTWYHKDBN";
const uint8_t strtobase[] = { 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
							  15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
							  15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
							  15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
							  15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
							  15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
							  15, 0, 15, 15, 15, 1, 14, 2, 13, 15,
							  15, 4, 11, 15, 15, 12, 15, 3, 15, 15,
							  15, 15, 5, 6, 8, 15, 7, 9, 15, 10,
							  15, 15, 15, 15, 15, 15, 15, 15, 15, 15 };

inline Base * base_init(uint32_t size)
{
	Base* result = (Base*)malloc(sizeof(Base));
	if (result == NULL) {
		fprintf(stderr, "base data malloc failed,the memory is full!\n");
		exit(1);
	}
	result->base = 3;
	result->flag = 0;
	result->refkmer = 0;
	result->count = 0;
	result->insert = NULL;
	result->data = seqlist_init(size, 1);
	result->score = seqlist_init(sizeof(Score), 1);
	return result;
}

inline void base_destory(Base * base)
{
	if (base->insert != NULL) {
		Base* p = base->insert->data;
		for (int i = 0; i < base->insert->length; i++, p++) {
			seqlist_destory(p->data);
			seqlist_destory(p->score);
		}
		seqlist_destory(base->insert);
	}
	seqlist_destory(base->data);
	seqlist_destory(base->score);
	free(base);
}

int32_t comparekmer(void * first, void * second)
{
	if (((Kmer*)first)->kmer == ((Kmer*)second)->kmer) {
		return 0;
	}
	else if (((Kmer*)first)->kmer < ((Kmer*)second)->kmer) {
		return -1;
	}
	return 1;
}

void base_add_data(Base * base, uint16_t kmer)
{
	Kmer temp = { .kmer = kmer,.count = 1 };
	Kmer* result = seqlist_find(base->data, &temp, comparekmer);
	if (result == NULL) {
		seqlist_append(base->data, &temp);
	}
	else {
		result->count++;
	}
	base->count++;
}

void base_clean_data(Base * base)
{
	base->data->length = 0;
	base->count = 0;
}

double base_get_coverage(Base * base, uint16_t kmer)
{
	Kmer* p = base->data->data;
	uint32_t count = 0;
	for (int i = 0; i < base->data->length; i++, p++) {
		if ((p->kmer & 0xf) == kmer) {
			count += p->count;
		}
	}
	return count / (double)base->count;
}

int32_t base_get_nlargest(Base * base, Kmer ** maxn, int32_t n)
{
	int32_t count = 0;
	if (base->data->length > 0) {
		Kmer* p = base->data->data;
		maxn[0] = p;
		p++;
		count++;
		int i = 1, j = 0;
		for (; i < base->data->length; i++, p++) {
			for (j = count - 1; j >= 0; j--) {
				if (p->count > maxn[j]->count) {
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

void base_merge_kmer(Base * base)
{
	const int N = 16;
	int32_t count = 0, index = 0;
	Kmer* temp[N];
	if (base->data->length > 0) {
		memset(temp, 0, sizeof(Kmer*) * N);
		Kmer *p = base->data->data, *q = p;
		for (int i = 0; i < base->data->length; i++, p++) {
			index = p->kmer & 0xf;
			if (temp[index] == 0) {
				count++;
				q->kmer = index;
				q->count = p->count;
				temp[index] = q;
				q++;
			}
			else {
				temp[index]->count += p->count;
			}
		}
		base->data->length = count;
	}
}

int32_t comparescore(void * first, void * second)
{
	if (((Score*)first)->base == ((Score*)second)->base) {
		return 0;
	}
	else if (((Score*)first)->base < ((Score*)second)->base) {
		return -1;
	}
	return 1;
}

void base_add_score(Base * base, uint16_t kmer, double score)
{
	Score temp = { .base = kmer & 0xf,.kmer = kmer,.score = score };
	Score* result = seqlist_find(base->score, &temp, comparescore);
	if (result == NULL) {
		seqlist_append(base->score, &temp);
	}
	else {
		seqlist_update(base->score, result, &temp);
	}
}

Score* base_get_score(Base * base, uint16_t kmer)
{
	if (kmer) {
		Score temp = { .base = kmer & 0xf,.kmer = kmer,.score = 0 };
		return seqlist_find(base->score, &temp, comparescore);
	}
	return base_max_score(base);
}

void base_clean_score(Base * base)
{
	base->score->length = 0;
}

inline Score * base_max_score(Base * base)
{
	Score* p = NULL, *q = NULL;
	if (base->score->length != 0) {
		q = p = base->score->data;
		for (int i = 0; i < base->score->length; i++, p++) {
			if (p->score > q->score) {
				q = p;
			}
		}
	}
	return q;
}

void base_set_flag(Base * base, uint8_t flag)
{
	base->flag |= flag;
}
