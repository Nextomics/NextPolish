#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdint.h>
#include "bseq.h"

uint8_t nt_table[128] = {
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

uint8_t bs_table[] = {
	65, 67, 71, 84
};

char **init_bases(int n){
	int i, j;
	char **arr = (char **)malloc(sizeof(char *) * n + sizeof(char) * 9 * n);
	char * const _ = (char *) (arr + n);
	for(i = 0; i < n; i++) {
		arr[i] = _ + i * 9;
		for(j = 0; j < 8; j++)
			arr[i][j] = "ACGT"[i >> (14 - (j << 1)) & 3];
		arr[i][8] = '\0';
	}
	return arr;
}

void destroy_bases(char **arr){
	free (arr);
}

int init_seq_mode(FILE * fp){
	uint8_t buf[2] = {0, 254};
	fwrite(buf, sizeof(uint8_t), 2, fp);
	return 2;
}

int find_seq_mode(gzFile f){
	int mode = 1;//1:fa/fq/gz; 2:bit

	uint8_t buf[2];
	int bread = gzread (f, buf, 2);
	if (bread == 2){
		if (buf[0] == 0 || buf[0] == 254) return 2;
		gzrewind (f);
	}else{
		fprintf(stderr, "[ERROR] failed to find seq file type or empty file\n");
		exit(EXIT_FAILURE);
	}
	return mode;
}

uint32_t seq2bit(FILE *fp, uint32_t name, uint32_t len, char *seq){// pls init_seq_mode first
	uint32_t start = 0;
	name = (uint32_t) name;

	uint32_t i, base_cnt = 0;
	uint64_t buffer = 0;
	fwrite(&name, sizeof(uint32_t), 1, fp);
	fwrite(&len, sizeof(uint32_t), 1, fp);
	start += 8;
	for(i = 0; i < len; i++) {
		buffer = buffer << 2 | nt_table[(uint8_t)seq[i]];
		++base_cnt;
		if (base_cnt == 16) {
			fwrite(&buffer, 4, 1, fp);
			start += 4;
			base_cnt = 0;
			buffer = 0;
		}
	}
	if (base_cnt) {
		buffer <<= (32 - (base_cnt << 1));
		fwrite(&buffer, 4, 1, fp);
		start += 4;
	}
	return start;
}

void seq2bit1(uint32_t *s, uint32_t len, char *seq){
	uint32_t i, j, base_cnt = 0;
	uint32_t buffer = 0;
	for(i = j = 0; i < len; i++) {
		buffer = buffer << 2 | nt_table[(uint8_t)seq[i]];
		++base_cnt;
		if (base_cnt == 16) {
			s[j++] = buffer;
			base_cnt = 0;
			buffer = 0;
		}
	}
	if (base_cnt) {
		buffer <<= (32 - (base_cnt << 1));
		s[j++] = buffer;
	}
}

void bit2seq1(uint32_t *s, uint32_t len, char *seq){
	
	uint32_t i, j = len - 1;
	int32_t l = ((len - 1) >> 4) + 1;
	int i_m, l_shift = ((l << 4) - len) << 1;
	while (l--){
		i = s[l];
		i_m = 32;
		if (l_shift){
			i_m -= l_shift;
			i = i >> l_shift;
			l_shift = 0;
		}
		for (; i_m; i_m -= 2){
			seq[j--] = "ACGT"[i & 3];
			i >>= 2;
		}
	}
	seq[len] = '\0';
}

void realloc_sbbuf(sbbuf *buf, uint32_t len){
	buf->m = len + 2;
	buf->seq = (char *) realloc(buf->seq, buf->m);
	buf->buffer = (char *) realloc(buf->buffer, buf->m + 32);
	buf->tmp = (uint32_t *) realloc(buf->tmp, (buf->m + 32) <<2);
}

sbbuf *init_sbbuf(uint32_t len){
	sbbuf *buf = (sbbuf *) malloc(sizeof(sbbuf));
	buf->bases_arr = init_bases(65536);
	buf->m = len + 2;
	buf->seq = (char *) malloc(buf->m);
	buf->buffer = (char *) malloc(buf->m + 32);
	buf->tmp = (uint32_t *) malloc((buf->m + 32) << 2);
	return buf;
}

void destroy_sbbbuf(sbbuf *buf){
	free (buf->seq);
	free (buf->buffer);
	free (buf->tmp);
	free (buf->bases_arr);
	free (buf);
}

void subfa(sbbuf *buf, FILE *fp, int rev, uint64_t offset, uint32_t start, uint32_t end)
{	
	uint32_t i;
	int l = end - start + 1;
	if (l + 1 >= buf->m) realloc_sbbuf(buf, end - start);
	fseeko(fp, offset + start, SEEK_SET);
	fread(buf->seq, 1, l, fp);
	if (rev) {
		char c1, c2;
		for(i = 0; i < l >> 1; ++i) {
			c1 = "ACGT"[3 - nt_table[(uint8_t)buf->seq[i]]];
			c2 = "ACGT"[3 - nt_table[(uint8_t)buf->seq[l - 1 - i]]];
			buf->seq[i] = c2; buf->seq[l - 1 - i] = c1;
		}
		if (l & 1)  buf->seq[l >> 1] = "ACGT"[3 - nt_table[(uint8_t)buf->seq[l >> 1]]];
	}
	buf->seq[end - start + 1] = '\0';
}

void subbit(sbbuf *buf, FILE *fp, int rev, uint64_t offset, uint32_t start, uint32_t end)
{	
	if (end - start  + 2 >= buf->m) realloc_sbbuf(buf, end - start);
	uint32_t start_align = (start >> 4) << 4;
	uint32_t end_align = ((end >> 4) + 1) << 4;
	uint32_t i, cnt = (end_align - start_align) >> 4;
	fseeko(fp, offset + (start_align >> 2), SEEK_SET);//here is the I/O bottleneck.
	fread(buf->tmp, 4, cnt, fp);
	for(i = 0; i < cnt; i++) {
		memcpy(buf->buffer + (i << 4), buf->bases_arr[buf->tmp[i] >> 16 & 0xFFFF], 8);
		memcpy(buf->buffer + (i << 4) + 8, buf->bases_arr[buf->tmp[i] & 0xFFFF], 8);
	}
	memcpy(buf->seq, buf->buffer + start - start_align, end - start + 1);
	if (rev) {
		char c1, c2;
		int l = end - start + 1;
		for(i = 0; i < l >> 1; ++i) {
			c1 = "ACGT"[3 - nt_table[(uint8_t)buf->seq[i]]];
			c2 = "ACGT"[3 - nt_table[(uint8_t)buf->seq[l - 1 - i]]];
			buf->seq[i] = c2; buf->seq[l - 1 - i] = c1;
		}
		if (l & 1)  buf->seq[l >> 1] = "ACGT"[3 - nt_table[(uint8_t)buf->seq[l >> 1]]];
	}
	buf->seq[end - start + 1] = '\0';
}

void subbit_(sbbuf *buf, uint32_t *fp, int rev, uint64_t offset, uint32_t start, uint32_t end)
{
	if (end - start  + 2 >= buf->m) realloc_sbbuf(buf, end - start);
	uint32_t start_align = (start >> 4) << 4;
	uint32_t end_align = ((end >> 4) + 1) << 4;
	uint32_t i, cnt = (end_align - start_align) >> 4;
	memcpy(buf->tmp, fp + ((offset + (start_align >> 2)) >> 2), 4 * cnt);
	for(i = 0; i < cnt; i++) {
		memcpy(buf->buffer + (i << 4), buf->bases_arr[buf->tmp[i] >> 16 & 0xFFFF], 8);
		memcpy(buf->buffer + (i << 4) + 8, buf->bases_arr[buf->tmp[i] & 0xFFFF], 8);
	}
	memcpy(buf->seq, buf->buffer + start - start_align, end - start + 1);
	if (rev) {
		char c1, c2;
		int l = end - start + 1;
		for(i = 0; i < l >> 1; ++i) {
			c1 = "ACGT"[3 - nt_table[(uint8_t)buf->seq[i]]];
			c2 = "ACGT"[3 - nt_table[(uint8_t)buf->seq[l - 1 - i]]];
			buf->seq[i] = c2; buf->seq[l - 1 - i] = c1;
		}
		if (l & 1)  buf->seq[l >> 1] = "ACGT"[3 - nt_table[(uint8_t)buf->seq[l >> 1]]];
	}
	buf->seq[end - start + 1] = '\0';
}

int kbit_read(kseq_t *seq){

	static char **bases_arr = NULL;
	if (!seq->seq.l){
		if (!seq->name.s){
			seq->name.s = (char *) calloc(1, 20); // only accept uint32_t name, sizeof(uint32_t) == 4
			seq->seq.m = 100000;
			seq->seq.s = (char *) malloc(seq->seq.m);
		}
	}
	if (!bases_arr) bases_arr = init_bases(65536);
	//Parse sequence name
	uint32_t name;
	if(gzread(seq->f->f, &name, 4) <= 0) {
		if (bases_arr) {
			destroy_bases(bases_arr);
			bases_arr = NULL;
		}
		seq->f->is_eof = 1;
		seq->seq.l = 0; //for multi-files share var.(seq)
		return EOF;
	}
	seq->offset += 4;
	seq->name.l = sprintf(seq->name.s, "%u", name);

	//Parse sequence length
	gzread(seq->f->f, &seq->seq.l, 4);
	if (seq->seq.l + 16 >= seq->seq.m){
		seq->seq.m = seq->seq.l + 17;
		seq->seq.s = (char*)realloc(seq->seq.s, seq->seq.m);
	}
	seq->offset += 4;
	//Parse sequence, convert 2bit to bases
	uint32_t end = (((seq->seq.l - 1) >> 4) + 1) << 4;
	uint32_t i, cnt = end >> 4;
	uint32_t *tmp = (uint32_t *)malloc(end);
	gzread(seq->f->f, tmp, cnt << 2);
	for(i = 0; i < cnt; i++) {
		memcpy(seq->seq.s + (i << 4), bases_arr[tmp[i] >> 16 & 0xFFFF], 8);
		memcpy(seq->seq.s + (i << 4) + 8, bases_arr[tmp[i] & 0xFFFF], 8);
	}
	seq->seq.s[seq->seq.l] = '\0';
	seq->offset += (end >> 2);
	free(tmp);
	return seq->seq.l;
}

int kseq_r(kseq_t *seq)
{
	if (!seq->fm) seq->fm = find_seq_mode(seq->f->f);
	if (seq->fm == 1) return kseq_read(seq);
	else return kbit_read(seq);
}