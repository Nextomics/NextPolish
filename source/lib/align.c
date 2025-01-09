#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef LGS_CORRECT
	#include "ctg_cns.h"
#else
	#include "nextcorrect.h"
#endif


void malloc_vd(int **V, uint8_t ***D, uint64_t max_mem_d){
	*V = malloc( max_mem_d * 2 * sizeof(int));
	int i;
	for (i = 0; i < 3 && ! *D; i++){
		*D = malloc( max_mem_d * sizeof(uint8_t *) + max_mem_d * (max_mem_d + 1)/2 * sizeof(uint8_t));
	}
	if (*D) {
		uint64_t d;
		uint8_t * const _ = (uint8_t *) (*D + max_mem_d);
		for (d = 0; d < max_mem_d; d ++ ) {
			(*D)[d] = _ + d * (d + 1)/2;
		}
	}
}

void reverse_str(char *str, int len)
{
	char tmp;
	char *p1 = str;
	char *p2 = str + len - 1;
	while (p1 < p2) {
		tmp = *p1;
		*p1++ = *p2;
		*p2-- = tmp;
	}
}

void align(char *query_seq, int q_len,
				char *target_seq, int t_len,
				alignment *align_rtn, int *V, uint8_t **D) {

	int max_d = (int) (0.4 * (q_len + t_len));

	float band_factor = q_len + t_len > 5000 ? 0.1 : 1;
	int band_size =  (int) (band_factor * (q_len + t_len));

	int x, new_min_k, new_max_k, k, k2, d;
	int y = x = 0;
	int kk = 0;
	int min_k = 0;
	int max_k = 0;
	int best_m = -1;
	bool aligned = false;
	int k_offset = max_d;

	for (d = 0; d < max_d && max_k - min_k <= band_size; d ++ ) {
		for (k = min_k; k <= max_k;  k += 2) {
			kk = k < 0 ? -1 * k - 1 : k;
			if ( (k == min_k) || ((k != max_k) && (V[ k - 1 + k_offset ] < V[ k + 1 + k_offset])) ) {
				x = V[ k + 1 + k_offset];
				D[d][kk] = 0;
			} else {
				x = V[ k - 1 + k_offset] + 1;
				D[d][kk] = 1;
			}

			y = x - k;
			while ( x < q_len && y < t_len && query_seq[x] == target_seq[y] ){
				x++;
				y++;
			}

			V[ k + k_offset] = x;

			if ( x + y > best_m) {
				best_m = x + y;
			}

			// if ( x >= q_len || y >= t_len) {//for prefix aln
			// 	aligned = true;
			// 	break;
			// }
			if ( x >= q_len && y >= t_len) { //for global aln
				aligned = true;
				break;
			}
		}
		// For banding
		new_min_k = max_k;
		new_max_k = min_k;

		k2 = min_k;
		while (k2 < new_min_k ){
			if (V[k2 + k_offset] * 2 - k2 >= best_m - 150) new_min_k = k2;
			k2 += 2;
		}
		
		k2 = max_k;
		while (k2 > new_max_k ){
			if (V[k2 + k_offset] * 2 - k2 >= best_m - 150) new_max_k = k2;
			k2 -= 2;
		}

		max_k = new_max_k + 1;
		min_k = new_min_k - 1;

		if (aligned == true){
			x--; // 0-based
			max_d = d;
			align_rtn->aln_t_e = align_rtn->aln_t_s + y - 1;// fix bug for tail with clip
			align_rtn->aln_t_len = y;
			align_rtn->aln_q_len = x + 1;

			int gap, aln_pos, pre_d, pre_k, pre_kk, pre_x, pre_y;
			gap = aln_pos = 0;

			while (true){
				while (x >= 0 && x >= k && query_seq[x] == target_seq[x - k]){
					align_rtn->t_aln_str[aln_pos] = align_rtn->q_aln_str[aln_pos] = query_seq[x];
					x --;
					aln_pos ++;
					gap = 0;
				}

				pre_d = d - 1;

				if (x < 0 && x - k < 0) break;
				if (D[d][kk]){
					pre_k = k - 1;
					pre_x = x - 1;
				}else{
					pre_k = k + 1;
					pre_x = x;
				}

				pre_y = pre_x - pre_k;
				pre_kk = pre_k < 0 ? -1 * pre_k - 1 : pre_k;

				if (pre_x == x && pre_y !=  x - k){ //advance in y

					if (x - k < 0) gap = 260;
					else{
						align_rtn->q_aln_str[aln_pos] = '-';
						align_rtn->t_aln_str[aln_pos++] = target_seq[x - k];
					}
					

				}else{ //advance in x
					
					if (x < 0) gap = 260;
					else{
						align_rtn->q_aln_str[aln_pos] = query_seq[x];
						align_rtn->t_aln_str[aln_pos ++] = '-';
					}
				}
				
				if (gap ++ > 250){//only allow the max length of a gap = 250
					aln_pos = 2;
					break;
				}

				d = pre_d;
				k = pre_k;
				kk = pre_kk;
				x = pre_x;
			}
			align_rtn->aln_len = aln_pos;
			reverse_str(align_rtn->t_aln_str, align_rtn->aln_len);
			reverse_str(align_rtn->q_aln_str, align_rtn->aln_len);
			align_rtn->t_aln_str[align_rtn->aln_len] = '\0';
			align_rtn->q_aln_str[align_rtn->aln_len] = '\0';
			break;
		}

	}
}



/////////////////////NM ALIGN/////////////////
static uint8_t MMH[] = {64, 16, 4, 1};
static uint8_t INS[] = {128, 32, 8, 2};
static uint8_t DEL[] = {192, 48, 12, 3};

typedef struct {
	int mas; //2
	int mis; //-4
	int gos; //-4
	int ges; //-2
} sco;

void align_nd_tb(const char *s1, const char *s2, uint32_t s1_l, uint32_t s2_l, uint8_t**d, alignment *aln){
	uint8_t i, j;
	aln->aln_len = 0;
	while(s1_l != 0 || s2_l != 0){
		i = s2_l&3;
		j = d[s1_l][s2_l>>2] & DEL[i];

		// printf("%d %d %d %d\n", s1_l, s2_l, i, aln->aln_len);
		if (j == MMH[i]){
				aln->t_aln_str[aln->aln_len] = s1[--s1_l];
				aln->q_aln_str[aln->aln_len] = s2[--s2_l];
				aln->aln_len ++;
		}else if (j == DEL[i]){
			aln->q_aln_str[aln->aln_len] = s2[--s2_l];
			aln->t_aln_str[aln->aln_len] = '-';
			aln->aln_len ++;
		}else{
			aln->t_aln_str[aln->aln_len] = s1[--s1_l];
			aln->q_aln_str[aln->aln_len] = '-';
			aln->aln_len ++;
		}
	}
	aln->t_aln_str[aln->aln_len] = aln->q_aln_str[aln->aln_len] = '\0';
	reverse_str(aln->t_aln_str, aln->aln_len);
	reverse_str(aln->q_aln_str, aln->aln_len);
}

void align_nd_dp(const char *s1, const char *s2, const uint32_t s1_l, const uint32_t s2_l, uint8_t**d){
	uint32_t i, j;
	sco score = {
		.mas = 2,
		.mis = -4,
		.gos = -4,
		.ges = -2
	};

	int32_t cs, ms, is, ds, *s = malloc((s2_l + 1) * sizeof(int32_t));
	s[0] = cs = 0;
	d[0][0] |= MMH[0];
	for (j = 1; j < s2_l + 1; j++){
		s[j] = s[j - 1] + ((d[0][(j - 1)>>2] & DEL[(j - 1)&3]) == DEL[(j - 1)&3] ? score.ges : score.gos);
		d[0][j>>2] |= DEL[j&3];
	}
	
	register int k;
	for (i = 1; i < s1_l + 1; i++){
		for (j = 0; j < s2_l + 1; j++){
			if (j == 0){
				cs = s[j] + ((d[i - 1][0]& DEL[0]) == INS[0] ? score.ges : score.gos);
				d[i][j] = INS[0];
			}else{
				k = j&3;
				ms = s[j - 1] + (s1[i - 1] == s2[j - 1] ? score.mas : score.mis);
				is = s[j] + ((d[i - 1][j>>2] & DEL[k]) == INS[k] ? score.ges : score.gos);
				ds = cs + ((d[i][(j - 1)>>2] & DEL[(j - 1)&3]) == DEL[(j - 1)&3] ? score.ges : score.gos);
				s[j - 1] = cs;

				d[i][j>>2] |= ms > is ? (ms > ds ? (cs = ms, MMH[k]) : (cs = ds, DEL[k])) : \
					(is > ds ? (cs = is, INS[k]) : (cs = ds, DEL[k]));
				if (j == s2_l) s[j] = cs;
			}
			// printf("%d %d %d %d\n",i,j,cs, d[i][j>>2] & DEL[j&3]);
		}
	}
	// printf("cs:%d\n",cs );
	free (s);
}

void align_nd(const char *s2, const uint32_t s2_l, const char *s1, const uint32_t s1_l, alignment *aln){
	uint64_t s2_size = (s2_l >> 2) + 1;
	uint64_t size = s2_size * (s1_l + 1) * sizeof(uint8_t) + (s1_l + 1) * sizeof(uint8_t *);

	uint8_t **d = calloc(1, size);
	// assert (d);
	d[0] = (uint8_t *)(d + s1_l + 1);
	uint32_t i;
	for (i = 1; i < s1_l + 1; i++ ) d[i] = d[i - 1] + s2_size;

	// alignment aln;
	// aln.t_aln_str = malloc(s1_l + s2_l);
	// aln.q_aln_str = malloc(s1_l + s2_l);
	align_nd_dp(s1, s2, s1_l, s2_l, d);
	align_nd_tb(s1, s2, s1_l, s2_l, d, aln);
	// printf("%s\n%s\n", aln.t_aln_str, aln.q_aln_str);
	free (d);
	// free (aln.t_aln_str);
	// free (aln.q_aln_str);
}
