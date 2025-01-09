#include "seqlist.h"

inline SeqList * seqlist_init(uint32_t stepsize, uint32_t max)
{
	SeqList* result = (SeqList*)calloc(sizeof(SeqList), 1);
	result->stepsize = stepsize;
	result->maxsize = max;
	if (max != 0) {
		result->data = calloc(sizeof(uint8_t), stepsize*max);
		if (result->data == NULL) {
			printf("seqlist malloc failed,the memory is full!\n");
			exit(1);
		}
	}
	return result;
}

inline void seqlist_destory(SeqList * seqlist)
{
	if (seqlist->data != NULL) {
		free(seqlist->data);
	}
	free(seqlist);
}

inline void seqlist_append(SeqList * seqlist, void * item)
{
	if (seqlist->length == seqlist->maxsize) {
		seqlist->maxsize += seqlist->maxsize / 2 + 1;
		if (seqlist->data == NULL) {
			seqlist->data = malloc(seqlist->maxsize*seqlist->stepsize);
		}
		else {
			seqlist->data = realloc(seqlist->data, seqlist->maxsize*seqlist->stepsize);
		}
		if (seqlist->data == NULL) {
			fprintf(stderr, "seqlist malloc failed,the memory is full!\n");
			exit(1);
		}
	}
	void* dest = (uint8_t*)seqlist->data + seqlist->stepsize*seqlist->length;
	memcpy(dest, item, seqlist->stepsize);
	seqlist->length++;
}

inline void seqlist_remove(SeqList * seqlist, uint32_t pos)
{
	if (seqlist->length > pos + 1) {
		uint8_t* dest = (uint8_t*)seqlist->data + seqlist->stepsize*pos;
		pos++;
		uint8_t* src = (uint8_t*)seqlist->data + seqlist->stepsize*pos;
		//memcpy(dest, src, seqlist->stepsize*(seqlist->length - pos - 1));
		for (; pos < seqlist->length; pos++) {
			memcpy(dest, src, seqlist->stepsize);
			dest += seqlist->stepsize;
			src += seqlist->stepsize;
		}
	}
	seqlist->length--;
}

inline void * seqlist_index(SeqList * seqlist, uint32_t pos)
{
	if (pos >= seqlist->length) {
		return NULL;
	}
	void* temp = (uint8_t*)seqlist->data + seqlist->stepsize*pos;
	return temp;
}

inline void seqlist_setitem(SeqList * seqlist, void * item, const uint32_t pos)
{
	if (seqlist->length > pos) {
		void* dest = (uint8_t*)seqlist->data + seqlist->stepsize*pos;
		memcpy(dest, item, seqlist->stepsize);
	}
}

void seqlist_update(SeqList * seqlist, void * dest, void * item)
{
	memcpy(dest, item, seqlist->stepsize);
}

inline void * seqlist_find(SeqList * seqlist, void * item, equalto comparefunc)
{
	void* result = NULL;
	for (uint16_t i = 0; i < seqlist->length; i++) {
		void* curptr = (uint8_t*)seqlist->data + seqlist->stepsize*i;
		if (comparefunc == NULL) {
			if (memcmp(curptr, item, seqlist->stepsize) == 0) {
				result = curptr;
				break;
			}
		}
		else if (comparefunc(curptr, item) == 0) {
			result = curptr;
			break;
		}
	}
	return result;
}

int32_t seqlist_get_index(SeqList * seqlist, void * item, equalto comparefunc)
{
	int32_t result = -1;
	for (uint16_t i = 0; i < seqlist->length; i++) {
		void* curptr = (uint8_t*)seqlist->data + seqlist->stepsize*i;
		if (comparefunc == NULL) {
			if (memcmp(curptr, item, seqlist->stepsize) == 0) {
				result = i;
				break;
			}
		}
		else if (comparefunc(curptr, item) == 0) {
			result = i;
			break;
		}
	}
	return result;
}
