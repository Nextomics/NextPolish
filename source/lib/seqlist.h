#ifndef SEQLIST_H
#define SEQLIST_H

#include<stdlib.h>
#include<stdint.h>
#include<stdio.h>
#include<malloc.h>
#include<string.h>

typedef struct {
	void *data;
	uint32_t maxsize;
	uint32_t length;
	uint32_t stepsize;
} SeqList;

SeqList* seqlist_init(uint32_t stepsize, uint32_t max);

void seqlist_destory(SeqList* seqlist);

void seqlist_append(SeqList* seqlist, void* item);

void seqlist_remove(SeqList* seqlist, uint32_t pos);

void* seqlist_index(SeqList* seqlist, uint32_t pos);

void seqlist_setitem(SeqList* seqlist, void* item, const uint32_t pos);
void seqlist_update(SeqList* seqlist, void* dest, void* item);

typedef int32_t(*equalto)(void* first, void* second);
void* seqlist_find(SeqList* seqlist, void* item, equalto comparefunc);
int32_t seqlist_get_index(SeqList* seqlist, void* item, equalto comparefunc);

#endif // !SEQLIST_H

