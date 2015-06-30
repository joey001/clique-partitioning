/*
 * fsqueue.cpp
 *
 *  Created on: Apr 9, 2014
 *      Author: zhou
 */
#include <vector>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "struct.h"

using namespace std;

void generateRandList(int *randlist, int len){
	int idx = 0;
	assert(randlist != NULL);

	for (idx = 0; idx < len; idx++){
		randlist[idx] = idx;
	}
	for (idx = 0; idx < len; idx++){
		int randid = rand() % len;
		int tmp = randlist[idx];
		randlist[idx] = randlist[randid];
		randlist[randid] = tmp;
	}
}
void swapAry(int *ary, int idx1, int idx2){
	int tmp = ary[idx1];
	ary[idx1] = ary[idx2];
	ary[idx2] = tmp;
}
/*******************定义StepMove*********************************/
StepMove make_1StepMove(int vertex, int target, int gain){
	StepMove sm;
	sm.mvertex = 1;
	sm.orderedVertices[0] = vertex;
	sm.orderedTarget[0] = target;
	sm.orderedTarget[1] = target;
	sm.inc = gain;
	return sm;
}
StepMove make_2StepMove(int vertex1, int vertex2, int target1, int target2, int gain){
	StepMove sm;
	sm.mvertex = 2;
	sm.orderedVertices[0] = vertex1;
	sm.orderedVertices[1] = vertex2;
	sm.orderedTarget[0] = target1;
	sm.orderedTarget[1] = target2;
	sm.inc = gain;
	return sm;
}

/************************利用堆实现的一个最大K队列***********************************/
static StepMove *m_pool;

FixedSizeQueue* newFixedSizeQueue(const int size){
	FixedSizeQueue *fq = new FixedSizeQueue;
	if (m_pool == NULL){
		m_pool = new StepMove[301];
	}
	assert(size > 0);
	fq->members = m_pool;
	fq->fixedsize = size;
	fq->cursize = 0;
	return fq;
}

void insertFixedSizeQueue(FixedSizeQueue *fq,StepMove &sm){
	int val = sm.inc;
	if (fq->cursize < fq->fixedsize){
		int idx = fq->cursize;
		for (; idx > 0 && (fq->members[(idx-1) / 2]).inc > val; idx = (idx - 1) / 2)
			fq->members[idx] = fq->members[(idx-1)/2];
		fq->members[idx] = sm;
		fq->cursize++;
	}else{
		if (val > fq->members[0].inc){
			fq->members[0] = sm;
			int curidx = 0;
			int smalleridx = curidx;
			int end = 0;
			while(!end){
				int left = 2 * curidx + 1;
				int right = 2 * curidx + 2;
				if (left < fq->cursize && (fq->members[left]).inc < (fq->members[smalleridx]).inc){
					smalleridx = left;
				}
				if (right < fq->cursize && (fq->members[right]).inc < (fq->members[smalleridx]).inc){
					smalleridx = right;
				}
				if (smalleridx != curidx){
					//swap the child node with parent node
					StepMove tmp = fq->members[curidx];
					fq->members[curidx] = fq->members[smalleridx];
					fq->members[smalleridx] = tmp;
					curidx = smalleridx;
				}else
					end = 1;
			}
		}
	}
}

StepMove randSelect(FixedSizeQueue *fq){
	int randidx = rand() % fq->cursize;
		return fq->members[randidx];
}
void disposeQueue(FixedSizeQueue *fq){
	delete fq;
}
/****************************定义random access list结构***************************/
RandAcessList* ral_init(int capacity){
	RandAcessList *ral = new RandAcessList;
	ral->vlist = new int[capacity];
	ral->vpos = new int[capacity];
	ral->vnum = 0;
	ral->capacity = capacity;
	for (int i = 0; i< capacity;i++){
		ral->vpos[i] = capacity;
	}
	return ral;
}

void ral_add(RandAcessList *ral, int vid){
//	assert(ral->vpos[vid] >= ral->vnum);
	ral->vlist[ral->vnum] = vid;
	ral->vpos[vid] = ral->vnum;
	ral->vnum++;
}
void ral_delete(RandAcessList *ral, int vid){
//	assert(ral->vpos[vid] < ral->vnum);
	int last_id = ral->vlist[ral->vnum - 1];
	int id_pos = ral->vpos[vid];
	ral->vlist[id_pos] = last_id;
	ral->vpos[last_id] = id_pos;
	ral->vnum--;
//	ral->vpos[vid] = ral->vnum; /*It is not obligatory*/
}

void ral_clear(RandAcessList *ral){
	ral->vnum = 0;
}
void ral_release(RandAcessList *ral){
	delete[] ral->vlist;
	delete[] ral->vpos;
	delete ral;
}
void ral_showList(RandAcessList *ral, FILE *f){
	fprintf(f, "%d: ",ral->vnum);
	for(int i = 0;i < ral->vnum; i++){
		fprintf(f, "%d ", ral->vlist[i]);
	}
	fprintf(f, "\n");
}
