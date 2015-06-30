/*
 * struct.h
 *
 *  Created on: Apr 9, 2014
 *      Author: zhou
 */

#ifndef STRUCT_H_
#define STRUCT_H_
#include <vector>
#include <stdio.h>
using namespace std;
typedef struct{
	int mvertex;
	int orderedVertices[1];
	int orderedTarget[1];
	int inc;
}StepMove;


extern StepMove make_1StepMove(int vertex, int target, int gain);
extern StepMove make_2StepMove(int vertex1, int vertex2, int target1, int target2, int gain);
extern void printMove(const StepMove *move);
extern void printOderedSteps(const vector<StepMove> *moves, int n);

typedef struct ST_FixedSizeQueue{
	StepMove *members;
	int cursize;
	int fixedsize;
}FixedSizeQueue;

extern void generateRandList(int *randlist, int len);
extern void swapAry(int *ary, int idx1, int idx2);


extern FixedSizeQueue* newFixedSizeQueue(const int size);
extern void insertFixedSizeQueue(FixedSizeQueue *fq,StepMove &sm);
extern StepMove randSelect(FixedSizeQueue *fq);
extern void disposeQueue(FixedSizeQueue *fq);

typedef struct ST_RandAccessList{
	int *vlist;
	int *vpos;
	int vnum;
	int capacity;
}RandAcessList;
extern RandAcessList* ral_init(int capacity);
extern void ral_add(RandAcessList *ral, int vid);
extern void ral_delete(RandAcessList *ral, int vid);
extern void ral_clear(RandAcessList *ral);
extern void ral_release(RandAcessList *ral);
extern void ral_showList(RandAcessList *ral, FILE *f);
#endif /* STRUCT_H_ */
