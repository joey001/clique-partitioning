/*
 * cpp_ls.cpp
 * perturbation with
 *
 *  Created on: Jan 6, 2014
 *      Author: zhou
 */
#include <iostream>
#include <sstream>
#include <fstream>
#include <list>
#include <vector>
#include <algorithm>
#include <set>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>
#include <unistd.h>
#include "struct.h"

using namespace std;

#define DEBUG 0
#define DETAIL 1
#define LOCAL_OPT_RECORD 1
#define MAX_LOCAL_REC_NUM 50

#define MIN(a,b) ((a)<(b)?(a):(b))
#define PARENT(idx) ((idx-1)/2)
#define LEFT(idx) (2*idx+1)
#define RIGHT(idx) (2*idx+2)
#define EMPTY_IDX 0
const int MAX_CHAR = 10000;
const int MAX_VAL = 999999999;
const float CONST_E = 2.71828f;

//Modified in the script
char param_filename[1000] = "/home/zhou/cpp/benchmarks/charon_busco/rand400-5.txt";
int param_knownbest = 12133;
int param_runcnt = 1;		// the running time for each instance
int param_time = 400;		//the max time for tabu search procedure, unit: second

float param_alpha = 0.1f; 	// the ratio for weak perturbation
float param_beta = 0.4f;	// the raito for strong perturbation
//int param_iter_tabu = 500;
int param_seed = 123;
int param_coftabu = 15;	// coefficient of tabu search
const unsigned int param_klen = 4;
int param_tabu_depth = 10000000;
FILE *frec = NULL;

vector<int*> analyse_pool;

//int param_max_queue_size = 10;
typedef struct STGammaRow{
	int *pos_heap;	//the position of each gamma value
	int *heap;	//the partition index are ordered as a heap according to the corresponding gamma value
	int *values; // the current row
	int size;
}GammaRow;

typedef struct STGammaData{
	int **gammatbl;
	GammaRow* rows;
	int *best_improve;
	int global_best;

	RandAcessList *bound_ral;
	int *isbound;
}GammaData;

typedef enum{DESCENT, TABU, DIRECTED} SearchType;

typedef struct ST_Stats{
	int best_obj_value;
	int best_part_num;
	SearchType found_searh;

	double best_found_time;
	double total_run_time;

	long long best_found_iter;
	long long total_run_itr;
#if(DETAIL)
	int *bestpat;
#endif
}Stats;

#if (LOCAL_OPT_RECORD)
typedef struct ST_LC_OPT{
	int *pv;
	int s;
//	bool operator < (const ST_LC_OPT &other){
//		return s > other.s;
//	}
}LC_opt;
LC_opt lc_recs[MAX_LOCAL_REC_NUM];
int lc_num = 0;
int lc_min = 0;
#endif


int nnode;	/*Number of vertices in graph*/
int **matrix; /*Adjacent matrix of graph*/
int max_eweight; /*The maxweight of all the edges*/

Stats finstats;	/*The recorder*/

/*-------------partition data, maintain the current solution--------------*/
int *ppos;			/*ppos[pid] is the position of "pid" in "pbkt"*/
int *pbkt;			/*Elements from pbkt[0] to pbkt[pbkt_size-1] include all the partition ids*/
int pbkt_size = 0;	/*The number of partitions*/
int *pcnt;	 		/*The number of vertices in each partition*/
int *pvertex;		/*The partition of each vertex, for example, v is in partition pvertex[v] */
GammaData *pgamma;	/*Heap organized gamma table for each vertex to each partition*/
int fcurrent = 0;
int wbound;
RandAcessList *wral;	/*Records all the move with more than -wbound gain*/

clock_t starttime;	/*The start time of */
int fbest = 0;
long long **tabutbl;
//long long tabu_itr = 0;
//long long gpass = 0; /*Count the number of phases*/
long long gitr = 0;
/*-------------------------------------------------------------------------*/

void loadCompletedGraph(char *filename){
	ifstream fin;
	max_eweight = -MAX_VAL;

	fin.open(filename);
	if (fin.fail()){
		cerr << "Can not open file " << filename << endl;
		exit(0);
	}
	fin >> nnode;
	if (fin.eof()){
		cerr << "Empty file" << filename << endl;
		exit(0);
	}
	matrix = new int*[nnode];
	for (int i = 0; i < nnode; i++)
		matrix[i] = new int[nnode];
	int ni = 0;
	int val;
	while (ni < nnode){
		for (int nj = ni; nj < nnode; nj++){
			fin >> val;
			matrix[ni][nj] = matrix[nj][ni] = -val;
			if (val > max_eweight){
				max_eweight = val;
			}
		}
		ni++;
	}
}
void showPartition(FILE *f, int *p){
	int *mark = new int[nnode+1];
	memset(mark, 0, sizeof(int) * (nnode+1));
	for (int i = 0; i < nnode; i++) mark[p[i]] = 1;
	int cnt = 1;
	for (int i = 0; i < nnode + 1; i++){
		if (mark[i] != 0)
			mark[i] = cnt++;
	}
	fprintf(f, "Solution in order:");
	for(int i = 0; i < nnode; i++){
		fprintf(f, "%d ",mark[p[i]]);
	}
	fprintf(f, "\n");
}

void printMove(const StepMove *move){
	if (move->mvertex == 1){
		printf("ITR %lld: %d->%d:%d \n",gitr, move->orderedVertices[0],move->orderedTarget[0], move->inc);
	}
	if (move->mvertex == 2){
		printf("ITR %lld: %d->%d,%d->%d:%d \n",gitr, move->orderedVertices[0],move->orderedTarget[0],
				move->orderedVertices[1],move->orderedTarget[1],
				move->inc);
	}
}

void printOderedSteps(const vector<StepMove> *moves, int n = 10){
	for (int i = 0; i < n; i++){
		printMove(&(moves->at(i)));
	}
	printf("\n");
}
/**
 * Calculate F by part
 */
int calculateSum(int *part){
	int sum = 0;
	for (int i = 0 ;i < nnode; i++)
		for (int j = i; j < nnode; j++)
			if (part[i] == part[j])
				sum += matrix[i][j];
	return sum;
}
/**
 * Assign each node to a singleton partition
 */
void allocatePartitionData(){
	ppos = new int[nnode+1];
	pbkt = new int[nnode+1];
	pbkt_size = 0;
	pcnt = new int[nnode+1];
	pvertex = new int[nnode];
	for (int i = 0; i < nnode+1; i++){
		ppos[i] = i;
		pbkt[i] = i;
	}
	pbkt_size = 1;
	memset(pcnt, 0, sizeof(int) * (nnode+1));
	memset(pvertex, -1, sizeof(int) * nnode);
}

void disposePartition(){
	delete[] ppos;
	delete[] pbkt;
	delete[] pcnt;
	delete[] pvertex;
}
void buildPartition(int *vpart){
	for (int i = 0; i < nnode; i++){
		int pid = vpart[i];
		if (pid == EMPTY_IDX || pid > nnode){
			printf("Unvalid parition for node %d\n",pid);
			exit(0);
		}
		if (ppos[pid] >= pbkt_size){ //
			//ASSERTATION
			if (pbkt_size == nnode+1){
				printf("The bucket is full, no new partition could be added\n");
				exit(0);
			}
			int end_pid = pbkt[pbkt_size];
			swapAry(pbkt, pbkt_size, ppos[pid]);
			swapAry(ppos, end_pid, pid);
			pbkt_size++;
			pcnt[pid]++;
		}else{
			pcnt[pid]++;
		}
	}
	memcpy(pvertex, vpart, sizeof(int)*nnode);
}
void updatePartition(int v, int target){
	int oldpartition = pvertex[v];
	assert(target !=EMPTY_IDX);
	if (oldpartition == target)
		printf("Meaningless move\n");
	pcnt[oldpartition]--;
	pcnt[target]++;
	//the old cluster of v is empty
	if (pcnt[oldpartition] == 0){
		int end_pid = pbkt[pbkt_size-1];
		swapAry(pbkt, pbkt_size-1, ppos[oldpartition]);
		swapAry(ppos, end_pid, oldpartition);
		pbkt_size--;
	}
	pvertex[v] = target;
}

int calculateDistance(int *p1, int *p2){
	int sum = 0;
	for (int i = 0; i < nnode; i++){
		for (int j = i+1; j < nnode; j++){
			if (p1[i] == p1[j]){
				if (p2[i] != p2[j])
					sum++;
			}else{
				if (p2[i] == p2[j])
					sum++;
			}
		}
	}
	return sum;
}


/*****************************Data Structure of Gamma Table*********************************************/
/**
 * The gamma table have already updated before invoking this function
 */
void adjustGammaRow(GammaRow *gr, int partid){
	int curidx = gr->pos_heap[partid];
	int *v = gr->values;
	assert(curidx < gr->size);
	//ascend
	int parent = PARENT(curidx);
	if (parent >= 0 && v[gr->heap[parent]] < v[gr->heap[curidx]]){
		do{
			swapAry(gr->pos_heap, gr->heap[parent], gr->heap[curidx]);
			swapAry(gr->heap, parent, curidx);
			curidx = parent;
			parent = PARENT(curidx);
		}while (parent >= 0 && v[gr->heap[parent]] < v[gr->heap[curidx]]);
	}else{ //descend
		int end = 0;
		int biggeridx = curidx;
		while (!end){
			int left = LEFT(curidx);
			int right = RIGHT(curidx);
			if (left < gr->size && v[gr->heap[left]] > v[gr->heap[biggeridx]]){
				biggeridx = left;
			}
			if (right < gr->size && v[gr->heap[right]] > v[gr->heap[biggeridx]]){
				biggeridx = right;
			}
			if (biggeridx != curidx){
				swapAry(gr->pos_heap, gr->heap[biggeridx], gr->heap[curidx]);
				swapAry(gr->heap, biggeridx, curidx);
				curidx = biggeridx;
			}else{
				end = 1;
			}
		}
	}
}

void buildGammaRowFromPartition(GammaRow *gr,int *values, int len){
	int curpart = 0;
	gr->values = values;
	gr->pos_heap = new int[len];
	gr->heap = new int[len];
	gr->size = 0;
	memset(gr->pos_heap, -1, sizeof(int) * len);
	memset(gr->heap, -1, sizeof(int) * len);

	while (gr->size < pbkt_size){
		curpart = pbkt[gr->size];
		gr->heap[gr->size] = curpart;
		gr->pos_heap[curpart] = gr->size;
		gr->size++;
		adjustGammaRow(gr, curpart);
//		printRow(gr);
	}
}
void dropGammaRow(GammaRow *gr, int partid){
	int curidx = gr->pos_heap[partid];

	swapAry(gr->pos_heap, partid, gr->heap[gr->size-1]);
	swapAry(gr->heap, curidx, gr->size-1);
	gr->pos_heap[partid] = -1;
	gr->size--;
	if (curidx != gr->size)
		adjustGammaRow(gr, gr->heap[curidx]);
}

void addGammaRow(GammaRow *gr, int partid){
	int curidx = gr->size;
	gr->heap[curidx] = partid;
	gr->pos_heap[partid] = curidx;
	gr->size++;
	adjustGammaRow(gr, partid);
}


GammaData* buildGammaData(){
	GammaData* gamma = new GammaData;
	gamma->gammatbl = new int*[nnode];
	gamma->best_improve = new int[nnode];
	gamma->rows = new GammaRow[nnode];

	gamma->isbound = new int[nnode];
	gamma->bound_ral = ral_init(nnode);
	memset(gamma->isbound, 0, sizeof(int) * nnode);

	memset(gamma->best_improve, -1, sizeof(int) * nnode);
	//for each node
	for(int i = 0; i < nnode; i++){
		gamma->gammatbl[i] = new int[nnode+1];
		memset(gamma->gammatbl[i], 0, sizeof(int) * (nnode+1));
		for (int j = 0; j < nnode; j++){
			gamma->gammatbl[i][pvertex[j]] += matrix[i][j];
		}
		buildGammaRowFromPartition(&(gamma->rows[i]), gamma->gammatbl[i], nnode+1);

		GammaRow *grow =  &(gamma->rows[i]);
		int top_part = grow->heap[0];
		int selft_part = pvertex[i];
	//ERROR: The nodes which does not meet the following condition will NOT be initialized
	//In order to compare, we did not remove this error
	//		if (gamma->gammatbl[i][best_part] > gamma->gammatbl[i][selft_part])
	//			gamma->best_improve[i] = best_part;
		if (selft_part != top_part)
			gamma->best_improve[i] = top_part;
		else{
			int left = grow->heap[1];
			gamma->best_improve[i] = left;
			if (grow->size > 2 ){
				if (grow->values[grow->heap[2]] > grow->values[left])
					gamma->best_improve[i] = grow->heap[2];
			}
		}
		int delta = gamma->gammatbl[i][gamma->best_improve[i]] - gamma->gammatbl[i][selft_part];
		if (delta >= -wbound && (!gamma->isbound[i])){
			ral_add(gamma->bound_ral, i);
			gamma->isbound[i] = 1;
		}else if (delta < -wbound && (gamma->isbound[i])){
			ral_delete(gamma->bound_ral, i);
			gamma->isbound[i] = 0;
		}
	}
	return gamma;
}

/**
 * Move vertex from src_part to dest_part, call this function before execute the move
 */
void updateGamma(GammaData *gamma, int vertex, int src_part, int dest_part){
	//move vertex i to partition
	assert(pvertex[vertex] == src_part);
	for (int i = 0; i < nnode; i++){
		GammaRow *grow = &(gamma->rows[i]);
		int *v = grow->values;
		//update the gamma table
		gamma->gammatbl[i][src_part] -= matrix[vertex][i];
		gamma->gammatbl[i][dest_part] += matrix[vertex][i];
		if (pcnt[src_part] == 1){
			dropGammaRow(grow, src_part);
		}else
			adjustGammaRow(grow, src_part);

//		printRow(&(gamma->rows[i]));
		if (pcnt[dest_part] == 0){
			addGammaRow(grow, dest_part);
		}else{
			adjustGammaRow(grow, dest_part);
		}

		//adjust the best
		int self_part = (i == vertex) ? dest_part : pvertex[i];
		int top_part = grow->heap[0];
		if (self_part != top_part)
			gamma->best_improve[i] = top_part;
		else{
			int left = grow->heap[1];
			gamma->best_improve[i] = left;
			if (grow->size > 2 ){
				if (v[grow->heap[2]] > v[left])
					gamma->best_improve[i] = grow->heap[2];
			}
		}
		int delta = gamma->gammatbl[i][gamma->best_improve[i]] - gamma->gammatbl[i][self_part];
		if (delta >= -wbound && (!gamma->isbound[i])){
			ral_add(gamma->bound_ral, i);
			gamma->isbound[i] = 1;
		}else if (delta < -wbound && (gamma->isbound[i])){
			ral_delete(gamma->bound_ral, i);
			gamma->isbound[i] = 0;
		}
	}
}
void disposeGammaRow(GammaRow *gr){
	delete[] gr->heap;
	delete[] gr->pos_heap;
}
void disposeGamma(GammaData *gamma){
	for (int i = 0; i < nnode; i++){
		delete[] gamma->gammatbl[i];
		disposeGammaRow(gamma->rows+i);
	}
	delete[] gamma->rows;
	delete[] gamma->gammatbl;
	delete[] gamma->best_improve;
}

/**************************Other useful funtion****************************/
/*if the two solution is the same, return 1, otherwise return 0*/
int compareSolution(int *sa, int *sb, int n){
	int *dict = new int[n];
	memset(dict, -1, sizeof(int) * n);
	for (int i = 0; i < n; i++){
		if (dict[sa[i]] == -1)
			dict[sa[i]] = sb[i];
		else
			if (sb[i] != dict[sa[i]])
				return 0;
	}
	delete[] dict;
	return 1;
}




int decideTarget(StepMove move){
	int target = move.orderedTarget[0];
	assert(move.mvertex == 1);
	assert(ppos[target] < pbkt_size);
	//Move v to a new cluster
	if (target == EMPTY_IDX){
		assert(pbkt_size > 0);
		target = pbkt[pbkt_size];
		pbkt_size++;
	}
	return target;
}
/**
 * update partition, partcnt, tabu,move freq
 */
void executeStepMove(StepMove &move){
	int target = decideTarget(move);
	int curvertex = move.orderedVertices[0];

	//update gamma table
	updateGamma(pgamma, curvertex, pvertex[curvertex], target);
	updatePartition(curvertex, target);
	fcurrent += move.inc;
}
void recordBest(SearchType srh){
	fbest = fcurrent;

	finstats.best_obj_value = fbest;
	finstats.best_part_num = pbkt_size - 1;

	finstats.best_found_time = (float)(clock() - starttime) / CLOCKS_PER_SEC;
	finstats.best_found_iter = gitr;
	finstats.found_searh = srh;
#if (DETAIL)
	memcpy(finstats.bestpat, pvertex, sizeof(int)*nnode);
#endif
}
void clearFinalStats(){
	finstats.best_found_time = 0.0f;
	finstats.total_run_time = 0.0f;

	finstats.best_found_iter = 0l;
	finstats.total_run_itr = 0l;

	finstats.best_part_num = -1;
	finstats.best_obj_value = -MAX_VAL;
#if(DETAIL)
	if (finstats.bestpat == NULL)
		finstats.bestpat = new int[nnode];
	memset(finstats.bestpat, 0, sizeof(int) * nnode);
#endif

}
void initDataStructure(){
	allocatePartitionData();
	clearFinalStats();
	int *initpart = new int[nnode];
	for (int i = 0; i < nnode; i++){
		initpart[i] = i+1;
	}
	buildPartition(initpart);
	fbest = fcurrent =0;
	tabutbl = new long long*[nnode];// for moving two nodes
	for (int i = 0; i < nnode; i++){
		tabutbl[i] = new long long[nnode+1];
		memset(tabutbl[i], -1, sizeof(long long) * (nnode+1));
	}
	pgamma = buildGammaData();

	starttime = clock();
	gitr = 0;
	wbound = max_eweight; /*set wbound as max_eweight*/
#if (LOCAL_OPT_RECORD)
	for (int i = 0; i < MAX_LOCAL_REC_NUM; i++){
		lc_recs[i].pv = NULL;
		lc_recs[i].s = 0;
	}
	lc_num = 0;
	lc_min = MAX_VAL;
#endif
}

/**
 * Dispose gamma, tabutbl, partitiondata
 */
void clearDataStructure(){
	disposeGamma(pgamma);
	disposePartition();
	for (int i = 0; i < nnode; i++){
			delete[] tabutbl[i];
	}
	delete[] tabutbl;
}
void descentSearch(){
	int improved = 1;

	while (improved == 1){
		improved = 0;
		for (int i = 0; i < pgamma->bound_ral->vnum; i++){
			int currentnode = pgamma->bound_ral->vlist[i];
			int best_part = pgamma->best_improve[currentnode];
			int gain = pgamma->gammatbl[currentnode][best_part] - pgamma->gammatbl[currentnode][pvertex[currentnode]];
			if (gain > 0){
				improved = 1;
				StepMove sm = make_1StepMove(currentnode, best_part, gain);
				executeStepMove(sm);
			}
		}
		gitr++;
	}
	//assert(calculateSum(pvertex) == fcurrent);
	if (fcurrent > fbest){
		recordBest(DESCENT);
	}
//	printf("After descent: current %d-%d \n",gitr, fcurrent);
}
int dbg_cnt=0;
int findBestMove(StepMove *chosenMove){
	StepMove maxmove;
	int bestinc = -MAX_VAL;
	int eqcnt = 0;
	dbg_cnt=0;
	//for (int nodeidx = 0; nodeidx < nnode; nodeidx++){
//	printf("bound len %d\n",pgamma->bound_ral->vnum);
	for (int i = 0; i < pgamma->bound_ral->vnum; i++){
		int nodeidx = pgamma->bound_ral->vlist[i];
		int best_part = pgamma->best_improve[nodeidx];
		assert(best_part != pvertex[nodeidx]);
		if (best_part == EMPTY_IDX && pcnt[pvertex[nodeidx]] == 1)
			continue; //meaningless move
		int gain = pgamma->gammatbl[nodeidx][best_part] - pgamma->gammatbl[nodeidx][pvertex[nodeidx]];
//		printf("node %d(%d): to part %d %d, gain %d\n ",nodeidx,pgamma->gammatbl[nodeidx][pvertex[nodeidx]],
//		                                                  best_part, pgamma->gammatbl[nodeidx][best_part],gain);
		StepMove curmove = make_1StepMove(nodeidx, best_part, gain);
		if (fcurrent + gain > fbest || tabutbl[nodeidx][best_part] < gitr){
			if (gain > bestinc){
				bestinc = gain;
				maxmove = curmove;
				eqcnt = 0;
			}else if (gain == bestinc){ // Find the best one from all the moves with negative improvement;
				eqcnt++;
				if (rand() % eqcnt == 0)
					maxmove = curmove;
			}
		}
		if (gain >= -wbound )
			dbg_cnt++;
	}
//	printf("count %d\n",dbg_cnt);
	if (bestinc == -MAX_VAL)
		return 0;
	memcpy(chosenMove, &maxmove, sizeof(StepMove));
	return 1;
}

void tabuExploreSearch(){
	int non_improve_cnt = 0;
#if (LOCAL_OPT_RECORD)
	int local_opt = 0;
	int *plocal_s = new int[nnode];
#endif
	int cycle_cnt = 0;
//	int hit_cnt = 0;
	while (non_improve_cnt < nnode ){
		StepMove chosenMove;
		int find = findBestMove(&chosenMove);
//		printf("Delta: %d\n",chosenMove.inc);
//		if (chosenMove.inc >= -wbound && chosenMove.inc <=wbound ){
//			hit_cnt++;
//		}
		if (find == 0){
//			printf("Iter %d None node could be moved\n",gitr);
			//exit(0);
			break;
		}
		int curvertex = chosenMove.orderedVertices[0];
		int oldpart = pvertex[curvertex];
		/*Update tabu table*/
		if (pcnt[oldpart] == 1){
			//tabutbl[curvertex][EMPTY_IDX] = gitr + MIN(param_coftabu, nnode / 4) + rand() % pbkt_size;
			tabutbl[curvertex][EMPTY_IDX] = gitr + 7 + rand() % pbkt_size;
		// if there are only two nodes in a cluster and an edge will be moved outside, the cluster will be clean out
		//else if (pcnt[oldpart] == 2 && move.mvertex == 2)
		//	tabutbl[curvertex][EMPTY_IDX] = tabutenure_base + rand() % pbkt_size;
		}else{
//			tabutbl[curvertex][oldpart] = gitr + MIN(param_coftabu, nnode / 4) + rand() % pbkt_size;
			tabutbl[curvertex][oldpart] = gitr + 7 + rand() % pbkt_size;
		}
		executeStepMove(chosenMove);
		if (fcurrent > fbest){
			recordBest(TABU);
			non_improve_cnt = 0;
			if (fbest >= param_knownbest){
				break;
			}
		}else
			non_improve_cnt++;
#if LOCAL_OPT_RECORD
		if (fcurrent > local_opt){
			local_opt = fcurrent;
			memcpy(plocal_s, pvertex, sizeof(int)*nnode);
		}
#endif
		gitr++;
		cycle_cnt++;
	}
//	if (cycle_cnt > 0){
//		printf("Iter %d: Hit probability %.2f\n",cycle_cnt, hit_cnt/(float)cycle_cnt);
//	}

//	tabu_itr = tabu_itr + MIN(param_coftabu, nnode / 4) + nnode; //jump D step to invalid the data in tabu table
//	printf("After tabu: current %d-%d \n",gitr, fcurrent);
	/*record local optimum*/
#if LOCAL_OPT_RECORD
	if (lc_num < MAX_LOCAL_REC_NUM){
		lc_recs[lc_num].pv = plocal_s;
		lc_recs[lc_num].s = local_opt;
		lc_num++;
	}
#endif
}
/**
 */
void directedPerturb(){
	int *moved = new int[nnode];
	memset(moved, 0, sizeof(int) * nnode);
//	int param_max_queue_size =  10 + random() % 10;
	int param_max_queue_size = 10 + rand() % (int)(0.05* nnode + 1);
#if (DEBUG)
	printf("PERTURBATION starts from %lldth\n",tabu_itr);
#endif
//	strength = param_alpha  * nnode + rand() % (int((param_beta - param_alpha) * nnode));
	//Ensure the strenth is less than the number of nodes;
	int strength = (int)(param_alpha * nnode) + rand() % (int)(param_beta* nnode);
	int itr_cnt = 0;
	while(itr_cnt < strength){
		/*A heap structure queue to keep the the param_max_queue_size maximum*/
		FixedSizeQueue* fq = newFixedSizeQueue(param_max_queue_size);
		for (int i = 0; i < nnode; i++){
			int best_part = pgamma->best_improve[i];
			if (moved[i] == 1 && (best_part != 0 || (best_part == 0 && pcnt[pvertex[i]]>1) ))
				continue;
			int gain = pgamma->gammatbl[i][best_part] - pgamma->gammatbl[i][pvertex[i]];
			StepMove sm = make_1StepMove(i, best_part, gain);
			insertFixedSizeQueue(fq, sm);
		}
		// chose a random move in the candidate queue
		StepMove choice = randSelect(fq);
		//execute the move
		executeStepMove(choice);
		if (fcurrent > fbest){
			recordBest(DIRECTED);
		}
#if (DEBUG)
		printMove(&choice);
#endif
		moved[choice.orderedVertices[0]] = 1;
		itr_cnt++;
		gitr++;
	}
#if (DEBUG)
	printf("PERTURBATION ends at %lldth\n",gpass);
#endif
	delete[] moved;
//	printf("After directed: current %d-%d \n",gitr, fcurrent);
}


void threePhaseMain(int index){
//	int perturb_itr = 0;
//	int end = 0;
//	//DEBUG time
//	struct timeval tabustart,tabuend;
//	struct timeval pstart,pend;
//	long long psum=0,tabusum=0;
	gitr = 0;
	int len = 0;
	initDataStructure();
	//gitr < 10000000
#if (LOCAL_OPT_RECORD)
	while ((double)(clock() - starttime) / CLOCKS_PER_SEC < (double)param_time || lc_num < MAX_LOCAL_REC_NUM){
#else
	while ((double)(clock() - starttime) / CLOCKS_PER_SEC < (double)param_time){
#endif
		descentSearch();
		tabuExploreSearch();
		if (fbest >= param_knownbest){
			break;
		}
		directedPerturb();
	}
	finstats.total_run_time = (double)(clock() - starttime) / CLOCKS_PER_SEC;
	finstats.total_run_itr = gitr;
	clearDataStructure();
}



int countPartition(int *p){
	int npart = 0;
	//check
	int *cnt = new int[nnode];
	memset(cnt, 0, sizeof(int) * nnode);
	for (int i = 0; i < nnode; i++)
		cnt[p[i]]++;
	for (int i = 0; i < nnode; i++)
		if (cnt[i] != 0)
			npart++;
	return npart;
}

void showUsage(){
	cerr << "usage: -f <filepath> [-s <cut off time>] [-a <alpha>] [-b <beta>] [-t <tabu tenure base>] [-r <run time>]\n" << endl;
}
void readParameters(int argc, char **argv){
	for (int i = 1; i < argc; i+=2){
		if (argv[i][0] != '-' || argv[i][2] != 0){
			showUsage();
			exit(0);
		}else if (argv[i][1] == 'f'){	/*The file name*/
			strncpy(param_filename, argv[i+1],1000);
		}else if (argv[i][1] == 's'){	/*The maximum time*/
			param_time = atoi(argv[i+1]);
		}else if(argv[i][1] == 'a'){
			param_alpha = atof(argv[i+1]);
		}else if (argv[i][1] == 'b'){
			param_beta =  atof(argv[i+1]);
		}else if (argv[i][1] == 't'){
			param_coftabu = atoi(argv[i+1]);
		}else if (argv[i][1] == 'r'){
			param_runcnt = atoi(argv[i+1]);
		}else if(argv[i][1] == 'v'){
			param_knownbest = atoi(argv[i+1]);
		}else if(argv[i][1] == 'g'){
			param_seed = atoi(argv[i+1]);
		}
	}
	/*check parameters*/
	if (strlen(param_filename) == 0){
		cerr << "No input data" << endl;
		exit(1);
	}
}


int setupRecordFile(){
	/*creat file in current direct*/
	char path_cwd[FILENAME_MAX];
	char *graph_name = basename(param_filename);
	char file_name[FILENAME_MAX];
	int randnum = rand();

	getcwd(path_cwd, FILENAME_MAX);
	sprintf(file_name, "%s/%s_%5d.rec",path_cwd, graph_name, randnum);

//	sprintf(file_name, "%s/%s_best",path_cwd, graph_name);
	frec = fopen(file_name, "a+");
	if (frec == NULL){
		return 0;
	}
	return 1;
}


int main(int argc, char **argv){
	/*display and read parameters*/
	readParameters(argc, argv);
	srand(param_seed);
	//set logging device
	if (0 == setupRecordFile()){
		cerr << "Failed in open log file" << endl;
		exit(2);
	}
	//frec = stdout;

	//load graph by jostle format
	loadCompletedGraph(param_filename);

	/*initialize some random parameter*/
//	srand((unsigned int)time(NULL));
	for (int i = 0; i < argc; i++){
		fprintf(frec, "%s ", argv[i]);
	}
	fprintf(frec, "\n");

	fprintf(frec, "idx\t best_v\t npat\t find_t\t find_i\t ttl_t\t  ttl_i\n");
	int cnt = 0;
	float sumtime = 0;
	int sumres = 0;
	long long sumiter = 0;

	int bestInAll = -MAX_VAL;
	int *bestInAlllPartition = new int[nnode];
	while (cnt < param_runcnt){
		threePhaseMain(cnt);
#if (DETAIL)
		int checked = calculateSum(finstats.bestpat);
		if (checked != finstats.best_obj_value){
			fprintf(stderr, "fatal error, result does not match %d %d", finstats.best_obj_value, checked);
			exit(1);
		}
#endif
		fprintf(frec, "%d\t %d\t %d\t %.2f\t %lld\t %.2f\t %lld\n",cnt+1, finstats.best_obj_value, finstats.best_part_num,
				finstats.best_found_time, finstats.best_found_iter, finstats.total_run_time,finstats.total_run_itr);
		if (finstats.best_obj_value > bestInAll){
			bestInAll = finstats.best_obj_value;
#if (DETAIL)
			memcpy(bestInAlllPartition, finstats.bestpat, sizeof(int) * nnode);
#endif
		}
		cnt++;
		sumtime += finstats.best_found_time;
		sumres += finstats.best_obj_value;
		sumiter += finstats.best_found_iter;
#if (LOCAL_OPT_RECORD)
		/*clear local optimum*/
		fprintf(frec, "local optimums number %d. gloabl best %d, Format(dis, fit)\n",lc_num, bestInAll);
		for (int i = 0; i < lc_num; i++){
			int di = calculateDistance(lc_recs[i].pv, finstats.bestpat);
			double nor_di = (double)di *2/ (nnode *(nnode-1));
			fprintf(frec, "%d, %.4f, %d\n", di, nor_di, lc_recs[i].s);
			delete lc_recs[i].pv;
		}
		lc_num = 0;
#endif
	}
	fprintf(frec,"best result: %d\n",bestInAll);
	fprintf(frec, "average time:%.2f\n", sumtime / param_runcnt);
	fprintf(frec, "average result:%.2f\n", (float)sumres/param_runcnt);
	fprintf(frec, "average best iteration: %d\n", sumiter/param_runcnt);
#if (DETAIL)
	showPartition(frec, bestInAlllPartition);
#endif
	fclose(frec);
	return 0;
}

