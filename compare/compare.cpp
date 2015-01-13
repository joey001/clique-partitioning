/*
 * compare.cpp
 *单纯比较各种领域选择方式效果如何
 *  Created on: Apr 16, 2014
 *      Author: zhou
 */
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "struct.h"
using namespace std;

#define MAX_VAL 999999
#define CHOICE_BEST 0
#define CHOICE_RAND 1
#define CHOICE_WORST 2
#define CHOICE_MIXED 3
int nnode;
int **matrix;
char param_graph_file[1000] = "/home/zhou/workspace/instances/cpp_instances/rand200-5.txt";
int param_run_cnt= 30;
int param_choice = CHOICE_MIXED;
unsigned int param_klen =8;

int best_obj;
int cur_obj;


#define EMPTY_IDX 0
int *ppos;
int *pbkt;
int pbkt_size = 0;
int *pcnt;
int *pvertex;

int **gamma_tbl;
/*-------------------partition-------------------------*/
/**
 * Assign each node to a singleton partition
 */

void disposePartition(){
	delete[] ppos;
	delete[] pbkt;
	delete[] pcnt;
	delete[] pvertex;
}
void swapAry(int *ary, int idx1, int idx2){
	int tmp = ary[idx1];
	ary[idx1] = ary[idx2];
	ary[idx2] = tmp;
}
void buildPartition(int *vpart){
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
 * Load the completed graph from file
 */
void loadCompletedGraph(char *filename){
	ifstream fin;
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
		}
		ni++;
	}
}

/*****************************Operation of Gamma Table*********************************************/
void buildGamma(){
	int i = 0;
	gamma_tbl = new int*[nnode];

	for (i = 0; i < nnode; i++){
		gamma_tbl[i] = new int[nnode+1];
		memset(gamma_tbl[i], 0, sizeof(int)*(nnode+1));
	}
	for (int i = 0; i < nnode-1; i++)
		for (int j = i+1; j < nnode; j++){
			gamma_tbl[i][pvertex[j]] += matrix[i][j];
			gamma_tbl[j][pvertex[i]] += matrix[i][j];
		}
}

/**
 * move node t to partition k
 */
void updateGamma(int t, int k){
	for (int i = 0; i < nnode; i++){
		gamma_tbl[i][pvertex[t]] -= matrix[t][i];
		gamma_tbl[i][k] += matrix[t][i];
	}
}
/****************************End of operation of Gamma Table*********************************************/

/**
 *Generate a list of integers with length len
 */
void generateRandList(int *randlist, int len){
	int idx = 0;
	assert(randlist != NULL);

	for (idx = 0; idx < len; idx++){
		randlist[idx] = idx;
	}
	for (idx = 0; idx < len; idx++){
		int randid = random() % len;
		int tmp = randlist[idx];
		randlist[idx] = randlist[randid];
		randlist[randid] = tmp;
	}
}


int decideTarget(StepMove move){
	int target = move.orderedTarget[0];
	assert(move.orderedTarget[0] == move.orderedTarget[1]);
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
	for (int i = 0; i < move.mvertex; i++){
		int curvertex = move.orderedVetexes[i];
		int oldpart = pvertex[curvertex];

		//update gamma table
		updateGamma(curvertex, target);
		updatePartition(curvertex, target);
	}
	cur_obj += move.inc;
	//check
//	for (int i = 0; i < nnode ;i++){
//		int p = pvertex[i];
//		if (ppos[p] >= pbkt_size){
//			printf("Miss parition %d\n",p);
//		}
//	}

}

int preEstimatePartCnt(){
	int *randlst = new int[nnode];
	int *initpart = new int[nnode];
	for (int i = 0; i < nnode; i++){
		initpart[i] = i+1;
	}
	int **gamma_tmp = new int*[nnode];
	for (int i = 0; i < nnode; i++){
		gamma_tmp[i] = new int[nnode+1];
		memset(gamma_tmp[i], 0, sizeof(int)*(nnode+1));
	}
	for (int i = 0; i < nnode-1; i++)
		for (int j = i+1; j < nnode; j++){
			gamma_tmp[i][initpart[j]] += matrix[i][j];
			gamma_tmp[j][initpart[i]] += matrix[i][j];
		}
	generateRandList(randlst, nnode);
	int improved = 1;
	while(improved){
		for (int i = 0; i < nnode; i++){
			int maxgain = -99999;
			int destpat = -1;
			int currentnode = randlst[i];
			for (int othernode = 0; othernode< nnode; othernode++){
				int gain = gamma_tmp[currentnode][initpart[othernode]] - gamma_tmp[currentnode][initpart[currentnode]];
				if (initpart[currentnode] != initpart[othernode] && gain > maxgain){
					maxgain = gain;
					destpat = initpart[othernode];
				}
			}
			if (maxgain > 0){ //Move currentnode to  destpat
				for (int i = 0; i < nnode; i++){
					gamma_tmp[i][initpart[currentnode]] -= matrix[currentnode][i];
					gamma_tmp[i][destpat] += matrix[currentnode][i];
				}
				initpart[currentnode] = destpat;
				improved = 1;
			}else
				improved = 0;
		}
	}
//	for (int i = 0; i < nnode; i++){
//		delete[] gamma_tmp[i];
//		printf("%d ", initpart[i]);
//	}
//	printf("\n");
	delete[] gamma_tmp;
	delete[] randlst;
	int *cnt_pat = new int[nnode+1];
	int cnt= 0;
	memset(cnt_pat, 0, sizeof(int) * (nnode+1));
	for (int i = 0; i < nnode; i++)
		cnt_pat[initpart[i]]++;
	for (int i = 1; i < nnode+1; i++){
		if (cnt_pat[i] > 0) cnt++;
	}
	delete[] initpart;
	delete[] cnt_pat;
	return cnt;
}
int generateInitSolution(){
	int number = preEstimatePartCnt();
	assert(number > 1);
	int nodes_per_part = nnode/(number-1);
	int *partition = new int[nnode];
	int *rlist = new int[nnode];
	int idx =0;
	int cid = 1;
	int cnt = 0;
	generateRandList(rlist, nnode);
	while (idx < nnode){
		if (cnt == nodes_per_part){
			cnt = 0;
			cid++;
		}
		partition[rlist[idx]] = cid;
		cnt++;
		idx++;
	}
	buildPartition(partition);
	buildGamma();
	delete[] partition;
	delete[] rlist;
	return calculateSum(pvertex);
}

void showUsage(){
	printf("Wrong format\n");
}

void bestSelectionDescent(int *iters){
	int local_opt = 0;
	*iters = 0;
	while (!local_opt){
		StepMove maxmove;
		int maxinc = -MAX_VAL;
		int eqcnt = 0;
		for (int nodeidx = 0; nodeidx < nnode; nodeidx++){
			for (int k = 0; k < pbkt_size; k++){
				int curp = pbkt[k];
				if (curp == pvertex[nodeidx]) continue;
				if (curp == EMPTY_IDX && pcnt[pvertex[nodeidx]] == 1) continue;
				int gain = gamma_tbl[nodeidx][curp] - gamma_tbl[nodeidx][pvertex[nodeidx]];
				StepMove curmove = make_1StepMove(nodeidx, curp, gain);
				if (gain > maxinc){
					maxinc = gain;
					maxmove = curmove;
					eqcnt = 0;
				}else if (gain == maxinc){
					eqcnt++;
					if (rand() % eqcnt == 0)
						maxmove = curmove;
				}
			}
		}
		if (maxinc > 0){
			executeStepMove(maxmove);
			(*iters)++;
		}
		else
			local_opt = 1;
	}
}

void randomSelectionDescent(int *iters){
	int local_opt = 0;
	*iters = 0;
	int *randlist = new int[nnode];
	while (!local_opt){
		StepMove maxmove;
		int maxgain = 0;
		generateRandList(randlist, nnode);
		for (int i = 0; i < nnode; i++){
			int nodeidx =randlist[i];
			for (int k = 0; k < pbkt_size; k++){
				int curp = pbkt[k];
				if (curp == pvertex[nodeidx]) continue;
				if (curp == EMPTY_IDX && pcnt[pvertex[nodeidx]] == 1) continue;
				int gain = gamma_tbl[nodeidx][curp] - gamma_tbl[nodeidx][pvertex[nodeidx]];
				StepMove curmove = make_1StepMove(nodeidx, curp, gain);
				if (gain > 0){
					maxgain = gain;
					maxmove = curmove;
					goto findpos;
				}
			}
		}
		findpos:
		if (maxgain > 0){
			executeStepMove(maxmove);
			(*iters)++;
		}else
			local_opt = 1;
	}
	delete[] randlist;
}
void worstSelectionDescent(int *iters){
	int local_opt = 0;
	*iters = 0;
	while (!local_opt){
		StepMove maxmove;
		int mininc = MAX_VAL;
		int eqcnt = 0;

		for (int nodeidx = 0; nodeidx < nnode; nodeidx++){
			for (int k = 0; k < pbkt_size; k++){
				int curp = pbkt[k];
				if (curp == pvertex[nodeidx]) continue;
				if (curp == EMPTY_IDX && pcnt[pvertex[nodeidx]] == 1) continue;
				int gain = gamma_tbl[nodeidx][curp] - gamma_tbl[nodeidx][pvertex[nodeidx]];
				StepMove curmove = make_1StepMove(nodeidx, curp, gain);
				if (gain > 0 && gain < mininc){
					mininc = gain;
					maxmove = curmove;
					eqcnt = 0;
				}else if (gain == mininc){
					eqcnt++;
					if (rand() % eqcnt == 0)
						maxmove = curmove;
				}
			}
		}
		if (mininc < MAX_VAL){
			executeStepMove(maxmove);
			(*iters)++;
		}
		else
			local_opt = 1;
	}
}
void mixedSelectionDescent(int *iters){
	int local_opt = 0;
	*iters = 0;
	int *randlist = new int[nnode];
	while (!local_opt){
		StepMove maxmove;
		vector<StepMove> kqueue;
		generateRandList(randlist, nnode);
		for (int i = 0; i < nnode; i++){
			int nodeidx =randlist[i];
			for (int k = 0; k < pbkt_size; k++){
				int curp = pbkt[k];
				if (curp == pvertex[nodeidx]) continue;
				if (curp == EMPTY_IDX && pcnt[pvertex[nodeidx]] == 1) continue;
				int gain = gamma_tbl[nodeidx][curp] - gamma_tbl[nodeidx][pvertex[nodeidx]];
				StepMove curmove = make_1StepMove(nodeidx, curp, gain);
				if (gain > 0 && kqueue.size() < param_klen){
					kqueue.push_back(curmove);
					if (kqueue.size() == param_klen)
						goto findpos;
				}
			}
		}
		findpos:
		if (kqueue.size() > 0){
			int mingain = MAX_VAL;
			StepMove sp;
			for (int i = 0; i < kqueue.size(); i++)
				if (kqueue[i].inc < mingain){
					mingain = kqueue[i].inc;
					sp = kqueue[i];
				}
			executeStepMove(sp);
			(*iters)++;
		}else
			local_opt = 1;
	}
	delete[] randlist;
}
typedef struct ST_Stats{
	int run_cnt;
	char method[200];

	int *objval;
	int best_objval;
	float ave_objval;

	float *runtime;
	float best_runtime;
	float ave_runtime;

	int *iterations;
	int best_iteration;
	float ave_iteration;
}RunStats;

void runCompare(int choice, int run_cnt, RunStats *p_RunSts){
	p_RunSts->run_cnt = run_cnt;
	p_RunSts->objval = new int[run_cnt];
	p_RunSts->runtime = new float[run_cnt];
	p_RunSts->iterations = new int[run_cnt];
	float sum_objval = 0;
	float sum_runtime = 0;
	int sum_iters = 0;
	int max_obj = -MAX_VAL;
	float min_time = MAX_VAL;
	int min_iter = MAX_VAL;
	for (int i = 0; i < run_cnt; i++){
		best_obj = cur_obj = generateInitSolution();
		clock_t start_time = clock();
		int iters = 0;
		switch (choice){
			case CHOICE_BEST:
				bestSelectionDescent(&iters);
				strcpy(p_RunSts->method, "Select best improvement");
				break;
			case CHOICE_RAND:
				randomSelectionDescent(&iters);
				strcpy(p_RunSts->method, "Select the first improvement");
				break;
			case CHOICE_WORST:
				worstSelectionDescent(&iters);
				strcpy(p_RunSts->method, "Select the worst improvement");
				break;
			case CHOICE_MIXED:
				mixedSelectionDescent(&iters);
				strcpy(p_RunSts->method, "Select the mixed improvement");
				break;

		}
		p_RunSts->runtime[i] = (double)(clock() - start_time) / CLOCKS_PER_SEC;
		p_RunSts->objval[i] = cur_obj;
		p_RunSts->iterations[i] = iters;
		if (p_RunSts->runtime[i] < min_time){
			p_RunSts->best_runtime = p_RunSts->runtime[i];
			min_time = p_RunSts->runtime[i];
		}
		if (p_RunSts->objval[i] > max_obj){
			p_RunSts->best_objval = p_RunSts->objval[i];
			max_obj = p_RunSts->objval[i];
		}
		if (p_RunSts->iterations[i] < min_iter){
			p_RunSts->best_iteration = p_RunSts->iterations[i];
			min_iter = p_RunSts->iterations[i];
		}
		sum_runtime += p_RunSts->runtime[i];
		sum_objval += cur_obj;
		sum_iters += iters;
		assert(cur_obj == calculateSum(pvertex));

		for (int i = 0; i < nnode;i++){
			delete[] gamma_tbl[i];
		}
		delete[] gamma_tbl;
		disposePartition();
	}
	p_RunSts->ave_runtime = sum_runtime / run_cnt;
	p_RunSts->ave_objval = (float)sum_objval / run_cnt;
	p_RunSts->ave_iteration = (float)sum_iters / run_cnt;
}
int main(int argc, char** argv){
	for (int i = 1; i < argc; i+=2){
		if (argv[i][0] != '-' || argv[i][2] != 0){
			showUsage();
			exit(0);
		}else if (argv[i][1] == 'f'){
			strncpy(param_graph_file, argv[i+1],1000);
		}else if (argv[i][1] == 'r'){
			param_run_cnt = atoi(argv[i+1]);
		}else if (argv[i][1] == 'c'){
			param_choice = atoi(argv[i+1]);
		}else if (argv[i][1] == 'k'){
			param_klen = atoi(argv[i+1]);
		}
	}
	if (strlen(param_graph_file) == 0){
		cerr << "No file indicated" << endl;
		exit(1);
	}
	srand((unsigned int)time(NULL));
	loadCompletedGraph(param_graph_file);
	RunStats rs;
	runCompare(param_choice, param_run_cnt,&rs);
	printf("\nChoice:%s\t |V|=%d\n", rs.method, nnode);
	printf("No.\tObjVal\ttime\titeration\t\n");
	for (int i = 0; i < rs.run_cnt; i++){
		printf("%d\t%d\t\%.2f\t%d\t\n",i, rs.objval[i], rs.runtime[i],rs.iterations[i]);
	}
	printf("Best objval:%d, Average Objeval:%f\n",rs.best_objval, rs.ave_objval);
	printf("Best time: %.2f Average time: %f\n",rs.best_runtime, rs.ave_runtime);
	printf("Best iteration: %d Average iteration: %f\n",rs.best_iteration, rs.ave_iteration);

}

