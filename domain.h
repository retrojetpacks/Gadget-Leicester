#ifndef ALLVARS_H
#include "allvars.h"
#endif



void domain_assign_load_or_work_balanced(int mode);
void domain_allocate(void);
void domain_Decomposition(void);
int domain_check_memory_bound(void);
int domain_compare_key(const void *a, const void *b);
int domain_compare_key(const void *a, const void *b);
int domain_compare_toplist(const void *a, const void *b);
double domain_particle_costfactor(int i);
int domain_countToGo(size_t nlimit);
int domain_decompose(void);
int domain_determineTopTree(void);
void domain_exchange(void);
void domain_findExchangeNumbers(int task, int partner, int sphflag, int *send, int *recv);
void domain_findExtent(void);
void domain_findSplit(int cpustart, int ncpu, int first, int last);
void domain_findSplit_balanced(int cpustart, int ncpu, int first, int last);
void domain_free(void);
void domain_shiftSplit(void);
void domain_sumCost(void);
int domain_topsplit(int node, peanokey startkey);
int domain_topsplit_local(int node, peanokey startkey, int mode);
int domain_topsplit_special(void);
int domain_compare_key(const void *a, const void *b);
int domain_check_for_local_refine(int i, double countlimit, double costlimit);
void domain_free_trick(void);
void domain_allocate_trick(void);
int domain_recursively_combine_topTree(int start, int ncpu);


