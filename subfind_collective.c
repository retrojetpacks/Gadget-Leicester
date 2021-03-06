#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"
#include "domain.h"

#ifdef SUBFIND
#include "subfind.h"
#include "fof.h"


#define TAG_GET_NGB_COUNT        200
#define TAG_POLLING_DONE         201
#define TAG_SET_ALL              202
#define TAG_SET_ALL_DATA         203
#define TAG_GET_NGB_INDICES      204
#define TAG_GET_TAILANDLEN       205
#define TAG_GET_TAILANDLEN_DATA  206
#define TAG_SET_TAILANDLEN       207
#define TAG_SET_TAILANDLEN_DATA  208
#define TAG_SET_HEADANDNEXT      209
#define TAG_SET_HEADANDNEXT_DATA 210
#define TAG_SET_NEXT             211
#define TAG_SET_NEXT_DATA        212
#define TAG_SETHEADGETNEXT       213
#define TAG_SETHEADGETNEXT_DATA  214
#define TAG_GET_NEXT             215
#define TAG_GET_NEXT_DATA        216
#define TAG_GET_HEAD             217
#define TAG_GET_HEAD_DATA        218
#define TAG_ADD_PARTICLE         219
#define TAG_ADDBOUND             220
#define TAG_ADDBOUND_DATA        221
#define TAG_NID                  222
#define TAG_NID_DATA             223
#define TAG_SETRANK              224
#define TAG_SETRANK_IN           225
#define TAG_SETRANK_OUT          226
#define TAG_GET_RANK             227
#define TAG_GET_RANK_DATA        228
#define TAG_MARK_PARTICLE        229
#define TAG_MARK_PARTICLE_DATA   230


#define MASK ((((long long)1)<< 32)-1)
#define HIGHBIT (1 << 30)

#define NEAREST(x) (((x)>boxhalf)?((x)-boxsize):(((x)<-boxhalf)?((x)+boxsize):(x)))




static long long *Head, *Next, *Tail;
static int *Len;
static int LocalLen;

static struct cand_dat
{
  long long head;
  long long rank;
  int len;
  int nsub;
  int subnr, parent;
  int bound_length;
}
 *candidates;



static struct unbind_data *ud;


static struct sort_density_data
{
  double density;
  long long index;		/* this will store the task in the upper word */
}
 *sd;


void subfind_unbind_independent_ones(int count_cand)
{
  int i, j, k, len, nsubs;

  /*subfind_loctree_treeallocate(All.TreeAllocFactor * NumPart, NumPart); */

  R2list = mymalloc(NumPart * sizeof(struct r2data));
  ud = mymalloc(NumPart * sizeof(struct unbind_data));

  qsort(candidates, count_cand, sizeof(struct cand_dat), subfind_compare_candidates_nsubs);

  for(k = 0, i = 0; k < count_cand; k++)
    if(candidates[k].parent == 0)
      {
	while(P[i].submark < candidates[k].nsub)
	  {
	    i++;
	    if(i >= NumPart)
	      endrun(13213);
	  }

	if(P[i].submark >= 0 && P[i].submark < HIGHBIT)
	  {
	    len = 0;
	    nsubs = P[i].submark;

	    if(nsubs != candidates[k].nsub)
	      {
		printf("TASK=%d i=%d k=%d nsubs=%d candidates[k].nsub=%d\n",
		       ThisTask, i, k, nsubs, candidates[k].nsub);
		endrun(13199);
	      }

	    while(i < NumPart)
	      {
		if(P[i].submark == nsubs)
		  {
		    P[i].submark = HIGHBIT;
		    if((P[i].origintask & HIGHBIT) == 0)
		      {
			ud[len].index = i;
			len++;
		      }
		    i++;
		  }
		else
		  break;
	      }

	    len = subfind_unbind(ud, len);

	    if(len >= All.DesLinkNgb)
	      {
		/* ok, we found a substructure */
		candidates[k].bound_length = len;

		for(j = 0; j < len; j++)
		  P[ud[j].index].submark = nsubs;	/* we use this to flag the substructures */
	      }
	    else
	      candidates[k].bound_length = 0;
	  }
      }

  myfree(ud);
  myfree(R2list);
}

void subfind_process_group_collectively(void)
{
  long long p, ss, head, head_attach, ngb_index1, ngb_index2, rank;
  long long prev, tail, tail_attach, tmp, next, index;
  int len, len_attach, totgrouplen1, totgrouplen2;
  int max_candidates, ncand, parent, totcand;
  int max_loc_length, max_length;
  int count, countall, *countlist, *offset;
  int i, j, k, grindex = 0, count_cand, nsubs, subnr;
  int count_leaves, tot_count_leaves;
  int master, ngbcount;
  double SubMass, SubPos[3], SubVel[3], SubCM[3], SubVelDisp, SubVmax, SubVmaxRad, SubSpin[3], SubHalfMass;
  struct cand_dat *tmp_candidates;
  MyIDType SubMostBoundID;
  double t0, t1;


  if(ThisTask == 0)
    printf("\ncollectively doing halo %d\n", GrNr);

  for(i = 0, NumPartGroup = 0; i < NumPart; i++)
    if(P[i].GrNr == GrNr)
      NumPartGroup++;

  MPI_Allreduce(&NumPartGroup, &totgrouplen1, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == ((GrNr - 1) % NTask))
    {
      for(grindex = 0; grindex < Ngroups; grindex++)
	if(Group[grindex].GrNr == GrNr)
	  break;
      if(grindex >= Ngroups)
	endrun(8);
      totgrouplen2 = Group[grindex].Len;
    }

  MPI_Bcast(&totgrouplen2, 1, MPI_INT, (GrNr - 1) % NTask, MPI_COMM_WORLD);

  if(totgrouplen1 != totgrouplen2)
    endrun(9);			/* inconsistency */


  qsort(P, NumPart, sizeof(struct particle_data), subfind_compare_P_GrNrGrNr);


  /* distribute this halo among the processors */
  t0 = second();

  All.DoDynamicUpdate = 0;

  domain_free_trick();

  domain_Decomposition();

  t1 = second();
  if(ThisTask == 0)
    printf("coldomain_Decomposition() took %g sec\n", timediff(t0, t1));

  qsort(P, NumPart, sizeof(struct particle_data), subfind_compare_P_GrNrGrNr);

  for(i = 0, NumPartGroup = 0; i < NumPart; i++)
    if(P[i].GrNr == GrNr)
      NumPartGroup++;

  subfind_loctree_copyExtent();	/* this will make sure that all the serial trees start from the same root node geometry */

  /* construct a tree for the halo */
  t0 = second();
  force_treebuild(NumPartGroup, NULL);
  t1 = second();
  if(ThisTask == 0)
    printf("force_treebuild() took %g sec\n", timediff(t0, t1));


  /* determine the radius that encloses a certain number of link particles */
  t0 = second();
  subfind_find_linkngb();
  t1 = second();
  if(ThisTask == 0)
    printf("find_linkngb() took %g sec\n", timediff(t0, t1));


  /* determine the indices of the nearest two denser neighbours within the link region */
  t0 = second();
  NgbLoc = mymalloc(NumPartGroup * sizeof(struct nearest_ngb_data));
  R2Loc = mymalloc(NumPartGroup * sizeof(struct nearest_r2_data));
  subfind_find_nearesttwo();
  myfree(R2Loc);
  t1 = second();
  if(ThisTask == 0)
    printf("find_nearesttwo() took %g sec\n", timediff(t0, t1));


  /* sort the densities */
  sd = mymalloc(NumPartGroup * sizeof(struct sort_density_data));
  for(i = 0; i < NumPartGroup; i++)
    {
      sd[i].density = P[i].u.DM_Density;
      sd[i].index = (((long long) ThisTask) << 32) + i;
    }
  t0 = second();
  parallel_sort(sd, NumPartGroup, sizeof(struct sort_density_data), subfind_compare_densities);
  t1 = second();
  if(ThisTask == 0)
    printf("parallel sort of densities done took %g sec\n", timediff(t0, t1));



  /* allocate and initialize distributed link list */
  Head = mymalloc(NumPartGroup * sizeof(long long));
  Next = mymalloc(NumPartGroup * sizeof(long long));
  Tail = mymalloc(NumPartGroup * sizeof(long long));
  Len = mymalloc(NumPartGroup * sizeof(int));

  for(i = 0; i < NumPartGroup; i++)
    {
      Head[i] = Next[i] = Tail[i] = -1;
      Len[i] = 0;
    }


  /* allocate a list to store subhalo candidates */
  max_candidates = (NumPartGroup / All.DesLinkNgb);
  candidates = mymalloc(max_candidates * sizeof(struct cand_dat));
  count_cand = 0;


  /* now find the subhalo candidates by building up link lists from high density to low density */
  t0 = second();
  for(master = 0; master < NTask; master++)
    {
      if(ThisTask != master)
	subfind_poll_for_requests();
      else
	{
	  for(k = 0; k < NumPartGroup; k++)
	    {
	      ngbcount = subfind_distlinklist_get_ngb_count(sd[k].index, &ngb_index1, &ngb_index2);

	      switch (ngbcount)	/* treat the different possible cases */
		{
		case 0:	/* this appears to be a lonely maximum -> new group */
		  subfind_distlinklist_set_all(sd[k].index, sd[k].index, sd[k].index, 1, -1);
		  break;

		case 1:	/* the particle is attached to exactly one group */
		  head = subfind_distlinklist_get_head(ngb_index1);
		  subfind_distlinklist_get_tailandlen(head, &tail, &len);
		  subfind_distlinklist_set_tailandlen(head, sd[k].index, len + 1);
		  subfind_distlinklist_set_headandnext(sd[k].index, head, -1);
		  subfind_distlinklist_set_next(tail, sd[k].index);
		  break;

		case 2:	/* the particle merges two groups together */
		  head = subfind_distlinklist_get_head(ngb_index1);
		  head_attach = subfind_distlinklist_get_head(ngb_index2);

		  if(head != head_attach)
		    {
		      subfind_distlinklist_get_tailandlen(head, &tail, &len);
		      subfind_distlinklist_get_tailandlen(head_attach, &tail_attach, &len_attach);

		      if(len_attach > len)	/* other group is longer, swap them */
			{
			  tmp = head;
			  head = head_attach;
			  head_attach = tmp;
			  tmp = tail;
			  tail = tail_attach;
			  tail_attach = tmp;
			  tmp = len;
			  len = len_attach;
			  len_attach = tmp;
			}

		      /* only in case the attached group is long enough we bother to register it 
		         as a subhalo candidate */

		      if(len_attach >= All.DesLinkNgb)
			{
			  if(count_cand < max_candidates)
			    {
			      candidates[count_cand].len = len_attach;
			      candidates[count_cand].head = head_attach;
			      count_cand++;
			    }
			  else
			    endrun(87);
			}

		      /* now join the two groups */
		      subfind_distlinklist_set_tailandlen(head, tail_attach, len + len_attach);
		      subfind_distlinklist_set_next(tail, head_attach);

		      ss = head_attach;
		      do
			{
			  ss = subfind_distlinklist_set_head_get_next(ss, head);
			}
		      while(ss >= 0);
		    }

		  /* finally, attach the particle to 'head' */
		  subfind_distlinklist_get_tailandlen(head, &tail, &len);
		  subfind_distlinklist_set_tailandlen(head, sd[k].index, len + 1);
		  subfind_distlinklist_set_headandnext(sd[k].index, head, -1);
		  subfind_distlinklist_set_next(tail, sd[k].index);
		  break;
		}
	    }

	  /* now tell the others to stop polling */
	  for(k = 0; k < NTask; k++)
	    if(k != ThisTask)
	      MPI_Send(&k, 1, MPI_INT, k, TAG_POLLING_DONE, MPI_COMM_WORLD);
	}
    }
  t1 = second();
  if(ThisTask == 0)
    printf("identification of primary candidates took %g sec\n", timediff(t0, t1));

  /* add the full thing as a subhalo candidate */
  t0 = second();
  for(master = 0, head = -1, prev = -1; master < NTask; master++)
    {
      if(ThisTask != master)
	subfind_poll_for_requests();
      else
	{
	  for(i = 0; i < NumPartGroup; i++)
	    {
	      index = (((long long) ThisTask) << 32) + i;

	      if(Head[i] == index)
		{
		  subfind_distlinklist_get_tailandlen(Head[i], &tail, &len);
		  next = subfind_distlinklist_get_next(tail);
		  if(next == -1)
		    {
		      if(prev < 0)
			head = index;

		      if(prev >= 0)
			subfind_distlinklist_set_next(prev, index);

		      prev = tail;
		    }
		}
	    }

	  /* now tell the others to stop polling */
	  for(k = 0; k < NTask; k++)
	    if(k != ThisTask)
	      MPI_Send(&k, 1, MPI_INT, k, TAG_POLLING_DONE, MPI_COMM_WORLD);
	}

      MPI_Bcast(&head, sizeof(head), MPI_BYTE, master, MPI_COMM_WORLD);
      MPI_Bcast(&prev, sizeof(prev), MPI_BYTE, master, MPI_COMM_WORLD);
    }

  if(ThisTask == NTask - 1)
    {
      if(count_cand < max_candidates)
	{
	  candidates[count_cand].len = totgrouplen1;
	  candidates[count_cand].head = head;
	  count_cand++;
	}
      else
	endrun(123123);
    }
  t1 = second();
  if(ThisTask == 0)
    printf("adding background as candidate took %g sec\n", timediff(t0, t1));

  /* go through the whole chain once to establish a rank order. For the rank we use Len[] */
  t0 = second();

  master = (head >> 32);

  if(ThisTask != master)
    subfind_poll_for_requests();
  else
    {
      p = head;
      rank = 0;

      while(p >= 0)
	{
	  p = subfind_distlinklist_setrank_and_get_next(p, &rank);
	}

      /* now tell the others to stop polling */
      for(i = 0; i < NTask; i++)
	if(i != master)
	  MPI_Send(&i, 1, MPI_INT, i, TAG_POLLING_DONE, MPI_COMM_WORLD);
    }

  MPI_Bcast(&rank, sizeof(rank), MPI_BYTE, master, MPI_COMM_WORLD);	/* just for testing */

  /* for each candidate, we now pull out the rank of its head */
  for(master = 0; master < NTask; master++)
    {
      if(ThisTask != master)
	subfind_poll_for_requests();
      else
	{
	  for(k = 0; k < count_cand; k++)
	    candidates[k].rank = subfind_distlinklist_get_rank(candidates[k].head);

	  /* now tell the others to stop polling */
	  for(i = 0; i < NTask; i++)
	    if(i != ThisTask)
	      MPI_Send(&i, 1, MPI_INT, i, TAG_POLLING_DONE, MPI_COMM_WORLD);
	}
    }

  t1 = second();
  if(ThisTask == 0)
    printf("establishing of rank order took %g sec  (p=%d)\n", timediff(t0, t1), (int) rank);


  /* establish total number of candidates */
  MPI_Allreduce(&count_cand, &totcand, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(ThisTask == 0)
    printf("\ntotal number of subhalo candidates=%d\n", totcand);


  for(i = 0; i < NumPartGroup; i++)
    Tail[i] = -1;

  for(i = 0; i < count_cand; i++)
    candidates[i].parent = 0;

  do
    {
      /* Let's see which candidates can be unbound independent from each other.
         We identify them with those candidates that have no embedded other candidate */
      t0 = second();
      if(ThisTask == 0)
	tmp_candidates = mymalloc(totcand * sizeof(struct cand_dat));

      count = count_cand;
      count *= sizeof(struct cand_dat);

      countlist = mymalloc(NTask * sizeof(int));
      offset = mymalloc(NTask * sizeof(int));

      MPI_Allgather(&count, 1, MPI_INT, countlist, 1, MPI_INT, MPI_COMM_WORLD);

      for(i = 1, offset[0] = 0; i < NTask; i++)
	offset[i] = offset[i - 1] + countlist[i - 1];

      MPI_Gatherv(candidates, countlist[ThisTask], MPI_BYTE,
		  tmp_candidates, countlist, offset, MPI_BYTE, 0, MPI_COMM_WORLD);

      if(ThisTask == 0)
	{
	  for(k = 0; k < totcand; k++)
	    {
	      tmp_candidates[k].nsub = k;
	      tmp_candidates[k].subnr = k;
	    }

	  qsort(tmp_candidates, totcand, sizeof(struct cand_dat), subfind_compare_candidates_rank);

	  for(k = 0; k < totcand; k++)
	    {
	      if(tmp_candidates[k].parent >= 0)
		{
		  tmp_candidates[k].parent = 0;

		  for(j = k + 1; j < totcand; j++)
		    {
		      if(tmp_candidates[j].rank > tmp_candidates[k].rank + tmp_candidates[k].len)
			break;

		      if(tmp_candidates[j].parent < 0)	/* ignore these */
			continue;

		      if(tmp_candidates[k].rank + tmp_candidates[k].len >=
			 tmp_candidates[j].rank + tmp_candidates[j].len)
			{
			  tmp_candidates[k].parent++;	/* we here count the number of subhalos that are enclosed */
			}
		      else
			{
			  printf("k=%d|%d has rank=%d and len=%d.  j=%d has rank=%d and len=%d bound=%d\n",
				 k, countall, (int) tmp_candidates[k].rank, (int) tmp_candidates[k].len,
				 (int) tmp_candidates[k].bound_length, (int) tmp_candidates[j].rank,
				 (int) tmp_candidates[j].len, (int) tmp_candidates[j].bound_length);
			  endrun(8812313);
			}
		    }
		}
	    }
	  qsort(tmp_candidates, totcand, sizeof(struct cand_dat), subfind_compare_candidates_subnr);
	}

      MPI_Scatterv(tmp_candidates, countlist, offset, MPI_BYTE,
		   candidates, countlist[ThisTask], MPI_BYTE, 0, MPI_COMM_WORLD);


      myfree(offset);
      myfree(countlist);

      if(ThisTask == 0)
	myfree(tmp_candidates);


      for(i = 0, count_leaves = 0, max_loc_length = 0; i < count_cand; i++)
	if(candidates[i].parent == 0)
	  {
	    count_leaves++;
	    if(candidates[i].len > max_loc_length)
	      max_loc_length = candidates[i].len;
	  }

      MPI_Allreduce(&count_leaves, &tot_count_leaves, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&max_loc_length, &max_length, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

      t1 = second();
      if(ThisTask == 0)
	printf
	  ("\nnumber of subhalo candidates that can be done independently=%d.\n(Largest size is %d, finding them took %g sec)\n",
	   tot_count_leaves, max_length, timediff(t0, t1));

      if(tot_count_leaves < NTask)	/* if there are only a few left, let's do them collectively */
	{
	  if(ThisTask == 0)
	    printf("too few, I do the rest collectively\n\n");
	  break;
	}

      if(max_length > 0.1 * All.TotNumPart / NTask)	/* seems large, let's rather do it collectively */
	{
	  if(ThisTask == 0)
	    printf("too big candidates, I do the rest collectively\n\n");
	  break;
	}

      for(i = 0; i < NumPart; i++)
	{
	  P[i].origintask = P[i].targettask = ThisTask;
	  P[i].submark = HIGHBIT;
	  if(i < NumPartGroup)
	    if(Tail[i] >= 0)	/* this means this particle is already bound to a substructure */
	      P[i].origintask |= HIGHBIT;
	}

      /* we now mark the particles that are in subhalo candidates that can be processed independently in parallel */
      nsubs = 0;
      t0 = second();
      for(master = 0; master < NTask; master++)
	{
	  ncand = count_cand;

	  MPI_Bcast(&ncand, sizeof(ncand), MPI_BYTE, master, MPI_COMM_WORLD);

	  for(k = 0; k < ncand; k++)
	    {
	      if(ThisTask == master)
		{
		  len = candidates[k].len;
		  parent = candidates[k].parent;	/* this is here actually the daughter count */
		}

	      MPI_Bcast(&len, sizeof(len), MPI_BYTE, master, MPI_COMM_WORLD);
	      MPI_Bcast(&parent, sizeof(parent), MPI_BYTE, master, MPI_COMM_WORLD);

	      if(parent == 0)
		{
		  if(ThisTask != master)
		    subfind_poll_for_requests();
		  else
		    {
		      for(i = 0, p = candidates[k].head; i < candidates[k].len; i++)
			{
			  subfind_distlinklist_mark_particle(p, master, nsubs);

			  if(p < 0)
			    {
			      printf("Bummer i=%d \n", i);
			      endrun(128);
			    }
			  p = subfind_distlinklist_get_next(p);
			}

		      /* now tell the others to stop polling */
		      for(i = 0; i < NTask; i++)
			if(i != ThisTask)
			  MPI_Send(&i, 1, MPI_INT, i, TAG_POLLING_DONE, MPI_COMM_WORLD);
		    }
		}

	      nsubs++;
	    }
	}
      t1 = second();
      if(ThisTask == 0)
	{
	  printf("particles are marked (took %g)\n", timediff(t0, t1));
	  fflush(stdout);
	}

      t0 = second();
      subfind_distribute_particles(1);	/* assemble the particles on individual processors */
      t1 = second();

      if(ThisTask == 0)
	{
	  printf("independent subhalos are assembled on individual CPUs for unbinding (%g sec)\n",
		 timediff(t0, t1));
	  fflush(stdout);
	}

      qsort(P, NumPart, sizeof(struct particle_data), subfind_compare_P_submark);	/* groups particles of the same canidate together */

      MPI_Barrier(MPI_COMM_WORLD);
      t0 = second();

      subfind_unbind_independent_ones(count_cand);

      MPI_Barrier(MPI_COMM_WORLD);
      t1 = second();

      if(ThisTask == 0)
	{
	  printf("unbinding of independent ones took %g sec\n", timediff(t0, t1));
	  fflush(stdout);
	}

      for(i = 0; i < NumPart; i++)
	P[i].origintask &= (HIGHBIT - 1);	/* clear high bit if set */

      t0 = second();
      subfind_distribute_particles(2);	/* bring them back to their original processor */
      t1 = second();


      if(ThisTask == 0)
	printf("particles have returned to their original processor (%g sec)\n", timediff(t0, t1));

      /* reestablish the original order */
      qsort(P, NumPart, sizeof(struct particle_data), subfind_compare_P_GrNrGrNr);


      /* now mark the bound particles */
      for(i = 0; i < NumPartGroup; i++)
	if(P[i].submark >= 0 && P[i].submark < nsubs)
	  Tail[i] = P[i].submark;	/* we use this to flag bound parts of substructures */

      for(i = 0; i < count_cand; i++)
	if(candidates[i].parent == 0)
	  candidates[i].parent = -1;
    }
  while(tot_count_leaves > 0);





  force_treebuild(NumPartGroup, NULL);	/* re construct the tree for the collective part */


  /**** now we do the collective unbinding of the subhalo candidates that contain other subhalo candidates ****/
  ud = mymalloc(NumPartGroup * sizeof(struct unbind_data));
  
  t0 = second();
  for(master = 0; master < NTask; master++)
    {
      ncand = count_cand;

      MPI_Bcast(&ncand, sizeof(ncand), MPI_BYTE, master, MPI_COMM_WORLD);

      for(k = 0; k < ncand; k++)
	{
	  if(ThisTask == master)
	    {
	      len = candidates[k].len;
	      nsubs = candidates[k].nsub;
	      parent = candidates[k].parent;	/* this is here actually the daughter count */
	    }

	  MPI_Bcast(&parent, sizeof(parent), MPI_BYTE, master, MPI_COMM_WORLD);

	  if(parent >= 0)
	    {
	      MPI_Bcast(&len, sizeof(len), MPI_BYTE, master, MPI_COMM_WORLD);
	      MPI_Bcast(&nsubs, sizeof(nsubs), MPI_BYTE, master, MPI_COMM_WORLD);

	      if(ThisTask == 0)
		printf("collective unbinding of cand=%d|%d on task=%d|%d of length=%d\n",
		       k, ncand, master, NTask, (int) len);

	      LocalLen = 0;

	      if(ThisTask != master)
		subfind_poll_for_requests();
	      else
		{
		  for(i = 0, p = candidates[k].head; i < candidates[k].len; i++)
		    {
		      subfind_distlinklist_add_particle(p);
		      if(p < 0)
			{
			  printf("Bummer i=%d \n", i);
			  endrun(123);

			}
		      p = subfind_distlinklist_get_next(p);
		    }

		  /* now tell the others to stop polling */
		  for(i = 0; i < NTask; i++)
		    if(i != ThisTask)
		      MPI_Send(&i, 1, MPI_INT, i, TAG_POLLING_DONE, MPI_COMM_WORLD);
		}

	      LocalLen = subfind_col_unbind(ud, LocalLen);

	      MPI_Allreduce(&LocalLen, &len, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	      if(len >= All.DesLinkNgb)
		{
		  /* ok, we found a substructure */
#ifdef VERBOSE
		  if(ThisTask == 0)
		    printf("substructure of len=%d found\n", (int) len);
#endif
		  for(i = 0; i < LocalLen; i++)
		    Tail[ud[i].index] = nsubs;	/* we use this to flag the substructures */

		  if(ThisTask == master)
		    {
		      candidates[k].bound_length = len;
		    }
		}
	      else
		{
		  if(ThisTask == master)
		    {
		      candidates[k].bound_length = 0;
		    }
		}
	    }
	}
    }
  t1 = second();
  if(ThisTask == 0)
    printf("the collective unbinding of remaining halos took %g sec\n", timediff(t0, t1));


  for(k = 0, count = 0; k < count_cand; k++)
    if(candidates[k].bound_length >= All.DesLinkNgb)
      {
	if(candidates[k].len < All.DesLinkNgb)
	  {
	    printf("candidates[k=%d].len=%d bound=%d\n", k, candidates[k].len, candidates[k].bound_length);
	    endrun(77);
	  }
	count++;
      }
  
  MPI_Allreduce(&count, &countall, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  
  if(ThisTask == 0)
    {
      printf("\nfound %d bound substructures in FoF group of length %d\n", countall, totgrouplen1);
      fflush(stdout);
    }



  /* now determine the parent subhalo for each candidate */
  t0 = second();
  parallel_sort(candidates, count_cand, sizeof(struct cand_dat), subfind_compare_candidates_boundlength);

  if(ThisTask == 0)
    tmp_candidates = mymalloc(totcand * sizeof(struct cand_dat));

  count = count_cand;
  count *= sizeof(struct cand_dat);

  countlist = mymalloc(NTask * sizeof(int));
  offset = mymalloc(NTask * sizeof(int));

  MPI_Allgather(&count, 1, MPI_INT, countlist, 1, MPI_INT, MPI_COMM_WORLD);

  for(i = 1, offset[0] = 0; i < NTask; i++)
    offset[i] = offset[i - 1] + countlist[i - 1];

  MPI_Gatherv(candidates, countlist[ThisTask], MPI_BYTE,
	      tmp_candidates, countlist, offset, MPI_BYTE, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      for(k = 0; k < totcand; k++)
	{
	  tmp_candidates[k].subnr = k;
	  tmp_candidates[k].parent = 0;
	}

      qsort(tmp_candidates, totcand, sizeof(struct cand_dat), subfind_compare_candidates_rank);

      for(k = 0; k < totcand; k++)  
	{
	  for(j = k + 1; j < totcand; j++)
	    {
	      if(tmp_candidates[j].rank > tmp_candidates[k].rank + tmp_candidates[k].len)
		break;

	      if(tmp_candidates[k].rank + tmp_candidates[k].len >=
		 tmp_candidates[j].rank + tmp_candidates[j].len)
		{
		  if(tmp_candidates[k].bound_length >= All.DesLinkNgb)
		    tmp_candidates[j].parent = tmp_candidates[k].subnr;
		}
	      else
		{
		  printf("k=%d|%d has rank=%d and len=%d.  j=%d has rank=%d and len=%d bound=%d\n",
			 k, countall, (int) tmp_candidates[k].rank, (int) tmp_candidates[k].len,
			 (int) tmp_candidates[k].bound_length, (int) tmp_candidates[j].rank,
			 (int) tmp_candidates[j].len, (int) tmp_candidates[j].bound_length);
		  endrun(1212313);
		}
	    }
	}

      qsort(tmp_candidates, totcand, sizeof(struct cand_dat), subfind_compare_candidates_subnr);
    }

  MPI_Scatterv(tmp_candidates, countlist, offset, MPI_BYTE,
	       candidates, countlist[ThisTask], MPI_BYTE, 0, MPI_COMM_WORLD);


  myfree(offset);
  myfree(countlist);

  if(ThisTask == 0)
    myfree(tmp_candidates);

  t1 = second();
  if(ThisTask == 0)
    printf("determination of parent subhalo took %g sec\n", timediff(t0, t1));


  /* Now let's save  some properties of the substructures */
  if(ThisTask == ((GrNr - 1) % NTask))
    {
      Group[grindex].Nsubs = countall;
    }

  t0 = second();
  for(master = 0, subnr = 0; master < NTask; master++)
    {
      ncand = count_cand;
      MPI_Bcast(&ncand, sizeof(int), MPI_INT, master, MPI_COMM_WORLD);

      for(k = 0; k < ncand; k++)
	{
	  if(ThisTask == master)
	    {
	      len = candidates[k].bound_length;
	      nsubs = candidates[k].nsub;
	      parent = candidates[k].parent;
	    }

	  MPI_Bcast(&len, sizeof(len), MPI_BYTE, master, MPI_COMM_WORLD);


	  if(len > 0)
	    {
	      MPI_Bcast(&nsubs, sizeof(nsubs), MPI_BYTE, master, MPI_COMM_WORLD);
	      MPI_Bcast(&parent, sizeof(parent), MPI_BYTE, master, MPI_COMM_WORLD);
	  
	      LocalLen = 0;

	      if(ThisTask != master)
		subfind_poll_for_requests();
	      else
		{
		  for(i = 0, p = candidates[k].head; i < candidates[k].len; i++)
		    {
		      subfind_distlinklist_add_bound_particles(p, nsubs);
		      p = subfind_distlinklist_get_next(p);
		    }

		  /* now tell the others to stop polling */
		  for(i = 0; i < NTask; i++)
		    if(i != ThisTask)
		      MPI_Send(&i, 1, MPI_INT, i, TAG_POLLING_DONE, MPI_COMM_WORLD);
		}

	      subfind_col_determine_sub_halo_properties(ud, LocalLen, &SubMass,
							&SubPos[0], &SubVel[0], &SubCM[0], &SubVelDisp,
							&SubVmax, &SubVmaxRad, &SubSpin[0], &SubMostBoundID,
							&SubHalfMass);

	      /* we have filled into ud the binding energy and the particle ID return */

	      if(((GrNr - 1) % NTask) == ThisTask)
		{
		  if(Nsubgroups >= MaxNsubgroups)
		    endrun(899);

		  if(subnr == 0)
		    {
		      for(j = 0; j < 3; j++)
			Group[grindex].Pos[j] = SubPos[j];
		    }

		  SubGroup[Nsubgroups].Len = len;
		  if(subnr == 0)
		    SubGroup[Nsubgroups].Offset = Group[grindex].Offset;
		  else
		    SubGroup[Nsubgroups].Offset =
		      SubGroup[Nsubgroups - 1].Offset + SubGroup[Nsubgroups - 1].Len;
		  SubGroup[Nsubgroups].GrNr = GrNr - 1;
		  SubGroup[Nsubgroups].SubNr = subnr;
		  SubGroup[Nsubgroups].SubParent = parent;
		  SubGroup[Nsubgroups].Mass = SubMass;
		  SubGroup[Nsubgroups].SubMostBoundID = SubMostBoundID;
		  SubGroup[Nsubgroups].SubVelDisp = SubVelDisp;
		  SubGroup[Nsubgroups].SubVmax = SubVmax;
		  SubGroup[Nsubgroups].SubVmaxRad = SubVmaxRad;
		  SubGroup[Nsubgroups].SubHalfMass = SubHalfMass;

		  for(j = 0; j < 3; j++)
		    {
		      SubGroup[Nsubgroups].Pos[j] = SubPos[j];
		      SubGroup[Nsubgroups].Vel[j] = SubVel[j];
		      SubGroup[Nsubgroups].CM[j] = SubCM[j];
		      SubGroup[Nsubgroups].Spin[j] = SubSpin[j];
		    }
		  
		  Nsubgroups++;
		}

	      /* Let's now assign the subgroup number */

	      for(i = 0; i < LocalLen; i++)
		{
		  P[ud[i].index].SubNr = subnr;
		}
	      
	      subnr++;
	    }
	}
    }
  t1 = second();
  if(ThisTask == 0)
    printf("determining substructure properties took %g sec\n", timediff(t0, t1));

  myfree(ud);
  myfree(candidates);
  myfree(Len);
  myfree(Tail);
  myfree(Next);
  myfree(Head);
  myfree(sd);
  myfree(NgbLoc);

  force_treefree();
  domain_free();
  domain_allocate_trick();
}





int subfind_col_unbind(struct unbind_data *d, int num)
{
  int iter = 0;
  int i, j, p, part_index, minindex, task;
  int unbound, totunbound, numleft, mincpu;
  int *npart, *offset, *nbu_count, count_bound_unbound, phaseflag;
  double s[3], dx[3], ddxx, v[3], dv[3], sloc[3], vloc[3], pos[3];
  double vel_to_phys, H_of_a;
  MyFloat minpot, *potlist;
  double boxsize, boxhalf;
  double mass, massloc, t0, t1;
  double *bnd_energy, energy_limit, energy_limit_local, weakly_bound_limit_local, weakly_bound_limit = 0;

  boxsize = All.BoxSize;
  boxhalf = 0.5 * All.BoxSize;

  vel_to_phys = 1.0 / All.Time;

  H_of_a = hubble_function(All.Time);

  phaseflag = 0;		/* this means we will recompute the potential for all particles */

  do
    {
      t0 = second();

      force_treebuild(num, d);

      /* let's compute the potential energy */

      subfind_potential_compute(num, d, phaseflag, weakly_bound_limit);

      if(phaseflag == 0)
	{
	  potlist = (MyFloat *) mymalloc(NTask * sizeof(MyFloat));

	  for(i = 0, minindex = -1, minpot = 1.0e30; i < num; i++)
	    {
	      if(P[d[i].index].u.DM_Potential < minpot || minindex == -1)
		{
		  minpot = P[d[i].index].u.DM_Potential;
		  minindex = d[i].index;
		}
	    }

	  MPI_Allgather(&minpot, sizeof(MyFloat), MPI_BYTE, potlist, sizeof(MyFloat), MPI_BYTE,
			MPI_COMM_WORLD);

	  for(i = 0, mincpu = -1, minpot = 1.0e30; i < NTask; i++)
	    if(potlist[i] < minpot)
	      {
		mincpu = i;
		minpot = potlist[i];
	      }

	  if(mincpu < 0)
	    endrun(112);

	  myfree(potlist);

	  if(ThisTask == mincpu)
	    {
	      for(j = 0; j < 3; j++)
		pos[j] = P[minindex].Pos[j];
	    }

	  MPI_Bcast(&pos[0], 3, MPI_DOUBLE, mincpu, MPI_COMM_WORLD);
	  /* pos[] now holds the position of minimum potential */
	  /* we take that as the center */
	}

      /* let's get bulk velocity and the center-of-mass */

      for(j = 0; j < 3; j++)
	sloc[j] = vloc[j] = 0;

      for(i = 0, massloc = 0; i < num; i++)
	{
	  part_index = d[i].index;

	  for(j = 0; j < 3; j++)
	    {
#ifdef PERIODIC
	      ddxx = NEAREST(P[part_index].Pos[j] - pos[j]);
#else
	      ddxx = P[part_index].Pos[j] - pos[j];
#endif
	      sloc[j] += P[part_index].Mass * ddxx;
	      vloc[j] += P[part_index].Mass * P[part_index].Vel[j];
	    }
	  massloc += P[part_index].Mass;
	}

      MPI_Allreduce(sloc, s, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(vloc, v, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&massloc, &mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      for(j = 0; j < 3; j++)
	{
	  s[j] /= mass;		/* center of mass */
	  v[j] /= mass;

	  s[j] += pos[j];

#ifdef PERIODIC
	  while(s[j] < 0)
	    s[j] += boxsize;
	  while(s[j] >= boxsize)
	    s[j] -= boxsize;
#endif
	}

      bnd_energy = mymalloc(num * sizeof(double));

      for(i = 0; i < num; i++)
	{
	  part_index = d[i].index;

	  for(j = 0; j < 3; j++)
	    {
	      dv[j] = vel_to_phys * (P[part_index].Vel[j] - v[j]);
#ifdef PERIODIC
	      dx[j] = All.Time * NEAREST(P[part_index].Pos[j] - s[j]);
#else
	      dx[j] = All.Time * (P[part_index].Pos[j] - s[j]);
#endif
	      dv[j] += H_of_a * dx[j];
	    }

	  P[part_index].v.DM_BindingEnergy =
	    P[part_index].u.DM_Potential + 0.5 * (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);

	  bnd_energy[i] = P[part_index].v.DM_BindingEnergy;
	}

      parallel_sort(bnd_energy, num, sizeof(double), subfind_compare_binding_energy);


      npart = mymalloc(NTask * sizeof(int));
      nbu_count = mymalloc(NTask * sizeof(int));
      offset = mymalloc(NTask * sizeof(int));

      MPI_Allgather(&num, 1, MPI_INT, npart, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Allreduce(&num, &numleft, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      for(i = 1, offset[0] = 0; i < NTask; i++)
	offset[i] = offset[i - 1] + npart[i - 1];

      j = (int)(0.25 * numleft);  /* index of limiting energy value */
 
      task = 0;
      while(j >= npart[task])
	{
	  j -= npart[task];
	  task++;
	}
	 
      if(ThisTask == task)
	energy_limit_local = bnd_energy[j];
      else
	energy_limit_local = 1.0e30; 
	
      MPI_Allreduce(&energy_limit_local, &energy_limit, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

      for(i = 0, count_bound_unbound = 0; i < num; i++)
	{
	  if(bnd_energy[i] > 0)
	    count_bound_unbound++;
	  else
	    count_bound_unbound--;
	}

      MPI_Allgather(&count_bound_unbound, 1, MPI_INT, nbu_count, 1, MPI_INT, MPI_COMM_WORLD);

      for(i = 0, count_bound_unbound = 0; i < ThisTask; i++)
	count_bound_unbound += nbu_count[i];

      for(i = 0; i < num - 1; i++)
	{
	  if(bnd_energy[i] > 0)
	    count_bound_unbound++;
	  else
	    count_bound_unbound--;
	  if(count_bound_unbound <= 0)
	    break;
	}

      if(num > 0 && count_bound_unbound <= 0)
	weakly_bound_limit_local = bnd_energy[i];
      else
	weakly_bound_limit_local = -1.0e30;

      MPI_Allreduce(&weakly_bound_limit_local, &weakly_bound_limit, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      for(i = 0, unbound = 0; i < num; i++)
	{
	  p= d[i].index;

	  if(P[p].v.DM_BindingEnergy > 0 && P[p].v.DM_BindingEnergy > energy_limit)
	    {
	      unbound++;

	      d[i] = d[num - 1];
	      num--;
	      i--;
	    }
	}

      myfree(offset);
      myfree(nbu_count);
      myfree(npart);
      myfree(bnd_energy);

      MPI_Allreduce(&unbound, &totunbound, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&num, &numleft, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      t1 = second();

      if(ThisTask == 0)
	printf("iter=%d phaseflag=%d unbound=%d numleft=%d  (took %g sec)\n", iter, phaseflag, totunbound,
	       numleft, timediff(t0, t1));

      if(phaseflag == 0)
	{
	  if(totunbound > 0)
	    phaseflag = 1;
	}
      else
	{
	  if(totunbound == 0)
	    {
	      phaseflag = 0;	/* this will make us repeat everything once more for all particles */
	      totunbound = 1;
	    }
	}

      iter++;
    }
  while(totunbound > 0 && numleft >= All.DesLinkNgb);

  return num;
}




void subfind_col_determine_sub_halo_properties(struct unbind_data *d, int num, double *totmass,
					       double *pos, double *vel, double *cm, double *veldisp,
					       double *vmax, double *vmaxrad, double *spin,
					       MyIDType * mostboundid, double *halfmassrad)
{
  int i, j, part_index, *npart, numtot, mbid, nhalf, offset;
  double s[3], sloc[3], v[3], vloc[3], max, vel_to_phys, H_of_a;
  double lx, ly, lz, dv[3], dx[3], disp;
  double loclx, locly, loclz, locdisp;
  double boxhalf, boxsize, ddxx;
  sort_r2list *loc_rr_list;
  int minindex, mincpu;
  double mass, maxrad, massloc, *masslist;
  MyFloat minpot, *potlist;


  boxsize = All.BoxSize;
  boxhalf = 0.5 * boxsize;

  vel_to_phys = 1.0 / All.Time;

  H_of_a = hubble_function(All.Time);

  potlist = (MyFloat *) mymalloc(NTask * sizeof(MyFloat));

  for(i = 0, minindex = -1, minpot = 1.0e30; i < num; i++)
    {
      if(P[d[i].index].u.DM_Potential < minpot || minindex == -1)
	{
	  minpot = P[d[i].index].u.DM_Potential;
	  minindex = d[i].index;
	}
    }

  MPI_Allgather(&minpot, sizeof(MyFloat), MPI_BYTE, potlist, sizeof(MyFloat), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0, mincpu = -1, minpot = 1.0e30; i < NTask; i++)
    if(potlist[i] < minpot)
      {
	mincpu = i;
	minpot = potlist[i];
      }

  if(mincpu < 0)
    {
      printf("ta=%d num=%d\n", ThisTask, num);
      endrun(121);
    }

  if(ThisTask == mincpu)
    {
      for(j = 0; j < 3; j++)
	s[j] = P[minindex].Pos[j];
    }

  MPI_Bcast(&s[0], 3, MPI_DOUBLE, mincpu, MPI_COMM_WORLD);

  /* s[] now holds the position of minimum potential */
  /* we take that as the center */
  for(j = 0; j < 3; j++)
    pos[j] = s[j];


  /* the ID of the most bound particle, we take the minimum binding energy */
  for(i = 0, minindex = -1, minpot = 1.0e30; i < num; i++)
    {
      if(P[d[i].index].v.DM_BindingEnergy < minpot || minindex == -1)
	{
	  minpot = P[d[i].index].v.DM_BindingEnergy;
	  minindex = d[i].index;
	}
    }

  MPI_Allgather(&minpot, sizeof(MyFloat), MPI_BYTE, potlist, sizeof(MyFloat), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0, minpot = 1.0e30; i < NTask; i++)
    if(potlist[i] < minpot)
      {
	mincpu = i;
	minpot = potlist[i];
      }

  if(ThisTask == mincpu)
    mbid = P[minindex].ID;

  MPI_Bcast(&mbid, sizeof(mbid), MPI_BYTE, mincpu, MPI_COMM_WORLD);

  myfree(potlist);

  *mostboundid = mbid;

  /* let's get bulk velocity and the center-of-mass */

  for(j = 0; j < 3; j++)
    sloc[j] = vloc[j] = 0;

  for(i = 0, massloc = 0; i < num; i++)
    {
      part_index = d[i].index;

      for(j = 0; j < 3; j++)
	{
#ifdef PERIODIC
	  ddxx = NEAREST(P[part_index].Pos[j] - pos[j]);
#else
	  ddxx = P[part_index].Pos[j] - pos[j];
#endif
	  sloc[j] += P[part_index].Mass * ddxx;
	  vloc[j] += P[part_index].Mass * P[part_index].Vel[j];
	}
      massloc += P[part_index].Mass;
    }

  MPI_Allreduce(sloc, s, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(vloc, v, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&massloc, &mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  *totmass = mass;

  for(j = 0; j < 3; j++)
    {
      s[j] /= mass;		/* center of mass */
      v[j] /= mass;

      vel[j] = vel_to_phys * v[j];

      s[j] += pos[j];

#ifdef PERIODIC
      while(s[j] < 0)
	s[j] += boxsize;
      while(s[j] >= boxsize)
	s[j] -= boxsize;
#endif

      cm[j] = s[j];
    }


  locdisp = loclx = locly = loclz = 0;

  loc_rr_list = mymalloc(sizeof(sort_r2list) * num);

  for(i = 0, massloc = 0; i < num; i++)
    {
      part_index = d[i].index;

      loc_rr_list[i].r = 0;
      loc_rr_list[i].mass = P[part_index].Mass;

      for(j = 0; j < 3; j++)
	{
#ifdef PERIODIC
	  ddxx = NEAREST(P[part_index].Pos[j] - s[j]);
#else
	  ddxx = P[part_index].Pos[j] - s[j];
#endif
	  dx[j] = All.Time * ddxx;
	  dv[j] = vel_to_phys * (P[part_index].Vel[j] - v[j]);
	  dv[j] += H_of_a * dx[j];

	  locdisp += P[part_index].Mass * dv[j] * dv[j];
	  /* for rotation curve computation, take minimum of potential as center */
#ifdef PERIODIC
	  ddxx = NEAREST(P[part_index].Pos[j] - pos[j]);
#else
	  ddxx = P[part_index].Pos[j] - pos[j];
#endif
	  ddxx = All.Time * ddxx;
	  loc_rr_list[i].r += ddxx * ddxx;
	}

      loclx += P[part_index].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
      locly += P[part_index].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
      loclz += P[part_index].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);

      loc_rr_list[i].r = sqrt(loc_rr_list[i].r);
    }

  MPI_Allreduce(&loclx, &lx, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&locly, &ly, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&loclz, &lz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&locdisp, &disp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  *veldisp = sqrt(disp / (3 * mass));	/* convert to 1d velocity dispersion */

  spin[0] = lx / mass;
  spin[1] = ly / mass;
  spin[2] = lz / mass;


  npart = (int *) mymalloc(NTask * sizeof(int));

  MPI_Allgather(&num, 1, MPI_INT, npart, 1, MPI_INT, MPI_COMM_WORLD);

  for(i = 0, numtot = 0; i < NTask; i++)
    numtot += npart[i];

  parallel_sort(loc_rr_list, num, sizeof(sort_r2list), subfind_compare_dist_rotcurve);

  nhalf = numtot / 2;
  mincpu = 0;

  while(nhalf >= npart[mincpu])
    {
      nhalf -= npart[mincpu];
      mincpu++;
    }

  if(ThisTask == mincpu)
    *halfmassrad = loc_rr_list[nhalf].r;

  MPI_Bcast(halfmassrad, sizeof(double), MPI_BYTE, mincpu, MPI_COMM_WORLD);


  /* compute cumulative mass */

  masslist = (double *) mymalloc(NTask * sizeof(double));

  for(i = 0, massloc = 0; i < num; i++)
    massloc += loc_rr_list[i].mass;

  MPI_Allgather(&massloc, 1, MPI_DOUBLE, masslist, 1, MPI_DOUBLE, MPI_COMM_WORLD);

  for(i = 1; i < NTask; i++)
    masslist[i] += masslist[i - 1];

  for(i = 1; i < num; i++)
    loc_rr_list[i].mass += loc_rr_list[i - 1].mass;

  if(ThisTask > 0)
    for(i = 0; i < num; i++)
      loc_rr_list[i].mass += masslist[ThisTask - 1];

  for(i = 0, offset = 0; i < ThisTask; i++)
    offset += npart[i];

  for(i = num - 1, max = 0, maxrad = 0; i + offset > 5 && i >=0 ; i--)
    if(loc_rr_list[i].mass / loc_rr_list[i].r > max)
      {
	max = loc_rr_list[i].mass / loc_rr_list[i].r;
	maxrad = loc_rr_list[i].r;
      }

  MPI_Allreduce(&max, vmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  if(max < *vmax)
    maxrad = 0;

  *vmax = sqrt(All.G * (*vmax));

  MPI_Allreduce(&maxrad, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  *vmaxrad = max;


  myfree(masslist);
  myfree(npart);
  myfree(loc_rr_list);
}



void subfind_poll_for_requests(void)
{
  int index, nsub, source, tag, ibuf[2], target, submark;
  long long head, next, rank, buf[4];
  MPI_Status status;

  do
    {
      MPI_Recv(&index, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

      source = status.MPI_SOURCE;
      tag = status.MPI_TAG;
      /*      MPI_Get_count(&status, MPI_BYTE, &count); */

      switch (tag)
	{
	case TAG_GET_NGB_COUNT:
	  MPI_Send(&NgbLoc[index].count, 1, MPI_INT, source, TAG_GET_NGB_COUNT, MPI_COMM_WORLD);
	  MPI_Send(&NgbLoc[index].index[0], 2 * sizeof(long long), MPI_BYTE, source, TAG_GET_NGB_INDICES,
		   MPI_COMM_WORLD);
	  break;
	case TAG_SET_ALL:
	  MPI_Recv(buf, 4 * sizeof(long long), MPI_BYTE, source, TAG_SET_ALL_DATA, MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);
	  Head[index] = buf[0];
	  Tail[index] = buf[1];
	  Len[index] = buf[2];
	  Next[index] = buf[3];
	  break;
	case TAG_GET_TAILANDLEN:
	  buf[0] = Tail[index];
	  buf[1] = Len[index];
	  MPI_Send(buf, 2 * sizeof(long long), MPI_BYTE, source, TAG_GET_TAILANDLEN_DATA, MPI_COMM_WORLD);
	  break;
	case TAG_SET_TAILANDLEN:
	  MPI_Recv(buf, 2 * sizeof(long long), MPI_BYTE, source, TAG_SET_TAILANDLEN_DATA, MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);
	  Tail[index] = buf[0];
	  Len[index] = buf[1];
	  break;
	case TAG_SET_HEADANDNEXT:
	  MPI_Recv(buf, 2 * sizeof(long long), MPI_BYTE, source, TAG_SET_HEADANDNEXT_DATA, MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);
	  Head[index] = buf[0];
	  Next[index] = buf[1];
	  break;
	case TAG_SET_NEXT:
	  MPI_Recv(buf, 1 * sizeof(long long), MPI_BYTE, source, TAG_SET_NEXT_DATA, MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);
	  Next[index] = buf[0];
	  break;
	case TAG_SETHEADGETNEXT:
	  MPI_Sendrecv(&Next[index], 1 * sizeof(long long), MPI_BYTE, source, TAG_SETHEADGETNEXT_DATA,
		       &head, 1 * sizeof(long long), MPI_BYTE, source, TAG_SETHEADGETNEXT_DATA,
		       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  Head[index] = head;
	  break;
	case TAG_GET_NEXT:
	  MPI_Send(&Next[index], 1 * sizeof(long long), MPI_BYTE, source, TAG_GET_NEXT_DATA, MPI_COMM_WORLD);
	  break;
	case TAG_GET_HEAD:
	  MPI_Send(&Head[index], 1 * sizeof(long long), MPI_BYTE, source, TAG_GET_HEAD_DATA, MPI_COMM_WORLD);
	  break;
	case TAG_ADD_PARTICLE:
	  if(Tail[index] < 0)	/* consider only particles not already in substructures */
	    {
	      ud[LocalLen].index = index;
	      if(index >= NumPartGroup)
		{
		  printf("What: index=%d NumPartGroup=%d\n", index, NumPartGroup);
		  endrun(199);
		}
	      LocalLen++;
	    }
	  break;
	case TAG_MARK_PARTICLE:
	  MPI_Recv(ibuf, 2, MPI_INT, source, TAG_MARK_PARTICLE_DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  target = ibuf[0];
	  submark = ibuf[1];

	  if(P[index].submark != HIGHBIT)
	    {
	      printf("TasK=%d i=%d P[i].submark=%d?\n", ThisTask, index, P[index].submark);
	      endrun(132);
	    }

	  P[index].targettask = target;
	  P[index].submark = submark;
	  break;
	case TAG_ADDBOUND:
	  MPI_Recv(&nsub, 1, MPI_INT, source, TAG_ADDBOUND_DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  if(Tail[index] == nsub)	/* consider only particles in this substructure */
	    {
	      ud[LocalLen].index = index;
	      LocalLen++;
	    }
	  break;
	case TAG_SETRANK:
	  MPI_Recv(&rank, sizeof(long long), MPI_BYTE, source, TAG_SETRANK_IN, MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);
	  do
	    {
	      Len[index] = rank++;
	      next = Next[index];
	      if(next < 0)
		break;
	      index = (next & MASK);
	    }
	  while((next >> 32) == ThisTask);
	  buf[0] = next;
	  buf[1] = rank;
	  MPI_Send(buf, 2 * sizeof(long long), MPI_BYTE, source, TAG_SETRANK_OUT, MPI_COMM_WORLD);
	  break;
	case TAG_GET_RANK:
	  rank = Len[index];
	  MPI_Send(&rank, 1 * sizeof(long long), MPI_BYTE, source, TAG_GET_RANK_DATA, MPI_COMM_WORLD);
	  break;

	case TAG_POLLING_DONE:
	  break;

	default:
	  endrun(1213);
	  break;
	}

    }
  while(tag != TAG_POLLING_DONE);



}

long long subfind_distlinklist_setrank_and_get_next(long long index, long long *rank)
{
  int task, i;
  long long next;
  long long buf[2];

  task = (index >> 32);
  i = (index & MASK);

  if(ThisTask == task)
    {
      Len[i] = *rank;
      *rank = *rank + 1;
      next = Next[i];
    }
  else
    {
      MPI_Send(&i, 1, MPI_INT, task, TAG_SETRANK, MPI_COMM_WORLD);
      MPI_Send(rank, sizeof(long long), MPI_BYTE, task, TAG_SETRANK_IN, MPI_COMM_WORLD);
      MPI_Recv(buf, 2 * sizeof(long long), MPI_BYTE, task, TAG_SETRANK_OUT, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      next = buf[0];
      *rank = buf[1];
    }
  return next;
}


long long subfind_distlinklist_set_head_get_next(long long index, long long head)
{
  int task, i;
  long long next;

  task = (index >> 32);
  i = (index & MASK);

  if(ThisTask == task)
    {
      Head[i] = head;
      next = Next[i];
    }
  else
    {
      MPI_Send(&i, 1, MPI_INT, task, TAG_SETHEADGETNEXT, MPI_COMM_WORLD);
      MPI_Sendrecv(&head, 1 * sizeof(long long), MPI_BYTE, task, TAG_SETHEADGETNEXT_DATA,
		   &next, 1 * sizeof(long long), MPI_BYTE, task, TAG_SETHEADGETNEXT_DATA, MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);

    }

  return next;
}




void subfind_distlinklist_set_next(long long index, long long next)
{
  int task, i;
  long long buf[1];

  task = (index >> 32);
  i = (index & MASK);

  if(ThisTask == task)
    {
      Next[i] = next;
    }
  else
    {
      MPI_Send(&i, 1, MPI_INT, task, TAG_SET_NEXT, MPI_COMM_WORLD);
      buf[0] = next;
      MPI_Send(buf, 1 * sizeof(long long), MPI_BYTE, task, TAG_SET_NEXT_DATA, MPI_COMM_WORLD);
    }
}

void subfind_distlinklist_add_particle(long long index)
{
  int task, i;

  task = (index >> 32);
  i = (index & MASK);

  if(ThisTask == task)
    {
      if(Tail[i] < 0)		/* consider only particles not already in substructures */
	{
	  ud[LocalLen].index = i;
	  if(i >= NumPartGroup)
	    {
	      printf("What: index=%d NumPartGroup=%d\n", i, NumPartGroup);
	      endrun(299);
	    }

	  LocalLen++;
	}
    }
  else
    {
      MPI_Send(&i, 1, MPI_INT, task, TAG_ADD_PARTICLE, MPI_COMM_WORLD);
    }
}

void subfind_distlinklist_mark_particle(long long index, int target, int submark)
{
  int task, i, ibuf[2];

  task = (index >> 32);
  i = (index & MASK);

  if(ThisTask == task)
    {
      if(P[i].submark != HIGHBIT)
	{
	  printf("Tas=%d i=%d P[i].submark=%d?\n", ThisTask, i, P[i].submark);
	  endrun(131);
	}

      P[i].targettask = target;
      P[i].submark = submark;
    }
  else
    {
      MPI_Send(&i, 1, MPI_INT, task, TAG_MARK_PARTICLE, MPI_COMM_WORLD);
      ibuf[0] = target;
      ibuf[1] = submark;
      MPI_Send(ibuf, 2, MPI_INT, task, TAG_MARK_PARTICLE_DATA, MPI_COMM_WORLD);
    }
}


void subfind_distlinklist_add_bound_particles(long long index, int nsub)
{
  int task, i;

  task = (index >> 32);
  i = (index & MASK);

  if(ThisTask == task)
    {
      if(Tail[i] == nsub)	/* consider only particles not already in substructures */
	{
	  ud[LocalLen].index = i;
	  LocalLen++;
	}
    }
  else
    {
      MPI_Send(&i, 1, MPI_INT, task, TAG_ADDBOUND, MPI_COMM_WORLD);
      MPI_Send(&nsub, 1, MPI_INT, task, TAG_ADDBOUND_DATA, MPI_COMM_WORLD);
    }
}


long long subfind_distlinklist_get_next(long long index)
{
  int task, i;
  long long next;

  task = (index >> 32);
  i = (index & MASK);

  if(ThisTask == task)
    {
      next = Next[i];
    }
  else
    {
      MPI_Send(&i, 1, MPI_INT, task, TAG_GET_NEXT, MPI_COMM_WORLD);
      MPI_Recv(&next, 1 * sizeof(long long), MPI_BYTE, task, TAG_GET_NEXT_DATA, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
    }

  return next;
}

long long subfind_distlinklist_get_rank(long long index)
{
  int task, i;
  long long rank;

  task = (index >> 32);
  i = (index & MASK);

  if(ThisTask == task)
    {
      rank = Len[i];
    }
  else
    {
      MPI_Send(&i, 1, MPI_INT, task, TAG_GET_RANK, MPI_COMM_WORLD);
      MPI_Recv(&rank, 1 * sizeof(long long), MPI_BYTE, task, TAG_GET_RANK_DATA, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
    }

  return rank;
}



long long subfind_distlinklist_get_head(long long index)
{
  int task, i;
  long long head;

  task = (index >> 32);
  i = (index & MASK);

  if(ThisTask == task)
    {
      head = Head[i];
    }
  else
    {
      MPI_Send(&i, 1, MPI_INT, task, TAG_GET_HEAD, MPI_COMM_WORLD);
      MPI_Recv(&head, 1 * sizeof(long long), MPI_BYTE, task, TAG_GET_HEAD_DATA, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
    }

  return head;
}



void subfind_distlinklist_set_headandnext(long long index, long long head, long long next)
{
  int task, i;
  long long buf[2];

  task = (index >> 32);
  i = (index & MASK);

  if(ThisTask == task)
    {
      Head[i] = head;
      Next[i] = next;
    }
  else
    {
      MPI_Send(&i, 1, MPI_INT, task, TAG_SET_HEADANDNEXT, MPI_COMM_WORLD);
      buf[0] = head;
      buf[1] = next;
      MPI_Send(buf, 2 * sizeof(long long), MPI_BYTE, task, TAG_SET_HEADANDNEXT_DATA, MPI_COMM_WORLD);
    }
}



void subfind_distlinklist_set_tailandlen(long long index, long long tail, int len)
{
  int task, i;
  long long buf[2];

  task = (index >> 32);
  i = (index & MASK);

  if(ThisTask == task)
    {
      Tail[i] = tail;
      Len[i] = len;
    }
  else
    {
      MPI_Send(&i, 1, MPI_INT, task, TAG_SET_TAILANDLEN, MPI_COMM_WORLD);
      buf[0] = tail;
      buf[1] = len;
      MPI_Send(buf, 2 * sizeof(long long), MPI_BYTE, task, TAG_SET_TAILANDLEN_DATA, MPI_COMM_WORLD);
    }
}




void subfind_distlinklist_get_tailandlen(long long index, long long *tail, int *len)
{
  int task, i;
  long long buf[2];

  task = (index >> 32);
  i = (index & MASK);

  if(ThisTask == task)
    {
      *tail = Tail[i];
      *len = Len[i];
    }
  else
    {
      MPI_Send(&i, 1, MPI_INT, task, TAG_GET_TAILANDLEN, MPI_COMM_WORLD);
      MPI_Recv(buf, 2 * sizeof(long long), MPI_BYTE, task, TAG_GET_TAILANDLEN_DATA, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      *tail = buf[0];
      *len = buf[1];
    }
}


void subfind_distlinklist_set_all(long long index, long long head, long long tail, int len, long long next)
{
  int task, i;
  long long buf[4];

  task = (index >> 32);
  i = (index & MASK);

  if(ThisTask == task)
    {
      Head[i] = head;
      Tail[i] = tail;
      Len[i] = len;
      Next[i] = next;
    }
  else
    {
      MPI_Send(&i, 1, MPI_INT, task, TAG_SET_ALL, MPI_COMM_WORLD);
      buf[0] = head;
      buf[1] = tail;
      buf[2] = len;
      buf[3] = next;
      MPI_Send(buf, 4 * sizeof(long long), MPI_BYTE, task, TAG_SET_ALL_DATA, MPI_COMM_WORLD);
    }
}

int subfind_distlinklist_get_ngb_count(long long index, long long *ngb_index1, long long *ngb_index2)
{
  int ngbcount, task, i;
  long long buf[2];

  task = (index >> 32);
  i = (index & MASK);


  if(ThisTask == task)
    {
      ngbcount = NgbLoc[i].count;
      *ngb_index1 = NgbLoc[i].index[0];
      *ngb_index2 = NgbLoc[i].index[1];
    }
  else
    {
      MPI_Send(&i, 1, MPI_INT, task, TAG_GET_NGB_COUNT, MPI_COMM_WORLD);
      MPI_Recv(&ngbcount, 1, MPI_INT, task, TAG_GET_NGB_COUNT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(buf, 2 * sizeof(long long), MPI_BYTE, task, TAG_GET_NGB_INDICES, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      *ngb_index1 = buf[0];
      *ngb_index2 = buf[1];
    }

  return ngbcount;
}




int subfind_compare_P_GrNrGrNr(const void *a, const void *b)
{
  if(abs(((struct particle_data *) a)->GrNr - GrNr) < abs(((struct particle_data *) b)->GrNr - GrNr))
    return -1;

  if(abs(((struct particle_data *) a)->GrNr - GrNr) > abs(((struct particle_data *) b)->GrNr - GrNr))
    return +1;

  if(((struct particle_data *) a)->ID < ((struct particle_data *) b)->ID)
    return -1;

  if(((struct particle_data *) a)->ID > ((struct particle_data *) b)->ID)
    return +1;

  return 0;
}

int subfind_compare_P_submark(const void *a, const void *b)
{
  if(((struct particle_data *) a)->submark < ((struct particle_data *) b)->submark)
    return -1;

  if(((struct particle_data *) a)->submark > ((struct particle_data *) b)->submark)
    return +1;

  return 0;
}


int subfind_compare_candidates_subnr(const void *a, const void *b)
{
  if(((struct cand_dat *) a)->subnr < ((struct cand_dat *) b)->subnr)
    return -1;

  if(((struct cand_dat *) a)->subnr > ((struct cand_dat *) b)->subnr)
    return +1;

  return 0;
}

int subfind_compare_candidates_nsubs(const void *a, const void *b)
{
  if(((struct cand_dat *) a)->nsub < ((struct cand_dat *) b)->nsub)
    return -1;

  if(((struct cand_dat *) a)->nsub > ((struct cand_dat *) b)->nsub)
    return +1;

  return 0;
}

int subfind_compare_candidates_boundlength(const void *a, const void *b)
{
  if(((struct cand_dat *) a)->bound_length > ((struct cand_dat *) b)->bound_length)
    return -1;

  if(((struct cand_dat *) a)->bound_length < ((struct cand_dat *) b)->bound_length)
    return +1;

  if(((struct cand_dat *) a)->rank < ((struct cand_dat *) b)->rank)
    return -1;

  if(((struct cand_dat *) a)->rank > ((struct cand_dat *) b)->rank)
    return +1;

  return 0;
}

int subfind_compare_candidates_rank(const void *a, const void *b)
{
  if(((struct cand_dat *) a)->rank < ((struct cand_dat *) b)->rank)
    return -1;

  if(((struct cand_dat *) a)->rank > ((struct cand_dat *) b)->rank)
    return +1;

  if(((struct cand_dat *) a)->len > ((struct cand_dat *) b)->len)
    return -1;

  if(((struct cand_dat *) a)->len < ((struct cand_dat *) b)->len)
    return +1;

  return 0;
}




int subfind_compare_dist_rotcurve(const void *a, const void *b)
{
  if(((sort_r2list *) a)->r < ((sort_r2list *) b)->r)
    return -1;

  if(((sort_r2list *) a)->r > ((sort_r2list *) b)->r)
    return +1;

  return 0;
}

int subfind_compare_binding_energy(const void *a, const void *b)
{
  if(*((double *) a) > *((double *) b))
    return -1;

  if(*((double *) a) < *((double *) b))
    return +1;

  return 0;
}


int subfind_compare_densities(const void *a, const void *b)	/* largest density first */
{
  if(((struct sort_density_data *) a)->density > (((struct sort_density_data *) b)->density))
    return -1;

  if(((struct sort_density_data *) a)->density < (((struct sort_density_data *) b)->density))
    return +1;

  return 0;
}



#endif
