#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

/*! \file bsmooth.c
 *  \brief SPH smoothing of the B fiold
 */

#if defined(MAGNETIC) && defined(BSMOOTH)


void bsmooth(void)
{
  long long ntot, ntotleft;
  int *noffset, *nbuffer, *nsend, *nsend_local, *ndonelist, *numlist;
  int i, j, n;
  int ndone;
  int maxfill, source;
  int level, ngrp, sendTask, recvTask;
  int place, nexport;
  double tstart, tend;
  double maxt, sumt, sumcomm;
  double timecomp = 0, timeimbalance = 0, timecommsumm = 0;
  double timengb, sumtimengb;
  MPI_Status status;

  int Smooth_Flag = 0;
  double dB;

  if(ThisTask == 0)
    printf("Flag_FullStep = %d, Main TimestepCounts = %d\n", Flag_FullStep, All.MainTimestepCounts);
  if(Flag_FullStep == 1)
    {
      if((All.MainTimestepCounts % All.BSmoothInt == 0) && (All.BSmoothInt >= 0))
	{
	  Smooth_Flag = 1;
	  if(ThisTask == 0)
	    printf("Smoothing B %d, %f\n", All.BSmoothInt, All.BSmoothFrac);
	}
      All.MainTimestepCounts++;
    }

  noffset = mymalloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = mymalloc(sizeof(int) * NTask);
  nsend_local = mymalloc(sizeof(int) * NTask);
  nsend = mymalloc(sizeof(int) * NTask * NTask);
  ndonelist = mymalloc(sizeof(int) * NTask);

  for(n = 0, NumSphUpdate = 0; n < N_gas; n++)
    {
#ifdef SFR
      if(P[n].Type == 0)
#endif
	{
	  SphP[n].Left = SphP[n].Right = 0;

	  if(P[n].Ti_endstep == All.Ti_Current)
	    NumSphUpdate++;
	}
    }

  numlist = mymalloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  myfree(numlist);

  i = 0;			/* beginn with this index */
  ntotleft = ntot;		/* particles left for all tasks together */

  while(ntotleft > 0)
    {
      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

      /* do local particles and prepare export list */
      tstart = second();
      for(nexport = 0, ndone = 0; i < N_gas && nexport < All.BunchSizeDensity - NTask; i++)
#ifdef SFR
	if(P[i].Type == 0)
#endif
	  if(P[i].Ti_endstep == All.Ti_Current)
	    {
	      ndone++;

	      for(j = 0; j < NTask; j++)
		Exportflag[j] = 0;

	      bsmooth_evaluate(i, 0);

	      for(j = 0; j < NTask; j++)
		{
		  if(Exportflag[j])
		    {
		      DensDataIn[nexport].Pos[0] = P[i].Pos[0];
		      DensDataIn[nexport].Pos[1] = P[i].Pos[1];
		      DensDataIn[nexport].Pos[2] = P[i].Pos[2];
		      DensDataIn[nexport].Hsml = PPP[i].Hsml;
		      DensDataIn[nexport].Index = i;
		      DensDataIn[nexport].Task = j;
		      nexport++;
		      nsend_local[j]++;
		    }
		}
	    }
      tend = second();
      timecomp += timediff(tstart, tend);

      qsort(DensDataIn, nexport, sizeof(struct densdata_in), dens_compare_key);

      for(j = 1, noffset[0] = 0; j < NTask; j++)
	noffset[j] = noffset[j - 1] + nsend_local[j - 1];

      tstart = second();

      MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

      tend = second();
      timeimbalance += timediff(tstart, tend);


      /* now do the particles that need to be exported */

      for(level = 1; level < (1 << PTask); level++)
	{
	  tstart = second();
	  for(j = 0; j < NTask; j++)
	    nbuffer[j] = 0;
	  for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
	      maxfill = 0;
	      for(j = 0; j < NTask; j++)
		{
		  if((j ^ ngrp) < NTask)
		    if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		      maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
	      if(maxfill >= All.BunchSizeForce)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&DensDataIn[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
				   recvTask, TAG_BSMTH_A,
				   &DensDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct densdata_in),
				   MPI_BYTE, recvTask, TAG_BSMTH_A, MPI_COMM_WORLD, &status);
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }
	  tend = second();
	  timecommsumm += timediff(tstart, tend);


	  tstart = second();
	  for(j = 0; j < nbuffer[ThisTask]; j++)
	    {
	      bsmooth_evaluate(j, 1);
	    }
	  tend = second();
	  timecomp += timediff(tstart, tend);

	  /* do a block to measure imbalance */
	  tstart = second();
	  MPI_Barrier(MPI_COMM_WORLD);
	  tend = second();
	  timeimbalance += timediff(tstart, tend);

	  /* get the result */
	  tstart = second();
	  for(j = 0; j < NTask; j++)
	    nbuffer[j] = 0;
	  for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
	      maxfill = 0;
	      for(j = 0; j < NTask; j++)
		{
		  if((j ^ ngrp) < NTask)
		    if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		      maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
	      if(maxfill >= All.BunchSizeForce)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&DensDataResult[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct densdata_out),
				   MPI_BYTE, recvTask, TAG_BSMTH_B,
				   &DensDataPartialResult[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct densdata_out),
				   MPI_BYTE, recvTask, TAG_BSMTH_B, MPI_COMM_WORLD, &status);

		      /* add the result to the particles */
		      for(j = 0; j < nsend_local[recvTask]; j++)
			{
			  source = j + noffset[recvTask];
			  place = DensDataIn[source].Index;

			  SphP[place].BSmooth[0] += DensDataPartialResult[source].BSmooth[0];
			  SphP[place].BSmooth[1] += DensDataPartialResult[source].BSmooth[1];
			  SphP[place].BSmooth[2] += DensDataPartialResult[source].BSmooth[2];
			  SphP[place].DensityNorm += DensDataPartialResult[source].DensityNorm;
			}
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }
	  tend = second();
	  timecommsumm += timediff(tstart, tend);

	  level = ngrp - 1;
	}

      MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);
      for(j = 0; j < NTask; j++)
	ntotleft -= ndonelist[j];
    }



  /* do final operations on results */
  tstart = second();
  for(i = 0; i < N_gas; i++)
    {
#ifdef SFR
      if(P[i].Type == 0)
#endif
	if(P[i].Ti_endstep == All.Ti_Current)
	  {
	    SphP[i].BSmooth[0] /= SphP[i].DensityNorm;
	    SphP[i].BSmooth[1] /= SphP[i].DensityNorm;
	    SphP[i].BSmooth[2] /= SphP[i].DensityNorm;
	    if(Smooth_Flag == 1)
	      {
		dB = All.BSmoothFrac * (SphP[i].BSmooth[0] - SphP[i].BPred[0]);
		SphP[i].BPred[0] += dB;
		SphP[i].B[0] += dB;
		dB = All.BSmoothFrac * (SphP[i].BSmooth[1] - SphP[i].BPred[1]);
		SphP[i].BPred[1] += dB;
		SphP[i].B[1] += dB;
		dB = All.BSmoothFrac * (SphP[i].BSmooth[2] - SphP[i].BPred[2]);
		SphP[i].BPred[2] += dB;
		SphP[i].B[2] += dB;
	      }
	  }
    }
  tend = second();
  timecomp += timediff(tstart, tend);

  myfree(ndonelist);
  myfree(nsend);
  myfree(nsend_local);
  myfree(nbuffer);
  myfree(noffset);


  /* collect some timing information */

  timengb = 0;

  MPI_Reduce(&timengb, &sumtimengb, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timecomp, &sumt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timecommsumm, &sumcomm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  timeimbalance += timecomp + timecommsumm;
  MPI_Reduce(&timeimbalance, &maxt, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      All.CPU_HydCompWalk += sumt / NTask;
      All.CPU_HydCommSumm += sumcomm / NTask;
      All.CPU_HydImbalance += maxt - (sumt + sumcomm) / NTask;
      All.CPU_EnsureNgb += sumtimengb / NTask;
    }
}



void bsmooth_evaluate(int target, int mode)
{
  int j, n;
  int startnode, numngb_inbox;
  double h, h2, hinv, hinv3;
  double wk;
  double dx, dy, dz, r, r2, u, mass_j;
  MyFloat pos[3];

  MyFloat BSmooth[3], DensityNorm, mwk_d;

  BSmooth[0] = BSmooth[1] = BSmooth[2] = DensityNorm = 0;

  if(mode == 0)
    {
      pos[0] = P[target].Pos[0];
      pos[1] = P[target].Pos[1];
      pos[2] = P[target].Pos[2];
      h = PPP[target].Hsml;
    }
  else
    {
      pos[0] = DensDataGet[target].Pos[0];
      pos[1] = DensDataGet[target].Pos[1];
      pos[2] = DensDataGet[target].Pos[2];
      h = DensDataGet[target].Hsml;
    }


  h2 = h * h;
  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif

  startnode = All.MaxPart;
  do
    {
      numngb_inbox = ngb_treefind_variable(&pos[0], h, &startnode);

      for(n = 0; n < numngb_inbox; n++)
	{
	  j = Ngblist[n];
	  /* MHD sould be solved as well for wind particles ! */

	  dx = pos[0] - P[j].Pos[0];
	  dy = pos[1] - P[j].Pos[1];
	  dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC			/*  now find the closest image in the given box size  */
	  if(dx > boxHalf_X)
	    dx -= boxSize_X;
	  if(dx < -boxHalf_X)
	    dx += boxSize_X;
	  if(dy > boxHalf_Y)
	    dy -= boxSize_Y;
	  if(dy < -boxHalf_Y)
	    dy += boxSize_Y;
	  if(dz > boxHalf_Z)
	    dz -= boxSize_Z;
	  if(dz < -boxHalf_Z)
	    dz += boxSize_Z;
#endif
	  r2 = dx * dx + dy * dy + dz * dz;

	  if(r2 < h2)
	    {
	      r = sqrt(r2);

	      u = r * hinv;

	      if(u < 0.5)
		{
		  wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
		}
	      else
		{
		  wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
		}

	      mass_j = P[j].Mass;

	      mwk_d = mass_j * wk / SphP[j].d.Density;

	      BSmooth[0] += SphP[j].BPred[0] * mwk_d;
	      BSmooth[1] += SphP[j].BPred[1] * mwk_d;
	      BSmooth[2] += SphP[j].BPred[2] * mwk_d;
	      DensityNorm += mwk_d;
	    }
	}
    }
  while(startnode >= 0);


  if(mode == 0)
    {
      SphP[target].BSmooth[0] = BSmooth[0];
      SphP[target].BSmooth[1] = BSmooth[1];
      SphP[target].BSmooth[2] = BSmooth[2];
      SphP[target].DensityNorm = DensityNorm;
    }
  else
    {
      DensDataResult[target].BSmooth[0] = BSmooth[0];
      DensDataResult[target].BSmooth[1] = BSmooth[1];
      DensDataResult[target].BSmooth[2] = BSmooth[2];
      DensDataResult[target].DensityNorm = DensityNorm;
    }
}






#endif
