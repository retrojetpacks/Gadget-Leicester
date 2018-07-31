#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

#if ( defined(CONDUCTION) || defined(CR_DIFFUSION) || defined(SMOOTH_PHI))



/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */
static struct smoothdata_in
{
  MyDouble Pos[3];
  MyFloat Hsml;
  int NodeList[NODELISTLENGTH];
}
 *SmoothDataIn, *SmoothDataGet;


static struct smoothdata_out
{
#ifdef SMOOTH_PHI
  MyFloat SmoothPhi;
#endif

#if defined(CONDUCTION) || defined(BLACK_HOLES) || defined (DUST)
  MyLongDouble SmoothedEntr;
#ifdef CONDUCTION_SATURATION
  MyFloat GradEntr[3];
#endif
#endif

#ifdef CR_DIFFUSION
  MyFloat CR_SmoothE0;
  MyFloat CR_Smoothn0;
#endif
#ifdef CR_DIFFUSION_GREEN
  MyFloat CR_WeightSum;
  MyFloat CR_WeightSum_egy;
#endif
}
 *SmoothDataResult, *SmoothDataOut;





void compute_smoothed_values(void)
{
  long long ntot, ntotleft;
  int *ndone_list, *send_offset, *send_count, *recv_count, *recv_offset, *sendcount_matrix;
  int ngrp, sendTask, recvTask, place, nexport, nimport, maxexport;
  int i, j, n, ndone, dummy;
  MPI_Status status;

  /* Display information message that this step is executed on Task 0 ... */
  if(ThisTask == 0)
    {
      printf("Updating SPH interpolants for:"
#ifdef CONDUCTION
	     " (temperature)"
#endif /* CONDUCTION */
#ifdef CR_DIFFUSION
	     " (CR diffusivity terms)"
#endif /* CR_DIFFUSION */
#ifdef SMOOTH_PHI
	     " (divB cleaning phi)"
#endif /* SMOOTH_PHI */
	     "\n");
    }


  for(n = 0, NumSphUpdate = 0; n < N_gas; n++)
    {
      if(density_isactive(n))
	{
	  if(P[n].Ti_endstep == All.Ti_Current)
	    NumSphUpdate++;
	}
    }

  sumup_large_ints(1, &NumSphUpdate, &ntot);


  Exportflag = (int *) mymalloc(NTask * sizeof(int));
  Exportindex = (int *) mymalloc(NTask * sizeof(int));
  Exportnodecount = (int *) mymalloc(NTask * sizeof(int));
  Ngblist = (int *) mymalloc(NumPart * sizeof(int));

  All.BunchSize =
    (int) ((0.5 * All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist)));
  DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));
  maxexport = All.BunchSize - 2 * NTopleaves / NODELISTLENGTH;

  send_count = (int *) mymalloc(sizeof(int) * NTask);
  send_offset = (int *) mymalloc(sizeof(int) * NTask);
  recv_count = (int *) mymalloc(sizeof(int) * NTask);
  recv_offset = (int *) mymalloc(sizeof(int) * NTask);
  ndone_list = (int *) mymalloc(sizeof(int) * NTask);
  sendcount_matrix = (int *) mymalloc(sizeof(int) * NTask * NTask);




  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */

  i = 0;			/* beginn with this index */
  ntotleft = ntot;		/* particles left for all tasks together */

  while(ntotleft > 0)
    {
      for(j = 0; j < NTask; j++)
	{
	  send_count[j] = 0;
	  Exportflag[j] = -1;
	}


      /* do local particles and prepare export list */
      for(nexport = 0, ndone = 0; i < N_gas && nexport < maxexport; i++)
	{
	  if(density_isactive(i))
	    if(P[i].Ti_endstep == All.Ti_Current)
	      {
		ndone++;
		compute_smoothed_evaluate(i, 0, &nexport, send_count);
	      }
	}


      qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);

      MPI_Allgather(send_count, NTask, MPI_INT, sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);

      for(j = 0, nimport = 0, recv_offset[0] = 0, send_offset[0] = 0; j < NTask; j++)
	{
	  recv_count[j] = sendcount_matrix[j * NTask + ThisTask];
	  nimport += recv_count[j];

	  if(j > 0)
	    {
	      send_offset[j] = send_offset[j - 1] + send_count[j - 1];
	      recv_offset[j] = recv_offset[j - 1] + recv_count[j - 1];
	    }
	}

      SmoothDataGet = (struct smoothdata_in *) mymalloc(nimport * sizeof(struct smoothdata_in));
      SmoothDataIn = (struct smoothdata_in *) mymalloc(nexport * sizeof(struct smoothdata_in));

      /* prepare particle data for export */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  SmoothDataIn[j].Pos[0] = P[place].Pos[0];
	  SmoothDataIn[j].Pos[1] = P[place].Pos[1];
	  SmoothDataIn[j].Pos[2] = P[place].Pos[2];
	  SmoothDataIn[j].Hsml = PPP[place].Hsml;

	  memcpy(SmoothDataIn[j].NodeList,
		 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
	}

      /* exchange particle data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      if(send_count[recvTask] > 0 || recv_count[recvTask] > 0)
		{
		  /* get the particles */
		  MPI_Sendrecv(&SmoothDataIn[send_offset[recvTask]],
			       send_count[recvTask] * sizeof(struct smoothdata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A,
			       &SmoothDataGet[recv_offset[recvTask]],
			       recv_count[recvTask] * sizeof(struct smoothdata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
		}
	    }
	}

      myfree(SmoothDataIn);
      SmoothDataResult = (struct smoothdata_out *) mymalloc(nimport * sizeof(struct smoothdata_out));
      SmoothDataOut = (struct smoothdata_out *) mymalloc(nexport * sizeof(struct smoothdata_out));


      /* now do the particles that were sent to us */
      for(j = 0; j < nimport; j++)
	compute_smoothed_evaluate(j, 1, &dummy, &dummy);

      MPI_Allgather(&ndone, 1, MPI_INT, ndone_list, 1, MPI_INT, MPI_COMM_WORLD);


      /* get the result */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;
	  if(recvTask < NTask)
	    {
	      if(send_count[recvTask] > 0 || recv_count[recvTask] > 0)
		{
		  /* send the results */
		  MPI_Sendrecv(&SmoothDataResult[recv_offset[recvTask]],
			       recv_count[recvTask] * sizeof(struct smoothdata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B,
			       &SmoothDataOut[send_offset[recvTask]],
			       send_count[recvTask] * sizeof(struct smoothdata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, &status);
		}
	    }
	}


      /* add the result to the local particles */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  if(P[place].Type == 0)
	    {

#ifdef CONDUCTION
	      SphP[place].SmoothedEntr += SmoothDataOut[j].SmoothedEntr;
#ifdef CONDUCTION_SATURATION
	      SphP[place].GradEntr[0] += SmoothDataOut[j].GradEntr[0];
	      SphP[place].GradEntr[1] += SmoothDataOut[j].GradEntr[1];
	      SphP[place].GradEntr[2] += SmoothDataOut[j].GradEntr[2];
#endif
#endif /* CONDUCTION */

#ifdef CR_DIFFUSION
	      SphP[place].CR_SmoothE0 += SmoothDataOut[j].CR_SmoothE0;
	      SphP[place].CR_Smoothn0 += SmoothDataOut[j].CR_Smoothn0;
#endif /* CR_DIFFUSION */

#ifdef SMOOTH_PHI
	      SphP[place].SmoothPhi += SmoothDataOut[j].SmoothPhi;
#endif /* SMOOTH_PHI */
	    }
	}


      for(j = 0; j < NTask; j++)
	ntotleft -= ndone_list[j];

      myfree(SmoothDataOut);
      myfree(SmoothDataResult);
      myfree(SmoothDataGet);
    }

  myfree(sendcount_matrix);
  myfree(ndone_list);
  myfree(recv_offset);
  myfree(recv_count);
  myfree(send_offset);
  myfree(send_count);

  myfree(DataNodeList);
  myfree(DataIndexTable);

  myfree(Ngblist);
  myfree(Exportnodecount);
  myfree(Exportindex);
  myfree(Exportflag);



  /* do final operations on results */
  for(i = 0; i < N_gas; i++)
    {
      if(density_isactive(i))
	if(P[i].Ti_endstep == All.Ti_Current)
	  {

#ifdef CONDUCTION
	    SphP[i].SmoothedEntr /= pow(SphP[i].d.Density, GAMMA);
#ifdef CONDUCTION_SATURATION
	    SphP[i].GradEntr[0] /= pow(SphP[i].d.Density, GAMMA);
	    SphP[i].GradEntr[1] /= pow(SphP[i].d.Density, GAMMA);
	    SphP[i].GradEntr[2] /= pow(SphP[i].d.Density, GAMMA);
#endif
#endif /* CONDUCTION */


#ifdef CR_DIFFUSION
	    SphP[i].CR_SmoothE0 /= SphP[i].d.Density;
	    SphP[i].CR_Smoothn0 /= SphP[i].d.Density;

	    /* Limit smoothed values so small-value particles
	     * will not be pulled into negative energy regimes
	     * by their neighbors due to far-too-high smoothed
	     * estimate */
	    if(SphP[i].CR_SmoothE0 > 2.0 * SphP[i].CR_E0)
	      SphP[i].CR_SmoothE0 = 2.0 * SphP[i].CR_E0;

	    if(SphP[i].CR_Smoothn0 > 2.0 * SphP[i].CR_n0)
	      SphP[i].CR_Smoothn0 = 2.0 * SphP[i].CR_n0;

#endif /* CR_DIFFUSION */

#ifdef SMOOTH_PHI
	    SphP[i].SmoothPhi /= SphP[i].d.Density;
#endif /* SMOOTH_PHI */
	  }
    }
}



/*! This function represents the core of the SPH density computation. The
*  target particle may either be local, or reside in the communication
*  buffer.
*/
void compute_smoothed_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int j, n, listindex = 0;
  int startnode, numngb_inbox;
  double h, h2, hinv, hinv3, hinv4;
  double wk, dwk;
  double dx, dy, dz, r, r2, u, mass_j;
  MyFloat *pos;

#ifdef CONDUCTION
  double smoothentr;

#ifdef CONDUCTION_SATURATION
  double gradentr[3];
#endif
#endif /* CONDUCTION */

#ifdef CR_DIFFUSION
  double rCR_SmoothE0 = 0.0;
  double rCR_Smoothn0 = 0.0;
#endif /* CR_DIFFUSION */

#ifdef SMOOTH_PHI
  double SmoothPhi = 0.0;
#endif /* SMOOTH_PHI */

#ifdef CONDUCTION
  smoothentr = 0;

#ifdef CONDUCTION_SATURATION
  gradentr[0] = gradentr[1] = gradentr[2] = 0;
#endif
#endif /* CONDUCTION */

  if(mode == 0)
    {
      pos = P[target].Pos;
      h = PPP[target].Hsml;
    }
  else
    {
      pos = SmoothDataGet[target].Pos;
      h = SmoothDataGet[target].Hsml;
    }


  h2 = h * h;
  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif
  hinv4 = hinv3 * hinv;



  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = SmoothDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  numngb_inbox = ngb_treefind_variable(&pos[0], h, target, &startnode, mode, nexport, nsend_local);

	  for(n = 0; n < numngb_inbox; n++)
	    {
	      j = Ngblist[n];

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
		      dwk = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
		    }
		  else
		    {
		      wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
		      dwk = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
		    }

		  mass_j = P[j].Mass;

#ifdef CONDUCTION
		  smoothentr += mass_j * wk * pow(SphP[j].d.Density, GAMMA_MINUS1) * SphP[j].Entropy;

#ifdef CONDUCTION_SATURATION
		  if(r > 0)
		    {
		      gradentr[0] +=
			mass_j * dwk * dx / r * SphP[j].Entropy * pow(SphP[j].d.Density, GAMMA_MINUS1);
		      gradentr[1] +=
			mass_j * dwk * dy / r * SphP[j].Entropy * pow(SphP[j].d.Density, GAMMA_MINUS1);
		      gradentr[2] +=
			mass_j * dwk * dz / r * SphP[j].Entropy * pow(SphP[j].d.Density, GAMMA_MINUS1);
		    }
#endif
#endif /* CONDUCTION */

#ifdef CR_DIFFUSION
		  rCR_SmoothE0 += mass_j * wk * SphP[j].CR_E0;
		  rCR_Smoothn0 += mass_j * wk * SphP[j].CR_n0;
#endif /* CR_DIFFUSION */

#ifdef SMOOTH_PHI
		  SmoothPhi += mass_j * wk * SphP[j].PhiPred;
#endif /* SMOOTH_PHI */
		}
	    }
	}
      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = SmoothDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }


  if(mode == 0)
    {

#ifdef CONDUCTION
      SphP[target].SmoothedEntr = smoothentr;
#ifdef CONDUCTION_SATURATION
      SphP[target].GradEntr[0] = gradentr[0];
      SphP[target].GradEntr[1] = gradentr[1];
      SphP[target].GradEntr[2] = gradentr[2];
#endif
#endif /* CONDUCTION */

#ifdef CR_DIFFUSION
      SphP[target].CR_SmoothE0 = rCR_SmoothE0;
      SphP[target].CR_Smoothn0 = rCR_Smoothn0;
#endif /* CR_DIFFUSION */

#ifdef SMOOTH_PHI
      SphP[target].SmoothPhi = SmoothPhi;
#endif /* SMOOTH_PHI */

    }
  else
    {

#ifdef CONDUCTION
      SmoothDataResult[target].SmoothedEntr = smoothentr;
#ifdef CONDUCTION_SATURATION
      SmoothDataResult[target].GradEntr[0] = gradentr[0];
      SmoothDataResult[target].GradEntr[1] = gradentr[1];
      SmoothDataResult[target].GradEntr[2] = gradentr[2];
#endif
#endif /* CONDUCTION */

#ifdef CR_DIFFUSION
      SmoothDataResult[target].CR_SmoothE0 = rCR_SmoothE0;
      SmoothDataResult[target].CR_Smoothn0 = rCR_Smoothn0;
#endif /* CR_DIFFUSION */

#ifdef SMOOTH_PHI
      SmoothDataResult[target].SmoothPhi = SmoothPhi;
#endif /* SMOOTH_PHI */
    }
}


#endif
