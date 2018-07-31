#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>


#include "allvars.h"
#include "proto.h"

#ifndef DEBUG
#define NDEBUG
#endif
#include <assert.h>

/*! \file smoothed_rad_accel.c
*  \brief Computation of SPH forces and rate of entropy generation
*
*  This file contains the "second SPH loop", where the SPH forces are
*  computed, and where the rate of change of entropy due to the shock heating
*  (via artificial viscosity) is computed.
*/


#ifdef SMOOTHED_RAD_ACCEL /* SHC */

struct hydrodata_in
{
  MyDouble Pos[3];
  MyFloat Hsml;
  MyFloat Mass;
  MyFloat Density;

  int NodeList[NODELISTLENGTH];
}
 *HydroDataIn, *HydroDataGet;


struct hydrodata_out
{
  MyLongDouble Acc[3];

}
 *HydroDataResult, *HydroDataOut;

static double hubble_a, atime, hubble_a2, fac_mu, fac_vsic_fix, a3inv, fac_egy;

/*! This function is the driver routine for the calculation of hydrodynamical
*  force and rate of change of entropy due to shock heating for all active
*  particles .
*/
void smoothed_rad_accel(void)
{
  int i, j, k, ngrp, ndone, ndone_flag, dummy;
  int sendTask, recvTask, nexport, nimport, place;
  double timecomp1 = 0 ;
  double tstart, tend;

  hubble_a = hubble_a2 = atime = fac_mu = fac_vsic_fix = a3inv = fac_egy = 1.0;

  /* allocate buffers to arrange communication */

  Ngblist = (int *) mymalloc(NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct hydrodata_in) +
					     sizeof(struct hydrodata_out) +
					     sizemax(sizeof(struct hydrodata_in),
						     sizeof(struct hydrodata_out))));
  DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));


  i = FirstActiveParticle;	/* first particle for this task */

  do
    {
      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */
      tstart = second();
      for(nexport = 0; i >= 0; i = NextActiveParticle[i])
	if(P[i].Type == 0)
	  {
	    if(sra_evaluate(i, 0, &nexport, Send_count) < 0)
	      break;
	  }
      tend = second();
      timecomp1 += timediff(tstart, tend);

#ifdef MYSORT
      mysort_dataindex(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#else
      qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#endif

      MPI_Allgather(Send_count, NTask, MPI_INT, Sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);

      for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	{
	  Recv_count[j] = Sendcount_matrix[j * NTask + ThisTask];
	  nimport += Recv_count[j];

	  if(j > 0)
	    {
	      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	    }
	}

      HydroDataGet = (struct hydrodata_in *) mymalloc(nimport * sizeof(struct hydrodata_in));
      HydroDataIn = (struct hydrodata_in *) mymalloc(nexport * sizeof(struct hydrodata_in));

      /* prepare particle data for export */

      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  for(k = 0; k < 3; k++)
	    {
	      HydroDataIn[j].Pos[k] = P[place].Pos[k];
	    }
	  HydroDataIn[j].Hsml = PPP[place].Hsml;
	  HydroDataIn[j].Mass = P[place].Mass;
	  HydroDataIn[j].Density = SphP[place].d.Density;

	  memcpy(HydroDataIn[j].NodeList,
		 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));


#ifdef PARTICLE_DEBUG
	  HydroDataIn[j].ID = P[place].ID;
#endif

#ifdef TIME_DEP_MAGN_DISP
	  HydroDataIn[j].Balpha = SphP[place].Balpha;
#endif
	}


      /* exchange particle data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* get the particles */
		  MPI_Sendrecv(&HydroDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct hydrodata_in), MPI_BYTE,
			       recvTask, TAG_HYDRO_A,
			       &HydroDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct hydrodata_in), MPI_BYTE,
			       recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}


      myfree(HydroDataIn);
      HydroDataResult = (struct hydrodata_out *) mymalloc(nimport * sizeof(struct hydrodata_out));
      HydroDataOut = (struct hydrodata_out *) mymalloc(nexport * sizeof(struct hydrodata_out));



      /* now do the particles that were sent to us */

      for(j = 0; j < nimport; j++)
	{
	  sra_evaluate(j, 1, &dummy, &dummy);
	}

      if(i < 0)
	ndone_flag = 1;
      else
	ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);


      /* get the result */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;
	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* send the results */
		  MPI_Sendrecv(&HydroDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct hydrodata_out),
			       MPI_BYTE, recvTask, TAG_HYDRO_B,
			       &HydroDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct hydrodata_out),
			       MPI_BYTE, recvTask, TAG_HYDRO_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}

      /* add the result to the local particles */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  for(k = 0; k < 3; k++)
	    {
	      SphP[place].ra.SmRadAccel[k] += HydroDataOut[j].Acc[k];
	    }

	}

      myfree(HydroDataOut);
      myfree(HydroDataResult);
      myfree(HydroDataGet);
    }
  while(ndone < NTask);


  myfree(DataNodeList);
  myfree(DataIndexTable);

  myfree(Ngblist);


}




/*! This function is the 'core' of the SPH force computation. A target
*  particle is specified which may either be local, or reside in the
*  communication buffer.
*/
int sra_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int startnode, numngb, listindex = 0;
  int j, k, n, timestep;
  MyDouble *pos, density;
  MyFloat h_i ;
  MyLongDouble acc[3], dtEntropy;

  double dx, dy, dz;
  double h_i2, hinv, hinv3;
  double h_j, wk_i, wk_j;
  double r, r2, u;

  if(mode == 0)
    {
	pos = P[target].Pos;
	h_i = PPP[target].Hsml;
	density = SphP[target].d.Density;
    }
  else
    {
	pos = HydroDataGet[target].Pos;
	h_i = HydroDataGet[target].Hsml;
	density = HydroDataGet[target].Density;
    }


  /* initialize variables before SPH loop is started */

  acc[0] = acc[1] = acc[2] = dtEntropy = 0;
  h_i2 = h_i * h_i;

  /* Now start the actual SPH computation for this particle */

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = HydroDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
#ifndef ONE_IN_KERNEL
	  numngb = ngb_treefind_pairs(pos, h_i, target, &startnode, mode, nexport, nsend_local);
#else
	  numngb = ngb_treefind_variable(pos, h_i, target, &startnode, mode, nexport, nsend_local);
#endif

	  if(numngb < 0)
	    return -1;

	  for(n = 0; n < numngb; n++)
	    {
	      j = Ngblist[n];

	      dx = pos[0] - P[j].Pos[0];
	      dy = pos[1] - P[j].Pos[1];
	      dz = pos[2] - P[j].Pos[2];
	      r2 = dx * dx + dy * dy + dz * dz;
#ifndef ONE_IN_KERNEL

	      h_j = PPP[j].Hsml;
	      if(r2 < h_i2 || r2 < h_j * h_j)
	      {
		  r = sqrt(r2);
		  if(r > 0)
		  {
		      if(r2 < h_i2)
		      {
			  hinv = 1.0 / h_i;
			  u = r * hinv;
			  hinv3 = hinv * hinv * hinv;
			  
			  if(u <= 0.5)
			      wk_i = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
			  else
			      wk_i = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
		      }
		      else
		      {
			  wk_i = 0;
		      }
		      
		      if(r2 < h_j * h_j)
		      {
			  hinv = 1.0 / h_j;
			  u = r * hinv;
			  hinv3 = hinv * hinv * hinv;
			  
			  if(u <= 0.5)
			      wk_j = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
			  else
			      wk_j = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
		      }
		      else
		      {
			  wk_j = 0;
		      }
		      
		      
		      acc[0] += P[j].Mass/SphP[j].d.Density * SphP[j].ra.RadAccel[0] * 0.5*(wk_i + wk_j);
		      acc[1] += P[j].Mass/SphP[j].d.Density * SphP[j].ra.RadAccel[1] * 0.5*(wk_i + wk_j);
		      acc[2] += P[j].Mass/SphP[j].d.Density * SphP[j].ra.RadAccel[2] * 0.5*(wk_i + wk_j);
		  }
	      }
#else
	      if(r2 < h_i2)
	      {
		  r = sqrt(r2);
		  hinv = 1.0 / h_i;
		  u = r * hinv;
		  hinv3 = hinv * hinv * hinv;
		  
		  if(u <= 0.5)
		      wk_i = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
		  else
		      wk_i = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
	      }
	      else
	      {
		  wk_i = 0;
	      }
	      
	      
	      acc[0] += P[j].Mass/density * SphP[j].ra.RadAccel[0] * wk_i;
	      acc[1] += P[j].Mass/density * SphP[j].ra.RadAccel[1] * wk_i;
	      acc[2] += P[j].Mass/density * SphP[j].ra.RadAccel[2] * wk_i;
#endif
	      
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = HydroDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }


  /* Now collect the result at the right place */
  if(mode == 0)
    {
      for(k = 0; k < 3; k++)
	SphP[target].ra.SmRadAccel[k] = acc[k];
    }
  else
    {
      for(k = 0; k < 3; k++)
	HydroDataResult[target].Acc[k] = acc[k];
    }

  return 0;
}

#endif /* SHC */
