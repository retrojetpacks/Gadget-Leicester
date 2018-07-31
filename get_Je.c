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

#ifdef RADTRANSFER

#define eV_to_Hz 2.41838e14
#define c_light 2.9979e10/All.UnitVelocity_in_cm_per_s // in code units

static struct stardata_in
{
  MyDouble Pos[3];
  MyFloat Hsml;
  int NodeList[NODELISTLENGTH];
}
 *StarDataIn, *StarDataGet;


static struct stardata_out
{
  MyLongDouble Je;
}
 *StarDataResult, *StarDataOut;

void star_density(void)
{
  int i, j, dummy;
  int ngrp, sendTask, recvTask, place, nexport = 0, nimport;

  /* allocate buffers to arrange communication */


  Ngblist = (int *) mymalloc(NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct stardata_in) + sizeof(struct stardata_out) +
					     sizemax(sizeof(struct stardata_in),
						     sizeof(struct stardata_out))));
  DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));

  for(j = 0; j < NTask; j++)
    {
      Send_count[j] = 0;
      Exportflag[j] = -1;
    }

  /* do local particles and prepare export list */
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)
	{
	  if(star_density_evaluate(i, 0, &nexport, Send_count) < 0)
	    break;
	}
    }

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

  StarDataGet = (struct stardata_in *) mymalloc(nimport * sizeof(struct stardata_in));
  StarDataIn = (struct stardata_in *) mymalloc(nexport * sizeof(struct stardata_in));

  /* prepare particle data for export */
  for(j = 0; j < nexport; j++)
    {
      place = DataIndexTable[j].Index;

      StarDataIn[j].Pos[0] = P[place].Pos[0];
      StarDataIn[j].Pos[1] = P[place].Pos[1];
      StarDataIn[j].Pos[2] = P[place].Pos[2];
      StarDataIn[j].Hsml = PPP[place].Hsml;

      memcpy(StarDataIn[j].NodeList,
	     DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));

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
	      MPI_Sendrecv(&StarDataIn[Send_offset[recvTask]],
			   Send_count[recvTask] * sizeof(struct stardata_in), MPI_BYTE,
			   recvTask, TAG_DENS_A,
			   &StarDataGet[Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(struct stardata_in), MPI_BYTE,
			   recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }

  myfree(StarDataIn);
  StarDataResult = (struct stardata_out *) mymalloc(nimport * sizeof(struct stardata_out));
  StarDataOut = (struct stardata_out *) mymalloc(nexport * sizeof(struct stardata_out));


  /* now do the particles that were sent to us */

  for(j = 0; j < nimport; j++)
    star_density_evaluate(j, 1, &dummy, &dummy);

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
	      MPI_Sendrecv(&StarDataResult[Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(struct stardata_out),
			   MPI_BYTE, recvTask, TAG_DENS_B,
			   &StarDataOut[Send_offset[recvTask]],
			   Send_count[recvTask] * sizeof(struct stardata_out),
			   MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}

    }

  /* add the result to the local particles */
  for(j = 0; j < nexport; j++)
    {
      place = DataIndexTable[j].Index;
      if(P[place].Type == 0)
	{
	  SphP[place].Je += StarDataOut[j].Je;
	}
    }

  myfree(StarDataOut);
  myfree(StarDataResult);
  myfree(StarDataGet);
  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);

}


/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
int star_density_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int j, n, numngb;
  int startnode, listindex = 0;
  double h, h2, hinv, hinv3;
  MyFloat Je = 0;
  double wk;
  double dx, dy, dz, r, r2, u, mass_j;
  MyDouble *pos;

  if(mode == 0)
    {
      pos = P[target].Pos;
      h = PPP[target].Hsml;
    }
  else
    {
      pos = StarDataGet[target].Pos;
      h = StarDataGet[target].Hsml;
    }


  h2 = h * h;
  hinv = 1.0 / h;
  hinv3 = hinv * hinv * hinv;

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = StarDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  numngb = ngb_treefind_stars(pos, h, target, &startnode, mode, nexport, nsend_local);

	  if(numngb < 0)
	    return -1;

	  for(n = 0; n < numngb; n++)
	    {
	      j = Ngblist[n];

	      dx = pos[0] - P[j].Pos[0];
	      dy = pos[1] - P[j].Pos[1];
	      dz = pos[2] - P[j].Pos[2];

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

		  Je += FLT(mass_j * wk);

		}
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = StarDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }

  if(mode == 0)
    {
      if(All.ComovingIntegrationOn)
	SphP[target].Je = Je * c_light / All.Time / 4.0 / M_PI;
      else
	SphP[target].Je = Je * c_light / 4.0 / M_PI;
    }
  else
    {
      if(All.ComovingIntegrationOn)
	StarDataResult[target].Je = Je * c_light / All.Time / 4.0 / M_PI;
      else
	StarDataResult[target].Je = Je * c_light / 4.0 / M_PI;
    }

  return 0;
}

#endif
