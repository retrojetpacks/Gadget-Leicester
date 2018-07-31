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

/* eddington tensor computation and arrangement of particle communication*/
void eddington(void)
{
  int i, j, k, ngrp, dummy, cost_edd;
  int sendTask, recvTask, nexport = 0, nimport, place;
  char buf[100];
  MyLongDouble trace;

  /* allocate buffers to arrange communication */

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct eddingtondata_in) +
					     sizeof(struct eddingtondata_out) +
					     sizemax(sizeof(struct eddingtondata_in),
						     sizeof(struct eddingtondata_out))));
  DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));

  for(j = 0; j < NTask; j++)
    {
      Send_count[j] = 0;
      Exportflag[j] = -1;
    }

  /* do local particles and prepare export list */
  for(i = 0; i < NumPart; i++)
    if(P[i].Type == 0)
      {
	if(eddington_treeevaluate(i, 0, &nexport, Send_count) < 0)
	  break;
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

  EddingtonDataGet = (struct eddingtondata_in *) mymalloc(nimport * sizeof(struct eddingtondata_in));
  EddingtonDataIn = (struct eddingtondata_in *) mymalloc(nexport * sizeof(struct eddingtondata_in));

  /* prepare particle data for export */

  for(j = 0; j < nexport; j++)
    {
      place = DataIndexTable[j].Index;

      for(k = 0; k < 3; k++)
	{
	  EddingtonDataIn[j].Pos[k] = P[place].Pos[k];
	}
      EddingtonDataIn[j].Hsml = PPP[place].Hsml;


      memcpy(EddingtonDataIn[j].NodeList,
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
	      MPI_Sendrecv(&EddingtonDataIn[Send_offset[recvTask]],
			   Send_count[recvTask] * sizeof(struct eddingtondata_in), MPI_BYTE,
			   recvTask, TAG_HYDRO_A,
			   &EddingtonDataGet[Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(struct eddingtondata_in), MPI_BYTE,
			   recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }

  myfree(EddingtonDataIn);
  EddingtonDataResult = (struct eddingtondata_out *) mymalloc(nimport * sizeof(struct eddingtondata_out));
  EddingtonDataOut = (struct eddingtondata_out *) mymalloc(nexport * sizeof(struct eddingtondata_out));



  /* now do the particles that were sent to us */
  for(j = 0; j < nimport; j++)
    {
      cost_edd = eddington_treeevaluate(j, 1, &dummy, &dummy);
    }

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
	      MPI_Sendrecv(&EddingtonDataResult[Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(struct eddingtondata_out),
			   MPI_BYTE, recvTask, TAG_HYDRO_B,
			   &EddingtonDataOut[Send_offset[recvTask]],
			   Send_count[recvTask] * sizeof(struct eddingtondata_out),
			   MPI_BYTE, recvTask, TAG_HYDRO_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }

  /* add the result to the local particles */
  for(j = 0; j < nexport; j++)
    {
      place = DataIndexTable[j].Index;

      for(k = 0; k < 6; k++)
	{
	  SphP[place].ET[k] += EddingtonDataOut[j].ET[k];
	}
    }

  myfree(EddingtonDataOut);
  myfree(EddingtonDataResult);
  myfree(EddingtonDataGet);

  /* do final operations divide by the trace */
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)
	{
	  trace = SphP[i].ET[0] + SphP[i].ET[1] + SphP[i].ET[2];

	  if(trace)
	    {
	      for(k = 0; k < 6; k++)
		{
		  SphP[i].ET[k] /= trace;
		}
	    }
	  else
	    {
	      for(k = 0; k < 6; k++)
		{
		  SphP[i].ET[k] = 0.0;
		}
	    }
	}
    }

#ifdef PRINT_ET
  sprintf(buf, "%s%s%i%s%f%s", All.OutputDir, "eddington_", ThisTask, "_", All.Time, ".txt");
  FdEddington = fopen(buf, "wa");
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)
	{
	  printf("%i %i \n", ThisTask, i);
	  fprintf(FdEddington, "%f %f %f {{%f, %f, %f}, {%f, %f, %f}, {%f, %f, %f}} \n ",
		  P[i].Pos[0], P[i].Pos[1], P[i].Pos[2],
		  SphP[i].ET[0], SphP[i].ET[3], SphP[i].ET[5],
		  SphP[i].ET[3], SphP[i].ET[1], SphP[i].ET[4], SphP[i].ET[5], SphP[i].ET[4], SphP[i].ET[2]);
	}
    }
  fclose(FdEddington);
#endif


  myfree(DataNodeList);
  myfree(DataIndexTable);

}


/*! This routine computes the eddington tensor ET for a given local
 *  particle, or for a particle in the communication buffer. Depending on
 *  the value of TypeOfOpeningCriterion, either the geometrical BH
 *  cell-opening criterion, or the `relative' opening criterion is used.
 */
int eddington_treeevaluate(int target, int mode, int *nexport, int *nsend_local)
{
  struct NODE *nop = 0;
  int k, no, nodesinlist, nexport_save, ninteractions, task, listindex = 0, ptype;
  double r4, r2, dx, dy, dz, stellar_mass, h, h_inv, h3_inv;
  double pos_x, pos_y, pos_z;
  MyLongDouble ET[6] = { 0, 0, 0, 0, 0, 0 };
  MyLongDouble fac[6] = { 0, 0, 0, 0, 0, 0 };

#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  double soft = 0;
#endif
#ifdef PERIODIC
  double boxsize, boxhalf;

  boxsize = All.BoxSize;
  boxhalf = 0.5 * All.BoxSize;
#endif
  nexport_save = *nexport;

  ninteractions = 0;
  nodesinlist = 0;
  ptype = 0;

  if(mode == 0)
    {
      pos_x = P[target].Pos[0];
      pos_y = P[target].Pos[1];
      pos_z = P[target].Pos[2];
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      soft = P[target].Hsml;
#endif
    }
  else				/* mode=1 */
    {
      pos_x = EddingtonDataGet[target].Pos[0];
      pos_y = EddingtonDataGet[target].Pos[1];
      pos_z = EddingtonDataGet[target].Pos[2];
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      soft = EddingtonDataGet[target].Hsml;
#endif
    }

#ifndef UNEQUALSOFTENINGS
  h = All.ForceSoftening[ptype];
  h_inv = 1.0 / h;
  h3_inv = h_inv * h_inv * h_inv;
#endif

  if(mode == 0)
    {
      no = All.MaxPart;		/* root node */
    }
  else
    {
      no = EddingtonDataGet[target].NodeList[0];
      no = Nodes[no].u.d.nextnode;	/* open it */
    }

  while(no >= 0)
    {
      while(no >= 0)
	{
	  if(no < All.MaxPart)	/* single particle */
	    {
	      /* the index of the node is the index of the particle */
	      /* observe the sign */

	      if(P[no].Type == 4)
		{
		  dx = P[no].Pos[0] - pos_x;
		  dy = P[no].Pos[1] - pos_y;
		  dz = P[no].Pos[2] - pos_z;
		  stellar_mass = P[no].Mass;
		}
	      else
		{
		  dx = 0;
		  dy = 0;
		  dz = 0;
		  stellar_mass = 0;
		}
	    }

	  else			/* not a single particle */
	    {
	      if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
		{
		  if(mode == 0)
		    {
		      if(Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
			{
			  Exportflag[task] = target;
			  Exportnodecount[task] = NODELISTLENGTH;
			}

		      if(Exportnodecount[task] == NODELISTLENGTH)
			{
			  if(*nexport >= All.BunchSize)
			    {
			      /* out of buffer space. Need to discard work for this particle and interrupt */
			      *nexport = nexport_save;
			      if(nexport_save == 0)
				endrun(17998);	/* in this case, the buffer is too small to process even a single particle */
			      for(task = 0; task < NTask; task++)
				nsend_local[task] = 0;
			      for(no = 0; no < nexport_save; no++)
				nsend_local[DataIndexTable[no].Task]++;
			      return -1;
			    }
			  Exportnodecount[task] = 0;
			  Exportindex[task] = *nexport;
			  DataIndexTable[*nexport].Task = task;
			  DataIndexTable[*nexport].Index = target;
			  DataIndexTable[*nexport].IndexGet = *nexport;
			  *nexport = *nexport + 1;
			  nsend_local[task]++;
			}

		      DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] =
			DomainNodeIndex[no - (All.MaxPart + MaxNodes)];

		      if(Exportnodecount[task] < NODELISTLENGTH)
			DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;
		    }
		  no = Nextnode[no - MaxNodes];
		  continue;
		}

	      nop = &Nodes[no];

	      if(mode == 1)
		{
		  if(nop->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		    {
		      no = -1;
		      continue;
		    }
		}

	      stellar_mass = nop->stellar_mass;

	      if(!(nop->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
		{
		  /* open cell */
		  if(stellar_mass)
		    {
		      no = nop->u.d.nextnode;
		      continue;
		    }
		}

	      dx = nop->stellar_s[0] - pos_x;
	      dy = nop->stellar_s[1] - pos_y;
	      dz = nop->stellar_s[2] - pos_z;
	    }

	  r2 = (dx * dx) + (dy * dy) + (dz * dz);
	  r4 = r2 * r2;

	  if(no < All.MaxPart)
	    {
#ifdef UNEQUALSOFTENINGS
	      h = All.ForceSoftening[ptype];
	      if(h < All.ForceSoftening[P[no].Type])
		h = All.ForceSoftening[P[no].Type];
#endif
	      no = Nextnode[no];
	    }
	  else			/* we have an  internal node. Need to check opening criterion */
	    {
	      if(All.ErrTolTheta)	/* check Barnes-Hut opening criterion */
		{
		  if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
		    {
		      /* open cell */
		      no = nop->u.d.nextnode;
		      continue;
		    }
		}
	      else		/* check in addition whether we lie inside the cell */
		{

		  if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
		    {
		      if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
			{
			  if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
			    {
			      no = nop->u.d.nextnode;
			      continue;
			    }
			}
		    }
		}

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	      h = All.ForceSoftening[ptype];
	      if(h < All.ForceSoftening[extract_max_softening_type(nop->u.d.bitflags)])
		{
		  h = All.ForceSoftening[extract_max_softening_type(nop->u.d.bitflags)];
		  if(r2 < h * h)
		    {
		      if(maskout_different_softening_flag(nop->u.d.bitflags))	/* signals that there are particles of different softening in the node */
			{
			  no = nop->u.d.nextnode;
			  continue;
			}
		    }
		}
#else
	      if(ptype == 0)
		h = soft;
	      else
		h = All.ForceSoftening[ptype];

	      if(h < nop->maxsoft)
		{
		  h = nop->maxsoft;
		  if(r2 < h * h)
		    {
		      no = nop->u.d.nextnode;
		      continue;
		    }
		}
#endif
#endif
	      no = nop->u.d.sibling;	/* ok, node can be used */

	    }

	  if(r4 > 0)
	    {
	      fac[0] = stellar_mass * dx * dx / r4;
	      fac[1] = stellar_mass * dy * dy / r4;
	      fac[2] = stellar_mass * dz * dz / r4;
	      fac[3] = stellar_mass * dx * dy / r4;
	      fac[4] = stellar_mass * dy * dz / r4;
	      fac[5] = stellar_mass * dz * dx / r4;
	    }
	  else
	    {
	      for(k = 0; k < 6; k++)
		{
		  fac[k] = 0;
		}
	    }

	  for(k = 0; k < 6; k++)
	    {
	      ET[k] += fac[k];
	    }

	  if(stellar_mass > 0)
	    ninteractions++;
	}
      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      no = EddingtonDataGet[target].NodeList[listindex];
	      if(no >= 0)
		{
		  nodesinlist++;
		  no = Nodes[no].u.d.nextnode;	/* open it */
		}
	    }
	}
    }


  /* store result at the proper place */
  if(mode == 0)
    {
      for(k = 0; k < 6; k++)
	{
	  SphP[target].ET[k] = ET[k];
	}
    }
  else
    {
      for(k = 0; k < 6; k++)
	{
	  EddingtonDataResult[target].ET[k] = ET[k];
	}
      *nexport = nodesinlist;
    }


  return ninteractions;
}


#endif
