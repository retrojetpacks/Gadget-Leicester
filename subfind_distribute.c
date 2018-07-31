#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#ifdef SUBFIND

#include "allvars.h"
#include "proto.h"
#include "domain.h"
#include "fof.h"
#include "subfind.h"

void subfind_distribute_groups(void)
{
  int i, nexport = 0, nimport = 0, target, ngrp, sendTask, recvTask;
  struct group_properties *send_Group;

  /* count how many we have of each task */
  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(i = 0; i < Ngroups; i++)
    {
      target = (Group[i].GrNr - 1) % NTask;
      if(target != ThisTask)
	Send_count[target]++;
    }
  MPI_Allgather(Send_count, NTask, MPI_INT, Sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);

  for(i = 0, Recv_offset[0] = Send_offset[0] = 0; i < NTask; i++)
    {
      Recv_count[i] = Sendcount_matrix[i * NTask + ThisTask];

      nimport += Recv_count[i];
      nexport += Send_count[i];

      if(i > 0)
	{
	  Send_offset[i] = Send_offset[i - 1] + Send_count[i - 1];
	  Recv_offset[i] = Recv_offset[i - 1] + Recv_count[i - 1];
	}
    }

  send_Group = (struct group_properties *) mymalloc(nexport * sizeof(struct group_properties));

  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(i = 0; i < Ngroups; i++)
    {
      target = (Group[i].GrNr - 1) % NTask;
      if(target != ThisTask)
	{
	  send_Group[Send_offset[target] + Send_count[target]] = Group[i];
	  Send_count[target]++;

	  Group[i] = Group[Ngroups - 1];
	  Ngroups--;
	  i--;
	}
    }

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
	{
	  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
	    {
	      /* get the group info */
	      MPI_Sendrecv(&send_Group[Send_offset[recvTask]],
			   Send_count[recvTask] * sizeof(struct group_properties), MPI_BYTE,
			   recvTask, TAG_DENS_A,
			   &Group[Ngroups + Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(struct group_properties), MPI_BYTE,
			   recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }

  Ngroups += nimport;

  myfree(send_Group);
}




void subfind_distribute_particles(int mode)
{
  int nexport = 0, nimport = 0;
  int i, n, ngrp, target, n_requests;
  struct particle_data *partBuf;
  MPI_Request *requests;

  requests = (MPI_Request *) mymalloc(10 * NTask * sizeof(MPI_Request));

  for(n = 0; n < NTask; n++)
    Send_count[n] = 0;

  for(n = 0; n < NumPart; n++)
    {
      if(mode == 0)
	{
	  if(P[n].GrNr > Ncollective && P[n].GrNr <= TotNgroups)	/* particle is in small group */
	    {
	      target = (P[n].GrNr - 1) % NTask;
	      if(target != ThisTask)
		Send_count[target] += 1;
	    }
	}
      else if(mode == 1)
	{
	  if(P[n].targettask != ThisTask)
	    Send_count[P[n].targettask] += 1;
	}
      else if(mode == 2)
	{
	  if(P[n].origintask != ThisTask)
	    Send_count[P[n].origintask] += 1;
	}
    }

  MPI_Allgather(Send_count, NTask, MPI_INT, Sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);

  for(i = 0, Send_offset[0] = Recv_offset[0] = 0; i < NTask; i++)
    {
      Recv_count[i] = Sendcount_matrix[i * NTask + ThisTask];

      nimport += Recv_count[i];
      nexport += Send_count[i];

      if(i > 0)
	{
	  Send_offset[i] = Send_offset[i - 1] + Send_count[i - 1];
	  Recv_offset[i] = Recv_offset[i - 1] + Recv_count[i - 1];
	}
    }


  if(NumPart - nexport + nimport > All.MaxPart)
    {
      printf
	("on task=%d the maximum particle number All.MaxPart=%d is reached (NumPart=%d togo=%d toget=%d)",
	 ThisTask, All.MaxPart, NumPart, nexport, nimport);
      endrun(8765);
    }

  partBuf = (struct particle_data *) mymalloc(nexport * sizeof(struct particle_data));

  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(n = 0; n < NumPart; n++)
    {
      if(mode == 0)
	{
	  if(!(P[n].GrNr > Ncollective && P[n].GrNr <= TotNgroups))	/* particle is in small group */
	    continue;

	  target = (P[n].GrNr - 1) % NTask;
	}
      else if(mode == 1)
	{
	  target = P[n].targettask;
	}
      else if(mode == 2)
	{
	  target = P[n].origintask;
	}
      
      if(target != ThisTask)
	{
	  partBuf[Send_offset[target] + Send_count[target]] = P[n];
	  Send_count[target]++;
	  
	  P[n] = P[NumPart - 1];
	  NumPart--;
	  n--;
	}
    }

  for(ngrp = 1, n_requests = 0; ngrp < (1 << PTask); ngrp++)
    {
      target = ThisTask ^ ngrp;

      if(target < NTask)
	{
	  if(Send_count[target] > 0)
	    MPI_Isend(partBuf + Send_offset[target], Send_count[target] * sizeof(struct particle_data),
		      MPI_BYTE, target, TAG_PDATA, MPI_COMM_WORLD, &requests[n_requests++]);
	}
    }


  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = ThisTask ^ ngrp;

      if(target < NTask)
	{
	  if(Recv_count[target] > 0)
	    MPI_Irecv(P + NumPart + Recv_offset[target], Recv_count[target] * sizeof(struct particle_data),
		      MPI_BYTE, target, TAG_PDATA, MPI_COMM_WORLD, &requests[n_requests++]);
	}
    }

  MPI_Waitall(n_requests, requests, MPI_STATUSES_IGNORE);

  NumPart += nimport;

  myfree(partBuf);
  myfree(requests);
}

#endif
