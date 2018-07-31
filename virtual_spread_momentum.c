#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

/*! \file virtual_spread_momentum.c \brief this routine transfer momentum of
 *  feedback/virtual particles to the SPH particles that contain it inside
 *  their smoothing radius. This routine should be called after the density at
 *  the virtual particle is determined (currently in virtaul_density.c), so
 *  that the correct normalisation of the momentum transferred is
 *  achieved. 
 */

#if defined(VIRTUAL)



/* "static" means that the structure is only used here and not passed outside
 * of this program */

static struct virtualdata_in
{
    MyDouble Pos[3];
    MyFloat Hsml;
    MyFloat NewDensity;
#if defined(REAL_EOS)
    MyFloat NewKappa;
#endif
    MyFloat DeltaPhotonMomentum;
#if defined(PHOTONS) && defined(VIRTUAL_HEATING)
    MyFloat DeltaHeat;
#endif
    MyFloat  WindOrPhoton;
    MyFloat Vel[3];
    MyIDType ID;
    int Index;
    int NodeList[NODELISTLENGTH];
}
*VirtualDataIn, *VirtualDataGet; /* These are things passed between processors */

static struct virtualdata_out
{
    MyLongDouble NewDensity;
    MyLongDouble DeltaPhotonMomentum;
}
*VirtualDataResult, *VirtualDataOut; /* These are structures local to the
					  * processor */

static int N_virt_gas=0;

static double hubble_a, ascale;

void virtual_spread_momentum(void)
{
  int i, j, k, n, bin;
  int ndone_flag, ndone;
  int ngrp, sendTask, recvTask, place, nexport, nimport, dummy;
  int Ntot_virt_gas=0;
  double rho;
  MPI_Status status;

//  if(ThisTask == 0)
//    {
//      printf("Determine which virtual particles have SPH neighbors\n");
//      fflush(stdout);
//    }

//  CPU_Step[CPU_MISC] += measure_time();

  if(All.ComovingIntegrationOn)
    {
      ascale = All.Time;
      hubble_a = hubble_function(All.Time);
    }
  else
    hubble_a = ascale = 1;

/* this calculation goes over all the active particles local to processor */

  /* allocate buffers to arrange communication */

  Ngblist = (int *) mymalloc(NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct virtualdata_in) +
					     sizeof(struct virtualdata_out) +
					     sizemax(sizeof(struct virtualdata_in),
						     sizeof(struct virtualdata_out))));
  DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));


  /** Let's first determine which particles may be swalled by whom */

  i = FirstActiveParticle;  /* first particle for this task */ 

  do
    {
      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */
      /* In particular, for each local black hole, go over the LOCAL neighbor
       * list and determine which particles can be accreted. This is mode = 0. */

      /* Then determine which of the local black holes have SPH neighbors in
       * other processors. This creates the export Table of something like

       particleID | host processor | processor containing SPH neighbors
      */

      /* The coordinates, smoothing radius, etc., of the black holes
       (corresponding to ParticleID) will be passed from the host (local)
       processor to the one with the SPH neighbors. The latter will then do
       the appropriate calculations and return the result to the processor in
       which black hole resides
       */

        for(nexport = 0; i >= 0; i = NextActiveParticle[i])
	  if(P[i].Type == 3)
	      if (P[i].NewDensity > 0 && P[i].DtJump > 0)
//	      if (P[i].NewDensity > 0)
	      {
//		  if (ThisTask == 0)
//		  {
//		      printf("Have photon i = %d, with density = %g, OldMomentum = %g \n", 
//			     i, P[i].NewDensity, P[i].OldPhotonMomentum); 
//		      fflush(stdout); 
//		  }


/* 0 is for local particles */
	      if(virtual_evaluate_onthespot(i, 0, &nexport, Send_count) < 0)
		  break;
	  }
/* break is needed if communication buffer is full. In that case it is first
 * emptied out and then calculation restarts */


/* The sort here re-arranges the export Table in a way to minimize the number
 * of MPI communication calls to improve efficiency of the run */

#ifdef MYSORT
      mysort_dataindex(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#else
      qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#endif

/* This passes the Export Table between the processes */
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

/* These are the things that need to be received == exported */
      VirtualDataGet = (struct virtualdata_in *) mymalloc(nimport * sizeof(struct virtualdata_in));
/* Same for things imported */
      VirtualDataIn = (struct virtualdata_in *) mymalloc(nexport * sizeof(struct virtualdata_in));

      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  for(k = 0; k < 3; k++)
	    {
		VirtualDataIn[j].Pos[k] = P[place].Pos[k];
		VirtualDataIn[j].Vel[k] = P[place].Vel[k];
	    }

	  VirtualDataIn[j].Hsml = PPP[place].Hsml;
	  VirtualDataIn[j].NewDensity = P[place].NewDensity;
#if defined(REAL_EOS)
          VirtualDataIn[j].NewKappa = P[place].NewKappa;
#endif
	  VirtualDataIn[j].WindOrPhoton = P[place].WindOrPhoton;
	  VirtualDataIn[j].DeltaPhotonMomentum = P[place].DeltaPhotonMomentum;
#if defined(PHOTONS) && defined(VIRTUAL_HEATING)
	  VirtualDataIn[j].DeltaHeat = P[place].DeltaHeat;
#endif
	  VirtualDataIn[j].ID = P[place].ID;

	  memcpy(VirtualDataIn[j].NodeList,
		 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
	}

/* This is a HORRIBLE MPI stuff. The communication between processors happen
 * pair-wise, i.e. pairs of processors identified that are allowed to exchange
 * data (both send and receive from each other). All combinations of these
 * pairs are considered to cover all the communications possible.
 */


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
		  MPI_Sendrecv(&VirtualDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct virtualdata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A,
			       &VirtualDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct virtualdata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
		}
	    }
	}

/* free up memory now that the data were passed around*/
      myfree(VirtualDataIn);

      VirtualDataResult = (struct virtualdata_out *) mymalloc(nimport * sizeof(struct virtualdata_out));
      VirtualDataOut = (struct virtualdata_out *) mymalloc(nexport * sizeof(struct virtualdata_out));

      /* now consider SPH particles that were sent to us. Thus mode=1 below */

      for(j = 0; j < nimport; j++)
      {
	  place = DataIndexTable[j].Index;
//	  if (P[place].NewDensity > 0)  virtual_evaluate_onthespot(j, 1, &dummy, &dummy);
	  virtual_evaluate_onthespot(j, 1, &dummy, &dummy);
      }
/* This has determined which non-local SPH particles will accrete which local
 * virtual particles */

      if(i < 0)
	ndone_flag = 1;
      else
	ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      /* Now return the results into the processors where the exported blackholes reside */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;
	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* send the results */
		  MPI_Sendrecv(&VirtualDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct virtualdata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B,
			       &VirtualDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct virtualdata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, &status);
		}
	    }

	}

      myfree(VirtualDataOut);
      myfree(VirtualDataResult);
      myfree(VirtualDataGet);
    }
  while(ndone < NTask);


  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);

  MPI_Reduce(&N_virt_gas, &Ntot_virt_gas, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
//  MPI_Reduce(&N_BH_swallowed, &Ntot_BH_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

/*  if(ThisTask == 0)
    {
      printf("Photon-SPH interaction done: %d gas particles interacted \n",
	     Ntot_virt_gas);
      fflush(stdout);
    }
*/

  if(ThisTask == 0)
    {
//      printf("passing photon momentum to SPH particles done .... \n");
//      fflush(stdout);
    }
}


/* this routine is to go through the SPH neighbors of the photon and pass them
 * all of the photon's momentum */

int virtual_evaluate_onthespot(int target, int mode, int *nexport, int *nSend_local)
{
  int startnode, numngb, j, k, n, index, id, listindex = 0;
  MyFloat *pos, *velocity, h_i, rho, h_virtual_search, h_host, mass, density, 
      delta_photon_momentum, delta_heat, windorphot;
  double dx, dy, dz, h_i2, r2, r, u, absol_vel, ddpm, dheat;
  double wk, hinv, hinv3;
#if defined(REAL_EOS)
  double opacity;
#endif

  if(mode == 0)
    {
      pos = P[target].Pos;
      h_i = PPP[target].Hsml;
      density = P[target].NewDensity;
#if defined(REAL_EOS)
      opacity =  P[target].NewKappa;
#endif
      delta_photon_momentum = P[target].DeltaPhotonMomentum;
#if defined(PHOTONS) && defined(VIRTUAL_HEATING)
      delta_heat = P[target].DeltaHeat;
#endif
      windorphot = P[target].WindOrPhoton;
      velocity = P[target].Vel;
      index = target;
      id = P[target].ID;

/*       printf("mode %d, Have phot i = %d, density = %g, delta_heat = %g, windorphot = %g \n", */
/* 	     mode, id, density, delta_heat, P[target].DeltaHeat);//windorphot ); */
/*       fflush(stdout); */

    }
  else
    {
      pos = VirtualDataGet[target].Pos;
      h_i = VirtualDataGet[target].Hsml;
      density = VirtualDataGet[target].NewDensity;
#if defined(REAL_EOS)
      opacity =  VirtualDataGet[target].NewKappa;
#endif
      delta_photon_momentum = VirtualDataGet[target].DeltaPhotonMomentum;
#if defined(PHOTONS) && defined(VIRTUAL_HEATING)
      delta_heat = VirtualDataGet[target].DeltaHeat;
#endif
      windorphot = VirtualDataGet[target].WindOrPhoton;
      velocity = VirtualDataGet[target].Vel;
      index = VirtualDataGet[target].Index;
      id = VirtualDataGet[target].ID;

/*       printf("mode %d, Have photon i = %d, with density = %g, delta_heat = %g, windorphot = %g \n", */
/* 	     mode, id, density,  delta_heat, VirtualDataGet[target].DeltaHeat);//windorphot ); */
/*       fflush(stdout); */

    }

  /* initialize variables before SPH loop is started */
  h_i2 = h_i * h_i;
//  h_virtual_search = All.FeedBackVelocity * dt * 2. + h_i + 0.1;
  h_virtual_search = h_i;

  /* Now start the actual SPH computation for this particle */
  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = VirtualDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  numngb = ngb_treefind_virtual_active(pos, h_virtual_search, target, &startnode, mode, nexport, nSend_local);

	  if(numngb < 0)
	    return -1;

	  absol_vel = sqrt(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2]);

	  for(n = 0; n < numngb; n++)
	    {
	      j = Ngblist[n];

	      if(P[j].Mass > 0)
	      {
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
		  u = sqrt(r2)/P[j].Hsml;

		  if (u < 1)
		  {
		      hinv = 1 /P[j].Hsml;
#ifndef  TWODIMS
		      hinv3 = hinv * hinv * hinv;
#else
		      hinv3 = hinv * hinv / boxSize_Z;
#endif
		      if(u < 0.5)
			  wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
		      else
			  wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);     
#ifndef REAL_EOS
                      ddpm = delta_photon_momentum* P[j].Mass* wk/density;
                      dheat = delta_heat * P[j].Mass* wk/density;
#else
                      ddpm = delta_photon_momentum* P[j].Mass*wk* SphP[j].Kappa/(density*opacity);
                      dheat = delta_heat * P[j].Mass*wk* SphP[j].Kappa/(density*opacity);
/*
                      fprintf (diag_shc,"VSM : %g %g %g %g %g\n",All.Time,density,opacity,P[target].NewKappa,VirtualDataGet[target].NewKappa);
                      fflush (diag_shc);
*/

#endif
#ifdef VIRTUAL_HEATING
		      
		      SphP[j].Injected_VIRTUAL_Energy += FLT(dheat *  windorphot);
#endif
		      for(k = 0; k < 3; k++)
			  SphP[j].i.Injected_BH_Momentum[k] += FLT(ddpm * velocity[k]/absol_vel);

		      N_virt_gas++;
		  }
		}
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = VirtualDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }

  return 0;
}

/*  Note: ngb_treefind_.... function used in this file is the same as the one in virtual_density.c
 */

#endif /* VIRTUAL */
