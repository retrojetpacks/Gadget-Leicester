#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

/*! \file virtual.c
 *  \brief SN: this routine Determine dust density & other quantities on SPH particles
 */

#if defined(DUST_TWO_POPULATIONS)


/* "static" means that the structure is only used here and not passed outside
 * of this program */

static struct virtualdata_in
{
    MyDouble Pos[3];
  MyFloat NewDensity, d_MicroDustMass; /* SN: density of pebbles and density of micron sized dust on SPH particles */
    MyFloat Hsml;
#ifdef NON_GRAY_RAD_TRANSFER
    MyFloat VirtualTemperature; /* Temperature of gas neighbors of the photon */
#endif
    MyFloat NewKappa;
#if defined(VIRTUAL_HEATING)
/* Despite the misleading name, DeltaPhotonEnergy here is the fraction of
 * energy (dimensionless, defined from 0 for complete absorption to 1 for
 * absorption and complete re-emission -- scattering) returned back to the
 * photon in an interaction. That fraction is kernel-averaged over photon
 * neighbors in virtual_density.c 
 */
    MyFloat DeltaPhotonEnergy;
#endif
    MyIDType ID;
    int Index;
    int NodeList[NODELISTLENGTH];
}
*VirtualDataIn, *VirtualDataGet; /* These are things passed between processors */

static struct virtualdata_out
{
  MyLongDouble NewDensity, d_MicroDustMass;
    MyLongDouble NewKappa;
#if defined(VIRTUAL_HEATING) 
    MyLongDouble DeltaPhotonEnergy;
#endif
#ifdef NON_GRAY_RAD_TRANSFER
    MyLongDouble VirtualTemperature;
#endif

}
*VirtualDataResult, *VirtualDataOut; /* These are structures local to the
					  * processor */


static int N_virt_gas=0;

static double hubble_a, ascale;

void dust_density_on_gas(void)
{
  int i, j, k, n, bin;
  int ndone_flag, ndone;
  int ngrp, sendTask, recvTask, place, nexport, nimport, dummy;
  int Ntot_virt_gas=0;//  int Ntot_virt_swallowed;
  double rho;
  MPI_Status status;

  if(ThisTask == 0)
    {
      printf("Determine dust density on SPH particles\n");
      fflush(stdout);
    }


//  CPU_Step[CPU_MISC] += measure_time();

  if(All.ComovingIntegrationOn)
    {
      ascale = All.Time;
      hubble_a = hubble_function(All.Time);
    }
  else
    hubble_a = ascale = 1;

/* this calculation goes over all the active particles local to processor */

  for(n = FirstActiveParticle; n >= 0; n = NextActiveParticle[n])
      {
	  if (P[n].Type != 0) continue;

	  P[n].NewDensity = 0;
          P[n].NewKappa = 0; 
	  //P[n].d_MicroDustMass = 0.;

#if defined(VIRTUAL_HEATING) 
	  P[n].DeltaPhotonEnergy = 0;
#endif
#ifdef NON_GRAY_RAD_TRANSFER
	  P[n].VirtualTemperature = 0;
#endif

      }

  /* Now let's invoke the functions that stochasticall swallow gas
   * and deal with black hole mergers.
   */



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

//  i = n;	/* first particle for this task */

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
	  if(P[i].Type == 0)
	  {
	    /* 	      printf("dust_density_on_gas.c: Have SPH i = %d, with dust density = %g \n", i, P[i].NewDensity); 
		      fflush(stdout); */

/* virtual_evaluate_select only determines which particles are to be swallowed by
 * which blackholes, but doesn't actually swallow any particles. This is done
 * to protect against same particle being swallowed by different
 * blackholes. Actual swallowing is done later in
 * blackhole_evaluate_swallow */


/* HERE WE SET hsml OF THE PHOTONS TO SOME NUMBER FOR A TEST OF NEIGHBOR FINDING ROUTINE */
//	      PPP[i].Hsml = 0.1 * (All.Time - P[i].StellarAge) * All.FeedBackVelocity + 0.2;
/* 0 is for local particles */
	      if(virtual_evaluate_select(i, 0, &nexport, Send_count) < 0)
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
//	      VirtualDataIn[j].Vel[k] = P[place].Vel[k];
	    }

	  VirtualDataIn[j].Hsml = P[place].Hsml;
	  VirtualDataIn[j].NewDensity = P[place].NewDensity;
          VirtualDataIn[j].NewKappa = P[place].NewKappa;
	  //          VirtualDataIn[j].d_MicroDustMass = P[place].d_MicroDustMass;

#if defined(VIRTUAL_HEATING)
	  VirtualDataIn[j].DeltaPhotonEnergy = P[place].DeltaPhotonEnergy;
#endif
#ifdef NON_GRAY_RAD_TRANSFER
	  VirtualDataIn[j].VirtualTemperature = P[place].VirtualTemperature;
#endif
//	  VirtualDataIn[j].Mass = P[place].Mass;

/*	  VirtualDataIn[j].Dt = 
	  (P[place].TimeBin ? (1 << P[place].TimeBin) : 0) * All.Timebase_interval / hubble_a; */
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
//	  PPP[place].Hsml = 0.1 * (All.Time - P[place].StellarAge) * All.FeedBackVelocity + 0.1;
	  virtual_evaluate_select(j, 1, &dummy, &dummy);
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

      /* add the result to the particles */
      for(j = 0; j < nexport; j++)
      {
	  place = DataIndexTable[j].Index;
	  P[place].NewDensity += VirtualDataOut[j].NewDensity;
	  //          P[place].d_MicroDustMass += VirtualDataOut[j].d_MicroDustMass;
#if defined(VIRTUAL_HEATING) 
          P[place].DeltaPhotonEnergy += VirtualDataOut[j].DeltaPhotonEnergy;
#endif
#ifdef NON_GRAY_RAD_TRANSFER
	  P[place].VirtualTemperature += VirtualDataOut[j].VirtualTemperature;
#endif
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

/*  if(ThisTask == 0)
    {
      printf("Selection of virtual particle with neighbors done");
      fflush(stdout);
    }
*/
}


int virtual_evaluate_select(int target, int mode, int *nexport, int *nSend_local)
{
  int startnode, numngb, j, k, n, index, id, listindex = 0, numngb_in_search_radius;
  MyFloat *pos, *velocity, h_i, rho, h_virtual_search, h_host, mass, NewDensity;
  double dx, dy, dz, h_i2, r2, r, u, vrel;
  double wk, hinv, hinv3;
  double NewKappa, d_MicroDustMass;

#if defined(VIRTUAL_HEATING) 
  MyFloat return_energy;
#endif
#ifdef NON_GRAY_RAD_TRANSFER
    MyFloat virtual_temperature, uold; 
    double dmax1, dmax2, a3inv, ascale, time_hubble_a, hubble_a;
  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = hubble_function(All.Time);
      time_hubble_a = All.Time * hubble_a;
      ascale = All.Time;
    }
  else
    a3inv = ascale = time_hubble_a = hubble_a = 1.;
#endif

  if(mode == 0)
    {
      pos = P[target].Pos;
      h_i = PPP[target].Hsml;
      NewDensity = P[target].NewDensity;
      NewKappa = P[target].NewKappa;
      //      d_MicroDustMass = P[target].d_MicroDustMass;

#if defined(VIRTUAL_HEATING)
      return_energy = P[target].DeltaPhotonEnergy;
#endif
#ifdef NON_GRAY_RAD_TRANSFER
      virtual_temperature = P[target].VirtualTemperature; 
#endif
      index = target;
      id = P[target].ID;
    }
  else
    {
      pos = VirtualDataGet[target].Pos;
      h_i = VirtualDataGet[target].Hsml;
      NewDensity = VirtualDataGet[target].NewDensity;
      NewKappa = VirtualDataGet[target].NewKappa;
      //      d_MicroDustMass = VirtualDataGet[target].d_MicroDustMass;
#ifdef NON_GRAY_RAD_TRANSFER
      virtual_temperature = VirtualDataGet[target].VirtualTemperature; 
#endif
#if defined(VIRTUAL_HEATING)
      return_energy = VirtualDataGet[target].DeltaPhotonEnergy;
#endif
      index = VirtualDataGet[target].Index;
      id = VirtualDataGet[target].ID;
    }

/* I'm trying to follow density.c here  */
  NewDensity = 0.;
  NewKappa = 0.;
  //  d_MicroDustMass = 0.;

#if defined(VIRTUAL_HEATING) 
  return_energy = 0.;
#endif
#ifdef NON_GRAY_RAD_TRANSFER
  virtual_temperature = 0.;
#endif

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

		       NewDensity += FLT(P[j].Mass * wk);
		       NewKappa += FLT(P[j].Mass* P[j].DustRadius * wk);
		       //		       Dlog_a_dt += FLT(P[j].Mass* P[j].Mass*P[j].LogDustRadius_by_dt * wk);
		       //d_MicroDustMass +=  FLT(3. * P[j].Mass * P[j].Mass * P[j].LogDustRadius_by_dt * wk);
		       /* Protect against dust mass goind negative */
		       //if (P[i].d_MicroDustMass <= 1.e-10*P[i].Mass) P[i].d_MicroDustMass = 1.e-10*P[i].Mass;
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

  /* Now collect the result at the right place */
  if(mode == 0)
    {
	P[target].NewDensity = NewDensity;
	P[target].NewKappa = NewKappa;
	//	P[target].d_MicroDustMass = d_MicroDustMass;

#if defined(VIRTUAL_HEATING)
	P[target].DeltaPhotonEnergy = return_energy;
#endif
#ifdef NON_GRAY_RAD_TRANSFER
	P[target].VirtualTemperature = virtual_temperature;
#endif
    }
  else
    {
	VirtualDataResult[target].NewDensity = NewDensity;
	VirtualDataResult[target].NewKappa = NewKappa;
	//	VirtualDataResult[target].d_MicroDustMass = d_MicroDustMass;
    }

  return 0;
}

/* Note that this function requires that the smoothing length is already
 * known. This is important to ensure the right number of neighbors. Thus it
 * cannot be called for a virtual particle whose neighbors were not previously
 * found...
 */
int ngb_treefind_virtual_active(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode,
			   int *nexport, int *nsend_local)
{
  int numngb, no, p, task, nexport_save;
  struct NODE *current;
  MyDouble dx, dy, dz, dist;

#ifdef PERIODIC
  MyDouble xtmp;
#endif
  nexport_save = *nexport;

  numngb = 0;
  no = *startnode;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

	  if(P[p].Type != 2) /* SN: I only care about dust particles here */
	    continue;

	  dist = hsml;
	  dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]);
	  if(dz > dist)
	    continue;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;

	  Ngblist[numngb++] = p;
//	  printf("neighbor type is %d", P[p].Type)
	}
      else
	{

	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(12312);

	      if(Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
		{
		  Exportflag[task] = target;
		  Exportnodecount[task] = NODELISTLENGTH;
		}

	      if(Exportnodecount[task] == NODELISTLENGTH)
		{
		  if(*nexport >= All.BunchSize)
		    {
		      *nexport = nexport_save;
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

	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	
/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
		  return numngb;
		}
	    }

	  no = current->u.d.sibling;	/* in case the node can be discarded */

	  dist = hsml + 0.5 * current->len;;
	  dx = NGB_PERIODIC_LONG_X(current->center[0] - searchcenter[0]);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(current->center[1] - searchcenter[1]);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(current->center[2] - searchcenter[2]);
	  if(dz > dist)
	    continue;
	  /* now test against the minimal sphere enclosing everything */
	  dist += FACT1 * current->len;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;

	  no = current->u.d.nextnode;	/* ok, we need to open the node */
	}
    }

  *startnode = -1;
  return numngb;
}

#endif /* VIRTUAL */
