#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

/*! \file blackhole.c
 *  \brief routines for gas accretion onto black holes, and black hole mergers
 */

#ifndef NO_BH_ACCRETION 
#ifdef BLACK_HOLES


/* "static" means that the structure is only used here and not passed outside
 * of this program */

static struct blackholedata_in
{
  MyDouble Pos[3];
  MyFloat Density;
  MyFloat Mdot;
  MyFloat Dt;
  MyFloat Hsml;
  MyFloat Mass;
  MyFloat BH_Mass;
  MyFloat Vel[3];
  MyFloat Csnd;
  MyIDType ID;
    int Index;
    int NodeList[NODELISTLENGTH];
#ifdef LIMITED_ACCRETION_AND_FEEDBACK
    MyFloat AccDisc_Mass;
#endif
#ifdef DUST
  MyFloat Dust_Mass;
#endif

}
*BlackholeDataIn, *BlackholeDataGet; /* These are things passed between processors */

static struct blackholedata_out
{
  MyLongDouble Mass;
  MyLongDouble BH_Mass;
  MyLongDouble AccretedMomentum[3];
#ifdef REPOSITION_ON_POTMIN
  MyFloat BH_MinPotPos[3];
  MyFloat BH_MinPot;
#endif
#ifdef LIMITED_ACCRETION_AND_FEEDBACK
  MyLongDouble AccDisc_Mass;
#endif
#ifdef DUST
  MyLongDouble Dust_Mass;
#endif
}
*BlackholeDataResult, *BlackholeDataOut; /* These are structures local to the
					  * processor */

static double hubble_a, ascale;

static int N_gas_swallowed, N_BH_swallowed, N_dust_swallowed;


void blackhole_accretion(void)
{
  int i, j, k, n, bin;
  int ndone_flag, ndone;
  int ngrp, sendTask, recvTask, place, nexport, nimport, dummy;
  int Ntot_gas_swallowed, Ntot_BH_swallowed, Ntot_dust_swallowed;
  double mdot, rho, bhvel, soundspeed, meddington, dt, mdot_in_msun_per_year;
  double mass_real, total_mass_real, total_mdoteddington;
  double mass_holes, total_mass_holes, total_mdot;
  MPI_Status status;



  if(ThisTask == 0)
    {
      printf("Beginning black-hole accretion\n");
      fflush(stdout);
    }

  CPU_Step[CPU_MISC] += measure_time();

  if(All.ComovingIntegrationOn)
    {
      ascale = All.Time;
      hubble_a = hubble_function(All.Time);
    }
  else
    hubble_a = ascale = 1;
/* this here just records the accretion information into the appropriate log
 * files. I will not need this */

  for(n = 0; n < TIMEBINS; n++)
    {
#ifdef DUST
      //      TimeBin_Dust_Mass[n] = 0;
#endif

      if(TimeBinActive[n])
	{
	  TimeBin_BH_mass[n] = 0;
	  TimeBin_BH_dynamicalmass[n] = 0;
	  TimeBin_BH_Mdot[n] = 0;
	}
    }

  /* Let's first compute the Mdot values */

/* this calculation goes over all the Active local black holes */

  for(n = FirstActiveParticle; n >= 0; n = NextActiveParticle[n])
    if(P[n].Type == 5)
      {
	mdot = 0;		/* if no accretion model is enabled, we have mdot=0 */

	dt = (P[n].TimeBin ? (1 << P[n].TimeBin) : 0) * All.Timebase_interval / hubble_a;



#ifndef NO_BH_ACCRETION 
	rho = P[n].b1.BH_Density;
#else
	P[n].b1.BH_Density = 1.e-30;
	rho =  P[n].b1.BH_Density;
#endif
	bhvel = sqrt(pow(P[n].Vel[0] - P[n].b3.BH_SurroundingGasVel[0], 2) +
		     pow(P[n].Vel[1] - P[n].b3.BH_SurroundingGasVel[1], 2) +
		     pow(P[n].Vel[2] - P[n].b3.BH_SurroundingGasVel[2], 2));

	if(All.ComovingIntegrationOn)
	  {
	    bhvel /= All.Time;
	    rho /= pow(All.Time, 3);
	  }

	soundspeed = sqrt(GAMMA * P[n].b2.BH_Entropy * pow(rho, GAMMA_MINUS1));

	/* Note: we take here a radiative efficiency of 0.1 for Eddington accretion */

#ifndef BH_MASS_REAL
	meddington = (4 * M_PI * GRAVITY * C * PROTONMASS / (0.1 * C * C * THOMPSON)) * P[n].BH_Mass
	    * All.UnitTime_in_s;
#ifdef BONDI
	mdot = 4. * M_PI * All.BlackHoleAccretionFactor * All.G * All.G *
	    P[n].BH_Mass * P[n].BH_Mass * rho / pow((pow(soundspeed, 2) + pow(bhvel, 2)), 1.5);
#endif


#else /* BH_MASS_REAL */
	meddington = (4 * M_PI * GRAVITY * C * PROTONMASS / (0.1 * C * C * THOMPSON)) * P[n].Mass
	    * All.UnitTime_in_s;
#ifdef BONDI
	mdot = 4. * M_PI * All.BlackHoleAccretionFactor * All.G * All.G *
	    P[n].Mass * P[n].Mass * rho / pow((pow(soundspeed, 2) + pow(bhvel, 2)), 1.5);
#endif

#endif /* BH_MASS_REAL */


#ifdef STAR_FROMATION_FEEDBACK
/* Eddington limit for a star of radius 10^11 cm, transformed into code units
*/
	meddington = (4 * M_PI * C * PROTONMASS * 1.e11/ THOMPSON)/All.UnitMass_in_g * All.UnitTime_in_s;
#endif

#ifdef LIMITED_ACCRETION_AND_FEEDBACK
// CBP (31/10/2008) -- Compute accretion rate using accretion disc mass.
//
// The idea here is quite simple. P[n].Mass tracks the mass of the black hole 
// plus reservoir (which includes the accretion disc). P[n].BH_Mass is the 
// mass of the black hole that will, once the black hole particle starts 
// swallowing gas particles, be smaller than P[n].Mass. A fraction of the mass
// of P[n].Mass will be accreted onto the black hole in a given timestep. 
// This fraction will be dictated by the viscous timescale -- the timescale 
// over which mass is transported through the accretion disc.
// 
// To determine what this rate is, we solve the differential equation
// 
//    dMdisc/dt = Mdot - min(Meddington, Mdisc/tvisc, ...)
//
// where Meddington is the simple Eddington limit computed above.

#ifdef ACCRETION_RADIUS


#ifndef ACCRETION_AT_EDDINGTON_RATE

	mdot = P[n].AccDisc_Mass/(dt + (All.ViscousTimescale/All.UnitTime_in_s));
#ifdef ENFORCE_EDDINGTON_LIMIT
	if(mdot > All.BlackHoleEddingtonFactor * meddington)
	  mdot = All.BlackHoleEddingtonFactor * meddington;
#endif
	P[n].AccDisc_Mass -= mdot*dt;
	P[n].BH_Mass += mdot*dt;
	
#else /* if ACCRETION_AT_EDDINGTON_RATE defined */
	mdot = All.BlackHoleEddingtonFactor * meddington;
	P[n].BH_Mass += mdot*dt;

#endif /* ACCRETION_AT_EDDINGTON_RATE */



#endif
// End of CBP (31/10/2008)

#else /* LIMITED_ACCRETION_AND_FEEDBACK */

#ifdef ENFORCE_EDDINGTON_LIMIT
	if(mdot > All.BlackHoleEddingtonFactor * meddington)
	  mdot = All.BlackHoleEddingtonFactor * meddington;
#endif

#endif /* LIMITED_ACCRETION_AND_FEEDBACK */



/*
*/
#ifdef LIMITED_ACCRETION_AND_FEEDBACK
// CBP (31/10/2008) -- Print variables to screen for debugging
/*	printf("\n time %lf, Acc rate/medd, BH_Mass, Mass, Disc M:: %lf %lf %lf %lf \n\n",All.Time,
	       mdot/meddington,P[n].BH_Mass,P[n].Mass, P[n].AccDisc_Mass);
	fflush(stdout);
*/
// End of CBP (31/10/2008)
#endif

	P[n].BH_Mdot = mdot;

#ifndef DUST /* This is the default reporting */
	if(P[n].BH_Mass > 0)
	  {
#ifndef LONGIDS
	    fprintf(FdBlackHolesDetails, "BH=%u %g %g %g %g %g\n",
		    P[n].ID, All.Time, P[n].BH_Mass, mdot, rho, soundspeed);
#else
	    fprintf(FdBlackHolesDetails, "BH=%llu %g %g %g %g %g\n",
		    P[n].ID, All.Time, P[n].BH_Mass, mdot, rho, soundspeed);
#endif
	  }
#else
	if(P[n].Mass < 0.95*All.SMBHmass) 
	  {
	    fprintf(FdBlackHolesDetails, "%u  %g  %g  %g  %g \n",
		    P[n].ID, All.Time, P[n].Mass, P[n].Total_Mass, P[n].Dust_Mass);
	  }
#endif

// hack by CBP (03/11/2008) -- moved this to start of loop

//	dt = (P[n].TimeBin ? (1 << P[n].TimeBin) : 0) * All.Timebase_interval / hubble_a;

#ifdef BH_DRAG
	/* add a drag force for the black-holes,
	   accounting for the accretion */
	double fac;

	if(P[n].BH_Mass > 0)
	  {
	    /*
	       fac = P[n].BH_Mdot * dt / P[n].BH_Mass;
	     */
	    fac = meddington * dt / P[n].BH_Mass;

	    if(fac > 1)
	      fac = 1;

	    if(dt > 0)
	      for(k = 0; k < 3; k++)
		P[n].g.GravAccel[k] +=
		  -ascale * ascale * fac / dt * (P[n].Vel[k] - P[n].b3.BH_SurroundingGasVel[k]) / ascale;
	  }
#endif

// SN: done above inside ACCRETION_RADIUS
//	P[n].BH_Mass += P[n].BH_Mdot * dt;
      }


  /* Now let's invoke the functions that stochasticall swallow gas
   * and deal with black hole mergers.
   */

  if(ThisTask == 0)
    {
      printf("Start swallowing of gas particles and black holes\n");
      fflush(stdout);
    }


  N_gas_swallowed = N_BH_swallowed = N_dust_swallowed = 0;


  /* allocate buffers to arrange communication */

  Ngblist = (int *) mymalloc(NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct blackholedata_in) +
					     sizeof(struct blackholedata_out) +
					     sizemax(sizeof(struct blackholedata_in),
						     sizeof(struct blackholedata_out))));
  DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));


  /** Let's first spread the feedback energy, and determine which particles may be swalled by whom */

  i = FirstActiveParticle;	/* first particle for this task */

  do
    {
      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */
      /* In particular, for each local black hole, go other the LOCAL neighbor
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
	if(P[i].Type == 5)
/* blackhole_evaluate only determines which particles are to be swallowed by
 * which blackholes, but doesn't actually swallow any particles. This is done
 * to protect against same particle being swallowed by different
 * blackholes. Actual swallowing is done later in
 * blackhole_evaluate_swallow */


/* 0 is for local particles */
	  if(blackhole_evaluate(i, 0, &nexport, Send_count) < 0)
	    break;
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
      BlackholeDataGet = (struct blackholedata_in *) mymalloc(nimport * sizeof(struct blackholedata_in));
/* Same for things imported */
      BlackholeDataIn = (struct blackholedata_in *) mymalloc(nexport * sizeof(struct blackholedata_in));

      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  for(k = 0; k < 3; k++)
	    {
	      BlackholeDataIn[j].Pos[k] = P[place].Pos[k];
	      BlackholeDataIn[j].Vel[k] = P[place].Vel[k];
	    }

	  BlackholeDataIn[j].Hsml = PPP[place].Hsml;
	  BlackholeDataIn[j].Mass = P[place].Mass;
	  BlackholeDataIn[j].BH_Mass = P[place].BH_Mass;
#ifdef DUST
	  BlackholeDataIn[j].Dust_Mass = P[place].Dust_Mass;
#endif
	  BlackholeDataIn[j].Density = P[place].b1.BH_Density;
	  BlackholeDataIn[j].Mdot = P[place].BH_Mdot;
	  BlackholeDataIn[j].Csnd =
	    sqrt(GAMMA * P[place].b2.BH_Entropy *
		 pow(P[place].b1.BH_Density / (ascale * ascale * ascale), GAMMA_MINUS1));
	  BlackholeDataIn[j].Dt =
	    (P[place].TimeBin ? (1 << P[place].TimeBin) : 0) * All.Timebase_interval / hubble_a;
	  BlackholeDataIn[j].ID = P[place].ID;

	  memcpy(BlackholeDataIn[j].NodeList,
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
		  MPI_Sendrecv(&BlackholeDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A,
			       &BlackholeDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
		}
	    }
	}

/* free up memory now that the data were passed around*/
      myfree(BlackholeDataIn);

      BlackholeDataResult = (struct blackholedata_out *) mymalloc(nimport * sizeof(struct blackholedata_out));
      BlackholeDataOut = (struct blackholedata_out *) mymalloc(nexport * sizeof(struct blackholedata_out));


      /* now consider the blackhole particles that were sent to us. Thus mode=1 below */

      for(j = 0; j < nimport; j++)
	blackhole_evaluate(j, 1, &dummy, &dummy);
/* This has determined which non-local blackholes will accrete which local SPH particles */

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
		  MPI_Sendrecv(&BlackholeDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct blackholedata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B,
			       &BlackholeDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct blackholedata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, &status);
		}
	    }

	}

      /* add the result to the particles */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

#ifdef REPOSITION_ON_POTMIN
	  if(P[place].BH_MinPot > BlackholeDataOut[j].BH_MinPot)
	    {
	      P[place].BH_MinPot = BlackholeDataOut[j].BH_MinPot;
	      for(k = 0; k < 3; k++)
		P[place].BH_MinPotPos[k] = BlackholeDataOut[j].BH_MinPotPos[k];
	    }
#endif
	}

      myfree(BlackholeDataOut);
      myfree(BlackholeDataResult);
      myfree(BlackholeDataGet);
    }
  while(ndone < NTask);





  /* Now do actual the swallowing of particles */

  i = FirstActiveParticle;	/* first particle for this task */

  do
    {
      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */

      for(nexport = 0; i >= 0; i = NextActiveParticle[i])
	if(P[i].Type == 5)
	  if(P[i].SwallowID == 0)
	    if(blackhole_evaluate_swallow(i, 0, &nexport, Send_count) < 0)
	      break;


      qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);

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

      BlackholeDataGet = (struct blackholedata_in *) mymalloc(nimport * sizeof(struct blackholedata_in));
      BlackholeDataIn = (struct blackholedata_in *) mymalloc(nexport * sizeof(struct blackholedata_in));

      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  for(k = 0; k < 3; k++)
	    BlackholeDataIn[j].Pos[k] = P[place].Pos[k];

	  BlackholeDataIn[j].Hsml = PPP[place].Hsml;
	  BlackholeDataIn[j].BH_Mass = P[place].BH_Mass;
#ifdef DUST
	  BlackholeDataIn[j].Dust_Mass = P[place].Dust_Mass;
#endif
	  BlackholeDataIn[j].ID = P[place].ID;

	  memcpy(BlackholeDataIn[j].NodeList,
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
		  MPI_Sendrecv(&BlackholeDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A,
			       &BlackholeDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
		}
	    }
	}


      myfree(BlackholeDataIn);
      BlackholeDataResult = (struct blackholedata_out *) mymalloc(nimport * sizeof(struct blackholedata_out));
      BlackholeDataOut = (struct blackholedata_out *) mymalloc(nexport * sizeof(struct blackholedata_out));


      /* now do the particles that were sent to us */

      for(j = 0; j < nimport; j++)
	blackhole_evaluate_swallow(j, 1, &dummy, &dummy);

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
		  MPI_Sendrecv(&BlackholeDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct blackholedata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B,
			       &BlackholeDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct blackholedata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, &status);
		}
	    }

	}

      /* add the result to the particles */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  P[place].b4.dBH_accreted_Mass += BlackholeDataOut[j].Mass;
	  P[place].b5.dBH_accreted_BHMass += BlackholeDataOut[j].BH_Mass;
#ifdef DUST
	  P[place].b5.dBH_accreted_DustMass += BlackholeDataOut[j].Dust_Mass;
#endif
	  for(k = 0; k < 3; k++)
	    P[place].b6.dBH_accreted_momentum[k] += BlackholeDataOut[j].AccretedMomentum[k];
	}

      myfree(BlackholeDataOut);
      myfree(BlackholeDataResult);
      myfree(BlackholeDataGet);
    }
  while(ndone < NTask);


  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);


  MPI_Reduce(&N_gas_swallowed, &Ntot_gas_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&N_BH_swallowed, &Ntot_BH_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&N_dust_swallowed, &Ntot_dust_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("Accretion done: %d gas particles swallowed, %d BH particles swallowed, %d Dust particles swallowed\n",
	     Ntot_gas_swallowed, Ntot_BH_swallowed, Ntot_dust_swallowed);
      fflush(stdout);
    }



#ifdef REPOSITION_ON_POTMIN
  for(n = FirstActiveParticle; n >= 0; n = NextActiveParticle[n])
    if(P[n].Type == 5)
      for(k = 0; k < 3; k++)
	P[n].Pos[k] = P[n].BH_MinPotPos[k];
#endif

  for(n = FirstActiveParticle; n >= 0; n = NextActiveParticle[n])
    if(P[n].Type == 5)
      {
#ifdef FLTROUNDOFFREDUCTION
	P[n].b4.BH_accreted_Mass = FLT(P[n].b4.dBH_accreted_Mass);
	P[n].b5.BH_accreted_BHMass = FLT(P[n].b5.dBH_accreted_BHMass);
#ifdef DUST
	P[n].b5.BH_accreted_DustMass = FLT(P[n].b5.dBH_accreted_DustMass);
#endif
	for(k = 0; k < 3; k++)
	  P[n].b6.BH_accreted_momentum[k] = FLT(P[n].b6.dBH_accreted_momentum[k]);
#endif
	if(P[n].b4.BH_accreted_Mass > 0)
	  {
	    for(k = 0; k < 3; k++)
	      P[n].Vel[k] =
		(P[n].Vel[k] * P[n].Mass + P[n].b6.BH_accreted_momentum[k]) /
		(P[n].Mass + P[n].b4.BH_accreted_Mass);

	    P[n].Mass += P[n].b4.BH_accreted_Mass;
	    P[n].BH_Mass += P[n].b5.BH_accreted_BHMass;
#ifdef DUST
	    P[n].Dust_Mass += P[n].b5.BH_accreted_DustMass;
#endif
// CBP (31/10/2008) -- Addition of mass to disc
#ifdef LIMITED_ACCRETION_AND_FEEDBACK
	    P[n].AccDisc_Mass += P[n].b4.BH_accreted_Mass;
#ifdef ACCRETION_AT_EDDINGTON_RATE
	    P[n].Mass = P[n].BH_Mass;
#endif
#endif
#ifdef DUST
	    P[n].Total_Mass += P[n].b4.BH_accreted_Mass;
#endif
	    P[n].b4.BH_accreted_Mass = 0;
	  }
	bin = P[n].TimeBin;
#ifndef PLANET_ACCRETION_REPORTING /* This is the default reporting */
	TimeBin_BH_mass[bin] += P[n].BH_Mass;
	TimeBin_BH_dynamicalmass[bin] += P[n].Mass;
	TimeBin_BH_Mdot[bin] += P[n].BH_Mdot;
#else
	if (P[n].Mass < 0.9*All.SMBHmass) 
	  {
#ifdef DUST
	    TimeBin_BH_mass[bin] += P[n].Dust_Mass;
	    TimeBin_BH_dynamicalmass[bin] += P[n].Total_Mass;
	    TimeBin_BH_Mdot[bin] += P[n].Mass; /* SN: gives total planet mass NOW */
#else
	    TimeBin_BH_mass[bin] += P[n].BH_Mass;
	    TimeBin_BH_dynamicalmass[bin] += P[n].Mass;
	    TimeBin_BH_Mdot[bin] += P[n].BH_Mdot;
#endif
	    //	    TimeBin_Dust_Mass[bin] += P[n].Dust_Mass;
	  }
#endif
      }

  mdot = 0;
  mass_holes = 0;
  mass_real = 0;

  for(bin = 0; bin < TIMEBINS; bin++)
    if(TimeBinCount[bin])
      {
	mass_holes += TimeBin_BH_mass[bin];
	mass_real += TimeBin_BH_dynamicalmass[bin];
	mdot += TimeBin_BH_Mdot[bin];
      }

  MPI_Reduce(&mass_holes, &total_mass_holes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&mass_real, &total_mass_real, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&mdot, &total_mdot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      /* convert to solar masses per yr */
      mdot_in_msun_per_year =
	total_mdot * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

      if (total_mass_holes > 0)
	  total_mdoteddington = total_mdot / ((4 * M_PI * GRAVITY * C * PROTONMASS /
					       (0.1 * C * C * THOMPSON)) * total_mass_holes * All.UnitTime_in_s);
      else
	total_mdoteddington = 0.;
#ifdef PLANET_ACCRETION_REPORTING 
      All.rad_pl = sqrt((All.xpl-All.xbh)*(All.xpl-All.xbh) + (All.ypl-All.ybh)*(All.ypl-All.ybh)  + (All.zpl-All.zbh)*(All.zpl-All.zbh));
      total_mdoteddington = All.rad_pl;/* Planet-Star distance */

      /* SN: the print-out's format is: time, Star-Planet separation,
	 Mplanet_total (accreted from t=0), Mpl_dust (accreted from t=0), total planet mass */      

      /*      fprintf(FdBlackHoles, "%g %g %g %g %g\n",
	      All.Time, All.rad_pl, total_mass_real, total_mass_holes, mdot_in_msun_per_year);*/

      fprintf(FdBlackHoles, "%g %g %g %g %g\n",
	      All.Time, All.rad_pl, total_mass_real, total_mass_holes, total_mdot);

#else
      fprintf(FdBlackHoles, "%g %d %g %g %g %g %g\n",
	      All.Time, All.TotBHs, total_mass_holes, total_mdot, mdot_in_msun_per_year,
	      total_mass_real, total_mdoteddington);
#endif
      
      
      fflush(FdBlackHoles);
    }


  fflush(FdBlackHolesDetails);

  CPU_Step[CPU_BLACKHOLES] += measure_time();
}


int blackhole_evaluate(int target, int mode, int *nexport, int *nSend_local)
{
  int startnode, numngb, j, k, n, index, id, listindex = 0;
  MyFloat *pos, *velocity, h_i, dt, mdot, rho, mass, bh_mass, csnd;
  double dx, dy, dz, h_i2, r2, r, u, hinv, hinv3, wk, vrel;

  // AH changes //

  double acc_boundary;
  
  // End of AH changes //

#ifdef BH_KINETICFEEDBACK
  /*  double deltavel; */
  double activetime, activeenergy;
#endif
#ifdef BH_THERMALFEEDBACK
  double energy;
#endif
#ifdef REPOSITION_ON_POTMIN
  MyFloat minpotpos[3] = { 0, 0, 0 }, minpot = 1.0e30;
#endif

  if(mode == 0)
    {
      pos = P[target].Pos;
      rho = P[target].b1.BH_Density;
      mdot = P[target].BH_Mdot;
      dt = (P[target].TimeBin ? (1 << P[target].TimeBin) : 0) * All.Timebase_interval / hubble_a;
      h_i = PPP[target].Hsml;
      mass = P[target].Mass;
      bh_mass = P[target].BH_Mass;
      velocity = P[target].Vel;
      csnd =
	sqrt(GAMMA * P[target].b2.BH_Entropy *
	     pow(P[target].b1.BH_Density / (ascale * ascale * ascale), GAMMA_MINUS1));
      index = target;
      id = P[target].ID;
#ifdef BH_KINETICFEEDBACK
      activetime = P[target].ActiveTime;
      activeenergy = P[target].ActiveEnergy;
#endif

    }
  else
    {
      pos = BlackholeDataGet[target].Pos;
      rho = BlackholeDataGet[target].Density;
      mdot = BlackholeDataGet[target].Mdot;
      dt = BlackholeDataGet[target].Dt;
      h_i = BlackholeDataGet[target].Hsml;
      mass = BlackholeDataGet[target].Mass;
      bh_mass = BlackholeDataGet[target].BH_Mass;
      velocity = BlackholeDataGet[target].Vel;
      csnd = BlackholeDataGet[target].Csnd;
      index = BlackholeDataGet[target].Index;
      id = BlackholeDataGet[target].ID;
#ifdef BH_KINETICFEEDBACK
      activetime = BlackholeDataGet[target].ActiveTime;
      activeenergy = BlackholeDataGet[target].ActiveEnergy;
#endif
    }

//  h_i += 0.01;

  /* initialize variables before SPH loop is started */
  h_i2 = h_i * h_i;

  /* Now start the actual SPH computation for this particle */
  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = BlackholeDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  numngb = ngb_treefind_blackhole(pos, h_i, target, &startnode, mode, nexport, nSend_local);

	  if(numngb < 0)
	    return -1;
	  
//	  printf("processor N = %d, in blackhole.c -- Nnb = %d \n", ThisTask, numngb);

	  for(n = 0; n < numngb; n++)
	    {
	      j = Ngblist[n];

	      if(P[j].Mass > 0)
		{
		  if(mass > 0)
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

		      if(r2 < h_i2)
			{
#ifdef REPOSITION_ON_POTMIN
			  /* if this option is switched on, we may also encounter dark matter particles or stars */
			  if(P[j].p.Potential < minpot)
			    {
			      /* compute relative velocities */

			      for(k = 0, vrel = 0; k < 3; k++)
				vrel += (P[j].Vel[k] - velocity[k]) * (P[j].Vel[k] - velocity[k]);

			      vrel = sqrt(vrel) / ascale;

			      if(vrel <= 0.25 * csnd)
				{
				  minpot = P[j].p.Potential;
				  for(k = 0; k < 3; k++)
				    minpotpos[k] = P[j].Pos[k];
				}
			    }
#endif

#ifndef NO_BH_MERGERS
			  if(P[j].Type == 5)	/* we have a black hole merger */
			    {
			      if(r2 > 0)
				{
				  /* compute relative velocity of BHs */

#ifndef BH_MERGERS_WITHIN_H

				  for(k = 0, vrel = 0; k < 3; k++)
				    vrel += (P[j].Vel[k] - velocity[k]) * (P[j].Vel[k] - velocity[k]);

				  vrel = sqrt(vrel) / ascale;

				  if(vrel > 0.5 * csnd)

#else //Corresponds to BH_MERGERS_WITHIN_H

				  // AH changes //

				  if(mass > 0.95*All.SMBHmass) {
				      
				    acc_boundary = All.InnerBoundary;
					  
				  }
				  else {

				    acc_boundary = All.SofteningBndry;
				    //  acc_boundary = All.SinkBoundary;

				  }

				  if(pow(r2, 0.5) > acc_boundary)
				  //if(pow(r2, 0.5) > All.SofteningBndry)

				  // End of AH changes //

#endif
				      {
/*
#ifndef LONGIDS
				      fprintf(FdBlackHolesDetails,
					      "ThisTask=%d, time=%g: id=%u would like to swallow %u, but vrel=%g csnd=%g\n",
					      ThisTask, All.Time, id, P[j].ID, vrel, csnd);
#else
				      fprintf(FdBlackHolesDetails,
					      "ThisTask=%d, time=%g: id=%llu would like to swallow %llu, but vrel=%g csnd=%g\n",
					      ThisTask, All.Time, id, P[j].ID, vrel, csnd);
#endif
*/
				      }
				      else
				      {
					//				      if(P[j].SwallowID != id)
/* SN: I wonder if this will work better for mergers ? */
					  if(P[j].Mass <= mass && P[j].SwallowID != id)
//				          if(P[j].Mass > 0)
					  {
					      P[j].SwallowID = id;
					  }
				      else
				      {
					/*				      fprintf(FdBlackHolesDetails,
					      "ThisTask=%d, time=%g: id=%u would like to swallow %u, but SwallowID=%u and separation =%g, M1/M2 = %g \n",
					      ThisTask, All.Time, id, P[j].ID, P[j].SwallowID, pow(r2, 0.5), P[j].Mass/mass);*/
				      }
				    }
				}
			    }
#endif

#ifdef DUST
			  if(P[j].Type == 2) { /* dust accretion on the BH */
			      if(r2 > 0) {
				  /* compute relative velocity */

				  for(k = 0, vrel = 0; k < 3; k++)
				    vrel += (P[j].Vel[k] - velocity[k]) * (P[j].Vel[k] - velocity[k]);

				  vrel = sqrt(vrel) / ascale;
				  /* SN: use InnerBoundary for SMBH,
				     but SinkBoundary for smaller
				     sinks */
				   if(mass > 0.95*All.SMBHmass) { 
				       acc_boundary = All.InnerBoundary; 
				   }
				   else {
				       acc_boundary = All.SinkBoundary; 
				   } 
				  
				  if(pow(r2, 0.5)<acc_boundary && P[j].Mass>0.) {
				    double etotal = vrel*vrel/2. - mass/(r + 1.e-20);
				  /*
				    fprintf(FdBlackHolesDetails,
					    "ThisTask=%d, time=%g: id=%u swallows dust %u, but SwallowID=%u and separation =%g\n",
					    ThisTask, All.Time, id, P[j].ID, P[j].SwallowID, pow(r2, 0.5));
				  */
				    if (etotal <= 0.) P[j].SwallowID = id;
				  }
				}
			  }
#endif



			  if(P[j].Type == 0)
			    {
			      /* here we have a gas particle */

			      r = sqrt(r2);
			      hinv = 1 / h_i;
#ifndef  TWODIMS
			      hinv3 = hinv * hinv * hinv;
#else
			      hinv3 = hinv * hinv / boxSize_Z;
#endif

			      u = r * hinv;

			      if(u < 0.5)
				wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
			      else
				wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

#ifdef SWALLOWGAS
#ifndef ACCRETION_RADIUS
			      /* compute accretion probability */
			      double p, w;

			      if((bh_mass - mass) > 0)
				p = (bh_mass - mass) * wk / rho;
			      else
				p = 0;

			      /* compute random number, uniform in [0,1] */
			      w = get_random_number(P[j].ID);
			      if(w < p)
				{
				  if(P[j].SwallowID < id)
				    P[j].SwallowID = id;
				}
#else
			      for(k = 0, vrel = 0; k < 3; k++)
				  vrel += (P[j].Vel[k] - velocity[k]) * (P[j].Vel[k] - velocity[k]);
			      vrel = sqrt(vrel) / ascale;
			      double etotal = vrel*vrel/2. - mass/(r + 1.e-20);

			      //#ifndef NO_ACCRETION || ACCRETION_OF_DUST_ONLY 				 
#ifndef NO_ACCRETION
			  /* SN: disregard gas particle accretion if this flag is defined */


			      if(mass > 0.95*All.SMBHmass) 
				{ //for the SMBH
				  if(r < All.InnerBoundary) {
				    double usph = SphP[j].Entropy /GAMMA_MINUS1 * pow(SphP[j].d.Density+1.e-40, GAMMA_MINUS1);
				    double etotal = vrel*vrel/2. + usph - mass/(All.InnerBoundary + 1.e-10);
			      //				      if (etotal < 0.) 
				    //				    {
				    if(P[j].SwallowID < id)
				      P[j].SwallowID = id;
				      //a				    }
				    /*				    printf("ThisTask=%d, time=%g: id=%llu swallows %llu (%g %g)\n",
								    ThisTask, All.Time, id, P[j].ID, bh_mass, P[j].Mass);*/
				  }
				}
			      else 
				{
#ifndef ACCRETION_OF_DUST_ONLY 				 /* don't accrete gas on the SMBH in this case */
				  if(r < All.SinkBoundary) 
				    {   //for the other sinks
#ifndef ACCRETION_DENSITY
				      if (etotal < 0.) 
#else
					if( SphP[j].d.Density >=
					    All.CritOverDensity * pow(All.UnitLength_in_cm, 3.)/ All.UnitMass_in_g)
#endif
					  {
					    if(P[j].SwallowID < id)
					      P[j].SwallowID = id;
					  }
				    }
#endif /*  ACCRETION_OF_DUST_ONLY */
				}
#endif

#endif
#endif

			      if(P[j].Mass > 0)
				{
#ifdef BH_THERMALFEEDBACK
				  energy = 0.*All.BlackHoleFeedbackFactor * 0.1 * mdot * dt *
				    pow(C / All.UnitVelocity_in_cm_per_s, 2);

#ifdef TMP_FEEDBACK
#ifdef FRACTION_OF_LSOLAR_FB
				  /* SN: All.BlackHoleFeedbackFactor
				     is actually feedback luminosity
				     in erg/s 
				   energy = All.BlackHoleFeedbackFactor * dt *All.UnitTime_in_s
				     /All.UnitEnergy_in_cgs;*/
				  energy = 0.;
				  if (mass < 0.95*All.SMBHmass) {
				    energy = All.BlackHoleFeedbackFactor * SOLAR_LUM * dt
				      /All.UnitEnergy_in_cgs * All.UnitTime_in_s;
				  }
#else
				  /* SN: in this case All.BlackHoleFeedbackFactor is GM_p/R_p */
				  /*				   energy = All.BlackHoleFeedbackFactor * mdot * dt * 
								   All.UnitMass_in_g/All.UnitEnergy_in_cgs;*/
				  energy = 0.;
				  /* Feedback only from smaller sinks */
				  if (mass < 0.95*All.SMBHmass) {
				    energy = All.BlackHoleFeedbackFactor * pow(mass *All.UnitMass_in_g, 0.6667) * 6.67e-8 * 
				      pow(4.*3.1415/3. * 5., 0.3333)/All.UnitEnergy_in_cgs * mdot * All.UnitMass_in_g * dt;
				  }
				     
#endif
#endif

				  SphP[j].i.dInjected_BH_Energy += FLT(energy * P[j].Mass * wk / rho);
#endif
				}

			    }
			}
		    }
		}
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = BlackholeDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }



  /* Now collect the result at the right place */
  if(mode == 0)
    {
#ifdef REPOSITION_ON_POTMIN
      for(k = 0; k < 3; k++)
	P[target].BH_MinPotPos[k] = minpotpos[k];
      P[target].BH_MinPot = minpot;
#endif
    }
  else
    {
#ifdef REPOSITION_ON_POTMIN
      for(k = 0; k < 3; k++)
	BlackholeDataResult[target].BH_MinPotPos[k] = minpotpos[k];
      BlackholeDataResult[target].BH_MinPot = minpot;
#endif
    }

  return 0;
}


int blackhole_evaluate_swallow(int target, int mode, int *nexport, int *nSend_local)
{
  int startnode, numngb, j, k, n, id, listindex = 0;
  MyLongDouble accreted_mass, accreted_dust_mass, accreted_BH_mass, accreted_momentum[3];
  MyFloat *pos, h_i, bh_mass;


  if(mode == 0)
    {
      pos = P[target].Pos;
      h_i = PPP[target].Hsml;
      id = P[target].ID;
      bh_mass = P[target].BH_Mass;
    }
  else
    {
      pos = BlackholeDataGet[target].Pos;
      h_i = BlackholeDataGet[target].Hsml;
      id = BlackholeDataGet[target].ID;
      bh_mass = BlackholeDataGet[target].BH_Mass;
    }

//  h_i += 0.01;
  accreted_mass = 0;
  accreted_dust_mass = 0;
  accreted_BH_mass = 0;
  accreted_momentum[0] = accreted_momentum[1] = accreted_momentum[2] = 0;


  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = BlackholeDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }


  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  numngb = ngb_treefind_blackhole(pos, h_i, target, &startnode, mode, nexport, nSend_local);

	  if(numngb < 0)
	    return -1;

	  for(n = 0; n < numngb; n++)
	    {
	      j = Ngblist[n];

//accretion or merger happens when id of the victim particle coincides with the BH particle
	      if(P[j].SwallowID == id)
		{
		  if(P[j].Type == 5)	/* we have a black hole merger */
		    {
#ifndef LONGIDS
		      fprintf(FdBlackHolesDetails,
			      "ThisTask=%d, time=%g: id=%u swallows %u (%g %g)\n",
			      ThisTask, All.Time, id, P[j].ID, bh_mass, P[j].BH_Mass);
#else
		      fprintf(FdBlackHolesDetails,
			      "ThisTask=%d, time=%g: id=%llu swallows %llu (%g %g)\n",
			      ThisTask, All.Time, id, P[j].ID, bh_mass, P[j].BH_Mass);
#endif

		      accreted_mass += FLT(P[j].Mass);
		      accreted_BH_mass += FLT(P[j].BH_Mass);

		      for(k = 0; k < 3; k++)
			accreted_momentum[k] += FLT(P[j].Mass * P[j].Vel[k]);

		      P[j].Mass = 0;
		      P[j].BH_Mass = 0;

		      N_BH_swallowed++;
		    }
		}

		if(P[j].SwallowID == id) {
		  if(P[j].Type == 2) {
		      accreted_mass += FLT(P[j].Mass);
		      accreted_dust_mass += FLT(P[j].Mass);

		      for(k = 0; k < 3; k++)
#ifndef BH_FIXED
		      for(k = 0; k < 3; k++)
			accreted_momentum[k] += FLT(P[j].Mass * P[j].Vel[k]);
#endif
		      P[j].Mass = 0;
		      N_dust_swallowed++;
		    }
            }

	      if(P[j].Type == 0)
		{
		  if(P[j].SwallowID == id)
		    {
		      accreted_mass += FLT(P[j].Mass);
#ifndef BH_FIXED
		      for(k = 0; k < 3; k++)
			accreted_momentum[k] += FLT(P[j].Mass * P[j].Vel[k]);
#endif

		      P[j].Mass = 0.;
		      N_gas_swallowed++;
		    }
		}
	    }
	}
      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = BlackholeDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }

  /* Now collect the result at the right place */
  if(mode == 0)
    {
      P[target].b4.dBH_accreted_Mass = accreted_mass;
      P[target].b5.dBH_accreted_BHMass = accreted_BH_mass;
      P[target].b5.dBH_accreted_DustMass = accreted_dust_mass;
      for(k = 0; k < 3; k++)
	P[target].b6.dBH_accreted_momentum[k] = accreted_momentum[k];
    }
  else
    {
      BlackholeDataResult[target].Mass = accreted_mass;
      BlackholeDataResult[target].BH_Mass = accreted_BH_mass;
#ifdef  DUST
      BlackholeDataResult[target].Dust_Mass = accreted_dust_mass;
#endif
      for(k = 0; k < 3; k++)
	BlackholeDataResult[target].AccretedMomentum[k] = accreted_momentum[k];
    }

  return 0;
}




int ngb_treefind_blackhole(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode,
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

#ifndef REPOSITION_ON_POTMIN
#ifdef DUST
	  if(P[p].Type != 0 && P[p].Type != 5 && P[p].Type != 2)
	    continue;
#else
	  if(P[p].Type != 0 && P[p].Type != 5)
	    continue;
#endif
#endif
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
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
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

#endif
#endif

#if defined(SGRA_POTENTIAL) || defined(CUSP_POTENTIAL) || defined(SIS_POTENTIAL) || defined(NFW_POTENTIAL) || defined(QUASAR_HEATING) || defined(OTHIN_ACCELERATOR)  || defined(FIND_SMBH)
/* shares some SMBH data (currently only one allowed! 02.08.09) between
 * the different processors */
void   FindQuasars()
{
    double xbh = 0.0, tot_xbh = 0., ybh = 0.0, tot_ybh = 0., zbh = 0.0, tot_zbh = 0.,
	mbh=0, tot_mbh=0, lum=0, tot_lum=0,CurrentLuminosity ;
    double L_Edd_dimensionless = (4 * M_PI * GRAVITY * C * PROTONMASS/ THOMPSON)
	/pow(All.UnitVelocity_in_cm_per_s, 2) *   All.UnitTime_in_s;

    int blackhole_number=0, tot_num = 0, i;
    //    for(i = 0; i < NumPart; i++)
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
	if (P[i].Type == 5 && P[i].Mass > 0.9*All.SMBHmass) 
	{
	    blackhole_number++;
	    xbh = P[i].Pos[0];
	    ybh = P[i].Pos[1];
	    zbh = P[i].Pos[2];
	    mbh = P[i].Mass;
//#if defined(NFW_POTENTIAL)
#ifdef  FRACTION_OF_LEDD_FEEDBACK
//	    lum = All.VirtualFeedBack * L_Edd_dimensionless * P[i].Mass;
	    lum = L_Edd_dimensionless * P[i].Mass;
#else
//	    lum = All.VirtualFeedBack * 0.1 * P[i].BH_Mdot *  pow(C/All.UnitVelocity_in_cm_per_s, 2);
/* SN -- Feedback appropriate for stellar radiation sources */
/* The first bit is the accretion luminosity, assuming GM/R_core = 0.01 Msun/ 1Rsun */
/* The second bit is the main sequence luminosity of the star */
	    CurrentLuminosity = (1.989e31/All.UnitMass_in_g)/(7.e10/All.UnitLength_in_cm) * P[i].BH_Mdot;
//              + SOLAR_LUM/All.UnitEnergy_in_cgs * All.UnitTime_in_s *
//              pow(P[i].Mass * All.UnitMass_in_g/SOLAR_MASS, 3.3);
	    /* Case of a fixed protostellar luminosity */
//	    CurrentLuminosity = 10. * SOLAR_LUM/All.UnitEnergy_in_cgs * All.UnitTime_in_s;

/* Model in which Msink/Rsink grows as Msink */
//	    CurrentLuminosity = P[i].Mass/(7.e10/All.UnitLength_in_cm) * P[i].BH_Mdot;
/*  Now limit that at Eddington luminosity */
	  if (CurrentLuminosity > L_Edd_dimensionless * P[i].Mass) 
	      CurrentLuminosity = L_Edd_dimensionless * P[i].Mass;
	  lum = All.VirtualFeedBack * CurrentLuminosity;
	  /*
	  printf("M of the black hole is %g (code units) and L is %g (phys units) \n", 
		 P[i].Mass,  CurrentLuminosity * All.UnitEnergy_in_cgs/All.UnitTime_in_s);
	  fflush(stdout);
	  */
#endif
//#endif
	}
    }
    MPI_Allreduce(&xbh, &tot_xbh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&ybh, &tot_ybh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&zbh, &tot_zbh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&mbh, &tot_mbh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//#ifdef NFW_POTENTIAL 
    MPI_Allreduce(&lum, &tot_lum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//#endif
    MPI_Allreduce(&blackhole_number, &tot_num, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    //    if (ThisTask == 0) 
    {
	if (tot_num == 1)
	{
	  All.xbh = tot_xbh;
	  All.ybh = tot_ybh;
	  All.zbh = tot_zbh;
	  All.mbh = tot_mbh;
	  //#ifdef NFW_POTENTIAL 
#ifdef SMBH_IRRADIATION
	  All.BH_Luminosity =  All.BH_LuminosityFactor * SOLAR_LUM;
#endif
	}
	/* else */
	/* { */
	/*     printf("Error: Found %d black holes, but assumed 1! Have to stop!\n", tot_num); */
	/*     fflush(stdout); */
	/*     exit(1); */
	/* } */
    }

/*     All.xbh = tot_xbh; */
/*     All.ybh = tot_ybh; */
/*     All.zbh = tot_zbh; */
/*     All.mbh = tot_mbh; */
/* //#ifdef NFW_POTENTIAL  */
/* #ifdef SMBH_IRRADIATION */
/*     All.BH_Luminosity =  All.BH_LuminosityFactor * SOLAR_LUM; */
/* #endif */
//#endif
}
#endif


/* SN: this is for preheating by the planet */
#if defined(PLANET_IRRADIATION) || defined(PLANET_ACCRETION_FEEDBACK) || defined(PLANET_ACCRETION_REPORTING)
/* find position of the planet -- a blackhole type particle that is
   less massive than SMBH */
void   FindPlanet()
{
    double xbh = 0.0, tot_xbh = 0., ybh = 0.0, tot_ybh = 0., zbh = 0.0, tot_zbh = 0.,
	mbh=0, tot_mbh=0, lum=0, tot_lum=0,CurrentLuminosity ;

    int blackhole_number=0, tot_num = 0, i;
    for(i = 0; i < NumPart; i++)
    {
	if (P[i].Type == 5 && P[i].Mass < 0.9*All.SMBHmass) 
	{
	    blackhole_number++;
	    xbh = P[i].Pos[0];
	    ybh = P[i].Pos[1];
	    zbh = P[i].Pos[2];
	    mbh = P[i].Mass;
#ifdef PLANET_ACCRETION_FEEDBACK
	    CurrentLuminosity = 6.67e-8 * mbh * All.UnitMass_in_g/All.OuterBoundary *
	      P[i].BH_Mdot *  All.UnitMass_in_g/All.UnitTime_in_s;
#endif
#ifdef PLANET_IRRADIATION
	    CurrentLuminosity =  All.BlackHoleFeedbackFactor * SOLAR_LUM;
#endif

	    	  printf("M of the PLANET is %g (code units) and Luminosity is %g (phys units) \n", 
		 mbh,  CurrentLuminosity );
		 fflush(stdout); 
	}
    }
    MPI_Allreduce(&xbh, &tot_xbh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&ybh, &tot_ybh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&zbh, &tot_zbh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&mbh, &tot_mbh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&blackhole_number, &tot_num, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (ThisTask == 0) 
    {
	if (tot_num == 1)
	{
	}
	else
	{
	    printf("Error: Found %d PLANETS, but assumed 1! Have to stop!\n", tot_num);
	    fflush(stdout);
	    exit(1);
	}
    }

    All.xpl = tot_xbh;
    All.ypl = tot_ybh;
    All.zpl = tot_zbh;
    All.mpl = tot_mbh;
    All.planet_luminosity = CurrentLuminosity;
}
#endif
