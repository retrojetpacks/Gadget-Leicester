#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

#include "allvars.h"
#include "proto.h"




/*! \file io.c
 *  \brief Output of a snapshot file to disk.
 */

static int n_type[6];
static long long ntot_type_all[6];



/*! This function writes a snapshot of the particle distribution to one or
 * several files using Gadget's default file format.  If
 * NumFilesPerSnapshot>1, the snapshot is distributed into several files,
 * which are written simultaneously. Each file contains data from a group of
 * processors of size roughly NTask/NumFilesPerSnapshot.
 */
void savepositions(int num)
{
  size_t bytes;
  char buf[500];
  int i, j, *temp, n, filenr, gr, ngroups, masterTask, lastTask;

  if(ThisTask == 0)
    printf("\nwriting snapshot file... \n");

  CPU_Step[CPU_MISC] += measure_time();

#if defined(SFR) || defined(BLACK_HOLES)
  rearrange_particle_sequence();
  /* ensures that new tree will be constructed */
  All.NumForcesSinceLastDomainDecomp = (long long) (1 + All.TreeDomainUpdateFrequency * All.TotNumPart);
#endif

#ifdef ORDER_SNAPSHOTS_BY_ID
  double t0, t1;

  t0 = second();
  for(i = 0; i < NumPart; i++)
    {
      P[i].GrNr = ThisTask;
      P[i].SubNr = i;
    }
  parallel_sort(P, NumPart, sizeof(struct particle_data), io_compare_P_ID);
  t1 = second();
  if(ThisTask == 0)
    printf("Reordering of particle-data in ID-sequence took = %g sec\n", timediff(t0, t1));
#endif


  if(!(CommBuffer = mymalloc(bytes = All.BufferSize * 1024 * 1024)))
    {
      printf("failed to allocate memory for `CommBuffer' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(2);
    }

  if(NTask < All.NumFilesPerSnapshot)
    {
      if(ThisTask == 0)
	printf("Fatal error.\nNumber of processors must be larger or equal than All.NumFilesPerSnapshot.\n");
      endrun(0);
    }
  if(All.SnapFormat < 1 || All.SnapFormat > 3)
    {
      if(ThisTask == 0)
	printf("Unsupported File-Format\n");
      endrun(0);
    }
#ifndef  HAVE_HDF5
  if(All.SnapFormat == 3)
    {
      if(ThisTask == 0)
	printf("Code wasn't compiled with HDF5 support enabled!\n");
      endrun(0);
    }
#endif


  /* determine global and local particle numbers */
  for(n = 0; n < 6; n++)
    n_type[n] = 0;

  for(n = 0; n < NumPart; n++)
    n_type[P[n].Type]++;

/* SHC for the photon dump */
//n_type[3] = n_type[3] * All.FractionOfPhotons;
/* END SHC */

  /* because ntot_type_all[] is of type `long long', we cannot do a simple
   * MPI_Allreduce() to sum the total particle numbers 
   */
  temp = (int *) mymalloc(NTask * 6 * sizeof(int));
  MPI_Allgather(n_type, 6, MPI_INT, temp, 6, MPI_INT, MPI_COMM_WORLD);
  for(i = 0; i < 6; i++)
    {
      ntot_type_all[i] = 0;
      for(j = 0; j < NTask; j++)
	ntot_type_all[i] += temp[j * 6 + i];
    }
  myfree(temp);



  /* assign processors to output files */
  distribute_file(All.NumFilesPerSnapshot, 0, 0, NTask - 1, &filenr, &masterTask, &lastTask);

  if(All.NumFilesPerSnapshot > 1)
    sprintf(buf, "%s%s_%03d.%d", All.OutputDir, All.SnapshotFileBase, num, filenr);
  else
    sprintf(buf, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, num);


  ngroups = All.NumFilesPerSnapshot / All.NumFilesWrittenInParallel;
  if((All.NumFilesPerSnapshot % All.NumFilesWrittenInParallel))
    ngroups++;

  for(gr = 0; gr < ngroups; gr++)
    {
      if((filenr / All.NumFilesWrittenInParallel) == gr)	/* ok, it's this processor's turn */
	write_file(buf, masterTask, lastTask);
      MPI_Barrier(MPI_COMM_WORLD);
    }

  myfree(CommBuffer);

#ifdef ORDER_SNAPSHOTS_BY_ID
  t0 = second();
  parallel_sort(P, NumPart, sizeof(struct particle_data), io_compare_P_GrNr_SubNr);
  t1 = second();
  if(ThisTask == 0)
    printf("Restoring order of particle-data took = %g sec\n", timediff(t0, t1));
#endif

  if(ThisTask == 0)
    printf("done with snapshot.\n");

  CPU_Step[CPU_SNAPSHOT] += measure_time();

#ifdef FOF
  if(ThisTask == 0)
    printf("\ncomputing group catalogue...\n");

  fof_fof(num);

  if(ThisTask == 0)
    printf("done with group catalogue.\n");

  CPU_Step[CPU_FOF] += measure_time();
#endif
}



/*! This function fills the write buffer with particle data. New output blocks can in
 *  principle be added here.
 */
void fill_write_buffer(enum iofields blocknr, int *startindex, int pc, int type)
{
  int n, k, pindex, dt_step;
  MyOutputFloat *fp;
  double dmax1, dmax2;
  MyIDType *ip;

#ifdef PERIODIC
  MyFloat boxSize;
#endif
#ifdef PMGRID
  double dt_gravkick_pm = 0;
#endif
  double dt_gravkick, dt_hydrokick, a3inv = 1, fac1, fac2;

#if defined(COOLING)
  double ne, nh0, nHeII;
#endif
#ifdef OUTPUTCOOLRATE
  double tcool, u;
#endif

  if(All.ComovingIntegrationOn)
    {
      a3inv = 1 / (All.Time * All.Time * All.Time);
      fac1 = 1 / (All.Time * All.Time);
      fac2 = 1 / pow(All.Time, 3 * GAMMA - 2);
    }
  else
    a3inv = fac1 = fac2 = 1;

#ifdef PMGRID
  if(All.ComovingIntegrationOn)
    dt_gravkick_pm =
      get_gravkick_factor(All.PM_Ti_begstep,
			  All.Ti_Current) -
      get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);
  else
    dt_gravkick_pm = (All.Ti_Current - (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2) * All.Timebase_interval;
#endif

  fp = (MyOutputFloat *) CommBuffer;
  ip = (MyIDType *) CommBuffer;

  pindex = *startindex;

  switch (blocknr)
    {
    case IO_POS:		/* positions */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      {
		fp[k] = P[pindex].Pos[k];
#ifdef PERIODIC
		boxSize = All.BoxSize;
#ifdef LONG_X
		if(k == 0)
		  boxSize = All.BoxSize * LONG_X;
#endif
#ifdef LONG_Y
		if(k == 1)
		  boxSize = All.BoxSize * LONG_Y;
#endif
#ifdef LONG_Z
		if(k == 2)
		  boxSize = All.BoxSize * LONG_Z;
#endif
		while(fp[k] < 0)
		  fp[k] += boxSize;
		while(fp[k] >= boxSize)
		  fp[k] -= boxSize;
#endif
	      }
	    n++;
	    fp += 3;
	  }
      break;

    case IO_VEL:		/* velocities */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    dt_step = (P[pindex].TimeBin ? (1 << P[pindex].TimeBin) : 0);

	    if(All.ComovingIntegrationOn)
	      {
		dt_gravkick =
		  get_gravkick_factor(P[pindex].Ti_begstep, All.Ti_Current) -
		  get_gravkick_factor(P[pindex].Ti_begstep, P[pindex].Ti_begstep + dt_step / 2);
		dt_hydrokick =
		  get_hydrokick_factor(P[pindex].Ti_begstep, All.Ti_Current) -
		  get_hydrokick_factor(P[pindex].Ti_begstep, P[pindex].Ti_begstep + dt_step / 2);
	      }
	    else
	      dt_gravkick = dt_hydrokick =
		(All.Ti_Current - (P[pindex].Ti_begstep + dt_step / 2)) * All.Timebase_interval;

	    for(k = 0; k < 3; k++)
	      {
		fp[k] = P[pindex].Vel[k] + P[pindex].g.GravAccel[k] * dt_gravkick;
		if(P[pindex].Type == 0)
                {
		  fp[k] += SphP[pindex].a.HydroAccel[k] * dt_hydrokick;
#ifdef RAD_ACCEL
		  fp[k] += SphP[pindex].ra.RadAccel[k] * dt_hydrokick;
#endif

#ifdef GAS_FIXED
		  fp[k] = 0.;
#endif
	        }
	      }
#ifdef PMGRID
	    for(k = 0; k < 3; k++)
	      fp[k] += P[pindex].GravPM[k] * dt_gravkick_pm;
#endif
	    for(k = 0; k < 3; k++)
	      fp[k] *= sqrt(a3inv);

	    n++;
	    fp += 3;
	  }
      break;

    case IO_ID:		/* particle ID */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *ip++ = P[pindex].ID;
	    n++;
	  }
      break;

    case IO_MASS:		/* particle mass */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
        {
	    *fp++ = P[pindex].Mass;
	    n++;
        }
      break;

    case IO_U:			/* internal energy */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ =
	      DMAX(All.MinEgySpec,
		   SphP[pindex].Entropy / GAMMA_MINUS1 * pow(SphP[pindex].d.Density * a3inv, GAMMA_MINUS1));
	    n++;
	  }
      break;

    case IO_RHO:		/* density */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].d.Density;
	    n++;
	  }
      break;

    case IO_NE:		/* electron abundance */
#ifdef COOLING
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Ne;
	    n++;
	  }
#endif
      break;

    case IO_NH:		/* neutral hydrogen fraction */
#ifdef COOLING
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    ne = SphP[pindex].Ne;

	    AbundanceRatios(DMAX(All.MinEgySpec,
				 SphP[pindex].Entropy / GAMMA_MINUS1 * pow(SphP[pindex].d.Density *
									   a3inv,
									   GAMMA_MINUS1)),
			    SphP[pindex].d.Density * a3inv, &ne, &nh0, &nHeII);

	    *fp++ = nh0;
	    n++;
	  }
#endif
      break;

#ifdef CHEMISTRY
    case IO_ELECT:		/* electron abundance */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].elec;
	    n++;
	  }
      break;

    case IO_HI:		/* neutral hydrogen abundance */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].HI;
	    n++;
	  }
      break;

    case IO_HII:		/* ionized hydrogen abundance */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].HII;
	    n++;
	  }
      break;

    case IO_HeI:		/* neutral Helium */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].HeI;
	    n++;
	  }
      break;

    case IO_HeII:		/* ionized Heluum */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].HeII;
	    n++;
	  }
      break;

    case IO_HeIII:		/* double ionised Helium */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].HeIII;
	    n++;
	  }
      break;

    case IO_H2I:		/* H2 molecule */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].H2I;
	    n++;
	  }
      break;

    case IO_H2II:		/* ionised H2 molecule */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].H2II;
	    n++;
	  }
      break;

    case IO_HM:		/* H minus */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].HM;
	    n++;
	  }
      break;
#else
    case IO_ELECT:
    case IO_HI:
    case IO_HII:
    case IO_HeI:
    case IO_HeII:
    case IO_HeIII:
    case IO_H2I:
    case IO_H2II:
    case IO_HM:
      break;
#endif

    case IO_HSML:		/* SPH smoothing length */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = PPP[pindex].Hsml;
	    n++;
	  }
      break;

    case IO_SFR:		/* star formation rate */
#ifdef SFR
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = get_starformation_rate(pindex);
	    n++;
	  }
#endif
      break;

    case IO_DUST:
#ifdef DUST
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].DustRadius;
	    //*fp++ = P[pindex].DustVcoll;
	    n++;

	    /* if (n%100 == 0){ */
	    /*   printf("Dust Radius %d %g %g %d %d\n",pindex,P[pindex].DustRadius,P[pindex].Mass,P[pindex].Type,type);  */
	    /*   fflush(stdout);  */
	    /* } */
	  }
#endif
      break;

    case IO_TWO_DUST:
#ifdef  DUST_TWO_POPULATIONS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].MicroDustMass;
	    n++;

	    /* if (n%100 == 0){ */
	    /*   printf("Dust Radius %d %g %g %d %d\n",pindex,P[pindex].DustRadius,P[pindex].Mass,P[pindex].Type,type);  */
	    /*   fflush(stdout);  */
	    /* } */
	  }
      break;
#endif

    case IO_VCOLL_DUST:
#ifdef DUST_REAL_PEBBLE_COLLISIONS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].DustVcoll;
	    n++;

	    /*	    if (n%3000 == 0 && P[pindex].Type == 2){
	      printf("Index %d Radius %g  Dust Vcoll== %g, dln(a)/dt= %g \n",pindex,P[pindex].DustRadius, P[pindex].DustVcoll, 
		     P[pindex].LogDustRadius_by_dt);
	      fflush(stdout);
	      }*/
	  }
      break;
#endif

    case IO_TGROW_DUST:
#ifdef DUST_REAL_PEBBLE_COLLISIONS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].LogDustRadius_by_dt; /* Grain growth time scale in this case*/
	    n++;

	    /* if (n%3000 == 0 && P[pindex].Type == 2){ */
	    /*   printf("Index %d Radius %g  Dust Vcoll== %g, dln(a)/dt= %g \n",pindex,P[pindex].DustRadius, P[pindex].DustVcoll,  */
	    /* 	     P[pindex].LogDustRadius_by_dt); */
	    /*   fflush(stdout); */
	    /* }*/
	  }
      break;
#endif

    case IO_AGE:		/* stellar formation time */
#ifdef STELLARAGE
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].StellarAge;
	    n++;
	  }
#endif
      break;

    case IO_Z:			/* gas and star metallicity */
#ifdef METALS
#ifndef SFR_METALS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].Metallicity;
	    n++;
	  }
#else
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    if(P[pindex].Zm[0] <= 0 || P[pindex].Zm[6] <= 0)
	      printf("io2 H=%7.3e, He=%7.3e\n", P[pindex].Zm[6], P[pindex].Zm[0]);

	    for(k = 0; k < 12; k++)
	      *fp++ = P[pindex].Zm[k];
	    n++;
	  }
#endif
#endif
      break;

    case IO_POT:		/* gravitational potential */
#ifdef OUTPUTPOTENTIAL
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].p.Potential;
	    n++;
	  }
#endif
      break;

    case IO_ACCEL:		/* acceleration */
#ifdef OUTPUTACCELERATION
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      fp[k] = fac1 * P[pindex].g.GravAccel[k];
#ifdef PMGRID
	    for(k = 0; k < 3; k++)
	      fp[k] += fac1 * P[pindex].GravPM[k];
#endif
	    if(P[pindex].Type == 0)
	      for(k = 0; k < 3; k++)
              {
		fp[k] += fac2 * SphP[pindex].a.HydroAccel[k];
#ifdef RAD_ACCEL
	        fp[k] += fac2 * SphP[pindex].ra.RadAccel[k];
#endif
              }
	    fp += 3;
	    n++;
	  }
#endif
      break;

    case IO_DTENTR:		/* rate of change of entropy */
#ifdef OUTPUTCHANGEOFENTROPY
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].e.DtEntropy;
	    n++;
	  }
#endif
      break;

    case IO_STRESSDIAG:	/* Diagonal components of viscous shear tensor */
#ifdef OUTPUTSTRESS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      fp[k] = SphP[pindex].u.s.StressDiag[k];
	    fp += 3;
	    n++;
	  }
#endif
      break;

    case IO_STRESSOFFDIAG:	/* Offdiagonal components of viscous shear tensor */
#ifdef OUTPUTSTRESS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      fp[k] = SphP[pindex].u.s.StressOffDiag[k];
	    fp += 3;
	    n++;
	  }
#endif
      break;

    case IO_STRESSBULK:	/* Viscous bulk tensor */
#ifdef OUTPUTBULKSTRESS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].u.s.StressBulk;
	    n++;
	  }
#endif
      break;

    case IO_SHEARCOEFF:	/* Shear viscosity coefficient */
#ifdef OUTPUTSHEARCOEFF
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = get_shear_viscosity(pindex) *
	      pow((SphP[pindex].Entropy * pow(SphP[pindex].d.Density * a3inv,
					      GAMMA_MINUS1) / GAMMA_MINUS1), 2.5);
	    n++;
	  }
#endif
      break;

    case IO_TSTP:		/* timestep  */
#ifdef OUTPUTTIMESTEP

      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = (P[pindex].TimeBin ? (1 << P[pindex].TimeBin) : 0) * All.Timebase_interval;
	    n++;
	  }
#endif
      break;

    case IO_BFLD:		/* magnetic field  */
#ifdef MAGNETIC
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      *fp++ = SphP[pindex].B[k];
	    n++;
	  }
#endif
      break;

    case IO_DBDT:		/* rate of change of magnetic field  */
#ifdef DBOUTPUT
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      *fp++ = SphP[pindex].DtB[k];
	    n++;
	  }
#endif
      break;

    case IO_DIVB:		/* divergence of magnetic field  */
#ifdef TRACEDIVB
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].divB;
	    n++;
	  }
#endif
      break;

    case IO_ABVC:		/* artificial viscosity of particle  */
#ifdef TIME_DEP_ART_VISC
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].alpha;
	    n++;
	  }
#endif
      break;

    case IO_AMDC:		/* artificial viscosity of particle  */
#ifdef TIME_DEP_MAGN_DISP
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Balpha;
	    n++;
	  }
#endif
      break;

    case IO_PHI:		/* divBcleaning fuction of particle  */
#ifdef DIVBCLEANING_DEDNER
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].PhiPred;
	    n++;
	  }
#endif
      break;

    case IO_COOLRATE:		/* current cooling rate of particle  */
#ifdef OUTPUTCOOLRATE
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    ne = SphP[pindex].Ne;

	    /* get cooling time */
	    u = SphP[pindex].Entropy / GAMMA_MINUS1 * pow(SphP[pindex].d.Density * a3inv, GAMMA_MINUS1);

	    tcool = GetCoolingTime(u, SphP[pindex].d.Density * a3inv, &ne);

	    /* convert cooling time with current thermal energy to du/dt */
	    if(tcool != 0)
	      *fp++ = u / tcool;
	    else
	      *fp++ = 0;
	    n++;
	  }
#endif
      break;

    case IO_CONDRATE:		/* current heating/cooling due to thermal conduction  */
      break;

    case IO_BSMTH:		/* smoothed magnetic field */
#ifdef OUTPUTBSMOOTH
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      *fp++ = SphP[pindex].BSmooth[k];
	    n++;
	  }
#endif
      break;

    case IO_DENN:		/* density normalization factor */
#ifdef OUTPUTDENSNORM
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].DensityNorm;
	    n++;
	  }
#endif
      break;

    case IO_CR_C0:
#ifdef COSMIC_RAYS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].CR_C0;
	    n++;
	  }
#endif
      break;

    case IO_CR_Q0:
#ifdef COSMIC_RAYS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].CR_q0;
	    n++;
	  }
#endif
      break;

    case IO_CR_P0:
#ifdef COSMIC_RAYS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = CR_Physical_Pressure(&SphP[pindex]);
	    n++;
	  }
#endif
      break;

    case IO_CR_E0:
#ifdef COSMIC_RAYS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].CR_E0;
	    n++;
	  }
#endif
      break;

    case IO_CR_n0:
#ifdef COSMIC_RAYS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].CR_n0;
	    n++;
	  }
#endif
      break;

    case IO_CR_ThermalizationTime:
#ifdef COSMIC_RAYS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ =
	      CR_Tab_GetThermalizationTimescale(SphP[pindex].CR_q0 *
						pow(SphP[pindex].d.Density * a3inv, 0.333333),
						SphP[pindex].d.Density * a3inv);
	    n++;
	  }
#endif
      break;

    case IO_CR_DissipationTime:
#ifdef COSMIC_RAYS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ =
	      CR_Tab_GetDissipationTimescale(SphP[pindex].CR_q0 *
					     pow(SphP[pindex].d.Density * a3inv, 0.333333),
					     SphP[pindex].d.Density * a3inv);
	    n++;
	  }
#endif
      break;

    case IO_BHMASS:
#ifdef BLACK_HOLES
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].BH_Mass;
	    n++;
	  }
#endif
      break;

    case IO_BHMDOT:
#ifdef BLACK_HOLES
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].BH_Mdot;
	    n++;
	  }
#endif
      break;

    case IO_ACCDISCMASS:
#ifdef BLACK_HOLES
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].AccDisc_Mass;
	    n++;
	  }
#endif
      break;

    case IO_BHMSEED:
#ifdef BLACK_HOLES
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].BH_Mseed;
	    n++;
	  }
#endif
      break;

    case IO_NEWDENS:
#ifdef VIRTUAL
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].NewDensity;
	    n++;
	  }
#endif
      break;

    case IO_DELPHOTONMOM:
#ifdef VIRTUAL
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].DeltaPhotonMomentum;
	    n++;
	  }
#endif
      break;

    case IO_OLDPHOTONMOM:
#ifdef VIRTUAL
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].OldPhotonMomentum;
	    n++;
	  }
#endif
      break;

    case IO_BIPHOTONEN:
#ifdef VIRTUAL
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].BirthPhotonEnergy;
	    n++;
	  }
#endif
      break;

    case IO_DELPHOTONEN:
#if defined(VIRTUAL_HEATING) && defined(ENERGY_IN_DENSITY)
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].DeltaPhotonEnergy;
	    n++;
	  }
#endif
      break;

    case IO_ACCMASS:
#ifdef BLACK_HOLES
#ifdef VIRTUAL
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].accreted_mass;
	    n++;
	  }
#endif
#endif
      break;

    case IO_ACCMOM:
#ifdef BLACK_HOLES
#ifdef VIRTUAL
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
		*fp++ = SphP[pindex].accreted_momentum[k];
	    n++;
	  }
#endif
#endif
      break;
 
    case IO_INJVIRTUALEN:
#ifdef VIRTUAL_HEATING
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Injected_VIRTUAL_Energy;
	    n++;
	  }
#endif
      break;

    case IO_EMVIRTUALEN:
#ifdef VIRTUAL_HEATING
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Emitted_VIRTUAL_Energy;
	    n++;
	  }
#endif
      break;

    case IO_RET_E_FRAC:
#ifdef VIRTUAL_HEATING
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Returned_E_Fraction;
	    n++;
	  }
#endif
      break;

    case IO_MACH:
#ifdef MACHNUM
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Shock_MachNumber;
	    n++;
	  }
#endif
      break;

    case IO_DTENERGY:
#ifdef MACHSTATISTIC
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Shock_DtEnergy;
	    n++;
	  }
#endif
      break;

    case IO_PRESHOCK_DENSITY:
#ifdef CR_OUTPUT_JUMP_CONDITIONS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].PreShock_PhysicalDensity;
	    n++;
	  }
#endif
      break;

    case IO_PRESHOCK_ENERGY:
#ifdef CR_OUTPUT_JUMP_CONDITIONS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].PreShock_PhysicalEnergy;
	    n++;
	  }
#endif
      break;

    case IO_PRESHOCK_XCR:
#ifdef CR_OUTPUT_JUMP_CONDITIONS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].PreShock_XCR;
	    n++;
	  }
#endif
      break;

    case IO_DENSITY_JUMP:
#ifdef CR_OUTPUT_JUMP_CONDITIONS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Shock_DensityJump;
	    n++;
	  }
#endif
      break;

    case IO_ENERGY_JUMP:
#ifdef CR_OUTPUT_JUMP_CONDITIONS
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Shock_EnergyJump;
	    n++;
	  }
#endif
      break;

    case IO_CRINJECT:
#if defined( COSMIC_RAYS ) && defined( CR_OUTPUT_INJECTION )
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].CR_Specific_SupernovaHeatingRate;
	    n++;
	  }
#endif
      break;

    case IO_TIDALTENSOR:
#ifdef OUTPUT_TIDALTENSOR
      for(n = 0; n < pc; pindex++)

	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 6; k++)
	      {
		fp[k] = P[pindex].tite.tidal_tensor[k];
	      }
	    n++;
	    fp += 6;
	  }
#endif
      break;

    case IO_DISTORTIONTENSOR:
#ifdef DISTORTIONTENSOR
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 9; k++)
	      {
		fp[k] = P[pindex].distortion_tensor[k];
	      }
	    n++;
	    fp += 9;
	  }
#endif
      break;


    case IO_LASTENTRY:
      endrun(213);
      break;
    }

  *startindex = pindex;
}




/*! This function tells the size of one data entry in each of the blocks
 *  defined for the output file.
 */
int get_bytes_per_blockelement(enum iofields blocknr, int mode)
{
  int bytes_per_blockelement = 0;

  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_BSMTH:
    case IO_ACCEL:
    case IO_BFLD:
    case IO_DBDT:
    case IO_STRESSDIAG:
    case IO_STRESSOFFDIAG:
      if(mode)
	bytes_per_blockelement = 3 * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = 3 * sizeof(MyOutputFloat);
      break;

    case IO_ID:
      bytes_per_blockelement = sizeof(MyIDType);
      break;

    case IO_MASS:
    case IO_U:
    case IO_RHO:
    case IO_NE:
    case IO_NH:
    case IO_ELECT:
    case IO_HI:
    case IO_HII:
    case IO_HeI:
    case IO_HeII:
    case IO_HeIII:
    case IO_H2I:
    case IO_H2II:
    case IO_HM:
    case IO_HSML:
    case IO_SFR:
    case IO_AGE:
    case IO_POT:
    case IO_DTENTR:
    case IO_STRESSBULK:
    case IO_SHEARCOEFF:
    case IO_TSTP:
    case IO_DIVB:
    case IO_ABVC:
    case IO_AMDC:
    case IO_PHI:
    case IO_COOLRATE:
    case IO_CONDRATE:
    case IO_DENN:
    case IO_CR_C0:
    case IO_CR_Q0:
    case IO_CR_P0:
    case IO_CR_E0:
    case IO_CR_n0:
    case IO_CR_ThermalizationTime:
    case IO_CR_DissipationTime:
    case IO_BHMASS:
    case IO_BHMDOT:
    case IO_ACCDISCMASS:
    case IO_BHMSEED:   
    case IO_NEWDENS:
    case IO_DELPHOTONMOM:
    case IO_OLDPHOTONMOM:   
    case IO_BIPHOTONEN:
    case IO_DELPHOTONEN:
    case IO_ACCMASS:   
    case IO_ACCMOM:
    case IO_INJVIRTUALEN:
    case IO_EMVIRTUALEN:   
    case IO_RET_E_FRAC:
    case IO_MACH:
    case IO_DTENERGY:
    case IO_PRESHOCK_DENSITY:
    case IO_PRESHOCK_ENERGY:
    case IO_PRESHOCK_XCR:
    case IO_DENSITY_JUMP:
    case IO_ENERGY_JUMP:
    case IO_CRINJECT:
      if(mode)
	bytes_per_blockelement = sizeof(MyInputFloat);
      else
	bytes_per_blockelement = sizeof(MyOutputFloat);
      break;

    case IO_Z:
      if(mode)
	bytes_per_blockelement = sizeof(MyInputFloat);
      else
	bytes_per_blockelement = sizeof(MyOutputFloat);
      break;

    case IO_TIDALTENSOR:
      if(mode)
	bytes_per_blockelement = 6 * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = 6 * sizeof(MyOutputFloat);
      break;
    case IO_DISTORTIONTENSOR:
      if(mode)
	bytes_per_blockelement = 9 * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = 9 * sizeof(MyOutputFloat);
      break;

    case IO_DUST:
      if(mode)
      	bytes_per_blockelement = sizeof(MyInputFloat);
      else
      	bytes_per_blockelement = sizeof(MyOutputFloat);
      break;

    case IO_TWO_DUST:
      if(mode)
      	bytes_per_blockelement = sizeof(MyInputFloat);
      else
      	bytes_per_blockelement = sizeof(MyOutputFloat);
      break;

    case IO_VCOLL_DUST:
      if(mode)
      	bytes_per_blockelement = sizeof(MyInputFloat);
      else
      	bytes_per_blockelement = sizeof(MyOutputFloat);
      break;

    case IO_TGROW_DUST:
      if(mode)
      	bytes_per_blockelement = sizeof(MyInputFloat);
      else
      	bytes_per_blockelement = sizeof(MyOutputFloat);
      break;

    case IO_LASTENTRY:
      endrun(214);
      break;
    }

  return bytes_per_blockelement;
}

int get_datatype_in_block(enum iofields blocknr)
{
  int typekey;

  switch (blocknr)
    {
    case IO_ID:
#ifdef LONGIDS
      typekey = 2;		/* native long long */
#else
      typekey = 0;		/* native int */
#endif
      break;

    default:
      typekey = 1;		/* native MyOutputFloat */
      break;
    }

  return typekey;
}



int get_values_per_blockelement(enum iofields blocknr)
{
  int values = 0;

  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_BSMTH:
    case IO_ACCEL:
    case IO_BFLD:
    case IO_DBDT:
    case IO_STRESSDIAG:
    case IO_STRESSOFFDIAG:
      values = 3;
      break;

    case IO_ID:
    case IO_MASS:
    case IO_U:
    case IO_RHO:
    case IO_NE:
    case IO_NH:
    case IO_ELECT:
    case IO_HI:
    case IO_HII:
    case IO_HeI:
    case IO_HeII:
    case IO_HeIII:
    case IO_H2I:
    case IO_H2II:
    case IO_HM:
    case IO_HSML:
    case IO_SFR:
    case IO_AGE:
    case IO_POT:
    case IO_DTENTR:
    case IO_STRESSBULK:
    case IO_SHEARCOEFF:
    case IO_TSTP:
    case IO_DIVB:
    case IO_ABVC:
    case IO_AMDC:
    case IO_PHI:
    case IO_COOLRATE:
    case IO_CONDRATE:
    case IO_DENN:
    case IO_CR_C0:
    case IO_CR_Q0:
    case IO_CR_P0:
    case IO_CR_E0:
    case IO_CR_n0:
    case IO_CR_ThermalizationTime:
    case IO_CR_DissipationTime:
    case IO_BHMASS:
    case IO_BHMDOT:
    case IO_ACCDISCMASS:
    case IO_BHMSEED:   
    case IO_NEWDENS:
    case IO_DELPHOTONMOM:
    case IO_OLDPHOTONMOM:   
    case IO_BIPHOTONEN:
    case IO_DELPHOTONEN:
    case IO_ACCMASS:   
    case IO_ACCMOM:
    case IO_INJVIRTUALEN:
    case IO_EMVIRTUALEN:   
    case IO_RET_E_FRAC:
    case IO_MACH:
    case IO_DTENERGY:
    case IO_PRESHOCK_DENSITY:
    case IO_PRESHOCK_ENERGY:
    case IO_PRESHOCK_XCR:
    case IO_DENSITY_JUMP:
    case IO_ENERGY_JUMP:
    case IO_CRINJECT:
    case IO_DUST:
    case IO_TWO_DUST:
    case IO_VCOLL_DUST:
    case IO_TGROW_DUST:
      values = 1;
      break;

    case IO_Z:
      values = 1;
      break;

    case IO_TIDALTENSOR:
      values = 6;
      break;
    case IO_DISTORTIONTENSOR:
      values = 9;
      break;

    case IO_LASTENTRY:
      endrun(215);
      break;
    }
  return values;
}




/*! This function determines how many particles there are in a given block,
 *  based on the information in the header-structure.  It also flags particle
 *  types that are present in the block in the typelist array.
 */
int get_particles_in_block(enum iofields blocknr, int *typelist)
{
  int i, nall, ntot_withmasses, ngas, nstars, ndust;

  nall = 0;
  ntot_withmasses = 0;

  for(i = 0; i < 6; i++)
    {
      typelist[i] = 0;

      if(header.npart[i] > 0)
	{
	  nall += header.npart[i];
	  typelist[i] = 1;
	}

      if(All.MassTable[i] == 0)
	ntot_withmasses += header.npart[i];
    }

  ngas = header.npart[0];
  ndust = header.npart[2];
  nstars = header.npart[4];


  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ACCEL:
    case IO_TSTP:
    case IO_ID:
    case IO_POT:
      return nall;
      break;

    case IO_MASS:
      for(i = 0; i < 6; i++)
	{
	  typelist[i] = 0;
	  if(All.MassTable[i] == 0 && header.npart[i] > 0)
	    typelist[i] = 1;
	}
      return ntot_withmasses;
      break;

    case IO_U:
    case IO_RHO:
    case IO_NE:
    case IO_NH:
    case IO_ELECT:
    case IO_HI:
    case IO_HII:
    case IO_HeI:
    case IO_HeII:
    case IO_HeIII:
    case IO_H2I:
    case IO_H2II:
    case IO_HM:
    case IO_HSML:
    case IO_SFR:
    case IO_DTENTR:
    case IO_STRESSDIAG:
    case IO_STRESSOFFDIAG:
    case IO_STRESSBULK:
    case IO_SHEARCOEFF:
    case IO_BSMTH:
    case IO_BFLD:
    case IO_DBDT:
    case IO_DIVB:
    case IO_ABVC:
    case IO_AMDC:
    case IO_PHI:
    case IO_COOLRATE:
    case IO_CONDRATE:
    case IO_DENN:
    case IO_CR_C0:
    case IO_CR_Q0:
    case IO_CR_P0:
    case IO_CR_E0:
    case IO_CR_n0:
    case IO_CR_ThermalizationTime:
    case IO_CR_DissipationTime:
    case IO_MACH:
    case IO_DTENERGY:
    case IO_PRESHOCK_DENSITY:
    case IO_PRESHOCK_ENERGY:
    case IO_PRESHOCK_XCR:
    case IO_DENSITY_JUMP:
    case IO_ENERGY_JUMP:
    case IO_CRINJECT:
      for(i = 1; i < 6; i++)
	typelist[i] = 0;
      return ngas;
      break;

    case IO_AGE:
      for(i = 0; i < 6; i++)
	if(i != 4)
	  typelist[i] = 0;
      return nstars;

    case IO_Z:
      for(i = 0; i < 6; i++)
	if(i != 0 && i != 4)
	  typelist[i] = 0;
      return ngas + nstars;
      break;

    case IO_BHMASS:
    case IO_BHMDOT:
      for(i = 0; i < 6; i++)
	if(i != 5)
	  typelist[i] = 0;
      return header.npart[5];
      break;

    case IO_ACCDISCMASS:
    case IO_BHMSEED:
      for(i = 0; i < 6; i++)
	if(i != 5)
	  typelist[i] = 0;
      return header.npart[5];
      break;
 
    case IO_NEWDENS:
    case IO_DELPHOTONMOM:
    case IO_OLDPHOTONMOM:
    case IO_BIPHOTONEN:
    case IO_DELPHOTONEN:
      for(i = 0; i < 6; i++)
	if(i != 3)
	  typelist[i] = 0;
      return header.npart[3];
      break;

    case IO_ACCMASS:
    case IO_ACCMOM:
      for(i = 0; i < 6; i++)
	if(i != 0)
	  typelist[i] = 0;
      return ngas;
      break;

    case IO_INJVIRTUALEN:
    case IO_EMVIRTUALEN:
    case IO_RET_E_FRAC:
      for(i = 0; i < 6; i++)
	if(i != 0)
	  typelist[i] = 0;
      return ngas;
      break;

    case IO_TIDALTENSOR:
    case IO_DISTORTIONTENSOR:
      return nall;
      break;

    case IO_DUST:
    case IO_TWO_DUST:
    case IO_VCOLL_DUST:
    case IO_TGROW_DUST:
      return nall;
      break;

    case IO_LASTENTRY:
      endrun(216);
      break;
    }

  endrun(212);
  return 0;
}



/*! This function tells whether a block in the output file is present or not.
 */
int blockpresent(enum iofields blocknr)
{
  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ID:
    case IO_MASS:
    case IO_U:
    case IO_RHO:
    case IO_HSML:
      return 1;			/* always present */

    case IO_NE:
    case IO_NH:
      if(All.CoolingOn == 0)
	return 0;
      else
	return 1;
      break;


    case IO_SFR:
    case IO_AGE:
    case IO_Z:
      if(All.StarformationOn == 0)
	return 0;
      else
	{
#ifdef SFR
	  if(blocknr == IO_SFR)
	    return 1;
#endif
#ifdef STELLARAGE
	  if(blocknr == IO_AGE)
	    return 1;
#endif
#ifdef METALS
	  if(blocknr == IO_Z)
	    return 1;
#endif
	}
      return 0;
      break;


    case IO_ELECT:
    case IO_HI:
    case IO_HII:
    case IO_HeI:
    case IO_HeII:
    case IO_HeIII:
    case IO_H2I:
    case IO_H2II:
    case IO_HM:
#ifdef CHEMISTRY
      return 1;
#else
      return 0;
#endif
      break;


    case IO_POT:
#ifdef OUTPUTPOTENTIAL
      return 1;
#else
      return 0;
#endif

    case IO_ACCEL:
#ifdef OUTPUTACCELERATION
      return 1;
#else
      return 0;
#endif
      break;


    case IO_DTENTR:
#ifdef OUTPUTCHANGEOFENTROPY
      return 1;
#else
      return 0;
#endif
      break;

    case IO_STRESSDIAG:
#ifdef OUTPUTSTRESS
      return 1;
#else
      return 0;
#endif
      break;

    case IO_STRESSOFFDIAG:
#ifdef OUTPUTSTRESS
      return 1;
#else
      return 0;
#endif
      break;

    case IO_STRESSBULK:
#ifdef OUTPUTBULKSTRESS
      return 1;
#else
      return 0;
#endif
      break;

    case IO_SHEARCOEFF:
#ifdef OUTPUTSHEARCOEFF
      return 1;
#else
      return 0;
#endif
      break;

    case IO_TSTP:
#ifdef IO_TSTP
      return 1;
#else
      return 0;
#endif


    case IO_BFLD:
#ifdef MAGNETIC
      return 1;
#else
      return 0;
#endif
      break;


    case IO_DBDT:
#ifdef DBOUTPUT
      return 1;
#else
      return 0;
#endif
      break;


    case IO_DIVB:
#ifdef TRACEDIVB
      return 1;
#else
      return 0;
#endif
      break;


    case IO_ABVC:
#ifdef TIME_DEP_ART_VISC
      return 1;
#else
      return 0;
#endif
      break;


    case IO_AMDC:
#ifdef TIME_DEP_MAGN_DISP
      return 1;
#else
      return 0;
#endif
      break;


    case IO_PHI:
#ifdef DIVBCLEANING_DEDNER
      return 1;
#else
      return 0;
#endif
      break;


    case IO_COOLRATE:
#ifdef OUTPUTCOOLRATE
      return 1;
#else
      return 0;
#endif
      break;


    case IO_CONDRATE:
#ifdef OUTPUTCONDRATE
      return 1;
#else
      return 0;
#endif
      break;


    case IO_BSMTH:
#ifdef OUTPUTBSMOOTH
      return 1;
#else
      return 0;
#endif
      break;


    case IO_DENN:
#ifdef OUTPUTDENSNORM
      return 1;
#else
      return 0;
#endif
      break;


    case IO_BHMASS:
    case IO_BHMDOT:
#ifdef BLACK_HOLES
      return 1;
#else
      return 0;
#endif
      break;

 case IO_ACCDISCMASS:
    case IO_BHMSEED:
#ifdef BLACK_HOLES
      return 1;
#else
      return 0;
#endif
      break;

    case IO_NEWDENS:
#ifdef VIRTUAL
      return 1;
#else
      return 0;
#endif
      break;

    case IO_DELPHOTONMOM:
    case IO_OLDPHOTONMOM:
    case IO_BIPHOTONEN:
#ifdef BLACK_HOLES
      return 1;
#else
      return 0;
#endif
      break;

    case IO_DELPHOTONEN:
#if defined(VIRTUAL_HEATING) && defined(ENERGY_IN_DENSITY)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_ACCMASS:
    case IO_ACCMOM:
#ifdef BLACK_HOLES
#ifdef VIRTUAL
      return 1;
#else
      return 0;
#endif
#endif
      break;

    case IO_INJVIRTUALEN:
    case IO_EMVIRTUALEN:
    case IO_RET_E_FRAC:
#ifdef VIRTUAL_HEATING
      return 1;
#else
      return 0;
#endif
      break;

    case IO_MACH:
#ifdef MACHNUM
      return 1;
#else
      return 0;
#endif
      break;


    case IO_DTENERGY:
#ifdef MACHSTATISTIC
      return 1;
#else
      return 0;
#endif
      break;


    case IO_CR_C0:
    case IO_CR_Q0:
    case IO_CR_P0:
    case IO_CR_E0:
    case IO_CR_n0:
    case IO_CR_ThermalizationTime:
    case IO_CR_DissipationTime:
#ifdef COSMIC_RAYS
      return 1;
#else
      return 0;
#endif
      break;


    case IO_PRESHOCK_DENSITY:
    case IO_PRESHOCK_ENERGY:
    case IO_PRESHOCK_XCR:
    case IO_DENSITY_JUMP:
    case IO_ENERGY_JUMP:
#ifdef CR_OUTPUT_JUMP_CONDITIONS
      return 1;
#else
      return 0;
#endif
      break;


    case IO_CRINJECT:
#ifdef CR_OUTPUT_INJECTION
      return 1;
#else
      return 0;
#endif
      break;

    case IO_TIDALTENSOR:
#ifdef OUTPUT_TIDALTENSOR
      return 1;
#else
      return 0;
#endif
    case IO_DISTORTIONTENSOR:
#ifdef DISTORTIONTENSOR
      return 1;
#else
      return 0;
#endif
      break;

    case IO_DUST: 
#ifdef DUST
      //      if(blocknr == IO_DUST)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_TWO_DUST: 
#ifdef DUST_TWO_POPULATIONS
      //      if(blocknr == IO_DUST)
      return 1;
#else
      return 0;
#endif
      break;

    case IO_VCOLL_DUST: 
#ifdef DUST_REAL_PEBBLE_COLLISIONS
      return 1;
#else
      return 0;
#endif
      break;

    case IO_TGROW_DUST:
#ifdef DUST_REAL_PEBBLE_COLLISIONS
      return 1;
#else
      return 0;
#endif
      break;


    case IO_LASTENTRY:
      return 0;			/* will not occur */
      break;
    }

  return 0;			/* default: not present */
}




/*! This function associates a short 4-character block name with each block number.
 *  This is stored in front of each block for snapshot FileFormat=2.
 */
void get_Tab_IO_Label(enum iofields blocknr, char *label)
{
  switch (blocknr)
    {
    case IO_POS:
      strncpy(label, "POS ", 4);
      break;
    case IO_VEL:
      strncpy(label, "VEL ", 4);
      break;
    case IO_ID:
      strncpy(label, "ID  ", 4);
      break;
    case IO_MASS:
      strncpy(label, "MASS", 4);
      break;
    case IO_U:
      strncpy(label, "U   ", 4);
      break;
    case IO_RHO:
      strncpy(label, "RHO ", 4);
      break;
    case IO_NE:
      strncpy(label, "NE  ", 4);
      break;
    case IO_NH:
      strncpy(label, "NH  ", 4);
      break;
    case IO_ELECT:
      strncpy(label, "elect ", 4);
      break;
    case IO_HI:
      strncpy(label, "HI ", 4);
      break;
    case IO_HII:
      strncpy(label, "HII ", 4);
      break;
    case IO_HeI:
      strncpy(label, "HeI ", 4);
      break;
    case IO_HeII:
      strncpy(label, "HeII ", 4);
      break;
    case IO_HeIII:
      strncpy(label, "HeIII ", 4);
      break;
    case IO_H2I:
      strncpy(label, "H2I ", 4);
      break;
    case IO_H2II:
      strncpy(label, "H2II ", 4);
      break;
    case IO_HM:
      strncpy(label, "HM ", 4);
      break;
    case IO_HSML:
      strncpy(label, "HSML", 4);
      break;
    case IO_SFR:
      strncpy(label, "SFR ", 4);
      break;
    case IO_AGE:
      strncpy(label, "AGE ", 4);
      break;
    case IO_Z:
      strncpy(label, "Z   ", 4);
      break;
    case IO_POT:
      strncpy(label, "POT ", 4);
      break;
    case IO_ACCEL:
      strncpy(label, "ACCE", 4);
      break;
    case IO_DTENTR:
      strncpy(label, "ENDT", 4);
      break;
    case IO_STRESSDIAG:
      strncpy(label, "STRD", 4);
      break;
    case IO_STRESSOFFDIAG:
      strncpy(label, "STRO", 4);
      break;
    case IO_STRESSBULK:
      strncpy(label, "STRB", 4);
      break;
    case IO_SHEARCOEFF:
      strncpy(label, "SHCO", 4);
      break;
    case IO_TSTP:
      strncpy(label, "TSTP", 4);
      break;
    case IO_BFLD:
      strncpy(label, "BFLD", 4);
      break;
    case IO_DBDT:
      strncpy(label, "DBDT", 4);
      break;
    case IO_DIVB:
      strncpy(label, "DIVB", 4);
      break;
    case IO_ABVC:
      strncpy(label, "ABVC", 4);
      break;
    case IO_AMDC:
      strncpy(label, "AMDC", 4);
      break;
    case IO_PHI:
      strncpy(label, "PHI ", 4);
      break;
    case IO_COOLRATE:
      strncpy(label, "COOR", 4);
      break;
    case IO_CONDRATE:
      strncpy(label, "CONR", 4);
      break;
    case IO_BSMTH:
      strncpy(label, "BFSM", 4);
      break;
    case IO_DENN:
      strncpy(label, "DENN", 4);
      break;
    case IO_CR_C0:
      strncpy(label, "CRC0", 4);
      break;
    case IO_CR_Q0:
      strncpy(label, "CRQ0", 4);
      break;
    case IO_CR_P0:
      strncpy(label, "CRP0", 4);
      break;
    case IO_CR_E0:
      strncpy(label, "CRE0", 4);
      break;
    case IO_CR_n0:
      strncpy(label, "CRn0", 4);
      break;
    case IO_CR_ThermalizationTime:
      strncpy(label, "CRco", 4);
      break;
    case IO_CR_DissipationTime:
      strncpy(label, "CRdi", 4);
      break;
    case IO_BHMASS:
      strncpy(label, "BHMA", 4);
      break;
    case IO_BHMDOT:
      strncpy(label, "BHMD", 4);
      break;
    case IO_ACCDISCMASS:
      strncpy(label, "ACDM", 4);
      break;
    case IO_BHMSEED:
      strncpy(label, "BHMS", 4);
      break;
    case IO_NEWDENS:
      strncpy(label, "NEWD", 4);
      break;
    case IO_DELPHOTONMOM:
      strncpy(label, "DPhM", 4);
      break;
    case IO_OLDPHOTONMOM:
      strncpy(label, "OPhM", 4);
      break;
    case IO_BIPHOTONEN:
      strncpy(label, "BPhE", 4);
      break;
    case IO_DELPHOTONEN:
      strncpy(label, "DPhE", 4);
      break;
    case IO_ACCMASS:
      strncpy(label, "ACMA", 4);
      break;
    case IO_ACCMOM:
      strncpy(label, "ACMO", 4);
      break;
    case IO_INJVIRTUALEN:
      strncpy(label, "INVE", 4);
      break;
    case IO_EMVIRTUALEN:
      strncpy(label, "EMVE", 4);
      break;
    case IO_RET_E_FRAC:
      strncpy(label, "REFR", 4);
      break;
    case IO_MACH:
      strncpy(label, "MACH", 4);
      break;
    case IO_DTENERGY:
      strncpy(label, "DTEG", 4);
      break;
    case IO_PRESHOCK_DENSITY:
      strncpy(label, "PSDE", 4);
      break;
    case IO_PRESHOCK_ENERGY:
      strncpy(label, "PSEN", 4);
      break;
    case IO_PRESHOCK_XCR:
      strncpy(label, "PSXC", 4);
      break;
    case IO_DENSITY_JUMP:
      strncpy(label, "DJMP", 4);
      break;
    case IO_ENERGY_JUMP:
      strncpy(label, "EJMP", 4);
      break;
    case IO_CRINJECT:
      strncpy(label, "CRDE", 4);
      break;
    case IO_TIDALTENSOR:
      strncpy(label, "TITE", 4);
      break;
    case IO_DISTORTIONTENSOR:
      strncpy(label, "DITE", 4);
      break;
    case IO_DUST: 
      strncpy(label, "DUST", 4);
      break;
    case IO_TWO_DUST: 
      strncpy(label, "DUS2", 4);
      break;
    case IO_VCOLL_DUST: 
      strncpy(label, "VCOL", 4);
      break;
    case IO_TGROW_DUST:
      strncpy(label, "TGRO", 4);
      break;

    case IO_LASTENTRY:
      endrun(217);
      break;
    }
}


void get_dataset_name(enum iofields blocknr, char *buf)
{
  strcpy(buf, "default");

  switch (blocknr)
    {
    case IO_POS:
      strcpy(buf, "Coordinates");
      break;
    case IO_VEL:
      strcpy(buf, "Velocities");
      break;
    case IO_ID:
      strcpy(buf, "ParticleIDs");
      break;
    case IO_MASS:
      strcpy(buf, "Masses");
      break;
    case IO_U:
      strcpy(buf, "InternalEnergy");
      break;
    case IO_RHO:
      strcpy(buf, "Density");
      break;
    case IO_NE:
      strcpy(buf, "ElectronAbundance");
      break;
    case IO_NH:
      strcpy(buf, "NeutralHydrogenAbundance");
      break;
    case IO_ELECT:
      strcpy(buf, "elect");
      break;
    case IO_HI:
      strcpy(buf, "HI");
      break;
    case IO_HII:
      strcpy(buf, "HII");
      break;
    case IO_HeI:
      strcpy(buf, "HeI");
      break;
    case IO_HeII:
      strcpy(buf, "HeII");
      break;
    case IO_HeIII:
      strcpy(buf, "HeIII");
      break;
    case IO_H2I:
      strcpy(buf, "H2I");
      break;
    case IO_H2II:
      strcpy(buf, "H2II");
      break;
    case IO_HM:
      strcpy(buf, "HM");
      break;
    case IO_HSML:
      strcpy(buf, "SmoothingLength");
      break;
    case IO_SFR:
      strcpy(buf, "StarFormationRate");
      break;
    case IO_AGE:
      strcpy(buf, "StellarFormationTime");
      break;
    case IO_Z:
      strcpy(buf, "Metallicity");
      break;
    case IO_POT:
      strcpy(buf, "Potential");
      break;
    case IO_ACCEL:
      strcpy(buf, "Acceleration");
      break;
    case IO_DTENTR:
      strcpy(buf, "RateOfChangeOfEntropy");
      break;
    case IO_STRESSDIAG:
      strcpy(buf, "DiagonalStressTensor");
      break;
    case IO_STRESSOFFDIAG:
      strcpy(buf, "OffDiagonalStressTensor");
      break;
    case IO_STRESSBULK:
      strcpy(buf, "BulkStressTensor");
      break;
    case IO_SHEARCOEFF:
      strcpy(buf, "ShearCoefficient");
      break;
    case IO_TSTP:
      strcpy(buf, "TimeStep");
      break;
    case IO_BFLD:
      strcpy(buf, "MagneticField");
      break;
    case IO_DBDT:
      strcpy(buf, "RateOfChangeOfMagneticField");
      break;
    case IO_DIVB:
      strcpy(buf, "DivergenceOfMagneticField");
      break;
    case IO_ABVC:
      strcpy(buf, "ArtificialViscosity");
      break;
    case IO_AMDC:
      strcpy(buf, "ArtificialMagneticDissipatio");
      break;
    case IO_PHI:
      strcpy(buf, "DivBcleaningFunctionPhi");
      break;
    case IO_COOLRATE:
      strcpy(buf, "CoolingRate");
      break;
    case IO_CONDRATE:
      strcpy(buf, "ConductionRate");
      break;
    case IO_BSMTH:
      strcpy(buf, "SmoothedMagneticField");
      break;
    case IO_DENN:
      strcpy(buf, "Denn");
      break;
    case IO_CR_C0:
      strcpy(buf, "CR_C0");
      break;
    case IO_CR_Q0:
      strcpy(buf, "CR_q0");
      break;
    case IO_CR_P0:
      strcpy(buf, "CR_P0");
      break;
    case IO_CR_E0:
      strcpy(buf, "CR_E0");
      break;
    case IO_CR_n0:
      strcpy(buf, "CR_n0");
      break;
    case IO_CR_ThermalizationTime:
      strcpy(buf, "CR_ThermalizationTime");
      break;
    case IO_CR_DissipationTime:
      strcpy(buf, "CR_DissipationTime");
      break;
    case IO_BHMASS:
      strcpy(buf, "BH_Mass");
      break;
    case IO_BHMDOT:
      strcpy(buf, "BH_Mdot");
      break;
    case IO_ACCDISCMASS:
      strcpy(buf, "AccDisc_Mass");
      break;
    case IO_BHMSEED:
      strcpy(buf, "BH_Mseed");
      break;
    case IO_NEWDENS:
      strcpy(buf, "NewDensity");
      break;
    case IO_DELPHOTONMOM:
      strcpy(buf, "DeltaPhotonMomentum");
      break;
    case IO_OLDPHOTONMOM:
      strcpy(buf, "OldPhotonMomentum");
      break;
    case IO_BIPHOTONEN:
      strcpy(buf, "BirthPhotonEnergy");
      break;
    case IO_DELPHOTONEN:
      strcpy(buf, "DeltaPhotonEnergy");
      break;
    case IO_ACCMASS:
      strcpy(buf, "accreted_mass");
      break;
    case IO_ACCMOM:
      strcpy(buf, "accreted_momentum");
      break;
    case IO_INJVIRTUALEN:
      strcpy(buf, "Injected_VIRTUAL_Energy");
      break;
    case IO_EMVIRTUALEN:
      strcpy(buf, "Emitted_VIRTUAL_Energy");
      break;
    case IO_RET_E_FRAC:
      strcpy(buf, "Returned_E_Fraction");
      break;
    case IO_MACH:
      strcpy(buf, "MachNumber");
      break;
    case IO_DTENERGY:
      strcpy(buf, "DtEnergy");
      break;
    case IO_PRESHOCK_DENSITY:
      strcpy(buf, "Preshock_Density");
      break;
    case IO_PRESHOCK_ENERGY:
      strcpy(buf, "Preshock_Energy");
      break;
    case IO_PRESHOCK_XCR:
      strcpy(buf, "Preshock_XCR");
      break;
    case IO_DENSITY_JUMP:
      strcpy(buf, "DensityJump");
      break;
    case IO_ENERGY_JUMP:
      strcpy(buf, "EnergyJump");
      break;
    case IO_CRINJECT:
      strcpy(buf, "CR_DtE");
      break;
    case IO_TIDALTENSOR:
      strcpy(buf, "TidalTensor");
      break;
    case IO_DISTORTIONTENSOR:
      strcpy(buf, "DistortionTensor");
      break;
    case IO_DUST:
      strcpy(buf, "DustRadius");
      break;
    case IO_TWO_DUST:
#ifdef DUST_TWO_POPULATIONS
      strcpy(buf, "MicroDustMass");
      break;
#endif

#ifdef DUST_REAL_PEBBLE_COLLISIONS
      strcpy(buf, "DustVcoll");
#endif
      break;

    case IO_LASTENTRY:
      endrun(218);
      break;
    }
}



/*! This function writes a snapshot file containing the data from processors
 *  'writeTask' to 'lastTask'. 'writeTask' is the one that actually writes.
 *  Each snapshot file contains a header first, then particle positions,
 *  velocities and ID's.  Then particle masses are written for those particle
 *  types with zero entry in MassTable.  After that, first the internal
 *  energies u, and then the density is written for the SPH particles.  If
 *  cooling is enabled, mean molecular weight and neutral hydrogen abundance
 *  are written for the gas particles. This is followed by the SPH smoothing
 *  length and further blocks of information, depending on included physics
 *  and compile-time flags.
 */
void write_file(char *fname, int writeTask, int lastTask)
{
  int type, bytes_per_blockelement, npart, nextblock, typelist[6];
  int n_for_this_task, ntask, n, p, pc, offset = 0, task;
  int blockmaxlen, ntot_type[6], nn[6];
  enum iofields blocknr;
  char label[8];
  int bnr;
  int blksize;
  MPI_Status status;
  FILE *fd = 0;

#ifdef HAVE_HDF5
  hid_t hdf5_file = 0, hdf5_grp[6], hdf5_headergrp = 0, hdf5_dataspace_memory, hdf5_paramgrp = 0;
  hid_t hdf5_datatype = 0, hdf5_dataspace_in_file = 0, hdf5_dataset = 0;
  herr_t hdf5_status;
  hsize_t dims[2], count[2], start[2];
  int rank = 0, pcsum = 0;
  char buf[500];
#endif

#define SKIP  {my_fwrite(&blksize,sizeof(int),1,fd);}

  /* determine particle numbers of each type in file */


  if(ThisTask == writeTask)
    {
      for(n = 0; n < 6; n++)
	ntot_type[n] = n_type[n];

      for(task = writeTask + 1; task <= lastTask; task++)
	{
	  MPI_Recv(&nn[0], 6, MPI_INT, task, TAG_LOCALN, MPI_COMM_WORLD, &status);
	  for(n = 0; n < 6; n++)
	    ntot_type[n] += nn[n];
	}

      for(task = writeTask + 1; task <= lastTask; task++)
	MPI_Send(&ntot_type[0], 6, MPI_INT, task, TAG_N, MPI_COMM_WORLD);
    }
  else
    {
      MPI_Send(&n_type[0], 6, MPI_INT, writeTask, TAG_LOCALN, MPI_COMM_WORLD);
      MPI_Recv(&ntot_type[0], 6, MPI_INT, writeTask, TAG_N, MPI_COMM_WORLD, &status);
    }



  /* fill file header */

  for(n = 0; n < 6; n++)
    {
      header.npart[n] = ntot_type[n];
      header.npartTotal[n] = (unsigned int) ntot_type_all[n];
      header.npartTotalHighWord[n] = (unsigned int) (ntot_type_all[n] >> 32);
    }

  for(n = 0; n < 6; n++)
    header.mass[n] = All.MassTable[n];

  header.time = All.Time;

  if(All.ComovingIntegrationOn)
    header.redshift = 1.0 / All.Time - 1;
  else
    header.redshift = 0;

  header.flag_sfr = 0;
  header.flag_feedback = 0;
  header.flag_cooling = 0;
  header.flag_stellarage = 0;
  header.flag_metals = 0;

#ifdef COOLING
  header.flag_cooling = 1;
#endif

#ifdef SFR
  header.flag_sfr = 1;
  header.flag_feedback = 1;
#ifdef STELLARAGE
  header.flag_stellarage = 1;
#endif
#ifdef METALS
  header.flag_metals = 1;
#endif
#endif

  header.num_files = All.NumFilesPerSnapshot;
  header.BoxSize = All.BoxSize;
  header.Omega0 = All.Omega0;
  header.OmegaLambda = All.OmegaLambda;
  header.HubbleParam = All.HubbleParam;

#ifdef OUTPUT_IN_DOUBLEPRECISION
  header.flag_doubleprecision = 1;
#else
  header.flag_doubleprecision = 0;
#endif

  /* open file and write header */

  if(ThisTask == writeTask)
    {
      if(All.SnapFormat == 3)
	{
#ifdef HAVE_HDF5
	  sprintf(buf, "%s.hdf5", fname);
	  hdf5_file = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	  hdf5_headergrp = H5Gcreate(hdf5_file, "/Header", 0);

	  for(type = 0; type < 6; type++)
	    {
	      if(header.npart[type] > 0)
		{
		  sprintf(buf, "/PartType%d", type);
		  hdf5_grp[type] = H5Gcreate(hdf5_file, buf, 0);
		}
	    }

	  write_header_attributes_in_hdf5(hdf5_headergrp);
#endif
	}
      else
	{
	  if(!(fd = fopen(fname, "w")))
	    {
	      printf("can't open file `%s' for writing snapshot.\n", fname);
	      endrun(123);
	    }

	  if(All.SnapFormat == 2)
	    {
	      blksize = sizeof(int) + 4 * sizeof(char);
	      SKIP;
	      my_fwrite((void *) "HEAD", sizeof(char), 4, fd);
	      nextblock = sizeof(header) + 2 * sizeof(int);
	      my_fwrite(&nextblock, sizeof(int), 1, fd);
	      SKIP;
	    }

	  blksize = sizeof(header);
	  SKIP;
	  my_fwrite(&header, sizeof(header), 1, fd);
	  SKIP;
	}
    }

  ntask = lastTask - writeTask + 1;

  for(bnr = 0; bnr < 1000; bnr++)
    {
      blocknr = (enum iofields) bnr;

      if(blocknr == IO_LASTENTRY)
	break;

      if(blockpresent(blocknr))
	{
	  bytes_per_blockelement = get_bytes_per_blockelement(blocknr, 0);

	  blockmaxlen = ((int) (All.BufferSize * 1024 * 1024)) / bytes_per_blockelement;

	  npart = get_particles_in_block(blocknr, &typelist[0]);

	  if(npart > 0)
	    {
	      if(ThisTask == writeTask)
		{

		  if(All.SnapFormat == 1 || All.SnapFormat == 2)
		    {
		      if(All.SnapFormat == 2)
			{
			  blksize = sizeof(int) + 4 * sizeof(char);
			  SKIP;
			  get_Tab_IO_Label(blocknr, label);
			  my_fwrite(label, sizeof(char), 4, fd);
			  nextblock = npart * bytes_per_blockelement + 2 * sizeof(int);
			  my_fwrite(&nextblock, sizeof(int), 1, fd);
			  SKIP;
			}

		      blksize = npart * bytes_per_blockelement;
		      SKIP;

		    }
		}

	      for(type = 0; type < 6; type++)
		{
		  if(typelist[type])
		    {
#ifdef HAVE_HDF5
		      if(ThisTask == writeTask && All.SnapFormat == 3 && header.npart[type] > 0)
			{
			  switch (get_datatype_in_block(blocknr))
			    {
			    case 0:
			      hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT);
			      break;
			    case 1:
#ifdef OUTPUT_IN_DOUBLEPRECISION
			      hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
			      hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
			      break;
			    case 2:
			      hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT64);
			      break;
			    }

			  dims[0] = header.npart[type];
			  dims[1] = get_values_per_blockelement(blocknr);
			  if(dims[1] == 1)
			    rank = 1;
			  else
			    rank = 2;

			  get_dataset_name(blocknr, buf);

			  hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);
			  hdf5_dataset =
			    H5Dcreate(hdf5_grp[type], buf, hdf5_datatype, hdf5_dataspace_in_file,
				      H5P_DEFAULT);

			  pcsum = 0;
			}
#endif

		      for(task = writeTask, offset = 0; task <= lastTask; task++)
			{
			  if(task == ThisTask)
			    {
			      n_for_this_task = n_type[type];

			      for(p = writeTask; p <= lastTask; p++)
				if(p != ThisTask)
				  MPI_Send(&n_for_this_task, 1, MPI_INT, p, TAG_NFORTHISTASK, MPI_COMM_WORLD);
			    }
			  else
			    MPI_Recv(&n_for_this_task, 1, MPI_INT, task, TAG_NFORTHISTASK, MPI_COMM_WORLD,
				     &status);

			  while(n_for_this_task > 0)
			    {
			      pc = n_for_this_task;

			      if(pc > blockmaxlen)
				pc = blockmaxlen;

			      if(ThisTask == task)
				fill_write_buffer(blocknr, &offset, pc, type);

			      if(ThisTask == writeTask && task != writeTask)
				MPI_Recv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, task,
					 TAG_PDATA, MPI_COMM_WORLD, &status);

			      if(ThisTask != writeTask && task == ThisTask)
				MPI_Ssend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, writeTask,
					  TAG_PDATA, MPI_COMM_WORLD);

			      if(ThisTask == writeTask)
				{
				  if(All.SnapFormat == 3)
				    {
#ifdef HAVE_HDF5
				      start[0] = pcsum;
				      start[1] = 0;

				      count[0] = pc;
				      count[1] = get_values_per_blockelement(blocknr);
				      pcsum += pc;

				      H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET,
							  start, NULL, count, NULL);

				      dims[0] = pc;
				      dims[1] = get_values_per_blockelement(blocknr);
				      hdf5_dataspace_memory = H5Screate_simple(rank, dims, NULL);

				      hdf5_status =
					H5Dwrite(hdf5_dataset, hdf5_datatype,
						 hdf5_dataspace_memory,
						 hdf5_dataspace_in_file, H5P_DEFAULT, CommBuffer);

				      H5Sclose(hdf5_dataspace_memory);
#endif
				    }
				  else
				    {
				      my_fwrite(CommBuffer, bytes_per_blockelement, pc, fd);
				    }
				}

			      n_for_this_task -= pc;
			    }
			}

#ifdef HAVE_HDF5
		      if(ThisTask == writeTask && All.SnapFormat == 3 && header.npart[type] > 0)
			{
			  if(All.SnapFormat == 3)
			    {
			      H5Dclose(hdf5_dataset);
			      H5Sclose(hdf5_dataspace_in_file);
			      H5Tclose(hdf5_datatype);
			    }
			}
#endif
		    }
		}

	      if(ThisTask == writeTask)
		{
		  if(All.SnapFormat == 1 || All.SnapFormat == 2)
		    SKIP;
		}
	    }
	}
    }

  if(ThisTask == writeTask)
    {
      if(All.SnapFormat == 3)
	{
#ifdef HAVE_HDF5
	  for(type = 5; type >= 0; type--)
	    if(header.npart[type] > 0)
	      H5Gclose(hdf5_grp[type]);
	  H5Gclose(hdf5_headergrp);
	  H5Fclose(hdf5_file);
#endif
	}
      else
	fclose(fd);
    }
}




#ifdef HAVE_HDF5
void write_header_attributes_in_hdf5(hid_t handle)
{
  hsize_t adim[1] = { 6 };
  hid_t hdf5_dataspace, hdf5_attribute;

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "NumPart_ThisFile", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npart);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "NumPart_Total", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotal);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "NumPart_Total_HighWord", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotalHighWord);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);


  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "MassTable", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, header.mass);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Time", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Redshift", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.redshift);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "BoxSize", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.BoxSize);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "NumFilesPerSnapshot", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.num_files);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Omega0", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.Omega0);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "OmegaLambda", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.OmegaLambda);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "HubbleParam", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.HubbleParam);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Sfr", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_sfr);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Cooling", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_cooling);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_StellarAge", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_stellarage);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Metals", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_metals);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Feedback", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_feedback);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_DoublePrecision", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_doubleprecision);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);
}
#endif





/*! This catches I/O errors occuring for my_fwrite(). In this case we
 *  better stop.
 */
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nwritten;

  if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
    {
      printf("I/O error (fwrite) on task=%d has occured: %s\n", ThisTask, strerror(errno));
      fflush(stdout);
      endrun(777);
    }
  return nwritten;
}


/*! This catches I/O errors occuring for fread(). In this case we
 *  better stop.
 */
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nread;

  if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
    {
      if(feof(stream))
	printf("I/O error (fread) on task=%d has occured: end of file\n", ThisTask);
      else
	printf("I/O error (fread) on task=%d has occured: %s\n", ThisTask, strerror(errno));
      fflush(stdout);
      endrun(778);
    }
  return nread;
}



#ifdef ORDER_SNAPSHOTS_BY_ID
int io_compare_P_ID(const void *a, const void *b)
{
  if(((struct particle_data *) a)->ID < (((struct particle_data *) b)->ID))
    return -1;

  if(((struct particle_data *) a)->ID > (((struct particle_data *) b)->ID))
    return +1;

  return 0;
}

int io_compare_P_GrNr_SubNr(const void *a, const void *b)
{
  if(((struct particle_data *) a)->GrNr < (((struct particle_data *) b)->GrNr))
    return -1;

  if(((struct particle_data *) a)->GrNr > (((struct particle_data *) b)->GrNr))
    return +1;

  if(((struct particle_data *) a)->SubNr < (((struct particle_data *) b)->SubNr))
    return -1;

  if(((struct particle_data *) a)->SubNr > (((struct particle_data *) b)->SubNr))
    return +1;

  return 0;
}
#endif
