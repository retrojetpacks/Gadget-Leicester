#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "allvars.h"
#include "proto.h"
#ifdef COSMIC_RAYS
#include "cosmic_rays.h"
#endif


void reconstruct_timebins(void)
{
  int i, n, prev, bin;

  for(bin = 0; bin < TIMEBINS; bin++)
    {
      TimeBinCount[bin] = 0;
      TimeBinCountSph[bin] = 0;
      FirstInTimeBin[bin] = -1;
      LastInTimeBin[bin] = -1;
#ifdef SFR
      TimeBinSfr[bin] = 0;
#endif
#ifdef BLACK_HOLES
      TimeBin_BH_mass[bin] = 0;
      TimeBin_BH_dynamicalmass[bin] = 0;
      TimeBin_BH_Mdot[bin] = 0;
#endif
#ifdef DUST
      //      TimeBin_Dust_Mass[bin]=0;
#endif
    }


  for(i = 0; i < NumPart; i++)
    {
      bin = P[i].TimeBin;

      if(TimeBinCount[bin] > 0)
	{
	  PrevInTimeBin[i] = LastInTimeBin[bin];
	  NextInTimeBin[i] = -1;
	  NextInTimeBin[LastInTimeBin[bin]] = i;
	  LastInTimeBin[bin] = i;
	}
      else
	{
	  FirstInTimeBin[bin] = LastInTimeBin[bin] = i;
	  PrevInTimeBin[i] = NextInTimeBin[i] = -1;
	}
      TimeBinCount[bin]++;
      if(P[i].Type == 0)
	TimeBinCountSph[bin]++;

#ifdef SFR
      if(P[i].Type == 0)
	TimeBinSfr[bin] += SphP[i].Sfr;
#endif
#if BLACK_HOLES
      if(P[i].Type == 5)
	{
	  /*	  TimeBin_BH_mass[bin] += P[i].BH_Mass;
	  TimeBin_BH_dynamicalmass[bin] += P[i].Mass;
	  TimeBin_BH_Mdot[bin] += P[i].BH_Mdot;*/
#ifdef DUST
	    TimeBin_BH_mass[bin] += P[i].Dust_Mass;
	    TimeBin_BH_dynamicalmass[bin] += P[i].Total_Mass;
	    TimeBin_BH_Mdot[bin] += P[i].Mass; /* SN: gives total planet mass NOW */
#else
	    TimeBin_BH_mass[bin] += P[i].BH_Mass;
	    TimeBin_BH_dynamicalmass[bin] += P[i].Mass;
	    TimeBin_BH_Mdot[bin] += P[i].BH_Mdot;
#endif
	}
#endif
    }

  FirstActiveParticle = -1;

  for(n = 0, prev = -1; n < TIMEBINS; n++)
    {
      if(TimeBinActive[n])
	for(i = FirstInTimeBin[n]; i >= 0; i = NextInTimeBin[i])
	  {
	    if(prev == -1)
	      FirstActiveParticle = i;

	    if(prev >= 0)
	      NextActiveParticle[prev] = i;

	    prev = i;
	  }
    }

  if(prev >= 0)
    NextActiveParticle[prev] = -1;


  int sum1, sum2;

  MPI_Allreduce(&NumForceUpdate, &sum1, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  for(i = FirstActiveParticle, NumForceUpdate = 0; i >= 0; i = NextActiveParticle[i])
    {
      NumForceUpdate++;
      if(i >= NumPart)
	{
	  printf("Bummer i=%d\n", i);
	  endrun(12);
	}

    }
  MPI_Allreduce(&NumForceUpdate, &sum2, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(ThisTask == 0)
    {
      printf("sum1=%d sum2=%d\n", sum1, sum2);
    }

   if(sum1 != sum2 && All.NumCurrentTiStep > 0)
    endrun(121); 

}



void drift_particle(int i, int time1)
{
  int j, time0, dt_step;
  double dt_drift, dt_gravkick, dt_hydrokick, dt_entr;

#ifdef VIRTUAL
  if (P[i].Type == 3) return;
#endif

#ifdef SOFTEREQS
  double a3inv;

  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;
#endif

  time0 = P[i].Ti_current;

  if(time1 < time0)
    {
      printf("i=%d time0=%d time1=%d\n", i, time0, time1);
      endrun(12);
    }

  if(time1 == time0)
    return;

  if(All.ComovingIntegrationOn)
    {
      dt_drift = get_drift_factor(time0, time1);
      dt_gravkick = get_gravkick_factor(time0, time1);
      dt_hydrokick = get_hydrokick_factor(time0, time1);
    }
  else
    {
      dt_drift = dt_gravkick = dt_hydrokick = (time1 - time0) * All.Timebase_interval;
    }


  for(j = 0; j < 3; j++)
    P[i].Pos[j] += P[i].Vel[j] * dt_drift;

#ifdef DISTORTIONTENSOR
  for(j = 0; j < 9; j++)
    P[i].distortion_tensor[j] += P[i].distortion_tensor_vel[j] * dt_drift;
#endif

#ifndef HPM
  if(P[i].Type == 0)
    {
#ifdef PMGRID
      for(j = 0; j < 3; j++)
	SphP[i].VelPred[j] +=
	  (P[i].g.GravAccel[j] + P[i].GravPM[j]) * dt_gravkick + SphP[i].a.HydroAccel[j] * dt_hydrokick;
#else
      for(j = 0; j < 3; j++)
	SphP[i].VelPred[j] += P[i].g.GravAccel[j] * dt_gravkick + SphP[i].a.HydroAccel[j] * dt_hydrokick;
#endif

#ifdef RAD_ACCEL
      for(j = 0; j < 3; j++)
          SphP[i].VelPred[j] += SphP[i].ra.RadAccel[j] * dt_hydrokick;
#endif

#ifdef DUST
      for(j = 0; j < 3; j++)
          SphP[i].VelPred[j] += SphP[i].da.DragAccel[j] * dt_hydrokick;
#endif

#ifdef GAS_FIXED
      for(j = 0; j < 3; j++)
          SphP[i].VelPred[j] = 0.;
#endif

      SphP[i].d.Density *= exp(-SphP[i].v.DivVel * dt_drift);
      PPP[i].Hsml *= exp(0.333333333333 * SphP[i].v.DivVel * dt_drift);

      if(PPP[i].Hsml < All.MinGasHsml)
	PPP[i].Hsml = All.MinGasHsml;


      dt_step = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0);

      dt_entr = (time1 - (P[i].Ti_begstep + dt_step / 2)) * All.Timebase_interval;

#ifndef MHM
#ifndef SOFTEREQS
      SphP[i].Pressure = (SphP[i].Entropy + SphP[i].e.DtEntropy * dt_entr) * pow(SphP[i].d.Density, GAMMA);
#else
      if(SphP[i].d.Density * a3inv >= All.PhysDensThresh)
	SphP[i].Pressure =
	  All.FactorForSofterEQS * (SphP[i].Entropy +
				    SphP[i].e.DtEntropy * dt_entr) * pow(SphP[i].d.Density,
									 GAMMA) + (1 -
										   All.
										   FactorForSofterEQS) *
	  GAMMA_MINUS1 * SphP[i].d.Density * All.InitGasU;
      else
	SphP[i].Pressure = (SphP[i].Entropy + SphP[i].e.DtEntropy * dt_entr) * pow(SphP[i].d.Density, GAMMA);
#endif
#else
      /* Here we use an isothermal equation of state */
      SphP[i].Pressure = GAMMA_MINUS1 * SphP[i].d.Density * All.InitGasU;
      SphP[i].Entropy = SphP[i].Pressure / pow(SphP[i].d.Density, GAMMA);
#endif

#ifdef COSMIC_RAYS
#if ( defined( CR_UPDATE_PARANOIA ) )
      CR_Particle_Update(SphP + i);
#endif
#ifndef CR_NOPRESSURE
      SphP[i].Pressure += CR_Comoving_Pressure(SphP + i);
#endif
#endif


#ifdef MAGNETIC
      for(j = 0; j < 3; j++)
	SphP[i].BPred[j] += SphP[i].DtB[j] * dt_entr;
#ifdef DIVBCLEANING_DEDNER
      SphP[i].PhiPred += SphP[i].DtPhi * dt_entr;
#endif
#endif

    }
#endif /* end of HPM */

  P[i].Ti_current = time1;
}



void move_particles(int time1)
{
  int i;

  if(ThisTask == 0)
    printf("MOVE\n");

  for(i = 0; i < NumPart; i++)
    drift_particle(i, time1);
}



/*! This function makes sure that all particle coordinates (Pos) are
 *  periodically mapped onto the interval [0, BoxSize].  After this function
 *  has been called, a new domain decomposition should be done, which will
 *  also force a new tree construction.
 */
#ifdef PERIODIC
void do_box_wrapping(void)
{
  int i, j;
  double boxsize[3];

  for(j = 0; j < 3; j++)
    boxsize[j] = All.BoxSize;

#ifdef LONG_X
  boxsize[0] *= LONG_X;
#endif
#ifdef LONG_Y
  boxsize[1] *= LONG_Y;
#endif
#ifdef LONG_Z
  boxsize[2] *= LONG_Z;
#endif

  for(i = 0; i < NumPart; i++)
    for(j = 0; j < 3; j++)
      {
	while(P[i].Pos[j] < 0)
	  P[i].Pos[j] += boxsize[j];

	while(P[i].Pos[j] >= boxsize[j])
	  P[i].Pos[j] -= boxsize[j];
      }
}
#endif



/*

#ifdef XXLINFO
#ifdef MAGNETIC
  double MeanB_part = 0, MeanB_sum;

#ifdef TRACEDIVB
  double MaxDivB_part = 0, MaxDivB_all;
  double dmax1, dmax2;
#endif
#endif
#ifdef TIME_DEP_ART_VISC
  double MeanAlpha_part = 0, MeanAlpha_sum;
#endif
#endif



#ifdef XXLINFO
      if(Flag_FullStep == 1)
        {
#ifdef MAGNETIC
          MeanB_part += sqrt(SphP[i].BPred[0] * SphP[i].BPred[0] +
                             SphP[i].BPred[1] * SphP[i].BPred[1] + SphP[i].BPred[2] * SphP[i].BPred[2]);
#ifdef TRACEDIVB
          MaxDivB_part = DMAX(MaxDivB, fabs(SphP[i].divB));
#endif
#endif
#ifdef TIME_DEP_ART_VISC
          MeanAlpha_part += SphP[i].alpha;
#endif
        }
#endif


#ifdef XXLINFO
  if(Flag_FullStep == 1)
    {
#ifdef MAGNETIC
      MPI_Reduce(&MeanB_part, &MeanB_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if(ThisTask == 0)
	MeanB = MeanB_sum / All.TotN_gas;
#ifdef TRACEDIVB
      MPI_Reduce(&MaxDivB_part, &MaxDivB_all, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      if(ThisTask == 0)
	MaxDivB = MaxDivB_all;
#endif
#endif
#ifdef TIME_DEP_ART_VISC
      MPI_Reduce(&MeanAlpha_part, &MeanAlpha_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if(ThisTask == 0)
	MeanAlpha = MeanAlpha_sum / All.TotN_gas;
#endif
    }
#endif
*/
