#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"
#ifdef COSMIC_RAYS
#include "cosmic_rays.h"
#endif

/*! \file timestep.c 
 *  \brief routines for 'kicking' particles in
 *  momentum space and assigning new timesteps */

static double fac1, fac2, fac3, hubble_a, atime, a3inv;
static double dt_displacement = 0;

#ifdef PMGRID
static double dt_gravkickA, dt_gravkickB;
#endif


/*! This function advances the system in momentum space, i.e. it does apply the 'kick' operation after the
 *  forces have been computed. Additionally, it assigns new timesteps to particles. At start-up, a
 *  half-timestep is carried out, as well as at the end of the simulation. In between, the half-step kick that
 *  ends the previous timestep and the half-step kick for the new timestep are combined into one operation.
 */
void advance_and_find_timesteps(void)
{
  int i, ti_step, ti_step_old, ti_min, tend, tstart, bin, binold, prev, next;
  double aphys;

#ifdef PMGRID
  int j, dt_step;
  double dt_gravkick, dt_hydrokick;
#endif
#ifdef MAKEGLASS
  double disp, dispmax, globmax, dmean, fac, disp2sum, globdisp2sum;
#endif
#if defined(TIME_DEP_ART_VISC) || defined(TIME_DEP_MAGN_DISP)
  double dmin1, dmin2;
#endif
#ifdef CHEMISTRY
  int ifunc, mode;
  double a_start, a_end;
#endif

  CPU_Step[CPU_MISC] += measure_time();

  if(All.ComovingIntegrationOn)
    {
      fac1 = 1 / (All.Time * All.Time);
      fac2 = 1 / pow(All.Time, 3 * GAMMA - 2);
      fac3 = pow(All.Time, 3 * (1 - GAMMA) / 2.0);
      hubble_a = hubble_function(All.Time);
      a3inv = 1 / (All.Time * All.Time * All.Time);
      atime = All.Time;
    }
  else
    fac1 = fac2 = fac3 = hubble_a = a3inv = atime = 1;

  if(Flag_FullStep || dt_displacement == 0)
    find_dt_displacement_constraint(hubble_a * atime * atime);

#ifdef PMGRID
  if(All.ComovingIntegrationOn)
    dt_gravkickB = get_gravkick_factor(All.PM_Ti_begstep, All.Ti_Current) -
      get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);
  else
    dt_gravkickB = (All.Ti_Current - (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2) * All.Timebase_interval;
#endif


#ifdef MAKEGLASS
  for(i = 0, dispmax = 0, disp2sum = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  P[i].g.GravAccel[j] *= -1;
#ifdef PMGRID
	  P[i].GravPM[j] *= -1;
	  P[i].g.GravAccel[j] += P[i].GravPM[j];
	  P[i].GravPM[j] = 0;
#endif
	}

      disp = sqrt(P[i].g.GravAccel[0] * P[i].g.GravAccel[0] +
		  P[i].g.GravAccel[1] * P[i].g.GravAccel[1] + P[i].g.GravAccel[2] * P[i].g.GravAccel[2]);

      disp *= 2.0 / (3 * All.Hubble * All.Hubble);

      disp2sum += disp * disp;

      if(disp > dispmax)
	dispmax = disp;
    }

  MPI_Allreduce(&dispmax, &globmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&disp2sum, &globdisp2sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  dmean = pow(P[0].Mass / (All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)), 1.0 / 3);

  if(globmax > dmean)
    fac = dmean / globmax;
  else
    fac = 1.0;

  if(ThisTask == 0)
    {
      printf("\nglass-making:  dmean= %g  global disp-maximum= %g  rms= %g\n\n",
	     dmean, globmax, sqrt(globdisp2sum / All.TotNumPart));
      fflush(stdout);
    }

  for(i = 0, dispmax = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  P[i].Vel[j] = 0;
	  P[i].Pos[j] += fac * P[i].g.GravAccel[j] * 2.0 / (3 * All.Hubble * All.Hubble);
	  P[i].g.GravAccel[j] = 0;
	}
    }
#endif


  //  All.total_heating = 0.;

  All.DoDynamicUpdate = ShouldWeDoDynamicUpdate();

  /* Now assign new timesteps and kick */

  if(All.DoDynamicUpdate)
    {
      GlobFlag++;
      DomainNumChanged = 0;
      DomainList = (int *) mymalloc(NTopleaves * sizeof(int));
      if(ThisTask == 0)
	printf("kicks will prepare for dynamic update of tree\n");
    }

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
  {
	  ti_step = get_timestep(i, &aphys, 0);
	  
	  /* make it a power 2 subdivision */
	  ti_min = TIMEBASE;
	  while(ti_min > ti_step)
		  ti_min >>= 1;
	  ti_step = ti_min;
	  
	  bin = get_timestep_bin(ti_step);
	  binold = P[i].TimeBin;
	  
	  if(bin > binold)		/* timestep wants to increase */
	  {
		  if(TimeBinActive[bin] == 0)	/* leave at old step if not synchronized */
		  {
			  bin = binold;
			  ti_step = bin ? (1 << bin) : 0;
		  }
	  }
	  
	  if(All.Ti_Current >= TIMEBASE)	/* we here finish the last timestep. */
	  {
		  ti_step = 0;
		  bin = 0;
	  }
	  
	  if((TIMEBASE - All.Ti_Current) < ti_step)	/* check that we don't run beyond the end */
	  {
		  endrun(888);		/* should not happen */
		  ti_step = TIMEBASE - All.Ti_Current;
		  ti_min = TIMEBASE;
		  while(ti_min > ti_step)
			  ti_min >>= 1;
		  ti_step = ti_min;
	  }
	  
	  if(bin != binold)
	  {
		  TimeBinCount[binold]--;
		  if(P[i].Type == 0)
		  {
			  TimeBinCountSph[binold]--;
#ifdef SFR
			  TimeBinSfr[binold] -= SphP[i].Sfr;
			  TimeBinSfr[bin] += SphP[i].Sfr;
#endif
		  }

#ifdef BLACK_HOLES
		  if(P[i].Type == 5)
		  {
#ifndef DUST
			  TimeBin_BH_mass[binold] -= P[i].BH_Mass;
			  TimeBin_BH_dynamicalmass[binold] -= P[i].Mass;
			  TimeBin_BH_Mdot[binold] -= P[i].BH_Mdot;
			  TimeBin_BH_mass[bin] += P[i].BH_Mass;
			  TimeBin_BH_dynamicalmass[bin] += P[i].Mass;
			  TimeBin_BH_Mdot[bin] += P[i].BH_Mdot;
#else
			  TimeBin_BH_mass[binold] -= P[i].Dust_Mass;
			  TimeBin_BH_dynamicalmass[binold] -= P[i].Total_Mass;
			  TimeBin_BH_Mdot[binold] -= P[i].Mass;
			  TimeBin_BH_mass[bin] += P[i].Dust_Mass;
			  TimeBin_BH_dynamicalmass[bin] += P[i].Total_Mass;
			  TimeBin_BH_Mdot[bin] += P[i].Mass;
#endif
		  }
#endif

		  prev = PrevInTimeBin[i];
		  next = NextInTimeBin[i];
		  
		  if(FirstInTimeBin[binold] == i)
			  FirstInTimeBin[binold] = next;
		  if(LastInTimeBin[binold] == i)
			  LastInTimeBin[binold] = prev;
		  if(prev >= 0)
			  NextInTimeBin[prev] = next;
		  if(next >= 0)
			  PrevInTimeBin[next] = prev;
		  
		  if(TimeBinCount[bin] > 0)
		  {
			  PrevInTimeBin[i] = LastInTimeBin[bin];
			  NextInTimeBin[LastInTimeBin[bin]] = i;
			  NextInTimeBin[i] = -1;
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
		  
		  P[i].TimeBin = bin;
	  }

		
	  ti_step_old = binold ? (1 << binold) : 0;
	  
	  tstart = P[i].Ti_begstep + ti_step_old / 2;	/* midpoint of old step */
	  tend = P[i].Ti_begstep + ti_step_old + ti_step / 2;	/* midpoint of new step */
	  
	  P[i].Ti_begstep += ti_step_old;
	  
#ifdef DUST
	  //	  All.total_heating += SphP[i].da.DragAccel[0]*P[i].Mass;
#endif

	  do_the_kick(i, tstart, tend, P[i].Ti_begstep);
  }
  if(All.DoDynamicUpdate)
  {
	  force_finish_kick_nodes();
	  myfree(DomainList);
  }
  
  /*  printf(" Total_heating ====== %g \n", All.total_heating/(All.Time + 1.e-35));
  fflush(stdout);
  */

#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current)	/* need to do long-range kick */
    {
      All.DoDynamicUpdate = 0;

      ti_step = TIMEBASE;
      while(ti_step > (dt_displacement / All.Timebase_interval))
	ti_step >>= 1;

      if(ti_step > (All.PM_Ti_endstep - All.PM_Ti_begstep))	/* PM-timestep wants to increase */
	{
	  /* we only increase if an integer number of steps will bring us to the end */
	  if(((TIMEBASE - All.PM_Ti_endstep) % ti_step) > 0)
	    ti_step = All.PM_Ti_endstep - All.PM_Ti_begstep;	/* leave at old step */
	}

      if(All.Ti_Current == TIMEBASE)	/* we here finish the last timestep. */
	ti_step = 0;

      tstart = (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2;
      tend = All.PM_Ti_endstep + ti_step / 2;

      if(All.ComovingIntegrationOn)
	dt_gravkick = get_gravkick_factor(tstart, tend);
      else
	dt_gravkick = (tend - tstart) * All.Timebase_interval;

      All.PM_Ti_begstep = All.PM_Ti_endstep;
      All.PM_Ti_endstep = All.PM_Ti_begstep + ti_step;

      if(All.ComovingIntegrationOn)
	dt_gravkickB = -get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);
      else
	dt_gravkickB =
	  -((All.PM_Ti_begstep + All.PM_Ti_endstep) / 2 - All.PM_Ti_begstep) * All.Timebase_interval;

      for(i = 0; i < NumPart; i++)
	{
	  for(j = 0; j < 3; j++)	/* do the kick */
#ifdef BH_FIXED
      if (P[i].Mass > 0.45*All.SMBHmass) 
      {
	  P[i].g.GravAccel[0]=P[i].g.GravAccel[1]=P[i].g.GravAccel[2]=0.;
	  P[i].Vel[j] = 0.;
      }
#endif
	    P[i].Vel[j] += P[i].GravPM[j] * dt_gravkick;

	  if(P[i].Type == 0)
	    {
	      dt_step = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0);

	      if(All.ComovingIntegrationOn)
		{
		  dt_gravkickA = get_gravkick_factor(P[i].Ti_begstep, All.Ti_Current) -
		    get_gravkick_factor(P[i].Ti_begstep, P[i].Ti_begstep + dt_step / 2);
		  dt_hydrokick = get_hydrokick_factor(P[i].Ti_begstep, All.Ti_Current) -
		    get_hydrokick_factor(P[i].Ti_begstep, P[i].Ti_begstep + dt_step / 2);
		}
	      else
		dt_gravkickA = dt_hydrokick =
		  (All.Ti_Current - (P[i].Ti_begstep + dt_step / 2)) * All.Timebase_interval;

	      for(j = 0; j < 3; j++)
		SphP[i].VelPred[j] = P[i].Vel[j]
		  + P[i].g.GravAccel[j] * dt_gravkickA
		  + SphP[i].a.HydroAccel[j] * dt_hydrokick + P[i].GravPM[j] * dt_gravkickB;

#ifdef RAD_ACCEL
	      for(j = 0; j < 3; j++)
		SphP[i].VelPred[j] = P[i].Vel[j] + SphP[i].ra.RadAccel[j] * dt_hydrokick;
#endif

#ifdef DUST
	      /*	      for(j = 0; j < 3; j++)
			      SphP[i].VelPred[j] = P[i].Vel[j] + SphP[i].da.DragAccel[j] * dt_hydrokick;*/
#endif

#ifdef GAS_FIXED
	      for(j = 0; j < 3; j++) {
	        P[i].Vel[j] = 0.;
		SphP[i].VelPred[j] = 0.0;
              }
#endif

	    }
	}
    }
#endif

  CPU_Step[CPU_TIMELINE] += measure_time();
}



void do_the_kick(int i, int tstart, int tend, int tcurrent)
{
  int j;
  MyFloat dv[3];
  double minentropy;
  double dt_entr, dt_gravkick, dt_hydrokick, dt_gravkick2, dt_hydrokick2, dt_entr2;

#ifdef DISTORTIONTENSOR
  MyFloat dv_distortion_tensor[9];
#endif

#ifdef VIRTUAL
  if (P[i].Type == 3) return;
#endif

  if(All.ComovingIntegrationOn)
    {
      dt_entr = (tend - tstart) * All.Timebase_interval;
      dt_entr2 = (tend - tcurrent) * All.Timebase_interval;
      dt_gravkick = get_gravkick_factor(tstart, tend);
      dt_hydrokick = get_hydrokick_factor(tstart, tend);
      dt_gravkick2 = get_gravkick_factor(tcurrent, tend);
      dt_hydrokick2 = get_hydrokick_factor(tcurrent, tend);
    }
  else
    {
      dt_entr = dt_gravkick = dt_hydrokick = (tend - tstart) * All.Timebase_interval;
      dt_gravkick2 = dt_hydrokick2 = dt_entr2 = (tend - tcurrent) * All.Timebase_interval;
    }


#ifdef DUST
  if (P[i].Type == 2) 
    {
      double dust_acc = 0.;
      for(j = 0; j < 3; j++){
	dust_acc += P[i].DeltaDustMomentum[j]*P[i].DeltaDustMomentum[j];
      }
	dust_acc = sqrt(dust_acc)/(P[i].Mass)/(dt_gravkick + dt_gravkick2);
    }
#endif /* DUST */

  /* do the kick */

  for(j = 0; j < 3; j++)
    {
      dv[j] = P[i].g.GravAccel[j] * dt_gravkick;

#ifdef DUST
      if (P[i].Type == 2) 
	{
	  dv[j] = 0.0;
	}
#endif /* DUST */

#ifdef BH_FIXED
      if (P[i].Mass > 0.45 * All.SMBHmass) dv[j] = 0.;
#endif
      P[i].Vel[j] += dv[j];
#ifdef GAS_FIXED
      if(P[i].Type == 0) P[i].Vel[j] = 0.;
#endif
    }

#ifdef DISTORTIONTENSOR
  /*my set of equations */
  dv_distortion_tensor[0] =
    (P[i].tite.tidal_tensor[0] * P[i].distortion_tensor[0] +
     P[i].tite.tidal_tensor[1] * P[i].distortion_tensor[3] +
     P[i].tite.tidal_tensor[2] * P[i].distortion_tensor[6]) * dt_gravkick;
  P[i].distortion_tensor_vel[0] += dv_distortion_tensor[0];

  dv_distortion_tensor[1] =
    (P[i].tite.tidal_tensor[0] * P[i].distortion_tensor[1] +
     P[i].tite.tidal_tensor[1] * P[i].distortion_tensor[4] +
     P[i].tite.tidal_tensor[2] * P[i].distortion_tensor[7]) * dt_gravkick;
  P[i].distortion_tensor_vel[1] += dv_distortion_tensor[1];

  dv_distortion_tensor[2] =
    (P[i].tite.tidal_tensor[0] * P[i].distortion_tensor[2] +
     P[i].tite.tidal_tensor[1] * P[i].distortion_tensor[5] +
     P[i].tite.tidal_tensor[2] * P[i].distortion_tensor[8]) * dt_gravkick;
  P[i].distortion_tensor_vel[2] += dv_distortion_tensor[2];

  dv_distortion_tensor[3] =
    (P[i].tite.tidal_tensor[1] * P[i].distortion_tensor[0] +
     P[i].tite.tidal_tensor[3] * P[i].distortion_tensor[3] +
     P[i].tite.tidal_tensor[4] * P[i].distortion_tensor[6]) * dt_gravkick;
  P[i].distortion_tensor_vel[3] += dv_distortion_tensor[3];

  dv_distortion_tensor[4] =
    (P[i].tite.tidal_tensor[1] * P[i].distortion_tensor[1] +
     P[i].tite.tidal_tensor[3] * P[i].distortion_tensor[4] +
     P[i].tite.tidal_tensor[4] * P[i].distortion_tensor[7]) * dt_gravkick;
  P[i].distortion_tensor_vel[4] += dv_distortion_tensor[4];

  dv_distortion_tensor[5] =
    (P[i].tite.tidal_tensor[1] * P[i].distortion_tensor[2] +
     P[i].tite.tidal_tensor[3] * P[i].distortion_tensor[5] +
     P[i].tite.tidal_tensor[4] * P[i].distortion_tensor[8]) * dt_gravkick;
  P[i].distortion_tensor_vel[5] += dv_distortion_tensor[5];

  dv_distortion_tensor[6] =
    (P[i].tite.tidal_tensor[2] * P[i].distortion_tensor[0] +
     P[i].tite.tidal_tensor[4] * P[i].distortion_tensor[3] +
     P[i].tite.tidal_tensor[5] * P[i].distortion_tensor[6]) * dt_gravkick;
  P[i].distortion_tensor_vel[6] += dv_distortion_tensor[6];

  dv_distortion_tensor[7] =
    (P[i].tite.tidal_tensor[2] * P[i].distortion_tensor[1] +
     P[i].tite.tidal_tensor[4] * P[i].distortion_tensor[4] +
     P[i].tite.tidal_tensor[5] * P[i].distortion_tensor[7]) * dt_gravkick;
  P[i].distortion_tensor_vel[7] += dv_distortion_tensor[7];

  dv_distortion_tensor[8] =
    (P[i].tite.tidal_tensor[2] * P[i].distortion_tensor[2] +
     P[i].tite.tidal_tensor[4] * P[i].distortion_tensor[5] +
     P[i].tite.tidal_tensor[5] * P[i].distortion_tensor[8]) * dt_gravkick;
  P[i].distortion_tensor_vel[8] += dv_distortion_tensor[8];
#endif

  if(P[i].Type == 0)		/* SPH stuff */
    {
      //      All.total_heating += SphP[i].da.DragAccel[0]*(dt_hydrokick+dt_hydrokick2);

      for(j = 0; j < 3; j++)
	{
	  dv[j] += SphP[i].a.HydroAccel[j] * dt_hydrokick;
	  P[i].Vel[j] += SphP[i].a.HydroAccel[j] * dt_hydrokick;

	  SphP[i].VelPred[j] =
	    P[i].Vel[j] - dt_gravkick2 * P[i].g.GravAccel[j] - dt_hydrokick2 * SphP[i].a.HydroAccel[j];

#ifdef RAD_ACCEL
	  dv[j] += SphP[i].ra.RadAccel[j] * dt_hydrokick;
	  P[i].Vel[j] += SphP[i].ra.RadAccel[j] * dt_hydrokick;
	  SphP[i].VelPred[j] = P[i].Vel[j] - dt_hydrokick2 * SphP[i].ra.RadAccel[j];
#endif

#ifdef DUST
	  /* dv[j] += SphP[i].da.DragAccel[j] * dt_hydrokick; */
	  /* P[i].Vel[j] += SphP[i].da.DragAccel[j] * dt_hydrokick; */
	  /* SphP[i].VelPred[j] = P[i].Vel[j] - dt_hydrokick2 * SphP[i].da.DragAccel[j]; */
	  /* SN: reset this to zero after the force has been applied */
	  SphP[i].da.DragAccel[j] = 0.;
#endif

#ifdef PMGRID
	  SphP[i].VelPred[j] += P[i].GravPM[j] * dt_gravkickB;
#endif

#ifdef GAS_FIXED
          P[i].Vel[j] = 0.;
          SphP[i].VelPred[j] = 0.;
#endif

#ifdef MAGNETIC
	  SphP[i].B[j] += SphP[i].DtB[j] * dt_entr;
	  SphP[i].BPred[j] = SphP[i].B[j] - SphP[i].DtB[j] * dt_entr2;
#endif
	}
#if defined(MAGNETIC) && defined(DIVBCLEANING_DEDNER)
      SphP[i].Phi += SphP[i].DtPhi * dt_entr;
      SphP[i].PhiPred = SphP[i].Phi - SphP[i].DtPhi * dt_entr2;
#endif
#ifdef TIME_DEP_ART_VISC
      SphP[i].alpha += SphP[i].Dtalpha * dt_entr;
      SphP[i].alpha = DMIN(SphP[i].alpha, All.ArtBulkViscConst);
      if(SphP[i].alpha < All.AlphaMin)
	SphP[i].alpha = All.AlphaMin;
#endif
#ifdef TIME_DEP_MAGN_DISP
      SphP[i].Balpha += SphP[i].DtBalpha * dt_entr;
      SphP[i].Balpha = DMIN(SphP[i].Balpha, All.ArtMagDispConst);
      if(SphP[i].Balpha < All.ArtMagDispMin)
	SphP[i].Balpha = All.ArtMagDispMin;
#endif
      /* In case of cooling, we prevent that the entropy (and
         hence temperature decreases by more than a factor 0.5 */

#ifndef DO_NOT_PROTECT

      if(SphP[i].e.DtEntropy * dt_entr > -0.5 * SphP[i].Entropy)
	SphP[i].Entropy += SphP[i].e.DtEntropy * dt_entr;
      else
	SphP[i].Entropy *= 0.5;
#else 
      if(SphP[i].e.DtEntropy * dt_entr > -0.9 * SphP[i].Entropy)
	SphP[i].Entropy += SphP[i].e.DtEntropy * dt_entr;
      else
	SphP[i].Entropy *= 0.1;
#endif



#ifdef CHEMISTRY
      /* update the chemical abundances for the new density and temperature */
      a_start = All.TimeBegin * exp(tstart * All.Timebase_interval);
      a_end = All.TimeBegin * exp(tend * All.Timebase_interval);

      /* time in cosmic expansion parameter */
      ifunc = compute_abundances(mode = 1, i, a_start, a_end);
#endif

      if(All.MinEgySpec)
	{
	  minentropy = All.MinEgySpec * GAMMA_MINUS1 / pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);
	  if(SphP[i].Entropy < minentropy)
	    {
	      SphP[i].Entropy = minentropy;
	      SphP[i].e.DtEntropy = 0;
	    }
	}

      /* In case the timestep increases in the new step, we
         make sure that we do not 'overcool'. */

      dt_entr = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) / 2 * All.Timebase_interval;
      if(SphP[i].Entropy + SphP[i].e.DtEntropy * dt_entr < 0.5 * SphP[i].Entropy)
	SphP[i].e.DtEntropy = -0.5 * SphP[i].Entropy / dt_entr;
    }

//#ifndef EXCLUDE_VIRTUAL_FROM_TREE
	 if(All.DoDynamicUpdate)
		 force_kick_node(i, dv);
/* #else
	 if (P[i].Type != 3)
		 if(All.DoDynamicUpdate)
			 force_kick_node(i, dv);
#endif
*/
  
}



/*! This function normally (for flag==0) returns the maximum allowed timestep of a particle, expressed in
 *  terms of the integer mapping that is used to represent the total simulated timespan. The physical
 *  acceleration is returned in aphys. The latter is used in conjunction with the PSEUDOSYMMETRIC integration
 *  option, which also makes of the second function of get_timestep. When it is called with a finite timestep
 *  for flag, it returns the physical acceleration that would lead to this timestep, assuming timestep
 *  criterion 0.
 */
int get_timestep(int p,		/*!< particle index */
		 double *aphys,	/*!< acceleration (physical units) */
		 int flag	/*!< either 0 for normal operation, or finite timestep to get corresponding
				   aphys */ )
{
  double ax, ay, az, ac;
  double csnd = 0, dt = 0, dt_courant = 0, rho;
  int ti_step;

#ifdef BLACK_HOLES
  double dt_accr, dt_a = 1.e10;
#endif

#ifdef VIRTUAL
  double dt_virt;
  double dt_thin, dt_thick, dt_abs;

#endif

#ifdef CONDUCTION
  double dt_cond;
#endif

#ifdef NS_TIMESTEP
  double dt_NS = 0;
#endif

#ifdef CR_DIFFUSION
  double dt_CR;
#endif

#ifdef NONEQUILIBRIUM
  double dt_cool, dt_elec;
#endif

  if(flag == 0)
    {
      ax = fac1 * P[p].g.GravAccel[0];
      ay = fac1 * P[p].g.GravAccel[1];
      az = fac1 * P[p].g.GravAccel[2];

#ifdef PMGRID
      ax += fac1 * P[p].GravPM[0];
      ay += fac1 * P[p].GravPM[1];
      az += fac1 * P[p].GravPM[2];
#endif

#ifdef DUST
      /* if(P[p].Type == 2) */
      /* 	{ */
      /* 	  ax += P[p].DeltaDustMomentum[0]/P[p].Mass/(dt_gravkick + dt_gravkick2); */
      /* 	  ay += P[p].DeltaDustMomentum[1]/P[p].Mass/(dt_gravkick + dt_gravkick2); */
      /* 	  az += P[p].DeltaDustMomentum[2]/P[p].Mass/(dt_gravkick + dt_gravkick2); */
      /* 	} */
#endif

      if(P[p].Type == 0)
	{
	  ax += fac2 * SphP[p].a.HydroAccel[0];
	  ay += fac2 * SphP[p].a.HydroAccel[1];
	  az += fac2 * SphP[p].a.HydroAccel[2];
#ifdef RAD_ACCEL
	  ax += fac2 * SphP[p].ra.RadAccel[0];
	  ay += fac2 * SphP[p].ra.RadAccel[1];
	  az += fac2 * SphP[p].ra.RadAccel[2];
#endif
#ifdef DUST
	  ax += fac2 * SphP[p].da.DragAccel[0];
	  ay += fac2 * SphP[p].da.DragAccel[1];
	  az += fac2 * SphP[p].da.DragAccel[2];
#endif
	}

      ac = sqrt(ax * ax + ay * ay + az * az);	/* this is now the physical acceleration */
#ifdef GAS_FIXED
      ax = 0.;
      ay = 0.;
      az = 0.;
      ac = 0.;
#endif
      *aphys = ac;
    }
  else
    ac = *aphys;

  if(ac == 0)
    ac = 1.0e-30;


  switch (All.TypeOfTimestepCriterion)
    {
    case 0:
      if(flag > 0)
	{
	  dt = flag * All.Timebase_interval;

	  dt /= hubble_a;	/* convert dloga to physical timestep  */

	  ac = 2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / (dt * dt);
	  *aphys = ac;
	  return flag;
	}
      dt = sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / ac);
#ifdef DUST_TIMESTEP
      if(P[p].Type == 0 && P[p].Mass > 0)
      	{
      	  ax += P[p].DeltaDustMomentum[0]/P[p].Mass/dt;
      	  ay += P[p].DeltaDustMomentum[1]/P[p].Mass/dt;
      	  az += P[p].DeltaDustMomentum[2]/P[p].Mass/dt;
      	  ac = sqrt(ax * ax + ay * ay + az * az);	/* this is now the physical acceleration */
      	  double dt_new;
	  if (ac > 0) dt_new = sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / ac);
	  /*      	  printf(" dt_new %g, dt_old  %g, acc =   %g \n",dt_new , dt, ac);
			  fflush(stdout);*/
	  	  if (dt_new < dt) dt = dt_new;
      	}
#endif

#ifdef DUST
      if(P[p].Type == 2) dt = dt/2; 
      // if(P[p].Type == 2) dt = dt/16; 
//SN: dust is integrated differently, in dust.c, so larger time steps should be possible
      /* if(P[p].Type == 2 && P[p].Mass > 0.) { */
      /* 	double dv_dust = sqrt(P[p].DeltaDustMomentum[0]*P[p].DeltaDustMomentum[0] + P[p].DeltaDustMomentum[1]*P[p].DeltaDustMomentum[1] +  */
      /* 			      P[p].DeltaDustMomentum[2]*P[p].DeltaDustMomentum[2])/P[p].Mass; */
      /* 	//	double dacc_dust = dv_dust/dt; */
      /* 	double v_dust =  sqrt(P[p].Vel[0]*P[p].Vel[0]+ P[p].Vel[1]*P[p].Vel[1] + P[p].Vel[2]*P[p].Vel[2]); */
      /* 	double dt_dust = All.ErrTolIntAccuracy*v_dust/sqrt(dv_dust*dv_dust + 1.d-20); */
      /* 	if (dt >= dt_dust) dt = dt_dust; */
      /* } */
#endif


#ifdef ADAPTIVE_GRAVSOFT_FORGAS
#ifdef ADAPTIVE_GRAVSOFT_FORGAS_HSML
      if(P[p].Type == 0)
	dt = sqrt(2 * All.ErrTolIntAccuracy * atime * PPP[p].Hsml / 2.8 / ac);
#else
      if(P[p].Type == 0)
	dt =
	  sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] *
	       pow(P[p].Mass / All.ReferenceGasMass, 1.0 / 3) / ac);
#endif
#endif
      break;

    default:
      endrun(888);
      break;
    }


  if(P[p].Type == 0)
    {
      csnd = sqrt(GAMMA * SphP[p].Pressure / SphP[p].d.Density);

#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
      double dmax1, dmax2;

      if(All.ComovingIntegrationOn)
	dt_courant = All.CourantFac * All.Time * DMAX(PPP[p].Hsml, All.SofteningTable[0]) / (fac3 * csnd);
      else
	dt_courant = All.CourantFac * DMAX(PPP[p].Hsml, All.SofteningTable[0]) / csnd;

      if(dt_courant > 2 * All.CourantFac * SphP[p].MinViscousDt)
	dt_courant = 2 * All.CourantFac * SphP[p].MinViscousDt;
#else
      if(All.ComovingIntegrationOn)
	dt_courant = 2 * All.CourantFac * All.Time * PPP[p].Hsml / (fac3 * SphP[p].MaxSignalVel);
      else
	dt_courant = 2 * All.CourantFac * PPP[p].Hsml / SphP[p].MaxSignalVel;
#endif

      if(dt_courant < dt)
	dt = dt_courant;

#ifdef MYFALSE
      dt_viscous = All.CourantFac * SphP[p].MaxViscStep / hubble_a;	/* to convert dloga to physical dt */

      if(dt_viscous < dt)
	dt = dt_viscous;
#endif

#ifdef NS_TIMESTEP
      if(fabs(SphP[p].ViscEntropyChange))
	{
	  dt_NS = VISC_TIMESTEP_PARAMETER * SphP[p].Entropy / SphP[p].ViscEntropyChange / hubble_a;

	  if(dt_NS < dt)
	    dt = dt_NS;
	}
#endif

#ifdef CONDUCTION
      if(fabs(SphP[p].CondEnergyChange))
	{
	  dt_cond =
	    COND_TIMESTEP_PARAMETER * SphP[p].Entropy / GAMMA_MINUS1 * pow(SphP[p].d.Density * a3inv,
									   GAMMA_MINUS1) /
	    fabs(SphP[p].CondEnergyChange);

	  /* convert from dloga back to dt */
	  dt_cond /= hubble_a;

	  if(dt_cond < dt)
	    dt = dt_cond;
	}
#endif

#if (defined(COSMIC_RAYS) && defined(COOLING) && defined(CR_TIMESTEP))
      if(SphP[p].CR_C0 != 0.0)
	{
	  CR_Particle_GetThermalizationRate(SphP + p, &dE_CR, &dN_CR);

	  dt_CR = 0.1 * (SphP[p].CR_E0 / dE_CR);


	  if((dt_CR < dt) && (dt_CR > 0.0))
	    {
	      dt = dt_CR;
	    }

	  dt_CR = 0.1 * (SphP[p].CR_n0 / dN_CR);

	  if((dt_CR < dt) && (dt_CR > 0.0))
	    {
	      dt = dt_CR;
	    }
	}
#endif


#ifdef CR_DIFFUSION
      /* For the diffusion timestep criterion, we look at the green's
       * function for diffusive effects
       * We find that in the exp, we have something like
       * r^2 / 4 alpha t
       * when alpha is the diffusivity.
       * Things should be safe if this factor is smaller than one.
       *
       * The approx. value of r is of the order of (m/rho)^1/3
       *
       */
      if(All.CR_DiffusionCoeff > 0.0)
	{
	  dt_CR =
	    All.CR_DiffusionTimeScale * pow(P[p].Mass / (SphP[p].d.Density * a3inv),
					    2.0 / 3.0) / All.CR_DiffusionCoeff;

	  if(All.CR_DiffusionDensScaling != 0.0)
	    {
	      dt_CR *=
		pow(SphP[p].d.Density * a3inv / All.CR_DiffusionDensZero, -All.CR_DiffusionDensScaling);
	    }

	  if(All.CR_DiffusionEntropyScaling != 0.0)
	    {
	      dt_CR *= pow(SphP[p].Entropy / All.CR_DiffusionEntropyZero, -All.CR_DiffusionEntropyScaling);
	    }


	  if(dt_CR < dt)
	    {
	      dt = dt_CR;

#ifdef CR_DIFFUSION_VERBOSE
	      printf("CR_DIFF timestep had to be applied for particle %i.\n"
		     "dt = %g -- rho_phys = %g -- entropy = %g -- a = %g\n",
		     dt_CR, SphP[p].d.Density * a3inv, SphP[p].Entropy, atime);
#endif
	    }

	}
#endif




    }

#ifdef VIRTUAL

  if (P[p].Type == 3)
   {

       dt_abs =  0.03 * All.OuterBoundary/C * All.UnitVelocity_in_cm_per_s * All.FeedBackVelocity;
       double dt_ff = 1.;
       if (P[p].NewDensity > 0.) dt_ff = 0.002/pow(rho, 0.5);
       if (dt_abs > dt_ff) dt_abs = dt_ff;
       if(dt > dt_abs) dt = dt_abs;

   }




//#endif


#endif

#ifdef BLACK_HOLES
  if(P[p].Type == 5)
    {
	dt_accr = 0.03 * (All.OuterBoundary/100.);
//	printf("Black hole mass = %g, ID = %d \n", P[p].Mass, P[p].ID);
//	fflush(stdout); //SN 12.12.2008: ask for a small timestep for all
//	black holes, not only those that are accreting, because one never
//	knows if accretion will start in the next timestep or a merger will
//	happen quickly
 // if(P[p].BH_Mdot > 0 && P[p].BH_Mass > 0)
	{
#ifdef ACCRETION_RADIUS
            if (P[p].Mass >= 0.45*All.SMBHmass && All.InnerBoundary > 0) 
		dt_a = All.AccDtBlackHole *pow((All.InnerBoundary + 0.5*P[p].Hsml),1.5)/pow(P[p].Mass, 0.5); 
            if (P[p].Mass < 0.45*All.SMBHmass && All.SinkBoundary > 0) 
		dt_a = All.AccDtBlackHole *pow((All.SinkBoundary + 0.5*P[p].Hsml),1.5)/pow(P[p].Mass, 0.5); 
	    if (dt_accr > dt_a) dt_accr = dt_a;
/*	    double dt_cross = 0.5 * All.InnerBoundary/ All.FeedBackVelocity;
	    dt_cross = dt_cross * dt_accr/(dt_cross + dt_accr);
	    double dt_cross = 0.002;
*/


//            if (P[p].Mass > 0.1*All.SMBHmass && dt_accr > 0.02) dt_accr = 0.02;

// Hack by CBP : assign a timestep to the black hole. 
//
// (a) Ensure that the black hole timestep is a small fraction of the timescale over which it accretes 
//     its mass, i.e.
//               dt_accr = All.AccDtBlackHole * P[p].BH_Mass/P[p].BH_Mdot;
//  (31st October 2008)
//	    dt_accr = All.AccDtBlackHole*P[p].BH_Mass/P[p].BH_Mdot;
#endif
	    //SN: this speeds up removal of particles too close to the
            //BH if they accumulate there too quickly
/*	  if(dt_cross < dt)
	    dt = dt_cross;
*/
	  if(dt_accr < dt)
	    dt = dt_accr;
	}

//#ifdef BH_FIXED
//	  dt = dt_accr;//SN: if BH is already fixed, the accretion time step is the only criterium
//#endif

#ifdef NO_BH_ACCRETION 
      dt = dt_accr;
#endif
    }
#endif

#ifdef NONEQUILIBRIUM
  /* another criterion given by the local cooling time */

  if(P[p].Type == 0)
    {
      dt_cool = fabs(SphP[p].t_cool);	/* still in yrs */
      dt_cool *= SEC_PER_YEAR;	/* in seconds */
      dt_cool /= All.UnitTime_in_s;
      dt_cool *= All.HubbleParam;	/* internal units */

      dt_cool = All.Epsilon * dt_cool;


      if(dt_cool > 0 && dt_cool < dt)
	dt = dt_cool;


      /* yet another criterion given by the electron number density change */

      dt_elec = fabs(SphP[p].t_elec);	/* still in yrs */
      dt_elec *= SEC_PER_YEAR;	/* in seconds */
      dt_elec /= All.UnitTime_in_s;
      dt_elec *= All.HubbleParam;	/* internal units */

      dt_elec = All.Epsilon * dt_elec;

      if(dt_elec > 0 && dt_elec < dt)
	dt = dt_elec;
    }
#endif



  /* convert the physical timestep to dloga if needed. Note: If comoving integration has not been selected,
     hubble_a=1.
   */
  dt *= hubble_a;

#ifdef ONLY_PM
  dt = All.MaxSizeTimestep;
#endif



  if(dt >= All.MaxSizeTimestep)
    dt = All.MaxSizeTimestep;


  if(dt >= dt_displacement)
    dt = dt_displacement;

  if(dt < All.MinSizeTimestep)
    {
#ifndef NOSTOP_WHEN_BELOW_MINTIMESTEP
      printf("warning: Timestep wants to be below the limit `MinSizeTimestep'\n");

      if(P[p].Type == 0)
	{
	  printf
	    ("Part-ID=%d  dt=%g dtc=%g ac=%g xyz=(%g|%g|%g)  hsml=%g  maxcsnd=%g dt0=%g eps=%g\n",
	     (int) P[p].ID, dt, dt_courant * hubble_a, ac, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], PPP[p].Hsml,
	     csnd,
	     sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / ac) * hubble_a,
	     All.SofteningTable[P[p].Type]);
#ifdef CONDUCTION
	  printf
	    ("Part-ID=%d  dt_cond=%g  A=%g  rho=%g  condench=%g  dtold=%g\n",
	     (int) P[p].ID, dt_cond * hubble_a, SphP[p].Entropy, SphP[p].d.Density,
	     SphP[p].CondEnergyChange, (P[p].TimeBin ? (1 << P[p].TimeBin) : 0) * All.Timebase_interval);
#endif

#ifdef NS_TIMESTEP
	  printf
	    ("Part-ID=%d  dt_NS=%g  A=%g  rho=%g  dotAvisc=%g  dtold=%g, meanpath=%g \n",
	     (int) P[p].ID, dt_NS * hubble_a, SphP[p].Entropy, SphP[p].d.Density,
	     SphP[p].ViscEntropyChange, (P[p].TimeBin ? (1 << P[p].TimeBin) : 0) * All.Timebase_interval,
	     All.IonMeanFreePath *
	     pow((SphP[p].Entropy * pow(SphP[p].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1),
		 2.0) / SphP[p].d.Density);

	  printf("Stressd=(%g|%g|%g) \n", SphP[p].u.s.StressDiag[0], SphP[p].u.s.StressDiag[1],
		 SphP[p].u.s.StressDiag[2]);
	  printf("Stressoffd=(%g|%g|%g) \n", SphP[p].u.s.StressOffDiag[0], SphP[p].u.s.StressOffDiag[1],
		 SphP[p].u.s.StressOffDiag[2]);
#endif


	}
      else
	{
	  printf("Part-ID=%d  dt=%g ac=%g xyz=(%g|%g|%g)\n", (int) P[p].ID, dt, ac, P[p].Pos[0], P[p].Pos[1],
		 P[p].Pos[2]);
	}
      fflush(stdout);
      endrun(888);
#endif
      dt = All.MinSizeTimestep;
    }

  ti_step = (int) (dt / All.Timebase_interval);

#ifdef CHEMISTRY
  if(ti_step == 0)
    {
      printf("\nError: A timestep of size zero was assigned on the integer timeline!\n"
	     "We better stop.\n"
	     "Task=%d Part-ID=%d dt=%g dt_elec=%g dt_cool=%g tibase=%g ti_step=%d ac=%g xyz=(%g|%g|%g)\n\n",
	     ThisTask, P[p].ID, dt, SphP[p].t_elec, SphP[p].t_cool, All.Timebase_interval, ti_step, ac,
	     P[p].Pos[0], P[p].Pos[1], P[p].Pos[2]);
      fflush(stdout);
      endrun(818);
    }
#endif


  if(!(ti_step > 0 && ti_step < TIMEBASE))
    {
      printf("\nError: A timestep of size zero was assigned on the integer timeline!\n"
	     "We better stop.\n"
	     "Task=%d Part-ID=%d dt=%g dtc=%g dtv=%g dtdis=%g tibase=%g ti_step=%d ac=%g xyz=(%g|%g|%g) tree=(%g|%g%g)\n\n",
	     ThisTask, (int) P[p].ID, dt, dt_courant, dt, dt_displacement,
//	     All.Timebase_interval, ti_step, ac,
	     All.Timebase_interval, ti_step, ac,
	     P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], P[p].g.GravAccel[0], P[p].g.GravAccel[1],
	     P[p].g.GravAccel[2]);
#ifdef PMGRID
      printf("pm_force=(%g|%g|%g)\n", P[p].GravPM[0], P[p].GravPM[1], P[p].GravPM[2]);
#endif
      if(P[p].Type == 0)
      {
	  printf("hydro-frc=(%g|%g|%g) dens=%g hsml=%g\n", SphP[p].a.HydroAccel[0], SphP[p].a.HydroAccel[1],
		 SphP[p].a.HydroAccel[2], SphP[p].d.Density, PPP[p].Hsml);
#ifdef RAD_ACCEL
	  printf("  RadAccel = (x,y,z) = (%g, %g, %g) \n", SphP[p].ra.RadAccel[0], SphP[p].ra.RadAccel[1], SphP[p].ra.RadAccel[2]);
#endif
#ifdef DUST
	  printf("  DragAccel = (x,y,z) = (%g, %g, %g) \n", SphP[p].da.DragAccel[0], SphP[p].da.DragAccel[1], SphP[p].da.DragAccel[2]);
#endif
	  printf(" Injected_BH_Momentum = (x,y,z) = (%g, %g, %g) \n",  SphP[p].i.Injected_BH_Momentum[0], 
		 SphP[p].i.Injected_BH_Momentum[1], SphP[p].i.Injected_BH_Momentum[2]);
	  printf(" Velocity = (x,y,z) = (%g, %g, %g) \n",  P[p].Vel[0], P[p].Vel[1], P[p].Vel[2]);
	  printf(" previous time step %g \n", (P[p].TimeBin ? (1 << P[p].TimeBin) : 0) * All.Timebase_interval);
//	  SphP[i].i.Injected_BH_Momentum[k] = 0.;
      }
#ifdef COSMIC_RAYS
      if(P[p].Type == 0)
	printf("Cosmic Ray Properties: C0: %g -- q0  : %g -- P  : %g\n"
	       "                       Rho: %g\n",
	       SphP[p].CR_C0, SphP[p].CR_q0, CR_Particle_Pressure(SphP + p), SphP[p].d.Density);
#endif

      fflush(stdout);
      endrun(818);
    }

  return ti_step;
}


/*! This function computes an upper limit ('dt_displacement') to the global timestep of the system based on
 *  the rms velocities of particles. For cosmological simulations, the criterion used is that the rms
 *  displacement should be at most a fraction MaxRMSDisplacementFac of the mean particle separation. Note that
 *  the latter is estimated using the assigned particle masses, separately for each particle type. If comoving
 *  integration is not used, the function imposes no constraint on the timestep.
 */
void find_dt_displacement_constraint(double hfac /*!<  should be  a^2*H(a)  */ )
{
  int i, type;
  int count[6];
  long long count_sum[6];
  double v[6], v_sum[6], mim[6], min_mass[6];
  double dt, dmean, asmth = 0;

  dt_displacement = All.MaxSizeTimestep;

  if(All.ComovingIntegrationOn)
    {
      for(type = 0; type < 6; type++)
	{
	  count[type] = 0;
	  v[type] = 0;
	  mim[type] = 1.0e30;
	}

      for(i = 0; i < NumPart; i++)
	{
	  v[P[i].Type] += P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2];
	  if(P[i].Mass > 0)
	    {
	      if(mim[P[i].Type] > P[i].Mass)
		mim[P[i].Type] = P[i].Mass;
	    }
	  count[P[i].Type]++;
	}

      MPI_Allreduce(v, v_sum, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(mim, min_mass, 6, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

      sumup_large_ints(6, count, count_sum);

#ifdef SFR
      /* add star and gas particles together to treat them on equal footing, using the original gas particle
         spacing. */
      v_sum[0] += v_sum[4];
      count_sum[0] += count_sum[4];
      v_sum[4] = v_sum[0];
      count_sum[4] = count_sum[0];
#ifdef BLACK_HOLES
      v_sum[0] += v_sum[5];
      count_sum[0] += count_sum[5];
      v_sum[5] = v_sum[0];
      count_sum[5] = count_sum[0];
      min_mass[5] = min_mass[0];
#endif
#endif

      for(type = 0; type < 6; type++)
	{
	  if(count_sum[type] > 0)
	    {
	      if(type == 0 || (type == 4 && All.StarformationOn))
		dmean =
		  pow(min_mass[type] / (All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
		      1.0 / 3);
	      else
		dmean =
		  pow(min_mass[type] /
		      ((All.Omega0 - All.OmegaBaryon) * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
		      1.0 / 3);

#ifdef BLACK_HOLES
	      if(type == 5)
		dmean =
		  pow(min_mass[type] / (All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
		      1.0 / 3);
#endif
	      dt = All.MaxRMSDisplacementFac * hfac * dmean / sqrt(v_sum[type] / count_sum[type]);

#ifdef PMGRID
	      asmth = All.Asmth[0];
#ifdef PLACEHIGHRESREGION
	      if(((1 << type) & (PLACEHIGHRESREGION)))
		asmth = All.Asmth[1];
#endif
	      if(asmth < dmean)
		dt = All.MaxRMSDisplacementFac * hfac * asmth / sqrt(v_sum[type] / count_sum[type]);
#endif

	      if(ThisTask == 0)
		printf("type=%d  dmean=%g asmth=%g minmass=%g a=%g  sqrt(<p^2>)=%g  dlogmax=%g\n",
		       type, dmean, asmth, min_mass[type], All.Time, sqrt(v_sum[type] / count_sum[type]), dt);


#ifdef NEUTRINOS
	      if(type != 2)	/* don't constrain the step to the neutrinos */
#endif
		if(dt < dt_displacement)
		  dt_displacement = dt;
	    }
	}

      if(ThisTask == 0)
	printf("displacement time constraint: %g  (%g)\n", dt_displacement, All.MaxSizeTimestep);
    }
}

int get_timestep_bin(int ti_step)
{
  int bin = -1;

  if(ti_step == 0)
    return 0;

  if(ti_step == 1)
    {
      printf("time-step of integer size 1 not allowed\n");
      endrun(112313);
    }

  while(ti_step)
    {
      bin++;
      ti_step >>= 1;
    }

  return bin;
}
