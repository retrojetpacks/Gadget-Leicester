#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"


/*! \file accel.c
 *  \brief driver routines to carry out force computation
 */


/*! This routine computes the accelerations for all active particles.  First, the gravitational forces are
 * computed. This also reconstructs the tree, if needed, otherwise the drift/kick operations have updated the
 * tree to make it fullu usable at the current time.
 *
 * If gas particles are presented, the `interior' of the local domain is determined. This region is guaranteed
 * to contain only particles local to the processor. This information will be used to reduce communication in
 * the hydro part.  The density for active SPH particles is computed next. If the number of neighbours should
 * be outside the allowed bounds, it will be readjusted by the function ensure_neighbours(), and for those
 * particle, the densities are recomputed accordingly. Finally, the hydrodynamical forces are added.
 */
void compute_accelerations(int mode)
{
#ifdef RADTRANSFER
  int iter = 0;
#endif

#if defined(BUBBLES) || defined(MULTI_BUBBLES)
  double hubble_a;
#endif

  if(ThisTask == 0)
    {
      printf("Start force computation...\n");
      fflush(stdout);
    }

#ifdef REIONIZATION
  heating();
#endif

  CPU_Step[CPU_MISC] += measure_time();

#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current)
    {
      long_range_force();

      CPU_Step[CPU_MESH] += measure_time();
    }
#endif


#ifndef ONLY_PM

  gravity_tree();		/* computes gravity accel. */

  if(All.TypeOfOpeningCriterion == 1 && All.Ti_Current == 0)
    gravity_tree();		/* For the first timestep, we redo it
				 * to allow usage of relative opening
				 * criterion for consistent accuracy.
				 */
#endif


#ifdef FORCETEST
  gravity_forcetest();
#endif


  if(All.TotN_gas > 0)
    {
      /***** density *****/
      if(ThisTask == 0)
	{
	  printf("Start density computation...\n");
	  fflush(stdout);
	}
      density();		/* computes density, and pressure */

#if defined(MAGNETIC) && defined(BSMOOTH)
      bsmooth();
#endif
#if (defined(CONDUCTION) || defined(CR_DIFFUSION) || defined(SMOOTH_PHI))
      compute_smoothed_values();
#endif



      /***** update smoothing lengths in tree *****/
      force_update_hmax();


      /***** hydro forces *****/
      if(ThisTask == 0)
	{
	  printf("Start hydro-force computation...\n");
	  fflush(stdout);
	}

      hydro_force();		/* adds hydrodynamical accelerations  and computes du/dt  */

#ifdef RADTRANSFER
      /***** compute eddington tensor *****/
      if(ThisTask == 0)
	{
	  printf("Start Eddington tensor computation...\n");
	  fflush(stdout);
	}

      eddington();

      if(ThisTask == 0)
	{
	  printf("%s\n", "done Eddington tensor!");
	  fflush(stdout);
	}

      star_density();

      /***** set simple initial conditions *****/
      if(All.Time == All.TimeBegin)
	{
	  if(ThisTask == 0)
	    {
	      printf("Setting simple inits...\n");
	      fflush(stdout);
	    }

	  set_simple_inits();

	  if(ThisTask == 0)
	    {
	      printf("%s\n", "done with simple inits!");
	      fflush(stdout);
	    }
	}

      /***** evolve the transport of radiation *****/
      if(ThisTask == 0)
	{
	  printf("start radtransfer...\n");
	  fflush(stdout);
	}

      do
	{
	  radiative_transfer();
	  iter++;
	  if(ThisTask == 0)
	    printf("%s %f\n", "the residue is ", All.Residue);
	}
      while(All.Residue > 0.0001);

      update_nHI();
      simple_output();

      if(ThisTask == 0)
	{
	  printf("%s \n %i %s\n", "done with radtransfer!", iter, "iterations in total!");
	  fflush(stdout);
	}
#endif

#ifdef MHM
      /***** kinetic feedback *****/
      kinetic_feedback_mhm();
#endif


#ifdef BLACK_HOLES
      /***** black hole accretion and feedback *****/
#ifdef FOF
      /* this will find new black hole seed halos */
      if(All.Time >= All.TimeNextBlackHoleCheck)
	{
	  fof_fof(-1);

	  if(All.ComovingIntegrationOn)
	    All.TimeNextBlackHoleCheck *= All.TimeBetBlackHoleSearch;
	  else
	    All.TimeNextBlackHoleCheck += All.TimeBetBlackHoleSearch;
	}
#endif
#ifndef NO_BH_ACCRETION 
      blackhole_accretion();
#endif
#ifdef DUST 
      dust_density();
#ifdef DUST_TWO_POPULATIONS
      dust_density_on_gas();
#endif
      dust_drag();
#ifdef DUST_SINK_ON_FLY
      dust_sink_formation();
#endif
#ifdef DUST_PEBBLES_BORN
      pebbles_born();
      // defined(DUST_TWO_POPULATIONS) && defined(DUST_PEBBLES_BORN)

#endif



#endif

#if defined(SGRA_POTENTIAL) || defined(CUSP_POTENTIAL) || defined(SIS_POTENTIAL) || defined(NFW_POTENTIAL) || defined(QUASAR_HEATING) || defined(OTHIN_ACCELERATOR) || defined(FIND_SMBH)
/* share position of the SMBH (currently only one allowed! 02.08.09) between
 * the different processors */
      FindQuasars();
#endif

#if defined(PLANET_IRRADIATION) || defined(PLANET_ACCRETION_FEEDBACK) || defined(PLANET_ACCRETION_REPORTING)
      FindPlanet();
#endif


/* This calls up rad transfer or momentum transfer via feedback particles */
#if defined(VIRTUAL) 

#if defined(FB_MOMENTUM)
/* wind particles */      
      momentum_feedback();
#endif

#ifdef PHOTONS
#ifdef OTHIN_ACCELERATOR
/* opticall thin heating assumption */
  sink_heating_othin();
#else 
/* photons via onte Carlo */
  fb_particles(); 
#endif
#endif /* PHOTONS */

/* #ifdef FB_MOMENTUM */
/* /\* This deals with momentum-transferring particles only, e.g., wind *\/ */
/*   fb_momentum();  */
/* #endif */

#endif /* VIRTUAL */

#endif /* BLACK_HOLES */


#ifdef COOLING
      /**** radiative cooling and star formation *****/
      cooling_and_starformation();

      CPU_Step[CPU_COOLINGSFR] += measure_time();
#endif


#ifdef CR_DIFFUSION_GREEN
      greenf_diffusion();
#endif


#ifdef BUBBLES
      /**** bubble feedback *****/
      if(All.Time >= All.TimeOfNextBubble)
	{
#ifdef FOF
	  fof_fof(-1);
	  bubble();
#else
	  bubble();
#endif
	  if(All.ComovingIntegrationOn)
	    {
              hubble_a = hubble_function(All.Time);
	      All.TimeOfNextBubble *= (1.0 + All.BubbleTimeInterval * hubble_a);
	    }
	  else
	    All.TimeOfNextBubble += All.BubbleTimeInterval / All.UnitTime_in_Megayears;

	  if(ThisTask == 0)
	    printf("Time of the bubble generation: %g\n", 1. / All.TimeOfNextBubble - 1.);
	}
#endif


#if defined(MULTI_BUBBLES) && defined(FOF)
      if(All.Time >= All.TimeOfNextBubble)
	{
	  fof_fof(-1);

	  if(All.ComovingIntegrationOn)
	    {
              hubble_a = Hubble_func(All.Time);
	      All.TimeOfNextBubble *= (1.0 + All.BubbleTimeInterval * hubble_a);
	    }
	  else
	    All.TimeOfNextBubble += All.BubbleTimeInterval / All.UnitTime_in_Megayears;

	  if(ThisTask == 0)
	    printf("Time of the bubble generation: %g\n", 1. / All.TimeOfNextBubble - 1.);
	}
#endif

    }

  if(ThisTask == 0)
    {
      printf("force computation done.\n");
      fflush(stdout);
    }
}
