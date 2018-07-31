#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"
#include "forcetree.h"


#if defined(VIRTUAL)  && defined(PHOTONS)

/* Aug 2007: SN
 * This routine creates new feedback particles on the base of sfr_eff.c
 */

void fb_particles(void)
/* cooling routine when star formation is enabled */
{
#ifdef DUST_OPACITY_FIT     
  double const5; /* kappa(T) = kappa0 (T/Tscale); kappa0 is from the insput file */
#endif




  unsigned int bits;
  int i, j, k, prev, bin, flag, stars_spawned, tot_spawned, new_virt, new_wind, N_iter,
     N_active, N_active_tot;
  struct state_of_system sys;
  double ascale = 1, hubble_a = 0, a3inv, correction_factor, wemit, correction_wind;
  double time_hubble_a, dir[3], theta, phi, norm, dt, dt_create, rseed, dtau, r, new_virt_exact, new_wind_exact;
  double CurrentLuminosity, EradFactor, RadEnergyDensity, u_eq, acmax, rho, 
      dtime, u_to_temp_fac, temp, uold, unew, dmax1, dmax2, energy_gamma, energy_wind_min, ax, ay, az, ac,
      updv, dudt_pdv, u1, u2, rhs, lhs, const4, diff, dr_jump, RadHeating;
  double AllowedJump = 0.5; /* fraction of lambda_orin that a photon can jump in Monte Carlo regime */

  double lambda_orig, desired_dr, desired_remainder, 
      dt_bh, dt_total, r2, dx, dy, dz, dt_tot_all;

  EradFactor = 1./THOMPSON/All.VirtualCrosSection * PROTONMASS  * (All.UnitVelocity_in_cm_per_s/C) *
      All.UnitVelocity_in_cm_per_s/All.UnitTime_in_s;

#ifdef CONSTANT_MEAN_MOLECULAR_WEIGHT
  u_to_temp_fac = All.MeanWeight  * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
#else
  u_to_temp_fac = 2.3  * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
#endif

/* SN: Multiplying this by particle mass in code units gives Eddington limit
 * in code units, defined using Thompson opacity */
  double L_Edd_dimensionless = (4 * M_PI * GRAVITY * C * PROTONMASS/ THOMPSON)/pow(All.UnitVelocity_in_cm_per_s, 2) *
      All.UnitTime_in_s;
/* NOTE THAT THIS IS DEFINED FOR EDD LIMIT FOR UNIT MASS!!! IF TOTAL MASS >> 1 THEN NEED TO RESET THIS */
#ifdef PHOTON_ENERGY_LEDD   
/* photon energy limited by the fraction of the eddington limit for the unit
 * mass (assumed = 1). Convenient when running star formation simulations
 * where the total blackhole/sink particle mass will be limited by the initial
 * gas mass */
  double energy_gamma_min = 1.* L_Edd_dimensionless * All.VirtualTime/(All.MaxPart);
  energy_wind_min  = 1. * energy_gamma_min; 
#else 
/* photon energy limited by the fraction of the SPH momentum */
  double energy_gamma_min = 10. * All.VirtualMomentum * C/All.UnitVelocity_in_cm_per_s;
  energy_wind_min  =  energy_gamma_min; 
#endif
  energy_gamma = 0.1 * All.VirtualMomentum * C/All.UnitVelocity_in_cm_per_s;
//L_Edd_dimensionless * All.VirtualTime/(All.MaxPart*NTask);

  CPU_Step[CPU_MISC] += measure_time();


/* number of steps over which photons are initialised at t=0 */

  if(ThisTask == 0)
    {
      printf("Beginning Feed-Back Particles\n");
      fflush(stdout);
    }

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

  stars_spawned = 0;

/* THIS BLOCK GOES OVER ALL ACTIVE BLACKHOLE PARTICLES. THE SPH PARTICLES AND
 * PHOTONS ARE TREATED IN ANOTHER BLOCK MARKED BY !!--!! BELOW
 */

  double added_energy_by_bh = 0, total_added_energy_by_bh = 0;

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
  {

/* Creation of new virtual/photon particles
 */	    
      if(P[i].Type == 5 && P[i].Mass > 0.5*All.SMBHmass && All.Time > All.VirtualStart) /* do this for the more massive BH only */
      {
	  flag = 1;
/* decide how many new virtual particles are needed */
	  dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / hubble_a;
	  P[i].dt1 = dt;
	  dt_create = dt;

#ifdef FRACTION_OF_LEDD_FEEDBACK
	  new_virt_exact = All.VirtualFeedBack * L_Edd_dimensionless * P[i].Mass * dt_create/energy_gamma_min;
	  new_wind_exact = All.VirtualFeedBack * L_Edd_dimensionless * P[i].Mass * dt_create/energy_wind_min;
#else
/* SN -- Feedback appropriate for stellar radiation sources */
/* The first bit is the accretion luminosity, assuming GM/R_core = Msink/ 3Rsun */
/* The second bit is the main sequence luminosity of the star */
	  CurrentLuminosity = (P[i].Mass-0.01)/(3.*6.96e10/All.UnitLength_in_cm) * P[i].BH_Mdot;
	  if (CurrentLuminosity <= 0) CurrentLuminosity = 0;
          CurrentLuminosity += SOLAR_LUM/All.UnitEnergy_in_cgs * All.UnitTime_in_s * pow(P[i].Mass * All.UnitMass_in_g/SOLAR_MASS, 3.3);
/*  Now limit that at Eddington luminosity */
//	  if (CurrentLuminosity > L_Edd_dimensionless * P[i].Mass) CurrentLuminosity = L_Edd_dimensionless * P[i].Mass;

//	  CurrentLuminosity = 2. * SOLAR_LUM/All.UnitEnergy_in_cgs * All.UnitTime_in_s;


	  new_virt_exact = All.VirtualFeedBack * CurrentLuminosity * dt_create/energy_gamma_min;
	  new_wind_exact = All.VirtualFeedBack * CurrentLuminosity * dt_create/energy_wind_min;

#endif

#ifndef FB_MOMENTUM 
	  new_wind_exact = 0.;
#endif

#ifndef VIRTUAL_RE_EMISSION /* don't want any photons in this case */
//	  new_virt_exact = 0.;
#endif


#if  defined(OTHIN_ACCELERATOR) 
	  new_virt_exact = 0.;
#endif
//	  new_virt_exact = 0.;

	  printf("Accretion Rate %g, Lum %g, new_virt_exact %g, new_wind_exact %g, Bh mass %g, dt %g \n", 
		 P[i].BH_Mdot, CurrentLuminosity*All.UnitEnergy_in_cgs/All.UnitTime_in_s, new_virt_exact, new_wind_exact, P[i].Mass, dt);
	  fflush(stdout);  

/* Deciding how many NEW PHOTONS excatly needed */
#ifdef STARTUP_LUMINOSITY
	  if (new_virt_exact > 1.e-4) 
#else
	  if (new_virt_exact > 0.) 
#endif
	  {
	      new_virt = (int) new_virt_exact + 1;
	  }
	  else
	  {
	      new_virt = 0;
	  }
	  if (new_virt > 0)
	  {
	      correction_factor = new_virt_exact/new_virt;
	  }
	  else
	  {
	      correction_factor = 1.;
	  }

/* Deciding how many NEW WIND PARTICLES excatly needed */
#ifdef STARTUP_LUMINOSITY
	  if (new_wind_exact > 1.e-3) 
#else
	  if (new_wind_exact > 0.) 
#endif
	  {
	      new_wind = (int) new_wind_exact + 1;
	  }
	  else
	  {
	      new_wind = 0;
	  }
	  if (new_wind > 0)
	  {
	      correction_wind = new_wind_exact/new_wind;
	  }
	  else
	  {
	      correction_wind = 1.;
	  }
	  
	  
 	  if (new_virt > 0 || new_wind > 0)
 	  {
/*  		printf("On Task=%d Need to spawn %d virtual particles, new_virt_exact is equal to %g \n", ThisTask, new_virt, new_virt_exact);
  		fflush(stdout);
*/
 	      flag=0;
 	  }
	  
	  if(flag == 0)		/* active virtual particle creation */
	  {
	      for(bits = 0; bits < new_virt + new_wind; bits++)
	      {
		  /* here we spawn a new star particle */
		  if(NumPart + stars_spawned >= All.MaxPart)
		  {
		      printf
			  ("On Task=%d with NumPart=%d we try!! to spawn %d photons. Sorry, no space left...(All.MaxPart=%d)\n",
			   ThisTask, NumPart, stars_spawned, All.MaxPart);
		      fflush(stdout);
		      endrun(8888);
		  }
		  
		  P[NumPart + stars_spawned] = P[i];
		  P[NumPart + stars_spawned].Type = 3;
		  P[NumPart + stars_spawned].NewDensity = P[i].b1.dBH_Density;
		  
#ifdef STELLARAGE
/* StellarAge is time when the particle was created  */
		  P[NumPart + stars_spawned].StellarAge = All.Time;
#endif
		  
//SN -- adding kicks to particles ala Cuadra et al SPH particles wind method
#ifdef ISOTROPIC
//		  theta = acos(2. * get_random_number(P[i].ID + 2 + bits + All.NumCurrentTiStep) - 1.);
		  theta= acos(2. *  gsl_rng_uniform(random_generator) - 1.);
#else
		  theta = acos(0.5 * get_random_number(P[i].ID + 2 + bits + All.NumCurrentTiStep) + 0.5);
#endif
//		  phi = 2 * M_PI * (get_random_number(P[i].ID + 3 + bits + All.NumCurrentTiStep) + All.Time);
		  phi = 2 * M_PI * gsl_rng_uniform(random_generator);
		  
		  dir[0] = sin(theta) * cos(phi);
		  dir[1] = sin(theta) * sin(phi);
		  dir[2] = cos(theta);

		  
		  /* this gives random velocity to the
		     winds */
		  for(j = 0, norm = 0; j < 3; j++)
		      norm += dir[j] * dir[j];
		  
		  norm = sqrt(norm);
/*		    if(get_random_number(P[i].ID + 5 + bits + All.NumCurrentTiStep ) < 0.5)
		    norm = -norm;
		    This is redundant given that I create partocles random in 4 pi angle already */ 		    
		  
		  if(norm != 0)
		  {
		      for(j = 0; j < 3; j++)
			  dir[j] /= norm;
		      
		      for(j = 0; j < 3; j++)
		      {
/* particles are kicked in real and velocity space in direction dir[j] */

//			  P[NumPart + stars_spawned].Vel[j] = All.FeedBackVelocity * dir[j];
			  P[NumPart + stars_spawned].Vel[j] = C/All.UnitVelocity_in_cm_per_s * dir[j]/All.FeedBackVelocity;
			  
/*	 initialise photons with the maximum Virtual Hsml = search radius for
 * sph particles  */
//			  P[NumPart + stars_spawned].Hsml = 3.;
			  P[NumPart + stars_spawned].Hsml += All.SofteningBndryMaxPhys;
#ifdef ACCRETION_RADIUS
			  P[NumPart + stars_spawned].Hsml += All.InnerBoundary;
#endif

			  PPP[NumPart + stars_spawned].Hsml = P[NumPart + stars_spawned].Hsml;
			  
			  
/* 			  P[NumPart + stars_spawned].Pos[j] += gsl_rng_uniform(random_generator) */
/* 			      * All.FeedBackVelocity * dt * dir[j] + (0. * All.SofteningBulgeMaxPhys + (0.25 + 0.25*gsl_rng_uniform(random_generator)) */
/* 				  * P[i].Hsml) * dir[j]; */
			      P[NumPart + stars_spawned].Pos[j] += 0.2*(0.25 + 0.25*gsl_rng_uniform(random_generator))
				  * P[i].Hsml * dir[j];

#ifdef VIRTUAL_FLY_THORUGH_EMPTY_SPACE
			  P[NumPart + stars_spawned].b1.BH_Density = 0.;
#endif
		      }
		  }
		  
		  NextActiveParticle[NumPart + stars_spawned] = FirstActiveParticle;
		  FirstActiveParticle = NumPart + stars_spawned;
		  NumForceUpdate++;
		  TimeBinCount[P[NumPart + stars_spawned].TimeBin]++;
		  
		  PrevInTimeBin[NumPart + stars_spawned] = i;
		  NextInTimeBin[NumPart + stars_spawned] = NextInTimeBin[i];
		  if(NextInTimeBin[i] >= 0)
		      PrevInTimeBin[NextInTimeBin[i]] = NumPart + stars_spawned;
		  NextInTimeBin[i] = NumPart + stars_spawned;
		  if(LastInTimeBin[P[i].TimeBin] == i)
		      LastInTimeBin[P[i].TimeBin] = NumPart + stars_spawned;
		  
		  P[NumPart + stars_spawned].ID += NumPart + stars_spawned;
		  P[NumPart + stars_spawned].NewDensity = 0.;
#ifdef NON_GRAY_RAD_TRANSFER
/* photons from the sink are born with high temperature, resulting in high opacity */
		  P[NumPart + stars_spawned].VirtualTemperature = 500.;
/* opacity is approximately kappa ~ kappa_0 (T/T0) */
		  P[NumPart + stars_spawned].Kappa = P[NumPart + stars_spawned].VirtualTemperature/All.EqTemp;
#endif
		  
/* DeltaPhotonMomentum us the change in the OldPhotonMomentum in a single SPH/photon interaction */
		  P[NumPart + stars_spawned].DeltaPhotonMomentum = 0.;
#if defined(VIRTUAL_HEATING) 
		  P[NumPart + stars_spawned].DeltaHeat = 0.;
		  P[NumPart + stars_spawned].DeltaPhotonEnergy = 0.;
#endif
/* initialise photon momentum */
//bits < new_virt + new_wind;
		  if (new_virt > 0 && bits < new_virt)
		  {
		      P[NumPart + stars_spawned].WindOrPhoton = 1; /* This marks it as a photon */
		      P[NumPart + stars_spawned].OldPhotonMomentum = energy_gamma_min/C* All.UnitVelocity_in_cm_per_s * correction_factor;
		      P[NumPart + stars_spawned].BirthPhotonEnergy = energy_gamma_min * correction_factor;
		  }
		  

		  if (new_wind > 0 && bits >= new_virt)
		  {
		      P[NumPart + stars_spawned].WindOrPhoton = 0; /* This marks it as a Wind particle */
		      P[NumPart + stars_spawned].OldPhotonMomentum = energy_wind_min/C* All.UnitVelocity_in_cm_per_s * correction_wind;
		      P[NumPart + stars_spawned].BirthPhotonEnergy = energy_wind_min * correction_wind;

		  }

                  added_energy_by_bh += P[NumPart + stars_spawned].BirthPhotonEnergy;
		  P[NumPart + stars_spawned].Mass = All.VirtualMass * P[NumPart + stars_spawned].OldPhotonMomentum;
#ifndef EXCLUDE_VIRTUAL_FROM_TREE
		  force_add_star_to_tree(i, NumPart + stars_spawned);
#endif
		  
		  stars_spawned++;
	      }
	  }
      }
  }				/* end of main loop over active particles */
  

  MPI_Allreduce(&stars_spawned, &tot_spawned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(tot_spawned > 0)
  {
      All.TotNumPart += tot_spawned;
      NumPart += stars_spawned;
  }
  stars_spawned = tot_spawned = 0.;

  
/*  This re-constructs the NextPhoton list */
  FirstPhoton = -1;  
  prev = -1;
//  for(i = 0; i < NumPart; i++) 
  for(i = FirstActiveParticle ; i >= 0; i = NextActiveParticle[i])
  { 
/* Note: photons that were just created are inactive until timestep.c assigns them a time step */
      if(P[i].Type == 3 &&  P[i].StellarAge < All.Time)
      {
	  if(prev == -1)
	  {
	      FirstPhoton = i;
	  }
	  
 	  if(prev >= 0) 
 	      NextPhoton[prev] = i; 
	  
	  prev = i;
      }
   } 
  
  if(prev >= 0)
  {
      NextPhoton[prev] = -1;
  }
  FirstActivePhoton = FirstPhoton;
  NextActivePhoton = NextPhoton;

  bits = 0;
  dt_total = 0.;
  for(i = FirstActivePhoton ; i >= 0; i = NextActivePhoton[i])
  {
/* This is the time the photon needs to be propagated for. */
      P[i].DtJump = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / hubble_a;
      dt_bh = P[i].DtJump;
      dt_total += dt_bh;
      P[i].dt2 = 0.; /* running time of the photon used in the iterative loop
			    * below. Runs from 0 to DtJump */
  }
  if (ThisTask == 0) 
  {
      printf("Step number %d, dt_total = %g \n", All.NumCurrentTiStep, dt_total);	     
      fflush(stdout);
  }
  N_iter = 0;
  int nendd = 1;
  MPI_Allreduce(&dt_total, &dt_tot_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  if (dt_tot_all <= 1.e-15 * All.VirtualTime * (All.MaxPart*NTask)) nendd = -1; 
  if (dt_tot_all <= 0.e-14) nendd = -1; /* don't even enter the loop in this case */

  if (nendd > 0 ) virtual_density();

  while (nendd >= 0)
  {

      N_iter++;

      for(i = FirstActivePhoton ; i >= 0; i = NextActivePhoton[i])
      {

/* remaining time for the photon to be active*/
	  dt_bh = P[i].DtJump;
/* if it is less than 0 the photon is inactive; move to the next photon */
	  if (dt_bh <= 0.) // || P[i].Mass == 0)
	  {
	      dt_bh = P[i].DtJump = 0.;
	      continue;
	  }
/* a counter for basic steps for all the photons in this loop */
	  bits++;
	  rho = P[i].NewDensity + 1.e-38;
/* original, i.e., physical, mean free path of the photon */
	  lambda_orig = 1./(All.VirtualCrosSection * THOMPSON/PROTONMASS * rho 
			    *  All.UnitDensity_in_cgs * All.UnitLength_in_cm);
	  if (rho*All.UnitDensity_in_cgs > 1.e-13) lambda_orig *= pow(10.,-5.2)*pow(rho*All.UnitDensity_in_cgs,-0.4);
#ifdef NON_GRAY_RAD_TRANSFER
	  if (P[i].VirtualTemperature > 0.)
	  lambda_orig /= P[i].Kappa;
#endif
/* Desired spatial jump in one step for the photon. Limited by spatial
 * resolution requirements */
	  desired_dr = 1./2. * 0.25 * P[i].Hsml; 
#ifdef VIRTUAL_FLY_THORUGH_EMPTY_SPACE
/* If there are no SPH particles containing the photon in its Hsml, but there
 * are SPH particles in the search radius of the photon, be careful. I the
 * photon desired_dr is too large it propagates too far in one step, skipping
 * lots of sph particles that it should have interacted with if it were
 * propagated by smaller desired_dr. BH_Density is actually the minimum
 * distance to an SPH particle within the Hsml, so this provides a convenient
 * way of allowing small enough steps */
//	  if (P[i].NewDensity <= 0 && P[i].b1.BH_Density <= P[i].Hsml)
	  if (desired_dr > 0.5 * P[i].b1.BH_Density + 2. * All.SofteningGas)
	  {
	      desired_dr = 1./2. * (0.5 * P[i].b1.BH_Density + 2. * All.SofteningGas);
	  }
#endif


/* if desired spatial jump takes longer than the requested time-jump to
 * inactive status, reduce dr accordingly */
	  if (desired_dr > C*dt_bh/All.UnitVelocity_in_cm_per_s/All.FeedBackVelocity) 
	      desired_dr = C*dt_bh/All.UnitVelocity_in_cm_per_s/All.FeedBackVelocity;

/* 	  if (ThisTask == 0)//  && P[i].DeltaHeat <= 0) */
/* 	  { */
/* 	      fprintf (photonSN,"%d, before %d, time %g, (x,y,z)  %g, %g, %g, den = %g, des_dr %g, hsml %g, DtJump %g \n", i, N_iter, */
/* 		       All.Time, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], P[i].NewDensity, desired_dr, P[i].Hsml, P[i].DtJump);/\*, */
/* 									  P[i].DeltaHeat, P[i].OldPhotonMomentum, P[i].DeltaPhotonEnergy);*\/ */
/* 	      fflush(photonSN); */
/* 	  } */




/* Now decide whether the jump will be in the standard linear regime or the
 * Accelerated Diffusion one */ 

	  if (desired_dr < All.FreeTravelLength * lambda_orig)
	  {
/* Standard case -- photon simply travels in straight line. The time
 * step is then simply the light crossing time of this distance  */
/* HOWEVER, the distance travelled must not exceed a fraction (AllowedJump) of
 * the mean free path or else the integration becomes inaccurate, thus: */
	      dr_jump = desired_dr;
	      if (dr_jump > AllowedJump * lambda_orig)
	      {
		  dr_jump = AllowedJump * lambda_orig;
	      }
/* Time it takes for this jump */
	      P[i].dt1 = dr_jump * All.UnitVelocity_in_cm_per_s/C * All.FeedBackVelocity;


/* Here I need to calculate delta tau, and the absorption of the photon's
 * momentum, and energy */

	      dtau = dr_jump/lambda_orig;
	      
#ifndef VIRTUAL_OTHIN
	      if (dtau > 1.e-5)
	      {
		  P[i].DeltaPhotonMomentum = P[i].OldPhotonMomentum * (1. - exp(-dtau));
	      }
	      else
	      {
		  P[i].DeltaPhotonMomentum = P[i].OldPhotonMomentum * dtau;
	      }     
#else
/* if the test is optically thin, this should be exactly dtau! */
	      P[i].DeltaPhotonMomentum = P[i].OldPhotonMomentum * dtau;
#endif
#if defined(VIRTUAL_HEATING)
	      P[i].DeltaHeat = P[i].DeltaPhotonMomentum * C/All.UnitVelocity_in_cm_per_s;
#endif
#ifndef VIRTUAL_OTHIN
/* in an optically thin test, photon momentum is not changed */
//	      P[i].OldPhotonMomentum -= P[i].DeltaPhotonMomentum;
#endif
/* P[i].Mass is equal to P[i].OldPhotonMomentum save for the constant in
 * front. So it can be used to track photon momentum in snapshot files.  */
	      P[i].Mass = All.VirtualMass * P[i].OldPhotonMomentum;

/* Now make the jump */
	      for (k = 0 ; k < 3; k++)
		  P[i].Pos[k] += P[i].dt1 * P[i].Vel[k];
		
	      if (ThisTask < -1)
	      {	  
		  dx =  P[i].Pos[0];
		  dy =  P[i].Pos[1];
		  dz =  P[i].Pos[2];
		  double dvx, dvy, dvz;
		  dvx =  P[i].Vel[0];
		  dvy =  P[i].Vel[1];
		  dvz =  P[i].Vel[2];

		  double rp2 = dx * dx + dy * dy + dz * dz;
		  double vp2 = dvx * dvx + dvy * dvy + dvz * dvz;
		  fprintf(photonSN,"%d, %d, time %g, desired_dr %g, dt_bh %g, r %g, dt1 %g, dr %g, lambda %g, density %g, Hsml = %g \n", i, N_iter,
			  All.Time, desired_dr, dt_bh, sqrt(rp2), 
			  P[i].dt1, sqrt(vp2)*P[i].dt1, lambda_orig, P[i].NewDensity, P[i].Hsml);
		  fflush(photonSN);
	      }



		  
/* reduce the amount of time remaining till the photon becomes inactive */
	      dt_bh -= P[i].dt1;
	      P[i].DtJump -= P[i].dt1;
	      dt_total -= P[i].dt1;	      
/* This is how long the photon has been active */
	      P[i].dt2 += P[i].dt1;
	  }
	  else
	  {
/* optically thick case: Acclerated diffusion here */	  
	      P[i].dt1 = desired_dr*desired_dr/(2*lambda_orig * C/All.UnitVelocity_in_cm_per_s);

/* CAREFULL: THERE SHOULD BE A PROTECTION CONDITION AGAINST THE PHOTON TAKING
 * TOO LONG TO MAKE THIS JUMP. IF dt1 is longer than dt_bh, the photon should
 * not be allowed to jump by this much.
 */

//	      desired_remainder = sqrt(2. * lambda_orig * C/All.UnitVelocity_in_cm_per_s * dt_bh);
/* The factor of 2 comes from the fact that after N scattering the
 * displacement r^2 = 2 N lambda^2 = 2 lambda c t */
	      
	      /* If photon is re-oriented randomly it means that its
	       * initial momentum is completely absorbed by the gas */
	      P[i].DeltaPhotonMomentum = P[i].OldPhotonMomentum;
#if defined(VIRTUAL_HEATING)
	      P[i].DeltaHeat = P[i].OldPhotonMomentum * C/All.UnitVelocity_in_cm_per_s * desired_dr*desired_dr/
		  (2*lambda_orig*lambda_orig); 
#endif
/* the energy of the new photon is now reduced if SPH particles absorb energy
 * from the radiation field */
//		  P[i].OldPhotonMomentum *= P[i].DeltaPhotonEnergy;


/* make the jump */
/* jenerate new random velocity for 2nd jump */
	      theta = acos(2. * gsl_rng_uniform(random_generator) - 1.);
	      phi = 2 * M_PI * gsl_rng_uniform(random_generator);
// 2 * M_PI * (get_random_number(P[i].ID + 3 + bits + All.NumCurrentTiStep) + All.Time);
	      dir[0] = sin(theta) * cos(phi);
	      dir[1] = sin(theta) * sin(phi);
	      dir[2] = cos(theta);
	      for(j = 0; j < 3; j++)
	      {
		  P[i].Pos[j] +=  desired_dr * dir[j] ;
		  P[i].Vel[j] =  C/All.UnitVelocity_in_cm_per_s * dir[j];
	      }

	      if (ThisTask < 0)
	      {
		  fprintf (photonSN,"%d, %d, %d, %g, %g, %g, %g, %g, density %g, Hsml = %g \n", i, N_iter, 1,
			   All.Time, dt_bh, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], P[i].NewDensity, P[i].Hsml);
		  fflush(photonSN);
	      }
	      dt_bh -= P[i].dt1;
	      P[i].DtJump -= P[i].dt1;
	      dt_total -= P[i].dt1;	 
    
	      P[i].dt1 = 0.;
	      P[i].dt2 = 0.;
//	      dt_bh = 0.;

	  }
	  
/* Re-define photon's Hsml. First find distance from an origin */
	  
	  dx =  P[i].Pos[0];
	  dy =  P[i].Pos[1];
	  dz =  P[i].Pos[2];
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
	  
	  if (P[i].NewDensity > 0) 
	  {
	      PPP[i].Hsml = 3. * pow((All.DesNumNgb * All.OriginalGasMass/(P[i].NewDensity + 1.e-38)), 0.33333) + All.SofteningGasMaxPhys;
	  }
	  else
	  {
	      PPP[i].Hsml *= 1.2;
	  }
	  
#ifdef VARIABLE_PHOTON_SEARCH_RADIUS
/* avoid unphysically large photon search radius */
	  if (PPP[i].Hsml > All.PhotonSearchRadius * sqrt(r2) + All.SofteningBulgeMaxPhys)
	      PPP[i].Hsml = All.PhotonSearchRadius * sqrt(r2) + All.SofteningBulgeMaxPhys;
#endif





	  P[i].Hsml = PPP[i].Hsml;
	  
      }
      
#if defined(VIRTUAL_HEATING) 
      for(i = FirstActivePhoton ; i >= 0; i = NextActivePhoton[i])
      {
	  if (P[i].NewDensity > 0) 
	  {
	      P[i].DeltaPhotonEnergy /= P[i].NewDensity;
#ifdef NON_GRAY_RAD_TRANSFER
/* this now gives average uold for SPH neighbors of this particle */
	      P[i].VirtualTemperature /= P[i].NewDensity; 
#endif
	  }
	  else
	  {
	      P[i].DeltaPhotonEnergy = 1.;
#ifdef NON_GRAY_RAD_TRANSFER
/*  SPH neighbor T is undefined in this case */
	      P[i].VirtualTemperature = 0.; 
#endif
	  }	      
      }
#endif

/* Now transfer the absorbed momentum to SPH particles, and calculation of the
 * heating (Injected_VIRTUAL_Energy) if it is enabled */
/* ------------------------------- */
      virtual_spread_momentum();
/* ------------------------------- */


/*  ------------------------------------
 EMISSION OF DAUGHTER PHOTONS BY PHOTONS.
----------------------------------------- */

#ifndef NO_DAUGHTER_PHOTONS
      
      for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
      {
	  if(P[i].Type != 3 || P[i].Mass <= 0) continue;
      
      
#if defined(VIRTUAL_HEATING)
	  if(P[i].Type == 3 && P[i].DeltaPhotonEnergy > 0 && P[i].DeltaPhotonMomentum > 0 && P[i].WindOrPhoton == 1) 
#else
          if(P[i].Type == 3 && P[i].DeltaPhotonMomentum > 0 && P[i].WindOrPhoton == 1) 
#endif
	  {
	      
/* toss a coin to decide weather a new photon should be emitted */
#if defined(VIRTUAL_HEATING)  

/* /\* DeltaPhotonEnergy is now the energy to be re-emitted from the location of the photon */
	      P[i].DeltaPhotonEnergy *= P[i].DeltaPhotonMomentum * C/All.UnitVelocity_in_cm_per_s;
/* we emit only one new daughter photon with the initial energy equal to the
 * original energy of the parent photon. This is done to prevent photon number
 * exploding exponentially.
 */

	      wemit = gsl_rng_uniform(random_generator);
	      if (wemit <  P[i].DeltaPhotonEnergy/P[i].BirthPhotonEnergy)
	      {
#else /* if we only track photon momentum but not the energy exchange with sph */
              wemit = gsl_rng_uniform(random_generator);
	      if (wemit <  P[i].DeltaPhotonMomentum * C/All.UnitVelocity_in_cm_per_s/P[i].BirthPhotonEnergy) 
	      {
#endif /* defined(VIRTUAL_HEATING) */
		  
/* now emit it */
		  
		  if(NumPart + stars_spawned >= All.MaxPart)
		  {
		      printf
			  ("On Task=%d with NumPart=%d we try! to spawn %d particles. Sorry, no space left...(All.MaxPart=%d)\n",
			   ThisTask, NumPart, stars_spawned, All.MaxPart);
		      printf(" DeltaEn = %g, BirthPhEn = %g \n",P[i].DeltaPhotonEnergy, P[i].BirthPhotonEnergy);
		      fflush(stdout);
		      endrun(8888);
		  }
		  
		  
		  P[NumPart + stars_spawned] = P[i];
		  
		  theta = acos(2. * gsl_rng_uniform(random_generator) - 1.);
		  phi = 2 * M_PI * gsl_rng_uniform(random_generator);
		  dir[0] = sin(theta) * cos(phi);
		  dir[1] = sin(theta) * sin(phi);
		  dir[2] = cos(theta);
		  
		  /* this gives random velocity to the particle  */
		  for(j = 0; j < 3; j++)
		  {
/* particles are kicked in real and velocity space in direction dir[j] */
		      P[NumPart + stars_spawned].Vel[j] = C/All.UnitVelocity_in_cm_per_s* dir[j]/All.FeedBackVelocity;
		  }
		  PPP[NumPart + stars_spawned].Hsml = P[NumPart + stars_spawned].Hsml;
		  NextActiveParticle[NumPart + stars_spawned] = FirstActiveParticle;
		  FirstActiveParticle = NumPart + stars_spawned;
		  NumForceUpdate++;
		  TimeBinCount[P[NumPart + stars_spawned].TimeBin]++;
	      
		  PrevInTimeBin[NumPart + stars_spawned] = i;
		  NextInTimeBin[NumPart + stars_spawned] = NextInTimeBin[i];
		  if(NextInTimeBin[i] >= 0)
		      PrevInTimeBin[NextInTimeBin[i]] = NumPart + stars_spawned;
		  NextInTimeBin[i] = NumPart + stars_spawned;
		  if(LastInTimeBin[P[i].TimeBin] == i)
		      LastInTimeBin[P[i].TimeBin] = NumPart + stars_spawned;
		  
		  P[NumPart + stars_spawned].ID += NumPart + stars_spawned;
//	      P[NumPart + stars_spawned].NewDensity = 0.;
		  
/* PhotonMomentum us the change in the OldPhotonMomentum in a single call to the SPH/photon interaction routine */
/* initialise photon momentum */
		  
		  P[NumPart + stars_spawned].BirthPhotonEnergy = P[i].BirthPhotonEnergy;
		  P[NumPart + stars_spawned].OldPhotonMomentum = P[NumPart + stars_spawned].BirthPhotonEnergy/C * All.UnitVelocity_in_cm_per_s;
		  P[NumPart + stars_spawned].Mass = All.VirtualMass * P[NumPart + stars_spawned].OldPhotonMomentum;
		  P[NumPart + stars_spawned].DeltaPhotonMomentum = 0;
		  P[NumPart + stars_spawned].WindOrPhoton = 1; /* This marks it as a photon */
#ifdef STELLARAGE
/* StellarAge is time when the particle was created  */
		  P[NumPart + stars_spawned].StellarAge = All.Time + P[i].dt2;
#endif
		  P[NumPart + stars_spawned].dt2 = 0.;
		  P[NumPart + stars_spawned].DtJump = P[i].DtJump;
		  P[NumPart + stars_spawned].dt1 = 0.;
/* add this new photon's time to the total propagation time of photons */
		  dt_total += P[i].DtJump;

#ifdef VIRTUAL_FLY_THORUGH_EMPTY_SPACE
		  P[NumPart + stars_spawned].b1.BH_Density = P[i].b1.BH_Density;
#endif

#if defined(VIRTUAL_HEATING) 
//		  P[NumPart + stars_spawned].DeltaPhotonEnergy = 0;
#endif
		  PPP[NumPart + stars_spawned].Hsml = P[NumPart + stars_spawned].Hsml;
#ifndef EXCLUDE_VIRTUAL_FROM_TREE
		  force_add_star_to_tree(i, NumPart + stars_spawned);
#endif
		  stars_spawned++;
		  P[i].DeltaPhotonMomentum = 0;
#if defined(VIRTUAL_HEATING)
		  P[i].DeltaPhotonEnergy = 0;
#endif
	      }
	      else
	      {
		  P[i].DeltaPhotonMomentum = 0;
#if defined(VIRTUAL_HEATING) 
		  P[i].DeltaPhotonEnergy = 0;
#endif
	      }

	  }
      }

#else


      for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
      {
	  if(P[i].Type != 3 || P[i].Mass <= 0) continue;
	  if(P[i].dt2 == 0) continue; /* no need to do anything in the smart diffusion regime */
            
#if defined(VIRTUAL_HEATING)
	  if(P[i].Type == 3 && P[i].DeltaPhotonEnergy > 0 && P[i].DeltaPhotonMomentum > 0 && P[i].WindOrPhoton == 1) 
#else
          if(P[i].Type == 3 && P[i].DeltaPhotonMomentum > 0 && P[i].WindOrPhoton == 1) 
#endif
	  {
	      
/* toss a coin to decide weather a new photon should be emitted */
#if defined(VIRTUAL_HEATING)  

/* /\* DeltaPhotonEnergy is now the energy to be re-emitted from the location of the photon */
	      P[i].DeltaPhotonEnergy *= P[i].DeltaPhotonMomentum * C/All.UnitVelocity_in_cm_per_s;

	      wemit = gsl_rng_uniform(random_generator);
//	      if (wemit <  P[i].DeltaPhotonEnergy*P[i].DeltaPhotonMomentum * C/All.UnitVelocity_in_cm_per_s/P[i].BirthPhotonEnergy)
	      if (wemit <  P[i].DeltaPhotonEnergy/P[i].BirthPhotonEnergy)
	      {
#else /* if we only track photon momentum but not the energy exchange with sph */
              wemit = gsl_rng_uniform(random_generator);
	      if (wemit <  P[i].DeltaPhotonMomentum * C/All.UnitVelocity_in_cm_per_s/P[i].BirthPhotonEnergy) 
	      {
#endif /* defined(VIRTUAL_HEATING) */
		  
/* now emit it */
		  theta = acos(2. * gsl_rng_uniform(random_generator) - 1.);
		  phi = 2 * M_PI * gsl_rng_uniform(random_generator);
		  dir[0] = sin(theta) * cos(phi);
		  dir[1] = sin(theta) * sin(phi);
		  dir[2] = cos(theta);
		  
		  /* this gives random velocity to the particle  */
		  for(j = 0; j < 3; j++)
		  {
/* particles are kicked in velocity space in direction dir[j] */
/*		      P[i].Vel[j] = (P[i].Vel[j] * P[i].BirthPhotonEnergy + 
				     C/All.UnitVelocity_in_cm_per_s* dir[j]/All.FeedBackVelocity * P[i].DeltaPhotonEnergy)/
				     (P[i].BirthPhotonEnergy + P[i].DeltaPhotonEnergy);*/
		      P[i].Vel[j] = C/All.UnitVelocity_in_cm_per_s* dir[j]/All.FeedBackVelocity;
		  }
#ifdef NON_GRAY_RAD_TRANSFER
		  P[i].Kappa = P[i].VirtualTemperature/All.EqTemp;
#endif

#ifdef STELLARAGE
/* StellarAge is time when the particle was created  */
		  P[i].StellarAge = All.Time + P[i].dt2;
#endif
		  P[i].dt2 = 0.;
		  P[i].dt1 = 0.;
/*
		  P[i].BirthPhotonEnergy += P[i].OldPhotonMomentum*C/All.UnitVelocity_in_cm_per_s;// + P[i].DeltaPhotonEnergy;
		  P[i].OldPhotonMomentum = P[i].BirthPhotonEnergy/C * All.UnitVelocity_in_cm_per_s ;
*/
		  P[i].Mass = All.VirtualMass * P[i].OldPhotonMomentum;
#ifdef VIRTUAL_FLY_THORUGH_EMPTY_SPACE
		  P[i].b1.BH_Density = P[i].b1.BH_Density;
#endif

/* 		  if (ThisTask == 0) */
/* 		  { */
/* 		      fprintf (photonSN,"%d, %d, %d,  %g, %g, %g, %g, Birth En %g, Delta Ph E = %g \n", i, N_iter, 1, */
/* 			       All.Time,  P[i].Vel[0], P[i].Vel[1], P[i].Pos[2], P[i].BirthPhotonEnergy, P[i].DeltaPhotonEnergy); */
/* 		      fflush(photonSN); */
/* 		  } */

		  P[i].DeltaPhotonMomentum = 0;
#if defined(VIRTUAL_HEATING) 
		  P[i].DeltaPhotonEnergy = 0;
#endif
		  
	      }
	      else
	      {
		  P[i].DeltaPhotonMomentum = 0;
#if defined(VIRTUAL_HEATING) 
		  P[i].DeltaPhotonEnergy = 0;
#endif
	      }


	  }
      }

#endif /* NO_DAUGHTER_PHOTONS */

/* ELIMINATION OF PHOTONS of too small an energy */

/* check if feedback/virtual particle has left an interesting sphere where I
 * track its motion. Or if it is "too old". In either case, remove it */


	  int count_deleted_photon = 0, count_deleted_photon_total=0;
	  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	  {
	      if (P[i].Type == 3)
	      {
		  r = sqrt(P[i].Pos[0]*P[i].Pos[0] + P[i].Pos[1]*P[i].Pos[1] + P[i].Pos[2]*P[i].Pos[2]);
//		  if (P[i].OldPhotonMomentum < 0.001 * P[i].BirthPhotonEnergy/C * All.UnitVelocity_in_cm_per_s)
		  if (P[i].OldPhotonMomentum < 0.001 * All.VirtualMomentum)
		  {
		      P[i].Mass = 0.;
		      P[i].DtJump = 0.;
		      count_deleted_photon++;
		      P[i].OldPhotonMomentum = 0;
		  }
	      }
	  }


   MPI_Allreduce(&count_deleted_photon, &count_deleted_photon_total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 

// This marks the photons that completed their time jumps as inactive 
       FirstPhoton = -1; 
       prev = -1; 
       N_active = N_active_tot = 0.;
       for(i = FirstActiveParticle ; i >= 0; i = NextActiveParticle[i]) 
       { 

 	  if(P[i].Type == 3 && P[i].DtJump > 1.e-40 && P[i].Mass > 0) 
 	  { 
	      N_active ++;
 	      if(prev == -1) 
 	      { 
 		  FirstPhoton = i; 
 	      } 
	      
 	      if(prev >= 0) 
 		  NextPhoton[prev] = i; 
	      
 	      prev = i; 
 	  } 
       }
      
       if(prev >= 0) 
       { 
 	  NextPhoton[prev] = -1; 
       } 
       FirstActivePhoton = FirstPhoton; 
       NextActivePhoton = NextPhoton; 
      
/* This is needed to share information between processors */
//      MPI_Allreduce(&dt_total, &dt_tot_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&N_active, &N_active_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	  
      if (N_active_tot == 0)   
      {
	  nendd = -1; 
      }
      
      
  MPI_Allreduce(&stars_spawned, &tot_spawned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(tot_spawned > 0) 
  { 
      All.TotNumPart += tot_spawned;
      NumPart += stars_spawned;      
  }  
  double rprnt = gsl_rng_uniform(random_generator); /* to avoid too much output in log file */
  if (rprnt < 0.05 && ThisTask == 0)
//  if (ThisTask == 0)
  {
//	  printf("bits = %d, Made iteration %d, total time remaining %g \n",
//		 bits, N_iter, dt_tot_all );

      printf("N_active = %d,  N_iter = %d, stars_spawned %d, NumPart %d ! \n", 
	     N_active, N_iter,stars_spawned, NumPart);
//      printf("Destroyed %d and Emitted %d photons, N_active = %d,  N_iter = %d ! \n", 
//	     count_deleted_photon_total, tot_spawned, N_active_tot, N_iter);
      fflush(stdout);
  }

  stars_spawned = 0;
  tot_spawned = 0;
/* update density at photon's location for use in the next loop*/
      virtual_density();
      
      
  }
      


  
#ifdef RAD_ACCEL /* SHC */
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
  {
      if(P[i].Type == 0 )
      {
	  if(All.ComovingIntegrationOn)
	      dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / time_hubble_a;
	  else
	      dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;
	  
#ifndef FB_MOMENTUM
	  for (k = 0 ; k < 3; k++)
	      SphP[i].ra.RadAccel[k] = 0.;
#endif

  	  if (dt != 0. && P[i].Mass > 0.)
  	  {
  	      for (k = 0 ; k < 3; k++)
  	      {
  		  SphP[i].ra.RadAccel[k] = SphP[i].i.Injected_BH_Momentum[k]/ dt / P[i].Mass;
  	      }
  	      ax =  SphP[i].ra.RadAccel[0];
  	      ay =  SphP[i].ra.RadAccel[1];
  	      az =  SphP[i].ra.RadAccel[2];
  	      ac = sqrt(ax * ax + ay * ay + az * az);
  	      acmax = 1.e3/dt;
  	      if (ac >= acmax)
  	      {
/*  		  printf("  RadAccel = (x,y,z) = (%g, %g, %g), dt = %g, mass = %g \n",  */
/*  			 SphP[i].ra.RadAccel[0], SphP[i].ra.RadAccel[1], SphP[i].ra.RadAccel[2], dt, P[i].Mass);  */
/*  		  fflush(stdout);  */
  		  for (k = 0 ; k < 3; k++)
  		  {
  		      SphP[i].ra.RadAccel[k] *= (acmax/ac);
  		  }
  	      }
  	  }
	  for (k = 0 ; k < 3; k++)
	  {
	      if (P[i].Mass > 0) P[i].Vel[k] = (P[i].Vel[k] * P[i].Mass + SphP[i].i.Injected_BH_Momentum[k])/P[i].Mass;
	      SphP[i].i.Injected_BH_Momentum[k] = 0.;
	  }
      }
  }
  
#endif /* RAD_ACCEL */

#ifdef SMART_HEATING

#ifdef VIRTUAL_HEATING
/*  THIS IS ENTERED WHEN SPH ENERGY BALANCE NEEDS TO BE SOLVED
 */
  double umin, umax;

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
  {
      dtime =  (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;
      if(P[i].Type == 0 && dtime > 0)
      {
		  if(P[i].Mass == 0)
		  {
		      SphP[i].Injected_VIRTUAL_Energy = 0;
		      continue;
		  }
/* this is for SPH particles to be removed by accretion on black holes or leaving
   the computational    domain*/
		  else
		  {
		      uold = DMAX(All.MinEgySpec, SphP[i].Entropy/
				  GAMMA_MINUS1 * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1));
		      updv = DMAX(All.MinEgySpec,
				  (SphP[i].Entropy + SphP[i].e.DtEntropy * dtime) /
				  GAMMA_MINUS1 * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1));
		      dudt_pdv = (updv - uold)/dtime;

/* RadEnergyDensity is radiation energy density at the SPH particle in physical units */
/* It is calculated using the amount of radiation heat absorbed by the SPH
 * particlew during time dtime. Using the known cross section it is easy to
 * reconstruct the radiation energy density then
 */

		      RadEnergyDensity = SphP[i].Injected_VIRTUAL_Energy/P[i].Mass/dtime * EradFactor;
		      RadHeating = SphP[i].Injected_VIRTUAL_Energy/P[i].Mass;
#ifdef DECREASE_POISSON_NOISE 
		      SphP[i].Injected_VIRTUAL_Energy_Old = (3*SphP[i].Injected_VIRTUAL_Energy_Old + SphP[i].Injected_VIRTUAL_Energy)/4.;
		      RadEnergyDensity = SphP[i].Injected_VIRTUAL_Energy_Old/P[i].Mass/dtime * EradFactor;
		      RadHeating = SphP[i].Injected_VIRTUAL_Energy_Old/P[i].Mass;
#endif

		      



#ifdef BACKGROUND_ILLUMINATION_TEMP
		      double T_bgr = All.EqTemp ;
//		      double u_equilibrium_prime = pow(RadEnergyDensity/RAD_CONST, 0.25)/u_to_temp_fac;
//		      RadEnergyDensity += pow(T_bgr, 4) * RAD_CONST;
#endif


		      double k_factor = 1./EradFactor;
/* using this k-factor ensures the correct units for the heating rate below in t_cool_estimate */
/* k_factor = THOMPSON * All.VirtualCrosSection/PROTONMASS * C *
   All.UnitTime_in_s/All.UnitVelocity_in_cm_per_s/All.UnitVelocity_in_cm_per_s;
*/

#ifdef REAL_EOS
/* COOLING by iterayion shc -- 01.07.09 */

/*
                      double shc_temp=4.;
                      double shc_tmp_den;
                      if (ThisTask == 0) shc_tmp_den=1.e-19;
                      if (ThisTask == 1) shc_tmp_den=1.e-13;
                      if (ThisTask == 2) shc_tmp_den=1.e-10;
                      if (ThisTask == 3) shc_tmp_den=1.e-7;
                      while (shc_temp < 2000.) {
                         fprintf (diag_shc,"%g %g %g\n",shc_temp,ene_per_gram(1,shc_tmp_den,shc_temp),SphP[1].Mu);
                         shc_temp++;
                      }       
                      exit (99999);
*/

                      double shcha_rho = SphP[i].d.Density * All.UnitDensity_in_cgs;
		      
//                      double shcha_T_left = (4. > uold*u_to_temp_fac/10.) ? 4.:uold*u_to_temp_fac/10.;
//                     double shcha_T_right = (2000. < uold*u_to_temp_fac*10.) ? 2000.:uold*u_to_temp_fac*10.;
                      double shcha_T_left =  4.;
		      double shcha_T_right = 2000.;
                      double shcha_T_center = (shcha_T_left + shcha_T_right) / 2.;

                      double shcha_u_left = ene_per_gram(i,shcha_rho,shcha_T_left)/All.UnitEnergy_in_cgs*All.UnitMass_in_g;
                      double shcha_u_right = ene_per_gram(i,shcha_rho,shcha_T_right)/All.UnitEnergy_in_cgs*All.UnitMass_in_g;
                      double shcha_u_center = ene_per_gram(i,shcha_rho,shcha_T_center)/All.UnitEnergy_in_cgs*All.UnitMass_in_g;

                      double shcha_F_left = shcha_u_left
                           - (uold+dtime*(dudt_pdv+k_factor*(RadEnergyDensity-RAD_CONST*pow(shcha_T_left,4))));
                      double shcha_F_right = shcha_u_right
                           - (uold+dtime*(dudt_pdv+k_factor*(RadEnergyDensity-RAD_CONST*pow(shcha_T_right,4))));
                      double shcha_F_center = shcha_u_center
                           - (uold+dtime*(dudt_pdv+k_factor*(RadEnergyDensity-RAD_CONST*pow(shcha_T_center,4))));

                      int shcha_icount = 0;
//                     while (fabs(shcha_F_center)/uold > 1.e-2 && shcha_icount < 50)
                      while (fabs(shcha_F_center)/uold > 1.e-2)
                      {
                           if (shcha_F_left * shcha_F_center < 0.)
                           {   
                               shcha_T_right = shcha_T_center;
                               shcha_F_right = shcha_F_center;
                           }
                           else
                           {
                               shcha_T_left = shcha_T_center;
                               shcha_F_left = shcha_F_center;
                           }
                           shcha_T_center = (shcha_T_left + shcha_T_right) / 2.;
                           shcha_u_center = ene_per_gram(i,shcha_rho,shcha_T_center)/All.UnitEnergy_in_cgs*All.UnitMass_in_g;
                           shcha_F_center = shcha_u_center
                         - (uold+dtime*(dudt_pdv+k_factor*(RadEnergyDensity-RAD_CONST*pow(shcha_T_center,4))));

                           shcha_icount++;
                      }
                      if (shcha_icount >= 50) {
                         unew = uold;
                      } else {
                         unew = shcha_u_center;
                         SphP[i].Temperature = shcha_T_center;
                      }

                      fprintf (diag_shc,"%d %g %g %g %g %g %g\n",i,unew,unew*All.UnitEnergy_in_cgs/All.UnitMass_in_g,shcha_rho,SphP[i].Temperature,SphP[i].Mu,SphP[i].Kappa);
                      fflush (diag_shc);
/*
*/

#else




/* EXPLICIT COOLING sn -- 26.05.09 */

/* 		      unew = uold + dtime * (dudt_pdv  + k_factor *  */
/* 			  (RadEnergyDensity - RAD_CONST*pow(uold * u_to_temp_fac ,4))); */


#ifdef DUST_OPACITY_FIT     
	      /* rad energy density of the central star at radius r */
		      rhs = uold + dtime * (dudt_pdv  + k_factor * pow(T_bgr, 5)/All.EqTemp * RAD_CONST +
					    k_factor * RadEnergyDensity * (200./All.EqTemp));
		      u_eq = pow((pow(T_bgr, 4) + 20.*RadEnergyDensity/RAD_CONST), 0.25)/u_to_temp_fac;
#else
//		      rhs = uold + dtime * (dudt_pdv  + k_factor * (pow(T_bgr, 4) * RAD_CONST + RadEnergyDensity));
//		      u_eq = pow((pow(T_bgr, 4) + RadEnergyDensity/RAD_CONST), 0.25)/u_to_temp_fac;

		      rhs = uold + dtime * (dudt_pdv  + k_factor * pow(T_bgr, 4) * RAD_CONST) + RadHeating;
		      u_eq = pow((pow(T_bgr, 4) + RadEnergyDensity/RAD_CONST), 0.25)/u_to_temp_fac;

#endif


#ifndef DUST_OPACITY_FIT 
//	      rhs = uold + dtime * (dudt_pdv  + k_factor * RadEnergyDensity);
	      if (SphP[i].d.Density*All.UnitDensity_in_cgs > 1.e-13) k_factor *= pow(10.,5.2)*pow(SphP[i].d.Density*All.UnitDensity_in_cgs,0.4);
	      const4 = k_factor * RAD_CONST*pow(u_to_temp_fac ,4) * dtime;
	      u1 = rhs;
	      u2 = pow(rhs/const4, 0.25);
	      if (u2 < umin) umin = u2;
	      umin = 0.03 * umin + All.MinGasTemp/u_to_temp_fac ;
	      umax = u1;
	      if (u2 > umax) umax = 2.*u2 + 200.;
	      double u1save = umin, u2save = umax;

	      unew = sqrt(umin * umax);
	      
	      lhs = unew + const4 * pow(unew, 4);
	      diff = fabs(lhs - rhs);
	      int count = 0;
	      while (diff > 0.001 * rhs)
	      {
		  if (lhs <= rhs)
		  {
		      umin = unew;
		  }
		  else
		  {
		      umax = unew;
		  }
		  unew = sqrt(umin * umax);
		  lhs = unew + const4 * pow(unew, 4);
		  diff = fabs(lhs - rhs);
		  count++;
		  if (count > 300)
		  {
		      printf("\n Cannot find solution for T in fb_particles.c ... \n");
		      printf("\n temp = %g,  t1 = %g, t2 = %g ... \n", u_to_temp_fac * unew, u_to_temp_fac * u1save, u_to_temp_fac * u2save);
		      exit(1);
		  }
	      }
#else
	      const5 = k_factor * RAD_CONST*pow(u_to_temp_fac ,5) * dtime/All.EqTemp;
	      u1 = rhs;
	      u2 = pow(rhs/const5, 0.2);
	      umin = u1;
	      if (u2 < umin) umin = u2;
	      umin = 0.1 * umin + All.MinGasTemp/u_to_temp_fac;
	      umax = u1;
	      if (u2 > umax) umax = 2.*u2 + 200.;
	      double u1save = umin, u2save = umax;
	      
	      unew = sqrt(umin * umax);
	      
	      lhs = unew + const5 * pow(unew, 5);
	      diff = fabs(lhs - rhs);
	      int count = 0;
	      while (diff > 0.0001 * rhs)
	      {
		  if (lhs < rhs)
		  {
		      umin = unew;
		  }
		  else
		  {
		      umax = unew;
		  }
		  unew = sqrt(umin * umax);
		  lhs = unew + const5 * pow(unew, 5);
		  diff = fabs(lhs - rhs);
		  count++;
		  if (count > 300)
		  {
		      printf("\n Cannot find solution for T in fb_particles.c ... \n");
		      printf("\n temp = %g,  t1 = %g, t2 = %g ... \n", u_to_temp_fac * unew, u_to_temp_fac * u1save, u_to_temp_fac * u2save);
		      exit(1);
		  }
	      } 
#endif

/* END of NEW energy equilibrium approach */


#endif



/* This definition of the energy to be emitted by the SPH particle makes sure that the energy is conserved
 */
/* #ifdef BACKGROUND_ILLUMINATION_TEMP */
/* //		      SphP[i].Emitted_VIRTUAL_Energy = (uold - unew_prime)* P[i].Mass + SphP[i].Injected_VIRTUAL_Energy * */

#ifndef RADIATIVE_EQUILIBRIUM

	      double tnew = u_to_temp_fac * unew;
	      double Reflect = 1.;
	      if (RadEnergyDensity > 0) 
	      {
		  Reflect = pow(tnew, 4) * RAD_CONST/RadEnergyDensity;
	      }
	      if (Reflect >= 1.) 
	      {
		  Reflect = 1.;
#ifdef DECREASE_POISSON_NOISE
		  SphP[i].Emitted_VIRTUAL_Energy = pow(tnew, 4) * RAD_CONST * k_factor*dtime*P[i].Mass - 
		      SphP[i].Injected_VIRTUAL_Energy_Old;
#else
		  SphP[i].Emitted_VIRTUAL_Energy = pow(tnew, 4) * RAD_CONST * k_factor*dtime*P[i].Mass - 
		      SphP[i].Injected_VIRTUAL_Energy;
#endif
	      }
	      else
	      {
		  SphP[i].Emitted_VIRTUAL_Energy = 0.;
	      }
	      SphP[i].Returned_E_Fraction = Reflect;

#else /* RADIATIVE_EQUILIBRIUM */

		  SphP[i].Emitted_VIRTUAL_Energy = (dudt_pdv  + pow(T_bgr, 4) * RAD_CONST * k_factor) * dtime * P[i].Mass;
		  SphP[i].Returned_E_Fraction = 1.;
#endif

	      SphP[i].e.DtEntropy = (unew * GAMMA_MINUS1/
                                     pow(SphP[i].d.Density * a3inv,
                                         GAMMA_MINUS1) - SphP[i].Entropy) / dtime;
/* this is how radiative heating cooling is passed to the rest of the code */
/* D Entropy/dt is calculated */

              SphP[i].Injected_VIRTUAL_Energy = 0.;
		  }
	      }
	  }


#ifdef SPH_EMIT


/* ------------------------------------------------------- //
   Creation of new virtual/photon particles by SPH; 
   Needed only when energy balance is solved ----------//
// ------------------------------------------------------- */	    

  double photon_energy_new = 0., energy_to_be_emitted = 0;
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
  {
      if(P[i].Type != 0) continue;

      if(SphP[i].Emitted_VIRTUAL_Energy > 0) 
      {
	  energy_to_be_emitted += SphP[i].Emitted_VIRTUAL_Energy;
/* decide how many new virtual particles are needed, and how energetic they should be */

	  dt_create = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;
/* Decide the energy and momentum of the new photons to be emitted by the SPH
 * particles. If chosen too low, it will lead to too many photon packets, and
 * the code will exceed the allowed maximum number of particles.  The choise
 * below ensures that the number of photon packets emitted is independent of
 * the amount of energy needed to be emitted by this SPH particle in the
 * timestep dt
 */
//	  energy_gamma = 5./All.PartAllocFactor * SphP[i].Emitted_VIRTUAL_Energy/dt_create * All.VirtualTime + 


	  energy_gamma = All.VirtualMomentum * C/All.UnitVelocity_in_cm_per_s;
/* this is the number of new photons emitted by the given SPH particle */	  
	  new_virt_exact = SphP[i].Emitted_VIRTUAL_Energy/energy_gamma;
	  
/* new_virt is the integer number of photons to be emitted */

	  new_virt = (int) new_virt_exact;
//	  if (new_virt_exact > 1.e-2) new_virt = (int) new_virt_exact + 1;
//	  new_virt = (int) new_virt_exact + 1;

	  for(bits = 0; bits < new_virt+1; bits++)
	  {
	      
/* here we spawn new photon particles. Notice that for bits <= new_virt,
 * photons emitted without questions, but for bits > new_virt, we emit a photon
 * with probability given by the excess of new_virt_exact over new_virt */


	      correction_factor = 1.;//new_virt_exact/new_virt;
//	      correction_factor = new_virt_exact/new_virt;

	      if (bits == new_virt)
	      {
		  wemit = gsl_rng_uniform(random_generator);
//		  wemit = get_random_number(P[i].ID + bits + All.NumCurrentTiStep);
		  if (wemit > (new_virt_exact - new_virt)) continue;
	      }

	      if(NumPart + stars_spawned >= All.MaxPart)
	      {
		  printf
		      ("On Task=%d with NumPart=%d we try to spawn %d particles. Sorry, no space left...(All.MaxPart=%d)\n",
		       ThisTask, NumPart, stars_spawned, All.MaxPart);
		  printf("  bits %d, new_virt_exact %g \n", bits, new_virt_exact);
		  updv = DMAX(All.MinEgySpec,
			      (SphP[i].Entropy + SphP[i].e.DtEntropy * dt_create) /
			      GAMMA_MINUS1 * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1));
		  printf("  gas T = %g, u %g, e_gamma %g \n", updv*u_to_temp_fac, updv*P[i].Mass, energy_gamma);
		  fflush(stdout);
		  endrun(8888);
	      }
	      
	      P[NumPart + stars_spawned] = P[i];
	      P[NumPart + stars_spawned].Type = 3; /* change the type of
						    * particle into photon
						    * particle */
#ifdef STELLARAGE
/* StellarAge is time when the particle was created  */
	      P[NumPart + stars_spawned].StellarAge = All.Time; /* present time */
#endif
	      
/*	      theta = acos(2. * get_random_number(P[i].ID + 2 + bits + All.NumCurrentTiStep) - 1.);
	      phi = 2 * M_PI * (get_random_number(P[i].ID + 3 + bits + All.NumCurrentTiStep) + All.Time);*/

	      theta = acos(2. * gsl_rng_uniform(random_generator)- 1.);
	      phi = 2 * M_PI * (gsl_rng_uniform(random_generator) + All.Time);	      

	      dir[0] = sin(theta) * cos(phi);
	      dir[1] = sin(theta) * sin(phi);
	      dir[2] = cos(theta);
	      
	      PPP[NumPart + stars_spawned].Hsml = 2.*P[NumPart + stars_spawned].Hsml;

#ifdef NON_GRAY_RAD_TRANSFER
	      uold = DMAX(All.MinEgySpec, SphP[i].Entropy/
			  GAMMA_MINUS1 * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1));
/* photons from the sink are born with high temperature, resulting in high opacity */
	      P[NumPart + stars_spawned].VirtualTemperature = uold *  u_to_temp_fac;
/* opacity is approximately kappa ~ kappa_0 (T/T0) */
	      P[NumPart + stars_spawned].Kappa = P[NumPart + stars_spawned].VirtualTemperature/All.EqTemp;
#endif

	      for(j = 0; j < 3; j++)
	      {
/* particles are kicked in real and velocity space in direction dir[j] */
		  P[NumPart + stars_spawned].Vel[j] = C/All.UnitVelocity_in_cm_per_s * dir[j]/All.FeedBackVelocity;
		  P[NumPart + stars_spawned].Pos[j] += (0.0 + 0.1*gsl_rng_uniform(random_generator))
		      * P[i].Hsml * dir[j];
/* This needs to be re-examined to make sure photons can escape from the initial position */
	      }
		  
#ifdef VIRTUAL_FLY_THORUGH_EMPTY_SPACE
	      P[NumPart + stars_spawned].b1.BH_Density = 0.;
#endif	      
	      NextActiveParticle[NumPart + stars_spawned] = FirstActiveParticle;
	      FirstActiveParticle = NumPart + stars_spawned;
	      NumForceUpdate++;
	      TimeBinCount[P[NumPart + stars_spawned].TimeBin]++;
	      
	      PrevInTimeBin[NumPart + stars_spawned] = i;
	      NextInTimeBin[NumPart + stars_spawned] = NextInTimeBin[i];
	      if(NextInTimeBin[i] >= 0)
		  PrevInTimeBin[NextInTimeBin[i]] = NumPart + stars_spawned;
	      NextInTimeBin[i] = NumPart + stars_spawned;
	      if(LastInTimeBin[P[i].TimeBin] == i)
		  LastInTimeBin[P[i].TimeBin] = NumPart + stars_spawned;
	      
	      P[NumPart + stars_spawned].ID += NumPart + stars_spawned;
	      P[NumPart + stars_spawned].NewDensity = 0.;
	      
//		    P[NumPart + stars_spawned].Mass = All.VirtualMass;
/* PhotonMomentum us the change in the OldPhotonMomentum in a single call to the SPH/photon interaction routine */
	      
	      
	      P[NumPart + stars_spawned].OldPhotonMomentum = energy_gamma *correction_factor/C * All.UnitVelocity_in_cm_per_s;
	      P[NumPart + stars_spawned].BirthPhotonEnergy = energy_gamma * correction_factor;
	      
	      P[NumPart + stars_spawned].DeltaPhotonMomentum = 0.;
	      P[NumPart + stars_spawned].DeltaPhotonEnergy = 0.;
	      P[NumPart + stars_spawned].DeltaHeat = 0.;
	      P[NumPart + stars_spawned].WindOrPhoton = 1; /* This marks it as a photon */
/* initialise photon momentum */
	      P[NumPart + stars_spawned].Mass = All.VirtualMass * P[NumPart + stars_spawned].OldPhotonMomentum;
//		    PPP[NumPart + stars_spawned].Hsml = All.InnerBoundary * P[i].Mass/2.;
//		    sum_mass_stars += P[NumPart + stars_spawned].Mass;
#ifndef EXCLUDE_VIRTUAL_FROM_TREE
	      force_add_star_to_tree(i, NumPart + stars_spawned);
#endif
	      stars_spawned++;
	      photon_energy_new += energy_gamma * correction_factor;
	  }
      }
  }
  tot_spawned = 0.;
  MPI_Allreduce(&stars_spawned, &tot_spawned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(tot_spawned > 0) 
  { 
      All.TotNumPart += tot_spawned;
      NumPart += stars_spawned;      
  }  
  if (ThisTask == 0 && tot_spawned > 0)  
  {
      printf
	  ("We spawn %d photon packets; energy to be emitted %g, new photon en %g \n",tot_spawned, energy_to_be_emitted, photon_energy_new);
      fflush(stdout);
  }

#endif /* SPH_EMIT */

#endif /* VIRTUAL_HEATING */

/* ELIMINATION OF PHOTONS */

/* check if feedback/virtual particle has left an interesting sphere where I
 * track its motion. Or if it is "too old". In either case, remove it */


/*   int count_deleted_photon = 0, count_deleted_photon_total=0; */
/* //  double deleted_momentum = 0, total_deleted_momentum = 0; */
/*   for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) */
/*   { */
/*       if (P[i].Type == 3 && P[i].DtJump >= 0)  */
/*       { */
/* 	  r = sqrt(P[i].Pos[0]*P[i].Pos[0] + P[i].Pos[1]*P[i].Pos[1] + P[i].Pos[2]*P[i].Pos[2]); */
/* 	  if (r > All.OuterBoundary || All.Time-P[i].StellarAge > All.VirtualTime)  */
/* 	  { */
/* 	      P[i].Mass = 0.; */
/* 	      count_deleted_photon++; */
/* 	      deleted_momentum += P[i].OldPhotonMomentum; */
/* 	  }  */
/* 	  else  if (P[i].OldPhotonMomentum < 0.001 * P[i].BirthPhotonEnergy/C * All.UnitVelocity_in_cm_per_s)  */
/* 	  { */
/* 	      P[i].Mass = 0.; */
/* 	      count_deleted_photon++; */
/* //            deleted_momentum += P[i].OldPhotonMomentum; */
/* 	  } */
/*       } */
/*   } */

/*   fflush(diag_shc_photon); */
  
/*   MPI_Allreduce(&count_deleted_photon, &count_deleted_photon_total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); */
/*   MPI_Allreduce(&deleted_momentum, &total_deleted_momentum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); */
/*   MPI_Allreduce(&added_energy_by_bh, &total_added_energy_by_bh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); */

/*   if(ThisTask == 0) */
/*   { */
/* /\* convert it to energy in cgs *\/ */
/*       total_deleted_momentum *= C/All.UnitVelocity_in_cm_per_s * All.UnitEnergy_in_cgs ; */

/*       fprintf (diag_shc_photon_total,"%g %i %g\n",All.Time, count_deleted_photon_total, total_deleted_momentum); */
/*       fflush(diag_shc_photon_total); */
/* /\* Also, print the number of new and deleted photons to the screen *\/ */
/*       printf("FB PARTICLES: spawned %d and destroyed %d photon particles \n", tot_spawned, count_deleted_photon_total); */
/*       fflush(stdout); */

/*       /\* SHC BEGIN *\/ */
/* //    SysState.EnergyRadDeleted += total_deleted_momentum * C/All.UnitVelocity_in_cm_per_s; */
/*       /\* SHC END *\/ */




#endif /* SMART_HEATING */
  
//#ifndef NO_NEW_PHOTONS
/* ELIMINATION OF PHOTONS */

/* check if feedback/virtual particle has left an interesting sphere where I
 * track its motion. Or if it is "too old". In either case, remove it */


  int count_deleted_photon = 0, count_deleted_photon_total=0;
  double deleted_momentum = 0, total_deleted_momentum = 0, total_photon_energy_new = 0;
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
  {
      if (P[i].Type == 3 && P[i].Mass != 0.) 
      {
	  r = sqrt(P[i].Pos[0]*P[i].Pos[0] + P[i].Pos[1]*P[i].Pos[1] + P[i].Pos[2]*P[i].Pos[2]);
//	  fprintf (photonSN,"%d, %g, %g, %g, %g, %g, %g \n",i , All.Time, P[i].dt1, P[i].dt2, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);

	  if (r > All.OuterBoundary || All.Time-P[i].StellarAge > All.VirtualTime) 
	  {
	      P[i].Mass = 0.;
	      count_deleted_photon++;
	      deleted_momentum += P[i].OldPhotonMomentum;
	      P[i].OldPhotonMomentum = 0;
	  } 
/*	  else  if (P[i].OldPhotonMomentum < 0.00001 * P[i].BirthPhotonEnergy/C * All.UnitVelocity_in_cm_per_s) 
	  {
	      P[i].Mass = 0.;
	      count_deleted_photon++;
//            deleted_momentum += P[i].OldPhotonMomentum;
	  }
*/
      }
  }
//  double background_heating = energy_to_be_emitted * All.UnitEnergy_in_cgs;
/*  double energy_deleted = photon_energy_new;
*/


//  fflush (photonSN);
   MPI_Allreduce(&count_deleted_photon, &count_deleted_photon_total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
   MPI_Allreduce(&deleted_momentum, &total_deleted_momentum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
   MPI_Allreduce(&photon_energy_new, &total_photon_energy_new, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 

   double energy_deleted = total_deleted_momentum * C/All.UnitVelocity_in_cm_per_s;
   double background_heating = pow(All.EqTemp, 4) * RAD_CONST * (All.TimeStep * All.UnitTime_in_s) * (0.2 * All.UnitMass_in_g) 
       * THOMPSON * All.VirtualCrosSection/PROTONMASS * C;
//   background_heating += CurrentLuminosity * All.UnitEnergy_in_cgs * All.TimeStep;


   SysState.EnergyRadDeleted += energy_deleted; 
   SysState.EnergyRadAdded += background_heating/All.UnitEnergy_in_cgs*All.UnitTime_in_s ; 


//   printf(" Energy eliminated %g, heating %g \n", energy_deleted, background_heating);

//#endif /* NO_NEW_PHOTONS */


  CPU_Step[CPU_RHD] += measure_time();

}


#endif /* closing FB_PARTCLES */
  
