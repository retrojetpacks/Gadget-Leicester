#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"
#include "forcetree.h"



#if defined(VIRTUAL)  && defined(FB_MOMENTUM)

/* This routine creates momentum feedback particles, propagates them and
 * transfer the momentum to SPH particles. No photons here -- they are in
 * fb_particles.c
 */

void momentum_feedback(void)
/* cooling routine when star formation is enabled */
{
  unsigned int bits;
  int i, j, k, prev, bin, flag, stars_spawned, tot_spawned, new_wind, N_iter, 
      N_active, N_active_tot;
  struct state_of_system sys;
  double ascale = 1, hubble_a = 0, a3inv, wemit, correction_wind;
  double time_hubble_a, dir[3], theta, phi, norm, dt, dt_create, rseed, dtau, r, new_wind_exact;
  double CurrentLuminosity, acmax, rho, dtau_max,
      dtime, dmax1, dmax2, energy_gamma, energy_wind_min, ax, ay, az, ac,
      dr_jump;
  double desired_dr, dt_bh, dt_total, r2, dx, dy, dz, dt_tot_all;


/* SN: Multiplying this by particle mass in code units gives Eddington limit
 * in code units, defined using Thompson opacity */
  double L_Edd_dimensionless = (4 * M_PI * GRAVITY * C * PROTONMASS/ THOMPSON)/pow(All.UnitVelocity_in_cm_per_s, 2) *
      All.UnitTime_in_s;

  energy_wind_min = All.VirtualMomentum * C/All.UnitVelocity_in_cm_per_s;


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
	  new_wind_exact = All.VirtualFeedBack * L_Edd_dimensionless * P[i].Mass * dt_create/energy_wind_min;

#ifdef ACCRETION_AT_EDDINGTON_RATE
	  new_wind_exact = All.VirtualFeedBack * L_Edd_dimensionless * P[i].BH_Mass * dt_create/energy_wind_min;
	  P[i].Mass = P[i].BH_Mass; /* CAUTION: these two are assumed equal in this case */
#endif
#else
/* SN -- Feedback appropriate for stellar radiation sources */
/* The first bit is the accretion luminosity, assuming GM/R_core = 0.01 Msun/ 1Rsun */
/* The second bit is the main sequence luminosity of the star */
//	  CurrentLuminosity = (1.99e31/All.UnitMass_in_g)/(1.e11/All.UnitLength_in_cm) * P[i].BH_Mdot;
//              + SOLAR_LUM/All.UnitEnergy_in_cgs * All.UnitTime_in_s * pow(P[i].Mass * All.UnitMass_in_g/SOLAR_MASS, 3.3);



	  CurrentLuminosity = All.VirtualFeedBack * 0.1 * P[i].BH_Mdot *  pow(C/All.UnitVelocity_in_cm_per_s, 2);

/*  Now limit that at Eddington luminosity */
	  if (CurrentLuminosity > L_Edd_dimensionless * P[i].Mass) CurrentLuminosity = L_Edd_dimensionless * P[i].Mass;

	  new_wind_exact = All.VirtualFeedBack * CurrentLuminosity * dt_create/energy_wind_min;

#endif

	  printf("Accretion Rate %g, new_wind_exact %g, Bh mass %g, dt %g \n", 
		 P[i].BH_Mdot, new_wind_exact, P[i].Mass, dt);
	  fflush(stdout);  

/* Deciding how many NEW WIND PARTICLES excatly needed */
#ifdef STARTUP_LUMINOSITY
	  if (new_wind_exact > 1.e-4) 
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
	  
	  
 	  if (new_wind > 0)
 	  {
/*  		printf("On Task=%d Need to spawn %d virtual particles, new_virt_exact is equal to %g \n", ThisTask, new_virt, new_virt_exact);
  		fflush(stdout);
*/
 	      flag=0;
 	  }
	  
	  if(flag == 0)		/* active virtual particle creation */
	  {
	      for(bits = 0; bits < new_wind; bits++)
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
		  theta = acos(2. * get_random_number(P[i].ID + 2 + bits + All.NumCurrentTiStep) - 1.);
#else
		  theta = acos(0.5 * get_random_number(P[i].ID + 2 + bits + All.NumCurrentTiStep) + 0.5);
#endif
		  phi = 2 * M_PI * (get_random_number(P[i].ID + 3 + bits + All.NumCurrentTiStep) + All.Time);
		  
		  dir[0] = sin(theta) * cos(phi);
		  dir[1] = sin(theta) * sin(phi);
		  dir[2] = cos(theta);

		  
		  P[NumPart + stars_spawned].Hsml += All.SofteningBndryMaxPhys;
#ifdef ACCRETION_RADIUS
		  P[NumPart + stars_spawned].Hsml += All.InnerBoundary;
#endif
		  PPP[NumPart + stars_spawned].Hsml = P[NumPart + stars_spawned].Hsml;

		  for(j = 0; j < 3; j++)
		  {
/* particles are kicked in real and velocity space in direction dir[j] */
		      
//			  P[NumPart + stars_spawned].Vel[j] = All.FeedBackVelocity * dir[j];
		      P[NumPart + stars_spawned].Vel[j] = C/All.UnitVelocity_in_cm_per_s * dir[j]/All.FeedBackVelocity;
		      
		      P[NumPart + stars_spawned].Pos[j] += 0.*(0.25 + 0.25*gsl_rng_uniform(random_generator))
			  * P[i].Hsml * dir[j];
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
		  
/* DeltaPhotonMomentum us the change in the OldPhotonMomentum in a single SPH/photon interaction */
		  P[NumPart + stars_spawned].DeltaPhotonMomentum = 0.;
#if defined(VIRTUAL_HEATING) 
		  P[NumPart + stars_spawned].DeltaPhotonEnergy = 0.;
#endif
//bits < new_virt + new_wind;
		  P[NumPart + stars_spawned].WindOrPhoton = 0; /* This marks it as a Wind particle */
		  P[NumPart + stars_spawned].OldPhotonMomentum = energy_wind_min/C* All.UnitVelocity_in_cm_per_s * correction_wind;
		  P[NumPart + stars_spawned].BirthPhotonEnergy = energy_wind_min * correction_wind;

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

  


/*  This NextPhoton list */

  FirstPhoton = -1;  
  prev = -1;
//  for(i = 0; i < NumPart; i++) 
  for(i = FirstActiveParticle ; i >= 0; i = NextActiveParticle[i])
  { 
/* Note: photons that were just created are inactive until timestep.c assigns them a time step */
      if(P[i].Type == 3 &&  P[i].StellarAge < All.Time  && P[i].WindOrPhoton == 0)
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
  for(i = FirstActivePhoton; i >= 0; i = NextActivePhoton[i])
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
  if (dt_tot_all <= 1.e-15 * All.VirtualTime * (All.MaxPart*NTask)) nendd = -1; 
/* don't even enter the loop in this case */

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
/* Desired spatial jump in one step for the photon. Limited by spatial
 * resolution requirements */
	  desired_dr = 0.2 * P[i].Hsml; 
#ifdef VIRTUAL_FLY_THORUGH_EMPTY_SPACE
/* If there are no SPH particles containing the photon in its Hsml, but there
 * are SPH particles in the search radius of the photon, be careful. I the
 * photon desired_dr is too large it propagates too far in one step, skipping
 * lots of sph particles that it should have interacted with if it were
 * propagated by smaller desired_dr. BH_Density is actually the minimum
 * distance to an SPH particle within the Hsml, so this provides a convenient
 * way of allowing small enough steps */
	  if (desired_dr > 0.5 * P[i].b1.BH_Density + 2. * All.SofteningGas)
	  {
	      desired_dr = 0.5 * P[i].b1.BH_Density + 2. * All.SofteningGas;
	  }
#endif

/* if desired spatial jump takes longer than the requested time-jump to
 * inactive status, reduce dr accordingly */
	  if (desired_dr > C*dt_bh/All.UnitVelocity_in_cm_per_s/All.FeedBackVelocity) 
	      desired_dr = C*dt_bh/All.UnitVelocity_in_cm_per_s/All.FeedBackVelocity;

	  P[i].dt1 = desired_dr * All.UnitVelocity_in_cm_per_s/C * All.FeedBackVelocity;

/* Let the wind particle be absorbed in several steps rather than one. This
 * should reduce statistical noise of the approach */

	  dtau_max = 10. * desired_dr/P[i].Hsml;
	  dtau = All.VirtualCrosSection * THOMPSON/PROTONMASS *
              P[i].NewDensity * desired_dr * All.UnitDensity_in_cgs *
              All.UnitLength_in_cm;
	  if (dtau > dtau_max) dtau = dtau_max; 

	  P[i].DeltaPhotonMomentum = P[i].OldPhotonMomentum * (1. - exp(-dtau));
	  P[i].OldPhotonMomentum -= P[i].DeltaPhotonMomentum;

/* P[i].Mass is equal to P[i].OldPhotonMomentum save for the constant in
 * front. So it can be used to track photon momentum in snapshot files.  */
	  P[i].Mass = All.VirtualMass * P[i].OldPhotonMomentum;

/* Now make the jump */
	  for (k = 0 ; k < 3; k++)
	      P[i].Pos[k] += P[i].dt1 * P[i].Vel[k];
				  
/* reduce the amount of time remaining till the photon becomes inactive */
	  dt_bh -= P[i].dt1;
	  P[i].DtJump -= P[i].dt1;
	  dt_total -= P[i].dt1;	      
/* This is how long the photon has been active */
	  P[i].dt2 += P[i].dt1;
	  
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



	 if (P[i].OldPhotonMomentum < 0.01 * P[i].BirthPhotonEnergy/C * All.UnitVelocity_in_cm_per_s) 
	  {
/* Hide this particle then, by stating that it has already propagated as long as needed to become inactive */
	      dt_total -= P[i].DtJump;	      
	      P[i].dt2 += P[i].DtJump;
	      P[i].DtJump = 0.e0;
	      P[i].dt1 = 0.;
	  } 


      }
      

/* update density at photon's location for use in the next loop*/
      virtual_density();


/* Now transfer the absorbed momentum to SPH particles*/
/* ------------------------------- */
      virtual_spread_momentum();
/* ------------------------------- */


// This marks the photons that completed their time jumps as inactive 
       FirstPhoton = -1; 
       prev = -1; 
       N_active = N_active_tot = 0.;
       for(i = FirstActiveParticle ; i >= 0; i = NextActiveParticle[i]) 
       { 
 	  if(P[i].Type == 3 && P[i].DtJump > 1.e-40 && P[i].WindOrPhoton == 0) 
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
      MPI_Allreduce(&dt_total, &dt_tot_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&N_active, &N_active_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      
      if (N_active_tot == 0) 
      {
	  nendd = -1; 
	  if(ThisTask == 0){   
	      printf("All photons processed: dt_tot_all = %g, local photon-iterations: %d \n", 
		     dt_tot_all, bits);
	      fflush(stdout);
	  }
      }

/*       if (bits >= 0 && ThisTask == 0) */
/*       { */
/* 	  printf("bits = %d, Made iteration %d, total time remaining %g \n", */
/* 		 bits, N_iter, dt_tot_all ); */
/* 	  fflush(stdout); */
/*       } */


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
	  
	  for (k = 0 ; k < 3; k++)
	      SphP[i].ra.RadAccel[k] = 0.;

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
 	      acmax = 1.e4/dt;
 	      if (ac >= acmax)
 	      {
 		  printf("  RadAccel = (x,y,z) = (%g, %g, %g), dt = %g, mass = %g \n",
 			 SphP[i].ra.RadAccel[0], SphP[i].ra.RadAccel[1], SphP[i].ra.RadAccel[2], dt, P[i].Mass);
 		  fflush(stdout);
 		  for (k = 0 ; k < 3; k++)
 		  {
 		      SphP[i].ra.RadAccel[k] *= (acmax/ac);
 		  }
 	      }
 	  }
	  for (k = 0 ; k < 3; k++)
	  {
//	      if (P[i].Mass > 0) P[i].Vel[k] = (P[i].Vel[k] * P[i].Mass + 0. * SphP[i].i.Injected_BH_Momentum[k])/P[i].Mass;
	      SphP[i].i.Injected_BH_Momentum[k] = 0.;
	  }
      }
  }
  
#endif /* RAD_ACCEL */
      
/* ELIMINATION OF PHOTONS */

/* check if feedback/virtual particle has left an interesting sphere where I
 * track its motion. Or if it is "too old". In either case, remove it */


  int count_deleted_photon = 0, count_deleted_photon_total=0;
  double deleted_momentum = 0, total_deleted_momentum = 0;
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
  {
      if (P[i].Type == 3 && P[i].Mass != 0.) 
      {
	  r = sqrt(P[i].Pos[0]*P[i].Pos[0] + P[i].Pos[1]*P[i].Pos[1] + P[i].Pos[2]*P[i].Pos[2]);
	  if (r > All.OuterBoundary || All.Time-P[i].StellarAge > All.VirtualTime) 
	  {
	      P[i].Mass = 0.;
	      count_deleted_photon++;
	      deleted_momentum += P[i].OldPhotonMomentum;
	  } 
	  else  if (P[i].OldPhotonMomentum < 0.01 * P[i].BirthPhotonEnergy/C * All.UnitVelocity_in_cm_per_s) 
	  {
	      P[i].Mass = 0.;
	      count_deleted_photon++;
//            deleted_momentum += P[i].OldPhotonMomentum;
	  }
      }
  }

  fflush(diag_shc_photon);
  
  MPI_Allreduce(&count_deleted_photon, &count_deleted_photon_total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&deleted_momentum, &total_deleted_momentum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
  {
/* convert it to energy in cgs */
      total_deleted_momentum *= C/All.UnitVelocity_in_cm_per_s * All.UnitEnergy_in_cgs ;

      fprintf (diag_shc_photon_total,"%g %i %g\n",All.Time, count_deleted_photon_total, total_deleted_momentum);
      fflush(diag_shc_photon_total);
/* Also, print the number of new and deleted photons to the screen */
      printf("FB PARTICLES: spawned %d and destroyed %d photon particles \n", tot_spawned, count_deleted_photon_total);
      fflush(stdout);

      /* SHC BEGIN */
//    SysState.EnergyRadDeleted += total_deleted_momentum * C/All.UnitVelocity_in_cm_per_s;
      SysState.EnergyRadDeleted += total_deleted_momentum/All.UnitEnergy_in_cgs;
      /* SHC END */

  }
  

  
//fprintf (diag_shc,"End fb_particle\n");
//fflush (diag_shc);
  
  CPU_Step[CPU_RHD] += measure_time();

}


#endif /* closing FB_PARTCLES */
  
