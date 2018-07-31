#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"
#include "forcetree.h"


/* SN: Feb 2009.  This routine is to used when radiation field can be
 * considered optically thin. In this case no photon packets is emitted of
 * course. 
 */



#if defined(VIRTUAL) && defined(OTHIN_ACCELERATOR)

void sink_heating_othin(void)
/* cooling routine when star formation is enabled */
{
  int i, k, prev;
  double ascale = 1, hubble_a = 0, a3inv;
  double time_hubble_a, r2, TotalLuminosity, umin, umax;
  double EradFactor, RadEnergyDensity, 
      dtime, u_to_temp_fac, temp, uold, unew, dmax1, dmax2, diff, dudt_pdv, rhs, lhs, const4, 
      updv, u1, u2;

  EradFactor = 1./C/THOMPSON/All.VirtualCrosSection * PROTONMASS * All.UnitVelocity_in_cm_per_s *
      All.UnitVelocity_in_cm_per_s/All.UnitTime_in_s;
#ifdef DUST_OPACITY_FIT     
  double Tscale = 300., ustar, const5; /* kappa(T) = kappa0 (T/Tscale); kappa0 is from the insput file */
#endif

  

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


  CPU_Step[CPU_MISC] += measure_time();

/*   FirstSinkParticle = -1; */
/*   prev = -1; */
/*   int bh_num = 0; */
/*   for(i = 0; i < NumPart; i++) */
/*   { */
/*       if(P[i].Type == 5) */
/*       { */
/* 	  bh_num ++; */
/*           if(prev == -1) */
/*           { */
/*               FirstSinkParticle = i; */
/*           } */
	  
/*           if(prev >= 0) */
/*               NextSinkParticle[prev] = i; */
	  
/*           prev = i; */
/*       } */
/*   } */
  
/*   if(prev >= 0) */
/*   { */
/*       NextSinkParticle[prev] = -1; */
/*   } */
  
/*   if(ThisTask == 0) */
/*     { */
/*       printf("Found %d Blackholes \n", bh_num); */
/*       fflush(stdout); */
/*     } */


 
//----------------

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
	      TotalLuminosity = All.BH_Luminosity * All.UnitEnergy_in_cgs/All.UnitTime_in_s;

	      r2 = (P[i].Pos[0]-All.xbh) * (P[i].Pos[0]-All.xbh) +
		       (P[i].Pos[1]-All.ybh) * (P[i].Pos[1]-All.ybh) +
		       (P[i].Pos[2]-All.zbh) * (P[i].Pos[2]-All.zbh) +
		       All.SofteningBndry*All.SofteningBndry;

	      RadEnergyDensity =  TotalLuminosity/(4 * M_PI *r2 * All.UnitLength_in_cm * All.UnitLength_in_cm)/C;


#ifdef BACKGROUND_ILLUMINATION_TEMP
	      double T_bgr = All.EqTemp ;
	      RadEnergyDensity += pow(T_bgr, 4) * RAD_CONST;
#endif

//	      temp = pow(RadEnergyDensity/RAD_CONST, 0.25);
/* This is the blackbody temperature corresponding to the above radiation
 * density field  */
/* This is the repective u for the sph particle */

	      double k_factor = 1./EradFactor;

/* The 100 below assumes that radiation from the first core has 20 times
 * higher opacity than that of the T ~ 10 K envelope. This is a temporary fix
 * to be updated later  */

#ifdef DUST_OPACITY_FIT     
	      /* rad energy density of the central star at radius r */
	      ustar = TotalLuminosity/(4 * M_PI *r2 * All.UnitLength_in_cm * All.UnitLength_in_cm * C);
	      rhs = uold + dtime * (dudt_pdv  + k_factor * pow(T_bgr, 5)/Tscale * RAD_CONST +
				    k_factor * ustar);
#else
	      rhs = uold + dtime * (dudt_pdv  + k_factor * (pow(T_bgr, 4) * RAD_CONST + 25. * TotalLuminosity/
				   (4 * M_PI *r2 * All.UnitLength_in_cm * All.UnitLength_in_cm)/C));
#endif


/* using this k-factor ensures the correct units for the heating rate below in t_cool_estimate */
/* k_factor = THOMPSON * All.VirtualCrosSection/PROTONMASS * C *
   All.UnitTime_in_s/All.UnitVelocity_in_cm_per_s/All.UnitVelocity_in_cm_per_s;
*/

/* 		      unew = uold + dtime * (dudt_pdv  + k_factor *  
/* 			  (RadEnergyDensity - RAD_CONST*pow(uold * u_to_temp_fac ,4))); 
/* EXPLICIT COOLING sn -- 26.05.09 */


#ifndef DUST_OPACITY_FIT 
//	      rhs = uold + dtime * (dudt_pdv  + k_factor * RadEnergyDensity);
	      const4 = k_factor * RAD_CONST*pow(u_to_temp_fac ,4) * dtime;
	      u1 = rhs;
	      u2 = pow(rhs/const4, 0.25);
	      umin = u1;
	      if (u2 < umin) umin = u2;
	      umin = 0.1 * umin;
	      umax = u1;
	      if (u2 > umax) umax = u2;
	      
	      unew = sqrt(umin * umax);
	      
	      lhs = unew + const4 * pow(unew, 4);
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
		  lhs = unew + const4 * pow(unew, 4);
		  diff = fabs(lhs - rhs);
		  count++;
		  if (count > 300)
		  {
		      printf("\n Cannot find solution for T in othin_accelerator.c ... \n");
		      exit(1);
		  }
	      }
#else
	      const5 = k_factor * RAD_CONST*pow(u_to_temp_fac ,5) * dtime/Tscale;
	      u1 = rhs;
	      u2 = pow(rhs/const5, 0.2);
	      umin = u1;
	      if (u2 < umin) umin = u2;
	      umin = 0.1 * umin;
	      umax = u1;
	      if (u2 > umax) umax = u2;
	      
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
		      printf("\n Cannot find solution for T in othin_accelerator.c ... \n");
		      exit(1);
		  }
	      } 
#endif

	      if(ThisTask == 0 && unew*u_to_temp_fac > 1.e6)
	      {   
		  printf("SPH particle %d temperature too high T = %g \n",
			 i, unew * u_to_temp_fac);
		  fflush(stdout);
	      }

//------------------

	      SphP[i].e.DtEntropy = 
		  (unew * GAMMA_MINUS1/pow(SphP[i].d.Density * a3inv,
					   GAMMA_MINUS1) - SphP[i].Entropy) / dtime;

	  }
      }
  }


  CPU_Step[CPU_RHD] += measure_time();

}


#endif /* closing FB_PARTCLES */
  
