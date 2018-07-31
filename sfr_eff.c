#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"
#include "forcetree.h"


#ifdef COOLING

/*
 * This routine does cooling and star formation for
 * the effective multi-phase model.
 */

#ifndef MHM
#ifndef SFR			/* normal cooling routine when star formation is disabled */
void cooling_and_starformation(void)
{
  int i;
  double dt, dtime, hubble_a = 0, a3inv, ne = 1;
  double time_hubble_a, unew, dmax1, dmax2;

  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = hubble_function(All.Time);
      time_hubble_a = All.Time * hubble_a;
    }
  else
    {
      a3inv = time_hubble_a = hubble_a = 1;
    }

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(P[i].Type == 0)
	{
	  dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;
	  /*  the actual time-step */

	  if(All.ComovingIntegrationOn)
	    dtime = All.Time * dt / time_hubble_a;
	  else
	    dtime = dt;

	  ne = SphP[i].Ne;	/* electron abundance (gives ionization state and mean molecular weight) */

	  unew = DoCooling(DMAX(All.MinEgySpec,
				(SphP[i].Entropy + SphP[i].e.DtEntropy * dt) /
				GAMMA_MINUS1 * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1)),
			   SphP[i].d.Density * a3inv, dtime, &ne);

	  SphP[i].Ne = ne;

	  if(P[i].TimeBin)	/* upon start-up, we need to protect against dt==0 */
	    {
	      if(dt > 0)
		{

#ifdef COSMIC_RAYS
		  unew += CR_Particle_ThermalizeAndDissipate(SphP + i, dtime);
#endif /* COSMIC_RAYS */

		  SphP[i].e.DtEntropy = (unew * GAMMA_MINUS1 /
					 pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1) - SphP[i].Entropy) / dt;

		  if(SphP[i].e.DtEntropy < -0.5 * SphP[i].Entropy / dt)
		    SphP[i].e.DtEntropy = -0.5 * SphP[i].Entropy / dt;

		}
	    }
	}
    }
}

#else

void cooling_and_starformation(void)
/* cooling routine when star formation is enabled */
{
  int i, k, bin, flag, stars_spawned, tot_spawned, stars_converted, tot_converted, number_of_stars_generated;
  int Pebbles_created;
  unsigned int bits;
  double dt, dtime, ascale = 1, hubble_a = 0, a3inv, ne = 1;
  double time_hubble_a, unew, mass_of_star;
  double sum_sm, total_sm, sm, rate, sum_mass_stars, total_sum_mass_stars;
  double p, prob;
  double cloudmass;
  double factorEVP;
  double tsfr, trelax;
  double egyhot, egyeff, egycurrent, tcool, x, y, rate_in_msunperyear;
  double sfrrate, totsfrrate, dmax1, dmax2;

  // AH changes //

  double rhocore;
  double tidalrho, tidalrho_core, jeansrho;
  double rhocrit;

  // End of AH changes //

#ifdef WINDS
  int j;
  double v;
  double norm, dir[3];

#ifdef ISOTROPICWINDS
  double theta, phi;
#endif
#endif
#ifdef METALS
  double w;
#endif
#ifdef COSMIC_RAYS
  double tinj, instant_reheat;
#endif


#if defined(QUICK_LYALPHA) || defined(BH_THERMALFEEDBACK) || defined (BH_KINETICFEEDBACK) || defined(VIRTUAL_HEATING)
  double temp, u_to_temp_fac;

#ifdef CONSTANT_MEAN_MOLECULAR_WEIGHT
  u_to_temp_fac = All.MeanWeight  * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
#else
  u_to_temp_fac = 2.3  * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
#endif


//  u_to_temp_fac = (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1
//    * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
#endif

#ifdef FLTROUNDOFFREDUCTION
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 0)
      SphP[i].i.Injected_BH_Energy = FLT(SphP[i].i.dInjected_BH_Energy);
#endif

  for(bin = 0; bin < TIMEBINS; bin++)
    if(TimeBinActive[bin])
      TimeBinSfr[bin] = 0;

  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = hubble_function(All.Time);
      time_hubble_a = All.Time * hubble_a;
      ascale = All.Time;
    }
  else
    a3inv = ascale = time_hubble_a = 1;


  Pebbles_created = 0.;
  stars_spawned = stars_converted = 0;
  sum_sm = sum_mass_stars = 0;

#ifndef BH_FORM
  for(bits = 0; GENERATIONS > (1 << bits); bits++);
#endif

  //  double total_heating = 0.;
  //All.total_heating = 0.;

#if defined(DUST) && (DUST_ENERGY_CONSERVATION)
  double du_corr = 0.;
  if (SysState.MassComp[0])
    du_corr = All.total_delta_energy/SysState.MassComp[0];

  if(ThisTask == 0)
	{
	  printf("dE/E ==== %g, du_corr == == %g \n", All.total_delta_energy/All.total_energy_t0, du_corr);
	  fflush(stdout);
	}
#endif

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {

      /* Remove dust particles of very small mass */
       if(P[i].Type == 2)
	 if (P[i].Mass <= 1.e-5*All.OriginalGasMass)  P[i].Mass = 0.;

      if(P[i].Type == 0)
	{
	  dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;
	  /*  the actual time-step */

	  if(All.ComovingIntegrationOn)
	    dtime = All.Time * dt / time_hubble_a;
	  else
	    dtime = dt;

#ifdef INCREASE_SPH_MASS
	  //	  P[i].Mass *= exp(dt/All.DeltaTimeGrow);
	  P[i].Mass += All.OriginalGasMass * dt/All.DeltaTimeGrow;
#ifdef INCREASE_SPH_MASS_LIMIT
	  if (P[i].Mass >= 1.3*All.OriginalGasMass) P[i].Mass = 1.3*All.OriginalGasMass;
#endif
#endif

	  /* check whether conditions for star formation are fulfilled.
	   *  
	   * f=1  normal cooling
	   * f=0  star formation
	   */
	  flag = 1;		/* default is normal cooling */

#ifdef BH_FORM
#if defined(SGRA_POTENTIAL) || defined(CUSP_POTENTIAL) || defined(NFW_POTENTIAL) || defined(QUASAR_HEATING) ||  defined(FIND_SMBH)
		 double dx = P[i].Pos[0]-All.xbh;
		 double dy = P[i].Pos[1]-All.ybh;
		 double dz = P[i].Pos[2]-All.zbh;
#else
		 double dx = P[i].Pos[0];
		 double dy = P[i].Pos[1];
		 double dz = P[i].Pos[2];
#endif
		 double r2 = dx * dx + dy * dy + dz * dz;

#if defined(PLANET_IRRADIATION) || defined(PLANET_ACCRETION_FEEDBACK)
		 /* SN: this is for preheating by the planet */
		 double dxp = P[i].Pos[0]-All.xpl;
		 double dyp = P[i].Pos[1]-All.ypl;
		 double dzp = P[i].Pos[2]-All.zpl;
		 double r2p = dxp * dxp + dyp * dyp + dzp * dzp;
#endif


#ifndef JEANS_MASS_SF
		 if( SphP[i].d.Density >=
		     (All.CritPhysDensity * pow(All.UnitLength_in_cm, 3.)/ All.UnitMass_in_g))
//				+ 3./(2.*M_PI) / pow(r2, 1.5)*All.CritOverDensity) )
		     flag = 0;  
#else
#ifndef VIRTUAL_HEATING
		     double u_to_temp_fac, temp;
#endif		     

#ifdef CONSTANT_MEAN_MOLECULAR_WEIGHT
		     u_to_temp_fac = All.MeanWeight  * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
		     temp = u_to_temp_fac * (SphP[i].Entropy + SphP[i].e.DtEntropy * dt) /
		     GAMMA_MINUS1 * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);

#if defined(TIDAL_JEANS_SF) 
#ifndef CUSP_POTENTIAL
		     printf("\n TIDAL_JEANS_SF is calculated for CUSP_POTENTIAL only. Change the potential option in the Makefile. Exiting. \n");
		     exit(1);
#endif

		     // AH changes - implementing tidal density for the halo i.e. bh + flat core + Jaffe cusp //

		     rhocore = (All.CuspMass/(4.*M_PI))*(All.CuspRadius/(pow(All.FlatRadius, 2)*pow((All.FlatRadius+All.CuspRadius), 2)));

		     jeansrho = ((pow(M_PI * BOLTZMANN * temp/GRAVITY/All.MeanWeight/PROTONMASS, 3)/pow(10.*All.OriginalGasMass * All.UnitMass_in_g, 2))
				 * pow(All.UnitLength_in_cm, 3.)/ All.UnitMass_in_g);

		     tidalrho_core = (3.*All.SMBHmass)/(2.*M_PI*pow(sqrt(r2), 3)) + rhocore;

		     tidalrho = ((3.*All.CuspMass*All.CuspRadius)/(2.*M_PI*pow(sqrt(r2), 3)*(All.FlatRadius+All.CuspRadius))
		                 - (3.*All.CuspMass*All.CuspRadius)/(2.*M_PI*pow(sqrt(r2), 3)*(sqrt(r2)+All.CuspRadius))
				 - (3.*All.CuspMass*All.CuspRadius)/(4.*M_PI*pow(sqrt(r2), 2)*pow((sqrt(r2)+All.CuspRadius), 2))
				 + (2.*pow(All.FlatRadius, 3)*rhocore)/pow(sqrt(r2), 3) 
		                 + (3.*All.SMBHmass)/(2.*M_PI*pow(sqrt(r2), 3)));
		     
		     if(sqrt(r2) < All.FlatRadius) { //inside core

			 if(All.CritPhysDensity*tidalrho_core >= All.CritOverDensity*jeansrho) {

			     rhocrit = All.CritPhysDensity*tidalrho_core;

			 }
			 else {

			     rhocrit = All.CritOverDensity*jeansrho;

			 }

		     }
		     else { //outside core

			 if(All.CritPhysDensity*tidalrho >= All.CritOverDensity*jeansrho) {

			     rhocrit = All.CritPhysDensity*tidalrho;

			 }
			 else {

			     rhocrit = All.CritOverDensity*jeansrho;

			 }

		     }

		     //Now impose the density threshold

		     if(SphP[i].d.Density >= rhocrit)

#else //just the Jeans constraint

		     if( SphP[i].d.Density >= 
		     (All.CritPhysDensity 
		      + All.CritOverDensity * pow(M_PI * BOLTZMANN * temp/GRAVITY/All.MeanWeight/PROTONMASS, 3)/pow(0.01 * SOLAR_MASS, 2))
//		     + All.CritOverDensity * pow(M_PI * BOLTZMANN * temp/GRAVITY/All.MeanWeight/PROTONMASS, 3)/pow(10.*All.OriginalGasMass * All.UnitMass_in_g, 2))
		     * pow(All.UnitLength_in_cm, 3.)/ All.UnitMass_in_g)

#endif

#else
		     u_to_temp_fac = 2.3  * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
		     temp = u_to_temp_fac * (SphP[i].Entropy + SphP[i].e.DtEntropy * dt) /
		     GAMMA_MINUS1 * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);
	
#ifdef TIDAL_JEANS_SF	
#ifndef CUSP_POTENTIAL
		     printf("\n TIDAL_JEANS_SF is calculated for CUSP_POTENTIAL only. Change the potential option in the Makefile. Exiting. \n");
		     exit(1);
#endif
                     // AH changes - implementing tidal density for the halo i.e. bh + flat core + Jaffe cusp //

		     rhocore = (All.CuspMass/(4.*M_PI))*(All.CuspRadius/(pow(All.FlatRadius, 2)*pow((All.FlatRadius+All.CuspRadius), 2)));

		     jeansrho = ((pow(M_PI * BOLTZMANN * temp/GRAVITY/2.3/PROTONMASS, 3)/pow(10.*All.OriginalGasMass * All.UnitMass_in_g, 2))
				 * pow(All.UnitLength_in_cm, 3.)/ All.UnitMass_in_g)

		     tidalrho_core = (3.*All.SMBHmass)/(2.*M_PI*pow(sqrt(r2), 3)) + rhocore;

		     tidalrho = ((3.*All.CuspMass*All.CuspRadius)/(2.*M_PI*pow(sqrt(r2), 3)*(All.FlatRadius+All.CuspRadius))
		                 - (3.*All.CuspMass*All.CuspRadius)/(2.*M_PI*pow(sqrt(r2), 3)*(sqrt(r2)+All.CuspRadius))
				 - (3.*All.CuspMass*All.CuspRadius)/(4.*M_PI*pow(sqrt(r2), 2)*pow((sqrt(r2)+All.CuspRadius), 2))
				 + (2.*pow(All.FlatRadius, 3)*rhocore)/pow(sqrt(r2), 3) 
		                 + (3.*All.SMBHmass)/(2.*M_PI*pow(sqrt(r2), 3)));
		     
		     if(sqrt(r2) < All.FlatRadius) { //inside core

			 if(All.CritPhysDensity*tidalrho_core >= All.CritOverDensity*jeansrho) {

			     rhocrit = All.CritPhysDensity*tidalrho_core;

			 }
			 else {

			     rhocrit = All.CritOverDensity*jeansrho;

			 }

		     }
		     else { //outside core

			 if(All.CritPhysDensity*tidalrho >= All.CritOverDensity*jeansrho) {

			     rhocrit = All.CritPhysDensity*tidalrho;

			 }
			 else {

			     rhocrit = All.CritOverDensity*jeansrho;

			 }

		     }

		     //Now impose the density threshold

		     if(SphP[i].d.Density >= rhocrit)

#else //just the Jeans constraint

		     if( SphP[i].d.Density >= 
			 (All.CritPhysDensity 
//		         + All.CritOverDensity * pow(M_PI * BOLTZMANN * temp/GRAVITY/2.3/PROTONMASS, 3)/pow(0.04 * SOLAR_MASS, 2))
			 + All.CritOverDensity * pow(M_PI * BOLTZMANN * temp/GRAVITY/2.3/PROTONMASS, 3)/pow(10.*All.OriginalGasMass * All.UnitMass_in_g, 2))
		         * pow(All.UnitLength_in_cm, 3.)/ All.UnitMass_in_g)

#endif

#endif

/* SN: at temperature T and Jeans mass Mj, the collapse density is
   (pi k_b T/mu/G)^3 Mj^-2
 */
		     flag = 0;  
#endif /* JEANS_MASS_SF */

#else



	  if(SphP[i].d.Density * a3inv >= All.PhysDensThresh)
	    flag = 0;

	  if(All.ComovingIntegrationOn)
	    if(SphP[i].d.Density < All.OverDensThresh)
	      flag = 1;

#endif /* BH_FORM */


/* #ifdef VIRTUAL */
/* #ifndef RAD_ACCEL */
/* 	  if (P[i].Mass > 0) */
/* 	  { */
/* /\* 	      double injected_mom1 = 0.; *\/ */
/* /\* 	      for (k = 0; k < 3; k++) *\/ */
/* /\* 	      { *\/ */
/* /\* 		  injected_mom1 += SphP[i].i.Injected_BH_Momentum[k]* SphP[i].i.Injected_BH_Momentum[k]; *\/ */
/* /\* 	      } *\/ */
/* /\* 	      injected_mom1 = sqrt(injected_mom1); *\/ */
/* /\* 	      if (injected_mom1 > 0)  *\/ */
/* /\* 	      { *\/ */
/* /\* 		  printf("SPH INJECTED MoM = %g\n",injected_mom1); *\/ */
/* /\* 		  fflush(stdout); *\/ */
/* /\* 	      } *\/ */
/* 	      for (k = 0; k < 3; k++) */
/* 	      { */
/* 		P[i].Vel[k] = (P[i].Vel[k] * P[i].Mass + SphP[i].i.Injected_BH_Momentum[k])/P[i].Mass; */
/* 		  /\* if (i%3000 == 0)   *\/ */
/* 		  /\*   { *\/ */
/* 		  /\*     printf("k = %d, SPH INJECTED V_kick = %g\n", k, SphP[i].i.Injected_BH_Momentum[k]/P[i].Mass);  *\/ */
/* 		  /\*     fflush(stdout); *\/ */
/* 		  /\*   } *\/ */
/* 		  SphP[i].i.Injected_BH_Momentum[k] = 0.; */
/* 	      } */
/* 	  } */
/* #endif */
/* #endif */


#ifdef BLACK_HOLES
	  if(P[i].Mass == 0)
	    flag = 1;
#endif

#ifdef WINDS
	  if(SphP[i].DelayTime > 0)
	    flag = 1;		/* only normal cooling for particles in the wind */

	  if(SphP[i].DelayTime > 0)
	    SphP[i].DelayTime -= dtime;

	  if(SphP[i].DelayTime > 0)
	    if(SphP[i].d.Density * a3inv < All.WindFreeTravelDensFac * All.PhysDensThresh)
	      SphP[i].DelayTime = 0;

	  if(SphP[i].DelayTime < 0)
	    SphP[i].DelayTime = 0;

#endif


#ifdef QUICK_LYALPHA
	  temp = u_to_temp_fac * (SphP[i].Entropy + SphP[i].e.DtEntropy * dt) /
	    GAMMA_MINUS1 * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);

	  if(SphP[i].d.Density > All.OverDensThresh && temp < 1.0e5)
	    flag = 0;
	  else
	    flag = 1;
#endif


#if !defined(NOISMPRESSURE) && !defined(QUICK_LYALPHA)
	  if(flag == 1)		/* normal implicit isochoric cooling */
#endif
	    {
	      SphP[i].Sfr = 0;
#if defined(COSMIC_RAYS) && defined(CR_OUTPUT_INJECTION)
	      SphP[i].CR_Specific_SupernovaHeatingRate = 0;
#endif
	      ne = SphP[i].Ne;	/* electron abundance (gives ionization state and mean molecular weight) */

	      unew = DMAX(All.MinEgySpec,
			  (SphP[i].Entropy + SphP[i].e.DtEntropy * dt) /
			  GAMMA_MINUS1 * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1));

#ifdef DUST
	      //	      if (All.NumCurrentTiStep = 0)
	      if (SphP[i].dh.DragHeating && P[i].Mass != 0)
	      	{
		  double du_dust= SphP[i].dh.DragHeating/P[i].Mass*dt;
#ifdef DUST_ENERGY_CONSERVATION
		  /* The energy conservation is thus expected to be re-established on the following time scale */
		  unew += du_dust - du_corr*dt/0.003/All.TimeMax;
#else
		  /* if (abs(du_dust) >= 0.* unew) */
		  /*   { */
		  /*     printf("unew %g, du_dust  %g \n",unew, du_dust); */
		  /*     fflush(stdout); */
		  /*   } */
		  unew += du_dust;
#endif
	      	  SphP[i].dh.DragHeating = 0.;
	      	}
#endif


#if defined(BH_THERMALFEEDBACK) || defined(BH_KINETICFEEDBACK) 
	      double dT = 0.;
	      if(SphP[i].i.Injected_BH_Energy)
		{
		  if(P[i].Mass == 0)
		    SphP[i].i.Injected_BH_Energy = 0;
		  else
		    unew += SphP[i].i.Injected_BH_Energy / P[i].Mass;

		  temp = u_to_temp_fac * unew;
		  /*		  dT = u_to_temp_fac * SphP[i].i.Injected_BH_Energy / P[i].Mass;
				  if (dT > 0.01) printf("T = %g, dT = %g  \n", temp, dT);*/


		  if(temp > 5.0e9)
		    unew = 5.0e9 / u_to_temp_fac;
#ifdef FLTROUNDOFFREDUCTION
		  SphP[i].i.dInjected_BH_Energy = 0;
#else
		  SphP[i].i.Injected_BH_Energy = 0;
#endif
		}
#endif


#ifndef VIRTUAL_HEATING /* VIRTUAL_HEATING: heating-cooling is calculated in fb_particles.c */


#if defined(PLANET_IRRADIATION) || defined(PLANET_ACCRETION_FEEDBACK)
	      unew = DoCooling(unew, SphP[i].d.Density * a3inv, dtime, &ne, r2, r2p);
#else
#ifndef DUST_TWO_POPULATIONS
	      unew = DoCooling(unew, SphP[i].d.Density * a3inv, dtime, &ne, r2);
#else
	      double zloc = P[i].NewDensity/SphP[i].d.Density;
	      double a_loc = 0.; /* mean pebble size at the SPH particle  */
	      if (zloc > 0)
		{
		  a_loc = P[i].NewKappa/P[i].NewDensity; /* mean pebble size at the SPH particle  */
		  //		  dMicro_dust_mass =  - P[i].d_MicroDustMass/P[i].NewDensity;
		}
	      //	      dMicro_dust_mass =  P[i].d_MicroDustMass;
	      P[i].MicroDustMass += P[i].d_MicroDustMass;
	      P[i].d_MicroDustMass = 0.;
	      double micro_dust_enhancement = P[i].MicroDustMass/P[i].Mass/All.InitialMicroDustZ;
	      /* Protect against dust mass goind negative */
	      if (P[i].MicroDustMass <= 1.e-10*P[i].Mass) 
		{
		  /*		  printf("i = %d, zloc   =  %g    Micro Mass/M = %g, z/R = %g z_micro_enhancement = %g dMicro/Mass = %g \n", i,  
			 zloc, P[i].MicroDustMass/P[i].Mass,  micro_dust_enhancement, 
			 dMicro_dust_mass/P[i].Mass);*/
		  P[i].MicroDustMass = 1.e-10*P[i].Mass;
		}


	      if (micro_dust_enhancement < 1./All.BetaCool) micro_dust_enhancement = 1./All.BetaCool;
	      /* if (ThisTask == 0 && i%100000000 == 0)  */
	      /* 	{ */
	      /* 	  printf("i = %d, Micro Mass = %g, z/R = %g z_micro_enhancement = %g dMicro/MMicro = %g R = %g \n", i,   */
	      /* 		 zloc, P[i].MicroDustMass,  micro_dust_enhancement,  */
	      /* 		 dMicro_dust_mass/P[i].Mass, pow(r2, 0.5)); */
	      /* 	} */

	      unew = DoCooling(unew, SphP[i].d.Density * a3inv, dtime, &ne, r2, micro_dust_enhancement); 



#endif /* end of DUST_TWO_POPULATIONS */
#endif



	      if(P[i].TimeBin)	/* upon start-up, we need to protect against dt==0 */
		{
		  /* note: the adiabatic rate has been already added in ! */

		  if(dt > 0)
		    {
#ifdef COSMIC_RAYS
		      unew += CR_Particle_ThermalizeAndDissipate(SphP + i, dtime);
#endif /* COSMIC_RAYS */

		      SphP[i].e.DtEntropy = (unew * GAMMA_MINUS1 /
					     pow(SphP[i].d.Density * a3inv,
						 GAMMA_MINUS1) - SphP[i].Entropy) / dt;
#ifndef DO_NOT_PROTECT
		      if(SphP[i].e.DtEntropy < -0.5 * SphP[i].Entropy / dt)
			SphP[i].e.DtEntropy = -0.5 * SphP[i].Entropy / dt;
#else
		      if(SphP[i].e.DtEntropy < -0.9 * SphP[i].Entropy / dt)
			  SphP[i].e.DtEntropy = -0.9 * SphP[i].Entropy / dt;
#endif

		    }
		}
#endif /* VIRTUAL_HEATING: heating-cooling is calculated in fb_particles.c */
	    }

	  if(flag == 0)		/* active star formation */
	    {
#ifdef BH_FORM
		/* here we turn the gas particle itself into a star */
		Stars_converted++;
		stars_converted++;
		
		sum_mass_stars += P[i].Mass;
		
		P[i].Type = 5;
/* SN: This negligible mass change helps to enforce mergers in blackhole.c
 */
		P[i].Mass *= (0.995 + 0.005 * gsl_rng_uniform(random_generator));
/* SN: Dec11.2008: make sure that BH mass is initialised equal to current mass
 * of the sink particle */
		P[i].BH_Mass = P[i].Mass;
#ifdef LIMITED_ACCRETION_AND_FEEDBACK
/* SN. BH_Mass describes the mass that already released the feedback, whereas
 * disc mass is the reservoir mass that will be accreted/turned into the
 * BH_Mass after viscous time. When sink particle is created, it should be
 * initialised with 0 BH mass so that feedback could be produced.
 */
		P[i].AccDisc_Mass = P[i].Mass;
		P[i].BH_Mass = 0;
#endif

//		P[i].SwallowID = 0;
/* end of SN changes */
		TimeBinCountSph[P[i].TimeBin]--;
		TimeBinSfr[P[i].TimeBin] -= SphP[i].Sfr;
#else /* BH_FORM */  
/* use Volker's star formation */
#if !defined(QUICK_LYALPHA)
	      tsfr = sqrt(All.PhysDensThresh / (SphP[i].d.Density * a3inv)) * All.MaxSfrTimescale;

	      factorEVP = pow(SphP[i].d.Density * a3inv / All.PhysDensThresh, -0.8) * All.FactorEVP;

	      egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

	      ne = SphP[i].Ne;

/* Addition by CBP. 27/07/2009. Added NFW_POTENTIAL. */

#if defined(SGRA_POTENTIAL) || defined(CUSP_POTENTIAL) || defined(SIS_POTENTIAL) || defined(NFW_POTENTIAL)
	      double dx = P[i].Pos[0]-All.xbh;
	      double dy = P[i].Pos[1]-All.ybh;
	      double dz = P[i].Pos[2]-All.zbh;
#else
	      double dx = P[i].Pos[0];
	      double dy = P[i].Pos[1];
	      double dz = P[i].Pos[2];
#endif
	      double r2 = dx * dx + dy * dy + dz * dz;

#ifdef DUST_TWO_POPULATIONS
	      tcool = 1.e30;
#else
	      tcool = GetCoolingTime(egyhot, SphP[i].d.Density * a3inv, &ne);
	      SphP[i].Ne = ne;
#endif


	      y =
		tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);

	      x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

	      egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

	      cloudmass = x * P[i].Mass;

	      if(tsfr < dtime)
		tsfr = dtime;

	      sm = (1 - All.FactorSN) * dtime / tsfr * cloudmass;	/* amount of stars expect to form */

	      p = sm / P[i].Mass;

	      sum_sm += P[i].Mass * (1 - exp(-p));


	      if(dt > 0)
		{
		  if(P[i].TimeBin)	/* upon start-up, we need to protect against dt==0 */
		    {
		      trelax = tsfr * (1 - x) / x / (All.FactorSN * (1 + factorEVP));
		      egycurrent =
			SphP[i].Entropy * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;

#ifdef COSMIC_RAYS
		      if(All.CR_SNEff > 0)
			{
			  tinj = SphP[i].CR_E0 / (p * All.FeedbackEnergy * All.CR_SNEff / dtime);

			  instant_reheat =
			    CR_Particle_SupernovaFeedback(&SphP[i], p * All.FeedbackEnergy * All.CR_SNEff,
							  tinj);
			}
		      else
			instant_reheat = 0;

#if defined(COSMIC_RAYS) && defined(CR_OUTPUT_INJECTION)
		      SphP[i].CR_Specific_SupernovaHeatingRate =
			(p * All.FeedbackEnergy * All.CR_SNEff - instant_reheat) / dtime;
#endif

		      egycurrent += instant_reheat;

		      egycurrent += CR_Particle_ThermalizeAndDissipate(SphP + i, dtime);
#endif /* COSMIC_RAYS */


#if defined(BH_THERMALFEEDBACK) || defined(BH_KINETICFEEDBACK) || defined(VIRTUAL)
		      if (P[i].Mass > 0)
		      {
			  for (k = 0, ; k < 3; k++)
			  {
			      P[i].Vel[k] = (P[i].Vel[k] * P[i].Mass + SphP[i].i.Injected_BH_Momentum[k]) /P[i].Mass;
			      SphP[i].i.Injected_BH_Momentum[k] = 0.;
			  }
		      }
#endif


#if !defined(NOISMPRESSURE)
		      SphP[i].Entropy =
			(egyeff +
			 (egycurrent -
			  egyeff) * exp(-dtime / trelax)) * GAMMA_MINUS1 /
			pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);

		      SphP[i].e.DtEntropy = 0;
#endif
		    }
		}



	      /* the upper bits of the gas particle ID store how man stars this gas
	         particle gas already generated */

	      if(bits == 0)
		number_of_stars_generated = 0;
	      else
		number_of_stars_generated = (P[i].ID >> (32 - bits));

	      mass_of_star = P[i].Mass / (GENERATIONS - number_of_stars_generated);


	      SphP[i].Sfr = (1 - All.FactorSN) * cloudmass / tsfr *
		(All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

	      TimeBinSfr[P[i].TimeBin] += SphP[i].Sfr;
#ifdef METALS
	      w = get_random_number(P[i].ID);
	      P[i].Metallicity += w * METAL_YIELD * (1 - exp(-p));
#endif

	      prob = P[i].Mass / mass_of_star * (1 - exp(-p));

#else /* belongs to ifndef(QUICK_LYALPHA) */

	      prob = 2.0;	/* this will always cause a star creation event */

	      if(bits == 0)
		number_of_stars_generated = 0;
	      else
		number_of_stars_generated = (P[i].ID >> (32 - bits));

	      mass_of_star = P[i].Mass / (GENERATIONS - number_of_stars_generated);

	      SphP[i].Sfr = 0;

#endif /* ends to QUICK_LYALPHA */

	      if(get_random_number(P[i].ID + 1) < prob)	/* ok, make a star */
		{
		  if(number_of_stars_generated == (GENERATIONS - 1))
		    {
		      /* here we turn the gas particle itself into a star */
		      Stars_converted++;
		      stars_converted++;

		      sum_mass_stars += P[i].Mass;

		      P[i].Type = 4;
		      TimeBinCountSph[P[i].TimeBin]--;
		      TimeBinSfr[P[i].TimeBin] -= SphP[i].Sfr;

#ifdef STELLARAGE
		      P[i].StellarAge = All.Time;
#endif
		    }
		  else
		    {
		      /* here we spawn a new star particle */

		      if(NumPart + stars_spawned >= All.MaxPart)
			{
			  printf
			    ("On Task=%d with NumPart=%d we try to spawn %d particles. Sorry, no space left...(All.MaxPart=%d)\n",
			     ThisTask, NumPart, stars_spawned, All.MaxPart);
			  fflush(stdout);
			  endrun(8888);
			}

		      P[NumPart + stars_spawned] = P[i];
		      P[NumPart + stars_spawned].Type = 4;

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

		      P[i].ID += (1 << (32 - bits));

		      P[NumPart + stars_spawned].Mass = mass_of_star;
//		      P[i].Mass -= P[NumPart + stars_spawned].Mass;
		      sum_mass_stars += P[NumPart + stars_spawned].Mass;
#ifdef STELLARAGE
		      P[NumPart + stars_spawned].StellarAge = All.Time;
#endif
		      force_add_star_to_tree(i, NumPart + stars_spawned);

		      stars_spawned++;
		    }
		}

#ifdef METALS
	      if(P[i].Type == 0)	/* to protect using a particle that has been turned into a star */
		P[i].Metallicity += (1 - w) * METAL_YIELD * (1 - exp(-p));
#endif



#ifdef WINDS
	      /* Here comes the wind model */

	      if(P[i].Type == 0)	/* to protect using a particle that has been turned into a star */
		{
		  p = All.WindEfficiency * sm / P[i].Mass;

		  prob = 1 - exp(-p);

		  if(get_random_number(P[i].ID + 2) < prob)	/* ok, make the particle go into the wind */
		    {
		      v =
			sqrt(2 * All.WindEnergyFraction * All.FactorSN *
			     All.EgySpecSN / (1 - All.FactorSN) / All.WindEfficiency);
#ifdef ISOTROPICWINDS
		      theta = acos(2 * get_random_number(P[i].ID + 3) - 1);
		      phi = 2 * M_PI * get_random_number(P[i].ID + 4);

		      dir[0] = sin(theta) * cos(phi);
		      dir[1] = sin(theta) * sin(phi);
		      dir[2] = cos(theta);
#else
		      dir[0] = P[i].g.GravAccel[1] * P[i].Vel[2] - P[i].g.GravAccel[2] * P[i].Vel[1];
		      dir[1] = P[i].g.GravAccel[2] * P[i].Vel[0] - P[i].g.GravAccel[0] * P[i].Vel[2];
		      dir[2] = P[i].g.GravAccel[0] * P[i].Vel[1] - P[i].g.GravAccel[1] * P[i].Vel[0];
#endif

		      for(j = 0, norm = 0; j < 3; j++)
			norm += dir[j] * dir[j];

		      norm = sqrt(norm);
		      if(get_random_number(P[i].ID + 5) < 0.5)
			norm = -norm;

		      if(norm != 0)
			{
			  for(j = 0; j < 3; j++)
			    dir[j] /= norm;

			  for(j = 0; j < 3; j++)
			    {
			      P[i].Vel[j] += v * ascale * dir[j];
			      SphP[i].VelPred[j] += v * ascale * dir[j];
			    }

			  SphP[i].DelayTime = All.WindFreeTravelLength / v;
			}
		    }
		}
#endif
#endif /* BH_FORM */
	    }
	}

    }				/* end of main loop over active particles */




  MPI_Allreduce(&stars_spawned, &tot_spawned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&stars_converted, &tot_converted, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(tot_spawned > 0 || tot_converted > 0)
    {
      if(ThisTask == 0)
	{
	  printf("SFR: spawned %d stars, converted %d gas particles into stars\n",
		 tot_spawned, tot_converted);
	  fflush(stdout);
	}


      All.TotNumPart += tot_spawned;
      All.TotN_gas -= tot_converted;
      NumPart += stars_spawned;

      /* Note: N_gas is only reduced once rearrange_particle_sequence is called */

      /* Note: New tree construction can be avoided because of  `force_add_star_to_tree()' */
    }

  for(bin = 0, sfrrate = 0; bin < TIMEBINS; bin++)
    if(TimeBinCount[bin])
      sfrrate += TimeBinSfr[bin];

  MPI_Allreduce(&sfrrate, &totsfrrate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Reduce(&sum_sm, &total_sm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sum_mass_stars, &total_sum_mass_stars, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ThisTask == 0)
    {
      if(All.TimeStep > 0)
	rate = total_sm / (All.TimeStep / time_hubble_a);
      else
	rate = 0;

      /* convert to solar masses per yr */

      rate_in_msunperyear = rate * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

      fprintf(FdSfr, "%g %g %g %g %g\n", All.Time, total_sm, totsfrrate, rate_in_msunperyear,
	      total_sum_mass_stars);
      fflush(FdSfr);
    }
}

double get_starformation_rate(int i)
{
  double rateOfSF;
  double a3inv;
  int flag;
  double tsfr;
  double factorEVP, egyhot, ne, tcool, y, x, cloudmass;



  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;


  flag = 1;			/* default is normal cooling */

  if(SphP[i].d.Density * a3inv >= All.PhysDensThresh)
    flag = 0;

  if(All.ComovingIntegrationOn)
    if(SphP[i].d.Density < All.OverDensThresh)
      flag = 1;

  if(flag == 1)
    return 0;

  tsfr = sqrt(All.PhysDensThresh / (SphP[i].d.Density * a3inv)) * All.MaxSfrTimescale;

  factorEVP = pow(SphP[i].d.Density * a3inv / All.PhysDensThresh, -0.8) * All.FactorEVP;

  egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

  ne = SphP[i].Ne;

//#if defined(QUASAR_HEATING)  || defined(BETA_COOLING) || defined(ISOTHERM)
#if defined(VIRTUAL)
//  exit(0);
#else
  tcool = GetCoolingTime(egyhot, SphP[i].d.Density * a3inv, &ne);
#endif

  y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);

  x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

  cloudmass = x * P[i].Mass;

  rateOfSF = (1 - All.FactorSN) * cloudmass / tsfr;

  /* convert to solar masses per yr */

  rateOfSF *= (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

  return rateOfSF;
}

#endif /* closing of SFR-conditional */



#endif /* closes SFR_MHM conditional */





#if defined(SFR) || defined(BLACK_HOLES)
void rearrange_particle_sequence(void)
{
  int i, j, flag = 0, flag_sum;
  struct particle_data psave;

#ifdef BLACK_HOLES
  int count_elim, count_gaselim, tot_elim, tot_gaselim;
#ifdef DUST
  int count_dustelim, tot_dustelim;
#endif
#endif

#ifdef SFR
  if(Stars_converted)
    {
      N_gas -= Stars_converted;
      Stars_converted = 0;

      for(i = 0; i < N_gas; i++)
	if(P[i].Type != 0)
	  {
	    for(j = N_gas; j < NumPart; j++)
	      if(P[j].Type == 0)
		break;

	    if(j >= NumPart)
	      endrun(181170);

	    psave = P[i];
	    P[i] = P[j];
	    SphP[i] = SphP[j];
	    P[j] = psave;
	  }
      flag = 1;
    }
#endif

#ifdef BLACK_HOLES
  count_elim = 0;
  count_gaselim = 0;
#ifdef DUST
  count_dustelim = 0;
#endif

  for(i = 0; i < NumPart; i++)
    if(P[i].Mass == 0)
      {
	TimeBinCount[P[i].TimeBin]--;

	if(TimeBinActive[P[i].TimeBin])
	  NumForceUpdate--;

	if(P[i].Type == 0)
	  {
	    TimeBinCountSph[P[i].TimeBin]--;

	    P[i] = P[N_gas - 1];
	    SphP[i] = SphP[N_gas - 1];

	    P[N_gas - 1] = P[NumPart - 1];

	    N_gas--;

	    count_gaselim++;
	  }
#ifdef DUST
	else
	  {
          if(P[i].Type == 2)  count_dustelim++;
	    P[i] = P[NumPart - 1];
	  }
#endif

	NumPart--;
	i--;

	count_elim++;
      }

  MPI_Allreduce(&count_elim, &tot_elim, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&count_gaselim, &tot_gaselim, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#ifdef DUST
  MPI_Allreduce(&count_dustelim, &tot_dustelim, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

  if(count_elim)
    flag = 1;

  if(ThisTask == 0)
    {
#ifdef DUST
      printf("Blackholes: Eliminated %d gas and %d dust particles and merged away %d black holes.\n",
	     tot_gaselim, tot_dustelim, tot_elim - tot_gaselim-tot_dustelim);
#else
      printf("Blackholes: Eliminated %d gas particles and merged away %d black holes.\n",
	     tot_gaselim, tot_elim - tot_gaselim);
#endif
      fflush(stdout);
    }

  All.TotNumPart -= tot_elim;
  All.TotN_gas -= tot_gaselim;
  All.TotBHs -= tot_elim - tot_gaselim;
#endif

  MPI_Allreduce(&flag, &flag_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);


  if(flag_sum)
    reconstruct_timebins();

}
#endif /* closing of SFR-conditional */



#if defined(SFR)
void init_clouds(void)
{
  double A0, dens, tcool, ne, coolrate, egyhot, x, u4, meanweight;
  double tsfr, y, peff, fac, neff, egyeff, factorEVP, sigma, thresholdStarburst;

  if(All.PhysDensThresh == 0)
    {
      A0 = All.FactorEVP;

      egyhot = All.EgySpecSN / A0;

      meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */

#ifdef CONSTANT_MEAN_MOLECULAR_WEIGHT
  meanweight = All.MeanWeight;
#endif

      u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
      u4 *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;


      if(All.ComovingIntegrationOn)
	dens = 1.0e6 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
      else
	dens = 1.0e6 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

      if(All.ComovingIntegrationOn)
	{
	  All.Time = 1.0;	/* to be guaranteed to get z=0 rate */
	  IonizeParams();
	}

      ne = 1.0;
      SetZeroIonization();
//#if defined(QUASAR_HEATING)  || defined(BETA_COOLING) || defined(ISOTHERM)
#ifdef VIRTUAL
//      exit(0);
#else
     tcool = GetCoolingTime(egyhot, dens, &ne);
#endif

      coolrate = egyhot / tcool / dens;

      x = (egyhot - u4) / (egyhot - All.EgySpecCold);

      All.PhysDensThresh =
	x / pow(1 - x,
		2) * (All.FactorSN * All.EgySpecSN - (1 -
						      All.FactorSN) * All.EgySpecCold) /
	(All.MaxSfrTimescale * coolrate);

      if(ThisTask == 0)
	{
	  printf("\nA0= %g  \n", A0);
	  printf("Computed: PhysDensThresh= %g  (int units)         %g h^2 cm^-3\n", All.PhysDensThresh,
		 All.PhysDensThresh / (PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs));
	  printf("EXPECTED FRACTION OF COLD GAS AT THRESHOLD = %g\n\n", x);
	  printf("tcool=%g dens=%g egyhot=%g\n", tcool, dens, egyhot);
	}

      dens = All.PhysDensThresh * 10;

      do
	{
	  tsfr = sqrt(All.PhysDensThresh / (dens)) * All.MaxSfrTimescale;
	  factorEVP = pow(dens / All.PhysDensThresh, -0.8) * All.FactorEVP;
	  egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

	  ne = 0.5;
//#if defined(QUASAR_HEATING)  || defined(BETA_COOLING) || defined(ISOTHERM)
#ifdef VIRTUAL
//      exit(0);
#else
     tcool = GetCoolingTime(egyhot, dens, &ne);
#endif
	  y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
	  x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
	  egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

	  peff = GAMMA_MINUS1 * dens * egyeff;

	  fac = 1 / (log(dens * 1.025) - log(dens));
	  dens *= 1.025;

	  neff = -log(peff) * fac;

	  tsfr = sqrt(All.PhysDensThresh / (dens)) * All.MaxSfrTimescale;
	  factorEVP = pow(dens / All.PhysDensThresh, -0.8) * All.FactorEVP;
	  egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

	  ne = 0.5;
//#if defined(QUASAR_HEATING)  || defined(BETA_COOLING) || defined(ISOTHERM)
#ifdef VIRTUAL
//      exit(0);
#else
     tcool = GetCoolingTime(egyhot, dens, &ne);
#endif
	  y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
	  x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
	  egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

	  peff = GAMMA_MINUS1 * dens * egyeff;

	  neff += log(peff) * fac;
	}
      while(neff > 4.0 / 3);

      thresholdStarburst = dens;

#ifdef MODIFIEDBONDI
      All.BlackHoleRefDensity = thresholdStarburst;
      All.BlackHoleRefSoundspeed = sqrt(GAMMA * GAMMA_MINUS1 * egyeff);
#endif


      if(ThisTask == 0)
	{
	  printf("Run-away sets in for dens=%g\n", thresholdStarburst);
	  printf("Dynamic range for quiescent star formation= %g\n", thresholdStarburst / All.PhysDensThresh);
	  fflush(stdout);
	}

      integrate_sfr();

      if(ThisTask == 0)
	{
	  sigma = 10.0 / All.Hubble * 1.0e-10 / pow(1.0e-3, 2);

	  printf("Isotherm sheet central density: %g   z0=%g\n",
		 M_PI * All.G * sigma * sigma / (2 * GAMMA_MINUS1) / u4,
		 GAMMA_MINUS1 * u4 / (2 * M_PI * All.G * sigma));
	  fflush(stdout);

	}

      if(All.ComovingIntegrationOn)
	{
	  All.Time = All.TimeBegin;
	  IonizeParams();
	}

#ifdef WINDS
      if(All.WindEfficiency > 0)
	if(ThisTask == 0)
	  printf("Windspeed: %g\n",
		 sqrt(2 * All.WindEnergyFraction * All.FactorSN * All.EgySpecSN / (1 - All.FactorSN) /
		      All.WindEfficiency));
#endif
    }
}

void integrate_sfr(void)
{
  double rho0, rho, rho2, q, dz, gam, sigma = 0, sigma_u4, sigmasfr = 0, ne, P1;
  double x = 0, y, P, P2, x2, y2, tsfr2, factorEVP2, egyhot2, tcool2, drho, dq;
  double meanweight, u4, z, tsfr, tcool, egyhot, factorEVP, egyeff, egyeff2;
  FILE *fd;


  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */
#ifdef CONSTANT_MEAN_MOLECULAR_WEIGHT
  meanweight = All.MeanWeight;
#endif
  u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
  u4 *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  if(All.ComovingIntegrationOn)
    {
      All.Time = 1.0;		/* to be guaranteed to get z=0 rate */
      IonizeParams();
    }

  if(ThisTask == 0)
    fd = fopen("eos.txt", "w");
  else
    fd = 0;

  for(rho = All.PhysDensThresh; rho <= 1000 * All.PhysDensThresh; rho *= 1.1)
    {
      tsfr = sqrt(All.PhysDensThresh / rho) * All.MaxSfrTimescale;

      factorEVP = pow(rho / All.PhysDensThresh, -0.8) * All.FactorEVP;

      egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

      ne = 1.0;
//#if defined(QUASAR_HEATING)  || defined(BETA_COOLING) || defined(ISOTHERM)
#ifdef VIRTUAL
//      exit(0);
#else
     tcool =  GetCoolingTime(egyhot, rho, &ne);
#endif
      
      y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
      x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

      egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

      P = GAMMA_MINUS1 * rho * egyeff;

      if(ThisTask == 0)
	{
	  fprintf(fd, "%g %g\n", rho, P);
	}
    }

  if(ThisTask == 0)
    fclose(fd);


  if(ThisTask == 0)
    fd = fopen("sfrrate.txt", "w");
  else
    fd = 0;

  for(rho0 = All.PhysDensThresh; rho0 <= 10000 * All.PhysDensThresh; rho0 *= 1.02)
    {
      z = 0;
      rho = rho0;
      q = 0;
      dz = 0.001;

      sigma = sigmasfr = sigma_u4 = 0;

      while(rho > 0.0001 * rho0)
	{
	  if(rho > All.PhysDensThresh)
	    {
	      tsfr = sqrt(All.PhysDensThresh / rho) * All.MaxSfrTimescale;

	      factorEVP = pow(rho / All.PhysDensThresh, -0.8) * All.FactorEVP;

	      egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

	      ne = 1.0;
//#if defined(QUASAR_HEATING)  || defined(BETA_COOLING) || defined(ISOTHERM)
#ifdef VIRTUAL
//	      exit(0);
#else
	      tcool = GetCoolingTime(egyhot, rho, &ne);
#endif
		   
	      y =
		tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
	      x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

	      egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

	      P = P1 = GAMMA_MINUS1 * rho * egyeff;

	      rho2 = 1.1 * rho;
	      tsfr2 = sqrt(All.PhysDensThresh / rho2) * All.MaxSfrTimescale;
	      factorEVP2 = pow(rho2 / All.PhysDensThresh, -0.8) * All.FactorEVP;
	      egyhot2 = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;
//#if defined(QUASAR_HEATING)  || defined(BETA_COOLING) || defined(ISOTHERM)
#ifdef VIRTUAL
//	      exit(0);
#else
	      tcool2 = GetCoolingTime(egyhot2, rho2, &ne);
#endif
	      y2 =
		tsfr2 / tcool2 * egyhot2 / (All.FactorSN * All.EgySpecSN -
					    (1 - All.FactorSN) * All.EgySpecCold);
	      x2 = 1 + 1 / (2 * y2) - sqrt(1 / y2 + 1 / (4 * y2 * y2));
	      egyeff2 = egyhot2 * (1 - x2) + All.EgySpecCold * x2;
	      P2 = GAMMA_MINUS1 * rho2 * egyeff2;

	      gam = log(P2 / P1) / log(rho2 / rho);
	    }
	  else
	    {
	      tsfr = 0;

	      P = GAMMA_MINUS1 * rho * u4;
	      gam = 1.0;


	      sigma_u4 += rho * dz;
	    }



	  drho = q;
	  dq = -(gam - 2) / rho * q * q - 4 * M_PI * All.G / (gam * P) * rho * rho * rho;

	  sigma += rho * dz;
	  if(tsfr > 0)
	    {
	      sigmasfr += (1 - All.FactorSN) * rho * x / tsfr * dz;
	    }

	  rho += drho * dz;
	  q += dq * dz;
	}


      sigma *= 2;		/* to include the other side */
      sigmasfr *= 2;
      sigma_u4 *= 2;


      if(ThisTask == 0)
	{
	  fprintf(fd, "%g %g %g %g\n", rho0, sigma, sigmasfr, sigma_u4);
	}
    }


  if(All.ComovingIntegrationOn)
    {
      All.Time = All.TimeBegin;
      IonizeParams();
    }

  if(ThisTask == 0)
    fclose(fd);
}

#endif /* closing of SFR-conditional */


#if defined(SFR)
void set_units_sfr(void)
{
  double meanweight;

#ifdef COSMIC_RAYS
  double feedbackenergyinergs;
#endif

  All.OverDensThresh =
    All.CritOverDensity * All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

  All.PhysDensThresh = All.CritPhysDensity * PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs;

  meanweight = 4 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: assuming NEUTRAL GAS */

#ifdef CONSTANT_MEAN_MOLECULAR_WEIGHT
  meanweight = All.MeanWeight;
#endif

  All.EgySpecCold = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempClouds;
  All.EgySpecCold *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */

#ifdef CONSTANT_MEAN_MOLECULAR_WEIGHT
  meanweight = All.MeanWeight;
#endif

  All.EgySpecSN = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempSupernova;
  All.EgySpecSN *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;


#ifdef COSMIC_RAYS
  if(All.CR_SNEff < 0.0)
    /* if CR_SNeff < 0.0, then substract CR Feedback energy from thermal
     * feedback energy
     */
    {
      if(ThisTask == 0)
	{
	  printf("%g percent of thermal feedback go into Cosmic Rays.\nRemaining ", -100.0 * All.CR_SNEff);
	}

      All.EgySpecSN *= (1.0 + All.CR_SNEff);
      All.CR_SNEff = -All.CR_SNEff / (1.0 + All.CR_SNEff);

    }

  All.FeedbackEnergy = All.FactorSN / (1 - All.FactorSN) * All.EgySpecSN;

  feedbackenergyinergs = All.FeedbackEnergy / All.UnitMass_in_g * (All.UnitEnergy_in_cgs * SOLAR_MASS);

  if(ThisTask == 0)
    {
      printf("Feedback energy per formed solar mass in stars= %g  ergs\n", feedbackenergyinergs);
      printf("OverDensThresh= %g\nPhysDensThresh= %g (internal units)\n", All.OverDensThresh,
	     All.PhysDensThresh);
    }

#endif
}


#endif /* closes SFR */

#endif /* closes COOLING */
