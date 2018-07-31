#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_sf_gamma.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"
#ifdef COSMIC_RAYS
#include "cosmic_rays.h"
#endif

#ifdef MACHNUM
#ifdef COSMIC_RAYS
#define h  All.HubbleParam
#define cm (h/All.UnitLength_in_cm)
#define s  (h/All.UnitTime_in_s)
#define LightSpeed (2.9979e10*cm/s)
#define c2   ( LightSpeed * LightSpeed )
#endif
#endif



/*! \file init.c
 *  \brief code for initialisation of a simulation from initial conditions
 */


/*! This function reads the initial conditions, and allocates storage for the
 *  tree(s). Various variables of the particle data are initialised and An
 *  intial domain decomposition is performed. If SPH particles are present,
 *  the inial SPH smoothing lengths are determined.
 */
void init(void)
{
  int i, j;
  double a3, atime;

#if defined(COSMIC_RAYS) && defined(MACHNUM)
  double Pth1, PCR1, rBeta, C_phys, q_phys;
#endif
#ifdef CR_INITPRESSURE
  double cr_pressure, q_phys, C_phys;
#endif
#ifdef BLACK_HOLES
  int count_holes = 0;
#endif
#ifdef CHEMISTRY
  int ifunc;
  double min_t_cool, max_t_cool;
  double min_t_elec, max_t_elec;
  double a_start, a_end;
#endif

#if defined(DISTORTIONTENSOR) || defined(OUTPUT_TIDALTENSOR)
  int i1;
#endif

  All.Time = All.TimeBegin;

  if(RestartFlag == 3 && RestartSnapNum < 0)
    {
      if(ThisTask == 0)
	printf("Need to give the snapshot number if FOF/SUBFIND is selected for output\n");
      endrun(0);
    }

  switch (All.ICFormat)
    {
    case 1:
    case 2:
    case 3:
      if(RestartFlag >= 2 && RestartSnapNum >= 0)
	{
	  char fname[1000];

	  sprintf(fname, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, RestartSnapNum);
	  read_ic(fname);
	}
      else
	read_ic(All.InitCondFile);
      break;
    case 4:
      if(RestartFlag == 2)
         read_ic(All.InitCondFile);
      else
         read_ic_cluster(All.InitCondFile);
      break;
    default:
      if(ThisTask == 0)
	printf("ICFormat=%d not supported.\n", All.ICFormat);
      endrun(0);
    }

  All.Time = All.TimeBegin;

#ifdef CR_DIFFUSION_GREEN
  All.TimeOfLastDiffusion = All.Time;
#endif

#ifdef COOLING
  IonizeParams();
#endif

#ifdef CHEMISTRY
  InitChem();
#endif

  if(All.ComovingIntegrationOn)
    {
      All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
      All.Ti_Current = 0;
      a3 = All.Time * All.Time * All.Time;
      atime = All.Time;
    }
  else
    {
      All.Timebase_interval = (All.TimeMax - All.TimeBegin) / TIMEBASE;
      All.Ti_Current = 0;
      a3 = 1;
      atime = 1;
    }

  set_softenings();

  All.NumCurrentTiStep = 0;	/* setup some counters */
  All.SnapshotFileCount = 0;
  if(RestartFlag == 2)
    {
      if(RestartSnapNum < 0)
	All.SnapshotFileCount = atoi(All.InitCondFile + strlen(All.InitCondFile) - 3) + 1;
      else
	All.SnapshotFileCount = RestartSnapNum + 1;
    }

#ifdef OUTPUTLINEOFSIGHT
  All.Ti_nextlineofsight = (int) (log(All.TimeFirstLineOfSight / All.TimeBegin) / All.Timebase_interval);
  if(RestartFlag == 2)
    endrun(78787);
#endif

  All.TotNumOfForces = 0;
  All.NumForcesSinceLastDomainDecomp = 0;
#if defined(MAGNETIC) && defined(BSMOOTH)
  All.MainTimestepCounts = 0;
#endif

  All.TopNodeAllocFactor = 0.008;
  All.TreeAllocFactor = 0.7;

  All.TopNodeAllocFactor = 0.015;//SN
  All.TreeAllocFactor = 1.5;//SN

  All.TopNodeAllocFactor = 0.02;//RJH
  All.TreeAllocFactor = 2.0;//RJH

  All.Cadj_Cost = 1.0e-30;
  All.Cadj_Cpu = 1.0e-3;

  if(All.ComovingIntegrationOn)
    if(All.PeriodicBoundariesOn == 1)
      check_omega();

  All.TimeLastStatistics = All.TimeBegin - All.TimeBetStatistics;
#ifdef BLACK_HOLES
  All.TimeNextBlackHoleCheck = All.TimeBegin;
#endif



#ifdef BUBBLES
  if(All.ComovingIntegrationOn)
    All.TimeOfNextBubble = 1. / (1. + All.FirstBubbleRedshift);
  else
    All.TimeOfNextBubble = All.TimeBegin + All.BubbleTimeInterval / All.UnitTime_in_Megayears;
  if(ThisTask == 0)
    printf("Initial time: %g and first bubble time %g \n", All.TimeBegin, All.TimeOfNextBubble);

  if(RestartFlag == 2 && All.TimeBegin > All.TimeOfNextBubble)
    {
      printf("Restarting from the snapshot file with the wrong FirstBubbleRedshift! \n");
      endrun(0);
    }
#endif

#ifdef MULTI_BUBBLES
  if(All.ComovingIntegrationOn)
    All.TimeOfNextBubble = 1. / (1. + All.FirstBubbleRedshift);
  else
    All.TimeOfNextBubble = All.TimeBegin + All.BubbleTimeInterval / All.UnitTime_in_Megayears;
  if(ThisTask == 0)
    printf("Initial time: %g and time of the first bubbles %g \n", All.TimeBegin, All.TimeOfNextBubble);
  if(RestartFlag == 2 && All.TimeBegin > All.TimeOfNextBubble)
    {
      printf("Restarting from the snapshot file with the wrong FirstBubbleRedshift! \n");
      endrun(0);
    }
#endif

  if(All.ComovingIntegrationOn)	/*  change to new velocity variable */
    {
      for(i = 0; i < NumPart; i++)
	for(j = 0; j < 3; j++)
	  P[i].Vel[j] *= sqrt(All.Time) * All.Time;
    }

  for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
      for(j = 0; j < 3; j++)
	P[i].g.GravAccel[j] = 0;
#if defined(DISTORTIONTENSOR) || defined(OUTPUT_TIDALTENSOR)
      for(i1 = 0; i1 < 6; i1++)
	{
	  P[i].tite.tidal_tensor[i1] = 0;
	}
#endif
#ifdef DISTORTIONTENSOR
      /*distortion tensor equals unity in the beginning */
      P[i].distortion_tensor[0] = 1;
      P[i].distortion_tensor[1] = 0;
      P[i].distortion_tensor[2] = 0;
      P[i].distortion_tensor[3] = 0;
      P[i].distortion_tensor[4] = 1;
      P[i].distortion_tensor[5] = 0;
      P[i].distortion_tensor[6] = 0;
      P[i].distortion_tensor[7] = 0;
      P[i].distortion_tensor[8] = 1;

      /*set 'velocity' of distortion tensor to zero in the beginning */
      P[i].distortion_tensor_vel[0] = 0;
      P[i].distortion_tensor_vel[1] = 0;
      P[i].distortion_tensor_vel[2] = 0;
      P[i].distortion_tensor_vel[3] = 0;
      P[i].distortion_tensor_vel[4] = 0;
      P[i].distortion_tensor_vel[5] = 0;
      P[i].distortion_tensor_vel[6] = 0;
      P[i].distortion_tensor_vel[7] = 0;
      P[i].distortion_tensor_vel[8] = 0;
      if(ThisTask == 0)
	printf("Distortion tensor init...\n");
#endif
#ifdef PMGRID
      for(j = 0; j < 3; j++)
	P[i].GravPM[j] = 0;
#endif
      P[i].Ti_begstep = 0;
      P[i].Ti_current = 0;
      P[i].TimeBin = 0;
      P[i].OldAcc = 0;
      P[i].GravCost = 1;
#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY)
      P[i].p.Potential = 0;
#endif
#ifdef STELLARAGE
      if(RestartFlag == 0)
	P[i].StellarAge = 0;
#endif

#ifdef VIRTUAL
    P[i].dt1 = P[i].dt2 = 0;
#ifdef NON_GRAY_RAD_TRANSFER
    P[i].Kappa = 1.; /* initially set it to 1 in code units */
#endif
#endif /* VIRTUAL */

#ifdef METALS
      if(RestartFlag == 0)
	P[i].Metallicity = 0;
#endif

#ifdef VIRTUAL_SET_MASS
/* sets All.OriginalGasMass to an SPH particle mass. CAUTION: it is assumed
 * that sph particles have same masses here */
      if (P[i].Type == 0) {
	  if(P[i].Mass>0.0) 
	      All.OriginalGasMass = P[i].Mass;
	  else
	      All.OriginalGasMass = All.MassTable[0];	      
      }
#endif

#ifdef BLACK_HOLES
      if(P[i].Type == 5)
	{
	  count_holes++;

	  if(RestartFlag == 0)
#ifndef SMBH_MASS_INIT
	    P[i].BH_Mass = All.SeedBlackHoleMass;

#if defined(PLANET_INITIAL_MASS) 
	  	  if (P[i].Mass < 0.9*All.SMBHmass && All.BlackHoleMass > 0) P[i].Mass = All.BlackHoleMass;
#endif
#ifdef BH_KICK
	  if (All.BlackHoleMass > 0) P[i].Mass = All.BlackHoleMass;
	  P[i].Vel[2] = 0.;
	  P[i].Vel[1] = All.BlackHoleKickVel * sin(All.BlackHoleKickAngle/180.*M_PI);
	  P[i].Vel[0] = All.BlackHoleKickVel * cos(All.BlackHoleKickAngle/180.*M_PI);
#endif


#else
#ifdef BH_KICK
#ifdef BH_FIXED
	  if (All.BlackHoleKickVel > 0.) 
	  {
	      printf("\n Blackhole should have kick vel = 0 if BH_FIXED enabled. set All.BlackHoleKickVel or #BH_FIXED \n");
	  exit(1);
	  }
#endif
	  if (All.BlackHoleMass > 0) P[i].Mass = All.BlackHoleMass;
	  P[i].Vel[2] = 0.;
	  P[i].Vel[1] = All.BlackHoleKickVel * sin(All.BlackHoleKickAngle/180.*M_PI);
	  P[i].Vel[0] = All.BlackHoleKickVel * cos(All.BlackHoleKickAngle/180.*M_PI);
#endif
	    P[i].BH_Mass = P[i].Mass;
// CBP (31/10/2008) -- Addition of Accretion Disc Mass -- used in blackhole.c -- 
//                     initialise to zero. Initialise BH mass accretion rate to zero.
	    P[i].BH_Mseed = P[i].BH_Mass;
	    P[i].AccDisc_Mass = 0.0;
// End of CBP
#endif
	}
#endif
    }

#ifdef BLACK_HOLES
  MPI_Allreduce(&count_holes, &All.TotBHs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

  for(i = 0; i < TIMEBINS; i++)
    TimeBinActive[i] = 1;

  reconstruct_timebins();

#ifdef PMGRID
  All.PM_Ti_endstep = All.PM_Ti_begstep = 0;
#endif

#ifdef DUST_POWERLAW
srand((unsigned)time(NULL));
#endif

#ifdef DUST
  All.total_heating = 0.;
  for(i = 0; i < NumPart; i++)	
    {
      P[i].DustVcoll = 0;
      if (RestartFlag == 0) 
	{
#ifdef DUST_POWERLAW
	  double r;
	  double logamin;
	  double logamax;
	  /*printf("Amax %d Amin %d", All.aDustMax, All.aDustMin);*/

	  logamin = log(All.aDustMin0);
	  logamax = log(All.aDustMax0);
	  r = rand();
	  P[i].DustRadius = exp(r/RAND_MAX * (logamax-logamin) + logamin);
	  /*	  printf("DUST_POWERLAW activate!");
	  printf("Particle %d, InitialRadius %g, New Radius %g \n",i,All.InitialDustRadius,P[i].DustRadius);
	  printf("logs %g, %g, %g, %d", logamin, logamax, r, RAND_MAX);*/

#else
	  P[i].DustRadius = All.InitialDustRadius;

#endif
	}
      if (P[i].Type == 2 && RestartFlag == 0) {
	P[i].Mass *= All.InitialDustMass;
	/*printf("particle %d, Dust Mass %g, mass = %g, type = %d \n",i,All.InitialDustMass, 
	  P[i].Mass, P[i].Type);*/
      }
    }
#endif

#ifdef FLEXSTEPS
  All.PresentMinStep = TIMEBASE;
  for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
      P[i].FlexStepGrp = (int) (TIMEBASE * get_random_number(P[i].ID));
    }
#endif


  for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
      for(j = 0; j < 3; j++)
	{
	  SphP[i].VelPred[j] = P[i].Vel[j];
	  SphP[i].a.HydroAccel[j] = 0;
#ifdef RAD_ACCEL
	  SphP[i].ra.RadAccel[j] = 0;
#endif
#ifdef DUST
	  SphP[i].da.DragAccel[j] = 0;
#endif
	}

      SphP[i].e.DtEntropy = 0;

#ifdef DUST
      SphP[i].dh.DragHeating = 0;
      SphP[i].dh.dDragHeating = 0;
#endif

#ifdef CHEMISTRY
      SphP[i].Gamma = GAMMA;	/* set universal value */
      SphP[i].t_cool = 0;
      SphP[i].t_elec = 0;
#endif

#ifdef SFR_DECOUPLING
      SphP[i].DensityOld = 0;
#endif
#ifdef SFR_PROMOTION
      SphP[i].DensityAvg = 0;
      SphP[i].EntropyAvg = 0;
      SphP[i].DensPromotion = 0;
      SphP[i].TempPromotion = 0;

#endif

      if(RestartFlag == 0)
	{
#ifndef READ_HSML
	  PPP[i].Hsml = 0;
#endif
	  SphP[i].d.Density = -1;
#ifdef COOLING
	  SphP[i].Ne = 1.0;
#endif
	  SphP[i].v.DivVel = 0;
	}
#ifdef WINDS
      SphP[i].DelayTime = 0;
#endif
#ifdef SFR
      SphP[i].Sfr = 0;
#endif
#ifdef MAGNETIC
      for(j = 0; j < 3; j++)
	{
	  SphP[i].DtB[j] = 0;
	  SphP[i].BPred[j] = SphP[i].B[j];
	}
#ifdef BINISET
      SphP[i].B[0] = All.BiniX;
      SphP[i].B[1] = All.BiniY;
      SphP[i].B[2] = All.BiniZ;
      SphP[i].BPred[0] = All.BiniX;
      SphP[i].BPred[1] = All.BiniY;
      SphP[i].BPred[2] = All.BiniZ;
#endif
#ifdef TIME_DEP_MAGN_DISP
#ifdef HIGH_MAGN_DISP_START
      SphP[i].Balpha = All.ArtMagDispConst;
#else
      SphP[i].Balpha = All.ArtMagDispMin;
#endif
      SphP[i].DtBalpha = 0.0;
#endif
#ifdef DIVBCLEANING_DEDNER
      SphP[i].Phi = SphP[i].PhiPred = SphP[i].DtPhi = 0;
#endif
#endif

#ifdef TIME_DEP_ART_VISC
#ifdef HIGH_ART_VISC_START
      if(HIGH_ART_VISC_START == 0)
	SphP[i].alpha = All.ArtBulkViscConst;
      if(HIGH_ART_VISC_START > 0)
	if(P[i].Pos[0] > HIGH_ART_VISC_START)
	  SphP[i].alpha = All.ArtBulkViscConst;
	else
	  SphP[i].alpha = All.AlphaMin;
      if(HIGH_ART_VISC_START < 0)
	if(P[i].Pos[0] < -HIGH_ART_VISC_START)
	  SphP[i].alpha = All.ArtBulkViscConst;
	else
	  SphP[i].alpha = All.AlphaMin;
#else
      SphP[i].alpha = All.AlphaMin;
#endif
      SphP[i].Dtalpha = 0.0;
#endif

#if defined(BH_THERMALFEEDBACK) || defined(BH_KINETICFEEDBACK)
      SphP[i].i.Injected_BH_Energy = 0;
#ifdef DECREASE_POISSON_NOISE 
      SphP[i].i.Injected_VIRTUAL_Energy_Old = 0;
#endif
#endif

#ifdef VIRTUAL_HEATING
      SphP[i].Injected_VIRTUAL_Energy = 0;
      SphP[i].Emitted_VIRTUAL_Energy = 0;
      SphP[i].Returned_E_Fraction = 1.; 
#endif
    }


  test_id_uniqueness();



  All.NumForcesSinceLastDomainDecomp = (long long) (1 + All.TotNumPart * All.TreeDomainUpdateFrequency);

  Flag_FullStep = 1;		/* to ensure that Peano-Hilber order is done */

  TreeReconstructFlag = 1;

  domain_Decomposition();	/* do initial domain decomposition (gives equal numbers of particles) */

  set_softenings();

  if(RestartFlag == 3)
    {
#ifdef FOF
      fof_fof(RestartSnapNum);
#endif
      endrun(0);
    }

  /* will build tree */
  ngb_treebuild();

  All.Ti_Current = 0;

  setup_smoothinglengths();

  /* at this point, the entropy variable actually contains the 
   * internal energy, read in from the initial conditions file. 
   * Once the density has been computed, we can convert to entropy.
   */

  for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
      if(header.flag_entropy_instead_u == 0)
	{
	  if(ThisTask == 0 && i == 0)
	    printf("Converting u -> entropy !\n");
	  SphP[i].Entropy = GAMMA_MINUS1 * SphP[i].Entropy / pow(SphP[i].d.Density / a3, GAMMA_MINUS1);
	}
      SphP[i].e.DtEntropy = 0;

      SphP[i].v.DivVel = 0;

#ifdef MACHNUM
      SphP[i].Shock_MachNumber = 1.0;
#ifdef COSMIC_RAYS
      Pth1 = SphP[i].Entropy * pow(SphP[i].d.Density / a3, GAMMA);

#ifdef CR_IC_PHYSICAL
      C_phys = SphP[i].CR_C0;
      q_phys = SphP[i].CR_q0;
#else
      C_phys = SphP[i].CR_C0 * pow(SphP[i].d.Density, (All.CR_Alpha - 1.0) / 3.0);
      q_phys = SphP[i].CR_q0 * pow(SphP[i].d.Density, 1.0 / 3.0);
#endif

      rBeta = gsl_sf_beta((All.CR_Alpha - 2.0) * 0.5, (3.0 - All.CR_Alpha) * 0.5) *
	gsl_sf_beta_inc((All.CR_Alpha - 2.0) * 0.5, (3.0 - All.CR_Alpha) * 0.5,
			1.0 / (1.0 + q_phys * q_phys));

      PCR1 = C_phys * c2 * SphP[i].d.Density * rBeta / 6.0;
      PCR1 *= pow(atime, -3.0 * GAMMA);
      SphP[i].PreShock_XCR = PCR1 / Pth1;

      SphP[i].PreShock_PhysicalDensity = SphP[i].d.Density / a3;
      SphP[i].PreShock_PhysicalEnergy =
	SphP[i].Entropy / GAMMA_MINUS1 * pow(SphP[i].d.Density / a3, GAMMA_MINUS1);

      SphP[i].Shock_DensityJump = 1.0001;
      SphP[i].Shock_EnergyJump = 1.0;
#endif /* COSMIC_RAYS */
#endif /* MACHNUM */

#ifdef REIONIZATION
      All.not_yet_reionized = 1;
#endif

#ifdef CR_IC_PHYSICAL
      /* Scale CR variables so that values from IC file are now the
       * physical values, not the adiabatic invariants
       */

      SphP[i].CR_C0 *= pow(SphP[i].d.Density, (1.0 - All.CR_Alpha) / 3.0);
      SphP[i].CR_q0 *= pow(SphP[i].d.Density, -1.0 / 3.0);
#endif

#ifdef CR_INITPRESSURE

      cr_pressure = CR_INITPRESSURE * SphP[i].Entropy * pow(SphP[i].d.Density / a3, GAMMA);
      SphP[i].Entropy *= (1 - CR_INITPRESSURE);
      q_phys = 1.685;
      C_phys =
	cr_pressure / (SphP[i].d.Density / a3 * CR_Tab_Beta(q_phys) * (C / All.UnitVelocity_in_cm_per_s) *
		       (C / All.UnitVelocity_in_cm_per_s) / 6.0);

      SphP[i].CR_C0 = C_phys * pow(SphP[i].d.Density, (1.0 - All.CR_Alpha) / 3.0);
      SphP[i].CR_q0 = q_phys * pow(SphP[i].d.Density, -1.0 / 3.0);
#endif
    }

#ifdef REAL_EOS
/* set Temperature, Mu and Kappa of gas particle */
  for(i = 0; i < N_gas; i++) {  
    SphP[i].Mu = All.MeanWeight;
    SphP[i].Temperature = SphP[i].Entropy * pow(SphP[i].d.Density, GAMMA_MINUS1) * SphP[i].Mu * PROTONMASS;
    SphP[i].Kappa = THOMPSON / PROTONMASS;
}
#endif


#ifdef CHEMISTRY

  if(ThisTask == 0)
    {
      printf("Initial abundances: \n");
      printf("HI=%g, HII=%g, HeI=%g, HeII=%g, HeIII=%g \n",
	     SphP[1].HI, SphP[1].HII, SphP[1].HeI, SphP[1].HeII, SphP[1].HeIII);

      printf("HM=%g, H2I=%g, H2II=%g, elec=%g, %d\n",
	     SphP[1].HM, SphP[1].H2I, SphP[1].H2II, SphP[1].elec, P[1].ID);

      printf("x=%g, y=%g, z=%g, vx=%g, vy=%g, vz=%g, density=%g, entropy=%g\n",
	     P[N_gas - 1].Pos[0], P[N_gas - 1].Pos[1], P[N_gas - 1].Pos[2], P[N_gas - 1].Vel[0],
	     P[N_gas - 1].Vel[1], P[N_gas - 1].Vel[2], SphP[N_gas - 1].Density, SphP[N_gas - 1].Entropy);
    }

  /* need predict the cooling time and elec_dot here */
  min_t_cool = min_t_elec = 1.0e30;
  max_t_cool = max_t_elec = -1.0e30;

  for(i = 0; i < N_gas; i++)
    {
      a_start = All.Time;
      a_end = All.Time + 0.001;	/* 0.001 as an arbitrary value */

      ifunc = compute_abundances(0, i, a_start, a_end);


      if(fabs(SphP[i].t_cool) < min_t_cool)
	min_t_cool = fabs(SphP[i].t_cool);
      if(fabs(SphP[i].t_cool) > max_t_cool)
	max_t_cool = fabs(SphP[i].t_cool);

      if(fabs(SphP[i].t_elec) < min_t_elec)
	min_t_elec = fabs(SphP[i].t_elec);
      if(fabs(SphP[i].t_elec) > max_t_elec)
	max_t_elec = fabs(SphP[i].t_elec);

    }

  fprintf(stdout, "PE %d t_cool min= %g, max= %g in yrs \n", ThisTask, min_t_cool, max_t_cool);
  fflush(stdout);
  fprintf(stdout, "PE %d t_elec min= %g, max= %g in yrs \n", ThisTask, min_t_elec, max_t_elec);
  fflush(stdout);

#endif



}


/*! This routine computes the mass content of the box and compares it to the
 * specified value of Omega-matter.  If discrepant, the run is terminated.
 */
void check_omega(void)
{
  double mass = 0, masstot, omega;
  int i;

  for(i = 0; i < NumPart; i++)
    mass += P[i].Mass;

  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  omega =
    masstot / (All.BoxSize * All.BoxSize * All.BoxSize) / (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G));

  if(fabs(omega - All.Omega0) > 1.0e-3)
    {
      if(ThisTask == 0)
	{
	  printf("\n\nI've found something odd!\n");
	  printf
	    ("The mass content accounts only for Omega=%g,\nbut you specified Omega=%g in the parameterfile.\n",
	     omega, All.Omega0);
	  printf("\nI better stop.\n");

	  fflush(stdout);
	}
      endrun(1);
    }
}



/*! This function is used to find an initial smoothing length for each SPH
 *  particle. It guarantees that the number of neighbours will be between
 *  desired_ngb-MAXDEV and desired_ngb+MAXDEV. For simplicity, a first guess
 *  of the smoothing length is provided to the function density(), which will
 *  then iterate if needed to find the right smoothing length.
 */
void setup_smoothinglengths(void)
{
  int i, no, p;

  if(RestartFlag == 0)
    {
      for(i = 0; i < N_gas; i++)
	{
	  no = Father[i];

	  while(10 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
	    {
	      p = Nodes[no].u.d.father;

	      if(p < 0)
		break;

	      no = p;
	    }

#ifndef READ_HSML
#ifndef TWODIMS
	  PPP[i].Hsml =
	    pow(3.0 / (4 * M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3) * Nodes[no].len;
#else
	  PPP[i].Hsml =
	    pow(1.0 / (M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 2) * Nodes[no].len;
#endif
#endif
	  if(PPP[i].Hsml > 10.0 * All.SofteningTable[0])
	    PPP[i].Hsml = All.SofteningTable[0];
	}
    }

#ifdef BLACK_HOLES
  if(RestartFlag == 0 || RestartFlag == 2)
    {
      for(i = 0; i < NumPart; i++)
	if(P[i].Type == 5)
	  PPP[i].Hsml = All.SofteningTable[5];
    }
#endif

#ifdef DUST
  if(RestartFlag == 0 || RestartFlag == 2)
    {
      for(i = 0; i < NumPart; i++)
      if(P[i].Type == 2) {
	//         P[i].DustRadius = All.InitialDustRadius;
        // PPP[i].DustRadius = All.InitialDustRadius;
//         P[i].DustRadius = 0.1;
//         PPP[i].DustRadius = 0.1;
         PPP[i].Hsml = All.SofteningTable[2];

//         PPP[i].HsmlDust = All.SofteningTable[2];
//         P[i].HsmlDust = All.SofteningTable[2];
      } else {
//        P[i].DustRadius = -1000.1;
      }
    }


#endif

#ifdef VIRTUAL
  if(RestartFlag == 0 || RestartFlag == 2)
    {
      for(i = 0; i < NumPart; i++)
	if(P[i].Type == 3)
	{
	  PPP[i].Hsml = All.SofteningTable[3];
#ifdef STELLAR_AGE
	  P[i].StellarAge = All.Time;
#endif
	}
    }
#endif  



  density();

#if defined(DUST) && defined(DUST_TWO_POPULATIONS)
  /* Initialise the micro dust densities */
  for(i = 0; i < NumPart; i++)	
    {
      if (RestartFlag == 0) 
	{
	  P[i].MicroDustMass = All.InitialMicroDustZ * P[i].Mass;
	  P[i].d_MicroDustMass = 0.;
	}
    }
#endif



#ifdef DUST
  dust_density();
#endif

#if defined(MAGNETIC) && defined(BFROMROTA)
  if(RestartFlag == 0)
    {
      if(ThisTask == 0)
	printf("Converting: Vector Potential -> Bfield\n");
      rot_a();
    }
#endif

#if defined(MAGNETIC) && defined(BSMOOTH)
  if(RestartFlag == 0)
    bsmooth();
#endif
}


void test_id_uniqueness(void)
{
  int i;
  MyIDType *ids, *ids_first;

  if(ThisTask == 0)
    {
      printf("Testing ID uniqueness...\n");
      fflush(stdout);
    }

  if(NumPart == 0)
    {
      printf("need at least one particle per cpu\n");
      endrun(8);
    }

  ids = mymalloc(NumPart * sizeof(MyIDType));
  ids_first = mymalloc(NTask * sizeof(MyIDType));

  for(i = 0; i < NumPart; i++)
    ids[i] = P[i].ID;

  parallel_sort(ids, NumPart, sizeof(MyIDType), compare_IDs);

#ifndef IGNORE_NON_UNIQUE
  for(i = 1; i < NumPart; i++)
    if(ids[i] == ids[i - 1])
      {
	printf("non-unique ID=%d found on task=%d\n", (int) ids[i], ThisTask);
	printf("particle TYPE =%d found on task=%d\n", (int) P[i].Type, ThisTask);
	endrun(12);
      }
#endif

  MPI_Allgather(&ids[0], sizeof(MyIDType), MPI_BYTE, ids_first, sizeof(MyIDType), MPI_BYTE, MPI_COMM_WORLD);

#ifndef IGNORE_NON_UNIQUE
  if(ThisTask < NTask - 1)
    if(ids[NumPart - 1] == ids_first[ThisTask + 1])
      {
	printf("non-unique ID=%d found on task=%d\n", (int) ids[NumPart - 1], ThisTask);
	endrun(13);
      }
#endif

  myfree(ids_first);
  myfree(ids);

  if(ThisTask == 0)
    {
      printf("success.\n");
      fflush(stdout);
    }
}

int compare_IDs(const void *a, const void *b)
{
  if(*((MyIDType *) a) < *((MyIDType *) b))
    return -1;

  if(*((MyIDType *) a) > *((MyIDType *) b))
    return +1;

  return 0;
}
