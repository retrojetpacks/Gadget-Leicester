#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>

#include "allvars.h"
#include "proto.h"
#ifdef COSMIC_RAYS
#include "cosmic_rays.h"
#endif


/*! \file begrun.c
 *  \brief initial set-up of a simulation run
 *
 *  This file contains various functions to initialize a simulation run. In
 *  particular, the parameterfile is read in and parsed, the initial
 *  conditions or restart files are read, and global variables are initialized
 *  to their proper values.
 */



/*! This function performs the initial set-up of the simulation. First, the
 *  parameterfile is set, then routines for setting units, reading
 *  ICs/restart-files are called, auxialiary memory is allocated, etc.
 */
void begrun(void)
{
  struct global_data_all_processes all;

  if(ThisTask == 0)
    {
      /*    printf("\nThis is P-Gadget, version `%s', svn-revision `%s'.\n", GADGETVERSION, svn_version()); */
      printf("\nThis is P-Gadget, version `%s'.\n", GADGETVERSION);
      printf("\nRunning on %d processors.\n", NTask);
    }

#if defined (NFW_POTENTIAL) && defined (CUSP_POTENTIAL)
  if(ThisTask==0) 
      printf("Cannot run with both NFW and Cusp potentials....\n");
  
  endrun(1);
#endif

#if defined(X86FIX) && defined(SOFTDOUBLEDOUBLE)
  x86_fix();			/* disable 80bit treatment of internal FPU registers in favour of proper IEEE 64bit double precision arithmetic */
#endif

#ifdef PEDANTIC_MEMORY_HANDLER
  mymalloc_init();
#endif

  read_parameter_file(ParameterFile);	/* ... read in parameters for this run */

#ifdef DEBUG
  write_pid_file();
  enable_core_dumps_and_fpu_exceptions();
#endif

  set_units();

#ifdef REGULATED_FB  
  // CBP (03/11/2008) -- Echo input parameters to screen
  if(ThisTask==0)
  {
      printf("Regulating black hole feedback!\n\n");
      printf("Viscous Timescale (Myrs) : %lf\n",All.ViscousTimescale);
      printf("Viscous Timescale (Codes Units) : %lf\n",All.ViscousTimescale/All.UnitTime_in_Megayears);
  }
#endif

/* Addition by CBP. 27/07/2009. */
#ifdef NFW_POTENTIAL
  All.VirialRadius = pow(2.*All.VirialMass*All.G/(All.DeltaVir*All.Hubble*All.Hubble),1./3.);
  All.ScaleRadius = All.VirialRadius/All.CNFW;

  if(ThisTask==0)
  {
      printf("Assuming a background NFW potential.\n");
      printf("Virial Mass & Radius (Code Units) : %g\t%g\n",All.VirialMass,All.VirialRadius);
      printf("Scale Radius (Code Units): %g\n",All.ScaleRadius);
  }
#endif
/* End of Addition by CBP. 27/07/2009. */

#ifdef COOLING
  All.Time = All.TimeBegin;
  InitCool();
#endif

#ifdef CHEMISTRY
  InitChem();
#endif

#if defined(SFR)
  init_clouds();
#endif

#ifdef PERIODIC
  ewald_init();
#endif

#ifdef PERIODIC
  boxSize = All.BoxSize;
  boxHalf = 0.5 * All.BoxSize;
#ifdef LONG_X
  boxHalf_X = boxHalf * LONG_X;
  boxSize_X = boxSize * LONG_X;
#endif
#ifdef LONG_Y
  boxHalf_Y = boxHalf * LONG_Y;
  boxSize_Y = boxSize * LONG_Y;
#endif
#ifdef LONG_Z
  boxHalf_Z = boxHalf * LONG_Z;
  boxSize_Z = boxSize * LONG_Z;
#endif
#endif

#ifdef DARKENERGY
#ifdef TIMEDEPDE
  fwa_init();
#endif
#endif

#ifdef TIME_DEP_ART_VISC
  All.ViscSource = All.ViscSource0 / log((GAMMA + 1) / (GAMMA - 1));
  All.DecayTime = 1 / All.DecayLength * sqrt((GAMMA - 1) / 2 * GAMMA);
#endif

  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  gsl_rng_set(random_generator, 42);	/* start-up seed */

#ifdef PMGRID
  if(RestartFlag != 3)
    long_range_init();
#endif

#ifdef SUBFIND
  GrNr = -1;
#endif

  All.TimeLastRestartFile = CPUThisRun;

  if(RestartFlag == 0 || RestartFlag == 2 || RestartFlag == 3)
    {
      set_random_numbers();

      init();			/* ... read in initial model */
    }
  else
    {
      all = All;		/* save global variables. (will be read from restart file) */

      restart(RestartFlag);	/* ... read restart file. Note: This also resets 
				   all variables in the struct `All'. 
				   However, during the run, some variables in the parameter
				   file are allowed to be changed, if desired. These need to 
				   copied in the way below.
				   Note:  All.PartAllocFactor is treated in restart() separately.  
				 */
      All.MinSizeTimestep = all.MinSizeTimestep;
      All.MaxSizeTimestep = all.MaxSizeTimestep;
      All.BufferSize = all.BufferSize;
      All.TimeLimitCPU = all.TimeLimitCPU;
      All.ResubmitOn = all.ResubmitOn;
      All.TimeBetSnapshot = all.TimeBetSnapshot;
      All.TimeBetStatistics = all.TimeBetStatistics;
      All.CpuTimeBetRestartFile = all.CpuTimeBetRestartFile;
      All.ErrTolIntAccuracy = all.ErrTolIntAccuracy;
      All.MaxRMSDisplacementFac = all.MaxRMSDisplacementFac;

      All.ErrTolForceAcc = all.ErrTolForceAcc;
      All.TypeOfTimestepCriterion = all.TypeOfTimestepCriterion;
      All.TypeOfOpeningCriterion = all.TypeOfOpeningCriterion;
      All.NumFilesWrittenInParallel = all.NumFilesWrittenInParallel;
      All.TreeDomainUpdateFrequency = all.TreeDomainUpdateFrequency;

      All.OutputListOn = all.OutputListOn;
      All.CourantFac = all.CourantFac;

      All.OutputListLength = all.OutputListLength;
      memcpy(All.OutputListTimes, all.OutputListTimes, sizeof(double) * All.OutputListLength);

#ifdef TIME_DEP_ART_VISC
      All.ViscSource = all.ViscSource;
      All.ViscSource0 = all.ViscSource0;
      All.DecayTime = all.DecayTime;
      All.DecayLength = all.DecayLength;
      All.AlphaMin = all.AlphaMin;
#endif

#ifdef MAGNETIC_DISSIPATION
      All.ArtMagDispConst = all.ArtMagDispConst;
#ifdef TIME_DEP_MAGN_DISP
      All.ArtMagDispMin = all.ArtMagDispMin;
      All.ArtMagDispSource = all.ArtMagDispSource;
      All.ArtMagDispTime = all.ArtMagDispTime;
#endif
#endif

#ifdef DIVBCLEANING_DEDNER
      All.DivBcleanParabolicSigma = all.DivBcleanParabolicSigma;
      All.DivBcleanHyperbolicSigma = all.DivBcleanHyperbolicSigma;
#endif

#ifdef DARKENERGY
      All.DarkEnergyParam = all.DarkEnergyParam;
#endif

      /* SN modifications
#ifdef VIRTUAL
    All.VirtualStart = all.VirtualStart;
    All.FeedBackVelocity= all.FeedBackVelocity;
    All.VirtualMomentum= all.VirtualMomentum;
    All.VirtualMass= all.VirtualMass;
    All.VirtualTime= all.VirtualTime;
    All.VirtualCrosSection= all.VirtualCrosSection;
    All.VirtualFeedBack= all.VirtualFeedBack;
    All.OuterBoundary= all.OuterBoundary;
    All.FreeTravelLength= all.FreeTravelLength;
#endif

#if defined(BACKGROUND_ILLUMINATION_TEMP) || defined(BETA_COOLING) || defined(ISOTHERM)
    All.EqTemp = all.EqTemp; 
#endif

#ifdef VIRTUAL_SET_MASS
    All.OriginalGasMass = all.OriginalGasMass;
#endif

#ifdef DUST
  All.InitialDustRadius = all.InitialDustRadius;
  All.InitialDustMass = all.InitialDustMass;
#endif
      */

      strcpy(All.ResubmitCommand, all.ResubmitCommand);
      strcpy(All.OutputListFilename, all.OutputListFilename);
      strcpy(All.OutputDir, all.OutputDir);
      strcpy(All.RestartFile, all.RestartFile);
      strcpy(All.EnergyFile, all.EnergyFile);
      strcpy(All.InfoFile, all.InfoFile);
      strcpy(All.CpuFile, all.CpuFile);
      strcpy(All.TimingsFile, all.TimingsFile);
      strcpy(All.SnapshotFileBase, all.SnapshotFileBase);

      if(All.TimeMax != all.TimeMax)
	readjust_timebase(All.TimeMax, all.TimeMax);

#ifdef NO_TREEDATA_IN_RESTART
      /* if this is not activated, the tree was stored in the restart-files,
         which also allocated the storage for it already */

      /* ensures that domain reconstruction will be done and new tree will be constructed */
      All.NumForcesSinceLastDomainDecomp = (long long) (1 + All.TotNumPart * All.TreeDomainUpdateFrequency);
#endif
    }

  open_outputfiles();

#ifdef PMGRID
  long_range_init_regionsize();
#endif

  reconstruct_timebins();

#ifdef COSMIC_RAYS
  CR_initialize_beta_tabs(All.CR_Alpha);
  CR_Tab_Initialize();
#ifdef COSMIC_RAY_TEST
  CR_test_routine();
#endif

#endif

// Addition by CBP (10/10/2008)
//
// All.SofteningBulgeMaxPhys should define an upper limit for the 
// photon search radius. 
#ifdef VARIABLE_PHOTON_SEARCH_RADIUS
  if(ThisTask==0) {
      printf("\nAllowing feedback particles to have variable search radius\n");
      printf("\n PhotonSearchRadius : %lf %lf\n",All.PhotonSearchRadius,All.SofteningBulgeMaxPhys);
      fflush(stdout);
  }

#endif
// End of CBP

#if defined(VIRTUAL) //&& defined(FB_MOMENTUM)
// Addition by CBP (10/11/2008)
//
// Set VirtualMomentum=1 in input file, 
  if(ThisTask==0) {
      printf("Resetting All.VirtualMomentum (%g) to All.OriginalGasMass (%g) : %g \n",All.VirtualMomentum,All.OriginalGasMass,All.VirtualMomentum*All.OriginalGasMass);
      fflush(stdout);
  }
  
  if(All.OriginalGasMass>0.0)
      All.VirtualMomentum *= All.OriginalGasMass;

// End of CBP
#endif

  if(All.ComovingIntegrationOn)
    init_drift_table();

  if(RestartFlag == 2)
    All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 100);
  else
    All.Ti_nextoutput = find_next_outputtime(All.Ti_Current);


  All.TimeLastRestartFile = CPUThisRun;
}




/*! Computes conversion factors between internal code units and the
 *  cgs-system.
 */
void set_units(void)
{
  double meanweight;

#ifdef CONDUCTION
#ifndef CONDUCTION_CONSTANT
  double coulomb_log;
#endif
#endif
#ifdef STATICNFW
  double Mtot;
#endif

  All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
  All.UnitTime_in_Megayears = All.UnitTime_in_s / SEC_PER_MEGAYEAR;

  if(All.GravityConstantInternal == 0)
    All.G = GRAVITY / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g * pow(All.UnitTime_in_s, 2);
  else
    All.G = All.GravityConstantInternal;
#ifdef TIMEDEPGRAV
  All.Gini=All.G;
  All.G = All.Gini * dGfak(All.TimeBegin);
#endif

  All.UnitDensity_in_cgs = All.UnitMass_in_g / pow(All.UnitLength_in_cm, 3);
  All.UnitPressure_in_cgs = All.UnitMass_in_g / All.UnitLength_in_cm / pow(All.UnitTime_in_s, 2);
  All.UnitCoolingRate_in_cgs = All.UnitPressure_in_cgs / All.UnitTime_in_s;
  All.UnitEnergy_in_cgs = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);

  /* convert some physical input parameters to internal units */

  All.Hubble = HUBBLE * All.UnitTime_in_s;

  if(ThisTask == 0)
    {
      printf("\nHubble (internal units) = %g\n", All.Hubble);
      printf("G (internal units) = %g\n", All.G);
      printf("UnitMass_in_g = %g \n", All.UnitMass_in_g);
      printf("UnitTime_in_s = %g \n", All.UnitTime_in_s);
      printf("UnitVelocity_in_cm_per_s = %g \n", All.UnitVelocity_in_cm_per_s);
      printf("UnitDensity_in_cgs = %g \n", All.UnitDensity_in_cgs);
      printf("UnitEnergy_in_cgs = %g \n", All.UnitEnergy_in_cgs);
      printf("\n");
    }

  meanweight = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: assuming NEUTRAL GAS */

#ifdef CONSTANT_MEAN_MOLECULAR_WEIGHT
  meanweight = All.MeanWeight;
#endif

  All.MinEgySpec = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.MinGasTemp;
  All.MinEgySpec *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

#if defined(SFR)
  set_units_sfr();
#endif


#define cm (All.HubbleParam/All.UnitLength_in_cm)
#define g  (All.HubbleParam/All.UnitMass_in_g)
#define s  (All.HubbleParam/All.UnitTime_in_s)
#define erg (g*cm*cm/(s*s))
#define keV (1.602e-9*erg)
#define deg 1.0
#define m_p (PROTONMASS * g)
#define k_B (BOLTZMANN * erg / deg)

#ifdef NAVIERSTOKES
  /* Braginskii-Spitzer shear viscosity parametrization */
  /* mu = 0.406 * m_p^0.5 * (k_b* T)^(5/2) / e^4 / logLambda  [g/cm/s] */
  /* eta = frac * mu */

  meanweight = 4.0 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* assuming full ionization */

#ifdef CONSTANT_MEAN_MOLECULAR_WEIGHT
  meanweight = All.MeanWeight;
#endif

#ifdef NAVIERSTOKES_CONSTANT
  All.NavierStokes_ShearViscosity = All.FractionSpitzerViscosity * 0.406 * pow(PROTONMASS, 0.5) * pow(BOLTZMANN * All.ShearViscosityTemperature, 5. / 2.) / pow(ELECTRONCHARGE, 4) / LOG_LAMBDA;	/*in cgs units */

  if(ThisTask == 0)
    printf("Constant shear viscosity in cgs units: eta = %g\n", All.NavierStokes_ShearViscosity);

  All.NavierStokes_ShearViscosity *= All.UnitTime_in_s * All.UnitLength_in_cm / All.UnitMass_in_g / All.HubbleParam;	/* in internal code units */

  if(ThisTask == 0)
    printf("Constant shear viscosity in internal code units: eta = %g\n", All.NavierStokes_ShearViscosity);

#else
  All.NavierStokes_ShearViscosity = All.FractionSpitzerViscosity * 0.406 * pow(PROTONMASS, 0.5) * pow((meanweight * PROTONMASS * GAMMA_MINUS1), 5. / 2.) / pow(ELECTRONCHARGE, 4) / LOG_LAMBDA;	/*in cgs units */
  /*T = mu*m_p*(gamma-1)/k_b * E * UnitEnergy/UnitMass */

  All.NavierStokes_ShearViscosity *= pow((All.UnitEnergy_in_cgs / All.UnitMass_in_g), 5. / 2.);	/* now energy can be multiplied later in the internal code units */
  All.NavierStokes_ShearViscosity *= All.UnitTime_in_s * All.UnitLength_in_cm / All.UnitMass_in_g / All.HubbleParam;	/* in internal code units */

  if(ThisTask == 0)
    printf("Variable shear viscosity in internal code units: eta = %g\n", All.NavierStokes_ShearViscosity);

#endif

#ifdef NAVIERSTOKES_BULK
  if(ThisTask == 0)
    printf("Costant bulk viscosity in internal code units: zeta = %g\n", All.NavierStokes_BulkViscosity);
#endif

#ifdef VISCOSITY_SATURATION
  /* calculate ion mean free path assuming complete ionization: 
     ion mean free path for hydrogen is similar to that of helium, 
     thus we calculate only for hydrogen */
  /* l_i = 3^(3/2)*(k*T)^2 / (4*\pi^(1/2)*ni*(Z*e)^4*lnL) */

  All.IonMeanFreePath = pow(3.0, 1.5) / (4.0 * sqrt(M_PI) * pow(ELECTRONCHARGE, 4) * LOG_LAMBDA);

  All.IonMeanFreePath *= pow(meanweight * PROTONMASS * GAMMA_MINUS1, 2) * pow((All.UnitEnergy_in_cgs / All.UnitMass_in_g), 2);	/*kT -> u */

  All.IonMeanFreePath /= (HYDROGEN_MASSFRAC / PROTONMASS) *
    (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam);
  /* n_H = rho * Hfr / mp *//* now is cgs units *///changed / to * in front of the unitdensity

  All.IonMeanFreePath *= All.HubbleParam / All.UnitLength_in_cm;
  /* in internal code units */
#endif

#endif

#ifdef CONDUCTION
#ifndef CONDUCTION_CONSTANT

  meanweight = m_p * 4.0 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));
#ifdef CONSTANT_MEAN_MOLECULAR_WEIGHT
  meanweight = All.MeanWeight;
#endif

  /* assuming full ionization */

  coulomb_log = 37.8;
  /* accordin1g to Sarazin's book */

  All.ConductionCoeff *=
    (1.84e-5 / coulomb_log * pow(meanweight / k_B * GAMMA_MINUS1, 2.5) * erg / (s * deg * cm));
  /* Kappa_Spitzer definition taken from Zakamska & Narayan 2003 
   * ( ApJ 582:162-169, Eq. (5) )
   */

  /* Note: Because we replace \nabla(T) in the conduction equation with
   * \nable(u), our conduction coefficient is not the usual kappa, but
   * rather kappa*(gamma-1)*mu/kB. We therefore need to multiply with 
   * another factor of (meanweight / k_B * GAMMA_MINUS1).
   */
  All.ConductionCoeff *= meanweight / k_B * GAMMA_MINUS1;

  /* The conversion of  ConductionCoeff between internal units and cgs
   * units involves one factor of 'h'. We take care of this here.
   */
  All.ConductionCoeff /= All.HubbleParam;

#ifdef CONDUCTION_SATURATION
  All.ElectronFreePathFactor = 8 * pow(3.0, 1.5) * pow(GAMMA_MINUS1, 2) / pow(3 + 5 * HYDROGEN_MASSFRAC, 2)
    / (1 + HYDROGEN_MASSFRAC) / sqrt(M_PI) / coulomb_log * pow(PROTONMASS, 3) / pow(ELECTRONCHARGE, 4)
    / (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam)
    * pow(All.UnitPressure_in_cgs / All.UnitDensity_in_cgs, 2);

  /* If the above value is multiplied with u^2/rho in code units (with rho being the physical density), then
   * one gets the electrong mean free path in centimeter. Since we want to compare this with another length
   * scale in code units, we now add an additional factor to convert back to code units.
   */
  All.ElectronFreePathFactor *= All.HubbleParam / All.UnitLength_in_cm;
#endif

#endif /* CONDUCTION_CONSTANT */
#endif /* CONDUCTION */

#if defined(CR_DIFFUSION) || defined(CR_DIFFUSION_GREEN)
  if(All.CR_DiffusionDensZero == 0.0)
    {
      /* Set reference density for CR Diffusion to rhocrit at z=0 */

      All.CR_DiffusionDensZero = 3.0 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
    }

  if(All.CR_DiffusionEntropyZero == 0.0)
    {
      All.CR_DiffusionEntropyZero = 1.0e4;
    }

  /* Set reference entropic function to correspond to 
     Reference Temperature @ ReferenceDensity */
  if(ThisTask == 0)
    {
      printf("CR Diffusion: T0 = %g\n", All.CR_DiffusionEntropyZero);
    }

  /* convert Temperature value in Kelvin to thermal energy per unit mass
     in internal units, and then to entropy */
  All.CR_DiffusionEntropyZero *=
    BOLTZMANN / (4.0 * PROTONMASS / (3.0 * HYDROGEN_MASSFRAC + 1.0)) *
    All.UnitMass_in_g / All.UnitEnergy_in_cgs / pow(All.CR_DiffusionDensZero, GAMMA_MINUS1);

  /* Change the density scaling, so that the temp scaling is mapped
     onto an entropy scaling that is numerically less expensive */
  All.CR_DiffusionDensScaling += GAMMA_MINUS1 * All.CR_DiffusionEntropyScaling;

  if(ThisTask == 0)
    {
      printf("CR Diffusion: Rho0 = %g -- A0 = %g\n", All.CR_DiffusionDensZero, All.CR_DiffusionEntropyZero);
    }

  if(All.CR_DiffusionTimeScale <= 0.0)
    {
      All.CR_DiffusionTimeScale = 0.1;
    }

#endif /* CR_DIFFUSION */

#ifdef STATICNFW
  R200 = pow(NFW_M200 * All.G / (100 * All.Hubble * All.Hubble), 1.0 / 3);
  Rs = R200 / NFW_C;
  Dc = 200.0 / 3 * NFW_C * NFW_C * NFW_C / (log(1 + NFW_C) - NFW_C / (1 + NFW_C));
  RhoCrit = 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
  V200 = 10 * All.Hubble * R200;
  if(ThisTask == 0)
    printf("V200= %g\n", V200);

  fac = 1.0;
  Mtot = enclosed_mass(R200);
  if(ThisTask == 0)
    printf("M200= %g\n", Mtot);
  fac = V200 * V200 * V200 / (10 * All.G * All.Hubble) / Mtot;
  Mtot = enclosed_mass(R200);
  if(ThisTask == 0)
    printf("M200= %g\n", Mtot);
#endif
}

#ifdef STATICNFW
/*! auxiliary function for static NFW potential
 */
double enclosed_mass(double R)
{
  /* Eps is in units of Rs !!!! */

  if(R > Rs * NFW_C)
    R = Rs * NFW_C;

  return fac * 4 * M_PI * RhoCrit * Dc *
    (-(Rs * Rs * Rs * (1 - NFW_Eps + log(Rs) - 2 * NFW_Eps * log(Rs) + NFW_Eps * NFW_Eps * log(NFW_Eps * Rs)))
     / ((NFW_Eps - 1) * (NFW_Eps - 1)) +
     (Rs * Rs * Rs *
      (Rs - NFW_Eps * Rs - (2 * NFW_Eps - 1) * (R + Rs) * log(R + Rs) +
       NFW_Eps * NFW_Eps * (R + Rs) * log(R + NFW_Eps * Rs))) / ((NFW_Eps - 1) * (NFW_Eps - 1) * (R + Rs)));
}
#endif



/*!  This function opens various log-files that report on the status and
 *   performance of the simulstion. On restart from restart-files
 *   (start-option 1), the code will append to these files.
 */
void open_outputfiles(void)
{
  char mode[2], buf[200];

  if(RestartFlag == 0)
    strcpy(mode, "w");
  else
    strcpy(mode, "a");

#ifdef BLACK_HOLES
  /* Note: This is done by everyone */
  sprintf(buf, "%sblackhole_details_%d.txt", All.OutputDir, ThisTask);
  if(!(FdBlackHolesDetails = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef REAL_EOS
    int i,j;
    double temp_shc;
    sprintf(buf, "KAPPA.DATA");
    data_kappa = fopen(buf, "r");
/*
    for (i=0;i<=2850;i++) {
        fscanf(data_kappa,"%le %le",&temp_shc,&kappa1[i]);
        fprintf(diag_shc,"%d %le\n",i,kappa1[i]);
    }
    for (i=0;i<=81;i++) {
    for (j=0;j<=1251;j++) {
        fscanf(data_kappa,"%le %le %le",&temp_shc,&temp_shc,&kappa2[i][j]);
        fprintf(diag_shc,"%d %d %le\n",i,j,kappa2[i][j]);
    }
    }
    fclose(data_kappa);
    fflush (diag_shc);
*/
/*  For Plank Opacity Krumholz et al 2007 
    for (i=0;i<=375;i++)    kappa_P[i] = 0.3+7.0*(double)(i/375);
    for (i=376;i<=575;i++)  kappa_P[i] = 7.3+0.7*(double)(i-375)/200.;
    for (i=576;i<=675;i++)  kappa_P[i] = 3.0+0.1*(double)(i-575)/100.;
    for (i=676;i<=960;i++)  kappa_P[i] = 2.8+0.3*(double)(i-675)/285.;
    for (i=961;i<=1100;i++) kappa_P[i] = 3.1-0.3*(double)(i-960)/140.;
*/

/*  For Plank Opacity Krumholz et al 2007 
    for (i=0;i<=350;i++)    kappa_R[i] = 0.1+4.4*(double)(i/350);
    for (i=351;i<=600;i++)  kappa_R[i] = 3.9;
    for (i=601;i<=700;i++)  kappa_R[i] = 0.7;
    for (i=701;i<=950;i++)  kappa_R[i] = 0.25; 
    for (i=951;i<=1000;i++) kappa_R[i] = 0.25-0.15*(double)(i-950)/50.;
*/
    
#endif
/* SHC */

  if(ThisTask != 0)		/* only the root processors writes to the log files */
    return;

  sprintf(buf, "%s%s", All.OutputDir, All.CpuFile);
  if(!(FdCPU = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.InfoFile);
  if(!(FdInfo = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.EnergyFile);
  if(!(FdEnergy = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.TimingsFile);
  if(!(FdTimings = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, "balance.txt");
  if(!(FdBalance = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  fprintf(FdBalance, "\n");
  fprintf(FdBalance, "Treewalk1      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWALK1],
	  CPU_SymbolImbalance[CPU_TREEWALK1]);
  fprintf(FdBalance, "Treewalk2      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWALK2],
	  CPU_SymbolImbalance[CPU_TREEWALK2]);
  fprintf(FdBalance, "Treewait1      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWAIT1],
	  CPU_SymbolImbalance[CPU_TREEWAIT1]);
  fprintf(FdBalance, "Treewait2      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWAIT2],
	  CPU_SymbolImbalance[CPU_TREEWAIT2]);
  fprintf(FdBalance, "Treesend       = '%c' / '%c'\n", CPU_Symbol[CPU_TREESEND],
	  CPU_SymbolImbalance[CPU_TREESEND]);
  fprintf(FdBalance, "Treerecv       = '%c' / '%c'\n", CPU_Symbol[CPU_TREERECV],
	  CPU_SymbolImbalance[CPU_TREERECV]);
  fprintf(FdBalance, "Treebuild      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEBUILD],
	  CPU_SymbolImbalance[CPU_TREEBUILD]);
  fprintf(FdBalance, "Treeupdate     = '%c' / '%c'\n", CPU_Symbol[CPU_TREEUPDATE],
	  CPU_SymbolImbalance[CPU_TREEUPDATE]);
  fprintf(FdBalance, "Treehmaxupdate = '%c' / '%c'\n", CPU_Symbol[CPU_TREEHMAXUPDATE],
	  CPU_SymbolImbalance[CPU_TREEHMAXUPDATE]);
  fprintf(FdBalance, "Treemisc =       '%c' / '%c'\n", CPU_Symbol[CPU_TREEMISC],
	  CPU_SymbolImbalance[CPU_TREEMISC]);
  fprintf(FdBalance, "Domain decomp  = '%c' / '%c'\n", CPU_Symbol[CPU_DOMAIN],
	  CPU_SymbolImbalance[CPU_DOMAIN]);
  fprintf(FdBalance, "Density compute= '%c' / '%c'\n", CPU_Symbol[CPU_DENSCOMPUTE],
	  CPU_SymbolImbalance[CPU_DENSCOMPUTE]);
  fprintf(FdBalance, "Density imbal  = '%c' / '%c'\n", CPU_Symbol[CPU_DENSWAIT],
	  CPU_SymbolImbalance[CPU_DENSWAIT]);
  fprintf(FdBalance, "Density commu  = '%c' / '%c'\n", CPU_Symbol[CPU_DENSCOMM],
	  CPU_SymbolImbalance[CPU_DENSCOMM]);
  fprintf(FdBalance, "Density misc   = '%c' / '%c'\n", CPU_Symbol[CPU_DENSMISC],
	  CPU_SymbolImbalance[CPU_DENSMISC]);
  fprintf(FdBalance, "Hydro compute  = '%c' / '%c'\n", CPU_Symbol[CPU_HYDCOMPUTE],
	  CPU_SymbolImbalance[CPU_HYDCOMPUTE]);
  fprintf(FdBalance, "Hydro imbalance= '%c' / '%c'\n", CPU_Symbol[CPU_HYDWAIT],
	  CPU_SymbolImbalance[CPU_HYDWAIT]);
  fprintf(FdBalance, "Hydro comm     = '%c' / '%c'\n", CPU_Symbol[CPU_HYDCOMM],
	  CPU_SymbolImbalance[CPU_HYDCOMM]);
  fprintf(FdBalance, "Hydro misc     = '%c' / '%c'\n", CPU_Symbol[CPU_HYDMISC],
	  CPU_SymbolImbalance[CPU_HYDMISC]);
  fprintf(FdBalance, "Drifts         = '%c' / '%c'\n", CPU_Symbol[CPU_DRIFT], CPU_SymbolImbalance[CPU_DRIFT]);
  fprintf(FdBalance, "Blackhole      = '%c' / '%c'\n", CPU_Symbol[CPU_BLACKHOLE],
	  CPU_SymbolImbalance[CPU_BLACKHOLE]);
  fprintf(FdBalance, "Kicks          = '%c' / '%c'\n", CPU_Symbol[CPU_TIMELINE],
	  CPU_SymbolImbalance[CPU_TIMELINE]);
  fprintf(FdBalance, "Potential      = '%c' / '%c'\n", CPU_Symbol[CPU_POTENTIAL],
	  CPU_SymbolImbalance[CPU_POTENTIAL]);
  fprintf(FdBalance, "PM             = '%c' / '%c'\n", CPU_Symbol[CPU_MESH], CPU_SymbolImbalance[CPU_MESH]);
  fprintf(FdBalance, "Peano-Hilbert  = '%c' / '%c'\n", CPU_Symbol[CPU_PEANO], CPU_SymbolImbalance[CPU_PEANO]);
  fprintf(FdBalance, "Cooling & SFR  = '%c' / '%c'\n", CPU_Symbol[CPU_COOLINGSFR],
	  CPU_SymbolImbalance[CPU_COOLINGSFR]);
  fprintf(FdBalance, "Snapshot dump  = '%c' / '%c'\n", CPU_Symbol[CPU_SNAPSHOT],
	  CPU_SymbolImbalance[CPU_SNAPSHOT]);
  fprintf(FdBalance, "FoF            = '%c' / '%c'\n", CPU_Symbol[CPU_FOF], CPU_SymbolImbalance[CPU_FOF]);
  fprintf(FdBalance, "Miscellaneous  = '%c' / '%c'\n", CPU_Symbol[CPU_MISC], CPU_SymbolImbalance[CPU_MISC]);
  fprintf(FdBalance, "\n");

#ifdef SFR
  sprintf(buf, "%s%s", All.OutputDir, "sfr.txt");
  if(!(FdSfr = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif


#ifdef BLACK_HOLES
  sprintf(buf, "%s%s", All.OutputDir, "blackholes.txt");
  if(!(FdBlackHoles = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif



#ifdef FORCETEST
  if(RestartFlag == 0)
    {
      sprintf(buf, "%s%s", All.OutputDir, "forcetest.txt");
      if(!(FdForceTest = fopen(buf, "w")))
	{
	  printf("error in opening file '%s'\n", buf);
	  endrun(1);
	}
      fclose(FdForceTest);
    }
#endif

#ifdef XXLINFO
  sprintf(buf, "%s%s", All.OutputDir, "xxl.txt");
  if(!(FdXXL = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  else
    {
      if(RestartFlag == 0)
	{
	  fprintf(FdXXL, "nstep time ");
#ifdef MAGNETIC
	  fprintf(FdXXL, "<|B|> ");
#ifdef TRACEDIVB
	  fprintf(FdXXL, "max(divB) ");
#endif
#endif
#ifdef TIME_DEP_ART_VISC
	  fprintf(FdXXL, "<alpha> ");
#endif
	  fprintf(FdXXL, "\n");
	  fflush(FdXXL);
	}
    }
#endif

#ifdef DARKENERGY
  sprintf(buf, "%s%s", All.OutputDir, "darkenergy.txt");
  if(!(FdDE = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  else
    {
      if(RestartFlag == 0)
	{
	  fprintf(FdDE, "nstep time H(a) ");
#ifndef TIMEDEPDE
	  fprintf(FdDE, "w0 Omega_L ");
#else
	  fprintf(FdDE, "w(a) Omega_L ");
#endif
#ifdef TIMEDEPGRAV
	  fprintf(FdDE, "dH dG ");
#endif
	  fprintf(FdDE, "\n");
	  fflush(FdDE);
	}
    }
#endif

}




/*!  This function closes the global log-files.
 */
void close_outputfiles(void)
{
#ifdef BLACK_HOLES
  fclose(FdBlackHolesDetails);	/* needs to be done by everyone */
#endif

  if(ThisTask != 0)		/* only the root processors writes to the log files */
    return;

  fclose(FdCPU);
  fclose(FdInfo);
  fclose(FdEnergy);
  fclose(FdTimings);
  fclose(FdBalance);

#ifdef SFR
  fclose(FdSfr);
#endif

#ifdef BLACK_HOLES
  fclose(FdBlackHoles);
#endif

#ifdef XXLINFO
  fclose(FdXXL);
#endif

#ifdef DARKENERGY
  fclose(FdDE); 
#endif
}





/*! This function parses the parameterfile in a simple way.  Each paramater is
 *  defined by a keyword (`tag'), and can be either of type douple, int, or
 *  character string.  The routine makes sure that each parameter appears
 *  exactly once in the parameterfile, otherwise error messages are
 *  produced that complain about the missing parameters.
 */
void read_parameter_file(char *fname)
{
#define REAL 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd, *fdout;
  char buf[200], buf1[200], buf2[200], buf3[400];
  int i, j, nt;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int pnum, errorFlag = 0;

  All.StarformationOn = 0;	/* defaults */


  if(sizeof(long long) != 8)
    {
      if(ThisTask == 0)
	printf("\nType `long long' is not 64 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(int) != 4)
    {
      if(ThisTask == 0)
	printf("\nType `int' is not 32 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(float) != 4)
    {
      if(ThisTask == 0)
	printf("\nType `float' is not 32 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(double) != 8)
    {
      if(ThisTask == 0)
	printf("\nType `double' is not 64 bit on this platform. Stopping.\n\n");
      endrun(0);
    }


  if(ThisTask == 0)		/* read parameter file on process 0 */
    {
      nt = 0;

      strcpy(tag[nt], "InitCondFile");
      addr[nt] = All.InitCondFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputDir");
      addr[nt] = All.OutputDir;
      id[nt++] = STRING;

      strcpy(tag[nt], "SnapshotFileBase");
      addr[nt] = All.SnapshotFileBase;
      id[nt++] = STRING;

      strcpy(tag[nt], "EnergyFile");
      addr[nt] = All.EnergyFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "CpuFile");
      addr[nt] = All.CpuFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "InfoFile");
      addr[nt] = All.InfoFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "TimingsFile");
      addr[nt] = All.TimingsFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "RestartFile");
      addr[nt] = All.RestartFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "ResubmitCommand");
      addr[nt] = All.ResubmitCommand;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListFilename");
      addr[nt] = All.OutputListFilename;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListOn");
      addr[nt] = &All.OutputListOn;
      id[nt++] = INT;

      strcpy(tag[nt], "Omega0");
      addr[nt] = &All.Omega0;
      id[nt++] = REAL;

      strcpy(tag[nt], "OmegaBaryon");
      addr[nt] = &All.OmegaBaryon;
      id[nt++] = REAL;

      strcpy(tag[nt], "OmegaLambda");
      addr[nt] = &All.OmegaLambda;
      id[nt++] = REAL;

      strcpy(tag[nt], "HubbleParam");
      addr[nt] = &All.HubbleParam;
      id[nt++] = REAL;

      strcpy(tag[nt], "BoxSize");
      addr[nt] = &All.BoxSize;
      id[nt++] = REAL;

      strcpy(tag[nt], "PeriodicBoundariesOn");
      addr[nt] = &All.PeriodicBoundariesOn;
      id[nt++] = INT;

      strcpy(tag[nt], "TimeOfFirstSnapshot");
      addr[nt] = &All.TimeOfFirstSnapshot;
      id[nt++] = REAL;

      strcpy(tag[nt], "CpuTimeBetRestartFile");
      addr[nt] = &All.CpuTimeBetRestartFile;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeBetStatistics");
      addr[nt] = &All.TimeBetStatistics;
      id[nt++] = REAL;

/*    SHC 09/04/2009 the fraction of photon dump */
      strcpy(tag[nt], "FractionOfPhotons");
      addr[nt] = &All.FractionOfPhotons;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeBegin");
      addr[nt] = &All.TimeBegin;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeMax");
      addr[nt] = &All.TimeMax;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeBetSnapshot");
      addr[nt] = &All.TimeBetSnapshot;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
      addr[nt] = &All.UnitVelocity_in_cm_per_s;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitLength_in_cm");
      addr[nt] = &All.UnitLength_in_cm;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitMass_in_g");
      addr[nt] = &All.UnitMass_in_g;
      id[nt++] = REAL;

      strcpy(tag[nt], "TreeDomainUpdateFrequency");
      addr[nt] = &All.TreeDomainUpdateFrequency;
      id[nt++] = REAL;

      strcpy(tag[nt], "ErrTolIntAccuracy");
      addr[nt] = &All.ErrTolIntAccuracy;
      id[nt++] = REAL;

      strcpy(tag[nt], "ErrTolTheta");
      addr[nt] = &All.ErrTolTheta;
      id[nt++] = REAL;

#ifdef SUBFIND
      strcpy(tag[nt], "ErrTolThetaSubfind");
      addr[nt] = &All.ErrTolThetaSubfind;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "ErrTolForceAcc");
      addr[nt] = &All.ErrTolForceAcc;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinGasHsmlFractional");
      addr[nt] = &All.MinGasHsmlFractional;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxSizeTimestep");
      addr[nt] = &All.MaxSizeTimestep;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinSizeTimestep");
      addr[nt] = &All.MinSizeTimestep;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxRMSDisplacementFac");
      addr[nt] = &All.MaxRMSDisplacementFac;
      id[nt++] = REAL;

      strcpy(tag[nt], "ArtBulkViscConst");
      addr[nt] = &All.ArtBulkViscConst;
      id[nt++] = REAL;

      strcpy(tag[nt], "CourantFac");
      addr[nt] = &All.CourantFac;
      id[nt++] = REAL;

      strcpy(tag[nt], "DesNumNgb");
      addr[nt] = &All.DesNumNgb;
      id[nt++] = INT;

#ifdef SUBFIND
      strcpy(tag[nt], "DesLinkNgb");
      addr[nt] = &All.DesLinkNgb;
      id[nt++] = INT;
#endif

      strcpy(tag[nt], "MaxNumNgbDeviation");
      addr[nt] = &All.MaxNumNgbDeviation;
      id[nt++] = REAL;

      strcpy(tag[nt], "ComovingIntegrationOn");
      addr[nt] = &All.ComovingIntegrationOn;
      id[nt++] = INT;

      strcpy(tag[nt], "ICFormat");
      addr[nt] = &All.ICFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "SnapFormat");
      addr[nt] = &All.SnapFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesPerSnapshot");
      addr[nt] = &All.NumFilesPerSnapshot;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesWrittenInParallel");
      addr[nt] = &All.NumFilesWrittenInParallel;
      id[nt++] = INT;

      strcpy(tag[nt], "ResubmitOn");
      addr[nt] = &All.ResubmitOn;
      id[nt++] = INT;

      strcpy(tag[nt], "CoolingOn");
      addr[nt] = &All.CoolingOn;
      id[nt++] = INT;

      strcpy(tag[nt], "StarformationOn");
      addr[nt] = &All.StarformationOn;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfTimestepCriterion");
      addr[nt] = &All.TypeOfTimestepCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfOpeningCriterion");
      addr[nt] = &All.TypeOfOpeningCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "TimeLimitCPU");
      addr[nt] = &All.TimeLimitCPU;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningHalo");
      addr[nt] = &All.SofteningHalo;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningDisk");
      addr[nt] = &All.SofteningDisk;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningBulge");
      addr[nt] = &All.SofteningBulge;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningGas");
      addr[nt] = &All.SofteningGas;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningStars");
      addr[nt] = &All.SofteningStars;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningBndry");
      addr[nt] = &All.SofteningBndry;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningHaloMaxPhys");
      addr[nt] = &All.SofteningHaloMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningDiskMaxPhys");
      addr[nt] = &All.SofteningDiskMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningBulgeMaxPhys");
      addr[nt] = &All.SofteningBulgeMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningGasMaxPhys");
      addr[nt] = &All.SofteningGasMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningStarsMaxPhys");
      addr[nt] = &All.SofteningStarsMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningBndryMaxPhys");
      addr[nt] = &All.SofteningBndryMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "BufferSize");
      addr[nt] = &All.BufferSize;
      id[nt++] = INT;

      strcpy(tag[nt], "PartAllocFactor");
      addr[nt] = &All.PartAllocFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "GravityConstantInternal");
      addr[nt] = &All.GravityConstantInternal;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitGasTemp");
      addr[nt] = &All.InitGasTemp;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinGasTemp");
      addr[nt] = &All.MinGasTemp;
      id[nt++] = REAL;

      strcpy(tag[nt], "MeanWeight");
      addr[nt] = &All.MeanWeight;
      id[nt++] = REAL;

#ifdef SGRA_POTENTIAL
      strcpy(tag[nt], "CuspMass");
      addr[nt] = &All.CuspMass;
      id[nt++] = REAL;

      strcpy(tag[nt], "CuspRadius");
      addr[nt] = &All.CuspRadius;
      id[nt++] = REAL;

      strcpy(tag[nt], "CuspAlpha1");
      addr[nt] = &All.CuspAlpha1;
      id[nt++] = REAL;

      strcpy(tag[nt], "CuspAlpha2");
      addr[nt] = &All.CuspAlpha2;
      id[nt++] = REAL;

#endif 

#if defined(CUSP_POTENTIAL) || defined(SIS_POTENTIAL) 
      strcpy(tag[nt], "CuspMass");
      addr[nt] = &All.CuspMass;
      id[nt++] = REAL;

      strcpy(tag[nt], "CuspRadius");
      addr[nt] = &All.CuspRadius;
      id[nt++] = REAL;

      strcpy(tag[nt], "FlatRadius");
      addr[nt] = &All.FlatRadius;
      id[nt++] = REAL;

      strcpy(tag[nt], "CuspGamma");
      addr[nt] = &All.CuspGamma;
      id[nt++] = REAL;
#endif

/* Addition by CBP. 27/07/2009. 
/ Input parameters needed for a Navarro, Frenk & White halo.
/ Note that the virial and scale radii are calculated in begrun.c -- as soon
/ as input parameters are read in.
*/

#ifdef NFW_POTENTIAL
      strcpy(tag[nt], "VirialMass");  // Virial mass of halo.
      addr[nt] = &All.VirialMass;       // mvir = (4PI/3)*DCrit*rhocrit*rvir**3
      id[nt++] = REAL;

      strcpy(tag[nt], "NFWConcentration"); // NFW concentration
      addr[nt] = &All.CNFW;
      id[nt++] = REAL;

      strcpy(tag[nt], "DeltaVir"); // Virial Overdensity Criterion
      addr[nt] = &All.DeltaVir;
      id[nt++] = REAL;

      strcpy(tag[nt], "FlatRadius");
      addr[nt] = &All.FlatRadius;
      id[nt++] = REAL;
#endif

/* End of Addition by CBP. 27/07/2009. */ 

#ifdef SMBH_IRRADIATION
      strcpy(tag[nt], "BH_LuminosityFactor");
      addr[nt] = &All.BH_LuminosityFactor;
      id[nt++] = REAL;
#endif



#ifdef SCALARFIELD
      strcpy(tag[nt], "ScalarBeta");
      addr[nt] = &All.ScalarBeta;
      id[nt++] = REAL;

      strcpy(tag[nt], "ScalarScreeningLength");
      addr[nt] = &All.ScalarScreeningLength;
      id[nt++] = REAL;
#endif

#ifdef OUTPUTLINEOFSIGHT
      strcpy(tag[nt], "TimeFirstLineOfSight");
      addr[nt] = &All.TimeFirstLineOfSight;
      id[nt++] = REAL;
#endif


#if defined(ADAPTIVE_GRAVSOFT_FORGAS) && !defined(ADAPTIVE_GRAVSOFT_FORGAS_HSML)
      strcpy(tag[nt], "ReferenceGasMass");
      addr[nt] = &All.ReferenceGasMass;
      id[nt++] = REAL;
#endif

#ifdef NAVIERSTOKES
      strcpy(tag[nt], "FractionSpitzerViscosity");
      addr[nt] = &All.FractionSpitzerViscosity;
      id[nt++] = REAL;
#endif

#ifdef NAVIERSTOKES_CONSTANT
      strcpy(tag[nt], "ShearViscosityTemperature");
      addr[nt] = &All.ShearViscosityTemperature;
      id[nt++] = REAL;
#endif

#ifdef NAVIERSTOKES_BULK
      strcpy(tag[nt], "NavierStokes_BulkViscosity");
      addr[nt] = &All.NavierStokes_BulkViscosity;
      id[nt++] = REAL;
#endif

#ifdef CHEMISTRY
      strcpy(tag[nt], "Epsilon");
      addr[nt] = &All.Epsilon;
      id[nt++] = REAL;
#endif


#ifdef CONDUCTION
      strcpy(tag[nt], "ConductionEfficiency");
      addr[nt] = &All.ConductionCoeff;
      id[nt++] = REAL;
#endif

#if defined(BUBBLES) || defined(MULTI_BUBBLES)
      strcpy(tag[nt], "BubbleDistance");
      addr[nt] = &All.BubbleDistance;
      id[nt++] = REAL;

      strcpy(tag[nt], "BubbleRadius");
      addr[nt] = &All.BubbleRadius;
      id[nt++] = REAL;

      strcpy(tag[nt], "BubbleTimeInterval");
      addr[nt] = &All.BubbleTimeInterval;
      id[nt++] = REAL;

      strcpy(tag[nt], "BubbleEnergy");
      addr[nt] = &All.BubbleEnergy;
      id[nt++] = REAL;

      strcpy(tag[nt], "FirstBubbleRedshift");
      addr[nt] = &All.FirstBubbleRedshift;
      id[nt++] = REAL;
#endif

#ifdef MULTI_BUBBLES
      strcpy(tag[nt], "MinFoFMassForNewSeed");
      addr[nt] = &All.MinFoFMassForNewSeed;
      id[nt++] = REAL;

      strcpy(tag[nt], "ClusterMass200");
      addr[nt] = &All.ClusterMass200;
      id[nt++] = REAL;

      strcpy(tag[nt], "massDMpart");
      addr[nt] = &All.massDMpart;
      id[nt++] = REAL;

#endif


#ifdef COSMIC_RAYS
      strcpy(tag[nt], "CR_SpectralIndex");
      addr[nt] = &All.CR_Alpha;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_SupernovaEfficiency");
      addr[nt] = &All.CR_SNEff;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_SupernovaSpectralIndex");
      addr[nt] = &All.CR_SNAlpha;
      id[nt++] = REAL;

#if defined(CR_DIFFUSION) || defined(CR_DIFFUSION_GREEN)
      strcpy(tag[nt], "CR_DiffusionCoeff");
      addr[nt] = &All.CR_DiffusionCoeff;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_DiffusionDensityScaling");
      addr[nt] = &All.CR_DiffusionDensScaling;
      id[nt++] = REAL;

      /* CR Diffusion scaling: reference density rho_0.
         If value is 0, then rho_crit @ z=0 is used. */
      strcpy(tag[nt], "CR_DiffusionReferenceDensity");
      addr[nt] = &All.CR_DiffusionDensZero;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_DiffusionTemperatureScaling");
      addr[nt] = &All.CR_DiffusionEntropyScaling;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_DiffusionReferenceTemperature");
      addr[nt] = &All.CR_DiffusionEntropyZero;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_DiffusionTimeScale");
      addr[nt] = &All.CR_DiffusionTimeScale;
      id[nt++] = REAL;
#endif /* CR_DIFFUSION */

#ifdef CR_SHOCK
      strcpy(tag[nt], "CR_ShockEfficiency");
      addr[nt] = &All.CR_ShockEfficiency;
      id[nt++] = REAL;
#if ( CR_SHOCK == 1 )		/* Constant Spectral Index method */
      strcpy(tag[nt], "CR_ShockSpectralIndex");
      addr[nt] = &All.CR_ShockAlpha;
      id[nt++] = REAL;
#else /* Mach-Number - Dependent method */
      strcpy(tag[nt], "CR_ShockCutoffFac");
      addr[nt] = &All.CR_ShockCutoff;
      id[nt++] = REAL;
#endif
#endif /* CR_SHOCK */

#ifdef FIX_QINJ
      strcpy(tag[nt], "Shock_Fix_Qinj");
      addr[nt] = &All.Shock_Fix_Qinj;
      id[nt++] = REAL;
#endif

#endif /* COSMIC_RAYS */


#ifdef MACHNUM
      strcpy(tag[nt], "Shock_LengthScale");
      addr[nt] = &All.Shock_Length;
      id[nt++] = REAL;

      strcpy(tag[nt], "Shock_DeltaDecayTimeMax");
      addr[nt] = &All.Shock_DeltaDecayTimeMax;
      id[nt++] = REAL;
#endif


#ifdef BLACK_HOLES
      strcpy(tag[nt], "TimeBetBlackHoleSearch");
      addr[nt] = &All.TimeBetBlackHoleSearch;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleAccretionFactor");
      addr[nt] = &All.BlackHoleAccretionFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleFeedbackFactor");
      addr[nt] = &All.BlackHoleFeedbackFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleEddingtonFactor");
      addr[nt] = &All.BlackHoleEddingtonFactor;
      id[nt++] = REAL;


      strcpy(tag[nt], "SeedBlackHoleMass");
      addr[nt] = &All.SeedBlackHoleMass;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinFoFMassForNewSeed");
      addr[nt] = &All.MinFoFMassForNewSeed;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleNgbFactor");
      addr[nt] = &All.BlackHoleNgbFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleActiveTime");
      addr[nt] = &All.BlackHoleActiveTime;
      id[nt++] = REAL;

      strcpy(tag[nt], "SMBHmass");
      addr[nt] = &All.SMBHmass;
      id[nt++] = REAL;
		

#if defined(PLANET_INITIAL_MASS) 
      strcpy(tag[nt], "BlackHoleMass");
      addr[nt] = &All.BlackHoleMass;
      id[nt++] = REAL;
#endif

#if defined(BH_KICK) 
      strcpy(tag[nt], "BlackHoleMass");
      addr[nt] = &All.BlackHoleMass;
      id[nt++] = REAL;
      strcpy(tag[nt], "BlackHoleKickVel");
      addr[nt] = &All.BlackHoleKickVel;
      id[nt++] = REAL;
      strcpy(tag[nt], "BlackHoleKickAngle");
      addr[nt] = &All.BlackHoleKickAngle;
      id[nt++] = REAL;
#endif

#ifdef ACCRETION_RADIUS
      strcpy(tag[nt], "InnerBoundary");
      addr[nt] = &All.InnerBoundary;
      id[nt++] = REAL;
#ifdef BH_FORM
      strcpy(tag[nt], "SinkBoundary");
      addr[nt] = &All.SinkBoundary;
      id[nt++] = REAL;
#endif
// Addition by CBP (30/10/2008)
      strcpy(tag[nt], "AccDtBlackHole");
      addr[nt] = &All.AccDtBlackHole;
      id[nt++] = REAL;
      strcpy(tag[nt], "ViscousTimescale");
      addr[nt] = &All.ViscousTimescale;
      id[nt++] = REAL;
      strcpy(tag[nt], "AccDiscSFR");
      addr[nt] = &All.AccDiscSFR;
      id[nt++] = REAL;
#endif

/* #ifdef VARIABLE_PHOTON_SEARCH_RADIUS */
/*       strcpy(tag[nt], "PhotonSearchRadius"); */
/*       addr[nt] = &All.PhotonSearchRadius; */
/*       id[nt++] = REAL; */
/* #endif      */ 

// Addition by CBP (14/11/2008)
// Parameters for beta-cooling model -- cooling occurs on a dynamical timescale.
#if defined(BETA_COOLING) || defined(EVAPORATION) || defined(EVAPORATION_RADIAL)
      strcpy(tag[nt], "BetaCool");   // Beta parameter
      addr[nt] = &All.BetaCool;
      id[nt++] = REAL;
#endif
// End of CBP      

#if defined(EVAPORATION) || defined(EVAPORATION_RADIAL)
      strcpy(tag[nt], "Evap_dens");   // Beta parameter
      addr[nt] = &All.Evap_dens;
      id[nt++] = REAL;
#endif

#if defined(EVAPORATION_RADIAL)
      strcpy(tag[nt], "Cool_ind");  
      addr[nt] = &All.Cool_ind;
      id[nt++] = REAL;
#endif 

#if defined(EVAPORATION_RADIAL)
      strcpy(tag[nt], "rho_cool_ind");  
      addr[nt] = &All.rho_cool_ind;
      id[nt++] = REAL;
#endif 

#if defined(BACKGROUND_ILLUMINATION_TEMP) || defined(BETA_COOLING) || defined(ISOTHERM) || defined(EVAPORATION) || defined(EVAPORATION_RADIAL)
      strcpy(tag[nt], "EquilibriumTemp");  // Equilibrium temperature 
      addr[nt] = &All.EqTemp;
      id[nt++] = REAL;
#endif

#ifdef POLYTROP //RJH 20/03/2018
      strcpy(tag[nt], "Poly_K");   // Polytropic cooling K const. in cgs
      addr[nt] = &All.Poly_K;
      id[nt++] = REAL;
#endif



#ifdef VIRTUAL
      strcpy(tag[nt], "VirtualStart");
      addr[nt] = &All.VirtualStart;
      id[nt++] = REAL;
      strcpy(tag[nt], "FeedBackVelocity");
      addr[nt] = &All.FeedBackVelocity;
      id[nt++] = REAL;
      strcpy(tag[nt], "VirtualMomentum");
      addr[nt] = &All.VirtualMomentum;
      id[nt++] = REAL;
      strcpy(tag[nt], "VirtualMass");
      addr[nt] = &All.VirtualMass;
      id[nt++] = REAL;
      strcpy(tag[nt], "VirtualTime");
      addr[nt] = &All.VirtualTime;
      id[nt++] = REAL;
      strcpy(tag[nt], "VirtualCrosSection");
      addr[nt] = &All.VirtualCrosSection;
      id[nt++] = REAL;
      strcpy(tag[nt], "VirtualFeedBack");
      addr[nt] = &All.VirtualFeedBack;
      id[nt++] = REAL;
      strcpy(tag[nt], "OuterBoundary");
      addr[nt] = &All.OuterBoundary;
      id[nt++] = REAL;
      strcpy(tag[nt], "FreeTravelLength");
      addr[nt] = &All.FreeTravelLength;
      id[nt++] = REAL;
#ifdef DUST_GROWTH
      strcpy(tag[nt], "FragmentationVelocity");
      addr[nt] = &All.FragmentationVelocity;
      id[nt++] = REAL;
#endif
#ifdef DUST_SINK_ON_FLY
      strcpy(tag[nt], "CritDustDensity");
      addr[nt] = &All.CritDustDensity;
      id[nt++] = REAL;
#endif
#endif


#ifdef FOF
      strcpy(tag[nt], "massDMpart");
      addr[nt] = &All.massDMpart;
      id[nt++] = REAL;
#endif
#endif


#ifdef SFR
      strcpy(tag[nt], "CritOverDensity");
      addr[nt] = &All.CritOverDensity;
      id[nt++] = REAL;

      strcpy(tag[nt], "CritPhysDensity");
      addr[nt] = &All.CritPhysDensity;
      id[nt++] = REAL;

      strcpy(tag[nt], "FactorSN");
      addr[nt] = &All.FactorSN;
      id[nt++] = REAL;

      strcpy(tag[nt], "FactorEVP");
      addr[nt] = &All.FactorEVP;
      id[nt++] = REAL;

      strcpy(tag[nt], "TempSupernova");
      addr[nt] = &All.TempSupernova;
      id[nt++] = REAL;

      strcpy(tag[nt], "TempClouds");
      addr[nt] = &All.TempClouds;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxSfrTimescale");
      addr[nt] = &All.MaxSfrTimescale;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindEfficiency");
      addr[nt] = &All.WindEfficiency;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindEnergyFraction");
      addr[nt] = &All.WindEnergyFraction;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindFreeTravelLength");
      addr[nt] = &All.WindFreeTravelLength;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindFreeTravelDensFac");
      addr[nt] = &All.WindFreeTravelDensFac;
      id[nt++] = REAL;
#endif

#ifdef SOFTEREQS
      strcpy(tag[nt], "FactorForSofterEQS");
      addr[nt] = &All.FactorForSofterEQS;
      id[nt++] = REAL;
#endif
#ifdef DARKENERGY
#ifndef TIMEDEPDE
      strcpy(tag[nt], "DarkEnergyParam");
      addr[nt] = &All.DarkEnergyParam;
      id[nt++] = REAL;
#endif
#endif

#ifdef RESCALEVINI
      strcpy(tag[nt], "VelIniScale");
      addr[nt] = &All.VelIniScale;
      id[nt++] = REAL;
#endif

#ifdef DARKENERGY
#ifdef TIMEDEPDE
      strcpy(tag[nt], "DarkEnergyFile");
      addr[nt] = All.DarkEnergyFile;
      id[nt++] = STRING;
#endif
#endif

#ifdef TIME_DEP_ART_VISC
      strcpy(tag[nt], "ViscositySourceScaling");
      addr[nt] = &All.ViscSource0;
      id[nt++] = REAL;

      strcpy(tag[nt], "ViscosityDecayLength");
      addr[nt] = &All.DecayLength;
      id[nt++] = REAL;

      strcpy(tag[nt], "ViscosityAlphaMin");
      addr[nt] = &All.AlphaMin;
      id[nt++] = REAL;
#endif

#ifdef MAGNETIC_DISSIPATION
      strcpy(tag[nt], "ArtificialMagneticDissipationConstant");
      addr[nt] = &All.ArtMagDispConst;
      id[nt++] = REAL;

#ifdef TIME_DEP_MAGN_DISP
      strcpy(tag[nt], "ArtificialMagneticDissipationMin");
      addr[nt] = &All.ArtMagDispMin;
      id[nt++] = REAL;

      strcpy(tag[nt], "ArtificialMagneticDissipationSource");
      addr[nt] = &All.ArtMagDispSource;
      id[nt++] = REAL;

      strcpy(tag[nt], "ArtificialMagneticDissipationDecaytime");
      addr[nt] = &All.ArtMagDispTime;
      id[nt++] = REAL;
#endif
#endif

#ifdef DIVBCLEANING_DEDNER
      strcpy(tag[nt], "DivBcleaningParabolicSigma");
      addr[nt] = &All.DivBcleanParabolicSigma;
      id[nt++] = REAL;

      strcpy(tag[nt], "DivBcleaningHyperbolicSigma");
      addr[nt] = &All.DivBcleanHyperbolicSigma;
      id[nt++] = REAL;
#endif

#ifdef MAGNETIC
#ifdef BINISET
      strcpy(tag[nt], "BiniX");
      addr[nt] = &All.BiniX;
      id[nt++] = REAL;

      strcpy(tag[nt], "BiniY");
      addr[nt] = &All.BiniY;
      id[nt++] = REAL;

      strcpy(tag[nt], "BiniZ");
      addr[nt] = &All.BiniZ;
      id[nt++] = REAL;
#endif

#ifdef BSMOOTH
      strcpy(tag[nt], "BSmoothInt");
      addr[nt] = &All.BSmoothInt;
      id[nt++] = INT;

      strcpy(tag[nt], "BSmoothFrac");
      addr[nt] = &All.BSmoothFrac;
      id[nt++] = REAL;
#endif
#endif

#ifdef DUST
      strcpy(tag[nt], "InitialDustRadius");
      addr[nt] = &All.InitialDustRadius;
      id[nt++] = REAL;
      strcpy(tag[nt], "InitialDustMass");
      addr[nt] = &All.InitialDustMass;
      id[nt++] = REAL;
#ifdef DUST_POWERLAW
      strcpy(tag[nt], "aDustMin0");
      addr[nt] = &All.aDustMin0;
      id[nt++] = REAL;
      strcpy(tag[nt], "aDustMax0");
      addr[nt] = &All.aDustMax0;
      id[nt++] = REAL;
#endif
#ifdef DUST_TWO_POPULATIONS
      strcpy(tag[nt], "InitialMicroDustZ");
      addr[nt] = &All.InitialMicroDustZ;
      id[nt++] = REAL;
#endif
#ifdef DUST_PEBBLES_BORN
      strcpy(tag[nt], "DustInjectionTime");
      addr[nt] = &All.DustInjectionTime;
      id[nt++] = REAL;
      strcpy(tag[nt], "PebbleInjectionMass");
      addr[nt] = &All.PebbleInjectionMass;
      id[nt++] = REAL;
      strcpy(tag[nt], "PebbleInjectionDens");
      addr[nt] = &All.PebbleInjectionDens;
      id[nt++] = REAL;
      strcpy(tag[nt], "PebbleInjectionSlab");
      addr[nt] = &All.PebbleInjectionSlab;
      id[nt++] = REAL;
#endif
#endif

#ifdef INCREASE_SPH_MASS
      /* time scale for SPH particle mass growth; SN */
      strcpy(tag[nt], "DeltaTimeGrow");
      addr[nt] = &All.DeltaTimeGrow;
      id[nt++] = REAL;
#endif

      if((fd = fopen(fname, "r")))
	{
	  sprintf(buf, "%s%s", fname, "-usedvalues");
	  if(!(fdout = fopen(buf, "w")))
	    {
	      printf("error opening file '%s' \n", buf);
	      errorFlag = 1;
	    }
	  else
	    {
	      while(!feof(fd))
		{
		  char *ret;

		  *buf = 0;
		  ret = fgets(buf, 200, fd);
		  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
		    continue;

		  if(buf1[0] == '%')
		    continue;

		  for(i = 0, j = -1; i < nt; i++)
		    if(strcmp(buf1, tag[i]) == 0)
		      {
			j = i;
			tag[i][0] = 0;
			break;
		      }

		  if(j >= 0)
		    {
		      switch (id[j])
			{
			case REAL:
			  *((double *) addr[j]) = atof(buf2);
			  fprintf(fdout, "%-35s%g\n", buf1, *((double *) addr[j]));
			  break;
			case STRING:
			  strcpy((char *) addr[j], buf2);
			  fprintf(fdout, "%-35s%s\n", buf1, buf2);
			  break;
			case INT:
			  *((int *) addr[j]) = atoi(buf2);
			  fprintf(fdout, "%-35s%d\n", buf1, *((int *) addr[j]));
			  break;
			}
		    }
		  else
		    {
#ifdef ALLOWEXTRAPARAMS
		      fprintf(stdout, "WARNING from file %s:   Tag '%s' ignored !\n", fname, buf1);
#else
		      fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n",
			      fname, buf1);
		      errorFlag = 1;
#endif
		    }
		}
	      fclose(fd);
	      fclose(fdout);

	      i = strlen(All.OutputDir);
	      if(i > 0)
		if(All.OutputDir[i - 1] != '/')
		  strcat(All.OutputDir, "/");

	      sprintf(buf1, "%s%s", fname, "-usedvalues");
	      sprintf(buf2, "%s%s", All.OutputDir, "parameters-usedvalues");
	      sprintf(buf3, "cp %s %s", buf1, buf2);
#ifndef NOCALLSOFSYSTEM
	      int ret;

	      ret = system(buf3);
#endif
	    }
	}
      else
	{
	  printf("Parameter file %s not found.\n", fname);
	  errorFlag = 1;
	}


      for(i = 0; i < nt; i++)
	{
	  if(*tag[i])
	    {
	      printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
	      errorFlag = 1;
	    }
	}

      if(All.OutputListOn && errorFlag == 0)
	errorFlag += read_outputlist(All.OutputListFilename);
      else
	All.OutputListLength = 0;
    }

  MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(errorFlag)
    {
      MPI_Finalize();
      exit(0);
    }

  /* now communicate the relevant parameters to the other processes */
  MPI_Bcast(&All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);



  for(pnum = 0; All.NumFilesWrittenInParallel > (1 << pnum); pnum++);

  if(All.NumFilesWrittenInParallel != (1 << pnum))
    {
      if(ThisTask == 0)
	printf("NumFilesWrittenInParallel MUST be a power of 2\n");
      endrun(0);
    }

  if(All.NumFilesWrittenInParallel > NTask)
    {
      if(ThisTask == 0)
	printf("NumFilesWrittenInParallel MUST be smaller than number of processors\n");
      endrun(0);
    }

#ifdef PERIODIC
  if(All.PeriodicBoundariesOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with periodic boundary conditions switched on.\n");
	  printf("You must set `PeriodicBoundariesOn=1', or recompile the code.\n");
	}
      endrun(0);
    }
#else
  if(All.PeriodicBoundariesOn == 1)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with periodic boundary conditions switched off.\n");
	  printf("You must set `PeriodicBoundariesOn=0', or recompile the code.\n");
	}
      endrun(0);
    }
#endif

#ifdef COOLING
  if(All.CoolingOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with cooling switched on.\n");
	  printf("You must set `CoolingOn=1', or recompile the code.\n");
	}
      endrun(0);
    }
#else
  if(All.CoolingOn == 1)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with cooling switched off.\n");
	  printf("You must set `CoolingOn=0', or recompile the code.\n");
	}
      endrun(0);
    }
#endif

  if(All.TypeOfTimestepCriterion >= 3)
    {
      if(ThisTask == 0)
	{
	  printf("The specified timestep criterion\n");
	  printf("is not valid\n");
	}
      endrun(0);
    }

#if defined(LONG_X) ||  defined(LONG_Y) || defined(LONG_Z)
#ifndef NOGRAVITY
  if(ThisTask == 0)
    {
      printf("Code was compiled with LONG_X/Y/Z, but not with NOGRAVITY.\n");
      printf("Stretched periodic boxes are not implemented for gravity yet.\n");
    }
  endrun(0);
#endif
#endif

#ifdef SFR

#ifndef MOREPARAMS
  if(ThisTask == 0)
    {
      printf("Code was compiled with SFR, but not with MOREPARAMS.\n");
      printf("This is not allowed.\n");
    }
  endrun(0);
#endif

  if(All.StarformationOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with star formation switched on.\n");
	  printf("You must set `StarformationOn=1', or recompile the code.\n");
	}
      endrun(0);
    }
  if(All.CoolingOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("You try to use the code with star formation enabled,\n");
	  printf("but you did not switch on cooling.\nThis mode is not supported.\n");
	}
      endrun(0);
    }
#else
  if(All.StarformationOn == 1)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with star formation switched off.\n");
	  printf("You must set `StarformationOn=0', or recompile the code.\n");
	}
      endrun(0);
    }
#endif


#ifdef METALS
#ifndef SFR
  if(ThisTask == 0)
    {
      printf("Code was compiled with METALS, but not with SFR.\n");
      printf("This is not allowed.\n");
    }
  endrun(0);
#endif
#endif

#ifndef MOREPARAMS
#ifdef TIME_DEP_ART_VISC
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with TIME_DEP_ART_VISC, but not with MOREPARAMS.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef DARKENERGY
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with DARKENERGY, but not with MOREPARAMS.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef TIMEDEPDE
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with TIMEDEPDE, but not with MOREPARAMS.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif
#endif

#ifdef TIMEDEPDE
#ifndef DARKENERGY
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with TIMEDEPDE, but not with DARKENERGY.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif
#endif

#ifndef MAGNETIC
#ifdef TRACEDIVB
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with TRACEDIVB, but not with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef DBOUTPUT
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with DBOUTPUT, but not with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef MAGFORCE
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with MAGFORCE, but not with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef BSMOOTH
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with BSMOOTH, but not with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef BFROMROTA
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with BFROMROTA, but not with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef MU0_UNITY
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with MU0_UNITY, but not with MAGNETIC.\n");
      fprintf(stdout, "This makes no sense.\n");
    }
  endrun(0);
#endif

#ifdef MAGNETIC_DISSIPATION
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with MAGNETIC_DISSIPATION, but not with MAGNETIC.\n");
      fprintf(stdout, "This makes no sense.\n");
    }
  endrun(0);
#endif

#ifdef TIME_DEP_MAGN_DISP
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with TIME_DEP_MAGN_DISP, but not with MAGNETIC.\n");
      fprintf(stdout, "This makes no sense.\n");
    }
  endrun(0);
#endif

#ifdef DIVBCLEANING_DEDNER
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with DIVBCLEANING_DEDNER, but not with MAGNETIC.\n");
      fprintf(stdout, "This makes no sense.\n");
    }
  endrun(0);
#endif

#endif

#if defined(NOWINDTIMESTEPPING) && defined(MAGNETIC)
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with NOWINDTIMESTEPPING and with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed, as it leads to inconsitent MHD for wind particles.\n");
    }
  endrun(0);
#endif

#ifndef MAGFORCE
#ifdef ARTBPRES
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with ARTBPRES, but not with MAGFORCE.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef DIVBFORCE
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with DIVBFORCE, but not with MAGFORCE.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif
#endif


#undef REAL
#undef STRING
#undef INT
#undef MAXTAGS


#ifdef COSMIC_RAYS
  if(ThisTask == 0)
    {
      printf("CR SN Efficiency: %g\n", All.CR_SNEff);
    }
#endif
}


/*! this function reads a table with a list of desired output times. The table
 *  does not have to be ordered in any way, but may not contain more than
 *  MAXLEN_OUTPUTLIST entries.
 */
int read_outputlist(char *fname)
{
  FILE *fd;

  if(!(fd = fopen(fname, "r")))
    {
      printf("can't read output list in file '%s'\n", fname);
      return 1;
    }

  All.OutputListLength = 0;
  do
    {
      if(fscanf(fd, " %lg ", &All.OutputListTimes[All.OutputListLength]) == 1)
	All.OutputListLength++;
      else
	break;
    }
  while(All.OutputListLength < MAXLEN_OUTPUTLIST);

  fclose(fd);

  printf("\nfound %d times in output-list.\n", All.OutputListLength);

  return 0;
}


/*! If a restart from restart-files is carried out where the TimeMax variable
 * is increased, then the integer timeline needs to be adjusted. The approach
 * taken here is to reduce the resolution of the integer timeline by factors
 * of 2 until the new final time can be reached within TIMEBASE.
 */
void readjust_timebase(double TimeMax_old, double TimeMax_new)
{
  int i;
  long long ti_end;

  if(sizeof(long long) != 8)
    {
      if(ThisTask == 0)
	printf("\nType 'long long' is not 64 bit on this platform\n\n");
      endrun(555);
    }

  if(ThisTask == 0)
    {
      printf("\nAll.TimeMax has been changed in the parameterfile\n");
      printf("Need to adjust integer timeline\n\n\n");
    }

  if(TimeMax_new < TimeMax_old)
    {
      if(ThisTask == 0)
	printf("\nIt is not allowed to reduce All.TimeMax\n\n");
      endrun(556);
    }

  if(All.ComovingIntegrationOn)
    ti_end = (long long) (log(TimeMax_new / All.TimeBegin) / All.Timebase_interval);
  else
    ti_end = (long long) ((TimeMax_new - All.TimeBegin) / All.Timebase_interval);

  while(ti_end > TIMEBASE)
    {
      All.Timebase_interval *= 2.0;

      ti_end /= 2;
      All.Ti_Current /= 2;

#ifdef PMGRID
      All.PM_Ti_begstep /= 2;
      All.PM_Ti_endstep /= 2;
#endif

      for(i = 0; i < NumPart; i++)
	{
	  P[i].Ti_begstep /= 2;
	  P[i].TimeBin--;
	  if(P[i].TimeBin == 0)
	    {
	      printf("Error in readjust_timebase(). Minimum Timebin for particle %d reached.\n", i);
	      endrun(8765);
	    }
	}

      All.Ti_nextlineofsight /= 2;
    }

  All.TimeMax = TimeMax_new;
}
