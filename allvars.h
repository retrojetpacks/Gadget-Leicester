/*! \file allvars.h
 *  \brief declares global variables.
 *
 *  This file declares all global variables. Further variables should be added here, and declared as
 *  'extern'. The actual existence of these variables is provided by the file 'allvars.c'. To produce
 *  'allvars.c' from 'allvars.h', do the following:
 *
 *     - Erase all #define statements
 *     - add #include "allvars.h" 
 *     - delete all keywords 'extern'
 *     - delete all struct definitions enclosed in {...}, e.g.
 *        "extern struct global_data_all_processes {....} All;"
 *        becomes "struct global_data_all_processes All;"
 */

#ifndef ALLVARS_H
#define ALLVARS_H

#include <mpi.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

#include "tags.h"
#ifdef CHEMISTRY
#include "chemistry.h"
#endif
#include "assert.h"



#define  GADGETVERSION   "3.0"	/*!< code version string */

#define  GENERATIONS     2	/*!< Number of star particles that may be created per gas particle */

#define  TIMEBINS        29

#define  TIMEBASE        (1<<TIMEBINS)	/*!< The simulated timespan is mapped onto the integer interval [0,TIMESPAN],
					 *   where TIMESPAN needs to be a power of 2. Note that (1<<28) corresponds
					 *   to 2^29
					 */
#ifndef  MULTIPLEDOMAINS
#define  MULTIPLEDOMAINS     1
#endif

#ifndef  TOPNODEFACTOR
#define  TOPNODEFACTOR       2.5
#endif

#define  NODELISTLENGTH      8

typedef unsigned long long peanokey;


#define  BITS_PER_DIMENSION 21	/* for Peano-Hilbert order. Note: Maximum is 10 to fit in 32-bit integer ! */
#define  PEANOCELLS (((peanokey)1)<<(3*BITS_PER_DIMENSION))




#define  GAMMA          7./5.   //5.0/3	/*!< adiabatic index of simulated gas */
#define  GAMMA_MINUS1  (GAMMA-1)

#define  HYDROGEN_MASSFRAC 0.76	/*!< mass fraction of hydrogen, relevant only for radiative cooling */

#define  METAL_YIELD       0.02	/*!< effective metal yield for star formation */

#define  MAX_REAL_NUMBER  1e37
#define  MIN_REAL_NUMBER  1e-37

//#define  RNDTABLE 8192
#define  RNDTABLE 262144

/* ... often used physical constants (cgs units) */

#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  C           2.9979e10
#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24
#define  PROTONMASS  1.6726e-24
#define  ELECTRONMASS 9.10953e-28
#define  THOMPSON     6.65245e-25
#define  ELECTRONCHARGE  4.8032e-10
#define  HUBBLE          3.2407789e-18	/* in h/sec */
#define  LYMAN_ALPHA      1215.6e-8	/* 1215.6 Angstroem */
#define  LYMAN_ALPHA_HeII  303.8e-8	/* 303.8 Angstroem */
#define  OSCILLATOR_STRENGTH       0.41615
#define  OSCILLATOR_STRENGTH_HeII  0.41615

#ifdef NAVIERSTOKES
#define  LOG_LAMBDA      37.8	/* logarithmic Coulomb factor */
#endif

#ifdef CHEMISTRY
#define  T_CMB0      2.728	/* present-day CMB temperature */
#endif

#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7


#ifndef FOF_PRIMARY_LINK_TYPES
#define FOF_PRIMARY_LINK_TYPES 2
#endif

#ifndef FOF_SECONDARY_LINK_TYPES
#define FOF_SECONDARY_LINK_TYPES 0
#endif


#ifndef ASMTH
/*! ASMTH gives the scale of the short-range/long-range force split in units of FFT-mesh cells */
#define ASMTH 1.25
#endif
#ifndef RCUT
/*! RCUT gives the maximum distance (in units of the scale used for the force split) out to which short-range
 * forces are evaluated in the short-range tree walk.
 */
#define RCUT  4.5
#endif

#define COND_TIMESTEP_PARAMETER 0.25
#define VISC_TIMESTEP_PARAMETER 0.25

#define MAXLEN_OUTPUTLIST 350	/*!< maxmimum number of entries in output list */

#define DRIFT_TABLE_LENGTH  1000	/*!< length of the lookup table used to hold the drift and kick factors */


#define MAXITER 150

#define LINKLENGTH 0.2
#define GROUP_MIN_LEN 32

#define MINRESTFAC 0.05


#ifndef LONGIDS
typedef unsigned int MyIDType;
#else
typedef unsigned long long MyIDType;
#endif


#ifndef DOUBLEPRECISION     /* default is single-precision */
typedef float  MyFloat;
typedef float  MyDouble;
#else
#if (DOUBLEPRECISION == 2)   /* mixed precision */
typedef float   MyFloat;
typedef double  MyDouble;
#else                        /* everything double-precision */
typedef double  MyFloat;
typedef double  MyDouble;
#endif
#endif

#ifdef OUTPUT_IN_DOUBLEPRECISION
typedef double MyOutputFloat;
#else
typedef float MyOutputFloat;
#endif

#ifdef INPUT_IN_DOUBLEPRECISION
typedef double MyInputFloat;
#else
typedef float MyInputFloat;
#endif

struct unbind_data
{
  int index;
};


#ifdef FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
extern MPI_Status mpistat;
#undef MPI_STATUS_IGNORE
#define MPI_STATUS_IGNORE &mpistat
#endif

#ifdef FLTROUNDOFFREDUCTION  
#define FLT(x) ((MyFloat)(x))
#ifdef SOFTDOUBLEDOUBLE      /* this requires a C++ compilation */
#include "dd.h"
typedef dd MyLongDouble;
#else
typedef long double MyLongDouble;
#endif
#else  /* not enabled */
#define FLT(x) (x)
typedef MyFloat MyLongDouble;
#endif  /* end FLTROUNDOFFREDUCTION */


#define CPU_ALL            0
#define CPU_TREEWALK1      1
#define CPU_TREEWALK2      2
#define CPU_TREEWAIT1      3
#define CPU_TREEWAIT2      4
#define CPU_TREESEND       5
#define CPU_TREERECV       6
#define CPU_TREEMISC       7
#define CPU_TREEBUILD      8
#define CPU_TREEUPDATE     9
#define CPU_TREEHMAXUPDATE 10
#define CPU_DOMAIN         11
#define CPU_DENSCOMPUTE    12
#define CPU_DENSWAIT       13
#define CPU_DENSCOMM       14
#define CPU_DENSMISC       15
#define CPU_HYDCOMPUTE     16
#define CPU_HYDWAIT        17
#define CPU_HYDCOMM        18
#define CPU_HYDMISC        19
#define CPU_DRIFT          20
#define CPU_BLACKHOLE      21
#define CPU_TIMELINE       22
#define CPU_POTENTIAL      23
#define CPU_MESH           24
#define CPU_PEANO          25
#define CPU_COOLINGSFR     26        
#define CPU_SNAPSHOT       27
#define CPU_FOF            28
#define CPU_BLACKHOLES     29
#define CPU_MISC           30
/*   SHC   */
#define CPU_RHD            31
#define CPU_PARTS          32  /* this gives the number of parts above (must be last) */
/*   SHC   */
/* #define CPU_PARTS          31  */ /* this gives the number of parts above (must be last) */

#define CPU_STRING_LEN 120


#ifndef  TWODIMS
#define  NUMDIMS 3		/*!< For 3D-normalized kernel */
#define  KERNEL_COEFF_1  2.546479089470	/*!< Coefficients for SPH spline kernel and its derivative */
#define  KERNEL_COEFF_2  15.278874536822
#define  KERNEL_COEFF_3  45.836623610466
#define  KERNEL_COEFF_4  30.557749073644
#define  KERNEL_COEFF_5  5.092958178941
#define  KERNEL_COEFF_6  (-15.278874536822)
#define  NORM_COEFF      4.188790204786	/*!< Coefficient for kernel normalization. Note:  4.0/3 * PI = 4.188790204786 */
#else
#define  NUMDIMS 2		/*!< For 2D-normalized kernel */
#define  KERNEL_COEFF_1  (5.0/7*2.546479089470)	/*!< Coefficients for SPH spline kernel and its derivative */
#define  KERNEL_COEFF_2  (5.0/7*15.278874536822)
#define  KERNEL_COEFF_3  (5.0/7*45.836623610466)
#define  KERNEL_COEFF_4  (5.0/7*30.557749073644)
#define  KERNEL_COEFF_5  (5.0/7*5.092958178941)
#define  KERNEL_COEFF_6  (5.0/7*(-15.278874536822))
#define  NORM_COEFF      M_PI	/*!< Coefficient for kernel normalization. */
#endif


#if defined (BLACK_HOLES) || defined (DUST)
#define PPP P
#else
#define PPP SphP
#endif


#define DMAX(a,b) (dmax1=(a),dmax2=(b),(dmax1>dmax2)?dmax1:dmax2)
#define DMIN(a,b) (dmin1=(a),dmin2=(b),(dmin1<dmin2)?dmin1:dmin2)
#define IMAX(a,b) (imax1=(a),imax2=(b),(imax1>imax2)?imax1:imax2)
#define IMIN(a,b) (imin1=(a),imin2=(b),(imin1<imin2)?imin1:imin2)

#ifdef PERIODIC
extern MyDouble boxSize, boxHalf;
#ifdef LONG_X
extern MyDouble boxSize_X, boxHalf_X;
#else
#define boxSize_X boxSize
#define boxHalf_X boxHalf
#endif
#ifdef LONG_Y
extern MyDouble boxSize_Y, boxHalf_Y;
#else
#define boxSize_Y boxSize
#define boxHalf_Y boxHalf
#endif
#ifdef LONG_Z
extern MyDouble boxSize_Z, boxHalf_Z;
#else
#define boxSize_Z boxSize
#define boxHalf_Z boxHalf
#endif
#endif

#ifdef PERIODIC
#define NGB_PERIODIC_LONG_X(x) (xtmp=fabs(x),(xtmp>boxHalf_X)?(boxSize_X-xtmp):xtmp)
#define NGB_PERIODIC_LONG_Y(x) (xtmp=fabs(x),(xtmp>boxHalf_Y)?(boxSize_Y-xtmp):xtmp)
#define NGB_PERIODIC_LONG_Z(x) (xtmp=fabs(x),(xtmp>boxHalf_Z)?(boxSize_Z-xtmp):xtmp)
#else
#define NGB_PERIODIC_LONG_X(x) fabs(x)
#define NGB_PERIODIC_LONG_Y(x) fabs(x)
#define NGB_PERIODIC_LONG_Z(x) fabs(x)
#endif

#define FACT1 0.366025403785	/* FACT1 = 0.5 * (sqrt(3)-1) */

#ifdef REAL_EOS
//extern double kappa1[2850],kappa2[81][1251];                /*!< SHC: Table of opacity from OPAL project */ 
//extern double kappa_P[1100],kappa_R[1000];                /*!< SHC Table of Krumholz et al opacity */
#endif

/*********************************************************/
/*  Global variables                                     */
/*********************************************************/


extern int FirstActiveParticle;
extern int *NextActiveParticle;

#ifdef OTHIN_ACCELERATOR
extern int FirstSinkParticle;
extern int *NextSinkParticle;
#endif

#ifdef VIRTUAL
extern int FirstPhoton;
extern int *NextPhoton;
extern int FirstActivePhoton;
extern int *NextActivePhoton;
#endif

extern int TimeBinCount[TIMEBINS];
extern int TimeBinCountSph[TIMEBINS];
extern int TimeBinActive[TIMEBINS];

extern int FirstInTimeBin[TIMEBINS];
extern int LastInTimeBin[TIMEBINS];
extern int *NextInTimeBin;
extern int *PrevInTimeBin;

#ifdef SFR
extern double TimeBinSfr[TIMEBINS];
#endif

#ifdef BLACK_HOLES
extern double TimeBin_BH_mass[TIMEBINS];
extern double TimeBin_BH_dynamicalmass[TIMEBINS];
extern double TimeBin_BH_Mdot[TIMEBINS];
#endif

#ifdef DUST
extern double TimeBin_Dust_Mass[TIMEBINS];
#endif


extern int ThisTask;		/*!< the number of the local processor  */
extern int NTask;		/*!< number of processors */
extern int PTask;		/*!< note: NTask = 2^PTask */

extern double CPUThisRun;	/*!< Sums CPU time of current process */

extern int NumForceUpdate;	/*!< number of active particles on local processor in current timestep  */
extern long long GlobNumForceUpdate;

extern int NumSphUpdate;	/*!< number of active SPH particles on local processor in current timestep  */

extern int MaxTopNodes;	        /*!< Maximum number of nodes in the top-level tree used for domain decomposition */

extern int RestartFlag;		/*!< taken from command line used to start code. 0 is normal start-up from
				   initial conditions, 1 is resuming a run from a set of restart files, while 2
				   marks a restart from a snapshot file. */
extern int RestartSnapNum;

extern int *Exportflag;	        /*!< Buffer used for flagging whether a particle needs to be exported to another process */
extern int *Exportnodecount;
extern int *Exportindex;

extern int *Send_offset, *Send_count, *Recv_count, *Recv_offset, *Sendcount_matrix;

extern size_t AllocatedBytes;
extern size_t FreeBytes;

extern double CPU_Step[CPU_PARTS];
extern char CPU_Symbol[CPU_PARTS];
extern char CPU_SymbolImbalance[CPU_PARTS];
extern char CPU_String[CPU_STRING_LEN + 1];

extern double WallclockTime;    /*!< This holds the last wallclock time measurement for timings measurements */

extern int Flag_FullStep;	/*!< Flag used to signal that the current step involves all particles */


extern int TreeReconstructFlag;
extern int GlobFlag;

extern int NumPart;		/*!< number of particles on the LOCAL processor */
extern int N_gas;		/*!< number of gas particles on the LOCAL processor  */
extern long long Ntype[6];	/*!< total number of particles of each type */
extern int NtypeLocal[6];	/*!< local number of particles of each type */

extern gsl_rng *random_generator;	/*!< the random number generator used */


#ifdef SFR
extern int Stars_converted;	/*!< current number of star particles in gas particle block */
 #ifdef DUST_PEBBLES_BORN
extern int Pebbles_created; 
#endif 
#endif


extern double TimeOfLastTreeConstruction;	/*!< holds what it says */

extern int *Ngblist;		/*!< Buffer to hold indices of neighbours retrieved by the neighbour search
				   routines */


extern double DomainCorner[3], DomainCenter[3], DomainLen, DomainFac;
extern int *DomainStartList, *DomainEndList;

extern double *DomainWork;
extern int *DomainCount;
extern int *DomainCountSph;
extern int *DomainTask;
extern int *DomainNodeIndex;
extern int *DomainList, DomainNumChanged;



extern peanokey *Key, *KeySorted;

extern struct topnode_data
{
  peanokey Size;
  peanokey StartKey;
  long long Count;
  MyFloat GravCost;
  int Daughter;
  int Pstart;
  int Blocks;
  int Leaf;
} *TopNodes;

extern int NTopnodes, NTopleaves;

extern double RndTable[RNDTABLE];


#ifdef SUBFIND
extern int GrNr;
extern int NumPartGroup;
#endif


/* variables for input/output , usually only used on process 0 */

extern FILE *diag_shc,*diag_shc_ene,*diag_shc_photon,*diag_shc_photon_total;		/*!< SHC */
#ifdef REAL_EOS
extern FILE *data_kappa;                /*!< SHC */
#endif
extern FILE *photonSN;		/*!< SN */

extern char ParameterFile[100];	/*!< file name of parameterfile used for starting the simulation */

extern FILE *FdInfo,		/*!< file handle for info.txt log-file. */
 *FdEnergy,			/*!< file handle for energy.txt log-file. */
 *FdTimings,			/*!< file handle for timings.txt log-file. */
 *FdBalance,			/*!< file handle for balance.txt log-file. */
 *FdCPU;			/*!< file handle for cpu.txt log-file. */

#ifdef SFR
extern FILE *FdSfr;		/*!< file handle for sfr.txt log-file. */
#endif

#ifdef RADTRANSFER
extern FILE *FdEddington;     /*!< file handle for eddington.txt log-file. */
extern FILE *FdRadtransfer;		/*!< file handle for radtransfer.txt log-file. */
#endif

#ifdef BLACK_HOLES
extern FILE *FdBlackHoles;	/*!< file handle for blackholes.txt log-file. */
extern FILE *FdBlackHolesDetails;
#endif


#ifdef FORCETEST
extern FILE *FdForceTest;	/*!< file handle for forcetest.txt log-file. */
#endif

#ifdef DARKENERGY
extern FILE *FdDE;  /*!< file handle for darkenergy.txt log-file. */
#endif

#ifdef XXLINFO
extern FILE *FdXXL;		/*!< file handle for xxl.txt log-file. */

#ifdef MAGNETIC
extern double MeanB;

#ifdef TRACEDIVB
extern double MaxDivB;
#endif
#endif
#ifdef TIME_DEP_ART_VISC
extern double MeanAlpha;
#endif
#endif

/*! table for the cosmological drift factors */
extern double DriftTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for gravitational forces */
extern double GravKickTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for hydrodynmical forces */
extern double HydroKickTable[DRIFT_TABLE_LENGTH];

extern void *CommBuffer;	/*!< points to communication buffer, which is used at a few places */

/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
extern struct global_data_all_processes
{
  long long TotNumPart;		/*!<  total particle numbers (global value) */
  long long TotN_gas;		/*!<  total gas particle number (global value) */

#ifdef BLACK_HOLES
  int TotBHs;
#endif

  int MaxPart;			/*!< This gives the maxmimum number of particles that can be stored on one
				   processor. */
  int MaxPartSph;		/*!< This gives the maxmimum number of SPH particles that can be stored on one
				   processor. */

  int ICFormat;			/*!< selects different versions of IC file-format */

  int SnapFormat;		/*!< selects different versions of snapshot file-formats */

  int DoDynamicUpdate;

  int NumFilesPerSnapshot;	/*!< number of files in multi-file snapshot dumps */
  int NumFilesWrittenInParallel;	/*!< maximum number of files that may be written simultaneously when
					   writing/reading restart-files, or when writing snapshot files */

  int BufferSize;		/*!< size of communication buffer in MB */
  int BunchSize;     	        /*!< number of particles fitting into the buffer in the parallel tree algorithm  */


  double PartAllocFactor;	/*!< in order to maintain work-load balance, the particle load will usually
				   NOT be balanced.  Each processor allocates memory for PartAllocFactor times
				   the average number of particles to allow for that */

  double TreeAllocFactor;	/*!< Each processor allocates a number of nodes which is TreeAllocFactor times
				   the maximum(!) number of particles.  Note: A typical local tree for N
				   particles needs usually about ~0.65*N nodes. */

  double TopNodeAllocFactor;	/*!< Each processor allocates a number of nodes which is TreeAllocFactor times
				   the maximum(!) number of particles.  Note: A typical local tree for N
				   particles needs usually about ~0.65*N nodes. */

#ifdef SCALARFIELD
  double ScalarBeta;
  double ScalarScreeningLength;
#endif

  /* some SPH parameters */

  int DesNumNgb;		/*!< Desired number of SPH neighbours */
#ifdef SUBFIND
  int DesLinkNgb;
  double ErrTolThetaSubfind;
#endif

  double MaxNumNgbDeviation;	/*!< Maximum allowed deviation neighbour number */

  double ArtBulkViscConst;	/*!< Sets the parameter \f$\alpha\f$ of the artificial viscosity */
  double InitGasTemp;		/*!< may be used to set the temperature in the IC's */
  double InitGasU;		/*!< the same, but converted to thermal energy per unit mass */
  double MinGasTemp;		/*!< may be used to set a floor for the gas temperature */
  double MinEgySpec;		/*!< the minimum allowed temperature expressed as energy per unit mass */
  double MeanWeight;		/*!< Constant mean molecular weight */

#ifdef INCREASE_SPH_MASS
  double DeltaTimeGrow; /* time scale for SPH particle mass growth; SN */
#endif

  /* some force counters  */

  long long TotNumOfForces;	/*!< counts total number of force computations  */

  long long NumForcesSinceLastDomainDecomp;	/*!< count particle updates since last domain decomposition */

  /* some variable for dynamic work-load adjustment based on CPU measurements */

  double Cadj_Cost;
  double Cadj_Cpu;

  /* system of units  */

  double UnitTime_in_s,		/*!< factor to convert internal time unit to seconds/h */
    UnitMass_in_g,		/*!< factor to convert internal mass unit to grams/h */
    UnitVelocity_in_cm_per_s,	/*!< factor to convert intqernal velocity unit to cm/sec */
    UnitLength_in_cm,		/*!< factor to convert internal length unit to cm/h */
    UnitPressure_in_cgs,	/*!< factor to convert internal pressure unit to cgs units (little 'h' still
				   around!) */
    UnitDensity_in_cgs,		/*!< factor to convert internal length unit to g/cm^3*h^2 */
    UnitCoolingRate_in_cgs,	/*!< factor to convert internal cooling rate to cgs units */
    UnitEnergy_in_cgs,		/*!< factor to convert internal energy to cgs units */
    UnitTime_in_Megayears,	/*!< factor to convert internal time to megayears/h */
    GravityConstantInternal,	/*!< If set to zero in the parameterfile, the internal value of the
				   gravitational constant is set to the Newtonian value based on the system of
				   units specified. Otherwise the value provided is taken as internal gravity
				   constant G. */
    G;				/*!< Gravity-constant in internal units */

  /* Cosmology */

  double Hubble;		/*!< Hubble-constant in internal units */
  double Omega0,		/*!< matter density in units of the critical density (at z=0) */
    OmegaLambda,		/*!< vaccum energy density relative to crictical density (at z=0) */
    OmegaBaryon,		/*!< baryon density in units of the critical density (at z=0) */
    HubbleParam;		/*!< little `h', i.e. Hubble constant in units of 100 km/s/Mpc.  Only needed to get absolute
				 * physical values for cooling physics
				 */

  double BoxSize;		/*!< Boxsize in case periodic boundary conditions are used */

  /* Code options */

  int ComovingIntegrationOn;	/*!< flags that comoving integration is enabled */
  int PeriodicBoundariesOn;	/*!< flags that periodic boundaries are enabled */
  int ResubmitOn;		/*!< flags that automatic resubmission of job to queue system is enabled */
  int TypeOfOpeningCriterion;	/*!< determines tree cell-opening criterion: 0 for Barnes-Hut, 1 for relative
				   criterion */
  int TypeOfTimestepCriterion;	/*!< gives type of timestep criterion (only 0 supported right now - unlike
				   gadget-1.1) */
  int OutputListOn;		/*!< flags that output times are listed in a specified file */
  int CoolingOn;		/*!< flags that cooling is enabled */
  int StarformationOn;		/*!< flags that star formation is enabled */


  /* parameters determining output frequency */

  int SnapshotFileCount;	/*!< number of snapshot that is written next */
  double TimeBetSnapshot,	/*!< simulation time interval between snapshot files */
    TimeOfFirstSnapshot,	/*!< simulation time of first snapshot files */
    CpuTimeBetRestartFile,	/*!< cpu-time between regularly generated restart files */
    TimeLastRestartFile,	/*!< cpu-time when last restart-file was written */
    TimeBetStatistics,		/*!< simulation time interval between computations of energy statistics */
    TimeLastStatistics;		/*!< simulation time when the energy statistics was computed the last time */
  int NumCurrentTiStep;		/*!< counts the number of system steps taken up to this point */
  /* SHC */
  double FractionOfPhotons;     /*!< Fraction of photon dump */
  int NumLastDumpStep;          /*!< counts the number of system steps taken up to this point */
  /* SHC */

  /* Current time of the simulation, global step, and end of simulation */

  double Time,			/*!< current time of the simulation */
    TimeBegin,			/*!< time of initial conditions of the simulation */
    TimeStep,			/*!< difference between current times of previous and current timestep */
    TimeMax;			/*!< marks the point of time until the simulation is to be evolved */

  /* variables for organizing discrete timeline */

  double Timebase_interval;	/*!< factor to convert from floating point time interval to integer timeline */
  int Ti_Current;		/*!< current time on integer timeline */
  int Ti_nextoutput;		/*!< next output time on integer timeline */

#ifdef FLEXSTEPS
  int PresentMinStep;		/*!< If FLEXSTEPS is used, particle timesteps are chosen as multiples of the present minimum timestep. */
  int PresentMaxStep;		/*!< If FLEXSTEPS is used, this is the maximum timestep in timeline units, rounded down to the next power 2 division */
#endif


#ifdef PMGRID
  int PM_Ti_endstep, PM_Ti_begstep;
  double Asmth[2], Rcut[2];
  double Corner[2][3], UpperCorner[2][3], Xmintot[2][3], Xmaxtot[2][3];
  double TotalMeshSize[2];
#endif

#ifdef CHEMISTRY
  double Epsilon;
#endif

  int Ti_nextlineofsight;
#ifdef OUTPUTLINEOFSIGHT
  double TimeFirstLineOfSight;
#endif

  /* variables that keep track of cumulative CPU consumption */

  double TimeLimitCPU;
  double CPU_Sum[CPU_PARTS];    /*!< sums wallclock time/CPU consumption in whole run */

  /* tree code opening criterion */

  double ErrTolTheta;		/*!< BH tree opening angle */
  double ErrTolForceAcc;	/*!< parameter for relative opening criterion in tree walk */


  /* adjusts accuracy of time-integration */

  double ErrTolIntAccuracy;	/*!< accuracy tolerance parameter \f$ \eta \f$ for timestep criterion. The
				   timesteps is \f$ \Delta t = \sqrt{\frac{2 \eta eps}{a}} \f$ */

  double MinSizeTimestep,	/*!< minimum allowed timestep. Normally, the simulation terminates if the
				   timestep determined by the timestep criteria falls below this limit. */
    MaxSizeTimestep;		/*!< maximum allowed timestep */

  double MaxRMSDisplacementFac;	/*!< this determines a global timestep criterion for cosmological simulations
				   in comoving coordinates.  To this end, the code computes the rms velocity
				   of all particles, and limits the timestep such that the rms displacement
				   is a fraction of the mean particle separation (determined from the
				   particle mass and the cosmological parameters). This parameter specifies
				   this fraction. */



  double CourantFac;		/*!< SPH-Courant factor */


  /* frequency of tree reconstruction/domain decomposition */


  double TreeDomainUpdateFrequency;	/*!< controls frequency of domain decompositions  */


  /* gravitational and hydrodynamical softening lengths (given in terms of an `equivalent' Plummer softening
   * length)
   *
   * five groups of particles are supported 0=gas,1=halo,2=disk,3=bulge,4=stars
   */
  double MinGasHsmlFractional,	/*!< minimum allowed SPH smoothing length in units of SPH gravitational
				   softening length */
    MinGasHsml;			/*!< minimum allowed SPH smoothing length */


  double SofteningGas,		/*!< for type 0 */
    SofteningHalo,		/*!< for type 1 */
    SofteningDisk,		/*!< for type 2 */
    SofteningBulge,		/*!< for type 3 */
    SofteningStars,		/*!< for type 4 */
    SofteningBndry;		/*!< for type 5 */

  double SofteningGasMaxPhys,	/*!< for type 0 */
    SofteningHaloMaxPhys,	/*!< for type 1 */
    SofteningDiskMaxPhys,	/*!< for type 2 */
    SofteningBulgeMaxPhys,	/*!< for type 3 */
    SofteningStarsMaxPhys,	/*!< for type 4 */
    SofteningBndryMaxPhys;	/*!< for type 5 */

  double SofteningTable[6];	/*!< current (comoving) gravitational softening lengths for each particle type */
  double ForceSoftening[6];	/*!< the same, but multiplied by a factor 2.8 - at that scale the force is Newtonian */


  /*! If particle masses are all equal for one type, the corresponding entry in MassTable is set to this
   *  value, * allowing the size of the snapshot files to be reduced
   */
  double MassTable[6];


  /* some filenames */
  char InitCondFile[100],
    OutputDir[100],
    SnapshotFileBase[100],
    EnergyFile[100],
    CpuFile[100],
    InfoFile[100], TimingsFile[100], RestartFile[100], ResubmitCommand[100], OutputListFilename[100];

  /*! table with desired output times */
  double OutputListTimes[MAXLEN_OUTPUTLIST];

  int OutputListLength;		/*!< number of times stored in table of desired output times */



#if defined(ADAPTIVE_GRAVSOFT_FORGAS) && !defined(ADAPTIVE_GRAVSOFT_FORGAS_HSML)
  double ReferenceGasMass;
#endif

#ifdef SGRA_POTENTIAL
    double CuspMass;
    double CuspRadius,CuspAlpha1, CuspAlpha2;
    double xbh, ybh, zbh, mbh; /* position of the SMBH */
#endif

#ifdef NFW_POTENTIAL 
    double VirialMass;   /* Virial Mass of Halo */
    double VirialRadius; /* Virial Radius of Halo */
    double ScaleRadius;  /* Scale Radius of Halo */
    double CNFW;         /* NFW concentration */
    double DeltaVir;     /* Virial overdensity criterion */
    double xbh, ybh, zbh, mbh; /* position of the SMBH */
    double FlatRadius;   /* Flattening radius */
    double BH_Luminosity; /* position of the SMBH */
#endif

#if defined(OTHIN_ACCELERATOR) || defined(SIS_POTENTIAL) || defined(CUSP_POTENTIAL) || defined(QUASAR_HEATING) || defined(FIND_SMBH)
    double CuspMass;
    double CuspRadius,CuspGamma,FlatRadius;
    double xbh, ybh, zbh, mbh; /* position of the SMBH */
    double BH_Luminosity; /* position of the SMBH */
    double BH_LuminosityFactor; /* A multiplicative factor introduced in the input file */
#endif

#if defined(PLANET_IRRADIATION) || defined(PLANET_ACCRETION_FEEDBACK) || defined(PLANET_ACCRETION_REPORTING)
/* SN: this is for preheating by the planet */
  double xpl, ypl, zpl, mpl, rad_pl, planet_luminosity; /* position of the planet */
#endif


#ifdef SFR		/* star formation and feedback sector */
  double CritOverDensity;
  double CritPhysDensity;
  double OverDensThresh;
  double PhysDensThresh;
  double EgySpecSN;
  double FactorSN;
  double EgySpecCold;
  double FactorEVP;
  double FeedbackEnergy;
  double TempSupernova;
  double TempClouds;
  double MaxSfrTimescale;
  double WindEfficiency;
  double WindEnergyFraction;
  double WindFreeTravelLength;
  double WindFreeTravelDensFac;
  double FactorForSofterEQS;
#endif

#ifdef DARKENERGY
  double DarkEnergyParam;	/*!< fixed w for equation of state */
#ifdef TIMEDEPDE
  char DarkEnergyFile[100];	/*!< tabelized w for equation of state */
#ifdef TIMEDEPGRAV
  double Gini;
#endif
#endif
#endif

#ifdef RESCALEVINI
  double VelIniScale;		/*!< Scale the initial velocities by this amount */
#endif

#ifdef TIME_DEP_ART_VISC
  double ViscSource0;		/*!< Given sourceterm in viscosity evolution */
  double DecayLength;		/*!< Number of h for the viscosity decay */
  double ViscSource;		/*!< Reduced sourceterm in viscosity evolution */
  double DecayTime;		/*!< Calculated decaytimescale */
  double AlphaMin;		/*!< Minimum of allowed viscosity parameter */
#endif

#ifdef CONDUCTION
  double ConductionCoeff;	/*!< Thermal Conductivity */
#ifdef CONDUCTION_SATURATION
  double ElectronFreePathFactor;	/*!< Factor to get electron mean free path */
#endif
#endif

#ifdef BINISET
  double BiniX, BiniY, BiniZ;	/*!< Initial values for B */
#endif

#ifdef BSMOOTH
  int BSmoothInt;
  double BSmoothFrac;
  int MainTimestepCounts;
#endif

#ifdef MAGNETIC_DISSIPATION
  double ArtMagDispConst;	/*!< Sets the parameter \f$\alpha\f$ of the artificial magnetic disipation */
#ifdef TIME_DEP_MAGN_DISP
  double ArtMagDispMin;
  double ArtMagDispSource;
  double ArtMagDispTime;
#endif
#endif

#ifdef DIVBCLEANING_DEDNER
  double DivBcleanParabolicSigma;
  double DivBcleanHyperbolicSigma;
#endif

#ifdef BLACK_HOLES
  double TimeNextBlackHoleCheck;
  double TimeBetBlackHoleSearch;
  double BlackHoleAccretionFactor;	/*!< Fraction of BH bondi accretion rate */
  double BlackHoleFeedbackFactor;	/*!< Fraction of the black luminosity feed into thermal feedback */
  double SeedBlackHoleMass;	/*!< Seed black hole mass */
  double MinFoFMassForNewSeed;	/*!< Halo mass required before new seed is put in */
  double BlackHoleNgbFactor;	/*!< Factor by which the normal SPH neighbour should be increased/decreased */
	double BlackHoleActiveTime;
	double SMBHmass;
  double BlackHoleEddingtonFactor;	/*! Factor above Eddington */

#ifdef ACCRETION_RADIUS
    double InnerBoundary;
// Addition by CBP (30/10/2008)
    double AccDtBlackHole;       /* Regulates size of black hole timestep 
				    -- fraction of orbital time at
				    Inner Boundary */
    double ViscousTimescale;     /* Viscous timescale of accretion disc 
				    surrounding black hole -- regulates
				    the black hole accretion rate mdot */
    double AccDiscSFR;           /* Star formation rate in the accretion 
				    disc */
// End of CBP
#ifdef  BH_FORM
    double SinkBoundary;
#endif
#endif

// Addition by CBP (10/10/2008)
#ifdef VARIABLE_PHOTON_SEARCH_RADIUS
    double PhotonSearchRadius;   /* Search radius for photons, used when
				    this is allowed to vary with distance 
				    from the source. */
#endif
// End of CBP (10/10/2008)

#if defined(BETA_COOLING) || defined(EVAPORATION) || defined(EVAPORATION_RADIAL)
    double BetaCool;             /* Beta parameter for cooling model */
#endif

#if defined(EVAPORATION) || defined(EVAPORATION_RADIAL)
  double Evap_dens; /* SN: Density g/cm^3 below which extra heating applied to evaporate gas clouds */
#endif

#if defined(EVAPORATION_RADIAL)
  double Cool_ind; /* RJH: Index for irradiation temp floor */  
  double rho_cool_ind; /* RJH: Index for density cooling switch */

#endif

#if defined(BACKGROUND_ILLUMINATION_TEMP) || defined(BETA_COOLING) || defined(ISOTHERM) || defined(EVAPORATION) || defined(EVAPORATION_RADIAL)
    double EqTemp;               /* Equilibrium temperature; SN -- also used
				  * as background illumination temperature in
				  * fb_partciles.c */				    
#endif

#ifdef POLYTROP  //RJH 20/03/2018
  double Poly_K; /*Polytropic K constant for cooling. In cgs*/
#endif

// Addition by CBP (06/11/2008)
#if defined(BH_KICK)
    double BlackHoleMass, BlackHoleKickVel, BlackHoleKickAngle;
#endif
// End of CBP

#if defined(PLANET_INITIAL_MASS)
  /* SN: redefine mass of the secondary BH using the mass in the input                                                                                                                                                                        
     file */
  double BlackHoleMass;
#endif


#ifdef VIRTUAL
    double VirtualStart;
    double FeedBackVelocity;
    double VirtualMomentum;
    double VirtualMass;
    double VirtualTime;
    double VirtualCrosSection;
    double VirtualFeedBack;
    double OuterBoundary;
    double FreeTravelLength;
#ifdef DUST_GROWTH
/* #ifndef DUST_VAPORIZE  */
/*   double FragmentationVelocity; */
/*   double FragmentationDensity, PebbleBirthTime, PebbleSizeGrTime, PebbleVelBreak, PebbleBirthMass; */
/* #endif */
  double FragmentationVelocity;
#endif
#endif

#ifdef VIRTUAL_SET_MASS
    double OriginalGasMass;
#endif


#ifdef FOF
  double massDMpart;
#endif
#ifdef MODIFIEDBONDI
  double BlackHoleRefDensity;
  double BlackHoleRefSoundspeed;
#endif
#endif

#ifdef COSMIC_RAYS
  double CR_Alpha;		/*!< Cosmic ray spectral index [2..3] */
  double CR_SNEff;		/*!< SN injection efficiency [0..1] */
  double CR_SNAlpha;		/*!< SN injection spectral index [2..3] */
  int bDebugFlag;		/*!< enables debug outputs after triggered */

#if defined(CR_DIFFUSION) || defined (CR_DIFFUSION_GREEN)
  double CR_DiffusionCoeff;	/*!< (temporary) fixed value for CR diffusivity */

  double CR_DiffusionDensScaling;	/*!< grade of density dependence of diffusivity */
  double CR_DiffusionDensZero;	/*!< Reference point density for diffusivity */

  double CR_DiffusionEntropyScaling;	/*!< grade of specific energy dependence of diffusivity */

  double CR_DiffusionEntropyZero;	/*!< Reference Entropic function for diffusivity */

  double CR_DiffusionTimeScale;	/*!< Parameter for Diffusion Timestep Criterion */

  double TimeOfLastDiffusion;
#endif				/* CR_DIFFUSION */

#if defined(CR_SHOCK)
#if (CR_SHOCK == 1)
  double CR_ShockAlpha;		/*!< spectral index to be used in shock injection */
#else
  double CR_ShockCutoff;	/*!< Cutoff factor x_inj for CR accel */
#endif
  double CR_ShockEfficiency;	/*!< energy fraction of shock energy fed into CR */
#endif				/* CR_SHOCK */

#ifdef FIX_QINJ
  double Shock_Fix_Qinj;	/*!< inject only CRps with threshold cutoff Shock_Fix_Qinj */
#endif

#endif				/* COSMIC_RAYS */

#ifdef MACHNUM
  double Shock_Length;		/*!< length scale on which the shock is smoothed out */
  double Shock_DeltaDecayTimeMax;	/*!< maximum time interval (Dloga) for which the 
					   Mach number is kept at its maximum */
#endif

#ifdef REIONIZATION
  int not_yet_reionized;	/*!< flag that makes sure that there is only one reionization */
#endif



#ifdef BUBBLES
  double BubbleDistance;
  double BubbleRadius;
  double BubbleTimeInterval;
  double BubbleEnergy;
  double TimeOfNextBubble;
  double FirstBubbleRedshift;
#ifdef FOF
  int BiggestGroupLen;
  float BiggestGroupCM[3];
  double BiggestGroupMass;
#endif
#endif

#if defined(MULTI_BUBBLES) && defined(FOF)
#ifndef BLACK_HOLES
  double MinFoFMassForNewSeed;	/*!< Halo mass required before new seed is put in */
  double massDMpart;
#endif
  double BubbleDistance;
  double BubbleRadius;
  double BubbleTimeInterval;
  double BubbleEnergy;
  double TimeOfNextBubble;
  double ClusterMass200;
  double FirstBubbleRedshift;
#endif

#ifdef NAVIERSTOKES
  double NavierStokes_ShearViscosity;
  double FractionSpitzerViscosity;
  double ShearViscosityTemperature;
#endif
#ifdef NAVIERSTOKES_BULK
  double NavierStokes_BulkViscosity;
#endif
#ifdef VISCOSITY_SATURATION
  double IonMeanFreePath;
#endif

#ifdef RADTRANSFER
  MyFloat kappaMean;
  MyFloat JMean;
  MyFloat RadIMean;
  MyDouble Residue;
  float TotVol;
#endif

#ifdef DUST
  double InitialDustRadius;
  double InitialDustMass;
#ifdef DUST_POWERLAW
  double aDustMin0;
  double aDustMax0;
#endif
#ifdef DUST_TWO_POPULATIONS 
  double InitialMicroDustZ;
#endif
  double total_heating;
  double total_energy_t0;
  double total_delta_energy;
#ifdef DUST_CONSTANT_TSTOP
  double DustTstop;
#endif
#ifdef DUST_PEBBLES_BORN
  double DustInjectionTime, PebbleInjectionMass, PebbleInjectionDens, PebbleInjectionSlab;
#endif
#ifdef DUST_SINK_ON_FLY
  double CritDustDensity;
#endif
#endif
}
All;




/*! This structure holds all the information that is
 * stored for each particle of the simulation.
 */
extern struct particle_data
{
  MyDouble Pos[3];   /*!< particle position at its current time */
  MyDouble Vel[3];   /*!< particle velocity at its current time */
  MyDouble Mass;     /*!< particle mass */
  MyIDType ID;

  union
  {
    MyFloat       GravAccel[3];		/*!< particle acceleration due to gravity */
    MyLongDouble dGravAccel[3];
  } g;
#ifdef PMGRID
  MyFloat GravPM[3];		/*!< particle acceleration due to long-range PM gravity force */
#endif
#ifdef FORCETEST
  MyFloat GravAccelDirect[3];	/*!< particle acceleration calculated by direct summation */
#endif
#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY) || defined(OUTPUTPOTENTIAL)
  union
  {
    MyFloat       Potential;		/*!< gravitational potential */
    MyLongDouble dPotential;
  } p;
#endif
#if defined(DISTORTIONTENSOR) || defined(OUTPUT_TIDALTENSOR)
  union
  {
    MyFloat       tidal_tensor[6];		/*!< tidal tensor (=second derivatives of grav. potential)*/ 
    MyLongDouble dtidal_tensor[6]; 
  } tite;
#endif
  MyFloat OldAcc;			/*!< magnitude of old gravitational force. Used in relative opening
                                          criterion */
#ifdef DISTORTIONTENSOR
   MyDouble distortion_tensor[9];       /*!< Distortion tensor for that particle*/
   MyDouble distortion_tensor_vel[9];   /*!< Distortion tensor 'velocity' for that particle*/   
#endif

#if defined(EVALPOTENTIAL) && defined(PMGRID)
  MyFloat PM_Potential;
#endif

#ifdef STELLARAGE
  MyFloat StellarAge;		/*!< formation time of star particle */
#endif
#ifdef METALS
  MyFloat Metallicity;		/*!< metallicity of gas or star particle */
#endif				/* closes METALS */

#if defined (BLACK_HOLES) || defined (DUST)
  MyFloat Hsml;

  union
  {
    MyFloat       NumNgb;
    MyLongDouble dNumNgb;
  } n;
#endif


#ifdef BLACK_HOLES
  int SwallowID;
  MyFloat BH_Mass;
  MyFloat BH_Mdot;
// CBP (31/10/2008) -- Addition of Accretion Disc Mass    
  MyFloat AccDisc_Mass;    
  MyFloat BH_Mseed;
// End of CBP
#ifdef DUST
  MyFloat DeltaDustMomentum[3]; /* Used to hold the new value of density for virtual particles */
  MyFloat NewDragAcc[3]; /* Used to hold the new value of density for virtual particles */
  MyFloat DeltaDragEnergy;
  MyFloat DustRadius; /* SN: Dust Radius for a pebble */
  MyFloat DustVcoll; /* SN: dust-gas velocity difference a pebble */
  MyFloat Dust_Mass, Total_Mass; /* SN: dust and total (dust+gas) mass accreted by a sink */
  MyFloat LogDustRadius_by_dt; /* SN: Time derivative of Dust radius defined on a dust particle */
#ifdef DUST_TWO_POPULATIONS
  /* SN -- in this case Dust Radius is ALSO calculated at SPH particle position */
  MyFloat MicroDustMass; /* Mass of micron sized grains tightly bound to the SPH particle */
  MyFloat d_MicroDustMass; /* addition to the above in a time step */
  MyFloat NewKappa, Dlog_a_dt; /* SN: Dust Radius, micron size dust density at SPH particle position */
  //#endif
#endif
#endif
#ifdef VIRTUAL
    int WindOrPhoton; /* A Flag to mark whether the feedback particle is a
		       * wind or photon particle. The former only transfer
		       * momentum and have a huge absorption cross section,
		       * whereas the latter are the usual photons -- SN March
		       * 09 */
    MyFloat NewDensity; /* Used to hold the new value of density for virtual particles */
    MyFloat dt1, dt2, DtJump;
    MyFloat DeltaPhotonMomentum;
    MyFloat OldPhotonMomentum;
    MyFloat BirthPhotonEnergy;
#ifdef NON_GRAY_RAD_TRANSFER
    MyFloat Kappa; /* individual opacity for photons */
    MyFloat VirtualTemperature; /* Temperature of gas neighbors of the photon */
#endif
/* #ifdef REAL_EOS */
/*     MyFloat NewKappa; /\* Used to hold the new value of density for virtual particles *\/ */
/* #endif */

#if defined(VIRTUAL_HEATING) 
    MyFloat DeltaPhotonEnergy;
    MyFloat DeltaHeat; /* energy transferred from photon to sph */
#endif

#endif /*  VIRTUAL */
  union
  {
    MyFloat BH_Density;
    MyLongDouble dBH_Density;
  } b1;
  union
  {
    MyFloat BH_Entropy;
    MyLongDouble dBH_Entropy;
  } b2;
  union
  {
    MyFloat BH_SurroundingGasVel[3];
    MyLongDouble dBH_SurroundingGasVel[3];
  } b3;
  union
  {
    MyFloat BH_accreted_Mass;
    MyLongDouble dBH_accreted_Mass;
  } b4;
  union
  {
    MyFloat BH_accreted_BHMass;
    MyLongDouble dBH_accreted_BHMass;
    MyFloat BH_accreted_DustMass;
    MyLongDouble dBH_accreted_DustMass;
  } b5;
  union
  {
    MyFloat BH_accreted_momentum[3];
    MyLongDouble dBH_accreted_momentum[3];
  } b6;
#if defined (DUST)
/*
    MyFloat HsmlDust;
    int    NumDustNgb;
  union
  {
  } nd;
*/
  union
  {
    MyFloat DUST_Density;
    MyLongDouble dDUST_Density;
  } d1;
  union
  {
    MyFloat DUST_Entropy;
    MyLongDouble dDUST_Entropy;
  } d2;
  union
  {
    MyFloat DUST_SurroundingGasVel[3];
    MyLongDouble dDUST_SurroundingGasVel[3];
  } d3;
  union
  {
    MyFloat DUST_accreted_Mass;
    MyLongDouble dDUST_accreted_Mass;
  } d4;
  union
  {
    MyFloat DUST_accreted_DUSTMass;
    MyLongDouble dDUST_accreted_DUSTMass;
  } d5;
  union
  {
    MyFloat DUST_accreted_momentum[3];
    MyLongDouble dDUST_accreted_momentum[3];
  } d6;
  union
  {
    MyFloat DUST_particle_density;
    MyLongDouble dDUST_particle_density;
  } d7;
  union
  {
    MyFloat DUST_particle_velocity[3];
    MyLongDouble dDUST_particle_velocity[3];
  } d9;
#ifdef DUST_TWO_POPULATIONS
  union
  {
    MyFloat DUST_micro_Density;
    MyLongDouble dDUST_micro_Density;
  } d8;
#endif
//  union
//  {
//    MyFloat       DragAccel[3];		/*!< acceleration due to radiation pressure */
//    MyLongDouble dDragAccel[3];
//  } da;
#endif
#ifdef REPOSITION_ON_POTMIN
  MyFloat BH_MinPotPos[3];
  MyFloat BH_MinPot;
#endif
#ifdef BH_KINETICFEEDBACK
  MyFloat ActiveTime;
  MyFloat ActiveEnergy;
#endif
#endif

#ifdef SUBFIND
  int GrNr;
  int SubNr;
  int DM_NumNgb;
  int targettask, origintask, submark;
  MyFloat DM_Hsml;
  union
  {
    MyFloat DM_Density;	
    MyFloat DM_Potential;
  } u;
  union
  {
    MyFloat DM_VelDisp;
    MyFloat DM_BindingEnergy;
  } v;
#endif

#if defined(ORDER_SNAPSHOTS_BY_ID) && !defined(SUBFIND)
  int     GrNr;
  int     SubNr;
#endif


  float GravCost;		/*!< weight factor used for balancing the work-load */

  int Ti_begstep;		/*!< marks start of current timestep of particle on integer timeline */
  int Ti_current;		/*!< current time of the particle */

  short int Type;		/*!< flags particle type.  0=gas, 1=halo, 2=disk, 3=bulge, 4=stars, 5=bndry */
  short int TimeBin;
}
 *P,				/*!< holds particle data on local processor */
 *DomainPartBuf;		/*!< buffer for particle data used in domain decomposition */



/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 * variables.
 */
extern struct sph_particle_data
{
  MyDouble Entropy;		/*!< current value of entropy (actually entropic function) of particle */
  MyFloat  Pressure;		/*!< current pressure */
  MyFloat  VelPred[3];		/*!< predicted SPH particle velocity at the current time */
#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
  MyFloat MinViscousDt;
#else
  MyFloat MaxSignalVel;           /*!< maximum signal velocity */
#endif
#ifdef BLACK_HOLES
#ifdef VIRTUAL
    MyFloat accreted_mass;
    MyFloat accreted_momentum[3];
#endif
#endif
#ifdef REAL_EOS
    MyDouble Temperature;       /*!< Temperature */
    MyDouble Mu;                /*!< Mean Molecular Weight */
    MyDouble Kappa;             /*!< Opacity */
#endif

  union
  {
    MyFloat       Density;		/*!< current baryonic mass density of particle */
    MyLongDouble dDensity;
  } d;
  union
  {
    MyFloat       DtEntropy;		/*!< rate of change of entropy */
    MyLongDouble dDtEntropy;
  } e;
  union
  {
    MyFloat       HydroAccel[3];	/*!< acceleration due to hydrodynamical force */
    MyLongDouble dHydroAccel[3];
  } a;
#ifdef RAD_ACCEL /* SHC */
  union
  {
    MyFloat       RadAccel[3];		/*!< acceleration due to radiation pressure */
    MyLongDouble dRadAccel[3];
#ifdef SMOOTHED_RAD_ACCEL /* SHC */
    MyFloat       SmRadAccel[3];	/*!< acceleration due to radiation pressure */
    MyLongDouble dSmRadAccel[3];
#endif  
} ra;
#endif
#ifdef DUST /* SHC */
  union
  {
    MyFloat       DragAccel[3];		/*!< acceleration due to radiation pressure */
    MyLongDouble dDragAccel[3];
  } da;
  union
  { 
    MyFloat       DragHeating;		/* Heating of SPH particles due to dust drag */
    MyLongDouble dDragHeating;
  } dh;
#endif
  union
  {
    MyFloat       DhsmlDensityFactor;	/*!< correction factor needed in entropy formulation of SPH */
    MyLongDouble dDhsmlDensityFactor;
  } h;
  union
  {
    MyFloat       DivVel;		/*!< local velocity divergence */
    MyLongDouble dDivVel;
  } v;
#ifndef NAVIERSTOKES
  union
  {
    MyFloat CurlVel;     	        /*!< local velocity curl */
    MyFloat       Rot[3];		/*!< local velocity curl */
    MyLongDouble dRot[3];
  } r;
#else
  union
  {
    MyFloat DV[3][3];
    struct
    {
      MyFloat DivVel;
      MyFloat CurlVel;
      MyFloat StressDiag[3];
      MyFloat StressOffDiag[3];
#ifdef NAVIERSTOKES_BULK
      MyFloat StressBulk;
#endif
    } s;
  } u;
#endif

#if !defined(BLACK_HOLES) || defined (DUST)
  MyFloat Hsml;			/*!< current smoothing length */
  MyFloat Left,			/*!< lower bound in iterative smoothing length search */
          Right;		/*!< upper bound in iterative smoothing length search */
  union
  {
    MyFloat       NumNgb;
    MyLongDouble dNumNgb;
  } n;
#endif

#ifdef VIRTUAL_HEATING
    MyFloat       Injected_VIRTUAL_Energy, Emitted_VIRTUAL_Energy, Returned_E_Fraction;
#ifdef DECREASE_POISSON_NOISE 
    MyFloat       Injected_VIRTUAL_Energy_Old;
#endif

#endif


#if defined(BH_THERMALFEEDBACK) || defined(BH_KINETICFEEDBACK) || defined(VIRTUAL)
  union
  {
    MyFloat       Injected_BH_Energy;
    MyLongDouble dInjected_BH_Energy;
    MyFloat       Injected_BH_Momentum[3];
    MyLongDouble dInjected_BH_Momentum[3];
  } i;
#endif

#ifdef COOLING
  MyFloat Ne;  /*!< electron fraction, expressed as local electron number
		    density normalized to the hydrogen number density. Gives
		    indirectly ionization state and mean molecular weight. */
#endif
#ifdef SFR
  MyFloat Sfr;
#endif
#ifdef WINDS
  MyFloat DelayTime;		/*!< remaining maximum decoupling time of wind particle */
#endif


#ifdef MAGNETIC
  MyFloat B[3], BPred[3];
  MyFloat DtB[3];
#if defined(TRACEDIVB) || defined(TIME_DEP_MAGN_DISP)
  MyFloat divB;
#endif
#if defined(BSMOOTH) || defined(BFROMROTA)
  MyFloat BSmooth[3];
  MyFloat DensityNorm;
#endif
#ifdef TIME_DEP_MAGN_DISP
  MyFloat Balpha, DtBalpha;
#endif
#ifdef DIVBCLEANING_DEDNER
  MyFloat Phi, PhiPred, DtPhi;
#ifdef SMOOTH_PHI
  MyFloat SmoothPhi;
#endif
#endif
#endif
#ifdef TIME_DEP_ART_VISC
  MyFloat alpha, Dtalpha;
#endif
#ifdef NS_TIMESTEP
  MyFloat ViscEntropyChange;
#endif
#ifdef CONDUCTION
  MyFloat CondEnergyChange;
  MyFloat SmoothedEntr;
#ifdef CONDUCTION_SATURATION
  MyFloat GradEntr[3];
#endif
#ifdef OUTPUTCOOLRATE
  MyFloat CondRate;
#endif
#endif


#ifdef MHM
  MyFloat FeedbackEnergy;
#endif

#ifdef COSMIC_RAYS
  MyFloat CR_C0;			/*!< Cosmic ray amplitude adiabatic invariable */
  MyFloat CR_q0;			/*!< Cosmic ray cutoff adiabatic invariable */
  MyFloat CR_E0;			/*!< Specific Energy at Rho0 */
  MyFloat CR_n0;			/*!< baryon fraction in cosmic rays */

  MyFloat CR_DeltaE;		/*!< Specific Energy growth during timestep */
  MyFloat CR_DeltaN;		/*!< baryon fraction growth during timestep */
#ifdef MACHNUM
  MyFloat CR_Gamma0;
#endif

#ifdef CR_OUTPUT_INJECTION
  MyFloat CR_Specific_SupernovaHeatingRate;
#endif

#ifdef CR_DIFFUSION
  MyFloat CR_SmoothE0;		/*!< SPH-smoothed interpolant of diffusion term */
  MyFloat CR_Smoothn0;		/*!< SPH-smoothed interpolant for diffusion term */
#endif				/* CR_DIFFUSION */
#ifdef CR_DIFFUSION_GREEN
  MyFloat CR_Kappa;
  MyFloat CR_Kappa_egy;
  MyFloat CR_WeightSum;
  MyFloat CR_WeightSum_egy;
#endif
#endif				/* COSMIC_RAYS */

#ifdef MACHNUM
  MyFloat Shock_MachNumber;	/*!< Mach number */
  MyFloat Shock_DecayTime;	/*!< Shock decay time */
#ifdef COSMIC_RAYS
  MyFloat Shock_DensityJump;	/*!< Density jump at the shock */
  MyFloat Shock_EnergyJump;	/*!< Energy jump at the shock */
  MyFloat PreShock_PhysicalDensity;	/*!< Specific energy in the preshock regime */
  MyFloat PreShock_PhysicalEnergy;	/*!< Density in the preshock regime */
  MyFloat PreShock_XCR;		/*!< XCR = PCR / Pth in the preshock regime */
#endif
#ifdef MACHSTATISTIC
  MyFloat Shock_DtEnergy;		/*!< Change of thermal specific energy at Shocks */
#endif
#endif				/* Mach number estimate */


#ifdef CHEMISTRY
  MyFloat elec;
  MyFloat HI;
  MyFloat HII;

  MyFloat HeI;
  MyFloat HeII;
  MyFloat HeIII;

  MyFloat H2I;
  MyFloat H2II;

  MyFloat HM;

  MyFloat Gamma;
  MyFloat t_elec, t_cool;
#endif

#ifdef RADTRANSFER
  MyFloat ET[6];                /* eddington tensor - symmetric -> only 6 elements needed */
  MyFloat RadJ;                 /* mean intensity */
  MyFloat P1, D1;               /* coefficients in GAuss-Seidel method */
  MyFloat Je;                   /* emissivity */
  MyFloat nH;
  MyFloat x;
  MyFloat nHI;                  /* HI number density */
  MyFloat nHII;                 /* HII number density */
  MyFloat n_elec;               /* electron number density */
#endif

}
  *SphP,				/*!< holds SPH particle data on local processor */
  *DomainSphBuf;			/*!< buffer for SPH particle data in domain decomposition */


extern peanokey *DomainKeyBuf;

/* global state of system 
*/
extern struct state_of_system
{
  double Mass,
/* SHC BEGIN */
    EnergyRadComp,EnergyRadAdded,EnergyRadDeleted,
#ifdef CHECK_ENERGY_CONSERVATION
    EnergySinkInt,
#endif
/* SHC END */
    EnergyKin,
    EnergyPot,
    EnergyInt,
    EnergyTot,
    Momentum[4],
    AngMomentum[4],
    CenterOfMass[4],
    MassComp[6],
    EnergyKinComp[6],
    EnergyPotComp[6],
    EnergyIntComp[6], EnergyTotComp[6], MomentumComp[6][4], AngMomentumComp[6][4], CenterOfMassComp[6][4];
}
SysState, SysStateAtStart, SysStateAtEnd;


/* Various structures for communication during the gravity computation.
 */

extern struct data_index
{
  int Task;
  int Index;
  int IndexGet;
}
 *DataIndexTable;		/*!< the particles to be exported are grouped
                                  by task-number. This table allows the
                                  results to be disentangled again and to be
                                  assigned to the correct particle */

extern struct data_nodelist
{
  int NodeList[NODELISTLENGTH];
}
*DataNodeList;		

extern struct gravdata_in
{
  MyFloat Pos[3];
#ifdef UNEQUALSOFTENINGS
  int Type;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  MyFloat Soft;
#endif
#endif
  MyFloat OldAcc;
  int NodeList[NODELISTLENGTH];
}
 *GravDataIn,			/*!< holds particle data to be exported to other processors */
 *GravDataGet;			/*!< holds particle data imported from other processors */


extern struct gravdata_out
{
  MyLongDouble Acc[3];
#ifdef EVALPOTENTIAL
  MyLongDouble Potential;
#endif
#if defined(DISTORTIONTENSOR) || defined(OUTPUT_TIDALTENSOR)
  MyLongDouble tidal_tensor[6];
#endif
  int Ninteractions;
}
 *GravDataResult,		/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *GravDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */


extern struct potdata_out
{
  MyLongDouble Potential;
}
 *PotDataResult,		/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *PotDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */




/*! Header for the standard file format.
 */
extern struct io_header
{
  int npart[6];			/*!< number of particles of each type in this file */
  double mass[6];		/*!< mass of particles of each type. If 0, then the masses are explicitly
				   stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;			/*!< time of snapshot file */
  double redshift;		/*!< redshift of snapshot file */
  int flag_sfr;			/*!< flags whether the simulation was including star formation */
  int flag_feedback;		/*!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[6];	/*!< total number of particles of each type in this snapshot. This can be
				   different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;		/*!< flags whether cooling was included  */
  int num_files;		/*!< number of files in multi-file snapshot */
  double BoxSize;		/*!< box-size of simulation in case periodic boundaries were used */
  double Omega0;		/*!< matter density in units of critical density */
  double OmegaLambda;		/*!< cosmological constant parameter */
  double HubbleParam;		/*!< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;		/*!< flags whether the file contains formation times of star particles */
  int flag_metals;		/*!< flags whether the file contains metallicity values for gas and star
				   particles */
  unsigned int npartTotalHighWord[6];	/*!< High word of the total number of particles of each type */
  int flag_entropy_instead_u;	/*!< flags that IC-file contains entropy instead of u */
  int flag_doubleprecision;	/*!< flags that snapshot contains double-precision instead of single precision */
  char fill[56];		/*!< fills to 256 Bytes */
}
header;				/*!< holds header for snapshot files */



enum iofields
{ IO_POS,
  IO_VEL,
  IO_ID,
  IO_MASS,
  IO_U,
  IO_RHO,
  IO_NE,
  IO_NH,
  IO_HSML,
  IO_DUST,
  IO_TWO_DUST,
  IO_VCOLL_DUST,
  IO_TGROW_DUST,
  IO_SFR,
  IO_AGE,
  IO_Z,
  IO_BHMASS,
  IO_BHMDOT,
  IO_ACCDISCMASS,
  IO_BHMSEED,
  IO_NEWDENS,
  IO_DELPHOTONMOM,
  IO_OLDPHOTONMOM,
  IO_BIPHOTONEN,
  IO_DELPHOTONEN,
  IO_ACCMASS,
  IO_ACCMOM,
  IO_INJVIRTUALEN,
  IO_EMVIRTUALEN,
  IO_RET_E_FRAC,
  IO_POT,
  IO_ACCEL,
  IO_CR_C0,
  IO_CR_Q0,
  IO_CR_P0,
  IO_CR_E0,
  IO_CR_n0,
  IO_CR_ThermalizationTime,
  IO_CR_DissipationTime,
  IO_ELECT,
  IO_HI,
  IO_HII,
  IO_HeI,
  IO_HeII,
  IO_HeIII,
  IO_H2I,
  IO_H2II,
  IO_HM,
  IO_DTENTR,
  IO_STRESSDIAG,
  IO_STRESSOFFDIAG,
  IO_STRESSBULK,
  IO_SHEARCOEFF,
  IO_TSTP,
  IO_BFLD,
  IO_DBDT,
  IO_DIVB,
  IO_ABVC,
  IO_AMDC,
  IO_PHI,
  IO_COOLRATE,
  IO_CONDRATE,
  IO_BSMTH,
  IO_DENN,
  IO_MACH,
  IO_DTENERGY,
  IO_PRESHOCK_DENSITY,
  IO_PRESHOCK_ENERGY,
  IO_PRESHOCK_XCR,
  IO_DENSITY_JUMP,
  IO_ENERGY_JUMP,
  IO_CRINJECT,
  IO_TIDALTENSOR,
  IO_DISTORTIONTENSOR,
  IO_LASTENTRY			/* This should be kept - it signals the end of the list */
};



/*
 * Variables for Tree
 * ------------------
 */

extern struct NODE
{
  MyFloat len;			/*!< sidelength of treenode */
  MyFloat center[3];		/*!< geometrical center of node */

#ifdef RADTRANSFER
  MyFloat stellar_mass;         /*!< mass in stars in the node*/
  MyFloat stellar_s[3];         /*!< enter of mass for the stars in the node*/
#endif

#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  MyFloat maxsoft;		/*!< hold the maximum gravitational softening of particles in the 
				   node if the ADAPTIVE_GRAVSOFT_FORGAS option is selected */
#endif
  union
  {
    int suns[8];		/*!< temporary pointers to daughter nodes */
    struct
    {
      MyFloat s[3];		/*!< center of mass of node */
      MyFloat mass;		/*!< mass of node */
      unsigned int bitflags;	/*!< flags certain node properties */
      int sibling;		/*!< this gives the next node in the walk in case the current node can be used */
      int nextnode;		/*!< this gives the next node in case the current node needs to be opened */
      int father;		/*!< this gives the parent node of each node (or -1 if we have the root node) */
    }
    d;
  }
  u;
#ifdef SCALARFIELD
  MyFloat s_dm[3];
  MyFloat mass_dm;
#endif
  int Ti_current;
}
 *Nodes_base,			/*!< points to the actual memory allocted for the nodes */
 *Nodes;			/*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart] 
				   gives the first allocated node */


extern struct extNODE
{
  MyLongDouble dp[3];
#ifdef SCALARFIELD
  MyLongDouble dp_dm[3];
  MyFloat vs_dm[3];
#endif
#ifdef FLTROUNDOFFREDUCTION
  MyFloat s_base[3];
  MyFloat len_base;
#ifdef SCALARFIELD
  MyFloat s_dm_base[3];
#endif
#endif
  MyFloat vs[3];
  MyFloat vmax;
  MyFloat divVmax;
  MyFloat hmax;			/*!< maximum SPH smoothing length in node. Only used for gas particles */
  int Ti_lastkicked;
  int Flag;
}
 *Extnodes, *Extnodes_base;


extern int MaxNodes;		/*!< maximum allowed number of internal nodes */
extern int Numnodestree;	/*!< number of (internal) nodes in each tree */


extern int *Nextnode;		/*!< gives next node in tree walk  (nodes array) */
extern int *Father;		/*!< gives parent node in tree (Prenodes array) */

#ifdef STATICNFW
extern double Rs, R200;
extern double Dc;
extern double RhoCrit, V200;
extern double fac;
#endif


#ifdef CHEMISTRY
/* ----- chemistry part ------- */

#define H_number_fraction 0.76
#define He_number_fraction 0.06

/* ----- Tables ------- */
extern double T[N_T], J0_nu[N_nu], J_nu[N_nu], nu[N_nu];
extern double k1a[N_T], k2a[N_T], k3a[N_T], k4a[N_T], k5a[N_T], k6a[N_T], k7a[N_T], k8a[N_T], k9a[N_T],
  k10a[N_T], k11a[N_T];
extern double k12a[N_T], k13a[N_T], k14a[N_T], k15a[N_T], k16a[N_T], k17a[N_T], k18a[N_T], k19a[N_T],
  k20a[N_T], k21a[N_T];
extern double ciHIa[N_T], ciHeIa[N_T], ciHeIIa[N_T], ciHeISa[N_T], reHIIa[N_T], brema[N_T];
extern double ceHIa[N_T], ceHeIa[N_T], ceHeIIa[N_T], reHeII1a[N_T], reHeII2a[N_T], reHeIIIa[N_T];

/* cross-sections */
#ifdef RADIATION
extern double sigma24[N_nu], sigma25[N_nu], sigma26[N_nu], sigma27[N_nu], sigma28[N_nu], sigma29[N_nu],
  sigma30[N_nu], sigma31[N_nu];
#endif
#endif

#endif


/*structures for eddington tensor*/
#ifdef RADTRANSFER
struct eddingtondata_in
{
  int NodeList[NODELISTLENGTH];
  MyDouble Pos[3];
  MyFloat Hsml;
  MyFloat ET[6];
  MyFloat nHI;
}
 *EddingtonDataIn, *EddingtonDataGet;


struct eddingtondata_out
{
  MyFloat ET[6];
  MyFloat P1, D1;
}
 *EddingtonDataResult, *EddingtonDataOut;
#endif



