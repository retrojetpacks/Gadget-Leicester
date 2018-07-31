#ifdef COSMIC_RAYS
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include "allvars.h"
#include "proto.h"
#include <stdio.h>

#include "cosmic_rays.h"

#include "auxiliary_functions.c"


/*! \file cosmic_rays.c
 *  \brief Computation of Cosmic Ray Physics
 *
 *  This file contains all the basic equations for cosmic ray physics
 *  It is divided into two sections. Internal computation functions and
 *  interfaces to Gadget itself that will plug into different parts of the
 *  code (hydrodynamics, cooling, star formation, ...)
 */


double loc_diss;
double loc_time;
double loc_tau;


double *rBeta1Tab;

double *rEkinTab;
double *rQTab;

#ifdef FIX_QINJ
double *rBetaTab;
double *rAlphaTab;
#endif

#ifdef CR_DIFFUSION
double *rBeta2Tab;
double *rBeta3Tab;
#endif

static gsl_integration_workspace *Integration_Workspace;


const double qthresh = 0.83;	/* Threshold value for hadronic interactions */


#define TabLength 1000
#define a_min 3.0
#define a_max 22.0

#define CR_CHANGE_FAC 0.05





#define QMIN  0.0001
#define QMAX  1000.0
#define TABLEN 1400
#define ALPHAS 40
#define ALPHA_STEP 0.05
#define ALPHA_MIN  2.05
#define ALPHA_MAX (ALPHA_MIN + ALPHAS * ALPHA_STEP)

double Tab_Q[TABLEN + 1];
double Tab_tau_Thermalization[TABLEN + 1];
double Tab_tau_Dissipation[TABLEN + 1];
double Tab_tau_Cooling[TABLEN + 1];
double Tab_beta[TABLEN + 1];

double Tab_alpha[ALPHAS + 1];
double Tab_injection_point[TABLEN + 1][ALPHAS + 1];
double Tab_energy_fraction[TABLEN + 1][ALPHAS + 1];
double Tab_mean_energy[TABLEN + 1][ALPHAS + 1];



double Tab_injection_point_maximum_value[ALPHAS + 1];
double Tab_injection_point_maximum_q[ALPHAS + 1];
int Tab_injection_point_maximum_index[ALPHAS + 1];

double CR_Tab_EnergyFraction(double q, double alpha);
double CR_Tab_MeanEnergy(double q, double alpha);



void CR_Tab_Initialize(void)
{
  int i, j;
  double dlogq, q, alpha;

  dlogq = log(QMAX / QMIN) / (TABLEN);

  for(j = 0; j <= ALPHAS; j++)
    {
      Tab_injection_point_maximum_value[j] = -1.0e60;
      Tab_alpha[j] = ALPHA_MIN + ALPHA_STEP * j;
    }

  for(i = 0; i <= TABLEN; i++)
    {
      Tab_Q[i] = log(QMIN) + i * dlogq;

      q = exp(Tab_Q[i]);

      Tab_tau_Thermalization[i] = log(CR_Particle_GetThermalizationTimescale(q, 1.0));	/* for unit density */
      Tab_tau_Dissipation[i] = log(CR_Particle_GetDissipationTimescale(q, 1.0));	/* for unit density */
      Tab_tau_Cooling[i] = log(CR_Particle_GetCoolingTimescale(q, 1.0));	/* for unit density */
      Tab_beta[i] = Beta((All.CR_Alpha - 2.0) * 0.5, (3.0 - All.CR_Alpha) * 0.5, 1.0 / (1.0 + q * q));

      for(j = 0; j <= ALPHAS; j++)
	{
	  alpha = ALPHA_MIN + ALPHA_STEP * j;

	  Tab_energy_fraction[i][j] = CR_EnergyFractionAfterShiftingQ(QMIN, q, alpha);

	  Tab_mean_energy[i][j] = log(CR_mean_kinetic_energy(q, alpha));

	  Tab_injection_point[i][j] = log(exp(Tab_tau_Cooling[i]) * Tab_energy_fraction[i][j]);

	  if(Tab_injection_point[i][j] > Tab_injection_point_maximum_value[j])
	    {
	      Tab_injection_point_maximum_value[j] = Tab_injection_point[i][j];
	      Tab_injection_point_maximum_q[j] = q;
	      Tab_injection_point_maximum_index[j] = i;
	    }
	}
    }
}

double CR_Tab_Beta(double q)
{
  double x;
  int i;

  if(q < QMIN)
    q = QMIN;

  if(q >= QMAX)
    q = 0.99 * QMAX;

  assert(q >= QMIN && q <= QMAX);

  x = (log(q / QMIN) / log(QMAX / QMIN)) * TABLEN;
  i = (int) x;
  x -= i;

  return Tab_beta[i] * (1 - x) + Tab_beta[i + 1] * x;
}

double CR_Tab_EnergyFraction(double q, double alpha)
{
  int i, j;
  double x;

  if(q < QMIN)
    return 1.0;
  if(q >= QMAX)
    return q = 0.99 * QMAX;

  assert(q >= QMIN && q <= QMAX);

  if(alpha < ALPHA_MIN)
    j = 0;
  else if(alpha >= ALPHA_MAX)
    j = ALPHAS;
  else
    j = (int) ((alpha - ALPHA_MIN + 1.0e-5) / ALPHA_STEP);

  assert(j >= 0 && j <= ALPHAS);

  x = (log(q / QMIN) / log(QMAX / QMIN)) * TABLEN;
  i = (int) x;
  x -= i;

  return Tab_energy_fraction[i][j] * (1 - x) + Tab_energy_fraction[i + 1][j] * x;
}


double CR_Tab_MeanEnergy(double q, double alpha)
{
  double x;
  int i, j;

  if(q < QMIN)
    q = QMIN;

  if(q >= QMAX)
    q = 0.99 * QMAX;

  assert(q >= QMIN && q <= QMAX);

  if(alpha < ALPHA_MIN)
    j = 0;
  else if(alpha >= ALPHA_MAX)
    j = ALPHAS;
  else
    j = (int) ((alpha - ALPHA_MIN + 1.0e-5) / ALPHA_STEP);

  assert(j >= 0 && j <= ALPHAS);

  x = (log(q / QMIN) / log(QMAX / QMIN)) * TABLEN;
  i = (int) x;
  x -= i;
  return exp(Tab_mean_energy[i][j] * (1 - x) + Tab_mean_energy[i + 1][j] * x);
}


double CR_Find_Qinj(double rho, double TimeScale, double alpha, double qstart)
{
  int j, i1, i2, im;
  double x, q;

  if(alpha < ALPHA_MIN)
    j = 0;
  else if(alpha >= ALPHA_MAX)
    j = ALPHAS;
  else
    j = (int) ((alpha - ALPHA_MIN + 1.0e-5) / ALPHA_STEP);

  assert(j >= 0 && j <= ALPHAS);

  TimeScale *= rho;

  if(TimeScale < exp(Tab_injection_point[0][j]))
    return qstart;

  TimeScale = log(TimeScale);


  if(Tab_injection_point_maximum_value[j] < TimeScale)
    return Tab_injection_point_maximum_q[j];

  i1 = 0;
  i2 = Tab_injection_point_maximum_index[j];

  while(i2 - i1 > 1)
    {
      im = i1 + (i2 - i1) / 2;

      if(Tab_injection_point[im][j] < TimeScale)
	i1 = im;
      else
	i2 = im;
    }

  x = (TimeScale - Tab_injection_point[i1][j]) / (Tab_injection_point[i2][j] - Tab_injection_point[i1][j]);

  q = exp(Tab_Q[i1] * (1 - x) + Tab_Q[i2] * x);

  return q;
}


double CR_Tab_GetThermalizationTimescale(double q, double rho)
{
  double tau, x;
  int i;

  if(q < QMIN)
    q = QMIN;

  if(q >= QMAX)
    q = 0.99 * QMAX;

  assert(q >= QMIN && q <= QMAX);

  x = (log(q / QMIN) / log(QMAX / QMIN)) * TABLEN;
  i = (int) x;
  x -= i;
  tau = exp(Tab_tau_Thermalization[i] * (1 - x) + Tab_tau_Thermalization[i + 1] * x);

  return tau / rho;
}

double CR_Tab_GetDissipationTimescale(double q, double rho)
{
  double tau, x;
  int i;

  if(q < QMIN)
    q = QMIN;

  if(q >= QMAX)
    q = 0.99 * QMAX;

  assert(q >= QMIN && q <= QMAX);

  x = (log(q / QMIN) / log(QMAX / QMIN)) * TABLEN;
  i = (int) x;
  x -= i;
  tau = exp(Tab_tau_Dissipation[i] * (1 - x) + Tab_tau_Dissipation[i + 1] * x);

  return tau / rho;
}

double CR_Tab_GetCoolingTimescale(double q, double rho)
{
  double tau, x;
  int i;

  if(q < QMIN)
    q = QMIN;

  if(q >= QMAX)
    q = 0.99 * QMAX;

  assert(q >= QMIN && q <= QMAX);

  x = (log(q / QMIN) / log(QMAX / QMIN)) * TABLEN;
  i = (int) x;
  x -= i;
  tau = exp(Tab_tau_Cooling[i] * (1 - x) + Tab_tau_Cooling[i + 1] * x);

  return tau / rho;
}



double Beta_Integrand(double x, void *Params)
{
  double *dParams = Params;

  return (pow(x, dParams[0] - 1.0) * pow(1.0 - x, dParams[1] - 1.0));
}



double Beta(double A, double B, double x)
{
  double rResult, rError;
  gsl_function Integrand;
  double Integration_Parameters[2];

#ifdef FIX_QINJ
  int iLeft, iRight;
  double rAlpha, rSlope;
#endif


  if(B > 0.0)
    /* The gsl_sf_beta( A, B ) is not defined for B <= 0.0
     * the coefficient A should never be a problem in our
     * case, it would require alpha < 2 to cause issues. */
    {

      rResult = gsl_sf_beta_inc(A, B, x) * gsl_sf_beta(A, B);

      return rResult;
    }

  /* Damn. We cannot use the nice optimized functions. Do some
   * nasty explicit numerical integration integration instead */

  /* Limit alpha_inj < 22 => M > 1.1: */
  if(A > 10.0)
    {
      A = 10.0;
      B = -9.5;
    }

#ifdef FIX_QINJ

  /* Use 1D interpolation table to speed up the routine for alpha > 3: */
#ifdef CR_DIFFUSION
  assert((All.CR_Diffusion_Gamma != 1.0) || (All.CR_Diffusion_Gamma != 2.0));
  /* In these cases, A + B = 0.5, and the diffusion Beta functions are 
     indistinguishable to the Beta functions of CR energy density and pressure,
     however the inferred alpha value is the wrong one! */
#endif

  if(fabs(A + B - 0.5) < 1e-6)
    {
      if(fabs(sqrt(1.0 / x - 1.0) - All.Shock_Fix_Qinj) < 1e-6)
	{

	  rAlpha = A - B + 2.5;

	  if(rAlpha != a_max)
	    {
	      iLeft = (int) ((rAlpha - a_min) / (a_max - a_min) * (double) (TabLength - 1) + 1.0);
	    }
	  else
	    {
	      iLeft = TabLength - 1;
	    }

	  iRight = iLeft + 1;

	  /* now linearly interpolate between the left and right cell */

	  rSlope = (rBetaTab[iRight] - rBetaTab[iLeft]) / (rAlphaTab[iRight] - rAlphaTab[iLeft]);

	  rResult = rBetaTab[iLeft] + rSlope * (rAlpha - rAlphaTab[iLeft]);

	  return rResult;

	}
    }

#endif /* FIX_QINJ */


  Integrand.function = &Beta_Integrand;
  Integrand.params = Integration_Parameters;

  Integration_Parameters[0] = A;
  Integration_Parameters[1] = B;

  if(Integration_Workspace == NULL)
    {
    }

  gsl_integration_qag(&Integrand,	/* Integrand function */
		      1.0e-20, x,	/* Range. Avoid 0^0 by starting nonzero */
		      0, 1e-3, 1000,	/* abs error, rel error, max. bins */
		      GSL_INTEG_GAUSS15,	/* Integration method */
		      Integration_Workspace, &rResult, &rError);

  return rResult;

}


double inline CR_physical_C(double C0, double Density,	/* (physical density) */
			    double Alpha)
     /* Obtain the actual value of the normalization parameter
      * C for the current locus and time. */
{
  return (C0 * pow(Density, (Alpha - 1.0) * 0.33333));
}



double inline CR_physical_q(double q0, double Density	/* (physical density) */
  )
     /* Obtain the actual value of the spectral cutoff parameter
      * q for the current locus and time. */
{
  return (q0 * pow(Density, 0.33333));
}




double CR_mean_kinetic_energy(double q, double Alpha)
     /* Return the mean kinetic energy of a CR particle for a
      * spectrum with a spectral index of alpha and a cutoff of q */
     /* Equation (5) */
{
  double rQFac;

  rQFac = 1.0 + q * q;

  return ((0.5 * pow(q, Alpha - 1.0) *
	   Beta((Alpha - 2.0) * 0.5, (3.0 - Alpha) * 0.5, 1.0 / rQFac) + (sqrt(rQFac) - 1.0)) * mpc2);
}


double inline CR_kinetic_energy(double P)
     /* Returns the kinetic energy of a particle with impulse P */
{
  return ((sqrt(1.0 + (P * P)) - 1.0) * mpc2);
}



double CR_q_from_mean_kinetic_energy(double T)
     /* Inversion of CR_mean_kinetic_energy()
      * this function looks up the q that corresponds
      * to a given mean kinetic energy at the global
      * spectral index for cosmic rays. */
     /* q(T), see section 5.1 */
{
  int iLeft, iRight;
  int iMed;

  double rFraction;
  double rX;

  double rTau;

  /*
     double rB;
     double rQa;
   */

  rTau = T / mpc2;

  /* Can we solve this with Equation (40) ? *//* << 1 */
  /*
     if(rTau < CR_CHANGE_FAC)   
     {
     rB = gsl_sf_beta((All.CR_Alpha - 2.0) * 0.5, (3.0 - All.CR_Alpha) * 0.5);

     rQa = pow(2.0 * rTau / rB, 1.0 / (All.CR_Alpha - 1.0));

     return (rQa + (pow(rQa, 4.0 - All.CR_Alpha) / ((3.0 - All.CR_Alpha) * rB)));
     }
   */

  /* Can we solve this with Equation (41) ? *//* >> 1 */
  /*
     if(rTau * CR_CHANGE_FAC > 1.0)     
     {
     return ((All.CR_Alpha - 2.0) / (All.CR_Alpha - 1.0) * (rTau + 1.0));
     }
   */

  /* the approximate formulae can have substantially larger error (up to 5-10%) than the table.
     For this reason, I switched them off */


  /* No asymptotic solution possible. Hence, do a logarithmic inverse-lookup
   *interpolation in the tabulated mean kinetic energy field. For this, we 
   * localize the cell of the current mean kinetic energy, then compute the
   * location in this cell.
   * Finally, we will convert from 1/(1+q^2) to q
   */

  iLeft = 1;
  iRight = TabLength;

  /* Do a binary cellfind */
  do
    {
      iMed = (iLeft + iRight) / 2;

      if(T < rEkinTab[iMed])
	{
	  iRight = iMed;
	}
      else
	{
	  iLeft = iMed;
	}
    }
  while(iRight > (iLeft + 1));


  /* now logarithmically interpolate between the left and right cell */
  rFraction = (log(T / rEkinTab[iLeft]) / log(rEkinTab[iRight] / rEkinTab[iLeft]));

  rX = rQTab[iLeft] * pow(rQTab[iRight] / rQTab[iLeft], rFraction);

  return rX;
}


int CR_initialize_beta_tabs(double Alpha)
{
  int i;
  double x, q;

#ifdef FIX_QINJ
  double a_tab, A, B;
  double rResult, rError;
  double Integration_Parameters[2];
  gsl_function Integrand;
#endif

  rEkinTab = calloc(TabLength + 1, sizeof(double));
  rQTab = calloc(TabLength + 1, sizeof(double));

  if((rEkinTab == NULL) || (rQTab == NULL))
    {
      printf("ERROR: Could not allocate memory for CR energy table.\n");
      endrun(1337);
    }


  q = 1.0e-3;
  x = pow(1.0e6, 1.0 / ((double) TabLength));

  /* Compute the mean kinetic energy table */
  for(i = 1; i <= TabLength; i++)
    /* NOTE : i = 0 gives us quite some trouble. q->infinity */
    {
      rEkinTab[i] = CR_mean_kinetic_energy(q, Alpha);
      rQTab[i] = q;

      q *= x;
    }

  Integration_Workspace = gsl_integration_workspace_alloc(1000);
  if(Integration_Workspace == NULL)
    {
      printf("ERROR: Could not allocate memory for integration workspace.\n");
      endrun(1338);
    }

#ifdef FIX_QINJ
  rBetaTab = calloc(TabLength + 1, sizeof(double));
  rAlphaTab = calloc(TabLength + 1, sizeof(double));

  if((rBetaTab == NULL) || (rAlphaTab == NULL))
    {
      printf("ERROR: Could not allocate memory for CR Beta table.\n");
      endrun(1337);
    }

  /* initialize tabular for Beta_x(A, B), alpha = (a_min, a_max): */
  Integrand.function = &Beta_Integrand;
  Integrand.params = Integration_Parameters;

  x = 1.0 / (1.0 + pow(All.Shock_Fix_Qinj, 2));

  for(i = 1; i <= TabLength; i++)
    {
      a_tab = a_min + (a_max - a_min) * ((double) (i - 1)) / ((double) (TabLength - 1));

      A = (a_tab - 2.0) / 2.0;
      B = (3.0 - a_tab) / 2.0;

      Integration_Parameters[0] = A;
      Integration_Parameters[1] = B;

      gsl_integration_qag(&Integrand,	/* Integrand function */
			  1.0e-20, x,	/* Range. Avoid 0^0 by starting nonzero */
			  0, 1e-3, 1000,	/* abs error, rel error, max. bins */
			  GSL_INTEG_GAUSS15,	/* Integration method */
			  Integration_Workspace, &rResult, &rError);

      rBetaTab[i] = rResult;
      rAlphaTab[i] = a_tab;
    }
#endif

  /* Finished. Return true */
  return 1;
}







void CR_free_beta_tabs(void)
{
  /* if ( rBeta1Tab != NULL ) myfree( rBeta1Tab ); */
  if(rEkinTab != NULL)
    myfree(rEkinTab);
  if(rQTab != NULL)
    myfree(rQTab);
  if(Integration_Workspace != NULL)
    myfree(Integration_Workspace);
#ifdef FIX_QINJ
  if(rBetaTab != NULL)
    myfree(rBetaTab);
  if(rAlphaTab != NULL)
    myfree(rAlphaTab);
#endif
}



/* ============================================================ */
/* ============ Interface functions to GADGET ================= */
/* ============================================================ */


double CR_MeanKineticEnergy(double q, double Alpha)
     /* Equation (5) */
{
  double rQFac;

  rQFac = 1.0 + q * q;

  return ((0.5 * pow(q, Alpha - 1.0) *
	   Beta((Alpha - 2.0) * 0.5, (3.0 - Alpha) * 0.5, 1.0 / rQFac) + (sqrt(rQFac) - 1.0)) * mpc2);
}


void CR_Particle_UpdateNonAdiabatics(SphParticle * Particle, double rCphys, double rQphys)
{
#ifndef CR_NO_CHANGE
  double rTq;
  double rQmeanKin;
  double rMeanKineticEnergy;
  double rQold, rCold, rho;
  double rEOld, rNOld;

  rEOld = Particle->CR_E0;
  rNOld = Particle->CR_n0;

  rQold = Particle->CR_q0;
  rCold = Particle->CR_C0;

  rho = Physical_Density(Particle);

  /* If the changes are small in comparison to the actual values of
   * E0 and n0, we can apply eq. (10) and (11) to compute their effects
   */
  if((fabs(Particle->CR_DeltaE) < CR_CHANGE_FAC * Particle->CR_E0) &&
     (fabs(Particle->CR_DeltaN) < CR_CHANGE_FAC * Particle->CR_n0))
    {
      /* Compute T_p(q) */
      rTq = CR_kinetic_energy(rQphys);

      /* Compute T_CR */
      rMeanKineticEnergy = Particle->CR_E0 * m_p / Particle->CR_n0;

      /* Equation (10) */
      Particle->CR_C0 *=
	((m_p * Particle->CR_DeltaE - rTq * Particle->CR_DeltaN) /
	 (m_p * Particle->CR_E0 - rTq * Particle->CR_n0)) + 1.0;


      /* Equation (11) */
      Particle->CR_q0 *=
	((m_p * Particle->CR_DeltaE -
	  rMeanKineticEnergy * Particle->CR_DeltaN) /
	 ((m_p * Particle->CR_E0 - rTq * Particle->CR_n0) * (All.CR_Alpha - 1.0))) + 1.0;

      /* Add the injected values */
      Particle->CR_E0 += Particle->CR_DeltaE;
      Particle->CR_n0 += Particle->CR_DeltaN;

      assert(Particle->CR_C0 > 0);
      assert(Particle->CR_q0 > 0);
      assert(Particle->CR_E0 > 0);
      assert(Particle->CR_n0 > 0);
    }
  else
    {
      /* We need to add energy and number density and
       * recompute C0 and q0 from it
       */

      /* If all cosmic ray particles or all cosmic ray energy is used
       * up, cancel all cosmic ray effects of particle
       */
      if((Particle->CR_E0 + Particle->CR_DeltaE) <= 0.0 || (Particle->CR_n0 + Particle->CR_DeltaN) <= 0.0)
	{
	  /*
	     printf("Warning! Strong Drop in CR energy/number. ID\n");
	     printf("EOld=%g DeltaE=%g Enew=%g\n",
	     Particle->CR_E0, Particle->CR_DeltaE, Particle->CR_E0 + Particle->CR_DeltaE);
	     printf("NOld=%g DeltaN=%g Nnew=%g\n",
	     Particle->CR_n0, Particle->CR_DeltaN, Particle->CR_n0 + Particle->CR_DeltaN);

	     printf("Dissip=%g %g %g\n", loc_diss, loc_time, loc_tau);
	     fflush(stdout);
	     endrun(777);
	   */

	  Particle->CR_E0 = 0.0;
	  Particle->CR_n0 = 0.0;

	  Particle->CR_q0 = 1.0e10;
	  Particle->CR_C0 = 0.0;

	  return;
	}

      Particle->CR_E0 += Particle->CR_DeltaE;
      Particle->CR_n0 += Particle->CR_DeltaN;


      rMeanKineticEnergy = Particle->CR_E0 * m_p / Particle->CR_n0;

      /* Compute a corresponding physical q for the given
       * mean kinetic energy
       *
       * For the density, we use CR_Rho0 instead of Density, because
       * CR_Rho0 contains the physical value
       */
      rQmeanKin = CR_q_from_mean_kinetic_energy(rMeanKineticEnergy);

      Particle->CR_q0 = rQmeanKin * pow(rho, -(1.0 / 3.0));

      /* Compute new C0 through the definition of n */
      Particle->CR_C0 = Particle->CR_n0 * (All.CR_Alpha - 1.0) * pow(Particle->CR_q0, All.CR_Alpha - 1.0);
    }

#endif /* CR_NO_CHANGE */

/* Zero out the "unconsidered changes" buffer */
  Particle->CR_DeltaE = 0.0;
  Particle->CR_DeltaN = 0.0;
}


void CR_Particle_Update(SphParticle * Particle)
{
  double rCphys;
  double rQphys;
  double rho;
  double rBeta;

  /* A value of q = 0.0 gives us serious trouble, because it messes up the
   * q^(1-alpha) term. So, set it to an arbitrary small value
   */
  if(Particle->CR_q0 == 0.0)
    {
      Particle->CR_q0 = 1.0e10;
    }

  rho = Physical_Density(Particle);

  /* If there have been non-adiabatic effects during the last timestep,
   * evaluate them and their effect on q0 and c0.
   */
  if((Particle->CR_DeltaE != 0.0) || (Particle->CR_DeltaN != 0.0))
    {
      /* Compute physical C */
      rCphys = CR_physical_C(Particle->CR_C0, rho, All.CR_Alpha);

      /* Compute physical q */
      rQphys = CR_physical_q(Particle->CR_q0, rho);

      /* Compute and cache Beta function */
      rBeta = CR_Tab_Beta(rQphys);

      /* Evaluate Non-Adiabatic process effect.
       * Cphys and qphys are passed so they need not be computed
       * again
       */
      CR_Particle_UpdateNonAdiabatics(Particle, rCphys, rQphys);
    }


  /* Compute physical C */
  rCphys = CR_physical_C(Particle->CR_C0, rho, All.CR_Alpha);

  /* Compute physical q */
  rQphys = CR_physical_q(Particle->CR_q0, rho);


  /* Compute Beta function */
  rBeta = CR_Tab_Beta(rQphys);

  /* Compute the specific energy  */
  Particle->CR_E0 =
    (rCphys * c2) / (All.CR_Alpha - 1.0) *
    (0.5 * rBeta + pow(rQphys, 1.0 - All.CR_Alpha) * (sqrt(1.0 + square(rQphys)) - 1.0));


  /* Compute the baryon fraction of CR particles at rho0. Note: This is an adiabatic invariant. */
  Particle->CR_n0 = (Particle->CR_C0 * pow(Particle->CR_q0, 1.0 - All.CR_Alpha) / (All.CR_Alpha - 1));


  /* Compute pressure at rho0 */
  /*
     Particle->CR_P0 = (rCphys * Particle->CR_Rho0 * rBeta * c2) / 6.0;
   */

  /* Compute adiabatic index at rho0 */
#ifdef MACHNUM
#if  ( CR_SHOCK == 2 )
  if(rBeta != 0.0)
    {
      Particle->CR_Gamma0 = ((All.CR_Alpha + 2.0) / 3.0 -
			     ((2.0 * pow(rQphys, 3.0 - All.CR_Alpha)) /
			      (3.0 * sqrt(1 + square(rQphys)) * rBeta)));
    }
  else
    Particle->CR_Gamma0 = 4. / 3.;	/* value does not matter, since if Beta = 0 -> P_CR = 0 and function 
					 * GetMachNumberCR() uses only 1D Mach finder since XCR < 0.01!      */
#endif
#endif

}


double CR_Physical_Pressure(SphParticle * Particle)
{
  double qphys, Cphys, rho;

  if(Particle->CR_C0 > 0)
    {
      rho = Physical_Density(Particle);

      qphys = Particle->CR_q0 * pow(rho, 0.333333);

      Cphys = Particle->CR_C0 * pow(rho, (All.CR_Alpha - 1.0) * 0.33333);

      return Cphys * rho * CR_Tab_Beta(qphys) * c2 / 6.0;
    }
  else
    return 0.0;
}


inline double CR_Comoving_Pressure(SphParticle * Particle)
{
  if(All.ComovingIntegrationOn)
    {
      return CR_Physical_Pressure(Particle) * pow(All.Time, 3.0 * GAMMA);
    }
  else
    {
      return CR_Physical_Pressure(Particle);
    }
}



inline double CR_Particle_Pressure(SphParticle * Particle)
     /* For compatiblity reasons with old version of code */
{
  return CR_Comoving_Pressure(Particle);
}


double CR_Physical_SpecificEnergy(SphParticle * Particle)
{
  double rCphys, rQphys, rBeta, rho;

  rho = Physical_Density(Particle);

  rCphys = CR_physical_C(Particle->CR_C0, rho, All.CR_Alpha);

  rQphys = CR_physical_q(Particle->CR_q0, rho);

  rBeta = CR_Tab_Beta(rQphys);

  /* Compute the specific energy at rho0 */
  return (rCphys * c2) / (All.CR_Alpha - 1.0) *
    (0.5 * rBeta + pow(rQphys, 1.0 - All.CR_Alpha) * (sqrt(1.0 + square(rQphys)) - 1.0));
}


inline double CR_Particle_SpecificEnergy(SphParticle * Particle)
     /* For compatibility reasons with old version of code */
{
  return CR_Physical_SpecificEnergy(Particle);
}



double inline CR_Particle_BaryonFraction(SphParticle * Particle)
     /* Output: Number of CR particles per baryon (scalefree) */
{
  return Particle->CR_n0;
}



double inline CR_Particle_SpecificNumber(SphParticle * Particle)
     /* Output: Number of CR particles per unit mass (physical units) */
{
  return Particle->CR_n0 / m_p;
}



double inline CR_Particle_MeanKineticEnergy(SphParticle * Particle)
     /* Output: Mean energy per CR particle (physical units) */
{
  return CR_Particle_SpecificEnergy(Particle) / Particle->CR_n0 * m_p;
}



void inline CR_Particle_Inject(SphParticle * Particle,
			       double SpecificEnergyChange, double BaryonFractionChange)
     /* Input: Change in CR energy per unit mass,
      *        Change in CR particles per baryon
      *
      * Output: None
      */
{
  Particle->CR_DeltaE += SpecificEnergyChange;
  Particle->CR_DeltaN += BaryonFractionChange;

#if ( defined( CR_UPDATE_PARANOIA ) && ( CR_UPDATE_PARANOIA >= 2 ) )
  /* If we're /very/ paranoid, we will just recompute the effects of
   * the CR injection instantly.
   * Note: This might result in a signifficant slowdown, as this function
   *       might be called several times per timestep
   */
  CR_Particle_Update(Particle);
#endif
}



void CR_Particle_GetThermalizationRate(SphParticle * Particle, double q0,
				       double *SpecEnergyChangeRate, double *BaryonFractionChangeRate)
     /* Basically, thermalization is supposed to be implemented like this:
      * get the thermalization rate with this function, then with the time-
      * step compute the energy and number change, and call the inject function
      * with that
      */
{
  double rMeanBetaP;
  double rPlasmaFrequency, rNe_plasma;
  double rQ, rX;
  double rC;
  double rNe;

  double rFactor1;		/* 4 pi e^4 / m_e c */
  double rFactor2;		/* 2 m_p c^2 / hbar */
  double rFactor3;		/* - 1/alpha * q^(2-alpha) * (1+q^2)^(-1/2) + 0.5 B  */
  double rBetaFunction;

  double rDensity;

  if(Particle->CR_C0 == 0.0)
    {
      *SpecEnergyChangeRate = 0.0;
      *BaryonFractionChangeRate = 0.0;

      return;
    }

  rDensity = Physical_Density(Particle);

  rQ = q0 * pow(rDensity, 0.3333);

  rC = CR_physical_C(Particle->CR_C0, rDensity, All.CR_Alpha);


  rX = 1.0 / (1.0 + square(rQ));


  rNe = Physical_ElectronDensity(Particle);


  rMeanBetaP = rQ * rQ / sqrt(1.0 + rQ * rQ);

  /* for the Plasma frequency, we put here a fiducial constant electron density, 
     corresponding to typical densities. Since this is only a logarithmic factor,
     this is unimportant, but makes the cooling timescale independent of density,
     instead of introducing an (extremely weak) dependence */

  rNe_plasma = 1.0 * pow(All.UnitLength_in_cm, 3);	/* corresponds to n = 1.0 cm^-3 */

  rPlasmaFrequency = sqrt(4.0 * M_PI * e2 * rNe_plasma / m_e);

  rFactor1 = 4.0 * M_PI * square(e2) / (m_e * LightSpeed);

  rBetaFunction = Beta(0.5 * (All.CR_Alpha - 1.0), 0.5 * (4.0 - All.CR_Alpha), rX);


  rFactor2 = 2.0 * mec2 / HBAR;

  rFactor3 = (-pow(rQ, 2.0 - All.CR_Alpha) / sqrt(1.0 + rQ * rQ) +
	      0.5 * rBetaFunction) / (2.0 - All.CR_Alpha);


  /* Equation (14) */
  *SpecEnergyChangeRate = rFactor1 * rC / m_p * rNe;

  *SpecEnergyChangeRate *=
    (log(rFactor2 * rMeanBetaP / rPlasmaFrequency) / All.CR_Alpha *
     (pow(rQ, -All.CR_Alpha) * sqrt(1.0 + square(rQ)) + rFactor3) - 0.5 * rFactor3);


  *BaryonFractionChangeRate = *SpecEnergyChangeRate * m_p / CR_kinetic_energy(rQ);

}



double CR_Particle_GetDissipationTimescale(double rQ, double rDensity)
{
  double rEnergy, rRate, elec_frac;
  double egy;
  double rNe;

#ifdef COOLING
  elec_frac = (1 + HYDROGEN_MASSFRAC) / (2 * HYDROGEN_MASSFRAC);
#else
  elec_frac = 1.0;
#endif

  rNe = rDensity * HYDROGEN_MASSFRAC * elec_frac / m_p;

  egy = CR_mean_kinetic_energy(rQ, All.CR_Alpha) * pow(rQ, 1 - All.CR_Alpha) / (All.CR_Alpha - 1) / m_p;

  if(rQ < qthresh)		/* rQ > q_th */
    {
      rQ = qthresh;
    }

  rEnergy = CR_mean_kinetic_energy(rQ, All.CR_Alpha);

  rRate = LightSpeed * rNe / (1.0 + HYDROGEN_MASSFRAC) * pow(rQ, 1 - All.CR_Alpha) * 0.019132 * (cm * cm / g) *	/* sigma_pp / m_p */
    rEnergy / (All.CR_Alpha - 1);

  return egy / rRate;
}


void CR_Particle_GetDissipationRate(SphParticle * Particle,
				    double *SpecEnergyChangeRate, double *BaryonFractionChangeRate)
     /* quantities output are physical, not comoving */
{
  double rQ;
  double rC;
  double rEnergy;
  double rDensity;
  double rNe;

  if(Particle->CR_C0 == 0.0)
    {
      *SpecEnergyChangeRate = 0.0;
      *BaryonFractionChangeRate = 0.0;

      return;
    }

  rDensity = Physical_Density(Particle);

  rNe = Physical_ElectronDensity(Particle);

  rQ = CR_physical_q(Particle->CR_q0, rDensity);

  rC = CR_physical_C(Particle->CR_C0, rDensity, All.CR_Alpha);

  if(rQ < qthresh)		/* rQ > q_th */
    rQ = qthresh;

  rEnergy = CR_mean_kinetic_energy(rQ, All.CR_Alpha);

  *SpecEnergyChangeRate = LightSpeed * rNe / (1.0 + HYDROGEN_MASSFRAC) * rC * pow(rQ, 1 - All.CR_Alpha) * 0.019132 * (cm * cm / g) *	/* sigma_pp / m_p */
    rEnergy / (All.CR_Alpha - 1);

  *BaryonFractionChangeRate = 0.0;
}




/* Get the Coulomb timescale for momentum q and density rho (both physical values)
 */
double CR_Particle_GetCoulombTimescale(double rQ, double rDensity)
{
  double rBeta, rPlasmaFrequency, rNe, prefac;
  double rLogFactor, elec_frac;

#ifdef COOLING
  elec_frac = (1 + HYDROGEN_MASSFRAC) / (2 * HYDROGEN_MASSFRAC);
#else
  elec_frac = 1.0;
#endif

  rNe = rDensity * HYDROGEN_MASSFRAC * elec_frac / m_p;

  rBeta = rQ / sqrt(1.0 + square(rQ));

  rPlasmaFrequency = sqrt(4.0 * M_PI * e2 * rNe / m_e);

  rLogFactor = log(2.0 * mec2 * rBeta * rQ / (HBAR * rPlasmaFrequency));

  prefac = (m_e * LightSpeed * LightSpeed * LightSpeed * m_p * rBeta) / (4 * M_PI * rNe * square(e2));

  return prefac * (sqrt(1 + square(rQ)) - 1) / (rLogFactor - 0.5 * rBeta * rBeta);
}


double CR_Particle_GetCoolingTimescale(double rQ, double rDensity)
{
  double t1, t2;

  t1 = CR_Particle_GetThermalizationTimescale(rQ, rDensity);
  t2 = CR_Particle_GetDissipationTimescale(rQ, rDensity);

  return 1.0 / (1.0 / t1 + 1.0 / t2);
}



/* Get the averaged Coulomb thermalization timescale for a spectrum with given cutoff q and density rho (both physical values)
 */
double CR_Particle_GetThermalizationTimescale(double rQ, double rDensity)
{
  double rMeanBetaP;
  double rPlasmaFrequency, rNe_plasma;
  double rX, elec_frac;
  double rNe;
  double prefac;
  double rLogFactor;
  double rBetaFunction0, rBetaFunction1, rBetaFunction2;

#ifdef COOLING
  elec_frac = (1 + HYDROGEN_MASSFRAC) / (2 * HYDROGEN_MASSFRAC);
#else
  elec_frac = 1.0;
#endif

  rNe = rDensity * HYDROGEN_MASSFRAC * elec_frac / m_p;

  /* for the Plasma frequency, we put here a fiducial constant electron density, 
     corresponding to typical densities. Since this is only a logarithmic factor,
     this is unimportant, but makes the cooling timescale independent of density,
     instead of introducing an (extremely weak) dependence */

  rNe_plasma = 1.0 * pow(All.UnitLength_in_cm, 3);	/* corresponds to n = 1.0 cm^-3 */

  rX = 1.0 / (1.0 + square(rQ));

  rBetaFunction0 = Beta(0.5 * (All.CR_Alpha - 2.0), 0.5 * (3 - All.CR_Alpha), rX);
  rBetaFunction1 = Beta(0.5 * (All.CR_Alpha - 1.0), -0.5 * All.CR_Alpha, rX);
  rBetaFunction2 = Beta(0.5 * (All.CR_Alpha - 1.0), -0.5 * (All.CR_Alpha - 2), rX);

  rMeanBetaP = rQ * rQ / sqrt(1.0 + rQ * rQ);

  rPlasmaFrequency = sqrt(4.0 * M_PI * e2 * rNe_plasma / m_e);

  rLogFactor = log(2.0 * mec2 * rMeanBetaP / (HBAR * rPlasmaFrequency));

  prefac =
    (m_e * LightSpeed * LightSpeed * LightSpeed * m_p) / (2 * M_PI * (All.CR_Alpha - 1.0) * rNe * square(e2));

  return prefac * (0.5 * rBetaFunction0 +
		   pow(rQ,
		       1 - All.CR_Alpha) * (sqrt(1 + square(rQ)) - 1)) / (rLogFactor * rBetaFunction1 -
									  0.5 * rBetaFunction2);
}


/* Compute the cutoff q at which the cooling timescale is equal to 
 * a time step.
 */
double CR_Particle_CoulombTimeScaleQ(SphParticle * Particle, double TimeScale)
{
#ifdef CR_THERMALIZATION

  double rQinj;
  double rT, rdT;
  double rNe, rPlasmaFrequency;
  double rBeta;

  if(TimeScale == 0.0)
    return 0.0;

  rQinj = 1.0;

  rT = TimeScale;
  rdT = TimeScale;


  /* Precompute electron density (we need this all the time) */
  rNe = Physical_ElectronDensity(Particle);

  /* Precompute plasma frequency */
  rPlasmaFrequency = sqrt(4.0 * M_PI * e2 * rNe / m_e);

  do
    {
      /* adjust the guess for Qinj according to the
       * difference between actual timescale and desired one
       */
      rQinj /= pow(rT / rdT, 0.3);

      rT = CR_kinetic_energy(rQinj);

      rBeta = rQinj / sqrt(1.0 + square(rQinj));

      rdT =
	4.0 * M_PI * square(e2) * rNe /
	(m_e * rBeta * LightSpeed) *
	(log(2.0 * mec2 * rBeta * rQinj / (HBAR * rPlasmaFrequency)) - square(rBeta) / 2.0);

      rdT *= TimeScale;

    }
  while((rT < 0.9 * rdT) || (rT * 0.9 > rdT));

  return rQinj;

#else /* ! CR_THERMALIZATION */

#ifdef FIX_QINJ
  /* inject only CRps with threshold cutoff for the hadronic CRp-p reaction: qinj = 0.828; */
  return All.Shock_Fix_Qinj;
#else
  return 0.0;
#endif

#endif
}

/* Compute the fraction of CR number density still contained in CR 
 *if q is shifted from Q1 to Q2 */
double CR_NumberFractionAfterShiftingQ(double Q1, double Q2, double SpectralIndex)
{
  return (pow(Q2 / Q1, (1.0 - SpectralIndex)));
}

/* Compute the fraction of CR energy still contained in CR 
 *if q is shifted from Q1 to Q2 */
double CR_EnergyFractionAfterShiftingQ(double Q1, double Q2, double SpectralIndex)
{
  return (CR_mean_kinetic_energy(Q2, SpectralIndex) /
	  CR_mean_kinetic_energy(Q1, SpectralIndex) * pow(Q2 / Q1, (1.0 - SpectralIndex)));
}




#ifdef SFR
double inline CR_Particle_SupernovaFeedback(SphParticle * Particle,
					    double SpecificEnergyInput, double TimeScale)
{
  double rho, rQcooling, rQinj, rInjectedEnergy, rTcooling, rBaryonFractionInput;

  rho = Physical_Density(Particle);

  /* pick a  supernova injection scale at T=10^7K (which is a fiducial thermal scale) */
  rQinj = sqrt(k_B * 1.0e7 / mpc2);

  rQcooling = CR_Find_Qinj(rho, TimeScale, All.CR_SNAlpha, rQinj);

  rInjectedEnergy = SpecificEnergyInput * CR_Tab_EnergyFraction(rQcooling, All.CR_SNAlpha);

  /* Compute mean kinetic energy per particle for cooling cutoff */
  rTcooling = CR_Tab_MeanEnergy(rQcooling, All.CR_SNAlpha);

  /* Coulomb cutoff takes place */
  rBaryonFractionInput = rInjectedEnergy * m_p / rTcooling;

  CR_Particle_Inject(Particle, rInjectedEnergy, rBaryonFractionInput);

  return (SpecificEnergyInput - rInjectedEnergy);
}
#endif


#if defined( CR_SHOCK )

#if ( CR_SHOCK == 1 )
double CR_Particle_ShockInject(SphParticle * Particle, double SpecificEnergyInput, double TimeStep)
{
  /* Return the amount of energy that remains in the CR */
  double rInjectedEnergy, rBaryonFractionInput, rho, tinj, rQcooling, rTcooling;

  rho = Physical_Density(Particle);

  tinj = Particle->CR_E0 / (All.CR_ShockEfficiency * SpecificEnergyInput / TimeStep);

  /* determine injection cut-off */

  rQcooling = CR_Find_Qinj(rho, tinj, All.CR_ShockAlpha, QMIN);

  rInjectedEnergy =
    All.CR_ShockEfficiency * SpecificEnergyInput * CR_Tab_EnergyFraction(rQcooling, All.CR_ShockAlpha);

  /* Compute mean kinetic energy per particle for cooling cutoff */
  rTcooling = CR_Tab_MeanEnergy(rQcooling, All.CR_ShockAlpha);

  /* Coulomb cutoff takes place */
  rBaryonFractionInput = rInjectedEnergy * m_p / rTcooling;

  assert(rInjectedEnergy > 0);
  assert(rBaryonFractionInput > 0);

  CR_Particle_Inject(Particle, rInjectedEnergy, rBaryonFractionInput);

  return rInjectedEnergy;
}

#else

double CR_Particle_ShockInject(SphParticle * Particle, double SpecificEnergyInput, double TimeStep)
{
  double q_inj, alpha_inj, eta_lin, zeta_lin, zeta_inj, egy_lin;
  double rInjectedEnergy, rQcooling, rTcooling, rBaryonFractionInput, rho, tinj;
  double M, temp_jump, dens_jump, egy, egy_diss;

#if ( CR_SHOCK == 3 )

  /* this is presently designed to work with CR_SHOCK == 3 */

  M = Particle->Shock_MachNumber;

  if(M < 1.25)			/* ignore very weak shocks */
    return 0;

  temp_jump =
    (2 * GAMMA * M * M - GAMMA_MINUS1) * (GAMMA_MINUS1 * M * M + 2) / ((GAMMA + 1) * (GAMMA + 1) * M * M);

  dens_jump = (GAMMA + 1) * M * M / (GAMMA_MINUS1 * M * M + 2);

  /* get current thermal energy of particle - will take this as proxy for pre-shock temperature */

  egy = Particle->Entropy / GAMMA_MINUS1 * pow(Physical_Density(Particle), GAMMA_MINUS1);

  alpha_inj = (4 - M * M * (1 - 3 * GAMMA)) / (2 * (M * M - 1));

#else /* this is  CR_SHOCK == 2  */

  M = Particle->Shock_MachNumber;

  if(M < 1.25)			/* ignore very weak shocks */
    return 0;

  temp_jump = Particle->Shock_EnergyJump;

  dens_jump = Particle->Shock_DensityJump;

  /* Shock_DensityJump = 1.0 leads to an infinite value of alpha: */
  if(Particle->Shock_DensityJump - 1.0 < 1.0e-6)
    Particle->Shock_DensityJump = 1.0 + 1.0e-6;

  egy = Particle->PreShock_PhysicalEnergy;

  alpha_inj = (Particle->Shock_DensityJump + 2.0) / (Particle->Shock_DensityJump - 1.0);

#endif

  /* Too small alpha values (for compression ratio = 4) will screw
   * up things in the beta function. So avoid that numerical effect by
   * limiting Alpha slightly larger than 2. */
  if(alpha_inj < ALPHA_MIN)
    alpha_inj = ALPHA_MIN;

  if(alpha_inj > ALPHA_MAX)
    alpha_inj = ALPHA_MAX;

  /* note:  CR_ShockCutoff ~ 3.5 */
  /* assume that we'll yet have to jump to the postshock temperature */
  q_inj = All.CR_ShockCutoff / LightSpeed * sqrt(2.0 * GAMMA_MINUS1 * temp_jump * egy);

  eta_lin = 4.0 * cube(All.CR_ShockCutoff) * exp(-square(All.CR_ShockCutoff)) / (1.77 /* sqrt(pi) */  *
										 (alpha_inj - 1.0));

  egy_lin = eta_lin * CR_Tab_MeanEnergy(q_inj, alpha_inj) / m_p;

  egy_diss = egy * (temp_jump - pow(dens_jump, GAMMA_MINUS1));

  if((egy_diss / egy) < 1e-10)
    egy_diss = egy * 1e-10;

  zeta_lin = egy_lin / egy_diss;

  zeta_inj = (1.0 - exp(-zeta_lin / All.CR_ShockEfficiency)) * All.CR_ShockEfficiency;

  rho = Physical_Density(Particle);

  tinj = Particle->CR_E0 / (zeta_inj * SpecificEnergyInput / TimeStep);

  /* determine injection cut-off */

  rQcooling = CR_Find_Qinj(rho, tinj, alpha_inj, q_inj);

  rInjectedEnergy = zeta_inj * SpecificEnergyInput * CR_Tab_EnergyFraction(rQcooling, alpha_inj);

  /* Compute mean kinetic energy per particle for cooling cutoff */
  rTcooling = CR_Tab_MeanEnergy(rQcooling, alpha_inj);

  /* Coulomb cutoff takes place */
  rBaryonFractionInput = rInjectedEnergy * m_p / rTcooling;

  CR_Particle_Inject(Particle, rInjectedEnergy, rBaryonFractionInput);

  return rInjectedEnergy;
}
#endif
#endif




/* Do thermalization and dissipation over a time of "Time",
 * return the amount of ReThermalized Energy to be added to
 * the gas internal energy reservoir
 */
double CR_Particle_ThermalizeAndDissipate(SphParticle * Particle, double Time)
{
  double rho;
  double rThermalizedEnergy = 0.0;

#ifdef CR_THERMALIZATION
  double q0, q1, qm, egy0, egy1, fac, tau_therm;
  double rEnergyRate;

  /*
     double rNumberRate;
   */
#endif
#ifdef CR_DISSIPATION
  double q, deltaE, tau_diss;
#endif

  CR_Particle_Update(Particle);

  rho = Physical_Density(Particle);

#ifdef CR_THERMALIZATION	/* Thermalization first */

  q0 = q1 = Particle->CR_q0 * pow(rho, 0.33333);

  egy0 = Particle->CR_E0;

  tau_therm = CR_Tab_GetThermalizationTimescale(q0, rho);


  /* note: it's better to use the implicit solver all the time. The relative change in number density
     can be large for coulomb cooling even if the change in energy is still relatively small */

  /*
     if(Time < CR_CHANGE_FAC * tau_therm)
     {
     rEnergyRate = egy0 / tau_therm;
     rNumberRate = rEnergyRate * m_p / CR_kinetic_energy(q0);

     assert(rEnergyRate * Time < egy0);

     CR_Particle_Inject(Particle, -rEnergyRate * Time, -rNumberRate * Time);
     }
     else
   */
  {
    /* implicit solver for the new spectral cut-off, amplitude CR_Q0 stays constant */

    fac = egy0 / CR_Tab_EnergyFraction(q0, All.CR_Alpha);

    /* first, find a bracket */
    do
      {
	q1 *= 2;
	egy1 = fac * CR_Tab_EnergyFraction(q1, All.CR_Alpha);
	rEnergyRate = egy1 / CR_Tab_GetThermalizationTimescale(q1, rho);

	if(q1 > QMAX)		/* ok, that's thermalize the rest of the energy */
	  {
	    Particle->CR_E0 = 0.0;
	    Particle->CR_n0 = 0.0;

	    Particle->CR_q0 = 1.0e10;
	    Particle->CR_C0 = 0.0;
	    q1 = q0 = Particle->CR_q0 * pow(rho, 0.33333);
	    break;
	  }
      }
    while(rEnergyRate * Time > egy0 - egy1);

    /* now that we bracketed the solution, let's find an accurate solution by bisection */
    while((q1 - q0) > 1.0e-6 * q1)
      {
	qm = (q1 + q0) / 2;
	egy1 = fac * CR_Tab_EnergyFraction(qm, All.CR_Alpha);
	rEnergyRate = egy1 / CR_Tab_GetThermalizationTimescale(qm, rho);

	if(rEnergyRate * Time - egy0 + egy1 > 0)
	  q0 = qm;
	else
	  q1 = qm;
      }

    Particle->CR_q0 = 0.5 * (q1 + q0) * pow(rho, -0.33333);
  }

  CR_Particle_Update(Particle);

  rThermalizedEnergy = egy0 - Particle->CR_E0;
#endif

#ifdef CR_DISSIPATION
  q = Particle->CR_q0 * pow(rho, 0.33333);

  tau_diss = CR_Tab_GetDissipationTimescale(q, rho);

  deltaE = Particle->CR_E0 * (1 - exp(-Time / tau_diss));

  loc_diss = deltaE;
  loc_time = Time;
  loc_tau = tau_diss;

  CR_Particle_Inject(Particle, -deltaE, 0.0);	/* note: no change of number density */

  CR_Particle_Update(Particle);
#endif

  return rThermalizedEnergy;
}


double CR_corresponding_q(double qinj, double alpha)
{
  double egy, qm, q1, q2, egym, egy1, egy2, x;
  int i1, i2, im;

  i1 = 0;
  i2 = TABLEN;

  egy = CR_Tab_MeanEnergy(qinj, alpha);

  while(i2 - i1 > 1)
    {
      im = i1 + (i2 - i1) / 2;
      qm = exp(Tab_Q[im]);

      egym =
	CR_kinetic_energy(qm) * (1 +
				 CR_Tab_GetThermalizationTimescale(qm,
								   1.0) / CR_Tab_GetDissipationTimescale(qm,
													 1.0));
      if(egym < egy)
	i1 = im;
      else
	i2 = im;
    }

  q1 = exp(Tab_Q[i1]);
  q2 = exp(Tab_Q[i2]);

  egy1 =
    CR_kinetic_energy(q1) * (1 +
			     CR_Tab_GetThermalizationTimescale(q1, 1.0) / CR_Tab_GetDissipationTimescale(q1,
													 1.0));
  egy2 =
    CR_kinetic_energy(q2) * (1 +
			     CR_Tab_GetThermalizationTimescale(q2, 1.0) / CR_Tab_GetDissipationTimescale(q2,
													 1.0));

  x = (egy - egy1) / (egy2 - egy1);

  return exp(Tab_Q[i1] * (1 - x) + Tab_Q[i2] * x);
}





#ifdef COSMIC_RAY_TEST
void CR_test_routine(void)
{
  double a3inv = 1, tsfr, factorEVP, egyhot, ne, tcool, y, x, egyeff, dtime, rQinj;
  double cloudmass, p, sm, utherm, tcr, th, tinj, instant_reheat, heatingrate_effective, coolrate,
    heatingrate;
  double trelax, egycurrent, coolrate_effective, qinj, fraction, egy_cr, bet, P_cr;
  int i, j, nsteps = 0;
  double temp, temp2, u_to_temp_fac, before, after, tdyn, q;
  double egy, egy_sn, meanweight;
  FILE *fd;


  u_to_temp_fac = (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1
    * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

  dtime = 0.001;

  i = 0;

  if(ThisTask == 0)
    {
      fd = fopen("CR_eqn_of_state.txt", "w");
      for(SphP[i].d.Density = 1.01 * All.PhysDensThresh; SphP[i].d.Density <= 10000.0 * All.PhysDensThresh;
	  SphP[i].d.Density *= 1.05)
	{
	  tsfr = sqrt(All.PhysDensThresh / (SphP[i].d.Density * a3inv)) * All.MaxSfrTimescale;

	  factorEVP = pow(SphP[i].d.Density * a3inv / All.PhysDensThresh, -0.8) * All.FactorEVP;

	  egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

	  ne = SphP[i].Ne;
	  tcool = GetCoolingTime(egyhot, SphP[i].d.Density * a3inv, &ne);
	  SphP[i].Ne = ne;

	  y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);

	  x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

	  cloudmass = x * P[i].Mass;

	  sm = (1 - All.FactorSN) * dtime / tsfr * cloudmass;	/* amount of stars expect to form */

	  p = sm / P[i].Mass;

	  egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

	  th = egyeff / ((1 - All.FactorSN) / tsfr * cloudmass / P[i].Mass * All.FeedbackEnergy);

	  trelax = tsfr * (1 - x) / x / (All.FactorSN * (1 + factorEVP));

	  tdyn = 1 / sqrt(4 * M_PI * All.G * SphP[i].d.Density);


	  heatingrate = (1 - All.FactorSN) / tsfr * cloudmass / P[i].Mass * All.FeedbackEnergy * All.CR_SNEff;


	  j = (int) ((All.CR_SNAlpha - ALPHA_MIN + 1.0e-5) / ALPHA_STEP);

	  qinj = Tab_injection_point_maximum_q[j];

	  q = CR_corresponding_q(qinj, All.CR_SNAlpha);

	  fraction = CR_Tab_EnergyFraction(qinj, All.CR_SNAlpha);

	  egy_cr = (p * All.FeedbackEnergy * All.CR_SNEff / dtime) * fraction
	    * CR_Tab_GetCoolingTimescale(q, SphP[i].d.Density);

	  bet = CR_Tab_Beta(q);

	  P_cr =
	    bet / 6.0 * (All.CR_Alpha - 1) * SphP[i].d.Density * egy_cr / (0.5 * bet +
									   pow(q,
									       1 - All.CR_Alpha) * (sqrt(1 +
													 square
													 (q))
												    - 1));

	  fprintf(fd, "%g  %g %g %g %g %g    %g   %g %g  %g %g  %g\n",
		  SphP[i].d.Density, tsfr, tcool, th, trelax, tdyn, egyeff, qinj, q, egy_cr, P_cr, fraction);
	}
      fclose(fd);

      SphP[i].d.Density = 1.01 * All.PhysDensThresh;

      fd = fopen("trcool.txt", "w");
      rQinj = sqrt(k_B * 1.0e7 / mpc2);
      for(q = 0.001; q <= 1000.0; q *= 1.025)
	{
	  egy = CR_mean_kinetic_energy(q, All.CR_Alpha) / CR_kinetic_energy(q);
	  egy_sn = CR_MeanKineticEnergy(q, All.CR_SNAlpha) / CR_kinetic_energy(q);

	  fprintf(fd, "%g  %g %g %g  %g %g  %g %g\n", q,
		  CR_Particle_GetThermalizationTimescale(q, SphP[i].d.Density),
		  CR_Particle_GetDissipationTimescale(q, SphP[i].d.Density),
		  CR_Particle_GetCoulombTimescale(q, SphP[i].d.Density), egy, egy_sn,
		  CR_Particle_GetCoolingTimescale(q, SphP[i].d.Density),
		  CR_EnergyFractionAfterShiftingQ(rQinj, q, All.CR_SNAlpha));
	}
      fclose(fd);

      fd = fopen("tcool_thermal.txt", "w");
      for(temp = 1.0e3; temp <= 1.0e9; temp *= 1.02)
	{
	  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */

#ifdef CONSTANT_MEAN_MOLECULAR_WEIGHT
  meanweight = All.MeanWeight;
#endif

	  egyhot = All.UnitMass_in_g / All.UnitEnergy_in_cgs / meanweight *
	    (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * temp;

	  tcool = GetCoolingTime(egyhot, SphP[i].d.Density * a3inv, &ne);

	  meanweight = 4.0 / (3 * HYDROGEN_MASSFRAC + 1 + 4 * HYDROGEN_MASSFRAC * ne);

#ifdef CONSTANT_MEAN_MOLECULAR_WEIGHT
  meanweight = All.MeanWeight;
#endif

	  temp2 = 1.0 / (All.UnitMass_in_g / All.UnitEnergy_in_cgs / meanweight *
			 (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS)) * egyhot;

	  fprintf(fd, "%g  %g\n", temp2, tcool);
	}
      fclose(fd);
    }

  SphP[i].d.Density = 1.01 * All.PhysDensThresh;

  if(ThisTask == 0)
    do
      {
	nsteps++;

	printf("\n");

	tsfr = sqrt(All.PhysDensThresh / (SphP[i].d.Density * a3inv)) * All.MaxSfrTimescale;

	factorEVP = pow(SphP[i].d.Density * a3inv / All.PhysDensThresh, -0.8) * All.FactorEVP;

	egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

	ne = SphP[i].Ne;
	tcool = GetCoolingTime(egyhot, SphP[i].d.Density * a3inv, &ne);
	SphP[i].Ne = ne;

	y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);

	x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

	egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

	trelax = tsfr * (1 - x) / x / (All.FactorSN * (1 + factorEVP));

	egycurrent = SphP[i].Entropy * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;

	temp = u_to_temp_fac * (SphP[i].Entropy) /
	  GAMMA_MINUS1 * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);

	SphP[i].Entropy =
	  (egyeff +
	   (egycurrent -
	    egyeff) * exp(-dtime / trelax)) * GAMMA_MINUS1 / pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);

	cloudmass = x * P[i].Mass;

	tcr =
	  CR_Particle_GetCoolingTimescale(SphP[i].CR_q0 * pow(SphP[i].d.Density, 0.3333), SphP[i].d.Density);

	th = egyeff / ((1 - All.FactorSN) / tsfr * cloudmass / P[i].Mass * All.FeedbackEnergy);

	tdyn = 1 / sqrt(4 * M_PI * All.G * SphP[i].d.Density);

	printf("tsfr=%g  tcool=%g  trelax=%g  tcr=%g  th=%g  tdyn=%g\n", tsfr, tcool, trelax, tcr, th, tdyn);


	sm = (1 - All.FactorSN) * dtime / tsfr * cloudmass;	/* amount of stars expect to form */

	p = sm / P[i].Mass;

	heatingrate = (1 - All.FactorSN) / tsfr * cloudmass / P[i].Mass * All.FeedbackEnergy * All.CR_SNEff;


	tinj = SphP[i].CR_E0 / (p * All.FeedbackEnergy * All.CR_SNEff / dtime);

	CR_Particle_Update(SphP + i);
	printf("[1] SphP[i].CR_q0=%g SphP[i].CR_C0=%g  SphP[i].CR_E0=%g\n", SphP[i].CR_q0, SphP[i].CR_C0,
	       SphP[i].CR_E0);
	before = SphP[i].CR_E0;
	instant_reheat = CR_Particle_SupernovaFeedback(&SphP[i], p * All.FeedbackEnergy * All.CR_SNEff, tinj);
	CR_Particle_Update(SphP + i);
	after = SphP[i].CR_E0;

	heatingrate_effective = (p * All.FeedbackEnergy * All.CR_SNEff - instant_reheat) / dtime;

	printf("SN: inject=%g  eff=%g  act=%g\n", p * All.FeedbackEnergy * All.CR_SNEff,
	       heatingrate_effective * dtime, after - before);

	heatingrate_effective = (after - before) / dtime;

	printf("[2] SphP[i].CR_q0=%g SphP[i].CR_C0=%g  SphP[i].CR_E0=%g\n", SphP[i].CR_q0, SphP[i].CR_C0,
	       SphP[i].CR_E0);

	coolrate = SphP[i].CR_E0 / CR_Tab_GetCoolingTimescale(SphP[i].CR_q0 * pow(SphP[i].d.Density, 0.3333),
							      SphP[i].d.Density);

	before = SphP[i].CR_E0;
	utherm = CR_Particle_ThermalizeAndDissipate(SphP + i, dtime);
	after = SphP[i].CR_E0;

	coolrate_effective = (before - after) / dtime;

	printf("thermalize: utherm=%g  %g\n", utherm, before - after);

	printf("[3] SphP[i].CR_q0=%g SphP[i].CR_C0=%g  SphP[i].CR_E0=%g\n", SphP[i].CR_q0, SphP[i].CR_C0,
	       SphP[i].CR_E0);

	printf
	  ("step=%d dt=%g egycurrent=%g(temp=%g)  CR_e0=%g(temp=%g) sn_raw=%g  cr_heat=%g  cr_cool=%g (%g)   Q0=%g \n",
	   nsteps, dtime,
	   egycurrent, temp,
	   SphP[i].CR_E0, u_to_temp_fac * SphP[i].CR_E0,
	   heatingrate, heatingrate_effective, coolrate_effective, coolrate,
	   SphP[i].CR_q0 * pow(SphP[i].d.Density, 0.3333));

      }
    while(nsteps < 5000);

  endrun(0);
}
#endif



#endif /* COSMIC RAYS */
