#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

/*! \file ene_per_gram.c
 *  \brief this routine calculates the internal energy per gram of a particle at a specific temperature
considering the all possible micro physics.
 */

#if defined(REAL_EOS)

/* "static" means that the structure is only used here and not passed outside
 * of this program */

double ene_per_gram(int idx, double rho, double temp)
{

double shc_chem[4],shc_dissene,shc_thetavib,shc_am,shc_temp0;
double shc_mu,shc_hatom,shc_uh,shc_uh2,shc_uh2tran,shc_uh2vib,shc_uh2diss,shc_uhe,shc_um;
double shc_dummy;
double shc_utotal;


/* For the simplest case */

/* double u_to_temp_fac = */
/*        All.MeanWeight  * PROTONMASS / BOLTZMANN * GAMMA_MINUS1; */
/*        shc_utotal = temp / u_to_temp_fac; */
/*        return shc_utotal; */




shc_chem[1] = 0.7;
shc_chem[2] = 0.22;
shc_chem[3] = 0.02;
shc_dissene = 4.4773 * 1.6021927e-12;
shc_thetavib = 6100.;
shc_am = 16.78;
shc_temp0 = 273.;




/* Find H/H2 fraction */
if (temp <= shc_temp0) {
     shc_dummy = 0.;
} else {
     shc_dummy = 2.11/(rho*shc_chem[1]) * exp(-52490./temp);
}
shc_hatom = 1./2.*(-1.*shc_dummy+sqrt(shc_dummy*(shc_dummy+4.)));

/* Energy of H */
shc_uh = 3./2. * shc_chem[1]*shc_hatom*BOLTZMANN/PROTONMASS * temp;

/* Energy of H2 : Translational */
shc_uh2tran = 3./2. * BOLTZMANN * temp;

/* Energy of H2 : Vibrational */
if (temp <= 50.) {
     shc_uh2vib = 0.;
} else {
     shc_uh2vib = BOLTZMANN*shc_thetavib / (exp(shc_thetavib/temp)-1.);
}

/* Energy of H2 : Total */
shc_uh2 = shc_uh2tran + shc_uh2vib + BOLTZMANN * shc_uh2rot(temp);
shc_uh2 = shc_chem[1] * (1.-shc_hatom)*shc_uh2 / (2.*PROTONMASS);

/* Energy of H2 : Dissociation */
shc_uh2diss = 1./2. * shc_chem[1] * shc_hatom * shc_dissene / PROTONMASS;

/* Energy of He */
shc_uhe = 3./8. * shc_chem[2] * BOLTZMANN / PROTONMASS * temp;

/* Energy of Metal */
shc_um = 3./2. * shc_chem[3] * BOLTZMANN / PROTONMASS * temp / shc_am;

/* Mean Molecular Weight, mu */
shc_mu = 1./(shc_chem[1]/2.*(1.+shc_hatom) + shc_chem[2]/4. + shc_chem[3]/shc_am);
SphP[idx].Mu = shc_mu;

/* Internal energy of a particle per unit gram */
shc_utotal = shc_uh + shc_uh2 + shc_uh2diss + shc_uhe + shc_um;

/*
fprintf (diag_shc,"ENE_PER_GRAM : %i %g %g %g\n",idx,rho,temp,shc_utotal);
fprintf (diag_shc,"2 : %g %g %g %g %g\n",shc_hatom,shc_uh,shc_uh2diss,shc_uhe,shc_um);
fprintf (diag_shc,"3 : %g %g %g %g %g\n",shc_uh,shc_uh2,shc_uh2diss,shc_uhe,shc_um);
fflush (diag_shc);
*/

/* For the simplest case */
/*
double u_to_temp_fac =
     All.MeanWeight  * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
//       All.MeanWeight  * PROTONMASS / BOLTZMANN * GAMMA_MINUS1;
//       shc_utotal = temp / u_to_temp_fac;

double shc_ene1 = temp / u_to_temp_fac;
fprintf (diag_shc,"ENE_PER_GRAM : %i %g %g %g %g\n",idx,rho,temp,shc_utotal,shc_ene1);
fflush (diag_shc);
*/

return shc_utotal;
}

double shc_uh2rot(double temp)
{

int j,k;
int shc_jmax;
double shc_zp,shc_fp,shc_zo,shc_fo,shc_zpi,shc_zoi,shc_v,shc_jjp1,shc_kkp1,shc_gj,shc_gk;
double shc_thetarot,shc_ratiop,shc_ratioo;
double shc_enerot;

shc_thetarot = 85.4;
shc_ratiop = 1.;
shc_ratioo = 3.;

shc_zp = 0.;
shc_fp = 0.;
shc_zo = 0.;
shc_fo = 0.;

shc_v = shc_thetarot / temp;
shc_jmax = ( ((int) (sqrt(200./shc_v))) < 20 ) ? ((int) (sqrt(200./shc_v))):20;
shc_v = -1.*shc_v;

for (j=0;j<=shc_jmax;j=j+2)
{
     k = j + 1;
     shc_jjp1 = (double) (j*(j+1));
     shc_kkp1 = (double) (k*(k+1));

     shc_gj = 2.*j + 1.;
     shc_gk = shc_gj + 2.;

     shc_zpi = shc_gj * exp(shc_jjp1*shc_v);
     shc_zp = shc_zp + shc_zpi;
     shc_fp = shc_fp + shc_jjp1 * shc_zpi;

     shc_zoi = shc_gk * exp(shc_kkp1*shc_v);
     shc_zo = shc_zo + shc_zoi;
     shc_fo = shc_fo + shc_kkp1 * shc_zoi;
}

shc_enerot = shc_thetarot*(shc_ratiop*shc_fp + shc_ratioo*shc_fo) / (shc_ratiop*shc_zp + shc_ratioo*shc_zo);

return shc_enerot;
}

#endif /* REAL_EOS */
