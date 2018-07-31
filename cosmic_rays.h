#ifdef COSMIC_RAYS
#ifndef COSMIC_RAYS_H
#define COSMIC_RAYS_H

/* Sanity checks for compiler switches in conjunction with cosmic
 * rays
 */

#if defined(CR_SHOCK) && !defined(COOLING) 
// #error Cannot compile with CR_SHOCK but without cooling.
#endif

typedef struct sph_particle_data SphParticle;

/* ============================================================ */
/* ============ Interface functions to GADGET ================= */
/* ============================================================ */

int extern CR_initialize_beta_tabs( double Alpha );
void extern CR_free_beta_tabs ( void );
double Beta(double A, double B, double x);

double CR_Tab_MeanEnergy(double q, double alpha);

void extern CR_Particle_Update( SphParticle* Particle );
double CR_Particle_GetCoolingTimescale(double q, double rho);
double CR_Physical_Pressure(SphParticle * Particle);

double extern CR_Particle_Pressure( SphParticle* Particle );
double extern CR_Comoving_Pressure( SphParticle* Particle );
double extern CR_Physical_Pressure( SphParticle* Particle );

double extern CR_get_energy(SphParticle * Particle, double q0);

double extern CR_Particle_SpecificEnergy( SphParticle* Particle );
double extern CR_Comoving_SpecificEnergy( SphParticle* Particle );
double extern CR_Physical_SpecificEnergy( SphParticle* Particle );

double extern CR_Particle_SpecificNumber( SphParticle* Particle );

double extern CR_Particle_BaryonFraction( SphParticle* Particle );
double extern CR_Particle_MeanKineticEnergy( SphParticle* Particle );

double CR_Particle_GetCooolingTimescale(double q0, double rDensity);

double CR_Particle_GetQforCoolingTimeScale(SphParticle * Particle, double TimeScale, double q0);
double CR_Particle_GetCoulombTimescale(double rQ, double rDensity);

double CR_Particle_GetThermalizationTimescale(double rQ, double rDensity);
double CR_Particle_GetDissipationTimescale(double rQ, double rDensity);
double CR_Particle_GetCooolingTimescale(double rQ, double rDensity);

double CR_EnergyFractionAfterShiftingQ(double Q1, double Q2, double SpectralIndex);
double CR_Particle_GetThermalizationTimescale(double rQ, double rDensity);
void CR_Tab_Initialize(void);
double CR_Tab_GetThermalizationTimescale(double q, double rho);
double CR_Tab_GetDissipationTimescale(double q, double rho);
double CR_Tab_GetCoolingTimescale(double q, double rho);
double  CR_Find_Qinj(double rho, double TimeScale, double alpha, double qstart);
double CR_corresponding_q(double qinj, double alpha);
double CR_Tab_Beta(double q);


void   CR_test_routine(void);


#ifdef COOLING
void extern CR_Particle_GetThermalizationRate( SphParticle* Particle, double q0, double* SpecEnergyChangeRate, double* SpecNumberChangeRate );
#endif

void extern CR_Particle_GetDissipationRate( SphParticle* Particle, double* SpecEnergyChangeRate, double* BaryonFractionChangeRate );

double CR_Particle_ThermalizeAndDissipate( SphParticle *Particle,
					   double Time );


#ifdef SFR
double extern CR_Particle_SupernovaFeedback( SphParticle* Particle, double SpecificEnergyInput, double TimeScale );
#endif

#if defined(CR_SHOCK) || defined(CR_SHOCK)
double extern CR_Particle_ShockInject( SphParticle* Particle,
				       double SpecificEnergyInput,
				       double TimeScale );
#endif

void extern CR_Particle_Inject( SphParticle* Particle, double DeltaE, double DeltaN );

double extern CR_q_from_mean_kinetic_energy( double T );
double extern CR_mean_kinetic_energy( double q, double Alpha );


double extern CR_Particle_TimeScaleQ( SphParticle *Particle, double TimeScale );

#endif
#endif




