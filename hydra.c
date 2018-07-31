#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>


#include "allvars.h"
#include "proto.h"
#ifdef COSMIC_RAYS
#include "cosmic_rays.h"
#endif
#ifdef MACH_NUM
#include "machfinder.h"
#endif

#ifndef DEBUG
#define NDEBUG
#endif
#include <assert.h>

/*! \file hydra.c
*  \brief Computation of SPH forces and rate of entropy generation
*
*  This file contains the "second SPH loop", where the SPH forces are
*  computed, and where the rate of change of entropy due to the shock heating
*  (via artificial viscosity) is computed.
*/


struct hydrodata_in
{
  MyDouble Pos[3];
  MyFloat Vel[3];
  MyFloat Hsml;
  MyFloat Mass;
  MyFloat Density;
  MyFloat Pressure;
  MyFloat F1;
  MyFloat DhsmlDensityFactor;
  int Timestep;

#ifdef PARTICLE_DEBUG
  MyIDType ID;			/*!< particle identifier */
#endif

#ifdef MAGNETIC
  MyFloat BPred[3];
#ifdef TIME_DEP_MAGN_DISP
  MyFloat Balpha;
#endif
#ifdef DIVBCLEANING_DEDNER
  MyFloat PhiPred;
#endif
#endif
#ifdef TIME_DEP_ART_VISC
  MyFloat alpha;
#endif

#if defined(CONDUCTION) || defined(CR_DIFFUSION) || defined(NAVIERSTOKES)
  MyFloat Entropy;
#endif

#ifdef CONDUCTION
  MyFloat SmoothedEntr;
#ifdef CONDUCTION_SATURATION
  MyFloat GradEntr[3];
#endif
#endif

#ifdef CR_DIFFUSION
  MyFloat CR_E0;		/*!< diffusion term for cosmic ray energy */
  MyFloat CR_n0;		/*!< diffusion term for cosmic ray baryon fraction */
  MyFloat CR_q0;
#endif

#ifdef NAVIERSTOKES
  MyFloat stressoffdiag[3];
  MyFloat stressdiag[3];
  MyFloat shear_viscosity;
#endif

#ifdef NAVIERSTOKES_BULK
  MyFloat divvel;
#endif

  int NodeList[NODELISTLENGTH];
}
 *HydroDataIn, *HydroDataGet;


struct hydrodata_out
{
  MyLongDouble Acc[3];
  MyLongDouble DtEntropy;
#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
  MyFloat MinViscousDt;
#else
  MyFloat MaxSignalVel;
#endif
#ifdef MAGNETIC
  MyFloat DtB[3];
#if defined(TRACEDIVB) || defined(TIME_DEP_MAGN_DISP)
  MyFloat divB;
#endif
#ifdef DIVBCLEANING_DEDNER
  MyFloat DtPhi;
#endif
#endif
#ifdef CONDUCTION
  MyFloat CondEnergyChange;
#ifdef OUTPUTCOOLRATE
  MyFloat CondRate;
#endif
#endif

#if  defined(CR_SHOCK) || defined(CR_DIFFUSION)
  MyFloat CR_EnergyChange;
  MyFloat CR_BaryonFractionChange;
#endif

}
 *HydroDataResult, *HydroDataOut;






#ifdef MACHNUM
double hubble_a, atime, hubble_a2, fac_mu, fac_vsic_fix, a3inv, fac_egy;
#else
static double hubble_a, atime, hubble_a2, fac_mu, fac_vsic_fix, a3inv, fac_egy;
#endif

#if defined(MAGNETIC) && defined(MAGFORCE) && defined(ARTBPRES)
static double u1;
#endif

/*! This function is the driver routine for the calculation of hydrodynamical
*  force and rate of change of entropy due to shock heating for all active
*  particles .
*/
void hydro_force(void)
{
  int i, j, k, ngrp, ndone, ndone_flag, dummy;
  int sendTask, recvTask, nexport, nimport, place;
  double soundspeed_i;
  double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 =
    0, timewait2 = 0;
  double timecomp, timecomm, timewait, tstart, tend, t0, t1;

#ifdef NAVIERSTOKES
  double fac;
#endif
#ifdef CONDUCTION
  double rEntropyFlow, dt, fac;
#endif

#if defined(CR_SHOCK)
  double rShockEnergy;
  double rNonRethermalizedEnergy;

#ifndef COOLING
  double utherm;
#endif
#endif

#ifdef WINDS
  double windspeed, hsml_c;
  double dmax1, dmax2;
#endif


#ifdef TIME_DEP_ART_VISC
  double f, cs_h;
#endif
#if defined(MAGNETIC) && defined(MAGFORCE)
#ifdef TIME_DEP_MAGN_DISP
  double mu0 = 1;
#endif
#ifdef ARTBPRES
  /* mean particle placing */
  /*
     u1 = pow(4. * M_PI / 3 / All.DesNumNgb, 1. / 3.);
   */
  u1 = 0.3;
#endif
#endif

  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
      hubble_a = hubble_function(All.Time);
      hubble_a2 = All.Time * All.Time * hubble_a;

      fac_mu = pow(All.Time, 3 * (GAMMA - 1) / 2) / All.Time;

      fac_egy = pow(All.Time, 3 * (GAMMA - 1));

      fac_vsic_fix = hubble_a * pow(All.Time, 3 * GAMMA_MINUS1);

      a3inv = 1 / (All.Time * All.Time * All.Time);
      atime = All.Time;
    }
  else
    hubble_a = hubble_a2 = atime = fac_mu = fac_vsic_fix = a3inv = fac_egy = 1.0;

#if defined(MAGFORCE) && defined(TIME_DEP_MAGN_DISP)
#ifndef MU0_UNITY
  mu0 *= (4 * M_PI);
  mu0 /= All.UnitTime_in_s * All.UnitTime_in_s *
    All.UnitLength_in_cm / (All.UnitMass_in_g * All.HubbleParam * All.HubbleParam);
#endif
#endif


#ifdef CONDUCTION
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(P[i].Type == 0)
	SphP[i].CondEnergyChange = 0.0;
    }
#endif


  /* allocate buffers to arrange communication */

  Ngblist = (int *) mymalloc(NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct hydrodata_in) +
					     sizeof(struct hydrodata_out) +
					     sizemax(sizeof(struct hydrodata_in),
						     sizeof(struct hydrodata_out))));
  DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));


  CPU_Step[CPU_HYDMISC] += measure_time();
  t0 = second();

  i = FirstActiveParticle;	/* first particle for this task */

  do
    {
      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */
      tstart = second();
      for(nexport = 0; i >= 0; i = NextActiveParticle[i])
	if(P[i].Type == 0)
	  {
	    if(hydro_evaluate(i, 0, &nexport, Send_count) < 0)
	      break;
	  }
      tend = second();
      timecomp1 += timediff(tstart, tend);

#ifdef MYSORT
      mysort_dataindex(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#else
      qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#endif

      tstart = second();

      MPI_Allgather(Send_count, NTask, MPI_INT, Sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);

      tend = second();
      timewait1 += timediff(tstart, tend);

      for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	{
	  Recv_count[j] = Sendcount_matrix[j * NTask + ThisTask];
	  nimport += Recv_count[j];

	  if(j > 0)
	    {
	      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	    }
	}

      HydroDataGet = (struct hydrodata_in *) mymalloc(nimport * sizeof(struct hydrodata_in));
      HydroDataIn = (struct hydrodata_in *) mymalloc(nexport * sizeof(struct hydrodata_in));

      /* prepare particle data for export */

      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  for(k = 0; k < 3; k++)
	    {
	      HydroDataIn[j].Pos[k] = P[place].Pos[k];
	      HydroDataIn[j].Vel[k] = SphP[place].VelPred[k];
	    }
	  HydroDataIn[j].Hsml = PPP[place].Hsml;
	  HydroDataIn[j].Mass = P[place].Mass;
	  HydroDataIn[j].DhsmlDensityFactor = SphP[place].h.DhsmlDensityFactor;
	  HydroDataIn[j].Density = SphP[place].d.Density;
	  HydroDataIn[j].Pressure = SphP[place].Pressure;
	  HydroDataIn[j].Timestep = (P[place].TimeBin ? (1 << P[place].TimeBin) : 0);

	  /* calculation of F1 */
#ifndef ALTVISCOSITY
	  soundspeed_i = sqrt(GAMMA * SphP[place].Pressure / SphP[place].d.Density);
#ifndef NAVIERSTOKES
	  HydroDataIn[j].F1 = fabs(SphP[place].v.DivVel) /
	    (fabs(SphP[place].v.DivVel) + SphP[place].r.CurlVel +
	     0.0001 * soundspeed_i / PPP[place].Hsml / fac_mu);
#else
	  HydroDataIn[j].F1 = fabs(SphP[place].v.DivVel) /
	    (fabs(SphP[place].v.DivVel) + SphP[place].u.s.CurlVel +
	     0.0001 * soundspeed_i / PPP[place].Hsml / fac_mu);
#endif

#else
	  HydroDataIn[j].F1 = SphP[place].v.DivVel;
#endif

	  memcpy(HydroDataIn[j].NodeList,
		 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));


#ifdef MAGNETIC
	  for(k = 0; k < 3; k++)
	    HydroDataIn[j].BPred[k] = SphP[place].BPred[k];
#ifdef DIVBCLEANING_DEDNER
#ifdef SMOOTH_PHI
	  HydroDataIn[j].PhiPred = SphP[place].SmoothPhi;
#else
	  HydroDataIn[j].PhiPred = SphP[place].PhiPred;
#endif
#endif
#endif


#if defined(CONDUCTION) || defined(CR_DIFFUSION) || defined(NAVIERSTOKES)
	  HydroDataIn[j].Entropy = SphP[place].Entropy;
#endif

#ifdef CONDUCTION
	  HydroDataIn[j].SmoothedEntr = SphP[place].SmoothedEntr;
#ifdef CONDUCTION_SATURATION
	  HydroDataIn[j].GradEntr[0] = SphP[place].GradEntr[0];
	  HydroDataIn[j].GradEntr[1] = SphP[place].GradEntr[1];
	  HydroDataIn[j].GradEntr[2] = SphP[place].GradEntr[2];
#endif
#endif

#ifdef TIME_DEP_ART_VISC
	  HydroDataIn[j].alpha = SphP[place].alpha;
#endif

#ifdef SFR_DECOUPLING
	  HydroDataIn[j].DensityOld = SphP[place].d.Density;
	  HydroDataIn[j].Entropy = SphP[place].Entropy;
#endif


#ifdef PARTICLE_DEBUG
	  HydroDataIn[j].ID = P[place].ID;
#endif

#ifdef CR_DIFFUSION
	  HydroDataIn[j].CR_E0 = SphP[place].CR_E0;
	  HydroDataIn[j].CR_n0 = SphP[place].CR_n0;
	  HydroDataIn[j].CR_q0 = SphP[place].CR_q0;
#endif

#ifdef NAVIERSTOKES
	  for(k = 0; k < 3; k++)
	    {
	      HydroDataIn[j].stressdiag[k] = SphP[i].u.s.StressDiag[k];
	      HydroDataIn[j].stressoffdiag[k] = SphP[i].u.s.StressOffDiag[k];
	    }
	  HydroDataIn[j].shear_viscosity = get_shear_viscosity(i);

#ifdef NAVIERSTOKES_BULK
	  HydroDataIn[j].divvel = SphP[i].u.s.DivVel;
#endif
#endif

#ifdef TIME_DEP_MAGN_DISP
	  HydroDataIn[j].Balpha = SphP[place].Balpha;
#endif
	}




      /* exchange particle data */
      tstart = second();
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* get the particles */
		  MPI_Sendrecv(&HydroDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct hydrodata_in), MPI_BYTE,
			       recvTask, TAG_HYDRO_A,
			       &HydroDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct hydrodata_in), MPI_BYTE,
			       recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}
      tend = second();
      timecommsumm1 += timediff(tstart, tend);


      myfree(HydroDataIn);
      HydroDataResult = (struct hydrodata_out *) mymalloc(nimport * sizeof(struct hydrodata_out));
      HydroDataOut = (struct hydrodata_out *) mymalloc(nexport * sizeof(struct hydrodata_out));



      /* now do the particles that were sent to us */

      tstart = second();
      for(j = 0; j < nimport; j++)
	{
	  hydro_evaluate(j, 1, &dummy, &dummy);
	}
      tend = second();
      timecomp2 += timediff(tstart, tend);

      if(i < 0)
	ndone_flag = 1;
      else
	ndone_flag = 0;

      tstart = second();
      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      tend = second();
      timewait2 += timediff(tstart, tend);


      /* get the result */
      tstart = second();
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;
	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* send the results */
		  MPI_Sendrecv(&HydroDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct hydrodata_out),
			       MPI_BYTE, recvTask, TAG_HYDRO_B,
			       &HydroDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct hydrodata_out),
			       MPI_BYTE, recvTask, TAG_HYDRO_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}
      tend = second();
      timecommsumm2 += timediff(tstart, tend);



      /* add the result to the local particles */
      tstart = second();
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  for(k = 0; k < 3; k++)
	    {
	      SphP[place].a.dHydroAccel[k] += HydroDataOut[j].Acc[k];
	    }

	  SphP[place].e.dDtEntropy += HydroDataOut[j].DtEntropy;

#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
	  if(SphP[place].MinViscousDt > HydroDataOut[j].MinViscousDt)
	    SphP[place].MinViscousDt = HydroDataOut[j].MinViscousDt;
#else
	  if(SphP[place].MaxSignalVel < HydroDataOut[j].MaxSignalVel)
	    SphP[place].MaxSignalVel = HydroDataOut[j].MaxSignalVel;
#endif

#ifdef CONDUCTION
	  SphP[place].CondEnergyChange += HydroDataOut[j].CondEnergyChange;
#ifdef OUTPUTCOOLRATE
	  SphP[place].CondRate += HydroDataOut[j].CondRate;
#endif
#endif
#ifdef MAGNETIC
	  for(k = 0; k < 3; k++)
	    SphP[place].DtB[k] += HydroDataOut[j].DtB[k];
#ifdef DIVBCLEANING_DEDNER
	  SphP[place].DtPhi += HydroDataOut[j].DtPhi;
#endif
#ifdef TRACEDIVB
	  SphP[place].divB += HydroDataOut[j].divB;
#endif
#endif
#if defined(CR_DIFFUSION)
	  CR_Particle_Inject(SphP + place,
			     HydroDataOut[j].CR_EnergyChange, HydroDataOut[j].CR_BaryonFractionChange);
#endif
	}
      tend = second();
      timecomp1 += timediff(tstart, tend);

      myfree(HydroDataOut);
      myfree(HydroDataResult);
      myfree(HydroDataGet);
    }
  while(ndone < NTask);


  myfree(DataNodeList);
  myfree(DataIndexTable);

  myfree(Ngblist);


  /* do final operations on results */


#ifdef FLTROUNDOFFREDUCTION
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 0)
      {
	SphP[i].e.DtEntropy = FLT(SphP[i].e.dDtEntropy);

	for(j = 0; j < 3; j++)
	  SphP[i].a.HydroAccel[j] = FLT(SphP[i].a.dHydroAccel[j]);
      }
#endif



  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 0)
      {
#ifdef CR_SHOCK
	/* state right here:
	 *
	 * _c denotes comoving quantities
	 * _p denotes physical quantities
	 *
	 *
	 * Delta u_p = rho_p^(gamma-1)/(gamma-1) Delta A
	 *
	 * Delta A = dA/dloga * Delta loga
	 *
	 * dA/dloga = DtE * (gamma-1) / ( H(a) a^2 rho_c^(gamma-1)
	 *
	 * => Delta u_p = DtE * dloga / ( H(a) a^2 a^(3(gamma-1)) )
	 */

	if(SphP[i].e.DtEntropy > 0.0)
	  {
	    rShockEnergy = SphP[i].e.DtEntropy *
	      (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / hubble_a2 / fac_egy;
	  }
	else
	  {
	    rShockEnergy = 0.0;
	  }

#endif /* CR_SHOCK */

	/* Translate energy change rate into entropy change rate */
	SphP[i].e.DtEntropy *= GAMMA_MINUS1 / (hubble_a2 * pow(SphP[i].d.Density, GAMMA_MINUS1));

#ifdef MACHNUM

	/* Estimates the Mach number of particle i for non-radiative runs,
	 * or the Mach number, density jump and specific energy jump
	 * in case of cosmic rays!
	 */
#if (CR_SHOCK == 2)
	GetMachNumberCR(SphP + i);
#else
	GetMachNumber(SphP + i);
#endif /* COSMIC_RAYS */
#endif /* MACHNUM */
#ifdef MACHSTATISTIC
	GetShock_DtEnergy(SphP + i);
#endif

#ifdef CR_SHOCK
	if(rShockEnergy > 0.0)
	  {
	    /* Feed fraction "All.CR_ShockEfficiency" into CR and see what
	     * amount of energy instantly gets rethermalized
	     *
	     * for this, we need the physical time step, which is
	     * Delta t_p = Delta t_c / hubble_a
	     */

	    rNonRethermalizedEnergy =
	      CR_Particle_ShockInject(SphP + i,
				      rShockEnergy,
				      (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval /
				      hubble_a);

	    /* Fraction of total energy that went and remained in CR is
	     * rNonRethermalizedEnergy / rShockEnergy,
	     * hence, we conserve energy if we do:
	     */
#ifndef CR_NO_CHANGE
	    SphP[i].e.DtEntropy *= (1.0 - rNonRethermalizedEnergy / rShockEnergy);
#endif /* CR_NO_CHANGE */

	    assert(rNonRethermalizedEnergy >= 0.0);

	    assert(rNonRethermalizedEnergy <= (rShockEnergy * All.CR_ShockEfficiency));


#ifndef COOLING
	    utherm = CR_Particle_ThermalizeAndDissipate(SphP + i, (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) *
							All.Timebase_interval / hubble_a);

	    /* we need to add this thermalized energy to the internal energy */

	    SphP[i].e.DtEntropy += GAMMA_MINUS1 * utherm * fac_egy / pow(SphP[i].d.Density, GAMMA_MINUS1) /
	      ((P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval);
#endif

	  }
#endif /* CR_SHOCK */


#ifdef NAVIERSTOKES
	/* sigma_ab * sigma_ab */
	for(k = 0, fac = 0; k < 3; k++)
	  {
	    fac += SphP[i].u.s.StressDiag[k] * SphP[i].u.s.StressDiag[k] +
	      2 * SphP[i].u.s.StressOffDiag[k] * SphP[i].u.s.StressOffDiag[k];
	  }

#ifndef NAVIERSTOKES_CONSTANT	/*entropy increase due to the shear viscosity */
#ifdef NS_TIMESTEP
	SphP[i].ViscEntropyChange = 0.5 * GAMMA_MINUS1 /
	  (hubble_a2 * pow(SphP[i].d.Density, GAMMA_MINUS1)) *
	  get_shear_viscosity(i) / SphP[i].d.Density * fac *
	  pow((SphP[i].Entropy * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1), 2.5);

	SphP[i].e.DtEntropy += SphP[i].ViscEntropyChange;
#else
	SphP[i].e.DtEntropy += 0.5 * GAMMA_MINUS1 /
	  (hubble_a2 * pow(SphP[i].d.Density, GAMMA_MINUS1)) *
	  get_shear_viscosity(i) / SphP[i].d.Density * fac *
	  pow((SphP[i].Entropy * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1), 2.5);
#endif

#else
	SphP[i].e.DtEntropy += 0.5 * GAMMA_MINUS1 /
	  (hubble_a2 * pow(SphP[i].d.Density, GAMMA_MINUS1)) *
	  get_shear_viscosity(i) / SphP[i].d.Density * fac;

#ifdef NS_TIMESTEP
	SphP[i].ViscEntropyChange = 0.5 * GAMMA_MINUS1 /
	  (hubble_a2 * pow(SphP[i].d.Density, GAMMA_MINUS1)) *
	  get_shear_viscosity(i) / SphP[i].d.Density * fac;
#endif

#endif

#ifdef NAVIERSTOKES_BULK	/*entropy increase due to the bulk viscosity */
	SphP[i].e.DtEntropy += GAMMA_MINUS1 /
	  (hubble_a2 * pow(SphP[i].d.Density, GAMMA_MINUS1)) *
	  All.NavierStokes_BulkViscosity / SphP[i].d.Density * pow(SphP[i].u.s.DivVel, 2);

#ifdef NS_TIMESTEP
	SphP[i].ViscEntropyChange = GAMMA_MINUS1 /
	  (hubble_a2 * pow(SphP[i].d.Density, GAMMA_MINUS1)) *
	  All.NavierStokes_BulkViscosity / SphP[i].d.Density * pow(SphP[i].u.s.DivVel, 2);
#endif

#endif

#endif /* these entropy increases directly follow from the general heat transfer equation */


#ifdef MAGNETIC
	/* take care of cosmological dilution */
	if(All.ComovingIntegrationOn)
	  for(k = 0; k < 3; k++)
	    SphP[i].DtB[k] -= 2.0 * SphP[i].BPred[k];
#endif

#ifdef WINDS
	/* if we have winds, we decouple particles briefly if delaytime>0 */

	if(SphP[i].DelayTime > 0)
	  {
	    for(k = 0; k < 3; k++)
	      SphP[i].a.HydroAccel[k] = 0;

	    SphP[i].e.DtEntropy = 0;

#ifdef NOWINDTIMESTEPPING
	    SphP[i].MaxSignalVel = 2 * sqrt(GAMMA * SphP[i].Pressure / SphP[i].d.Density);
#else
	    windspeed = sqrt(2 * All.WindEnergyFraction * All.FactorSN *
			     All.EgySpecSN / (1 - All.FactorSN) / All.WindEfficiency) * All.Time;
	    windspeed *= fac_mu;
	    hsml_c = pow(All.WindFreeTravelDensFac * All.PhysDensThresh /
			 (SphP[i].d.Density * a3inv), (1. / 3.));
	    SphP[i].MaxSignalVel = hsml_c * DMAX((2 * windspeed), SphP[i].MaxSignalVel);
#endif
	  }
#endif

#ifdef SPH_BND_PARTICLES
	if(P[i].ID == 0)
	  {
	    SphP[i].e.DtEntropy = 0;
	    //	    SphP[i].ViscEntropyChange = 0;
	    for(k = 0; k < 3; k++)
	      SphP[i].a.HydroAccel[k] = 0;
	  }
#endif
#ifdef TIME_DEP_ART_VISC
	cs_h = sqrt(GAMMA * SphP[i].Pressure / SphP[i].d.Density) / PPP[i].Hsml;
	f = fabs(SphP[i].u.s.DivVel) /
	  (fabs(SphP[i].u.s.DivVel) + SphP[i].u.s.CurlVel + 0.0001 * cs_h / fac_mu);
	SphP[i].Dtalpha = -(SphP[i].alpha - All.AlphaMin) * All.DecayTime *
	  0.5 * SphP[i].MaxSignalVel / (PPP[i].Hsml * fac_mu)
	  + f * All.ViscSource * DMAX(0.0, -SphP[i].u.s.DivVel);
	if(All.ComovingIntegrationOn)
	  SphP[i].Dtalpha /= (hubble_a * All.Time * All.Time);
#endif
#ifdef MAGNETIC
#ifdef TIME_DEP_MAGN_DISP
	SphP[i].DtBalpha = -(SphP[i].Balpha - All.ArtMagDispMin) * All.ArtMagDispTime *
	  0.5 * SphP[i].MaxSignalVel / (PPP[i].Hsml * fac_mu)
	  + All.ArtMagDispSource * fabs(SphP[i].divB) / sqrt(mu0 * SphP[i].d.Density);
#endif
#ifdef DIVBCLEANING_DEDNER
	SphP[place].DtPhi -=
	  SphP[i].PhiPred * All.DivBcleanParabolicSigma * 0.5 * SphP[i].MaxSignalVel / (PPP[i].Hsml * fac_mu);
#endif
#endif
      }


#ifdef CONDUCTION		/* note: this needs to be done for all particles, not just active ones */
  for(i = 0; i < N_gas; i++)
    if(P[i].Type == 0)
      {
	if(All.ComovingIntegrationOn)
	  SphP[i].CondEnergyChange *= All.Time / hubble_a;

	if(P[i].TimeBin)
	  {
	    rEntropyFlow =
	      SphP[i].CondEnergyChange / P[i].Mass * GAMMA_MINUS1 * pow(SphP[i].d.Density * a3inv,
									-GAMMA_MINUS1);
	    if(fabs(rEntropyFlow) > 0.5 * SphP[i].Entropy)
	      {
		printf
		  ("NOTE: step=%d task=%d large entropy change for particle i=%d (ID=%d) due to conduction: A=%g  dA=%g rho=%g dt=%g\n",
		   All.NumCurrentTiStep, ThisTask, i, (int) P[i].ID, SphP[i].Entropy, rEntropyFlow,
		   SphP[i].d.Density, (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval);
		fflush(stdout);

		fac = 0.5 * SphP[i].Entropy / fabs(rEntropyFlow);

		rEntropyFlow *= fac;
		SphP[i].CondEnergyChange *= fac;
	      }

	    SphP[i].Entropy += rEntropyFlow;

	    if(TimeBinActive[P[i].TimeBin])
	      {
		dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;
		SphP[i].CondEnergyChange /= (P[i].Mass * dt);
	      }
	  }
	else
	  {
	    SphP[i].CondEnergyChange /= P[i].Mass;
	  }
      }
#endif

  /* collect some timing information */

  t1 = WallclockTime = second();
  timeall += timediff(t0, t1);

  timecomp = timecomp1 + timecomp2;
  timewait = timewait1 + timewait2;
  timecomm = timecommsumm1 + timecommsumm2;

  CPU_Step[CPU_HYDCOMPUTE] += timecomp;
  CPU_Step[CPU_HYDWAIT] += timewait;
  CPU_Step[CPU_HYDCOMM] += timecomm;
  CPU_Step[CPU_HYDMISC] += timeall - (timecomp + timewait + timecomm);
}




/*! This function is the 'core' of the SPH force computation. A target
*  particle is specified which may either be local, or reside in the
*  communication buffer.
*/
int hydro_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int startnode, numngb, listindex = 0;
  int j, k, n, timestep;
  MyDouble *pos;
  MyFloat *vel;
  MyFloat mass, h_i, dhsmlDensityFactor, rho, pressure, f1, f2;
  MyLongDouble acc[3], dtEntropy;

#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
  MyFloat minViscousDt;
#else
  MyFloat maxSignalVel;
#endif
  double dx, dy, dz, dvx, dvy, dvz;
  double h_i2, hinv, hinv4;
  double p_over_rho2_i, p_over_rho2_j, soundspeed_i, soundspeed_j;
  double hfc, dwk_i, vdotr, vdotr2, visc, mu_ij, rho_ij, vsig;
  double h_j, dwk_j;
  double r, r2, u;
  double hfc_visc;
  double dmin1, dmin2;
  double BulkVisc_ij;
  int imax1, imax2;

#ifdef NAVIERSTOKES
  double faci, facj;
  MyFloat *stressdiag;
  MyFloat *stressoffdiag;
  MyFloat shear_viscosity;

#ifdef VISCOSITY_SATURATION
  double VelLengthScale_i, VelLengthScale_j;
  double IonMeanFreePath_i, IonMeanFreePath_j;
#endif
#ifdef NAVIERSTOKES_BULK
  double facbi, facbj;
  MyFloat divvel;
#endif
#endif

#if defined(CONDUCTION) || defined(CR_DIFFUSION) || defined(NAVIERSTOKES)
  double Entropy;
#endif

#ifdef CONDUCTION
  double condEnergyChange;
  double rEnergy_i, rEnergy_j, rEnergy_i_ns, rEnergy_j_ns;
  double rKappa_i, rKappa_j;	/* Thermal Conduction Coefficients */
  double rEnergyFlow;
  double smoothentr;
  double dtcond;
  double rEntropy2Energy_i;

#ifdef CONDUCTION_SATURATION
  double gradentr[3];
  double electron_free_path, temp_scale_length;
#endif
#endif

#ifdef CR_DIFFUSION
  double CR_E0_i, CR_n0_i;
  double CR_q_i, CR_q_j;
  double rKappaCR_i, rKappaCR_j;
  double cr_efac_i, cr_efac_j;
  double rCR_DiffusionTerm;
  double rCR_EnergyFlow;
  double rCR_MassFlow;
  double dtdiff;
#endif


#ifdef TIME_DEP_ART_VISC
  MyFloat alpha;
#endif

#ifdef ALTVISCOSITY
  double mu_i, mu_j;
#endif

#ifndef NOVISCOSITYLIMITER
  double dt;
#endif

#ifdef MAGNETIC
  MyFloat *bpred;
  double dtB[3];
  double dBx, dBy, dBz;
  double magfac, magfac_i, magfac_j, magfac_i_base;
  double mu0_1;

#ifdef MAGFORCE
  double mm_i[3][3], mm_j[3][3];
  double b2_i, b2_j;
  int l;

#ifdef ARTBPRES
  double wk_i, wk_j, hinv3;
  double w1_i, w1_j, R_abp_i, R_abp_j;
  double mm_abp_i[3][3], mm_abp_j[3][3];
#endif
#endif
#if defined(MAGNETIC_DISSIPATION) || defined(DIVBCLEANING_DEDNER)
  double magfac_sym;
#endif
#ifdef MAGNETIC_DISSIPATION
  double dTu_diss_b, Balpha_ij;

#ifdef MAGDISSIPATION_PERPEN
  double mft, mvt[3];
#endif
#ifdef TIME_DEP_MAGN_DISP
  double Balpha;
#endif
#endif
#ifdef DIVBCLEANING_DEDNER
  double PhiPred, DtPhi, phifac;
#endif
#ifdef TRACEDIVB
  double divB;
#endif
#ifdef MAGNETIC_SIGNALVEL
  double magneticspeed_i, magneticspeed_j, vcsa2_i, vcsa2_j, Bpro2_i, Bpro2_j;
#endif
#endif

#ifdef PARTICLE_DEBUG
  MyIDType ID;			/*!< particle identifier */
#endif

#ifdef CONVENTIONAL_VISCOSITY
  double c_ij, h_ij;
#endif



  if(mode == 0)
    {
      pos = P[target].Pos;
      vel = SphP[target].VelPred;
      h_i = PPP[target].Hsml;
      mass = P[target].Mass;
      dhsmlDensityFactor = SphP[target].h.DhsmlDensityFactor;
      rho = SphP[target].d.Density;
      pressure = SphP[target].Pressure;
      timestep = (P[target].TimeBin ? (1 << P[target].TimeBin) : 0);
      soundspeed_i = sqrt(GAMMA * pressure / rho);

#ifndef ALTVISCOSITY
#ifndef NAVIERSTOKES
      f1 = fabs(SphP[target].v.DivVel) /
	(fabs(SphP[target].v.DivVel) + SphP[target].r.CurlVel +
	 0.0001 * soundspeed_i / PPP[target].Hsml / fac_mu);
#else
      f1 = fabs(SphP[target].v.DivVel) /
	(fabs(SphP[target].v.DivVel) + SphP[target].u.s.CurlVel +
	 0.0001 * soundspeed_i / PPP[target].Hsml / fac_mu);
#endif
#else
      f1 = SphP[target].v.DivVel;
#endif
#ifdef MAGNETIC
      bpred = SphP[target].BPred;
#ifdef DIVBCLEANING_DEDNER
#ifdef SMOOTH_PHI
      PhiPred = SphP[target].SmoothPhi;
#else
      PhiPred = SphP[target].PhiPred;
#endif
#endif
#endif
#ifdef TIME_DEP_ART_VISC
      alpha = SphP[target].alpha;
#endif

#if defined(CONDUCTION) || defined(CR_DIFFUSION)|| defined(NAVIERSTOKES)
      Entropy = SphP[target].Entropy;
#endif

#ifdef CONDUCTION
      smoothentr = SphP[target].SmoothedEntr;
#ifdef CONDUCTION_SATURATION
      gradentr[0] = SphP[target].GradEntr[0];
      gradentr[1] = SphP[target].GradEntr[1];
      gradentr[2] = SphP[target].GradEntr[2];
#endif
#endif

#ifdef PARTICLE_DEBUG
      ID = P[target].ID;
#endif

#if CR_DIFFUSION
      CR_E0_i = SphP[target].CR_E0;
      CR_n0_i = SphP[target].CR_n0;
      CR_q_i = SphP[target].CR_q0 * pow(SphP[target].d.Density * a3inv, 0.33333);
      cr_efac_i = CR_Tab_MeanEnergy(CR_q_i, All.CR_Alpha - 0.3333) / CR_Tab_MeanEnergy(CR_q_i, All.CR_Alpha);
#endif

#ifdef NAVIERSTOKES
      stressdiag = SphP[target].u.s.StressDiag;
      stressoffdiag = SphP[target].u.s.StressOffDiag;
      shear_viscosity = get_shear_viscosity(target);
#ifdef NAVIERSTOKES_BULK
      divvel = SphP[target].u.s.DivVel;
#endif
#endif

#ifdef TIME_DEP_MAGN_DISP
      Balpha = SphP[target].Balpha;
#endif
    }
  else
    {
      pos = HydroDataGet[target].Pos;
      vel = HydroDataGet[target].Vel;
      h_i = HydroDataGet[target].Hsml;
      mass = HydroDataGet[target].Mass;
      dhsmlDensityFactor = HydroDataGet[target].DhsmlDensityFactor;
      rho = HydroDataGet[target].Density;
      pressure = HydroDataGet[target].Pressure;
      timestep = HydroDataGet[target].Timestep;
      soundspeed_i = sqrt(GAMMA * pressure / rho);
      f1 = HydroDataGet[target].F1;
#ifdef MAGNETIC
      bpred = HydroDataGet[target].BPred;
#ifdef DIVBCLEANING_DEDNER
      PhiPred = HydroDataGet[target].PhiPred;
#endif
#endif
#ifdef TIME_DEP_ART_VISC
      alpha = HydroDataGet[target].alpha;
#endif

#if defined(CONDUCTION) || defined(CR_DIFFUSION) || defined(NAVIERSTOKES)
      Entropy = HydroDataGet[target].Entropy;
#endif

#ifdef CONDUCTION
      smoothentr = HydroDataGet[target].SmoothedEntr;
#ifdef CONDUCTION_SATURATION
      gradentr[0] = HydroDataGet[target].GradEntr[0];
      gradentr[1] = HydroDataGet[target].GradEntr[1];
      gradentr[2] = HydroDataGet[target].GradEntr[2];
#endif
#endif


#ifdef PARTICLE_DEBUG
      ID = HydroDataGet[target].ID;
#endif

#if CR_DIFFUSION
      CR_E0_i = HydroDataGet[target].CR_E0;
      CR_n0_i = HydroDataGet[target].CR_n0;
      CR_q_i = HydroDataGet[target].CR_q0 * pow(rho * a3inv, 0.33333);
      cr_efac_i = CR_Tab_MeanEnergy(CR_q_i, All.CR_Alpha - 0.3333) / CR_Tab_MeanEnergy(CR_q_i, All.CR_Alpha);
#endif

#ifdef NAVIERSTOKES
      stressdiag = HydroDataGet[target].stressdiag;
      stressoffdiag = HydroDataGet[target].stressoffdiag;
      shear_viscosity = HydroDataGet[target].shear_viscosity;
#endif
#ifdef NAVIERSTOKES
      stressdiag = HydroDataGet[target].stressdiag;
      stressoffdiag = HydroDataGet[target].stressoffdiag;
      shear_viscosity = HydroDataGet[target].shear_viscosity;
#ifdef NAVIERSTOKES_BULK
      divvel = HydroDataGet[target].divvel;
#endif
#endif
#ifdef TIME_DEP_MAGN_DISP
      Balpha = HydroDataGet[target].Balpha;
#endif
    }


  /* initialize variables before SPH loop is started */

  acc[0] = acc[1] = acc[2] = dtEntropy = 0;


#ifdef CONDUCTION
  condEnergyChange = 0;
  dtcond = timestep * All.Timebase_interval;
  if(timestep == 0)
    dtcond = 1;
#endif

#ifdef CR_DIFFUSION
  dtdiff = timestep * All.Timebase_interval;
  if(All.ComovingIntegrationOn)
    dtdiff *= All.Time / hubble_a;
#endif


#ifdef MAGNETIC
  for(k = 0; k < 3; k++)
    dtB[k] = 0;
#ifdef DIVBCLEANING_DEDNER
  DtPhi = 0;
#endif
#ifdef TRACEDIVB
  divB = 0;
#endif
#ifdef MAGFORCE
  magfac_i_base = 1 / (rho * rho);
  mu0_1 = 1;
#ifndef MU0_UNITY
  magfac_i_base /= (4 * M_PI);
  mu0_1 /= (4 * M_PI);
  mu0_1 *= All.UnitTime_in_s * All.UnitTime_in_s *
    All.UnitLength_in_cm / (All.UnitMass_in_g * All.HubbleParam * All.HubbleParam);
#endif
#ifdef CORRECTBFRC
  magfac_i_base *= dhsmlDensityFactor;
#endif
  for(k = 0, b2_i = 0; k < 3; k++)
    {
      b2_i += bpred[k] * bpred[k];
      for(l = 0; l < 3; l++)
	mm_i[k][l] = bpred[k] * bpred[l];
    }
  for(k = 0; k < 3; k++)
    mm_i[k][k] -= 0.5 * b2_i;
#ifdef MAGNETIC_SIGNALVEL
  vcsa2_i = soundspeed_i * soundspeed_i + mu0_1 * b2_i / rho;
#endif
#endif /* end of MAGFORCE */
#endif /* end of MAGNETIC */

  p_over_rho2_i = pressure / (rho * rho);
  p_over_rho2_i *= dhsmlDensityFactor;
  h_i2 = h_i * h_i;

#if defined(CONDUCTION)
  rEntropy2Energy_i = pow(rho * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;
#endif

#ifdef CONDUCTION
  /* this gives the thermal energy per unit mass for particle i */
  rEnergy_i = smoothentr * rEntropy2Energy_i;
  rEnergy_i_ns = Entropy * rEntropy2Energy_i;

#  ifdef CONDUCTION_CONSTANT
  rKappa_i = All.ConductionCoeff;
#  else
  rKappa_i = All.ConductionCoeff * pow(rEnergy_i_ns, 2.5);
#    ifdef CONDUCTION_SATURATION
  electron_free_path = All.ElectronFreePathFactor * rEnergy_i * rEnergy_i / (rho * a3inv);
  temp_scale_length =
    atime * fabs(smoothentr) / sqrt(gradentr[0] * gradentr[0] +
				    gradentr[1] * gradentr[1] + gradentr[2] * gradentr[2]);
  rKappa_i /= (1 + 4.2 * electron_free_path / temp_scale_length);
#    endif
#  endif
#  ifdef SFR
  if(rho * a3inv >= All.PhysDensThresh)
    rKappa_i = 0;
#  endif
#endif



#ifdef CR_DIFFUSION
  rKappaCR_i = All.CR_DiffusionCoeff;

  if(All.CR_DiffusionDensScaling != 0.0)
    {
      rKappaCR_i *= pow(rho * a3inv / All.CR_DiffusionDensZero, All.CR_DiffusionDensScaling);
    }

  if(All.CR_DiffusionEntropyScaling != 0.0)
    {
      rKappaCR_i *= pow(Entropy / All.CR_DiffusionEntropyZero, All.CR_DiffusionEntropyScaling);
    }
#endif



#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
  minViscousDt = 1.0e32;
#else
  maxSignalVel = 0;
#endif


  /* Now start the actual SPH computation for this particle */

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = HydroDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  numngb = ngb_treefind_pairs(pos, h_i, target, &startnode, mode, nexport, nsend_local);

	  if(numngb < 0)
	    return -1;

	  for(n = 0; n < numngb; n++)
	    {
	      j = Ngblist[n];

#ifdef BLACK_HOLES
	      if(P[j].Mass == 0)
		continue;
#endif

#ifdef NOWINDTIMESTEPPING
#ifdef WINDS
	      if(P[j].Type == 0)
		if(SphP[j].DelayTime > 0)	/* ignore the wind particles */
		  continue;
#endif
#endif
	      dx = pos[0] - P[j].Pos[0];
	      dy = pos[1] - P[j].Pos[1];
	      dz = pos[2] - P[j].Pos[2];
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
	      h_j = PPP[j].Hsml;
	      if(r2 < h_i2 || r2 < h_j * h_j)
		{
		  r = sqrt(r2);
		  if(r > 0)
		    {
		      p_over_rho2_j = SphP[j].Pressure / (SphP[j].d.Density * SphP[j].d.Density);
		      soundspeed_j = sqrt(GAMMA * p_over_rho2_j * SphP[j].d.Density);
		      dvx = vel[0] - SphP[j].VelPred[0];
		      dvy = vel[1] - SphP[j].VelPred[1];
		      dvz = vel[2] - SphP[j].VelPred[2];
		      vdotr = dx * dvx + dy * dvy + dz * dvz;

		      if(All.ComovingIntegrationOn)
			vdotr2 = vdotr + hubble_a2 * r2;
		      else
			vdotr2 = vdotr;

		      if(r2 < h_i2)
			{
			  hinv = 1.0 / h_i;
#ifndef  TWODIMS
			  hinv4 = hinv * hinv * hinv * hinv;
#else
			  hinv4 = hinv * hinv * hinv / boxSize_Z;
#endif
			  u = r * hinv;
			  if(u < 0.5)
			    dwk_i = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
			  else
			    dwk_i = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
#if defined(MAGNETIC) && defined(MAGFORCE) && defined(ARTBPRES)
#ifndef  TWODIMS
			  hinv3 = hinv * hinv * hinv;
#else
			  hinv3 = hinv * hinv / boxSize_Z;
#endif
			  if(u <= 0.5)
			    wk_i = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
			  else
			    wk_i = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
			  if(u1 <= 0.5)
			    w1_i = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u1 - 1) * u1 * u1);
			  else
			    w1_i = hinv3 * KERNEL_COEFF_5 * (1.0 - u1) * (1.0 - u1) * (1.0 - u1);
#endif
			}
		      else
			{
			  dwk_i = 0;
#if defined(MAGNETIC) && defined(MAGFORCE) && defined(ARTBPRES)
			  wk_i = 0;
			  w1_i = 1;
#endif
			}

		      if(r2 < h_j * h_j)
			{
			  hinv = 1.0 / h_j;
#ifndef  TWODIMS
			  hinv4 = hinv * hinv * hinv * hinv;
#else
			  hinv4 = hinv * hinv * hinv / boxSize_Z;
#endif
			  u = r * hinv;
			  if(u < 0.5)
			    dwk_j = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
			  else
			    dwk_j = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
#if defined(MAGNETIC) && defined(MAGFORCE) && defined(ARTBPRES)
#ifndef  TWODIMS
			  hinv3 = hinv * hinv * hinv;
#else
			  hinv3 = hinv * hinv / boxSize_Z;
#endif
			  if(u <= 0.5)
			    wk_j = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
			  else
			    wk_j = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
			  if(u1 <= 0.5)
			    w1_j = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u1 - 1) * u1 * u1);
			  else
			    w1_j = hinv3 * KERNEL_COEFF_5 * (1.0 - u1) * (1.0 - u1) * (1.0 - u1);
#endif
			}
		      else
			{
			  dwk_j = 0;
#if defined(MAGNETIC) && defined(MAGFORCE) && defined(ARTBPRES)
			  wk_j = 0;
			  w1_j = 1.;
#endif
			}

#ifdef MAGNETIC
#if defined(MAGFORCE) && defined(ARTBPRES)
		      R_abp_i = 0.3 * pow(wk_i / w1_i, 4);
		      R_abp_j = 0.3 * pow(wk_j / w1_j, 4);
#endif
		      dBx = bpred[0] - SphP[j].BPred[0];
		      dBy = bpred[1] - SphP[j].BPred[1];
		      dBz = bpred[2] - SphP[j].BPred[2];
		      magfac = P[j].Mass / r;	/* we moved 'dwk_i / rho' down ! */
#ifdef TRACEDIVB
		      divB += magfac / rho * dwk_i * (dBx * dx + dBy * dy + dBz * dz);
#endif
		      if(All.ComovingIntegrationOn)
			magfac *= 1. / (hubble_a * All.Time * All.Time);
		      /* last factor takes care of all cosmological prefactor */
#ifdef CORRECTDB
		      magfac *= dhsmlDensityFactor;
#endif
#if defined(MAGNETIC_DISSIPATION) || defined(DIVBCLEANING_DEDNER)
		      magfac_sym = 2 * magfac * rho * (dwk_i + dwk_j) /
			((rho + SphP[j].d.Density) * (rho + SphP[j].d.Density));
#endif
#ifdef MAGNETIC_DISSIPATION
#ifdef TIME_DEP_MAGN_DISP
		      Balpha_ij = 0.5 * (Balpha + SphP[j].Balpha);
#else
		      Balpha_ij = All.ArtMagDispConst;
#endif
#endif
		      magfac *= dwk_i / rho;
		      dtB[0] +=
			magfac * ((bpred[0] * dvy - bpred[1] * dvx) * dy +
				  (bpred[0] * dvz - bpred[2] * dvx) * dz);
		      dtB[1] +=
			magfac * ((bpred[1] * dvz - bpred[2] * dvy) * dz +
				  (bpred[1] * dvx - bpred[0] * dvy) * dx);
		      dtB[2] +=
			magfac * ((bpred[2] * dvx - bpred[0] * dvz) * dx +
				  (bpred[2] * dvy - bpred[1] * dvz) * dy);
#ifdef DIVBINDUCTION
		      dtB[0] += magfac * vel[0] * (dBx * dx + dBy * dy + dBz * dz);
		      dtB[1] += magfac * vel[1] * (dBx * dx + dBy * dy + dBz * dz);
		      dtB[2] += magfac * vel[2] * (dBx * dx + dBy * dy + dBz * dz);
#endif
#ifdef DIVBCLEANING_DEDNER
#ifdef SMOOTH_PHI
		      phifac = magfac_sym * (PhiPred - SphP[j].SmoothPhi);
#else
		      phifac = magfac_sym * (PhiPred - SphP[j].PhiPred);
#endif
		      dtB[0] -= phifac * dx;
		      dtB[1] -= phifac * dy;
		      dtB[2] -= phifac * dz;
		      phifac = magfac_sym;
#endif
#ifdef MAGFORCE
		      magfac_j = 1 / (SphP[j].d.Density * SphP[j].d.Density);
#ifndef MU0_UNITY
		      magfac_j /= (4 * M_PI);
#endif
#ifdef CORRECTBFRC
		      magfac_j *= dwk_j * SphP[j].a3.DhsmlDensityFactor;
		      magfac_i = dwk_i * magfac_i_base;
#else
		      magfac_i = magfac_i_base;
#endif
		      for(k = 0, b2_j = 0; k < 3; k++)
			{
			  b2_j += SphP[j].BPred[k] * SphP[j].BPred[k];
			  for(l = 0; l < 3; l++)
			    mm_j[k][l] = SphP[j].BPred[k] * SphP[j].BPred[l];
			}
		      for(k = 0; k < 3; k++)
			mm_j[k][k] -= 0.5 * b2_j;
#ifdef MAGNETIC_SIGNALVEL
		      vcsa2_j = soundspeed_j * soundspeed_j + mu0_1 * b2_j / SphP[j].d.Density;
		      Bpro2_j = (SphP[j].BPred[0] * dx + SphP[j].BPred[1] * dy + SphP[j].BPred[1] * dy) / r;
		      Bpro2_j *= Bpro2_j;
		      magneticspeed_j = sqrt(vcsa2_j +
					     sqrt(DMAX((vcsa2_j * vcsa2_j -
							4 * soundspeed_j * soundspeed_j * Bpro2_j
							* mu0_1 / SphP[j].d.Density), 0))) / 1.4142136;
		      Bpro2_i = (bpred[0] * dx + bpred[1] * dy + bpred[1] * dy) / r;
		      Bpro2_i *= Bpro2_i;
		      magneticspeed_i = sqrt(vcsa2_i +
					     sqrt(DMAX((vcsa2_i * vcsa2_i -
							4 * soundspeed_i * soundspeed_i * Bpro2_i
							* mu0_1 / rho), 0))) / 1.4142136;
#endif
#ifdef MAGNETIC_DISSIPATION
#ifdef MAGDISSIPATION_PERPEN
		      dTu_diss_b = -magfac_sym * Balpha_ij * ((b2_i - Bpro2_i) - (b2_j - Bpro2_j));
#else
		      dTu_diss_b = -magfac_sym * Balpha_ij * (dBx * dBx + dBy * dBy + dBz * dBz);
#endif
#endif
#ifdef CORRECTBFRC
		      magfac = P[j].Mass / r;
#else
		      magfac = P[j].Mass * 0.5 * (dwk_i + dwk_j) / r;
#endif
		      if(All.ComovingIntegrationOn)
			magfac *= pow(All.Time, 3 * GAMMA);
		      /* last factor takes care of all cosmological prefactor */
#ifndef MU0_UNITY
		      magfac *= All.UnitTime_in_s * All.UnitTime_in_s *
			All.UnitLength_in_cm / (All.UnitMass_in_g * All.HubbleParam * All.HubbleParam);
		      /* take care of B unit conversion into GADGET units ! */
#endif
		      for(k = 0; k < 3; k++)
			acc[k] +=
			  magfac * ((mm_i[k][0] * magfac_i + mm_j[k][0] * magfac_j) * dx +
				    (mm_i[k][1] * magfac_i + mm_j[k][1] * magfac_j) * dy +
				    (mm_i[k][2] * magfac_i + mm_j[k][2] * magfac_j) * dz);
#ifdef ARTBPRES
		      for(k = 0; k < 3; k++)
			for(l = 0; l < 3; l++)
			  mm_abp_j[k][l] = SphP[j].BPred[k] * SphP[j].BPred[l] * R_abp_j / 2.0;
		      for(k = 0; k < 3; k++)
			for(l = 0; l < 3; l++)
			  mm_abp_i[k][l] = bpred[k] * bpred[l] * R_abp_i / 2.0;
		      for(k = 0; k < 3; k++)
			acc[k] -= magfac * ((mm_abp_i[k][0] * magfac_i +
					     mm_abp_j[k][0] * magfac_j) * dx +
					    (mm_abp_i[k][1] * magfac_i +
					     mm_abp_j[k][1] * magfac_j) * dy +
					    (mm_abp_i[k][2] * magfac_i + mm_abp_j[k][2] * magfac_j) * dz);
#endif
#ifdef DIVBFORCE
		      for(k = 0; k < 3; k++)
			acc[k] -=
			  magfac * (((bpred[k] * bpred[0]) * magfac_i
				     + (bpred[k] * SphP[j].BPred[0]) * magfac_j) * dx
				    + ((bpred[k] * bpred[1]) * magfac_i
				       + (bpred[k] * SphP[j].BPred[1]) * magfac_j) * dy
				    + ((bpred[k] * bpred[2]) * magfac_i +
				       (bpred[k] * SphP[j].BPred[2]) * magfac_j) * dz);
#endif
#endif
#endif /* end of MAGNETIC */


#ifndef MAGNETIC_SIGNALVEL
		      vsig = soundspeed_i + soundspeed_j;
#else
		      vsig = magneticspeed_i + magneticspeed_j;
#endif
#ifndef ALTERNATIVE_VISCOUS_TIMESTEP
		      if(vsig > maxSignalVel)
			maxSignalVel = vsig;
#endif
		      if(vdotr2 < 0)	/* ... artificial viscosity */
			{
#ifndef ALTVISCOSITY
#ifndef CONVENTIONAL_VISCOSITY
			  mu_ij = fac_mu * vdotr2 / r;	/* note: this is negative! */
#else
			  c_ij = 0.5 * (soundspeed_i + soundspeed_j);
			  h_ij = 0.5 * (h_i + h_j);
			  mu_ij = fac_mu * h_ij * vdotr2 / (r2 + 0.0001 * h_ij * h_ij);
#endif
			  vsig -= 3 * mu_ij;

#ifndef ALTERNATIVE_VISCOUS_TIMESTEP
			  if(vsig > maxSignalVel)
			    maxSignalVel = vsig;
#endif

			  rho_ij = 0.5 * (rho + SphP[j].d.Density);
#ifndef NAVIERSTOKES
			  f2 =
			    fabs(SphP[j].v.DivVel) / (fabs(SphP[j].v.DivVel) + SphP[j].r.CurlVel +
						      0.0001 * soundspeed_j / fac_mu / PPP[j].Hsml);
#else
			  f2 =
			    fabs(SphP[j].v.DivVel) / (fabs(SphP[j].v.DivVel) + SphP[j].u.s.CurlVel +
						      0.0001 * soundspeed_j / fac_mu / PPP[j].Hsml);
#endif

#ifdef NO_SHEAR_VISCOSITY_LIMITER
			  f1 = f2 = 1;
#endif
#ifdef TIME_DEP_ART_VISC
			  BulkVisc_ij = 0.5 * (alpha + SphP[j].alpha);
#else
			  BulkVisc_ij = All.ArtBulkViscConst;
#endif
#ifndef CONVENTIONAL_VISCOSITY
			  visc = 0.25 * BulkVisc_ij * vsig * (-mu_ij) / rho_ij * (f1 + f2);
#else
			  visc =
			    (-BulkVisc_ij * mu_ij * c_ij + 2 * BulkVisc_ij * mu_ij * mu_ij) /
			    rho_ij * (f1 + f2) * 0.5;
#endif

#else /* start of ALTVISCOSITY block */
			  if(f1 < 0)
			    mu_i = h_i * fabs(f1);	/* f1 hold here the velocity divergence of particle i */
			  else
			    mu_i = 0;
			  if(SphP[j].v.DivVel < 0)
			    mu_j = h_j * fabs(SphP[j].v.DivVel);
			  else
			    mu_j = 0;
			  visc = All.ArtBulkViscConst * ((soundspeed_i + mu_i) * mu_i / rho +
							 (soundspeed_j + mu_j) * mu_j / SphP[j].d.Density);
#endif /* end of ALTVISCOSITY block */


			  /* .... end artificial viscosity evaluation */
			  /* now make sure that viscous acceleration is not too large */
#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
			  if(visc > 0)
			    {
			      dt = fac_vsic_fix * vdotr2 /
				(0.5 * (mass + P[j].Mass) * (dwk_i + dwk_j) * r * visc);

			      dt /= hubble_a;

			      if(dt < minViscousDt)
				minViscousDt = dt;
			    }
#endif

#ifndef NOVISCOSITYLIMITER
			  dt =
			    2 * IMAX(timestep,
				     (P[j].TimeBin ? (1 << P[j].TimeBin) : 0)) * All.Timebase_interval;
			  if(dt > 0 && (dwk_i + dwk_j) < 0)
			    {
#ifdef BLACK_HOLE
			      if((mass + P[j].Mass) > 0)
#endif
				visc = DMIN(visc, 0.5 * fac_vsic_fix * vdotr2 /
					    (0.5 * (mass + P[j].Mass) * (dwk_i + dwk_j) * r * dt));
			    }
#endif
			}
		      else
			{
			  visc = 0;
			}

		      p_over_rho2_j *= SphP[j].h.DhsmlDensityFactor;

		      hfc_visc = 0.5 * P[j].Mass * visc * (dwk_i + dwk_j) / r;
		      /* Formulation derived from the Lagrangian */
		      hfc = hfc_visc + P[j].Mass * (p_over_rho2_i * dwk_i + p_over_rho2_j * dwk_j) / r;
#ifdef WINDS
		      if(P[j].Type == 0)
			if(SphP[j].DelayTime > 0)	/* No force by wind particles */
			  {
			    hfc = hfc_visc = 0;
			  }
#endif

#ifndef NOACCEL
		      acc[0] += FLT(-hfc * dx);
		      acc[1] += FLT(-hfc * dy);
		      acc[2] += FLT(-hfc * dz);
#endif

		      dtEntropy += FLT(0.5 * hfc_visc * vdotr2);


#ifdef NAVIERSTOKES
		      faci = mass * shear_viscosity / (rho * rho) * dwk_i / r;

#ifndef NAVIERSTOKES_CONSTANT
		      faci *= pow((Entropy * pow(rho * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1), 2.5);	/*multiplied by E^5/2 */
#endif
		      facj = P[j].Mass * get_shear_viscosity(j) /
			(SphP[j].d.Density * SphP[j].d.Density) * dwk_j / r;

#ifndef NAVIERSTOKES_CONSTANT
		      facj *= pow((SphP[j].Entropy * pow(SphP[j].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1), 2.5);	/*multiplied by E^5/2 */
#endif

#ifdef NAVIERSTOKES_BULK
		      facbi = mass * All.NavierStokes_BulkViscosity / (rho * rho) * dwk_i / r;
		      facbj = P[j].Mass * All.NavierStokes_BulkViscosity /
			(SphP[j].d.Density * SphP[j].d.Density) * dwk_j / r;
#endif

#ifdef WINDS
		      if(P[j].Type == 0)
			if(SphP[j].DelayTime > 0)	/* No visc for wind particles */
			  {
			    faci = facj = 0;
#ifdef NAVIERSTOKES_BULK
			    facbi = facbj = 0;
#endif
			  }
#endif

#ifdef VISCOSITY_SATURATION
		      IonMeanFreePath_i = All.IonMeanFreePath * pow((Entropy * pow(rho * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1), 2.0) / rho;	/* u^2/rho */

		      IonMeanFreePath_j = All.IonMeanFreePath * pow((SphP[j].Entropy * pow(SphP[j].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1), 2.0) / SphP[j].d.Density;	/* u^2/rho */

		      for(k = 0, VelLengthScale_i = 0, VelLengthScale_j = 0; k < 3; k++)
			{
			  if(fabs(stressdiag[k]) > 0)
			    {
			      VelLengthScale_i = 2 * soundspeed_i / fabs(stressdiag[k]);

			      if(VelLengthScale_i < IonMeanFreePath_i && VelLengthScale_i > 0)
				{
				  stressdiag[k] = stressdiag[k] * (VelLengthScale_i / IonMeanFreePath_i);

				}
			    }
			  if(fabs(SphP[j].u.s.StressDiag[k]) > 0)
			    {
			      VelLengthScale_j = 2 * soundspeed_j / fabs(SphP[j].u.s.StressDiag[k]);

			      if(VelLengthScale_j < IonMeanFreePath_j && VelLengthScale_j > 0)
				{
				  SphP[j].u.s.StressDiag[k] = SphP[j].u.s.StressDiag[k] *
				    (VelLengthScale_j / IonMeanFreePath_j);

				}
			    }
			  if(fabs(stressoffdiag[k]) > 0)
			    {
			      VelLengthScale_i = 2 * soundspeed_i / fabs(stressoffdiag[k]);

			      if(VelLengthScale_i < IonMeanFreePath_i && VelLengthScale_i > 0)
				{
				  stressoffdiag[k] =
				    stressoffdiag[k] * (VelLengthScale_i / IonMeanFreePath_i);
				}
			    }
			  if(fabs(SphP[j].u.s.StressOffDiag[k]) > 0)
			    {
			      VelLengthScale_j = 2 * soundspeed_j / fabs(SphP[j].u.s.StressOffDiag[k]);

			      if(VelLengthScale_j < IonMeanFreePath_j && VelLengthScale_j > 0)
				{
				  SphP[j].u.s.StressOffDiag[k] = SphP[j].u.s.StressOffDiag[k] *
				    (VelLengthScale_j / IonMeanFreePath_j);
				}
			    }
			}
#endif

		      /* Acceleration due to the shear viscosity */
		      acc[0] += faci * (stressdiag[0] * dx + stressoffdiag[0] * dy + stressoffdiag[1] * dz)
			+ facj * (SphP[j].u.s.StressDiag[0] * dx + SphP[j].u.s.StressOffDiag[0] * dy +
				  SphP[j].u.s.StressOffDiag[1] * dz);

		      acc[1] += faci * (stressoffdiag[0] * dx + stressdiag[1] * dy + stressoffdiag[2] * dz)
			+ facj * (SphP[j].u.s.StressOffDiag[0] * dx + SphP[j].u.s.StressDiag[1] * dy +
				  SphP[j].u.s.StressOffDiag[2] * dz);

		      acc[2] += faci * (stressoffdiag[1] * dx + stressoffdiag[2] * dy + stressdiag[2] * dz)
			+ facj * (SphP[j].u.s.StressOffDiag[1] * dx + SphP[j].u.s.StressOffDiag[2] * dy +
				  SphP[j].u.s.StressDiag[2] * dz);

		      /*Acceleration due to the bulk viscosity */
#ifdef NAVIERSTOKES_BULK
#ifdef VISCOSITY_SATURATION
		      VelLengthScale_i = 0;
		      VelLengthScale_j = 0;

		      if(fabs(divvel) > 0)
			{
			  VelLengthScale_i = 3 * soundspeed_i / fabs(divvel);

			  if(VelLengthScale_i < IonMeanFreePath_i && VelLengthScale_i > 0)
			    {
			      divvel = divvel * (VelLengthScale_i / IonMeanFreePath_i);
			    }
			}

		      if(fabs(SphP[j].u.s.DivVel) > 0)
			{
			  VelLengthScale_j = 3 * soundspeed_j / fabs(SphP[j].u.s.DivVel);

			  if(VelLengthScale_j < IonMeanFreePath_j && VelLengthScale_j > 0)
			    {
			      SphP[j].u.s.DivVel = SphP[j].u.s.DivVel *
				(VelLengthScale_j / IonMeanFreePath_j);

			    }
			}
#endif


		      acc[0] += facbi * divvel * dx + facbj * SphP[j].u.s.DivVel * dx;
		      acc[1] += facbi * divvel * dy + facbj * SphP[j].u.s.DivVel * dy;
		      acc[2] += facbi * divvel * dz + facbj * SphP[j].u.s.DivVel * dz;
#endif
#endif /* end NAVIERSTOKES */


#ifdef MAGNETIC
#ifdef MAGNETIC_DISSIPATION
		      magfac_sym *= vsig * 0.5 * Balpha_ij;
		      dtEntropy += dTu_diss_b * 0.25 * vsig * mu0_1 * r2;
#ifdef MAGDISSIPATION_PERPEN
		      mft = (dBx * dx + dBy * dy + dBz * dz) / r;
		      mvt[0] = dBx - mft * dx / r;
		      mvt[1] = dBy - mft * dy / r;
		      mvt[2] = dBz - mft * dz / r;
		      mft = (mvt[0] * dx + mvt[1] * dy + mvt[2] * dy) / r;
		      dtB[0] += magfac_sym * mft * dx;
		      dtB[1] += magfac_sym * mft * dy;
		      dtB[2] += magfac_sym * mft * dz;
#else
		      magfac_sym *= (dBx * dx + dBy * dy + dBz * dz) / r;
		      dtB[0] += magfac_sym * dx;
		      dtB[1] += magfac_sym * dy;
		      dtB[2] += magfac_sym * dz;
#endif
#endif
#ifdef DIVBCLEANING_DEDNER
		      DtPhi -= All.DivBcleanHyperbolicSigma * phifac * 0.25 * vsig * vsig *
			(dBx * dx + dBy * dy + dBz * dz) / (fac_mu * fac_mu);
#endif
#endif


#ifdef CONDUCTION
		      rEnergy_j =
			SphP[j].SmoothedEntr * pow(SphP[j].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;
		      rEnergy_j_ns =
			SphP[j].Entropy * pow(SphP[j].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;
#ifdef CONDUCTION_CONSTANT
		      rKappa_j = All.ConductionCoeff;
#else
		      rKappa_j = All.ConductionCoeff * pow(rEnergy_j_ns, 2.5);
#ifdef CONDUCTION_SATURATION
		      electron_free_path =
			All.ElectronFreePathFactor * rEnergy_j * rEnergy_j / (SphP[j].d.Density * a3inv);
		      temp_scale_length =
			atime * fabs(smoothentr) / sqrt(SphP[j].GradEntr[0] * SphP[j].GradEntr[0] +
							SphP[j].GradEntr[1] * SphP[j].GradEntr[1] +
							SphP[j].GradEntr[2] * SphP[j].GradEntr[2]);
		      rKappa_j /= (1 + 4.2 * electron_free_path / temp_scale_length);
#endif
#endif
#ifdef SFR
		      if(SphP[j].d.Density * a3inv >= All.PhysDensThresh)
			rKappa_j = 0;
#endif
		      if((rKappa_i + rKappa_j) > 0)
			{
			  rEnergyFlow =
			    (mass * P[j].Mass) /
			    (rho * SphP[j].d.Density) *
			    4 * (rKappa_i * rKappa_j) /
			    (rKappa_i + rKappa_j) * (rEnergy_j - rEnergy_i) * 0.5 * (dwk_i +
										     dwk_j) / r * dtcond;
			}
		      else
			{
			  rEnergyFlow = 0;
			}
#ifdef WINDS
		      if(P[j].Type == 0)
			if(SphP[j].DelayTime > 0)	/* No cond. by wind particles */
			  {
			    rEnergyFlow = 0;
			  }
#endif
		      if(timestep > 0)
			SphP[j].CondEnergyChange += 0.5 * rEnergyFlow;
		      condEnergyChange += -0.5 * rEnergyFlow;
#endif /* end of CONDUCTION */


#ifdef CR_DIFFUSION
		      /* *************** COSMIC RAY DIFFUSION ******************* */ */
			/* In the following, the cosmic ray diffusion is computed.
			 * First, the total energy/mass transfer during particle i's
			 * timestep are evaluated, then the corresponding amount of
			 * specific CR energy and baryon fraction are adjusted in
			 * both particle i and particle j
			 */
			rKappaCR_j = All.CR_DiffusionCoeff;

		      if(All.CR_DiffusionDensScaling != 0.0)
			{
			  rKappaCR_j *= pow(SphP[j].d.Density * a3inv /
					    All.CR_DiffusionDensZero, All.CR_DiffusionDensScaling);
			}

		      if(All.CR_DiffusionEntropyScaling != 0.0)
			{
			  rKappaCR_j *= pow(SphP[j].Entropy /
					    All.CR_DiffusionEntropyZero, All.CR_DiffusionEntropyScaling);
			}

		      /* Compute the common term element that is typical for
		       * conduction/diffusion effects. So we only have to do
		       * some of the computations once. 
		       */

		      if(rKappaCR_i * rKappaCR_j != 0.0)
			{
			  rCR_DiffusionTerm = (All.CR_Alpha - 1) / (All.CR_Alpha - 1.33333) *
			    (mass * P[j].Mass) /
			    (rho * SphP[j].d.Density) *
			    4 * (rKappaCR_i * rKappaCR_j) /
			    (rKappaCR_i + rKappaCR_j) * 0.5 * (dwk_i + dwk_j) / r * dtdiff;

			  /* Compute the flow of total CR energy from particle j to particle i
			   * with smoothed SPH interpolants
			   * (conserved quantity) 
			   */
			  CR_q_j = SphP[j].CR_q0 * pow(SphP[j].d.Density * a3inv, 0.333333);

			  cr_efac_j =
			    CR_Tab_MeanEnergy(CR_q_j, All.CR_Alpha - 0.33333) / CR_Tab_MeanEnergy(CR_q_j,
												  All.
												  CR_Alpha);

			  rCR_EnergyFlow =
			    rCR_DiffusionTerm * (SphP[j].CR_SmoothE0 * SphP[j].d.Density *
						 pow(CR_q_j,
						     0.33333) * cr_efac_j - CR_E0_i * rho * pow(CR_q_i,
												0.33333) *
						 cr_efac_i);

			  /* Compute the CR particle flow 
			   */
			  rCR_MassFlow =
			    rCR_DiffusionTerm * (SphP[j].CR_Smoothn0 * SphP[j].d.Density *
						 pow(CR_q_j, 0.33333) - CR_n0_i * rho * pow(CR_q_i, 0.33333));


			  if(fabs(rCR_EnergyFlow) > 0.1 * DMAX(SphP[j].CR_SmoothE0, CR_E0_i))
			    {
			      rCR_EnergyFlow *=
				0.1 * DMAX(SphP[j].CR_SmoothE0, CR_E0_i) / fabs(rCR_EnergyFlow);
			    }

			  if(fabs(rCR_MassFlow) > 0.1 * DMAX(SphP[j].CR_Smoothn0, CR_n0_i))
			    {
			      rCR_MassFlow *= 0.1 * DMAX(SphP[j].CR_Smoothn0, CR_n0_i) / fabs(rCR_MassFlow);
			    }


			  CR_Particle_Inject(SphP + j,
					     0.5 * rCR_EnergyFlow / P[j].Mass,
					     0.5 * rCR_MassFlow / P[j].Mass);

			  rCR_EnergyChange_i += -0.5 * rCR_EnergyFlow / mass;
			  rCR_BaryonFractionChange_i += -0.5 * rCR_MassFlow / mass;

			}
#endif /* ************ COSMIC RAY DIFFUSION - END ***************** */
		    }
		}
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = HydroDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }


  /* Now collect the result at the right place */
  if(mode == 0)
    {
      for(k = 0; k < 3; k++)
	SphP[target].a.dHydroAccel[k] = acc[k];
      SphP[target].e.dDtEntropy = dtEntropy;
#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
      SphP[target].MinViscousDt = minViscousDt;
#else
      SphP[target].MaxSignalVel = maxSignalVel;
#endif
#ifdef CONDUCTION
      SphP[target].CondEnergyChange += condEnergyChange;
#ifdef OUTPUTCOOLRATE
      SphP[target].CondRate = 2 * atime * condEnergyChange / (dtcond * mass);
#endif
#endif
#ifdef MAGNETIC
      for(k = 0; k < 3; k++)
	SphP[target].DtB[k] = dtB[k];
#ifdef DIVBCLEANING_DEDNER
      SphP[target].DtPhi = DtPhi;
#endif
#ifdef TRACEDIVB
      SphP[target].divB = divB;
#endif
#endif

#if defined(CR_DIFFUSION)
      CR_Particle_Inject(SphP + target, rCR_EnergyChange_i, rCR_BaryonFractionChange_i);
#endif
    }
  else
    {
      for(k = 0; k < 3; k++)
	HydroDataResult[target].Acc[k] = acc[k];
      HydroDataResult[target].DtEntropy = dtEntropy;
#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
      HydroDataResult[target].MinViscousDt = minViscousDt;
#else
      HydroDataResult[target].MaxSignalVel = maxSignalVel;
#endif
#ifdef CONDUCTION
      HydroDataResult[target].CondEnergyChange = condEnergyChange;
#ifdef OUTPUTCOOLRATE
      HydroDataResult[target].CondRate = 2 * atime * condEnergyChange / (dtcond * mass);
#endif
#endif
#ifdef MAGNETIC
      for(k = 0; k < 3; k++)
	HydroDataResult[target].DtB[k] = dtB[k];
#ifdef DIVBCLEANING_DEDNER
      HydroDataResult[target].DtPhi = DtPhi;
#endif
#ifdef TRACEDIVB
      HydroDataResult[target].divB = divB;
#endif
#endif
#if defined(CR_DIFFUSION)
      HydroDataResult[target].CR_EnergyChange = rCR_EnergyChange_i;
      HydroDataResult[target].CR_BaryonFractionChange = rCR_BaryonFractionChange_i;
#endif
    }

  return 0;
}
