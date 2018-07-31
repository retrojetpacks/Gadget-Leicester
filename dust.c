#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

/*! \file dust.c
 *  \brief routines for the descriptions for dust
 */

/* "static" means that the structure is only used here and not passed outside
 * of this program */

#if defined(DUST)

static struct dustdata_in {
  MyDouble Pos[3];
  MyDouble Vel[3];
  MyFloat dustmass;
  MyFloat density;
  MyFloat dust_density;
  MyFloat Hsml;
  MyFloat DeltaDustMomentum[3];
  MyFloat DeltaDragEnergy;
#ifdef DUST_TWO_POPULATIONS
  MyFloat LogDustRadius_by_dt;
#endif
#ifdef DUST_REAL_PEBBLE_COLLISIONS
  MyFloat DustVel[3];
#endif
  MyIDType ID;
    int Index;
    int NodeList[NODELISTLENGTH];

}
*DustDataIn, *DustDataGet; /* These are things passed between processors */

static struct dustdata_out {
  MyLongDouble dustmass;
  MyLongDouble AccretedMomentum[3];
  MyLongDouble density;
  MyLongDouble dust_particle_density;
  MyLongDouble DeltaDustMomentum[3];
#ifdef DUST_REAL_PEBBLE_COLLISIONS
  MyLongDouble DustVel[3];
#endif
 MyLongDouble DeltaDragEnergy;
}
*DustDataResult, *DustDataOut; /* These are structures local to the
                                * processor */


static double hubble_a, ascale;

/* Calculate dust density */
void dust_density(void) {

   int i,j,k;
   int ndone_flag, ndone;
   int ngrp, recvTask, place, nexport, nimport, dummy;
   //   double a3inv;
   //   double rho, delta_vel, dust_vel, soundspeed, dt;
   MPI_Status status;

   if(ThisTask == 0) {
      printf("Beginning dust density estimation\n");
      fflush(stdout);
   }

   CPU_Step[CPU_MISC] += measure_time();

   if(All.ComovingIntegrationOn) {
      ascale = All.Time;
      hubble_a = hubble_function(All.Time);
      //      a3inv = 1 / (All.Time * All.Time * All.Time);
   } else {
      hubble_a = ascale = 1;
      //      a3inv = 1 / (All.Time * All.Time * All.Time);
   }


   /* allocate buffers to arrange communication */
   Ngblist = (int *) mymalloc(NumPart * sizeof(int));

   All.BunchSize =
   (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                            sizeof(struct dustdata_in) + sizeof(struct dustdata_out) +
                                    sizemax(sizeof(struct dustdata_in), sizeof(struct dustdata_out))));
   DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
   DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));

/* For Dust-Gas interactyion */
   i = FirstActiveParticle;  /* first particle for this task */ 


   do { 

      for (j = 0; j < NTask; j++) {   
        Send_count[j] = 0;
        Exportflag[j] = -1; 
      }   

      /* do local particles and prepare export list */
      for(nexport = 0; i >= 0; i = NextActiveParticle[i])
         if(P[i].Type == 2)
            if(dust_evaluate_density(i, 0, &nexport, Send_count) < 0) break;

      /* The sort here re-arranges the export Table */
#ifdef MYSORT
      mysort_dataindex(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#else
      qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#endif

      /* This passes the Export Table between the processes */
      MPI_Allgather(Send_count, NTask, MPI_INT, Sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);

      for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++) {   
         Recv_count[j] = Sendcount_matrix[j * NTask + ThisTask];
         nimport += Recv_count[j];

         if(j > 0) {
            Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1]; 
            Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1]; 
         }
      }   

      /* These are the things that need to be received == exported */
      DustDataGet = (struct dustdata_in *) mymalloc(nimport * sizeof(struct dustdata_in));
      /* Same for things imported */
      DustDataIn = (struct dustdata_in *) mymalloc(nexport * sizeof(struct dustdata_in));

      for(j = 0; j < nexport; j++) {   
         place = DataIndexTable[j].Index;

         for(k = 0; k < 3; k++) {
            DustDataIn[j].Pos[k] = P[place].Pos[k]; 
            DustDataIn[j].Vel[k] = P[place].Vel[k]; 
         }

         DustDataIn[j].dustmass = P[place].Mass;
         DustDataIn[j].dust_density = P[place].d7.dDUST_particle_density;
#ifdef DUST_REAL_PEBBLE_COLLISIONS
         DustDataIn[j].DustVel[0] = P[place].d9.dDUST_particle_velocity[0];
         DustDataIn[j].DustVel[1] = P[place].d9.dDUST_particle_velocity[1];
         DustDataIn[j].DustVel[2] = P[place].d9.dDUST_particle_velocity[2];
#endif
         DustDataIn[j].Hsml = PPP[place].Hsml;
         DustDataIn[j].ID = P[place].ID;

         memcpy(DustDataIn[j].NodeList,DataNodeList[DataIndexTable[j].IndexGet].NodeList,NODELISTLENGTH*sizeof(int));
      }

      /* This is a HORRIBLE MPI stuff */

      /* exchange particle data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++) {
         recvTask = ThisTask ^ ngrp;

         if(recvTask < NTask) {
            if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0) {

              /* get the particles */
               MPI_Sendrecv(&DustDataIn[Send_offset[recvTask]],
                            Send_count[recvTask] * sizeof(struct dustdata_in), MPI_BYTE,
                            recvTask, TAG_DENS_A,
                            &DustDataGet[Recv_offset[recvTask]],
                            Recv_count[recvTask] * sizeof(struct dustdata_in), MPI_BYTE,
                            recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
            }
         }
      }

      /* free up memory now that the data were passed around*/
      myfree(DustDataIn);

      DustDataResult = (struct dustdata_out *) mymalloc(nimport * sizeof(struct dustdata_out));
      DustDataOut = (struct dustdata_out *) mymalloc(nexport * sizeof(struct dustdata_out));

      /* now consider SPH particles that were sent to us. Thus mode=1 below */
      for(j = 0; j < nimport; j++) {
         place = DataIndexTable[j].Index;
         dust_evaluate_density(j, 1, &dummy, &dummy);
/*
         P[place].d7.dDUST_particle_density = 0.;
         PPP[place].d7.dDUST_particle_density = 0.;
   fprintf(diag_shc,"Before %d %d %g %g\n",place,j,P[place].d7.dDUST_particle_density,DustDataResult[j].dust_particle_density);
   fflush(diag_shc);
*/
      }

      if(i < 0)
         ndone_flag = 1;
      else
         ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      /* Now return the results into the processors where the exported blackholes reside */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++) {
         recvTask = ThisTask ^ ngrp;
         if(recvTask < NTask) {
            if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0) {

               /* send the results */
               MPI_Sendrecv(&DustDataResult[Recv_offset[recvTask]],
                            Recv_count[recvTask] * sizeof(struct dustdata_out),
                            MPI_BYTE, recvTask, TAG_DENS_B,
                            &DustDataOut[Send_offset[recvTask]],
                            Send_count[recvTask] * sizeof(struct dustdata_out),
                            MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, &status);
            }
         }
      }   

      /* add the result to the particles */
      for(j = 0; j < nexport; j++) {   
         place = DataIndexTable[j].Index;
         P[place].d7.dDUST_particle_density += DustDataOut[j].dust_particle_density;
#ifdef DUST_REAL_PEBBLE_COLLISIONS
	 for (k=0; k <3; k++) P[place].d9.dDUST_particle_velocity[k] += DustDataOut[j].DustVel[k];
#endif
/*
   fprintf(diag_shc,"Time=%g place=%d j=%d Den=%g Out=%g Result=%g\n",All.Time,place,j,P[place].d7.dDUST_particle_density,DustDataOut[j].dust_particle_density,DustDataResult[j].dust_particle_density);
   fflush(diag_shc);
*/
      }   

      myfree(DustDataOut);
      myfree(DustDataResult);
      myfree(DustDataGet);

   }

   while(ndone < NTask);

   for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
     {
       if(P[i].Type == 2) {
	 P[i].d7.DUST_particle_density = FLT(P[i].d7.dDUST_particle_density);
	 PPP[i].d7.dDUST_particle_density = P[i].d7.dDUST_particle_density;
	 PPP[i].d7.DUST_particle_density = FLT(P[i].d7.dDUST_particle_density);
	 for (k=0; k <3; k++) P[i].d9.DUST_particle_velocity[k] = FLT(P[i].d9.dDUST_particle_velocity[k]);
	 /*	 for (k=0; k <3; k++) PPP[i].d9.dDUST_particle_velocity[k] = P[i].d9.dDUST_particle_velocity[k];
		 for (k=0; k <3; k++) PPP[i].d9.DUST_particle_velocity[k] = FLT(P[i].d9.dDUST_particle_velocity[k]);*/
       }
     }

   myfree(DataNodeList);
   myfree(DataIndexTable);
   myfree(Ngblist);

}




/* Calculate dust-gas drag */
void dust_drag(void) {

   int n,i,j,k;
   int ndone_flag, ndone;
   int ngrp, sendTask, recvTask, place, nexport, nimport, dummy;
   double a3inv, latent_heat_const = 1.e11; // 1.e11 for rocks in erg/g //4.e10 for water ice
#ifdef DUST_FE_AND_ICE_GRAINS
   double latent_heat_const_ice = 4.e10 ;
#endif
   double rho, delta_vel, dust_vel, soundspeed, dt, acc, ts, vold, adot, dust_temperature, u_to_temp_fac, 
     lambda_h2, rey, C_drag, vnew_nog, f_pebble = 0.;
   double dE_latent = 0., pvap, delta_peb_vel = 0;
   MPI_Status status;
#ifdef DUST_GROWTH
   double adust_min = 0.1, adust_max =  1.e5, t_coll, v_coll;
#endif

#ifdef CONSTANT_MEAN_MOLECULAR_WEIGHT
  u_to_temp_fac = All.MeanWeight  * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
#else
  u_to_temp_fac = 2.3  * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
#endif


   if(ThisTask == 0) {
      printf("Beginning dust management\n");
      fflush(stdout);
   }

   CPU_Step[CPU_MISC] += measure_time();

   if(All.ComovingIntegrationOn) {
      ascale = All.Time;
      hubble_a = hubble_function(All.Time);
      a3inv = 1 / (All.Time * All.Time * All.Time);
   } else {
      hubble_a = ascale = 1;
      a3inv = 1 / (All.Time * All.Time * All.Time);
   }

/* For Dust only */
   for(n = FirstActiveParticle; n >= 0; n = NextActiveParticle[n]) {

      if(P[n].Type != 2) continue;

      for (k=0;k<=2;k++) {
         P[n].NewDragAcc[k] = 0.;
         P[n].DeltaDustMomentum[k] = 0.;
#ifdef DUST_REAL_PEBBLE_COLLISIONS
	 //         P[n].d9.dDUST_particle_velocity[k] = 0.;
#endif
      }
      P[n].DeltaDragEnergy = 0.;
#ifdef DUST_TWO_POPULATIONS
      P[n].LogDustRadius_by_dt = 0.;
#endif
#ifdef DUST_GROWTH_FIXED_SIZE
      PPP[n].DustRadius = All.InitialDustRadius;
      P[n].DustRadius = All.InitialDustRadius;
#endif
   }


   /* allocate buffers to arrange communication */
   Ngblist = (int *) mymalloc(NumPart * sizeof(int));

   All.BunchSize =
   (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                            sizeof(struct dustdata_in) + sizeof(struct dustdata_out) +
                                    sizemax(sizeof(struct dustdata_in), sizeof(struct dustdata_out))));
   DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
   DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));

/* For Dust-Gas interactyion */
   i = FirstActiveParticle;  /* first particle for this task */ 

   do { 

      for (j = 0; j < NTask; j++) {   
        Send_count[j] = 0;
        Exportflag[j] = -1; 
      }   

      /* do local particles and prepare export list */
      for(nexport = 0; i >= 0; i = NextActiveParticle[i])
      if(P[i].Type == 2) {

         /* Find drag acceleration */
         dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / hubble_a;
         rho = P[i].d1.DUST_Density;

#ifdef DUST_MDUST_GROW
	 if (All.Time <= 0.46) P[i].Mass *= exp( dt/All.VirtualStart);
#endif


         soundspeed = sqrt(8./M_PI * P[i].d2.DUST_Entropy * pow(rho, GAMMA_MINUS1));
         dust_vel = sqrt(pow(P[i].Vel[0], 2) + pow(P[i].Vel[1], 2) + pow(P[i].Vel[2], 2));
         delta_vel = sqrt(pow(P[i].Vel[0] - P[i].d3.DUST_SurroundingGasVel[0], 2) +
                          pow(P[i].Vel[1] - P[i].d3.DUST_SurroundingGasVel[1], 2) +
                          pow(P[i].Vel[2] - P[i].d3.DUST_SurroundingGasVel[2], 2));
	 
	 delta_peb_vel = 0.;
	 if (P[i].d7.DUST_particle_density > 0.)
	   {
	     P[i].d9.DUST_particle_velocity[0] = P[i].d9.DUST_particle_velocity[0]/P[i].d7.DUST_particle_density;
	     P[i].d9.DUST_particle_velocity[1] = P[i].d9.DUST_particle_velocity[1]/P[i].d7.DUST_particle_density;
	     P[i].d9.DUST_particle_velocity[2] = P[i].d9.DUST_particle_velocity[2]/P[i].d7.DUST_particle_density;
	     delta_peb_vel = sqrt(pow(P[i].Vel[0] - P[i].d9.DUST_particle_velocity[0], 2) + 
				  pow(P[i].Vel[1] - P[i].d9.DUST_particle_velocity[1], 2) +
				  pow(P[i].Vel[2] - P[i].d9.DUST_particle_velocity[2], 2));
	     P[i].DustVcoll = delta_peb_vel*All.UnitVelocity_in_cm_per_s/1.e2 + 1.e-30;
	   }

	 /* SN: grain material density. Assuming something between rock (3.5 g/cm^3) and iron (9 g/cm^3) */
	 double rho_dust = 3., rho_tidal = 0., rdist = 0., sigma_disc;
	 rdist = pow((P[i].Pos[0] - All.xbh)*(P[i].Pos[0] - All.xbh) + (P[i].Pos[1] - All.ybh)*(P[i].Pos[1] - All.ybh) + 
		     (P[i].Pos[2] - All.zbh)*(P[i].Pos[2] - All.zbh), 0.5);
	 sigma_disc = 20./pow(rdist, 1.5); /* disc surface desnity profile assumed */


	 /* SN -- DYNAMICS OF THE GRAIN  */	 
         if(dt > 0)
	   {
	     /*  Calculate H2 mean free path -- assuming cross section 10^{-15} cm^2 per molecule, Reynolds number
	      */
	     lambda_h2 = All.MeanWeight  * PROTONMASS/(All.UnitDensity_in_cgs *rho)/1.e-15/All.UnitLength_in_cm; /* in code units */
	     
	     rey = 6 * delta_vel * PPP[i].DustRadius/All.UnitLength_in_cm/(lambda_h2 *soundspeed);
	     
#ifndef DUST_EPSTEIN
	     if (3./2*lambda_h2*All.UnitLength_in_cm >= PPP[i].DustRadius) /* Epstein drag law */
	       {
		 ts = 1./(rho*soundspeed/(rho_dust*PPP[i].DustRadius) * All.UnitMass_in_g/All.UnitLength_in_cm/All.UnitLength_in_cm);
	       }
	     else
	       {
		 if (delta_vel > 0)
		   {
		     /* Expressions for C_drag from Weidenschilling 77 */
		     if (rey >= 800.) C_drag = 0.44;
		     if (rey < 800. && rey >= 1.) C_drag = 24.*pow(rey,-0.6);
		     if (rey < 1.) C_drag = 24./rey;
		     
		     ts = rho_dust*PPP[i].DustRadius/(rho*delta_vel)/All.UnitMass_in_g * All.UnitLength_in_cm * All.UnitLength_in_cm;
		     ts *= 8./3./C_drag;
		   }
		 else
		   {
		     ts = 0.66667/(rho*soundspeed/(rho_dust*PPP[i].DustRadius))/All.UnitMass_in_g*All.UnitLength_in_cm * PPP[i].DustRadius/lambda_h2;
		   }
	       }
#else /* Epstein drag only */
	     ts = 1./(rho*soundspeed/(rho_dust*PPP[i].DustRadius) * All.UnitMass_in_g/All.UnitLength_in_cm/All.UnitLength_in_cm);
#endif
	     double vhalf[3] ;
	     for(k = 0; k < 3; k++) {
               vold =  P[i].Vel[k];
	       /* vsteady is average velocity of gas and dust  */
	       double vsteady = (vold * P[i].d7.DUST_particle_density + P[i].d3.DUST_SurroundingGasVel[k]*rho)/(P[i].d7.DUST_particle_density + rho + 1.e-30);
	       //	       double vsteady = P[i].d3.DUST_SurroundingGasVel[k];
	       /* printf("DUST density = %g gas density = %g \n",P[i].d7.DUST_particle_density,rho ); */
	       /* fflush(stdout); */

	       P[i].Vel[k] = vsteady + (vold-vsteady)*exp(-dt/ts) + P[i].g.GravAccel[k]*ts * (1. - exp(-dt/ts));
	       //vhalf[k] = vsteady + (vold-vsteady)*exp(-dt/ts) + P[i].g.GravAccel[k]*ts/2. * (1. - exp(-dt*2./ts));
	       vnew_nog = vold + P[i].g.GravAccel[k] * dt;

	       /* Pass the momentum back to the gas to satisfy 2nd Newton's law. */
	       P[i].DeltaDustMomentum[k] -=  P[i].Mass*(vnew_nog - P[i].Vel[k]);
#ifndef DUST_NO_FRICTION_HEATING
	       P[i].DeltaDragEnergy +=  P[i].Mass *pow((P[i].Vel[k] - vsteady), 2)* (1. - exp(-2.*dt/ts))/2. ;
#endif
	     }
	     P[i].DustVcoll = pow(P[i].g.GravAccel[0]*P[i].g.GravAccel[0] + 
				  P[i].g.GravAccel[1]*P[i].g.GravAccel[1] + 
				  P[i].g.GravAccel[2]*P[i].g.GravAccel[2], 0.5)* 
	                          ts*All.UnitVelocity_in_cm_per_s/1.e2;
	       /*	     P[i].DustVcoll = pow((P[i].Vel[0] - vhalf[0])*(P[i].Vel[0] - vhalf[0]) + 
				  (P[i].Vel[1] - vhalf[1])*(P[i].Vel[1] - vhalf[1]) + 
				  (P[i].Vel[2] - vhalf[2])*(P[i].Vel[2] - vhalf[2]), 0.5) *All.UnitVelocity_in_cm_per_s/1.e2; */
	     P[i].DustVcoll += 0.2; /* SN -- added Brownian motion velocity of 50 cm/s */
	   }
	 /* SN -- END -- DYNAMICS OF THE GRAIN  */


	 /* SN -- Begin dust growth, vaporisation and resulting energy exchange with gas */

	 adot = 0.;
#if defined(DUST_GROWTH)
	 if (All.Time > 0 && All.Time > All.VirtualTime)
	   {
#ifdef DUST_TWO_POPULATIONS
	     adot = All.FeedBackVelocity*(0.01/All.InitialMicroDustZ)*
	       P[i].d8.DUST_micro_Density/(4.*rho_dust/All.UnitDensity_in_cgs)  * (delta_vel+1./All.UnitVelocity_in_cm_per_s);
	     /* Tapper off dust growth at too high collision velocity */
	     adot *= pow((1.e2/(1.e2+delta_vel*All.UnitVelocity_in_cm_per_s)),2);
#endif 
	     if (P[i].d7.DUST_particle_density > 0. && All.Time > 0)
	       {
		 t_coll = 4. * (rho_dust/All.UnitDensity_in_cgs)*(PPP[i].DustRadius/All.UnitLength_in_cm)
		   /P[i].d7.DUST_particle_density/(P[i].DustVcoll*1.e2/All.UnitVelocity_in_cm_per_s + 1.e-20);
		 P[i].LogDustRadius_by_dt = t_coll;
		 
#ifdef DUST_TWO_POPULATIONS
		 adot -= PPP[i].DustRadius/All.UnitLength_in_cm/t_coll*(0.01/All.InitialMicroDustZ);
#else

#ifndef DUST_REAL_PEBBLE_COLLISIONS
		 adot = PPP[i].DustRadius/All.UnitLength_in_cm/3./t_coll;
#else
		 /* SN: an extra factor of 2 comes from the fact that
		    self-collisions of grains are considered, so the
		    grain growth is "shared" between the two particles
		    -- otherwise we'd have to kill one of them */
		 adot = PPP[i].DustRadius/All.UnitLength_in_cm/3./t_coll;
		 if (All.FragmentationVelocity < 0.1)
		   {
		     adot = 0.;
		   }
		 else
		   {
		 /* SN: this turns growth into fragmentation at high speed collisions */
		     double xcoll = P[i].DustVcoll/All.FragmentationVelocity;
		     adot *= (1- pow(xcoll, 2))/(1+ pow(xcoll, 2));
		     /*		     if (ThisTask == 0 && i%100 == 0) printf("i = %d, t_coll = %g, V_coll = %g \n", 
				     i,  t_coll*All.UnitTime_in_s/3.15e7, P[i].DustVcoll);*/
		   }
#endif
#endif
	       }
	   }
#endif

	   double old_a = PPP[i].DustRadius;
#ifdef DUST_VAPORIZE
/* Dust vaporization */

	 if (rho > 0.*P[i].d7.DUST_particle_density) {
	   dust_temperature = M_PI/8.* pow(soundspeed * All.UnitVelocity_in_cm_per_s, 2)*All.MeanWeight  * PROTONMASS / BOLTZMANN;
	   /* SN -- the extra factor of M_PI/8. is due to how the 'sound speed' is defined -- it is (8/pi) * kb T/mu */
	   //c$$$  CHON from Helled et al 2008
	   //pvap = 5.53e7 * exp(-1.e4/dust_temperature);
#ifndef DUST_FE_AND_ICE_GRAINS
	   pvap = pow(10. , (-24605./dust_temperature+13.176)); //Rocks from Podolak et al 1988
#else
	 if (P[i].ID%2 == 0)
	   {
	     /* Fe grains */
	     pvap = pow(10. , (-24605./dust_temperature+13.176)); //Rocks from Podolak et al 1988
	   }
	 else
	   {
	     /* Water Ice grains */
	   if (dust_temperature <= 600.) // Water Ice from Podolak 88
	     {
	       pvap = pow(10.,(11.6 - 2104./dust_temperature));
	     }
	   else
	     {
	       pvap = 5. + 5.2e-3 * dust_temperature;
	     }
	   }
#endif

	   adot -= 1./(rho_dust*sqrt(2.*3.1415)) / soundspeed/All.UnitVelocity_in_cm_per_s * pvap/All.UnitVelocity_in_cm_per_s;
	 }
	 /*	 if (ThisTask == 0 && i%10 == 0) printf("i = %d, a = %g, T = %g da = %g \n",
		 i,  old_a, dust_temperature, adot*dt*All.UnitLength_in_cm);*/

	 double old_latent_heat = 0;
	 dE_latent = 0.;
#ifndef DUST_FE_AND_ICE_GRAINS
	 old_latent_heat = latent_heat_const * (old_a - adust_min)/(All.InitialDustRadius-adust_min);
#else
	 if (P[i].ID%2 == 0)
	   {
	     /* Fe grains */
	     old_latent_heat = latent_heat_const * (old_a - adust_min)/(All.InitialDustRadius-adust_min);
	   }
	 else
	   {
	     /* Water Ice grains */
	     old_latent_heat = latent_heat_const_ice * (old_a - adust_min)/(All.InitialDustRadius-adust_min);
	   }
#endif
#endif
#ifndef DUST_GROWTH_FIXED_SIZE
         PPP[i].DustRadius = PPP[i].DustRadius + (adot*dt)*All.UnitLength_in_cm;
#else
	 if (adot > 0.)
	   {
	     P[i].LogDustRadius_by_dt = old_a/(adot*All.UnitLength_in_cm+1.e-30);//adot*dt*All.UnitLength_in_cm/old_a;
	   }
	 else
	   {
	     P[i].LogDustRadius_by_dt = 1.e30;
	   }
	 /* if (ThisTask == 0 && i%1000 == 0) printf("i = %d, Pebble Dlog_by_dt = %g, rho_p/rho_gas = %g adot = %g \n", */
	 /* 					  i,  P[i].LogDustRadius_by_dt, P[i].d7.DUST_particle_density/rho, adot); */
#endif

#ifdef DUST_GROWTH /*RJH*/
	 if (PPP[i].DustRadius < adust_min) PPP[i].DustRadius = adust_min;
	 if (PPP[i].DustRadius > adust_max) PPP[i].DustRadius = adust_max;
#endif
#ifdef DUST_TWO_POPULATIONS
	 if (dt > 0.) 
	   {
	     P[i].LogDustRadius_by_dt = adot*dt*All.UnitLength_in_cm/old_a;
	     if (ThisTask == 0 && i%1000000 == 0) printf("i = %d, Pebble Dlog_by_dt = %g, rho_p/rho_gas = %g \n", 
							 i,  P[i].LogDustRadius_by_dt/dt, P[i].d7.DUST_particle_density/rho);
	   }
	 else
	   {
	     P[i].LogDustRadius_by_dt = -1.e-30;
	   }
	 double delta_mass =  P[i].LogDustRadius_by_dt, dm_max = 0.1;
	 if (abs(delta_mass) > dm_max) 
	   {
	     P[i].LogDustRadius_by_dt = dm_max * P[i].LogDustRadius_by_dt/abs(delta_mass);
	   }
	 P[i].Mass += P[i].Mass * 3.* P[i].LogDustRadius_by_dt;
	 if (P[i].Mass < 1.e-15) P[i].Mass = 1.e-15;
#endif

         P[i].DustRadius = PPP[i].DustRadius;
#ifdef DUST_VAPORIZE
	 double new_latent_heat = 0.;
#ifndef DUST_FE_AND_ICE_GRAINS
	 /* same composition grains */
	 new_latent_heat = latent_heat_const * ( P[i].DustRadius - adust_min)/(All.InitialDustRadius-adust_min); 
#else
	 if (P[i].ID%2 == 0)
	   {
	     /* Fe grains */
	     new_latent_heat = latent_heat_const * ( P[i].DustRadius - adust_min)/(All.InitialDustRadius-adust_min); 
	   }
	 else
	   {
	     /* Water Ice grains */
	     new_latent_heat = latent_heat_const_ice * ( P[i].DustRadius - adust_min)/(All.InitialDustRadius-adust_min); 
	   }
#endif
	 dE_latent = new_latent_heat - old_latent_heat;
	 P[i].DeltaDragEnergy += dE_latent * P[i].Mass* All.UnitMass_in_g/All.UnitEnergy_in_cgs; //this takes heat away from gas if grains are vaporised
#endif

	 

         if(dust_evaluate_select(i, 0, &nexport, Send_count) < 0)
            break;
      }

      /* The sort here re-arranges the export Table */
#ifdef MYSORT
      mysort_dataindex(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#else
      qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#endif

      /* This passes the Export Table between the processes */
      MPI_Allgather(Send_count, NTask, MPI_INT, Sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);

      for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++) {   
         Recv_count[j] = Sendcount_matrix[j * NTask + ThisTask];
         nimport += Recv_count[j];

         if(j > 0) {
            Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1]; 
            Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1]; 
         }
      }   

      /* These are the things that need to be received == exported */
      DustDataGet = (struct dustdata_in *) mymalloc(nimport * sizeof(struct dustdata_in));
      /* Same for things imported */
      DustDataIn = (struct dustdata_in *) mymalloc(nexport * sizeof(struct dustdata_in));

      for(j = 0; j < nexport; j++) {   
         place = DataIndexTable[j].Index;

         for(k = 0; k < 3; k++) {
            DustDataIn[j].Pos[k] = P[place].Pos[k]; 
            DustDataIn[j].Vel[k] = P[place].Vel[k]; 
            DustDataIn[j].DeltaDustMomentum[k] = P[place].DeltaDustMomentum[k];
#ifdef DUST_REAL_PEBBLE_COLLISIONS
	    DustDataIn[j].DustVel[k] = P[place].d9.dDUST_particle_velocity[k];
#endif
         }

         DustDataIn[j].dustmass = P[place].Mass;
         DustDataIn[j].density = P[place].d1.DUST_Density;
         DustDataIn[j].Hsml = P[place].Hsml;
         DustDataIn[j].DeltaDragEnergy = P[place].DeltaDragEnergy;
#ifdef DUST_TWO_POPULATIONS
         DustDataIn[j].LogDustRadius_by_dt = P[place].LogDustRadius_by_dt;
#endif
         DustDataIn[j].ID = P[place].ID;

         memcpy(DustDataIn[j].NodeList,DataNodeList[DataIndexTable[j].IndexGet].NodeList,NODELISTLENGTH*sizeof(int));
      }

      /* This is a HORRIBLE MPI stuff */

      /* exchange particle data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++) {
         sendTask = ThisTask;
         recvTask = ThisTask ^ ngrp;

         if(recvTask < NTask) {
            if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0) {

              /* get the particles */
               MPI_Sendrecv(&DustDataIn[Send_offset[recvTask]],
                            Send_count[recvTask] * sizeof(struct dustdata_in), MPI_BYTE,
                            recvTask, TAG_DENS_A,
                            &DustDataGet[Recv_offset[recvTask]],
                            Recv_count[recvTask] * sizeof(struct dustdata_in), MPI_BYTE,
                            recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
            }
         }
      }

      /* free up memory now that the data were passed around*/
      myfree(DustDataIn);

      DustDataResult = (struct dustdata_out *) mymalloc(nimport * sizeof(struct dustdata_out));
      DustDataOut = (struct dustdata_out *) mymalloc(nexport * sizeof(struct dustdata_out));

      /* now consider SPH particles that were sent to us. Thus mode=1 below */
      for(j = 0; j < nimport; j++) {
         place = DataIndexTable[j].Index;
         dust_evaluate_select(j, 1, &dummy, &dummy);
      }

      if(i < 0)
         ndone_flag = 1;
      else
         ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      /* Now return the results into the processors where the exported blackholes reside */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++) {
         sendTask = ThisTask;
         recvTask = ThisTask ^ ngrp;
         if(recvTask < NTask) {
            if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0) {

               /* send the results */
               MPI_Sendrecv(&DustDataResult[Recv_offset[recvTask]],
                            Recv_count[recvTask] * sizeof(struct dustdata_out),
                            MPI_BYTE, recvTask, TAG_DENS_B,
                            &DustDataOut[Send_offset[recvTask]],
                            Send_count[recvTask] * sizeof(struct dustdata_out),
                            MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, &status);
            }
         }
      }   

      /* add the result to the particles */
      for(j = 0; j < nexport; j++) {   
         place = DataIndexTable[j].Index;
         for (k=0;k<=2;k++) P[place].DeltaDustMomentum[k] += DustDataOut[j].DeltaDustMomentum[k];
         P[place].DeltaDragEnergy += DustDataOut[j].DeltaDragEnergy;
      }   

      myfree(DustDataOut);
      myfree(DustDataResult);
      myfree(DustDataGet);

   }

   while(ndone < NTask);

   myfree(DataNodeList);
   myfree(DataIndexTable);
   myfree(Ngblist);

}


/* SN: finds dust particle neighbours of a dust particle, using HSML
   found from iterations on GAS PARTICLE neighbours in density.c */
int dust_evaluate_density(int target, int mode, int *nexport, int *nSend_local) {

  MyFloat *pos, h_i, *DeltaDustMomentum,dustmass, dust_density,density,dust_particle_density;
   double dx, dy, dz, dt, h_i2, r2, u;
   double wk, hinv, hinv3, a3inv, time_hubble_a;
   int index, id, n, j, k;
   int startnode, numngb, listindex = 0;
   MyFloat return_energy;
   double acc_boundary;
#ifdef DUST_REAL_PEBBLE_COLLISIONS
   MyFloat dust_vel[3];
#endif


   if(All.ComovingIntegrationOn) {
      /* Factors for comoving integration of hydro */
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = hubble_function(All.Time);
      time_hubble_a = All.Time * hubble_a;
      ascale = All.Time;
   } else
      a3inv = ascale = time_hubble_a = hubble_a = 1.;

   if(mode == 0) {
      pos = P[target].Pos;
      h_i = PPP[target].Hsml;
      dustmass = PPP[target].Mass;
      density = PPP[target].d1.DUST_Density;
      //dust_density = P[target].d7.dDUST_particle_density;
      /*      DeltaDustMomentum = P[target].DeltaDustMomentum;
	      return_energy = P[target].DeltaDragEnergy;*/
      index = target;
      id = P[target].ID;
   } else {
      pos = DustDataGet[target].Pos;
      h_i = DustDataGet[target].Hsml;
      dustmass = DustDataGet[target].dustmass;
      density = DustDataGet[target].density;
      //dust_density = DustDataGet[target].dust_density;
      /*      DeltaDustMomentum = DustDataGet[target].DeltaDustMomentum;
	      return_energy = DustDataGet[target].DeltaDragEnergy;*/
      index = DustDataGet[target].Index;
      id = DustDataGet[target].ID;
   }

   dust_particle_density = 0.;
#ifdef DUST_REAL_PEBBLE_COLLISIONS
   for (k=0; k <3; k++) dust_vel[k] = 0.;
#endif

   /* Now start the actual SPH reaction due to the dust action */
   if(mode == 0) {
      startnode = All.MaxPart;      /* root node */
   } else {
      startnode = DustDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;      /* open it */
   }

   while(startnode >= 0) {
      while(startnode >= 0) {
         numngb = ngb_treefind_dust_active(pos, h_i, target, &startnode, mode, nexport, nSend_local, 2);

         if(numngb < 0) return -1;

         for(n = 0; n < numngb; n++) {
            j = Ngblist[n];

/*
            fprintf (diag_shc,"%d %d %d %d %g\n",target,P[target].Type,j,P[j].Type,P[j].Mass);
            fflush(diag_shc);
*/

            if(P[j].Mass > 0) {
               dx = pos[0] - P[j].Pos[0];
               dy = pos[1] - P[j].Pos[1];
               dz = pos[2] - P[j].Pos[2];
#ifdef PERIODIC /*  now find the closest image in the given box size  */
               if(dx > boxHalf_X) dx -= boxSize_X;
               if(dx < -boxHalf_X) dx += boxSize_X;
               if(dy > boxHalf_Y) dy -= boxSize_Y;
               if(dy < -boxHalf_Y) dy += boxSize_Y;
               if(dz > boxHalf_Z) dz -= boxSize_Z;
               if(dz < -boxHalf_Z) dz += boxSize_Z;
#endif
               r2 = dx * dx + dy * dy + dz * dz;

//             u = sqrt(r2)/P[j].Hsml;
               u = sqrt(r2)/h_i;

               if (u < 1) {
//                hinv = 1 /P[j].Hsml;
                  hinv = 1 /h_i;
#ifndef  TWODIMS
                  hinv3 = hinv * hinv * hinv;
#else
                  hinv3 = hinv * hinv / boxSize_Z;
#endif
                  if(u < 0.5) wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                  else wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

                  /* density estimation treatment */
                  dust_particle_density += FLT(dustmass * wk);
#ifdef DUST_REAL_PEBBLE_COLLISIONS
		  for (k=0; k <3; k++) dust_vel[k] += FLT(dustmass * wk * P[j].Vel[k]);
#endif

               }

            }
         }
      }

      if(mode == 1) {
         listindex++;
         if(listindex < NODELISTLENGTH) {
            startnode = DustDataGet[target].NodeList[listindex];
            if(startnode >= 0) startnode = Nodes[startnode].u.d.nextnode;      /* open it */
         }
      }
   }

   if(mode == 0){
      P[target].d7.dDUST_particle_density = dust_particle_density;
#ifdef DUST_REAL_PEBBLE_COLLISIONS
      P[target].d9.dDUST_particle_velocity[0] = dust_vel[0];
      P[target].d9.dDUST_particle_velocity[1] = dust_vel[1];
      P[target].d9.dDUST_particle_velocity[2] = dust_vel[2];
#endif
   } else {
      DustDataResult[target].dust_particle_density = dust_particle_density;
#ifdef DUST_REAL_PEBBLE_COLLISIONS
      for (k=0; k <3; k++) DustDataResult[target].DustVel[k] = dust_vel[k];
#endif
   }


   return 0;
}



int dust_evaluate_select(int target, int mode, int *nexport, int *nSend_local) {

  MyFloat *pos, h_i, *DeltaDustMomentum,dustmass,density, u_old, u_new, u_inc, entropy;
  double dx, dy, dz, dt, h_i2, r2, u, dmax1, dmax2;
   double wk, hinv, hinv3, a3inv, time_hubble_a;
   int index, id, n, j, k;
   int startnode, numngb, listindex = 0;
   MyFloat return_energy;
   double acc_boundary;
#ifdef DUST_TWO_POPULATIONS
   MyFloat d_micro_dust_mass;
#endif

   if(All.ComovingIntegrationOn) {
      /* Factors for comoving integration of hydro */
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = hubble_function(All.Time);
      time_hubble_a = All.Time * hubble_a;
      ascale = All.Time;
   } else
      a3inv = ascale = time_hubble_a = hubble_a = 1.;

   if(mode == 0) {
      pos = P[target].Pos;
      h_i = PPP[target].Hsml;
      dustmass = PPP[target].Mass;
      density = PPP[target].d1.DUST_Density;
      DeltaDustMomentum = P[target].DeltaDustMomentum;
      return_energy = P[target].DeltaDragEnergy;
#ifdef DUST_TWO_POPULATIONS
      d_micro_dust_mass =  dustmass * 3.* P[target].LogDustRadius_by_dt;
#endif
      index = target;
      id = P[target].ID;
   } else {
      pos = DustDataGet[target].Pos;
      h_i = DustDataGet[target].Hsml;
      dustmass = DustDataGet[target].dustmass;
      density = DustDataGet[target].density;
      DeltaDustMomentum = DustDataGet[target].DeltaDustMomentum;
      return_energy = DustDataGet[target].DeltaDragEnergy;
#ifdef DUST_TWO_POPULATIONS
      d_micro_dust_mass =  3.*dustmass * DustDataGet[target].LogDustRadius_by_dt;
#endif
      index = DustDataGet[target].Index;
      id = DustDataGet[target].ID;
   }

   /* initialize variables before SPH loop is started */
   h_i2 = h_i * h_i;

   /* Now start the actual SPH reaction due to the dust action */
   if(mode == 0) {
      startnode = All.MaxPart;      /* root node */
   } else {
      startnode = DustDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;      /* open it */
   }

   while(startnode >= 0) {
      while(startnode >= 0) {
         numngb = ngb_treefind_dust_active(pos, h_i, target, &startnode, mode, nexport, nSend_local, 0);

         if(numngb < 0) return -1;

         for(n = 0; n < numngb; n++) {
            j = Ngblist[n];

            if(P[j].Mass > 0) {
               dx = pos[0] - P[j].Pos[0];
               dy = pos[1] - P[j].Pos[1];
               dz = pos[2] - P[j].Pos[2];
#ifdef PERIODIC /*  now find the closest image in the given box size  */
               if(dx > boxHalf_X) dx -= boxSize_X;
               if(dx < -boxHalf_X) dx += boxSize_X;
               if(dy > boxHalf_Y) dy -= boxSize_Y;
               if(dy < -boxHalf_Y) dy += boxSize_Y;
               if(dz > boxHalf_Z) dz -= boxSize_Z;
               if(dz < -boxHalf_Z) dz += boxSize_Z;
#endif
               r2 = dx * dx + dy * dy + dz * dz;

               u = sqrt(r2)/h_i;

               if (u < 1) {
                  hinv = 1 /h_i;
#ifndef  TWODIMS
                  hinv3 = hinv * hinv * hinv;
#else
                  hinv3 = hinv * hinv / boxSize_Z;
#endif
                  if(u < 0.5) wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                  else wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

                  if(All.ComovingIntegrationOn)
                     dt = (P[j].TimeBin ? (1 << P[j].TimeBin) : 0) * All.Timebase_interval / time_hubble_a;
                  else
                     dt = (P[j].TimeBin ? (1 << P[j].TimeBin) : 0) * All.Timebase_interval;


		  /* SN. Note that the energy and momentum are passed
		     from dust to gas explicitly when dust particles
		     are active, not SPH. There is thus a certain
		     mismatch between time stepping of dust and gas
		     particles but conservation is more important to
		     me. */
                  if (dt>0. && density>0.)
		    {
		      //		      for (k=0;k<=2;k++) SphP[j].da.DragAccel[k] -= DeltaDustMomentum[k]/density*FLT(wk) / dt;
		      for (k=0;k<=2;k++) P[j].Vel[k] -= DeltaDustMomentum[k]/density*FLT(wk);

		       u_old = DMAX(All.MinEgySpec, SphP[j].Entropy/ GAMMA_MINUS1 * pow(density, GAMMA_MINUS1));
		       u_new = u_old + FLT(return_energy*wk/density);
		       u_inc = u_new/u_old;
		       if (u_inc > 1.5) u_inc = 1.5; /* don't allow to much heating in one interaction -- should never happen */
		       SphP[j].Entropy *= u_inc;

		       SphP[j].dh.DragHeating += 1.e-20*FLT(return_energy*wk/density*P[j].Mass/dt); 
#ifdef DUST_TWO_POPULATIONS
		       P[j].d_MicroDustMass -= FLT(d_micro_dust_mass *wk/density*P[j].Mass);
#endif
		       //if (P[j].MicroDustMass <= 1.e-12*P[j].Mass) P[j].MicroDustMass = 1.e-12*P[j].Mass;
		    }

               }
            }
         }
      }

      if(mode == 1) {
         listindex++;
         if(listindex < NODELISTLENGTH) {
            startnode = DustDataGet[target].NodeList[listindex];
            if(startnode >= 0) startnode = Nodes[startnode].u.d.nextnode;      /* open it */
         }
      }
   }

   return 0;
}


int dust_evaluate(int target, int mode, int *nexport, int *nSend_local) {

   MyFloat *pos, h_i, *DeltaDustMomentum,dustmass,density;
   double dx, dy, dz, dt, h_i2, r2, u;
   double wk, hinv, hinv3, a3inv, time_hubble_a;
   int index, id, n, j, k;
   int startnode, numngb, listindex = 0;
   MyFloat return_energy;
   double acc_boundary;

   if(All.ComovingIntegrationOn) {
      /* Factors for comoving integration of hydro */
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = hubble_function(All.Time);
      time_hubble_a = All.Time * hubble_a;
      ascale = All.Time;
   } else
      a3inv = ascale = time_hubble_a = hubble_a = 1.;

   if(mode == 0) {
      pos = P[target].Pos;
      h_i = PPP[target].Hsml;
      dustmass = PPP[target].Mass;
      density = PPP[target].d1.DUST_Density;
      DeltaDustMomentum = P[target].DeltaDustMomentum;
      return_energy = P[target].DeltaDragEnergy;
      index = target;
      id = P[target].ID;
   } else {
      pos = DustDataGet[target].Pos;
      h_i = DustDataGet[target].Hsml;
      dustmass = DustDataGet[target].dustmass;
      density = DustDataGet[target].density;
      DeltaDustMomentum = DustDataGet[target].DeltaDustMomentum;
      return_energy = DustDataGet[target].DeltaDragEnergy;
      index = DustDataGet[target].Index;
      id = DustDataGet[target].ID;
   }

/*  Delete these later
   *DeltaDustMomentum = 0.;
   return_energy = 0.;
*/

   /* initialize variables before SPH loop is started */
   h_i2 = h_i * h_i;

   /* Now start the actual SPH reaction due to the dust action */
   if(mode == 0) {
      startnode = All.MaxPart;      /* root node */
   } else {
      startnode = DustDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;      /* open it */
   }

   while(startnode >= 0) {
      while(startnode >= 0) {
         numngb = ngb_treefind_dust_active(pos, h_i, target, &startnode, mode, nexport, nSend_local, 2);

         if(numngb < 0) return -1;

         for(n = 0; n < numngb; n++) {
            j = Ngblist[n];

            if(P[j].Type == 2 && P[j].Mass > 0) {
               dx = pos[0] - P[j].Pos[0];
               dy = pos[1] - P[j].Pos[1];
               dz = pos[2] - P[j].Pos[2];
#ifdef PERIODIC /*  now find the closest image in the given box size  */
               if(dx > boxHalf_X) dx -= boxSize_X;
               if(dx < -boxHalf_X) dx += boxSize_X;
               if(dy > boxHalf_Y) dy -= boxSize_Y;
               if(dy < -boxHalf_Y) dy += boxSize_Y;
               if(dz > boxHalf_Z) dz -= boxSize_Z;
               if(dz < -boxHalf_Z) dz += boxSize_Z;
#endif
               r2 = dx * dx + dy * dy + dz * dz;

               acc_boundary = All.InnerBoundary;

               if(pow(r2,0.5)<acc_boundary && P[j].SwallowID<id) P[j].SwallowID = id;

            }
         }
      }

      if(mode == 1) {
         listindex++;
         if(listindex < NODELISTLENGTH) {
            startnode = DustDataGet[target].NodeList[listindex];
            if(startnode >= 0) startnode = Nodes[startnode].u.d.nextnode;      /* open it */
         }
      }
   }

   return 0;
}

int dust_evaluate_swallow(int target, int mode, int *nexport, int *nSend_local) {

   MyFloat *pos, h_i, *DeltaDustMomentum,dustmass,density;
   double dx, dy, dz, dt, h_i2, r2, u;
   double wk, hinv, hinv3, a3inv, time_hubble_a;
   int index, id, n, j, k;
   int startnode, numngb, listindex = 0;
   MyFloat accreted_mass, accreted_DUST_mass, accreted_momentum[3];
   MyFloat return_energy;

   if(All.ComovingIntegrationOn) {
      /* Factors for comoving integration of hydro */
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = hubble_function(All.Time);
      time_hubble_a = All.Time * hubble_a;
      ascale = All.Time;
   } else
      a3inv = ascale = time_hubble_a = hubble_a = 1.;

   if(mode == 0) {
      pos = P[target].Pos;
      h_i = PPP[target].Hsml;
      dustmass = PPP[target].Mass;
      density = PPP[target].d1.DUST_Density;
      DeltaDustMomentum = P[target].DeltaDustMomentum;
      return_energy = P[target].DeltaDragEnergy;
      index = target;
      id = P[target].ID;
   } else {
      pos = DustDataGet[target].Pos;
      h_i = DustDataGet[target].Hsml;
      dustmass = DustDataGet[target].dustmass;
      density = DustDataGet[target].density;
      DeltaDustMomentum = DustDataGet[target].DeltaDustMomentum;
      return_energy = DustDataGet[target].DeltaDragEnergy;
      index = DustDataGet[target].Index;
      id = DustDataGet[target].ID;
   }

   accreted_mass = 0;
   accreted_DUST_mass = 0;
   accreted_momentum[0] = accreted_momentum[1] = accreted_momentum[2] = 0;

   /* initialize variables before SPH loop is started */
   h_i2 = h_i * h_i;

   /* Now start the actual SPH computation for this particle */
   if(mode == 0) {
      startnode = All.MaxPart;      /* root node */
   } else {
      startnode = DustDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;      /* open it */
   }

   while(startnode >= 0) {
      while(startnode >= 0) {
         numngb = ngb_treefind_dust_active(pos, h_i, target, &startnode, mode, nexport, nSend_local, 2);

         if(numngb < 0) return -1;

         for(n = 0; n < numngb; n++) {
            j = Ngblist[n];

/*          Dust-Dust merging */
            if(P[j].Type == 2 && P[j].SwallowID == id) {

               accreted_mass += FLT(P[j].Mass);
               for(k = 0; k < 3; k++) accreted_momentum[k] += FLT(P[j].Mass * P[j].Vel[k]);
               P[j].Mass = 0;

               printf("DUST-DUST merging, %d and %d, Mass is %g\n",target,j,accreted_mass);
               fflush(stdout);

            }
         }
      }

      if(mode == 1) {
         listindex++;
         if(listindex < NODELISTLENGTH) {
            startnode = DustDataGet[target].NodeList[listindex];
            if(startnode >= 0) startnode = Nodes[startnode].u.d.nextnode;      /* open it */
         }
      }
   }

  /* Now collect the result at the right place */
  if(mode == 0)
    {    
      P[target].d4.dDUST_accreted_Mass = accreted_mass;
//      P[target].d5.dDUST_accreted_DUSTMass = accreted_DUST_mass;
      for(k = 0; k < 3; k++) 
      P[target].d6.dDUST_accreted_momentum[k] = accreted_momentum[k];
    }    
  else 
    {    
      DustDataResult[target].dustmass = accreted_mass;
//      DustDataResult[target].DUST_Mass = accreted_DUST_mass;
      for(k = 0; k < 3; k++) 
      DustDataResult[target].AccretedMomentum[k] = accreted_momentum[k];
    }    

   return 0;
}


int ngb_treefind_dust_active(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode, int *nexport, int *nsend_local, int ptype_want) {

   int numngb, no, p, task, nexport_save;
   struct NODE *current;
   MyDouble dx, dy, dz, dist;

#ifdef PERIODIC
   MyDouble xtmp;
#endif
   nexport_save = *nexport;

  numngb = 0;
  no = *startnode;

   while(no >= 0) {
      if(no < All.MaxPart) {   /* single particle */
         p = no;
         no = Nextnode[no];

         if(P[p].Type != ptype_want) continue; /* only Sph particles here */

         dist = hsml;
         dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]);
         if(dx > dist) continue;
         dy = NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]);
         if(dy > dist) continue;
         dz = NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]);
         if(dz > dist) continue;
         if(dx * dx + dy * dy + dz * dz > dist * dist) continue;

         Ngblist[numngb++] = p;

      } else {

         if(no >= All.MaxPart + MaxNodes) { /* pseudo particle */
            if(mode == 1)
            endrun(12312);

            if(Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target) {
              Exportflag[task] = target;
              Exportnodecount[task] = NODELISTLENGTH;
            }

            if(Exportnodecount[task] == NODELISTLENGTH) {
              if(*nexport >= All.BunchSize) {
                 *nexport = nexport_save;
                 for(task = 0; task < NTask; task++)
                 nsend_local[task] = 0;
                 for(no = 0; no < nexport_save; no++)
                 nsend_local[DataIndexTable[no].Task]++;
                 return -1;
              }
              Exportnodecount[task] = 0;
              Exportindex[task] = *nexport;
              DataIndexTable[*nexport].Task = task;
              DataIndexTable[*nexport].Index = target;
              DataIndexTable[*nexport].IndexGet = *nexport;
              *nexport = *nexport + 1;
              nsend_local[task]++;
            }

            DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] = DomainNodeIndex[no - (All.MaxPart + MaxNodes)];

            if(Exportnodecount[task] < NODELISTLENGTH) DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;

            no = Nextnode[no - MaxNodes];
            continue;
         }

         current = &Nodes[no];

         if(mode == 1) {
            if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL)) {
               *startnode = -1;
               return numngb;
            }
         }

         no = current->u.d.sibling;  /* in case the node can be discarded */

         dist = hsml + 0.5 * current->len;;
         dx = NGB_PERIODIC_LONG_X(current->center[0] - searchcenter[0]);
         if(dx > dist) continue;
         dy = NGB_PERIODIC_LONG_Y(current->center[1] - searchcenter[1]);
         if(dy > dist) continue;
         dz = NGB_PERIODIC_LONG_Z(current->center[2] - searchcenter[2]);
         if(dz > dist) continue;
         /* now test against the minimal sphere enclosing everything */
         dist += FACT1 * current->len;
         if(dx * dx + dy * dy + dz * dz > dist * dist) continue;

         no = current->u.d.nextnode; /* ok, we need to open the node */
      }
   }

   *startnode = -1;
   return numngb;
}

#if defined(DUST_PEBBLES_BORN) /* block to create pebble/dust particles */
void pebbles_born(void)
{
  int i, j, k, stars_spawned = 0, tot_spawned = 0;
  double dt, rand_peb, probability_of_birth, peb_init_mass, dir[3], grav;

#ifdef DUST_TWO_POPULATIONS
  peb_init_mass = 0.01 * All.OriginalGasMass *All.PebbleBirthMass; 
#else
  //  double total_mgas=0., local_mgas=0., total_mgas_give_birth, total_mgas_all;
  peb_init_mass = All.OriginalGasMass *All.PebbleInjectionMass; 
#endif

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
  {
    dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;
#if defined(DUST_TWO_POPULATIONS)
    if(P[i].Type == 0 && P[i].MicroDustMass > 0.)
      {
	rand_peb = get_random_number(P[i].ID + 300*All.NumCurrentTiStep);
	probability_of_birth = P[i].MicroDustMass/All.PebbleBirthTime *dt/peb_init_mass;
	P[i].MicroDustMass /= (1. + dt/All.PebbleBirthTime);
	if (rand_peb < probability_of_birth)
	  {
	    P[NumPart + stars_spawned] = P[i];
	    P[NumPart + stars_spawned].Type = 2;
	    for(j = 0; j < 3; j++)
	      {
		/* add a small random shift to the initial pebble particle position */
		P[NumPart + stars_spawned].Pos[j] += 0.02 *(-2 + gsl_rng_uniform(random_generator))
		  * P[i].Hsml;
	      }

	    
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
	    P[NumPart + stars_spawned].Mass = peb_init_mass;
	    P[NumPart + stars_spawned].DustRadius = All.InitialDustRadius;
	    PPP[NumPart + stars_spawned].DustRadius = P[NumPart + stars_spawned].DustRadius;
	    force_add_star_to_tree(i, NumPart + stars_spawned);
	    stars_spawned++;

	  }
      }
#else
    if(P[i].Type == 0)
      {
//	total_mgas += P[i].Mass;
	rand_peb = get_random_number(P[i].ID + 300*All.NumCurrentTiStep);
	probability_of_birth = dt/All.DustInjectionTime;
	/*	printf("\n  PebbleInjectionDens %g  PebbleInjectionSlab %g .... \n", 
		All.PebbleInjectionDens, All.PebbleInjectionSlab);*/
	grav = sqrt(P[i].g.GravAccel[0] * P[i].g.GravAccel[0] +
		  P[i].g.GravAccel[1] * P[i].g.GravAccel[1] + P[i].g.GravAccel[2] * P[i].g.GravAccel[2]);
	for (k =0; k < 3; k++)
	  {
	    dir[k] = -P[i].g.GravAccel[k]/(grav + 1.e-40);
	  }
	if (SphP[i].d.Density*All.UnitDensity_in_cgs > All.PebbleInjectionDens) 
	  {
	    probability_of_birth = 0.;
	  }
	else
	  {
	    if (fabs(P[i].Pos[2]) > All.PebbleInjectionSlab) probability_of_birth = 0.;
//	    local_mgas += P[i].Mass; /* mass of gas satisfying pebble birth criteria */
	  }
	if (rand_peb < probability_of_birth)
	  {
	    P[NumPart + stars_spawned] = P[i];
	    P[NumPart + stars_spawned].Type = 2;
	    for(j = 0; j < 3; j++)
	      {
		/* add a small random shift to the initial pebble particle position */
		/*	P[NumPart + stars_spawned].Pos[j] += 0.02 *(-2 + gsl_rng_uniform(random_generator))
		 * P[i].Hsml;*/
			P[NumPart + stars_spawned].Pos[j] += dir[j] * 0.004;
	      }

	    
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
	    P[NumPart + stars_spawned].Mass = peb_init_mass;
	    P[NumPart + stars_spawned].DustRadius = All.InitialDustRadius;
	    PPP[NumPart + stars_spawned].DustRadius = P[NumPart + stars_spawned].DustRadius;
	    force_add_star_to_tree(i, NumPart + stars_spawned);
	    stars_spawned++;
	  }
      }
#endif
  }
  MPI_Allreduce(&stars_spawned, &tot_spawned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
/*  MPI_Reduce(&local_mgas, &total_mgas_give_birth, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&total_mgas, &total_mgas_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (ThisTask == 0) printf("\n total_mgas_give_birth %g total_mgas_all %g ... \n", total_mgas, total_mgas_all); */
  if(tot_spawned > 0)
  {
      All.TotNumPart += tot_spawned;
      NumPart += stars_spawned;
      //      printf("NEW PEBBLES BORN = %d \n", tot_spawned);
  }


}

#endif /* end of DUST_PEBBLES_BORN */


#ifdef DUST_SINK_ON_FLY
void dust_sink_formation(void) {
  int i;
  double flag, crit_dens, rdist;
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if (P[i].Type == 2) {

      flag = 1.;
      /* rdist = pow((P[i].Pos[0] - All.xbh)*(P[i].Pos[0] - All.xbh) + (P[i].Pos[1] - All.ybh)*(P[i].Pos[1] - All.ybh) +  */
      /* 		     (P[i].Pos[2] - All.zbh)*(P[i].Pos[2] - All.zbh), 0.5); */
      /* crit_dens = All.CritPhysDensity * All.SMBHmass/(2.*M_PI*pow(rdist, 3) + 1.e-30); */

      crit_dens = All.CritDustDensity/All.UnitDensity_in_cgs;

      /*      printf("DUST PARTICLE i %d, density %g    while crit dens = %g \n", i, P[i].d7.dDUST_particle_density, crit_dens);
	      fflush(stdout);*/

      if (P[i].d7.dDUST_particle_density >= crit_dens) flag = 0.;

      /* SINK FORMATION */
      if (flag == 0)
	{
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


	  /*	  printf("DUST PARTICLE i %d, CONVERTED INTO A SINK PARTICLE!!!!  Mass %g \n", i, P[i].Mass);
		  fflush(stdout);*/
	}
  }

}
#endif //end of DUST_SINK_ON_FLY



#endif /* DUST */
