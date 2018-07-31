#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>


#include "allvars.h"
#include "proto.h"

#ifndef DEBUG
#define NDEBUG
#endif
#include <assert.h>

#ifdef RADTRANSFER

#define f_ij -10.0
#define f_ii -10.0/3.0
#define sigma 6.3e-18*All.UnitLength_in_cm/All.UnitLength_in_cm
#define alpha 2.59e-13
#define eV_to_erg 1.60184e-12
#define eV_to_Hz 2.41838e14
#define c_light 2.9979e10/All.UnitVelocity_in_cm_per_s


/* this function loops over all particles and computes the mean intensity for the particles 
that do not need to be exported to other processors. the particle that need to be exportes 
are then exported and their thei mean intensity is calculated as well */
void radiative_transfer(void)
{
  int i, j, k, ngrp, dummy;
  int sendTask, recvTask, nexport = 0, nimport, place;
  double residue = 0, sumresidue = 0;

  /* allocate buffers to arrange communication */

  Ngblist = (int *) mymalloc(NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct eddingtondata_in) +
					     sizeof(struct eddingtondata_out) +
					     sizemax(sizeof(struct eddingtondata_in),
						     sizeof(struct eddingtondata_out))));
  DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));

  for(j = 0; j < NTask; j++)
    {
      Send_count[j] = 0;
      Exportflag[j] = -1;
    }




  /* do local particles and prepare export list */
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)
	{
	  if(radtransfer_evaluate(i, 0, &nexport, Send_count) < 0)
	    break;
	  if(nexport == 0)
	    {
	      residue += evaluate(i);
	    }
	}
    }

#ifdef MYSORT
  mysort_dataindex(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#else
  qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#endif

  MPI_Allgather(Send_count, NTask, MPI_INT, Sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);

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

  EddingtonDataGet = (struct eddingtondata_in *) mymalloc(nimport * sizeof(struct eddingtondata_in));
  EddingtonDataIn = (struct eddingtondata_in *) mymalloc(nexport * sizeof(struct eddingtondata_in));

  /* prepare particle data for export */

  for(j = 0; j < nexport; j++)
    {
      place = DataIndexTable[j].Index;

      for(k = 0; k < 3; k++)
	{
	  EddingtonDataIn[j].Pos[k] = P[place].Pos[k];
	  EddingtonDataIn[j].ET[k] = SphP[place].ET[k];
	  EddingtonDataIn[j].ET[k + 3] = SphP[place].ET[k + 3];
	}
      EddingtonDataIn[j].Hsml = PPP[place].Hsml;
      EddingtonDataIn[j].nHI = SphP[place].nHI;

      memcpy(EddingtonDataIn[j].NodeList,
	     DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
    }

  /* exchange particle data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
	{
	  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
	    {
	      /* get the particles */
	      MPI_Sendrecv(&EddingtonDataIn[Send_offset[recvTask]],
			   Send_count[recvTask] * sizeof(struct eddingtondata_in), MPI_BYTE,
			   recvTask, TAG_HYDRO_A,
			   &EddingtonDataGet[Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(struct eddingtondata_in), MPI_BYTE,
			   recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }

  myfree(EddingtonDataIn);
  EddingtonDataResult = (struct eddingtondata_out *) mymalloc(nimport * sizeof(struct eddingtondata_out));
  EddingtonDataOut = (struct eddingtondata_out *) mymalloc(nexport * sizeof(struct eddingtondata_out));

  /* now do the particles that were sent to us */
  for(j = 0; j < nimport; j++)
    {
      radtransfer_evaluate(j, 1, &dummy, &dummy);
    }

  /* get the result */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
	{
	  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
	    {
	      /* send the results */
	      MPI_Sendrecv(&EddingtonDataResult[Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(struct eddingtondata_out),
			   MPI_BYTE, recvTask, TAG_HYDRO_B,
			   &EddingtonDataOut[Send_offset[recvTask]],
			   Send_count[recvTask] * sizeof(struct eddingtondata_out),
			   MPI_BYTE, recvTask, TAG_HYDRO_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }

  /* add the result to the local particles */
  for(j = 0; j < nexport; j++)
    {
      place = DataIndexTable[j].Index;
      SphP[place].P1 += EddingtonDataOut[j].P1;
      SphP[place].D1 += EddingtonDataOut[j].D1;
    }

  for(j = 0; j < nexport; j++)
    {
      place = DataIndexTable[j].Index;
      residue += evaluate(j);
    }

  myfree(EddingtonDataOut);
  myfree(EddingtonDataResult);
  myfree(EddingtonDataGet);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);

  MPI_Allreduce(&residue, &sumresidue, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  All.Residue = sumresidue / FLT(All.TotN_gas);

}

/* this function evaluates the coefficients  D1 and P1 needed for the Gauss-Seidel method integration */
int radtransfer_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int startnode, numngb, listindex = 0;
  int j, k, n;
  MyFloat ET[6] = { 0, 0, 0, 0, 0, 0 }, ET_j[6] =
  {0, 0, 0, 0, 0, 0};
  MyFloat kappa_i, kappa_j, h_i, RadJ_j, kappa_ij;
  MyDouble *pos;
  MyFloat mass, rho;
  MyFloat D1 = 0.0, P1 = 0.0, a;

  double dx, dy, dz;
  double h_j, hinv, hinv4;
  double dwk_i;
  double r, r2, r4, u;

  if(All.ComovingIntegrationOn)
    a = All.Time;
  else
    a = 1.0;

  if(mode == 0)
    {
      for(k = 0; k < 6; k++)
	{
	  ET[k] = SphP[target].ET[k];
	}
      pos = P[target].Pos;
      h_i = PPP[target].Hsml;
      kappa_i =
	(a / c_light) * (SphP[target].nHI * sigma - All.kappaMean + (All.JMean / All.RadIMean));

    }
  else
    {
      for(k = 0; k < 6; k++)
	{
	  ET[k] = EddingtonDataGet[target].ET[k];
	}
      pos = EddingtonDataGet[target].Pos;
      h_i = EddingtonDataGet[target].Hsml;
      kappa_i =
	(a / c_light) * (EddingtonDataGet[target].nHI * sigma - All.kappaMean +
				(All.JMean / All.RadIMean));
    }

  hinv = 1.0 / h_i;
  hinv4 = hinv * hinv * hinv * hinv;

  if(mode == 0)
    {
      startnode = All.MaxPart;
    }
  else
    {
      startnode = EddingtonDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;
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
	      dx = pos[0] - P[j].Pos[0];
	      dy = pos[1] - P[j].Pos[1];
	      dz = pos[2] - P[j].Pos[2];
	      r2 = dx * dx + dy * dy + dz * dz;
	      r = sqrt(r2);
	      if(r < h_i)
		{
		  r4 = r2 * r2;
		  h_j = PPP[j].Hsml;
		  mass = P[j].Mass;
		  rho = SphP[j].d.Density;
		  RadJ_j = SphP[j].RadJ;
		  kappa_j =
		    (a / c_light) * (SphP[j].nHI * sigma - All.kappaMean + (All.JMean / All.RadIMean));

		  for(k = 0; k < 6; k++)
		    {
		      ET_j[k] = SphP[j].ET[k];
		    }
		  kappa_ij = 1.0 / kappa_i + 1.0 / kappa_j;

		  if(r)
		    {
		      u = r * hinv;

		      if(u < 0.5)
			dwk_i = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
		      else
			dwk_i = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
		    }
		  else
		    dwk_i = 0;

		  if(r4)
		    {
		      D1 += f_ii * (mass / rho) * ET[0] * kappa_ij * r * dwk_i * dx * dx / r4
			+ f_ii * (mass / rho) * ET[1] * kappa_ij * r * dwk_i * dy * dy / r4
			+ f_ii * (mass / rho) * ET[2] * kappa_ij * r * dwk_i * dz * dz / r4
			+ 2.0 * f_ij * (mass / rho) * ET[3] * kappa_ij * r * dwk_i * dx * dy / r4
			+ 2.0 * f_ij * (mass / rho) * ET[4] * kappa_ij * r * dwk_i * dy * dz / r4
			+ 2.0 * f_ij * (mass / rho) * ET[5] * kappa_ij * r * dwk_i * dz * dx / r4;
		      
		      P1 += f_ii * (mass / rho) * RadJ_j * ET_j[0] * kappa_ij * r * dwk_i * dx * dx / r4
			+ f_ii * (mass / rho) * RadJ_j * ET_j[1] * kappa_ij * r * dwk_i * dy * dy / r4
			+ f_ii * (mass / rho) * RadJ_j * ET_j[2] * kappa_ij * r * dwk_i * dz * dz / r4
			+ 2.0 * f_ij * (mass / rho) * RadJ_j * ET_j[3] * kappa_ij * r * dwk_i * dx * dy / r4
			+ 2.0 * f_ij * (mass / rho) * RadJ_j * ET_j[4] * kappa_ij * r * dwk_i * dy * dz / r4
			+ 2.0 * f_ij * (mass / rho) * RadJ_j * ET_j[5] * kappa_ij * r * dwk_i * dz * dx / r4;
		      
		    }
		}
	    }
	}
    }


  if(mode == 1)
    {
      listindex++;
      if(listindex < NODELISTLENGTH)
	{
	  startnode = EddingtonDataGet[target].NodeList[listindex];
	  if(startnode >= 0)
	    startnode = Nodes[startnode].u.d.nextnode;
	}
    }


  if(mode == 0)
    {
      SphP[target].P1 = (c_light / 4.0 / a) * P1;
      SphP[target].D1 = (c_light / 4.0 / a) * D1;
    }
  else
    {
      EddingtonDataResult[target].P1 = (c_light / 4.0 / a) * P1;
      EddingtonDataResult[target].D1 = (c_light / 4.0 / a) * D1;
    }

  return 0;
}


/* this function evaluates the radiation mean intensity using the Gauss-Seidel method */
double evaluate(int target)
{
  MyFloat RadJOld, dt, P2, D2;
  double residue = 0;

  dt = All.TimeStep;
  D2 = 0.5 * (SphP[target].nHI * sigma - All.kappaMean + (All.JMean / All.RadIMean));
  P2 = 0.5 * (SphP[target].Je / All.RadIMean);

  RadJOld = SphP[target].RadJ;
  SphP[target].RadJ = (SphP[target].RadJ + dt * (SphP[target].P1 + P2)) / (1.0 + dt * (SphP[target].D1 + D2));

  residue = fabs(SphP[target].RadJ - RadJOld);

  All.RadIMean += All.RadIMean * (SphP[target].RadJ - RadJOld) / FLT(All.TotN_gas);
  All.kappaMean += SphP[target].nHI * sigma * (SphP[target].RadJ - RadJOld) / FLT(All.TotN_gas);

  //printf("%g %g %i \n", SphP[target].RadJ + dt * (SphP[target].P1 + P2), 1.0 + dt * (SphP[target].D1 + D2), target);

  return residue;
}

/* this function calculates the mean radiation intensity, absorption coefficient and emission coefficient */
void mean(void)
{
  int i;
  double RadI = 0, kappa = 0, je = 0;
  double sumRadI = 0, sumkappa = 0, sumJ = 0;
  //double Vol = 0, sumVol = 0, sumVolAll;

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)
	{
	  //Vol = P[i].Mass / SphP[i].d.Density;
	  RadI += SphP[i].RadJ * All.RadIMean;
	  kappa += SphP[i].nHI * sigma * SphP[i].RadJ;
	  je += SphP[i].Je;
	  //sumVol += Vol;
	}
    }


  MPI_Allreduce(&RadI, &sumRadI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&je, &sumJ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&kappa, &sumkappa, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  //MPI_Allreduce(&sumVol, &sumVolAll, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  //All.TotVol = sumVolAll;
  All.RadIMean = sumRadI / FLT(All.TotN_gas);
  All.JMean = sumJ / FLT(All.TotN_gas);
  All.kappaMean = sumkappa / FLT(All.TotN_gas);

  //printf("%f %f %f \n", All.RadIMean, All.JMean, All.kappaMean);

}

void set_simple_inits(void)
{
  int i;

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)
	{
	  /* in code units */
	  if(All.ComovingIntegrationOn)
	    SphP[i].nHI = SphP[i].d.Density * All.UnitMass_in_g / All.HubbleParam / PROTONMASS /
	      All.Time / All.Time / All.Time;
	  else
	    SphP[i].nHI = 1.e-3 * All.UnitLength_in_cm * All.UnitLength_in_cm * All.UnitLength_in_cm; 
	  //SphP[i].d.Density * All.UnitMass_in_g / PROTONMASS;
	  SphP[i].nH = SphP[i].nHI;
	  SphP[i].x = 0.0;
	  SphP[i].nHII = 0.0;
	  SphP[i].n_elec = 0.0;
	  SphP[i].RadJ = 1.0 / FLT(All.TotN_gas);
	}
    }

  mean();
  All.RadIMean = 1.0;
}

void simple_output(void)
{
  char buf[100];
  int i;
  double fac, x, y, z, dx, dy, dz, r2, nHII, xx;
  MyFloat dt, fact, facv, a;

  dt = All.TimeStep;

  if(All.ComovingIntegrationOn)
    {
      fact = All.UnitTime_in_s / All.HubbleParam;
      facv = All.UnitLength_in_cm * All.UnitLength_in_cm * All.UnitLength_in_cm / 
	All.HubbleParam / All.HubbleParam / All.HubbleParam;
    }
  else
    {
      fact = All.UnitTime_in_s;
      facv = All.UnitLength_in_cm * All.UnitLength_in_cm * All.UnitLength_in_cm;
    }
  
  if(All.ComovingIntegrationOn)
    a = All.Time;
  else
    a = 1.0;

  sprintf(buf, "%s%s%i%s%f%s", All.OutputDir, "radtransfer_", ThisTask, "_", All.Time, ".txt");
  FdRadtransfer = fopen(buf, "wa");
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)
	fprintf(FdRadtransfer, "%f %f %f %f \n", P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], SphP[i].x);
    }

  
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)
	{
	  dx = x - P[i].Pos[0];
	  dy = y - P[i].Pos[1];
	  dz = z - P[i].Pos[2];
	  r2 = dx*dx + dy*dy + dz*dz;
	  fac = (a/c_light)*SphP[i].Je/All.RadIMean/4.0/M_PI/r2;

	  nHII = SphP[i].nHII + dt * All.UnitEnergy_in_cgs * fac * All.RadIMean / 13.6 / eV_to_erg *
	    sigma * SphP[i].nHI 
	    - alpha * fact * dt * SphP[i].nHII * SphP[i].nHII / facv;;
	  if(nHII > SphP[i].nH)
	    nHII = SphP[i].nH;
	  xx = (SphP[i].nH-nHII)/SphP[i].nH;
	  fprintf(FdRadtransfer, "%f %f %f %f \n", P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], xx);
	}
    }
  
  fclose(FdRadtransfer);
}


void update_nHI(void)
{
  int i;
  MyFloat dt, fact, facv;

  dt = All.TimeStep;

  if(All.ComovingIntegrationOn)
    {
      fact = All.UnitTime_in_s / All.HubbleParam;
      facv = All.UnitLength_in_cm * All.UnitLength_in_cm * All.UnitLength_in_cm / 
	All.HubbleParam / All.HubbleParam / All.HubbleParam;
    }
  else
    {
 
      fact = All.UnitTime_in_s;
      facv = All.UnitLength_in_cm * All.UnitLength_in_cm * All.UnitLength_in_cm;
    }

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)
	{
	  /* hydrogen after photoionization and recombination - monochromatic light */
	  SphP[i].nHII +=
	    dt * All.UnitEnergy_in_cgs * SphP[i].RadJ * All.RadIMean / 13.6 / eV_to_erg *
	    sigma * SphP[i].nHI 
	    - alpha * fact * dt * SphP[i].n_elec * SphP[i].n_elec / facv;
	  
	  if(SphP[i].nHII > SphP[i].nH)
	    SphP[i].nHII = SphP[i].nH;
	  
	  if(SphP[i].nHII < 0.0)
	    SphP[i].nHII = 0;
	  
	  SphP[i].nHI = SphP[i].nH - SphP[i].nHII;
	  SphP[i].n_elec = SphP[i].nHII;
	  SphP[i].x = SphP[i].nHI / SphP[i].nH;
	
	  
	  printf("%g %g %g %g\n", SphP[i].nHII, SphP[i].nHI, SphP[i].nH/facv, SphP[i].x);
	  printf("%g \n", (1.0/(alpha*SphP[i].nH/facv))/3600/24/365/1.e6);
      
	}

    }
  
}

#endif
