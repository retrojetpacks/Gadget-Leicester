#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"







/* This routine allocates memory for 
 * particle storage, both the collisionless and the SPH particles.
 * The memory for the ordered binary tree of the timeline
 * is also allocated.
 */
void allocate_memory(void)
{
  size_t bytes;
  double bytes_tot = 0;


  Exportflag = (int *) mymalloc(NTask * sizeof(int));
  Exportindex = (int *) mymalloc(NTask * sizeof(int));
  Exportnodecount = (int *) mymalloc(NTask * sizeof(int));

  Send_count = (int *) mymalloc(sizeof(int) * NTask);
  Send_offset = (int *) mymalloc(sizeof(int) * NTask);
  Recv_count = (int *) mymalloc(sizeof(int) * NTask);
  Recv_offset = (int *) mymalloc(sizeof(int) * NTask);
  Sendcount_matrix = (int *) mymalloc(sizeof(int) * NTask * NTask);

  NextActiveParticle = (int *) mymalloc(bytes = All.MaxPart * sizeof(int));
  bytes_tot += bytes;

#ifdef  VIRTUAL
  NextPhoton = (int *) mymalloc(bytes = All.MaxPart * sizeof(int));
  bytes_tot += bytes;

#ifdef OTHIN_ACCELERATOR
  NextSinkParticle = (int *) mymalloc(bytes = All.MaxPart * sizeof(int));
  bytes_tot += bytes;
#endif

#endif


  NextInTimeBin = (int *) mymalloc(bytes = All.MaxPart * sizeof(int));
  bytes_tot += bytes;

  PrevInTimeBin = (int *) mymalloc(bytes = All.MaxPart * sizeof(int));
  bytes_tot += bytes;


  if(All.MaxPart > 0)
    {
      if(!(P = (struct particle_data *) mymalloc(bytes = All.MaxPart * sizeof(struct particle_data))))
	{
	  printf("failed to allocate memory for `P' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(1);
	}
      bytes_tot += bytes;

      if(ThisTask == 0)
	printf("\nAllocated %g MByte for particle storage.\n\n", bytes_tot / (1024.0 * 1024.0));
    }

  if(All.MaxPartSph > 0)
    {
      bytes_tot = 0;

      if(!
	 (SphP =
	  (struct sph_particle_data *) mymalloc(bytes = All.MaxPartSph * sizeof(struct sph_particle_data))))
	{
	  printf("failed to allocate memory for `SphP' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(1);
	}
      bytes_tot += bytes;

      if(ThisTask == 0)
	printf("Allocated %g MByte for storage of SPH data.\n\n", bytes_tot / (1024.0 * 1024.0));

    }
}
