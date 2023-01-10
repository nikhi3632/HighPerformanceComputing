
#include <stdio.h>

/**************************************************************/
void init_myid(size_t * j_st, size_t * j_en, size_t * j_dim, size_t ny, size_t nproc) 
{ 
  size_t  ip,ndim,end_dim ;

  ndim=ny/nproc;

  for(ip = 0; ip <= nproc-1; ip++) 
  {
    j_dim[ip] = ndim;
  }

  if (ny%nproc != 0) 
  {
    for(ip = 0; ip < ny%nproc; ip++) 
    {
      j_dim[ip] = j_dim[ip]+1;
    }
  }

  j_st[0]=0;
  j_en[0]=j_dim[0]-1;

  if (nproc > 1)
  {
    for(ip = 1; ip <= nproc-1; ip++) 
    {
      j_st[ip] = j_en[ip-1]+1;
      j_en[ip] = j_st[ip]+j_dim[ip]-1;
    }
  }
} 
/**************************************************************/
