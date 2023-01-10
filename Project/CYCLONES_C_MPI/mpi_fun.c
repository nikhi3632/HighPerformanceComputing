
#include "mpi.h"

/**************************************************************/
void init_myid(int myid,int * j_st, int * j_en, int * my, int * dims,int * disps,int ny,int nproc)
{ 
  int ip ;

  *my = ny/nproc ;
  if (myid < ny%nproc)
  {
    *my=*my+1;
  }

  MPI_Allgather(my,1,MPI_INT,dims,1,MPI_INT,MPI_COMM_WORLD);

  disps[0]=0;

  if(nproc > 1)  
  {
    for (ip=1;ip<=nproc-1;ip++)
    {
      disps[ip]=disps[ip-1]+dims[ip-1];
    }
  }

  *j_st = disps[myid];
  *j_en = *j_st + *my -1;
} 
/**************************************************************/
/**************************************************************/
void distru2(double ** ru, int myid, int nx, int my, int nproc)                            
{
  double rin_s[2][nx], rin_r[2][nx] ;
  MPI_Request request[2];
  MPI_Status status[2];
  int ix , pro_s,pro_r;
  //============
  //------------   send_up
  pro_s = myid + 1;
  pro_r = myid - 1;
  if (myid == nproc-1) 
  {
    pro_s = 0;
  }
  if (myid == 0) 
  {
    pro_r = nproc-1;
  }

  for(ix = 0; ix <= nx-1; ix++)
  {
    rin_s[1][ix] = ru[my+1][ix];
    rin_s[0][ix] = ru[my][ix];
  }

  MPI_Issend(rin_s,nx*2,MPI_DOUBLE, pro_s,99,MPI_COMM_WORLD,&request[0]);
  MPI_Irecv (rin_r,nx*2,MPI_DOUBLE, pro_r,99,MPI_COMM_WORLD,&request[1]);
  MPI_Waitall(2,request,status);

  for(ix = 0; ix <= nx-1; ix++)
  {
    ru[1][ix] = rin_r[1][ix];
    ru[0][ix] = rin_r[0][ix];
  }

  //------------   send_down
  pro_s = myid - 1;
  pro_r = myid + 1;
  if (myid == nproc-1)
  {
    pro_r = 0;
  }
  if (myid == 0) 
  { 
    pro_s = nproc-1;
  }

  for(ix = 0; ix <= nx-1; ix++)
  {
    rin_s[0][ix] = ru[2][ix];
    rin_s[1][ix] = ru[3][ix];
  }

  MPI_Issend(rin_s,nx*2,MPI_DOUBLE, pro_s,99,MPI_COMM_WORLD,&request[0]);
  MPI_Irecv (rin_r,nx*2,MPI_DOUBLE, pro_r,99,MPI_COMM_WORLD,&request[1]);
  MPI_Waitall(2,request,status);

  for(ix = 0; ix <= nx-1; ix++)
  {
    ru[my+2][ix] = rin_r[0][ix];
    ru[my+3][ix] = rin_r[1][ix];
  }
}
/**************************************************************/
/**************************************************************/
