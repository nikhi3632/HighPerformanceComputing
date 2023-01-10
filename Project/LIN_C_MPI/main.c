/*This code solves the linear shallow water equations */

#include <stdlib.h>                                                             
#include <stdio.h>                                                              
#include <math.h>                                                               
#include <time.h>                                                               
#include "flux.h"
#include "tendency.h"
#include "mpi_fun.h"
#include "input_output.h"
#include "mpi.h"

/*******************************************/
int main(int argc, char *argv[])
{
  const int        nt     = 1000 ;
  const int        nrk3   = 3 ;
  const int        nx     = 600 ;
  const int        ny     = nx ;
  const double     g      = 9.81 ;
  const double     pi     =  acos(-1.0) ;
  const double     c      = 30.0 ;
  const double     h0     = c*c/g ;
  const double     cour   = 0.40 ;
  const double     dx     = 100.0 ;
  const double     dy     = dx ;
  const double     dt     = cour*dx/c/sqrt(2.0) ;
  const double     omega  = 4.e-3 ;
  const double     amp    = 10.0 ;
  /************************************************/
  int it,ir,ix,iy,ixs,iys;
  double xx,yy,at;
  int irec = 0;
  clock_t tic , toc;                                                          
  /************************************************/
  int  nproc,myid,j_st,j_en,my,iyy,iy_p;
  /************************************************/
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  int *dims   = (int *) malloc( nproc* sizeof(int));                 
  int *disps  = (int *) malloc( nproc* sizeof(int));                 

  init_myid(myid,&j_st,&j_en,&my,dims,disps,ny,nproc);

  printf(" %6d, %6d, %6d, %6d\n", myid,j_st,j_en,my);

  /************************************************/
  /************************************************/
  double **h_a = (double **) malloc((my+4) * sizeof(double *));                  
  h_a[0] = (double *) malloc((my+4) * nx * sizeof(double));                      
  for(iy = 1; iy < (my+4); iy++) h_a[iy] = h_a[iy-1] + nx;                       
                                                                            
  double **h_b = (double **) malloc((my+4) * sizeof(double *));                  
  h_b[0] = (double *) malloc((my+4) * nx * sizeof(double));                      
  for(iy = 1; iy < (my+4); iy++) h_b[iy] = h_b[iy-1] + nx;                       
                                                                            
  double **u_a = (double **) malloc((my+4) * sizeof(double *));                  
  u_a[0] = (double *) malloc((my+4) * nx * sizeof(double));                      
  for(iy = 1; iy < (my+4); iy++) u_a[iy] = u_a[iy-1] + nx;                       
                                                                            
  double **u_b = (double **) malloc((my+4) * sizeof(double *));                  
  u_b[0] = (double *) malloc((my+4) * nx * sizeof(double));                      
  for(iy = 1; iy < (my+4); iy++) u_b[iy] = u_b[iy-1] + nx;                       
                                                                            
  double **v_a = (double **) malloc((my+4) * sizeof(double *));                  
  v_a[0] = (double *) malloc((my+4) * nx * sizeof(double));                      
  for(iy = 1; iy < (my+4); iy++) v_a[iy] = v_a[iy-1] + nx;                       
                                                                            
  double **v_b = (double **) malloc((my+4) * sizeof(double *));                  
  v_b[0] = (double *) malloc((my+4) * nx * sizeof(double));                      
  for(iy = 1; iy < (my+4); iy++) v_b[iy] = v_b[iy-1] + nx;                       
                                                                            
  double **fl_x = (double **) malloc((my+4) * sizeof(double *));                 
  fl_x[0] = (double *) malloc((my+4) * nx * sizeof(double));                     
  for(iy = 1; iy < (my+4); iy++) fl_x[iy] = fl_x[iy-1] + nx;                     
                                                                            
  double **fl_y = (double **) malloc((my+4) * sizeof(double *));                 
  fl_y[0] = (double *) malloc((my+4) * nx * sizeof(double));                     
  for(iy = 1; iy < (my+4); iy++) fl_y[iy] = fl_y[iy-1] + nx;                     
                                                                            
  double **dudx = (double **) malloc((my+4) * sizeof(double *));                 
  dudx[0] = (double *) malloc((my+4) * nx * sizeof(double));                     
  for(iy = 1; iy < (my+4); iy++) dudx[iy] = dudx[iy-1] + nx;                     
                                                                            
  double **dvdy = (double **) malloc((my+4) * sizeof(double *));                 
  dvdy[0] = (double *) malloc((my+4) * nx * sizeof(double));                     
  for(iy = 1; iy < (my+4); iy++) dvdy[iy] = dvdy[iy-1] + nx;                     
                                                                            
  double **dhdx = (double **) malloc((my+4) * sizeof(double *));                 
  dhdx[0] = (double *) malloc((my+4) * nx * sizeof(double));                     
  for(iy = 1; iy < (my+4); iy++) dhdx[iy] = dhdx[iy-1] + nx;                     
                                                                            
  double **dhdy = (double **) malloc((my+4) * sizeof(double *));                 
  dhdy[0] = (double *) malloc((my+4) * nx * sizeof(double));                     
  for(iy = 1; iy < (my+4); iy++) dhdy[iy] = dhdy[iy-1] + nx;                     
                                                                            
  double **rforce = (double **) malloc((my+4) * sizeof(double *));               
  rforce[0] = (double *) malloc((my+4) * nx * sizeof(double));                   
  for(iy = 1; iy < (my+4); iy++) rforce[iy] = rforce[iy-1] + nx; 
                                                                                
  /************************************************/
  /************************************************/
  MPI_Barrier(MPI_COMM_WORLD);
  tic = clock();
  ixs = 300; iys = 300;

  iy_p = 0;
  int id_print = 0 ;
  if (iys >= j_st && iys <= j_en)
  {  
    id_print = 1;
    iy_p = iys-j_st+2;
    printf("%8d, %8d, %8d, %8d\n", myid,id_print,iy_p,iys);
  }

  for(iy = j_st; iy <= j_en; iy++)
  { 
    iyy = iy - j_st +2 ;
    for(ix = 0; ix <= nx-1; ix++)
    { 
      xx = (double) ix - (double) ixs ;
      yy = (double) iy - (double) iys ;
      rforce[iyy][ix]= exp(-(xx*xx+yy*yy)/2.);
    }
  }

  for(iy = 0; iy <= my+3; iy++)
  { 
    for(ix = 0; ix <= nx-1; ix++)
    { 
      h_a[iy][ix]=0. ; h_b[iy][ix]=0. ;
      u_a[iy][ix]=0. ; u_b[iy][ix]=0. ;
      v_a[iy][ix]=0. ; v_b[iy][ix]=0. ;
      dudx[iy][ix]=0.; dvdy[iy][ix]=0.;
      dhdx[iy][ix]=0.; dhdy[iy][ix]=0.;
    }
  }

  ///*********************** Begin rk3 loop *************************/
  for(it = 1; it <= nt; it++) /* start rk3 integration */
  {
    for(ir = 1; ir <= nrk3; ir++)  /* start rk3 substeb */
    {
      switch(ir)
      {
        case 1:
          at=dt/3.;
          break;
        case 2:
          at=dt/2.;
          break;
        case 3:
          at=dt;
          break;
        default:
          printf("invalid rk3 substep");
          MPI_Abort(MPI_COMM_WORLD, 99);
          break;
      }

      distru2(h_b,myid,nx,my,nproc);
      distru2(u_b,myid,nx,my,nproc);
      distru2(v_b,myid,nx,my,nproc);

      flux_x(u_b,fl_x,nx,my);
      dudx_at_h(fl_x,dudx,dx,nx,my);

      flux_y(v_b,fl_y,nx,my);
      dvdy_at_h(fl_y,dvdy,dy,nx,my);

      flux_x(h_b,fl_x,nx,my);
      dhdx_at_u(fl_x,dhdx,dx,nx,my);

      flux_y(h_b,fl_y,nx,my);
      dhdy_at_v(fl_y,dhdy,dy,nx,my);

      for(iy = 2; iy <= my+1; iy++)
      { 
        for(ix = 0; ix <= nx-1; ix++)
        { 
          h_b[iy][ix] = h_a[iy][ix] -at*h0*(dudx[iy][ix]+dvdy[iy][ix]);
          u_b[iy][ix] = u_a[iy][ix] -at*g*dhdx[iy][ix] ;
          v_b[iy][ix] = v_a[iy][ix] -at*g*dhdy[iy][ix] ;
        }
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }  /* end rk3 substerp */

    for(iy = 2; iy <= my+1; iy++)
    { 
      for(ix = 0; ix <= nx-1; ix++)
      { 
        h_b[iy][ix] = h_b[iy][ix] + rforce[iy][ix]*amp*sin(2.*pi*omega*dt*(double)(it-1)) ;
        u_a[iy][ix] = u_b[iy][ix];
        v_a[iy][ix] = v_b[iy][ix];
        h_a[iy][ix] = h_b[iy][ix];
      }
    }

    if( id_print == 1)
    {
      printf("%10d,%12.6f\n", it,h_b[iy_p][ixs]);                              
    }

    if ((it)%100 == 0)
    {
      irec++ ;
      MPI_Barrier(MPI_COMM_WORLD);
      write_to_file(h_b,myid,nx,my,ny,irec,dims,nproc) ;
    }
  }  /* end  rk3 integration */
    
  ///*********************** END rk3 loop *************************/
  free (dims) ;                                                            
  free (disps) ;                                                           
  toc = clock();
  if(myid ==0) 
  {
    printf("Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
  }     

  free(h_a[0]); free(h_a);                                                   
  free(h_b[0]); free(h_b);                                                   
                                                                            
  free(u_a[0]); free(u_a);                                                   
  free(u_b[0]); free(u_b);                                                   
                                                                            
  free(v_a[0]); free(v_a);                                                   
  free(v_b[0]); free(v_b);                                                   
                                                                            
  free(dudx[0]); free(dudx);                                                 
  free(dvdy[0]); free(dvdy);                                                 
                                                                            
  free(dhdx[0]); free(dhdx);                                                 
  free(dhdy[0]); free(dhdy);                                                 
                                                                            
  free(fl_x[0]); free(fl_x);                                                 
  free(fl_y[0]); free(fl_y);                                                 
                                                                            
  free(rforce[0]); free(rforce);                                             

  MPI_Finalize();
  return 0;
}
/******************************************************/
/******************************************************/
