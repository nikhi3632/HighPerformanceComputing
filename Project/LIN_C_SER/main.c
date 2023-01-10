/*This code solves the linear shallow water equations */

#include <stdlib.h>                                                             
#include <stdio.h>                                                              
#include <math.h>                                                               
#include <time.h>                                                               
#include "flux.h"
#include "tendency.h"
#include "input_output.h"

/*******************************************/
int main()
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
  FILE  *file_out;                                                            
  int it,ir,ix,iy,ixs,iys;
  double xx,yy,at;
  int irec = 0;
  clock_t tic , toc;
/************************************************/
  double **h_a = (double **) malloc(ny * sizeof(double *));                  
  h_a[0] = (double *) malloc(ny * nx * sizeof(double));                      
  for(iy = 1; iy < ny; iy++) h_a[iy] = h_a[iy-1] + nx;                       
                                                                            
  double **h_b = (double **) malloc(ny * sizeof(double *));                  
  h_b[0] = (double *) malloc(ny * nx * sizeof(double));                      
  for(iy = 1; iy < ny; iy++) h_b[iy] = h_b[iy-1] + nx;                       
                                                                            
  double **u_a = (double **) malloc(ny * sizeof(double *));                  
  u_a[0] = (double *) malloc(ny * nx * sizeof(double));                      
  for(iy = 1; iy < ny; iy++) u_a[iy] = u_a[iy-1] + nx;                       
                                                                            
  double **u_b = (double **) malloc(ny * sizeof(double *));                  
  u_b[0] = (double *) malloc(ny * nx * sizeof(double));                      
  for(iy = 1; iy < ny; iy++) u_b[iy] = u_b[iy-1] + nx;                       
                                                                            
  double **v_a = (double **) malloc(ny * sizeof(double *));                  
  v_a[0] = (double *) malloc(ny * nx * sizeof(double));                      
  for(iy = 1; iy < ny; iy++) v_a[iy] = v_a[iy-1] + nx;                       
                                                                            
  double **v_b = (double **) malloc(ny * sizeof(double *));                  
  v_b[0] = (double *) malloc(ny * nx * sizeof(double));                      
  for(iy = 1; iy < ny; iy++) v_b[iy] = v_b[iy-1] + nx;                       
                                                                    
  double **fl_x = (double **) malloc(ny * sizeof(double *));                 
  fl_x[0] = (double *) malloc(ny * nx * sizeof(double));                     
  for(iy = 1; iy < ny; iy++) fl_x[iy] = fl_x[iy-1] + nx;                     
                                                                            
  double **fl_y = (double **) malloc(ny * sizeof(double *));                 
  fl_y[0] = (double *) malloc(ny * nx * sizeof(double));                     
  for(iy = 1; iy < ny; iy++) fl_y[iy] = fl_y[iy-1] + nx;                     
                                                                            
  double **dudx = (double **) malloc(ny * sizeof(double *));                 
  dudx[0] = (double *) malloc(ny * nx * sizeof(double));                     
  for(iy = 1; iy < ny; iy++) dudx[iy] = dudx[iy-1] + nx;                     
                                                                            
  double **dvdy = (double **) malloc(ny * sizeof(double *));                 
  dvdy[0] = (double *) malloc(ny * nx * sizeof(double));                     
  for(iy = 1; iy < ny; iy++) dvdy[iy] = dvdy[iy-1] + nx;                     
                                                                            
  double **dhdx = (double **) malloc(ny * sizeof(double *));                 
  dhdx[0] = (double *) malloc(ny * nx * sizeof(double));                     
  for(iy = 1; iy < ny; iy++) dhdx[iy] = dhdx[iy-1] + nx;                     
                                                                            
  double **dhdy = (double **) malloc(ny * sizeof(double *));                 
  dhdy[0] = (double *) malloc(ny * nx * sizeof(double));                     
  for(iy = 1; iy < ny; iy++) dhdy[iy] = dhdy[iy-1] + nx;                     
                                                                            
  double **rforce = (double **) malloc(ny * sizeof(double *));               
  rforce[0] = (double *) malloc(ny * nx * sizeof(double));                   
  for(iy = 1; iy < ny; iy++) rforce[iy] = rforce[iy-1] + nx;                 

  /************************************************/
  /************************************************/
  tic = clock();
  ixs = 300; iys = 300;

  for(iy = 0; iy <= ny-1; iy++) 
  { 
    for(ix = 0; ix <= nx-1; ix++)
    { 
      xx = (double) ix - (double) ixs;
      yy = (double) iy - (double) iys;
      rforce[iy][ix]= exp(-(xx*xx+yy*yy)/2.);
    }
  }

  for(iy = 0; iy <= ny-1; iy++)
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

  /*********************** Begin rk3 loop *************************/
  for(it = 1; it <= nt; it++) /* start rk3 integration */
  {
    for(ir = 1; ir <= nrk3; ir++)  /* start rk3 substeb */
    {
      switch (ir)
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
        exit(1);
        break;
      }

      flux_x(u_b,fl_x,nx,ny);
      dudx_at_h(fl_x,dudx,dx,nx,ny) ;

      flux_y(v_b,fl_y,nx,ny);
      dvdy_at_h(fl_y,dvdy,dy,nx,ny) ;

      flux_x(h_b,fl_x,nx,ny);
      dhdx_at_u(fl_x,dhdx,dx,nx,ny) ;

      flux_y(h_b,fl_y,nx,ny);
      dhdy_at_v(fl_y,dhdy,dy,nx,ny) ;

      for(iy = 0; iy <= ny-1; iy++)
      { 
        for(ix = 0; ix <= nx-1; ix++)
        { h_b[iy][ix] = h_a[iy][ix] -at*h0*(dudx[iy][ix]+dvdy[iy][ix]);
          u_b[iy][ix] = u_a[iy][ix] -at*g*dhdx[iy][ix] ;
          v_b[iy][ix] = v_a[iy][ix] -at*g*dhdy[iy][ix] ;
        }
      }
    }  /* end rk3 substerp */

    for(iy = 0; iy <= ny-1; iy++)
    { 
      for(ix = 0; ix <= nx-1; ix++)
      { 
        h_b[iy][ix] = h_b[iy][ix] + rforce[iy][ix]*amp*sin(2.*pi*omega*dt*(double)(it-1)) ;
        u_a[iy][ix] = u_b[iy][ix];
        v_a[iy][ix] = v_b[iy][ix];
        h_a[iy][ix] = h_b[iy][ix];
      }
    }

    printf("%10d,%12.6f\n", it,h_b[iys][ixs]);                              
    /*********************** END rk3 loop *************************/
    if ((it)%100 == 0) 
    {   
      irec++ ;
      write_to_file(h_b,irec,nx,ny) ;
    }
  }  /* end  rk3 integration */
     
  toc = clock();                                                             
  printf("Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

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
                                                                            
  free(rforce[0]); free(rforce);                                             

  return 0;
}
/******************************************************/
/******************************************************/

