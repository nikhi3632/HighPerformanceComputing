/*This code solves the nonlinear shallow water equations in a rotating frame*/

#include <stdlib.h>                                                             
#include <stdio.h>                                                              
#include <math.h>                                                               
#include <time.h>                                                               
#include "flux.h"
#include "tendency.h"
#include "init.h"
#include "adv_diff_cor.h"
#include "input_output.h"
/*******************************************/
int main()
{
  const int        nt     = 5001 ;
  const int        nrk3   = 3 ;
  const int        nx     = 600 ;
  const int        ny     = nx ;
  const double     g      = 9.81 ;
  const double     c      = 60.0 ;
  const double     h0     = c*c/g ;
  const double     cour   = 0.40 ;
  const double     dx     = 25.e3 ;
  const double     dy     = dx ;
  const double     dt     = cour*dx/c/sqrt(2.0) ;
  const double     gam    = 0.02/dt ;
  /************************************************/
  int it,ir,ix,iy,ixs,iys;
  double at;
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

  double **fl_x = (double **) malloc(ny * sizeof(double *));                 
  fl_x[0] = (double *) malloc(ny * nx * sizeof(double));                     
  for(iy = 1; iy < ny; iy++) fl_x[iy] = fl_x[iy-1] + nx;                     
                                                                            
  double **fl_y = (double **) malloc(ny * sizeof(double *));                 
  fl_y[0] = (double *) malloc(ny * nx * sizeof(double));                     
  for(iy = 1; iy < ny; iy++) fl_y[iy] = fl_y[iy-1] + nx;                     

  double **adv_u = (double **) malloc(ny * sizeof(double *));                 
  adv_u[0] = (double *) malloc(ny * nx * sizeof(double));                     
  for(iy = 1; iy < ny; iy++) adv_u[iy] = adv_u[iy-1] + nx;                     

  double **adv_v = (double **) malloc(ny * sizeof(double *));                 
  adv_v[0] = (double *) malloc(ny * nx * sizeof(double));                     
  for(iy = 1; iy < ny; iy++) adv_v[iy] = adv_v[iy-1] + nx;                     

  double **adv_h = (double **) malloc(ny * sizeof(double *));                 
  adv_h[0] = (double *) malloc(ny * nx * sizeof(double));                     
  for(iy = 1; iy < ny; iy++) adv_h[iy] = adv_h[iy-1] + nx;                     

  double **diff_u = (double **) malloc(ny * sizeof(double *));                 
  diff_u[0] = (double *) malloc(ny * nx * sizeof(double));                     
  for(iy = 1; iy < ny; iy++) diff_u[iy] = diff_u[iy-1] + nx;                     

  double **diff_v = (double **) malloc(ny * sizeof(double *));                 
  diff_v[0] = (double *) malloc(ny * nx * sizeof(double));                     
  for(iy = 1; iy < ny; iy++) diff_v[iy] = diff_v[iy-1] + nx;                     

  double **diff_h = (double **) malloc(ny * sizeof(double *));                 
  diff_h[0] = (double *) malloc(ny * nx * sizeof(double));                     
  for(iy = 1; iy < ny; iy++) diff_h[iy] = diff_h[iy-1] + nx;                     

  double **tend_u = (double **) malloc(ny * sizeof(double *));                 
  tend_u[0] = (double *) malloc(ny * nx * sizeof(double));                     
  for(iy = 1; iy < ny; iy++) tend_u[iy] = tend_u[iy-1] + nx;                     

  double **tend_v = (double **) malloc(ny * sizeof(double *));                 
  tend_v[0] = (double *) malloc(ny * nx * sizeof(double));                     
  for(iy = 1; iy < ny; iy++) tend_v[iy] = tend_v[iy-1] + nx;                     

  double **tend_h = (double **) malloc(ny * sizeof(double *));                 
  tend_h[0] = (double *) malloc(ny * nx * sizeof(double));                     
  for(iy = 1; iy < ny; iy++) tend_h[iy] = tend_h[iy-1] + nx;                     

  double **sig_u = (double **) malloc(ny * sizeof(double *));                 
  sig_u[0] = (double *) malloc(ny * nx * sizeof(double));                     
  for(iy = 1; iy < ny; iy++) sig_u[iy] = sig_u[iy-1] + nx;                     

  double **sig_v = (double **) malloc(ny * sizeof(double *));                 
  sig_v[0] = (double *) malloc(ny * nx * sizeof(double));                     
  for(iy = 1; iy < ny; iy++) sig_v[iy] = sig_v[iy-1] + nx;                     

  double **sig_h = (double **) malloc(ny * sizeof(double *));                 
  sig_h[0] = (double *) malloc(ny * nx * sizeof(double));                     
  for(iy = 1; iy < ny; iy++) sig_h[iy] = sig_h[iy-1] + nx;                     

  double * f_u  = (double *) malloc(ny * sizeof(double));                            
  double * f_v  = (double *) malloc(ny * sizeof(double));                            

  /************************************************/
  /************************************************/
  absorb(sig_u,sig_v,sig_h,dt,nx,ny);
  etat0(u_a,v_a,h_a,f_u,f_v,dy,nx,ny);

  tic = clock();
  ixs = 300; iys = 300;
  
  for(iy = 0; iy <= ny-1; iy++)
  { 
    for(ix = 0; ix <= nx-1; ix++)
      { 
        h_b[iy][ix] = h_a[iy][ix];
        u_b[iy][ix] = u_a[iy][ix];
        v_b[iy][ix] = v_a[iy][ix];
        adv_u[iy][ix] = 0.;
        adv_v[iy][ix] = 0.;
        adv_h[iy][ix] = 0.;
        diff_u[iy][ix] = 0.;
        diff_v[iy][ix] = 0.;
        diff_h[iy][ix] = 0.;
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
        dudx_at_h(fl_x,dudx,dx,nx,ny);

        flux_y(v_b,fl_y,nx,ny);
        dvdy_at_h(fl_y,dvdy,dy,nx,ny);

        flux_x(h_b,fl_x,nx,ny);
        dhdx_at_u(fl_x,dhdx,dx,nx,ny);

        flux_y(h_b,fl_y,nx,ny);
        dhdy_at_v(fl_y,dhdy,dy,nx,ny);

        diffusion(h_b,diff_h,nx,ny);
        diffusion(u_b,diff_u,nx,ny);
        diffusion(v_b,diff_v,nx,ny);

        advection_h(u_b,v_b,h_b,adv_h,dx,dy,nx,ny);
        advection_u(u_b,v_b,adv_u,dx,dy,nx,ny);
        advection_v(u_b,v_b,adv_v,dx,dy,nx,ny);

        for(iy = 0; iy <= ny-1; iy++)
        {
          for(ix = 0; ix <= nx-1; ix++)
          { 
            tend_u[iy][ix] = g*dhdx[iy][ix] + adv_u[iy][ix] + gam*diff_u[iy][ix];
            tend_v[iy][ix] = g*dhdy[iy][ix] + adv_v[iy][ix] + gam*diff_v[iy][ix];
            tend_h[iy][ix] = (h0+h_b[iy][ix])*(dudx[iy][ix] + dvdy[iy][ix])+ adv_h[iy][ix] + gam*diff_h[iy][ix];
          }
        }

        coriolis(tend_u,tend_v,u_b,v_b,f_u,f_v,nx,ny);

        for(iy = 0; iy <= ny-1; iy++)
        { 
          for(ix = 0; ix <= nx-1; ix++)
            { 
              h_b[iy][ix]=((1.-sig_h[iy][ix]*at*0.5)*h_a[iy][ix]-at*tend_h[iy][ix])/(1.+sig_h[iy][ix]*at*0.5);
              u_b[iy][ix]=((1.-sig_u[iy][ix]*at*0.5)*u_a[iy][ix]-at*tend_u[iy][ix])/(1.+sig_u[iy][ix]*at*0.5);
              v_b[iy][ix]=((1.-sig_v[iy][ix]*at*0.5)*v_a[iy][ix]-at*tend_v[iy][ix])/(1.+sig_v[iy][ix]*at*0.5);

            }
        }
      }  /* end rk3 substerp */

      for(iy = 0; iy <= ny-1; iy++)
      { 
        for(ix = 0; ix <= nx-1; ix++)
          { 
            u_a[iy][ix] = u_b[iy][ix];
            v_a[iy][ix] = v_b[iy][ix];
            h_a[iy][ix] = h_b[iy][ix];
          }
      }

      printf("%10d,%12.6f\n", it,h_b[iys][ixs]);                              
      /*********************** END rk3 loop *************************/
      if ((it-1)%200 == 0) 
      {   
        irec++ ;
        write_to_file(u_b,v_b,h_b,irec,nx,ny);
      }
    }/* end  rk3 integration */
    
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

  free(adv_u[0]); free(adv_u);                                                 
  free(adv_v[0]); free(adv_v);                                                 
  free(adv_h[0]); free(adv_h);                                                 

  free(diff_u[0]); free(diff_u);                                                 
  free(diff_v[0]); free(diff_v);                                                 
  free(diff_h[0]); free(diff_h);                                                 

  free(tend_u[0]); free(tend_u);                                                 
  free(tend_v[0]); free(tend_v);                                                 
  free(tend_h[0]); free(tend_h);                                                 

  free(fl_x[0]); free(fl_x);                                                 
  free(fl_y[0]); free(fl_y);                                                 

  free(sig_u[0]); free(sig_u);                                                 
  free(sig_v[0]); free(sig_v);                                                 
  free(sig_h[0]); free(sig_h);                                                 
  free(f_u); free(f_v);                                                 

  return 0;
}
/******************************************************/
/******************************************************/

