#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "flux.h"
#include "tendency.h"
#include "bound.h"
#include <omp.h>

/*******************************************/
int main()
{
  const size_t     nt     = 1000 ;
  const size_t     nrk3   = 3 ;
  const size_t     nx     = 600 ;
  const size_t     ny     = nx ;
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
  size_t it,ir,ix,iy,ixs,iys,nproc,ip;
  double xx,yy,at;
  int irec ;
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

  float *rdat_s = (float *) malloc(ny * nx* sizeof(double));
  /************************************************/
  /************************************************/

  //  Include an openmp pragma so that "nproc" can be picked up from 
  //  the setting of the environment OMP_NUM_THREADS 

  nproc = omp_get_max_threads();

  //nproc = 1; // needs to be commented after parallelezation

  size_t *j_st  = (size_t *) malloc( nproc* sizeof(size_t));
  size_t *j_en  = (size_t *) malloc( nproc* sizeof(size_t));
  size_t *j_dim = (size_t *) malloc( nproc* sizeof(size_t));

  // The following routine computes the outer loop decomposition 
  // for load-balanced parallelization

  init_myid(j_st,j_en,j_dim,ny,nproc) ;

  for(ip = 0; ip <= nproc-1; ip++)
  {
    printf(" %6zu, %6zu, %6zu, %6zu, %6zu\n", ip,j_st[ip],j_en[ip],j_dim[ip],ny);
  }

  tic = clock();
  /************************************************/
  /************************************************/
  /************************************************/

  irec = 0; ixs = 300; iys = 300;

  //  Add a pragma to parallelize the following loop
  #pragma omp parallel for private(xx,yy,ix,iy)
  for(iy = 0; iy <= ny-1; iy++) 
  { 
    for(ix = 0; ix <= nx-1; ix++)
    { 
      xx = (double) ix - (double) ixs ;
      yy = (double) iy - (double) iys ;
      rforce[iy][ix]= exp(-(xx*xx+yy*yy)/2.);
    }
  }

  //  Add a pragma to parallelize the following loop
  #pragma omp parallel for private(ix,iy)
  for(iy = 0; iy <= ny-1; iy++)
  { 
    for(ix = 0; ix <= nx-1; ix++)
    { 
      h_a[iy][ix]=0. ; h_b[iy][ix]=0. ;
      u_a[iy][ix]=0. ; u_b[iy][ix]=0. ;
      v_a[iy][ix]=0. ; v_b[iy][ix]=0. ;
      dudx[iy][ix]=0. ; dvdy[iy][ix]=0. ;
      dhdx[iy][ix]=0. ; dhdy[iy][ix]=0. ;
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
      /* Write a pragma to begin a parallel region
      Write pragmas to parallelize all the following code 
      up to the end of rk3 substep loop */
      // ip=0 ; // this statement needs to be commented and replaced by loops over ip: 1--> nproc
      #pragma omp parallel
      {
        #pragma omp for private(ip)
        for(ip = 0; ip <= nproc -1 ; ip++)
        {
          flux_x(u_b,fl_x,nx,ny,j_st[ip],j_en[ip]);
          flux_y(v_b,fl_y,nx,ny,j_st[ip],j_en[ip]);
        }

        #pragma omp for private(ip)
        for(ip = 0; ip <= nproc -1 ; ip++)
        {
          dudx_at_h(fl_x,dudx,dx,nx,ny,j_st[ip],j_en[ip]);
          dvdy_at_h(fl_y,dvdy,dy,nx,ny,j_st[ip],j_en[ip]);
        }

        #pragma omp for private(ip)
        for(ip = 0; ip <= nproc -1 ; ip++)
        {
          flux_x(h_b,fl_x,nx,ny,j_st[ip],j_en[ip]);
          flux_y(h_b,fl_y,nx,ny,j_st[ip],j_en[ip]);
        }

        #pragma omp for private(ip)
        for(ip = 0; ip <= nproc -1 ; ip++)
        {
          dhdx_at_u(fl_x,dhdx,dx,nx,ny,j_st[ip],j_en[ip]);
          dhdy_at_v(fl_y,dhdy,dy,nx,ny,j_st[ip],j_en[ip]);
        }

        //  Add a pragma to parallelize the following loop
        #pragma omp for private(ix,iy)                            
        for(iy = 0; iy <= ny-1; iy++)
        { 
          for(ix = 0; ix <= nx-1; ix++)
            { 
              h_b[iy][ix] = h_a[iy][ix] -at*h0*(dudx[iy][ix]+dvdy[iy][ix]);
              u_b[iy][ix] = u_a[iy][ix] -at*g*dhdx[iy][ix] ;
              v_b[iy][ix] = v_a[iy][ix] -at*g*dhdy[iy][ix] ;
            }
        }
      }  //   End the above parallel region
    }  /* end rk3 substerp */

    //  Add a pragma to parallelize the following loop
    #pragma omp parallel for private(ix,iy)
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

    printf("%10zu,%12.6f\n", it,h_b[iys][ixs]);

    if ((it)%100 == 0) 
    {  
      irec++ ;
      printf(" Archiving data irec = ,%10d\n", irec);
      for(iy = 0; iy <= ny-1; iy++)
      { 
        for(ix = 0; ix <= nx-1; ix++)
        { 
          rdat_s[nx*iy+ix] = (float)h_b[iy][ix];
        }
      }
      if (irec == 1) 
      { 
        file_out=fopen("shallow.bin","wb");
      } 
      else
      { 
        file_out=fopen("shallow.bin","ab");
      }
      fwrite(rdat_s, sizeof(float), nx * ny, file_out);
      fclose(file_out); 
    }  
  }  /* end  rk3 integration */
  /*********************** END rk3 loop *************************/
  toc = clock();
  // The following code computes the elapsed time. Use it to compare the speedup
  // gained by running the code on 1, 2, 4 , 6 and 8 threads 
  printf("Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC/nproc);

  free(j_st);
  free(j_en);
  free(j_dim);

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
  free(rdat_s);

  return 0;
}
/******************************************************/
/******************************************************/

