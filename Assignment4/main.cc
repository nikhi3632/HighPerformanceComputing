
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "flux.h"
#include "tendency.h"

using namespace std;

/*******************************************/
int main()
{
  const size_t     nt     = 1000 ;  // total number of time steps
  const size_t     nrk3   = 3 ;
  const size_t     nx     = 600 ;   // x-number of grid points
  const size_t     ny     = nx ;    // y-number of grid points
  const double     g      = 9.81 ;
  const double     pi     = acos(-1.0) ;
  const double     c      = 30.0 ;
  const double     h0     = c*c/g ;
  const double     cour   = 0.40 ;
  const double     dx     = 100.0 ;  // x-grid spacing
  const double     dy     = dx ;     // y-grid spacing
  const double     dt     = cour*dx/c/sqrt(2.0) ;
  const double     omega  = 4.e-3 ;   // frequency of the forcing
  const double     amp    = 10.0 ;    // amplitude of the forcing
/************************************************/
  FILE* file_out;
  size_t it,ir,ix,iy,ixs,iys;
  double xx,yy,at;
  int irec = 0;
  cout.precision(6);
/************************************************/
// allocate the fields needed for computation 
/************************************************/
  vector<vector<double> > h_a(ny, vector<double>(nx));
  vector<vector<double> > h_b(ny, vector<double>(nx));

  vector<vector<double> > u_a(ny, vector<double>(nx));
  vector<vector<double> > u_b(ny, vector<double>(nx));

  vector<vector<double> > v_a(ny, vector<double>(nx));
  vector<vector<double> > v_b(ny, vector<double>(nx));

  vector<vector<double> > dudx(ny, vector<double>(nx));
  vector<vector<double> > dvdy(ny, vector<double>(nx));

  vector<vector<double> > dhdx(ny, vector<double>(nx));
  vector<vector<double> > dhdy(ny, vector<double>(nx));

  vector<vector<double> > fl_x(ny, vector<double>(nx));
  vector<vector<double> > fl_y(ny, vector<double>(nx));

  vector<vector<double> > rforce(ny, vector<double>(nx));

  float * rdat_s = new float[nx*ny];
/************************************************/
/************************************************/
  ixs = 300; iys = 300;   //x,y center of the forcing

  for(iy = 0; iy <= ny-1; iy++) 
  { 
    for(ix = 0; ix <= nx-1; ix++)
    { 
      xx = (double) ix - (double) ixs;
      yy = (double) iy - (double) iys;

      rforce[iy][ix]= exp(-(xx*xx+yy*yy)/2.);
    }
  }
  // initialaze h, u and v fields
  for(iy = 0; iy <= ny-1; iy++)
  { 
    for(ix = 0; ix <= nx-1; ix++)
    { 
      h_a[iy][ix]=0.; h_b[iy][ix]=0.;
      u_a[iy][ix]=0.; u_b[iy][ix]=0.;
      v_a[iy][ix]=0.; v_b[iy][ix]=0.;

      dudx[iy][ix]=0.; dvdy[iy][ix]=0.;
      dhdx[iy][ix]=0.; dhdy[iy][ix]=0.;
    }
  }
  /*********************** Begin rk3 loop *************************/
  for(it = 1; it <= nt; it++) /* start rk3 integration */
  {
    for(ir = 1; ir <= nrk3; ir++)  /* start rk3 substep */
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
      dudx_at_h(fl_x,dudx,dx,nx,ny) ;  // compute du/dx at h points

      flux_y(v_b,fl_y,nx,ny);
      dvdy_at_h(fl_y,dvdy,dy,nx,ny) ;  // compute dv/dy at h points

      flux_x(h_b,fl_x,nx,ny);
      dhdx_at_u(fl_x,dhdx,dx,nx,ny) ;  // compute dh/dx at u points

      flux_y(h_b,fl_y,nx,ny);
      dhdy_at_v(fl_y,dhdy,dy,nx,ny) ;  // compute dh/dy at v points

      //compute u_b, v_b, h_b
      for(iy = 0; iy <= ny-1; iy++)
      { 
        for(ix = 0; ix <= nx-1; ix++)
        { 
          u_b[iy][ix] = u_a[iy][ix] - at*g*dhdx[iy][ix];
          v_b[iy][ix] = v_a[iy][ix] - at*g*dhdy[iy][ix];
          h_b[iy][ix] = h_a[iy][ix] - at*h0*(dudx[iy][ix]+dvdy[iy][ix]);
        }
      }
    } 
    /* end rk3 substep */

    // the following add the forcing to h and overwrite the fields at (n) by those computed at (n+1)
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

    cout << fixed << setw(12) << it  << fixed << setw(12) << h_b[iys][ixs] << endl;
    /*********************** END rk3 loop *************************/
    // the following code  opens a binary file shallow.bin to write  h for each 100 time steps

    if ((it)%100 == 0) 
    {   
      irec++;
      cout << fixed << setw(12) << "Archiving data irec ="  << 
              fixed << setw(12) <<  irec << endl;
      for(iy = 0; iy <= ny-1; iy++)
      { 
        for(ix = 0; ix <= nx-1; ix++)
        { 
          rdat_s[nx*iy+ix] = (float) h_b[iy][ix];
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
  delete [] rdat_s;
  return 0;
}
/******************************************************/
/******************************************************/

