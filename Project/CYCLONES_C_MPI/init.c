#include <math.h>                                                               
#include "input_output.h"                                                       
/*****************************************************************/
void absorb(double ** sig_u,double ** sig_v,double ** sig_h,                   
            double dt, int j_st, int j_en, int nx, int ny, int my)                                        
{ 
  const int np  = 10;
  double sig_p[np];
  double sigx[nx], sigy[ny];
  double val;
  int i;
  int ix , iy;
  int iyy;
  int ixm1, iym1 ;

  for(i = 0; i <= np-1; i++) 
  { 
    sig_p[i] = 0.8*pow((double)(np-i)/(double)(np),3.);
  }

  for(ix = 0; ix <= nx-1; ix++) 
  { 
    sigx[ix] = 0.;
  }

  for(ix = 0; ix <= np-1; ix++) 
  { 
    sigx[ix] = sig_p[ix];
  }

  for(ix = nx-np+1; ix <= nx-1; ix++) 
  {  
    sigx[ix] = sigx[nx-ix]; 
  }

  for(iy = 0; iy <= ny-1; iy++) 
  { 
    sigy[iy] = 0.;
  }

  for(iy = 0; iy <= np-1; iy++) 
  { 
    sigy[iy] = sig_p[iy]; 
  }

  for(iy = ny-np+1; iy <= ny-1; iy++) 
  {  
    sigy[iy] = sigy[ny-iy];
  }
  /*****************************************************************/

  for(iy = j_st; iy <= j_en; iy++) 
  {
    iyy = iy - j_st + 2;
    for(ix = 0; ix <= nx-1; ix++) 
    {
      val = sigx[ix];
      if (sigy[iy] > sigx[ix]) 
      {
        val = sigy[iy];
      }
      sig_h[iyy][ix]= val/dt;
    }
  }

  iy = j_st - 1;
  iyy = iy - j_st + 2;
  if(j_st == 0)
  {
    iy = ny - 1;
  }
  for(ix = 0; ix <= nx -1; ix++)
  {
    val = sigx[ix];
    if(sigy[iy] > sigx[ix])
    {
      val = sigy[iy];
    }
    sig_h[iyy][ix] = val/dt;
  }

  for(iy = 2; iy <= my+1; iy++) 
  {
    for(ix = 0; ix <= nx-1; ix++) 
    {
      ixm1=ix-1; 
      if(ixm1 < 0) 
      {
        ixm1 = ixm1+nx;
      }
      sig_u[iy][ix]= (sig_h[iy][ix] + sig_h[iy][ixm1])*0.5;
    }
  }

  for(iy = 2; iy <= my+1; iy++) 
  {
    iym1 = iy-1; 
    if(iym1 < 0) 
    {
      iym1 = iym1+ny;
    }
    for(ix = 0; ix <= nx-1; ix++) 
    {
      sig_v[iy][ix]= (sig_h[iy][ix] + sig_h[iym1][ix])*0.5;
    }
  }
} 
/**************************************************************/
void etat0(double ** u_a,double ** v_a,double ** h_a,                         
           double * f_u, double * f_v,                                        
           double dy, int j_st, int j_en, int myid,
           int nx, int my, int ny, int * dims, int nproc)                                        
{ 
  const double  radius   = 6357.e3;
  const double  omega    = 7.2921e-5;
  const double  beta     = 2.0*omega/radius;
  const double  box      = 7500.e3;

  double fcor[ny], ry[ny];
  int iy,iyp,iyy;

  for(iy = 0; iy <= ny-1; iy++) 
  {
    ry[iy] = -box + dy * (double) (iy);
    fcor[iy] = beta*ry[iy]; 
  }
  /********************************************************/
  for(iy = j_st; iy <= j_en; iy++)
  {
    iyy = iy - j_st + 2;
    iyp = iy + 1;
    if(iy == ny-1) 
    {
      iyp = 0;
    }
    f_v[iyy]= fcor[iy];
    f_u[iyy]= 0.5*(fcor[iyp]+fcor[iy]) ;
  }
  read_from_file(u_a,v_a,h_a,myid,nx,my,ny,dims,nproc);
} 
/**************************************************************/
