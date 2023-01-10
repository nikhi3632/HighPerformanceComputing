#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <valarray>
#include "tendency.h"

using namespace std;

/*******************************************/
int main()
{
  const size_t     nrk3   =  3;
  const size_t     nx     =  160;
  const double     g      =  9.81;
  const double     c      =  30.;
  const double     h0     =  c*c/g;
  const double     cour   =  0.5;
  const double     dx     =  100.;
  const double     dd     =  ((double)nx)*dx;
  const double     dt     =  cour*dx/c;
  const double     dl     =  1000.;
/************************************************/
  valarray<double> u_a(nx),u_b(nx),h_a(nx),h_b(nx); //  _a =n, _b=n+1
  valarray<double> tend_u(nx), tend_h(nx); // tendency for u and h
  valarray<double> h_exa(nx), rx(nx);  // rx=x

/************************************************/
  ofstream file_out;
  file_out.precision(6);
  cout.precision(6);
/************************************************/
  size_t nt, it, ix, ir, ncycle;
  double at, time, xx;
/************************************************/
  // nt = (size_t)(2.10*dd/(c*dt));   /* total number of time steps */ 
  nt = (size_t)((2.10)*(nx/cour));    /* equivalent to the above nt */
  for(ix = 0; ix <= nx-1; ix++) 
  {  
    rx[ix]= -8000. + 0.5 * dx + dx * (double)ix;     /* position x(i) at h points */
    h_a[ix]= exp(-(rx[ix]/dl)*(rx[ix]/dl)) ;
  }  /* initialize h(t=0) */

// the following code computes  the conditions for velocity
  for(ix = 1; ix <= nx-1; ix++) 
  {
    u_a[ix] = 0.5*g/c*(h_a[ix]+h_a[ix-1]);
  }
  u_a[0] = 0.5*g/c*(h_a[0]+h_a[nx-1]);
//   the following code initializes h_b and u_b needed for the first rk3 step
  for(ix = 1; ix <= nx-1; ix++) 
  { 
    h_b[ix] = h_a[ix]; 
    u_b[ix] = u_a[ix];
  }
////
/*********************** Begin rk3 integration *************************/
  for(it = 1; it <= nt; it++) /* begin rk3 loop */
  {
    for(ir = 1; ir <= nrk3; ir++)  /* begin rk3 substep loop*/
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

      tendency_h(h_b, u_b, tend_h, dx, h0, nx);   /* compute h tendency dh/dt */
      tendency_u(h_b, u_b, tend_u, dx, g,  nx);   /* compute u tendency du/dt */

      for(ix = 0; ix <= nx-1; ix++) 
      {  /* update u(s+1), h(s+1) */
        u_b[ix]= u_a[ix]-(at*tend_u[ix]); 
        h_b[ix]= h_a[ix]-(at*tend_h[ix]);
      }

    }  /* end rk3 substep loop */
    for(ix = 0; ix <= nx-1; ix++) 
    {
      h_a[ix] = h_b[ix];   // ovewrite h(n) by h(n+1)
      u_a[ix] = u_b[ix];   // ovewrite u(n) by u(n+1)
    }  
    cout << fixed << setw(12) << it  << fixed << setw(12) << u_a[nx/2] << fixed << setw(12) << h_a[nx/2] << endl;
  }  /* end  rk3 loop */
/*********************** END rk3 integration *************************/

  // write code to compute the exact solution of h  : h_exa at t =nt*dt 
  // fill in
  time = dt*nt;
  ncycle = (int)(c*time/dd);
  for(ix = 0; ix <= nx-1; ix++)
  {
    xx = rx[ix] - c*time + dd*double(ncycle);
    if (xx < -0.5*dd)
    {
      xx = xx + dd;
    }
    h_exa[ix] = exp(-(xx/dl)*(xx/dl));
  }
  // the following code  opens the file shallow.dat to write rx in km, h_a, h_exa
  file_out.open("shallow.dat");  
  for(ix = 0; ix <= nx-1; ix++) 
  {
    file_out <<  fixed << setw(12) << rx[ix]/1000. << fixed << setw(12) << h_a[ix] << fixed << setw(12) << h_exa[ix] << endl;
  } 

  file_out.close();

  return 0;
};
/******************************************************/
/******************************************************/
