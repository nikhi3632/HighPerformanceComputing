
#include <iostream>
#include <valarray>
#include "flux.h"

void tendency_h(const std::valarray<double> &h_b, const std::valarray<double> &u_b,
                 std::valarray<double> &tend_h, const double dx, const double h0, const size_t nx)
{ 
  double fl_x[nx];
  size_t ix;

  flux_x(u_b,fl_x,nx);

  for(ix = 0; ix <= nx-2; ix++)
  { 
    tend_h[ix] = h0 * (fl_x[ix+1] - fl_x[ix])/dx;
  }

  tend_h[nx-1] = h0 * (fl_x[0] - fl_x[nx-1])/dx;
} 

///


//  write a procedure to compute the tendency of u
//void   tendency_u(fill in)
void tendency_u(const std::valarray<double> &h_b, const std::valarray<double> &u_b,
                 std::valarray<double> &tend_u, const double dx, const double g, const size_t nx)
{
  double fl_x[nx];
  size_t ix;

  flux_x(h_b,fl_x,nx);

  for(ix = 1; ix <= nx-1; ix++)
  { 
    tend_u[ix] = g * (fl_x[ix] - fl_x[ix-1])/dx;
  }

  tend_u[0] = g * (fl_x[0] - fl_x[nx-1])/dx;
}



