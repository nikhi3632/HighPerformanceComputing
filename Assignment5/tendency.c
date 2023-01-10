
#include <stdio.h>
/*****************************************************************/
// The functions dudx_at_h, dvdy_at_h,dhdx_at_u,dhdy_at_v                         
/*****************************************************************/

/*****************************************************************/

void dudx_at_h(double ** fl_x, double ** restrict  dudx, double dx, 
              size_t nx, size_t ny, size_t j_st ,size_t j_en)
{ 
  size_t ix , iy;

  for(iy = j_st; iy <= j_en; iy++)
  { 
    for(ix = 0; ix <= nx-2; ix++)
    { 
      dudx[iy][ix]=(fl_x[iy][ix+1]-fl_x[iy][ix])/dx;
    }
  }

  for(iy = j_st; iy <= j_en; iy++)
  { 
    dudx[iy][nx-1]=(fl_x[iy][0]-fl_x[iy][nx-1])/dx;
  }
}

/*****************************************************************/

void dvdy_at_h(double ** fl_y, double ** restrict  dvdy, double dy,
              size_t nx, size_t ny, size_t j_st ,size_t j_en)
{
  size_t ix , iy , j_start,j_end;
  j_start = j_st;
  if(j_st < 1) 
  {
    j_start = 1;
  }

  j_end = j_en; 
  if(j_en > ny-2) 
  {
    j_end = ny-2;
  }

  for(iy = j_start; iy <= j_end; iy++)
  { 
    for(ix = 0; ix <= nx-1; ix++)
    {
      dvdy[iy][ix]=(fl_y[iy+1][ix]-fl_y[iy][ix])/dy;
    }
  }
  
  if(j_en == (ny -1))
  {
    for(ix = 0; ix <= nx-1; ix++)
    {
      dvdy[ny-1][ix]=(fl_y[0][ix]-fl_y[ny-1][ix])/dy;
    }
  }
}

/*****************************************************************/

void dhdx_at_u(double ** fl_x, double ** restrict  dhdx, double dx,
              size_t nx, size_t ny, size_t j_st ,size_t j_en)
{
  size_t ix , iy;
  
  for(iy = j_st; iy <= j_en; iy++)
  { 
    for(ix = 1; ix <= nx-1;ix++)
    {
      dhdx[iy][ix]=(fl_x[iy][ix]-fl_x[iy][ix-1])/dx;
    }
  }

  for(iy = j_st; iy <= j_en; iy++)
  {
    dhdx[iy][0]=(fl_x[iy][0]-fl_x[iy][nx-1])/dx;
  }
}

/*****************************************************************/

void dhdy_at_v(double ** fl_y,  double ** restrict  dhdy, double dy,
              size_t nx, size_t ny, size_t j_st ,size_t j_en)
{
  size_t ix , iy , j_start,j_end;
  j_start = j_st;
  if(j_st < 1) 
  {
    j_start = 1;
  }

  j_end = j_en; 
  if(j_en > ny-2) 
  {
    j_end = ny-2;
  }

  for(iy = j_start; iy <= j_end; iy++)
  { 
    for(ix = 0; ix <= nx-1; ix++)
    {
      dhdy[iy][ix]=(fl_y[iy][ix]-fl_y[iy-1][ix])/dy;
    }
  }

  if(j_st == 0)
  {
    for(ix = 0; ix <= nx-1; ix++)
    {
      dhdy[0][ix]=(fl_y[0][ix]-fl_y[ny-1][ix])/dy;
    }
  }
}

/*****************************************************************/
