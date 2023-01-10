
#include <stdio.h>

/**************************************************************/
// The functions flux_x, flux_y
/**************************************************************/
void flux_x(double ** h_b, double ** restrict fl_x,
            size_t nx, size_t ny, size_t j_st ,size_t j_en)
{ 
  size_t ix , iy;

  for(iy = j_st; iy <= j_en; iy++) 
  {
    for(ix = 1; ix <= nx-2; ix++) 
    {
      fl_x[iy][ix]=(-h_b[iy][ix-1]+26.*h_b[iy][ix]-h_b[iy][ix+1])/24.;
    }
  }

  for(iy = j_st; iy <= j_en; iy++) 
  {
    fl_x[iy][0]=(-h_b[iy][nx-1]+26.*h_b[iy][0]-h_b[iy][1])/24.;
    fl_x[iy][nx-1]=(-h_b[iy][nx-2]+26.*h_b[iy][nx-1]-h_b[iy][0])/24.;
  }
}

///**************************************************************/

void flux_y(double ** h_b, double ** restrict fl_y, 
            size_t nx, size_t ny, size_t j_st ,size_t j_en )   
{
  size_t ix , iy, j_start,j_end;

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
      fl_y[iy][ix]=(-h_b[iy-1][ix]+26.*h_b[iy][ix]-h_b[iy+1][ix])/24.;
    }
  }

  if(j_st == 0) 
  {
    for(ix = 0; ix <= nx-1; ix++)
    {
      fl_y[0][ix]=(-h_b[ny-1][ix]+26.*h_b[0][ix]-h_b[1][ix])/24.;
    }
  }

  if(j_en == ny-1) 
  {
    for(ix = 0; ix <= nx-1; ix++)
    {
      fl_y[ny-1][ix]=(-h_b[ny-2][ix]+26.*h_b[ny-1][ix]-h_b[0][ix])/24.;
    }
  }
}
/**************************************************************/

