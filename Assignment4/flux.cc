
#include <iostream>
#include <vector>
/*****************************************************************/
void flux_x(const std::vector<std::vector<double> > &h_b,
                  std::vector < std::vector<double> > &fl_x,
                  size_t nx, size_t ny)
{ 
  size_t ix , iy;

  for(iy = 0; iy <= ny-1; iy++) 
  {
    for(ix = 1; ix <= nx-2; ix++) 
    {
      fl_x[iy][ix]=(-h_b[iy][ix-1]+26.*h_b[iy][ix]-h_b[iy][ix+1])/24.;
    }
  }

  for(iy = 0; iy <= ny-1; iy++) 
  {
      fl_x[iy][0]=(-h_b[iy][nx-1]+26.*h_b[iy][0]-h_b[iy][1])/24.;
      fl_x[iy][nx-1]=(-h_b[iy][nx-2]+26.*h_b[iy][nx-1]-h_b[iy][0])/24.;
  }
}
///**************************************************************/
void flux_y(const std::vector<std::vector<double> > &h_b,
                  std::vector < std::vector<double> > &fl_y,
                  size_t nx, size_t ny)
{
  size_t ix , iy;

  for(iy = 1; iy <= ny-2; iy++)
  {
    for(ix = 0; ix <= nx-1; ix++)
    {
      fl_y[iy][ix]=(-h_b[iy-1][ix]+26.*h_b[iy][ix]-h_b[iy+1][ix])/24.;
    }
  }

  for(ix = 0; ix <= nx-1; ix++)
  {
    fl_y[0][ix]=(-h_b[ny-1][ix]+26.*h_b[0][ix]-h_b[1][ix])/24.;
    fl_y[ny-1][ix]=(-h_b[ny-2][ix]+26.*h_b[ny-1][ix]-h_b[0][ix])/24.;
  }
}
/**************************************************************/

