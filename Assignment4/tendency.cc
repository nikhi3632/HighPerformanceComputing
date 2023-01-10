
#include <iostream>
#include <vector>

/*****************************************************************/
void dudx_at_h(const std::vector<std::vector<double> >&fl_x,
                      std::vector<std::vector<double> >&dudx,
                      double dx, size_t nx, size_t ny)
{ 
  size_t ix , iy;

  for(iy = 0; iy <= ny-1; iy++)
  { 
    for(ix = 0; ix <= nx-2; ix++)
    { 
      dudx[iy][ix]=(fl_x[iy][ix+1]-fl_x[iy][ix])/dx;
    }
  }

  for(iy = 0; iy <= ny-1; iy++)
  { 
    dudx[iy][nx-1]=(fl_x[iy][0]-fl_x[iy][nx-1])/dx;
  }
} 
/*****************************************************************/
void dvdy_at_h(const std::vector<std::vector<double> >&fl_y,
                      std::vector < std::vector<double> >&dvdy,
                      double dy, size_t nx, size_t ny)
{
  size_t ix , iy;

  for(iy = 0; iy <= ny-2; iy++)
  { 
    for(ix = 0; ix <= nx-1; ix++)
    {
      dvdy[iy][ix]=(fl_y[iy+1][ix]-fl_y[iy][ix])/dy;
    }
  }

  for(ix = 0; ix <= nx-1; ix++)
  {
    dvdy[ny-1][ix]=(fl_y[0][ix]-fl_y[ny-1][ix])/dy;
  }
}
/*****************************************************************/
void dhdx_at_u(const std::vector<std::vector<double> >&fl_x,
                      std::vector<std::vector<double> >&dhdx,
                      double dx, size_t nx, size_t ny)
{
  size_t ix , iy;

  for(iy = 0; iy <= ny-1; iy++)
  { 
    for(ix = 1; ix <= nx-1; ix++)
    { 
      dhdx[iy][ix]=(fl_x[iy][ix]-fl_x[iy][ix-1])/dx;
    }
  }

  for(iy = 0; iy <= ny-1; iy++)
  { 
    dhdx[iy][0]=(fl_x[iy][0]-fl_x[iy][nx-1])/dx;
  }
}
/*****************************************************************/
void dhdy_at_v(const std::vector<std::vector<double> >&fl_y,
                      std::vector<std::vector<double> >&dhdy,
                      double dy, size_t nx, size_t ny)
{
  size_t ix , iy;

  for(iy = 1; iy <= ny-1; iy++)
  { 
    for(ix = 0; ix <= nx-1; ix++)
    {
      dhdy[iy][ix]=(fl_y[iy][ix]-fl_y[iy-1][ix])/dy;
    }
  }

  for(ix = 0; ix <= nx-1; ix++)
  {
    dhdy[0][ix]=(fl_y[0][ix]-fl_y[ny-1][ix])/dy;
  }
}
/*****************************************************************/
