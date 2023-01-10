
/*****************************************************************/
void flux_x(double ** h_b, double ** restrict fl_x, int nx, int my)
{ 
  int ix,iy;
  int ixp1,ixm1;

  for(iy = 0; iy <= my+3; iy++)
  {
    for(ix = 0; ix <= nx-1; ix++)
    {
      ixp1 = ix+1;
      ixm1 = ix-1;
      if(ixp1 > nx-1) 
      {
        ixp1 = ixp1-nx;
      }
      if(ixm1 < 0) 
      {
        ixm1 = ixm1+nx;
      }
      fl_x[iy][ix]=(-h_b[iy][ixm1]+26.*h_b[iy][ix]-h_b[iy][ixp1])/24.;
    }
  }
} 
/**************************************************************/
void flux_y(double ** h_b, double ** restrict fl_y, int nx ,int my)
{
  int ix,iy;

  for(iy = 1; iy <= my+3-1; iy++)
  {
    for(ix = 0; ix <= nx-1; ix++)
    {
      fl_y[iy][ix]=(-h_b[iy-1][ix]+26.*h_b[iy][ix]-h_b[iy+1][ix])/24.;
    }
  }
}
/**************************************************************/
