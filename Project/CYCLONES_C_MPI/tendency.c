/*****************************************************************/
void dudx_at_h(double ** fl_x, double ** dudx,                                
                 double dx, int nx, int my)                                    
{ 
  int ix,iy;
  int ixp1;

  for(iy = 0; iy <= my+3; iy++)
  { 
    for(ix = 0; ix <= nx-1; ix++)
    { 
      ixp1 = ix+1 ;
      if(ixp1 > nx-1) 
      {
        ixp1 = ixp1-nx;
      }
      dudx[iy][ix]=(fl_x[iy][ixp1]-fl_x[iy][ix])/dx;
    }
  }
} 
/*****************************************************************/
void dvdy_at_h(double ** fl_y, double ** dvdy,                                
                 double dy, int nx, int my)                                    
{
  int ix,iy;

  for(iy = 0; iy <= my+3-1; iy++)
  {
    for(ix = 0; ix <= nx-1; ix++)
    {
      dvdy[iy][ix]=(fl_y[iy+1][ix]-fl_y[iy][ix])/dy;
    }
  }
}
/*****************************************************************/
void  dhdx_at_u( double ** fl_x, double ** dhdx,                                
                 double dx, int nx, int my)                                    
{
  int ix,iy;
  int ixm1;

  for(iy = 0; iy <= my+3; iy++)
  { 
    for(ix = 0; ix <= nx-1; ix++)
    {
      ixm1 = ix-1;
      if(ixm1 < 0) 
      {
        ixm1 = ixm1+nx;
      }
      dhdx[iy][ix]=(fl_x[iy][ix]-fl_x[iy][ixm1])/dx;
    }
  }
}
/*****************************************************************/
void dhdy_at_v( double ** fl_y, double ** dhdy,                                
                 double dy, int nx, int my)                                    
{
  int ix,iy;

  for(iy = 1; iy <= my+3; iy++)
  {
    for(ix = 0; ix <= nx-1; ix++)
    {
      dhdy[iy][ix]=(fl_y[iy][ix]-fl_y[iy-1][ix])/dy;
    }
  }
}
/*****************************************************************/
