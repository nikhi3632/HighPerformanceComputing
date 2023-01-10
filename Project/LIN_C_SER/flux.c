
/*****************************************************************/
void flux_x(double ** h_b, double ** restrict fl_x, int nx, int ny)    
{ 
  int ixp1,ixm1;
  int ix,iy;

  for(iy = 0; iy <= ny-1; iy++) 
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
        ixm1=ixm1+nx;
      }
      fl_x[iy][ix]=(-h_b[iy][ixm1]+26.*h_b[iy][ix]-h_b[iy][ixp1])/24.;
    }

  }
}

///**************************************************************/
void flux_y(double ** h_b, double ** restrict fl_y, int nx ,int ny)
{
  int iyp1,iym1;
  int ix,iy;
  
  for(iy = 0; iy <= ny-1; iy++)
  {
    iyp1 = iy+1 ;
    iym1 = iy-1 ;
    if(iyp1 > ny-1) 
    {
      iyp1 = iyp1-ny;
    }
    if(iym1 < 0) 
    {
      iym1=iym1+ny;
    }
    for(ix = 0; ix <= nx-1; ix++)
    {
      fl_y[iy][ix]=(-h_b[iym1][ix]+26.*h_b[iy][ix]-h_b[iyp1][ix])/24.;
    }
  }
}
/**************************************************************/

