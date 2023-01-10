/*****************************************************************/
/*****************************************************************/
void advection_h(double ** u, double ** v,double ** h, double ** adv,          
                 double dx, double dy, int nx, int my)                        
{ 
  double uu,vv;
  int ixp1,ixp2,ixm1,ixm2;
  int iyp1,iyp2,iym1,iym2;
  int ix,iy;

  for(iy = 2; iy <= my+1; iy++)
  { 
    for(ix = 0; ix <= nx-1; ix++)
    { 
      ixp1 = ix+1;
      ixm1 = ix-1;
      ixp2 = ix+2;
      ixm2 = ix-2;
      if(ixp1 > nx-1) 
      {
        ixp1 = ixp1-nx;
      }
      if(ixm1 < 0) 
      {
        ixm1 = ixm1+nx;
      }
      if(ixp2 > nx-1) 
      {
        ixp2 = ixp2-nx;
      }
      if(ixm2 < 0) 
      {
        ixm2 = ixm2+nx;
      }
      uu = 0.5*(u[iy][ixp1]+u[iy][ix]) ;
      adv[iy][ix] = uu*(8.*(h[iy][ixp1]-h[iy][ixm1])-(h[iy][ixp2]-h[iy][ixm2]))/(12.*dx) ;
    }
  }

  for(iy = 2; iy <= my+1; iy++)
  { 
    iyp1 = iy+1;
    iym1 = iy-1;
    iyp2 = iy+2;
    iym2 = iy-2;

    for(ix = 0; ix <= nx-1; ix++)
    { 
      vv=0.5*(v[iy][ix]+v[iyp1][ix]);
      adv[iy][ix]= adv[iy][ix]+vv*(8.*(h[iyp1][ix]-h[iym1][ix])-(h[iyp2][ix]-h[iym2][ix]))/(12.*dy);
    }
  } 
} 
/*****************************************************************/
/*****************************************************************/
void advection_u(double ** u, double ** v, double ** adv,                      
                 double dx, double dy, int nx, int my)                        
{ 
  double uu,vv ;
  int ixp1,ixp2,ixm1,ixm2 ;
  int iyp1,iyp2,iym1,iym2 ;
  int ix,iy;

  for(iy = 2; iy <= my+1; iy++)
  { 
    for(ix = 0; ix <= nx-1; ix++)
    { 
      ixp1 = ix+1;
      ixm1 = ix-1;
      ixp2 = ix+2;
      ixm2 = ix-2;
      if(ixp1 > nx-1) 
      {
        ixp1 = ixp1-nx;
      }
      if(ixm1 < 0) 
      {
        ixm1 = ixm1+nx;
      }
      if(ixp2 > nx-1) 
      {
        ixp2 = ixp2-nx;
      }
      if(ixm2 < 0) 
      {
        ixm2 = ixm2+nx;
      }
      uu = u[iy][ix];
      adv[iy][ix]= uu*(8.*(u[iy][ixp1]-u[iy][ixm1])-(u[iy][ixp2]-u[iy][ixm2]))/(12.*dx) ;
    }
  }

  for(iy = 2; iy <= my+1; iy++)
  { 
    iyp1 = iy+1;
    iym1 = iy-1;
    iyp2 = iy+2;
    iym2 = iy-2;

    for(ix = 0; ix <= nx-1;ix ++)
    { 
      ixm1 = ix-1 ;
      if(ixm1 < 0) 
      {
        ixm1 = ixm1+nx;
      }
      vv = 0.25*(v[iy][ix] + v[iyp1][ix]+v[iy][ixm1]+v[iyp1][ixm1]);
      adv[iy][ix]= adv[iy][ix] + vv*(8.*(u[iyp1][ix]-u[iym1][ix])-(u[iyp2][ix]-u[iym2][ix]))/(12.*dy) ;
    }
  }
} 
/*****************************************************************/
/*****************************************************************/
void advection_v(double ** u, double ** v, double ** adv ,                      
                 double dx, double dy, int nx, int my)                        
{ 
  double uu,vv ;
  int ixp1,ixp2,ixm1,ixm2 ;
  int iyp1,iyp2,iym1,iym2 ;
  int ix,iy;

  for(iy = 2; iy <= my+1; iy++)
  { 
    iym1 = iy-1;
    for(ix = 0; ix <= nx-1; ix++)
    { 
      ixp1 = ix+1 ;
      ixm1 = ix-1 ;
      ixp2 = ix+2 ;
      ixm2 = ix-2 ;
      if(ixp1 > nx-1) 
      {
        ixp1 = ixp1-nx;
      }
      if(ixm1 < 0) 
      {
        ixm1 = ixm1+nx;
      }
      if(ixp2 > nx-1) 
      {
        ixp2 = ixp2-nx;
      }
      if(ixm2 < 0) 
      {
        ixm2 = ixm2+nx;
      }
      uu=0.25*(u[iy][ix]+u[iym1][ix]+u[iy][ixp1]+u[iym1][ixp1]) ;
      adv[iy][ix]= uu*(8.*(v[iy][ixp1]-v[iy][ixm1])-(v[iy][ixp2]-v[iy][ixm2]))/(12.*dx) ;
    }
  }

  for(iy = 2; iy <= my+1; iy++)
  { 
    iyp1 = iy+1;
    iym1 = iy-1;
    iyp2 = iy+2;
    iym2 = iy-2;

    for(ix = 0; ix <= nx-1; ix++)
    { 
      vv = v[iy][ix];
      adv[iy][ix] = adv[iy][ix] + vv*(8.*(v[iyp1][ix]-v[iym1][ix])-(v[iyp2][ix]-v[iym2][ix]))/(12.*dy);
    }
  }
} 
/*****************************************************************/
/*****************************************************************/
void diffusion(double ** rh,double ** diff, int nx, int my)                                                
{ 
  int ixp1,ixp2,ixm1,ixm2;
  int iyp1,iyp2,iym1,iym2;
  int ix,iy;

  for(iy = 2; iy <= my+1; iy++)
  { 
    for(ix = 0; ix <= nx-1; ix++)
    { 
      ixp1 = ix+1;
      ixm1 = ix-1;
      ixp2 = ix+2;
      ixm2 = ix-2;
      if(ixp1 > nx-1) 
      {
        ixp1 = ixp1-nx;
      }
      if(ixm1 < 0) 
      {
        ixm1 = ixm1+nx;
      }
      if(ixp2 > nx-1) 
      {
        ixp2 = ixp2-nx;
      }
      if(ixm2 < 0) 
      {
        ixm2 = ixm2+nx;
      }
      diff[iy][ix] = rh[iy][ixm2]-4.*rh[iy][ixm1]+6.*rh[iy][ix]-4.*rh[iy][ixp1]+rh[iy][ixp2] ;
    }
  }

  for(iy = 2; iy <= my+1; iy++)
  { 
    iyp1 = iy+1;
    iym1 = iy-1;
    iyp2 = iy+2;
    iym2 = iy-2;

    for(ix = 0; ix <= nx-1; ix++)
    { 
      diff[iy][ix]= diff[iy][ix]+rh[iym2][ix]-4.*rh[iym1][ix]+6.*rh[iy][ix]-4.*rh[iyp1][ix]+rh[iyp2][ix] ;
    }
  } 
} 
/*****************************************************************/
/*****************************************************************/
void coriolis(double ** tend_u,double ** tend_v,double ** u_b,                 
              double ** v_b, double * f_u, double * f_v,                       
              int nx, int my)                                                 
{ 
  int ixp1,ixm1;
  int iyp1,iym1;
  int ix,iy;

  for(iy = 2; iy <= my+1; iy++)
  { 
    iyp1 = iy+1;
    iym1 = iy-1;

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
      tend_u[iy][ix] = tend_u[iy][ix]-f_u[iy]*(v_b[iy][ix]+v_b[iyp1][ix]+v_b[iy][ixm1]+v_b[iyp1][ixm1])*0.25;
      tend_v[iy][ix] = tend_v[iy][ix]+f_v[iy]*(u_b[iy][ix]+u_b[iym1][ix]+u_b[iy][ixp1]+u_b[iym1][ixp1])*0.25;
    }
  } 
} 
/*****************************************************************/
/*****************************************************************/
