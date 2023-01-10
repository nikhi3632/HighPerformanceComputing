#include "head.h"

QP evaluate_f_qp(QP t,QP om,QP om0) 
{
  QP v;
  v = (cosl(om * t)- cosl(om0 * t))/((om * om) - (om0 * om0));
  return v;
}

DP evaluate_f_dp(DP t,DP om,DP om0) 
{
  DP v;
  v = (cosl(om * t)- cosl(om0 * t))/((om * om) - (om0 * om0));
  return v;
}

SP evaluate_f_sp(SP t,SP om,SP om0) 
{
  SP v;
  v = (cosl(om * t)- cosl(om0 * t))/((om * om) - (om0 * om0));
  return v;
}

SP taylor_sp(SP t, SP om, SP om0, SP pi2 ,  size_t nd)
{
  SP v, a, a0, fac, tay;
  a   = (om * t); 
  a0  = (om0 * t);
  fac = 1.0f;
  tay = 0.0f;
  size_t j;
  for(j = 1; j <= nd; j++) 
  {
    fac = (fac * j);
    tay = tay + t * t * powf((a - a0), (SP)(j - 1)) * cosf(a0 + pi2*j)/fac/(a + a0);
  }
  v = tay;
  return v;
}
