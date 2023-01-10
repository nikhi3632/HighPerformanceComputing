#include<stdio.h>
#include<math.h>
#include "head.h"

int main(int argc, char **argv)
{
  if(argc != 2)
  {
    printf("Incorrect Usage - Correct Usage is: ./main <epsilon_value>\n");
    exit(1);
  }

  QP epsilon = atof(argv[1]);

  const size_t m = 3 ;

  QP pi_qp, pi2_qp, om0_qp, om_qp;
  QP x_qp, val_f_qp;

  DP om0_dp, om_dp;
  DP x_dp, val_f_dp, rel_f_dp;

  SP pi2_sp, om0_sp, om_sp;
  SP x_sp, val_t_sp, rel_t_sp, val_f_sp, rel_f_sp;

  pi_qp   = acosl(-1.0L);
  pi2_qp  = pi_qp * 0.5L  ;
  om0_qp  = 0.16L       ;
  x_qp    = pi_qp/8.0L  ;
  om_qp   = om0_qp + epsilon ;

  om0_dp  =  om0_qp ;
  x_dp    =  x_qp   ;
  om_dp   =  om_qp  ;

  pi2_sp  =  pi2_qp ;
  om0_sp  =  om0_qp ;
  x_sp    =  x_qp   ;
  om_sp   =  om_qp  ;

  val_f_qp = evaluate_f_qp(x_qp, om_qp, om0_qp);
  val_f_dp = evaluate_f_dp(x_dp, om_dp, om0_dp);
  val_f_sp = evaluate_f_sp(x_sp, om_sp, om0_sp);
  
  rel_f_dp = fabsl((val_f_qp - (QP)val_f_dp)/val_f_qp);
  rel_f_sp = fabsl((val_f_qp - (QP)val_f_sp)/val_f_qp);

  val_t_sp = taylor_sp(x_sp, om_sp, om0_sp, pi2_sp,m);

  rel_t_sp = fabsl((val_f_qp - (QP)val_t_sp)/val_f_qp);

  printf("    quad prec,    double prec,    single prec,    rel double prec,    rel single prec,     taylor: single prec,    taylor: rel single prec\n" );                
            
  printf("%18.8Le%18.8e%18.8e%18.8e%18.8e%18.8e%18.8e\n",
        val_f_qp, val_f_dp, val_f_sp, rel_f_dp, rel_f_sp, val_t_sp, rel_t_sp);
  printf("successful execution \n")  ;

  return 0;
}
