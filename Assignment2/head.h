
#include <stdio.h>
#include <math.h>

typedef long double QP;
typedef double DP;
typedef float SP;

QP evaluate_f_qp(QP , QP , QP); 
DP evaluate_f_dp(DP , DP , DP); 
SP evaluate_f_sp(SP , SP , SP);
SP taylor_sp(SP , SP , SP , SP  ,  size_t) ;
