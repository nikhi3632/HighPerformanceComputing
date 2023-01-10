/*****************************************************************/
void absorb(double ** sig_u,double ** sig_v,double ** sig_h,                   
            double dt, int j_st, int j_en, int nx, int ny, int my);

void etat0(double ** u_a,double ** v_a,double ** h_a,                         
           double * f_u, double * f_v,                                        
           double dy, int j_st, int j_en, int myid,
           int nx, int my, int ny, int * dims, int nproc);
/**************************************************************/
/**************************************************************/
