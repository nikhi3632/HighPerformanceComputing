
/*****************************************************************/
/*****************************************************************/
void advection_h(double ** u, double ** v,double ** h, double ** adv,
                 double dx, double dy, int nx, int my);
/*****************************************************************/
/*****************************************************************/
void advection_u(double ** u, double ** v, double ** adv,
                 double dx, double dy, int nx, int my);
/*****************************************************************/
/*****************************************************************/
void advection_v(double ** u, double ** v, double ** adv,
                 double dx, double dy, int nx, int my);
/*****************************************************************/
/*****************************************************************/
void diffusion(double ** rh,double ** diff,
               int nx, int my);
/*****************************************************************/
/*****************************************************************/
void coriolis(double ** tend_u,double ** tend_v,double ** u_b,
              double ** v_b, double * f_u, double * f_v, int nx, int my);
/*****************************************************************/
/*****************************************************************/