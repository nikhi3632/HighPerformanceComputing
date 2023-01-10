
/*****************************************************************/
void dudx_at_h(double ** fl_x, double ** dudx, 
               double dx, int nx, int my);
/*****************************************************************/
void dvdy_at_h(double ** fl_y, double ** dvdy,
               double dy, int nx, int my);
/*****************************************************************/
void dhdx_at_u(double ** fl_x, double ** dhdx,
               double dx, int nx, int my);
/*****************************************************************/
void dhdy_at_v(double ** fl_y, double ** dhdy,
               double dy, int nx, int my);
/*****************************************************************/
