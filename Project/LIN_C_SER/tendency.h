void dudx_at_h(double ** fl_x, double ** restrict dudx, double dx, int nx, int ny);
void dvdy_at_h(double ** fl_y, double ** restrict dvdy, double dy, int nx, int ny);
void dhdx_at_u(double ** fl_x, double ** restrict dhdx, double dx, int nx, int ny);
void dhdy_at_v(double ** fl_y, double ** restrict dhdy, double dy, int nx, int ny);
