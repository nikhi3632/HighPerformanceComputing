/*****************************************************************/
void dudx_at_h(const std::vector<std::vector<double> >&fl_x,
                     std::vector < std::vector<double> >&dudx,
                     double dx, size_t nx, size_t ny);
/*****************************************************************/
void dvdy_at_h(const std::vector<std::vector<double> >&fl_y,
                     std::vector<std::vector<double> >&dvdy,
                     double dy, size_t nx, size_t ny);
/*****************************************************************/
//write a function prototypes  for dhdx_at_u and dhdy_at_v
void dhdx_at_u(const std::vector<std::vector<double> >&fl_x,
                     std::vector < std::vector<double> >&dhdx,
                     double dx, size_t nx, size_t ny);
/*****************************************************************/
void dhdy_at_v(const std::vector<std::vector<double> >&fl_y,
                     std::vector < std::vector<double> >&dhdy,
                     double dy, size_t nx, size_t ny);
/*****************************************************************/

