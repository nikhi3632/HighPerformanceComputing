
void tendency_h(const std::valarray<double> &h_b, const  std::valarray<double> &u_b,
                 std::valarray<double>   &tend_h, const double dx, const double h0, const size_t nx);
///
//write a function prototype for tendency_u
//void tendency_u (fill in)

void tendency_u(const std::valarray<double> &h_b, const  std::valarray<double> &u_b,
                 std::valarray<double>   &tend_u, const double dx, const double g, const size_t nx);
