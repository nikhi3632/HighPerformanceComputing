
// The fuction prototypes  dudx_at_h, dvdy_at_h,dhdx_at_u,dhdy_at_v for parallelization.    

void  dudx_at_h( double ** fl_x, double ** restrict  dudx, double dx,
                 size_t nx, size_t ny, size_t j_st ,size_t j_en) ;

void  dvdy_at_h( double ** fl_y,  double ** restrict  dvdy, double dy,
                 size_t nx, size_t ny, size_t j_st ,size_t j_en) ;

void  dhdx_at_u( double ** fl_x,  double ** restrict  dhdx, double dx,
                 size_t nx, size_t ny, size_t j_st ,size_t j_en) ;

void  dhdy_at_v( double ** fl_y,  double ** restrict  dhdy, double dy,
                 size_t nx, size_t ny, size_t j_st ,size_t j_en) ;

