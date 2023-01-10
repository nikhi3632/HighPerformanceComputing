
#include <iostream>
#include <valarray>

void flux_x(const std::valarray<double> &h_b, double * fl_x, const size_t nx) 
 { 
    size_t ix;

    for(ix = 1; ix <= nx-2; ix++) 
    {
        fl_x[ix] = (-h_b[ix-1]+26.*h_b[ix]-h_b[ix+1])/24.;
    }

    fl_x[0] = (-h_b[nx-1]+26.*h_b[0]-h_b[1])/24.;
    fl_x[nx-1] = (-h_b[nx-2]+26.*h_b[nx-1]-h_b[0])/24.;
 } 

