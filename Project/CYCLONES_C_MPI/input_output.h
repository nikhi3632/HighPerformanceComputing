
void read_from_file(double ** u_a, double ** v_a, double ** h_a,               
                    int myid, int nx, int my, int ny,
                    int * dims, int nproc);                   
/**************************************************************/
void write_to_file(double ** u_b, double ** v_b, double ** h_b,
                   int myid, int nx, int my, int ny,              
                   int irec, const int * dims, int nproc);                   
/**************************************************************/
