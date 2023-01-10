
#include <stdio.h>                                                              
#include "mpi.h"

/**************************************************************/
/**************************************************************/
void write_to_file(double **h_b, int myid,int nx,int my, int ny, 
                    int irec, const int * dims, int nproc)                          
{
     FILE  *file_out;                                                            
     int ix ,iy,ip;
     int offsets[nproc], counts[nproc];
     float dat_loc[my * nx], dat_glo[ny * nx];

     for(ip = 0; ip <= nproc-1; ip++)
     {
          counts[ip] = dims[ip]*nx;
     }

     offsets[0]=0 ;

     if (nproc > 1) 
     {
          for(ip = 1; ip <= nproc-1; ip++)
          {
             offsets[ip]=offsets[ip-1]+counts[ip-1] ;
          }
     }

     if(myid == 0) 
     {
          printf(" Archiving data irec = ,%10d\n", irec);                      
     }

     //==============    archiving h
     for(iy = 2; iy <= my+1; iy++)
     { 
          for(ix = 0; ix <= nx-1; ix++)
          { 
               dat_loc[nx*(iy-2)+ix] = (float) h_b[iy][ix];
          }
     }

     MPI_Gatherv(dat_loc,nx*my,MPI_FLOAT,dat_glo,counts,offsets,MPI_FLOAT,0,MPI_COMM_WORLD);

     if (myid == 0)
     {
          if (irec == 1)
          { 
               file_out=fopen("shallow.bin","wb");
          }                               
          else
          { 
               file_out=fopen("shallow.bin","ab");
          }                               
          fwrite(dat_glo, sizeof(float), nx * ny, file_out);                       
          fclose(file_out);                                                        
     }
}
/**************************************************************/
/**************************************************************/

