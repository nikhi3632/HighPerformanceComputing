
#include <stdlib.h>                                                             
#include <stdio.h>
#include "mpi.h"
                                                            
/**************************************************************/
/**************************************************************/
void read_from_file(double ** u_a, double ** v_a, double ** h_a,               
                    int myid, int nx, int my, int ny,
                    int * dims, int nproc)                                  
{
     FILE *file_in;
     int ix,iy,ip;
     int offsets[nproc], counts[nproc];
     float dat_loc[my * nx], dat_glo[ny * nx];

     for(ip = 0; ip <= nproc-1; ip++)
     {
          counts[ip] = dims[ip]*nx;
     }
     
     offsets[0] = 0;

     if(nproc > 1)
     {
          for(ip = 1; ip <= nproc -1 ; ip++)
          {
               offsets[ip] = offsets[ip-1] + counts[ip-1];
          }
     }

     //==============    reading h
     if(myid == 0)
     {
          printf("Reading initial data for h \n");                      
          file_in=fopen("shallow_init_h.bin","rb");                               
          fread(dat_glo, sizeof(float), nx * ny, file_in);                       
          fclose(file_in);
     }
                                                             
     MPI_Barrier(MPI_COMM_WORLD);
     MPI_Scatterv(dat_glo,counts,offsets,MPI_REAL,dat_loc,
                  nx*my,MPI_REAL,0,MPI_COMM_WORLD);

     for(iy = 2; iy <= my+1; iy++)
     { 
          for(ix = 0; ix <= nx-1; ix++)
          { 
               h_a[iy][ix] = (double)dat_loc[nx*(iy-2)+ix];
          }
     }
     //==============    reading u
     if(myid == 0)
     {
          printf("Reading initial data for u \n");                      
          file_in=fopen("shallow_init_u.bin","rb");                               
          fread(dat_glo, sizeof(float), nx * ny, file_in);                       
          fclose(file_in); 
     }
                                                             
     MPI_Barrier(MPI_COMM_WORLD);
     MPI_Scatterv(dat_glo,counts,offsets,MPI_REAL,dat_loc,
                  nx*my,MPI_REAL,0,MPI_COMM_WORLD);
     
     for(iy = 2; iy <= my+1; iy++)
     { 
          for(ix = 0; ix <= nx-1; ix++)
          { 
               u_a[iy][ix] = (double)dat_loc[nx*(iy-2)+ix];
          }
     }
     //==============    reading v
     if(myid == 0)
     {
          printf("Reading initial data for v \n");                      
          file_in = fopen("shallow_init_v.bin","rb");                               
          fread(dat_glo, sizeof(float), nx * ny, file_in);                       
          fclose(file_in);
     }

     MPI_Barrier(MPI_COMM_WORLD);
     MPI_Scatterv(dat_glo,counts,offsets,MPI_REAL,dat_loc,
                  nx*my,MPI_REAL,0,MPI_COMM_WORLD);
                                                             
     for(iy = 2; iy <= my+1; iy++)
     { 
          for(ix = 0; ix <= nx-1; ix++)
          { 
               v_a[iy][ix] = (double)dat_loc[nx*(iy-2)+ix];
          }
     }
}
/**************************************************************/
/**************************************************************/
void write_to_file(double ** u_b, double ** v_b, double ** h_b,
                   int myid, int nx, int my, int ny,              
                   int irec, const int * dims, int nproc)                                  

{
     FILE *file_out;
     int ix,iy,ip;
     int offsets[nproc], counts[nproc];
     float dat_loc[my * nx], dat_glo[ny * nx];

     for(ip = 0; ip <= nproc -1 ; ip++)
     {
          counts[ip] = dims[ip]*nx;
     }
     
     offsets[0] = 0;

     if(nproc > 1)
     {
          for(ip = 1; ip <= nproc -1 ; ip++)
          {
               offsets[ip] = offsets[ip-1] + counts[ip-1];
          }
     }

     if(myid == 0)
     {
          printf(" Archiving data irec=,%10d\n", irec);                      
     }

     //==============    archiving h
     for(iy = 2; iy <= my+1; iy++)
     { 
          for(ix = 0; ix <= nx-1; ix++)
          { 
               dat_loc[nx*(iy-2)+ix] = (float)h_b[iy][ix];
          }
     }

     MPI_Gatherv(dat_loc,ny*my,MPI_FLOAT,dat_glo,
                 counts,offsets,MPI_FLOAT,0,MPI_COMM_WORLD);

     if(myid == 0)
     {
          if (irec == 1)
          { 
               file_out=fopen("shallow_h.bin","wb");
          }                               
          else
          { 
               file_out=fopen("shallow_h.bin","ab");
          }                               
          fwrite(dat_glo, sizeof(float), nx * ny, file_out);                       
          fclose(file_out);
     }                                                        

     //==============    archiving u
     for(iy = 2; iy <= my+1; iy++)
     { 
          for(ix = 0; ix <= nx-1; ix++)
          { 
               dat_loc[nx*(iy-2)+ix] = (float) u_b[iy][ix];
          }
     }

     MPI_Gatherv(dat_loc,ny*my,MPI_FLOAT,dat_glo,
                 counts,offsets,MPI_FLOAT,0,MPI_COMM_WORLD);

     if(myid == 0)
     {
          if (irec == 1)
          { 
               file_out=fopen("shallow_u.bin","wb");
          }                               
          else
          { 
               file_out=fopen("shallow_u.bin","ab");
          }                               
          fwrite(dat_glo, sizeof(float), nx * ny, file_out);                       
          fclose(file_out);
     }
                                                             
     //==============    archiving v
     for(iy = 2; iy <= my+1; iy++)
     { 
          for(ix = 0; ix <= nx-1; ix++)
          { 
               dat_loc[nx*(iy-2)+ix] = (float)v_b[iy][ix];
          }
     }

     MPI_Gatherv(dat_loc,ny*my,MPI_FLOAT,dat_glo,
                 counts,offsets,MPI_FLOAT,0,MPI_COMM_WORLD);

     if(myid == 0)
     {
          if (irec == 1)
          { 
               file_out=fopen("shallow_v.bin","wb");
          }                               
          else
          { 
               file_out=fopen("shallow_v.bin","ab");
          }                               
          fwrite(dat_glo, sizeof(float), nx * ny, file_out);                       
          fclose(file_out);
     }                                                                 
}
/**************************************************************/
/**************************************************************/

