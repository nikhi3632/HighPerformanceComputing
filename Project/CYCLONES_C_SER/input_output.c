
#include <stdlib.h>                                                             
#include <stdio.h>                                                              
/**************************************************************/
/**************************************************************/
void read_from_file(double ** u_a, double ** v_a, double ** h_a,               
                    int nx, int ny)                                  
{
     FILE  *file_in;
     int ix,iy;
     float dat_glo[ny * nx];

     //==============    reading h
     printf("Reading initial data for h \n");                      

     file_in=fopen("shallow_init_h.bin","rb");                               
     fread(dat_glo, sizeof(float), nx * ny, file_in);                       
     fclose(file_in);                                                        

     for(iy = 0; iy <= ny-1; iy++)
     { 
          for(ix = 0; ix <= nx-1; ix++)
          { 
               h_a[iy][ix] = (double)dat_glo[nx*iy+ix];
          }
     }
     //==============    reading u
     printf("Reading initial data for u \n");                      

     file_in=fopen("shallow_init_u.bin","rb");                               
     fread(dat_glo, sizeof(float), nx * ny, file_in);                       
     fclose(file_in);                                                        

     for(iy = 0; iy <= ny-1; iy++)
     { 
          for(ix = 0; ix <= nx-1; ix++)
          { 
               u_a[iy][ix] = (double)dat_glo[nx*iy+ix];
          }
     }
     //==============    reading v
     printf("Reading initial data for v \n");                      

     file_in = fopen("shallow_init_v.bin","rb");                               
     fread(dat_glo, sizeof(float), nx * ny, file_in);                       
     fclose(file_in);                                                        

     for(iy = 0; iy <= ny-1; iy++)
     { 
          for(ix = 0; ix <= nx-1; ix++)
          { 
               v_a[iy][ix] = (double)dat_glo[nx*iy+ix];
          }
     }
}
/**************************************************************/
/**************************************************************/
void write_to_file(double ** u_b, double ** v_b, double ** h_b,               
                   int irec,int nx, int ny)                                  

{
     FILE  *file_out;
     int ix,iy;
     float dat_glo[ny * nx];

     printf(" Archiving data irec=,%10d\n", irec);                      
     //==============    archiving h
     for(iy = 0; iy <= ny-1; iy++)
     { 
          for(ix = 0; ix <= nx-1; ix++)
          { 
               dat_glo[nx*iy+ix] = (float)h_b[iy][ix];
          }
     }

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

     //==============    archiving u
     for(iy = 0; iy <= ny-1; iy++)
     { 
          for(ix = 0; ix <= nx-1; ix++)
          { 
               dat_glo[nx*iy+ix] = (float) u_b[iy][ix];
          }
     }
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
     //==============    archiving v
     for(iy = 0; iy <= ny-1; iy++)
     { 
          for(ix = 0; ix <= nx-1; ix++)
          { 
               dat_glo[nx*iy+ix] = (float)v_b[iy][ix];
          }
     }
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
/**************************************************************/
/**************************************************************/

