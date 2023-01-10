#include <stdio.h>                                                              

/**************************************************************/
void write_to_file(double ** h_b, int irec,int nx, int ny)                    
{
    FILE  *file_out;                                                            
    int ix , iy ;
    float dat_glo[ny * nx];

    printf(" Archiving data irec = ,%10d\n", irec);                      
    //==============    archiving h
    for(iy = 0; iy <= ny-1; iy++)
    {  
        for(ix = 0; ix <= nx-1; ix++)
        { 
            dat_glo[nx*iy+ix] = (float) h_b[iy][ix];
        }
    }

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
/**************************************************************/
/**************************************************************/
