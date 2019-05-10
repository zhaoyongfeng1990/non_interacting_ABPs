#include <sstream>
#include <cstring>
#include <ctime>
#include "particles.h"

int main(int argc, char *argv[])
{
    MPI_Init(NULL, NULL);
    int pid=0;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    int NP=0;
    MPI_Comm_size(MPI_COMM_WORLD, &NP);
    uint64_t localN=NumSample/NP;

    uint64_t seed = 0;
    clock_t it, ft;
    it=clock();
    ran.init(pid+it,pid+it+1);
    
    MPI_Type_contiguous(4*3, MPI_DOUBLE, &MPI_PARTICLE);
    MPI_Type_commit(&MPI_PARTICLE);

    particle pList[NumSample];

    double totaltime = 0;
    double time = 0;
    for (uint64_t i = 0; i < localN; ++i)
    {
        pList[i].x = ran.random4d()*Lx;
        pList[i].y = ran.random4d()*Ly;
        pList[i].theta = ran.random4d()*PI;
    }
    
    collectingResult(pList);

    FILE *outfile = NULL;
    int timepoint = 0;
    if (pid == 0)
    {
        char filename[26]; //file name will be like p0000.txt
        char buffer[13];
        sprintf(buffer, "%06d", timepoint);
        strcpy(filename, "p");
        strcat(filename, buffer);
        strcat(filename, ".txt");
        outfile = fopen(filename, "w");
        for (uint64_t i = 0; i < NumSample; i++)
        {  
            for(int j=0; j<4; ++j)
                fprintf(outfile, "%e %e %e ", pList[i].x[j], pList[i].y[j], pList[i].theta[j]);
        }
        //for (uint64_t i = 0; i < NLocalBoxes; i++)
        //{
        //    for(uint64_t j=0; j<bList[i].NumParticles; j++)
        //        fprintf(outfile, "%e %e %e %e ", bList[i].partList[j].x, bList[i].partList[j].y, bList[i].partList[j].theta);
        //}
        fprintf(outfile, "\n");
        fflush(outfile);
        fclose(outfile);
        ++timepoint;
    }

    uint64_t idx=1;
    uint64_t count=0;
    while (time < TotalTime)
    {
        for (uint64_t i = 0; i < localN; ++i)
        {
            step(&pList[i]);
        }
        time += dt;
        
        if(time >= idx)
        {
            collectingResult(pList);
            if (pid == 0)
            {
                char filename[26]; //file name will be like p0000.txt
                char buffer[13];
                sprintf(buffer, "%06d", timepoint);
                strcpy(filename, "p");
                strcat(filename, buffer);
                strcat(filename, ".txt");
                outfile = fopen(filename, "w");
                for (uint64_t i = 0; i < NumSample; i++)
                {  
                    for(int j=0; j<4; ++j)
                        fprintf(outfile, "%e %e %e ", pList[i].x[j], pList[i].y[j], pList[i].theta[j]);
                }
                //for (uint64_t i = 0; i < NLocalBoxes; i++)
                //{
                //    for(uint64_t j=0; j<bList[i].NumParticles; j++)
                //        fprintf(outfile, "%e %e %e %e ", bList[i].partList[j].x, bList[i].partList[j].y, bList[i].partList[j].theta);
                //}
                fprintf(outfile, "\n");
                fflush(outfile);
                fclose(outfile);
                ++timepoint;
            }
            ++idx;
        }
    }

    MPI_Type_free(&MPI_PARTICLE);
    ft = clock();
    printf("%f\n", (double)(ft - it) / CLOCKS_PER_SEC);
    MPI_Finalize();
    return 0;
}