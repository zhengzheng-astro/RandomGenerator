#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "random.h"
#include <time.h>

#define N 1000000

int main(int argc, char **argv)
{
 mt_state_t vstate;
 int     i,seed,option;
 double  T,a,x,uz;
 double  TotalTime;
 clock_t begin,end;

 if(argc!=4) {
   printf("run time test for drawing uz\n");
   printf("Usage: %s -seed option T\n",argv[0]);
   printf("   seed   - random seed\n");
   printf("   option - 0: ZM02\n");
   printf("            1: LZ26 (optimized)\n");
   printf("            2: Semelin2007\n");
   printf("            3: RASCAS2020\n");
   printf("            4: Smith2015\n");
   printf("   T      - Temperature (K)\n");
   //printf("   a      - Voigt parameter, Delta nu_L/(2 Delta nu_D)=4.7e-4(T/10^4K)^{-1/2}\n");
   exit(-1);
 }

 seed  =atoi(argv[1]);
 if(seed<0) seed=-seed;
 option=atoi(argv[2]);
 T     =atof(argv[3]);
 a     =4.70e-4/sqrt(T/1e4);

 i=N;
 printf("# seed = %d ,  method = %d ,  T = %11.4e ( a = %e ) , N = %8d\n",seed,option,T,a,i);
 printf("# x     time_per_draw (sec)\n");

 seed_ran_mt(&vstate, seed);
 for(x=0;x<21;x=x+0.1) {
    switch(option) {
       case 0:
            begin=clock();
            for(i=0;i<N;i++) get_uz0(&uz,a,x,&vstate);
            end=clock();
            break;
       case 1:
            begin=clock();
            for(i=0;i<N;i++) get_uz_opt(&uz,a,x,&vstate);
            end=clock();
            break;
       case 2:
            begin=clock();
            for(i=0;i<N;i++) get_uz_Semelin(&uz,a,x,&vstate);
            end=clock();
            break;
       case 3:
            begin=clock();
            for(i=0;i<N;i++) get_uz_RASCAS(&uz,a,x,&vstate);
            end=clock();
            break;
       default:
            begin=clock();
            for(i=0;i<N;i++) get_uz_Smith(&uz,a,x,&vstate);
            end=clock();
    }
  TotalTime=(double)(end-begin)/CLOCKS_PER_SEC;
  printf("%7.4f %e\n",x,TotalTime/N);
  fflush(stdin);
 }

}
