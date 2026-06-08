#ifndef RANDOM_H
#define RANDOM_H

typedef struct
  {
    unsigned long mt[624];
    int mti;
  }
mt_state_t;

double ran_mt(void *vstate);
void seed_ran_mt(void *vstate, unsigned long int s);

void get_uz0(double *uz, double a, double x, void *vstate);
double thres_expu2(void *vstate, double thres);
double bin_expu2(void *vstate, double thres1, double thres2);
double bin2_expu2(void *vstate, double thres1, double thres2);
void get_uz1(double *uz, double a, double xt, void *vstate);
void get_uz_opt(double *uz, double a, double xt, void *vstate);

void get_uz_Semelin(double *uz, double a, double x, void *vstate);
void get_uz_RASCAS(double *uz, double a, double x, void *vstate);
void get_uz_Smith(double *uz, double a, double x, void *vstate);

#endif
