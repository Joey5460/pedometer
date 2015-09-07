#ifndef _MATH_ECG_H_
#define _MATH_ECG_H_
void bf_max_min(const int * buf, unsigned int size, int *max, int *min);
double bf_mean(const int *buf, unsigned int size);
double bf_log2(double x);
void filtering(int * dataY,unsigned int size,int sr);
int down_sample(int vin, int *vout, int ifreq, int ofreq);

unsigned short step_counter(unsigned short data, unsigned short sr);
void reset_pmd(void);
#endif

