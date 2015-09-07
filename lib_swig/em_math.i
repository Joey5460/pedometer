%module em_math

%{
#define SWIG_FILE_WITH_INT
#include "em_math.h"
%}
unsigned short step_counter(unsigned short data, unsigned short sr);
void reset_pmd(void);

%include "carrays.i"
%array_class(int, intArray);
void filtering(int * dataY,unsigned int size,int sr);
%include "typemaps.i"
%apply int *OUTPUT {int *vout};
int down_sample(int vin, int *vout, int ifreq, int ofreq);

