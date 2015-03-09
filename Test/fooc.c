#include <math.h>
void fooc(double c[],double a[],double *b,int *m,int *n)
{
   int i;
   for ( i =0 ; i < (*m)*(*n) ; i++) 
     c[i] = sin(a[i]) + *b; 
}
