c//Example 1 with  a simple C code
            f1=['#include <math.h>'
            'void fooc(double c[],double a[],double *b,int *m,int *n)'
            '{'
            '   int i;'
            '   for ( i =0 ; i < (*m)*(*n) ; i++) '
            '     c[i] = sin(a[i]) + *b; '
            '}'];
            
            mputl(f1,'fooc.c')
            
            //creating the shared library (a gateway, a Makefile and a loader are
            //generated.
            
            ilib_for_link('fooc','fooc.c',[],"c")
            // load the shared library
            
            exec loader.sce
            
            //using the new primitive
            a=[1,2,3;4,5,6];b= %pi;
            [m,n]=size(a);
            
            // Inputs:
            // a is in position 2 and double
            // b                3     double
            // n                4     integer
            // m                5     integer
            // Outputs:
            // c is in position 1 and double with size [m,n]
            c=call("fooc",a,2,"d",b,3,"d",m,4,"i",n,5,"i","out",[m,n],1,"d");
            
