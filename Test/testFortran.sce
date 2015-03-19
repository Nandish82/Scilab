//Example 2 with  a simple Fortran code
            f1=['      subroutine foof(c,a,b,n,m)'
            '      integer n,m'
            '      double precision a(*),b,c(*)'
            '      do 10 i=1,m*n '
            '        c(i) = sin(a(i))+b'
            '   10 continue'
            '      end'];
            mputl(f1,'foof.f')
            
            //creating the shared library (a gateway, a Makefile and a loader are
            //generated.
            
            ilib_for_link('foof','foof.f',[],"f")
            
            // load the shared library
            
            exec loader.sce
            
            //using the new primitive
            a=[1,2,3;4,5,6];b= %pi;
            [m,n]=size(a);
            c=call("foof",a,2,"d",b,3,"d",m,4,"i",n,5,"i","out",[m,n],1,"d");
