      subroutine foof(c,a,b,n,m)
      integer n,m
      double precision a(*),b,c(*)
      do 10 i=1,m*n 
        c(i) = sin(a(i))+b
   10 continue
      end
