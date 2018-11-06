c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2018 Lawrence Livermore National Security, LLC
c
      subroutine setexactandrhs3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  exact,rhs,dx,xlower)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
c variables in 1d axis indexed
c
      double precision 
     &     dx(0:3-1),
     &     xlower(0:3-1)
c variables in 3d cell indexed         
      double precision
     &     exact(ifirst0-1:ilast0+1,
     &          ifirst1-1:ilast1+1,
     &          ifirst2-1:ilast2+1),
     &     rhs(ifirst0:ilast0,
     &          ifirst1:ilast1,
     &          ifirst2:ilast2)
c
c***********************************************************************     
c
      integer ic0,ic1,ic2
      double precision x, y, z, sinsin, pi

      pi=3.141592654

c     write(6,*) "In fluxcorrec()"
c     ******************************************************************
      do ic2=ifirst2,ilast2
         z = xlower(2) + dx(2)*(ic2-ifirst2+0.5)
         do ic1=ifirst1,ilast1
            y = xlower(1) + dx(1)*(ic1-ifirst1+0.5)
            do ic0=ifirst0,ilast0
               x = xlower(0) + dx(0)*(ic0-ifirst0+0.5)
               sinsin = sin(pi*x) * sin(pi*y) * sin(pi*z)
               exact(ic0,ic1,ic2) = sinsin
               rhs(ic0,ic1,ic2) = -3*pi*pi*sinsin
            enddo
         enddo
      enddo

      return
      end   
c
