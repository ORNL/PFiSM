c Adapted from test example in SAMRAI distribution.
c
define(NDIM,3)dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine comprhs3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  ngc0,ngc1,ngc2,
     &  dx,
     &  y,
     &  diff0,
     &  diff1,
     &  diff2,
     &  rhs)
c***********************************************************************
      implicit none
      double precision one
      parameter(one=1.d0)
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ngc0,ngc1,ngc2
      double precision dt, dx(0:NDIM-1)
      double precision
     &  y(CELL3dVECG(ifirst,ilast,ngc)),
     &  diff0(SIDE3d0(ifirst,ilast,0)),
     &  diff1(SIDE3d1(ifirst,ilast,0)),
     &  diff2(SIDE3d2(ifirst,ilast,0))
c output arrays:
      double precision
     &  rhs(CELL3d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer ic0,ic1,ic2
      double precision dgrade0p, dgrade0m, 
     &                 dgrade1p, dgrade1m,
     &                 dgrade2p, dgrade2m
      double precision invdx2(0:NDIM-1)

      invdx2(0)=1./(dx(0)*dx(0))
      invdx2(1)=1./(dx(1)*dx(1))
      invdx2(2)=1./(dx(2)*dx(2))

c
c  Computes RHS for 1 eqn diffusion
c
c      RHS = div ( D * grad(y) )
c
      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            do ic0=ifirst0,ilast0

c        compute  D(E)grad(E) in X, Y, and Z

            dgrade0p = diff0(ic0+1,ic1,ic2) *
     &                 (y(ic0+1,ic1,ic2) - y(ic0,ic1,ic2))
            dgrade0m = diff0(ic0,ic1,ic2) *
     &                 (y(ic0,ic1,ic2) - y(ic0-1,ic1,ic2))
            dgrade1p = diff1(ic0,ic1+1,ic2) *
     &                 (y(ic0,ic1+1,ic2) - y(ic0,ic1,ic2))
            dgrade1m = diff1(ic0,ic1,ic2) *
     &                 (y(ic0,ic1,ic2) - y(ic0,ic1-1,ic2))
            dgrade2p = diff2(ic0,ic1,ic2+1) *
     &                 (y(ic0,ic1,ic2+1) - y(ic0,ic1,ic2))
            dgrade2m = diff1(ic0,ic1,ic2) *
     &                 (y(ic0,ic1,ic2) - y(ic0,ic1,ic2-1))

c        compute  RHS

            rhs(ic0,ic1,ic2) = (dgrade0p - dgrade0m)*invdx2(0) +
     &                         (dgrade1p - dgrade1m)*invdx2(1) +
     &                         (dgrade2p - dgrade2m)*invdx2(2)

            enddo
         enddo
      enddo
c
      return
      end
c
