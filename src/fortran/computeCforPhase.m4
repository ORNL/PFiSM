define(NDIM,3)dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine compcforphase(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  phi, ngphi,
     &  mobility, dwell_height, gamma,
     &  cfield, ngc)
c***********************************************************************
      implicit none
      double precision one
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ngphi, ngc
      double precision mobility
      double precision gamma
      double precision dwell_height
      double precision
     &  phi(CELL3d(ifirst,ilast,ngphi))
c output arrays:
      double precision
     &  cfield(CELL3d(ifirst,ilast,ngc))
c
c***********************************************************************
c
      integer ic0,ic1,ic2
      double precision factor
      double precision dwell2ndderiv
      double precision gamma_inv
c
      factor = mobility*32.d0*dwell_height
      gamma_inv = 1.d0/gamma
      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            do ic0=ifirst0,ilast0
               dwell2ndderiv = ( 1.d0 + 6.d0 * phi(ic0,ic1,ic2) 
     &                                  * ( phi(ic0,ic1,ic2) - 1.d0 ) )
               cfield(ic0,ic1,ic2) = gamma_inv 
     &                             + factor*dwell2ndderiv
            enddo
         enddo
      enddo
c
      return
      end
c

