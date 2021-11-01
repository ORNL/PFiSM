c Adapted from test example in SAMRAI distribution.
c
define(NDIM,3)dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine comprhs3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  dx,
     &  temperature, ngt,
     &  diffusivity,
     &  dphidt,
     &  latent_heat, cp,
     &  rhs)
c***********************************************************************
      implicit none
      double precision one
      parameter(one=1.d0)
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ngt
      double precision dx(0:NDIM-1)
      double precision diffusivity, latent_heat, cp
      double precision
     &  temperature(CELL3d(ifirst,ilast,ngt)),
     &  dphidt(CELL3d(ifirst,ilast,0))
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
      double precision factor

      invdx2(0)=1./(dx(0)*dx(0))
      invdx2(1)=1./(dx(1)*dx(1))
      invdx2(2)=1./(dx(2)*dx(2))

      factor = latent_heat/cp
c
c  Computes RHS for 1 eqn diffusion
c
c      RHS = div ( D * grad(y) )
c
      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            do ic0=ifirst0,ilast0

               dgrade0p = (temperature(ic0+1,ic1,ic2)
     &                    -temperature(ic0,ic1,ic2))
               dgrade0m = (temperature(ic0,ic1,ic2)
     &                    -temperature(ic0-1,ic1,ic2))
               dgrade1p = (temperature(ic0,ic1+1,ic2)
     &                    -temperature(ic0,ic1,ic2))
               dgrade1m = (temperature(ic0,ic1,ic2)
     &                    -temperature(ic0,ic1-1,ic2))
               dgrade2p = (temperature(ic0,ic1,ic2+1)
     &                    -temperature(ic0,ic1,ic2))
               dgrade2m = (temperature(ic0,ic1,ic2)
     &                    -temperature(ic0,ic1,ic2-1))

c        compute  RHS

            rhs(ic0,ic1,ic2) = (dgrade0p - dgrade0m)*invdx2(0) +
     &                         (dgrade1p - dgrade1m)*invdx2(1) +
     &                         (dgrade2p - dgrade2m)*invdx2(2)
            rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2)*diffusivity
            rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2)
     &                       + factor*dphidt(ic0,ic1,ic2)
            enddo
         enddo
      enddo
c
      return
      end
c
c

      subroutine comprhsphase3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  phi, ngp,
     &  temp, ngt,
     &  dx,
     &  mobility, well_height, epsilon,
     &  latent_heat, tmelting,
     &  rhs)
c***********************************************************************
      implicit none
      double precision one
      parameter(one=1.d0)
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ngp, ngt
      double precision dx(0:NDIM-1)
      double precision mobility, well_height, epsilon
      double precision latent_heat, tmelting
      double precision
     &  phi(CELL3d(ifirst,ilast,ngp))
      double precision
     &  temp(CELL3d(ifirst,ilast,ngt))
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
      double precision phil, diffusion
      double precision invdx2(0:NDIM-1)

      invdx2(0)=1./(dx(0)*dx(0))
      invdx2(1)=1./(dx(1)*dx(1))
      invdx2(2)=1./(dx(2)*dx(2))

      diffusion = epsilon*epsilon
c
c  Computes RHS for 1 eqn diffusion
c
c      RHS = div ( D * grad(y) )
c
      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            do ic0=ifirst0,ilast0

c        compute  D(E)grad(E) in X, Y, and Z

            dgrade0p = phi(ic0+1,ic1,ic2) - phi(ic0,ic1,ic2)
            dgrade0m = phi(ic0,ic1,ic2) - phi(ic0-1,ic1,ic2)
            dgrade1p = phi(ic0,ic1+1,ic2) - phi(ic0,ic1,ic2)
            dgrade1m = phi(ic0,ic1,ic2) - phi(ic0,ic1-1,ic2)
            dgrade2p = phi(ic0,ic1,ic2+1) - phi(ic0,ic1,ic2)
            dgrade2m = phi(ic0,ic1,ic2) - phi(ic0,ic1,ic2-1)

c        compute  RHS

            rhs(ic0,ic1,ic2) = (dgrade0p - dgrade0m)*invdx2(0) +
     &                         (dgrade1p - dgrade1m)*invdx2(1) +
     &                         (dgrade2p - dgrade2m)*invdx2(2)
            rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2)*diffusion

            phil = phi(ic0,ic1,ic2)
            rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2)
     &         -32.d0*well_height*phil*(1.d0-phil)*(1.d0-2.d0*phil)
            rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2)
     &         -6.d0*latent_heat*(temp(ic0,ic1,ic2)-tmelting)
     &                          *phil*(1.d0-phil)/tmelting

            rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2)*mobility
            enddo
         enddo
      enddo
c
      return
      end
c

