c Adapted from test example in SAMRAI distribution.
c
c Includes frunctions               for im_ex=0 or 1:
c                     comprhs3d and comprhsphase3d
c                                   for im_ex=2:
c                     comprhsimpphase implicit
c                     comprhsex3d and comprhsexphase explicit
c
define(NDIM,3)dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine comprhs3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  temperature, ngt,
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
      double precision latent_heat, cp
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
      double precision factor

      factor = latent_heat/cp
c
c  Computes RHS
c
      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            do ic0=ifirst0,ilast0

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
      subroutine addvel2rhs3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  dx,
     &  field, ng,
     &  frame_velocity,
     &  rhs)
c***********************************************************************
      implicit none
      double precision one
      parameter(one=1.d0)
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ng
      double precision dx(0:NDIM-1)
      double precision frame_velocity
      double precision
     &  field(CELL3d(ifirst,ilast,ng))
c output arrays:
      double precision
     &  rhs(CELL3d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer ic0,ic1,ic2
      double precision invdx

      invdx=1./dx(2)

      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            do ic0=ifirst0,ilast0
           rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2)+
     &         0.5d0*frame_velocity*invdx*(field(ic0,ic1,ic2+1)
     &         -field(ic0,ic1,ic2-1))
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
     &  mobility, well_height,
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
      double precision mobility, well_height
      double precision latent_heat, tmelting
      double precision phi(CELL3d(ifirst,ilast,ngp))
      double precision temp(CELL3d(ifirst,ilast,ngt))
c output arrays:
      double precision rhs(CELL3d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer ic0,ic1,ic2
      double precision phil, diffusion, tmp
c
      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            do ic0=ifirst0,ilast0

c        compute  RHS
            phil = phi(ic0,ic1,ic2)
            tmp =
     &         -32.d0*well_height*phil*(1.d0-phil)*(1.d0-2.d0*phil)
            tmp = tmp
     &         -6.d0*latent_heat*(temp(ic0,ic1,ic2)-tmelting)
     &                          *phil*(1.d0-phil)/tmelting

            rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2) + tmp*mobility
            enddo
         enddo
      enddo
c
      return
      end
c
c
c
      subroutine comprhsex3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  dx,
     &  field, ng,
     &  frame_velocity,
     &  rhs)
c***********************************************************************
      implicit none
      double precision one
      parameter(one=1.d0)
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ng
      double precision dx(0:NDIM-1)
      double precision frame_velocity
      double precision
     &  field(CELL3d(ifirst,ilast,ng))
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
      double precision invdx2(0:3-1)
      double precision factor
      double precision vel

      vel = frame_velocity

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

c        compute  RHS

           rhs(ic0,ic1,ic2) =
     &         0.5d0*vel/dx(2)*(field(ic0,ic1,ic2+1)
     &         -field(ic0,ic1,ic2-1))

            enddo
         enddo
      enddo
c
      return
      end
c
c
      subroutine comprhsdiffusion3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  dx,
     &  field, ng,
     &  diffusivity,
     &  rhs)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ng
      double precision dx(0:2)
      double precision diffusivity
      double precision field(CELL3d(ifirst,ilast,ng))
c output arrays:
      double precision rhs(CELL3d(ifirst,ilast,0))
c***********************************************************************
c
      integer ic0,ic1,ic2
      double precision dgrade0p, dgrade0m,
     &                 dgrade1p, dgrade1m,
     &                 dgrade2p, dgrade2m
      double precision invdx2(0:2)

      invdx2(0)=1./(dx(0)*dx(0))
      invdx2(1)=1./(dx(1)*dx(1))
      invdx2(2)=1./(dx(2)*dx(2))
c
c  Computes RHS for diffusion eqn
c
c      RHS = div ( D * grad(y) )
c
      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            do ic0=ifirst0,ilast0
               dgrade0p = (field(ic0+1,ic1,ic2)
     &                    -field(ic0,ic1,ic2))
               dgrade0m = (field(ic0,ic1,ic2)
     &                    -field(ic0-1,ic1,ic2))
               dgrade1p = (field(ic0,ic1+1,ic2)
     &                    -field(ic0,ic1,ic2))
               dgrade1m = (field(ic0,ic1,ic2)
     &                    -field(ic0,ic1-1,ic2))
               dgrade2p = (field(ic0,ic1,ic2+1)
     &                    -field(ic0,ic1,ic2))
               dgrade2m = (field(ic0,ic1,ic2)
     &                    -field(ic0,ic1,ic2-1))

c        compute  RHS

            rhs(ic0,ic1,ic2) = (dgrade0p - dgrade0m)*invdx2(0) +
     &                         (dgrade1p - dgrade1m)*invdx2(1) +
     &                         (dgrade2p - dgrade2m)*invdx2(2)
            rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2)*diffusivity
            enddo
         enddo
      enddo
c
      return
      end
