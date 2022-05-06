#include "defines.h"

  module kernels_fluid_cg
    use dimensions_m
    use kernels_fluid
    implicit none
    
   private :: equil
   private :: equilCG
   private :: tensor_product
   private :: linear,linear2
   
  contains

    attributes(device) pure function linear(i,j,k)
      integer, intent(in) :: i,j,k
      integer :: linear

      linear = i*1000**2 + j*1000 + k
    end function linear

    attributes(device) pure function linear2(i,j,k)
      integer, intent(in) :: i,j,k      
      integer :: linear2

      linear2 = (50+i)*1000**2 + (50+j)*1000 + (50+k)
    end function linear2
    
    attributes(device) function equil(rho,u,v,w, l)
     real, intent(in) :: rho,u,v,w
     integer, intent(in) :: l
     real :: equil
     real :: uv

     uv = (1.0/cssq) * (u*ex_d(l) + v*ey_d(l) + w*ez_d(l))
     equil = rho * p_d(l)*(ONE+uv+HALF*(uv*uv)-(HALF/cssq) * (u*u + v*v + w*w))
    end function equil
    
    attributes(device) function equilCG(rho,u,v,w, k, l)
     real, intent(in) :: rho,u,v,w
     integer, intent(in) :: k,l
     real :: equilCG
     real :: uv

     uv = (1.0/cssq) * (u*ex_d(l) + v*ey_d(l) + w*ez_d(l))
     equilCG = rho * ( phi_d(l) + varphi_d(l)*alphaCG_d(k) + &
      p_d(l)*(uv+HALF*(uv*uv)-(HALF/cssq) * (u*u + v*v + w*w)))
    end function equilCG
    
    attributes(device) pure function tensor_product(v,w)

     implicit none

     real, intent(in), dimension(3) :: v,w
     real, dimension(3,3) :: tensor_product
  
     tensor_product(1,1)=v(1)*w(1)
     tensor_product(2,1)=v(2)*w(1)
     tensor_product(3,1)=v(3)*w(1)
     
     tensor_product(1,2)=v(1)*w(2)
     tensor_product(2,2)=v(2)*w(2)
     tensor_product(3,2)=v(3)*w(2)
     
     tensor_product(1,3)=v(1)*w(3)
     tensor_product(2,3)=v(2)*w(3)
     tensor_product(3,3)=v(3)*w(3)
     
     return

    end function tensor_product
    



    attributes(device) function fcut(r, r1, r2)
     real, intent(in) :: r, r1, r2
     real fcut

     if ( r <= r1 ) then
        fcut = ONE
     elseif ( r > r2 ) then
        fcut = ZERO
     else
        fcut = HALF * cos((r-r1)*Pi/(r2-r1)) + HALF
     endif
    end function fcut
    
    attributes(global) subroutine setup1(vx,vy,vz)
      real, value :: vx,vy,vz
      integer :: i,j,k, l
      integer :: gli,glj,glk
      real    :: rhoR,rhoB, u,v,w, eqR,eqB, distx,disty,distz, rdist,tempr
      real,parameter    :: radius = glz/8
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z

      gli = i + offset_d(1)
      glj = j + offset_d(2)
      glk = k + offset_d(3)

      if (1==i*j*k+offset_d(3)) write(*,*) 'CUDA setup2] offset:', gli, glj, glk

      distx = gli - (glx/2 + 0.5)
      disty = glj - (gly/2 + 0.5)
      distz = glk - (glz/2 + 0.5)

      !rdist = sqrt(distx*distx + disty*disty + distz*distz)
      tempr = 1.0 !fcut(rdist, radius, radius+0.1)
#ifdef DIFFDENS
      rhoR = tempr*densR_d
      rhoB = (1.0 - tempr)*densB_d
#else
      rhoR = tempr
      rhoB = (1.0 - tempr)
#endif
      u = vx
      v = vy
      w = vz
  
      rhoR_d(i,j,k) = rhoR
      rhoB_d(i,j,k) = rhoB
      do l = 0, npops-1
#ifdef DIFFDENS
	      eqR = equilCG(rhoR, u,v,w, 1, l)
	      eqB = equilCG(rhoB, u,v,w, 2, l)
#else
          eqR = equil(rhoR, u,v,w, l)
	      eqB = equil(rhoB, u,v,w, l)
#endif
        popsR_d(i,j,k, l, 1) = eqR
        popsB_d(i,j,k, l, 1) = eqB

        ! Clears other stuff
        popsR_d(i,j,k, l, 2) = 0.0
        popsB_d(i,j,k, l, 2) = 0.0
      end do

      myfluid_d(i,j,k, 1) = fluid_fluid
      myfluid_d(i,j,k, 2) = fluid_fluid
      debugfluid_d(i,j,k) = 0

      countmk_d = 0
      countrm_d = 0
      countn2p_d = 0
      countp2n1_d = 0
      countp2n2_d = 0
    end subroutine setup1

    attributes(global) subroutine setup2(vx,vy,vz)
      real, value :: vx,vy,vz
      integer :: i,j,k, l
      integer :: gli,glj,glk
      real    :: rhoR,rhoB, u,v,w, eqR,eqB, distx,disty,distz, rdist,tempr
      real,parameter    :: radius = glz/8
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z

      gli = i + offset_d(1)
      glj = j + offset_d(2)
      glk = k + offset_d(3)

      if (1==i*j*k+offset_d(3)) write(*,*) 'CUDA setup2] offset:', gli, glj, glk

      distx = gli - (glx/2 + 0.5)
      disty = glj - (gly/2 + 0.5)
      distz = glk - (glz/2 + 0.5)

      rdist = sqrt(distx*distx + disty*disty + distz*distz)
      tempr = fcut(rdist, radius, radius+0.1)
#ifdef DIFFDENS
      rhoR = tempr*densR_d
      rhoB = (1.0 - tempr)*densB_d
#else
      rhoR = tempr
      rhoB = (1.0 - tempr)
#endif
      u = vx
      v = vy
      w = vz
  
      rhoR_d(i,j,k) = rhoR
      rhoB_d(i,j,k) = rhoB
      do l = 0, npops-1
#ifdef DIFFDENS
	      eqR = equilCG(rhoR, u,v,w, 1, l)
	      eqB = equilCG(rhoB, u,v,w, 2, l)
#else
          eqR = equil(rhoR, u,v,w, l)
	      eqB = equil(rhoB, u,v,w, l)
#endif
        popsR_d(i,j,k, l, 1) = eqR
        popsB_d(i,j,k, l, 1) = eqB

        ! Clears other stuff
        popsR_d(i,j,k, l, 2) = 0.0
        popsB_d(i,j,k, l, 2) = 0.0
      end do

      myfluid_d(i,j,k, 1) = fluid_fluid
      myfluid_d(i,j,k, 2) = fluid_fluid
      debugfluid_d(i,j,k) = 0

      countmk_d = 0
      countrm_d = 0
      countn2p_d = 0
      countp2n1_d = 0
      countp2n2_d = 0
    end subroutine setup2
    
    
    attributes(global) subroutine setup2_sphere(vx,vy,vz,cx,cy,cz,radius)
      real, value :: vx,vy,vz,cx,cy,cz,radius
      integer :: i,j,k, l
      integer :: gli,glj,glk
      real    :: rhoR,rhoB, u,v,w, eqR,eqB, distx,disty,distz, rdist,tempr
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z

      gli = i + offset_d(1)
      glj = j + offset_d(2)
      glk = k + offset_d(3)

      if (1==i*j*k+offset_d(3)) write(*,*) 'CUDA setup2] offset:', gli, glj, glk

      distx = gli - cx
      disty = glj - cy
      distz = glk - cz

      rdist = sqrt(distx*distx + disty*disty + distz*distz)
      tempr = fcut(rdist, radius, radius+0.1)
#ifdef DIFFDENS
      rhoR = tempr*densR_d
      rhoB = (1.0 - tempr)*densB_d
#else
      rhoR = tempr
      rhoB = (1.0 - tempr)
#endif
      u = vx
      v = vy
      w = vz
  
      rhoR_d(i,j,k) = rhoR
      rhoB_d(i,j,k) = rhoB
      do l = 0, npops-1
#ifdef DIFFDENS
	      eqR = equilCG(rhoR, u,v,w, 1, l)
	      eqB = equilCG(rhoB, u,v,w, 2, l)
#else
          eqR = equil(rhoR, u,v,w, l)
	      eqB = equil(rhoB, u,v,w, l)
#endif
        popsR_d(i,j,k, l, 1) = eqR
        popsB_d(i,j,k, l, 1) = eqB

        ! Clears other stuff
        popsR_d(i,j,k, l, 2) = 0.0
        popsB_d(i,j,k, l, 2) = 0.0
      end do

      myfluid_d(i,j,k, 1) = fluid_fluid
      myfluid_d(i,j,k, 2) = fluid_fluid
      debugfluid_d(i,j,k) = 0

      countmk_d = 0
      countrm_d = 0
      countn2p_d = 0
      countp2n1_d = 0
      countp2n2_d = 0
    end subroutine setup2_sphere
    
    attributes(global) subroutine setup2_cylinder(vx,vy,vz,cx,cy,cz, &
     radius,isel)
      real, value :: vx,vy,vz,cx,cy,cz,radius
      integer, value :: isel
      integer :: i,j,k, l
      integer :: gli,glj,glk
      real    :: rhoR,rhoB, u,v,w, eqR,eqB, distx,disty,distz, rdist,tempr
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
      gli = i + offset_d(1)
      glj = j + offset_d(2)
      glk = k + offset_d(3)


      distx = gli - cx
      disty = glj - cy
      distz = glk - cz
      select case(isel)
      case(1)
        rdist = sqrt(disty*disty + distz*distz)
      case(2)
        rdist = sqrt(distx*distx + distz*distz)
      case(3)
        rdist = sqrt(distx*distx + disty*disty)
      end select
      
      tempr = fcut(rdist, radius, radius+0.1)
#ifdef DIFFDENS
      rhoR = tempr*densR_d
      rhoB = (1.0 - tempr)*densB_d
#else
      rhoR = tempr
      rhoB = (1.0 - tempr)
#endif
      u = vx
      v = vy
      w = vz
  
      rhoR_d(i,j,k) = rhoR
      rhoB_d(i,j,k) = rhoB
      do l = 0, npops-1
#ifdef DIFFDENS
	      eqR = equilCG(rhoR, u,v,w, 1, l)
	      eqB = equilCG(rhoB, u,v,w, 2, l)
#else
          eqR = equil(rhoR, u,v,w, l)
	      eqB = equil(rhoB, u,v,w, l)
#endif
        popsR_d(i,j,k, l, 1) = eqR
        popsB_d(i,j,k, l, 1) = eqB

        ! Clears other stuff
        popsR_d(i,j,k, l, 2) = 0.0
        popsB_d(i,j,k, l, 2) = 0.0
      end do

      myfluid_d(i,j,k, 1) = fluid_fluid
      myfluid_d(i,j,k, 2) = fluid_fluid
      debugfluid_d(i,j,k) = 0

      countmk_d = 0
      countrm_d = 0
      countn2p_d = 0
      countp2n1_d = 0
      countp2n2_d = 0
    end subroutine setup2_cylinder
    
    attributes(global) subroutine setup2_collid(vx,vy,vz)
      real, value :: vx,vy,vz
      integer :: i,j,k, l
      integer :: gli,glj,glk
      real    :: rhoR,rhoB, u,v,w, eqR,eqB, distx,disty,distz, rdist,tempr
      real ::  distx2,disty2,distz2,gfi,gfj,gfk,rdist2
      real,parameter    :: radius = glz/6
      real,parameter    :: myvel=0.05
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z

      gli = i + offset_d(1)
      glj = j + offset_d(2)
      glk = k + offset_d(3)
      
      gfi=real(gli)
      gfj=real(glj)
      gfk=real(glk)

      if (1==i*j*k+offset_d(3)) write(*,*) 'CUDA setup2] offset:', gli, glj, glk

      distx = gfi - (real(glx)/4.0 + 0.5)
      disty = gfj - (real(gly)/2.0 + 0.5)
      distz = gfk - (real(glz)/2.0 + 0.5)
      
      distx2 = gfi - (3.0*real(glx)/4.0 + 0.5)
      disty2 = gfj - (real(gly)/2.0 + 0.5)
      distz2 = gfk - (real(glz)/2.0 + 0.5)

      rdist = sqrt(distx*distx + disty*disty + distz*distz)
      rdist2 = sqrt(distx2*distx2 + disty2*disty2 + distz2*distz2)
      tempr = fcut(rdist, radius, radius+0.1)+fcut(rdist2, radius, radius+0.1)
#ifdef DIFFDENS
      rhoR = tempr*densR_d
      rhoB = (1.0 - tempr)*densB_d
#else
      rhoR = tempr
      rhoB = (1.0 - tempr)
#endif
      u = vx+myvel*fcut(rdist, radius, radius+0.1)-myvel*fcut(rdist2, radius, radius+0.1)
      v = vy
      w = vz
  
      rhoR_d(i,j,k) = rhoR
      rhoB_d(i,j,k) = rhoB
      do l = 0, npops-1
#ifdef DIFFDENS
	      eqR = equilCG(rhoR, u,v,w, 1, l)
	      eqB = equilCG(rhoB, u,v,w, 2, l)
#else
          eqR = equil(rhoR, u,v,w, l)
	      eqB = equil(rhoB, u,v,w, l)
#endif
        popsR_d(i,j,k, l, 1) = eqR
        popsB_d(i,j,k, l, 1) = eqB

        ! Clears other stuff
        popsR_d(i,j,k, l, 2) = 0.0
        popsB_d(i,j,k, l, 2) = 0.0
      end do

      myfluid_d(i,j,k, 1) = fluid_fluid
      myfluid_d(i,j,k, 2) = fluid_fluid
      debugfluid_d(i,j,k) = 0

      countmk_d = 0
      countrm_d = 0
      countn2p_d = 0
      countp2n1_d = 0
      countp2n2_d = 0
    end subroutine setup2_collid

    attributes(global) subroutine setup2_zplanes(vx,vy,vz)
      real, value :: vx,vy,vz
      integer :: i,j,k, l
      integer :: gli,glj,glk
      real    :: rhoR,rhoB, u,v,w, eqR,eqB, tempr
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z

      gli = i + offset_d(1)
      glj = j + offset_d(2)
      glk = k + offset_d(3)

      if (1==i*j*k+offset_d(3)) write(*,*) 'CUDA setup2_zplanes] offset:', gli, glj, glk

      ! B @ 1/8-> 5/8 z-planes
      if (glz*1/8<glk .and. glk<glz*5/8) then
        tempr = 0.0
      else
        tempr = 1.0
      endif

      rhoR = tempr
      rhoB = 1.0 - tempr

      u = vx
      v = vy
      w = vz
      ! write(*,*) 'CUDA setup2',rho,i,j,k
  
      rhoR_d(i,j,k) = rhoR
      rhoB_d(i,j,k) = rhoB
      do l = 0, npops-1
	      eqR = equil(rhoR, u,v,w, l)
        popsR_d(i,j,k, l, 1) = eqR
	      eqB = equil(rhoB, u,v,w, l)
        popsB_d(i,j,k, l, 1) = eqB

        ! Clears other stuff
        popsR_d(i,j,k, l, 2) = 0.0
        popsB_d(i,j,k, l, 2) = 0.0
      end do

      myfluid_d(i,j,k, 1) = fluid_fluid
      myfluid_d(i,j,k, 2) = fluid_fluid
      debugfluid_d(i,j,k) = 0

      countmk_d = 0
      countrm_d = 0
      countn2p_d = 0
      countp2n1_d = 0
      countp2n2_d = 0
    end subroutine setup2_zplanes
    
    attributes(global) subroutine setup2_xplanes(vx,vy,vz)
      real, value :: vx,vy,vz
      integer :: i,j,k, l
      integer :: gli,glj,glk,hx
      real    :: rhoR,rhoB, u,v,w, eqR,eqB, tempr
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z

      gli = i + offset_d(1)
      glj = j + offset_d(2)
      glk = k + offset_d(3)

      if (1==i*j*k+offset_d(3)) write(*,*) 'CUDA setup2_zplanes] offset:', gli, glj, glk
      hx = nint(real(glx)*0.5)
      ! B @ 1/8-> 5/8 x-planes
      if (gli>hx) then
        tempr = 0.0
      else
        tempr = 1.0
      endif

      rhoR = tempr
      rhoB = 1.0 - tempr

      u = vx
      v = vy
      w = vz
      ! write(*,*) 'CUDA setup2',rho,i,j,k
  
      rhoR_d(i,j,k) = rhoR
      rhoB_d(i,j,k) = rhoB
      do l = 0, npops-1
	      eqR = equil(rhoR, u,v,w, l)
        popsR_d(i,j,k, l, 1) = eqR
	      eqB = equil(rhoB, u,v,w, l)
        popsB_d(i,j,k, l, 1) = eqB

        ! Clears other stuff
        popsR_d(i,j,k, l, 2) = 0.0
        popsB_d(i,j,k, l, 2) = 0.0
      end do
#ifndef APPLYBC
      myfluid_d(i,j,k, 1) = fluid_fluid
      myfluid_d(i,j,k, 2) = fluid_fluid
#endif
      debugfluid_d(i,j,k) = 0

      countmk_d = 0
      countrm_d = 0
      countn2p_d = 0
      countp2n1_d = 0
      countp2n2_d = 0
    end subroutine setup2_xplanes

    attributes(global) subroutine setupCyl(vx,vy,vz)
      real, value :: vx,vy,vz
      integer :: i,j,k, l
      integer :: gli,glj,glk
      real    :: rhoR,rhoB, u,v,w, eqR,eqB, distx,disty,distz, rdist,tempr
      real,parameter    :: radius = glz/8
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z

      gli = i + offset_d(1)
      glj = j + offset_d(2)
      glk = k + offset_d(3)

      if (1==i*j*k+offset_d(3)) write(*,*) 'CUDA setup2] offset:', gli, glj, glk

      distx = gli - (glx/2 + 0.5)
      disty = glj - (gly/2 + 0.5)
      distz = glk - (glz/4 + 0.5)

      rdist = sqrt(distx*distx + disty*disty)
      tempr = fcut(rdist, radius, radius+0.1)

      rhoR = tempr
      rhoB = 1.0 - tempr

      u = vx
      v = vy
      w = vz
  
      rhoR_d(i,j,k) = rhoR
      rhoB_d(i,j,k) = rhoB
      do l = 0, npops-1
	      eqR = equil(rhoR, u,v,w, l)
        popsR_d(i,j,k, l, 1) = eqR
	      eqB = equil(rhoB, u,v,w, l)
        popsB_d(i,j,k, l, 1) = eqB

        ! Clears other stuff
        popsR_d(i,j,k, l, 2) = 0.0
        popsB_d(i,j,k, l, 2) = 0.0
      end do

      myfluid_d(i,j,k, 1) = fluid_fluid
      myfluid_d(i,j,k, 2) = fluid_fluid
      debugfluid_d(i,j,k) = 0

      countmk_d = 0
      countrm_d = 0
      countn2p_d = 0
      countp2n1_d = 0
      countp2n2_d = 0
    end subroutine setupCyl

    attributes(global) subroutine setupPops(u,v,w)
      real, value :: u,v,w
      integer :: i,j,k, l
      integer :: gli,glj,glk
      real    :: rhoR,rhoB, eqR,eqB

  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z

      rhoR = rhoR_d(i,j,k)
      rhoB = rhoB_d(i,j,k)
      do l = 0, npops-1
	      eqR = equil(rhoR, u,v,w, l)
        popsR_d(i,j,k, l, 1) = eqR
	      eqB = equil(rhoB, u,v,w, l)
        popsB_d(i,j,k, l, 1) = eqB

        ! Clears other stuff
        popsR_d(i,j,k, l, 2) = 0.0
        popsB_d(i,j,k, l, 2) = 0.0
      end do

      myfluid_d(i,j,k, 1) = fluid_fluid
      myfluid_d(i,j,k, 2) = fluid_fluid
      debugfluid_d(i,j,k) = 0

      countmk_d = 0
      countrm_d = 0
      countn2p_d = 0
      countp2n1_d = 0
      countp2n2_d = 0
    end subroutine setupPops

    !!!!!!!!!!!!!!!!!!! Color gradient !!!!!!!!!!!!!!!!!!!
#ifdef D3Q27
#define manual_rho(pops) pops(i,j,k,0,flip) + pops(i,j,k,1,flip) + pops(i,j,k,2, flip) +   \
       pops(i,j,k,3, flip) + pops(i,j,k,4, flip) + pops(i,j,k,5, flip) + pops(i,j,k,6, flip) + pops(i,j,k,7, flip) +                     \
       pops(i,j,k,8, flip) + pops(i,j,k,9, flip) + pops(i,j,k,10, flip) + pops(i,j,k,11, flip) + pops(i,j,k,12, flip) +                  \
       pops(i,j,k,13, flip) + pops(i,j,k,14, flip) + pops(i,j,k,15, flip) + pops(i,j,k,16, flip) + pops(i,j,k,17, flip) +                \
       pops(i,j,k,18, flip) + pops(i,j,k,19, flip) + pops(i,j,k,20, flip) + pops(i,j,k,21, flip) + pops(i,j,k,22, flip) +                \
       pops(i,j,k,23, flip) + pops(i,j,k,24, flip) + pops(i,j,k,25, flip) + pops(i,j,k,26, flip)
#else
#define manual_rho(pops) pops(i,j,k,0,flip) + pops(i,j,k,1,flip) + pops(i,j,k,2, flip) +   \
       pops(i,j,k,3, flip) + pops(i,j,k,4, flip) + pops(i,j,k,5, flip) + pops(i,j,k,6, flip) + pops(i,j,k,7, flip) +                     \
       pops(i,j,k,8, flip) + pops(i,j,k,9, flip) + pops(i,j,k,10, flip) + pops(i,j,k,11, flip) + pops(i,j,k,12, flip) +                  \
       pops(i,j,k,13, flip) + pops(i,j,k,14, flip) + pops(i,j,k,15, flip) + pops(i,j,k,16, flip) + pops(i,j,k,17, flip) + pops(i,j,k,18, flip)
#endif
    attributes(global) subroutine init_rho_CG(step, flip)
      integer, value :: step,flip
      real    :: rhoR,rhoB
      integer :: i,j,k
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z

      if (step==1 .and. 1==i*j*k+offset_d(3)) write(*,*) 'CUDA init_rho_CG]'

      ! Only on fluid nodes
      if (myfluid_d(i,j,k, flip) == fluid_fluid) then
        rhoR = manual_rho(popsR_d)		! compute_rho(popsR_d,i,j,k,flip)
        rhoB = manual_rho(popsB_d)		! compute_rho(popsB_d,i,j,k,flip)

        if (rhoR<minRho .or. rhoR>maxRho) then
          write(*,*) 'init_rho_CG]Range error rhoR', step, linear(i,j,k), rhoR
          stop_d = __LINE__
        endif
        if (rhoB<minRho .or. rhoB>maxRho) then
          write(*,*) 'init_rho_CG]Range error rhoB', step, linear(i,j,k), rhoB
          stop_d = __LINE__
        endif

        if (rhoR < MINDENS) rhoR = MINDENS
        if (rhoB < MINDENS) rhoB = MINDENS

        rhoR_d(i,j,k) = rhoR
        rhoB_d(i,j,k) = rhoB
      else
        rhoR_d(i,j,k) = MINDENS
        rhoB_d(i,j,k) = MINDENS
      endif
    end subroutine init_rho_CG
    
    attributes(global) subroutine init_rhoR_CG(step, flip)
      integer, value :: step,flip
      integer :: i,j,k
      real    :: rhoR

      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z

      if (step==1 .and. 1==i*j*k+offset_d(3)) write(*,*) 'CUDA init_rhoR_CG]'

      ! Only on fluid nodes
      if (myfluid_d(i,j,k, flip) == fluid_fluid) then
#ifdef D3Q27
      rhoR = popsR_d(i,j,k,0, flip) + popsR_d(i,j,k,1, flip) + &
               popsR_d(i,j,k,2, flip) + popsR_d(i,j,k,3, flip) + &
               popsR_d(i,j,k,4, flip) + popsR_d(i,j,k,5, flip) + &
               popsR_d(i,j,k,6, flip) + popsR_d(i,j,k,7, flip) + &
               popsR_d(i,j,k,8, flip) + popsR_d(i,j,k,9, flip) + &
               popsR_d(i,j,k,10,flip) + popsR_d(i,j,k,11,flip) + &
               popsR_d(i,j,k,12,flip) + popsR_d(i,j,k,13,flip) + &
               popsR_d(i,j,k,14,flip) + popsR_d(i,j,k,15,flip) + &
               popsR_d(i,j,k,16,flip) + popsR_d(i,j,k,17,flip) + &
               popsR_d(i,j,k,18,flip) + popsR_d(i,j,k,19,flip) + &
               popsR_d(i,j,k,20,flip) + popsR_d(i,j,k,21,flip) + &
               popsR_d(i,j,k,22,flip) + popsR_d(i,j,k,23,flip) + &
               popsR_d(i,j,k,24,flip) + popsR_d(i,j,k,25,flip) + &
               popsR_d(i,j,k,26,flip)
#else
        rhoR = popsR_d(i,j,k,0, flip) + popsR_d(i,j,k,1, flip) + &
               popsR_d(i,j,k,2, flip) + popsR_d(i,j,k,3, flip) + &
               popsR_d(i,j,k,4, flip) + popsR_d(i,j,k,5, flip) + &
               popsR_d(i,j,k,6, flip) + popsR_d(i,j,k,7, flip) + &
               popsR_d(i,j,k,8, flip) + popsR_d(i,j,k,9, flip) + &
               popsR_d(i,j,k,10,flip) + popsR_d(i,j,k,11,flip) + &
               popsR_d(i,j,k,12,flip) + popsR_d(i,j,k,13,flip) + &
               popsR_d(i,j,k,14,flip) + popsR_d(i,j,k,15,flip) + &
               popsR_d(i,j,k,16,flip) + popsR_d(i,j,k,17,flip) + &
               popsR_d(i,j,k,18,flip)
#endif
        if (rhoR<minRho .or. rhoR>maxRho) then
          write(*,*) 'init_rhoR_CG]Range error rhoR', step, linear(i,j,k), rhoR
          stop_d = __LINE__
        endif

        if (rhoR < MINDENS) then
          ! write(*,*) 'init_rhoR_CG] rhoR fixed from', step, linear(i,j,k), rhoR
          rhoR = MINDENS
        endif

        rhoR_d(i,j,k) = rhoR
      else
        if(fixdenswall)then
          rhoR_d(i,j,k) = denswallR_d
        else
          rhoR_d(i,j,k) = MINDENS
        endif
        if(store_vel_d)then
          vel_d(1,i,j,k) = 0.0
          vel_d(2,i,j,k) = 0.0
          vel_d(3,i,j,k) = 0.0
        endif
      endif

    end subroutine init_rhoR_CG
    
    attributes(global) subroutine init_rhoB_CG(step, flip)
      integer, value :: step,flip
      integer :: i,j,k
      real    :: rhoB

      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z

      if (step==1 .and. 1==i*j*k+offset_d(3)) write(*,*) 'CUDA init_rhoB_CG]'

      ! Only on fluid nodes
      if (myfluid_d(i,j,k, flip) == fluid_fluid) then
#ifdef D3Q27
        rhoB = popsB_d(i,j,k,0, flip) + popsB_d(i,j,k,1, flip) + &
               popsB_d(i,j,k,2, flip) + popsB_d(i,j,k,3, flip) + &
               popsB_d(i,j,k,4, flip) + popsB_d(i,j,k,5, flip) + &
               popsB_d(i,j,k,6, flip) + popsB_d(i,j,k,7, flip) + &
               popsB_d(i,j,k,8, flip) + popsB_d(i,j,k,9, flip) + &
               popsB_d(i,j,k,10,flip) + popsB_d(i,j,k,11,flip) + &
               popsB_d(i,j,k,12,flip) + popsB_d(i,j,k,13,flip) + &
               popsB_d(i,j,k,14,flip) + popsB_d(i,j,k,15,flip) + &
               popsB_d(i,j,k,16,flip) + popsB_d(i,j,k,17,flip) + &
               popsB_d(i,j,k,18,flip) + popsB_d(i,j,k,19,flip) + &
               popsB_d(i,j,k,20,flip) + popsB_d(i,j,k,21,flip) + &
               popsB_d(i,j,k,22,flip) + popsB_d(i,j,k,23,flip) + &
               popsB_d(i,j,k,24,flip) + popsB_d(i,j,k,25,flip) + &
               popsB_d(i,j,k,26,flip)
#else
        rhoB = popsB_d(i,j,k,0, flip) + popsB_d(i,j,k,1, flip) + &
               popsB_d(i,j,k,2, flip) + popsB_d(i,j,k,3, flip) + &
               popsB_d(i,j,k,4, flip) + popsB_d(i,j,k,5, flip) + &
               popsB_d(i,j,k,6, flip) + popsB_d(i,j,k,7, flip) + &
               popsB_d(i,j,k,8, flip) + popsB_d(i,j,k,9, flip) + &
               popsB_d(i,j,k,10,flip) + popsB_d(i,j,k,11,flip) + &
               popsB_d(i,j,k,12,flip) + popsB_d(i,j,k,13,flip) + &
               popsB_d(i,j,k,14,flip) + popsB_d(i,j,k,15,flip) + &
               popsB_d(i,j,k,16,flip) + popsB_d(i,j,k,17,flip) + &
               popsB_d(i,j,k,18,flip)
#endif
        if (rhoB<minRho .or. rhoB>maxRho) then
          write(*,*) 'init_rhoB_CG]Range error rhoB', step, linear(i,j,k), rhoB
          stop_d = __LINE__
        endif

        if (rhoB < MINDENS) then
          ! write(*,*) 'init_rhoB_CG] rhoB fixed from', step, linear(i,j,k), rhoB
          rhoB = MINDENS
        endif
        rhoB_d(i,j,k) = rhoB
      else
        if(fixdenswall)then
          rhoB_d(i,j,k) = denswallB_d
        else
          rhoB_d(i,j,k) = MINDENS
        endif
      endif

    end subroutine init_rhoB_CG
    
    attributes(global) subroutine init_nearsel_CG(step, flip)
      integer, value :: step,flip
      integer :: i,j,k
      real    :: phase_field

      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z

      if (step==1 .and. 1==i*j*k+offset_d(3)) write(*,*) 'CUDA init_rhoR_CG]'

      ! Only on fluid nodes
      if (myfluid_d(i,j,k, flip) == fluid_fluid) then
#ifdef DIFFDENS
        phase_field=(rhoR_d(i,j,k)/densR_d-rhoB_d(i,j,k)/densB_d)/ &
         (rhoR_d(i,j,k)/densR_d+rhoB_d(i,j,k)/densB_d)
#else
        phase_field=(rhoR_d(i,j,k)-rhoB_d(i,j,k))/(rhoR_d(i,j,k)+rhoB_d(i,j,k))
#endif
        if(abs(phase_field)<0.5)then
          nearsel_d(i,j,k) = 1
        else
          nearsel_d(i,j,k) = 0
        endif
      else
        nearsel_d(i,j,k) = 0
      endif

    end subroutine init_nearsel_CG
    
    attributes(global) subroutine init_gradnearsel_CG(step, flip)
      integer, value :: step,flip
      integer :: i,j,k,li,lj,lk,l,i1,j1,k1
      real,shared    :: loc_rhoR(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real,shared    :: loc_rhoB(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real :: rhoR_shifted,rhoB_shifted,psi_shifted,phase_field
      real :: psix,psiy,psiz,psinorm,psinorm_sq

      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z

      if (step==1 .and. 1==i*j*k+offset_d(3)) write(*,*) 'CUDA init_rhoR_CG]'
      
      li = threadIdx%x
      lj = threadIdx%y
      lk = threadIdx%z

      loc_rhoR(li,lj,lk) = rhoR_d(i,j,k)

      ! Halo Faces
      if (li==1) then
	loc_rhoR(li-1,lj,lk) = rhoR_d(i-1,j,k)
      endif
      if (li==TILE_DIMx) then
	loc_rhoR(li+1,lj,lk) = rhoR_d(i+1,j,k)
      endif

      if (lj==1) then
	loc_rhoR(li,lj-1,lk) = rhoR_d(i,j-1,k)
      endif
      if (lj==TILE_DIMy) then
	loc_rhoR(li,lj+1,lk) = rhoR_d(i,j+1,k)
      endif

      if (lk==1) then
	loc_rhoR(li,lj,lk-1) = rhoR_d(i,j,k-1)
      endif
      if (lk==TILE_DIMz) then
	loc_rhoR(li,lj,lk+1) = rhoR_d(i,j,k+1)
      endif

      ! Halo edges
      if (li==1 .and. lj==1) loc_rhoR(li-1,lj-1,lk) = rhoR_d(i-1,j-1,k)
      if (li==1 .and. lj==TILE_DIMy) loc_rhoR(li-1,lj+1,lk) = rhoR_d(i-1,j+1,k)

      if (li==1 .and. lk==1) loc_rhoR(li-1,lj,lk-1) = rhoR_d(i-1,j,k-1)
      if (li==1 .and. lk==TILE_DIMz) loc_rhoR(li-1,lj,lk+1) = rhoR_d(i-1,j,k+1)

      if (li==TILE_DIMx .and. lj==1) loc_rhoR(li+1,lj-1,lk) = rhoR_d(i+1,j-1,k)
      if (li==TILE_DIMx .and. lj==TILE_DIMy) loc_rhoR(li+1,lj+1,lk) = rhoR_d(i+1,j+1,k)

      if (li==TILE_DIMx .and. lk==1) loc_rhoR(li+1,lj,lk-1) = rhoR_d(i+1,j,k-1)
      if (li==TILE_DIMx .and. lk==TILE_DIMz) loc_rhoR(li+1,lj,lk+1) = rhoR_d(i+1,j,k+1)

      if (lj==1 .and. lk==1) loc_rhoR(li,lj-1,lk-1) = rhoR_d(i,j-1,k-1)
      if (lj==1 .and. lk==TILE_DIMz) loc_rhoR(li,lj-1,lk+1) = rhoR_d(i,j-1,k+1)

      if (lj==TILE_DIMy .and. lk==1) loc_rhoR(li,lj+1,lk-1) = rhoR_d(i,j+1,k-1)
      if (lj==TILE_DIMy .and. lk==TILE_DIMz) loc_rhoR(li,lj+1,lk+1) = rhoR_d(i,j+1,k+1)
             
      ! Halo corner
      if (li==1 .and. lj==1 .and. lk==1) loc_rhoR(li-1,lj-1,lk-1) = rhoR_d(i-1,j-1,k-1) 
      if (li==TILE_DIMx .and. lj==1 .and. lk==1) loc_rhoR(li+1,lj-1,lk-1) = rhoR_d(i+1,j-1,k-1)
      if (li==1 .and. lj==TILE_DIMy .and. lk==1) loc_rhoR(li-1,lj+1,lk-1) = rhoR_d(i-1,j+1,k-1) 
      if (li==1 .and. lj==1 .and. lk==TILE_DIMz) loc_rhoR(li-1,lj-1,lk+1) = rhoR_d(i-1,j-1,k+1) 
      if (li==1 .and. lj==TILE_DIMy .and. lk==TILE_DIMz) loc_rhoR(li-1,lj+1,lk+1) = rhoR_d(i-1,j+1,k+1) 
      if (li==TILE_DIMx .and. lj==1 .and. lk==TILE_DIMz) loc_rhoR(li+1,lj-1,lk+1) = rhoR_d(i+1,j-1,k+1) 
      if (li==TILE_DIMx .and. lj==TILE_DIMy .and. lk==1) loc_rhoR(li+1,lj+1,lk-1) = rhoR_d(i+1,j+1,k-1) 
      if (li==TILE_DIMx .and. lj==TILE_DIMy .and. lk==TILE_DIMz) loc_rhoR(li+1,lj+1,lk+1) = rhoR_d(i+1,j+1,k+1) 
           
      ! Ripeto per il blue
      loc_rhoB(li,lj,lk) = rhoB_d(i,j,k)

      ! Halo Faces
      if (li==1) then
	loc_rhoB(li-1,lj,lk) = rhoB_d(i-1,j,k)
      endif
      if (li==TILE_DIMx) then
	loc_rhoB(li+1,lj,lk) = rhoB_d(i+1,j,k)
      endif

      if (lj==1) then
	loc_rhoB(li,lj-1,lk) = rhoB_d(i,j-1,k)
      endif
      if (lj==TILE_DIMy) then
	loc_rhoB(li,lj+1,lk) = rhoB_d(i,j+1,k)
      endif

      if (lk==1) then
	loc_rhoB(li,lj,lk-1) = rhoB_d(i,j,k-1)
      endif
      if (lk==TILE_DIMz) then
	loc_rhoB(li,lj,lk+1) = rhoB_d(i,j,k+1)
      endif

      ! Halo edges
      if (li==1 .and. lj==1) loc_rhoB(li-1,lj-1,lk) = rhoB_d(i-1,j-1,k)
      if (li==1 .and. lj==TILE_DIMy) loc_rhoB(li-1,lj+1,lk) = rhoB_d(i-1,j+1,k)

      if (li==1 .and. lk==1) loc_rhoB(li-1,lj,lk-1) = rhoB_d(i-1,j,k-1)
      if (li==1 .and. lk==TILE_DIMz) loc_rhoB(li-1,lj,lk+1) = rhoB_d(i-1,j,k+1)

      if (li==TILE_DIMx .and. lj==1) loc_rhoB(li+1,lj-1,lk) = rhoB_d(i+1,j-1,k)
      if (li==TILE_DIMx .and. lj==TILE_DIMy) loc_rhoB(li+1,lj+1,lk) = rhoB_d(i+1,j+1,k)

      if (li==TILE_DIMx .and. lk==1) loc_rhoB(li+1,lj,lk-1) = rhoB_d(i+1,j,k-1)
      if (li==TILE_DIMx .and. lk==TILE_DIMz) loc_rhoB(li+1,lj,lk+1) = rhoB_d(i+1,j,k+1)

      if (lj==1 .and. lk==1) loc_rhoB(li,lj-1,lk-1) = rhoB_d(i,j-1,k-1)
      if (lj==1 .and. lk==TILE_DIMz) loc_rhoB(li,lj-1,lk+1) = rhoB_d(i,j-1,k+1)

      if (lj==TILE_DIMy .and. lk==1) loc_rhoB(li,lj+1,lk-1) = rhoB_d(i,j+1,k-1)
      if (lj==TILE_DIMy .and. lk==TILE_DIMz) loc_rhoB(li,lj+1,lk+1) = rhoB_d(i,j+1,k+1)
         
      ! Halo corner
      if (li==1 .and. lj==1 .and. lk==1) loc_rhoB(li-1,lj-1,lk-1) = rhoB_d(i-1,j-1,k-1) 
      if (li==TILE_DIMx .and. lj==1 .and. lk==1) loc_rhoB(li+1,lj-1,lk-1) = rhoB_d(i+1,j-1,k-1)
      if (li==1 .and. lj==TILE_DIMy .and. lk==1) loc_rhoB(li-1,lj+1,lk-1) = rhoB_d(i-1,j+1,k-1) 
      if (li==1 .and. lj==1 .and. lk==TILE_DIMz) loc_rhoB(li-1,lj-1,lk+1) = rhoB_d(i-1,j-1,k+1) 
      if (li==1 .and. lj==TILE_DIMy .and. lk==TILE_DIMz) loc_rhoB(li-1,lj+1,lk+1) = rhoB_d(i-1,j+1,k+1) 
      if (li==TILE_DIMx .and. lj==1 .and. lk==TILE_DIMz) loc_rhoB(li+1,lj-1,lk+1) = rhoB_d(i+1,j-1,k+1) 
      if (li==TILE_DIMx .and. lj==TILE_DIMy .and. lk==1) loc_rhoB(li+1,lj+1,lk-1) = rhoB_d(i+1,j+1,k-1) 
      if (li==TILE_DIMx .and. lj==TILE_DIMy .and. lk==TILE_DIMz) loc_rhoB(li+1,lj+1,lk+1) = rhoB_d(i+1,j+1,k+1) 
      
      call syncthreads
      


      ! Only on fluid nodes
      if (myfluid_d(i,j,k, flip) == fluid_fluid) then
      
        ! Psi calc
        psix = ZERO
        psiy = ZERO
        psiz = ZERO
#ifdef HIGHGRAD
        do l = 1, linksd3q27-1
          i1 = li + exd3q27_d(l)
          j1 = lj + eyd3q27_d(l)
          k1 = lk + ezd3q27_d(l)

          rhoR_shifted = loc_rhoR(i1,j1,k1)
          rhoB_shifted = loc_rhoB(i1,j1,k1)
          psi_shifted = (rhoR_shifted/densR_d - rhoB_shifted/densB_d) / &
           (rhoR_shifted/densR_d + rhoB_shifted/densB_d)
          psix = psix + ad3q27_d(l)*exd3q27_d(l)* psi_shifted
          psiy = psiy + ad3q27_d(l)*eyd3q27_d(l)* psi_shifted
          psiz = psiz + ad3q27_d(l)*ezd3q27_d(l)* psi_shifted
        enddo
#else
        do l = 1, npops-1
          i1 = li + ex_d(l)
          j1 = lj + ey_d(l)
          k1 = lk + ez_d(l)

          rhoR_shifted = loc_rhoR(i1,j1,k1)
          rhoB_shifted = loc_rhoB(i1,j1,k1)
          
#ifdef DIFFDENS
          psi_shifted = (rhoR_shifted/densR_d - rhoB_shifted/densB_d) / &
           (rhoR_shifted/densR_d + rhoB_shifted/densB_d)
#else
          psi_shifted = (rhoR_shifted - rhoB_shifted) / (rhoR_shifted + rhoB_shifted)
#endif
          psix = psix + a_d(l)*ex_d(l)* psi_shifted
          psiy = psiy + a_d(l)*ey_d(l)* psi_shifted
          psiz = psiz + a_d(l)*ez_d(l)* psi_shifted
        enddo
#endif
        psinorm_sq = psix*psix + psiy*psiy + psiz*psiz
        psinorm = sqrt(psinorm_sq)
         
#ifdef DIFFDENS
        phase_field=(rhoR_d(i,j,k)/densR_d-rhoB_d(i,j,k)/densB_d)/ &
         (rhoR_d(i,j,k)/densR_d+rhoB_d(i,j,k)/densB_d)
#else
        phase_field=(rhoR_d(i,j,k)-rhoB_d(i,j,k))/(rhoR_d(i,j,k)+rhoB_d(i,j,k))
#endif
        if(psinorm>=0.55 .and. abs(phase_field)<0.5)then
          nearsel_d(i,j,k) = 1
        else
          nearsel_d(i,j,k) = 0
        endif
      else
        nearsel_d(i,j,k) = 0
      endif

    end subroutine init_gradnearsel_CG
    
    attributes(global) subroutine applybchalo_meandenswall(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l, flop, lopp, itemp, mytype
      integer :: i1,j1,k1
      integer :: i2,j2,k2
      integer :: li,lj,lk
      integer :: istat
      logical :: alltrue
      real :: rho1,rho2,rtemp
      real,shared    :: loc_rhoR(-1:TILE_DIMx+2,-1:TILE_DIMy+2,-1:TILE_DIMz+2)
      real,shared    :: loc_rhoB(-1:TILE_DIMx+2,-1:TILE_DIMy+2,-1:TILE_DIMz+2)
      
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x - 1
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y - 1
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z - 1
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA stream_BGK] i,k:',i,j,k
      
      flop = 3 - flip
      
      if (i>nx+1) return
      if (j>ny+1) return
      if (k>nz_d+1) return
      
      li = threadIdx%x-1
      lj = threadIdx%y-1
      lk = threadIdx%z-1
      
      alltrue = .false.
      if( allthreads(myfluid_d(i,j,k, flip)==1) )alltrue = .true.
      call syncthreads
      if(alltrue)return
      
      loc_rhoR(li,lj,lk) = rhoR_d(i,j,k)

      ! Halo Faces
      if (li==0) then
	    loc_rhoR(li-1,lj,lk) = rhoR_d(i-1,j,k)
      endif
      if (li==TILE_DIMx+1) then
	    loc_rhoR(li+1,lj,lk) = rhoR_d(i+1,j,k)
      endif

      if (lj==0) then
	    loc_rhoR(li,lj-1,lk) = rhoR_d(i,j-1,k)
      endif
      if (lj==TILE_DIMy+1) then
	    loc_rhoR(li,lj+1,lk) = rhoR_d(i,j+1,k)
      endif

      if (lk==0) then
	    loc_rhoR(li,lj,lk-1) = rhoR_d(i,j,k-1)
      endif
      if (lk==TILE_DIMz+1) then
	    loc_rhoR(li,lj,lk+1) = rhoR_d(i,j,k+1)
      endif

      ! Halo edges
      if (li==0 .and. lj==0) loc_rhoR(li-1,lj-1,lk) = rhoR_d(i-1,j-1,k)
      if (li==0 .and. lj==TILE_DIMy+1) loc_rhoR(li-1,lj+1,lk) = rhoR_d(i-1,j+1,k)

      if (li==0 .and. lk==0) loc_rhoR(li-1,lj,lk-1) = rhoR_d(i-1,j,k-1)
      if (li==0 .and. lk==TILE_DIMz+1) loc_rhoR(li-1,lj,lk+1) = rhoR_d(i-1,j,k+1)

      if (li==TILE_DIMx+1 .and. lj==0) loc_rhoR(li+1,lj-1,lk) = rhoR_d(i+1,j-1,k)
      if (li==TILE_DIMx+1 .and. lj==TILE_DIMy+1) loc_rhoR(li+1,lj+1,lk) = rhoR_d(i+1,j+1,k)

      if (li==TILE_DIMx+1 .and. lk==0) loc_rhoR(li+1,lj,lk-1) = rhoR_d(i+1,j,k-1)
      if (li==TILE_DIMx+1 .and. lk==TILE_DIMz+1) loc_rhoR(li+1,lj,lk+1) = rhoR_d(i+1,j,k+1)

      if (lj==0 .and. lk==0) loc_rhoR(li,lj-1,lk-1) = rhoR_d(i,j-1,k-1)
      if (lj==0 .and. lk==TILE_DIMz+1) loc_rhoR(li,lj-1,lk+1) = rhoR_d(i,j-1,k+1)

      if (lj==TILE_DIMy+1 .and. lk==0) loc_rhoR(li,lj+1,lk-1) = rhoR_d(i,j+1,k-1)
      if (lj==TILE_DIMy+1 .and. lk==TILE_DIMz+1) loc_rhoR(li,lj+1,lk+1) = rhoR_d(i,j+1,k+1)
             
      ! Halo corner
      if (li==0 .and. lj==0 .and. lk==0) loc_rhoR(li-1,lj-1,lk-1) = rhoR_d(i-1,j-1,k-1) 
      if (li==TILE_DIMx+1 .and. lj==0 .and. lk==0) loc_rhoR(li+1,lj-1,lk-1) = rhoR_d(i+1,j-1,k-1)
      if (li==0 .and. lj==TILE_DIMy+1 .and. lk==0) loc_rhoR(li-1,lj+1,lk-1) = rhoR_d(i-1,j+1,k-1) 
      if (li==0 .and. lj==0 .and. lk==TILE_DIMz+1) loc_rhoR(li-1,lj-1,lk+1) = rhoR_d(i-1,j-1,k+1) 
      if (li==0 .and. lj==TILE_DIMy+1 .and. lk==TILE_DIMz+1) loc_rhoR(li-1,lj+1,lk+1) = rhoR_d(i-1,j+1,k+1) 
      if (li==TILE_DIMx+1 .and. lj==0 .and. lk==TILE_DIMz) loc_rhoR(li+1,lj-1,lk+1) = rhoR_d(i+1,j-1,k+1) 
      if (li==TILE_DIMx+1 .and. lj==TILE_DIMy+1 .and. lk==0) loc_rhoR(li+1,lj+1,lk-1) = rhoR_d(i+1,j+1,k-1) 
      if (li==TILE_DIMx+1 .and. lj==TILE_DIMy+1 .and. lk==TILE_DIMz+1) loc_rhoR(li+1,lj+1,lk+1) = rhoR_d(i+1,j+1,k+1) 
           
      ! Ripeto per il blue
      loc_rhoB(li,lj,lk) = rhoB_d(i,j,k)

      ! Halo Faces
      if (li==0) then
	    loc_rhoB(li-1,lj,lk) = rhoB_d(i-1,j,k)
      endif
      if (li==TILE_DIMx+1) then
	    loc_rhoB(li+1,lj,lk) = rhoB_d(i+1,j,k)
      endif

      if (lj==0) then
	    loc_rhoB(li,lj-1,lk) = rhoB_d(i,j-1,k)
      endif
      if (lj==TILE_DIMy+1) then
	    loc_rhoB(li,lj+1,lk) = rhoB_d(i,j+1,k)
      endif

      if (lk==0) then
	    loc_rhoB(li,lj,lk-1) = rhoB_d(i,j,k-1)
      endif
      if (lk==TILE_DIMz+1) then
	    loc_rhoB(li,lj,lk+1) = rhoB_d(i,j,k+1)
      endif

      ! Halo edges
      if (li==0 .and. lj==0) loc_rhoB(li-1,lj-1,lk) = rhoB_d(i-1,j-1,k)
      if (li==0 .and. lj==TILE_DIMy+1) loc_rhoB(li-1,lj+1,lk) = rhoB_d(i-1,j+1,k)

      if (li==0 .and. lk==0) loc_rhoB(li-1,lj,lk-1) = rhoB_d(i-1,j,k-1)
      if (li==0 .and. lk==TILE_DIMz+1) loc_rhoB(li-1,lj,lk+1) = rhoB_d(i-1,j,k+1)

      if (li==TILE_DIMx+1 .and. lj==0) loc_rhoB(li+1,lj-1,lk) = rhoB_d(i+1,j-1,k)
      if (li==TILE_DIMx+1 .and. lj==TILE_DIMy+1) loc_rhoB(li+1,lj+1,lk) = rhoB_d(i+1,j+1,k)

      if (li==TILE_DIMx+1 .and. lk==0) loc_rhoB(li+1,lj,lk-1) = rhoB_d(i+1,j,k-1)
      if (li==TILE_DIMx+1 .and. lk==TILE_DIMz+1) loc_rhoB(li+1,lj,lk+1) = rhoB_d(i+1,j,k+1)

      if (lj==0 .and. lk==0) loc_rhoB(li,lj-1,lk-1) = rhoB_d(i,j-1,k-1)
      if (lj==0 .and. lk==TILE_DIMz+1) loc_rhoB(li,lj-1,lk+1) = rhoB_d(i,j-1,k+1)

      if (lj==TILE_DIMy+1 .and. lk==0) loc_rhoB(li,lj+1,lk-1) = rhoB_d(i,j+1,k-1)
      if (lj==TILE_DIMy+1 .and. lk==TILE_DIMz+1) loc_rhoB(li,lj+1,lk+1) = rhoB_d(i,j+1,k+1)
         
      ! Halo corner
      if (li==0 .and. lj==0 .and. lk==0) loc_rhoB(li-1,lj-1,lk-1) = rhoB_d(i-1,j-1,k-1) 
      if (li==TILE_DIMx+1 .and. lj==0 .and. lk==0) loc_rhoB(li+1,lj-1,lk-1) = rhoB_d(i+1,j-1,k-1)
      if (li==0 .and. lj==TILE_DIMy+1 .and. lk==0) loc_rhoB(li-1,lj+1,lk-1) = rhoB_d(i-1,j+1,k-1) 
      if (li==0 .and. lj==0 .and. lk==TILE_DIMz+1) loc_rhoB(li-1,lj-1,lk+1) = rhoB_d(i-1,j-1,k+1) 
      if (li==0 .and. lj==TILE_DIMy+1 .and. lk==TILE_DIMz+1) loc_rhoB(li-1,lj+1,lk+1) = rhoB_d(i-1,j+1,k+1) 
      if (li==TILE_DIMx+1 .and. lj==0 .and. lk==TILE_DIMz+1) loc_rhoB(li+1,lj-1,lk+1) = rhoB_d(i+1,j-1,k+1) 
      if (li==TILE_DIMx+1 .and. lj==TILE_DIMy+1 .and. lk==0) loc_rhoB(li+1,lj+1,lk-1) = rhoB_d(i+1,j+1,k-1) 
      if (li==TILE_DIMx+1 .and. lj==TILE_DIMy+1 .and. lk==TILE_DIMz+1) loc_rhoB(li+1,lj+1,lk+1) = rhoB_d(i+1,j+1,k+1) 
      
      call syncthreads
      
      
      ! compute_mean     
      if (myfluid_d(i,j,k, flip) <= fluid_wall) then
        rho1=0.0
        rho2=0.0
        rtemp=0.0
        do l = 1, npops-1
          i2 = i + ex_d(l)
          j2 = j + ey_d(l)
          k2 = k + ez_d(l)
          if (myfluid_d(i2,j2,k2, flip) <= fluid_fluid) then
            i1 = li + ex_d(l)
            j1 = lj + ey_d(l)
            k1 = lk + ez_d(l)
            rho1 = rho1 + a_d(l)* loc_rhoR(i1,j1,k1)
            rho2 = rho2 + a_d(l)* loc_rhoB(i1,j1,k1)
            rtemp = rtemp + a_d(l)
          endif
        enddo
        rhoR_d(i,j,k)= rho1/rtemp
        rhoB_d(i,j,k)= rho2/rtemp
      endif
      
      return
      
    end subroutine applybchalo_meandenswall
    
    attributes(global) subroutine applybc_meandenswall(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l, flop, lopp, itemp, mytype
      integer :: i1,j1,k1
      integer :: i2,j2,k2
      integer :: li,lj,lk
      integer :: istat
      integer :: itest=0
      integer,shared :: isf(1:TILE_DIMx,1:TILE_DIMy,1:TILE_DIMz)
      real :: rho1,rho2,rtemp
      real,shared    :: loc_rhoR(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real,shared    :: loc_rhoB(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA stream_BGK] i,k:',i,j,k
      
      flop = 3 - flip
      
      if (i>nx+1) return
      if (j>ny+1) return
      if (k>nz_d+1) return
      
      li = threadIdx%x
      lj = threadIdx%y
      lk = threadIdx%z
      
      isf(li,lj,lk) = 0
      
      if (myfluid_d(i,j,k, flip) <= fluid_wall) then
        isf(li,lj,lk) = 1
      endif
      
      
      istat=atomicor(itest,isf(li,lj,lk))
      
      call syncthreads
      
      if(itest==0)return
      
      loc_rhoR(li,lj,lk) = rhoR_d(i,j,k)

      ! Halo Faces
      if (li==1) then
	loc_rhoR(li-1,lj,lk) = rhoR_d(i-1,j,k)
      endif
      if (li==TILE_DIMx) then
	loc_rhoR(li+1,lj,lk) = rhoR_d(i+1,j,k)
      endif

      if (lj==1) then
	loc_rhoR(li,lj-1,lk) = rhoR_d(i,j-1,k)
      endif
      if (lj==TILE_DIMy) then
	loc_rhoR(li,lj+1,lk) = rhoR_d(i,j+1,k)
      endif

      if (lk==1) then
	loc_rhoR(li,lj,lk-1) = rhoR_d(i,j,k-1)
      endif
      if (lk==TILE_DIMz) then
	loc_rhoR(li,lj,lk+1) = rhoR_d(i,j,k+1)
      endif

      ! Halo edges
      if (li==1 .and. lj==1) loc_rhoR(li-1,lj-1,lk) = rhoR_d(i-1,j-1,k)
      if (li==1 .and. lj==TILE_DIMy) loc_rhoR(li-1,lj+1,lk) = rhoR_d(i-1,j+1,k)

      if (li==1 .and. lk==1) loc_rhoR(li-1,lj,lk-1) = rhoR_d(i-1,j,k-1)
      if (li==1 .and. lk==TILE_DIMz) loc_rhoR(li-1,lj,lk+1) = rhoR_d(i-1,j,k+1)

      if (li==TILE_DIMx .and. lj==1) loc_rhoR(li+1,lj-1,lk) = rhoR_d(i+1,j-1,k)
      if (li==TILE_DIMx .and. lj==TILE_DIMy) loc_rhoR(li+1,lj+1,lk) = rhoR_d(i+1,j+1,k)

      if (li==TILE_DIMx .and. lk==1) loc_rhoR(li+1,lj,lk-1) = rhoR_d(i+1,j,k-1)
      if (li==TILE_DIMx .and. lk==TILE_DIMz) loc_rhoR(li+1,lj,lk+1) = rhoR_d(i+1,j,k+1)

      if (lj==1 .and. lk==1) loc_rhoR(li,lj-1,lk-1) = rhoR_d(i,j-1,k-1)
      if (lj==1 .and. lk==TILE_DIMz) loc_rhoR(li,lj-1,lk+1) = rhoR_d(i,j-1,k+1)

      if (lj==TILE_DIMy .and. lk==1) loc_rhoR(li,lj+1,lk-1) = rhoR_d(i,j+1,k-1)
      if (lj==TILE_DIMy .and. lk==TILE_DIMz) loc_rhoR(li,lj+1,lk+1) = rhoR_d(i,j+1,k+1)
             
      ! Halo corner
      if (li==1 .and. lj==1 .and. lk==1) loc_rhoR(li-1,lj-1,lk-1) = rhoR_d(i-1,j-1,k-1) 
      if (li==TILE_DIMx .and. lj==1 .and. lk==1) loc_rhoR(li+1,lj-1,lk-1) = rhoR_d(i+1,j-1,k-1)
      if (li==1 .and. lj==TILE_DIMy .and. lk==1) loc_rhoR(li-1,lj+1,lk-1) = rhoR_d(i-1,j+1,k-1) 
      if (li==1 .and. lj==1 .and. lk==TILE_DIMz) loc_rhoR(li-1,lj-1,lk+1) = rhoR_d(i-1,j-1,k+1) 
      if (li==1 .and. lj==TILE_DIMy .and. lk==TILE_DIMz) loc_rhoR(li-1,lj+1,lk+1) = rhoR_d(i-1,j+1,k+1) 
      if (li==TILE_DIMx .and. lj==1 .and. lk==TILE_DIMz) loc_rhoR(li+1,lj-1,lk+1) = rhoR_d(i+1,j-1,k+1) 
      if (li==TILE_DIMx .and. lj==TILE_DIMy .and. lk==1) loc_rhoR(li+1,lj+1,lk-1) = rhoR_d(i+1,j+1,k-1) 
      if (li==TILE_DIMx .and. lj==TILE_DIMy .and. lk==TILE_DIMz) loc_rhoR(li+1,lj+1,lk+1) = rhoR_d(i+1,j+1,k+1) 
           
      ! Ripeto per il blue
      loc_rhoB(li,lj,lk) = rhoB_d(i,j,k)

      ! Halo Faces
      if (li==1) then
	loc_rhoB(li-1,lj,lk) = rhoB_d(i-1,j,k)
      endif
      if (li==TILE_DIMx) then
	loc_rhoB(li+1,lj,lk) = rhoB_d(i+1,j,k)
      endif

      if (lj==1) then
	loc_rhoB(li,lj-1,lk) = rhoB_d(i,j-1,k)
      endif
      if (lj==TILE_DIMy) then
	loc_rhoB(li,lj+1,lk) = rhoB_d(i,j+1,k)
      endif

      if (lk==1) then
	loc_rhoB(li,lj,lk-1) = rhoB_d(i,j,k-1)
      endif
      if (lk==TILE_DIMz) then
	loc_rhoB(li,lj,lk+1) = rhoB_d(i,j,k+1)
      endif

      ! Halo edges
      if (li==1 .and. lj==1) loc_rhoB(li-1,lj-1,lk) = rhoB_d(i-1,j-1,k)
      if (li==1 .and. lj==TILE_DIMy) loc_rhoB(li-1,lj+1,lk) = rhoB_d(i-1,j+1,k)

      if (li==1 .and. lk==1) loc_rhoB(li-1,lj,lk-1) = rhoB_d(i-1,j,k-1)
      if (li==1 .and. lk==TILE_DIMz) loc_rhoB(li-1,lj,lk+1) = rhoB_d(i-1,j,k+1)

      if (li==TILE_DIMx .and. lj==1) loc_rhoB(li+1,lj-1,lk) = rhoB_d(i+1,j-1,k)
      if (li==TILE_DIMx .and. lj==TILE_DIMy) loc_rhoB(li+1,lj+1,lk) = rhoB_d(i+1,j+1,k)

      if (li==TILE_DIMx .and. lk==1) loc_rhoB(li+1,lj,lk-1) = rhoB_d(i+1,j,k-1)
      if (li==TILE_DIMx .and. lk==TILE_DIMz) loc_rhoB(li+1,lj,lk+1) = rhoB_d(i+1,j,k+1)

      if (lj==1 .and. lk==1) loc_rhoB(li,lj-1,lk-1) = rhoB_d(i,j-1,k-1)
      if (lj==1 .and. lk==TILE_DIMz) loc_rhoB(li,lj-1,lk+1) = rhoB_d(i,j-1,k+1)

      if (lj==TILE_DIMy .and. lk==1) loc_rhoB(li,lj+1,lk-1) = rhoB_d(i,j+1,k-1)
      if (lj==TILE_DIMy .and. lk==TILE_DIMz) loc_rhoB(li,lj+1,lk+1) = rhoB_d(i,j+1,k+1)
         
      ! Halo corner
      if (li==1 .and. lj==1 .and. lk==1) loc_rhoB(li-1,lj-1,lk-1) = rhoB_d(i-1,j-1,k-1) 
      if (li==TILE_DIMx .and. lj==1 .and. lk==1) loc_rhoB(li+1,lj-1,lk-1) = rhoB_d(i+1,j-1,k-1)
      if (li==1 .and. lj==TILE_DIMy .and. lk==1) loc_rhoB(li-1,lj+1,lk-1) = rhoB_d(i-1,j+1,k-1) 
      if (li==1 .and. lj==1 .and. lk==TILE_DIMz) loc_rhoB(li-1,lj-1,lk+1) = rhoB_d(i-1,j-1,k+1) 
      if (li==1 .and. lj==TILE_DIMy .and. lk==TILE_DIMz) loc_rhoB(li-1,lj+1,lk+1) = rhoB_d(i-1,j+1,k+1) 
      if (li==TILE_DIMx .and. lj==1 .and. lk==TILE_DIMz) loc_rhoB(li+1,lj-1,lk+1) = rhoB_d(i+1,j-1,k+1) 
      if (li==TILE_DIMx .and. lj==TILE_DIMy .and. lk==1) loc_rhoB(li+1,lj+1,lk-1) = rhoB_d(i+1,j+1,k-1) 
      if (li==TILE_DIMx .and. lj==TILE_DIMy .and. lk==TILE_DIMz) loc_rhoB(li+1,lj+1,lk+1) = rhoB_d(i+1,j+1,k+1) 
      
      
      call syncthreads
      
      
      ! compute_mean     
      if (myfluid_d(i,j,k, flip) <= fluid_wall) then
        rho1=0.0
        rho2=0.0
        rtemp=0.0
        do l = 1, npops-1
          i2 = i + ex_d(l)
          j2 = j + ey_d(l)
          k2 = k + ez_d(l)
          if (myfluid_d(i2,j2,k2, flip) <= fluid_fluid) then
            i1 = li + ex_d(l)
            j1 = lj + ey_d(l)
            k1 = lk + ez_d(l)
            rho1 = rho1 + a_d(l)* loc_rhoR(i1,j1,k1)
            rho2 = rho2 + a_d(l)* loc_rhoB(i1,j1,k1)
            rtemp = rtemp + a_d(l)
          endif
        enddo
        rhoR_d(i,j,k)= rho1/rtemp
        rhoB_d(i,j,k)= rho2/rtemp
      endif
      
      return
      
    end subroutine applybc_meandenswall
    
    attributes(global) subroutine applybc_meandenswall_x(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l, flop
      integer :: i1,j1,k1
      real :: rho1,rho2,rtemp
      
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1

      if (j>ny+1) return
      if (k>nz_d+1) return

      flop = 3 - flip

      ! do x=0
      i = 0
      ! compute_mean     
      call meandenswall(i,j,k,flip)

      ! do x=nx+1
      i = nx+1
      ! compute_mean     
      call meandenswall(i,j,k,flip)
      
    end subroutine applybc_meandenswall_x

    attributes(global) subroutine applybc_meandenswall_y(step, flip)
    integer, value :: step,flip
    integer :: i,j,k, l, flop
    integer :: i1,j1,k1
    real :: rho1,rho2,rtemp

    i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
    k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
    
    if (i>nx+1) return
    if (k>nz_d+1) return

    flop = 3 - flip

    ! do y=0
    j = 0
    ! compute_mean     
    call meandenswall(i,j,k,flip)

    ! do y=ny+1
    j = ny+1
    ! compute_mean     
    call meandenswall(i,j,k,flip)
    
  end subroutine applybc_meandenswall_y

  attributes(global) subroutine applybc_meandenswall_z(step, flip)
    integer, value :: step,flip
    integer :: i,j,k, l, flop
    integer :: i1,j1,k1
    real :: rho1,rho2,rtemp

    i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
    j = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
    
    if (i>nx+1) return
    if (j>ny+1) return

    flop = 3 - flip

    ! do z=0
    k = 0
    ! compute_mean     
    call meandenswall(i,j,k,flip)

    ! do z=nz_d+1
    k = nz_d+1
    ! compute_mean     
    call meandenswall(i,j,k,flip)
    
  end subroutine applybc_meandenswall_z
  
  attributes(device) subroutine meandenswall(i,j,k,flip)
    
    implicit none
    
    integer, intent(in) :: i,j,k,flip
    integer :: l,i1,j1,k1
    real :: rho1,rho2,rtemp
    
     if (myfluid_d(i,j,k, flip) <= fluid_wall) then
        rho1=0.0
        rho2=0.0
        rtemp=0.0
        do l = 1, npops-1
          i1 = i + ex_d(l)
          j1 = j + ey_d(l)
          k1 = k + ez_d(l)
          if (myfluid_d(i1,j1,k1, flip) <= fluid_fluid) then
            rho1 = rho1 + a_d(l)* rhoR_d(i1,j1,k1)
            rho2 = rho2 + a_d(l)* rhoB_d(i1,j1,k1)
            rtemp = rtemp + a_d(l)
          endif
        enddo
        rhoR_d(i,j,k)= rho1/rtemp
        rhoB_d(i,j,k)= rho2/rtemp
      endif
      
    end subroutine meandenswall
    
    attributes(global) subroutine applybchalo_fixdenswall(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l, flop, lopp, itemp, mytype
      integer :: i1,j1,k1
      integer :: i2,j2,k2
      integer :: li,lj,lk
      integer :: istat
      integer :: itest=0
      
      real :: rho1,rho2,rtemp
      
      
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x - 1
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y - 1
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z - 1
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA stream_BGK] i,k:',i,j,k
      
      flop = 3 - flip
      
      if (i>nx+1) return
      if (j>ny+1) return
      if (k>nz_d+1) return
      
      ! set_fix   
      if (myfluid_d(i,j,k, flip) <= fluid_wall) then
        rhoR_d(i,j,k)= denswallR_d
        rhoB_d(i,j,k)= denswallB_d
      endif
      
      return
      
    end subroutine applybchalo_fixdenswall
    
#ifdef D3Q27
#define manual_u(pops)   invrho * ( pops(1) - pops(2) + pops(7) - pops(8) - pops(9) + \
       pops(10) + pops(11) - pops(12) - pops(13) + pops(14) + pops(19) - pops(20) + \
       pops(21) - pops(22) - pops(23) + pops(24) + pops(25) - pops(26) )

#define manual_v(pops)   invrho * ( pops(3) - pops(4) + pops(7) - pops(8) + pops(9) - \
       pops(10) + pops(15) - pops(16) - pops(17) + pops(18) + pops(19) - pops(20) + \
       pops(21) - pops(22) + pops(23) - pops(24) - pops(25) + pops(26) )
#define manual_w(pops)   invrho * ( pops(5) - pops(6) + pops(11) - pops(12) + pops(13) - \
       pops(14) + pops(15) - pops(16) + pops(17) - pops(18) + pops(19) - pops(20) - \
       pops(21) + pops(22) + pops(23) - pops(24) + pops(25) - pops(26) )
#else
#define manual_u(pops)   invrho * ( pops(1) - pops(2) + pops(7) - pops(8) - pops(9) + \
       pops(10) + pops(11) - pops(12) - pops(13) + pops(14) )

#define manual_v(pops)   invrho * ( pops(3) - pops(4) + pops(7) - pops(8) + pops(9) - \
       pops(10) + pops(15) - pops(16) - pops(17) + pops(18) )

#define manual_w(pops)   invrho * ( pops(5) - pops(6) + pops(11) - pops(12) + pops(13) - \
       pops(14) + pops(15) - pops(16) + pops(17) - pops(18) )
#endif
    attributes(global) subroutine time_step_CG(step, flip)
      integer, value :: step,flip
      real    :: rhoR,rhoB,rhosum,invrho, u,v,w
      real    :: viscavg,omega,oneminusomega
      real    :: rhoR_shifted,rhoB_shifted,rhosum_shifted
!      real    :: grad_rhoRx,grad_rhoRy,grad_rhoRz, grad_rhoBx,grad_rhoBy,grad_rhoBz
      real    :: psi_shifted,phase_field
      real    :: psix,psiy,psiz,psinorm_sq,psinorm,acoeff,e_dot_psi,myrhoR
#ifdef DIFFDENS
      real :: v_psi_R_dot,v_psi_B_dot,tens_contr
      real, dimension(3) :: psi_R,psi_B,myvel
      real, dimension(3,3) :: mat_R,mat_B,mat_TR,mat_TB,myG_R,myG_B
#endif
      real    :: temp,temp1, cosphi, feq,dnx,dny,dnz
      real    :: FRx,FRy,FRz,FBx,FBy,FBz,dist_rep,ushR,vshR,wshR,ushB,vshB,wshB
      real, parameter :: mylimit=1.e-20
      
!      real,shared    :: loc_tempR(TILE_DIMx,TILE_DIMy,TILE_DIMz,0:npops-1)
!      real,shared    :: loc_tempB(TILE_DIMx,TILE_DIMy,TILE_DIMz,0:npops-1)
!      real,shared    :: loc_fdum(TILE_DIMx,TILE_DIMy,TILE_DIMz,0:npops-1)

!      real,shared    :: loc_popsR(TILE_DIMx,TILE_DIMy,TILE_DIMz, 0:18)
!      real,shared    :: loc_popsB(TILE_DIMx,TILE_DIMy,TILE_DIMz, 0:18)
      real,shared    :: loc_rhoR(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real,shared    :: loc_rhoB(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)

#ifdef D3Q27
      real :: loc_tempR(0:26),loc_tempB(0:26),loc_fdum(0:26)
#else
      real :: loc_tempR(0:18),loc_tempB(0:18),loc_fdum(0:18)
#endif
      integer :: i,j,k, i1,j1,k1,l
      integer :: li,lj,lk
      integer :: gli,glj,glk
      integer :: inx,iny,inz,llx,lly,llz,ss,istatR,istatB

      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z

      if (step==1 .and. 1==i*j*k+offset_d(3)) write(*,*) 'CUDA time_step_CG]'

      li = threadIdx%x
      lj = threadIdx%y
      lk = threadIdx%z

      loc_rhoR(li,lj,lk) = rhoR_d(i,j,k)

      ! Halo Faces
      if (li==1) then
	loc_rhoR(li-1,lj,lk) = rhoR_d(i-1,j,k)
      endif
      if (li==TILE_DIMx) then
	loc_rhoR(li+1,lj,lk) = rhoR_d(i+1,j,k)
      endif

      if (lj==1) then
	loc_rhoR(li,lj-1,lk) = rhoR_d(i,j-1,k)
      endif
      if (lj==TILE_DIMy) then
	loc_rhoR(li,lj+1,lk) = rhoR_d(i,j+1,k)
      endif

      if (lk==1) then
	loc_rhoR(li,lj,lk-1) = rhoR_d(i,j,k-1)
      endif
      if (lk==TILE_DIMz) then
	loc_rhoR(li,lj,lk+1) = rhoR_d(i,j,k+1)
      endif

      ! Halo edges
      if (li==1 .and. lj==1) loc_rhoR(li-1,lj-1,lk) = rhoR_d(i-1,j-1,k)
      if (li==1 .and. lj==TILE_DIMy) loc_rhoR(li-1,lj+1,lk) = rhoR_d(i-1,j+1,k)

      if (li==1 .and. lk==1) loc_rhoR(li-1,lj,lk-1) = rhoR_d(i-1,j,k-1)
      if (li==1 .and. lk==TILE_DIMz) loc_rhoR(li-1,lj,lk+1) = rhoR_d(i-1,j,k+1)

      if (li==TILE_DIMx .and. lj==1) loc_rhoR(li+1,lj-1,lk) = rhoR_d(i+1,j-1,k)
      if (li==TILE_DIMx .and. lj==TILE_DIMy) loc_rhoR(li+1,lj+1,lk) = rhoR_d(i+1,j+1,k)

      if (li==TILE_DIMx .and. lk==1) loc_rhoR(li+1,lj,lk-1) = rhoR_d(i+1,j,k-1)
      if (li==TILE_DIMx .and. lk==TILE_DIMz) loc_rhoR(li+1,lj,lk+1) = rhoR_d(i+1,j,k+1)

      if (lj==1 .and. lk==1) loc_rhoR(li,lj-1,lk-1) = rhoR_d(i,j-1,k-1)
      if (lj==1 .and. lk==TILE_DIMz) loc_rhoR(li,lj-1,lk+1) = rhoR_d(i,j-1,k+1)

      if (lj==TILE_DIMy .and. lk==1) loc_rhoR(li,lj+1,lk-1) = rhoR_d(i,j+1,k-1)
      if (lj==TILE_DIMy .and. lk==TILE_DIMz) loc_rhoR(li,lj+1,lk+1) = rhoR_d(i,j+1,k+1)
             
      ! Halo corner
      if (li==1 .and. lj==1 .and. lk==1) loc_rhoR(li-1,lj-1,lk-1) = rhoR_d(i-1,j-1,k-1) 
      if (li==TILE_DIMx .and. lj==1 .and. lk==1) loc_rhoR(li+1,lj-1,lk-1) = rhoR_d(i+1,j-1,k-1)
      if (li==1 .and. lj==TILE_DIMy .and. lk==1) loc_rhoR(li-1,lj+1,lk-1) = rhoR_d(i-1,j+1,k-1) 
      if (li==1 .and. lj==1 .and. lk==TILE_DIMz) loc_rhoR(li-1,lj-1,lk+1) = rhoR_d(i-1,j-1,k+1) 
      if (li==1 .and. lj==TILE_DIMy .and. lk==TILE_DIMz) loc_rhoR(li-1,lj+1,lk+1) = rhoR_d(i-1,j+1,k+1) 
      if (li==TILE_DIMx .and. lj==1 .and. lk==TILE_DIMz) loc_rhoR(li+1,lj-1,lk+1) = rhoR_d(i+1,j-1,k+1) 
      if (li==TILE_DIMx .and. lj==TILE_DIMy .and. lk==1) loc_rhoR(li+1,lj+1,lk-1) = rhoR_d(i+1,j+1,k-1) 
      if (li==TILE_DIMx .and. lj==TILE_DIMy .and. lk==TILE_DIMz) loc_rhoR(li+1,lj+1,lk+1) = rhoR_d(i+1,j+1,k+1) 
           
      ! Ripeto per il blue
      loc_rhoB(li,lj,lk) = rhoB_d(i,j,k)

      ! Halo Faces
      if (li==1) then
	loc_rhoB(li-1,lj,lk) = rhoB_d(i-1,j,k)
      endif
      if (li==TILE_DIMx) then
	loc_rhoB(li+1,lj,lk) = rhoB_d(i+1,j,k)
      endif

      if (lj==1) then
	loc_rhoB(li,lj-1,lk) = rhoB_d(i,j-1,k)
      endif
      if (lj==TILE_DIMy) then
	loc_rhoB(li,lj+1,lk) = rhoB_d(i,j+1,k)
      endif

      if (lk==1) then
	loc_rhoB(li,lj,lk-1) = rhoB_d(i,j,k-1)
      endif
      if (lk==TILE_DIMz) then
	loc_rhoB(li,lj,lk+1) = rhoB_d(i,j,k+1)
      endif

      ! Halo edges
      if (li==1 .and. lj==1) loc_rhoB(li-1,lj-1,lk) = rhoB_d(i-1,j-1,k)
      if (li==1 .and. lj==TILE_DIMy) loc_rhoB(li-1,lj+1,lk) = rhoB_d(i-1,j+1,k)

      if (li==1 .and. lk==1) loc_rhoB(li-1,lj,lk-1) = rhoB_d(i-1,j,k-1)
      if (li==1 .and. lk==TILE_DIMz) loc_rhoB(li-1,lj,lk+1) = rhoB_d(i-1,j,k+1)

      if (li==TILE_DIMx .and. lj==1) loc_rhoB(li+1,lj-1,lk) = rhoB_d(i+1,j-1,k)
      if (li==TILE_DIMx .and. lj==TILE_DIMy) loc_rhoB(li+1,lj+1,lk) = rhoB_d(i+1,j+1,k)

      if (li==TILE_DIMx .and. lk==1) loc_rhoB(li+1,lj,lk-1) = rhoB_d(i+1,j,k-1)
      if (li==TILE_DIMx .and. lk==TILE_DIMz) loc_rhoB(li+1,lj,lk+1) = rhoB_d(i+1,j,k+1)

      if (lj==1 .and. lk==1) loc_rhoB(li,lj-1,lk-1) = rhoB_d(i,j-1,k-1)
      if (lj==1 .and. lk==TILE_DIMz) loc_rhoB(li,lj-1,lk+1) = rhoB_d(i,j-1,k+1)

      if (lj==TILE_DIMy .and. lk==1) loc_rhoB(li,lj+1,lk-1) = rhoB_d(i,j+1,k-1)
      if (lj==TILE_DIMy .and. lk==TILE_DIMz) loc_rhoB(li,lj+1,lk+1) = rhoB_d(i,j+1,k+1)
         
      ! Halo corner
      if (li==1 .and. lj==1 .and. lk==1) loc_rhoB(li-1,lj-1,lk-1) = rhoB_d(i-1,j-1,k-1) 
      if (li==TILE_DIMx .and. lj==1 .and. lk==1) loc_rhoB(li+1,lj-1,lk-1) = rhoB_d(i+1,j-1,k-1)
      if (li==1 .and. lj==TILE_DIMy .and. lk==1) loc_rhoB(li-1,lj+1,lk-1) = rhoB_d(i-1,j+1,k-1) 
      if (li==1 .and. lj==1 .and. lk==TILE_DIMz) loc_rhoB(li-1,lj-1,lk+1) = rhoB_d(i-1,j-1,k+1) 
      if (li==1 .and. lj==TILE_DIMy .and. lk==TILE_DIMz) loc_rhoB(li-1,lj+1,lk+1) = rhoB_d(i-1,j+1,k+1) 
      if (li==TILE_DIMx .and. lj==1 .and. lk==TILE_DIMz) loc_rhoB(li+1,lj-1,lk+1) = rhoB_d(i+1,j-1,k+1) 
      if (li==TILE_DIMx .and. lj==TILE_DIMy .and. lk==1) loc_rhoB(li+1,lj+1,lk-1) = rhoB_d(i+1,j+1,k-1) 
      if (li==TILE_DIMx .and. lj==TILE_DIMy .and. lk==TILE_DIMz) loc_rhoB(li+1,lj+1,lk+1) = rhoB_d(i+1,j+1,k+1) 
      
      
      do l = 0, npops-1
        loc_tempR(l) = popsR_d(i,j,k, l,flip)
        loc_tempB(l) = popsB_d(i,j,k, l,flip)
      enddo
      call syncthreads

      if (myfluid_d(i,j,k, flip) == fluid_fluid) then
      
        rhoR = loc_rhoR(li,lj,lk)
        rhoB = loc_rhoB(li,lj,lk)
        rhosum = rhoR + rhoB

        if(lcompute_totrho)then
          if(mod(step,niterVTK_d)==0)then
            istatR=atomicadd(totrhobuff_d(1),rhoR_d(i,j,k)+rhoB_d(i,j,k))
          endif
        endif
#ifdef DIFFDENS
        phase_field = (rhoR/densR_d - rhoB/densB_d) / &
           (rhoR/densR_d + rhoB/densB_d)
        viscavg = 1.0 / (0.5*(1.0 + phase_field)*(1.0/viscR_d) & 
			            +0.5*(1.0 - phase_field)*(1.0/viscB_d))
#else
        viscavg = rhoR/rhosum*viscR_d + rhoB/rhosum*viscB_d
#endif
        omega  = 1.0 / ( viscavg/cssq  + 0.5)
        oneminusomega = 1.0 - omega

        invrho = 1.0 / rhosum
        u   = manual_u(loc_tempR) + manual_u(loc_tempB)
        v   = manual_v(loc_tempR) + manual_v(loc_tempB)
        w   = manual_w(loc_tempR) + manual_w(loc_tempB)
        
        vel_d(1,i,j,k) = u
        vel_d(2,i,j,k) = v
        vel_d(3,i,j,k) = w

        FRx= f_cost_d(1)
        FRy= f_cost_d(2)
        FRz= f_cost_d(3)
        FBx= f_cost_d(1)
        FBy= f_cost_d(2)
        FBz= f_cost_d(3)

        ! Psi calc
        psix = ZERO
        psiy = ZERO
        psiz = ZERO
#ifdef DIFFDENS
        psi_R = ZERO
        psi_B = ZERO
        myvel(1) = u
        myvel(2) = v
        myvel(3) = w
#endif
#ifdef HIGHGRAD
        do l = 1, linksd3q27-1
          i1 = li + exd3q27_d(l)
          j1 = lj + eyd3q27_d(l)
          k1 = lk + ezd3q27_d(l)

          rhoR_shifted = loc_rhoR(i1,j1,k1)
          rhoB_shifted = loc_rhoB(i1,j1,k1)
          psi_shifted = (rhoR_shifted/densR_d - rhoB_shifted/densB_d) / &
           (rhoR_shifted/densR_d + rhoB_shifted/densB_d)
          psix = psix + ad3q27_d(l)*exd3q27_d(l)* psi_shifted
          psiy = psiy + ad3q27_d(l)*eyd3q27_d(l)* psi_shifted
          psiz = psiz + ad3q27_d(l)*ezd3q27_d(l)* psi_shifted
#ifdef DIFFDENS
          psi_R(1) = psi_R(1) + ad3q27_d(l)*exd3q27_d(l)* rhoR_shifted
          psi_R(2) = psi_R(2) + ad3q27_d(l)*eyd3q27_d(l)* rhoR_shifted
          psi_R(3) = psi_R(3) + ad3q27_d(l)*ezd3q27_d(l)* rhoR_shifted
          psi_B(1) = psi_B(1) + ad3q27_d(l)*exd3q27_d(l)* rhoB_shifted
          psi_B(2) = psi_B(2) + ad3q27_d(l)*eyd3q27_d(l)* rhoB_shifted
          psi_B(3) = psi_B(3) + ad3q27_d(l)*ezd3q27_d(l)* rhoB_shifted
#endif
        enddo
#else
        do l = 1, npops-1
          i1 = li + ex_d(l)
          j1 = lj + ey_d(l)
          k1 = lk + ez_d(l)

          rhoR_shifted = loc_rhoR(i1,j1,k1)
          rhoB_shifted = loc_rhoB(i1,j1,k1)
          
#ifdef DIFFDENS
          psi_shifted = (rhoR_shifted/densR_d - rhoB_shifted/densB_d) / &
           (rhoR_shifted/densR_d + rhoB_shifted/densB_d)
          psi_R(1) = psi_R(1) + a_d(l)*ex_d(l)* rhoR_shifted
          psi_R(2) = psi_R(2) + a_d(l)*ey_d(l)* rhoR_shifted
          psi_R(3) = psi_R(3) + a_d(l)*ez_d(l)* rhoR_shifted
          psi_B(1) = psi_B(1) + a_d(l)*ex_d(l)* rhoB_shifted
          psi_B(2) = psi_B(2) + a_d(l)*ey_d(l)* rhoB_shifted
          psi_B(3) = psi_B(3) + a_d(l)*ez_d(l)* rhoB_shifted
#else
          psi_shifted = (rhoR_shifted - rhoB_shifted) / (rhoR_shifted + rhoB_shifted)
#endif
          psix = psix + a_d(l)*ex_d(l)* psi_shifted
          psiy = psiy + a_d(l)*ey_d(l)* psi_shifted
          psiz = psiz + a_d(l)*ez_d(l)* psi_shifted
        enddo
#endif
        psinorm_sq = psix*psix + psiy*psiy + psiz*psiz
        
#ifdef NEARCONTACT
        !near contact step
        psinorm = sqrt(psinorm_sq)
        if(nearsel_d(i,j,k)==fluid_fluid)then
        
          if(psinorm .ne. 0.0)then
            dnx=-psix/psinorm
            dny=-psiy/psinorm
            dnz=-psiz/psinorm
          else
            dnx=0.0
            dny=0.0
            dnz=0.0
          endif
          inx=nint(dnx)
          iny=nint(dny)
          inz=nint(dnz)
          do ss=3,6 
            llx=i+ss*inx 
            lly=j+ss*iny
            llz=k+ss*inz
            if(nearsel_d(llx,lly,llz)==fluid_fluid)then
              dist_rep=sqrt((real(ss*inx))**2.0 + (real(ss*iny))**2.0 + &
               (real(ss*inz))**2.0)
              if(dist_rep .ne. 0.0)then
                FRx = FRx - A_rep_d*dnx/dist_rep  !**2 
                FRy = FRy - A_rep_d*dny/dist_rep  !**2 
                FRz = FRz - A_rep_d*dnz/dist_rep  !**2 
              endif
            endif
          enddo
        endif
#endif

        !relaxation step
        ushR = u+FRx/rhosum
        vshR = v+FRy/rhosum
        wshR = w+FRz/rhosum
        ushB = u+FBx/rhosum
        vshB = v+FBy/rhosum
        wshB = w+FBz/rhosum
#ifdef DIFFDENS
        v_psi_R_dot = dot_product(myvel,psi_R)
        v_psi_B_dot = dot_product(myvel,psi_B)
        mat_R=tensor_product(myvel,psi_R)
        mat_B=tensor_product(myvel,psi_B)
        mat_TR=transpose(mat_R)
        mat_TB=transpose(mat_B)
        
        myG_R= mat_R + mat_TR
        myG_B= mat_B + mat_TB
        
        !bgk step
        do l = 0, npops-1
            tens_contr = myG_R(1,1)*cmat_d(1,1,l)+ &
             myG_R(2,2)*cmat_d(2,2,l)+myG_R(3,3)*cmat_d(3,3,l)
            loc_tempR(l) = loc_tempR(l)*oneminusomega + &
             (viscavg*(psi_d(l)*v_psi_R_dot + xi_d(l)*tens_contr) )*omega - &
             equilCG(rhoR,u,v,w,1,l)*oneminusomega + equilCG(rhoR,ushR,vshR,wshR,1,l)
        enddo
        do l = 0, npops-1
            tens_contr = myG_B(1,1)*cmat_d(1,1,l)+ &
             myG_B(2,2)*cmat_d(2,2,l)+myG_B(3,3)*cmat_d(3,3,l)
            loc_tempB(l) = loc_tempB(l)*oneminusomega + &
             (viscavg*(psi_d(l)*v_psi_B_dot + xi_d(l)*tens_contr) )*omega - &
             equilCG(rhoB,u,v,w,2,l)*oneminusomega + equilCG(rhoB,ushB,vshB,wshB,2,l)
        enddo
#else
        !bgk step
        do l = 0, npops-1
            loc_tempR(l) = loc_tempR(l)*oneminusomega - &
             equil(rhoR,u,v,w,l)*oneminusomega + equil(rhoR,ushR,vshR,wshR,l)
        enddo
        do l = 0, npops-1
            loc_tempB(l) = loc_tempB(l)*oneminusomega - &
             equil(rhoB,u,v,w,l)*oneminusomega + equil(rhoB,ushB,vshB,wshB,l)
        enddo
#endif        
        
        !perturbation step
        
        if (psinorm_sq>mylimit) then
#ifndef NEARCONTACT
          psinorm = sqrt(psinorm_sq)
#endif
          acoeff =  NINE/FOUR * omega * sigma_cg_d

          do l = 0, npops-1
              e_dot_psi = ex_d(l)*psix + ey_d(l)*psiy + ez_d(l)*psiz
              temp = psinorm*(p_d(l)*(e_dot_psi*e_dot_psi)/psinorm_sq - b_l_d(l))
              !if(isnan(temp)) temp=ZERO

              loc_fdum(l) = loc_tempR(l) + loc_tempB(l) + acoeff*temp
          enddo

          !recoloring step
          do l=0, npops-1
#ifdef DIFFDENS
              feq = equilCG(rhoR,ZERO,ZERO,ZERO,1,l) + equilCG(rhoB,ZERO,ZERO,ZERO,2,l)
#else
              feq = equil(rhoR,ZERO,ZERO,ZERO,l) + equil(rhoB,ZERO,ZERO,ZERO,l)
#endif
              e_dot_psi = ex_d(l)*psix + ey_d(l)*psiy + ez_d(l)*psiz
              temp1 = rec_fact_d(l) * psinorm

              if (temp1<=mylimit) then
                cosphi = ZERO
              else
                cosphi = e_dot_psi/temp1
              endif

              temp = beta_CG_d * rhoR * rhoB * cosphi/(rhosum*rhosum)

              popsR_d(i,j,k,l, flip) = loc_fdum(l)*rhoR/rhosum + temp*feq
              popsB_d(i,j,k,l, flip) = loc_fdum(l)*rhoB/rhosum - temp*feq
          enddo
        else
        ! Store to mem if not at the interface
          do l = 0, npops-1
            popsR_d(i,j,k,l, flip) = loc_tempR(l)
          enddo
          do l = 0, npops-1
            popsB_d(i,j,k,l, flip) = loc_tempB(l)
          enddo
        endif
      endif
    end subroutine time_step_CG


    attributes(global) subroutine flipflopPop0_CG(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, flop

      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      if (step==1 .and. 1==i*j*k+offset_d(3)) write(*,*) 'CUDA flipflopPop0_CG]'

      flop = 3 - flip
      
      popsR_d(i,j,k,0, flop) = popsR_d(i,j,k,0, flip)
      popsB_d(i,j,k,0, flop) = popsB_d(i,j,k,0, flip)
    end subroutine flipflopPop0_CG

    attributes(global) subroutine flipflopRPop0_CG(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, flop

      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      if (step==1 .and. 1==i*j*k+offset_d(3)) write(*,*) 'CUDA flipflopRPop0_CG]'

      flop = 3 - flip
      
      popsR_d(i,j,k,0, flop) = popsR_d(i,j,k,0, flip)
    end subroutine flipflopRPop0_CG
    
    attributes(global) subroutine flipflopBPop0_CG(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, flop

      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      if (step==1 .and. 1==i*j*k+offset_d(3)) write(*,*) 'CUDA flipflopBPop0_CG]'

      flop = 3 - flip
      
      popsB_d(i,j,k,0, flop) = popsB_d(i,j,k,0, flip)
    end subroutine flipflopBPop0_CG    


    attributes(global) subroutine stream_CG(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l, flop
      integer :: i1,j1,k1
      integer :: atm_i,atm_j,atm_k, atm_st,atm_en, i_atm,i_list


      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      if (step==1 .and. 1==i*j*k+offset_d(3)) write(*,*) 'CUDA stream_CG]'

      flop = 3 - flip


      ! Streaming      
      if (myfluid_d(i,j,k, flip) < fluid_dead) then
        do l = 1, npops-1
          i1 = i + ex_d(l)
          j1 = j + ey_d(l)
          k1 = k + ez_d(l)
          popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
        enddo

        do l = 1, npops-1
          i1 = i + ex_d(l)
          j1 = j + ey_d(l)
          k1 = k + ez_d(l)
          popsB_d(i1,j1,k1,l, flop) = popsB_d(i,j,k,l, flip)
        enddo
      endif
    
    end subroutine stream_CG
    
    attributes(global) subroutine streamR_CG(step, flip)
      integer, value :: step,flip
      integer :: i,j,k,flop

      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      if (step==1 .and. 1==i*j*k+offset_d(3)) write(*,*) 'CUDA streamR_CG]'

      flop = 3 - flip
      if(lcompute_totrho)then
        if(i==1 .and. j==1 .and. k==1)totrhobuff_d(1)=0.0
      endif

      ! Streaming
      if (myfluid_d(i,j,k, flip) < fluid_dead) then
          popsR_d(i  ,j  ,k  , 0, flop) = popsR_d(i,j,k, 0,flip)
          popsR_d(i+1,j  ,k  , 1, flop) = popsR_d(i,j,k, 1,flip)
          popsR_d(i-1,j  ,k  , 2, flop) = popsR_d(i,j,k, 2,flip)
          popsR_d(i  ,j+1,k  , 3, flop) = popsR_d(i,j,k, 3,flip)
          popsR_d(i  ,j-1,k  , 4, flop) = popsR_d(i,j,k, 4,flip)
          popsR_d(i  ,j  ,k+1, 5, flop) = popsR_d(i,j,k, 5,flip)
          popsR_d(i  ,j  ,k-1, 6, flop) = popsR_d(i,j,k, 6,flip)
          popsR_d(i+1,j+1,k  , 7, flop) = popsR_d(i,j,k, 7,flip)
          popsR_d(i-1,j-1,k  , 8, flop) = popsR_d(i,j,k, 8,flip)
          popsR_d(i-1,j+1,k  , 9, flop) = popsR_d(i,j,k, 9,flip)
          popsR_d(i+1,j-1,k  ,10, flop) = popsR_d(i,j,k,10,flip)
          popsR_d(i+1,j  ,k+1,11, flop) = popsR_d(i,j,k,11,flip)
          popsR_d(i-1,j  ,k-1,12, flop) = popsR_d(i,j,k,12,flip)
          popsR_d(i-1,j  ,k+1,13, flop) = popsR_d(i,j,k,13,flip)
          popsR_d(i+1,j  ,k-1,14, flop) = popsR_d(i,j,k,14,flip)
          popsR_d(i  ,j+1,k+1,15, flop) = popsR_d(i,j,k,15,flip)
          popsR_d(i  ,j-1,k-1,16, flop) = popsR_d(i,j,k,16,flip)
          popsR_d(i  ,j-1,k+1,17, flop) = popsR_d(i,j,k,17,flip)
          popsR_d(i  ,j+1,k-1,18, flop) = popsR_d(i,j,k,18,flip)
#ifdef D3Q27
          popsR_d(i+1,j+1,k+1,19, flop) = popsR_d(i,j,k,19,flip)
          popsR_d(i-1,j-1,k-1,20, flop) = popsR_d(i,j,k,20,flip)
          popsR_d(i+1,j+1,k-1,21, flop) = popsR_d(i,j,k,21,flip)
          popsR_d(i-1,j-1,k+1,22, flop) = popsR_d(i,j,k,22,flip)
          popsR_d(i-1,j+1,k+1,23, flop) = popsR_d(i,j,k,23,flip)
          popsR_d(i+1,j-1,k-1,24, flop) = popsR_d(i,j,k,24,flip)
          popsR_d(i+1,j-1,k+1,25, flop) = popsR_d(i,j,k,25,flip)
          popsR_d(i-1,j+1,k-1,26, flop) = popsR_d(i,j,k,26,flip)
#endif
      endif
    end subroutine streamR_CG

    attributes(global) subroutine streamB_CG(step, flip)
      integer, value :: step,flip
      integer :: i,j,k,flop

      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      if (step==1 .and. 1==i*j*k+offset_d(3)) write(*,*) 'CUDA streamB_CG]'

      flop = 3 - flip

      ! Streaming
      if (myfluid_d(i,j,k, flip) < fluid_dead) then
          popsB_d(i  ,j  ,k  , 0, flop) = popsB_d(i,j,k, 0,flip)
          popsB_d(i+1,j  ,k  , 1, flop) = popsB_d(i,j,k, 1,flip)
          popsB_d(i-1,j  ,k  , 2, flop) = popsB_d(i,j,k, 2,flip)
          popsB_d(i  ,j+1,k  , 3, flop) = popsB_d(i,j,k, 3,flip)
          popsB_d(i  ,j-1,k  , 4, flop) = popsB_d(i,j,k, 4,flip)
          popsB_d(i  ,j  ,k+1, 5, flop) = popsB_d(i,j,k, 5,flip)
          popsB_d(i  ,j  ,k-1, 6, flop) = popsB_d(i,j,k, 6,flip)
          popsB_d(i+1,j+1,k  , 7, flop) = popsB_d(i,j,k, 7,flip)
          popsB_d(i-1,j-1,k  , 8, flop) = popsB_d(i,j,k, 8,flip)
          popsB_d(i-1,j+1,k  , 9, flop) = popsB_d(i,j,k, 9,flip)
          popsB_d(i+1,j-1,k  ,10, flop) = popsB_d(i,j,k,10,flip)
          popsB_d(i+1,j  ,k+1,11, flop) = popsB_d(i,j,k,11,flip)
          popsB_d(i-1,j  ,k-1,12, flop) = popsB_d(i,j,k,12,flip)
          popsB_d(i-1,j  ,k+1,13, flop) = popsB_d(i,j,k,13,flip)
          popsB_d(i+1,j  ,k-1,14, flop) = popsB_d(i,j,k,14,flip)
          popsB_d(i  ,j+1,k+1,15, flop) = popsB_d(i,j,k,15,flip)
          popsB_d(i  ,j-1,k-1,16, flop) = popsB_d(i,j,k,16,flip)
          popsB_d(i  ,j-1,k+1,17, flop) = popsB_d(i,j,k,17,flip)
          popsB_d(i  ,j+1,k-1,18, flop) = popsB_d(i,j,k,18,flip)
#ifdef D3Q27
          popsB_d(i+1,j+1,k+1,19, flop) = popsB_d(i,j,k,19,flip)
          popsB_d(i-1,j-1,k-1,20, flop) = popsB_d(i,j,k,20,flip)
          popsB_d(i+1,j+1,k-1,21, flop) = popsB_d(i,j,k,21,flip)
          popsB_d(i-1,j-1,k+1,22, flop) = popsB_d(i,j,k,22,flip)
          popsB_d(i-1,j+1,k+1,23, flop) = popsB_d(i,j,k,23,flip)
          popsB_d(i+1,j-1,k-1,24, flop) = popsB_d(i,j,k,24,flip)
          popsB_d(i+1,j-1,k+1,25, flop) = popsB_d(i,j,k,25,flip)
          popsB_d(i-1,j+1,k-1,26, flop) = popsB_d(i,j,k,26,flip)
#endif
      endif
    end subroutine streamB_CG

    attributes(global) subroutine stream_CG_x(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l, flop
      integer :: i1,j1,k1
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1

      if (j>ny+1) return
      if (k>nz_d+1) return

      flop = 3 - flip

      ! Stream x=0
      i = 0
      if (myfluid_d(i,j,k, flip) < fluid_dead) then
        do l = 1, npops-1
          i1 = i + ex_d(l)
          j1 = j + ey_d(l)
          k1 = k + ez_d(l)
          if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz_d) then
            popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
            popsB_d(i1,j1,k1,l, flop) = popsB_d(i,j,k,l, flip)
          endif
        end do
      endif

      ! Stream x=nx+1
      i = nx+1
      if (myfluid_d(i,j,k, flip) < fluid_dead) then
        do l = 1, npops-1
          i1 = i + ex_d(l)
          j1 = j + ey_d(l)
          k1 = k + ez_d(l)
          if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz_d) then
            popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
            popsB_d(i1,j1,k1,l, flop) = popsB_d(i,j,k,l, flip)
          endif
        end do
      endif
    end subroutine stream_CG_x

    attributes(global) subroutine stream_CG_y(step, flip)
    integer, value :: step,flip
    integer :: i,j,k, l, flop
    integer :: i1,j1,k1

    i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
    k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
    
    if (i>nx+1) return
    if (k>nz_d+1) return

    flop = 3 - flip

    ! Stream y=0
    j = 0
    if (myfluid_d(i,j,k, flip) < fluid_dead) then
      do l = 1, npops-1
        i1 = i + ex_d(l)
        j1 = j + ey_d(l)
        k1 = k + ez_d(l)
        if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz_d) then
          popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
          popsB_d(i1,j1,k1,l, flop) = popsB_d(i,j,k,l, flip)
        endif
      end do
    endif

    ! Stream y=ny+1
    j = ny+1
    if (myfluid_d(i,j,k, flip) < fluid_dead) then
      do l = 1, npops-1
        i1 = i + ex_d(l)
        j1 = j + ey_d(l)
        k1 = k + ez_d(l)
        if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz_d) then
          popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
          popsB_d(i1,j1,k1,l, flop) = popsB_d(i,j,k,l, flip)
        endif
      end do
    endif
  end subroutine stream_CG_y

  attributes(global) subroutine stream_CG_z(step, flip)
    integer, value :: step,flip
    integer :: i,j,k, l, flop
    integer :: i1,j1,k1

    i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
    j = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
    
    if (i>nx+1) return
    if (j>ny+1) return

    flop = 3 - flip

    ! Stream z=0
    k = 0
    if (myfluid_d(i,j,k, flip) < fluid_dead) then
      do l = 1, npops-1
        i1 = i + ex_d(l)
        j1 = j + ey_d(l)
        k1 = k + ez_d(l)
        if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz_d) then
          popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
          popsB_d(i1,j1,k1,l, flop) = popsB_d(i,j,k,l, flip)
        endif
      end do
    endif

    ! Stream z=nz_d+1
    k = nz_d+1
    if (myfluid_d(i,j,k, flip) < fluid_dead) then
      do l = 1, npops-1
        i1 = i + ex_d(l)
        j1 = j + ey_d(l)
        k1 = k + ez_d(l)
        if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz_d) then
          popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
          popsB_d(i1,j1,k1,l, flop) = popsB_d(i,j,k,l, flip)
        endif
      end do
    endif
  end subroutine stream_CG_z


    ! Periodic BC 	-----------------------------------------------
    attributes(global) subroutine bc_per_x(step, flip)
      integer, value :: step,flip
      integer :: j,k !, l
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
      ! if (1==step .and. 1==j*k) write(*,*) 'CUDA bc_per_x] i,j:',j,k

      ! Apply BC also on virtual nodes...
      if (j<=ny+1 .and. k<=nz_d+1) then
        popsR_d(1, j,k, 1, flip) = popsR_d(nx+1,j,k, 1, flip)
        popsR_d(1, j,k, 7, flip) = popsR_d(nx+1,j,k, 7, flip)
        popsR_d(1, j,k,10, flip) = popsR_d(nx+1,j,k,10, flip)
        popsR_d(1, j,k,11, flip) = popsR_d(nx+1,j,k,11, flip)
        popsR_d(1, j,k,14, flip) = popsR_d(nx+1,j,k,14, flip)

        popsR_d(nx,j,k, 2, flip) = popsR_d(0,j,k, 2, flip)
        popsR_d(nx,j,k, 8, flip) = popsR_d(0,j,k, 8, flip)
        popsR_d(nx,j,k, 9, flip) = popsR_d(0,j,k, 9, flip)
        popsR_d(nx,j,k,12, flip) = popsR_d(0,j,k,12, flip)
        popsR_d(nx,j,k,13, flip) = popsR_d(0,j,k,13, flip)

        ! do l = 1, npops-1
        !      if (ex_d(l)>0) then
        !       if (1==step .and. 1==j*k) write(*,*) 'CUDA bc_per_x] 1<-nx+1 l',l
        !       popsR_d(1, j,k, l, flip) = popsR_d(nx+1,j,k, l, flip)
        !      endif
        !      if (ex_d(l)<0) then
        !       if (1==step .and. 1==j*k) write(*,*) 'CUDA bc_per_x] nx<-0 l',l
        !       popsR_d(nx,j,k, l, flip) = popsR_d(0,j,k, l, flip)
        !      endif
        ! end do
      endif
    end subroutine bc_per_x


    attributes(global) subroutine bc_per_y(step, flip)
      integer, value :: step,flip
      integer :: i,k, l
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
      ! if (1==step .and. 1==i*k) write(*,*) 'CUDA bc_per_y] i,k:',i,k

      if (0<=i .and. i<=nx+1 .and. 0<=k .and. k<=nz_d+1) then
        do l = 1, npops-1
             if (ey_d(l)>0) popsR_d(i,1,k, l, flip) = popsR_d(i,ny+1,k, l, flip)
             if (ey_d(l)<0) popsR_d(i,ny,k, l, flip) = popsR_d(i,0,k, l, flip)
        end do
      endif
    end subroutine bc_per_y


    attributes(global) subroutine bc_per_z(step, flip)
      integer, value :: step,flip
      integer :: i,j, l
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      j = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
      ! if (1==step .and. 1==i*j) write(*,*) 'CUDA bc_per_z] i,j:',i,j

      if (0<=j .and. j<=ny+1 .and. 0<=i .and. i<=nx+1) then
        do l = 1, npops-1
             if (ez_d(l)>0) popsR_d(i,j,1, l, flip) = popsR_d(i,j,nz_d+1, l, flip)
             if (ez_d(l)<0) popsR_d(i,j,nz_d, l, flip) = popsR_d(i,j,0, l, flip)
        end do
      endif
    end subroutine bc_per_z

    ! Periodic BC 	2 fluids ----------------------------------------------
    attributes(global) subroutine bc_per_x2(step, flip, faiPops,faiRho)
      integer, value :: step,flip
      logical, value :: faiPops, faiRho
      integer :: j,k, l
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
      ! if (1==step .and. 1==j*k) write(*,*) 'CUDA bc_per_x2] i,j:',j,k

      if (1<=j .and. j<=ny .and. 1<=k .and. k<=nz_d) then
        if (faiPops) then
          do l = 1, npops-1
            popsR_d(nx+1, j,k, l, flip) = popsR_d(1,j,k, l, flip)
            popsB_d(nx+1, j,k, l, flip) = popsB_d(1,j,k, l, flip)
          
            popsR_d(0,j,k, l, flip) = popsR_d(nx,j,k, l, flip)
            popsB_d(0,j,k, l, flip) = popsB_d(nx,j,k, l, flip)
          end do
        endif

        if (faiRho) then
          rhoR_d (nx+1,j,k) = rhoR_d (1,j,k)
          rhoB_d (nx+1,j,k) = rhoB_d (1,j,k)

          rhoR_d (0,j,k) = rhoR_d (nx,j,k)
          rhoB_d (0,j,k) = rhoB_d (nx,j,k)
        endif
      endif
    end subroutine bc_per_x2


    attributes(global) subroutine bc_per_y2(step, flip, faiPops,faiRho)
      integer, value :: step,flip
      logical, value :: faiPops, faiRho
      integer :: i,k, l
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
      ! if (1==step .and. 1==i*k) write(*,*) 'CUDA bc_per_y2] i,k:',i,k

      if (1<=i .and. i<=nx .and. 1<=k .and. k<=nz_d) then
        if (faiPops) then
          do l = 1, npops-1
            popsR_d(i,ny+1,k, l, flip) = popsR_d(i,1,k, l, flip)
            popsB_d(i,ny+1,k, l, flip) = popsB_d(i,1,k, l, flip)
          
            popsR_d(i,0,k, l, flip) = popsR_d(i,ny,k, l, flip)
            popsB_d(i,0,k, l, flip) = popsB_d(i,ny,k, l, flip)
          end do
        endif

        if (faiRho) then
          rhoR_d(i,ny+1,k) = rhoR_d(i,1,k)
          rhoB_d(i,ny+1,k) = rhoB_d(i,1,k)

          rhoR_d(i,0,k) = rhoR_d(i,ny,k)
          rhoB_d(i,0,k) = rhoB_d(i,ny,k)
        endif
      endif
    end subroutine bc_per_y2


    attributes(global) subroutine bc_per_z2(step, flip, faiPops,faiRho)
      integer, value :: step,flip
      logical, value :: faiPops, faiRho
      integer :: i,j, l
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      j = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
      ! if (1==step .and. 1==i*j) write(*,*) 'CUDA bc_per_z2] i,j:',i,j

      if (1<=j .and. j<=ny .and. 1<=i .and. i<=nx) then
        if (faiPops) then
          do l = 1, npops-1        
            popsR_d(i,j,nz_d+1, l, flip) = popsR_d(i,j,1, l, flip)
            popsB_d(i,j,nz_d+1, l, flip) = popsB_d(i,j,1, l, flip)
          
            popsR_d(i,j,0, l, flip) = popsR_d(i,j,nz_d, l, flip)
            popsB_d(i,j,0, l, flip) = popsB_d(i,j,nz_d, l, flip)
          end do
        endif

        if (faiRho) then
          rhoR_d(i,j,nz_d+1) = rhoR_d(i,j,1)
          rhoB_d(i,j,nz_d+1) = rhoB_d(i,j,1)

          rhoR_d(i,j,0) = rhoR_d(i,j,nz_d)
          rhoB_d(i,j,0) = rhoB_d(i,j,nz_d)
        endif
      endif
    end subroutine bc_per_z2


    attributes(global) subroutine bc_edge_z2(step, flip, faiPops,faiRho)
      integer, value :: step,flip
      logical, value :: faiPops, faiRho
      integer :: k, l
  
      k = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1

      if (0<=k .and. k<=nz_d+1) then
        if (faiPops) then
          do l = 1, npops-1
            popsR_d(0,0,k, l, flip) = popsR_d(nx,ny,k, l, flip)
            popsB_d(0,0,k, l, flip) = popsB_d(nx,ny,k, l, flip)

            popsR_d(0,ny+1,k, l, flip) = popsR_d(nx,1,k, l, flip)
            popsB_d(0,ny+1,k, l, flip) = popsB_d(nx,1,k, l, flip)

            popsR_d(nx+1,0,k, l, flip) = popsR_d(1,ny,k, l, flip)
            popsB_d(nx+1,0,k, l, flip) = popsB_d(1,ny,k, l, flip)

            popsR_d(nx+1,ny+1,k, l, flip) = popsR_d(1,1,k, l, flip)
            popsB_d(nx+1,ny+1,k, l, flip) = popsB_d(1,1,k, l, flip)
          end do
        endif

        if (faiRho) then
          rhoR_d (0,0,k) = rhoR_d (nx,ny,k)
          rhoB_d (0,0,k) = rhoB_d (nx,ny,k)

          rhoR_d (0,ny+1,k) = rhoR_d (nx,1,k)
          rhoB_d (0,ny+1,k) = rhoB_d (nx,1,k)

          rhoR_d (nx+1,0,k) = rhoR_d (1,ny,k)
          rhoB_d (nx+1,0,k) = rhoB_d (1,ny,k)

          rhoR_d (nx+1,ny+1,k) = rhoR_d (1,1,k)
          rhoB_d (nx+1,ny+1,k) = rhoB_d (1,1,k)
        endif
      endif
    end subroutine bc_edge_z2

    attributes(global) subroutine bc_edge_y2(step, flip, faiPops,faiRho)
      integer, value :: step,flip
      logical, value :: faiPops, faiRho
      integer :: j, l
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1

      if (1<=j .and. j<=ny) then
        if (faiPops) then
          do l = 1, npops-1
            popsR_d(0,j,0, l, flip) = popsR_d(nx,j,nz_d, l, flip)
            popsB_d(0,j,0, l, flip) = popsB_d(nx,j,nz_d, l, flip)

            popsR_d(0,j,nz_d+1, l, flip) = popsR_d(nx,j,1, l, flip)
            popsB_d(0,j,nz_d+1, l, flip) = popsB_d(nx,j,1, l, flip)

            popsR_d(nx+1,j,0, l, flip) = popsR_d(1,j,nz_d, l, flip)
            popsB_d(nx+1,j,0, l, flip) = popsB_d(1,j,nz_d, l, flip)

            popsR_d(nx+1,j,nz_d+1, l, flip) = popsR_d(1,j,1, l, flip)
            popsB_d(nx+1,j,nz_d+1, l, flip) = popsB_d(1,j,1, l, flip)
          end do
        endif

        if (faiRho) then
          rhoR_d (0,j,0) = rhoR_d (nx,j,nz_d)
          rhoB_d (0,j,0) = rhoB_d (nx,j,nz_d)

          rhoR_d (0,j,nz_d+1) = rhoR_d (nx,j,1)
          rhoB_d (0,j,nz_d+1) = rhoB_d (nx,j,1)

          rhoR_d (nx+1,j,0) = rhoR_d (1,j,nz_d)
          rhoB_d (nx+1,j,0) = rhoB_d (1,j,nz_d)

          rhoR_d (nx+1,j,nz_d+1) = rhoR_d (1,j,1)
          rhoB_d (nx+1,j,nz_d+1) = rhoB_d (1,j,1)
        endif
      endif
    end subroutine bc_edge_y2

    attributes(global) subroutine bc_edge_x2(step, flip, faiPops,faiRho)
      integer, value :: step,flip
      logical, value :: faiPops, faiRho
      integer :: i, l
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1

      if (1<=i .and. i<=nx) then
        if (faiPops) then
          do l = 1, npops-1
            popsR_d(i,0,0, l, flip) = popsR_d(i,ny,nz_d, l, flip)
            popsB_d(i,0,0, l, flip) = popsB_d(i,ny,nz_d, l, flip)

            popsR_d(i,0,nz_d+1, l, flip) = popsR_d(i,ny,1, l, flip)
            popsB_d(i,0,nz_d+1, l, flip) = popsB_d(i,ny,1, l, flip)

            popsR_d(i,ny+1,0, l, flip) = popsR_d(i,1,nz_d, l, flip)
            popsB_d(i,ny+1,0, l, flip) = popsB_d(i,1,nz_d, l, flip)

            popsR_d(i,ny+1,nz_d+1, l, flip) = popsR_d(i,1,1, l, flip)
            popsB_d(i,ny+1,nz_d+1, l, flip) = popsB_d(i,1,1, l, flip)
          end do
        endif

        if (faiRho) then
          rhoR_d (i,0,0) = rhoR_d (i,ny,nz_d)
          rhoB_d (i,0,0) = rhoB_d (i,ny,nz_d)

          rhoR_d (i,0,nz_d+1) = rhoR_d (i,ny,1)
          rhoB_d (i,0,nz_d+1) = rhoB_d (i,ny,1)

          rhoR_d (i,ny+1,0) = rhoR_d (i,1,nz_d)
          rhoB_d (i,ny+1,0) = rhoB_d (i,1,nz_d)

          rhoR_d (i,ny+1,nz_d+1) = rhoR_d (i,1,1)
          rhoB_d (i,ny+1,nz_d+1) = rhoB_d (i,1,1)
        endif
      endif
    end subroutine bc_edge_x2


    attributes(global) subroutine bc_corners(step, flip)
      integer, value :: step,flip
      integer :: j,k, l
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
    end subroutine bc_corners
    
    attributes(global) subroutine applybchalo_CG_R(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l, flop, lopp, itemp, mytype
      integer :: i1,j1,k1
      real :: u,v,w,uv,rho

      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x - 1
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y - 1
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z - 1
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA stream_BGK] i,k:',i,j,k

      flop = 3 - flip
      
      if (i>nx+1) return
      if (j>ny+1) return
      if (k>nz_d+1) return
      
      call applybc_R(i,j,k,flip) 
       
    end subroutine applybchalo_CG_R
    
    attributes(global) subroutine applybc_CG_R(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l, flop, lopp, itemp, mytype
      integer :: i1,j1,k1
      real :: u,v,w,uv,rho

      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA stream_BGK] i,k:',i,j,k

      flop = 3 - flip
      
      if (i>nx) return
      if (j>ny) return
      if (k>nz_d) return
        
      call applybc_R(i,j,k,flip)  
       
    end subroutine applybc_CG_R
    
    attributes(global) subroutine applybc_CG_R_x(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l, flop
      integer :: i1,j1,k1
      real :: rho1,rho2,rtemp
      
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1

      if (j>ny+1) return
      if (k>nz_d+1) return

      flop = 3 - flip

      ! do x=0
      i = 0
      call applybc_R(i,j,k,flip)     
      

      ! do x=nx+1
      i = nx+1
      call applybc_R(i,j,k,flip)      
      
    end subroutine applybc_CG_R_x

    attributes(global) subroutine applybc_CG_R_y(step, flip)
    integer, value :: step,flip
    integer :: i,j,k, l, flop
    integer :: i1,j1,k1
    real :: rho1,rho2,rtemp

    i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
    k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
    
    if (i>nx+1) return
    if (k>nz_d+1) return

    flop = 3 - flip

    ! do y=0
    j = 0
    call applybc_R(i,j,k,flip)  

    ! do y=ny+1
    j = ny+1
    call applybc_R(i,j,k,flip)  
    
  end subroutine applybc_CG_R_y

  attributes(global) subroutine applybc_CG_R_z(step, flip)
    integer, value :: step,flip
    integer :: i,j,k, l, flop
    integer :: i1,j1,k1
    real :: rho1,rho2,rtemp

    i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
    j = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
    
    if (i>nx+1) return
    if (j>ny+1) return

    flop = 3 - flip

    ! do z=0
    k = 0  
    call applybc_R(i,j,k,flip)  

    ! do z=nz_d+1
    k = nz_d+1
    call applybc_R(i,j,k,flip)  
    
  end subroutine applybc_CG_R_z
    
    attributes(global) subroutine applybchalo_CG_B(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l, flop, lopp, itemp, mytype
      integer :: i1,j1,k1
      real :: u,v,w,uv,rho

      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x - 1
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y - 1
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z - 1
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA stream_BGK] i,k:',i,j,k

      flop = 3 - flip
      
      if (i>nx+1) return
      if (j>ny+1) return
      if (k>nz_d+1) return
      
      call applybc_B(i,j,k,flip)
       
    end subroutine applybchalo_CG_B
    
     attributes(global) subroutine applybc_CG_B(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l, flop, lopp, itemp, mytype
      integer :: i1,j1,k1
      real :: u,v,w,uv,rho

      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA stream_BGK] i,k:',i,j,k

      flop = 3 - flip
      
      if (i>nx) return
      if (j>ny) return
      if (k>nz_d) return
      
      call applybc_B(i,j,k,flip)
       
    end subroutine applybc_CG_B
    
    attributes(global) subroutine applybc_CG_B_x(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l, flop
      integer :: i1,j1,k1
      real :: rho1,rho2,rtemp
      
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1

      if (j>ny+1) return
      if (k>nz_d+1) return

      flop = 3 - flip

      ! do x=0
      i = 0
      call applybc_B(i,j,k,flip)     
      

      ! do x=nx+1
      i = nx+1
      call applybc_B(i,j,k,flip)      
      
    end subroutine applybc_CG_B_x

    attributes(global) subroutine applybc_CG_B_y(step, flip)
    integer, value :: step,flip
    integer :: i,j,k, l, flop
    integer :: i1,j1,k1
    real :: rho1,rho2,rtemp

    i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
    k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
    
    if (i>nx+1) return
    if (k>nz_d+1) return

    flop = 3 - flip

    ! do y=0
    j = 0
    call applybc_B(i,j,k,flip)  

    ! do y=ny+1
    j = ny+1
    call applybc_B(i,j,k,flip)  
    
  end subroutine applybc_CG_B_y

  attributes(global) subroutine applybc_CG_B_z(step, flip)
    integer, value :: step,flip
    integer :: i,j,k, l, flop
    integer :: i1,j1,k1
    real :: rho1,rho2,rtemp

    i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
    j = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
    
    if (i>nx+1) return
    if (j>ny+1) return

    flop = 3 - flip

    ! do z=0
    k = 0  
    call applybc_B(i,j,k,flip)  

    ! do z=nz_d+1
    k = nz_d+1
    call applybc_B(i,j,k,flip)  
    
  end subroutine applybc_CG_B_z
    
    
    attributes(device) subroutine applybc_R(i,j,k,flip)
    
    implicit none
    
    integer, intent(in) :: i,j,k,flip
    integer :: l,i1,j1,k1,lopp,itemp,mytype
    real :: rho,u,v,w,uv
    
    if (myfluid_d(i,j,k, flip) <= fluid_wall) then
        itemp=-int(myfluid_d(i,j,k, flip))
        if(itemp>0)then
          mytype=bctype_d(1,itemp)
          if(mytype==1)then
            rho = bcrho_d(1,itemp)
            do l = 1, npops-1
              i1 = i + ex_d(l)
              j1 = j + ey_d(l)
              k1 = k + ez_d(l)
              if (myfluid_d(i1,j1,k1, flip) == fluid_fluid) then
                lopp=opp_d(l)
                u = vel_d(1,i1,j1,k1)
                v = vel_d(2,i1,j1,k1)
                w = vel_d(3,i1,j1,k1)
                
                uv = (ONE/cssq) * (u*ex_d(lopp) + v*ey_d(lopp) + w*ez_d(lopp))
#ifdef DIFFDENS 
                popsR_d(i,j,k,l, flip) = -popsR_d(i1,j1,k1,lopp, flip) + &
                 TWO* rho*(phi_d(lopp) + varphi_d(lopp)*alphaCG_d(1)+ &
                 p_d(lopp)*(HALF*(uv*uv)-(HALF/cssq)*(u*u + v*v + w*w)))           
#else               
                popsR_d(i,j,k,l, flip) = -popsR_d(i1,j1,k1,lopp, flip) + &
                 TWO* rho*p_d(lopp)*(ONE+HALF*(uv*uv)-(HALF/cssq)*(u*u + v*v + w*w))
#endif
              endif
            enddo
          elseif(mytype==2)then
            u = bcvel_d(1,itemp)
            v = bcvel_d(2,itemp)
            w = bcvel_d(3,itemp)
            do l = 1, npops-1
              i1 = i + ex_d(l)
              j1 = j + ey_d(l)
              k1 = k + ez_d(l)
              if (myfluid_d(i1,j1,k1, flip) == fluid_fluid) then
                lopp=opp_d(l)
                rho = rhoR_d(i1,j1,k1)
                
                uv = (ONE/cssq) * (u*ex_d(lopp) + v*ey_d(lopp) + w*ez_d(lopp))
               
                popsR_d(i,j,k,l, flip) = popsR_d(i1,j1,k1,lopp, flip) - &
                 TWO* rho*p_d(lopp)*uv
              endif
            enddo
          elseif(mytype==3)then
            rho = bcrho_d(1,itemp)
            u = bcvel_d(1,itemp)
            v = bcvel_d(2,itemp)
            w = bcvel_d(3,itemp)
            do l = 1, npops-1
              i1 = i + ex_d(l)
              j1 = j + ey_d(l)
              k1 = k + ez_d(l)
              if (myfluid_d(i1,j1,k1, flip) == fluid_fluid) then
                lopp=opp_d(l)
                
                uv = (ONE/cssq) * (u*ex_d(lopp) + v*ey_d(lopp) + w*ez_d(lopp))
               
                popsR_d(i,j,k,l, flip) = popsR_d(i1,j1,k1,lopp, flip) - &
                 TWO* rho*p_d(lopp)*uv
              endif
            enddo
          else
            do l = 1, npops-1
              i1 = i + ex_d(l)
              j1 = j + ey_d(l)
              k1 = k + ez_d(l)
              if (myfluid_d(i1,j1,k1, flip) == fluid_fluid) then
                lopp=opp_d(l)
                popsR_d(i,j,k,l, flip) = popsR_d(i1,j1,k1,lopp, flip)
              endif
            enddo
          endif
        else
          do l = 1, npops-1
            i1 = i + ex_d(l)
            j1 = j + ey_d(l)
            k1 = k + ez_d(l)
            if (myfluid_d(i1,j1,k1, flip) == fluid_fluid) then
              lopp=opp_d(l)
              popsR_d(i,j,k,l, flip) = popsR_d(i1,j1,k1,lopp, flip)
            endif
          enddo
        endif
      endif  
      
    end subroutine applybc_R
    
    attributes(device) subroutine applybc_B(i,j,k,flip)
    
    implicit none
    
    integer, intent(in) :: i,j,k,flip
    integer :: l,i1,j1,k1,lopp,itemp,mytype
    real :: rho,u,v,w,uv
    
    if (myfluid_d(i,j,k, flip) <= fluid_wall) then
        itemp=-int(myfluid_d(i,j,k, flip))
        if(itemp>0)then
          mytype=bctype_d(2,itemp)
          if(mytype==1)then
            rho = bcrho_d(2,itemp)
            do l = 1, npops-1
              i1 = i + ex_d(l)
              j1 = j + ey_d(l)
              k1 = k + ez_d(l)
              if (myfluid_d(i1,j1,k1, flip) == fluid_fluid) then
                lopp=opp_d(l)
                u = vel_d(1,i1,j1,k1)
                v = vel_d(2,i1,j1,k1)
                w = vel_d(3,i1,j1,k1)
                
                uv = (ONE/cssq) * (u*ex_d(lopp) + v*ey_d(lopp) + w*ez_d(lopp))
#ifdef DIFFDENS 
                popsB_d(i,j,k,l, flip) = -popsB_d(i1,j1,k1,lopp, flip) + &
                 TWO* rho*(phi_d(lopp) + varphi_d(lopp)*alphaCG_d(2)+ &
                 p_d(lopp)*(HALF*(uv*uv)-(HALF/cssq)*(u*u + v*v + w*w)))           
#else                  
                popsB_d(i,j,k,l, flip) = -popsB_d(i1,j1,k1,lopp, flip) + &
                 TWO* rho*p_d(lopp)*(ONE+HALF*(uv*uv)-(HALF/cssq)*(u*u + v*v + w*w))
#endif
              endif
            enddo
          elseif(mytype==2)then
            u = bcvel_d(1,itemp)
            v = bcvel_d(2,itemp)
            w = bcvel_d(3,itemp)
            do l = 1, npops-1
              i1 = i + ex_d(l)
              j1 = j + ey_d(l)
              k1 = k + ez_d(l)
              if (myfluid_d(i1,j1,k1, flip) == fluid_fluid) then
                lopp=opp_d(l)
                rho = rhoB_d(i1,j1,k1)
                
                uv = (ONE/cssq) * (u*ex_d(lopp) + v*ey_d(lopp) + w*ez_d(lopp))
               
                popsB_d(i,j,k,l, flip) = popsB_d(i1,j1,k1,lopp, flip) - &
                 TWO* rho*p_d(lopp)*uv
              endif
            enddo
          elseif(mytype==3)then
            rho = bcrho_d(2,itemp)
            u = bcvel_d(1,itemp)
            v = bcvel_d(2,itemp)
            w = bcvel_d(3,itemp)
            do l = 1, npops-1
              i1 = i + ex_d(l)
              j1 = j + ey_d(l)
              k1 = k + ez_d(l)
              if (myfluid_d(i1,j1,k1, flip) == fluid_fluid) then
                lopp=opp_d(l)
                
                uv = (ONE/cssq) * (u*ex_d(lopp) + v*ey_d(lopp) + w*ez_d(lopp))
               
                popsB_d(i,j,k,l, flip) = popsB_d(i1,j1,k1,lopp, flip) - &
                 TWO* rho*p_d(lopp)*uv
              endif
            enddo
          else
            do l = 1, npops-1
              i1 = i + ex_d(l)
              j1 = j + ey_d(l)
              k1 = k + ez_d(l)
              if (myfluid_d(i1,j1,k1, flip) == fluid_fluid) then
                lopp=opp_d(l)
                popsB_d(i,j,k,l, flip) = popsB_d(i1,j1,k1,lopp, flip)
              endif
            enddo
          endif
        else
          do l = 1, npops-1
            i1 = i + ex_d(l)
            j1 = j + ey_d(l)
            k1 = k + ez_d(l)
            if (myfluid_d(i1,j1,k1, flip) == fluid_fluid) then
              lopp=opp_d(l)
              popsB_d(i,j,k,l, flip) = popsB_d(i1,j1,k1,lopp, flip)
            endif
          enddo
        endif
      endif  
      
    end subroutine applybc_B

    !!! Periodic treat myfluid
    attributes(global) subroutine isfluid_per_x(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y

      if (j<=ny .and. k<=nz_d) then
        do i=1,nbuff
          myfluid_d(nx+i, j,k, flip) = myfluid_d( i,j,k, flip)
          myfluid_d(1-i, j,k, flip) = myfluid_d(nx+1-i,j,k, flip)
        enddo
      endif
    end subroutine isfluid_per_x


    attributes(global) subroutine isfluid_per_y(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y

      if (i<=nx .and. k<=nz_d) then
        do j=1,nbuff
          myfluid_d(i,ny+j,k, flip) = myfluid_d(i, j,k, flip)
          myfluid_d(i,1-j,k, flip) = myfluid_d(i,ny+1-j,k, flip)
        enddo
      endif
    end subroutine isfluid_per_y


    attributes(global) subroutine isfluid_per_z(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIM + threadIdx%y

      if (j<=ny .and. i<=nx) then
        do k=1,nbuff
          myfluid_d(i,j,nz_d+k, flip) = myfluid_d(i,j, k, flip)        
          myfluid_d(i,j,1-k, flip) = myfluid_d(i,j,nz_d+1-k, flip)
        enddo
      endif
    end subroutine isfluid_per_z


    attributes(global) subroutine isfluid_edge_z(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l
  
      k = (blockIdx%x-1) * TILE_DIM + threadIdx%x - nbuff

      if (k<=nz_d+nbuff) then
        do j=1,nbuff
          do i=1,nbuff
             myfluid_d (1-i,1-j,k, flip) = myfluid_d (nx+1-i,ny+1-j,k, flip)
             myfluid_d (1-i,ny+j,k, flip) = myfluid_d (nx+1-i,j,k, flip)
             myfluid_d (nx+i,1-j,k, flip) = myfluid_d ( i,ny+1-j,k, flip)
             myfluid_d (nx+i,ny+j,k, flip) = myfluid_d ( i, j,k, flip)
           enddo
         enddo
      endif
    end subroutine isfluid_edge_z

    attributes(global) subroutine isfluid_edge_y(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x

      if (j<=ny) then
        do k=1,nbuff
          do i=1,nbuff
            myfluid_d (1-i,j,1-k, flip) = myfluid_d (nx+1-i,j, nz_d+1-k, flip)
            myfluid_d (1-i,j,nz_d+k, flip) = myfluid_d (nx+1-i,j,  k, flip)
            myfluid_d (nx+i,j,1-k, flip) = myfluid_d ( i,j, nz_d+1-k, flip)
            myfluid_d (nx+i,j,nz_d+k, flip) = myfluid_d ( i,j,  k, flip)
          enddo
        enddo
      endif
    end subroutine isfluid_edge_y

    attributes(global) subroutine isfluid_edge_x(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x

      if (i<=nx) then
        do k=1,nbuff
          do j=1,nbuff
            myfluid_d (i,1-j,1-k, flip) = myfluid_d (i,ny+1-j,nz_d+1-k, flip)
            myfluid_d (i,1-j,nz_d+k, flip) = myfluid_d (i,ny+1-j,k, flip)
            myfluid_d (i,ny+j,1-k, flip) = myfluid_d (i,j,nz_d+1-k, flip)
            myfluid_d (i,ny+j,nz_d+k, flip) = myfluid_d (i, j, k, flip)
          enddo
        enddo
      endif
    end subroutine isfluid_edge_x
    
    
    !nearsl
    
        !!! Periodic treat myfluid
    attributes(global) subroutine nearsel_per_x(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y

      if (j<=ny .and. k<=nz_d) then
        do i=1,nbuff
          nearsel_d(nx+i, j,k) = nearsel_d( i,j,k)
          nearsel_d(1-i, j,k) = nearsel_d(nx+1-i,j,k)
        enddo
      endif
    end subroutine nearsel_per_x


    attributes(global) subroutine nearsel_per_y(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y

      if (i<=nx .and. k<=nz_d) then
        do j=1,nbuff
          nearsel_d(i,ny+j,k) = nearsel_d(i, j,k)
          nearsel_d(i,1-j,k) = nearsel_d(i,ny+1-j,k)
        enddo
      endif
    end subroutine nearsel_per_y


    attributes(global) subroutine nearsel_per_z(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIM + threadIdx%y

      if (j<=ny .and. i<=nx) then
        do k=1,nbuff
          nearsel_d(i,j,nz_d+k) = nearsel_d(i,j, k)        
          nearsel_d(i,j,1-k) = nearsel_d(i,j,nz_d+1-k)
        enddo
      endif
    end subroutine nearsel_per_z


    attributes(global) subroutine nearsel_edge_z(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l
  
      k = (blockIdx%x-1) * TILE_DIM + threadIdx%x - nbuff

      if (k<=nz_d+nbuff) then
        do j=1,nbuff
          do i=1,nbuff
             nearsel_d (1-i,1-j,k) = nearsel_d (nx+1-i,ny+1-j,k)
             nearsel_d (1-i,ny+j,k) = nearsel_d (nx+1-i,j,k)
             nearsel_d (nx+i,1-j,k) = nearsel_d ( i,ny+1-j,k)
             nearsel_d (nx+i,ny+j,k) = nearsel_d ( i, j,k)
           enddo
         enddo
      endif
    end subroutine nearsel_edge_z

    attributes(global) subroutine nearsel_edge_y(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x

      if (j<=ny) then
        do k=1,nbuff
          do i=1,nbuff
            nearsel_d (1-i,j,1-k) = nearsel_d (nx+1-i,j, nz_d+1-k)
            nearsel_d (1-i,j,nz_d+k) = nearsel_d (nx+1-i,j,  k)
            nearsel_d (nx+i,j,1-k) = nearsel_d ( i,j, nz_d+1-k)
            nearsel_d (nx+i,j,nz_d+k) = nearsel_d ( i,j,  k)
          enddo
        enddo
      endif
    end subroutine nearsel_edge_y

    attributes(global) subroutine nearsel_edge_x(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x

      if (i<=nx) then
        do k=1,nbuff
          do j=1,nbuff
            nearsel_d (i,1-j,1-k) = nearsel_d (i,ny+1-j,nz_d+1-k)
            nearsel_d (i,1-j,nz_d+k) = nearsel_d (i,ny+1-j,k)
            nearsel_d (i,ny+j,1-k) = nearsel_d (i,j,nz_d+1-k)
            nearsel_d (i,ny+j,nz_d+k) = nearsel_d (i, j, k)
          enddo
        enddo
      endif
    end subroutine nearsel_edge_x


    attributes(global) subroutine bc_periodic_ext_x(step, flip)
      integer, value :: step,flip
      integer :: j,k, l
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
      ! if (1==step .and. 1==j*k) write(*,*) 'CUDA bc_periodic_ext_x] i,j:',j,k

      ! Copy on virtual nodes...
      if (0<=j .and. j<=ny+1 .and. 0<=k .and. k<=nz_d+1) then
        do l = 0, npops-1
          popsR_d(0   ,j,k, l, flip) = popsR_d(nx,j,k, l, flip)
          popsR_d(nx+1,j,k, l, flip) = popsR_d(1 ,j,k, l, flip)
        end do
      endif
    end subroutine bc_periodic_ext_x

    attributes(global) subroutine bc_periodic_ext_y(step, flip)
      integer, value :: step,flip
      integer :: i,k, l
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
      ! if (1==step .and. 1==i*k) write(*,*) 'CUDA bc_periodic_ext_y] i,k:',i,k

      if (0<=i .and. i<=nx+1 .and. 0<=k .and. k<=nz_d+1) then
        do l = 0, npops-1
          popsR_d(i,0   ,k, l, flip) = popsR_d(i,ny,k, l, flip)
          popsR_d(i,ny+1,k, l, flip) = popsR_d(i,1 ,k, l, flip)
        end do
      endif
    end subroutine bc_periodic_ext_y

    attributes(global) subroutine bc_periodic_ext_z(step, flip)
      integer, value :: step,flip
      integer :: i,j, l
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      j = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
      ! if (1==step .and. 1==i*j) write(*,*) 'CUDA bc_periodic_ext_z] i,j:',i,j

      if (0<=j .and. j<=ny+1 .and. 0<=i .and. i<=nx+1) then
        do l = 0, npops-1
          popsR_d(i,j,0   , l, flip) = popsR_d(i,j,nz_d, l, flip)
          popsR_d(i,j,nz_d+1, l, flip) = popsR_d(i,j,1 , l, flip)
        end do
      endif
    end subroutine bc_periodic_ext_z


    ! Bounce-Back BC
    attributes(global) subroutine bc_bb_x(step, flip)
      real :: rho1,rho2
      integer, value :: step,flip
      integer :: j,k
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y
      ! if (1==step .and. 1==j*k) write(*,*) 'CUDA bc_bb_x] j,k:',j,k

      if (0<=j .and. j<=ny+1 .and. 0<=k .and. k<=nz_d+1) then
	      !rho1 = compute_rho(popsR_d,1,j,k,flip)

        popsR_d(1,j,k, 1, flip) = popsR_d(0, j,   k,    2, flip)
        popsR_d(1,j,k, 7, flip) = popsR_d(0, j-1, k,    8, flip)
        popsR_d(1,j,k,10, flip) = popsR_d(0, j+1, k,    9, flip)
        popsR_d(1,j,k,11, flip) = popsR_d(0, j  , k-1, 12, flip)
        popsR_d(1,j,k,14, flip) = popsR_d(0, j  , k+1, 13, flip)

	      !rho2 = compute_rho(popsR_d,nx,j,k,flip)

        popsR_d(nx,j,k, 2, flip) = popsR_d(nx+1, j,   k,    1, flip)
        popsR_d(nx,j,k, 8, flip) = popsR_d(nx+1, j+1, k,    7, flip)
        popsR_d(nx,j,k, 9, flip) = popsR_d(nx+1, j-1, k,   10, flip)
        popsR_d(nx,j,k,12, flip) = popsR_d(nx+1, j,   k+1, 11, flip)
        popsR_d(nx,j,k,13, flip) = popsR_d(nx+1, j,   k-1, 14, flip)
      endif
    end subroutine bc_bb_x


    attributes(global) subroutine bc_bb_y(step, flip)
      real :: rho1,rho2
      integer, value :: step,flip
      integer :: i,k
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y
      ! if (1==step .and. 1==i*k) write(*,*) 'CUDA bc_bb_y] i,k:',i,k

      if (0<=i .and. i<=nx+1 .and. 0<=k .and. k<=nz_d+1) then
          popsR_d(i,1,k, 3, flip) = popsR_d(i,   0, k,    4, flip)
          popsR_d(i,1,k, 7, flip) = popsR_d(i-1, 0, k,    8, flip)
          popsR_d(i,1,k, 9, flip) = popsR_d(i+1, 0, k,   10, flip)
          popsR_d(i,1,k,15, flip) = popsR_d(i,   0, k-1, 16, flip)
          popsR_d(i,1,k,18, flip) = popsR_d(i,   0, k+1, 17, flip)

          popsR_d(i,ny,k, 4, flip) = popsR_d(i,   ny+1, k,    3, flip)
          popsR_d(i,ny,k, 8, flip) = popsR_d(i+1, ny+1, k,    7, flip)
          popsR_d(i,ny,k,10, flip) = popsR_d(i-1, ny+1, k,    9, flip)
          popsR_d(i,ny,k,16, flip) = popsR_d(i,   ny+1, k+1, 15, flip)
          popsR_d(i,ny,k,17, flip) = popsR_d(i,   ny+1, k-1, 18, flip)
      endif
    end subroutine bc_bb_y


    attributes(global) subroutine bc_bb_z(step, flip)
      real :: rho1,rho2
      integer, value :: step,flip
      integer :: i,j
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIM + threadIdx%y
      ! if (1==step .and. 1==i*j) write(*,*) 'CUDA bc_bb_z] i,j:',i,j

      if (0<=j .and. j<=ny+1 .and. 0<=i .and. i<=nx+1) then
        !rho1 = compute_rho(popsR_d,i,j,1,flip)

        popsR_d(i,j,1, 5, flip) = popsR_d(i,   j,   0,  6, flip)
        popsR_d(i,j,1,11, flip) = popsR_d(i-1, j,   0, 12, flip) + p_d(12)*buz_d*rho1 * wall_z0_d(1)
        popsR_d(i,j,1,13, flip) = popsR_d(i+1, j,   0, 14, flip) - p_d(14)*buz_d*rho1 * wall_z0_d(1)
        popsR_d(i,j,1,15, flip) = popsR_d(i,   j-1, 0, 16, flip)
        popsR_d(i,j,1,17, flip) = popsR_d(i,   j+1, 0, 18, flip)

        !rho2 = compute_rho(popsR_d,i,j,nz_d,flip)

        popsR_d(i,j,nz_d, 6, flip) = popsR_d(i,   j,   nz_d+1,  5, flip)
        popsR_d(i,j,nz_d,12, flip) = popsR_d(i+1, j,   nz_d+1, 11, flip) - p_d(11)*buz_d*rho2 * wall_z1_d(1)
        popsR_d(i,j,nz_d,14, flip) = popsR_d(i-1, j,   nz_d+1, 13, flip) + p_d(13)*buz_d*rho2 * wall_z1_d(1)
        popsR_d(i,j,nz_d,16, flip) = popsR_d(i,   j+1, nz_d+1, 15, flip)
        popsR_d(i,j,nz_d,18, flip) = popsR_d(i,   j-1, nz_d+1, 17, flip)
      endif
    end subroutine bc_bb_z

    ! Bounce-Back BC 2 fluids
    attributes(global) subroutine bc_bb_x2(step, flip)
      integer, value :: step,flip
      integer :: j,k
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y
      ! if (1==step .and. 1==j*k) write(*,*) 'CUDA bc_bb_x] j,k:',j,k

      if (0<=j .and. j<=ny+1 .and. 0<=k .and. k<=nz_d+1) then	    
        popsR_d(1,j,k, 1, flip) = popsR_d(0, j,   k,    2, flip)
        popsR_d(1,j,k, 7, flip) = popsR_d(0, j-1, k,    8, flip)
        popsR_d(1,j,k,10, flip) = popsR_d(0, j+1, k,    9, flip)
        popsR_d(1,j,k,11, flip) = popsR_d(0, j  , k-1, 12, flip)
        popsR_d(1,j,k,14, flip) = popsR_d(0, j  , k+1, 13, flip)

        popsR_d(nx,j,k, 2, flip) = popsR_d(nx+1, j,   k,    1, flip)
        popsR_d(nx,j,k, 8, flip) = popsR_d(nx+1, j+1, k,    7, flip)
        popsR_d(nx,j,k, 9, flip) = popsR_d(nx+1, j-1, k,   10, flip)
        popsR_d(nx,j,k,12, flip) = popsR_d(nx+1, j,   k+1, 11, flip)
        popsR_d(nx,j,k,13, flip) = popsR_d(nx+1, j,   k-1, 14, flip)

        popsB_d(1,j,k, 1, flip) = popsB_d(0, j,   k,    2, flip)
        popsB_d(1,j,k, 7, flip) = popsB_d(0, j-1, k,    8, flip)
        popsB_d(1,j,k,10, flip) = popsB_d(0, j+1, k,    9, flip)
        popsB_d(1,j,k,11, flip) = popsB_d(0, j  , k-1, 12, flip)
        popsB_d(1,j,k,14, flip) = popsB_d(0, j  , k+1, 13, flip)

        popsB_d(nx,j,k, 2, flip) = popsB_d(nx+1, j,   k,    1, flip)
        popsB_d(nx,j,k, 8, flip) = popsB_d(nx+1, j+1, k,    7, flip)
        popsB_d(nx,j,k, 9, flip) = popsB_d(nx+1, j-1, k,   10, flip)
        popsB_d(nx,j,k,12, flip) = popsB_d(nx+1, j,   k+1, 11, flip)
        popsB_d(nx,j,k,13, flip) = popsB_d(nx+1, j,   k-1, 14, flip)
      endif
    end subroutine bc_bb_x2


    attributes(global) subroutine bc_bb_y2(step, flip)
      integer, value :: step,flip
      integer :: i,k
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y
      ! if (1==step .and. 1==i*k) write(*,*) 'CUDA bc_bb_y] i,k:',i,k

      if (0<=i .and. i<=nx+1 .and. 0<=k .and. k<=nz_d+1) then
          popsR_d(i,1,k, 3, flip) = popsR_d(i,   0, k,    4, flip)
          popsR_d(i,1,k, 7, flip) = popsR_d(i-1, 0, k,    8, flip)
          popsR_d(i,1,k, 9, flip) = popsR_d(i+1, 0, k,   10, flip)
          popsR_d(i,1,k,15, flip) = popsR_d(i,   0, k-1, 16, flip)
          popsR_d(i,1,k,18, flip) = popsR_d(i,   0, k+1, 17, flip)

          popsR_d(i,ny,k, 4, flip) = popsR_d(i,   ny+1, k,    3, flip)
          popsR_d(i,ny,k, 8, flip) = popsR_d(i+1, ny+1, k,    7, flip)
          popsR_d(i,ny,k,10, flip) = popsR_d(i-1, ny+1, k,    9, flip)
          popsR_d(i,ny,k,16, flip) = popsR_d(i,   ny+1, k+1, 15, flip)
          popsR_d(i,ny,k,17, flip) = popsR_d(i,   ny+1, k-1, 18, flip)

          popsB_d(i,1,k, 3, flip) = popsB_d(i,   0, k,    4, flip)
          popsB_d(i,1,k, 7, flip) = popsB_d(i-1, 0, k,    8, flip)
          popsB_d(i,1,k, 9, flip) = popsB_d(i+1, 0, k,   10, flip)
          popsB_d(i,1,k,15, flip) = popsB_d(i,   0, k-1, 16, flip)
          popsB_d(i,1,k,18, flip) = popsB_d(i,   0, k+1, 17, flip)

          popsB_d(i,ny,k, 4, flip) = popsB_d(i,   ny+1, k,    3, flip)
          popsB_d(i,ny,k, 8, flip) = popsB_d(i+1, ny+1, k,    7, flip)
          popsB_d(i,ny,k,10, flip) = popsB_d(i-1, ny+1, k,    9, flip)
          popsB_d(i,ny,k,16, flip) = popsB_d(i,   ny+1, k+1, 15, flip)
          popsB_d(i,ny,k,17, flip) = popsB_d(i,   ny+1, k-1, 18, flip)
      endif
    end subroutine bc_bb_y2


    attributes(global) subroutine bc_bb_z2(step, flip)
      real :: rho1,rho2, rho3,rho4
      integer, value :: step,flip
      integer :: i,j
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIM + threadIdx%y
      ! if (1==step .and. 1==i*j) write(*,*) 'CUDA bc_bb_z] i,j:',i,j

      if (0<=j .and. j<=ny+1 .and. 0<=i .and. i<=nx+1) then
        !rho1 = compute_rho(popsR_d,i,j,1,flip)
        popsR_d(i,j,1, 5, flip) = popsR_d(i,   j,   0,  6, flip)
        popsR_d(i,j,1,11, flip) = popsR_d(i-1, j,   0, 12, flip) + p_d(12)*buz_d*rho1 * wall_z0_d(1)
        popsR_d(i,j,1,13, flip) = popsR_d(i+1, j,   0, 14, flip) - p_d(14)*buz_d*rho1 * wall_z0_d(1)
        popsR_d(i,j,1,15, flip) = popsR_d(i,   j-1, 0, 16, flip)
        popsR_d(i,j,1,17, flip) = popsR_d(i,   j+1, 0, 18, flip)

        !rho2 = compute_rho(popsR_d,i,j,nz_d,flip)
        popsR_d(i,j,nz_d, 6, flip) = popsR_d(i,   j,   nz_d+1,  5, flip)
        popsR_d(i,j,nz_d,12, flip) = popsR_d(i+1, j,   nz_d+1, 11, flip) - p_d(11)*buz_d*rho2 * wall_z1_d(1)
        popsR_d(i,j,nz_d,14, flip) = popsR_d(i-1, j,   nz_d+1, 13, flip) + p_d(13)*buz_d*rho2 * wall_z1_d(1)
        popsR_d(i,j,nz_d,16, flip) = popsR_d(i,   j+1, nz_d+1, 15, flip)
        popsR_d(i,j,nz_d,18, flip) = popsR_d(i,   j-1, nz_d+1, 17, flip)

        !rho3 = compute_rho(popsB_d,i,j,1,flip)
        popsB_d(i,j,1, 5, flip) = popsB_d(i,   j,   0,  6, flip)
        popsB_d(i,j,1,11, flip) = popsB_d(i-1, j,   0, 12, flip) + p_d(12)*buz_d*rho3 * wall_z0_d(1)
        popsB_d(i,j,1,13, flip) = popsB_d(i+1, j,   0, 14, flip) - p_d(14)*buz_d*rho3 * wall_z0_d(1)
        popsB_d(i,j,1,15, flip) = popsB_d(i,   j-1, 0, 16, flip)
        popsB_d(i,j,1,17, flip) = popsB_d(i,   j+1, 0, 18, flip)

        !rho4 = compute_rho(popsB_d,i,j,nz_d,flip)
        popsB_d(i,j,nz_d, 6, flip) = popsB_d(i,   j,   nz_d+1,  5, flip)
        popsB_d(i,j,nz_d,12, flip) = popsB_d(i+1, j,   nz_d+1, 11, flip) - p_d(11)*buz_d*rho4 * wall_z1_d(1)
        popsB_d(i,j,nz_d,14, flip) = popsB_d(i-1, j,   nz_d+1, 13, flip) + p_d(13)*buz_d*rho4 * wall_z1_d(1)
        popsB_d(i,j,nz_d,16, flip) = popsB_d(i,   j+1, nz_d+1, 15, flip)
        popsB_d(i,j,nz_d,18, flip) = popsB_d(i,   j-1, nz_d+1, 17, flip)
      endif
    end subroutine bc_bb_z2

  end module kernels_fluid_cg
