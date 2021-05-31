#include "defines.h"

  module kernels_fluid
    use dimensions_m
    implicit none

    ! device arrays
    ! real(4), device :: popsR_d(0:nx+1,0:ny+1,0:nz+1, 0:npops-1, 2), popsB_d(0:nx+1,0:ny+1,0:nz+1, 0:npops-1, 2)
    ! real(4), device :: force_d(3, nx,ny,nz)
    ! real(4), device :: rhoR_d(1:nx,1:ny,1:nz),rhoB_d(1:nx,1:ny,1:nz)
    real(4), allocatable, device :: popsR_d(:,:,:,:,:), popsB_d(:,:,:,:,:)
    real(4), allocatable, device :: force_d(:, :,:,:)
    real(4), allocatable, device :: rhoR_d(:,:,:),rhoB_d(:,:,:)

    real(4), device :: pos_d(13, numAtoms,2)
    real(4), device :: forceAtoms_d(12, numAtoms,2)
    real(8), device :: myf_d(3, numAtoms), myt_d(3, numAtoms)
    real(8), device :: myf2_d(3, numAtoms), myt2_d(3, numAtoms)
    integer(2), device :: listAtoms_d(numAtoms * nx/TILE_DIMx* ny/TILE_DIMy* nz/TILE_DIMz)
    integer(2), device :: findAtoms_d(2,nx/TILE_DIMx, ny/TILE_DIMy, nz/TILE_DIMz)

    integer(1), allocatable, device :: issub_d(:, :, :)
    integer, device :: rmax_issub_d, sphereMax_d
    integer(1), allocatable, device :: myfluid_d(:,:,:, :)
    
    integer, device :: ndouble
    integer, device, allocatable :: exdouble(:), eydouble(:), ezdouble(:)

    integer, device :: stop_d, countmk_d, countrm_d, partVol_d, oldpartVol_d
    integer, device :: countn2p_d

    real, device :: minRho = -0.1, minPops = -0.1

    ! Debug
    real, device :: debugline_d(numAtoms, 21,4)
    real, device :: debugn2pf_d(2000,5)
    real, device :: debugn2pt_d(2000,5)
    integer(2), allocatable, device :: debugfluid_d(:,:,:)
   
  contains
    attributes(global) subroutine debug_issub
      integer :: i,j,k, id, val

      id = threadIdx%x

      do k = -rmax_issub_d,rmax_issub_d
        do j = -rmax_issub_d,rmax_issub_d
          do i = -rmax_issub_d,rmax_issub_d
            val = issub_d(i,j,k)
            write(*,*) 'CUDA-issub', i,j,k, val
          enddo
        enddo
      enddo

    end subroutine debug_issub

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

!!!!!!!!!!!!! Helpers for recovering phys vars
    attributes(device) function compute_rho(pops, i,j,k, flip)
     integer, intent(in) :: i,j,k, flip
     real, intent(in) :: pops(0:nx+1,0:ny+1,0:nz+1, 0:npops-1, 2)
     real :: compute_rho

     compute_rho = pops(i,j,k,0, flip) + pops(i,j,k,1, flip) + pops(i,j,k,2, flip) + &
       pops(i,j,k,3, flip) + pops(i,j,k,4, flip) + &
       pops(i,j,k,5, flip) + pops(i,j,k,6, flip) + pops(i,j,k,7, flip) + &
       pops(i,j,k,8, flip) + pops(i,j,k,9, flip) + &
       pops(i,j,k,10, flip) + pops(i,j,k,11, flip) + pops(i,j,k,12, flip) + &
       pops(i,j,k,13, flip) + pops(i,j,k,14, flip) + &
       pops(i,j,k,15, flip) + pops(i,j,k,16, flip) + &
       pops(i,j,k,17, flip) + pops(i,j,k,18, flip)
    end function compute_rho

    attributes(device) function compute_u_1fl(pops, i,j,k, flip, invrho)
     integer, intent(in) :: i,j,k, flip
     real, intent(in) :: pops(0:nx+1,0:ny+1,0:nz+1, 0:npops-1, 2), invrho
     real :: compute_u_1fl

     compute_u_1fl   = invrho * ( pops(i,j,k,1, flip) - pops(i,j,k,2, flip) + pops(i,j,k,7, flip) - &
       pops(i,j,k,8, flip) - pops(i,j,k,9, flip) + &
       pops(i,j,k,10, flip) + pops(i,j,k,11, flip) - pops(i,j,k,12, flip) - &
       pops(i,j,k,13, flip) + pops(i,j,k,14, flip) )
    end function compute_u_1fl

    attributes(device) function compute_v_1fl(pops, i,j,k, flip, invrho)
     integer, intent(in) :: i,j,k, flip
     real, intent(in) :: pops(0:nx+1,0:ny+1,0:nz+1, 0:npops-1, 2), invrho
     real :: compute_v_1fl

     compute_v_1fl    = invrho * ( pops(i,j,k,3, flip) - pops(i,j,k,4, flip) + pops(i,j,k,7, flip) - &
       pops(i,j,k,8, flip) + pops(i,j,k,9, flip) - &
       pops(i,j,k,10, flip) + pops(i,j,k,15, flip) - pops(i,j,k,16, flip) - &
       pops(i,j,k,17, flip) + pops(i,j,k,18, flip) )
    end function compute_v_1fl

    attributes(device) function compute_w_1fl(pops, i,j,k, flip, invrho)
     integer, intent(in) :: i,j,k, flip
     real, intent(in) :: pops(0:nx+1,0:ny+1,0:nz+1, 0:npops-1, 2), invrho
     real :: compute_w_1fl

     compute_w_1fl    = invrho * ( pops(i,j,k,5, flip) - pops(i,j,k,6, flip) + pops(i,j,k,11, flip) - &
       pops(i,j,k,12, flip) + pops(i,j,k,13, flip) - &
       pops(i,j,k,14, flip) + pops(i,j,k,15, flip) - pops(i,j,k,16, flip) + &
       pops(i,j,k,17, flip) - pops(i,j,k,18, flip) )

    end function compute_w_1fl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	 	END of Helpers for recovering phys vars


    ! Setup
    attributes(global) subroutine setup(vx,vy,vz)
      real, value :: vx,vy,vz
      integer :: i,j,k, l
      real    :: rho, u,v,w, eq, x,y,z
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z

      ! if (k>nz/2) then
      !   rho = 1.01
      ! else
      !   rho = 0.99
      ! endif

      x = 3.14159265358979 / nx * i
      y = 3.14159265358979 / ny * j
      z = 3.14159265358979 / nz * k
      rho = 1 + 0.1*sin(x)*sin(x) * sin(y)*sin(y) * sin(z)*sin(z)

      u = vx
      v = vy
      w = vz
      ! write(*,*) 'CUDA setup',rho,i,j,k
  
      rhoR_d(i,j,k) = rho
      do l = 0, npops-1
        eq = equil(rho, u,v,w, l)
        popsR_d(i,j,k, l, 1) = eq
      end do

      myfluid_d(i,j,k, 1) = fluid_fluid
      myfluid_d(i,j,k, 2) = fluid_fluid
    end subroutine setup

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

    attributes(global) subroutine setup2(vx,vy,vz)
      real, value :: vx,vy,vz
      integer :: i,j,k, l
      real    :: rhoR,rhoB, u,v,w, eqR,eqB, distx,disty,distz, rdist,tempr
      real,parameter    :: radius = nx/8
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z

      distx = i - (nx/2 + 0.5)
      disty = j - (ny/2 + 0.5)
      distz = k - (nz/2 + 0.5)
      rdist = sqrt(distx*distx + disty*disty + distz*distz)
      tempr = fcut(rdist, radius, radius+0.1)
      rhoR = tempr
      rhoB = 1.0 - tempr

      ! distx = i - (nx/2 + 0.5)
      ! disty = j - (ny/2 + 0.5)
      ! distz = k - (nz/4 + 0.5)
      ! rdist = sqrt(distx*distx + disty*disty + distz*distz)
      ! tempr = fcut(rdist, radius, radius+0.1)
      ! rhoR = rhoR + 0.0 + (1.0 - 0.0) * tempr
      ! rhoB = rhoB + 1.0 + (0.0 - 1.0) * tempr

      ! distx = i - (nx/2 + 0.5)
      ! disty = j - (ny/2 + 0.5)
      ! distz = k - (3*nz/4 + 0.5)
      ! rdist = sqrt(distx*distx + disty*disty + distz*distz)
      ! tempr = fcut(rdist, radius, radius+0.1)
      ! rhoR = rhoR + 0.0 + (1.0 - 0.0) * tempr
      ! rhoB = rhoB + 1.0 + (0.0 - 1.0) * tempr

      ! rhoR = 1.0
      ! rhoB = 0.0

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
      partVol_d = 0
      oldpartVol_d = 0
      countn2p_d = 0
    end subroutine setup2

#ifdef BGK
    ! Compute forces on lattice
    attributes(global) subroutine lb_force()
      integer :: i,j,k, l
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z

      force_d(1, i,j,k) = 0.0
      force_d(2, i,j,k) = 0.0
      force_d(3, i,j,k) = 0.0
    end subroutine


    ! Time Step
    attributes(global) subroutine time_step(step, flip, omega,oneminusomega)
      integer, value :: step,flip
      real, value :: omega,oneminusomega
      real    :: rho,invrho, u,v,w
      integer :: i,j,k, l, flop
      integer :: i1,j1,k1
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA time_step] i,k:',i,j,k

      rho = compute_rho(popsR_d,i,j,k,flip)
      if (rho < 1.0E-7) then
       write(*,*) 'Rho error:',rho,i,j,k
       stop_d = __LINE__
      endif

      invrho = 1.0 / rho
      u   = compute_u_1fl(popsR_d,i,j,k,flip, invrho)
      v   = compute_v_1fl(popsR_d,i,j,k,flip, invrho)
      w   = compute_w_1fl(popsR_d,i,j,k,flip, invrho)
  
      flop = 3 - flip
      do l = 0, npops-1
        i1 = i + ex_d(l)
        j1 = j + ey_d(l)
        k1 = k + ez_d(l)
        popsR_d(i1,j1,k1, l, flop) = popsR_d(i,j,k, l, flip)*oneminusomega + equil(rho, u,v,w, l)*omega
      end do
    end subroutine time_step

    attributes(global) subroutine time_step_force_cost(step, flip, omega,oneminusomega)
      integer, value :: step,flip
      real, value :: omega,oneminusomega
      real    :: rho,invrho, u,v,w, fx,fy,fz
      real    :: oldflu,feq,feqshift
      integer :: i,j,k, l, flop
      integer :: i1,j1,k1
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      ! if (1==step.and.1==i*j*k) write(*,*) 'CUDA time_step_force_cost] i,k:',i,j,k

      rho = compute_rho(popsR_d,i,j,k,flip)
      if (rho < 1.0E-7) then
       write(*,*) 'Rho error:',rho,i,j,k
       stop_d = __LINE__
      endif

      invrho = 1.0 / rho
      u   = compute_u_1fl(popsR_d,i,j,k,flip, invrho)
      v   = compute_v_1fl(popsR_d,i,j,k,flip, invrho)
      w   = compute_w_1fl(popsR_d,i,j,k,flip, invrho)

      fx = fx_d * invrho
      fy = fy_d * invrho
      fz = fz_d * invrho
      if (abs(fx) > 5.0E-1 .or. abs(fy) > 5.0E-1 .or. abs(fz) > 5.0E-1) then
       write(*,*) 'Big fx,fy,fz:',fx,fy,fz
       stop_d = __LINE__
      endif
  
      flop = 3 - flip
      do l = 0, npops-1
        i1 = i + ex_d(l)
        j1 = j + ey_d(l)
        k1 = k + ez_d(l)
      
        !!!!!
        ! THIS is not like Krueger...same as ML
        !!!!!
        oldflu   = popsR_d(i,j,k, l, flip)
        feq      = equil(rho, u,   v,   w, l)
        feqshift = equil(rho, u+fx,v+fy,w+fz, l)

        popsR_d(i1,j1,k1, l, flop) = oldflu*oneminusomega  - feq*oneminusomega + feqshift
      end do
    end subroutine time_step_force_cost

    attributes(global) subroutine time_step_force(step, flip, omega,oneminusomega)
      integer, value :: step,flip
      real, value :: omega,oneminusomega
      real    :: rho,invrho, u,v,w, fx,fy,fz
      real    :: oldflu,feq,feqshift
      integer :: i,j,k, l, flop
      integer :: i1,j1,k1
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      ! if (1==step.and.1==i*j*k) write(*,*) 'CUDA time_step_force] i,k:',i,j,k

      rho = compute_rho(popsR_d,i,j,k,flip)
      if (rho < 1.0E-7) then
       write(*,*) 'Rho error:',rho,i,j,k
       stop_d = __LINE__
      endif

      invrho = 1.0 / rho
      u   = compute_u_1fl(popsR_d,i,j,k,flip, invrho)
      v   = compute_v_1fl(popsR_d,i,j,k,flip, invrho)
      w   = compute_w_1fl(popsR_d,i,j,k,flip, invrho)

      fx = force_d(1, i,j,k) * invrho
      fy = force_d(2, i,j,k) * invrho
      fz = force_d(3, i,j,k) * invrho
      if (abs(fx) > 5.0E-1 .or. abs(fy) > 5.0E-1 .or. abs(fz) > 5.0E-1) then
       write(*,*) 'Big fx,fy,fz:',fx,fy,fz
       stop_d = __LINE__
      endif
  
      flop = 3 - flip
      do l = 0, npops-1
        i1 = i + ex_d(l)
        j1 = j + ey_d(l)
        k1 = k + ez_d(l)
      
        !!!!!
        ! THIS is not like Krueger...same as ML
        !!!!!
        oldflu   = popsR_d(i,j,k, l, flip)
        feq      = equil(rho, u,   v,   w, l)
        feqshift = equil(rho, u+fx,v+fy,w+fz, l)
        popsR_d(i1,j1,k1, l, flop) = oldflu*oneminusomega - feq*oneminusomega + feqshift
      end do
    end subroutine time_step_force
#endif

    attributes(global) subroutine init_rho_isfluid_BGK(step, flip)
    integer, value :: step,flip
    real    :: rhoR
    integer :: i,j,k, l

    i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
    j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
    k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
    ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA init_rho_isfluid_BGK] i,k:',i,j,k

    ! Only on fluid nodes
    if (myfluid_d(i,j,k, flip) == fluid_fluid) then
      rhoR = compute_rho(popsR_d,i,j,k,flip)

      if (rhoR<minRho .or. rhoR>10) then
        ! write(*,*) 'init_rho_isfluid_BGK]Range error rhoR', step, linear(i,j,k), rhoR
        stop_d = __LINE__
      endif

      rhoR_d(i,j,k) = rhoR
    else
      rhoR_d(i,j,k) = MINDENS
    endif
    end subroutine init_rho_isfluid_BGK


    attributes(global) subroutine time_step_BGK(step, flip, omega,oneminusomega)
      integer, value :: step,flip
      real, value :: omega,oneminusomega
      real    :: rhoR, invrho, u,v,w  
      integer :: i,j,k, l
      
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA time_step_CG] i,k:',i,j,k

      if (myfluid_d(i,j,k, flip) == fluid_fluid) then
        rhoR = rhoR_d(i,j,k)

        invrho = 1.0 / rhoR
        u   = compute_u_1fl(popsR_d,i,j,k,flip, invrho)
        v   = compute_v_1fl(popsR_d,i,j,k,flip, invrho)
        w   = compute_w_1fl(popsR_d,i,j,k,flip, invrho)

        !bgk step
        do l = 0, npops-1
          popsR_d(i,j,k,l, flip) = popsR_d(i,j,k, l, flip)*oneminusomega + equil(rhoR,u,v,w, l)*omega
        enddo
      endif
    end subroutine time_step_BGK

    attributes(global) subroutine copy_BGK(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l, flop

      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA copy_BGK] i,k:',i,j,k

      flop = 3 - flip
      
      do l = 0, npops-1
        popsR_d(i,j,k,l, flop) = popsR_d(i,j,k,l, flip)
      enddo
    end subroutine copy_BGK


    attributes(global) subroutine stream_BGK(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l, flop
      integer :: i1,j1,k1


      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA stream_BGK] i,k:',i,j,k

      flop = 3 - flip

      ! Streaming      
      if (myfluid_d(i,j,k, flip) /= fluid_spheredead) then
        do l = 1, npops-1
          i1 = i + ex_d(l)
          j1 = j + ey_d(l)
          k1 = k + ez_d(l)
          popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
        enddo
      endif    
    end subroutine stream_BGK

    attributes(global) subroutine stream_BGK_x(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l, flop
      integer :: i1,j1,k1
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1

      if (j>ny+1) return
      if (k>nz+1) return

      flop = 3 - flip

      ! Stream x=0
      i = 0
      if (myfluid_d(i,j,k, flip) /= fluid_spheredead) then
        do l = 1, npops-1
          i1 = i + ex_d(l)
          j1 = j + ey_d(l)
          k1 = k + ez_d(l)
          if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz) then
            popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
          endif
        end do
      endif

      ! Stream x=nx+1
      i = nx+1
      if (myfluid_d(i,j,k, flip) /= fluid_spheredead) then
        do l = 1, npops-1
          i1 = i + ex_d(l)
          j1 = j + ey_d(l)
          k1 = k + ez_d(l)
          if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz) then
            popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
          endif
        end do
      endif
    end subroutine stream_BGK_x

    attributes(global) subroutine stream_BGK_y(step, flip)
    integer, value :: step,flip
    integer :: i,j,k, l, flop
    integer :: i1,j1,k1

    i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
    k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
    
    if (i>nx+1) return
    if (k>nz+1) return

    flop = 3 - flip

    ! Stream y=0
    j = 0
    if (myfluid_d(i,j,k, flip) /= fluid_spheredead) then
      do l = 1, npops-1
        i1 = i + ex_d(l)
        j1 = j + ey_d(l)
        k1 = k + ez_d(l)
        if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz) then
          popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
        endif
      end do
    endif

    ! Stream y=ny+1
    j = ny+1
    if (myfluid_d(i,j,k, flip) /= fluid_spheredead) then
      do l = 1, npops-1
        i1 = i + ex_d(l)
        j1 = j + ey_d(l)
        k1 = k + ez_d(l)
        if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz) then
          popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
        endif
      end do
    endif
  end subroutine stream_BGK_y

  attributes(global) subroutine stream_BGK_z(step, flip)
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
    if (myfluid_d(i,j,k, flip) /= fluid_spheredead) then
      do l = 1, npops-1
        i1 = i + ex_d(l)
        j1 = j + ey_d(l)
        k1 = k + ez_d(l)
        if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz) then
          popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
        endif
      end do
    endif

    ! Stream z=nz+1
    k = nz+1
    if (myfluid_d(i,j,k, flip) /= fluid_spheredead) then
      do l = 1, npops-1
        i1 = i + ex_d(l)
        j1 = j + ey_d(l)
        k1 = k + ez_d(l)
        if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz) then
          popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
        endif
      end do
    endif
  end subroutine stream_BGK_z    

#ifdef CG
    !!!!!!!!!!!!!!!!!!! Color gradient !!!!!!!!!!!!!!!!!!!
    attributes(global) subroutine init_isfluid_CG(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l
      integer :: atm_i,atm_j,atm_k, atm_st,atm_en, i_atm,i_list
      logical :: updated
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA init_rho_isfluid_CG] i,k:',i,j,k

      ! Set myfluid_d
      if (withParticles) then
        atm_i = blockIdx%x
        atm_j = blockIdx%y
        atm_k = blockIdx%z
        atm_st = findAtoms_d(1, atm_i,atm_j,atm_k)
        atm_en = findAtoms_d(2, atm_i,atm_j,atm_k)
        ! if (step==1) write(*,*) 'findAtms',linear(i,j,k), atm_st,atm_en, linear(atm_i,atm_j,atm_k)

        ! Update myfluid_d
        updated = .false.
        debugfluid_d(i,j,k) = 0
        if (atm_st>0) then
          do i_list = atm_st, atm_en
            i_atm = listAtoms_d(i_list)
            ! if (step==1) write(*,*) 'CUDA-init_rho_isfluid_CG',linear(i,j,k), i_atm

            if (.not. updated) then
              updated = setmyfluid(i,j,k, i_atm, step, flip)
            endif
          enddo
        endif
        if (.not. updated) myfluid_d(i,j,k, flip) = fluid_fluid
      endif
    end subroutine init_isfluid_CG


    attributes(global) subroutine del_fluid_CG(step, flip)
      integer, value :: step,flip
      integer :: i,j,k
      integer :: atm_i,atm_j,atm_k, atm_st,atm_en, i_atm,i_list
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z

      if (withParticles) then
        atm_i = blockIdx%x 
        atm_j = blockIdx%y 
        atm_k = blockIdx%z 
        atm_st = findAtoms_d(1, atm_i,atm_j,atm_k)
        atm_en = findAtoms_d(2, atm_i,atm_j,atm_k)

        ! del_fluid 1st (old rho), if needed
        if (atm_st>0) then
          do i_list = atm_st, atm_en
            i_atm = listAtoms_d(i_list)

            ! Add/delete fluid if particle moved
            if (particle_moved(i_atm) == 1) call del_fluid(step, flip, i_atm,i,j,k)
          enddo
        endif
      endif
    end subroutine del_fluid_CG

    attributes(global) subroutine init_rho_CG(step, flip)
      integer, value :: step,flip
      real    :: rhoR,rhoB
      integer :: i,j,k
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z

      ! Only on fluid nodes
      if (myfluid_d(i,j,k, flip) == fluid_fluid) then
        rhoR = compute_rho(popsR_d,i,j,k,flip)
        rhoB = compute_rho(popsB_d,i,j,k,flip)

        if (rhoR<minRho .or. rhoR>10) then
          write(*,*) 'init_rho_CG]Range error rhoR', step, linear(i,j,k), rhoR
          stop_d = __LINE__
        endif
        if (rhoB<minRho .or. rhoB>10) then
          write(*,*) 'init_rho_CG]Range error rhoB', step, linear(i,j,k), rhoB
          stop_d = __LINE__
        endif

        rhoR_d(i,j,k) = rhoR
        rhoB_d(i,j,k) = rhoB
      else
        rhoR_d(i,j,k) = MINDENS
        rhoB_d(i,j,k) = MINDENS
      endif
    end subroutine init_rho_CG

    attributes(global) subroutine make_fluid_CG(step, flip)
      integer, value :: step,flip
      integer :: i,j,k
      integer :: atm_i,atm_j,atm_k, atm_st,atm_en, i_atm,i_list
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z

      if (withParticles) then
        atm_i = blockIdx%x 
        atm_j = blockIdx%y 
        atm_k = blockIdx%z 
        atm_st = findAtoms_d(1, atm_i,atm_j,atm_k)
        atm_en = findAtoms_d(2, atm_i,atm_j,atm_k)

        ! make_fluid, if needed
        if (atm_st>0) then
          do i_list = atm_st, atm_en
            i_atm = listAtoms_d(i_list)

            ! Add/delete fluid if particle moved
            if (particle_moved(i_atm) == 1) call make_fluid(step, flip, i_atm,i,j,k)
          enddo
        endif
      endif
    end subroutine make_fluid_CG

    attributes(device) pure function getDeltai(i, xxx)
    integer, intent(in) :: i, xxx
    integer :: getDeltai
    
      getDeltai = i - xxx
      if (xperiodic) then
        if (xxx-rdim < 1 .or. xxx+rdim > nx) then
          if (xxx+rdim > nx .and. i < 2*rdim) then
            getDeltai = i - xxx + nx
          else if (xxx-rdim < 1 .and. i > nx - 2*rdim) then
            getDeltai = i - xxx - nx
          endif
        endif
      endif
    end function getDeltai

    attributes(device) pure function getDeltaj(j, yyy)
    integer, intent(in) :: j, yyy
    integer :: getDeltaj
    
      getDeltaj = j - yyy
      if (yperiodic) then
        if (yyy-rdim < 1 .or. yyy+rdim > ny) then
          if (yyy+rdim > ny .and. j < 2*rdim) then
            getDeltaj = j - yyy + ny
          else if (yyy-rdim < 1 .and. j > ny - 2*rdim) then
            getDeltaj = j - yyy - ny
          endif
        endif
      endif
    end function getDeltaj

    attributes(device) pure function getDeltak(k, zzz)
    integer, intent(in) :: k, zzz
    integer :: getDeltak
    
      getDeltak = k - zzz
      if (zperiodic) then
        if (zzz-rdim < 1 .or. zzz+rdim > nz) then
          if (zzz+rdim > nz .and. k < 2*rdim) then
            getDeltak = k - zzz + nz
          else if (zzz-rdim < 1 .and. k > nz - 2*rdim) then
            getDeltak = k - zzz - nz
          endif
        endif
      endif
    end function getDeltak

    attributes(device) function setmyfluid(i,j,k, iatm, step, flip)
      integer, intent(in) :: i,j,k, iatm, step, flip
      logical :: setmyfluid
      integer :: deltai,deltaj,deltak
      integer :: xxx,yyy,zzz
      integer(1) :: val
      integer    :: tempInt
  
      setmyfluid = .false.

      xxx = nint(pos_d(1,iatm, flip))
      yyy = nint(pos_d(2,iatm, flip))
      zzz = nint(pos_d(3,iatm, flip))

      deltai = getDeltai(i, xxx)
      deltaj = getDeltaj(j, yyy)
      deltak = getDeltak(k, zzz)
      
      ! Check if inside particle      
      if (abs(deltai)<=rmax_issub_d.and.abs(deltaj)<=rmax_issub_d.and.abs(deltak)<=rmax_issub_d) then
        
        val = issub_d( deltai,deltaj,deltak )
        if (val /= fluid_fluid) then
          setmyfluid = .true.
          myfluid_d(i,j,k, flip) = val
          tempInt = atomicadd( partVol_d, 1)
          debugfluid_d(i,j,k) = iatm
          ! write(*,*) 'isfl',step, i,j,k, myfluid_d(i,j,k, flip)+0
        endif
      endif
    end function setmyfluid


    attributes(global) subroutine compute_densities_wall(step, flip)
      integer, value :: step,flip
      real    :: rhoR,rhoB, rhosum,invrho,  Rsum,Bsum, DsumR
      integer :: Dsum
      integer :: i,j,k
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
    
      if(myfluid_d(i,j,k, flip)==fluid_wall .or. myfluid_d(i,j,k, flip)==fluid_spherelist) then
        call compute_onebelt_density_2fl_weights(i,j,k, flip, Rsum, Bsum, DsumR)

        if(DsumR==ZERO .and. myfluid_d(i,j,k, flip)==fluid_spherelist) then
          call compute_secbelt_density_twofluids(i,j,k, flip, Rsum, Bsum,Dsum)

          if(Dsum==ZERO) then
            call fix_onebelt_density_twofluids(i,j,k, flip, Rsum, Bsum, Dsum)
            ! errorCount = errorCount + nint(Dsum)
          else
            ! secbelt /= 0
            rhoR_d(i,j,k) = Rsum
            rhoB_d(i,j,k) = Bsum
          endif
        else
          ! onebelt /= 0
          rhoR_d(i,j,k) = Rsum
          rhoB_d(i,j,k) = Bsum
        endif
      endif
    end subroutine compute_densities_wall

    ! pos(1:3) -> xxx,yyy,zzz
    ! pos(4:6) -> vxx,vyy,vzz
    ! pos(7:9) -> oxx,oyy,ozz
    ! pos(10:13) -> q0,q1,q2,q3
    ! force(1:3) -> fxx,fyy,fzz
    ! force(4:6) -> fxb,fyb,fzb
    ! force(7:9) -> tqx,tqy,tqz
    ! force(10:12) -> tqb,tqb,tqb

    attributes(device) pure function particle_moved(iatm)
    integer, intent(in) :: iatm
    integer :: particle_moved
    integer :: xxx,yyy,zzz, xxo,yyo,zzo

      xxx = nint(pos_d(1,iatm, 1))
      xxo = nint(pos_d(1,iatm, 2))
      yyy = nint(pos_d(2,iatm, 1))
      yyo = nint(pos_d(2,iatm, 2))    
      zzz = nint(pos_d(3,iatm, 1))    
      zzo = nint(pos_d(3,iatm, 2))

      if ( (xxx /= xxo) .or. (yyy /= yyo) .or. (zzz /= zzo) ) then
        particle_moved = 1
      else
        particle_moved = 0
      endif
    end function particle_moved

    
    attributes(device) subroutine del_fluid(step, flip, iatm, i,j,k )
      integer, value    :: step, flip
      integer, value    :: iatm, i,j,k
      logical           :: lmove
      real              :: rtemp(3), ftemp(3),ttemp(3), modr, tempFoo
      real              :: rhoR,rhoB,rhosum,invrho, myu,myv,myw
      real              :: aaa,bbb,ccc
      real(8)           :: acc(3), accFoo1,accFoo2,accFoo3
      integer           :: Dsum, tempInt
      integer :: deltai,deltaj,deltak, atmDist
      integer :: xxx,yyy,zzz

      ! if (1==step) write(*,*) 'CUDA-mk-del-fl]', linear(i,j,k), iatm

      ! apply periodic conditions if necessary
      ! i=pimage(ixpbc,i,nx)
      ! j=pimage(iypbc,j,ny)
      ! k=pimage(izpbc,k,nz)

      xxx = nint(pos_d(1,iatm, flip))
      yyy = nint(pos_d(2,iatm, flip))
      zzz = nint(pos_d(3,iatm, flip))

      deltai = getDeltai(i, xxx)
      deltaj = getDeltaj(j, yyy)
      deltak = getDeltak(k, zzz)

      atmDist = deltai*deltai + deltaj*deltaj + deltak*deltak
      ! only remove fluid if near THIS atom...
      if (atmDist > sphereMax_d) return


      ! On SphereList points...
        
      !!!!!!!! Delete Fluid
      if(myfluid_d(i,j,k, 3-flip)==fluid_fluid .and. myfluid_d(i,j,k, flip) == fluid_spherelist) then
        !fluid node is trasformed to solid
        !formula taken from eq. 18 of PRE 83, 046707 (2011)
        rhoR = rhoR_d(i,j,k)
        rhoB = rhoB_d(i,j,k)
        rhosum = rhoR + rhoB
        invrho = 1.0 / rhosum

        myu   = compute_u_1fl(popsR_d,i,j,k,flip, invrho) + compute_u_1fl(popsB_d,i,j,k,flip, invrho)
        myv   = compute_v_1fl(popsR_d,i,j,k,flip, invrho) + compute_v_1fl(popsB_d,i,j,k,flip, invrho)
        myw   = compute_w_1fl(popsR_d,i,j,k,flip, invrho) + compute_w_1fl(popsB_d,i,j,k,flip, invrho)

        ftemp(1) = (rhoR + rhoB) * myu
        ftemp(2) = (rhoR + rhoB) * myv
        ftemp(3) = (rhoR + rhoB) * myw
        
        ! Double
        acc(1) = ftemp(1)
        acc(2) = ftemp(2)
        acc(3) = ftemp(3)
        accFoo1 = atomicadd( myf2_d(1,iatm), acc(1) )
        accFoo2 = atomicadd( myf2_d(2,iatm), acc(2) )
        accFoo3 = atomicadd( myf2_d(3,iatm), acc(3) )

        ! naive atomic
        tempFoo = atomicadd( forceAtoms_d(1,iatm,flip), ftemp(1) )
        tempFoo = atomicadd( forceAtoms_d(2,iatm,flip), ftemp(2) )
        tempFoo = atomicadd( forceAtoms_d(3,iatm,flip), ftemp(3) )

#ifdef DEBUG_MKRM        
        tempInt = atomicadd( countrm_d, 1)
        write(*,*) 'rm-fl]', linear(i,j,k), ftemp(1),ftemp(2),ftemp(3)
        ! write(*,*) 'myrm]', linear(i,j,k), myu,myv,myw
        ! write(*,*) 'myrm]', linear(i,j,k), rhoR,rhoB
#endif        

        if (lrotate) then
          rtemp(1) = i - pos_d(1,iatm,flip)
          rtemp(2) = j - pos_d(2,iatm,flip)
          rtemp(3) = k - pos_d(3,iatm,flip)

          aaa = ONE/real(nx)
          bbb = ONE/real(ny)
          ccc = ONE/real(nz)
          
          rtemp(1) = rtemp(1) - real(nx)*nint(aaa*rtemp(1))
          rtemp(2) = rtemp(2) - real(ny)*nint(bbb*rtemp(2))
          rtemp(3) = rtemp(3) - real(nz)*nint(ccc*rtemp(3))
          
          modr = sqrt( rtemp(1)*rtemp(1) + rtemp(2)*rtemp(2) + rtemp(3)*rtemp(3))
          rtemp(1:3) = rdim / modr * rtemp(1:3)
          ttemp(1) = xcross(rtemp,ftemp)
          ttemp(2) = ycross(rtemp,ftemp)
          ttemp(3) = zcross(rtemp,ftemp)

          ! Double
          acc(1) = ttemp(1)
          acc(2) = ttemp(2)
          acc(3) = ttemp(3)
          accFoo1 = atomicadd( myt2_d(1,iatm), acc(1) )
          accFoo2 = atomicadd( myt2_d(2,iatm), acc(2) )
          accFoo3 = atomicadd( myt2_d(3,iatm), acc(3) )

          ! naive atomic
          tempFoo = atomicadd( forceAtoms_d(7,iatm,flip), ttemp(1) )
          tempFoo = atomicadd( forceAtoms_d(8,iatm,flip), ttemp(2) )
          tempFoo = atomicadd( forceAtoms_d(9,iatm,flip), ttemp(3) )
        endif
      endif
        
      ! On SphereDead points...
      
      !!!!!!!! Delete Fluid
      if(myfluid_d(i,j,k, 3-flip)==fluid_fluid .and. myfluid_d(i,j,k, flip)==fluid_spheredead) then
        !fluid node is trasformed to solid
        !formula taken from eq. 18 of PRE 83, 046707 (2011)
        
        rhoR = rhoR_d(i,j,k)
        rhoB = rhoB_d(i,j,k)
        rhosum = rhoR + rhoB
        invrho = 1.0 / rhosum

        myu   = compute_u_1fl(popsR_d,i,j,k,flip, invrho) + compute_u_1fl(popsB_d,i,j,k,flip, invrho)
        myv   = compute_v_1fl(popsR_d,i,j,k,flip, invrho) + compute_v_1fl(popsB_d,i,j,k,flip, invrho)
        myw   = compute_w_1fl(popsR_d,i,j,k,flip, invrho) + compute_w_1fl(popsB_d,i,j,k,flip, invrho)

        ftemp(1) = (rhoR + rhoB) * myu
        ftemp(2) = (rhoR + rhoB) * myv
        ftemp(3) = (rhoR + rhoB) * myw

        ! Double
        acc(1) = ftemp(1)
        acc(2) = ftemp(2)
        acc(3) = ftemp(3)
        accFoo1 = atomicadd( myf2_d(1,iatm), acc(1) )
        accFoo2 = atomicadd( myf2_d(2,iatm), acc(2) )
        accFoo3 = atomicadd( myf2_d(3,iatm), acc(3) )
        
        ! naive atomic
        tempFoo = atomicadd( forceAtoms_d(1,iatm,flip), ftemp(1) )
        tempFoo = atomicadd( forceAtoms_d(2,iatm,flip), ftemp(2) )
        tempFoo = atomicadd( forceAtoms_d(3,iatm,flip), ftemp(3) )

#ifdef DEBUG_MKRM
        tempInt = atomicadd( countrm_d, 1)
        write(*,*) 'rm-fl2]', linear(i,j,k), ftemp(1),ftemp(2),ftemp(3)
        ! write(*,*) 'myrm2]', linear(i,j,k), myu,myv,myw
        ! write(*,*) 'myrm2]', linear(i,j,k), rhoR,rhoB
#endif

        if (lrotate) then
          rtemp(1) = i - pos_d(1,iatm,flip)
          rtemp(2) = j - pos_d(2,iatm,flip)
          rtemp(3) = k - pos_d(3,iatm,flip)
          modr = sqrt( rtemp(1)*rtemp(1) + rtemp(2)*rtemp(2) + rtemp(3)*rtemp(3))
          rtemp(1:3) = rdim / modr * rtemp(1:3)
          ttemp(1) = xcross(rtemp,ftemp)
          ttemp(2) = ycross(rtemp,ftemp)
          ttemp(3) = zcross(rtemp,ftemp)

          ! Double
          acc(1) = ttemp(1)
          acc(2) = ttemp(2)
          acc(3) = ttemp(3)
          accFoo1 = atomicadd( myt2_d(1,iatm), acc(1) )
          accFoo2 = atomicadd( myt2_d(2,iatm), acc(2) )
          accFoo3 = atomicadd( myt2_d(3,iatm), acc(3) )

          ! naive atomic
          tempFoo = atomicadd( forceAtoms_d(7,iatm,flip), ttemp(1) )
          tempFoo = atomicadd( forceAtoms_d(8,iatm,flip), ttemp(2) )
          tempFoo = atomicadd( forceAtoms_d(9,iatm,flip), ttemp(3) )
        endif
      endif    
    end subroutine del_fluid

    attributes(device) subroutine make_fluid(step, flip, iatm, i,j,k )
      integer, value    :: step, flip
      integer, value    :: iatm, i,j,k
      logical           :: lmove
      real              :: rtemp(3), ftemp(3),ttemp(3), modr, tempFoo
      real              :: rhoR,rhoB,rhosum,invrho, myu,myv,myw, Rsum, Bsum
      real              :: aaa,bbb,ccc
      real(8)           :: acc(3), accFoo1,accFoo2,accFoo3
      integer           :: Dsum, tempInt
      integer :: deltai,deltaj,deltak, atmDist
      integer :: xxx,yyy,zzz

      ! if (1==step) write(*,*) 'CUDA-mk-del-fl]', linear(i,j,k), iatm

      ! Here i,j,k are already periodic-wrapped

      xxx = nint(pos_d(1,iatm, flip))
      yyy = nint(pos_d(2,iatm, flip))
      zzz = nint(pos_d(3,iatm, flip))

      deltai = getDeltai(i, xxx)
      deltaj = getDeltaj(j, yyy)
      deltak = getDeltak(k, zzz)

      atmDist = deltai*deltai + deltaj*deltaj + deltak*deltak
      ! only remove fluid if near THIS atom...
      if (atmDist > sphereMax_d) return


      ! On SphereList points...
      if(myfluid_d(i,j,k, 3-flip)==fluid_spherelist .and. myfluid_d(i,j,k, flip) == fluid_fluid) then
        call compute_onebelt_density_twofluids(i,j,k, flip, Rsum, Bsum, Dsum)
        if(Dsum==ZERO)then
          call compute_secbelt_density_twofluids(i,j,k, flip, Rsum, Bsum,Dsum)
          if(Dsum==ZERO)then
            call fix_onebelt_density_twofluids(i,j,k, flip, Rsum, Bsum, Dsum)
            ! errorCount = errorCount + nint(Dsum)
          endif
        endif

        myu = pos_d(4,iatm,flip)
        myv = pos_d(5,iatm,flip)
        myw = pos_d(6,iatm,flip)

        !formula taken from eq. 25 of PRE 83, 046707 (2011)
        call initialize_newnode_fluids(i,j,k, Rsum,Bsum, myu,myv,myw, flip)

        !formula taken from eq. 26 of PRE 83, 046707 (2011)
        ftemp(1) = -(Rsum+Bsum) * myu
        ftemp(2) = -(Rsum+Bsum) * myv
        ftemp(3) = -(Rsum+Bsum) * myw

        ! Double
        acc(1) = ftemp(1)
        acc(2) = ftemp(2)
        acc(3) = ftemp(3)
        accFoo1 = atomicadd( myf2_d(1,iatm), acc(1) )
        accFoo2 = atomicadd( myf2_d(2,iatm), acc(2) )
        accFoo3 = atomicadd( myf2_d(3,iatm), acc(3) )

        ! naive atomic
        tempFoo = atomicadd( forceAtoms_d(1,iatm,flip), ftemp(1) )
        tempFoo = atomicadd( forceAtoms_d(2,iatm,flip), ftemp(2) )
        tempFoo = atomicadd( forceAtoms_d(3,iatm,flip), ftemp(3) )

#ifdef DEBUG_MKRM
        tempInt = atomicadd( countmk_d, 1)
        write(*,*) 'mk-fl]', linear(i,j,k), ftemp(1),ftemp(2),ftemp(3)
#endif        

        if (lrotate) then
          rtemp(1) = i - pos_d(1,iatm,flip)
          rtemp(2) = j - pos_d(2,iatm,flip)
          rtemp(3) = k - pos_d(3,iatm,flip)

          aaa = ONE/real(nx)
          bbb = ONE/real(ny)
          ccc = ONE/real(nz)
          
          rtemp(1) = rtemp(1) - real(nx)*nint(aaa*rtemp(1))
          rtemp(2) = rtemp(2) - real(ny)*nint(bbb*rtemp(2))
          rtemp(3) = rtemp(3) - real(nz)*nint(ccc*rtemp(3))

          modr = sqrt( rtemp(1)*rtemp(1) + rtemp(2)*rtemp(2) + rtemp(3)*rtemp(3))
          rtemp(1:3) = rdim / modr * rtemp(1:3)
          ttemp(1) = xcross(rtemp,ftemp)
          ttemp(2) = ycross(rtemp,ftemp)
          ttemp(3) = zcross(rtemp,ftemp)

          ! Double
          acc(1) = ttemp(1)
          acc(2) = ttemp(2)
          acc(3) = ttemp(3)
          accFoo1 = atomicadd( myt2_d(1,iatm), acc(1) )
          accFoo2 = atomicadd( myt2_d(2,iatm), acc(2) )
          accFoo3 = atomicadd( myt2_d(3,iatm), acc(3) )

          ! naive atomic
          tempFoo = atomicadd( forceAtoms_d(7,iatm,flip), ttemp(1) )
          tempFoo = atomicadd( forceAtoms_d(8,iatm,flip), ttemp(2) )
          tempFoo = atomicadd( forceAtoms_d(9,iatm,flip), ttemp(3) )
        endif
      endif

      
      ! On SphereDead points...
      if(myfluid_d(i,j,k, 3-flip)==fluid_spheredead .and. myfluid_d(i,j,k, flip)==fluid_fluid) then
        call compute_onebelt_density_twofluids(i,j,k, flip, Rsum, Bsum, Dsum)
        if(Dsum==ZERO)then
          call compute_secbelt_density_twofluids(i,j,k, flip, Rsum, Bsum,Dsum)
          if(Dsum==ZERO)then
            call fix_onebelt_density_twofluids(i,j,k, flip, Rsum, Bsum, Dsum)
            ! errorCount = errorCount + nint(Dsum)
          endif
        endif

        myu = pos_d(4,iatm,flip)
        myv = pos_d(5,iatm,flip)
        myw = pos_d(6,iatm,flip)
        
        !formula taken from eq. 25 of PRE 83, 046707 (2011)
        call initialize_newnode_fluids(i,j,k, Rsum,Bsum, myu,myv,myw, flip)

        !formula taken from eq. 26 of PRE 83, 046707 (2011)
        ftemp(1) = -(Rsum+Bsum) * myu
        ftemp(2) = -(Rsum+Bsum) * myv
        ftemp(3) = -(Rsum+Bsum) * myw

        ! Double
        acc(1) = ftemp(1)
        acc(2) = ftemp(2)
        acc(3) = ftemp(3)
        accFoo1 = atomicadd( myf2_d(1,iatm), acc(1) )
        accFoo2 = atomicadd( myf2_d(2,iatm), acc(2) )
        accFoo3 = atomicadd( myf2_d(3,iatm), acc(3) )

        ! naive atomic
        tempFoo = atomicadd( forceAtoms_d(1,iatm,flip), ftemp(1) )
        tempFoo = atomicadd( forceAtoms_d(2,iatm,flip), ftemp(2) )
        tempFoo = atomicadd( forceAtoms_d(3,iatm,flip), ftemp(3) )

#ifdef DEBUG_MKRM
        tempInt = atomicadd( countmk_d, 1)
        write(*,*) 'mk-fl2]', linear(i,j,k), ftemp(1),ftemp(2),ftemp(3)
#endif        

        if (lrotate) then
          rtemp(1) = i - pos_d(1,iatm,flip)
          rtemp(2) = j - pos_d(2,iatm,flip)
          rtemp(3) = k - pos_d(3,iatm,flip)
          modr = sqrt( rtemp(1)*rtemp(1) + rtemp(2)*rtemp(2) + rtemp(3)*rtemp(3))
          rtemp(1:3) = rdim / modr * rtemp(1:3)
          ttemp(1) = xcross(rtemp,ftemp)
          ttemp(2) = ycross(rtemp,ftemp)
          ttemp(3) = zcross(rtemp,ftemp)

          ! Double
          acc(1) = ttemp(1)
          acc(2) = ttemp(2)
          acc(3) = ttemp(3)
          accFoo1 = atomicadd( myt2_d(1,iatm), acc(1) )
          accFoo2 = atomicadd( myt2_d(2,iatm), acc(2) )
          accFoo3 = atomicadd( myt2_d(3,iatm), acc(3) )

          ! naive atomic
          tempFoo = atomicadd( forceAtoms_d(7,iatm,flip), ttemp(1) )
          tempFoo = atomicadd( forceAtoms_d(8,iatm,flip), ttemp(2) )
          tempFoo = atomicadd( forceAtoms_d(9,iatm,flip), ttemp(3) )
        endif
      endif
    
    end subroutine make_fluid

    attributes(device) pure function cross(a,b)
      implicit none
      real :: cross(3)
      real, dimension(3), intent(in) :: a, b
    
      cross(1) = a(2) * b(3) - a(3) * b(2)
      cross(2) = a(3) * b(1) - a(1) * b(3)
      cross(3) = a(1) * b(2) - a(2) * b(1)
     end function cross

    attributes(device) pure function xcross(a,b)
      implicit none
      real :: xcross
      real, dimension(3), intent(in) :: a, b

      xcross = a(2) * b(3) - a(3) * b(2)
    end function xcross

    attributes(device) pure function ycross(a,b)
      implicit none
      real :: ycross
      real, dimension(3), intent(in) :: a, b

      ycross = a(3) * b(1) - a(1) * b(3)
    end function ycross

    attributes(device) pure function zcross(a,b)
      implicit none
      real :: zcross
      real, dimension(3), intent(in) :: a, b

      zcross = a(1) * b(2) - a(2) * b(1)
    end function zcross

    attributes(device) subroutine compute_onebelt_density_twofluids(i,j,k, flip, Rsum, Bsum, Dsum)
      implicit none
      integer, intent(in) :: i,j,k, flip
      real, intent(out) :: Rsum,Bsum
      integer, intent(out) :: Dsum
      integer :: i1,j1,k1, l
    
      Rsum=ZERO; Bsum=ZERO; Dsum=0
      !compute mean density value
      do l = 1, linksd3q27
        i1 = i + exd3q27(l)
        j1 = j + eyd3q27(l)
        k1 = k + ezd3q27(l)
        
        ! CYCLE_OUT_INTERVAL(ishift, minx-nbuff, maxx+nbuff)
        ! CYCLE_OUT_INTERVAL(jshift, miny-nbuff, maxy+nbuff)
        ! CYCLE_OUT_INTERVAL(kshift, minz-nbuff, maxz+nbuff)
        if (i1<0 .or. i1>nx+1 .or. j1<0 .or. j1>ny+1 .or. k1<0 .or. k1>nz+1) then
          stop_d = __LINE__
          write(*,*) '1belt] RangeErr i,j,k=',linear(i,j,k),linear(i1,j1,k1)
        endif

        if ( myfluid_d(i1,j1,k1, flip)==fluid_fluid .and. &
              myfluid_d(i1,j1,k1, 3-flip)==fluid_fluid) then
            Rsum = Rsum + rhoR_d(i1,j1,k1)
            Bsum = Bsum + rhoB_d(i1,j1,k1)
            Dsum = Dsum + 1            
        endif
      enddo

      if (Dsum /= ZERO) then
        Rsum = Rsum/Dsum
        Bsum = Bsum/Dsum
      endif
    end subroutine compute_onebelt_density_twofluids

    attributes(device) subroutine compute_secbelt_density_twofluids(i,j,k, flip, Rsum, Bsum, Dsum)
      implicit none
      integer, intent(in) :: i,j,k, flip
      real, intent(out) :: Rsum,Bsum
      integer, intent(out) :: Dsum
      integer :: i1,j1,k1, l

      Rsum=ZERO; Bsum=ZERO; Dsum=0
      !compute mean density value
      do l = 1, ndouble
        i1 = pimage(i + exdouble(l), nx)
        j1 = pimage(j + eydouble(l), ny)
        k1 = pimage(k + ezdouble(l), nz)

        if (i1<0 .or. i1>nx+1 .or. j1<0 .or. j1>ny+1 .or.k1<0 .or. k1>nz+1) then
          stop_d = __LINE__
          write(*,*) '2belt] RangeErr i,j,k=',linear(i,j,k),linear(i1,j1,k1)
        endif

        if ( myfluid_d(i1,j1,k1, flip)==fluid_fluid .and. &
              myfluid_d(i1,j1,k1, 3-flip)==fluid_fluid) then
            Rsum = Rsum + rhoR_d(i1,j1,k1)
            Bsum = Bsum + rhoB_d(i1,j1,k1)
            Dsum = Dsum + 1
            if (linear(i,j,k)== 270503 ) write(*,*) "2ndbelt", linear(i1,j1,k1), Dsum
        endif
      enddo

      if (Dsum /= ZERO) then
        Rsum = Rsum/Dsum
        Bsum = Bsum/Dsum
      endif
    end subroutine compute_secbelt_density_twofluids

    attributes(device) subroutine fix_onebelt_density_twofluids(i,j,k, flip, Rsum, Bsum, Dsum)
      implicit none
      integer, intent(in) :: i,j,k, flip
      real, intent(out) :: Rsum,Bsum
      integer, intent(out) :: Dsum

      ! TODO
      write (*,*) "fix_onebelt_density_twofluids]TO BE IMPL'ED"
    end subroutine fix_onebelt_density_twofluids

    attributes(device) subroutine compute_onebelt_density_2fl_weights(i,j,k, flip, Rsum, Bsum, Dsum)
      implicit none
      integer, intent(in) :: i,j,k, flip
      real, intent(out) :: Rsum,Bsum
      real, intent(out) :: Dsum
      integer :: i1,j1,k1, l
    
      Rsum=ZERO; Bsum=ZERO; Dsum=ZERO
      !compute mean density value
      do l = 1, linksd3q27
        i1 = i + exd3q27(l)
        j1 = j + eyd3q27(l)
        k1 = k + ezd3q27(l)

        if (i1<0 .or. i1>nx+1 .or. j1<0 .or. j1>ny+1 .or. k1<0 .or. k1>nz+1) then
          stop_d = __LINE__
          write(*,*) '1beltWghts] RangeErr i,j,k=',linear(i,j,k),linear(i1,j1,k1)
        endif
        
        if ( myfluid_d(i1,j1,k1, flip)==fluid_fluid) then
            Rsum = Rsum + pd3q27(l) * rhoR_d(i1,j1,k1)
            Bsum = Bsum + pd3q27(l) * rhoB_d(i1,j1,k1)
            Dsum = Dsum + pd3q27(l)            
        endif
      enddo

      if (Dsum /= ZERO) then
        Rsum = Rsum/Dsum
        Bsum = Bsum/Dsum
      endif
    end subroutine compute_onebelt_density_2fl_weights
    
    attributes(device) subroutine initialize_newnode_fluids(i,j,k, Rsum,Bsum, u,v,w, flip)
      implicit none
      integer, intent(in) :: i,j,k,flip
      real, intent(in) :: Rsum,Bsum, u,v,w
      real :: eqR,eqB
      integer :: l

      rhoR_d(i,j,k) = Rsum
      rhoB_d(i,j,k) = Bsum

      do l = 0, npops-1
	      eqR = equil(Rsum, u,v,w, l)
        popsR_d(i,j,k, l, flip) = eqR

	      eqB = equil(Bsum, u,v,w, l)
        popsB_d(i,j,k, l, flip) = eqB
      end do
    end subroutine initialize_newnode_fluids
    

    attributes(global) subroutine time_step_CG(step, flip)
      integer, value :: step,flip
      real    :: rhoR,rhoB,rhosum,invrho, u,v,w
      real    :: rhoavg,omega,oneminusomega
      real    :: rhoR_shifted,rhoB_shifted,rhosum_shifted
!      real    :: grad_rhoRx,grad_rhoRy,grad_rhoRz, grad_rhoBx,grad_rhoBy,grad_rhoBz
      real    :: psi_shifted
      real    :: psix,psiy,psiz,psinorm_sq,psinorm,acoeff,e_dot_psi
      real    :: temp,temp1, cosphi, feq
      real, parameter :: mylimit=1.e-20
      real    :: tempR(0:npops-1), tempB(0:npops-1), fdum(0:npops-1)

      integer :: i,j,k, i1,j1,k1,l
      
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA time_step_CG] i,k:',i,j,k

      if (myfluid_d(i,j,k, flip) == fluid_fluid) then
        rhoR = rhoR_d(i,j,k)
        rhoB = rhoB_d(i,j,k)
        rhosum = rhoR + rhoB

        rhoavg = rhoR/rhosum*viscR + rhoB/rhosum*viscB
        omega  = 1.0 / ( rhoavg/cssq  + 0.5)
        oneminusomega = 1.0 - omega

        invrho = 1.0 / rhosum
        u   = compute_u_1fl(popsR_d,i,j,k,flip, invrho) + compute_u_1fl(popsB_d,i,j,k,flip, invrho)
        v   = compute_v_1fl(popsR_d,i,j,k,flip, invrho) + compute_v_1fl(popsB_d,i,j,k,flip, invrho)
        w   = compute_w_1fl(popsR_d,i,j,k,flip, invrho) + compute_w_1fl(popsB_d,i,j,k,flip, invrho)


        !bgk step
        do l = 0, npops-1
            tempR(l) = popsR_d(i,j,k, l, flip)*oneminusomega + equil(rhoR,u,v,w, l)*omega
            tempB(l) = popsB_d(i,j,k, l, flip)*oneminusomega + equil(rhoB,u,v,w, l)*omega
        enddo


        ! Psi calc
        psix = ZERO
        psiy = ZERO
        psiz = ZERO
        do l = 1, npops-1
          i1 = mod(i + ex_d(l) +nx-1, nx) + 1
          j1 = mod(j + ey_d(l) +ny-1, ny) + 1
          k1 = mod(k + ez_d(l) +nz-1, nz) + 1
          if (i1<1 .or. i1>nx .or. j1<1 .or. j1>ny .or.k1<1 .or. k1>nz) then
            stop_d = __LINE__
          endif

          rhoR_shifted = rhoR_d(i1,j1,k1)
          rhoB_shifted = rhoB_d(i1,j1,k1)
          psi_shifted = (rhoR_shifted - rhoB_shifted) / (rhoR_shifted + rhoB_shifted)
          psix = psix + a_d(l)*ex_d(l)* psi_shifted
          psiy = psiy + a_d(l)*ey_d(l)* psi_shifted
          psiz = psiz + a_d(l)*ez_d(l)* psi_shifted
        enddo

        !perturbation step
        psinorm_sq = psix*psix + psiy*psiy + psiz*psiz
        if (psinorm_sq>mylimit) then
          psinorm = sqrt(psinorm_sq)
          acoeff =  NINE/FOUR * omega * sigma_cg

          do l = 0, npops-1
              e_dot_psi = ex_d(l)*psix + ey_d(l)*psiy + ez_d(l)*psiz
              temp = psinorm*(p_d(l)*(e_dot_psi*e_dot_psi)/psinorm_sq - b_l_d(l))
              !if(isnan(temp)) temp=ZERO

              fdum(l) = tempR(l) + tempB(l) + acoeff*temp
          enddo

          !recoloring step
          do l=0, npops-1
              feq = equil(rhoR,ZERO,ZERO,ZERO, l) + equil(rhoB,ZERO,ZERO,ZERO, l)
              e_dot_psi = ex_d(l)*psix + ey_d(l)*psiy + ez_d(l)*psiz
              temp1 = rec_fact_d(l) * psinorm

              if (temp1<=mylimit) then
                cosphi = ZERO
              else
                cosphi = e_dot_psi/temp1
              endif

              temp = beta_CG * rhoR * rhoB * cosphi/(rhosum*rhosum)

              tempR(l) = fdum(l)*rhoR/rhosum + temp*feq
              tempB(l) = fdum(l)*rhoB/rhosum - temp*feq
          enddo
        endif

        ! Store to mem
        do l = 0, npops-1
          popsR_d(i,j,k,l, flip) = tempR(l)
          popsB_d(i,j,k,l, flip) = tempB(l)
        enddo
        
      endif
    end subroutine time_step_CG

    attributes(global) subroutine partBB_CG(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l
      integer :: atm_i,atm_j,atm_k, atm_st,atm_en, i_atm,i_list
      integer :: deltai,deltaj,deltak, atmDist,okDist
      integer :: xxx,yyy,zzz


      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA partBB_CG] i,k:',i,j,k

      atm_i = blockIdx%x 
      atm_j = blockIdx%y 
      atm_k = blockIdx%z 
      atm_st = findAtoms_d(1, atm_i,atm_j,atm_k)
      atm_en = findAtoms_d(2, atm_i,atm_j,atm_k)
      if (atm_st>0) then
        do i_list = atm_st, atm_en
          i_atm = listAtoms_d(i_list)
          
          if ( myfluid_d(i,j,k, flip)==fluid_spherelist ) then
            xxx = nint(pos_d(1,i_atm, flip))
            yyy = nint(pos_d(2,i_atm, flip))
            zzz = nint(pos_d(3,i_atm, flip))

            deltai = getDeltai(i, xxx)
            deltaj = getDeltaj(j, yyy)
            deltak = getDeltak(k, zzz)
            
            atmDist = deltai*deltai + deltaj*deltaj + deltak*deltak
            ! only call n2p_bb if near THIS atom...
            if (atmDist < sphereMax_d) then
              if (debugfluid_d(i,j,k) /= i_atm) then
                xxx = nint(pos_d(1,debugfluid_d(i,j,k), flip))
                yyy = nint(pos_d(2,debugfluid_d(i,j,k), flip))
                zzz = nint(pos_d(3,debugfluid_d(i,j,k), flip))

                deltai = getDeltai(i, xxx)
                deltaj = getDeltaj(j, yyy)
                deltak = getDeltak(k, zzz)
                okDist = deltai*deltai + deltaj*deltaj + deltak*deltak
                write(*,*) 'partBB_CG]Wrong atom', linear(i,j,k), &
                    i_atm+0.001*atmDist, debugfluid_d(i,j,k)+0.001*okDist
                ! stop_d = __LINE__
              endif
              call node_to_particle_bounce_back_bc2(step, flip, i_atm, i,j,k)
            endif
          endif
        enddo
      endif
    end subroutine partBB_CG

    attributes(global) subroutine partBB_CG_p2n(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l


      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA partBB_CG] i,k:',i,j,k

      if ( myfluid_d(i,j,k, flip)==fluid_spherelist ) then
        call particle_to_node_bounce_back_bc2_phase1(step, flip, i,j,k)
      endif
    end subroutine partBB_CG_p2n

    attributes(global) subroutine partBB_CG_p2n2(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l
      integer :: atm_i,atm_j,atm_k, atm_st,atm_en, i_atm,i_list
      integer :: deltai,deltaj,deltak, atmDist
      integer :: xxx,yyy,zzz


      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA partBB_CG] i,k:',i,j,k

      atm_i = blockIdx%x 
      atm_j = blockIdx%y 
      atm_k = blockIdx%z 
      atm_st = findAtoms_d(1, atm_i,atm_j,atm_k)
      atm_en = findAtoms_d(2, atm_i,atm_j,atm_k)
      if (atm_st>0) then
        do i_list = atm_st, atm_en
          i_atm = listAtoms_d(i_list)
          
          if ( myfluid_d(i,j,k, flip)==fluid_spherelist ) then
            xxx = nint(pos_d(1,i_atm, flip))
            yyy = nint(pos_d(2,i_atm, flip))
            zzz = nint(pos_d(3,i_atm, flip))

            deltai = getDeltai(i, xxx)
            deltaj = getDeltaj(j, yyy)
            deltak = getDeltak(k, zzz)

            atmDist = deltai*deltai + deltaj*deltaj + deltak*deltak
            ! only call P2N_bb if near THIS atom...
            if (atmDist < sphereMax_d) then
              if (debugfluid_d(i,j,k) /= i_atm) then
                write(*,*) 'partBB_CG_p2n2]Wrong atom', linear(i,j,k), &
                    i_atm+0.001*atmDist, debugfluid_d(i,j,k)+0
                ! stop_d = __LINE__
              endif
              call particle_to_node_bounce_back_bc2_phase2(step, flip, i_atm, i,j,k)
            endif
          endif
        enddo
      endif
    end subroutine partBB_CG_p2n2

    attributes(global) subroutine copy_CG(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l, flop

      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA copy_CG] i,k:',i,j,k

      flop = 3 - flip
      
      do l = 0, npops-1
        popsR_d(i,j,k,l, flop) = popsR_d(i,j,k,l, flip)
        popsB_d(i,j,k,l, flop) = popsB_d(i,j,k,l, flip)
      enddo
    end subroutine copy_CG


    attributes(global) subroutine stream_CG(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l, flop
      integer :: i1,j1,k1
      integer :: atm_i,atm_j,atm_k, atm_st,atm_en, i_atm,i_list


      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA stream_CG] i,k:',i,j,k

      flop = 3 - flip


      ! Streaming      
      if (myfluid_d(i,j,k, flip) /= fluid_spheredead) then
        do l = 1, npops-1
          i1 = i + ex_d(l)
          j1 = j + ey_d(l)
          k1 = k + ez_d(l)
          popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
          popsB_d(i1,j1,k1,l, flop) = popsB_d(i,j,k,l, flip)
        enddo
      endif
    
    end subroutine stream_CG

    attributes(global) subroutine stream_CG_x(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l, flop
      integer :: i1,j1,k1
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1

      if (j>ny+1) return
      if (k>nz+1) return

      flop = 3 - flip

      ! Stream x=0
      i = 0
      if (myfluid_d(i,j,k, flip) /= fluid_spheredead) then
        do l = 1, npops-1
          i1 = i + ex_d(l)
          j1 = j + ey_d(l)
          k1 = k + ez_d(l)
          if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz) then
            popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
            popsB_d(i1,j1,k1,l, flop) = popsB_d(i,j,k,l, flip)
          endif
        end do
      endif

      ! Stream x=nx+1
      i = nx+1
      if (myfluid_d(i,j,k, flip) /= fluid_spheredead) then
        do l = 1, npops-1
          i1 = i + ex_d(l)
          j1 = j + ey_d(l)
          k1 = k + ez_d(l)
          if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz) then
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
    if (k>nz+1) return

    flop = 3 - flip

    ! Stream y=0
    j = 0
    if (myfluid_d(i,j,k, flip) /= fluid_spheredead) then
      do l = 1, npops-1
        i1 = i + ex_d(l)
        j1 = j + ey_d(l)
        k1 = k + ez_d(l)
        if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz) then
          popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
          popsB_d(i1,j1,k1,l, flop) = popsB_d(i,j,k,l, flip)
        endif
      end do
    endif

    ! Stream y=ny+1
    j = ny+1
    if (myfluid_d(i,j,k, flip) /= fluid_spheredead) then
      do l = 1, npops-1
        i1 = i + ex_d(l)
        j1 = j + ey_d(l)
        k1 = k + ez_d(l)
        if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz) then
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
    if (myfluid_d(i,j,k, flip) /= fluid_spheredead) then
      do l = 1, npops-1
        i1 = i + ex_d(l)
        j1 = j + ey_d(l)
        k1 = k + ez_d(l)
        if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz) then
          popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
          popsB_d(i1,j1,k1,l, flop) = popsB_d(i,j,k,l, flip)
        endif
      end do
    endif

    ! Stream z=nz+1
    k = nz+1
    if (myfluid_d(i,j,k, flip) /= fluid_spheredead) then
      do l = 1, npops-1
        i1 = i + ex_d(l)
        j1 = j + ey_d(l)
        k1 = k + ez_d(l)
        if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz) then
          popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
          popsB_d(i1,j1,k1,l, flop) = popsB_d(i,j,k,l, flip)
        endif
      end do
    endif
  end subroutine stream_CG_z

    

    ! PARTICLES 	-----------------------------------------------
    ! pos(1:3) -> xxx,yyy,zzz
    ! pos(4:6) -> vxx,vyy,vzz
    ! pos(7:9) -> oxx,oyy,ozz
    ! pos(10:13) -> q0,q1,q2,q3
    ! force(1:3) -> fxx,fyy,fzz
    ! force(4:6) -> fxb,fyb,fzb
    ! force(7:9) -> tqx,tqy,tqz
    ! force(10:12) -> tqb,tqb,tqb

    attributes(global) subroutine time_step_force_particles(step, flip)
      integer, value    :: step,flip
      integer           :: i,j
      integer :: flop
  
      i = (blockIdx%x-1) * TILE_DIMPART + threadIdx%x
      if (i>numAtoms) return

      if (1==step .and. 1==i) write(*,*) 'CUDA time_step_force_particles] i:',i

      if (withSidewall) call add_sidewall(step, flip, i)

      if (withHz) then
        do j = i+1, numAtoms	! loop over succ parts
          call addHertzian(step, flip, i,j)
        enddo
      endif

      call addForce_particle_bounce_back(step, flip, i)
      call nve_lf(step, flip, i)

      ! call zero_force_part(step, flip, i)
      flop = 3 - flip

      forceAtoms_d(1,i,flop) = ext_fxx_d
      forceAtoms_d(2,i,flop) = ext_fyy_d
      forceAtoms_d(3,i,flop) = ext_fzz_d
      forceAtoms_d(4,i,flop) = 0.0
      forceAtoms_d(5,i,flop) = 0.0
      forceAtoms_d(6,i,flop) = 0.0

      myf_d(1,i) = 0.0
      myf_d(2,i) = 0.0
      myf_d(3,i) = 0.0
      myf2_d(1,i) = ext_fxx_d
      myf2_d(2,i) = ext_fyy_d
      myf2_d(3,i) = ext_fzz_d

      if (lrotate) then
        forceAtoms_d(7,i,flop) = ext_tqx_d
        forceAtoms_d(8,i,flop) = ext_tqy_d
        forceAtoms_d(9,i,flop) = ext_tqz_d
        forceAtoms_d(10,i,flop) = 0.0
        forceAtoms_d(11,i,flop) = 0.0
        forceAtoms_d(12,i,flop) = 0.0

        myt_d(1,i) = 0.0
        myt_d(2,i) = 0.0
        myt_d(3,i) = 0.0
        myt2_d(1,i) = ext_tqx_d
        myt2_d(2,i) = ext_tqy_d
        myt2_d(3,i) = ext_tqz_d
      endif

      if (i == 1) then
        if (countmk_d /= countrm_d) then
          write(*,*) 'countmk_d /= countrm_d ->', countmk_d, countrm_d
          stop_d = __LINE__
        endif

        if (oldpartVol_d>0 .and. partVol_d /= oldpartVol_d) then        
          write(*,*) 'partVol_d /= oldpartVol_d ->', partVol_d, oldpartVol_d
          stop_d = __LINE__      
        endif
        
        countmk_d = 0
        countrm_d = 0
        oldpartVol_d = partVol_d
        partVol_d = 0
        countn2p_d = 0
      endif
    end subroutine time_step_force_particles


    attributes(device) subroutine add_sidewall(step, flip, iatm)
      integer, value    :: step, flip, iatm
      real, parameter	  :: inflimit = sidewall_rdist + HALF, rlimit = ONE+HALF, gmin = FIVE*HALF * sidewall_k
      real, parameter	  :: suplimitx = nx - inflimit
      real, parameter	  :: suplimity = ny - inflimit
      real, parameter	  :: suplimitz = nz - inflimit
      real, parameter	  :: suprcapx = nx - 1.5
      real, parameter	  :: suprcapy = ny - 1.5
      real, parameter	  :: suprcapz = nz - 1.5
      real              :: rrr,ggg

      ! vmin = sidewall_k*(ONE**(FIVE*HALF))    
      if (1==step.and.1==iatm) write(*,*) 'sw]', inflimit, rlimit, gmin

      if(.not. xperiodic) then
        if (pos_d(1,iatm,flip)-rdim <= inflimit) then
          rrr = pos_d(1,iatm,flip)-rdim
          if (rrr<=rlimit) then
            ggg = gmin
          else
            ggg = FIVE*HALF*sidewall_k*(inflimit-rrr)**(THREE*HALF)
          endif
          
          forceAtoms_d(1,iatm,flip) = forceAtoms_d(1,iatm,flip) + ggg
          write(*,*) 'sw] xi:', iatm, ggg
        endif
  
        if( pos_d(1,iatm,flip)+rdim >= suplimitx) then
          rrr = pos_d(1,iatm,flip) + rdim
          if (rrr>=suprcapx) then
            ggg = gmin
          else
            ggg = FIVE*HALF * sidewall_k * (rrr-suplimitx)**(THREE*HALF)
          endif
          
          forceAtoms_d(1,iatm,flip) = forceAtoms_d(1,iatm,flip) - ggg
          write(*,*) 'sw] xi:', iatm, ggg
        endif
      endif

      if(.not. yperiodic) then
        if (pos_d(2,iatm,flip)-rdim <= inflimit) then
          rrr = pos_d(2,iatm,flip)-rdim
          if (rrr<=rlimit) then
            ggg = gmin
          else
            ggg = FIVE*HALF*sidewall_k*(inflimit-rrr)**(THREE*HALF)
          endif
          
          forceAtoms_d(2,iatm,flip) = forceAtoms_d(2,iatm,flip) + ggg
          write(*,*) 'sw] yi:', iatm, ggg
        endif
  
        if( pos_d(2,iatm,flip)+rdim >= suplimity) then
          rrr = pos_d(2,iatm,flip) + rdim
          if (rrr>=suprcapy) then
            ggg = gmin
          else
            ggg = FIVE*HALF * sidewall_k * (rrr-suplimity)**(THREE*HALF)
          endif
          
          forceAtoms_d(2,iatm,flip) = forceAtoms_d(2,iatm,flip) - ggg
          write(*,*) 'sw] yi:', iatm, ggg
        endif
      endif

      if(.not. zperiodic) then
        if (pos_d(3,iatm,flip)-rdim <= inflimit) then
          rrr = pos_d(3,iatm,flip)-rdim
          if (rrr<=rlimit) then
            ggg = gmin
          else
            ggg = FIVE*HALF*sidewall_k*(inflimit-rrr)**(THREE*HALF)
          endif
          
          forceAtoms_d(3,iatm,flip) = forceAtoms_d(3,iatm,flip) + ggg
          write(*,*) 'sw] zi:', iatm, ggg
        endif
  
        if( pos_d(3,iatm,flip)+rdim >= suplimitz) then
          rrr = pos_d(3,iatm,flip) + rdim
          if (rrr>=suprcapz) then
            ggg = gmin
          else
            ggg = FIVE*HALF * sidewall_k * (rrr-suplimitz)**(THREE*HALF)
          endif
          
          forceAtoms_d(3,iatm,flip) = forceAtoms_d(3,iatm,flip) - ggg
          write(*,*) 'sw] zi:', iatm, ggg
        endif
      endif
    end subroutine add_sidewall


    attributes(device) subroutine addForce_particle_bounce_back(step, flip, i)
      integer, value    :: step, flip, i
      integer :: flop
      real :: fxx,fyy,fzz
      real :: fxb,fyb,fzb, fxbo,fybo,fzbo
      real :: tqx,tqy,tqz
      real :: txb,tyb,tzb, txbo,tybo,tzbo

      ! if (1==step) write(*,*) 'CUDA addForce_particle_bounce_back] i:',i

      flop = 3 - flip

      !Double prec accum
      forceAtoms_d(1,i,flip) = myf2_d(1,i)
      forceAtoms_d(2,i,flip) = myf2_d(2,i)
      forceAtoms_d(3,i,flip) = myf2_d(3,i)

      fxx = forceAtoms_d(1,i,flip)
      fyy = forceAtoms_d(2,i,flip)
      fzz = forceAtoms_d(3,i,flip)
#ifdef DEBUG_FORCE
      debugline_d(i, 8,1) = pos_d(1,i,flip)
      debugline_d(i, 8,2) = pos_d(2,i,flip)
      debugline_d(i, 8,3) = pos_d(3,i,flip)

      debugline_d(i, 8,4) = partVol_d
      debugline_d(i, 9,1) = countmk_d
      debugline_d(i, 9,2) = countrm_d

      debugline_d(i, 1,1) = fxx
      debugline_d(i, 1,2) = fyy
      debugline_d(i, 1,3) = fzz
#endif      

      !Double prec accum
      forceAtoms_d(4,i,flip) = myf_d(1,i)
      forceAtoms_d(5,i,flip) = myf_d(2,i)
      forceAtoms_d(6,i,flip) = myf_d(3,i)
      
      fxb = forceAtoms_d(4,i,flip)
      fyb = forceAtoms_d(5,i,flip)
      fzb = forceAtoms_d(6,i,flip)
      
      fxbo = forceAtoms_d(4,i,flop)
      fybo = forceAtoms_d(5,i,flop)
      fzbo = forceAtoms_d(6,i,flop)

! f(t) = (f(t+1/2)+f(t-1/2))/2
      fxx = fxx + (fxb+fxbo) * HALF
      fyy = fyy + (fyb+fybo) * HALF
      fzz = fzz + (fzb+fzbo) * HALF
      
      forceAtoms_d(1,i,flip) = fxx
      forceAtoms_d(2,i,flip) = fyy
      forceAtoms_d(3,i,flip) = fzz

#ifdef DEBUG_FORCE
      debugline_d(i, 2,1) = fxx
      debugline_d(i, 2,2) = fyy
      debugline_d(i, 2,3) = fzz

      debugline_d(i, 3,1) = fxb
      debugline_d(i, 3,2) = fyb
      debugline_d(i, 3,3) = fzb

      debugline_d(i, 4,1) = fxbo
      debugline_d(i, 4,2) = fybo
      debugline_d(i, 4,3) = fzbo
#endif      
      
      if (.not. lrotate) return

      ! Get in Double prec 
      forceAtoms_d(7,i,flip) = myt2_d(1,i)
      forceAtoms_d(8,i,flip) = myt2_d(2,i)
      forceAtoms_d(9,i,flip) = myt2_d(3,i)

      forceAtoms_d(10,i,flip) = myt_d(1,i)
      forceAtoms_d(11,i,flip) = myt_d(2,i)
      forceAtoms_d(12,i,flip) = myt_d(3,i)

      tqx = forceAtoms_d(7,i,flip)
      tqy = forceAtoms_d(8,i,flip)
      tqz = forceAtoms_d(9,i,flip)

      txb = forceAtoms_d(10,i,flip)
      tyb = forceAtoms_d(11,i,flip)
      tzb = forceAtoms_d(12,i,flip)

      txbo = forceAtoms_d(10,i,flop)
      tybo = forceAtoms_d(11,i,flop)
      tzbo = forceAtoms_d(12,i,flop)

! f(t) = (f(t+1/2)+f(t-1/2))/2
      tqx = tqx + (txb+txbo) * HALF
      tqy = tqy + (tyb+tybo) * HALF
      tqz = tqz + (tzb+tzbo) * HALF

      forceAtoms_d(7,i,flip) = tqx
      forceAtoms_d(8,i,flip) = tqy
      forceAtoms_d(9,i,flip) = tqz

#ifdef DEBUG_FORCE
      debugline_d(i, 5,1) = tqx
      debugline_d(i, 5,2) = tqy
      debugline_d(i, 5,3) = tqz

      debugline_d(i, 6,1) = txb
      debugline_d(i, 6,2) = tyb
      debugline_d(i, 6,3) = tzb

      debugline_d(i, 7,1) = txbo
      debugline_d(i, 7,2) = tybo
      debugline_d(i, 7,3) = tzbo
#endif
    end subroutine addForce_particle_bounce_back


    attributes(device) subroutine addHertzian(step, flip, iatm,jatm)
      integer, value    :: step, flip, iatm,jatm
      real		:: xi(3), xj(3), rsq, mxrsqcut, rparcap, xdf,ydf,zdf, rrr,ggg
      real    :: ux,uy,uz, dotuv, visc, ftemp(3), tempFoo1

      real, parameter	:: rsqcut = hz_rmin**TWO
      real, parameter	:: gmin = FIVE*HALF * hz_k * (hz_rmin-hz_rcap)**(THREE*HALF)      
      real, parameter	:: rpar = 2*rdim
      real, parameter	:: lubfactor = lubric_k * SIX*Pi*  (rdim*rdim)**TWO / (rpar*rpar)
      real, parameter	:: rparcut = rpar + lubric_rmin
      real, parameter	:: rsrparcut = rparcut**TWO

      mxrsqcut = max(rsrparcut,rsqcut)
      rparcap = rpar + lubric_rcap

      if (1==step) write(*,*) 'Hz]', mxrsqcut, rparcap, gmin

      xi(1) = pos_d(1, iatm, flip)
      xi(2) = pos_d(2, iatm, flip)
      xi(3) = pos_d(3, iatm, flip)
      
      xj(1) = pos_d(1, jatm, flip)
      xj(2) = pos_d(2, jatm, flip)
      xj(3) = pos_d(3, jatm, flip)

      xdf = xj(1) - xi(1)
      ydf = xj(2) - xi(2)
      zdf = xj(3) - xi(3)

      ! write(*,*) 'Hz] xdf:', xdf,ydf,zdf
      if (xperiodic) xdf = xdf - nx*nint(xdf / nx)
      if (yperiodic) ydf = ydf - ny*nint(ydf / ny)
      if (zperiodic) zdf = zdf - nz*nint(zdf / nz)

      rsq = xdf*xdf + ydf*ydf + zdf*zdf
      if (rsq <= mxrsqcut) then
      	
        rrr = sqrt(rsq)
        if (rrr <= hz_rmin) then
          if (rrr<=hz_rcap) then
            ggg = gmin
          else
            ggg = FIVE*HALF * hz_k * (hz_rmin-rrr)**(THREE*HALF)
          endif
          
          ftemp(1) = ggg*xdf/rrr
          ftemp(2) = ggg*ydf/rrr
          ftemp(3) = ggg*zdf/rrr

          ! Atomic
          tempFoo1 = atomicadd( forceAtoms_d(1,iatm,flip), -ftemp(1) )
          tempFoo1 = atomicadd( forceAtoms_d(1,jatm,flip),  ftemp(1) )

          tempFoo1 = atomicadd( forceAtoms_d(2,iatm,flip), -ftemp(2) )
          tempFoo1 = atomicadd( forceAtoms_d(2,jatm,flip),  ftemp(2) )

          tempFoo1 = atomicadd( forceAtoms_d(3,iatm,flip), -ftemp(3) )
          tempFoo1 = atomicadd( forceAtoms_d(3,jatm,flip),  ftemp(3) )

          ! forceAtoms_d(1,iatm,flip) = forceAtoms_d(1,iatm,flip) - ggg*xdf/rrr
          ! forceAtoms_d(1,jatm,flip) = forceAtoms_d(1,jatm,flip) + ggg*xdf/rrr          

          ! forceAtoms_d(2,iatm,flip) = forceAtoms_d(2,iatm,flip) - ggg*ydf/rrr
          ! forceAtoms_d(2,jatm,flip) = forceAtoms_d(2,jatm,flip) + ggg*ydf/rrr

          ! forceAtoms_d(3,iatm,flip) = forceAtoms_d(3,iatm,flip) - ggg*zdf/rrr
          ! forceAtoms_d(3,jatm,flip) = forceAtoms_d(3,jatm,flip) + ggg*zdf/rrr
          write(*,*) 'hz] i:', iatm, ggg*xdf/rrr, ggg*ydf/rrr, ggg*zdf/rrr
        endif

        if (llubrication) then
          if (rrr<=rparcut) then
            ux = xdf/rrr
            uy = ydf/rrr
            uz = zdf/rrr
            
            if (rrr<=rparcap) rrr=rparcap
            
            if (1==step) write(*,*) '******************************  CUDA verify visco in HZ  ******************************'
            visc = 1
            ! if (.not. lunique_omega) then
              ! xcm = nint( xxx(iatm) + xdf*HALF )
              ! ycm = nint( yyy(iatm) + ydf*HALF )
              ! zcm = nint( zzz(iatm) + zdf*HALF )
              ! xcm = pimage(ixpbc,xcm,nx)
              ! ycm = pimage(iypbc,ycm,ny)
              ! zcm = pimage(izpbc,zcm,nz)
              ! visc=omega_to_viscosity(omega(xcm,ycm,zcm))
            ! endif
            
            dotuv = ux*(pos_d(4, iatm, flip)-pos_d(4, jatm, flip)) + uy*(pos_d(5, iatm, flip)-pos_d(5, jatm, flip)) &
                + uz*(pos_d(6, iatm, flip)-pos_d(6, jatm, flip))

            forceAtoms_d(1,iatm,flip) = forceAtoms_d(1,iatm,flip) - visc*lubfactor*ux* dotuv*(ONE/(rrr-rpar)-ONE)
            forceAtoms_d(1,jatm,flip) = forceAtoms_d(1,jatm,flip) + visc*lubfactor*ux* dotuv*(ONE/(rrr-rpar)-ONE)
                                                        
            forceAtoms_d(2,iatm,flip) = forceAtoms_d(2,iatm,flip) - visc*lubfactor*uy* dotuv*(ONE/(rrr-rpar)-ONE)
            forceAtoms_d(2,jatm,flip) = forceAtoms_d(2,jatm,flip) + visc*lubfactor*uy* dotuv*(ONE/(rrr-rpar)-ONE)

            forceAtoms_d(3,iatm,flip) = forceAtoms_d(3,iatm,flip) - visc*lubfactor*uz* dotuv*(ONE/(rrr-rpar)-ONE)
            forceAtoms_d(3,jatm,flip) = forceAtoms_d(3,jatm,flip) + visc*lubfactor*uz* dotuv*(ONE/(rrr-rpar)-ONE)
          endif
        endif
      endif
    end subroutine addHertzian


    attributes(device) subroutine nve_lf(step, flip, i)
      integer, value    :: step,flip, i
      integer :: flop
      real :: fxx,fyy,fzz
      real :: bxx,byy,bzz, vxx,vyy,vzz
      real :: xxx,yyy,zzz, xxo,yyo,zzo
      logical :: lfix_moment = .false.

      ! if (1==step) write(*,*) 'CUDA nve_lf] i:', i

      flop = 3 - flip

      xxo = pos_d(1,i,flip)
      yyo = pos_d(2,i,flip)
      zzo = pos_d(3,i,flip)
      ! write(*,*) 'x] i:', i, xxo,yyo,zzo
      vxx = pos_d(4,i,flip)
      vyy = pos_d(5,i,flip)
      vzz = pos_d(6,i,flip)
      ! write(*,*) 'v] i:', i, vxx,vyy,vzz

      fxx = forceAtoms_d(1,i,flip)
      fyy = forceAtoms_d(2,i,flip)
      fzz = forceAtoms_d(3,i,flip)
      if (abs(fxx)*rmass > 5.0E-1 .or. abs(fyy)*rmass > 5.0E-1 .or. abs(fzz)*rmass > 5.0E-1) then
        write(*,*) 'CUDA nve_lf] i:', i
        write(*,*) 'Big force:',fxx,fyy,fzz
        stop_d = __LINE__
       endif

!     update velocities
      bxx = vxx + rmass * fxx
      byy = vyy + rmass * fyy
      bzz = vzz + rmass * fzz
      
      if(lfix_moment)then
      endif

!     update positions
      xxx = xxo + bxx
      yyy = yyo + byy
      zzz = zzo + bzz

! calculate full timestep velocity
      vxx = HALF * (vxx + bxx)
      vyy = HALF * (vyy + byy)
      vzz = HALF * (vzz + bzz)

! calculate kinetic energy at time step n
      call getkin

! restore free atom half step velocity
      vxx = bxx
      vyy = byy
      vzz = bzz

! periodic boundary condition
      if (xperiodic) then
        xxx = pbc_images_centeredx(xxx, cx)

        ! xxx = xxx - cx
        ! aaa = ONE/real(nx)
        ! write(*,*) 'nvea]xxx, nint(..) ', xxx, nint(aaa*xxx)
        ! xxx = xxx - real(nx)*nint(aaa*xxx) + cx
        ! write(*,*) 'nvea]xxx, step ', xxx, step

        if (xxx< 0.5 .or. xxx>nx+0.5) then
          write(*,*) 'nve]x step,i,pos ', step, i, xxx
          stop_d = __LINE__
        endif
      endif
      if (yperiodic) then
        yyy = pbc_images_centeredy(yyy, cy)
        if (yyy< 0.5 .or. yyy>ny+0.5) then
          write(*,*) 'nve]y step,i,pos ', step, i, yyy
          stop_d = __LINE__
        endif
      endif
      if (zperiodic) then
        zzz = pbc_images_centeredz(zzz, cz)
        if (zzz< 0.5 .or. zzz>nz+0.5) then
          write(*,*) 'nve]z step,i,pos', step, i, zzz
          stop_d = __LINE__
        endif
      endif

      if (lrotate) then
	      call rot_midstep_impl_lf(step, flip, i)
      endif

! Write to mem
      pos_d(1,i,flop) = xxx
      pos_d(2,i,flop) = yyy
      pos_d(3,i,flop) = zzz
      ! write(*,*) 'xn] i:', i, xxx,yyy,zzz
      pos_d(4,i,flop) = vxx
      pos_d(5,i,flop) = vyy
      pos_d(6,i,flop) = vzz
      ! write(*,*) 'vn] i:', i, bxx,byy,bzz
    end subroutine nve_lf


    attributes(device) subroutine getkin
    end subroutine getkin


    attributes(device) subroutine rot_midstep_impl_lf(step, flip, i)
      integer, value    :: step, flip, i
      integer :: flop, itq
      integer, parameter :: mxquat = 20
! working variables
      real, parameter :: onefive = 1.5d0, sqquattol = 1.0d-24
      real :: delx,dely,delz
      real :: oqx,oqy,oqz
      real :: opx,opy,opz
      real :: oxx,oyy,ozz
      real :: tqx,tqy,tqz
      real :: txb,tyb,tzb, txbo,tybo,tzbo
      real :: trx,try,trz
      real :: eps,rnorm,myrot(9)
      real :: qn0,qn1,qn2,qn3
      real :: qn0a,qn1a,qn2a,qn3a
      real :: dq0,dq1,dq2,dq3
      real :: q0,q1,q2,q3      
      real :: engrotsum


      engrotsum = 0

      flop = 3 - flip
      
  
      oxx = pos_d(7,i,flip)
      oyy = pos_d(8,i,flip)
      ozz = pos_d(9,i,flip)

      q0 = pos_d(10,i,flip)
      q1 = pos_d(11,i,flip)
      q2 = pos_d(12,i,flip)
      q3 = pos_d(13,i,flip)
      
      tqx = forceAtoms_d(7,i,flip)
      tqy = forceAtoms_d(8,i,flip)
      tqz = forceAtoms_d(9,i,flip)


      !current rotational matrix 
      myrot = quat_2_rotmat(q0,q1,q2,q3)
      
!     store current angular velocity which are still at n-1/2
      opx = oxx
      opy = oyy
      opz = ozz

#ifdef DEBUG_ROT      
      debugline_d(i, 15,1) = tqx
      debugline_d(i, 15,2) = tqy
      debugline_d(i, 15,3) = tqz

      debugline_d(i, 16,1) = q0
      debugline_d(i, 16,2) = q1
      debugline_d(i, 16,3) = q2
      debugline_d(i, 16,4) = q3

      debugline_d(i, 17,1) = oxx
      debugline_d(i, 17,2) = oyy
      debugline_d(i, 17,3) = ozz

      debugline_d(i, 18,1) = myrot(1)
      debugline_d(i, 18,2) = myrot(2)
      debugline_d(i, 18,3) = myrot(3)
      
      debugline_d(i, 19,1) = myrot(4)
      debugline_d(i, 19,2) = myrot(5)
      debugline_d(i, 19,3) = myrot(6)

      debugline_d(i, 20,1) = myrot(7)
      debugline_d(i, 20,2) = myrot(8)
      debugline_d(i, 20,3) = myrot(9)

      debugline_d(i, 21,1) = rotinx
      debugline_d(i, 21,2) = rotiny
      debugline_d(i, 21,3) = rotinz
#endif      

!     compute angular velocity increment of one step
      trx = (tqx*myrot(1) + tqy*myrot(4) + tqz*myrot(7)) /rotinx + (rotiny-rotinz)*opy*opz/rotinx
      try = (tqx*myrot(2) + tqy*myrot(5) + tqz*myrot(8)) /rotiny + (rotinz-rotinx)*opz*opx/rotiny
      trz = (tqx*myrot(3) + tqy*myrot(6) + tqz*myrot(9)) /rotinz + (rotinx-rotiny)*opx*opy/rotinz
      
#ifdef DEBUG_ROT      
      debugline_d(i, 11,1) = trx
      debugline_d(i, 11,2) = try
      debugline_d(i, 11,3) = trz
#endif      

      delx = trx
      dely = try
      delz = trz
      
!     angular velocity at time step n
      opx = oxx + delx*HALF
      opy = oyy + dely*HALF
      opz = ozz + delz*HALF
      
      
!     add its contribution to the total rotational kinetic energy at time step n
      engrotsum = engrotsum + (rotinx*opx*opx + rotiny*opy*opy + rotinz*opz*opz)
       
!     compute derivative at n
      dq0 = (-q1*opx-q2*opy-q3*opz)*HALF
      dq1 = ( q0*opx-q3*opy+q2*opz)*HALF
      dq2 = ( q3*opx+q0*opy-q1*opz)*HALF
      dq3 = (-q2*opx+q1*opy+q0*opz)*HALF

#ifdef DEBUG_ROT
      debugline_d(i, 12,1) = dq0
      debugline_d(i, 12,2) = dq1
      debugline_d(i, 12,3) = dq2
      debugline_d(i, 12,4) = dq3
#endif
      
!     update at t+1/2
      qn0 = q0 + dq0*HALF
      qn1 = q1 + dq1*HALF
      qn2 = q2 + dq2*HALF
      qn3 = q3 + dq3*HALF
      
      rnorm=ONE/sqrt(qn0*qn0 + qn1*qn1 + qn2*qn2 + qn3*qn3)
      qn0 = qn0 * rnorm
      qn1 = qn1 * rnorm
      qn2 = qn2 * rnorm
      qn3 = qn3 * rnorm
      
!     angular velocity at timestep n+1/2
      oxx = oxx + delx
      oyy = oyy + dely
      ozz = ozz + delz
      
!     store angular velocity at timestep n+1/2
      pos_d(7,i,flop) = oxx
      pos_d(8,i,flop) = oyy
      pos_d(9,i,flop) = ozz

      opx = oxx
      opy = oyy
      opz = ozz
      
!     assign new quaternions
      
!     iteration of new quaternions (lab fixed)
  
      itq = 0
      eps = 1.0d9
  
      do while((itq<mxquat).and.(eps>sqquattol))
        itq=itq+1
        
!       compute derivative at n+1/2
        dq0 = (-qn1*opx-qn2*opy-qn3*opz)*HALF
        dq1 = ( qn0*opx-qn3*opy+qn2*opz)*HALF
        dq2 = ( qn3*opx+qn0*opy-qn1*opz)*HALF
        dq3 = (-qn2*opx+qn1*opy+qn0*opz)*HALF
        
!       update at t+1/2
        qn0a = qn0 + dq0*HALF
        qn1a = qn1 + dq1*HALF
        qn2a = qn2 + dq2*HALF
        qn3a = qn3 + dq3*HALF
        
        rnorm = ONE / sqrt(qn0a*qn0a + qn1a*qn1a + qn2a*qn2a + qn3a*qn3a)
        qn0a = qn0a * rnorm
        qn1a = qn1a * rnorm
        qn2a = qn2a * rnorm
        qn3a = qn3a*rnorm
        
!      convergence test 
        eps = (qn0a-qn0)*(qn0a-qn0)+(qn1a-qn1)*(qn1a-qn1)+(qn2a-qn2)*(qn2a-qn2)+ (qn3a-qn3)*(qn3a-qn3)
              
        qn0 = qn0a
        qn1 = qn1a
        qn2 = qn2a
        qn3 = qn3a
      enddo
      
! store new quaternions
      q0 = q0 + dq0
      q1 = q1 + dq1
      q2 = q2 + dq2
      q3 = q3 + dq3
      
      rnorm = ONE / sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3)
      q0 = q0 * rnorm
      q1 = q1 * rnorm
      q2 = q2 * rnorm
      q3 = q3 * rnorm

! Write to mem
      pos_d(10,i,flop) = q0
      pos_d(11,i,flop) = q1
      pos_d(12,i,flop) = q2
      pos_d(13,i,flop) = q3

#ifdef DEBUG_ROT
      debugline_d(i, 13,1) = q0
      debugline_d(i, 13,2) = q1
      debugline_d(i, 13,3) = q2
      debugline_d(i, 13,4) = q3

      debugline_d(i, 14,1) = itq
      debugline_d(i, 14,2) = eps
#endif      
    end subroutine rot_midstep_impl_lf


    attributes(device) pure function quat_2_rotmat(q0,q1,q2,q3)
     implicit none
     real, intent(in) :: q0,q1,q2,q3
     real, dimension(9) :: quat_2_rotmat
         
     quat_2_rotmat(1) = q0*q0 + q1*q1 - q2*q2 - q3*q3
     quat_2_rotmat(2) = TWO*(q1*q2 - q0*q3)
     quat_2_rotmat(3) = TWO*(q1*q3 + q0*q2)
     quat_2_rotmat(4) = TWO*(q1*q2 + q0*q3)
     quat_2_rotmat(5) = q0*q0 - q1*q1 + q2*q2 - q3*q3
     quat_2_rotmat(6) = TWO*(q2*q3-q0*q1)
     quat_2_rotmat(7) = TWO*(q1*q3-q0*q2)
     quat_2_rotmat(8) = TWO*(q2*q3+q0*q1)
     quat_2_rotmat(9) = q0*q0 - q1*q1 - q2*q2 + q3*q3
    end function quat_2_rotmat

    attributes(device) pure function qconj(qs)
      implicit none      
      real, intent(in), dimension(0:3) :: qs      
      real, dimension(0:3) :: qconj
      
      qconj(0) = qs(0)
      qconj(1:3) = -qs(1:3)
     end function qconj

    attributes(device) pure function qmult(qsa,qsb)
       implicit none
       real, intent(in), dimension(0:3) :: qsa,qsb
       real, dimension(0:3) :: qmult
       
       qmult(0) = qsa(0)*qsb(0) - qsa(1)*qsb(1) - qsa(2)*qsb(2) - qsa(3)*qsb(3)
       qmult(1:3) = qsa(0)*qsb(1:3) + qsb(0)*qsa(1:3) + cross(qsa(1:3),qsb(1:3))
    end function qmult
      
    attributes(device) pure function qtrimult(qsa,qsb,qsc)
       implicit none
       real, intent(in), dimension(0:3) :: qsa,qsb,qsc
       real, dimension(0:3) :: qtrimult,qtemps
       
       qtemps = qmult(qsb,qsc)
       qtrimult = qmult(qsa,qtemps)
    end function qtrimult


    attributes(device) subroutine node_to_particle_bounce_back_bc2(step, flip, iatm, i,j,k)
      integer, value    :: step, flip
      integer, value    :: iatm, i,j,k
      real              :: rtemp(3), ftemp(3),ttemp(3), rversor(3), otemp(3), modr, f2p,f2pR,f2pB
      real              :: qtemp(0:3),qversor(0:3),tempconj(0:3), oat(0:3)
      real              :: vx,vy,vz, vxs,vys,vzs, tempFoo1,tempFoo2,tempFoo3
      real              :: aaa,bbb,ccc
      real(8)           :: acc(3), accFoo1,accFoo2,accFoo3
      integer           :: iloop, ii,jj,kk, io,jo,ko, indlow,indhig
      integer           :: xxx,yyy,zzz, deltai,deltaj,deltak, isOk
      integer    :: tempInt


      xxx = nint(pos_d(1,iatm, flip))
      yyy = nint(pos_d(2,iatm, flip))
      zzz = nint(pos_d(3,iatm, flip))

      vxs = pos_d(4,iatm,flip)
      vys = pos_d(5,iatm,flip)
      vzs = pos_d(6,iatm,flip)

      if (lrotate) then
        qtemp(0) = pos_d(10,iatm,flip)
        qtemp(1) = pos_d(11,iatm,flip)
        qtemp(2) = pos_d(12,iatm,flip)
        qtemp(3) = pos_d(13,iatm,flip)
        qversor(0) = ZERO
        qversor(1) = pos_d(7,iatm,flip)
        qversor(2) = pos_d(8,iatm,flip)
        qversor(3) = pos_d(9,iatm,flip)
        tempconj = qconj(qtemp)
        oat = qtrimult(qtemp,qversor,tempconj)
        otemp = oat(1:3)
        
        rversor(1) = i - pos_d(1,iatm,flip)
        rversor(2) = j - pos_d(2,iatm,flip)
        rversor(3) = k - pos_d(3,iatm,flip)

        ! PBC treatment
        aaa = ONE/real(nx)
        bbb = ONE/real(ny)
        ccc = ONE/real(nz)
        rversor(1) = rversor(1) - real(nx)*nint(aaa*rversor(1))
        rversor(2) = rversor(2) - real(ny)*nint(bbb*rversor(2))
        rversor(3) = rversor(3) - real(nz)*nint(ccc*rversor(3))

        modr = sqrt( rversor(1)*rversor(1) + rversor(2)*rversor(2) + rversor(3)*rversor(3) )
        rversor(1:3) = rdim / modr * rversor(1:3)
      endif

      do iloop = 1, 9
        indlow = iloop*2 - 1
        indhig = iloop*2
      
        ii = pimage(i + ex(indlow), nx)
        jj = pimage(j + ey(indlow), ny)
        kk = pimage(k + ez(indlow), nz)

        io = pimage(i + ex(indhig), nx)
        jo = pimage(j + ey(indhig), ny)
        ko = pimage(k + ez(indhig), nz)
      
        if (myfluid_d(ii,jj,kk, flip)==fluid_fluid .and. myfluid_d(io,jo,ko, flip)/=fluid_fluid) then
          deltai = getDeltai(ii, xxx)
          deltaj = getDeltaj(jj, yyy)
          deltak = getDeltak(kk, zzz)
          if (abs(deltai)<=rmax_issub_d.and.abs(deltaj)<=rmax_issub_d.and.abs(deltak)<=rmax_issub_d) then
            isOk = 1
          else
            ! write(*,*) 'n2p deltaErr:', linear(ii,jj,kk), deltai,deltaj,deltak
            ! stop_d = __LINE__
          endif
          
          if (lrotate) then
            ! rtemp = rversor
            rtemp(1) = rversor(1) + HALF*ex_d(indlow)
            rtemp(2) = rversor(2) + HALF*ey_d(indlow)
            rtemp(3) = rversor(3) + HALF*ez_d(indlow)
            vx = vxs + xcross(otemp,rtemp)
            vy = vys + ycross(otemp,rtemp)
            vz = vzs + zcross(otemp,rtemp)
          endif
    
          f2pR = TWO * popsR_d(ii,jj,kk,indhig,flip) - &
            p_d(indhig)*buz_d*rhoR_d(ii,jj,kk)* (ex_d(indhig)*vx + ey_d(indhig)*vy + ez_d(indhig)*vz)
          f2pB = TWO * popsB_d(ii,jj,kk,indhig,flip) - &
            p_d(indhig)*buz_d*rhoB_d(ii,jj,kk)* (ex_d(indhig)*vx + ey_d(indhig)*vy + ez_d(indhig)*vz)
          f2p = f2pR + f2pB
    
          ftemp(1) = f2p * ex_d(indhig)
          ftemp(2) = f2p * ey_d(indhig)
          ftemp(3) = f2p * ez_d(indhig)
          
          ! Double
          acc(1) = ftemp(1)
          acc(2) = ftemp(2)
          acc(3) = ftemp(3)
          accFoo1 = atomicadd( myf_d(1,iatm), acc(1) )
          accFoo2 = atomicadd( myf_d(2,iatm), acc(2) )
          accFoo3 = atomicadd( myf_d(3,iatm), acc(3) )

          ! naive atomic
          tempFoo1 = atomicadd( forceAtoms_d(4,iatm,flip), ftemp(1) )
          tempFoo2 = atomicadd( forceAtoms_d(5,iatm,flip), ftemp(2) )
          tempFoo3 = atomicadd( forceAtoms_d(6,iatm,flip), ftemp(3) )
#ifdef DEBUG_N2P
          tempInt = 1 + atomicadd( countn2p_d, 1)
          debugn2pf_d(tempInt, 1) = linear(ii,jj,kk)
          debugn2pf_d(tempInt, 2) = indhig
          debugn2pf_d(tempInt, 3) = ftemp(1)
          debugn2pf_d(tempInt, 4) = ftemp(2)
          debugn2pf_d(tempInt, 5) = ftemp(3)
          ! write(*,*) 'f1 R+B n2px]', linear(ii,jj,kk), indhig, ftemp(1)
          ! write(*,*) 'f1 R+B n2py]', linear(ii,jj,kk), indhig, ftemp(2)
          ! write(*,*) 'f1 R+B n2pz]', linear(ii,jj,kk), indhig, ftemp(3)
#endif          
          
          if(lrotate)then
            ttemp(1) = xcross(rtemp,ftemp) * abs(ex_d(indhig))
            ttemp(2) = ycross(rtemp,ftemp) * abs(ey_d(indhig))
            ttemp(3) = zcross(rtemp,ftemp) * abs(ez_d(indhig))

            ! Double
            acc(1) = ttemp(1)
            acc(2) = ttemp(2)
            acc(3) = ttemp(3)
            accFoo1 = atomicadd( myt_d(1,iatm), acc(1) )
            accFoo2 = atomicadd( myt_d(2,iatm), acc(2) )
            accFoo3 = atomicadd( myt_d(3,iatm), acc(3) )

            ! naive atomic
            tempFoo1 = atomicadd( forceAtoms_d(10,iatm,flip), ttemp(1) )
            tempFoo2 = atomicadd( forceAtoms_d(11,iatm,flip), ttemp(2) )
            tempFoo3 = atomicadd( forceAtoms_d(12,iatm,flip), ttemp(3) )
#ifdef DEBUG_N2P
            debugn2pt_d(tempInt, 1) = linear(ii,jj,kk)
            debugn2pt_d(tempInt, 2) = indhig
            debugn2pt_d(tempInt, 3) = ttemp(1)
            debugn2pt_d(tempInt, 4) = ttemp(2)
            debugn2pt_d(tempInt, 5) = ttemp(3)
            ! write(*,*) 't1 R+B n2px]', linear(ii,jj,kk), indhig, ttemp(1)
            ! write(*,*) 't1 R+B n2py]', linear(ii,jj,kk), indhig, ttemp(2)
            ! write(*,*) 't1 R+B n2pz]', linear(ii,jj,kk), indhig, ttemp(3)
#endif            
          endif
        endif

        if (myfluid_d(ii,jj,kk, flip)/=fluid_fluid .and. myfluid_d(io,jo,ko, flip)==fluid_fluid) then
          deltai = getDeltai(io, xxx)
          deltaj = getDeltaj(jo, yyy)
          deltak = getDeltak(ko, zzz)
          if (abs(deltai)<=rmax_issub_d.and.abs(deltaj)<=rmax_issub_d.and.abs(deltak)<=rmax_issub_d) then
            isOk = 1
          else
            ! write(*,*) 'n2p deltaErr:', linear(ii,jj,kk), deltai,deltaj,deltak
            ! stop_d = __LINE__
          endif

          if (lrotate) then
            ! rtemp = rversor
            rtemp(1) = rversor(1) + HALF*ex_d(indhig)
            rtemp(2) = rversor(2) + HALF*ey_d(indhig)
            rtemp(3) = rversor(3) + HALF*ez_d(indhig)
            vx = vxs + xcross(otemp,rtemp)
            vy = vys + ycross(otemp,rtemp)
            vz = vzs + zcross(otemp,rtemp)
          endif
    
          f2pR = TWO * popsR_d(io,jo,ko,indlow,flip) - &
            p_d(indlow)*buz_d*rhoR_d(io,jo,ko) * (ex_d(indlow)*vx + ey_d(indlow)*vy + ez_d(indlow)*vz)
          f2pB = TWO * popsB_d(io,jo,ko,indlow,flip) - &
            p_d(indlow)*buz_d*rhoB_d(io,jo,ko) * (ex_d(indlow)*vx + ey_d(indlow)*vy + ez_d(indlow)*vz)
          f2p = f2pR + f2pB
    
          ftemp(1) = f2p * ex_d(indlow)
          ftemp(2) = f2p * ey_d(indlow)
          ftemp(3) = f2p * ez_d(indlow)

          ! Double
          acc(1) = ftemp(1)
          acc(2) = ftemp(2)
          acc(3) = ftemp(3)
          accFoo1 = atomicadd( myf_d(1,iatm), acc(1) )
          accFoo2 = atomicadd( myf_d(2,iatm), acc(2) )
          accFoo3 = atomicadd( myf_d(3,iatm), acc(3) )
          
          ! naive atomic
          tempFoo1 = atomicadd( forceAtoms_d(4,iatm,flip), ftemp(1) )
          tempFoo2 = atomicadd( forceAtoms_d(5,iatm,flip), ftemp(2) )
          tempFoo3 = atomicadd( forceAtoms_d(6,iatm,flip), ftemp(3) )
#ifdef DEBUG_N2P
          tempInt = 1 + atomicadd( countn2p_d, 1)
          debugn2pf_d(tempInt, 1) = linear(io,jo,ko)
          debugn2pf_d(tempInt, 2) = indlow
          debugn2pf_d(tempInt, 3) = ftemp(1)
          debugn2pf_d(tempInt, 4) = ftemp(2)
          debugn2pf_d(tempInt, 5) = ftemp(3)
          ! write(*,*) 'f2 R+B n2px]', linear(io,jo,ko), indlow, ftemp(1)
          ! write(*,*) 'f2 R+B n2py]', linear(io,jo,ko), indlow, ftemp(2)
          ! write(*,*) 'f2 R+B n2pz]', linear(io,jo,ko), indlow, ftemp(3)
#endif          
                    
          if(lrotate)then
            ttemp(1) = xcross(rtemp,ftemp) * abs(ex_d(indlow))
            ttemp(2) = ycross(rtemp,ftemp) * abs(ey_d(indlow))
            ttemp(3) = zcross(rtemp,ftemp) * abs(ez_d(indlow))

            ! Double
            acc(1) = ttemp(1)
            acc(2) = ttemp(2)
            acc(3) = ttemp(3)
            accFoo1 = atomicadd( myt_d(1,iatm), acc(1) )
            accFoo2 = atomicadd( myt_d(2,iatm), acc(2) )
            accFoo3 = atomicadd( myt_d(3,iatm), acc(3) )
            
            ! naive atomic
            tempFoo1 = atomicadd( forceAtoms_d(10,iatm,flip), ttemp(1) )
            tempFoo2 = atomicadd( forceAtoms_d(11,iatm,flip), ttemp(2) )
            tempFoo3 = atomicadd( forceAtoms_d(12,iatm,flip), ttemp(3) )
#ifdef DEBUG_N2P
            debugn2pt_d(tempInt, 1) = linear(io,jo,ko)
            debugn2pt_d(tempInt, 2) = indlow
            debugn2pt_d(tempInt, 3) = ttemp(1)
            debugn2pt_d(tempInt, 4) = ttemp(2)
            debugn2pt_d(tempInt, 5) = ttemp(3)
            ! write(*,*) 't2 R+B n2px]', linear(io,jo,ko), indlow, ttemp(1)
            ! write(*,*) 't2 R+B n2py]', linear(io,jo,ko), indlow, ttemp(2)
            ! write(*,*) 't2 R+B n2pz]', linear(io,jo,ko), indlow, ttemp(3)
#endif
          endif
        endif
      enddo
    end subroutine node_to_particle_bounce_back_bc2

    attributes(device) pure function pimage(i, nMax)
      integer, intent(in) :: i, nMax
      integer :: pimage

      pimage = mod(i +nMax-1, nMax) + 1
    end function


    attributes(device) pure function pbc_images_centeredx(x, atomx)
      real, intent(in) :: x, atomx
      real :: pbc_images_centeredx, xxs, aaa

      xxs = x - atomx
      aaa = ONE/real(nx)

      pbc_images_centeredx = xxs - real(nx)*nint(aaa*xxs) + atomx
    end function

    attributes(device) pure function pbc_images_centeredy(y, atomy)
      real, intent(in) :: y, atomy
      real :: pbc_images_centeredy, yys, bbb

      yys = y - atomy
      bbb = ONE/real(ny)

      pbc_images_centeredy = yys - real(ny)*nint(bbb*yys) + atomy
    end function

    attributes(device) pure function pbc_images_centeredz(z, atomz)
      real, intent(in) :: z, atomz
      real :: pbc_images_centeredz, zzs, ccc

      zzs = z - atomz
      ccc = ONE/real(nz)

      pbc_images_centeredz = zzs - real(nz)*nint(ccc*zzs) + atomz
    end function

    attributes(device) subroutine particle_to_node_bounce_back_bc2_phase1(step, flip, i,j,k)
      integer, value    :: step, flip
      integer, value    :: i,j,k
      integer           :: iloop, ii,jj,kk, io,jo,ko, indlow,indhig
      
      do iloop = 1, 9
        indlow = iloop*2 - 1
        indhig = iloop*2

        ii = pimage(i + ex(indlow), nx)
        jj = pimage(j + ey(indlow), ny)
        kk = pimage(k + ez(indlow), nz)

        io = pimage(i + ex(indhig), nx)
        jo = pimage(j + ey(indhig), ny)
        ko = pimage(k + ez(indhig), nz)

        if(myfluid_d(ii,jj,kk, flip)==fluid_fluid .and. myfluid_d(io,jo,ko, flip)/=fluid_fluid)then
#ifdef DEBUG_P2N
          write(*,*) 'p2n1]',linear(ii,jj,kk), indhig, linear(i,j,k), indlow
#endif
          popsR_d(i,j,k,indlow, flip) = popsR_d(ii,jj,kk,indhig, flip)
          popsB_d(i,j,k,indlow, flip) = popsB_d(ii,jj,kk,indhig, flip)
        endif

        if(myfluid_d(ii,jj,kk, flip)/=fluid_fluid .and. myfluid_d(io,jo,ko, flip)==fluid_fluid)then
#ifdef DEBUG_P2N
          write(*,*) 'p2n1]',linear(io,jo,ko), indlow, linear(i,j,k), indhig
#endif
          popsR_d(i,j,k,indhig, flip) = popsR_d(io,jo,ko,indlow, flip)
          popsB_d(i,j,k,indhig, flip) = popsB_d(io,jo,ko,indlow, flip)
        endif
      enddo
    end subroutine particle_to_node_bounce_back_bc2_phase1
    
    attributes(device) subroutine particle_to_node_bounce_back_bc2_phase2(step, flip, iatm, i,j,k)
      integer, value    :: step, flip
      integer, value    :: iatm, i,j,k
      real              :: rtemp(3), ftemp(3),ttemp(3), rversor(3), otemp(3), modr, f2pR,f2pB
      real              :: qtemp(0:3),qversor(0:3),tempconj(0:3), oat(0:3)
      real              :: vx,vy,vz, vxs,vys,vzs
      real              :: aaa,bbb,ccc
      integer           :: iloop, ii,jj,kk, io,jo,ko, indlow,indhig
      
      vxs = pos_d(4,iatm,flip)
      vys = pos_d(5,iatm,flip)
      vzs = pos_d(6,iatm,flip)

      vx = vxs
      vy = vys
      vz = vzs
      
      if (lrotate) then
        qtemp(0) = pos_d(10,iatm,flip)
        qtemp(1) = pos_d(11,iatm,flip)
        qtemp(2) = pos_d(12,iatm,flip)
        qtemp(3) = pos_d(13,iatm,flip)
        qversor(0) = ZERO
        qversor(1) = pos_d(7,iatm,flip)
        qversor(2) = pos_d(8,iatm,flip)
        qversor(3) = pos_d(9,iatm,flip)
        tempconj = qconj(qtemp)
        oat = qtrimult(qtemp,qversor,tempconj)
        otemp = oat(1:3)
        
        rversor(1) = i - pos_d(1,iatm,flip)
        rversor(2) = j - pos_d(2,iatm,flip)
        rversor(3) = k - pos_d(3,iatm,flip)
        ! PBC treatment
        aaa = ONE/real(nx)
        bbb = ONE/real(ny)
        ccc = ONE/real(nz)
        rversor(1) = rversor(1) - real(nx)*nint(aaa*rversor(1))
        rversor(2) = rversor(2) - real(ny)*nint(bbb*rversor(2))
        rversor(3) = rversor(3) - real(nz)*nint(ccc*rversor(3))
        
        modr = sqrt( rversor(1)*rversor(1) + rversor(2)*rversor(2) + rversor(3)*rversor(3) )
        rversor(1:3) = rdim / modr * rversor(1:3)
      endif

      do iloop = 1, 9
        indlow = iloop*2 - 1
        indhig = iloop*2
      
        ii = pimage(i + ex(indlow), nx)
        jj = pimage(j + ey(indlow), ny)
        kk = pimage(k + ez(indlow), nz)

        io = pimage(i + ex(indhig), nx)
        jo = pimage(j + ey(indhig), ny)
        ko = pimage(k + ez(indhig), nz)

        if (myfluid_d(ii,jj,kk, flip)==fluid_fluid .and. myfluid_d(io,jo,ko, flip)/=fluid_fluid)then
          if (lrotate) then
            ! rtemp = rversor
            rtemp(1) = rversor(1) + HALF*ex_d(indlow)
            rtemp(2) = rversor(2) + HALF*ey_d(indlow)
            rtemp(3) = rversor(3) + HALF*ez_d(indlow)
            vx = vxs + xcross(otemp,rtemp)
            vy = vys + ycross(otemp,rtemp)
            vz = vzs + zcross(otemp,rtemp)
          endif
          f2pR = p_d(indhig)*buz_d*rhoR_d(ii,jj,kk) * (ex_d(indhig)*vx + ey_d(indhig)*vy + ez_d(indhig)*vz)
          popsR_d(i,j,k,indlow, flip) = popsR_d(i,j,k,indlow, flip) - f2pR
          
          f2pB = p_d(indhig)*buz_d*rhoB_d(ii,jj,kk) * (ex_d(indhig)*vx + ey_d(indhig)*vy + ez_d(indhig)*vz)
          popsB_d(i,j,k,indlow, flip) = popsB_d(i,j,k,indlow, flip) - f2pB
#ifdef DEBUG_P2N
          write(*,*) 'R+B p2nx]', linear(ii,jj,kk), indhig, f2pR, f2pB
#endif
        endif

        if (myfluid_d(ii,jj,kk, flip)/=fluid_fluid .and. myfluid_d(io,jo,ko, flip)==fluid_fluid)then
          if (lrotate) then
            ! rtemp = rversor
            rtemp(1) = rversor(1) + HALF*ex_d(indhig)
            rtemp(2) = rversor(2) + HALF*ey_d(indhig)
            rtemp(3) = rversor(3) + HALF*ez_d(indhig)
            vx = vxs + xcross(otemp,rtemp)
            vy = vys + ycross(otemp,rtemp)
            vz = vzs + zcross(otemp,rtemp)
          endif
          f2pR = p_d(indlow)*buz_d*rhoR_d(io,jo,ko) * (ex_d(indlow)*vx + ey_d(indlow)*vy + ez_d(indlow)*vz)
          popsR_d(i,j,k,indhig, flip)= popsR_d(i,j,k,indhig, flip) - f2pR
          
          f2pB = p_d(indlow)*buz_d*rhoB_d(io,jo,ko) * (ex_d(indlow)*vx + ey_d(indlow)*vy + ez_d(indlow)*vz)
          popsB_d(i,j,k,indhig, flip)= popsB_d(i,j,k,indhig, flip) - f2pB
#ifdef DEBUG_P2N
          write(*,*) 'R+B p2nx]', linear(io,jo,ko), indlow, f2pR, f2pB
#endif
        endif
      enddo
    end subroutine particle_to_node_bounce_back_bc2_phase2
    #endif    

    ! Periodic BC 	-----------------------------------------------
    attributes(global) subroutine bc_per_x(step, flip)
      integer, value :: step,flip
      integer :: j,k !, l
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
      ! if (1==step .and. 1==j*k) write(*,*) 'CUDA bc_per_x] i,j:',j,k

      ! Apply BC also on virtual nodes...
      if (j<=ny+1 .and. k<=nz+1) then
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

      if (0<=i .and. i<=nx+1 .and. 0<=k .and. k<=nz+1) then
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
             if (ez_d(l)>0) popsR_d(i,j,1, l, flip) = popsR_d(i,j,nz+1, l, flip)
             if (ez_d(l)<0) popsR_d(i,j,nz, l, flip) = popsR_d(i,j,0, l, flip)
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

      if (1<=j .and. j<=ny .and. 1<=k .and. k<=nz) then
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

      if (1<=i .and. i<=nx .and. 1<=k .and. k<=nz) then
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
            popsR_d(i,j,nz+1, l, flip) = popsR_d(i,j,1, l, flip)
            popsB_d(i,j,nz+1, l, flip) = popsB_d(i,j,1, l, flip)
          
            popsR_d(i,j,0, l, flip) = popsR_d(i,j,nz, l, flip)
            popsB_d(i,j,0, l, flip) = popsB_d(i,j,nz, l, flip)
          end do
        endif

        if (faiRho) then
          rhoR_d(i,j,nz+1) = rhoR_d(i,j,1)
          rhoB_d(i,j,nz+1) = rhoB_d(i,j,1)

          rhoR_d(i,j,0) = rhoR_d(i,j,nz)
          rhoB_d(i,j,0) = rhoB_d(i,j,nz)
        endif
      endif
    end subroutine bc_per_z2


    attributes(global) subroutine bc_edge_z2(step, flip, faiPops,faiRho)
      integer, value :: step,flip
      logical, value :: faiPops, faiRho
      integer :: k, l
  
      k = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1

      if (1<=k .and. k<=nz) then
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
            popsR_d(0,j,0, l, flip) = popsR_d(nx,j,nz, l, flip)
            popsB_d(0,j,0, l, flip) = popsB_d(nx,j,nz, l, flip)

            popsR_d(0,j,nz+1, l, flip) = popsR_d(nx,j,1, l, flip)
            popsB_d(0,j,nz+1, l, flip) = popsB_d(nx,j,1, l, flip)

            popsR_d(nx+1,j,0, l, flip) = popsR_d(1,j,nz, l, flip)
            popsB_d(nx+1,j,0, l, flip) = popsB_d(1,j,nz, l, flip)

            popsR_d(nx+1,j,nz+1, l, flip) = popsR_d(1,j,1, l, flip)
            popsB_d(nx+1,j,nz+1, l, flip) = popsB_d(1,j,1, l, flip)
          end do
        endif

        if (faiRho) then
          rhoR_d (0,j,0) = rhoR_d (nx,j,nz)
          rhoB_d (0,j,0) = rhoB_d (nx,j,nz)

          rhoR_d (0,j,nz+1) = rhoR_d (nx,j,1)
          rhoB_d (0,j,nz+1) = rhoB_d (nx,j,1)

          rhoR_d (nx+1,j,0) = rhoR_d (1,j,nz)
          rhoB_d (nx+1,j,0) = rhoB_d (1,j,nz)

          rhoR_d (nx+1,j,nz+1) = rhoR_d (1,j,1)
          rhoB_d (nx+1,j,nz+1) = rhoB_d (1,j,1)
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
            popsR_d(i,0,0, l, flip) = popsR_d(i,ny,nz, l, flip)
            popsB_d(i,0,0, l, flip) = popsB_d(i,ny,nz, l, flip)

            popsR_d(i,0,nz+1, l, flip) = popsR_d(i,ny,1, l, flip)
            popsB_d(i,0,nz+1, l, flip) = popsB_d(i,ny,1, l, flip)

            popsR_d(i,ny+1,0, l, flip) = popsR_d(i,1,nz, l, flip)
            popsB_d(i,ny+1,0, l, flip) = popsB_d(i,1,nz, l, flip)

            popsR_d(i,ny+1,nz+1, l, flip) = popsR_d(i,1,1, l, flip)
            popsB_d(i,ny+1,nz+1, l, flip) = popsB_d(i,1,1, l, flip)
          end do
        endif

        if (faiRho) then
          rhoR_d (i,0,0) = rhoR_d (i,ny,nz)
          rhoB_d (i,0,0) = rhoB_d (i,ny,nz)

          rhoR_d (i,0,nz+1) = rhoR_d (i,ny,1)
          rhoB_d (i,0,nz+1) = rhoB_d (i,ny,1)

          rhoR_d (i,ny+1,0) = rhoR_d (i,1,nz)
          rhoB_d (i,ny+1,0) = rhoB_d (i,1,nz)

          rhoR_d (i,ny+1,nz+1) = rhoR_d (i,1,1)
          rhoB_d (i,ny+1,nz+1) = rhoB_d (i,1,1)
        endif
      endif
    end subroutine bc_edge_x2


    attributes(global) subroutine bc_corners(step, flip)
      integer, value :: step,flip
      integer :: j,k, l
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
    end subroutine bc_corners

    !!! Periodic treat myfluid
    attributes(global) subroutine isfluid_per_x2(step, flip)
      integer, value :: step,flip
      integer :: j,k, l
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1

      if (1<=j .and. j<=ny .and. 1<=k .and. k<=nz) then
        myfluid_d(nx+1, j,k, flip) = myfluid_d( 1,j,k, flip)
        myfluid_d(   0, j,k, flip) = myfluid_d(nx,j,k, flip)
      endif
    end subroutine isfluid_per_x2


    attributes(global) subroutine isfluid_per_y2(step, flip)
      integer, value :: step,flip
      integer :: i,k, l
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1

      if (1<=i .and. i<=nx .and. 1<=k .and. k<=nz) then
        myfluid_d(i,ny+1,k, flip) = myfluid_d(i, 1,k, flip)
        myfluid_d(i,   0,k, flip) = myfluid_d(i,ny,k, flip)
      endif
    end subroutine isfluid_per_y2


    attributes(global) subroutine isfluid_per_z2(step, flip)
      integer, value :: step,flip
      integer :: i,j, l
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      j = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1

      if (1<=j .and. j<=ny .and. 1<=i .and. i<=nx) then
        myfluid_d(i,j,nz+1, flip) = myfluid_d(i,j, 1, flip)        
        myfluid_d(i,j,   0, flip) = myfluid_d(i,j,nz, flip)
      endif
    end subroutine isfluid_per_z2


    attributes(global) subroutine isfluid_edge_z2(step, flip)
      integer, value :: step,flip
      integer :: k, l
  
      k = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1

      if (1<=k .and. k<=nz) then
        myfluid_d (   0,   0,k, flip) = myfluid_d (nx,ny,k, flip)
        myfluid_d (   0,ny+1,k, flip) = myfluid_d (nx, 1,k, flip)
        myfluid_d (nx+1,   0,k, flip) = myfluid_d ( 1,ny,k, flip)
        myfluid_d (nx+1,ny+1,k, flip) = myfluid_d ( 1, 1,k, flip)
      endif
    end subroutine isfluid_edge_z2

    attributes(global) subroutine isfluid_edge_y2(step, flip)
      integer, value :: step,flip
      integer :: j, l
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1

      if (1<=j .and. j<=ny) then
        myfluid_d (   0,j,   0, flip) = myfluid_d (nx,j, nz, flip)
        myfluid_d (   0,j,nz+1, flip) = myfluid_d (nx,j,  1, flip)
        myfluid_d (nx+1,j,   0, flip) = myfluid_d ( 1,j, nz, flip)
        myfluid_d (nx+1,j,nz+1, flip) = myfluid_d ( 1,j,  1, flip)
      endif
    end subroutine isfluid_edge_y2

    attributes(global) subroutine isfluid_edge_x2(step, flip)
      integer, value :: step,flip
      integer :: i, l
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1

      if (1<=i .and. i<=nx) then
        myfluid_d (i,   0,   0, flip) = myfluid_d (i,ny,nz, flip)
        myfluid_d (i,   0,nz+1, flip) = myfluid_d (i,ny, 1, flip)
        myfluid_d (i,ny+1,   0, flip) = myfluid_d (i, 1,nz, flip)
        myfluid_d (i,ny+1,nz+1, flip) = myfluid_d (i, 1, 1, flip)
      endif
    end subroutine isfluid_edge_x2


    attributes(global) subroutine bc_periodic_ext_x(step, flip)
      integer, value :: step,flip
      integer :: j,k, l
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
      ! if (1==step .and. 1==j*k) write(*,*) 'CUDA bc_periodic_ext_x] i,j:',j,k

      ! Copy on virtual nodes...
      if (0<=j .and. j<=ny+1 .and. 0<=k .and. k<=nz+1) then
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

      if (0<=i .and. i<=nx+1 .and. 0<=k .and. k<=nz+1) then
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
          popsR_d(i,j,0   , l, flip) = popsR_d(i,j,nz, l, flip)
          popsR_d(i,j,nz+1, l, flip) = popsR_d(i,j,1 , l, flip)
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

      if (0<=j .and. j<=ny+1 .and. 0<=k .and. k<=nz+1) then
	      rho1 = compute_rho(popsR_d,1,j,k,flip)

        popsR_d(1,j,k, 1, flip) = popsR_d(0, j,   k,    2, flip)
        popsR_d(1,j,k, 7, flip) = popsR_d(0, j-1, k,    8, flip)
        popsR_d(1,j,k,10, flip) = popsR_d(0, j+1, k,    9, flip)
        popsR_d(1,j,k,11, flip) = popsR_d(0, j  , k-1, 12, flip)
        popsR_d(1,j,k,14, flip) = popsR_d(0, j  , k+1, 13, flip)

	      rho2 = compute_rho(popsR_d,nx,j,k,flip)

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

      if (0<=i .and. i<=nx+1 .and. 0<=k .and. k<=nz+1) then
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
        rho1 = compute_rho(popsR_d,i,j,1,flip)

        popsR_d(i,j,1, 5, flip) = popsR_d(i,   j,   0,  6, flip)
        popsR_d(i,j,1,11, flip) = popsR_d(i-1, j,   0, 12, flip) + p_d(12)*buz_d*rho1 * wall_z0_d(1)
        popsR_d(i,j,1,13, flip) = popsR_d(i+1, j,   0, 14, flip) - p_d(14)*buz_d*rho1 * wall_z0_d(1)
        popsR_d(i,j,1,15, flip) = popsR_d(i,   j-1, 0, 16, flip)
        popsR_d(i,j,1,17, flip) = popsR_d(i,   j+1, 0, 18, flip)

        rho2 = compute_rho(popsR_d,i,j,nz,flip)

        popsR_d(i,j,nz, 6, flip) = popsR_d(i,   j,   nz+1,  5, flip)
        popsR_d(i,j,nz,12, flip) = popsR_d(i+1, j,   nz+1, 11, flip) - p_d(11)*buz_d*rho2 * wall_z1_d(1)
        popsR_d(i,j,nz,14, flip) = popsR_d(i-1, j,   nz+1, 13, flip) + p_d(13)*buz_d*rho2 * wall_z1_d(1)
        popsR_d(i,j,nz,16, flip) = popsR_d(i,   j+1, nz+1, 15, flip)
        popsR_d(i,j,nz,18, flip) = popsR_d(i,   j-1, nz+1, 17, flip)
      endif
    end subroutine bc_bb_z

    ! Bounce-Back BC 2 fluids
    attributes(global) subroutine bc_bb_x2(step, flip)
      integer, value :: step,flip
      integer :: j,k
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y
      ! if (1==step .and. 1==j*k) write(*,*) 'CUDA bc_bb_x] j,k:',j,k

      if (0<=j .and. j<=ny+1 .and. 0<=k .and. k<=nz+1) then	    
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

      if (0<=i .and. i<=nx+1 .and. 0<=k .and. k<=nz+1) then
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
        rho1 = compute_rho(popsR_d,i,j,1,flip)
        popsR_d(i,j,1, 5, flip) = popsR_d(i,   j,   0,  6, flip)
        popsR_d(i,j,1,11, flip) = popsR_d(i-1, j,   0, 12, flip) + p_d(12)*buz_d*rho1 * wall_z0_d(1)
        popsR_d(i,j,1,13, flip) = popsR_d(i+1, j,   0, 14, flip) - p_d(14)*buz_d*rho1 * wall_z0_d(1)
        popsR_d(i,j,1,15, flip) = popsR_d(i,   j-1, 0, 16, flip)
        popsR_d(i,j,1,17, flip) = popsR_d(i,   j+1, 0, 18, flip)

        rho2 = compute_rho(popsR_d,i,j,nz,flip)
        popsR_d(i,j,nz, 6, flip) = popsR_d(i,   j,   nz+1,  5, flip)
        popsR_d(i,j,nz,12, flip) = popsR_d(i+1, j,   nz+1, 11, flip) - p_d(11)*buz_d*rho2 * wall_z1_d(1)
        popsR_d(i,j,nz,14, flip) = popsR_d(i-1, j,   nz+1, 13, flip) + p_d(13)*buz_d*rho2 * wall_z1_d(1)
        popsR_d(i,j,nz,16, flip) = popsR_d(i,   j+1, nz+1, 15, flip)
        popsR_d(i,j,nz,18, flip) = popsR_d(i,   j-1, nz+1, 17, flip)

        rho3 = compute_rho(popsB_d,i,j,1,flip)
        popsB_d(i,j,1, 5, flip) = popsB_d(i,   j,   0,  6, flip)
        popsB_d(i,j,1,11, flip) = popsB_d(i-1, j,   0, 12, flip) + p_d(12)*buz_d*rho3 * wall_z0_d(1)
        popsB_d(i,j,1,13, flip) = popsB_d(i+1, j,   0, 14, flip) - p_d(14)*buz_d*rho3 * wall_z0_d(1)
        popsB_d(i,j,1,15, flip) = popsB_d(i,   j-1, 0, 16, flip)
        popsB_d(i,j,1,17, flip) = popsB_d(i,   j+1, 0, 18, flip)

        rho4 = compute_rho(popsB_d,i,j,nz,flip)
        popsB_d(i,j,nz, 6, flip) = popsB_d(i,   j,   nz+1,  5, flip)
        popsB_d(i,j,nz,12, flip) = popsB_d(i+1, j,   nz+1, 11, flip) - p_d(11)*buz_d*rho4 * wall_z1_d(1)
        popsB_d(i,j,nz,14, flip) = popsB_d(i-1, j,   nz+1, 13, flip) + p_d(13)*buz_d*rho4 * wall_z1_d(1)
        popsB_d(i,j,nz,16, flip) = popsB_d(i,   j+1, nz+1, 15, flip)
        popsB_d(i,j,nz,18, flip) = popsB_d(i,   j-1, nz+1, 17, flip)
      endif
    end subroutine bc_bb_z2
  end module kernels_fluid
