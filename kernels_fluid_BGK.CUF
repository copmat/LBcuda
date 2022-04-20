#include "defines.h"

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
#ifdef D3Q27
#define manual_u(pops)   invrho * ( pops(i,j,k,1,flip) - pops(i,j,k,2,flip) + pops(i,j,k,7,flip) - pops(i,j,k,8,flip) - pops(i,j,k,9,flip) + \
       pops(i,j,k,10,flip) + pops(i,j,k,11,flip) - pops(i,j,k,12,flip) - pops(i,j,k,13,flip) + pops(i,j,k,14,flip) + pops(i,j,k,19,flip) - pops(i,j,k,20,flip) + \
       pops(i,j,k,21,flip) - pops(i,j,k,22,flip) - pops(i,j,k,23,flip) + pops(i,j,k,24,flip) + pops(i,j,k,25,flip) - pops(i,j,k,26,flip) )

#define manual_v(pops)   invrho * ( pops(i,j,k,3,flip) - pops(i,j,k,4,flip) + pops(i,j,k,7,flip) - pops(i,j,k,8,flip) + pops(i,j,k,9,flip) - \
       pops(i,j,k,10,flip) + pops(i,j,k,15,flip) - pops(i,j,k,16,flip) - pops(i,j,k,17,flip) + pops(i,j,k,18,flip) + pops(i,j,k,19,flip) - pops(i,j,k,20,flip) + \
       pops(i,j,k,21,flip) - pops(i,j,k,22,flip) + pops(i,j,k,23,flip) - pops(i,j,k,24,flip) - pops(i,j,k,25,flip) + pops(i,j,k,26,flip) )
#define manual_w(pops)   invrho * ( pops(i,j,k,5,flip) - pops(i,j,k,6,flip) + pops(i,j,k,11,flip) - pops(i,j,k,12,flip) + pops(i,j,k,13,flip) - \
       pops(i,j,k,14,flip) + pops(i,j,k,15,flip) - pops(i,j,k,16,flip) + pops(i,j,k,17,flip) - pops(i,j,k,18,flip) + pops(i,j,k,19,flip) - pops(i,j,k,20,flip) - \
       pops(i,j,k,21,flip) + pops(i,j,k,22,flip) + pops(i,j,k,23,flip) - pops(i,j,k,24,flip) + pops(i,j,k,25,flip) - pops(i,j,k,26,flip) )
#else
#define manual_u(pops)   invrho * ( pops(i,j,k,1,flip) - pops(i,j,k,2,flip) + pops(i,j,k,7,flip) - pops(i,j,k,8,flip) - pops(i,j,k,9,flip) + \
       pops(i,j,k,10,flip) + pops(i,j,k,11,flip) - pops(i,j,k,12,flip) - pops(i,j,k,13,flip) + pops(i,j,k,14,flip) )

#define manual_v(pops)   invrho * ( pops(i,j,k,3,flip) - pops(i,j,k,4,flip) + pops(i,j,k,7,flip) - pops(i,j,k,8,flip) + pops(i,j,k,9,flip) - \
       pops(i,j,k,10,flip) + pops(i,j,k,15,flip) - pops(i,j,k,16,flip) - pops(i,j,k,17,flip) + pops(i,j,k,18,flip) )

#define manual_w(pops)   invrho * ( pops(i,j,k,5,flip) - pops(i,j,k,6,flip) + pops(i,j,k,11,flip) - pops(i,j,k,12,flip) + pops(i,j,k,13,flip) - \
       pops(i,j,k,14,flip) + pops(i,j,k,15,flip) - pops(i,j,k,16,flip) + pops(i,j,k,17,flip) - pops(i,j,k,18,flip) )
#endif
  module kernels_bgk
    use dimensions_m
    use kernels_fluid
    implicit none
    
!    private :: equil,compute_u_1fl,compute_v_1fl, &
!     compute_w_1fl
!    private :: compute_rho
  contains
   
    attributes(device) function equil(rho,u,v,w, l)
     real, intent(in) :: rho,u,v,w
     integer, intent(in) :: l
     real :: equil
     real :: uv

     uv = (1.0/cssq) * (u*ex_d(l) + v*ey_d(l) + w*ez_d(l))
     equil = rho * p_d(l)*(ONE+uv+HALF*(uv*uv)-(HALF/cssq) * (u*u + v*v + w*w))
    end function equil
 
    ! Setup
    attributes(global) subroutine setup(vx,vy,vz)
      real, value :: vx,vy,vz
      integer :: i,j,k, l,gli,glj,glk
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
      z = 3.14159265358979 / nz_d * k
      rho = densR_d !1 + 0.1*sin(x)*sin(x) * sin(y)*sin(y) * sin(z)*sin(z)
      
      gli = i + offset_d(1)
      glj = j + offset_d(2)
      glk = k + offset_d(3)
      
      
      u = vx
      v = vy
      w = vz
      
      if(store_vel_d)then
        vel_d(1,i,j,k) = u
        vel_d(2,i,j,k) = v
        vel_d(3,i,j,k) = w
      endif
      
      ! write(*,*) 'CUDA setup',rho,i,j,k
  
      rhoR_d(i,j,k) = rho
      do l = 0, npops-1
        eq = equil(rho, u,v,w, l)
        popsR_d(i,j,k, l, 1) = eq
      end do
#ifndef APPLYBC
      myfluid_d(i,j,k, 1) = fluid_fluid
      myfluid_d(i,j,k, 2) = fluid_fluid
#endif
    end subroutine setup
    
    attributes(global) subroutine setup_taylorgreen(vx,vy,vz)
      real, value :: vx,vy,vz
      integer :: i,j,k, l,gli,glj,glk
      real    :: rho, u,v,w, eq, x,y,z, xs,ys,zs
      real, parameter :: fact=1.0/sqrt(3.0)
  
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
      z = 3.14159265358979 / nz_d * k
      rho = densR_d !1 + 0.1*sin(x)*sin(x) * sin(y)*sin(y) * sin(z)*sin(z)
      
      gli = i + offset_d(1)
      glj = j + offset_d(2)
      glk = k + offset_d(3)
      
      xs = real(gli)/real(glx)*2.0*Pi
      ys = real(glj)/real(gly)*2.0*Pi
      zs = real(glk)/real(glz)*2.0*Pi
      
      
      u = -fact*sin(xs)*cos(ys)*cos(zs)*0.1
      v = -fact*cos(xs)*sin(ys)*cos(zs)*0.1
      w = 2.0*fact*cos(xs)*cos(ys)*sin(zs)*0.1
      
      if(store_vel_d)then
        vel_d(1,i,j,k) = u
        vel_d(2,i,j,k) = v
        vel_d(3,i,j,k) = w
      endif
      
      ! write(*,*) 'CUDA setup',rho,i,j,k
  
      rhoR_d(i,j,k) = rho
      do l = 0, npops-1
        eq = equil(rho, u,v,w, l)
        popsR_d(i,j,k, l, 1) = eq
      end do

      myfluid_d(i,j,k, 1) = fluid_fluid
      myfluid_d(i,j,k, 2) = fluid_fluid
    end subroutine setup_taylorgreen
    
    attributes(global) subroutine setup_p(vx,vy,vz)
      real, value :: vx,vy,vz
      integer :: i,j,k, l,gli,glj,glk
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
      z = 3.14159265358979 / nz_d * k
      rho = densR_d !1 + 0.1*sin(x)*sin(x) * sin(y)*sin(y) * sin(z)*sin(z)
      
      gli = i + offset_d(1)
      glj = j + offset_d(2)
      glk = k + offset_d(3)
      
      if(gli==glx/4 .and. glj==gly/2 .and. glk==glz/2)then
        rho = densR_d*(ONE+1.e-2)
      endif
      
      u = vx
      v = vy
      w = vz
      
      if(store_vel_d)then
        vel_d(1,i,j,k) = u
        vel_d(2,i,j,k) = v
        vel_d(3,i,j,k) = w
      endif
      
      ! write(*,*) 'CUDA setup',rho,i,j,k
  
      rhoR_d(i,j,k) = rho
      do l = 0, npops-1
        eq = equil(rho, u,v,w, l)
        popsR_d(i,j,k, l, 1) = eq
      end do

      myfluid_d(i,j,k, 1) = fluid_fluid
      myfluid_d(i,j,k, 2) = fluid_fluid
    end subroutine setup_p


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
    attributes(global) subroutine time_step_mom_BGK(step, flip, omega)
      integer, value :: step,flip
      real, value :: omega
      real    :: rho,invrho, u,v,w,oneminusomega,omegaminusone,ush,vsh,wsh
      integer :: i,j,k, l, flop
      integer :: i1,j1,k1
      
      oneminusomega = 1.0 -omega
      omegaminusone = omega - 1.0
      
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA time_step] i,k:',i,j,k
      
      if (myfluid_d(i,j,k, flip) == fluid_fluid) then
      
        rho = manual_rho(popsR_d) !compute_rho(popsR_d,i,j,k,flip)
        if (rho<minRho .or. rho>10) then
         write(*,*) 'Rho error:',rho,i,j,k
         stop_d = __LINE__
        endif
        rhoR_d(i,j,k) = rho

    

        invrho = 1.0 / rho
        u   = manual_u(popsR_d) !compute_u_1fl(popsR_d,i,j,k,flip, invrho)
        v   = manual_v(popsR_d) !compute_v_1fl(popsR_d,i,j,k,flip, invrho)
        w   = manual_w(popsR_d) !compute_w_1fl(popsR_d,i,j,k,flip, invrho)
        if(store_vel_d) then
          vel_d(1,i,j,k) = u
          vel_d(2,i,j,k) = v
          vel_d(3,i,j,k) = w
        endif
        ush = u + fx_d * invrho
        vsh = v + fy_d * invrho
        wsh = w + fz_d * invrho
  
        flop = 3 - flip
        do l = 0, npops-1
          !i1 = i + ex_d(l)
          !j1 = j + ey_d(l)
          !k1 = k + ez_d(l)
          !popsR_d(i1,j1,k1, l, flop) = popsR_d(i,j,k, l, flip)*oneminusomega + equil(rho, u,v,w, l)*omega
          popsR_d(i,j,k,l,flip) = popsR_d(i,j,k,l,flip)*oneminusomega + &
           equil(rho, u,v,w, l)*omegaminusone + equil(rho, ush,vsh,wsh, l)
        end do
      
      else
      
        rhoR_d(i,j,k) = MINDENS
      endif
      
    end subroutine time_step_mom_BGK

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

      rho = manual_rho(popsR_d) !compute_rho(popsR_d,i,j,k,flip)
      if (rho < 1.0E-7) then
       write(*,*) 'Rho error:',rho,i,j,k
       stop_d = __LINE__
      endif

      invrho = 1.0 / rho
      u   = manual_u(popsR_d) !compute_u_1fl(popsR_d,i,j,k,flip, invrho)
      v   = manual_v(popsR_d) !compute_v_1fl(popsR_d,i,j,k,flip, invrho)
      w   = manual_w(popsR_d) !compute_w_1fl(popsR_d,i,j,k,flip, invrho)

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

      rho = manual_rho(popsR_d) !compute_rho(popsR_d,i,j,k,flip)
      if (rho < 1.0E-7) then
       write(*,*) 'Rho error:',rho,i,j,k
       stop_d = __LINE__
      endif

      invrho = 1.0 / rho
      u   = manual_u(popsR_d) !compute_u_1fl(popsR_d,i,j,k,flip, invrho)
      v   = manual_v(popsR_d) !compute_v_1fl(popsR_d,i,j,k,flip, invrho)
      w   = manual_w(popsR_d) !compute_w_1fl(popsR_d,i,j,k,flip, invrho)

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
        ! THIS is not like Kruger...same as ML
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
      rhoR = manual_rho(popsR_d) !compute_rho(popsR_d,i,j,k,flip)

      if (rhoR<minRho .or. rhoR>10) then
        ! write(*,*) 'init_rho_isfluid_BGK]Range error rhoR', step, linear(i,j,k), rhoR
        stop_d = __LINE__
      endif

      rhoR_d(i,j,k) = rhoR
    else
      rhoR_d(i,j,k) = MINDENS
    endif
    end subroutine init_rho_isfluid_BGK
    
 



    attributes(global) subroutine time_step_BGK(step, flip, omega)
      integer, value :: step,flip
      real, value :: omega
      real    :: rhoR, invrho, u,v,w  ,oneminusomega,omegaminusone,ush,vsh,wsh
      integer :: i,j,k, l
      
      oneminusomega = 1.0 - omega
      omegaminusone = omega - 1.0
      
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA time_step_CG] i,k:',i,j,k

      if (myfluid_d(i,j,k, flip) == fluid_fluid) then
        rhoR = rhoR_d(i,j,k)

        invrho = 1.0 / rhoR
        u   = manual_u(popsR_d) !compute_u_1fl(popsR_d,i,j,k,flip, invrho)
        v   = manual_v(popsR_d) !compute_v_1fl(popsR_d,i,j,k,flip, invrho)
        w   = manual_w(popsR_d) !compute_w_1fl(popsR_d,i,j,k,flip, invrho)
        if(store_vel_d) then
          vel_d(1,i,j,k) = u
          vel_d(2,i,j,k) = v
          vel_d(3,i,j,k) = w
        endif
        ush = u + fx_d * invrho
        vsh = v + fy_d * invrho
        wsh = w + fz_d * invrho
        
        !bgk step
        do l = 0, npops-1
          popsR_d(i,j,k,l, flip) = popsR_d(i,j,k, l, flip)*oneminusomega + &
           equil(rhoR, u,v,w, l)*omegaminusone + equil(rhoR, ush,vsh,wsh, l)
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
    
    attributes(global) subroutine bc_BGK_per_x(step, flip, faiPops,faiRho)
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
            
            popsR_d(0,j,k, l, flip) = popsR_d(nx,j,k, l, flip)
          end do
        endif

        if (faiRho) then
          rhoR_d (nx+1,j,k) = rhoR_d (1,j,k)

          rhoR_d (0,j,k) = rhoR_d (nx,j,k)

        endif

      endif
      
    end subroutine bc_BGK_per_x
    
    attributes(global) subroutine bc_per_x_hvar(step,faiRho,faiRhoB,faiVel)
      integer, value :: step
      logical, value :: faiRho,faiRhoB,faiVel
      integer :: i,j,k, l,m
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y
      ! if (1==step .and. 1==j*k) write(*,*) 'CUDA bc_per_x2] i,j:',j,k
      
      if (j<=ny .and. k<=nz_d) then

        if (faiRho) then
          do i=1,nbuff
            rhoR_d (nx+i,j,k) = rhoR_d (i,j,k)

            rhoR_d (1-i,j,k) = rhoR_d (nx+1-i,j,k)
          enddo
        endif
        if (faiRhoB) then
          do i=1,nbuff
            rhoB_d (nx+i,j,k) = rhoB_d (i,j,k)

            rhoB_d (1-i,j,k) = rhoB_d (nx+1-i,j,k)
          enddo
        endif
        if (faiVel) then
          do i=1,nbuff
            do m=1,3
              vel_d (m,nx+i,j,k) = vel_d (m,i,j,k)

              vel_d (m,1-i,j,k) = vel_d (m,nx+1-i,j,k)
            enddo
          enddo
        endif

      endif
      
    end subroutine bc_per_x_hvar
    
    
    attributes(global) subroutine bc_per_x_backf(step,faiF)
      integer, value :: step
      logical, value :: faiF
      integer :: j,k, l,m
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
      ! if (1==step .and. 1==j*k) write(*,*) 'CUDA bc_per_x2] i,j:',j,k
      
      if (1<=j .and. j<=ny .and. 1<=k .and. k<=nz_d) then
        
        if (faiF) then
          do m=1,3
            force_d (m,1,j,k) = force_d (m,nx+1,j,k) + force_d (m,1,j,k)

            force_d (m,nx,j,k) = force_d (m,0,j,k) + force_d (m,nx,j,k)
          enddo
        endif

      endif
      
    end subroutine bc_per_x_backf

    attributes(global) subroutine bc_BGK_per_y(step, flip, faiPops,faiRho)
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
            
            popsR_d(i,0,k, l, flip) = popsR_d(i,ny,k, l, flip)
          end do
        endif

        if (faiRho) then
          rhoR_d(i,ny+1,k) = rhoR_d(i,1,k)

          rhoR_d(i,0,k) = rhoR_d(i,ny,k)
        endif
        
      endif

      
      
      
    end subroutine bc_BGK_per_y
    
    attributes(global) subroutine bc_per_y_hvar(step,faiRho,faiRhoB,faiVel)
      integer, value :: step
      logical, value :: faiRho,faiRhoB,faiVel
      integer :: i,j,k, l,m
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y
      ! if (1==step .and. 1==i*k) write(*,*) 'CUDA bc_per_y2] i,k:',i,k
      
      if (i<=nx .and. k<=nz_d) then

          

        if (faiRho) then
          do j=1,nbuff
            rhoR_d(i,ny+j,k) = rhoR_d(i,j,k)

            rhoR_d(i,1-j,k) = rhoR_d(i,ny+1-j,k)
          enddo
        endif
        if (faiRhoB) then
          do j=1,nbuff
            rhoB_d(i,ny+j,k) = rhoB_d(i,j,k)

            rhoB_d(i,1-j,k) = rhoB_d(i,ny+1-j,k)
          enddo
        endif
        if (faiVel) then
          do j=1,nbuff
            do m=1,3
              vel_d(m,i,ny+j,k) = vel_d(m,i,j,k)

              vel_d(m,i,1-j,k) = vel_d(m,i,ny+1-j,k)
            enddo
          enddo
        endif
        
      endif

      
      
      
    end subroutine bc_per_y_hvar
    
    attributes(global) subroutine bc_per_y_backf(step,faiF)
      integer, value :: step
      logical, value :: faiF
      integer :: i,k, l,m
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
      ! if (1==step .and. 1==i*k) write(*,*) 'CUDA bc_per_y2] i,k:',i,k
      
      if (1<=i .and. i<=nx .and. 1<=k .and. k<=nz_d) then
        
        if (faiF) then
          do m=1,3
            force_d(m,i,1,k) = force_d(m,i,ny+1,k) + force_d(m,i,1,k)

            force_d(m,i,ny,k) = force_d(m,i,0,k) + force_d(m,i,ny,k)
          enddo
        endif
        
      endif

      
      
      
    end subroutine bc_per_y_backf
    
    attributes(global) subroutine bc_BGK_per_z(step, flip, faiPops,faiRho)
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
          
            popsR_d(i,j,0, l, flip) = popsR_d(i,j,nz_d, l, flip)
          end do
        endif

        if (faiRho) then
          rhoR_d(i,j,nz_d+1) = rhoR_d(i,j,1)

          rhoR_d(i,j,0) = rhoR_d(i,j,nz_d)
        endif
      endif
    end subroutine bc_BGK_per_z
    
    attributes(global) subroutine bc_per_z_hvar(step,faiRho,faiRhoB,faiVel)
      integer, value :: step
      logical, value :: faiRho,faiRhoB,faiVel
      integer :: i,j,k, l,m
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIM + threadIdx%y
      ! if (1==step .and. 1==i*j) write(*,*) 'CUDA bc_per_z2] i,j:',i,j

      if (j<=ny .and. i<=nx) then

        if (faiRho) then
          do k=1,nbuff
            rhoR_d(i,j,nz_d+k) = rhoR_d(i,j,k)

            rhoR_d(i,j,1-k) = rhoR_d(i,j,nz_d+1-k)
          enddo
        endif
        if (faiRhoB) then
          do k=1,nbuff
            rhoB_d(i,j,nz_d+k) = rhoB_d(i,j,k)

            rhoB_d(i,j,1-k) = rhoB_d(i,j,nz_d+1-k)
          enddo
        endif
        if (faiVel) then
          do k=1,nbuff
            do m=1,3
              vel_d(m,i,j,nz_d+k) = vel_d(m,i,j,k)

              vel_d(m,i,j,1-k) = vel_d(m,i,j,nz_d+1-k)
            enddo
          enddo
        endif
      endif
    end subroutine bc_per_z_hvar
    
    attributes(global) subroutine bc_per_z_backf(step,faiF)
      integer, value :: step
      logical, value :: faiF
      integer :: i,j, l,m
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      j = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
      ! if (1==step .and. 1==i*j) write(*,*) 'CUDA bc_per_z2] i,j:',i,j

      if (1<=j .and. j<=ny .and. 1<=i .and. i<=nx) then
        
        if (faiF) then
          do m=1,3
            force_d(m,i,j,1) = force_d(m,i,j,nz_d+1) + force_d(m,i,j,1)

            force_d(m,i,j,nz_d) = force_d(m,i,j,0) + force_d(m,i,j,nz_d)
          enddo
        endif
      endif
    end subroutine bc_per_z_backf
    
    attributes(global) subroutine merge_per_z_backf(step,faiF)
      integer, value :: step
      logical, value :: faiF
      integer :: i,j, l,ii,m
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      j = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
      ! if (1==step .and. 1==i*j) write(*,*) 'CUDA bc_per_z2] i,j:',i,j

      if (0<=j .and. j<=ny+1 .and. 0<=i .and. i<=nx+1) then
        
        if (faiF) then
          do m=1,3
            force_d(m,i,j,1) = force_buf_up_d(m,i,j,1) + force_d(m,i,j,1)

            force_d(m,i,j,nz_d) = force_buf_down_d(m,i,j,1) + force_d(m,i,j,nz_d)
          enddo
        endif
        
      endif
    end subroutine merge_per_z_backf
    
    attributes(global) subroutine bc_BGK_edge_z(step, flip, faiPops,faiRho)
      integer, value :: step,flip
      logical, value :: faiPops, faiRho
      integer :: k, l
  
      k = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1

      if (0<=k .and. k<=nz_d+1) then
        if (faiPops) then
          do l = 1, npops-1
            popsR_d(0,0,k, l, flip) = popsR_d(nx,ny,k, l, flip)

            popsR_d(0,ny+1,k, l, flip) = popsR_d(nx,1,k, l, flip)

            popsR_d(nx+1,0,k, l, flip) = popsR_d(1,ny,k, l, flip)

            popsR_d(nx+1,ny+1,k, l, flip) = popsR_d(1,1,k, l, flip)
          end do
        endif

        if (faiRho) then
          rhoR_d (0,0,k) = rhoR_d (nx,ny,k)

          rhoR_d (0,ny+1,k) = rhoR_d (nx,1,k)

          rhoR_d (nx+1,0,k) = rhoR_d (1,ny,k)

          rhoR_d (nx+1,ny+1,k) = rhoR_d (1,1,k)
        endif
      endif
    end subroutine bc_BGK_edge_z
    
    attributes(global) subroutine bc_edge_z_hvar(step,faiRho,faiRhoB,faiVel)
      integer, value :: step
      logical, value :: faiRho,faiRhoB,faiVel
      integer :: k,i,j, l,m
  
      k = (blockIdx%x-1) * TILE_DIM + threadIdx%x  - nbuff

      if (k<=nz_d+nbuff) then

        if (faiRho) then
          do j=1,nbuff
            do i=1,nbuff
               rhoR_d(1-i,1-j,k) = rhoR_d(nx+1-i,ny+1-j,k)
               rhoR_d(1-i,ny+j,k) = rhoR_d(nx+1-i,j,k)
               rhoR_d(nx+i,1-j,k) = rhoR_d( i,ny+1-j,k)
               rhoR_d(nx+i,ny+j,k) = rhoR_d( i, j,k)
             enddo
           enddo
        endif
        if (faiRhoB) then
          do j=1,nbuff
            do i=1,nbuff
               rhoB_d(1-i,1-j,k) = rhoB_d(nx+1-i,ny+1-j,k)
               rhoB_d(1-i,ny+j,k) = rhoB_d(nx+1-i,j,k)
               rhoB_d(nx+i,1-j,k) = rhoB_d( i,ny+1-j,k)
               rhoB_d(nx+i,ny+j,k) = rhoB_d( i, j,k)
             enddo
           enddo
        endif
        if (faiVel) then
          do j=1,nbuff
            do i=1,nbuff
              do m=1,3
                vel_d(m,1-i,1-j,k) = vel_d(m,nx+1-i,ny+1-j,k)
                vel_d(m,1-i,ny+j,k) = vel_d(m,nx+1-i,j,k)
                vel_d(m,nx+i,1-j,k) = vel_d(m, i,ny+1-j,k)
                vel_d(m,nx+i,ny+j,k) = vel_d(m, i, j,k)
              enddo
            enddo
          enddo
        endif
      endif
    end subroutine bc_edge_z_hvar
    
    attributes(global) subroutine bc_edge_z_backf(step,faiF)
      integer, value :: step
      logical, value :: faiF
      integer :: k, l,m
  
      k = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1

      if (0<=k .and. k<=nz_d+1) then

        if (faiF) then
          do m=1,3
            force_d(m,nx,ny,k) = force_d(m,0,0,k) + force_d(m,nx,ny,k)

            force_d(m,nx,1,k) = force_d(m,0,ny+1,k) + force_d(m,nx,1,k)

            force_d(m,1,ny,k) = force_d(m,nx+1,0,k) + force_d(m,1,ny,k)

            force_d(m,1,1,k) = force_d(m,nx+1,ny+1,k) + force_d(m,1,1,k)
          enddo
        endif
      endif
    end subroutine bc_edge_z_backf
    
    attributes(global) subroutine bc_BGK_edge_y(step, flip, faiPops,faiRho)
      integer, value :: step,flip
      logical, value :: faiPops, faiRho
      integer :: j, l
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1

      if (1<=j .and. j<=ny) then
        if (faiPops) then
          do l = 1, npops-1
            popsR_d(0,j,0, l, flip) = popsR_d(nx,j,nz_d, l, flip)

            popsR_d(0,j,nz_d+1, l, flip) = popsR_d(nx,j,1, l, flip)

            popsR_d(nx+1,j,0, l, flip) = popsR_d(1,j,nz_d, l, flip)

            popsR_d(nx+1,j,nz_d+1, l, flip) = popsR_d(1,j,1, l, flip)
          end do
        endif

        if (faiRho) then
          rhoR_d (0,j,0) = rhoR_d (nx,j,nz_d)

          rhoR_d (0,j,nz_d+1) = rhoR_d (nx,j,1)

          rhoR_d (nx+1,j,0) = rhoR_d (1,j,nz_d)

          rhoR_d (nx+1,j,nz_d+1) = rhoR_d (1,j,1)
        endif
      endif
    end subroutine bc_BGK_edge_y
    
    attributes(global) subroutine bc_edge_y_hvar(step,faiRho,faiRhoB,faiVel)
      integer, value :: step
      logical, value :: faiRho,faiRhoB,faiVel
      integer :: j,i,k, l,m
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x

      if (j<=ny) then

        if (faiRho) then
          do k=1,nbuff
            do i=1,nbuff
              rhoR_d(1-i,j,1-k) = rhoR_d(nx+1-i,j, nz_d+1-k)
              rhoR_d(1-i,j,nz_d+k) = rhoR_d(nx+1-i,j,  k)
              rhoR_d(nx+i,j,1-k) = rhoR_d( i,j, nz_d+1-k)
              rhoR_d(nx+i,j,nz_d+k) = rhoR_d( i,j,  k)
            enddo
          enddo
        endif
        if (faiRhoB) then
          do k=1,nbuff
            do i=1,nbuff
              rhoB_d(1-i,j,1-k) = rhoB_d(nx+1-i,j, nz_d+1-k)
              rhoB_d(1-i,j,nz_d+k) = rhoB_d(nx+1-i,j,  k)
              rhoB_d(nx+i,j,1-k) = rhoB_d( i,j, nz_d+1-k)
              rhoB_d(nx+i,j,nz_d+k) = rhoB_d( i,j,  k)
            enddo
          enddo
        endif
        if (faiVel) then
          do k=1,nbuff
            do i=1,nbuff
              do m=1,3
                vel_d(m,1-i,j,1-k) = vel_d(m,nx+1-i,j, nz_d+1-k)
                vel_d(m,1-i,j,nz_d+k) = vel_d(m,nx+1-i,j,  k)
                vel_d(m,nx+i,j,1-k) = vel_d(m, i,j, nz_d+1-k)
                vel_d(m,nx+i,j,nz_d+k) = vel_d(m, i,j,  k)
              enddo
            enddo
          enddo
        endif
      endif
    end subroutine bc_edge_y_hvar
    
    attributes(global) subroutine bc_edge_y_backf(step,faiF)
      integer, value :: step
      logical, value :: faiF
      integer :: j, l,m
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1

      if (1<=j .and. j<=ny) then

        if (faiF) then
          do m=1,3
            force_d(m,nx,j,nz_d) = force_d(m,0,j,0) + force_d(m,nx,j,nz_d)

            force_d(m,nx,j,1) = force_d(m,0,j,nz_d+1) + force_d(m,nx,j,1)

            force_d(m,1,j,nz_d) = force_d(m,nx+1,j,0) + force_d(m,1,j,nz_d)

            force_d(m,1,j,1) = force_d(m,nx+1,j,nz_d+1) + force_d(m,1,j,1)
          enddo
        endif
      endif
    end subroutine bc_edge_y_backf
    
    attributes(global) subroutine bc_BGK_edge_x(step, flip, faiPops,faiRho)
      integer, value :: step,flip
      logical, value :: faiPops, faiRho
      integer :: i, l
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1

      if (1<=i .and. i<=nx) then
        if (faiPops) then
          do l = 1, npops-1
            popsR_d(i,0,0, l, flip) = popsR_d(i,ny,nz_d, l, flip)

            popsR_d(i,0,nz_d+1, l, flip) = popsR_d(i,ny,1, l, flip)

            popsR_d(i,ny+1,0, l, flip) = popsR_d(i,1,nz_d, l, flip)

            popsR_d(i,ny+1,nz_d+1, l, flip) = popsR_d(i,1,1, l, flip)
          end do
        endif

        if (faiRho) then
          rhoR_d (i,0,0) = rhoR_d (i,ny,nz_d)

          rhoR_d (i,0,nz_d+1) = rhoR_d (i,ny,1)

          rhoR_d (i,ny+1,0) = rhoR_d (i,1,nz_d)

          rhoR_d (i,ny+1,nz_d+1) = rhoR_d (i,1,1)
        endif
      endif
    end subroutine bc_BGK_edge_x
    
    attributes(global) subroutine bc_edge_x_hvar(step,faiRho,faiRhoB,faiVel)
      integer, value :: step
      logical, value :: faiRho,faiRhoB,faiVel
      integer :: i,j,k, l,m
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x

      if (i<=nx) then

        if (faiRho) then
          do k=1,nbuff
            do j=1,nbuff
              rhoR_d(i,1-j,1-k) = rhoR_d(i,ny+1-j,nz_d+1-k)
              rhoR_d(i,1-j,nz_d+k) = rhoR_d(i,ny+1-j,k)
              rhoR_d(i,ny+j,1-k) = rhoR_d(i,j,nz_d+1-k)
              rhoR_d(i,ny+j,nz_d+k) = rhoR_d(i, j, k)
            enddo
          enddo
        endif
        if (faiRhoB) then
          do k=1,nbuff
            do j=1,nbuff
              rhoB_d(i,1-j,1-k) = rhoB_d(i,ny+1-j,nz_d+1-k)
              rhoB_d(i,1-j,nz_d+k) = rhoB_d(i,ny+1-j,k)
              rhoB_d(i,ny+j,1-k) = rhoB_d(i,j,nz_d+1-k)
              rhoB_d(i,ny+j,nz_d+k) = rhoB_d(i, j, k)
            enddo
          enddo
        endif
        if (faiVel) then
          do k=1,nbuff
            do j=1,nbuff
              do m=1,3
                vel_d(m,i,1-j,1-k) = vel_d(m,i,ny+1-j,nz_d+1-k)
                vel_d(m,i,1-j,nz_d+k) = vel_d(m,i,ny+1-j,k)
                vel_d(m,i,ny+j,1-k) = vel_d(m,i,j,nz_d+1-k)
                vel_d(m,i,ny+j,nz_d+k) = vel_d(m,i, j, k)
              enddo
            enddo
          enddo
        endif
      endif
    end subroutine bc_edge_x_hvar
    
     attributes(global) subroutine bc_edge_x_backf(step,faiF)
      integer, value :: step
      logical, value :: faiF
      integer :: i, l,m
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1

      if (1<=i .and. i<=nx) then
        
        if (faiF) then
          do m=1,3
            force_d(m,i,ny,nz_d) = force_d(m,i,0,0) + force_d(m,i,ny,nz_d)

            force_d(m,i,ny,1) = force_d(m,i,0,nz_d+1) + force_d(m,i,ny,1)

            force_d(m,i,1,nz_d) = force_d(m,i,ny+1,0) + force_d(m,i,1,nz_d)

            force_d(m,i,1,1) = force_d(m,i,ny+1,nz_d+1) + force_d(m,i,1,1)
          enddo
        endif
      endif
    end subroutine bc_edge_x_backf
    
    attributes(global) subroutine bc_BGK_corners(step, flip)
      integer, value :: step,flip
      integer :: j,k, l
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
    end subroutine bc_BGK_corners
    
    attributes(global) subroutine applybchalo_BGK(step, flip)
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
      
      ! Streaming      
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
               
                popsR_d(i,j,k,l, flip) = -popsR_d(i1,j1,k1,lopp, flip) + &
                 TWO* rho*p_d(lopp)*(ONE+HALF*(uv*uv)-(HALF/cssq)*(u*u + v*v + w*w))
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
       
    end subroutine applybchalo_BGK

    attributes(global) subroutine stream_BGK_old(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l, flop
      integer :: i1,j1,k1


      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA stream_BGK] i,k:',i,j,k

      flop = 3 - flip

      ! Streaming      
      if (myfluid_d(i,j,k, flip) < fluid_dead) then
        do l = 1, npops-1
          i1 = i + ex_d(l)
          j1 = j + ey_d(l)
          k1 = k + ez_d(l)
          popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
        enddo
      endif    
    end subroutine stream_BGK_old

    attributes(global) subroutine stream_BGK_x(step, flip)
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
    if (myfluid_d(i,j,k, flip) < fluid_dead) then
      do l = 1, npops-1
        i1 = i + ex_d(l)
        j1 = j + ey_d(l)
        k1 = k + ez_d(l)
        if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz_d) then
          popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
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
        endif
      end do
    endif
  end subroutine stream_BGK_z   
  
  attributes(global) subroutine stream_BGK_newbc(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l, flop, lopp
      integer :: i1,j1,k1


      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA stream_BGK] i,k:',i,j,k

      flop = 3 - flip

      ! Streaming      
      if (myfluid_d(i,j,k, flip) == fluid_fluid) then
        do l = 1, npops-1
          i1 = i + ex_d(l)
          j1 = j + ey_d(l)
          k1 = k + ez_d(l)
          if (myfluid_d(i1,j1,k1, flip) == fluid_fluid) then
            popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
          elseif(myfluid_d(i1,j1,k1, flip) == fluid_wall)then
            lopp=opp_d(l)
            popsR_d(i,j,k,lopp, flop) = popsR_d(i,j,k,l, flip)
          endif
        enddo
      endif    
    end subroutine stream_BGK_newbc

    attributes(global) subroutine stream_BGK_x_newbc(step, flip)
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
      if (myfluid_d(i,j,k, flip) == fluid_fluid) then
        do l = 1, npops-1
          i1 = i + ex_d(l)
          j1 = j + ey_d(l)
          k1 = k + ez_d(l)
          if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz_d) then
            if (myfluid_d(i1,j1,k1, flip) == fluid_fluid) then
              popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
            endif
          endif
        end do
      endif

      ! Stream x=nx+1
      i = nx+1
      if (myfluid_d(i,j,k, flip) == fluid_fluid) then
        do l = 1, npops-1
          i1 = i + ex_d(l)
          j1 = j + ey_d(l)
          k1 = k + ez_d(l)
          if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz_d) then
            if (myfluid_d(i1,j1,k1, flip) == fluid_fluid) then
              popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
            endif
          endif
        end do
      endif
    end subroutine stream_BGK_x_newbc

    attributes(global) subroutine stream_BGK_y_newbc(step, flip)
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
    if (myfluid_d(i,j,k, flip) == fluid_fluid) then
      do l = 1, npops-1
        i1 = i + ex_d(l)
        j1 = j + ey_d(l)
        k1 = k + ez_d(l)
        if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz_d) then
          if (myfluid_d(i1,j1,k1, flip) == fluid_fluid) then
            popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
          endif
        endif
      end do
    endif

    ! Stream y=ny+1
    j = ny+1
    if (myfluid_d(i,j,k, flip) == fluid_fluid) then
      do l = 1, npops-1
        i1 = i + ex_d(l)
        j1 = j + ey_d(l)
        k1 = k + ez_d(l)
        if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz_d) then
          if (myfluid_d(i1,j1,k1, flip) == fluid_fluid) then
            popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
          endif
        endif
      end do
    endif
  end subroutine stream_BGK_y_newbc

  attributes(global) subroutine stream_BGK_z_newbc(step, flip)
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
    if (myfluid_d(i,j,k, flip) == fluid_fluid) then
      do l = 1, npops-1
        i1 = i + ex_d(l)
        j1 = j + ey_d(l)
        k1 = k + ez_d(l)
        if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz_d) then
          if (myfluid_d(i1,j1,k1, flip) == fluid_fluid) then
            popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
          endif
        endif
      end do
    endif

    ! Stream z=nz_d+1
    k = nz_d+1
    if (myfluid_d(i,j,k, flip) == fluid_fluid) then
      do l = 1, npops-1
        i1 = i + ex_d(l)
        j1 = j + ey_d(l)
        k1 = k + ez_d(l)
        if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz_d) then
          if (myfluid_d(i1,j1,k1, flip) == fluid_fluid) then
            popsR_d(i1,j1,k1,l, flop) = popsR_d(i,j,k,l, flip)
          endif
        endif
      end do
    endif
  end subroutine stream_BGK_z_newbc   
  
  attributes(global) subroutine stream_BGK(step, flip)
      integer, value :: step,flip
      integer :: i,j,k,flop

      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      if (step==1 .and. 1==i*j*k+offset_d(3)) write(*,*) 'CUDA streamR_CG]'

      flop = 3 - flip

      ! Streaming
      if (myfluid_d(i,j,k, flip) < fluid_dead) then
          popsR_d(i  ,j  ,k  , 0, flop) = popsR_d(i,j,k, 0, flip)
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
    end subroutine stream_BGK 
    
    attributes(global) subroutine flipflopPop0_BGK(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, flop

      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      if (step==1 .and. 1==i*j*k+offset_d(3)) write(*,*) 'CUDA flipflopRPop0_CG]'

      flop = 3 - flip
      if(myfluid_d(i,j,k,flip) < fluid_dead)then
        popsR_d(i,j,k,0, flop) = popsR_d(i,j,k,0, flip)
      endif
    end subroutine flipflopPop0_BGK

  end module kernels_bgk
