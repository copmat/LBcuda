#include "defines.h"

#ifdef SCPD3Q27
#define manual_rho(pops) pops(m,i,j,k,0,flip) + pops(m,i,j,k,1,flip) + pops(m,i,j,k,2, flip) +   \
       pops(m,i,j,k,3, flip) + pops(m,i,j,k,4, flip) + pops(m,i,j,k,5, flip) + pops(m,i,j,k,6, flip) + pops(m,i,j,k,7, flip) +           \
       pops(m,i,j,k,8, flip) + pops(m,i,j,k,9, flip) + pops(m,i,j,k,10,flip) + pops(m,i,j,k,11,flip) + pops(m,i,j,k,12,flip) +           \
       pops(m,i,j,k,13, flip) + pops(m,i,j,k,14, flip) + pops(m,i,j,k,15, flip) + pops(m,i,j,k,16, flip) + pops(m,i,j,k,17, flip) +      \
       pops(m,i,j,k,18, flip) + pops(m,i,j,k,19, flip) + pops(m,i,j,k,20, flip) + pops(m,i,j,k,21, flip) + pops(m,i,j,k,22, flip) +      \
       pops(m,i,j,k,23, flip) + pops(m,i,j,k,24, flip) + pops(m,i,j,k,25, flip) + pops(m,i,j,k,26, flip)
#elif defined SCPD3Q19
#define manual_rho(pops) pops(m,i,j,k,0,flip) + pops(m,i,j,k,1,flip) + pops(m,i,j,k,2, flip) +   \
       pops(m,i,j,k,3, flip) + pops(m,i,j,k,4, flip) + pops(m,i,j,k,5, flip) + pops(m,i,j,k,6, flip) + pops(m,i,j,k,7, flip) +           \
       pops(m,i,j,k,8, flip) + pops(m,i,j,k,9, flip) + pops(m,i,j,k,10,flip) + pops(m,i,j,k,11,flip) + pops(m,i,j,k,12,flip) +           \
       pops(m,i,j,k,13, flip) + pops(m,i,j,k,14, flip) + pops(m,i,j,k,15, flip) + pops(m,i,j,k,16, flip) + pops(m,i,j,k,17, flip) +      \
       pops(m,i,j,k,18, flip)
#else
#define manual_rho(pops) pops(m,i,j,k,0,flip) + pops(m,i,j,k,1,flip) + pops(m,i,j,k,2, flip) +   \
       pops(m,i,j,k,3, flip) + pops(m,i,j,k,4, flip) + pops(m,i,j,k,5, flip) + pops(m,i,j,k,6, flip)
#endif

  module kernels_scp
    use dimensions_m
    use kernels_fluid
    implicit none
    
!    private :: equil,compute_u_1fl,compute_v_1fl, &
!     compute_w_1fl
!    private :: compute_rho
  contains
   
    attributes(device) function equil_SCP(rho,u,v,w, l,m)
     real, intent(in) :: rho(numscp),u,v,w
     integer, intent(in) :: l,m
     real :: equil_SCP
     real :: uv

     uv = (1.0/cssq) * (u*ex_d(l) + v*ey_d(l) + w*ez_d(l))
     equil_SCP = rho(m) * p_scp(l)*(ONE+uv+HALF*(uv*uv)-(HALF/cssq) * (u*u + v*v + w*w))
    end function equil_SCP
 
    ! Setup
    attributes(global) subroutine setup_SCP(vx,vy,vz)
      real, value :: vx,vy,vz
      integer :: i,j,k, l,m,gli,glj,glk
      real    :: rho(numscp), u,v,w, eq, x,y,z
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
      if(i>nx .or. j>ny .or. k>nz_d)return

      ! if (k>nz/2) then
      !   rho = 1.01
      ! else
      !   rho = 0.99
      ! endif

      x = 3.14159265358979 / nx * i
      y = 3.14159265358979 / ny * j
      z = 3.14159265358979 / nz_d * k
      rho = 1.0 !densR_d !1 + 0.1*sin(x)*sin(x) * sin(y)*sin(y) * sin(z)*sin(z)
      
      gli = i + offset_d(1)
      glj = j + offset_d(2)
      glk = k + offset_d(3)
      
      
      u = vx
      v = vy
      w = vz
      
      if(store_vel_d)then
        u = vel_d(1,i,j,k)
        v = vel_d(2,i,j,k)
        w = vel_d(3,i,j,k)
      endif
      
      ! write(*,*) 'CUDA setup',rho,i,j,k
      do m=1,numscp
        scalar_d(m,i,j,k) = rho(m)
      enddo
      do l = 0, npops_scp-1
        do m = 1, numscp
          eq = equil_SCP(rho, u,v,w, l,m)
          popsSCP_d(m,i,j,k, l, 1) = eq
        enddo
      end do
      
      
      
    end subroutine setup_SCP
    
    attributes(global) subroutine setup_SCP_p(vx,vy,vz)
      real, value :: vx,vy,vz
      integer :: i,j,k, l,m,gli,glj,glk
      real    :: rho(numscp), u,v,w, eq, x,y,z
  
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
      ! write(*,*) 'CUDA setup',rho,i,j,k
      do m=1,numscp
        scalar_d(m,i,j,k) = rho(m)
      enddo
      do l = 0, npops_scp-1
        do m = 1, numscp
          eq = equil_SCP(rho, u,v,w, l,m)
          popsSCP_d(m,i,j,k, l, 1) = eq
        enddo
      end do

     
    end subroutine setup_SCP_p




    ! Time Step
    attributes(global) subroutine time_step_mom_SCP(step, flip, omega)
      integer, value :: step,flip
      real, value :: omega
      real    :: invrho, u,v,w,oneminusomega,omegaminusone,ush,vsh,wsh, &
       myscalar(numscp)
      integer :: i,j,k, l,m, flop
      integer :: i1,j1,k1
      
      oneminusomega = 1.0 -omega
      omegaminusone = omega - 1.0
      
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA time_step] i,k:',i,j,k
      
      if (myfluid_d(i,j,k, flip) == fluid_fluid) then
        invrho = 1.0 / rhoR_d(i,j,k)
        do m=1,numscp
          myscalar(m) = manual_rho(popsSCP_d)
          if (myscalar(m)<minSCP .or. myscalar(m)>maxSCP) then
           write(*,*) 'myscalar(m) error:',myscalar(m),i,j,k,m
           stop_d = __LINE__
          endif
          scalar_d(m,i,j,k) = myscalar(m)
        enddo
    

        
        u = vel_d(1,i,j,k)
        v = vel_d(2,i,j,k)
        w = vel_d(3,i,j,k)
        
        ush = u + fx_d * invrho
        vsh = v + fy_d * invrho
        wsh = w + fz_d * invrho
  
        flop = 3 - flip
        do l = 0, npops_scp-1
          !i1 = i + ex_d(l)
          !j1 = j + ey_d(l)
          !k1 = k + ez_d(l)
          do m=1,numscp
            !popsSCP_d(m,i1,j1,k1, l, flop) = popsSCP_d(m,i,j,k, l, flip)*oneminusomega + equil_SCP(myscalar, u,v,w, l,m)*omega
            popsSCP_d(m,i,j,k,l,flip) = popsSCP_d(m,i,j,k,l,flip)*oneminusomega + &
             equil_SCP(myscalar, u,v,w, l,m)*omegaminusone + equil_SCP(myscalar, ush,vsh,wsh, l,m)
          enddo
        end do
      
      else
        do m=1,numscp
          scalar_d(m,i,j,k) = MINDENS
        enddo
      endif
      
    end subroutine time_step_mom_SCP

  

    attributes(global) subroutine init_rho_isfluid_SCP(step, flip)
    integer, value :: step,flip
    real    :: myscalar(numscp)
    integer :: i,j,k, l,m

    i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
    j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
    k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
    ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA init_rho_isfluid_SCP] i,k:',i,j,k

    ! Only on fluid nodes
    if (myfluid_d(i,j,k, flip) == fluid_fluid) then
      do m=1,numscp
        myscalar(m) = manual_rho(popsSCP_d)

        if (myscalar(m)<minSCP .or. myscalar(m)>maxSCP) then
          ! write(*,*) 'init_rho_isfluid_SCP]Range error scalar', step, linear(i,j,k),m, myscalar(m)
          stop_d = __LINE__
        endif

        scalar_d(m,i,j,k) = myscalar(m)
      enddo
    else
      do m=1,numscp
        scalar_d(m,i,j,k) = MINDENS
      enddo
    endif
    end subroutine init_rho_isfluid_SCP
    
 



    attributes(global) subroutine time_step_SCP(step, flip, omega)
      integer, value :: step,flip
      real, value :: omega
      real    :: myscalar(numscp), invrho, u,v,w  ,oneminusomega,omegaminusone,ush,vsh,wsh
      integer :: i,j,k, l,m
      
      oneminusomega = 1.0 - omega
      omegaminusone = omega - 1.0
      
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA time_step_CG] i,k:',i,j,k

      if (myfluid_d(i,j,k, flip) == fluid_fluid) then
        invrho = 1.0 / rhoR_d(i,j,k)
        do m=1,numscp
          myscalar(m) = scalar_d(m,i,j,k)
        enddo

        u = vel_d(1,i,j,k)
        v = vel_d(2,i,j,k)
        w = vel_d(3,i,j,k)
        
        ush = u + fx_d * invrho
        vsh = v + fy_d * invrho
        wsh = w + fz_d * invrho
        
        !scp step
        do l = 0, npops_scp-1
          do m=1,numscp
            popsSCP_d(m,i,j,k,l, flip) = popsSCP_d(m,i,j,k, l, flip)*oneminusomega + &
             equil_SCP(myscalar, u,v,w, l,m)*omegaminusone + equil_SCP(myscalar, ush,vsh,wsh, l,m)
          enddo
        enddo
        
      endif
    end subroutine time_step_SCP

    attributes(global) subroutine copy_SCP(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l, m,flop

      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA copy_SCP] i,k:',i,j,k

      flop = 3 - flip
      
      do l = 0, npops_scp-1
        do m=1,numscp
          popsSCP_d(m,i,j,k,l, flop) = popsSCP_d(m,i,j,k,l, flip)
        enddo
      enddo
    end subroutine copy_SCP
    
    attributes(global) subroutine bc_SCP_per_x(step, flip, faiPops,faiRho)
      integer, value :: step,flip
      logical, value :: faiPops, faiRho
      integer :: j,k, l,m
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
      ! if (1==step .and. 1==j*k) write(*,*) 'CUDA bc_per_x2] i,j:',j,k
      
      if (1<=j .and. j<=ny .and. 1<=k .and. k<=nz_d) then

         
        if (faiPops) then
          do l = 1, npops_scp-1
            do m=1,numscp
              popsSCP_d(m,nx+1, j,k, l, flip) = popsSCP_d(m,1,j,k, l, flip)
            
              popsSCP_d(m,0,j,k, l, flip) = popsSCP_d(m,nx,j,k, l, flip)
            enddo
          end do
        endif

        if (faiRho) then
          do m=1,numscp
            scalar_d (m,nx+1,j,k) = scalar_d (m,1,j,k)

            scalar_d (m,0,j,k) = scalar_d (m,nx,j,k)
          enddo
        endif

      endif
      
    end subroutine bc_SCP_per_x
    
    attributes(global) subroutine bc_per_x_svar(step,faiRho)
      integer, value :: step
      logical, value :: faiRho
      integer :: i,j,k, l,m
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y
      ! if (1==step .and. 1==j*k) write(*,*) 'CUDA bc_per_x2] i,j:',j,k
      
      if (j<=ny .and. k<=nz_d) then

        if (faiRho) then
          do i=1,nbuff
            do m=1,numscp
              scalar_d (m,nx+i,j,k) = scalar_d (m,i,j,k)

              scalar_d (m,1-i,j,k) = scalar_d (m,nx+1-i,j,k)
            enddo
          enddo
        endif

      endif
      
    end subroutine bc_per_x_svar

    attributes(global) subroutine bc_SCP_per_y(step, flip, faiPops,faiRho)
      integer, value :: step,flip
      logical, value :: faiPops, faiRho
      integer :: i,k, l,m
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
      ! if (1==step .and. 1==i*k) write(*,*) 'CUDA bc_per_y2] i,k:',i,k
      
      if (1<=i .and. i<=nx .and. 1<=k .and. k<=nz_d) then

          
        if (faiPops) then
          do l = 1, npops_scp-1
            do m=1,numscp
              popsSCP_d(m,i,ny+1,k, l, flip) = popsSCP_d(m,i,1,k, l, flip)
            
              popsSCP_d(m,i,0,k, l, flip) = popsSCP_d(m,i,ny,k, l, flip)
            enddo
          end do
        endif

        if (faiRho) then
          do m=1,numscp
            scalar_d(m,i,ny+1,k) = scalar_d(m,i,1,k)

            scalar_d(m,i,0,k) = scalar_d(m,i,ny,k)
          enddo
        endif
        
      endif

      
      
      
    end subroutine bc_SCP_per_y
    
    attributes(global) subroutine bc_per_y_svar(step,faiRho)
      integer, value :: step
      logical, value :: faiRho
      integer :: i,j,k, l,m
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y
      ! if (1==step .and. 1==i*k) write(*,*) 'CUDA bc_per_y2] i,k:',i,k
      
      if (i<=nx .and. k<=nz_d) then

          

        if (faiRho) then
          do j=1,nbuff
            scalar_d (m,i,ny+j,k) = scalar_d (m,i,j,k)

            scalar_d (m,i,1-j,k) = scalar_d (m,i,ny+1-j,k)
          enddo
        endif
        
      endif

      
      
      
    end subroutine bc_per_y_svar
    
    attributes(global) subroutine bc_SCP_per_z(step, flip, faiPops,faiRho)
      integer, value :: step,flip
      logical, value :: faiPops, faiRho
      integer :: i,j, l,m
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      j = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
      ! if (1==step .and. 1==i*j) write(*,*) 'CUDA bc_per_z2] i,j:',i,j

      if (1<=j .and. j<=ny .and. 1<=i .and. i<=nx) then
        if (faiPops) then
          do l = 1, npops_scp-1   
            do m=1,numscp     
              popsSCP_d(m,i,j,nz_d+1, l, flip) = popsSCP_d(m,i,j,1, l, flip)
          
              popsSCP_d(m,i,j,0, l, flip) = popsSCP_d(m,i,j,nz_d, l, flip)
            enddo
          end do
        endif

        if (faiRho) then
          do m=1,numscp
            scalar_d(m,i,j,nz_d+1) = scalar_d(m,i,j,1)

            scalar_d(m,i,j,0) = scalar_d(m,i,j,nz_d)
          enddo
        endif
      endif
    end subroutine bc_SCP_per_z
    
    attributes(global) subroutine bc_per_z_svar(step,faiRho)
      integer, value :: step
      logical, value :: faiRho
      integer :: i,j,k, l,m
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIM + threadIdx%y
      ! if (1==step .and. 1==i*j) write(*,*) 'CUDA bc_per_z2] i,j:',i,j

      if (j<=ny .and. i<=nx) then

        if (faiRho) then
          do k=1,nbuff
            do m=1,numscp
              scalar_d (m,i,j,nz_d+k) = scalar_d (m,i,j,k)

              scalar_d (m,i,j,1-k) = scalar_d (m,i,j,nz_d+1-k)
            enddo
          enddo
        endif

      endif
    end subroutine bc_per_z_svar
    
    attributes(global) subroutine bc_SCP_edge_z(step, flip, faiPops,faiRho)
      integer, value :: step,flip
      logical, value :: faiPops, faiRho
      integer :: k, l,m
  
      k = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1

      if (0<=k .and. k<=nz_d+1) then
        if (faiPops) then
          do l = 1, npops_scp-1
            do m=1,numscp
              popsSCP_d(m,0,0,k, l, flip) = popsSCP_d(m,nx,ny,k, l, flip)

              popsSCP_d(m,0,ny+1,k, l, flip) = popsSCP_d(m,nx,1,k, l, flip)

              popsSCP_d(m,nx+1,0,k, l, flip) = popsSCP_d(m,1,ny,k, l, flip)

              popsSCP_d(m,nx+1,ny+1,k, l, flip) = popsSCP_d(m,1,1,k, l, flip)
            enddo
          end do
        endif

        if (faiRho) then
          do m=1,numscp
            scalar_d (m,0,0,k) = scalar_d (m,nx,ny,k)

            scalar_d (m,0,ny+1,k) = scalar_d (m,nx,1,k)

            scalar_d (m,nx+1,0,k) = scalar_d (m,1,ny,k)

            scalar_d (m,nx+1,ny+1,k) = scalar_d (m,1,1,k)
          enddo
        endif
      endif
    end subroutine bc_SCP_edge_z
    
    attributes(global) subroutine bc_edge_z_svar(step,faiRho)
      integer, value :: step
      logical, value :: faiRho
      integer :: k,i,j, l,m
  
      k = (blockIdx%x-1) * TILE_DIM + threadIdx%x  - nbuff

      if (k<=nz_d+nbuff) then

        if (faiRho) then
          do j=1,nbuff
            do i=1,nbuff
              do m=1,numscp
                scalar_d(m,1-i,1-j,k) = scalar_d(m,nx+1-i,ny+1-j,k)
                scalar_d(m,1-i,ny+j,k) = scalar_d(m,nx+1-i,j,k)
                scalar_d(m,nx+i,1-j,k) = scalar_d(m, i,ny+1-j,k)
                scalar_d(m,nx+i,ny+j,k) = scalar_d(m, i, j,k)
              enddo
            enddo
          enddo
        endif

      endif
    end subroutine bc_edge_z_svar
    
    attributes(global) subroutine bc_SCP_edge_y(step, flip, faiPops,faiRho)
      integer, value :: step,flip
      logical, value :: faiPops, faiRho
      integer :: j, l,m
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1

      if (1<=j .and. j<=ny) then
        if (faiPops) then
          do l = 1, npops_scp-1
            do m=1,numscp
              popsSCP_d(m,0,j,0, l, flip) = popsSCP_d(m,nx,j,nz_d, l, flip)

              popsSCP_d(m,0,j,nz_d+1, l, flip) = popsSCP_d(m,nx,j,1, l, flip)

              popsSCP_d(m,nx+1,j,0, l, flip) = popsSCP_d(m,1,j,nz_d, l, flip)

              popsSCP_d(m,nx+1,j,nz_d+1, l, flip) = popsSCP_d(m,1,j,1, l, flip)
            enddo
          end do
        endif

        if (faiRho) then
          do m=1,numscp
            scalar_d(m,0,j,0) = scalar_d(m,nx,j,nz_d)

            scalar_d(m,0,j,nz_d+1) = scalar_d(m,nx,j,1)

            scalar_d(m,nx+1,j,0) = scalar_d(m,1,j,nz_d)

            scalar_d(m,nx+1,j,nz_d+1) = scalar_d(m,1,j,1)
          enddo
        endif
      endif
    end subroutine bc_SCP_edge_y
    
    attributes(global) subroutine bc_edge_y_svar(step,faiRho)
      integer, value :: step
      logical, value :: faiRho
      integer :: j,i,k, l,m
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x

      if (j<=ny) then

        if (faiRho) then
          do k=1,nbuff
            do i=1,nbuff
              do m=1,numscp
                scalar_d(m,1-i,j,1-k) = scalar_d(m,nx+1-i,j, nz_d+1-k)
                scalar_d(m,1-i,j,nz_d+k) = scalar_d(m,nx+1-i,j,  k)
                scalar_d(m,nx+i,j,1-k) = scalar_d(m, i,j, nz_d+1-k)
                scalar_d(m,nx+i,j,nz_d+k) = scalar_d(m, i,j,  k)
              enddo
            enddo
          enddo
        endif

      endif
    end subroutine bc_edge_y_svar
    
    
    attributes(global) subroutine bc_SCP_edge_x(step, flip, faiPops,faiRho)
      integer, value :: step,flip
      logical, value :: faiPops, faiRho
      integer :: i, l,m
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1

      if (1<=i .and. i<=nx) then
        if (faiPops) then
          do l = 1, npops_scp-1
            do m=1,numscp
              popsSCP_d(m,i,0,0, l, flip) = popsSCP_d(m,i,ny,nz_d, l, flip)

              popsSCP_d(m,i,0,nz_d+1, l, flip) = popsSCP_d(m,i,ny,1, l, flip)

              popsSCP_d(m,i,ny+1,0, l, flip) = popsSCP_d(m,i,1,nz_d, l, flip)

              popsSCP_d(m,i,ny+1,nz_d+1, l, flip) = popsSCP_d(m,i,1,1, l, flip)
            enddo
          end do
        endif

        if (faiRho) then
          do m=1,numscp
            scalar_d(m,i,0,0) = scalar_d(m,i,ny,nz_d)

            scalar_d(m,i,0,nz_d+1) = scalar_d(m,i,ny,1)

            scalar_d(m,i,ny+1,0) = scalar_d(m,i,1,nz_d)

            scalar_d(m,i,ny+1,nz_d+1) = scalar_d(m,i,1,1)
          enddo
        endif
      endif
    end subroutine bc_SCP_edge_x
    
    attributes(global) subroutine bc_edge_x_svar(step,faiRho)
      integer, value :: step
      logical, value :: faiRho
      integer :: i,j,k, l,m
  
      i = (blockIdx%x-1) * TILE_DIM + threadIdx%x

      if (i<=nx) then

        if (faiRho) then
          do k=1,nbuff
            do j=1,nbuff
              do m=1,numscp
                scalar_d(m,i,1-j,1-k) = scalar_d(m,i,ny+1-j,nz_d+1-k)
                scalar_d(m,i,1-j,nz_d+k) = scalar_d(m,i,ny+1-j,k)
                scalar_d(m,i,ny+j,1-k) = scalar_d(m,i,j,nz_d+1-k)
                scalar_d(m,i,ny+j,nz_d+k) = scalar_d(m,i, j, k)
              enddo
            enddo
          enddo
        endif

      endif
    end subroutine bc_edge_x_svar
    
    
    attributes(global) subroutine bc_SCP_corners(step, flip)
      integer, value :: step,flip
      integer :: j,k, l
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
    end subroutine bc_SCP_corners
    
    attributes(global) subroutine applybc_SCP(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l, flop, lopp, itemp, mytype,m
      integer :: i1,j1,k1
      real :: u,v,w,uv,rho(numscp)

      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA stream_SCP] i,k:',i,j,k

      flop = 3 - flip
#if 0
      ! Streaming      
      if (myfluid_d(i,j,k, flip) <= fluid_wall) then
        itemp=-int(myfluid_d(i,j,k, flip))
        if(itemp>0)then
          mytype=bcscptype_d(itemp)
          if(mytype==1)then
            do m=1,numscp
              rho(m) = bcscp_d(m,itemp)
            enddo
            do l = 1, npops_scp-1
              i1 = i + ex_d(l)
              j1 = j + ey_d(l)
              k1 = k + ez_d(l)
              if (myfluid_d(i1,j1,k1, flip) == fluid_fluid) then
                lopp=opp_d(l)
                u = vel_d(1,i1,j1,k1)
                v = vel_d(2,i1,j1,k1)
                w = vel_d(3,i1,j1,k1)
                
                uv = (ONE/cssq) * (u*ex_d(lopp) + v*ey_d(lopp) + w*ez_d(lopp))
                do m=1,numscp
                  popsSCP_d(m,i,j,k,l, flip) = -popsSCP_d(m,i1,j1,k1,lopp, flip) + &
                   TWO* rho(m)*p_d(lopp)*(ONE+HALF*(uv*uv)-(HALF/cssq)*(u*u + v*v + w*w))
                enddo
              endif
            enddo
          elseif(mytype==2)then
            u = bcvel_d(1,itemp)
            v = bcvel_d(2,itemp)
            w = bcvel_d(3,itemp)
            do l = 1, npops_scp-1
              i1 = i + ex_d(l)
              j1 = j + ey_d(l)
              k1 = k + ez_d(l)
              if (myfluid_d(i1,j1,k1, flip) == fluid_fluid) then
                lopp=opp_d(l)
                do m=1,numscp
                  rho(m) = scalar_d(m,i1,j1,k1)
                enddo
                uv = (ONE/cssq) * (u*ex_d(lopp) + v*ey_d(lopp) + w*ez_d(lopp))
                do m=1,numscp
                  popsSCP_d(m,i,j,k,l, flip) = popsSCP_d(m,i1,j1,k1,lopp, flip) - &
                   TWO* rho(m)*p_d(lopp)*uv
                enddo
              endif
            enddo
          elseif(mytype==3)then
            do m=1,numscp
              rho(m) = bcscp_d(m,itemp)
            enddo
            u = bcvel_d(1,itemp)
            v = bcvel_d(2,itemp)
            w = bcvel_d(3,itemp)
            do l = 1, npops_scp-1
              i1 = i + ex_d(l)
              j1 = j + ey_d(l)
              k1 = k + ez_d(l)
              if (myfluid_d(i1,j1,k1, flip) == fluid_fluid) then
                lopp=opp_d(l)
                
                uv = (ONE/cssq) * (u*ex_d(lopp) + v*ey_d(lopp) + w*ez_d(lopp))
                do m=1,numscp
                  popsSCP_d(m,i,j,k,l, flip) = popsSCP_d(m,i1,j1,k1,lopp, flip) - &
                   TWO* rho(m)*p_d(lopp)*uv
                enddo
              endif
            enddo
          else
            do l = 1, npops_scp-1
              i1 = i + ex_d(l)
              j1 = j + ey_d(l)
              k1 = k + ez_d(l)
              if (myfluid_d(i1,j1,k1, flip) == fluid_fluid) then
                lopp=opp_d(l)
                do m=1,numscp
                  popsSCP_d(m,i,j,k,l, flip) = popsSCP_d(m,i1,j1,k1,lopp, flip)
                enddp
              endif
            enddo
          endif
        else
          do l = 1, npops_scp-1
            i1 = i + ex_d(l)
            j1 = j + ey_d(l)
            k1 = k + ez_d(l)
            if (myfluid_d(i1,j1,k1, flip) == fluid_fluid) then
              lopp=opp_d(l)
              do m=1,numscp
                popsSCP_d(m,i,j,k,l, flip) = popsSCP_d(m,i1,j1,k1,lopp, flip)
              enddo
            endif
          enddo
        endif
      endif  
#else
      if (myfluid_d(i,j,k, flip) == fluid_wall) then
        do l = 1, npops_scp-1
            i1 = i + ex_d(l)
            j1 = j + ey_d(l)
            k1 = k + ez_d(l)
            if (myfluid_d(i1,j1,k1, flip) == fluid_fluid) then
              lopp=opp_d(l)
              do m=1,numscp
                popsSCP_d(m,i,j,k,l, flip) = popsSCP_d(m,i1,j1,k1,lopp, flip)
              enddo
            endif
          enddo
      endif
#endif  
    end subroutine applybc_SCP
    
    attributes(global) subroutine applybchalo_SCP(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l, flop, lopp, itemp, mytype,m
      integer :: i1,j1,k1
      real :: u,v,w,uv,rho

      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x - 1
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y - 1
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z - 1
      ! if (step==1 .and. 1==i*j*k) write(*,*) 'CUDA stream_SCP] i,k:',i,j,k

      flop = 3 - flip
      
      if (i>nx+1) return
      if (j>ny+1) return
      if (k>nz_d+1) return
 
      if (myfluid_d(i,j,k, flip) == fluid_wall) then
        do l = 1, npops_scp-1
            i1 = i + ex_d(l)
            j1 = j + ey_d(l)
            k1 = k + ez_d(l)
            if(i1<1 .or. i1>nx)cycle
            if(j1<1 .or. j1>ny)cycle
            if(k1<1 .or. k1>nz_d)cycle
            if (myfluid_d(i1,j1,k1, flip) == fluid_fluid) then
              lopp=opp_d(l)
              do m=1,numscp
                popsSCP_d(m,i,j,k,l, flip) = popsSCP_d(m,i1,j1,k1,lopp, flip)
              enddo
            endif
          enddo
      endif
 
    end subroutine applybchalo_SCP

    attributes(global) subroutine stream_SCP_x(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l, flop,m
      integer :: i1,j1,k1
  
      j = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
      k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1

      if (j>ny+1) return
      if (k>nz_d+1) return

      flop = 3 - flip

      ! Stream x=0
      i = 0
      if (myfluid_d(i,j,k, flip) < fluid_dead) then
        do l = 1, npops_scp-1
          i1 = i + ex_d(l)
          j1 = j + ey_d(l)
          k1 = k + ez_d(l)
          if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz_d) then
            do m=1,numscp
              popsSCP_d(m,i1,j1,k1,l, flop) = popsSCP_d(m,i,j,k,l, flip)
            enddo
          endif
        end do
      endif

      ! Stream x=nx+1
      i = nx+1
      if (myfluid_d(i,j,k, flip) < fluid_dead) then
        do l = 1, npops_scp-1
          i1 = i + ex_d(l)
          j1 = j + ey_d(l)
          k1 = k + ez_d(l)
          if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz_d) then
            do m=1,numscp
              popsSCP_d(m,i1,j1,k1,l, flop) = popsSCP_d(m,i,j,k,l, flip)
            enddo
          endif
        end do
      endif
    end subroutine stream_SCP_x

    attributes(global) subroutine stream_SCP_y(step, flip)
    integer, value :: step,flip
    integer :: i,j,k, l, flop,m
    integer :: i1,j1,k1

    i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
    k = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
    
    if (i>nx+1) return
    if (k>nz_d+1) return

    flop = 3 - flip

    ! Stream y=0
    j = 0
    if (myfluid_d(i,j,k, flip) < fluid_dead) then
      do l = 1, npops_scp-1
        i1 = i + ex_d(l)
        j1 = j + ey_d(l)
        k1 = k + ez_d(l)
        if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz_d) then
          do m=1,numscp
            popsSCP_d(m,i1,j1,k1,l, flop) = popsSCP_d(m,i,j,k,l, flip)
          enddo
        endif
      end do
    endif

    ! Stream y=ny+1
    j = ny+1
    if (myfluid_d(i,j,k, flip) < fluid_dead) then
      do l = 1, npops_scp-1
        i1 = i + ex_d(l)
        j1 = j + ey_d(l)
        k1 = k + ez_d(l)
        if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz_d) then
          do m=1,numscp
            popsSCP_d(m,i1,j1,k1,l, flop) = popsSCP_d(m,i,j,k,l, flip)
          enddo
        endif
      end do
    endif
  end subroutine stream_SCP_y

  attributes(global) subroutine stream_SCP_z(step, flip)
    integer, value :: step,flip
    integer :: i,j,k, l, flop,m
    integer :: i1,j1,k1

    i = (blockIdx%x-1) * TILE_DIM + threadIdx%x - 1
    j = (blockIdx%y-1) * TILE_DIM + threadIdx%y - 1
    
    if (i>nx+1) return
    if (j>ny+1) return

    flop = 3 - flip

    ! Stream z=0
    k = 0
    if (myfluid_d(i,j,k, flip) < fluid_dead) then
      do l = 1, npops_scp-1
        i1 = i + ex_d(l)
        j1 = j + ey_d(l)
        k1 = k + ez_d(l)
        if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz_d) then
          do m=1,numscp
            popsSCP_d(m,i1,j1,k1,l, flop) = popsSCP_d(m,i,j,k,l, flip)
          enddo
        endif
      end do
    endif

    ! Stream z=nz_d+1
    k = nz_d+1
    if (myfluid_d(i,j,k, flip) < fluid_dead) then
      do l = 1, npops_scp-1
        i1 = i + ex_d(l)
        j1 = j + ey_d(l)
        k1 = k + ez_d(l)
        if (1<=i1 .and. i1<=nx .and. 1<=j1 .and. j1<=ny .and. 1<=k1 .and. k1<=nz_d) then
          do m=1,numscp
            popsSCP_d(m,i1,j1,k1,l, flop) = popsSCP_d(m,i,j,k,l, flip)
          enddo
        endif
      end do
    endif
  end subroutine stream_SCP_z   
  
  attributes(global) subroutine stream_SCP(step, flip)
      integer, value :: step,flip
      integer :: i,j,k,flop,m

      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      if (step==1 .and. 1==i*j*k+offset_d(3)) write(*,*) 'CUDA streamR_CG]'

      flop = 3 - flip

      ! Streaming
      if (myfluid_d(i,j,k, flip) < fluid_dead) then
        do m=1,numscp
#ifdef SCPD3Q27
          popsSCP_d(m,i  ,j  ,k  , 0, flop) = popsSCP_d(m,i,j,k, 0, flip)
          popsSCP_d(m,i+1,j  ,k  , 1, flop) = popsSCP_d(m,i,j,k, 1,flip)
          popsSCP_d(m,i-1,j  ,k  , 2, flop) = popsSCP_d(m,i,j,k, 2,flip)
          popsSCP_d(m,i  ,j+1,k  , 3, flop) = popsSCP_d(m,i,j,k, 3,flip)
          popsSCP_d(m,i  ,j-1,k  , 4, flop) = popsSCP_d(m,i,j,k, 4,flip)
          popsSCP_d(m,i  ,j  ,k+1, 5, flop) = popsSCP_d(m,i,j,k, 5,flip)
          popsSCP_d(m,i  ,j  ,k-1, 6, flop) = popsSCP_d(m,i,j,k, 6,flip)
          popsSCP_d(m,i+1,j+1,k  , 7, flop) = popsSCP_d(m,i,j,k, 7,flip)
          popsSCP_d(m,i-1,j-1,k  , 8, flop) = popsSCP_d(m,i,j,k, 8,flip)
          popsSCP_d(m,i-1,j+1,k  , 9, flop) = popsSCP_d(m,i,j,k, 9,flip)
          popsSCP_d(m,i+1,j-1,k  ,10, flop) = popsSCP_d(m,i,j,k,10,flip)
          popsSCP_d(m,i+1,j  ,k+1,11, flop) = popsSCP_d(m,i,j,k,11,flip)
          popsSCP_d(m,i-1,j  ,k-1,12, flop) = popsSCP_d(m,i,j,k,12,flip)
          popsSCP_d(m,i-1,j  ,k+1,13, flop) = popsSCP_d(m,i,j,k,13,flip)
          popsSCP_d(m,i+1,j  ,k-1,14, flop) = popsSCP_d(m,i,j,k,14,flip)
          popsSCP_d(m,i  ,j+1,k+1,15, flop) = popsSCP_d(m,i,j,k,15,flip)
          popsSCP_d(m,i  ,j-1,k-1,16, flop) = popsSCP_d(m,i,j,k,16,flip)
          popsSCP_d(m,i  ,j-1,k+1,17, flop) = popsSCP_d(m,i,j,k,17,flip)
          popsSCP_d(m,i  ,j+1,k-1,18, flop) = popsSCP_d(m,i,j,k,18,flip)
          popsSCP_d(m,i+1,j+1,k+1,19, flop) = popsSCP_d(m,i,j,k,19,flip)
          popsSCP_d(m,i-1,j-1,k-1,20, flop) = popsSCP_d(m,i,j,k,20,flip)
          popsSCP_d(m,i+1,j+1,k-1,21, flop) = popsSCP_d(m,i,j,k,21,flip)
          popsSCP_d(m,i-1,j-1,k+1,22, flop) = popsSCP_d(m,i,j,k,22,flip)
          popsSCP_d(m,i-1,j+1,k+1,23, flop) = popsSCP_d(m,i,j,k,23,flip)
          popsSCP_d(m,i+1,j-1,k-1,24, flop) = popsSCP_d(m,i,j,k,24,flip)
          popsSCP_d(m,i+1,j-1,k+1,25, flop) = popsSCP_d(m,i,j,k,25,flip)
          popsSCP_d(m,i-1,j+1,k-1,26, flop) = popsSCP_d(m,i,j,k,26,flip)
#elif defined SCPD3Q19
          popsSCP_d(m,i  ,j  ,k  , 0, flop) = popsSCP_d(m,i,j,k, 0, flip)
          popsSCP_d(m,i+1,j  ,k  , 1, flop) = popsSCP_d(m,i,j,k, 1,flip)
          popsSCP_d(m,i-1,j  ,k  , 2, flop) = popsSCP_d(m,i,j,k, 2,flip)
          popsSCP_d(m,i  ,j+1,k  , 3, flop) = popsSCP_d(m,i,j,k, 3,flip)
          popsSCP_d(m,i  ,j-1,k  , 4, flop) = popsSCP_d(m,i,j,k, 4,flip)
          popsSCP_d(m,i  ,j  ,k+1, 5, flop) = popsSCP_d(m,i,j,k, 5,flip)
          popsSCP_d(m,i  ,j  ,k-1, 6, flop) = popsSCP_d(m,i,j,k, 6,flip)
          popsSCP_d(m,i+1,j+1,k  , 7, flop) = popsSCP_d(m,i,j,k, 7,flip)
          popsSCP_d(m,i-1,j-1,k  , 8, flop) = popsSCP_d(m,i,j,k, 8,flip)
          popsSCP_d(m,i-1,j+1,k  , 9, flop) = popsSCP_d(m,i,j,k, 9,flip)
          popsSCP_d(m,i+1,j-1,k  ,10, flop) = popsSCP_d(m,i,j,k,10,flip)
          popsSCP_d(m,i+1,j  ,k+1,11, flop) = popsSCP_d(m,i,j,k,11,flip)
          popsSCP_d(m,i-1,j  ,k-1,12, flop) = popsSCP_d(m,i,j,k,12,flip)
          popsSCP_d(m,i-1,j  ,k+1,13, flop) = popsSCP_d(m,i,j,k,13,flip)
          popsSCP_d(m,i+1,j  ,k-1,14, flop) = popsSCP_d(m,i,j,k,14,flip)
          popsSCP_d(m,i  ,j+1,k+1,15, flop) = popsSCP_d(m,i,j,k,15,flip)
          popsSCP_d(m,i  ,j-1,k-1,16, flop) = popsSCP_d(m,i,j,k,16,flip)
          popsSCP_d(m,i  ,j-1,k+1,17, flop) = popsSCP_d(m,i,j,k,17,flip)
          popsSCP_d(m,i  ,j+1,k-1,18, flop) = popsSCP_d(m,i,j,k,18,flip)
#else
          popsSCP_d(m,i  ,j  ,k  , 0, flop) = popsSCP_d(m,i,j,k, 0, flip)
          popsSCP_d(m,i+1,j  ,k  , 1, flop) = popsSCP_d(m,i,j,k, 1,flip)
          popsSCP_d(m,i-1,j  ,k  , 2, flop) = popsSCP_d(m,i,j,k, 2,flip)
          popsSCP_d(m,i  ,j+1,k  , 3, flop) = popsSCP_d(m,i,j,k, 3,flip)
          popsSCP_d(m,i  ,j-1,k  , 4, flop) = popsSCP_d(m,i,j,k, 4,flip)
          popsSCP_d(m,i  ,j  ,k+1, 5, flop) = popsSCP_d(m,i,j,k, 5,flip)
          popsSCP_d(m,i  ,j  ,k-1, 6, flop) = popsSCP_d(m,i,j,k, 6,flip)
#endif
        enddo
      endif
    end subroutine stream_SCP 
    

  end module kernels_scp