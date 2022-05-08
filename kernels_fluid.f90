#include "defines.h"
  module kernels_fluid
    use dimensions_m
    implicit none

    integer(4), constant :: nz_d
    real,constant        :: cz_d

    ! device arrays
    real(4), allocatable, device :: popsR_d(:,:,:,:,:), popsB_d(:,:,:,:,:)
    real(4), allocatable, device :: popsSCP_d(:,:,:,:,:,:)
    real(4), allocatable, device :: force_d(:,:,:,:)
    real(4), allocatable, device :: rhoR_d(:,:,:),rhoB_d(:,:,:)
    integer(1), allocatable, device :: nearsel_d(:,:,:)
    real(4), allocatable, device :: vel_d(:,:,:,:)
    
    real(4), allocatable, device :: scalar_d(:,:,:,:)

    real(4), allocatable, device :: pos_d(:,:,:)
    real(4), allocatable, device :: forceAtoms_d(:,:,:)
    real(8), allocatable, device :: myf_d(:,:), myt_d(:,:)
    real(8), allocatable, device :: myf2_d(:,:), myt2_d(:,:)
    
    real(4), allocatable, device :: force_buf_up_d(:,:,:,:)
    real(4), allocatable, device :: force_buf_down_d(:,:,:,:)
    
    ! integer(2), allocatable, device :: listAtoms_d(:)
    ! integer(2), allocatable, device :: findAtoms_d(:,:,:,:)

    integer, allocatable, device :: matrixAtoms_d(:,:,:,:)
    integer, allocatable, device :: dim_matAtm_d(:,:,:)
    integer, parameter :: maxAtomInBlock = 1000

    integer(1), allocatable, device :: issub_d(:, :, :)
    integer, device :: rmax_issub_d, sphereMax_d
    integer(1), allocatable, device :: myfluid_d(:,:,:, :)
    
    integer, device :: ndouble
    integer, device, allocatable :: exdouble(:), eydouble(:), ezdouble(:)

    integer, device :: stop_d, countmk_d, countrm_d
    integer, device, allocatable :: partVol_d(:), oldpartVol_d(:)
    integer, device :: countn2p_d, countp2n1_d,countp2n2_d

    real, device :: minRho = -10.0, maxRho = 10.0
    real, device :: max_vel_d = 0.0
    real, device :: minSCP = -10.0, maxSCP = 10.0

    ! Debug
    real, allocatable, device :: debugline_d(:,:,:)
    real, device :: debugn2pf_d(20000,5)
    real, device :: debugn2pt_d(20000,5)

    real, device :: debugp2n1_d(20000,4)
    real, device :: debugp2n2_d(20000,4)

    real, device :: debug_rm_d(20000,5)
    real, device :: debug_mk_d(20000,5)
    integer(2), allocatable, device :: debugfluid_d(:,:,:)
    public :: compute_rho
  contains
  
    attributes(device) function compute_rho(pops, i,j,k, flip)
     integer, intent(in) :: i,j,k, flip
     real, intent(in) :: pops(0:nx+1,0:ny+1,0:nz_d+1, 0:npops-1, 2)
     real :: compute_rho
#ifdef D3Q27
     compute_rho = pops(i,j,k,0, flip) + pops(i,j,k,1, flip) + pops(i,j,k,2, flip) + &
       pops(i,j,k,3, flip) + pops(i,j,k,4, flip) + &
       pops(i,j,k,5, flip) + pops(i,j,k,6, flip) + pops(i,j,k,7, flip) + &
       pops(i,j,k,8, flip) + pops(i,j,k,9, flip) + &
       pops(i,j,k,10, flip) + pops(i,j,k,11, flip) + pops(i,j,k,12, flip) + &
       pops(i,j,k,13, flip) + pops(i,j,k,14, flip) + &
       pops(i,j,k,15, flip) + pops(i,j,k,16, flip) + &
       pops(i,j,k,17, flip) + pops(i,j,k,18, flip) + &
       pops(i,j,k,19, flip) + pops(i,j,k,20, flip) + &
       pops(i,j,k,21, flip) + pops(i,j,k,22, flip) + &
       pops(i,j,k,23, flip) + pops(i,j,k,24, flip) + &
       pops(i,j,k,25, flip) + pops(i,j,k,26, flip)
#else
     compute_rho = pops(i,j,k,0, flip) + pops(i,j,k,1, flip) + pops(i,j,k,2, flip) + &
       pops(i,j,k,3, flip) + pops(i,j,k,4, flip) + &
       pops(i,j,k,5, flip) + pops(i,j,k,6, flip) + pops(i,j,k,7, flip) + &
       pops(i,j,k,8, flip) + pops(i,j,k,9, flip) + &
       pops(i,j,k,10, flip) + pops(i,j,k,11, flip) + pops(i,j,k,12, flip) + &
       pops(i,j,k,13, flip) + pops(i,j,k,14, flip) + &
       pops(i,j,k,15, flip) + pops(i,j,k,16, flip) + &
       pops(i,j,k,17, flip) + pops(i,j,k,18, flip)
#endif       
    end function compute_rho  
    
    attributes(global) subroutine copy_isfluid(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, flop
      real    :: rhoR

      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x - nbuff
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y - nbuff
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z - nbuff
      
      if (i>nx+nbuff) return
      if (j>ny+nbuff) return
      if (k>nz_d+nbuff) return
      
      flop = 3 - flip
      
      myfluid_d(i,j,k, flop) = myfluid_d(i,j,k, flip)


    end subroutine copy_isfluid
    
    
!    attributes(device) function compute_u_1fl(pops, i,j,k, flip, invrho)
!     integer, intent(in) :: i,j,k, flip
!     real, intent(in) :: pops(0:nx+1,0:ny+1,0:nz_d+1, 0:npops-1, 2), invrho
!     real :: compute_u_1fl
!
!     compute_u_1fl   = invrho * ( pops(i,j,k,1, flip) - pops(i,j,k,2, flip) + pops(i,j,k,7, flip) - &
!       pops(i,j,k,8, flip) - pops(i,j,k,9, flip) + &
!       pops(i,j,k,10, flip) + pops(i,j,k,11, flip) - pops(i,j,k,12, flip) - &
!       pops(i,j,k,13, flip) + pops(i,j,k,14, flip) )
!    end function compute_u_1fl

!    attributes(device) function compute_v_1fl(pops, i,j,k, flip, invrho)
!     integer, intent(in) :: i,j,k, flip
!     real, intent(in) :: pops(0:nx+1,0:ny+1,0:nz_d+1, 0:npops-1, 2), invrho
!     real :: compute_v_1fl
!
!     compute_v_1fl    = invrho * ( pops(i,j,k,3, flip) - pops(i,j,k,4, flip) + pops(i,j,k,7, flip) - &
!       pops(i,j,k,8, flip) + pops(i,j,k,9, flip) - &
!       pops(i,j,k,10, flip) + pops(i,j,k,15, flip) - pops(i,j,k,16, flip) - &
!       pops(i,j,k,17, flip) + pops(i,j,k,18, flip) )
!    end function compute_v_1fl

!    attributes(device) function compute_w_1fl(pops, i,j,k, flip, invrho)
!     integer, intent(in) :: i,j,k, flip
!     real, intent(in) :: pops(0:nx+1,0:ny+1,0:nz_d+1, 0:npops-1, 2), invrho
!     real :: compute_w_1fl
!
!     compute_w_1fl    = invrho * ( pops(i,j,k,5, flip) - pops(i,j,k,6, flip) + pops(i,j,k,11, flip) - &
!       pops(i,j,k,12, flip) + pops(i,j,k,13, flip) - &
!       pops(i,j,k,14, flip) + pops(i,j,k,15, flip) - pops(i,j,k,16, flip) + &
!       pops(i,j,k,17, flip) - pops(i,j,k,18, flip) )
!
!    end function compute_w_1fl
    
  end module