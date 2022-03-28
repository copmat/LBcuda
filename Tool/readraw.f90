#include "defines.h"
 program generate_isfluid
 
 real(kind=8), parameter :: & 
  Pi=3.141592653589793238462643383279502884d0
 integer :: i,j,k,itemp,ii,jj,kk,ll
 real(8) :: dist(3),cm(3),radius,width,cell(9),rdist
 integer, parameter :: inx=1
 integer, parameter :: iny=1
 integer, parameter :: inz=1
 integer, parameter :: nx=8
 integer, parameter :: ny=8
 integer, parameter :: nz=8
 character(len=*), parameter :: filename = 'output/isfluid.00000000.raw'
 logical :: lexist
 integer(1), dimension(inx:nx,iny:ny,inz:nz) :: isfluid
 logical :: lpbc(3)
 integer :: ndim,izero,itre
 integer(1), parameter :: THREE=int(3,kind=1)
 integer(1), parameter :: ZERO=int(0,kind=1)
 real(4),allocatable :: data1(:,:,:)
 real(4),allocatable :: data3(:,:,:,:)
 
#ifdef D3Q27
    
    ! Q3D17 related
    integer(4), parameter :: npops = 27
    
    
    character(len=*), parameter :: lattice='D3Q27'
    real, parameter :: zeta = 9.0 / 19.0
    real, parameter :: p0 = ( 2.0 / 3.0 )**3.0
    real, parameter :: p1 = ( 2.0 / 3.0 )**2.0 * (1.0/6.0)
    real, parameter :: p2 = ( 1.0 / 6.0 )**2.0 * (2.0/3.0)
    real, parameter :: p3 = ( 1.0 / 6.0 )**3.0
    real, dimension(0:npops-1), parameter, public :: p = &
     (/ p0,p1,p1,p1,p1,p1,p1,p2,p2,p2,p2,p2,p2,p2,p2, &
      p2,p2,p2,p2,p3,p3,p3,p3,p3,p3,p3,p3 /)
    !lattice vectors
    integer, dimension(0:npops-1), parameter :: &
   	!       0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26
  	ex = (/ 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 1,-1/)
    integer, dimension(0:npops-1), parameter:: &
  	ey = (/ 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1, 1,-1,-1, 1/)
    integer, dimension(0:npops-1), parameter:: &
  	ez = (/ 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1, 1,-1, 1,-1/)
    integer, dimension(0:npops-1), parameter:: &
  	opp =(/ 0, 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13,16,15,18,17,20,19,22,21,24,23,26,25/)
  	
#else

    ! Q3D19 related
    integer(4), parameter :: npops = 19
    
   
    character(len=*), parameter :: lattice='D3Q19'
    real, parameter :: zeta = 1.0 / 2.0
    real, parameter :: p0 = 1.0 / 3.0
    real, parameter :: p1 = 1.0 / 18.0
    real, parameter :: p2 = 1.0 / 36.0
    real, dimension(0:npops-1), parameter :: p = (/p0,p1,p1,p1,p1,p1,p1,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2/)
    
    !lattice vectors
    integer, dimension(0:npops-1), parameter :: &
   	!       0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18
  	ex = (/ 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1, 0, 0, 0, 0/)
    integer, dimension(0:npops-1), parameter:: &
  	ey = (/ 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1/)
    integer, dimension(0:npops-1), parameter:: &
  	ez = (/ 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1/)
    integer, dimension(0:npops-1), parameter:: &
  	opp =(/ 0, 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13,16,15,18,17/)
    
#endif
  
 lpbc=.true.
 isfluid=int(3,kind=1)
 cm(1)=dble(nx)*0.5d0
 cm(2)=dble(ny)*0.5d0
 cm(3)=dble(nz)*0.5d0
 radius=12.d0
 
  cell(1)=real(nx,kind=8)
  cell(2)=0.d0
  cell(3)=0.d0
  cell(4)=0.d0
  cell(5)=real(ny,kind=8)
  cell(6)=0.d0
  cell(7)=0.d0
  cell(8)=0.d0
  cell(9)=real(nz,kind=8)
 
 ndim=1
 inquire(file=filename,exist=lexist)
 if(.not. lexist)then
   write(6,*)'file ',filename,' not found!'
   stop
 endif
 open(unit=15,file=filename,access='stream',status='old')
 
 
 if(ndim==1)then
   allocate(data1(nx,ny,nz))
   read(15)data1(1:nx,1:ny,1:nz)
   
   do k=1,nz
     do j=1,ny
       do i=1,nx
         write(6,*)i,j,k,nint(data1(i,j,k))
       enddo
     enddo
   enddo
 endif
 

 
 write(6,*)'Everything is ok'
 
 contains
 
  subroutine pbc_images_onevec(imcons,cells, dists)
  implicit none
  
  integer, intent(in) :: imcons
  real(kind=8), intent(in) :: cells(9)
  real(kind=8) :: dists(3)

  integer :: i
  real(kind=8) aaa,bbb,ccc
  
  select case(imcons)
  case(0)
    return
  case(1)
    aaa=1.d0/cells(1)
      dists(1)=dists(1)-cell(1)*nint(aaa*dists(1))
  case(2)
    bbb=1.d0/cells(5)
      dists(2)=dists(2)-cell(5)*nint(bbb*dists(2))
  case(3)
    aaa=1.d0/cells(1)
    bbb=1.d0/cells(5)
      dists(1)=dists(1)-cell(1)*nint(aaa*dists(1))
      dists(2)=dists(2)-cell(5)*nint(bbb*dists(2))
  case(4)
    ccc=1.d0/cells(9)
      dists(3)=dists(3)-cell(9)*nint(ccc*dists(3))
  case(5)
    aaa=1.d0/cells(1)
    ccc=1.d0/cells(9)
      dists(1)=dists(1)-cell(1)*nint(aaa*dists(1))
      dists(3)=dists(3)-cell(9)*nint(ccc*dists(3))
  case(6)
    bbb=1.d0/cells(5)
    ccc=1.d0/cells(9)
      dists(2)=dists(2)-cell(5)*nint(bbb*dists(2))
      dists(3)=dists(3)-cell(9)*nint(ccc*dists(3))
  case(7)
    aaa=1.d0/cells(1)
    bbb=1.d0/cells(5)
    ccc=1.d0/cells(9)
      dists(1)=dists(1)-cell(1)*nint(aaa*dists(1))
      dists(2)=dists(2)-cell(5)*nint(bbb*dists(2))
      dists(3)=dists(3)-cell(9)*nint(ccc*dists(3))
  end select
  
 end subroutine pbc_images_onevec
 
  pure function fcut(r,inner_cut,outer_cut)

!***********************************************************************
!
!     LBsoft function for fading an observable (r) within a given
!     interval
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification January 2018
!
!***********************************************************************

  implicit none

  real(kind=8), intent(in) :: r,inner_cut,outer_cut
  real(kind=8) :: fcut

  if ( r <= inner_cut ) then
    fcut = 1.d0
  elseif ( r > outer_cut ) then
      fcut = 0.d0
  else
      fcut = 0.5d0*cos((r-inner_cut)*Pi/(outer_cut-inner_cut))+0.5d0
  endif

  return

 end function fcut
 
 pure function pimage(ipbcsub,i,nssub)
 
!***********************************************************************
!     
!     LBsoft sfunction to impose the pbc 
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification November 2018
!     
!***********************************************************************
 
  implicit none
  logical, intent(in) :: ipbcsub
  integer, intent(in) :: i,nssub
  integer :: pimage
  
  pimage=i
  
  if(ipbcsub)then
    if(i<1) then
      pimage=i+nssub
    endif
    if(i>nssub) then
      pimage=i-nssub
    endif
  endif
  
  return
  
 end function pimage

 end program generate_isfluid
