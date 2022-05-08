#include "defines.h"

  module dimensions_m
    implicit none

    integer(4), parameter :: glx = MYDIMESION, gly = glx, glz = glx
    ! integer(4), parameter :: glx = 128, gly = glx, glz = glx
    integer(4), parameter :: nx = glx, ny = gly
    !integer(4) :: nz
#ifdef NEARCONTACT
    integer, parameter :: nbuff=4
#else
    integer, parameter :: nbuff=2
#endif
    logical, parameter :: lprintraw=.true.
    logical, parameter :: lprintvel=.true.
    logical, parameter :: lprintrhoB=.false.
    logical, parameter :: lprintphase=.false.
    logical, parameter :: lpoptransf=.false.
    logical, parameter :: lcompute_totrho=.true.
    
    logical, save :: diagnostic=.true.
    integer, save :: tdiagnostic=1
    real, save :: totrho(1)
    
!
! mpi stuff
    integer:: nprocs, myrank, ierr
    integer:: xyplane, xyplane_int1, xyplane3, xyplanebuff, xyplanesvar
    integer :: xyplanescp
    integer, parameter::  mpid=3
    integer:: front(2),rear(2),right(2)
    integer:: up(2),down(2),left(2)
    logical remdims(mpid)
    logical Tperiodic(mpid)
    logical rreorder
    integer:: lbecomm,localcomm
    integer:: mydev 
    integer:: prgrid(mpid)
    integer:: mpicoords(mpid)
    integer:: gsizes(3),lsizes(3),start_idx(3)
    integer:: offset(3)
    integer, save :: numbc=0
    integer, allocatable :: bctype(:,:)
    real, allocatable :: bcvel(:,:)
    real, allocatable :: bcrho(:,:)
#ifdef FIXDENSWALL    
    logical, parameter :: fixdenswall=.true.
#else
    logical, parameter :: fixdenswall=.false.
#endif    
    integer, allocatable :: bcscptype(:)
    real, allocatable :: bcscp(:,:)
    
!
#ifdef PBC
    logical, parameter :: xperiodic = .true.
    logical, parameter :: yperiodic = .true.
    logical, parameter :: zperiodic = .true.
#else
    logical, parameter :: xperiodic = .false.
    logical, parameter :: yperiodic = .false.
    logical, parameter :: zperiodic = .false.
#endif
    logical, parameter :: uniqueOmega = .true.
    logical, parameter :: wantOut = .true.
    logical, parameter :: textVTK = .false.
    logical, save :: forced = .false.
    logical, save :: const_forced = .false.
    logical, save :: store_vel = .true.
#ifdef APPLYBC
    logical, save :: lreadisfluid=.true.
    logical, save :: lwriteisf=.true.
#else
    logical, save :: lreadisfluid=.false.
    logical, save :: lwriteisf=.false.
#endif
    logical, save :: lwriteisf_every=.false.
    logical, save :: lbcforce=.false.
    real, dimension(3), save :: f_cost = (/ 0.0, 0.0, 0.0 /)
    real, save    :: A_rep=0.0

    ! Moving boundaries..
    real, dimension(3), parameter :: wall_x0 = (/ 0.0, 0.0, 0.0 /), wall_x1 = (/ 0.0, 0.0, 0.0 /)
    real, dimension(3), parameter :: wall_y0 = (/ 0.0, 0.0, 0.0 /), wall_y1 = (/ 0.0, 0.0, 0.0 /)
    real, dimension(3), parameter :: wall_z0 = (/ 0.0, 0.0, 0.0 /), wall_z1 = (/ 0.0, 0.0, 0.0 /)
    ! real, dimension(3), parameter :: wall_z0 = (/ -0.1, 0.0, 0.0 /), wall_z1 = (/ 0.1, 0.0, 0.0 /)

    ! CG
#ifdef CG    
    logical, parameter :: withCG = .true.
#else
    logical, parameter :: withCG = .false.
#endif

    real, save :: sigma_CG = 0.1, beta_CG = 0.999

    ! Particles
    logical :: withParticles = .false.
    logical :: lrotate = .false.

    integer :: numAtoms = 1

    logical, parameter :: atomsFromFile = .true.
    real, parameter    :: rdim  = 5.5
    real, parameter    :: rmass = 1.0 / 1500.0   ! 1/mass for particles

    logical, parameter :: lfix_moment = .false.
    real               :: max_vel = 0.0
    
    real    :: densR,densB
    real    :: denswallR,denswallB
    
    ! External particle forces
    real    :: ext_fxx, ext_fyy, ext_fzz
    real    :: ext_tqx, ext_tqy, ext_tqz

    ! Side Wall
    logical, parameter :: withSidewall = .false.
    real, parameter    :: sidewall_rdist = 5.5
    real, parameter    :: sidewall_k = 4

    ! Hertzian
    logical, parameter :: withHz = .true.
    real, parameter    :: hz_k = 10
    real, parameter    :: hz_rmin = 10.5
    real, parameter    :: hz_rcap = 9.0
    logical, parameter :: llubrication = .false.
    real, parameter    :: lubric_k = 10
    real, parameter    :: lubric_rmin = 10.0
    real, parameter    :: lubric_rcap = 9.0
    
    ! Precision dependant
    real, parameter :: ZERO = 0.0
    real, parameter :: ONE  = 1.0
    real, parameter :: HALF = 0.5
    real, parameter :: TWO  = 2.0
    real, parameter :: THREE = 3.0
    real, parameter :: FOUR = 4.0
    real, parameter :: FIVE = 5.0
    real, parameter :: SIX = 6.0
    real, parameter :: EIGHT = 8.0
    real, parameter :: NINE = 9.0
    real, parameter :: TEN = 10.0
    real, parameter :: TWELVE = 12.0
    real, parameter :: TWENTYFOUR = 24.0
    real, parameter :: TWENTYSEVEN = 27.0
    real, parameter :: FIFTYFOUR = 54.0
    
    real, parameter :: cssq = 1.0 / 3.0
    real, parameter :: pref_bouzidi = TWO / cssq
    real, parameter :: Pi = real(3.141592653589793238462643383279502884d0,kind=4)
    
    real :: tauR = 1.0
    real :: tauB = 1.0
    real :: omega = 1.0
    
    ! From tau
    real :: viscR = cssq*(ONE-HALF)
    real :: viscB = cssq*(ONE-HALF)

    ! SCP
    integer, parameter :: numscp=9
#ifdef SCP    
    logical, parameter :: withSCP = .true.
#else
    logical, parameter :: withSCP = .false.
#endif

    !!!!!!!!!!!!!!!!!!!!! END OF DIRECT PARAMETERS !!!!!!!!!!!!!!!!!!!!!

#ifdef D3Q27
    
    ! Q3D17 related
    integer(4), parameter :: npops = 27
    integer(8), parameter :: lup = (glx+2)*(gly+2)*(glz+2)
    integer(8), parameter :: dataSz = 4 * lup * npops
    character(len=*), parameter :: lattice='D3Q27'
    real, parameter :: zeta = 9.0 / 19.0
    real, parameter :: p0 = ( TWO / THREE )**THREE
    real, parameter :: p1 = ( TWO / THREE )**TWO * (ONE/SIX)
    real, parameter :: p2 = ( ONE / SIX )**TWO * (TWO/THREE)
    real, parameter :: p3 = ( ONE / SIX )**THREE
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
  	
  	! CG Shimpei Saito
    real, parameter :: b0 = -(TEN/TWENTYSEVEN)
    real, parameter :: b1 = TWO/TWENTYSEVEN
    real, parameter :: b2 = ONE/FIFTYFOUR
    real, parameter :: b3 = ONE/216.0
    real, dimension(0:npops-1), parameter, public :: &
      b_l = (/b0,b1,b1,b1,b1,b1,b1,b2,b2,b2,b2,b2,b2,b2,b2,b2,b2,b2,b2, &
       b3,b3,b3,b3,b3,b3,b3,b3/) 
    
    real, parameter :: phi0 = ZERO
    real, parameter :: phi1 = TWO/19.0
    real, parameter :: phi2 = ONE/38.0
    real, parameter :: phi3 = ONE/152.0
    
    real, dimension(0:npops-1), parameter, public :: phi = &
     (/phi0,phi1,phi1,phi1,phi1,phi1,phi1, &
      phi2,phi2,phi2,phi2,phi2,phi2,phi2,phi2,phi2,phi2,phi2,phi2, &
      phi3,phi3,phi3,phi3,phi3,phi3,phi3,phi3 /) 
    
    real, parameter :: varphi0 = ONE
    real, parameter :: varphi1 = -TWO/19.0
    real, parameter :: varphi2 = -ONE/38.0
    real, parameter :: varphi3 = -ONE/152.0
    
    real, dimension(0:npops-1), parameter, public :: varphi = &
     (/varphi0,varphi1,varphi1,varphi1,varphi1,varphi1,varphi1, &
      varphi2,varphi2,varphi2,varphi2,varphi2,varphi2, &
      varphi2,varphi2,varphi2,varphi2,varphi2,varphi2, &
      varphi3,varphi3,varphi3,varphi3,varphi3,varphi3,varphi3,varphi3 /) 
    
    real, parameter :: psi0 = -20.0/NINE
    real, parameter :: psi1 = -TWO/NINE
    real, parameter :: psi2 = ONE/36.0
    real, parameter :: psi3 = ONE/36.0
    
    real, dimension(0:npops-1), parameter, public :: psi = &
     (/psi0,psi1,psi1,psi1,psi1,psi1,psi1, &
      psi2,psi2,psi2,psi2,psi2,psi2,psi2,psi2,psi2,psi2,psi2,psi2, &
      psi3,psi3,psi3,psi3,psi3,psi3,psi3,psi3 /) 
    
    real, parameter :: xi0 = ZERO
    real, parameter :: xi1 = ONE/THREE
    real, parameter :: xi2 = ONE/TWELVE
    real, parameter :: xi3 = ONE/48.0
    
    real, dimension(0:npops-1), parameter, public :: xi = &
     (/xi0,xi1,xi1,xi1,xi1,xi1,xi1, &
      xi2,xi2,xi2,xi2,xi2,xi2,xi2,xi2,xi2,xi2,xi2,xi2, &
      xi3,xi3,xi3,xi3,xi3,xi3,xi3,xi3 /) 
  	
  	real, dimension(0:npops-1), parameter, public :: &
      a = (/ ZERO,p1/cssq,p1/cssq,p1/cssq,p1/cssq,p1/cssq,p1/cssq,p2/cssq, &
      p2/cssq,p2/cssq,p2/cssq,p2/cssq,p2/cssq,p2/cssq,p2/cssq,p2/cssq, &
      p2/cssq,p2/cssq,p2/cssq,p3/cssq,p3/cssq,p3/cssq,p3/cssq, &
      p3/cssq,p3/cssq,p3/cssq,p3/cssq /)
    
#else

    ! Q3D19 related
    integer(4), parameter :: npops = 19
    integer(8), parameter :: lup = (glx+2)*(gly+2)*(glz+2)
    integer(8), parameter :: dataSz = 4 * lup * npops
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
  	
  	! CG
    real, parameter :: b_xi = 0.5 !Leclaire
    real, parameter :: b0 = -(TWO+TWO*b_xi)/(THREE*b_xi+TWELVE)
    real, parameter :: b1 = (b_xi)/(SIX*b_xi+TWENTYFOUR)
    real, parameter :: b2 = ONE/(SIX*b_xi+TWENTYFOUR)
    real, dimension(0:npops-1), parameter, public :: &
      b_l = (/b0,b1,b1,b1,b1,b1,b1,b2,b2,b2,b2,b2,b2,b2,b2,b2,b2,b2,b2/) 
    
    real, parameter :: phi0 = ZERO
    real, parameter :: phi1 = ONE/TWELVE
    real, parameter :: phi2 = ONE/TWENTYFOUR
    
    real, dimension(0:npops-1), parameter, public :: phi = &
     (/phi0,phi1,phi1,phi1,phi1,phi1,phi1, &
      phi2,phi2,phi2,phi2,phi2,phi2,phi2,phi2,phi2,phi2,phi2,phi2/) 
    
    real, parameter :: varphi0 = ONE
    real, parameter :: varphi1 = -ONE/TWELVE
    real, parameter :: varphi2 = -ONE/TWENTYFOUR
    
    real, dimension(0:npops-1), parameter, public :: varphi = &
     (/varphi0,varphi1,varphi1,varphi1,varphi1,varphi1,varphi1, &
      varphi2,varphi2,varphi2,varphi2,varphi2,varphi2, &
      varphi2,varphi2,varphi2,varphi2,varphi2,varphi2/) 
    
    real, parameter :: psi0 = -FIVE/TWO
    real, parameter :: psi1 = -ONE/SIX
    real, parameter :: psi2 = ONE/TWENTYFOUR
    
    real, dimension(0:npops-1), parameter, public :: psi = &
     (/psi0,psi1,psi1,psi1,psi1,psi1,psi1, &
      psi2,psi2,psi2,psi2,psi2,psi2,psi2,psi2,psi2,psi2,psi2,psi2/) 
    
    real, parameter :: xi0 = ZERO
    real, parameter :: xi1 = ONE/FOUR
    real, parameter :: xi2 = ONE/EIGHT
    
    real, dimension(0:npops-1), parameter, public :: xi = &
     (/xi0,xi1,xi1,xi1,xi1,xi1,xi1, &
      xi2,xi2,xi2,xi2,xi2,xi2,xi2,xi2,xi2,xi2,xi2,xi2/) 
    
    real, dimension(0:npops-1), parameter, public :: &
      a = (/ZERO,p1/cssq,p1/cssq,p1/cssq,p1/cssq,p1/cssq,p1/cssq,p2/cssq, &
      p2/cssq,p2/cssq,p2/cssq,p2/cssq,p2/cssq,p2/cssq,p2/cssq,p2/cssq, &
      p2/cssq,p2/cssq,p2/cssq/)
    
    
#endif

    real, dimension(2) ::  alphaCG

    ! Q3D27 for onebelt
    integer, parameter, public :: linksd3q27 = 27
    integer, dimension(0:linksd3q27-1), parameter, public :: & 
    !           0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26
     exd3q27=(/ 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 1,-1/)       
    integer, dimension(0:linksd3q27-1), parameter, public :: &
     eyd3q27=(/ 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1, 1,-1,-1, 1/)
    integer, dimension(0:linksd3q27-1), parameter, public :: & 
	 ezd3q27=(/ 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1, 1,-1, 1,-1/)
	integer, dimension(0:linksd3q27-1), parameter, public :: & 
	 oppd3q27=(/ 0, 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13,16,15,18,17,20,19,22,21,24,23,26,25/)
    real, parameter :: p0d3q27 = ( TWO / THREE )**THREE
    real, parameter :: p1d3q27 = ( TWO / THREE )**TWO * (ONE/SIX)
    real, parameter :: p2d3q27 = ( ONE / SIX )**TWO * (TWO/THREE)
    real, parameter :: p3d3q27 = ( ONE / SIX )**THREE
    real, dimension(0:linksd3q27-1), parameter, public :: pd3q27 = &
     (/ p0d3q27,p1d3q27,p1d3q27,p1d3q27,p1d3q27,p1d3q27,p1d3q27, &
      p2d3q27,p2d3q27,p2d3q27,p2d3q27,p2d3q27,p2d3q27,p2d3q27,p2d3q27, &
      p2d3q27,p2d3q27,p2d3q27,p2d3q27,p3d3q27,p3d3q27,p3d3q27,p3d3q27, &
      p3d3q27,p3d3q27,p3d3q27,p3d3q27 /)
    real, dimension(0:linksd3q27-1), parameter, public :: ad3q27 = &
     (/ ZERO,p1d3q27/cssq,p1d3q27/cssq,p1d3q27/cssq,p1d3q27/cssq,p1d3q27/cssq,p1d3q27/cssq, &
      p2d3q27/cssq,p2d3q27/cssq,p2d3q27/cssq,p2d3q27/cssq,&
      p2d3q27/cssq,p2d3q27/cssq,p2d3q27/cssq,p2d3q27/cssq,&
      p2d3q27/cssq,p2d3q27/cssq,p2d3q27/cssq,p2d3q27/cssq,&
      p3d3q27/cssq,p3d3q27/cssq,p3d3q27/cssq,p3d3q27/cssq,&
      p3d3q27/cssq,p3d3q27/cssq,p3d3q27/cssq,p3d3q27/cssq /)



    ! isfluid
    integer, parameter :: fluid_wall = 0
    integer, parameter :: fluid_fluid = 1
    integer, parameter :: fluid_spherelist = 2
    integer, parameter :: fluid_dead = 3
    integer, parameter :: fluid_spheredead = 4
    integer, parameter :: fluid_particleCM = 5

  
    ! center PBC
    real, parameter :: cx=real(nx+1) * HALF
    real, parameter :: cy=real(ny+1) * HALF
    real, parameter :: cz=real(glz+1) * HALF

    ! Device constants...
    logical, constant :: store_vel_d
    integer, constant :: numbc_d
    integer, constant :: numscp_d
    integer, constant :: nbuff_d
    integer, dimension(0:npops-1), constant :: ex_d,ey_d,ez_d,opp_d
    integer, dimension(0:linksd3q27-1), constant :: exd3q27_d,eyd3q27_d,ezd3q27_d
    real, constant :: fx_d,fy_d,fz_d, buz_d
    real, constant :: viscR_d,viscB_d
    real, constant :: densR_d,densB_d
    real, constant :: denswallR_d,denswallB_d
    real, dimension(3), constant :: wall_x0_d, wall_x1_d
    real, dimension(3), constant :: wall_y0_d, wall_y1_d
    real, dimension(3), constant :: wall_z0_d, wall_z1_d
    real, dimension(2), constant ::  alphaCG_d
    real, dimension(0:npops-1), constant :: p_d, a_d, b_l_d
    real, dimension(0:npops-1), constant :: phi_d, varphi_d, psi_d, xi_d
    real, dimension(0:linksd3q27-1), constant :: ad3q27_d
    real, constant :: sigma_CG_d, beta_CG_d
    real, dimension(0:npops-1), constant :: rec_fact_d
    
    real, dimension(1:3,1:3,0:npops-1), constant :: cmat_d
    
    real, dimension(3), constant :: f_cost_d 
    real, constant    :: A_rep_d

    real, constant    :: ext_fxx_d, ext_fyy_d, ext_fzz_d
    real, constant    :: ext_tqx_d, ext_tqy_d, ext_tqz_d

    integer, parameter :: TILE_DIM = 16
    integer, parameter :: TILE_DIMx = TILE1
    integer, parameter :: TILE_DIMy = TILE2
    integer, parameter :: TILE_DIMz = TILE3
    integer, parameter :: TILE_DIMPART = 32
    
    integer, parameter :: NHALO = 2

    ! MPI
    integer, constant  :: offset_d(3)
    
    !BC
    integer, allocatable, device :: bctype_d(:,:)
    real, allocatable, device :: bcvel_d(:,:)
    real, allocatable, device :: bcrho_d(:,:)
    real, device :: totrho_d(1)
    real, device :: totrhobuff_d(1)
    integer, constant :: niterVTK_d
    
    integer, allocatable, device :: bcscptype_d(:)
    real, allocatable, device :: bcscp_d(:,:)

    ! Particles
    real, parameter    :: atmWeight = 1.0 / rmass
    real, parameter    :: rotinx = TWO/FIVE * atmWeight * rdim*rdim
    real, parameter    :: rotiny = TWO/FIVE * atmWeight * rdim*rdim
    real, parameter    :: rotinz = TWO/FIVE * atmWeight * rdim*rdim
	
    ! PARTICLES 	-----------------------------------------------
    ! pos(1:3) -> xxx,yyy,zzz
    ! pos(4:6) -> vxx,vyy,vzz
    ! pos(7:9) -> oxx,oyy,ozz
    ! pos(10:13) -> q0,q1,q2,q3
    ! force(1:3) -> fxx,fyy,fzz
    ! force(4:6) -> fxb,fyb,fzb
    ! force(7:9) -> tqx,tqy,tqz
    ! force(10:12) -> tqb,tqb,tqb
    
    
    !SCP
    
#ifdef SCPD3Q27
    
    ! Q3D27 related
    real, parameter :: cssq_scp = 1.0 / 3.0
    integer(4), parameter :: npops_scp = 27
    integer(8), parameter :: lup_scp = (glx+2)*(gly+2)*(glz+2)
    integer(8), parameter :: dataSz_scp = 4 * lup_scp * npops_scp
    character(len=*), parameter :: lattice_scp='D3Q27'
    real, parameter :: p0_scp = ( TWO / THREE )**THREE
    real, parameter :: p1_scp = ( TWO / THREE )**TWO * (ONE/SIX)
    real, parameter :: p2_scp = ( ONE / SIX )**TWO * (TWO/THREE)
    real, parameter :: p3_scp = ( ONE / SIX )**THREE
    real, dimension(0:npops_scp-1), parameter, public :: p_scp = &
     (/ p0_scp,p1_scp,p1_scp,p1_scp,p1_scp,p1_scp,p1_scp,p2_scp,p2_scp,p2_scp,p2_scp,p2_scp,p2_scp,p2_scp,p2_scp, &
      p2_scp,p2_scp,p2_scp,p2_scp,p3_scp,p3_scp,p3_scp,p3_scp,p3_scp,p3_scp,p3_scp,p3_scp /)
    !lattice vectors
    integer, dimension(0:npops_scp-1), parameter :: &
   	!       0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26
  	ex_scp = (/ 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 1,-1/)
    integer, dimension(0:npops_scp-1), parameter:: &
  	ey_scp = (/ 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1, 1,-1,-1, 1/)
    integer, dimension(0:npops_scp-1), parameter:: &
  	ez_scp = (/ 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1, 1,-1, 1,-1/)
    integer, dimension(0:npops_scp-1), parameter:: &
  	opp_scp =(/ 0, 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13,16,15,18,17,20,19,22,21,24,23,26,25/)
  	
  	

    
#elif defined SCPD3Q19
    
    ! Q3D19 related
    real, parameter :: cssq_scp = 1.0 / 3.0
    integer(4), parameter :: npops_scp = 19
    integer(8), parameter :: lup_scp = (glx+2)*(gly+2)*(glz+2)
    integer(8), parameter :: dataSz_scp = 4 * lup_scp * npops_scp
    character(len=*), parameter :: lattice_scp='D3Q19'
    real, parameter :: p0_scp = 1.0 / 3.0
    real, parameter :: p1_scp = 1.0 / 18.0
    real, parameter :: p2_scp = 1.0 / 36.0
    real, dimension(0:npops_scp-1), parameter :: p_scp = (/p0_scp, &
     p1_scp,p1_scp,p1_scp,p1_scp,p1_scp,p1_scp,p2_scp,p2_scp,p2_scp,p2_scp, &
     p2_scp,p2_scp,p2_scp,p2_scp,p2_scp,p2_scp,p2_scp,p2_scp/)
    
    !lattice vectors
    integer, dimension(0:npops_scp-1), parameter :: &
   	!       0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18
  	ex_scp = (/ 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1, 0, 0, 0, 0/)
    integer, dimension(0:npops_scp-1), parameter:: &
  	ey_scp = (/ 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1/)
    integer, dimension(0:npops_scp-1), parameter:: &
  	ez_scp = (/ 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1/)
    integer, dimension(0:npops_scp-1), parameter:: &
  	opp_scp =(/ 0, 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13,16,15,18,17/)
  	

#else

    ! Q3D7 related
    real, parameter :: cssq_scp = 1.0 / 3.0
    integer(4), parameter :: npops_scp = 7
    integer(8), parameter :: lup_scp = (glx+2)*(gly+2)*(glz+2)
    integer(8), parameter :: dataSz_scp = 4 * lup_scp * npops_scp
    character(len=*), parameter :: lattice_scp='D3Q7'
    real, parameter :: p0_scp = 1.0 - 3.0 * cssq_scp
    real, parameter :: p1_scp = 0.5 * cssq_scp
    real, dimension(0:npops_scp-1), parameter :: p_scp = (/p0_scp, &
     p1_scp,p1_scp,p1_scp,p1_scp,p1_scp,p1_scp/)
    
    !lattice vectors
    integer, dimension(0:npops_scp-1), parameter :: &
   	!       0, 1, 2, 3, 4, 5, 6, 
  	ex_scp = (/ 0, 1,-1, 0, 0, 0, 0 /)
    integer, dimension(0:npops_scp-1), parameter:: &
  	ey_scp = (/ 0, 0, 0, 1,-1, 0, 0 /)
    integer, dimension(0:npops_scp-1), parameter:: &
  	ez_scp = (/ 0, 0, 0, 0, 0, 1,-1 /)
    integer, dimension(0:npops_scp-1), parameter:: &
  	opp_scp =(/ 0, 2, 1, 4, 3, 6, 5 /)
    
#endif
     
#ifdef SCPD3Q27GRAD      
      integer(4), parameter :: npops_grd = 27
      real, parameter :: cssq_grd = 1.0 / 3.0
      real, parameter :: p0_grd = ( TWO / THREE )**THREE
      real, parameter :: p1_grd = ( TWO / THREE )**TWO * (ONE/SIX)
      real, parameter :: p2_grd = ( ONE / SIX )**TWO * (TWO/THREE)
      real, parameter :: p3_grd = ( ONE / SIX )**THREE
      real, dimension(0:npops_grd-1), parameter, public :: &
      a_scp = (/ ZERO,p1_grd/cssq_grd,p1_grd/cssq_grd,p1_grd/cssq_grd,p1_grd/cssq_grd,p1_grd/cssq_grd,p1_grd/cssq_grd,p2_grd/cssq_grd, &
      p2_grd/cssq_grd,p2_grd/cssq_grd,p2_grd/cssq_grd,p2_grd/cssq_grd,p2_grd/cssq_grd,p2_grd/cssq_grd,p2_grd/cssq_grd,p2_grd/cssq_grd, &
      p2_grd/cssq_grd,p2_grd/cssq_grd,p2_grd/cssq_grd,p3_grd/cssq_grd,p3_grd/cssq_grd,p3_grd/cssq_grd,p3_grd/cssq_grd, &
      p3_grd/cssq_grd,p3_grd/cssq_grd,p3_grd/cssq_grd,p3_grd/cssq_grd /)
      
      !lattice vectors
    integer, dimension(0:npops_grd-1), parameter :: &
   	!       0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26
  	ex_grd = (/ 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 1,-1/)
    integer, dimension(0:npops_scp-1), parameter:: &
  	ey_grd = (/ 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1, 1,-1,-1, 1/)
    integer, dimension(0:npops_scp-1), parameter:: &
  	ez_grd = (/ 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1, 1,-1, 1,-1/)
    integer, dimension(0:npops_scp-1), parameter:: &
  	opp_grd =(/ 0, 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13,16,15,18,17,20,19,22,21,24,23,26,25/)
      
#elif defined SCPD3Q19GRAD       
      integer(4), parameter :: npops_grd = 19
      real, parameter :: cssq_grd = 1.0 / 3.0
      real, parameter :: p0_grd = 1.0 / 3.0
      real, parameter :: p1_grd = 1.0 / 18.0
      real, parameter :: p2_grd = 1.0 / 36.0 
      real, dimension(0:npops_grd-1), parameter, public :: &
      a_scp = (/ZERO,p1_grd/cssq_grd,p1_grd/cssq_grd,p1_grd/cssq_grd,p1_grd/cssq_grd,p1_grd/cssq_grd,p1_grd/cssq_grd,p2_grd/cssq_grd, &
      p2_grd/cssq_grd,p2_grd/cssq_grd,p2_grd/cssq_grd,p2_grd/cssq_grd,p2_grd/cssq_grd,p2_grd/cssq_grd,p2_grd/cssq_grd,p2_grd/cssq_grd, &
      p2_grd/cssq_grd,p2_grd/cssq_grd,p2_grd/cssq_grd /)
      
      !lattice vectors
    integer, dimension(0:npops_grd-1), parameter :: &
   	!       0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18
  	ex_grd = (/ 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1, 0, 0, 0, 0/)
    integer, dimension(0:npops_grd-1), parameter:: &
  	ey_grd = (/ 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1/)
    integer, dimension(0:npops_grd-1), parameter:: &
  	ez_grd = (/ 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1/)
    integer, dimension(0:npops_grd-1), parameter:: &
  	opp_grd =(/ 0, 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13,16,15,18,17/)
      
      
#else
      integer(4), parameter :: npops_grd = 7
      real, parameter :: cssq_grd = 1.0 / 3.0
      real, parameter :: p0_grd = 1.0 - 3.0 * cssq_grd
      real, parameter :: p1_grd = 0.5 * cssq_grd
      real, dimension(0:npops_grd-1), parameter, public :: &
      a_scp = (/ZERO,p1_grd/cssq_grd,p1_grd/cssq_grd,p1_grd/cssq_grd,p1_grd/cssq_grd,p1_grd/cssq_grd,p1_grd/cssq_grd /)
      
      !lattice vectors
    integer, dimension(0:npops_grd-1), parameter :: &
   	!       0, 1, 2, 3, 4, 5, 6, 
  	ex_grd = (/ 0, 1,-1, 0, 0, 0, 0 /)
    integer, dimension(0:npops_grd-1), parameter:: &
  	ey_grd = (/ 0, 0, 0, 1,-1, 0, 0 /)
    integer, dimension(0:npops_grd-1), parameter:: &
  	ez_grd = (/ 0, 0, 0, 0, 0, 1,-1 /)
    integer, dimension(0:npops_grd-1), parameter:: &
  	opp_grd =(/ 0, 2, 1, 4, 3, 6, 5 /)
      
#endif
    
  end module dimensions_m