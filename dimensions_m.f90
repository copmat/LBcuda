#include "defines.h"

  module dimensions_m
    implicit none

    ! integer(4), parameter :: nx = 32, ny = nx, nz = nx
    integer(4), parameter :: nx = 128, ny = nx, nz = nx
    logical, parameter :: xperiodic = .true.
    logical, parameter :: yperiodic = .true.
    logical, parameter :: zperiodic = .true.
    logical, parameter :: uniqueOmega = .true.
    logical, parameter :: wantOut = .true.
    logical, parameter :: textVTK = .false.
    logical, parameter :: forced = .false.
    logical, parameter :: const_forced = .false.
    real, dimension(3), parameter :: f_cost = (/ 0.0, 0.0, -0.000981 /)

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

    real, parameter :: sigma_CG = 0.1, beta_CG = 0.999

    ! Particles
#ifdef PARTICLES
    logical, parameter :: withParticles = .true.
#else
    logical, parameter :: withParticles = .false.
#endif    
    
    logical, parameter :: lrotate = .true.
    integer, parameter :: numAtoms = 100
    logical, parameter :: atomsFromFile = .true.
    real, parameter    :: rdim  = 4.5
    real, parameter    :: rmass = 1.0 / 1500.0   ! 1/mass for particles

    logical, parameter :: lfix_moment = .false.
    
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

    !!!!!!!!!!!!!!!!!!!!! END OF DIRECT PARAMETERS !!!!!!!!!!!!!!!!!!!!!

    ! Q3D19 related
    integer(4), parameter :: npops = 19
    integer(8), parameter :: lup = (nx+2)*(ny+2)*(nz+2)
    integer(8), parameter :: dataSz = 4 * lup * npops

    real, parameter :: cssq = 1.0 / 3.0
    real, parameter :: p0 = 1.0 / 3.0
    real, parameter :: p1 = 1.0 / 18.0
    real, parameter :: p2 = 1.0 / 36.0
    real, dimension(0:npops-1), parameter :: p = (/p0,p1,p1,p1,p1,p1,p1,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2/)

    real, parameter :: tauR = 1.0
    real, parameter :: tauB = 1.0


    ! Precision dependant
    real, parameter :: ZERO = 0.0
    real, parameter :: ONE  = 1.0
    real, parameter :: HALF = 0.5
    real, parameter :: TWO  = 2.0
    real, parameter :: THREE = 3.0
    real, parameter :: FOUR = 4.0
    real, parameter :: FIVE = 5.0
    real, parameter :: SIX = 6.0
    real, parameter :: NINE = 9.0
    real, parameter :: TWELVE = 12.0
    real, parameter :: TWENTYFOUR = 24.0
    real, parameter :: pref_bouzidi = TWO / cssq
    real, parameter :: Pi = 3.141592653589793238462643383279502884d0

    ! From tau
    real, parameter :: viscR = cssq*(tauR-HALF)
    real, parameter :: viscB = cssq*(tauB-HALF)

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

    ! Q3D27 for onebelt
    integer, parameter, public :: linksd3q27 = 26
    integer, dimension(0:linksd3q27), parameter, public :: &
      !           0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26
     exd3q27 = (/ 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1 /)
    integer, dimension(0:linksd3q27), parameter, public :: &
     eyd3q27 = (/ 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1, 1,-1 /)
    integer, dimension(0:linksd3q27), parameter, public :: &
     ezd3q27 = (/ 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1 /)

    real, parameter :: p0d3q27 = ( TWO / THREE )**THREE
    real, parameter :: p1d3q27 = ( TWO / THREE )**TWO * (ONE/SIX)
    real, parameter :: p2d3q27 = ( ONE / SIX )**TWO * (TWO/THREE)
    real, parameter :: p3d3q27 = ( ONE / SIX )**THREE
    real, dimension(0:linksd3q27), parameter, public :: pd3q27 = (/ p0d3q27,p1d3q27,p1d3q27,p1d3q27,p1d3q27,p1d3q27,p1d3q27, &
      p2d3q27,p2d3q27,p2d3q27,p2d3q27,p2d3q27,p2d3q27,p2d3q27,p2d3q27, &
      p2d3q27,p2d3q27,p2d3q27,p2d3q27,p3d3q27,p3d3q27,p3d3q27,p3d3q27, &
      p3d3q27,p3d3q27,p3d3q27,p3d3q27 /)

    ! CG
    real, parameter :: b_xi = 0.5 !Leclaire
    real, parameter :: b0 = -(TWO+TWO*b_xi)/(THREE*b_xi+TWELVE)
    real, parameter :: b1 = (b_xi)/(SIX*b_xi+TWENTYFOUR)
    real, parameter :: b2 = ONE/(SIX*b_xi+TWENTYFOUR)
    real, dimension(0:npops-1), parameter, public :: &
      b_l = (/b0,b1,b1,b1,b1,b1,b1,b2,b2,b2,b2,b2,b2,b2,b2,b2,b2,b2,b2/) 

    real, dimension(0:npops-1), parameter, public :: &
      a = (/ZERO,p1/cssq,p1/cssq,p1/cssq,p1/cssq,p1/cssq,p1/cssq,p2/cssq, &
      p2/cssq,p2/cssq,p2/cssq,p2/cssq,p2/cssq,p2/cssq,p2/cssq,p2/cssq, &
      p2/cssq,p2/cssq,p2/cssq/)

    ! isfluid
    integer, parameter :: fluid_wall = 0
    integer, parameter :: fluid_fluid = 1
    integer, parameter :: fluid_spherelist = 2
    integer, parameter :: fluid_spheredead = 4
    integer, parameter :: fluid_particleCM = 5

  
    ! Particles PBC
    real, parameter :: cx=real(nx+1) * HALF
    real, parameter :: cy=real(ny+1) * HALF
    real, parameter :: cz=real(nz+1) * HALF

    ! Device constants...
    integer, dimension(0:npops-1), constant :: ex_d,ey_d,ez_d
    real, constant :: fx_d,fy_d,fz_d, buz_d
    real, dimension(3), constant :: wall_x0_d, wall_x1_d
    real, dimension(3), constant :: wall_y0_d, wall_y1_d
    real, dimension(3), constant :: wall_z0_d, wall_z1_d
    real, dimension(0:npops-1), constant :: p_d, a_d, b_l_d
    real, constant :: sigma_CG_d, beta_CG_d
    real, dimension(0:npops-1), constant :: rec_fact_d

    real, constant    :: ext_fxx_d, ext_fyy_d, ext_fzz_d
    real, constant    :: ext_tqx_d, ext_tqy_d, ext_tqz_d

    integer, parameter :: TILE_DIM = 8
    integer, parameter :: TILE_DIMx = 4
    integer, parameter :: TILE_DIMy = 4
    integer, parameter :: TILE_DIMz = 4
    integer, parameter :: TILE_DIMPART = 16

    ! MPI
    integer idrank, mxrank

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

  end module dimensions_m
