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
  module kernels_fluid_part
    use dimensions_m
    use kernels_fluid
    implicit none
    
    private :: linear
    private :: equil
    
    contains


    attributes(device) pure function linear(i,j,k)
      integer, intent(in) :: i,j,k
      integer :: linear

      linear = i*1000**2 + j*1000 + k
    end function linear


    attributes(device) function compute_u_1fl(pops, i,j,k, flip, invrho)
     integer, intent(in) :: i,j,k, flip
     real, intent(in) :: pops(0:nx+1,0:ny+1,0:nz_d+1, 0:npops-1, 2), invrho
     real :: compute_u_1fl

     compute_u_1fl   = invrho * ( pops(i,j,k,1, flip) - pops(i,j,k,2, flip) + pops(i,j,k,7, flip) - &
       pops(i,j,k,8, flip) - pops(i,j,k,9, flip) + &
       pops(i,j,k,10, flip) + pops(i,j,k,11, flip) - pops(i,j,k,12, flip) - &
       pops(i,j,k,13, flip) + pops(i,j,k,14, flip) )
    end function compute_u_1fl

    attributes(device) function compute_v_1fl(pops, i,j,k, flip, invrho)
     integer, intent(in) :: i,j,k, flip
     real, intent(in) :: pops(0:nx+1,0:ny+1,0:nz_d+1, 0:npops-1, 2), invrho
     real :: compute_v_1fl

     compute_v_1fl    = invrho * ( pops(i,j,k,3, flip) - pops(i,j,k,4, flip) + pops(i,j,k,7, flip) - &
       pops(i,j,k,8, flip) + pops(i,j,k,9, flip) - &
       pops(i,j,k,10, flip) + pops(i,j,k,15, flip) - pops(i,j,k,16, flip) - &
       pops(i,j,k,17, flip) + pops(i,j,k,18, flip) )
    end function compute_v_1fl

    attributes(device) function compute_w_1fl(pops, i,j,k, flip, invrho)
     integer, intent(in) :: i,j,k, flip
     real, intent(in) :: pops(0:nx+1,0:ny+1,0:nz_d+1, 0:npops-1, 2), invrho
     real :: compute_w_1fl

     compute_w_1fl    = invrho * ( pops(i,j,k,5, flip) - pops(i,j,k,6, flip) + pops(i,j,k,11, flip) - &
       pops(i,j,k,12, flip) + pops(i,j,k,13, flip) - &
       pops(i,j,k,14, flip) + pops(i,j,k,15, flip) - pops(i,j,k,16, flip) + &
       pops(i,j,k,17, flip) - pops(i,j,k,18, flip) )
    end function compute_w_1fl

    attributes(device) function equil(rho,u,v,w, l)
     real, intent(in) :: rho,u,v,w
     integer, intent(in) :: l
     real :: equil
     real :: uv

     uv = (1.0/cssq) * (u*ex_d(l) + v*ey_d(l) + w*ez_d(l))
     equil = rho * p_d(l)*(ONE+uv+HALF*(uv*uv)-(HALF/cssq) * (u*u + v*v + w*w))
    end function equil


#ifdef PARTICLES
    !!!!!!!!!!!!!!!!!!! Color gradient !!!!!!!!!!!!!!!!!!!
    attributes(global) subroutine init_isfluid_CG(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l
      integer :: atm_i,atm_j,atm_k, atm_st,atm_en, i_atm,i_list
      logical :: updated
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      if (step==0 .and. 1==i*j*k+offset_d(3)) write(*,*) 'CUDA init_rho_isfluid_CG] i,k:',i,j,k
      
      atm_i = blockIdx%x
      atm_j = blockIdx%y
      atm_k = blockIdx%z
      
      atm_st = 1
      atm_en = dim_matAtm_d(atm_i,atm_j,atm_k)

      ! Update myfluid_d
      updated = .false.
      debugfluid_d(i,j,k) = 0
      if (atm_st>0) then
        do i_list = atm_st, atm_en
          i_atm = matrixAtoms_d(i_list, atm_i,atm_j,atm_k)

          if (.not. updated) then
            updated = setmyfluid(i,j,k, i_atm, step, flip)
          endif
        enddo
      endif
      if (.not. updated) myfluid_d(i,j,k, flip) = fluid_fluid
    
    end subroutine init_isfluid_CG


    attributes(global) subroutine del_fluid_CG(step, flip, lrotate)
      integer, value :: step,flip
      logical, value :: lrotate
      integer :: i,j,k
      integer :: atm_i,atm_j,atm_k, atm_st,atm_en, i_atm,i_list
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z

      atm_i = blockIdx%x 
      atm_j = blockIdx%y 
      atm_k = blockIdx%z 
      ! atm_st = findAtoms_d(1, atm_i,atm_j,atm_k)
      ! atm_en = findAtoms_d(2, atm_i,atm_j,atm_k)
      atm_st = 1
      atm_en = dim_matAtm_d(atm_i,atm_j,atm_k)

      ! del_fluid 1st (old rho), if needed
      if (atm_st>0) then
        do i_list = atm_st, atm_en
          ! i_atm = listAtoms_d(i_list)
          i_atm = matrixAtoms_d(i_list, atm_i,atm_j,atm_k)

          ! Add/delete fluid if particle moved
          if (particle_moved(i_atm) == 1) call del_fluid(step, flip, i_atm,i,j,k, lrotate)
        enddo
      endif
      
    end subroutine del_fluid_CG


    attributes(global) subroutine make_fluid_CG(step, flip, lrotate)
      integer, value :: step,flip
      logical, value :: lrotate
      integer :: i,j,k
      integer :: atm_i,atm_j,atm_k, atm_st,atm_en, i_atm,i_list
  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z

      atm_i = blockIdx%x 
      atm_j = blockIdx%y 
      atm_k = blockIdx%z 
      ! atm_st = findAtoms_d(1, atm_i,atm_j,atm_k)
      ! atm_en = findAtoms_d(2, atm_i,atm_j,atm_k)
      atm_st = 1
      atm_en = dim_matAtm_d(atm_i,atm_j,atm_k)

      ! make_fluid, if needed
      if (atm_st>0) then
        do i_list = atm_st, atm_en
          ! i_atm = listAtoms_d(i_list)
          i_atm = matrixAtoms_d(i_list, atm_i,atm_j,atm_k)

          ! Add/delete fluid if particle moved
          if (particle_moved(i_atm) == 1) call make_fluid(step, flip, i_atm,i,j,k, lrotate)
        enddo
      endif
    end subroutine make_fluid_CG

    attributes(device) pure function getDeltai(i, xxx)
    integer, intent(in) :: i, xxx
    integer :: getDeltai
    
      getDeltai = i - xxx
      if (xperiodic) then
        if (xxx-rdim < 1 .or. xxx+rdim > nx) then
          if (xxx+rdim > nx .and. i < rdim) then
            getDeltai = i - xxx + nx
          else if (xxx-rdim < 1 .and. i > nx - rdim) then
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
          if (yyy+rdim > ny .and. j < rdim) then
            getDeltaj = j - yyy + ny
          else if (yyy-rdim < 1 .and. j > ny - rdim) then
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
        if (zzz-rdim < 1 .or. zzz+rdim > glz) then
          if (zzz+rdim > glz .and. k < rdim) then
            getDeltak = k - zzz + glz
          else if (zzz-rdim < 1 .and. k > glz - rdim) then
            getDeltak = k - zzz - glz
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
      deltak = getDeltak(k+ offset_d(3), zzz)   ! MPI along z
      
      ! Check if inside particle      
      if (abs(deltai)<=rmax_issub_d.and.abs(deltaj)<=rmax_issub_d.and.abs(deltak)<=rmax_issub_d) then
        
        val = issub_d( deltai,deltaj,deltak )
        if (val /= fluid_fluid) then
          setmyfluid = .true.
          myfluid_d(i,j,k, flip) = val
          tempInt = atomicadd( partVol_d(iatm), 1)
          debugfluid_d(i,j,k) = iatm
          ! write(*,*) 'sf1',linear(i,j,k), val+0
        else
          ! write(*,*) 'sf2',linear(i,j,k)
        endif
      else
        ! write(*,*) 'sf3',linear(i,j,k), deltai,deltaj,deltak
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

    attributes(global) subroutine compute_virtual_rho(step, flip)
      integer, value :: step,flip
      integer :: i,j,k, l
      integer :: iatm
      real              :: rtemp(3), ftemp(3),ttemp(3), modr, tempFoo
      real              :: rhoR,rhoB, angleCos
      real              :: aaa,bbb,ccc
      real, dimension(0:3), parameter :: qversor = (/ ZERO , ONE , ZERO , ZERO /)
      real, dimension(0:3) :: qtemp,qtemp2, tempconj

  
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      if (step==1 .and. 1==i*j*k+offset_d(3)) write(*,*) 'CUDA compute_virtual_rho]'
      
      if(myfluid_d(i,j,k, flip)==fluid_spherelist) then
        iatm = debugfluid_d(i,j,k)

        rtemp(1) = i - pos_d(1,iatm,flip)
        rtemp(2) = j - pos_d(2,iatm,flip)
        rtemp(3) = k+ offset_d(3) - pos_d(3,iatm,flip)

        aaa = ONE/real(nx)
        bbb = ONE/real(ny)
        ccc = ONE/real(glz)
        
        rtemp(1) = rtemp(1) - real(nx)*nint(aaa*rtemp(1))
        rtemp(2) = rtemp(2) - real(ny)*nint(bbb*rtemp(2))
        rtemp(3) = rtemp(3) - real(glz)*nint(ccc*rtemp(3))

        modr = sqrt( rtemp(1)*rtemp(1) + rtemp(2)*rtemp(2) + rtemp(3)*rtemp(3))
        rtemp(1:3) = rtemp(1:3)/ modr

        ! Get particle x-versor
        qtemp(0) = pos_d(10,iatm,flip)
        qtemp(1) = pos_d(11,iatm,flip)
        qtemp(2) = pos_d(12,iatm,flip)
        qtemp(3) = pos_d(13,iatm,flip)
        
        tempconj = qconj(qtemp)
        qtemp2 = qtrimult(qtemp,qversor,tempconj)

        ! Calc angle cosine
        angleCos = qtemp2(1)*rtemp(1)+ qtemp2(2)*rtemp(2)+ qtemp2(3)*rtemp(3)
        if (angleCos<0) then
          rhoR_d(i,j,k) = 0.0
          rhoB_d(i,j,k) = 1.0
        else
          rhoR_d(i,j,k) = 1.0
          rhoB_d(i,j,k) = 0.0
        endif
      endif
    end subroutine compute_virtual_rho

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

    
    attributes(device) subroutine del_fluid(step, flip, iatm, i,j,k, lrotate )
      integer, value    :: step, flip
      integer, value    :: iatm, i,j,k
      logical, value    :: lrotate
      logical           :: lmove
      real              :: rtemp(3), ftemp(3),ttemp(3), modr, tempFoo
      real              :: rhoR,rhoB,rhosum,invrho, myu,myv,myw
      real              :: aaa,bbb,ccc
      real(8)           :: acc(3), accFoo1,accFoo2,accFoo3
      integer           :: Dsum, tempInt
      integer :: deltai,deltaj,deltak, atmDist
      integer :: xxx,yyy,zzz

      ! if (1==step) write(*,*) 'CUDA-mk-del-fl]', linear(i,j,k), iatm

      ! Here i,j,k are already periodic-wrapped

      ! We check the distance from old atom position
      xxx = nint(pos_d(1,iatm, 3-flip))
      yyy = nint(pos_d(2,iatm, 3-flip))
      zzz = nint(pos_d(3,iatm, 3-flip))

      deltai = getDeltai(i, xxx)
      deltaj = getDeltaj(j, yyy)
      deltak = getDeltak(k+ offset_d(3), zzz)     ! MPI along z

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
              
        myu   = manual_u(popsR_d) + manual_u(popsB_d) !compute_u_1fl(popsR_d,i,j,k,flip, invrho) + compute_u_1fl(popsB_d,i,j,k,flip, invrho)
        myv   = manual_v(popsR_d) + manual_v(popsB_d) !compute_v_1fl(popsR_d,i,j,k,flip, invrho) + compute_v_1fl(popsB_d,i,j,k,flip, invrho)
        myw   = manual_w(popsR_d) + manual_w(popsB_d) !compute_w_1fl(popsR_d,i,j,k,flip, invrho) + compute_w_1fl(popsB_d,i,j,k,flip, invrho)

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
        tempInt = 1 + atomicadd( countrm_d, 1)
        debug_rm_d (tempInt, 1) = linear(i,j,k)
        debug_rm_d (tempInt, 2) = ftemp(1)
        debug_rm_d (tempInt, 3) = ftemp(2)
        debug_rm_d (tempInt, 4) = ftemp(3)
        debug_rm_d (tempInt, 5) = iatm
        ! write(*,*) 'rm-fl]', linear(i,j,k), ftemp(1),ftemp(2),ftemp(3)
        ! write(*,*) 'myrm]', linear(i,j,k), myu,myv,myw
        ! write(*,*) 'myrm]', linear(i,j,k), rhoR,rhoB
#endif        

        if (lrotate) then
          rtemp(1) = i - pos_d(1,iatm,flip)
          rtemp(2) = j - pos_d(2,iatm,flip)
          rtemp(3) = k+ offset_d(3) - pos_d(3,iatm,flip)

          aaa = ONE/real(nx)
          bbb = ONE/real(ny)
          ccc = ONE/real(glz)
          
          rtemp(1) = rtemp(1) - real(nx)*nint(aaa*rtemp(1))
          rtemp(2) = rtemp(2) - real(ny)*nint(bbb*rtemp(2))
          rtemp(3) = rtemp(3) - real(glz)*nint(ccc*rtemp(3))
          
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
        
        myu   = manual_u(popsR_d) + manual_u(popsB_d) !compute_u_1fl(popsR_d,i,j,k,flip, invrho) + compute_u_1fl(popsB_d,i,j,k,flip, invrho)
        myv   = manual_v(popsR_d) + manual_v(popsB_d) !compute_v_1fl(popsR_d,i,j,k,flip, invrho) + compute_v_1fl(popsB_d,i,j,k,flip, invrho)
        myw   = manual_w(popsR_d) + manual_w(popsB_d) !compute_w_1fl(popsR_d,i,j,k,flip, invrho) + compute_w_1fl(popsB_d,i,j,k,flip, invrho)
        
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
        tempInt = 1 + atomicadd( countrm_d, 1)
        debug_rm_d (tempInt, 1) = linear(i,j,k)
        debug_rm_d (tempInt, 2) = ftemp(1)
        debug_rm_d (tempInt, 3) = ftemp(2)
        debug_rm_d (tempInt, 4) = ftemp(3)
        debug_rm_d (tempInt, 5) = iatm
        ! write(*,*) 'rm-fl2]', linear(i,j,k), ftemp(1),ftemp(2),ftemp(3)
        ! write(*,*) 'myrm2]', linear(i,j,k), myu,myv,myw
        ! write(*,*) 'myrm2]', linear(i,j,k), rhoR,rhoB
#endif

        if (lrotate) then
          rtemp(1) = i - pos_d(1,iatm,flip)
          rtemp(2) = j - pos_d(2,iatm,flip)
          rtemp(3) = k+ offset_d(3) - pos_d(3,iatm,flip)
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

    attributes(device) subroutine make_fluid(step, flip, iatm, i,j,k, lrotate )
      integer, value    :: step, flip
      integer, value    :: iatm, i,j,k
      logical, value    :: lrotate
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

      ! We check the distance from old atom position
      xxx = nint(pos_d(1,iatm, 3-flip))
      yyy = nint(pos_d(2,iatm, 3-flip))
      zzz = nint(pos_d(3,iatm, 3-flip))

      deltai = getDeltai(i, xxx)
      deltaj = getDeltaj(j, yyy)
      deltak = getDeltak(k+ offset_d(3), zzz)       ! MPI over z

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
        tempInt = 1 + atomicadd( countmk_d, 1)
        debug_mk_d (tempInt, 1) = linear(i,j,k)
        debug_mk_d (tempInt, 2) = ftemp(1)
        debug_mk_d (tempInt, 3) = ftemp(2)
        debug_mk_d (tempInt, 4) = ftemp(3)
        debug_mk_d (tempInt, 5) = iatm
        ! write(*,*) 'mk-fl]', linear(i,j,k), ftemp(1),ftemp(2),ftemp(3)
#endif        

        if (lrotate) then
          rtemp(1) = i - pos_d(1,iatm,flip)
          rtemp(2) = j - pos_d(2,iatm,flip)
          rtemp(3) = k+ offset_d(3) - pos_d(3,iatm,flip)

          aaa = ONE/real(nx)
          bbb = ONE/real(ny)
          ccc = ONE/real(glz)
          
          rtemp(1) = rtemp(1) - real(nx)*nint(aaa*rtemp(1))
          rtemp(2) = rtemp(2) - real(ny)*nint(bbb*rtemp(2))
          rtemp(3) = rtemp(3) - real(glz)*nint(ccc*rtemp(3))

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
        tempInt = 1 + atomicadd( countmk_d, 1)
        debug_mk_d (tempInt, 1) = linear(i,j,k)
        debug_mk_d (tempInt, 2) = ftemp(1)
        debug_mk_d (tempInt, 3) = ftemp(2)
        debug_mk_d (tempInt, 4) = ftemp(3)
        debug_mk_d (tempInt, 5) = iatm
        ! write(*,*) 'mk-fl2]', linear(i,j,k), ftemp(1),ftemp(2),ftemp(3)
#endif        

        if (lrotate) then
          rtemp(1) = i - pos_d(1,iatm,flip)
          rtemp(2) = j - pos_d(2,iatm,flip)
          rtemp(3) = k+ offset_d(3) - pos_d(3,iatm,flip)
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
      do l = 1, linksd3q27-1
        i1 = i + exd3q27(l)
        j1 = j + eyd3q27(l)
        k1 = k + ezd3q27(l)
        
        ! CYCLE_OUT_INTERVAL(ishift, minx-nbuff, maxx+nbuff)
        ! CYCLE_OUT_INTERVAL(jshift, miny-nbuff, maxy+nbuff)
        ! CYCLE_OUT_INTERVAL(kshift, minz-nbuff, maxz+nbuff)
        if (i1<0 .or. i1>nx+1 .or. j1<0 .or. j1>ny+1 .or. k1<0 .or. k1>nz_d+1) then
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
        k1 = pimage(k + ezdouble(l), glz)
        if (k1<0 .or. k1>nz_d+1) cycle

        if (i1<0 .or. i1>nx+1 .or. j1<0 .or. j1>ny+1 .or.k1<0 .or. k1>nz_d+1) then
          stop_d = __LINE__
          write(*,*) '2belt] RangeErr i,j,k=',linear(i,j,k),linear(i1,j1,k1)
        endif

        if (myfluid_d(i1,j1,k1, flip)==fluid_fluid .and. &
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
      ! write (*,*) "fix_onebelt_density_twofluids]TO BE IMPL'ED"
    end subroutine fix_onebelt_density_twofluids

    attributes(device) subroutine compute_onebelt_density_2fl_weights(i,j,k, flip, Rsum, Bsum, Dsum)
      implicit none
      integer, intent(in) :: i,j,k, flip
      real, intent(out) :: Rsum,Bsum
      real, intent(out) :: Dsum
      integer :: i1,j1,k1, l
    
      Rsum=ZERO; Bsum=ZERO; Dsum=ZERO
      !compute mean density value
      do l = 1, linksd3q27-1
        i1 = i + exd3q27(l)
        j1 = j + eyd3q27(l)
        k1 = k + ezd3q27(l)

        if (i1<0 .or. i1>nx+1 .or. j1<0 .or. j1>ny+1 .or. k1<0 .or. k1>nz_d+1) then
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
    

    attributes(global) subroutine partBB_CG(step, flip, lrotate)
      integer, value :: step,flip
      logical, value :: lrotate
      integer :: i,j,k, l
      integer :: atm_i,atm_j,atm_k, atm_st,atm_en, i_atm,i_list
      integer :: deltai,deltaj,deltak, atmDist,okDist
      integer :: xxx,yyy,zzz


      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      if (step==1 .and. 1==i*j*k+offset_d(3)) write(*,*) 'CUDA partBB_CG] i,k:',i,j,k

      atm_i = blockIdx%x 
      atm_j = blockIdx%y 
      atm_k = blockIdx%z 
      ! atm_st = findAtoms_d(1, atm_i,atm_j,atm_k)
      ! atm_en = findAtoms_d(2, atm_i,atm_j,atm_k)
      atm_st = 1
      atm_en = dim_matAtm_d(atm_i,atm_j,atm_k)
      if (atm_st>0) then
        do i_list = atm_st, atm_en
          ! i_atm = listAtoms_d(i_list)
          i_atm = matrixAtoms_d(i_list, atm_i,atm_j,atm_k)
          
          if ( myfluid_d(i,j,k, flip)==fluid_spherelist ) then
            xxx = nint(pos_d(1,i_atm, flip))
            yyy = nint(pos_d(2,i_atm, flip))
            zzz = nint(pos_d(3,i_atm, flip))

            deltai = getDeltai(i, xxx)
            deltaj = getDeltaj(j, yyy)
            deltak = getDeltak(k+ offset_d(3), zzz)
            
            atmDist = deltai*deltai + deltaj*deltaj + deltak*deltak
            ! only call n2p_bb if near THIS atom...
            if (atmDist < sphereMax_d) then
#ifdef CHECK_WRONGATOM
              if (debugfluid_d(i,j,k) /= i_atm) then
                xxx = nint(pos_d(1,debugfluid_d(i,j,k), flip))
                yyy = nint(pos_d(2,debugfluid_d(i,j,k), flip))
                zzz = nint(pos_d(3,debugfluid_d(i,j,k), flip))

                deltai = getDeltai(i, xxx)
                deltaj = getDeltaj(j, yyy)
                deltak = getDeltak(k+ offset_d(3), zzz)
                okDist = deltai*deltai + deltaj*deltaj + deltak*deltak
                write(*,*) 'partBB_CG]Wrong atom', linear(i,j,k), i_atm+0.001*atmDist, debugfluid_d(i,j,k)+0.001*okDist
                ! stop_d = __LINE__
              endif
#endif
              call node_to_particle_bounce_back_bc2(step, flip, i_atm, i,j,k, lrotate)
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
      if (step==1 .and. 1==i*j*k+offset_d(3)) write(*,*) 'CUDA partBB_CG_p2n] i,k:',i,j,k

      if ( myfluid_d(i,j,k, flip)==fluid_spherelist ) then
        call particle_to_node_bounce_back_bc2_phase1(step, flip, i,j,k)
      endif
    end subroutine partBB_CG_p2n

    attributes(global) subroutine partBB_CG_p2n2(step, flip, lrotate)
      integer, value :: step,flip
      logical, value :: lrotate
      integer :: i,j,k, l
      integer :: atm_i,atm_j,atm_k, atm_st,atm_en, i_atm,i_list
      integer :: deltai,deltaj,deltak, atmDist
      integer :: xxx,yyy,zzz


      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      if (step==1 .and. 1==i*j*k+offset_d(3)) write(*,*) 'CUDA partBB_CG_p2n2] i,k:',i,j,k

      atm_i = blockIdx%x 
      atm_j = blockIdx%y 
      atm_k = blockIdx%z 
      ! atm_st = findAtoms_d(1, atm_i,atm_j,atm_k)
      ! atm_en = findAtoms_d(2, atm_i,atm_j,atm_k)
      atm_st = 1
      atm_en = dim_matAtm_d(atm_i,atm_j,atm_k)
      if (atm_st>0) then
        do i_list = atm_st, atm_en
          ! i_atm = listAtoms_d(i_list)
          i_atm = matrixAtoms_d(i_list, atm_i,atm_j,atm_k)
          
          if ( myfluid_d(i,j,k, flip)==fluid_spherelist ) then
            xxx = nint(pos_d(1,i_atm, flip))
            yyy = nint(pos_d(2,i_atm, flip))
            zzz = nint(pos_d(3,i_atm, flip))

            deltai = getDeltai(i, xxx)
            deltaj = getDeltaj(j, yyy)
            deltak = getDeltak(k+ offset_d(3), zzz)

            atmDist = deltai*deltai + deltaj*deltaj + deltak*deltak
            ! only call P2N_bb if near THIS atom...
            if (atmDist < sphereMax_d) then
#ifdef CHECK_WRONGATOM
              if (debugfluid_d(i,j,k) /= i_atm) then
                write(*,*) 'partBB_CG_p2n2]Wrong atom', linear(i,j,k), &
                    i_atm+0.001*atmDist, debugfluid_d(i,j,k)+0
                ! stop_d = __LINE__
              endif
#endif
              call particle_to_node_bounce_back_bc2_phase2(step, flip, i_atm, i,j,k, lrotate)
            endif
          endif
        enddo
      endif
    end subroutine partBB_CG_p2n2

    ! PARTICLES 	-----------------------------------------------
    ! pos(1:3) -> xxx,yyy,zzz
    ! pos(4:6) -> vxx,vyy,vzz
    ! pos(7:9) -> oxx,oyy,ozz
    ! pos(10:13) -> q0,q1,q2,q3
    ! force(1:3) -> fxx,fyy,fzz
    ! force(4:6) -> fxb,fyb,fzb
    ! force(7:9) -> tqx,tqy,tqz
    ! force(10:12) -> tqb,tqb,tqb

    attributes(global) subroutine inter_force_particles(step, flip, numAtoms)
      integer, value    :: step,flip,numAtoms
      integer           :: i,j

      i = (blockIdx%x-1) * TILE_DIMPART + threadIdx%x
      if (i>numAtoms) return

      if (offset_d(3)==0 .and. 1==step .and. 1==i) write(*,*) 'CUDA inter_force_particles] i:',i

      if (withSidewall) call add_sidewall(step, flip, i)

      if (withHz) then
        do j = i+1, numAtoms	! loop over succ parts
          call addHertzian(step, flip, i,j)
        enddo
      endif
    end subroutine inter_force_particles

    attributes(device) subroutine zero_forces_part(step, i, lrotate)
    integer, value    :: step, i
    logical, value    :: lrotate

    !Zero forces

    ! Queste sono BB
    myf_d(1,i) = 0.0
    myf_d(2,i) = 0.0
    myf_d(3,i) = 0.0
    
    ! myf2_d(1,i) = 0.0
    ! myf2_d(2,i) = 0.0
    ! myf2_d(3,i) = 0.0

    if (lrotate) then

      ! Queste sono BB per torque
      myt_d(1,i) = 0.0
      myt_d(2,i) = 0.0
      myt_d(3,i) = 0.0
      
      ! myt2_d(1,i) = 0.0
      ! myt2_d(2,i) = 0.0
      ! myt2_d(3,i) = 0.0
    endif
    end subroutine zero_forces_part

    attributes(global) subroutine time_step_force_particles(step, flip, lrotate, numAtoms)
      integer, value    :: step,flip,numAtoms
      logical, value    :: lrotate
      integer           :: i
      integer :: flop
  
      i = (blockIdx%x-1) * TILE_DIMPART + threadIdx%x
      if (i>numAtoms) return

      if (offset_d(3)==0 .and. 1==step .and. 1==i) write(*,*) 'CUDA time_step_force_particles] i:',i

      ! Only for DEBUG!!!!
      ! call zero_forces_part(step, i, lrotate)

      call addForce_particle_bounce_back(step, flip, i, lrotate)
      call nve_lf(step, flip, i, lrotate)

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

      if (oldpartVol_d(i)>0 .and. partVol_d(i) /= oldpartVol_d(i)) then        
        ! write(*,*) 'partVol_d /= oldpartVol_d ->Atom', i, partVol_d(i), oldpartVol_d(i)
        ! stop_d = __LINE__      
      endif
      oldpartVol_d(i) = partVol_d(i)
      partVol_d(i) = 0

      if (i == 1) then
        if (countmk_d /= countrm_d) then
          ! write(*,*) 'countmk_d /= countrm_d ->', countmk_d, countrm_d
          ! stop_d = __LINE__
        endif

        countmk_d = 0
        countrm_d = 0
        countn2p_d = 0
        countp2n1_d = 0
        countp2n2_d = 0
      endif
    end subroutine time_step_force_particles


    attributes(device) subroutine add_sidewall(step, flip, iatm)
      integer, value    :: step, flip, iatm
      real, parameter	  :: inflimit = sidewall_rdist + HALF, rlimit = ONE+HALF, gmin = FIVE*HALF * sidewall_k
      real, parameter	  :: suplimitx = nx - inflimit
      real, parameter	  :: suplimity = ny - inflimit
      real	  :: suplimitz
      real, parameter	  :: suprcapx = nx - 1.5
      real, parameter	  :: suprcapy = ny - 1.5
      real	  :: suprcapz
      real              :: rrr,ggg

      suplimitz = nz_d - inflimit
      suprcapz = nz_d - 1.5

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


    attributes(device) subroutine addForce_particle_bounce_back(step, flip, i, lrotate)
      integer, value    :: step, flip, i
      logical, value    :: lrotate
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

      debugline_d(i, 8,4) = partVol_d(i)
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
      real(8) :: acc(3), accFoo1,accFoo2,accFoo3

      real, parameter	:: rsqcut = hz_rmin**TWO
      real, parameter	:: gmin = FIVE*HALF * hz_k * (hz_rmin-hz_rcap)**(THREE*HALF)      
      real, parameter	:: rpar = 2*rdim
      real, parameter	:: lubfactor = lubric_k * SIX*Pi*  (rdim*rdim)**TWO / (rpar*rpar)
      real, parameter	:: rparcut = rpar + lubric_rmin
      real, parameter	:: rsrparcut = rparcut**TWO

      mxrsqcut = max(rsrparcut,rsqcut)
      rparcap = rpar + lubric_rcap

      if (offset_d(3)==0 .and. 1==step.and.2==iatm*jatm) write(*,*) 'Hz]', mxrsqcut, rparcap, gmin

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
      if (zperiodic) zdf = zdf - glz*nint(zdf / glz)

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

          ! Double prec
          acc(1) = ftemp(1)
          acc(2) = ftemp(2)
          acc(3) = ftemp(3)          
          accFoo1 = atomicadd( myf2_d(1,iatm), -acc(1) )
          accFoo2 = atomicadd( myf2_d(2,iatm), -acc(2) )
          accFoo3 = atomicadd( myf2_d(3,iatm), -acc(3) )
          accFoo1 = atomicadd( myf2_d(1,jatm), acc(1) )
          accFoo2 = atomicadd( myf2_d(2,jatm), acc(2) )
          accFoo3 = atomicadd( myf2_d(3,jatm), acc(3) )
          ! Done double

          ! write(*,*) '<=hz_rmin]', linear(0,iatm,jatm), xdf,ydf,zdf, rrr
          ! write(*,*) '<=hz_rmin]', linear(0,iatm,jatm), ftemp(1),ftemp(2),ftemp(3)
          ! write(*,*) '<=hz_rmin]', iatm, linear( nint(xi(1)),nint(xi(2)),nint(xi(3)) )
          ! write(*,*) '<=hz_rmin]', jatm, linear( nint(xj(1)),nint(xj(2)),nint(xj(3)) )
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
              ! zcm = pimage(izpbc,zcm,nz_d)
              ! visc=omega_to_viscosity(omega(xcm,ycm,zcm))
            ! endif
            
            dotuv = ux*(pos_d(4, iatm, flip)-pos_d(4, jatm, flip)) + uy*(pos_d(5, iatm, flip)-pos_d(5, jatm, flip)) &
                + uz*(pos_d(6, iatm, flip)-pos_d(6, jatm, flip))

            ftemp(1) = visc*lubfactor*ux* dotuv*(ONE/(rrr-rpar)-ONE)
            ftemp(2) = visc*lubfactor*uy* dotuv*(ONE/(rrr-rpar)-ONE)
            ftemp(3) = visc*lubfactor*uz* dotuv*(ONE/(rrr-rpar)-ONE)

            forceAtoms_d(1,iatm,flip) = forceAtoms_d(1,iatm,flip) - ftemp(1)
            forceAtoms_d(1,jatm,flip) = forceAtoms_d(1,jatm,flip) + ftemp(1)
                                                        
            forceAtoms_d(2,iatm,flip) = forceAtoms_d(2,iatm,flip) - ftemp(2)
            forceAtoms_d(2,jatm,flip) = forceAtoms_d(2,jatm,flip) + ftemp(2)

            forceAtoms_d(3,iatm,flip) = forceAtoms_d(3,iatm,flip) - ftemp(3)
            forceAtoms_d(3,jatm,flip) = forceAtoms_d(3,jatm,flip) + ftemp(3)

            write(*,*) 'hz<=rparcut] i:', iatm, ftemp(1),ftemp(2),ftemp(3)
          endif
        endif
      endif
    end subroutine addHertzian


    attributes(device) subroutine nve_lf(step, flip, i, lrotate)
      integer, value    :: step,flip, i
      logical, value    :: lrotate
      integer :: flop
      real :: fxx,fyy,fzz
      real :: bxx,byy,bzz, vxx,vyy,vzz
      real :: xxx,yyy,zzz, xxo,yyo,zzo
      logical :: lfix_moment = .false.
      real :: velMag

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

      velMag = sqrt(bxx*bxx + byy*byy + bzz*bzz)
      if (velMag > max_vel_d) then
        write(*,*) 'CUDA nve_lf] i:', i
        write(*,*) 'Extreme vel:', velMag, bxx,byy,bzz
        stop_d = __LINE__
      endif
      
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
        if (zzz< 0.5 .or. zzz>glz+0.5) then
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


    attributes(device) subroutine node_to_particle_bounce_back_bc2(step, flip, iatm, i,j,k, lrotate)
      integer, value    :: step, flip
      integer, value    :: iatm, i,j,k
      logical, value    :: lrotate
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
        rversor(3) = k+ offset_d(3) - pos_d(3,iatm,flip)

        ! PBC treatment
        aaa = ONE/real(nx)
        bbb = ONE/real(ny)
        ccc = ONE/real(glz)
        rversor(1) = rversor(1) - real(nx)*nint(aaa*rversor(1))
        rversor(2) = rversor(2) - real(ny)*nint(bbb*rversor(2))
        rversor(3) = rversor(3) - real(glz)*nint(ccc*rversor(3))

        modr = sqrt( rversor(1)*rversor(1) + rversor(2)*rversor(2) + rversor(3)*rversor(3) )
        rversor(1:3) = rdim / modr * rversor(1:3)
      endif

      do iloop = 1, 9
        indlow = iloop*2 - 1
        indhig = iloop*2
      
        ii = i + ex(indlow)
        jj = j + ey(indlow)
        kk = k + ez(indlow)
        if (kk<0 .or. kk>nz_d+1) cycle

        io = i + ex(indhig)
        jo = j + ey(indhig)
        ko = k + ez(indhig)
        if (ko<0 .or. ko>nz_d+1) cycle
      
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
      ccc = ONE/real(glz)

      pbc_images_centeredz = zzs - real(glz)*nint(ccc*zzs) + atomz
    end function

    attributes(device) subroutine particle_to_node_bounce_back_bc2_phase1(step, flip, i,j,k)
      integer, value    :: step, flip
      integer, value    :: i,j,k
      integer           :: iloop, ii,jj,kk, io,jo,ko, indlow,indhig
      integer    :: tempInt
      
      do iloop = 1, 9
        indlow = iloop*2 - 1
        indhig = iloop*2

        ii = i + ex(indlow)
        jj = j + ey(indlow)
        kk = k + ez(indlow)
        if (kk<0 .or. kk>nz_d+1) cycle

        io = i + ex(indhig)
        jo = j + ey(indhig)
        ko = k + ez(indhig)
        if (ko<0 .or. ko>nz_d+1) cycle

        if(myfluid_d(ii,jj,kk, flip)==fluid_fluid .and. myfluid_d(io,jo,ko, flip)/=fluid_fluid)then
#ifdef DEBUG_P2N
          tempInt = 1 + atomicadd( countp2n1_d, 1)
          debugp2n1_d(tempInt, 1) = linear(ii,jj,kk+ offset_d(3))
          debugp2n1_d(tempInt, 2) = indhig
          debugp2n1_d(tempInt, 3) = linear(i,j,k+ offset_d(3))
          debugp2n1_d(tempInt, 4) = indlow
          ! write(*,*) 'p2n1]',linear(ii,jj,kk), indhig, linear(i,j,k), indlow
#endif
          popsR_d(i,j,k,indlow, flip) = popsR_d(ii,jj,kk,indhig, flip)
          popsB_d(i,j,k,indlow, flip) = popsB_d(ii,jj,kk,indhig, flip)
        endif

        if(myfluid_d(ii,jj,kk, flip)/=fluid_fluid .and. myfluid_d(io,jo,ko, flip)==fluid_fluid)then
#ifdef DEBUG_P2N
          tempInt = 1 + atomicadd( countp2n1_d, 1)
          debugp2n1_d(tempInt, 1) = linear(io,jo,ko+ offset_d(3))
          debugp2n1_d(tempInt, 2) = indlow
          debugp2n1_d(tempInt, 3) = linear(i,j,k+ offset_d(3))
          debugp2n1_d(tempInt, 4) = indhig
          ! write(*,*) 'p2n1]',linear(io,jo,ko), indlow, linear(i,j,k), indhig
#endif
          popsR_d(i,j,k,indhig, flip) = popsR_d(io,jo,ko,indlow, flip)
          popsB_d(i,j,k,indhig, flip) = popsB_d(io,jo,ko,indlow, flip)
        endif
      enddo
    end subroutine particle_to_node_bounce_back_bc2_phase1
    
    attributes(device) subroutine particle_to_node_bounce_back_bc2_phase2(step, flip, iatm, i,j,k, lrotate)
      integer, value    :: step, flip
      integer, value    :: iatm, i,j,k
      logical, value    :: lrotate
      real              :: rtemp(3), ftemp(3),ttemp(3), rversor(3), otemp(3), modr, f2pR,f2pB
      real              :: qtemp(0:3),qversor(0:3),tempconj(0:3), oat(0:3)
      real              :: vx,vy,vz, vxs,vys,vzs
      real              :: aaa,bbb,ccc
      integer           :: iloop, ii,jj,kk, io,jo,ko, indlow,indhig
      integer    :: tempInt
      
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
        rversor(3) = k+ offset_d(3) - pos_d(3,iatm,flip)
        ! PBC treatment
        aaa = ONE/real(nx)
        bbb = ONE/real(ny)
        ccc = ONE/real(glz)
        rversor(1) = rversor(1) - real(nx)*nint(aaa*rversor(1))
        rversor(2) = rversor(2) - real(ny)*nint(bbb*rversor(2))
        rversor(3) = rversor(3) - real(glz)*nint(ccc*rversor(3))
        
        modr = sqrt( rversor(1)*rversor(1) + rversor(2)*rversor(2) + rversor(3)*rversor(3) )
        rversor(1:3) = rdim / modr * rversor(1:3)
      endif

      do iloop = 1, 9
        indlow = iloop*2 - 1
        indhig = iloop*2
      
        ii = i + ex(indlow)
        jj = j + ey(indlow)
        kk = k + ez(indlow)
        if (kk<0 .or. kk>nz_d+1) cycle

        io = i + ex(indhig)
        jo = j + ey(indhig)
        ko = k + ez(indhig)
        if (ko<0 .or. ko>nz_d+1) cycle

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
          tempInt = 1 + atomicadd( countp2n2_d, 1)
          debugp2n2_d(tempInt, 1) = linear(ii,jj,kk+ offset_d(3))
          debugp2n2_d(tempInt, 2) = indhig
          debugp2n2_d(tempInt, 3) = f2pR
          debugp2n2_d(tempInt, 4) = f2pB
          ! write(*,*) 'R+B p2nx]', linear(ii,jj,kk), indhig, f2pR, f2pB
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
          tempInt = 1 + atomicadd( countp2n2_d, 1)
          debugp2n2_d(tempInt, 1) = linear(io,jo,ko+ offset_d(3))
          debugp2n2_d(tempInt, 2) = indlow
          debugp2n2_d(tempInt, 3) = f2pR
          debugp2n2_d(tempInt, 4) = f2pB
          ! write(*,*) 'R+B p2nx]', linear(io,jo,ko), indlow, f2pR, f2pB
#endif
        endif
      enddo
    end subroutine particle_to_node_bounce_back_bc2_phase2
#endif    

    attributes(global) subroutine zeroListAtomsGPU(step, flip)
      integer, value    :: step,flip
      integer           :: i,j,k


      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
      if (i>nx/TILE_DIMx) return
      if (j>ny/TILE_DIMy) return
      if (k>nz_d/TILE_DIMz) return

      dim_matAtm_d(i,j,k) = 0
    end subroutine zeroListAtomsGPU

    attributes(global) subroutine makeListAtomsGPU(step, flip, numAtoms)
      integer, value    :: step,flip,numAtoms
      real           :: x,y,z
      integer           :: i,j,k, l
      integer           :: tempInt
      integer           :: boxatom(3,4)
      logical           :: splitted(3)
      integer           :: loopi,loopj,loopk
      integer           :: mini,minj,mink, maxi,maxj,maxk


      l = (blockIdx%x-1) * TILE_DIMPART + threadIdx%x
      if (l>numAtoms) return

        i = ceiling( pos_d(1,l, flip) )
        j = ceiling( pos_d(2,l, flip) )
        k = ceiling( pos_d(3,l, flip) )
        x = pos_d(1,l, flip)
        y = pos_d(2,l, flip)
        z = pos_d(3,l, flip)
        boxatom(:,1) = (/ x-rdim, y-rdim, z-rdim /)
        boxatom(:,2) = (/ x+rdim, y+rdim, z+rdim /)

        ! x-component
        if (xperiodic .and. boxatom(1,1)< 1) then 
          splitted(1) = .true.
          boxatom(1,1) = 1
          boxatom(1,2) = i+rdim
          boxatom(1,3) = i-rdim + nx
          boxatom(1,4) = nx
        elseif (xperiodic .and. boxatom(1, 2) > nx) then 
          splitted(1) = .true.
          boxatom(1,1) = 1
          boxatom(1,2) = i+rdim - nx
          boxatom(1,3) = i-rdim
          boxatom(1,4) = nx
        else
          splitted(1) = .false.
          boxatom(1,3) = TILE_DIMx	! This avoid loopi==1
          boxatom(1,4) = -TILE_DIMx
        endif

        ! y-component
        if (yperiodic .and. boxatom(2,1)< 1) then 
          splitted(2) = .true.
          boxatom(2,1) = 1
          boxatom(2,2) = j+rdim
          boxatom(2,3) = j-rdim + ny
          boxatom(2,4) = ny
        elseif (yperiodic .and. boxatom(2, 2) > ny) then 
          splitted(2) = .true.
          boxatom(2,1) = 1
          boxatom(2,2) = j+rdim - ny
          boxatom(2,3) = j-rdim
          boxatom(2,4) = ny
        else
          splitted(2) = .false.
          boxatom(2,3) = TILE_DIMy	! This avoid loopj==1
          boxatom(2,4) = -TILE_DIMy
        endif

        ! z-component
        if (zperiodic .and. boxatom(3,1)< 1) then 
          splitted(3) = .true.
          boxatom(3,1) = 1
          boxatom(3,2) = k+rdim
          boxatom(3,3) = k-rdim + glz
          boxatom(3,4) = glz
        elseif (zperiodic .and. boxatom(3, 2) > glz) then 
          splitted(3) = .true.
          boxatom(3,1) = 1
          boxatom(3,2) = k+rdim - glz
          boxatom(3,3) = k-rdim
          boxatom(3,4) = glz
        else
          splitted(3) = .false.
          boxatom(3,3) = TILE_DIMz	! This avoid loopk==1
          boxatom(3,4) = -TILE_DIMz
        endif

      do loopk=0, 1
	mink = ( boxatom(3,2*loopk+1)-1 ) / TILE_DIMz + 1
	maxk = ( boxatom(3,2*loopk+2)-1 ) / TILE_DIMz + 1
  
	do k = mink- offset_d(3)/TILE_DIMz, maxk- offset_d(3)/TILE_DIMz
	  if (k<1 .or. k > nz_d/TILE_DIMz) cycle

      	  do loopj=0, 1
      	    minj = (boxatom(2,2*loopj+1)-1) / TILE_DIMy + 1
      	    maxj = (boxatom(2,2*loopj+2)-1) / TILE_DIMy + 1
            
      	    do j = minj, maxj
	      if (j<1 .or. j > ny/TILE_DIMy) cycle

      	      do loopi=0, 1
      	        mini = (boxatom(1,2*loopi+1)-1) / TILE_DIMx + 1
      	        maxi = (boxatom(1,2*loopi+2)-1) / TILE_DIMx + 1

      	        do i = mini, maxi
	          if (i<1 .or. i > nx/TILE_DIMx) cycle

      	          ! put atom l in list(i,j,k)
		  tempInt = 1 + atomicadd( dim_matAtm_d(i,j,k), 1)
      ! write(*,*) 'makeListAtomsGPU', step, offset_d(3), l,  linear(i*TILE_DIMx,j*TILE_DIMy,k*TILE_DIMz)

		  if (tempInt > maxAtomInBlock) then
            	     write(*,*) 'Too many atoms in CUDA block:', linear(i,j,k), tempInt, maxAtomInBlock
            	     stop_d = __LINE__
		  endif

		  matrixAtoms_d(tempInt, i,j,k) = l
      	        enddo
      	      enddo

      	    enddo
      	  enddo

	enddo
      enddo

    end subroutine makeListAtomsGPU

  end module kernels_fluid_part