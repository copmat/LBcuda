program TestLB
#include "defines.h"

  use cudafor
  use dimensions_m
  use profiling_m,   only : timer_init,itime_start, &
                        startPreprocessingTime,print_timing_partial, &
                        reset_timing_partial,printSimulationTime, &
                        print_timing_final,itime_counter,idiagnostic, &
                        ldiagnostic,start_timing2,end_timing2, &
                        set_value_ldiagnostic,set_value_idiagnostic, &
                        startSimulationTime,print_memory_registration, &
                        get_memory_cuda,print_memory_registration_cuda, &
                        get_memory
  use kernels_fluid
  use kernels_bgk
  use kernels_scp
  use kernels_fluid_cg
  use kernels_fluid_part
  use write_output

#ifdef SERIAL
! do nothing  
#else
  use mpi
#endif

  implicit none

  integer(4) :: nz


  ! host arrays
  real(4), allocatable, pinned :: pop_pinned(:,:,:,:)  ! pinned (aka page-locked) host memory, must be allocatable
  real(4), allocatable, pinned :: popB_pinned(:,:,:,:)  ! pinned (aka page-locked) host memory, must be allocatable
  real(4), allocatable, pinned :: popSCP_pinned(:,:,:,:,:)
  real(4), allocatable :: rhoR(:,:,:), rhoB(:,:,:), vel(:,:,:,:), phase(:,:,:)
  real(4), allocatable :: scalar(:,:,:,:)
  real(4), allocatable :: x_atm(:,:), v_atm(:,:), q(:,:)
  real(8), allocatable :: buffer_for(:,:),buffer_forall(:,:)
  real(4), allocatable :: pos(:,:,:)
  real(4), allocatable :: forceAtoms(:,:,:)
  integer(1), allocatable :: myfluid(:,:,:,:)
  integer(1) :: prova1
  integer, allocatable, dimension(:) :: inidom,findom
  integer :: nsphere=0,nspheredead=0, totSphere
  integer, allocatable :: spherelist(:,:),spherelistdead(:,:)

  ! events for timing
  type (cudaEvent) :: startEvent, stopEvent
  integer(8)	   :: clock_rate, c_stopAtom,c_startAtom, diff_listatomGPU

  ! misc
  type (dim3) :: dimGrid,dimBlock, dimGridx,dimGridy,dimGridz,dimBlock2
  type (dim3) :: dimGridhalo,dimBlockhalo 
  type (dim3) :: dimGridhalox,dimGridhaloy,dimGridhaloz
  type (dim3) :: dimGridsidex,dimGridsidey,dimGridsidez
  
  integer     :: dimGridAtm, dimBlockAtm
  type (dim3) :: dimGridSubd, dimBlockSubd
  type (cudaDeviceProp) :: prop
  real(4) :: mytime,mymemory,totmemory,myramhost
  integer :: istat, flipflop, nIter,nIterOut,nIterVTK,step,nIterVTKmin
  integer :: nIterVTK2d=0
  integer :: numVTK2d=0
  integer, allocatable, dimension(:) :: ivecVTK2d,iselVTK2d
  integer :: i,j,k
  logical :: pinnedFlag = .true.

  ! Phys params
  ! Initial fluid velocity
  real :: vx,vy,vz
  ! Input type
  integer :: initType = -1
  
  ! Debug
  real,allocatable :: debugline(:,:,:)
  integer,allocatable :: partVol(:)
  logical :: debug1 = .false.


  step = 0
  flipflop = 1

  stop_d = 0
  
  
  gsizes(1)=glx
  gsizes(2)=gly
  gsizes(3)=glz
#ifdef SERIAL
  myrank = 0
  nprocs = 1
  nz = glz
  nz_d = glz
  offset(1)= 0
  offset(2)= 0
  offset(3)= 0
#else
  call mpi_init(ierr)
  call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  nz = glz / nprocs
  if(mod(glz,nprocs).ne.0)then
    write (*,*) 'Global gnz is not divisible for the number of MPI tasks'
    write (*,*) 'mod(glz,nprocs) = ',mod(glz,nprocs)
    call mystop
  endif
  call setupMPI(nz)
  nz_d = glz / nprocs
#endif
  lsizes(1)=nx
  lsizes(2)=ny
  lsizes(3)=nz
  start_idx(1)=1
  start_idx(2)=1
  start_idx(3)=offset(3)+1
  
  allocate(inidom(0:nprocs-1))
  allocate(findom(0:nprocs-1))
  inidom(:)=0
  findom(:)=-1
  inidom(0)=1
  findom(0)=nz+inidom(0)-1
  do i=1,nprocs-1
    inidom(i)=findom(i-1)+1
    findom(i)=inidom(i)+nz-1
  enddo
  
  
  call readParameters

  call setupGPUgridAndAlloc
#ifdef PARTICLES
  if (withParticles) then
	  call SetupParts
  endif
#endif
  call makeOutput
  call SetConstants
  
  call set_value_ldiagnostic(diagnostic)
  call set_value_idiagnostic(tdiagnostic)
  ! start diagnostic if requested
  if(ldiagnostic)then
    call timer_init()
    call startPreprocessingTime()  
  endif
  
  
#ifdef PARTICLES
  call spherical_template
#endif

  call setupRho_Pops
  
  if (withCG) then
    call fixPeriodic(.true., .false., "after_setup")
    call fixPeriodic_hvar(.true.,.true.,.false.,"after_setup")
  else
    call fixPeriodic_BGK(.true., .false., "after_setup")
    call fixPeriodic_hvar(.true.,.false.,.false.,"after_setup")
  endif
  call fixPeriodic_isfluid("after_setup")
  call copy_isfluid<<<dimGridhalo,dimBlockhalo>>>(step, flipflop)
!  myfluid = myfluid_d
!  istat = cudaDeviceSynchronize
!  open(unit=16,file='mio.data',status='replace')
!  do k=1-nbuff,nz+nbuff
!    do j=1-nbuff,ny+nbuff
!      do i=1-nbuff,nx+nbuff
!        write(16,'(5i4)')i,j,k,myfluid(i,j,k,1:2)
!      enddo
!    enddo
!  enddo
  
  istat = cudaMemcpy(myfluid,myfluid_d, (nx+2*nbuff)*(ny+2*nbuff)*(nz+2*nbuff)*2 )
  if (istat/=0) write(*,*) 'status after myfluid:',istat
  istat = cudaDeviceSynchronize
  if (istat/=0) write(*,*) 'status after myfluid',istat

  if (wantOut) call OutputVTK("step_", 0, 1)
!  istat = cudaDeviceSynchronize
!  if (myrank== 0) write(*,*) 'ciaone'
!  call MPI_FINALIZE(ierr)
!  stop
  
  step = 1

  if (myrank== 0) then
    call printParams
    write(*,*) 'Run...'
  endif
  
  ! start diagnostic if requested
  if(ldiagnostic)then
    !call print_timing_partial(1,1,itime_start,6)
    !call reset_timing_partial()
    call startSimulationTime()
    call get_memory_cuda(mymemory,totmemory)
    call print_memory_registration_cuda(6,'DEVICE memory occupied at the start', &
       'total DEVICE memory',mymemory,totmemory)
    call get_memory(myramhost)
    call print_memory_registration(6,'HOST memory occupied at the start', &
     myramhost)
  endif
  
  
  istat = cudaEventRecord(startEvent, 0)


  do step=1, nIter

    if (withCG) then
#ifdef PARTICLES
      if (withParticles) then
        call init_isfluid_CG<<<dimGrid, dimBlock>>>(step, flipflop)
        call fixPeriodic_isfluid("step")
        call del_fluid_CG<<<dimGrid, dimBlock>>>(step, flipflop, lrotate)

#ifdef CHECK_VOLP
        partVol(1:numAtoms) = partVol_d(1:numAtoms)

        call dumpPartvol(step)
#endif
      endif


#ifdef FLUIDSPLIT 
      if(ldiagnostic)call start_timing2("LB","moments")
      call init_rhoB_CG<<<dimGrid, dimBlock>>>(step, flipflop) 
      call init_rhoR_CG<<<dimGrid, dimBlock>>>(step, flipflop) 
      if(ldiagnostic)call end_timing2("LB","moments")
#else
      call init_rho_CG<<<dimGrid, dimBlock>>>(step, flipflop)
#endif

      call fixPeriodic(.false., .true., "init_rho_CG")
      if (debug1) call OutputVTK("rho_isfluid_",step, flipflop)


      if (withParticles) then
        call make_fluid_CG<<<dimGrid, dimBlock>>>(step, flipflop, lrotate)
        call fixPeriodic(.false., .true., "make_fluid_CG")

#ifdef DEBUG_MKRM
        call dumpMkRm(step)
#endif

      endif

      if (withParticles) then
        ! call compute_densities_wall<<<dimGrid, dimBlock>>>(step, flipflop)
        call compute_virtual_rho<<<dimGrid, dimBlock>>>(step, flipflop)
        
        call fixPeriodic(.false., .true., "compute_densities_wall")
        if (debug1) call OutputVTK("wall_",step, flipflop)
      endif


      call time_step_CG<<<dimGrid, dimBlock>>>(step, flipflop)
      if (debug1) call OutputVTK("cg_",step, flipflop)

      call fixPeriodic(.true., .false., "time_step_CG")
      

      if (withParticles) then
        call partBB_CG<<<dimGrid, dimBlock>>>(step, flipflop, lrotate)
        call abortOnLastErrorAndSync("partBB_CG",step, flipflop)

        ! call OutputVTK("partBB_n2p_",step, flipflop)

        call partBB_CG_p2n<<<dimGrid, dimBlock>>>(step, flipflop)          
        ! call OutputVTK("partBB_p2n1_",step, flipflop)

        call partBB_CG_p2n2<<<dimGrid, dimBlock>>>(step, flipflop, lrotate)          
        if (debug1) call OutputVTK("partBB_",step, flipflop)

        call fixPeriodic(.true., .false., "partBB_CG_p2n2")

#ifdef DEBUG_N2P
        call dumpTerms_n2p(step)
#endif
#ifdef DEBUG_P2N
        call dumpTerms_p2n(step)
#endif
      endif


      ! Boundary Conditions with new flipflop
      if (debug1) call OutputVTK("bb_",step, flipflop)
      
! 
#ifdef FLUIDSPLIT 
      !call flipflopRPop0_CG<<<dimGrid, dimBlock>>>(step, flipflop)
      call streamR_CG<<<dimGrid, dimBlock >>>(step, flipflop)

      !call flipflopBPop0_CG<<<dimGrid, dimBlock>>>(step, flipflop)
      call streamB_CG<<<dimGrid, dimBlock >>>(step, flipflop) 
! 
#else 
      call flipflopPop0_CG<<<dimGrid, dimBlock>>>(step, flipflop)
      call stream_CG<<<dimGrid, dimBlock>>>(step, flipflop) 
#endif 

      call stream_CG_x<<<dimGridx, dimBlock2>>>(step, flipflop)
      call stream_CG_y<<<dimGridy, dimBlock2>>>(step, flipflop)
      call stream_CG_z<<<dimGridz, dimBlock2>>>(step, flipflop)
      

      if (withParticles) then
        call inter_force_particles<<<dimGridAtm, dimBlockAtm>>>(step, flipflop, numAtoms)
# ifdef SERIAL
# else
        call abortOnLastErrorAndSync("inter_force_particles",step, flipflop)
        call mpiAddForces()
# endif        
        call time_step_force_particles<<<dimGridAtm, dimBlockAtm>>>(step, flipflop, lrotate, numAtoms)
        
#ifdef DEBUG_FORCE
        call dumpForces(step)
#endif

#ifdef DEBUG_ROT
        if (lrotate) call dumpRot(step, flipflop)
#endif

        call updateListAtoms(3 - flipflop)  ! Use new pos
      endif
! if not def PARTICLE      
#else


      if(ldiagnostic)call start_timing2("LB","moments")
      call init_rhoB_CG<<<dimGrid, dimBlock>>>(step, flipflop) 
      call init_rhoR_CG<<<dimGrid, dimBlock>>>(step, flipflop) 
      if(ldiagnostic)call end_timing2("LB","moments")
      
      if(ldiagnostic)call start_timing2("LB","fixPhvar")
      call fixPeriodic_hvar(.true.,.true.,.false.,"init_rho_CG")
      if(ldiagnostic)call end_timing2("LB","fixPhvar")
      
#ifdef NEARCONTACT
      if(ldiagnostic)call start_timing2("LB","nearsel")
      call init_gradnearsel_CG<<<dimGrid, dimBlock>>>(step, flipflop)
      call fixPeriodic_nearsel("nearsel_CG")
      if(ldiagnostic)call end_timing2("LB","nearsel")
#endif      
      
      
      if (debug1) call OutputVTK("rho_isfluid_",step, flipflop)

#ifdef APPLYBC 
      if(.not. fixdenswall)then
        if(ldiagnostic)call start_timing2("LB","fixWall")
        call compute_meandenswall("fixWall")
        if(ldiagnostic)call end_timing2("LB","fixWall")
      endif
#endif 


      if(ldiagnostic)call start_timing2("LB","time_step")
      call time_step_CG<<<dimGrid, dimBlock>>>(step, flipflop)
      if(ldiagnostic)call end_timing2("LB","time_step")
      
      if (debug1) call OutputVTK("cg_",step, flipflop)
      
      if(ldiagnostic)call start_timing2("LB","fixPpop")
      call fixPeriodic(.true., .false., "time_step_CG")
      if(ldiagnostic)call end_timing2("LB","fixPpop")

#ifdef APPLYBC  
      if(ldiagnostic)call start_timing2("LB","applybc")
      call fixPeriodic_hvar(.false.,.false.,.true.,"applybchalo")
      call compute_applybc("applybchalo") 
      if(ldiagnostic)call end_timing2("LB","applybc")
#endif      


      ! Boundary Conditions with new flipflop
      if (debug1) call OutputVTK("bb_",step, flipflop)
      
! 
      if(ldiagnostic)call start_timing2("LB","streaming")
      call streamR_CG<<<dimGrid, dimBlock >>>(step, flipflop)
      call streamB_CG<<<dimGrid, dimBlock >>>(step, flipflop) 
! 

      call stream_CG_x<<<dimGridx, dimBlock2>>>(step, flipflop)
      call stream_CG_y<<<dimGridy, dimBlock2>>>(step, flipflop)
      call stream_CG_z<<<dimGridz, dimBlock2>>>(step, flipflop)
      if(ldiagnostic)call end_timing2("LB","streaming")
      
#endif
    else
#ifdef BGK
      ! if (forced) then
      !   if (const_forced) then
      !     call time_step_force_cost<<<dimGrid, dimBlock>>>(step, flipflop, omega, 1.0 - omega)
      !   else
      !     call lb_force<<<dimGrid, dimBlock>>>()
      !     call time_step_force<<<dimGrid, dimBlock>>>(step, flipflop, omega, 1.0 - omega)
      !   endif
      ! else
      !   call time_step<<<dimGrid, dimBlock>>>(step, flipflop, omega, 1.0 - omega)
      ! endif
      
      
      
      ! call OutputVTK("bb_",step, flipflop)
      
#ifdef MOMBGK      
      if(ldiagnostic)call start_timing2("LB","moments")
      call init_rho_isfluid_BGK<<<dimGrid, dimBlock>>>(step, flipflop)
      if(ldiagnostic)call end_timing2("LB","moments")
      
      if(ldiagnostic)call start_timing2("LB","fixPhvar")
      call fixPeriodic_hvar(.true.,.false.,.false.,"init_rho")
      if(ldiagnostic)call end_timing2("LB","fixPhvar")
      
      if(ldiagnostic)call start_timing2("LB","timestep")
      call time_step_BGK<<<dimGrid, dimBlock>>>(step, flipflop, omega)
      if(ldiagnostic)call end_timing2("LB","timestep")
#else
      if(ldiagnostic)call start_timing2("LB","timestep")
      call time_step_mom_BGK<<<dimGrid, dimBlock>>>(step, flipflop, omega)
      if(ldiagnostic)call end_timing2("LB","timestep")
      
      if(ldiagnostic)call start_timing2("LB","fixPhvar")
      call fixPeriodic_hvar(.true.,.false.,.false.,"init_rho")
      if(ldiagnostic)call end_timing2("LB","fixPhvar")
      
      if(withSCP)then
        if(ldiagnostic)call start_timing2("LB","timestepSCP")
        call time_step_mom_SCP<<<dimGrid, dimBlock>>>(step,flipflop,omega)
        if(ldiagnostic)call end_timing2("LB","timestepSCP")
        
        if(ldiagnostic)call start_timing2("LB","fixPsvar")
        call fixPeriodic_SCP(.false.,.true.,"init_scp")
        if(ldiagnostic)call end_timing2("LB","fixPsvar")
      endif
      
#endif
      

      if(ldiagnostic)call start_timing2("LB","fixPpops")
      call fixPeriodic_BGK(.true., .false., "time_step_BGK")
      if(ldiagnostic)call end_timing2("LB","fixPpops")
      
      if(withSCP)then
        if(ldiagnostic)call start_timing2("LB","fixPscp")
        call fixPeriodic_SCP(.true.,.false.,"fixPscp")
        if(ldiagnostic)call end_timing2("LB","fixPscp")
      endif
        
#ifdef APPLYBC  
      if(ldiagnostic)call start_timing2("LB","applybc")
      call fixPeriodic_hvar(.false.,.false.,.true.,"applybchalo_BGK")
      call applybchalo_BGK<<<dimGridhalo,dimBlockhalo >>>(step, flipflop) 
      if(ldiagnostic)call end_timing2("LB","applybc")
#endif
      
      
      if(ldiagnostic)call start_timing2("LB","stream")
      call stream_BGK<<<dimGrid, dimBlock >>>(step, flipflop)
      call stream_BGK_x<<<dimGridx, dimBlock2>>>(step, flipflop)
      call stream_BGK_y<<<dimGridy, dimBlock2>>>(step, flipflop)
      call stream_BGK_z<<<dimGridz, dimBlock2>>>(step, flipflop)
      if(ldiagnostic)call end_timing2("LB","stream")
      
      if(withSCP)then
        if(ldiagnostic)call start_timing2("LB","streamscp")
        call stream_SCP<<<dimGrid, dimBlock >>>(step, flipflop)
        call stream_SCP_x<<<dimGridx, dimBlock2>>>(step, flipflop)
        call stream_SCP_y<<<dimGridy, dimBlock2>>>(step, flipflop)
        call stream_SCP_z<<<dimGridz, dimBlock2>>>(step, flipflop)
        if(ldiagnostic)call end_timing2("LB","streamscp")
      endif
      
#endif
    endif

    call abortOnLastErrorAndSync("time step",step, flipflop)


    ! Work on new vals...
    flipflop = 3 - flipflop

    istat = cudaEventRecord(stopEvent, 0)
    istat = cudaEventSynchronize(stopEvent)
    istat = cudaEventElapsedTime(mytime, startEvent, stopEvent)


    if (mod(step,niterOut)==0) then
      if (myrank== 0) then
        write(*,*) ''
        write(*,fmt=100)  step, mytime*0.001, mytime/niterOut, lup * niterOut / mytime / 1000.
100     format('Iter:',I10,'   Time(s):', G10.4, '   ms/it:',G9.3, '   Mlup/s:',F8.0)
#ifdef PARTICLES
        if (withParticles) then
          if (diff_listatomGPU/real(clock_rate)>0.5) write(*,fmt=101)  step, diff_listatomGPU/real(clock_rate)
101       format('Iter:',I10, '   ListAtomGPU Time(s):', F10.3)
        
          diff_listatomGPU = 0
        endif
#endif
      endif

       istat = cudaEventRecord(startEvent, 0)
    endif

    if (wantOut) then
      ! if (mod(step,niterVTK)==0) call DumpBorderPlanez("step_",step, flipflop)      
      if (mod(step,nIterVTKmin)==0)then
        if(ldiagnostic)call start_timing2("IO","OutVTK")
        call OutputVTK("step_",step, flipflop)
        if(ldiagnostic)call end_timing2("IO","OutVTK")  
      endif
    endif
    
  enddo
  
  ! finalize and print the diagnostic data
  if(ldiagnostic)then
    call printSimulationTime()
    call print_timing_final(idiagnostic,itime_counter, &
     itime_start,1,1,6)
    call get_memory_cuda(mymemory,totmemory)
    call print_memory_registration_cuda(6,'DEVICE memory occupied at the end', &
       'total DEVICE memory',mymemory,totmemory)
    call get_memory(myramhost)
    call print_memory_registration(6,'HOST memory occupied at the end', &
     myramhost)
  endif


#ifdef SERIAL
! do nothing  
#else
  call MPI_FINALIZE(ierr)
#endif

  contains
  
   function GET_RANK_POINT(ii,jj,kk)
     
      implicit none

      integer, intent(in) :: ii,jj,kk
     
      integer :: k
      integer :: GET_RANK_POINT
  
        do k=0,nprocs-1
          if(kk<=findom(i))then
            GET_RANK_POINT=k
            exit
          endif
        enddo   
      
    end function GET_RANK_POINT
  
    pure function viscosity_to_omega(visc)
      use dimensions_m
      implicit none
      real, intent(in) :: visc
      real :: viscosity_to_omega
      
      viscosity_to_omega = ONE / ( visc /cssq  + HALF)
    end function viscosity_to_omega


    subroutine readParameters
#ifdef SERIAL
! do nothing  
#else
      use mpi
#endif

      implicit none
      integer, parameter :: maxlen=120
      integer, parameter :: inputio=25
      character(len = maxlen)  :: string,inipFile,arg
      integer               :: i, itest, found, max_vel_found,count
      logical :: lcheck,lexist
      real :: mydist
#ifdef NEWINPUT      
      namelist /simulation/ nIter,nIterOut,nIterVTK,initType, &
       lreadisfluid,lwriteisf,store_vel,numbc,lwriteisf_every,diagnostic, &
       nIterVTK2d,numVTK2d      !,numscp
      namelist /fluid/ vx,vy,vz,densR,densB,tauR,tauB,f_cost,const_forced,forced,&
       bctype,bcvel,bcrho,denswallR,denswallB,A_rep,sigma_CG,beta_CG
      namelist /passive/ bcscptype,bcscp
      namelist /particle/ withParticles,numAtoms,lrotate,ext_fxx,ext_fyy,&
      ext_fzz,ext_tqx,ext_tqy,ext_tqz,max_vel, &
       ivecVTK2d,iselVTK2d
     
#endif      
      nIter = 0.2*1000*1000*1000 / lup * 10
      nIterOut = nIter / 10
      nIter = 5
      nIterOut = 1
      nIterVTK = 100000000 !nIterOut
      ext_fxx=ZERO;ext_fyy=ZERO;ext_fzz=ZERO
      ext_tqx=ZERO;ext_tqy=ZERO;ext_tqz=ZERO
      withParticles=.false.
      lrotate=.false.
      
#ifdef NEWINPUT
      lcheck=.false.
      count = command_argument_count()
      if(count>=1)then
        i=1
        arg=repeat(' ',maxlen)
        call getarg(i, arg)
        inipFile=repeat(' ',maxlen)
        inipFile=trim(arg)
        if(myrank==0)write(6,*) 'file  = ',trim(inipFile)
        lcheck=.true.
      else
        do i = 1, count
          arg=repeat(' ',maxlen)
          call getarg(i, arg)
          if(lcheck)then
            inipFile=repeat(' ',maxlen)
            inipFile=trim(arg)
            if(myrank==0)write(6,*) 'file  = ',trim(inipFile)
            exit
          endif
          if(trim(arg)=='./lbCUDA')lcheck=.true.
        enddo
      endif
      if(len(trim(inipFile))==0)lcheck=.false. 
      if(.not. lcheck)then
        if(myrank==0)then
          write(6,*) 'error!'
          write(6,*) 'the command line should be'
          write(6,*) '[executable] [file]'
          write(6,*) 'file  = name of input file'
          write(6,*) 'STOP!'
        endif
        call mystop
      endif     
      inquire(file=trim(inipFile),exist=lexist)
      if(.not. lexist)then
        if(myrank==0)then
          write(6,*)'ERROR: file ',trim(inipFile),' does not exist!'
          write(6,*) 'the command line should be'
          write(6,*) '[executable] [file]'
          write(6,*) 'file  = name of input file'
        endif
        call mystop
      endif
      open(unit=inputio,file=trim(inipFile),status='old')
      
      read(inputio,nml=simulation)
      
      if(numbc>0)then
        allocate(bcvel(3,numbc), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of bcvel failed'
        endif
        bcvel(1:3,1:numbc) = 0.0
        if(withCG)then
          allocate(bctype(2,numbc), STAT=istat)
          if (istat /= 0) then
            write(*,*) 'Allocation of bctype failed'
          endif
          bctype(1:2,1:numbc) = 0
          
          allocate(bcrho(2,numbc), STAT=istat)
          if (istat /= 0) then
            write(*,*) 'Allocation of bcrho failed'
          endif
          bcrho(1:2,1:numbc) = 0.0
        else
          allocate(bctype(1,numbc), STAT=istat)
          if (istat /= 0) then
            write(*,*) 'Allocation of bctype failed'
          endif
          bctype(1,1:numbc) = 0
        
          allocate(bcrho(1,numbc), STAT=istat)
          if (istat /= 0) then
            write(*,*) 'Allocation of bcrho failed'
          endif
          bcrho(1,1:numbc) = 0.0
        endif
      endif
        if(withSCP)then
          allocate(bcscptype(numbc), STAT=istat)
          if (istat /= 0) then
            write(*,*) 'Allocation of particle bcscptype failed'
          endif
          bcscptype(1:numbc) = 0
          
          allocate(bcscp(numscp,numbc), STAT=istat)
          if (istat /= 0) then
            write(*,*) 'Allocation of particle bcscp failed'
          endif
          bcscp(1:numscp,1:numbc) = 0.0
        endif

      if(numVTK2d>0)then
        allocate(ivecVTK2d(numVTK2d), STAT=istat)
	    if (istat /= 0) then
          write(*,*) 'Allocation of ivecVTK2d failed'
        endif
        ivecVTK2d=0
        
        allocate(iselVTK2d(numVTK2d), STAT=istat)
	    if (istat /= 0) then
          write(*,*) 'Allocation of iselVTK2d failed'
        endif
        
      endif
      
      read(inputio,nml=fluid)
      read(inputio,nml=passive)
      read(inputio,nml=particle)
      close(inputio)
      
      if(numVTK2d>0)then
        nIterVTKmin=min(nIterVTK,nIterVTK2d)
        if(any(ivecVTK2d>3) .or. any(ivecVTK2d<1))then
          write (*,*) 'when ivecVTK2d should be from 1 to 3'
          write (*,*) 'now ivecVTK2d = ',ivecVTK2d
          call mystop
        endif
      else
        nIterVTKmin=nIterVTK
      endif
      
      if(numbc>0)then
        if(.not. store_vel)then
          write (*,*) 'when numbc>0 store_vel should be true'
          write (*,*) 'now store_vel = ',store_vel
          call mystop
        endif
      endif
#ifndef PBC
      if(.not. lreadisfluid)then
        write (*,*) 'when PBC is false lreadisfluid should be true'
        write (*,*) 'lreadisfluid = ',lreadisfluid
        call mystop
      endif
#endif

      if(withSCP)then
        if(.not. store_vel)then
          write (*,*) 'when withSCP is true store_vel should be true'
          write (*,*) 'now store_vel = ',store_vel
          call mystop
        endif
      endif
      
#ifdef APPLYBC
      if(nbuff<2)then
        write (*,*) 'when APPLYBC is true nbuff should be at least 2'
        write (*,*) 'now nbuff = ',nbuff
        call mystop
      endif
#endif      
      
#else
      
      open(unit=11, file='input.dat', status='old', form='FORMATTED', action='read', iostat=itest)
      if (itest/=0) then
        write (*,*) "input.dat] Error opening file. STOP"
        call mystop
      endif

      found = 0
      ! Header line
      ! read(11, '(a120)', end=120) string

      ! Number of iterations
      read(11, *, end=120) nIter
      found = 1
      
      ! Number of iterations for terminal output
      read(11, *, end=120) nIterOut
      found = 2
      
      ! Number of iterations for VTK output
      read(11, *, end=120) nIterVTK
      found = 3
      
      ! Number of iterations for input Type
      read(11, *, end=120) initType
      found = 4
      
      read(11, *, end=120) vx,vy,vz
      found = 5

      read(11, *, end=120) string, numAtoms
      found = 6
      if (string(1:3) == 'YES') then
        withParticles = .true.
      else if (string(1:3) == 'NO') then
        withParticles = .false.
      else
        found = 5
      endif

      if (withParticles) then
        max_vel_found = 0
        read(11, *, end=120) ext_fxx,ext_fyy,ext_fzz
        found = 7

        read(11, *, end=120) string, max_vel
        found = 8
        if (string(1:3) == 'MAX_PARTICLE_VEL') then
          max_vel_found = 1
        else
          found = 7
        endif

        read(11, *, end=120) string
        found = 9
        if (string(1:3) == 'YES') then
          lrotate = .true.
        else if (string(1:3) == 'NO') then
          lrotate = .false.
        else
          found = 8
        endif

        if (lrotate) then
          read(11, *, end=120) ext_tqx,ext_tqy,ext_tqz
          found = 10
        endif
      endif

      120 continue
      close(11)
      
      ! write (*,*) 'found',found, withParticles, lrotate
      
      
      
      if ((lrotate .and. found /= 10) .or. &
          (.not. lrotate .and. withParticles .and. found /= 9) .or. &
          .not. withParticles .and. found /= 6) then
        if (found<1) write (*,*) 'Number of iterations not found'
        if (found<2) write (*,*) 'nIterOut not found'
        if (found<3) write (*,*) 'nIterVTK not found'
        if (found<4) write (*,*) 'initType not found'
        if (found<5) write (*,*) 'vx,vy,vz,densR,densB not found'
        if (found<6) write (*,*) 'withParticles not found'
        if (withParticles .and. found<7) write (*,*) 'ext_fxx,ext_fyy,ext_fzz not found'
        if (withParticles .and. max_vel_found<1) write (*,*) 'MAX_PARTICLE_VEL not found'        
        if (withParticles .and. found<9) write (*,*) 'lrotate not found'
        if (lrotate .and. found<10) write (*,*) 'ext_tqx,ext_tqy,ext_tqz not found'
        call mystop
      endif
      
#endif

      
    end subroutine 


    subroutine fixPeriodic_isfluid(msg)
#ifdef SERIAL
      implicit none
      character(len=*), intent(in) :: msg

      
        if(xperiodic)&
         call isfluid_per_x<<<dimGridsidex, dimBlock2>>>(step, flipflop)
        if(yperiodic)&
         call isfluid_per_y<<<dimGridsidey, dimBlock2>>>(step, flipflop)
        if(zperiodic)&
         call isfluid_per_z<<<dimGridsidez, dimBlock2>>>(step, flipflop)

        if(xperiodic)&
         call isfluid_edge_x<<<(nx+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop)
        if(yperiodic)&
         call isfluid_edge_y<<<(ny+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop)
        if(zperiodic)&
         call isfluid_edge_z<<<(nz+2*nbuff+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop)
    
#else
      use mpi
      implicit none
      character(len=*), intent(in) :: msg

      integer      :: reqs_right(5)
      integer      :: status_up(MPI_STATUS_SIZE)
      integer      :: reqs_up(5)
      integer      :: status_down(MPI_STATUS_SIZE)
      integer      :: reqs_down(5)
      integer      :: tag

        
       if(xperiodic)& 
        call isfluid_per_x<<<dimGridsidex, dimBlock2>>>(step, flipflop)
       if(yperiodic)&
        call isfluid_per_y<<<dimGridsidey, dimBlock2>>>(step, flipflop)

        call abortOnLastErrorAndSync(msg,step, flipflop)

        ! if it is the same task, there's no need to send/recv data
        if(up(2) == myrank) then
         if(zperiodic)&
          call isfluid_per_z<<<dimGridsidez, dimBlock2>>>(step, flipflop)
        else
          ! 1) isend to up task 
          tag = 405
          call mpi_isend(myfluid_d(1-nbuff,1-nbuff,nz+1-nbuff, flipflop), 1, xyplane_int1, up(2), tag, &
                            lbecomm, reqs_up(1), ierr)
          !
          ! 2) recv from down task
          tag = 405
          call  mpi_recv(myfluid_d(1-nbuff,1-nbuff,1-nbuff, flipflop), 1, xyplane_int1, down(2), tag, &
                            lbecomm, status_up, ierr)
          !
          ! 3) isend to down task
          tag = 406
          call mpi_isend(myfluid_d(1-nbuff,1-nbuff, 1, flipflop), 1, xyplane_int1, down(2), tag, &
                            lbecomm, reqs_down(1), ierr)
          !
          ! 4) recv from up task
          tag = 406
          call  mpi_recv(myfluid_d(1-nbuff,1-nbuff,nz+1, flipflop), 1, xyplane_int1, up(2), tag, &
                            lbecomm, status_down, ierr)
          !
          ! 5) syncronize
          call mpi_wait(reqs_up(1), status_up, ierr)
          call mpi_wait(reqs_down(1), status_down, ierr)
        endif

        ! if it is the same task, there's no need to send/recv data
        if(up(2) == myrank) then
         if(xperiodic)&
          call isfluid_edge_x<<<(nx+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop)
         if(yperiodic)&
          call isfluid_edge_y<<<(ny+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop)
         if(zperiodic)&
          call isfluid_edge_z<<<(nz+2*nbuff+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop)
        else
         if(zperiodic)&
          call isfluid_edge_z<<<(nz+2*nbuff+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop)
        endif
      
#endif      
    end subroutine fixPeriodic_isfluid
    
    
    subroutine compute_meandenswall(msg)

      
      implicit none
      character(len=*), intent(in) :: msg

      call applybc_meandenswall<<<dimGrid,dimBlock>>>(step, flipflop)
      call applybc_meandenswall_x<<<dimGridx, dimBlock2>>>(step, flipflop)
      call applybc_meandenswall_y<<<dimGridy, dimBlock2>>>(step, flipflop)
      call applybc_meandenswall_z<<<dimGridz, dimBlock2>>>(step, flipflop)
      
    end subroutine compute_meandenswall
    
    subroutine compute_applybc(msg)

      
      implicit none
      
      character(len=*), intent(in) :: msg
      
      call applybc_CG_R<<<dimGrid,dimBlock>>>(step, flipflop) 
      call applybc_CG_B<<<dimGrid,dimBlock>>>(step, flipflop) 
      
      call applybc_CG_R_x<<<dimGridx, dimBlock2>>>(step, flipflop)
      call applybc_CG_R_y<<<dimGridy, dimBlock2>>>(step, flipflop)
      call applybc_CG_R_z<<<dimGridz, dimBlock2>>>(step, flipflop)
      
      call applybc_CG_B_x<<<dimGridx, dimBlock2>>>(step, flipflop)
      call applybc_CG_B_y<<<dimGridy, dimBlock2>>>(step, flipflop)
      call applybc_CG_B_z<<<dimGridz, dimBlock2>>>(step, flipflop)
      
    end subroutine compute_applybc
    
    subroutine fixPeriodic_nearsel(msg)
#ifdef SERIAL
      implicit none
      character(len=*), intent(in) :: msg

      
        if(xperiodic)&
         call nearsel_per_x<<<dimGridsidex, dimBlock2>>>(step, flipflop)
        if(yperiodic)&
         call nearsel_per_y<<<dimGridsidey, dimBlock2>>>(step, flipflop)
        if(zperiodic)&
         call nearsel_per_z<<<dimGridsidez, dimBlock2>>>(step, flipflop)

        if(xperiodic)&
         call nearsel_edge_x<<<(nx+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop)
        if(yperiodic)&
         call nearsel_edge_y<<<(ny+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop)
        if(zperiodic)&
         call nearsel_edge_z<<<(nz+2*nbuff+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop)
    
#else
      use mpi
      implicit none
      character(len=*), intent(in) :: msg

      integer      :: reqs_right(5)
      integer      :: status_up(MPI_STATUS_SIZE)
      integer      :: reqs_up(5)
      integer      :: status_down(MPI_STATUS_SIZE)
      integer      :: reqs_down(5)
      integer      :: tag

        
       if(xperiodic)& 
        call nearsel_per_x<<<dimGridsidex, dimBlock2>>>(step, flipflop)
       if(yperiodic)&
        call nearsel_per_y<<<dimGridsidey, dimBlock2>>>(step, flipflop)

        call abortOnLastErrorAndSync(msg,step, flipflop)

        ! if it is the same task, there's no need to send/recv data
        if(up(2) == myrank) then
         if(zperiodic)&
          call nearsel_per_z<<<dimGridsidez, dimBlock2>>>(step, flipflop)
        else
          ! 1) isend to up task 
          tag = 405
          call mpi_isend(nearsel_d(1-nbuff,1-nbuff,nz+1-nbuff), 1, xyplane_int1, up(2), tag, &
                            lbecomm, reqs_up(1), ierr)
          !
          ! 2) recv from down task
          tag = 405
          call  mpi_recv(nearsel_d(1-nbuff,1-nbuff,1-nbuff), 1, xyplane_int1, down(2), tag, &
                            lbecomm, status_up, ierr)
          !
          ! 3) isend to down task
          tag = 406
          call mpi_isend(nearsel_d(1-nbuff,1-nbuff, 1), 1, xyplane_int1, down(2), tag, &
                            lbecomm, reqs_down(1), ierr)
          !
          ! 4) recv from up task
          tag = 406
          call  mpi_recv(nearsel_d(1-nbuff,1-nbuff,nz+1), 1, xyplane_int1, up(2), tag, &
                            lbecomm, status_down, ierr)
          !
          ! 5) syncronize
          call mpi_wait(reqs_up(1), status_up, ierr)
          call mpi_wait(reqs_down(1), status_down, ierr)
        endif

        ! if it is the same task, there's no need to send/recv data
        if(up(2) == myrank) then
         if(xperiodic)&
          call nearsel_edge_x<<<(nx+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop)
         if(yperiodic)&
          call nearsel_edge_y<<<(ny+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop)
         if(zperiodic)&
          call nearsel_edge_z<<<(nz+2*nbuff+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop)
        else
         if(zperiodic)&
          call nearsel_edge_z<<<(nz+2*nbuff+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop)
        endif
      
#endif      
    end subroutine fixPeriodic_nearsel

    subroutine fixPeriodic(faiPops,faiRho, msg)
#ifdef SERIAL
      implicit none
      character(len=*), intent(in) :: msg
#else
      use mpi
      implicit none
      character(len=*), intent(in) :: msg
#ifdef D3Q27
      integer      :: reqs_up(9)
      integer      :: reqs_down(9)
#else
      integer      :: reqs_up(5)
      integer      :: reqs_down(5)
#endif
      
      integer      :: status_up(MPI_STATUS_SIZE)
      integer      :: status_down(MPI_STATUS_SIZE)
      integer      :: tag

#endif
      logical :: faiPops, faiRho

      if (withCG) then
#ifdef SERIAL
          if(xperiodic)call bc_per_x2<<<dimGridx, dimBlock2>>>(step, flipflop, faiPops,faiRho)
          if(yperiodic)call bc_per_y2<<<dimGridy, dimBlock2>>>(step, flipflop, faiPops,faiRho)
          if(zperiodic)call bc_per_z2<<<dimGridz, dimBlock2>>>(step, flipflop, faiPops,faiRho)

          call bc_edge_x2<<<(nx+2+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop, faiPops,faiRho)
          call bc_edge_y2<<<(ny+2+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop, faiPops,faiRho)
          call bc_edge_z2<<<(nz+2+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop, faiPops,faiRho)

          !call bc_corners<<<dimGridz, dimBlock2>>>(step, flipflop)
#else
          if(xperiodic)call bc_per_x2<<<dimGridx, dimBlock2>>>(step, flipflop, faiPops,faiRho)
          if(yperiodic)call bc_per_y2<<<dimGridy, dimBlock2>>>(step, flipflop, faiPops,faiRho)
          
          
          
          call abortOnLastErrorAndSync(msg,step, flipflop)

! if it is the same task, there's no need to send/recv data
          if(up(2) == myrank) then
             if(zperiodic)call bc_per_z2<<<dimGridz, dimBlock2>>>(step, flipflop, faiPops,faiRho)
          else 
            if (faiRho) then
! 
! I've to exchange pop with other task
! -----------------------------------------------------------------------------------
! send RED rho
! 1) isend to up task 
              tag = 405
              call mpi_isend(rhoR_d(0,0,  nz), 1, xyplane, up(2), tag, &
                                lbecomm, reqs_up(1), ierr)
!
! 2) recv from down task
              tag = 405
              call  mpi_recv(rhoR_d(0,0,   0), 1, xyplane, down(2), tag, &
                                lbecomm, status_up, ierr)
!
! 3) isend to down task
              tag = 406
              call mpi_isend(rhoR_d(0,0,   1), 1, xyplane, down(2), tag, &
                                lbecomm, reqs_down(1), ierr)
!
! 4) recv from up task
              tag = 406
              call  mpi_recv(rhoR_d(0,0,nz+1), 1, xyplane, up(2), tag, &
                                lbecomm, status_down, ierr)
!
! 5) syncronize
              call mpi_wait(reqs_up(1), status_up, ierr)
              call mpi_wait(reqs_down(1), status_down, ierr)
            endif
!
! -----------------------------------------------------------------------------------
! Send RED pops
! 1) isend to up task
            if (faiPops) then
              tag = 005
              call mpi_isend(popsR_d(0,0,nz  , 5,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, reqs_up(1), ierr)
              tag = 011
              call mpi_isend(popsR_d(0,0,nz  ,11,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, reqs_up(2), ierr)
              tag = 013
              call mpi_isend(popsR_d(0,0,nz  ,13,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, reqs_up(3), ierr)
              tag = 015
              call mpi_isend(popsR_d(0,0,nz  ,15,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, reqs_up(4), ierr)
              tag = 017
              call mpi_isend(popsR_d(0,0,nz  ,17,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, reqs_up(5), ierr)
#ifdef D3Q27
              tag = 019
              call mpi_isend(popsR_d(0,0,nz  ,19,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, reqs_up(6), ierr)
              tag = 022
              call mpi_isend(popsR_d(0,0,nz  ,22,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, reqs_up(7), ierr)
              tag = 023
              call mpi_isend(popsR_d(0,0,nz  ,23,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, reqs_up(8), ierr)
              tag = 025
              call mpi_isend(popsR_d(0,0,nz  ,25,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, reqs_up(9), ierr)
#endif

!
! 2) recv from down task
              tag = 005
              call mpi_recv(popsR_d(0,0,0  , 5,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 011
              call mpi_recv(popsR_d(0,0, 0  ,11,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 013
              call mpi_recv(popsR_d(0,0, 0  ,13,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 015
              call mpi_recv(popsR_d(0,0, 0  ,15,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 017
              call mpi_recv(popsR_d(0,0, 0  ,17,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, status_up, ierr)
#ifdef D3Q27
              tag = 019
              call mpi_recv(popsR_d(0,0, 0  ,19,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 022
              call mpi_recv(popsR_d(0,0, 0  ,22,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 023
              call mpi_recv(popsR_d(0,0, 0  ,23,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 025
              call mpi_recv(popsR_d(0,0, 0  ,25,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, status_up, ierr)
#endif
!
! 3) isend to down task
              tag = 006
              call mpi_isend(popsR_d(0,0, 1  , 6,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, reqs_down(1), ierr)
              tag = 012
              call mpi_isend(popsR_d(0,0, 1  ,12,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, reqs_down(2), ierr)
              tag = 014
              call mpi_isend(popsR_d(0,0, 1  ,14,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, reqs_down(3), ierr)
              tag = 016
              call mpi_isend(popsR_d(0,0, 1  ,16,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, reqs_down(4), ierr)
              tag = 018
              call mpi_isend(popsR_d(0,0, 1  ,18,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, reqs_down(5), ierr)
#ifdef D3Q27
              tag = 020
              call mpi_isend(popsR_d(0,0, 1  ,20,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, reqs_down(6), ierr)
              tag = 021
              call mpi_isend(popsR_d(0,0, 1  ,21,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, reqs_down(7), ierr)
              tag = 024
              call mpi_isend(popsR_d(0,0, 1  ,24,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, reqs_down(8), ierr)
              tag = 026
              call mpi_isend(popsR_d(0,0, 1  ,26,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, reqs_down(9), ierr)
#endif
!
! 4) recv from up task
              tag = 006
              call mpi_recv(popsR_d( 0,0,nz+1, 6,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 012
              call mpi_recv(popsR_d( 0,0,nz+1,12,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 014
              call mpi_recv(popsR_d( 0,0,nz+1,14,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 016
              call mpi_recv(popsR_d( 0,0,nz+1,16,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 018
              call mpi_recv(popsR_d( 0,0,nz+1,18,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, status_down, ierr)
#ifdef D3Q27
              tag = 020
              call mpi_recv(popsR_d( 0,0,nz+1,20,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 021
              call mpi_recv(popsR_d( 0,0,nz+1,21,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 024
              call mpi_recv(popsR_d( 0,0,nz+1,24,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 026
              call mpi_recv(popsR_d( 0,0,nz+1,26,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, status_down, ierr)
#endif
!
! 5) sincronizza
              call mpi_wait(reqs_up(1), status_up, ierr)
              call mpi_wait(reqs_down(1), status_down, ierr)
              call mpi_wait(reqs_up(2), status_up, ierr)
              call mpi_wait(reqs_down(2), status_down, ierr)
              call mpi_wait(reqs_up(3), status_up, ierr)
              call mpi_wait(reqs_down(3), status_down, ierr)
              call mpi_wait(reqs_up(4), status_up, ierr)
              call mpi_wait(reqs_down(4), status_down, ierr)
              call mpi_wait(reqs_up(5), status_up, ierr)
              call mpi_wait(reqs_down(5), status_down, ierr)
#ifdef D3Q27
              call mpi_wait(reqs_up(6), status_up, ierr)
              call mpi_wait(reqs_down(6), status_down, ierr)
              call mpi_wait(reqs_up(7), status_up, ierr)
              call mpi_wait(reqs_down(7), status_down, ierr)
              call mpi_wait(reqs_up(8), status_up, ierr)
              call mpi_wait(reqs_down(8), status_down, ierr)
              call mpi_wait(reqs_up(9), status_up, ierr)
              call mpi_wait(reqs_down(9), status_down, ierr)
#endif
            endif
!
! -----------------------------------------------------------------------------------
! send BLUE rho
! 1) isend to up task
           if (faiRho) then
            tag = 205 
            call mpi_isend(rhoB_d(0,0,  nz), 1, xyplane, up(2), tag, &
                              lbecomm, reqs_up(1), ierr)
!
! 2) recv from down task
            tag = 205
            call  mpi_recv(rhoB_d(0,0,   0), 1, xyplane, down(2), tag, &
                              lbecomm, status_up, ierr)
!
! 3) isend to down task
            tag = 206
            call mpi_isend(rhoB_d(0,0,   1), 1, xyplane, down(2), tag, &
                              lbecomm, reqs_down(1), ierr)
!
! 4) recv from up task
            tag = 206
            call  mpi_recv(rhoB_d(0,0,nz+1), 1, xyplane, up(2), tag, &
                              lbecomm, status_down, ierr)
!
! 5) sincronizza
            call mpi_wait(reqs_up(1), status_up, ierr)
            call mpi_wait(reqs_down(1), status_down, ierr)
           endif
!
! -----------------------------------------------------------------------------------
! send BLUE popos
! 1) isend to up task
            if (faiPops) then
            tag = 105
            call mpi_isend(popsB_d(0,0,nz  , 5,flipflop), 1, xyplane, up(2), tag, &
                              lbecomm, reqs_up(1), ierr)
            tag = 111
            call mpi_isend(popsB_d(0,0,nz  ,11,flipflop), 1, xyplane, up(2), tag, &
                              lbecomm, reqs_up(2), ierr)
            tag = 113
            call mpi_isend(popsB_d(0,0,nz  ,13,flipflop), 1, xyplane, up(2), tag, &
                              lbecomm, reqs_up(3), ierr)
            tag = 115
            call mpi_isend(popsB_d(0,0,nz  ,15,flipflop), 1, xyplane, up(2), tag, &
                              lbecomm, reqs_up(4), ierr)
            tag = 117
            call mpi_isend(popsB_d(0,0,nz  ,17,flipflop), 1, xyplane, up(2), tag, &
                              lbecomm, reqs_up(5), ierr)
#ifdef D3Q27
            tag = 119
            call mpi_isend(popsB_d(0,0,nz  ,19,flipflop), 1, xyplane, up(2), tag, &
                              lbecomm, reqs_up(6), ierr)
            tag = 122
            call mpi_isend(popsB_d(0,0,nz  ,22,flipflop), 1, xyplane, up(2), tag, &
                              lbecomm, reqs_up(7), ierr)
            tag = 123
            call mpi_isend(popsB_d(0,0,nz  ,23,flipflop), 1, xyplane, up(2), tag, &
                              lbecomm, reqs_up(8), ierr)
            tag = 125
            call mpi_isend(popsB_d(0,0,nz  ,25,flipflop), 1, xyplane, up(2), tag, &
                              lbecomm, reqs_up(9), ierr)
#endif
!
! 2) recv from down task
            tag = 105
            call mpi_recv(popsB_d( 0,0,0  , 5,flipflop), 1, xyplane, down(2), tag, &
                              lbecomm, status_up, ierr)
            tag = 111
            call mpi_recv(popsB_d(0,0, 0  ,11,flipflop), 1, xyplane, down(2), tag, &
                              lbecomm, status_up, ierr)
            tag = 113
            call mpi_recv(popsB_d(0,0, 0  ,13,flipflop), 1, xyplane, down(2), tag, &
                              lbecomm, status_up, ierr)
            tag = 115
            call mpi_recv(popsB_d(0,0, 0  ,15,flipflop), 1, xyplane, down(2), tag, &
                              lbecomm, status_up, ierr)
            tag = 117
            call mpi_recv(popsB_d(0,0, 0  ,17,flipflop), 1, xyplane, down(2), tag, &
                              lbecomm, status_up, ierr)
#ifdef D3Q27
            tag = 119
            call mpi_recv(popsB_d(0,0, 0  ,19,flipflop), 1, xyplane, down(2), tag, &
                              lbecomm, status_up, ierr)
            tag = 122
            call mpi_recv(popsB_d(0,0, 0  ,22,flipflop), 1, xyplane, down(2), tag, &
                              lbecomm, status_up, ierr)
            tag = 123
            call mpi_recv(popsB_d(0,0, 0  ,23,flipflop), 1, xyplane, down(2), tag, &
                              lbecomm, status_up, ierr)
            tag = 125
            call mpi_recv(popsB_d(0,0, 0  ,25,flipflop), 1, xyplane, down(2), tag, &
                              lbecomm, status_up, ierr)
#endif
!
! 3) isend to down task
            tag = 106
            call mpi_isend(popsB_d(0,0, 1  , 6,flipflop), 1, xyplane, down(2), tag, &
                              lbecomm, reqs_down(1), ierr)
            tag = 112
            call mpi_isend(popsB_d(0,0, 1  ,12,flipflop), 1, xyplane, down(2), tag, &
                              lbecomm, reqs_down(2), ierr)
            tag = 114
            call mpi_isend(popsB_d(0,0, 1  ,14,flipflop), 1, xyplane, down(2), tag, &
                              lbecomm, reqs_down(3), ierr)
            tag = 116
            call mpi_isend(popsB_d(0,0, 1  ,16,flipflop), 1, xyplane, down(2), tag, &
                              lbecomm, reqs_down(4), ierr)
            tag = 118
            call mpi_isend(popsB_d(0,0, 1  ,18,flipflop), 1, xyplane, down(2), tag, &
                              lbecomm, reqs_down(5), ierr)
#ifdef D3Q27
            tag = 120
            call mpi_isend(popsB_d(0,0, 1  ,20,flipflop), 1, xyplane, down(2), tag, &
                              lbecomm, reqs_down(6), ierr)
            tag = 121
            call mpi_isend(popsB_d(0,0, 1  ,21,flipflop), 1, xyplane, down(2), tag, &
                              lbecomm, reqs_down(7), ierr)
            tag = 124
            call mpi_isend(popsB_d(0,0, 1  ,24,flipflop), 1, xyplane, down(2), tag, &
                              lbecomm, reqs_down(8), ierr)
            tag = 126
            call mpi_isend(popsB_d(0,0, 1  ,26,flipflop), 1, xyplane, down(2), tag, &
                              lbecomm, reqs_down(9), ierr)
#endif
!
! 4) recv from up task
            tag = 106
            call mpi_recv(popsB_d( 0,0,nz+1, 6,flipflop), 1, xyplane, up(2), tag, &
                              lbecomm, status_down, ierr)
            tag = 112
            call mpi_recv(popsB_d( 0,0,nz+1,12,flipflop), 1, xyplane, up(2), tag, &
                              lbecomm, status_down, ierr)
            tag = 114
            call mpi_recv(popsB_d( 0,0,nz+1,14,flipflop), 1, xyplane, up(2), tag, &
                              lbecomm, status_down, ierr)
            tag = 116
            call mpi_recv(popsB_d( 0,0,nz+1,16,flipflop), 1, xyplane, up(2), tag, &
                              lbecomm, status_down, ierr)
            tag = 118
            call mpi_recv(popsB_d( 0,0,nz+1,18,flipflop), 1, xyplane, up(2), tag, &
                              lbecomm, status_down, ierr)
#ifdef D3Q27
            tag = 120
            call mpi_recv(popsB_d( 0,0,nz+1,20,flipflop), 1, xyplane, up(2), tag, &
                              lbecomm, status_down, ierr)
            tag = 121
            call mpi_recv(popsB_d( 0,0,nz+1,21,flipflop), 1, xyplane, up(2), tag, &
                              lbecomm, status_down, ierr)
            tag = 124
            call mpi_recv(popsB_d( 0,0,nz+1,24,flipflop), 1, xyplane, up(2), tag, &
                              lbecomm, status_down, ierr)
            tag = 126
            call mpi_recv(popsB_d( 0,0,nz+1,26,flipflop), 1, xyplane, up(2), tag, &
                              lbecomm, status_down, ierr)
#endif
!
! 5) sincronizza
            call mpi_wait(reqs_up(1), status_up, ierr)
            call mpi_wait(reqs_down(1), status_down, ierr)
            call mpi_wait(reqs_up(2), status_up, ierr)
            call mpi_wait(reqs_down(2), status_down, ierr)
            call mpi_wait(reqs_up(3), status_up, ierr)
            call mpi_wait(reqs_down(3), status_down, ierr)
            call mpi_wait(reqs_up(4), status_up, ierr)
            call mpi_wait(reqs_down(4), status_down, ierr)
            call mpi_wait(reqs_up(5), status_up, ierr)
            call mpi_wait(reqs_down(5), status_down, ierr)
#ifdef D3Q27
            call mpi_wait(reqs_up(6), status_up, ierr)
            call mpi_wait(reqs_down(6), status_down, ierr)
            call mpi_wait(reqs_up(7), status_up, ierr)
            call mpi_wait(reqs_down(7), status_down, ierr)
            call mpi_wait(reqs_up(8), status_up, ierr)
            call mpi_wait(reqs_down(8), status_down, ierr)
            call mpi_wait(reqs_up(9), status_up, ierr)
            call mpi_wait(reqs_down(9), status_down, ierr)
#endif
            endif

          endif

          if(up(2) == myrank) then
            if(xperiodic)&
             call bc_edge_x2<<<(nx+2+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop, faiPops,faiRho)
            if(yperiodic)&
             call bc_edge_y2<<<(ny+2+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop, faiPops,faiRho)  
            if(zperiodic)&
             call bc_edge_z2<<<(nz+2+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop, faiPops,faiRho)
          else
            ! NON VA FATTO
              !call abortOnLastErrorAndSync(msg,step, flipflop)
              
!             mpi rho+pop myrank <--> up(2)
             if(zperiodic)&
              call bc_edge_z2<<<(nz+2+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop, faiPops,faiRho)
          endif

          !call bc_corners<<<dimGridz, dimBlock2>>>(step, flipflop)

#endif
      
      
      endif
      
      
      
    end subroutine fixPeriodic
    
        subroutine fixPeriodic_BGK(faiPops,faiRho, msg)
#ifdef SERIAL
      implicit none
      character(len=*), intent(in) :: msg
#else
      use mpi
      implicit none
      character(len=*), intent(in) :: msg
#ifdef D3Q27
      integer      :: reqs_up(9)
      integer      :: reqs_down(9)
#else
      integer      :: reqs_up(5)
      integer      :: reqs_down(5)
#endif
      
      integer      :: status_up(MPI_STATUS_SIZE)
      integer      :: status_down(MPI_STATUS_SIZE)
      integer      :: tag

#endif
      logical :: faiPops, faiRho

      

          if(xperiodic)call bc_BGK_per_x<<<dimGridx, dimBlock2>>>(step, flipflop, faiPops,faiRho)
          if(yperiodic)call bc_BGK_per_y<<<dimGridy, dimBlock2>>>(step, flipflop, faiPops,faiRho)
          
          
          call abortOnLastErrorAndSync(msg,step, flipflop)

! if it is the same task, there's no need to send/recv data
          if(up(2) == myrank) then
             if(zperiodic)call bc_BGK_per_z<<<dimGridz, dimBlock2>>>(step, flipflop, faiPops,faiRho)
          else 
            if (faiRho) then
! 
! I've to exchange pop with other task
! -----------------------------------------------------------------------------------
! send RED rho
! 1) isend to up task 
              tag = 405
              call mpi_isend(rhoR_d(0,0,  nz), 1, xyplane, up(2), tag, &
                                lbecomm, reqs_up(1), ierr)
!
! 2) recv from down task
              tag = 405
              call  mpi_recv(rhoR_d(0,0,   0), 1, xyplane, down(2), tag, &
                                lbecomm, status_up, ierr)
!
! 3) isend to down task
              tag = 406
              call mpi_isend(rhoR_d(0,0,   1), 1, xyplane, down(2), tag, &
                                lbecomm, reqs_down(1), ierr)
!
! 4) recv from up task
              tag = 406
              call  mpi_recv(rhoR_d(0,0,nz+1), 1, xyplane, up(2), tag, &
                                lbecomm, status_down, ierr)
!
! 5) syncronize
              call mpi_wait(reqs_up(1), status_up, ierr)
              call mpi_wait(reqs_down(1), status_down, ierr)
            endif
!
! -----------------------------------------------------------------------------------
! Send RED pops
! 1) isend to up task
            if (faiPops) then
              tag = 005
              call mpi_isend(popsR_d(0,0,nz  , 5,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, reqs_up(1), ierr)
              tag = 011
              call mpi_isend(popsR_d(0,0,nz  ,11,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, reqs_up(2), ierr)
              tag = 013
              call mpi_isend(popsR_d(0,0,nz  ,13,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, reqs_up(3), ierr)
              tag = 015
              call mpi_isend(popsR_d(0,0,nz  ,15,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, reqs_up(4), ierr)
              tag = 017
              call mpi_isend(popsR_d(0,0,nz  ,17,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, reqs_up(5), ierr)
#ifdef D3Q27
              tag = 019
              call mpi_isend(popsR_d(0,0,nz  ,19,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, reqs_up(6), ierr)
              tag = 022
              call mpi_isend(popsR_d(0,0,nz  ,22,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, reqs_up(7), ierr)
              tag = 023
              call mpi_isend(popsR_d(0,0,nz  ,23,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, reqs_up(8), ierr)
              tag = 025
              call mpi_isend(popsR_d(0,0,nz  ,25,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, reqs_up(9), ierr)
#endif

!
! 2) recv from down task
              tag = 005
              call mpi_recv(popsR_d(0,0,0  , 5,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 011
              call mpi_recv(popsR_d(0,0, 0  ,11,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 013
              call mpi_recv(popsR_d(0,0, 0  ,13,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 015
              call mpi_recv(popsR_d(0,0, 0  ,15,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 017
              call mpi_recv(popsR_d(0,0, 0  ,17,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, status_up, ierr)
#ifdef D3Q27
              tag = 019
              call mpi_recv(popsR_d(0,0, 0  ,19,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 022
              call mpi_recv(popsR_d(0,0, 0  ,22,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 023
              call mpi_recv(popsR_d(0,0, 0  ,23,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 025
              call mpi_recv(popsR_d(0,0, 0  ,25,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, status_up, ierr)
#endif
!
! 3) isend to down task
              tag = 006
              call mpi_isend(popsR_d(0,0, 1  , 6,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, reqs_down(1), ierr)
              tag = 012
              call mpi_isend(popsR_d(0,0, 1  ,12,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, reqs_down(2), ierr)
              tag = 014
              call mpi_isend(popsR_d(0,0, 1  ,14,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, reqs_down(3), ierr)
              tag = 016
              call mpi_isend(popsR_d(0,0, 1  ,16,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, reqs_down(4), ierr)
              tag = 018
              call mpi_isend(popsR_d(0,0, 1  ,18,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, reqs_down(5), ierr)
#ifdef D3Q27
              tag = 020
              call mpi_isend(popsR_d(0,0, 1  ,20,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, reqs_down(6), ierr)
              tag = 021
              call mpi_isend(popsR_d(0,0, 1  ,21,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, reqs_down(7), ierr)
              tag = 024
              call mpi_isend(popsR_d(0,0, 1  ,24,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, reqs_down(8), ierr)
              tag = 026
              call mpi_isend(popsR_d(0,0, 1  ,26,flipflop), 1, xyplane, down(2), tag, &
                                lbecomm, reqs_down(9), ierr)
#endif
!
! 4) recv from up task
              tag = 006
              call mpi_recv(popsR_d( 0,0,nz+1, 6,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 012
              call mpi_recv(popsR_d( 0,0,nz+1,12,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 014
              call mpi_recv(popsR_d( 0,0,nz+1,14,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 016
              call mpi_recv(popsR_d( 0,0,nz+1,16,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 018
              call mpi_recv(popsR_d( 0,0,nz+1,18,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, status_down, ierr)
#ifdef D3Q27
              tag = 020
              call mpi_recv(popsR_d( 0,0,nz+1,20,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 021
              call mpi_recv(popsR_d( 0,0,nz+1,21,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 024
              call mpi_recv(popsR_d( 0,0,nz+1,24,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 026
              call mpi_recv(popsR_d( 0,0,nz+1,26,flipflop), 1, xyplane, up(2), tag, &
                                lbecomm, status_down, ierr)
#endif
!
! 5) sincronizza
              call mpi_wait(reqs_up(1), status_up, ierr)
              call mpi_wait(reqs_down(1), status_down, ierr)
              call mpi_wait(reqs_up(2), status_up, ierr)
              call mpi_wait(reqs_down(2), status_down, ierr)
              call mpi_wait(reqs_up(3), status_up, ierr)
              call mpi_wait(reqs_down(3), status_down, ierr)
              call mpi_wait(reqs_up(4), status_up, ierr)
              call mpi_wait(reqs_down(4), status_down, ierr)
              call mpi_wait(reqs_up(5), status_up, ierr)
              call mpi_wait(reqs_down(5), status_down, ierr)
#ifdef D3Q27
              call mpi_wait(reqs_up(6), status_up, ierr)
              call mpi_wait(reqs_down(6), status_down, ierr)
              call mpi_wait(reqs_up(7), status_up, ierr)
              call mpi_wait(reqs_down(7), status_down, ierr)
              call mpi_wait(reqs_up(8), status_up, ierr)
              call mpi_wait(reqs_down(8), status_down, ierr)
              call mpi_wait(reqs_up(9), status_up, ierr)
              call mpi_wait(reqs_down(9), status_down, ierr)
#endif
            endif
!
! -----------------------------------------------------------------------------------



          endif

          if(up(2) == myrank) then
            if(xperiodic)&
             call bc_BGK_edge_x<<<(nx+2+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop, faiPops,faiRho)
            if(yperiodic)&
             call bc_BGK_edge_y<<<(ny+2+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop, faiPops,faiRho)  
            if(zperiodic)&
             call bc_BGK_edge_z<<<(nz+2+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop, faiPops,faiRho)
          else
            ! NON VA FATTO
              !call abortOnLastErrorAndSync(msg,step, flipflop)
              
!             mpi rho+pop myrank <--> up(2)
             if(zperiodic)&
              call bc_BGK_edge_z<<<(nz+2+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop, faiPops,faiRho)
          endif

          !call bc_corners<<<dimGridz, dimBlock2>>>(step, flipflop)


      
      
    
      
      
      
    end subroutine fixPeriodic_BGK
    
    subroutine fixPeriodic_hvar(faiRho,faiRhoB,faiVelsub, msg)
#ifdef SERIAL
      implicit none
      character(len=*), intent(in) :: msg
#else
      use mpi
      implicit none
      character(len=*), intent(in) :: msg

      integer      :: reqs_up(5)
      integer      :: reqs_down(5)

      
      integer      :: status_up(MPI_STATUS_SIZE)
      integer      :: status_down(MPI_STATUS_SIZE)
      integer      :: tag

#endif
      logical :: faiRho,faiRhoB,faiVelsub
      logical :: faivel
          
          
          
          faivel = (faivelsub .and. store_vel)
          
         if(xperiodic)&
          call bc_per_x_hvar<<<dimGridsidex, dimBlock2>>>(step,faiRho,faiRhoB,faiVel)
         if(yperiodic)&
          call bc_per_y_hvar<<<dimGridsidey, dimBlock2>>>(step,faiRho,faiRhoB,faiVel)
          
          
          call abortOnLastErrorAndSync(msg,step, flipflop)

! if it is the same task, there's no need to send/recv data
          if(up(2) == myrank) then
           if(zperiodic)&
             call bc_per_z_hvar<<<dimGridsidez, dimBlock2>>>(step,faiRho,faiRhoB,faiVel)
          else 
            if (faiRho) then
! 
! I've to exchange pop with other task
! -----------------------------------------------------------------------------------
! send RED rho
! 1) isend to up task 
              tag = 405
              call mpi_isend(rhoR_d(1-nbuff,1-nbuff,nz+1-nbuff), 1, xyplanebuff, up(2), tag, &
                                lbecomm, reqs_up(1), ierr)
!
! 2) recv from down task
              tag = 405
              call  mpi_recv(rhoR_d(1-nbuff,1-nbuff,1-nbuff), 1, xyplanebuff, down(2), tag, &
                                lbecomm, status_up, ierr)
!
! 3) isend to down task
              tag = 406
              call mpi_isend(rhoR_d(1-nbuff,1-nbuff, 1), 1, xyplanebuff, down(2), tag, &
                                lbecomm, reqs_down(1), ierr)
!
! 4) recv from up task
              tag = 406
              call  mpi_recv(rhoR_d(1-nbuff,1-nbuff,nz+1), 1, xyplanebuff, up(2), tag, &
                                lbecomm, status_down, ierr)
!
! 5) syncronize
              call mpi_wait(reqs_up(1), status_up, ierr)
              call mpi_wait(reqs_down(1), status_down, ierr)
            endif
            if (faiRhoB) then
! 
! I've to exchange pop with other task
! -----------------------------------------------------------------------------------
! send BLUE rho
! 1) isend to up task 
              tag = 205
              call mpi_isend(rhoB_d(1-nbuff,1-nbuff,nz+1-nbuff), 1, xyplanebuff, up(2), tag, &
                                lbecomm, reqs_up(1), ierr)
!
! 2) recv from down task
              tag = 205
              call  mpi_recv(rhoB_d(1-nbuff,1-nbuff,1-nbuff), 1, xyplanebuff, down(2), tag, &
                                lbecomm, status_up, ierr)
!
! 3) isend to down task
              tag = 206
              call mpi_isend(rhoB_d(1-nbuff,1-nbuff, 1), 1, xyplanebuff, down(2), tag, &
                                lbecomm, reqs_down(1), ierr)
!
! 4) recv from up task
              tag = 206
              call  mpi_recv(rhoB_d(1-nbuff,1-nbuff,nz+1), 1, xyplanebuff, up(2), tag, &
                                lbecomm, status_down, ierr)
!
! 5) syncronize
              call mpi_wait(reqs_up(1), status_up, ierr)
              call mpi_wait(reqs_down(1), status_down, ierr)
            endif
            
            if (faiVel) then
! 
! I've to exchange pop with other task
! -----------------------------------------------------------------------------------
! send velocity rho
! 1) isend to up task 
              tag = 605
              call mpi_isend(vel_d(1,1-nbuff,1-nbuff,nz+1-nbuff), 1, xyplane3, up(2), tag, &
                                lbecomm, reqs_up(1), ierr)
!
! 2) recv from down task
              tag = 605
              call  mpi_recv(vel_d(1,1-nbuff,1-nbuff,1-nbuff), 1, xyplane3, down(2), tag, &
                                lbecomm, status_up, ierr)
!
! 3) isend to down task
              tag = 606
              call mpi_isend(vel_d(1,1-nbuff,1-nbuff, 1), 1, xyplane3, down(2), tag, &
                                lbecomm, reqs_down(1), ierr)
!
! 4) recv from up task
              tag = 606
              call  mpi_recv(vel_d(1,1-nbuff,1-nbuff,nz+1), 1, xyplane3, up(2), tag, &
                                lbecomm, status_down, ierr)
!
! 5) syncronize
              call mpi_wait(reqs_up(1), status_up, ierr)
              call mpi_wait(reqs_down(1), status_down, ierr)
            endif

          endif

          if(up(2) == myrank) then
           if(xperiodic)&
            call bc_edge_x_hvar<<<(nx+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step,faiRho,faiRhoB,faiVel)
           if(yperiodic)&
            call bc_edge_y_hvar<<<(ny+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step,faiRho,faiRhoB,faiVel)
           if(zperiodic)&
            call bc_edge_z_hvar<<<(nz+2*nbuff+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step,faiRho,faiRhoB,faiVel)
          else
            ! NON VA FATTO
              !call abortOnLastErrorAndSync(msg,step, flipflop)
              
!             mpi rho+pop myrank <--> up(2)
            if(zperiodic)&
             call bc_edge_z_hvar<<<(nz+2*nbuff+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step,faiRho,faiRhoB,faiVel)
          endif

          

      
    end subroutine fixPeriodic_hvar
    
    subroutine fixPeriodic_SCP(faiPops,faiRho, msg)
#ifdef SERIAL
      implicit none
      character(len=*), intent(in) :: msg
#else
      use mpi
      implicit none
      character(len=*), intent(in) :: msg
#ifdef D3Q27
      integer      :: reqs_up(9)
      integer      :: reqs_down(9)
#else
      integer      :: reqs_up(5)
      integer      :: reqs_down(5)
#endif
      
      integer      :: status_up(MPI_STATUS_SIZE)
      integer      :: status_down(MPI_STATUS_SIZE)
      integer      :: tag

#endif
      logical :: faiPops, faiRho

      

          if(xperiodic)call bc_SCP_per_x<<<dimGridx, dimBlock2>>>(step, flipflop, faiPops,faiRho)
          if(yperiodic)call bc_SCP_per_y<<<dimGridy, dimBlock2>>>(step, flipflop, faiPops,faiRho)
          
          
          call abortOnLastErrorAndSync(msg,step, flipflop)

! if it is the same task, there's no need to send/recv data
          if(up(2) == myrank) then
             if(zperiodic)call bc_SCP_per_z<<<dimGridz, dimBlock2>>>(step, flipflop, faiPops,faiRho)
          else 
            if (faiRho) then
! 
! I've to exchange pop with other task
! -----------------------------------------------------------------------------------
! send RED rho
! 1) isend to up task 
              tag = 405
              call mpi_isend(scalar_d(1,0,0,  nz), 1, xyplanescp, up(2), tag, &
                                lbecomm, reqs_up(1), ierr)
!
! 2) recv from down task
              tag = 405
              call  mpi_recv(scalar_d(1,0,0,   0), 1, xyplanescp, down(2), tag, &
                                lbecomm, status_up, ierr)
!
! 3) isend to down task
              tag = 406
              call mpi_isend(scalar_d(1,0,0,   1), 1, xyplanescp, down(2), tag, &
                                lbecomm, reqs_down(1), ierr)
!
! 4) recv from up task
              tag = 406
              call  mpi_recv(scalar_d(1,0,0,nz+1), 1, xyplanescp, up(2), tag, &
                                lbecomm, status_down, ierr)
!
! 5) syncronize
              call mpi_wait(reqs_up(1), status_up, ierr)
              call mpi_wait(reqs_down(1), status_down, ierr)
            endif
!
! -----------------------------------------------------------------------------------
! Send RED pops
! 1) isend to up task
            if (faiPops) then
#ifdef SCPD3Q27
              tag = 005
              call mpi_isend(popsSCP_d(1,0,0,nz  , 5,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, reqs_up(1), ierr)
              tag = 011
              call mpi_isend(popsSCP_d(1,0,0,nz  ,11,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, reqs_up(2), ierr)
              tag = 013
              call mpi_isend(popsSCP_d(1,0,0,nz  ,13,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, reqs_up(3), ierr)
              tag = 015
              call mpi_isend(popsSCP_d(1,0,0,nz  ,15,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, reqs_up(4), ierr)
              tag = 017
              call mpi_isend(popsSCP_d(1,0,0,nz  ,17,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, reqs_up(5), ierr)
              tag = 019
              call mpi_isend(popsSCP_d(1,0,0,nz  ,19,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, reqs_up(6), ierr)
              tag = 022
              call mpi_isend(popsSCP_d(1,0,0,nz  ,22,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, reqs_up(7), ierr)
              tag = 023
              call mpi_isend(popsSCP_d(1,0,0,nz  ,23,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, reqs_up(8), ierr)
              tag = 025
              call mpi_isend(popsSCP_d(1,0,0,nz  ,25,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, reqs_up(9), ierr)
#elif defined SCPD3Q19
              tag = 005
              call mpi_isend(popsSCP_d(1,0,0,nz  , 5,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, reqs_up(1), ierr)
              tag = 011
              call mpi_isend(popsSCP_d(1,0,0,nz  ,11,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, reqs_up(2), ierr)
              tag = 013
              call mpi_isend(popsSCP_d(1,0,0,nz  ,13,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, reqs_up(3), ierr)
              tag = 015
              call mpi_isend(popsSCP_d(1,0,0,nz  ,15,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, reqs_up(4), ierr)
              tag = 017
              call mpi_isend(popsSCP_d(1,0,0,nz  ,17,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, reqs_up(5), ierr)
#else
              tag = 005
              call mpi_isend(popsSCP_d(1,0,0,nz  , 5,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, reqs_up(1), ierr)
#endif

!
! 2) recv from down task
#ifdef SCPD3Q27
              tag = 005
              call mpi_recv(popsSCP_d(1,0,0,0  , 5,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 011
              call mpi_recv(popsSCP_d(1,0,0, 0  ,11,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 013
              call mpi_recv(popsSCP_d(1,0,0, 0  ,13,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 015
              call mpi_recv(popsSCP_d(1,0,0, 0  ,15,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 017
              call mpi_recv(popsSCP_d(1,0,0, 0  ,17,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 019
              call mpi_recv(popsSCP_d(1,0,0, 0  ,19,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 022
              call mpi_recv(popsSCP_d(1,0,0, 0  ,22,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 023
              call mpi_recv(popsSCP_d(1,0,0, 0  ,23,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 025
              call mpi_recv(popsSCP_d(1,0,0, 0  ,25,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, status_up, ierr)
#elif defined SCPD3Q19
              tag = 005
              call mpi_recv(popsSCP_d(1,0,0,0  , 5,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 011
              call mpi_recv(popsSCP_d(1,0,0, 0  ,11,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 013
              call mpi_recv(popsSCP_d(1,0,0, 0  ,13,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 015
              call mpi_recv(popsSCP_d(1,0,0, 0  ,15,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, status_up, ierr)
              tag = 017
              call mpi_recv(popsSCP_d(1,0,0, 0  ,17,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, status_up, ierr)
#else
              tag = 005
              call mpi_recv(popsSCP_d(1,0,0,0  , 5,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, status_up, ierr)
#endif
!
! 3) isend to down task
#ifdef SCPD3Q27
              tag = 006
              call mpi_isend(popsSCP_d(1,0,0, 1  , 6,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, reqs_down(1), ierr)
              tag = 012
              call mpi_isend(popsSCP_d(1,0,0, 1  ,12,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, reqs_down(2), ierr)
              tag = 014
              call mpi_isend(popsSCP_d(1,0,0, 1  ,14,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, reqs_down(3), ierr)
              tag = 016
              call mpi_isend(popsSCP_d(1,0,0, 1  ,16,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, reqs_down(4), ierr)
              tag = 018
              call mpi_isend(popsSCP_d(1,0,0, 1  ,18,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, reqs_down(5), ierr)
              tag = 020
              call mpi_isend(popsSCP_d(1,0,0, 1  ,20,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, reqs_down(6), ierr)
              tag = 021
              call mpi_isend(popsSCP_d(1,0,0, 1  ,21,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, reqs_down(7), ierr)
              tag = 024
              call mpi_isend(popsSCP_d(1,0,0, 1  ,24,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, reqs_down(8), ierr)
              tag = 026
              call mpi_isend(popsSCP_d(1,0,0, 1  ,26,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, reqs_down(9), ierr)
#elif defined SCPD3Q19
              tag = 006
              call mpi_isend(popsSCP_d(1,0,0, 1  , 6,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, reqs_down(1), ierr)
              tag = 012
              call mpi_isend(popsSCP_d(1,0,0, 1  ,12,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, reqs_down(2), ierr)
              tag = 014
              call mpi_isend(popsSCP_d(1,0,0, 1  ,14,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, reqs_down(3), ierr)
              tag = 016
              call mpi_isend(popsSCP_d(1,0,0, 1  ,16,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, reqs_down(4), ierr)
              tag = 018
              call mpi_isend(popsSCP_d(1,0,0, 1  ,18,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, reqs_down(5), ierr)
#else
              tag = 006
              call mpi_isend(popsSCP_d(1,0,0, 1  , 6,flipflop), 1, xyplanescp, down(2), tag, &
                                lbecomm, reqs_down(1), ierr)
#endif
!
! 4) recv from up task
#ifdef SCPD3Q27
              tag = 006
              call mpi_recv(popsSCP_d(1, 0,0,nz+1, 6,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 012
              call mpi_recv(popsSCP_d(1, 0,0,nz+1,12,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 014
              call mpi_recv(popsSCP_d(1, 0,0,nz+1,14,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 016
              call mpi_recv(popsSCP_d(1, 0,0,nz+1,16,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 018
              call mpi_recv(popsSCP_d(1, 0,0,nz+1,18,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 020
              call mpi_recv(popsSCP_d(1, 0,0,nz+1,20,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 021
              call mpi_recv(popsSCP_d(1, 0,0,nz+1,21,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 024
              call mpi_recv(popsSCP_d(1, 0,0,nz+1,24,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 026
              call mpi_recv(popsSCP_d(1, 0,0,nz+1,26,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, status_down, ierr)
#elif defined SCPD3Q19
              tag = 006
              call mpi_recv(popsSCP_d(1, 0,0,nz+1, 6,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 012
              call mpi_recv(popsSCP_d(1, 0,0,nz+1,12,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 014
              call mpi_recv(popsSCP_d(1, 0,0,nz+1,14,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 016
              call mpi_recv(popsSCP_d(1, 0,0,nz+1,16,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, status_down, ierr)
              tag = 018
              call mpi_recv(popsSCP_d(1, 0,0,nz+1,18,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, status_down, ierr)
#else
              tag = 006
              call mpi_recv(popsSCP_d(1, 0,0,nz+1, 6,flipflop), 1, xyplanescp, up(2), tag, &
                                lbecomm, status_down, ierr)
#endif
!
! 5) sincronizza
#ifdef SCPD3Q27
              call mpi_wait(reqs_up(1), status_up, ierr)
              call mpi_wait(reqs_down(1), status_down, ierr)
              call mpi_wait(reqs_up(2), status_up, ierr)
              call mpi_wait(reqs_down(2), status_down, ierr)
              call mpi_wait(reqs_up(3), status_up, ierr)
              call mpi_wait(reqs_down(3), status_down, ierr)
              call mpi_wait(reqs_up(4), status_up, ierr)
              call mpi_wait(reqs_down(4), status_down, ierr)
              call mpi_wait(reqs_up(5), status_up, ierr)
              call mpi_wait(reqs_down(5), status_down, ierr)
              call mpi_wait(reqs_up(6), status_up, ierr)
              call mpi_wait(reqs_down(6), status_down, ierr)
              call mpi_wait(reqs_up(7), status_up, ierr)
              call mpi_wait(reqs_down(7), status_down, ierr)
              call mpi_wait(reqs_up(8), status_up, ierr)
              call mpi_wait(reqs_down(8), status_down, ierr)
              call mpi_wait(reqs_up(9), status_up, ierr)
              call mpi_wait(reqs_down(9), status_down, ierr)
#elif defined SCPD3Q19
              call mpi_wait(reqs_up(1), status_up, ierr)
              call mpi_wait(reqs_down(1), status_down, ierr)
              call mpi_wait(reqs_up(2), status_up, ierr)
              call mpi_wait(reqs_down(2), status_down, ierr)
              call mpi_wait(reqs_up(3), status_up, ierr)
              call mpi_wait(reqs_down(3), status_down, ierr)
              call mpi_wait(reqs_up(4), status_up, ierr)
              call mpi_wait(reqs_down(4), status_down, ierr)
              call mpi_wait(reqs_up(5), status_up, ierr)
              call mpi_wait(reqs_down(5), status_down, ierr)
#else
              call mpi_wait(reqs_up(1), status_up, ierr)
              call mpi_wait(reqs_down(1), status_down, ierr)
#endif
            endif
!
! -----------------------------------------------------------------------------------



          endif

          if(up(2) == myrank) then
            if(xperiodic)&
             call bc_SCP_edge_x<<<(nx+2+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop, faiPops,faiRho)
            if(yperiodic)&
             call bc_SCP_edge_y<<<(ny+2+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop, faiPops,faiRho)  
            if(zperiodic)&
             call bc_SCP_edge_z<<<(nz+2+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop, faiPops,faiRho)
          else
            ! NON VA FATTO
              !call abortOnLastErrorAndSync(msg,step, flipflop)
              
!             mpi rho+pop myrank <--> up(2)
             if(zperiodic)&
              call bc_SCP_edge_z<<<(nz+2+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step, flipflop, faiPops,faiRho)
          endif      
      
    end subroutine fixPeriodic_SCP
    
    subroutine fixPeriodic_svar(faiRho, msg)
#ifdef SERIAL
      implicit none
      character(len=*), intent(in) :: msg
#else
      use mpi
      implicit none
      character(len=*), intent(in) :: msg

      integer      :: reqs_up(5)
      integer      :: reqs_down(5)

      
      integer      :: status_up(MPI_STATUS_SIZE)
      integer      :: status_down(MPI_STATUS_SIZE)
      integer      :: tag

#endif
      logical :: faiRho

          
          
          
  
          
         if(xperiodic)&
          call bc_per_x_svar<<<dimGridsidex, dimBlock2>>>(step,faiRho)
         if(yperiodic)&
          call bc_per_y_svar<<<dimGridsidey, dimBlock2>>>(step,faiRho)
          
          
          call abortOnLastErrorAndSync(msg,step, flipflop)

! if it is the same task, there's no need to send/recv data
          if(up(2) == myrank) then
            if(zperiodic)&
             call bc_per_z_svar<<<dimGridsidez, dimBlock2>>>(step,faiRho)
          else 
            if (faiRho) then

! 
! I've to exchange pop with other task
! -----------------------------------------------------------------------------------
! send BLUE rho
! 1) isend to up task 
              tag = 805
              call mpi_isend(scalar_d(1,1-nbuff,1-nbuff,nz+1-nbuff), 1, xyplanesvar, up(2), tag, &
                                lbecomm, reqs_up(1), ierr)
!
! 2) recv from down task
              tag = 805
              call  mpi_recv(scalar_d(1,1-nbuff,1-nbuff,1-nbuff), 1, xyplanesvar, down(2), tag, &
                                lbecomm, status_up, ierr)
!
! 3) isend to down task
              tag = 806
              call mpi_isend(scalar_d(1,1-nbuff,1-nbuff, 1), 1, xyplanesvar, down(2), tag, &
                                lbecomm, reqs_down(1), ierr)
!
! 4) recv from up task
              tag = 806
              call  mpi_recv(scalar_d(1,1-nbuff,1-nbuff,nz+1), 1, xyplanesvar, up(2), tag, &
                                lbecomm, status_down, ierr)
!
! 5) syncronize
              call mpi_wait(reqs_up(1), status_up, ierr)
              call mpi_wait(reqs_down(1), status_down, ierr)
            endif

          endif

          if(up(2) == myrank) then
           if(xperiodic)&
            call bc_edge_x_svar<<<(nx+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step,faiRho)
           if(yperiodic)&
            call bc_edge_y_svar<<<(ny+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step,faiRho)  
           if(zperiodic)&
            call bc_edge_z_svar<<<(nz+2*nbuff+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step,faiRho)
          else
            ! NON VA FATTO
              !call abortOnLastErrorAndSync(msg,step, flipflop)
              
!             mpi rho+pop myrank <--> up(2)
            if(zperiodic)&
             call bc_edge_z_svar<<<(nz+2*nbuff+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step,faiRho)
          endif

          

      
    end subroutine fixPeriodic_svar
    
    subroutine fixPeriodic_back_force(faiF, msg)
#ifdef SERIAL
      implicit none
      character(len=*), intent(in) :: msg
#else
      use mpi
      implicit none
      character(len=*), intent(in) :: msg

      integer      :: reqs_up(5)
      integer      :: reqs_down(5)

      
      integer      :: status_up(MPI_STATUS_SIZE)
      integer      :: status_down(MPI_STATUS_SIZE)
      integer      :: tag

#endif
      logical :: faiF

          if(.not. lbcforce)return
 
          if(xperiodic)&
           call bc_per_x_backf<<<dimGridx, dimBlock2>>>(step,faiF)
          if(yperiodic)&
           call bc_per_y_backf<<<dimGridy, dimBlock2>>>(step,faiF)
          
          
          call abortOnLastErrorAndSync(msg,step, flipflop)

! if it is the same task, there's no need to send/recv data
          if(up(2) == myrank) then
             if(zperiodic)&
              call bc_per_z_backf<<<dimGridz, dimBlock2>>>(step,faiF)
          else 
            
            
            if (faiF) then
! 
! I've to exchange pop with other task
! -----------------------------------------------------------------------------------
! send BLUE rho
! 1) isend to up task 
              tag = 605
              call mpi_isend(force_d(1,0,0,   0), 1, xyplane3, up(2), tag, &
                                lbecomm, reqs_up(1), ierr)
!
! 2) recv from down task
              tag = 605
              call  mpi_recv(force_buf_down_d(1,0,0,   0), 1, xyplane3, down(2), tag, &
                                lbecomm, status_up, ierr) !force_d(1,0,0,  nz)
!
! 3) isend to down task
              tag = 606
              call mpi_isend(force_d(1,0,0,nz+1), 1, xyplane3, down(2), tag, &
                                lbecomm, reqs_down(1), ierr)
!
! 4) recv from up task
              tag = 606
              call  mpi_recv(force_buf_up_d(1,0,0,1), 1, xyplane3, up(2), tag, &
                                lbecomm, status_down, ierr) !force_d(1,0,0,   1)
!
! 5) syncronize
              call mpi_wait(reqs_up(1), status_up, ierr)
              call mpi_wait(reqs_down(1), status_down, ierr)
              
            endif
            
            call merge_per_z_backf<<<dimGridz, dimBlock2>>>(step,faiF)

          endif

          if(up(2) == myrank) then
            if(xperiodic)&
             call bc_edge_x_backf<<<(nx+2+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step,faiF)
            if(yperiodic)&
             call bc_edge_y_backf<<<(ny+2+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step,faiF)  
            if(zperiodic)&
             call bc_edge_z_backf<<<(nz+2+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step,faiF)
          else
            ! NON VA FATTO
              !call abortOnLastErrorAndSync(msg,step, flipflop)
              
!             mpi rho+pop myrank <--> up(2)
             if(zperiodic)&
              call bc_edge_z_backf<<<(nz+2+TILE_DIM-1)/TILE_DIM, TILE_DIM>>>(step,faiF)
          endif


      
    end subroutine fixPeriodic_back_force



    subroutine setVirtualNodes
      implicit none
      if (xperiodic) then
          call bc_periodic_ext_x<<<dimGridx, dimBlock2>>>(step, flipflop)
          call bc_periodic_ext_x<<<dimGridx, dimBlock2>>>(step, flipflop)
      endif
      if (yperiodic) then
          call bc_periodic_ext_y<<<dimGridy, dimBlock2>>>(step, flipflop)
          call bc_periodic_ext_y<<<dimGridy, dimBlock2>>>(step, flipflop)
      endif
      if (zperiodic) then
          call bc_periodic_ext_z<<<dimGridz, dimBlock2>>>(step, flipflop)
          call bc_periodic_ext_z<<<dimGridz, dimBlock2>>>(step, flipflop)
      endif
    end subroutine setVirtualNodes


    subroutine SetConstants
      implicit none
      integer l, i,j,k
      real, dimension(3) :: v1,v2
#ifdef DIFFDENS      
      if(densR.ne.ONE)then
        densR_d=densR
      else
        densR_d=ONE
      endif
      if(densB.ne.ONE)then
        densB_d=densB
      else
        densB_d=ONE
      endif
      
      if (densR.ne.1.0 .and. densB.ne.1.0) then
        write (*,*) 'the density of at leat one of the two fluids should be one', &
         densR,densB
        call mystop
      endif
#else
      densR = ONE
      densB = ONE
      densR_d = densR
      densB_d = densB
#endif  
      denswallR_d = denswallR
      denswallB_d = denswallB    
      
      viscR = cssq*(tauR-HALF)
      viscB = cssq*(tauB-HALF)
      omega = ONE / tauR
      
      viscR_d = viscR
      viscB_d = viscB
      
      sigma_CG_d = sigma_CG
      beta_CG_d = beta_CG
      
      nbuff_d = nbuff
      
      f_cost_d = f_cost
      A_rep_d = A_rep

      numscp_d = numscp
      
      store_vel_d = store_vel
      
      totrho = 0.0
      totrho_d = totrho
      totrhobuff_d = totrho
      
      
      
      niterVTK_d = nIterVTKmin
      
      if(densR.eq.1.0)then
        alphaCG(1)=p(0)
        alphaCG(2)=(-ONE+alphaCG(1))*(densR/densB)+ONE
      else
        alphaCG(2)=p(0)
        alphaCG(1)=(ONE/(-ONE+alphaCG(2)))*(densB/densR)+ONE
      endif
      
      do l=1,2
        alphaCG_d(l)=alphaCG(l)
      enddo
      
      cz_d = real(glz+1) * HALF

      do l=0,npops-1
        p_d(l) = p(l)
        ex_d(l) = ex(l)
        ey_d(l) = ey(l)
        ez_d(l) = ez(l)
        opp_d(l) = opp(l)
        a_d(l) = a(l)
        phi_d(l) = phi(l)
        varphi_d(l) = varphi(l)
        psi_d(l) = psi(l)
        xi_d(l) = xi(l)
        b_l_d(l) = b_l(l)
        rec_fact_d(l) = sqrt( real(ex(l)*ex(l) + ey(l)*ey(l) + ez(l)*ez(l)) )
      enddo
      
      do l=0,npops-1
		v1(1) = ex(l)
		v1(2) = ey(l)
		v1(3) = ez(l)
	   	v2(1:3) = v1(1:3)
        cmat_d(1:3,1:3,l) = tensor_product(v1,v2)
      enddo
      
      do l=0,linksd3q27
        exd3q27_d(l) = exd3q27(l)
        eyd3q27_d(l) = eyd3q27(l)
        ezd3q27_d(l) = ezd3q27(l)
        ad3q27_d(l)=ad3q27(l)
      enddo
      
      
      fx_d = f_cost(1)
      fy_d = f_cost(2)
      fz_d = f_cost(3)

      ext_fxx_d = ext_fxx
      ext_fyy_d = ext_fyy
      ext_fzz_d = ext_fzz
      if (lrotate) then
        ext_tqx_d = ext_tqx
        ext_tqy_d = ext_tqy
        ext_tqz_d = ext_tqz
      endif

      buz_d = pref_bouzidi
      wall_x0_d = wall_x0
      wall_x1_d = wall_x1
      wall_y0_d = wall_y0
      wall_y1_d = wall_y1
      wall_z0_d = wall_z0
      wall_z1_d = wall_z1

      ! Max speed
      max_vel_d = max_vel

      ! Double belt stuff..
      ndouble = 98
      
      allocate(exdouble(ndouble),eydouble(ndouble),ezdouble(ndouble))
      l = 0
      do k = -2,2
        do j = -2,2
          do i = -2,2
            if(i==-2 .or. i==2 .or. j==-2 .or. j==2 .or. k==-2 .or. k==2)then
               l = l + 1
               exdouble(l) = i
               eydouble(l) = j
               ezdouble(l) = k
               ! write (*,*) l, "e[x,y,z]double=", i,j,k
            endif
          enddo
        enddo
      enddo
    end subroutine SetConstants
    
    pure function tensor_product(v,w)

     implicit none

     real, intent(in),dimension(3) :: v,w
     real,dimension(3,3) :: tensor_product
  
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
    
    subroutine SetupParts
      implicit none
      real(4) dx, dist_3d(3), dist
      integer i,j, iterOut


      if (nx < 2*rdim .or.ny < 2*rdim .or. nz < 2*rdim) then
        write (*,*) 'Particle does not fit the lattice:rdim=',rdim, 'Box=', nx,ny,nz
        call mystop
      endif

      if (myrank== 0) write(*,*) 'Initial setup of particles...'

      call read_input_atom

      call zeroListAtomsGPU<<<dimGridSubd, dimBlockSubd>>>(step, flipflop)
      call makeListAtomsGPU<<<dimGridAtm, dimBlockAtm>>>(step, flipflop, numAtoms)

      ! Sanity checks...
      iterOut = numAtoms / 10
      do i=1,numAtoms
        if (x_atm(1,i) < 1 .or. x_atm(1,i)>glx .or. x_atm(2,i)<1 .or. x_atm(2,i)>gly .or. x_atm(3,i)<1 .or. x_atm(3,i)>glz) then
          write (*,*) 'Particle #',i, 'out of bounds', x_atm(:,i)
          call mystop
        endif

        do j=i+1,numAtoms
          dist_3d = x_atm(:,i) - x_atm(:,j)
          dist = norm2(dist_3d)
          if (dist < 2*rdim+1) then
            write (*,fmt=200) 'Dist less than 2*rdim+1',i,j, dist, '<',rdim*2+1
            200     format(A, 2I8, F8.4, A, F8.4)
            write (*,*) 'Atom',i, x_atm(:,i)
            write (*,*) 'Atom',j, x_atm(:,j)
            call mystop
          else
            ! write (*,fmt=200) 'Dist more than 2*rdim+1',i,j, dist, '>=',rdim*2+1
          endif
        enddo
      enddo
    end subroutine SetupParts


    subroutine read_input_atom
      implicit none
      real                  :: x,y,z, vx,vy,vz, eul(3), qs(4)
      character(len = 120)  :: string, directive
      integer               :: i, itest, l

      open(unit=11, file='input.xyz', status='old', form='FORMATTED', action='read', iostat=itest)
      if (itest/=0) then
        write (*,*) "input.xyz] Error opening file. STOP"
        call mystop
      endif

      ! Number of atoms
      read(11, *, end=120) i
      if (myrank==0) write (*,*) 'input.xyz] Using:', numAtoms, ' particles from a file with:', i
      if (i < numAtoms) then        
        write (*,*) 'Input file:input.xyz contains less particles than needed. STOP'
        call mystop
      endif

      allocate(x_atm(3, numAtoms), STAT=istat)
      if (istat /= 0) then
        write(*,*) 'Allocation of particle positions failed'
      end if

      allocate(v_atm(3, numAtoms), STAT=istat)
      if (istat /= 0) then
        write(*,*) 'Allocation of particle veloc failed'
      end if

      if(lrotate) then
        allocate(q(4, numAtoms), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of particle q failed'
        endif
      endif

      ! Format line
      read(11, '(a120)', end=120) string

      do l = 1, numAtoms
        if (.not. lrotate) then
          read(11,*) string, x,y,z, vx,vy,vz
        else
          read(11,*) string, x,y,z, vx,vy,vz, eul(1),eul(2),eul(3)
        endif

        x_atm(1,l) = x
        x_atm(2,l) = y
        x_atm(3,l) = z
        v_atm(1,l) = vx
        v_atm(2,l) = vy
        v_atm(3,l) = vz
        if (lrotate) then
          call eul2q(eul(1),eul(3),eul(2),q(1:4, l))
          ! write(*,*) 'read_input_atom]q=  ',  q(1:4, l)
          ! write(*,*) 'read_input_atom]eul=',  eul(1),eul(2),eul(3)
        endif
      enddo

      120 continue
      close(11)

      do l = 1, numAtoms
        do i = 1, 3
          pos(i,l,1) = x_atm(i,l)   ! New & old pos
          pos(i,l,2) = x_atm(i,l)
          pos(3+i,l,1) = v_atm(i,l) ! New & old veloc
          pos(3+i,l,2) = v_atm(i,l)
        enddo            

        do i = 7, 9
          pos(i,l,1) = 0.0
        enddo

        forceAtoms(1,l,1) = ext_fxx
        forceAtoms(2,l,1) = ext_fyy
        forceAtoms(3,l,1) = ext_fzz
        forceAtoms(4,l,1) = 0.0
        forceAtoms(5,l,1) = 0.0
        forceAtoms(6,l,1) = 0.0

        forceAtoms(1,l,2) = ext_fxx
        forceAtoms(2,l,2) = ext_fyy
        forceAtoms(3,l,2) = ext_fzz
        forceAtoms(4,l,2) = 0.0
        forceAtoms(5,l,2) = 0.0
        forceAtoms(6,l,2) = 0.0

        if (lrotate) then
          forceAtoms(7,l,1) = ext_tqx
          forceAtoms(8,l,1) = ext_tqy
          forceAtoms(9,l,1) = ext_tqz
          forceAtoms(10,l,1) = 0.0
          forceAtoms(11,l,1) = 0.0
          forceAtoms(12,l,1) = 0.0

          forceAtoms(7,l,2) = ext_tqx
          forceAtoms(8,l,2) = ext_tqy
          forceAtoms(9,l,2) = ext_tqz
          forceAtoms(10,l,2) = 0.0
          forceAtoms(11,l,2) = 0.0
          forceAtoms(12,l,2) = 0.0
          
            do i = 1, 4
            pos(9+i,l,1) = q(i,l)
            pos(9+i,l,2) = q(i,l)
          enddo
        endif

      enddo

      pos_d = pos
      forceAtoms_d = forceAtoms

      if (myrank==0) then        
        write(6,'(a)') "File input.xyz correctly closed"
        write(6,*) "*******************************************************************************"
        write(6,*) "                                                                               "
      endif
    end subroutine read_input_atom

    subroutine eul2q(phis,psis,thetas,qs)
        implicit none
        
        real, intent(in) :: phis,psis,thetas
        real, intent(out), dimension(0:3) :: qs        
        real :: cy,sy,cp,sp,cr,sr
        
        cy = cos(phis * HALF)
        sy = sin(phis * HALF)
        cp = cos(psis * HALF)
        sp = sin(psis * HALF)
        cr = cos(thetas * HALF)
        sr = sin(thetas * HALF)
        
        qs(0) = cy * cp * cr + sy * sp * sr
        qs(1) = cy * cp * sr - sy * sp * cr
        qs(2) = sy * cp * sr + cy * sp * cr
        qs(3) = sy * cp * cr - cy * sp * sr
       end subroutine eul2q

    subroutine updateListAtoms(flip)
      implicit none
      integer, intent(in)   :: flip

      ! On th GPU
      call system_clock(c_startAtom)

      call zeroListAtomsGPU<<<dimGridSubd, dimBlockSubd>>>(1, flip)
      call makeListAtomsGPU<<<dimGridAtm, dimBlockAtm>>>(1, flip, numAtoms)

      call system_clock(c_stopAtom)
      diff_listatomGPU = diff_listatomGPU + c_stopAtom-c_startAtom
    end subroutine updateListAtoms

    subroutine printParams
      implicit none
      real volatms
      integer i

      write(*,*) '************************'
      write(*,fmt=103) 'Lattice box', glx,gly,glz
      write(*,fmt=100) 'nIters',nIter
      write(*,fmt=100) 'nIterOut',nIterOut
      write(*,fmt=100) 'initType',initType
      write(*,fmt=102) 'lreadisfluid',lreadisfluid
      write(*,fmt=102) 'lwriteisf',lwriteisf
      write(*,fmt=102) 'lwriteisf_every',lwriteisf_every
      write(*,fmt=108)lattice
      
      if (withCG) then
        write(*,fmt=100) 'Fluid components', 2
        write(*,fmt=105) 'Color gradient'
        write(*,fmt=101) 'Color gradient sigma', sigma_CG
        write(*,fmt=101) 'Color gradient beta', beta_CG
        
        write(*,fmt=109) 'Initial density ',densR,densB
        write(*,fmt=109) 'Set fluid wall density',denswallR,denswallB
        write(*,fmt=104) 'Initial velocity ',vx,vy,vz

#ifdef NEARCONTACT
        write(*,fmt=105) 'Near contact'
        write(*,fmt=101) 'Near contact A_rep', A_rep
#endif
        if (withParticles) then
          write(*,fmt=103) 'With Particles'
          write(*,fmt=104) 'External force', ext_fxx,ext_fyy,ext_fzz
          if (lrotate) then
            write(*,fmt=103) 'With Rotation'
            write(*,fmt=104) 'External Torque', ext_tqx,ext_tqy,ext_tqz
          else
            write(*,fmt=103) 'No Rotation'
          endif
        else
          write(*,fmt=103) 'No Particles'
        endif
      else
        write(*,fmt=100) 'Fluid components', 1
        write(*,fmt=105) 'BGK'
        write(*,fmt=104) 'Initial fluid vel', vx,vy,vz
      endif
      
      if (withSCP) then
        write(*,fmt=100) 'Passive scalar', numscp
      endif  
      
      write(*,fmt=102) 'X axis is periodic', xperiodic
      write(*,fmt=102) 'Y axis is periodic', yperiodic
      write(*,fmt=102) 'Z axis is periodic', zperiodic
      write(*,fmt=101) 'Fluid tau', tauR
      write(*,fmt=101) 'Fluid viscosity',  viscR
      write(*,fmt=101) 'Fluid blue tau', tauB
      write(*,fmt=101) 'Fluid blue viscosity',  viscB
      write(*,fmt=102) 'Constant omega mode', uniqueOmega
      write(*,fmt=101) 'Omega',  omega
      write(*,fmt=102) 'Forced', forced
      write(*,fmt=102) 'const_forced', const_forced
      if (forced .and. const_forced) write(*,fmt=104) 'Const Forced', (f_cost(i), i=1,3)
      write(*,fmt=100) 'numbc', numbc 
      if(numbc>0)then
        do i=1,numbc
          if(withCG)then
            write(*,'(a,i4,a,2i4)')'bc number: ',i,'   type: ',bctype(1:2,i)
            write(*,'(a,2f8.4)')'bc rho: ',bcrho(1:2,i)
          else
            write(*,'(a,i4,a,i4)')'bc number: ',i,'   type: ',bctype(1,i)
            write(*,'(a,1f8.4)')'bc rho: ',bcrho(1:1,i)
          endif
          write(*,'(a,3f8.4)')'bc vel: ',bcvel(1:3,i)
          if(withSCP)then
            write(*,'(a,9f8.4)')'bc scp: ',bcscp(1:numscp,i)
          endif
        enddo
      endif
      write(*,fmt=102) 'store_vel', store_vel
      write(*,fmt=102) 'VTK output', wantOut
      if (wantOut) then
        write(*,fmt=100) 'VTK every', nIterVTK 
        write(*,fmt=102) 'VTK ASCII', textVTK 
      endif
      write(*,*) '************************'
      if (withParticles) then
        write(*,fmt=100) 'Particles', numAtoms
        write(*,fmt=101) 'Radius', rdim
        write(*,fmt=101) 'Weight', atmWeight
        ! write(*,fmt=100) 'Size of findAtoms_d', dimGrid%x*dimGrid%y*dimGrid%z
        ! write(*,fmt=100) 'Size of listAtoms_d', numAtoms* dimGrid%x*dimGrid%y*dimGrid%z
        volatms = 4.0/3.0 * Pi * rdim**THREE
        write(*,fmt=106) 'Vol. fraction (%)', 100*volatms*numAtoms / real(glx*gly*glz)
      endif
      if(numVTK2d>0)then
        
        write(*,fmt=100) '2D VTK every', nIterVTK2d 
        write(*,fmt=100)'numVTK2d',numVTK2d 
        do i=1,numVTK2d
          select case(ivecVTK2d(i))
          case(1)
	        write(*,'(a)')'VTK2d plane: YZ'
	        write(*,'(a,8i)')'VTK2d at X : ',iselVTK2d(i)
	      case(2)
	        write(*,'(a)')'VTK2d plane: XZ'
	        write(*,'(a,8i)')'VTK2d at Y : ',iselVTK2d(i)
	      case(3)
	        write(*,'(a)')'VTK2d plane: XY'
	        write(*,'(a,8i)')'VTK2d at Z : ',iselVTK2d(i)
	      end select
        enddo
      endif
      write(*,*) '************************'

      write(*,*)

100    format(a20, ':', I8)
101    format(a20, ':', F20.8)
102    format(a20, ':', L)
103    format(a20, 3I8 )
104    format(a20, ':', 3F20.8)
105    format('              Method:', a20)
106    format(a20, ':', F20.2)
107    format(a20, ':', 5F20.8)
108    format('             Lattice:', a5)
109    format(a20, ':', 2F20.8)
    end subroutine printParams


    subroutine makeOutput
      implicit none
      character(len=120) :: makedirectory
      logical lexist

      makedirectory=repeat(' ',255)
      makedirectory = 'output'
      inquire(file=trim(makedirectory),exist=lexist)
        
      if(.not. lexist)then
        if(myrank==0) then
          makedirectory=repeat(' ',255)
          makedirectory = 'mkdir output'
          call system(makedirectory)
        endif
      endif
      
      call test_little_endian(lelittle)
      
    end subroutine makeOutput

    
    subroutine dumpTerms_n2p(step)
      implicit none
      integer, intent(in) :: step
      character(len=120) :: mynamefile
      integer :: l, int1,int2
      integer :: countn2p
      real :: debugn2pf(20000,5)
      real :: debugn2pt(20000,5)

      mynamefile=repeat(' ',120)
      mynamefile='output/debugN2pf_'//write_fmtnumb(myrank) // '.' // write_fmtnumb(step)//'.txt'
      open(unit=111, file=trim(mynamefile), status='replace')

      mynamefile='output/debugN2pt_'//write_fmtnumb(myrank) // '.' // write_fmtnumb(step)//'.txt'
      open(unit=112, file=trim(mynamefile), status='replace')

      countn2p = countn2p_d
      debugn2pf = debugn2pf_d
      debugn2pt = debugn2pt_d

      do l= 1, countn2p
        int1 = debugn2pf(l,1)
        int2 = debugn2pf(l,2)
        write(111,fmt=200) int1,int2,debugn2pf(l,3), debugn2pf(l,4), debugn2pf(l,5)

        int1 = debugn2pt(l,1)
        int2 = debugn2pt(l,2)
        write(112,fmt=200) int1,int2,debugn2pt(l,3), debugn2pt(l,4), debugn2pt(l,5)
      enddo

200    format(2I12,' ', 3F20.15)

      close(111)
      close(112)
    end subroutine dumpTerms_n2p

    subroutine dumpTerms_p2n(step)
      implicit none
      integer, intent(in) :: step
      character(len=120) :: mynamefile
      integer :: l, int1,int2,int3,int4
      integer :: countp2n1, countp2n2
      real :: debugp2n1(20000,4)
      real :: debugp2n2(20000,4)

      mynamefile=repeat(' ',120)
      mynamefile='output/debugP2n1_'//write_fmtnumb(myrank) // '.' // write_fmtnumb(step)//'.txt'
      open(unit=111, file=trim(mynamefile), status='replace')

      countp2n1 = countp2n1_d
      debugp2n1 = debugp2n1_d

      do l= 1, countp2n1
        int1 = debugp2n1(l,1)
        int2 = debugp2n1(l,2)
        int3 = debugp2n1(l,3)
        int4 = debugp2n1(l,4)
        write(111,fmt=200) int1,int2, int3,int4
200     format(4I12)
      enddo
      close(111)


      mynamefile='output/debugP2n2_'//write_fmtnumb(myrank) // '.' // write_fmtnumb(step)//'.txt'
      open(unit=111, file=trim(mynamefile), status='replace')

      countp2n2 = countp2n2_d      
      debugp2n2 = debugp2n2_d

      do l= 1, countp2n2
        int1 = debugp2n2(l,1)
        int2 = debugp2n2(l,2)
        write(111,fmt=201) int1,int2, debugp2n2(l,3), debugp2n2(l,4)
201     format(2I12,' ', 2F20.15)
      enddo
      close(111)
    end subroutine dumpTerms_p2n    

    subroutine dumpRot(step, flip)
      implicit none
      integer, intent(in) :: step, flip
      character(len=120) :: mynamefile
      integer :: iatm

      ! Only master dumps...
      if(myrank/=0) return

      mynamefile='output/debugRot_'//write_fmtnumb(step)//'.txt'
      open(unit=111, file=trim(mynamefile), status='replace')

      debugline = debugline_d
      
      100 format (A,3G20.7)
      101 format (A,4G20.7)
      102 format (A,F6.0,G20.7)
      103 format (A,I8, 4G20.7)
      200 format (A,4F20.8)

      iatm = 1

      write(111,103) 'rot1]tqx', step, debugline(iatm, 15,1:3)
      write(111,101) 'rot1]q',   debugline(iatm, 16,1:4)
      write(111,100) 'rot1]oxx', debugline(iatm, 17,1:3)
      write(111,100) 'rot1]rot13',  debugline(iatm, 18,1:3)
      write(111,100) 'rot1]rot46',  debugline(iatm, 19,1:3)
      write(111,100) 'rot1]rot79',  debugline(iatm, 20,1:3)
      write(111,100) 'rot1]rotinx', debugline(iatm, 21,1:3)

      write(111,100) 'rot1]tr',  debugline(iatm, 11,1:3)
      write(111,101) 'rot1]dq',  debugline(iatm, 12,1:4)
      write(111,101) 'rot2]q ',  debugline(iatm, 13,1:4)
      write(111,102) 'rot2]it',  debugline(iatm, 14,1:2)

      q = pos_d(10:13,:, flip)
      write(111,200) 'rot2]x ',take_rotversorx(q(1,1),q(2,1),q(3,1),q(4,1))
      close(111)
    end subroutine dumpRot

    subroutine dumpForces(step)
      implicit none
      integer, intent(in) :: step
      character(len=120) :: mynamefile
      integer :: iatm, int1,int2,int3


      ! Only master dumps...
      if(myrank/=0) return

      debugline = debugline_d
      v_atm = pos_d(4:6,:, 3-flipflop)

      100 format (A,I8, 3G20.7)
      101 format (A,I8, 3G20.7, I8, 3G20.7)
      102 format (A, 3I8)

      do iatm = 1, numAtoms
        mynamefile = 'output/' // 'debugForces_atm'// trim(write_fmtnumb(iatm)) &
          // '_step' // trim(write_fmtnumb(step)) // '.txt'
        open(unit=111,file=trim(mynamefile),status='replace',action='write')

        write(111,100) 'fxBx] step:', step, debugline(iatm, 1,1:3)
        write(111,100) 'fxx ] step:', step, debugline(iatm, 2,1:3)
        write(111,100) 'fxb ] step:', step, debugline(iatm, 3,1:3)
        write(111,100) 'fxbo] step:', step, debugline(iatm, 4,1:3)

        write(111,100) 'tqx ] step:', step, debugline(iatm, 5,1:3)
        write(111,100) 'txb ] step:', step, debugline(iatm, 6,1:3)
        write(111,100) 'txbo] step:', step, debugline(iatm, 7,1:3)

        int1 = debugline(iatm, 8,4)
        write(111,101) 'partPos:', step, debugline(iatm, 8,1:3), int1, v_atm(1, iatm), v_atm(2, iatm), v_atm(3, iatm)

        int2 = debugline(iatm, 9,1)
        int3 = debugline(iatm, 9,2)
        if (int2+int3>0) write(111,102) 'countmk, countrm=', step, int2,int3
        
        close(111)
      enddo
    end subroutine dumpForces

    subroutine mpiAddForces()
#ifdef SERIAL
      implicit none
#else
      use mpi
      implicit none
      integer      :: istat4(4)

      ! istat4(1) = cudaMemcpy(buffer_for(1,  1), myf_d , 3*numAtoms )
      ! istat4(2) = cudaMemcpy(buffer_for(4,  1), myf2_d, 3*numAtoms )
      ! if (lrotate) then
      !   istat4(3) = cudaMemcpy(buffer_for(7, 1), myt_d, 3*numAtoms )
      !   istat4(4) = cudaMemcpy(buffer_for(10,1), myt2_d, 3*numAtoms )
      ! endif
      buffer_for(1:3,  1:numAtoms) = myf_d
      buffer_for(4:6,  1:numAtoms) = myf2_d
      if (lrotate) then
        buffer_for(7:9,  1:numAtoms) = myt_d
        buffer_for(10:12,  1:numAtoms) = myt2_d
      endif

      istat = cudaDeviceSynchronize
      if (istat/=0) then
        write(*,*) 'status after mpiAddForces1 Sync:', cudaGetErrorString(istat)
        write(*,*) 'status after mpiAddForces1 Sync:', istat4
        write(*,*) 'Exiting ....'
        call mystop
      endif

      if (lrotate) then
        call mpi_allreduce(buffer_for,buffer_forall, 12*numAtoms,MPI_REAL8, &
          MPI_SUM, lbecomm, ierr)
      else
        call mpi_allreduce(buffer_for,buffer_forall, 6*numAtoms,MPI_REAL8, &
          MPI_SUM, lbecomm, ierr)
      endif
      
      ! istat4(1) = cudaMemcpy(myf_d,  buffer_forall(1,  1), 3*numAtoms )
      ! istat4(2) = cudaMemcpy(myf2_d, buffer_forall(4,  1), 3*numAtoms )
      ! if (lrotate) then
      !   istat4(3) = cudaMemcpy(myt_d , buffer_forall(7, 1), 3*numAtoms )
      !   istat4(4) = cudaMemcpy(myt2_d, buffer_forall(10,1), 3*numAtoms )
      ! endif

      myf_d = buffer_forall(1:3,  1:numAtoms)
      myf2_d = buffer_forall(4:6,  1:numAtoms)
      if (lrotate) then
        myt_d = buffer_forall(7:9,  1:numAtoms)
        myt2_d = buffer_forall(10:12,  1:numAtoms)
      endif

      istat = cudaDeviceSynchronize
      if (istat/=0) then
        write(*,*) 'status after mpiAddForces2 Sync:', cudaGetErrorString(istat)
        write(*,*) 'status after mpiAddForces2 Sync:', istat4
        write(*,*) 'Exiting ....'
        call mystop
      endif
#endif
    end subroutine mpiAddForces

    subroutine print_rho_pops2(filenam,step,flip)
      implicit none
      integer, intent(in) :: step, flip
      character(len=*), intent(in) :: filenam
      character(len=120) :: mynamefile
      integer :: i,j,k,l, iosub = 55

!       mynamefile=repeat(' ',120)
!       mynamefile='output/'//trim(filenam)//write_fmtnumb(step)//'.'//write_fmtnumb(myrank)//'.dat'
!       open(unit=iosub, file=trim(mynamefile), status='replace')
        
!       do k=0,nz+1
!         do j=0,ny+1
!           do i=0,nx+1
!               do l=0,npops-1
!     !!                              1,2,3,4,    5,		        6,         
!                 write(iosub,fmt=200)i,j,k,l,pop_pinned(i,j,k,l), popB_pinned(i,j,k,l)
! 200    format(4I5,' ', 2F20.15)
!               enddo
!           enddo
!         enddo
!       enddo    
!       close(iosub)
    
      mynamefile=repeat(' ',120)
      mynamefile='output/'//trim(filenam)//write_fmtnumb(step)//'.'//write_fmtnumb(myrank)//'.dat1'
      open(unit=iosub+1,file=trim(mynamefile),status='replace')
    
      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1            
    !!                               1,2,3,      4,          5,       !    6
              write(iosub+1,fmt=201) i,j,offset(3)+k, rhoR(i,j,k), rhoB(i,j,k)  !, myfluid(i,j,k, flip)
201    format(3I5,' ', 2F20.12, I2)
          enddo
        enddo
      enddo
    
      close(iosub+1)
   end subroutine print_rho_pops2

   subroutine dumpPartvol(step)
    implicit none
    integer, intent(in) :: step
    character(len=120) :: mynamefile
    integer :: l


    mynamefile=repeat(' ',120)
    mynamefile='output/partVol_'//write_fmtnumb(myrank) // '.' // write_fmtnumb(step)//'.txt'
    open(unit=111, file=trim(mynamefile), status='replace')

    do l= 1, numAtoms
      write(111,fmt=200) l, partVol(l)
      200    format(I8,' ', I8)
    enddo

    close(111)
  end subroutine dumpPartvol

  subroutine dumpMkRm(step)
    implicit none
    integer, intent(in) :: step
    character(len=120) :: mynamefile
    integer :: l, int1,int2, int5
    integer :: countmk, countrm
    real :: debug_rm(20000,5)
    real :: debug_mk(20000,5)


    countmk = countmk_d
    debug_mk = debug_mk_d

    if (countmk>0) then
      mynamefile=repeat(' ',120)
      mynamefile='output/debugMK-f_'//write_fmtnumb(myrank) // '.' // write_fmtnumb(step)//'.txt'
      open(unit=111, file=trim(mynamefile), status='replace')

      do l = 1, countmk
        int1 = debug_mk(l,1)
        int5 = debug_mk(l,5)
        write(111,fmt=200) int1, debug_mk(l,2), debug_mk(l,3), debug_mk(l,4), int5
200    format(I12,' ', 3F20.15, I12)
      enddo

      close(111)
    endif

    countrm = countrm_d
    debug_rm = debug_rm_d

    if (countrm>0) then
      mynamefile='output/debugRM-f_'//write_fmtnumb(myrank) // '.' // write_fmtnumb(step)//'.txt'
      open(unit=111, file=trim(mynamefile), status='replace')

      do l = 1, countrm
        int1 = debug_rm(l,1)
        int5 = debug_rm(l,5)
        write(111,fmt=200) int1, debug_rm(l,2), debug_rm(l,3), debug_rm(l,4), int5
      enddo

      close(111)
    endif
  end subroutine dumpMkRm

   subroutine OutputVTK(fname, step, flip)
      implicit none
      integer, intent(in) :: step, flip
      character(len=*), intent(in) :: fname
      integer(8) :: c_start,c_stop,diff_xfer,c_start2, pops_time
      character(len=120) :: mynamefile, mynamefile0, mynamefile1
      real(4)     :: plane0(nx+2,ny+2), plane1(nx+2,ny+2)
      logical, save :: lfirst=.true.
      integer, parameter :: myunit=865
      integer :: isel,i,mydest
      character(len=13) :: string13a
      character(len=14) :: string14a
      
      if(lfirst)then
        lfirst=.false.
        if(myrank== 0)then
          open(unit=myunit,file='bcvalue.dat',status='replace',action='write')
          write(myunit,'(a)')'# totrho'
        endif
      endif
      
      myfluid = myfluid_d
      istat = cudaDeviceSynchronize
      if (istat/=0) write(*,*) 'status after myfluid',istat
      
      ! if (fname(1:5) /= 'step_' .and. fname(1:4)/= 'stop') return
      
      call system_clock(c_start)
      if(lpoptransf)then
        istat = cudaMemcpy(pop_pinned, popsR_d(:,:,:,:,flip), (nx+2)*(ny+2)*(nz+2)*npops )
        if (istat/=0) write(*,*) 'status after copy R:',istat, 'flipflop',flipflop
        istat = cudaDeviceSynchronize
        if (istat/=0) write(*,*) 'status after deviceSync',istat
      endif
      if (.not. withCG) then
        if(lpoptransf)then
          call moments1Fl(nz, pop_pinned,rhoR, vel, myfluid, flip)
        else
          istat = cudaMemcpy(rhoR, rhoR_d, (nx+2*nbuff)*(ny+2*nbuff)*(nz+2*nbuff) )
          if (istat/=0) write(*,*) 'status after rho R:',istat
          istat = cudaDeviceSynchronize
          if (istat/=0) write(*,*) 'status after rho R',istat
          if(lprintvel)then
            istat = cudaMemcpy(vel, vel_d, (nx+2*nbuff)*(ny+2*nbuff)*(nz+2*nbuff)*3 )
            if (istat/=0) write(*,*) 'status after vel:',istat
            istat = cudaDeviceSynchronize
            if (istat/=0) write(*,*) 'status after vel',istat
          endif
        endif
      endif
      call system_clock(c_stop)
      diff_xfer = c_stop-c_start
      
      if (withCG) then
        call system_clock(c_start2)
        if(lpoptransf)then
          istat = cudaMemcpy(popB_pinned, popsB_d(:,:,:,:,flip), (nx+2)*(ny+2)*(nz+2)*npops )
          if (istat/=0) write(*,*) 'status after copy B:',istat, 'flipflop',flipflop
          istat = cudaDeviceSynchronize
          if (istat/=0) write(*,*) 'status after deviceSync B',istat
         
          call moments2Fl(nz, pop_pinned,rhoR, popB_pinned,rhoB, vel, phase, myfluid, flip)
        else
          istat = cudaMemcpy(rhoR, rhoR_d, (nx+2*nbuff)*(ny+2*nbuff)*(nz+2*nbuff) )
          if (istat/=0) write(*,*) 'status after rho R:',istat
          istat = cudaDeviceSynchronize
          if (istat/=0) write(*,*) 'status after rho R',istat

          istat = cudaMemcpy(rhoB, rhoB_d, (nx+2*nbuff)*(ny+2*nbuff)*(nz+2*nbuff) )
          if (istat/=0) write(*,*) 'status after rho B:',istat
          istat = cudaDeviceSynchronize
          if (istat/=0) write(*,*) 'status after rho B',istat
          if(lprintvel)then
            istat = cudaMemcpy(vel, vel_d, (nx+2*nbuff)*(ny+2*nbuff)*(nz+2*nbuff)*3 )
            if (istat/=0) write(*,*) 'status after vel:',istat
            istat = cudaDeviceSynchronize
            if (istat/=0) write(*,*) 'status after vel',istat
          endif
        endif
         call system_clock(c_stop)
         diff_xfer = diff_xfer + c_stop-c_start2

         call system_clock(c_start2)
         if (debug1) call print_rho_pops2(fname,step, flip)
         call system_clock(c_stop)
         pops_time = c_stop-c_start2

         mynamefile = repeat(' ',120)
         if (fname(1:5) == 'step_') then
          mynamefile = 'rhoB'
         else
          mynamefile = trim(fname)//'rhoB'
         endif
         
         if(lprintrhoB)then
           do k=1,nz
             do j=1,ny
               do i=1,nx
                 if(myfluid(i,j,k,flip) /= fluid_fluid )then
                   rhoB(i,j,k) = ZERO
                 endif
               enddo
             enddo
           enddo
         endif
         if(lprintphase)then
           do k=1,nz
             do j=1,ny
               do i=1,nx
                 if(myfluid(i,j,k,flip) == fluid_fluid )then
                   phase(i,j,k) = (rhoR(i,j,k)-rhoB(i,j,k))/(rhoR(i,j,k)+rhoB(i,j,k))
                 else
                   phase(i,j,k) = ZERO
                 endif
               enddo
             enddo
           enddo
         endif
         
      endif
      
      do k=1,nz
         do j=1,ny
           do i=1,nx
             if(myfluid(i,j,k,flip) /= fluid_fluid )then
               rhoR(i,j,k) = ZERO
             endif
           enddo
         enddo
       enddo
      
      if(lprintvel)then
        do k=1,nz
          do j=1,ny
            do i=1,nx
              if(myfluid(i,j,k,flip) /= fluid_fluid )then
                vel(1:3,i,j,k) = ZERO
              endif
            enddo
          enddo
        enddo
      endif
      
      
      
      if(lwriteisf .and. (lwriteisf_every .or. step==0))then
        mynamefile = repeat(' ',120)
        if (fname(1:5) == 'step_') then
           mynamefile = 'isfluid'
        else
           mynamefile = trim(fname)//'_isfluid'
        endif
        call writeImageDataVTI_isfluid_nopart(nz, mynamefile,step, myfluid, flip, textVTK)
        if(numVTK2d>0)then
          do i=1,numVTK2d
            isel=iselVTK2d(i)
            string13a='isf2d_'//write_fmtnumb(i)//'_'
            select case(ivecVTK2d(i))
            case(1)
              call writeImageDataVTI_YZ(isel,nz,string13a,step, 'isfluid', rhoR, &
               textVTK, myfluid, flip,.true.)
            case(2)
              call writeImageDataVTI_XZ(isel,nz,string13a,step, 'isfluid', rhoR, &
               textVTK, myfluid, flip,.true.)
            case(3)
              mydest=GET_RANK_POINT(nx/2,ny/2,isel)
              call writeImageDataVTI_XY(isel,mydest,nz,string13a,step,'isfluid', rhoR, &
               textVTK, myfluid, flip,.true.)
            end select
          enddo
        endif
      endif


      
      if (myrank== 0) then

          
          if(lcompute_totrho)then
            istat = cudaMemcpy(totrho,totrho_d,1)
            if (istat/=0) write(*,*) 'status after copy totrho:',istat, &
             'flipflop',flipflop
            istat = cudaDeviceSynchronize
            if (istat/=0) write(*,*) 'status after deviceSync',istat
          endif
          
          write(myunit,'(f20.10)')totrho
        
      endif
      
      if(numVTK2d>0)then
        if(mod(step,niterVTK2d).eq.0)then
          do i=1,numVTK2d
            isel=iselVTK2d(i)
            string13a='rho2d_'//write_fmtnumb(i)//'_'
            select case(ivecVTK2d(i))
            case(1)
              call writeImageDataVTI_YZ(isel,nz,string13a,step, 'densR', rhoR, &
               textVTK, myfluid, flip,.false.)
               
              if(lprintrhoB)then
               string14a='rhoB2d_'//write_fmtnumb(i)//'_'
               call writeImageDataVTI_YZ(isel,nz,string14a,step, 'densB', &
               rhoB, textVTK, myfluid, flip,.false.)
              endif
               
              if(lprintphase)then
               string14a='phas2d_'//write_fmtnumb(i)//'_'
               call writeImageDataVTI_YZ(isel,nz,string14a,step, 'phase', &
               phase, textVTK, myfluid, flip,.false.)
              endif
              
              if(lprintvel)then
                string13a='vel2d_'//write_fmtnumb(i)//'_'
                call writeImageDataVTI_YZ_3d(isel,nz,string13a,step, 'vel', &
                 vel, textVTK, myfluid, flip)
              endif
               
            case(2)
              call writeImageDataVTI_XZ(isel,nz,string13a,step, 'densR', rhoR, &
               textVTK, myfluid, flip,.false.)
               
              if(lprintrhoB)then
               string14a='rhoB2d_'//write_fmtnumb(i)//'_'
               call writeImageDataVTI_XZ(isel,nz,string14a,step, 'densB', &
               rhoB, textVTK, myfluid, flip,.false.)
              endif
               
              if(lprintphase)then
               string14a='phas2d_'//write_fmtnumb(i)//'_'
               call writeImageDataVTI_XZ(isel,nz,string14a,step, 'phase', &
               phase, textVTK, myfluid, flip,.false.)
              endif
              
              if(lprintvel)then
                string13a='vel2d_'//write_fmtnumb(i)//'_'
                call writeImageDataVTI_XZ_3d(isel,nz,string13a,step, 'vel', &
                 vel, textVTK, myfluid, flip)
              endif
               
            case(3)
              mydest=GET_RANK_POINT(nx/2,ny/2,isel)
              call writeImageDataVTI_XY(isel,mydest,nz,string13a,step, 'densR', &
               rhoR, textVTK, myfluid, flip,.false.)
               
              if(lprintrhoB)then
               string14a='rhoB2d_'//write_fmtnumb(i)//'_'
               call writeImageDataVTI_XY(isel,mydest,nz,string14a,step, 'densB', &
               rhoB, textVTK, myfluid, flip,.false.)
              endif
               
              if(lprintphase)then
               string14a='phas2d_'//write_fmtnumb(i)//'_'
               call writeImageDataVTI_XY(isel,mydest,nz,string14a,step, 'phase', &
               phase, textVTK, myfluid, flip,.false.)
              endif
               
              if(lprintvel)then
                string13a='vel2d_'//write_fmtnumb(i)//'_'
                call writeImageDataVTI_XY_3d(isel,mydest,nz,string13a,step, 'vel', &
                 vel, textVTK, myfluid, flip)
              endif
               
            end select
          enddo
          call system_clock(c_stop)
          if (myrank== 0) then
            write(*,fmt=201) (c_stop-c_start)/real(clock_rate), diff_xfer/real(clock_rate), pops_time/real(clock_rate)
201         format('Time for 2DVTK:', G10.2, ' secs (tranfer time:', G10.2, ' secs, pops:',G10.2, ')')
          endif
        endif
      endif
      
      
      if(mod(step,niterVTK).ne.0)return
      
      call writeImageDataVTI(nz,'rho',step, 'densR', rhoR, textVTK, myfluid, flip)
      if (withCG) then
        if(lprintrhoB)then
          call writeImageDataVTI(nz,'rhoB',step, 'rhoB', rhoB, textVTK, myfluid, flip)
        endif
        if(lprintphase)then
          call writeImageDataVTI(nz, "phase",step, 'phase', phase, textVTK, myfluid, flip)
        endif
      endif
      if(lprintvel)call writeImageDataVTI_3d(nz, "vel",step,"velocity",vel, textVTK, myfluid, flip)
      

      if (withParticles) then
        call system_clock(c_start2)

        x_atm = pos_d(1:3,:, flip)
        v_atm = pos_d(4:6,:, flip)
        q = pos_d(10:13,:, flip)

        call system_clock(c_stop)
        diff_xfer = diff_xfer + c_stop-c_start2

        mynamefile = repeat(' ',120)
        if (fname(1:5) == 'step_') then
          mynamefile = 'particle'
        else
          mynamefile = trim(fname)//'particle'
        endif
        call writeParticleVTK(mynamefile,step, x_atm,v_atm,q, textVTK)

        ! mynamefile = repeat(' ',120)
        ! if (fname(1:5) == 'step_') then
        !   mynamefile = 'isfluid'
        ! else
        !   mynamefile = trim(fname)//'_isfluid'
        ! endif
        ! if (fname(1:5) == 'step_' .and. step/=0) then
        !   ! After time step, use old value
        !   call writeImageDataVTI_isfluid(nz, mynamefile,step, myfluid, 3-flip, textVTK, totSphere)
        ! else
        !   call writeImageDataVTI_isfluid(nz, mynamefile,step, myfluid, flip, textVTK, totSphere)
        ! endif
     endif

      call system_clock(c_stop)
      if (myrank== 0) then
        write(*,fmt=301) (c_stop-c_start)/real(clock_rate), diff_xfer/real(clock_rate), pops_time/real(clock_rate)
301     format('Time for 3DVTK:', G10.2, ' secs (tranfer time:', G10.2, ' secs, pops:',G10.2, ')')
      endif
      
   end subroutine OutputVTK


   subroutine DumpBorderPlanez(fname, step, flip)
    implicit none
    integer, intent(in) :: step, flip
    character(len=*), intent(in) :: fname
    integer(8) :: c_start,c_stop,diff_xfer,c_start2, pops_time
    character(len=120) :: mynamefile0, mynamefile1
    integer     :: l
    real(4)     :: plane0(nx+2,ny+2), plane1(nx+2,ny+2)
    

    call system_clock(c_start)

    do l = 0, -1
      istat = cudaMemcpy(plane0, popsR_d(:,:,0,l,flip), (nx+2)*(ny+2) )
      istat = cudaMemcpy(plane1, popsR_d(:,:,nz+1,l,flip), (nx+2)*(ny+2) )
      if (istat/=0) write(*,*) 'status after copy R:',istat, 'flipflop',flipflop

      istat = cudaDeviceSynchronize
      if (istat/=0) write(*,*) 'status after deviceSync',istat

      if (fname(1:5) == 'step_') then
        mynamefile0 = 'plane0_popR' // trim(write_fmtnumb0(l,1))
        mynamefile1 = 'plane1_popR' // trim(write_fmtnumb0(l,1))
      else
        mynamefile0 = trim(fname)//'plane0_popR' // trim(write_fmtnumb0(l,1))
        mynamefile1 = trim(fname)//'plane1_popR' // trim(write_fmtnumb0(l,1))
      endif

      call writeImageDataVTI_2d(mynamefile0,step, 'pop'// trim(write_fmtnumb0(l,1)), plane0, 0, textVTK)
      call writeImageDataVTI_2d(mynamefile1,step, 'pop'// trim(write_fmtnumb0(l,1)), plane1, nz+1, textVTK)
    enddo

    call system_clock(c_stop)
    diff_xfer = c_stop-c_start

    istat = cudaMemcpy(plane0, rhoR_d(:,:,0), (nx+2)*(ny+2) )
    if (istat/=0) write(*,*) 'status after copy R0:',istat
    istat = cudaDeviceSynchronize
    if (istat/=0) write(*,*) 'status after R0 deviceSync',istat

    istat = cudaMemcpy(plane1, rhoR_d(:,:,nz+1), (nx+2)*(ny+2) )
    if (istat/=0) write(*,*) 'status after copy R1:',istat
    istat = cudaDeviceSynchronize
    if (istat/=0) write(*,*) 'status after R1 deviceSync',istat

    if (fname(1:5) == 'step_') then
      mynamefile0 = 'plane0_rhoR'
      mynamefile1 = 'plane1_rhoR'
    else
      mynamefile0 = trim(fname)//'plane0_rhoR'
      mynamefile1 = trim(fname)//'plane1_rhoR'
    endif

    call writeImageDataVTI_2d(mynamefile0,step, 'rho1', plane0, 0, textVTK)
    call writeImageDataVTI_2d(mynamefile1,step, 'rho1', plane1, nz+1, textVTK)
  end subroutine


   subroutine spherical_template
    implicit none
    
    integer :: i,j,k,l,m,rmax,rmin,ishift,jshift,kshift
    integer :: ioppshift,joppshift,koppshift, istat
    real :: rdist,sqrcut
    integer(1), allocatable, dimension(:,:,:) :: issub
    character(len=120) :: extent
        
    sqrcut = rdim**TWO
    
    rmin = floor(rdim)
    rmax = rdim + 1    
    if (myrank==0) write(*,fmt=10) rmin,rmax, sqrcut,rmax*rmax
    10 format('spherical_template] rmin,rmax=', 2I5, ' sqrcut in',F12.4,I5)


    allocate(issub(-rmax:rmax,-rmax:rmax,-rmax:rmax))
    issub(:,:,:) = fluid_fluid
        
    do k = -rmax,rmax
      do j = -rmax,rmax
        do i = -rmax,rmax
          rdist = i**TWO + j**TWO + k**TWO
          if (rdist<=sqrcut) then
            issub(i,j,k) = fluid_spheredead
            if (i>rmin .or. i<-rmin .or. j>rmin .or. j<-rmin .or. k>rmin .or. k<-rmin) then
              write (*,*) 'Strange error here. Stopping'
              call mystop
            endif
          endif
        enddo
      enddo
    enddo
    
    do k = -rmin,rmin
      do j = -rmin,rmin
        do i = -rmin,rmin
          if (issub(i,j,k)==fluid_spheredead) then

            do l = 1, npops-1
              if (issub(i+ex(l),j+ey(l),k+ez(l)) == fluid_fluid) then
                issub(i,j,k) = fluid_spherelist
                exit
              endif
            enddo

          endif
        enddo
      enddo
    enddo
    
    l = 0
    do k = -rmax,rmax
      do j = -rmax,rmax
        do i = -rmax,rmax
          if (issub(i,j,k)==fluid_spherelist) then
            l = l+1
          endif
        enddo
      enddo
    enddo
    
    issub(0,0,0) = fluid_particleCM
    
    !!!!!!!!!!! Now calc spherelist, spherelistdead
    nsphere = l
    allocate(spherelist(3,nsphere))

    l = 0
    do k = -rmax,rmax
      do j = -rmax,rmax
        do i = -rmax,rmax
          if (issub(i,j,k)==fluid_spherelist) then
            l = l+1
            spherelist(1,l) = i
            spherelist(2,l) = j
            spherelist(3,l) = k
          endif
        enddo
      enddo
    enddo
  
  
    l = 0
    do k = -rmax,rmax
      do j = -rmax,rmax
        do i = -rmax,rmax
          if (issub(i,j,k)==fluid_spheredead) then
            l = l+1
          endif
        enddo
      enddo
    enddo
    
    nspheredead = l
    allocate(spherelistdead(3,nspheredead))

    totSphere = nsphere + nspheredead + 1
    if (myrank==0) write (*,fmt=11) nsphere,nspheredead, totSphere
    11 format('nsphere,nspheredead=', 2I5,' totSphere=', I5)

    
    l = 0
    do k = -rmax,rmax
      do j = -rmax,rmax
        do i = -rmax,rmax
          if (issub(i,j,k)==fluid_spheredead) then
            l = l+1
            spherelistdead(1,l) = i
            spherelistdead(2,l) = j
            spherelistdead(3,l) = k
          endif
        enddo
      enddo
    enddo

    allocate(issub_d(-rmax:rmax,-rmax:rmax,-rmax:rmax))
    rmax_issub_d = rmax
    sphereMax_d = rmax*rmax
    ! issub_d(-rdim:rdim, -rdim:rdim, -rdim:rdim) = issub(-rdim:rdim, -rdim:rdim, -rdim:rdim)
    issub_d = issub
    
    ! do k = -rmax,rmax
    !   do j = -rmax,rmax
    !     do i = -rmax,rmax
    !       write(21,*) i,j,k, issub(i,j,k)
    !     enddo
    !   enddo
    ! enddo

    ! call debug_issub<<<1,1>>>

    istat = cudaGetLastError()
    if (istat/=0) then
      write(*,*) 'status after spherical_template:',cudaGetErrorString(istat)
      write(*,*) 'Exiting ....'
      call mystop
    endif

   end subroutine spherical_template


   subroutine setupGPUgridAndAlloc
    implicit none
    integer*8 gpuMEM, numElems

    if (mod(nx, TILE_DIMx)/= 0) then
        write(*,*) 'nx must be a multiple of TILE_DIM'
        call mystop
    end if
    if (mod(ny, TILE_DIMy) /= 0) then
        write(*,*) 'ny must be a multiple of TILE_DIMy'
        call mystop
    end if
    if (mod(nz, TILE_DIMz) /= 0) then
        write(*,*) 'nz must be a multiple of TILE_DIMz'
        call mystop
    end if
    
  
    ! if (withParticles .and. mod(numAtoms, TILE_DIMPART) /= 0) then
    !    write(*,*) 'numAtoms(',numAtoms,') must be a multiple of TILE_DIMPART(',TILE_DIMPART,')'
    !    call mystop
    ! end if
  
    dimGrid  = dim3(nx/TILE_DIMx, ny/TILE_DIMy, nz/TILE_DIMz)
    dimBlock = dim3(TILE_DIMx, TILE_DIMy, TILE_DIMz)
    !(ceiling((ny+2)/TILE_DIM)) 
    dimGridx  = dim3((ny+2+TILE_DIM-1)/TILE_DIM, (nz+2+TILE_DIM-1)/TILE_DIM, 1)
    dimGridy  = dim3((nx+2+TILE_DIM-1)/TILE_DIM, (nz+2+TILE_DIM-1)/TILE_DIM, 1)
    dimGridz  = dim3((nx+2+TILE_DIM-1)/TILE_DIM, (ny+2+TILE_DIM-1)/TILE_DIM, 1)
    
    dimBlock2 = dim3(TILE_DIM, TILE_DIM, 1)
    
    
    dimGridhalox = dim3((ny+2*nbuff+TILE_DIM-1)/TILE_DIM, (nz+2*nbuff+TILE_DIM-1)/TILE_DIM, 1)
    dimGridhaloy = dim3((nx+2*nbuff+TILE_DIM-1)/TILE_DIM, (nz+2*nbuff+TILE_DIM-1)/TILE_DIM, 1)
    dimGridhaloz = dim3((nx+2*nbuff+TILE_DIM-1)/TILE_DIM, (ny+2*nbuff+TILE_DIM-1)/TILE_DIM, 1)
    
    dimGridsidex = dim3((ny+TILE_DIM-1)/TILE_DIM, (nz+TILE_DIM-1)/TILE_DIM, 1)
    dimGridsidey = dim3((nx+TILE_DIM-1)/TILE_DIM, (nz+TILE_DIM-1)/TILE_DIM, 1)
    dimGridsidez = dim3((nx+TILE_DIM-1)/TILE_DIM, (ny+TILE_DIM-1)/TILE_DIM, 1)
    
    dimGridhalo = dim3((nx+2*nbuff+TILE_DIMx-1)/TILE_DIMx, &
     (ny+2*nbuff+TILE_DIMy-1)/TILE_DIMy,(nz+2*nbuff+TILE_DIMz-1)/TILE_DIMz)
    dimBlockhalo = dim3(TILE_DIMx, TILE_DIMy, TILE_DIMz)
  
    dimGridAtm  = (numAtoms+TILE_DIMPART-1)/TILE_DIMPART
    dimBlockAtm = TILE_DIMPART

    dimGridSubd  = dim3((dimGrid%x+TILE_DIMx-1) / TILE_DIMx, (dimGrid%y+TILE_DIMy-1) / TILE_DIMy, (dimGrid%z+TILE_DIMz-1) / TILE_DIMz)
    dimBlockSubd = dim3(TILE_DIMx, TILE_DIMy, TILE_DIMz)
  
    ! For timing
    istat = cudaEventCreate(startEvent)
    istat = cudaEventCreate(stopEvent)
    call system_clock(count_rate=clock_rate)
  
    ! write parameters
  
    if (myrank== 0) then
	    write(*,*)
	    write(*,fmt=10) nx,ny,nz,npops
	    write(*,fmt=11) dimGrid%x, dimGrid%y, dimGrid%z, dimBlock%x, dimBlock%y, dimBlock%z
	    write(*,fmt=12) dimGridx%x, dimGridx%y, dimGridx%z, dimBlock2%x, dimBlock2%y, dimBlock2%z
	    write(*,fmt=13) dimGridy%x, dimGridy%y, dimGridy%z
	    write(*,fmt=14) dimGridz%x, dimGridz%y, dimGridz%z
	    write(*,fmt=15) dimGridAtm, dimBlockAtm
	    write(*,fmt=20) dimGridSubd%x,dimGridSubd%y,dimGridSubd%z, dimBlockSubd%x,dimBlockSubd%y,dimBlockSubd%z
    endif
  
    10 format('Matrix size:', 4i5)
    11 format('dimGrid:', 3i4, ',   dimBlock:', 3i4)
    12 format('dimGridx:', 3i4, ',   dimBlock2:', 3i4)
    13 format('dimGridy:', 3i4)
    14 format('dimGridz:', 3i4)
    15 format('dimGridAtm:', i4, ',   dimBlockAtm:', i4)
    20 format('dimGridSubd:', 3i4, ',   dimBlockSubd:', 3i4)

    if (myrank== 0) write(*,*) 'Mem size:', dataSz/1024./1024.0, 'MB'

    if (myrank== 0) write(*,*) 'Allocation of HOST pops array, pinned:', pinnedFlag
    allocate(pop_pinned(0:nx+1,0:ny+1,0:nz+1,0:npops-1), STAT=istat, PINNED=pinnedFlag)
    if (istat /= 0) then
       write(*,*) 'Allocation of pop_pinned failed'
       call mystop
    end if
    pop_pinned = 0.0
  
    allocate(rhoR(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff), STAT=istat)
    if (istat /= 0) then
       write(*,*) 'Allocation of rhoR failed'
       call mystop
    end if
    allocate(vel(3,1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff), STAT=istat)
    if (istat /= 0) then
       write(*,*) 'Allocation of velocity failed'
       call mystop
    end if

    allocate(myfluid(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff, 2), STAT=istat)
    if (istat /= 0) then
      write(*,*) 'Allocation of particle myfluid failed'
      call mystop
    end if

    if (myrank== 0) write(*,*) 'Allocation of GPU pops:'
    allocate(popsR_d(0:nx+1,0:ny+1,0:nz+1,0:npops-1, 2), STAT=istat)
    if (istat /= 0) then
       write(*,*) 'Allocation of popsR_d failed'
       call mystop
    end if
  
    allocate(rhoR_d(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff), STAT=istat)
    if (istat /= 0) then
       write(*,*) 'Allocation of rhoR_d failed'
       call mystop
    end if
    
    if(withSCP)then
      allocate(popSCP_pinned(1:numscp,0:nx+1,0:ny+1,0:nz+1,0:npops_scp-1), STAT=istat, PINNED=pinnedFlag)
      if (istat /= 0) then
         write(*,*) 'Allocation of popSCP_pinned failed'
         call mystop
      end if
      popSCP_pinned = 0.0
    
      allocate(popsSCP_d(1:numscp,0:nx+1,0:ny+1,0:nz+1,0:npops_scp-1, 2), STAT=istat)
      if (istat /= 0) then
         write(*,*) 'Allocation of popsSCP_d failed'
         call mystop
      end if
      
      allocate(scalar(1:numscp,1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff), STAT=istat)
      if (istat /= 0) then
         write(*,*) 'Allocation of scalar failed'
         call mystop
      end if
      
      allocate(scalar_d(1:numscp,1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff), STAT=istat)
      if (istat /= 0) then
         write(*,*) 'Allocation of scalar_d failed'
         call mystop
      end if
    
    endif
    
    if(store_vel) then
      allocate(vel_d(1:3,1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff), STAT=istat)
      if (istat /= 0) then
        write(*,*) 'Allocation of vel_d failed'
      endif
      vel_d = 0.0
    endif

    if (.not. withCG .and. forced .and. .not. const_forced) then
      allocate(force_d(3,1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff), STAT=istat)
      if (istat /= 0) then
        write(*,*) 'Allocation of force_d failed'
      endif
      if(lbcforce)then
        allocate(force_buf_up_d(3, 0:nx+1,0:ny+1,1), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of force_buf_up_d failed'
        endif
        
        allocate(force_buf_down_d(3, 0:nx+1,0:ny+1,1), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of force_buf_down_d failed'
        endif
      endif
    endif
    
    numbc_d = numbc
    
    if(numbc>0)then
      
      allocate(bcvel_d(3,numbc), STAT=istat)
      if (istat /= 0) then
        write(*,*) 'Allocation of bcvel_d failed'
      endif
      bcvel_d = bcvel
      
      if(withCG)then
        allocate(bctype_d(2,numbc), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of bctype_d failed'
        end if
        bctype_d = bctype
      
        allocate(bcrho_d(2,numbc), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of bcrho_d failed'
        endif
      else
        allocate(bctype_d(1,numbc), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of bctype_d failed'
        end if
        bctype_d = bctype
        
        allocate(bcrho_d(1,numbc), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of bcrho_d failed'
        endif
      endif
      bcrho_d = bcrho
      
      if(withSCP)then
        allocate(bcscptype_d(numbc), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of particle bcscptype_d failed'
        end if
        bcscptype_d = bcscptype
        
        allocate(bcscp_d(numscp,numbc), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of particle bcscp_d failed'
        endif
        bcscp_d = bcscp
      endif
    endif

    allocate(myfluid_d(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff,2), STAT=istat)
    if (istat /= 0) then
      write(*,*) 'Allocation of myfluid_d failed'
    end if
    allocate(debugfluid_d(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff), STAT=istat)
    if (istat /= 0) then
      write(*,*) 'Allocation of debugfluid_d failed'
    end if
    
  
    if (withCG) then
      if (myrank== 0) write(*,*) 'Allocation of HOST other pops:'
      allocate(popB_pinned(0:nx+1,0:ny+1,0:nz+1,0:npops-1), STAT=istat, PINNED=pinnedFlag)
      if (istat /= 0) then
        write(*,*) 'Allocation of popB_pinned failed'
        call mystop
      end if
      popB_pinned = 0.0
  
      allocate(rhoB(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff), STAT=istat)
      if (istat /= 0) then
        write(*,*) 'Allocation of rhoB failed'
      end if
  
      allocate(phase(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff), STAT=istat)
      if (istat /= 0) then
        write(*,*) 'Allocation of phase failed'
      end if

      if (myrank== 0) write(*,*) 'Allocation of GPU other pops:'
      allocate(popsB_d(0:nx+1,0:ny+1,0:nz+1,0:npops-1, 2), STAT=istat)
      if (istat /= 0) then
        write(*,*) 'Allocation of popsB_d failed'
        call mystop
      end if
      popsB_d = 0.0
  
      allocate(rhoB_d(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff), STAT=istat)
      if (istat /= 0) then
        write(*,*) 'Allocation of rhoB_d failed'
      end if
      
      allocate(nearsel_d(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff), STAT=istat)
      if (istat /= 0) then
        write(*,*) 'Allocation of nearsel_d failed'
      endif
      nearsel_d = 0

      if (withParticles) then
        if (myrank== 0) write(*,fmt=16) numAtoms
        16 format('Allocating ', i4, ' particles on GPU...')

        allocate(pos(13, numAtoms,2), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of pos failed'
        end if

        allocate(forceAtoms(12, numAtoms,2), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of forceAtoms failed'
        end if

        allocate(pos_d(13, numAtoms,2), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of pos_d failed'
        end if

        allocate(forceAtoms_d(12, numAtoms,2), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of forceAtoms_d failed'
        end if

        allocate(myf_d(3, numAtoms), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of myf_d failed'
        end if

        allocate(myf2_d(3, numAtoms), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of myf2_d failed'
        end if

        allocate(myt_d(3, numAtoms), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of myt_d failed'
        end if

        allocate(myt2_d(3, numAtoms), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of myt2_d failed'
        end if

        allocate(buffer_for(12, numAtoms), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of buffer_for failed'
        end if
        allocate(buffer_forall(12, numAtoms), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of buffer_forall failed'
        end if

        allocate(matrixAtoms_d(maxAtomInBlock, nx/TILE_DIMx, ny/TILE_DIMy, nz/TILE_DIMz), STAT=istat)
        numElems = maxAtomInBlock * nx/TILE_DIMx * ny/TILE_DIMy * nz/TILE_DIMz
        if (myrank== 0) write(*,*) 'Mem (MB) for matrixAtoms_d=', numElems * 4 / 1024**2
        if (istat /= 0) then
          write(*,*) 'Allocation of matrixAtoms_d failed'
        end if

        allocate(dim_matAtm_d(nx/TILE_DIMx, ny/TILE_DIMy, nz/TILE_DIMz), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of dim_matAtm_d failed'
        end if

        ! DEBUG_FORCE
        allocate(debugline(numAtoms, 21,4), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of CPU debugline failed'
        end if

        allocate(debugline_d(numAtoms, 21,4), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of debugline_d failed'
        end if

        ! DEBUG CHECK_VOLP
        allocate(partVol(numAtoms), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of partVol failed'
        end if

        allocate(partVol_d(numAtoms), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of partVol_d failed'
        end if
        allocate(oldpartVol_d(numAtoms), STAT=istat)
        if (istat /= 0) then
          write(*,*) 'Allocation of oldpartVol_d failed'
        end if
        
      endif
    endif
 
  
    ! output device info and transfer size
    istat = cudaGetDeviceProperties(prop, 0)
    if (istat/=0) then
      write(*,*) 'status after cudaGetDeviceProperties:',cudaGetErrorString(istat)
      write(*,*) 'Exiting ....'
      call mystop
    endif
  
    if (myrank== 0) then
      gpuMEM = (nx+2)*(ny+2)*(nz+2)
      gpuMEM = gpuMEM *(npops+1)*2* 2* 4/1024./1024.
      write(*,*)
      write(*,*) 'Device: ', trim(prop%name)

      write(*,*) '************************'
      write(*,*) 'GPU LB (R+B) size (MB): ', gpuMEM
      write(*,*) 'Transfer size (MB): ', gpuMEM / 4 * npops/(npops+1)
      istat = cudaEventRecord(startEvent, 0)

      istat = cudaMemcpy(pop_pinned, popsR_d, (nx+2)*(ny+2)*(nz+2)*npops )
      if (istat/=0) write(*,*) 'status after pops-gpuSync:',istat
      istat = cudaEventRecord(stopEvent, 0)
      istat = cudaEventSynchronize(stopEvent)

      istat = cudaEventElapsedTime(mytime, startEvent, stopEvent)
      write(*,*) '  Device to Host bandwidth (GiB/s): ', dataSz/mytime*(1.e+3/1024**3)
      write(*,*)
    endif
   end subroutine setupGPUgridAndAlloc


   subroutine abortOnLastErrorAndSync(msg, step, flipflop)
    implicit none
    integer, intent(in) :: step, flipflop
    character(len=*), intent(in) :: msg
    integer :: gpustop, gpustop_forall, ierr,istat0
    
    istat0 = cudaGetLastError()
    
    istat = cudaDeviceSynchronize
    if (istat0/=0) then
      write(*,*) 'status after ',msg,':', cudaGetErrorString(istat)
      write(*,*) 'Exiting ....'
      call mystop
    endif
    if (istat/=0) then
      write(*,*) 'status after ',msg,' Sync:', cudaGetErrorString(istat)
      write(*,*) 'Exiting ....'
      call mystop
    endif

    gpustop = stop_d
    call mpi_allreduce(gpustop, gpustop_forall, 1,MPI_INTEGER, &
          MPI_SUM, lbecomm, ierr)

    if (gpustop_forall/=0) then
      if (gpustop/=0) then
        write(*,*) 'status after ',msg,' mystop:', gpustop
        write(*,*) 'Exiting ....'
      endif
      call OutputVTK("stop-",step, flipflop)
      call mystop
    endif
   end subroutine

   subroutine read_and_setRed_Blue
    implicit none
    integer :: minx,maxx, miny,maxy, minz,maxz

    open(unit=11, file='rhoR.init.bin', form='unformatted',action='read',status='old')
    read(11) minx,maxx, miny,maxy, minz,maxz
    if (minx/=1 .or. miny/=1 .or. minz/=1 .or. maxx/=nx .or. maxy/=ny .or. maxz/=nz) then
      write(*,fmt=10) 'Wrong start rhoR file', minx,maxx, miny,maxy, minz,maxz
10  format(a, 6i8)          
      call mystop
    endif
    read(11) rhoR(1:nx, 1:ny, 1:nz)
    close(11)
    
    open(unit=11, file='rhoB.init.bin', form='unformatted',action='read',status='old')
    read(11) minx,maxx, miny,maxy, minz,maxz
    if (minx/=1 .or. miny/=1 .or. minz/=1 .or. maxx/=nx .or. maxy/=ny .or. maxz/=nz) then
      write(*,fmt=10) 'Wrong start rhoB file', minx,maxx, miny,maxy, minz,maxz
      call mystop
    endif        
    read(11) rhoB(1:nx, 1:ny, 1:nz)
    close(11)
    write(*,*) ''

    istat = cudaMemcpy(rhoR_d, rhoR(:,:,:), (nx+2)*(ny+2)*(nz+2) )
    if (istat/=0) write(*,*) 'status after copy rho R:',istat
    istat = cudaDeviceSynchronize
    if (istat/=0) write(*,*) 'status after sync rho R',istat

    istat = cudaMemcpy(rhoB_d, rhoB(:,:,:), (nx+2)*(ny+2)*(nz+2) )
    if (istat/=0) write(*,*) 'status after copy rho B:',istat
    istat = cudaDeviceSynchronize
    if (istat/=0) write(*,*) 'status after sync rho B',istat

    myfluid_d = fluid_fluid
    debugfluid_d = 0

    countmk_d = 0
    countrm_d = 0
    countn2p_d = 0
    countp2n1_d = 0
    countp2n2_d = 0

   end subroutine read_and_setRed_Blue

   subroutine read_and_setPops
    implicit none
    integer :: minx,maxx, miny,maxy, minz,maxz, l
    character(len=120) :: mynamefile

    do l = 0, npops-1
      mynamefile='popR.' // write_fmtnumb(l) //'.init.bin'
      open(unit=11, file=mynamefile, form='unformatted',action='read',status='old')
      read(11) minx,maxx, miny,maxy, minz,maxz
      if (minx/=1 .or. miny/=1 .or. minz/=1 .or. maxx/=nx .or. maxy/=ny .or. maxz/=nz) then
        write(*,fmt=10) 'Wrong start popsR file', minx,maxx, miny,maxy, minz,maxz
  10  format(a, 6i8)          
        call mystop
      endif
      read(11) pop_pinned(1:nx, 1:ny, 1:nz, l)
      close(11)
      
      mynamefile='popB.' // write_fmtnumb(l) //'.init.bin'
      open(unit=11, file=mynamefile, form='unformatted',action='read',status='old')
      read(11) minx,maxx, miny,maxy, minz,maxz
      if (minx/=1 .or. miny/=1 .or. minz/=1 .or. maxx/=nx .or. maxy/=ny .or. maxz/=nz) then
        write(*,fmt=10) 'Wrong start popsB file', minx,maxx, miny,maxy, minz,maxz
        call mystop
      endif        
      read(11) popB_pinned(1:nx, 1:ny, 1:nz, l)
      close(11)
    enddo

    istat = cudaMemcpy(popsR_d(:,:,:,:,flipflop), pop_pinned, (nx+2)*(ny+2)*(nz+2)*npops )
    if (istat/=0) write(*,*) 'status after copy pops R:',istat

    istat = cudaDeviceSynchronize
    if (istat/=0) write(*,*) 'status after sync pops  R',istat

    istat = cudaMemcpy(popsB_d(:,:,:,:,flipflop), popB_pinned, (nx+2)*(ny+2)*(nz+2)*npops )
    if (istat/=0) write(*,*) 'status after copy pops  B:',istat
    istat = cudaDeviceSynchronize
    if (istat/=0) write(*,*) 'status after sync pops  B',istat

    write(*,*) ''
   end subroutine read_and_setPops

   subroutine setupRho_Pops()
    implicit none
    integer :: minx,maxx, miny,maxy, minz,maxz
    integer(kind=1), allocatable, dimension(:,:,:) :: myis
    
    logical :: lexist
    integer :: i,ii,j,jj,k,kk,itemp,mytype
    
    if (myrank== 0) write(*,*) 'Initial setup...'
    rhoR = 0
    rhoB = 0

    if (withCG) then
      if (myrank== 0) write(*,*) 'Initial setup red & blue pops...'
      if (initType == 0) then
        call setup1<<<dimGrid, dimBlock>>>(vx,vy,vz)
      else if (initType == 1) then
        call setup2<<<dimGrid, dimBlock>>>(vx,vy,vz)

      else if (initType == 2) then
          call setup2_zplanes<<<dimGrid, dimBlock>>>(vx,vy,vz)

      else if (initType == 3) then
        if (nprocs /= 1) then
          write(*,*) 'initType == 2 only in serial...'
          call mystop
        endif

        call read_and_setRed_Blue()
        call setupPops<<<dimGrid, dimBlock>>>(vx,vy,vz)

      else if (initType == 4) then
        if (nprocs /= 1) then
          write(*,*) 'initType == 3 only in serial...'
          call mystop
        endif

        call read_and_setRed_Blue
        call read_and_setPops

      else if (initType == 5) then
        call setupCyl<<<dimGrid, dimBlock>>>(vx,vy,vz)

      else if (initType == 6) then 
        call set_random_dens_fluids

        istat = cudaMemcpy(rhoR_d, rhoR(:,:,:), (nx+2)*(ny+2)*(nz+2) )
        if (istat/=0) write(*,*) 'status after copy rho R:',istat
        istat = cudaDeviceSynchronize
        if (istat/=0) write(*,*) 'status after sync rho R',istat

        istat = cudaMemcpy(rhoB_d, rhoB(:,:,:), (nx+2)*(ny+2)*(nz+2) )
        if (istat/=0) write(*,*) 'status after copy rho B:',istat
        istat = cudaDeviceSynchronize
        if (istat/=0) write(*,*) 'status after sync rho B',istat

        call setupPops<<<dimGrid, dimBlock>>>(vx,vy,vz)
      else if (initType == 7) then 
        call setup2_xplanes<<<dimGrid, dimBlock>>>(vx,vy,vz)
      else if (initType == 8) then 
        call setup2_collid<<<dimGrid, dimBlock>>>(vx,vy,vz)
      else if (initType == 9) then
        call setup2_sphere<<<dimGrid, dimBlock>>>(vx,vy,vz,64.0,64.0,16.0,10.0)
      else if (initType == 10) then
        call setup2_cylinder<<<dimGrid, dimBlock>>>(vx,vy,vz,&
         152.0,16.0,25.0,20.0,2)
      else
        write(*,*) 'Wront init type=', initType
        call mystop
      endif
#ifdef PARTICLES
      if (withParticles) then
        partVol_d(1:numAtoms) = 0

        call init_isfluid_CG<<<dimGrid, dimBlock>>>(step, flipflop)
        call fixPeriodic_isfluid("init_isfluid_CG_step0")

#ifdef CHECK_VOLP
        partVol(1:numAtoms) = partVol_d(1:numAtoms)
        call dumpPartvol(step)
#endif
        partVol_d(1:numAtoms) = 0
      endif
#endif
    else
      if (initType == 0) then
        call setup<<<dimGrid, dimBlock>>>(vx,vy,vz)
      elseif (initType == 1) then
        call setup_p<<<dimGrid, dimBlock>>>(vx,vy,vz)
      elseif (initType == 2) then
        call setup_taylorgreen<<<dimGrid, dimBlock>>>(vx,vy,vz)
      endif  
    end if
    
  
  

    if(withSCP)then
      if (myrank== 0) write(*,*) 'Initial setup SCP...'
      call setup_SCP<<<dimGrid, dimBlock>>>(vx,vy,vz)    
    endif
      
    if(lreadisfluid)then
      if(myrank==0)write(6,'(a)')'Reading isfluid.init.bin'
      allocate(myis(glx,gly,glz))
      inquire(file='isfluid.init.bin',exist=lexist)
      if(.not. lexist)then
        write(*,'(a)') 'File isfluid.init.bin not flound!'
        call mystop
      endif
      open(unit=11, file='isfluid.init.bin', form='unformatted',action='read',status='old')
      read(11) minx,maxx, miny,maxy, minz,maxz
      if (minx/=1 .or. miny/=1 .or. minz/=1 .or. maxx/=glx .or. maxy/=gly .or. maxz/=glz) then
        write(*,fmt=10) 'Wrong start isfluid file', minx,maxx, miny,maxy, minz,maxz
10      format(a, 6i8)          
        call mystop
      endif
      read(11) myis(1:glx, 1:gly, 1:glz)
      close(11)
      
      do k=1,glz
        do j=1,gly
          do i=1,glx
            if(myis(i,j,k)<0)then
              itemp=-int(myis(i,j,k))
              if(itemp>numbc)then
                if(myrank==0)write(*,*)'BC type larger than numbc : ', &
                 itemp,numbc
                call mystop
              endif
            endif
          enddo
        enddo
      enddo
      
      myfluid(0:nx+1,0:ny+1,0:nz+1,1:2)=int(3,kind=1)
      minx=start_idx(1)
      miny=start_idx(2)
      minz=start_idx(3)
      maxx=minx+lsizes(1)-1
      maxy=miny+lsizes(2)-1
      maxz=minz+lsizes(3)-1
      do k=1,glz
        if(k<minz .or. k>maxz)cycle
        kk=k-minz+1
        do j=1,gly
          if(j<miny .or. j>maxy)cycle
          jj=j-miny+1
          do i=1,glx
            if(i<minx .or. i>maxx)cycle
            ii=i-minx+1
            myfluid(ii,jj,kk,1:2)=myis(i,j,k) !i+j+k!nint(rand_noseeded(i,j,k,1)*120.0) 
          enddo
        enddo
      enddo
      
      myfluid_d = myfluid
      
    else
    
      myfluid_d = fluid_fluid
      
    endif
    
   
    
    
    call abortOnLastErrorAndSync("setupRho_Pops",0, 1)
    
    if(lreadisfluid)deallocate(myis)

 
  end subroutine

  subroutine set_random_dens_fluids
    implicit none
    integer :: i,j,k, k1
    logical :: linit_seed = .false.
    real, parameter :: meanR = 0.5, stdevR = 1.0e-4
    real, parameter :: meanB = 0.5, stdevB = 1.0e-4
    real :: localRhor, myRandom
    INTEGER :: randomSZ
    integer, allocatable :: seed(:)
    
    if(linit_seed)then
      do k=1, nz
        do j=1, ny
          do i=1, nx
            rhoR(i,j,k) = meanR + stdevR*gauss()
          enddo
        enddo
      enddo
    else
      call RANDOM_SEED (size=randomSZ)
      write (*, *) 'Random seed size=', randomSZ
      allocate(seed(randomSZ))
      call random_seed(get=seed)
      write (*, *) seed
      
      do k1=1, glz
        do j=1, ny
          do i=1, nx
            ! localRhor  = meanR + stdevR*gauss_noseeded(i,j,k,1)

            call RANDOM_NUMBER (myRandom)
            localRhor  = meanR + stdevR*myRandom

            ! localRhor = 1000 * k1 + j + i*0.001

            if (offset(3)+1<=k1 .and. k1<=offset(3)+nz) then
              k = k1 - offset(3)
              rhoR(i,j,k) = localRhor
              if (rhoR(i,j,k)<minRho .or. rhoR(i,j,k)>maxRho) then
                write(*,fmt=20) 'rhoR', i,j,k, rhoR(i,j,k)
  20  format('set_random_dens_fluids]Fixing range error ', a, 3i8, f10.8)
                rhoR(i,j,k) = meanR   ! call mystop
              endif
              rhoB(i,j,k) = 1.0 - rhoR(i,j,k)
            else
              if (j==1 .and. i==1) write (*,*) 'skipping k=', k1
            endif

          enddo
        enddo
      enddo  
    endif

    ! if(linit_seed)then
    !   do k=1, nz
    !     do j=1, ny
    !       do i=1, nx
    !         rhoB(i,j,k) = meanB + stdevB*gauss()
    !       enddo
    !     enddo
    !   enddo
    ! else  
    !   do k=1, nz
    !     do j=1, ny
    !       do i=1, nx
    !         rhoB(i,j,k) = meanB + stdevB*gauss_noseeded(i,j,k,100)
    !         if (rhoB(i,j,k)<minRho .or. rhoB(i,j,k)>10) then
    !           write(*,fmt=20) 'rhoB', i,j,k, rhoB(i,j,k)
    !           rhoB(i,j,k) = meanB   ! call mystop
    !         endif
    !       enddo
    !     enddo
    !   enddo  
    ! endif
   end subroutine set_random_dens_fluids
  
  
   function gauss()
    implicit none
    real :: gauss
    real :: dtemp1,dtemp2
    real, parameter :: mylimit=1.e-32
  
    
    call random_number(dtemp1)
    call random_number(dtemp2)
    
    dtemp1=dtemp1*(ONE-mylimit)+mylimit
  
  ! Box-Muller transformation
    gauss=sqrt(- TWO *log(dtemp1))*cos(TWO*pi*dtemp2)
   end function gauss
  
   
   pure function gauss_noseeded(i,j,k,l)
    implicit none
    integer, intent(in) :: i,j,k,l
    integer :: kk
    real :: gauss_noseeded
    real :: dtemp1,dtemp2
    real, parameter :: mylimit=1.e-32
    real, parameter :: FIFTY = 50
    
    dtemp1=rand_noseeded(i,j,k,l)
    kk=nint(dtemp1*FIFTY)
    dtemp2=rand_noseeded(i,j,k,l+kk)
    
    dtemp1=dtemp1*(ONE-mylimit)+mylimit
    
  ! Box-Muller transformation
    gauss_noseeded=sqrt(- TWO *log(dtemp1))*cos(TWO*pi*dtemp2)
   end function gauss_noseeded
  
  
  pure function rand_noseeded(i,j,k,l)
    implicit none
    integer, intent(in) :: i,j,k,l
    integer :: isub,jsub,ksub,lsub,msub
    integer ::ii,jj
    real(4) :: s,t,u33,u97,csub,uni
    real :: rand_noseeded
    
    real(4), parameter :: c =  362436.0/16777216.0
    real(4), parameter :: cd= 7654321.0/16777216.0
    real(4), parameter :: cm=16777213.0/16777216.0
  
    
  ! initial values of i,j,k must be in range 1 to 178 (not all 1)
  ! initial value of l must be in range 0 to 168.      
    isub=mod(i,178)
    isub=isub+1
    
    jsub=mod(j,178)
    jsub=jsub+1
    
    ksub=mod(k,178)
    ksub=ksub+1
    
    lsub=mod(l,169)
    
  ! initialization on fly  
    ii=97
    s=0.0
    t=0.5
    do jj=1,24
      msub=mod(mod(isub*jsub,179)*ksub,179)
      isub=jsub
      jsub=ksub
      ksub=msub
      lsub=mod(53*lsub+1,169)
      if(mod(lsub*msub,64).ge.32)s=s+t
      t=0.5*t
    enddo
    u97=s
    
    ii=33
    s=0.0
    t=0.5
    do jj=1,24
      msub=mod(mod(isub*jsub,179)*ksub,179)
      isub=jsub
      jsub=ksub
      ksub=msub
      lsub=mod(53*lsub+1,169)
      if(mod(lsub*msub,64).ge.32)s=s+t
      t=0.5*t
    enddo
    u33=s
    uni=u97-u33
    if (uni.lt.0.0) uni = uni + 1.0
    csub = c-cd
    if (csub.lt.0.0) csub = csub + cm
    uni = uni-csub
    if (uni.lt.0.0) uni = uni+1.0
    rand_noseeded = real(uni)
   end function rand_noseeded
end program