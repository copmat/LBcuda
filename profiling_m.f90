#include "defines.h"
#define MPI
#define PRC 4
#define IOOUT 6
#define ZERO 0.0
 module profiling_m
 
!***********************************************************************
!     
!     LBsoft module containing subroutines for profiling code 
!     originally written in MUPHY by S. Melchionna et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     modified by: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 use dimensions_m, only : myrank,nprocs
 use cudafor
 
 implicit none
#ifdef MPI    
    include 'mpif.h'
#endif  
 
 integer, save :: ier
 
 private
 
 
 integer, parameter :: FUNCTION_NAME_LEN=40
 
 integer, save :: nsection = 0
 
 integer, parameter :: mxsection = 20
 integer, parameter :: mxroutine = 50
  
 type single_routine_timer_t
   real(kind=PRC) :: timing !< cumulated timing
   integer :: icall !< number of calls
   character(len=FUNCTION_NAME_LEN) :: aroutine=' ' !< registered routine names
   character(len=FUNCTION_NAME_LEN) :: asection=' ' !< parent section

 end type
  
 type routine_timer_t
   real(kind=PRC) :: timing(mxroutine) !< cumulated timing
   integer :: icall(mxroutine) !< number of calls
   character(len=FUNCTION_NAME_LEN) :: secname=' ' !< parent section
 end type
  
 type section_t
   character(len=FUNCTION_NAME_LEN) :: name_sec=' ' !< registered section name
   real(kind=PRC) :: timing             !< cumulated section timing
   integer :: icall                   !< number of section calls
   character(len=FUNCTION_NAME_LEN) :: aroutine(mxroutine)=' ' !< registered routine names
   integer :: nroutine = 0                                    !< number of routines in list
   type(routine_timer_t) :: timer_routine           !< routine timer
 end type

 type(section_t), save :: section_cum(mxsection)
 type(section_t), save :: section_blk(mxsection)
  
 type(single_routine_timer_t), save :: timer_standalone(mxroutine)
  
 integer :: nroutine_global=0
 integer :: id_routine_global(mxsection*1000)
 
 type (cudaEvent) :: timing_routine_cuda_i(mxroutine)
 type (cudaEvent) :: timing_routine_cuda_e(mxroutine)
 logical, dimension(mxroutine), save :: ltiming_routine_cuda_i=.false.
 logical, dimension(mxroutine), save :: ltiming_routine_cuda_e=.false.
 
 real(kind=PRC) :: timing_routine(mxroutine)
 character(len=FUNCTION_NAME_LEN) :: routine_name_global(mxsection*1000)=''
  
 logical :: lactive_routine=.false.
 integer :: active_level = 2
  
 type (cudaEvent), save :: tstart, tend, tstart0
  
 integer :: iotiming
    
 integer, save ::            n_io_unit = 0  ! counter for extra units
 integer, save ::            io_unit = 14   ! starting no. for extra units
  
 integer, protected, public, save :: itime_start = 1                    !< start of time counter
 integer, protected, public, save :: itime_counter = 0                          !< time counter
 integer, protected, public, save :: idiagnostic = 50 !< make diagnostic every
 logical, protected, public, save :: ldiagnostic = .false.
 
 integer, save :: clock_max_param=0 
 integer, save :: clock_rate_param=0
 
 real(kind=PRC), save :: clock_huge=real(0.d0,kind=PRC)
 
INTERFACE
    FUNCTION get_mem ( ) bind ( C, name = "get_mem" )
      USE ISO_C_BINDING, ONLY : c_double
      REAL(kind=c_double) :: get_mem
    END FUNCTION get_mem
END INTERFACE
 
 
 public :: timer_init,startPreprocessingTime,reset_timing_partial
 public :: print_timing_partial,printSimulationTime
 public :: print_timing_final,start_timing2,end_timing2
 public :: set_value_idiagnostic,set_value_ldiagnostic,startSimulationTime
 public :: print_memory_registration,print_memory_registration_cuda
 public :: get_memory_cuda,get_memory
 
 contains
 
 subroutine set_value_idiagnostic(itemp1)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the nstep interval for the 
!     diagnostic profiling
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: itemp1
  
  idiagnostic = itemp1
  
  return
  
 end subroutine set_value_idiagnostic
 
 subroutine set_value_ldiagnostic(ltemp1)
 
!***********************************************************************
!     
!     LBsoft subroutine for setting the nstep interval for the 
!     diagnostic profiling
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(in) :: ltemp1
  
  ldiagnostic = ltemp1
  
  return
  
 end subroutine set_value_ldiagnostic
 
 
  FUNCTION new_io_unit()
  
!***********************************************************************
!     
!     LBsoft subroutine for driving the timing services
!     originally written in MUPHY by S. Melchionna et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     modified by: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************

   INTEGER :: new_io_unit

   n_io_unit = n_io_unit + 1
   io_unit = io_unit + 1
   new_io_unit = io_unit

  END FUNCTION new_io_unit


  SUBROUTINE startPreprocessingTime()

!***********************************************************************
!     
!     LBsoft subroutine for driving the timing services
!     originally written in MUPHY by S. Melchionna et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     modified by: M. Lauricella
!     last modification June 2018
!     
!***********************************************************************
     integer :: istat
     !tstart0 = current_time()
     istat = cudaEventCreate(tstart0)
     istat = cudaEventRecord(tstart0, 0)
     
  ENDSUBROUTINE

  SUBROUTINE startSimulationTime()

!***********************************************************************
!     
!     LBsoft subroutine for driving the timing services
!     originally written in MUPHY by S. Melchionna et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     modified by: M. Lauricella
!     last modification June 2018
!     
!***********************************************************************
     integer :: istat
     !tstart = current_time()
     istat = cudaEventCreate(tstart)
     istat = cudaEventRecord(tstart, 0)
  ENDSUBROUTINE
  
  SUBROUTINE timer_init()
  
!***********************************************************************
!     
!     LBsoft subroutine for driving the timing services
!     originally written in MUPHY by S. Melchionna et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     modified by: M. Lauricella
!     last modification June 2018
!     
!***********************************************************************
  
  INTEGER :: isec


  iotiming = new_io_unit()
  IF(myrank==0) OPEN(iotiming,FILE='time.dat')
  IF(myrank==0) WRITE(iotiming,*) '# Number of MPI Processes:',nprocs

  timing_routine(:) = 0

  DO isec=1, SIZE(section_cum)

     section_cum(isec)%timer_routine%timing(:) = 0
     section_cum(isec)%timer_routine%icall(:) = 0

     section_blk(isec)%timer_routine%timing(:) = 0
     section_blk(isec)%timer_routine%icall(:) = 0

  ENDDO

  
 

  END SUBROUTINE timer_init
  
  SUBROUTINE printSimulationTime()
  
!***********************************************************************
!     
!     LBsoft subroutine for driving the timing services
!     originally written in MUPHY by S. Melchionna et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     modified by: M. Lauricella
!     last modification June 2018
!     
!***********************************************************************
      integer :: istat
      
      REAL(kind=PRC) :: ttt0
      REAL(kind=PRC), DIMENSION(1) :: ttt1,ttt2
  
      !ttt = current_time()
      istat = cudaEventCreate(tend)
      istat = cudaEventRecord(tend, 0)
      
      istat = cudaEventSynchronize(tend)
      istat = cudaEventElapsedTime(ttt0, tstart0, tend)
      ttt1=ttt0/1000.0
      

     IF(ttt1(1)< ZERO ) THEN
#ifdef WALLCLOCK
       if(abs(ttt1(1))>(clock_huge*0.5d0))then
         ttt1 = ttt - tstart0 + clock_huge
         if(ttt1(1)<0)ttt1 = ZERO
       else
         ttt1 = ZERO
       endif
#else
       ttt1 = ZERO
#endif
     ENDIF
     
     istat = cudaEventElapsedTime(ttt0, tstart, tend)
     ttt2=ttt0/1000.0
     IF(ttt2(1)< ZERO ) THEN
#ifdef WALLCLOCK
       if(abs(ttt2(1))>(clock_huge*0.5d0))then
         ttt2 = ttt - tstart + clock_huge
         if(ttt2(1)<0)ttt2 = ZERO
       else
         ttt2 = ZERO
       endif
#else
       ttt2 = ZERO
#endif
     ENDIF
 

     
     IF(myrank==0) THEN
        WRITE(IOOUT,'(/a/)') '             --------- END OF TIME CYCLE --------- '
        WRITE(IOOUT,'(a,f16.6,1x,a)') 'TOTAL TIME (PREPARE+SIM)',ttt1(1),'Seconds'
        WRITE(IOOUT,'(a,f16.6,1x,a)') 'TOTAL TIME (SIM)        ',ttt2(1),'Seconds'
     ENDIF
  ENDSUBROUTINE
  
  SUBROUTINE reset_timing_partial()
  
!***********************************************************************
!     
!     LBsoft subroutine for driving the timing services
!     originally written in MUPHY by S. Melchionna et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     modified by: M. Lauricella
!     last modification June 2018
!     
!***********************************************************************

  INTEGER :: isec

    ! reset counter and cumulators
    nsection = 0 ! number of sections called

    DO isec=1, SIZE(section_cum)

       section_cum(isec)%timing = 0
       section_cum(isec)%nroutine = 0 ! number of routines called
       section_cum(isec)%icall = 0

       section_cum(isec)%timer_routine%icall(:) = 0 ! number of routine calls
       section_cum(isec)%timer_routine%timing(:) = 0 ! timings of routines

       section_blk(isec)%timing = 0
       section_blk(isec)%nroutine = 0 ! number of routines called
       section_blk(isec)%icall = 0

       section_blk(isec)%timer_routine%icall(:) = 0 ! number of routine calls
       section_blk(isec)%timer_routine%timing(:) = 0 ! timings of routines

    ENDDO

  END SUBROUTINE reset_timing_partial
  
  SUBROUTINE print_timing_partial(ifreq_l, itime, itime_start, ioout_l)

!***********************************************************************
!     
!     LBsoft subroutine for driving the timing services
!     originally written in MUPHY by S. Melchionna et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     modified by: M. Lauricella
!     last modification June 2018
!     
!***********************************************************************

INTEGER,INTENT(in) :: ifreq_l, itime, itime_start, ioout_l

CHARACTER(len=30),SAVE :: strtmp
INTEGER,SAVE :: iotempo_partial
LOGICAL,SAVE :: newjob=.true.
INTEGER :: i,isec
REAL(kind=PRC) :: total_time

IF(newjob) THEN
  iotempo_partial = new_io_unit()
  IF(myrank==0)WRITE(strtmp,'(a,i6)') 'tempo_',myrank
  CALL compress_blanks(strtmp)
  newjob=.false.
ENDIF

!OPEN(iotempo_partial,FILE=strtmp)

!DO isec=1, nsection
!   WRITE(iotempo_partial,*)  TRIM(ADJUSTL(section_blk(isec)%name_sec)),section_blk(isec)%timing
!ENDDO

!CLOSE(iotempo_partial)

IF(myrank==0) THEN
  WRITE(ioout_l,'(/a,i0,a)') '***** BLOCK TIMINGS (every ',ifreq_l,' timesteps) *****'
ENDIF

total_time = SUM(section_blk(1:nsection)%timing)

DO isec=1, nsection


   IF(myrank==0) THEN

      WRITE(ioout_l,1) 100.*section_blk(isec)%timing / MAX(1.e-6,total_time), &
                     TRIM(section_blk(isec)%name_sec), &
                     section_blk(isec)%timing,' sec/call,', &
                     section_blk(isec)%icall,' calls,', &
                     section_blk(isec)%icall / (1.*ifreq_l),' calls/step'

1     FORMAT(/'(',f5.1,'%)',1x,a,2x,g10.3,a,1x,i0,a,1x,g10.3,a)

   ENDIF

   !
   ! print statistics for block timing

   CALL print_routine_timing(section_blk(isec),ifreq_l,ioout_l)

ENDDO
!IF(myrank==0) THEN
!   WRITE(ioout_l,'(/a,f12.3/)') 'Sum of section timings:',SUM(section_cum(1:nsection)%timing)
!ENDIF

DO isec=1, nsection

   ! for block timers set accumulators to zero
   section_blk(isec)%timer_routine%timing(:) = 0
   section_blk(isec)%timer_routine%icall(:) = 0
   section_blk(isec)%timing = 0
   section_blk(isec)%icall = 0

ENDDO

END SUBROUTINE print_timing_partial
  
  SUBROUTINE print_timing_final(ifreq_l,itime,itime_start,natms_l,natms_tot_l,ioout_l)

!***********************************************************************
!     
!     LBsoft subroutine for driving the timing services
!     originally written in MUPHY by S. Melchionna et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     modified by: M. Lauricella
!     last modification June 2018
!     
!***********************************************************************

INTEGER,INTENT(in) :: ifreq_l,itime,itime_start,natms_l,natms_tot_l,ioout_l
CHARACTER(len=30),SAVE :: strtmp
INTEGER :: nmax,nmin,nave,i,j,isec,nrout,itemp,itempa(1)
REAL(kind=PRC) :: total_time,dtemp,dtempa(1)
LOGICAL,SAVE :: newjob=.true.
REAL(kind=PRC), allocatable, dimension(:) :: mytime_sec,dtemp_temp
INTEGER, allocatable, dimension(:) :: myorder_sec,tot_calls
LOGICAL :: lswapped
LOGICAL :: ltemp

IF(myrank==0) THEN
   WRITE(ioout_l,'(/a/)') '********** GLOBAL TIMINGS (full simulation) **********'
   WRITE(iotiming,'(/a/)')'********** GLOBAL TIMINGS (full simulation) **********'
ENDIF

total_time = 0
nrout = 0
DO i=1, mxroutine
   IF(timer_standalone(i)%icall<=0) CYCLE
   total_time = total_time + timer_standalone(i)%timing
   nrout = nrout + 1
ENDDO

 


CALL print_standalone(ifreq_l,ioout_l,nrout,total_time)

total_time = SUM(section_cum(1:nsection)%timing)


if(myrank==0)WRITE(ioout_l,*)
if(myrank==0)WRITE(iotiming,*)

allocate(mytime_sec(nsection))
allocate(myorder_sec(nsection))

DO i=1, nsection
  mytime_sec(i)=section_cum(i)%timing
  myorder_sec(i)=i
ENDDO


  
DO j=nsection-1,1,-1
  lswapped=.false.
  DO i=1,j
    IF(mytime_sec(i)>mytime_sec(i+1))THEN
      dtemp=mytime_sec(i)
      itemp=myorder_sec(i)
      mytime_sec(i)=mytime_sec(i+1)
      mytime_sec(i+1)=dtemp
      myorder_sec(i)=myorder_sec(i+1)
      myorder_sec(i+1)=itemp
      lswapped=.true.
    ENDIF
  ENDDO
  IF(.not.lswapped)EXIT
ENDDO


DO i=nsection,1,-1
   
   isec=myorder_sec(i)

   IF(myrank==0) THEN

      WRITE(ioout_l,1) 100.*section_cum(isec)%timing / MAX(1.e-6,total_time), &
                     TRIM(section_cum(isec)%name_sec), &
                     section_cum(isec)%timing,' time/call[s],', &
                     section_cum(isec)%icall,' calls'!, &
                     !section_cum(isec)%icall / (real(ifreq_l,kind=PRC)),' calls/step'

      WRITE(iotiming,1) 100.*section_cum(isec)%timing / MAX(1.e-6,total_time), &
                     TRIM(section_cum(isec)%name_sec), &
                     section_cum(isec)%timing,' time/call[s],', &
                     section_cum(isec)%icall,' calls'!, &
                     !section_cum(isec)%icall / (real(ifreq_l,kind=PRC)),' calls/step'

1     FORMAT('(',f5.1,'%)',1x,a,2x,g10.3,a,1x,i0,a,1x,g10.3,a)
   ENDIF


   !
   ! print statistics for block timing
   !CALL print_routine_timing(section_cum(isec),itime-itime_start+1,ioout_l)

ENDDO
IF(myrank==0) THEN
   WRITE(ioout_l,'(/a,f12.3/)') 'Sum of section timings:',SUM(section_cum(1:nsection)%timing)
   WRITE(iotiming,'(/a,f12.3/)') 'Sum of section timings:',SUM(section_cum(1:nsection)%timing)
ENDIF

deallocate(mytime_sec)
deallocate(myorder_sec)


END SUBROUTINE print_timing_final

SUBROUTINE print_standalone(ifreq,ioout_l,nroutine,total_time)

!***********************************************************************
!     
!     LBsoft subroutine for driving the timing services
!     originally written in MUPHY by S. Melchionna et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     modified by: M. Lauricella
!     last modification June 2018
!     
!***********************************************************************

  INTEGER,INTENT(in) :: ifreq,ioout_l,nroutine
  REAL(kind=PRC),INTENT(in) :: total_time

  REAL(kind=PRC),DIMENSION(nroutine+2) :: t_routine_max,t_routine_min,t_routine_ave
  INTEGER,DIMENSION(nroutine+2) :: ishell
  REAL(kind=PRC),DIMENSION(nroutine+2) :: tt1
  REAL(kind=PRC),DIMENSION(nroutine+2) :: tt2

  REAL(kind=PRC) :: tot_section_ave,tot_section_min
  REAL(kind=PRC) :: f,ff,tres,rate
  INTEGER :: i,j

  t_routine_max(1:nroutine) = timer_standalone(1:nroutine)%timing
  t_routine_min(1:nroutine) = timer_standalone(1:nroutine)%timing
  t_routine_ave(1:nroutine) = timer_standalone(1:nroutine)%timing

  f = SUM(t_routine_ave(1:nroutine))

  ! residual time of untimed regions
  t_routine_max(nroutine+1) = total_time - f
  t_routine_min(nroutine+1) = total_time - f
  t_routine_ave(nroutine+1) = total_time - f

  IF(t_routine_max(nroutine+1) < 0) THEN
      write(ioout_l,*)myrank,'negative time for untimed routine',total_time,f 
  ENDIF

  ! total timing
  t_routine_max(nroutine+2) = total_time
  t_routine_min(nroutine+2) = total_time
  t_routine_ave(nroutine+2) = total_time


  DO j=1,nroutine+1
     ishell(j) = j
  ENDDO
  tt2(:) = t_routine_max(:)

  CALL shellsort_r(nroutine+1,tt2,ishell)

  IF(myrank==0) THEN

      write(ioout_l,'(a10,4x,a10,1x,a30,5x,20(a12))') &
              'Time Share', &
              'Section', &
              'Routine', &
              'Time/Call[s]', &
              '#Calls', &
              'MinTime[s]', &
              'MaxTime[s]'

     f = 0
     ff = 0
     DO j=nroutine+1,1,-1 ! ascending order

        i=ishell(j)

        IF(timer_standalone(i)%icall==0 .OR. timer_standalone(i)%icall==0) CYCLE

        rate = timer_standalone(i)%icall/(1.*ifreq)

        IF(timer_standalone(i)%icall<0) THEN
           rate = 0
        ENDIF

        IF(timer_standalone(i)%icall>0) ff = ff + t_routine_max(i)

!        WRITE(ioout_l,1) t_routine_max(i)*100. / t_routine_max(nroutine+2), &
!                       TRIM(YELLOW(TRIM(timer_standalone(i)%asection))), &
!                       TRIM(RED(TRIM(timer_standalone(i)%aroutine))), &
!                       t_routine_max(i) / timer_standalone(i)%icall, &
!                       timer_standalone(i)%icall, &
!!                       rate, &
!                       t_routine_min(i)/timer_standalone(i)%icall, &
!                       t_routine_max(i)/timer_standalone(i)%icall
        WRITE(ioout_l,1) t_routine_max(i)*100. / t_routine_max(nroutine+2), &
                       TRIM(timer_standalone(i)%asection), &
                       TRIM(timer_standalone(i)%aroutine), &
                       t_routine_max(i) / timer_standalone(i)%icall, &
                       timer_standalone(i)%icall, &
!                       rate, &
                       t_routine_min(i)/timer_standalone(i)%icall, &
                       t_routine_max(i)/timer_standalone(i)%icall
        WRITE(iotiming,1) t_routine_max(i)*100. / t_routine_max(nroutine+2), &
                       TRIM(timer_standalone(i)%asection), &
                       TRIM(timer_standalone(i)%aroutine), &
                       t_routine_max(i) / timer_standalone(i)%icall, &
                       timer_standalone(i)%icall, &
!                       rate, &
                       t_routine_min(i)/timer_standalone(i)%icall, &
                       t_routine_max(i)/timer_standalone(i)%icall
     ENDDO

  ENDIF

1 FORMAT(3x,'(',f5.1,'%)',2x, &
         a20, &
         a40, 2x,&
         g15.3, &
         i12, &
         g12.3, &
         g12.3)

END SUBROUTINE print_standalone

SUBROUTINE print_routine_timing(section,ifreq,ioout_l)

!***********************************************************************
!     
!     LBsoft subroutine for driving the timing services
!     originally written in MUPHY by S. Melchionna et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     modified by: M. Lauricella
!     last modification June 2018
!     
!***********************************************************************

  TYPE(section_t) :: section
  INTEGER :: ifreq,ioout_l

  REAL(kind=PRC),DIMENSION(section%nroutine+2) :: t_routine_max,t_routine_min,t_routine_ave
  INTEGER,DIMENSION(section%nroutine+2) :: ishell
  REAL(kind=PRC) :: tot_section_ave,tot_section_min
  REAL(kind=PRC) :: f,ff,tres,rate
  REAL(kind=PRC),DIMENSION(section%nroutine+2) :: tt1
  REAL(kind=PRC),DIMENSION(section%nroutine+2) :: tt2
  INTEGER :: i,j

  t_routine_max(1:section%nroutine) = &
   section%timer_routine%timing(1:section%nroutine)
  t_routine_min(1:section%nroutine) = &
   section%timer_routine%timing(1:section%nroutine)
  t_routine_ave(1:section%nroutine) = &
   section%timer_routine%timing(1:section%nroutine)

  f = SUM(t_routine_ave(1:section%nroutine))

  t_routine_max(section%nroutine+1) = section%timing - f
  t_routine_min(section%nroutine+1) = section%timing - f
  t_routine_ave(section%nroutine+1) = section%timing - f

  t_routine_max(section%nroutine+2) = section%timing
  t_routine_min(section%nroutine+2) = section%timing
  t_routine_ave(section%nroutine+2) = section%timing


  section%aroutine(section%nroutine+1) = 'Untimed zones'
  section%timer_routine%icall(section%nroutine+1) = -1

  DO j=1,section%nroutine+1
     ishell(j) = j
  ENDDO
  tt2(:) = t_routine_max(:)
  CALL shellsort_r(section%nroutine+1,tt2,ishell)

  IF(myrank==0) THEN

     f = 0
     ff = 0
     DO j=section%nroutine+1,1,-1 ! ascending order

        i=ishell(j)

        IF(section%timer_routine%icall(i)==0 .OR. section%icall==0) CYCLE

        section%timer_routine%icall(i) = ABS(section%timer_routine%icall(i))

        rate = section%timer_routine%icall(i)/(1.*ifreq)

        IF(section%timer_routine%icall(i)<0) THEN
           rate = 0
        ENDIF

        f = f + t_routine_max(i)*100. / t_routine_max(section%nroutine+2)

        IF(section%timer_routine%icall(i)>0) ff = ff + t_routine_max(i)

        WRITE(ioout_l,1) t_routine_max(i)*100. / t_routine_max(section%nroutine+2), &
                       section%aroutine(i), &
                       t_routine_max(i) / section%timer_routine%icall(i), &
                       section%timer_routine%icall(i), &
                       rate, &
                       t_routine_min(i)/section%timer_routine%icall(i), &
                       t_routine_max(i)/section%timer_routine%icall(i)
        WRITE(iotiming,1) t_routine_max(i)*100. / t_routine_max(section%nroutine+2), &
                       section%aroutine(i), &
                       t_routine_max(i) / section%timer_routine%icall(i), &
                       section%timer_routine%icall(i), &
                       rate, &
                       t_routine_min(i)/section%timer_routine%icall(i), &
                       t_routine_max(i)/section%timer_routine%icall(i)
     ENDDO

  ENDIF

1 FORMAT(3x,'(',f5.1,'%)',1x, &
         a20,1x, &
         g10.3,1x,'sec/call,', &
         i10,1x,'calls,' &
         g10.3,1x,'calls/step,',1x, &
         '(min/max ',g10.3,',',g10.3,' sec)')

END SUBROUTINE print_routine_timing

     SUBROUTINE shellsort_r(n,rlist,ilist)

!***********************************************************************
!     
!     LBsoft subroutine for driving the timing services
!     originally written in MUPHY by S. Melchionna et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     modified by: M. Lauricella
!     last modification June 2018
!     
!***********************************************************************

       IMPLICIT NONE

       REAL(kind=PRC) :: rlist(n)
       INTEGER :: ilist(n)
       REAL(kind=PRC) :: rimax
       INTEGER :: n,nl,nn,i,ix,j,iimax

       !
       !     set up sort

       IF(n.GT.1) THEN

          !     number of lists
          nl = n/2

          !     iterate shell sort

10        DO nn = 1,nl
             !
             !     begin insertion sort on nnth list

             DO i = nn+nl,n,nl

                rimax = rlist(i)
                iimax = ilist(i)
                ix = i
                !
                !     find location for insertion

                j = i
100             j = j-nl

                IF(j.LT.1) GOTO 110

                IF (rlist(j).GT.rimax) THEN

                   ix = j

                ELSE

                   j =1

                ENDIF

                GOTO 100
110             CONTINUE

                !
                !     insert in index array

                DO j = i,ix+nl,-nl

                   rlist(j) = rlist(j-nl)
                   ilist(j) = ilist(j-nl)

                ENDDO

                rlist(ix) = rimax
                ilist(ix) = iimax

             ENDDO

          ENDDO

          nl = nl/2
          IF(nl.GT.0) GOTO 10

       ENDIF

     END SUBROUTINE shellsort_r
  
SUBROUTINE start_timing2(secname,subname)

!***********************************************************************
!     
!     LBsoft subroutine for driving the timing services
!     originally written in MUPHY by S. Melchionna et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     modified by: M. Lauricella
!     last modification June 2018
!     
!***********************************************************************

    CHARACTER(len=*) :: secname,subname
    INTEGER :: i,j
!$omp master
    j = INDEX(subname,'(')
    IF(j>0) THEN
       CALL start_timing(secname,subname(:j-1))
    ELSE
       CALL start_timing(secname,subname)
    ENDIF
!$omp end master
  END SUBROUTINE start_timing2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE end_timing2(secname,subname)
  
!***********************************************************************
!     
!     LBsoft subroutine for driving the timing services
!     originally written in MUPHY by S. Melchionna et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     modified by: M. Lauricella
!     last modification June 2018
!     
!***********************************************************************

    CHARACTER(len=*) :: secname,subname
    INTEGER :: i,j
!$omp master
    j = INDEX(subname,'(')
    IF(j>0) THEN
       CALL end_timing(secname,subname(:j-1))
    ELSE
       CALL end_timing(secname,subname)
    ENDIF
!$omp end master
  END SUBROUTINE end_timing2

  SUBROUTINE start_timing(secname,subname)
  
!***********************************************************************
!     
!     LBsoft subroutine for driving the timing services
!     originally written in MUPHY by S. Melchionna et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     modified by: M. Lauricella
!     last modification June 2018
!     
!***********************************************************************

    CHARACTER(len=*) :: secname,subname
    INTEGER :: isub,j,isec,sec_ilenght,sub_ilenght,myend,istat
    CHARACTER(len=FUNCTION_NAME_LEN) :: reconame
    CHARACTER(len=FUNCTION_NAME_LEN) :: sec_temp,sub_temp
    
    sec_ilenght=LEN_TRIM(secname)
    sub_ilenght=LEN_TRIM(subname) 
    sec_temp=repeat(' ',FUNCTION_NAME_LEN)
    sec_temp(1:sec_ilenght)=trim(secname)
    isec = id_timing_sec(sec_temp,nsection,section_cum(:)%name_sec)
    
    IF(isec<0) THEN ! section not in list. It is registered now.

       nsection = nsection + 1

       IF(nsection>mxsection) THEN
          write(IOOUT,*)'timing limit reached:'//TRIM(secname)
       ENDIF
       
       section_cum(nsection)%name_sec = sec_temp
       section_blk(nsection)%name_sec = sec_temp

       isec = nsection

    ENDIF
    
    sub_temp = repeat(' ',FUNCTION_NAME_LEN)
    sub_temp = trim(subname(1:sub_ilenght))
    myend= FUNCTION_NAME_LEN - sec_ilenght
    reconame = repeat(' ',FUNCTION_NAME_LEN)
    reconame = secname(1:sec_ilenght)//sub_temp(1:myend)

    isub = id_timing_sub(reconame)

    IF(isub<0) THEN ! routine not in list. It is registered now.

       nroutine_global = nroutine_global + 1

       section_cum(isec)%nroutine = section_cum(isec)%nroutine + 1
       section_blk(isec)%nroutine = section_blk(isec)%nroutine + 1

       IF(section_cum(isec)%nroutine>mxroutine) THEN
          write(IOOUT,*)'timing limit reached:'//TRIM(subname)
          stop
       ENDIF

       section_cum(isec)%aroutine(section_cum(isec)%nroutine) = subname
       section_blk(isec)%aroutine(section_blk(isec)%nroutine) = subname

       section_cum(isec)%timer_routine%icall(section_cum(isec)%nroutine) = 0
       section_blk(isec)%timer_routine%icall(section_blk(isec)%nroutine) = 0

       isub = nroutine_global

       id_routine_global(nroutine_global) = nroutine_global
       routine_name_global(nroutine_global) = reconame
       !write(IOOUT,'(a,i2,a,i2,3a)')'sec ',isec,' sub ',isub,' name: "', &
       ! trim(reconame),'"'

    ENDIF

    IF(timing_routine(isub)>ZERO) THEN
       write(IOOUT,*)'timer already running for routine ',TRIM(subname)
       stop
    ENDIF


    section_cum(isec)%timer_routine%icall(isub) = &
     section_cum(isec)%timer_routine%icall(isub) + 1
    section_blk(isec)%timer_routine%icall(isub) = &
     section_blk(isec)%timer_routine%icall(isub) + 1
     
    section_cum(isec)%icall = section_cum(isec)%icall + 1
    section_blk(isec)%icall = section_blk(isec)%icall + 1

    lactive_routine = .true.
    active_level = active_level + 1
    if(.not. ltiming_routine_cuda_i(isub))then
      ltiming_routine_cuda_i(isub)=.true.
      istat=cudaEventCreate(timing_routine_cuda_i(isub))
    endif
    
    istat = cudaEventRecord(timing_routine_cuda_i(isub), 0)
    !istat = cudaDeviceSynchronize
    !timing_routine(isub)=current_time()
    
    RETURN

  END SUBROUTINE start_timing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE end_timing(secname,subname)
  
!***********************************************************************
!     
!     LBsoft subroutine for driving the timing services
!     originally written in MUPHY by S. Melchionna et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     modified by: M. Lauricella
!     last modification June 2018
!     
!***********************************************************************

    CHARACTER(len=*) :: secname,subname
    INTEGER :: isub,j,isec,sec_ilenght,sub_ilenght,myend
    REAL(kind=PRC) :: ttt
    CHARACTER(len=FUNCTION_NAME_LEN) :: reconame
    CHARACTER(len=FUNCTION_NAME_LEN) :: sec_temp,sub_temp
    
    integer :: istat
    
    
    sec_ilenght=LEN_TRIM(secname)
    sub_ilenght=LEN_TRIM(subname) 
    sec_temp=repeat(' ',FUNCTION_NAME_LEN)
    sec_temp(1:sec_ilenght)=trim(secname)
    
    isec = id_timing_sec(sec_temp,nsection,section_cum(:)%name_sec)

    IF(isec<0) THEN
       write(IOOUT,*)'end_timing: section:'//TRIM(sec_temp)//' not started for timing'
       stop
    ENDIF
    
    sub_temp = repeat(' ',FUNCTION_NAME_LEN)
    sub_temp = trim(subname(1:sub_ilenght))
    myend= FUNCTION_NAME_LEN - sec_ilenght
    reconame = repeat(' ',FUNCTION_NAME_LEN)
    reconame = secname(1:sec_ilenght)//sub_temp(1:myend)

    isub = id_timing_sub(reconame)

    IF(isub<0) THEN ! should never happen that the routine is not found
       write(IOOUT,*)'end_timing: section:'//TRIM(sec_temp)// &
        ' routine not registered : '//TRIM(sub_temp)
       stop
    ENDIF
    
    !ttt = current_time()
    if(.not. ltiming_routine_cuda_e(isub))then
      ltiming_routine_cuda_e(isub)=.true.
      istat=cudaEventCreate(timing_routine_cuda_e(isub))
    endif
    
    istat = cudaEventRecord(timing_routine_cuda_e(isub), 0)
    istat = cudaDeviceSynchronize
    istat = cudaEventSynchronize(timing_routine_cuda_e(isub))
    istat = cudaEventElapsedTime(ttt,timing_routine_cuda_i(isub), &
     timing_routine_cuda_e(isub))
    ttt=ttt/1000.0

    IF(ttt< ZERO ) THEN
       write(IOOUT,*)'Warning. Negative timing for routine ', &
        TRIM(reconame),timing_routine(isub),ttt
#ifdef WALLCLOCK
       if(abs(ttt)>(clock_huge*0.5d0))then
         timing_routine(isub) = ttt + clock_huge
         if(timing_routine(isub)<0)timing_routine(isub) = ZERO
       else
         timing_routine(isub) = ZERO
       endif
#else
       timing_routine(isub) = ZERO
#endif
    else
      timing_routine(isub) = ttt
    ENDIF

    ! accumulate times

    section_cum(isec)%timer_routine%timing(isub) =  &
     section_cum(isec)%timer_routine%timing(isub) + timing_routine(isub)
    section_blk(isec)%timer_routine%timing(isub) =  &
     section_blk(isec)%timer_routine%timing(isub) + timing_routine(isub)

    timer_standalone(isub)%icall = timer_standalone(isub)%icall + 1
    timer_standalone(isub)%timing = timer_standalone(isub)%timing + &
     timing_routine(isub)
    timer_standalone(isub)%aroutine = TRIM(subname)
    timer_standalone(isub)%asection = TRIM(secname)
    
    section_cum(isec)%timing = section_cum(isec)%timing + timing_routine(isub)
    section_blk(isec)%timing = section_blk(isec)%timing + timing_routine(isub)


    timing_routine(isub) = ZERO ! -99.0

    lactive_routine = .false.
    active_level = active_level - 1

  END SUBROUTINE end_timing
  
  FUNCTION id_timing_sec(secname,n,atiming_l)
    CHARACTER(len=FUNCTION_NAME_LEN) :: secname
    INTEGER :: n
    CHARACTER(len=FUNCTION_NAME_LEN) :: atiming_l(:)
    INTEGER :: id_timing_sec
    INTEGER :: i

    id_timing_sec = -99
    DO i=1, n
       IF(trim(atiming_l(i))==trim(secname)) THEN
          id_timing_sec = i
          RETURN
       ENDIF
    ENDDO

  END FUNCTION id_timing_sec
  
  FUNCTION id_timing_sub(subname)
  
!***********************************************************************
!     
!     LBsoft subroutine for driving the timing services
!     originally written in MUPHY by S. Melchionna et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     modified by: M. Lauricella
!     last modification June 2018
!     
!***********************************************************************

    CHARACTER(len=FUNCTION_NAME_LEN) :: subname
    INTEGER :: id_timing_sub
    INTEGER :: i,j,isec,isub

    isec = -99
    id_timing_sub = -99

    DO i=1, nroutine_global
       IF(trim(subname) == trim(routine_name_global(i))) THEN
          id_timing_sub = id_routine_global(i)
          RETURN
       ENDIF
    ENDDO

  END FUNCTION id_timing_sub
  
  FUNCTION current_time() RESULT(myout)
  REAL(kind=PRC) :: myout
  INTEGER, SAVE :: clock_start, clock_stop,c1
  LOGICAL, SAVE :: newjob=.true.

  IF(newjob) THEN
  
     CALL system_clock(clock_stop, clock_rate_param, clock_max_param)
     clock_huge = &
      real(clock_max_param,kind=PRC)/real(clock_rate_param,kind=PRC)
     newjob=.false.
  ENDIF
  
  !CALL SYSTEM_CLOCK(c1)
  !out=real(c1,kind=PRC)/real(clock_rate,kind=PRC)

  
#ifdef WALLCLOCK

#ifdef MPI
    myout=MPI_Wtime()
#else
    myout=wtime()
#endif

#else


    CALL CPU_TIME(myout)
#endif

  
  RETURN
  
  END FUNCTION current_time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE compress_blanks(text)
  
!***********************************************************************
!     
!     LBsoft subroutine for driving the timing services
!     originally written in MUPHY by S. Melchionna et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     modified by: M. Lauricella
!     last modification June 2018
!     
!***********************************************************************

  ! compress all single blanks

    IMPLICIT NONE

    CHARACTER (len=*) :: text
    INTEGER :: i

    DO
      i = INDEX(TRIM(text), " ")
      IF(i==0) EXIT
      text(i:) = text(i+1:)
    ENDDO

  END SUBROUTINE compress_blanks
  
  function wtime()

!***********************************************************************
!     
!     LBsoft subroutine for computing the wall-clock time
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     modified by: M. Lauricella
!     last modification June 2018
!     
!***********************************************************************
  
  implicit none

  integer :: clock_max
  integer :: clock_rate
  integer :: clock_reading
  real(kind=PRC) wtime

  call system_clock ( clock_reading, clock_rate, clock_max )

  wtime = real ( clock_reading, kind = PRC ) &
        / real ( clock_rate, kind = PRC )

  return
  
 end function wtime
 
 subroutine print_memory_registration(iu,mybanner,mymemory)
 
!***********************************************************************
!     
!     LBsoft subroutine for printing the memory registration
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: iu
  character(len=*), intent(in) :: mybanner
  real(kind=PRC), intent(in) :: mymemory
  real(kind=PRC) :: mymemoryt,mymega,mykilo
  
  character(len=12) :: r_char,r_char2
  
  character(len=*),parameter :: of='(a)'
  
  real(kind=PRC), parameter :: convert1=real(1024.d0,kind=PRC)
  real(kind=PRC), parameter :: convert2=real(1000.d0,kind=PRC)
  
  mymemoryt=mymemory/convert1
  mymega=floor(mymemoryt)
  mykilo=(mymemoryt-mymega)*convert2
  
  if(myrank/=0)return
  write(iu,of)"                                                                               "
  write(iu,of)"********************************MEMORY MONITOR*********************************"
  write(iu,of)"                                                                               "
  write (r_char,'(i12)')nint(mymemoryt)
  write (r_char2,'(i12)')nint(mykilo)
  write(iu,'(6a)')trim(mybanner)," = ",trim(adjustl(r_char)),".",trim(adjustl(r_char2))," (mb)"
  write(iu,of)"                                                                               "
  write(iu,of)"*******************************************************************************"
  write(iu,of)"                                                                               "
  
  return
  
 end subroutine print_memory_registration
 
 subroutine print_memory_registration_cuda(iu,mybanner,mybanner2,&
  mymemory,totmem)
 
!***********************************************************************
!     
!     LBcuda subroutine for printing the memory registration
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification April 2022
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: iu
  character(len=*), intent(in) :: mybanner,mybanner2
  real(kind=PRC), intent(in) :: mymemory,totmem
  
  character(len=12) :: r_char,r_char2
  
  character(len=*),parameter :: of='(a)'
  
  
  
  if(myrank/=0)return
  write (r_char,'(f12.4)')mymemory
  write (r_char2,'(f12.4)')totmem
  write(iu,of)"                                                                               "
  write(iu,of)"******************************GPU MEMORY MONITOR*******************************"
  write(iu,of)"                                                                               "
  write(iu,'(4a)')trim(mybanner)," = ",trim(adjustl(r_char))," (GB)"
  write(iu,'(4a)')trim(mybanner2)," = ",trim(adjustl(r_char2))," (GB)"
  write(iu,of)"                                                                               "
  write(iu,of)"*******************************************************************************"
  write(iu,of)"                                                                               "
  
  return
  
 end subroutine print_memory_registration_cuda
 
 subroutine get_memory_cuda(fout,fout2)

!***********************************************************************
!     
!     LBsoft subroutine for register the memory usage
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     modified by: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************  
  
  implicit none
  
  real(kind=PRC), intent(out) :: fout,fout2
  real(kind=PRC) :: myd(2),myd2(2)
  integer(kind=cuda_count_kind) :: free, total
  integer :: istat
  
  istat = cudaMemGetInfo( free, total )
  fout = real(total-free,kind=4)/(1024.0**3.0)
  fout2 = real(total,kind=4)/(1024.0**3.0)
  
#ifdef MPI
  myd(1)=fout
  myd(2)=fout2
  call MPI_ALLREDUCE(myd,myd2,2,MPI_REAL, &
     MPI_SUM,MPI_COMM_WORLD,ier)
  fout=myd2(1)
  fout2=myd2(2)
#endif  
  
  return
  
 end subroutine get_memory_cuda
 
  subroutine get_memory(fout)
 
!***********************************************************************
!     
!     LBsoft subroutine for register the memory usage
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     modified by: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  use iso_c_binding
  
  implicit none
  
  real(kind=PRC), intent(out) :: fout
  real(kind=PRC) :: myd(1),myd2(1)

  fout = real( get_mem() ,kind=PRC)
  
#ifdef MPI
  myd(1)=fout
  call MPI_ALLREDUCE(myd,myd2,1,MPI_REAL, &
     MPI_SUM,MPI_COMM_WORLD,ier)
  fout=myd2(1)
#endif  
  
  return
  
 end subroutine get_memory
  
 end module profiling_m
