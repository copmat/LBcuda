subroutine initialize_fluids
  implicit none
  integer :: idistselect = 1
 

  select case(idistselect)
  case(1)
    call set_random_dens_fluids
  case(2)
    call set_uniform_dens_fluids
  case default
    call set_initial_dens_fluids
  end select
  
  call driver_set_initial_pop_fluids
 end subroutine initialize_fluids


 subroutine set_random_dens_fluids
  implicit none
  integer :: i,j,k
  
  if(linit_seed)then

    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          rhoR(i,j,k)=meanR+stdevR*gauss()
        enddo
      enddo
    enddo

  else

    do k=minz,maxz
      do j=miny,maxy
        do i=minx,maxx
          rhoR(i,j,k)=meanR+stdevR*gauss_noseeded(i,j,k,1)
        enddo
      enddo
    enddo

  endif
 end subroutine set_random_dens_fluids


 subroutine set_initial_pop_fluids(rho, pops)
  implicit none
  real(kind=PRC), allocatable, dimension(:,:,:)  :: rho
  real(kind=PRC), allocatable, dimension(:,:,:,:)  :: pops
  integer :: i,j,k,l
  
  do k=minz-1,maxz+1
    do j=miny-1,maxy+1
      do i=minx-1,maxx+1
        fluidsub(i,j,k,0)=equil_pop00(rhosub(i,j,k),usub(i,j,k), &
         vsub(i,j,k),wsub(i,j,k))
          enddo
      enddo
    enddo
 end subroutine set_initial_pop_fluids
  

 function gauss()
  implicit none
  real(kind=PRC) :: gauss
  real(kind=PRC) :: dtemp1,dtemp2
  real(kind=PRC), parameter :: mylimit=1.e-50

  
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
  real(kind=PRC) :: gauss_noseeded
  real(kind=PRC) :: dtemp1,dtemp2
  real(kind=PRC), parameter :: mylimit=1.e-50
  
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
  real(kind=PRC) :: rand_noseeded
  
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
  if(uni.lt.0.0)uni=uni+1.0
  csub=c-cd
  if(csub.lt.0.0)csub=csub+cm
  uni=uni-csub
  if(uni.lt.0.0)uni=uni+1.0
  rand_noseeded=real(uni,kind=PRC)
 end function rand_noseeded

