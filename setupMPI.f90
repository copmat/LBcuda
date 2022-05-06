
      subroutine setupMPI(nz)

      use dimensions_m
      use cudafor

      use mpi
      implicit none
      integer(4), intent(in) :: nz
      character*15 file_name5

      integer:: i, istat, ndev, num_gpu

!      call mpi_init(ierr)
!      call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)
!      call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
!
! building virtual topology
!
      rreorder=.false.
!
      Tperiodic(1) = xperiodic !.true.
      Tperiodic(2) = yperiodic !.true.
      Tperiodic(3) = zperiodic !.true.
!
      prgrid(1) = 1
      prgrid(2) = 1
      prgrid(3) = nprocs
!
      call MPI_cart_create(MPI_COMM_WORLD, mpid, prgrid, &
                              Tperiodic,rreorder,lbecomm,ierr)

      call MPI_comm_rank(lbecomm, myrank, ierr)
      call MPI_cart_coords(lbecomm, myrank, mpid, &
                            mpicoords, ierr)
!
      offset(1) = mpicoords(1)*nx
      offset(2) = mpicoords(2)*ny
      offset(3) = mpicoords(3)*nz
!
      call mpi_barrier(lbecomm,ierr)
!
! mpidata type (to remove for performance)
! xy plane is contiguous arrays (stride.eq.1)
      call MPI_type_contiguous((ny+2)*(nx+2), MPI_REAL, xyplane, ierr)
      call MPI_type_commit(xyplane,ierr)
#ifdef MPI_DEBUG
      if(myrank.eq.0) then
         write(6,*) 'MPI-DEBUG) INFO: xyplane (KB)-->', (ny+2)*(nx+2) *4 / 1024
      endif
#endif
      
      call MPI_type_contiguous((ny+2*nbuff)*(nx+2*nbuff)*nbuff, MPI_REAL, xyplanebuff, ierr)
      call MPI_type_commit(xyplanebuff,ierr)
#ifdef MPI_DEBUG
      if(myrank.eq.0) then
         write(6,*) 'MPI-DEBUG) INFO: xyplanebuff (KB)-->', (ny+2*nbuff)*(nx+2*nbuff)*nbuff *4 / 1024
      endif
#endif

      call MPI_type_contiguous((ny+2*nbuff)*(nx+2*nbuff)*nbuff, MPI_INTEGER1, xyplane_int1, ierr)
      call MPI_type_commit(xyplane_int1,ierr)
#ifdef MPI_DEBUG
      if(myrank.eq.0) then
         write(6,*) 'MPI-DEBUG) INFO: xyplane_int1 (KB)-->', (ny+2*nbuff)*(nx+2*nbuff)*nbuff / 1024
      endif
#endif

! xy plane is contiguous arrays 
      call MPI_type_contiguous((ny+2*nbuff)*(nx+2*nbuff)*nbuff*3, MPI_REAL, xyplane3, ierr)
      call MPI_type_commit(xyplane3,ierr)
#ifdef MPI_DEBUG
      if(myrank.eq.0) then
         write(6,*) 'MPI-DEBUG) INFO: xyplane3 (KB)-->', (ny+2*nbuff)*(nx+2*nbuff)*nbuff*3 *4 / 1024
      endif
#endif

! xy plane is contiguous arrays 
      call MPI_type_contiguous((ny+2*nbuff)*(nx+2*nbuff)*numscp, MPI_REAL, xyplanesvar, ierr)
      call MPI_type_commit(xyplanesvar,ierr)
#ifdef MPI_DEBUG
      if(myrank.eq.0) then
         write(6,*) 'MPI-DEBUG) INFO: xyplanesvar (KB)-->', (ny+2*nbuff)*(nx+2*nbuff)*numscp *4 / 1024
      endif
#endif

      call MPI_type_contiguous((ny+2)*(nx+2)*numscp, MPI_REAL, xyplanescp, ierr)
      call MPI_type_commit(xyplanescp,ierr)
#ifdef MPI_DEBUG
      if(myrank.eq.0) then
         write(6,*) 'MPI-DEBUG) INFO: xyplanescp (KB)-->', (ny+2)*(nx+2)*numscp *4 / 1024
      endif
#endif

! x dir  & y dir
      call MPI_cart_shift(lbecomm, 0, 1, rear(2), front(2), ierr)
      call MPI_cart_shift(lbecomm, 1, 1, left(2), right(2), ierr)
      call MPI_cart_shift(lbecomm, 2, 1, down(2), up(2), ierr)

      call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, &
                               MPI_INFO_NULL, localcomm, ierr)
      call MPI_comm_size(localcomm, ndev, ierr)
      call MPI_Comm_rank(localcomm, mydev, ierr)

      istat = cudaGetDeviceCount(num_gpu)
      mydev = mod(mydev, num_gpu)
      ! mydev = 0
      istat = cudaSetDevice(mydev)
      if (istat/=0) then
            write(*,*) 'status after cudaSetDevice:', cudaGetErrorString(istat)
            write(*,*) 'Exiting ....'
            call mystop
      endif

      offset_d = offset
#ifdef MPI_DEBUG
      write(6,*) 'MPI-DEBUG',myrank,'CPU offset = ', offset
      write(6,*) 'MPI-DEBUG',myrank,':tot  tasks ', nprocs
      write(6,*) 'MPI-DEBUG',myrank,':Device     ', mydev, ndev
      write(6,*) 'MPI-DEBUG',myrank,':num_gpu    ', num_gpu
      write(6,*) 'MPI-DEBUG',myrank,':task mpi-z ', mpicoords(3)
      write(6,*) 'MPI-DEBUG',myrank,':up    task ', up(2)
      write(6,*) 'MPI-DEBUG',myrank,':down  task ', down(2)
#endif
      end subroutine setupMPI
      
      subroutine MY_CART_SHIFT(mycomm_cart,myndims,mydims,myperiods, &
      mycoords,myshift,getfrom,giveto,lfromdo,ltodo,lfromloc,ltoloc, &
      ldoipdc)
      
      use dimensions_m
      use mpi
      
      implicit none

      integer, intent(in) :: mycomm_cart,myndims
      integer, intent(in), dimension(myndims) :: mydims,mycoords,myshift
      logical, intent(in), dimension(myndims) :: myperiods
      integer, intent(out) :: getfrom,giveto
      logical, intent(out) :: lfromdo,ltodo
      !if true the bounce back condition should be treated in the same rank
      logical, intent(out) :: lfromloc,ltoloc
      !if true the periodic condition should be treated in the same rank
      logical, intent(out) :: ldoipdc
      
      integer :: i,ierrs
      integer, dimension(myndims) :: cordfrom,cordto
      logical, dimension(myndims) :: mytoshift,myfromshift
      logical :: ldoto,ldofrom, lfromipdc,ltoipdc
      
      mytoshift(1:myndims)=.true.
      myfromshift(1:myndims)=.true.
      cordto(1:myndims)=mycoords(1:myndims)+myshift(1:myndims)
      cordfrom(1:myndims)=mycoords(1:myndims)-myshift(1:myndims)
      do i=1,3
        !if I am moving and I am not periodic
        if(myshift(i)/=0 .and. (.not. myperiods(i)))then
          !If I am going out do not nothing
          if(cordto(i)>=mydims(i) .or. cordto(i)<0)then
            mytoshift(i)=.false.
          endif
          !If I am going out do not nothing
          if(cordfrom(i)>=mydims(i) .or. cordfrom(i)<0)then
            myfromshift(i)=.false.
          endif
        endif
      enddo
      
      ldoto=.false.
      if(all(mytoshift))ldoto=.true.
      
      ldofrom=.false.
      if(all(myfromshift))ldofrom=.true.
      
      !write(6,*)'ldoto ',ldoto
      !write(6,*)'ldofrom ',ldofrom
      lfromloc=.false.
      ltoloc=.false.
      lfromdo=.false.
      ltodo=.false.
      lfromipdc=.false.
      ltoipdc=.false.
      if(ldoto)then
        call MPI_Cart_rank(mycomm_cart,cordto,giveto,ierrs)
        if(myrank.ne.giveto)then
          ltodo=.true.
        else
          !if I should send to me so the problem is internal periodic
          ltoipdc=.true.
        endif
      else
        giveto=myrank
        ltoloc=.true.
      endif
      
      if(ldofrom)then
        call MPI_Cart_rank(mycomm_cart,cordfrom,getfrom,ierrs)
        if(myrank.ne.getfrom)then
          lfromdo=.true.
        else
          !if I should send to me so the problem is internal periodic
          lfromipdc=.true.
        endif
      else
        getfrom=myrank
        lfromloc=.true.
      endif
      
      ldoipdc=.false.
      if(ltoipdc)then
        if(lfromipdc)then
          ldoipdc=.true.
        else
          write(6,*)'WARNING in MY_CART_SHIFT ',ltoipdc,lfromipdc
        endif
      endif
        
      
      
     end subroutine

      subroutine mystop
       use mpi
       implicit none
       integer:: ierr
        
#ifdef SERIAL
       stop
#else
       call mpi_abort(MPI_COMM_WORLD, 10, ierr)
#endif
      end subroutine mystop
      
      subroutine open_file_vtk_par(iotest,nn,myname,e_io)
 
!***********************************************************************
!     
!     LBsoft subroutine for opening the vtk legacy file
!     in parallel IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  use mpi
  implicit none
  
  integer :: iotest,nn
  character(len=nn) :: myname
  integer, intent(out) :: e_io
  
  call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(myname), &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL,iotest,e_io)
  return
  
 endsubroutine open_file_vtk_par
 
 subroutine print_header_vtk_par(iotest,offsetsub,nn,header,e_io)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing the header part of
!     in VTK legacy file in parallel IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  use mpi
  implicit none
  
  integer, intent(in) :: iotest,offsetsub,nn
  character(len=500) :: header
  integer, intent(out) :: e_io
  
  integer(kind=MPI_OFFSET_KIND) :: myoffset
  
  myoffset=int(offsetsub,kind=MPI_OFFSET_KIND)
  
  call MPI_File_write_at(iotest,myoffset,header(1:nn),nn, &
   MPI_CHARACTER,MPI_STATUS_IGNORE,e_io)
  
  return
  
 endsubroutine print_header_vtk_par
 
 subroutine print_footer_vtk_par(iotest,offsetsub,footer,e_io)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing the footer part of
!     in VTK legacy file in parallel IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  use mpi
  implicit none
  
  integer, intent(in) :: iotest,offsetsub
  character(len=30) :: footer
  integer, intent(out) :: e_io
  
  integer(kind=MPI_OFFSET_KIND) :: myoffset
  
  myoffset=int(offsetsub,kind=MPI_OFFSET_KIND)
  
  call MPI_File_write_at(iotest,myoffset,footer(1:30),30, &
   MPI_CHARACTER,MPI_STATUS_IGNORE,e_io)
  
  return
  
 end subroutine print_footer_vtk_par
 
 subroutine close_file_vtk_par(iotest,e_io)
 
!***********************************************************************
!     
!     LBsoft subroutine for closing the vtk legacy file
!     in parallel IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  use mpi
  implicit none
  
  integer :: iotest
  integer, intent(out) :: e_io
  
  call MPI_FILE_CLOSE(iotest, e_io)
  
  return
  
 endsubroutine close_file_vtk_par
 
 
 subroutine print_binary_1d_vtk_col(iotest,headoff,nbyte,nn,myrank, &
  globalDims,lDims,mymin,is,ie,myvar,e_io)
 
!***********************************************************************
!     
!     LBsoft subroutine for writing a scalar field with single
!     precision in VTK legacy binary format in parallel IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  use mpi
  implicit none
  
  integer, intent(in) :: iotest,headoff,nbyte,nn,myrank
  integer, dimension(3), intent(in) :: is,ie
  integer, intent(in), dimension(3) :: globalDims,lDims,mymin
  real(4), dimension(is(1):ie(1),is(2):ie(2),is(3):ie(3)), intent(in) :: myvar
  !dimension(is(1):ie(1),is(2):ie(2),is(3):ie(3))
  integer, intent(out) :: e_io
  
  integer(kind=MPI_OFFSET_KIND) :: myoffset
  
  integer, parameter :: byteint=4
  
  integer :: filetypesub
  integer, parameter :: nbuffsub=0
  
  integer :: memDims(3),memOffs(3),mystarts(3),imemtype
  
  myoffset=int(headoff,kind=MPI_OFFSET_KIND)
  
  mystarts = mymin - 1
  
  
  if(myrank==0)call MPI_File_write_at(iotest,myoffset,int(nbyte,kind=4),1, &
     MPI_INTEGER,MPI_STATUS_IGNORE,e_io)
  myoffset = myoffset + 1 * byteint
  
  call MPI_Type_create_subarray(3,globalDims,lDims,mystarts, &
    MPI_ORDER_FORTRAN,MPI_REAL4,filetypesub,e_io)
        
  call MPI_Type_commit(filetypesub, e_io)
   
  call MPI_File_Set_View(iotest,myoffset,MPI_REAL4,filetypesub, &
    "native",MPI_INFO_NULL,e_io)
  ! We need full local sizes: memDims
  memDims = lDims + 2*nbuffsub
  memOffs = [ nbuffsub, nbuffsub, nbuffsub ]
  

  call MPI_TYPE_CREATE_SUBARRAY(3,memDims,lDims,memOffs, &
   MPI_ORDER_FORTRAN,MPI_REAL4,imemtype,e_io)

  call MPI_TYPE_COMMIT(imemtype,e_io)

  call MPI_FILE_WRITE_ALL(iotest,myvar,1,imemtype,MPI_STATUS_IGNORE,e_io)
  
  
  return
  
 end subroutine print_binary_1d_vtk_col
 
 subroutine print_binary_3d_vtk_col(iotest,headoff,nbyte,nn,myrank, &
  globalDims,lDims,mymin,is,ie,myvar,e_io)
  
!***********************************************************************
!     
!     LBsoft subroutine for writing a vector field with single
!     precision in VTK legacy binary format in parallel IO
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
  use mpi
  implicit none
  
  integer, intent(in) :: iotest,headoff,nbyte,nn,myrank
  integer, dimension(4), intent(in) :: is,ie
  integer, intent(in), dimension(3) :: globalDims,lDims,mymin
  real(4), intent(in), &
   dimension(is(1):ie(1),is(2):ie(2),is(3):ie(3),is(4):ie(4)) :: myvar
  integer, intent(out) :: e_io
  
  integer(kind=MPI_OFFSET_KIND) :: myoffset
  
  integer, parameter :: byteint=4
  
  integer :: filetypesub
  integer, parameter :: nbuffsub=0
  
  integer :: imemtype
  
  integer, dimension(4) :: velglobalDims,velldims,velmystarts, &
   velmemDims,velmemOffs
  
  myoffset=int(headoff,kind=MPI_OFFSET_KIND)
  
  velglobalDims(1)=3
  velglobalDims(2:4)=globalDims(1:3)
  velldims(1)=3
  velldims(2:4)=lDims(1:3)
  velmystarts(1) = 0
  velmystarts(2:4) = mymin(1:3)-1
  
  if(myrank==0)call MPI_File_write_at(iotest,myoffset,int(nbyte,kind=4),1, &
     MPI_INTEGER,MPI_STATUS_IGNORE,e_io)
  myoffset = myoffset + 1 * byteint
  
  call MPI_Type_create_subarray(4,velglobalDims,velldims,velmystarts, &
    MPI_ORDER_FORTRAN,MPI_REAL4,filetypesub,e_io)
        
  call MPI_Type_commit(filetypesub, e_io)
   
  call MPI_File_Set_View(iotest,myoffset,MPI_REAL4,filetypesub, &
    "native",MPI_INFO_NULL,e_io)
  ! We need full local sizes: memDims
  velmemDims(1) = vellDims(1)
  velmemDims(2:4) = vellDims(2:4) + 2*nbuffsub
  velmemOffs = [ 0, nbuffsub, nbuffsub, nbuffsub ]

  call MPI_TYPE_CREATE_SUBARRAY(4,velmemDims,velldims,velmemOffs, &
   MPI_ORDER_FORTRAN,MPI_REAL4,imemtype,e_io)

  call MPI_TYPE_COMMIT(imemtype,e_io)

  call MPI_FILE_WRITE_ALL(iotest,myvar,1,imemtype,MPI_STATUS_IGNORE,e_io)
  
  return
  
 end subroutine print_binary_3d_vtk_col
