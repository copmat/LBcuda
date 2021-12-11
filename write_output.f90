 module write_output
    use dimensions_m
    use, intrinsic ::  iso_c_binding
    implicit none
    
    logical, save :: lelittle

    ! integer(1)  :: debugInt1(0:nx+1,0:ny+1,0:nz+1)
    ! real        :: debugReal(0:nx+1,0:ny+1,0:nz+1)
 contains
 
 subroutine test_little_endian(ltest)
 
!***********************************************************************
!     
!     LBsoft subroutine for checking if the computing architecture
!     is working in little-endian or big-endian
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none 
  integer, parameter :: ik1 = selected_int_kind(2) 
  integer, parameter :: ik4 = selected_int_kind(9) 
   
  logical, intent(out) :: ltest
   
  if(btest(transfer(int((/1,0,0,0/),ik1),1_ik4),0)) then 
    !it is little endian
    ltest=.true.
  else 
    !it is big endian
    ltest=.false.
  end if 
   
  return
   
 end subroutine test_little_endian 
 
 function space_fmtnumb(inum)
 
!***********************************************************************
!     
!     LBsoft function for returning the string of six characters 
!     with integer digits and leading spaces to the left
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none

  integer,intent(in) :: inum
  character(len=6) :: space_fmtnumb
  integer :: numdigit,irest
  real(kind=8) :: tmp
  character(len=22) :: cnumberlabel

  numdigit=dimenumb(inum)
  irest=6-numdigit
  if(irest>0)then
    write(cnumberlabel,"(a,i8,a,i8,a)")"(a",irest,",i",numdigit,")"
    write(space_fmtnumb,fmt=cnumberlabel)repeat(' ',irest),inum
  else
    write(cnumberlabel,"(a,i8,a)")"(i",numdigit,")"
    write(space_fmtnumb,fmt=cnumberlabel)inum
  endif
  
  return

 end function space_fmtnumb
 
 function space_fmtnumb12(inum)
 
!***********************************************************************
!     
!     LBsoft function for returning the string of six characters 
!     with integer digits and leading TWELVE spaces to the left
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none

  integer,intent(in) :: inum
  character(len=12) :: space_fmtnumb12
  integer :: numdigit,irest
  real(kind=8) :: tmp
  character(len=22) :: cnumberlabel

  numdigit=dimenumb(inum)
  irest=12-numdigit
  if(irest>0)then
    write(cnumberlabel,"(a,i8,a,i8,a)")"(a",irest,",i",numdigit,")"
    write(space_fmtnumb12,fmt=cnumberlabel)repeat(' ',irest),inum
  else
    write(cnumberlabel,"(a,i8,a)")"(i",numdigit,")"
    write(space_fmtnumb12,fmt=cnumberlabel)inum
  endif
  
  return

 end function space_fmtnumb12

 subroutine writeImageDataVTI(nz, fname, step, varname, myvar, textual)
  use dimensions_m
  
  implicit none
  integer(4), intent(in) :: nz
  character(len=*),intent(in) :: fname, varname
  integer,intent(in) :: step
  real(4),allocatable,intent(in) :: myvar(:,:,:)
  logical,intent(in) :: textual
  character(len=120) :: fnameFull,extent
  integer i,j,k, iotest
  integer(4)         :: length
  
  integer :: e_io
  character(1), parameter:: end_rec = char(10)
  character(1) :: string1
  character(len=500) :: headervtk
  character(len=30) :: footervtk
  integer :: iend,iini,nele,toffset,new_toffset,bytechar,byter4,byter8
  integer :: indent,ioffset,vtkoffset,ndatavtk,nn,nheadervtk,endoff
  integer :: new_myoffset,nnloc,byteint
  character(len=8) :: namevar
  integer, parameter :: ncomp=1
  integer, parameter :: ndimvtk=1
  character(len=*),parameter :: topology='ImageData' 
  real(kind=4), allocatable, dimension(:,:,:) :: service1
  
  iini=0
  bytechar=kind(end_rec)
  byteint=kind(iini)
  byter4  = 4
  byter8  = 8
  
  nn=glx*gly*glz
  
  iotest = 55
  
  if (nprocs==1) then
  
  fnameFull = 'output/' // trim(fname) // '_' // trim(write_fmtnumb(myrank)) // '.' //trim(write_fmtnumb(step)) // '.vti'
  open(unit=iotest,file=trim(fnameFull),status='replace',action='write')

  extent =  trim(write_fmtnumb(1)) // ' ' // trim(write_fmtnumb(nx)) // ' ' &
        // trim(write_fmtnumb(1)) // ' ' // trim(write_fmtnumb(ny)) // ' ' &
        // trim(write_fmtnumb(offset(3) + 1)) // ' ' // trim(write_fmtnumb(offset(3) +nz))

  write(iotest,*) '<VTKFile type="ImageData" version="1.0" byte_order="LittleEndian" >'
  write(iotest,*) ' <ImageData WholeExtent="' // trim(extent) // '" >'
  write(iotest,*) ' <Piece Extent="' // trim(extent) // '">'
  write(iotest,*) '   <PointData>'

  if (textual) then
    write(iotest,*) '    <DataArray type="Float32" Name="',trim(varname),'" format="ascii" >'

    do k=1,nz
      do j=1,ny
        do i=1,nx
          write(iotest,fmt='("     ", F20.8)') myvar(i,j,k)
        enddo
      enddo
    enddo

    write(iotest,*) '    </DataArray>'
    write(iotest,*) '   </PointData>'
    write(iotest,*) ' </Piece>'
    write(iotest,*) ' </ImageData>'
  else
    write(iotest,*) '    <DataArray type="Float32" Name="',trim(varname),'" format="appended" offset="0" />'
    write(iotest,*) '   </PointData>'
    write(iotest,*) ' </Piece>'
    write(iotest,*) ' </ImageData>'
    write(iotest,*) ' <AppendedData encoding="raw">'
    close(iotest)

    open(unit=iotest,file=trim(fnameFull),status='old',position='append', &
                access='stream',form='unformatted',action='write')
    write(iotest) '_'
    length = 4*(nx+0)*(ny+0)*(nz+0)
    write(iotest) length, myvar(1:nx, 1:ny, 1:nz)
    close(iotest)

    open(unit=iotest,file=trim(fnameFull),status='old',position='append', &
                action='write')
    write(iotest,*) ''
    write(iotest,*) ' </AppendedData>'
  endif

  write(iotest,*) '</VTKFile >'
  close(iotest)
  
  else
  
  fnameFull=repeat(' ',120)
  fnameFull = 'output/' // trim(fname) // '.' //trim(write_fmtnumb8(step)) // '.vti'
  
  call open_file_vtk_par(iotest,120,fnameFull,e_io)
  
  
  
    
  namevar=repeat(' ',8)
  namevar=trim(fname)
    
  extent =  space_fmtnumb(1) // ' ' // space_fmtnumb(glx) // ' ' &
      // space_fmtnumb(1) // ' ' // space_fmtnumb(gly) // ' ' &
      // space_fmtnumb(1) // ' ' // space_fmtnumb(glz)
    
  indent=0
    
  toffset=0
  headervtk=repeat(' ',500)
  iend=0
  
  iini=iend+1
  nele=22
  iend=iend+nele
  headervtk(iini:iend)='<?xml version="1.0"?>'//end_rec
  new_toffset=toffset
  new_toffset = new_toffset + nele * bytechar
    
  iini=iend+1
  nele=67
  iend=iend+nele
  if(lelittle)then  
    headervtk(iini:iend) = '<VTKFile type="'//trim(topology)// &
     '" version="0.1" byte_order="LittleEndian">'//end_rec
  else
    headervtk(iini:iend) = '<VTKFile type="'//trim(topology)// &
     '" version="0.1" byte_order="BigEndian">   '//end_rec
  endif
  
  new_toffset = new_toffset + 67 * bytechar
  
  indent = indent + 2
  iini=iend+1
  nele=70
  iend=iend+nele
  headervtk(iini:iend) = repeat(' ',indent)//'<'//trim(topology)//' WholeExtent="'//&
               trim(extent)//'">'//end_rec
  

  new_toffset = new_toffset + 70 * bytechar
    
  indent = indent + 2
  iini=iend+1
  nele=63
  iend=iend+nele
  headervtk(iini:iend) = repeat(' ',indent)//'<Piece Extent="'//trim(extent)//'">'//end_rec
  
  new_toffset = new_toffset + 63 * bytechar
    
  ! initializing offset pointer
  ioffset = 0 
  
  indent = indent + 2
  iini=iend+1
  nele=18
  iend=iend+nele
  headervtk(iini:iend)=repeat(' ',indent)//'<PointData>'//end_rec
  
  new_toffset = new_toffset + 18 * bytechar
  
  indent = indent + 2
  iini=iend+1
  nele=115
  iend=iend+nele
  
    
  
  write(string1,'(i1)')ncomp
   headervtk(iini:iend)=repeat(' ',indent)//'<DataArray type="Float32" Name="'// &
   namevar//'" NumberOfComponents="'//string1// '" '//&
   'format="appended" offset="'//space_fmtnumb12(ioffset)//'"/>'//end_rec
  
  new_toffset = new_toffset + 115 * bytechar
 
  
  indent = indent - 2
  iini=iend+1
  nele=19
  iend=iend+nele
  headervtk(iini:iend)=repeat(' ',indent)//'</PointData>'//end_rec
  
  new_toffset = new_toffset + 19 * bytechar
  
  
  indent = indent - 2
  iini=iend+1
  nele=13
  iend=iend+nele
  headervtk(iini:iend)=repeat(' ',indent)//'</Piece>'//end_rec
  
  
  new_toffset = new_toffset + 13 * bytechar
  
  
  indent = indent - 2
  iini=iend+1
  nele=15
  iend=iend+nele
  headervtk(iini:iend)=repeat(' ',indent)//'</'//trim(topology)//'>'//end_rec
  
  new_toffset = new_toffset + 15 * bytechar
   
  
  iini=iend+1
  nele=32
  iend=iend+nele
  headervtk(iini:iend)=repeat(' ',indent)//'<AppendedData encoding="raw">'//end_rec
  
  new_toffset = new_toffset + 32 * bytechar
  
  iini=iend+1
  nele=1
  iend=iend+nele
  headervtk(iini:iend)='_'
  
  new_toffset = new_toffset + 1 * bytechar
  
  vtkoffset=new_toffset
  toffset=new_toffset+byteint+ndimvtk*nn*byter4
  ndatavtk=ndimvtk*nn*byter4
  nheadervtk=iend
  
  if(myrank==0)then
    call print_header_vtk_par(iotest,0,nheadervtk,headervtk,e_io)  
  endif
  
  footervtk=repeat(' ',30)
  iini=0
  iend=iini
  
  iini=iend+1
  nele=1
  iend=iend+nele
  footervtk(iini:iend)=end_rec
  
  new_myoffset = toffset
  new_myoffset = new_myoffset + 1 * bytechar
  
  iini=iend+1
  nele=18
  iend=iend+nele
  footervtk(iini:iend)=repeat(' ',indent)//'</AppendedData>'//end_rec
  
  new_myoffset = new_myoffset + 18 * bytechar
  
  iini=iend+1
  nele=11
  iend=iend+nele
  footervtk(iini:iend)='</VTKFile>'//end_rec
  
  endoff=vtkoffset+ndatavtk+byteint
  if(myrank==0)call print_footer_vtk_par(iotest,endoff,footervtk,e_io)
  nnloc=nx*ny*nz
  allocate(service1(1:nx,1:ny,1:nz))
  service1(1:nx,1:ny,1:nz)=0.0
  do k=1,nz
    do j=1,ny
      do i=1,nx
        !if(myfluid(i,j,kk)==3)cycle
        service1(i,j,k)=real(myvar(i,j,k),kind=4)
      enddo
    enddo
  enddo
  
  
  call print_binary_1d_vtk_col(iotest,vtkoffset,ndatavtk,nnloc,myrank, &
     gsizes,lsizes,start_idx,lbound(service1),ubound(service1),service1,e_io)
  
  call close_file_vtk_par(iotest,e_io)
  
  deallocate(service1)
  
  endif
  
 end subroutine writeImageDataVTI
 
 

 subroutine writeImageDataVTI_3d(nz, fname, step, velsub, textual)
  use dimensions_m
  implicit none
  integer(4), intent(in) :: nz
  character(len=*),intent(in) :: fname
  integer,intent(in) :: step
  real(4),allocatable,intent(in) :: velsub(:,:,:,:)
  logical,intent(in) :: textual
  character(len=120) :: fnameFull,extent
  integer i,j,k, iotest
  integer(4)         :: length
  
  integer :: e_io
  character(1), parameter:: end_rec = char(10)
  character(1) :: string1
  character(len=500) :: headervtk
  character(len=30) :: footervtk
  integer :: iend,iini,nele,toffset,new_toffset,bytechar,byter4,byter8
  integer :: indent,ioffset,vtkoffset,ndatavtk,nn,nheadervtk,endoff
  integer :: new_myoffset,nnloc,byteint
  character(len=8) :: namevar
  integer, parameter :: ncomp=3
  integer, parameter :: ndimvtk=3
  character(len=*),parameter :: topology='ImageData' 
  real(kind=4), allocatable, dimension(:,:,:,:) :: service4
  
  iini=0
  bytechar=kind(end_rec)
  byteint=kind(iini)
  byter4  = 4
  byter8  = 8
  
  nn=glx*gly*glz
  

  iotest = 55
  
  if (nprocs==1) then
  
  fnameFull = 'output/' // trim(fname) // '_' // trim(write_fmtnumb(myrank)) // '.' //trim(write_fmtnumb(step)) // '.vti'

  open(unit=iotest,file=trim(fnameFull),status='replace',action='write')

  extent =  trim(write_fmtnumb(1)) // ' ' // trim(write_fmtnumb(nx)) // ' ' &
        // trim(write_fmtnumb(1)) // ' ' // trim(write_fmtnumb(ny)) // ' ' &
        // trim(write_fmtnumb(offset(3) + 1)) // ' ' // trim(write_fmtnumb(offset(3) + nz))

  write(iotest,*) '<VTKFile type="ImageData" version="1.0">'
  write(iotest,*) ' <ImageData WholeExtent="' // trim(extent) // '" >'
  write(iotest,*) ' <Piece Extent="' // trim(extent) // '">'
  write(iotest,*) '   <PointData>'

  if (textual) then
    write(iotest,*) '    <DataArray type="Float32" Name="vel" NumberOfComponents="3" format="ascii" >'

    do k=1,nz
      do j=1,ny
        do i=1,nx
          write(iotest,fmt='("     ", 3F20.8)') velsub(1:3,i,j,k)
        enddo
      enddo
    enddo

    write(iotest,*) '    </DataArray>'
    write(iotest,*) '   </PointData>'
    write(iotest,*) ' </Piece>'
    write(iotest,*) ' </ImageData>'
  else
    write(iotest,*) '    <DataArray type="Float32" Name="vel" NumberOfComponents="3" format="appended" offset="0" />'
    write(iotest,*) '   </PointData>'
    write(iotest,*) ' </Piece>'
    write(iotest,*) ' </ImageData>'
    write(iotest,*) ' <AppendedData encoding="raw">'
    close(iotest)

    open(unit=iotest,file=trim(fnameFull),status='old',position='append', &
                access='stream',form='unformatted',action='write')
    write(iotest) '_'
    length = 4*nx*ny*nz*3
    write(iotest) length, velsub(1:3, 1:nx,1:ny,1:nz)
    close(iotest)

    open(unit=iotest,file=trim(fnameFull),status='old',position='append', &
                action='write')
    write(iotest,*) ''
    write(iotest,*) ' </AppendedData>'
  endif

  write(iotest,*) '</VTKFile >'
  close(iotest)
  
  else
  
  fnameFull=repeat(' ',120)
  fnameFull = 'output/' // trim(fname) // '.' //trim(write_fmtnumb8(step)) // '.vti'
  
  call open_file_vtk_par(iotest,120,fnameFull,e_io)
  
  
  
    
  namevar=repeat(' ',8)
  namevar=trim(fname)
    
  extent =  space_fmtnumb(1) // ' ' // space_fmtnumb(glx) // ' ' &
      // space_fmtnumb(1) // ' ' // space_fmtnumb(gly) // ' ' &
      // space_fmtnumb(1) // ' ' // space_fmtnumb(glz)
    
  indent=0
    
  toffset=0
  headervtk=repeat(' ',500)
  iend=0
  
  iini=iend+1
  nele=22
  iend=iend+nele
  headervtk(iini:iend)='<?xml version="1.0"?>'//end_rec
  new_toffset=toffset
  new_toffset = new_toffset + nele * bytechar
    
  iini=iend+1
  nele=67
  iend=iend+nele
  if(lelittle)then  
    headervtk(iini:iend) = '<VTKFile type="'//trim(topology)// &
     '" version="0.1" byte_order="LittleEndian">'//end_rec
  else
    headervtk(iini:iend) = '<VTKFile type="'//trim(topology)// &
     '" version="0.1" byte_order="BigEndian">   '//end_rec
  endif
  
  new_toffset = new_toffset + 67 * bytechar
  
  indent = indent + 2
  iini=iend+1
  nele=70
  iend=iend+nele
  headervtk(iini:iend) = repeat(' ',indent)//'<'//trim(topology)//' WholeExtent="'//&
               trim(extent)//'">'//end_rec
  

  new_toffset = new_toffset + 70 * bytechar
    
  indent = indent + 2
  iini=iend+1
  nele=63
  iend=iend+nele
  headervtk(iini:iend) = repeat(' ',indent)//'<Piece Extent="'//trim(extent)//'">'//end_rec
  
  new_toffset = new_toffset + 63 * bytechar
    
  ! initializing offset pointer
  ioffset = 0 
  
  indent = indent + 2
  iini=iend+1
  nele=18
  iend=iend+nele
  headervtk(iini:iend)=repeat(' ',indent)//'<PointData>'//end_rec
  
  new_toffset = new_toffset + 18 * bytechar
  
  indent = indent + 2
  iini=iend+1
  nele=115
  iend=iend+nele
  
    
  
  write(string1,'(i1)')ncomp
   headervtk(iini:iend)=repeat(' ',indent)//'<DataArray type="Float32" Name="'// &
   namevar//'" NumberOfComponents="'//string1// '" '//&
   'format="appended" offset="'//space_fmtnumb12(ioffset)//'"/>'//end_rec
  
  new_toffset = new_toffset + 115 * bytechar
 
  
  indent = indent - 2
  iini=iend+1
  nele=19
  iend=iend+nele
  headervtk(iini:iend)=repeat(' ',indent)//'</PointData>'//end_rec
  
  new_toffset = new_toffset + 19 * bytechar
  
  
  indent = indent - 2
  iini=iend+1
  nele=13
  iend=iend+nele
  headervtk(iini:iend)=repeat(' ',indent)//'</Piece>'//end_rec
  
  
  new_toffset = new_toffset + 13 * bytechar
  
  
  indent = indent - 2
  iini=iend+1
  nele=15
  iend=iend+nele
  headervtk(iini:iend)=repeat(' ',indent)//'</'//trim(topology)//'>'//end_rec
  
  new_toffset = new_toffset + 15 * bytechar
   
  
  iini=iend+1
  nele=32
  iend=iend+nele
  headervtk(iini:iend)=repeat(' ',indent)//'<AppendedData encoding="raw">'//end_rec
  
  new_toffset = new_toffset + 32 * bytechar
  
  iini=iend+1
  nele=1
  iend=iend+nele
  headervtk(iini:iend)='_'
  
  new_toffset = new_toffset + 1 * bytechar
  
  vtkoffset=new_toffset
  toffset=new_toffset+byteint+ndimvtk*nn*byter4
  ndatavtk=ndimvtk*nn*byter4
  nheadervtk=iend
  
  if(myrank==0)then
    call print_header_vtk_par(iotest,0,nheadervtk,headervtk,e_io)  
  endif
  
  footervtk=repeat(' ',30)
  iini=0
  iend=iini
  
  iini=iend+1
  nele=1
  iend=iend+nele
  footervtk(iini:iend)=end_rec
  
  new_myoffset = toffset
  new_myoffset = new_myoffset + 1 * bytechar
  
  iini=iend+1
  nele=18
  iend=iend+nele
  footervtk(iini:iend)=repeat(' ',indent)//'</AppendedData>'//end_rec
  
  new_myoffset = new_myoffset + 18 * bytechar
  
  iini=iend+1
  nele=11
  iend=iend+nele
  footervtk(iini:iend)='</VTKFile>'//end_rec
  
  endoff=vtkoffset+ndatavtk+byteint
  if(myrank==0)call print_footer_vtk_par(iotest,endoff,footervtk,e_io)
  nnloc=nx*ny*nz
  allocate(service4(1:3,1:nx,1:ny,1:nz))
  
  do k=1,nz
    do j=1,ny
      do i=1,nx
        !if(myfluid(i,j,kk)==3)cycle
        service4(1:3,i,j,k)=real(velsub(1:3,i,j,k),kind=4)
      enddo
    enddo
  enddo
  
  
  call print_binary_3d_vtk_col(iotest,vtkoffset,ndatavtk,nnloc,myrank, &
     gsizes,lsizes,start_idx,lbound(service4),ubound(service4),service4,e_io)
  
  call close_file_vtk_par(iotest,e_io)
  
  deallocate(service4)
  
  endif
  
 end subroutine

 subroutine writeParticleVTK(fname, step, x_atm,v_atm, q, textual)
  use dimensions_m
  implicit none
  character(len=*),intent(in) :: fname
  integer,intent(in) :: step
  real(4), allocatable :: x_atm(:,:), v_atm(:,:), q(:,:)
  logical,intent(in) :: textual
  character(len=120) :: fnameFull
  integer iatm, iotest
  


  iotest = 55
  fnameFull = 'output/' // trim(fname) // '_' // trim(write_fmtnumb(step)) // '.vtk'
  open(unit=iotest,file=trim(fnameFull),status='replace',action='write')

  write(iotest,fmt=100) '# vtk DataFile Version 2.0'
  write(iotest,fmt=100) 'Field Emission Device - Charge Density Plot'
  write(iotest,fmt=100) 'ASCII'
  write(iotest,fmt=100) 'DATASET POLYDATA'
  write(iotest,fmt=101) 'POINTS ',numAtoms,' float'

  do iatm = 1,numAtoms
    write(iotest,108) x_atm(1, iatm), x_atm(2, iatm), x_atm(3, iatm)
  enddo

  write(iotest,fmt=100) ''
  write(iotest,fmt=102) 'POINT_DATA ', numAtoms
  write(iotest,fmt=100) 'SCALARS Type float 1'
  write(iotest,fmt=100) 'LOOKUP_TABLE default'
  
  do iatm = 1 ,numAtoms
    write(iotest,fmt=103) iatm !Only Spheres...
  enddo

  write(iotest,fmt=100) ''
  write(iotest,fmt=100) 'VECTORS vel float'
  do iatm = 1 ,numAtoms
    write(iotest,*) v_atm(1, iatm), v_atm(2, iatm), v_atm(3, iatm)
  enddo

  if (lrotate) then
    write(iotest,fmt=100) ''
    write(iotest,fmt=100) 'VECTORS Vectx float'
    do iatm = 1 ,numAtoms
      write(iotest,201)take_rotversorx(q(1,iatm),q(2,iatm),q(3,iatm),q(4,iatm))
    enddo
    
    write(iotest,fmt=100) ''
    write(iotest,fmt=100) 'VECTORS Vecty float'
    do iatm = 1 ,numAtoms
      write(iotest,201)take_rotversory(q(1,iatm),q(2,iatm),q(3,iatm),q(4,iatm))
    enddo
    
    write(iotest,fmt=100) ''
    write(iotest,fmt=100) 'VECTORS Vectz float'
    do iatm = 1 ,numAtoms
      write(iotest,201)take_rotversorz(q(1,iatm),q(2,iatm),q(3,iatm),q(4,iatm))
    enddo
  endif


  close(iotest)

  100 format (A)
  101 format (A,I12,A)
  102 format (A,I12)
  103 format (I12)
  104 format (3F20.12)
  108 format (ES16.5,ES16.5,ES16.5)
  200 format (A,4F20.8)
  201 format (4F20.8)
 end subroutine writeParticleVTK

 pure function qconj(qs)
  implicit none      
  real, intent(in), dimension(0:3) :: qs      
  real, dimension(0:3) :: qconj
  
  qconj(0) = qs(0)
  qconj(1:3) = -qs(1:3)
 end function qconj

 pure function cross(a,b) 
   implicit none
   real, dimension(3) :: cross
   real, dimension(3), intent(in) :: a, b
 
   cross(1) = a(2) * b(3) - a(3) * b(2)
   cross(2) = a(3) * b(1) - a(1) * b(3)
   cross(3) = a(1) * b(2) - a(2) * b(1)   
  end function cross

 pure function qmult(qsa,qsb)
  implicit none
  real, intent(in), dimension(0:3) :: qsa,qsb
  real, dimension(0:3) :: qmult
  
  qmult(0) = qsa(0)*qsb(0) - qsa(1)*qsb(1) - qsa(2)*qsb(2) - qsa(3)*qsb(3)
  qmult(1:3) = qsa(0)*qsb(1:3) + qsb(0)*qsa(1:3) + cross(qsa(1:3),qsb(1:3))
 end function qmult
    
 pure function qtrimult(qsa,qsb,qsc)
  implicit none
  real, intent(in), dimension(0:3) :: qsa,qsb,qsc
  real, dimension(0:3) :: qtrimult,qtemps

  qtemps = qmult(qsb,qsc)
  qtrimult = qmult(qsa,qtemps)
 end function qtrimult

 function take_rotversorx(qa,qb,qc,qd)
  implicit none
  real, intent(in) :: qa,qb,qc,qd
  real :: q0,q1,q2,q3
  real, dimension(3) :: take_rotversorx
  real, dimension(0:3), parameter :: qversor = (/ ZERO , ONE , ZERO , ZERO /)
  real, dimension(0:3) :: qtemp,qtemp2, tempconj
  
  qtemp(0) = qa
  qtemp(1) = qb
  qtemp(2) = qc
  qtemp(3) = qd
  
  tempconj = qconj(qtemp)
  qtemp2 = qtrimult(qtemp,qversor,tempconj)
  take_rotversorx = qtemp2(1:3)
 end function take_rotversorx
  
 function take_rotversory(qa,qb,qc,qd)
  implicit none
  real, intent(in) :: qa,qb,qc,qd
  real :: q0,q1,q2,q3
  real, dimension(3) :: take_rotversory
  real, dimension(0:3), parameter :: qversor = (/ ZERO , ZERO , ONE , ZERO /)
  real, dimension(0:3) :: qtemp,qtemp2,tempconj

  qtemp(0) = qa
  qtemp(1) = qb
  qtemp(2) = qc
  qtemp(3) = qd

  tempconj = qconj(qtemp)
  qtemp2 = qtrimult(qtemp,qversor,tempconj)
  take_rotversory = qtemp2(1:3)
 end function take_rotversory
  
 function take_rotversorz(qa,qb,qc,qd)
  implicit none
  real, intent(in) :: qa,qb,qc,qd
  real :: q0,q1,q2,q3
  real, dimension(3) :: take_rotversorz
  real, dimension(0:3), parameter :: qversor = (/ ZERO , ZERO , ZERO , ONE /)
  real, dimension(0:3) :: qtemp,qtemp2,tempconj
  
  qtemp(0) = qa
  qtemp(1) = qb
  qtemp(2) = qc
  qtemp(3) = qd
  
  tempconj = qconj(qtemp)
  qtemp2 = qtrimult(qtemp,qversor,tempconj)
  take_rotversorz = qtemp2(1:3)
 end function take_rotversorz


 function dimenumb(inum)
  implicit none
  integer,intent(in) :: inum
  integer :: dimenumb, i
  real :: tmp

  i=1
  tmp = real(inum)
  do
    if (tmp<10.0) exit
    i = i+1
    tmp = tmp / 10.0
  enddo

  dimenumb = i
 end function dimenumb


 function write_fmtnumb(inum)
  implicit none
  integer,intent(in) :: inum
  character(len=6) :: write_fmtnumb
  integer :: numdigit,irest
  real :: tmp
  character(len=22) :: cnumberlabel

  numdigit = dimenumb(inum)
  irest=6-numdigit
  if(irest>0)then
    write(cnumberlabel,"(a,i8,a,i8,a)")"(a",irest,",i",numdigit,")"
    write(write_fmtnumb,fmt=cnumberlabel)repeat('0',irest),inum
  else
    write(cnumberlabel,"(a,i8,a)")"(i",numdigit,")"
    write(write_fmtnumb,fmt=cnumberlabel)inum
  endif
 end function write_fmtnumb

 
 function write_fmtnumb8(inum)
  implicit none
  integer,intent(in) :: inum
  character(len=8) :: write_fmtnumb8
  integer :: numdigit,irest
  real :: tmp
  character(len=22) :: cnumberlabel

  numdigit=dimenumb(inum)
  irest=8-numdigit
  if(irest>0)then
    write(cnumberlabel,"(a,i8,a,i8,a)")"(a",irest,",i",numdigit,")"
    write(write_fmtnumb8,fmt=cnumberlabel)repeat('0',irest),inum
  else
    write(cnumberlabel,"(a,i8,a)")"(i",numdigit,")"
    write(write_fmtnumb8,fmt=cnumberlabel)inum
  endif
 end function write_fmtnumb8

 function write_fmtnumb0(inum, sz)
  implicit none
  integer,intent(in) :: inum, sz
  character(len=6) :: write_fmtnumb0
  integer :: numdigit,irest
  real :: tmp
  character(len=22) :: cnumberlabel

  numdigit = dimenumb(inum)
  irest = sz - numdigit
  if(irest>0)then
    write(cnumberlabel,"(a,i8,a,i8,a)")"(a",irest,",i",numdigit,")"
    write(write_fmtnumb0,fmt=cnumberlabel)repeat('0',irest),inum
  else
    write(cnumberlabel,"(a,i8,a)")"(i",numdigit,")"
    write(write_fmtnumb0,fmt=cnumberlabel)inum
  endif
 end function write_fmtnumb0


 subroutine moments1Fl(nz, pops, rhoR,vel)
  use dimensions_m
  implicit none
  integer(4), intent(in) :: nz
  real(4), allocatable, pinned,intent(in) :: pops(:,:,:,:)
  real(4), allocatable :: rhoR(:,:,:), vel(:,:,:,:)
  integer i,j,k, l
  real(4) locrho,invrho, locu,locv,locw
  real(4), parameter :: MINDENS = 0.0


  do k=1,nz
    do j=1,ny
      do i=1,nx
       locrho = pops(i,j,k,0) + pops(i,j,k,1) + pops(i,j,k,2) + &
        pops(i,j,k,3) + pops(i,j,k,4) + &
        pops(i,j,k,5) + pops(i,j,k,6) + pops(i,j,k,7) + &
        pops(i,j,k,8) + pops(i,j,k,9) + &
        pops(i,j,k,10) + pops(i,j,k,11) + pops(i,j,k,12) + &
        pops(i,j,k,13) + pops(i,j,k,14) + &
        pops(i,j,k,15) + pops(i,j,k,16) + &
        pops(i,j,k,17) + pops(i,j,k,18)

        invrho = ONE / locrho

        locu   = invrho * ( pops(i,j,k,1) - pops(i,j,k,2) + pops(i,j,k,7) - &
         pops(i,j,k,8) - pops(i,j,k,9) + &
         pops(i,j,k,10) + pops(i,j,k,11) - pops(i,j,k,12) - &
         pops(i,j,k,13) + pops(i,j,k,14) )

        locv    = invrho * ( pops(i,j,k,3) - pops(i,j,k,4) + pops(i,j,k,7) - &
         pops(i,j,k,8) + pops(i,j,k,9) - &
         pops(i,j,k,10) + pops(i,j,k,15) - pops(i,j,k,16) - &
         pops(i,j,k,17) + pops(i,j,k,18) )

        locw    = invrho * ( pops(i,j,k,5) - pops(i,j,k,6) + pops(i,j,k,11) - &
         pops(i,j,k,12) + pops(i,j,k,13) - &
         pops(i,j,k,14) + pops(i,j,k,15) - pops(i,j,k,16) + &
         pops(i,j,k,17) - pops(i,j,k,18) )

        if (locrho > MINDENS) then
          rhoR(i,j,k) = locrho
          vel(1,i,j,k) = locu
          vel(2,i,j,k) = locv
          vel(3,i,j,k) = locw
        else
          rhoR(i,j,k) = MINDENS
          vel(1,i,j,k) = 0.0
          vel(2,i,j,k) = 0.0
          vel(3,i,j,k) = 0.0
        endif

      enddo
    enddo
  enddo
  ! write (*,*) 1,1,1, rhoR(1,1,1), (pops(1,1,1,l), l=0,18)
 end subroutine moments1Fl

 subroutine moments2Fl(nz, pops,rhoR, popsB,rhoB, vel, phase, myfluid, flip)
  use dimensions_m
  implicit none
  integer(4), intent(in) :: nz
  real(4), allocatable, pinned,intent(in) :: pops(:,:,:,:), popsB(:,:,:,:)
  real(4), allocatable :: rhoR(:,:,:), rhoB(:,:,:), vel(:,:,:,:), phase(:,:,:)
  integer(1), allocatable, intent(in) :: myfluid(:,:,:,:)
  integer, intent(in) :: flip
  integer i,j,k, l
  real(4) locrhor,locrhob,rhosum,invrho, locu,locv,locw
  real(4), parameter :: MINDENS = 0.0


  do k=1,nz
    do j=1,ny
      do i=1,nx
        locrhor = pops(i,j,k,0) + pops(i,j,k,1) + pops(i,j,k,2) + pops(i,j,k,3) + pops(i,j,k,4) + &
          pops(i,j,k,5) + pops(i,j,k,6) + pops(i,j,k,7) + pops(i,j,k,8) + pops(i,j,k,9) + &
          pops(i,j,k,10) + pops(i,j,k,11) + pops(i,j,k,12) + pops(i,j,k,13) + pops(i,j,k,14) + &
          pops(i,j,k,15) + pops(i,j,k,16) + pops(i,j,k,17) + pops(i,j,k,18)

        locrhob = popsB(i,j,k,0) + popsB(i,j,k,1) + popsB(i,j,k,2) + popsB(i,j,k,3) + popsB(i,j,k,4) + &
          popsB(i,j,k,5) + popsB(i,j,k,6) + popsB(i,j,k,7) + popsB(i,j,k,8) + popsB(i,j,k,9) + &
          popsB(i,j,k,10) + popsB(i,j,k,11) + popsB(i,j,k,12) + popsB(i,j,k,13) + popsB(i,j,k,14) + &
          popsB(i,j,k,15) + popsB(i,j,k,16) + popsB(i,j,k,17) + popsB(i,j,k,18)

        rhosum = locrhor + locrhob
        invrho = ONE / rhosum

        locu   = invrho * ( &
         pops(i,j,k,1) - pops(i,j,k,2) + pops(i,j,k,7) - pops(i,j,k,8) - pops(i,j,k,9) + &
         pops(i,j,k,10) + pops(i,j,k,11) - pops(i,j,k,12) - pops(i,j,k,13) + pops(i,j,k,14) + &
         popsB(i,j,k,1) - popsB(i,j,k,2) + popsB(i,j,k,7) - popsB(i,j,k,8) - popsB(i,j,k,9) + &
         popsB(i,j,k,10) + popsB(i,j,k,11) - popsB(i,j,k,12) - popsB(i,j,k,13) + popsB(i,j,k,14) )

        locv    = invrho * ( &
         pops(i,j,k,3) - pops(i,j,k,4) + pops(i,j,k,7) - pops(i,j,k,8) + pops(i,j,k,9) - &
         pops(i,j,k,10) + pops(i,j,k,15) - pops(i,j,k,16) - pops(i,j,k,17) + pops(i,j,k,18) + & 
         popsB(i,j,k,3) - popsB(i,j,k,4) + popsB(i,j,k,7) - popsB(i,j,k,8) + popsB(i,j,k,9) - &
         popsB(i,j,k,10) + popsB(i,j,k,15) - popsB(i,j,k,16) - popsB(i,j,k,17) + popsB(i,j,k,18) )

        locw    = invrho * ( &
         pops(i,j,k,5) - pops(i,j,k,6) + pops(i,j,k,11) - pops(i,j,k,12) + pops(i,j,k,13) - &
         pops(i,j,k,14) + pops(i,j,k,15) - pops(i,j,k,16) + pops(i,j,k,17) - pops(i,j,k,18) + & 
         popsB(i,j,k,5) - popsB(i,j,k,6) + popsB(i,j,k,11) - popsB(i,j,k,12) + popsB(i,j,k,13) - &
         popsB(i,j,k,14) + popsB(i,j,k,15) - popsB(i,j,k,16) + popsB(i,j,k,17) - popsB(i,j,k,18) )

        if (myfluid(i,j,k, flip) == fluid_fluid .and. rhosum > MINDENS) then
          rhoR(i,j,k) = locrhor
          rhoB(i,j,k) = locrhob
          vel(1,i,j,k) = locu
          vel(2,i,j,k) = locv
          vel(3,i,j,k) = locw
          phase(i,j,k) = (locrhor - locrhob) / (locrhor + locrhob)
        else
          rhoR(i,j,k) = MINDENS
          rhoB(i,j,k) = MINDENS
          vel(1,i,j,k) = 0.0
          vel(2,i,j,k) = 0.0
          vel(3,i,j,k) = 0.0
          phase(i,j,k) = 10.0
        endif

      enddo
    enddo
  enddo
 end subroutine moments2Fl

 subroutine writeImageDataVTI_isfluid(nz, fname, step, myfluid, flip, textual, totSphere)
  use dimensions_m
  use kernels_fluid
  implicit none
  integer(4), intent(in) :: nz
  character(len=*),intent(in) :: fname
  integer,intent(in) :: step, flip, totSphere
  integer(1), allocatable, intent(in) :: myfluid(:,:,:,:)
  logical,intent(in) :: textual
  character(len=120) :: fnameFull,extent
  integer i,j,k, iotest
  integer(4)         :: length
  real(8)            :: totSumR
  integer(8)         :: totSumI, totCenters, isFlCent(3,numAtoms)

  iotest = 55

  ! fnameFull = 'output/' // 'partisfluid_' // trim(write_fmtnumb(step)) // '.txt'
  ! open(unit=iotest,file=trim(fnameFull),status='replace',action='write')

  totSumI = 0
  totCenters = 0
  do k=1,nz
    do j=1,ny
      do i=1,nx
        if (myfluid(i,j,k, flip) /= fluid_fluid) then
          totSumI = totSumI + 1
          if (myfluid(i,j,k, flip) == fluid_particleCM) then
            totCenters = totCenters + 1
            isFlCent(1,totCenters) = i
            isFlCent(2,totCenters) = j
            isFlCent(3,totCenters) = k + offset_d(3)
          endif
        endif
      enddo
    enddo
  enddo
  write(*,fmt=200) step, myrank, totSumI,totCenters
  200     format('isfluid] step=', I8, ' rank=', I6,' myfluid=', I10,' CM=', I10)

  if (nprocs==1) then
    if (myrank==0 .and. totCenters /= numAtoms) then
      write(*,fmt=201) step, totCenters, numAtoms
      201     format('isfluid] step=', I8, ' ERROR!!!!!!!!!!   particle number diff:',I8,' vs ',I8)

      if (step<2) call mystop
    endif

    if (myrank==0 .and. totSumI /= totSphere*numAtoms) then
      write(*,fmt=202) step, totSumI, totSphere*numAtoms
      202     format('isfluid] step=', I8, ' ERROR!!!!!!!!!!   particle volume diff:',I8,' vs ',I8)
      if (step<2) call mystop
    endif
  else
    ! write(*,fmt=201) step, totCenters, numAtoms
    ! write(*,fmt=202) step, totSumI, totSphere*numAtoms
  endif
  
  
  fnameFull = 'output/' // trim(fname) // '_' // trim(write_fmtnumb(myrank)) // '_' // trim(write_fmtnumb(step)) // '.vti'
  open(unit=iotest,file=trim(fnameFull),status='replace',action='write')

  extent =  trim(write_fmtnumb(1)) // ' ' // trim(write_fmtnumb(nx)) // ' ' &
        // trim(write_fmtnumb(1)) // ' ' // trim(write_fmtnumb(ny)) // ' ' &
        // trim(write_fmtnumb(offset(3) + 1)) // ' ' // trim(write_fmtnumb(offset(3) + nz))

  write(iotest,*) '<VTKFile type="ImageData" version="1.0" byte_order="LittleEndian" >'
  write(iotest,*) ' <ImageData WholeExtent="' // trim(extent) // '" >'
  write(iotest,*) ' <Piece Extent="' // trim(extent) // '">'
  write(iotest,*) '   <PointData>'

  if (textual) then
    write(iotest,*) '    <DataArray type="UInt8" Name="',trim(fname),'" format="ascii" >'

    do k=1,nz
      do j=1,ny
        do i=1,nx
          write(iotest,fmt='("     ", I3)') myfluid(i,j,k, flip)
        enddo
      enddo
    enddo

    write(iotest,*) '    </DataArray>'
    write(iotest,*) '   </PointData>'
    write(iotest,*) ' </Piece>'
    write(iotest,*) ' </ImageData>'
  else
    write(iotest,*) '    <DataArray type="UInt8" Name="',trim(fname),'" format="appended" offset="0" />'
    write(iotest,*) '   </PointData>'
    write(iotest,*) ' </Piece>'
    write(iotest,*) ' </ImageData>'
    write(iotest,*) ' <AppendedData encoding="raw">'
    close(iotest)

    open(unit=iotest,file=trim(fnameFull),status='old',position='append', &
                access='stream',form='unformatted',action='write')
    write(iotest) '_'
    length = 1*nx*ny*nz
    write(iotest) length, myfluid(1:nx,1:ny,1:nz, flip)
    close(iotest)

    open(unit=iotest,file=trim(fnameFull),status='old',position='append', &
                action='write')
    write(iotest,*) ''
    write(iotest,*) ' </AppendedData>'
  endif

  write(iotest,*) '</VTKFile >'
  close(iotest)


  ! debugInt1 = debugmkdel_d
  ! totSumI = 0.0
  ! do k=1,nz
  !   do j=1,ny
  !     do i=1,nx
  !       totSumI = totSumI + debugInt1(i,j,k)
  !     enddo
  !   enddo
  ! enddo
  ! write(*,fmt='(A, I10)') 'debugmkdel_d=', totSumI


  ! fnameFull = 'output/debugmkdel_' // trim(write_fmtnumb(step)) // '.vti'
  ! open(unit=iotest,file=trim(fnameFull),status='replace',action='write')

  ! extent =  trim(write_fmtnumb(1)) // ' ' // trim(write_fmtnumb(nx)) // ' ' &
  !       // trim(write_fmtnumb(1)) // ' ' // trim(write_fmtnumb(ny)) // ' ' &
  !       // trim(write_fmtnumb(1)) // ' ' // trim(write_fmtnumb(nz))

  ! write(iotest,*) '<VTKFile type="ImageData" version="1.0" byte_order="LittleEndian" >'
  ! write(iotest,*) ' <ImageData WholeExtent="' // trim(extent) // '" >'
  ! write(iotest,*) ' <Piece Extent="' // trim(extent) // '">'
  ! write(iotest,*) '   <PointData>'

  ! if (textual) then
  !   write(iotest,*) '    <DataArray type="UInt8" Name="debugmkdel" format="ascii" >'

  !   do k=1,nz
  !     do j=1,ny
  !       do i=1,nx
  !         write(iotest,fmt='("     ", I3)') debugInt1(i,j,k)
  !       enddo
  !     enddo
  !   enddo

  !   write(iotest,*) '    </DataArray>'
  !   write(iotest,*) '   </PointData>'
  !   write(iotest,*) ' </Piece>'
  !   write(iotest,*) ' </ImageData>'
  ! else
  !   write(iotest,*) '    <DataArray type="UInt8" Name="debugmkdel" format="appended" offset="0" />'
  !   write(iotest,*) '   </PointData>'
  !   write(iotest,*) ' </Piece>'
  !   write(iotest,*) ' </ImageData>'
  !   write(iotest,*) ' <AppendedData encoding="raw">'
  !   close(iotest)

  !   open(unit=iotest,file=trim(fnameFull),status='old',position='append', &
  !               access='stream',form='unformatted',action='write')
  !   write(iotest) '_'
  !   length = 1*nx*ny*nz
  !   write(iotest) length, debugInt1(1:nx,1:ny,1:nz)
  !   close(iotest)

  !   open(unit=iotest,file=trim(fnameFull),status='old',position='append', &
  !               action='write')
  !   write(iotest,*) ''
  !   write(iotest,*) ' </AppendedData>'
  ! endif

  ! write(iotest,*) '</VTKFile >'
  ! close(iotest)


  ! debugInt1 = debugbb_d

  ! fnameFull = 'output/debugbb_' // trim(write_fmtnumb(step)) // '.vti'
  ! open(unit=iotest,file=trim(fnameFull),status='replace',action='write')

  ! extent =  trim(write_fmtnumb(1)) // ' ' // trim(write_fmtnumb(nx)) // ' ' &
  !       // trim(write_fmtnumb(1)) // ' ' // trim(write_fmtnumb(ny)) // ' ' &
  !       // trim(write_fmtnumb(1)) // ' ' // trim(write_fmtnumb(nz))

  ! write(iotest,*) '<VTKFile type="ImageData" version="1.0" byte_order="LittleEndian" >'
  ! write(iotest,*) ' <ImageData WholeExtent="' // trim(extent) // '" >'
  ! write(iotest,*) ' <Piece Extent="' // trim(extent) // '">'
  ! write(iotest,*) '   <PointData>'

  ! if (textual) then
  !   write(iotest,*) '    <DataArray type="UInt8" Name="debugbb" format="ascii" >'

  !   do k=1,nz
  !     do j=1,ny
  !       do i=1,nx
  !         write(iotest,fmt='("     ", I3)') debugInt1(i,j,k)
  !       enddo
  !     enddo
  !   enddo

  !   write(iotest,*) '    </DataArray>'
  !   write(iotest,*) '   </PointData>'
  !   write(iotest,*) ' </Piece>'
  !   write(iotest,*) ' </ImageData>'
  ! else
  !   write(iotest,*) '    <DataArray type="UInt8" Name="debugbb" format="appended" offset="0" />'
  !   write(iotest,*) '   </PointData>'
  !   write(iotest,*) ' </Piece>'
  !   write(iotest,*) ' </ImageData>'
  !   write(iotest,*) ' <AppendedData encoding="raw">'
  !   close(iotest)

  !   open(unit=iotest,file=trim(fnameFull),status='old',position='append', &
  !               access='stream',form='unformatted',action='write')
  !   write(iotest) '_'
  !   length = 1*nx*ny*nz
  !   write(iotest) length, debugInt1(1:nx,1:ny,1:nz)
  !   close(iotest)

  !   open(unit=iotest,file=trim(fnameFull),status='old',position='append', &
  !               action='write')
  !   write(iotest,*) ''
  !   write(iotest,*) ' </AppendedData>'
  ! endif

  ! write(iotest,*) '</VTKFile >'
  ! close(iotest)
 end subroutine writeImageDataVTI_isfluid


 subroutine writeImageDataVTI_2d(fname, step, varname, plane, z, textual)
    use dimensions_m
    implicit none
    integer(4), intent(in) :: z
    character(len=*),intent(in) :: fname, varname
    integer,intent(in) :: step
    real(4)            :: plane(nx+2,ny+2)
    logical,intent(in) :: textual
    character(len=120) :: fnameFull,extent
    integer i,j,k, iotest
    integer(4)         :: length
  
  
    iotest = 55
    fnameFull = 'output/' // trim(fname) // '_' // trim(write_fmtnumb(myrank)) // '.' //trim(write_fmtnumb(step)) // '.vti'
    open(unit=iotest,file=trim(fnameFull),status='replace',action='write')
  
    extent =  trim(write_fmtnumb(0)) // ' ' // trim(write_fmtnumb(nx+1)) // ' ' &
          // trim(write_fmtnumb(0)) // ' ' // trim(write_fmtnumb(ny+1)) // ' ' &
          // trim(write_fmtnumb(offset(3) + z)) // ' ' // trim(write_fmtnumb(offset(3) + z))
  
    write(iotest,*) '<VTKFile type="ImageData" version="1.0" byte_order="LittleEndian" >'
    write(iotest,*) ' <ImageData WholeExtent="' // trim(extent) // '" >'
    write(iotest,*) ' <Piece Extent="' // trim(extent) // '">'
    write(iotest,*) '   <PointData>'
  
    if (textual) then
      write(iotest,*) '    <DataArray type="Float32" Name="',trim(varname),'" format="ascii" >'
  
      
      do j=0,ny+1
        do i=0,nx+1
          write(iotest,fmt='("     ", F20.8)') plane(i,j)
        enddo
      enddo
  
      write(iotest,*) '    </DataArray>'
      write(iotest,*) '   </PointData>'
      write(iotest,*) ' </Piece>'
      write(iotest,*) ' </ImageData>'
    else
      write(iotest,*) '    <DataArray type="Float32" Name="',trim(varname),'" format="appended" offset="0" />'
      write(iotest,*) '   </PointData>'
      write(iotest,*) ' </Piece>'
      write(iotest,*) ' </ImageData>'
      write(iotest,*) ' <AppendedData encoding="raw">'
      close(iotest)
  
      open(unit=iotest,file=trim(fnameFull),status='old',position='append', &
                  access='stream',form='unformatted',action='write')
      write(iotest) '_'
      length = 4*(nx+2)*(ny+2)
      write(iotest) length, plane(0:nx+1, 0:ny+1)
      close(iotest)
  
      open(unit=iotest,file=trim(fnameFull),status='old',position='append', &
                  action='write')
      write(iotest,*) ''
      write(iotest,*) ' </AppendedData>'
    endif
  
    write(iotest,*) '</VTKFile >'
    close(iotest)
   end subroutine writeImageDataVTI_2d

end module write_output

