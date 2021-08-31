 module write_output
    use dimensions_m
    implicit none

    ! integer(1)  :: debugInt1(0:nx+1,0:ny+1,0:nz+1)
    ! real        :: debugReal(0:nx+1,0:ny+1,0:nz+1)
 contains

 subroutine writeImageDataVTI(nz, fname, step, varname, rhoR, textual)
  use dimensions_m
  implicit none
  integer(4), intent(in) :: nz
  character(len=*),intent(in) :: fname, varname
  integer,intent(in) :: step
  real(4),allocatable,intent(in) :: rhoR(:,:,:)
  logical,intent(in) :: textual
  character(len=120) :: fnameFull,extent
  integer i,j,k, iotest
  integer(4)         :: length


  iotest = 55
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
          write(iotest,fmt='("     ", F20.8)') rhoR(i,j,k)
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
    length = 4*nx*ny*nz
    write(iotest) length, rhoR(1:nx,1:ny,1:nz)
    close(iotest)

    open(unit=iotest,file=trim(fnameFull),status='old',position='append', &
                action='write')
    write(iotest,*) ''
    write(iotest,*) ' </AppendedData>'
  endif

  write(iotest,*) '</VTKFile >'
  close(iotest)
 end subroutine writeImageDataVTI


 subroutine writeImageDataVTI_3d(nz, fname, step, vel, textual)
  use dimensions_m
  implicit none
  integer(4), intent(in) :: nz
  character(len=*),intent(in) :: fname
  integer,intent(in) :: step
  real(4),allocatable,intent(in) :: vel(:,:,:,:)
  logical,intent(in) :: textual
  character(len=120) :: fnameFull,extent
  integer i,j,k, iotest
  integer(4)         :: length


  iotest = 55
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
          write(iotest,fmt='("     ", F20.8)') vel(1,i,j,k),vel(2,i,j,k),vel(3,i,j,k)
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
    write(iotest) length, vel(1:3, 1:nx,1:ny,1:nz)
    close(iotest)

    open(unit=iotest,file=trim(fnameFull),status='old',position='append', &
                action='write')
    write(iotest,*) ''
    write(iotest,*) ' </AppendedData>'
  endif

  write(iotest,*) '</VTKFile >'
  close(iotest)
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

  ! write(iotest,fmt=100) ''
  ! write(iotest,fmt=100) 'VECTORS vel float'
  ! do iatm = 1 ,numAtoms
  !   write(iotest,*) v_atm(1, iatm), v_atm(2, iatm), v_atm(3, iatm)
  ! enddo

  if (lrotate) then
    write(iotest,fmt=100) 'VECTORS Vectx float'
    do iatm = 1 ,numAtoms
      write(iotest,201)take_rotversorx(q(1,iatm),q(2,iatm),q(3,iatm),q(4,iatm))
    enddo
    
    write(iotest,fmt=100) 'VECTORS Vecty float'
    do iatm = 1 ,numAtoms
      write(iotest,201)take_rotversory(q(1,iatm),q(2,iatm),q(3,iatm),q(4,iatm))
    enddo
    
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

 subroutine writeImageDataVTI_isfluid(nz, fname, step, myfluid, flip, textual)
  use dimensions_m
  use kernels_fluid
  implicit none
  integer(4), intent(in) :: nz
  character(len=*),intent(in) :: fname
  integer,intent(in) :: step, flip
  integer(1), allocatable, intent(in) :: myfluid(:,:,:,:)
  logical,intent(in) :: textual
  character(len=120) :: fnameFull,extent
  integer i,j,k, iotest
  integer(4)         :: length
  real(8)            :: totSumR
  integer(8)         :: totSumI

  iotest = 55

  ! fnameFull = 'output/' // 'partisfluid_' // trim(write_fmtnumb(step)) // '.txt'
  ! open(unit=iotest,file=trim(fnameFull),status='replace',action='write')

  totSumI = 0
  do k=1,nz
    do j=1,ny
      do i=1,nx
        if (myfluid(i,j,k, flip) /= fluid_fluid) then
          totSumI = totSumI + 1
          ! write(iotest,fmt='(3I5, I3)') i,j,k, myfluid(i,j,k, flip)
        endif
      enddo
    enddo
  enddo
  write(*,fmt='(A,I5, A,I10)') 'writeImageDataVTI] step=', step, ' myfluid=', totSumI

  ! close(iotest)
  
  
  fnameFull = 'output/' // trim(fname) // '_' // trim(write_fmtnumb(step)) // '.vti'
  open(unit=iotest,file=trim(fnameFull),status='replace',action='write')

  extent =  trim(write_fmtnumb(1)) // ' ' // trim(write_fmtnumb(nx)) // ' ' &
        // trim(write_fmtnumb(1)) // ' ' // trim(write_fmtnumb(ny)) // ' ' &
        // trim(write_fmtnumb(1)) // ' ' // trim(write_fmtnumb(nz))

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

end module write_output

