#define timer(func, store) call system_clock(start_t, rate); call func; call system_clock(end_t); store  = store + real(end_t - start_t)/real(rate);
program test_smoothing
  use util_mod
  implicit none
  integer :: nelem, npoin
  integer, parameter :: loops = 1
  integer, allocatable, dimension(:,:) :: inpoel
  type(Vector), allocatable, dimension(:) :: coord
  logical, allocatable, dimension(:) :: fixed
  integer :: ipoin, ielem, IOstatus, dummy
  character(len=32) :: arg1!, arg2
  integer :: start_t, end_t, rate
  integer :: i
  real*4 :: time
  integer, dimension(:), allocatable :: points
  type(Body) :: body0
  type(IndexedVector), dimension(:), allocatable :: ds
  type(Triangle) :: a
  !%%%%%%%%%%%%%%%%%%%%
  time = 0.d0
  call getarg(iargc(), arg1)
  open(1, file =trim(arg1), status = 'old')
  read(1,*)
  read(1,*) nelem, npoin
  allocate(inpoel(3,nelem))
  allocate(coord(npoin))
  allocate(fixed(npoin))
  allocate(points(npoin))
  allocate(ds(npoin))
  do ipoin = 1, npoin
     !SE MUEVEN TODOS POR DEFECTO
     fixed(ipoin) = .false.
  end do
  !LECTURA DE LAS COORDENADAS
  read(1,*)
  do ipoin = 1, npoin
     read(1,*) dummy, coord(ipoin)%x, coord(ipoin)%y
  end do
  !LECTURA DE INPOEL
  read(1,*)
  do ielem = 1, nelem
     read(1,*) dummy, inpoel(1,ielem), inpoel(2,ielem), inpoel(3,ielem)
  end do
  !LECTURA DE BORDES (FIXED)
  read(1,*)
  do 
     read(1,*, IOSTAT = IOstatus) ipoin, dummy
     if(IOstatus /= 0) exit
     fixed(ipoin) = .true.
  end do
  !%%%%%%%%%%%%%%%%%%%%
  do ipoin = 1, npoin
     points(ipoin) = ipoin
  end do

  call body0%init(points)
  
  do ipoin = 1, npoin
     coord(ipoin)%x = coord(ipoin)%x + ds(ipoin)%x
     coord(ipoin)%y = coord(ipoin)%y + ds(ipoin)%y
  end do
     
  ! print*, 'Time:', time/loops
  call IOprint(npoin, nelem, coord(:)%x, coord(:)%y, inpoel, "remeshing.msh", 0.d0)

end program test_smoothing

subroutine IOprint(npoin, nelem, X, Y, inpoel, filename, setoff)
  integer, intent(in) :: npoin, nelem, inpoel(3,nelem)
  real*8, intent(in) :: X(npoin), Y(npoin), setoff
  character(len=*), intent(in) :: filename
  !%%%%%%%%%%%%%%%%%%%%
  integer :: ipoin, ielem
  !%%%%%%%%%%%%%%%%%%%%
  open(2, file=trim(filename), status='unknown')
  write(2, '(A)') 'BackgroundMesh V 1.0'
  write(2, '(A)') 'MESH dimension 2 ElemType Triangle Nnode 3'
  write(2, '(A)') 'Coordinates'
  do ipoin = 1, npoin
     write(2, '(I7, 3E14.6)') ipoin, X(ipoin) + setoff, Y(ipoin)
  end do
  write(2, '(A)') 'end Coordinates'

  write(2, '(A)') 'Elements'
  do ielem = 1, nelem
     write(2, '(5I7)') ielem, inpoel(1, ielem), inpoel(2, ielem), inpoel(3, ielem)
  end do
  write(2, '(A)') 'end Elements'
  write(2, '(A)') 'DesiredSize Nodes'

  ! do ipoin = 1, npoin
  !    write(2, '(I7, E14.6)') ipoin, HH_NEW(ipoin)
  ! end do
  write(2, '(A)') 'end DesiredSize Nodes'
  write(*, *)
  write(*, *)'****-------> FIN DEL CALCULO <-------****'
  write(*, *)
endsubroutine IOprint
