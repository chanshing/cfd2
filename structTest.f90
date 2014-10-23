program structTest
  use, intrinsic :: iso_fortran_env  
  !$ use omp_lib
  
  implicit none
!  integer, parameter :: cp = REAL32
  integer, parameter :: cp = REAL64
!  integer, parameter :: cp = REAL128

  !===========================
  ! STRUCTURE  
  type :: structQuad
    real(cp)  :: coord1(3)
    real(cp)  :: coord2(3)
    real(cp)  :: coord3(3)
    real(cp)  :: coord4(3)
    real(cp)  :: n(3)
    real(cp)  :: area
    !---
    integer   :: node_numbers(4)
    integer   :: id
    character :: type
    logical   :: is_triangle
    complex(cp)  :: some_large_array(100)
  end type

  type(structQuad),allocatable  :: q(:)
  !===========================
  ! STATIC ARRAYS 
  real(cp),allocatable  :: q_coord1(:,:)
  real(cp),allocatable  :: q_coord2(:,:)
  real(cp),allocatable  :: q_coord3(:,:)
  real(cp),allocatable  :: q_coord4(:,:)
  real(cp),allocatable  :: q_n(:,:)
  real(cp),allocatable  :: q_area(:)
  !---
  integer,allocatable   :: q_node_numbers(:,:)
  integer,allocatable   :: q_id(:)
  character,allocatable :: q_type(:)
  logical,allocatable   :: q_is_triangle(:)
  complex(cp),allocatable  :: q_some_large_array(:,:)
  !===========================
  ! STRUCTURE: HOT AND COLD
  type :: structQuadHot
    real(cp)  :: coord1(3)
    real(cp)  :: coord2(3)
    real(cp)  :: coord3(3)
    real(cp)  :: coord4(3)
    real(cp)  :: n(3)
    real(cp)  :: area
  end type
  !---
  type :: structQuadCold
    integer   :: node_numbers(4)
    integer   :: id
    character :: type
    logical   :: is_triangle
    complex(cp)  :: some_large_array(100)
  end type

  type(structQuadHot),allocatable,target  :: qHot(:)
  type(structQuadHot), pointer     :: qHotPtr
  type(structQuadCold),allocatable,target  :: qCold(:)
  type(structQuadCold), pointer     :: qColdPtr
  !===========================

!$ real(REAL64) :: t1, t2, t3
  integer           :: stat
  integer,parameter :: Nx=600, Ny=600
  integer,parameter :: NOQUAD=Nx*Ny
  real(cp),parameter:: X=1.e0_cp, Y=1.e0_cp
  real(cp),parameter:: dx=X/Nx, dy=Y/Ny
  real(cp)          :: a(3), b(3)
  ! Loop counters:
  integer           :: i, ii, c
  
  !===
  
  !===
  ! 1 Static arrays
  !===
  !$ t1 = omp_get_wtime()
  allocate( q_coord1(3,NOQUAD), &
            q_coord2(3,NOQUAD), &
            q_coord3(3,NOQUAD), &
            q_coord4(3,NOQUAD), &
            q_n(3,NOQUAD), &
            q_area(NOQUAD), &
            q_node_numbers(4,NOQUAD), &
            q_id(NOQUAD), &
            q_type(NOQUAD), &
            q_is_triangle(NOQUAD), &
            q_some_large_array(100,NOQUAD), &
            stat=stat ) 
  if ( stat /= 0 ) stop 'Cannot allocate memory!'
  
  ! Fill cold part
  !$omp parallel workshare
  q_id = 0
  q_type = 'q'
  q_is_triangle = .false.
  q_node_numbers = 0
  q_some_large_array = ( 0.e0_cp, 1.e0_cp )
  !$omp end parallel workshare
  
  ! Fill hot part
  !$omp parallel do private(i,ii,c)
  do i=1,Nx
    do ii=1,Ny
      c = ii+(i-1)*Ny
      ! Coordinates
      q_coord1(:,c) = [ real(i-1,kind=cp)*dx, real(ii-1,kind=cp)*dy, 0.e0_cp ]
      q_coord2(:,c) = [ real(i  ,kind=cp)*dx, real(ii-1,kind=cp)*dy, 0.e0_cp ]
      q_coord3(:,c) = [ real(i  ,kind=cp)*dx, real(ii  ,kind=cp)*dy, 0.e0_cp ]
      q_coord4(:,c) = [ real(i-1,kind=cp)*dx, real(ii  ,kind=cp)*dy, 0.e0_cp ]
    enddo ! ii
  enddo ! i
  !$omp end parallel do
  
  !$ t2 = omp_get_wtime()
  
  ! Calculate
  !$omp parallel do private(c, a, b)
  do c=1,NOQUAD  
    ! Normal vector
    a = q_coord2(:,c) - q_coord1(:,c)
    b = q_coord3(:,c) - q_coord1(:,c)
    q_n(1,c) = a(2) * b(3) - a(3) * b(2)
    q_n(2,c) = a(3) * b(1) - a(1) * b(3)
    q_n(3,c) = a(1) * b(2) - a(2) * b(1)
    ! Normalize
    q_area(c) = sqrt(dot_product(q_n(:,c), q_n(:,c)))
    q_n(:,c) = q_n(:,c) / q_area(c)
  enddo ! c
  !$omp end parallel do
  
  !$ t3 = omp_get_wtime()
  
  !$ write(*,'(A,F10.3,A,F10.3)') 'STATIC ARRAYS:     Setup ',t2-t1, ' Calculations ',t3-t1
  if ( allocated(q_coord1) ) deallocate (q_coord1)
  if ( allocated(q_coord2) ) deallocate (q_coord2)
  if ( allocated(q_coord3) ) deallocate (q_coord3)
  if ( allocated(q_coord4) ) deallocate (q_coord4)
  if ( allocated(q_n) ) deallocate(q_n)
  if ( allocated(q_area) ) deallocate(q_area)
  if ( allocated(q_node_numbers) ) deallocate(q_node_numbers)
  if ( allocated(q_id) ) deallocate(q_id)
  if ( allocated(q_type) ) deallocate(q_type)
  if ( allocated(q_is_triangle) ) deallocate(q_is_triangle)
  if ( allocated(q_some_large_array) ) deallocate(q_some_large_array)

  !===
  ! 2 Structure
  !===
  !$ t1 = omp_get_wtime()
  allocate( q(NOQUAD), stat=stat ) 
  if ( stat /= 0 ) stop 'Cannot allocate memory!'
  
  ! Fill cold part
  !$omp parallel do private(c)
  do c=1,NOQUAD
    q(c)%id = 0
    q(c)%type = 'q'
    q(c)%is_triangle = .false.
    q(c)%node_numbers = 0
    q(c)%some_large_array = ( 0.e0_cp, 1.e0_cp )
  enddo
  !$omp end parallel do

  ! Fill hot part
  !$omp parallel do private(i,ii,c)
  do i=1,Nx
    do ii=1,Ny
      c = ii+(i-1)*Ny
      ! Coordinates
      q(c)%coord1 = [ real(i-1,kind=cp)*dx, real(ii-1,kind=cp)*dy, 0.e0_cp ]
      q(c)%coord2 = [ real(i  ,kind=cp)*dx, real(ii-1,kind=cp)*dy, 0.e0_cp ]
      q(c)%coord3 = [ real(i  ,kind=cp)*dx, real(ii  ,kind=cp)*dy, 0.e0_cp ]
      q(c)%coord4 = [ real(i-1,kind=cp)*dx, real(ii  ,kind=cp)*dy, 0.e0_cp ]
    enddo ! ii
  enddo ! i
  !$omp end parallel do
  !$ t2 = omp_get_wtime()
  ! Calculate
  !$omp parallel do private(c, a, b)
  do c=1,NOQUAD  
    ! Normal vector
    a = q(c)%coord2 - q(c)%coord1
    b = q(c)%coord3 - q(c)%coord1
    q(c)%n(1) = a(2) * b(3) - a(3) * b(2)
    q(c)%n(2) = a(3) * b(1) - a(1) * b(3)
    q(c)%n(3) = a(1) * b(2) - a(2) * b(1)
    ! Normalize
    q(c)%area = sqrt(dot_product(q(c)%n, q(c)%n))
    q(c)%n = q(c)%n / q(c)%area
  enddo ! c
  !$omp end parallel do
  !$ t3 = omp_get_wtime()
  !$ write(*,'(A,F10.3,A,F10.3)') 'STRUCTURE:         Setup ',t2-t1, ' Calculations ',t3-t1
   
  !===
  ! 3 Structure (HOT & COLD)
  !===
  !$ t1 = omp_get_wtime()
  allocate( qHot(NOQUAD), stat=stat ) 
  allocate( qCold(NOQUAD), stat=stat ) 
  if ( stat /= 0 ) stop 'Cannot allocate memory!'
  
  ! Fill cold part
  !$omp parallel do private(c)
  do c=1,NOQUAD
    qCold(c)%id = 0
    qCold(c)%type = 'q'
    qCold(c)%is_triangle = .false.
    qCold(c)%node_numbers = 0
    qCold(c)%some_large_array = ( 0.e0_cp, 1.e0_cp )
  enddo
  !$omp end parallel do
  ! Fill hot part
  !$omp parallel do private(c)
  do i=1,Nx
    do ii=1,Ny
      c = ii+(i-1)*Ny
      ! Coordinates
      qHot(c)%coord1 = [ real(i-1,kind=cp)*dx, real(ii-1,kind=cp)*dy, 0.e0_cp ]
      qHot(c)%coord2 = [ real(i  ,kind=cp)*dx, real(ii-1,kind=cp)*dy, 0.e0_cp ]
      qHot(c)%coord3 = [ real(i  ,kind=cp)*dx, real(ii  ,kind=cp)*dy, 0.e0_cp ]
      qHot(c)%coord4 = [ real(i-1,kind=cp)*dx, real(ii  ,kind=cp)*dy, 0.e0_cp ]
    enddo ! ii
  enddo ! i
  !$omp end parallel do
  !$ t2 = omp_get_wtime()
  ! Calculate
  !$omp parallel do private(c,a,b)
  do c=1,NOQUAD  
    ! Normal vector
    a = qHot(c)%coord2 - qHot(c)%coord1
    b = qHot(c)%coord3 - qHot(c)%coord1
    qHot(c)%n(1) = a(2) * b(3) - a(3) * b(2)
    qHot(c)%n(2) = a(3) * b(1) - a(1) * b(3)
    qHot(c)%n(3) = a(1) * b(2) - a(2) * b(1)
    ! Normalize
    qHot(c)%area = sqrt(dot_product(qHot(c)%n, qHot(c)%n))
    qHot(c)%n = qHot(c)%n / qHot(c)%area
  enddo ! c
  !$omp end parallel do
  !$ t3 = omp_get_wtime()
  !$ write(*,'(A,F10.3,A,F10.3)') 'STRUCTURE HOT    : Setup ',t2-t1, ' Calculations ',t3-t1
  if ( allocated(qHot) ) deallocate(qHot)
  if ( allocated(qCold) ) deallocate(qCold)

  !===
  ! 4 Structure (HOT & COLD) + PTR
  !===
  !$ t1 = omp_get_wtime()
  allocate( qHot(NOQUAD), stat=stat ) 
  allocate( qCold(NOQUAD), stat=stat ) 
  if ( stat /= 0 ) stop 'Cannot allocate memory!'
  
  ! Fill cold part
  !$omp parallel do private(c,qColdPtr)
  do c=1,NOQUAD
    qColdPtr => qCold(c)
    qColdPtr%id = 0
    qColdPtr%type = 'q'
    qColdPtr%is_triangle = .false.
    qColdPtr%node_numbers = 0
    qColdPtr%some_large_array = ( 0.e0_cp, 1.e0_cp )
  enddo !c
  !$omp end parallel do
  ! Fill hot part
!  c = 0 ! Init counter
  !$omp parallel do private(i,ii,c,qHotPtr)
  do i=1,Nx
    do ii=1,Ny
!      c = c + 1 ! Increment counter
      c = ii+(i-1)*Ny
      qHotPtr => qHot(c)
      ! Coordinates
      qHotPtr%coord1 = [ real(i-1,kind=cp)*dx, real(ii-1,kind=cp)*dy, 0.e0_cp ]
      qHotPtr%coord2 = [ real(i  ,kind=cp)*dx, real(ii-1,kind=cp)*dy, 0.e0_cp ]
      qHotPtr%coord3 = [ real(i  ,kind=cp)*dx, real(ii  ,kind=cp)*dy, 0.e0_cp ]
      qHotPtr%coord4 = [ real(i-1,kind=cp)*dx, real(ii  ,kind=cp)*dy, 0.e0_cp ]
    enddo ! ii
  enddo ! i
  !$omp end parallel do
  !$ t2 = omp_get_wtime()
  ! Calculate
  !$omp parallel do private(c, a, b, qHotPtr)
  do c=1,NOQUAD  
    qHotPtr => qHot(c)
    ! Normal vector
    a = qHotPtr%coord2 - qHotPtr%coord1
    b = qHot(c)%coord3 - qHotPtr%coord1
    qHotPtr%n(1) = a(2) * b(3) - a(3) * b(2)
    qHotPtr%n(2) = a(3) * b(1) - a(1) * b(3)
    qHotPtr%n(3) = a(1) * b(2) - a(2) * b(1)
    ! Normalize
    qHotPtr%area = sqrt(dot_product(qHotPtr%n, qHotPtr%n))
    qHotPtr%n = qHotPtr%n / qHotPtr%area
  enddo ! c
  !$omp end parallel do
  !$ t3 = omp_get_wtime()
  !$ write(*,'(A,F10.3,A,F10.3)') 'STRUCTURE HOT,PTR: Setup ',t2-t1, ' Calculations ',t3-t1
  if ( allocated(qHot) ) deallocate(qHot)
  if ( allocated(qCold) ) deallocate(qCold)

end program
