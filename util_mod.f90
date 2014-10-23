module util_mod
  implicit none
  private :: &
       wallType_addEdge,&
       wallType_updateNormals

  type IndexedScalar
     integer :: point
     real(kind=8) :: val
  end type IndexedScalar

  type Vector
     real(kind=8) :: x
     real(kind=8) :: y
  end type Vector
  
  type, extends(Vector) :: IndexedVector
     integer :: point
  end type IndexedVector

  type, extends(IndexedVector) :: ValuedIndexedVector
     real(kind=8) :: val
  end type ValuedIndexedVector

  type Edge
     integer, dimension(2) :: points
  end type Edge

  type, extends(Edge) :: IndexedEdge
     integer :: element
  end type IndexedEdge

  type Triangle
     integer, dimension(3) :: points
     real(kind=8) :: area
  end type Triangle

  type WallType
     integer, dimension(:), allocatable :: points
     type(IndexedEdge), dimension(:), allocatable :: edges
     type(ValuedIndexedVector), dimension(:), allocatable :: normals
   contains
     procedure, public :: addEdge => wallType_addEdge
     procedure, public :: updateNormals => wallType_updateNormals
  end type WallType

  type CSR
     real(kind=8), dimension(:), allocatable :: val
     real(kind=8), dimension(:), allocatable :: col_idx
     real(kind=8), dimension(:), allocatable :: row_ptr
  end type CSR

  type Link
     real(kind=8), dimension(:), allocatable :: val
     real(kind=8), dimension(:), allocatable :: ptr
  end type Link

  interface assignment(=)
     procedure assign_real_to_vector
     procedure assign_real_to_indexedVector
     procedure assign_real_to_valuedIndexedVector
     procedure assign_vector_to_indexedVector
     procedure assign_vector_to_valuedIndexedVector
     procedure assign_indexedVector_to_vector
     procedure assign_indexedVector_to_valuedIndexedvector
     procedure assign_valuedIndexedVector_to_vector
     procedure assign_valuedIndexedVector_to_indexedvector
     procedure assign_indexedEdge_to_edge
     procedure assign_edge_to_indexedEdge
  end interface assignment(=)

  interface operator(+)
     procedure add_real_vector
     procedure add_real_indexedVector
     procedure add_real_valuedIndexedVector
     procedure add_vector_vector
     procedure add_vector_real
     procedure add_vector_indexedVector
     procedure add_vector_valuedIndexedVector
     procedure add_indexedVector_indexedVector
     procedure add_indexedVector_real
     procedure add_indexedVector_vector
     procedure add_indexedVector_valuedIndexedVector
     procedure add_valuedIndexedVector_valuedIndexedVector
     procedure add_valuedIndexedVector_real
     procedure add_valuedIndexedVector_vector
     procedure add_valuedIndexedVector_indexedVector
  end interface operator(+)

  interface operator(-)
     procedure substract_real_vector
     procedure substract_real_indexedVector
     procedure substract_real_valuedIndexedVector
     procedure substract_vector_vector
     procedure substract_vector_real
     procedure substract_vector_indexedVector
     procedure substract_vector_valuedIndexedVector
     procedure substract_indexedVector_indexedVector
     procedure substract_indexedVector_real
     procedure substract_indexedVector_vector
     procedure substract_indexedVector_valuedIndexedVector
     procedure substract_valuedIndexedVector_valuedIndexedVector
     procedure substract_valuedIndexedVector_real
     procedure substract_valuedIndexedVector_vector
     procedure substract_valuedIndexedVector_indexedVector
     procedure substract_vectorArray_vectorArray
  end interface operator(-)

  interface operator(*)
     procedure dot_real_vector
     procedure dot_real_indexedVector
     procedure dot_real_valuedIndexedVector
     procedure dot_vector_real
     procedure dot_vector_vector
     procedure dot_vector_indexedVector
     procedure dot_vector_valuedIndexedVector
     procedure dot_indexedVector_real
     procedure dot_indexedVector_vector
     procedure dot_indexedVector_indexedVector
     procedure dot_indexedVector_valuedIndexedVector
     procedure dot_valuedIndexedVector_real
     procedure dot_valuedIndexedVector_vector
     procedure dot_valuedIndexedVector_indexedVector
     procedure dot_valuedIndexedVector_valuedIndexedVector
  end interface operator(*)

  interface operator(/)
     procedure div_vector_real
     procedure div_indexedVector_real
     procedure div_valuedIndexedVector_real
  end interface operator(/)

contains

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%  ASSIGNMENT OVERLOADING %%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  pure subroutine assign_real_to_vector(v1,a)
    type(Vector), intent(out) :: v1
    real(kind=8), intent(in) :: a
    v1%x = a
    v1%y = a
  end subroutine assign_real_to_vector

  pure subroutine assign_real_to_indexedVector(v1,a)
    type(IndexedVector), intent(inout) :: v1
    real(kind=8), intent(in) :: a
    v1%x = a
    v1%y = a
  end subroutine assign_real_to_indexedVector

  pure subroutine assign_real_to_valuedIndexedVector(v1,a)
    type(ValuedIndexedVector), intent(inout) :: v1
    real(kind=8), intent(in) :: a
    v1%x = a
    v1%y = a
  end subroutine assign_real_to_valuedIndexedVector

  pure subroutine assign_vector_to_indexedVector(v1,v2)
    type(IndexedVector), intent(inout) :: v1
    type(Vector), intent(in) :: v2
    v1%x = v2%x
    v1%y = v2%y
  end subroutine assign_vector_to_indexedVector

  pure subroutine assign_vector_to_valuedIndexedVector(v1,v2)
    type(ValuedIndexedVector), intent(inout) :: v1
    type(Vector), intent(in) :: v2
    v1%x = v2%x
    v1%y = v2%y
  end subroutine assign_vector_to_valuedIndexedVector

  pure subroutine assign_indexedVector_to_vector(v1,v2)
    type(Vector), intent(out) :: v1
    type(IndexedVector), intent(in) :: v2
    v1%x = v2%x
    v1%y = v2%y
  end subroutine assign_indexedVector_to_vector

  pure subroutine assign_indexedVector_to_valuedIndexedvector(v1,v2)
    type(ValuedIndexedVector), intent(inout) :: v1
    type(IndexedVector), intent(in) :: v2
    v1%x = v2%x
    v1%y = v2%y
  end subroutine assign_indexedVector_to_valuedIndexedvector

  pure subroutine assign_valuedIndexedVector_to_vector(v1,v2)
    type(Vector), intent(out) :: v1
    type(ValuedIndexedVector), intent(in) :: v2
    v1%x = v2%x
    v1%y = v2%y
  end subroutine assign_valuedIndexedVector_to_vector

  pure subroutine assign_valuedIndexedVector_to_indexedVector(v1,v2)
    type(IndexedVector), intent(inout) :: v1
    type(ValuedIndexedVector), intent(in) :: v2
    v1%x = v2%x
    v1%y = v2%y
  end subroutine assign_valuedIndexedVector_to_indexedVector

  pure subroutine assign_indexedEdge_to_edge(e1,e2)
    type(Edge), intent(out) :: e1
    type(IndexedEdge), intent(in) :: e2
    e1%points(:) = e2%points(:)
  end subroutine assign_indexedEdge_to_edge

  pure subroutine assign_edge_to_indexedEdge(e1,e2)
    type(IndexedEdge), intent(inout) :: e1
    type(Edge), intent(in) :: e2
    e1%points(:) = e2%points(:)
  end subroutine assign_edge_to_indexedEdge

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%% OPERATOR(+) %%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  pure type(Vector) function add_real_vector(a,v2)
    real(kind=8), intent(in) :: a
    type(Vector), intent(in) :: v2
    add_real_vector%x = a + v2%x
    add_real_vector%y = a + v2%y
  end function add_real_vector

  pure type(IndexedVector) function add_real_indexedVector(a,v2)
    real(kind=8), intent(in) :: a
    type(IndexedVector), intent(in) :: v2
    add_real_indexedVector%point = v2%point
    add_real_indexedVector%x = a + v2%x
    add_real_indexedVector%y = a + v2%y
  end function add_real_indexedVector

  pure type(ValuedIndexedVector) function add_real_valuedIndexedVector(a,v2)
    real(kind=8), intent(in) :: a
    type(ValuedIndexedVector), intent(in) :: v2
    add_real_valuedIndexedVector%point = v2%point
    add_real_valuedIndexedVector%val = v2%val
    add_real_valuedIndexedVector%x = a + v2%x
    add_real_valuedIndexedVector%y = a + v2%y
  end function add_real_valuedIndexedVector
  
  pure type(Vector) function add_vector_vector(v1,v2)
    type(Vector), intent(in) :: v1, v2
    add_vector_vector%x = v1%x + v2%x
    add_vector_vector%y = v1%y + v2%y
  end function add_vector_vector

  pure type(Vector) function add_vector_real(v1,a)
    type(Vector), intent(in) :: v1
    real(kind=8), intent(in) :: a
    add_vector_real%x = v1%x + a
    add_vector_real%y = v1%y + a
  end function add_vector_real

  pure type(IndexedVector) function add_vector_indexedVector(v1,v2)
    type(Vector), intent(in) :: v1
    type(IndexedVector), intent(in) :: v2
    add_vector_indexedVector%point = v2%point
    add_vector_indexedVector%x = v1%x + v2%x
    add_vector_indexedVector%y = v1%y + v2%y
  end function add_vector_indexedVector

  pure type(ValuedIndexedVector) function add_vector_valuedIndexedVector(v1,v2)
    type(Vector), intent(in) :: v1
    type(ValuedIndexedVector), intent(in) :: v2
    add_vector_valuedIndexedVector%point = v2%point
    add_vector_valuedIndexedVector%val = v2%val
    add_vector_valuedIndexedVector%x = v1%x + v2%x
    add_vector_valuedIndexedVector%y = v1%y + v2%y
  end function add_vector_valuedIndexedVector

  pure type(IndexedVector) function add_indexedVector_indexedVector(v1,v2)
    type(IndexedVector), intent(in) :: v1, v2
    add_indexedVector_indexedVector%point = v1%point
    add_indexedVector_indexedVector%x = v1%x + v2%x
    add_indexedVector_indexedVector%y = v1%y + v2%y
  end function add_indexedVector_indexedVector

  pure type(IndexedVector) function add_indexedVector_real(v1,a)
    type(IndexedVector), intent(in) :: v1
    real(kind=8), intent(in) :: a
    add_indexedVector_real%point = v1%point
    add_indexedVector_real%x = v1%x + a
    add_indexedVector_real%y = v1%y + a
  end function add_indexedVector_real

  pure type(IndexedVector) function add_indexedVector_vector(v1,v2)
    type(IndexedVector), intent(in) :: v1
    type(Vector), intent(in) :: v2
    add_indexedVector_vector%point = v1%point
    add_indexedVector_vector%x = v1%x + v2%x
    add_indexedVector_vector%y = v1%y + v2%y
  end function add_indexedVector_vector

  pure type(ValuedIndexedVector) function add_indexedVector_valuedIndexedVector(v1,v2)
    type(IndexedVector), intent(in) :: v1
    type(ValuedIndexedVector), intent(in) :: v2
    add_indexedVector_valuedIndexedVector%point = v2%point
    add_indexedVector_valuedIndexedVector%val = v2%val
    add_indexedVector_valuedIndexedVector%x = v1%x + v2%x
    add_indexedVector_valuedIndexedVector%y = v1%y + v2%y
  end function add_indexedVector_valuedIndexedVector

  pure type(ValuedIndexedVector) function add_valuedIndexedVector_valuedIndexedVector(v1,v2)
    type(ValuedIndexedVector), intent(in) :: v1
    type(ValuedIndexedVector), intent(in) :: v2
    add_valuedIndexedVector_valuedIndexedVector%point = v1%point
    add_valuedIndexedVector_valuedIndexedVector%val = v1%val
    add_valuedIndexedVector_valuedIndexedVector%x = v1%x + v2%x
    add_valuedIndexedVector_valuedIndexedVector%y = v1%y + v2%y
  end function add_valuedIndexedVector_valuedIndexedVector

  pure type(ValuedIndexedVector) function add_valuedIndexedVector_real(v1,a)
    type(ValuedIndexedVector), intent(in) :: v1
    real(kind=8), intent(in) :: a
    add_valuedIndexedVector_real%point = v1%point
    add_valuedIndexedVector_real%val = v1%val
    add_valuedIndexedVector_real%x = v1%x + a
    add_valuedIndexedVector_real%y = v1%y + a
  end function add_valuedIndexedVector_real

  pure type(ValuedIndexedVector) function add_valuedIndexedVector_vector(v1,v2)
    type(ValuedIndexedVector), intent(in) :: v1
    type(Vector), intent(in) :: v2
    add_valuedIndexedVector_vector%point = v1%point
    add_valuedIndexedVector_vector%val = v1%val
    add_valuedIndexedVector_vector%x = v1%x + v2%x
    add_valuedIndexedVector_vector%y = v1%y + v2%y
  end function add_valuedIndexedVector_vector

  pure type(ValuedIndexedVector) function add_valuedIndexedVector_indexedVector(v1,v2)
    type(ValuedIndexedVector), intent(in) :: v1
    type(IndexedVector), intent(in) :: v2
    add_valuedIndexedVector_indexedVector%point = v1%point
    add_valuedIndexedVector_indexedVector%val = v1%val
    add_valuedIndexedVector_indexedVector%x = v1%x + v2%x
    add_valuedIndexedVector_indexedVector%y = v1%y + v2%y
  end function add_valuedIndexedVector_indexedVector

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%% OPERATOR(-) %%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  pure type(Vector) function substract_real_vector(a,v2)
    real(kind=8), intent(in) :: a
    type(Vector), intent(in) :: v2
    substract_real_vector%x = a - v2%x
    substract_real_vector%y = a - v2%y
  end function substract_real_vector

  pure type(IndexedVector) function substract_real_indexedVector(a,v2)
    real(kind=8), intent(in) :: a
    type(IndexedVector), intent(in) :: v2
    substract_real_indexedVector%point = v2%point
    substract_real_indexedVector%x = a - v2%x
    substract_real_indexedVector%y = a - v2%y
  end function substract_real_indexedVector

  pure type(ValuedIndexedVector) function substract_real_valuedIndexedVector(a,v2)
    real(kind=8), intent(in) :: a
    type(ValuedIndexedVector), intent(in) :: v2
    substract_real_valuedIndexedVector%point = v2%point
    substract_real_valuedIndexedVector%val = v2%val
    substract_real_valuedIndexedVector%x = a - v2%x
    substract_real_valuedIndexedVector%y = a - v2%y
  end function substract_real_valuedIndexedVector
  
  pure type(Vector) function substract_vector_vector(v1,v2)
    type(Vector), intent(in) :: v1, v2
    substract_vector_vector%x = v1%x - v2%x
    substract_vector_vector%y = v1%y - v2%y
  end function substract_vector_vector

  pure type(Vector) function substract_vector_real(v1,a)
    type(Vector), intent(in) :: v1
    real(kind=8), intent(in) :: a
    substract_vector_real%x = v1%x - a
    substract_vector_real%y = v1%y - a
  end function substract_vector_real

  pure type(IndexedVector) function substract_vector_indexedVector(v1,v2)
    type(Vector), intent(in) :: v1
    type(IndexedVector), intent(in) :: v2
    substract_vector_indexedVector%point = v2%point
    substract_vector_indexedVector%x = v1%x - v2%x
    substract_vector_indexedVector%y = v1%y - v2%y
  end function substract_vector_indexedVector

  pure type(ValuedIndexedVector) function substract_vector_valuedIndexedVector(v1,v2)
    type(Vector), intent(in) :: v1
    type(ValuedIndexedVector), intent(in) :: v2
    substract_vector_valuedIndexedVector%point = v2%point
    substract_vector_valuedIndexedVector%val = v2%val
    substract_vector_valuedIndexedVector%x = v1%x - v2%x
    substract_vector_valuedIndexedVector%y = v1%y - v2%y
  end function substract_vector_valuedIndexedVector

  pure type(IndexedVector) function substract_indexedVector_indexedVector(v1,v2)
    type(IndexedVector), intent(in) :: v1, v2
    substract_indexedVector_indexedVector%point = v1%point
    substract_indexedVector_indexedVector%x = v1%x - v2%x
    substract_indexedVector_indexedVector%y = v1%y - v2%y
  end function substract_indexedVector_indexedVector

  pure type(IndexedVector) function substract_indexedVector_real(v1,a)
    type(IndexedVector), intent(in) :: v1
    real(kind=8), intent(in) :: a
    substract_indexedVector_real%point = v1%point
    substract_indexedVector_real%x = v1%x - a
    substract_indexedVector_real%y = v1%y - a
  end function substract_indexedVector_real

  pure type(IndexedVector) function substract_indexedVector_vector(v1,v2)
    type(IndexedVector), intent(in) :: v1
    type(Vector), intent(in) :: v2
    substract_indexedVector_vector%point = v1%point
    substract_indexedVector_vector%x = v1%x - v2%x
    substract_indexedVector_vector%y = v1%y - v2%y
  end function substract_indexedVector_vector

  pure type(ValuedIndexedVector) function substract_indexedVector_valuedIndexedVector(v1,v2)
    type(IndexedVector), intent(in) :: v1
    type(ValuedIndexedVector), intent(in) :: v2
    substract_indexedVector_valuedIndexedVector%point = v2%point
    substract_indexedVector_valuedIndexedVector%val = v2%val
    substract_indexedVector_valuedIndexedVector%x = v1%x - v2%x
    substract_indexedVector_valuedIndexedVector%y = v1%y - v2%y
  end function substract_indexedVector_valuedIndexedVector

  pure type(ValuedIndexedVector) function substract_valuedIndexedVector_valuedIndexedVector(v1,v2)
    type(ValuedIndexedVector), intent(in) :: v1
    type(ValuedIndexedVector), intent(in) :: v2
    substract_valuedIndexedVector_valuedIndexedVector%point = v1%point
    substract_valuedIndexedVector_valuedIndexedVector%val = v1%val
    substract_valuedIndexedVector_valuedIndexedVector%x = v1%x - v2%x
    substract_valuedIndexedVector_valuedIndexedVector%y = v1%y - v2%y
  end function substract_valuedIndexedVector_valuedIndexedVector

  pure type(ValuedIndexedVector) function substract_valuedIndexedVector_real(v1,a)
    type(ValuedIndexedVector), intent(in) :: v1
    real(kind=8), intent(in) :: a
    substract_valuedIndexedVector_real%point = v1%point
    substract_valuedIndexedVector_real%val = v1%val
    substract_valuedIndexedVector_real%x = v1%x - a
    substract_valuedIndexedVector_real%y = v1%y - a
  end function substract_valuedIndexedVector_real

  pure type(ValuedIndexedVector) function substract_valuedIndexedVector_vector(v1,v2)
    type(ValuedIndexedVector), intent(in) :: v1
    type(Vector), intent(in) :: v2
    substract_valuedIndexedVector_vector%point = v1%point
    substract_valuedIndexedVector_vector%val = v1%val
    substract_valuedIndexedVector_vector%x = v1%x - v2%x
    substract_valuedIndexedVector_vector%y = v1%y - v2%y
  end function substract_valuedIndexedVector_vector

  pure type(ValuedIndexedVector) function substract_valuedIndexedVector_indexedVector(v1,v2)
    type(ValuedIndexedVector), intent(in) :: v1
    type(IndexedVector), intent(in) :: v2
    substract_valuedIndexedVector_indexedVector%point = v1%point
    substract_valuedIndexedVector_indexedVector%val = v1%val
    substract_valuedIndexedVector_indexedVector%x = v1%x - v2%x
    substract_valuedIndexedVector_indexedVector%y = v1%y - v2%y
  end function substract_valuedIndexedVector_indexedVector

  pure function substract_vectorArray_vectorArray(arr1,arr2) result(arr)
    type(Vector), dimension(:), intent(in) :: arr1, arr2
    type(Vector), dimension(size(arr1)) :: arr
    integer :: i
    forall(i = 1:size(arr1)) arr(i) = arr1(i) + arr2(i)
  end function substract_vectorArray_vectorArray

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%% OPERATOR(*) %%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  pure type(Vector) function dot_real_vector(a,v2)
    real(kind=8), intent(in) :: a
    type(Vector), intent(in) :: v2
    dot_real_vector%x = v2%x*a
    dot_real_vector%y = v2%y*a
  end function dot_real_vector

  pure type(IndexedVector) function dot_real_indexedVector(a,v2)
    real(kind=8), intent(in) :: a
    type(IndexedVector), intent(in) :: v2
    dot_real_indexedVector%point = v2%point
    dot_real_indexedVector%x = v2%x*a
    dot_real_indexedVector%y = v2%y*a
  end function dot_real_indexedVector

  pure type(ValuedIndexedVector) function dot_real_valuedIndexedVector(a,v2)
    real(kind=8), intent(in) :: a
    type(ValuedIndexedVector), intent(in) :: v2
    dot_real_valuedIndexedVector%point = v2%point
    dot_real_valuedIndexedVector%val = v2%val
    dot_real_valuedIndexedVector%x = v2%x*a
    dot_real_valuedIndexedVector%y = v2%y*a
  end function dot_real_valuedIndexedVector

  pure type(Vector) function dot_vector_real(v1,a)
    type(Vector), intent(in) :: v1
    real(kind=8), intent(in) :: a
    dot_vector_real%x = v1%x*a
    dot_vector_real%y = v1%y*a
  end function dot_vector_real

  pure real(kind=8) function dot_vector_vector(v1,v2)
    type(Vector), intent(in) :: v1, v2
    dot_vector_vector = v1%x*v2%x + v1%y*v2%y
  end function dot_vector_vector

  pure real(kind=8) function dot_vector_indexedVector(v1,v2)
    type(Vector), intent(in) :: v1
    type(IndexedVector), intent(in) :: v2
    dot_vector_indexedVector = v1%x*v2%x + v1%y*v2%y
  end function dot_vector_indexedVector

  pure real(kind=8) function dot_vector_valuedIndexedVector(v1,v2)
    type(Vector), intent(in) :: v1
    type(ValuedIndexedVector), intent(in) :: v2
    dot_vector_valuedIndexedVector = v1%x*v2%x + v1%y*v2%y
  end function dot_vector_valuedIndexedVector

  pure type(IndexedVector) function dot_indexedVector_real(v1,a)
    type(IndexedVector), intent(in) :: v1
    real(kind=8), intent(in) :: a
    dot_indexedVector_real%point = v1%point
    dot_indexedVector_real%x = v1%x*a
    dot_indexedVector_real%y = v1%y*a
  end function dot_indexedVector_real

  pure real(kind=8) function dot_indexedVector_vector(v1,v2)
    type(IndexedVector), intent(in) :: v1
    type(Vector), intent(in) :: v2
    dot_indexedVector_vector = v1%x*v2%x + v1%y*v2%y
  end function dot_indexedVector_vector

  pure real(kind=8) function dot_indexedVector_indexedVector(v1,v2)
    type(IndexedVector), intent(in) :: v1
    type(IndexedVector), intent(in) :: v2
    dot_indexedVector_indexedVector = v1%x*v2%x + v1%y*v2%y
  end function dot_indexedVector_indexedVector

  pure real(kind=8) function dot_indexedVector_valuedIndexedVector(v1,v2)
    type(IndexedVector), intent(in) :: v1
    type(ValuedIndexedVector), intent(in) :: v2
    dot_indexedVector_valuedIndexedVector = v1%x*v2%x + v1%y*v2%y
  end function dot_indexedVector_valuedIndexedVector

  pure type(ValuedIndexedVector) function dot_valuedIndexedVector_real(v1,a)
    type(ValuedIndexedVector), intent(in) :: v1
    real(kind=8), intent(in) :: a
    dot_valuedIndexedVector_real%point = v1%point
    dot_valuedIndexedVector_real%val = v1%val
    dot_valuedIndexedVector_real%x = v1%x*a
    dot_valuedIndexedVector_real%y = v1%y*a
  end function dot_valuedIndexedVector_real

  pure real(kind=8) function dot_valuedIndexedVector_vector(v1,v2)
    type(ValuedIndexedVector), intent(in) :: v1
    type(Vector), intent(in) :: v2
    dot_valuedIndexedVector_vector = v1%x*v2%x + v1%y*v2%y
  end function dot_valuedIndexedVector_vector

  pure real(kind=8) function dot_valuedIndexedVector_indexedVector(v1,v2)
    type(ValuedIndexedVector), intent(in) :: v1
    type(IndexedVector), intent(in) :: v2
    dot_valuedIndexedVector_indexedVector = v1%x*v2%x + v1%y*v2%y
  end function dot_valuedIndexedVector_indexedVector

  pure real(kind=8) function dot_valuedIndexedVector_valuedIndexedVector(v1,v2)
    type(ValuedIndexedVector), intent(in) :: v1
    type(valuedIndexedVector), intent(in) :: v2
    dot_valuedIndexedVector_valuedIndexedVector = v1%x*v2%x + v1%y*v2%y
  end function dot_valuedIndexedVector_valuedIndexedVector

  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%% OPERATOR(/) %%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  pure type(Vector) function div_vector_real(v1,a)
    type(Vector), intent(in) :: v1
    real(kind=8), intent(in) :: a
    div_vector_real%x = v1%x/a
    div_vector_real%y = v1%y/a
  end function div_vector_real

  pure type(IndexedVector) function div_indexedVector_real(v1,a)
    type(IndexedVector), intent(in) :: v1
    real(kind=8), intent(in) :: a
    div_indexedVector_real%point = v1%point
    div_indexedVector_real%x = v1%x/a
    div_indexedVector_real%y = v1%y/a
  end function div_indexedVector_real

  pure type(ValuedIndexedVector) function div_valuedIndexedVector_real(v1,a)
    type(ValuedIndexedVector), intent(in) :: v1
    real(kind=8), intent(in) :: a
    div_valuedIndexedVector_real%point = v1%point
    div_valuedIndexedVector_real%val = v1%val
    div_valuedIndexedVector_real%x = v1%x/a
    div_valuedIndexedVector_real%y = v1%y/a
  end function div_valuedIndexedVector_real

  pure real(kind=8) function norm(v)
    type(Vector), intent(in) :: v
    norm = dsqrt(v%x*v%x + v%y*v%y)
  end function norm

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%% TYPE PROCEDURES %%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  pure subroutine wallType_addEdge(this, iedge)
    class(WallType), intent(inout) :: this
    type(IndexedEdge), intent(in) :: iedge

    integer :: nedge, npoin, i
    type(IndexedEdge), dimension(:), allocatable :: tmp_edges
    integer, dimension(:), allocatable :: tmp_points
    if(allocated(this%edges)) then
       nedge = size(this%edges) + 1
       allocate(tmp_edges(nedge))
       tmp_edges(1:nedge-1) = this%edges
       call move_alloc(tmp_edges, this%edges)
       this%edges(nedge) = iedge
    else
       allocate(this%edges(1))
       this%edges(1) = iedge
    end if
    ! UPDATE LIST OF POINTS
    if(allocated(this%points)) then
       do i = 1, size(iedge%points)
          if(.not.any(this%points == iedge%points(i))) then
             npoin = size(this%points) + 1
             allocate(tmp_points(npoin))
             tmp_points(1:npoin-1) = this%points
             call move_alloc(tmp_points, this%points)
             this%points(npoin) = iedge%points(i)
          end if
       end do
    else
       allocate(this%points(size(iedge%points)))
       this%points = iedge%points
    end if
  end subroutine wallType_addEdge

  pure subroutine wallType_updateNormals(this, coord)
    class(WallType), intent(inout) :: this
    type(Vector), dimension(:), intent(in) :: coord

    integer :: ipoi1, ipoi2, i, j
    type(Vector) :: v
    type(Edge) :: iedge
    if(.not.allocated(this%normals)) then
       allocate(this%normals(size(this%points)))
       this%normals(:)%point = this%points(:)
    end if
    forall(i = 1:size(this%normals)) &
         this%normals(i) = Vector(0.d0,0.d0)
    do i = 1, size(this%edges)
       iedge = this%edges(i)
       ipoi1 = iedge%points(1)
       ipoi2 = iedge%points(2)
       v%x = coord(ipoi2)%y - coord(ipoi1)%y
       v%y = coord(ipoi1)%x - coord(ipoi2)%x
       v = v/dsqrt(v*v)
       do j = 1, size(this%normals)
          if(this%normals(j)%point == ipoi1.or.&
               this%normals(j)%point == ipoi2) then
             this%normals(j) = this%normals(j) + v
          end if
       end do
    end do
    ! NORMALIZE VECTORS
    do i = 1, size(this%normals)
       this%normals(i) =&
            this%normals(i)/dsqrt(this%normals(i)*this%normals(i))
    end do
  end subroutine wallType_updateNormals
  
end module util_mod
