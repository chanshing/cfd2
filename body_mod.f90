module body_mod
  ! %%%%%%%%%% TIPO DERIVADO "Body" %%%%%%%%%%
  ! === METODOS PUBLICOS ===: 
  ! * initialize(points[, rotationCenter, t])
  ! * rotate(angle, coord)
  ! * translate(r, coord)
  ! * computeBoundaryDisplacement(dt, coord, ds)
  ! === OVERRIDABLE ===:
  ! * r_t(t)
  ! * alpha_t(t)
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  use util_mod
  implicit none
  private
  public :: Body

  type Body
     integer, dimension(:), allocatable :: points
     type(Vector) :: rotationCenter = Vector(0,0)
     real(kind=8) :: t = 0.d0
     integer :: nsize
   contains
     procedure, public :: rotate
     procedure, public :: translate
     procedure, private, nopass :: r_t
     procedure, private, nopass :: alpha_t
     procedure, public :: computeBoundaryDisplacement
     procedure, public :: initialize
  end type Body

contains

  subroutine initialize(b, points, rotationCenter, t)
    class(Body) :: b
    integer, dimension(:), intent(in) :: points
    type(Vector), intent(in), optional :: rotationCenter
    real(kind=8), intent(in), optional :: t
    
    b%nsize = size(points)
    allocate(b%points(b%nsize))
    b%points(:) = points(:)
    if(present(rotationCenter)) b%rotationCenter = rotationCenter
    if(present(t)) b%t = t
  end subroutine initialize

  subroutine rotate(b, angle, coord)
    class(Body) :: b
    real(kind=8), intent(in) :: angle
    type(Vector), dimension(:), intent(inout) :: coord

    type(Vector) :: r
    integer :: i
    do i = 1, b%nsize
       r = coord(b%points(i)) - b%rotationCenter
       coord(b%points(i))%x = b%rotationCenter%x + &
            r%x*dcos(angle) - r%y*dsin(angle)
       coord(b%points(i))%y = b%rotationCenter%y + &
            r%x*dsin(angle) + r%y*dcos(angle)
    end do
  end subroutine rotate

  subroutine translate(b, r, coord)
    class(Body) :: b
    type(Vector), intent(in) :: r
    type(Vector), dimension(:), intent(inout) :: coord

    integer :: i
    do i = 1, b%nsize
       coord(b%points(i)) = coord(b%points(i)) + r
    end do
  end subroutine translate

  subroutine computeBoundaryDisplacement(b, dt, coord, ds)
    class(Body) :: b
    real(kind=8), intent(in) :: dt
    type(Vector), dimension(:), intent(in) :: coord
    type(IndexedVector), dimension(b%nsize), intent(out) :: ds

    real(kind=8) :: t
    type(Vector) :: r, r_new
    real(kind=8) :: dalpha
    type(Vector) :: dr
    integer :: i

    t = b%t

    ! DESPLAZAMIENTO DE LOS BORDES POR ROTACION
    dalpha = b%alpha_t(t + dt) - b%alpha_t(t)
    do i = 1, b%nsize
       r = coord(b%points(i)) - b%rotationCenter
       r_new%x = r%x*dcos(dalpha) - r%y*dsin(dalpha)
       r_new%y = r%x*dsin(dalpha) + r%y*dcos(dalpha)
       ds(i) = IndexedVector(&
            idx = b%points(i), &
            x = r_new%x - r%x, &
            y = r_new%y - r%y)
    end do

    ! DESPLAZAMIENTO DE LOS BORDES POR TRANSLACION
    dr = b%r_t(t + dt) - b%r_t(t)
    do i = 1, b%nsize
       ds(i) = ds(i) + dr
    end do

    ! UPDATE TIME
    b%t = t + dt
  end subroutine computeBoundaryDisplacement

  real(kind=8) function alpha_t(t)
    real(kind=8), intent(in) :: t
    alpha_t = 0
  end function alpha_t

  type(Vector) function r_t(t)
    real(kind=8), intent(in) :: t
    r_t = Vector(0,0)
  end function r_t

end module body_mod
