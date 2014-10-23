module mesh_mod
  use util_mod
  implicit none
  public
  private :: meshType_update,&
       body_alpha_t,&
       body_r_t

  type, extends(Triangle) :: ElementType
     type(Vector), dimension(3) :: dshape
  end type ElementType

  type, extends(WallType) :: Body
     type(Vector) :: rotationCenter = Vector(0,0)
     real(kind=8) :: time = 0.d0
   contains
     procedure, private, nopass :: r_t => body_r_t
     procedure, private, nopass :: alpha_t => body_alpha_t
  end type Body

  type BodyContainer
     class(Body), pointer :: obj
  end type BodyContainer

  type MeshType
     integer :: nelem, npoin
     integer, dimension(:), allocatable :: fixedPoints
     type(Vector), dimension(:), allocatable :: coord
     type(Vector), dimension(:), allocatable :: vel
     type(ElementType), dimension(:), allocatable :: elements
     type(BodyContainer), dimension(:), allocatable :: bodies
   contains
     procedure, public :: update => meshType_update
  end type MeshType

  type DistantConditionsType
     type(Vector) :: velocity
     real(kind=8) :: temperature
     real(kind=8) :: density
     real(kind=8) :: mach
     real(kind=8) :: pressure
     real(kind=8) :: sound_speed
  end type DistantConditionsType

  type BoundaryConditionsType
     type(IndexedVector), dimension(:), allocatable :: &
          velocities
     type(IndexedScalar), dimension(:), allocatable :: &
          temperatures
     type(IndexedScalar), dimension(:), allocatable :: &
          densities
     type(WallType) :: &
          normalVelocities
  end type BoundaryConditionsType

  type SolverParamsType
     integer :: max_iterations
     integer :: isMoving
     integer :: isLocalTimeStep
     real(kind=8) :: shock_factor
     real(kind=8) :: timestep_factor
  end type SolverParamsType

  type ReferenceValuesType
     type(Vector) :: gravity
     real(kind=8) :: viscosity
     real(kind=8) :: gas_constant
     real(kind=8) :: thermal_conductivity
     real(kind=8) :: specific_heat_capacity
     real(kind=8) :: heat_capacity_ratio
     type(Vector) :: velocity
     real(kind=8) :: temperature
     real(kind=8) :: density
     real(kind=8) :: mach
     real(kind=8) :: pressure
     real(kind=8) :: sound_speed
     real(kind=8) :: total_temperature
     real(kind=8) :: total_density
     real(kind=8) :: total_pressure
  end type ReferenceValuesType

  type PrimitiveVarsType
     real(kind=8) :: density
     type(Vector) :: velocity
     real(kind=8) :: energy
     real(kind=8) :: pressure
     real(kind=8) :: temperature
     real(kind=8) :: mach
     real(kind=8) :: sound_speed
  end type PrimitiveVarsType

contains

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%% TYPE PROCEDURES %%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  pure real(kind=8) function body_alpha_t(t)
    real(kind=8), intent(in) :: t
    body_alpha_t = 0
  end function body_alpha_t

  pure type(Vector) function body_r_t(t)
    real(kind=8), intent(in) :: t
    body_r_t = Vector(0,0)
  end function body_r_t

  subroutine meshType_update(this)
    class(MeshType), intent(inout) :: this
    integer :: i
    integer :: ipoi1, ipoi2, ipoi3
    type(Vector) :: p1, p2, p3
    real(kind=8) :: area

    do i = 1, this%nelem
       ipoi1 = this%elements(i)%points(1)
       ipoi2 = this%elements(i)%points(2)
       ipoi3 = this%elements(i)%points(3)
       p1 = this%coord(ipoi1)
       p2 = this%coord(ipoi2)
       p3 = this%coord(ipoi3)
       ! AREA DEL ELEMENTO
       area = 0.5*(&
            (p2%x*p3%y + p3%x*p1%y + p1%x*p2%y) -&
            (p2%x*p1%y + p3%x*p2%y + p1%x*p3%y))
       this%elements(i)%area = area

       ! DERIVADAS DE LAS FUNCIONES DE FORMA
       this%elements(i)%dshape(1) = &
            Vector(p2%y - p3%y, p3%x - p2%x)/(2*area)
       this%elements(i)%dshape(2) = &
            Vector(p3%y - p1%y, p1%x - p3%x)/(2*area)
       this%elements(i)%dshape(3) = &
            Vector(p1%y - p2%y, p2%x - p1%x)/(2*area)
    end do

    do i = 1, size(this%bodies)
       ! ACTUALIZA NORMALES AL CUERPO
       call this%bodies(i)%obj%updateNormals(this%coord)
    end do

  end subroutine meshType_update

end module mesh_mod
