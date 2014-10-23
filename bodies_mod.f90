module bodies_mod
  use mesh_mod
  private
  public :: bodies
  public :: init_bodies

  ! EDIT: NUMERO DE CUERPOS
  integer, parameter :: nbodies = 2
  type(BodyContainer), dimension(nbodies) :: bodies 

  ! EDIT: DECLARACION DE LOS CUERPOS
  type, extends(Body) :: BodyType1
   contains
     procedure, nopass :: alpha_t => alpha_t1
  end type BodyType1
  type(BodyType1), save :: body1

  type, extends(Body) :: BodyType2
   contains
     procedure, nopass, private :: alpha_t => alpha_t2
  end type BodyType2
  type(BodyType2), save :: body2

contains

  ! EDIT: INICIALIZAR ARREGLO DE CUERPOS
  subroutine init_bodies()
    allocate(bodies(1)%obj, source=body1)
    allocate(bodies(2)%obj, source=body2)
  end subroutine init_bodies

  ! EDIT: DESCRIPCION DEL MOVIMIENTO
  pure real(kind=8) function alpha_t1(t)
    real(kind=8), intent(in) :: t
    alpha_t1 = 0
  end function alpha_t1
  
  pure real(kind=8) function alpha_t2(t)
    real(kind=8), intent(in) :: t
    alpha_t2 = 0
  end function alpha_t2
  
end module bodies_mod
