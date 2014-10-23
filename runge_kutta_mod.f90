module runge_kutta_mod
  implicit none
  private
  real(kind=8), dimension(:,:), allocatable, save :: U_old
  real(kind=8), dimension(:,:), allocatable, save :: rhs
  public:: rk4, rk_config
contains
  subroutine rk4(U, f, dt)
    real(kind=8), dimension(:,:), intent(inout) :: U
    real(kind=8), intent(in) :: dt
    interface
       subroutine f(rhs, U)
         real(kind=8), dimension(:,:), intent(in) :: U
         real(kind=8), dimension(:,:), intent(out) :: rhs
       end subroutine f
    end interface
    integer :: stage
    real(kind=8) :: alpha
    call error_check()
    U_old = U
    do stage = 1, 4
       alpha = 1.d0/(5 - stage)
       call f(rhs, U)
       U = U_old + alpha*dt*rhs
    end do
  contains
    subroutine error_check()
      if(.not.allocated(U_old).or..not.allocated(rhs)) then
         stop "U_old o rhs NO ASIGNADOS"
      end if
      if(any(shape(U_old)/=shape(U))&
           .or.any(shape(rhs)/=shape(U))) then
         stop "U_old o rhs NO CONFORMA"
      end if
    end subroutine error_check
  end subroutine rk4

  subroutine rk_config(U)
    real(kind=8), dimension(:,:), intent(in) :: U
    integer :: m, n
    m = size(U,1)
    n = size(U,2)
    if(allocated(U_old)) deallocate(U_old)
    allocate(U_old(m,n))
    if(allocated(rhs)) deallocate(rhs)
    allocate(rhs(m,n))
  end subroutine rk_config
end module runge_kutta_mod
