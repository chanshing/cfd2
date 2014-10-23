module functions
  implicit none
  type fargs
     real :: a, b, c
  end type fargs
contains
  real function f(args, x)
    class(*), intent(in) :: args
    real, intent(in) :: x

    real :: a, b, c
    select type(args)
    type is (fargs)
       a = args%a
       b = args%b
       c = args%c
       f = a*x**2 + b*x + c
    end select
  end function f
end module functions

module printer_mod
  implicit none
contains
  subroutine printer(f, args)
    interface
       real function f(args)
         class(*), intent(in) :: args
       end function f
    end interface
    class(*) :: args
    print*, f(args)
  end subroutine printer
end module printer_mod

module integrator_mod
  implicit none
contains
  real function integ(f, args, x0, x1, dx)
    interface
       real function f(args, x)
         class(*), intent(in) :: args
         real, intent(in) :: x
       end function f
    end interface
    class(*), intent(in) :: args
    real, intent(in) :: x0, x1, dx

    real :: x
    integ = 0
    x = x0
    do while(x < x1)
       integ = integ + f(args, x)*dx
       x = x + dx
    end do
  end function integ
end module integrator_mod

program tester
  use functions
  use printer_mod
  use integrator_mod
  type(fargs) :: args
  args = fargs(1, 0, 0)
  print*, integ(f, args, 0.0, 1.0, 0.00001)

end program tester
  
