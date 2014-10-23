module lala
  private
  public :: one

  type one
   contains
     procedure, public :: fun
  end type one

contains

  subroutine fun(this)
    class(one) :: this
    call pfun()
  end subroutine fun

  subroutine pfun()
    print*, "Yeah"
  end subroutine pfun

end module lala

module lele
  use lala
  private
  public :: lalar
contains
  subroutine lalar(wan)
    type(one) :: wan
    call wan%fun
  end subroutine lalar
end module lele

program lol
  use lala
  use lele
  type(one) :: wan
  call lalar(wan)
end program lol
