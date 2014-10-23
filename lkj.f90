program lkj
  use util_mod
  type(Vector) :: v 
  v = Vector(1.d0,2.d0)
  print*, norm(v)
  print*, dabs(v%x*v%x + v%y*v%y)
end program lkj
