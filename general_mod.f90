module general_mod
  use mesh_mod
  implicit none
  private

  public ::&
       U_to_prim,&
       restoreBC,&
       initCond

contains

  pure subroutine U_to_prim(pVars, U, refVal)
    type(PrimitiveVarsType), dimension(:), intent(out) :: pVars
    type(ReferenceValuesType), intent(in) :: refVal
    real(kind=8), dimension(:,:), intent(in) :: U
    integer :: ipoin
    type(Vector) :: v
    real(kind=8) :: rho, e, P, T, mach, c, v_sq
    real(kind=8) :: gamma, R
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gamma = refVal%heat_capacity_ratio
    R = refVal%gas_constant
    do ipoin = 1, size(U,2)
       rho = U(1,ipoin)
       v%x = U(2,ipoin)/rho
       v%y = U(3,ipoin)/rho
       e = U(4,ipoin)/rho
       v_sq = (v%x**2 + v%y**2)
       P = rho*(gamma - 1.d0)*(e - .5d0*v_sq)
       T = P/(rho*R)
       mach = dsqrt(v_sq/(T*gamma*R))
       c = dsqrt(gamma*R*T)
       pVars(ipoin) = PrimitiveVarsType(&
            density=rho,&
            velocity=v,&
            energy=e,&
            pressure=P,&
            temperature=T,&
            mach=mach,&
            sound_speed=c)
    end do
  end subroutine U_to_prim

  pure subroutine restoreBC(U, bcond, mesh, refVal)
    real(kind=8), dimension(:,:), intent(inout) :: U
    type(BoundaryConditionsType), intent(in) :: bcond
    type(MeshType), intent(in) :: mesh
    type(ReferenceValuesType), intent(in) :: refVal
    integer :: i, j, ipoin
    real(kind=8) :: rho, v1, v2, e, T, Cv
    type(Vector) :: v, w, n
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! DENSITY
    do i = 1, size(bcond%densities)
       ipoin = bcond%densities(i)%point
       rho = U(1,ipoin)
       v1 = U(2,ipoin)/rho
       v2 = U(3,ipoin)/rho
       e = U(4,ipoin)/rho
       rho = bcond%densities(i)%val
       U(1,ipoin) = rho
       U(2,ipoin) = rho*v1
       U(3,ipoin) = rho*v2
       U(4,ipoin) = rho*e
    end do

    ! VELOCITY
    do i = 1, size(bcond%velocities)
       ipoin = bcond%velocities(i)%point
       rho = U(1,ipoin)
       v1 = U(2,ipoin)/rho
       v2 = U(3,ipoin)/rho
       v1 = bcond%velocities(i)%x
       v2 = bcond%velocities(i)%y
       U(2,ipoin) = rho*v1
       U(3,ipoin) = rho*v2
    end do

    ! TEMPERATURE
    Cv = refVal%specific_heat_capacity
    do i = 1, size(bcond%temperatures)
       ipoin = bcond%temperatures(i)%point
       T = bcond%temperatures(i)%val
       rho = U(1,ipoin)
       v1 = U(2,ipoin)/rho
       v2 = U(3,ipoin)/rho
       e = Cv*T + 0.5d0*(v1**2 + v2**2)
       U(4,ipoin) = rho*e
    end do

    ! VELOCITY ON BODY SURFACE
    do i = 1, size(mesh%bodies)
       do j = 1, size(mesh%bodies(i)%obj%normals)
          ipoin = mesh%bodies(i)%obj%normals(j)%point
          n = mesh%bodies(i)%obj%normals(j)
          rho = U(1,ipoin)
          v1 = U(2,ipoin)/rho
          v2 = U(3,ipoin)/rho
          v = Vector(v1,v2)
          w = mesh%vel(ipoin)
          v = v - (v*n)*n + (w*n)*n
          U(2,ipoin) = rho*v%x
          U(3,ipoin) = rho*v%y
       end do
    end do

    ! VELOCITY ON BODY SURFACE
    ! do i = 1, size(mesh%bodies)
    !    do j = 1, size(mesh%bodies(i)%obj%normals)
    !       ipoin = mesh%bodies(i)%obj%normals(j)%point
    !       n = mesh%bodies(i)%obj%normals(j)
    !       rho = U(1,ipoin)
    !       w = mesh%vel(ipoin)
    !       v = w
    !       U(2,ipoin) = rho*v%x
    !       U(3,ipoin) = rho*v%y
    !    end do
    ! end do
  end subroutine restoreBC

  subroutine initCond(U, bCond, mesh, refVal)
    real(kind=8), dimension(:,:), intent(out) :: U
    type(MeshType), intent(inout) :: mesh
    type(BoundaryConditionsType), intent(in) :: bCond
    type(ReferenceValuesType), intent(in) :: refVal
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    real(kind=8) :: rho, v1, v2, e, Cv, T
    rho = refVal%density
    v1 = refVal%velocity%x
    v2 = refVal%velocity%y
    T = refVal%temperature
    Cv = refVal%specific_heat_capacity
    e = Cv*T + 0.5d0*(v1**2 + v2**2)
    U(1,:) = rho
    U(2,:) = rho*v1
    U(3,:) = rho*v2
    U(4,:) = rho*e
    call restoreBC(U, bCond, mesh, refVal)

  end subroutine initCond

  pure real(kind=8) function density_by_mach(mach, refVal)
    real(kind=8), intent(in) :: mach
    type(ReferenceValuesType), intent(in) :: refVal
    real(kind=8) :: gamma, rho_tot
    gamma = refVal%heat_capacity_ratio
    rho_tot = refVal%total_density
    density_by_mach = &
         rho_tot/((1 + 0.5*(gamma - 1)*mach**2)**(1/(gamma - 1)))
  end function density_by_mach

  pure real(kind=8) function pressure_by_mach(mach, refVal)
    real(kind=8), intent(in) :: mach
    type(ReferenceValuesType), intent(in) :: refVal
    real(kind=8) :: gamma, P_tot
    gamma = refVal%heat_capacity_ratio
    P_tot = refVal%total_pressure
    pressure_by_mach = &
         P_tot/((1 + 0.5*(gamma - 1)*mach**2)**(gamma/(gamma - 1)))
  end function pressure_by_mach

  pure real(kind=8) function temperature_by_mach(mach, refVal)
    real(kind=8), intent(in) :: mach
    type(ReferenceValuesType), intent(in) :: refVal
    real(kind=8) :: gamma, T_tot
    gamma = refVal%heat_capacity_ratio
    T_tot = refVal%total_temperature
    temperature_by_mach = &
         T_tot/((1 + 0.5*(gamma - 1)*mach**2))
  end function temperature_by_mach

end module general_mod
