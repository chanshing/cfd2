module solver_mod
  use mesh_mod
  use general_mod
  implicit none
  private
  save

  public :: &
       compute_rhs,&
       compute_rhs_config

  logical :: isMeshMove = .false.
  real(kind=8), dimension(:), allocatable :: mass
  type(PrimitiveVarsType), dimension(:), allocatable :: pVars

  ! ARGUMENTOS
  type(MeshType), pointer :: mesh
  type(ReferenceValuesType), pointer :: refVal
  real(kind=8), pointer :: dt

contains

  subroutine compute_rhs_config &
       (mesh0, refVal0, dt0, isMeshMove0)
    type(MeshType), intent(in), target :: mesh0
    type(ReferenceValuesType), intent(in), target :: refVal0
    real(kind=8), intent(in), target :: dt0
    logical, intent(in), optional :: isMeshMove0
    mesh => mesh0
    refVal => refVal0
    dt => dt0
    if(present(isMeshMove0)) isMeshMove = isMeshMove0
    if(.not.allocated(pVars)) allocate(pVars(mesh%npoin))
    if(.not.allocated(mass)) allocate(mass(mesh%npoin))
    call compute_mass()
  end subroutine compute_rhs_config

  subroutine compute_rhs(rhs, U)
    real(kind=8), dimension(:,:), intent(out) :: rhs
    real(kind=8), dimension(:,:), intent(in) :: U
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    real(kind=8), dimension(4,3) :: rhs_local
    real(kind=8), dimension(4) :: U_local
    real(kind=8) :: rho, v1, v2, e, w1, w2
    real(kind=8) :: area
    real(kind=8) :: tau_r, tau_v, tau_e, nu_shoc
    real(kind=8) :: gamma, Cv, k, mu, T
    integer :: ielem, ipoi1, ipoi2, ipoi3, i
    real(kind=8), dimension(3) :: Nx, Ny
    real(kind=8), dimension(4) :: Ux, Uy
    real(kind=8), dimension(1:4) :: AiUi
    real(kind=8), dimension(1:4) :: A1AiUi
    real(kind=8), dimension(1:4) :: A2AiUi
    real(kind=8), parameter, dimension(3,3) :: N = &
         reshape((/&
         0.d0, .5d0, .5d0, &
         .5d0, 0.d0, .5d0, &
         .5d0, .5d0, 0.d0 &
         /), (/3,3/))

    if(isMeshMove) call compute_mass()

    call U_to_prim(pVars, U, refVal)

    rhs(:,:) = 0.d0
    gamma = refVal%heat_capacity_ratio
    Cv = refVal%specific_heat_capacity
    do ielem = 1, mesh%nelem
       rhs_local(:,:) = 0.d0
       area = mesh%elements(ielem)%area
       ipoi1 = mesh%elements(ielem)%points(1)
       ipoi2 = mesh%elements(ielem)%points(2)
       ipoi3 = mesh%elements(ielem)%points(3)
       Nx(1:3) = mesh%elements(ielem)%dshape(:)%x
       Ny(1:3) = mesh%elements(ielem)%dshape(:)%y
       Ux(1:4) = &
            U(1:4,ipoi1)*Nx(1) + &
            U(1:4,ipoi2)*Nx(2) + &
            U(1:4,ipoi3)*Nx(3)
       Uy(1:4) = &
            U(1:4,ipoi1)*Ny(1) + &
            U(1:4,ipoi2)*Ny(2) + &
            U(1:4,ipoi3)*Ny(3)
       T = (pVars(ipoi1)%temperature + &
            pVars(ipoi2)%temperature + &
            pVars(ipoi3)%temperature)/3.d0
       mu = get_mu(refVal, T)
       k = get_k(refVal, T)

       call compute_stability(tau_r, tau_v, tau_e, nu_shoc, &
            ielem, mesh, pVars, refVal, mu, dt)

       do i = 1, 3
          U_local(1:4) = &
               N(1,i)*U(1:4,ipoi1) +&
               N(2,i)*U(1:4,ipoi2) +&
               N(3,i)*U(1:4,ipoi3)
          rho = U_local(1)
          v1 = U_local(2)/rho
          v2 = U_local(3)/rho
          e = U_local(4)/rho
          w1 = &
               N(1,i)*mesh%vel(ipoi1)%x +&
               N(2,i)*mesh%vel(ipoi2)%x +&
               N(3,i)*mesh%vel(ipoi3)%x
          w2 = &
               N(1,i)*mesh%vel(ipoi1)%y +&
               N(2,i)*mesh%vel(ipoi2)%y +&
               N(3,i)*mesh%vel(ipoi3)%y 

          AiUi(1) = Ux(2) + Uy(3)
          AiUi(2) = -1.0d0/2.0d0*Ux(1)*(2*v1**2 - (gamma - 1)*(v1**2 + &
               v2**2)) - Ux(2)*v1*(gamma - 3) - Ux(3)*v2*(gamma - 1) + &
               Ux(4)*(gamma - 1) - Uy(1)* v1*v2 + Uy(2)*v2 + Uy(3)*v1
          AiUi(3) = -Ux(1)*v1*v2 + Ux(2)*v2 + Ux(3)*v1 -               &
               1.0d0/2.0d0*Uy(1)*(2*v2**2 - ( gamma - 1)*(v1**2 +      &
               v2**2)) - Uy(2)*v1*(gamma - 1) - Uy(3)*v2*(gamma - 3) + &
               Uy(4)*(gamma - 1)
          AiUi(4) = -Ux(1)*v1*(e*gamma - (gamma - 1)*(v1**2 + v2**2))  &
               - 1.0d0/ 2.0d0*Ux(2)*(-2*e*gamma + 2*v1**2*(gamma - 1)  &
               + (gamma - 1)*(v1**2 + v2**2)) - Ux(3)*v1*v2*(gamma -   &
               1) + Ux(4)*gamma*v1 - Uy(1)*v2*(e* gamma - (gamma -     &
               1)*(v1**2 + v2**2)) - Uy(2)*v1*v2*(gamma - 1) -         &
               1.0d0/2.0d0*Uy(3)*(-2*e*gamma + 2*v2**2*(gamma - 1) +   &
               (gamma - 1)*( v1**2 + v2**2)) + Uy(4)*gamma*v2

          A1AiUi(1) = -w1*AiUi(1) + AiUi(2)
          A1AiUi(2) = -v2*(gamma - 1)*AiUi(3) + (-w1 + v1*(-gamma +    &
               3))*AiUi (2) + (gamma - 1)*AiUi(4) + (-v1**2 +          &
               (1.0d0/2.0d0)*(gamma - 1)*(v1**2 + v2**2))*AiUi(1)
          A1AiUi(3) = -v1*v2*AiUi(1) + v2*AiUi(2) + (-w1 + v1)*AiUi(3) 
          A1AiUi(4) = -v1*v2*(gamma - 1)*AiUi(3) + v1*(-e*gamma +      &
               (gamma - 1 )*(v1**2 + v2**2))*AiUi(1) + (-w1 +          &
               gamma*v1)*AiUi(4) + (e* gamma - v1**2*(gamma - 1) -     &
               1.0d0/2.0d0*(gamma - 1)*(v1**2 + v2** 2))*AiUi(2)

          A2AiUi(1) = -w2*AiUi(1) + AiUi(3)
          A2AiUi(2) = -v1*v2*AiUi(1) + v1*AiUi(3) + (-w2 + v2)*AiUi(2)
          A2AiUi(3) = -v1*(gamma - 1)*AiUi(2) + (-w2 + v2*(-gamma +    &
               3))*AiUi (3) + (gamma - 1)*AiUi(4) + (-v2**2 +          &
               (1.0d0/2.0d0)*(gamma - 1)*(v1**2 + v2**2))*AiUi(1)
          A2AiUi(4) = -v1*v2*(gamma - 1)*AiUi(2) + v2*(-e*gamma +      &
               (gamma - 1 )*(v1**2 + v2**2))*AiUi(1) + (-w2 +          &
               gamma*v2)*AiUi(4) + (e* gamma - v2**2*(gamma - 1) -     &
               1.0d0/2.0d0*(gamma - 1)*(v1**2 + v2** 2))*AiUi(3)

          A1AiUi(1) = A1AiUi(1)*tau_r
          A1AiUi(2:3) = A1AiUi(2:3)*tau_v
          A1AiUi(4) = A1AiUi(4)*tau_e
          A2AiUi(1) = A2AiUi(1)*tau_r
          A2AiUi(2:3) = A2AiUi(2:3)*tau_v
          A2AiUi(4) = A2AiUi(4)*tau_e

          rhs_local(:,1) = rhs_local(:,1) &
               + N(1,i)*AiUi &
               + (Nx(1)*A1AiUi + Ny(1)*A2AiUi) &
               + nu_shoc*(Nx(1)*Ux + Ny(1)*Uy) &
               - N(1,i)*(w1*Ux + w2*Uy)
          rhs_local(:,2) = rhs_local(:,2) &
               + N(2,i)*AiUi &
               + (Nx(2)*A1AiUi + Ny(2)*A2AiUi) &
               + nu_shoc*(Nx(2)*Ux + Ny(2)*Uy) &
               - N(2,i)*(w1*Ux + w2*Uy)
          rhs_local(:,3) = rhs_local(:,3) &
               + N(3,i)*AiUi &
               + (Nx(3)*A1AiUi + Ny(3)*A2AiUi) &
               + nu_shoc*(Nx(3)*Ux + Ny(3)*Uy) &
               - N(3,i)*(w1*Ux + w2*Uy)

          ! call rhs_diffusive(rhs_local, &
          !      Nx, Ny, Ux, Uy, &
          !      rho, v1, v2, e, &
          !      Cv, k, mu)
       end do

       rhs_local = rhs_local*area/3.d0

       do i = 1, 4
          rhs(i,ipoi1) = rhs(i,ipoi1) + rhs_local(i,1)
          rhs(i,ipoi2) = rhs(i,ipoi2) + rhs_local(i,2)
          rhs(i,ipoi3) = rhs(i,ipoi3) + rhs_local(i,3)
       end do
       
    end do

    do i = 1, mesh%npoin
       rhs(:,i) = -rhs(:,i)/mass(i)
    end do

  end subroutine compute_rhs

  pure subroutine rhs_diffusive(rhs_local, &
       Nx, Ny, Ux, Uy, &
       rho, v1, v2, e, &
       Cv, k, mu)
    real(kind=8), dimension(4,3), intent(inout) :: rhs_local
    real(kind=8), dimension(3), intent(in) :: Nx, Ny
    real(kind=8), dimension(4), intent(in) :: Ux, Uy
    real(kind=8), intent(in) :: rho, v1, v2, e
    real(kind=8), intent(in) :: Cv, k, mu
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    real(kind=8), dimension(2:4) :: K1jUj
    real(kind=8), dimension(2:4) :: K2jUj
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! K1jUj(1) = 0
    K1jUj(2) = (2.0d0/3.0d0)*mu*(-2*Ux(1)*v1 + 2*Ux(2) + Uy(1)*v2 - Uy(3))/rho
    K1jUj(3) = mu*(-Ux(1)*v2 + Ux(3) - Uy(1)*v1 + Uy(2))/rho
    K1jUj(4) = (1.0d0/3.0d0)*(Cv*mu*(-Uy(1)*v1*v2 + 3*Uy(2)*v2 - 2*Uy(3)*v1) &
         - Ux(1)*(Cv*mu*(3*(v1**2 + v2**2) + v1**2) - 3*k*((v1**2 + v2**2) - e)) &
         + Ux(2)*v1*(4*Cv*mu - 3*k) + 3*Ux(3)*v2*(Cv*mu - k) + 3*Ux(4)*k)/(Cv*rho)

    ! K2jUj(1) = 0
    K2jUj(2) = mu*(-Ux(1)*v2 + Ux(3) - Uy(1)*v1 + Uy(2))/rho
    K2jUj(3) = (2.0d0/3.0d0)*mu*(Ux(1)*v1 - Ux(2) - 2*Uy(1)*v2 + 2*Uy(3))/rho
    K2jUj(4) = (1.0d0/3.0d0)*(Cv*mu*(-Ux(1)*v1*v2 - 2*Ux(2)*v2 + 3*Ux(3)*v1) &
         - Uy(1)*(Cv*mu*(3*(v1**2 + v2**2) + v2**2) - 3*k*((v1**2 + v2**2) - e)) &
         + 3*Uy(2)*v1*(Cv*mu - k) + Uy(3)*v2*(4*Cv*mu - 3*k) + 3*Uy(4)*k)/(Cv*rho)

    rhs_local(2:4,1) = rhs_local(2:4,1) + (Nx(1)*K1jUj(2:4) + Ny(1)*K2jUj(2:4))
    rhs_local(2:4,2) = rhs_local(2:4,2) + (Nx(2)*K1jUj(2:4) + Ny(2)*K2jUj(2:4))
    rhs_local(2:4,3) = rhs_local(2:4,3) + (Nx(3)*K1jUj(2:4) + Ny(3)*K2jUj(2:4))
  end subroutine rhs_diffusive

  pure real(kind=8) function get_mu(refVal, T)
    type(ReferenceValuesType), intent(in) :: refVal
    real(kind=8), intent(in) :: T
    real(kind=8) :: mu_ref, T_ref
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mu_ref = refVal%viscosity
    T_ref = refVal%temperature
    get_mu = mu_ref*(T/T_ref)**1.5d0*(T_ref + 110)/(T + 110)
  end function get_mu

  pure real(kind=8) function get_k(refVal, T)
    type(ReferenceValuesType), intent(in) :: refVal
    real(kind=8), intent(in) :: T
    real(kind=8) :: k_ref, T_ref
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k_ref = refVal%thermal_conductivity
    T_ref = refVal%temperature
    get_k = k_ref*(T/T_ref)**1.5d0*(T_ref + 194)/(T + 194)
  end function get_k

  pure subroutine compute_stability(tau_r, tau_v, tau_e, nu_shoc, &
       ielem, mesh, pVars, refVal, mu, dt)
    real(kind=8), intent(out) :: tau_r, tau_v, tau_e, nu_shoc
    integer, intent(in) :: ielem
    type(MeshType), intent(in) :: mesh
    type(PrimitiveVarsType), dimension(:), intent(in) :: pVars
    type(ReferenceValuesType), intent(in) :: refVal
    real(kind=8), intent(in) :: mu, dt
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    integer :: points(3)
    real(kind=8) :: rho(3), T(3), c(3), v_ref, rho_ref
    type(Vector) :: v(3)
    type(Vector) :: grad_N(3)
    type(Vector) :: j, grad_rho, r, r_e
    real(kind=8) :: inv_tau1, inv_tau2, inv_tau3_v, inv_tau3_e, tau_shoc
    real(kind=8) :: h_rgn, h_rgn_e, h_shoc
    real(kind=8), parameter :: beta = 1
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    points(:) = mesh%elements(ielem)%points(:)
    rho(:) = pVars(points(:))%density
    T(:) = pVars(points(:))%temperature
    c(:) = pVars(points(:))%sound_speed
    v(1) = pVars(points(1))%velocity - mesh%vel(points(1))
    v(2) = pVars(points(2))%velocity - mesh%vel(points(2))
    v(3) = pVars(points(3))%velocity - mesh%vel(points(3))
    v_ref = norm(&
         refVal%velocity - &
         (mesh%vel(points(1)) + &
         mesh%vel(points(2)) + &
         mesh%vel(points(3)))/3.d0)
    rho_ref = refVal%density
    grad_N(:) = mesh%elements(ielem)%dshape(:)
    
    ! 1/TAU1
    grad_rho = &
         rho(1)*grad_N(1) + &
         rho(2)*grad_N(2) + &
         rho(3)*grad_N(3)
    if(norm(grad_rho) > tiny(0d0)) then
       j = grad_rho/norm(grad_rho)
       inv_tau1 = &
            c(1)*dabs(j*grad_N(1)) + dabs(v(1)*grad_N(1)) + &
            c(2)*dabs(j*grad_N(2)) + dabs(v(2)*grad_N(2)) + &
            c(3)*dabs(j*grad_N(3)) + dabs(v(3)*grad_N(3))
    else
       ! inv_tau1 = &
       !      dabs(v(1)*grad_N(1)) + &
       !      dabs(v(2)*grad_N(2)) + &
       !      dabs(v(3)*grad_N(3))
       inv_tau1 = 0.d0
    end if

    ! 1/TAU2
    inv_tau2 = 2/dt

    ! 1/TAU3_V
    ! r = &
    !      norm(v(1))*grad_N(1) + &
    !      norm(v(2))*grad_N(2) + &
    !      norm(v(3))*grad_N(3)
    ! if(norm(r) > tiny(0d0)) then
    !    r = r/norm(r)
    !    h_rgn = &
    !         dabs(r*grad_N(1)) + &
    !         dabs(r*grad_N(2)) + &
    !         dabs(r*grad_N(3))
    !    h_rgn = 2.d0/h_rgn
    !    inv_tau3_v = 2/(0.25d0*h_rgn**2/mu)
    ! else
    !    inv_tau3_v = 0.d0
    ! end if
    ! DEBUGGG
    inv_tau3_v = 0.d0

    ! 1/TAU3_E
    ! r_e = &
    !      T(1)*grad_N(1) + &
    !      T(2)*grad_N(2) + &
    !      T(3)*grad_N(3)
    ! if(norm(r_e) > tiny(0d0)) then
    !    r_e = r_e/norm(r_e)
    !    h_rgn_e = &
    !         dabs(r_e*grad_N(1)) + &
    !         dabs(r_e*grad_N(2)) + &
    !         dabs(r_e*grad_N(3))
    !    h_rgn_e = 2.d0/h_rgn_e
    !    inv_tau3_e = 1/(0.25d0*h_rgn_e**2/mu)
    ! else
    !    inv_tau3_e = 0.d0
    ! end if
    ! DEBUGGG
    inv_tau3_e = 0.d0

    tau_r = (1/(inv_tau1**2 + inv_tau2**2))**0.5
    tau_v = (1/(inv_tau1**2 + inv_tau2**2 + inv_tau3_v**2))**0.5
    tau_e = (1/(inv_tau1**2 + inv_tau2**2 + inv_tau3_e**2))**0.5

    ! NU_SHOC
    if(norm(grad_rho) > tiny(0d0)) then
       h_shoc = &
            dabs(j*grad_N(1)) + &
            dabs(j*grad_N(2)) + &
            dabs(j*grad_N(3))
       h_shoc = 2/h_shoc
       tau_shoc = h_shoc*0.5/v_ref*(norm(grad_rho)*h_shoc/rho_ref)**beta
       nu_shoc = tau_shoc*(v_ref)**2
    else
       nu_shoc = 0.d0
    end if
  end subroutine compute_stability

  subroutine compute_mass()
    ! real(kind=8), dimension(:), intent(out) :: mass
    ! type(MeshType), intent(in) :: mesh
    integer :: ipoi1, ipoi2, ipoi3, ielem
    real(kind=8) :: area
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mass(:) = 0.d0
    do ielem = 1, mesh%nelem
       area = mesh%elements(ielem)%area
       ipoi1 = mesh%elements(ielem)%points(1)
       ipoi2 = mesh%elements(ielem)%points(2)
       ipoi3 = mesh%elements(ielem)%points(3)
       mass(ipoi1) = mass(ipoi1) + area/3.d0
       mass(ipoi2) = mass(ipoi2) + area/3.d0
       mass(ipoi3) = mass(ipoi3) + area/3.d0
    end do
  end subroutine compute_mass

end module solver_mod
