module iter_solver_mod
  use mesh_mod
  use general_mod
  implicit none
  private
  save

  public :: &
       iterator,&
       iterator_config,&
       compute_config
  
  real(kind=8), dimension(:), allocatable :: mass
  type(PrimitiveVarsType), dimension(:), allocatable :: pVars
  real(kind=8), dimension(:,:), allocatable :: dU_old, dU, r, U_old, l, U_theta, dU_theta

  ! ARGUMENTOS
  type(MeshType), pointer :: mesh
  type(ReferenceValuesType), pointer :: refVal
  real(kind=8), pointer :: dt
contains

  subroutine iterator_config(U)
    real(kind=8), dimension(:,:), intent(in) :: U
    integer :: m, n
    m = size(U,1)
    n = size(U,2)
    if(allocated(dU_old)) deallocate(dU_old)
    allocate(dU_old(m,n))
    if(allocated(dU)) deallocate(dU)
    allocate(dU(m,n))
    if(allocated(dU_theta)) deallocate(dU_theta)
    allocate(dU_theta(m,n))
    if(allocated(r)) deallocate(r)
    allocate(r(m,n))
    if(allocated(l)) deallocate(l)
    allocate(l(m,n))
    if(allocated(U_old)) deallocate(U_old)
    allocate(U_old(m,n))
    if(allocated(U_theta)) deallocate(U_theta)
    allocate(U_theta(m,n))
  end subroutine iterator_config

  subroutine iterator(U, dt)
    real(kind=8), dimension(:,:), intent(inout) :: U
    real(kind=8), intent(in) :: dt
    real(kind=8), parameter :: theta = 0.5
    integer :: i, ipoin

    call compute_r(r,U)
    do ipoin = 1, mesh%npoin
       dU(:,ipoin) = dt*r(:,ipoin)/mass(ipoin)
    end do

    do i = 1, 3
       dU_theta = (1-theta)*dU_old + theta*dU
       call compute_l(l, dU_theta, U)
       do ipoin = 1, mesh%npoin
          dU(:,ipoin) = dU_old(:,ipoin) + (dt*r(:,ipoin) - l(:,ipoin))/mass(ipoin)
       end do
    end do

    U = U + dU

    ! call compute_r(r,U_old)
    ! do ipoin = 1, mesh%npoin
    !    U(:,ipoin) = U_old(:,ipoin) + dt*r(:,ipoin)/mass(ipoin)
    ! end do

    ! do i = 1, 3
    !    U_theta = (1 - theta)*U_old + theta*U
    !    call compute_r(r,U_theta)
    !    do ipoin = 1, mesh%npoin
    !       U(:,ipoin) = U_old(:,ipoin) + dt*r(:,ipoin)/mass(ipoin)
    !    end do
    ! end do

  end subroutine iterator

  subroutine compute_config(mesh0, refVal0, dt0)
    type(MeshType), intent(in), target :: mesh0
    type(ReferenceValuesType), intent(in), target :: refVal0
    real(kind=8), intent(in), target :: dt0
    mesh => mesh0
    refVal => refVal0
    dt => dt0
    if(.not.allocated(pVars)) allocate(pVars(mesh%npoin))
    if(.not.allocated(mass)) allocate(mass(mesh%npoin))
    call compute_mass()
  end subroutine compute_config

  subroutine compute_l(l, dU, U)
    real(kind=8), dimension(:,:), intent(out) :: l
    real(kind=8), dimension(:,:), intent(in) :: dU
    real(kind=8), dimension(:,:), intent(in) :: U
    real(kind=8) :: U_local(4)
    real(kind=8) :: dU_local(4)
    real(kind=8) :: l_local(4,3)
    real(kind=8), dimension(3) :: Nx, Ny
    real(kind=8) :: rho, v1, v2, e, w1, w2
    real(kind=8) :: tau_r, tau_v, tau_e, nu_shoc
    real(kind=8) :: area, gamma
    real(kind=8) :: mu, T, k
    real(kind=8) :: tauA1dU(4), tauA2dU(4)
    integer :: i, ipoi1, ipoi2, ipoi3, ielem
    real(kind=8), parameter, dimension(3,3) :: N = &
         reshape((/&
         0.d0, .5d0, .5d0, &
         .5d0, 0.d0, .5d0, &
         .5d0, .5d0, 0.d0 &
         /), (/3,3/))

    l(:,:) = 0.d0
    gamma = refVal%heat_capacity_ratio
    do ielem = 1, mesh%nelem
       l_local(:,:) = 0.d0
       area = mesh%elements(ielem)%area
       ipoi1 = mesh%elements(ielem)%points(1)
       ipoi2 = mesh%elements(ielem)%points(2)
       ipoi3 = mesh%elements(ielem)%points(3)
       Nx(1:3) = mesh%elements(ielem)%dshape(:)%x
       Ny(1:3) = mesh%elements(ielem)%dshape(:)%y
       T = (pVars(ipoi1)%temperature + &
            pVars(ipoi2)%temperature + &
            pVars(ipoi3)%temperature)/3.d0
       mu = get_mu(refVal, T)
       k = get_k(refVal, T)

       call compute_stability(tau_r, tau_v, tau_e, nu_shoc, &
            ielem, mesh, pVars, refVal, mu, dt)

       do i = 1, 3
          ! U_local(1:4) = &
          !      N(1,i)*U(1:4,ipoi1) +&
          !      N(2,i)*U(1:4,ipoi2) +&
          !      N(3,i)*U(1:4,ipoi3)
          U_local(1:4) = (&
               U(1:4,ipoi1) +&
               U(1:4,ipoi2) +&
               U(1:4,ipoi3))/3.d0
          dU_local(1:4) = &
               N(1,i)*dU(1:4,ipoi1) +&
               N(2,i)*dU(1:4,ipoi2) +&
               N(3,i)*dU(1:4,ipoi3)
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
          tauA1dU(1) = tau_r*(-w1*dU_local(1) + dU_local(2))
          tauA1dU(2) = (1.0d0/2.0d0)*tau_v*(-dU_local(1)*(2*v1**2 - (gamma - 1)*(v1**2 &
               + v2**2)) - 2*dU_local(2)*(w1 + v1*(gamma - 3)) - 2*dU_local(3)*v2*(gamma - 1) + &
               2*dU_local(4)*(gamma - 1))
          tauA1dU(3) = tau_v*(-dU_local(1)*v1*v2 + dU_local(2)*v2 - dU_local(3)*(w1 - v1))
          tauA1dU(4) = -1.0d0/2.0d0*tau_e*(2*dU_local(1)*v1*(e*gamma - (gamma - 1)*(v1 &
               **2 + v2**2)) + dU_local(2)*(-2*e*gamma + 2*v1**2*(gamma - 1) + (gamma - &
               1)*(v1**2 + v2**2)) + 2*dU_local(3)*v1*v2*(gamma - 1) + 2*dU_local(4)*(w1 - gamma &
               *v1))
          tauA2dU(1) = tau_r*(-w2*dU_local(1) + dU_local(3))
          tauA2dU(2) = tau_v*(-dU_local(1)*v1*v2 - dU_local(2)*(w2 - v2) + dU_local(3)*v1)
          tauA2dU(3) = (1.0d0/2.0d0)*tau_v*(-dU_local(1)*(2*v2**2 - (gamma - 1)*(v1**2 &
               + v2**2)) - 2*dU_local(2)*v1*(gamma - 1) - 2*dU_local(3)*(w2 + v2*(gamma - 3)) + &
               2*dU_local(4)*(gamma - 1))
          tauA2dU(4) = -1.0d0/2.0d0*tau_e*(2*dU_local(1)*v2*(e*gamma - (gamma - 1)*(v1 &
               **2 + v2**2)) + 2*dU_local(2)*v1*v2*(gamma - 1) + dU_local(3)*(-2*e*gamma + 2*v2 &
               **2*(gamma - 1) + (gamma - 1)*(v1**2 + v2**2)) + 2*dU_local(4)*(w2 - &
               gamma*v2))

          l_local(:,1) = l_local(:,1) + N(1,i)*dU_local + (Nx(1)*tauA1dU + Ny(1)*tauA2dU) 
          l_local(:,2) = l_local(:,2) + N(2,i)*dU_local + (Nx(2)*tauA1dU + Ny(2)*tauA2dU) 
          l_local(:,3) = l_local(:,3) + N(3,i)*dU_local + (Nx(3)*tauA1dU + Ny(3)*tauA2dU) 

       end do
       l_local = l_local*area/3.d0

       do i = 1, 4
          l(i,ipoi1) = l(i,ipoi1) + l_local(i,1)
          l(i,ipoi2) = l(i,ipoi2) + l_local(i,2)
          l(i,ipoi3) = l(i,ipoi3) + l_local(i,3)
       end do
    end do

  end subroutine compute_l
  
  subroutine compute_r(r, U)
    real(kind=8), dimension(:,:), intent(out) :: r
    real(kind=8), dimension(:,:), intent(in) :: U
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    real(kind=8), dimension(4,3) :: r_local
    real(kind=8), dimension(4) :: U_local
    real(kind=8) :: rho, v1, v2, e, w1, w2
    real(kind=8) :: area
    real(kind=8) :: tau_r, tau_v, tau_e, nu_shoc
    real(kind=8) :: gamma, Cv, k, mu, T
    integer :: ielem, ipoi1, ipoi2, ipoi3, i
    real(kind=8), dimension(3) :: Nx, Ny
    real(kind=8), dimension(4) :: Ux, Uy
    real(kind=8), dimension(4) :: AiUi
    real(kind=8), dimension(4) :: A1AiUi
    real(kind=8), dimension(4) :: A2AiUi
    real(kind=8), parameter, dimension(3,3) :: N = &
         reshape((/&
         0.d0, .5d0, .5d0, &
         .5d0, 0.d0, .5d0, &
         .5d0, .5d0, 0.d0 &
         /), (/3,3/))

    call U_to_prim(pVars, U, refVal)

    r(:,:) = 0.d0
    gamma = refVal%heat_capacity_ratio
    Cv = refVal%specific_heat_capacity
    do ielem = 1, mesh%nelem
       r_local(:,:) = 0.d0
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
          ! U_local(1:4) = &
          !      N(1,i)*U(1:4,ipoi1) +&
          !      N(2,i)*U(1:4,ipoi2) +&
          !      N(3,i)*U(1:4,ipoi3)
          U_local(1:4) = (&
               U(1:4,ipoi1) +&
               U(1:4,ipoi2) +&
               U(1:4,ipoi3))/3.d0
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

          r_local(:,1) = r_local(:,1) &
               + N(1,i)*AiUi &
               + (Nx(1)*A1AiUi + Ny(1)*A2AiUi) &
               + nu_shoc*(Nx(1)*Ux + Ny(1)*Uy) &
               - N(1,i)*(w1*Ux + w2*Uy)
          r_local(:,2) = r_local(:,2) &
               + N(2,i)*AiUi &
               + (Nx(2)*A1AiUi + Ny(2)*A2AiUi) &
               + nu_shoc*(Nx(2)*Ux + Ny(2)*Uy) &
               - N(2,i)*(w1*Ux + w2*Uy)
          r_local(:,3) = r_local(:,3) &
               + N(3,i)*AiUi &
               + (Nx(3)*A1AiUi + Ny(3)*A2AiUi) &
               + nu_shoc*(Nx(3)*Ux + Ny(3)*Uy) &
               - N(3,i)*(w1*Ux + w2*Uy)

          ! call rhs_diffusive(r_local, &
          !      Nx, Ny, Ux, Uy, &
          !      rho, v1, v2, e, &
          !      Cv, k, mu)

       end do

       r_local = r_local*area/3.d0

       do i = 1, 4
          r(i,ipoi1) = r(i,ipoi1) + r_local(i,1)
          r(i,ipoi2) = r(i,ipoi2) + r_local(i,2)
          r(i,ipoi3) = r(i,ipoi3) + r_local(i,3)
       end do
       
    end do

    r = -r

  end subroutine compute_r

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

  subroutine compute_stability(tau_r, tau_v, tau_e, nu_shoc, &
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
    real(kind=8) :: inv_tau1, inv_tau2, inv_tau3_v, inv_tau3_e
    real(kind=8) :: tau_shoc1, tau_shoc2, tau_shoc
    real(kind=8) :: h_rgn, h_rgn_e, h_shoc
    real(kind=8), parameter :: beta = 1
    real(kind=8), parameter :: tol = 1d-10
    real(kind=8) :: dummy
    integer, save :: counter = 0
    logical :: doit
    integer :: elem
    type(Vector) :: v_el
    real(kind=8) :: c_el
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elem = 19910
    counter = counter + 1
    if(mod(counter,4) == 0) then
       doit = .true.
    else
       doit = .false.
    end if
    points(:) = mesh%elements(ielem)%points(:)
    rho(:) = pVars(points(:))%density
    T(:) = pVars(points(:))%temperature
    c(:) = pVars(points(:))%sound_speed
    v(:) = pVars(points(:))%velocity - mesh%vel(points(:))
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
    dummy = norm(grad_rho)
    ! if(norm(grad_rho) > tol) then
    !    j = grad_rho/norm(grad_rho)
    !    inv_tau1 = &
    !         c(1)*dabs(j*grad_N(1)) + dabs(v(1)*grad_N(1)) + &
    !         c(2)*dabs(j*grad_N(2)) + dabs(v(2)*grad_N(2)) + &
    !         c(3)*dabs(j*grad_N(3)) + dabs(v(3)*grad_N(3))
    ! else
    !    inv_tau1 = &
    !         dabs(v(1)*grad_N(1)) + &
    !         dabs(v(2)*grad_N(2)) + &
    !         dabs(v(3)*grad_N(3))
    ! end if
    v_el = (v(1) + v(2) + v(3))/3.d0
    c_el = (c(1) + c(2) + c(3))/3.d0
    if(norm(grad_rho) < tol) grad_rho = Vector(1.d0,0.d0)
    j = grad_rho/norm(grad_rho)
    inv_tau1 = &
         c_el*dabs(j*grad_N(1)) + dabs(v_el*grad_N(1)) + &
         c_el*dabs(j*grad_N(2)) + dabs(v_el*grad_N(2)) + &
         c_el*dabs(j*grad_N(3)) + dabs(v_el*grad_N(3))
    if(ielem == elem) then
       if(doit) then
          ! print*, grad_rho
          ! print*, rho(:)
          write(101,*) dummy, inv_tau1
       end if
    end if

    ! 1/TAU2
    inv_tau2 = 2/dt

    ! 1/TAU3_V
    r = &
         norm(v(1))*grad_N(1) + &
         norm(v(2))*grad_N(2) + &
         norm(v(3))*grad_N(3)
    dummy = norm(r)
    if(norm(r) > tiny(0d0)) then
       r = r/norm(r)
       h_rgn = &
            dabs(r*grad_N(1)) + &
            dabs(r*grad_N(2)) + &
            dabs(r*grad_N(3))
       h_rgn = 2.d0/h_rgn
       inv_tau3_v = 1/(0.25d0*h_rgn**2/mu)
    else
       inv_tau3_v = 0.d0
    end if
    if(ielem == elem) then
       if(doit) then
          write(102,*) dummy, inv_tau3_v
       end if
    end if
    ! %%%%%%%%%% DEBUG %%%%%%%%%%
    inv_tau3_v = 0.d0

    ! 1/TAU3_E
    r_e = &
         T(1)*grad_N(1) + &
         T(2)*grad_N(2) + &
         T(3)*grad_N(3)
    dummy = norm(r_e)
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
    if(norm(r_e) < tiny(0d0)) r_e = Vector(1.d0,1.d0)
    r_e = r_e/norm(r_e)
    h_rgn_e = &
         dabs(r_e*grad_N(1)) + &
         dabs(r_e*grad_N(2)) + &
         dabs(r_e*grad_N(3))
    h_rgn_e = 2.d0/h_rgn_e
    inv_tau3_e = 1/(0.25d0*h_rgn_e**2/mu)
    if(ielem == elem) then
       if(doit) then
          write(103,*) dummy, inv_tau3_e
       end if
    end if
    ! %%%%%%%%%% DEBUG %%%%%%%%%%
    inv_tau3_e = 0

    tau_r = (1/(inv_tau1**2 + inv_tau2**2))**0.5
    tau_v = (1/(inv_tau1**2 + inv_tau2**2 + inv_tau3_v**2))**0.5
    tau_e = (1/(inv_tau1**2 + inv_tau2**2 + inv_tau3_e**2))**0.5

    ! NU_SHOC
    dummy = norm(grad_rho)
    if(norm(grad_rho) > tol) then
       h_shoc = &
            dabs(j*grad_N(1)) + &
            dabs(j*grad_N(2)) + &
            dabs(j*grad_N(3))
       h_shoc = 2/h_shoc
       tau_shoc1 = h_shoc*0.5/v_ref*(norm(grad_rho)*h_shoc/rho_ref)
       tau_shoc2 = h_shoc*0.5/v_ref*(norm(grad_rho)*h_shoc/rho_ref)**2
       tau_shoc = 0.5*(tau_shoc1 + tau_shoc2)
       nu_shoc = tau_shoc*(v_ref)**2
    else
       nu_shoc = 0.d0
    end if
    if(ielem == elem) then
       if(doit) then
          write(104,*) dummy, nu_shoc
       end if
    end if
  end subroutine compute_stability

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
end module iter_solver_mod
