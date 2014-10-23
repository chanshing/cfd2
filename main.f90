program main
  use mesh_mod
  use general_mod
  use IO_mod
  ! use solver_mod
  use iter_solver_mod
  use runge_kutta_mod
  implicit none
  type(SolverParamsType) :: solverParams
  type(ReferenceValuesType) :: refVal
  type(BoundaryConditionsType) :: bCond
  type(MeshType) :: mesh
  type(PrimitiveVarsType), dimension(:), allocatable :: pVars
  real(kind=8), dimension(:,:), allocatable :: U
  real(kind=8) :: time, dt
  real(kind=8) :: progress_time
  integer :: start_t, end_t, rate
  integer :: ipoin, iter, niter

  call IOread(&
       refVal,&
       solverParams,&
       mesh,&
       bcond)

  allocate(U(4,mesh%npoin))
  allocate(pVars(mesh%npoin))

  call initCond(U, bCond, mesh, refVal)

  dt = 2d-7*solverParams%timestep_factor

  ! call compute_rhs_config(mesh, refVal, dt)
  call compute_config(mesh, refVal, dt)
  call iterator_config(U)

  ! call rk_config(U)

  open(101, file='inv_tau1.mat', status='unknown')
  open(102, file='inv_tau3_v.mat', status='unknown')
  open(103, file='inv_tau3_e.mat', status='unknown')
  open(104, file='nu_shoc.mat', status='unknown')
  write(101,*) "# name: inv_tau1"
  write(101,*) "# type: matrix"
  write(101,*) "# rows: 1000"
  write(101,*) "# columns: 2"
  write(102,*) "# name: inv_tau3_v"
  write(102,*) "# type: matrix"
  write(102,*) "# rows: 1000"
  write(102,*) "# columns: 2"
  write(103,*) "# name: inv_tau3_e"
  write(103,*) "# type: matrix"
  write(103,*) "# rows: 1000"
  write(103,*) "# columns: 2"
  write(104,*) "# name: nu_shoc"
  write(104,*) "# type: matrix"
  write(104,*) "# rows: 1000"
  write(104,*) "# columns: 2"

  time = 0
  progress_time = 0
  niter = solverParams%max_iterations
  do iter = 1, niter
     call system_clock(start_t, rate)
     time = time + dt
     ! call rk4(U, compute_rhs, dt)
     call iterator(U, dt)
     call restoreBC(U, bCond, mesh, refVal)
     call system_clock(end_t)
     progress_time = progress_time + dble(end_t - start_t)/dble(rate)
     call printProgress(iter, niter, progress_time, time)
  end do

  call U_to_prim(pVars, U, refVal)
  open(2, file=trim(filename)//'.flavia.res', status='unknown')

  !CCCC----> ESCRITURA DE VELOCIDADES
  !CCCC------------------------------
  write(2, '(A15, 5(I8, 2X))') 'VELOCITY', 2, iter, 2, 1, 1
  write(2, '(A)') 'VEL_X'
  write(2, '(A)') 'VEL_Y'
  do ipoin = 1, mesh%npoin
     write(2, '(I8, 3E13.4)') ipoin, pVars(ipoin)%velocity%x, pVars(ipoin)%velocity%y
  end do

  !CCCC----> ESCRITURA DE LAS DENSIDADES
  !CCCC---------------------------------
  write(2, '(A15, 5(I8, 2X))') 'DENSITY', 2, ITER, 1, 1, 1
  write(2, '(A)') 'DENSITY'
  do ipoin = 1, mesh%npoin
     write(2, '(I8, 1E13.4)') ipoin, pVars(ipoin)%density
  end do

  close(2)

end program main
