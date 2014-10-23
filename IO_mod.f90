module IO_mod
  use mesh_mod
  use general_mod
  implicit none
  private
  public :: IOread, &
       IOParamsType, &
       printProgress, &
       filename

  type IOParamsType
     integer :: isRestart
     integer :: PRINT_INTERVAL
     integer :: output_movie
     character(len=4) :: output_density
     character(len=4) :: output_velocity
     character(len=4) :: output_mach
     character(len=4) :: output_pressure
     character(len=4) :: output_temperature
     character(len=4) :: output_energy
     character(len=4) :: output_position
  end type IOParamsType

  type(IOParamsType) :: IOParams

  character(len=80) :: filename
  
contains

  subroutine IOread(&
       refVal,&
       solverParams,&
       mesh,&
       bCond)
    type(ReferenceValuesType) :: refVal
    type(SolverParamsType) :: solverParams
    type(MeshType) :: mesh
    type(BoundaryConditionsType) :: bCond

    call readInputData(&
         refVal,&
         solverParams)
    call readMeshData(&
         mesh, &
         bCond, &
         refVal)

  end subroutine IOread

  subroutine readInputData(&
       refVal,&
       solverParams)
    type(ReferenceValuesType) :: refVal
    type(SolverParamsType) :: solverParams
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    type(Vector) :: V
    real(kind=8) :: T, rho, mach, P, c, gamma, R
    real(kind=8) :: dummy
    open(1, file='EULER.DAT', status='old')
    read(1, '(A)') filename
    close(1)
    open(1, file=trim(filename)//'-1.dat', status='old')

    read(1,*)
    read(1,*) IOParams%isRestart, &
         solverParams%max_iterations, &
         IOParams%print_interval, &
         solverParams%isLocalTimeStep
    read(1,*)
    read(1,*) solverParams%timestep_factor, &
         refVal%velocity%x, &
         refVal%velocity%y, &
         refVal%mach, &
         refVal%temperature, &
         refVal%density, &
         refVal%pressure
    read(1,*)
    read(1,*) refVal%viscosity, &
         refVal%gravity%x, &
         refVal%gravity%y, &
         dummy
    read(1,*)
    read(1,*) refVal%thermal_conductivity, &
         refVal%gas_constant, &
         refVal%specific_heat_capacity, &
         refVal%heat_capacity_ratio, &
         dummy
    read(1,*)
    read(1,*) solverParams%shock_factor
    read(1,*)
    read(1,*)
    read(1,*) solverParams%isMoving, dummy, dummy
    read(1,*)
    read(1,*)
    read(1,*) IOParams%output_density, &
         IOParams%output_velocity, &
         IOParams%output_mach, &
         IOParams%output_pressure, &
         IOParams%output_temperature, &
         IOParams%output_energy, &
         IOParams%output_position
    close(1)

    ! DEDUCIR VALORES FALTANTES
    solverParams%shock_factor = 1.d0/solverParams%shock_factor
    V = refVal%velocity
    T = refVal%temperature
    rho = refVal%density
    mach = refVal%mach
    P = refVal%pressure
    gamma = refVal%heat_capacity_ratio
    R = refVal%gas_constant

    if(T < tiny(0.d0)) then
       T = P/(R*rho)
       refVal%temperature = T
    end if
    if(P < tiny(0.d0)) then
       P = R*rho*T
       refVal%pressure = P
    end if
    if(rho < tiny(0.d0)) then
       rho = P/(R*T)
       refVal%density = rho
    end if
    c = dsqrt(gamma*R*T)
    refVal%sound_speed = c
    if(norm(V) < tiny(0.d0)) then
       V%y = 0.d0
       V%x = c*mach
       refVal%velocity = V
    end if

    refVal%total_temperature = &
         T*(1 + 0.5*(gamma - 1)*mach**2)
    refVal%total_pressure = &
         P*(1 + 0.5*(gamma - 1)*mach**2)**(gamma/(gamma - 1))
    refVal%total_density = &
         rho*(1 + 0.5*(gamma - 1)*mach**2)**(1/(gamma - 1))

  end subroutine readInputData

  subroutine readMeshData(mesh, bCond, refVal)
    use bodies_mod
    implicit integer(i,j,k,l,m,n)
    type(MeshType) :: mesh
    type(BoundaryConditionsType) :: bCond
    type(ReferenceValuesType) :: refVal
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    type(IndexedEdge) :: iedge
    real(kind=8) :: factor, factorX, factorY

    open(1, file=trim(filename)//'.dat', status='old')
    read(1,*) ! NUMERO DE NODOS, NUMERO DE ELEMENTOS
    read(1,*) npoin, nelem
    read(1,*) 
    read(1,*) nfixrho, nfixv, nnoslip, nwall, &
         nfixt, nsets, nmaster, nslave, &
         nfixpoint, nmove
    read(1,*)
    read(1,*)
    read(1,*)

    mesh%npoin = npoin
    mesh%nelem = nelem
    allocate(mesh%coord(npoin), &
         mesh%elements(nelem), &
         mesh%vel(npoin), &
         mesh%fixedPoints(nfixpoint))
    allocate(bCond%velocities(nfixv), &
         bCond%temperatures(nfixt), &
         bCond%densities(nfixrho))

    !CCCC-----> COORDENADAS DE LOS NODOS
    read(1,*)
    do ipoin = 1, npoin
       read(1,*) i, &
            mesh%coord(i)%x, &
            mesh%coord(i)%y
    end do

    !CCCC----> CONECTIVIDADES DE LOS ELEMENTOS
    read(1,*)
    do ielem = 1, nelem
       read(1,*) i, &
            mesh%elements(i)%points(1), &
            mesh%elements(i)%points(2), &
            mesh%elements(i)%points(3)
    end do

    !CCCC----> PUNTOS CON DENSIDAD PRESCRITA
    read(1,*)
    do i = 1, nfixrho
       read(1,*) &
            bCond%densities(i)%point, &
            factor
       bCond%densities(i)%val = factor*refVal%density
    end do

    !CCCC----> NODOS CON VELOCIDAD FIJA
    read(1,*)
    do i = 1, nfixv
       read(1,*) &
            bCond%velocities(i)%point, &
            factorX, &
            factorY
       bCond%velocities(i)%x = factorX*refVal%velocity%x
       bCond%velocities(i)%y = factorY*refVal%velocity%y
    end do
    
    !CCCC----> NODOS CON NO SLIP
    read(1,*)
    do i = 1, nnoslip
       read(1,*)
    end do

    !CCCC----> PAREDES CON VELOCIDAD NORMAL NULA
    read(1,*)
    do i = 1, nwall
       read(1,*) iedge%points(1), iedge%points(2)
        call bCond%normalVelocities%addEdge(iedge)
    end do
    call bCond%normalVelocities%updateNormals(mesh%coord)

    !CCCC----> NODOS CON TEMPERATURA PRESCRITA
    read(1,*)
    do i = 1, nfixt
       read(1,*) &
            bCond%temperatures(i)%point, &
            factor
       bCond%temperatures(i)%val = factor*refVal%temperature
    end do

    !CCCC----> LECTURA DE LOS SETS 
    call init_bodies()
    read(1,*)
    do i = 1, nsets
       read(1,*) &
            iedge%element, &
            iedge%points(1), &
            iedge%points(2), &
            ibody
       call bodies(ibody)%obj%addEdge(iedge)
    end do
    allocate(mesh%bodies(size(bodies)))
    mesh%bodies = bodies
    
    !CCCC----> PERIODICAL MASTER
    read(1,*)
    do i = 1, nmaster
       read(1,*)
    end do
    read(1,*)
    do i = 1, nmaster
       read(1,*)
    end do

    !CCCC----> NODOS FIJOS
    read(1,*)
    do i = 1, nfixpoint
       read(1,*) mesh%fixedPoints(i), idummy
    end do

    !CCCC----> NODOS CON MOVIMIENTO
    read(1,*)
    do i = 1, nmove
       read(1,*)
    end do
    close(1)
    
    call mesh%update()
    mesh%vel(:)%x = 0.d0
    mesh%vel(:)%y = 0.d0

  end subroutine readMeshData

  subroutine printProgress(iter, niter, progress_time, time)
    integer, intent(in) :: iter, niter
    real(kind=8), intent(in) :: progress_time
    real(kind=8), intent(in) :: time
    integer :: time_sec
    real(kind=8) :: time_ms
    integer :: estimated_time, hr, min, sec
    character(len=30) :: Format0, Format1, Format2
    Format0 = "(A, I7.2, A, I7.2)"
    Format1 = "(A, I2, A, F7.2, A)"
    Format2 = "(A, I3.2, A, I3.2, A, I3.2, A)"
    if(mod(iter,IOParams%PRINT_INTERVAL) == 0) then
       write(*,Format0) "Iteraciones:", iter, "/", niter

       time_sec = int(time)
       time_ms = (time - time_sec)*1000
       write(*,Format1) "t =", time_sec, "s", time_ms, "ms"

       estimated_time = int(estimatedTime(iter, niter, progress_time))
       hr = int(estimated_time/3600)
       min = int(estimated_time/60)
       sec = mod(int(estimated_time),60)
       write(*,Format2) "Tiempo restante:", &
            hr, "h", min, "m", sec, "s"
    end if
  end subroutine printProgress

  pure real(kind=8) function estimatedTime&
       (iter, niter, progress_time)
    integer, intent(in) :: iter, niter
    real(kind=8), intent(in) :: progress_time
    estimatedTime = dble(niter - iter)/dble(iter)*progress_time
  end function estimatedTime
    
end module IO_mod
