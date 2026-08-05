program fortran_reference_benchmarks
  implicit none

  call diffusion_case()
  call river_and_deposition_case()
  call glacial_case()
  call rigid_rotation_case()

contains

  subroutine diffusion_case()
    integer, parameter :: nx=65, ny=65, nn=nx*ny, nstep=10
    integer :: i, j, node, step
    double precision, parameter :: length=100.d3, dt=1.d5, diffusivity=100.d0
    double precision :: x, y, pi
    double precision :: h(nn), basement(nn), area(nn), erosion(nn)
    double precision :: uplift(nn), river(nn), hillslope(nn)

    pi = acos(-1.d0)
    call FastScape_Init()
    call FastScape_Set_NX_NY(nx,ny)
    call FastScape_Setup()
    call FastScape_Set_XL_YL(length,length)
    call FastScape_Set_DT(dt)
    call FastScape_Set_BC(1111)

    do j=1,ny
      y = dble(j-1)*length/dble(ny-1)
      do i=1,nx
        x = dble(i-1)*length/dble(nx-1)
        node = i+(j-1)*nx
        h(node) = 1000.d0*sin(pi*x/length)*sin(pi*y/length)
      enddo
    enddo
    uplift = 0.d0
    river = 0.d0
    hillslope = diffusivity
    call FastScape_Init_H(h)
    call FastScape_Set_U(uplift)
    call FastScape_Set_Erosional_Parameters(river,-1.d0,0.4d0,1.d0, &
                                             hillslope,-1.d0,0.d0,0.d0,-2.d0)
    basement = h
    area = 0.d0
    erosion = 0.d0
    call write_grid('output/fortran_diffusion_initial.csv',nx,ny,length,length, &
                    h,basement,area,erosion)
    do step=1,nstep
      call FastScape_Execute_Step()
    enddo
    call FastScape_Copy_H(h)
    call FastScape_Copy_Basement(basement)
    call FastScape_Copy_Drainage_Area(area)
    call FastScape_Copy_Erosion_Rate(erosion)
    call write_grid('output/fortran_diffusion_final.csv',nx,ny,length,length, &
                    h,basement,area,erosion)
    call FastScape_Destroy()
  end subroutine diffusion_case

  subroutine river_and_deposition_case()
    integer, parameter :: nx=81, ny=65, nn=nx*ny, nstep=80
    integer :: i, j, node, step
    double precision, parameter :: xlength=100.d3, ylength=80.d3, dt=500.d0
    double precision :: x, y, pi
    double precision :: h(nn), basement(nn), area(nn), erosion(nn)
    double precision :: uplift(nn), river(nn), hillslope(nn)

    pi = acos(-1.d0)
    call FastScape_Init()
    call FastScape_Set_NX_NY(nx,ny)
    call FastScape_Setup()
    call FastScape_Set_XL_YL(xlength,ylength)
    call FastScape_Set_DT(dt)
    call FastScape_Set_BC(1000)
    do j=1,ny
      y = dble(j-1)*ylength/dble(ny-1)
      do i=1,nx
        x = dble(i-1)*xlength/dble(nx-1)
        node = i+(j-1)*nx
        h(node) = 80.d0 + 900.d0/(1.d0+exp(-(y-0.52d0*ylength)/2500.d0)) &
                  + 12.d0*sin(2.d0*pi*x/17.d3)*sin(2.d0*pi*y/13.d3)
      enddo
    enddo
    uplift = 0.d0
    river = 7.d-6
    hillslope = 5.d-2
    call FastScape_Init_H(h)
    call FastScape_Set_U(uplift)
    call FastScape_Set_Erosional_Parameters(river,1.d-5,0.45d0,1.d0, &
                                             hillslope,7.5d-2,1.d0,1.d0,1.d0)
    basement = h
    area = 0.d0
    erosion = 0.d0
    call write_grid('output/fortran_river_deposition_initial.csv',nx,ny, &
                    xlength,ylength,h,basement,area,erosion)
    do step=1,nstep
      call FastScape_Execute_Step()
    enddo
    call FastScape_Copy_H(h)
    call FastScape_Copy_Basement(basement)
    call FastScape_Copy_Drainage_Area(area)
    call FastScape_Copy_Erosion_Rate(erosion)
    call write_grid('output/fortran_river_deposition_final.csv',nx,ny, &
                    xlength,ylength,h,basement,area,erosion)
    call FastScape_Destroy()
  end subroutine river_and_deposition_case

  subroutine glacial_case()
    integer, parameter :: nx=65, ny=65, nn=nx*ny
    integer :: i, j, node
    double precision, parameter :: length=100.d3, dt=500.d0
    double precision :: x, y, shape
    double precision :: h(nn), initial_h(nn), basement(nn), area(nn), erosion(nn)
    double precision :: ice(nn), sliding(nn), glacial_coefficient(nn)

    call FastScape_Init()
    call FastScape_Set_NX_NY(nx,ny)
    call FastScape_Setup()
    call FastScape_Set_XL_YL(length,length)
    call FastScape_Set_DT(dt)
    call FastScape_Set_BC(1111)
    do j=1,ny
      y = dble(j-1)*length/dble(ny-1)
      do i=1,nx
        x = dble(i-1)*length/dble(nx-1)
        node = i+(j-1)*nx
        shape = exp(-((x-0.5d0*length)/22.d3)**6 &
                    -((y-0.5d0*length)/15.d3)**6)
        h(node) = 1200.d0 + 400.d0*(x/length) + 100.d0*abs(y-0.5d0*length)/length
        ice(node) = 1000.d0*shape
        sliding(node) = 100.d0*shape
        glacial_coefficient(node) = 1.d-4
      enddo
    enddo
    ice(1:nx)=0.d0
    ice(nn-nx+1:nn)=0.d0
    ice(1:nn:nx)=0.d0
    ice(nx:nn:nx)=0.d0
    initial_h = h
    call FastScape_Init_H(h)
    call FastScape_Set_Glacial_Parameters(glacial_coefficient,1.d0,1.d0,10.d0)
    call FastScape_Set_Ice_Thickness(ice)
    call FastScape_Set_Basal_Sliding_Velocity(sliding)
    basement = h
    area = ice
    erosion = sliding
    call write_grid('output/fortran_glacial_initial.csv',nx,ny,length,length, &
                    h,basement,area,erosion)
    call FastScape_Execute_Step()
    call FastScape_Copy_H(h)
    call FastScape_Copy_Basement(basement)
    call FastScape_Copy_Drainage_Area(area)
    call FastScape_Copy_Glacial_Erosion_Rate(erosion)
    call write_grid('output/fortran_glacial_final.csv',nx,ny,length,length, &
                    h,basement,area,erosion)
    call FastScape_Destroy()
  end subroutine glacial_case

  subroutine rigid_rotation_case()
    integer, parameter :: nx=51, ny=51, nn=nx*ny
    integer :: i, j, node, step, nstep
    double precision, parameter :: length=100.d3, period=1.d6
    double precision :: x, y, dx, radius, omega, maximum_speed, dt, pi
    double precision :: h(nn), basement(nn), area(nn), erosion(nn)
    double precision :: uplift(nn), river(nn), hillslope(nn), vx(nn), vy(nn)

    pi = acos(-1.d0)
    dx = length/dble(nx-1)
    radius = sqrt(2.d0)*0.5d0*length
    omega = 2.d0*pi/period
    maximum_speed = omega*radius
    ! Advect_TVD enforces a directional Courant number of 0.2 internally.
    ! Stay below that value so Execute_Step does not silently shorten the
    ! modeled interval while this driver still counts one requested step.
    nstep = ceiling(period/(0.18d0*dx/maximum_speed))
    dt = period/dble(nstep)
    call FastScape_Init()
    call FastScape_Set_Advection_Scheme(2)
    call FastScape_Set_NX_NY(nx,ny)
    call FastScape_Setup()
    call FastScape_Set_XL_YL(length,length)
    call FastScape_Set_DT(dt)
    call FastScape_Set_BC(1111)
    do j=1,ny
      y = dble(j-1)*length/dble(ny-1)
      do i=1,nx
        x = dble(i-1)*length/dble(nx-1)
        node = i+(j-1)*nx
        h(node) = 900.d0*exp(-((x-0.62d0*length)/11.d3)**2 &
                              -((y-0.48d0*length)/8.d3)**2) &
                  +300.d0*exp(-((x-0.42d0*length)/7.d3)**2 &
                              -((y-0.64d0*length)/13.d3)**2)
        vx(node) = -omega*(y-0.5d0*length)
        vy(node) =  omega*(x-0.5d0*length)
      enddo
    enddo
    vx(1:nx)=0.d0; vx(nn-nx+1:nn)=0.d0
    vx(1:nn:nx)=0.d0; vx(nx:nn:nx)=0.d0
    vy(1:nx)=0.d0; vy(nn-nx+1:nn)=0.d0
    vy(1:nn:nx)=0.d0; vy(nx:nn:nx)=0.d0
    uplift = 0.d0
    river = 0.d0
    hillslope = 0.d0
    call FastScape_Init_H(h)
    call FastScape_Set_U(uplift)
    call FastScape_Set_V(vx,vy)
    call FastScape_Set_Erosional_Parameters(river,-1.d0,0.4d0,1.d0, &
                                             hillslope,-1.d0,0.d0,0.d0,-2.d0)
    basement = h
    area = 0.d0
    erosion = 0.d0
    call write_grid('output/fortran_rotation_initial.csv',nx,ny,length,length, &
                    h,basement,area,erosion)
    do step=1,nstep
      call FastScape_Execute_Step()
    enddo
    call FastScape_Copy_H(h)
    call FastScape_Copy_Basement(basement)
    call FastScape_Copy_Drainage_Area(area)
    call FastScape_Copy_Erosion_Rate(erosion)
    call write_grid('output/fortran_rotation_final.csv',nx,ny,length,length, &
                    h,basement,area,erosion)
    call FastScape_Destroy()
  end subroutine rigid_rotation_case

  subroutine write_grid(filename,nx,ny,xlength,ylength,h,basement,area,erosion)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nx, ny
    double precision, intent(in) :: xlength, ylength
    double precision, intent(in) :: h(nx*ny), basement(nx*ny)
    double precision, intent(in) :: area(nx*ny), erosion(nx*ny)
    integer :: unit, i, j, node
    double precision :: x, y

    open(newunit=unit,file=filename,status='replace',action='write')
    write(unit,'(a)') 'x_m,y_m,elevation_m,basement_m,drainage_area_m2,erosion_rate_m_per_year'
    do j=1,ny
      y = dble(j-1)*ylength/dble(ny-1)
      do i=1,nx
        x = dble(i-1)*xlength/dble(nx-1)
        node = i+(j-1)*nx
        write(unit,'(*(g0,:,","))') x,y,h(node),basement(node),area(node),erosion(node)
      enddo
    enddo
    close(unit)
  end subroutine write_grid

end program fortran_reference_benchmarks
