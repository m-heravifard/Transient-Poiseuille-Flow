!======================================================================
! 1D unsteady channel flow with periodic forcing
!
! Time-dependent solution of a 1D diffusion-type equation using
! a theta-method (Crank–Nicolson: theta = 0.5) and a tridiagonal solver.
!
! Outputs:
!   iteration_error.plt : iteration vs. error
!   u_periodic.plt      : snapshots of u(y) at selected times
!   u_center_iter.plt   : u at channel center vs. time step
!   u_y.plt             : final numerical and analytical u(y)
!
! Author  : Mohammad E. Heravifard
!======================================================================
program Project1
    implicit none

    !-------------------------------
    ! Variable declarations
    !-------------------------------
    integer  :: iter, i, k, N, max_iteration
    real(8)  :: rho_star, H_star, U_star, mu_star
    real(8)  :: y_star, Re_H, y, dy, dt, r, error, theta
    real(8)  :: sigma, temp, pi, time, a, b, end_time
    real(8)  :: A_coef, B_coef, C_coef, D_coef, F_coef
    real(8), allocatable :: u_old(:), u_new(:), y_space(:)
    real(8), allocatable :: du(:), Ondiag(:), up(:), low(:)
    real(8), allocatable :: rhs(:), u_exact(:)

    !-------------------------------
    ! Constants and parameters
    !-------------------------------
    pi = acos(-1.0d0)

    rho_star = 1.0d0
    H_star   = 1.0d0
    U_star   = 1.0d0

    Re_H = 20.0d0
    mu_star = rho_star * U_star * H_star / Re_H

    N = 41
    allocate(u_old(N), u_new(N), y_space(N), du(N), Ondiag(N), up(N), low(N), rhs(N), u_exact(N))

    y_star = 2.0d0
    y      = y_star / H_star

    dy = y / real(N-1, kind=8)

    dt = 0.01d0                      ! time step (choose as needed)
    r  = dt / (Re_H * dy * dy)
    write(*,*) "r =", r

    ! Theta-method: 0 => explicit, 0.5 => Crank–Nicolson, 1 => fully implicit
    theta = 0.5d0

    ! Periodic forcing parameters (example values)
    a = 1.0d0      ! forcing amplitude
    b = 1.0d0      ! forcing frequency (non-dimensional)

    ! Coefficients for the more detailed analytical solution
    ! TODO: set physically consistent values or formulas for:
    B_coef = 0.0d0
    C_coef = 0.0d0
    D_coef = 0.0d0
    F_coef = 0.0d0
    A_coef = U_star * mu_star * a / (H_star*H_star * rho_star * b * pi)

    !-------------------------------
    ! Grid generation
    !-------------------------------
    do i = 1, N
        y_space(i) = -1.0d0 + real(i-1,kind=8) * dy
    end do

    !-------------------------------
    ! Initial condition
    !-------------------------------
    u_old = 0.0d0
    du    = 0.0d0

    !-------------------------------
    ! Files for output
    !-------------------------------
    open(5050, file="iteration_error.plt")
    write(5050,*) 'variables="Iteration", "Error"'

    end_time      = 100.0d0
    max_iteration = nint(end_time / dt)

    open(8181, file="u_periodic.plt")
    open(8182, file="u_center_iter.plt")
    write(8182,*) 'variables="iter", "u_center"'

    !-------------------------------
    ! Time-stepping loop
    !-------------------------------
    do iter = 1, max_iteration
        time = real(iter-1,kind=8) * dt

        ! Discretization (theta-method)
        do i = 2, N-1
            Ondiag(i) = 1.0d0 + 2.0d0 * r * theta
            up(i)     = -r * theta
            low(i)    = -r * theta

            rhs(i) = r * dy*dy                                   &
     &              + r * dy*dy * a * (                           &
     &                  theta      * cos(b*pi*(time+dt))          &
     &                + (1.0d0-theta)* cos(b*pi*time)             &
     &              )                                             &
     &              + r * (u_old(i-1) - 2.0d0*u_old(i) + u_old(i+1))
        end do

        ! Solve tridiagonal system A * du = rhs
        call TRIDAG(2, N-1, low, Ondiag, up, rhs)
        du = rhs

        ! Boundary conditions (homogeneous Dirichlet)
        du(1) = 0.0d0
        du(N) = 0.0d0

        ! Update solution
        u_new = u_old + du
        u_old = u_new

        ! Store centerline value
        write(8182,*) iter, u_new(N/2+1)

        ! Error (L2 norm of increment)
        error = 0.0d0
        do i = 2, N-1
            error = error + du(i)*du(i)
        end do
        error = sqrt(error / real(N-2,kind=8))

        write(5050,*) iter, error
        write(*,*)    iter, error

        ! Write snapshots at selected iterations (example indices)
        if (iter == 92502 .or. iter == 93501 .or. iter == 94502) then
            write(8181,*) 'variables="u", "y"'
            write(8181,*) 'zone T="', iter, '"'
            do i = 1, N
                write(8181,*) u_new(i), y_space(i)
            end do
        end if

        ! Optional early stop
        if (error < 1.0d-8) then
            exit
        end if
    end do

    close(5050)
    close(8181)
    close(8182)

    !-------------------------------
    ! Analytical solution (base part)
    !-------------------------------
    sigma = 0.0d0
    do i = 1, N
        sigma = 0.0d0
        do k = 1, 10000
            temp = ((-1.0d0)**(k+1) / ((2*k-1)**3)) *               &
     &             cos( (2*k-1)*pi*y_space(i)/2.0d0 ) *             &
     &             exp( - ((2*k-1)**2) * pi*pi * time / (4.0d0*Re_H) )
            sigma = sigma + temp
        end do

        ! Base Poiseuille + transient series
        u_exact(i) = 0.5d0 * ( 1.0d0 - y_space(i)*y_space(i)        &
     &              - (32.0d0/(pi*pi*pi))*sigma )

        ! NOTE:
        ! The more detailed oscillatory analytical term that depends on
        ! A_coef, B_coef, C_coef, D_coef, F_coef is problem-specific
        ! and has been omitted here. You can re-introduce it once
        ! the correct formulas for these coefficients are defined.
        !
        ! Example placeholder (commented):
        ! u_exact(i) = u_exact(i)                                    &
        !   - sin(b*pi*time) * ( ... )                               &
        !   + cos(b*pi*time) * ( ... )

    end do

    !-------------------------------
    ! Final numerical + analytical profile
    !-------------------------------
    open(100, file="u_y.plt")
    write(100,*) 'variables="u", "y"'
    write(100,*) 'zone T="Numerical"'
    do i = 1, N
        write(100,*) u_new(i), y_space(i)
    end do

    write(100,*) 'zone T="Analytical"'
    do i = 1, N
        write(100,*) u_exact(i), y_space(i)
    end do
    close(100)

    stop
end program Project1

!======================================================================
! TRIDAG : Thomas algorithm for tridiagonal systems
!
! Solves A * d = d in-place on [ilo:ihi], where:
!   low(i)  : subdiagonal (i>ilo)
!   diag(i) : main diagonal
!   up(i)   : superdiagonal (i<ihi)
!   d(i)    : RHS on input, solution on output
!======================================================================
subroutine TRIDAG(ilo, ihi, low, diag, up, d)
    implicit none
    integer, intent(in)    :: ilo, ihi
    real(8), intent(inout) :: low(:), diag(:), up(:), d(:)

    integer :: i
    real(8) :: m

    ! Forward elimination
    do i = ilo+1, ihi
        m       = low(i) / diag(i-1)
        diag(i) = diag(i) - m * up(i-1)
        d(i)    = d(i)    - m * d(i-1)
    end do

    ! Back substitution
    d(ihi) = d(ihi) / diag(ihi)
    do i = ihi-1, ilo, -1
        d(i) = (d(i) - up(i) * d(i+1)) / diag(i)
    end do
end subroutine TRIDAG

