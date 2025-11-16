!===============================================================
! 1D Channel Flow (Re_H = 20)
! 
! Steady solution of a 1D diffusion-like equation using
! an implicit theta-method and tridiagonal solver (TRIDAG).
!
! Outputs:
!   - iteration_error.plt : iteration vs. error
!   - u_y.plt             : numerical and analytical u(y)
!
! Author : Mohammad E. Heravifard
!===============================================================
program Project1
    implicit none

    integer            :: iter, i, k, N, it
    real(8)            :: rho_star, H_star, U_star, mu_star
    real(8)            :: y_star, Re_H, y, dy, dt, r, error
    real(8)            :: theta, sigma, temp, pi, t_star, time
    real(8), allocatable :: u_old(:), u_new(:), y_space(:)
    real(8), allocatable :: du(:), Ondiag(:), up(:), low(:)
    real(8), allocatable :: rhs(:), u_exact(:)

    !-------------------------------
    ! Physical and numerical params
    !-------------------------------
    rho_star = 1.0d0
    H_star   = 1.0d0
    U_star   = 1.0d0

    Re_H = 20.0d0
    mu_star = rho_star * U_star * H_star / Re_H

    N = 41
    allocate(u_old(N), u_new(N), y_space(N), du(N), Ondiag(N), up(N), low(N), rhs(N), u_exact(N))

    y_star = 2.0d0               ! physical width (nondimensional)
    y      = y_star / H_star

    dy = y / real(N-1, kind=8)

    dt = 0.05d0                  ! time step
    r  = dt / (Re_H * dy * dy)   ! diffusion-like parameter

    ! theta-method (0: explicit, 0.5: Crank-Nicolson, 1: implicit)
    theta = 1.0d0                ! fully implicit

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
    ! Iteration log file
    !-------------------------------
    open(5050, file="iteration_error.plt")
    write(5050,*) 'variables="Iteration", "Error"'

    !-------------------------------
    ! Main iteration loop
    !-------------------------------
    error = 1.0d10
    it    = 300

    do iter = 1, it

        ! Discretization (combined theta-method)
        do i = 2, N-1
            Ondiag(i) = 1.0d0 + 2.0d0 * r * theta
            up(i)     = -r * theta
            low(i)    = -r * theta

            rhs(i) = r * dy * dy + r * (u_old(i-1) - 2.0d0*u_old(i) + u_old(i+1))
        end do

        ! Solve tridiagonal system for du: A * du = rhs
        call TRIDAG(2, N-1, low, Ondiag, up, rhs)
        du = rhs

        ! Boundary conditions
        du(1) = 0.0d0
        du(N) = 0.0d0

        ! Update solution
        u_new = u_old + du
        u_old = u_new

        ! Compute error (L2 norm of du)
        error = 0.0d0
        do i = 2, N-1
            error = error + du(i)*du(i)
        end do
        error = sqrt(error / real(N-2,kind=8))

        write(5050,*) iter, error

        ! Optional steady-state criterion
        if (error < 1.0d-8) then
            exit
        end if

    end do

    close(5050)

    !-------------------------------
    ! Analytical solution
    !-------------------------------
    pi   = acos(-1.0d0)
    time = 0.0d0   ! set physical time here if needed

    do i = 1, N
        sigma = 0.0d0
        do k = 1, 10000
            temp = ((-1.0d0)**(k+1) / ( (2*k-1)**3 )) * &
                   cos( (2*k-1)*pi*y_space(i)/2.0d0 ) * &
                   exp( - ( (2*k-1)**2 * pi*pi * time / (4.0d0*Re_H) ) )
            sigma = sigma + temp
        end do
        u_exact(i) = 0.5d0 * ( 1.0d0 - y_space(i)*y_space(i) - (32.0d0/(pi*pi*pi))*sigma )
    end do

    !-------------------------------
    ! Output numerical + analytical
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

!===============================================================
! TRIDAG : simple Thomas algorithm for tridiagonal systems
!
! Solves for d in-place: A * d = d,
! where A has:
!   - low(i)   on subdiagonal (i>ilo)
!   - diag(i)  on main diagonal
!   - up(i)    on superdiagonal (i<ihi)
!
! System is assumed defined on indices [ilo:ihi].
!===============================================================
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
