subroutine loop(n, z, cz, l, n1, theta, intLog, p, q, mu, tau, a, d, w, rho)

  ! compile with:
  !    f2py -c -m loop loop.f90 -llapack -lblas --opt='-O3' (--fcompiler=gfortran)

  ! mu - Discharge coefficient (Molin p. 3)
  ! tau - porosity parameter, 0 is solid, 1 is open
  ! w - frequency omega of oscillations
  ! a - cylinder radius
  ! d - amplitude of motion

  integer                        :: n, ok
  integer, dimension(2*n)        :: pivot
  integer, parameter             :: dp = kind(0.d0)
  real(8)                        :: pi = 4 * atan(1.0_dp) ! or 3.1415926535897932
  real(8), dimension(n)          :: n1, l, tmp
  real(8), dimension(n, n)       :: theta, intLog
  real(8)                        :: b, c, mu, tau, a, d, w, rho, K, KC
  complex(8), dimension(2*n,2*n) :: p
  complex(8), dimension(2*n)     :: q
  complex(8), dimension(n+1)     :: z
  complex(8), dimension(n)       :: cz

  ! The special !f2py comment line, specifies input arguments and
  ! objects to be returned from this routine.
  !
  !f2py intent(in) n, z, cz, l, n1, n2, theta, intLog, p, mu, tau, a, d, w, rho
  !f2py intent(in, out) q
  !f2py depend(n) z, cz, l, n1, n2, theta, intLog, p, q
  
  ! Calculating the opening angle using the law of cosines.
  ! Calculating the integral of log using the trapezoidal rule
  do i=1, n
     do j=1, n
        b = abs( cz(i) - z(j+1) )
        c = abs( cz(i) - z(j) )
        theta(i, j) = -acos( (b**2 + c**2 - l(j)**2) / (2.0_dp*b*c) )
        intLog(i, j) = 0.5_dp * ( log(b) + log(c) ) * l(j)
     end do
  end do
  
  ! Put pi on the diagonal
  do i = 1, n
     theta(i, i) = pi
  end do

  ! Parameters for fluid and porous media
  K = (1.0_dp - tau) / (mu*tau**2)  ! Resistance coefficient (Molin p. 3)
  KC = d / (2.0_dp*a)               ! Keulegan-Carpenter number
  b = 1.0_dp / (K * KC)             ! Porosity coefficient from Chwang

  !
  p(n+1:, :n) = -theta
  p(:n, n+1:) = theta
  p(n+1:, n+1:) = p(n+1:, n+1:) - cmplx(0.0_dp, 2.0_dp*rho*w*b*intLog, dp)

  !
  do i = 1, n
     tmp(i) = dot_product( intLog(i, :), n1 )
  end do
  q(n+1:) = cmplx(0.0, 2.0_dp * w * d * tmp, dp)
  
  ! Calculating the velocity potential for each panel
  ! using the method: ZGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO)
  ! from the external library LAPACK. B contains the new solution
  call zgesv(2*n, 1, p, 2*n, pivot, q, 2*n, ok)
  
end subroutine loop
