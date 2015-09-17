subroutine loop(n, z, cz, l, n1, n2, p, q)

  ! compile with:
  !    f2py -c -m loop loop.f90 (--fcompiler=gfortran --opt='-O3')

  integer                       :: n
  integer, parameter            :: dp = kind(0.d0)
  real(8)                       :: pi = 4 * atan(1.0_dp) ! or 3.1415926535897932
  real(8), dimension(n)         :: n1, n2, l
  real(8), dimension(n, n)      :: p
  real(8), dimension(n, 3)      :: q
  real(8)                       :: log_r, b, c
  complex(8), dimension(n+1)    :: z
  complex(8), dimension(n)      :: cz

  ! The special !f2py comment line, specifies input arguments and
  ! objects to be returned from the loop routine.
  !
  !f2py intent(in) n, z, cz, l, n1, n2
  !f2py intent(in, out) p, q
  !f2py depend(n) p, q, z, cz, l, n1, n2
  
  ! Calculating the opening angle using the law of cosines.
  ! Calculating the rhs integral using the trapezoidal rule
  do i=1, n
     do j=1, n
        b = abs( cz(i) - z(j+1) )
        c = abs( cz(i) - z(j) )
        p(i, j) = -acos( (b**2 + c**2 - l(j)**2) / (2*b*c) )
        log_r =  0.5_dp * ( log(b) + log(c) )
        q(i, 1) = q(i, 1) + n1(j) * log_r * l(j)
        q(i, 2) = q(i, 2) + n2(j) * log_r * l(j)
        q(i, 3) = q(i, 3) + (real(cz(j))*n2(j) - aimag(cz(j))*n1(j)) * l(j) * log_r
     end do
  end do

  ! Replace NaN with 0
  do i = 1, n
     do j = 1, n
        if ( isnan(p(i, j)) ) then
           p(i, j) = 0
        end if
     end do
  end do
  
  ! Put -pi on the diagonal
  do i = 1, n
     p(i, i) = -pi
  end do
  
end subroutine loop
