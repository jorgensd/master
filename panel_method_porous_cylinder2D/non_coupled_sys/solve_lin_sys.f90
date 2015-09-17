subroutine solve_lin_sys(N, P, Q)

  ! compile with:
  !    f2py -c -m solve_lin_sys solve_lin_sys.f90

  integer                          :: N
  real(kind=8), dimension(N, N)    :: P
  real(kind=8), dimension(N, 3)    :: Q
  integer, dimension(N)            :: pivot
  integer                          :: ok

  ! The special !f2py comment line, specifies input arguments and
  ! objects to be returned from the loop routine.
  !
  !f2py intent(in) P, Q, N
  !f2py intent(out) Q
  !f2py depend(N) P, Q

  ! Calculating the velocity potential for each panel
  ! using the method: DGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO)
  ! from the external library LAPACK. B contains the new solution
  call dgesv(N-1, 3, P, N-1, pivot, Q, N-1, ok)
  
end subroutine solve_lin_sys
