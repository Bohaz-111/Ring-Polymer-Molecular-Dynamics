program qm_kubo_v2
  use iso_fortran_env
  integer, parameter :: dp = real64

  real(dp), parameter :: pi = 3.141592653589793238_dp
  real(dp), parameter :: on2m = 0.5_dp
  real(dp), parameter :: beta = 8.0_dp
  integer,  parameter :: Ngrid = 401
  real(dp), parameter :: Lbox  = 3.0_dp      ! domain [-Lbox, Lbox]
  real(dp), parameter :: mLbox  = -Lbox
  integer,  parameter :: Nkeep = 5         ! eigenstates used in thermal sums

  real(dp), parameter :: tmax = 20.0_dp
  real(dp), parameter :: dt   = 0.05_dp
  integer,  parameter :: Nt   = int(tmax/dt) + 1

  real(dp) :: x(Ngrid), w(Ngrid), Hmat(Ngrid,Ngrid), Tmat(Ngrid,Ngrid)
  real(dp) :: eval(Ngrid), evec(Ngrid,Ngrid)
  real(dp) :: Ccorr(Nt), time_grid(Nt), Zpart, dE, xmn, wmn, hstep
  integer  :: igrid, jgrid, istate, jstate, it

  call sindvr(on2m,mLbox,Lbox,Ngrid,Hmat,w,x)

  do igrid = 1, Ngrid
     Hmat(igrid,igrid) = Hmat(igrid,igrid) + 0.25_dp * x(igrid)**4
  end do

  call jacobi_symeig(Ngrid, Hmat, eval, evec)

  Zpart = 0.0_dp
  do jstate = 1, Nkeep
     Zpart = Zpart + exp(-beta*eval(jstate))
  end do

  ! Time grid
  do it = 1, Nt
     time_grid(it) = dt * real(it-1, dp)
  end do
  Ccorr = 0.0_dp

  ! Accumulate Kubo-transformed C_xx(t)
  ! C~(t) = (1/Z) sum_{m != n} [(e^{-βE_m} - e^{-βE_n})/(β(E_n-E_m))] |x_mn|^2 cos[(E_n-E_m)t]
  do istate = 1, Nkeep
     do jstate = 1, Nkeep
        if (istate == jstate) cycle
        dE = eval(jstate) - eval(istate)
        xmn = x_matrix_element(Ngrid, evec(:,istate), evec(:,jstate), x)
        wmn = (exp(-beta*eval(istate)) - exp(-beta*eval(jstate))) / (beta * dE)
        do it = 1, Nt
           Ccorr(it) = Ccorr(it) + wmn * (xmn*xmn) * cos(dE * time_grid(it))
        end do
     end do
  end do
  Ccorr = Ccorr / Zpart

  call write_two_col('QM_beta8_quartic_Kubo_v2.dat', time_grid, Ccorr, Nt)
  print *, 'Wrote: QM_beta8_quartic_Kubo_v2.dat'

contains

  subroutine sindvr(on2m, a, b, n, t, w, x)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: on2m, a, b
    real(dp), intent(out) :: t(n,n), w(n), x(n)

    integer  :: i, j, k, m
    real(dp) :: alfa, beta, cosa, sina, cosb, sinb
    real(dp) :: dx, wt, t0, temp
    real(dp), allocatable :: s(:)

    m    = n + 1
    alfa = 0.5_dp*on2m*(pi/(b-a))**2 
    beta = pi/(2.0_dp*m)

    cosa = 1.0_dp
    sina = 0.0_dp
    cosb = cos(beta)
    sinb = sin(beta)

    dx = (b-a)/real(m,dp)
    wt = sqrt(dx)
    t0 = alfa*(2.0_dp*m*m + 1.0_dp)/3.0_dp

    allocate(s(2*n))
    do k = 1, 2*n
      alfa = -alfa
      temp = cosa*sinb + sina*cosb
      cosa = cosa*cosb - sina*sinb
      sina = temp
      s(k) = alfa/(sina*sina)
    end do

    t = 0.0_dp
    do j = 2, n
      do i = 1, j-1
        t(i,j) = s(j-i) - s(j+i)
        t(j,i) = t(i,j)
      end do
    end do

    do j = 1, n
      t(j,j) = t0 - s(2*j)
      w(j)   = wt
      x(j)   = a + real(j,dp)*dx
    end do

    deallocate(s)
  end subroutine sindvr

!===================== Symmetric Jacobi eigensolver =====================
  ! Diagonalizes real symmetric A (in/out). Returns eigenvalues (ascending)
  ! in w and eigenvectors in columns of V.
  subroutine jacobi_symeig(n, A, w, V)
    integer,  intent(in)    :: n
    real(dp), intent(inout) :: A(n,n)
    real(dp), intent(out)   :: w(n)
    real(dp), intent(out)   :: V(n,n)

    integer :: p, q, i, j, sweep, max_sweeps
    real(dp) :: app, aqq, apq, tau, t, c, s, g, thresh, offmax
    real(dp) :: Akkp, Akkq

    V = 0.0_dp
    do i = 1, n
      V(i,i) = 1.0_dp
    end do

    max_sweeps = 100

    do sweep = 1, max_sweeps
      offmax = max_offdiag(n, A)
      if (offmax < 1.0e-12_dp) exit

      do p = 1, n-1
        do q = p+1, n
          apq = A(p,q)
          if (abs(apq) <= 1.0e-20_dp) cycle
          app = A(p,p)
          aqq = A(q,q)

          tau = (aqq - app) / (2.0_dp * apq)
          if (tau >= 0.0_dp) then
            t = 1.0_dp / (tau + sqrt(1.0_dp + tau*tau))
          else
            t = 1.0_dp / (tau - sqrt(1.0_dp + tau*tau))
          end if
          c = 1.0_dp / sqrt(1.0_dp + t*t)
          s = t * c

          do i = 1, n
            if (i /= p .and. i /= q) then
              Akkp   = A(i,p)
              Akkq   = A(i,q)
              A(i,p) = c*Akkp - s*Akkq
              A(p,i) = A(i,p)
              A(i,q) = s*Akkp + c*Akkq
              A(q,i) = A(i,q)
            end if
          end do

          A(p,p) = c*c*app - 2.0_dp*s*c*apq + s*s*aqq
          A(q,q) = s*s*app + 2.0_dp*s*c*apq + c*c*aqq
          A(p,q) = 0.0_dp
          A(q,p) = 0.0_dp

          do i = 1, n
            Akkp = V(i,p)
            Akkq = V(i,q)
            V(i,p) = c*Akkp - s*Akkq
            V(i,q) = s*Akkp + c*Akkq
          end do
        end do
      end do
    end do

    do i = 1, n
      w(i) = A(i,i)
    end do

    call sort_eigs(n, w, V)
  end subroutine jacobi_symeig

  pure real(dp) function max_offdiag(n, A) result(val)
    integer,  intent(in) :: n
    real(dp), intent(in) :: A(n,n)
    integer :: i, j
    val = 0.0_dp
    do i = 1, n-1
      do j = i+1, n
        val = max(val, abs(A(i,j)))
      end do
    end do
  end function max_offdiag

  subroutine sort_eigs(n, w, V)
    integer,  intent(in)    :: n
    real(dp), intent(inout) :: w(n)
    real(dp), intent(inout) :: V(n,n)
    integer :: i, j, k
    real(dp) :: wk
    real(dp), allocatable :: tmpcol(:)

    allocate(tmpcol(n))

    do i = 1, n-1
      k = i
      wk = w(i)
      do j = i+1, n
        if (w(j) < wk) then
          k = j
          wk = w(j)
        end if
      end do
      if (k /= i) then
        w(k) = w(i)
        w(i) = wk
        tmpcol = V(:,i)
        V(:,i) = V(:,k)
        V(:,k) = tmpcol
      end if
    end do

    deallocate(tmpcol)
  end subroutine sort_eigs

  pure function x_matrix_element(N, v, w, x) result(val)
    integer,  intent(in) :: N
    real(dp), intent(in) :: v(N), w(N), x(N)
    real(dp) :: val
    integer :: i
    val = 0.0_dp
    do i = 1, N
       val = val + v(i) * x(i) * w(i)
    end do
  end function x_matrix_element

  subroutine write_two_col(fname, xv, yv, n)
    character(*), intent(in) :: fname
    real(dp),     intent(in) :: xv(:), yv(:)
    integer,      intent(in) :: n
    integer :: u,i
    open(newunit=u, file=fname, status='replace', action='write')
      write(u,'(A)') '# t    Ctilde_xx(t)'
      do i = 1, n
        write(u,'(2ES24.15)') xv(i), yv(i)
      end do
    close(u)
  end subroutine write_two_col

end program qm_kubo_v2
