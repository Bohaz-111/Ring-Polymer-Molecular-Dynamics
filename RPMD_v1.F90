program RPMD_v1
   use iso_fortran_env
   implicit none

   integer,  parameter :: dp    = real64
   integer,  parameter :: n     = 32              ! Number of beads
   real(dp), parameter :: pi    = 3.14159265358979323846264338327950288419_dp
   real(dp), parameter :: temp  = 0.125_dp
   real(dp), parameter :: tau   = 5.0_dp             ! Thermostat time constant
   real(dp), parameter :: dt    = 0.1_dp     
   integer,  parameter :: steps = 1200                 
   integer,  parameter :: equil_steps = 10000    
   integer,  parameter :: ntraj = 20000
   real(dp), parameter :: beta = 1.0/temp
   real(dp), parameter :: beta_n = beta/real(n, dp)

   real(dp)    :: q(n), p(n), f(n)
   real(dp)    :: qacf(steps), qacf_sum(steps), time_lag(steps)
   real(dp),    allocatable :: q_series(:, :)
   integer     :: i, j, k, traj


   call random_seed()
   
   qacf = 0.0_dp
   qacf_sum = 0.0_dp

   do traj = 1, ntraj
      print *, 'Trajectory', traj, '/', ntraj

      allocate(q_series(steps, n))

      call init_config(q, p)
      do i = 1, equil_steps
         call step_vv(q, p, f, dt)
         call thermostat(p, dt, tau, 1.0_dp/beta_n)
      end do

      do i = 1, steps
         call step_vv(q, p, f, dt)
         q_series(i, :) = q(:)
      end do

      call compute_vacf(q_series, steps, n, time_lag, qacf, steps, dt)
      qacf_sum = qacf_sum + qacf

      deallocate(q_series)
   end do

   qacf = qacf_sum / real(ntraj, dp)
   call write_two_col('n32_beta8_RPMD_v1.dat', time_lag, qacf, steps)

contains

   subroutine init_config(q, p)
      real(dp), intent(out) :: q(:), p(:)
      integer :: i, nn
      nn = size(q)
      q = 0.0_dp
      do i = 1, nn
         p(i) = sqrt(1.0_dp/beta_n) * randn()
      end do
   end subroutine init_config

   pure subroutine force(q, f)
      real(dp), intent(in)  :: q(:)
      real(dp), intent(out) :: f(:)
      integer :: jp, jm, j, nn
      real(dp) :: delta_q(size(q)), delta_qp, delta_qm
      nn = size(q)
    
      do j = 1, nn
         jm = j-1; if (jm < 1)  jm = nn
         jp = j+1; if (jp > nn) jp = 1
       
         f(j) = -1.0_dp/beta_n**2 * (2.0_dp*q(j) - q(jp) - q(jm)) - q(j)**3
      end do
   end subroutine force

   pure subroutine step_vv(q, p, f, dt)
      real(dp), intent(inout) :: q(:), p(:), f(:)
      real(dp), intent(in)    :: dt
      real(dp) :: halfdt
      integer  :: i, nn
      nn = size(q)
      halfdt = 0.5_dp*dt
      call force(q, f)
      p = p + halfdt * f
      q = q + dt * p
      call force(q, f)
      p = p + halfdt * f
   end subroutine step_vv

   subroutine thermostat(p, dt, tau, T)
      real(dp), intent(inout) :: p(:)
      real(dp), intent(in)    :: dt, tau, T

      integer  :: np, dof, kshape
      real(dp) :: K, Kbar, c, s, r1, sum_r2, alpha2, alpha, factor

      np   = size(p)
      dof  = np
      K    = 0.5_dp * sum(p*p)
      Kbar = 0.5_dp * real(dof,dp) * T

      c = exp(-dt/tau)
      s = 1.0_dp - c

      r1     = randn()
      kshape = (dof-2)/2
      sum_r2 = rand_gamma(kshape) + randn()**2

      factor = (Kbar / (real(dof,dp)*max(K, tiny(1.0_dp))))
      alpha2 = c                                            &
            + factor * s * (r1*r1 + sum_r2)               &
            + 2.0_dp * exp(-0.5_dp*dt/tau)                &
               * sqrt( factor * s ) * r1

      alpha2 = max(alpha2, 0.0_dp)
      alpha  = sqrt(alpha2)
      p      = alpha * p
   end subroutine thermostat

   real(dp) function rand_gamma(k)
      integer, intent(in) :: k
      real(dp) :: v(k)
      integer  :: i
      call random_number(v)
      do i = 1, k
         v(i) = max(v(i), 1.0e-12_dp)
      end do
      v = -log(v)
      rand_gamma = 2.0_dp * sum(v)
   end function rand_gamma

   real(dp) function randn()
      real(dp) :: u1, u2
      call random_number(u1)
      call random_number(u2)
      u1    = max(u1, 1.0e-12_dp)
      randn = sqrt(-2.0_dp*log(u1)) * cos(2.0_dp*pi*u2)
   end function randn

   subroutine compute_vacf(q_series, nsteps, nk, time_lag, vacf_out, ntau_out, dt)
      real(dp), intent(in)  :: q_series(:, :)
      integer, intent(in)   :: nsteps, nk, ntau_out
      real(dp), intent(in)  :: dt
      real(dp), intent(out) :: vacf_out(:), time_lag(:)
      integer :: i, t0, norig
      real(dp), allocatable :: xcent(:)

      allocate(xcent(nsteps))
      do i = 1, nsteps
         xcent(i) = sum(q_series(i,1:nk)) / real(nk,dp)
      end do

      do i = 1, nsteps
         time_lag(i) = (i-1)*dt
         norig = nsteps - (i-1)
         vacf_out(i) = 0.0_dp
         do t0 = 1, norig
            vacf_out(i) = vacf_out(i) + xcent(t0)*xcent(t0 + i - 1)
         end do
         vacf_out(i) = vacf_out(i) / real(norig,dp)
      end do

      deallocate(xcent)
   end subroutine compute_vacf

   subroutine write_two_col(fname, x, y, n)
      character(*), intent(in) :: fname
      real(dp),     intent(in) :: x(:), y(:)
      integer,      intent(in) :: n
      integer :: u,i
      open(newunit=u, file=fname, status='replace', action='write')
         write(u,'(A)') '# Time    Value'
         do i=1,n
               write(u,*) x(i), y(i)
         end do
      close(u)
      print *, 'saved to:  ', fname
   end subroutine write_two_col

end program RPMD_v1
