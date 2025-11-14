program CM_v1
   use iso_fortran_env
   implicit none

   integer,  parameter :: dp    = real64
   real(dp), parameter :: pi    = 3.14159265358979323846264338327950288419_dp
   real(dp), parameter :: temp  = 0.125_dp
   real(dp), parameter :: tau   = 5.0_dp             ! Thermostat time constant
   real(dp), parameter :: dt    = 0.1_dp     
   integer,  parameter :: steps = 1200                 
   integer,  parameter :: equil_steps = 15000    
   integer,  parameter :: ntraj = 10000
   real(dp), parameter :: beta = 1.0/temp

   real(dp)    :: q, p, f
   real(dp)    :: qacf(steps), qacf_sum(steps), time_lag(steps)
   real(dp),    allocatable :: q_series(:)
   integer     :: i, j, k, traj


   call random_seed()

   qacf = 0.0_dp
   qacf_sum = 0.0_dp

   do traj = 1, ntraj
      print *, 'Trajectory', traj, '/', ntraj

      allocate(q_series(steps))

      call init_config(q, p)
      do i = 1, equil_steps
         call step_vv(q, p, f, dt)
         call thermostat(p, dt, tau, temp)
      end do

      do i = 1, steps
         call step_vv(q, p, f, dt)
         q_series(i) = q
      end do

      call compute_vacf(q_series, steps, time_lag, qacf, steps, dt)
      qacf_sum = qacf_sum + qacf

      deallocate(q_series)
   end do

   qacf = qacf_sum / real(ntraj, dp)
   call write_two_col('n1_beta8_CMD_v1.dat', time_lag, qacf, steps)

contains

   subroutine init_config(q, p)
      real(dp), intent(out) :: q, p
      integer :: i
      q = 0.0_dp
      p = sqrt(temp) * randn()
   end subroutine init_config

   pure subroutine force(q, f)
      real(dp), intent(in)  :: q
      real(dp), intent(out) :: f
      real(dp), parameter :: c2 = 0.52291785536284252_dp
      real(dp), parameter :: c4 = 0.12921633977845648_dp
      real(dp), parameter :: c6 = 0.019357304963654978_dp

      !f = -q**3 !Classical MD
      f = -(2.0_dp*c2*q + 4.0_dp*c4*q**3 + 6.0_dp*c6*q**5) !Centroid MD

   end subroutine force

   pure subroutine step_vv(q, p, f, dt)
      real(dp), intent(inout) :: q, p, f
      real(dp), intent(in)    :: dt
      real(dp) :: halfdt
      integer  :: i
      halfdt = 0.5_dp*dt
      call force(q, f)
      p = p + halfdt * f
      q = q + dt * p
      call force(q, f)
      p = p + halfdt * f
   end subroutine step_vv

   subroutine thermostat(p, dt, tau, T)
      real(dp), intent(inout) :: p
      real(dp), intent(in)    :: dt, tau, T

      integer  :: dof
      real(dp) :: K, Kbar, c, s, r1, sum_r2, alpha2, alpha, factor

      dof  = 1
      K    = 0.5_dp * p*p
      Kbar = 0.5_dp * real(dof,dp) * T

      c = exp(-dt/tau)
      s = 1.0_dp - c

      r1     = randn()

      factor = (Kbar / (real(dof,dp)*max(K, tiny(1.0_dp))))
      alpha2 = c                                            &
            + factor * s * (r1*r1)               &
            + 2.0_dp * exp(-0.5_dp*dt/tau)                &
               * sqrt( factor * s ) * r1

      alpha2 = max(alpha2, 0.0_dp)
      alpha  = sqrt(alpha2)
      p      = alpha * p
   end subroutine thermostat

   real(dp) function randn()
      real(dp) :: u1, u2
      call random_number(u1)
      call random_number(u2)
      u1    = max(u1, 1.0e-12_dp)
      randn = sqrt(-2.0_dp*log(u1)) * cos(2.0_dp*pi*u2)
   end function randn

   subroutine compute_vacf(q_series, nsteps, time_lag, vacf_out, ntau_out, dt)
      real(dp), intent(in)  :: q_series(:)
      integer, intent(in)   :: nsteps, ntau_out
      real(dp), intent(in)  :: dt
      real(dp), intent(out) :: vacf_out(:), time_lag(:)
      integer :: i, t0, norig


      do i = 1, nsteps
         time_lag(i) = (i-1)*dt
         norig = nsteps - (i-1)
         vacf_out(i) = 0.0_dp
         do t0 = 1, norig
            vacf_out(i) = vacf_out(i) + q_series(t0)*q_series(t0 + i - 1)
         end do
         vacf_out(i) = vacf_out(i) / real(norig,dp)
      end do

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

end program CM_v1
