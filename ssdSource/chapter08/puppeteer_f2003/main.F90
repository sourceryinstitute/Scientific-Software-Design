! Copyright (c) 2011, Damian Rouson, Jim Xia, and Xiaofeng Xu.
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!     * Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright
!       notice, this list of conditions and the following disclaimer in the
!       documentation and/or other materials provided with the distribution.
!     * Neither the names of Damian Rouson, Jim Xia, and Xiaofeng Xu nor the
!       names of any other contributors may be used to endorse or promote products
!       derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL DAMIAN ROUSON, JIM XIA, and XIAOFENG XU BE LIABLE 
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
! ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

program main
  use air_module  ,only : air,air_
  use cloud_module  ,only :cloud,cloud_
  use ground_module  ,only : ground,ground_
  use atmosphere_module ,only : atmosphere,atmosphere_
  use global_parameters_module ,only :debugging !print call tree if true

  implicit none ! Prevent implicit typing

  ! This code integrates the Lorenz equations over time using separate
  ! abstractions for equation and hiding the coupling of those 
  ! abstractions inside an abstraction that follows the Puppeteer design
  ! pattern of Rouson, Adalsteinsson and Xia (ACM TOMS 37:1, 2010).

  type(air)     ,allocatable :: sky   !puppet for 1st Lorenz equation
  type(cloud)   ,allocatable :: puff  !puppet for 2nd Lorenz equation
  type(ground)  ,allocatable :: earth !puppet for 3rd Lorenz equation
  type(atmosphere)           :: boundary_layer  ! Puppeteer
  integer                    :: step            ! time step
  integer ,parameter         :: num_steps=1000  ! total time steps 
  real    ,parameter         :: x=1.,y=1.,z=1.  ! initial conditions
  real                       :: t               ! time coordinate
                                ! Lorenz parameters
  real    ,parameter         :: sigma=10.,rho=28.,beta=8./3.,dt=.02 

! variables for testing
  real ,dimension(:), allocatable  :: results_check
  real ,parameter ,dimension(4)  &
            :: results= (/-10.1632519,10.0000000,-15.8904438,20.3053379/)
  real ,parameter                  :: epsilon=1.0e-6

  if (debugging) print *,'main: start'
  allocate (sky, puff, earth)
  sky = air_(x,sigma)
  puff = cloud_(y,rho)
  earth = ground_(z,beta)

  ! transfer allocations into puppeteer
  boundary_layer = atmosphere_(sky,puff,earth) 
  ! all puppets have now been deallocated

  t=0.
  write(*,'(f10.4)',advance='no') t
  print *,boundary_layer%state_vector()
  do step=1,num_steps
    call integrate(boundary_layer,dt)
    t = t + dt
    write(*,'(f10.4)',advance='no') t
    print *,boundary_layer%state_vector()
    if (step==num_steps) then
      results_check=boundary_layer%state_vector()
      if (abs(sum(results_check-results))<=1.0e-6) print *, 'Test passed'
    end if      
  end do
  if (debugging) print *,'main: end'

contains
  ! abstract Trapezoidal rule integration
  subroutine integrate(integrand,dt)
    type(atmosphere) ,intent(inout) :: integrand
    real             ,intent(in)    :: dt
    type(atmosphere) ,allocatable   :: integrand_estimate,residual
    integer ,parameter :: num_iterations=5
    integer :: newton_iteration,num_equations,num_statevars,i,j
    real ,dimension(:,:) ,allocatable  :: dRHS_dState,jacobian,identity

    if (debugging) print *,'  integrate: start'
    allocate(integrand_estimate, source=integrand)
    allocate(residual)
    do newton_iteration=1,num_iterations
      dRHS_dState   = integrand_estimate%dRHS_dV()
      num_equations = size(dRHS_dState,1)
      num_statevars = size(dRHS_dState,2)
      if (num_equations /= num_statevars) &
        stop 'integrate: ill-posed problem.'
      identity = reshape( &
        source=(/((0.,i=1,num_equations),j=1,num_statevars)/), &
        shape=(/num_equations,num_statevars/) )
      forall(i=1:num_equations) identity(i,i)=1.
      jacobian = identity - 0.5*dt*dRHS_dState
      residual = integrand_estimate - &
        ( integrand + (integrand%t() + integrand_estimate%t())*(0.5*dt))
      integrand_estimate = integrand_estimate - &
        (jacobian .inverseTimes. residual)
    end do
    integrand = integrand_estimate
    if (debugging) print *,'  integrate: end'
  end subroutine
end program main
