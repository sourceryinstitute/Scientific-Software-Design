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
  use lorenz_module ,only : lorenz,lorenz_
  use timed_lorenz_module ,only : timed_lorenz,timed_lorenz_
  use explicit_euler_module ,only : explicit_euler
  use runge_kutta_2nd_module ,only : runge_kutta_2nd
  implicit none ! Prevent implicit typing

  type(explicit_euler)  :: lorenz_integrator       !Integration strategy
  type(runge_kutta_2nd) :: timed_lorenz_integrator !Integration strategy
  type(lorenz)          :: attractor ! Lorenz equation/state abstraction
  type(timed_lorenz)    :: timed_attractor ! Time-stamped abstraction
  integer               :: step            ! Time step counter
  integer ,parameter    :: num_steps=2000,space_dimension=3

  ! Lorenz parameters and step size
  real    ,parameter    :: sigma=10.,rho=28.,beta=8./3.,dt=0.01 
  real    ,parameter ,dimension(space_dimension) &
                        :: initial=(/1.,1.,1./)

  ! variables for testing
  real ,dimension(:), allocatable  :: results_check
  real ,parameter ,dimension(space_dimension+1)  &
            :: results= (/20.0003624,-5.19188404,1.74187779,31.4894981/)
  real ,parameter                  :: epsilon=1.0e-6                            

  ! Initialize and choose strategy
  attractor = &
    lorenz_(initial,sigma,rho,beta,lorenz_integrator) 
  print *,'lorenz attractor:' 
  print *,attractor%output()

  ! Run explicit Euler at increased resolution for comparison to RK2
  do step=1,4*num_steps  
    call attractor%integrate(dt/4.)
    print *,attractor%output()
  end do
  ! Re-initialize, choose new strategy
  timed_attractor = &
    timed_lorenz_(initial,sigma,rho,beta,timed_lorenz_integrator)
  print *,''
  print *,'timed_lorenz attractor:'
  print *,timed_attractor%output()
  do step=1,num_steps
    call timed_attractor%integrate(dt)
    print *,timed_attractor%output()
    if (step==num_steps) then
      results_check=timed_attractor%output()
      if (abs(sum(results_check-results))<=1.0e-6) print *, 'Test passed'
    end if 
  end do
end program main
