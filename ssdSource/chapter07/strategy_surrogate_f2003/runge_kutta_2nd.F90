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

module runge_kutta_2nd_module 
  use surrogate_module,only : surrogate
  use strategy_module ,only : strategy   ! time integration strategy
  use integrand_module,only : integrand  ! abstract integrand

  implicit none                          ! Prevent implicit typing
  private                                ! Hide everything by default

  ! 2nd-order Runge-Kutta time integration
  type, extends(strategy) ,public :: runge_kutta_2nd
  contains
    procedure, nopass :: integrate       ! integration procedure
  end type

contains

  ! Time integrator
  subroutine integrate(this,dt)
    class(surrogate) ,intent(inout) :: this      ! integrand
    real             ,intent(in)    :: dt        ! time step size
    class(integrand) ,allocatable   :: this_half ! function evaluation
                                                 ! at interval t+dt/2.

    select type (this)
      class is (integrand)
        allocate(this_half,source=this)
        this_half = this + this%t()*(0.5*dt)     ! predictor step
        this      = this + this_half%t()*dt      ! corrector step
      class default
        stop 'integrate: unsupported class'
    end select
  end subroutine
end module 
