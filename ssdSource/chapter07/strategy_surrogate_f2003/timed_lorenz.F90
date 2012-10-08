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

module timed_lorenz_module
  use lorenz_module  ,only :lorenz,lorenz_,integrand ! Parent and grandparent
  use strategy_module,only :strategy         ! Time integration strategy
  implicit none              ! Prevent implicit typing
  private                    ! Hide everything by default
  public :: integrand        ! Expose abstract integrand and integrator
  public :: timed_lorenz,timed_lorenz_

  type ,extends(lorenz) :: timed_lorenz
    private
    real :: time                               ! time stamp
  contains
    procedure ,public :: t => dTimed_lorenz_dt ! time deriv. (eval. RHS)
    procedure ,public :: add       => add_timed_lorenz
    procedure ,public :: multiply  => multiply_timed_lorenz
    procedure ,public :: assign    => assign_timed_lorenz
    procedure ,public :: output                ! accessor: return state
  end type timed_lorenz

  interface timed_lorenz_
    procedure constructor
  end interface

contains

  ! Constructor: allocate and initialize state
  type(timed_lorenz) function constructor(initial,s,r,b,this_strategy)
    real ,dimension(:)  ,intent(in)  :: initial
    real                ,intent(in)  :: s ,r ,b ! Lorenz parameters:
                                                ! sigma, rho and beta
    class(strategy),intent(in) :: this_strategy ! marching algorithm
    constructor%lorenz = lorenz_ &
      (initial_state=initial,s=s,r=r,b=b,this_strategy=this_strategy)
    constructor%time   = 0.
  end function

  ! time derivative (expresses evolution equations)
  function dTimed_lorenz_dt(this) result(dState_dt)
    class(timed_lorenz),intent(in)  :: this
    class(integrand)   ,allocatable :: dState_dt 
    type(timed_lorenz) ,allocatable :: local_dState_dt

    allocate(local_dState_dt)
    local_dState_dt%time   = 1.                 ! dt/dt = 1.
    local_dState_dt%lorenz = this%lorenz%t()    ! delegate to parent
    ! avoid unnecessary memory allocation
    call move_alloc(local_dState_dt,dState_dt)
  end function

  ! add two instances of timed_lorenz
  function add_timed_lorenz(lhs,rhs) result(sum)
    class(timed_lorenz)     ,intent(in)  :: lhs
    class(integrand)        ,intent(in)  :: rhs
    class(integrand)        ,allocatable :: sum
    type(timed_lorenz)      ,allocatable :: local_sum

    select type(rhs)
      class is (timed_lorenz)
        allocate(local_sum)
        local_sum%time   = lhs%time   + rhs%time
        local_sum%lorenz = lhs%lorenz + rhs%lorenz
      class default
        stop 'add_timed_lorenz: type not supported'
    end select

    ! avoid unnecessary memory allocation
    call move_alloc(local_sum,sum)
  end function

  ! multiply one instance of timed_lorenz by a real scalar
  function multiply_timed_lorenz(lhs,rhs) result(product)
    class(timed_lorenz)     ,intent(in)  :: lhs
    real                    ,intent(in)  :: rhs
    class(integrand)        ,allocatable :: product
    type(timed_lorenz)      ,allocatable :: local_product

    allocate(local_product) 
    local_product%time   = lhs%time  * rhs
    local_product%lorenz = lhs%lorenz* rhs

    ! transfer allocation from local_product to result product
    call move_alloc(local_product,product)
  end function

  ! assign one instance of timed_lorenz to another
  subroutine assign_timed_lorenz(lhs,rhs)
    class(timed_lorenz)     ,intent(inout) :: lhs
    class(integrand)        ,intent(in)    :: rhs
    select type(rhs)
      class is (timed_lorenz)
        lhs%time   = rhs%time
        lhs%lorenz = rhs%lorenz
      class default
        stop 'assign_timed_lorenz: type not supported'
    end select
  end subroutine

  ! return state
  function output(this) result(coordinates)
    class(timed_lorenz) ,intent(in)  :: this
    real ,dimension(:)  ,allocatable :: coordinates

    ! assignment automatically allocates coordinates
    coordinates = [ this%time, this%lorenz%output() ]
  end function
end module timed_lorenz_module
