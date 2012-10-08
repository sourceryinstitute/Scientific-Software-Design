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

module lorenz_module
  use strategy_module ,only : strategy   ! time integration strategy
  use integrand_module ,only : integrand ! Abstract integrand
  implicit none ! Prevent implicit typing
  private       ! Hide everything by default

  public :: integrand ! Expose integrand type
  public :: lorenz,lorenz_    ! Expose lorenz type and constructor

  type ,extends(integrand) :: lorenz
    private
    real ,dimension(:) ,allocatable :: state   ! solution vector
    real :: sigma ,rho ,beta                   ! Lorenz parameters
  contains
    procedure ,public :: t  => dlorenz_dt ! Time deriv. (eval. RHS)
    procedure ,public :: add      => add_lorenz
    procedure ,public :: multiply => multiply_lorenz
    procedure ,public :: assign   => assign_lorenz
    procedure ,public :: output  ! Accessor: return state
  end type

  interface lorenz_
    module procedure constructor
  end interface

contains

  ! Constructor: allocate and initialize
  type(lorenz) function constructor(initial_state,s,r,b,this_strategy)
    real ,dimension(:) ,intent(in)  :: initial_state
    real               ,intent(in)  :: s ,r ,b
    class(strategy)    ,intent(in)  :: this_strategy
    ! constructor%state is automatically allocated by the assignment
    constructor%state=initial_state
    constructor%sigma=s
    constructor%rho=r
    constructor%beta=b
    call constructor%set_quadrature(this_strategy)
  end function

  ! Time derivative (specifies evolution equations)
  function dLorenz_dt(this) result(dState_dt)
    class(lorenz)    ,intent(in)  :: this 
    class(integrand) ,allocatable :: dState_dt
    type(lorenz)     ,allocatable :: local_dState_dt
    allocate(local_dState_dt)
    call local_dState_dt%set_quadrature(this%get_quadrature())
    allocate(local_dState_dt%state(size(this%state)))
    ! 1st Lorenz equation
    local_dState_dt%state(1) = this%sigma*(this%state(2)-this%state(1))
    ! 2nd Lorenz equation
    local_dState_dt%state(2) = this%state(1)*(this%rho-this%state(3)) &
                              -this%state(2)
    ! 3rd Lorenz equation
    local_dState_dt%state(3) = this%state(1)*this%state(2)            &
                              -this%beta*this%state(3)
    local_dState_dt%sigma = 0.
    local_dState_dt%rho   = 0.
    local_dState_dt%beta  = 0.
    ! transfer the allocation from local_dState_dt to dState_dt
    call move_alloc(local_dState_dt,dState_dt)
  end function

  ! Add two instances of type lorenz
  function add_lorenz(lhs,rhs) result(sum)
    class(lorenz)    ,intent(in)  :: lhs
    class(integrand) ,intent(in)  :: rhs
    class(integrand) ,allocatable :: sum
    type(lorenz)     ,allocatable :: local_sum
    select type(rhs)
      class is (lorenz)
        allocate(local_sum)
        call local_sum%set_quadrature(lhs%get_quadrature())
        ! local_sum%state is automatically allocated by assignment
        local_sum%state = lhs%state + rhs%state
        local_sum%sigma = lhs%sigma + rhs%sigma
        local_sum%rho   = lhs%rho   + rhs%rho  
        local_sum%beta  = lhs%beta  + rhs%beta 
      class default
        stop 'assign_lorenz: unsupported class'
    end select
    ! no additional allocation needed by using move_alloc
    call move_alloc(local_sum,sum)
  end function

  ! Multiply an instance of lorenz by a real scalar
  function multiply_lorenz(lhs,rhs) result(product)
    class(lorenz)    ,intent(in)  :: lhs
    real             ,intent(in)  :: rhs
    class(integrand) ,allocatable :: product
    type(lorenz)     ,allocatable :: local_product
    allocate(local_product)
    call local_product%set_quadrature(lhs%get_quadrature())
    ! local_product%state is automatically allocated by assignment
    local_product%state = lhs%state* rhs
    local_product%sigma = lhs%sigma* rhs
    local_product%rho   = lhs%rho  * rhs
    local_product%beta  = lhs%beta * rhs
    ! avoid unnecessary memory allocation by using move_alloc
    call move_alloc(local_product,product)
  end function

  ! Assign one instance to another
  subroutine assign_lorenz(lhs,rhs)
    class(lorenz)    ,intent(inout) :: lhs
    class(integrand) ,intent(in)    :: rhs
    select type(rhs)
      class is (lorenz)
        ! let assignment automatically allocate lhs%state
        lhs%state = rhs%state
        lhs%sigma = rhs%sigma
        lhs%rho   = rhs%rho  
        lhs%beta  = rhs%beta 
        call lhs%set_quadrature(rhs%get_quadrature())
      class default
        stop 'assign_lorenz: unsupported class'
    end select
  end subroutine

  ! Accessor: return state
  function output(this) result(coordinates)
    class(lorenz)      ,intent(in)  :: this
    real ,dimension(:) ,allocatable :: coordinates
    ! assignment allocates coordinates automatically
    coordinates = this%state
  end function 
end module lorenz_module
