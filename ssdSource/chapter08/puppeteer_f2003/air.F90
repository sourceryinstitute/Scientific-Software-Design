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

module air_module
  use global_parameters_module ,only : debugging !print call treeif true
  implicit none ! Prevent implicit typing
  private       ! Hide everything by default
  public :: air,air_

  ! Number of evolution equations/variables exposed to the outside world
  ! (via Jacobian sub-block shapes):
  integer ,parameter :: num_eqs=2,num_vars=2 

  ! This type tracks the evolution of the first state variable in the
  ! Lorenz system according to the first Lorenz equation.  It also
  ! tracks the corresponding paramater (sigma) according to the
  ! differential equation d(sigma)/dt=0. For illustrative purposes, this
  ! implementation exposes the number of state variables (2) to the
  ! puppeteer without providing direct access to them or exposing
  ! anything about their layout, storage location or identifiers (x and
  ! sigma).  Their existence is apparent in the rank (2) of the matrix 
  ! d_dAir() returns as its diagonal Jacobian element contribution. 

  type air
    private 
    real :: x,sigma !1st Lorenz equation solution variable and parameter
  contains
    procedure :: coordinate! accessor: return phase-space coordinates
    procedure :: t ! time derivative ( evaluate RHS of Lorenz equation)
    procedure :: d_dAir ! contribution to diagonal Jacobian element
    procedure :: d_dy   ! contribution to off-diagonal Jacobian element
    procedure ,private :: add ! add two instances
    procedure ,private :: subtract ! subtract one instance from another
    procedure ,private :: multiply ! multiply instance by a real scalar
    ! map defined operators to corresponding procedures
    generic   :: operator(+) => add
    generic   :: operator(-) => subtract
    generic   :: operator(*) => multiply
  end type 
  
  interface air_
    module procedure constructor
  end interface

contains
  ! constructor: allocate and initialize components
  type(air) function constructor(x_initial,s)
    real           ,intent(in)  :: x_initial
    real           ,intent(in)  :: s
    if (debugging) print *,'      air%construct: start'
    constructor%x = x_initial
    constructor%sigma  = s
    if (debugging) print *,'      air%construct: end'
  end function

  ! accessor: returns phase-space coordinates
  function coordinate(this) result(return_x)
    class(air)          ,intent(in)  :: this
    real  ,dimension(:) ,allocatable :: return_x
    return_x = [ this%x ,this%sigma ]
  end function

  ! time derivative (evolution equations)
  function t(this,y) result(dx_dt)
    class(air)         ,intent(in) :: this
    real ,dimension(:) ,intent(in) :: y
    type(air)                      :: dx_dt 
    if (debugging) print *,'      air%t: start'
    dx_dt%x=this%sigma*(y(1)-this%x)
    dx_dt%sigma=0.
    if (debugging) print *,'      air%t: end'
  end function

  ! contribution to diagonal Jacobian element
  function d_dAir(this,y) result(dRHS_dx)
    class(air)  ,intent(in)                :: this
    real        ,dimension(:,:) ,allocatable :: dRHS_dx
    real        ,dimension(:) ,allocatable :: y
    if (debugging) print *,'      air%d_dAir: start'
   !allocate(dRHS_dx(num_eqs,num_vars))
   !dRHS_dx = [ d{sigma*(y-x)}/dx   d{sigma*(y-x)}/dsigma ]
   !          [ d{0}/dx             d{0}/dsigma           ]
    if (size(y) /= 1) stop 'd_dAir: invalid y size'
    dRHS_dx = reshape(source=(/-this%sigma,0.,y(1)-this%x,0./),      &
               shape=(/num_eqs,num_vars/))
    if (debugging) print *,'      air%d_dAir: end'
  end function

  ! contribution to off-diagonal Jacobian element
  function d_dy(this,y) result(dRHS_dy)
    class(air)  ,intent(in)                  :: this
    real        ,dimension(:,:) ,allocatable :: dRHS_dy
    real        ,dimension(:)   ,allocatable :: y

    if (debugging) print *,'      air%d_dy: start'
    allocate(dRHS_dy(num_eqs,size(y))) 
   !dRHS_dy = [ d{sigma*(y(1)-x(1))}/dy(1)  0  ... 0]
   !          [ d{0}/dy(1)                  0  ... 0]
    dRHS_dy = 0.
    dRHS_dy(1,1) = this%sigma
    if (debugging) print *,'      air%d_dy: end'
  end function

  function add(lhs,rhs) result(sum) ! add two instances
    class(air) ,intent(in) :: lhs,rhs
    type(air)              :: sum

    if (debugging) print *,'      air%add: start'
    sum%x     = lhs%x     + rhs%x
    sum%sigma = lhs%sigma + rhs%sigma
    if (debugging) print *,'      air%add: end'
  end function

  ! subtract one instance from another
  function subtract(lhs,rhs) result(difference)
    class(air) ,intent(in) :: lhs,rhs
    type(air)              :: difference
    if (debugging) print *,'      air%subtract: start'
    difference%x     = lhs%x     - rhs%x
    difference%sigma = lhs%sigma - rhs%sigma
    if (debugging) print *,'      air%subtract: end'
  end function

  ! multiply an instance by a real scalar
  function multiply(lhs,rhs) result(product)
    class(air) ,intent(in) :: lhs
    real       ,intent(in) :: rhs
    type(air)              :: product
    if (debugging) print *,'      air%multiply: start'
    product%x     = lhs%x    *rhs
    product%sigma = lhs%sigma*rhs
    if (debugging) print *,'      air%multiply: end'
  end function
end module air_module
