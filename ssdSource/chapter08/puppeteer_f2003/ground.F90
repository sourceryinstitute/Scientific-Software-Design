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

module ground_module
  use global_parameters_module ,only :debugging !print call tree if true

  implicit none ! Prevent implicit typing
  private       ! Hide everything by default
  public :: ground,ground_

  ! This type tracks the evolution of the third state variable in the
  ! Lorenz system according to the third Lorenz equation.  It also
  ! tracks the corresponding paramater (beta) according to the
  ! differential equation d(beta)/dt=0. For illustrative purposes,
  ! this implementation does not expose the number of state variables(2)
  ! to the puppeteer because no iteration is required and the need for
  ! arithmetic operations on beta is therefore an internal concern.  The
  ! rank of the matrix d_dGround() returns is thus 1 to reflect the
  ! only variable on which the puppeteer needs to iterate when handling
  ! nonlinear couplings in implicit solvers.

  ! Number of evolution equations/variables exposed to the outside
  ! world (via Jacobian sub-block shapes):
  integer ,parameter :: num_eqs=1,num_vars=1

  type ground
    private 
    real :: z,beta ! 3rd Lorenz equation solution variable and parameter
  contains
    procedure :: set_coordinate ! accessor set phase-space coordinate
    procedure :: coordinate ! accessor: return phase-space coordinate
    procedure :: t          ! time derivative
    procedure :: d_dGround  ! contribution to diagonal Jacobian element
    procedure :: d_dx   ! contribution to off-diagonal Jacobian element
    procedure :: d_dy   ! contribution to off-diagonal Jacobian element
    procedure ,private :: add ! add two instances
    procedure ,private :: subtract ! subtract one instance from another
    procedure ,private :: multiply ! multiply instance by a real scalar
    ! map defined operators to corresponding procedures
    generic ,public  :: operator(+) => add
    generic ,public  :: operator(-) => subtract
    generic ,public  :: operator(*) => multiply
  end type

  interface ground_
    module procedure constructor
  end interface

contains

  ! constructor: allocate and initialize components
  type(ground) function constructor(z_initial,b)
    real ,intent(in) :: z_initial
    real ,intent(in) :: b
    if (debugging) print *,'      ground%construct: start'
    constructor%z = z_initial
    constructor%beta = b
    if (debugging) print *,'      ground%construct: end'
  end function

  ! accessor (returns phase-space coordinate)
  function coordinate(this) result(return_z)
    class(ground)       ,intent(in)  :: this
    real ,dimension(:) ,allocatable :: return_z
    return_z = [ this%z ]
  end function

  ! accessor (returns phase-space coordinate)
  subroutine set_coordinate(this,z_update) 
    class(ground) ,intent(inout)  :: this
    real ,intent(in) :: z_update
    this%z = z_update 
  end subroutine

  ! time derivative (evolution equations)
  function t(this,x,y) result(dz_dt)
    class(ground)      ,intent(in) :: this
    real ,dimension(:) ,intent(in) :: x,y
    type(ground)                   :: dz_dt 
    if (debugging) print *,'      ground%t: start'
    dz_dt%z = x(1)*y(1) - this%beta*this%z
    dz_dt%beta = 0.
    if (debugging) print *,'      ground%t: end'
  end function

  ! contribution to diagonal Jacobian element
  function d_dGround(this,x_ignored,y_ignored) result(dRHS_dz)
    class(ground)                     ,intent(in) :: this
    real ,dimension(:)   ,allocatable ,intent(in) :: x_ignored,y_ignored
    real ,dimension(:,:) ,allocatable             :: dRHS_dz
    if (debugging) print *,'      ground%d_dGround: start'
   !dRHS_dz = [ d{x(1)*y(1) - beta*z}/dz(1) ]
    allocate(dRHS_dz(num_eqs,num_vars))
    dRHS_dz(1,1) = -this%beta
    if (debugging) print *,'      ground%d_dGround: end'
  end function

  ! contribution to off-diagonal Jacobian element
  function d_dx(this,x,y) result(dRHS_dx)
    class(ground)                     ,intent(in) :: this
    real ,dimension(:)   ,allocatable ,intent(in) :: x,y
    real ,dimension(:,:) ,allocatable             :: dRHS_dx
    if (debugging) print *,'      ground%d_dx: start'
    allocate(dRHS_dx(num_eqs,size(x)))
   !dRHS_dz = [ d{x(1)*y(1) - beta*z(1)}/dx(1)   0  ...  0 ]
    dRHS_dx=0.
    dRHS_dx(1,1) = y(1)
    if (debugging) print *,'      ground%d_dx: end'
  end function

  ! contribution to off-diagonal Jacobian element
  function d_dy(this,x,y) result(dRHS_dy)
    class(ground)                     ,intent(in) :: this
    real ,dimension(:)   ,allocatable ,intent(in) :: x,y
    real ,dimension(:,:) ,allocatable             :: dRHS_dy
    if (debugging) print *,'      ground%d_dy: start'
    allocate(dRHS_dy(num_eqs,size(y)))
   !dRHS_dz = [ d{x(1)*y(1) - beta*z(1)}/dy(1)   0  ...  0 ]
    dRHS_dy = 0.
    dRHS_dy(1,1) = x(1)
    if (debugging) print *,'      ground%d_dy: end'
  end function

  function add(lhs,rhs) result(sum) ! add two instances
    class(ground) ,intent(in) :: lhs,rhs
    type(ground)              :: sum
    if (debugging) print *,'      ground%add: start'
    sum%z    = lhs%z    + rhs%z
    sum%beta = lhs%beta + rhs%beta
    if (debugging) print *,'      ground%add: end'
  end function

  ! subtract one instance from another
  function subtract(lhs,rhs) result(difference)
    class(ground) ,intent(in) :: lhs,rhs
    type(ground)              :: difference
    if (debugging) print *,'      ground%subtract: start'
    difference%z    = lhs%z    - rhs%z
    difference%beta = lhs%beta - rhs%beta
    if (debugging) print *,'      ground%subtract: end'
  end function 

  ! multiply an instance by a real scalar
  function multiply(lhs,rhs) result(product)
    class(ground) ,intent(in) :: lhs
    real      ,intent(in)  :: rhs
    type(ground)              :: product
    if (debugging) print *,'      ground%multiply: start'
    product%z    = lhs%z   * rhs
    product%beta = lhs%beta* rhs
    if (debugging) print *,'      ground%multiply: end'
  end function 
end module ground_module
