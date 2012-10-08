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

module cloud_module
  use global_parameters_module ,only :debugging !print call tree if true

  implicit none  ! Prevent implicit typing
  private        ! Hide everything by default
  public :: cloud,cloud_

  ! This type tracks the evolution of the second state variable in the
  ! Lorenz system according to the second Lorenz equation.  It also 
  ! tracks the corresponding paramater (rho) according to the
  ! differential equation d(rho)/dt=0. For illustrative purposes, this
  ! implementation does not expose the number of state variables (2) to
  ! the puppeteer because no iteration is required and the need for
  ! arithmetic operations on rho is therefore an internal concern.  The
  ! rank of the matrix d_dCloud() returns is thus 1 to reflect the only
  ! variable on which the puppeteer needs to iterate when handling
  ! nonlinear couplings in implicit solvers.

  ! Number of evolution equations/variables exposed to the outside
  ! world (via Jacobian sub-block shapes):
  integer ,parameter :: num_eqs=1,num_vars=1 

  type cloud
    private 
    real :: y,rho ! 2nd Lorenz equation solution variable and parameter
  contains
    procedure :: set_coordinate ! accessor: set phase-space coordinate
    procedure :: coordinate! accessor: return phase-space coordinate
    procedure :: t         ! time derivative
    procedure :: d_dCloud  ! contribution to diagonal Jacobian element
    procedure :: d_dx  ! contribution to off-diagonal Jacobian element
    procedure :: d_dz  ! contribution to off-diagonal Jacobian element
    procedure ,private :: add      ! add two instances
    procedure ,private :: subtract ! subtract one instance from another
    procedure ,private :: multiply ! multiply instance by a real scalar
    ! map defined operators to corresponding procedures
    generic   ,public  :: operator(+) => add
    generic   ,public  :: operator(-) => subtract
    generic   ,public  :: operator(*) => multiply
  end type cloud

  interface cloud_
    module procedure constructor
  end interface

contains
  ! constructor: allocate and initialize components
  type(cloud) function constructor(y_initial,r)
    real ,intent(in)  :: y_initial
    real ,intent(in)  :: r
    if (debugging) print *,'      cloud: start'
    constructor%y = y_initial
    constructor%rho = r
    if (debugging) print *,'      cloud: end'
  end function

  ! accessor: set component
  subroutine set_coordinate(this,y_update)
    class(cloud) ,intent(inout) :: this
    real ,intent(in) :: y_update
    this%y = y_update
  end subroutine

  ! accessor (returns phase-space coordinate)
  function coordinate(this) result(return_y)
    class(cloud)       ,intent(in)  :: this
    real ,dimension(:) ,allocatable :: return_y
    return_y = [this%y]
  end function

  ! contribution to diagonal Jacobian element
  function d_dCloud(this,x_ignored,z_ignored) result(dRHS_dy)
    class(cloud)                      ,intent(in) :: this
    real ,dimension(:)   ,allocatable ,intent(in) :: x_ignored,z_ignored
    real ,dimension(:,:) ,allocatable             :: dRHS_dy
    if (debugging) print *,'      cloud%d_dCloud: start'
    allocate(dRHS_dy(num_eqs,num_vars))
   !dRHS_dy(1) = [ d{x(1)*(rho-z(1))-y(1)}/dy(1) ]
    dRHS_dy(1,1) = -1.
    if (debugging) print *,'      cloud%d_dCloud: end'
  end function

  ! contribution to off-diagonal Jacobian element
  function d_dx(this,x,z) result(dRHS_dx)
    class(cloud)                      ,intent(in) :: this
    real ,dimension(:)   ,allocatable ,intent(in) :: x,z
    real ,dimension(:,:) ,allocatable             :: dRHS_dx
    if (debugging) print *,'      cloud%d_dx: start'
    allocate(dRHS_dx(num_eqs,size(x)))
   !dRHS_dx = [ d{x(1)*(rho-z(1))-y}/dx(1)  0  ... 0]
    dRHS_dx = 0.
    dRHS_dx(1,1) = this%rho-z(1)
    if (debugging) print *,'      cloud%d_dx: end'
  end function

  ! contribution to off-diagonal Jacobian element
  function d_dz(this,x,z) result(dRHS_dz)
    class(cloud)                       ,intent(in) :: this
    real ,dimension(:)   ,allocatable ,intent(in) :: x,z
    real ,dimension(:,:) ,allocatable             :: dRHS_dz
    if (debugging) print *,'      cloud%d_dz: start'
    allocate(dRHS_dz(num_eqs,size(z)))
   !dRHS_dz = [ d{x(1)*(rho-z(1))-y(1)}/dz(1)  0  ... 0]
    dRHS_dz = 0.
    dRHS_dz(1,1) = -x(1)
    if (debugging) print *,'      cloud%d_dz: end'
  end function

  ! time derivative (evolution equations)
  function t(this,x,z) result(dy_dt)
    class(cloud)       ,intent(in) :: this
    real ,dimension(:) ,intent(in) :: x,z
    type(cloud)                    :: dy_dt 
    if (debugging) print *,'      cloud%t: start'
    dy_dt%y = x(1)*(this%rho-z(1))-this%y
    dy_dt%rho = 0.
    if (debugging) print *,'      cloud%t: end'
  end function

  function add(lhs,rhs) result(sum) ! add two instances
    class(cloud) ,intent(in) :: lhs,rhs
    type(cloud)              :: sum
    if (debugging) print *,'      cloud%add: start'
    sum%y   = lhs%y   + rhs%y
    sum%rho = lhs%rho + rhs%rho
    if (debugging) print *,'      cloud%add: end'
  end function 

  ! subtract one instance from another
  function subtract(lhs,rhs) result(difference)
    class(cloud) ,intent(in) :: lhs,rhs
    type(cloud)              :: difference
    if (debugging) print *,'      cloud%subtract: start'
    difference%y   = lhs%y   - rhs%y
    difference%rho = lhs%rho - rhs%rho
    if (debugging) print *,'      cloud%subtract: end'
  end function

  ! multiply an instance by a real scalar
  function multiply(lhs,rhs) result(product)
    class(cloud) ,intent(in) :: lhs
    real         ,intent(in) :: rhs
    type(cloud)              :: product
    if (debugging) print *,'      cloud%multiply: start'
    product%y   = lhs%y* rhs
    product%rho = lhs%rho* rhs
    if (debugging) print *,'      cloud%multiply: end'
  end function
end module cloud_module
