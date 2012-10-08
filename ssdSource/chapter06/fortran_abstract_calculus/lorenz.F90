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
  use integrand_module ,only : integrand,integrate
  implicit none

  private             ! Hide everything by default
  public :: integrate ! Expose time integration procedure
  public :: lorenz

  type ,extends(integrand) :: lorenz
    private 
    real ,dimension(:) ,allocatable :: state   ! solution vector
    real :: sigma ,rho ,beta                   ! Lorenz parameters
  contains
     procedure ,public :: t        => dLorenz_dt
     procedure ,public :: add      => add_lorenz
     procedure ,public :: multiply => multiply_lorenz
     procedure ,public :: assign   => assign_lorenz
     procedure ,public :: output
  end type

  interface lorenz
    procedure constructor
  end interface

contains

  type(lorenz) function constructor(initial_state,s,r,b)
    real ,dimension(:) ,intent(in)  :: initial_state 
    real               ,intent(in)  :: s ,r ,b  ! passed values for
                                                ! sigma, rho and beta
    constructor%state=initial_state
    constructor%sigma=s
    constructor%rho=r
    constructor%beta=b
  end function

  function output(this) result(coordinates)
    class(lorenz)      ,intent(in)  :: this
    real ,dimension(:) ,allocatable :: coordinates
    coordinates = this%state
  end function output

  ! time derivative: encapsulates Lorenz equations
  function dLorenz_dt(this) result(dState_dt)
    class(lorenz)           ,intent(in)  :: this
    class(integrand) ,allocatable :: dState_dt 
    type(lorenz)            ,allocatable :: delta

    allocate(delta) 
    allocate(delta%state(size(this%state)))
    ! 1st lorenz equation
    delta%state(1)=this%sigma*( this%state(2) -this%state(1))
    ! 2nd lorenz equation
    delta%state(2)=this%state(1)*(this%rho-this%state(3))-this%state(2)
    ! 3rd lorenz equation
    delta%state(3)=this%state(1)*this%state(2)-this%beta*this%state(3)
    ! hold Lorenz parameters constant over time
    delta%sigma=0.
    delta%rho=0.
    delta%beta=0.
    call move_alloc (delta, dState_dt)
  end function

  function add_Lorenz(lhs,rhs) result(sum) ! add two Lorenz objects
    class(lorenz)           ,intent(in)  :: lhs
    class(integrand) ,intent(in)  :: rhs
    class(integrand) ,allocatable :: sum
    type(lorenz)            ,allocatable :: local_sum

    allocate (lorenz :: local_sum)
    select type(rhs)
      class is (lorenz)
        local_sum%state = lhs%state + rhs%state
        local_sum%sigma = lhs%sigma + rhs%sigma
        local_sum%rho   = lhs%rho   + rhs%rho
        local_sum%beta  = lhs%beta  + rhs%beta
      class default
        stop 'add_Lorenz: rhs argument type not supported'
    end select
    call move_alloc(local_sum, sum)
  end function

  ! multiply a Lorenz object by a real scalar
  function multiply_Lorenz(lhs,rhs) result(product)
    class(lorenz) ,intent(in)  :: lhs
    real          ,intent(in)  :: rhs
    class(integrand) ,allocatable :: product
    type(lorenz)            ,allocatable :: local_product

    allocate (local_product)
    local_product%state = lhs%state*rhs
    local_product%sigma = lhs%sigma*rhs
    local_product%rho   = lhs%rho  *rhs
    local_product%beta  = lhs%beta *rhs
    call move_alloc(local_product, product)
  end function

  ! assign one lorenz object to another
  subroutine assign_lorenz(lhs,rhs)
    class(lorenz)           ,intent(inout) :: lhs
    class(integrand) ,intent(in)    :: rhs

    select type(rhs)
      class is (lorenz)
        lhs%state = rhs%state 
        lhs%sigma = rhs%sigma
        lhs%rho   = rhs%rho  
        lhs%beta  = rhs%beta
      class default
        stop 'assign_lorenz: rhs argument type not supported'
    end select
  end subroutine
end module lorenz_module
