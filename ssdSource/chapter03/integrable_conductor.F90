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

module integrable_conductor_class
  use kind_parameters ,only : ikind ,rkind
  use field_class     ,only : field, distribution
  implicit none
  private
  public :: integrable_conductor,integrable_conductor_
  type integrable_conductor
    private
    type(field) :: T_field ! temperatures
    real(rkind) :: alpha   ! thermal diffusivity
  contains
    procedure            :: t
    procedure            :: rhs_operator 
    procedure            :: rhs_operator_size
    procedure            :: temperature
    procedure            :: time_derivative
    procedure ,private   :: product
    generic              :: operator(*) => product
    procedure ,private   :: total
    generic              :: operator(+) => total
    procedure ,private ,pass(rhs) :: inverseTimes
    generic              :: operator(.inverseTimes.) => inverseTimes
  end type 
  interface integrable_conductor_
    module procedure constructor
  end interface
contains
  type(integrable_conductor) function constructor(spec,T_distribution)
    use problem_class ,only : problem
    type(problem) :: spec ! problem specification
    procedure(distribution) T_distribution
    constructor%T_field = field(spec,T_distribution)
    constructor%alpha   = spec%diffusivity()
  end function

  pure integer(ikind) function rhs_operator_size(this)
    class(integrable_conductor) ,intent(in) :: this
    rhs_operator_size = this%T_field%field_size()
  end function

  function rhs_operator(this)
    class(integrable_conductor) ,intent(in)  :: this
    real(rkind) ,dimension(:,:) ,allocatable :: rhs_operator
    rhs_operator = this%T_field%xx_matrix()
  end function

  function temperature(this)
    class(integrable_conductor) ,intent(in)  :: this
    real(rkind) ,dimension(:) ,allocatable :: temperature
    temperature = this%T_field%nodal_values()
  end function

  pure type(integrable_conductor) function product(this,factor)
    class(integrable_conductor) ,intent(in) :: this
    real(rkind)                 ,intent(in) :: factor
    product%T_field = this%T_field * factor
    product%alpha   = this%alpha   * factor
  end function

  pure type(integrable_conductor) function total(lhs,rhs)
    class(integrable_conductor) ,intent(in) :: lhs 
    type(integrable_conductor)  ,intent(in) :: rhs
    total%T_field = lhs%T_field + rhs%T_field
    total%alpha   = lhs%alpha   + rhs%alpha
  end function

  pure function t(this) result(d_dt)
    class(integrable_conductor) ,intent(in) :: this
    type(integrable_conductor)              :: d_dt
    d_dt%T_field = this%T_field%xx()*this%alpha
    d_dt%alpha   = 0.0_rkind
  end function

  pure real(rkind) function time_derivative(this) 
    class(integrable_conductor) ,intent(in) :: this
    time_derivative = this%alpha*this%T_field%xx_boundary()
  end function

  type(integrable_conductor) function inverseTimes(lhs,rhs)
    real(rkind) ,dimension(:,:) ,intent(in) :: lhs
    class(integrable_conductor)       ,intent(in) :: rhs
    inverseTimes%T_field = lhs .inverseTimes. rhs%T_field
    inverseTimes%alpha = rhs%alpha
  end function
end module integrable_conductor_class
