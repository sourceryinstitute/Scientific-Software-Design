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

module integrand_module
  implicit none       ! Prevent implicit typing
  private             ! Hide everything by default
  public :: integrate ! expose time integration procedure

  type ,abstract ,public :: integrand
  contains
    procedure(time_derivative      ) ,deferred :: t
    procedure( symmetric_operator  ) ,deferred :: add
    procedure( symmetric_assignment) ,deferred :: assign
    procedure(asymmetric_operator  ) ,deferred :: multiply

    generic :: operator(+) => add ! Map operators to proceures
    generic :: operator(*) => multiply
    generic :: assignment(=) => assign
  end type

  abstract interface
    function time_derivative(this) result(dState_dt)
      import :: integrand 
      class(integrand) ,intent(in)  :: this
      class(integrand) ,allocatable :: dState_dt 
    end function time_derivative
    function symmetric_operator(lhs,rhs) result(operator_result)
      import :: integrand 
      class(integrand) ,intent(in)  :: lhs,rhs
      class(integrand) ,allocatable :: operator_result
    end function symmetric_operator
    function asymmetric_operator(lhs,rhs) result(operator_result)
      import :: integrand 
      class(integrand) ,intent(in)  :: lhs
      class(integrand) ,allocatable :: operator_result
      real                    ,intent(in)  :: rhs
    end function asymmetric_operator
    subroutine symmetric_assignment(lhs,rhs) 
      import :: integrand 
      class(integrand) ,intent(in)    :: rhs
      class(integrand) ,intent(inout) :: lhs
    end subroutine symmetric_assignment
  end interface

contains
  subroutine integrate(model,dt)    
    class(integrand) :: model
    real ,intent(in) :: dt       ! time step size 
    model = model + model%t()*dt ! Explicit Euler formula
  end subroutine 
end module integrand_module
