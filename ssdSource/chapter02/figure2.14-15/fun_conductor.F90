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

module fun_conductor_class
  use kind_parameters ,only : ikind,rkind
  use conductor_class ,only : conductor,conductor_
  implicit none
  private               
  public::fun_conductor,fun_conductor_
  type ,extends(conductor) :: fun_conductor
    private
    real(rkind) :: dt
  contains
    procedure :: time_step
  end type
  interface fun_conductor_
    module procedure constructor
  end interface 
contains
   type(fun_conductor) function constructor(file,T_distribution) 
    use differentiator_class ,only : differentiator,differentiator_
    use problem_class        ,only : problem,problem_
    integer(ikind) ,intent(in)             :: file
    type(differentiator)                   :: diff
    type(problem)                          :: prob
    real(rkind), dimension(:) ,allocatable :: T
    real(rkind)                            :: dx
    integer(ikind)                         :: i
    interface
      pure real(rkind) function T_distribution(x)
        use kind_parameters ,only : rkind
        real(rkind) ,intent(in) :: x
      end function
    end interface
    prob = problem_(file)
    diff = differentiator_(prob)
    dx   = prob%spacing()
    constructor%dt = prob%time_step()
    allocate(T(prob%nodes()))
    forall(i=1:prob%nodes()) T(i) = T_distribution(i*dx)
    constructor%conductor = conductor_(diff,prob,T)
  end function
  real(rkind) function time_step(this) 
    class(fun_conductor) :: this
    time_step = this%dt
  end function 
end module
