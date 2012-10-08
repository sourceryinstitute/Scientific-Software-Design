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

module conductor_class
  use kind_parameters ,only : ikind,rkind
  use differentiator_class, only : differentiator
  implicit none
  private           ! Hide everything by default
  public::conductor,conductor_ ! Expose type and type-bound procedures
  integer(ikind), parameter :: end_points=2 
  type conductor
    private
    type(differentiator) :: diff
    real(rkind) :: diffusivity
    real(rkind) ,dimension(end_points) :: boundary
    real(rkind) ,dimension(:), allocatable :: temp 
  contains                                   !^ internal temperatures
    procedure :: heat_for
    procedure :: temperature
    procedure :: time_derivative
  end type
  interface conductor_
    module procedure constructor
  end interface
contains
  type(conductor) function constructor(diff,prob,T_init) 
    use conduction_module ,only : read_physics
    use problem_class     ,only : problem
    type(differentiator) ,intent(in)    :: diff
    type(problem)        ,intent(in)    :: prob
    real(rkind), dimension(:) ,optional :: T_init
    constructor%diff = diff
    constructor%diffusivity = prob%diffusivity()
    constructor%boundary = prob%boundary_vals()
    if (present(T_init)) then
      if (size(T_init)/=prob%nodes())stop 'In conductor: size mismatch.'
      constructor%temp = T_init
    else; allocate(constructor%temp(prob%nodes()))
          constructor%temp = constructor%boundary(1)
    end if
  end function
  function temperature(this)
    class(conductor), intent(in) :: this
    real(rkind), dimension(:), allocatable :: temperature
    temperature = (/ this%boundary(1), this%temp, this%boundary(2) /)
  end function 
  subroutine heat_for(this,dt)
    class(conductor) ,intent(inout) :: this
    real(rkind)      ,intent(in)    :: dt
    real(rkind) ,dimension(:) ,allocatable :: T_xx !2nd derivative
    allocate(T_xx(size(this%temp)))
    T_xx = this%diff%laplacian(this%temp,this%boundary)
    this%temp = this%temp + dt*this%diffusivity*T_xx
  end subroutine
  real(rkind) function time_derivative(this)
    class(conductor)   :: this
    real(rkind) ,dimension(:) ,allocatable :: T_xx_end
    real(rkind)                            :: T_2
    if (size(this%temp)>1) then 
      T_2=this%temp(2)
    else; T_2=this%boundary(2)
    end if
    T_xx_end = this%diff%laplacian(this%temp,(/this%boundary(1),T_2/) )
    time_derivative = this%diffusivity*T_xx_end(1)
  end function
end module
