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

module problem_class
  use kind_parameters ,only: ikind,rkind
  implicit none
  private
  public :: problem,problem_
  type problem
    private  ! Default values:
    integer(ikind) :: num_nodes=4
    real(rkind):: air_temp=1.,chip_temp=1.!boundary temperatures
    real(rkind):: alpha=1.,length=1.      !physical parameters
    real(rkind):: dt=0.1,dx               !numerical parameters 
  contains              !^ constructor computes default 
    procedure :: boundary_vals
    procedure :: diffusivity
    procedure :: nodes
    procedure :: time_step
    procedure :: spacing 
  end type
  interface problem_
    module procedure spec ! constructor
  end interface
contains
  type(problem) function spec(file)
    use conduction_module &
      ,only : read_physics,read_numerics,stable_dt
    integer(ikind) :: file
    integer(ikind) :: elements_default

    elements_default = spec%num_nodes+1 
    spec%dx = spec%length/elements_default ! default element size

    if (.not. read_physics(file                              &
      ,spec%alpha,spec%length,spec%chip_temp,spec%air_temp)) &
      print *,'In problem constructor: Using default physics spec.'

    if (.not. read_numerics(file                    &
      ,spec%dt,spec%num_nodes,spec%dx,spec%length)) &
      print *,'In problem constructor: Using default numerics spec.'

    spec%dt = min(spec%dt,stable_dt(spec%dx,spec%alpha))
  end function
  pure function boundary_vals(this)
    class(problem) ,intent(in) :: this
    integer(ikind) ,parameter :: end_points=2
    real(rkind) ,dimension(end_points) :: boundary_vals
    boundary_vals = (/this%chip_temp,this%air_temp/)
  end function
  pure real(rkind) function diffusivity(this)
    class(problem) ,intent(in) :: this
    diffusivity = this%alpha
  end function
  pure integer(ikind) function nodes(this)
    class(problem) ,intent(in) :: this
    nodes = this%num_nodes
  end function
  pure real(rkind) function time_step(this)
    class(problem) ,intent(in) :: this
    time_step = this%dt
  end function
  pure real(rkind) function spacing(this)
    class(problem) ,intent(in) :: this
    spacing = this%dx 
  end function
end module
