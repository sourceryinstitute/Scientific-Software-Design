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

module vortex_module
  use hermetic_module ,only : hermetic
  implicit none
  private
  public :: vortex,vortex_

  integer ,parameter :: space_dim=3 
  type ,extends(hermetic) :: vortex
    private
    real ,dimension(space_dim) :: point
    type(vortex) ,pointer :: next=>null()
  contains
    procedure :: force_finalization => manual_finalizer
    procedure :: output
    procedure :: assert
    final :: finalize
  end type

  interface vortex_
    module procedure ring
  end interface

contains
  function ring(radius,num_points) 
    real ,intent(in) :: radius 
    integer ,intent(in) :: num_points 
    type(vortex) ,pointer :: ring,current
    integer :: i
    real ,parameter :: pi = acos(-1.)
    real :: delta,theta
    delta = 2.*pi*radius/real(num_points) ! arclength increment
    allocate(ring)
    ring%point =  [radius,0.,0.] 
    current => ring
    do i=1,num_points-1
      allocate(current%next)
      theta = i*delta/radius
      current%next%point = [radius*cos(theta),radius*sin(theta),0.]
      current => current%next 
      if (.not.associated(current)) exit ! not in the original printing
    end do
    current%next => ring ! point tail to head
  end function

  recursive subroutine finalize (this) 
    type(vortex), intent(inout) :: this 
    type(vortex) ,pointer :: next
    next => this%next
    this%next=>null() ! break the chain so recursion terminates
    if (associated(next)) deallocate (next) ! automatic recursion
  end subroutine 

  subroutine manual_finalizer(this)
    class(vortex) ,intent(inout) :: this
    type(vortex) ,pointer :: current,next
    current => this%next
    this%next=>null() ! break the chain so loop terminates
    do while(associated(current%next))
      next => current%next ! start after break in chain
      deallocate(current)
      current=>next
    end do
  end subroutine

  subroutine output(this) 
    class(vortex) ,intent(in) ,target :: this
    type(vortex) ,pointer :: current
    current => this
    do 
      print *,current%point
      current => current%next
      if (associated(current,this)) exit
    end do
  end subroutine

  subroutine assert(this, node_index, check_value)
    class(vortex) ,intent(in) ,target :: this
    integer                   ,intent(in) :: node_index
    real         ,dimension(:),intent(in) :: check_value
    type(vortex) ,pointer    :: current
    real         ,parameter  :: epsilon=1.0e-6
    integer                   :: i         
    current => this
    i=1
    do
      current => current%next 
      i=i+1
      if (i==node_index) then
        if (abs(sum(current%point-check_value)) <=epsilon) then
          print *, 'Test passed'
          exit
        end if
      end if
      if (associated(current,this)) exit
    end do
  end subroutine
end module
