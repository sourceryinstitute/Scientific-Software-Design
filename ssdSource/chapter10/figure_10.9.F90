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

subroutine precondition(this,constructor)
  use assertion_utility ,only : assert
  implicit none
  type(field) ,intent(in) :: this
  logical ,intent(in) :: constructor
  logical :: both_allocated ,at_least_one_allocated
  both_allocated= allocated(this%fourier) .and. allocated(this%physical)
  call assert( [.not. both_allocated] &
  ,[error_message('redundant or inconsistent argument')])
  at_least_one_allocated = &
  allocated(this%fourier) .or. allocated(this%physical)
  if (constructor) then
    call assert( [.not. at_least_one_allocated] &
    ,[error_message('constructor argument pre-allocated')])
  else; call assert( [at_least_one_allocated] &
    ,[error_message('argument data missing')])
  end if
end subroutine  

subroutine postcondition(this,public_operator,deletable)
  use assertion_utility ,only : assert
  implicit none
  type(field) ,intent(in) :: this
  logical ,intent(in) :: public_operator,deletable
  logical :: both_allocated,at_least_one_allocated
  both_allocated= allocated(this%fourier) .and. allocated(this%physical)
  call assert( [.not. both_allocated] &
  ,[error_message('redundant or inconsistent result')])
  at_least_one_allocated = &
  allocated(this%fourier) .or. allocated(this%physical)
  if (public_operator .and. deletable) then
     if (this%temporary) then
       call assert([.not. at_least_one_allocated] &
       ,[error_message('temporary argument persisting']))
     else; call assert([at_least_one_allocated] &
       ,[error_message('persistent argument deleted')])
     end if
  else 
     call assert([at_least_one_allocated] &
     ,[error_message('invalid result']))
  end if
end subroutine
