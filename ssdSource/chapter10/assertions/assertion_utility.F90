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

module assertion_utility
  use iso_fortran_env ,only : error_unit  
  implicit none
  private
  public :: error_message,assert,assert_identical
  type error_message
    character(:) ,allocatable :: string
  end type
contains
  subroutine assert(assertion,text)
    logical ,dimension(:) ,intent(in) :: assertion
    type(error_message) ,dimension(:) ,intent(in) :: text
    integer :: i
    logical :: any_failures 
    call assert_identical( [size(assertion),size(text)] )
    any_failures=.false.
    do i=1,size(assertion)
      if (.not. assertion(i)) then
        any_failures=.true.
        write(error_unit,*) 'Assertion failed with message: '
        if (allocated(text(i)%string)) then
          write(error_unit,*) text(i)%string
        else
          write(error_unit,*) '(no message provided).'
        end if
      end if
    end do
    if (any_failures) stop 'Execution halted on failed assertion(s)!'
  end subroutine
  subroutine assert_identical(integers)
    integer ,dimension(:) ,intent(in) :: integers
    integer :: i
    logical :: any_mismatches
    any_mismatches = .false.
    do i=2,size(integers)
      if (integers(i) /= integers(1)) then
        any_mismatches = .true.
        write(error_unit,*) &
        'Value ',i,' does not match expected value ',integers(1)
      end if
    end do
    if (any_mismatches) stop 'Execution halted on failed assertion!'
  end subroutine
end module
