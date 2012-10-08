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

module field_module
  use kind_parameters ,only : rkind, ikind
  implicit none
  private
  public :: field
  type :: field 
    real(rkind), allocatable :: f(:)
  contains
    procedure :: add            => total
    procedure :: subtract       => difference      
    procedure :: multiply_field => product_
    procedure :: multiply_real  => multiple
    procedure :: assign         => assign_field_f
    procedure :: copy           => copy_filed
    procedure :: state          => field_values
    generic   :: operator(+)    => add
    generic   :: operator(-)    => subtract
    generic   :: operator(*)    => multiply_real,multiply_field
    generic   :: assignment(=)  => assign, copy
  end type
    
contains

  function field_values (this)
    class(field), intent(in) :: this
    real(rkind), allocatable :: field_values(:)
    field_values = this%f
  end function

  subroutine assign_field_f (lhs, rhs)
    class(field), intent(inout) :: lhs
    real(rkind), intent(in) :: rhs(:)
    lhs%f = rhs
  end subroutine

  subroutine copy_filed (lhs, rhs)
    class(field), intent(inout) :: lhs
    class(field), intent(in) :: rhs
    lhs%f = rhs%f
  end subroutine

  function total(lhs,rhs)
    class(field) ,intent(in) :: lhs
    class(field) ,intent(in) :: rhs
    class(field) ,allocatable :: total
    allocate (total)
    total%f = lhs%f + rhs%f
  end function

  function difference(lhs,rhs)
    class(field) ,intent(in) :: lhs
    class(field) ,intent(in)  :: rhs
    class(field) ,allocatable :: difference
    allocate (difference)
    difference%f = lhs%f - rhs%f
  end function

  function product_(lhs,rhs)
    class(field) ,intent(in) :: lhs
    class(field) ,intent(in)  :: rhs
    class(field) ,allocatable :: product_
    allocate(product_)
    product_%f = lhs%f * rhs%f
  end function

  function multiple(lhs,rhs)
    class(field) ,intent(in) :: lhs
    real(rkind) ,intent(in)  :: rhs
    class(field) ,allocatable :: multiple
    allocate(multiple)
    multiple%f = lhs%f * rhs
  end function
end module
