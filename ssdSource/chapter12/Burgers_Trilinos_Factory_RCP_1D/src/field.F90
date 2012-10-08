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
  use iso_c_binding, only: c_int, c_double
  use FEpetra_Comm, only:Epetra_Comm
  implicit none
  private
  public :: field, initial_field
  
  type ,abstract :: field 
  contains
    procedure(field_op_field) ,deferred :: add
    procedure(field_eq_field) ,deferred :: assign         
    procedure(field_op_field) ,deferred :: subtract    
    procedure(field_op_field) ,deferred :: multiply_field 
    procedure(field_op_real)  ,deferred :: multiply_real  
    procedure(real_to_real)  ,deferred :: runge_kutta_2nd_step
    procedure(derivative) ,deferred :: x  ! 1st derivative
    procedure(derivative) ,deferred :: xx ! 2nd derivative
    procedure(output_interface) ,deferred  :: output
    generic :: operator(+)   => add 
    generic :: operator(-)   => subtract
    generic :: operator(*)   => multiply_real,multiply_field
    generic :: assignment(=) => assign
    procedure(force_finalize_interface) ,deferred :: force_finalize
  end type
    
  abstract interface
    real(c_double) pure function initial_field(x)
      import :: c_double
      real(c_double) ,intent(in) :: x
    end function
    function field_op_field(lhs,rhs)
      import :: field
      class(field) ,intent(in) :: lhs,rhs
      class(field) ,allocatable :: field_op_field
    end function
    function field_op_real(lhs,rhs)
      import :: field,c_double
      class(field) ,intent(in)  :: lhs 
      real(c_double) ,intent(in) :: rhs 
      class(field) ,allocatable :: field_op_real
    end function
    real(c_double) function real_to_real(this,nu,grid_resolution)
      import :: field,c_double,c_int
      class(field) ,intent(in)  :: this
      real(c_double) ,intent(in) :: nu
      integer(c_int),intent(in) :: grid_resolution
    end function
    function derivative(this)
      import :: field
      class(field) ,intent(in)  :: this
      class(field) ,allocatable :: derivative
    end function
    subroutine field_eq_field(lhs,rhs)
      import :: field
      class(field) ,intent(inout) :: lhs
      class(field) ,intent(in) :: rhs
    end subroutine
    subroutine output_interface(this,comm)
      import :: field,Epetra_Comm
      class(field) ,intent(in) :: this
      class(Epetra_Comm), intent(in) :: comm
    end subroutine
    subroutine force_finalize_interface(this)
      import :: field
      class(field), intent(inout)  :: this
    end subroutine
  end interface
end module
