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

module differentiator_class
  use kind_parameters ,only: ikind,rkind
  implicit none
  private                  ! Hide everything by default.
  public :: differentiator,differentiator_ ! Expose type/constructor/type-bound procs.
  type differentiator
    private
    real(rkind),dimension(:,:),allocatable::diff_matrix
  contains
    procedure :: laplacian  ! return Laplacian
    procedure :: lap_matrix ! return Laplacian matrix operator
  end type
  interface differentiator_
    module procedure constructor
  end interface 
contains
  type(differentiator) function constructor(spec) 
    use problem_class ,only : problem
    use conduction_module ,only : differencing_stencil
    type(problem) ,intent(in) :: spec
    integer(ikind) ,parameter :: diagonals=3
    allocate(constructor%diff_matrix(spec%nodes(),diagonals))
    constructor%diff_matrix = &
      differencing_stencil(spec%spacing(),spec%nodes())
  end function
  pure function laplacian(this,T,Tboundaries)
    use conduction_module, only : differentiate
    class(differentiator), intent(in) :: this
    real(rkind) ,dimension(:) ,intent(in) :: T
    real(rkind) ,dimension(:) ,allocatable :: laplacian
    integer(ikind) ,parameter :: end_points=2
    real(rkind), dimension(end_points), intent(in) :: Tboundaries
    allocate(laplacian(size(T)))
    laplacian = differentiate(& !compute derivative at internal points
      this%diff_matrix,T,Tboundaries(1),Tboundaries(2) )
  end function
  pure function lap_matrix(this)
    class(differentiator), intent(in) :: this
    real(rkind) ,allocatable ,dimension(:,:) :: lap_matrix
    lap_matrix = this%diff_matrix
  end function
end module 
