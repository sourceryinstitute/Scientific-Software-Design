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

module linear_solve_module
  use kind_parameters ,only : rkind,ikind
  implicit none
contains
  function gaussian_elimination(lhs,rhs) result(x)
    real(rkind) ,dimension(:,:) ,intent(in) :: lhs
    real(rkind) ,dimension(:) ,allocatable ,intent(in) :: rhs
    real(rkind) ,dimension(:)   ,allocatable :: x,b ! Linear system:
    real(rkind) ,dimension(:,:) ,allocatable :: A   ! Ax = b
    real(rkind)            :: factor
    real(rkind) ,parameter :: pivot_tolerance=1.0E-02
    integer(ikind)         :: row,col,n,p ! p=pivot row/col
    n=size(lhs,1)
    b = rhs ! Copy rhs side to preserve required intent
    if ( n /= size(lhs,2) .or. n /= size(b)) &
      stop 'gaussian_elimination: ill-posed system'
    allocate(x(n))
    A = lhs            ! Copy lhs side to preserve required intent
    !______ Gaussian elimination _______________________________________
    do p=1,n-1         ! Forward elimination
      if (abs(A(p,p))<pivot_tolerance) &
        stop 'gaussian_elimination: use pivoting'
      do row=p+1,n
        factor=A(row,p)/A(p,p)
        forall(col=p:n) A(row,col) = A(row,col) - A(p,col)*factor
        b(row) = b(row) - b(p)*factor
      end do
    end do
    x(n) = b(n)/A(n,n) ! Back substitution 
    do row=n-1,1,-1
      x(row) = (b(row) - sum(A(row,row+1:n)*x(row+1:n)))/A(row,row)
    end do
  end function
end module
