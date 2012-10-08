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

  function thomasTimes(lhs,rhs)
    real(rkind) ,dimension(:,:) ,allocatable ,intent(in) :: lhs
    class(field)                             ,intent(in) :: rhs
    type(field)                 ,allocatable :: tomasTimes
    real(rkind) ,dimension(:)   ,allocatable :: x,b ! Linear system:
    real(rkind) ,dimension(:,:) ,allocatable :: A   ! Ax = b
    real(rkind)            :: factor
    integer(ikind)         :: row,n 
    n=size(lhs,1)
    b = rhs%node    ! Copy rhs side to preserve required intent
    if ( n /= size(lhs,2) .or. n /= size(b)) &
      stop 'thomasTimes: ill-posed system'
    allocate(x(n))
    A = lhs            ! Copy lhs side to preserve required intent
    !______ Establish upper triangular matrix _______________________
    do row=2,n         
      factor=A(row,row-1)/A(row-1,row-1)
      A(row,row) = A(row,row) - factor*A(row-1,row)
      b(row) = b(row) - factor*b(row-1)
    end do
    !______ Back supstitution  _______________________________________
    x(n) = b(n)/A(n,n) ! Back substitution 
    do row=n-1,1,-1
      x(row) = (b(row) - A(row,row+1)*x(row+1))/A(row,row)
    end do
    allocate(thomasTimes)
    thomasTimes%node = x
  end function
