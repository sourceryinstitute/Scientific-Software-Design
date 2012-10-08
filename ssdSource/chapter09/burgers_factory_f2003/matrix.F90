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

module matrix_module
  use kind_parameters ,only : rkind, ikind
  implicit none
  private
  public :: matrix
  public :: new_matrix 
  type  :: matrix 
    integer(ikind), allocatable :: pivot(:)
    real(rkind), allocatable :: lu(:,:) 
  contains
    procedure  :: back_substitute            
    procedure  :: is_built 
    procedure  :: matrix_eq_matrix 
    generic    :: operator(.inverseTimes.) => back_substitute 
    generic    :: assignment(=) => matrix_eq_matrix 
  end type
    
  interface new_matrix
    module procedure constructor
  end interface 

  contains
    logical function is_built(this)
      class(matrix),intent(in)   :: this 
      is_built = allocated (this%pivot)
    end function
    
    function constructor(A) result(new_matrix)
      class(matrix) ,allocatable               :: new_matrix
      real(rkind), intent(in) ,dimension(:,:)  :: A 
      integer                                  :: n,info
      
      n=size(A,1)
      allocate(new_matrix)
      allocate(new_matrix%pivot(n), new_matrix%lu(n,n))
      new_matrix%lu=A
      call dgetrf(n,n,new_matrix%lu,n,new_matrix%pivot,info)  
    end function
      
    subroutine matrix_eq_matrix(lhs,rhs)
      class(matrix) ,intent(in)  :: rhs
      class(matrix) ,intent(out) :: lhs
      lhs%pivot = rhs%pivot
      lhs%lu    = rhs%lu
    end subroutine

    function back_substitute(this,b) result(x)
      class(matrix) ,intent(in)                  :: this
      real(rkind)  ,intent(in)  ,dimension(:)   :: b
      real(rkind)  ,allocatable ,dimension(:)   :: x, temp_x
      real(rkind)  ,allocatable ,dimension(:,:) :: lower, upper 
      real(rkind)  ,allocatable ,dimension(:)   :: local_b 
      integer(ikind)                            :: n,i,j 
      real(rkind)                               :: temp 
       
      n=size(this%lu,1)
      allocate(lower(n,n), upper(n,n))
      lower=0.0_rkind
      upper=0.0_rkind
      do i=1,n
        do j=i,n
          upper(i,j)=this%lu(i,j)  
        end do
      end do
      do i=1,n
        lower(i,i)=1.0_rkind 
      end do
      do i=2,n
        do j=1,i-1 
          lower(i,j)=this%lu(i,j)
        end do
      end do
      allocate(local_b(n))
      local_b=b
      do i=1,n
        if (this%pivot(i)/=i) then 
          temp=local_b(i) 
          local_b(i)=local_b(this%pivot(i))
          local_b(this%pivot(i))=temp
        end if
      end do
      allocate(temp_x(n))
      temp_x(1)=local_b(1)
      do i=2,n
        temp_x(i)=local_b(i)-sum(lower(i,1:i-1)*temp_x(1:i-1)) 
      end do
      allocate(x(n))
      x(n)=temp_x(n)/upper(n,n)
      do i=n-1,1,-1
        x(i)=(temp_x(i)-sum(upper(i,i+1:n)*x(i+1:n)))/upper(i,i)
      end do
      do i=n,1,-1
        if (this%pivot(i)/=i) then
          temp=x(i)
          x(i)=x(this%pivot(i))
          x(this%pivot(i))=temp
        end if
      end do 
    end function
end module
