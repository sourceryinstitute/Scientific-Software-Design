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

module field_class
  use kind_parameters      ,only : ikind, rkind
  use differentiator_class ,only : differentiator,differentiator_
  implicit none 
  private
  abstract interface
    pure function distribution (x) result(ret)
        import
        real(rkind), intent(in) :: x
        real(rkind) ret
    end function
  end interface
  type field
    private
    real(rkind) ,allocatable ,dimension(:) :: node ! internal points
  contains
    procedure           :: nodal_values
    procedure           :: field_size
    procedure ,private  :: product
    generic             :: operator(*) => product
    procedure ,private  :: ratio
    generic             :: operator(/) => ratio
    procedure ,private  :: total
    generic             :: operator(+) => total
    procedure ,private ,pass(rhs) :: inverseTimes
    generic             :: operator(.inverseTimes.) => inverseTimes
    procedure           :: xx          ! 2nd-order spatial derivative
    procedure           :: xx_boundary ! 2nd-order boundary derivative
    procedure           :: xx_matrix   ! matrix derivative operator 
  end type
  interface field
    module procedure constructor
  end interface
  public :: field, distribution
  type(differentiator) ,protected    :: stencil
  integer(ikind)       ,parameter    :: end_points=2
  real(rkind) ,dimension(end_points) :: boundary
contains
  type(field) function constructor(spec,sample)
    use problem_class ,only : problem
    type(problem)  :: spec
    procedure (distribution) sample
    real(rkind)    :: dx
    integer(ikind) :: i
    stencil = differentiator_(spec)
    dx      = spec%spacing()
    allocate(constructor%node(spec%nodes()))
    forall(i=1:spec%nodes()) constructor%node(i) = sample(i*dx)
    boundary = spec%boundary_vals()
  end function

  pure type(field) function xx(this) 
    class(field) ,intent(in) :: this
    xx%node = stencil%laplacian(this%node,(/boundary(1),boundary(2)/))
  end function

  pure real(rkind) function xx_boundary(this) 
    class(field) ,intent(in) :: this
    type(field)              :: this_xx
    this_xx = this%xx()
    xx_boundary = this_xx%node(1)
  end function

  pure function xx_matrix(this)
    class(field) ,intent(in) :: this
    real(rkind) ,dimension(:,:) ,allocatable :: xx_matrix
    xx_matrix = stencil%lap_matrix()
  end function

  pure integer(ikind) function field_size(this) 
    class(field) ,intent(in) :: this
    field_size = size(this%node)
  end function

  pure function nodal_values(this) 
    class(field) ,intent(in) :: this
    real(rkind) ,allocatable ,dimension(:) :: nodal_values
    nodal_values = this%node
  end function

  pure type(field) function ratio(this,denominator) 
    class(field) ,intent(in) :: this
    real(rkind)  ,intent(in) :: denominator
    ratio%node = this%node/denominator
  end function

  pure type(field) function product(this,factor) 
    class(field) ,intent(in) :: this
    real(rkind)  ,intent(in) :: factor
    product%node = this%node*factor
  end function

  pure type(field) function total(lhs,rhs) 
    class(field) ,intent(in) :: lhs
    type(field)  ,intent(in) :: rhs
    total%node = lhs%node + rhs%node
  end function

  function inverseTimes(lhs,rhs)
    use linear_solve_module ,only : gaussian_elimination ! See Appendix A
    real(rkind) ,dimension(:,:) ,intent(in) :: lhs
    class(field)                             ,intent(in) :: rhs
    type(field)                 ,allocatable :: inverseTimes
    allocate(inverseTimes)
    inverseTimes%node = gaussian_elimination(lhs,rhs%node)
  end function
end module field_class
