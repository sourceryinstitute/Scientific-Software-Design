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

module periodic_6th_order_module
  use field_module ,only : field,initial_field
  use kind_parameters ,only : rkind,ikind
  use matrix_module       ,only : matrix, new_matrix 
  implicit none
  private
  public :: periodic_6th_order,periodic_6th_order_
  type ,extends(field) :: periodic_6th_order
    private
    real(rkind) ,dimension(:) ,allocatable :: f
  contains
    procedure :: add => total
    procedure :: assign => copy
    procedure :: subtract => difference
    procedure :: multiply_field => product
    procedure :: multiply_real => multiple
    procedure :: runge_kutta_2nd_step => rk2_dt
    procedure :: x  => df_dx   ! 1st derivative w.r.t. x
    procedure :: xx => d2f_dx2 ! 2nd derivative w.r.t. x
    procedure :: output
    procedure :: assert
  end type

  real(rkind) ,dimension(:) ,allocatable :: x_node

  interface periodic_6th_order_
    module procedure constructor
  end interface 

  real(rkind) ,parameter :: pi=acos(-1._rkind)
  real(rkind) ,parameter :: first_coeff_6th(5)= (/1.0_rkind/3.0_rkind,&
    0.0_rkind,14.0_rkind/9.0_rkind,1.0_rkind/9.0_rkind,0.0_rkind/)
  real(rkind) ,parameter :: second_coeff_6th(5)=(/2.0_rkind/11.0_rkind,&
    0.0_rkind,12.0_rkind/11.0_rkind,3.0_rkind/11.0_rkind,0.0_rkind/) 

contains
  function constructor(initial,grid_resolution)
    type(periodic_6th_order) ,pointer :: constructor
    procedure(initial_field) ,pointer :: initial
    integer(ikind) ,intent(in) :: grid_resolution
    integer :: i
    allocate(constructor)
    allocate(constructor%f(grid_resolution))
    if (.not. allocated(x_node)) x_node = grid()
    forall (i=1:size(x_node)) constructor%f(i)=initial(x_node(i))  
  contains
    pure function grid()
      integer(ikind) :: i
      real(rkind) ,dimension(:) ,allocatable :: grid
      allocate(grid(grid_resolution))
      forall(i=1:grid_resolution) &
        grid(i)  = 2.*pi*real(i-1,rkind)/real(grid_resolution,rkind)  
    end function
  end function

  real(rkind) function rk2_dt(this,nu, grid_resolution) 
    class(periodic_6th_order) ,intent(in) :: this
    real(rkind) ,intent(in) :: nu
    integer(ikind) ,intent(in) :: grid_resolution
    real(rkind)             :: dx, CFL, k_max 
    dx=2.0*pi/grid_resolution
    k_max=grid_resolution/2.0_rkind
    CFL=2.0/(24.0*(1-cos(k_max*dx))/11.0/(1.0+4.0/11.0*cos(k_max*dx))+ &
        3.0*(1.0-cos(2.0*k_max*dx))/22.0/(1.0+4.0/11.0*cos(k_max*dx)))
    rk2_dt = CFL*dx**2/nu 
  end function

  function total(lhs,rhs)
    class(periodic_6th_order) ,intent(in) :: lhs
    class(field) ,intent(in) :: rhs
    class(field) ,allocatable :: total
    type(periodic_6th_order) ,allocatable :: local_total
    select type(rhs)
      class is (periodic_6th_order)
        allocate(local_total)
        local_total%f = lhs%f + rhs%f
        call move_alloc(local_total,total)
      class default
        stop 'periodic_6th_order%total: unsupported rhs class.'
    end select
  end function

  function difference(lhs,rhs)
    class(periodic_6th_order) ,intent(in) :: lhs
    class(field) ,intent(in)  :: rhs
    class(field) ,allocatable :: difference
    type(periodic_6th_order) ,allocatable :: local_difference
    select type(rhs)
      class is (periodic_6th_order)
        allocate(local_difference)
        local_difference%f = lhs%f - rhs%f
        call move_alloc(local_difference,difference)
      class default
        stop 'periodic_6th_order%difference: unsupported rhs class.'
    end select
  end function

  function product(lhs,rhs)
    class(periodic_6th_order) ,intent(in) :: lhs
    class(field) ,intent(in)  :: rhs
    class(field) ,allocatable :: product
    type(periodic_6th_order) ,allocatable :: local_product
    select type(rhs)
      class is (periodic_6th_order)
        allocate(local_product)
        local_product%f = lhs%f * rhs%f
        call move_alloc(local_product,product)
      class default
        stop 'periodic_6th_order%product: unsupported rhs class.'
    end select
  end function

  function multiple(lhs,rhs)
    class(periodic_6th_order) ,intent(in) :: lhs
    real(rkind) ,intent(in)  :: rhs
    class(field) ,allocatable :: multiple
    type(periodic_6th_order) ,allocatable :: local_multiple
    allocate(local_multiple)
    local_multiple%f = lhs%f * rhs
    call move_alloc(local_multiple,multiple)
  end function

  subroutine copy(lhs,rhs)
    class(field) ,intent(in) :: rhs
    class(periodic_6th_order) ,intent(inout) :: lhs
    select type(rhs)
      class is (periodic_6th_order)
        lhs%f = rhs%f
      class default
        stop 'periodic_6th_order%copy: unsupported copy class.'
    end select
  end subroutine

  function df_dx(this) 
    class(periodic_6th_order) ,intent(in) :: this
    class(field) ,allocatable  :: df_dx
    class(matrix),allocatable,save  :: lu_matrix 
    integer(ikind) :: i,nx, x_east, x_west
    integer(ikind) :: x_east_plus1,x_east_plus2
    integer(ikind) :: x_west_minus1,x_west_minus2
    real(rkind) ,dimension(:,:) ,allocatable :: A
    real(rkind) ,dimension(:)   ,allocatable :: b,coeff
    real(rkind) :: dx
    class(periodic_6th_order) ,allocatable :: df_dx_local

    nx=size(x_node)
    dx=2.*pi/real(nx,rkind)
    coeff = first_coeff_6th 
    if (.NOT. allocated(lu_matrix)) allocate(lu_matrix) 
    if (.NOT. lu_matrix%is_built()) then
      allocate(A(nx,nx))
      !__________Initialize coeffecient matrix A _____
      A=0.0_rkind
      do i=1, nx
        x_east = mod(i,nx)+1
        x_west = nx-mod(nx+1-i,nx)
        if (i==2) then
          x_east_plus1=x_east+1; x_west_minus1=nx
          x_east_plus2=x_east+2; x_west_minus2=nx-1
        else if (i==3) then
          x_east_plus1=x_east+1; x_west_minus1=1
          x_east_plus2=x_east+2; x_west_minus2=nx
        else if (i==nx-1) then
          x_east_plus1=1; x_west_minus1=x_west-1 
          x_east_plus2=2; x_west_minus2=x_west-2
        else if (i==nx-2) then
          x_east_plus1=nx; x_west_minus1=x_west-1
          x_east_plus2=1; x_west_minus2=x_west-2
        else
          x_east_plus1=x_east+1; x_west_minus1=x_west-1
          x_east_plus2=x_east+2; x_west_minus2=x_west-2
        end if
        A(i,x_west_minus1) =coeff(2)
        A(i,x_west)   =coeff(1)
        A(i,i)        =1.0_rkind
        A(i,x_east)   =coeff(1)
        A(i,x_east_plus1) =coeff(2)
      end do
      lu_matrix=new_matrix(A)
      deallocate(A)
    end if 
    allocate(b(nx))
    b=0.0
    do i=1,nx
      x_east = mod(i,nx)+1
      x_west = nx-mod(nx+1-i,nx)
      if (i==2) then
        x_east_plus1=x_east+1; x_west_minus1=nx
        x_east_plus2=x_east+2; x_west_minus2=nx-1
      else if (i==3) then
        x_east_plus1=x_east+1; x_west_minus1=1
        x_east_plus2=x_east+2; x_west_minus2=nx
      else if (i==nx-1) then
        x_east_plus1=1; x_west_minus1=x_west-1
        x_east_plus2=2; x_west_minus2=x_west-2
      else if (i==nx-2) then
        x_east_plus1=nx; x_west_minus1=x_west-1
        x_east_plus2=1; x_west_minus2=x_west-2
      else
        x_east_plus1=x_east+1; x_west_minus1=x_west-1
        x_east_plus2=x_east+2; x_west_minus2=x_west-2
      end if

      b(i)=(0.25*coeff(4)*(this%f(x_east_plus1)-this%f(x_west_minus1))+&
          0.5*coeff(3)*(this%f(x_east)-this%f(x_west))+ &
          coeff(5)/6.0*(this%f(x_east_plus2)-this%f(x_west_minus2)))/dx
    end do
    allocate(df_dx_local)
    df_dx_local%f=lu_matrix .inverseTimes. b
    call move_alloc(df_dx_local, df_dx)
  end function

  function d2f_dx2(this)
    class(periodic_6th_order)  ,intent(in)   :: this
    class(field) ,allocatable :: d2f_dx2
    class(matrix),allocatable, save      :: lu_matrix 
    integer(ikind) :: i,nx,x_east,x_west
    integer(ikind) :: x_east_plus1,x_east_plus2
    integer(ikind) :: x_west_minus1,x_west_minus2 
    real(rkind) ,dimension(:,:) ,allocatable :: A
    real(rkind) ,dimension(:)   ,allocatable :: coeff,b
    real(rkind)                              :: dx
    class(periodic_6th_order)  ,allocatable  :: d2f_dx2_local 

    nx=size(this%f)
    dx=2.*pi/real(nx,rkind)
    coeff = second_coeff_6th
    if (.NOT. allocated(lu_matrix)) allocate(lu_matrix) 
    if (.NOT. lu_matrix%is_built()) then
      allocate(A(nx,nx))

      !__________Initialize coeffecient matrix A _____
      A=0.0_rkind
      do i=1, nx
        x_east = mod(i,nx)+1
        x_west = nx-mod(nx+1-i,nx)
        if (i==2) then
          x_east_plus1=x_east+1; x_west_minus1=nx
          x_east_plus2=x_east+2; x_west_minus2=nx-1
        else if (i==3) then
          x_east_plus1=x_east+1; x_west_minus1=1
          x_east_plus2=x_east+2; x_west_minus2=nx
        else if (i==nx-1) then
          x_east_plus1=1; x_west_minus1=x_west-1
          x_east_plus2=2; x_west_minus2=x_west-2
        else if (i==nx-2) then
          x_east_plus1=nx; x_west_minus1=x_west-1
          x_east_plus2=1; x_west_minus2=x_west-2
        else
          x_east_plus1=x_east+1; x_west_minus1=x_west-1
          x_east_plus2=x_east+2; x_west_minus2=x_west-2
        end if
        A(i,x_west_minus1) =coeff(2)
        A(i,x_west)        =coeff(1)
        A(i,i)             =1.0_rkind
        A(i,x_east)        =coeff(1)
        A(i,x_east_plus1)  =coeff(2)
      end do
      lu_matrix=new_matrix(A)
      deallocate(A)
    end if
    allocate(b(nx))
    do i=1, nx
      x_east = mod(i,nx)+1
      x_west = nx-mod(nx+1-i,nx)
      if (i==2) then
        x_east_plus1=x_east+1; x_west_minus1=nx
        x_east_plus2=x_east+2; x_west_minus2=nx-1
      else if (i==3) then
        x_east_plus1=x_east+1; x_west_minus1=1
        x_east_plus2=x_east+2; x_west_minus2=nx
      else if (i==nx-1) then
        x_east_plus1=1; x_west_minus1=x_west-1
        x_east_plus2=2; x_west_minus2=x_west-2
      else if (i==nx-2) then
        x_east_plus1=nx; x_west_minus1=x_west-1
        x_east_plus2=1; x_west_minus2=x_west-2
      else
        x_east_plus1=x_east+1; x_west_minus1=x_west-1
        x_east_plus2=x_east+2; x_west_minus2=x_west-2
      end if
      b(i)=(0.25*coeff(4)* &
        (this%f(x_east_plus1)-2.0*this%f(i)+this%f(x_west_minus1))+ &
        coeff(3)*(this%f(x_east)-2.0*this%f(i)+this%f(x_west))+ &
        coeff(5)/9.0* &
        (this%f(x_east_plus2)-this%f(i)+this%f(x_west_minus2)))/dx**2
    end do
    allocate(d2f_dx2_local)
    d2f_dx2_local%f=lu_matrix .inverseTimes. b
    call move_alloc(d2f_dx2_local, d2f_dx2)
  end function
  
  subroutine output(this)
    class(periodic_6th_order) ,intent(in) :: this
    integer(ikind) :: i
    do i=1,size(x_node) 
      print *, x_node(i), this%f(i)
    end do
  end subroutine
 
  subroutine assert(this, node_index, check_value)
    class(periodic_6th_order) ,intent(in) :: this 
    integer                   ,intent(in) :: node_index
    real(rkind)               ,intent(in) :: check_value
    real(rkind)               ,parameter  :: epsilon=1.0e-6    
  
    if (abs(this%f(node_index)-check_value) <=epsilon) print *, 'Test passed'
  end subroutine

end module
