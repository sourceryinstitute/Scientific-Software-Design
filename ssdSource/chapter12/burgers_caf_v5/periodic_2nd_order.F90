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

module periodic_2nd_order_module
  use kind_parameters ,only : rkind, ikind
  use field_module ,only : field

  implicit none
  private
  public :: periodic_2nd_order, initial_field
  type periodic_2nd_order
    private
    real(rkind), allocatable :: global_f(:)[:]
  contains
    procedure :: construct
    procedure :: assign        => copy
    procedure :: add           => add_field
    procedure :: multiply      => multiply_field
    procedure :: x             => df_dx
    procedure :: xx            => d2f_dx2
    procedure :: output
    procedure :: runge_kutta_2nd_step => rk2_dt
    procedure :: assert
    generic   :: assignment(=) => assign
    generic   :: operator(+)   => add
    generic   :: operator(*)   => multiply
    ! The following procedures were not included in the textbook:
    procedure, nopass :: this_image_contains_midpoint
    procedure :: has_a_zero_at
  end type

  real(rkind) ,parameter :: pi=acos(-1._rkind)
  ! The following variable was not included in the textbook:
  real, allocatable :: local_grid(:)

  abstract interface
    real(rkind) pure function initial_field(x)
      import :: rkind
      real(rkind) ,intent(in) :: x
    end function
  end interface

  contains

  subroutine output (this)
    class(periodic_2nd_order), intent(in) :: this
    integer(ikind) :: i
    do i = 1, size(this%global_f)
        print *, (this_image()-1)*size(this%global_f) + i, &
                 this%global_f(i)
    end do
  end subroutine

  subroutine assert(this, node_index, check_value)
    class(periodic_2nd_order) ,intent(in) :: this
    integer                   ,intent(in) :: node_index
    real(rkind)               ,intent(in) :: check_value
    integer                               :: i
    real(rkind)               ,parameter  :: epsilon=1.0e-6
 
    do i=1, size(this%global_f)
      if (((this_image()-1)*size(this%global_f)+i)==node_index) then
        if (abs(this%global_f(i)-check_value) <=epsilon) print *, 'Test passed'
      end if
    end do
  end subroutine

  subroutine construct (this, initial, num_grid_pts)
    class(periodic_2nd_order), intent(inout) :: this
    procedure(initial_field) ,pointer :: initial
    integer(ikind) ,intent(in) :: num_grid_pts
    integer :: i, local_grid_size

    !<-- assume mod(num_grid_pts, num_images()) == 0
    local_grid_size = num_grid_pts / num_images()

    ! set up the global grid points
    allocate (this%global_f(local_grid_size)[*])

    local_grid = grid() ! This line was not in the textbook 
    this%global_f(:) = local_grid ! The text book version directly assigns grid() to this%globa_f(:)

    do i = 1, local_grid_size
        this%global_f(i) = initial(this%global_f(i))
    end do

    sync all

  contains
    pure function grid()
      integer(ikind) :: i
      real(rkind) ,dimension(:) ,allocatable :: grid
      allocate(grid(local_grid_size))
      do i=1,local_grid_size
        grid(i)  = 2.*pi*(local_grid_size*(this_image()-1)+i-1) &
                   /real(num_grid_pts,rkind)  
      end do
    end function
  end subroutine

  real(rkind) function rk2_dt(this,nu, num_grid_pts) 
    class(periodic_2nd_order) ,intent(in) :: this
    real(rkind) ,intent(in) :: nu
    integer(ikind) ,intent(in) :: num_grid_pts
    real(rkind)             :: dx, CFL, k_max 
    dx=2.0*pi/num_grid_pts
    k_max=num_grid_pts/2.0_rkind
    CFL=1.0/(1.0-cos(k_max*dx))
    rk2_dt = CFL*dx**2/nu 
  end function

  ! this is the assignment
  subroutine copy(lhs,rhs)
    class(periodic_2nd_order) ,intent(inout) :: lhs
    class(field) ,intent(in) :: rhs

    ! update global field
    lhs%global_f(:) = rhs%state()
    sync all
  end subroutine

  function add_field (this, rhs)
    class(periodic_2nd_order), intent(in) :: this
    class(field), intent(in) :: rhs
    class(field), allocatable :: add_field

    allocate (add_field)
    add_field = rhs%state()+this%global_f(:)
  end function

  function multiply_field (this, rhs)
    class(periodic_2nd_order), intent(in) :: this, rhs
    class(field), allocatable :: multiply_field

    allocate (multiply_field)
    multiply_field = this%global_f(:)*rhs%global_f(:)
  end function

  function df_dx(this)
    class(periodic_2nd_order), intent(in) :: this
    class(field) ,allocatable  :: df_dx
    integer(ikind) :: i,nx, me, east, west
    real(rkind) :: dx
    real(rkind), allocatable :: tmp_field_array(:)

    nx=size(this%global_f)
    dx=2.*pi/(real(nx,rkind)*num_images())

    allocate(df_dx,tmp_field_array(nx))

    me = this_image()

    if (me == 1) then 
        west = num_images()
        east = merge(1,2,num_images()==1)
    else if (me == num_images()) then
        west = me - 1
        east = 1
    else
        west = me - 1
        east = me + 1
    end if

    tmp_field_array(1) = &
       0.5*(this%global_f(2)-this%global_f(nx)[west])/dx

    tmp_field_array(nx) = &
       0.5*(this%global_f(1)[east]-this%global_f(nx-1))/dx

    do i=2,nx-1
      tmp_field_array(i)=&
        0.5*(this%global_f(i+1)-this%global_f(i-1))/dx
    end do

    df_dx = tmp_field_array
  end function

  function d2f_dx2(this)
    class(periodic_2nd_order), intent(in) :: this
    class(field) ,allocatable  :: d2f_dx2
    integer(ikind) :: i,nx, me, east, west
    real(rkind) :: dx
    real(rkind), allocatable :: tmp_field_array(:)

    nx=size(this%global_f)
    dx=2.*pi/(real(nx,rkind)*num_images())

    allocate(d2f_dx2,tmp_field_array(nx))

    me = this_image()

    if (me == 1) then 
        west = num_images()
        east = merge(1,2,num_images()==1)
    else if (me == num_images()) then
        west = me - 1
        east = 1
    else
        west = me - 1
        east = me + 1
    end if

    tmp_field_array(1) = &
       (this%global_f(2)-2.0*this%global_f(1)+this%global_f(nx)[west])&
       /dx**2

    tmp_field_array(nx) =&
       (this%global_f(1)[east]-2.0*this%global_f(nx)+this%global_f(nx-1))&
       /dx**2

    do i=2,nx-1
      tmp_field_array(i)=&
        (this%global_f(i+1)-2.0*this%global_f(i)+this%global_f(i-1))&
        /dx**2
    end do

    d2f_dx2 = tmp_field_array
  end function

  ! The following procedures were not included in the initial publication of the textbook:

  pure function has_a_zero_at(this, expected_location) result(zero_at_expected_location)
    class(periodic_2nd_order) ,intent(in) :: this
    real(rkind) ,intent(in) :: expected_location
    real(rkind), parameter :: tolerance = 1.0E-06_rkind
    integer :: nearest_grid_point
    logical :: zero_at_expected_location
    nearest_grid_point = minloc(abs(local_grid-expected_location),dim=1)
    zero_at_expected_location = merge(.true.,.false., abs(this%global_f(nearest_grid_point)) < tolerance  )
  end function

  function this_image_contains_midpoint() result(within_bounds)
    logical within_bounds
    !<-- assume mod(num_grid_pts, num_images()) == 0
    within_bounds = merge(.true.,.false., (this_image()==num_images()/2+1) )
  end function

end module
