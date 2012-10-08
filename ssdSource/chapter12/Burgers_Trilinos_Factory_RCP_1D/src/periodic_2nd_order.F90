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
  use FEpetra_Comm, only:Epetra_Comm
  use FEpetra_Map!, only: Epetra_Map
  use FEpetra_Vector!, only:Epetra_Vector
  use FEpetra_CrsMatrix, only:Epetra_CrsMatrix
  use ForTrilinos_enum_wrappers
  use ForTrilinos_error 
  use field_module ,only : field,initial_field
  use iso_c_binding ,only : c_double,c_int
  use ForTrilinos_assertion_utility ,only: &
    error_message,assert,assert_identical
  implicit none
  private
  public :: periodic_2nd_order

  type ,extends(field) :: periodic_2nd_order
    private
    type(Epetra_Vector) :: f
  contains
    procedure :: add => total
    procedure :: subtract => difference
    procedure :: assign => copy
    procedure :: multiply_field => product
    procedure :: multiply_real => multiple
    procedure :: runge_kutta_2nd_step => rk2_dt
    procedure :: x  => df_dx     ! 1st derivative w.r.t. x
    procedure :: xx => d2f_dx2   ! 2nd derivative w.r.t. x
    procedure :: output
    procedure :: force_finalize
  end type
   
  real(c_double) ,parameter                 :: pi=acos(-1._c_double)
  real(c_double) ,dimension(:) ,allocatable :: x_node
  type(Epetra_Map)             ,allocatable :: map

  interface periodic_2nd_order
    procedure constructor
  end interface 
contains
  function constructor(initial,grid_resolution,comm) result(this)
    type(periodic_2nd_order) ,pointer  :: this
    procedure(initial_field) ,pointer :: initial
    integer(c_int) ,intent(in) :: grid_resolution
    integer(c_int) :: i,j
    class(Epetra_Comm), intent(in) :: comm
    integer(c_int) :: NumGlobalElements
    integer(c_int),dimension(:),allocatable :: MyGlobalElements
    integer(c_int)      :: NumMyElements,IndexBases=1,status
    real(c_double) ,dimension(:) ,allocatable :: f_v
    type(error) :: ierr
    NumGlobalElements=grid_resolution
    allocate(this)
    if (.not. allocated(x_node)) x_node = grid()
    if (.not. allocated(map)) then
      allocate(map,stat=status)
      ierr=error(status,'periodic_2nd_order: create map')
      call ierr%check_allocation()
      map = Epetra_Map(NumGlobalElements,IndexBases,comm)
    end if
    NumMyElements= map%NumMyElements()
    allocate(MyGlobalElements(NumMyElements))
    MyGlobalElements = map%MyGlobalElements()
    allocate(f_v(NumMyElements))
    forall(i=1:NumMyElements)f_v(i)=initial(x_node(MyGlobalElements(i)))
    this%f=Epetra_Vector(map,zero_initial=.true.)
    call this%f%ReplaceGlobalValues(NumMyElements,f_v,MyGlobalElements)
  contains
    pure function grid()
      integer(c_int) :: i
      real(c_double) ,dimension(:) ,allocatable :: grid
      allocate(grid(grid_resolution))
      forall(i=1:grid_resolution) &
        grid(i)= 2.*pi*real(i-1,c_double)/real(grid_resolution,c_double)
    end function
  end function

  subroutine copy(lhs,rhs)
    class(field) ,intent(in) :: rhs
    class(periodic_2nd_order) ,intent(inout) :: lhs
    select type(rhs)
      class is (periodic_2nd_order)
        lhs%f = rhs%f
      class default
        stop 'periodic_2nd_order%copy: unsupported copy class.'
    end select
  end subroutine

  real(c_double) function rk2_dt(this,nu, grid_resolution)
    class(periodic_2nd_order) ,intent(in) :: this
    real(c_double) ,intent(in) :: nu
    integer(c_int) ,intent(in) :: grid_resolution
    real(c_double)             :: dx, CFL, k_max
    dx=2.0*pi/grid_resolution
    k_max=grid_resolution/2.0_c_double
    CFL=1.0/(1.0-cos(k_max*dx))
    rk2_dt = CFL*dx**2/nu
  end function

  function total(lhs,rhs)
    class(periodic_2nd_order) ,intent(in) :: lhs
    class(field) ,intent(in) :: rhs
    class(field) ,allocatable :: total
    type(periodic_2nd_order) ,allocatable :: local_total
    select type(rhs)
      class is (periodic_2nd_order)
        allocate(periodic_2nd_order::local_total)
        local_total%f=Epetra_Vector(map,zero_initial=.true.)
        call local_total%f%Update( &
          1._c_double,lhs%f,1._c_double,rhs%f,0._c_double)
        call move_alloc(local_total,total)
      class default
        stop 'periodic_2nd_order%total: unsupported rhs class.'
    end select
  end function

  function difference(lhs,rhs)
    class(periodic_2nd_order) ,intent(in) :: lhs
    class(field) ,intent(in)  :: rhs
    class(field) ,allocatable :: difference
   type(periodic_2nd_order) ,allocatable :: local_difference
    select type(rhs)
      class is (periodic_2nd_order)
        allocate(periodic_2nd_order::local_difference)
        local_difference%f=Epetra_Vector(map,zero_initial=.true.)
        call local_difference%f%Update(&
          1._c_double,lhs%f,-1._c_double,rhs%f,0._c_double)
        call move_alloc(local_difference,difference)
      class default
        stop 'periodic_2nd_order%difference: unsupported rhs class.'
    end select
  end function

 function product(lhs,rhs)
   class(periodic_2nd_order) ,intent(in) :: lhs
   class(field) ,intent(in)  :: rhs
   class(field) ,allocatable :: product
   type(periodic_2nd_order) ,allocatable :: local_product
   select type(rhs)
    class is (periodic_2nd_order)
      allocate(periodic_2nd_order::local_product)
      local_product%f=Epetra_Vector(map,zero_initial=.true.)
      call local_product%f%Multiply(1._c_double,lhs%f,rhs%f,0._c_double)
      call move_alloc(local_product,product)
     class default
      stop 'periodic_2nd_order%product: unsupported rhs class.'
   end select
  end function

  function multiple(lhs,rhs)
    class(periodic_2nd_order) ,intent(in) :: lhs
    real(c_double) ,intent(in)  :: rhs
    class(field) ,allocatable :: multiple
    type(periodic_2nd_order) ,allocatable :: local_multiple
    allocate(periodic_2nd_order::local_multiple)
    local_multiple%f=Epetra_Vector(map,zero_initial=.true.)
    call local_multiple%f%Scale(rhs,lhs%f)
   call move_alloc(local_multiple,multiple)
  end function

  function df_dx(this) 
    class(periodic_2nd_order) ,intent(in) :: this
    class(field) ,allocatable  :: df_dx
    type(Epetra_Vector) :: x
    type(Epetra_CrsMatrix) :: A 
    type(error) :: err
    real(c_double) ,dimension(:)   ,allocatable :: c
    real(c_double) :: dx
    integer(c_int) :: nx
    type(periodic_2nd_order), allocatable :: df_dx_local
    integer(c_int),dimension(:),allocatable :: MyGlobalElements
    integer(c_int),dimension(:),allocatable :: MyGlobalElements_diagonal
    integer(c_int),dimension(:),allocatable :: NumNz
    integer(c_int) :: NumGlobalElements,NumMyElements,i
    integer(c_int) :: indices(2), NumEntries
    real(c_double) ::values(2)
    real(c_double),parameter :: zero =0.0
    integer(c_int),parameter :: diagonal=1

    ! Executable code
    nx=size(x_node)
    dx=2.*pi/real(nx,c_double)
    NumGlobalElements = nx

! Get update list and number of local equations from given Map
    NumMyElements = map%NumMyElements()
    call assert_identical( [NumGlobalElements,map%NumGlobalElements()] )
    allocate(MyGlobalElements(NumMyElements))
    MyGlobalElements = map%MyGlobalElements()

! Create an integer vector NumNz that is used to build the Epetra Matrix
! NumNz(i) is the number of non-zero elements for the ith global eqn.
! on this processor
  allocate(NumNz(NumMyElements))

! We are building a tridiagonal matrix where each row has (-1 0 1)
! So we need 2 off-diagonal terms (except for the first and last eqn.)
  NumNz = 3

! Create a Epetra_Matrix
  A = Epetra_CrsMatrix(FT_Epetra_DataAccess_E_Copy,map,NumNz)

! Add rows one at a time
! Need some vectors to help
! off diagonal values will always be -1 and 1
  values(1) = -1.0/(2.0*dx)
  values(2) = 1.0/(2.0*dx)
  do i=1,NumMyElements
    if (MyGlobalElements(i)==1) then
      indices(1) = NumGlobalElements 
      indices(2) = 2
      NumEntries = 2
    else if(MyGlobalElements(i)==NumGlobalElements) then
      indices(1) = NumGlobalElements-1
      indices(2) = 1
      NumEntries = 2
    else
      indices(1) = MyGlobalElements(i)-1
      indices(2) = MyGlobalElements(i)+1
      NumEntries = 2
    end if
     call A%InsertGlobalValues(&
       MyGlobalElements(i),NumEntries,values,indices,err)
     call assert( [err%error_code()==0_c_int] , &
       [error_message('A%InsertGlobalValues: failed')] )
  !Put in the diaogonal entry
     MyGlobalElements_diagonal=MyGlobalElements+i-1
     call A%InsertGlobalValues(MyGlobalElements(i), &
       diagonal,[zero],MyGlobalElements_diagonal,err)
     call assert( [err%error_code()==0_c_int] , &
       [error_message('A%InsertGlobalValues: failed')] )
  end do

  !Finish up
    call A%FillComplete()
  !create vector x 
    x=Epetra_Vector(A%RowMap())
    Call A%Multiply_Vector(.false.,this%f,x)
    allocate(c(NumMyElements))
    c=x%ExtractCopy()
 !create vector of df_dx
    allocate(periodic_2nd_order::df_dx_local)
    df_dx_local%f=Epetra_Vector(map,zero_initial=.true.)
    call df_dx_local%f%ReplaceGlobalValues(&
      NumMyElements,c,MyGlobalElements)
    call move_alloc(df_dx_local, df_dx)
  end function
  
  function d2f_dx2(this) 
    class(periodic_2nd_order) ,intent(in) :: this
    class(field) ,allocatable  :: d2f_dx2
    type(Epetra_Vector) :: x
    type(Epetra_CrsMatrix) :: A 
    type(error) :: err
    real(c_double) ,dimension(:)   ,allocatable :: c
    real(c_double) :: dx
    integer(c_int) :: nx
    type(periodic_2nd_order) ,allocatable :: d2f_dx2_local
    integer(c_int),dimension(:),allocatable :: MyGlobalElements,NumNz
    integer(c_int),dimension(:),allocatable :: MyGlobalElements_diagonal
    integer(c_int) :: NumGlobalElements,NumMyElements,i
    integer(c_int) :: indices(2),NumEntries
    real(c_double) :: values(2),two_dx2  
    integer(c_int),parameter :: diagonal=1

  ! Executable code
   nx=size(x_node)
   dx=2.*pi/real(nx,c_double)
   NumGlobalElements = nx

! Get update list and number of local equations from given Map
  NumMyElements = map%NumMyElements()
  call assert_identical( [NumGlobalElements,map%NumGlobalElements()] )
  allocate(MyGlobalElements(NumMyElements))
  MyGlobalElements = map%MyGlobalElements()

! Create an integer vector NumNz that is used to build the Epetra Matrix
! NumNz(i) is the number of non-zero elements for the ith global eqn.
! on this processor
  allocate(NumNz(NumMyElements))

! We are building a tridiagonal matrix where each row has (1 -2  1)
! So we need 2 off-diagonal terms (except for the first and last eqn.)
  NumNz = 3

! Create a Epetra_Matrix
  A = Epetra_CrsMatrix(FT_Epetra_DataAccess_E_Copy,map,NumNz)

! Add rows one at a time
! Need some vectors to help
! off diagonal values will always be 1 and 1
  values(1) = 1.0/(dx*dx)
  values(2) = 1.0/(dx*dx)
  two_dx2   =-2.0/(dx*dx)
  do i=1,NumMyElements
    if (MyGlobalElements(i)==1) then
      indices(1) = NumGlobalElements 
      indices(2) = 2
      NumEntries = 2
    else if(MyGlobalElements(i)==NumGlobalElements) then
      indices(1) = NumGlobalElements-1
      indices(2) = 1
      NumEntries = 2
    else
      indices(1) = MyGlobalElements(i)-1
      indices(2) = MyGlobalElements(i)+1
      NumEntries = 2
    end if
     call A%InsertGlobalValues( &
       MyGlobalElements(i),NumEntries,values,indices,err)
     call assert( [err%error_code()==0_c_int] , &
       [error_message('A%InsertGlobalValues: failed')] )
  !Put in the diaogonal entry
     MyGlobalElements_diagonal=MyGlobalElements+i-1
     call A%InsertGlobalValues( MyGlobalElements(i) &
       ,diagonal,[two_dx2],MyGlobalElements_diagonal,err)
     call assert( [err%error_code()==0_c_int] , &
       [error_message('A%InsertGlobalValues: failed')] )
  end do

  !Finish up
    call A%FillComplete()
  !create vector x 
    x=Epetra_Vector(A%RowMap())
    Call A%Multiply_Vector(.false.,this%f,x)
    allocate(c(NumMyElements))
    c=x%ExtractCopy()
  !create vector of df_dx
    allocate(periodic_2nd_order::d2f_dx2_local)
    d2f_dx2_local%f=Epetra_Vector(map,zero_initial=.true.)
    call d2f_dx2_local%f%ReplaceGlobalValues( &
      NumMyElements,c,MyGlobalElements)
    call move_alloc(d2f_dx2_local, d2f_dx2)
  end function

  subroutine output(this,comm)
    class(periodic_2nd_order) ,intent(in) :: this
    class(Epetra_Comm),intent(in) ::comm
    integer(c_int) :: i,NumMyElements,NumGlobalElements
    integer(c_int), dimension(:), allocatable :: MyGlobalElements
    real(c_double), dimension(:), allocatable :: f_v
    real(c_double), dimension(:), allocatable :: f
    NumGlobalElements=map%NumGlobalElements()
    NumMyElements=map%NumMyElements()
    allocate(MyGlobalElements(NumMyElements))
    MyGlobalElements=map%MyGlobalElements()
    allocate(f_v(NumMyElements))
    f_v=this%f%ExtractCopy()
    allocate(f(NumGlobalElements))
    call comm%GatherAll(f_v,f,NumMyElements)
    do i=1,NumGlobalElements
      if (comm%MyPID()==0) write(20,'(2(E20.12,1x))') x_node(i),f(i)
    enddo
  end subroutine
  
  subroutine force_finalize(this)
    class(periodic_2nd_order), intent(inout) :: this
    call this%f%force_finalize
  end subroutine
end module
