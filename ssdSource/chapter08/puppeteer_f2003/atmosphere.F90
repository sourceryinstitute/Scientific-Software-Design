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

module atmosphere_module
  use air_module ,only : air,air_ ! puppet for 1st Lorenz eq. and
                                  ! corresonding state variable
  use cloud_module ,only : cloud,cloud_ ! puppet for 2nd Lorenz eq. and 
                                        ! corresonding state variable
  use ground_module ,only : ground,ground_ ! puppet for 3rd Lorenz eq. and
                                           ! corresonding state variable
  use global_parameters_module ,only:debugging !print call tree if true

 !implicit none        ! Prevent implicit typing
  private              ! Hide everything by default
  public :: atmosphere,atmosphere_

  type atmosphere ! Puppeteer
    private 
    type(air)    ,allocatable :: air_puppet
    type(cloud)  ,allocatable :: cloud_puppet
    type(ground) ,allocatable :: ground_puppet
  contains
    procedure :: t ! time derivative
    procedure :: dRHS_dV ! Jacobian contribution (dR/dV)
    procedure :: state_vector ! return atmosphere solution vector
    procedure :: add ! add two atmospheres
    procedure :: subtract ! subtract one atmosphere from another
    procedure :: multiply ! multiply an atmosphere by a real scalar
    procedure ,pass(rhs) :: inverseTimes !abstract Gaussian elimination
    procedure :: assign ! assign one  atmosphere to another
    generic :: operator(+) => add
    generic :: operator(-) => subtract
    generic :: operator(*) => multiply
    generic :: assignment(=) => assign
    generic :: operator(.inverseTimes.) => inverseTimes
  end type atmosphere

  interface atmosphere_
    module procedure constructor
  end interface

contains

  ! constructor for an atmosphere object
  type(atmosphere) function constructor &
    (air_target,cloud_target,ground_target)
    type(air)    ,allocatable ,intent(inout) :: air_target
    type(cloud)  ,allocatable ,intent(inout) :: cloud_target
    type(ground) ,allocatable ,intent(inout) :: ground_target
    if (debugging) print *,'    atmosphere%construct(): start'
    ! transfer allocations from puppets to Puppeteer
    call move_alloc(air_target, constructor%air_puppet)
    call move_alloc(ground_target, constructor%ground_puppet)
    call move_alloc(cloud_target, constructor%cloud_puppet)
    if (debugging) print *,'    atmosphere%construct(): end'
  end function 

  ! time derivative (evolution equations)
  function t(this) result(dState_dt)
    class(atmosphere) ,intent(in)  :: this
    class(atmosphere) ,allocatable :: dState_dt 
    type(atmosphere)  ,allocatable :: delta
    if (debugging) print *,'    atmosphere%t(): start'
    allocate(delta) 
    delta%air_puppet= this%air_puppet%t(this%cloud_puppet%coordinate())
    delta%cloud_puppet = this%cloud_puppet%t( &
      this%air_puppet%coordinate(),this%ground_puppet%coordinate())
    delta%ground_puppet = this%ground_puppet%t( &
      this%air_puppet%coordinate(),this%cloud_puppet%coordinate())
    call move_alloc(delta, dState_dt)
    if (debugging) print *,'    atmosphere%t(): end'
  end function 

  ! atmosphere contribution to Jacobian
  function dRHS_dV(this) result(dRHS_dState)
    class(atmosphere)    ,intent(in)  :: this

    ! Sub-blocks of dR/dV array
    real ,dimension(:,:) ,allocatable ::                             &
                          dAir_dAir ,dAir_dCloud ,dAir_dGround,      &
                          dCloud_dAir ,dCloud_dCloud ,dCloud_dGround,&
                          dGround_dAir ,dGround_dCloud ,dGround_dGround

    real ,dimension(:,:) ,allocatable :: dRHS_dState

    real ,dimension(:)   ,allocatable ::                             &
                   air_coordinate, cloud_coordinate, ground_coordinate

    integer :: air_eqs ,air_vars,cloud_eqs ,cloud_vars , &
               ground_eqs ,ground_vars ,i ,j ,rows ,cols

    if (debugging) print *,'atmosphere%dRHS_dV(): start'

       ! Calculate matrices holding partial derivative of puppet
       ! evolution equation right-hand sides with respect to the
       ! dependent variables of each puppet.

         air_coordinate  = this%air_puppet%coordinate()
         cloud_coordinate= this%cloud_puppet%coordinate()
         ground_coordinate=this%ground_puppet%coordinate()

         ! compute diagonal block submatrix
         dAir_dAir       = this%air_puppet%d_dAir(cloud_coordinate)
         dCloud_dCloud   = this%cloud_puppet%d_dCloud(air_coordinate,&
                                                      ground_coordinate)
         dGround_dGround = this%ground_puppet%d_dGround(air_coordinate,&
                                                       cloud_coordinate)
         air_eqs     = size(dAir_dAir,1)        ! submatrix rows
         air_vars    = size(dAir_dAir,2)        ! submatrix columns
         cloud_eqs   = size(dCloud_dCloud,1)    ! submatrix rows
         cloud_vars  = size(dCloud_dCloud,2)    ! submatrix columns
         ground_eqs  = size(dGround_dGround,1)  ! submatrix rows
         ground_vars = size(dGround_dGround,2)  ! submatrix columns

         ! compute off-diagonal values
         dAir_dCloud    = this%air_puppet%d_dy(cloud_coordinate)
         dAir_dGround = reshape(source=(/(0.,i=1,air_eqs*ground_vars)/)&
                                        ,shape=(/air_eqs,ground_vars/))
         dCloud_dAir    = this%cloud_puppet%d_dx(air_coordinate,     &
                                                 ground_coordinate)
         dCloud_dGround = this%cloud_puppet%d_dz(air_coordinate,     &
                                                 ground_coordinate)
         dGround_dAir   = this%ground_puppet%d_dx(air_coordinate,    &
                                                  cloud_coordinate)
         dGround_dCloud = this%ground_puppet%d_dy(air_coordinate,    &
                                                  cloud_coordinate)
         rows=air_eqs+cloud_eqs+ground_eqs
         cols=air_vars+cloud_vars+ground_vars
         allocate(dRHS_dState(rows,cols))

         ! Begin result assembly
         dRHS_dState(1:air_eqs, 1:air_vars)                     =    &
                       dAir_dAir
         dRHS_dState(1:air_eqs, air_vars+1:air_vars+cloud_vars) =    &
                       dAir_dCloud
         dRHS_dState(1:air_eqs, air_vars+cloud_vars+1:cols)     =    &
                       dAir_dGround
         dRHS_dState(air_eqs+1:air_eqs+cloud_eqs, 1:air_vars)   =    &
                       dCloud_dAir
         dRHS_dState(air_eqs+1:air_eqs+cloud_eqs,                    &
                       air_vars+1:air_vars+cloud_vars)          =    &
                       dCloud_dCloud
         dRHS_dState(air_eqs+1:air_eqs+cloud_eqs,                    &
                       air_vars+cloud_vars+1:cols)              =    &
                       dCloud_dGround
         dRHS_dState(air_eqs+cloud_eqs+1:rows, 1:air_vars)      =    &
                       dGround_dAir
         dRHS_dState(air_eqs+cloud_eqs+1:rows,                       &
                       air_vars+1:air_vars+cloud_vars)          =    &
                       dGround_dCloud
         dRHS_dState(air_eqs+cloud_eqs+1:rows,                       &
                       air_vars+cloud_vars+1:cols)              =    &
                       dGround_dGround
      ! result assembly finishes
    if (debugging) print *,'    atmosphere%dRHS_dV(): end'
  end function dRHS_dV

  ! assemble and return solution vector
  function state_vector(this) result(phase_space)
    class(atmosphere)  ,intent(in)  :: this
    real ,dimension(:) ,allocatable :: state 
    real ,dimension(:) ,allocatable :: x,y,z,phase_space
    integer                         :: x_start,y_start,z_start
    integer                         :: x_end  ,y_end  ,z_end
    x = this%air_puppet%coordinate()
    x_start=1
    x_end=x_start+size(x)-1
    y = this%cloud_puppet%coordinate()
    y_start=x_end+1
    y_end=y_start+size(y)-1
    z = this%ground_puppet%coordinate()
    z_start=y_end+1
    z_end=z_start+size(z)-1
    allocate(phase_space(size(x)+size(y)+size(z)))
    phase_space(x_start:x_end) = x
    phase_space(y_start:y_end) = y
    phase_space(z_start:z_end) = z
  end function

  function add(lhs,rhs) result(sum)
    class(atmosphere) ,intent(in)  :: lhs
    class(atmosphere) ,intent(in)  :: rhs
    class(atmosphere) ,allocatable :: sum
    if (debugging) print *,'    atmosphere%add(): start'
    allocate(sum)
    sum%air_puppet    = lhs%air_puppet    + rhs%air_puppet
    sum%cloud_puppet  = lhs%cloud_puppet  + rhs%cloud_puppet
    sum%ground_puppet = lhs%ground_puppet + rhs%ground_puppet
    if (debugging) print *,'    atmosphere%add(): end'
  end function 

  function subtract(lhs,rhs) result(difference)
    class(atmosphere)       ,intent(in)  :: lhs
    class(atmosphere) ,intent(in)  :: rhs
    class(atmosphere) ,allocatable :: difference
    if (debugging) print *,' atmosphere%subtract(): start'
    allocate(difference) 
    difference%air_puppet = lhs%air_puppet - rhs%air_puppet
    difference%cloud_puppet = lhs%cloud_puppet - rhs%cloud_puppet
    difference%ground_puppet = lhs%ground_puppet - rhs%ground_puppet
    if (debugging) print *,'    atmosphere%subtract(): end'
  end function

  ! Solve linear system Ax=b by Gaussian elimination
  function inverseTimes(lhs,rhs) result(product)
    class(atmosphere) ,intent(in) :: rhs
    class(atmosphere) ,allocatable :: product
    real ,dimension(:,:) ,allocatable ,intent(in) :: lhs 
    real ,dimension(:)   ,allocatable :: x,b
    real ,dimension(:,:) ,allocatable :: A
    real ,parameter  :: pivot_tolerance=1.0E-02
    integer :: row,col,n,p ! p=pivot row/col 
    real :: factor

    n=size(lhs,1)
    b = rhs%state_vector()
    if (n /= size(lhs,2) .or. n /= size(b)) &
      stop 'atmosphere.f03: ill-posed matrix problem in inverseTimes()'
    allocate(x(n))
    A = lhs 
    do p=1,n-1         ! Forward elimination
      if (abs(A(p,p))<pivot_tolerance) &
        stop 'invert: use an algorithm with pivoting'
      do row=p+1,n
        factor=A(row,p)/A(p,p)
        forall(col=p:n)
          A(row,col) = A(row,col) - A(p,col)*factor
        end forall
        b(row) = b(row) - b(p)*factor
      end do
    end do
    x(n) = b(n)/A(n,n) ! Back substitution
    do row=n-1,1,-1
      x(row) = (b(row) - sum(A(row,row+1:n)*x(row+1:n)))/A(row,row)
    end do
    allocate(product,source=rhs) 
    product%air_puppet = air_(x(1),x(2))
    call product%cloud_puppet%set_coordinate(x(3))
    call product%ground_puppet%set_coordinate(x(4))
  end function 

  ! multiply atmosphere object by a real scalar
  function multiply(lhs,rhs) result(product)
    class(atmosphere) ,intent(in)  :: lhs
    real              ,intent(in)  :: rhs
    class(atmosphere) ,allocatable :: product
    if (debugging) print *,'    atmosphere%multiply(): start'
    allocate(product) 
    product%air_puppet    = lhs%air_puppet    * rhs
    product%cloud_puppet  = lhs%cloud_puppet  * rhs
    product%ground_puppet = lhs%ground_puppet * rhs
    if (debugging) print *,'    atmosphere%multiply(): end'
  end function 

  ! assign one atmosphere object to another
  subroutine assign(lhs,rhs)
    class(atmosphere) ,intent(inout) :: lhs
    class(atmosphere) ,intent(in)    :: rhs
    if (debugging) print *,'    atmosphere%assign(): start'
    lhs%air_puppet    = rhs%air_puppet    ! automatic allocation
    lhs%cloud_puppet  = rhs%cloud_puppet  ! automatic allocation
    lhs%ground_puppet = rhs%ground_puppet ! automatic allocation
    if (debugging) print *,'    atmosphere%assign(): end'
  end subroutine 
end module atmosphere_module
