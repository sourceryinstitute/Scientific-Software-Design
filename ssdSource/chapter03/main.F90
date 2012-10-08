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

program integrable_fin_test 
  use iso_fortran_env            ,only : input_unit
  use kind_parameters            ,only : rkind,ikind
  use initializer                ,only : linear
  use problem_class              ,only : problem,problem_
  use integrable_conductor_class ,only : integrable_conductor,integrable_conductor_
  implicit none
  real(rkind) ,dimension(:,:) ,allocatable :: I ! identity matrix
  real(rkind) ,parameter    :: tolerance=1.0E-06
  real(rkind)               :: dt=0.1           ! default time step 
  type(integrable_conductor):: fin              ! heat equation solver
  integer(ikind)            :: scheme           ! quadrature choice
  type(problem)             :: specs            ! problem specifications
  namelist /test_suite/ scheme

  enum ,bind(c) 
    enumerator forward_euler,backward_euler
  end enum
  if (.not. get_scheme_from(input_unit)) scheme=forward_euler ! default
  specs = problem_(input_unit)
  fin   = integrable_conductor_(specs,linear)
  print '(a,5g9.2)','initial temperature = ',fin%temperature()

  dt = specs%time_step()
  select case (scheme)
    case(forward_euler)
      fin = fin + fin%t()*dt
    case(backward_euler)
      I = identity(fin%rhs_operator_size())
      fin = ( I - fin%rhs_operator()*dt ) .inverseTimes. fin
    case default; stop 'In main: no method specified.'
  end select

  print '(a,5g9.2)','final temperature   = ',fin%temperature()
  if (abs(fin%time_derivative())<tolerance) then
    print '(2(a,es9.3))','|dT/dt|=',fin%time_derivative(),'<',tolerance
    print *,'In main: test passed. :)'
  else
    print '(2(a,es9.3))','|dT/dt|=',fin%time_derivative(),'>',tolerance
    print *,'In main: test failed. :('
  end if
contains
  logical function get_scheme_from(file)
    integer(ikind) ,intent(in) :: file
    integer(ikind) ,parameter  :: found=0
    integer(ikind)             :: namelist_status
    read(file,nml=test_suite,iostat=namelist_status)
    if (namelist_status == found) then
      get_scheme_from = .true.
    else
      get_scheme_from = .false.
    end if
    rewind(file) ! future reads start from the file head
  end function
   
   pure function identity(n)
      integer(ikind) ,intent(in) :: n
      integer(ikind)             :: i
      real(rkind) ,dimension(:,:) ,allocatable :: identity
      allocate(identity(n,n))
      identity = 0.
      forall(i=1:n) identity(i,i) = 1.
  end function
end program
