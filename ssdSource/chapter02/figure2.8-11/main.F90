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

program fin_test 
  use iso_fortran_env ,only : input_unit
  use kind_parameters ,only : rkind
  use problem_class   ,only : problem,problem_
  use differentiator_class ,only : differentiator,differentiator_
  use conductor_class      ,only : conductor,conductor_
  implicit none

  real(rkind), parameter :: tolerance=1.0E-06
  type(problem)   :: specification
  type(conductor) :: fin
  type(differentiator) :: finite_difference

  specification      = problem_(input_unit)
  finite_difference  = differentiator_(specification)
  fin = conductor_(finite_difference,specification)
  print '(a,5g9.3)','initial temperature = ',fin%temperature()
  call fin%heat_for(specification%time_step())
  print '(a,5g9.3)','final temperature   = ',fin%temperature()

  if (abs(fin%time_derivative())<tolerance) then
    print '(2(a,es9.3))','|dT/dt|=',fin%time_derivative(),'<',tolerance
    print *,'In main: test passed. :)'
  else
    print '(2(a,es9.3))','|dT/dt|=',fin%time_derivative(),'>',tolerance
    print *,'In main: test failed. :('
  end if
end program
