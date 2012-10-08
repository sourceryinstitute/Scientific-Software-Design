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

program main
  use iso_fortran_env ,only : input_unit
  use kind_parameters ,only : rkind,ikind
  use conduction_module
  implicit none
  real(rkind)    ,parameter :: tolerance=1.0E-06
  integer(ikind) ,parameter :: chip_end=1
  integer(ikind)            :: elements=4,intervals           
  real(rkind) :: air_temp=1.,chip_temp=1. & !boundary temperatures
    ,diffusivity=1.,fin_length=1.,time=0.  & !physical parameters
    ,dt=0.01,dx                              !time step, grid spacing
  real(rkind),dimension(:)  ,allocatable:: T,T_xx !temperature,2nd deriv
  real(rkind),dimension(:,:),allocatable:: laplacian !differential op
  dx=fin_length/elements                     !default element size
  if (.not. &
    read_physics(input_unit,diffusivity,fin_length,chip_temp,air_temp))&
    print *,'In main: Using default physics specification.'
  if (.not. read_numerics(input_unit,dt,elements,dx,fin_length)) &
    print *,'In main: Using default numerical specification.'
  allocate(T(elements),laplacian(elements,diagonals),T_xx(elements)) 
  T         = air_temp
  laplacian = differencing_stencil(dx,elements)
  print *,'time,temperature=',time,T
  dt = min(dt,stable_dt(dx,diffusivity))
  do
    time = time + dt 
    T_xx = differentiate(laplacian,T,chip_temp,air_temp) 
    T = T + dt*diffusivity*T_xx
    print *,'time,temperature=',time,T
    if (T_xx(chip_end)<tolerance) then
      print *,'Test passed.'
      exit
    end if
  end do
  print *,'steady state reached at time ',time
end program 
