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
#include "ForTrilinos_config.h"
#ifdef HAVE_MPI
  use mpi
  use FEpetra_MpiComm, only:Epetra_MpiComm
#else
  use FEpetra_SerialComm, only:Epetra_SerialComm
#endif
  use ForTrilinos_utils, only : valid_kind_parameters
  use iso_c_binding, only : c_int,c_double
  use field_module ,only : field,initial_field
  use field_factory_module ,only : field_factory
  use periodic_2nd_order_factory_module,only: periodic_2nd_order_factory
  use initializer ,only : u_initial,zero
  implicit none
#ifdef HAVE_MPI
  type(Epetra_MpiComm) :: comm
#else
  type(Epetra_SerialComm) :: comm
#endif
  class(field), pointer :: u,half_uu,u_half
  class(field_factory), allocatable :: field_creator
  procedure(initial_field) ,pointer :: initial 
  real(c_double) :: dt,half=0.5,t=0.,t_final=0.6,nu=1.
  real(c_double) :: t_start,t_end
  integer(c_int) :: tstep
  integer(c_int), parameter :: grid_resolution=(64**3)/4
  integer(c_int)      :: MyPID, NumProc
  logical             :: verbose
  integer :: rc,ierr 
 
  if (.not. valid_kind_parameters()) &
    stop 'C interoperability not supported on this platform.'
#ifdef HAVE_MPI
  call MPI_INIT(ierr) 
  t_start=MPI_Wtime()
  comm = Epetra_MpiComm(MPI_COMM_WORLD)
#else
  call cpu_time(t_start) 
  comm = Epetra_SerialComm()
#endif
  allocate(periodic_2nd_order_factory :: field_creator)
  initial => u_initial
  u => field_creator%create(initial,grid_resolution,comm)
  initial => zero
  half_uu => field_creator%create(initial,grid_resolution,comm)
  u_half => field_creator%create(initial,grid_resolution,comm)
  do tstep=1,20 !2nd-order Runge-Kutta:
   dt = u%runge_kutta_2nd_step(nu ,grid_resolution)
    half_uu = u*u*half
    u_half = u + (u%xx()*nu - half_uu%x())*dt*half ! first substep
    half_uu = u_half*u_half*half
    u  = u + (u_half%xx()*nu - half_uu%x())*dt ! second substep
    t = t + dt
    if (comm%MyPID()==0) print *,'timestep=',tstep
  end do
  if (comm%MyPID()==0) write(10,*) 'u at t=',t
#ifdef HAVE_MPI
  t_end= MPI_Wtime()
#else
  call cpu_time(t_end)
#endif
  if (comm%MyPID()==0) write(10,*) 'Elapsed CPU time=',t_end-t_start
  call u%output(comm)
  call half_uu%force_finalize
  call u_half%force_finalize
  call u%force_finalize
  call comm%force_finalize
#ifdef HAVE_MPI
  call MPI_FINALIZE(rc)
#endif
end program
