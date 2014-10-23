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

module conduction_module
  use kind_parameters ,only : rkind,ikind
  implicit none            ! Tri-diagonal array positions:
  integer(ikind) ,parameter::diagonals=3,low_diag=1,diag=2,up_diag=3 
contains
  pure logical function found(namelist_io)
    integer(ikind) ,intent(in) ::namelist_io
    if (namelist_io==0)  then 
      found=.true.
    else
      found=.false.
    end if
  end function

  logical function read_physics(unit,alpha,L_fin,T_chip,T_air)
    real(rkind)    ,intent(out) :: alpha,L_fin,T_chip,T_air 
    integer(ikind) ,intent(in)  :: unit
    integer(ikind)              :: physics_io
    namelist/physics/ alpha,L_fin,T_chip,T_air
    read(unit,nml=physics,iostat=physics_io)
    read_physics = found(physics_io)
  end function 
 
  logical function read_numerics(unit,dt,nodes,dx,L_fin)
    real(rkind)    ,intent(out) :: dt,dx
    integer(ikind) ,intent(out) :: nodes
    real(rkind)    ,intent(in)  :: L_fin
    integer(ikind) ,intent(in)  :: unit
    integer(ikind)              :: numerics_io,elements
    namelist/numerics/ dt,nodes
    read(unit,nml=numerics,iostat=numerics_io)
    read_numerics = found(numerics_io)
    elements = nodes+1
    dx  = L_fin/elements
  end function
 
  pure real(rkind) function stable_dt(dx,alpha) 
    real(rkind) ,intent(in) :: dx,alpha
    real(rkind) ,parameter  :: safety_factor=0.9
    stable_dt = safety_factor*(dx**2/alpha)
  end function

  pure function differencing_stencil(dx,nodes) result(centralDiff)
    real(rkind)    ,intent(in) :: dx
    integer(ikind) ,intent(in) :: nodes
    real(rkind),dimension(:,:),allocatable :: centralDiff
    allocate(centralDiff(nodes,diagonals))
    centralDiff(:,low_diag)=  1./dx**2
    centralDiff(:,    diag)= -2./dx**2
    centralDiff(:, up_diag)=  1./dx**2 
  end function

  pure function differentiate(finiteDiff,T,T1st,Tlast)result(T_xx)
    real(rkind),dimension(:)  ,intent(in) ::T   
    real(rkind),dimension(:,:),intent(in) ::finiteDiff!differentiation op
    real(rkind)               ,intent(in) :: T1st,Tlast
    real(rkind),dimension(:)  ,allocatable:: T_xx
    integer(ikind) :: nodes,i
    nodes = size(T)
    allocate(T_xx(nodes))
    T_xx(1) =  finiteDiff(1,low_diag)*T1st &
              +finiteDiff(1,    diag)*T(1) &
              +finiteDiff(1, up_diag)*T(2)
    forall(i=2:nodes-1)
      T_xx(i) =  finiteDiff(i,low_diag)*T(i-1) &
                +finiteDiff(i,    diag)*T(i  ) &
                +finiteDiff(i, up_diag)*T(i+1) 
    end forall
    T_xx(nodes) =  finiteDiff(nodes,low_diag)*T(nodes-1) &
                  +finiteDiff(nodes,    diag)*T(nodes  ) &
                  +finiteDiff(nodes, up_diag)*Tlast
  end function 
end module 
