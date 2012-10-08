!*********************************************************************
! ForTrilinos: Object-Oriented Fortran 2003 interface to Trilinos
!                Copyright 2010 Sandia Corporation
!
! Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
! the U.S. Government retains certain rights in this software.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright
!    notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright
!    notice, this list of conditions and the following disclaimer in the
!    documentation and/or other materials provided with the distribution.
!
! 3. Neither the name of the Corporation nor the names of the
!    contributors may be used to endorse or promote products derived from
!    this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! Questions? Contact Karla Morris  (knmorri@sandia.gov)
!                    Damian Rouson (rouson@sandia.gov)
!*********************************************************************
module ForTrilinos_universal

  ! This module implements a base type that all ForTrilinos derived 
  ! types (except 'hermetic') extend.  It provides a universal dummy 
  ! argument class to which any actual argument can be passed in an 
  ! Epetra type-bound procedure.  The deferred binding "generalize" 
  ! ensures that each Epetra derived  type implements a type-bound 
  ! procedure that can be invoked to create an equivalent general entity
  ! of derived type ForTrilinos_Object_ID_t, which can then be converted
  ! to any other Epetra derived type.

  use ForTrilinos_hermetic ,only : hermetic
  use ForTrilinos_enums ,only : ForTrilinos_Object_ID_t
  implicit none
  type ,abstract ,public ,extends(hermetic) :: universal
  contains
    procedure(generalize_interface) ,deferred :: generalize 
  end type

  ! Implementations of this procedure will take in a specific struct and
  ! return a generic struct that can be used to call a casting function.

  abstract interface
    type(ForTrilinos_Object_ID_t) function generalize_interface(this)
      import :: universal,ForTrilinos_Object_ID_t
      class(universal) ,intent(in) ,target :: this
    end function
  end interface
end module
