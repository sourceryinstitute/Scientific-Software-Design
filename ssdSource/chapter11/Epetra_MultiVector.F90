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

module Epetra_MultiVector_module ! ---------- Excerpt --------------
  use ForTrilinos_enums ,only: FT_Epetra_MultiVector_ID_t,&
      FT_Epetra_Map_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_universal,only:universal
  use ForTrilinos_error
  use FEpetra_BlockMap  ,only: Epetra_BlockMap
  use iso_c_binding     ,only: c_int,c_double,c_char
  use forepetra
  implicit none
  private                      ! Hide everything by default
  public :: Epetra_MultiVector ! Expose type/constructors/methods
  type ,extends(universal)           :: Epetra_MultiVector !"shell"
    private
    type(FT_Epetra_MultiVector_ID_t) :: MultiVector_id
  contains
     procedure         :: cpp_delete => ctrilinos_delete_EpetraMultiVector
     procedure         :: get_EpetraMultiVector_ID 
     procedure ,nopass :: alias_EpetraMultiVector_ID
     procedure         :: generalize 
     procedure         :: Random
     procedure         :: Norm2
     procedure         :: NumVectors
  end type

   interface Epetra_MultiVector ! constructors
     module procedure from_scratch,duplicate,from_struct
   end interface

contains

  type(Epetra_MultiVector) function from_struct(id)
     type(FT_Epetra_MultiVector_ID_t) ,intent(in) :: id
     from_struct%MultiVector_id = id  
     call from_struct%register_self
  end function

  type(Epetra_MultiVector) function from_scratch( &
     BlockMap,Num_Vectors,zero &
  )
   use ForTrilinos_enums ,only: FT_boolean_t,FT_TRUE,FT_FALSE
   use iso_c_binding     ,only: c_int
   class(Epetra_BlockMap) ,intent(in) :: BlockMap
   integer(c_int)         ,intent(in) :: Num_Vectors
   logical                ,intent(in) :: zero 
   integer(FT_boolean_t)              :: zero_in 
   type(FT_Epetra_MultiVector_ID_t)   :: from_scratch_id
   if (zero) zero_in=FT_TRUE
   if (.not.zero) zero_in=FT_FALSE
   from_scratch_id = Epetra_MultiVector_Create( &
     BlockMap%get_EpetraBlockMap_ID(),Num_Vectors,zero_in &
   )
   from_scratch = from_struct(from_scratch_id)
  end function

  type(Epetra_MultiVector) function duplicate(this)
    type(Epetra_MultiVector) ,intent(in) :: this
    type(FT_Epetra_MultiVector_ID_t) :: duplicate_id
    duplicate_id = Epetra_MultiVector_Duplicate(this%MultiVector_id)
    duplicate = from_struct(duplicate_id)
  end function

  type(FT_Epetra_MultiVector_ID_t) function get_EpetraMultiVector_ID( &
      this &
  )
    class(Epetra_MultiVector) ,intent(in) :: this 
    get_EpetraMultiVector_ID=this%MultiVector_id
  end function
  
  type(FT_Epetra_MultiVector_ID_t) function alias_EpetraMultiVector_ID( &
    generic_id &
  )
    use ForTrilinos_table_man,only: CT_Alias
    use iso_c_binding        ,only: c_loc,c_int
    use ForTrilinos_enums    ,only: ForTrilinos_Universal_ID_t, &
                                    FT_Epetra_MultiVector_ID
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,allocatable ,target :: alias_id
    integer(c_int) :: status
    type(error) :: ierr
    allocate(alias_id,source=CT_Alias( & 
      generic_id,FT_Epetra_MultiVector_ID),stat=status &
    )
    ierr=error(status,'FEpetra_MultiVector:alias_EpetraMultiVector_ID')
    call ierr%check_success()
    alias_EpetraMultiVector_ID=degeneralize_EpetraMultiVector( &
      c_loc(alias_id) &
    )
  end function

  type(ForTrilinos_Universal_ID_t) function generalize(this)
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding     ,only : c_loc
   class(Epetra_MultiVector) ,intent(in) ,target :: this
   generalize = generalize_all(c_loc(this%MultiVector_id))
  end function

 type(FT_Epetra_MultiVector_ID_t) function degeneralize_EpetraMultiVector( &
    generic_id) bind(C)
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t, &
                                  FT_Epetra_MultiVector_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)                      ,value   :: generic_id
    type(FT_Epetra_MultiVector_ID_t) ,pointer :: local_ptr=>null()
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraMultiVector = local_ptr
  end function
 
  subroutine Random(this,err)
   class(Epetra_MultiVector) ,intent(inout) :: this
   type(error)     ,optional ,intent(out)   :: err
   integer(c_int)                           :: error_out
   error_out = Epetra_MultiVector_Random (this%MultiVector_id)
   if (present(err)) err=error(error_out)
  end subroutine

  function Norm2(this,err) result(Norm2_val)
    class(Epetra_MultiVector)   ,intent(in)  :: this
    type(error)    ,optional    ,intent(out) :: err
    real(c_double) ,dimension(:),allocatable :: Norm2_val 
    integer(c_int)                           :: error_out
    integer(c_int)                           :: status
    type(error)                              :: ierr
    if (.not.allocated(Norm2_val)) then
      allocate(Norm2_val(this%NumVectors()),stat=status)
      ierr=error(status,'FEpetra_MultiVector:Norm2')
      call ierr%check_success()
    endif
    error_out = Epetra_MultiVector_Norm2(this%MultiVector_id,Norm2_val)
    if (present(err)) err=error(error_out)
  end function
 
  integer(c_int) function NumVectors(this)
    class(Epetra_MultiVector) ,intent(in) :: this
    NumVectors=Epetra_MultiVector_NumVectors(this%MultiVector_id)
  end function 

  subroutine ctrilinos_delete_EpetraMultiVector(this)
    class(Epetra_MultiVector),intent(inout) :: this
    call Epetra_MultiVector_Destroy( this%MultiVector_id ) 
  end subroutine

end module 

