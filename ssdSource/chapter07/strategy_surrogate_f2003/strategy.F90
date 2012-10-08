module strategy_module
  ! Substitute for integrand (avoiding circular references)
  use surrogate_module ,only : surrogate

  implicit none
  private

  ! Abstract time integration strategy
  type ,abstract ,public :: strategy
  contains
    ! Abstract integration procedure interface
    procedure(integrator_interface) ,nopass ,deferred :: integrate
  end type

  abstract interface
    subroutine integrator_interface(this,dt)
      import :: surrogate
      class(surrogate) ,intent(inout) :: this ! integrand
      real             ,intent(in)    :: dt ! time step size
    end subroutine
  end interface
end module
