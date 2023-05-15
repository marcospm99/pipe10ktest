! module type
! implicit none

! type tipo
!     double precision, allocatable :: im(:) 
!     double precision, allocatable :: re(:)
! end type tipo

! contains

!     subroutine allocata()
!     implicit none

!         type(tipo) :: A

!     allocate(A%re(10),A%im(10))
!     write(*,*) 'Matriz allocateada', A%re(0)

!     end subroutine

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     subroutine manipula()
!     implicit none

!         type(tipo) :: a

!         a%re(:) = 12
!         a%im(:) = 123
!         write(*,*) 'ha funcionado', a%re(1), a%im(1)

!     end subroutine

! end module

! !!!!!!!!!!!!!!!!!!!!!

! program main

! use type

! type(tipo) :: A

! call allocata()
! call manipula()

! end program

module type
  implicit none

  type tipo
      double precision, allocatable :: im(:) 
      double precision, allocatable :: re(:)
  end type tipo

contains

  function allocata() result(A)
    implicit none

    type(tipo) :: A

    allocate(A%re(10),A%im(10))
    write(*,*) 'Matriz allocateada', A%re(1)

  end function

  subroutine manipula(A)
    implicit none

    type(tipo) :: A
    A%re(:) = 12
    A%im(:) = 123
    write(*,*) 'ha funcionado', A%re(1), A%im(1)

  end subroutine

end module

program main
  use type
  type(tipo) :: A

  A = allocata()
  call manipula(A)

end program
