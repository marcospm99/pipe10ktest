


!*************************
!
!
!*************************
#include "../parallel.h"
 module mpif
 
 use mpi
!*************************
   implicit none
   save

   integer :: mpi_er, mpi_tg, mpi_rq(0:2*_Np), mpi_st(mpi_status_size)
   integer :: mpi_rnk, mpi_sze

 contains

!------------------------------------------------------------------------
!  initialise 
!------------------------------------------------------------------------
   subroutine mpi_precompute()
      
      call mpi_init(mpi_er)
      call mpi_comm_rank(mpi_comm_world, mpi_rnk, mpi_er)
      call mpi_comm_size(mpi_comm_world, mpi_sze, mpi_er)
      if(mpi_sze /= _Np) stop 'mpi_precompute: incorrect num procs'
      
   end subroutine mpi_precompute

!*************************
 end module mpif
!*************************