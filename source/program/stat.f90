#include "../parallel.h"
!*************************************************************************
   module std
!*************************************************************************
   use io
   implicit none
   double precision :: mean_ur(i_N), stdv_ur(i_N)
   double precision :: mean_ut(i_N), stdv_ut(i_N)
   double precision :: mean_uz(i_N), stdv_uz(i_N), stdv_rz(i_N)
   double precision :: d1, d(i_N)
   integer :: n,n_
   
   call vel_sta

   do n = 1, mes_D%pN
      n_ = mes_D%pNi + n - 1
      mean_ur(n_) = mean_ur(n_) + sum(vel_r%Re(:,:,n))
      stdv_ur(n_) = stdv_ur(n_) + sum(vel_r%Re(:,:,n)**2)
      mean_ut(n_) = mean_ut(n_) + sum(vel_t%Re(:,:,n))
      stdv_ut(n_) = stdv_ut(n_) + sum(vel_t%Re(:,:,n)**2)
      mean_uz(n_) = mean_uz(n_) + sum(vel_z%Re(:,:,n))
      stdv_uz(n_) = stdv_uz(n_) + sum(vel_z%Re(:,:,n)**2)
      stdv_rz(n_) = stdv_rz(n_) + sum(vel_z%Re(:,:,n)*vel_r%Re(:,:,n))
   end do

   call mpi_allreduce(mean_ur, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   mean_ur = d
   call mpi_allreduce(stdv_ur, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   stdv_ur = d
   call mpi_allreduce(mean_ut, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   mean_ut = d
   call mpi_allreduce(stdv_ut, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   stdv_ut = d
   call mpi_allreduce(mean_uz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   mean_uz = d
   call mpi_allreduce(stdv_uz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   stdv_uz = d

   call mpi_allreduce(stdv_rz, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   stdv_rz = d
   
   call mpi_barrier(mpi_comm_world, mpi_er)
   call mpi_finalize(mpi_er)

subroutine initializeSTD()
   mean_ur = 0d0
   stdv_ur = 0d0
   mean_ut = 0d0
   stdv_ut = 0d0
   mean_uz = 0d0
   stdv_uz = 0d0
   stdv_rz = 0d0
end


!*************************************************************************
 END MODULE
!*************************************************************************

