!**************************************************************************
!  IN/OUT 
!
!**************************************************************************
#include "../parallel.h"
 module sta
!**************************************************************************
   use velocity
   use h5lt
   use hdf5

   implicit none
   save


   double precision, private :: uclm, utaum,ucl, utau, Ub

! ------------------------- stats  -------------------------------

   double precision :: mean_ur(i_N), stdv_ur(i_N)
   double precision :: mean_ut(i_N), stdv_ut(i_N)
   double precision :: mean_uz(i_N), stdv_uz(i_N), stdv_rz(i_N)

   double precision :: d(i_N) ! auxiliary mem
   integer :: csta

! ------------------------- HDF5 -------------------------------

   integer:: info,ierr
   integer(hid_t):: fid,pid, dset_id,dspace_id
   integer:: h5err
   integer(hsize_t),dimension(3):: dims
   integer(hsize_t),dimension(1):: hdims,maxdims
   integer(hsize_t),dimension(2):: hdims2
   
   integer(hid_t) :: header_id,sta_id ! Group identifiers
   integer(size_t)::siever
   parameter (siever = 4*1024*1024)


! --------------------- Program -------------------------------------

 contains
! Compute statistics
   subroutine compute_sta()
   implicit none
   integer:: n, n_

   ! Compute friction velocity

   if (mpi_rnk ==0 ) then
      ucl = 1d0 + dot_product(vel_uz%Re(1:1+i_KL,0),mes_D%dr0(:,0))
      utau = dot_product(vel_uz%Re(i_N-i_KL:i_N,0),mes_D%dr1(:,1))
      utau = dsqrt(dabs((Utau-2d0)/d_Re))

      ! add statistics
      uclm = uclm + ucl
      utaum = utaum + utau 
   endif

   call vel_sta() ! Compute vel_r
   ! call vel_diss

!  call diss_sta1()
!  call diss_sta2()
!  call diss_sta3()

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
   csta = csta + 1

   call mpi_barrier(mpi_comm_world, mpi_er)

end subroutine

subroutine initialiseSTD()
implicit none
   csta = 0
   utaum   = 0d0
   uclm    = 0d0
   mean_ur = 0d0
   stdv_ur = 0d0
   mean_ut = 0d0
   stdv_ut = 0d0
   mean_uz = 0d0
   stdv_uz = 0d0
   stdv_rz = 0d0
end

subroutine saveStats(fnameima)
implicit none

    character(len = 256):: fnameima
    call mpi_reduce(mean_ur, d, i_N, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    mean_ur = d
    call mpi_reduce(stdv_ur, d, i_N, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    stdv_ur = d
    call mpi_reduce(mean_ut, d, i_N, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    mean_ut = d
    call mpi_reduce(stdv_ut, d, i_N, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    stdv_ut = d
    call mpi_reduce(mean_uz, d, i_N, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    mean_uz = d
    call mpi_reduce(stdv_uz, d, i_N, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    stdv_uz = d
    call mpi_reduce(stdv_rz, d, i_N, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    stdv_rz = d

    if(mpi_rnk==0)  then

       call h5fcreate_f(trim(fnameima),H5F_ACC_TRUNC_F,fid,h5err)
       call h5gcreate_f(fid, '/header', header_id, h5err)
       call h5gcreate_f(fid, '/sta'   , sta_id   , h5err)

       hdims = (/1/)
       call h5ltmake_dataset_double_f(header_id,"time",1,hdims,(/tim_t/),h5err)
       call h5ltmake_dataset_double_f(header_id,"Re",1,hdims,(/d_Re/),h5err)
       call h5ltmake_dataset_double_f(header_id,"alpha",1,hdims,(/d_alpha/),h5err)

       call h5ltmake_dataset_int_f(header_id,"N" ,1,hdims,(/i_N/),h5err)
       call h5ltmake_dataset_int_f(header_id,"num" ,1,hdims,(/csta/),h5err)
       call h5ltmake_dataset_int_f(header_id,"M" ,1,hdims,(/i_M/),h5err)
       call h5ltmake_dataset_int_f(header_id,"K" ,1,hdims,(/i_K/),h5err)
       call h5ltmake_dataset_int_f(header_id,"Mp",1,hdims,(/i_Mp/),h5err)

       call h5ltmake_dataset_double_f(sta_id,"utau",1,hdims,utaum,h5err)
       call h5ltmake_dataset_double_f(sta_id,"ucl",1,hdims,uclm,h5err)

       call h5ltmake_dataset_double_f(header_id,"dt"   ,1,hdims,(/tim_dt/),h5err)

       hdims = (/i_N/)
       call h5ltmake_dataset_double_f(header_id,"r"   ,1,hdims,mes_D%r(1:i_N,1),h5err)

       call h5ltmake_dataset_double_f(sta_id,"mean_ur",1,hdims,mean_ur,h5err)
       call h5ltmake_dataset_double_f(sta_id,"mean_uz",1,hdims,mean_uz,h5err)
       call h5ltmake_dataset_double_f(sta_id,"mean_ut",1,hdims,mean_ut,h5err)

       call h5ltmake_dataset_double_f(sta_id,"stdv_ur",1,hdims,stdv_ur,h5err)
       call h5ltmake_dataset_double_f(sta_id,"stdv_ut",1,hdims,stdv_ut,h5err)
       call h5ltmake_dataset_double_f(sta_id,"stdv_uz",1,hdims,stdv_uz,h5err)
       call h5ltmake_dataset_double_f(sta_id,"stdv_rz",1,hdims,stdv_rz,h5err)

       call h5gclose_f(header_id,h5err)
       call h5gclose_f(sta_id,h5err)

       call h5fclose_f(fid,h5err)
   endif
   ! Once saved, initialise e verything to 0 again. 
   call initialiseSTD()

end subroutine

end module sta

