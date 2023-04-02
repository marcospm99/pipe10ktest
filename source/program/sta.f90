!**************************************************************************
!  IN/OUT 
!
!**************************************************************************
#include "../parallel.h"
 module sta
!**************************************************************************
   use variables
   use velocity
   use h5lt
   use hdf5
   use transform

   implicit none
   save


   double precision, private :: uclm, utaum,ucl, utau, Ub
   double precision, private :: aa, bb, cc

! ------------------------- stats  -------------------------------

   double precision :: mean_ur(i_N), stdv_ur(i_N)
   double precision :: mean_ut(i_N), stdv_ut(i_N)
   double precision :: mean_uz(i_N), stdv_uz(i_N), stdv_rz(i_N), stdv_rt(i_N), stdv_tz(i_N)
   double precision :: mom_ur(i_N,10)
   double precision :: mom_ut(i_N,10)
   double precision :: mom_uz(i_N,10)
   double precision :: urf(i_N), utf(i_N), uzf(i_N)
   double precision :: pfpF(i_N)
   double precision :: or(i_N), ot(i_N), oz(i_N)
   

   double precision :: dissr(i_N,3),disst(i_N,3),dissz(i_N,3), diss(i_N,3), dzduzsq(i_N), dzduzcub(i_N)
   double precision :: tmpr1(i_N),tmpr2(i_N),tmpt1(i_N),tmpt2(i_N),tmpz1(i_N),tmpz2(i_N)
   double precision :: factor

   double precision :: d(i_N),dd(i_n,10) ! auxiliary mem
   integer :: csta

   !type (coll), private  :: c1,c2,c3 ! Three colls are defined here. Why! They are really big. 
                                     ! They are defined as private. They cannot be used anywhere else
                                     ! Remove in future versions. We have to pass the routine some workarray 

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
   integer:: n, n_, ii, jj, kk 

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

   call staFFT() ! Compute stdv_[u,t,z]
   call var_coll_dissp() ! Compute dissipation
   call vort()
   call pressure()





      

   do n = 1, mes_D%pN
      n_ = mes_D%pNi + n - 1
      mean_ur(n_) = mean_ur(n_) + sum(vel_r%Re(:,:,n))
      stdv_ur(n_) = stdv_ur(n_) + sum(vel_r%Re(:,:,n)**2)
      mean_ut(n_) = mean_ut(n_) + sum(vel_t%Re(:,:,n))
      stdv_ut(n_) = stdv_ut(n_) + sum(vel_t%Re(:,:,n)**2)
      mean_uz(n_) = mean_uz(n_) + sum(vel_z%Re(:,:,n))
      stdv_uz(n_) = stdv_uz(n_) + sum(vel_z%Re(:,:,n)**2)
      stdv_rz(n_) = stdv_rz(n_) + sum(vel_z%Re(:,:,n)*vel_r%Re(:,:,n))
      stdv_rt(n_) = stdv_rz(n_) + sum(vel_r%Re(:,:,n)*vel_t%Re(:,:,n))
      stdv_tz(n_) = stdv_rz(n_) + sum(vel_t%Re(:,:,n)*vel_z%Re(:,:,n))

      do kk =  0,i_Th-1
         do jj =  0,i_pZ-1
            aa = 1d0
            bb = 1d0
            cc = 1d0
            do ii = 1,10
               aa = aa * vel_r%Re(jj,kk,n)
               bb = bb * vel_t%Re(jj,kk,n)
               cc = cc * vel_z%Re(jj,kk,n)
               mom_ur(n_,ii) =  mom_ur(n_,ii) + aa
               mom_ut(n_,ii) =  mom_ut(n_,ii) + bb
               mom_uz(n_,ii) =  mom_uz(n_,ii) + cc
            enddo
         enddo
      enddo
   end do
   csta = csta + 1

   call mpi_barrier(mpi_comm_world, mpi_er)

end subroutine compute_sta


! !------------------------------------------------------------------------
! !  Pressure field:
! !      Results: 
! !      
! !------------------------------------------------------------------------
subroutine pressure()
implicit none

type (coll) :: pfp !pressure field perturbation
double precision :: factor
_loop_km_vars

call vel_adjPPE(1) ! compute p'
pfp = c4


! fourier interp

_loop_km_begin

       tmpr1 =  pfp%Re(:,nh)
       tmpr2 =  pfp%Im(:,nh)

factor = 2d0
   if (m==0) factor = 1d0

          pfpF(:) = pfpF(:) + factor*(tmpr1(:)**2+tmpr2(:)**2)
          
_loop_km_end
  

end subroutine pressure



! !------------------------------------------------------------------------
! !  Turbulent budgets:
! !      Dissipation, % optimization pending
! !      Rest of the derivatives: budgets in postproc
! !------------------------------------------------------------------------
 subroutine var_coll_dissp()
    implicit none

       double precision :: factor
       integer :: n
       _loop_km_vars


       !!! Empezamos para uz !!!


       call var_coll_copy(vel_uz,c1)
       call var_coll_meshmult(0,mes_D%dr(1),c1, c2) 

       _loop_km_begin
       tmpr1 =  c2%Re(:,nh)
       tmpr2 =  c2%Im(:,nh)

       tmpz1 = -c1%Im(:,nh)*ad_k1a1(k)
       tmpz2 =  c1%Re(:,nh)*ad_k1a1(k)

       tmpt1 = -c1%Im(:,nh)*ad_m1r1(:,m)
       tmpt2 =  c1%Re(:,nh)*ad_m1r1(:,m)
     
       factor = 2d0
       if (m==0) factor = 1d0

          dissr(:,3) = dissr(:,3) + factor*(tmpr1(:)**2+tmpr2(:)**2)
          dissz(:,3) = dissz(:,3) + factor*(tmpz1(:)**2+tmpz2(:)**2)
          disst(:,3) = disst(:,3) + factor*(tmpt1(:)**2+tmpt2(:)**2)
     
       _loop_km_end
  


!       ! Lo mismo para theta

      call var_coll_copy(vel_ut,c1)
      call var_coll_copy(vel_ur,c3)
      call var_coll_meshmult(1,mes_D%dr(1),c1, c2) 

      _loop_km_begin

      tmpr1 =  c2%Re(:,nh)
      tmpr2 =  c2%Im(:,nh)

      tmpz1 = -c1%Im(:,nh)*ad_k1a1(k)
      tmpz2 =  c1%Re(:,nh)*ad_k1a1(k)

      tmpt1 = -c1%Im(:,nh)*ad_m1r1(:,m) + c3%Re(:,nh)*mes_D%r(:,-1) !+
      tmpt2 =  c1%Re(:,nh)*ad_m1r1(:,m) + c3%Im(:,nh)*mes_D%r(:,-1) !+

      factor = 2d0
      if (m==0) factor = 1d0
      
         dissr(:,2) = dissr(:,2) + factor*(tmpr1(:)**2+tmpr2(:)**2)
         dissz(:,2) = dissz(:,2) + factor*(tmpz1(:)**2+tmpz2(:)**2)
         disst(:,2) = disst(:,2) + factor*(tmpt1(:)**2+tmpt2(:)**2)

      _loop_km_end

      
      ! Lo mismo para r

      call var_coll_copy(vel_ur,c1)
      call var_coll_copy(vel_ut,c3)
      call var_coll_meshmult(1,mes_D%dr(1),c1, c2)

      _loop_km_begin
      tmpr1 =  c2%Re(:,nh)
      tmpr2 =  c2%Im(:,nh)

      tmpz1 = -c1%Im(:,nh)*ad_k1a1(k)
      tmpz2 =  c1%Re(:,nh)*ad_k1a1(k)

      tmpt1 = -c1%Im(:,nh)*ad_m1r1(:,m) - c3%Re(:,nh)*mes_D%r(:,-1) !-
      tmpt2 =  c1%Re(:,nh)*ad_m1r1(:,m) - c3%Im(:,nh)*mes_D%r(:,-1) !-


      factor = 2d0
      if (m==0) factor = 1d0

         dissr(:,1) = dissr(:,1) + factor*(tmpr1(:)**2+tmpr2(:)**2)
         dissz(:,1) = dissz(:,1) + factor*(tmpz1(:)**2+tmpz2(:)**2)
         disst(:,1) = disst(:,1) + factor*(tmpt1(:)**2+tmpt2(:)**2)

   _loop_km_end


   !dzduzsq

   call var_coll_copy(vel_uz,c1)

   _loop_km_begin
   tmpr1 = (c1%Re(:,nh)**2 + c1%Im(:,nh)**2) *ad_k1a1(k)
   tmpr2 = (2*c1%Re(:,nh)*c1%Im(:,nh)) *ad_k1a1(k)

         factor = 2d0
      if (m==0) factor = 1d0

      dzduzsq(:) = dzduzsq(:) + factor*(tmpr1(:)**2+tmpr2(:)**2)


   _loop_km_end
   
     !dzduzcub

   call var_coll_copy(vel_uz,c1)

   _loop_km_begin
   tmpr1 = (c1%Re(:,nh)**3 - 3*c1%Im(:,nh)**2*c1%Re(:,nh)) *ad_k1a1(k)
   tmpr2 = (3*c1%Im(:,nh)*c1%Re(:,nh)**2 - c1%Im(:,nh)**3) *ad_k1a1(k)

         factor = 2d0
      if (m==0) factor = 1d0

      dzduzcub(:) = dzduzcub(:) + factor*(tmpr1(:)**2+tmpr2(:)**2)


   _loop_km_end 


 end subroutine var_coll_dissp


!------------------------------------------------------------------------
!  Velocities using fourier, % optimization pending
!------------------------------------------------------------------------
 subroutine staFFT()
    implicit none

    double precision :: factor
    integer :: n, comp, n_
    _loop_km_vars

    !call var_coll_copy(vel_ur,c1)

    _loop_km_begin

     factor = 2d0
     if (m==0) then;
        factor = 1d0
     end if

     urf(:) = urf(:) + factor*(vel_ur%Re(:,nh)**2+vel_ur%Im(:,nh)**2)
     utf(:) = utf(:) + factor*(vel_ut%Re(:,nh)**2+vel_ut%Im(:,nh)**2)
     uzf(:) = uzf(:) + factor*(vel_uz%Re(:,nh)**2+vel_uz%Im(:,nh)**2)

     _loop_km_end

     end subroutine staFFT

   !!!!!!!!!!!!!
   ! Vorticity !
   !!!!!!!!!!!!!

     subroutine vort()

      implicit none
       _loop_km_vars
 
      call var_coll_curl(vel_ur,vel_ut,vel_uz,c1,c2,c3)

      _loop_km_begin

     factor = 2d0
     if (m==0) factor = 1d0

     or(:) = or(:) + sqrt(factor*(c1%Re(:,nh)**2+c1%Im(:,nh)**2))
     ot(:) = ot(:) + sqrt(factor*(c2%Re(:,nh)**2+c2%Im(:,nh)**2))
     oz(:) = oz(:) + sqrt(factor*(c3%Re(:,nh)**2+c3%Im(:,nh)**2))

      _loop_km_end

     end subroutine vort







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
   stdv_rt = 0d0
   stdv_tz = 0d0

   mom_ur = 0d0
   mom_uz = 0d0
   mom_ut = 0d0

   dissr = 0d0
   disst = 0d0
   dissz = 0d0

   urf = 0d0
   pfpF = 0d0

   diss = 0d0
   dzduzsq = 0d0
   dzduzcub = 0d0
end subroutine initialiseSTD


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
    call mpi_reduce(stdv_rt, d, i_N, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    stdv_rt = d
    call mpi_reduce(stdv_tz, d, i_N, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    stdv_tz = d
    call mpi_reduce(pfpF, d, i_N, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    pfpF = d

   !! Dissipation
       diss(:,1) = dissr(:,1) + disst(:,1) + dissz(:,1)
       diss(:,2) = dissr(:,2) + disst(:,2) + dissz(:,2)
       diss(:,3) = dissr(:,3) + disst(:,3) + dissz(:,3)

   ! call mpi_reduce(diss, dissr, 3*i_N, mpi_double_precision,  & ! Usammos dissr que ya está usada y el no. de elementos a escribir es 3*i_N
   !     mpi_sum, 0, mpi_comm_world, mpi_er)
   !  diss = dissr

   call mpi_reduce(diss(1,1), d, i_N, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    diss(:,1) = d

    call mpi_reduce(diss(1,2), d, i_N, mpi_double_precision,  &
    mpi_sum, 0, mpi_comm_world, mpi_er)
 diss(:,2) = d
    call mpi_reduce(diss(1,3), d, i_N, mpi_double_precision,  &
    mpi_sum, 0, mpi_comm_world, mpi_er)
 diss(:,3) = d

   ! call mpi_reduce(urf, d, i_N, mpi_double_precision,  &
   !     mpi_sum, 0, mpi_comm_world, mpi_er)
   !  urf = d

   ! or, ot, oz: vort's
   call mpi_reduce(or, d, i_N, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    or = d
   call mpi_reduce(ot, d, i_N, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    ot = d
   call mpi_reduce(oz, d, i_N, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    oz = d

   call mpi_reduce(dzduzsq, d, i_N, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    dzduzsq = d
       call mpi_reduce(dzduzcub, d, i_N, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    dzduzcub = d

    call mpi_reduce(mom_ur, dd, i_N*10, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    mom_ur = dd
    call mpi_reduce(mom_uz, dd, i_N*10, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    mom_uz = dd
    call mpi_reduce(mom_ut, dd, i_N*10, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    mom_ut = dd

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
       call h5ltmake_dataset_double_f(sta_id,"stdv_rt",1,hdims,stdv_rt,h5err)
       call h5ltmake_dataset_double_f(sta_id,"stdv_tz",1,hdims,stdv_tz,h5err)

       call h5ltmake_dataset_double_f(sta_id,"urf",1,hdims,urf,h5err)
       call h5ltmake_dataset_double_f(sta_id,"pfpF",1,hdims,pfpF,h5err)

       
       call h5ltmake_dataset_double_f(sta_id,"disstt",1,hdims,diss(:,2),h5err)
       call h5ltmake_dataset_double_f(sta_id,"dissrr",1,hdims,diss(:,1),h5err)
       call h5ltmake_dataset_double_f(sta_id,"disszz",1,hdims,diss(:,3),h5err)       
       
       call h5ltmake_dataset_double_f(sta_id,"vortr",1,hdims,or,h5err)
       call h5ltmake_dataset_double_f(sta_id,"vortt",1,hdims,ot,h5err)
       call h5ltmake_dataset_double_f(sta_id,"vortz",1,hdims,oz,h5err)  

       call h5ltmake_dataset_double_f(sta_id,"dzduzsq",1,hdims,dzduzsq,h5err)
       call h5ltmake_dataset_double_f(sta_id,"dzduzcub",1,hdims,dzduzcub,h5err)
         

       hdims2 = (/i_N,10/)

       call h5ltmake_dataset_double_f(sta_id,"mom_ur",2,hdims2,mom_ur,h5err)
       call h5ltmake_dataset_double_f(sta_id,"mom_uz",2,hdims2,mom_uz,h5err)
       call h5ltmake_dataset_double_f(sta_id,"mom_ut",2,hdims2,mom_ut,h5err)

       call h5gclose_f(header_id,h5err)
       call h5gclose_f(sta_id,h5err)

       call h5fclose_f(fid,h5err)
   endif
   ! Once saved, initialise e verything to 0 again. 
   call initialiseSTD()

end subroutine saveStats

end module sta

