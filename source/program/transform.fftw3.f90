!*************************************************************************
! Transformation using the FFT, de-aliased.
! Written for FFTW version 3 -- documentation at  www.fftw.org
! See FFTW reference, sections 2.3, 4.3.3, 4.4.1-2
!
! \sum_{-Ma}^{Ma} == \sum_0^{3M-1}
! logical size of DFT n==3M=Th
!
!*************************************************************************
#include "../parallel.h"
 module transform
!*************************************************************************
   
   use wksp
   use parameters
   use variables
   
   implicit none
   save

   integer, parameter, private :: i_3K = 3*i_K
   integer, parameter, private :: i_3M = 3*i_M
   integer, parameter, private :: i_Ma = (3*i_M)/2
   double complex,     private :: X(0:i_3K-1, 0:_Ms-1)   ! Mem = 512*3*128/8*72*16/1024**3 no preocupante, 0.026 GB
   double complex,     private :: Y(0:i_3K-1, 0:_Ms-1)   ! Idem
   double complex,     private :: Xs(0:i_pZ-1, 0:i_Ma)   ! Mem = 3*512/8*3/2*128*16*72/1024**3, no preocupante, 0.038 GB
   double precision,   private :: Ys(0:i_pZ-1, 0:i_3M-1) ! Mem = 3*512/8*3*128*8*72/1024**3, no preocupante, 0.038 GB
   integer*8,          private :: plan_c2cf, plan_c2cb, plan_r2c, plan_c2r
      integer :: nh, n,m, m_
      double precision :: scale_
   

   double complex,     allocatable :: Taux(:,:,:) !Taux(0:i_3K-1, 0:_Ms-1, i_pN) ! Mem = 3*512*128/8*72*384/9*16/1024**3, Preocupante, 1.12 GB total
   
   
   !  double complex,     private :: Ts(0:i_pZ-1, 0:i_M1, i_pN) ! Mem = 2*8*154*127*40*80/1024**3, Preocupante, 0.93 GB total
  

 contains


!------------------------------------------------------------------------
! Setup plans for transforms.  
!------------------------------------------------------------------------
   subroutine tra_precompute()
      implicit none
      integer, parameter :: flag=32 !=FFTW_PATIENT see fftw3.f
      integer :: sgn, n(1), howmany, inembed(1), onembed(1)

         ! if(.not.allocated(s1%Re))  allocate(s1%Re(0:_Hs1, i_pN))
         ! if(.not.allocated(s1%Im))  allocate(s1%Im(0:_Hs1, i_pN))

         if(.not.allocated(Taux))  allocate(Taux(0:i_3K-1, 0:_Ms-1, i_pN))

      
      n = (/i_3K/)            
      howmany = _Ms
      inembed = (/i_3K*_Ms/)
      onembed = (/i_3K*_Ms/)
      sgn = 1
      call dfftw_plan_many_dft(plan_c2cf, 1, n, howmany,  &
         X, inembed, 1, i_3K,  Y, onembed, 1, i_3K,  sgn, flag)
      sgn = -1
      call dfftw_plan_many_dft(plan_c2cb, 1, n, howmany,  &
         Y, onembed, 1, i_3K,  X, inembed, 1, i_3K,  sgn, flag)

      n = (/i_3M/)
      howmany = i_pZ
      inembed = (/i_pZ*(i_Ma+1)/)
      onembed = (/i_pZ*i_3M/)
      call dfftw_plan_many_dft_c2r(plan_c2r, 1, n, howmany,  &
         Xs, inembed, i_pZ, 1,  Ys, onembed, i_pZ, 1,  flag)
      call dfftw_plan_many_dft_r2c(plan_r2c, 1, n, howmany,  &
         Ys, onembed, i_pZ, 1,  Xs, inembed, i_pZ, 1,  flag)
      
   end subroutine tra_precompute


!-------------------------------------------------------------------------
!  convert collocated -> physical space
!-------------------------------------------------------------------------
   subroutine tra_coll2phys(c1,p1, c2,p2, c3,p3)
   
      implicit none
      
      type (coll), intent(inout)  :: c1,c2,c3
      type (phys), intent(inout) :: p1,p2,p3

      call var_coll2spec(c1,s1) 

      call tra_spec2phys(s1,p1)

      call var_coll2spec(c2,s1)
      call tra_spec2phys(s1,p2)

      call var_coll2spec(c3,s1)
      call tra_spec2phys(s1,p3)
      
   end subroutine tra_coll2phys


   subroutine tra_coll2phys1d(c,p)
   type (coll), intent(inout)  :: c
   type (phys), intent(inout) :: p

   call var_coll2spec(c,s1)
   call tra_spec2phys(s1,p)

   end subroutine tra_coll2phys1d


!-------------------------------------------------------------------------
!  convert collocated -> physical space
!-------------------------------------------------------------------------
   subroutine tra_phys2coll(p,c, p2,c2, p3,c3)
   
      implicit none
      type (phys), intent(in)  :: p,p2,p3
      type (coll), intent(inout) :: c,c2,c3

      call tra_phys2spec(p,s1)
      call var_spec2coll(s1,c)

      call tra_phys2spec(p2,s1)
      call var_spec2coll(s1,c2)

      call tra_phys2spec(p3, s1)
      call var_spec2coll(s1,c3)

      
   end subroutine tra_phys2coll

   subroutine tra_phys2coll1d(p,c)
      type (coll), intent(inout)  :: C
      type (phys), intent(in) :: p

      call tra_phys2spec(p,  s1)
      
      call var_spec2coll(s1, c)

   end subroutine tra_phys2coll1d


!------------------------------------------------------------------------
!  Convert spectral to real space
!------------------------------------------------------------------------
   subroutine tra_spec2phys(s, p)

      implicit none
      type (spec), intent(in)  :: s
      type (phys), intent(inout) :: p

      integer :: nh, n,m,m_
      				! for each r_n ...   
   



      do n = 1, mes_D%pN
         if(mpi_rnk/_Nr==0) then
            X(0:i_K1,   0) = dcmplx(s%Re(0:i_K1,n),s%Im(0:i_K1,n))
            X(i_K:2*i_K,0) = 0d0
            X(2*i_K+1:, 0) = dcmplx(s%Re(i_K1:1:-1,n),-s%Im(i_K1:1:-1,n))
            m_ = 1
            nh = 2*i_K-1
         else
            m_ = 0
            nh = i_K1
         end if
         do m = m_, _Ms1
            X(0:i_K1,   m) = dcmplx(s%Re(nh:nh+i_K1,n),s%Im(nh:nh+i_K1,n))
            X(i_K:2*i_K,m) = 0d0
            X(2*i_K+1:, m) = dcmplx(s%Re(nh-i_K1:nh-1,n),s%Im(nh-i_K1:nh-1,n))
            nh = nh + 2*i_K-1
         end do
         call dfftw_execute(plan_c2cf)
         
         Taux(:,:,n) = Y
      end do
      call tra_T2Ts()
      
      do n = 1, mes_D%pN
         Xs(:,0:i_M1) = Ts(:,:,n) 
                 
! #endif
         Xs(:,i_M:) = 0d0
         call dfftw_execute(plan_c2r)
         p%Re(:,:,n) = Ys
      end do
   
   end subroutine tra_spec2phys


!------------------------------------------------------------------------
!  Convert real to spectral space
!------------------------------------------------------------------------
   subroutine tra_phys2spec(p, s)

      implicit none
      type (phys), intent (in)  :: p
      type (spec), intent (inout) :: s ! Había que ponerlo en inout
      integer :: nh, n,m, m_
      double precision :: scale_


      ! if(mpi_rnk==0) write(*,*) 'master: el tamaño es ', size(s1%Re,1),  'x', size(s1%Re,2)
         			! scale, FFTW 4.7.2
      scale_ = 1d0 / dble(i_3K*i_3M)

      
      
      do n = 1, mes_D%pN
         Ys = scale_ * p%Re(:,:,n)
         call dfftw_execute(plan_r2c)

         Ts(:,:,n) = Xs(:,0:i_M1)
      end do
      
      
      call tra_Ts2T()
      
      
      do n = 1, mes_D%pN
         Y = Taux(:,:,n)
         
         call dfftw_execute(plan_c2cb)
         if(mpi_rnk/_Nr==0) then
            
            s%Re(0:i_K1,n) =  dble(X(0:i_K1,0))
            s%Im(0:i_K1,n) = dimag(X(0:i_K1,0))
            m_ = 1
            nh = 2*i_K-1
            
            ! write(*,*) 'master: el tamaño es ', size(s1%Re,1),  'x', size(s1%Re,2)
            
            
         else
            
            m_ = 0
            nh = i_K1
            ! if (nh > size(s1%Re, 1) .or. nh > size(s1%Im, 1)) write(*,*) 'Error: nh out of bounds'
            !  write(*,*) 'slave: el tamaño es ', size(s1%Re,1),  'x', size(s1%Re,2)
         end if

         
         
         do m = m_,_Ms1
            ! if (m==m_+1) write(*,*) 'hasta aquí' ! No llega a la segunda iteración de m, problema en nh
            ! if (m==m_) write(*,*) 'hola'  

            s%Re(nh:nh+i_K1,n) =  dble(X(0:i_K1,m))
            s%Im(nh:nh+i_K1,n) = dimag(X(0:i_K1,m))
            s%Re(nh-i_K1:nh-1,n) =  dble(X(2*i_K+1:,m))
            s%Im(nh-i_K1:nh-1,n) = dimag(X(2*i_K+1:,m))
            
            
            nh = nh + 2*i_K-1

            
            

         end do
         

      end do

      
         
   end subroutine tra_phys2spec


!------------------------------------------------------------------------
! transposes
!------------------------------------------------------------------------
! #if _Ns != 1
   subroutine tra_T2Ts()
 
           implicit none

      integer :: stp, dst,src, l,j, rnk,rko
      integer :: n,m, pm0,jz0
 

      rnk = mpi_rnk/_Nr
      rko = modulo(mpi_rnk,_Nr)

        

      do stp = 0, _Ns-1
         src  = modulo(mpi_sze-stp+rnk, _Ns)*_Nr + rko
         mpi_tg = stp + rko
         call mpi_irecv( brecv(1,stp), 2*mes_D%pN*_Ms*i_pZ,  &
            mpi_double_precision, src, mpi_tg, mpi_comm_world,  &
            mpi_rq(stp), mpi_er)
      end do

      do stp = 0, _Ns-1
         dst  = modulo(stp+rnk, _Ns)*_Nr + rko
         jz0  = (dst/_Nr)*i_pZ
         l = 1
         do n = 1, mes_D%pN
            do m = 0, _Ms1
               do j = jz0, jz0+i_pZ-1
                  bsend(l,  stp) =  dble(Taux(j,m,n))
                  bsend(l+1,stp) = dimag(Taux(j,m,n))
                  l = l + 2
               end do
            end do
         end do

         mpi_tg = stp + rko
         call mpi_isend( bsend(1,stp), 2*mes_D%pN*_Ms*i_pZ,  &
            mpi_double_precision, dst, mpi_tg, mpi_comm_world,  &
            mpi_rq(mpi_sze+stp), mpi_er)
      end do
 
      do stp = 0, _Ns-1
         src  = modulo(mpi_sze-stp+rnk, _Ns)*_Nr + rko
         call mpi_wait( mpi_rq(stp), mpi_st, mpi_er)
         pm0 = (src/_Nr)*_Ms
         l = 1
         do n = 1, mes_D%pN
            do m = pm0, pm0+_Ms1
               do j = 0, i_pZ-1
                  Ts(j,m,n) = dcmplx(brecv(l,stp),brecv(l+1,stp))
                  l = l + 2
               end do
            end do
         end do
      end do

      do stp = 0, _Ns-1
         call mpi_wait( mpi_rq(mpi_sze+stp), mpi_st, mpi_er)
      end do

    
   end subroutine tra_T2Ts

!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
   subroutine tra_Ts2T()

      integer :: stp, dst,src, l,j, rnk,rko
      integer :: n,m, pm0,jz0

      rnk = mpi_rnk/_Nr
      rko = modulo(mpi_rnk,_Nr)

      do stp = 0, _Ns-1
         src  = modulo(mpi_sze-stp+rnk, _Ns)*_Nr + rko
         mpi_tg = stp + rko
         call mpi_irecv( brecv(1,stp), 2*i_pN*_Ms*i_pZ,  &
            mpi_double_precision, src, mpi_tg, mpi_comm_world,  &
            mpi_rq(stp), mpi_er)
      end do

      do stp = 0, _Ns-1
         dst  = modulo(stp+rnk, _Ns)*_Nr + rko
         pm0 = (dst/_Nr)*_Ms
         l = 1
         do n = 1, mes_D%pN
            do m = pm0, pm0+_Ms1
               do j = 0, i_pZ-1
                  bsend(l,  stp) =  dble(Ts(j,m,n))
                  bsend(l+1,stp) = dimag(Ts(j,m,n))

                  l = l + 2
               end do
            end do
         end do
         mpi_tg = stp + rko
         call mpi_isend( bsend(1,stp), 2*i_pN*_Ms*i_pZ,  &
            mpi_double_precision, dst, mpi_tg, mpi_comm_world,  &
            mpi_rq(mpi_sze+stp), mpi_er)
      end do

      do stp = 0, _Ns-1
         src  = modulo(mpi_sze-stp+rnk, _Ns)*_Nr + rko
         call mpi_wait( mpi_rq(stp), mpi_st, mpi_er)
         jz0  = (src/_Nr)*i_pZ
         l = 1

         do n = 1, mes_D%pN
            do m = 0, _Ms1
               do j = jz0, jz0+i_pZ-1
                  Taux(j,m,n) = dcmplx(brecv(l,stp),brecv(l+1,stp))
                  l = l + 2
               end do
            end do
         end do
      end do

      do stp = 0, _Ns-1
         call mpi_wait( mpi_rq(mpi_sze+stp), mpi_st, mpi_er)
      end do
            
   end subroutine tra_Ts2T


! #endif

!*************************************************************************
 end module transform
!*************************************************************************

