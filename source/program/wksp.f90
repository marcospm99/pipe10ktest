#include "../parallel.h"
 module wksp
use mpif
use type

implicit none
! save
   

   double complex, allocatable, dimension(:,:,:) :: Ts
   ! double complex :: Ts(0:i_pZ-1, 0:i_M1, i_pN) 
   

! Velocity

   double precision :: vel_nu
   double precision :: vel_Pr0
   double precision :: vel_U(i_N)
   double precision :: vel_Up(i_N)


   type (lumesh)         :: LDp(0:i_pH1), LDm(0:i_pH1)
   type (lumesh)         :: LDz(0:i_pH1), LNp(0:i_pH1)


   type (mesh)           :: Ltp(0:i_pH1), Ltm(0:i_pH1)
   type (mesh)           :: Ltz(0:i_pH1)


   type (coll)           :: c1,c2,c3,c4
   type (coll)           :: Nr_,Nt_,Nz_
   type (coll)           :: ur_,ut_,uz_ 
   type (coll)           :: vel_Nr, vel_Nt, vel_Nz 
   type (coll)           :: vel_ur, vel_ut, vel_uz  ! Imprescindibles



  !! Esenciales, no se pueden tocar
  !type (coll) :: r,t,z ! SOLO PARA vel_corrector y no sé qué pasa

   type (phys) :: vel_r, vel_t,vel_z 
   type (phys) :: p1, p2, p3, p4

   !  double precision :: bsend(2*i_pN*_Ms*i_pZ,0:_Ns-1) ! -> 1.16 GB
   !  double precision :: brecv(2*i_pN*_Ms*i_pZ,0:_Ns-1) ! -> 1.16 GB

      ! Usar estas, que son más grandes

       double precision, allocatable, dimension(:,:) :: bsend
       double precision, allocatable, dimension(:,:) :: brecv

      ! double precision :: bsend(2*i_pN*(i_pH1+1)*3,0:_Nr-1)
      ! double precision :: brecv(2*i_pN*(i_pH1+1)*3,0:_Nr-1)
    

       type (spec)      :: s1 !-> 2.3 GB each one



   ! Memory chunks, very little BUT

   !  double precision :: inRe(1-i_KL:i_N), inIm(1-i_KL:i_N)
   !  double precision :: inRe(i_N), inIm(i_N)
    double precision :: d(i_N)
    double precision, save :: dt,dz, r_(i_N)
    double precision :: re(i_N), im(i_N)


      contains


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine allospec()

   implicit none
   integer :: dimension, ii

   ! chunks

      allocate(Ts(0:i_pZ-1, 0:i_M1, i_pN)) 
      allocate(bsend(2*i_pN*(i_pH1+1)*3,0:_Nr-1))
      allocate(brecv(2*i_pN*(i_pH1+1)*3,0:_Nr-1))


   ! spec
      if(.not.allocated(s1%Re))   allocate(s1%Re(0:_Hs1, i_pN))
      if(.not.allocated(s1%Im))   allocate(s1%Im(0:_Hs1, i_pN))


   ! ! phys
      if(.not.allocated(p1%Re))      allocate(p1%Re   (0:i_pZ-1, 0:i_Th-1, i_pN))
      if(.not.allocated(p2%Re))      allocate(p2%Re   (0:i_pZ-1, 0:i_Th-1, i_pN))
      if(.not.allocated(p3%Re))      allocate(p3%Re   (0:i_pZ-1, 0:i_Th-1, i_pN))
      if(.not.allocated(p4%Re))      allocate(p4%Re   (0:i_pZ-1, 0:i_Th-1, i_pN))
      if(.not.allocated(vel_r%Re))   allocate(vel_r%Re(0:i_pZ-1, 0:i_Th-1, i_pN))
      if(.not.allocated(vel_t%Re))   allocate(vel_t%Re(0:i_pZ-1, 0:i_Th-1, i_pN))
      if(.not.allocated(vel_z%Re))   allocate(vel_z%Re(0:i_pZ-1, 0:i_Th-1, i_pN))

   ! colls 

! if(.not.allocated(c1%Re))       allocate(c1%Re(i_N, 0:i_pH1))
! if(.not.allocated(c1%Im))       allocate(c1%Im(i_N, 0:i_pH1))
! if(.not.allocated(c2%Re))       allocate(c2%Re(i_N, 0:i_pH1))
! if(.not.allocated(c2%Im))       allocate(c2%Im(i_N, 0:i_pH1))
! if(.not.allocated(c3%Re))       allocate(c3%Re(i_N, 0:i_pH1))
! if(.not.allocated(c3%Im))       allocate(c3%Im(i_N, 0:i_pH1))
! if(.not.allocated(c4%Re))       allocate(c4%Re(i_N, 0:i_pH1))
! if(.not.allocated(c4%Im))       allocate(c4%Im(i_N, 0:i_pH1))

! if(.not.allocated(Nr_%Re))      allocate(Nr_%Re(i_N, 0:i_pH1))
! if(.not.allocated(Nr_%Im))      allocate(Nr_%Im(i_N, 0:i_pH1))
! if(.not.allocated(Nt_%Re))      allocate(Nt_%Re(i_N, 0:i_pH1))
! if(.not.allocated(Nt_%Im))      allocate(Nt_%Im(i_N, 0:i_pH1))
! if(.not.allocated(Nz_%Re))      allocate(Nz_%Re(i_N, 0:i_pH1))
! if(.not.allocated(Nz_%Im))      allocate(Nz_%Im(i_N, 0:i_pH1))
! if(.not.allocated(ur_%Re))      allocate(ur_%Re(i_N, 0:i_pH1))
! if(.not.allocated(ur_%Im))      allocate(ur_%Im(i_N, 0:i_pH1))
! if(.not.allocated(ut_%Re))      allocate(ut_%Re(i_N, 0:i_pH1))
! if(.not.allocated(ut_%Im))      allocate(ut_%Im(i_N, 0:i_pH1))
! if(.not.allocated(uz_%Re))      allocate(uz_%Re(i_N, 0:i_pH1))
! if(.not.allocated(uz_%Im))      allocate(uz_%Im(i_N, 0:i_pH1))
! if(.not.allocated(vel_Nr%Re))   allocate(vel_Nr%Re(i_N, 0:i_pH1))
! if(.not.allocated(vel_Nr%Im))   allocate(vel_Nr%Im(i_N, 0:i_pH1))
! if(.not.allocated(vel_Nt%Re))   allocate(vel_Nt%Re(i_N, 0:i_pH1))
! if(.not.allocated(vel_Nt%Im))   allocate(vel_Nt%Im(i_N, 0:i_pH1))
! if(.not.allocated(vel_Nz%Re))   allocate(vel_Nz%Re(i_N, 0:i_pH1))
! if(.not.allocated(vel_Nz%Im))   allocate(vel_Nz%Im(i_N, 0:i_pH1))
! if(.not.allocated(vel_ur%Re))   allocate(vel_ur%Re(i_N, 0:i_pH1))
! if(.not.allocated(vel_ur%Im))   allocate(vel_ur%Im(i_N, 0:i_pH1))
! if(.not.allocated(vel_ut%Re))   allocate(vel_ut%Re(i_N, 0:i_pH1))
! if(.not.allocated(vel_ut%Im))   allocate(vel_ut%Im(i_N, 0:i_pH1))
! if(.not.allocated(vel_uz%Re))   allocate(vel_uz%Re(i_N, 0:i_pH1))
! if(.not.allocated(vel_uz%Im))   allocate(vel_uz%Im(i_N, 0:i_pH1))

   !  c1%Re(:,:)=0d0
   !  c1%Im(:,:)=0d0
   !  c2%Re(:,:)=0d0
   !  c2%Im(:,:)=0d0
   !  c3%Re(:,:)=0d0
   !  c3%Im(:,:)=0d0
   !  c4%Re(:,:)=0d0
   !  c4%Im(:,:)=0d0
   !  Nr_%Re(:,:)=0d0
   !  Nr_%Im(:,:)=0d0
   !  Nt_%Re(:,:)=0d0
   !  Nt_%Im(:,:)=0d0
   !  Nz_%Re(:,:)=0d0
   !  Nz_%Im(:,:)=0d0
   !  ur_%Re(:,:)=0d0
   !  ur_%Im(:,:)=0d0
   !  ut_%Re(:,:)=0d0
   !  ut_%Im(:,:)=0d0
   !  uz_%Re(:,:)=0d0
   !  uz_%Im(:,:)=0d0
   !  vel_Nr%Re(:,:)=0d0
   !  vel_Nr%Im(:,:)=0d0
   !  vel_Nt%Re(:,:)=0d0
   !  vel_Nt%Im(:,:)=0d0
   !  vel_Nz%Re(:,:)=0d0
   !  vel_Nz%Im(:,:)=0d0
   !  vel_ur%Re(:,:)=0d0
   !  vel_ur%Im(:,:)=0d0
   !  vel_ut%Re(:,:)=0d0
   !  vel_ut%Im(:,:)=0d0
   !  vel_uz%Re(:,:)=0d0
   !  vel_uz%Im(:,:)=0d0

   ! write(*,*) 'tamaño de los coll', size(c1%re,1), 'x', size(c1%re,2)




















!!!!!!!!!!!!!

! lumesh y mesh

dimension = 3*i_KL+1

do ii =0,i_pH1
if(.not.allocated(LDp(ii)%M)) allocate(LDp(ii)%M(dimension, i_N),LDp(ii)%ipiv(i_N))
if(.not.allocated(LDm(ii)%M)) allocate(LDm(ii)%M(dimension, i_N),LDm(ii)%ipiv(i_N))
if(.not.allocated(LDz(ii)%M)) allocate(LDz(ii)%M(dimension, i_N),LDz(ii)%ipiv(i_N))
if(.not.allocated(LNp(ii)%M)) allocate(LNp(ii)%M(dimension, i_N),LNp(ii)%ipiv(i_N))

! mesh
! if(.not.allocated(Ltp(ii)%M)) allocate(Ltp(ii)%M(2*i_KL+1, 1-i_KL:i_N))
! if(.not.allocated(Ltm(ii)%M)) allocate(Ltm(ii)%M(2*i_KL+1, 1-i_KL:i_N))
! if(.not.allocated(Ltz(ii)%M)) allocate(Ltz(ii)%M(2*i_KL+1, 1-i_KL:i_N))

enddo
   

   
   end subroutine

   subroutine deallospec()
   implicit none
   integer :: ii

   ! memory chunks
      deallocate(Ts) 
      deallocate(bsend)
      deallocate(brecv) 
   ! spec
   deallocate(s1%Re,s1%Im)

   ! colls
   ! deallocate(c1%Re,c1%Im,c2%Re,c2%Im,c3%Re,c3%Im,c4%Re,c4%Im)
   ! deallocate(Nr_%Re,Nr_%Im,Nt_%Re,Nt_%Im,Nz_%Re,Nz_%Im)
   ! deallocate(ur_%Re,ur_%Im,ut_%Re,ut_%Im,uz_%Re,uz_%Im)
   ! deallocate(vel_Nr%Re,vel_Nr%Im,vel_Nt%Re,vel_Nt%Im,vel_Nz%Re,vel_Nz%Im)
   ! deallocate(vel_ur%Re,vel_ur%Im,vel_ut%Re,vel_ut%Im,vel_uz%Re,vel_uz%Im)


   ! ! phys

   deallocate(p1%Re,p2%Re,p3%Re,p4%Re,vel_r%Re,vel_t%Re,vel_z%Re)

   !lumesh y mesh

   do ii =0,i_pH1

   ! lumesh
   deallocate(LDp(ii)%M, LDp(ii)%ipiv)
   deallocate(LDm(ii)%M, LDm(ii)%ipiv)
   deallocate(LDz(ii)%M, LDz(ii)%ipiv)
   deallocate(LNp(ii)%M, LNp(ii)%ipiv)

   ! ! mesh
   ! deallocate(Ltp(ii)%M)
   ! deallocate(Ltm(ii)%M)
   ! deallocate(Ltz(ii)%M)
   enddo

   end subroutine


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine setzero()
implicit none
    s1%Re(:,:)=0d0
    s1%Im(:,:)=0d0
end subroutine
   
 end module wksp