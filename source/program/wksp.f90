!*************************************************************************
#include "../parallel.h"
 module wksp
!*************************************************************************
   use type
   use mpif


      ! DISCLAIMER: hay un pto debil, al ser types allocatables, es necesario ponerlos como inout
   implicit none
   save

   ! Phys
   type (phys)    :: vel_r
   type (phys)    :: vel_t
   type (phys)    :: vel_z
   type (phys)    :: p1,p2,p3,p4

   ! Coll
   type (coll) :: vel_ur,vel_ut,vel_uz
   type (coll) :: vel_Nr,vel_Nt,vel_Nz
   type (coll) :: Nr_,Nt_,Nz_
   type (coll) :: ur_,ut_,uz_
   type (coll) :: c1,c2,c3,c4
   type (coll) :: cp

   ! Spec
   type (spec)    :: s1

   ! Lumesh
   type (lumesh) :: LDp(0:i_pH1), LDm(0:i_pH1)
   type (lumesh) :: LDz(0:i_pH1), LNp(0:i_pH1)
   
   ! Mesh
   type (mesh)   :: Ltp(0:i_pH1), Ltm(0:i_pH1), Ltz(0:i_pH1)


   ! FFT y variables (allocatear el valor máximo de los tamaños)

   double precision, allocatable, dimension(:,:) :: bsend
   double precision, allocatable, dimension(:,:) :: brecv

   double complex, allocatable, dimension(:,:,:) :: Ts
   double complex, allocatable, dimension(:,:,:) :: T1
   
   
   ! Hay que deallocatear T1 también!!!!!!!!!!!


   contains

   subroutine initmem()
   implicit none

   ! miscelaneous
   if(.not.allocated(Ts)) allocate(Ts(0:i_pZ-1, 0:i_M1, i_pN))
   if(.not.allocated(T1))  allocate(T1(0:i_3K-1, 0:_Ms-1, i_pN))

   if(.not.allocated(bsend)) allocate(bsend(max(2*i_pN*(i_pH1+1)*3,2*i_pN*_Ms*i_pZ),0:max(_Nr-1,_Ns-1)))
   if(.not.allocated(brecv)) allocate(brecv(max(2*i_pN*(i_pH1+1)*3,2*i_pN*_Ms*i_pZ),0:max(_Nr-1,_Ns-1)))

   ! spec
   if(.not.allocated(s1%Re)) allocate(s1%Re(0:_Hs1, i_pN))
   if(.not.allocated(s1%Im)) allocate(s1%Im(0:_Hs1, i_pN))

   ! phys

   if(.not.allocated(p1%Re)) allocate(p1%Re(0:i_pZ-1, 0:i_Th-1, i_pN))
   if(.not.allocated(p2%Re)) allocate(p2%Re(0:i_pZ-1, 0:i_Th-1, i_pN))
   if(.not.allocated(p3%Re)) allocate(p3%Re(0:i_pZ-1, 0:i_Th-1, i_pN))
   if(.not.allocated(p4%Re)) allocate(p4%Re(0:i_pZ-1, 0:i_Th-1, i_pN))

   if(.not.allocated(vel_r%Re)) allocate(vel_r%Re(0:i_pZ-1, 0:i_Th-1, i_pN))
   if(.not.allocated(vel_t%Re)) allocate(vel_t%Re(0:i_pZ-1, 0:i_Th-1, i_pN))
   if(.not.allocated(vel_z%Re)) allocate(vel_z%Re(0:i_pZ-1, 0:i_Th-1, i_pN))


   ! coll
   if (.not.allocated(cp%Re)) allocate(cp%Re(i_N, 0:i_pH1))
   if (.not.allocated(cp%Im)) allocate(cp%Im(i_N, 0:i_pH1))

   if (.not.allocated(c1%Re)) allocate(c1%Re(i_N, 0:i_pH1))
   if (.not.allocated(c1%Im)) allocate(c1%Im(i_N, 0:i_pH1))
   if (.not.allocated(c2%Re)) allocate(c2%Re(i_N, 0:i_pH1))
   if (.not.allocated(c2%Im)) allocate(c2%Im(i_N, 0:i_pH1))
   if (.not.allocated(c3%Re)) allocate(c3%Re(i_N, 0:i_pH1))
   if (.not.allocated(c3%Im)) allocate(c3%Im(i_N, 0:i_pH1))
   if (.not.allocated(c4%Re)) allocate(c4%Re(i_N, 0:i_pH1))
   if (.not.allocated(c4%Im)) allocate(c4%Im(i_N, 0:i_pH1))

   if (.not.allocated(vel_ur%Re)) allocate(vel_ur%Re(i_N, 0:i_pH1))
   if (.not.allocated(vel_ur%Im)) allocate(vel_ur%Im(i_N, 0:i_pH1))
   if (.not.allocated(vel_ut%Re)) allocate(vel_ut%Re(i_N, 0:i_pH1))
   if (.not.allocated(vel_ut%Im)) allocate(vel_ut%Im(i_N, 0:i_pH1))
   if (.not.allocated(vel_uz%Re)) allocate(vel_uz%Re(i_N, 0:i_pH1))
   if (.not.allocated(vel_uz%Im)) allocate(vel_uz%Im(i_N, 0:i_pH1))

   if (.not.allocated(ur_%Re)) allocate(ur_%Re(i_N, 0:i_pH1))
   if (.not.allocated(ur_%Im)) allocate(ur_%Im(i_N, 0:i_pH1))
   if (.not.allocated(ut_%Re)) allocate(ut_%Re(i_N, 0:i_pH1))
   if (.not.allocated(ut_%Im)) allocate(ut_%Im(i_N, 0:i_pH1))
   if (.not.allocated(uz_%Re)) allocate(uz_%Re(i_N, 0:i_pH1))
   if (.not.allocated(uz_%Im)) allocate(uz_%Im(i_N, 0:i_pH1))

   if (.not.allocated(Nr_%Re)) allocate(Nr_%Re(i_N, 0:i_pH1))
   if (.not.allocated(Nr_%Im)) allocate(Nr_%Im(i_N, 0:i_pH1))
   if (.not.allocated(Nt_%Re)) allocate(Nt_%Re(i_N, 0:i_pH1))
   if (.not.allocated(Nt_%Im)) allocate(Nt_%Im(i_N, 0:i_pH1))
   if (.not.allocated(Nz_%Re)) allocate(Nz_%Re(i_N, 0:i_pH1))
   if (.not.allocated(Nz_%Im)) allocate(Nz_%Im(i_N, 0:i_pH1))

   if (.not.allocated(vel_Nr%Re)) allocate(vel_Nr%Re(i_N, 0:i_pH1))
   if (.not.allocated(vel_Nr%Im)) allocate(vel_Nr%Im(i_N, 0:i_pH1))
   if (.not.allocated(vel_Nt%Re)) allocate(vel_Nt%Re(i_N, 0:i_pH1))
   if (.not.allocated(vel_Nt%Im)) allocate(vel_Nt%Im(i_N, 0:i_pH1))
   if (.not.allocated(vel_Nz%Re)) allocate(vel_Nz%Re(i_N, 0:i_pH1))
   if (.not.allocated(vel_Nz%Im)) allocate(vel_Nz%Im(i_N, 0:i_pH1))

   end subroutine

   subroutine deallocmem()
   implicit none

   ! miscelaneous
   deallocate(T1)
   deallocate(Ts)

   deallocate(bsend)
   deallocate(brecv)

   ! spec
   deallocate(s1%Re)
   deallocate(s1%Im)

   ! phys
   deallocate(p1%Re)
   deallocate(p2%Re)
   deallocate(p3%Re)
   deallocate(p4%Re)
   deallocate(vel_r%Re)
   deallocate(vel_t%Re)
   deallocate(vel_z%Re)

   ! coll
   deallocate(cp%Re)
   deallocate(cp%Im)

   deallocate(c1%Re)
   deallocate(c1%Im)
   deallocate(c2%Re)
   deallocate(c2%Im)
   deallocate(c3%Re)
   deallocate(c3%Im)
   deallocate(c4%Re)
   deallocate(c4%Im)

   deallocate(vel_ur%Re)
   deallocate(vel_ur%Im)
   deallocate(vel_ut%Re)
   deallocate(vel_ut%Im)
   deallocate(vel_uz%Re)
   deallocate(vel_uz%Im)

   deallocate(ur_%Re)
   deallocate(ur_%Im)
   deallocate(ut_%Re)
   deallocate(ut_%Im)
   deallocate(uz_%Re)
   deallocate(uz_%Im)

   deallocate(Nr_%Re)
   deallocate(Nr_%Im)
   deallocate(Nt_%Re)
   deallocate(Nt_%Im)
   deallocate(Nz_%Re)
   deallocate(Nz_%Im)

   deallocate(vel_Nr%Re)
   deallocate(vel_Nr%Im)
   deallocate(vel_Nt%Re)
   deallocate(vel_Nt%Im)
   deallocate(vel_Nz%Re)
   deallocate(vel_Nz%Im)

   end subroutine

end module