!*************************************************************************
#include "../parallel.h"
 module wksp
!*************************************************************************
   use type

   implicit none
   save

   type (phys) :: vel_r
   type (phys) :: vel_t
   type (phys) :: vel_z
   type (coll) :: vel_ur
   type (coll) :: vel_ut
   type (coll) :: vel_uz
   type (coll) :: vel_Nr
   type (coll) :: vel_Nt
   type (coll) :: vel_Nz


   type (lumesh) :: LDp(0:i_pH1), LDm(0:i_pH1)
   type (lumesh) :: LDz(0:i_pH1), LNp(0:i_pH1)
   type (mesh)   :: Ltp(0:i_pH1), Ltm(0:i_pH1)
   type (mesh)   :: Ltz(0:i_pH1)
   type (coll)   :: Nr_,Nt_,Nz_,ur_,ut_,uz_

   type (coll) :: c1,c2,c3,c4
   type (phys) :: p1,p2,p3,p4

   double complex :: Ts(0:i_pZ-1, 0:i_M1, i_pN)
   type (spec)    :: s1


   ! FFT y variables (allocatear el valor máximo de los tamaños)

   ! double precision :: bsend(2*i_pN*_Ms*i_pZ,0:_Ns-1)
   ! double precision :: brecv(2*i_pN*_Ms*i_pZ,0:_Ns-1)

   double precision :: bsend(2*(i_pH1+1)*i_pN*3,0:_Nr-1)
   double precision :: brecv(2*(i_pH1+1)*i_pN*3,0:_Nr-1)
   
   ! Hay que deallocatear T1 también!!!!!!!!!!!

end module