#include "../parallel.h"
 module wksp

use type

implicit none
save


   
    type (coll) :: c1,c2,c3,c4          ! Three colls are defined here. Why! They are really big. 
                                     ! They are defined as private. They cannot be used anywhere else
                                     ! Remove in future versions. We have to pass the routine some workarray 
   type (phys) :: p1,p2,p3,p4 
  
! Velocity

   double precision :: vel_nu
   double precision :: vel_Pr0
   double precision :: vel_U(i_N)
   double precision :: vel_Up(i_N)

   type (lumesh)     :: LDp(0:i_pH1), LDm(0:i_pH1)
   type (lumesh)     :: LDz(0:i_pH1), LNp(0:i_pH1)
   type (mesh)       :: Ltp(0:i_pH1), Ltm(0:i_pH1)
   type (mesh)       :: Ltz(0:i_pH1)
   type (coll)       :: Nr_,Nt_,Nz_,ur_,ut_,uz_ ! Mem_tot = 6*2*8*384*1819*72/1024**3


  

  !! Esenciales, no se pueden tocar
  type (coll) :: r,t,z ! SOLO PARA vel_corrector y no sé qué pasa

   type (phys) :: vel_r
   type (phys) :: vel_t
   type (phys) :: vel_z
   type (coll) :: vel_ur
   type (coll) :: vel_ut
   type (coll) :: vel_uz
   type (coll) :: vel_Nr
   type (coll) :: vel_Nt
   type (coll) :: vel_Nz

   !  double precision :: bsend(2*i_pN*_Ms*i_pZ,0:_Ns-1) ! -> 1.16 GB
   !  double precision :: brecv(2*i_pN*_Ms*i_pZ,0:_Ns-1) ! -> 1.16 GB
    type (spec) :: s1,s2,s3 ! -> 2.3 GB each one

    double precision :: wk(i_N, 0:i_pH1)
 end module wksp