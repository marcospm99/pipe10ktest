!*************************************************************************
#include "../parallel.h"
 module type
!*************************************************************************
   use parameters

   implicit none
   save

   type harm
      integer              :: pH0,pH1, pH0_(0:_Np-1),pH1_(0:_Np-1)
   end type harm

   type spec
      double precision     :: Re(0:_Hs1, i_pN)
      double precision     :: Im(0:_Hs1, i_pN)
   end type spec

!   integer,          parameter :: i_pH1 = (_Nr+_Hs1)/_Nr-1
   type coll
      double precision     :: Re(i_N, 0:i_pH1)
      double precision     :: Im(i_N, 0:i_pH1)
   end type coll

   type phys
      double precision     :: Re(0:i_pZ-1, 0:i_Th-1, i_pN)
   end type phys


end module type