!*************************************************************************
#include "../parallel.h"
 module type
!*************************************************************************

   use parameters
!    use meshs
   implicit none
   save

   type harm
      integer              :: pH0,pH1, pH0_(0:_Np-1),pH1_(0:_Np-1)
   end type harm

   ! type spec
   !    double precision     :: Re(0:_Hs1, i_pN)
   !    double precision     :: Im(0:_Hs1, i_pN)
   ! end type spec

   type spec
      double precision, allocatable     :: Re(:,:)
      double precision, allocatable     :: Im(:,:)
   end type spec

!   integer,          parameter :: i_pH1 = (_Nr+_Hs1)/_Nr-1

   type coll
      double precision     :: Re(i_N, 0:i_pH1)
      double precision     :: Im(i_N, 0:i_pH1)
   end type coll

   ! type coll
   !    double precision, allocatable     :: Re(:,:)
   !    double precision, allocatable     :: Im(:,:)
   ! end type coll


   ! type phys
   !    double precision     :: Re(0:i_pZ-1, 0:i_Th-1, i_pN)
   ! end type phys


   type phys
      double precision, allocatable, dimension(:,:,:)     :: Re
   end type phys

                             ! M(KL+1+n-j, j) = A(n,j)
   type mesh
      double precision :: M(2*i_KL+1, 1-i_KL:i_N)
   end type mesh
                                ! M(2*KL+1+n-j, j) = A(n,j)
                                ! see lapack dgbtrf
   type lumesh
      integer          :: ipiv(i_N)
      double precision :: M(3*i_KL+1, i_N)
   end type lumesh


   type rdom
      integer          :: N, pNi,pN, pNi_(0:_Np-1),pN_(0:_Np-1)
      double precision :: r(i_N,-3:3)
      double precision :: intrdr(i_N)
      double precision :: dr0(1:1+i_KL,0:i_KL), dr1(1:1+i_KL,0:i_KL)
      type (mesh)      :: dr(i_KL)
      type (mesh)      :: radLap
   end type rdom


   type (rdom) :: mes_D

   end module