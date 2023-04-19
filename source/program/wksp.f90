#include "../parallel.h"
 module wksp

use type
!use variables

implicit none
save



   type (coll) :: c1,c2,c3,c4          ! Three colls are defined here. Why! They are really big. 
                                     ! They are defined as private. They cannot be used anywhere else
                                     ! Remove in future versions. We have to pass the routine some workarray 
   type (phys) :: p1,p2,p3,p4 
  
  type (coll) :: r,t,z ! SOLO PARA vel_corrector y no sé qué pasa

  type (spec) :: s1,s2,s3


 end module wksp