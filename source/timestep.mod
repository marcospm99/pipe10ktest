  ,&  [   k820309    �
          2021.4.0    "ǈd                                                                                                          
       ./program/timestep.f90 TIMESTEP                                                     
                         @                               ' 7                   #M                 �   �                                 �                        
  p ��������  p 	         & p ��������p �           p 	         p �                                           @                                ' Q                   #IPIV    #M                 �                                    �                           p          p �           p �                                      �                                    �	                       
  p          p          p �           p          p �                                           @                                ' �L                  #RE    #IM 	                �   �                                 ��                       
  p          p �         & p         p d          p �         p e                                     �   �                            	     ��      ^&               
  p          p �         & p         p d          p �         p e                                          @                           
     '�                   #PH0    #PH1    #PH0_    #PH1_                 �                                                               �                                                              �   �                                 P                          p           & p         p O           p P                                      �   �                                 P       H                  p           & p         p O           p P                                           @                                '�G                  #N    #PNI    #PN    #PNI_    #PN_    #R    #INTRDR    #DR0    #DR1    #DR    #RADLAP                 �                                                               �                                                              �                                                              �   �                                 P                          p           & p         p O           p P                                      �   �                                 P       L                  p           & p         p O           p P                                      �   �                                 @      �                
  p ��������  p �         & p ��������p            p �         p                                       �                                    �       �,                
  p          p �           p �                                      �   �                                        �2                
  p          & p        p          & p         p            p          p                                       �   �                                        X3             	   
  p          & p        p          & p         p            p          p                                       �                                            4       7      
      #MESH    p          p            p                                       �                                     7      �            #MESH                                                      
                 
                                 0D0                                                  
                   
                        �                                                                                               d                                                                                                  �               192                                                 �      #HARM 
                                                                                           �                                                           !                                                      512                                            "     �G     #RDOM                                                 #                                                      1                                           $     
                                                    %                                                      4                                             &     
                 
                 R���Q�?        0.51D0                                             '                                
        L            1275070495                                             (                                
         X            1476395009                                             )                                
          D            1140850688                                             *     
                 
                 �h㈵��>        1D-5                                            +                                                         ,                                          ��������                                                     -                                                      6                                             .     
                 
                       �?        0.5D0                                             /     
                 
                 .����B}T        1D99           @                                0     
                  @@                               1     
                  @@                               2     
                  @                                3                       @@                               4                       @                                5     
                  @                                6     
                  @                                7                       @                                8                       @                                9            #         @                                   :                     #         @                                   ;                    #PM <   #C1 =   #C2 >   #A ?             
                                  <                     
                                  =     
                
                                  >     
                D     �                            ?     e       7            p           & p         p d          p e                        #MESH    #         @                                   @                    #PM A   #BC B   #C1 C   #C2 D   #A E             
  @                               A                     
                                  B                     
                                  C     
                
                                  D     
                D @   �                            E     e       Q            p           & p         p d          p e                        #LUMESH    #         @                                   F                    #BC G   #A H   #B I             
                                  G                     
@ @   �                            H     e       Q           p           & p         p d          p e                        #LUMESH              
D @                               I      �L             #COLL    #         @                                   J                    #S K   #A L   #B M   #C N   #D O             
                                  K                     
      �                            L     e       7           p           & p         p d          p e                        #MESH              
                                 M      �L             #COLL              
                                 N      �L             #COLL              
D                                 O      �L             #COLL    #         @                                   P                    #A Q             
D                                 Q      �L             #COLL    #         @                                   R                    #N_ S   #N T             
                                 S      �L             #COLL              
D                                 T      �L             #COLL    #         @                                   U                    #C1 V   #C2 W   #C3 X             
                                 V      �L             #COLL              
                                 W      �L             #COLL              
                                 X      �L             #COLL    #         @                                   Y                     #         @                                   Z                        �   (      fn#fn    �   @   J   VARIABLES      W       MESH+TYPE    _  �   a   MESH%M+TYPE    +  a       LUMESH+TYPE !   �  �   a   LUMESH%IPIV+TYPE    (  �   a   LUMESH%M+TYPE    �  `       COLL+TYPE    D  �   a   COLL%RE+TYPE      �   a   COLL%IM+TYPE    �  v       HARM+TYPE    R  H   a   HARM%PH0+TYPE    �  H   a   HARM%PH1+TYPE    �  �   a   HARM%PH0_+TYPE    �  �   a   HARM%PH1_+TYPE    :  �       RDOM+TYPE    �  H   a   RDOM%N+TYPE    6	  H   a   RDOM%PNI+TYPE    ~	  H   a   RDOM%PN+TYPE    �	  �   a   RDOM%PNI_+TYPE    r
  �   a   RDOM%PN_+TYPE      �   a   RDOM%R+TYPE !   �  �   a   RDOM%INTRDR+TYPE    �  �   a   RDOM%DR0+TYPE    b  �   a   RDOM%DR1+TYPE    >  �   a   RDOM%DR+TYPE !   �  Z   a   RDOM%RADLAP+TYPE "   >  s       D_TIME+PARAMETERS &   �  p       D_TIMESTEP+PARAMETERS !   !  p       I_PH1+PARAMETERS    �  s       I_N+PARAMETERS       J       VAR_H+VARIABLES     N  p       I_K1+PARAMETERS    �  s       I_K+PARAMETERS    1  J       MES_D+TYPE     {  q       I_MP+PARAMETERS #   �  @       D_ALPHA+PARAMETERS     ,  q       I_KL+PARAMETERS &   �  v       D_IMPLICIT+PARAMETERS *     z       MPI_DOUBLE_PRECISION+MPIF    �  z       MPI_MAX+MPIF $     z       MPI_COMM_WORLD+MPIF #   �  t       D_DTERR+PARAMETERS    �  @       MPI_RNK+MPIF &   5  p       I_MAXTSTEP+PARAMETERS (   �  q       I_SAVE_RATE2+PARAMETERS %     u       D_COURANT+PARAMETERS #   �  t       D_MAXDT+PARAMETERS    �  @       TIM_T    ?  @       TIM_DT      @       TIM_DTERR    �  @       TIM_IT    �  @       TIM_STEP    ?  @       TIM_CORR_DT      @       TIM_CFL_DT    �  @       TIM_CFL_DIR    �  @       TIM_NEW_DT    ?  @       TIM_PCN      H       TIM_PRECOMPUTE    �  g       TIM_MESH_INIT !   .  @   a   TIM_MESH_INIT%PM !   n  @   a   TIM_MESH_INIT%C1 !   �  @   a   TIM_MESH_INIT%C2     �  �   a   TIM_MESH_INIT%A     �  o       TIM_LUMESH_INIT #     @   a   TIM_LUMESH_INIT%PM #   K  @   a   TIM_LUMESH_INIT%BC #   �  @   a   TIM_LUMESH_INIT%C1 #   �  @   a   TIM_LUMESH_INIT%C2 "     �   a   TIM_LUMESH_INIT%A "   �  ^       TIM_LUMESH_INVERT %     @   a   TIM_LUMESH_INVERT%BC $   Y  �   a   TIM_LUMESH_INVERT%A $   	   R   a   TIM_LUMESH_INVERT%B    [   k       TIM_MESHMULT    �   @   a   TIM_MESHMULT%S    !  �   a   TIM_MESHMULT%A    �!  R   a   TIM_MESHMULT%B    "  R   a   TIM_MESHMULT%C    X"  R   a   TIM_MESHMULT%D    �"  O       TIM_ZEROBC    �"  R   a   TIM_ZEROBC%A    K#  W       TIM_NLINCORR     �#  R   a   TIM_NLINCORR%N_    �#  R   a   TIM_NLINCORR%N     F$  `       TIM_MEASURECORR #   �$  R   a   TIM_MEASURECORR%C1 #   �$  R   a   TIM_MEASURECORR%C2 #   J%  R   a   TIM_MEASURECORR%C3    �%  H       TIM_CHECK_CGCE    �%  H       TIM_NEW_TSTEP 