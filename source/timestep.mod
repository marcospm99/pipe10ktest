  .&  [   k820309    ļ
          2021.4.0    ú	ôc                                                                                                          
       ./program/timestep.f90 TIMESTEP                      @                              
                         @                               ' m                   #M                                                     ¤                        
  p æ˙˙˙˙˙˙˙  p 	         & p ũ˙˙˙˙˙˙˙p           p 	         p                                           @                                ' ĸ                   #IPIV    #M                                                                               p          p           p                                                                                                 
  p          p          p           p          p                                           @                                ' Ē                  #RE    #IM 	                                                    ¨
                       
  p          p         & p         p           p         p                                                                     	     ¨
      DU               
  p          p         & p         p           p         p                                           @                           
     '                  #N    #PNI    #PN    #PNI_    #PN_    #R    #INTRDR    #DR0    #DR1    #DR    #RADLAP                                                                                                                                                                                                                                                H                          p           & p         p G           p H                                                                          H       ,                  p           & p         p G           p H                                                                          
      P                
  p û˙˙˙˙˙˙  p         & p ũ˙˙˙˙˙˙˙p            p         p                                                                                 PV                
  p          p           p                                                                                 Pb                
  p          & p        p          & p         p            p          p                                                                                  c             	   
  p          & p        p          & p         p            p          p                                                                                  āc       m      
      #MESH    p          p            p                                                                            m      `            #MESH                      @                                'H                   #PH0    #PH1    #PH0_    #PH1_                                                                                                                                                                                  H                          p           & p         p G           p H                                                                          H       (                  p           & p         p G           p H                                                                           
                 
                                 0D0                                                  
                   
                        đŋ                                                                                                                                                                                                               384                                                 H      #HARM                                                                                            ˙                                                           !                                                      512                                            "          #RDOM 
                                                #                                                      1                                           $     
                                                    %                                                      4                                             &     
                 
                 R¸ëQā?        0.51D0                                             '                                
        L            1275070495                                             (                                
         X            1476395009                                             )                                
          D            1140850688                                             *     
                 
                 ņhãĩøä>        1D-5                                            +                                                         ,                                          ˙˙˙˙˙˙˙˙                                                     -                                       ô              500                                             .     
                 
                       ā?        0.5D0                                             /     
                 
                 .ĸŽB}T        1D99           @                                0     
                  @@                               1     
                  @@                               2     
                  @                                3                       @@                               4                       @                                5     
                  @                                6     
                  @                                7                       @                                8                       @                                9            #         @                                   :                     #         @                                   ;                    #PM <   #C1 =   #C2 >   #A ?             
                                  <                     
                                  =     
                
                                  >     
                D                                 ?            m            p           & p         p           p                         #MESH    #         @                                   @                    #PM A   #BC B   #C1 C   #C2 D   #A E             
  @                               A                     
                                  B                     
                                  C     
                
                                  D     
                D @                               E            ĸ            p           & p         p           p                         #LUMESH    #         @                                   F                    #BC G   #A H   #B I             
                                  G                     
@ @                               H            ĸ           p           & p         p           p                         #LUMESH              
D @                               I      Ē             #COLL    #         @                                   J                    #S K   #A L   #B M   #C N   #D O             
                                  K                     
                                  L            m           p           & p         p           p                         #MESH              
                                  M      Ē            #COLL              
                                  N      Ē            #COLL              D                                 O      Ē             #COLL    #         @                                   P                    #A Q             D                                 Q      Ē             #COLL    #         @                                   R                    #N_ S   #N T             
                                  S      Ē            #COLL              
D                                 T      Ē             #COLL    #         @                                   U                    #C1 V   #C2 W   #C3 X             
                                  V      Ē            #COLL              
                                  W      Ē            #COLL              
                                  X      Ē            #COLL    #         @                                   Y                     #         @                                   Z                            (      fn#fn    Č   @   J   VARIABLES      W       MESH+MESHS    _  Ė   a   MESH%M+MESHS    +  a       LUMESH+MESHS "        a   LUMESH%IPIV+MESHS    (  ŧ   a   LUMESH%M+MESHS    ä  `       COLL+VARIABLES "   D  Ė   a   COLL%RE+VARIABLES "     Ė   a   COLL%IM+VARIABLES    Ü  ´       RDOM+MESHS      H   a   RDOM%N+MESHS    Ø  H   a   RDOM%PNI+MESHS       H   a   RDOM%PN+MESHS     h  Ŧ   a   RDOM%PNI_+MESHS      Ŧ   a   RDOM%PN_+MESHS    Ā  Ė   a   RDOM%R+MESHS "   	     a   RDOM%INTRDR+MESHS    (
  Ü   a   RDOM%DR0+MESHS      Ü   a   RDOM%DR1+MESHS    ā  Ļ   a   RDOM%DR+MESHS "     Z   a   RDOM%RADLAP+MESHS    ā  v       HARM+VARIABLES #   V  H   a   HARM%PH0+VARIABLES #     H   a   HARM%PH1+VARIABLES $   æ  Ŧ   a   HARM%PH0_+VARIABLES $     Ŧ   a   HARM%PH1_+VARIABLES "   >  s       D_TIME+PARAMETERS &   ą  p       D_TIMESTEP+PARAMETERS !   !  p       I_PH1+PARAMETERS      s       I_N+PARAMETERS       J       VAR_H+VARIABLES     N  p       I_K1+PARAMETERS    ž  s       I_K+PARAMETERS    1  J       MES_D+MESHS     {  q       I_MP+PARAMETERS #   ė  @       D_ALPHA+PARAMETERS     ,  q       I_KL+PARAMETERS &     v       D_IMPLICIT+PARAMETERS *     z       MPI_DOUBLE_PRECISION+MPIF      z       MPI_MAX+MPIF $     z       MPI_COMM_WORLD+MPIF #     t       D_DTERR+PARAMETERS    õ  @       MPI_RNK+MPIF &   5  p       I_MAXTSTEP+PARAMETERS (   Ĩ  s       I_SAVE_RATE2+PARAMETERS %     u       D_COURANT+PARAMETERS #     t       D_MAXDT+PARAMETERS      @       TIM_T    A  @       TIM_DT      @       TIM_DTERR    Á  @       TIM_IT      @       TIM_STEP    A  @       TIM_CORR_DT      @       TIM_CFL_DT    Á  @       TIM_CFL_DIR      @       TIM_NEW_DT    A  @       TIM_PCN      H       TIM_PRECOMPUTE    É  g       TIM_MESH_INIT !   0  @   a   TIM_MESH_INIT%PM !   p  @   a   TIM_MESH_INIT%C1 !   °  @   a   TIM_MESH_INIT%C2     đ  Ž   a   TIM_MESH_INIT%A       o       TIM_LUMESH_INIT #     @   a   TIM_LUMESH_INIT%PM #   M  @   a   TIM_LUMESH_INIT%BC #     @   a   TIM_LUMESH_INIT%C1 #   Í  @   a   TIM_LUMESH_INIT%C2 "     °   a   TIM_LUMESH_INIT%A "   Ŋ  ^       TIM_LUMESH_INVERT %     @   a   TIM_LUMESH_INVERT%BC $   [  °   a   TIM_LUMESH_INVERT%A $      R   a   TIM_LUMESH_INVERT%B    ]   k       TIM_MESHMULT    Č   @   a   TIM_MESHMULT%S    !  Ž   a   TIM_MESHMULT%A    ļ!  R   a   TIM_MESHMULT%B    "  R   a   TIM_MESHMULT%C    Z"  R   a   TIM_MESHMULT%D    Ŧ"  O       TIM_ZEROBC    û"  R   a   TIM_ZEROBC%A    M#  W       TIM_NLINCORR     ¤#  R   a   TIM_NLINCORR%N_    ö#  R   a   TIM_NLINCORR%N     H$  `       TIM_MEASURECORR #   ¨$  R   a   TIM_MEASURECORR%C1 #   ú$  R   a   TIM_MEASURECORR%C2 #   L%  R   a   TIM_MEASURECORR%C3    %  H       TIM_CHECK_CGCE    æ%  H       TIM_NEW_TSTEP 