
  �6  Y   k820309    a          16.0        Z˭W                                                                                                           
       ../utility.F90 W90_UTILITY              UTILITY_INV3 UTILITY_INV2 UTILITY_RECIP_LATTICE UTILITY_METRIC UTILITY_COMPAR UTILITY_CART_TO_FRAC UTILITY_FRAC_TO_CART UTILITY_STRING_TO_COORD UTILITY_LOWERCASE UTILITY_STRIP UTILITY_ZGEMM UTILITY_TRANSLATE_HOME UTILITY_ROTATE UTILITY_MATMUL_DIAG UTILITY_ROTATE_DIAG UTILITY_COMMUTATOR_DIAG UTILITY_RE_TR UTILITY_IM_TR W0GAUSS WGAUSS UTILITY_DIAGONALIZE                      @                              
       DP                                                                                                       #         @                                                       #A    #B    #DET              
                                      	              
    p          p          p            p          p                                    D                                     	              
     p          p          p            p          p                                    D                                     
       #         @                                                       #A    #B 	   #DET 
             
                                                    
    p          p          p            p          p                                    D                                	                   
     p          p          p            p          p                                    D                                
     
       #         @                                                       #REAL_LAT    #RECIP_LAT    #VOLUME                                                                                                             
                                      	              
 	   p          p          p            p          p                                    D                                     	              
 
    p          p          p            p          p                                    D @                                   
       #         @                                                       #REAL_LAT    #RECIP_LAT    #REAL_METRIC    #RECIP_METRIC              
                                      	              
    p          p          p            p          p                                    
                                      	              
    p          p          p            p          p                                    D                                     	              
     p          p          p            p          p                                    D                                     	              
     p          p          p            p          p                          #         @                                                       #A    #B    #IFPOS    #IFNEG                                       
                                                    
    p          p            p                                    
                                                    
    p          p            p                                    D                                                       D                                             #         @                                                      #CART    #FRAC    #RECIP_LAT                                              
                                                    
    p          p            p                                    D                                                   
     p          p            p                                    
                                      	              
    p          p          p            p          p                          #         @                                                      #FRAC    #CART    #REAL_LAT               
                                                    
    p          p            p                                    D                                                   
     p          p            p                                    
                                       	              
    p          p          p            p          p                          #         @                                   !                    #STRING_TMP "   #OUTVEC #                                                 
  @                              "     x                                D                                #                   
     p          p            p                          $         @                                $     x                      #STRING %                                                   
  @                              %                    1 $         @                                &     x                      #STRING '                                               
  @                              '                    1 #         @                                   (                    #C )   #A +   #TRANSA ,   #B -   #TRANSB .   #N *                                                                  D @                              )                           p        5 � p        r *   p          5 � p        r *     5 � p        r *       5 � p        r *     5 � p        r *                              
@ @                              +                          p        5 � p        r *   p          5 � p        r *     5 � p        r *       5 � p        r *     5 � p        r *                               
@ @                              ,                                    
@ @                              -                          p        5 � p        r *   p          5 � p        r *     5 � p        r *       5 � p        r *     5 � p        r *                               
@ @                              .                                     
@ @                               *           #         @                                   /                    #VEC 0   #REAL_LAT 1   #RECIP_LAT 2             
D @                              0                   
     p          p            p                                    
  @                              1     	              
    p          p          p            p          p                                    
  @                              2     	              
    p          p          p            p          p                          (        `                                3                    (                    #MAT 4   #ROT 6   #DIM 5     p        5 O p        p          5 O p          5 O p            5 O p          5 O p                                                                                          4                     &      p        5 � p        r 5   p          5 � p        r 5     5 � p        r 5       5 � p        r 5     5 � p        r 5                               @                              6                     '      p        5 � p        r 5   p          5 � p        r 5     5 � p        r 5       5 � p        r 5     5 � p        r 5                                                                5            (        `                               7                    ,                    #MAT1 8   #MAT2 :   #DIM 9   p          5 O p            5 O p                                                                                                                                8                     *      p        5 � p        r 9   p          5 � p        r 9     5 � p        r 9       5 � p        r 9     5 � p        r 9                                                              :                     +      p        5 � p        r 9   p          5 � p        r 9     5 � p        r 9       5 � p        r 9     5 � p        r 9                                                                9            (        `                                ;                    0                    #MAT <   #ROT >   #DIM =   p          5 O p            5 O p                                                                                               <                     .      p        5 � p        r =   p          5 � p        r =     5 � p        r =       5 � p        r =     5 � p        r =                              D @                              >                     /      p        5 � p        r =   p          5 � p        r =     5 � p        r =       5 � p        r =     5 � p        r =                               D @                               =            (        `                                ?                    4                    #MAT1 @   #MAT2 B   #DIM A   p          5 O p            5 O p                                                                   D @                              @                     2      p        5 � p        r A   p          5 � p        r A     5 � p        r A       5 � p        r A     5 � p        r A                              D @                              B                     3      p        5 � p        r A   p          5 � p        r A     5 � p        r A       5 � p        r A     5 � p        r A                               D @                               A            %         @                                C                    
       #MAT D                                                                                         @                              D                    5              &                   &                                           %         @                                E                    
       #MAT F                                                              @                              F                    6              &                   &                                           %         @                                G                    
       #X H   #N I                                              @                              H     
                                                  I            %         @                                J                    
       #X K   #N L                                                                           K     
                                                  L            #         @                                   M                    #MAT N   #DIM O   #EIG P   #ROT Q                                                                         
                                 N                          p        5 � p        r O   p          5 � p        r O     5 � p        r O       5 � p        r O     5 � p        r O                               
@ @                               O                    D @                              P                    
     p          5 � p        r O       5 � p        r O                              D @                              Q                           p        5 � p        r O   p          5 � p        r O     5 � p        r O       5 � p        r O     5 � p        r O                     %         @   @                           R                    
       #X S                            
  @                              S     
      %         @   @                           T                    
       #X U                                
                                 U     
      %         @   @                           V                    
       #X W                             
  @                              W     
         �   #      fn#fn !   �   s  b   uapp(W90_UTILITY    6  C   J  W90_CONSTANTS !   y  p       DP+W90_CONSTANTS    �  _       UTILITY_INV3    H  �   a   UTILITY_INV3%A    �  �   a   UTILITY_INV3%B !   �  @   a   UTILITY_INV3%DET    �  _       UTILITY_INV2    O  �   a   UTILITY_INV2%A      �   a   UTILITY_INV2%B !   �  @   a   UTILITY_INV2%DET &   �  �       UTILITY_RECIP_LATTICE /   �  �   a   UTILITY_RECIP_LATTICE%REAL_LAT 0   {  �   a   UTILITY_RECIP_LATTICE%RECIP_LAT -   /	  @   a   UTILITY_RECIP_LATTICE%VOLUME    o	  �       UTILITY_METRIC (   �	  �   a   UTILITY_METRIC%REAL_LAT )   �
  �   a   UTILITY_METRIC%RECIP_LAT +   _  �   a   UTILITY_METRIC%REAL_METRIC ,     �   a   UTILITY_METRIC%RECIP_METRIC    �  �       UTILITY_COMPAR !   L  �   a   UTILITY_COMPAR%A !   �  �   a   UTILITY_COMPAR%B %   t  @   a   UTILITY_COMPAR%IFPOS %   �  @   a   UTILITY_COMPAR%IFNEG %   �  �       UTILITY_CART_TO_FRAC *     �   a   UTILITY_CART_TO_FRAC%CART *     �   a   UTILITY_CART_TO_FRAC%FRAC /   �  �   a   UTILITY_CART_TO_FRAC%RECIP_LAT %   [  j       UTILITY_FRAC_TO_CART *   �  �   a   UTILITY_FRAC_TO_CART%FRAC *   Y  �   a   UTILITY_FRAC_TO_CART%CART .   �  �   a   UTILITY_FRAC_TO_CART%REAL_LAT (   �  �       UTILITY_STRING_TO_COORD 3   )  P   a   UTILITY_STRING_TO_COORD%STRING_TMP /   y  �   a   UTILITY_STRING_TO_COORD%OUTVEC "     �       UTILITY_LOWERCASE )   �  L   a   UTILITY_LOWERCASE%STRING    �  ~       UTILITY_STRIP %   Y  L   a   UTILITY_STRIP%STRING    �  �       UTILITY_ZGEMM     W  $  a   UTILITY_ZGEMM%C     {  $  a   UTILITY_ZGEMM%A %   �  P   a   UTILITY_ZGEMM%TRANSA     �  $  a   UTILITY_ZGEMM%B %     P   a   UTILITY_ZGEMM%TRANSB     c  @   a   UTILITY_ZGEMM%N '   �  n       UTILITY_TRANSLATE_HOME +     �   a   UTILITY_TRANSLATE_HOME%VEC 0   �  �   a   UTILITY_TRANSLATE_HOME%REAL_LAT 1   Y  �   a   UTILITY_TRANSLATE_HOME%RECIP_LAT      >      UTILITY_ROTATE #   K  $  a   UTILITY_ROTATE%MAT #   o   $  a   UTILITY_ROTATE%ROT #   �!  @   a   UTILITY_ROTATE%DIM $   �!        UTILITY_MATMUL_DIAG )   �"  $  a   UTILITY_MATMUL_DIAG%MAT1 )   $  $  a   UTILITY_MATMUL_DIAG%MAT2 (   )%  @   a   UTILITY_MATMUL_DIAG%DIM $   i%  �       UTILITY_ROTATE_DIAG (   T&  $  a   UTILITY_ROTATE_DIAG%MAT (   x'  $  a   UTILITY_ROTATE_DIAG%ROT (   �(  @   a   UTILITY_ROTATE_DIAG%DIM (   �(  �       UTILITY_COMMUTATOR_DIAG -   �)  $  a   UTILITY_COMMUTATOR_DIAG%MAT1 -   �*  $  a   UTILITY_COMMUTATOR_DIAG%MAT2 ,   ,  @   a   UTILITY_COMMUTATOR_DIAG%DIM    U,  �       UTILITY_RE_TR "   �,  �   a   UTILITY_RE_TR%MAT    �-  �       UTILITY_IM_TR "   (.  �   a   UTILITY_IM_TR%MAT    �.  ~       W0GAUSS    J/  @   a   W0GAUSS%X    �/  @   a   W0GAUSS%N    �/  |       WGAUSS    F0  @   a   WGAUSS%X    �0  @   a   WGAUSS%N $   �0  �       UTILITY_DIAGONALIZE (   o1  $  a   UTILITY_DIAGONALIZE%MAT (   �2  @   a   UTILITY_DIAGONALIZE%DIM (   �2  �   a   UTILITY_DIAGONALIZE%EIG (   �3  $  a   UTILITY_DIAGONALIZE%ROT    �4  f       QE_ERF    5  @   a   QE_ERF%X    Q5  j       GAUSS_FREQ    �5  @   a   GAUSS_FREQ%X    �5  g       QE_ERFC    b6  @   a   QE_ERFC%X 