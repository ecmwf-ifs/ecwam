      SUBROUTINE CIGETDEAC

! ----------------------------------------------------------------------

!**** *CIGETDEAC* - DEFINE THE DIMENSIONLESS ENERGY ATTENUATION COEFFICIENT 

!*    PURPOSE.
!     --------

!       CIGETDEAC DEFINES THE DIMENSIONLESS ENERGY ATTENUATION COEFFICIENT

!**   INTERFACE.
!     ----------

!       *CALL* *CIGETDEAC

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCES.                                                       
!     -----------  

!     KOHOUT AND MEYLAN, 2008: JGR, 113, doi:10.1029/2007JC004424

! ----------------------------------------------------------------------

      USE YOWICE   , ONLY : NICT    ,NICH    ,TICMIN  ,HICMIN   ,
     &              DTIC   ,DHIC    ,CIDEAC
      USE YOWTEST  , ONLY : IU06     ,ITEST

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: IT, IH
      REAL :: DHI, DCI
      REAL :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('CIGETDEAC',0,ZHOOK_HANDLE)
#endif

! I have kindly received the values from Figure 6 of Kohout and Meylan
! from Alison Kohout
! I have also extrapolated the values for thickness=0.2 and
! for periods between 1 and 5 sec by assuming that the values asymptotes to -3 to -1
! for 1 sec. linear dependence on ice thickness

!     ice thickness discretisation
      NICH=36
      HICMIN=0.2
      DHIC=0.1
      
!     wave period discretisation
      NICT=16
      TICMIN=1.0
      DTIC=1.0

      IF(ALLOCATED(CIDEAC)) DEALLOCATE(CIDEAC)
      ALLOCATE(CIDEAC(NICT,NICH))

!     assumption:
      CIDEAC( 1, 1)=-2.00
      CIDEAC( 1,NICH)=-1.00
      DHI=CIDEAC( 1,NICH)-CIDEAC( 1, 1)
      DO IH=2,NICH-1
        CIDEAC( 1,IH)=CIDEAC( 1, 1)+(IH-1)*DHI/(NICH-1)
      ENDDO

      !!!! all values for h=0.20 and 0.30 were extrapolated
! h =    0.2000000    
      CIDEAC( 6, 1)=  -5.50
      CIDEAC( 7, 1)=  -7.20
      CIDEAC( 8, 1)=  -8.90
      CIDEAC( 9, 1)=  -10.10
      CIDEAC(10, 1)=  -10.75
      CIDEAC(11, 1)=  -11.10
      CIDEAC(12, 1)=  -11.30
      CIDEAC(13, 1)=  -11.50
      CIDEAC(14, 1)= -11.70
      CIDEAC(15, 1)= -12.35
      CIDEAC(16, 1)= -13.00

! h =    0.3000000    
      CIDEAC( 6, 2)=  -4.84
      CIDEAC( 7, 2)=  -6.35
      CIDEAC( 8, 2)=  -7.99
      CIDEAC( 9, 2)=  -9.46
      CIDEAC(10, 2)=  -10.43
      CIDEAC(11, 2)=  -10.93
      CIDEAC(12, 2)=  -11.15
      CIDEAC(13, 2)=  -11.37
      CIDEAC(14, 2)= -11.59
      CIDEAC(15, 2)= -12.00
      CIDEAC(16, 2)= -13.44

! h =    0.4000000    
      CIDEAC( 6, 3)=   -4.179189    
      CIDEAC( 7, 3)=   -5.494833    
      CIDEAC( 8, 3)=   -7.073360    
      CIDEAC( 9, 3)=   -8.809628    
      CIDEAC(10, 3)=   -10.11117    
      CIDEAC(11, 3)=   -10.75490    
      CIDEAC(12, 3)=   -11.01079    
      CIDEAC(13, 3)=   -11.23960    
      CIDEAC(14, 3)=   -11.46686    
      CIDEAC(15, 3)=   -11.64766    
      CIDEAC(16, 3)=   -11.88710    
 
! h =    0.5000000    
      CIDEAC( 6, 4)=   -3.671164    
      CIDEAC( 7, 4)=   -4.688528    
      CIDEAC( 8, 4)=   -6.063766    
      CIDEAC( 9, 4)=   -7.556777    
      CIDEAC(10, 4)=   -9.092730    
      CIDEAC(11, 4)=   -10.21499    
      CIDEAC(12, 4)=   -10.62448    
      CIDEAC(13, 4)=   -10.84617    
      CIDEAC(14, 4)=   -11.03255    
      CIDEAC(15, 4)=   -11.24719    
      CIDEAC(16, 4)=   -11.49634    
 
! h =    0.6000000    
      CIDEAC( 6, 5)=   -3.279724    
      CIDEAC( 7, 5)=   -4.138198    
      CIDEAC( 8, 5)=   -5.303435    
      CIDEAC( 9, 5)=   -6.634662    
      CIDEAC(10, 5)=   -8.086823    
      CIDEAC(11, 5)=   -9.434999    
      CIDEAC(12, 5)=   -10.24225    
      CIDEAC(13, 5)=   -10.53463    
      CIDEAC(14, 5)=   -10.77608    
      CIDEAC(15, 5)=   -10.94032    
      CIDEAC(16, 5)=   -11.14464    
 
! h =    0.7000000    
      CIDEAC( 6, 6)=   -2.986051    
      CIDEAC( 7, 6)=   -3.809261    
      CIDEAC( 8, 6)=   -4.791394    
      CIDEAC( 9, 6)=   -5.896924    
      CIDEAC(10, 6)=   -7.261892    
      CIDEAC(11, 6)=   -8.680384    
      CIDEAC(12, 6)=   -9.695383    
      CIDEAC(13, 6)=   -10.25629    
      CIDEAC(14, 6)=   -10.57562    
      CIDEAC(15, 6)=   -10.75792    
      CIDEAC(16, 6)=   -10.88391    
 
! h =    0.8000000    
      CIDEAC( 6, 7)=   -2.802249    
      CIDEAC( 7, 7)=   -3.473820    
      CIDEAC( 8, 7)=   -4.312738    
      CIDEAC( 9, 7)=   -5.404909    
      CIDEAC(10, 7)=   -6.585416    
      CIDEAC(11, 7)=   -7.936078    
      CIDEAC(12, 7)=   -9.055226    
      CIDEAC(13, 7)=   -9.891881    
      CIDEAC(14, 7)=   -10.28646    
      CIDEAC(15, 7)=   -10.50197    
      CIDEAC(16, 7)=   -10.69553    

! h =    0.9000000    
      CIDEAC( 6, 8)=   -2.649911    
      CIDEAC( 7, 8)=   -3.208950    
      CIDEAC( 8, 8)=   -3.984321    
      CIDEAC( 9, 8)=   -4.919455    
      CIDEAC(10, 8)=   -6.056835    
      CIDEAC(11, 8)=   -7.323155    
      CIDEAC(12, 8)=   -8.486164    
      CIDEAC(13, 8)=   -9.498883    
      CIDEAC(14, 8)=   -10.05186    
      CIDEAC(15, 8)=   -10.26166    
      CIDEAC(16, 8)=   -10.44009    
 
! h =     1.000000    
      CIDEAC( 6, 9)=   -2.494633    
      CIDEAC( 7, 9)=   -3.058290    
      CIDEAC( 8, 9)=   -3.697340    
      CIDEAC( 9, 9)=   -4.548424    
      CIDEAC(10, 9)=   -5.556311    
      CIDEAC(11, 9)=   -6.724046    
      CIDEAC(12, 9)=   -7.960689    
      CIDEAC(13, 9)=   -8.995591    
      CIDEAC(14, 9)=   -9.717956    
      CIDEAC(15, 9)=   -10.12654    
      CIDEAC(16, 9)=   -10.34523    
 
! h =     1.100000    
      CIDEAC( 6,10)=   -2.399249    
      CIDEAC( 7,10)=   -2.882891    
      CIDEAC( 8,10)=   -3.501052    
      CIDEAC( 9,10)=   -4.253071    
      CIDEAC(10,10)=   -5.189890    
      CIDEAC(11,10)=   -6.225914    
      CIDEAC(12,10)=   -7.413078    
      CIDEAC(13,10)=   -8.483686    
      CIDEAC(14,10)=   -9.360973    
      CIDEAC(15,10)=   -9.923790    
      CIDEAC(16,10)=   -10.21649    
 
! h =     1.200000    
      CIDEAC( 6,11)=   -2.314447    
      CIDEAC( 7,11)=   -2.792313    
      CIDEAC( 8,11)=   -3.336744    
      CIDEAC( 9,11)=   -4.020543    
      CIDEAC(10,11)=   -4.842176    
      CIDEAC(11,11)=   -5.899529    
      CIDEAC(12,11)=   -6.954353    
      CIDEAC(13,11)=   -8.057040    
      CIDEAC(14,11)=   -8.976917    
      CIDEAC(15,11)=   -9.661180    
      CIDEAC(16,11)=   -10.01354    
 
! h =     1.300000    
      CIDEAC( 6,12)=   -2.243207    
      CIDEAC( 7,12)=   -2.647049    
      CIDEAC( 8,12)=   -3.217140    
      CIDEAC( 9,12)=   -3.814464    
      CIDEAC(10,12)=   -4.590047    
      CIDEAC(11,12)=   -5.466197    
      CIDEAC(12,12)=   -6.520925    
      CIDEAC(13,12)=   -7.642193    
      CIDEAC(14,12)=   -8.617725    
      CIDEAC(15,12)=   -9.381130    
      CIDEAC(16,12)=   -9.848615    
 
! h =     1.400000    
      CIDEAC( 6,13)=   -2.169400    
      CIDEAC( 7,13)=   -2.576199    
      CIDEAC( 8,13)=   -3.010994    
      CIDEAC( 9,13)=   -3.624323    
      CIDEAC(10,13)=   -4.324925    
      CIDEAC(11,13)=   -5.187501    
      CIDEAC(12,13)=   -6.182012    
      CIDEAC(13,13)=   -7.190685    
      CIDEAC(14,13)=   -8.241459    
      CIDEAC(15,13)=   -9.042115    
      CIDEAC(16,13)=   -9.622642    

! h =     1.500000    
      CIDEAC( 6,14)=   -2.062091    
      CIDEAC( 7,14)=   -2.452876    
      CIDEAC( 8,14)=   -2.952611    
      CIDEAC( 9,14)=   -3.506958    
      CIDEAC(10,14)=   -4.177629    
      CIDEAC(11,14)=   -4.920455    
      CIDEAC(12,14)=   -5.839159    
      CIDEAC(13,14)=   -6.873662    
      CIDEAC(14,14)=   -7.827965    
      CIDEAC(15,14)=   -8.873965    
      CIDEAC(16,14)=   -9.413249    

! h =     1.600000    
      CIDEAC( 6,15)=   -1.996322    
      CIDEAC( 7,15)=   -2.392122    
      CIDEAC( 8,15)=   -2.823834    
      CIDEAC( 9,15)=   -3.359782    
      CIDEAC(10,15)=   -3.958120    
      CIDEAC(11,15)=   -4.667717    
      CIDEAC(12,15)=   -5.541727    
      CIDEAC(13,15)=   -6.549209    
      CIDEAC(14,15)=   -7.520937    
      CIDEAC(15,15)=   -8.429287    
      CIDEAC(16,15)=   -9.221933    
 
! h =     1.700000    
      CIDEAC( 6,16)=   -1.955446    
      CIDEAC( 7,16)=   -2.320779    
      CIDEAC( 8,16)=   -2.715107    
      CIDEAC( 9,16)=   -3.211328    
      CIDEAC(10,16)=   -3.808961    
      CIDEAC(11,16)=   -4.510885    
      CIDEAC(12,16)=   -5.320290    
      CIDEAC(13,16)=   -6.226068    
      CIDEAC(14,16)=   -7.185133    
      CIDEAC(15,16)=   -8.094091    
      CIDEAC(16,16)=   -8.907555    

! h =     1.800000    
      CIDEAC( 6,17)=   -1.919889    
      CIDEAC( 7,17)=   -2.260443    
      CIDEAC( 8,17)=   -2.674223    
      CIDEAC( 9,17)=   -3.108331    
      CIDEAC(10,17)=   -3.647704    
      CIDEAC(11,17)=   -4.258914    
      CIDEAC(12,17)=   -5.110578    
      CIDEAC(13,17)=   -5.947557    
      CIDEAC(14,17)=   -6.857039    
      CIDEAC(15,17)=   -7.801445    
      CIDEAC(16,17)=   -8.697256    
 
! h =     1.900000    
      CIDEAC( 6,18)=   -1.866967    
      CIDEAC( 7,18)=   -2.277320    
      CIDEAC( 8,18)=   -2.620277    
      CIDEAC( 9,18)=   -3.053688    
      CIDEAC(10,18)=   -3.530046    
      CIDEAC(11,18)=   -4.170623    
      CIDEAC(12,18)=   -4.870700    
      CIDEAC(13,18)=   -5.723269    
      CIDEAC(14,18)=   -6.651424    
      CIDEAC(15,18)=   -7.524860    
      CIDEAC(16,18)=   -8.409969    
 
! h =     2.000000    
      CIDEAC( 6,19)=   -1.822846    
      CIDEAC( 7,19)=   -2.088338    
      CIDEAC( 8,19)=   -2.531257    
      CIDEAC( 9,19)=   -2.942987    
      CIDEAC(10,19)=   -3.422226    
      CIDEAC(11,19)=   -4.036695    
      CIDEAC(12,19)=   -4.701579    
      CIDEAC(13,19)=   -5.505970    
      CIDEAC(14,19)=   -6.318645    
      CIDEAC(15,19)=   -7.290503    
      CIDEAC(16,19)=   -8.141773    

! h =     2.100000    
      CIDEAC( 6,20)=   -1.779950    
      CIDEAC( 7,20)=   -2.108976    
      CIDEAC( 8,20)=   -2.475516    
      CIDEAC( 9,20)=   -2.862222    
      CIDEAC(10,20)=   -3.349712    
      CIDEAC(11,20)=   -3.842627    
      CIDEAC(12,20)=   -4.558055    
      CIDEAC(13,20)=   -5.257342    
      CIDEAC(14,20)=   -6.072487    
      CIDEAC(15,20)=   -7.056084    
      CIDEAC(16,20)=   -7.865724    
 
! h =     2.200000    
      CIDEAC( 6,21)=   -1.758808    
      CIDEAC( 7,21)=   -2.097306    
      CIDEAC( 8,21)=   -2.397066    
      CIDEAC( 9,21)=   -2.773985    
      CIDEAC(10,21)=   -3.249429    
      CIDEAC(11,21)=   -3.731577    
      CIDEAC(12,21)=   -4.400121    
      CIDEAC(13,21)=   -5.056503    
      CIDEAC(14,21)=   -5.947149    
      CIDEAC(15,21)=   -6.758241    
      CIDEAC(16,21)=   -7.665355    
 
! h =     2.300000    
      CIDEAC( 6,22)=   -1.742336    
      CIDEAC( 7,22)=   -2.030318    
      CIDEAC( 8,22)=   -2.380843    
      CIDEAC( 9,22)=   -2.738259    
      CIDEAC(10,22)=   -3.130276    
      CIDEAC(11,22)=   -3.666918    
      CIDEAC(12,22)=   -4.254450    
      CIDEAC(13,22)=   -4.941644    
      CIDEAC(14,22)=   -5.660053    
      CIDEAC(15,22)=   -6.561189    
      CIDEAC(16,22)=   -7.462294    
 
! h =     2.400000    
      CIDEAC( 6,23)=   -1.682265    
      CIDEAC( 7,23)=   -1.974037    
      CIDEAC( 8,23)=   -2.300742    
      CIDEAC( 9,23)=   -2.670091    
      CIDEAC(10,23)=   -3.057821    
      CIDEAC(11,23)=   -3.594969    
      CIDEAC(12,23)=   -4.102533    
      CIDEAC(13,23)=   -4.751281    
      CIDEAC(14,23)=   -5.549400    
      CIDEAC(15,23)=   -6.337473    
      CIDEAC(16,23)=   -7.156771    
 
! h =     2.500000    
      CIDEAC( 6,24)=   -1.648687    
      CIDEAC( 7,24)=   -1.903361    
      CIDEAC( 8,24)=   -2.280064    
      CIDEAC( 9,24)=   -2.629137    
      CIDEAC(10,24)=   -3.024701    
      CIDEAC(11,24)=   -3.445110    
      CIDEAC(12,24)=   -3.987420    
      CIDEAC(13,24)=   -4.637731    
      CIDEAC(14,24)=   -5.365729    
      CIDEAC(15,24)=   -6.134120    
      CIDEAC(16,24)=   -6.950096    

! h =     2.600000    
      CIDEAC( 6,25)=   -1.641218    
      CIDEAC( 7,25)=   -1.956000    
      CIDEAC( 8,25)=   -2.195282    
      CIDEAC( 9,25)=   -2.547474    
      CIDEAC(10,25)=   -2.952879    
      CIDEAC(11,25)=   -3.382897    
      CIDEAC(12,25)=   -3.878702    
      CIDEAC(13,25)=   -4.551814    
      CIDEAC(14,25)=   -5.228397    
      CIDEAC(15,25)=   -5.965115    
      CIDEAC(16,25)=   -6.785750    
 
! h =     2.700000    
      CIDEAC( 6,26)=   -1.600150    
      CIDEAC( 7,26)=   -1.844258    
      CIDEAC( 8,26)=   -2.257075    
      CIDEAC( 9,26)=   -2.519011    
      CIDEAC(10,26)=   -2.883503    
      CIDEAC(11,26)=   -3.288872    
      CIDEAC(12,26)=   -3.823005    
      CIDEAC(13,26)=   -4.376071    
      CIDEAC(14,26)=   -5.016023    
      CIDEAC(15,26)=   -5.787023    
      CIDEAC(16,26)=   -6.621093    
 
! h =     2.800000    
      CIDEAC( 6,27)=   -1.586459    
      CIDEAC( 7,27)=   -1.823506    
      CIDEAC( 8,27)=   -2.191002    
      CIDEAC( 9,27)=   -2.443355    
      CIDEAC(10,27)=   -2.839071    
      CIDEAC(11,27)=   -3.269969    
      CIDEAC(12,27)=   -3.717166    
      CIDEAC(13,27)=   -4.261279    
      CIDEAC(14,27)=   -4.884097    
      CIDEAC(15,27)=   -5.582966    
      CIDEAC(16,27)=   -6.413867    
 
! h =     2.900000    
      CIDEAC( 6,28)=   -1.530613    
      CIDEAC( 7,28)=   -1.881014    
      CIDEAC( 8,28)=   -2.084314    
      CIDEAC( 9,28)=   -2.427499    
      CIDEAC(10,28)=   -2.796773    
      CIDEAC(11,28)=   -3.155616    
      CIDEAC(12,28)=   -3.628541    
      CIDEAC(13,28)=   -4.160691    
      CIDEAC(14,28)=   -4.748525    
      CIDEAC(15,28)=   -5.441663    
      CIDEAC(16,28)=   -6.239977    
 
! h =     3.000000    
      CIDEAC( 6,29)=   -1.520659    
      CIDEAC( 7,29)=   -1.827379    
      CIDEAC( 8,29)=   -2.099156    
      CIDEAC( 9,29)=   -2.403214    
      CIDEAC(10,29)=   -2.764410    
      CIDEAC(11,29)=   -3.155527    
      CIDEAC(12,29)=   -3.555711    
      CIDEAC(13,29)=   -4.072697    
      CIDEAC(14,29)=   -4.656217    
      CIDEAC(15,29)=   -5.335746    
      CIDEAC(16,29)=   -6.086100    
 
! h =     3.100000    
      CIDEAC( 6,30)=   -1.468517    
      CIDEAC( 7,30)=   -1.802037    
      CIDEAC( 8,30)=   -2.039742    
      CIDEAC( 9,30)=   -2.369589    
      CIDEAC(10,30)=   -2.647677    
      CIDEAC(11,30)=   -3.087644    
      CIDEAC(12,30)=   -3.489749    
      CIDEAC(13,30)=   -3.964992    
      CIDEAC(14,30)=   -4.532036    
      CIDEAC(15,30)=   -5.117778    
      CIDEAC(16,30)=   -5.921745    
 
! h =     3.200000    
      CIDEAC( 6,31)=   -1.462488    
      CIDEAC( 7,31)=   -1.748595    
      CIDEAC( 8,31)=   -1.976536    
      CIDEAC( 9,31)=   -2.330469    
      CIDEAC(10,31)=   -2.633314    
      CIDEAC(11,31)=   -2.967257    
      CIDEAC(12,31)=   -3.411971    
      CIDEAC(13,31)=   -3.881270    
      CIDEAC(14,31)=   -4.440731    
      CIDEAC(15,31)=   -5.080121    
      CIDEAC(16,31)=   -5.782035    
 
! h =     3.300000    
      CIDEAC( 6,32)=   -1.404058    
      CIDEAC( 7,32)=   -1.682052    
      CIDEAC( 8,32)=   -1.989879    
      CIDEAC( 9,32)=   -2.303284    
      CIDEAC(10,32)=   -2.609248    
      CIDEAC(11,32)=   -2.970038    
      CIDEAC(12,32)=   -3.335601    
      CIDEAC(13,32)=   -3.827919    
      CIDEAC(14,32)=   -4.264277    
      CIDEAC(15,32)=   -4.969507    
      CIDEAC(16,32)=   -5.624526    
 
! h =     3.400000    
      CIDEAC( 6,33)=   -1.421658    
      CIDEAC( 7,33)=   -1.694868    
      CIDEAC( 8,33)=   -1.939824    
      CIDEAC( 9,33)=   -2.256738    
      CIDEAC(10,33)=   -2.563078    
      CIDEAC(11,33)=   -2.854945    
      CIDEAC(12,33)=   -3.300406    
      CIDEAC(13,33)=   -3.748025    
      CIDEAC(14,33)=   -4.209355    
      CIDEAC(15,33)=   -4.798466    
      CIDEAC(16,33)=   -5.484742    
 
! h =     3.500000    
      CIDEAC( 6,34)=   -1.433080    
      CIDEAC( 7,34)=   -1.658592    
      CIDEAC( 8,34)=   -1.953619    
      CIDEAC( 9,34)=   -2.189728    
      CIDEAC(10,34)=   -2.483670    
      CIDEAC(11,34)=   -2.876348    
      CIDEAC(12,34)=   -3.227396    
      CIDEAC(13,34)=   -3.643594    
      CIDEAC(14,34)=   -4.174257    
      CIDEAC(15,34)=   -4.713649    
      CIDEAC(16,34)=   -5.365061    

! h =     3.600000    
      CIDEAC( 6,35)=   -1.404438    
      CIDEAC( 7,35)=   -1.603982    
      CIDEAC( 8,35)=   -1.976813    
      CIDEAC( 9,35)=   -2.186624    
      CIDEAC(10,35)=   -2.514024    
      CIDEAC(11,35)=   -2.848864    
      CIDEAC(12,35)=   -3.172881    
      CIDEAC(13,35)=   -3.602279    
      CIDEAC(14,35)=   -4.101719    
      CIDEAC(15,35)=   -4.649209    
      CIDEAC(16,35)=   -5.305425    
 
! h =     3.700000    
      CIDEAC( 6,36)=   -1.392771    
      CIDEAC( 7,36)=   -1.646687    
      CIDEAC( 8,36)=   -1.921211    
      CIDEAC( 9,36)=   -2.149506    
      CIDEAC(10,36)=   -2.463910    
      CIDEAC(11,36)=   -2.780142    
      CIDEAC(12,36)=   -3.117836    
      CIDEAC(13,36)=   -3.524015    
      CIDEAC(14,36)=   -4.022181    
      CIDEAC(15,36)=   -4.536423    
      CIDEAC(16,36)=   -5.194000    
 
!     extrapolation:
      DO IH=1,NICH
        DCI=CIDEAC(6,IH)-CIDEAC(1,IH)
        DO IT=2,5
          CIDEAC(IT,IH)=CIDEAC(1,IH)+DCI*(IT-1)*DTIC/(5*DTIC)
        ENDDO
      ENDDO

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('CIGETDEAC',1,ZHOOK_HANDLE)
#endif
      RETURN
      END SUBROUTINE CIGETDEAC
