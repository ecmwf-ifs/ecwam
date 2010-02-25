      PROGRAM intwaminput

! ----------------------------------------------------------------------

!**** *INTWAMINPUT* -

!     J. BIDLOT    ECMWF  JULY 2000 

!*    PURPOSE.
!     --------
!     TO TRANSFORM A GRIB INPUT WAVE FIELD (SPECTRA OR OTHERS) INTO
!     A CORRESPONDING FIELD WITH GRID, DIRECTION AND FREQUENCY AS
!     SPECIFIED BY THE MODEL INPUT GRID FILE (gridglou).
!
!     IF THERE IS NOT NEED FOR INTERPOLATION, THE INPUT FILE WILL
!     SIMPLY BE RECOPIED TO THE OUTPUT DESTINATION. NO CODING/DECODING
!     SHOULD OCCUR.

!     IN CASE THE INPUT FIELD IS PARAMETER 251 (WAVE SPECTRA) A FILE
!     NAMED wave_spectral_resolution WILL BE OUTPUT CONTAINING THE
!     WAVE SPECTRAL RESOLUTION (number of directions and frequencies)
!     AND WHETHER OR NOT THE FIELD HAS TO BE INTERPOLATED.

!**   INTERFACE.
!     ----------
!     gridglou must be present and connected to unit 7 (i.e fort.7)
!     The input should be in input_field
!     The output will be placed in output_field

!     METHOD.
!     -------
!     EXTERNALS.
!     ----------
!     REFERENCES
!     ----------
!       NONE.
! ----------------------------------------------------------------------

      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRIBHD, ONLY : NKSEK0   ,NKSEK1   ,NKSEK2   ,NKSEK3   ,
     &            NKSEK4   ,NPSEK2   ,NPSEK3   ,KSEK0    ,KSEK1    ,
     &            KSEK2    ,KSEK3    ,KSEK4    ,PSEK2    ,PSEK3    ,
     &            PPEPS    ,PPREC    ,NTENCODE ,NGRBRESI ,
     &            NGRBRESS ,HOPERI   ,HOPERS   ,LGRHDIFS ,LNEWLVTP
      USE YOWGRID  , ONLY : DELPHI   ,NLONRGG  ,IJL
      USE YOWMAP   , ONLY : IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP     ,
     &            KTAG     ,NPRECR   ,NPRECI
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY      ,
     &            NBLO     ,NIBLO    ,CLDOMAIN
      USE YOWSTAT  , ONLY : MARSTYPE ,YCLASS   ,YEXPVER  ,ISTREAM  ,
     &            IREFRA   ,NENSFNB  ,NTOTENS
      USE YOWSPEC  , ONLY : NSTART   ,NEND
      USE YOWTEST  , ONLY : IU06     ,ITEST

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, PARAMETER :: NBIT = 180000

      INTEGER :: I4(2), IUIN, KRET, KPLENG, ISIZE, KLEN
      INTEGER :: KLENG, KLENP, KWORD, NC, NR, NTOT, I, J, JRGG, IR
      INTEGER :: IREPR, IPARAM, KANG, KFRE, IALLFLD, JSN
      INTEGER :: IFR1, KFR1, IFR, IC, MMSHIFT, MLAST 
      INTEGER :: IFORP, KZLEV, K, M, KKK, MMM, JPARAM, KKKPR, NGYFULL
      INTEGER :: IPERIODIC, NRFULL, JSTREAM, ISTART, ISTOP, IDUM
      INTEGER :: IU07, ILEVTYPE, KSTREAM
      INTEGER :: KAMONOP, KAMOEAP, KAMOSOP, KAMOWEP
      INTEGER :: KRMONOP, KRMOEAP, KRMOSOP, KRMOWEP 
      INTEGER :: KSEC0(NKSEK0), KSEC2(NKSEK2),
     &           KSEC3(NKSEK3), KSEC4(NKSEK4)
      INTEGER :: I_GET_UNIT

      INTEGER, ALLOCATABLE :: KSEC1(:), KGRIB(:), KLONRGG(:)

      REAL X4(2)
      REAL :: RMONOP, RMOSOP, RMOEAP, RMOWEP
      REAL ::  DELLO
      REAL :: ONETHIRD, TWOTHIRD
      REAL ::  PSEC2(NPSEK2),PSEC3(NPSEK3)
      REAL, ALLOCATABLE :: PSEC4(:)
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FIELD, TEMP

      CHARACTER*  2 MARSFCTYPE
      CHARACTER*  4 CSTREAM 
      CHARACTER* 14 CDATE
      CHARACTER* 11 IFILENAME
      CHARACTER* 12 OFILENAME

      LOGICAL LLINTERPOL, LLNONWAVE, LASTREAM
      LOGICAL LLFR1OK

! ----------------------------------------------------------------------

      IU06 = 6
      IU07 = I_GET_UNIT(IU06, 'wam_grid_tables', 'r', 'u', 0)

      WRITE(IU06,*) ' INTWAMINPUT STARTED '

      CALL INIWCST

!     1.1 INITIALISATION OF VARIABLES WITH DEFAULT VALUES
!         ---------------------------------------------------
      ITEST  =   0 
      LMESSPASS = .FALSE.
      MARSTYPE = 'an'
      YCLASS   = 'od'
      YEXPVER  = '0001' 
      NENSFNB = 0
      NTOTENS = 0
      ISTREAM = 1045 !!! is changed to an ifs stream also change LNEWLVTP 
      LGRHDIFS=.FALSE.
      LNEWLVTP=.FALSE.

      IRANK=1
      NPROC=1
      X4=1.
      NPRECR = KIND(X4)
      I4=1
      NPRECI = KIND(I4)
      KTAG=100

      ONETHIRD = 1./3.
      TWOTHIRD = 2./3.

      ALLOCATE (NSTART(NPROC),NEND(NPROC))

! ----------------------------------------------------------------------

!*    2. READ PREPROC OUTPUT.
!        --------------------

      CALL READPRE (IU07)

      NINF=1
      NSUP=NIBLO

      NSTART=1
      NEND=IJL(1)

!     SET GRIB HEADERS FOR INPUTS/OUTPUTS

      CALL SETWGRIBHD

      ALLOCATE(KSEC1(NKSEK1))

! ----------------------------------------------------------------------

!     3. READ GRIB INPUT.
!        ---------------

       IFILENAME='input_field'
       OFILENAME='output_field'
  

       CALL PBOPEN(IUIN,IFILENAME,'r',KRET)
       IF(KRET.LT.0) THEN
         IF(KRET.EQ.-1) THEN
            WRITE(IU06,*) '***********************************'
            WRITE(IU06,*) '*  COULD NOT OPEN '//IFILENAME
            WRITE(IU06,*) '***********************************'
            CALL ABORT1
         ENDIF
         IF(KRET.EQ.-2) THEN
           WRITE(IU06,*) '***********************************'
           WRITE(IU06,*) '*  INVALID FILENAME'//IFILENAME
           WRITE(IU06,*) '***********************************'
           CALL ABORT1
         ENDIF
       ENDIF


       ISIZE=NBIT
1001   KPLENG=ISIZE*KIND(ISIZE)
       IF(.NOT.ALLOCATED(KGRIB)) ALLOCATE(KGRIB(ISIZE))
       CALL PBGRIB(IUIN,KGRIB,KPLENG,KLEN,KRET)
       IF(KRET.EQ.-1) THEN
          WRITE(IU06,*) '***********************************'
          WRITE(IU06,*) '*  PBGRIB :  END OF FILE ENCOUNTED'
          WRITE(IU06,*) '***********************************'
          CALL ABORT1
       ENDIF
       IF(KRET.EQ.-2) THEN
         WRITE(IU06,*) '***********************************'
         WRITE(IU06,*) '*  PBGRIB : FILE HANDLING ERROR'
         WRITE(IU06,*) '***********************************'
         CALL ABORT1
       ENDIF
       IF(KRET.EQ.-3) THEN
         CALL KGRIBSIZE(IU06, IUIN, KLEN, ISIZE, 'INTWAMINPUT')
         DEALLOCATE(KGRIB)
         GOTO 1001
       ENDIF

!      GET SIZE OF THE UNPACKED DATA
       ALLOCATE(PSEC4(1))
       KLENG = ISIZE
       KLENP = 1
       KRET  = 1
       KSEC3 = KSEK3
       PSEC3 = PSEK3

!      CALL GRSDBG(1)
       CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &              PSEC4, KLENP, KGRIB, KLENG , KWORD, 'I', KRET)
       IF(KRET.GT.0) THEN
         WRITE(IU06,*) '***************************************'
         WRITE(IU06,*) '*                                     *'
         WRITE(IU06,*) '*  GRIBEX FAILED WITH ERROR CODE ',KRET
         WRITE(IU06,*) '*                                     *'
         WRITE(IU06,*) '*     THE PROGRAM ABORTS              *'
         WRITE(IU06,*) '***************************************'
         CALL ABORT1
       ENDIF

       NR=KSEC2(3)
       NRFULL=NR
       IF(NR.GT.NKSEK2-22) THEN
         WRITE(IU06,*) '***************************************'
         WRITE(IU06,*) '*                                     *'
         WRITE(IU06,*) '*  THE DIMENSION FOR KSEC2 ',NR+22
         WRITE(IU06,*) '*  IS LARGER THAN THE PRESCRIBED ONE ',NKSEK2
         WRITE(IU06,*) '*                                     *'
         WRITE(IU06,*) '*     THE PROGRAM ABORTS              *'
         WRITE(IU06,*) '***************************************'
         CALL ABORT1
       ENDIF
       JRGG = KSEC2(17)

       IF(JRGG.EQ.1) THEN
         NC=0
         NTOT=0
         DO J=1,NR
           NC = MAX(NC,KSEC2(22+J))
           NTOT=NTOT+KSEC2(22+J)
         ENDDO

         IR=0
         DO J=1,NR
           IF(KSEC2(22+J).NE.0) IR=IR+1
         ENDDO
         NR=IR
         KLEN=NTOT
       ELSEIF(JRGG.EQ.0) THEN
         NC=KSEC2(2)
         KLEN = NC*NR
       ELSE
         WRITE(IU06,*)
     &'   STRUCTURE OF THE FIELD NOT KNOWN'
         CALL ABORT1
       ENDIF

       IREPR=KSEC2(1)

       IF(IREPR.NE.0.AND.IREPR.NE.4) THEN
         WRITE(IU06,*) '***************************************'
         WRITE(IU06,*) '*                                     *'
         WRITE(IU06,*) '*  UNKNOWN GRID REPRESENTATION = ',IREPR
         WRITE(IU06,*) '*  CAN ONLY DEAL WITH                 *'
         WRITE(IU06,*) '*  LATITUDE/LONGITUDE GRID (IREPR=0)  *'
         WRITE(IU06,*) '*   OR GAUSSIAN (IREPR=4)             *'
         WRITE(IU06,*) '*                                     *'
         WRITE(IU06,*) '*     THE PROGRAM ABORTS              *'
         WRITE(IU06,*) '***************************************'
         CALL ABORT1
       ENDIF

!*     DETERMINE CODE FOR DATA FIELD TYPE.

       IPARAM = KSEC1(6)

!      DETERMINE GRID PARAMETERS.

       ALLOCATE(KLONRGG(NR))
       KLONRGG=0


       RMONOP = REAL(KSEC2(4)/1000)+0.1*(MOD(KSEC2(4),1000)/100)+
     &          0.01*(MOD(KSEC2(4),100)/10)+0.001*MOD(KSEC2(4),10)

       RMOWEP = REAL(KSEC2(5)/1000)+0.1*(MOD(KSEC2(5),1000)/100)+
     &          0.01*(MOD(KSEC2(5),100)/10)+0.001*MOD(KSEC2(5),10)

       RMOSOP = REAL(KSEC2(7)/1000)+0.1*(MOD(KSEC2(7),1000)/100)+
     &          0.01*(MOD(KSEC2(7),100)/10)+0.001*MOD(KSEC2(7),10)

!!! THERE IS A DANGER THAT THE DEFINITON FOR RMOEAP MIGHT VARY DUE TO
!!! THE AMBIGOUS DEFINITION FOR IRREGULAR GRIDS. FOR NON WAVE FIELDS,
!!! A GAUSSIAN GRID IMPLIES THAT THE GRID IS GLOBAL, THEREFORE
!!! RMOEAP IS IMPLICITLY KNOWN.
       CSTREAM='****'
       ILEVTYPE=KSEC1(7)
       JSTREAM=KSEC1(40)
       CALL WSTREAM_STRG(JSTREAM,CSTREAM,NENSFNB,NTOTENS,MARSFCTYPE,
     &                   KSTREAM,LASTREAM)

       IF(CSTREAM.EQ.'****' .OR.
     &   (LASTREAM .AND. ILEVTYPE.NE.209 .AND. ILEVTYPE.NE.212 )) THEN 
         LLNONWAVE=.TRUE.
       ELSE
         LLNONWAVE=.FALSE.
       ENDIF

       IF(IREPR.EQ.4 .AND. LLNONWAVE) THEN
         DELLO = 360./MAX(1,NC)
         RMOEAP = RMOWEP+360. - DELLO
         IPERIODIC = 1
       ELSE
         RMOEAP = REAL(KSEC2(8)/1000)+0.1*(MOD(KSEC2(8),1000)/100)+
     &           0.01*(MOD(KSEC2(8),100)/10)+0.001*MOD(KSEC2(8),10)
         IF(JRGG.EQ.1) CALL ADJUST (RMOWEP, RMOEAP)
         IPERIODIC = 0
         DELLO=(RMOEAP-RMOWEP)/MAX(1,NC-1)
         IF (RMOEAP-RMOWEP+1.5*DELLO.GE.360.) IPERIODIC = 1
       ENDIF

       IF(JRGG.EQ.1) THEN
          ISTART=0
         DO WHILE(KSEC2(23+ISTART).EQ.0)
           ISTART=ISTART+1
         ENDDO
         ISTOP=0
         DO WHILE(KSEC2(23+NRFULL-1-ISTOP).EQ.0)
           ISTOP=ISTOP+1
         ENDDO
         DO J=1,NR
           JSN=NR-J+1
           KLONRGG(JSN) = KSEC2(22+J+ISTART)
         ENDDO
         IDUM = KSEC2(4)-ISTART*KSEC2(10)
         RMONOP = REAL(IDUM/1000)+0.1*(MOD(IDUM,1000)/100)+
     &            0.01*(MOD(IDUM,100)/10)+0.001*MOD(IDUM,10)
         IDUM = KSEC2(7)+ISTOP*KSEC2(10)
         RMOSOP = REAL(IDUM/1000)+0.1*(MOD(IDUM,1000)/100)+
     &          0.01*(MOD(IDUM,100)/10)+0.001*MOD(IDUM,10)
        ELSEIF(JRGG.EQ.0) THEN
          KLONRGG=NC
        ELSE
          WRITE(IU06,*)
     &    ' REPRESENTATION OF THE FIELD NOT KNOWN'
        CALL ABORT1
       ENDIF

!      FIND WHETHER INTERPOLATION IS NEEDED

       LLINTERPOL=.TRUE.
       LLFR1OK=.TRUE.

       NGYFULL=INT(180./XDELLA)+1

       KAMONOP=NINT(AMONOP*100.)
       KRMONOP=NINT(RMONOP*100.)
       KAMOSOP=NINT(AMOSOP*100.)
       KRMOSOP=NINT(RMOSOP*100.)
       KAMOWEP=NINT(AMOWEP*100.)
       KRMOWEP=NINT(RMOWEP*100.)
       KAMOEAP=NINT(AMOEAP*100.)
       KRMOEAP=NINT(RMOEAP*100.)
       IF((KAMONOP.EQ.KRMONOP.OR.KRMONOP.EQ.9000).AND.
     &    (KAMOSOP.EQ.KRMOSOP.OR.KRMOSOP.EQ.-9000).AND.
     &    KAMOWEP.EQ.KRMOWEP.AND.KAMOEAP.EQ.KRMOEAP ) THEN
          IF(JRGG.EQ.IRGG.AND.
     &       (NR.EQ.NGY.OR.NR.EQ.NGYFULL)) THEN
             IF(IRGG.EQ.1) THEN
               LLINTERPOL=.FALSE.
               DO J=1,NGY
                 IF(KLONRGG(J).NE.NLONRGG(J)) THEN
                   LLINTERPOL=.TRUE.
                   EXIT
                 ENDIF
               ENDDO
             ELSE
               IF(NC.EQ.NGX) LLINTERPOL=.FALSE.
             ENDIF
          ENDIF
       ENDIF
       IF(LLINTERPOL) THEN
           WRITE(IU06,*) ' SPATIAL INTERPOLATION REQUIRED '
           WRITE(IU06,*) ' RMONOP = ',RMONOP,' -> AMONOP = ',AMONOP
           WRITE(IU06,*) ' RMOSOP = ',RMOSOP,' -> AMOSOP = ',AMOSOP
           WRITE(IU06,*) ' RMOWEP = ',RMOWEP,' -> AMOWEP = ',AMOWEP
           WRITE(IU06,*) ' RMOEAP = ',RMOEAP,' -> AMOEAP = ',AMOEAP
           WRITE(IU06,*) ' JRGG = ',JRGG,' -> IRGG = ',IRGG
           WRITE(IU06,*) ' NR = ',NR,' -> NGY = ',NGY
           WRITE(IU06,*) ' NC = ',NC,' -> NGX = ',NGX
       ENDIF

       IF(IPARAM.EQ.251) THEN
         KANG=KSEC1(46)
         KFRE=KSEC1(47)
!        input data first frequency
         KFR1=KSEC1(50+KANG)
!        model first frequency
         IFR1=NINT(KSEC1(49)*FR(1))
         IF(KANG.NE.NANG .OR. KFRE.NE.NFRE .OR. IFR1.NE.KFR1 ) THEN
!        it's assumed here that we will not change the initial direction
!        without changing the total number of directions
           LLINTERPOL=.TRUE.
         ENDIF

         IF(IFR1.NE.KFR1 ) THEN
           LLFR1OK=.FALSE.
           IF(IFR1.LT.KFR1) THEN
             IC=1
             IFR=IFR1
             DO WHILE (IFR.LT.KFR1 .AND. IC.LT.NFRE) 
               IC=IC+1
               IFR=NINT(KSEC1(49)*FR(IC))
             ENDDO
             MMSHIFT=IC-1
           ELSE
             IC=1
             IFR=KFR1
             DO WHILE (IFR.LT.IFR1 .AND. IC.LT.NFRE) 
               IC=IC+1
               IFR=KSEC1(49+KANG+IC)
             ENDDO
             MMSHIFT=-(IC-1)
           ENDIF
         ELSE
           LLFR1OK=.TRUE.
           MMSHIFT=0
         ENDIF


         open(17,FILE='wave_spectral_resolution')
         WRITE(17,*) 'NANG ',KANG
         WRITE(17,*) 'NFRE ',KFRE
         WRITE(17,*) 'INTERPOL ',LLINTERPOL
         close(17)
       ENDIF

       IF(.NOT.LLINTERPOL) THEN
!      THERE IS NO NEED FOR INTERPOLATION THE INPUT FILE
!      WILL SIMPLY BE RECOPIED TO THE OUTPUT DESTINATION.
         WRITE(IU06,*) ' THERE IS NO NEED FOR INTERPOLATION'
         WRITE(IU06,*) ' THE INPUT FILE WILL SIMPLY BE RECOPIED'
         WRITE(IU06,*) ' TO THE OUTPUT DESTINATION.'

         CALL SYSTEM('cp '//IFILENAME//' '//OFILENAME)

         GOTO 8888 
       ENDIF

       DEALLOCATE(PSEC4)


       IF(IPARAM.EQ.250) THEN
         WRITE(IU06,*) '***************************************'
         WRITE(IU06,*) '*                                     *'
         WRITE(IU06,*) '* THE FIELDS OF PARAMETER 250 MUST BE *'
         WRITE(IU06,*) '* CONVERTED TO PARAMETER 251 USING    *'
         WRITE(IU06,*) '* convert_grbspec                     *' 
         WRITE(IU06,*) '*                                     *'
         WRITE(IU06,*) '***************************************'
         CALL ABORT1

       ELSEIF(IPARAM.EQ.251) THEN

         KANG=KSEC1(46)
         KFRE=KSEC1(47)

         WRITE(IU06,*) ' '
         IF(KANG.EQ.NANG .AND. KFRE.EQ.NFRE .AND. LLFR1OK ) THEN
           WRITE(IU06,*) ' NO DIRECTIONAL OR FREQUENCY INTERPOLATION',
     &                   ' REQUIRED. '
         ELSEIF(KANG.EQ.NANG .AND. .NOT.LLFR1OK ) THEN
           WRITE(IU06,*) ' NO DIRECTIONAL INTERPOLATION REQUIRED'
           WRITE(IU06,*) ' BUT FREQUENCY INTERPOLATION IS REQUIRED'
           WRITE(IU06,*) ' BECAUSE OF DIFFERENT FIRST FREQUENCY'
           WRITE(IU06,*) ' INPUT FIRST FREQUENCY: ',KFR1
           WRITE(IU06,*) ' MODEL FIRST FREQUENCY: ',IFR1
         ELSEIF( KANG.EQ.2*NANG .OR. 2*KANG.EQ.NANG .OR.
     &           KANG.EQ.3*NANG .OR. 3*KANG.EQ.NANG .OR.
     &         2*KANG.EQ.3*NANG .OR. 3*KANG.EQ.2*NANG) THEN
           WRITE(IU06,*) ' DIRECTIONAL INTERPOLATION REQUIRED.'
           WRITE(IU06,*) ' INPUT NUMBER OF DIRECTIONS: ',KANG 
           WRITE(IU06,*) ' MODEL NUMBER OF DIRECTIONS: ',NANG 
           WRITE(IU06,*) ' '
           IF(KFRE.EQ.NFRE .AND. LLFR1OK ) THEN
             WRITE(IU06,*) ' BUT NO FREQUENCY INTERPOLATION REQUIRED.'
           ELSEIF(.NOT.LLFR1OK .OR. KFRE.NE.NFRE ) THEN
             WRITE(IU06,*) ' AND FREQUENCY INTERPOLATION ALSO REQUIRED'
             WRITE(IU06,*) ' INPUT NUMBER OF FREQUENCIES: ',KFRE
             WRITE(IU06,*) ' MODEL NUMBER OF FREQUENCIES: ',NFRE
             WRITE(IU06,*) ' INPUT FIRST FREQUENCY: ',KFR1
             WRITE(IU06,*) ' MODEL FIRST FREQUENCY: ',IFR1
           ELSE
             WRITE(IU06,*) '*******************************************'
             WRITE(IU06,*) ' FREQUENCY REQUIRED.'
             WRITE(IU06,*) ' BUT THE PROGRAM DOES NOT KNOW WHAT TO DO !'
             WRITE(IU06,*) ' IT ABORTS '
             WRITE(IU06,*) ' INPUT NUMBER OF FREQUENCIES: ',KFRE
             WRITE(IU06,*) ' MODEL NUMBER OF FREQUENCIES: ',NFRE
             WRITE(IU06,*) '*******************************************'
             CALL ABORT1
           ENDIF
         ELSEIF( KANG.EQ.NANG .AND. KFRE.NE.NFRE) THEN
           WRITE(IU06,*) ' NO DIRECTIONAL INTERPOLATION REQUIRED.'
           IF(.NOT.LLFR1OK .OR. KFRE.NE.NFRE ) THEN
             WRITE(IU06,*) ' BUT FREQUENCY INTERPOLATION IS REQUIRED.'
             WRITE(IU06,*) ' INPUT NUMBER OF FREQUENCIES: ',KFRE
             WRITE(IU06,*) ' MODEL NUMBER OF FREQUENCIES: ',NFRE
             WRITE(IU06,*) ' INPUT FIRST FREQUENCY: ',KFR1
             WRITE(IU06,*) ' MODEL FIRST FREQUENCY: ',IFR1
           ELSE
             WRITE(IU06,*) '*******************************************'
             WRITE(IU06,*) ' FREQUENCY REQUIRED.'
             WRITE(IU06,*) ' BUT THE PROGRAM DOES NOT KNOW WHAT TO DO !'
             WRITE(IU06,*) ' IT ABORTS '
             WRITE(IU06,*) ' INPUT NUMBER OF FREQUENCIES: ',KFRE
             WRITE(IU06,*) ' MODEL NUMBER OF FREQUENCIES: ',NFRE
             WRITE(IU06,*) '*******************************************'
             CALL ABORT1
           ENDIF
         ELSE
           WRITE(IU06,*) '*********************************************'
           WRITE(IU06,*) ' DIRECTIONAL REQUIRED.'
           WRITE(IU06,*) ' BUT THE PROGRAM DOES NOT KNOW HOW TO PROCEED'
           WRITE(IU06,*) ' IT ABORTS ' 
           WRITE(IU06,*) ' INPUT NUMBER OF DIRECTIONS: ',KANG 
           WRITE(IU06,*) ' MODEL NUMBER OF DIRECTIONS: ',NANG 
           WRITE(IU06,*) ' INPUT NUMBER OF FREQUENCIES: ',KFRE
           WRITE(IU06,*) ' MODEL NUMBER OF FREQUENCIES: ',NFRE
           WRITE(IU06,*) '*********************************************'
           CALL ABORT1
         ENDIF
         WRITE(IU06,*) ' '

         ALLOCATE(FIELD(NGX,NGY))
         ALLOCATE(TEMP(NGX,NGY))

         KKKPR=0

         IF(.NOT. LLFR1OK) THEN
!        SET MODEL VALUES TO MISSING FOR THE FIRST MMSHIFT FREQEUNCY BINS
!        (WHEN MMSHIFT>0) SINCE NO INPUT DATA AVAILABLE FOR THOSE BINS.
           FIELD=PSEK3(2)
           DO M=1,MMSHIFT
             DO K=1,NANG
               CALL OUTWAMINPUT(FIELD, KSEC1, K, M, OFILENAME)
             ENDDO
           ENDDO
         ENDIF

         DO IALLFLD=1,KANG*KFRE

!          GET THE DATA AND PERFORM INTERPOLATION TO WAVE MODEL GRID.

           CALL GRB2WGRD (IU06, ITEST, NKSEK0, NKSEK1, NKSEK2, NKSEK3,
     &                    NKSEK4, NPSEK2, NPSEK3, KSEK3, PSEK3,
     &                    KGRIB, ISIZE,
     &                    NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                    AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &                    PPREC, PPEPS,
     &                    CDATE, IFORP, JPARAM, KZLEV, KKK, MMM, FIELD)

           DEALLOCATE(KGRIB)

!          FREQUENCY INTERPOLATION
!          -----------------------
           M=MMM+MMSHIFT
           MLAST=M
           IF(M.LE.0 .OR. M.GT.NFRE) GOTO 2000

!          DIRECTION INTERPOLATION
!          -----------------------
           IF( KANG.EQ.NANG ) THEN
!            NO SPECTRAL INTERPOLATION REQUIRED
             K=KKK
             CALL OUTWAMINPUT(FIELD, KSEC1, K, M, OFILENAME)
           ELSEIF( 2*KANG.EQ.NANG) THEN
!            DOUBLING THE NUMBER OF DIRECTIONS
             K=2*KKK-1
             CALL OUTWAMINPUT(FIELD, KSEC1, K, M, OFILENAME)
             K=2*KKK
             CALL OUTWAMINPUT(FIELD, KSEC1, K, M, OFILENAME)

           ELSEIF( 3*KANG.EQ.NANG) THEN
!            TRIPLING THE NUMBER OF DIRECTIONS
             K=3*KKK-2
             CALL OUTWAMINPUT(FIELD, KSEC1, K, M, OFILENAME)
             K=3*KKK-1
             CALL OUTWAMINPUT(FIELD, KSEC1, K, M, OFILENAME)
             K=3*KKK
             CALL OUTWAMINPUT(FIELD, KSEC1, K, M, OFILENAME)

           ELSEIF( 3*KANG.EQ.2*NANG) THEN
!            3/2 OF THE NUMBER OF DIRECTIONS
             IF(MOD(KKK,2).EQ.1) THEN

!              KEEP INPUT (!! TO AVERAGE OUT FOR NEXT OUTPUT BIN)
               DO J=1,NGY
                 JSN=NGY-J+1
                 DO I=1,NLONRGG(JSN)
                   TEMP(I,J)=0.5*FIELD(I,J)
                 ENDDO
               ENDDO

               K=3*KKK/2
               CALL OUTWAMINPUT(FIELD, KSEC1, K, M, OFILENAME)

             ELSE
!              AVERAGE OUT THE INPUT WITH PREVIOUSLY SAVED ONE 
               DO J=1,NGY
                 JSN=NGY-J+1
                 DO I=1,NLONRGG(JSN)
                   IF(FIELD(I,J).NE.PSEK3(2)) THEN
                     IF(TEMP(I,J).NE.PSEK3(2)) THEN
                       TEMP(I,J)=TEMP(I,J)+0.5*FIELD(I,J)
                     ELSE
                       TEMP(I,J)=0.5*FIELD(I,J)
                     ENDIF
                   ELSE
                     IF(TEMP(I,J).NE.PSEK3(2)) THEN
                       TEMP(I,J)=0.5*TEMP(I,J)
                     ENDIF
                   ENDIF
                 ENDDO
               ENDDO
               K=3*KKK/2-1
               CALL OUTWAMINPUT(TEMP, KSEC1, K, M, OFILENAME)

               K=3*KKK/2
               CALL OUTWAMINPUT(FIELD, KSEC1, K, M, OFILENAME)
             ENDIF

           ELSEIF( KANG.EQ.2*NANG) THEN
!            HALVING THE NUMBER OF DIRECTIONS
             IF(MOD(KKK,2).EQ.1) THEN
!              SAVE DIRECTIONAL CONTRIBUTION TO BE COMBINED
!              WITH THE NEXT ONE.
               DO J=1,NGY
                JSN=NGY-J+1
                DO I=1,NLONRGG(JSN)
                  IF(FIELD(I,J).NE.PSEK3(2)) THEN
                    TEMP(I,J)=FIELD(I,J)
                  ELSE
                    TEMP(I,J)=0. 
                  ENDIF
                ENDDO
               ENDDO

               KKKPR=KKK
             ELSEIF(KKK.EQ.KKKPR+1) THEN
               DO J=1,NGY
                JSN=NGY-J+1
                DO I=1,NLONRGG(JSN)
                  IF(FIELD(I,J).NE.PSEK3(2)) THEN
                    FIELD(I,J)=0.5*(FIELD(I,J)+TEMP(I,J))
                  ELSE
                    IF(TEMP(I,J).NE.0.) THEN
                      FIELD(I,J)=0.5*TEMP(I,J) 
                    ENDIF
                  ENDIF
                ENDDO
               ENDDO

               K=KKK/2
               CALL OUTWAMINPUT(FIELD, KSEC1, K, M, OFILENAME)
             ELSE
               WRITE(IU06,*) '************************************'
               WRITE(IU06,*) ' FOR THIS TYPE OF'
               WRITE(IU06,*) ' DIRECTIONAL INTERPOLATION,'
               WRITE(IU06,*) ' THE DIRECTIONS IN THE INPUT FIELD ' 
               WRITE(IU06,*) ' MUST BE CONCECUTIVE, STARTING WITH'
               WRITE(IU06,*) ' THE FIRST ONE !!!!' 
               WRITE(IU06,*) ' THE PROGRAM ABORTS ' 
               WRITE(IU06,*) '*************************************'
               CALL ABORT1
             ENDIF
           ELSEIF( KANG.EQ.3*NANG) THEN
!            1/3 OF THE NUMBER OF DIRECTIONS
             IF(MOD(KKK,3).EQ.1) THEN
               DO J=1,NGY
                 JSN=NGY-J+1
                 DO I=1,NLONRGG(JSN)
                   TEMP(I,J)=0.
                 ENDDO
               ENDDO
             ENDIF
             IF(MOD(KKK,3).NE.0) THEN
!              SAVE DIRECTIONAL CONTRIBUTION TO BE COMBINED
!              WITH THE NEXT ONE.
               DO J=1,NGY
                JSN=NGY-J+1
                DO I=1,NLONRGG(JSN)
                  IF(FIELD(I,J).NE.PSEK3(2)) THEN
                    TEMP(I,J)=TEMP(I,J)+ONETHIRD*FIELD(I,J)
                  ENDIF
                ENDDO
               ENDDO

             ELSE
               DO J=1,NGY
                JSN=NGY-J+1
                DO I=1,NLONRGG(JSN)
                  IF(FIELD(I,J).NE.PSEK3(2)) THEN
                    FIELD(I,J)=ONETHIRD*FIELD(I,J)+TEMP(I,J)
                  ELSE
                    IF(TEMP(I,J).NE.0.) THEN
                      FIELD(I,J)=TEMP(I,J)
                    ENDIF
                  ENDIF
                ENDDO
               ENDDO

               K=KKK/3
               CALL OUTWAMINPUT(FIELD, KSEC1, K, M, OFILENAME)
             ENDIF

           ELSEIF( 2*KANG.EQ.3*NANG) THEN
!            2/3 OF THE NUMBER OF DIRECTIONS
             IF(MOD(KKK,3).EQ.1) THEN
!              KEEP INPUT
               DO J=1,NGY
                 JSN=NGY-J+1
                 DO I=1,NLONRGG(JSN)
                   TEMP(I,J)=FIELD(I,J)
                 ENDDO
               ENDDO

             ELSEIF(MOD(KKK,3).EQ.2) THEN
!              ADD CONTRIBUTION OF PREVIOUSLY SAVED TO INPUT 
               DO J=1,NGY
                 JSN=NGY-J+1
                 DO I=1,NLONRGG(JSN)
                   IF(FIELD(I,J).NE.PSEK3(2)) THEN
                     IF(TEMP(I,J).NE.PSEK3(2)) THEN
                       TEMP(I,J)=TWOTHIRD*TEMP(I,J)+ONETHIRD*FIELD(I,J)
                     ELSE
                       TEMP(I,J)=ONETHIRD*FIELD(I,J)
                     ENDIF
                   ELSE
                     TEMP(I,J)=TWOTHIRD*TEMP(I,J)
                   ENDIF
                 ENDDO
               ENDDO
               K=2*(KKK/3)+1
               CALL OUTWAMINPUT(TEMP, KSEC1, K, M, OFILENAME)

!              KEEP INPUT
               DO J=1,NGY
                 JSN=NGY-J+1
                 DO I=1,NLONRGG(JSN)
                   TEMP(I,J)=FIELD(I,J)
                 ENDDO
               ENDDO

             ELSE
!              ADD INPUT TO PREVIOUSLY SAVED
               DO J=1,NGY
                 JSN=NGY-J+1
                 DO I=1,NLONRGG(JSN)
                   IF(FIELD(I,J).NE.PSEK3(2)) THEN
                     IF(TEMP(I,J).NE.PSEK3(2)) THEN
                       TEMP(I,J)=ONETHIRD*TEMP(I,J)+TWOTHIRD*FIELD(I,J)
                     ELSE
                       TEMP(I,J)=TWOTHIRD*FIELD(I,J)
                     ENDIF
                   ELSE
                     IF(TEMP(I,J).NE.PSEK3(2)) THEN
                       TEMP(I,J)=ONETHIRD*TEMP(I,J)
                     ENDIF
                   ENDIF
                 ENDDO
               ENDDO
               K=2*(KKK/3)
               CALL OUTWAMINPUT(TEMP, KSEC1, K, M, OFILENAME)

             ENDIF

           ELSE
             WRITE(IU06,*) '************************************'
             WRITE(IU06,*) ' NOT READY FOR THIS TYPE OF'
             WRITE(IU06,*) ' DIRECTIONAL INTERPOLATION,'
             WRITE(IU06,*) ' INPUT NUMBER OF DIRECTIONS: ',KANG 
             WRITE(IU06,*) ' MODEL NUMBER OF DIRECTIONS: ',NANG 
             WRITE(IU06,*) '*************************************'
             CALL ABORT1
           ENDIF



2000       CONTINUE
!          CONTINUE READING AND DECODING UNTIL LAST FIELD IS READ.
!          -------------------------------------------------------
           IF(IALLFLD.EQ.KANG*KFRE) EXIT

2001       KPLENG=ISIZE*KIND(ISIZE)
           IF(.NOT.ALLOCATED(KGRIB)) ALLOCATE(KGRIB(ISIZE))

           CALL PBGRIB(IUIN,KGRIB,KPLENG,KLEN,KRET)
           IF(KRET.EQ.-1) THEN
             IF(IALLFLD.EQ.1) THEN
               WRITE(IU06,*) ' ONLY ONE SPECTRAL FIELD HAS BEEN READ !'
               EXIT
             ELSE
               WRITE(IU06,*) '**********************************'
               WRITE(IU06,*) '*  PBGRIB:  END OF FILE ENCOUNTED'
               WRITE(IU06,*) '**********************************'
               CALL ABORT1
             ENDIF
           ENDIF
           IF(KRET.EQ.-2) THEN
             WRITE(IU06,*) '**********************************'
             WRITE(IU06,*) '*  PBGRIB: FILE HANDLING ERROR'
             WRITE(IU06,*) '**********************************'
             CALL ABORT1
           ENDIF
           IF(KRET.EQ.-3) THEN
             CALL KGRIBSIZE(IU06, IUIN, KLEN, ISIZE, 'INTWAMINPUT')
             DEALLOCATE(KGRIB)
             GOTO 2001
           ENDIF

         ENDDO

!        SET MODEL VALUES TO MISSING IF NO VALUES WERE PROVIDED 
!        BY INPUT DATA 
         FIELD=PSEK3(2)
         DO M=MLAST+1,NFRE
           DO K=1,NANG
             CALL OUTWAMINPUT(FIELD, KSEC1, K, M, OFILENAME)
           ENDDO
         ENDDO

       ELSE
!      ALL THE OTHER WAVE PARAMETERS

!        GET THE DATA AND PERFORM INTERPOLATION TO WAVE MODEL GRID.

         ALLOCATE(FIELD(NGX,NGY))
         CALL GRB2WGRD (IU06, ITEST, NKSEK0, NKSEK1, NKSEK2, NKSEK3,
     &                  NKSEK4, NPSEK2, NPSEK3, KSEK3, PSEK3,
     &                  KGRIB, ISIZE,
     &                  NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                  AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &                  PPREC, PPEPS,
     &                  CDATE, IFORP, JPARAM, KZLEV, KKK, MMM, FIELD)

         DEALLOCATE(KGRIB)


!        ENCODE AND OUTPUT INTEGRATED PARAMETER
!        --------------------------------------

         KKK=0
         MMM=0

         CALL OUTWAMINPUT(FIELD, KSEC1, KKK, MMM, OFILENAME)

       ENDIF

8888   WRITE (IU06,*) ' PROGRAM INTWAMINPUT: ALL DONE'

      END
