      PROGRAM intwaminput

! ----------------------------------------------------------------------

!**** *INTWAMINPUT* -

!     J. BIDLOT    ECMWF  JULY 2000 
!     J. BIDLOT    ECMWF  MARCH 2010 : adapt to use grib_api 

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
      USE YOWGRIB_HANDLES , ONLY :NGRIB_HANDLE_WAM_I,NGRIB_HANDLE_WAM_S
      USE YOWGRIBHD, ONLY : PPEPS    ,PPREC    ,NGRBRESI ,PPMIN_RESET,
     &            NGRBRESS ,HOPERI   ,HOPERS   ,LGRHDIFS ,LNEWLVTP ,
     &            LPADPOLES
      USE YOWGRID  , ONLY : DELPHI   ,NLONRGG  ,IJL
      USE YOWMAP   , ONLY : IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO   ,IQGAUSS
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP     ,
     &            KTAG     ,NPRECR   ,NPRECI
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY      ,
     &            NBLO     ,NIBLO    ,CLDOMAIN
      USE YOWPCONS , ONLY : ZMISS    ,EPSMIN
      USE YOWSTAT  , ONLY : MARSTYPE ,YCLASS   ,YEXPVER  ,ISTREAM  ,
     &            NLOCGRB , IREFRA   ,NENSFNB  ,NTOTENS  ,NPROMA_WAM
      USE YOWSPEC  , ONLY : NSTART   ,NEND
      USE YOWUNPOOL ,ONLY : LLUNSTR
      USE YOWWIND  , ONLY : NXFF     ,NYFF
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE MPL_MODULE,ONLY : MPL_INIT, MPL_END

      USE GRIB_API_INTERFACE


! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: NBIT = 1600000

      INTEGER :: I4(2), KRET, KPLENG, ISIZE, KLEN
      INTEGER :: ISYSTEM_STAT, ISTAT
      INTEGER :: KLENG, KLENP, KWORD, NC, NR, NTOT, I, J, JRGG, IR
      INTEGER :: IREPR, IPARAM, ITABLE, KANG, KFRE, IALLFLD, JSN
      INTEGER :: IFR1, KFR1, IFR, IC, JC, MMSHIFT, MLAST
      INTEGER :: IFORP, KZLEV, K, M, KKK, MMM, JPARAM, KKKPR, NGYFULL
      INTEGER :: IPERIODIC, NRFULL, JSTREAM, ISTART, ISTOP, IDUM
      INTEGER :: IU07, ILEVTYPE, KSTREAM
      INTEGER :: KAMONOP, KAMOEAP, KAMOSOP, KAMOWEP
      INTEGER :: KRMONOP, KRMOEAP, KRMOSOP, KRMOWEP
      INTEGER :: I_GET_UNIT
      INTEGER :: LFILE, KFILE_HANDLE, KGRIB_HANDLE, IGRIB_LEN
      INTEGER :: IYYYYMMDD, IHHMM
      INTEGER :: IVAL
      INTEGER :: IRET
      INTEGER :: IPLPRESENT, NB_PL
      INTEGER :: IDIRSCALING, IFRESCALING
      INTEGER :: IUOUT, KFDB
      INTEGER, DIMENSION(:), ALLOCATABLE :: PL
      INTEGER, ALLOCATABLE :: KGRIB_BUFR(:)
      INTEGER, ALLOCATABLE ::  NLONRGG_LOC(:)
      INTEGER, ALLOCATABLE ::  KLONRGG(:)
      INTEGER(KIND=JPKSIZE_T) :: KBYTES

      REAL :: X4(2)
      REAL :: STEP, START_STEP, END_STEP
      REAL :: RMONOP, RMOSOP, RMOEAP, RMOWEP
      REAL :: PRPLRADI
      REAL ::  DELLO, DELLA
      REAL :: ONETHIRD, TWOTHIRD
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FIELD, TEMP
      REAL, ALLOCATABLE, DIMENSION(:,:) :: XLON, YLAT
      REAL, ALLOCATABLE :: SCFR(:), SCTH(:)


      CHARACTER(LEN=1) :: CFDBSF
      CHARACTER(LEN=2) :: MARSFCTYPE
      CHARACTER(LEN=4) :: CSTREAM
      CHARACTER(LEN=8) :: CSTEPTYPE
      CHARACTER(LEN=14) :: CDATE
      CHARACTER(LEN=11) :: IFILENAME
      CHARACTER(LEN=12) :: OFILENAME
      CHARACTER(LEN=12) :: CGRIDTYPE

      LOGICAL :: LLINTERPOL, LLNONWAVE, LASTREAM
      LOGICAL :: LLFR1OK
      LOGICAL :: LLEXIST
      LOGICAL :: LFDB, LFDBOPEN
      LOGICAL :: LLWAIT

! ----------------------------------------------------------------------

      CALL MPL_INIT()
      IU06 = 6
      IU07 = I_GET_UNIT(IU06, 'wam_grid_tables', 'r', 'u', 0)

      WRITE(IU06,*) ' INTWAMINPUT STARTED '

!!!!!  because this program will write its output in grib, it is not yet ready for unsctructured grid
      LLUNSTR=.FALSE.

      PRPLRADI=1.0
      CALL INIWCST(PRPLRADI)

!     1.1 INITIALISATION OF VARIABLES WITH DEFAULT VALUES
!         ---------------------------------------------------

      ITEST = 1 
      LMESSPASS = .TRUE.
      MARSTYPE = 'an'
      YCLASS   = 'od'
      YEXPVER  = '0001' 
      NENSFNB = 0
      NTOTENS = 0
      ISTREAM = 1045 !!! is changed to an ifs stream also change LNEWLVTP 
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
      ! reset PPMIN to avoid imposing a minmimum value in *WGRIBENCODE*
      PPMIN_RESET=LOG10(PPEPS)+ABS(PPREC)+EPSMIN 
      LPADPOLES=.FALSE. ! do not pad poles in *WGRIBENCODE*

      KFDB=0
      LFDB=.FALSE.
      LFDBOPEN=.FALSE.
      ALLOCATE (NSTART(NPROC),NEND(NPROC))

      NPROMA_WAM=1
      NPROMA_WAM=HUGE(NPROMA_WAM)/2

      KGRIB_HANDLE=0
! ----------------------------------------------------------------------

!*    2. READ PREPROC OUTPUT.
!        --------------------

      CALL READPRE (IU07)

      NINF=1
      NSUP=NIBLO

      NSTART=1
      NEND=IJL(1)

      IF(LLUNSTR) THEN
        WRITE(*,*) 'NOT YET READY FOR UNSTRUCTURED GRID '
        CALL ABORT1
      ELSE
        NXFF=NGX
        NYFF=NGY
        ALLOCATE(NLONRGG_LOC(NYFF))
      ENDIF

      IF(.NOT.ALLOCATED(XLON)) ALLOCATE(XLON(NXFF,NYFF))
      XLON(:,:)=ZMISS
      IF(.NOT.ALLOCATED(YLAT)) ALLOCATE(YLAT(NXFF,NYFF))
      YLAT(:,:)=ZMISS

      IF(LLUNSTR) THEN
        WRITE(*,*) 'NOT YET READY FOR UNSTRUCTURED GRID '
!!!!    still too many NGX, NGY
        CALL ABORT1
      ELSE
        NLONRGG_LOC(:)=NLONRGG(:)
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(J,I,JSN)
        DO J=1,NYFF
          JSN=NYFF-J+1
          DO I=1,NLONRGG(JSN)
            XLON(I,J) = AMOWEP + (I-1)*ZDELLO(JSN)
            YLAT(I,J) = AMOSOP + (JSN-1)*XDELLA
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO
      ENDIF


      IFILENAME='input_field'
      OFILENAME='output_field'
  
!     CONNECT TO INPUT FILE
      LFILE=0
      LLEXIST=.FALSE.
      IF (IFILENAME.NE. ' ') LFILE=LEN_TRIM(IFILENAME)
      INQUIRE(FILE=IFILENAME(1:LFILE),EXIST=LLEXIST)
      IF(LLEXIST) THEN
        CALL IGRIB_OPEN_FILE(KFILE_HANDLE,IFILENAME(1:LFILE),'r')
      ELSE
        WRITE(*,*)'****************************'
        WRITE(*,*)'*                          *'
        WRITE(*,*)'*GRIB DATA NOT FOUND IN *'
        WRITE(*,*)  IFILENAME 
        WRITE(*,*)'*PROGRAM WILL ABORT        *'
        WRITE(*,*)'*                          *'
        WRITE(*,*)'****************************'
        CALL ABORT1
      ENDIF

!     LOAD THE DATA
1021  ISIZE=NBIT
      KBYTES=ISIZE*NPRECI
      ALLOCATE(KGRIB_BUFR(ISIZE))

      CALL IGRIB_READ_FROM_FILE(KFILE_HANDLE,KGRIB_BUFR,KBYTES,IRET)

      IF(IRET.EQ.JPGRIB_BUFFER_TOO_SMALL) THEN
!!!     *IGRIB_READ_FROM_FILE* does not read through the file if
!!!     the size is too small, so figure out the size and read again.
        CALL KGRIBSIZE(IU06, KBYTES, NBIT, 'INTWAMINPUT')
        DEALLOCATE(KGRIB_BUFR)
        GOTO 1021
      ELSEIF(IRET.EQ.JPGRIB_END_OF_FILE) THEN
        WRITE(IU06,*) '*************************************'
        WRITE(IU06,*) '* INTWAMINPUT: END OF FILE ENCOUNTED'
        WRITE(IU06,*) '*************************************'
        CALL ABORT1
      ELSEIF(IRET.NE.JPGRIB_SUCCESS) THEN
        WRITE(IU06,*) '*************************************'
        WRITE(IU06,*) '* INTWAMINPUT: FILE HANDLING ERROR'
        WRITE(IU06,*) '*************************************'
        CALL ABORT1
      ENDIF

      KGRIB_HANDLE=-99
      CALL IGRIB_NEW_FROM_MESSAGE(KGRIB_HANDLE,KGRIB_BUFR)

      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'Nj',NR)

      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'gridType', CGRIDTYPE)
      IF(CGRIDTYPE(1:10) == 'regular_gg') THEN
        JRGG=1
        IREPR=4
      ELSEIF(CGRIDTYPE(1:10) == 'reduced_gg') THEN
        JRGG=1
        IREPR=4
      ELSEIF(CGRIDTYPE(1:7) == 'regular') THEN
        JRGG=0
        IREPR=0
      ELSEIF(CGRIDTYPE(1:7) == 'reduced') THEN
        JRGG=1
        IREPR=0
      ELSE
        WRITE(IU06,*) '***********************************'
        WRITE(IU06,*) '*  GRID TYPE NOT RECOGNIZED !!!'
        WRITE(IU06,*) '   gridType = ', CGRIDTYPE 
        WRITE(IU06,*) '***********************************'
        CALL ABORT1
      ENDIF

      IF(JRGG.EQ.1) THEN
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'PLPresent',IPLPRESENT)
        IF (IPLPRESENT == 1) THEN
          CALL IGRIB_GET_VALUE(KGRIB_HANDLE,
     &                         'numberOfPointsAlongAMeridian',NB_PL)
          ALLOCATE(PL(NB_PL))
          CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'pl',PL)
        ELSE
          WRITE(IU06,*) 'NUMBER OF POINTS PER LATITUDE MISSING !!!'
          CALL ABORT1
        ENDIF
        NC=0
        DO J=1,NB_PL
          NC = MAX(NC,PL(J))
        ENDDO
        IR=0
        DO J=1,NB_PL
          IF(PL(J).NE.0) IR=IR+1
        ENDDO
        NR=IR
      ELSEIF(JRGG.EQ.0) THEN
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'Ni',IVAL)
        NC=IVAL
      ELSE
        WRITE(IU06,*)
     &  '   STRUCTURE OF THE FIELD NOT KNOWN'
        CALL ABORT1
      ENDIF

!*    DETERMINE CODE FOR DATA FIELD TYPE.
      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'paramId',IVAL)
      ITABLE=IVAL/1000
      IPARAM=IVAL-ITABLE*1000

!     DATE. 
      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'dataDate',IYYYYMMDD)
      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'time',IHHMM)
      WRITE(CDATE(1:12),'(I8.8,I4.4)') IYYYYMMDD,IHHMM 
      CDATE(13:14)='00'

      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'stepType',CSTEPTYPE)
!     FORECAST STEP (in seconds)
      IF(CSTEPTYPE(1:7) .EQ. 'instant') THEN
        CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stepUnits','s')
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'step',STEP)
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'startStep',START_STEP)
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'endStep',END_STEP)
!       THE DATA ARE VALID BETWEEN TWO TIMES. TAKE THE MIDDLE POINT
        IF(START_STEP.NE.END_STEP) THEN
          STEP=(END_STEP-START_STEP)/2
        ENDIF
        IFORP=STEP
      ELSE
        WRITE(*,*) 'UNKNOWN DEFINITION OF FORECAST STEP TYPE !!!'
        WRITE(*,*) 'stepType = ',CSTEPTYPE
        CALL ABORT1
      ENDIF


!     DETERMINE GRID PARAMETERS.

      IF(.NOT.ALLOCATED(KLONRGG)) ALLOCATE(KLONRGG(NR))
      KLONRGG=0

      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,
     &                     'latitudeOfFirstGridPointInDegrees',RMONOP)
      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,
     &                     'latitudeOfLastGridPointInDegrees',RMOSOP)

      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,
     &                     'longitudeOfFirstGridPointInDegrees',RMOWEP)

!!!   THERE IS A DANGER THAT THE DEFINITON FOR RMOEAP MIGHT VARY DUE TO
!!!   THE AMBIGOUS DEFINITION FOR IRREGULAR GRIDS. FOR NON WAVE FIELDS,
!!!   A GAUSSIAN GRID IMPLIES THAT THE GRID IS GLOBAL, THEREFORE
!!!   RMOEAP IS IMPLICITLY KNOWN.
      CSTREAM='****'
      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'levtype',ILEVTYPE)
      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'stream',JSTREAM)
      CALL WSTREAM_STRG(JSTREAM,CSTREAM,NENSFNB,NTOTENS,MARSFCTYPE,
     &                  KSTREAM,LASTREAM)

      IF(CSTREAM.EQ.'****' .OR.
     &  (LASTREAM .AND. ILEVTYPE.NE.209 .AND. ILEVTYPE.NE.212 )) THEN 
        LLNONWAVE=.TRUE.
      ELSE
        LLNONWAVE=.FALSE.
      ENDIF

      IF(IREPR.EQ.4 .AND. LLNONWAVE) THEN
        DELLO = 360./MAX(1,NC)
        RMOEAP = RMOWEP+360. - DELLO
        IPERIODIC = 1
      ELSE
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,
     &                      'longitudeOfLastGridPointInDegrees',RMOEAP)
        IF(JRGG.EQ.1) CALL ADJUST (RMOWEP, RMOEAP)
        IPERIODIC = 0
        DELLO=(RMOEAP-RMOWEP)/MAX(1,NC-1)
        IF (RMOEAP-RMOWEP+1.5*DELLO.GE.360.) IPERIODIC = 1
      ENDIF

      IF(JRGG.EQ.1) THEN
        ISTART=1
        DO WHILE(PL(ISTART).EQ.0 .AND. ISTART.LT.NB_PL)
          ISTART=ISTART+1
        ENDDO
        ISTART=ISTART-1

        ISTOP=0
        DO WHILE(PL(NB_PL-ISTOP).EQ.0 .AND. ISTOP.LT.NB_PL)
          ISTOP=ISTOP+1
        ENDDO

        DO J=1,NR-ISTART
          JSN=NR-J+1
          KLONRGG(JSN) = PL(J+ISTART)
        ENDDO
        DEALLOCATE(PL)

        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,
     &                       'jDirectionIncrementInDegrees',DELLA)
        RMONOP =  RMONOP-ISTART*DELLA
        RMOSOP =  RMOSOP+ISTOP*DELLA

      ELSEIF(JRGG.EQ.0) THEN
        KLONRGG=NC
      ELSE
        WRITE(IU06,*) ' REPRESENTATION OF THE FIELD NOT KNOWN'
        CALL ABORT1
      ENDIF

      IF(IPARAM.EQ.251) THEN
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'numberOfDirections',KANG)
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'numberOfFrequencies',KFRE)

!       input data first frequency
        IF(.NOT.ALLOCATED(SCFR)) ALLOCATE(SCFR(KFRE))
        IF(.NOT.ALLOCATED(SCTH)) ALLOCATE(SCTH(KANG))

        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,
     &                       'directionScalingFactor',IDIRSCALING)
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'scaledDirections',SCTH)

        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,
     &                      'frequencyScalingFactor',IFRESCALING)
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'scaledFrequencies',SCFR)

        KFR1=NINT(SCFR(1))
!       model first frequency
        IFR1=NINT(IFRESCALING*FR(1))

        IF(IFR1.NE.KFR1 ) THEN
          LLFR1OK=.FALSE.
          IF(IFR1.LT.KFR1) THEN
            IC=1
            IFR=IFR1
            DO WHILE (IFR.LT.KFR1 .AND. IC.LT.NFRE) 
              IC=IC+1
              IFR=NINT(IFRESCALING*FR(IC))
            ENDDO
            MMSHIFT=IC-1
          ELSE
            IC=1
            IFR=KFR1
            DO WHILE (IFR.LT.IFR1 .AND. IC.LT.NFRE) 
              IC=IC+1
              IFR=NINT(SCFR(IC))
            ENDDO
            MMSHIFT=-(IC-1)
          ENDIF
        ELSE
          LLFR1OK=.TRUE.
          MMSHIFT=0
        ENDIF

      ENDIF

!      FIND WHETHER INTERPOLATION IS NEEDED
!      ------------------------------------

       LLINTERPOL=.TRUE.

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
           WRITE(IU06,111) RMONOP,AMONOP
           WRITE(IU06,*) ' RMOSOP = ',RMOSOP,' -> AMOSOP = ',AMOSOP
           WRITE(IU06,111) RMOSOP,AMOSOP
           WRITE(IU06,*) ' RMOWEP = ',RMOWEP,' -> AMOWEP = ',AMOWEP
           WRITE(IU06,111) RMOWEP,AMOWEP
           WRITE(IU06,*) ' RMOEAP = ',RMOEAP,' -> AMOEAP = ',AMOEAP
           WRITE(IU06,111) RMOEAP,AMOEAP
           WRITE(IU06,*) ' JRGG = ',JRGG,' -> IRGG = ',IRGG
           WRITE(IU06,*) ' NR = ',NR,' -> NGY = ',NGY
           WRITE(IU06,*) ' NC = ',NC,' -> NGX = ',NGX
       ENDIF
111     FORMAT(4x,'HEX: ',3(Z16.16,2x))

       IF(IPARAM.EQ.251) THEN
         IF(KANG.NE.NANG .OR. KFRE.NE.NFRE .OR. IFR1.NE.KFR1 ) THEN
!        it's assumed here that we will not change the initial direction
!        without changing the total number of directions
           LLINTERPOL=.TRUE.
         ENDIF

         open(17,FILE='wave_spectral_resolution')
         WRITE(17,*) 'NANG ',KANG
         WRITE(17,*) 'NFRE ',KFRE
         WRITE(17,*) 'INTERPOL ',LLINTERPOL
         close(17)
       ENDIF

       IF(.NOT.LLINTERPOL) THEN
!      THERE IS NO NEED FOR INTERPOLATION, THE INPUT FILE
!      WILL SIMPLY BE RECOPIED TO THE OUTPUT DESTINATION.
         WRITE(IU06,*) ' THERE IS NO NEED FOR INTERPOLATION'
         WRITE(IU06,*) ' THE INPUT FILE WILL SIMPLY BE RECOPIED'
         WRITE(IU06,*) ' TO THE OUTPUT DESTINATION.'

         LLWAIT=.TRUE.
         ISYSTEM_STAT=-1
         ISTAT=-1
#if defined(__PGI)
! Not supported by the PGI compiler yet
#else
         CALL EXECUTE_COMMAND_LINE('cp '//IFILENAME//' '//OFILENAME,
     &        WAIT=LLWAIT, EXITSTAT=ISYSTEM_STAT,CMDSTAT=ISTAT)
#endif
         IF(ISTAT.NE.0) THEN
            WRITE(IU06,*) '*************************************'
            WRITE(IU06,*) '* INTWAMINPUT: FILE COPYING ERROR !'
            WRITE(IU06,*) '* ISTAT = ',ISTAT
            WRITE(IU06,*) '* ISYSTEM_STAT = ',ISYSTEM_STAT
            WRITE(IU06,*) '*************************************'
            CALL ABORT1
         ENDIF

         GOTO 8888 
       ELSE
!        PREPARE OUTPUT
         CALL IGRIB_GET_VALUE(KGRIB_HANDLE,
     &                        'localDefinitionNumber',NLOCGRB)
!        GRIB HANDLES FOR OUTPUT
         LGRHDIFS=.FALSE.
!        FOR INTEGRATED PARAMETERS
         CALL PRESET_WGRIB_TEMPLATE("I",NGRIB_HANDLE_WAM_I)
!        FOR SPECTRA
         CALL PRESET_WGRIB_TEMPLATE("S",NGRIB_HANDLE_WAM_S)

         CALL IGRIB_OPEN_FILE(IUOUT,OFILENAME,'w')
       ENDIF

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

         IF(.NOT.ALLOCATED(FIELD)) ALLOCATE(FIELD(NGX,NGY))
         IF(.NOT.ALLOCATED(TEMP)) ALLOCATE(TEMP(NGX,NGY))

         KKKPR=0

         IF(.NOT. LLFR1OK) THEN
!        SET MODEL VALUES TO MISSING FOR THE FIRST MMSHIFT FREQUENCY BINS
!        (WHEN MMSHIFT>0) SINCE NO INPUT DATA AVAILABLE FOR THOSE BINS.
           FIELD=ZMISS
           
           DO M=1,MMSHIFT
             DO K=1,NANG
               CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,
     &                       ITABLE, IPARAM, 0, K , M,
     &                       CDATE, IFORP, MARSTYPE,
     &                       LFDB, CFDBSF, KFDB, LFDBOPEN, IUOUT)
             ENDDO
           ENDDO
         ENDIF

         DO IALLFLD=1,KANG*KFRE

!          GET THE DATA AND PERFORM INTERPOLATION TO WAVE MODEL GRID.

           CALL GRIB2WGRID (IU06, ITEST, NPROMA_WAM, 
     &                      KGRIB_HANDLE, KGRIB_BUFR, ISIZE,
     &                      LLUNSTR,
     &                      NXFF, NYFF, NLONRGG_LOC,
     &                      IRGG, XDELLA, ZDELLO,
     &                      XLON, YLAT,
     &                      ZMISS, PPREC, PPEPS,
     &                      CDATE, IFORP, JPARAM, KZLEV,KKK,MMM,FIELD)


           DEALLOCATE(KGRIB_BUFR)

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
             CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,
     &                     ITABLE, IPARAM, 0, K , M,
     &                     CDATE, IFORP, MARSTYPE,
     &                     LFDB, CFDBSF, KFDB, LFDBOPEN, IUOUT)
           ELSEIF( 2*KANG.EQ.NANG) THEN
!            DOUBLING THE NUMBER OF DIRECTIONS
             K=2*KKK-1
             CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,
     &                     ITABLE, IPARAM, 0, K , M,
     &                     CDATE, IFORP, MARSTYPE,
     &                     LFDB, CFDBSF, KFDB, LFDBOPEN, IUOUT)
             K=2*KKK
             CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,
     &                     ITABLE, IPARAM, 0, K , M,
     &                     CDATE, IFORP, MARSTYPE,
     &                     LFDB, CFDBSF, KFDB, LFDBOPEN, IUOUT)
           ELSEIF( 3*KANG.EQ.NANG) THEN
!            TRIPLING THE NUMBER OF DIRECTIONS
             K=3*KKK-2
             CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,
     &                     ITABLE, IPARAM, 0, K , M,
     &                     CDATE, IFORP, MARSTYPE,
     &                     LFDB, CFDBSF, KFDB, LFDBOPEN, IUOUT)
             K=3*KKK-1
             CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,
     &                     ITABLE, IPARAM, 0, K , M,
     &                     CDATE, IFORP, MARSTYPE,
     &                     LFDB, CFDBSF, KFDB, LFDBOPEN, IUOUT)
             K=3*KKK
             CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,
     &                     ITABLE, IPARAM, 0, K , M,
     &                     CDATE, IFORP, MARSTYPE,
     &                     LFDB, CFDBSF, KFDB, LFDBOPEN, IUOUT)
           ELSEIF( 3*KANG.EQ.2*NANG) THEN
!            3/2 OF THE NUMBER OF DIRECTIONS
             IF(MOD(KKK,2).EQ.1) THEN

!              KEEP INPUT (!! TO AVERAGE OUT FOR NEXT OUTPUT BIN)
               DO J=1,NGY
                 JSN=NGY-J+1
                 DO I=1,NLONRGG(JSN)
                   IF(FIELD(I,J).NE.ZMISS) THEN
                     TEMP(I,J)=0.5*FIELD(I,J)
                   ENDIF
                 ENDDO
               ENDDO

               K=3*KKK/2
               CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,
     &                       ITABLE, IPARAM, 0, K , M,
     &                       CDATE, IFORP, MARSTYPE,
     &                       LFDB, CFDBSF, KFDB, LFDBOPEN, IUOUT)
             ELSE
!              AVERAGE OUT THE INPUT WITH PREVIOUSLY SAVED ONE 
               DO J=1,NGY
                 JSN=NGY-J+1
                 DO I=1,NLONRGG(JSN)
                   IF(FIELD(I,J).NE.ZMISS) THEN
                     IF(TEMP(I,J).NE.ZMISS) THEN
                       TEMP(I,J)=TEMP(I,J)+0.5*FIELD(I,J)
                     ELSE
                       TEMP(I,J)=0.5*FIELD(I,J)
                     ENDIF
                   ELSE
                     IF(TEMP(I,J).NE.ZMISS) THEN
                       TEMP(I,J)=0.5*TEMP(I,J)
                     ENDIF
                   ENDIF
                 ENDDO
               ENDDO
               K=3*KKK/2-1
               CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,
     &                       ITABLE, IPARAM, 0, K , M,
     &                       CDATE, IFORP, MARSTYPE,
     &                       LFDB, CFDBSF, KFDB, LFDBOPEN, IUOUT)
               K=3*KKK/2
               CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,
     &                       ITABLE, IPARAM, 0, K , M,
     &                       CDATE, IFORP, MARSTYPE,
     &                       LFDB, CFDBSF, KFDB, LFDBOPEN, IUOUT)
             ENDIF

           ELSEIF( KANG.EQ.2*NANG) THEN
!            HALVING THE NUMBER OF DIRECTIONS
             IF(MOD(KKK,2).EQ.1) THEN
!              SAVE DIRECTIONAL CONTRIBUTION TO BE COMBINED
!              WITH THE NEXT ONE.
               DO J=1,NGY
                JSN=NGY-J+1
                DO I=1,NLONRGG(JSN)
                  IF(FIELD(I,J).NE.ZMISS) THEN
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
                  IF(FIELD(I,J).NE.ZMISS) THEN
                    FIELD(I,J)=0.5*(FIELD(I,J)+TEMP(I,J))
                  ELSE
                    IF(TEMP(I,J).NE.0.) THEN
                      FIELD(I,J)=0.5*TEMP(I,J) 
                    ENDIF
                  ENDIF
                ENDDO
               ENDDO

               K=KKK/2
               CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,
     &                       ITABLE, IPARAM, 0, K , M,
     &                       CDATE, IFORP, MARSTYPE,
     &                       LFDB, CFDBSF, KFDB, LFDBOPEN, IUOUT)
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
                  IF(FIELD(I,J).NE.ZMISS) THEN
                    TEMP(I,J)=TEMP(I,J)+ONETHIRD*FIELD(I,J)
                  ENDIF
                ENDDO
               ENDDO

             ELSE
               DO J=1,NGY
                JSN=NGY-J+1
                DO I=1,NLONRGG(JSN)
                  IF(FIELD(I,J).NE.ZMISS) THEN
                    FIELD(I,J)=ONETHIRD*FIELD(I,J)+TEMP(I,J)
                  ELSE
                    IF(TEMP(I,J).NE.0.) THEN
                      FIELD(I,J)=TEMP(I,J)
                    ENDIF
                  ENDIF
                ENDDO
               ENDDO

               K=KKK/3
               CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,
     &                       ITABLE, IPARAM, 0, K , M,
     &                       CDATE, IFORP, MARSTYPE,
     &                       LFDB, CFDBSF, KFDB, LFDBOPEN, IUOUT)
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
                   IF(FIELD(I,J).NE.ZMISS) THEN
                     IF(TEMP(I,J).NE.ZMISS) THEN
                       TEMP(I,J)=TWOTHIRD*TEMP(I,J)+ONETHIRD*FIELD(I,J)
                     ELSE
                       TEMP(I,J)=ONETHIRD*FIELD(I,J)
                     ENDIF
                   ELSE
                     IF(TEMP(I,J).NE.ZMISS) THEN
                       TEMP(I,J)=TWOTHIRD*TEMP(I,J)
                     ENDIF
                   ENDIF
                 ENDDO
               ENDDO
               K=2*(KKK/3)+1
               CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, TEMP,
     &                       ITABLE, IPARAM, 0, K , M,
     &                       CDATE, IFORP, MARSTYPE,
     &                       LFDB, CFDBSF, KFDB, LFDBOPEN, IUOUT)

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
                   IF(FIELD(I,J).NE.ZMISS) THEN
                     IF(TEMP(I,J).NE.ZMISS) THEN
                       TEMP(I,J)=ONETHIRD*TEMP(I,J)+TWOTHIRD*FIELD(I,J)
                     ELSE
                       TEMP(I,J)=TWOTHIRD*FIELD(I,J)
                     ENDIF
                   ELSE
                     IF(TEMP(I,J).NE.ZMISS) THEN
                       TEMP(I,J)=ONETHIRD*TEMP(I,J)
                     ENDIF
                   ENDIF
                 ENDDO
               ENDDO
               K=2*(KKK/3)
               CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, TEMP,
     &                       ITABLE, IPARAM, 0, K , M,
     &                       CDATE, IFORP, MARSTYPE,
     &                       LFDB, CFDBSF, KFDB, LFDBOPEN, IUOUT)

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

           CALL IGRIB_RELEASE(KGRIB_HANDLE)
           KGRIB_HANDLE=0

2021       ISIZE=NBIT
           KBYTES=ISIZE*NPRECI
           ALLOCATE(KGRIB_BUFR(ISIZE))
           CALL IGRIB_READ_FROM_FILE(KFILE_HANDLE,KGRIB_BUFR,
     &                               KBYTES,IRET)
           IF(IRET.EQ.JPGRIB_BUFFER_TOO_SMALL) THEN
!!!          *IGRIB_READ_FROM_FILE* does not read through the file if
!!!          the size is too small, so figure out the size and read again.
             CALL KGRIBSIZE(IU06, KBYTES, NBIT, 'INTWAMINPUT 2')
             DEALLOCATE(KGRIB_BUFR)
             GOTO 2021
           ELSEIF(IRET.EQ.JPGRIB_END_OF_FILE) THEN
             IF(IALLFLD.EQ.1) THEN
               WRITE(IU06,*) ' ONLY ONE SPECTRAL FIELD HAS BEEN READ!'
               EXIT
             ELSE
               WRITE(IU06,*) '**********************************'
               WRITE(IU06,*) '*   END OF FILE ENCOUNTED'
               WRITE(IU06,*) '**********************************'
               CALL ABORT1
             ENDIF
           ELSEIF(IRET.NE.JPGRIB_SUCCESS) THEN
             WRITE(IU06,*) '*************************************'
             WRITE(IU06,*) '* INTWAMINPUT: FILE HANDLING ERROR'
             WRITE(IU06,*) '*************************************'
             CALL ABORT1
           ENDIF

           KGRIB_HANDLE=-99
           CALL IGRIB_NEW_FROM_MESSAGE(KGRIB_HANDLE,KGRIB_BUFR)

         ENDDO


!        SET MODEL VALUES TO MISSING IF NO VALUES WERE PROVIDED 
!        BY INPUT DATA 
         IF(IALLFLD.GT.1) THEN
           FIELD=ZMISS
           DO M=MLAST+1,NFRE
             DO K=1,NANG
               CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,
     &                       ITABLE, IPARAM, 0, K , M,
     &                       CDATE, IFORP, MARSTYPE,
     &                       LFDB, CFDBSF, KFDB, LFDBOPEN, IUOUT)
             ENDDO
           ENDDO
         ENDIF

       ELSE
!      ALL THE OTHER WAVE PARAMETERS

!        GET THE DATA AND PERFORM INTERPOLATION TO WAVE MODEL GRID.

         IF(.NOT.ALLOCATED(FIELD)) ALLOCATE(FIELD(NGX,NGY))
         CALL GRIB2WGRID (IU06, ITEST, NPROMA_WAM, 
     &                    KGRIB_HANDLE, KGRIB_BUFR, ISIZE,
     &                    LLUNSTR,
     &                    NXFF, NYFF, NLONRGG_LOC,
     &                    IRGG, XDELLA, ZDELLO,
     &                    XLON, YLAT,
     &                    ZMISS, PPREC, PPEPS,
     &                    CDATE, IFORP, JPARAM, KZLEV,KKK,MMM,FIELD)

         DEALLOCATE(KGRIB_BUFR)


!        ENCODE AND OUTPUT INTEGRATED PARAMETER
!        --------------------------------------

         KKK=0
         MMM=0

         CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,
     &                 ITABLE, IPARAM, KZLEV, KKK , MMM,
     &                 CDATE, IFORP, MARSTYPE,
     &                 LFDB, CFDBSF, KFDB, LFDBOPEN, IUOUT)

       ENDIF

8888  CONTINUE

      CALL IGRIB_RELEASE(KGRIB_HANDLE)
      CALL IGRIB_CLOSE_FILE(KFILE_HANDLE)
      IF(LLINTERPOL) CALL IGRIB_CLOSE_FILE(IUOUT)
 
      CALL MPL_END()
      WRITE (IU06,*) ' PROGRAM INTWAMINPUT: ALL DONE'

      END
