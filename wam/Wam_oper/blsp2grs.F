      PROGRAM BLSP2GRS
C
C ======================================================================
C
C**** blsp2grs  convert blsp file in grs file
C****           where a blsp file contains blocked WAM spectra 
C****           and a grs file contains gridded WAM spectra.
C
C        B. HANSEN   *ECMWF*   17/05/93
C
C     Interface.
C     ----------
C         INPUT.
C         ------
C           *jpuin*     standard user input to read name of blsp 
C                       and grs file.
C           *jpubls*    input unit opened for blsp file.
C           *cpugrd*    input unit for file GRID.
C           *cpupar*    input unit for WAModel parameter file.
C
C         OUTPUT.
C         -------
C           *iu_grs*    output unit opened for grs file.
C
C     Method.
C     -------
C       NONE.
C
C     Modification.
C     -------------
C       NONE.
C
C ----------------------------------------------------------------------
C
      CHARACTER* 4 cpugrd
      CHARACTER* 6 cpupar
      PARAMETER (
     . jpuin  =         5, ! input unit for user information.
     . jpuso  =         6, ! output unit for printer output.
     . jpubls =        10, ! input unit opened for blsp file.
     . cpugrd =    'GRID', ! input unit for file GRID.
     . cpupar =  'PARWAM'  ! input unit for WAModel parameter file.
     .          )
C
      INTEGER GETCLO, GETCLA, OPTVAL, system
      CHARACTER*  2 clclass
      CHARACTER*  3 cll1
      CHARACTER*  4 clexpver
      CHARACTER*  6 cldomain ; DATA cldomain /'global'/
      CHARACTER* 12 clfmt
      CHARACTER* 28 clopts
      CHARACTER*128 clarg
      CHARACTER     cloptlet
      DATA clopts/'b;g;k;l;a;c;s;f;C;V;D;h;m;F;'/
      DATA clclass,clexpver / 'RD', 'xxxx' /
C
      CHARACTER* 28 clct
      DATA clct / '__COMPILED_ON__' /
C
      CHARACTER*128 clfile(2) ! For filenames.
      CHARACTER*  5 cltask    ! Informs BLSPCON about what to do.
      CHARACTER*  6 clstat    ! For the status of a file to be opened.
C
      LOGICAL llexist, lloutspp, lllog, llfdb, llhelp
      DATA  lloutspp, lllog,  llfdb,   llhelp /
     .     .FALSE.,  .TRUE., .FALSE., .FALSE. /
C
      INTEGER isect1_3
      DATA isect1_3 / -999 /   ! Model version number
C
C     FOR TIMING:
      REAL wam_user_clock, time0, time
C
C ----------------------------------------------------------------------
C
C*       COMMON MAPSIZE DEFINES LENGTH OF PRIMARY BITMAP
C
      COMMON /MAPSIZE/ JPMAP1, JPMAP2
C
C ----------------------------------------------------------------------
C
C*       PARAMETER NAMELIST AS IN COMMON PARAM OF THE WAM-MODEL.
C        -------------------------------------------------------
C
      LOGICAL llper          !  the land sea mask is periodic (T)
      NAMELIST /PARWAM/
     .  iang, ifre, igx, igy, iblo, iiblo, iover, ioutp, ioutt,
     .  imaxc, imaxf, ibinp, iibl1, iibld, ibld, iiblc, iblc ,
     .  itaumax, iumax, iustar, ialpha, idepth, irefra, llper
C
C ----------------------------------------------------------------------
C
C* 0.    COMMON OUTPUT CONTROL.
C  ----------------------------
C
#include <outctrl.h>
      DATA csubna /"BLSP2GR"/
C
C ----------------------------------------------------------------------
C
C*    1. INITIALIZATIONS.
C        ----------------
 1000 CONTINUE
      time0=-wam_user_clock()
      WRITE(0,*)' BLSP2GRS v3.0  ', clct
      msub = 1
      mtest = 0
      msuvi = 0
      csuvi = ' '
      ctest = ' '
      llsuvi = .FALSE.
      lltest = .FALSE.
      is = 0
      it = 0
      ispecs=0
 
      clfile=' '
C
C        DEFAULT LENGTH OF PRIMARY BITMAP AND SECONDARY BITMAP.
C        ------------------------------------------------------
      jpmap1   = 1 908 000
      jpmap2   =       300
C
C ----------------------------------------------------------------------
C
C* 2.0   READ FILENMAES AND OPEN.
C        ------------------------
C
 2000 CONTINUE
      OPTVAL=GETCLO(clopts,clarg)
 
      cloptlet=CHAR(OPTVAL)
      IF ( OPTVAL .GT. 0 ) THEN
        IF (lltest) WRITE(jpuso,*) ' OPTION = ', cloptlet
 
C       GETS VARIABLE ARGUMENT FOR OPTION
        MORARG=GETCLA(clarg)
        IF (lltest) WRITE(jpuso,*) ' MORARG : ', MORARG
        IF(MORARG.NE.0) THEN
           IF (lltest)
     .     WRITE(jpuso,*)  ' ARGUMENT PRESENT =', clarg
           IF ( cloptlet .EQ. 'a' ) THEN
 2001         is=is+1
              csuvi(is)=clarg(1:7)
              MORARG=GETCLA(clarg)
              IF (lltest) WRITE(jpuso,*) ' MORARG : ', MORARG
              IF(MORARG.NE.0) THEN
                 GOTO  2001
              ELSE
                 CALL  uraopc (csubna, msub, llsuvi, lltest)
              ENDIF
           ELSEIF ( cloptlet .EQ. 'c' ) THEN
 2002         it=it+1
              ctest(it)=clarg(1:7)
              MORARG=GETCLA(clarg)
              IF(MORARG.NE.0) THEN
                 GOTO 2002
              ELSE
                 CALL  uraopc (csubna, msub, llsuvi, lltest)
              ENDIF
           ELSEIF ( cloptlet .EQ. 'k' ) THEN
              l1=len_trim(clarg)
              WRITE(cll1,'(i3)')l1
              clfmt='(i'//cll1//')'
              READ(clarg(1:l1),FMT=clfmt)msuvi
              CALL  uraopc (csubna, msub, llsuvi, lltest)
              IF (llsuvi) WRITE(jpuso,*)'      ',csubna,': '
           ELSEIF ( cloptlet .EQ. 'l' ) THEN
              l1=len_trim(clarg)
              WRITE(cll1,'(i3)')l1
              clfmt='(i'//cll1//')'
              READ(clarg(1:l1),FMT=clfmt)mtest
              CALL  uraopc (csubna, msub, llsuvi, lltest)
           ELSEIF ( cloptlet .EQ. 'b' ) THEN
              clfile(1)=clarg
           ELSEIF ( cloptlet .EQ. 'g' ) THEN
              clfile(2)=clarg
           ELSEIF ( cloptlet .EQ. 's' ) THEN
              lloutspp=.TRUE.
              l1=len_trim(clarg)
              WRITE(cll1,'(i3)')l1
              clfmt='(i'//cll1//')'
              READ(clarg(1:l1),FMT=clfmt)ispecs
           ELSEIF ( cloptlet .EQ. 'C' ) THEN
              clclass=clarg
           ELSEIF ( cloptlet .EQ. 'V' ) THEN
              largmnt=len_trim(clarg)
              clexpver(5-largmnt:4)=clarg
           ELSEIF ( cloptlet .EQ. 'm' ) THEN
              l1=len_trim(clarg)
              WRITE(cll1,'(i3)')l1
              clfmt='(i'//cll1//')'
              READ(clarg(1:l1),FMT=clfmt)isect1_3
           ELSEIF ( cloptlet .EQ. 'D' ) THEN
              IF ( clarg(1:6) .EQ. 'medite' ) THEN
                 cldomain=clarg(1:6)
C                isect1_3=204
              ELSEIF ( clarg(1:6) .EQ. 'global' ) THEN
                 cldomain=clarg(1:6)
C                isect1_3=1
              ELSEIF ( clarg(1:5) .EQ. 'north' ) THEN
                 cldomain=clarg(1:5)
C                isect1_3=304
              ELSEIF ( clarg(1:5) .EQ. 'south' ) THEN
                 cldomain=clarg(1:5)
C                isect1_3=305
              ELSEIF ( clarg(1:4) .EQ. 'both' ) THEN
                 cldomain=clarg(1:4)
C                isect1_3=306
              ELSE
                 CALL BLSP_USAGE ( llhelp )
              ENDIF
           ELSE
              CALL BLSP_USAGE ( llhelp )
           ENDIF
        ELSEIF ( cloptlet .EQ. 'F' ) THEN
          lflow=.TRUE.
        ELSEIF ( cloptlet .EQ. 'f' ) THEN
          llfdb=.TRUE.
          IF (llsuvi) THEN
            WRITE(jpuso,*)' llfdb changed to: ', llfdb
            WRITE(jpuso,*)' '
          ENDIF
        ELSEIF ( cloptlet .EQ. 'h' ) THEN
           llhelp = .TRUE.
           CALL BLSP_USAGE ( llhelp )
        ENDIF
 
        GOTO 2000
      ELSE
        IF (lltest) WRITE(jpuso,*) '   OPTVAL .LE. 0 '
      ENDIF
      IF ( llsuvi ) THEN
         WRITE(jpuso,'(a7)') ' CSUVI: '
         WRITE(jpuso, *    )'>',(CSUVI(jc),jc=1,is),'<'
         WRITE(jpuso,'(a7)') ' CTEST: '
         WRITE(jpuso, *    )'>',(CTEST(jc),jc=1,it),'<'
      ENDIF
 
      IF (isect1_3 .EQ. -999) THEN
        WRITE(0,'(/,"  NO MODEL VERSION NUMBER GIVEN!!",/,
     &                   "  *******************************",/)')
        CALL BLSP_USAGE ( llhelp )
      ENDIF

 2100 CONTINUE
C
      l1=0
      l2=0
      IF ( clfile(1) .NE. ' ') l1=len_trim(clfile(1))
      IF ( clfile(2) .NE. ' ') l2=len_trim(clfile(2))
      IF ( l1 .EQ. 0 .OR. l2 .EQ. 0 ) THEN
         WRITE(jpuso,*)'      FILENAMES MISSING:'
         WRITE(jpuso,*)'      clfile(1) : "',clfile(1)(1:17),'"'
         WRITE(jpuso,*)'      clfile(2) : "',clfile(2)(1:17),'"'
         WRITE(jpuso,*)'              '
         CALL BLSP_USAGE (llhelp)
      ELSEIF ( l1 .LT. 15 .OR. l2 .LT. 15 ) THEN
         WRITE(jpuso,*)'      FILENAMES PROPABLY',
     .                 ' NOT ACCORDING CONVENTION'
         WRITE(jpuso,*)'      clfile(1) : "',clfile(1)(1:17),'"'
         WRITE(jpuso,*)'      clfile(2) : "',clfile(2)(1:17),'"'
         WRITE(jpuso,*)'              '
         CALL BLSP_USAGE (llhelp)
      ENDIF
C
      clstat='NEW   '
      cltask='BLOCK'
      INQUIRE(FILE=clfile(1)(1:l1),EXIST=llexist)
      IF (llexist) THEN
        ifu1 = system ("[ -s "//clfile(1)(1:l1)//" ] " )
        IF ( ifu1 .EQ. 0 ) THEN
          clstat='OLD'
          cltask='GRID'
          READ (clfile(1)(4:17),'(i4,3i2,i4)') iyyyy, im, id, ih, ifcst
          IF (llsuvi) THEN
            WRITE(jpuso,   *)'IYYYY, IM, ID, IH, IFCST'
            WRITE(jpuso,2200) iyyyy, im, id, ih, ifcst 
          ENDIF
 2200     FORMAT(i5, 3i4, i7)
        ELSE
          clstat='UNKNOWN'
          llexist = .FALSE.
        ENDIF
      ENDIF
      IF(lltest)WRITE(jpuso,*) ' OPEN unit=',jpubls,' FILE=',
     .   clfile(1)(1:l1),' STATUS=',clstat,' FORM=UNFORMATTED'
      OPEN (jpubls,FILE=clfile(1)(1:l1), STATUS=clstat,
     .        FORM='UNFORMATTED')
      IF(lltest)WRITE(jpuso,*)' DONE '
C
      INQUIRE(FILE=clfile(2)(1:l2),EXIST=llexist)
      ifu2 = system ("[ -s "//clfile(2)(1:l2)//" ] " )
      IF (llexist .AND. clstat .EQ. 'OLD   ' .AND. ifu2 .EQ. 0) THEN
         WRITE(jpuso,*)'              '
         WRITE(jpuso,*)'      ++++++++++++++++++++++++++'
         WRITE(jpuso,*)'      + BOTH FILES ARE PRESENT +'
         WRITE(jpuso,*)'      ++++++++++++++++++++++++++'
         WRITE(jpuso,*)'              '
         WRITE(jpuso,*)'FILE 1: ',clfile(1)(1:l1)
         WRITE(jpuso,*)'FILE 2: ',clfile(2)(1:l1)
         WRITE(jpuso,*)'              '
         CALL BLSP_USAGE ( llhelp )
      ELSEIF (.NOT.llexist .AND. 
     .       (clstat .EQ. 'NEW   ' .OR. clstat .EQ. 'UNKNOWN' )) THEN
         WRITE(jpuso,*)'              '
         WRITE(jpuso,*)'      ++++++++++++++++++++++++++'
         WRITE(jpuso,*)'      + NO FILE IS PRESENT     +'
         WRITE(jpuso,*)'      ++++++++++++++++++++++++++'
         WRITE(jpuso,*)'              '
         WRITE(jpuso,*)'FILE 1: ',clfile(1)(1:l1)
         WRITE(jpuso,*)'FILE 2: ',clfile(2)(1:l1)
         WRITE(jpuso,*)'              '
         CALL BLSP_USAGE ( llhelp )
      ENDIF
      IF(llexist) THEN
        clstat='r'
      ELSE
         clstat='w'
      ENDIF
      CALL PBOPEN(iu_grs,clfile(2)(1:l1),clstat,iret)
      IF(lltest)WRITE(jpuso,*) ' PBOPEN unit=',iu_grs,' FILE=',
     .   clfile(2)(1:l1),' ACCESS=',clstat,' FORM=UNFORMATTED'
C
      IF     (iret.EQ.-1) THEN
        WRITE (jpuso,*) ' CANNOT OPEN FILE ',clfile(2)(1:l1)
      ELSEIF (iret .EQ. -1) THEN
        WRITE (jpuso,*) ' INVALID FILE NAME ',clfile(2)(1:l1)
      ELSEIF (iret .EQ. -3) THEN
        WRITE (jpuso,*) ' INVALID OPEN MODE ',clstat,' FOR FILE: ',
     .  clfile(2)(1:l1)
      ENDIF
      IF(iret .NE. 0) CALL abort
C
C ----------------------------------------------------------------------
C
C*    3.0 READ GRID INFORMATION FILE AND CONVERT BLSP FILE
C         INTO GRS FILE.
C         ------------------------------------------------
C
 3000 CONTINUE
      iu_parwam =  i_get_unit (-1, cpupar, 's', 'f', 0)
      REWIND(iu_parwam)
      READ(iu_parwam,parwam)
      IF ( lltest ) THEN
        WRITE(jpuso,*) ' READ PARWAM: '
        WRITE(jpuso, PARWAM)
      ENDIF
      CLOSE (iu_parwam)
C
      JPMAP1=iang*ifre*igx*igy
      JPMAP2=iang*ifre
C
      iu_grd =  i_get_unit (-1, cpugrd, 's', 'u', 0)
      CALL BLSPCON ( jpuso, iu_grd, jpubls, iu_grs, cltask,
     .  iang, ifre, igx, igy, iblo, iiblo, iover, ioutp, ioutt,
     .  imaxc, imaxf, ibinp, iibl1, iibld, ibld, iiblc, iblc,
     .  itaumax, iumax, iustar, ialpha, idepth, irefra,
     .  iyyyy, im, id, ih, lloutspp, ispecs, lllog, llfdb,
     .  ifcst, clclass, clexpver, cldomain, isect1_3 )
 4000 CONTINUE
C
      IF (llsuvi) THEN
        time=time0+wam_user_clock()
        time=time*1E-06
        WRITE (jpuso,*) ' ++++++++++++++++++++++++++++++'
        WRITE (jpuso,*) ' + TOTAL USER TIME IN SECONDS +'
        WRITE (jpuso,*) ' +                            +'
        WRITE (jpuso,*) ' + ', time
        WRITE (jpuso,*) ' +                            +'
        WRITE (jpuso,*) ' ++++++++++++++++++++++++++++++'
      ENDIF
C
      END
