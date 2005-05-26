      SUBROUTINE WSTREAM_STRG(ISTREAM,CSTREAM,NENSFNB,NTOTENS,
     &                        MARSFCTYPE) 

! ----------------------------------------------------------------------

!**** *WSTREAM_STRG* -

!     J. BIDLOT     ECMWF    AUGUST 2000.


!*    PURPOSE.
!     --------
!     TO RETURN A 4 CHARACTER STRING WHICH CORRESPONDS TO THE VALUE
!     ISTREAM WHEN DEALING WITH GRIB DATA AND THE FDB.
!     IF NO MATCH IS FOUND, CSTREAM IS SET TO '****'
!     BASED ON ISTREAM, NENSFNB (the ensemble number)
!                   and NTOTENS (total number of ensemble members)
!     FORECAST TYPE WILL ALSO BE RETURNED. 

!**   INTERFACE.
!     ----------
!     *CALL WSTREAM_STRG*(ISTREAM,CSTREAM)

! ----------------------------------------------------------------------
      INTEGER :: ISTREAM, NENSFNB, NTOTENS
      CHARACTER*2 :: MARSFCTYPE
      CHARACTER*4 CSTREAM

      IF(ISTREAM.EQ.1045) THEN
!       DETERMINISTIC WAVE FORECAST
        CSTREAM = 'wave'
        MARSFCTYPE = 'fc'
      ELSE IF(ISTREAM.EQ.1080) THEN
!       DETERMINISTIC WAVE MONTHLY MEANS AND CLIMATOLOGY
        CSTREAM = 'wamo'
        MARSFCTYPE = 'fc'
      ELSE IF(ISTREAM.EQ.1081) THEN
!       ENSEMBLE WAVE FORECAST
        CSTREAM = 'waef'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'pf'
        ENDIF
      ELSE IF(ISTREAM.EQ.1082) THEN
!       WAVE SEASONAL FORECAST
        CSTREAM = 'wasf'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
      ELSE IF(ISTREAM.EQ.1083) THEN
!       WAVE MULTI ANALYSIS
        CSTREAM = 'mawv'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
      ELSE IF(ISTREAM.EQ.1084) THEN
!       ENSEMBLE WAVE HINDCAST
        CSTREAM = 'ewhc'
        MARSFCTYPE = 'fc'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'pf'
        ENDIF
      ELSE IF(ISTREAM.EQ.1085) THEN
!       WAVE HINDCAST
        CSTREAM = 'wvhc'
        MARSFCTYPE = 'fc'
      ELSE IF(ISTREAM.EQ.1086) THEN
!       WAVE ENSEMBLE FORECAST OVERLAP 
        CSTREAM = 'weov'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'pf'
        ENDIF
      ELSE IF(ISTREAM.EQ.1027) THEN
!       WAVE SHORT CUT-OFF FORECAST
        CSTREAM = 'scwv'
        MARSFCTYPE = 'fc'
      ELSE IF(ISTREAM.EQ.1029) THEN
!       WAVE DELAYED CUT-OFF FORECAST
        CSTREAM = 'dcwv'
        MARSFCTYPE = 'fc'
      ELSE IF(ISTREAM.EQ.1092) THEN
!       WAVE SEASONAL MONTHLY MEANS
        CSTREAM = 'swmm'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
      ELSE IF(ISTREAM.EQ.1095) THEN
!       WAVE MONTHLY FORECAST (old)
        CSTREAM = 'wamf'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
      ELSE IF(ISTREAM.EQ.1096) THEN
!       WAVE MONTHLY FORECAST MEANS (old)
        CSTREAM = 'wmfm'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
      ELSE IF(ISTREAM.EQ.1203) THEN
!       MONTHLY FORECAST (new)
        CSTREAM = 'mnfw'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
      ELSE IF(ISTREAM.EQ.1204) THEN
!       MONTHLY FORECAST HINDCAST
        CSTREAM = 'mfhw'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
      ELSE IF(ISTREAM.EQ.1205) THEN
!       MONTHLY FORECAST ANOMALIES
        CSTREAM = 'mfaw'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
      ELSE IF(ISTREAM.EQ.1209) THEN
!       MONTHLY FORECAST MEANS
        CSTREAM = 'mfwm'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
      ELSE IF(ISTREAM.EQ.1210) THEN
!       MONTHLY FORECAST HINDCAST MEANS
        CSTREAM = 'mhwm'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
      ELSE IF(ISTREAM.EQ.1211) THEN
!       MONTHLY FORECAST HINDCAST ANOMALIES
        CSTREAM = 'mawm'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
      ELSE IF(ISTREAM.EQ.1222) THEN
!       MULTI MODEL WAVE SEASONAL FORECAST
        CSTREAM = 'wams'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
      ELSE IF(ISTREAM.EQ.1223) THEN
!       MULTI MODEL WAVE SEASONAL FORECAST MONTHLY MEANS
        CSTREAM = 'mswm'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
      ELSE
        CSTREAM='****'
        MARSFCTYPE = 'fc'
      ENDIF

      RETURN
      END SUBROUTINE WSTREAM_STRG
