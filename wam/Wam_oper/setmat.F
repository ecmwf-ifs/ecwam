      SUBROUTINE SETMAT (XM, NOBS, P, D, NDIM2)                         

! --------------------------------------------------------------------- 

!     DETERMINE MATRIX P+D                                              

      DIMENSION XM(NOBS,NOBS), P(NDIM2,NDIM2), D(NDIM2,NDIM2)           

      DO JOBS=1,NOBS                                               
        DO IOBS=1,NOBS                                            
          XM(IOBS,JOBS)=P(IOBS,JOBS)+D(IOBS,JOBS)                     
        ENDDO
      ENDDO

      RETURN                                                            
      END SUBROUTINE SETMAT




