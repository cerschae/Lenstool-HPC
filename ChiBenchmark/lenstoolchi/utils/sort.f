c=======================================================================
c Numerical recipes sorting routines.
c=======================================================================

      SUBROUTINE SORT(RA,N)
	integer n
	real RA(N)
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            goto 30
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
 30   call flip(ra,n)	
      return
	END
	
c=======================================================================

      SUBROUTINE SORT2(N,RA,RB)
      DIMENSION RA(N),RB(N),WKSP(N),IWKSP(N)
      CALL INDEXX(N,RA,IWKSP)
      DO 11 J=1,N
        WKSP(J)=RA(J)
 11   CONTINUE
      DO 12 J=1,N
        RA(J)=WKSP(IWKSP(J))
 12   CONTINUE
      DO 13 J=1,N
        WKSP(J)=RB(J)
 13   CONTINUE
      DO 14 J=1,N
        RB(J)=WKSP(IWKSP(J))
 14   CONTINUE
      call flip(ra,n)	
      call flip(rb,n)
      RETURN
      END

c=======================================================================

      SUBROUTINE SORT3(N,RA,RB,RC)
      DIMENSION RA(N),RB(N),RC(N),WKSP(N),IWKSP(N)
      CALL INDEXX(N,RA,IWKSP)
      DO 11 J=1,N
        WKSP(J)=RA(J)
11    CONTINUE
      DO 12 J=1,N
        RA(J)=WKSP(IWKSP(J))
12    CONTINUE
      DO 13 J=1,N
        WKSP(J)=RB(J)
13    CONTINUE
      DO 14 J=1,N
        RB(J)=WKSP(IWKSP(J))
14    CONTINUE
      DO 15 J=1,N
        WKSP(J)=RC(J)
15    CONTINUE
      DO 16 J=1,N
        RC(J)=WKSP(IWKSP(J))
16    CONTINUE
      call flip(ra,n)	
      call flip(rb,n)
      call flip(rc,n)
      RETURN
      END

c=======================================================================

      SUBROUTINE INDEXX(N,ARRIN,INDX)
      DIMENSION ARRIN(N),INDX(N)
      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
      END

c=======================================================================
c 
c       SUBROUTINE SORT2(N,RA,RB)
c       DIMENSION RA(N),RB(N)
c       L=N/2+1
c       IR=N
c 10    CONTINUE
c         IF(L.GT.1)THEN
c           L=L-1
c           RRA=RA(L)
c           RRB=RB(L)
c         ELSE
c           RRA=RA(IR)
c           RRB=RB(IR)
c           RA(IR)=RA(1)
c           RB(IR)=RB(1)
c           IR=IR-1
c           IF(IR.EQ.1)THEN
c             RA(1)=RRA
c             RB(1)=RRB
c             RETURN
c           ENDIF
c         ENDIF
c         I=L
c         J=L+L
c 20      IF(J.LE.IR)THEN
c           IF(J.LT.IR)THEN
c             IF(RA(J).LT.RA(J+1))J=J+1
c           ENDIF
c           IF(RRA.LT.RA(J))THEN
c             RA(I)=RA(J)
c             RB(I)=RB(J)
c             I=J
c             J=J+J
c           ELSE
c             J=IR+1
c           ENDIF
c         GO TO 20
c         ENDIF
c         RA(I)=RRA
c         RB(I)=RRB
c       GO TO 10
c  30   call flip(ra,n)	
c       call flip(rb,n)
c       return	
c       END
c 		
c=======================================================================

      SUBROUTINE iSORT(IA,N)
	integer n,IA(N)
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          IIA=IA(L)
        ELSE
          IIA=IA(IR)
          IA(IR)=IA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            IA(1)=IIA
            goto 30
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(IA(J).LT.IA(J+1))J=J+1
          ENDIF
          IF(IIA.LT.IA(J))THEN
            IA(I)=IA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        IA(I)=IIA
      GO TO 10
 30   call iflip(IA,n)	
      return
	END
	
c=======================================================================

      subroutine flip(x,n)
      integer i,j,n
      real x(n),temp
      do i=1,n/2
        j = n-(i-1)
        temp = x(i)
        x(i) = x(j)
        x(j) = temp
      enddo
      return
      end

c=======================================================================

      subroutine iflip(x,n)
      integer i,j,n,x(n),temp
      do i=1,n/2
        j = n-(i-1)
        temp = x(i)
        x(i) = x(j)
        x(j) = temp
      enddo
      return
      end

c=======================================================================

