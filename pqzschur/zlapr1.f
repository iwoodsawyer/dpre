      SUBROUTINE ZLAPR1( BASE, K, S, A, INCA, ALPHA, BETA, SCAL )
      IMPLICIT NONE
C
C     PURPOSE
C
C     Computes the general product of K complex scalars trying to avoid
C     over- and underflow.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     BASE    (input)  DOUBLE PRECISION
C             Machine base.
C
C     K       (input)  INTEGER
C             The number of scalars.  K >= 0.
C
C     S       (input)  INTEGER array, dimension (K)
C             The signature array. Each entry of S must be 1 or -1.
C
C     A       (input)  COMPLEX*16 array, dimension (K)
C             Vector of real scalars.
C
C     INCA    (input)  INTEGER
C             Increment for the array A. incA <> 0.
C
C     ALPHA   (output)  COMPLEX*16
C             ALPHA is a real scalar with 1.0 <= ABS(ALPHA) < BASE such
C             that
C
C                  ALPHA / BETA * BASE**(SCAL)
C
C             is the general product of A. 
C
C     BETA    (output)  COMPLEX*16
C             BETA is either 0.0 or 1.0.
C             See also the description of ALPHA.
C
C     SCAL    (output)  INTEGER
C             Scaling factor, see ALPHA.
C
C     CONTRIBUTOR
C
C     D. Kressner, Technical Univ. Berlin, Germany, Dec. 2002.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      COMPLEX*16        CONE, CZERO
      PARAMETER         ( CONE = ( 1.0D+0, 0.0D+0 ),
     $                  CZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Array Arguments ..
      INTEGER           S(*)
      COMPLEX*16        A(*)
C     .. Scalar Arguments ..
      INTEGER           INCA, K, SCAL
      DOUBLE PRECISION  BASE
      COMPLEX*16        ALPHA, BETA
C     .. Local Scalars ..
      INTEGER           I, INDA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DCMPLX
C
C     .. Executable Statements ..
C
      ALPHA = CONE
      BETA = CONE
      SCAL = 0
      INDA = 1
      DO 50  I = 1, K
         IF ( S(I).EQ.1 ) THEN
            ALPHA = ALPHA * A(INDA)
         ELSE
            IF ( A(INDA).EQ.CZERO ) THEN
               BETA = CZERO
            ELSE
               ALPHA = ALPHA / A(INDA)
            END IF
         END IF
         IF ( ABS( ALPHA ).EQ.ZERO ) THEN
            SCAL = 0
            ALPHA = CZERO
         ELSE
   10       IF ( ABS( ALPHA ).GE.ONE )  GO TO 20
               ALPHA = DCMPLX( BASE, ZERO )*ALPHA
               SCAL = SCAL - 1
            GO TO 10
   20       CONTINUE
   30       IF ( ABS( ALPHA ).LT.BASE ) GO TO 40
               ALPHA = ALPHA / DCMPLX( BASE, ZERO )
               SCAL = SCAL + 1
            GO TO 30
   40       CONTINUE
         END IF
         INDA = INDA + INCA
   50 CONTINUE
C
      RETURN
C
C *** Last line of ZLAPR1 ***
      END
