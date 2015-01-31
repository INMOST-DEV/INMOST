C***********************************************************
C    Perform    Y:= BETA * Y + ALPHA * A * X
C    ALPHA, BETA are real numbers,
C    A is the matrix in  sparse row format:
C N          order of matrix
C A          non-zero entries of the matrix in  sparse row format
C jA         indexes of columns for entries 
C iA(1,..,n+1) pointers to jA for the first entry of each row 
C***********************************************************
      SUBROUTINE matvec( IMATVEC, ALPHA, X, BETA, Y, IA,JA,A )
      IMPLICIT NONE

      INTEGER   IMATVEC(*), N
      REAL*8    X(*), Y(*), ALPHA, BETA
      Integer   iA(*),jA(*)
      REAL*8    A(*)

      Integer  i,j

      N = IMATVEC(1)

      do i = 1, N
         Y(i) = BETA*Y(i) 
         do j = IA(i),IA(i+1)-1
            Y(i) = Y(i) + A(j)*X(JA(j))*ALPHA
         end do
      end do

      return
      end

