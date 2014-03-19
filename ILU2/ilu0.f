c-----------------------------------------------------------------------
c This routine solves the system (LU) y = x
c-----------------------------------------------------------------------
        SUBROUTINE prevec0( iprevec, dummy, x, y, iwork, dwork)
        INTEGER   iprevec(*), dummy, iwork(*), n
        REAL*8    x(*),y(*),dwork(*)
        n = iprevec(1)
        call lusol(n, x, y, dwork(1), iwork(n+2), iwork(1))
        return
        end

	subroutine ilu0(n, a, ja, ia, alu, jlu, ju, iw, ierr)
	implicit real*8 (a-h,o-z)
	real*8 a(*), alu(*)
        integer ja(*), ia(*), ju(*), jlu(*), iw(*)
c input:
c---------
c n       = dimension of matrix
c a, ja,
c ia      = original matrix in compressed sparse row storage.
c
c output:   
c alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
c           the L and U factors together. The diagonal (stored in
c           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
c           contains the i-th row of L (excluding the diagonal entry=1)
c           followed by the i-th row of U.
c
c ju	  = pointer to the diagonal elements in alu, jlu.
c
c ierr	  = integer indicating error code on return
c	     ierr = 0 --> normal return
c	     ierr = k --> code encountered a zero pivot at step k.
c work arrays:
c-------------
c iw	    = integer work array of length n.
c-----------------------------------------------------------------------
      integer maxbuf, nbuf
      parameter (maxbuf=1000) ! max number of non-zeros in a row
      integer jbuf(maxbuf)
      real*8  abuf(maxbuf)
c-----------------------------------------------------------------------
        ju0 = n+2
        jlu(1) = ju0
c
c initialize work vector to zero's
c
	do 31 i=1, n
           iw(i) = 0
 31     continue
c
c main loop
c
	do 500 ii = 1, n
           js = ju0
c
c sort column indexes
c
           nbuf=0
           do j=ia(ii),ia(ii+1)-1
            nbuf=nbuf+1
            abuf(nbuf)=a(j)
            jbuf(nbuf)=ja(j)
           end do
           if (nbuf.gt.maxbuf) stop 'ilu0: maxbuf too small'
           call dsortilu (jbuf, abuf, nbuf, 2, ierr)

c generating row number ii of L and U.
c
           do 100 j=1,nbuf  ! j=ia(ii),ia(ii+1)-1
c
c     copy row ii of a, ja, ia into row ii of alu, jlu (L/U) matrix.
c
              jcol = jbuf(j)  ! ja(j)
              if (jcol .eq. ii) then
                 alu(ii) = abuf(j)  ! a(j)
                 iw(jcol) = ii
                 ju(ii)  = ju0
              else
                 alu(ju0) = abuf(j)  ! a(j)
                 jlu(ju0) = jbuf(j)  ! ja(j)
                 iw(jcol) = ju0
                 ju0 = ju0+1
              endif
 100       continue
           jlu(ii+1) = ju0
           jf = ju0-1
           jm = ju(ii)-1
c
c     exit if diagonal element is reached.
c
           do 150 j=js, jm
              jrow = jlu(j)
              tl = alu(j)*alu(jrow)
              alu(j) = tl
c
c     perform  linear combination
c
              do 140 jj = ju(jrow), jlu(jrow+1)-1
                 jw = iw(jlu(jj))
                 if (jw .ne. 0) alu(jw) = alu(jw) - tl*alu(jj)
 140          continue
 150       continue
c
c     invert  and store diagonal element.
c
           if (alu(ii) .eq. 0.0d0) goto 600
           alu(ii) = 1.0d0/alu(ii)
c
c     reset pointer iw to zero
c
           iw(ii) = 0
           do 201 i = js, jf
 201          iw(jlu(i)) = 0
 500       continue
           ierr = 0
           return
c
c     zero pivot :
c
 600       ierr = ii
c
           return
           end

        subroutine lusol(n, y, x, alu, jlu, ju)
        real*8 x(n), y(n), alu(*)
        integer n, jlu(*), ju(*)
c-----------------------------------------------------------------------
c
c This routine solves the system (LU) x = y,
c given an LU decomposition of a matrix stored in (alu, jlu, ju)
c modified sparse row format
c
c-----------------------------------------------------------------------
c on entry:
c n   = dimension of system
c y   = the right-hand-side vector
c alu, jlu, ju
c     = the LU matrix as provided from the ILU routines.
c
c on return
c x   = solution of LU x = y.
c-----------------------------------------------------------------------
c
c Note: routine is in place: call lusol (n, x, x, alu, jlu, ju)
c       will solve the system with rhs x and overwrite the result on x .
c
c-----------------------------------------------------------------------
c local variables
c
        integer i,k
c
c forward solve
c
        do 40 i = 1, n
           x(i) = y(i)
           do 41 k=jlu(i),ju(i)-1
              x(i) = x(i) - alu(k)* x(jlu(k))
 41        continue
 40     continue
c
c     backward solve.
c
        do 90 i = n, 1, -1
           do 91 k=ju(i),jlu(i+1)-1
              x(i) = x(i) - alu(k)*x(jlu(k))
 91        continue
           x(i) = alu(i)*x(i)
 90     continue
c
        return
        end
