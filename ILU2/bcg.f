c---->------------------------------------------------------------------<
c  Preconditioned bi-conjugate gradient stabilized method
c  Templates for the solution if linear Systems...
c  http://www.netlib.org
c---->------------------------------------------------------------------<
      SUBROUTINE slpbcgs(
     >  prevec, IPREVEC, iW,rW, 
     >  matvec, IMATVEC, IA,JA,A,
     >  WORK, MW, NW, 
     >  N, RHS, SOL,
     >  ITER, RESID,
     >  INFO, NUNIT )
c---->
      IMPLICIT NONE
c---->------------------------------------------------------------------<
c  Argument types:
c
      EXTERNAL  matvec, prevec
      INTEGER   IMATVEC(*), IPREVEC(*), IA(*), JA(*), iW(*)
      INTEGER   N, MW, NW, ITER, INFO, NUNIT
      REAL*8    RESID, A(*)
      REAL*8    RHS(*), SOL(*), WORK(MW,NW), rW(*)
c---->
c  Argument Descriptions:
c 
c  prevec   : extern : Precondition-vector routine
c  IPREVEC  : input  : Configuration data for 'prevec'
c  matvec   : extern : Matrix-vector multiply routine
c  IMATVEC  : input  : Configuration data for 'matvec'
c
c  WORK     : work   : Workspace (MW,NW)
c  MW       : input  : leading  dimension of workspace >= N
c  NW       : input  : trailing dimension of workspace >= 8
c
c  N        : input  : Length of vectors
c  RHS      : input  : RHS vector
c  SOL      : in/out : Initial guess / iterated solution
c  ITER     : in/out : Maximum iterations / actual iterations
c  RESID    : in/out : Convergence target / Norm of final residual
c  INFO     : output : = 0, converged
c                    : = 1, did not converge
c                    : = 2, could not continue since OMEGA = ZERO
c                    : = 3, could not continue since || S || is too small
c                    : = 4, could not continue since RHO = ZERO
c                    : < 0, error with input
c---->
c  External routine specifications:
c
c    matvec( IMATVEC, A, X, B, Y )  <=>  Y = A * Mat * X + B * Y
c    prevec( IPREVEC, i, X, Y )  <=>  Y = (MatP_{i})^{-1} * X
c      where MatP is the approximation of Mat
c---->------------------------------------------------------------------<
c  Local Parameters
c
      REAL*8    ZERO,ONE
      PARAMETER ( ZERO = 0.0 , ONE = 1.0 )
c---->------------------------------------------------------------------<
c  Local Variables:
c
      INTEGER MAXIT
      INTEGER JR, JP, JS, JR1, JP1, JS1, JV, JT, i
      REAL*8  RHO, RHOPREV, OMEGA, ALPHA, BETA, TOL, TMP(1), TMP2
     & ,minSol,maxSol
c
c---->------------------------------------------------------------------<
c  External BLAS, etc.:
c
      EXTERNAL  ddot,daxpy,dcopy,dscal
      REAL*8    ddot
      INTRINSIC sqrt, min
c---->------------------------------------------------------------------<
c
c    Test the input parameters.
c
      INFO = 0
c
      if ( N .eq. 0 ) then
         return
      else if ( N .lt. 0 ) then
         INFO = -10
      else if ( MW .lt. N ) then
         INFO = -20
      else if ( NW .lt. 8 ) then
         INFO = -30
      else if ( ITER .le. 0 ) then
         INFO = -40
      endif
c
      if ( INFO .ne. 0 ) return
c---->------------------------------------------------------------------<
c  Save input iteration limit and convergence tolerance
c
      MAXIT = ITER
      TOL   = RESID
c---->
c  Alias workspace columns.
c
      JV  = 1
      JT  = JV + 1
      JR  = JT + 1
      JP  = JR + 1
      JS  = JP + 1
      JR1 = JS + 1
      JP1 = JR1+ 1
      JS1 = JP1+ 1
c---->
c  Set initial residual 
c
      call dcopy( N, RHS, 1, WORK(1,JR), 1 )
c
      TMP2 = ddot( N, SOL, 1, SOL, 1 )
      if ( TMP2 .ne. ZERO ) then
        call matvec( IMATVEC, -ONE, SOL, ONE, WORK(1,JR), IA,JA,A )
c       call matvec( IMATVEC, -ONE, SOL, ZERO, WORK(1,JR) )
c       call daxpy( N, ONE, RHS, 1, WORK(1,JR), 1 )
      endif
c---->
      TMP2 = ddot( N, WORK(1,JR), 1, WORK(1,JR), 1 )
      RESID = sqrt( TMP2 )
c---->
      ITER = 0
      RHO = ZERO
      if ( RESID .lt. TOL ) GOTO 20

      call dcopy( N, WORK(1,JR), 1, WORK(1,JR1), 1 )
c---->------------------------------------------------------------------<
c  PBCGS  iteration point
c---->--
   10   continue
c
          ITER = ITER + 1
c---->----
          RHOPREV = RHO
          RHO = ddot( N, WORK(1,JR), 1, WORK(1,JR1), 1 )

          if ( RHO .eq. ZERO ) then
             if ( NUNIT .gt. 0 ) then
                write(NUNIT,*) "PBCGS: Bad rho_tilde: method fails" 
             end if
             INFO = 4
             goto 20
          end if

          IF (ITER.eq.1) THEN
             call dcopy( N, WORK(1,JR), 1, WORK(1,JP), 1 )
             TMP(1) = ZERO
             call dcopy( N, TMP, 0, WORK(1,JV), 1 )
             call dcopy( N, TMP, 0, WORK(1,JT), 1 )
          ELSE
             BETA = ( RHO / RHOPREV ) * ( ALPHA / OMEGA ) 
             call daxpy( N, -OMEGA, WORK(1,JV), 1, WORK(1,JP), 1 )
             call dscal( N, BETA, WORK(1,JP), 1 )
             call daxpy( N, ONE, WORK(1,JR), 1, WORK(1,JP), 1 )
          END IF

          call prevec( IPREVEC, 0, WORK(1,JP), WORK(1,JP1), iW,rW )
          call matvec( IMATVEC, ONE, WORK(1,JP1), ZERO, WORK(1,JV),
     &                 IA,JA,A )
 
          TMP2 = ddot( N, WORK(1,JV), 1, WORK(1,JR1), 1 )

          ALPHA = RHO / TMP2

          call dcopy( N, WORK(1,JR), 1, WORK(1,JS), 1 )
          call daxpy( N, -ALPHA, WORK(1,JV), 1, WORK(1,JS), 1 )

          TMP2 = ddot( N, WORK(1,JS), 1, WORK(1,JS), 1 )

          if ( sqrt( TMP2 ) .lt. 1d-16 ) then
            call daxpy( N, ALPHA, WORK(1,JP1), 1, SOL, 1 )

            call matvec( IMATVEC, -ONE, SOL, ZERO, WORK(1,JR), IA,JA,A )
            call daxpy( N, ONE, RHS, 1, WORK(1,JR), 1 )
        
            TMP2 = ddot( N, WORK(1,JR), 1, WORK(1,JR), 1 )
            RESID = sqrt( TMP2 )
            INFO = 3
            goto 20 
          end if

          call prevec( IPREVEC, 0, WORK(1,JS), WORK(1,JS1), iW,rW )
          call matvec( IMATVEC, ONE, WORK(1,JS1), ZERO, WORK(1,JT), 
     &                 IA,JA,A )

          OMEGA = ddot( N, WORK(1,JT), 1, WORK(1,JS), 1 )

          TMP2 = ddot( N, WORK(1,JT), 1, WORK(1,JT), 1 )

          OMEGA = OMEGA / TMP2

          call daxpy( N, ALPHA, WORK(1,JP1), 1, SOL, 1 )
          call daxpy( N, OMEGA, WORK(1,JS1), 1, SOL, 1 )

          call dcopy( N, WORK(1,JS), 1, WORK(1,JR), 1 )
          call daxpy( N, -OMEGA, WORK(1,JT), 1, WORK(1,JR), 1 )

c---->----
c  Check convergence
          TMP2 = ddot( N, WORK(1,JR), 1, WORK(1,JR), 1 )
          RESID = sqrt( TMP2 )
c---->------------------------------------------------------------------<
c  Continue BPCGS loop while:
c    1)  Less than maximum iteration
c    2)  Have not converged 
c         print*,'pbcgs: ',ITER, RESID
       if (NUNIT.gt.0)    write(NUNIT,*) ITER, RESID , '[A'
c  For continuation it is necessary that OMEGA .ne. ZERO
c
          if ( ITER .lt. MAXIT .and. RESID .ge. TOL 
     >         .and. OMEGA .ne. ZERO ) go to 10
c---->--
c
c  Convergence failure?
c
        if ( ITER .ge. MAXIT .and. RESID .ge. TOL ) INFO = 1
        if ( OMEGA .eq. ZERO )                      INFO = 2
c---->------------------------------------------------------------------<
c  Output
c
  20    continue
        TMP2 = ddot( N , SOL, 1, SOL, 1 )
        TMP2 = sqrt( TMP2 )
        if ( NUNIT .gt. 0 ) then
          minSol=1d30
          maxSol=-1d30
          do i=1,n
             minSol=min(SOL(i),minSol)
             maxSol=max(SOL(i),maxSol)
          end do
          WRITE(NUNIT,9000) ITER,RESID,TMP2,minSol,maxSol
 9000     FORMAT('ILU: iters:',I4,', residual=',E9.3,' (|SOL|= ',E10.4,
     &           ', ',E11.4,' < SOL < ',E11.4,')')
        end if
c---->------------------------------------------------------------------<
      return
      end
c---->------------------------------------------------------------------<
