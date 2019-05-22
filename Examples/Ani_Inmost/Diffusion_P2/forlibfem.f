C ================================================================
      Subroutine FEM3Dext(XY1, XY2, XY3, XY4,
     &                    lbE, lbF, lbR, lbP, DDATA, IDATA, iSYS,
     &                    LDA, A, F, nRow, nCol,
     &                    templateR, templateC)
C ================================================================
      Include 'fem3Dtet.fd'
      Include 'assemble.fd'

C ================================================================
      Real*8  XY1(*), XY2(*), XY3(*), XY4(*)

      Integer lbE, lbF(4), lbR(6), lbP(4)
      Real*8  DDATA(*)
      Integer IDATA(*), iSYS(*)

      Real*8  A(LDA, *), F(*)
      Integer templateR(*), templateC(*)

C Local variables
      Integer  Ddiff, Dconv, Drhs, Dbc
      External Ddiff, Dconv, Drhs, Dbc

      Real*8   XYP(3, 4)
      Real*8   x, y, z, eBC(1)

C ================================================================
      nRow = 10
      nCol = 10

c ... set up templates
      Do i = 1, 4
         templateR(i) = Vdof   
      End do
      Do i = 1, 6
         templateR(i+4) = Rdof   
      End do

      Do k = 1, nCol
         templateC(k) = templateR(k)
      End do

c ... compute the stiffness matrix 
      Do i = 1, 3
         XYP(i, 1) = XY1(i)
         XYP(i, 2) = XY2(i)
         XYP(i, 3) = XY3(i)
         XYP(i, 4) = XY4(i)
      End do


      label = lbE

      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              GRAD, FEM_P2, GRAD, FEM_P2,
     &              label, Ddiff, DDATA, IDATA, iSYS, 5,
     &              LDA, A, ir, ic)






c ... compute right hand side
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              IDEN, FEM_P0, IDEN, FEM_P2,
     &              lbE, Drhs, DDATA, IDATA, iSYS, 5,
     &              LDA, F, ir, ic)

c ... impose boundary conditions (assume nRow = nCol) at nodes of triangle
      Do k = 1, 4
         If(lbP(k).ne.0) Then
            If(k.EQ.1) Then
               x = XY1(1)
               y = XY1(2)
               z = XY1(3)
            ElseIf(k.EQ.2) Then
               x = XY2(1)
               y = XY2(2)
               z = XY2(3)
            ElseIf(k.EQ.3) Then
               x = XY3(1)
               y = XY3(2)
               z = XY3(3)
            ElseIf(k.EQ.4) Then
               x = XY4(1)
               y = XY4(2)
               z = XY4(3)
            End if

c  ...  Dbc may change iSYS
            ibc = Dbc(x, y, z, lbP(k), DDATA, IDATA, iSYS, eBC)

            If(ibc.EQ.BC_DIRICHLET) Then
               Call applyDIR(LDA, nRow, A, F, k, eBC(1))
              !  write(*,*) iSYS(3),iSYS(5)
!               do i=1,10 
    !           write(*,*) (A(i,k),i=1,10), F(k)

            End if
         End if
      End do
!edge bc
      k = 0
      Do i = 1, 3
         Do j = i + 1, 4  
            k = k + 1

            If( (lbP(i).GT.0 .AND. lbP(j).GT.0) ) Then
               x = (XYP(1, i) + XYP(1, j)) / 2
               y = (XYP(2, i) + XYP(2, j)) / 2
               z = (XYP(3, i) + XYP(3, j)) / 2

c  ...  Dbc will change iSYS(1:2)
               lbc = min( lbP(i), lbP(j) )
!               If (lbc.eq.9.and.x.le.0d0) lbc = 3 ! reset label for edge on inhomogeneous Dirichlet face
               ibc = Dbc(x, y, z, lbc, DDATA, IDATA, iSYS, eBC)

               If(ibc.EQ.BC_DIRICHLET) Then
                  Call applyDIR(LDA, nRow, A, F, k+4,  eBC(1))
                End if
            End if
         End do
      End do





      Return
      End



C ======================================================================
C  Diffusion tensor             
C ======================================================================
      Integer Function Ddiff(x, y, z, label, DDATA, IDATA, iSYS, Coef)
      include 'fem3Dtet.fd'

      Real*8  x, y, z, DDATA(*), Coef(1, 1)
      Integer label, IDATA(*), iSYS(*)
C ======================================================================
      iSYS(1) = 1
      iSYS(2) = 1

       Coef(1, 1) = 1.0
      Ddiff = TENSOR_SCALAR

      Return
      End





C ======================================================================
C Boundary condition
C ======================================================================
      Integer Function Dbc(x, y, z, label, DDATA, IDATA, iSYS, eBC)
      Include 'assemble.fd'
C ======================================================================
      Real*8  x, y, z, DDATA(*), eBC(*)
      Integer label, IDATA(*), iSYS(*)

C ======================================================================
      real*8 M_PI
      M_PI = 3.14159265358979323846
      iSYS(1) = 1
      iSYS(2) = 1

      Dbc = BC_DIRICHLET
      eBC(1)  = sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z) !/(3*M_PI*M_PI)

      Return
      End



C ======================================================================
C Right hand side
C ======================================================================
      Integer Function Drhs(x, y, z, label, DDATA, IDATA, iSYS, F)
      Include 'fem3Dtet.fd'
C ======================================================================
      Real*8  x, y, z, DDATA(*), F(*)
      Integer label, IDATA(*), iSYS(*)

C ======================================================================
      Real*8 M_PI
       M_PI = 3.14159265358979323846
      iSYS(1) = 1
      iSYS(2) = 1

      F(1) = sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z)*(3*M_PI*M_PI)
      Drhs = TENSOR_SCALAR

      Return
      End


