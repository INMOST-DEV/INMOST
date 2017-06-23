C ================================================================
      Subroutine FEM3Dext(
C ================================================================
     &           XY1, XY2, XY3, XY4,
     &           lbE, lbF, lbR, lbP, DDATA, IDATA, iSYS,
     &           LDA, A, F, nRow, nCol,
     &           templateR, templateC)
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
      Integer  Ddiff, Drhs, Dbc
      External Ddiff, Drhs, Dbc

      Real*8   x, y, z, eBC(3)

C ================================================================
      nRow = 12
      nCol = 12

c ... set up templates
      Do i = 1, nRow
         templateR(i) = Vdof   
      End do

      Do k = 1, nCol
         templateC(k) = templateR(k)
      End do

c ... compute the stiffness matrix 
      label = lbE

      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              GRAD, FEM_P1vector, GRAD, FEM_P1vector,
     &              label, Ddiff, DDATA, IDATA, iSYS, 1,
     &              LDA, A, ir, ic)

c ... compute right hand side
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              IDEN, FEM_P0, IDEN, FEM_P1vector,
     &              lbE, Drhs, DDATA, IDATA, iSYS, 1,
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
               Call applyDIR(LDA, nRow, A, F, k,   eBC(1))
               Call applyDIR(LDA, nRow, A, F, k+4, eBC(2))
               Call applyDIR(LDA, nRow, A, F, k+8, eBC(3))
            End if
         End if
      End do

      Return
      End



C ======================================================================
C  Diffusion tensor             
C ======================================================================
      Integer Function Ddiff(x, y, z, label, DDATA, IDATA, iSYS, Coef)
      include 'fem3Dtet.fd'

      Real*8  x, y, z, DDATA(*), Coef(9, 9), M, L, LM
      Integer label, IDATA(*), iSYS(*)
C ======================================================================
      iSYS(1) = 9
      iSYS(2) = 9

      M = DDATA(1)
      L = DDATA(2)

      Do i = 1, 9
         Do j = 1, 9
            Coef(i, j) = 0D0
         End do
         Coef(i, i) = M
      End do

      LM = M + L
      Coef(1, 1) = Coef(1, 1) + LM
      Coef(5, 5) = Coef(5, 5) + LM
      Coef(9, 9) = Coef(9, 9) + LM

      Coef(1, 5) = LM
      Coef(5, 1) = LM

      Coef(1, 9) = LM
      Coef(9, 1) = LM

      Coef(5, 9) = LM
      Coef(9, 5) = LM

      Ddiff = TENSOR_SYMMETRIC

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
      iSYS(1) = 3
      iSYS(2) = 1

      Dbc = BC_DIRICHLET
      If(label.EQ.1) Then
         eBC(1) = 0
         eBC(2) = 0
         eBC(3) = z
      Else
         eBC(1) = 0
         eBC(2) = 0
         eBC(3) = z
      End if

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
      iSYS(1) = 3
      iSYS(2) = 1

      F(1) = 0D0
      F(2) = 0D0
      F(3) = 0D0
      Drhs = TENSOR_GENERAL

      Return
      End


