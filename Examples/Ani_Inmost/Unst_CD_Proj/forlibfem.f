C ======================================================================
      Subroutine FEM3Dext(XY1, XY2, XY3, XY4,
     &                        lbE, lbF, lbR, lbP, DDATA, IDATA, iSYS,
     &                        LDA, A, F, nRow, nCol,
     &             templateR, templateC)
C ======================================================================
      implicit none
      Include 'fem3Dtet.fd'
      Include 'assemble.fd'
C ======================================================================
      Real*8   XY1(*), XY2(*), XY3(*), XY4(*)

      Integer  lbE, lbF(4), lbR(6), lbP(4), LDA
      Real*8   DDATA(*)
      Integer  IDATA(*), iSYS(*)

      Real*8   A(LDA, *), F(*)
      Integer  nRow, nCol, templateR(*), templateC(*)


C Local variables
      Integer  Ddiff, Dconv, Drhs, Dbc
      External Ddiff, Dconv, Drhs, Dbc

      Integer  ANI_Dnull
      Real*8   calVol, calEdge, tet_diameter
      External calVol, calEdge, tet_diameter, ANI_Dnull

      Integer  i,j,k,n, idum, label, ibc, ir, ic
      Real*8   B(4, 4), C(4, 4), FF, vol, DiffCoef(9,9), ConvCoef(9,9)
      Real*8   x, y, z, eBC(1), Dbb, bb, bnrm, h
      Real*8   Pe, delta, dT, ConcAvr

!====================================================
       Real*8   D(4, 4),u_loc(4), DDATASIMP(4), res(4)
      Integer  IDATASIMP(4)
      Integer  ANI_D1x1_scalar
      External ANI_D1x1_scalar

C ======================================================================
      nRow = 4
      nCol = 4

      n       = iSYS(3)   ! index of triangle
      dT      = DDATA(5)   ! time step
      ConcAvr = DDATA(5+n) ! mean value of 2u^{n-1}-u^{n}/2

c      write(*,*) DDATA(1),DDATA(2),DDATA(3),DDATA(4)
c ... set up templates
      Do i = 1, nRow
         templateR(i) = Vdof   
         templateC(i) = Vdof   

      End do


!!!==============================================
      do i=1,4
      u_loc(i) = DDATA(5+iSYS(21)+iSYS(3+i))
      enddo
      DDATASIMP(1) = 1.5/dT
      IDATASIMP(1) = 1

c ... compute the stiffness matrix 
      label = lbE
!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              IDEN, FEM_P1, IDEN, FEM_P1,
     &          label, ANI_D1x1_scalar, DDATASIMP, IDATASIMP, iSYS, 2,
     &              4, D, ir, ic)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              GRAD, FEM_P1, GRAD, FEM_P1,
     &              label, Ddiff, DDATA, IDATA, iSYS, 1,
     &              4, A, ir, ic)

c ... compute the convection matrix 
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              GRAD, FEM_P1, IDEN, FEM_P1,
     &              label, Dconv, DDATA, IDATA, iSYS, 2,
     &              4, B, ir, ic)

      Do i = 1, 4
         Do j = 1, 4
            A(i, j) = A(i, j) + B(i, j)
         End do
      End do



c ... compute the mass matrix 
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              IDEN, FEM_P1, IDEN, FEM_P1,
     &              label, ANI_Dnull, DDATA, IDATA, iSYS, 2,
     &              4, B, ir, ic)
c ... add scaled mass matrix by 1.5/DeltaT
      Do i = 1, 4
         Do j = 1, 4
            A(i, j) = A(i, j) + B(i, j) * 1.5d0 / dT
         End do
      End do

c ... compute right hand side
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              IDEN, FEM_P0, IDEN, FEM_P1,
     &              lbE, Drhs, DDATA, IDATA, iSYS, 1,
     &              LDA, F, ir, ic)

c ... SUPG contribution
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              IDEN, FEM_P0, IDEN, FEM_P0,
     &              lbE, Drhs, DDATA, IDATA, iSYS, 1,
     &              1, FF, ir, ic)

      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              GRAD, FEM_P1, IDEN, FEM_P0,
     &              label, Dconv, DDATA, IDATA, iSYS, 1,
     &              4, C, ir, ic)

      vol = dabs(calVol(XY1, XY2, XY3, XY4))
      Call tet_center(XY1, XY2, XY3, XY4, x, y, z)

      idum = Ddiff(x, y, z, label, DDATA, IDATA, iSYS, DiffCoef)
      idum = Dconv(x, y, z, label, DDATA, IDATA, iSYS, ConvCoef)
      
      Dbb = 0D0
      Do i = 1, 3
         Do j = 1, 3
            Dbb = Dbb + DiffCoef(i,j) * ConvCoef(1,i) * ConvCoef(1,j)
         End do
      End do

      bb  = ConvCoef(1,1)**2 + ConvCoef(1,2)**2 + ConvCoef(1,3)**2
      bnrm = dsqrt(bb)

      h = tet_diameter(XY1, XY2, XY3, XY4)

      If(bb.NE.0D0) Then
          Dbb = Dbb / bb
          Pe = h * bnrm / Dbb
      Else
          Pe = 0D0
      End if
     
      delta = 0D0
      If(Pe.GE.1) Then
c        delta = (h / (2*bnrm)) * (1D0 - 1D0/Pe)
         delta = 1D-2
      End if 

      Do i = 1, 4
         Do j = 1, 4
            A(i, j) = A(i, j) + delta * C(1,i) * 
     &              ( C(1, j) / vol + 0.5d0 / dT )
         End do
      End do
      Do i = 1, 4
         F(i) = F(i) + delta * C(1, i) * ( FF / vol + ConcAvr / dT )
      End do

      Do i = 1, 4
         res(i) = 0.0
 	do j=1,4
		res(i) = res(i) + D(i,j)*u_loc(j)
	enddo
!         if(iSYS(3+i).eq.1) then
!           write(*,*) F(i),res(i)
! 	endif
!         F(i) = F(i) + res(i)
      End do
	


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

            ibc = Dbc(x, y, z, lbP(k), DDATA, IDATA, iSYS, eBC)

            If(ibc.EQ.BC_DIRICHLET) Then
               Call applyDIR(LDA, nRow, A, F, k, eBC(1))
            End if
         End if
      End do

      Do i = 1, 4
         res(i) = 0.0
 	do j=1,4
		res(i) = res(i) + D(i,j)*u_loc(j)
	enddo
         F(i) = F(i) + res(i)
      End do
      Return
      End







C ======================================================================
C  Diffusion tensor             
C ======================================================================
      Integer Function Ddiff(x, y, z, label, DDATA, IDATA, iSYS, Coef)
      include 'fem3Dtet.fd'

      Real*8  x, y, z, DDATA(*), Coef(9, 9)
      Integer label, IDATA(*), iSYS(*)
C ======================================================================
      iSYS(1) = 3
      iSYS(2) = 3

      Do i = 1, 3
         Do j = 1, 3
            Coef(i, j) = 0D0
         End do
         Coef(i, i) = DDATA(1)
      End do

      Ddiff = TENSOR_SYMMETRIC

      Return
      End



C ======================================================================
C  Convection tensor             
C ======================================================================
      Integer Function Dconv(x, y, z, label, DDATA, IDATA, iSYS, Conv)
      include 'fem3Dtet.fd'
C ======================================================================
      Real*8  x, y, z, DDATA(*), Conv(9, 9)
      Integer label, IDATA(*), iSYS(*)
C ======================================================================
      iSYS(1) = 1
      iSYS(2) = 3

      Conv(1,1) = DDATA(2)
      Conv(1,2) = DDATA(3)
      Conv(1,3) = DDATA(4)

      Dconv = TENSOR_GENERAL  

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
      iSYS(1) = 1
      iSYS(2) = 1

      If(label.EQ.1) Then  
         Dbc = BC_DIRICHLET
         eBC(1) = 0D0
      Else If(label.EQ.2) Then 
         Dbc = BC_DIRICHLET
         eBC(1) = 1D0
      Else
         Write(*,*) 'Dbc: wrong label=', label
         Stop
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
      iSYS(1) = 1
      iSYS(2) = 1

      F(1) = 0D0
      Drhs = TENSOR_SCALAR

      Return
      End


