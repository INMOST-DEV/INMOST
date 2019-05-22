C ======================================================================
C  Diffusion tensor             
C ======================================================================
      Integer Function Ddiff(x, y, z, label, DDATA, IDATA, iSYS, Coef)
      include 'fem3Dtet.fd'
C ======================================================================
C  Diffusion tensor is passed through DDATA array. DDATA is filled 
C  in main.f using subroutine Drotate at the end of this file.
C ======================================================================
      Real*8  x, y, z, DDATA(3,*), Coef(9, 9)
      Integer label, IDATA(*), iSYS(*)
C ======================================================================
      iSYS(1) = 3
      iSYS(2) = 3

c HMFEM needs inverted diffusion tensor; we use Lapack for this
      Do i = 1, 3
         Do j = i, 3
            Coef(i, j) = DDATA(i, j)
         End do
      End do

      Call invertSPDmatrix(3, Coef, 9)

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
      iSYS(1) = 1
      iSYS(2) = 1

      If (label.ge.1.and.label.le.6) then
          eBC(1) = 0D0
          Dbc = BC_DIRICHLET
      Else if (label.ge.7.and.label.le.12) then
          eBC(1) = 2D0
          Dbc = BC_DIRICHLET
      Else
          stop 'wrong labelF'
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



C ======================================================================
c Templated routine for elemental matrix
C ======================================================================
      Subroutine FEM3Dext(
     &           XY1, XY2, XY3, XY4,
     &           lbE, lbF, lbR, lbP, DDATA, IDATA, iSYS,
     &           LDA, A, F, nRow, nCol,
     &           templateR, templateC)
      Implicit none
      Include 'fem3Dtet.fd'
C ======================================================================
      Real*8  XY1(*), XY2(*), XY3(*), XY4(*)

      Integer lbE, lbF(4), lbR(6), lbP(4)
      Real*8  DDATA(*)
      Integer IDATA(*), iSYS(*)

      Integer LDA, nRow, nCol
      Real*8  A(LDA, *), F(*)
      Integer templateR(*), templateC(*)

C Local variables
      Integer  Ddiff, Drhs, Dbc
      External Ddiff, Drhs, Dbc
      Real*8   tri_area 
      External tri_area

      Real*8   Q(4, 4), B(4), L(4), G(1), det, alpha, s
      Real*8   x, y, z, eBC(1)
      Integer  i, j, k, ir,ic, INFO
C ======================================================================
      nRow = 4
      nCol = 4

c ... set up templates (unknowns are located on faces)
      Do i = 1, 4
         templateR(i) = Fdof
         templateC(i) = Fdof
      End do

c ... compute the mass matrix (Q)
c     pass a copy of iSYS which is changed inside fem3Dtet
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              IDEN, FEM_RT0, IDEN, FEM_RT0,
     &              lbE, Ddiff, DDATA, IDATA, iSYS, 2,
     &              4, Q, ir, ic)

c ... compute right hand side
c     pass a copy of iSYS which is changed inside fem3Dtet
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              IDEN, FEM_P0, IDEN, FEM_P0,
     &              lbE, Drhs, DDATA, IDATA, iSYS, 2,
     &              1, G, ir, ic)

c ... compute the constraint matrix
      B(1) = tri_area(XY1, XY2, XY3)
      B(2) = tri_area(XY2, XY3, XY4)
      B(3) = tri_area(XY3, XY4, XY1)
      B(4) = tri_area(XY4, XY1, XY2)

c ... elliminate fluxes: invert the mass matrix first
      call invertSPDmatrix(4,Q,4)

      Do i = 1, 4
         Do j = 1, 4
            Q(i, j) = Q(i, j) * B(j)
         End do
      End do

      Do i = 1, 4
         L(i) = Q(1, i) * B(1) + Q(2, i) * B(2) 
     &        + Q(3, i) * B(3) + Q(4, i) * B(4)
      End do

c ... elliminate central pressure
      alpha = L(1) + L(2) + L(3) + L(4)

      Do i = 1, 4
         Do j = i, 4
            A(i, j) = B(i) * Q(i, j) - L(i) * L(j) / alpha
            A(j, i) = A(i, j)
         End do
         F(i) = L(i) * G(1) / alpha
      End do

c ... impose boundary conditions (assume nRow = nCol)
      Do k = 1, 4
         If(lbF(k).GT.0) Then
            If(k.EQ.1) Then
               x = (XY1(1) + XY2(1) + XY3(1)) / 3
               y = (XY1(2) + XY2(2) + XY3(2)) / 3
               z = (XY1(3) + XY2(3) + XY3(3)) / 3
            ElseIf(k.EQ.2) Then
               x = (XY2(1) + XY3(1) + XY4(1)) / 3
               y = (XY2(2) + XY3(2) + XY4(2)) / 3
               z = (XY2(3) + XY3(3) + XY4(3)) / 3
            ElseIf(k.EQ.3) Then
               x = (XY3(1) + XY4(1) + XY1(1)) / 3
               y = (XY3(2) + XY4(2) + XY1(2)) / 3
               z = (XY3(3) + XY4(3) + XY1(3)) / 3
            ElseIf(k.EQ.4) Then
               x = (XY4(1) + XY1(1) + XY2(1)) / 3
               y = (XY4(2) + XY1(2) + XY2(2)) / 3
               z = (XY4(3) + XY1(3) + XY2(3)) / 3
            End if

c           pass of copy of iSYS which is changed inside Dbc
            i = Dbc(x, y, z, lbF(k), DDATA, IDATA, iSYS, eBC)

            Do i = 1, nRow
               F(i) = F(i) - A(i, k) * eBC(1)
               A(k, i) = 0
               A(i, k) = 0
            End do

            A(k, k) = 1
            F(k) = eBC(1)
         End if
      End do

      Return
      End



C ======================================================================
C  The routine takes an input matrix D and rotates it about i-th
C  coordinate axis by angle theta.
C ======================================================================
      Subroutine Drotate(D, i, theta)
C ======================================================================
      Real*8  D(3, 3), theta, c, s, a, b
C ======================================================================
      j = i + 1
      If( j.GT.3 ) j = 1

      c = dcos(theta)
      s = dsin(theta)

      Do k = 1, 3
         a = D(k, i)
         b = D(k, j)

         D(k, i) = a * c - b * s
         D(k, j) = a * s + b * c
      End do

      Do k = 1, 3
         a = D(i, k)
         b = D(j, k)

         D(i, k) = a * c - b * s
         D(j, k) = a * s + b * c
      End do

      Return
      End

c ==========================================================
      Subroutine HMFEMrecoverUP(nv, vrt, nt, tet, material, DDATA,IDATA,
     &                          SOL, SOLatTets, VECatTets, iW, MaxWi)
c ==========================================================
c  Recovers pressures SOLatTets at cell centers and velocities 
c  VECatTets at cells centers from Lagrange multipliers SOL 
c  defined on mesh faces
c ==========================================================
      Implicit none
      include 'fem3Dtet.fd'

      Integer  nv, nt, MaxWi
      Real*8   vrt(3,*), DDATA(*)
      Integer  tet(4,*), material(*), IDATA(*),  iW(*)
      Real*8   SOL(*), SOLatTets(*), VECatTets(3,*)

      EXTERNAL FEM3Dext, Ddiff, Drhs
      Integer  Ddiff, Drhs

      Real*8    tri_area, calVol
      External  tri_area, calVol

C=====LOCAL VARIABLES======
      Integer   n,i,j,k, i1,i2,i3,i4, ir,ic, ifc, info
      Integer   nf, iIFE, inEP, iIEP, iiW, iSYS(5)

      Real*8    localA(5, 5), localF(5), xyzc(3), rW(100), s
      Real*8    face_area(4), vol

      Integer   Face2OpPnt(4)
      Data      Face2OpPnt / 4, 1, 2, 3 /

c ==========================================================
c     make a list E -> F
      iIFE = 1
      inEP = iIFE + 4 * nt
      iIEP = inEP + nv
      iiW  = iIEP + 4 * nt
      If(iiW .gt.MaxWi) Stop 'Please increase MaxWi'

      Call listE2F(nv, nf, nt, tet, iW(iIFE), iW(inEP), iW(iIEP))

c     loop over mesh cells
      Do n = 1, nt
         i1 = tet(1, n)
         i2 = tet(2, n)
         i3 = tet(3, n)
         i4 = tet(4, n)

c        compute mass matrix using material properties
c        we use the user-specified routine Ddiff()
         iSYS(1) = n
         iSYS(2) = i1
         iSYS(3) = i2
         iSYS(4) = i3
         iSYS(5) = i4

         Call fem3Dtet(
     &        vrt(1,i1), vrt(1,i2), vrt(1,i3), vrt(1,i4),
     &        IDEN, FEM_RT0, IDEN, FEM_RT0,
     &        material(n), Ddiff, DDATA, IDATA, iSYS, 2,
     &        5, localA, ir, ic)

c        add mass conservation equation (constraint matrix)
c        we use the divergent theorem here as we did in forlibfem.f
         localA(1, 5) = tri_area(vrt(1,i1), vrt(1,i2), vrt(1,i3))
         localA(2, 5) = tri_area(vrt(1,i2), vrt(1,i3), vrt(1,i4))
         localA(3, 5) = tri_area(vrt(1,i3), vrt(1,i4), vrt(1,i1))
         localA(4, 5) = tri_area(vrt(1,i4), vrt(1,i1), vrt(1,i2))

         localA(5, 5) = 0D0

         Do i = 1, 4
            face_area(i) = localA(i, 5)
         End do

c        compute right hand side (Lagrange muptipliers are in SOL)
  !       if(n.eq.227) then

!         write(*,*)iW(iIFE + 4*(n-1) + 0),iW(iIFE + 4*(n-1) + 1),
!     & iW(iIFE + 4*(n-1) + 2),iW(iIFE + 4*(n-1) + 3)
 !        write(*,*) tet(1,n),tet(2,n),tet(3,n),tet(4,n)
 !       endif

         Do i = 1, 4
            ifc = iW(iIFE + 4*(n-1) + i-1)
            

            localF(i) = localA(i, 5) * SOL(ifc)
         End do

c        we use the user-specified routine Drhs()
         Call fem3Dtet(
     &        vrt(1,i1), vrt(1,i2), vrt(1,i3), vrt(1,i4),
     &        IDEN, FEM_P0, IDEN, FEM_P0,
     &        material(n), Drhs, DDATA, IDATA, iSYS, 2,
     &        1, localF(5), ir, ic)

c        solve the 5x5 system localA x = localF. The first four 
c        entries in the solution x are velocities on tetrahedron faces.
c        The 5-th entry is the solution p at the mass center of the
c        tetrahedron. The user may use the LAPACK routine dsysv.

         Call dsysv('U', 5, 1, localA, 5, iW, localF, 5, rW, 100, info)
         If (info.ne.0) stop 'error in dsysv'

c        recover pressure
         SOLatTets(n) = localF(5)

c        recover velocity vector 
         vol = calVol(vrt(1,i1), vrt(1,i2), vrt(1,i3), vrt(1,i4))
         vol = dabs(vol)

         Do i = 1, 3
            xyzc(i) = (vrt(i,i1) + vrt(i,i2) + vrt(i,i3) + vrt(i,i4))/4
            VECatTets(i,n) = 0d0
         End do

         Do j = 1, 4 ! loop over faces
            k = Face2OpPnt(j)

            s = face_area(k) / (3*vol) * localF(j)

            i1 = tet(k, n)
            Do i = 1, 3
               VECatTets(i,n) = VECatTets(i,n) + s * (xyzc(i)-vrt(i,i1))
            End do
         End do
      End do

      Return
      End


