      subroutine  refine(
     & number_of_refines,
     & nv,nvmax,nt,ntmax,nb,nbmax,
     & vrt,tet,bnd ,material,labelF,
     & MapMtr,Ref2MapMtr,iW,MaxWi)
c ==========================================================
c This program loads a trivial mesh (data/cub6.out) consisting 
c of six tetrahedra, refines it uniformly 3 times, refines it
c locally 4 times and then  saves the final mesh in GMV format.
c ==========================================================
      implicit none
      
      integer number_of_refines
      integer nvmax,ntmax,nbmax

c ... nvmax - maximum number of mesh nodes
c ... ntmax - maximum number of mesh tetrahedra
c ... nbmax - maximum number of boundary faces
c      parameter(nvmax = 150 000, ntmax = 2*nvmax, nbmax = 100 000)


      Integer   MaxWi
c ... standard mesh arrays (see doc/aft_guide.pdf for more detail)
c ... number of points, tets, and boundary faces
      Integer  nv, nt, nb 

c ... work memory
      Integer  iW(*)
c ... array  keeps mappings of each coarse tetrahedron to make it equilateral
      Real*8 MapMtr(3,3,*)
c ... array keeps references of each fine cell to MapMtr
      Integer Ref2MapMtr(*)



c ... coordinates of mesh points 
      Real*8   vrt(3,*)

c ... connectivity table for tets and their labels
      Integer  tet(4,*),material(*)

c ... connectivity table for boundary faces, and their labels
      Integer  bnd(3,*),labelF(*)

c ... additional mesh arrays are needed to call mesh_metric() from libmba3D.a 
c     ivfix(nvfix) is array of fixed (never touched) points
c     ibfix(nbfix) is array of fixed (never touched) faces
c     itfix(ntfix) is array of fixed (never touched) elements
      Integer   nvfix, ivfix(6), nbfix, ibfix(1), ntfix, itfix(1) 


c Local variables
      Integer   i, j, k, ii
      Real*8    xc,yc,zc
c		write(*,*) "enter refine ",nv,nt,nb

c      write(*,*) "begin init refine"
c ... Initialize the refinement (filling MapMtr, Ref2MapMtr)
      Call initializeRefinement(
     &      nv, nt, vrt, tet,
     &      MapMtr, Ref2MapMtr)
c      write(*,*) "begin refine ",nv,nt
c ... Refine the mesh uniformly 3 times
      Do i = 1, number_of_refines
c         write(*,*) "refine ",i
         Call uniformRefinement(
     &      nv, nvmax, nb, nbmax, nt, ntmax,
     &      vrt, bnd, tet, labelF, material,
     &      MapMtr, Ref2MapMtr,
     &      iW, MaxWi)
c       write(*,*) "end refine ",nv,nt
      End do

      end subroutine
