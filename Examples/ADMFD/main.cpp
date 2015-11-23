#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "inmost.h"
#include "matrix.hpp"
using namespace INMOST;

#ifndef M_PI
#define M_PI 3.141592653589
#endif

#if defined(USE_MPI)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);
#else
#define BARRIER
#endif

//shortcuts
typedef Storage::real real;
typedef Storage::integer integer;
typedef Storage::enumerator enumerator;
typedef Storage::real_array real_array;
typedef Storage::var_array var_array;

const real reg_abs = 1.0e-12; //regularize abs(x) as sqrt(x*x+reg_abs)
const real reg_div = 1.0e-15; //regularize (|x|+reg_div)/(|x|+|y|+2*reg_div) to reduce to 1/2 when |x| ~= |y| ~= 0




int main(int argc,char ** argv)
{
	Solver::Initialize(&argc,&argv,""); // Initialize the solver and MPI activity
#if defined(USE_PARTITIONER)
	Partitioner::Initialize(&argc,&argv); // Initialize the partitioner activity
#endif
	if( argc > 1 )
	{
    double ttt; // Variable used to measure timing
    bool repartition = false; // Is it required to redistribute the mesh?
		Mesh * m = new Mesh(); // Create an empty mesh
    { // Load the mesh
		  ttt = Timer();
		  m->SetCommunicator(INMOST_MPI_COMM_WORLD); // Set the MPI communicator for the mesh
      if( m->GetProcessorRank() == 0 ) std::cout << "Processors: " << m->GetProcessorsNumber() << std::endl;
		  if( m->isParallelFileFormat(argv[1]) ) //The format is 
		  {
			  m->Load(argv[1]); // Load mesh from the parallel file format
			  repartition = true; // Ask to repartition the mesh
		  }
		  else if( m->GetProcessorRank() == 0 ) m->Load(argv[1]); // Load mesh from the serial file format
		  BARRIER
		  if( m->GetProcessorRank() == 0 ) std::cout << "Load the mesh: " << Timer()-ttt << std::endl;
    }
		
		
#if defined(USE_PARTITIONER)
    if (m->GetProcessorsNumber() > 1 && !repartition) // Currently only non-distributed meshes are supported by Inner_RCM partitioner
    {
      { // Compute mesh partitioning
			  ttt = Timer();
			  Partitioner p(m); //Create Partitioning object
			  p.SetMethod(Partitioner::Inner_RCM,repartition ? Partitioner::Repartition : Partitioner::Partition); // Specify the partitioner
			  p.Evaluate(); // Compute the partitioner and store new processor ID in the mesh
			  BARRIER
			  if( m->GetProcessorRank() == 0 ) std::cout << "Evaluate: " << Timer()-ttt << std::endl;
      }

      { //Distribute the mesh
			  ttt = Timer();
			  m->Redistribute(); // Redistribute the mesh data
			  m->ReorderEmpty(CELL|FACE|EDGE|NODE); // Clean the data after reordring
			  BARRIER
			  if( m->GetProcessorRank() == 0 ) std::cout << "Redistribute: " << Timer()-ttt << std::endl;
      }
		}
#endif

    { // prepare geometrical data on the mesh
      ttt = Timer();
		  Mesh::GeomParam table;
		  table[CENTROID]    = CELL | FACE; //Compute averaged center of mass
		  table[NORMAL]      = FACE;        //Compute normals
		  table[ORIENTATION] = FACE;        //Check and fix normal orientation
		  table[MEASURE]     = CELL | FACE; //Compute volumes and areas
		  //table[BARYCENTER]  = CELL | FACE; //Compute volumetric center of mass
		  m->PrepareGeometricData(table); //Ask to precompute the data
		  BARRIER
		  if( m->GetProcessorRank() == 0 ) std::cout << "Prepare geometric data: " << Timer()-ttt << std::endl;
    }
		
    // data tags for
		Tag tag_P;  // Pressure
		Tag tag_K;  // Diffusion tensor
    Tag tag_F;  // Forcing term
    Tag tag_BC; // Boundary conditions
    Tag tag_M;  // Stiffness matrix
    Tag tag_I;  // Temporary local indexation
    Tag tag_H;  // Harmonic points
    Tag tag_W;  // Gradient matrix acting on harmonic points on faces and returning gradient on faces
    Tag tag_D;  // Entries for scaling matrix D

    if( m->GetProcessorsNumber() > 1 ) //skip for one processor job
    { // Exchange ghost cells
		  ttt = Timer();
		  m->ExchangeGhost(1,FACE); // Produce layer of ghost cells
      m->ExchangeData(tag_H,FACE,0); //Synchronize harmonic points
		  BARRIER
		  if( m->GetProcessorRank() == 0 ) std::cout << "Exchange ghost: " << Timer()-ttt << std::endl;
    }
    
    { //initialize data
		  if( m->HaveTag("PERM") ) // is diffusion tensor already defined on the mesh? (PERM from permeability)
        tag_K = m->GetTag("PERM"); // get the diffusion tensor
      
      if( !tag_K.isValid() || !tag_K.isDefined(CELL) ) // diffusion tensor was not initialized or was not defined on cells.
      {
        tag_K = m->CreateTag("PERM",DATA_REAL,CELL,NONE,6); // create a new tag for symmetric diffusion tensor K
		    for( Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell ) // loop over mesh cells
        {
          real_array K = cell->RealArray(tag_K);
          // assign a symmetric positive definite tensor K
				  K[0] = 1.0; //XX
          K[1] = 0.0; //XY
          K[2] = 0.0; //XZ
          K[3] = 1.0; //YY
          K[4] = 0.0; //YZ
          K[5] = 1.0; //ZZ
        }

        m->ExchangeData(tag_K,CELL,0); //Exchange diffusion tensor
      }

      if( m->HaveTag("PRESSURE") ) //Is there a pressure on the mesh?
        tag_P = m->GetTag("PRESSURE"); //Get the pressure
      
      if( !tag_P.isValid() || !tag_P.isDefined(CELL) ) // Pressure was not initialized or was not defined on nodes
      {
        srand(1); // Randomization
        tag_P = m->CreateTag("PRESSURE",DATA_REAL,CELL|FACE,FACE,1); // Create a new tag for the pressure
        for(Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell) //Loop over mesh cells
          cell->Real(tag_P) = (rand()*1.0)/(RAND_MAX*1.0); // Prescribe random value in [0,1]
      }

      
      if( m->HaveTag("BOUNDARY_CONDITION") ) //Is there boundary condition on the mesh?
      {
        tag_BC = m->GetTag("BOUNDARY_CONDITION");
        
        //initialize unknowns at boundary
      }

      m->ExchangeData(tag_P,CELL|FACE,0); //Synchronize initial solution with boundary unknowns
      

      //run in a loop to identify boundary pressures so that they enter as primary variables

      //3 entries are coordinates and last entry is the coefficient for FrontCell
      tag_H = m->CreateTag("HARMONIC_POINT",DATA_REAL,FACE,NONE,4);
      // Assemble coefficients for harmonic averaging matrix H
      for(Mesh::iteratorFace fKL = m->BeginFace(); fKL != m->EndFace(); ++fKL) //go over faces
      {
        real_array H = fKL->RealArrayDF(tag_H); //access data structure
        Cell cK = fKL->BackCell(); //get cell K from the back of the normal
        Cell cL = fKL->FrontCell(); //get cell L to the front of the normal
        if( !cL.isValid() ) //this is boundary face
        {
          fKL->Centroid(H.data()); //set center of face as harmonic point
          H[3] = 0.0; //set coefficient 
          continue;
        }

        real dK, dL, lK, lL, D;
        real coefK, coefL, coefQ, coefDiv;
        
        rMatrix y(3,1); //harmonic point
        rMatrix yK(3,1); //projection of the center of cell K onto surface
        rMatrix yL(3,1); //projection of the center of cell L onto surface
        rMatrix xK(3,1); //center of cell K
        rMatrix xL(3,1); //center of cell L
        rMatrix xKL(3,1); //point on face
        rMatrix nKL(3,1); //normal to face
        rMatrix lKs(3,1); //reminder of the co-normal vector Kn in cell K when projection to normal is substracted, i.e. (Kn - n n.Kn)
        rMatrix lLs(3,1); //reminder of the co-normal vector Kn in cell L when projection to normal is substracted, i.e. (Kn - n n.Kn)
        rMatrix KK = FromTensor(cK->RealArray(tag_K)).Transpose(); //diffusion tensor in cell K
        rMatrix KL = FromTensor(cL->RealArray(tag_K)).Transpose(); //diffusion tensor in cell L

        cK->Centroid(xK.data()); //retrive center of cell K
        cL->Centroid(xL.data()); //retrive center of cell L
        fKL->Centroid(xKL.data()); //retrive center of face
        fKL->OrientedUnitNormal(cK,nKL.data()); //retrive unit normal to face

        lKs = KK*nKL; //get conormal in cell K
        lLs = KL*nKL; //get conormal in cell L
        lK = nKL.DotProduct(lKs); //find projection of conormal onto normal in cell K
        lL = nKL.DotProduct(lLs); //find projection of conormal onto normal in cell L
        lKs -= nKL*lK; //obtain reminder in cell K
        lLs -= nKL*lL; //obtain reminder in cell L

        D = -nKL.DotProduct(xKL); // Compute constant in equation for plane (nx,ny,nz,D)
        dK = nKL.DotProduct(xK)+D; //compute distance to the center of cell K
        dL = nKL.DotProduct(xL)+D; //compute distance to the center of cell L

        if( dK*dL > 0 ) //check consistency of geometry
        {
          std::cout << "Cell centers are on the same side from the face" << std::endl;
        }

        dK = fabs(dK); //signs should be encorporated automatically
        dL = fabs(dL);

        if( dK < 1.0e-9 || dL < 1.0e-9 ) std::cout << "Cell center is located on the face" << std::endl;

        yK = xK + nKL*dK; //compute projection of the center of cell K onto face
        yL = xL - nKL*dL; //compute projection of the center of cell L onto face

        //compute coefficients for harmonic point
        coefDiv = (lL*dK+lK*dL);
        coefL = lL*dK/coefDiv;
        coefK = lK*dL/coefDiv;
        coefQ = dK*dL/coefDiv;

        y = yK*coefK + yL*coefL + (lKs - lLs)*coefQ; //calculate position of harmonic point

        std::copy(y.data(),y.data()+3,H.data()); //store position of harmonic point
        H[3] = coefL; //store multiplier for pressure in FrontCell, the other is coefK = 1 - coefL
      } //end of loop over faces

      tag_M = m->CreateTag("INNER_PRODUCT",DATA_REAL,CELL,NONE);
      //assemble inner product matrix M acting on faces for each cell
      for(Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell)
      {
        real xP[3]; //center of the cell
        real nF[3]; //normal to the face
        real aF, vP = cell->Volume(); //area of the face
        cell->Centroid(xP); //obtain cell center
        ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell
        enumerator NF = faces.size(), k; //number of faces;

        rMatrix IP(NF,NF), N(NF,3), R(NF,3); //matrices for inner product

        //assemble matrices R and N
        k = 0;
        for(ElementArray<Face>::iterator face = faces.begin(); face != faces.end(); ++face)
        {
          aF = face->Area();
          real_array yF = face->RealArrayDF(tag_H); //point on face is harmonic point
          face->OrientedUnitNormal(cell->self(),nF);
          // assemble R matrix, follows chapter 3.4.3 in the book "Mimetic Finite Difference Method for Elliptic Problems" by Lourenco et al
          R(k,0) = (yF[0]-xP[0])*vP;
          R(k,1) = (yF[1]-xP[1])*vP;
          R(k,2) = (yF[2]-xP[2])*vP;
          // assemble N matrix, follow chapter 3.4. in the book and projection operator definition from
          // chapter 2.3.1 in "Mimetic scalar products of differential forms" by Brezzi et al
          N(k,0) = nF[0]*aF;
          N(k,1) = nF[1]*aF;
          N(k,2) = nF[2]*aF;
          k++;
        } //end of loop over faces

        //inner product definition from Corollary 4.2,  "Mimetic scalar products of differential forms" by Brezzi et al
        IP = R*(R.Transpose()*N).Invert()*R.Transpose(); // Consistency part
        IP += (rMatrix::Unit(NF) - N*(N.Transpose()*N).Invert()*N.Transpose())*(2.0/static_cast<real>(NF)*IP.Trace()); //Stability part
        
        //assert(IP.isSymmetric()); //test positive definitiness as well!
        /*
        if( !IP.isSymmetric() ) 
        {
          std::cout << "unsymmetric" << std::endl;
          IP.Print();
          std::cout << "check" << std::endl;
          (IP-IP.Transpose()).Print();
        }
        */
        real_array M = cell->RealArrayDV(tag_M); //access data structure for inner product matrix in mesh
        M.resize(NF*NF); //resize variable array
        std::copy(IP.data(),IP.data()+NF*NF,M.data()); //write down the inner product matrix
      } //end of loop over cells

      tag_W = m->CreateTag("GRAD",DATA_REAL,CELL,NONE);
      tag_D = m->CreateTag("DIAG",DATA_VARIABLE,CELL,NONE);
      //Assemble all-positive gradient matrix W on cells
      for(Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell)
      {
        real xP[3]; //center of the cell
        real nF[3]; //normal to the face
        cell->Centroid(xP);
        ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell
        enumerator NF = faces.size(), k; //number of faces;
        dynarray<real,64> NG(NF,0), G(NF,0); // count how many times gradient was calculated for face
        rMatrix K = FromTensor(cell->RealArrayDF(tag_K)); //get permeability for the cell
        rMatrix GRAD(NF,NF), NK(NF,3), R(NF,3); //big gradient matrix, co-normals, directions
        rMatrix GRADF(NF,NF), NKF(NF,3), RF(NF,3); //corresponding matrices with masked negative values when multiplied by individual co-normals
        rMatrix GRADFstab(NF,NF);
        GRAD.Zero();
        k = 0;
        for(ElementArray<Face>::iterator face = faces.begin(); face != faces.end(); ++face) //loop over faces
        {
          real_array yF = face->RealArrayDF(tag_H); //point on face is harmonic point
          face->OrientedUnitNormal(cell->self(),nF);
          // assemble matrix of directions
          R(k,0) = (yF[0]-xP[0]);
          R(k,1) = (yF[1]-xP[1]);
          R(k,2) = (yF[2]-xP[2]);
          // assemble matrix of co-normals 
          rMatrix nK = FromVector(nF).Transpose()*K;
          NK(k,0) = nK(0,0);
          NK(k,1) = nK(0,1);
          NK(k,2) = nK(0,2);
          k++;
        } //end of loop over faces

        GRAD = NK*(NK.Transpose()*R).Invert()*NK.Transpose(); //stability part
        //std::cout << "consistency" << std::endl;
        //GRADF.Print();

        GRAD += (rMatrix::Unit(NF) - R*(R.Transpose()*R).Invert()*R.Transpose())*(2.0/NF*GRADF.Trace());
        
        /*
        std::cout << "stability" << std::endl;
        GRADFstab.Print();
        GRADF += GRADFstab;
        std::cout << "consistency+stability" << std::endl;
        GRADF.Print();

        k = 0;
        for(ElementArray<Face>::iterator face = faces.begin(); face != faces.end(); ++face) //loop over faces
        {
          real nK[3] = {NK(k,0),NK(k,1),NK(k,2)};
          real N = 0.0;
          for(enumerator l = 0; l < NF; ++l) //run through faces
          {
            if( nK[0]*NK(l,0) + nK[1]*NK(l,1) + nK[2]*NK(l,2) > 0.0 && nK[0]*R(l,0) + nK[1]*R(l,1) + nK[2]*R(l,2) > 0.0 ) //all entries are positive
            {
              NKF(l,0) = NK(l,0);
              NKF(l,1) = NK(l,1);
              NKF(l,2) = NK(l,2);
              RF(l,0) = R(l,0);
              RF(l,1) = R(l,1);
              RF(l,2) = R(l,2);
              G[l] = 1; //mask active rows
              NG[l]++; //count number of equations
              N++; //count number of active rows
            }
            else
            {
              NKF(l,0) = NKF(l,1) = NKF(l,2) = 0.0;
              RF(l,0) = RF(l,1) = RF(l,2) = 0.0;
              G[l] = 0; //mask active rows
            }
          }
          //NKF.Print();
          //RF.Print();
          GRADF.Zero();
          GRADF = NKF*(NKF.Transpose()*RF).Invert()*NKF.Transpose(); //stability part
          std::cout << "consistency" << std::endl;
          GRADF.Print();
          GRADFstab.Zero();
          GRADFstab = (FromDiagonal(G.data(),NF) - RF*(RF.Transpose()*RF).Invert()*RF.Transpose())*(2.0/N*GRADF.Trace());
          real alpha = 1.0;
          for(enumerator i = 0; i < NF; ++i)
          {
            for(enumerator j = 0; j < NF; ++j)
              if( GRADF(i,j) + alpha*GRADFstab(i,j) < 0.0 )
              {
                alpha = std::min(alpha,-GRADF(i,j)/GRADFstab(i,j));
              }
          }
          std::cout << "stability" << std::endl;
          GRADFstab.Print();
          std::cout << "alpha " << alpha << std::endl;
          GRADF += GRADFstab*alpha; //consistency part
          std::cout << "consistency+stability" << std::endl;
          GRADF.Print();
          GRAD += GRADF;

          //GRAD.Print();
          k++;
        } //end of loop over faces
        GRAD = GRAD*FromDiagonalInverse(NG.data(),NF); //normalize gradients if we computed them multiple times

        //std::cout << "Gradient matrix: " << std::endl;
        //GRAD.Print();
        */
        real_array W = cell->RealArrayDV(tag_W); //access data structure for gradient matrix in mesh
        W.resize(NF*NF); //resize the structure
        std::copy(GRAD.data(),GRAD.data()+NF*NF,W.data()); //write down the gradient matrix

        cell->VariableArrayDV(tag_D).resize(NF); //resize scaling matrix D for the future use
      } //end of loop over cells

      if( m->HaveTag("FORCE") ) //Is there force on the mesh?
      {
        Tag tag_F0 = m->GetTag("FORCE"); //initial force
        assert(tag_F0.isDefined(CELL)); //assuming it was defined on cells
        tag_F = m->CreateTag("RECOMPUTED_FORCE",DATA_REAL,CELL,NONE,1); //recomputed force on face
        // compute expression H^T M H f, where H is harmonic averaging and M is inner product over faces
        for(Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell)
        {
          ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell
          enumerator NF = faces.size(), k = 0;
          Cell cK = cell->self();
          real gK = cK->Real(tag_F0), gL; //force on current cell and on adjacent cell
          rMatrix M(cK->RealArrayDF(tag_M).data(),NF,NF); //inner product matrix
          rMatrix gF(NF,1); //vector of forces on faces
          rMatrix MgF(NF,1); //vector of forces multiplied by inner product matrix
          //first average force onto faces with harmonic mean
          k = 0;
          for(ElementArray<Face>::iterator face = faces.begin(); face != faces.end(); ++face) //go over faces
          {
            real_array H = face->RealArrayDF(tag_H); //access structure for harmonic averaging
            if( !face->Boundary() ) //internal face
            {
              Cell cL = cK->Neighbour(face->self());
              real aL = (cK == face->BackCell()) ? H[3] : 1-H[3]; //harmonic coefficient for neighbour cell 
              real aK = 1 - aL; //harmonic coefficient for current cell
              gL = cL->Real(tag_F0); //retrive force on adjacent cell
              gF(k,0) = gK*aK + gL*aL; //harmonic averaging
            }
            else gF(k,0) = gK;
            k++;
          } //end of loop over faces
          MgF = M*gF; //multiply by inner product matrix
          //now apply H^T for current cell only
          real & rgK = cell->Real(tag_F); //access force in current cell
          rgK = 0.0;
          k = 0;
          for(ElementArray<Face>::iterator face = faces.begin(); face != faces.end(); ++face) //go over faces
          {
            real_array H = face->RealArrayDF(tag_H); //access structure for harmonic averaging
            if( !face->Boundary() ) //internal face
            {
              Cell cL = cell->Neighbour(face->self());
              real aL = (cK == face->BackCell()) ? H[3] : 1-H[3]; //harmonic coefficient for neighbour cell
              real aK = 1 - aL; //harmonic coefficient for current cell
              real & rgL = face->FrontCell()->RealDF(tag_F); // access force in adjacent cell
              //harmonic distribution
              rgK += MgF(k,0)*aK; // application of H^T with inner product in current cell
              rgL += MgF(k,0)*aL; // application of H^T with inner product in current cell
            }
            else rgK += MgF(k,0);
            k++;
          } //end of loop over faces
        } //end of loop over cells
        m->ExchangeData(tag_F,CELL,0); //synchronize force
      } // end of force

    } //end of initialize data

   
		
    

    { //Main loop for problem solution
      Automatizator aut(m); // declare class to help manage unknowns
      Sparse::RowMerger & merger = aut.GetMerger(); //get structure that helps matrix assembly
      dynamic_variable P(aut,aut.RegisterDynamicTag(tag_P,CELL|FACE)); //register pressure as primary unknown
      enhanced_variable calc(aut); //declare variable that helps calculating the value with variations
      aut.EnumerateDynamicTags(); //enumerate all primary variables

     
      
      Sparse::Matrix Jacobian("",aut.GetFirstIndex(),aut.GetLastIndex()); //matrix for jacobian
      Sparse::Vector Update("",aut.GetFirstIndex(),aut.GetLastIndex()); //vector for update
      Sparse::Vector Residual("",aut.GetFirstIndex(),aut.GetLastIndex()); //vector for residual
      real Residual_norm = 0.0;

      do
      {
        std::fill(Residual.Begin(),Residual.End(),0.0); //clean up the residual

        //First we need to evaluate the gradient at each cell for scaling matrix D
        for(Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell) //loop over cells
        {
          ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell
          enumerator NF = faces.size(), k = 0;
          Cell cK = cell->self();
          variable pK = P(cK), pL; //pressure for back and front cells
          rMatrix GRAD(cK->RealArrayDV(tag_W).data(),NF,NF); //Matrix for gradient
          vMatrix pF(NF,1); //vector of pressures on faces
          vMatrix FLUX(NF,1); //computed flux on faces
          for(ElementArray<Face>::iterator face = faces.begin(); face != faces.end(); ++face)
          {
            real_array H = face->RealArrayDF(tag_H); //access structure for harmonic averaging
            if( !face->Boundary() ) //internal face
            {
              Cell cL = cK->Neighbour(face->self()); //
              real aL = (cK == face->BackCell()) ? H[3] : 1-H[3]; //harmonic coefficient for neighbour cell
              real aK = 1 - aL; //harmonic coefficient for current cell
              pL = P(face->FrontCell()); //retrive force on adjacent cell
              pF(k,0) = pK*aK + pL*aL;
            }
            else if( face->HaveData(tag_P) ) //boundary condition
              pF(k,0) = P(face->self()); //get pressure on boundary
            else
              pF(k,0) = pK; //get pressure in current cell

            pF(k,0) -= pK; //substract pressure in the center to get differences
            k++;
          }
          FLUX = GRAD*pF; //fluxes on faces
          var_array D = cell->VariableArrayDV(tag_D); //access data structure on mesh for variations
          for(enumerator l = 0; l < NF; ++l) D[l] = FLUX(l,0); //copy the computed flux with variations into mesh
        } //end of loop over cells

        //Now we need to assemble and transpose nonlinear gradient matrix
        for(Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell) //loop over cells
        {
          ElementArray<Face> facesK = cell->getFaces(), facesL; //obtain faces of the cell
          enumerator NF = facesK.size(), k, l;
          Cell cK = cell->self();
          variable pK = P(cK), pL; //pressure for back and front cells
          rMatrix GRAD(cK->RealArrayDV(tag_W).data(),NF,NF); //Matrix for gradient
          vMatrix D(NF,NF); //Matrix for diagonal
          vMatrix pF(NF,1); //vector of pressures on faces
          D.Zero();
          //assemble diagonal matrix, loop through faces and access corresponding entries
          k = 0;
          for(ElementArray<Face>::iterator faceK = facesK.begin(); faceK != facesK.end(); ++faceK) //loop over faces of current cell
          {
            if( !faceK->Boundary() ) //face is on boundary
            {
              variable DK(cK->VariableArrayDV(tag_D)[k]);
              Cell cL = cK->Neighbour(faceK->self()); //retrive neighbour
              facesL = cL->getFaces(); //retrive faces of other cells
              l = 0; //count face number 
              for(ElementArray<Face>::iterator faceL = facesL.begin(); faceL != facesL.end(); ++faceL) //loop over faces of other cell
              {
                if( faceK->self() == faceL->self() ) //faces match - get flux entry
                {
                  variable DL(cL->VariableArrayDV(tag_D)[l]); //retrive flux entry as variable
                  variable absDK = sqrt(DK*DK+reg_abs) + reg_div;
                  variable absDL = sqrt(DL*DL+reg_abs) + reg_div;
                  D(k,k) = absDL/(absDK+absDL); // compute absolute value and write to D matrix
                  break;
                }
                l++; //advance the face number
              }
            }
            else
            {
              D(k,k) = 1; //should somehow use the gradient defined on boundary
            }

            real_array H = faceK->RealArrayDF(tag_H); //access structure for harmonic averaging
            if( !faceK->Boundary() ) //internal face
            {
              Cell cL = cK->Neighbour(faceK->self()); //
              real aL = (cK == faceK->BackCell()) ? H[3] : 1-H[3]; //harmonic coefficient for neighbour cell
              real aK = 1 - aL; //harmonic coefficient for current cell
              pL = P(faceK->FrontCell()); //retrive force on adjacent cell
              pF(k,0) = pK*aK + pL*aL;
            }
            else if( faceK->HaveData(tag_P) ) //boundary condition
              pF(k,0) = P(faceK->self()); //get pressure on boundary
            else
              pF(k,0) = pK; //get pressure in current cell

            pF(k,0) -= pK; //substract pressure in the center to get differences

            k++;
          }
          vMatrix DGRAD = D*GRAD; //multiply the matrix with diagonal
          rMatrix M(cK->RealArrayDV(tag_M).data(),NF,NF); //inner product matrix
          vMatrix DIVGRADHp = DGRAD.Transpose()*M*DGRAD*pF; //cell-wise divgradp
          //redistribute divgradp onto cells by applying H^T
          enumerator indK = P.Index(cK->self()),indL;
          Storage::real & ResidK = Residual[indK];
          Sparse::Row & JacK = Jacobian[indK]; //access force in current cell
          k = 0;
          for(ElementArray<Face>::iterator face = facesK.begin(); face != facesK.end(); ++face) //go over faces
          {
            real_array H = face->RealArrayDF(tag_H); //access structure for harmonic averaging
            if( !face->Boundary() ) //internal face
            {
              Cell cL = cell->Neighbour(face->self());
              real aL = (cK == face->BackCell()) ? H[3] : 1-H[3]; //harmonic coefficient for neighbour cell
              real aK = 1 - aL; //harmonic coefficient for current cell
              indL = P.Index(face->FrontCell());
              Storage::real & ResidL = Residual[indL];
              Sparse::Row & JacL = Jacobian[indL]; // access force in adjacent cell
              //harmonic distribution
              merger.Merge(JacK,1.0,JacK,-aK,DIVGRADHp(k,0).GetRow());// application of H^T with inner product in current cell
              ResidK += DIVGRADHp(k,0).GetValue()*aK;
              merger.Merge(JacL,1.0,JacL,-aL,DIVGRADHp(k,0).GetRow());// application of H^T with inner product in current cell
              ResidL += DIVGRADHp(k,0).GetValue()*aL;
            }
            else 
            {
              indL = P.Index(face->self());
              Storage::real & ResidL = Residual[indL];
              Sparse::Row & JacL = Jacobian[indL]; // access force in adjacent cell
              merger.Merge(JacK,1.0,JacK,-1.0,DIVGRADHp(k,0).GetRow());
              ResidK += DIVGRADHp(k,0).GetValue();
              merger.Merge(JacL,1.0,JacL,-1.0,DIVGRADHp(k,0).GetRow());
              ResidL += DIVGRADHp(k,0).GetValue();
            }
            k++;
          } //end of loop over faces
        }

        if( tag_F.isValid() )
        {
          for(Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell)
            Residual[P.Index(cell->self())] += cell->RealDF(tag_F);
        }
        

        Residual_norm = 0.0;
        for(Sparse::Vector::iterator it = Residual.Begin(); it != Residual.End(); ++it) Residual_norm += (*it)*(*it); //compute residual norm on current iteration

        std::cout << "Nonlinear residual: " << Residual_norm << std::endl;

        Solver S(Solver::INNER_MPTILUC);
        S.SetMatrix(Jacobian);
        S.SetParameterReal("drop_tolerance",1.0e-3);
        S.SetParameterReal("reuse_tolerance",1.0e-6);
        S.SetParameterEnum("gmres_substeps",4);
        if( S.Solve(Residual,Update) )
        {
          for(Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell)
            cell->Real(tag_P) += Update[P.Index(cell->self())];
        }
        else
        {
          std::cout << "Unable to solve: " << S.GetReason() << std::endl;
        }
      
      } while( Residual_norm > 1.0e-4 ); //check the residual norm
    } 

		delete m; //clean up the mesh
	}
	else
	{
		std::cout << argv[0] << " mesh_file" << std::endl;
	}

#if defined(USE_PARTITIONER)
	Partitioner::Finalize(); // Finalize the partitioner activity
#endif
	Solver::Finalize(); // Finalize solver and close MPI activity
	return 0;
}
