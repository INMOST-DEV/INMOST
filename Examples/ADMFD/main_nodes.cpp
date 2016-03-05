#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "inmost.h"
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
    Tag tag_W;  // Gradient matrix acting on harmonic points on faces and returning gradient on faces
    Tag tag_S;  // Stabilization part for DIVKGRAD
    Tag tag_D;  // Entries for scaling matrix D
    Tag tag_LF; // Coefficients of linearly computed fluxes on faces, 2 of them per internal face

    if( m->GetProcessorsNumber() > 1 ) //skip for one processor job
    { // Exchange ghost cells
		  ttt = Timer();
		  m->ExchangeGhost(1,FACE); // Produce layer of ghost cells
		  BARRIER
		  if( m->GetProcessorRank() == 0 ) std::cout << "Exchange ghost: " << Timer()-ttt << std::endl;
    }
    
    { //initialize data
		  if( m->HaveTag("PERM") ) // is diffusion tensor already defined on the mesh? (PERM from permeability)
        tag_K = m->GetTag("PERM"); // get the diffusion tensor
      
      if( !tag_K.isValid() || !tag_K.isDefined(CELL) ) // diffusion tensor was not initialized or was not defined on cells.
      {
        tag_K = m->CreateTag("PERM",DATA_REAL,CELL,NONE,6); // create a new tag for symmetric diffusion tensor K
        for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) ) // loop over mesh cells
        {
          Cell cell = m->CellByLocalID(q);
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
        tag_P = m->CreateTag("PRESSURE",DATA_REAL,CELL,NONE,1); // Create a new tag for the pressure
        for(Mesh::iteratorElement e = m->BeginElement(CELL); e != m->EndElement(); ++e) //Loop over mesh cells
          e->Real(tag_P) = (rand()*1.0)/(RAND_MAX*1.0); // Prescribe random value in [0,1]
      }

      if( !tag_P.isDefined(NODE) )
      {
        tag_P = m->CreateTag("PRESSURE",DATA_REAL,NODE,NONE,1);
        for(Mesh::iteratorElement e = m->BeginElement(NODE); e != m->EndElement(); ++e) //Loop over mesh cells
          e->Real(tag_P) = (rand()*1.0)/(RAND_MAX*1.0); // Prescribe random value in [0,1]
      }


      
      if( m->HaveTag("BOUNDARY_CONDITION") ) //Is there boundary condition on the mesh?
      {
        tag_BC = m->GetTag("BOUNDARY_CONDITION");
        
        //initialize unknowns at boundary
      }
      m->ExchangeData(tag_P,CELL|NODE,0); //Synchronize initial solution with boundary unknowns
      //run in a loop to identify boundary pressures so that they enter as primary variables
      tag_M = m->CreateTag("INNER_PRODUCT",DATA_REAL,CELL,NONE);
      //assemble inner product matrix M acting on faces for each cell
      for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
      {
        Cell cell = m->CellByLocalID(q);
        real xP[3]; //center of the cell
        real nF[3]; //normal to the face
        real yF[3]; //center of the face
        real aF; //area of the face
        real vP = cell->Volume(); //volume of the cell
        cell->Centroid(xP); //obtain cell center
        ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell
        int NF = (int)faces.size(); //number of faces;
        rMatrix IP(NF,NF), N(NF,3), R(NF,3); //matrices for inner product
        rMatrix K = rMatrix::FromTensor(cell->RealArrayDF(tag_K).data(),cell->RealArrayDF(tag_K).size()); //get permeability for the cell
        //assemble matrices R and N
        //follows chapter 5.1.4
        //in the book "Mimetic Finite Difference Method for Elliptic Problems" 
        //by Lourenco et al
        for(int k = 0; k < NF; ++k)
        {
          aF = faces[k]->Area();
          faces[k]->Centroid(yF); //point on face
          faces[k]->OrientedUnitNormal(cell->self(),nF);
          // assemble R matrix, formula (5.29)
          R(k,0) = (yF[0]-xP[0])*aF;
          R(k,1) = (yF[1]-xP[1])*aF;
          R(k,2) = (yF[2]-xP[2])*aF;
          // assemble N marix, formula (5.25)
          //N(k,0) = nF[0]*aF;
          //N(k,1) = nF[1]*aF;
          //N(k,2) = nF[2]*aF;
          rMatrix nK = rMatrix::FromVector(nF,3).Transpose()*K;
          N(k,0) = nK(0,0);
          N(k,1) = nK(0,1);
          N(k,2) = nK(0,2);
        } //end of loop over faces
        // formula (5.31)
        IP = R*(R.Transpose()*N).Invert().first*R.Transpose(); // Consistency part
        // formula (5.32)
        IP += (rMatrix::Unit(NF) - N*(N.Transpose()*N).Invert().first*N.Transpose())*(2.0/(static_cast<real>(NF)*vP)*(R*K.Invert().first*R.Transpose()).Trace()); //Stability part
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
      tag_S = m->CreateTag("STAB",DATA_REAL,CELL,NONE);
      //Assemble gradient matrix W on cells
      Tag node_index = m->CreateTag("TMP_NODE_INDEX",DATA_INTEGER,NODE,NONE,1);
      for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
      {
        Cell cell = m->CellByLocalID(q);
        real xP[3]; //center of the cell
        real yN[3]; //center of the node
        real nF[3]; //normal to the face
        real aF; //area of the face
        real vP = cell->Volume(); //volume of the cell
        cell->Centroid(xP);
        ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell
        ElementArray<Node> nodes = cell->getNodes(); //obtain faces of the cell
        int NF = (int)faces.size(); //number of faces;
        int NN = (int)nodes.size(); //number of faces;
        rMatrix K = rMatrix::FromTensor(cell->RealArrayDF(tag_K).data(),cell->RealArrayDF(tag_K).size()); //get permeability for the cell
        //rMatrix U,S,V;
        //K0.SVD(U,S,V);
        //for(int k = 0; k < 3; ++k) S(k,k) = sqrt(S(k,k));
        //rMatrix K = U*S*V;
        rMatrix A(NF,NN); //Averaging matrix
        rMatrix GRAD(NF,NN), NK(NF,3), R(NN,3), STAB(NN,NN); //big gradient matrix, co-normals, directions
        GRAD.Zero();
        for(int k = 0; k < NF; ++k) //loop over faces
        {
          aF = faces[k]->Area();
          faces[k]->OrientedUnitNormal(cell->self(),nF);
          // assemble matrix of co-normals 
          rMatrix nK = rMatrix::FromVector(nF,3).Transpose()*K;
          NK(k,0) = nK(0,0);
          NK(k,1) = nK(0,1);
          NK(k,2) = nK(0,2);
        } //end of loop over faces
        for(int k = 0; k < NN; ++k)
        {
          nodes[k]->Centroid(yN);
          nodes[k]->Integer(node_index) = k;
          // assemble matrix of directions
          R(k,0) = (yN[0]-xP[0]);
          R(k,1) = (yN[1]-xP[1]);
          R(k,2) = (yN[2]-xP[2]);
        }
        A.Zero();
        for(int k = 0; k < NF; ++k)
        {
            ElementArray<Node> face_nodes = faces[k].getNodes();
            real fraction = faces[k].Area() / static_cast<real>(face_nodes.size());
            for(int l = 0; l < (int)face_nodes.size(); ++l)
                A(k,face_nodes[l].Integer(node_index)) += fraction;
        }
        //std::cout << "A:" << std::endl;
        //A.Print();
        GRAD = NK*(NK.Transpose()*A*R).Invert().first*NK.Transpose()*A; //stability part
        //GRAD += (rMatrix::Unit(NF) - R*(R.Transpose()*R).Invert().first*R.Transpose())*(1.0/(static_cast<real>(NF)*vP)*(NK*K.Invert().first*NK.Transpose()).Trace());
        GRAD += (rMatrix::Unit(NF) - A*R*((A*R).Transpose()*A*R).Invert().first*(A*R).Transpose())*A*(2.0/(static_cast<real>(NF)*vP)*(NK*K.Invert(true).first*NK.Transpose()).Trace());
        real_array W = cell->RealArrayDV(tag_W); //access data structure for gradient matrix in mesh
        W.resize(NF*NN); //resize the structure
        std::copy(GRAD.data(),GRAD.data()+NF*NN,W.data()); //write down the gradient matrix
        cell->VariableArrayDV(tag_D).resize(NF); //resize scaling matrix D for the future use

        //now compute stabilization part
        /*
        for(int k = 0; k < NN; ++k)
        {
            aF = 0;
            ElementArray<Face> node_faces = nodes[k].getFaces();
            for(int q = 0; q < (int)node_faces.size(); ++q)
                aF += node_faces[q].Area() / static_cast<real>(node_faces[q].nbAdjElements(NODE));
            R(k,0) *= aF;
            R(k,1) *= aF;
            R(k,2) *= aF;
        }
        */
        STAB = (rMatrix::Unit(NN) - R*(R.Transpose()*R).Invert().first*R.Transpose())*2.0/(static_cast<real>(NN)*vP);
        real_array S = cell->RealArrayDV(tag_S);
        S.resize(NN*NN);
        std::copy(STAB.data(),STAB.data()+NN*NN,S.data());
      } //end of loop over cells

      if( m->HaveTag("FORCE") ) //Is there force on the mesh?
      {
        tag_F = m->GetTag("FORCE"); //initial force
        assert(tag_F.isDefined(CELL)); //assuming it was defined on cells
      } // end of force

      tag_LF = m->CreateTag("LINEAR_FLUX",DATA_VARIABLE,FACE,NONE,2);

    } //end of initialize data

   
		
    

    { //Main loop for problem solution
      Automatizator aut(m); // declare class to help manage unknowns
      Automatizator::MakeCurrent(&aut);
      dynamic_variable P(aut,aut.RegisterDynamicTag(tag_P,CELL|NODE)); //register pressure as primary unknown
      variable calc; //declare variable that helps calculating the value with variations
      aut.EnumerateDynamicTags(); //enumerate all primary variables

     
      Residual R("",aut.GetFirstIndex(),aut.GetLastIndex());
      Sparse::Vector Update  ("",aut.GetFirstIndex(),aut.GetLastIndex()); //vector for update
      /*
      {//Annotate matrix
        for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
        {
          Cell cell = m->CellByLocalID(q);
          R.GetJacobian().Annotation(P.Index(cell)) = "Cell-centered pressure value";
        }
        for( int q = 0; q < m->FaceLastLocalID(); ++q ) if( m->isValidFace(q) )
        {
          Face face = m->FaceByLocalID(q);
          if( tag_BC.isValid() && face.HaveData(tag_BC) )
            R.GetJacobian().Annotation(P.Index(face)) = "Pressure guided by boundary condition";
          else
            R.GetJacobian().Annotation(P.Index(face)) = "Interface pressure";
        }
      }
      */
      int iter = 0;
      do
      {
        R.Clear(); //clean up the residual
        //First we need to evaluate the gradient at each cell for scaling matrix D
        for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) ) //loop over cells
        {
          Cell cell = m->CellByLocalID(q);
          ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell
          ElementArray<Node> nodes = cell->getNodes(); //obtain faces of the cell
          int NF = (int)faces.size();
          int NN = (int)nodes.size();
          rMatrix GRAD(cell->RealArrayDV(tag_W).data(),NF,NN); //Matrix for gradient
          vMatrix pF(NN,1); //vector of pressure differences on faces
          vMatrix FLUX(NF,1); //computed flux on faces
          for(int k = 0; k < NN; ++k)
            pF(k,0) = P(nodes[k]) - P(cell);
          FLUX = GRAD*pF; //fluxes on faces
          for(int k = 0; k < NF; ++k) //copy the computed flux value with variations into mesh
            faces[k]->VariableArray(tag_LF)[(faces[k]->BackCell() == cell)? 0 : 1] = FLUX(k,0);
        } //end of loop over cells

        //Now we need to assemble and transpose nonlinear gradient matrix
        for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) ) //loop over cells
        {
          const real eps1 = 1.0e-7;
          const real eps2 = 1.0e-3;
          Cell cell = m->CellByLocalID(q);
          ElementArray<Face> faces = cell->getFaces(); //obtain faces of the cell
          ElementArray<Node> nodes = cell->getNodes(); //obtain faces of the cell
          int NF = (int)faces.size();
          int NN = (int)nodes.size();
          rMatrix K = rMatrix::FromTensor(cell->RealArrayDF(tag_K).data(),cell->RealArrayDF(tag_K).size()); //get permeability for the cell
          rMatrix GRAD(cell->RealArrayDV(tag_W).data(),NF,NN); //Matrix for gradient
          rMatrix STAB(cell->RealArrayDV(tag_S).data(),NN,NN);
          rMatrix M(cell->RealArrayDV(tag_M).data(),NF,NF); //inner product matrix
          vMatrix D(NF,NF); //Matrix for diagonal
          vMatrix FLUX(NF,1); //computed flux on faces
          vMatrix pF(NN,1); //vector of pressure differences on faces
          D.Zero();
          for(int k = 0; k < NN; ++k)
            pF(k,0) = P(nodes[k]) - P(cell);
          //assemble diagonal matrix, loop through faces and access corresponding entries
          for(int k = 0; k < NF; ++k) //loop over faces of current cell
          {
            var_array LF = faces[k]->VariableArray(tag_LF);
            variable & u = LF[(faces[k]->BackCell() == cell) ? 0 : 1]; //my flux value
            if( faces[k].Boundary() ) 
            {
              D(k,k) = 1.0; // no flux balancing on the boundary
              FLUX(k,0) = u; //restore matrix of fluxes
            }
            else
            {
              variable & v = LF[(faces[k]->BackCell() == cell) ? 1 : 0]; //neighbour flux value
              v = -v;
              //single flux definition
              FLUX(k,0) = u;

              //Van Albada
              variable mu = (u*v+v*v+eps1)/(u*u+v*v+2*eps1);
              variable sgn = 1;
              //Van Leer
              //variable mu = soft_fabs(v,eps1)/(soft_fabs(u,eps1)+soft_fabs(v,eps1));
              //variable sgn = 1+soft_sign(u*v,eps1);
              D(k,k) = variation(mu*sgn,0.0);//+eps2;
              //D(k,k) = 1;
              //if( fabs(D(k,k)) < 1.0e-3 ) std::cout << "diagonal multiplier is small: " << get_value(D(k,k)) << " u: " << get_value(u) << " v: " << get_value(v) << " f: " << get_value(D(k,k)*u) << std::endl;
              //dual flux definition
              //FLUX(k,0) = u*D(k,k);
            }
          }
          vMatrix DIVKGRAD = -(D*GRAD).Transpose()*M*(D*GRAD); //cell-wise div
          DIVKGRAD += STAB / DIVKGRAD.Trace();
          /*
          {
              rMatrix check = DIVKGRAD, U,S,V;
              check.SVD(U,S,V);
              int rank = 0;
              for(int k = 0; k < (int)S.Rows(); ++k)
                  if( fabs(S(k,k)) > 1.0e-5 )
                      rank++;
              if( rank != NN ) std::cout << "Matrix is not full rank: " << rank << " should be " << NN << std::endl;
          }
          */
          /*
          std::cout << "D" << std::endl;
          D.Print();
          std::cout << "GRAD" << std::endl;
          GRAD.Print();
          std::cout << "M" << std::endl;
          M.Print();
          std::cout << "DIV" << std::endl;
          DIV.Print();
          */
          vMatrix DIVKGRADp = DIVKGRAD*pF;//*cell->Volume();

          

          //std::cout << "DIVKGRADp" << std::endl;
          //DIVKGRAD.Print();

          for(int k = 0; k < NN; ++k) //loop over faces of current cell
            R[P.Index(cell)] -= DIVKGRADp(k,0);
          for(int k = 0; k < NN; ++k) //loop over faces of current cell
          {
            int index = P.Index(nodes[k]);
            if( tag_BC.isValid() && nodes[k].HaveData(tag_BC) )
            {
              real_array BC = nodes[k].RealArray(tag_BC);
              R[index] += BC[0]*P(nodes[k]) + BC[1]*DIVKGRADp(k,0) - BC[2];
            }
            else
              R[index] += DIVKGRADp(k,0);
          }
        }

        if( tag_F.isValid() )
        {
          for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
          {
            Cell cell = m->CellByLocalID(q);
            if( cell->HaveData(tag_F) ) R[P.Index(cell)] += cell->Real(tag_F);//*cell->Volume();
          }
        }

        R.Rescale();
        
        //R.GetJacobian().Save("jacobian.mtx");
        //R.GetResidual().Save("residual.mtx");


        std::cout << "Nonlinear residual: " << R.Norm() << std::endl;

        if( R.Norm() < 1.0e-4 ) break;

        //Solver S(Solver::INNER_MPTILU2);
        Solver S(Solver::INNER_MPTILUC);
        S.SetMatrix(R.GetJacobian());
        S.SetParameterReal("relative_tolerance", 1.0e-11);
        S.SetParameterReal("absolute_tolerance", 1.0e-9);
        S.SetParameterReal("drop_tolerance", 1.0e-2);
        S.SetParameterReal("reuse_tolerance", 1.0e-4);
        S.SetParameterEnum("gmres_substeps", 2);
                
        //std::fill(Update.Begin(),Update.End(),0.0);
        bool flag = S.Solve(R.GetResidual(),Update);
        //if(  )
        {
          for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
          {
            Cell cell = m->CellByLocalID(q);
            cell->Real(tag_P) -= Update[P.Index(cell)];
          }
          for( int q = 0; q < m->NodeLastLocalID(); ++q ) if( m->isValidNode(q) )
          {
            Node node = m->NodeByLocalID(q);
            node->Real(tag_P) -= Update[P.Index(node)];
          }
          {
            std::stringstream file;
            file << "iter" << iter << ".vtk";
            m->Save(file.str());
          }
        }
        if( !flag )
        {
          std::cout << "Unable to solve: " << S.GetReason() << std::endl;
          break;
        }
        iter++;
      } while( R.Norm() > 1.0e-4 && iter < 20 ); //check the residual norm
    } 
    
    if( m->HaveTag("REFERENCE_SOLUTION") )
    {
      Tag tag_E = m->CreateTag("ERRROR",DATA_REAL,CELL,NONE,1);
      Tag tag_R = m->GetTag("REFERENCE_SOLUTION");
      real C, L2, volume;
      C = L2 = volume = 0.0;
      for( int q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
      {
        Cell cell = m->CellByLocalID(q);
        real err = cell->Real(tag_P) - cell->Real(tag_R);
        real vol = cell->Volume();
        if( C < fabs(err) ) C = fabs(err);
        L2 += err*err*vol;
        volume += vol;
        cell->Real(tag_E) = err;
      }
      L2 = sqrt(L2/volume);
      std::cout << "Error on cells, C-norm " << C << " L2-norm " << L2 << std::endl;
      C = L2 = volume = 0.0;
      if( tag_R.isDefined(NODE) )
      {
        tag_E = m->CreateTag("ERRROR",DATA_REAL,NODE,NONE,1);
        for( int q = 0; q < m->NodeLastLocalID(); ++q ) if( m->isValidNode(q) )
        {
          Node node = m->NodeByLocalID(q);
          real err = node->Real(tag_P) - node->Real(tag_R);
          real vol = 0;
          ElementArray<Cell> cells = node->getCells();
          for(int k = 0; k < (int)cells.size(); ++k)
              vol += cells[k].Volume();
          vol /= static_cast<real>(cells.size());
          if( C < fabs(err) ) C = fabs(err);
          L2 += err*err*vol;
          volume += vol;
          node->Real(tag_E) = err;
        }
        L2 = sqrt(L2/volume);
        std::cout << "Error on nodes, C-norm " << C << " L2-norm " << L2 << std::endl;
      }
      else std::cout << "Reference solution was not defined on nodes" << std::endl;
    }

    m->Save("out.gmv");
    m->Save("out.vtk");

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