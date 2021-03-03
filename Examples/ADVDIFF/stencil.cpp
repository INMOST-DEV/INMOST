#include "stencil.h"

using namespace INMOST;

//#define HARMONIC_VERSION

//shortcuts
typedef Storage::real            real;
typedef Storage::real_array      real_array;
typedef Storage::integer         integer;
typedef Storage::reference_array ref_array;



bool predicate(real a, real b)
{
	return fabs(a-b) < 1.0e-10*(fabs(a)+fabs(b));
}

bool find_stencils(Cell cK, 
	               std::vector<Stencil> & compute, 
				   Tag tag_BC, 
				   Tag tag_K,
				   Tag tag_iT, 
				   Tag tag_iC,
				   //Tag tag_U,
				   MarkerType boundary_marker,
				   ElementType bridge_layers,
				   real regularization,
				   int max_layers)
{
	//integer cKid = cK.LocalID();
	Mesh & m = *cK.GetMeshLink(); //link to the mesh
	integer layer_wgt = -1; //weight of current layer of cells
	ElementArray<Cell> layer(&m), //currently considered layer of cells
		               next_layer(&m), //layer to be considered further
					   following_layer(&m), 
					   todo_layer(&m), //all the scheduled cells
					   done_cells(&m); //all the visited cells (for data cleanup)
	MarkerType scheduled, //cells to be visited in the following sweep
		       done, //already visited cells
			   todo; //cells to be considered in current sweep

	//remember condition number
	std::vector<real> mincond(compute.size());
	std::vector<real> minsum(compute.size());
	std::vector<integer> minwgt(compute.size());

	//fill data
	for(integer k = 0; k < (integer)compute.size(); ++k)
	{
		minwgt[k] = max_layers*3;
		mincond[k] = 1.0e+100;
		minsum[k] = 1.0e+100;
		compute[k].computed = false;
	}

	//allocate markers
	scheduled = m.CreatePrivateMarker();
	done = m.CreatePrivateMarker();
	todo = m.CreatePrivateMarker();

	//fill zero correction
	real_array correction = cK.RealArray(tag_iC), correctionN;
	correction.resize(3);
	std::fill(correction.begin(),correction.end(),0.0);

	//fill unit tensor
	real_array tensor = cK.RealArray(tag_iT), tensorN;
	tensor.resize(9);
	std::fill(tensor.begin(),tensor.end(),0.0);
	tensor[0] = tensor[4] = tensor[8] = 1.0;

	//prepare arrays
	following_layer.push_back(cK);
	done_cells.push_back(cK);
	cK.SetPrivateMarker(done); //avoid visiting myself

	
	//store all the precomputed interpolation directions
	std::vector<real> dirs; //all directions
	std::vector<real> bndrhs; //right hand sides coming out of boundary conditions
	std::vector<real> bndmlt; //multiplier for original cell
	std::vector<integer> weights; //layer number
	std::vector<HandleType> handles; //elements corresponding to each direction

	//variables used for interpolation calculation
	rMatrix xK(1,3), //current cell position
			xN(1,3), //neighbouring cell position
			yN(1,3), //projection of neighbouring cell position onto interface
			yS(1,3), //position of harmonic point
			xF(1,3), //interface position
			xD(1,3), //direction
			nF(1,3) //boundary face normal
			;
	integer num_tensors, num_tensorsN; //number of variants of interpolation tensor
	real lambdaN, //projection of co-normal onto normal in neighbour cell
		 lambdaC, //projection of co-normal onto normal in current cell
		 lambdaK,
		 rN,
		 rK,
	     rQ,
		 mult; //boundary condition multiplier
	(void)lambdaC,(void)rN,(void)rK,(void)rQ;
	//variables used for approximation
	rMatrix A(3,3), U(3,3), S(3,3), V(3,3), Sinv(3,3), I(3,3), v(1,3); //gradient matrix and SVD storage
	rMatrix coef(1,3); //row of coefficients for triplet
	real //velocity = 0.0, //velocity value 
		 cond, //current condition number of gradient matrix
		 Smin, //minimal singular value
		 Smax, //maximal singular value
		 coef_sum; //sum of coefficients

	//fill with zeroes
	Sinv.Zero();
	
	Cell cC, //current cell
		 cN, //neighbouring cell
		 cI; //element used in interpolation
	rMatrix KK(3,3), //diffusion tensor at first cell, for boundary condition calculation
			KC(3,3), //diffusion tensor in current cell, for interpolation calculation
			KN(3,3), //diffusion tensor in neighbouring cell, for interpolation calculation
			KD(3,3), //difference
			iTN(3,3), //interpolation tensor for the neighbouring cell
			iCN(1,3), //interpolation correction for the neighbouring cell
			iQ(1,3) //correction due to tensor for harmonic point
			;
	integer total_vectors; //total number of vectors
	integer combination_wgt; //weight of combination
	bool selected_combination; //is there a combination of vectors
	bool have_similar; //is there similar correction tensor
	//retrive tensor at current cell
	if( tag_K.isValid() )
		KK = rMatrix::FromTensor(cK.RealArray(tag_K).data(),cK.RealArray(tag_K).size());
	else
	{
		KK = rMatrix::Unit(3);
		KC = rMatrix::Unit(3);
		KN = rMatrix::Unit(3);
		KD.Zero();
	}

	integer C[3] = {0,1,2}; //unknowns for combinatorial search
	integer have_stencils = 0; //record how many stencils we have found
	
	
	cK.Centroid(xK.data()); //retrive current cell center

	do
	{
		//advance weight according to the layer number
		layer_wgt++;
		//First mark all the elements that I'm going to visit
		for(integer k = 0; k < (integer)following_layer.size(); ++k)
		{
			ElementArray<Element> bridge = following_layer[k].getAdjElements(bridge_layers); //retrive faces of the current cell
			for(integer q = 0; q < (int)bridge.size(); ++q)
			{
				ElementArray<Cell> cells = bridge[q].getCells(done | todo,true); //retrive cells that were not yet considered or done
				cells.SetPrivateMarker(todo); //mark cells as already done
				todo_layer.Unite(cells); //push cells into array
			}
		}
		//Make following layer as current
		layer.swap(following_layer);
		//Make todo layer as next following layer
		following_layer.swap(todo_layer);
		//Clear todo layer for the next iteration
		todo_layer.clear();
		//Visit boundary conditions only for the first layer
		bool first_layer = true;
		//Now grow layer by layer until all todo cells are done.
		//This loop is required in order to traverse all the cells
		//neighbouring current cell over bridge elements.
		//Say, we have to consider all the neighbouring cells sharing
		//edges with the current. But interpolation algorithm allows
		//to consider only elements over faces. Thus we traverse cells
		//over faces until all the cells over edge were considered.
		while(!layer.empty())
		{
			for(integer k = 0; k < (integer)layer.size(); ++k)
			{
				tensor = layer[k].RealArray(tag_iT); //retrive interpolation tensors
				correction = layer[k].RealArray(tag_iC); //retrive interpolation correctors
				num_tensors = tensor.size()/9; //total number of variants
				for(integer it_tensor = 0; it_tensor < num_tensors; ++it_tensor)
				{
					//convert into matrices
					rMatrix iT(&tensor[it_tensor*9],3,3), iC(&correction[it_tensor*3],1,3);
					ElementArray<Face> faces = layer[k].getFaces(); //retrive faces of the cell
					for(integer q = 0; q < (integer)faces.size(); ++q)
					{
						const double eps = 1.0e-30;
						//if( faces[q].Boundary() ) //this is boundary face
						if( faces[q].GetMarker(boundary_marker) ) //this is boundary face
						{
							if( first_layer )
							{
								faces[q].Centroid(xF.data()); //center of face
								faces[q].OrientedUnitNormal(layer[k],nF.data()); //normal to face
								real bcconds[3] = {0.0,1.0,0.0}; //pure neumann boundary condition
								if( tag_BC.isValid() && faces[q].HaveData(tag_BC) ) //are there boundary conditions on face?
								{
									//retrive boundary conditions
									real_array bc = faces[q].RealArray(tag_BC);
									bcconds[0] = bc[0];
									bcconds[1] = bc[1];
									bcconds[2] = bc[2];
								}
								//calculate direction corresponding to boundary condition
								xD.Zero();
								//Dirichlet part
								if( fabs(bcconds[0]) > 1.0e-12 ) xD += bcconds[0]*((xF - xK)*iT - iC);
								//Neumann part
								if( fabs(bcconds[1]) > 1.0e-12 ) 
								{
									lambdaK = nF.DotProduct(nF*KK); //project conormal onto normal
									xD += bcconds[1]*nF*(KK*lambdaK/(lambdaK+eps)+rMatrix::Unit(3)*eps/(lambdaK+eps))*iT;
								}
								//Store information about boundary condition
								dirs.insert(dirs.end(),xD.data(),xD.data()+3);
								bndrhs.push_back(bcconds[2]);
								bndmlt.push_back(bcconds[0]);
								handles.push_back(InvalidHandle());
								weights.push_back(layer_wgt);
							}
						}
						else if( !faces[q].Boundary() ) //treat neighbouring cell
						{
							cC = layer[k]; //current cell
							cN = layer[k].Neighbour(faces[q]); //neighbour cell
							if( cN.isValid() && cN.GetPrivateMarker(todo) )
							{
								//retrive tensors for interpolation
								if( tag_K.isValid() )
								{
									KC = rMatrix::FromTensor(cC.RealArray(tag_K).data(),cC.RealArray(tag_K).size());
									KN = rMatrix::FromTensor(cN.RealArray(tag_K).data(),cN.RealArray(tag_K).size());
									KD = KC-KN;
								}
								cN.Centroid(xN.data()); //neighbour center
								faces[q].Centroid(xF.data()); //center of face
								faces[q].OrientedUnitNormal(layer[k],nF.data()); //normal to face
								//projections
								lambdaC = nF.DotProduct(nF*KC);
								lambdaN = nF.DotProduct(nF*KN)+regularization; //project conormal onto normal
								lambdaK = nF.DotProduct(nF*KK)+regularization; //project conormal onto normal
								mult = 1;
								//distances
#if defined(HARMONIC_VERSION)
								rK = nF.DotProduct(xF-xK);
								rN = nF.DotProduct(xN-xF);
								//correction due to tensor
								iQ = -iC * iT.Invert();
								rQ = nF.DotProduct(iQ);
								//harmonic point
								yS = (lambdaC*rN*(xK - iQ - rQ*nF) + lambdaN*(rK+rQ)*xN + rN*(rK+rQ)*nF*KD)/(lambdaC*rN + lambdaN*(rK+rQ));
								mult = lambdaN*(rK+rQ)/(lambdaC*rN + lambdaN*(rK+rQ));
								xD = (yS - xK)*iT - iC;
								cI = cN;
#else //HARMONIC_VERSION
								if( lambdaN > 0.0 || KD.FrobeniusNorm() == 0.0 ) //neighbour cell is permiable or there is no tensor jump
#endif
								{
									if( lambdaN ) KD /= lambdaN;
									//tensor and corrector for the current interface
									iTN = nF.Transpose()*nF*KD;
									iCN = (xF-xK)*iTN;
									iTN += rMatrix::Unit(3);
									//fit into global chain
									iTN = iTN*iT; //products of tensors
									iCN = iCN*iT + iC; //sum of products with correctors
									//record tensor for neighbouring cell
									tensorN = cN.RealArray(tag_iT); //access interpolation tensor
									correctionN = cN.RealArray(tag_iC); //access interpolation correction
									//search for duplicates
									have_similar = false;
									num_tensorsN = tensorN.size()/9;
									for(integer it_tensorN = 0; it_tensorN < num_tensorsN && !have_similar; ++it_tensorN)
									{
										//compare tensor and corrector
										//TODO: floating point comparison
										if( std::equal(correctionN.data()+it_tensorN*3,correctionN.data()+(it_tensorN+1)*3,iCN.data(),predicate) &&
											std::equal(tensorN.data()+it_tensorN*9,tensorN.data()+(it_tensorN+1)*9,iTN.data(),predicate) )
											have_similar = true; //this is duplicate
									}
									if( !have_similar ) //insert new interpolation information
									{
										tensorN.insert(tensorN.end(),iTN.data(),iTN.data()+9);
										correctionN.insert(correctionN.end(),iCN.data(),iCN.data()+3);
									}
#if !defined(HARMONIC_VERSION)
									xD = (xN-xK)*iTN - iCN;
									cI = cN;
#endif //HARMONIC_VERSION
								}
#if !defined(HARMONIC_VERSION) 								
								else //cell is impermiable, using Neumann boundary condition on interface
								{
									xD = nF*(KK*lambdaK/(lambdaK+eps)+rMatrix::Unit(3)*eps/(lambdaK+eps))*iT;
									cI = cK;
									mult = 0;
								}
#endif //HARMONIC_VERSION

								//Store information about boundary condition
								dirs.insert(dirs.end(),xD.data(),xD.data()+3);
								bndrhs.push_back(0.0);
								bndmlt.push_back(mult);
								handles.push_back(cI.GetHandle());
								weights.push_back(layer_wgt);
								//schedule permiable neighbouring cell for consideration
								if( !cN->GetPrivateMarker(scheduled) && mult > 0 )
								{
									next_layer.push_back(cN);
									cN.SetPrivateMarker(scheduled);
								}
							} //is neighbour valid?
						} // if,else for boundary face
					} //loop over faces of current cell
				} //loop over interpolation tensors
			} //loop over cells in the layer
			//clean-up markers
			next_layer.RemPrivateMarker(scheduled);
			next_layer.RemPrivateMarker(todo);
			//identify cells that were done
			next_layer.SetPrivateMarker(done);
			//cells that were already considered
			done_cells.insert(done_cells.end(),next_layer.begin(),next_layer.end());
			//select next layer as current
			layer.swap(next_layer);
			//clear next layer for refill at the next sweep
			next_layer.clear();
			//do not visit boundary conditions of subsequent cells
			first_layer = false; 
		} //all the cells over bridge elements were considered
		total_vectors = (integer)dirs.size()/3; //current set of elements
		//start considering combinations from the beginning
		C[0] = 0;
		C[1] = 1;
		C[2] = 2;
		do
		{
			//gradient matrix for current combinatorial selection
			combination_wgt = 0;
			for(integer q = 0; q < 3; ++q)
			{
				A(q,0) = dirs[C[q]*3+0];
				A(q,1) = dirs[C[q]*3+1];
				A(q,2) = dirs[C[q]*3+2];
				combination_wgt += weights[C[q]];
			}
			//calculate inverse with singular value decomposition
			if( A.SVD(U,S,V) )
			{
				Smax = Smin = S(0,0);
				//invert S matrix
				for(integer q = 0; q < 3; ++q)
				{
					if( S(q,q) > 0.0 )
					{
						Smin = S(q,q);
						Sinv(q,q) = 1.0/S(q,q);
					}
					else Sinv(q,q) = 0.0;
				}
				//condition number
				cond = Smax/Smin;
				//this should be unit
				I = V*Sinv*U.Transpose()*U*S*V.Transpose();
				if( false )
				{
					std::cout << "C: " << C[0] << "," << C[1] << "," << C[2] << std::endl;
					std::cout << "A: " << std::endl; 
					A.Print();
					std::cout << "U: " << std::endl; 
					U.Print();
					std::cout << "S: " << std::endl; 
					S.Print();
					std::cout << "V: " << std::endl; 
					V.Print();
					std::cout << "cond: " << cond << std::endl;
					std::cout << "I: " << std::endl; 
					I.Print();
					for(integer l = 0; l < (integer)compute.size(); ++l)
					{
						std::cout << "vector[" << l << "]: " << compute[l].v[0] << "," << compute[l].v[1] << "," << compute[l].v[2] << std::endl;
						std::cout << "coefs: " << std::endl;
						(rMatrix(compute[l].v,1,3)*V*Sinv*U.Transpose()).Print();
					}
				}
				
				//check factorization
				//if( (I - rMatrix::Unit(3)).FrobeniusNorm() < 1.0e-8 )
				{
					for(integer l = 0; l < (integer)compute.size(); ++l) //if( cond < mincond[l] )//&& combination_wgt < minwgt[l] )
					{
						//current vector
						v = rMatrix(compute[l].v,1,3);
						//compute coefficients
						coef = v*V*Sinv*U.Transpose();
						//check coefficients
						if( (v-coef*A).FrobeniusNorm() < 1.0e-8 )
						{
							//check positivity
							bool nonnegative = true;
							coef_sum = 0.0;
							for(integer q = 0; q < 3; ++q)
								coef_sum += fabs(coef(0,q));
								
							
							for(integer q = 0; q < 3; ++q)
							{
								if( coef(0,q) < 0.0 )
								{
									if( coef(0,q)/coef_sum < -1.0e-5 )
										nonnegative = false;
									//else //truncate small values
									//	coef(0,q) = 0.0;
								}
							}
							//it's positive and condition number is improved
							if( nonnegative && coef_sum < minsum[l] )
							{
								/*
								if( (I - rMatrix::Unit(3)).FrobeniusNorm() > 1.0e-8 )
								{
#pragma omp critical
									{
										std::cout << "deviation from unit: " << std::endl;
										(I - rMatrix::Unit(3)).Print();
										std::cout << "deviation from conormal: " << std::endl;
										(v-coef*A).Print();
									}
								}
								*/
								//remember stencil
								if( !compute[l].computed )
								{
									have_stencils++;
									compute[l].computed = 1;
								}
								//resize variable-sized arrays (do not touch fixed size)
								if( compute[l].coefs.size() != 4 ) compute[l].coefs.resize(4);
								if( compute[l].elems.size() != 4 ) compute[l].elems.resize(4);
								compute[l].coefs[0] = coef(0,0)*compute[l].sign*bndmlt[C[0]];
								compute[l].coefs[1] = coef(0,1)*compute[l].sign*bndmlt[C[1]];
								compute[l].coefs[2] = coef(0,2)*compute[l].sign*bndmlt[C[2]];
								compute[l].coefs[3] = -compute[l].sign*(bndmlt[C[0]]*coef(0,0) + bndmlt[C[1]]*coef(0,1) + bndmlt[C[2]]*coef(0,2));
								compute[l].elems.at(0) = handles[C[0]];
								compute[l].elems.at(1) = handles[C[1]];
								compute[l].elems.at(2) = handles[C[2]];
								compute[l].elems.at(3) = cK.GetHandle();
								compute[l].nonnegative = nonnegative;
								*compute[l].rhs = compute[l].sign*(bndrhs[C[0]]*coef(0,0) + bndrhs[C[1]]*coef(0,1) + bndrhs[C[2]]*coef(0,2));
								mincond[l] = cond;
								minsum[l] = coef_sum;
								minwgt[l] = combination_wgt;
							} //if nonnegative
						} //check coefficients
					} //loop over stencils
				} //test for decomposition
			}
			else std::cout << __FILE__ << ":" << __LINE__ << " svd failed" << std::endl;
			//proceed to the next combination of vectors
			selected_combination = false;
			for(integer q = 2; q >= 0 && !selected_combination; q--)
			{
				while( C[q] < total_vectors && !selected_combination )
				{
					C[q]++;
					selected_combination = C[q] < total_vectors;
					for(integer j = q+1; j < 3; ++j)
					{
						C[j] = C[j-1]+1;
						if( C[j] >= total_vectors ) selected_combination = false;
					}
				}
			}
		}
		while(selected_combination); //there is a combination of vectors
	} while( have_stencils != (integer)compute.size() && layer_wgt < max_layers );
	//unmark all visited cells
	done_cells.RemPrivateMarker(done);
	//clean up used memory
	for(integer k = 0; k < (integer)done_cells.size(); ++k)
	{
		done_cells[k].RealArray(tag_iT).clear();
		done_cells[k].RealArray(tag_iC).clear();
	}
	//deallocate markers
	m.ReleasePrivateMarker(done);
	m.ReleasePrivateMarker(todo);
	m.ReleasePrivateMarker(scheduled);
	return have_stencils == (integer)compute.size();
}
