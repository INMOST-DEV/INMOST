#include "conv_diff.h"
#include "stencil.h"


using namespace INMOST;




//shortcuts
typedef Storage::bulk            bulk;
typedef Storage::real            real;
typedef Storage::integer         integer;
typedef Storage::enumerator      enumerator;
typedef Storage::real_array      real_array;
typedef Storage::var_array       var_array;
typedef Storage::reference_array ref_array;

//todo
extern bool split_diffusion;
extern double degenerate_diffusion_regularization;
extern int max_layers; //maximum number of layers to be considered in interpolation
extern ElementType bridge_layers; //used as bridge elements for layers


ConvectionDiffusion::ConvectionDiffusion(Mesh * _m, Tag _tag_U, Tag _tag_K, Tag _tag_BC, MarkerType _boundary_face, bool correct_convection, bool correct_diffusion) 
		: m(_m), tag_U(_tag_U), tag_K(_tag_K), tag_BC(_tag_BC), boundary_face(_boundary_face)
{
	perform_correction_diffusion = correct_convection;
	perform_correction_convection = correct_diffusion;
	if( tag_U.isValid() )
	{
		//4 entries - triplet for upwind corrector, including the cell
		//back cell corrector
		tag_CONV_CB = m->CreateTag("CONV_BACK_COEFS"   ,DATA_REAL     ,FACE,NONE,4);
		tag_CONV_EB = m->CreateTag("CONV_BACK_ELEMS"   ,DATA_REFERENCE,FACE,NONE,4);
		tag_CONV_RB = m->CreateTag("CONV_BACK_RHS"     ,DATA_REAL     ,FACE,NONE,1);
		//front cell corrector
		tag_CONV_CF = m->CreateTag("CONV_FRONT_COEFS"  ,DATA_REAL     ,FACE,NONE,4);
		tag_CONV_EF = m->CreateTag("CONV_FRONT_ELEMS"  ,DATA_REFERENCE,FACE,NONE,4);
		tag_CONV_RF = m->CreateTag("CONV_FRONT_RHS"    ,DATA_REAL     ,FACE,NONE,1);
		//upstream
		tag_CONV_CU = m->CreateTag("CONV_UPSTREAM_COEF",DATA_REAL     ,FACE,NONE,1);
		tag_CONV_EU = m->CreateTag("CONV_UPSTREAM_ELEM",DATA_REFERENCE,FACE,NONE,1);
		tag_CONV_RU = m->CreateTag("CONV_UPSTREAM_RHS" ,DATA_REAL     ,FACE,NONE,1);
		tag_CONV_VU = m->CreateTag("CONV_UPSTREAM_CORR",DATA_REAL     ,FACE,NONE,6);
		//flag
		tag_CONV_F = m->CreateTag("CONV_FLAG"          ,DATA_BULK     ,FACE,NONE,1);
	}
	if( tag_K.isValid() )
	{
		//4 entries - triplet for transversal corrector, including the cell
		//back cell corrector
		tag_DIFF_CB = m->CreateTag("DIFF_BACK_COEFS"   ,DATA_REAL     ,FACE,NONE,4);
		tag_DIFF_EB = m->CreateTag("DIFF_BACK_ELEMS"   ,DATA_REFERENCE,FACE,NONE,4);
		tag_DIFF_RB = m->CreateTag("DIFF_BACK_RHS"     ,DATA_REAL     ,FACE,NONE,1);
		//front cell corrector
		tag_DIFF_CF = m->CreateTag("DIFF_FRONT_COEFS"  ,DATA_REAL     ,FACE,NONE,4);
		tag_DIFF_EF = m->CreateTag("DIFF_FRONT_ELEMS"  ,DATA_REFERENCE,FACE,NONE,4);
		tag_DIFF_RF = m->CreateTag("DIFF_FRONT_RHS"    ,DATA_REAL     ,FACE,NONE,1);
		//two-point transmissibility
		tag_DIFF_CT = m->CreateTag("DIFF_TP_COEF"      ,DATA_REAL     ,FACE,NONE,1);
		tag_DIFF_RT = m->CreateTag("DIFF_TP_RHS"       ,DATA_REAL     ,FACE,NONE,1);
		tag_DIFF_VT = m->CreateTag("DIFF_TP_CORR"      ,DATA_REAL     ,FACE,NONE,6);
		//flag
		tag_DIFF_F = m->CreateTag("DIFF_FLAG"          ,DATA_BULK     ,FACE,NONE,1);
	}

	build_flux = 0;
		
	if( m->GetProcessorsNumber() > 1 )
	{
		build_flux = m->CreateMarker();
#if defined(USE_OMP)
#pragma omp parallel for
#endif
		for(integer q = 0; q < m->FaceLastLocalID(); ++q ) if( m->isValidFace(q) )
		{
			Face fKL = m->FaceByLocalID(q);
			Element::Status stat = fKL->BackCell()->GetStatus();
			if( fKL->FrontCell().isValid() ) stat |= fKL->FrontCell()->GetStatus();
			if( stat & (Element::Shared | Element::Owned) )
				fKL->SetMarker(build_flux);
		}	
	}

	initialization_failure = false; //there is a failure during intialization
	//Precompute two-point part, transversal correction direction,
	// upstream and upstream correction directions and flags
#if defined(USE_OMP)
#pragma omp parallel
#endif
	{
		//private memory
		rMatrix	xK(1,3), //center of cell K
				yK(1,3), //projection of xK onto interface
				xL(1,3), //center of cell L
				yL(1,3), //projection of xL onto interface
				xKL(1,3), //center of face common to K and L
				nKL(1,3), //normal vector to face
				KKn(1,3), //co-normal at cell K
				KLn(1,3), //co-normal at cell L
				v(1,3), //vector for correction
				uK(1,3), //vector for upwind correction in back cell
				uL(1,3), //vector for upwind correction in front cell
				r(1,3), //path to reach upstream concentration from downstream concentration
				r0(1,3),
				KK(3,3), //tenosr at cell K
				KL(3,3), //tensor at cell L
				KD(3,3), //difference of tensors
				gammaK(1,3), //transversal part of co-normal at cell K
				gammaL(1,3),  //transversal part of co-normal at cell L
				gamma(1,3), //common transversal part for flux
				iT(3,3), //heterogeneous interpolation tensor
				iC(1,3) //heterogeneous interpolation correction
				;
		real	A, //area of the face
				U, //normal component of the velocity
				C, //coefficient for upstream cell
				T = 0, //two-point transmissibility
				R = 0, //right hand side
				dK = 0, //distance from center to interface at cell K
				dL = 0, //distance from center to interface at cell L
				lambdaK = 0, //projection of co-normal onto normal at cell K
				lambdaL = 0, //projection of co-normal onto normal at cell L
				div //divider
				;
		const real eps = degenerate_diffusion_regularization;
		Cell cK, cL, cU;
		Face fKL;
		bulk flag_DIFF = 0, flag_CONV = 0;
		KK.Zero();
		KL.Zero();
		KD.Zero();
		U = 0.0;
#if defined(USE_OMP)
#pragma omp for
#endif
		for(integer q = 0; q < m->FaceLastLocalID(); ++q ) if( m->isValidFace(q) )
		{
			fKL = m->FaceByLocalID(q);
			if( !BuildFlux(fKL) ) continue;

			fKL.Centroid(xKL.data());
			fKL.UnitNormal(nKL.data());
			A = fKL.Area();
			if( tag_U.isValid() ) U = fKL.Real(tag_U);
					
			cK = fKL.BackCell();
			cL = fKL.FrontCell();
			assert(cK.isValid());
					
			cK.Centroid(xK.data());
			dK = nKL.DotProduct(xKL-xK);
			yK = xK + dK*nKL;
			if( tag_K.isValid() )
				KK = rMatrix::FromTensor(cK.RealArray(tag_K).data(),cK.RealArray(tag_K).size());//.Transpose();
			KKn = nKL*KK;
			lambdaK = nKL.DotProduct(KKn);
			if( lambdaK < 0 ) std::cout << __FILE__ << ":" << __LINE__ << " lambdaK=" << lambdaK << std::endl;
			gammaK = KKn - lambdaK*nKL;
					
			//Diffusion part
			uK.Zero();
			uL.Zero();
			if( cL.isValid() ) //internal, both cells are present
			{
				cL.Centroid(xL.data());
				dL = nKL.DotProduct(xL-xKL);
				yL = xL - dL*nKL;
				if( tag_K.isValid() )
					KL = rMatrix::FromTensor(cL.RealArray(tag_K).data(),cL.RealArray(tag_K).size());//.Transpose();
				KLn = nKL*KL;
				lambdaL = nKL.DotProduct(KLn);
				if( lambdaL < 0 ) std::cout << __FILE__ << ":" << __LINE__ << " lambdaL=" << lambdaL << std::endl;
				gammaL = KLn - lambdaL*nKL;
						
						
				//don't forget the area!
				R = 0;
				if( split_diffusion )
				{
					div = (dK*lambdaL+dL*lambdaK);
					if( div > 0 )
					{
						T = A * lambdaK*lambdaL/div;
						gamma = A * (lambdaK*lambdaL*(yK-yL) + lambdaL*dK*gammaK + lambdaK*dL*gammaL)/div;
						uK = gamma;
						uL = gamma;
					}
					else if( div < 0 )
					{
						T = -A * lambdaK*lambdaL/div;
						gamma = A * (lambdaK*lambdaL*(yK-yL) + lambdaL*dK*gammaK + lambdaK*dL*gammaL)/div;
						uK = 2*A*KKn - gamma;
						uL = 2*A*KLn - gamma;
					}
					else
					{
						uK = A * KKn;
						uL = A * KLn;
						T = 0;
					}
				}
				else //un-splitted diffusion
				{
					uK = A * KKn;
					uL = A * KLn;
					T = 0;
				}
				
				if( T < 0 ) std::cout << __FILE__ << ":" << __LINE__ << " T=" << T << std::endl;

						
				//internal
				if( uK.FrobeniusNorm() > 0 || uL.FrobeniusNorm() > 0 )
					flag_DIFF = 3;
				else
					flag_DIFF = 2;
			}
			else if( fKL.GetMarker(boundary_face) ) //boundary interface
			{
				//don't forget the area!
				if( split_diffusion )
				{
					if( dK > 0 )
					{
						T = lambdaK / dK;
						gamma = KKn - ((xKL - xK) * lambdaK / dK);
						uK = gamma;
					}
					else if( dK < 0 )
					{
						T = -lambdaK / dK;
						gamma = KKn - ((xKL - xK) * lambdaK / dK);
						uK = 2*KKn - gamma;
					}
					else
					{
						T = 0;
						uK = KKn;
					}
				}
				else //un-splitted diffusion
				{
					uK = KKn;
					T = 0;
				}
				
				if( T < 0 ) std::cout << __FILE__ << ":" << __LINE__ << " T=" << T << std::endl;
						
				real bcconds[3] = {0.0,1.0,0.0}; //pure neumann boundary condition
				if( tag_BC.isValid() && fKL.HaveData(tag_BC) ) //are there boundary conditions on face?
				{
					//retrive boundary conditions
					real_array bc = fKL.RealArray(tag_BC);
					bcconds[0] = bc[0];
					bcconds[1] = bc[1];
					bcconds[2] = bc[2];
				}

				//account for boundary conditions
				R  = A*T*bcconds[2]/(bcconds[0] + bcconds[1]*T);
				uK *=A*  bcconds[0]/(bcconds[0] + bcconds[1]*T);
				T  *=A*  bcconds[0]/(bcconds[0] + bcconds[1]*T);
						
				//on BC
				if( uK.FrobeniusNorm() > 0 )
					flag_DIFF = 1;
				else
					flag_DIFF = 0;
			}
			else std::cout << "No adjacent cell on non-boundary face" << std::endl;
			//record data for diffusion part
			if( tag_K.isValid() )
			{
				fKL.Bulk(tag_DIFF_F) = flag_DIFF;
				fKL.Real(tag_DIFF_RT) = R;
				fKL.Real(tag_DIFF_CT) = T;
				real_array VT = fKL.RealArray(tag_DIFF_VT);
				VT[0] = uK(0,0);
				VT[1] = uK(0,1);
				VT[2] = uK(0,2);
				VT[3] = -uL(0,0);
				VT[4] = -uL(0,1);
				VT[5] = -uL(0,2);
			}
			//Advection part
			uK.Zero();
			uL.Zero();
			C = 0.0;
			cU = InvalidCell();
			R = 0;
			
			if( U > 0.0 ) //flow out of back cell to front cell
			{
				cU = cK;
				C = U*A;
				
				//upstream corrector
				uK = (xK - xKL)*U*A;
						
				//downstream corrector
				if( cL.isValid() ) //internal face
				{
					r0 = r = xK - xL;
					if( tag_K.isValid() ) //heterogeneous media
					{
						KD = KL - KK;
						if( KD.FrobeniusNorm() > 0.0 ) 
						{
							KD /= lambdaK + eps;
							iT = (rMatrix::Unit(3) + nKL.Transpose()*nKL*KD);
							iC = -nKL*dL*KD;
							r = (r*iT - iC)*(lambdaK/(lambdaK+eps)) + r*(eps/(lambdaK+eps));
							//r = (r*iT - iC)*(lambdaK/(lambdaK+eps)) + (r + nKL*(KL-KK)/(U+eps))*(eps/(lambdaK+eps));
						}
					}
					uL = (xKL - xL - r)*U*A;
							
					flag_CONV = 2; //internal
				}
				else if( fKL.GetMarker(boundary_face) ) flag_CONV = 1; //on outlet BC
				else std::cout << "No adjacent cell on non-boundary face" << std::endl;
			}
			else if( U < 0.0 ) //flow out of front cell to back cell
			{
				if( cL.isValid() ) //internal face
				{
					cU = cL;
					C = U*A;
							
					uL = (xKL - xL)*U*A;
							
					r = xL - xK;
					if( tag_K.isValid() ) //heterogeneous media
					{
						KD = KK - KL;
						if( KD.FrobeniusNorm() > 0.0 ) 
						{
							KD /= lambdaL + eps;
							iT = (rMatrix::Unit(3) + nKL.Transpose()*nKL*KD);
							iC = nKL*dK*KD;
							r = (r*iT - iC)*(lambdaL/(lambdaL+eps)) + r*(eps/(lambdaL+eps));
							//r = (r*iT - iC)*(lambdaL/(lambdaL+eps)) + (r + nKL*(KK-KL)/(U-eps))*(eps/(lambdaL+eps));
						}
					}
					uK = (r + xK - xKL)*U*A;
							
					flag_CONV = 2; //internal
				}
				else if( fKL.GetMarker(boundary_face) )//boundary face
				{
					real bcconds[3] = {0.0,1.0,0.0}; //pure neumann boundary condition
					if( tag_BC.isValid() && fKL.HaveData(tag_BC) ) //are there boundary conditions on face?
					{
						//retrive boundary conditions
						real_array bc = fKL.RealArray(tag_BC);
						bcconds[0] = bc[0];
						bcconds[1] = bc[1];
						bcconds[2] = bc[2];
					}
					if( bcconds[0] > 0 && fabs(bcconds[1]) < 1.0e-7 ) //Dirichlet BC
					{
						C = 0;
						R = U*A*bcconds[2]/bcconds[0];
								
						flag_CONV = 0; //on inlet BC
					}
					else if( tag_K.isValid() ) //get expression out of diffusion
					{
						v = - A * (KKn - (xKL - xK) * lambdaK / dK);
						T =   A * lambdaK / dK;
								
						cU = cK;
						C = U*A*T*bcconds[1]/(bcconds[0] + bcconds[1]*T + 1.0e-100);
						R = U*A  *bcconds[2]/(bcconds[0] + bcconds[1]*T + 1.0e-100);
						uK =U*A*v*bcconds[1]/(bcconds[0] + bcconds[1]*T + 1.0e-100);
								
						flag_CONV = 1; //treat as outlet BC
					}
					else
					{
						if( bcconds[0] > 0 && bcconds[1] > 0 ) //Mixed BC
						{
							initialization_failure = true;
							std::cout << "Mixed BC type on inlet boundary face " << fKL->LocalID() << " is inconsistent for pure advection" << std::endl;
						}
						else  if( bcconds[1] > 0 ) // Neumann BC
						{
							initialization_failure = true;
							std::cout << "Neumann BC type on inlet boundary face " << fKL->LocalID() << " is inconsistent for pure advection" << std::endl;
						}
					}
				}
				else std::cout << "No adjacent cell on non-boundary face" << std::endl;
			}
			else flag_CONV = 0; //otherwise velocity is zero, no flow through face
			//record data for advection part
			if( tag_U.isValid() )
			{
				fKL.Bulk(tag_CONV_F) = flag_CONV;
				fKL.Real(tag_CONV_CU) = C;
				fKL.Reference(tag_CONV_EU) = cU.GetHandle();
				fKL.Real(tag_CONV_RU) = R;
				real_array VU = fKL.RealArray(tag_CONV_VU);
				VU[0] = uK(0,0);
				VU[1] = uK(0,1);
				VU[2] = uK(0,2);
				VU[3] = uL(0,0);
				VU[4] = uL(0,1);
				VU[5] = uL(0,2);
			} // tag_U is defined
		} // for q
	} //openmp
	Tag failure = m->CreateTag("FAILURE",DATA_INTEGER,CELL,NONE,1);
	//Compute triplets
	integer negative = 0, total = 0;
#if defined(USE_OMP)
#pragma omp parallel reduction(+:negative,total)
#endif
	{
		Cell cK; //current cell
		Face fKL;
		ElementArray<Face> faces; //faces of the cell
		std::vector<Stencil> compute; //approximations to be computed
		Tag tag_iC, tag_iT; //temporary data used for interpolation tensors
		{
			std::stringstream name;
			name << "INTERPOLATION_TENSOR_" << m->GetLocalProcessorRank();
			tag_iT = m->CreateTag(name.str(),DATA_REAL,CELL,NONE);
			name.clear();
			name << "INTERPOLATION_CORRECTION_" << m->GetLocalProcessorRank();
			tag_iC = m->CreateTag(name.str(),DATA_REAL,CELL,NONE);
		}
		bulk flag;
		real_array v;
#if defined(USE_OMP)
#pragma omp for schedule(dynamic)
#endif
		for(integer q = 0; q < m->CellLastLocalID(); ++q ) if( m->isValidCell(q) )
		{
			cK = m->CellByLocalID(q);
			compute.clear();
			faces = cK.getFaces();
			for(INMOST_DATA_ENUM_TYPE k = 0; k < faces.size(); ++k)
			{
				fKL = faces[k];

				if( !BuildFlux(fKL) ) continue;
				//Diffusion stencils
				if( tag_K.isValid() )
				{
					flag = fKL.Bulk(tag_DIFF_F);
					if( flag == 3 || flag == 1 )
					{
						Stencil s;
						s.type = DIFFUSION;
						if( fKL.FaceOrientedOutside(cK) )
						{
							s.coefs = fKL.RealArray(tag_DIFF_CB);
							s.elems = fKL.ReferenceArray(tag_DIFF_EB);
							s.rhs = &fKL.Real(tag_DIFF_RB);
							v = fKL.RealArray(tag_DIFF_VT);
							s.v[0] = v[0];
							s.v[1] = v[1];
							s.v[2] = v[2];
							s.sign = 1;
						}
						else
						{
							s.coefs = fKL.RealArray(tag_DIFF_CF);
							s.elems = fKL.ReferenceArray(tag_DIFF_EF);
							s.rhs = &fKL.Real(tag_DIFF_RF);
							v = fKL.RealArray(tag_DIFF_VT);
							s.v[0] = v[3];
							s.v[1] = v[4];
							s.v[2] = v[5];
							s.sign = -1;
						}
						if( s.v[0]*s.v[0] + s.v[1]*s.v[1] + s.v[2]*s.v[2] > 0.0 )
							compute.push_back(s);
					}
				}
				//Advection stencils
				if( tag_U.isValid() )
				{
					flag = fKL.Bulk(tag_CONV_F);
					if( flag >= 2 || flag == 1 )
					{
						Stencil s;
						s.type = CONVECTION;
						if( fKL.FaceOrientedOutside(cK) )
						{
							s.coefs = fKL.RealArray(tag_CONV_CB);
							s.elems = fKL.ReferenceArray(tag_CONV_EB);
							s.rhs = &fKL.Real(tag_CONV_RB);
							v = fKL.RealArray(tag_CONV_VU);
							s.v[0] = v[0];
							s.v[1] = v[1];
							s.v[2] = v[2];
							s.sign = -1;
						}
						else
						{
							s.coefs = fKL.RealArray(tag_CONV_CF);
							s.elems = fKL.ReferenceArray(tag_CONV_EF);
							s.rhs = &fKL.Real(tag_CONV_RF);
							v = fKL.RealArray(tag_CONV_VU);
							s.v[0] = v[3];
							s.v[1] = v[4];
							s.v[2] = v[5];
							s.sign = 1;
						}
						if( s.v[0]*s.v[0] + s.v[1]*s.v[1] + s.v[2]*s.v[2] > 0.0 )
							compute.push_back(s);
					}
				}
			}
			//Gathered all the stencils at once
			if( !compute.empty() )
			{
				if( !find_stencils(cK,compute,tag_BC,tag_K,tag_iT,tag_iC,boundary_face,bridge_layers,degenerate_diffusion_regularization,max_layers) )
				{
					initialization_failure = true;
					cK.Integer(failure) = 1;
					std::cout << "Cannot find stencil for cell " << cK->LocalID() << std::endl;
					for(size_t it = 0; it < compute.size(); ++it)
					{
						if( !compute[it].computed )
						{
							std::cout << (compute[it].type == DIFFUSION ? "Diffusion" : "Advection");
							std::cout << " vec (" << compute[it].v[0] << "," << compute[it].v[1] << "," << compute[it].v[2] << ")";
							std::cout << std::endl;
						}
					}
					std::cout << "Total: " << compute.size() << std::endl;
					find_stencils(cK,compute,tag_BC,tag_K,tag_iT,tag_iC,boundary_face,bridge_layers,degenerate_diffusion_regularization,max_layers);
				}
				else for(size_t it = 0; it < compute.size(); ++it)
					if( !compute[it].nonnegative ) negative++;
				total += (integer)compute.size();
			}
		} //end of loop over cells
		//release interpolation tags
		m->DeleteTag(tag_iC);
		m->DeleteTag(tag_iT);
	} //openmp

	negative = m->Integrate(negative);
	total = m->Integrate(total);
	if( m->GetProcessorRank() == 0 )
			std::cout << " negative flux approximation: " << negative << "/" << total << std::endl;

	{ //Synchronize initialization_failure
		integer sync = (initialization_failure ? 1 : 0);
		sync = m->Integrate(sync);
		if( sync ) initialization_failure = true;
	}
	if( !initialization_failure )
		m->DeleteTag(failure);
}

bool ConvectionDiffusion::Failure() const 
{
	return initialization_failure;
}


ConvectionDiffusion::~ConvectionDiffusion()
{
	if( tag_U.isValid() )
	{
		m->DeleteTag(tag_CONV_CB);
		m->DeleteTag(tag_CONV_EB);
		m->DeleteTag(tag_CONV_RB);
		m->DeleteTag(tag_CONV_CF);
		m->DeleteTag(tag_CONV_EF);
		m->DeleteTag(tag_CONV_RF);
		m->DeleteTag(tag_CONV_CU);
		m->DeleteTag(tag_CONV_EU);
		m->DeleteTag(tag_CONV_RU);
		m->DeleteTag(tag_CONV_VU);
		m->DeleteTag(tag_CONV_F);
	}
	if( tag_K.isValid() )
	{
		m->DeleteTag(tag_DIFF_CB);
		m->DeleteTag(tag_DIFF_EB);
		m->DeleteTag(tag_DIFF_RB);
		m->DeleteTag(tag_DIFF_CF);
		m->DeleteTag(tag_DIFF_EF);
		m->DeleteTag(tag_DIFF_RF);
		m->DeleteTag(tag_DIFF_CT);
		m->DeleteTag(tag_DIFF_RT);
		m->DeleteTag(tag_DIFF_VT);
		m->DeleteTag(tag_DIFF_F);
	}
	
	if( build_flux )
		m->ReleaseMarker(build_flux,FACE);
}

void ConvectionDiffusion::DiffusionFluxes(abstract_variable & param, Face fKL, variable & fluxT, variable & corrK, variable & corrL, variable & corrK_cK, variable & corrL_cL) const
{
	if( !BuildFlux(fKL) ) return;

	get_variable<variable> Func(param);
	//diffusion part
	if( tag_K.isValid() )
	{
		Cell cK = fKL->BackCell();
		Cell cL = fKL->FrontCell();
		//diffusion data
		real       coefT  = fKL.Real           (tag_DIFF_CT);
		real       rhsT   = fKL.Real           (tag_DIFF_RT);
		//corrector in back cell
		real_array coefsK = fKL->RealArray     (tag_DIFF_CB);
		ref_array  elemsK = fKL->ReferenceArray(tag_DIFF_EB);
		real       rhsK   = fKL->Real          (tag_DIFF_RB);
		//corrector in front cell
		real_array coefsL = fKL->RealArray     (tag_DIFF_CF);
		ref_array  elemsL = fKL->ReferenceArray(tag_DIFF_EF);
		real       rhsL   = fKL->Real          (tag_DIFF_RF);
		//flag
		bulk       flag   = fKL->Bulk          (tag_DIFF_F);

		//two-point part
		fluxT -= rhsT;
		if( flag > 1 ) fluxT  -= coefT*(Func(cL) - Func(cK));
		else fluxT  += coefT*Func(cK);

		if( flag > 1 ) //full nonlinear part
		{
			if( perform_correction_diffusion )
			{
				corrK -= rhsK;
				corrL -= rhsL;
				for(enumerator k = 0; k < elemsK.size(); ++k)
					if( elemsK[k].isValid() ) corrK -= Func(elemsK[k])*coefsK[k];
				for(enumerator k = 0; k < elemsL.size(); ++k)
					if( elemsL[k].isValid() ) corrL -= Func(elemsL[k])*coefsL[k];
				//if( scheme_type == NTPFA || scheme_type == NTPFA_PICARD )
				{
					if( elemsK.back().isValid() ) corrK_cK -= Func(elemsK.back())*coefsK.back();
					if( elemsL.back().isValid() ) corrL_cL -= Func(elemsL.back())*coefsL.back();
				}
			}
		}
		else //one-sided corrector
		{
			if( flag == 1 && perform_correction_diffusion ) // only fluxK is present
			{
				corrK -= rhsK;
				for(enumerator k = 0; k < elemsK.size(); ++k)
					if( elemsK[k].isValid() ) corrK -= Func(elemsK[k])*coefsK[k];
			} //no corrections otherwise
		} //flag
	} //diffusion part
}


bool ConvectionDiffusion::BuildFlux(Face f) const
{
	if( build_flux && !f.GetMarker(build_flux) ) return false;
	return true;
}

void ConvectionDiffusion::AdvectionFluxes(abstract_variable & param, Face fKL, variable & fluxT, variable & corrK, variable & corrL, variable & corrK_cK, variable & corrL_cL) const
{
	if( !BuildFlux(fKL) ) return;

	get_variable<variable> Func(param);
	//advection part
	if( tag_U.isValid() )
	{
		//upwind data
		real       coefU  = fKL.Real            (tag_CONV_CU);
		real       rhsU   = fKL.Real            (tag_CONV_RU);
		Cell       cU     = Cell(m,fKL.Reference(tag_CONV_EU));
		//corrector in back cell
		real_array coefsK = fKL->RealArray      (tag_CONV_CB);
		ref_array  elemsK = fKL->ReferenceArray (tag_CONV_EB);
		real       rhsK   = fKL->Real           (tag_CONV_RB);
		//corrector in front cell
		real_array coefsL = fKL->RealArray      (tag_CONV_CF);
		ref_array  elemsL = fKL->ReferenceArray (tag_CONV_EF);
		real       rhsL   = fKL->Real           (tag_CONV_RF);
		//flag
		bulk       flag   = fKL->Bulk           (tag_CONV_F);

		fluxT += rhsU;
		if( cU.isValid() ) fluxT += Func(cU)*coefU;

		if( flag == 2 ) //both fluxes present
		{
			if( perform_correction_convection )
			{
				corrK += rhsK;
				corrL += rhsL;
				for(enumerator k = 0; k < elemsK.size(); ++k)
					if( elemsK[k].isValid() ) corrK += Func(elemsK[k])*coefsK[k];
				for(enumerator k = 0; k < elemsL.size(); ++k)
					if( elemsL[k].isValid() ) corrL += Func(elemsL[k])*coefsL[k];

				//if( scheme_type == NTPFA || scheme_type == NTPFA_PICARD )
				{
					if( elemsK.back().isValid() ) corrK_cK += Func(elemsK.back())*coefsK.back();
					if( elemsL.back().isValid() ) corrL_cL += Func(elemsL.back())*coefsL.back();
				}
			}
		}
		else
		{
			if( flag == 1 && perform_correction_convection ) // only corrK is present
			{
				corrK += rhsK;
				for(enumerator k = 0; k < elemsK.size(); ++k)
					if( elemsK[k].isValid() ) corrK += Func(elemsK[k])*coefsK[k];
			} //no corrections otherwise
		}
		
	} //advection part
}

void ConvectionDiffusion::Averaging(SchemeType scheme_type, INMOST_DATA_REAL_TYPE regularization, const variable & fluxT, const variable & corrK, const variable & corrL, const variable & corrK_cK, const variable & corrL_cL, variable & fluxK, variable & fluxL, bool boundary) const
{
	if( boundary )
	{
		//compute flux with correction
		fluxK = fluxT + corrK;
	}
	else
	{
		//perform nonlinear weighting
		if( scheme_type == NTPFA || scheme_type == NTPFA_PICARD )
			LimitedAverageNonlinearTwoPoint(corrK,corrL,corrK_cK,corrL_cL,fluxK,fluxL,regularization,scheme_type);
		else
			LimitedAverageNonlinearMultiPoint(corrK,corrL,fluxK,fluxL,regularization,scheme_type);
		//add upwind/two-point part
		fluxK += fluxT;
		fluxL += fluxT;
	}
}
