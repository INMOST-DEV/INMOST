#include "inmost_model.h"
#if defined(USE_AUTODIFF) && defined(USE_MESH) && defined(USE_SOLVER)

namespace INMOST
{
	
	void Model::AddEntry(std::string name, AbstractEntry & entry)
	{
		assert( !isInitialized() ); // do not add new data after initialization
		assert( GetEntry(name) == NULL ); // do not overwrite entry with the same name
		Entries.push_back(std::make_pair(name,&entry));
	}
	
	void Model::AddFirstEntry(std::string name, AbstractEntry & entry)
	{
		assert( !isInitialized() ); // do not add new data after initialization
		assert( GetEntry(name) == NULL ); // do not overwrite entry with the same name
		Entries.insert(Entries.begin(),std::make_pair(name,&entry));
	}
	
	void Model::AddAfterEntry(std::string name, AbstractEntry & entry, std::string after)
	{
		assert( !isInitialized() ); // do not add new data after initialization
		assert( GetEntry(name) == NULL ); // do not overwrite entry with the same name
		std::vector< std::pair< std::string, AbstractEntry *> >::iterator find = Entries.end();
		for(std::vector< std::pair< std::string, AbstractEntry *> >::iterator it = Entries.begin(); it != Entries.end(); ++it)
			if( it->first == after )
			{
				find = ++it;
				break;
			}
		Entries.insert(find,std::make_pair(name,&entry));
	}
	
	void Model::AddMesh(std::string name, Mesh & m)
	{
		assert( !isInitialized() ); // do not add new data after initialization
		assert( GetMesh(name) == NULL ); // do not overwrite mesh with the same name
		Meshes.push_back(std::make_pair(name,&m));
	}
	
	void Model::AddSubModel(std::string name, AbstractSubModel & submodel)
	{
		assert( !isInitialized() ); // do not add new data after initialization
		assert( GetSubModel(name) == NULL ); // do not overwrite submodel with the same name
		SubModels.push_back(std::make_pair(name,&submodel));
	}

	void Model::AddCouplingTerm(std::string submodel, std::string name, AbstractCouplingTerm& term)
	{
		assert(!isInitialized()); // do not add new data after initialization
		assert(GetCouplingTerm(submodel, name) == NULL); //do not overwrite existing coupling term
		AbstractSubModel* subm = GetSubModel(submodel);
		if (subm == NULL && submodel.size() != 0)
			std::cout << "Warning: submodel " << submodel << " corresponding to coupling term " << name << " does not exist!" << std::endl; //TODO: warnings in parallel
		CouplingTerms[subm].push_back(std::make_pair(name, &term));
	}

	void Model::AddScalarFunction(std::string name, AbstractScalarFunction& func)
	{
		assert(!isInitialized()); // do not add new data after initialization
		assert(GetScalarFunction(name) == NULL); //do not overwrite existing coupling term
		ScalarFunctions.push_back(std::make_pair(name, &func));
	}

	void Model::AddMatrixFunction(std::string name, AbstractMatrixFunction& func)
	{
		assert(!isInitialized()); // do not add new data after initialization
		assert(GetMatrixFunction(name) == NULL); //do not overwrite existing coupling term
		MatrixFunctions.push_back(std::make_pair(name, &func));
	}
		
	void Model::AddOperator(std::string name, AbstractOperator & op)
	{
		assert( !isInitialized() ); // do not add new data after initialization
		assert( GetOperator(name) == NULL ); // do not overwrite submodel with the same name
		Operators.push_back(std::make_pair(name,&op));
	}
	
	Mesh * Model::GetMesh(std::string name)
	{
		for(std::vector< std::pair<std::string,Mesh *> >::iterator it = Meshes.begin();
			it != Meshes.end(); ++it)
			if( it->first == name )
				return it->second;
		return NULL;
	}
	
	const Mesh * Model::GetMesh(std::string name) const
	{
		for(std::vector< std::pair<std::string,Mesh *> >::const_iterator it = Meshes.begin();
			it != Meshes.end(); ++it)
			if( it->first == name )
				return it->second;
		return NULL;
	}
	
	std::vector<std::string> Model::GetMeshesNames() const
	{
		std::vector<std::string> ret;
		for(std::vector< std::pair<std::string,Mesh *> >::const_iterator it = Meshes.begin();
			it != Meshes.end(); ++it)
			ret.push_back(it->first);
		return ret;
	}
	
	AbstractEntry * Model::GetEntry(std::string name)
	{
		for(std::vector< std::pair<std::string,AbstractEntry *> >::iterator it = Entries.begin();
			it != Entries.end(); ++it)
			if( it->first == name )
				return it->second;
		return NULL;
	}
	
	const AbstractEntry * Model::GetEntry(std::string name) const
	{
		for(std::vector< std::pair<std::string,AbstractEntry *> >::const_iterator it = Entries.begin();
			it != Entries.end(); ++it)
			if( it->first == name )
				return it->second;
		return NULL;
	}
	
	std::vector<std::string> Model::GetEntriesNames() const
	{
		std::vector<std::string> ret;
		for(std::vector< std::pair<std::string,AbstractEntry *> >::const_iterator it = Entries.begin();
			it != Entries.end(); ++it)
			ret.push_back(it->first);
		return ret;
	}
	
	AbstractSubModel * Model::GetSubModel(std::string name)
	{
		for(std::vector< std::pair<std::string,AbstractSubModel *> >::iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			if( it->first == name )
				return it->second;
		return NULL;
	}
	
	const AbstractSubModel * Model::GetSubModel(std::string name) const
	{
		for(std::vector< std::pair<std::string,AbstractSubModel *> >::const_iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			if( it->first == name )
				return it->second;
		return NULL;
	}

	std::string Model::GetSubModelName(const AbstractSubModel* ref) const
	{
		if (ref == NULL) return "";
		for (std::vector< std::pair<std::string, AbstractSubModel*> >::const_iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			if (it->second == ref)
				return it->first;
		std::cout << "Warning: the sub-model with address " << (void*)ref << " was not associated with the model!" << std::endl; //TODO: warnings in parallel
		return "";
	}
	
	std::vector<std::string> Model::GetSubModelsNames() const
	{
		std::vector<std::string> ret;
		for(std::vector< std::pair<std::string,AbstractSubModel *> >::const_iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			ret.push_back(it->first);
		return ret;
	}

	AbstractCouplingTerm* Model::GetCouplingTerm(std::string submodel, std::string name)
	{
		AbstractSubModel* subm = GetSubModel(submodel);
		std::map< AbstractSubModel*, std::vector< std::pair< std::string, AbstractCouplingTerm*> > >::iterator terms = CouplingTerms.find(subm);
		for (std::vector< std::pair<std::string, AbstractCouplingTerm*> >::iterator it = terms->second.begin();
			it != terms->second.end(); ++it)
			if (it->first == name)
				return it->second;
		return NULL;
	}

	const AbstractCouplingTerm* Model::GetCouplingTerm(std::string submodel, std::string name) const
	{
		const AbstractSubModel* subm = GetSubModel(submodel);
		std::map< AbstractSubModel*, std::vector< std::pair< std::string, AbstractCouplingTerm*> > >::const_iterator terms = CouplingTerms.find(const_cast<AbstractSubModel*>(subm));
		for (std::vector< std::pair<std::string, AbstractCouplingTerm*> >::const_iterator it = terms->second.begin();
			it != terms->second.end(); ++it)
			if (it->first == name)
				return it->second;
		return NULL;
	}

	std::vector< std::pair< std::string, std::vector<std::string> > > Model::GetCouplingTermNames() const
	{
		std::vector< std::pair< std::string, std::vector<std::string> > > ret;
		std::map< AbstractSubModel*, std::vector< std::pair< std::string, AbstractCouplingTerm*> > >::const_iterator terms;
		for (terms = CouplingTerms.begin(); terms != CouplingTerms.end(); ++terms)
		{
			ret.push_back(std::make_pair(GetSubModelName(terms->first), std::vector<std::string>()));
			for (std::vector< std::pair<std::string, AbstractCouplingTerm*> >::const_iterator it = terms->second.begin();
				it != terms->second.end(); ++it)
				ret.back().second.push_back(it->first);
		}	
		return ret;
	}

	AbstractScalarFunction* Model::GetScalarFunction(std::string name)
	{
		for (std::vector< std::pair<std::string, AbstractScalarFunction*> >::iterator it = ScalarFunctions.begin();
			it != ScalarFunctions.end(); ++it)
			if (it->first == name)
				return it->second;
		return NULL;
	}

	const AbstractScalarFunction* Model::GetScalarFunction(std::string name) const
	{
		for (std::vector< std::pair<std::string, AbstractScalarFunction*> >::const_iterator it = ScalarFunctions.begin();
			it != ScalarFunctions.end(); ++it)
			if (it->first == name)
				return it->second;
		return NULL;
	}



	std::vector<std::string> Model::GetScalarFunctionNames() const
	{
		std::vector<std::string> ret;
		for (std::vector< std::pair<std::string, AbstractMatrixFunction*> >::const_iterator it = MatrixFunctions.begin();
			it != MatrixFunctions.end(); ++it)
			ret.push_back(it->first);
		return ret;
	}

	AbstractMatrixFunction* Model::GetMatrixFunction(std::string name)
	{
		for (std::vector< std::pair<std::string, AbstractMatrixFunction*> >::iterator it = MatrixFunctions.begin();
			it != MatrixFunctions.end(); ++it)
			if (it->first == name)
				return it->second;
		return NULL;
	}

	const AbstractMatrixFunction* Model::GetMatrixFunction(std::string name) const
	{
		for (std::vector< std::pair<std::string, AbstractMatrixFunction*> >::const_iterator it = MatrixFunctions.begin();
			it != MatrixFunctions.end(); ++it)
			if (it->first == name)
				return it->second;
		return NULL;
	}



	std::vector<std::string> Model::GetMatrixFunctionNames() const
	{
		std::vector<std::string> ret;
		for (std::vector< std::pair<std::string, AbstractMatrixFunction*> >::const_iterator it = MatrixFunctions.begin();
			it != MatrixFunctions.end(); ++it)
			ret.push_back(it->first);
		return ret;
	}


	
	AbstractOperator * Model::GetOperator(std::string name)
	{
		for(std::vector< std::pair<std::string,AbstractOperator *> >::iterator it = Operators.begin();
			it != Operators.end(); ++it)
			if( it->first == name )
				return it->second;
		return NULL;
	}
	
	const AbstractOperator * Model::GetOperator(std::string name) const
	{
		for(std::vector< std::pair<std::string,AbstractOperator *> >::const_iterator it = Operators.begin();
			it != Operators.end(); ++it)
			if( it->first == name )
				return it->second;
		return NULL;
	}
	
	std::vector<std::string> Model::GetOperatorsNames() const
	{
		std::vector<std::string> ret;
		for(std::vector< std::pair<std::string,AbstractOperator *> >::const_iterator it = Operators.begin();
			it != Operators.end(); ++it)
			ret.push_back(it->first);
		return ret;
	}
	
	bool Model::PrepareEntries()
	{
		bool success = true;
		//first initialize submodels
		for(std::vector< std::pair<std::string, AbstractSubModel *> >::iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			success &= it->second->PrepareEntries(*this);
		return success;
	}

	bool Model::Initialize()
	{
		bool success = true;
		//first initialize submodels
		for(std::vector< std::pair<std::string, AbstractSubModel *> >::iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			success &= it->second->Initialize(*this);
		//second initialize operators
		for(std::vector< std::pair<std::string, AbstractOperator *> >::iterator it = Operators.begin();
			it != Operators.end(); ++it)
			success &= it->second->Initialize(*this);
		//third register all the entries
		for(std::vector< std::pair<std::string, AbstractEntry *> >::iterator it = Entries.begin();
			it != Entries.end(); ++it)
		{
			it->second->reg_index = aut.RegisterEntry(*it->second);
			it->second->SetOffsetTag(aut.GetEntry(it->second->reg_index).GetOffsetTag());
		}
		//initialize coupling terms
		for (std::map< AbstractSubModel*, std::vector< std::pair<std::string, AbstractCouplingTerm*> > >::iterator it = CouplingTerms.begin();
			it != CouplingTerms.end(); ++it)
		{
			for (std::vector< std::pair<std::string, AbstractCouplingTerm*> >::iterator jt = it->second.begin(); jt != it->second.end(); ++jt)
				success &= jt->second->Initialize(*this);
		}
		//initialize scalar functions
		for (std::vector< std::pair<std::string, AbstractScalarFunction*> >::iterator it = ScalarFunctions.begin();
			it != ScalarFunctions.end(); ++it)
			success &= it->second->Initialize(*this);
		//initialize matrix functions
		for (std::vector< std::pair<std::string, AbstractMatrixFunction*> >::iterator it = MatrixFunctions.begin();
			it != MatrixFunctions.end(); ++it)
			success &= it->second->Initialize(*this);
		initialized = success;
		return success;
	}
	
	bool Model::PrepareIterations()
	{
		bool success = true;
		for(std::vector< std::pair<std::string, AbstractSubModel *> >::iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			success &= it->second->PrepareIterations();
		//for operator
		for (std::vector< std::pair<std::string, AbstractOperator*> >::iterator it = Operators.begin();
			it != Operators.end(); ++it)
			success &= it->second->PrepareIterations();
		//for coupling terms
		for (std::map< AbstractSubModel*, std::vector< std::pair<std::string, AbstractCouplingTerm*> > >::iterator it = CouplingTerms.begin();
			it != CouplingTerms.end(); ++it)
		{
			for (std::vector< std::pair<std::string, AbstractCouplingTerm*> >::iterator jt = it->second.begin(); jt != it->second.end(); ++jt)
				success &= jt->second->PrepareIterations();
		}
		//for scalar functions
		for (std::vector< std::pair<std::string, AbstractScalarFunction*> >::iterator it = ScalarFunctions.begin();
			it != ScalarFunctions.end(); ++it)
			success &= it->second->PrepareIterations();
		//for matrix functions
		for (std::vector< std::pair<std::string, AbstractMatrixFunction*> >::iterator it = MatrixFunctions.begin();
			it != MatrixFunctions.end(); ++it)
			success &= it->second->PrepareIterations();
		return success;
	}
	
	bool Model::FillResidual(Residual & R) const
	{
		bool success = true;
		R.Clear();
		Automatizator::MakeCurrent(&aut);
		double total_time = Timer(), model_time = 0;
		for(std::vector< std::pair<std::string, AbstractSubModel *> >::const_iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
		{
			model_time = Timer();
			success &= it->second->FillResidual(R);
			model_time = Timer() - model_time;
			//~ std::cout << it->first << ":" << model_time << " ";
		}
		for (std::map< AbstractSubModel*, std::vector< std::pair<std::string, AbstractCouplingTerm*> > >::const_iterator it = CouplingTerms.begin();
			it != CouplingTerms.end(); ++it)
		{
			for (std::vector< std::pair<std::string, AbstractCouplingTerm*> >::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt)
			{
				model_time = Timer();
				success &= jt->second->FillResidual(R);
				model_time = Timer() - model_time;
				//~ std::cout << it->first << ":" << jt->first << ":" << model_time << " ";
			}
		}
		total_time = Timer() - total_time;
		//~ std::cout << "total:" << total_time << std::endl;
		Automatizator::RemoveCurrent();
		return success;
	}
	
	bool Model::UpdateSolution(const Sparse::Vector & sol, double alpha)
	{
		bool success = true;
		for(std::vector< std::pair<std::string, AbstractSubModel *> >::const_iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			success &= it->second->UpdateSolution(sol,alpha);
		return success;
	}
	
	bool Model::UpdateTimeStep()
	{
		bool success = true;
		for(std::vector< std::pair<std::string, AbstractSubModel *> >::iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			success &= it->second->UpdateTimeStep();
		for (std::map< AbstractSubModel*, std::vector< std::pair<std::string, AbstractCouplingTerm*> > >::iterator it = CouplingTerms.begin();
			it != CouplingTerms.end(); ++it)
		{
			for (std::vector< std::pair<std::string, AbstractCouplingTerm*> >::iterator jt = it->second.begin(); jt != it->second.end(); ++jt)
				success &= jt->second->UpdateTimeStep();
		}
		for (std::vector< std::pair<std::string, AbstractScalarFunction*> >::iterator it = ScalarFunctions.begin();
			it != ScalarFunctions.end(); ++it)
			success &= it->second->UpdateTimeStep();
		for (std::vector< std::pair<std::string, AbstractMatrixFunction*> >::iterator it = MatrixFunctions.begin();
			it != MatrixFunctions.end(); ++it)
			success &= it->second->UpdateTimeStep();
		return success;
	}
	
	bool Model::SetTimeStep(double dt)
	{
		bool success = true;
		for(std::vector< std::pair<std::string, AbstractSubModel *> >::const_iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			success &= it->second->SetTimeStep(dt);

		for (std::map< AbstractSubModel*, std::vector< std::pair<std::string, AbstractCouplingTerm*> > >::iterator it = CouplingTerms.begin();
			it != CouplingTerms.end(); ++it)
		{
			for (std::vector< std::pair<std::string, AbstractCouplingTerm*> >::iterator jt = it->second.begin(); jt != it->second.end(); ++jt)
				success &= jt->second->SetTimeStep(dt);
		}
		for (std::vector< std::pair<std::string, AbstractScalarFunction*> >::const_iterator it = ScalarFunctions.begin();
			it != ScalarFunctions.end(); ++it)
			success &= it->second->SetTimeStep(dt);
		for (std::vector< std::pair<std::string, AbstractMatrixFunction*> >::const_iterator it = MatrixFunctions.begin();
			it != MatrixFunctions.end(); ++it)
			success &= it->second->SetTimeStep(dt);
		return success;
		
	}
	
	bool Model::RestoreTimeStep()
	{
		bool success = true;
		for(std::vector< std::pair<std::string, AbstractSubModel *> >::iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			success &= it->second->RestoreTimeStep();
		for (std::map< AbstractSubModel*, std::vector< std::pair<std::string, AbstractCouplingTerm*> > >::iterator it = CouplingTerms.begin();
			it != CouplingTerms.end(); ++it)
		{
			for (std::vector< std::pair<std::string, AbstractCouplingTerm*> >::iterator jt = it->second.begin(); jt != it->second.end(); ++jt)
				success &= jt->second->RestoreTimeStep();
		}
		for (std::vector< std::pair<std::string, AbstractScalarFunction*> >::const_iterator it = ScalarFunctions.begin();
			it != ScalarFunctions.end(); ++it)
			success &= it->second->RestoreTimeStep();
		for (std::vector< std::pair<std::string, AbstractMatrixFunction*> >::const_iterator it = MatrixFunctions.begin();
			it != MatrixFunctions.end(); ++it)
			success &= it->second->RestoreTimeStep();
		return success;
	}
	
	double Model::UpdateMultiplier(const Sparse::Vector & sol) const
	{
		double alpha = 1;
		for(std::vector< std::pair<std::string, AbstractSubModel *> >::const_iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			alpha = std::min(alpha,it->second->UpdateMultiplier(sol));

		for (std::map< AbstractSubModel*, std::vector< std::pair<std::string, AbstractCouplingTerm*> > >::const_iterator it = CouplingTerms.begin();
			it != CouplingTerms.end(); ++it)
		{
			for (std::vector< std::pair<std::string, AbstractCouplingTerm*> >::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt)
				alpha = std::min(alpha,jt->second->UpdateMultiplier(sol));
		}
		for (std::vector< std::pair<std::string, AbstractScalarFunction*> >::const_iterator it = ScalarFunctions.begin();
			it != ScalarFunctions.end(); ++it)
			alpha = std::min(alpha, it->second->UpdateMultiplier(sol));
		for (std::vector< std::pair<std::string, AbstractMatrixFunction*> >::const_iterator it = MatrixFunctions.begin();
			it != MatrixFunctions.end(); ++it)
			alpha = std::min(alpha, it->second->UpdateMultiplier(sol));

#if defined(USE_MPI)
		double tmp = alpha;
		MPI_Allreduce(&tmp,&alpha,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#endif
		return alpha;
	}
	
	double Model::AdjustTimeStep(double dt) const
	{
		double ret = dt;
		for(std::vector< std::pair<std::string, AbstractSubModel *> >::const_iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			ret = std::min(ret,it->second->AdjustTimeStep(dt));

		for (std::map< AbstractSubModel*, std::vector< std::pair<std::string, AbstractCouplingTerm*> > >::const_iterator it = CouplingTerms.begin();
			it != CouplingTerms.end(); ++it)
		{
			for (std::vector< std::pair<std::string, AbstractCouplingTerm*> >::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt)
				ret = std::min(ret, jt->second->AdjustTimeStep(dt));
		}
		for (std::vector< std::pair<std::string, AbstractScalarFunction*> >::const_iterator it = ScalarFunctions.begin();
			it != ScalarFunctions.end(); ++it)
			ret = std::min(ret, it->second->AdjustTimeStep(dt));
		for (std::vector< std::pair<std::string, AbstractMatrixFunction*> >::const_iterator it = MatrixFunctions.begin();
			it != MatrixFunctions.end(); ++it)
			ret = std::min(ret, it->second->AdjustTimeStep(dt));
#if defined(USE_MPI)
		double tmp = ret;
		MPI_Allreduce(&tmp,&ret,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#endif
		return ret;
	}
	
	
	void Model::ReportErrors(const Residual & R) const
	{
		for(std::vector< std::pair< std::string, AbstractEntry *> >::const_iterator it = Entries.begin(); it != Entries.end(); ++it)
		{
			ElementType etypes = it->second->GetElementType();
			Mesh * m = it->second->GetMeshLink();
			//if( m->GetProcessorRank() == 0 )
			//	std::cout << "Errors for entry " << it->first << std::endl;
			for(ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype) ) if ( etype & etypes )
			{
				if( !m->HaveGlobalID(etype) ) m->AssignGlobalID(etype); //to report element number
				//account for variable size!
				std::vector<INMOST_DATA_REAL_TYPE> err_comp(it->second->Size(),0.0);
				//std::vector<Element> err_elem(it->second->Size());
				std::vector<INMOST_DATA_REAL_TYPE> err_int(it->second->Size(),0.0);
				INMOST_DATA_REAL_TYPE max_err = 0, int_err = 0;
				Element e;
				for(Mesh::iteratorElement jt = m->BeginElement(etype); jt != m->EndElement(); ++jt) if( jt->GetStatus() != Element::Ghost )
				{
					rMatrix err = R.Value(it->second->Index(jt->self()));
					double block_err = err.FrobeniusNorm();
					if( block_err > max_err )
					{
						e = jt->self();
						max_err = block_err;
					}
					int_err += block_err*block_err;
					INMOST_DATA_ENUM_TYPE N = err.Rows()*err.Cols();
					if( N > it->second->Size() ) continue; //No account for variable size
					//~ std::cout << jt->LocalID() << " block err " << block_err << " comp err ";
					for(INMOST_DATA_ENUM_TYPE k = 0; k < N; ++k)
					{
						//~ std::cout << err.data()[k] << " ";
						err_int[k] += err.data()[k]*err.data()[k];
						if( fabs(err.data()[k]) > err_comp[k] )
						{
							err_comp[k] = fabs(err.data()[k]);
							//err_elem[k] = jt->self();
						}
					}
					//~ std::cout << std::endl;
				}
				m->Integrate(&err_int[0],it->second->Size());
				m->AggregateMax(&err_comp[0],it->second->Size());
				int_err = m->Integrate(int_err);
				max_err = m->AggregateMax(max_err); //we need MPI_Allreduce to zero proc to know index of element with largest error
				if( m->GetProcessorRank() == 0 )
				{
					std::cout << std::setw(30) << it->first << " element type " << ElementTypeName(etype) << " error integral " << std::setw(12) << sqrt(int_err) << " maximal " << std::setw(12) << max_err << std::endl;// " on element " << e.GlobalID() << std::endl;
					for(int k = 0; k < (int)err_int.size(); ++k)
					{
						std::cout << "\t\t" << std::setw(2) << k << " error integral " << std::setw(12) << sqrt(err_int[k]) << " maximal " << std::setw(12) << err_comp[k] << " tag " << it->second->GetValueTag(k).GetTagName() << std::endl;
					}
				}
			}
		}
	}


	/*
	
	void Model::PrepareAdaptation(Mesh &m)
	{
		//first initialize submodels
		for (std::vector< std::pair<std::string, AbstractSubModel*> >::iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			it->second->PrepareAdaptation(m);
	}
	
	
	void Model::Adaptation(Mesh & m) const
	{
		//~ array<HandleType> old_cells;
		//~ for(int k = 0; k < m.CellLastLocalID(); ++k) if( m.isValidCell(k) )
		//~ {
			//~ Cell c = m.CellByLocalID(k);
			//~ if( c.Hidden() ) old_cells.push_back(c.GetHandle());
		//~ }
		//std::cout << "old cells: " << old_cells.size() << std::endl;
		//~ SearchKDTree tree(&m,old_cells.data(),old_cells.size());
		for(std::vector< std::pair<std::string, AbstractSubModel *> >::const_iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			it->second->Adaptation(m);
	}

	void Model::NewNode(Cell& c, Node &n, Storage::reference_array cell_hanging_nodes)
	{
		for (std::vector< std::pair<std::string, AbstractSubModel*> >::const_iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			it->second->NewNode(c, n, cell_hanging_nodes);
	}
	void Model::NewNode(Face& f, Node& n, Storage::reference_array face_hanging_nodes)
	{
		for (std::vector< std::pair<std::string, AbstractSubModel*> >::const_iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			it->second->NewNode(f, n, face_hanging_nodes);
	}
	void Model::NewNode(Edge& e, Node& n)
	{
		for (std::vector< std::pair<std::string, AbstractSubModel*> >::const_iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			it->second->NewNode(e, n);
	}

	void Model::NewEdge(Cell& c, Edge& e)
	{
		for (std::vector< std::pair<std::string, AbstractSubModel*> >::const_iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			it->second->NewEdge(c, e);
	}
	void Model::NewEdge(Face& f, Edge& e)
	{
		for (std::vector< std::pair<std::string, AbstractSubModel*> >::const_iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			it->second->NewEdge(f, e);
	}
	void Model::NewFace(Cell& c, Face& f)
	{
		for (std::vector< std::pair<std::string, AbstractSubModel*> >::const_iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			it->second->NewFace(c, f);
	}
	
	void Model::CellRefinement(Cell & old_cell, ElementArray<Cell> & new_cells, ElementSet & new_cell_set)
	{
		for(std::vector< std::pair<std::string, AbstractSubModel *> >::const_iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			it->second->CellRefinement(old_cell,new_cells,new_cell_set);
	}
	void Model::FaceRefinement(Face & old_face, ElementArray<Face> & new_faces)
	{
		for(std::vector< std::pair<std::string, AbstractSubModel *> >::const_iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			it->second->FaceRefinement(old_face,new_faces);
	}
	void Model::EdgeRefinement(Edge & old_edge, ElementArray<Edge> & new_edges)
	{
		for(std::vector< std::pair<std::string, AbstractSubModel *> >::const_iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			it->second->EdgeRefinement(old_edge,new_edges);
	}

	
	void Model::CellCoarsening(ElementArray<Cell> & old_cells, Cell & new_cell, ElementSet & old_cells_set)
	{
		for(std::vector< std::pair<std::string, AbstractSubModel *> >::const_iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			it->second->CellCoarsening(old_cells,new_cell,old_cells_set);
	}
	void Model::FaceCoarsening(ElementArray<Face> & old_faces, Face & new_face)
	{
		for(std::vector< std::pair<std::string, AbstractSubModel *> >::const_iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			it->second->FaceCoarsening(old_faces,new_face);
	}
	void Model::EdgeCoarsening(ElementArray<Edge> & old_edges, Edge & new_edge)
	{
		for(std::vector< std::pair<std::string, AbstractSubModel *> >::const_iterator it = SubModels.begin();
			it != SubModels.end(); ++it)
			it->second->EdgeCoarsening(old_edges,new_edge);
	}
	*/
}
#endif
