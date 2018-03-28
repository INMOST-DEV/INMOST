#ifndef INMOST_MODEL_INCLUDED
#define INMOST_MODEL_INCLUDED

#include "inmost_common.h"
#include "inmost_autodiff.h"
#include "inmost_residual.h"
#include "inmost_operator.h"

#if defined(USE_AUTODIFF) && defined(USE_MESH) && defined(USE_SOLVER)

namespace INMOST
{
	
	class Model;
	class AbstractOperator;
	
	/// A class to manage a submodel within a model.
	/// Each submodel is responsible to define unknowns
	/// of the model on the mesh.
	class AbstractSubModel
	{
	public:
		/// Let the submodel introduce it's unknowns
		virtual bool PrepareEntries(Model & P) = 0;
		/// Initialize coupling and dependent unknowns.
		virtual bool Initialize(Model & P) = 0;
		/// Initialize data needed for FillResidual.
		/// Called once before nonlinear iterations.
		virtual bool PrepareIterations() {return true;}
		/// Fill part of the residual related to my unknowns.
		virtual bool FillResidual(Residual & R) const = 0;
		/// Update solution.
		virtual bool UpdateSolution(const Sparse::Vector & sol, double alpha) = 0;
		/// Update time step.
		virtual bool UpdateTimeStep() = 0;
		/// Provide time step.
		virtual bool SetTimeStep(double dt) = 0;
		/// Roll back to previous step.
		virtual bool RestoreTimeStep() = 0;
		/// Calculate multiplier for update for this model. Can simply return 1.
                virtual double UpdateMultiplier(const Sparse::Vector & sol) const {(void)sol; return 1;}
		/// Calculate time step for this model. Can simply return dt.
		virtual double AdjustTimeStep(double dt) const {return dt;}
	};
	
	/// A class to organize a model.
	class Model
	{
		Automatizator & aut; ///< Automatizator that is used to manage all unknowns of the model.
		//todo: decide later on how it should be stored
		std::vector< std::pair<std::string, AbstractSubModel *> > SubModels; ///< A set of submodels of the model.
		std::vector< std::pair< std::string, AbstractEntry *> > Entries; ///< A set of entries for blocks of unknowns of the model.
		std::vector< std::pair<std::string, Mesh *> > Meshes; ///< A set of meshes of the model.
		std::vector< std::pair<std::string, AbstractOperator *> > Operators; ///< A set of operators used by submodels
		bool initialized; ///< Indicates whether a model was initialized.
	public:
		Model(Automatizator & aut) : aut(aut), initialized(false) {}
		//todo:
		//Model(const Model & b) aut(b.aut) {}
		//todo:
		//Model & operator =(Model const & b);
		virtual ~Model() {}
		/// Add an entry of block unknowns to a model.
		/// The model stores a link to the entry and may modify it contents.
		/// The intries should be added from Model::Initialize function,
		/// either by model or by any of the submodels.
		void AddEntry(std::string name, AbstractEntry & entry);
		/// Add a mesh to a model.
		/// The model stores a link to the provided mesh, so it should not
		/// be deallocated. The meshes are provided by the user from outside
		/// before Model::Initialize function was called.
		/// The meshes are requested by name by each submodel.
		/// Same mesh can be added with different names for submodels.
		void AddMesh(std::string name, Mesh & m);
		/// Add a submodel to a model.
		/// Submodels are added by the user from outside.
		/// All submodels are initialized on Model::Initialize function.
		void AddSubModel(std::string name, AbstractSubModel & submodel);
		/// Add an operator to a model.
		void AddOperator(std::string name, AbstractOperator & op);
		/// Retrive an automatizator.
		Automatizator & GetAutomatizator() {return aut;}
		/// Retrive an automatizator.
		const Automatizator & GetAutomatizator() const {return aut;}
		/// Retrive a mesh by name.
		Mesh * GetMesh(std::string);
		/// Retrive a mesh by name.
		const Mesh * GetMesh(std::string) const;
		/// Retrive all names of meshes.
		std::vector<std::string> GetMeshesNames() const;
		/// Retrive an entry that describe unknowns of the model by name.
		AbstractEntry * GetEntry(std::string name);
		/// Retrive an entry that describe unknowns of the model by name.
		const AbstractEntry * GetEntry(std::string name) const;
		/// Retrive all names of entries.
		std::vector<std::string> GetEntriesNames() const;
		/// Retrive a submodel of the model by name.
		AbstractSubModel * GetSubModel(std::string name);
		/// Retrive a submodel of the model by name.
		const AbstractSubModel * GetSubModel(std::string name) const;
		/// Retrive all names of submodules.
		std::vector<std::string> GetSubModelsNames() const;
		/// Retrive an operator of the model by name.
		AbstractOperator * GetOperator(std::string name);
		/// Retrive an operator of the model by name.
		const AbstractOperator * GetOperator(std::string name) const;
		/// Retrive all names of operators.
		std::vector<std::string> GetOperatorsNames() const;
		/// Each submodel introduces it's unknowns into the model
		/// so that later it can be accessed
		bool PrepareEntries();
		/// Initialze all entries and submodels.
		bool Initialize();
		/// Initialize data needed for FillResidual.
		/// Called once before nonlinear iterations.
		bool PrepareIterations();
		/// Compute the residual of the model.
		bool FillResidual(Residual & R) const;
		/// Update solution.
		/// alpha is the parameter that scales the update solution.
		/// Usually calculated with UpdateMultiplier.
		bool UpdateSolution(const Sparse::Vector & sol, double alpha);
		/// Move to the next time step
		bool UpdateTimeStep();
		/// Provide new time step.
		bool SetTimeStep(double dt);
		/// Roll back to previous time step
		bool RestoreTimeStep();
		/// Check was the model initialized.
		bool isInitialized() const {return initialized;}
		/// Update variables  contained in all block of automatizator on ghost elements of the grid.
		/// For synchronization of data in individual blocks see AbstractEntry::SynhronizeData.
		void SynchronizeData() { aut.SynchronizeData(); }
		/// Calculate multiplier for update.
		double UpdateMultiplier(const Sparse::Vector & sol) const;
		/// Calculate optimal time step for submodels.
		double AdjustTimeStep(double dt) const;
	};
}

#endif //USE_AUTODIFF && USE_MESH && USE_SOLVER

#endif //INMOST_MODEL_INCLUDED

