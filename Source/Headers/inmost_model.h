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
		virtual bool PrepareEntries(Model& P) = 0;
		/// Initialize coupling and dependent unknowns.
		virtual bool Initialize(Model& P) = 0;
		/// Initialize data needed for FillResidual.
		/// Called once before nonlinear iterations.
		virtual bool PrepareIterations() { return true; }
		/// Fill part of the residual related to my unknowns.
		virtual bool FillResidual(Residual& R) const = 0;
		/// Update solution.
		virtual bool UpdateSolution(const Sparse::Vector& sol, double alpha) = 0;
		/// Update time step.
		virtual bool UpdateTimeStep() = 0;
		/// Provide time step.
		virtual bool SetTimeStep(double dt) = 0;
		/// Roll back to previous step.
		virtual bool RestoreTimeStep() = 0;
		/// Calculate multiplier for update for this model. Can simply return 1.
		virtual double UpdateMultiplier(const Sparse::Vector& sol) const { (void)sol; return 1; }
		/// Calculate time step for this model. Can simply return dt.
		virtual double AdjustTimeStep(double dt) const { return dt; }
		/// Adapt the data of the model after the mesh refinement/coarsement.
		/// No algorithm by default
		/// If this submodel depends on provided adapted mesh, it should update it's data
		virtual void Adaptation(Mesh& m) const { (void)m; }

		virtual void PrepareAdaptation(Mesh& m) { (void)m; }

		virtual void CellRefinement(Cell& old_cell, ElementArray<Cell>& new_cells) { (void)old_cell; (void)new_cells; }
		virtual void FaceRefinement(Face& old_face, ElementArray<Face>& new_faces) { (void)old_face; (void)new_faces; }
		virtual void EdgeRefinement(Edge& old_edge, ElementArray<Edge>& new_edges) { (void)old_edge; (void)new_edges; }

		virtual void CellCoarsening(ElementArray<Cell>& old_cells, Cell& new_cell) { (void)old_cells; (void)new_cell; }
		virtual void FaceCoarsening(ElementArray<Face>& old_faces, Face& new_face) { (void)old_faces; (void)new_face; }
		virtual void EdgeCoarsening(ElementArray<Edge>& old_edges, Edge& new_edge) { (void)old_edges; (void)new_edge; }
	};

	/// A class to manage residual fill-in of the coupling between models.
	/// Use this to introduce coupling terms that depends on the data between the models.
	/// Coupling term does not introduce it's own unknowns and writes to the residual.
	class AbstractCouplingTerm
	{
	public:
		/// Initialize coupling and dependent unknowns.
		virtual bool Initialize(Model& P) = 0;
		/// Fill part of the residual related to my unknowns.
		virtual bool FillResidual(Residual& R) const = 0;
		
		///
		/// This functions might be not necessary:
		/// 
		
		/// Initialize data needed for FillResidual.
		/// Called once before nonlinear iterations.
		virtual bool PrepareIterations() { return true; }
		/// Update time step.
		virtual bool UpdateTimeStep() { return true; }
		/// Provide time step.
		virtual bool SetTimeStep(double dt) { return true; }
		/// Roll back to previous step.
		virtual bool RestoreTimeStep() { return true; }
		/// Calculate multiplier for update for this model. Can simply return 1.
		virtual double UpdateMultiplier(const Sparse::Vector& sol) const { (void)sol; return 1; }
		/// Calculate time step for this model. Can simply return dt.
		virtual double AdjustTimeStep(double dt) const { return dt; }
	};

	/// This can be used to abstract implementation of particular scalar functions, 
	/// that depend on mesh data and variables of the models.
	class AbstractScalarFunction
	{
	public:
		/// Initialize coupling and dependent unknowns.
		virtual bool Initialize(Model& P) = 0;
		/// The main function provided by operator.
		/// @param e Element on which operator is evaluated.
		/// @return scalar value of the coupling function.
		virtual variable Evaluate(const Storage& e) const = 0;

		///
		/// This functions might be not necessary:
		/// 
		
		/// Initialize data needed for FillResidual.
		/// Called once before nonlinear iterations.
		virtual bool PrepareIterations() { return true; }
		/// Update time step.
		virtual bool UpdateTimeStep() { return true; }
		/// Provide time step.
		virtual bool SetTimeStep(double dt) { return true; }
		/// Roll back to previous step.
		virtual bool RestoreTimeStep() { return true; }
		/// Calculate multiplier for update for this function. Can simply return 1.
		virtual double UpdateMultiplier(const Sparse::Vector& sol) const { (void)sol; return 1; }
		/// Calculate time step for this model. Can simply return dt.
		virtual double AdjustTimeStep(double dt) const { return dt; }
	};

	/// This can be used to abstract implementation of particular matrix functions, 
	/// that depend on mesh data and variables of the models.
	class AbstractMatrixFunction
	{
	public:
		/// Initialize coupling and dependent unknowns.
		virtual bool Initialize(Model& P) = 0;
		/// The main function provided by operator.
		/// @param e Element on which operator is evaluated.
		/// @return matrix value of the coupling function.
		virtual vMatrix Evaluate(const Storage& e) const = 0;

		///
		/// This functions might be not necessary:
		/// 

		/// Initialize data needed for FillResidual.
		/// Called once before nonlinear iterations.
		virtual bool PrepareIterations() { return true; }
		/// Update time step.
		virtual bool UpdateTimeStep() { return true; }
		/// Provide time step.
		virtual bool SetTimeStep(double dt) { return true; }
		/// Roll back to previous step.
		virtual bool RestoreTimeStep() { return true; }
		/// Calculate multiplier for update for this function. Can simply return 1.
		virtual double UpdateMultiplier(const Sparse::Vector& sol) const { (void)sol; return 1; }
		/// Calculate time step for this model. Can simply return dt.
		virtual double AdjustTimeStep(double dt) const { return dt; }
	};
	
	/// A class to organize a model.
	class Model
	{
		Automatizator & aut; ///< Automatizator that is used to manage all unknowns of the model.
		//todo: decide later on how it should be stored
		std::vector< std::pair<std::string, AbstractSubModel *> > SubModels; ///< A set of submodels of the model.
		std::vector< std::pair<std::string, AbstractEntry *> > Entries; ///< A set of entries for blocks of unknowns of the model. The naming convention is "SubModelName:EntryName"
		std::vector< std::pair<std::string, Mesh *> > Meshes; ///< A set of meshes of the model.
		std::map< AbstractSubModel*, std::vector< std::pair<std::string, AbstractCouplingTerm*> > > CouplingTerms; ///< A set of coupling terms, tied to each submodel
		//std::map< AbstractSubModel*, std::vector< std::pair<std::string, AbstractEntry*> > > SubModelEntries;
		std::vector< std::pair<std::string, AbstractScalarFunction*> > ScalarFunctions; ///< A set of scalar functions, present in the model
		std::vector< std::pair<std::string, AbstractMatrixFunction*> > MatrixFunctions; ///< A set of matrix functions, present in the model
		std::vector< std::pair<std::string, AbstractOperator *> > Operators; ///< A set of operators used by submodels
		bool initialized; ///< Indicates whether the model was initialized.
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
		void AddEntry(std::string name, AbstractEntry& entry);
		/// Add an entry of block unknowns to a model after an entry with certain name.
		/// The model stores a link to the entry and may modify it contents.
		/// The intries should be added from Model::Initialize function,
		/// either by model or by any of the submodels.
		void AddAfterEntry(std::string name, AbstractEntry& entry, std::string after);
		/// Add an entry of block unknowns to a model as a first entry.
		/// The model stores a link to the entry and may modify it contents.
		/// The intries should be added from Model::Initialize function,
		/// either by model or by any of the submodels.
		void AddFirstEntry(std::string name, AbstractEntry& entry);
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
		/// \warning A pointer to the original object is stored, thus do not destroy the original object.
		void AddSubModel(std::string name, AbstractSubModel & submodel);
		/// Add coupling term to a submodel.
		/// An empty submodel name will add a global coupling term.
		/// Coupling terms linked to submodels are not evaluated for non-active submodels.
		/// Global coupling terms are always evaluated and should monitor 
		/// the state of the unknown to avoid writing to respective residual.
		/// \warning A pointer to the original object is stored, thus do not destroy the original object.
		void AddCouplingTerm(std::string submodel, std::string name, AbstractCouplingTerm& term);
		/// Add scalar coupling function to a model.
		/// \warning A pointer to the original object is stored, thus do not destroy the original object.
		void AddScalarFunction(std::string name, AbstractScalarFunction& func);
		/// Add matrix coupling function to a model.
		/// \warning A pointer to the original object is stored, thus do not destroy the original object.
		void AddMatrixFunction(std::string name, AbstractMatrixFunction& func);
		/// Add an operator to a model.
		/// \warning A pointer to the original object is stored, thus do not destroy the original object.
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
		/// Retrive the name by submodel address
		std::string GetSubModelName(const AbstractSubModel* model) const;
		/// Retrive all names of submodules.
		std::vector<std::string> GetSubModelsNames() const;
		/// Retrive a coupling term of the model by name.
		AbstractCouplingTerm* GetCouplingTerm(std::string submodel, std::string name);
		/// Retrive a coupling term of the model by name.
		const AbstractCouplingTerm* GetCouplingTerm(std::string submodel, std::string name) const;
		/// Retrive all names of coupling terms.
		std::vector< std::pair< std::string, std::vector<std::string> > > GetCouplingTermNames() const;
		/// Retrive a scalar function of the model by name.
		AbstractScalarFunction* GetScalarFunction(std::string name);
		/// Retrive a scalar function of the model by name.
		const AbstractScalarFunction* GetScalarFunction(std::string name) const;
		/// Retrive all names of scalar functions.
		std::vector<std::string> GetScalarFunctionNames() const;
		/// Retrive a matrix functon of the model by name.
		AbstractMatrixFunction* GetMatrixFunction(std::string name);
		/// Retrive a matrix function of the model by name.
		const AbstractMatrixFunction* GetMatrixFunction(std::string name) const;
		/// Retrive all names of matrix functions.
		std::vector<std::string> GetMatrixFunctionNames() const;
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
		/// Adapt the data of the model after the mesh refinement/coarsement.
		/// Those model that use the adapted mesh should update their data
		void Adaptation(Mesh & m) const;
		
		void PrepareAdaptation(Mesh& m);
		
		void CellRefinement(Cell & old_cell, ElementArray<Cell> & new_cells);
		void FaceRefinement(Face & old_face, ElementArray<Face> & new_faces);
		void EdgeRefinement(Edge & old_edge, ElementArray<Edge> & new_edges);
		
		void CellCoarsening(ElementArray<Cell> & old_cells, Cell & new_cell);
		void FaceCoarsening(ElementArray<Face> & old_faces, Face & new_face);
		void EdgeCoarsening(ElementArray<Edge> & old_edges, Edge & new_edge);
		
		void ReportErrors(const Residual & R) const;
	};
}

#endif //USE_AUTODIFF && USE_MESH && USE_SOLVER

#endif //INMOST_MODEL_INCLUDED

