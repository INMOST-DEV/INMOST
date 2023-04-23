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
		/// Let the submodel introduce its unknowns
		virtual bool PrepareEntries(Model& P) = 0;
		/// Initialize coupling and dependent unknowns.
		virtual bool Initialize(Model& P) = 0;
		/// Setup coupling with unknowns of other models
		virtual bool SetupCoupling(Model& P) { return true; }
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
		virtual bool SetTimeStep(double dt) {(void)dt; return true;}
		/// Provide current time.
		virtual bool SetTime(double t) {(void)t; return true;}
		/// Roll back to previous step.
		virtual bool RestoreTimeStep() = 0;
		/// Calculate multiplier for update for this model. Can simply return 1.
		virtual double UpdateMultiplier(const Sparse::Vector& sol) const { (void)sol; return 1; }
		/// Calculate time step for this model. Can simply return dt.
		virtual double AdjustTimeStep(double dt) const { return dt; }
		/*
		/// Adapt the data of the model after the mesh refinement/coarsement.
		/// No algorithm by default
		/// If this submodel depends on provided adapted mesh, it should update its data
		virtual void Adaptation(Mesh& m) const { (void)m; }

		virtual void PrepareAdaptation(Mesh& m) { (void)m; }

		//new hanging node was introduced, no connectivity at function call
		virtual void NewNode(Cell& c, Node& n, Storage::reference_array cell_hanging_nodes) { (void)c; (void)n; (void)cell_hanging_nodes; }
		virtual void NewNode(Face& f, Node& n, Storage::reference_array face_hanging_nodes) { (void)f; (void)n; (void)face_hanging_nodes; }
		virtual void NewNode(Edge& e, Node& n) { (void)e; (void)n; }
		//new hanging edge was introduced, no connectivity at function call
		virtual void NewEdge(Cell& c, Edge& e) { (void)c; (void)e; }
		virtual void NewEdge(Face& f, Edge& e) { (void)f; (void)e; }
		//new hanging face was introduced, no connectivity at function call
		virtual void NewFace(Cell& c, Face& f) { (void)c; (void)f; }
		//element was refined into new elements
		virtual void CellRefinement(Cell& old_cell, ElementArray<Cell>& new_cells, ElementSet& new_cell_set) { (void)old_cell; (void)new_cells; (void)new_cell_set; }
		virtual void FaceRefinement(Face& old_face, ElementArray<Face>& new_faces) { (void)old_face; (void)new_faces; }
		virtual void EdgeRefinement(Edge& old_edge, ElementArray<Edge>& new_edges) { (void)old_edge; (void)new_edges; }
		//old elements were connected into an old element
		virtual void CellCoarsening(ElementArray<Cell>& old_cells, Cell& new_cell, ElementSet& old_cells_set) { (void)old_cells; (void)new_cell; (void)old_cells_set; }
		virtual void FaceCoarsening(ElementArray<Face>& old_faces, Face& new_face) { (void)old_faces; (void)new_face; }
		virtual void EdgeCoarsening(ElementArray<Edge>& old_edges, Edge& new_edge) { (void)old_edges; (void)new_edge; }
		*/
		/// Prepare submodel data for output
		virtual void PrepareOutput() {}
		/// Remove submodel data for output
		virtual void FinishOutput() {}
	};

	/// A class to manage residual fill-in of the coupling between models.
	/// Use this to introduce coupling terms that depends on the data between the models.
	/// Coupling term does not introduce its own unknowns and writes to the residual.
	class AbstractCouplingTerm
	{
	public:
		/// Initialize coupling and dependent unknowns.
		virtual bool Initialize(Model& P) = 0;
		/// Setup coupling with unknowns of otheer models
		virtual bool SetupCoupling(Model& P) { return true; }
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
		/// Provide current time.
		virtual bool SetTime(double t) { return true; }
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
		/// Setup coupling with unknowns of otheer models
		virtual bool SetupCoupling(Model& P) { return true; }
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
		/// Provide current time.
		virtual bool SetTime(double t) { return true; }
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
		/// Setup coupling with unknowns of otheer models
		virtual bool SetupCoupling(Model& P) { return true; }
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
		/// Provide current time.
		virtual bool SetTime(double t) { return true; }
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
		std::map< std::string, bool> activeEntries, activeSubModels; ///< Active entries and models

		bool initialized; ///< Indicates whether a model was initialized.
		bool setentries; ///< Indicates whether entries were set.
	public:
		Model(Automatizator & aut) : aut(aut), initialized(false), setentries(false) {}
		//todo:
		//Model(const Model & b) aut(b.aut) {}
		//todo:
		//Model & operator =(Model const & b);
		virtual ~Model() {}
		void ActivateSubModel(std::string name) { activeSubModels[name] = true; }
		void DeactivateSubModel(std::string name) { activeSubModels[name] = false; }
		void ToggleEntryState();
		void ActivateEntry(std::string name) { activeEntries[name] = true; }
		void DeactivateEntry(std::string name) { activeEntries[name] = false; }
		/// Add an entry of block unknowns to a model.
		/// The model stores a link to the entry and may modify its contents.
		/// The intries should be added from Model::Initialize function,
		/// either by model or by any of the submodels.
		void AddEntry(std::string name, AbstractEntry& entry);
		/// Add an entry of block unknowns to a model after an entry with certain name.
		/// The model stores a link to the entry and may modify its contents.
		/// The intries should be added from Model::Initialize function,
		/// either by model or by any of the submodels.
		void AddAfterEntry(std::string name, AbstractEntry& entry, std::string after);
		/// Add an entry of block unknowns to a model as a first entry.
		/// The model stores a link to the entry and may modify its contents.
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
		/// Retrieve an automatizator.
		Automatizator & GetAutomatizator() {return aut;}
		/// Retrieve an automatizator.
		const Automatizator & GetAutomatizator() const {return aut;}
		/// Retrieve a mesh by name.
		Mesh * GetMesh(std::string);
		/// Retrieve a mesh by name.
		const Mesh * GetMesh(std::string) const;
		/// Retrieve all names of meshes.
		std::vector<std::string> GetMeshesNames() const;
		/// Retrieve an entry that describe unknowns of the model by name.
		AbstractEntry * GetEntry(std::string name);
		/// Retrieve an entry that describe unknowns of the model by name.
		const AbstractEntry * GetEntry(std::string name) const;
		/// Retrieve all names of entries.
		std::vector<std::string> GetEntriesNames() const;
		/// Retrieve a submodel of the model by name.
		AbstractSubModel * GetSubModel(std::string name);
		/// Retrieve a submodel of the model by name.
		const AbstractSubModel * GetSubModel(std::string name) const;
		/// Retrieve the name by submodel address
		std::string GetSubModelName(const AbstractSubModel* model) const;
		/// Retrieve all names of submodules.
		std::vector<std::string> GetSubModelsNames() const;
		/// Retrieve a coupling term of the model by name.
		AbstractCouplingTerm* GetCouplingTerm(std::string submodel, std::string name);
		/// Retrieve a coupling term of the model by name.
		const AbstractCouplingTerm* GetCouplingTerm(std::string submodel, std::string name) const;
		/// Retrieve all names of coupling terms.
		std::vector< std::pair< std::string, std::vector<std::string> > > GetCouplingTermNames() const;
		/// Retrieve a scalar function of the model by name.
		AbstractScalarFunction* GetScalarFunction(std::string name);
		/// Retrieve a scalar function of the model by name.
		const AbstractScalarFunction* GetScalarFunction(std::string name) const;
		/// Retrieve all names of scalar functions.
		std::vector<std::string> GetScalarFunctionNames() const;
		/// Retrieve a matrix functon of the model by name.
		AbstractMatrixFunction* GetMatrixFunction(std::string name);
		/// Retrieve a matrix function of the model by name.
		const AbstractMatrixFunction* GetMatrixFunction(std::string name) const;
		/// Retrieve all names of matrix functions.
		std::vector<std::string> GetMatrixFunctionNames() const;
		/// Retrieve an operator of the model by name.
		AbstractOperator * GetOperator(std::string name);
		/// Retrieve an operator of the model by name.
		const AbstractOperator * GetOperator(std::string name) const;
		/// Retrieve all names of operators.
		std::vector<std::string> GetOperatorsNames() const;
		/// Each submodel introduces its unknowns into the model
		/// so that later it can be accessed
		bool PrepareEntries();
		/// Initialize all entries and submodels.
		bool Initialize();
		/// Link to all couplings after PrepareEntries and Initialize
		bool SetupCoupling();
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
		/// Provide current time.
		bool SetTime(double t);
		/// Roll back to previous time step
		bool RestoreTimeStep();
		/// Check if the model was initialized.
		bool isInitialized() const {return initialized;}
		/// Check if the entries were set.
		bool areEntriesSet() const { return setentries; }
		/// Update variables  contained in all block of automatizator on ghost elements of the grid.
		/// For synchronization of data in individual blocks see AbstractEntry::SynchronizeData.
		void SynchronizeData() { aut.SynchronizeData(); }
		/// Calculate multiplier for update.
		double UpdateMultiplier(const Sparse::Vector & sol) const;
		/// Calculate optimal time step for submodels.
		double AdjustTimeStep(double dt) const;
		/*
		/// Adapt the data of the model after the mesh refinement/coarsement.
		/// Those model that use the adapted mesh should update their data
		void Adaptation(Mesh & m) const;
		
		void PrepareAdaptation(Mesh& m);

		void NewNode(Cell& c, Node& n, Storage::reference_array cell_hanging_nodes);
		void NewNode(Face& f, Node& n, Storage::reference_array face_hanging_nodes);
		void NewNode(Edge& e, Node& n);

		void NewEdge(Cell& c, Edge& e);
		void NewEdge(Face& f, Edge& e);

		void NewFace(Cell& c, Face& f);
		
		void CellRefinement(Cell & old_cell, ElementArray<Cell> & new_cells, ElementSet & new_cell_set);
		void FaceRefinement(Face & old_face, ElementArray<Face> & new_faces);
		void EdgeRefinement(Edge & old_edge, ElementArray<Edge> & new_edges);
		
		void CellCoarsening(ElementArray<Cell> & old_cells, Cell & new_cell, ElementSet & old_cells_set);
		void FaceCoarsening(ElementArray<Face> & old_faces, Face & new_face);
		void EdgeCoarsening(ElementArray<Edge> & old_edges, Edge & new_edge);
		*/
		void ReportErrors(const Residual & R) const;
		void PrepareOutput();
		void FinishOutput();
	};
}

#endif //USE_AUTODIFF && USE_MESH && USE_SOLVER

#endif //INMOST_MODEL_INCLUDED

