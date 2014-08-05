#pragma once
#ifndef INMOST_PARTITIONER_H_INCLUDED
#define INMOST_PARTITIONER_H_INCLUDED
#include "inmost.h"

#if defined(USE_PARTITIONER) && !defined(USE_MESH)
#warning "USE_PARTITIONER require USE_MESH"
#undef USE_PARTITIONER
#endif


#if defined(USE_PARTITIONER)

namespace INMOST
{
	/// Main class to modify or improve the mesh distribution for better load balancing.
	class Partitioner
	{
	public:
		/// Type of the Partitioner can be currently used in this version of INMOST.
		/// @see Zoltan: Parallel Partitioning, Load Balancing and Data-Management Services. http://www.cs.sandia.gov/Zoltan/
		/// @see ParMETIS: Parallel Graph Partitioning and Fill-reducing Matrix Ordering. http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview
		enum Type
		{
			Zoltan_Parmetis, ///< Parmetis partitioner with the Zoltan package interface.
			Zoltan_Scotch,   ///< Scotch partitioner with the Zoltan package interface.
			Zoltan_PHG,      ///< Zoltan topology-based method using Partitioning of HyperGraph.
			Zoltan_RCB,      ///< Zoltan geometry-based method using Recursive Coordinate Bisection.
			Zoltan_RIB,      ///< Zoltan geometry-based method using Recursive Inertial Bisection.
			Zoltan_HSFC,     ///< Zoltan geometry-based method using Hilbert Space-Filling Curve partitioning.
			Parmetis,        ///< Parmetis partitioner with the original interface.
			Inner_RCM        ///< Internal serial only partitioner based on the Reverse Cuthill–McKee algorithm ordering.
		};
		enum Action
		{
			Partition,       ///< Partition "from scratch", not taking into account the current mesh distribution.
			Repartition,     ///< Repartition the existing partition but try to stay close to the current mesh distribution.
			Refine           ///< Refine the current partition assuming only small changes of mesh distribution.
		};
	private:
		enum Type pt;
		enum Action pa;
		void * pzz;
		Tag weight_tag;
		Mesh * m;
	public:
		/// Initialize the use of partitioner.
		/// @param argc The number of arguments transmitted to the function main.
		/// @param argv The pointer to arguments transmitted to the function main.
		/// The shortest call to this function with the default solver parameters is the following: Initialize(NULL,NULL);
		/// @see Partitioner::SetMethod
		/// @see Partitioner::Finalize
		static void Initialize(int * argc, char *** argv);
		/// Finalize the use of partitioner.
		/// @see Partitioner::Initialize
		static void Finalize();
		/// The default constructor of the partitioner for the specified mesh.
		Partitioner(Mesh * m);
		Partitioner(const Partitioner & other);
		Partitioner & operator =(Partitioner const & other);
		~Partitioner();
		/// Evaluate the earlier specified partitioner.
		/// @see Partitioner::SetMethod
		void Evaluate();
		/// Set the partitioner method to be used.
		/// @param t The concrete Type of the partitioner from the selected package.
		/// @param a The partitioner Action, the default is Repartition.
		/// @see Partitioner::Evaluate
		void SetMethod(enum Type t, enum Action a = Repartition);
		/// Compute the specific weights for the selected partitioner.
		void SetWeight(Tag weight);
		/// Reset the computed weights for the partitioner.
		void ResetWeight();
		
		/// Get the Mesh pointer for the current partitioner.
		Mesh * GetMesh();
		/// Get the Tag of the computed weights for the current partitioner.
		Tag GetWeight();
	};
}
#endif //USE_PARTITIONER


#endif //INMOST_PARTITIONER_H_INCLUDED
