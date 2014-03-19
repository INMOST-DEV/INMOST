#pragma once
#ifndef INMOST_PARTITIONER_H_INCLUDED
#define INMOST_PARTITIONER_H_INCLUDED
#include "inmost.h"

#if defined(USE_PARTITIONER) && !defined(USE_MESH)
#warning "USE_PARITIONER require USE_MESH"
#undef USE_PARTITIONER
#endif


#if defined(USE_PARTITIONER)

namespace INMOST
{
	class Partitioner
	{
	public:
		enum Type{Zoltan_Parmetis, Zoltan_Scotch, Zoltan_PHG, Zoltan_RCB, Zoltan_RIB, Zoltan_HSFC, Parmetis, Inner_RCM};
		enum Action{Partition, Repartition, Refine};
	private:
		enum Type pt;
		enum Action pa;
		void * pzz;
		Tag weight_tag;
		Mesh * m;
	public:
		static void Initialize(int * argc, char *** argv);
		static void Finalize();
		Partitioner(Mesh * m);
		Partitioner(const Partitioner & other);
		Partitioner & operator =(Partitioner const & other);
		~Partitioner();
		void Evaluate();
		void SetMethod(enum Type t, enum Action a = Repartition);
		void SetWeight(Tag weight);
		void ResetWeight();
		
		Mesh * GetMesh();
		Tag GetWeight();
	};
}
#endif //USE_PARTITIONER


#endif //INMOST_PARTITIONER_H_INCLUDED
