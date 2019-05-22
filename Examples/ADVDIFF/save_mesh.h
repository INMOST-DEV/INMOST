#ifndef _SAVE_MESH_H
#define _SAVE_MESH_H


#include "inmost.h"


/// Specification of output format.
enum OutputMeshFormat
{
	PMF, //< Internal binary format.
	GMV, //< Compatible with General Mesh Viewer and Visit.
	VTK, //< Compatible with Paraview and Visit.
	XML //< Internal user-friendly ascii format based on xml.
};

/// Save provided mesh in chosen format with provided prefix.
/// @param m Mesh to be saved.
/// @param prefix Name of the mesh without format specification.
/// @param format Format of the mesh.
std::string SaveMesh(INMOST::Mesh * m, std::string prefix, OutputMeshFormat format);
/// Save provided mesh in chosen format with provided prefix and counter.
/// @param m Mesh to be saved.
/// @param prefix Name of the mesh without counter and format specification.
/// @param counter The number will be attached to prefix.
/// @param format Format of the mesh.
std::string SaveMesh(INMOST::Mesh * m, std::string prefix, int counter, OutputMeshFormat format);

#endif //_SAVE_MESH_H