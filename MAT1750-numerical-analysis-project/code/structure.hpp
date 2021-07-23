/************************ structure.hpp ************************

 * Storage and processing of all geometric and electrical
 * information pertaining to the structure being simulated.
 *
 * Author: Shashwat Sharma
 * Created on: Nov 22, 2018

 ***************************************************************/

#ifndef STRUCTURE_H
#define STRUCTURE_H

#include <complex>
#include <vector>
#include <map>
#include <numeric>
#include <utility>
#include <set>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <algorithm>

#include "coordinate_system.hpp"
#include "mesh.hpp"
// #include "layer_manager.hpp"
// #include "GDS_processor.hpp"

#define MSH 11
#define GDS 12


/*! \brief Parent structure class.*/
class Structure
{
 public:

	// ====== Geometry and mesh data ======

 	// Main structure
 	Mesh *mesh;

 	// Objects in the structure
 	std::vector<Mesh *> objects;

 	// Number of Objects
 	int N_obj;

 	// Base mesh characteristic length, applicable for GDS files meshed with Gmsh
	int lc = 1000;

	// Switch to toggle whether the dual mesh should also be generated
	bool dual_mesh = false;
	

	// ====== Mesh import ======

	std::string MeshFile;			///< Mesh file name
	std::string TechFile;        	///< Tech file name
	std::string PortFile;        	///< Port definition file name

	int mesh_filetype = MSH;        ///< File-type through which the mesh is to be imported

	void PickMeshImporter();
	void ExtractMeshData_Msh();						/// Read mesh from Gmsh mesh file (.msh)
	void ExtractMeshData_Gmsh();                 	/// Read mesh from processed gmsh object
	void ImportFromGDS();							/// Run GDS importer and mesher


	// ====== Mesh processing ======
	

	void TestForContacts_Conformal();
	void TestForContacts_NonConformal();

	
};


#endif
