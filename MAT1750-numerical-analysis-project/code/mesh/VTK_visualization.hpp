/********************* VTK_visualization.hpp ********************

 * Toolkit for mesh visualization via VTK.
 *
 * Author: Shashwat Sharma
 * Created on: Nov 22, 2018

 ***************************************************************/

#ifndef VTKVIS_H
#define VTKVIS_H

#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <fstream>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolygon.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>

#include <vtkDataArray.h>;
#include <vtkLookupTable.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>

//#include <vtkPolyDataMapper.h>
//#include <vtkDataArrayTemplate.h>
//#include <vtkActor.h>
//#include <vtkCamera.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderWindow.h>


#include "coordinate_system.hpp"
#include "mesh.hpp"


#define VTK_RED 1
#define VTK_GREEN 2
#define VTK_BLUE 3


class VTK_visualization
{
 public:

	char VTK_red [3];
    char VTK_green [3];
    char VTK_blue [3];

	Mesh *mesh;

	std::vector<int> all_nodeIds, all_triangleIds;
	std::vector<int> dual_nodeIds, dual_triangleIds;
        
	// ====== Functions ======
        
	void Initialize(Mesh *_mesh);

	void VisualizeField(std::string filename, std::vector<int> &nodeIds, std::vector<int> &triangleIds, double *scalar_field_nodes = NULL, double *scalar_field_faces = NULL, std::complex<double> *vector_field_nodes = NULL, std::complex<double> *vector_field_faces = NULL, bool isDual = false);

	void VisualizeNodeScalarField(vtkSmartPointer<vtkPolyData> polyData, double *edgeFieldValues, std::string fieldName);
	void VisualizeFaceScalarField(vtkSmartPointer<vtkPolyData> polyData, double *edgeFieldValues, std::string fieldName);
	
	void VisualizeNodeVectorField(vtkSmartPointer<vtkPolyData> polyData, std::complex<double> *vectorField, std::string fieldName);
	void VisualizeFaceVectorField(vtkSmartPointer<vtkPolyData> polyData, std::complex<double> *vectorField, std::string fieldName);

	void VisualizePoints(std::vector<vect3D> &points, std::string fieldName, int colour = VTK_RED);

	
	double findMax(std::vector<double> list);
		
	std::vector<double> approximateFieldAtNodes(std::vector<double> *edgeFieldValues);
	void VisualizeField(std::vector<int> &nodeIds, std::vector<int> &triangleIds, std::vector<double> *edgeFieldValues);
	void VisualizeFaceScalarField(vtkSmartPointer<vtkPolyData> polyData, std::vector<double> *edgeFieldValues);

	void findMaxMin(double *list, int len, double &minval, double &maxval);
	

};

#endif

