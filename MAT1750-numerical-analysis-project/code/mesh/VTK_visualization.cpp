/********************* vtk_visualization.cpp ********************

 * Toolkit for mesh visualization via VTK.
 *
 * Author: Shashwat Sharma
 * Created on: Nov 22, 2018

 ***************************************************************/

#include "VTK_visualization.hpp"



/*! \brief Function to initialize the mesh nodes and triangles for visualization with VTK. Also initializes the dual mesh, if one exists.*/
void VTK_visualization::Initialize (Mesh *_mesh)
{
    mesh = _mesh;

	// Generate and store index lists for exporting the entire mesh
    for (int ii = 0; ii < mesh->Nodes.size(); ii++)
        all_nodeIds.push_back(ii);

    for (int ii = 0; ii < mesh->Triangles.size(); ii++)
        all_triangleIds.push_back(ii);


	if (mesh->dual_mesh)
	{
		for (int ii = 0; ii < mesh->DualTriangles.size(); ii++)
			dual_triangleIds.push_back(ii);
	}

    // Set colours
    VTK_red = {255, 0, 0};
    VTK_green = {0, 255, 0};
    VTK_blue = {0, 0, 255};
	

	return;
	
}


/*! \brief Function to export a set of nodes and triangles to VTK, and optionally a scalar field which resides on nodes (can be passed as NULL to not visualize scalar field).*/
void VTK_visualization::VisualizeField(std::string filename, std::vector<int> &nodeIds, std::vector<int> &triangleIds, double *scalar_field_nodes, double *scalar_field_faces, std::complex<double> *vector_field_nodes, std::complex<double> *vector_field_faces, bool isDual)
{
	
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> faces = vtkSmartPointer<vtkCellArray>::New(); //, pointFaces = vtkSmartPointer<vtkCellArray>::New();
    
    for (unsigned int ii = 0; ii < nodeIds.size(); ii++)
    {
        // Get the mesh nodes for each elementary entity
        vtkIdType id = points->InsertNextPoint(mesh->Nodes[nodeIds[ii]]->p);
        //pointFaces->InsertNextCell(1, &id);
    }

	if (!isDual)
	{
		for (unsigned int ii = 0; ii < triangleIds.size(); ii++)
		{
			// The triangle faces with the 3 points are stored here for VTK/Paraview
			vtkIdType face0 [3];
		
			face0[0] = mesh->Triangles[triangleIds[ii]]->v[0];
			face0[1] = mesh->Triangles[triangleIds[ii]]->v[1];
			face0[2] = mesh->Triangles[triangleIds[ii]]->v[2];
			faces->InsertNextCell(3, face0); //all faces are added to vtk in the correct order   
		}
	}
	else
	{
		for (unsigned int ii = 0; ii < triangleIds.size(); ii++)
		{
			// The triangle faces with the 3 points are stored here for VTK/Paraview
			vtkIdType face0 [3];
		
			face0[0] = mesh->DualTriangles[triangleIds[ii]]->v[0];
			face0[1] = mesh->DualTriangles[triangleIds[ii]]->v[1];
			face0[2] = mesh->DualTriangles[triangleIds[ii]]->v[2];
			faces->InsertNextCell(3, face0); //all faces are added to vtk in the correct order   
		}
	}
	

    // Create PolyData
    vtkSmartPointer<vtkPolyData> polygonPolyData = vtkSmartPointer<vtkPolyData>::New();
    polygonPolyData->SetPoints(points);
    polygonPolyData->SetPolys(faces);
    
    // ------ Apply fields ------

    if (scalar_field_nodes != NULL)
        vtkVisualizeNodeScalarField(polygonPolyData, scalar_field_nodes, filename + " - Node scalar field");
	if (scalar_field_faces != NULL)
		vtkVisualizeFaceScalarField(polygonPolyData, scalar_field_faces, filename + " - Face scalar field");
	if (vector_field_nodes != NULL)
        vtkVisualizeNodeVectorField(polygonPolyData, vector_field_nodes, filename + " - Node vector field");
	if (vector_field_faces != NULL)
		vtkVisualizeFaceVectorField(polygonPolyData, vector_field_faces, filename + " - Face vector field");
	
	
    // Writes the file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

    // std::string fileName = std::string(mesh->MeshFile) + "MeshedVTP.vtp";
    filename += ".vtp";
    const char *cfileName = filename.c_str();
    writer->SetFileName(cfileName);
#if VTK_MAJOR_VERSION <= 5
	writer->SetInput(polygonPolyData);
#else
	writer->SetInputData(polygonPolyData);
#endif
    writer->Write();	 // writes the .vtp file
    std::cout << "Done writing '" << filename << "'" << std::endl;

	return;
	
}


/*! \brief Export a vector field residing on mesh nodes.*/
void VTK_visualization::VisualizeNodeVectorField(vtkSmartPointer<vtkPolyData> polyData, std::complex<double> *vectorField, string fieldName)
{
	
    vtkSmartPointer <vtkDoubleArray> field = vtkSmartPointer <vtkDoubleArray>::New();
    field->SetNumberOfComponents(3);
    field->SetNumberOfTuples(polyData->GetNumberOfPoints());
    field->SetName(fieldName.c_str());
	
    for (int ii = 0; ii < polyData->GetNumberOfPoints(); ii++)
	{
        double fieldValue [3];

        fieldValue[0] = std::real(vectorField[ii*3]);
        fieldValue[1] = std::real(vectorField[ii*3+1]);
		fieldValue[2] = std::real(vectorField[ii*3+2]);
        // field->InsertNextTuple(fieldValue);
		field->SetTuple3(ii, fieldValue[0], fieldValue[1], fieldValue[2]);
    }
	
    polyData->GetPointData()->SetVectors(field); 
    
    return;

}


/*! \brief Export a scalar field residing on mesh nodes.*/
void VTK_visualization::VisualizeNodeScalarField(vtkSmartPointer<vtkPolyData> polyData, double *scalarField, std::string fieldName)
{
    
    vtkSmartPointer <vtkDoubleArray> field = vtkSmartPointer <vtkDoubleArray>::New();
    field->SetNumberOfValues(polyData->GetNumberOfPoints());
    field->SetName(fieldName.c_str());
	
    for (int ii = 0; ii < polyData->GetNumberOfPoints(); ii++)
        field->SetValue(ii, scalarField[ii]);
	
    polyData->GetPointData()->SetScalars(field);

	
    return;
	
}


/*! \brief Export a scalar field residing on mesh faces.*/
void VTK_visualization::VisualizeFaceScalarField(vtkSmartPointer<vtkPolyData> polyData, double *scalarField, std::string fieldName)
{
	
	vtkSmartPointer <vtkDoubleArray> field = vtkSmartPointer <vtkDoubleArray>::New();
	field->SetNumberOfValues(polyData->GetPolys()->GetNumberOfCells());
    field->SetName(fieldName.c_str());
	
    for (int ii = 0; ii < polyData->GetPolys()->GetNumberOfCells(); ii++)
        field->SetValue(ii, scalarField[ii]);
		// field->SetValue(ii, ii);
	
    polyData->GetCellData()->SetScalars(field); 


	return;

}


/*! \brief Export a vector field residing on mesh faces.*/
void VTK_visualization::VisualizeFaceVectorField(vtkSmartPointer<vtkPolyData> polyData, std::complex<double> *vectorField, std::string fieldName)
{
	
    vtkSmartPointer <vtkDoubleArray> field = vtkSmartPointer <vtkDoubleArray>::New();
    field->SetNumberOfComponents(3);
    field->SetNumberOfTuples(polyData->GetPolys()->GetNumberOfCells());
    field->SetName(fieldName.c_str());
	
    for (int ii = 0; ii < polyData->GetPolys()->GetNumberOfCells(); ii++)
	{
		double fieldValue [3];

		fieldValue[0] = std::real(vectorField[ii*3]);
        fieldValue[1] = std::real(vectorField[ii*3+1]);
		fieldValue[2] = std::real(vectorField[ii*3+2]);
        // field->InsertNextTuple(fieldValue);
		field->SetTuple3(ii, fieldValue[0], fieldValue[1], fieldValue[2]);
    }
    polyData->GetCellData()->SetVectors(field); 

	
    return;
	
}


/*! \brief Export a point cloud.*/
void VTK_visualization::VisualizePoints(std::vector<vect3D> &points, std::string fieldName, int colour)
{

    vtkSmartPointer<vtkPoints> vtk_points = vtkSmartPointer<vtkPoints>::New();

    for (int ii = 0; ii < points.size(); ii++)
        vtk_points->InsertNextPoint (points[ii].x, points[ii].y, points[ii].z);

    vtkSmartPointer<vtkPolyData> pointsPolydata = vtkSmartPointer<vtkPolyData>::New();
    pointsPolydata->SetPoints(vtk_points);


    vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    
    #if VTK_MAJOR_VERSION <= 5
        vertexFilter->SetInputConnection(pointsPolydata->GetProducerPort());
    #else
        vertexFilter->SetInputData(pointsPolydata);
    #endif
        vertexFilter->Update();
     
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->ShallowCopy(vertexFilter->GetOutput());

    vtkSmartPointer<vtkUnsignedCharArray> vtk_colours = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(1);
    colors->SetName ("Colours");

    switch colour
    {
        case VTK_RED:
        {
            vtk_colours->InsertNextTupleValue(VTK_red);
            break;
        }
        case VTK_GREEN:
        {
            vtk_colours->InsertNextTupleValue(VTK_green);
            break;
        }
        case VTK_BLUE:
        {
            vtk_colours->InsertNextTupleValue(VTK_blue);
            break;
        }
    }
     
    polydata->GetPointData()->SetScalars(vtk_colours);

    // Write the file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

    filename += ".vtp";
    const char *cfileName = filename.c_str();
    writer->SetFileName(cfileName);
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(polygonPolyData);
#else
    writer->SetInputData(polygonPolyData);
#endif
    writer->Write();     // writes the .vtp file
    std::cout << "Done writing '" << filename << "'" << std::endl;


    return;

}















// Approximates a scalar field given on edge midpoints onto nodes of a mesh by averaging the field values of connected edges 
std::vector<double> VTK_visualization::approximateFieldAtNodes(std::vector<double> *edgeFieldValues)
{
    std::vector<double> approximatedPointField (mesh->Nodes.size()); //store the scalar field approximations at the mesh nodes
    std::vector<double> numEdgesConnected (mesh->Nodes.size()); //store the number of edges connected to each point
    
    // Iterate through the edges and add the field values to the connected nodes, while storing in numEdgesConnected the number of times a certain node has been added to
    for (int ii = 0; ii < mesh->Edges.size(); ii++)
	{
        int node1 = mesh->Edges[ii]->Nodes[0];
        int node2 = mesh->Edges[ii]->Nodes[1];
        approximatedPointField[node1] += (*edgeFieldValues)[ii];
        approximatedPointField[node2] += (*edgeFieldValues)[ii];
        numEdgesConnected[node1] += 1.0;
        numEdgesConnected[node2] += 1.0;
    }
    // Divide each of the node field values by the number of edges connected to get the average field values
    for (int ii = 0; ii < mesh->Nodes.size(); ii++)
	{
        approximatedPointField[ii] /= numEdgesConnected[ii];
    }
    
    return approximatedPointField;
}


// Exports a set of nodes and triangles to VTK, as well as a scalar field which resides on edge midpoints (can be passed as 0 to not visualize scalar field)
void VTK_visualization::VisualizeField(std::vector<int> &nodeIds, std::vector<int> &triangleIds, std::vector<double> *edgeFieldValues)
{
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> faces = vtkSmartPointer<vtkCellArray>::New();
    
    for (unsigned int ii = 0; ii < nodeIds.size(); ii++)
	{
        // Get the mesh nodes for each elementary entity
        points->InsertNextPoint(mesh->Nodes[nodeIds[ii]]->p);
    }
    
    for (unsigned int ii = 0; ii < triangleIds.size(); ii++)
	{
        // The triangle faces with the 3 points are stored here for VTK/Paraview
        vtkIdType face0 [3];
        face0[0] = mesh->Triangles[triangleIds[ii]]->v[0];
        face0[1] = mesh->Triangles[triangleIds[ii]]->v[1];
        face0[2] = mesh->Triangles[triangleIds[ii]]->v[2];
        faces->InsertNextCell(3, face0); //all faces are added to vtk in the correct order 
    }

    // Create PolyData
    vtkSmartPointer<vtkPolyData> polygonPolyData = vtkSmartPointer<vtkPolyData>::New();
    polygonPolyData->SetPoints(points);
    polygonPolyData->SetPolys(faces);
    
    //Colour the scalar field
    if (edgeFieldValues != NULL)
        vtkVisualizeFaceScalarField(polygonPolyData, edgeFieldValues);

    // Writes the file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    std::string fileName = std::string(mesh->MeshFile) + "MeshedVTP.vtp";
    const char *cfileName = fileName.c_str();
    writer->SetFileName(cfileName);
    #if VTK_MAJOR_VERSION <= 5
        writer->SetInput(polygonPolyData);
    #else
    	writer->SetInputData(polygonPolyData);
    #endif
        writer->Write();	 //writes the .vtp file
    cout << "Done writing '" << fileName << "'" << endl;
}

//Find the maximum value in an unsorted list
double VTK_visualization::findMax(std::vector<double> list)
{
    double max = 0.0;
    for (int ii = 0; ii < list.size(); ii++)
	{
        if (list[ii] > max)
            max = list[ii];
    }
    return max;
}


/*! \brief Export a scalar field residing on edge midpoints onto the mesh nodes of a vtkPolyData structure using approximations.*/
void VTK_visualization::VisualizeFaceScalarField(vtkSmartPointer<vtkPolyData> polyData, vector<double> *edgeFieldValues)
{
    //Create the scalar field approximation that will be used at the mesh nodes
    std::vector<double> scalarField = approximateFieldAtNodes(edgeFieldValues);

	
    //for (int ii = 0; ii < scalarField.size(); ii++){
    //    cout << scalarField[ii] << endl;
    //}
    
    //Range for the colour table
    double minScalar = 0, maxScalar = findMax(scalarField); //polyData->GetPolys()->GetNumberOfCells(); //polyData->GetNumberOfPoints()
    std::cout << "Maxscalar: " << maxScalar << std::endl;
    
    // Create the color map
    vtkSmartPointer<vtkLookupTable> colorLookupTable = vtkSmartPointer<vtkLookupTable>::New();
    colorLookupTable->SetTableRange(minScalar, maxScalar);
    colorLookupTable->Build();
    
    //Generate the colors for each point based on the scalar field
    vtkSmartPointer <vtkUnsignedCharArray> colors = vtkSmartPointer <vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetName("Colors");
   
    for (int ii = 0; ii < polyData->GetNumberOfPoints(); ii++){ //polyData->GetNumberOfPoints() //polyData->GetPolys()->GetNumberOfCells()
        double dcolor [3];
        colorLookupTable->GetColor(scalarField[ii], dcolor);
        //std::cout << "dcolor: " << dcolor[0] << " " << dcolor[1] << " " << dcolor[2] << std::endl;
        float color [3];
        for(unsigned int j = 0; j < 3; j++)
		{
            color[j] = static_cast<float>(255.0 * dcolor[j]);
        }
        //cout << "Color " << ii << ": (" << color[0] << "," << color[1] << "," << color[2] << ")" << endl;
        //std::cout << "color: " << (int)color[0] << " " << (int)color[1] << " " << (int)color[2] << std::endl;
        colors->InsertNextTuple(color);
    }
    
    polyData->GetPointData()->SetScalars(colors);
    //polyData->GetCellData()->SetScalars(colors);

	return;
	
}


std::string getNameOfFile (string fileName)
{
    std::size_t found = fileName.find_last_of("/\\");
    return fileName.substr(found+1);
}




// Find the maximum value in an unsorted list
void VTK_visualization::findMaxMin(double *list, int len, double &minval, double &maxval)
{
    maxval = 1.0e-50, minval = 1.0e50;
    for (int ii = 0; ii < len; ii++)
    {
        if (list[ii] > maxval)
            maxval = list[ii];
        if (list[ii] < minval)
            minval = list[ii];
    }
    return;
}


