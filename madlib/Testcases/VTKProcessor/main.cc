#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPointLocator.h>
#include <vtkDataArray.h>

int main(int argc, char* argv[])
{

    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " path/to/input.vtu path/to/output.vtu" << std::endl;
        return EXIT_FAILURE;
    }

    std::string inputFilename = argv[1];
    std::string outputFilename = argv[2];

    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(inputFilename.c_str());
    reader->Update();

    // initialize unstructured grid from reader
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = reader->GetOutput();
    
    // print some information
    vtkSmartPointer<vtkPoints> points = unstructuredGrid->GetPoints();
    std::cout << "INFO -- number of points: " << points->GetNumberOfPoints() << "\n";
    vtkSmartPointer<vtkPointData> pointData = unstructuredGrid->GetPointData();
    std::cout << "INFO -- number of data arrays: " << pointData->GetNumberOfArrays() << "\n";

    // Get the velocity vector array
    vtkDataArray* velocityArray = pointData->GetVectors("Velocity");
    if (!velocityArray)
    {
        std::cerr << "Error: Could not find vector field named 'Velocity'." << std::endl;
        return EXIT_FAILURE;
    }

    // Create array for velocity magnitude
    vtkNew<vtkDoubleArray> magnitudeArray;
    magnitudeArray->SetName("VelocityMagnitude");
    magnitudeArray->SetNumberOfComponents(1);
    magnitudeArray->SetNumberOfTuples(velocityArray->GetNumberOfTuples());

    for (vtkIdType i = 0; i < velocityArray->GetNumberOfTuples(); ++i)
    {
        double v[3];
        velocityArray->GetTuple(i, v);
        double mag = vtkMath::Norm(v);
        magnitudeArray->SetTuple1(i, mag);
    }

    // Add the new array to the point data
    pointData->AddArray(magnitudeArray);

    // Write to .vtu
    vtkNew<vtkXMLUnstructuredGridWriter> writer;
    writer->SetFileName(outputFilename.c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();

    std::cout << "Successfully wrote velocity magnitude to " << outputFilename << std::endl;

    return EXIT_SUCCESS;
}