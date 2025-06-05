#include "MAdModel.h"
#include "MeshDataBase.h"
#include "MeshDataBaseInterface.h"
#include "PWLinearSField.h"
#include "LogSizeField.h"
#include "MAdSolution.h"
#include "AnisoMeshSize.h"
#include "MeshQualityManager.h"
#include "MAdQuadrature.h"
#include "CMeshAdapter.h"
#include "MappingInterface.h"

#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkPointLocator.h>
#include <vtkDataArray.h>

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

#ifdef _HAVE_MKL_
#include "mkl_service.h"
#endif

using namespace MAd;

void fixBoundaryEdges(pMesh msh, CMeshAdapter* mshAdp)
{
  // --- travel mesh edges and constraint boundaries
  EIter iteEdg = M_edgeIter(msh);
  while ( pEdge edg = EIter_next(iteEdg) ) {
    if ( edg->nbrFaces() == 1 ) mshAdp->setConstraint(edg);
  }
  EIter_delete(iteEdg);
}

int main(int argc, char* argv[]){
#ifdef _HAVE_MKL_
  atexit(mkl_free_buffers);
#endif
#ifdef _HAVE_MKL_
  atexit(mkl_free_buffers);
#endif
  
  int result = 0;
  
  // initialize empty mesh and solution field
  pGModel mod = NULL;
  GM_create(&mod);
  pMesh  msh = M_new(mod);
  LogSField * logFld = new LogSField(msh);
  
  // load mesh
  std::string mshFil = argv[1];
  M_load(msh, mshFil.c_str());

  // build mesh adapter
  CMeshAdapter *mshAdp = new CMeshAdapter(msh, logFld, false);
  mshAdp->setNbrCollapse(3);
  mshAdp->setNbrSplit(3);
  mshAdp->setNbrSwap(3);
  mshAdp->setMaxIterationsNumber(5);
  // build solution field
  logFld->setMaxGrad(1.5);
  logFld->setMaxAni(1000.);

  // fix boundary edges
  fixBoundaryEdges(msh, mshAdp);
  
  // combine mesh and field
  SolAtVertices *sol = new SolAtVertices(msh, logFld, 1, "Solution");

  // read solution file with vtu loader
  if (argc != 4)
    {
        std::cerr << "Usage: " << argv[0] << " path/to/mesh.msh path/to/sol.vtu path/to/output.msh" << std::endl;
        return EXIT_FAILURE;
    }

  // build vtk reader from the .vtu file
  std::string filename = argv[2];
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  // initialize unstructured grid from reader
  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = reader->GetOutput();
  
  // print some information
  vtkSmartPointer<vtkPoints> points = unstructuredGrid->GetPoints();
  std::cout << "Number of points: " << points->GetNumberOfPoints() << std::endl;
  vtkSmartPointer<vtkCellArray> cells = unstructuredGrid->GetCells();
  std::cout << "Number of cells: " << cells->GetNumberOfCells() << std::endl;
  vtkSmartPointer<vtkPointData> pointData = unstructuredGrid->GetPointData();
  std::cout << "Number of data arrays: " << pointData->GetNumberOfArrays() << std::endl;

  vtkSmartPointer<vtkPointLocator> locator = vtkSmartPointer<vtkPointLocator>::New();
  locator->SetDataSet(unstructuredGrid);
  locator->BuildLocator();

  // associate data to vertex
  VIter iteVer = M_vertexIter(msh);
  vtkSmartPointer<vtkDataArray> array = pointData->GetArray(0);
  if (array->GetNumberOfComponents() != 1) array = pointData->GetArray(1);
  int ver_id = 0;
  while ( pVertex ver = iteVer->next() ) {
    // get coordinates of ver
    double xyz[3];
    V_coord(ver, xyz);
    // find its ID in the vtk datatree
    vtkIdType nearestPointId = locator->FindClosestPoint(xyz);
    double nearestPoint[3];
    points->GetPoint(nearestPointId, nearestPoint);
    // attach value to the vertex
    doubleVector *val = new doubleVector(1);
    (*val)(0) = array->GetTuple1(nearestPointId);
    sol->setValue(ver, val);
    ver_id++;
  }
  
  double fstCmp = logFld->getMeshComplexity();
  
  std::cout << "  initial complexity " << fstCmp << std::endl;
  std::cout << "  final complexity   " << fstCmp << std::endl;

  printf("\n");
  printf("  -----------------------------\n");
  printf("  -         Adaptation        -\n");
  printf("  -----------------------------\n");
  printf("\n");

  logFld->computeAdaptationMet(fstCmp, 2);

  mshAdp->run();

  std::cout << "  dof " << msh->nbrVertices() << std::endl;

  M_writeMsh(msh, argv[3], 2);
  
  if ( !(mshAdp->checkTheMesh(1)) ) {	
    std::cout << "  WARNING - unsuccessful mesh check " << std::endl;
  }
  
  msh->reorderPoints();
  
  double maxAni = 1., meanAni = 0.;
  int idx = 0;
  iteVer->reset();
  while ( pVertex ver = iteVer->next() ) {
    doubleMatrix eigVec(2,2);
    doubleVector eigVal(2);
    MAdMetric met = logFld->findSize(ver)->getMetric();
    met.eig(eigVec, eigVal);
    if ( std::max( eigVal(0) / eigVal(1), eigVal(1) / eigVal(0) ) > maxAni ) {
      idx = ver->getId();
      maxAni = std::max( eigVal(0) / eigVal(1), eigVal(1) / eigVal(0) );
    }
    meanAni += sqrt(std::max( eigVal(0) / eigVal(1), eigVal(1) / eigVal(0) ));
  }
  
  maxAni   = sqrt(maxAni);
  meanAni /= msh->nbrVertices();
  
  std::cout << "  maximum anisotropy " << maxAni << std::endl;
  std::cout << "  idx " << idx << std::endl;
  std::cout << "  average anisotropy " << meanAni << std::endl;
  
  double meanShape, bestShape, worstShape;
  
  logFld->qualityInfo(meanShape, bestShape, worstShape);
  
  if ( meanShape < 0.75 || bestShape < 0.9 || worstShape < 1e-3 ) {
    std::cout << "  WARNING - meanShape " << meanShape << ", bestShape " << bestShape;
    std::cout << ", worstShape " << worstShape << std::endl;
  }
  
  VIter_delete(iteVer);
  
  delete logFld;
  delete sol;
  delete mshAdp;
  GM_delete(mod);
  M_delete(msh);

  return result;
}
