#include "CMeshAdapter.h"
#include "LogSizeField.h"
#include "MAdSolution.h"

#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkPointLocator.h>
#include <vtkDataArray.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

#ifdef _HAVE_MKL_
#include "mkl_service.h"
#endif

using namespace MAd;


void printHelp() {
    std::cout << "Usage: program_name [options]\n";
    std::cout << "Options:\n";
    std::cout << "  --help, -h          show this help message\n";
    std::cout << "  -m, --mesh <val>  input mesh (required)\n";
    std::cout << "  -f, --filelist <val>  list of solution files (required)\n";
    std::cout << "  -o, --output <val>  path to the adapted mesh (required)\n";
    std::cout << "  -c, --cmp <val>  complexity (default='-1')\n";
    std::cout << "  -q, --field <val>  field used for adaptation (default='VelocityMagnitude')\n";
    std::cout << "  -grad, --maxgrad <val>  max gradient (default='1.5')\n";
    std::cout << "  -ani, --maxani <val>  max anisotropy (default='1000.')\n";
    std::cout << "  -xlen, --maxlen <val>  max length (default='-1')\n";
    std::cout << "  -nlen, --minlen <val>  min length (default='-1')\n";
    std::cout << "  -met, --metric <val>  metric combination method"
                 " {'intersection', 'mean'} (default='intersection')\n";
}

//std::pair<bool, std::unordered_map<std::string, std::string>> parseArgs(int argc, char* argv[]) {
std::unordered_map<std::string, std::string> parseArgs(int argc, char* argv[]) {
  /*
  Parses the input arguments and returns the extracted options.
  */
  std::unordered_map<std::string, std::string> options;

  // required options
  bool mesh_bool = false, file_bool = false, out_bool = false;
  // default options
  options["cmp"] = "-1";
  options["field"] = "VelocityMagnitude";
  options["maxgrad"] = "1.5";
  options["maxani"] = "1000";
  options["maxlength"] = "-1";
  options["minlength"] = "-1";
  options["metric"] = "intersection";
  // arg parsing
  for (int i = 1; i < argc; i++) {
    std::string arg = argv[i];
    if ((arg == "-h" || arg == "--help")) {
      throw std::invalid_argument("INFO -- show help message\n\n");
    }
    else if ((arg == "-m" || arg == "--mesh") && i + 1 < argc) {
      options["mesh"] = argv[++i];
      std::cout << "INFO -- mesh file: " << options["mesh"] << "\n";
      mesh_bool = true;
    }
    else if ((arg == "-f" || arg == "--filelist") && i + 1 < argc) {
      options["filelist"] = argv[++i];
      std::cout << "INFO -- file list: " << options["filelist"] << "\n";
      file_bool = true;
    }
    else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
      options["output"] = argv[++i];
      std::cout << "INFO -- output: " << options["output"] << "\n";
      out_bool = true;
    }
    else if ((arg == "-c" || arg == "--cmp") && i + 1 < argc) {
      options["cmp"] = argv[++i];
      std::cout << "INFO -- complexity: " << options["cmp"] << "\n";
    }
    else if ((arg == "-q" || arg == "--field") && i + 1 < argc) {
      options["field"] = argv[++i];
      std::cout << "INFO -- field: " << options["field"] << "\n";
    }
    else if ((arg == "-grad" || arg == "--maxgrad") && i + 1 < argc) {
      options["maxgrad"] = argv[++i];
      std::cout << "INFO -- max gradient: " << options["maxgrad"] << "\n";
    }
    else if ((arg == "-ani" || arg == "--maxani") && i + 1 < argc) {
      options["maxani"] = argv[++i];
      std::cout << "INFO -- max anisotropy: " << options["maxani"] << "\n";
    }
    else if ((arg == "-xlen" || arg == "--maxlength") && i + 1 < argc) {
      options["maxlength"] = argv[++i];
      std::cout << "INFO -- max length: " << options["maxlength"] << "\n";
    }
    else if ((arg == "-nlen" || arg == "--minlength") && i + 1 < argc) {
      options["minlength"] = argv[++i];
      std::cout << "INFO -- min length: " << options["minlength"] << "\n";
    }
    else if ((arg == "-met" || arg == "--metric") && i + 1 < argc) {
      options["metric"] = argv[++i];
      if (options["metric"] != "intersection" && options["metric"] != "mean")
        options["metric"] = "intersection";
      std::cout << "INFO -- metric: " << options["metric"] << "\n";
    }
    else {
      std::ostringstream errorMsg;
      errorMsg<< "ERRROR -- unknown argument: " << arg << "\n\n";
      throw std::invalid_argument(errorMsg.str());
    }
  }
  if (mesh_bool && file_bool && out_bool)
    return options;
  else {
    throw std::invalid_argument("ERROR -- missing required arguments\n\n");
  }

}


LogSField* buildLogSField(pMesh &msh, double maxgrad, double maxani, double maxlength, double minlength) {
  /*
  Creates and returns a LogSField object.
  */
  LogSField *logFld = new LogSField(msh);
  // build solution field
  logFld->setMaxGrad(maxgrad);
  logFld->setMaxAni(maxani);
  if (maxlength > 0.)
    logFld->setMaxLength(maxlength);
  if (minlength > 0.)
    logFld->setMinLength(minlength);
  return logFld;
}


CMeshAdapter* buildCMeshAdapter(pMesh &msh, LogSField *logFld) {
  /*
  Creates and returns a CMeshAdapter object.
  */
  bool curved = (msh->getOrder() > 1);
  CMeshAdapter *mshAdp = new CMeshAdapter(msh, logFld, curved);
  mshAdp->setNbrCollapse(3);
  mshAdp->setNbrSplit(3);
  mshAdp->setNbrSwap(3);
  mshAdp->setMaxIterationsNumber(5);
  return mshAdp;
}


std::tuple<vtkSmartPointer<vtkPoints>, vtkSmartPointer<vtkDataArray>, vtkSmartPointer<vtkPointLocator>> buildVtkObjects(std::string filename, std::string field) {
  /*
  Builds several vtk objects (points, array and locator). 
  */
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  // initialize unstructured grid from reader
  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = reader->GetOutput();
  
  // print some information
  vtkSmartPointer<vtkPoints> points = unstructuredGrid->GetPoints();
  std::cout << "INFO -- number of points: " << points->GetNumberOfPoints() << "\n";
  vtkSmartPointer<vtkPointData> pointData = unstructuredGrid->GetPointData();
  std::cout << "INFO -- number of data arrays: " << pointData->GetNumberOfArrays() << "\n";

  vtkSmartPointer<vtkDataArray> array;
  array = pointData->GetArray(field.c_str());

  vtkSmartPointer<vtkPointLocator> locator = vtkSmartPointer<vtkPointLocator>::New();
  locator->SetDataSet(unstructuredGrid);
  locator->BuildLocator();

  return std::make_tuple(points, array, locator);
}


void combineMeshSolAtVertices(std::string filename, std::string field, pMesh &msh, SolAtVertices *sol) {
  /*
  Combines solution at mesh points.
  */
  vtkSmartPointer<vtkPoints> points;
  vtkSmartPointer<vtkDataArray> array;
  vtkSmartPointer<vtkPointLocator> locator;

  std::tie(points, array, locator) = buildVtkObjects(filename, field);

  // associate data to vertex
  VIter iteVer = M_vertexIter(msh);
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
  }
}


void combineMeshSolAtElements(std::string filename, std::string field, pMesh &msh, SolAtElements *sol, const double* geo) {
  /*
  Combines solution at mesh points (HO version).
  */
  vtkSmartPointer<vtkPoints> points;
  vtkSmartPointer<vtkDataArray> array;
  vtkSmartPointer<vtkPointLocator> locator;

  std::tie(points, array, locator) = buildVtkObjects(filename, field);

  // associate data to vertex
  FIter iteFac = M_faceIter(msh);
  int nbrNod = ( msh->getOrder() + 1 ) * ( msh->getOrder() + 2 ) / 2;
  double* val = new double[nbrNod];
  while ( pFace fac = iteFac->next() ) {
    double xyz[3];
    pVertex facVer[3];
    for (int iVer = 0; iVer < 3; iVer++ ) {
      facVer[iVer] = fac->getVertex(iVer);
    }
    for (int iVer = 0; iVer < nbrNod; iVer++ ) {
      // get coordinates of HO nodes
      xyz[0] = geo[3*iVer] * facVer[0]->seeCoord()[0] + geo[3*iVer+1] * facVer[1]->seeCoord()[0] + geo[3*iVer+2] * facVer[2]->seeCoord()[0];
      xyz[1] = geo[3*iVer] * facVer[0]->seeCoord()[1] + geo[3*iVer+1] * facVer[1]->seeCoord()[1] + geo[3*iVer+2] * facVer[2]->seeCoord()[1];
      xyz[2] = geo[3*iVer] * facVer[0]->seeCoord()[2] + geo[3*iVer+1] * facVer[1]->seeCoord()[2] + geo[3*iVer+2] * facVer[2]->seeCoord()[2];

      vtkIdType nearestPointId = locator->FindClosestPoint(xyz);
      double nearestPoint[3];
      // cellCentersOutput->GetPoint(nearestPointId, nearestPoint);
      points->GetPoint(nearestPointId, nearestPoint);

      // attach value to the vertex
      val[iVer] = array->GetTuple1(nearestPointId);
    }
    sol->setValue(fac, val);
  }
}


void updateAdaptationMetIntersection(pMesh &refMsh, LogSField *logFld, pMesh &tmpMsh, LogSField *tmpLogFld) {
  /*
  Intersects the metrics of the given meshes and updates the reference metric.
  */
  pVertex ver, tmpVer;
  VIter iteVer = M_vertexIter(refMsh);
  VIter tmpIteVer = M_vertexIter(tmpMsh);
  while ( (ver = iteVer->next()) && (tmpVer =  tmpIteVer->next())) {
    MAdMetric met = logFld->findSize(ver)->getMetric();
    MAdMetric tmpMet = tmpLogFld->findSize(tmpVer)->getMetric();
    MAdMetric newMet;
    if (M_dim(refMsh) == 2) {
      newMet = intersection2D(met, tmpMet);
    }
    else {
      newMet = intersection(met, tmpMet);
    }
    // update metric
    AnisoMeshSize *siz = (AnisoMeshSize*) logFld->findSize(ver);
    siz->setMetric(newMet);
  }
}


void updateAdaptationMetMean(pMesh &refMsh, LogSField *logFld, pMesh &tmpMsh, LogSField *tmpLogFld, int ite) {
  /*
  Computes the iterative mean of the metrics of the given meshes and updates the reference metric.
  Note: ite is the number of samples used to compute the mean:
        see P. Pébay (2008): https://www.osti.gov/servlets/purl/1028931
  */
  pVertex ver, tmpVer;
  VIter iteVer = M_vertexIter(refMsh);
  VIter tmpIteVer = M_vertexIter(tmpMsh);
  while ( (ver = iteVer->next()) && (tmpVer =  tmpIteVer->next())) {
    MAdMetric met = logFld->findSize(ver)->getMetric();
    MAdMetric tmpMet = tmpLogFld->findSize(tmpVer)->getMetric();
    // get metric matrices
    double matMet[3][3], tmpMatMet[3][3];
    met.getMat(matMet);
    tmpMet.getMat(tmpMatMet);
    // iterative mean
    for (int i=0;i<3;i++){
      for (int j=0;j<3;j++){
        matMet[i][j] += (tmpMatMet[i][j] - matMet[i][j]) / ite;
      }
    }
    // update metric
    met.setMat(matMet);
    AnisoMeshSize *siz = (AnisoMeshSize*) logFld->findSize(ver);
    siz->setMetric(met);
  }
}


double computeErrorBound(pMesh &msh, LogSField *logFld) {
  /*
  Computes and returns an error estimate based on the mesh complexity, order and dimension:
  see equation 5.1 of O. Coulaud et al. (2024): https://doi.org/10.1016/j.jcp.2024.112774
  */
  int order = M_maxOrder(msh);
  std::cout << "INFO -- mesh order: " << order << "\n";
  int dim =  M_dim(msh);
  std::cout << "INFO -- mesh dim: " << dim << "\n";
  double cmp = logFld->getMeshComplexity();
  std::cout << "INFO -- mesh cmp: " << cmp << "\n";
  return pow(cmp, -(order + 1) / dim);
}


double computeError(pMesh &msh, LogSField *logFld) {
  /*
  Computes and returns an error estimate based on the metric Frobenius norm.
  */
  double totErr = 0;
  pVertex ver;
  VIter iteVer = M_vertexIter(msh);
  while ( (ver = iteVer->next())) {
    MAdMetric met = logFld->findSize(ver)->getMetric();
    totErr += met.frobNormSq();
  }
  return sqrt(totErr);
}


double computeMinEdge(pMesh &msh, LogSField *logFld) {
  /*
  Computes and returns the minimal edge length of the mesh.
  */
  double minEdgeLength = __DBL_MAX__, edgeLength;
  EIter iteEdg = M_edgeIter(msh);
  while ( pEdge edg = EIter_next(iteEdg) ) {
    edgeLength = logFld->SF_VV_lengthSq(edg->getVertex(0), edg->getVertex(1));
    if (edgeLength < minEdgeLength)
      minEdgeLength = edgeLength;
  }
  return sqrt(minEdgeLength);
}


int main(int argc, char* argv[]) {
  #ifdef _HAVE_MKL_
    atexit(mkl_free_buffers);
  #endif
  #ifdef _HAVE_MKL_
    atexit(mkl_free_buffers);
  #endif

  // parse options
  std::unordered_map<std::string, std::string> options;
  try {
    options = parseArgs(argc, argv);
  }
  catch(const std::exception& e) {
    std::cout << e.what();
    printHelp();
    return EXIT_FAILURE;
  }

  // read solution files list
  int nb_lines = 0;
  std::string line;
  std::ifstream myfile(options["filelist"]);
  std::vector<std::string> files;

  while (std::getline(myfile, line)) {
    files.push_back(line);
    ++nb_lines;
  }

  if (nb_lines <1) {
    std::cout << "ERROR -- empty files list\n";
    return EXIT_FAILURE;
  }
  std::cout << "INFO -- number of files: " << nb_lines << "\n";
  
  // initialize empty meshes
  pGModel mod = NULL;
  GM_create(&mod);
  pMesh  msh = M_new(mod);

  // load mesh
  std::string mshFil = options["mesh"];
  M_load(msh, mshFil.c_str());
  
  // first solution
  std::cout << "INFO -- processing solution 0: " << files[0] << "..\n";

  // mesh adapter options
  double cmp = stof(options["cmp"]);
  double maxgrad = stof(options["maxgrad"]);
  double maxani = stof(options["maxani"]);
  double maxlength = stof(options["maxlength"]);
  double minlength = stof(options["minlength"]);

  // build mesh adapter and logSField
  LogSField *logFld = buildLogSField(msh, maxgrad, maxani, maxlength, minlength);
  CMeshAdapter *mshAdp = buildCMeshAdapter(msh, logFld);
  std::cout << "DEBUG -- built logsfield and cmeshadapter\n";

  // build sol
  double *geo = nullptr;
  SolAtVertices *solAV = nullptr;
  SolAtElements *solAE = nullptr;
  int ord = msh->getOrder();
  if (ord == 1) {
    std::cout << "INFO -- order: " << ord << " => Hessian based adaptation\n";
    solAV = new SolAtVertices(msh, logFld, 1, "Solution");
    std::cout << "DEBUG -- built solAtVertices\n";
    combineMeshSolAtVertices(files[0], options["field"], msh, solAV);
    std::cout << "DEBUG -- combined mesh and sol\n";
  }
  else {
    int nbrNod = ( ord + 1 ) * ( ord + 2 ) / 2;
    geo = new double[3*nbrNod];
    memset(geo, 0, sizeof(double)*3*nbrNod);
    geo[0] = 1.;
    geo[4] = 1.;
    geo[8] = 1.;
    
    int iVer = 3;
    for (int i = 0; i <= ord-1; i++ ) {
      for (int j = 0; j <= ord-i; j++ ) {
        int k = ord-i-j;
        if ( j == ord  || k == ord ) continue;
        geo[3*iVer  ] = (double) i / (double) ord;
        geo[3*iVer+1] = (double) j / (double) ord;
        geo[3*iVer+2] = (double) k / (double) ord;
        iVer++;
      }
    }
    solAE = new SolAtElements(msh, logFld, 1, ord, geo, "Solution");
    std::cout << "DEBUG -- built solAtElements\n";
    combineMeshSolAtElements(files[0], options["field"], msh, solAE, geo);
    std::cout << "DEBUG -- combined mesh and sol\n";
  }

  // compute complexity if not imposed
  if (cmp < 0) {
    cmp = logFld->getMeshComplexity();
    std::string cmp_file = options["output"].substr(0, options["output"].find_last_of("/")) \
                           + "/cmp.txt";
    std::ofstream cmpFil(cmp_file, std::ios::trunc);
    cmpFil << cmp << std::endl;
    std::cout << "INFO -- complexity: " << cmp << " printed to: " << cmp_file << "\n";
  }

  // compute reference metric
  logFld->computeAdaptationMet(cmp, 2);

  // compute error estimate
  std::string error_file = options["output"].substr(0, options["output"].find_last_of("/")) \
                           + "/error.txt";
  std::ofstream errFil(error_file, std::ios::trunc);
  double totErr = computeError(msh, logFld);
  errFil << totErr << std::endl;
  std::cout << "INFO -- error " << totErr << " printed to: " << error_file << "\n";
  totErr = computeErrorBound(msh, logFld);

  // for each next solution file:
  // 0. build temporary mesh adapter and logSField
  // 1. combine temporary mesh and solution
  // 2. compute temporary metric
  // 3. update the reference metric
  for (int ii=1; ii<nb_lines; ii++) {
    std::cout << "INFO -- processing solution " << ii << ": " << files[ii] << "..\n";
    // create temporary mesh
    pGModel tmpMod = NULL;
    GM_create(&tmpMod);
    pMesh tmpMsh = M_new(tmpMod);
    M_load(tmpMsh, mshFil.c_str());
    // build temporary mesh adapter and logSField
    LogSField *tmpLogFld = buildLogSField(tmpMsh, maxgrad, maxani, maxlength, minlength);
    // combine tomparary mesh and solution
    SolAtVertices *tmpSolAV = nullptr;
    SolAtElements *tmpSolAE = nullptr;
    if (ord == 1) {
      tmpSolAV = new SolAtVertices(tmpMsh, tmpLogFld, 1, "tmpSolution");
      combineMeshSolAtVertices(files[ii], options["field"], tmpMsh, tmpSolAV);
      std::cout << "INFO -- combined temporary mesh and solution\n";
    }
    else {
      tmpSolAE = new SolAtElements(tmpMsh, tmpLogFld, 1, ord, geo, "tmpSolution");
      combineMeshSolAtElements(files[ii], options["field"], tmpMsh, tmpSolAE, geo);
      std::cout << "INFO -- combined temporary mesh and solution\n";
    }
    // compute metric
    tmpLogFld->computeAdaptationMet(cmp, 2);
    std::cout << "INFO -- computed temporary metric\n";
    // update metric (iterative mean or intersection)
    if (options["metric"] == "intersection")
      updateAdaptationMetIntersection(msh, logFld, tmpMsh, tmpLogFld);
    else
      updateAdaptationMetMean(msh, logFld, tmpMsh, tmpLogFld, ii + 1);
    std::cout << "INFO -- updated reference metric\n";
    // destruct all temporary objects
    delete tmpLogFld;
    if (ord == 1)
      delete tmpSolAV;
    else
      delete tmpSolAE;
    GM_delete(tmpMod);
    M_delete(tmpMsh);
  }

  printf("\n");
  printf("  -----------------------------\n");
  printf("  -         Adaptation        -\n");
  printf("  -----------------------------\n");
  printf("\n");

  mshAdp->run();

  std::cout << "\nINFO -- dof " << msh->nbrVertices() << "\n";

  M_writeMsh(msh, options["output"].c_str(), 2);
  std::cout << "INFO -- adapted mesh written to: " << options["output"] << "\n";
  
  // check mesh
  if ( !(mshAdp->checkTheMesh(1)) ) {	
    std::cerr << "WARNING -- unsuccessful mesh check " << "\n\n";
  }
  
  // compute mesh anisotropy statistics
  msh->reorderPoints();
  double maxAni = 1., meanAni = 0.;
  VIter iteVer = M_vertexIter(msh);
  while ( pVertex ver = iteVer->next() ) {
    doubleMatrix eigVec(2,2);
    doubleVector eigVal(2);
    MAdMetric met = logFld->findSize(ver)->getMetric();
    met.eig(eigVec, eigVal);
    if ( std::max( eigVal(0) / eigVal(1), eigVal(1) / eigVal(0) ) > maxAni ) {
      maxAni = std::max( eigVal(0) / eigVal(1), eigVal(1) / eigVal(0) );
    }
    meanAni += sqrt(std::max( eigVal(0) / eigVal(1), eigVal(1) / eigVal(0) ));
  }
  maxAni   = sqrt(maxAni);
  meanAni /= msh->nbrVertices();
  std::cout << "INFO -- maximum anisotropy: " << maxAni << "\n";
  std::cout << "INFO -- average anisotropy: " << meanAni << "\n";
  
  // compute mesh quality
  double meanShape, bestShape, worstShape;
  logFld->qualityInfo(meanShape, bestShape, worstShape);
  if ( meanShape < 0.75 || bestShape < 0.9 || worstShape < 1e-3 ) {
    std::cout << "INFO -- meanShape: " << meanShape << ", bestShape: " << bestShape;
    std::cout << ", worstShape: " << worstShape << "\n\n";
  }

  // compute min edge length
  double minEdgeLength = computeMinEdge(msh, logFld);
  std::cout << "\nINFO -- minimum edge length: " << minEdgeLength << "\n";
  
  VIter_delete(iteVer);
  delete logFld;
  if (ord == 1)
    delete solAV;
  else
    delete solAE;
  delete mshAdp;
  GM_delete(mod);
  M_delete(msh);

  return EXIT_SUCCESS;
}
