/***************************************************************************
* pAspect.cpp
*
* Project: pRPL, v 1.0
* Purpose: Demonstration program for pRPL. It derives aspects from DEM data
*
* Usage:   pAspect input-demFilename output-aspFilename cellsize [scale] dcmpMethod
*
*          cellsize -   the length of the cell's border, assuming the cell 
*                       is a square. The cellsize for the test data, i.e., usa_nw.dem,
*                       is 0.008333 in geographic latitude-longitude degree
*          scale -      the horizontal-to-vertical unit transformation scale. 
*                       The scale for the test data is 111120 (latitude-longitude 
*                       degrees to meters, 1 degree roughly equals to 111120 meters
*                       on the equator)
*          dcmpMethod - the method for decomposing the data. 
*                       1: ROW-WISE
*                       2: COLUMN-WISE
*                       3: BLOCK-WISE
*
* Example: mpirun -np 4 ./pAspect ./data/usa_nw.dem ./data/usa_nw.asp 0.008333 111120 1
*          Use 4 processors, and ROW-WISE decomposition
*
* Author:  Qingfeng (Gene) Guan
* E-mail:  guanqf {at} gmail.com
****************************************************************************
* Copyright (c) 2008, Qingfeng Guan
* 
****************************************************************************/

#include "layer.h"
#include "aspectTransition.h"
#include <iostream>
#include <string>
#include <fstream>
#include <time.h>
#include <sstream>

using namespace std;
using namespace pRPL;

int main(int argc, char *argv[]) {
  const string usage("usage: pAspect input-demFilename output-aspFilename cellsize [scale] dcmpMethod");
  
  /* Declare processes, and initialize them */
  PRProcess aspPrc(MPI_COMM_WORLD);
  aspPrc.init(argc, argv);

  /* Handle arguments */  
  if(argc != 5 && argc != 6) {
    if(0 == aspPrc.id()) {
      cout << usage << endl;
    }
    aspPrc.finalize();
    return -1;
  }
  string demFilename, aspectFilename;
  demFilename.assign(argv[1]);
  aspectFilename.assign(argv[2]);
  double cellSize = atof(argv[3]);
  double scale;
  DomDcmpMethod dcmpMethod;
  int nSubSpcs = aspPrc.nPrcs();
  if(argc == 5) {
    scale = 1.0;
    dcmpMethod = static_cast<DomDcmpMethod>(atoi(argv[4]));
  }
  else {
    scale = atof(argv[4]);
    dcmpMethod = static_cast<DomDcmpMethod>(atoi(argv[5]));
  }
  
  /* Declare layer and transition objects */
  Layer<double> aspectLayer(aspPrc, "aspect");
  Layer<short> demLayer(aspPrc, "dem");
  AspectTransition calcAspect(cellSize, scale);

  /* Load dem data and neighborhood configuration */
  if(demLayer.isMaster()) {
    demLayer.newCellSpace();
    fstream demFile(demFilename.c_str(), ios::in);
    demFile >> *(demLayer.cellSpace());
    demFile.close();
    SpaceDims dims = (demLayer.cellSpace())->dims();
    aspectLayer.newCellSpace(dims);

    demLayer.newNbrhood();
    fstream nbrFile("./nbrhoods/moore.nbr", ios::in);
    nbrFile >> *(demLayer.nbrhood());
    nbrFile.close();
    
    aspectLayer.newNbrhood(*(demLayer.nbrhood()));
  }
  
  /* Record the start time */
  double timeStart;
  aspPrc.sync();
  if(aspPrc.isMaster()) {
    timeStart = MPI_Wtime(); 
  }
  
  /* Decompose and distribute layers */
  if(!demLayer.smplDcmpDstrbt(dcmpMethod, nSubSpcs, nSubSpcs)) {
    aspPrc.abort();
    return -1;
  }
  cout << aspPrc.id() << " distributing dem ok" << endl;
  
  if(!aspectLayer.distribute(false)) {
    cout << "layer[" << aspectLayer.title() 
         << "] error during distributing" << endl;
    aspPrc.abort();
  }
  cout << aspPrc.id() << " distributing aspect ok" << endl;
  
  /* Calculate aspects */
  cout << aspPrc.id() << "updating aspect layer" << endl;
  calcAspect.demLayer(demLayer);
  if(!aspectLayer.update(calcAspect)) {
    cout << "layer[" << aspectLayer.title() 
         << "] error during update" << endl;
    aspPrc.abort();
    return -1;
  }
  
  /* Gather aspects */
  cout << aspPrc.id() << "gathering aspect layer" << endl;
  if(!aspectLayer.gatherCellSpace()) {
    cout << "layer[" << aspectLayer.title() 
         << "] error during gathering cellspace" << endl;
    aspPrc.abort();
    return -1;
  }
  
  /* Record the end time, log computing time */
  double timeEnd;
  aspPrc.sync();
  if(aspPrc.isMaster()) {
    timeEnd = MPI_Wtime();

    ofstream timeFile;
    timeFile.open("./data/time.log", ios::app);
    if(!timeFile) {
      cerr << "Error: unable to open the time log file" << endl;
    }
    timeFile << demFilename << "\t" \
             << aspectFilename << "\t" \
             << dcmpMethod << "\t" \
             << aspPrc.nPrcs() << "\t" \
             << timeEnd - timeStart \
             << endl;
    timeFile.close();
  }
  
  /* Output aspects */
  if(aspectLayer.isMaster()) {
    fstream slopeFile((aspectFilename).c_str(), ios::out);
    slopeFile << *(aspectLayer.cellSpace());
    slopeFile.close();
  }
  
  /* Release gathering types and finalize processes */
  aspectLayer.freeGatherTypes();
  aspPrc.finalize();
  
  return 0;
}
