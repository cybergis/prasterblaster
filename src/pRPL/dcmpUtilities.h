#ifndef DCMPUTILITIES_H
#define DCMPUTILITIES_H

/***************************************************************************
* dcmpUtilities.h
*
* Project: pRPL, v 1.0
* Purpose: Header file for utility funcitons for decomposition
* Author:  Qingfeng (Gene) Guan
* E-mail:  guanqf {at} gmail.com
****************************************************************************
* Copyright (c) 2008, Qingfeng Guan
* NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
* purposes, NO COMMERCIAL usages are allowed unless the author is 
* contacted and a permission is granted
* 
****************************************************************************/

#include "basicTypes.h"
#include "cellSpace.h"
#include "neighborhood.h"

namespace pRPL {
  static const int nbrLocations[16] = {-1, 0,
                                       -1, 1,
                                       0, 1,
                                       1, 1,
                                       1, 0,
                                       1, -1,
                                       0, -1,
                                       -1, -1};

  template <class elemType>
  bool checkDcmpParms(const SpaceDims &cellSpaceDims,
                      const Neighborhood<elemType> &nbrhood);
};

template <class elemType>
inline bool pRPL::
checkDcmpParms(const SpaceDims &cellSpaceDims,
               const Neighborhood<elemType> &nbrhood) {
  bool parmsValid = true;
  if(!cellSpaceDims.valid()) {
    cerr << __FILE__ << " " << __FUNCTION__ 
         << " Error: unable to decompose an empty CellSpace" << endl;
    parmsValid = false;
  }
  else if(nbrhood.empty()) {
    cerr << __FILE__ << " " << __FUNCTION__ 
         << " Error: unable to decompose CellSpace using" \
         << " an empty Neighborhood" << endl;
    parmsValid = false;
  }
  else if(cellSpaceDims.nRows() < nbrhood.nRows() ||
          cellSpaceDims.nCols() < nbrhood.nCols()) {
    cerr << __FILE__ << " " << __FUNCTION__ 
         << " Error: the CellSpace's dimensions have to be greater" \
         << " than that of the Neighborhood" << endl;
    parmsValid = false;
  }
  return parmsValid;
}

#endif
