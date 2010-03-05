#include "subIDMap.h"

/***************************************************************************
* subIDMap.cpp
*
* Project: pRPL, v 1.0
* Purpose: Implementation for class pRPL::SubIDMap
* Author:  Qingfeng (Gene) Guan
* E-mail:  guanqf {at} gmail.com
****************************************************************************
* Copyright (c) 2008, Qingfeng Guan
* NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
* purposes, NO COMMERCIAL usages are allowed unless the author is 
* contacted and a permission is granted
* 
****************************************************************************/

void pRPL::SubIDMap::
toIntVect(IntVect &vIDPack) const {
  vIDPack.erase(vIDPack.begin(), vIDPack.end());
  vIDPack.push_back(size());
  SubIDMap::const_iterator iVal = begin();
  while(iVal != end()) {
    vIDPack.push_back((*iVal).first);
    vIDPack.push_back((*iVal).second);
    iVal++;
  }
}

bool pRPL::SubIDMap::
fromIntVect(const IntVect &vIDPack) {
  if(vIDPack.empty()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to load an empty vector to a SubIDMap" \
         << endl;
    return false;
  }
  int nVals = vIDPack[0];
  if(nVals <= 0) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid size of SubIDMap (" \
         << nVals << "). It should be greater than ZERO"
         << endl;
    return false;
  }
  if(vIDPack.size() != nVals*2 + 1) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid size of the vector (" \
         << vIDPack.size() << "). It should equal to" \
         << nVals*2 + 1 << endl;
    return false;
  }
  for(int iVal = 1; iVal < vIDPack.size()-1; iVal+=2) {
    int subSpcID = vIDPack[iVal];
    int prcID = vIDPack[iVal+1];
    if(subSpcID < 0 || prcID < 0) {
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: invalid SubSpcID (" \
           << subSpcID << "), or invalid PrcID (" \
           << prcID << ")" \
           << endl;
      return false;
    }
    insert(SubIDMap::value_type(subSpcID, prcID));
  }
  return true;
}

int pRPL::SubIDMap::
maxID() const {
  map<int, int>::const_iterator iMax
    = std::max_element(begin(), end());
  return (*iMax).first;
}

int pRPL::SubIDMap::
minID() const {
  map<int, int>::const_iterator iMin
    = std::min_element(begin(), end());
  return (*iMin).first;
}
