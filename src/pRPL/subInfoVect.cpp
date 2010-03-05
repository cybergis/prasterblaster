#include "subInfoVect.h"

/***************************************************************************
* subInfoVect.cpp
*
* Project: pRPL, v 1.0
* Purpose: Implementation for class pRPL::SubInfoVect
* Author:  Qingfeng (Gene) Guan
* E-mail:  guanqf {at} gmail.com
****************************************************************************
* Copyright (c) 2008, Qingfeng Guan
* NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
* purposes, NO COMMERCIAL usages are allowed unless the author is 
* contacted and a permission is granted
* 
****************************************************************************/

void pRPL::SubInfoVect::
clear() {
  SubInfoVect::iterator iVal = begin();
  while(iVal != end()) {
    if((*iVal).first) {
      delete (*iVal).first;
      (*iVal).first = 0;
    }
    iVal++;
  }
  vector<InfoWLPair>::clear();
}

pRPL::SubInfoVect::iterator pRPL::SubInfoVect::
erase(SubInfoVect::iterator iVal) {
  if((*iVal).first) {
    delete (*iVal).first;
    (*iVal).first = 0;
  }
  return vector<InfoWLPair>::erase(iVal);
}

pRPL::SubInfoVect::iterator pRPL::SubInfoVect::
erase(SubInfoVect::iterator iBegin,
      SubInfoVect::iterator iEnd) {
  SubInfoVect::iterator iVal = iBegin;
  while(iVal != iEnd) {
    if((*iVal).first) {
      delete (*iVal).first;
      (*iVal).first = 0;
    }
    iVal++;
  }
  return vector<InfoWLPair>::erase(iBegin, iEnd);
}

const pRPL::SubSpaceInfo* pRPL::SubInfoVect::
findSubInfo(int subID) const {
  const SubSpaceInfo* pSubInfo = 0;
  SubInfoVect::const_iterator iSubInfo = begin();
  while(iSubInfo != end()) {
    if((*iSubInfo).first->id() == subID) {
      pSubInfo = (*iSubInfo).first;
      break;
    }
    iSubInfo++;
  }
  return pSubInfo;
}
