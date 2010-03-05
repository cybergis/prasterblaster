#include "subSpaceMap.h"

/***************************************************************************
* subSpaceMap.cpp
*
* Project: pRPL, v 1.0
* Purpose: Implementation for class pRPL::SubSpaceMap
* Author:  Qingfeng (Gene) Guan
* E-mail:  guanqf {at} gmail.com
****************************************************************************
* Copyright (c) 2008, Qingfeng Guan
* NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
* purposes, NO COMMERCIAL usages are allowed unless the author is 
* contacted and a permission is granted
* 
****************************************************************************/

void pRPL::SubSpaceMap::
_changePrcID(int &iPrc, bool &prcUp,
             int stepSize) {
  if(prcUp) {
    iPrc++;
    if(iPrc == _nPrcs) {
      prcUp = false;
      iPrc -= stepSize;
    }
  }
  else {
    iPrc--;
    if(iPrc == -1) {
      prcUp = true;
      iPrc += stepSize;
    }
  }
}

bool pRPL::SubSpaceMap::
mapping(SubIDMap &subIDMap) {
  if(_mapped) {
    return false;
  }
  /*
  subIDMap.erase(subIDMap.begin(), subIDMap.end());
  _mPrcInfos.erase(_mPrcInfos.begin(), _mPrcInfos.end());
  */
  _mWLs.clear();
  
  int iPrc;
  for(iPrc = 0; iPrc < _nPrcs; iPrc++) {
    _mPrcInfos.insert(PIMapValue(iPrc, vector<SubSpaceInfo *>()));
    _mWLs[iPrc] = 0;
  }
  
  iPrc = 0;
  bool prcUp = true;
  for(int iInfo = 0; iInfo < _pvpSubSpcInfos->size(); iInfo++) {
    /*
    if(iPrc == _masterID) {
      _changePrcID(iPrc, prcUp, 2);
    }
    */
    SubSpaceInfo *pInfo = (_pvpSubSpcInfos->at(iInfo)).first;
    _mPrcInfos[iPrc].push_back(pInfo);
    _mWLs[iPrc] += (_pvpSubSpcInfos->at(iInfo)).second;
    subIDMap.insert(SubIDMap::value_type(pInfo->id(), iPrc));
    _changePrcID(iPrc, prcUp, 1);
  }
  _mapped = true;
  return true;
}
