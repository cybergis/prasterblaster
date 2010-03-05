#ifndef SUBINFOVECT_H
#define SUBINFOVECT_H

/***************************************************************************
* subInfoVect.h
*
* Project: pRPL, v 1.0
* Purpose: Header file for class pRPL::InfoWLPair, and pRPL::SubInfoVect
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
#include "subSpaceInfo.h"

namespace pRPL {
  class InfoWLPair: public pair<SubSpaceInfo *, int> {
    public:
    InfoWLPair() {}
    InfoWLPair(SubSpaceInfo* pInfo, int wl)
      :pair<SubSpaceInfo*, int>(pInfo, wl) {}
    InfoWLPair(const InfoWLPair &rhs)
      :pair<SubSpaceInfo*, int>(rhs.first, rhs.second) {}
    ~InfoWLPair() {}

    InfoWLPair &operator=(const InfoWLPair &rhs) {
      if(this != &rhs) {
        this->first = rhs.first;
        this->second = rhs.second;
      }
      return *this;
    }
    bool operator==(const InfoWLPair &rhs) const {
      return (this->second == rhs.second);
    }
    bool operator!=(const InfoWLPair &rhs) const {
      return (this->second != rhs.second);
    }
    bool operator>(const InfoWLPair &rhs) const {
      return (this->second > rhs.second);
    }
    bool operator>=(const InfoWLPair &rhs) const {
      return (this->second >= rhs.second);
    }
    bool operator<(const InfoWLPair &rhs) const {
      return (this->second < rhs.second);
    }
    bool operator<=(const InfoWLPair &rhs) const {
      return (this->second <= rhs.second);
    }
  };

  class SubInfoVect: public vector<InfoWLPair> {
    public:
    SubInfoVect() {}
    SubInfoVect(const SubInfoVect &rhs)
      :vector<InfoWLPair>(rhs) {}
    ~SubInfoVect() {
      clear();
    }

    void clear();
    SubInfoVect::iterator erase(SubInfoVect::iterator iVal);
    SubInfoVect::iterator erase(SubInfoVect::iterator iBegin,
                                SubInfoVect::iterator iEnd);
    const SubSpaceInfo* findSubInfo(int subID) const;
  };
};

#endif
