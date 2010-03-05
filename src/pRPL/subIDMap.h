#ifndef SUBIDMAP_H
#define SUBIDMAP_H

/***************************************************************************
* subIDMap.h
*
* Project: pRPL, v 1.0
* Purpose: Header file for class pRPL::SubIDMap
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

namespace pRPL {
  class SubIDMap: public map<int, int> {
    public:
      SubIDMap() 
        :map<int, int>() {}
      SubIDMap(const SubIDMap &rhs)
        :map<int, int>(rhs) {}
      SubIDMap(const IntVect &vIDPack);
      ~SubIDMap() {}

      void toIntVect(IntVect &vIDPack) const;
      bool fromIntVect(const IntVect &vIDPack);

      int maxID() const;
      int minID() const;
  };
};

inline pRPL::SubIDMap::
SubIDMap(const IntVect &vIDPack)
  :map<int, int>() {
  if(!fromIntVect(vIDPack)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to construct a SubIDMap" \
         << endl;
    exit(-1);
  }
}

#endif
