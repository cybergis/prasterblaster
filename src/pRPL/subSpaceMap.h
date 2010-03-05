#ifndef SUBSPACEMAP_H
#define SUBSPACEMAP_H

/***************************************************************************
* subSpaceMap.h
*
* Project: pRPL, v 1.0
* Purpose: Header file for class pRPL::SubSpaceMap
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
#include "subIDMap.h"
#include "subInfoVect.h"

namespace pRPL {
  typedef map<int, vector<SubSpaceInfo *> > PrcInfoMap;
  typedef PrcInfoMap::value_type PIMapValue;
  typedef PrcInfoMap::iterator PIMapItr;

  class SubSpaceMap {
    public:
      SubSpaceMap()
        :_pvpSubSpcInfos(0),
         _nPrcs(0),
         _masterID(-1),
         _mapped(false) {}
      SubSpaceMap(const SubInfoVect &vpSubSpcInfos,
                  int nPrcs, int masterID);
      SubSpaceMap(const SubSpaceMap &rhs)
        :_pvpSubSpcInfos(rhs._pvpSubSpcInfos),
         _nPrcs(rhs._nPrcs),
         _masterID(rhs._masterID),
         _mPrcInfos(rhs._mPrcInfos),
         _mapped(rhs._mapped) {}

      ~SubSpaceMap() {}

      const vector<SubSpaceInfo *>& operator[](int iPrc) {
        return _mPrcInfos[iPrc];
      }
      
      bool mapping(SubIDMap &subIDMap);
      bool mapped() const {
        return _mapped;
      }
      int nPrcs() const {
        return _nPrcs;
      }
      int masterID() const {
        return _masterID;
      }
      const map<int, int>& workloadMap() const {
        return _mWLs;
      }
      
    private:
      void _changePrcID(int &iPrc, bool &prcUp,
                        int stepSize);

    private:
      const SubInfoVect *_pvpSubSpcInfos;
      int _nPrcs;
      int _masterID;
      PrcInfoMap _mPrcInfos;
      map<int, int> _mWLs;
      bool _mapped;
  };
};

inline pRPL::SubSpaceMap::
SubSpaceMap(const SubInfoVect &vpSubSpcInfos,
            int nPrcs, int masterID)
  :_mapped(false) {
  if(vpSubSpcInfos.empty()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to construct a SubSpaceMap" \
         << " with an empty vector of SubSpaceInfo" \
         << endl;
    exit(-1);
  }
  _pvpSubSpcInfos = &vpSubSpcInfos;

  /*
  if(nPrcs <= 1) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to construct a SubSpaceMap." \
         << " The number of processes should be greater than ONE" \
         << endl;
    exit(-1);
  }
  */
  _nPrcs = nPrcs;
  
  if(masterID < 0 ||
     masterID >= nPrcs) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to construct a SubSpaceMap." \
         << " Invalid master node ID (" << masterID<< ")" \
         << endl;
    exit(-1);
  }
  _masterID = masterID;
}

#endif
