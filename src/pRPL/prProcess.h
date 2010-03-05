#ifndef PRPROCESS_H
#define PRPROCESS_H

/***************************************************************************
* prProcess.h
*
* Project: pRPL, v 1.0
* Purpose: Header file for class pRPL::PRProcess
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
#include "subInfoVect.h"
#include "subSpaceMap.h"
#include "subIDMap.h"
#include "mpi.h"

namespace pRPL {
  class PRProcess {
    public:
      PRProcess();
      PRProcess(MPI_Comm comm,
                int groupID = -1);
      ~PRProcess();

      bool initialized() const;
      bool active() const;
      bool init(int argc = 0, 
                char *argv[] = 0);
      void abort() const;
      void finalize() const;
      void sync() const;

      bool set(MPI_Comm &comm,
               int groupID = -1);
      bool grouping(int nGroups,
                    bool incldMaster,
                    PRProcess *pGrpedPrc,
                    PRProcess *pGrpMaster = 0) const;

      const MPI_Comm &comm() const {
        return _comm;
      }
      int id() const {
        return _id;
      }
      int groupID() const {
        return _grpID;
      }
      int nPrcs() const {
        return _nTotalPrcs; 
      }
      int masterID() const {
        return _masterID;
      }
      bool isMaster() const;

      void cleanPrcInfoMap();
      void cleanSubIDMap();
      void cleanSubInfos();

      bool hasPrcInfoMap() const;
      bool hasSubIDMap() const;

      bool mapping();
      bool bcastSubInfos();

      void bcastStr(string &str) const;

      template <class elemType>
      void bcastVal(elemType &val) const;

      template <class elemType>
      void bcastVct(vector<elemType> &vect) const;

      template <class elemType>
      void gatherVal(elemType &val,
                     vector<elemType> &vect) const;

      template <class elemType>
      void gatherVct(vector<elemType> &vSend,
                     vector<elemType> &vRecv) const;

      template <class elemType>
      void allGatherVal(elemType &val,
                        vector<elemType> &vect) const;

      template <class elemType>
      void allGatherVct(vector<elemType> &vSend,
                        vector<elemType> &vRecv) const;

      SubInfoVect &glbSubInfos() {
        return _vpGlbSubSpcInfos;
      }
      const SubInfoVect &glbSubInfos() const {
        return _vpGlbSubSpcInfos;
      }
      int nTotSubSpcs() const {
        return _vpGlbSubSpcInfos.size();
      }

      SubSpaceMap *prcInfoMap() {
        return _pPrcInfoMap;
      }
      const SubSpaceMap *prcInfoMap() const {
        return _pPrcInfoMap;
      }

      SubIDMap *subIDMap() {
        return _pSubIDMap;
      }
      const SubIDMap *subIDMap() const {
        return _pSubIDMap;
      }

      const vector<SubSpaceInfo *> &lclSubInfos() const {
        return (*_pPrcInfoMap)[_id];
      }
      int nLclSubSpcs() const {
        return (*_pPrcInfoMap)[_id].size();
      }

    private:
      MPI_Comm _comm;
      int _id;
      int _grpID;
      int _nTotalPrcs;
      int _masterID;

      SubInfoVect _vpGlbSubSpcInfos;
      SubSpaceMap *_pPrcInfoMap;

      SubIDMap *_pSubIDMap;
  };
};

inline pRPL::PRProcess::
PRProcess()
  :_comm(MPI_COMM_NULL),
   _id(-1),
   _grpID(-1),
   _nTotalPrcs(-1),
   _masterID(-1),
   _pPrcInfoMap(0),
   _pSubIDMap(0) {}

inline pRPL::PRProcess::
PRProcess(MPI_Comm comm,
          int groupID)
  :_comm(comm),
   _id(-1),
   _grpID(groupID),
   _nTotalPrcs(0),
   _masterID(0),
   _pPrcInfoMap(0),
   _pSubIDMap(0) {}

inline pRPL::PRProcess::
~PRProcess() {
  cleanPrcInfoMap();
  cleanSubIDMap();
  cleanSubInfos();
}

inline void pRPL::PRProcess::
abort() const {
  MPI_Abort(_comm, -1);
}

inline void pRPL::PRProcess::
finalize() const {
  MPI_Finalize();
}

inline void pRPL::PRProcess::
sync() const {
  MPI_Barrier(_comm);
}

inline bool pRPL::PRProcess::
isMaster() const {
  return (_id == _masterID);
}

template <class elemType>
void pRPL::PRProcess::
bcastVal(elemType &val) const {
  if(active()) {
    MPI_Bcast(&val, sizeof(elemType), MPI_CHAR, _masterID, _comm);
  }
}

template <class elemType>
void pRPL::PRProcess::
bcastVct(vector<elemType> &vect) const {
  if(!active()) {
    return;
  }
  int vectSize = 0;
  if(isMaster()) {
    vectSize = vect.size();
  }
  else {
    vect.clear();
  }
  MPI_Bcast(&vectSize, 1, MPI_INT, _masterID, _comm);

  if(!isMaster()) {
    vect.resize(vectSize);
  }

  MPI_Bcast(&(vect[0]), vectSize*sizeof(elemType), MPI_CHAR, _masterID, _comm);
}

template <class elemType>
void pRPL::PRProcess::
gatherVal(elemType &val,
          vector<elemType> &vect) const {
  if(!active()) {
    return;
  }
  if(isMaster()) {
    vect.clear();
    vect.resize(_nTotalPrcs);
  }
  MPI_Gather(&val, sizeof(elemType), MPI_CHAR,
             &(vect[0]), sizeof(elemType), MPI_CHAR,
             _masterID, _comm);
}
                           
template <class elemType>
void pRPL::PRProcess::
gatherVct(vector<elemType> &vSend,
          vector<elemType> &vRecv) const {
  if(!active()) {
    return;
  }
  int size = 0;
  int *aSizes = 0;
  int *aDspls = 0;

  if(isMaster()) {
    aSizes = (int *)malloc(_nTotalPrcs * sizeof(int));
    aDspls = (int *)malloc(_nTotalPrcs * sizeof(int));
  }
  size = vSend.size() * sizeof(elemType);
  
  MPI_Gather(&size, 1, MPI_INT,
             aSizes, 1, MPI_INT,
             _masterID, _comm);
  
  if(isMaster()) {
    int totSize = 0;
    for(int i = 0; i < _nTotalPrcs; i++) {
      aDspls[i] = totSize;
      totSize += aSizes[i];
    }
    vRecv.resize(totSize / sizeof(elemType));
  }
  
  MPI_Gatherv(&(vSend[0]), size, MPI_CHAR,
              &(vRecv[0]), aSizes, aDspls, MPI_CHAR,
              _masterID, _comm);
  
  if(isMaster()) {
    free(aSizes);
    free(aDspls);
  }
}

template <class elemType>
void pRPL::PRProcess::
allGatherVal(elemType &val,
             vector<elemType> &vect) const {
  if(!active()) {
    return;
  }
  vect.clear();
  vect.resize(_nTotalPrcs);
  MPI_Allgather(&val, sizeof(elemType), MPI_CHAR,
                &(vect[0]), sizeof(elemType), MPI_CHAR,
                _comm);
}
                        
template <class elemType>
void pRPL::PRProcess::
allGatherVct(vector<elemType> &vSend,
             vector<elemType> &vRecv) const {
  if(!active()) {
    return;
  }
  int *aSizes = (int *)malloc(_nTotalPrcs * sizeof(int));
  int *aDspls = (int *)malloc(_nTotalPrcs * sizeof(int));
  int size = vSend.size() * sizeof(elemType);

  MPI_Allgather(&size, 1, MPI_INT,
                aSizes, 1, MPI_INT, _comm);
  
  int totSize = 0;
  for(int i = 0; i < _nTotalPrcs; i++) {
    aDspls[i] = totSize;
    totSize += aSizes[i];
  }
  vRecv.resize(totSize / sizeof(elemType));
  
  MPI_Allgatherv(&(vSend[0]), size, MPI_CHAR,
                 &(vRecv[0]), aSizes, aDspls, MPI_CHAR,
                 _comm);
  
  free(aSizes);
  free(aDspls);
}

#endif
