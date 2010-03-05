#include "prProcess.h"

/***************************************************************************
* prProcess.cpp
*
* Project: pRPL, v 1.0
* Purpose: Implementation for class pRPL::PRProcess
* Author:  Qingfeng (Gene) Guan
* E-mail:  guanqf {at} gmail.com
****************************************************************************
* Copyright (c) 2008, Qingfeng Guan
* NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
* purposes, NO COMMERCIAL usages are allowed unless the author is 
* contacted and a permission is granted
* 
****************************************************************************/

bool pRPL::PRProcess::
initialized() const{
  int mpiStarted;
  MPI_Initialized(&mpiStarted);
  return static_cast<bool>(mpiStarted);
}

bool pRPL::PRProcess::
active() const {
  return (_comm != MPI_COMM_NULL &&
          _id != -1 &&
          _nTotalPrcs != -1);
}

bool pRPL::PRProcess::
init(int argc,
     char *argv[]) {
  bool done = true;
  if(_comm == MPI_COMM_NULL) {
    return done;
  }
  
  if(!initialized()) {
    MPI_Init(&argc, &argv);
  }
  MPI_Comm_rank(_comm, &_id);
  MPI_Comm_size(_comm, &_nTotalPrcs);
  
  if(_masterID < 0 ||
     _masterID >= _nTotalPrcs) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid master ID (" \
         << _masterID << "). There are totally " \
         << _nTotalPrcs << " processes." \
         << endl;
    done = false;
  }
  return done;
}

bool pRPL::PRProcess::
set(MPI_Comm &comm,
    int groupID) {
  _comm = comm;
  _masterID = 0;
  _grpID = groupID;
  return(init());
}

bool pRPL::PRProcess::
grouping(int nGroups,
         bool incldMaster,
         PRProcess *pGrpedPrc,
         PRProcess *pGrpMaster) const {
  if(!initialized()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: PRProcess has NOT been initialized," \
         << " unable to be grouped" << endl;
    return false;
  }
  
  if(!active()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: inactive PRProcess," \
         << " unable to group a Null communicator." \
         << " id = " << _id << " nTotPrcs = " << _nTotalPrcs << endl;
    return false;
  }

  if(nGroups <= 0 ||
     nGroups > _nTotalPrcs) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid number of groups (" \
         << nGroups << ") as the total number of processes is " \
         << _nTotalPrcs << endl;
    return false;
  }
  
  if(!incldMaster && _nTotalPrcs <= 1) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error:  " << _nTotalPrcs << " processes can NOT" \
         << " be grouped without the master process" << endl;
    return false;
  }

  MPI_Group glbGrp;
  MPI_Comm glbComm = _comm;
  MPI_Comm_group(glbComm, &glbGrp);
  int myID = -1;
  int grpID = -1;
  MPI_Comm grpComm = MPI_COMM_NULL;
  
  if(incldMaster) {
    myID = _id;
    grpID = myID % nGroups; 
    MPI_Comm_split(glbComm, grpID, myID, &grpComm);
    if(!pGrpedPrc->set(grpComm, grpID)) {
      return false;
    }
    if(pGrpMaster != 0) {
      MPI_Group masterGrp= MPI_GROUP_NULL;
      MPI_Comm masterComm = MPI_COMM_NULL;
      int grpMasterRange[1][3] = {{0, nGroups-1, 1}};
      MPI_Group_range_incl(glbGrp, 1, grpMasterRange, &masterGrp);
      MPI_Comm_create(glbComm, masterGrp, &masterComm);
      if(!pGrpMaster->set(masterComm)) {
        return false;
      }
    }
  }
  else {
    int excldRanks[1] = {_masterID};
    MPI_Group glbGrp2 = MPI_GROUP_NULL;
    MPI_Group_excl(glbGrp, 1, excldRanks, &glbGrp2);
    MPI_Comm_create(_comm, glbGrp2, &glbComm);
    glbGrp = glbGrp2;
    if(!isMaster()) {
      MPI_Comm_rank(glbComm, &myID);
      grpID = myID % nGroups;
      MPI_Comm_split(glbComm, grpID, myID, &grpComm);
      if(!pGrpedPrc->set(grpComm, grpID)) {
        return false;
      }
      if(pGrpMaster != 0) {
        MPI_Group masterGrp= MPI_GROUP_NULL;
        MPI_Comm masterComm = MPI_COMM_NULL;
        int grpMasterRange[1][3] = {{0, nGroups-1, 1}};
        MPI_Group_range_incl(glbGrp, 1, grpMasterRange, &masterGrp);
        MPI_Comm_create(glbComm, masterGrp, &masterComm);
        if(!pGrpMaster->set(masterComm)) {
          return false;
        }
      }
    }
  }
  
  return true;
}

void pRPL::PRProcess::
cleanPrcInfoMap() {
  if(_pPrcInfoMap) {
    delete _pPrcInfoMap;
    _pPrcInfoMap = 0;
  }
}

void pRPL::PRProcess::
cleanSubIDMap() {
  if(_pSubIDMap) {
    delete _pSubIDMap;
    _pSubIDMap = 0;
  }
}

void pRPL::PRProcess::
cleanSubInfos() {
  _vpGlbSubSpcInfos.clear();
}

bool pRPL::PRProcess::
hasPrcInfoMap() const {
  bool hasIt = true;
    if(!_pPrcInfoMap ||
       !_pPrcInfoMap->mapped()) {
      hasIt = false;
    }
  return hasIt;
}

bool pRPL::PRProcess::
hasSubIDMap() const {
  bool hasIt = true;
  if(!_pSubIDMap) {
    hasIt = false;
  }
  return hasIt;
}

bool pRPL::PRProcess::
mapping() {
  if(!active()) {
    return true;
  }
  
  if(_vpGlbSubSpcInfos.empty()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: the vector of SubSpaceInfos on the" \
         << " master process[" << _id << "] is empty" \
         << endl;
    return false;
  }

  cleanPrcInfoMap();
  cleanSubIDMap();
  _pPrcInfoMap = new SubSpaceMap(_vpGlbSubSpcInfos, 
                                 _nTotalPrcs, _masterID);
  _pSubIDMap = new SubIDMap();
  if(!_pPrcInfoMap->mapping(*_pSubIDMap)) {
    return false;
  }
  /*
  if(isMaster()) {
    cout << *_pSubIDMap << endl;
  }
  */
  return true;
}

bool pRPL::PRProcess::
bcastSubInfos() {
  if(!active()) {
    return true;
  }
  
  IntVect vInfoPack;
  int packSize;
  if(isMaster()) { /* Master process */
    for(int iSubInfo = 0; iSubInfo < _vpGlbSubSpcInfos.size(); iSubInfo++) {
      (_vpGlbSubSpcInfos[iSubInfo].first)->add2IntVect(vInfoPack);
      vInfoPack.push_back(_vpGlbSubSpcInfos[iSubInfo].second);
    }
    packSize = vInfoPack.size();
  }
  else {
    cleanSubInfos();
  }

  MPI_Bcast(&packSize, 1, MPI_INT, _masterID, _comm);
  if(!isMaster()) {
    vInfoPack.resize(packSize);
  }

  MPI_Bcast(&(vInfoPack[0]), packSize, MPI_INT, _masterID, _comm);

  if(!isMaster()) {
    IntVctItr iVal = vInfoPack.begin();
    while(iVal != vInfoPack.end()) {
      _vpGlbSubSpcInfos.push_back(InfoWLPair(new SubSpaceInfo(vInfoPack, iVal),
                                             *iVal));
      iVal++;
    }
  }

  return true;
}

void pRPL::PRProcess::
bcastStr(string &str) const {
  if(!active()) {
    return;
  }
  
  int strSize = 0;
  if(isMaster()) {
    strSize = str.size();
  }
  MPI_Bcast(&strSize, 1, MPI_INT, _masterID, _comm);
  if(!isMaster()) {
    str.resize(strSize);
  }
  MPI_Bcast(&(str[0]), strSize, MPI_CHAR, _masterID, _comm);
}
