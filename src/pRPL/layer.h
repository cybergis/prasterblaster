#ifndef LAYER_H
#define LAYER_H

/***************************************************************************
* layer.h
*
* Project: pRPL, v 1.0
* Purpose: Header file for class pRPL::Layer
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
#include "prProcess.h"
#include "cellSpace.h"
#include "subSpace.h"
#include "subSpaceMap.h"
#include "subIDMap.h"
#include "neighborhood.h"
#include "smplDcmp.h"
#include "quadtreeDcmp.h"
#include "subInfoVect.h"
#include "mpi.h"

namespace pRPL{
  template <class elemType>
  class Layer {
    public:
      Layer();
      Layer(PRProcess &prPrc,
            const string layerName = "Untitled");
      Layer(const Layer<elemType> &rhs);
      ~Layer();

      Layer<elemType>& operator=(const Layer<elemType> &rhs);

      const string &name() const;
      void name(const string &layerName);
      bool hasPRPrc() const;
      PRProcess *const prPrc() const;
      int id() const;
      const string title() const;
      bool isMaster() const;

      void cleanCellSpace();
      void cleanNbrhood();
      void cleanPrcInfoMap();
      void cleanSubIDMap();
      void cleanSubSpcs();

      bool hasCellSpace() const;
      bool hasNbrhood() const;

      void prPrc(PRProcess &prPrc);

      bool newCellSpace();
      bool newCellSpace(const SpaceDims& dims);
      bool newCellSpace(const SpaceDims& dims,
                        const elemType &initVal);
      bool newCellSpace(int nRows, int nCols);
      bool newCellSpace(int nRows, int nCols,
                        const elemType &initVal);

      bool newNbrhood();
      bool newNbrhood(const vector<CellCoord> &vNbrCoords,
                      double weight = 1.0);
      bool newNbrhood(const vector<CellCoord> &vNbrCoords,
                      const vector<double> &vNbrWeights);
                      
      template<class elemType2>
      bool newNbrhood(const Neighborhood<elemType2> &nbr);

      bool smplDcmp(DomDcmpMethod dcmpMethod,
                    int nSubSpcs1,
                    int nSubSpcs2 = 1);
      bool quadDcmp(Transition<elemType> &transition,
                    int maxNumLeaves,
                    int minWorkload = 0);

      bool bcastNbrhood();
      bool dstrbtCellSpace(bool dstrbtData = true);
      bool bcastCellSpace(bool bcastData = true);
      bool distribute(bool dstrbtData = true);
      bool broadcast(bool bcastData = true);

      bool smplDcmpDstrbt(DomDcmpMethod dcmpMethod,
                          int nSubSpcs1,
                          int nSubSpcs2 = 1,
                          bool dstrbtData = true);
      bool quadDcmpDstrbt(Transition<elemType> &transition,
                          int maxNumLeaves,
                          int minWorkload = 0,
                          bool dstrbtData = true);

      CellSpace<elemType> *cellSpace();
      const CellSpace<elemType> *cellSpace() const;
      Neighborhood<elemType> *nbrhood();
      const Neighborhood<elemType> *nbrhood() const;
      SubInfoVect *glbSubInfos();

      SubIDMap *subIDMap();

      int nLclSubSpcs() const;
      SubSpace<elemType> *operator[](int iSubSpc);
      const SubSpace<elemType> *operator[](int iSubSpc) const;
      SubSpace<elemType> *findSubSpc(int subID);
      const SubSpace<elemType> *findSubSpc(int subID) const;
      double sizeRatio() const;

      void resetExchange();
      bool exchangeBegin(SubSpace<elemType> *pSubSpc);
      bool exchangeSize();
      bool exchangeEnd(Transition<elemType> *pTransition);

      bool update(Transition<elemType> &transition);
      bool update(map<int, IntVect> &mLclIdxs2Eval,
                  Transition<elemType> &transition);
      bool update(IntVect &vGlbIdxs2Eval,
                  Transition<elemType> &transition);

      bool buildGatherTypes();
      void freeGatherTypes();
      bool gatherCellSpace();

    protected:
      PRProcess *_pPRPrc;
      string _name;

      CellSpace<elemType> *_pCellSpace; // Master process
      Neighborhood<elemType> *_pNbrhood; // Master & Slave
      map<int, MPI_Datatype> _mGatherTypes; // Master
      vector<SubSpace<elemType> *> _vpSubSpcs; // Slave process

      vector<MPI_Request> _vInfoReqs; // Slave process
      vector<MPI_Request> _vUpdtReqs; // Slave process
      int _iUpdtInfo; // Slave process
      int _iUpdtPack; // Slave process
      map<int, char* > _mSendLocs; // Slave process
      map<int, IntPair> _mSendSizes; // Slave process
      map<int, char* > _mRecvLocs; // Slave process
      map<int, pair<SubSpace<elemType> *, IntPair> > _mRecvSizes; // Slave process
  };
};

template <class elemType>
inline pRPL::Layer<elemType>::
Layer()
  :_pPRPrc(0),
   _name("Untitled"),
   _pCellSpace(0),
   _pNbrhood(0){}

template <class elemType>
inline pRPL::Layer<elemType>::
Layer(PRProcess &prPrc,
      const string layerName)
  :_name(layerName),
   _pCellSpace(0),
   _pNbrhood(0) {
  if(!prPrc.initialized()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: PRProcess has NOT been initialized yet." \
         << " Unable to construct a process." << endl;
    exit(-1);
  }
  _pPRPrc = &prPrc;
}

template <class elemType>
inline pRPL::Layer<elemType>::
Layer(const Layer<elemType> &rhs)
  :_pPRPrc(rhs._pPRPrc),
   _name(rhs._name),
   _pCellSpace(0),
   _pNbrhood(0),
   _mGatherTypes(rhs._mGatherTypes){
  if(rhs._pCellSpace) {
    _pCellSpace = new CellSpace<elemType>(*(rhs._pCellSpace));
  }
  if(rhs._pNbrhood) {
    _pNbrhood = new Neighborhood<elemType>(*(rhs._pNbrhood));
  }
  for(int iSubSpc = 0; iSubSpc < rhs._vpSubSpcs.size(); iSubSpc++) {
    _vpSubSpcs.push_back(new SubSpace<elemType>(*(rhs._vpSubSpcs[iSubSpc])));
  }
  _vInfoReqs.resize(rhs._vInfoReqs.size());
  _vUpdtReqs.resize(rhs._vUpdtReqs.size());
}

template <class elemType>
inline pRPL::Layer<elemType>::
~Layer() {
  cleanCellSpace();
  cleanNbrhood();
  cleanSubSpcs();
}

template <class elemType>
inline pRPL::Layer<elemType>& pRPL::Layer<elemType>::
operator=(const Layer<elemType> &rhs) {
  if(this == &rhs) {
    return *this;
  }
  
  _pPRPrc = rhs._pPRPrc;
  /*
  if(_name == "Untitled") {
    _name = rhs._name + "_copy";
  }
  */
  if(rhs._pCellSpace) {
    if(_pCellSpace) {
      *(_pCellSpace) = *(rhs._pCellSpace);
    }
    else {
      _pCellSpace = new CellSpace<elemType>(*(rhs._pCellSpace));
    }
  }
  else {
    cleanCellSpace();
  }
  
  if(rhs._pNbrhood) {
    if(_pNbrhood) {
      *(_pNbrhood) = *(rhs._pNbrhood);
    }
    else {
      _pNbrhood = new Neighborhood<elemType>(*(rhs._pNbrhood));
    }
  }
  else {
    cleanNbrhood();
  }
  
  if(nLclSubSpcs() != rhs.nLclSubSpcs()) {
    cleanSubSpcs();
    for(int iSubSpc = 0; iSubSpc < rhs.nLclSubSpcs(); iSubSpc++) {
      _vpSubSpcs.push_back(new SubSpace<elemType>(*(rhs._vpSubSpcs[iSubSpc])));
    }
  }
  else {
    for(int iSubSpc = 0; iSubSpc < rhs.nLclSubSpcs(); iSubSpc++) {
      *(_vpSubSpcs[iSubSpc]) = *(rhs._vpSubSpcs[iSubSpc]);
    }
  }
  
  _mGatherTypes = rhs._mGatherTypes;
  _vInfoReqs.resize(rhs._vInfoReqs.size());
  _vUpdtReqs.resize(rhs._vUpdtReqs.size());
  
  return *this;
}

template <class elemType>
inline const string& pRPL::Layer<elemType>::
name() const {
  return _name;
}

template <class elemType>
inline void pRPL::Layer<elemType>::
name(const string &layerName) {
  _name = layerName;
}

template <class elemType>
inline bool pRPL::Layer<elemType>::
hasPRPrc() const {
  bool hasIt = true;
  if(!_pPRPrc) {
    cerr << __FILE__ << __FUNCTION__ \
         << " Error: no prPrc were found on the process" \
         << endl;
    hasIt = false;
  }
  return hasIt;
}

template <class elemType>
inline pRPL::PRProcess *const pRPL::Layer<elemType>::
prPrc() const {
  return _pPRPrc;
}

template <class elemType>
inline int pRPL::Layer<elemType>::
id() const {
  int myID = -1;
  if(hasPRPrc()) {
    myID = _pPRPrc->id();
  }
  return myID;
}

template <class elemType>
inline const string pRPL::Layer<elemType>::
title() const {
  ostringstream myTitle;
  myTitle << _name << id();
  return myTitle.str();
}

template <class elemType>
inline bool pRPL::Layer<elemType>::
isMaster() const {
  return (hasPRPrc() &&
          _pPRPrc->isMaster());
}

template <class elemType>
void pRPL::Layer<elemType>::
cleanCellSpace() {
  if(_pCellSpace) {
    delete _pCellSpace;
    _pCellSpace = 0;
  }
}

template <class elemType>
void pRPL::Layer<elemType>::
cleanNbrhood() {
  if(_pNbrhood) {
    delete _pNbrhood;
    _pNbrhood = 0;
  }
}

template <class elemType>
void pRPL::Layer<elemType>::
cleanSubSpcs() {
  int nSubSpcs = _vpSubSpcs.size();
  for(int iSub = 0; iSub < nSubSpcs; iSub++) {
    if(_vpSubSpcs[iSub]) {
      delete _vpSubSpcs[iSub];
    }
  }
  _vpSubSpcs.clear();
}

template <class elemType>
bool pRPL::Layer<elemType>::
hasCellSpace() const {
  bool hasIt = true;
  if(!_pCellSpace ||
     _pCellSpace->empty()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: no CellSpace associated with layer[" \
         << title() << "]" << endl;
    hasIt = false;
  }
  return hasIt;
}

template <class elemType>
bool pRPL::Layer<elemType>::
hasNbrhood() const {
  bool hasIt = true;
  if(!_pNbrhood) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: no neighborhood associated with layer[" \
         << title() << "]" << endl;
    hasIt = false;
  }
  return hasIt;
}

template <class elemType>
void pRPL::Layer<elemType>::
prPrc(PRProcess &prPrc) {
  _pPRPrc = &prPrc;
}

template <class elemType>
bool pRPL::Layer<elemType>::
newCellSpace() {
  if(!hasPRPrc()) {
    return false;
  }
  cleanCellSpace();
  if(isMaster()) {
    _pCellSpace = new CellSpace<elemType>();
    if(!_pCellSpace) {
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: unable to new a CellSpace" \
           << endl;
      return false;
    }
  }
  return true;
}

template <class elemType>
bool pRPL::Layer<elemType>::
newCellSpace(const SpaceDims& dims) {
  if(!hasPRPrc()) {
    return false;
  }
  cleanCellSpace();
  if(isMaster()) {
    _pCellSpace = new CellSpace<elemType>(dims);
    if(!_pCellSpace) {
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: unable to new a CellSpace" \
           << endl;
      return false;
    }
  }
  return true;
}

template <class elemType>
bool pRPL::Layer<elemType>::
newCellSpace(const SpaceDims& dims,
             const elemType &initVal) {
  if(!hasPRPrc()) {
    return false;
  }
  cleanCellSpace();
  if(isMaster()) {
    _pCellSpace = new CellSpace<elemType>(dims, initVal);
    if(!_pCellSpace) {
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: unable to new a CellSpace" \
           << endl;
      return false;
    }
  }
  return true;
}

template <class elemType>
bool pRPL::Layer<elemType>::
newCellSpace(int nRows, int nCols) {
  return newCellSpace(SpaceDims(nRows, nCols));
}

template <class elemType>
bool pRPL::Layer<elemType>::
newCellSpace(int nRows, int nCols,
             const elemType &initVal) {
  return newCellSpace(SpaceDims(nRows, nCols), initVal);
}

template <class elemType>
bool pRPL::Layer<elemType>::
newNbrhood() {
  if(!hasPRPrc()) {
    return false;
  }
  
  cleanNbrhood();
  _pNbrhood = new Neighborhood<elemType>();
  if(!_pNbrhood) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to new a Neighborhood on process[" \
         << title() << "]" << endl;
    return false;
  }
  return true;
}

template <class elemType>
bool pRPL::Layer<elemType>::
newNbrhood(const vector<CellCoord> &vNbrCoords,
           double weight) {
  if(!hasPRPrc()) {
    return false;
  }
  
  cleanNbrhood();
  _pNbrhood = new Neighborhood<elemType>(vNbrCoords, weight);
  if(!_pNbrhood) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to new a Neighborhood on process[" \
         << title() << "]" << endl;
    return false;
  }
  return true;
}

template <class elemType>
bool pRPL::Layer<elemType>::
newNbrhood(const vector<CellCoord> &vNbrCoords,
           const vector<double> &vNbrWeights) {
  if(!hasPRPrc()) {
    return false;
  }
  
  cleanNbrhood();
  _pNbrhood = new Neighborhood<elemType>(vNbrCoords, vNbrWeights);
  if(!_pNbrhood) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to new a Neighborhood on process[" \
         << title() << "]" << endl;
    return false;
  }
  return true;
}

template<class elemType> template<class elemType2>
bool pRPL::Layer<elemType>::
newNbrhood(const Neighborhood<elemType2> &nbr) {
  if(!hasPRPrc()) {
    return false;
  }
  
  cleanNbrhood();
  _pNbrhood = new Neighborhood<elemType>(nbr);
  if(!_pNbrhood) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to new a Neighborhood on process[" \
         << title() << "]" << endl;
    return false;
  }
  return true;
}
      
template <class elemType>
bool pRPL::Layer<elemType>::
smplDcmp(DomDcmpMethod dcmpMethod,
         int nSubSpcs1, int nSubSpcs2) {
  if(!hasPRPrc()) {
    return false;
  }
  if(!_pPRPrc->active()) {
    return true;
  }
  if(!isMaster()) {
    return true;
  }

  if(!hasCellSpace() || !hasNbrhood()) {
    return false;
  }
  const SpaceDims &glbDims = _pCellSpace->dims();

  SmplDcmp<elemType> dcmp(glbDims, *_pNbrhood);
  SubInfoVect &glbSubInfos = _pPRPrc->glbSubInfos();
  switch(dcmpMethod) {
    case SMPL_ROW: /* Row-wise */
      if(!dcmp.rowDcmp(glbSubInfos, nSubSpcs1)) {
        return false;
      }
      break;
    case SMPL_COL: /* Column-wise */
      if(!dcmp.colDcmp(glbSubInfos, nSubSpcs1)) {
        return false;
      }
      break;
    case SMPL_BLK: /* Block */
      if(!dcmp.blockDcmp(glbSubInfos, nSubSpcs1, nSubSpcs2)) {
        return false;
      }
      break;
    default:
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: invalid decomposition method(" << dcmpMethod << ")." \
           << endl;
      return false;
  }
  return true;
}

template <class elemType>
bool pRPL::Layer<elemType>::
quadDcmp(Transition<elemType> &transition,
         int maxNumLeaves,
         int minWorkload) {
  if(!hasPRPrc()) {
    return false;
  }
  if(!_pPRPrc->active()) {
    return true;
  }
  if(!isMaster()) {
    return true;
  }
  
  if(!hasCellSpace() || !hasNbrhood()) {
    return false;
  }
  
  QuadTree<elemType> quadtree(*_pCellSpace, *_pNbrhood, transition,
                              maxNumLeaves, minWorkload);
  bool keepGrowing = true;
  while(keepGrowing) {
    keepGrowing = quadtree.grow();
  }
  quadtree.updateNbrs();
  quadtree.dcmp(_pPRPrc->glbSubInfos());
  /*
  for(int iLeaf = 0; iLeaf < quadtree.nLeaves(); iLeaf++) {
    cout << "leave[" << quadtree.leaf(iLeaf)->id() \
         << "] workload = " << quadtree.leaf(iLeaf)->workload() << endl;
  }
  */
  return true;
}

template <class elemType>
bool pRPL::Layer<elemType>::
bcastNbrhood() {
  if(!hasPRPrc()) {
    return false;
  }
  if(!_pPRPrc->active()) {
    return true;
  }
  
  const MPI_Comm &comm =  _pPRPrc->comm();
  int masterID = _pPRPrc->masterID();

  vector<CellCoord> vNbrCoords;
  vector<double> vNbrWeights;
  double weight = 0.0;
  double info[3]; /* 0 - number of nbrs
                     1 - if euqally weighted
                     2 - weight */

  if(isMaster()) {
    if(!_pNbrhood){
      info[0] = 0.0;
      info[1] = 0.0;
      info[2] = 0.0;
    }
    else {
      if(_pNbrhood->isEquallyWeighted(weight)) {
        _pNbrhood->toVect(vNbrCoords);
        info[1] = 1.0;
        info[2] = weight;
      }
      else {
        _pNbrhood->toVects(vNbrCoords, vNbrWeights);
        info[1] = 0.0;
        info[2] = 0.0;
      }
      info[0] = (double)vNbrCoords.size();
    }
  } // End of Master process
  else {
    cleanNbrhood();
  }

  MPI_Bcast(info, 3, MPI_DOUBLE, masterID, comm);

  int nNbrs = (int)(info[0]);
  if(nNbrs == 0) {
    return true;
  }
  int isEquallyWeighted = (int)(info[1]);

  if(!isMaster()) {
    vNbrCoords.resize(nNbrs);
    if(isEquallyWeighted == 0) {
      vNbrWeights.resize(nNbrs);
    }
  }

  MPI_Bcast(&(vNbrCoords[0]), nNbrs*sizeof(CellCoord), MPI_CHAR,
            masterID, comm);

  if(isEquallyWeighted == 0) {
    MPI_Bcast(&(vNbrWeights[0]), nNbrs, MPI_DOUBLE, 
              masterID, comm);
  }

  if(!isMaster()) {
    if(isEquallyWeighted == 1) {
      weight = info[2];
      if(!newNbrhood(vNbrCoords, weight)) {
        return false;
      }
    }
    else if(isEquallyWeighted == 0){
      if(!newNbrhood(vNbrCoords, vNbrWeights)) {
        return false;
      }
    }
  }

  return true;
}

template <class elemType>
bool pRPL::Layer<elemType>::
dstrbtCellSpace(bool dstrbtData) {
  if(!hasPRPrc()) {
    return false;
  }
  if(!_pPRPrc->active()) {
    return true;
  }
  
  const MPI_Comm &comm = _pPRPrc->comm();
  int masterID = _pPRPrc->masterID();
  
  if(!_pPRPrc->hasPrcInfoMap()) {
    cerr << __FILE__ << __FUNCTION__ \
         << " Error: no PrcInfoMap on process[" \
         << _pPRPrc->id() << "], or the PrcInfoMap has NOT been mapped yet" \
         << endl;
    MPI_Abort(comm, -1);
    return false;
  }
  SubSpaceMap &prcInfoMap = *(_pPRPrc->prcInfoMap());
  const vector<SubSpaceInfo *> &lclSubInfos = _pPRPrc->lclSubInfos();
  int nTotalNbrs = 0;
  cleanSubSpcs();
  for(int iSub = 0; iSub < lclSubInfos.size(); iSub++){
    _vpSubSpcs.push_back(new SubSpace<elemType>(*(lclSubInfos[iSub])));
    nTotalNbrs += lclSubInfos[iSub]->nTotNbrs();
  }
  _vInfoReqs.resize(2 * nTotalNbrs);
  _vUpdtReqs.resize(2 * nTotalNbrs);
  int nSubSpcs = _vpSubSpcs.size();
  
  vector<MPI_Request> vRequests;
  vector<MPI_Status> vStatus;

  if(isMaster()) { /* Master process */
    if(!dstrbtData) {
      return true;
    }

    if(!hasCellSpace()) {
      MPI_Abort(comm, -1);
      return false;
    }
    CellSpace<elemType> &cellSpace = *_pCellSpace;
    int nGlbCols = cellSpace.nCols();

    int nPrcs = _pPRPrc->nPrcs();
    int nTotSubSpcs = _pPRPrc->nTotSubSpcs();
    vRequests.resize(nTotSubSpcs-nSubSpcs);
    MPI_Datatype aNewTypes[nTotSubSpcs-nSubSpcs];
    int iReq = 0;

    for(int iPrc = 0; iPrc < nPrcs; iPrc++) {
      const vector<SubSpaceInfo *> &vpInfos = prcInfoMap[iPrc];
      for(int iSub = 0; iSub < vpInfos.size(); iSub++) {
        SubSpaceInfo &subInfo = *(vpInfos[iSub]);
        if(cellSpace.dims() != subInfo.glbDims()) {
          cerr << __FILE__ << " " << __FUNCTION__ \
               << " Error: the dimensions of the cellSpace (" \
               << cellSpace.dims() << ") within layer[" 
               << title() << "] do NOT match the subInfo[" << subInfo.id() \
               << "]'s global dimensions ("\
               << subInfo.glbDims() << ")" << endl;
          MPI_Abort(comm, -1);
          return false;
        }
        int nRows = subInfo.nRows();
        int nCols = subInfo.nCols();
        int iRowBegin = subInfo.iRowBegin();
        int iColBegin = subInfo.iColBegin();
        
        if(iPrc != id()) {
          MPI_Type_vector(nRows, nCols*sizeof(elemType), 
                          nGlbCols*sizeof(elemType), MPI_CHAR, &(aNewTypes[iReq]));
          MPI_Type_commit(&(aNewTypes[iReq]));
          MPI_Isend(&(cellSpace[iRowBegin][iColBegin]), 1, aNewTypes[iReq],
                    iPrc, iPrc+iSub+2, comm, &(vRequests[iReq]));
          iReq++;
        }
        else {
          SubSpace<elemType> &subSpace = *(_vpSubSpcs[iSub]);
          int iLclRow, iGlbRow;
          for(iLclRow = 0, iGlbRow = iRowBegin; iLclRow < nRows; iLclRow++, iGlbRow++) {
            memcpy(subSpace[iLclRow], 
                   &(cellSpace[iGlbRow][iColBegin]),
                   nCols*sizeof(elemType));
          }
        }
      } // End of iSub loop
    } // End of iPrc loop
    
    vStatus.resize(iReq);
    MPI_Waitall(iReq, &(vRequests[0]), &(vStatus[0]));

    for(int iType = 0; iType < iReq; iType++) {
      MPI_Type_free(&(aNewTypes[iType]));
    }
  }
  else { /* Slave process */
    if(nSubSpcs > 0 && dstrbtData) {
      vRequests.resize(nSubSpcs);
      vStatus.resize(nSubSpcs);
      for(int iSubSpc = 0; iSubSpc < nSubSpcs; iSubSpc++) {
        SubSpace<elemType> &subSpace = *(_vpSubSpcs[iSubSpc]);
        int subSize = subSpace.size();
        MPI_Irecv(subSpace[0], subSize*sizeof(elemType), MPI_CHAR,
                  masterID, id()+iSubSpc+2, comm, &(vRequests[iSubSpc]));
      }
      MPI_Waitall(nSubSpcs, &(vRequests[0]), &(vStatus[0]));
    }
  }

  return true;
}

template <class elemType>
bool pRPL::Layer<elemType>::
bcastCellSpace(bool bcastData) {
  if(!hasPRPrc()) {
    return false;
  }
  if(!_pPRPrc->active()) {
    return true;
  }
  
  const MPI_Comm &comm = _pPRPrc->comm();
  int masterID = _pPRPrc->masterID();
  
  int spcDims[2];
  if(isMaster()) {
    if(!hasCellSpace()) {
      MPI_Abort(comm, -1);
      return false;
    }
    spcDims[0] = _pCellSpace->nRows();
    spcDims[1] = _pCellSpace->nCols();
  }
  else {
    cleanCellSpace();
  }

  MPI_Bcast(spcDims, 2, MPI_INT, masterID, comm);

  if(!isMaster()) {
    _pCellSpace = new CellSpace<elemType>(spcDims[0],
                                          spcDims[1]);
  }

  if(bcastData) {
    CellSpace<elemType> &cellSpc = *_pCellSpace;
    MPI_Bcast(cellSpc[0], cellSpc.size() * sizeof(elemType),
              MPI_CHAR, masterID, comm);
  }

  return true;
}

template <class elemType>
bool pRPL::Layer<elemType>::
distribute(bool dstrbtData) {
  return (bcastNbrhood() &&
          dstrbtCellSpace(dstrbtData));
}

template <class elemType>
bool pRPL::Layer<elemType>::
broadcast(bool bcastData) {
  return (bcastNbrhood() &&
          bcastCellSpace(bcastData));
}

template <class elemType>
bool pRPL::Layer<elemType>::
smplDcmpDstrbt(DomDcmpMethod dcmpMethod,
               int nSubSpcs1,
               int nSubSpcs2,
               bool dstrbtData) {
  if(!hasPRPrc()) {
    return false;
  }
  if(!_pPRPrc->active()) {
    return true;
  }
  
  if(!_pPRPrc->hasPrcInfoMap() &&
     !_pPRPrc->hasSubIDMap()) {
    if(!smplDcmp(dcmpMethod, nSubSpcs1, nSubSpcs2)) {
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: failed to decompose layer[" << title() << "]" \
           << endl;
      _pPRPrc->abort();
      return false;
    }
    if(!_pPRPrc->bcastSubInfos()) {
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: failed to broadcast the SubSpaceInfos onto process[" \
           << _pPRPrc->id() << "]" << endl;
      _pPRPrc->abort();
      return false;
    }
    if(!_pPRPrc->mapping()) {
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: failed to map the SubSpaceInfos on process[" \
           << _pPRPrc->id() << "]" << endl;
      _pPRPrc->abort();
      return false;
    }
  }
  
  if(!bcastNbrhood()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: failed to broadcast the Neighborhood onto process[" \
         << _pPRPrc->id() << "]" << endl;
      _pPRPrc->abort();
      return false;
  }
  if(!dstrbtCellSpace(dstrbtData)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: failed to distribute the SubSpaces onto process[" \
         << _pPRPrc->id() << "]" << endl;
      _pPRPrc->abort();
      return false;
  }
          
  _pPRPrc->sync();
  return true;
}

template <class elemType>
bool pRPL::Layer<elemType>::
quadDcmpDstrbt(Transition<elemType> &transition,
               int maxNumLeaves,
               int minWorkload,
               bool dstrbtData) {
  if(!hasPRPrc()) {
    return false;
  }
  if(!_pPRPrc->active()) {
    return true;
  }
  
  if(!_pPRPrc->hasPrcInfoMap() &&
     !_pPRPrc->hasSubIDMap()) {
    if(!quadDcmp(transition, maxNumLeaves, minWorkload)) {
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: failed to decompose layer[" << title() << "]" \
           << endl;
      _pPRPrc->abort();
      return false;
    }
    if(!_pPRPrc->bcastSubInfos()) {
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: failed to broadcast the SubSpaceInfos onto process[" \
           << _pPRPrc->id() << "]" << endl;
      _pPRPrc->abort();
      return false;
    }
    if(!_pPRPrc->mapping()) {
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: failed to map the SubSpaceInfos on process[" \
           << _pPRPrc->id() << "]" << endl;
      _pPRPrc->abort();
      return false;
    }
  }
  
  if(!bcastNbrhood()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: failed to broadcast the Neighborhood onto process[" \
         << _pPRPrc->id() << "]" << endl;
      _pPRPrc->abort();
      return false;
  }
  if(!dstrbtCellSpace(dstrbtData)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: failed to distribute the SubSpaces onto process[" \
         << _pPRPrc->id() << "]" << endl;
      _pPRPrc->abort();
      return false;
  }
          
  _pPRPrc->sync();
  return true;
}

template <class elemType>
inline pRPL::CellSpace<elemType>* pRPL::Layer<elemType>::
cellSpace() {
  return _pCellSpace;
}

template <class elemType>
inline const pRPL::CellSpace<elemType>* pRPL::Layer<elemType>::
cellSpace() const {
  return _pCellSpace;
}

template <class elemType>
inline pRPL::Neighborhood<elemType>* pRPL::Layer<elemType>::
nbrhood() {
  return _pNbrhood;
}

template <class elemType>
inline const pRPL::Neighborhood<elemType>* pRPL::Layer<elemType>::
nbrhood() const {
  return _pNbrhood;
}

template <class elemType>
inline int pRPL::Layer<elemType>::
nLclSubSpcs() const {
  return _vpSubSpcs.size();
}

template <class elemType>
pRPL::SubSpace<elemType>* pRPL::Layer<elemType>::
operator[](int iSubSpc) {
  return _vpSubSpcs.at(iSubSpc);
}

template <class elemType>
const pRPL::SubSpace<elemType>* pRPL::Layer<elemType>::
operator[](int iSubSpc) const {
  return _vpSubSpcs.at(iSubSpc);
}

template <class elemType>
pRPL::SubSpace<elemType>* pRPL::Layer<elemType>::
findSubSpc(int subID) {
  SubSpace<elemType> *pSub = 0;
  for(int iSub = 0; iSub < nLclSubSpcs(); iSub++) {
    if(_vpSubSpcs[iSub]->id() == subID) {
      pSub = _vpSubSpcs[iSub];
      break;
    }
  }
  return pSub;
}

template <class elemType>
const pRPL::SubSpace<elemType>* pRPL::Layer<elemType>::
findSubSpc(int subID) const {
  const SubSpace<elemType> *pSub = 0;
  for(int iSub = 0; iSub < nLclSubSpcs(); iSub++) {
    if(_vpSubSpcs[iSub]->id() == subID) {
      pSub = _vpSubSpcs[iSub];
      break;
    }
  }
  return pSub;
}

template <class elemType>
double pRPL::Layer<elemType>::
sizeRatio() const {
  double ratio = 0.0;
  for(int iSub = 0; iSub < nLclSubSpcs(); iSub++) {
    ratio += _vpSubSpcs[iSub].sizeRatio();
  }
  return ratio;
}

template <class elemType>
void pRPL::Layer<elemType>::
resetExchange() {
  _mSendSizes.clear();
  _mSendLocs.clear();
  _mRecvSizes.clear();
  _mRecvLocs.clear();
  _iUpdtInfo = 0;
  _iUpdtPack = 0;
}

template <class elemType>
bool pRPL::Layer<elemType>::
exchangeBegin(SubSpace<elemType> *pSubSpc) {
  if(!hasPRPrc()) {
    return false;
  }
  if(!_pPRPrc->active()) {
    return true;
  }
  
  const MPI_Comm &comm = _pPRPrc->comm();

  if(!hasNbrhood()) {
    MPI_Abort(comm, -1);
    return false;
  }

  if(!pSubSpc) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: missing pSubSpc" \
         << " on layer[" << title() << "]" << endl;
    MPI_Abort(comm, -1);
    return false;
  }
  SubSpace<elemType> &subSpc = *pSubSpc;

  if(!_pPRPrc->hasSubIDMap()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: missing subIDMap" \
         << " on process[" << _pPRPrc->id() << "]" << endl;
    MPI_Abort(comm, -1);
    return false;
  }
  SubIDMap &subIDMap = *(_pPRPrc->subIDMap());
  int maxSubID = subIDMap.maxID();
  
  int nNbrDirs = subSpc.nNbrDirs();
  for(int iDir = 0; iDir < nNbrDirs; iDir++) {
    int oppDir = subSpc.oppositeDir(iDir);
    MeshDir oppDirM = subSpc.spcDir2MeshDir(oppDir);
    /* Send updated cells */
    if(subSpc.hasNbrs(iDir) &&
       _pNbrhood->hasNbrs(oppDirM)) {
      const IntVect& nbrSpcIDs = subSpc.nbrSubSpcIDs(iDir);
      for(int iNbr = 0; iNbr < nbrSpcIDs.size(); iNbr++) {
        int nbrSpcID = nbrSpcIDs[iNbr];
        int nbrPrcID = subIDMap[nbrSpcID];
        CharVect& stream2Send = subSpc.stream2Send(iDir, iNbr);
        _mSendSizes[_iUpdtInfo] = make_pair(id(), stream2Send.size());
        MPI_Isend(&(_mSendSizes[_iUpdtInfo]), 2, MPI_INT,
                  nbrPrcID, nbrSpcID,
                  comm, &(_vInfoReqs[_iUpdtInfo]));
        if(_mSendSizes[_iUpdtInfo].second > 0) {
          _mSendLocs[_iUpdtPack] = &(stream2Send[0]);
          MPI_Isend(_mSendLocs[_iUpdtPack],
                    _mSendSizes[_iUpdtInfo].second,
                    MPI_CHAR, nbrPrcID, nbrSpcID+maxSubID+1,
                    comm, &(_vUpdtReqs[_iUpdtPack]));
          _iUpdtPack++;
        }
        _iUpdtInfo++;
      }
    }
    
    /* Receive updated cells */
    if(subSpc.hasNbrs(oppDir) &&
       _pNbrhood->hasNbrs(oppDirM)) {
      const IntVect& nbrSpcIDs = subSpc.nbrSubSpcIDs(oppDir);
      for(int iNbr = 0; iNbr < nbrSpcIDs.size(); iNbr++) {
        int nbrSpcID = nbrSpcIDs[iNbr];
        int nbrPrcID = subIDMap[nbrSpcID];
        _mRecvSizes[_iUpdtInfo] = make_pair(pSubSpc, IntPair()); 
        MPI_Irecv(&(_mRecvSizes[_iUpdtInfo].second), 2, MPI_INT, 
                  nbrPrcID, subSpc.id(),
                  comm, &(_vInfoReqs[_iUpdtInfo]));
        _iUpdtInfo++;
      }
    }
  } // End of iDir loop
  
  return true;
}

template <class elemType>
bool pRPL::Layer<elemType>::
exchangeSize() {
  if(!hasPRPrc()) {
    return false;
  }
  if(!_pPRPrc->active()) {
    return true;
  }
  
  vector<MPI_Status> vInfoStatus(_vInfoReqs.size());
  MPI_Waitall(_iUpdtInfo, &(_vInfoReqs[0]), &(vInfoStatus[0]));

  const MPI_Comm &comm = _pPRPrc->comm();
  if(!_pPRPrc->hasSubIDMap()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: missing subIDMap" \
         << " on process[" << _pPRPrc->id() << "]" << endl;
    MPI_Abort(comm, -1);
    return false;
  }
  SubIDMap &subIDMap = *(_pPRPrc->subIDMap());
  int maxSubID = subIDMap.maxID();

  typename map<int, pair<SubSpace<elemType>*, IntPair> >
    ::iterator iRecv = _mRecvSizes.begin();
  while(iRecv != _mRecvSizes.end()) {
    SubSpace<elemType> *pSubSpc = (*iRecv).second.first;
    int packSize = (*iRecv).second.second.second;
    if(packSize > 0) {
      int crntBufSize = pSubSpc->recvBufSize();
      pSubSpc->recvBufSize(crntBufSize + packSize);
    }
    iRecv++;
  }

  map<int, int> mRecved;
  iRecv = _mRecvSizes.begin();
  while(iRecv != _mRecvSizes.end()) {
    SubSpace<elemType> *pSubSpc = (*iRecv).second.first;
    int subID = pSubSpc->id();
    int nbrPrcID = (*iRecv).second.second.first;
    int packSize = (*iRecv).second.second.second;
    if(packSize > 0) {
      int nRecved = mRecved[subID];
      CharVect &streamRecved = pSubSpc->streamRecved();
      _mRecvLocs[_iUpdtPack] = &(streamRecved[nRecved]);
      MPI_Irecv(_mRecvLocs[_iUpdtPack], packSize,
                MPI_CHAR, nbrPrcID, subID+maxSubID+1,
                comm, &(_vUpdtReqs[_iUpdtPack]));
      mRecved[subID] += packSize;
      _iUpdtPack++;
    }
    iRecv++;
  }

  return true;
}

template <class elemType>
bool pRPL::Layer<elemType>::
exchangeEnd(Transition<elemType> *pTransition) {
  if(!hasPRPrc()) {
    return false;
  }
  if(!_pPRPrc->active()) {
    return true;
  }
  
  vector<MPI_Status> vUpdtStatus(_vUpdtReqs.size());
  MPI_Waitall(_iUpdtPack, &(_vUpdtReqs[0]), &(vUpdtStatus[0]));

  for(int iSubSpc = 0; iSubSpc < nLclSubSpcs(); iSubSpc++) {
    SubSpace<elemType> &subSpc = *(_vpSubSpcs[iSubSpc]);
    if(!subSpc.loadStreamRecved(pTransition, _pNbrhood)) {
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: SubSpace[" << subSpc.id() \
           << "] on process[" << title() \
           << "] failed to load the received cells" << endl;
      return false;
    }
    subSpc.cleanStreamBuf();
  }

  return true;
}

template <class elemType>
bool pRPL::Layer<elemType>::
update(Transition<elemType> &transition) {
  if(!hasPRPrc()) {
    return false;
  }
  if(!_pPRPrc->active()) {
    return true;
  }
  
  const MPI_Comm &comm = _pPRPrc->comm();

  Transition<elemType> *pTrans = transition.needFinalize() ? 
                                 &transition : 0;
  if(transition.needExchange()) { // Exchange needed
    resetExchange();
    for(int iSubSpc = 0; iSubSpc < nLclSubSpcs(); iSubSpc++) {
      SubSpace<elemType> *pSubSpc = _vpSubSpcs[iSubSpc];
      if(transition.edgesFirst()) {
        if(!pSubSpc->calcBRs(&transition,
                             _pNbrhood,
                             &(_pPRPrc->glbSubInfos())) ||
           !pSubSpc->updateEdges(&transition, _pNbrhood)) {
          MPI_Abort(comm, -1);
          return false;
        }
      }
      else {
        if(!pSubSpc->updateAll(&transition, _pNbrhood)) {
          MPI_Abort(comm, -1);
          return false;
        }
      }
      if(!exchangeBegin(pSubSpc)) {
        MPI_Abort(comm, -1);
        return false;
      }
    }

    if(!exchangeSize()) {
      MPI_Abort(comm, -1);
      return false;
    }

    for(int iSubSpc = 0; iSubSpc < nLclSubSpcs(); iSubSpc++) {
      SubSpace<elemType> *pSubSpc = _vpSubSpcs[iSubSpc];
      if(transition.edgesFirst()) {
        if(!pSubSpc->updateInterior(&transition, _pNbrhood)) {
          MPI_Abort(comm, -1);
          return false;
        }
      }
      if(!pSubSpc->updateFinalize(pTrans, _pNbrhood)) {
        MPI_Abort(comm, -1);
        return false;
      }
    }

    if(!exchangeEnd(pTrans)) {
      MPI_Abort(comm, -1);
      return false;
    }
    resetExchange();
  }
  else { // No exchange
    for(int iSubSpc = 0; iSubSpc < nLclSubSpcs(); iSubSpc++) {
      SubSpace<elemType> *pSubSpc = _vpSubSpcs[iSubSpc];
      if(!pSubSpc->updateAll(&transition, _pNbrhood) ||
         !pSubSpc->updateFinalize(pTrans, _pNbrhood)) {
        MPI_Abort(comm, -1);
        return false;
      }
    }
  }
  
  return true;
}

template <class elemType>
bool pRPL::Layer<elemType>::
update(map<int, IntVect> &mLclIdxs2Eval,
       Transition<elemType> &transition) {
  if(!hasPRPrc()) {
    return false;
  }
  if(!_pPRPrc->active()) {
    return true;
  }
  
  const MPI_Comm &comm = _pPRPrc->comm();

  Transition<elemType> *pTrans = transition.needFinalize() ?
                                 &transition : 0;
  if(transition.needExchange()) { // Exchange needed
    resetExchange();
    for(int iSubSpc = 0; iSubSpc < nLclSubSpcs(); iSubSpc++) {
      SubSpace<elemType> *pSubSpc = _vpSubSpcs[iSubSpc];
      int subID = pSubSpc->id();
      if(transition.edgesFirst()) {
        if(!pSubSpc->calcBRs(&transition,
                             _pNbrhood,
                             &(_pPRPrc->glbSubInfos())) ||
           !pSubSpc->updateEdges(mLclIdxs2Eval[subID],
                                 &transition, _pNbrhood)) {
          MPI_Abort(comm, -1);
          return false;
        }
      }
      else {
        if(!pSubSpc->updateAll(mLclIdxs2Eval[subID],
                               &transition, _pNbrhood)) {
          return false;
        }
      }
      if(!exchangeBegin(pSubSpc)) {
        MPI_Abort(comm, -1);
        return false;
      }
    }

    if(!exchangeSize()) {
      MPI_Abort(comm, -1);
      return false;
    }
    
    for(int iSubSpc = 0; iSubSpc < nLclSubSpcs(); iSubSpc++) {
      SubSpace<elemType> *pSubSpc = _vpSubSpcs[iSubSpc];
      if(transition.edgesFirst()) {
        int subID = pSubSpc->id();
        if(!pSubSpc->updateInterior(mLclIdxs2Eval[subID],
                                    &transition, _pNbrhood)) {
          MPI_Abort(comm, -1);
          return false;
        }
      }
      if(!pSubSpc->updateFinalize(pTrans, _pNbrhood)) {
        MPI_Abort(comm, -1);
        return false;
      }
    }

    if(!exchangeEnd(pTrans)) {
      MPI_Abort(comm, -1);
      return false;
    }
    resetExchange();
  }
  else { // No exchange
    for(int iSubSpc = 0; iSubSpc < nLclSubSpcs(); iSubSpc++) {
      SubSpace<elemType> *pSubSpc = _vpSubSpcs[iSubSpc];
      int subID = pSubSpc->id();
      if(!pSubSpc->updateAll(mLclIdxs2Eval[subID],
                             &transition, _pNbrhood) ||
         !pSubSpc->updateFinalize(pTrans, _pNbrhood)) {
        MPI_Abort(comm, -1);
        return false;
      }
    }
  }

  return true;
}

template <class elemType>
bool pRPL::Layer<elemType>::
update(IntVect &vGlbIdxs2Eval,
       Transition<elemType> &transition) {
  if(!hasPRPrc()) {
    return false;
  }
  if(!_pPRPrc->active()) {
    return true;
  }
  
  map<int, IntVect> mLclIdxs2Eval;
  for(int iSubSpc = 0; iSubSpc < nLclSubSpcs(); iSubSpc++) {
    SubSpace<elemType> *pSubSpc = _vpSubSpcs[iSubSpc];
    const CoordBR &workBR = pSubSpc->workBR();
    int subID = pSubSpc->id();
    IntVect::iterator iIdx = vGlbIdxs2Eval.begin();
    while(iIdx != vGlbIdxs2Eval.end()) {
      int lclIdx = pSubSpc->glbIdx2lclIdx(*iIdx);
      CellCoord lclCoord = pSubSpc->idx2coord(lclIdx);
      if(workBR.contain(lclCoord)) {
        mLclIdxs2Eval[subID].push_back(lclIdx);
        vGlbIdxs2Eval.erase(iIdx);
      }
      else {
        iIdx++;
      }
    }
  }
  
  return (update(mLclIdxs2Eval, transition));
}

template <class elemType>
bool pRPL::Layer<elemType>::
buildGatherTypes() {
  if(!hasPRPrc()) {
    return false;
  }
  if(!_pPRPrc->active()) {
    return true;
  }
  
  freeGatherTypes();
  if(isMaster()) { // Master process
    if(!hasCellSpace()) {
      return false;
    }
    CellSpace<elemType> &cellSpace = *_pCellSpace;
    int nGlbCols = cellSpace.nCols();

    int nTotSubSpcs = _pPRPrc->nTotSubSpcs();
    const SubInfoVect &glbSubInfos = _pPRPrc->glbSubInfos();
    for(int iInfo = 0; iInfo < nTotSubSpcs; iInfo++) {
      const SubSpaceInfo &subInfo = *(glbSubInfos[iInfo].first);
      int subID = subInfo.id();
      const CoordBR &workBR = subInfo.workBR();
      int nWorkRows = workBR.nRows();
      int nWorkCols = workBR.nCols();
      MPI_Type_vector(nWorkRows, nWorkCols*sizeof(elemType), 
                      nGlbCols*sizeof(elemType), MPI_CHAR,
                      &(_mGatherTypes[subID]));
      MPI_Type_commit(&(_mGatherTypes[subID]));
    }
  }
  else {
    int nSubSpcs = nLclSubSpcs();
    for(int iSubSpc = 0; iSubSpc < nSubSpcs; iSubSpc++) {
      SubSpace<elemType> &subSpc = *(_vpSubSpcs[iSubSpc]);
      int subID = subSpc.id();
      int nTotCols = subSpc.nCols();
      const CoordBR &workBR = subSpc.workBR();
      int nWorkRows = workBR.nRows();
      int nWorkCols = workBR.nCols();
      MPI_Type_vector(nWorkRows, nWorkCols*sizeof(elemType), 
                      nTotCols*sizeof(elemType), MPI_CHAR,
                      &(_mGatherTypes[subID]));
      MPI_Type_commit(&(_mGatherTypes[subID]));
    }
  }
  return true;
}

template <class elemType>
void pRPL::Layer<elemType>::
freeGatherTypes() {
  map<int, MPI_Datatype>::iterator iType = _mGatherTypes.begin();
  while(iType != _mGatherTypes.end()) {
    MPI_Type_free(&(iType->second));
    iType++;
  }
  _mGatherTypes.clear();
}

template <class elemType>
bool pRPL::Layer<elemType>::
gatherCellSpace() {
  if(!hasPRPrc()) {
    return false;
  }
  if(!_pPRPrc->active()) {
    return true;
  }
  
  const MPI_Comm &comm = _pPRPrc->comm();
  int masterID = _pPRPrc->masterID();

  vector<MPI_Request> vReqs;
  vector<MPI_Status> vStatus;

  if(isMaster()) { // Master process
    if(!hasCellSpace()) {
      MPI_Abort(comm, -1);
      return false;
    }
    CellSpace<elemType> &cellSpace = *_pCellSpace;

    if(!_pPRPrc->hasSubIDMap()) {
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: missing subIDMap" \
           << " on process[" << _pPRPrc->id() << "]" << endl;
      MPI_Abort(comm, -1);
      return false;
    }
    SubIDMap &subIDMap = *(_pPRPrc->subIDMap());

    int nTotSubSpcs = _pPRPrc->nTotSubSpcs();
    if(nTotSubSpcs != subIDMap.size()) {
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: on the master process, the size of the global" \
           << " SubSpaceInfo vector ("\
           << nTotSubSpcs << ") is NOT the same as that of" \
           << " that of the SubSpace ID map (" \
           << subIDMap.size() << ")" << endl;
      MPI_Abort(comm, -1);
      return false;
    }
    if(nTotSubSpcs != _mGatherTypes.size()) {
      buildGatherTypes();
    }

    vReqs.resize(nTotSubSpcs);
    const SubInfoVect &glbSubInfos = _pPRPrc->glbSubInfos();
    int iReq = 0;
    for(int iSubSpc = 0; iSubSpc < nTotSubSpcs; iSubSpc++) {
      const SubSpaceInfo &subInfo = *(glbSubInfos[iSubSpc].first);
      int subID = subInfo.id();
      int prcID = subIDMap[subID];
      const CoordBR &workBR = subInfo.workBR();
      CellCoord workBeginCoord = subInfo.lclCoord2glbCoord(workBR.nwCorner());
      if(prcID != _pPRPrc->id()) {
        MPI_Irecv(&(cellSpace[workBeginCoord.iRow()][workBeginCoord.iCol()]),
                  1, _mGatherTypes[subID], prcID, subID,
                  comm, &(vReqs[iReq]));
        iReq++;
      }
      else {
        SubSpace<elemType> *pSubSpc = findSubSpc(subID);
        if(pSubSpc == 0) {
          cerr << __FILE__ << " " << __FUNCTION__ \
               << " Error: unable to find SubSpace[" << subID \
               << "] on Layer[" << title() << "]" << endl;
          MPI_Abort(comm, -1);
          return false;
        }
        SubSpace<elemType> &subSpc = *pSubSpc;
        int nWorkRows = workBR.nRows();
        int nWorkCols = workBR.nCols();
        int iLclRow, iGlbRow;
        for(iLclRow = workBR.minIRow(), iGlbRow = workBeginCoord.iRow();
            iLclRow <= workBR.maxIRow(); iLclRow++, iGlbRow++) {
          memcpy(&(cellSpace[iGlbRow][workBeginCoord.iCol()]),
                 &(subSpc[iLclRow][workBR.minICol()]),
                 nWorkCols*sizeof(elemType));
        }
      }
    }
    vStatus.resize(iReq);
    MPI_Waitall(iReq, &(vReqs[0]), &(vStatus[0]));
  }
  else { // Slave process
    int nSubSpcs = _vpSubSpcs.size();
    if(nSubSpcs != _mGatherTypes.size()) {
      buildGatherTypes();
    }

    vReqs.resize(nSubSpcs);
    vStatus.resize(nSubSpcs);
    for(int iSubSpc = 0; iSubSpc < nSubSpcs; iSubSpc++) {
      SubSpace<elemType> &subSpc = *(_vpSubSpcs[iSubSpc]);
      int subID = subSpc.id();
      //int iRowWB = subSpc.iRowWorkBegin();
      int iRowWB = subSpc.workBR().minIRow();
      //int iColWB = subSpc.iColWorkBegin();
      int iColWB = subSpc.workBR().minICol();
      MPI_Isend(&(subSpc[iRowWB][iColWB]), 1, _mGatherTypes[subID],
                masterID, subID, comm, &(vReqs[iSubSpc]));
    }
    MPI_Waitall(nSubSpcs, &(vReqs[0]), &(vStatus[0]));
  }

  return true;
}

#endif
