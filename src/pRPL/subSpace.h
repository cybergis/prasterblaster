#ifndef SUBSPACE_H
#define SUBSPACE_H

/***************************************************************************
* subSpace.h
*
* Project: pRPL, v 1.0
* Purpose: Header file for class pRPL::SubSpace
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
#include "basicCell.h"
#include "cellSpace.h"
#include "neighborhood.h"
#include "subSpaceInfo.h"
#include "subInfoVect.h"
#include "mpi.h"

namespace pRPL {
  template <class elemType>
  class SubSpace: public CellSpace<elemType> {
    public:
      // Constructor
      SubSpace();
      SubSpace(const CellSpace<elemType> &cellSpace,
               const SubSpaceInfo& subSpaceInfo);
      SubSpace(const SubSpaceInfo& subSpaceInfo);
      SubSpace(const SubSpace<elemType> &rhs);

      // Deconstructor
      ~SubSpace();
      
      SubSpace<elemType>& operator=(const SubSpace<elemType> &rhs);
      
      bool hasInfo() const;
      int id() const;
      DomDcmpType domDcmpType() const;
      const SpaceDims& glbDims() const;
      const CoordBR& MBR() const;
      const CoordBR& workBR() const;
      double sizeRatio() const;
      
      int nNbrDirs() const;
      int nTotNbrs() const;
      bool hasNbrs(int iDir) const;
      int nNbrs(int iDir) const;
      const IntVect& nbrSubSpcIDs(int iDir) const;
      int nbrDir(int nbrID) const;
      int oppositeDir(int iDir) const;
      MeshDir spcDir2MeshDir(int iDir) const;

      const CellCoord lclCoord2glbCoord(const CellCoord& lclCoord) const;
      const CellCoord glbCoord2lclCoord(const CellCoord& glbCoord) const;
      int lclCoord2glbIdx(int iRowLcl, int iColLcl) const;
      int lclCoord2glbIdx(const CellCoord &lclCoord) const;
      const CellCoord glbIdx2lclCoord(int glbIdx) const;
      int glbIdx2lclIdx(int glbIdx) const;
      int glbCoord2glbIdx(int iRowGlb, int iColGlb) const;
      int glbCoord2glbIdx(const CellCoord &glbCoord) const;

      bool makeStream2Send();
      bool loadStreamRecved(Transition<elemType> *pTransition = 0,
                            Neighborhood<elemType> *pNbrhood = 0);
      void cleanStreamBuf();

      CharVect& stream2Send(int iDir, int iNbr);
      CharVect& streamRecved();
      int recvBufSize() const;
      bool recvBufSize(int nCells);

      bool calcBRs(const Transition<elemType> *pTransition,
                   const Neighborhood<elemType> *pNbrhood,
                   const SubInfoVect *pvSubInfos);
      bool updateEdges(Transition<elemType> *pTransition,
                       Neighborhood<elemType> *pNbrhood);
      bool updateInterior(Transition<elemType> *pTransition,
                          Neighborhood<elemType> *pNbrhood);
      bool updateAll(Transition<elemType> *pTransition,
                     Neighborhood<elemType> *pNbrhood);

      bool updateEdges(IntVect &vIdxs2Eval,
                       Transition<elemType> *pTransition,
                       Neighborhood<elemType> *pNbrhood);
      bool updateInterior(IntVect &vIdxs2Eval,
                          Transition<elemType> *pTransition,
                          Neighborhood<elemType> *pNbrhood);
      bool updateAll(IntVect &vIdxs2Eval,
                     Transition<elemType> *pTransition,
                     Neighborhood<elemType> *pNbrhood);

    protected:
     void _initBRs();
     void _rowwiseBRs(const Transition<elemType> *pTransition,
                      const Neighborhood<elemType> *pNbrhood);
     void _colwiseBRs(const Transition<elemType> *pTransition,
                      const Neighborhood<elemType> *pNbrhood);
     bool _blockwiseBRs(const Transition<elemType> *pTransition,
                        const Neighborhood<elemType> *pNbrhood,
                        const SubInfoVect *pvSubInfos);

    protected:
      const SubSpaceInfo *_pInfo;

      vector<vector<CoordBR> > _mSendBRs;
      vector<CoordBR> _vEdgeBRs;
      CoordBR _interiorBR;
      
      vector<vector<CharVect> > _mStream2Send;
      CharVect _vStreamRecved;
  };
};

/****************************************
*             Protected Methods         *
*****************************************/
template <class elemType>
void pRPL::SubSpace<elemType>::
_initBRs() {
  _mSendBRs.clear();
  _mSendBRs.resize(_pInfo->nNbrDirs());
  _vEdgeBRs.clear();
  _vEdgeBRs.resize(_pInfo->nEdges());
  _interiorBR = CoordBR();
}

template <class elemType>
void pRPL::SubSpace<elemType>::
_rowwiseBRs(const Transition<elemType> *pTransition,
            const Neighborhood<elemType> *pNbrhood) {
  int maxICol = pNbrhood->maxICol();
  int minICol = pNbrhood->minICol();
  int maxIRow = pNbrhood->maxIRow();
  int minIRow = pNbrhood->minIRow();
  
  int iRowWorkBegin = _pInfo->iRowWorkBegin();
  int iColWorkBegin = _pInfo->iColWorkBegin();
  int iRowWorkEnd = _pInfo->iRowWorkEnd();
  int iColWorkEnd = _pInfo->iColWorkEnd();
  const SpaceDims &dims = CellSpace<elemType>::_dims;
  
  if(_pInfo->hasNbrs(UPPER_DIR)) {
    _mSendBRs[UPPER_DIR].push_back(CoordBR());
  }
  if(_pInfo->hasNbrs(UPPER_DIR) && maxIRow > 0) {
    if(pTransition->onlyUpdtCtrCell()) {
      _mSendBRs[UPPER_DIR][0].nwCorner(iRowWorkBegin, iColWorkBegin);
      _mSendBRs[UPPER_DIR][0].seCorner(iRowWorkBegin - 1 + maxIRow,
                                       iColWorkEnd);
      _vEdgeBRs[UPPER_DIR].nwCorner(iRowWorkBegin, iColWorkBegin);
      _vEdgeBRs[UPPER_DIR].seCorner(iRowWorkBegin - 1 + maxIRow,
                                    iColWorkEnd);
      _interiorBR.nwCorner(iRowWorkBegin + maxIRow,
                           iColWorkBegin);
    }
    else {
      _mSendBRs[UPPER_DIR][0].nwCorner(0, 0);
      _mSendBRs[UPPER_DIR][0].seCorner(iRowWorkBegin - 1 + maxIRow,
                                       dims.nCols() - 1);
      _vEdgeBRs[UPPER_DIR].nwCorner(iRowWorkBegin, iColWorkBegin);
      _vEdgeBRs[UPPER_DIR].seCorner(iRowWorkBegin - 1 + maxIRow - minIRow,
                                    iColWorkEnd);
      _interiorBR.nwCorner(iRowWorkBegin + maxIRow - minIRow,
                           iColWorkBegin);
    }
  }
  else {
    _interiorBR.nwCorner(iRowWorkBegin, iColWorkBegin);
  }
    
  if(_pInfo->hasNbrs(LOWER_DIR)) {
    _mSendBRs[LOWER_DIR].push_back(CoordBR());
  }
  if(_pInfo->hasNbrs(LOWER_DIR) & minIRow < 0) {
    if(pTransition->onlyUpdtCtrCell()) {
      _mSendBRs[LOWER_DIR][0].nwCorner(iRowWorkEnd + 1 + minIRow,
                                       iColWorkBegin);
      _mSendBRs[LOWER_DIR][0].seCorner(iRowWorkEnd, iColWorkEnd);
      _vEdgeBRs[LOWER_DIR].nwCorner(iRowWorkEnd + 1 + minIRow,
                                    iColWorkBegin);
      _vEdgeBRs[LOWER_DIR].seCorner(iRowWorkEnd, iColWorkEnd);
      _interiorBR.seCorner(iRowWorkEnd + minIRow,
                           iColWorkEnd);
    }
    else {
      _mSendBRs[LOWER_DIR][0].nwCorner(iRowWorkEnd + 1 + minIRow,
                                       0);
      _mSendBRs[LOWER_DIR][0].seCorner(dims.nRows() - 1,
                                       dims.nCols() - 1);
      _vEdgeBRs[LOWER_DIR].nwCorner(iRowWorkEnd + 1 + minIRow - maxIRow,
                                    iColWorkBegin);
      _vEdgeBRs[LOWER_DIR].seCorner(iRowWorkEnd, iColWorkEnd);
      _interiorBR.seCorner(iRowWorkEnd + minIRow - maxIRow,
                           iColWorkEnd);
    }
  }
  else {
    _interiorBR.seCorner(iRowWorkEnd, iColWorkEnd);
  }
}

template <class elemType>
void pRPL::SubSpace<elemType>::
_colwiseBRs(const Transition<elemType> *pTransition,
            const Neighborhood<elemType> *pNbrhood) {
  int maxICol = pNbrhood->maxICol();
  int minICol = pNbrhood->minICol();
  int maxIRow = pNbrhood->maxIRow();
  int minIRow = pNbrhood->minIRow();
  
  int iRowWorkBegin = _pInfo->iRowWorkBegin();
  int iColWorkBegin = _pInfo->iColWorkBegin();
  int iRowWorkEnd = _pInfo->iRowWorkEnd();
  int iColWorkEnd = _pInfo->iColWorkEnd();
  const SpaceDims &dims = CellSpace<elemType>::_dims;
  
  if(_pInfo->hasNbrs(LEFT_DIR)) {
    _mSendBRs[LEFT_DIR].push_back(CoordBR());
  }
  if(_pInfo->hasNbrs(LEFT_DIR) && maxICol > 0) {
    if(pTransition->onlyUpdtCtrCell()) {
      _mSendBRs[LEFT_DIR][0].nwCorner(iRowWorkBegin, iColWorkBegin);
      _mSendBRs[LEFT_DIR][0].seCorner(iRowWorkEnd,
                                      iColWorkBegin - 1 + maxICol);
      _vEdgeBRs[LEFT_DIR].nwCorner(iRowWorkBegin, iColWorkBegin);
      _vEdgeBRs[LEFT_DIR].seCorner(iRowWorkEnd,
                                   iColWorkBegin - 1 + maxICol);
      _interiorBR.nwCorner(iRowWorkBegin,
                           iColWorkBegin + maxICol);
    }
    else {
      _mSendBRs[LEFT_DIR][0].nwCorner(0, 0);
      _mSendBRs[LEFT_DIR][0].seCorner(dims.nRows() - 1,
                                      iColWorkBegin - 1 + maxICol);
      _vEdgeBRs[LEFT_DIR].nwCorner(iRowWorkBegin, iColWorkBegin);
      _vEdgeBRs[LEFT_DIR].seCorner(iRowWorkEnd,
                                   iColWorkBegin - 1 + maxICol - minICol);
      _interiorBR.nwCorner(iRowWorkBegin,
                           iColWorkBegin + maxICol - minICol);
    }
  }
  else {
    _interiorBR.nwCorner(iRowWorkBegin, iColWorkBegin);
  }
    
  if(_pInfo->hasNbrs(RIGHT_DIR)) {
    _mSendBRs[RIGHT_DIR].push_back(CoordBR());
  }
  if(_pInfo->hasNbrs(RIGHT_DIR) && minICol < 0) {
    if(pTransition->onlyUpdtCtrCell()) {
      _mSendBRs[RIGHT_DIR][0].nwCorner(iRowWorkBegin,
                                       iColWorkEnd + 1 + minICol);
      _mSendBRs[RIGHT_DIR][0].seCorner(iRowWorkEnd, iColWorkEnd);
      _vEdgeBRs[RIGHT_DIR].nwCorner(iRowWorkBegin,
                                    iColWorkEnd + 1 + minICol);
      _vEdgeBRs[RIGHT_DIR].seCorner(iRowWorkEnd, iColWorkEnd);
      _interiorBR.seCorner(iRowWorkEnd,
                           iColWorkEnd + minICol);
    }
    else {
      _mSendBRs[RIGHT_DIR][0].nwCorner(0,
                                       iColWorkEnd + 1 + minICol);
      _mSendBRs[RIGHT_DIR][0].seCorner(dims.nRows() - 1,
                                       dims.nCols() - 1);
      _vEdgeBRs[RIGHT_DIR].nwCorner(iRowWorkBegin,
                                    iColWorkEnd + 1 + minICol - maxICol);
      _vEdgeBRs[RIGHT_DIR].seCorner(iRowWorkEnd, iColWorkEnd);
      _interiorBR.seCorner(iRowWorkEnd,
                           iColWorkEnd + minICol - maxICol);
    }
  }
  else {
    _interiorBR.seCorner(iRowWorkEnd, iColWorkEnd);
  }
}

template <class elemType>
bool pRPL::SubSpace<elemType>::
_blockwiseBRs(const Transition<elemType> *pTransition,
              const Neighborhood<elemType> *pNbrhood,
              const SubInfoVect *pvSubInfos) {
  int maxICol = pNbrhood->maxICol();
  int minICol = pNbrhood->minICol();
  int maxIRow = pNbrhood->maxIRow();
  int minIRow = pNbrhood->minIRow();
  
  int iRowWorkBegin = _pInfo->iRowWorkBegin();
  int iColWorkBegin = _pInfo->iColWorkBegin();
  int iRowWorkEnd = _pInfo->iRowWorkEnd();
  int iColWorkEnd = _pInfo->iColWorkEnd();
  
  int iRowBegin = _pInfo->iRowBegin();
  int iColBegin = _pInfo->iColBegin();
  int iRowEnd = _pInfo->iRowEnd();
  int iColEnd = _pInfo->iColEnd();

  int nRows = CellSpace<elemType>::nRows();
  int nCols = CellSpace<elemType>::nCols();

  for(int iDir = 0; iDir < 8; iDir++) {
    _mSendBRs[iDir].resize(_pInfo->nNbrs(iDir));
  }
  
  int iRowNIntr, iRowSIntr, iColWIntr, iColEIntr;
  if(_pInfo->hasNbrs(NORTH_DIR) && maxIRow > 0) {
    const IntVect &vNbrIDs = _pInfo->nbrSubSpcIDs(NORTH_DIR);
    for(int iNbr = 0; iNbr < vNbrIDs.size(); iNbr++) {
      int nbrSpcID = vNbrIDs[iNbr];
      const SubSpaceInfo *pNbrInfo = pvSubInfos->findSubInfo(nbrSpcID);
      if(!pNbrInfo) {
        cerr << __FILE__ << " " << __FUNCTION__ \
             << " Error: unable to find neighboring SubSpace [" \
             << nbrSpcID << "]" << endl;
        return false;
      }
      int iNbrColBegin = pNbrInfo->iColBegin();
      iNbrColBegin = (iNbrColBegin <= iColBegin) ?
                      0 : (iNbrColBegin - iColBegin);
      _mSendBRs[NORTH_DIR][iNbr].nwCorner(0,
                                          iNbrColBegin);
      int iNbrColEnd = pNbrInfo->iColEnd();
      iNbrColEnd = (iNbrColEnd >= iColEnd) ?
                    nCols - 1 : (iNbrColEnd - iColBegin);
      _mSendBRs[NORTH_DIR][iNbr].seCorner(iRowWorkBegin - 1 + maxIRow,
                                          iNbrColEnd);
    }
    if(pTransition->onlyUpdtCtrCell()) {
      _vEdgeBRs[NORTH_PDIR].nwCorner(iRowWorkBegin, iColWorkBegin);
      _vEdgeBRs[NORTH_PDIR].seCorner(iRowWorkBegin - 1 + maxIRow,
                                     iColWorkEnd);
      iRowNIntr = iRowWorkBegin + maxIRow;
    }
    else {
      _vEdgeBRs[NORTH_PDIR].nwCorner(iRowWorkBegin, iColWorkBegin);
      _vEdgeBRs[NORTH_PDIR].seCorner(iRowWorkBegin - 1 + maxIRow - minIRow,
                                     iColWorkEnd);
      iRowNIntr = iRowWorkBegin + maxIRow - minIRow;
    }
  }
  else {
    iRowNIntr = iRowWorkBegin;
  }
  
  if(_pInfo->hasNbrs(SOUTH_DIR) && minIRow < 0) {
    const IntVect &vNbrIDs = _pInfo->nbrSubSpcIDs(SOUTH_DIR);
    for(int iNbr = 0; iNbr < vNbrIDs.size(); iNbr++) {
      int nbrSpcID = vNbrIDs[iNbr];
      const SubSpaceInfo *pNbrInfo = pvSubInfos->findSubInfo(nbrSpcID);
      if(!pNbrInfo) {
        cerr << __FILE__ << " " << __FUNCTION__ \
             << " Error: unable to find neighboring SubSpace [" \
             << nbrSpcID << "]" << endl;
        return false;
      }
      int iNbrColBegin = pNbrInfo->iColBegin();
      iNbrColBegin = (iNbrColBegin <= iColBegin) ?
                      0 : (iNbrColBegin - iColBegin);
      _mSendBRs[SOUTH_DIR][iNbr].nwCorner(iRowWorkEnd + 1 + minIRow,
                                          iNbrColBegin);
      int iNbrColEnd = pNbrInfo->iColEnd();
      iNbrColEnd = (iNbrColEnd >= iColEnd) ?
                    nCols - 1 : (iNbrColEnd - iColBegin);
      _mSendBRs[SOUTH_DIR][iNbr].seCorner(nRows - 1,
                                          iNbrColEnd);
    }
    if(pTransition->onlyUpdtCtrCell()) {
      _vEdgeBRs[SOUTH_PDIR].nwCorner(iRowWorkEnd + 1 + minIRow,
                                     iColWorkBegin);
      _vEdgeBRs[SOUTH_PDIR].seCorner(iRowWorkEnd, iColWorkEnd);
      iRowSIntr = iRowWorkEnd + minIRow;
    }
    else {
      _vEdgeBRs[SOUTH_PDIR].nwCorner(iRowWorkEnd + 1 + minIRow - maxIRow,
                                     iColWorkBegin);
      _vEdgeBRs[SOUTH_PDIR].seCorner(iRowWorkEnd, iColWorkEnd);
      iRowSIntr = iRowWorkEnd + minIRow - maxIRow;
    }
  }
  else {
    iRowSIntr = iRowWorkEnd;
  }
  
  if(_pInfo->hasNbrs(WEST_DIR) && maxICol > 0) {
    const IntVect &vNbrIDs = _pInfo->nbrSubSpcIDs(WEST_DIR);
    for(int iNbr = 0; iNbr < vNbrIDs.size(); iNbr++) {
      int nbrSpcID = vNbrIDs[iNbr];
      const SubSpaceInfo *pNbrInfo = pvSubInfos->findSubInfo(nbrSpcID);
      if(!pNbrInfo) {
        cerr << __FILE__ << " " << __FUNCTION__ \
             << " Error: unable to find neighboring SubSpace [" \
             << nbrSpcID << "]" << endl;
        return false;
      }
      int iNbrRowBegin = pNbrInfo->iRowBegin();
      iNbrRowBegin = (iNbrRowBegin <= iRowBegin) ?
                      0 : (iNbrRowBegin - iRowBegin);
      _mSendBRs[WEST_DIR][iNbr].nwCorner(iNbrRowBegin,
                                         0);
      int iNbrRowEnd = pNbrInfo->iRowEnd();
      iNbrRowEnd = (iNbrRowEnd >= iRowEnd) ?
                    nRows - 1 : (iNbrRowEnd - iRowBegin);
      _mSendBRs[WEST_DIR][iNbr].seCorner(iNbrRowEnd,
                                         iColWorkBegin - 1 + maxICol);
    }
    if(pTransition->onlyUpdtCtrCell()) {
      _vEdgeBRs[WEST_PDIR].nwCorner(iRowNIntr, iColWorkBegin);
      _vEdgeBRs[WEST_PDIR].seCorner(iRowSIntr,
                                    iColWorkBegin - 1 + maxICol);
      iColWIntr = iColWorkBegin + maxICol;
    }
    else {
      _vEdgeBRs[WEST_PDIR].nwCorner(iRowNIntr, iColWorkBegin);
      _vEdgeBRs[WEST_PDIR].seCorner(iRowSIntr,
                                    iColWorkBegin - 1 + maxICol - minICol);
      iColWIntr = iColWorkBegin + maxICol - minICol;
    }
  }
  else {
    iColWIntr = iColWorkBegin;
  }
  
  if(_pInfo->hasNbrs(EAST_DIR) && minICol < 0) {
    const IntVect &vNbrIDs = _pInfo->nbrSubSpcIDs(EAST_DIR);
    for(int iNbr = 0; iNbr < vNbrIDs.size(); iNbr++) {
      int nbrSpcID = vNbrIDs[iNbr];
      const SubSpaceInfo *pNbrInfo = pvSubInfos->findSubInfo(nbrSpcID);
      if(!pNbrInfo) {
        cerr << __FILE__ << " " << __FUNCTION__ \
             << " Error: unable to find neighboring SubSpace [" \
             << nbrSpcID << "]" << endl;
        return false;
      }
      int iNbrRowBegin = pNbrInfo->iRowBegin();
      iNbrRowBegin = (iNbrRowBegin <= iRowBegin) ?
                      0 : (iNbrRowBegin - iRowBegin);
      _mSendBRs[EAST_DIR][iNbr].nwCorner(iNbrRowBegin,
                                         iColWorkEnd + 1 + minICol);
      int iNbrRowEnd = pNbrInfo->iRowEnd();
      iNbrRowEnd = (iNbrRowEnd >= iRowEnd) ?
                    nRows - 1 : (iNbrRowEnd - iRowBegin);
      _mSendBRs[EAST_DIR][iNbr].seCorner(iNbrRowEnd,
                                         nCols - 1);
    }
    if(pTransition->onlyUpdtCtrCell()) {
      _vEdgeBRs[EAST_PDIR].nwCorner(iRowNIntr,
                                    iColWorkEnd + 1 + minICol);
      _vEdgeBRs[EAST_PDIR].seCorner(iRowSIntr, iColWorkEnd);
      iColEIntr = iColWorkEnd + minICol;
    }
    else {
      _vEdgeBRs[EAST_PDIR].nwCorner(iRowNIntr,
                                    iColWorkEnd + 1 + minICol - maxICol);
      _vEdgeBRs[EAST_PDIR].seCorner(iRowSIntr, iColWorkEnd);
      iColEIntr = iColWorkEnd + minICol - maxICol;
    }
  }
  else {
    iColEIntr = iColWorkEnd;
  }
  
  _interiorBR.nwCorner(iRowNIntr, iColWIntr);
  _interiorBR.seCorner(iRowSIntr, iColEIntr);
  
  if(_pInfo->hasNbrs(NORTHEAST_DIR) &&
     pNbrhood->hasNbrs(SOUTHWEST_DIR)) {
    _mSendBRs[NORTHEAST_DIR][0].nwCorner(0,
                                         iColWorkEnd + minICol + 1);
    _mSendBRs[NORTHEAST_DIR][0].seCorner(iRowWorkBegin + maxIRow -1,
                                         nCols - 1);
  }
  if(_pInfo->hasNbrs(SOUTHEAST_DIR) &&
     pNbrhood->hasNbrs(NORTHWEST_DIR)) {
    _mSendBRs[SOUTHEAST_DIR][0].nwCorner(iRowWorkEnd + minIRow + 1,
                                         iColWorkEnd + minICol + 1);
    _mSendBRs[SOUTHEAST_DIR][0].seCorner(nRows - 1, nCols - 1);
  }
  if(_pInfo->hasNbrs(SOUTHWEST_DIR) &&
     pNbrhood->hasNbrs(NORTHEAST_DIR)) {
    _mSendBRs[SOUTHWEST_DIR][0].nwCorner(iRowWorkEnd + minIRow + 1,
                                         0);
    _mSendBRs[SOUTHWEST_DIR][0].seCorner(nRows - 1,
                                         iColWorkBegin + maxICol - 1);
  }
  if(_pInfo->hasNbrs(NORTHWEST_DIR) &&
     pNbrhood->hasNbrs(SOUTHEAST_DIR)) {
    _mSendBRs[NORTHWEST_DIR][0].nwCorner(0, 0);
    _mSendBRs[NORTHWEST_DIR][0].seCorner(iRowWorkBegin + maxIRow - 1,
                                         iColWorkBegin + maxICol - 1);
  }
  
  return true;
}


/****************************************
*             Public Methods            *
*****************************************/
template <class elemType>
inline pRPL::SubSpace<elemType>::
SubSpace()
  :CellSpace<elemType>(),
   _pInfo(0){}

template <class elemType>
inline pRPL::SubSpace<elemType>::
SubSpace(const pRPL::CellSpace<elemType> &cellSpace,
         const pRPL::SubSpaceInfo& subSpaceInfo)
  :_pInfo(&subSpaceInfo){
  if(cellSpace.empty()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to construct a SubSpace from"
         << " an empty CellSpace" << endl;
    exit(-1);
  }

  if(!(subSpaceInfo.glbDims() == cellSpace.dims())) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: the subSpaceInfo does NOT match the CellSpace" << endl;
    exit(-1);
  }

  SpaceDims dims = subSpaceInfo.dims();
  if(dims.valid()) {
    CellSpace<elemType>::_dims = dims;
    CellSpace<elemType>::_initMem(dims);
  }
  else {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid subSpaceInfo dimensions" << endl;
    exit(-1);
  }
  
  int iElemLcl;
  int iRowGlb = _pInfo->iRowBegin();
  int iColGlb = _pInfo->iColBegin();
  for(iElemLcl = 0; iElemLcl < CellSpace<elemType>::size(); iElemLcl++) {
    CellSpace<elemType>::_matrix[iElemLcl] = cellSpace[iRowGlb][iColGlb];
    iColGlb++;
    if(iColGlb > _pInfo->iColEnd()) {
      iRowGlb++;
      iColGlb = _pInfo->iColBegin();
    }
  } // End of iElemLcl loop 
}

template <class elemType>
inline pRPL::SubSpace<elemType>::
SubSpace(const SubSpaceInfo& subSpaceInfo)
  :CellSpace<elemType>(subSpaceInfo.dims()),
   _pInfo(&subSpaceInfo) {}

template <class elemType>
inline pRPL::SubSpace<elemType>::
SubSpace(const SubSpace<elemType> &rhs)
  :CellSpace<elemType>(CellSpace<elemType>(rhs)),
   _pInfo(rhs._pInfo),
   _mSendBRs(rhs._mSendBRs),
   _vEdgeBRs(rhs._vEdgeBRs),
   _interiorBR(rhs._interiorBR) {}

template <class elemType>
inline pRPL::SubSpace<elemType>::
~SubSpace() {}

template <class elemType>
inline pRPL::SubSpace<elemType>& pRPL::SubSpace<elemType>::
operator=(const SubSpace<elemType> &rhs) {
  if(this != &rhs) {
    CellSpace<elemType>(*this) 
      = CellSpace<elemType>::operator=(CellSpace<elemType>(rhs));
    _pInfo = rhs._pInfo;
    _mSendBRs = rhs._mSendBRs;
    _vEdgeBRs = rhs._vEdgeBRs;
    _interiorBR = rhs._interiorBR;
  }
  return *this;
}

template <class elemType>
inline bool pRPL::SubSpace<elemType>::
hasInfo() const {
  bool infoExists = true;
  if(!_pInfo) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: SubSpaceInfo is missing" << endl;
    infoExists = false;
  }
  return infoExists;
}

template <class elemType>
inline int pRPL::SubSpace<elemType>::
id() const {
  return _pInfo->id();
}

template <class elemType>
inline pRPL::DomDcmpType pRPL::SubSpace<elemType>::
domDcmpType() const {
  return _pInfo->domDcmpType();
}

template <class elemType>
inline const pRPL::SpaceDims& pRPL::SubSpace<elemType>::
glbDims() const {
  return _pInfo->glbDims();
}

template <class elemType>
const pRPL::CoordBR& pRPL::SubSpace<elemType>::
MBR() const {
  return _pInfo->MBR();
}

template <class elemType>
inline const pRPL::CoordBR& pRPL::SubSpace<elemType>::
workBR() const {
  return _pInfo->workBR();
}

template <class elemType>
inline double pRPL::SubSpace<elemType>::
sizeRatio() const {
  return _pInfo->sizeRatio();
}

template <class elemType>
const pRPL::CellCoord pRPL::SubSpace<elemType>::
lclCoord2glbCoord(const CellCoord& lclCoord) const {
  return _pInfo->lclCoord2glbCoord(lclCoord);
}
template <class elemType>
const pRPL::CellCoord pRPL::SubSpace<elemType>::
glbCoord2lclCoord(const CellCoord& glbCoord) const {
  return _pInfo->glbCoord2lclCoord(glbCoord);
}

template <class elemType>
inline int pRPL::SubSpace<elemType>::
lclCoord2glbIdx(int iRowLcl, int iColLcl) const {
  return _pInfo->lclCoord2glbIdx(iRowLcl,
                                 iColLcl);
}

template <class elemType>
inline int pRPL::SubSpace<elemType>::
lclCoord2glbIdx(const CellCoord &lclCoord) const {
  return _pInfo->lclCoord2glbIdx(lclCoord);
}

template <class elemType>
inline const pRPL::CellCoord pRPL::SubSpace<elemType>::
glbIdx2lclCoord(int glbIdx) const {
  return _pInfo->glbIdx2lclCoord(glbIdx);
}

template <class elemType>
inline int pRPL::SubSpace<elemType>::
glbIdx2lclIdx(int glbIdx) const {
  return _pInfo->glbIdx2lclIdx(glbIdx);
}

template <class elemType>
inline int pRPL::SubSpace<elemType>::
glbCoord2glbIdx(int iRowGlb, int iColGlb) const {
  return _pInfo->glbCoord2glbIdx(iRowGlb,
                                 iColGlb);
}

template <class elemType>
inline int pRPL::SubSpace<elemType>::
glbCoord2glbIdx(const CellCoord &glbCoord) const {
  return _pInfo->glbCoord2glbIdx(glbCoord);
}

template <class elemType>
inline bool pRPL::SubSpace<elemType>::
hasNbrs(int iDir) const {
  return _pInfo->hasNbrs(iDir);
}

template <class elemType>
inline int pRPL::SubSpace<elemType>::
nNbrs(int iDir) const {
  return _pInfo->nNbrs(iDir);
}

template <class elemType>
inline const pRPL::IntVect& pRPL::SubSpace<elemType>::
nbrSubSpcIDs(int iDir) const {
  return _pInfo->nbrSubSpcIDs(iDir);
}

template <class elemType>
inline int pRPL::SubSpace<elemType>::
nbrDir(int nbrID) const {
  return _pInfo->nbrDir(nbrID);
}

template <class elemType>
inline int pRPL::SubSpace<elemType>::
oppositeDir(int iDir) const {
  return _pInfo->oppositeDir(iDir);
}

template <class elemType>
inline pRPL::MeshDir pRPL::SubSpace<elemType>::
spcDir2MeshDir(int iDir) const {
  return _pInfo->spcDir2MeshDir(iDir);
}

template <class elemType>
inline int pRPL::SubSpace<elemType>::
nNbrDirs() const {
  return _pInfo->nNbrDirs();
}

template <class elemType>
inline int pRPL::SubSpace<elemType>::
nTotNbrs() const {
  return _pInfo->nTotNbrs();
}

template <class elemType>
bool pRPL::SubSpace<elemType>::
makeStream2Send() {
  if(CellSpace<elemType>::empty() ||
     !hasInfo()) {
    return false;
  }

  int elemSize = sizeof(elemType);
  int intSize = sizeof(int);

  int numNbrDirs = nNbrDirs();
  _mStream2Send.resize(numNbrDirs);
  for(int iDir = 0; iDir < numNbrDirs; iDir++) {
    int numNbrs = nNbrs(iDir);
    _mStream2Send[iDir].resize(numNbrs);
    for(int iNbr = 0; iNbr < numNbrs; iNbr++) {
      const CoordBR &sendBR = _mSendBRs[iDir][iNbr];
      if(!sendBR.valid(CellSpace<elemType>::_dims)) {
        continue;
      }

      typename map<elemType, IntVect>::const_iterator iCellGroup = 
        CellSpace<elemType>::_mUpdtCells.begin();
      while(iCellGroup != CellSpace<elemType>::_mUpdtCells.end()) {
        const elemType &newVal = (*iCellGroup).first;
        const IntVect &vIdxs = (*iCellGroup).second;
        IntVect vIdxs2Send;
        for(int iIdx = 0; iIdx < vIdxs.size(); iIdx++) {
          CellCoord coord(vIdxs[iIdx],
                          CellSpace<elemType>::_dims);
          if(sendBR.contain(coord)) {
            vIdxs2Send.push_back(lclCoord2glbIdx(coord));
          }
        }
        int nIdxs2Send = vIdxs2Send.size();
        if(nIdxs2Send > 0) {
          int streamSize = _mStream2Send[iDir][iNbr].size();
          _mStream2Send[iDir][iNbr].resize(streamSize
                                           + elemSize
                                           + (nIdxs2Send+1)*intSize);
          memcpy(&(_mStream2Send[iDir][iNbr][streamSize]),
                 &newVal, elemSize);
          memcpy(&(_mStream2Send[iDir][iNbr][streamSize+elemSize]),
                 &nIdxs2Send, intSize);
          memcpy(&(_mStream2Send[iDir][iNbr][streamSize+elemSize+intSize]),
                 &(vIdxs2Send[0]), intSize*nIdxs2Send);
        }
        iCellGroup++;
      } // End of iCellGroup loop
    } // End of iNbr loop
  } // End of iDir loop
  return true;
}

template <class elemType>
bool pRPL::SubSpace<elemType>::
loadStreamRecved(Transition<elemType> *pTransition,
                 Neighborhood<elemType> *pNbrhood) {
  if(CellSpace<elemType>::empty() ||
     !hasInfo()) {
    return false;
  }
  if(pTransition) {
    if(!pTransition->cellSpace(this) ||
       !pTransition->nbrhood(pNbrhood)) {
      return false;
    }
  }

  int elemSize = sizeof(elemType);
  int intSize = sizeof(int);

  int iChar = 0;
  while(iChar < _vStreamRecved.size()) {
    elemType newVal;
    memcpy(&newVal, &(_vStreamRecved[iChar]), elemSize);
    iChar += elemSize;

    int nIdxs;
    memcpy(&nIdxs, &(_vStreamRecved[iChar]), intSize);
    iChar += intSize;

    for(int iIdx = 0; iIdx < nIdxs; iIdx++) {
      int glbIdx;
      memcpy(&glbIdx, &(_vStreamRecved[iChar]), intSize);
      iChar += intSize;

      int lclIdx = glbIdx2lclIdx(glbIdx);
      if(!CellSpace<elemType>::validIdx(lclIdx)) {
        cerr << __FILE__ " " << __FUNCTION__ \
             << " Error: unable to load global idx[" << glbIdx \
             << "] on SubSpace[" << id() << "]" << endl;
        return false;
      }
      if(pTransition) {
        CellCoord lclCoord(lclIdx, CellSpace<elemType>::_dims);
        if(!pTransition->finalize(newVal, lclCoord)) {
          return false;
        }
      }
      else {
        CellSpace<elemType>::_matrix[lclIdx] = newVal;
      }
    }
  }

  return true;
}

template <class elemType>
void pRPL::SubSpace<elemType>::
cleanStreamBuf() {
  _mStream2Send.clear();
  _vStreamRecved.clear();
}

template <class elemType>
pRPL::CharVect& pRPL::SubSpace<elemType>::
stream2Send(int iDir, int iNbr) {
  return _mStream2Send.at(iDir).at(iNbr);
}

template <class elemType>
pRPL::CharVect& pRPL::SubSpace<elemType>::
streamRecved() {
  return _vStreamRecved;
}

template <class elemType>
inline int pRPL::SubSpace<elemType>::
recvBufSize() const {
  return _vStreamRecved.size();
}

template <class elemType>
bool pRPL::SubSpace<elemType>::
recvBufSize(int nChars) {
  if(nChars <= 0) {
    cerr << __FILE__ << __FUNCTION__ \
         << " Error: unable to re-allocate the receiving buffer" 
         <<" for SubSpace[" << id() \
         << "] with a non-positive number (" << nChars \
         << ")" << endl;
    return false;
  }

  _vStreamRecved.resize(nChars);

  return true;
}

template <class elemType>
bool pRPL::SubSpace<elemType>::
calcBRs(const Transition<elemType> *pTransition,
        const Neighborhood<elemType> *pNbrhood,
        const SubInfoVect *pvSubInfos) {
  if(CellSpace<elemType>::empty() ||
     !hasInfo()) {
    return false;
  }
  if(!pTransition ||
     !pNbrhood ||
     !pvSubInfos) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to calcuate Bounding Rectangles" \
         << " using a NULL pointer to Transition (" << pTransition << "), or" \
         << " using a NULL pointer to Neighborhood (" << pNbrhood << "), or" \
         << " using a NULL pointer to SubInfoVect ("  << pvSubInfos << ")" << endl;
    return false;
  }
  
  _initBRs();
  bool done = true;
  if(pTransition->needExchange()) {
    switch(_pInfo->domDcmpType()) {
      case NON_DCMP:
        break;
      case ROWWISE_DCMP:
        _rowwiseBRs(pTransition, pNbrhood);
        break;
      case COLWISE_DCMP:
        _colwiseBRs(pTransition, pNbrhood);
        break;
      case BLOCK_DCMP:
        done = _blockwiseBRs(pTransition, pNbrhood, pvSubInfos);
        break;
      default:
        cerr << __FILE__ << __FUNCTION__ \
             << " Error: invalid decomposition type (" \
             << _pInfo->domDcmpType() << ")" << endl;
        done = false;
        break;
    }
  }
  
  return done;
}

template <class elemType>
bool pRPL::SubSpace<elemType>::
updateEdges(Transition<elemType> *pTransition,
            Neighborhood<elemType> *pNbrhood) {
  if(!pTransition) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to update SubSpace[" << id()\
         << "] using a NULL pointer to Transition" \
         << endl;
    return false;
  }
  if(!pTransition->onlyUpdtCtrCell() &&
     (!pNbrhood || pNbrhood->empty())) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to update SubSpace[" << id()\
         << "] using a NULL pointer to Neighborhood or" \
         << " an empty Neighborhood while intending to" \
         << " update Neighboring cells" \
         << endl;
    return false;
  }
  if(CellSpace<elemType>::empty()||
     !hasInfo()) {
    return false;
  }
  if(!pTransition->cellSpace(this) ||
     !pTransition->nbrhood(pNbrhood)) {
    return false;
  }
  
  vector<pair<int, elemType> > vUpdtedCells;
  for(int iEdge = 0; iEdge < _pInfo->nEdges(); iEdge++) {
    const CoordBR &edgeBR = _vEdgeBRs[iEdge];
    if(!edgeBR.valid(CellSpace<elemType>::_dims)) {
      continue;
    }
    for(int iRow = edgeBR.minIRow(); iRow <= edgeBR.maxIRow(); iRow++) {
      for(int iCol = edgeBR.minICol(); iCol <= edgeBR.maxICol(); iCol++) {
        CellCoord coord(iRow, iCol);
        vUpdtedCells.clear();
        if(!pTransition->evaluate(vUpdtedCells, coord)) {
          return false;
        }
        if(!vUpdtedCells.empty()) {
          for(int iCell = 0; iCell < vUpdtedCells.size(); iCell++) {
            int &iElem = vUpdtedCells[iCell].first;
            if(!CellSpace<elemType>::validIdx(iElem)) {
              cerr << __FILE__ << " " << __FUNCTION__ \
                   << " Error: unable to update coord[" \
                   << coord << "] on SubSpace[" << id() \
                   << "]" << endl;
              return false;
            }
            elemType &val = vUpdtedCells[iCell].second;
            if(pTransition->needExchange() ||
               pTransition->needFinalize()) {
              CellSpace<elemType>::_mUpdtCells[val].push_back(iElem);
            }
            if(!pTransition->needFinalize()) {
              CellSpace<elemType>::_matrix[iElem] = val;
            }
          }
        }
      } // End of iCol loop
    } // End of iRow loop
  } // End of iEdge loop

  if(pTransition->needExchange()) {
    makeStream2Send();
  }
  if(!pTransition->needFinalize() &&
     pTransition->needExchange()) {
    CellSpace<elemType>::_mUpdtCells.clear();
  }
  return true;
}

template <class elemType>
bool pRPL::SubSpace<elemType>::
updateInterior(Transition<elemType> *pTransition,
               Neighborhood<elemType> *pNbrhood) {
  if(!pTransition) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to update SubSpace[" << id()\
         << "] using a NULL pointer to Transition" \
         << endl;
    return false;
  }
  if(!pTransition->onlyUpdtCtrCell() &&
     (!pNbrhood || pNbrhood->empty())) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to update SubSpace[" << id()\
         << "] using a NULL pointer to Neighborhood or" \
         << " an empty Neighborhood while intending to" \
         << " update Neighboring cells" \
         << endl;
    return false;
  }
  if(CellSpace<elemType>::empty() ||
     !hasInfo()) {
    return false;
  }
  if(!pTransition->cellSpace(this) ||
     !pTransition->nbrhood(pNbrhood)) {
    return false;
  }

  if(!_interiorBR.valid(CellSpace<elemType>::_dims)) {
    return true;
  }

  vector<pair<int, elemType> > vUpdtedCells;
  for(int iRow = _interiorBR.minIRow(); iRow <= _interiorBR.maxIRow(); iRow++) {
    for(int iCol = _interiorBR.minICol(); iCol <= _interiorBR.maxICol(); iCol++) {
      CellCoord coord(iRow, iCol);
      vUpdtedCells.clear();
      if(!pTransition->evaluate(vUpdtedCells, coord)) {
        return false;
      }
      if(!vUpdtedCells.empty()) {
        for(int iCell = 0; iCell < vUpdtedCells.size(); iCell++) {
          int &iElem = vUpdtedCells[iCell].first;
          if(!CellSpace<elemType>::validIdx(iElem)) {
            cerr << __FILE__ << " " << __FUNCTION__ \
                 << " Error: unable to update coord[" \
                 << coord << "] on SubSpace[" << id() \
                 << "]" << endl;
            return false;
          }
          elemType &val = vUpdtedCells[iCell].second;
          if(pTransition->needFinalize()) {
            CellSpace<elemType>::_mUpdtCells[val].push_back(iElem);
          }
          else {
            CellSpace<elemType>::_matrix[iElem] = val;
          }
        }
      }
    }
  }

  return true;
}

template <class elemType>
bool pRPL::SubSpace<elemType>::
updateAll(Transition<elemType> *pTransition,
          Neighborhood<elemType> *pNbrhood) {
  if(!pTransition) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to update SubSpace[" << id()\
         << "] using a NULL pointer to Transition" \
         << endl;
    return false;
  }
  if(CellSpace<elemType>::empty() ||
     !hasInfo()) {
    return false;
  }
  const CoordBR &workBR = _pInfo->workBR();
  if(pTransition->needExchange()) {
    bool ndFnlz = pTransition->needFinalize();
    pTransition->needFinalize(true);
    if(!CellSpace<elemType>::update(pTransition, pNbrhood,
                                    &workBR) ||
       !makeStream2Send()) {
      return false;
    }
    pTransition->needFinalize(ndFnlz);
    if(!pTransition->needFinalize()) {
      if(!CellSpace<elemType>::updateFinalize()) {
        return false;
      }
    }
  }
  else {
    if(!CellSpace<elemType>::update(pTransition, pNbrhood,
                                    &workBR)) {
      return false;
    }
  }
  return true;
}

template <class elemType>
bool pRPL::SubSpace<elemType>::
updateEdges(IntVect &vIdxs2Eval,
            Transition<elemType> *pTransition,
            Neighborhood<elemType> *pNbrhood) {
  if(!pTransition) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to update SubSpace[" << id()\
         << "] using a void pointer to Transition" \
         << endl;
    return false;
  }
  if(!pTransition->onlyUpdtCtrCell() &&
     (!pNbrhood || pNbrhood->empty())) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to update SubSpace[" << id()\
         << "] using a void pointer to Neighborhood or" \
         << " an empty Neighborhood while intending to" \
         << " update Neighboring cells" \
         << endl;
    return false;
  }
  if(CellSpace<elemType>::empty()||
     !hasInfo()) {
    return false;
  }
  if(!pTransition->cellSpace(this) ||
     !pTransition->nbrhood(pNbrhood)) {
    return false;
  }
  
  vector<pair<int, elemType> > vUpdtedCells;
  for(int iEdge = 0; iEdge < _pInfo->nEdges(); iEdge++) {
    if(!vIdxs2Eval.empty()) {
      const CoordBR &edgeBR = _vEdgeBRs[iEdge];
      if(!edgeBR.valid(CellSpace<elemType>::_dims)) {
        continue;
      }
      IntVctItr iIdx = vIdxs2Eval.begin();
      while(iIdx != vIdxs2Eval.end()) {
        int iElem = *iIdx;
        CellCoord coord(iElem, CellSpace<elemType>::_dims);
        if(edgeBR.contain(coord)) {
          vUpdtedCells.clear();
          if(!pTransition->evaluate(vUpdtedCells, coord)) {
            return false;
          }
          if(!vUpdtedCells.empty()) {
            for(int iCell = 0; iCell < vUpdtedCells.size(); iCell++) {
              int &iElem = vUpdtedCells[iCell].first;
              if(!CellSpace<elemType>::validIdx(iElem)) {
                cerr << __FILE__ << " " << __FUNCTION__ \
                     << " Error: unable to update coord[" \
                     << coord << "] on SubSpace[" << id() \
                     << "]" << endl;
                return false;
              }
              elemType &val = vUpdtedCells[iCell].second;
              if(pTransition->needExchange() ||
                 pTransition->needFinalize()) {
                CellSpace<elemType>::_mUpdtCells[val].push_back(iElem);
              }
              if(!pTransition->needFinalize()) {
                CellSpace<elemType>::_matrix[iElem] = val;
              }
            }
          }
          vIdxs2Eval.erase(iIdx);
        }
        else {
          iIdx++;
        }
      } // End of iIdx loop
    } // if vIdxs2Eval is not empty
  } // End of iEdge loop

  if(pTransition->needExchange()) {
    makeStream2Send();
  }
  if(!pTransition->needFinalize() &&
     pTransition->needExchange()) {
    CellSpace<elemType>::_mUpdtCells.clear();
  }
  return true;
}

template <class elemType>
bool pRPL::SubSpace<elemType>::
updateInterior(IntVect &vIdxs2Eval,
               Transition<elemType> *pTransition,
               Neighborhood<elemType> *pNbrhood) {
  if(!pTransition) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to update SubSpace[" << id()\
         << "] using a void pointer to Transition" \
         << endl;
    return false;
  }
  if(!pTransition->onlyUpdtCtrCell() &&
     (!pNbrhood || pNbrhood->empty())) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to update SubSpace[" << id()\
         << "] using a void pointer to Neighborhood or" \
         << " an empty Neighborhood while intending to" \
         << " update Neighboring cells" \
         << endl;
    return false;
  }
  if(CellSpace<elemType>::empty() ||
     !hasInfo()) {
    return false;
  }
  if(!pTransition->cellSpace(this) ||
     !pTransition->nbrhood(pNbrhood)) {
    return false;
  }

  if(!vIdxs2Eval.empty() &&
     _interiorBR.valid(CellSpace<elemType>::_dims)) {
    vector<pair<int, elemType> > vUpdtedCells;
    IntVctItr iIdx = vIdxs2Eval.begin();
    while(iIdx != vIdxs2Eval.end()) {
      int iElem = *iIdx;
      CellCoord coord(iElem, CellSpace<elemType>::_dims);
      if(_interiorBR.contain(coord)) {
        vUpdtedCells.clear();
        if(!pTransition->evaluate(vUpdtedCells, coord)) {
          return false;
        }
        if(!vUpdtedCells.empty()) {
          for(int iCell = 0; iCell < vUpdtedCells.size(); iCell++) {
            int &iElem = vUpdtedCells[iCell].first;
            if(!CellSpace<elemType>::validIdx(iElem)) {
              cerr << __FILE__ << " " << __FUNCTION__ \
                   << " Error: unable to update coord[" \
                   << coord << "] on SubSpace[" << id() \
                   << "]" << endl;
              return false;
            }
            elemType &val = vUpdtedCells[iCell].second;
            if(pTransition->needFinalize()) {
              CellSpace<elemType>::_mUpdtCells[val].push_back(iElem);
            }
            else {
              CellSpace<elemType>::_matrix[iElem] = val;
            }
          }
        }
        vIdxs2Eval.erase(iIdx);
      }
      else {
        iIdx++;
      }
    } // End of iIdx loop
  }
  return true;
}

template <class elemType>
bool pRPL::SubSpace<elemType>::
updateAll(IntVect &vIdxs2Eval,
          Transition<elemType> *pTransition,
          Neighborhood<elemType> *pNbrhood) {
  if(!pTransition) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to update SubSpace[" << id()\
         << "] using a void pointer to Transition" \
         << endl;
    return false;
  }
  if(CellSpace<elemType>::empty() ||
     !hasInfo()) {
    return false;
  }

  const CoordBR &workBR = _pInfo->workBR();
  if(pTransition->needExchange()) {
    bool ndFnlz = pTransition->needFinalize();
    pTransition->needFinalize(true);
    if(!CellSpace<elemType>::update(vIdxs2Eval,
                                    pTransition, pNbrhood,
                                    &workBR) ||
       !makeStream2Send()) {
      return false;
    }
    pTransition->needFinalize(ndFnlz);
    if(!pTransition->needFinalize()) {
      if(!CellSpace<elemType>::updateFinalize()) {
        return false;
      }
    }
  }
  else {
    if(!CellSpace<elemType>::update(vIdxs2Eval,
                                    pTransition, pNbrhood,
                                    &workBR)) {
      return false;
    }
  }
  return true;
}

#endif
