#ifndef SMPLDCMP_H
#define SMPLDCMP_H

/***************************************************************************
* smplDcmp.h
*
* Project: pRPL, v 1.0
* Purpose: Header file for class pRPL::SmplDcmp
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
#include "cellSpace.h"
#include "neighborhood.h"
#include "subSpaceInfo.h"
#include "dcmpUtilities.h"
#include "subInfoVect.h"

namespace pRPL {
  template <class elemType>
  class SmplDcmp {
    public:
      SmplDcmp()
        :_pNbrhood(0) {}

      SmplDcmp(const SpaceDims &cellSpaceDims,
               const Neighborhood<elemType> &nbrhood);
      ~SmplDcmp() {}

      bool rowDcmp(SubInfoVect &vpSubSpcInfos,
                   int nSubSpcs) const;
      bool colDcmp(SubInfoVect &vpSubSpcInfos,
                   int nSubSpcs) const;
      bool blockDcmp(SubInfoVect &vpSubSpcInfos,
                     int nRowSubSpcs,
                     int nColSubSpcs) const;

    private:
      void _smpl1DDcmp(int &subBegin, int &subEnd,
                       int glbBegin, int glbEnd,
                       int nSubSpcs, int iSubSpc) const;

    private:
      SpaceDims _glbDims;
      CoordBR _glbWorkBR;
      const Neighborhood<elemType> *_pNbrhood;
  };
};

/****************************************
*             Private Methods           *
*****************************************/

template <class elemType>
void pRPL::SmplDcmp<elemType>::
_smpl1DDcmp(int &subBegin, int &subEnd,
            int glbBegin, int glbEnd,
            int nSubSpcs, int iSubSpc) const {
  int glbRange = glbEnd - glbBegin + 1;
  int remainder = glbRange % nSubSpcs;
  int blockSize;
  if(remainder == 0) {
    blockSize = glbRange / nSubSpcs;
    subBegin = glbBegin + iSubSpc * blockSize;
    subEnd = subBegin + blockSize - 1;
  }
  else {
    blockSize = glbRange / nSubSpcs + 1;
    if(iSubSpc < remainder) {
      subBegin = glbBegin + iSubSpc * blockSize;
      subEnd = subBegin + blockSize - 1;
    }
    else {
      if(iSubSpc == remainder) {
        subBegin = glbBegin+ iSubSpc * blockSize;
        subEnd = subBegin + blockSize - 2;	
      }
      else {
        subBegin = glbBegin 
                + remainder * blockSize 
                + (iSubSpc - remainder)*(blockSize-1);
        subEnd = subBegin + blockSize - 2;
      }
    }
  }
}

/****************************************
*             Public Methods            *
*****************************************/

template <class elemType>
pRPL::SmplDcmp<elemType>::
SmplDcmp(const SpaceDims &cellSpaceDims,
         const Neighborhood<elemType> &nbrhood) {
  if(!checkDcmpParms(cellSpaceDims, nbrhood)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid parameters to construct a SmplDcmp" \
         << endl;
    exit(-1);
  }
  if(!nbrhood.calcWorkBR(_glbWorkBR, cellSpaceDims)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid Neighborhood to construct a SmplDcmp" \
         << endl;
    exit(-1);
  }
  _glbDims = cellSpaceDims;
  _pNbrhood = &nbrhood;
}

template <class elemType>
bool pRPL::SmplDcmp<elemType>::
rowDcmp(SubInfoVect &vpSubSpcInfos,
        int nSubSpcs) const {
  vpSubSpcInfos.clear();
  if(nSubSpcs < 1 ||
     nSubSpcs > _glbWorkBR.nRows()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid number of SubSpaces (" << nSubSpcs << ")" \
         << " considering the global working bounding rectangle (" << _glbWorkBR << ")" \
         << endl;
    return false;
  }

  DomDcmpType dcmpType = ROWWISE_DCMP;
  int glbBegin = _glbWorkBR.nwCorner().iRow();
  int glbEnd = _glbWorkBR.seCorner().iRow();

  for(int iSubSpc = 0; iSubSpc < nSubSpcs; iSubSpc++) {
    int subBegin, subEnd;
    _smpl1DDcmp(subBegin, subEnd,
                glbBegin, glbEnd,
                nSubSpcs, iSubSpc);
    CellCoord nwCorner(subBegin + _pNbrhood->minIRow(),
                       0);
    CellCoord seCorner(subEnd + _pNbrhood->maxIRow(),
                       _glbDims.nCols() - 1);
    CoordBR subMBR(nwCorner, seCorner);
    SpaceDims dims(subMBR.nRows(), subMBR.nCols());
    CoordBR workBR;
    if(!_pNbrhood->calcWorkBR(workBR, dims)) {
      return false;
    }

    vector<IntVect> mNbrSpcIDs(2);
    if(iSubSpc != 0) {
      mNbrSpcIDs[UPPER_DIR].push_back(iSubSpc - 1);
    }
    if(iSubSpc != nSubSpcs-1) {
      mNbrSpcIDs[LOWER_DIR].push_back(iSubSpc + 1);
    }

    int workload = workBR.size();
    SubInfoVect::iterator iInsert = 
      upper_bound(vpSubSpcInfos.begin(),
                  vpSubSpcInfos.end(),
                  InfoWLPair(0, workload),
                  greater<InfoWLPair>());
    vpSubSpcInfos.insert(iInsert,
      InfoWLPair(new SubSpaceInfo(iSubSpc, dcmpType, _glbDims,
                                  subMBR, workBR, mNbrSpcIDs),
                 workload));
  }
  return true;
}

template <class elemType>
bool pRPL::SmplDcmp<elemType>::
colDcmp(SubInfoVect &vpSubSpcInfos,
        int nSubSpcs) const {
  vpSubSpcInfos.clear();
  if(nSubSpcs < 1 ||
     nSubSpcs > _glbWorkBR.nCols()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid number of SubSpaces (" << nSubSpcs << ")" \
         << " considering the global working bounding rectangle (" << _glbWorkBR << ")" \
         << endl;
    return false;
  }

  DomDcmpType dcmpType = COLWISE_DCMP;
  int glbBegin = _glbWorkBR.nwCorner().iCol();
  int glbEnd = _glbWorkBR.seCorner().iCol();

  for(int iSubSpc = 0; iSubSpc < nSubSpcs; iSubSpc++) {
    int subBegin, subEnd;
    _smpl1DDcmp(subBegin, subEnd,
                glbBegin, glbEnd,
                nSubSpcs, iSubSpc);
    CellCoord nwCorner(0, subBegin + _pNbrhood->minICol());
    CellCoord seCorner(_glbDims.nRows() - 1, 
                       subEnd + _pNbrhood->maxICol());
    CoordBR subMBR(nwCorner, seCorner);
    SpaceDims dims(subMBR.nRows(), subMBR.nCols());
    CoordBR workBR;
    if(!_pNbrhood->calcWorkBR(workBR, dims)) {
      return false; 
    }

    vector<IntVect> mNbrSpcIDs(2);
    if(iSubSpc != 0) {
      mNbrSpcIDs[LEFT_DIR].push_back(iSubSpc - 1);
    }
    if(iSubSpc != nSubSpcs-1) {
      mNbrSpcIDs[RIGHT_DIR].push_back(iSubSpc + 1);
    }

    int workload = workBR.size();
    SubInfoVect::iterator iInsert = 
      upper_bound(vpSubSpcInfos.begin(),
                  vpSubSpcInfos.end(),
                  InfoWLPair(0, workload),
                  greater<InfoWLPair>());
    vpSubSpcInfos.insert(iInsert,
      InfoWLPair(new SubSpaceInfo(iSubSpc, dcmpType, _glbDims,
                                  subMBR, workBR, mNbrSpcIDs),
                 workload));
  }
  return true;
}

template <class elemType>
bool pRPL::SmplDcmp<elemType>::
blockDcmp(SubInfoVect &vpSubSpcInfos,
          int nRowSubSpcs,
          int nColSubSpcs) const {
  vpSubSpcInfos.clear();
  int nSubSpcs = nRowSubSpcs * nColSubSpcs;
  if(nSubSpcs < 1 ||
     nRowSubSpcs > _glbWorkBR.nRows() ||
     nColSubSpcs > _glbWorkBR.nCols()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid number of subspaces" \
         << endl;
    return false;
  }

  DomDcmpType dcmpType = BLOCK_DCMP;
  int glbRowBegin = _glbWorkBR.nwCorner().iRow();
  int glbRowEnd = _glbWorkBR.seCorner().iRow();
  int glbColBegin = _glbWorkBR.nwCorner().iCol();
  int glbColEnd = _glbWorkBR.seCorner().iCol();
  
  for(int iRowSubSpc = 0; iRowSubSpc < nRowSubSpcs; iRowSubSpc++) {
    int subRowBegin, subRowEnd;
    _smpl1DDcmp(subRowBegin, subRowEnd,
                glbRowBegin, glbRowEnd,
                nRowSubSpcs, iRowSubSpc);
    for(int iColSubSpc = 0; iColSubSpc < nColSubSpcs; iColSubSpc++) {
      int subColBegin, subColEnd;
      _smpl1DDcmp(subColBegin, subColEnd,
                  glbColBegin, glbColEnd,
                  nColSubSpcs, iColSubSpc);

      CellCoord nwCorner(subRowBegin + _pNbrhood->minIRow(), 
                         subColBegin + _pNbrhood->minICol());
      CellCoord seCorner(subRowEnd + _pNbrhood->maxIRow(), 
                         subColEnd + _pNbrhood->maxICol());

      CoordBR subMBR(nwCorner, seCorner);
      SpaceDims dims(subMBR.nRows(), subMBR.nCols());
      CoordBR workBR;
      if(!_pNbrhood->calcWorkBR(workBR, dims)) {
        return false; 
      }

      vector<IntVect> mNbrSpcIDs(8);
      for(int iDir = 0; iDir < 8; iDir++) {
        int iRowNbrSpc = iRowSubSpc + *(nbrLocations + 2*iDir);
        int iColNbrSpc = iColSubSpc + *(nbrLocations + 2*iDir + 1);
        if(iRowNbrSpc >= 0 && iRowNbrSpc < nRowSubSpcs &&
           iColNbrSpc >= 0 && iColNbrSpc < nColSubSpcs) {
          int nbrSpcID = iRowNbrSpc * nColSubSpcs + iColNbrSpc;
          mNbrSpcIDs[iDir].push_back(nbrSpcID);
        }
      }

      int iSubSpc = iRowSubSpc * nColSubSpcs + iColSubSpc;
      int workload = workBR.size();
      SubInfoVect::iterator iInsert = 
        upper_bound(vpSubSpcInfos.begin(),
                    vpSubSpcInfos.end(),
                    InfoWLPair(0, workload),
                    greater<InfoWLPair>());
      vpSubSpcInfos.insert(iInsert,
        InfoWLPair(new SubSpaceInfo(iSubSpc, dcmpType, _glbDims,
                                    subMBR, workBR, mNbrSpcIDs),
                   workload));
    } // End of iColSubSpc loop
  } // End of iRowSubSpc loop
  return true;
}

#endif
