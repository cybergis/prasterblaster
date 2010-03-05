#ifndef SUBSPACEINFO_H
#define SUBSPACEINFO_H

/***************************************************************************
* subSpaceInfo.h
*
* Project: pRPL, v 1.0
* Purpose: Header file for class pRPL::SubSpaceInfo
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
  class SubSpaceInfo {
    public:
      SubSpaceInfo();
      SubSpaceInfo(int id,
                   DomDcmpType domDcmpType,
                   const SpaceDims &glbDims,
                   const CoordBR &MBR,
                   const CoordBR &workBR,
                   const vector<IntVect> &mNbrSpcIDs);
      SubSpaceInfo(const SubSpaceInfo &rhs);
      SubSpaceInfo(const IntVect &vInfoPack,
                   IntVctItr &iVal);

      ~SubSpaceInfo() {}

      SubSpaceInfo& operator=(const SubSpaceInfo &rhs);
      bool operator==(const SubSpaceInfo &rhs) const;
      bool operator!=(const SubSpaceInfo &rhs) const;

      int id() const;
      DomDcmpType domDcmpType() const;
      int nGlbRows() const;
      int nGlbCols() const;
      const SpaceDims& glbDims() const;
      int nRows() const;
      int nCols() const;
      const SpaceDims& dims() const;
      int iRowBegin() const;
      int iColBegin() const;
      int iRowEnd() const;
      int iColEnd() const;
      const CoordBR& MBR() const;
      int iRowWorkBegin() const;
      int iColWorkBegin() const;
      int iRowWorkEnd() const;
      int iColWorkEnd() const;
      const CoordBR& workBR() const;
      double sizeRatio() const;
      
      bool validGlbIdx(int glbIdx) const;
      bool validGlbCoord(const CellCoord &glbCoord) const;
      bool validGlbCoord(int iRowGlb, int iColGlb) const;
      const CellCoord glbIdx2glbCoord(int glbIdx) const;
      int glbCoord2glbIdx(const CellCoord &glbCoord) const;
      int glbCoord2glbIdx(int iRowGlb, int iColGlb) const;

      bool validIdx(int idx) const;
      bool validCoord(const CellCoord &coord) const;
      bool validCoord(int iRow, int iCol) const;
      int coord2idx(const CellCoord &coord) const;
      int coord2idx(int iRow, int iCol) const;
      const CellCoord idx2coord(int idx) const;

      const CellCoord lclCoord2glbCoord(const CellCoord& lclCoord) const;
      const CellCoord lclCoord2glbCoord(int iRowLcl, int iColLcl) const;
      const CellCoord glbCoord2lclCoord(const CellCoord& glbCoord) const;
      const CellCoord glbCoord2lclCoord(int iRowGlb, int iColGlb) const;
      int lclCoord2glbIdx(const CellCoord &lclCoord) const;
      int lclCoord2glbIdx(int iRowLcl, int iColLcl) const;
      const CellCoord glbIdx2lclCoord(int glbIdx) const;
      int glbIdx2lclIdx(int glbIdx) const;
      int lclIdx2glbIdx(int lclIdx) const;

      int nNbrDirs() const;
      int nEdges() const;
      int nTotNbrs() const;
      bool validSpcDir(int iDir) const;
      bool hasNbrs(int iDir) const;
      int nNbrs(int iDir) const;
      const IntVect& nbrSubSpcIDs(int iDir) const;
      int nbrDir(int nbrID) const;
      MeshDir spcDir2MeshDir(int iDir) const;
      int oppositeDir(int iDir) const;

      void add2IntVect(IntVect &vInfoPack) const;
      bool fromIntVect(const IntVect &vInfoPack,
                       IntVctItr &iVal);

    protected:
      int _id;
      DomDcmpType _domDcmpType;
      SpaceDims _glbDims;
      CoordBR _MBR; /* MBR is in global coordinates */
      SpaceDims _dims;
      CoordBR _workBR; /* workBR is in local coordinates */
      vector<IntVect> _mNbrSpcIDs;
  };
};

inline pRPL::SubSpaceInfo::
SubSpaceInfo()
  :_id(ERROR_SPCID),
   _domDcmpType(NON_DCMP) {}

inline pRPL::SubSpaceInfo::
SubSpaceInfo(int id,
             DomDcmpType domDcmpType,
             const SpaceDims &glbDims,
             const CoordBR &MBR,
             const CoordBR &workBR,
             const vector<IntVect> &mNbrSpcIDs)
 :_domDcmpType(domDcmpType) {
  if(id < 0) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid SubSpace ID (" \
         << id << ")" << endl;
    exit(-1);
  }
  _id = id;

  if(!glbDims.valid()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: SubSpaceInfo[" << _id \
         << "] invalid glbDims (" << glbDims 
         << ")" << endl;
    exit(-1);
  }
  _glbDims = glbDims;

  if(!MBR.valid(glbDims)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: SubSpaceInfo[" << _id \
         << "] invalid MBR (" \
         << MBR << ") with the global dimensions (" \
         << glbDims << ")" << endl;
    exit(-1);
  }
  _MBR = MBR;
  _dims = SpaceDims(MBR.nRows(), MBR.nCols());

  if(!workBR.valid(_dims)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: SubSpaceInfo[" << _id \
         << "] invalid working range (" \
         << workBR << ") within the local dimensions (" \
         << _dims << ")" << endl;
    exit(-1);
  }
  _workBR = workBR;
  
  int numNbrDirs = nNbrDirs();
  if(numNbrDirs != mNbrSpcIDs.size()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: SubSpaceInfo[" << _id \
         << "] invalid number of neighbor ID directions (" \
         << mNbrSpcIDs.size() << ")" << endl;
    exit(-1);
  }
  _mNbrSpcIDs = mNbrSpcIDs;
}

inline pRPL::SubSpaceInfo::
SubSpaceInfo(const SubSpaceInfo &rhs)
  :_id(rhs._id),
   _domDcmpType(rhs._domDcmpType),
   _glbDims(rhs._glbDims),
   _dims(rhs._dims),
   _MBR(rhs._MBR),
   _workBR(rhs._workBR),
   _mNbrSpcIDs(rhs._mNbrSpcIDs) {}

inline pRPL::SubSpaceInfo::
SubSpaceInfo(const IntVect &vInfoPack,
             IntVctItr &iVal) {
  if(!fromIntVect(vInfoPack, iVal)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to construct a SubSpaceInfo" << endl;
    exit(-1);
  }
  /*
  cout << "spc[" << _id << "]" << endl \
       << "domDcmpType " << _domDcmpType << endl \
       << "glbDims " << _glbDims << endl \
       << "dims " << _dims << endl \
       << "MBR " << _MBR << endl \
       << "workBR " << _workBR << endl;
  for(int iDir = 0; iDir < nNbrDirs(); iDir++) {
    cout << "dir[" << iDir << "]" << endl;
    for(int iNbr = 0; iNbr < nNbrs(iDir); iNbr++) {
      cout << "nbr id " << _mNbrSpcIDs[iDir][iNbr] << endl;
    }
  }
  */ 
}

inline int pRPL::SubSpaceInfo::
id() const {
  return _id;
}

inline pRPL::DomDcmpType pRPL::SubSpaceInfo::
domDcmpType() const {
  return _domDcmpType;  
}

inline int pRPL::SubSpaceInfo::
nGlbRows() const {
  return _glbDims.nRows();
}

inline int pRPL::SubSpaceInfo::
nGlbCols() const {
  return _glbDims.nCols();
}

inline const pRPL::SpaceDims& pRPL::SubSpaceInfo::
glbDims() const {
  return _glbDims;
}

inline int pRPL::SubSpaceInfo::
nRows() const {
  return _dims.nRows();
}

inline int pRPL::SubSpaceInfo::
nCols() const {
  return _dims.nCols();
}

inline const pRPL::SpaceDims& pRPL::SubSpaceInfo::
dims() const {
  return _dims;
}

inline int pRPL::SubSpaceInfo::
iRowBegin() const {
  return _MBR.minIRow();
}

inline int pRPL::SubSpaceInfo::
iColBegin() const {
  return _MBR.minICol();
}

inline int pRPL::SubSpaceInfo::
iRowEnd() const {
  return _MBR.maxIRow();
}

inline int pRPL::SubSpaceInfo::
iColEnd() const {
  return _MBR.maxICol();
}

inline const pRPL::CoordBR& pRPL::SubSpaceInfo::
MBR() const {
  return _MBR;  
}

inline int pRPL::SubSpaceInfo::
iRowWorkBegin() const {
  return _workBR.minIRow();
}

inline int pRPL::SubSpaceInfo::
iColWorkBegin() const {
  return _workBR.minICol();
}

inline int pRPL::SubSpaceInfo::
iRowWorkEnd() const {
  return _workBR.maxIRow();
}

inline int pRPL::SubSpaceInfo::
iColWorkEnd() const {
  return _workBR.maxICol();
}

inline const pRPL::CoordBR& pRPL::SubSpaceInfo::
workBR() const {
  return _workBR;
}

#endif
