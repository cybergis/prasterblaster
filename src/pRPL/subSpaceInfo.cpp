#include "subSpaceInfo.h"

/***************************************************************************
* subSpaceInfo.cpp
*
* Project: pRPL, v 1.0
* Purpose: Implementation for class pRPL::SubSpaceInfo
* Author:  Qingfeng (Gene) Guan
* E-mail:  guanqf {at} gmail.com
****************************************************************************
* Copyright (c) 2008, Qingfeng Guan
* NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
* purposes, NO COMMERCIAL usages are allowed unless the author is 
* contacted and a permission is granted
* 
****************************************************************************/

pRPL::SubSpaceInfo& pRPL::SubSpaceInfo::
operator=(const SubSpaceInfo &rhs) {
  if(this != &rhs) {
    _id = rhs._id;
    _domDcmpType = rhs._domDcmpType;
    _glbDims = rhs._glbDims;
    _dims = rhs._dims;
    _MBR = rhs._MBR;
    _workBR = rhs._workBR;
    _mNbrSpcIDs = rhs._mNbrSpcIDs;
  }
  return *this;
}

bool pRPL::SubSpaceInfo::
operator==(const SubSpaceInfo &rhs) const {
  return (_id == rhs._id &&
          _domDcmpType == rhs._domDcmpType &&
          _glbDims == rhs._glbDims &&
          _dims == rhs._dims &&
          _MBR == rhs._MBR &&
          _workBR == rhs._workBR &&
          _mNbrSpcIDs == rhs._mNbrSpcIDs);
}

bool pRPL::SubSpaceInfo::
operator!=(const SubSpaceInfo &rhs) const {
  return !(operator==(rhs));
}

double pRPL::SubSpaceInfo::
sizeRatio() const {
  return (double)(_dims.size()) / (double)(_glbDims.size());
}

bool pRPL::SubSpaceInfo::
validGlbIdx(int glbIdx) const {
  bool valid = true;
  if(!_glbDims.validIdx(glbIdx)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: in SubSpace[" << _id \
         << "], glbIdx [" << glbIdx \
         << "] out of global boundary [" \
         << _glbDims << "]" << endl;
    valid = false;
  }
  return valid;
}

bool pRPL::SubSpaceInfo::
validGlbCoord(const CellCoord &glbCoord) const {
  bool valid = true;
  if(!glbCoord.valid(_glbDims)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: in SubSpace[" << _id \
         << "], glbCoord [" \
         << glbCoord << "] out of global boundary [" \
         << _glbDims << "]" << endl;
    valid = false;  
  }
  return valid; 
}

bool pRPL::SubSpaceInfo::
validGlbCoord(int iRowGlb, int iColGlb) const {
  return validGlbCoord(CellCoord(iRowGlb, iColGlb));
}

const pRPL::CellCoord pRPL::SubSpaceInfo::
glbIdx2glbCoord(int glbIdx) const {
  return CellCoord(glbIdx, _glbDims);
}

int pRPL::SubSpaceInfo::
glbCoord2glbIdx(const CellCoord &glbCoord) const {
  return glbCoord.toIdx(_glbDims);
}

int pRPL::SubSpaceInfo::
glbCoord2glbIdx(int iRowGlb, int iColGlb) const {
  return glbCoord2glbIdx(CellCoord(iRowGlb, iColGlb));
}

bool pRPL::SubSpaceInfo::
validIdx(int idx) const {
  bool valid = true;
  if(!_dims.validIdx(idx)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: in SubSpace[" << _id \
         << "], index [" << idx \
         << "] out of boundary [" \
         << _dims << "]" << endl;
    valid = false;
  }
  return valid;
}

bool pRPL::SubSpaceInfo::
validCoord(const CellCoord &coord) const {
  bool valid = true;
  if(!coord.valid(_dims)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: in SubSpace[" << _id \
         << "], coordinate [" \
         << coord << "] out of boundary [" \
         << _dims << "]" << endl;
    valid = false;
  }
  return valid;
}

bool pRPL::SubSpaceInfo::
validCoord(int iRow, int iCol) const {
  return validCoord(CellCoord(iRow, iCol));
}

int pRPL::SubSpaceInfo::
coord2idx(const CellCoord &coord) const {
  return coord.toIdx(_dims);
}

int pRPL::SubSpaceInfo::
coord2idx(int iRow, int iCol) const {
  return coord2idx(CellCoord(iRow, iCol));
}

const pRPL::CellCoord pRPL::SubSpaceInfo::
idx2coord(int idx) const {
  return CellCoord(idx, _dims);
}

const pRPL::CellCoord pRPL::SubSpaceInfo::
lclCoord2glbCoord(const CellCoord& lclCoord) const {
  return (lclCoord + _MBR.nwCorner()); 
}

const pRPL::CellCoord pRPL::SubSpaceInfo::
lclCoord2glbCoord(int iRowLcl, int iColLcl) const {
  return lclCoord2glbCoord(CellCoord(iRowLcl, iColLcl));  
}

const pRPL::CellCoord pRPL::SubSpaceInfo::
glbCoord2lclCoord(const CellCoord& glbCoord) const {
  return (glbCoord - _MBR.nwCorner());
}

const pRPL::CellCoord pRPL::SubSpaceInfo::
glbCoord2lclCoord(int iRowGlb, int iColGlb) const {
  return glbCoord2lclCoord(CellCoord(iRowGlb, iColGlb));
}

int pRPL::SubSpaceInfo::
lclCoord2glbIdx(const CellCoord &lclCoord) const {
  return glbCoord2glbIdx(lclCoord2glbCoord(lclCoord));
}

int pRPL::SubSpaceInfo::
lclCoord2glbIdx(int iRowLcl, int iColLcl) const {
  return lclCoord2glbIdx(CellCoord(iRowLcl, iColLcl));
}

const pRPL::CellCoord pRPL::SubSpaceInfo::
glbIdx2lclCoord(int glbIdx) const {
  return glbCoord2lclCoord(CellCoord(glbIdx, _glbDims));
}

int pRPL::SubSpaceInfo::
glbIdx2lclIdx(int glbIdx) const {
  return coord2idx(glbIdx2lclCoord(glbIdx));
}

int pRPL::SubSpaceInfo::
lclIdx2glbIdx(int lclIdx) const {
  return lclCoord2glbIdx(idx2coord(lclIdx));
}

int pRPL::SubSpaceInfo::
nNbrDirs() const {
  int nNbrDirects;
  switch(_domDcmpType) {
    case NON_DCMP:
      nNbrDirects = 0;
      break;
    case ROWWISE_DCMP:
    case COLWISE_DCMP:
      nNbrDirects = 2;
      break;
    case BLOCK_DCMP:
      nNbrDirects = 8;
      break;
    default:
      cerr << __FILE__ << __FUNCTION__ \
           << " Error: invalid decomposition type (" \
           << _domDcmpType << ")" << endl;
      nNbrDirects = -1;
      break;
  }
  return nNbrDirects;
}

int pRPL::SubSpaceInfo::
nEdges() const {
  int numEdges;
  switch(_domDcmpType) {
    case NON_DCMP:
      numEdges = 0;
      break;
    case ROWWISE_DCMP:
    case COLWISE_DCMP:
      numEdges = 2;
      break;
    case BLOCK_DCMP:
      numEdges = 4;
      break;
    default:
      cerr << __FILE__ << __FUNCTION__ \
           << " Error: invalid decomposition type (" \
           << _domDcmpType << ")" << endl;
      numEdges = -1;
      break;
  }
  return numEdges;
}

int pRPL::SubSpaceInfo::
nTotNbrs() const {
  int numNbrs = 0;
  vector<IntVect>::const_iterator iNbrMap = _mNbrSpcIDs.begin();
  while(iNbrMap != _mNbrSpcIDs.end()) {
    numNbrs += (*iNbrMap).size();
    iNbrMap++;
  }
  return numNbrs;
}

bool pRPL::SubSpaceInfo::
validSpcDir(int iDir) const {
  bool valid = true;
  if(iDir < 0 ||
     iDir >= nNbrDirs()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invliad direction (" << iDir \
         << ") of neighboring SubSpaces on SubSpace[" \
         << _id << "]" << endl;
    valid = false;
  }
  return valid;
}

bool pRPL::SubSpaceInfo::
hasNbrs(int iDir) const {
  bool hasIt = true;
  if(!validSpcDir(iDir) ||
     _mNbrSpcIDs[iDir].size() == 0) {
    hasIt = false;
  }
  return hasIt;
}

int pRPL::SubSpaceInfo::
nNbrs(int iDir) const {
  int numNbrs = 0;
  if(validSpcDir(iDir)) {
    numNbrs = _mNbrSpcIDs[iDir].size();
  }
  return numNbrs;
}

const pRPL::IntVect& pRPL::SubSpaceInfo::
nbrSubSpcIDs(int iDir) const {
  return _mNbrSpcIDs.at(iDir);
}

int pRPL::SubSpaceInfo::
nbrDir(int nbrID) const {
  int dir = -1;
  for(int iDir = 0; iDir < nNbrDirs(); iDir++) {
    const IntVect &vNbrs = nbrSubSpcIDs(iDir);
    if(std::find(vNbrs.begin(), vNbrs.end(), nbrID) != vNbrs.end()) {
      dir = iDir;
      break;
    }
  }
  return dir;
}

pRPL::MeshDir pRPL::SubSpaceInfo::
spcDir2MeshDir(int iDir) const {
  MeshDir opDir = ERROR_DIR;
  if(validSpcDir(iDir)) {
    switch(_domDcmpType) {
      case NON_DCMP:
        break;
      case ROWWISE_DCMP:
        opDir = (iDir == UPPER_DIR) ? NORTH_DIR : SOUTH_DIR;
        break;
      case COLWISE_DCMP:
        opDir = (iDir == LEFT_DIR) ? WEST_DIR : EAST_DIR;
        break;
      case BLOCK_DCMP:
        opDir = static_cast<MeshDir>(iDir);
        break;
      default:
        break;
    }
  }
  return opDir;
}

int pRPL::SubSpaceInfo::
oppositeDir(int iDir) const {
  int opDir = -1;
  if(validSpcDir(iDir)) {
    switch(_domDcmpType) {
      case NON_DCMP:
        break;
      case ROWWISE_DCMP:
      case COLWISE_DCMP:
        opDir = (iDir == 0) ? 1 : 0;
        break;
      case BLOCK_DCMP:
        opDir = static_cast<int>(pRPL::oppositeDir(static_cast<MeshDir>(iDir)));
        break;
      default:
        break;
    }
  }
  return opDir;
}

void pRPL::SubSpaceInfo::
add2IntVect(IntVect &vInfoPack) const {
  vInfoPack.push_back(_id);
  vInfoPack.push_back(_domDcmpType);
  vInfoPack.push_back(_glbDims.nRows());
  vInfoPack.push_back(_glbDims.nCols());
  vInfoPack.push_back(_MBR.nwCorner().toIdx(_glbDims));
  vInfoPack.push_back(_MBR.seCorner().toIdx(_glbDims));
  vInfoPack.push_back(_workBR.nwCorner().toIdx(_dims));
  vInfoPack.push_back(_workBR.seCorner().toIdx(_dims));

  for(int iDir = 0; iDir < nNbrDirs(); iDir++) {
    const IntVect &vNbrIDs = _mNbrSpcIDs[iDir];
    int nNbrs = vNbrIDs.size();
    vInfoPack.push_back(nNbrs);
    for(int iNbr = 0; iNbr < nNbrs; iNbr++) {
      vInfoPack.push_back(vNbrIDs[iNbr]);
    }
  }
}

bool pRPL::SubSpaceInfo::
fromIntVect(const IntVect &vInfoPack,
            IntVctItr &iVal) {
  if(vInfoPack.end() - iVal < 8 ||
     iVal < vInfoPack.begin()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid iterator of the integer vector" \
         << " to construct a SubSpaceInfo" << endl;
    return false;
  }

  _id = *iVal;
  iVal++;
  if(_id < 0) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid SubSpace ID (" \
         << _id << ")" << endl;
    return false;
  }
  
  _domDcmpType = static_cast<DomDcmpType>(*iVal);
  iVal++;

  _glbDims.nRows(*iVal);
  iVal++;
  _glbDims.nCols(*iVal);
  iVal++;
  if(!_glbDims.valid()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: SubSpace[" << _id \
         << "] invalid global SpaceDims (" \
         << _glbDims << ")" << endl;
    return false;
  }

  _MBR.nwCorner(CellCoord(*iVal, _glbDims));
  iVal++;
  _MBR.seCorner(CellCoord(*iVal, _glbDims));
  iVal++;
  if(!_MBR.valid(_glbDims)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: SubSpace[" << _id \
         << "] invalid MBR (" << _MBR \
         << ")" << endl;
    return false;
  }
  _dims = SpaceDims(_MBR.nRows(), _MBR.nCols());

  _workBR.nwCorner(CellCoord(*iVal, _dims));
  iVal++;
  _workBR.seCorner(CellCoord(*iVal, _dims));
  iVal++;
  if(!_workBR.valid(_dims)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: SubSpace[" << _id \
         << "] invalid workBR (" \
         << _workBR << ")" << endl;
    return false;
  }
  
  int numNbrDirs = nNbrDirs();
  _mNbrSpcIDs.resize(numNbrDirs);
  for(int iDir = 0; iDir < numNbrDirs; iDir++) {
    int nNbrs = *iVal;
    iVal++;
    
    for(int iNbr = 0; iNbr < nNbrs; iNbr++) {
      _mNbrSpcIDs[iDir].push_back(*iVal);
      iVal++;
    }
  }

  return true;
}
