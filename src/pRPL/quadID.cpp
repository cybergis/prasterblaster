#include "quadID.h"

/***************************************************************************
* quadID.cpp
*
* Project: pRPL, v 1.0
* Purpose: Implementation for class pRPL::Quadrant, and pRPL::QuadID
* Author:  Qingfeng (Gene) Guan
* E-mail:  guanqf {at} gmail.com
****************************************************************************
*
* NOTE: The Neighbor-finding algorithm is based on Parthajit 
* Bhattacharya's (2001) thesis "Efficient Neighbor Finding Algorithms 
* in Quadtree and Octree"
*
****************************************************************************
* Copyright (c) 2008, Qingfeng Guan
* NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
* purposes, NO COMMERCIAL usages are allowed unless the author is 
* contacted and a permission is granted
* 
****************************************************************************/


/*********************************
*            Quadrant            *
**********************************/
const int pRPL::Quadrant::_tableSize;
const int pRPL::Quadrant::
_nbrPrstTable[_tableSize][_tableSize] 
  = {{1, 1, 0, 0},
     {1, 0, 1, 0},
     {0, 1, 0, 1},
     {0, 0, 1, 1}};

const pRPL::QuadDir pRPL::Quadrant::
_nbrIdxTable[_tableSize][_tableSize]
  = {{SW_QUAD, NE_QUAD, NE_QUAD, SW_QUAD},
     {SE_QUAD, NW_QUAD, NW_QUAD, SE_QUAD},
     {NW_QUAD, SE_QUAD, SE_QUAD, NW_QUAD},
     {NE_QUAD, SW_QUAD, SW_QUAD, NE_QUAD}};

/*
  upper -> SPECIAL1_DIR
  same -> SPECIAL2_DIR
*/
const pRPL::MeshDir pRPL::Quadrant::
_vtxNbrTable[_tableSize][_tableSize]
  = {{SPECIAL1_DIR, NORTH_DIR, WEST_DIR, SPECIAL2_DIR},
     {NORTH_DIR, SPECIAL1_DIR, SPECIAL2_DIR, EAST_DIR},
     {WEST_DIR, SPECIAL2_DIR, SPECIAL1_DIR, SOUTH_DIR},
     {SPECIAL2_DIR, EAST_DIR, SOUTH_DIR, SPECIAL1_DIR}};

int pRPL::Quadrant::
_prmDir2tblIdx(MeshDir dir) const {
  int iDir;
  switch(dir) {
    case NORTH_DIR:
      iDir = 0;
      break;
    case WEST_DIR:
      iDir = 1;
      break;
    case EAST_DIR:
      iDir = 2;
      break;
    case SOUTH_DIR:
      iDir = 3;
      break;
    default:
      cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid direction [" \
         << dir << "]" << endl;
      iDir = -1;
      break;
  }
  return iDir;
}

int pRPL::Quadrant::
_dglDir2tblIdx(MeshDir dir) const {
  int iDir;
  switch(dir) {
    case NORTHWEST_DIR:
      iDir = 0;
      break;
    case NORTHEAST_DIR:
      iDir = 1;
      break;
    case SOUTHWEST_DIR:
      iDir = 2;
      break;
    case SOUTHEAST_DIR:
      iDir = 3;
      break;
    default:
      cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid direction [" \
         << dir << "]" << endl;
      iDir = -1;
      break;
  }
  return iDir;
}

bool pRPL::Quadrant::
flip() {
  bool valid = true;
  switch(_quad) {
    case NW_QUAD:
      _quad = SE_QUAD;
      break;
    case NE_QUAD:
      _quad = SW_QUAD;
      break;
    case SW_QUAD:
      _quad = NE_QUAD;
      break;
    case SE_QUAD:
      _quad = NW_QUAD;
      break;
    case ERROR_QUAD:
      valid = false;
      break;
  }
  return valid;
}

int pRPL::Quadrant::
nbrPresent(MeshDir dir) const {
  int iDir = _prmDir2tblIdx(dir);
  if(iDir < 0 || iDir > 3) {
    return -1;
  }
  return _nbrPrstTable[_quad][iDir];
}

pRPL::QuadDir pRPL::Quadrant::
edgeNbrIdx(MeshDir dir) const {
  int iDir = _prmDir2tblIdx(dir);
  if(iDir < 0 || iDir > 3) {
    return ERROR_QUAD;
  }
  return _nbrIdxTable[_quad][iDir];
}

pRPL::MeshDir pRPL::Quadrant::
vtxNbrIdx(MeshDir dir) const {
  int iDir = _dglDir2tblIdx(dir);
  if(iDir < 0 || iDir > 3) {
    return ERROR_DIR;
  }
  return _vtxNbrTable[_quad][iDir];
}

/*********************************
*             QuadID             *
**********************************/
int pRPL::QuadID::
toInt() const {
  int iID = 0;
  for(int iQuad = 0; iQuad < size(); iQuad++) {
    iID = iID*10 + at(iQuad).dir() + 1;
  }
  return iID;
}

bool pRPL::QuadID::
geEdgeNbrID(QuadID &nbrID,
            int &count,
            MeshDir dir) const {
  nbrID.erase(nbrID.begin(), nbrID.end());
  nbrID.resize(size());
  count = 1;
  int i = size();
  while(i > 0) {
    int iNbrPrsnt = at(i-1).nbrPresent(dir);
    switch(iNbrPrsnt) {
      case 1:
        nbrID[i-1].set(at(i-1).edgeNbrIdx(dir));
        i--;
        count++;
        break;
      case 0:
        nbrID[i-1].set(at(i-1).edgeNbrIdx(dir));
        i--;
        while(i > 0) {
          nbrID[i-1] = at(i-1);
          i--;
        }
        return true;
      case -1:
        cerr << __FILE__ << " " << __FUNCTION__ \
             << " Error: invalid direction" << endl;
        return false;
    }
  }
  return false;
}

bool pRPL::QuadID::
geVtxNbrID(QuadID &nbrID, 
           int &count,
           MeshDir dir) const{
  nbrID.erase(nbrID.begin(), nbrID.end());
  nbrID.resize(size());
  count = 1;
  int i = size();
  while(i > 0) {
    nbrID[i-1] = at(i-1);
    if(!nbrID[i-1].flip()) {
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: invalid quadrant" << endl;
      return false;
    }
    MeshDir vtxDir = at(i-1).vtxNbrIdx(dir);
    switch(vtxDir) {
      case ERROR_DIR:
        cerr << __FILE__ << " " << __FUNCTION__ \
             << " Error: invalid direction" << endl;
        return false;
      case SPECIAL1_DIR: /* GO UPPER PARENTS */
        count++;
        i--;
        break;
      case SPECIAL2_DIR: /* SAME PARENT*/
        i--;
        while(i > 0) {
          nbrID[i-1] = at(i-1);
          i--;
        }
        return true;
      default:
        i--;
        QuadID tmpNbrID;
        QuadID tmpCurrentID(begin(), begin()+i);
        int tmpCount;
        if(!tmpCurrentID.geEdgeNbrID(tmpNbrID, tmpCount, vtxDir)) {
          return false;
        }
        while(i > 0) {
          nbrID[i-1] = tmpNbrID[i-1];
          i--;
        }
        count += tmpCount;
        return true;
    }
  }
  return false;
}

ostream& pRPL::
operator<<(ostream &os, const QuadID &quadID) {
  for(int iVal = 0; iVal < quadID.size(); iVal++) {
    os << quadID[iVal].dir();
    if(iVal != quadID.size()-1) {
      os << "-";
    }
  }
  return os;
}
