#include "quadtreeNode.h"

/***************************************************************************
* quadtreeNode.cpp
*
* Project: pRPL, v 1.0
* Purpose: Implementation for class pRPL::QuadTreeNode
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


/****************************************
*             Private Methods           *
*****************************************/
bool pRPL::QuadTreeNode::
_leavesInDir(vector<const QuadTreeNode *> &vpLeaves,
             MeshDir dir) const {
  if(_vpChildren.empty()) {
    return false;
  }
  QuadTreeNode *nbr;
  switch(dir) {
    case NORTH_DIR:
      nbr = _vpChildren[NW_QUAD];
      if(!nbr->_leavesInDir(vpLeaves, dir)) {
        vpLeaves.push_back(nbr);
      }
      nbr = _vpChildren[NE_QUAD];
      if(!nbr->_leavesInDir(vpLeaves, dir)) {
        vpLeaves.push_back(nbr);
      }
      break;
    case WEST_DIR:
      nbr = _vpChildren[NW_QUAD];
      if(!nbr->_leavesInDir(vpLeaves, dir)) {
        vpLeaves.push_back(nbr);
      }
      nbr = _vpChildren[SW_QUAD];
      if(!nbr->_leavesInDir(vpLeaves, dir)) {
        vpLeaves.push_back(nbr);
      }
      break;
    case EAST_DIR:
      nbr = _vpChildren[NE_QUAD];
      if(!nbr->_leavesInDir(vpLeaves, dir)) {
        vpLeaves.push_back(nbr);
      }
      nbr = _vpChildren[SE_QUAD];
      if(!nbr->_leavesInDir(vpLeaves, dir)) {
        vpLeaves.push_back(nbr);
      }
      break;
    case SOUTH_DIR:
      nbr = _vpChildren[SW_QUAD];
      if(!nbr->_leavesInDir(vpLeaves, dir)) {
        vpLeaves.push_back(nbr);
      }
      nbr = _vpChildren[SE_QUAD];
      if(!nbr->_leavesInDir(vpLeaves, dir)) {
        vpLeaves.push_back(nbr);
      }
      break;
    case NORTHWEST_DIR:
      nbr = _vpChildren[NW_QUAD];
      if(!nbr->_leavesInDir(vpLeaves, dir)) {
        vpLeaves.push_back(nbr);
      }
      break;
    case NORTHEAST_DIR:
      nbr = _vpChildren[NE_QUAD];
      if(!nbr->_leavesInDir(vpLeaves, dir)) {
        vpLeaves.push_back(nbr);
      }
      break;
    case SOUTHWEST_DIR:
      nbr = _vpChildren[SW_QUAD];
      if(!nbr->_leavesInDir(vpLeaves, dir)) {
        vpLeaves.push_back(nbr);
      }
      break;
    case SOUTHEAST_DIR:
      nbr = _vpChildren[SE_QUAD];
      if(!nbr->_leavesInDir(vpLeaves, dir)) {
        vpLeaves.push_back(nbr);
      }
      break;
    default:
      cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid direction" << endl;
      return false;
  }
  return true;
}

const pRPL::QuadTreeNode* pRPL::QuadTreeNode::
_getNbr(const QuadID &nbrID, int count) const {
  const QuadTreeNode *nbr = this;
  int depth = level();

  if(count > depth/2) {
    nbr = _pRoot;
    for(int i = 0; i < depth; i++) {
      if(nbr->_vpChildren.empty()) {
        return nbr;
      }
      nbr = nbr->_vpChildren[nbrID[i].dir()];
    }
  }
  else {
    for(int i = 0; i < count; i++) {
      nbr = nbr->_pParent;
    }
    for(int i = depth-count; i < depth; i++) {
      if(nbr->_vpChildren.empty()) {
        return nbr;
      }
      nbr = nbr->_vpChildren[nbrID[i].dir()];
    }
  }

  return nbr;
}

/****************************************
*             Public Methods            *
*****************************************/
pRPL::QuadTreeNode& pRPL::QuadTreeNode::
operator=(const QuadTreeNode& rhs) {
  if(this != &rhs) {
    _id = rhs._id;
    _pRoot = rhs._pRoot;
    _pParent = rhs._pParent;
    _vpChildren = rhs._vpChildren;
    _mpNbrs = rhs._mpNbrs;
    _glbDims = rhs._glbDims;
    _MBR = rhs._MBR;
    _dims = rhs._dims;
    _workload = rhs._workload;
  }
  return *this;
}

pRPL::QuadTreeNode* pRPL::QuadTreeNode::
child(QuadDir quad) const {
  if(!hasChildren()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: QuadTreeNode [" << _id \
         << "] has no children"
         << endl;
    return 0;
  }
  if(quad < NW_QUAD || quad > SE_QUAD) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid direction (" \
         << quad << ")" \
         << endl;
    return 0;
  }
  return _vpChildren[quad];
}

bool pRPL::QuadTreeNode::
dividable() const {
  return (_dims.nRows() >= 2 &&
          _dims.nCols() >= 2);
}

bool pRPL::QuadTreeNode::
divide() {
  if(!_vpChildren.empty()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: QuadTreeNode [" << _id \
         << "] has been divided already"
         << endl;
    return false;
  }
  if(!dividable()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: QuadTreeNode [" << _id \
         << "] can NOT be divided." \
         << " Dims [" << _dims << "]" \
         << endl;
    return false;
  }

  int nHalfRows = _dims.nRows() / 2;
  int nHalfCols = _dims.nCols() / 2;
  int cntrIRow = _MBR.minIRow() + nHalfRows - 1;
  int cntrICol = _MBR.minICol() + nHalfCols - 1;

  QuadID nwID(_id); nwID.push_back(Quadrant(NW_QUAD));
  QuadID neID(_id); neID.push_back(Quadrant(NE_QUAD));
  QuadID swID(_id); swID.push_back(Quadrant(SW_QUAD));
  QuadID seID(_id); seID.push_back(Quadrant(SE_QUAD));

  CoordBR nwMBR, neMBR, swMBR, seMBR;
  nwMBR.nwCorner(_MBR.nwCorner());
  nwMBR.seCorner(CellCoord(cntrIRow, cntrICol));
  neMBR.nwCorner(CellCoord(_MBR.minIRow(), cntrICol+1));
  neMBR.seCorner(CellCoord(cntrIRow, _MBR.maxICol()));
  swMBR.nwCorner(CellCoord(cntrIRow+1, _MBR.minICol()));
  swMBR.seCorner(CellCoord(_MBR.maxIRow(), cntrICol));
  seMBR.nwCorner(CellCoord(cntrIRow+1, cntrICol+1));
  seMBR.seCorner(_MBR.seCorner());

  const QuadTreeNode *pRoot;
  if(_pRoot) {
    pRoot = _pRoot;
  }
  else {
    pRoot = this;
  }

  _vpChildren.push_back(new QuadTreeNode(nwID, pRoot, this,
                                         _glbDims, nwMBR));
  _vpChildren.push_back(new QuadTreeNode(neID, pRoot, this,
                                         _glbDims, neMBR));
  _vpChildren.push_back(new QuadTreeNode(swID, pRoot, this,
                                         _glbDims, swMBR));
  _vpChildren.push_back(new QuadTreeNode(seID, pRoot, this,
                                         _glbDims, seMBR));

  return true;
}

void pRPL::QuadTreeNode::
updateNbrs() {
  _mpNbrs.clear();
  _mpNbrs.resize(8);
  for(int iDir = 0; iDir < 8; iDir++) {
    MeshDir dir = static_cast<MeshDir>(iDir);
    QuadID geNbrID;
    int count;
    bool gotID;
    if(iDir % 2) {
      gotID = _id.geVtxNbrID(geNbrID, count, dir);  
    }
    else {
      gotID = _id.geEdgeNbrID(geNbrID, count, dir);
    }
    if(gotID) {
      const QuadTreeNode *pNbr = _getNbr(geNbrID, count);
      if(!(pNbr->_vpChildren.empty())) {
        MeshDir opDir = oppositeDir(dir);
        pNbr->_leavesInDir(_mpNbrs[dir], opDir);
      }
      else {
        _mpNbrs[dir].push_back(pNbr);
      }
    }
  }

  for(int iDir = 1; iDir < 9; iDir+=2) {
    vector<const QuadTreeNode *> &vpNbrs = _mpNbrs[iDir];
    int iPrmDir1 = iDir - 1;
    int iPrmDir2 = (iDir == 7) ? 0 : (iDir + 1);
    vector<const QuadTreeNode *> &prmNbrs1 = _mpNbrs[iPrmDir1];
    vector<const QuadTreeNode *> &prmNbrs2 = _mpNbrs[iPrmDir2];

    vector<const QuadTreeNode *>::iterator iNbr = vpNbrs.begin();
    while(iNbr != vpNbrs.end()) {
      if(std::find(prmNbrs1.begin(), prmNbrs1.end(), *iNbr) != prmNbrs1.end() ||
         std::find(prmNbrs2.begin(), prmNbrs2.end(), *iNbr) != prmNbrs2.end()) {
        vpNbrs.erase(iNbr);
      }
      else {
        iNbr++;
      }
    }
  }
}

bool pRPL::QuadTreeNode::
leaves(vector<QuadTreeNode *> &vpLeaves) const {
  if(_vpChildren.empty()) {
    return false;
  }
  for(int iQuad = 0; iQuad < 4; iQuad++) {
    QuadTreeNode *nbr = _vpChildren[iQuad];
    if(!nbr->leaves(vpLeaves)) {
      vpLeaves.push_back(nbr);
    }
  }
  return true;
}

void pRPL::QuadTreeNode::
nbrIDsMap(vector<IntVect> &mNbrIDs) const {
  mNbrIDs.clear();
  mNbrIDs.resize(8);
  if(!_mpNbrs.empty()) {
    for(int iDir = 0; iDir < 8; iDir++) {
      const vector<const QuadTreeNode *> &vpNbrs = _mpNbrs[iDir];
      for(int iNbr = 0; iNbr < vpNbrs.size(); iNbr++) {
        mNbrIDs[iDir].push_back(vpNbrs[iNbr]->intID());
      }
    }
  }
}

