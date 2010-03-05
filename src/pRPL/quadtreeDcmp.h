#ifndef QUADTREEDCMP_H
#define QUADTREEDCMP_H

/***************************************************************************
* quadtreeDcmp.h
*
* Project: pRPL, v 1.0
* Purpose: Header file for class pRPL::QuadTree
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

#include "basicTypes.h"
#include "cellSpace.h"
#include "subSpaceInfo.h"
#include "neighborhood.h"
#include "quadtreeNode.h"
#include "dcmpUtilities.h"
#include "subInfoVect.h"

namespace pRPL {
  template <class elemType>
  class QuadTree {
    public:
      QuadTree()
        :_pRoot(0),
         _maxNumLeaves(0),
         _minWorkload(0),
         _pCellSpace(0),
         _pNbrhood(0),
         _pTrans(0) {}

      QuadTree(CellSpace<elemType> &cellspace,
               Neighborhood<elemType> &nbrhood,
               Transition<elemType> &transition,
               int maxNumLeaves,
               int minWorkload = 0);
      ~QuadTree();

      const QuadTreeNode *root() const{
        return _pRoot;
      }
      const CellSpace<elemType> *cellSpace() const {
        return _pCellSpace;  
      }
      const Neighborhood<elemType> *nbrhood() const {
        return _pNbrhood;
      }

      int maxNumLeaves() const {
        return _maxNumLeaves;
      }

      int nLeaves() const {
        return _vpLeaves.size();
      }
      QuadTreeNode* leaf(int iLeaf) const;

      void findAllLeaves();
      bool updateWorkloads();

      void updateNbrs() const;
      bool grow();

      bool dcmp(SubInfoVect &vpSubSpcInfos) const;

    private:
      void _updateLeaves(int iLeaf);
      CoordBR _quadBR2subBR(const CoordBR &quadBR) const;
      bool _makeSubInfo(CoordBR &subMBR,
                        CoordBR &workBR,
                        int iLeaf) const;

    private:
      QuadTreeNode *_pRoot;
      vector< QuadTreeNode *> _vpLeaves;
      int _maxNumLeaves;
      int _minWorkload;

      CellSpace<elemType> *_pCellSpace;
      Neighborhood<elemType> *_pNbrhood;
      Transition<elemType> *_pTrans;
  };
};

/****************************************
*             Private Methods           *
*****************************************/

template <class elemType>
void pRPL::QuadTree<elemType>::
_updateLeaves(int iLeaf) {
  QuadTreeNode *pLeaf = _vpLeaves[iLeaf];
  if(pLeaf->hasChildren()) {
    for(int iQuad = 0; iQuad < 4; iQuad++) {
      _vpLeaves.push_back(pLeaf->child(static_cast<QuadDir>(iQuad)));
    }
    _vpLeaves.erase(_vpLeaves.begin() + iLeaf);
  }
}

template <class elemType>
pRPL::CoordBR pRPL::QuadTree<elemType>::
_quadBR2subBR(const CoordBR &quadBR) const {
  CellCoord nwCorner(quadBR.nwCorner().iRow() + _pNbrhood->minIRow(),
                     quadBR.nwCorner().iCol() + _pNbrhood->minICol());
  CellCoord seCorner(quadBR.seCorner().iRow() + _pNbrhood->maxIRow(),
                     quadBR.seCorner().iCol() + _pNbrhood->maxICol());
  return(CoordBR(nwCorner, seCorner));
}

template <class elemType>
bool pRPL::QuadTree<elemType>::
_makeSubInfo(CoordBR &subMBR,
             CoordBR &workBR,
             int iLeaf) const {
  QuadTreeNode *pLeaf = _vpLeaves[iLeaf];
  const vector<vector<const QuadTreeNode *> > &mpNbrs = pLeaf->nbrsMap();
  if(mpNbrs.empty()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: leaf[" << pLeaf->id() << "] has NOT built" \
         << " the neighbor map yet" << endl;
    return false;
  }

  const CoordBR &quadBR = pLeaf->MBR();
  subMBR = _quadBR2subBR(quadBR);
  SpaceDims dims(subMBR.nRows(), subMBR.nCols());
  if(!_pNbrhood->calcWorkBR(workBR, dims)) {
    return false;
  }

  return true;
}

/****************************************
*             Public Methods            *
*****************************************/

template <class elemType>
inline pRPL::QuadTree<elemType>::
QuadTree(CellSpace<elemType> &cellspace,
         Neighborhood<elemType> &nbrhood,
         Transition<elemType> &transition,
         int maxNumLeaves,
         int minWorkload)
  :_pTrans(&transition) {
  if(!checkDcmpParms(cellspace.dims(), nbrhood)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid parameters to construct a QuadTree" << endl;
    exit(-1);
  }
  _pCellSpace = &cellspace;
  _pNbrhood = &nbrhood;

  if(maxNumLeaves < 4) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: at least 4 leaves to construct a QuadTree" << endl;
    exit(-1);
  }
  _maxNumLeaves = maxNumLeaves;
  
  if(minWorkload < 0) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: the minimum workload should be non-negative" << endl;
    exit(-1);
  }
  _minWorkload = minWorkload;
  
  QuadID rootID;
  SpaceDims glbDims(cellspace.dims());
  CoordBR rootMBR;
  if(!nbrhood.calcWorkBR(rootMBR, glbDims)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid Neighborhood to construct a QuadTree" << endl;
    exit(-1);
  }
  
  _pRoot = new QuadTreeNode(rootID, 0, 0,
                            glbDims, rootMBR);
  if(!(_pRoot->divide())) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to divide the root" \
         << endl;
    delete _pRoot;
    exit(-1);
  }
  findAllLeaves();
  updateWorkloads();
}

template <class elemType>
inline pRPL::QuadTree<elemType>::
~QuadTree() {
  if(!_pRoot) {
    delete _pRoot;
  }
}

template <class elemType>
pRPL::QuadTreeNode* pRPL::QuadTree<elemType>::
leaf(int iLeaf) const {
  if(iLeaf < 0 || iLeaf >= _vpLeaves.size()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid leaf index[" \
         << iLeaf << "]. there are only " \
         << _vpLeaves.size() << " leaves" << endl;
    return 0;
  }
  return _vpLeaves[iLeaf];
}

template <class elemType>
void pRPL::QuadTree<elemType>::
findAllLeaves() {
  _pRoot->leaves(_vpLeaves);
}

template <class elemType>
bool pRPL::QuadTree<elemType>::
updateWorkloads() {
  if(!_pTrans->cellSpace(_pCellSpace) ||
     !_pTrans->nbrhood(_pNbrhood)) {
    return false;
  }

  for(int iLeaf = 0; iLeaf < _vpLeaves.size(); iLeaf++) {
    if(_vpLeaves[iLeaf]->workload() < 0) {
      const CoordBR &mbr = _vpLeaves[iLeaf]->MBR();
      int wl = _pTrans->workload(mbr);
      _vpLeaves[iLeaf]->workload(wl); 
    }
  }
  return true;
}

template <class elemType>
void pRPL::QuadTree<elemType>::
updateNbrs() const {
  for(int iLeaf = 0; iLeaf < _vpLeaves.size(); iLeaf++) {
    _vpLeaves[iLeaf]->updateNbrs();
  }
}

template <class elemType>
bool pRPL::QuadTree<elemType>::
grow() {
  if(_maxNumLeaves-_vpLeaves.size() < 3) {
    cout << "QuadTree growth stops. " \
         << "The maximum number of leaves has been reached. " \
         << nLeaves() << " leaves have been built" << endl;
    return false;
  }
  vector<QuadTreeNode *>::iterator iLeaf = _vpLeaves.begin();
  int maxLeaf = 0;
  int maxWL = (*iLeaf)->workload();
  for(; iLeaf != _vpLeaves.end(); iLeaf++) {
    if(maxWL < (*iLeaf)->workload()) {
      maxLeaf = iLeaf - _vpLeaves.begin();
      maxWL = (*iLeaf)->workload();
    }
  }
  if(maxWL <= _minWorkload) {
    cout << "QuadTree growth stops. " \
         << "The minimum workload has been reached. " \
         << nLeaves() << " leaves have been built" << endl;
    return false;
  }

  QuadTreeNode *pMaxLeaf = _vpLeaves[maxLeaf];
  if(!(pMaxLeaf->divide())) {
    return false;
  }
  _updateLeaves(maxLeaf);
  
  if(!updateWorkloads()) {
    return false;
  }

  return true;
}

template <class elemType>
bool pRPL::QuadTree<elemType>::
dcmp(SubInfoVect &vpSubSpcInfos) const {
  vpSubSpcInfos.clear();
  DomDcmpType dcmpType = BLOCK_DCMP;
  const SpaceDims &glbDims = _pRoot->glbDims();
  for(int iLeaf = 0; iLeaf < _vpLeaves.size(); iLeaf++) {
    CoordBR subMBR, workBR;
    if(!_makeSubInfo(subMBR, workBR, iLeaf)) {
      return false;
    }

    int subID = _vpLeaves[iLeaf]->intID();
    vector<IntVect> mNbrIDs;
    _vpLeaves[iLeaf]->nbrIDsMap(mNbrIDs);

    int workload = _vpLeaves[iLeaf]->workload();
    SubInfoVect::iterator iInsert = upper_bound(vpSubSpcInfos.begin(),
                                                vpSubSpcInfos.end(),
                                                InfoWLPair(0, workload),
                                                greater<InfoWLPair>());
    vpSubSpcInfos.insert(iInsert,
      InfoWLPair(new SubSpaceInfo(subID, dcmpType, glbDims, 
                                  subMBR, workBR, mNbrIDs),
                 workload));
  }
  return true;
}

#endif
