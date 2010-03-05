#ifndef QUADTREENODE_H
#define QUADTREENODE_H

/***************************************************************************
* quadtreeNode.h
*
* Project: pRPL, v 1.0
* Purpose: Header file for class pRPL::QuadTreeNode
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
#include "quadID.h"

namespace pRPL {
  class QuadTreeNode {
    public:
      QuadTreeNode()
        :_pRoot(0),
         _pParent(0),
         _workload(-1) {}
      QuadTreeNode(const QuadID &id,
                   const QuadTreeNode *pRoot,
                   const QuadTreeNode *pParent,
                   const SpaceDims &glbDims,
                   const CoordBR &MBR,
                   int workload = -1);
      QuadTreeNode(const QuadTreeNode& rhs);
      ~QuadTreeNode();

      QuadTreeNode& operator=(const QuadTreeNode& rhs);

      bool dividable() const;
      bool divide();

      const QuadID& id() const {
        return _id;
      }
      int intID() const {
        return _id.toInt();  
      }
      int level() const {
        return _id.size();
      }

      const QuadTreeNode *root() const {
        return _pRoot;
      }
      const QuadTreeNode *parent() const {
        return _pParent;
      }

      bool hasChildren() const {
        return !(_vpChildren.empty());
      }
      QuadTreeNode* child(QuadDir quad) const;

      const vector<const QuadTreeNode *> &nbrs(MeshDir dir) const {
        return _mpNbrs[dir];
      }
      
      const SpaceDims& glbDims() const {
        return _glbDims;  
      }
      const CoordBR& MBR() const {
        return _MBR;
      }
      const SpaceDims& dims() const {
        return _dims;
      }

      void workload(int workLoad) {
        _workload = workLoad;  
      }
      int workload() const {
        return _workload;
      }
      
      const vector<vector<const QuadTreeNode *> >& nbrsMap() const {
        return _mpNbrs;
      }

      bool leaves(vector<QuadTreeNode *> &vpLeaves) const;
      void updateNbrs();
      void nbrIDsMap(vector<IntVect> &mNbrIDs) const;

    protected:
      bool _leavesInDir(vector<const QuadTreeNode *> &vpLeaves,
                        MeshDir dir) const;
      const QuadTreeNode* _getNbr(const QuadID &nbrID, int count) const;

    protected:
      QuadID _id;
      const QuadTreeNode *_pRoot;
      const QuadTreeNode *_pParent;
      vector< QuadTreeNode* > _vpChildren;

      vector<vector<const QuadTreeNode *> > _mpNbrs;

      SpaceDims _glbDims;
      CoordBR _MBR;
      SpaceDims _dims;
      int _workload;
  };
};

inline pRPL::QuadTreeNode::
QuadTreeNode(const QuadID &id,
             const QuadTreeNode *pRoot,
             const QuadTreeNode *pParent,
             const SpaceDims &glbDims,
             const CoordBR &MBR,
             int workload)
  :_id(id),
   _pRoot(pRoot),
   _pParent(pParent),
   _workload(workload) {
  if(!glbDims.valid()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid SpaceDims [" \
         << glbDims << "]" \
         << endl;
    exit(-1);
  }
  _glbDims = glbDims;

  if(!MBR.valid(glbDims)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invliad MBR [" << MBR \
         << "] over CellSpace ["
         << glbDims << "]"\
         << endl;
    exit(-1);
  }
  _MBR = MBR;
  _dims = SpaceDims(_MBR.nRows(), _MBR.nCols());
}

inline pRPL::QuadTreeNode::
QuadTreeNode(const QuadTreeNode& rhs)
  :_id(rhs._id),
   _pRoot(rhs._pRoot),
   _pParent(rhs._pParent),
   _vpChildren(rhs._vpChildren),
   _mpNbrs(rhs._mpNbrs),
   _glbDims(rhs._glbDims),
   _MBR(rhs._MBR),
   _dims(rhs._dims),
   _workload(rhs._workload) {}

inline pRPL::QuadTreeNode::
~QuadTreeNode() {
  if(!_vpChildren.empty()) {
    for(int iChild = 0; iChild < _vpChildren.size(); iChild++) {
      delete _vpChildren[iChild];
    }
  }
}

#endif
