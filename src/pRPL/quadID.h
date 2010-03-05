#ifndef QUADID_H
#define QUADID_H

/***************************************************************************
* quadID.h
*
* Project: pRPL, v 1.0
* Purpose: Header file for class pRPL::QuadDir, pRPL::Quadrant, 
*          and pRPL::QuadID
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

namespace pRPL {
  enum QuadDir{ERROR_QUAD = -1,
               NW_QUAD = 0,
               NE_QUAD,
               SW_QUAD,
               SE_QUAD};

  /* Quadrant */
  class Quadrant {
    public:
      Quadrant()
        :_quad(ERROR_QUAD) {}
      Quadrant(QuadDir quad)
        :_quad(quad) {}
      Quadrant(const Quadrant &rhs)
        :_quad(rhs._quad) {}
      ~Quadrant() {}

      Quadrant& operator=(const Quadrant &rhs) {
        if(this != &rhs) {
          _quad = rhs._quad;  
        }
        return *this;
      }
      bool operator==(const Quadrant &rhs) const {
        return (_quad == rhs._quad);  
      }

      QuadDir dir() const {
        return _quad;
      }
      void set(QuadDir quad) {
        _quad = quad;
      }

      bool flip();
      int nbrPresent(MeshDir dir) const;
      QuadDir edgeNbrIdx(MeshDir dir) const;
      MeshDir vtxNbrIdx(MeshDir dir) const;

    private:
      int _prmDir2tblIdx(MeshDir dir) const;
      int _dglDir2tblIdx(MeshDir dir) const;

    private:
      QuadDir _quad;
      static const int _tableSize = 4;
      static const int _nbrPrstTable[_tableSize][_tableSize];
      static const QuadDir _nbrIdxTable[_tableSize][_tableSize];
      static const MeshDir _vtxNbrTable[_tableSize][_tableSize];
  };

  /* QuadID */
  class QuadID: public vector<Quadrant> {
    public:
      QuadID()
        :vector<Quadrant>() {}
      QuadID(vector<Quadrant>::const_iterator begin,
             vector<Quadrant>::const_iterator end)
        :vector<Quadrant>(begin, end) {}
      ~QuadID() {}

      int toInt() const;

      bool geEdgeNbrID(QuadID &nbrID, 
                       int &count,
                       MeshDir dir) const;
      bool geVtxNbrID(QuadID &nbrID, 
                      int &count,
                      MeshDir dir) const;
  };

  ostream& operator<<(ostream &os, const QuadID &quadID);
};

#endif
