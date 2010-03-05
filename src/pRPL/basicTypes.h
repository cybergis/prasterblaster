#ifndef BASICTYPES_H
#define BASICTYPES_H

/***************************************************************************
* basicTypes.h
*
* Project: pRPL, v 1.0
* Purpose: Header file for class pRPL::DomDcmpType, pRPL::DomDcmpMethod,
*          pRPL::SglColDir, pRPL::SglRowDir, pRPL::MeshDir, pRPL::PrimeDir,
*          pRPL::SpaceDims, pRPL::CellCoord, and pRPL::CoordBR
* Author:  Qingfeng (Gene) Guan
* E-mail:  guanqf {at} gmail.com
****************************************************************************
* Copyright (c) 2008, Qingfeng Guan
* NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
* purposes, NO COMMERCIAL usages are allowed unless the author is 
* contacted and a permission is granted
* 
****************************************************************************/

#include <utility>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string.h>
#include <cmath>

using namespace std;

namespace pRPL {
  enum DomDcmpType{NON_DCMP = 0,
                   ROWWISE_DCMP,
                   COLWISE_DCMP,
                   BLOCK_DCMP};

  enum DomDcmpMethod{SINGLE_PRC = 0,
                     SMPL_ROW,
                     SMPL_COL,
                     SMPL_BLK,
                     CMPX_QUAD};

  ostream& operator<<(ostream &os, const DomDcmpMethod &dcmpMethod);

  enum SglColDir{UPPER_DIR = 0,
                 LOWER_DIR};

  enum SglRowDir{LEFT_DIR = 0,
                 RIGHT_DIR};

  enum MeshDir{ERROR_DIR = -1,
               NORTH_DIR = 0,
               NORTHEAST_DIR,
               EAST_DIR,
               SOUTHEAST_DIR,
               SOUTH_DIR,
               SOUTHWEST_DIR,
               WEST_DIR,
               NORTHWEST_DIR,
               SPECIAL1_DIR,
               SPECIAL2_DIR};

  enum PrimeDir{NORTH_PDIR = 0,
                EAST_PDIR,
                SOUTH_PDIR,
                WEST_PDIR};

  MeshDir oppositeDir(MeshDir dir);

  typedef vector<int> IntVect;
  typedef vector<int>::iterator IntVctItr;
  ostream& operator<<(ostream &os, const IntVect &vInt);

  typedef vector<double> DblVect;
  typedef vector<double>::iterator DblVectItr;
  ostream& operator<<(ostream &os, const DblVect &vDbl);

  typedef vector<char> CharVect;
  typedef vector<char>::iterator CharVctItr;
  ostream& operator<<(ostream &os, const CharVect &vChar);

  typedef pair<int, int> IntPair;
  ostream& operator<<(ostream &os, const IntPair &pair);

  typedef map<int, int> IntMap;
  typedef map<int, int>::value_type IntMapValType;
  typedef map<int, int>::iterator IntMapItr;
  ostream& operator<<(ostream &os, const IntMap &mInt);

  static const int ERROR_SPCID = -1;
  static const int ERROR_COORD = -9999;

  /* class SpaceDims */
  class SpaceDims {
    public:
      SpaceDims()
        :_nRows(0),
         _nCols(0),
         _size(0) {}
      SpaceDims(int nRows, int nCols)
        :_nRows(nRows),
         _nCols(nCols),
         _size(nRows*nCols) {}
      SpaceDims(const SpaceDims& rhs)
        :_nRows(rhs._nRows),
         _nCols(rhs._nCols),
         _size(rhs._size) {}
      ~SpaceDims() {}

      SpaceDims& operator=(const SpaceDims& rhs) {
        if(this != &rhs) {
          _nRows = rhs._nRows;
          _nCols = rhs._nCols;
          _size = rhs._size;
        }
        return *this;
      }
      bool operator==(const SpaceDims& rhs) const {
        return (_nRows == rhs._nRows &&
                _nCols == rhs._nCols);
      }
      bool operator!=(const SpaceDims &rhs) const {
        return (_nRows != rhs._nRows ||
                _nCols != rhs._nCols);
      }
      void nRows(int numRows) {
        _nRows = numRows;
        _size = _nRows * _nCols;
      }
      void nCols(int numCols) {
        _nCols = numCols;
        _size = _nRows * _nCols;
      }
      int nRows() const {
        return _nRows;
      }
      int nCols() const {
        return _nCols;
      }
      int size() const {
        return _size;
      }
      bool valid() const {
        return (_nRows >= 0 && _nCols >= 0);
      }
      bool isNone() const {
        return (_size == 0);
      }
      bool validIdx(int idx) const {
        return (idx >= 0 && idx < _size);
      }

    private:
      int _nRows;
      int _nCols;
      int _size;
  };

  /* class CellCoord */
  class CellCoord {
    public:
      CellCoord()
        :_iRow(ERROR_COORD),
         _iCol(ERROR_COORD) {}
      CellCoord(int iRow, int iCol)
        :_iRow(iRow),
         _iCol(iCol) {}
      CellCoord(const CellCoord& rhs)
        :_iRow(rhs._iRow),
         _iCol(rhs._iCol) {}
      CellCoord(int idx,
                const SpaceDims& dims) {
        fromIdx(idx, dims);
      }
      ~CellCoord() {}

      CellCoord& operator=(const CellCoord &rhs) {
        if(this != &rhs) {
          _iRow = rhs._iRow;
          _iCol = rhs._iCol;
        }
        return *this;
      }
      void operator()(const CellCoord &rhs) {
        _iRow = rhs._iRow;
        _iCol = rhs._iCol;
      }
      void operator()(int iRow, int iCol) {
        _iRow = iRow;
        _iCol = iCol;
      }
      bool operator==(const CellCoord &rhs) const {
        return (_iRow == rhs._iRow &&
                _iCol == rhs._iCol);
      }
      bool operator!=(const CellCoord &rhs) const {
        return (_iRow != rhs._iRow ||
                _iCol != rhs._iCol);
      }
      CellCoord operator+(const CellCoord &rhs) const{
        return CellCoord(_iRow + rhs._iRow,
                         _iCol + rhs._iCol);
      }
      CellCoord operator-(const CellCoord &rhs) const{
        return CellCoord(_iRow - rhs._iRow,
                         _iCol - rhs._iCol);
      }
      void iRow(int idxRow) {
        _iRow = idxRow;
      }
      void iCol(int idxCol) {
        _iCol = idxCol;
      }
      int iRow() const {
        return _iRow;
      }
      int iCol() const {
        return _iCol;
      }

      bool valid(const SpaceDims& dims) const {
        return (_iRow >= 0 &&
                _iRow < dims.nRows() &&
                _iCol >= 0 &&
                _iCol < dims.nCols());
      }

      int toIdx(const SpaceDims& dims) const {
        int idx;
        if(!valid(dims)) {
          idx = ERROR_COORD;
        }
        else {
          idx = _iRow*dims.nCols() + _iCol;
        }
        return idx;
      }

      void fromIdx(int idx,
                   const SpaceDims& dims) {
        if(!dims.validIdx(idx)) {
          _iRow = ERROR_COORD;
          _iCol = ERROR_COORD;
        }
        else {
          _iRow = idx / dims.nCols();
          _iCol = idx - _iRow*dims.nCols();
          if(!valid(dims)) {
            _iRow = ERROR_COORD;
            _iCol = ERROR_COORD;
          }
        }
      }

    private:
      int _iRow;
      int _iCol;
  };

  /* class CoordBR */
  class CoordBR {
    public:
      CoordBR() {}
      CoordBR(const CellCoord& nw,
              const CellCoord& se)
        :_nwCorner(nw),
         _seCorner(se) {}
      CoordBR(const CoordBR& rhs)
        :_nwCorner(rhs._nwCorner),
         _seCorner(rhs._seCorner) {}
      ~CoordBR() {}

      CoordBR& operator=(const CoordBR &rhs) {
        if(this != &rhs) {
          _nwCorner = rhs._nwCorner;
          _seCorner = rhs._seCorner;
        }
        return *this;
      }
      bool operator==(const CoordBR &rhs) const {
        return (_nwCorner == rhs._nwCorner &&
                _seCorner == rhs._seCorner);
      }
      bool operator!=(const CoordBR &rhs) const {
        return (_nwCorner != rhs._nwCorner ||
                _seCorner != rhs._seCorner);
      }
      void nwCorner(const CellCoord &nw) {
        _nwCorner = nw;
      }
      void nwCorner(int iRow, int iCol) {
        _nwCorner(iRow, iCol);
      }
      void seCorner(const CellCoord &se) {
        _seCorner = se;
      }
      void seCorner(int iRow, int iCol) {
        _seCorner(iRow, iCol);
      }
      
      const CellCoord& nwCorner() const {
        return _nwCorner;  
      }
      const CellCoord& seCorner() const {
        return _seCorner;  
      }
      
      int minIRow() const {
        return _nwCorner.iRow();
      }
      int minICol() const {
        return _nwCorner.iCol();
      }
      int maxIRow() const {
        return _seCorner.iRow();
      }
      int maxICol() const {
        return _seCorner.iCol();
      }

      int nRows() const {
        return _seCorner.iRow() - _nwCorner.iRow() + 1;
      }
      int nCols() const {
        return _seCorner.iCol() - _nwCorner.iCol() + 1;
      }
      int size() const {
        return nRows() * nCols();
      }

      bool valid() const {
        return (_nwCorner.iRow() <= _seCorner.iRow() &&
                _nwCorner.iCol() <= _seCorner.iCol());
      }
      bool valid(const SpaceDims& dims) const {
        return (_nwCorner.valid(dims) &&
                _seCorner.valid(dims) &&
                valid());
      }
      bool contain(const CellCoord &coord) const {
        return (coord.iRow() >= _nwCorner.iRow() &&
                coord.iRow() <= _seCorner.iRow() &&
                coord.iCol() >= _nwCorner.iCol() &&
                coord.iCol() <= _seCorner.iCol());
      }
      bool contain(const CoordBR &rhs) const {
        return (rhs.valid() &&
                contain(rhs._nwCorner) &&
                contain(rhs._seCorner));
      }
      bool within(const CoordBR &rhs) const {
        return (rhs.contain(*this));
      }

      CoordBR overlap(const CoordBR &rhs) const {
        CoordBR olBR;
        if(valid() && rhs.valid()) {
          int iWRow = max(minIRow(), rhs.minIRow());
          int iERow = min(maxIRow(), rhs.maxIRow());
          int iNCol = max(minICol(), rhs.minICol());
          int iSCol = min(maxICol(), rhs.maxICol());
          olBR.nwCorner(iWRow, iNCol);
          olBR.seCorner(iERow, iSCol);
          if(!olBR.valid()) {
            olBR = CoordBR();
          }
        }
        return olBR;
      }

      int overlapArea(const CoordBR &rhs) const {
        CoordBR olBR = overlap(rhs);
        if(olBR == CoordBR()) {
          return 0;
        }
        return olBR.size();
      }

    private:
      CellCoord _nwCorner;
      CellCoord _seCorner;
  };

  ostream& operator<<(ostream &os, const SpaceDims &dims);

  ostream& operator<<(ostream &os, const CellCoord &coord);

  ostream& operator<<(ostream &os, const CoordBR &br);
};

inline pRPL::MeshDir pRPL::
oppositeDir(MeshDir dir) {
  MeshDir opDir = ERROR_DIR;
  if(dir >= NORTH_DIR &&
     dir <= NORTHWEST_DIR) {
    int tmpDir = dir + 4;
    if(tmpDir >= 8) {
      tmpDir -= 8;
    }
    opDir = static_cast<MeshDir>(tmpDir);
  }
  return opDir;
}

inline ostream& pRPL::
operator<<(ostream &os, const IntVect &vInt) {
  for(int iVal = 0; iVal < vInt.size(); iVal++) {
    os << vInt[iVal] << " ";
  }
  return os;
}

inline ostream& pRPL::
operator<<(ostream &os, const DblVect &vDbl) {
  for(int iVal = 0; iVal < vDbl.size(); iVal++) {
    os << vDbl[iVal] << " ";
  }
  return os;
}

inline ostream& pRPL::
operator<<(ostream &os, const CharVect &vChar) {
  for(int iVal = 0; iVal < vChar.size(); iVal++) {
    os << vChar[iVal] << " ";
  }
  return os;
}

inline ostream& pRPL::
operator<<(ostream &os, const IntPair &pair) {
  os << pair.first << " " << pair.second;
  return os;
}

inline ostream& pRPL::
operator<<(ostream &os, const IntMap &mInt) {
  IntMap::const_iterator iVal = mInt.begin();
  while(iVal != mInt.end()) {
    os << "{" << (*iVal).first << " " << (*iVal).second << "} ";
    iVal++;
  }
  return os;
}

inline ostream& pRPL::
operator<<(ostream &os, const SpaceDims &dims) {
  os << dims.nRows() << " " << dims.nCols();
  return os;
}

inline ostream& pRPL::
operator<<(ostream &os, const CellCoord &coord) {
  os << coord.iRow() << " " << coord.iCol();
  return os;
}

inline ostream& pRPL::
operator<<(ostream &os, const CoordBR &br) {
  os << br.nwCorner() << "; " << br.seCorner();
  return os;
}

inline ostream& pRPL::
operator<<(ostream &os, const DomDcmpMethod &dcmpMethod) {
  switch(dcmpMethod) {
    case SINGLE_PRC:
      os << "NON_DCMP";
      break;
    case SMPL_ROW:
      os << "ROW_WISE";
      break;
    case SMPL_COL:
      os << "COL_WISE";
      break;
    case SMPL_BLK:
      os << "BLK_WISE";
      break;
    case CMPX_QUAD:
      os << "QUAD_TREE";
      break;
  }
  return os;
}

#endif
