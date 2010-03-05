#ifndef NEIGHBORHOOD_H
#define NEIGHBORHOOD_H

/***************************************************************************
* neighborhood.h
*
* Project: pRPL, v 1.0
* Purpose: Header file for class pRPL::Neighborhood
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
#include "weightedCell.h"
#include "cellSpace.h"

namespace pRPL {
  template <class elemType>
  class Neighborhood {
    public:
      // Constructor
      Neighborhood();
      Neighborhood(const vector<CellCoord> &vNbrCoords,
                   double weight = 1.0);
      Neighborhood(const vector<CellCoord> &vNbrCoords,
                   const vector<double> &vNbrWeights);

      Neighborhood(const Neighborhood<elemType> &rhs);
      
      template<class elemType2>
      Neighborhood(const Neighborhood<elemType2> &rhs);

      // Deconstructor
      ~Neighborhood() {}

      WeightedCell<elemType>& operator[](int iNbr);
      const WeightedCell<elemType>& operator[](int iNbr) const;
      Neighborhood<elemType>& operator=(const Neighborhood<elemType> &rhs);

      bool empty() const;
      int size() const;
      bool isEquallyWeighted(double &weight) const;
      void clear();

      int minIRow() const;
      int minICol() const;
      int maxIRow() const;
      int maxICol() const;
      int nRows() const;
      int nCols() const;
      const CoordBR& getMBR() const;
      bool hasNbrs(MeshDir dir) const;
      const IntVect* nbrIDs(MeshDir dir) const;
      bool calcWorkBR(CoordBR &workBR,
                      const SpaceDims &dims) const;

      void toVect(vector<CellCoord> &vNbrCoords) const;
      void toVects(vector<CellCoord> &vNbrCoords,
                   vector<double> &vNbrWeights) const;
      bool fromVect(const vector<CellCoord> &vNbrCoords,
                    double weight = 1.0);
      bool fromVects(const vector<CellCoord> &vNbrCoords,
                     const vector<double> &vNbrWeights);

      bool specify(int cntrIRow, int cntrICol,
                   const CellSpace<elemType> &cellSpace);
      bool specify(const CellCoord &cntrCoord,
                   const CellSpace<elemType> &cellSpace);
      bool specify(const BasicCell<elemType> &cntrCell,
                   const CellSpace<elemType> &cellSpace);
      
      bool totalVal(elemType &totalVal, 
                    bool includingCtr = true) const;
      
      template<class Predicate>
      int count(Predicate pred, 
                bool includingCtr = true) const;
      
    private:
      bool _initCoords(const vector<CellCoord> &vNbrCoords,
                       double weight = 1.0);
      bool _initWeights(const vector<double> &vNbrWeights);

      bool _validNumDirects() const;

    private:
      vector< WeightedCell<elemType> > _vNbrs;
      CoordBR _MBR;
      vector<IntVect> _mNbrIDMap;
      bool _specified;
  };

  template<class elemType>
  ostream& operator<<(ostream &os,
                      const Neighborhood<elemType> &nbrhood);
  template <class elemType>
  istream& operator>>(istream &is,
                      Neighborhood<elemType> &nbrhood);
};

/****************************************
*             Private Methods           *
*****************************************/

template <class elemType>
inline bool pRPL::Neighborhood<elemType>::
_initCoords(const vector<CellCoord> &vNbrCoords,
            double weight) {
  if(vNbrCoords.empty()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: the coordinate vector is empty" << endl;
    return false;
  }
  if(find(vNbrCoords.begin(),
          vNbrCoords.end(),
          CellCoord(0, 0)) == vNbrCoords.end()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: the coordinate vector has to have"
         << " a coordinate [0, 0]" << endl;
    return false;
  }

  _vNbrs.erase(_vNbrs.begin(), _vNbrs.end());
  _mNbrIDMap.erase(_mNbrIDMap.begin(), _mNbrIDMap.end());
  _mNbrIDMap.resize(8);

  CellCoord nwCorner, seCorner;
  vector<CellCoord> addedCoords;
  for(int iNbr = 0; iNbr < vNbrCoords.size(); iNbr++) {
    if(find(addedCoords.begin(),
            addedCoords.end(),
            vNbrCoords[iNbr]) != addedCoords.end()) {
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: a duplicated coordinate [" \
           << vNbrCoords[iNbr] \
           << "] has been found" << endl;
      return false;
    }
    else {
      addedCoords.push_back(vNbrCoords[iNbr]);
    }

    _vNbrs.push_back(WeightedCell<elemType>(vNbrCoords[iNbr], weight));

    // calc _MBR
    if(0 == iNbr) {
      nwCorner = vNbrCoords[iNbr];
      seCorner = vNbrCoords[iNbr];
    }
    else{
      if(nwCorner.iRow() > vNbrCoords[iNbr].iRow()) {
         nwCorner.iRow(vNbrCoords[iNbr].iRow());
      }
      if(nwCorner.iCol() > vNbrCoords[iNbr].iCol()) {
         nwCorner.iCol(vNbrCoords[iNbr].iCol());
      }
      if(seCorner.iRow() < vNbrCoords[iNbr].iRow()) {
         seCorner.iRow(vNbrCoords[iNbr].iRow());
      }
      if(seCorner.iCol() < vNbrCoords[iNbr].iCol()) {
         seCorner.iCol(vNbrCoords[iNbr].iCol());
      }
    }

    // calc _mNbrIDMap
    if(vNbrCoords[iNbr].iRow() < 0) {
      if(vNbrCoords[iNbr].iCol() < 0) {
        _mNbrIDMap[NORTHWEST_DIR].push_back(iNbr);
      }
      else if(vNbrCoords[iNbr].iCol() > 0) {
        _mNbrIDMap[NORTHEAST_DIR].push_back(iNbr);
      }
      _mNbrIDMap[NORTH_DIR].push_back(iNbr);
    }
    if(vNbrCoords[iNbr].iRow() > 0) {
      if(vNbrCoords[iNbr].iCol() < 0) {
        _mNbrIDMap[SOUTHWEST_DIR].push_back(iNbr);
      }
      else if(vNbrCoords[iNbr].iCol() > 0) {
        _mNbrIDMap[SOUTHEAST_DIR].push_back(iNbr);
      }
      _mNbrIDMap[SOUTH_DIR].push_back(iNbr);
    }
    if(vNbrCoords[iNbr].iCol() < 0) {
      _mNbrIDMap[WEST_DIR].push_back(iNbr);
    }
    if(vNbrCoords[iNbr].iCol() > 0) {
      _mNbrIDMap[EAST_DIR].push_back(iNbr);
    }
  } // End of iNbr loop

  _MBR.nwCorner(nwCorner);
  _MBR.seCorner(seCorner);

  return true;
}

template <class elemType>
inline bool pRPL::Neighborhood<elemType>::
_initWeights(const vector<double> &vNbrWeights) {
  if(_vNbrs.size() != vNbrWeights.size()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: the number of coordinates does"
         << " NOT match the number of weights" << endl;
    return false;
  }
  for(int iNbr = 0; iNbr < vNbrWeights.size(); iNbr++) {
    _vNbrs[iNbr].weight(vNbrWeights[iNbr]);
  }
  return true;
}

template <class elemType>
inline bool pRPL::Neighborhood<elemType>::
_validNumDirects() const {
  bool numDirectsValid = true;
  if(_mNbrIDMap.size() != 8) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid _mNbrIDMap size" << endl;
    numDirectsValid = false;
  }
  return numDirectsValid;
}

/****************************************
*             Public Methods            *
*****************************************/

template <class elemType>
inline pRPL::Neighborhood<elemType>::
Neighborhood()
  :_specified(false){} 

template <class elemType>
inline pRPL::Neighborhood<elemType>::
Neighborhood(const vector<CellCoord> &vNbrCoords,
             double weight)
  :_specified(false) {
  if(!_initCoords(vNbrCoords, weight)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to construct a Neighborhood object" << endl;
    exit(-1);
  }
}

template <class elemType>
inline pRPL::Neighborhood<elemType>::
Neighborhood(const vector<CellCoord> &vNbrCoords,
             const vector<double> &vNbrWeights)
  :_specified(false) {
  if(!_initCoords(vNbrCoords) ||
     !_initWeights(vNbrWeights)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to construct a Neighborhood object" << endl;
    exit(-1);
  }
}

template <class elemType>
inline pRPL::Neighborhood<elemType>::
Neighborhood(const Neighborhood<elemType> &rhs)
  :_specified(false) {
  _vNbrs = rhs._vNbrs;
  _MBR = rhs._MBR;
  _mNbrIDMap = rhs._mNbrIDMap;
}

template<class elemType> template<class elemType2>
inline pRPL::Neighborhood<elemType>::
Neighborhood(const Neighborhood<elemType2> &rhs)
  :_specified(false) {
  vector<CellCoord> vCoords;
  vector<double> vWeights;
  rhs.toVects(vCoords, vWeights);
  if(!_initCoords(vCoords) ||
     !_initWeights(vWeights)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to construct a Neighborhood object" << endl;
    exit(-1);
  }
}

template <class elemType>
inline pRPL::WeightedCell<elemType>& pRPL::Neighborhood<elemType>::
operator[](int iNbr) {
  return _vNbrs.at(iNbr);  
}

template <class elemType>
inline const pRPL::WeightedCell<elemType>& pRPL::Neighborhood<elemType>::
operator[](int iNbr) const {
  return _vNbrs.at(iNbr);  
}

template <class elemType>
inline pRPL::Neighborhood<elemType>& pRPL::Neighborhood<elemType>::
operator=(const Neighborhood<elemType> &rhs) {
  if(this != &rhs) {
    _vNbrs = rhs._vNbrs;
    _MBR = rhs._MBR;
    _mNbrIDMap = rhs._mNbrIDMap;
    _specified = rhs._specified;
  }
  return *this;
}

template <class elemType>
inline bool pRPL::Neighborhood<elemType>::
empty() const{
  return _vNbrs.empty();
}

template <class elemType>
inline int pRPL::Neighborhood<elemType>::
size() const{
  return _vNbrs.size();  
}

template <class elemType>
bool pRPL::Neighborhood<elemType>::
isEquallyWeighted(double &weight) const {
  bool equal = true;
  if(_vNbrs.empty()) {
    return equal;
  }
  weight = _vNbrs[0].weight();
  for(int iNbr = 1; iNbr < _vNbrs.size(); iNbr++) {
    if(weight != _vNbrs[iNbr].weight()) {
      equal = false;
      break;
    }
  }
  return equal;
}

template <class elemType>
void pRPL::Neighborhood<elemType>::
clear() {
  _vNbrs.erase(_vNbrs.begin(), _vNbrs.end());
  _MBR = CoordBR();
  _mNbrIDMap.erase(_mNbrIDMap.begin(),
                   _mNbrIDMap.end());
  _specified = false;
}

template <class elemType>
inline int pRPL::Neighborhood<elemType>::
minIRow() const {
  return _MBR.minIRow();
}

template <class elemType>
inline int pRPL::Neighborhood<elemType>::
minICol() const {
  return _MBR.minICol();
}

template <class elemType>
inline int pRPL::Neighborhood<elemType>::
maxIRow() const {
  return _MBR.maxIRow();
}

template <class elemType>
inline int pRPL::Neighborhood<elemType>::
maxICol() const {
  return _MBR.maxICol();
}

template <class elemType>
inline int pRPL::Neighborhood<elemType>::
nRows() const {
  return _MBR.nRows();
}

template <class elemType>
inline int pRPL::Neighborhood<elemType>::
nCols() const {
  return _MBR.nCols();
}

template <class elemType>
inline const pRPL::CoordBR& pRPL::Neighborhood<elemType>::
getMBR() const {
  return _MBR;
}

template <class elemType>
inline bool pRPL::Neighborhood<elemType>::
hasNbrs(MeshDir dir) const {
  bool nbrExist;
  if(!_validNumDirects()) {
    nbrExist = false; 
  }
  else {
    if(!_mNbrIDMap[dir].empty()) {
      nbrExist = true;
    }
    else {
      nbrExist = false;
    }
  }
  return nbrExist;
}

template <class elemType>
inline const pRPL::IntVect* pRPL::Neighborhood<elemType>::
nbrIDs(MeshDir dir) const {
  return &(_mNbrIDMap.at(dir));
}

template <class elemType>
bool pRPL::Neighborhood<elemType>::
calcWorkBR(CoordBR &workBR,
           const SpaceDims &dims) const {
  if(empty()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to calculate workBR" \
         << " using an empty Neighborhood" << endl;
    return false;
  }
  if(dims.isNone()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to calculate workBR" \
         << " using an empty SpaceDims" << endl;
    return false;
  }
  if(dims.nRows() < nRows() ||
     dims.nCols() < nCols()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: the dimensions (" << dims << ") is smaller than" \
         << " the dimensions of the Neighborhood (" << nRows() \
         << " " << nCols() << ")" << endl;
    return false;
  }
  workBR.nwCorner(-minIRow(), -minICol());
  workBR.seCorner(dims.nRows() - maxIRow() - 1,
                  dims.nCols() - maxICol() - 1);
  if(!workBR.valid(dims)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid workBR (" << workBR \
         << ") has been produced" << endl;
    return false;
  }
  
  if((maxIRow() > 0 && workBR.nRows() < maxIRow()) ||
     (minIRow() < 0 && workBR.nRows() < -minIRow()) ||
     (maxICol() > 0 && workBR.nCols() < maxICol()) ||
     (minICol() < 0 && workBR.nCols() < -minICol())) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: workBR (" << workBR \
         << ") is too small to accomidate the Neighborhood" << endl;
    return false;
  }
  return true;
}

template <class elemType>
void pRPL::Neighborhood<elemType>::
toVect(vector<CellCoord> &vNbrCoords) const {
  vNbrCoords.erase(vNbrCoords.begin(), vNbrCoords.end());
  for(int iNbr = 0; iNbr < _vNbrs.size(); iNbr++) {
    vNbrCoords.push_back(_vNbrs[iNbr].coord());
  }
}

template <class elemType>
void pRPL::Neighborhood<elemType>::
toVects(vector<CellCoord> &vNbrCoords,
        vector<double> &vNbrWeights) const {
  vNbrCoords.erase(vNbrCoords.begin(), vNbrCoords.end());
  vNbrWeights.erase(vNbrWeights.begin(), vNbrWeights.end());

  for(int iNbr = 0; iNbr < _vNbrs.size(); iNbr++) {
    vNbrCoords.push_back(_vNbrs[iNbr].coord());
    vNbrWeights.push_back(_vNbrs[iNbr].weight());
  }
}

template <class elemType>
bool pRPL::Neighborhood<elemType>::
fromVect(const vector<CellCoord> &vNbrCoords,
         double weight) {
  _specified = false;
  if(!_initCoords(vNbrCoords, weight)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to initialize a Neighborhood object" \
         << " from a CellCoord vector" << endl;
    return false;
  }
  return true;
}

template <class elemType>
bool pRPL::Neighborhood<elemType>::
fromVects(const vector<CellCoord> &vNbrCoords,
          const vector<double> &vNbrWeights) {
  _specified = false;
  if(!_initCoords(vNbrCoords) ||
     !_initWeights(vNbrWeights)) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to initialize a Neighborhood object" \
         << " from a CellCoord vector and a weight vector" << endl;
    return false;
  }
  return true;
}

template <class elemType>
bool pRPL::Neighborhood<elemType>::
specify(int cntrIRow, int cntrICol,
        const CellSpace<elemType> &cellSpace) {
  if(empty()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to specify an empty Neighborhood" << endl;
    return false;
  }
  if(cellSpace.empty()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to specify a Neighborhood"
         << " on an empty CellSpace" << endl;
    return false;
  }

  if(cntrIRow < -minIRow() ||
     cntrIRow > cellSpace.nRows()-1-maxIRow() ||
     cntrICol < -minICol() ||
     cntrICol > cellSpace.nCols()-1-maxICol()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid center cell's coordinate[" \
         << cntrIRow << " " << cntrICol << "]" \
         << " on dimensions (" << cellSpace.dims() \
         << ")"<< endl;
    return false;
  }

  for(int iNbr = 0; iNbr < _vNbrs.size(); iNbr++) {
    int iRow = cntrIRow + _vNbrs[iNbr].iRow();
    int iCol = cntrICol + _vNbrs[iNbr].iCol();
    _vNbrs[iNbr].val(cellSpace[iRow][iCol]);
  }
  _specified = true;

  return true;
}

template <class elemType>
bool pRPL::Neighborhood<elemType>::
specify(const CellCoord &cntrCoord,
        const CellSpace<elemType> &cellSpace) {
  return specify(cntrCoord.iRow(),
                 cntrCoord.iCol(),
                 cellSpace);
}

template <class elemType>
bool pRPL::Neighborhood<elemType>::
specify(const BasicCell<elemType> &cntrCell,
        const CellSpace<elemType> &cellSpace) {
  return specify(cntrCell.iRow(),
                 cntrCell.iCol(),
                 cellSpace);
}

template <class elemType>
bool pRPL::Neighborhood<elemType>::
totalVal(elemType &totalVal,
         bool includingCtr) const {
  if(!_specified) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to calculate the total value"
         << " of an unspecified Neighborhood" << endl;
    return false;  
  }

  totalVal = static_cast<elemType>(0.0);
  for(int iNbr = 0; iNbr < _vNbrs.size(); iNbr++) {
    if(_vNbrs[iNbr].coord() == CellCoord(0, 0) &&
       !includingCtr) {
      continue;
    }
    totalVal += static_cast<elemType>(_vNbrs[iNbr].val() *
                                      _vNbrs[iNbr].weight());
  }

  return true;
}

template<class elemType> template<class Predicate>
int pRPL::Neighborhood<elemType>::
count(Predicate pred,
      bool includingCtr) const {
  if(!_specified) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to calculate the total value"
         << " of an unspecified Neighborhood" << endl;
    return 0;  
  }
  vector<elemType> vVals;
  for(int iNbr = 0; iNbr < _vNbrs.size(); iNbr++) {
    if(_vNbrs[iNbr].coord() == CellCoord(0, 0) &&
       !includingCtr) {
      continue;
    }
    vVals.push_back(_vNbrs[iNbr].val());
  }
  return std::count_if(vVals.begin(), vVals.end(), pred);
}

template<class elemType>
ostream& pRPL::
operator<<(ostream &os,
           const Neighborhood<elemType> &nbrhood) {
  int nNbrs = nbrhood.size();
  os << nNbrs << endl;
  for(int iNbr = 0; iNbr < nNbrs; iNbr++) {
    const WeightedCell<elemType> &nbr = nbrhood[iNbr];
    os << nbr.iRow() << "\t" \
       << nbr.iCol() << "\t" \
       << nbr.weight() << endl;
  }
  os << endl;
  return os;
}

template <class elemType>
istream& pRPL::
operator>>(istream &is,
           Neighborhood<elemType> &nbrhood) {
  vector<CellCoord> vNbrCoords;
  vector<double> vNbrWeights;
  int nNbrs;
  int iRow, iCol;
  double weight;
  
  is >> nNbrs;
  for(int iNbr = 0; iNbr < nNbrs; iNbr++) {
    is >> iRow;
    is >> iCol;
    is >> weight;
    vNbrCoords.push_back(CellCoord(iRow, iCol));
    vNbrWeights.push_back(weight);
  }
  nbrhood.fromVects(vNbrCoords, vNbrWeights);
  
  return is;
}

#endif
