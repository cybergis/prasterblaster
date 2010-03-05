#include "aspectTransition.h"

/***************************************************************************
* slopeTransition.cpp
*
* Project: pRPL, v 1.0
* Purpose: Implementation for class AspectTransition
* Author:  Qingfeng (Gene) Guan
* E-mail:  guanqf {at} gmail.com
****************************************************************************
* Copyright (c) 2008, Qingfeng Guan
* 
****************************************************************************/

const double AspectTransition::_PI;
const double AspectTransition::_EPSINON;

void AspectTransition::
cellSize(double size) {
  _cellSize = size;
}

void AspectTransition::
scale(double h2vScale) {
  _scale = h2vScale;
}

void AspectTransition::
demLayer(Layer<short> &layerD) {
  _pDEMLayer = &layerD;
  _pDEMNbrhood = layerD.nbrhood();
}

bool AspectTransition::
cellSpace(CellSpace<double> *pSubSpc) {
  Transition<double>::_pCellSpace = pSubSpc;
  if(!pSubSpc) {
    return false;
  }
  if(_pDEMLayer) {
    Layer<short> &layerD = *(_pDEMLayer);
    int subID = ((SubSpace<double> *)pSubSpc)->id();
    SubSpace<short> *pSubSpcD = layerD.findSubSpc(subID);
    if(pSubSpcD) {
      _pDEMSubSpc = pSubSpcD;
    }
  }
  return true;
}

bool AspectTransition::
evaluate(vector<pair<int, double> > &vUpdtedCells,
         const CellCoord &coord) {
  if(_cellSize <= 0.0 || _scale <= 0.0) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid cellsize (" << _cellSize \
         << "), or invalid scale (" << _scale << ")" \
         << endl;
    return false;
  }
  if(!Transition<double>::_pCellSpace ||
      Transition<double>::_pCellSpace->empty()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: no CellSpace or an empty CellSpace" \
         << " to operate on" << endl;
    return false;
  }
  CellSpace<double> &cellSpc = *(Transition<double>::_pCellSpace);

  if(!_pDEMSubSpc ||
      _pDEMSubSpc->empty()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: no extra CellSpace or an empty extra CellSpace" \
         << " to operate on" << endl;
    return false;
  }
  CellSpace<short> &cellSpcD = *(_pDEMSubSpc);

  if(!_pDEMNbrhood ||
     _pDEMNbrhood->empty()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: no extra Neighborhood or an empty extra Neighborhood" \
         << " to operate on" << endl;
    return false;
  }
  Neighborhood<short> &nbrhoodD = *(_pDEMNbrhood);

  if(!nbrhoodD.specify(coord, cellSpcD)) {
    return false;
  }
  
  double aspect;
  if(nbrhoodD.count(bind2nd(equal_to<short>(), _noData)) > 0) {
    aspect = 0.0;
  }
  else {
    double cellSize = _cellSize * _scale;
    double demDiff;
    demDiff = static_cast<double>((nbrhoodD[1].val() + nbrhoodD[4].val() + nbrhoodD[4].val()+ nbrhoodD[6].val()) - 
                                  (nbrhoodD[3].val() + nbrhoodD[5].val() + nbrhoodD[5].val()+ nbrhoodD[8].val()));
    double az = demDiff / (2.0 * cellSize);
    demDiff = static_cast<double>((nbrhoodD[6].val() + nbrhoodD[7].val() + nbrhoodD[7].val()+ nbrhoodD[8].val()) - 
                                  (nbrhoodD[1].val() + nbrhoodD[2].val() + nbrhoodD[2].val()+ nbrhoodD[3].val()));
    double bz = demDiff / (2.0 * cellSize);
    double slope = pow((pow(az, 2.0) + pow(bz, 2.0)), 0.5);
    slope = atan(slope) * 180.0 / _PI;
  
    /* Determine the quadrant & calculate aspect */
    if(az < -_EPSINON) {
      if(bz < -_EPSINON) {
        aspect = (_PI / 2.0) - atan(fabs((bz / az)));
      }
      else {
        aspect = (_PI / 2.0) + atan(fabs((bz / az)));
      }
    } /* End of if(az < 0.0) */
    else {
      if(az > _EPSINON) {
        if(bz < -_EPSINON) {
          aspect = (3.0 * _PI / 2.0) + atan(fabs((bz / az)));
        }
        else {
          aspect = (3.0 * _PI / 2.0) - atan(fabs((bz / az)));
        }
      } /* End of if(az > 0.0) */
      else {
        /* az equals zero */
        if (bz < -_EPSINON) {
          aspect = 0.0;
        }
        else {
          aspect = _PI;
        }
      } /* End of if(az == 0.0) */
    } /* End of if(az >= 0.0) */

    /* Convert from radians to degrees */
    aspect = 360.0 * (aspect / (2.0 * _PI));
    /* Correct to north */
    aspect += 180.0;
    if(aspect >= 360.0) {
      aspect -= 360.0;
    }
    
    if(slope >= -_EPSINON &&
       slope <= _EPSINON) {
      aspect = 0.0;
    }
    else {
      aspect = 1.0 + aspect;
    }
  }
  
  int idx = cellSpc.coord2idx(coord);
  vUpdtedCells.push_back(make_pair(idx, aspect));

  return true;
}

int AspectTransition::
workload(const CoordBR &workBR) {
  int wLoad = 0;
  return wLoad;
}
