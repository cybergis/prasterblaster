#ifndef ASPECTTRANSITION_H
#define ASPECTTRANSITION_H

/***************************************************************************
* aspectTransition.h
*
* Project: pRPL, v 1.0
* Purpose: Header file for class AspectTransition
* Description:
*          _cellSize - the length of the cell's border, assuming the cell 
*                      is a square. The cellsize for the test data, i.e., usa_nw.dem,
*                      is 0.008333 in geographic latitude-longitude degree
*          _scale -    the horizontal-to-vertical unit transformation scale. 
*                      The scale for the test data is 111120 (latitude-longitude 
*                      degrees to meters, 1 degree roughly equals to 111120 meters
*                      on the equator)
*         _noData -    the value indicating non-data within the data. The default
*                      value is set to be -32767
*
* Author:  Qingfeng (Gene) Guan
* E-mail:  guanqf {at} gmail.com
****************************************************************************
* Copyright (c) 2008, Qingfeng Guan
* 
****************************************************************************/

#include "cellSpace.h"
#include "neighborhood.h"
#include "subSpace.h"
#include "layer.h"
#include <cmath>
#include <functional>

using namespace pRPL;

class AspectTransition : public Transition<double> {
  public:
    AspectTransition()
      :Transition<double>(true, false, false, false),
       _cellSize(0.0),
       _scale(1.0),
       _noData(-32767),
       _pDEMLayer(0),
       _pDEMSubSpc(0),
       _pDEMNbrhood(0) {}
    AspectTransition(double cellSize, double scale)
      :Transition<double>(true, false, false, false),
       _cellSize(cellSize),
       _scale(scale),
       _noData(-32767),
       _pDEMLayer(0),
       _pDEMSubSpc(0),
       _pDEMNbrhood(0) {}
    ~AspectTransition() {}

    void cellSize(double size);
    void scale(double h2vScale);
    void demLayer(Layer<short> &layerD);

    virtual bool cellSpace(CellSpace<double> *pSubSpc);

    virtual bool evaluate(vector<pair<int, double> > &vUpdtedCells,
                          const CellCoord &coord);

    virtual int workload(const CoordBR &workBR);

  protected:
    double _cellSize;
    double _scale;
    short _noData;
    
    Layer<short> *_pDEMLayer;
    SubSpace<short> *_pDEMSubSpc;
    Neighborhood<short> *_pDEMNbrhood;
    
  private:
    static const double _PI = 3.14159265;
    static const double _EPSINON = 0.000001;
};

#endif
