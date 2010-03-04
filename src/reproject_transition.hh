/* 
   Programmer: David Mattli
*/

#ifndef REPROJECT_TRANSITION_HH
#define REPROJECT_TRANSITION_HH

#include "pRPL/cellSpace.h"
#include "pRPL/neighborhood.h"
#include "pRPL/subSpace.h"
#include "pRPL/layer.h"

using namespace pRPL;

class ReprojectTransition : public pRPL::Transition<unsigned char> 
{
public:
	ReprojectTransition();
	~ReprojectTransition();
	bool inputLayer(Layer<unsigned char> *__inputLayer);
	virtual bool cellSpace(CellSpace<unsigned char> *cellspace);
       

protected: 
	Layer<unsigned char> *_inputLayer;
	SubSpace<unsigned char> *inputCellSpace;

};

#endif // REPROJECT_TRANSITION_HH
