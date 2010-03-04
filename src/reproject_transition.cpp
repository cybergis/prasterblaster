

#include "reproject_transition.hh"

ReprojectTransition::ReprojectTransition()
	: _Transition<unsigned char>(true, true, true, true)
{

	return;
}

bool ReprojectTransition::inputLayer(Layer<unsigned char> *__inputLayer)
{
	if (__inputLayer == 0) {
		return false;
	}
	
	_inputLayer = __inputLayer;
	
	return true;
}
								     
bool ReprojectTransition::cellSpace(CellSpace<unsigned char> *cellspace)
{
	if (pCellSpc == 0 || _inputLayer == 0) {
		return false;
	}

	

	return true;
}



