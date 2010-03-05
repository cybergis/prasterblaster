
#include <vector>

#include "projection.h"

class ProjectionFactory 
{
public:
  ProjectionFactory();
  ~ProjectionFactory();
  Projection fromProj4Specification(vector<string> args);


};
