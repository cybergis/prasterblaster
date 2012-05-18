

#ifndef QUADTREE_HH
#define QUADTREE_HH

#include <cstddef>
#include <vector>

#include "projectedraster.hh"

using std::vector;

class QuadTree 
{
public:
	QuadTree(Area boundry, size_t maximum_partition_size);
	QuadTree(size_t rows, size_t columns, size_t maximum_partition_size);
	~QuadTree();

	void subdivide();

	vector<Area> collectLeaves();

	QuadTree *northWest;
	QuadTree *northEast;
	QuadTree *southWest;
	QuadTree *southEast;

	size_t max_partition;
	Area boundry;
};


#endif // QUADTREE_HH
