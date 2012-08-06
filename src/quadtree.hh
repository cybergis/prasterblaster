

#ifndef QUADTREE_HH
#define QUADTREE_HH

#include <cstddef>
#include <vector>

#include "projectedraster.hh"

using std::vector;


struct QuadNode
{
	QuadNode() {
		northWest = NULL;
		northEast = NULL;
		southWest = NULL;
		southEast = NULL;
		return;
	}
	
	QuadNode(Area _boundry) {
		northWest = NULL;
		northEast = NULL;
		southWest = NULL;
		southEast = NULL;

		boundry = _boundry;
		return;
	}

	~QuadNode() {
		if (northWest != NULL) {
			delete northWest;
		}

		if (northEast != NULL) {
			delete northEast;
		}
		
		if (southWest != NULL) {
			delete southWest;
		}
		
		if (southEast != NULL) {
			delete southEast;
		}
		
	}
	
	Area boundry;

	QuadNode *northWest;
	QuadNode *northEast;
	QuadNode *southWest;
	QuadNode *southEast;
};

class QuadTree 
{
public:
	QuadTree(Area boundry, size_t maximum_partition_size);
	QuadTree(size_t rows, size_t columns, size_t maximum_partition_size);
	~QuadTree();

	void subdivide();

	vector<Area> collectLeaves();

	QuadNode *rootNode;

	size_t max_partition;
};


#endif // QUADTREE_HH
