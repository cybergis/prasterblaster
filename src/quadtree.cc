

#include "quadtree.hh"

QuadTree::QuadTree(Area _boundry, size_t  maximum_partition_size)
{
	max_partition = maximum_partition_size;

	rootNode = new QuadNode();
	rootNode->boundry = _boundry;
	rootNode->depth = 0;

	subdivide();

	return;
}

QuadTree::QuadTree(size_t rows, size_t columns, size_t maximum_partition_size)
{
	max_partition = maximum_partition_size;

	rootNode = new QuadNode();
	rootNode->boundry = Area(0, rows - 1, columns - 1, 0);
	rootNode->depth = 0;

	subdivide();
	return;
}

QuadTree::~QuadTree()
{
	if (rootNode != NULL) {
		delete rootNode;
	}
}
size_t getWidth(QuadNode *node) 
{
	return node->boundry.lr.x - node->boundry.ul.x + 1;
}

size_t getHeight(QuadNode *node)
{
	return node->boundry.ul.y - node->boundry.lr.y + 1;
}

void QuadTree::subdivide()
{
	int middle_width = 0;
	int middle_height = 0;

	int east_begin = 0;
	int east_end = 0;
	int west_begin = 0;
	int west_end = 0;
	int north_begin = 0;
	int north_end = 0;
	int south_begin = 0;
	int south_end = 0;

	vector<QuadNode*> stack;
	stack.push_back(rootNode);

	while(stack.size() > 0) {
		// Pop element off stack
		QuadNode *n = stack.back();
		stack.pop_back();
		Area boundry = n->boundry;

		// Calculate parameters of QuadNode from top of stack
		middle_width = ((boundry.lr.x - boundry.ul.x) / 2);
		middle_height = ((boundry.ul.y - boundry.lr.y ) / 2);

		east_begin = boundry.ul.x + middle_width + 1;
		east_end = boundry.lr.x;
		west_begin = boundry.ul.x;
		west_end = boundry.ul.x + middle_width;
		north_begin =  boundry.lr.y + middle_height + 1;
		north_end = boundry.ul.y;
		south_begin = boundry.lr.y;
		south_end = boundry.lr.y + middle_height;

		

		size_t area_size = getWidth(n) * getHeight(n);
		if (area_size > max_partition) {
			if (getWidth(n) > 1) {
				n->southWest = new QuadNode(Area(west_begin,
								 south_end,
								 west_end,
								 south_begin));
				n->southWest->depth = n->depth + 1;
				n->southEast = new QuadNode(Area(east_begin,
								 south_end,
								 east_end,
								 south_begin));
				n->southEast->depth = n->depth + 1;

				stack.push_back(n->southWest);
				stack.push_back(n->southEast);
			}

			if (getHeight(n) > 1) {
				n->northWest = new QuadNode(Area(west_begin,
								 north_end,
								 west_end,
								 north_begin));
				n->northWest->depth = n->depth + 1;
				n->northEast = new QuadNode(Area(east_begin,
								 north_end,
								 east_end,
								 north_begin));
				n->northEast->depth = n->depth + 1;
				stack.push_back(n->northWest);
				stack.push_back(n->northEast);

			}
		}
	}

	return;
}

vector<Area> QuadTree::collectLeaves()
{
	vector<Area> r;
	vector<QuadNode*> s;

	if (rootNode == NULL) {
		return r;
	}

	s.push_back(rootNode);

	while(s.size() > 0) {
		QuadNode *n = s.back();
		s.pop_back();

		if (n->northWest != NULL) {
			s.push_back(n->northWest);
		}

		if (n->northEast != NULL) {
			s.push_back(n->northEast);
		}

		if (n->southWest != NULL) {
			s.push_back(n->southWest);
		}

		if (n->southEast != NULL) {
			s.push_back(n->southEast);
		}

		// Check if we are a leaf node
		if (n->northWest == NULL &&
		    n->northEast == NULL &&
		    n->southWest == NULL &&
		    n->southEast == NULL) {
			// We are a leaf, push the boundry onto the stack
			r.push_back(n->boundry);
		}
	}

	return r;
}
