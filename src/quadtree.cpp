

#include "quadtree.hh"

QuadTree::QuadTree(Area _boundry, size_t  maximum_partition_size)
{
	max_partition = maximum_partition_size;
	boundry = _boundry;

	northWest = NULL;
	northEast = NULL;
	southWest = NULL;
	southEast = NULL;
	
	subdivide();

	return;
}

QuadTree::QuadTree(size_t rows, size_t columns, size_t maximum_partition_size)
{
	max_partition = maximum_partition_size;
	boundry = Area(0, rows - 1, columns - 1, 0);
	northWest = NULL;
	northEast = NULL;
	southWest = NULL;
	southEast = NULL;
	
	subdivide();
	return;
}

QuadTree::~QuadTree()
{
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

void QuadTree::subdivide()
{
	size_t width = boundry.lr.x - boundry.ul.x + 1;
	size_t height = boundry.ul.y - boundry.lr.y + 1;
	size_t area_size = width * height;

	const int middle_width = ((boundry.lr.x - boundry.ul.x) / 2);
	const int middle_height = ((boundry.ul.y - boundry.lr.y) / 2);

	const int east_begin = boundry.ul.x + middle_width + 1;
	const int east_end = boundry.lr.x;
	const int west_begin = boundry.ul.x;
	const int west_end = boundry.ul.x + middle_width;
	const int north_begin =  boundry.lr.y + middle_height + 1;
	const int north_end = boundry.ul.y;
	const int south_begin = boundry.lr.y;
	const int south_end = boundry.lr.y + middle_height;


	if (area_size > max_partition)
	{
		if (width > 1) {
			
			southWest = new QuadTree(Area(west_begin,
						      south_end,
						      west_end,
						      south_begin),
						 max_partition);
			
			southEast = new QuadTree(Area(east_begin,
						      south_end,
						      east_end,
						      south_begin),
						 max_partition);
		
		}
		
		if (height > 1) {
			northWest = new QuadTree(Area(west_begin,
						      north_end,
						      west_end,
						      north_begin),
						 max_partition);
			
			northEast = new QuadTree(Area(east_begin,
						      north_end,
						      east_end,
						      north_begin),
						 max_partition);
		}
		
	}
	
	return;
}

vector<Area> QuadTree::collectLeaves()
{
	vector<Area> r;
	vector<Area> tmp;

	if ((northWest == NULL) 
	    && (northEast == NULL)
	    && (southWest == NULL)
	    && (southEast == NULL)) {
		r.push_back(boundry);
		return r;
	}

	if (northWest != NULL) {
		tmp = northWest->collectLeaves();
		r.insert(r.end(), tmp.begin(), tmp.end());
	}
	
	if (northEast != NULL) {
		tmp = northEast->collectLeaves();
		r.insert(r.end(), tmp.begin(), tmp.end());
	}
	
	if (southWest != NULL) {
		tmp = southWest->collectLeaves();
		r.insert(r.end(), tmp.begin(), tmp.end());
	}
	
	if (southEast != NULL) {
		tmp = southEast->collectLeaves();
		r.insert(r.end(), tmp.begin(), tmp.end());
	}

	return r;
}
