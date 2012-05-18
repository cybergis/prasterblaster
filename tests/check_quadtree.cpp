
#include <gtest/gtest.h>

#include <memory>

#include "projectedraster.hh"
#include "sharedptr.hh"
#include "quadtree.hh"

TEST(QUADTREE, quadtree_init) 
{
	Area a(0, 99, 99, 0);
	QuadTree tree(a, 25*25);

	vector<Area> leaves = tree.collectLeaves();

	for (int i = 0; i < leaves.size(); ++i) {
		printf("Leaf %d: %f %f %f %f\n", i, 
		       leaves.at(i).ul.x, leaves.at(i).ul.y, 
		       leaves.at(i).lr.x, leaves.at(i).lr.y);

	}

	return;
}
