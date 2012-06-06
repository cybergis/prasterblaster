
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

	}

	return;
}
