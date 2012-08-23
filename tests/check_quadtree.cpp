
#include <gtest/gtest.h>

#include <memory>
#include <cstdlib>

#include "projectedraster.hh"
#include "sharedptr.hh"
#include "quadtree.hh"

TEST(QUADTREE, quadtree_init) 
{
	Area a(0, 99, 99, 0);
	QuadTree tree(a, 50*50);

	vector<Area> leaves = tree.collectLeaves();
	printf("%zd\n", leaves.size());
	for (int i = 0; i < leaves.size(); ++i) {
		Area t = leaves.at(i);
		printf("Node: %d {%f %f %f %f} \n",
		       i,
		       t.ul.x, t.ul.y, t.lr.x, t.lr.y);
	}

	return;
}

TEST(QUADTREE, quadtree_big)
{

	QuadTree tree(161190, 104424, 50*50);

	vector<Area> leaves = tree.collectLeaves();
	printf("%zd leaves \n", leaves.size());
	return;
}
