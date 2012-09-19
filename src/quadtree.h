//
// Copyright 0000 <Nobody>
// @file
// @author David Matthew Mattli <dmattli@usgs.gov>
//
// @section LICENSE
//
// This software is in the public domain, furnished "as is", without
// technical support, and with no warranty, express or implied, as to
// its usefulness for any purpose.
//
// @section DESCRIPTION
//
//
//

#ifndef SRC_QUADTREE_H_
#define SRC_QUADTREE_H_

#include <cstddef>
#include <vector>

#include "projectedraster.h"

using std::vector;

namespace librasterblaster {
struct QuadNode {
  QuadNode() {
    northWest = NULL;
    northEast = NULL;
    southWest = NULL;
    southEast = NULL;
    return;
  }

  explicit QuadNode(Area _boundry) {
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

  size_t depth;
};

class QuadTree {
 public:
  QuadTree(Area boundry, size_t maximum_partition_size);
  QuadTree(size_t rows, size_t columns, size_t maximum_partition_size);
  QuadTree(size_t rows, 
           size_t columns, 
           size_t maximum_partition_size, 
           size_t maximum_height, 
           size_t maximum_width);
  ~QuadTree();

  void subdivide();

  vector<Area> collectLeaves();

  QuadNode *rootNode;

  size_t max_partition;
  size_t max_height;
  size_t max_width;
};
}

#endif  // SRC_QUADTREE_H_
