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

#include "src/utils.h"

using std::vector;

namespace librasterblaster {
/**
 * \class QuadNode "quadtree.h"
 *
 * A struct representing a quadtree node with up to four children

 This struct is used with the QuadTree class to partition an area. It contains a
 member called 'boundry' that represents the inclusive area represented by the
 node. It also has four QuadTree* members: northWest, northEast, southWest,
 southEast. These are either NULL or point to child nodes. If all four QuadTree*
 pointers are NULL it mean the QuadNode is a child node.

 */
/* ! A constructor

   This constructor creates an empty node with four NULL pointers for children.
  
 */
struct QuadNode {
  QuadNode() {
    northWest = NULL;
    northEast = NULL;
    southWest = NULL;
    southEast = NULL;
    return;
  }

  /* ! A constructor
     
    This constructor takes an Area argument that represents its inclusive
    boundry. The four QuadTree* pointers are set to NULL.
   */
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

/**
 * \class QuadTree "quadtree.h"
 
 This class partitions a space using a quadtree.

 */
class QuadTree {
 public:
  /*! A constructor 

     This contructor initializes a QuadTree with a given boundry and maximum
     partition size.

     @param boundry This parameter is an inclusive Area that is partitioned.

     @param maximum_partition_size This size_t specifies the maximum size of the
     generated partitions. The partitions may be smaller than specified.

   */
  QuadTree(Area boundry, size_t maximum_partition_size);

  /*! A constructor
     
     This constructor initializes a QuadTree that partitions an area equal to
     the given number of rows and columns. 

     @param rows Number of rows in the area to be partitioned

     @param columns Number of columns in the area to be partitioned

     @param maximum_partition_size Maximum size of the partitions to be created

   */
  QuadTree(size_t rows, size_t columns, size_t maximum_partition_size);
  ~QuadTree();
  /*! A normal member function

     This function returns a vector<Area> that lists the child nodes of the quad
     tree. These Areas partition the boundry specified in the constructor.

   */
  vector<Area> collectLeaves();

 private:
  void subdivide();

  QuadNode *rootNode;

  size_t max_partition;
};
}

#endif  // SRC_QUADTREE_H_
