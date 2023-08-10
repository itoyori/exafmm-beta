#ifndef build_tree_tbb_h
#define build_tree_tbb_h
#include "logger.h"
#include "thread.h"
#include "types.h"

namespace EXAFMM_NAMESPACE {
  class BuildTree {
    template <typename T>
    using global_ptr = ityr::ori::global_ptr<T>;

    typedef vec<8,int> ivec8;                                   //!< Vector of 8 integer types

  public:
    //! Binary tree is used for counting number of bodies with a recursive approach
    struct BinaryTreeNode {
      ivec8                      NBODY;                                   //!< Number of descendant bodies
      global_ptr<BinaryTreeNode> LEFT;                                    //!< Pointer to left child
      global_ptr<BinaryTreeNode> RIGHT;                                   //!< Pointer to right child
      global_ptr<BinaryTreeNode> BEGIN;                                   //!< Pointer to begining of memory space
      global_ptr<BinaryTreeNode> END;                                     //!< Pointer to end of memory space
    };

  private:
    //! Octree is used for building the FMM tree structure as "nodes", then transformed to "cells" data structure
    struct OctreeNode {
      int                    IBODY;                                       //!< Index offset for first body in node
      int                    NBODY;                                       //!< Number of descendant bodies
      int                    NNODE;                                       //!< Number of descendant nodes
      global_ptr<OctreeNode> CHILD[8];                                    //!< Pointer to child node
      vec3                   X;                                           //!< Coordinate at center
    };

    const int              ncrit;                                         //!< Number of bodies per leaf cell
    const int              nspawn;                                        //!< Threshold of NBODY for spawning new threads
    int                    numLevels;                                     //!< Number of levels in tree
    GB_iter                B0;                                            //!< Iterator of first body
    global_ptr<OctreeNode> N0;                                            //!< Pointer to octree root node

  private:
    static global_ptr<OctreeNode> alloc_octree_node() {
      return ityr::ori::malloc<OctreeNode>(1);
    }

    static void free_octree_node(global_ptr<OctreeNode> node) {
      ityr::ori::free(node, 1);
    }

    //! Recursive functor for counting bodies in each octant using binary tree
    struct CountBodies {
      GBodies bodies;                                          //!< Vector of bodies
      int begin;                                                //!< Body begin index
      int end;                                                  //!< Body end index
      vec3 X;                                                   //!< Coordinate of node center
      global_ptr<BinaryTreeNode> binNode;                                 //!< Pointer to binary tree node
      int nspawn;                                               //!< Threshold of NBODY for spawning new threads
      CountBodies(GBodies _bodies, int _begin, int _end, vec3 _X,// Constructor
		  global_ptr<BinaryTreeNode> _binNode, int _nspawn) :
	bodies(_bodies), begin(_begin), end(_end), X(_X),       // Initialize variables
	binNode(_binNode), nspawn(_nspawn) {}
      //! Get number of binary tree nodes for a given number of bodies
      inline int getNumBinNode(int n) const {
	if (n <= nspawn) return 1;                              // If less then threshold, use only one node
	else return 4 * ((n - 1) / nspawn) - 1;                 // Else estimate number of binary tree nodes
      }
      void operator() () const {                                // Overload operator()
	/* assert(getNumBinNode(end - begin) <= binNode->END - binNode->BEGIN + 1); */
        auto cs = ityr::make_checkout(binNode, 1, ityr::checkout_mode::read_write);
        auto binNode_ = &cs[0];

	if (end - begin <= nspawn) {                            //  If number of bodies is less than threshold
          for (int i=0; i<8; i++) binNode_->NBODY[i] = 0;        //   Initialize number of bodies in octant
          binNode_->LEFT = binNode_->RIGHT = nullptr;                //   Initialize pointers to left and right child node

          ityr::for_each(
              body_seq_policy,
              ityr::count_iterator<int>(begin),
              ityr::count_iterator<int>(end),
              ityr::make_global_iterator(bodies.begin() + begin, ityr::checkout_mode::read),
              [&](int i, const auto& B) {
                vec3 x = B.X;                                   //  Coordinates of body
                if (B.ICELL < 0) {                                //  If using residual index
                  auto mp_X = static_cast<vec3 Body::*>(&Source::X);
                  x = (&bodies[i+B.ICELL])->*(mp_X);                      //   Use coordinates of first body in residual group
                }
                int octant = (x[0] > X[0]) + ((x[1] > X[1]) << 1) + ((x[2] > X[2]) << 2);// Which octant body belongs to
                binNode_->NBODY[octant]++;                           //    Increment body count in octant
              });

	} else {                                                //  Else if number of bodies is larger than threshold
	  int mid = (begin + end) / 2;                          //   Split range of bodies in half
	  while ((&bodies[mid])->*(&Body::ICELL) < 0) mid++;                  //   Don't split residual groups
	  int numLeftNode = getNumBinNode(mid - begin);         //   Number of binary tree nodes on left branch
	  int numRightNode = getNumBinNode(end - mid);          //   Number of binary tree nodes on right branch
	  /* assert(numLeftNode + numRightNode <= binNode->END - binNode->BEGIN);// Bounds checking for node count */

          global_ptr<BinaryTreeNode> binNode_left, binNode_right;

          binNode_left = binNode_->BEGIN;
          binNode_right = binNode_->BEGIN + numLeftNode;
          binNode_->LEFT = binNode_left;
          binNode_->RIGHT = binNode_right;

          auto [csl, csr] = ityr::make_checkouts(binNode_left , 1, ityr::checkout_mode::write,
                                                 binNode_right, 1, ityr::checkout_mode::write);
          auto binNode_left_  = &csl[0];
          auto binNode_right_ = &csr[0];

          binNode_left_->BEGIN = binNode_left + 1;             //   Assign next memory address to left begin pointer
          binNode_left_->END = binNode_left + numLeftNode;     //   Keep track of last memory address used by left
          binNode_right_->BEGIN = binNode_right + 1;           //   Assign next memory address to right begin pointer
          binNode_right_->END = binNode_right + numRightNode;  //   Keep track of last memory address used by right

          CountBodies leftBranch(bodies, begin, mid, X, binNode_left, nspawn);// Recursion for left branch
          CountBodies rightBranch(bodies, mid, end, X, binNode_right, nspawn);// Recursion for right branch

          cs.checkin();
          csl.checkin();
          csr.checkin();

          ityr::parallel_invoke(
            [=]() { leftBranch(); },
            [=]() { rightBranch(); }
          );

	  binNode->*(&BinaryTreeNode::NBODY) =
            ivec8(binNode_left->*(&BinaryTreeNode::NBODY)) + ivec8(binNode_right->*(&BinaryTreeNode::NBODY));// Sum contribution from both branches
	}                                                       //  End if for number of bodies
      }                                                         // End overload operator()
    };

    //! Recursive functor for sorting bodies according to octant (Morton order)
    struct MoveBodies {
      GBodies bodies;                                          //!< Vector of bodies
      GBodies buffer;                                          //!< Buffer for bodies
      int begin;                                                //!< Body begin index
      int end;                                                  //!< Body end index
      global_ptr<BinaryTreeNode> binNode;                                 //!< Pointer to binary tree node
      mutable ivec8 octantOffset;                               //!< Offset of octant
      vec3 X;                                                   //!< Coordinates of node center
      MoveBodies(GBodies _bodies, GBodies _buffer, int _begin, int _end,// Constructor
		 global_ptr<BinaryTreeNode> _binNode, ivec8 _octantOffset, vec3 _X) :
	bodies(_bodies), buffer(_buffer), begin(_begin), end(_end),// Initialize variables
	binNode(_binNode), octantOffset(_octantOffset), X(_X) {}
      void operator() () const {                                // Overload operator()
        global_ptr<BinaryTreeNode> binNode_left = binNode->*(&BinaryTreeNode::LEFT);
	if (binNode_left == nullptr) {                            //  If there are no more child nodes
          ivec8 counts = binNode->*(&BinaryTreeNode::NBODY);
          Body* buffers[8];
          ivec8 offsets;
          for (int i = 0; i < 8; i++) {
            if (counts[i] > 0) {
              buffers[i] = ityr::ori::checkout(
                  &buffer[octantOffset[i]], counts[i], ityr::checkout_mode::write);
              offsets[i] = 0;
            }
          }

          ityr::for_each(
              body_seq_policy,
              ityr::count_iterator<int>(begin),
              ityr::count_iterator<int>(end),
              ityr::make_global_iterator(bodies.begin() + begin, ityr::checkout_mode::read),
              [&](int i, const auto& B) {
                vec3 x = B.X;                                   //  Coordinates of body
                if (B.ICELL < 0) {                                //  If using residual index
                  auto mp_X = static_cast<vec3 Body::*>(&Source::X);
                  x = (&bodies[i+B.ICELL])->*(mp_X);                      //   Use coordinates of first body in residual group
                }
                int octant = (x[0] > X[0]) + ((x[1] > X[1]) << 1) + ((x[2] > X[2]) << 2);// Which octant body belongs to`
                /* buffer[octantOffset[octant]] = B;               //   Permute bodies out-of-place according to octant */
                /* octantOffset[octant]++;                                 //  Increment body count in octant */
                buffers[octant][offsets[octant]] = B;
                offsets[octant]++;
              });

          for (int i = 0; i < 8; i++) {
            if (counts[i] > 0) {
              ityr::ori::checkin(buffers[i], counts[i], ityr::checkout_mode::write);
            }
          }
	} else {                                                //  Else if there are child nodes
	  int mid = (begin + end) / 2;                          //   Split range of bodies in half
	  while ((&bodies[mid])->*(&Body::ICELL) < 0) mid++;                  //   Don't split residual groups

          MoveBodies leftBranch(bodies, buffer, begin, mid, binNode_left, octantOffset, X);// Recursion for left branch

          global_ptr<BinaryTreeNode> binNode_right = binNode->*(&BinaryTreeNode::RIGHT);
          ivec8 octantOffset_ = octantOffset + binNode_left->*(&BinaryTreeNode::NBODY);
          MoveBodies rightBranch(bodies, buffer, mid, end, binNode_right, octantOffset_, X);// Recursion for right branch

          ityr::parallel_invoke(
            [=]() { leftBranch(); },
            [=]() { rightBranch(); }
          );
	}                                                       //  End if for child existance
      }                                                         // End overload operator()
    };

    //! Recursive functor for building nodes of an octree adaptively using a top-down approach
    struct BuildNodes {
      GBodies bodies;                                          //!< Vector of bodies
      GBodies buffer;                                          //!< Buffer for bodies
      int begin;                                                //!< Body begin index
      int end;                                                  //!< Body end index
      global_ptr<BinaryTreeNode> binNode;                                 //!< Pointer to binary tree node
      vec3 X;                                                   //!< Coordinate of node center
      real_t R0;                                                //!< Radius of root cell
      int ncrit;                                                //!< Number of bodies per leaf cell
      int nspawn;                                               //!< Threshold of NBODY for spawning new threads
      int level;                                                //!< Current tree level
      bool direction;                                           //!< Direction of buffer copying
      //! Constructor
      BuildNodes(GBodies _bodies,
		 GBodies _buffer, int _begin, int _end, global_ptr<BinaryTreeNode> _binNode,
		 vec3 _X, real_t _R0, int _ncrit, int _nspawn, int _level=0, bool _direction=false) :
	bodies(_bodies), buffer(_buffer),    // Initialize variables
	begin(_begin), end(_end), binNode(_binNode), X(_X), R0(_R0),
	ncrit(_ncrit), nspawn(_nspawn), level(_level), direction(_direction) {}

      //! Create an octree node
      global_ptr<OctreeNode> makeOctNode(bool nochild) const {
        global_ptr<OctreeNode> octNode = alloc_octree_node();
        auto cs = ityr::make_checkout(octNode, 1, ityr::checkout_mode::write);
        OctreeNode* octNode_ = new (&cs[0]) OctreeNode();                             // Allocate memory for single node
        octNode_->IBODY = begin;                                 // Index of first body in node
        octNode_->NBODY = end - begin;                           // Number of bodies in node
        octNode_->NNODE = 1;                                     // Initialize counter for decendant nodes
        octNode_->X = X;                                         // Center coordinates of node
        if (nochild) {                                          // If node has no children
          for (int i=0; i<8; i++) octNode_->CHILD[i] = nullptr;     //  Initialize pointers to children
        }                                                       // End if for node children
	return octNode;                                         // Return node
      }

      //! Exclusive scan with offset
      inline ivec8 exclusiveScan(ivec8 input, int offset) const {
	ivec8 output;                                           // Output vector
	for (int i=0; i<8; i++) {                               // Loop over elements
	  output[i] = offset;                                   //  Set value
	  offset += input[i];                                   //  Increment offset
	}                                                       // End loop over elements
	return output;                                          // Return output vector
      }

      //! Get maximum number of binary tree nodes for a given number of bodies
      inline int getMaxBinNode(int n) const {
	return (4 * n) / nspawn;                                // Conservative estimate of number of binary tree nodes
      }

      global_ptr<OctreeNode> operator() () {                                // Overload operator()
	/* double tic = logger::get_time(); */
	/* assert(getMaxBinNode(end - begin) <= binNode->END - binNode->BEGIN);// Bounds checking for node range */
	if (begin == end) {                                     //  If no bodies are left
	  return nullptr;                                               //   End buildNodes()
	}                                                       //  End if for no bodies
	if (end - begin <= ncrit) {                             //  If number of bodies is less than threshold
          if (direction) {                                          //  If direction of data is from bodies to buffer
            auto [b_src, b_dest] =
              ityr::make_checkouts(bodies.begin() + begin, end - begin, ityr::checkout_mode::read,
                                   buffer.begin() + begin, end - begin, ityr::checkout_mode::write);
            for (int i = 0; i < end - begin; i++) {
              b_dest[i] = b_src[i];
            }
          }
          return makeOctNode(true);        //  Create an octree node and assign it's pointer
	}                                                       //  End if for number of bodies
	global_ptr<OctreeNode> octNode = makeOctNode(false);                           //  Create an octree node with child nodes
	/* double toc = logger::get_time(); */
	/* timer["Make node"] += toc - tic; */
	CountBodies countBodies(bodies, begin, end, X, binNode, nspawn);// Instantiate recursive functor
	countBodies();                                          //  Count bodies in each octant using binary recursion
	/* tic = logger::get_time(); */
	/* timer["Count bodies"] += tic - toc; */
        ivec8 counts = binNode->*(&BinaryTreeNode::NBODY);
	ivec8 octantOffset = exclusiveScan(counts, begin);//  Exclusive scan to obtain offset from octant count
	/* toc = logger::get_time(); */
	/* timer["Exclusive scan"] += toc - tic; */ 
	MoveBodies moveBodies(bodies, buffer, begin, end, binNode, octantOffset, X);// Instantiate recursive functor
	moveBodies();                                           //  Sort bodies according to octant
	/* tic = logger::get_time(); */
	/* timer["Move bodies"] += tic - toc; */
	global_ptr<BinaryTreeNode> binNodeOffset = binNode->*(&BinaryTreeNode::BEGIN);        //  Initialize pointer offset for binary tree nodes

        global_ptr<BinaryTreeNode> binNodeOffsets[8];
        for (int i = 0; i < 8; i++) {
	  int maxBinNode = getMaxBinNode(counts[i]);    //   Get maximum number of binary tree nodes
          binNodeOffsets[i] = binNodeOffset;
	  binNodeOffset += maxBinNode;                          //   Increment offset for binNode memory address
        }

        auto children = global_ptr<global_ptr<OctreeNode>>(&(octNode->*(&OctreeNode::CHILD)));
        global_vec<BinaryTreeNode> binNodeChild_vec(8);              //  Allocate new root for this branch
        global_ptr<BinaryTreeNode> binNodeChild = binNodeChild_vec.begin();
        ityr::for_each(
            ityr::execution::par,
            ityr::count_iterator<int>(0),
            ityr::count_iterator<int>(8),
            [=, *this](int i) {
              /* toc = logger::get_time(); */
              int maxBinNode = getMaxBinNode(counts[i]);    //   Get maximum number of binary tree nodes
              /* assert(binNodeOffset + maxBinNode <= binNode->END);   //    Bounds checking for node count */
              vec3 Xchild = X;                                      //    Initialize center coordinates of child node
              real_t r = R0 / (1 << (level + 1));                   //    Radius of cells for child's level
              for (int d=0; d<3; d++) {                             //    Loop over dimensions
                Xchild[d] += r * (((i & 1 << d) >> d) * 2 - 1);     //     Shift center coordinates to that of child node
              }                                                     //    End loop over dimensions
              (&binNodeChild[i])->*(&BinaryTreeNode::BEGIN) = binNodeOffsets[i];                //    Assign first memory address from offset
              (&binNodeChild[i])->*(&BinaryTreeNode::END) = binNodeOffsets[i] + maxBinNode;     //    Keep track of last memory address
              /* tic = logger::get_time(); */
              /* timer["Get node range"] += tic - toc; */
              BuildNodes buildNodes(buffer, bodies,//    Instantiate recursive functor
            		            octantOffset[i], octantOffset[i] + counts[i],
            		            &binNodeChild[i], Xchild, R0, ncrit, nspawn, level+1, !direction);
              children[i] = buildNodes();
            });

	for (int i=0; i<8; i++) {                               //  Loop over children
          if (global_ptr<OctreeNode>(children[i]))
            octNode->*(&OctreeNode::NNODE) += global_ptr<OctreeNode>(children[i])->*(&OctreeNode::NNODE);// If child exists increment child node count
	}                                                       //  End loop over chlidren
        return octNode;
      }                                                         // End overload operator()
    };

    //! Recursive functor for creating cell data structure from nodes
    struct Nodes2cells {
      global_ptr<OctreeNode> octNode;                                     //!< Pointer to octree node
      GB_iter B0;                                                //!< Iterator of first body
      GC_iter C;                                                 //!< Iterator of current cell
      GC_iter C0;                                                //!< Iterator of first cell
      mutable GC_iter CN;                                        //!< Iterator of cell counter
      vec3 X0;                                                  //!< Coordinate of root cell center
      real_t R0;                                                //!< Radius of root cell
      int nspawn;                                               //!< Threshold of NNODE for spawning new threads
      int level;                                                //!< Current tree level
      int iparent;                                              //!< Index of parent cell
      Nodes2cells(global_ptr<OctreeNode> _octNode, GB_iter _B0, GC_iter _C, // Constructor
		  GC_iter _C0, GC_iter _CN, vec3 _X0, real_t _R0,
		  int _nspawn, int _level=0, int _iparent=0) :
	octNode(_octNode), B0(_B0), C(_C), C0(_C0), CN(_CN),    // Initialize variables
	X0(_X0), R0(_R0), nspawn(_nspawn), level(_level), iparent(_iparent) {}
      //! Get cell index
      uint64_t getKey(vec3 X, vec3 Xmin, real_t diameter) const {
	int iX[3] = {0, 0, 0};                                  // Initialize 3-D index
	for (int d=0; d<3; d++) iX[d] = int((X[d] - Xmin[d]) / diameter);// 3-D index
	uint64_t index = ((1 << 3 * level) - 1) / 7;            // Levelwise offset
	for (int l=0; l<level; l++) {                           // Loop over levels
	  for (int d=0; d<3; d++) index += (iX[d] & 1) << (3 * l + d); // Interleave bits into Morton key
	  for (int d=0; d<3; d++) iX[d] >>= 1;                  //  Bitshift 3-D index
	}                                                       // End loop over levels
	return index;                                           // Return Morton key
      }
      int operator() () const {                                // Overload operator()
        int numLevels;
        int nchild = 0;
        global_ptr<OctreeNode> children[8];
        GC_iter CNs[8];
        GC_iter Ci = CN;                                       //   CN points to the next free memory address

        auto [cso, csc] =
          ityr::make_checkouts(octNode, 1, ityr::checkout_mode::read,
                               C      , 1, ityr::checkout_mode::write);
        auto o = &cso[0];

        Cell* c = new (&csc[0]) Cell();
        c->IPARENT = iparent;                                   //  Index of parent cell
        c->R       = R0 / (1 << level);                         //  Cell radius
        c->X       = o->X;                                //  Cell center
        c->NBODY   = o->NBODY;                            //  Number of decendant bodies
        c->IBODY   = o->IBODY;                            //  Index of first body in cell
        c->BODY    = B0 + c->IBODY;                             //  Iterator of first body in cell
        c->ICELL   = getKey(c->X, X0-R0, 2*c->R);               //  Get Morton key
        if (o->NNODE == 1) {                              //  If node has no children
          c->ICHILD = 0;                                        //   Set index of first child cell to zero
          c->NCHILD = 0;                                        //   Number of child cells
          assert(c->NBODY > 0);                                 //   Check for empty leaf cells
          numLevels = level;
        } else {                                                //  Else if node has children
          int octants[8];                                       //   Map of child index to octants
          for (int i=0; i<8; i++) {                             //   Loop over octants
            if (o->CHILD[i]) {                            //    If child exists for that octant
              octants[nchild] = i;                              //     Map octant to child index
              nchild++;                                         //     Increment child cell counter
            }                                                   //    End if for child existance
          }                                                     //   End loop over octants
          c->ICHILD = Ci - C0;                                  //   Set Index of first child cell
          c->NCHILD = nchild;                                   //   Number of child cells
          assert(c->NCHILD > 0);                                //   Check for childless non-leaf cells
          CN += nchild;                                         //   Increment next free memory address

          for (int i=0; i<nchild; i++) {
            int octant = octants[i];                            //    Get octant from child index
            children[i] = o->CHILD[octant];
            CNs[i] = CN;
            CN += children[i]->*(&OctreeNode::NNODE) - 1;            //    Increment next free memory address
          }                                                     //   End loop over octants
        }

        if (nchild > 0) {
          cso.checkin();
          csc.checkin();

          int numLevels_ = ityr::transform_reduce(
              ityr::execution::par,
              ityr::count_iterator<int>(0),
              ityr::count_iterator<int>(nchild),
              ityr::reducer::max<int>{},
              [=, *this](int i) {
                Nodes2cells nodes2cells(children[i],     //    Instantiate recursive functor
                                        B0, Ci+i, C0, CNs[i], X0, R0, nspawn, level+1, C-C0);
                return nodes2cells();
              });                                                     //   End loop over children

          for (int i=0; i<nchild; i++) {                        //   Loop over children
            free_octree_node(children[i]);
          }                                                     //   End loop over children
          numLevels = std::max(numLevels_, level+1);             //   Update maximum level of tree
        }

        return numLevels;
      }                                                         // End overload operator()
    };

    //! Transform Xmin & Xmax to X (center) & R (radius)
    Box bounds2box(Bounds bounds) {
      vec3 Xmin = bounds.Xmin;                                  // Set local Xmin
      vec3 Xmax = bounds.Xmax;                                  // Set local Xmax
      Box box;                                                  // Bounding box
      for (int d=0; d<3; d++) box.X[d] = (Xmax[d] + Xmin[d]) / 2; // Calculate center of domain
      box.R = 0;                                                // Initialize localRadius
      for (int d=0; d<3; d++) {                                 // Loop over dimensions
	box.R = std::max(box.X[d] - Xmin[d], box.R);            //  Calculate min distance from center
	box.R = std::max(Xmax[d] - box.X[d], box.R);            //  Calculate max distance from center
      }                                                         // End loop over dimensions
      box.R *= 1.00001;                                         // Add some leeway to radius
      return box;                                               // Return box.X and box.R
    }

    //! Grow tree structure top down
    void growTree(GBodies bodies, GBodies buffer, Box box, global_vec<BinaryTreeNode>& bin_nodes_vec) {
      assert(box.R > 0);                                        // Check for bounds validity
      if (ityr::is_master()) {
        logger::startTimer("Grow tree");                          // Start timer
      }
      B0 = bodies.begin();                                      // Bodies iterator
      int maxBinNode = (4 * bodies.size()) / nspawn;            // Get maximum size of binary tree

      bin_nodes_vec.resize(maxBinNode);

      BinaryTreeNode root_node;
      root_node.BEGIN = bin_nodes_vec.begin();
      root_node.END = bin_nodes_vec.end();

      N0 = ityr::root_exec([=, *this] {
        global_vec<BinaryTreeNode> root_bin_node_vec(1, root_node);

        BuildNodes buildNodes(bodies, buffer, 0, bodies.size(),
                              root_bin_node_vec.begin(), box.X, box.R, ncrit, nspawn);// Instantiate recursive functor
        return buildNodes();                                             // Recursively build octree nodes
      });

#if 0
      logger::printTitle("Grow tree");
      std::cout << std::setw(logger::stringLength) << std::left
		<< "Make node" << " : " << timer["Make node"] << " s\n"
		<< std::setw(logger::stringLength) << std::left
		<< "Count bodies" << " : " << timer["Count bodies"] << " s\n"
		<< std::setw(logger::stringLength) << std::left
		<< "Exclusive scan" << " : " << timer["Exclusive scan"] << " s\n"
		<< std::setw(logger::stringLength) << std::left
		<< "Move bodies" << " : " << timer["Move bodies"] << " s\n"
		<< std::setw(logger::stringLength) << std::left
		<< "Get node range" << " : " << timer["Get node range"] << " s\n"
		<< std::setw(logger::stringLength) << std::left
		<< "Total grow tree" << " : " << timer["Make node"] +
	timer["Count bodies"] + timer["Exclusive scan"] +
	timer["Move bodies"] + timer["Get node range"] << " s" << std::endl;
#endif
      if (ityr::is_master()) {
        logger::stopTimer("Grow tree");                           // Stop timer
      }
    }

    //! Link tree structure
    void linkTree(Box box, global_vec<Cell>& cells_vec) {
      if (ityr::is_master()) {
        logger::startTimer("Link tree");                          // Start timer
      }

      if (N0 != nullptr) {                                         // If the node tree is not empty
        std::size_t ncells = N0->*(&OctreeNode::NNODE);

	cells_vec.resize(ncells);                                //  Allocate cells array
	GC_iter C0 = cells_vec.begin();                              //  Cell begin iterator

        numLevels = ityr::root_exec([=, *this]() {
          Nodes2cells nodes2cells(N0, B0, C0, C0, C0+1, box.X, box.R, nspawn);// Instantiate recursive functor
          auto ret = nodes2cells();                                          //  Convert nodes to cells recursively
          free_octree_node(N0);
          return ret;
        });
      }                                                         // End if for empty node tree
      if (ityr::is_master()) {
        logger::stopTimer("Link tree");                           // Stop timer
      }
    }

  public:
    BuildTree(int _ncrit, int _nspawn) : ncrit(_ncrit), nspawn(_nspawn), numLevels(0) {}

    //! Build tree structure top down
    void buildTree(GBodies bodies, GBodies buffer, Bounds bounds,
                   global_vec<Cell>& cells_vec, global_vec<BinaryTreeNode>& bin_node_vec) {
      Box box = bounds2box(bounds);                             // Get box from bounds
      if (bodies.empty()) {                                     // If bodies vector is empty
	N0 = nullptr;                                              //  Reinitialize N0 with NULL
      } else {                                                  // If bodies vector is not empty
#if 0
        if (bodies.size() > buffer.size()) buffer.resize(bodies.size());// Enlarge buffer if necessary
#else
        assert(bodies.size() <= buffer.size());
#endif
	growTree(bodies, buffer, box, bin_node_vec);                          //  Grow tree from root
      }                                                         // End if for empty root
      linkTree(box, cells_vec);                                     // Form parent-child links in tree
    }

    //! Print tree structure statistics
    void printTreeData(GCells cells) {
      if (logger::verbose && !cells.empty()) {                  // If verbose flag is true
	logger::printTitle("Tree stats");                       //  Print title
        int nbody = cells.begin()->*(static_cast<int Cell::*>(&CellBase::NBODY));
	std::cout  << std::setw(logger::stringLength) << std::left//  Set format
		   << "Bodies"     << " : " << nbody << std::endl// Print number of bodies
		   << std::setw(logger::stringLength) << std::left//  Set format
		   << "Cells"      << " : " << cells.size() << std::endl// Print number of cells
		   << std::setw(logger::stringLength) << std::left//  Set format
		   << "Tree depth" << " : " << numLevels << std::endl;//  Print number of levels
      }                                                         // End if for verbose flag
    }
  };
}
#endif
