/**
 * @file dtree.hpp
 * @author Parikshit Ram (pram@cc.gatech.edu)
 *
 * Density Estimation Tree class
 *
 * This file is part of MLPACK 1.0.8.
 *
 * MLPACK is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * MLPACK is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 * details (LICENSE.txt).
 *
 * You should have received a copy of the GNU General Public License along with
 * MLPACK.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __MLPACK_METHODS_DET_DTREE_HPP
#define __MLPACK_METHODS_DET_DTREE_HPP

#include <mlpack/core.hpp>

namespace mlpack {
namespace det /** Density Estimation Trees */ {

/**
 * A density estimation tree is similar to both a decision tree and a space
 * partitioning tree (like a kd-tree).  Each leaf represents a constant-density
 * hyper-rectangle.  The tree is constructed in such a way as to minimize the
 * integrated square error between the probability distribution of the tree and
 * the observed probability distribution of the data.  Because the tree is
 * similar to a decision tree, the density estimation tree can provide very fast
 * density estimates for a given point.
 *
 * For more information, see the following paper:
 *
 * @code
 * @incollection{ram2011,
 *   author = {Ram, Parikshit and Gray, Alexander G.},
 *   title = {Density estimation trees},
 *   booktitle = {{Proceedings of the 17th ACM SIGKDD International Conference
 *       on Knowledge Discovery and Data Mining}},
 *   series = {KDD '11},
 *   year = {2011},
 *   pages = {627--635}
 * }
 * @endcode
 */
class DTree
{
 public:
  /**
   * Create an empty density estimation tree.
   */
  DTree();

  /**
   * Create a density estimation tree with the given bounds and the given number
   * of total points.  Children will not be created.
   *
   * @param maxVals Maximum values of the bounding box.
   * @param minVals Minimum values of the bounding box.
   * @param totalPoints Total number of points in the dataset.
   */
  DTree(const arma::vec& maxVals,
        const arma::vec& minVals,
        const size_t totalPoints);

  /**
   * Create a density estimation tree on the given data.  Children will be
   * created following the procedure outlined in the paper.  The data will be
   * modified; it will be reordered similar to the way BinarySpaceTree modifies
   * datasets.
   *
   * @param data Dataset to build tree on.
   */
  DTree(arma::mat& data);

  /**
   * Create a child node of a density estimation tree given the bounding box
   * specified by maxVals and minVals, using the size given in start and end and
   * the specified error.  Children of this node will not be created
   * recursively.
   *
   * @param maxVals Upper bound of bounding box.
   * @param minVals Lower bound of bounding box.
   * @param start Start of points represented by this node in the data matrix.
   * @param end End of points represented by this node in the data matrix.
   * @param error log-negative error of this node.
   */
  DTree(const arma::vec& maxVals,
        const arma::vec& minVals,
        const size_t start,
        const size_t end,
        const double logNegError);

  /**
   * Create a child node of a density estimation tree given the bounding box
   * specified by maxVals and minVals, using the size given in start and end,
   * and calculating the error with the total number of points given.  Children
   * of this node will not be created recursively.
   *
   * @param maxVals Upper bound of bounding box.
   * @param minVals Lower bound of bounding box.
   * @param start Start of points represented by this node in the data matrix.
   * @param end End of points represented by this node in the data matrix.
   */
  DTree(const arma::vec& maxVals,
        const arma::vec& minVals,
        const size_t totalPoints,
        const size_t start,
        const size_t end);

  //! Clean up memory allocated by the tree.
  ~DTree();

  /**
   * Greedily expand the tree.  The points in the dataset will be reordered
   * during tree growth.
   *
   * @param data Dataset to build tree on.
   * @param oldFromNew Mappings from old points to new points.
   * @param useVolReg If true, volume regularization is used.
   * @param maxLeafSize Maximum size of a leaf.
   * @param minLeafSize Minimum size of a leaf.
   */
  double Grow(arma::mat& data,
              arma::Col<size_t>& oldFromNew,
              const bool useVolReg = false,
              const size_t maxLeafSize = 10,
              const size_t minLeafSize = 5);

  /**
   * Perform alpha pruning on a tree.  Returns the new value of alpha.
   *
   * @param oldAlpha Old value of alpha.
   * @param points Total number of points in dataset.
   * @param useVolReg If true, volume regularization is used.
   * @return New value of alpha.
   */
  double PruneAndUpdate(const double oldAlpha,
                        const size_t points,
                        const bool useVolReg = false);

  /**
   * Compute the logarithm of the density estimate of a given query point.
   *
   * @param query Point to estimate density of.
   */
  double ComputeValue(const arma::vec& query) const;

  /**
   * Print the tree in a depth-first manner (this function is called
   * recursively).
   *
   * @param fp File to write the tree to.
   * @param level Level of the tree (should start at 0).
   */
  void WriteTree(FILE *fp, const size_t level = 0) const;

  /**
   * Index the buckets for possible usage later; this results in every leaf in
   * the tree having a specific tag (accessible with BucketTag()).  This
   * function calls itself recursively.
   *
   * @param tag Tag for the next leaf; leave at 0 for the initial call.
   */
  int TagTree(const int tag = 0);

  /**
   * Return the tag of the leaf containing the query.  This is useful for
   * generating class memberships.
   *
   * @param query Query to search for.
   */
  int FindBucket(const arma::vec& query) const;

  /**
   * Compute the variable importance of each dimension in the learned tree.
   *
   * @param importances Vector to store the calculated importances in.
   */
  void ComputeVariableImportance(arma::vec& importances) const;

  /**
   * Compute the log-negative-error for this point, given the total number of
   * points in the dataset.
   *
   * @param totalPoints Total number of points in the dataset.
   */
  double LogNegativeError(const size_t totalPoints) const;

  /**
   * Return whether a query point is within the range of this node.
   */
  bool WithinRange(const arma::vec& query) const;

 private:
  // The indices in the complete set of points
  // (after all forms of swapping in the original data
  // matrix to align all the points in a node
  // consecutively in the matrix. The 'old_from_new' array
  // maps the points back to their original indices.

  //! The index of the first point in the dataset contained in this node (and
  //! its children).
  size_t start;
  //! The index of the last point in the dataset contained in this node (and its
  //! children).
  size_t end;

  //! Upper half of bounding box for this node.
  arma::vec maxVals;
  //! Lower half of bounding box for this node.
  arma::vec minVals;

  //! The splitting dimension for this node.
  size_t splitDim;

  //! The split value on the splitting dimension for this node.
  double splitValue;

  //! log-negative-L2-error of the node.
  double logNegError;

  //! Sum of the error of the leaves of the subtree.
  double subtreeLeavesLogNegError;

  //! Number of leaves of the subtree.
  size_t subtreeLeaves;

  //! If true, this node is the root of the tree.
  bool root;

  //! Ratio of the number of points in the node to the total number of points.
  double ratio;

  //! The logarithm of the volume of the node.
  double logVolume;

  //! The tag for the leaf, used for hashing points.
  int bucketTag;

  //! Upper part of alpha sum; used for pruning.
  double alphaUpper;

  //! The left child.
  DTree* left;
  //! The right child.
  DTree* right;

 public:
  //! Return the starting index of points contained in this node.
  size_t Start() const { return start; }
  //! Return the first index of a point not contained in this node.
  size_t End() const { return end; }
  //! Return the split dimension of this node.
  size_t SplitDim() const { return splitDim; }
  //! Return the split value of this node.
  double SplitValue() const { return splitValue; }
  //! Return the log negative error of this node.
  double LogNegError() const { return logNegError; }
  //! Return the log negative error of all descendants of this node.
  double SubtreeLeavesLogNegError() const { return subtreeLeavesLogNegError; }
  //! Return the number of leaves which are descendants of this node.
  size_t SubtreeLeaves() const { return subtreeLeaves; }
  //! Return the ratio of points in this node to the points in the whole
  //! dataset.
  double Ratio() const { return ratio; }
  //! Return the inverse of the volume of this node.
  double LogVolume() const { return logVolume; }
  //! Return the left child.
  DTree* Left() const { return left; }
  //! Return the right child.
  DTree* Right() const { return right; }
  //! Return whether or not this is the root of the tree.
  bool Root() const { return root; }
  //! Return the upper part of the alpha sum.
  double AlphaUpper() const { return alphaUpper; }

  //! Return the maximum values.
  const arma::vec& MaxVals() const { return maxVals; }
  //! Modify the maximum values.
  arma::vec& MaxVals() { return maxVals; }

  //! Return the minimum values.
  const arma::vec& MinVals() const { return minVals; }
  //! Modify the minimum values.
  arma::vec& MinVals() { return minVals; }

 private:

  // Utility methods.

  /**
   * Find the dimension to split on.
   */
  bool FindSplit(const arma::mat& data,
                 size_t& splitDim,
                 double& splitValue,
                 double& leftError,
                 double& rightError,
                 const size_t maxLeafSize = 10,
                 const size_t minLeafSize = 5) const;

  /**
   * Split the data, returning the number of points left of the split.
   */
  size_t SplitData(arma::mat& data,
                   const size_t splitDim,
                   const double splitValue,
                   arma::Col<size_t>& oldFromNew) const;

};

}; // namespace det
}; // namespace mlpack

#endif // __MLPACK_METHODS_DET_DTREE_HPP
