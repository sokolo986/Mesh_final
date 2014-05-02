#ifndef CS207_GRAPH_HPP
#define CS207_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include "CS207/Util.hpp"
#include "Point.hpp"

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <map>


using namespace std;
/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
struct internal_Node;
 struct internal_Edge;
 private:

 public:
  // PUBLIC TYPE DEFINITIONS

  typedef V node_value_type;
  typedef E edge_value_type;
  /** Type of this graph. */
  typedef Graph graph_type;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  typedef Node node_type;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  typedef Edge edge_type;

  /** Type of indexes and sizes. Return type of Node::index() and
      Graph::num_nodes(), argument type of Graph::node. */
  typedef unsigned size_type;

  /** Type of node iterators, which iterate over all graph nodes. */
  class node_iterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class edge_iterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class incident_iterator;

  // CONSTRUCTOR AND DESTRUCTOR

  /** Construct an empty graph. */
  Graph():NodePoint_(), edges_()  {
  }
  /** Default destructor */
  ~Graph() = default;

  // NODES

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node :private totally_ordered<Node> {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Node x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node() {
    }

    /** Return this node's position. */
    Point position() const {

      return graph_->NodePoint_[uid_].point;//graph_->NodePoint_[uid_].point;
    }

    /* set position when given a point */
    void set_position(const Point & p)
    {
        graph_->NodePoint_[uid_].point = p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->NodePoint_[uid_].index_;//uid_;//graph_->NodePoint_[uid_].index_;//uid_;
    }

    // HW1 #1: YOUR CODE HERE
    // Supply definitions and specifications for:
    node_value_type& value()
    {
        return graph_->NodePoint_[uid_].value_;
    }


    const node_value_type& value() const
    {
        return graph_->NodePoint_[uid_].value_;
    }


    bool operator==(const Node& n) const{
        return graph_==n.graph_ && this->uid_== n.uid_;
    }
    bool operator< (const Node& n) const{
        if (graph_!=n.graph_)
            return true;
        return (this->uid_)< n.uid_;

    }



    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    size_type degree() const
    {
        if (graph_->adjmap_.count(*this)>0)
            return graph_->adjmap_[*this].size();
        else
            return 0;
    }

    incident_iterator edge_begin() const
    {
        return incident_iterator(this->graph_, *this, 0);
    }
    incident_iterator edge_end() const
    {
            return incident_iterator(this->graph_, *this, graph_->adjmap_[*this].size());
    }


   private:
    // Only Graph can access our private members
    friend class Graph;

    Graph* graph_;
    size_type uid_;
    Node(const Graph* graph,size_type uid):graph_(const_cast<Graph*>(graph)), uid_(uid){};

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return i2u_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new size() == old size() + 1
   * @post result_node.index() == old size()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position,const node_value_type & val = node_value_type ()) {
    internal_Node newNode;
    newNode.point = position;
    newNode.index_ =  i2u_.size() ;
    newNode.value_  = val;

    NodePoint_.push_back(newNode);
    i2u_.push_back(NodePoint_.size()-1);

    return Node(this, NodePoint_.size()-1);
  }



  /** Return the node with index @a i.
   * @pre 0 <= @a i < size()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < size());
    return Node(this,i2u_[i]);

  }

  /** Remove a node from the graph.
   * @param[in] n Node to be removed
   * @pre @a n is a valid node of this graph.
   * @post new size() == old size() - 1
   *
   * Can invalidate outstanding iterators. @a n becomes invalid, as do any
   * other Node objects equal to @a n. All other Node objects remain valid.
   *
   * Complexity: Polynomial in size().
   */
  void remove_node(const Node& n) {
    // HW1 #1: YOUR CODE HERE
    for (unsigned i = 0; i < i2u_.size(); ++i)
    {
        if (has_edge(n, node(i)))
		{
			remove_edge(n, node(i));
		}
    }
    for (auto it = i2u_.begin()+n.index()+1; it < i2u_.end(); ++it )
    {
		--NodePoint_[(*it)].index_ ;
    }
    i2u_.erase (i2u_.begin()+n.index());
    // now it should be empty for the adj vector as remove_edge function helped us to remove all the edges
    this->adjmap_.erase(n);


  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    NodePoint_.clear();
    edges_.clear();
    i2u_.clear();
    i2e_.clear();
    adjmap_.clear();
  }


  // EDGES

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge :private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->edges_[edgeId_].node_a;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph_->edges_[edgeId_].node_b;
    }

	size_type index() const {
      return graph_->edges_[edgeId_].index_;
    }

    double length() const
    {
        // to be implemented, return 0 for now
        return norm(node1().position() - node2().position());
    }
    // HW1 #1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    bool operator==(const Edge& ed) const
    {
        return (this->node1() == ed.node1() && this->node2() == ed.node2());
    }
    bool operator<(const Edge& ed ) const
    {
        return (this->node1() < ed.node1() && this->node2() < ed.node2());
    }
    edge_value_type& value()
    {
        return graph_->edges_[edgeId_].value_;
    }


    const edge_value_type& value() const
    {
        return graph_->edges_[edgeId_].value_;
    }
   private:
    // Only Graph can access our private members
    friend class Graph;

    Graph* graph_;
    size_type edgeId_;
    Edge(const Graph* graph,size_type edgeId):graph_(const_cast<Graph*>(graph)), edgeId_(edgeId){};

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return i2e_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges() && i >= 0);
    return Edge(this, i2e_[i]);        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return true if, for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (auto it = i2e_.begin(); it < i2e_.end(); ++it )
    {
        if ((edges_[(*it)].node_a == a &&  edges_[(*it)].node_b == b) ||
            (edges_[(*it)].node_b == a &&  edges_[(*it)].node_a  == b))
            return true;
    }

    return false;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b, const edge_value_type & val = edge_value_type ()) {
    // HW0: YOUR CODE HERE
    // test if we have an edge first.
    for (auto it = i2e_.begin(); it < i2e_.end(); ++it )
    {
        if ((edges_[(*it)].node_a == a &&  edges_[(*it)].node_b == b) ||
            (edges_[(*it)].node_b == a &&  edges_[(*it)].node_a  == b))
            return Edge(this, it - i2e_.begin());;
    }

    internal_Edge newEdge;
    newEdge.node_a = a;
    newEdge.node_b = b ;
	newEdge.index_ = i2e_.size() ;
	newEdge.value_  = val;
    edges_.push_back(newEdge);
	i2e_.push_back(edges_.size()-1);

    // add edge information to the hash_map,adjmap_, here we add the node a as the key, and the edge id to the vector
    adjmap_[a].push_back(newEdge.index_);
    adjmap_[b].push_back(newEdge.index_);

    return Edge(this,edges_.size()-1);

  }

  /** Remove an edge, if any, returning the number of edges removed.
   * @param[in] a,b The nodes potentially defining an edge to be removed.
   * @return 1 if old has_edge(@a a, @a b), 0 otherwise
   * @pre @a a and @a b are valid nodes of this graph
   * @post !has_edge(@a a, @a b)
   * @post new num_edges() == old num_edges() - result
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Can invalidate all edge and incident iterators.
   * Invalidates any edges equal to Edge(@a a, @a b). Must not invalidate
   * other outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type remove_edge(const Node& a, const Node& b) {
    // HW1 #1: YOUR CODE HERE
    if (! has_edge(a, b))
		return 0;


    // delete the item in i2e_ vector and decrement all the relavent index_ by 1
    for (auto it = i2e_.begin(); it < i2e_.end(); ++it )
    {
        if ((edges_[(*it)].node_a.index() == a.index() &&  edges_[(*it)].node_b.index() == b.index()) ||
            (edges_[(*it)].node_b.index() == a.index() &&  edges_[(*it)].node_a.index() == b.index()))
		{

		    			// here we delete the index from adj vector and update the respective index
			for (auto a_i = adjmap_[a].begin(); a_i < adjmap_[a].end(); ++a_i )
			{
               if ( (*a_i) == edges_[(*it)].index_)
                {
                    for (auto a_it = a_i+1; a_it < adjmap_[a].end(); ++a_it ) //adjmap_[a].begin()+ (a_i- adjmap_[a].begin())
                    {
                        --(*a_it) ;
                    }
                    adjmap_[a].erase(a_i);
                    break;
                }
			}

			for (auto b_i = adjmap_[b].begin(); b_i < adjmap_[b].end(); ++b_i )
			{
                if ((*b_i) == edges_[(*it)].index_)
                {
                    for (auto b_it = b_i+1; b_it < adjmap_[b].end(); ++b_it )
                    {
                        --(*b_it) ;
                    }
                    adjmap_[b].erase(b_i);
                    break;
                }
			}



            // here we update the actual index in the edge object

			for (auto i = i2e_.begin()+edges_[(*it)].index_+1; i < i2e_.end(); ++i )
			{
				--edges_[(*i)].index_ ;
			}
			i2e_.erase (i2e_.begin()+edges_[(*it)].index_);


			return 1;

		}

    }


  }

  /** Remove an edge, if any, returning the number of edges removed.
   * @param[in] e The edge to remove
   * @pre @a e is a valid edge of this graph
   * @pre has_edge(@a e.node1(), @a e.node2())
   * @post !has_edge(@a e.node1(), @a e.node2())
   * @post new num_edges() == old num_edges() - 1
   *
   * This is a synonym for remove_edge(@a e.node1(), @a e.node2()), but its
   * implementation can assume that @a e is definitely an edge of the graph.
   * This might allow a faster implementation.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Can invalidate all edge and incident iterators.
   * Invalidates any edges equal to Edge(@a a, @a b). Must not invalidate
   * other outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type remove_edge(const Edge& e) {
    // HW1 #1: YOUR CODE HERE

			// here we delete the index from adj vector and update the respective index
			Node a = e.node1;
			Node b = e.node2;
			for (auto a_i = adjmap_[a].begin(); a_i < adjmap_[a].end(); ++a_i )
			{
                if ((*a_i) == e.index())
                {
                    for (auto a_it = adjmap_[a].begin()+a_i+1; a_it < adjmap_[a].end(); ++a_it )
                    {
                        --(*a_it) ;
                    }
                    adjmap_[a].erase(a_i);
                    break;

                }
			}

			for (auto b_i = adjmap_[b].begin(); b_i < adjmap_[b].end(); ++b_i )
			{
                if ((*b_i) == e.index())
                {
                    for (auto b_it = adjmap_[b].begin()+b_i+1; b_it < adjmap_[b].end(); ++b_it )
                    {
                        --(*b_it) ;
                    }
                    adjmap_[b].erase(b_i);
                    break;
                }
			}




	for (auto it = i2e_.begin()+e.index()+1; it < i2e_.end(); ++it )
	{
        --edges_[(*it)].index_ ;
	}
	i2e_.erase (i2e_.begin()+e.index());

	return 0;
  }

  // ITERATORS

  /** @class Graph::node_iterator
   * @brief Iterator class for nodes. A forward iterator. */
  class node_iterator :private totally_ordered<node_iterator>{
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Node value_type;
    /** Type of pointers to elements. */
    typedef Node* pointer;
    /** Type of references to elements. */
    typedef Node& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid node_iterator. */
    node_iterator() {
      // HW1 #2: YOUR CODE HERE
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const
    {
        auto it = Node(graph_,graph_->i2u_[nIteratorId_]);
        return it;
    }

    node_iterator& operator++()
    {
        nIteratorId_ ++ ;
        return *this;
    }
    bool operator==(const node_iterator& nit) const
    {
        return (this->graph_ == nit.graph_ && this->nIteratorId_ == nit.nIteratorId_);
    }

   private:
    friend class Graph;
    friend class edge_iterator;
    Graph* graph_;
    size_type nIteratorId_;
    node_iterator(const Graph* graph,size_type nIteratorId):graph_(const_cast<Graph*>(graph)), nIteratorId_(nIteratorId){};

    // HW1 #3: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const
  {
      return node_iterator(this, 0);
  }
  node_iterator node_end() const
  {
      return node_iterator(this, this->i2u_.size());
  }


  /** @class Graph::edge_iterator
   * @brief Iterator class for edges. A forward iterator. */
  class edge_iterator  :private totally_ordered<edge_iterator>{
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid edge_iterator. */
    edge_iterator() {
      // HW1 #3: YOUR CODE HERE
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const
    {
        return Edge(graph_,graph_->i2e_[eIteratorId_]);
    }
    edge_iterator& operator++()
    {
        eIteratorId_ ++ ;
        return *this;
    }
    bool operator==(const edge_iterator& eit) const
    {
        return (this->graph_ == eit.graph_ && this->eIteratorId_ == eit.eIteratorId_);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE

    Graph* graph_;
    size_type eIteratorId_;
    edge_iterator(const Graph* graph,size_type eIteratorId):graph_(const_cast<Graph*>(graph)), eIteratorId_(eIteratorId){};



  };

  // HW1 #3: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const
  {
    return edge_iterator(this, 0);
  }
   edge_iterator edge_end() const
   {
    return edge_iterator(this, this->i2e_.size());

   }


  /** @class Graph::incident_iterator
   * @brief Iterator class for edges incident to a given node. A forward
   * iterator. */
  class incident_iterator :private totally_ordered<incident_iterator>{
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid incident_iterator. */
    incident_iterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const
    {
       size_type edgeIndex  =  this->graph_->adjmap_[thisnode_][adjedgeId_];

        return graph_->edge(edgeIndex);
    }
    incident_iterator& operator++()
    {
        ++adjedgeId_;
        return *this;
    }
    bool operator==(const incident_iterator& it) const
    {
        return (this->graph_ == it.graph_ && this->thisnode_ == it.thisnode_ && this->adjedgeId_ == it.adjedgeId_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE

    Graph* graph_;
    size_type adjedgeId_;
    Node thisnode_;
    incident_iterator(const Graph* graph,Node thisnode, size_type adjedgeId):graph_(const_cast<Graph*>(graph)),thisnode_(thisnode), adjedgeId_(adjedgeId){};

  };

 private:
     struct internal_Node {
      Point point;   // The text held by an element
      node_value_type value_;
      size_type index_;      // The unique identifcation for an element
    };

    vector<internal_Node> NodePoint_;
    vector<size_type> i2u_;

     struct internal_Edge {
      node_type node_a;
      node_type node_b;
      edge_value_type value_;
	  size_type index_;
    };
    vector<internal_Edge> edges_;
	vector<size_type> i2e_;

    map<Node, vector<size_type> > adjmap_; // store the node object as the key, and a vector of Edge index as the value

};

#endif
