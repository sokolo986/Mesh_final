/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Graph.hpp"
#include "Point.hpp"
#include <map>
#include <algorithm>

// Gravity in meters/sec^2
static constexpr double grav = 9.81;

//Define Graph Type
//Create a type nodeValueType_vm, containing the velocity of the node, initial length of the node and the mass of the node
struct nodeValueType_vm
{
    Point vel;
    double mass;
};

struct EdgeValueType_Kl
{
/** Edge vale type, it contains value K, spring constant and length, which is distance between node
**/
    double k;
    double length;
};

typedef Graph<nodeValueType_vm, EdgeValueType_Kl> GraphType;
typedef GraphType::Node Node;
typedef GraphType::Edge Edge;

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g Graph
 * @param[in] t The current time (useful for time-dependent forces)
 * @param[in] dt The time step
 * @param[in] force Function object defining the force per node
 * @pre G::node_value_type supports ???????? YOU CHOOSE
 * @return the next time step (usually @a t + @a dt)
 *
 * @a force is called as @a force(n, @a t), where n is a node of the graph
 * and @a t is the current time parameter. @a force must return a Point
 * representing the node's force at time @a t.
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {

  auto NodeNumber = g.num_nodes();
  double mass = 1.0/NodeNumber;

 // Calculate new positions and set the positions
  for (auto it = g.node_begin(); it != g.node_end(); ++it)
  {
    Point orig = (*it).position();
    auto vel = (*it).value().vel;
    Point newLocation = orig + vel * dt;
    (*it).set_position(newLocation);
   }

  // apply constraint
  constraint(g, t);

  // Calculate velocity and update velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it)
  {
    auto vel = (*it).value().vel;
    auto f = force((*it),t);
    (*it).value().vel = vel + f* (dt /mass);
  }
  return t + dt;
}



struct GravityForce {
  /** Struct/Class of gravity force
  * @param[in] node object with value function.
  * @param[out] Point object representing Gravity force
  * @post return a Point object representing Gravity force, Point should be (0, 0, - mass of the node * gravity constant)
  *
  */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point Gforce = Point(0,0,-n.value().mass*grav);
    return Gforce;
    }
};


struct MassSpringForce {
  /** Struct/Class of Mass Spring force
  * @param[in] node object with value function.
  * @param[out] Point object representing mass spring force
  * @post return a Point object representing mass spring force based on the calculation : difference =  (current position - adjacent nodes position)
  * @post mass spring force = \sum -K * difference / Norm (difference) * ( norm(difference) - length), length is provided as initial length between nodes before movement
  *
  */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point Sforce = Point(0,0,0);
    for (auto it = n.edge_begin(); it != n.edge_end();++it)
    {
        auto adj = n;
        if (n == (*it).node1())
            adj = (*it).node2();
        else
            adj = (*it).node1();
        Point diff = n.position() - adj.position();
        Sforce += (-(*it).value().k) * diff/norm(diff) * ( norm(diff) - (*it).value().length);
    }
    return Sforce;
    }
};

struct DampingForce {
  /** Struct/Class of Damping force
  * @param[in] node object with value function.
  * @param[out] Point object representing Damping force
  * @post return a Point object representing Damping force based on the calculation : -1/N * velocity, where N is the total number of nodes. Therefore here it could be defined as -mass * velocity
  */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return -n.value().vel * n.value().mass;
    }
};


struct constraint_fixPoint{
  /** Struct/Class of constraint on fixed points. If node's position is (0,0,0) or (1,0,0), their position dont change.
  * @param[in] graph object G and time t
  * @param[out] No output
  * @post if the node is either (0,0,0) or (1,0,0), set the node velocity to 0.
  */
  Node node0_;
  Node node1_;
  constraint_fixPoint(Node node0, Node node1): node0_(node0), node1_(node1){}
  template <typename G>
  void operator()(G& g, double t){
    node0_.value().vel = Point(0,0,0);
    node0_.set_position(Point(0,0,0));
    node1_.value().vel = Point(0,0,0);
    node1_.set_position(Point(1,0,0));

  }

};



struct constraint_z{
  /** Struct/Class of constraint z = -0.75
  * @param[in] graph object G and time t
  * @param[out] No output
  * @post if the node violate innerproduct (x, (0,0,1)) < -0.75, set the node position to nearest point on the plane, and velocity on z direction to 0
  */
  template <typename G>
  void operator()(G& g, double t){
    Point constPoint(0,0,1);
    double boundary = -0.75;
    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      if (inner_prod((*it).position(), constPoint) < boundary){
        auto position = (*it).position();
        position.z = -0.75;
        (*it).set_position(position);
        (*it).value().vel.z = 0;
      }

    }
  }
};

struct constraint_sphere{
  /** Struct/Class of constraint on a sphere. Sphere center = (0.5, 0.5, -0.5), radius = 0.15
  * @param[in] graph object G and time t
  * @param[out] No output, updated location and velocity of graph object
  * @post if node violates constraint |x - c| < r, where x is location of node, c is center, r is radius
  * @post We will fix the node by set the position to nearest point on surface
  * @post we will also set the velocity to v_i = v_i - innerproduct(v_i, R_i)*R_i, R_i is (x-c)/|x-c|
  * We will use the algorithm first calculate l = distance to center, and use r/l as the ratio to apply on each direction
  * for example x_new = (r/l *(x-x_0) + x_0, x_0 is the center x coordinate
  */
  template <typename G>
  void operator()(G& g, double t, Point center=Point(0.5,0.5,-0.5), double radius = 0.15){

    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      double l = norm((*it).position()-center);
      if (l<radius && l !=0 ){
        double ratio = radius/l;
        auto position = (*it).position();
        auto R = (position - center) / norm(position - center);
        position.x = ratio * (position.x - center.x)+center.x;
        position.y = ratio * (position.y - center.y)+center.y;
        position.z = ratio * (position.z - center.z)+center.z;
        (*it).set_position(position);
        // update velocity now.
        (*it).value().vel -= inner_prod((*it).value().vel, R) * R;
      }

    }
  }
};


struct MyComparator {
  /** Struct/Class of comparator to compare edge length
  * @param[in] two edge objects
  * @param[out] boolean, true if the first edge length is smaller than 2nd edge length.
  * @post return true if the first edge length is smaller than 2nd edge length.
  */
   template <typename EDGE>
   bool operator()(const EDGE& e1, const EDGE& e2) const {
    return e1.length() < e2.length();
  }
}MyComparator;

struct constraint_selfCollision{
  /** Struct/Class of constraint on self collision
  * @param[in] graph object G and time t
  * @param[out] No output, updated location and velocity of graph object
  * @post if node violates constraint |x - c| < r, where x is location of node, c is the other nodes' location, r is radius of the other node
  * @post We will fix the node by set the position to nearest point on surface of other node
  * @post we will also set the velocity to v_i = v_i - innerproduct(v_i, R_i)*R_i, R_i is (x-c)/|x-c|
  * We will use the algorithm first calculate l = distance to center, and use r/l as the ratio to apply on each direction
  * for example x_new = (r/l *(x-x_0) + x_0, x_0 is the center x coordinate
  * Initial value of the radius is set to the minimum of the edge length of the node
  * Note: I have implemented the radius as 0.9 * minimum edge length of the node, so it makes more sense, otherwise the node couldnt move smoothly.
  */
  template <typename G>
  void operator()(G& g, double t){

    std::map<Node, double> NodeRadius;
    for (auto it = g.node_begin(); it !=g.node_end(); ++it){
      auto minEdge = *std::min_element((*it).edge_begin(), (*it).edge_end(), MyComparator);
      NodeRadius[*it] = 0.9*minEdge.length();
    }

    for (auto mi = NodeRadius.begin(); mi != NodeRadius.end(); ++ mi){
      auto cs = constraint_sphere();
      cs(g, t, (*mi).first.position(), (*mi).second);
    }
  }
};


template<class constr1, class constr2, class constr3, class constr4>
struct make_combined_constraint_combiner {
 /** Struct/Class of constraint combiner. Combines 4 different constraint to 1 function
  * @param[in] constraint 1, 2, 3,4. They need to take in parameter g and t, where g is graph, t is time
  * @param[out] No output, will call each constraint in sequence
  * @post All constraint functions will be called
  */

  constr1 _first;
  constr2 _second;
  constr3 _third;
  constr4 _fourth;
  make_combined_constraint_combiner(constr1& c1,  constr2& c2,  constr3& c3,  constr4& c4):_first(c1), _second(c2),_third(c3),_fourth(c4){}

  template<typename G>
  void operator()(G& g, double t)
  {
      _first(g, t);
      _second(g, t);
      _third(g, t);
      _fourth(g,t);
  }
};


 /** This is a wrapper function,
  * @param[in] constraint 1, 2, 3, 4. They need to take in parameter g and t, where g is graph, t is time
  * @param[out] No output, will call each constraint in sequence
  * @post All constraint functions will be called
  */
template<class _T1, class _T2, class _T3, class _T4>
make_combined_constraint_combiner<_T1, _T2, _T3, _T4> make_combined_constraint(_T1&& _first, _T2&& _second, _T3&& _third, _T4&& _fourth)
{
    return make_combined_constraint_combiner<_T1, _T2, _T3, _T4>(_first, _second, _third, _fourth);
}



template<class force1, class force2>
struct make_combined_force_combiner {
 /** Struct/Class of force combiner. Combines 2 different force to 1
  * @param[in] force 1, 2. They need to take in parameter g and t, where g is graph, t is time
  * @param[out] No output, will call each force in sequence
  * @post All force functions will be called
  */
  force1 _first;
  force2 _second;
  make_combined_force_combiner(force1& f1,  force2& f2):_first(f1), _second(f2){}

  template<typename Node>
  Point operator()(Node n, double t)
  {
      return _first(n,  t)+_second(n,  t);
  }
};


template<class force1, class force2, class force3>
struct make_combined_force_combiner2 {
 /** Struct/Class of force combiner. Combines 3 different force to 1
  * @param[in] force 1, 2, 3. They need to take in parameter g and t, where g is graph, t is time
  * @param[out] No output, will call each force in sequence
  * @post All force functions will be called
  */
  force1 _first;
  force2 _second;
  force3 _third;
  make_combined_force_combiner2(force1& f1,  force2& f2,  force3& f3):_first(f1), _second(f2),_third(f3){}

  template<typename Node>
  Point operator()(Node n, double t)
  {
      return _first(n,t)+_second(n,t) + _third(n,t);
  }
};


 /** This is a wrapper function, to wrap on make_combined_force_combiner
  * @param[in] force 1, 2. They need to take in parameter g and t, where g is graph, t is time
  * @param[out] No output, will call each force in sequence
  * @post All force functions will be called
  */
template<class _T1, class _T2>
make_combined_force_combiner<_T1, _T2> make_combined_force(_T1&& _first, _T2&& _second)
{
    return make_combined_force_combiner<_T1, _T2>(_first, _second);
}

 /** This is a wrapper function, to wrap on make_combined_force_combiner
  * @param[in] force 1, 2,3. They need to take in parameter g and t, where g is graph, t is time
  * @param[out] No output, will call each force in sequence
  * @post All force functions will be called
  */
template<class _T1, class _T2, class _T3>
make_combined_force_combiner2<_T1, _T2, _T3> make_combined_force(_T1&& _first, _T2&& _second, _T3&& _third)
{
    return make_combined_force_combiner2<_T1, _T2, _T3>(_first, _second, _third);
}


/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    if(n.position() ==Point(0,0,0) ||n.position() ==Point(1,0,0))
        return Point(0,0,0);

    //calculate the gravity force first
    Point Gforce = Point(0,0,-n.value().mass*grav);

    // Now calculate the spring force
    Point Sforce = Point(0,0,0);
    for (auto it = n.edge_begin(); it != n.edge_end();++it)
    {
        auto adj = n;
        if (n == (*it).node1())
            adj = (*it).node2();
        else
            adj = (*it).node1();

        Point diff = n.position() - adj.position();
        Sforce += (-(*it).value().k) * diff/norm(diff) * ( norm(diff) - (*it).value().length);
    }

    return Sforce + Gforce;
    }
};


int main(int argc, char* argv[]) {
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  GraphType graph;
  std::vector<typename GraphType::node_type> nodes;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));


std::cout << graph.num_nodes() << std::endl;
  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CS207::getline_parsed(tets_file, t)) {
    for (unsigned i = 1; i < t.size(); ++i) {
      graph.add_edge(nodes[t[0]], nodes[t[1]]);
      graph.add_edge(nodes[t[0]], nodes[t[2]]);
//#if 0
      // Diagonal edges: include as of HW2 #2
      graph.add_edge(nodes[t[0]], nodes[t[3]]);
      graph.add_edge(nodes[t[1]], nodes[t[2]]);
//#endif
      graph.add_edge(nodes[t[1]], nodes[t[3]]);
      graph.add_edge(nodes[t[2]], nodes[t[3]]);
    }
  }

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  viewer.launch();

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();

  // Begin the mass-spring simulation
  // Define some parameters, k is the spring constant
  double dt = 0.0001;
  double t_start = 0;
  double t_end = 5.0;
  double k = 100;
// update the speed to 0, mass to 1/Node Numbers and initial length to the length() function defined as Graph.hpp
  auto NodeNumber = graph.num_nodes();
  double mass = 1.0/NodeNumber;
// initialize two nodes, their position are fixed.
  Node node000, node001;
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it)
  {
      (*it).value().vel = Point(0,0,0);
      (*it).value().mass = mass;
      if ((*it).position() == Point(0,0,0))
        node000=(*it);
      if ((*it).position() == Point(1,0,0))
        node001 = (*it);
  }
  // Initialize fixpoint constraint where these two nodes are not changeing the position.
  constraint_fixPoint fixpoint(node000, node001);

  // Initilize the length to be length at time 0
  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it )
  {
      (*it).value().length = (*it).length();
      (*it).value().k = k;
  }

  for (double t = t_start; t < t_end; t += dt) {
    symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce()) //
                    ,make_combined_constraint(constraint_fixPoint(node000, node001), constraint_z(), constraint_sphere(), constraint_selfCollision()) );
    viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
    viewer.set_label(t);

    // These lines slow down the animation for small graphs, like grid0_*.
    if (graph.size() < 100)
      CS207::sleep(0.001);
  }

  return 0;
}
