/**
 * @file shallow_water.cpp
 * Implementation of a shallow water system using Mesh
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D point list (one per line) defined by three doubles
 * Second file: Triangles (one per line) defined by 3 indices into the point list
 */

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"
#include "Point.hpp"
#include <fstream>
#include "omp.h"
#include "Mesh.hpp"
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <chrono>

// Standard gravity (average gravity at Earth's surface) in meters/sec^2
static constexpr double grav = 9.80665;

/** Water column characteristics */
typedef struct QVar {
  double h;	  // Height of fluid
  double hu;	// Height times average x velocity of column
  double hv;	// Height times average y velocity of column

  /** Default constructor.
   *
   * A default water column is 1 unit high with no velocity. */
  QVar()
    : h(1), hu(0), hv(0) {
  }
  /** Construct the given water column. */
  QVar(double h_, double hu_, double hv_)
    : h(h_), hu(hu_), hv(hv_) {
  }


  QVar operator+(const QVar& q){
	double hnew = h+q.h;
	double hunew = hu+q.hu;
	double hvnew = hv+q.hv;
	//return QVar(h+q.h, hu + q.hu, hv + q.hv);
	return QVar(hnew, hunew, hvnew);
  }

  QVar operator-(const QVar& q){
	return QVar(h-q.h, hu - q.hu, hv - q.hv);
  }



  bool operator==(const QVar& q) const {
    if (h==q.h && hu==q.hu && hv == q.hv)
		return false;
	else
		return true;
  }
  bool operator!=(const QVar& q) const {
    return !(*this == q);
  }

  QVar& operator=(const QVar& q){
	h = q.h;
	hu = q.hu;
	hv = q.hv;
	return *this;
  }

  QVar operator* (double n){
	return QVar(h*n, hu*n, hv*n);
  }

  QVar operator/ (double n){
	return QVar(h/n, hu/n, hv/n);

  }

  QVar& operator+=(const QVar& b) {
    h+=b.h;
	hu += b.hu;
	hv += b.hv;
    return *this;
  }

  QVar& operator-=(const QVar& b) {
    h -=b.h;
	hu -= b.hu;
	hv -= b.hv;
    return *this;
  }


}QVar;


QVar operator*(double n, QVar& q){
	return QVar(n*q.h,n*q.hu,n*q.hv);
}

QVar operator/(double n,QVar& q) {
	return QVar(n/q.h, n/q.hu, n/q.hv);
}


/** Function object for calculating shallow-water flux.
 *          |e
 *   T_k    |---> n = (nx,ny)   T_m
 *   QBar_k |                   QBar_m
 *          |
 * @param[in] nx, ny Defines the 2D outward normal vector n = (@a nx, @a ny)
 *            from triangle T_k to triangle T_m. The length of n is equal to the
 *            the length of the edge, |n| = |e|.
 * @param[in] dt The time step taken by the simulation. Used to compute the
 *               Lax-Wendroff dissipation term.
 * @param[in] qk The values of the conserved variables on the left of the edge.
 * @param[in] qm The values of the conserved variables on the right of the edge.
 * @return The flux of the conserved values across the edge e
 */
struct EdgeFluxCalculator {
  QVar operator()(double nx, double ny, double dt,
                  const QVar& qk, const QVar& qm) {
    double e_length = sqrt(nx*nx + ny*ny);
    nx /= e_length;
    ny /= e_length;

    // The velocities normal to the edge
    double wm = (qm.hu*nx + qm.hv*ny) / qm.h;
    double wk = (qk.hu*nx + qk.hv*ny) / qk.h;

    // Lax-Wendroff local dissipation coefficient
    double vm = sqrt(grav*qm.h) + sqrt(qm.hu*qm.hu + qm.hv*qm.hv) / qm.h;
    double vk = sqrt(grav*qk.h) + sqrt(qk.hu*qk.hu + qk.hv*qk.hv) / qk.h;
    double a  = dt * std::max(vm*vm, vk*vk);

    // Helper values
    double scale = 0.5 * e_length;
    double gh2   = 0.5 * grav * (qm.h*qm.h + qk.h*qk.h);

    // Simple flux with dissipation for stability
    return QVar(scale * (wm*qm.h  + wk*qk.h)           - a * (qm.h  - qk.h),
                scale * (wm*qm.hu + wk*qk.hu + gh2*nx) - a * (qm.hu - qk.hu),
                scale * (wm*qm.hv + wk*qk.hv + gh2*ny) - a * (qm.hv - qk.hv));
  }
};

/** Node position function object for use in the SDLViewer. */
struct NodePosition {
  template <typename NODE>
  Point operator()(const NODE& n) {
	// return the height stored in node value as the z direction
    return Point(n.position().x, n.position().y, n.value().h);
  }
};

// Define NodeData, EdgeData, TriData, etc
typedef Mesh<QVar,double,QVar> MeshType;


/** Integrate a hyperbolic conservation law defined over the mesh m
 * with flux functor f by dt in time.
 */
template <typename MESH, typename FLUX>
/*
 * implementation of euler approximation of the shallow water PDE
 * @pre: a valid Mesh class instance @mesh
 * @post: all triangle in the mesh class have new value() = old value - dt/area() * total flux,
		  where total flux is calculated by all three edges of the triangle
   @return: return total time t+dt
*/
double hyperbolic_step(MESH& mesh, FLUX& f, double t, double dt) {
  // Step the finite volume model in time by dt.
  // Implement Equation 7 from your pseudocode here.




  /*
  provide interface for parallel computing on all the node
  */

  struct step{
    double dt,t;
    FLUX f;
    MESH mesh;
    step(double dt, double t, FLUX& f, MESH& mesh): dt(dt),t(t), f(f), mesh(mesh) {}
    void operator()(MeshType::tri_iterator it){
	// value function will return the flux
	QVar total_flux=QVar(0,0,0);
	QVar qm = QVar(0,0,0);
	// iterate through 3 edges of a triangle
	auto edgetemp = (*it).edge1();
	for (int num = 0; num < 3; num++)
	{
		if (num ==0)
			edgetemp= (*it).edge1();
		else if (num==1)
			edgetemp = (*it).edge2();
		else
			edgetemp = (*it).edge3();

		if (  mesh.has_neighbor(edgetemp.index()) ) // it has a common triangle
		{
			auto nx =  ((*it).norm_vector(edgetemp)).x;
			auto ny =  ((*it).norm_vector(edgetemp)).y;

			// find the neighbour of a common edge
			for (auto i = mesh.tri_edge_begin(edgetemp.index()); i != mesh.tri_edge_end(edgetemp.index()); ++i){
				if (!(*i==*it))
					qm = (*i).value();
			}
			// calculat the total flux
			total_flux += f(nx, ny, dt, (*it).value(), qm);
		}
		else{
			// when it doesnt have a neighbour shared with this edge
			auto nx =  ((*it).norm_vector(edgetemp)).x;
			auto ny =  ((*it).norm_vector(edgetemp)).y;
			qm = QVar((*it).value().h, 0, 0 ); // approximation

			total_flux += f(nx, ny, dt, (*it).value(), qm);
		}
	}

	(*it).value() +=  total_flux * (- dt / (*it).area());

    }
  };

  step step1(dt, t, f, mesh);
  applytoall(mesh.tri_begin(),mesh.tri_end(), step1, 4);







/*
omp_set_num_threads(8);
#pragma omp parallel
{
  for (auto it = mesh.tri_begin(); it!=mesh.tri_end() ; ++it)
  {
  #pragma omp single nowait
  {
	// value function will return the flux
	QVar total_flux=QVar(0,0,0);
	QVar qm = QVar(0,0,0);
	// iterate through 3 edges of a triangle
	auto edgetemp = (*it).edge1();
	for (int num = 0; num < 3; num++)
	{
		if (num ==0)
			edgetemp= (*it).edge1();
		else if (num==1)
			edgetemp = (*it).edge2();
		else
			edgetemp = (*it).edge3();

		if (  mesh.has_neighbor(edgetemp.index()) ) // it has a common triangle
		{
			auto nx =  ((*it).norm_vector(edgetemp)).x;
			auto ny =  ((*it).norm_vector(edgetemp)).y;

			// find the neighbour of a common edge
			for (auto i = mesh.tri_edge_begin(edgetemp.index()); i != mesh.tri_edge_end(edgetemp.index()); ++i){
				if (!(*i==*it))
					qm = (*i).value();
			}
			// calculat the total flux
			total_flux += f(nx, ny, dt, (*it).value(), qm);
		}
		else{
			// when it doesnt have a neighbour shared with this edge
			auto nx =  ((*it).norm_vector(edgetemp)).x;
			auto ny =  ((*it).norm_vector(edgetemp)).y;
			qm = QVar((*it).value().h, 0, 0 ); // approximation

			total_flux += f(nx, ny, dt, (*it).value(), qm);
		}
	}

	(*it).value() +=  total_flux * (- dt / (*it).area());
	}
  }
 }
*/
  return t + dt;
}

  struct post{

    MeshType m;
    post(MeshType& m): m(m) {}
    void operator()(MeshType::node_iterator it){
	QVar sum = QVar(0,0,0);
	double sumTriArea = 0;
	// for each node, iterate through its adjacent triangles

	for (auto adji = m.vertex_begin((*it).index()); adji !=  m.vertex_end((*it).index()); ++ adji)
	{
		auto tri = (*adji);
		sum += tri.area() * tri.value();
		sumTriArea += tri.area();
	}

	(*it).value() = sum/sumTriArea; // update nodes value
    }
  };

/** Convert the triangle-averaged values to node-averaged values for viewing. */
template <typename MESH>
/*
 * post process of a mesh instance to update the values of all the nodes values
 * @pre: a valid mesh instance
 * @post: update the nodes values based on approximation of the average of neighbours values
		  refer to equation 9 on the notes
*/
void post_process(MESH& m) {
  // Translate the triangle-averaged values to node-averaged values
  // Implement Equation 8 from your pseudocode here
  	// iterate through all the nodes



omp_set_num_threads(4);
#pragma omp parallel
{
  for ( auto it = m.node_begin(); it!= m.node_end(); ++it)
  {
  #pragma omp single nowait
  {
	QVar sum = QVar(0,0,0);
	double sumTriArea = 0;
	// for each node, iterate through its adjacent triangles

	for (auto adji = m.vertex_begin((*it).index()); adji !=  m.vertex_end((*it).index()); ++ adji)
	{
		auto tri = (*adji);
		sum += tri.area() * tri.value();
		sumTriArea += tri.area();
	}

	(*it).value() = sum/sumTriArea; // update nodes value
  }
  }

}





  //post post(m);
  //applytoall(m.node_begin(),m.node_end(), post, 4);




}


struct EdgeComparator {
  /** Struct/Class of comparator to compare edge length
  * @param[in] two triangle objects
  * @param[out] boolean, true if the first triangle has the smallest edge length than the second triangle
  */
   template <typename Edge>
   bool operator()(const Edge& t1, const Edge& t2) const {
	return t1.length() < t2.length();
  }
}EdgeComparator;

struct HeightComparator {
  /** Struct/Class of comparator to compare height value stored in Node
  * @param[in] two Node objects
  * @param[out] boolean, true if the height in first Node value  is smaller than  height in Second Node value
  */
   template <typename Node>
   bool operator()(const Node& t1, const Node& t2) const {
	return t1.value().h < t2.value().h;
  }
}HeightComparator;



int main(int argc, char* argv[])
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: shallow_water NODES_FILE TRIS_FILE\n";
    exit(1);
  }

  MeshType mesh;
  // HW4B: Need node_type before this can be used!
  omp_set_num_threads(4);

auto start = std::chrono::high_resolution_clock::now();


  std::vector<typename MeshType::node_type> mesh_node;

  // Read all Points and add them to the Mesh
  std::ifstream nodes_file(argv[1]);
  Point p;
  while (CS207::getline_parsed(nodes_file, p)) {
    mesh_node.push_back(mesh.add_node(p));
  }

  // Read all mesh triangles and add them to the Mesh
  std::ifstream tris_file(argv[2]);
  std::array<int,3> t;
  while (CS207::getline_parsed(tris_file, t)) {
    // HW4B: Need to implement add_triangle before this can be used!
    mesh.add_triangle(mesh_node[t[0]], mesh_node[t[1]], mesh_node[t[2]]);
  }

  // Print out the stats
  std::cout << mesh.num_nodes() << " "
            << mesh.num_edges() << " "
            << mesh.num_triangles() << std::endl;

  // HW4B Initialization
  // Set the initial conditions
  int wave = 0, peddle=0, dam=1;
  if (wave){
	for ( auto it = mesh.node_begin(); it!= mesh.node_end(); ++it){
		auto x = (*it).position().x;
		auto y = (*it).position().y;
		double h = 1-0.75 * exp(-80 * ( (x-0.75)*(x-0.75) + y*y ));
		mesh.value((*it),QVar( h,0,0));
	}
  }
  else if (peddle){
	for ( auto it = mesh.node_begin(); it!= mesh.node_end(); ++it){
		auto x = (*it).position().x;
		auto y = (*it).position().y;
		double h = (x-0.75)*(x-0.75) + y*y -0.15*0.15 ;
		int H =0;
		if (h < 0)
			H = 1;
		mesh.value((*it), QVar(1+0.75*H,0,0));
	}
  }
  else if (dam){

omp_set_num_threads(8);
#pragma omp parallel
{
	for ( auto it = mesh.node_begin(); it!= mesh.node_end(); ++it){
	#pragma omp single nowait
	{
		auto x = (*it).position().x;
		int H =0;
		if (x < 0)
			H = 1;
		mesh.value((*it), QVar(1+0.75*H,0,0));
	}
	}
}
  }


  // Perform any needed precomputation
// initialize triangle
//omp_set_num_threads(4);
#pragma omp parallel
{
  for (auto it = mesh.tri_begin(); it < mesh.tri_end(); ++it ) {
  	#pragma omp single nowait
  {
	(*it).value() = ((*it).node1().value() + (*it).node2().value() + (*it).node3().value())/3.0;
  }
  }
}

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch();

  // HW4B: Need to define Mesh::node_type and node/edge iterator
  // before these can be used!

  auto node_map = viewer.empty_node_map(mesh);
  viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
                   CS207::DefaultColor(), NodePosition(), node_map);
  viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);

  viewer.center_view();

  // HW4B: Timestep
  // CFL stability condition requires dt <= dx / max|velocity|
  // For the shallow water equations with u = v = 0 initial conditions
  //   we can compute the minimum edge length and maximum original water height
  //   to set the time-step
  // Compute the minimum edge length and maximum water height for computing dt
  auto min_length = *std::min_element(mesh.edge_begin(), mesh.edge_end(), EdgeComparator);

  auto max_h = *std::max_element(mesh.node_begin(), mesh.node_end(), HeightComparator);

  double dt = 0.25 * min_length.length() / (sqrt(grav * max_h.value().h));
  double t_start = 0;
  double t_end = 0.1;

  // Preconstruct a Flux functor
  EdgeFluxCalculator f;

  // Begin the time stepping
  for (double t = t_start; t < t_end; t += dt) {
    // Step forward in time with forward Euler
    hyperbolic_step(mesh, f, t, dt);

    // Update node values with triangle-averaged values
    post_process(mesh);

    // Update the viewer with new node positions
    // HW4B: Need to define node_iterators before these can be used!
    viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
                     CS207::DefaultColor(), NodePosition(), node_map);
    viewer.set_label(t);

    // These lines slow down the animation for small meshes.
    // Feel free to remove them or tweak the constants.
    if (mesh.num_nodes() < 100)
      CS207::sleep(0.05);
  }
auto elapsed = std::chrono::high_resolution_clock::now() - start;
long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();

cout << microseconds << endl;

  return 0;
}
