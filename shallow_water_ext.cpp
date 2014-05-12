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
#include <limits>

#include "Mesh.hpp"

namespace shallow_water{
// Standard gravity (average gravity at Earth's surface) in meters/sec^2
static constexpr double grav = 9.80665;

// global variables for the boat
double boat_height = 0;
double systime = 0;
double vx = 0.05;
double vy = 0.05;

// bathymetry
double function_b(Point& pos){
	return std::sin(pos.x*pos.x+pos.y*pos.y*5)/5;
}
// partial derivatives of b(pos)
double function_dbdx(Point& pos){
	return 2*pos.x*std::cos(pos.x*pos.x+pos.y*pos.y*5)/5;
}

double function_dbdy(Point& pos){
	return 2*pos.y*std::cos(pos.x*pos.x+pos.y*pos.y*5)/5;
}



/** Water column characteristics */
typedef struct QVar {
	double h;   // Height of fluid
	double hu;  // Height times average x velocity of column
	double hv;  // Height times average y velocity of column

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
	// More operators?
	QVar operator+(const QVar& other){
		double h_ = h + other.h;
		double hu_ = hu + other.hu;
		double hv_ = hv + other.hv;
		return QVar(h_,hu_,hv_);
	}

	QVar operator-(const QVar& other){
		double h_ = h - other.h;
		double hu_ = hu - other.hu;
		double hv_ = hv - other.hv;
		return QVar(h_,hu_,hv_);
	}

	QVar& operator=(const QVar& other){
		h = other.h;
		hu = other.hu;
		hv = other.hv;
		return *this;
	}

}QVar;

QVar operator*(const double& lhs, const QVar& rhs){
	double h_ = lhs * rhs.h;
	double hu_ = lhs * rhs.hu;
	double hv_ = lhs * rhs.hv;
	return QVar(h_,hu_,hv_);
}

QVar operator/(const QVar& lhs, const double& rhs) {
	double h_ = lhs.h / rhs;
	double hu_ = lhs.hu / rhs;
	double hv_ = lhs.hv / rhs;
	return QVar(h_, hu_, hv_);
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

	// allow arbitrary force to be applied to each triangle
	QVar operator()(double nx, double ny, double dt,
			const QVar& qk, const QVar& qm, double pressure, double rho) {
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
		double gh2   = 0.5 * grav * (qm.h*qm.h + qk.h*qk.h) + 0.5*(qm.h+qk.h)*pressure/rho;
		/*
		if (pressure > 1){
			std::cout << scale * (wm*qm.hu + wk*qk.hu + gh2*nx) - a * (qm.hu - qk.hu) << "\t"
			<< scale * (wm*qm.hu + wk*qk.hu + 0.5*(qm.h+qk.h)*pressure/rho*nx) - a * (qm.hu - qk.hu) << std::endl;
		}
		*/


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
			// HW4B: You may change this to plot something other than the
			// positions of the nodes
			//return n.position();
			return Point(n.position().x, n.position().y, n.value().h);
		}
};

struct BoatNodePosition {
	template <typename NODE>
		Point operator()(const NODE& n) {
			// HW4B: You may change this to plot something other than the
			// positions of the nodes
			//return n.position();
			//std::cout << n.position().x << " " << n.position().y  << " " << n.position().z + boat_height + 0.2 << std::endl;
			return Point(n.position().x + vx*systime, n.position().y+vy*systime, n.position().z + boat_height);

		}
};

// HW4B: Placeholder for Mesh Type!
// Define NodeData, EdgeData, TriData, etc
typedef Mesh<QVar,double,QVar> MeshType;
typedef unsigned size_type;


template <typename MESH, typename FLUX>
double hyperbolic_step_ext(MESH& m, FLUX& f, double t, double dt, double rho,
                           std::vector<double>& boat_loc, double pressure,
                           int* launchBall, int* reset_ball,  Point Ballcenter) {

	double pressure_applied = 0;
	boat_height = -1e20;

	for (auto it = m.tri_begin(); it!=m.tri_end(); ++it){

        if ((*it).PointInTriangle(Ballcenter))
        {
            if (Ballcenter.z <= (*it).node1().value().h || Ballcenter.z <= (*it).node2().value().h || Ballcenter.z <= (*it).node3().value().h)
                {
                    *launchBall=0;
                    (*reset_ball)=1;

                }
        }

		Point tri_center = ((*it).node1().position() + (*it).node2().position() + (*it).node3().position())/3.0;

		if (isUnderPressure(*it, boat_loc)){
			pressure_applied = pressure;
			//std::cout << "triangle pressure applied\n";
		}
		else{
			pressure_applied = 0;
		}


	QVar total_flux=QVar(0,0,0);
	QVar qm = QVar(0,0,0);
			auto edgetemp = (*it).edge1();
	for (int num = 0; num < 3; num++)
	{
		if (num ==0)
			edgetemp= (*it).edge1();
		else if (num==1)
			edgetemp = (*it).edge2();
		else
			edgetemp = (*it).edge3();

		if (  m.has_neighbor(edgetemp.index()) ) // it has a common triangle
		{
			auto nx =  ((*it).norm_vector(edgetemp)).x;
			auto ny =  ((*it).norm_vector(edgetemp)).y;

			// find the neighbour of a common edge
			for (auto i = m.tri_edge_begin(edgetemp.index()); i != m.tri_edge_end(edgetemp.index()); ++i){
				if (!(*i==*it))
					qm = (*i).value();
			}
			// calculat the total flux
			total_flux = total_flux+f(nx, ny, dt, (*it).value(), qm, pressure_applied, rho) ;
		}
		else{
			// when it doesnt have a neighbour shared with this edge
			auto nx =  ((*it).norm_vector(edgetemp)).x;
			auto ny =  ((*it).norm_vector(edgetemp)).y;
			qm = QVar((*it).value().h, 0, 0 ); // approximation

			total_flux = total_flux+ f(nx, ny, dt, (*it).value(), qm, pressure_applied, rho) ;
		}
	}


		//step through
		(*it).value() = (*it).value() - dt/(*it).area()*total_flux
			+ dt*QVar(0, -grav*((*it).value().h)*function_dbdx(tri_center), -grav*((*it).value().h)*function_dbdy(tri_center));

		if (isUnderPressure(*it, boat_loc)){
			if ((*it).value().h >  boat_height){
				boat_height = (*it).value().h;
			}
			//std::cout << "boat_height value found\n";
		}
	}
	return t + dt;
}

/** Convert the triangle-averaged values to node-averaged values for viewing. */
template <typename MESH>
void post_process(MESH& m) {
  // Translate the triangle-averaged values to node-averaged values
  // Implement Equation 8 from your pseudocode here
  	// iterate through all the nodes
  for ( auto it = m.node_begin(); it!= m.node_end(); ++it)
  {

	QVar sum = QVar(0,0,0);
	double sumTriArea = 0;
	// for each node, iterate through its adjacent triangles

	for (auto adji = m.vertex_begin((*it).index()); adji !=  m.vertex_end((*it).index()); ++ adji)
	{
		auto tri = (*adji);
		sum = sum+ tri.area() * tri.value();
		sumTriArea = sumTriArea+tri.area();
	}

	(*it).value() = sum/sumTriArea; // update nodes value
  }
}

void boat_step(double dt, std::vector<double>& boat_loc) {

	boat_loc[0] += dt*vx;
	boat_loc[1] += dt*vx;
	boat_loc[2] += dt*vy;
	boat_loc[3] += dt*vy;

}

bool isUnderPressure(MeshType::Triangle tri, std::vector<double> boat_loc){


	auto loopnode  = tri.node1();
	for (int i=0;i<3;++i){
        if (i ==0)
			loopnode= tri.node1();
		else if (i==1)
			loopnode = tri.node2();
		else
			loopnode = tri.node3();


		auto pos = loopnode.position();
		if (pos.x > boat_loc[0] && pos.x < boat_loc[1] && pos.y > boat_loc[2] && pos.y < boat_loc[3]){
			return true;
		}
	}
	return false;
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


struct Node_Color{
    int Length_;
    Node_Color(const int Length): Length_(Length){};

    template <typename NODE>
    CS207::Color operator()(NODE& n) //Node_Color(Graph<int>::Node n)
     {
      double c = (double(n.value().h))/1000;
      if (c > 1) c = 1;
      else if (c< 0) c=0;
    return CS207::Color::make_heat(1);
}
};

}
/*




int main(int argc, char* argv[])
{
	// Check arguments
	if (argc < 3) {
		std::cerr << "Usage: shallow_water NODES_FILE TRIS_FILE\n";
		exit(1);
	}

	MeshType mesh;
	// HW4B: Need node_type before this can be used!
	std::vector<typename MeshType::node_type> mesh_node;

	// Read all Points and add them to the Mesh
	std::ifstream nodes_file(argv[1]);
	Point p;
	while (CS207::getline_parsed(nodes_file, p)) {
		// HW4B: Need to implement add_node before this can be used!
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

	// Set initial conditions
	// If dam* used, default to the 3rd initial condition
	// otherwise, allow user to choose between 1st and 2nd initial condition
	std::string input = argv[1];
	std::string dam = "data/dam";
	int initial;
	if (input.compare(0,8,dam) == 0) {
		initial = 3;
	}
	else {
		// ask user choose initial coniditions 1 or 2
		std::cout << "Specify Initial Condition: (Enter 1 or 2)" << std::endl;
		std::cin >> initial;
	}

	// Initial condition 1
	if (initial == 1) {
		for (auto n = mesh.node_begin(); n != mesh.node_end(); ++n) {
			auto x = (*n).position().x;
			auto y = (*n).position().y;
			double h = 1.0 - 0.75 * exp(-80.0*(pow(x-0.75, 2.0) + pow(y, 2.0)));
			mesh.value((*n),QVar( h,0,0));
		}
	}

	// Initial condition 2
	else if (initial == 2) {
		for (auto n = mesh.node_begin(); n != mesh.node_end(); ++n) {
			auto x = (*n).position().x;
			auto y = (*n).position().y;
			double temp = pow(x-0.75, 2.0) + pow(y, 2.0) - pow(0.15, 2.0);
			if (temp < 0) {
                mesh.value((*n),QVar(1.0 + 0.75, 0, 0));


			}
			else {
                mesh.value((*n),QVar(1.0, 0, 0));
				//mesh.update_node_value(QVar(1.0, 0, 0), (*n));
			}
		}
	}

	// Initial condition 3
	else if (initial == 3) {
		for (auto n = mesh.node_begin(); n != mesh.node_end(); ++n) {
			auto x = (*n).position().x;
			if (x < 0) {
                mesh.value((*n),QVar(1.0 + 0.75, 0, 0));

			}
			else {
                mesh.value((*n),QVar(1.0, 0, 0));

				//mesh.update_node_value(QVar(1.0, 0, 0), (*n));
			}
		}
	}

	// Set initial triangle values

  for (auto it = mesh.tri_begin(); it != mesh.tri_end(); ++it ) {
	(*it).value() = ((*it).node1().value() + (*it).node2().value() + (*it).node3().value())/3.0;
  }


	// introduce a boat
	std::vector<double> boat_loc (4,0.0);
	boat_loc[0] = 0.1;
	boat_loc[1] = 0.2;
	boat_loc[2] = 0.1;
	boat_loc[3] = 0.2;
	double boat_weight = 100.0;
	double pressure = boat_weight/((boat_loc[3]-boat_loc[2])*(boat_loc[1]-boat_loc[0]));

	MeshType mesh_boat;
	std::vector<typename MeshType::node_type> mesh_node_boat;
	std::ifstream nodes_boat_file(argv[3]);
	Point p_boat;
	while (CS207::getline_parsed(nodes_boat_file, p_boat)) {
		mesh_node_boat.push_back(mesh_boat.add_node(p_boat));
	}

	std::ifstream tris_boat_file(argv[4]);
	std::array<int,3> t_boat;
	while (CS207::getline_parsed(tris_boat_file, t_boat)) {
		mesh_boat.add_triangle(mesh_node_boat[t_boat[0]], mesh_node_boat[t_boat[1]], mesh_node_boat[t_boat[2]]);
	}




	std::cout << "boat mesh: " << mesh_boat.num_nodes() << " "
		<< mesh_boat.num_edges() << " "
		<< mesh_boat.num_triangles() << std::endl;



	// Perform any needed precomputation
	// Launch the SDLViewer
	CS207::SDLViewer viewer;
	viewer.launch();

	// HW4B: Need to define Mesh::node_type and node/edge iterator
	// before these can be used!

	  auto min_length = *std::min_element(mesh.edge_begin(), mesh.edge_end(), EdgeComparator);

      auto max_h = *std::max_element(mesh.node_begin(), mesh.node_end(), HeightComparator);

      double dt = 0.25 * min_length.length() / (sqrt(grav * max_h.value().h));
      double t_start = 0;
      double t_end = 20;
        std::cout << "dt = " << dt << std::endl;


	auto node_map = viewer.empty_node_map(mesh);
	auto boat_node_map = viewer.empty_node_map(mesh_boat);
	viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
			Node_Color(max_h.value().h), NodePosition(), node_map);
	viewer.add_nodes(mesh_boat.node_begin(), mesh_boat.node_end(),
			CS207::DefaultColor(), BoatNodePosition(), boat_node_map);
	//viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);
	viewer.add_triangles(mesh.tri_begin(), mesh.tri_end(), node_map);
	//viewer.add_edges(mesh_boat.edge_begin(), mesh_boat.edge_end(), boat_node_map);
    viewer.add_triangles(mesh_boat.tri_begin(), mesh_boat.tri_end(), boat_node_map);
	viewer.center_view();



	// HW4B: Timestep
	// CFL stability condition requires dt <= dx / max|velocity|
	// For the shallow water equations with u = v = 0 initial conditions
	//   we can compute the minimum edge length and maximum original water height
	//   to set the time-step
	// Compute the minimum edge length and maximum water height for computing dt
/*	double min_edge_length = std::numeric_limits<double>::max();
	double max_height = 0;
	for (auto it = mesh.edge_begin(); it!= mesh.edge_end(); ++it){
		if (min_edge_length > (*it).value()){
			min_edge_length = (*it).value();
		}
	}
	for (auto it = mesh.tri_begin(); it != mesh.tri_end(); ++it){
		if (max_height < (*it).value().h){
			max_height = (*it).value().h;
		}
	}

	double dt = 0.25 * min_edge_length / (sqrt(grav * max_height));
	std::cout << "dt = " << dt << std::endl;
*/

/*
	// Preconstruct a Flux functor
	EdgeFluxCalculator f;




	// Begin the time stepping
	for (double t = t_start; t < t_end; t += dt) {
		// Step forward in time with forward Euler
		systime = t;
		hyperbolic_step_ext(mesh, f, t, dt, 1000.0, boat_loc, pressure);

		// Update node values with triangle-averaged values
		post_process(mesh);
		// Update the viewer with new node positions
		// HW4B: Need to define node_iterators before these can be used!
		viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
				Node_Color(max_h.value().h), NodePosition(), node_map);
		viewer.add_nodes(mesh_boat.node_begin(), mesh_boat.node_end(),
				CS207::DefaultColor(), BoatNodePosition(), boat_node_map);
		viewer.set_label(t);

		boat_step(dt, boat_loc);

		// These lines slow down the animation for small meshes.
		// Feel free to remove them or tweak the constants.
		if (mesh.num_nodes() < 100)
			CS207::sleep(0.05);

	}

	return 0;
}

*/
