#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"
#include "Point.hpp"
#include <fstream>
#include <limits>

#include "Mesh.hpp"
#include "shallow_water_ext.cpp"
#include "mass_spring.cpp"
using namespace shallow_water;

/**
    SDL listener, listen to the event of keyboard and mouse.
    @param[in] mouse and keyboard event
    @param[in] initialVel, launchBall.
    @post initiailVel will change on the x, y, z direction based on keyboard event, to increase or decrease the speed.
    @post space or right click will set launchBall to 1
**/

template<typename MeshType>
struct my_listener : public CS207::SDL_Listener {
  void handle(SDL_Event e) {
    switch (e.type) {
      case SDL_MOUSEBUTTONDOWN: {
        if (e.button.button == SDL_BUTTON_RIGHT ) {
          std::cout <<"cannonball launched...approaching targets" <<endl;
          *launchBall = 1;
        }
        }
      case SDL_KEYDOWN: {
        if (e.key.keysym.sym == SDLK_x )
         {
            std::cout << "press x to increase the speed in x direction. speed is" <<(*initialVel).x << endl;
            (*initialVel).x = (*initialVel).x+0.1;
         }
         if (e.key.keysym.sym == SDLK_c)
         {
             std::cout << "press c to decrease the speed in x direction. speed is" <<(*initialVel).x << endl;
            (*initialVel).x = (*initialVel).x-0.1;
         }
         if (e.key.keysym.sym == SDLK_d)
         {
             std::cout << "press d to increase the speed in y direction. speed is" <<(*initialVel).y << endl;
            (*initialVel).y = (*initialVel).y+0.1;
         }
         if (e.key.keysym.sym == SDLK_e)
         {
             std::cout << "press e to decrease the speed in y direction. speed is" <<(*initialVel).y << endl;
            (*initialVel).y = (*initialVel).y-0.1;
         }
         if (e.key.keysym.sym == SDLK_f)
         {
             std::cout << "press f to increase the speed in z direction. speed is" <<(*initialVel).z << endl;
            (*initialVel).z = (*initialVel).z+0.1;
         }
         if (e.key.keysym.sym == SDLK_g)
         {
             std::cout << "press g to increase the speed in z direction. speed is" <<(*initialVel).z << endl;
            (*initialVel).z = (*initialVel).z-0.1;
         }
        if (e.key.keysym.sym == SDLK_UP)
         {
             std::cout << "press up to increase the speed in z direction. speed is" <<(*initialVel).z << endl;
            (*initialVel).z = (*initialVel).z+0.1;
         }
        if (e.key.keysym.sym == SDLK_DOWN)
         {
             std::cout << "press down to decrease the speed in z direction. speed is" <<(*initialVel).z << endl;
            (*initialVel).z = (*initialVel).z-0.1;
         }
        if (e.key.keysym.sym == SDLK_LEFT)
         {
             std::cout << "press left to decrease the speed in x direction and increase y direction. x speed is"
                                    <<(*initialVel).x << " y speed is " << (*initialVel).y << endl;
            (*initialVel).x = (*initialVel).x-0.1;
            (*initialVel).y = (*initialVel).y+0.1;
         }
         if (e.key.keysym.sym == SDLK_RIGHT)
         {
             std::cout << "press right to increase the speed in x direction and decrease y direction. x speed is"
                                    <<(*initialVel).x << " y speed is " << (*initialVel).y << endl;
            (*initialVel).x = (*initialVel).x + 0.1;
            (*initialVel).y = (*initialVel).y-0.1;
         }
        if (e.key.keysym.sym == SDLK_SPACE)
         {
          std::cout <<"cannonball launched...approaching targets" <<endl;
          *launchBall = 1;
         }

      }


       break;
    }
  }

  // constructor
  my_listener(CS207::SDLViewer& viewer, MeshType& mesh, int& launchBall, Point& initialVel) :
                                        viewer_(viewer), mesh_(mesh), launchBall(&launchBall),initialVel(&initialVel)  {};
  private:
   CS207::SDLViewer& viewer_;
   MeshType& mesh_;
   int* launchBall;
   Point* initialVel;
};

struct NodeComparatorX {
  /** Struct/Class of comparator to compare x value of a Node
  * @param[in] two Node objects
  * @param[out] boolean, true if the first Node has the smaller position in x direction than the second node
  */
   template <typename NODE>
   bool operator()(const NODE& t1, const NODE& t2) const {
	return t1.position().x < t2.position().x;
  }
}NodeComparatorX;
struct NodeComparatorY {
 /** Struct/Class of comparator to compare y value of a Node
  * @param[in] two Node objects
  * @param[out] boolean, true if the first Node has the smaller position in y direction than the second node
  */
   template <typename NODE>
   bool operator()(const NODE& t1, const NODE& t2) const {
	return t1.position().y < t2.position().y;
  }
}NodeComparatorY;


int main(int argc, char* argv[])
{
	// Check arguments
	if (argc < 3) {
		std::cerr << "Usage: combine Water_NODES_FILE Water_TRIS_FILE Boat_NODES_FILE Boat_TRIS_FILE Ball_NODES_FILE Ball_TRIS_FILE\n";
		exit(1);
	}

	MeshType mesh; // this mesh is the shallow water mesh
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

    // to calculate the time step for euler step
    auto min_length = *std::min_element(mesh.edge_begin(), mesh.edge_end(), EdgeComparator);

    auto max_h = *std::max_element(mesh.node_begin(), mesh.node_end(), HeightComparator);

    double dt = 0.25 * min_length.length() / (sqrt(grav * max_h.value().h));
    double t_start = 0;
    double t_end = 20;


	auto node_map = viewer.empty_node_map(mesh);
	auto boat_node_map = viewer.empty_node_map(mesh_boat);
	viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
			Node_Color(max_h.value().h), NodePosition(), node_map);
	viewer.add_nodes(mesh_boat.node_begin(), mesh_boat.node_end(),
			CS207::DefaultColor(), BoatNodePosition(), boat_node_map);

	viewer.add_triangles(mesh.tri_begin(), mesh.tri_end(), node_map);

    viewer.add_triangles(mesh_boat.tri_begin(), mesh_boat.tri_end(), boat_node_map);


	// Preconstruct a Flux functor
	EdgeFluxCalculator f;

/**********************************************************
working on the ball now
*************************************************************/
  mass_spring::MeshType ballmesh;
  std::vector<typename mass_spring::MeshType::node_type> ballnodes;

  // Create a nodes_file from the first input argument
  std::ifstream ball_nodes_file(argv[5]);

  // Interpret each line of the nodes_file as a 3D Point and add to the Mesh
  Point ball_p;


  while (CS207::getline_parsed(ball_nodes_file, ball_p))
  {
     ballnodes.push_back(ballmesh.add_node(ball_p/50.0));
  }

  // Calculate the shift point

  auto minx =  *std::min_element(mesh.node_begin(), mesh.node_end(), NodeComparatorX);
  auto miny =  *std::min_element(mesh.node_begin(), mesh.node_end(), NodeComparatorY);
  auto min_h = *std::min_element(mesh.node_begin(), mesh.node_end(), HeightComparator);

  Point Ballcenter = mass_spring::get_center(ballmesh);//Ballcenter/ballmesh.num_nodes();

  Point Balltarget(minx.position().x+.5, miny.position().y-.5, min_h.value().h+0.1);
  Point shift = Balltarget - Ballcenter;
  std::cout << "shift is " << shift.x << shift.y << shift.z << endl;
  for(auto it=ballmesh.node_begin(); it != ballmesh.node_end(); ++ it)
  {
    (*it).set_position((*it).position() + shift) ;
  }

  auto ballmeshOrig = ballmesh; // save a copy to reset the ball

  Point initialVel(2, 4, 3);

  // Create a trianges_file from the second input argument
  std::ifstream triangles_file(argv[6]);

  // Interpret each line of the tets_file as three ints which refer to nodes
  std::array<int,3> ball_t;


  // add in the triangles
  while (CS207::getline_parsed(triangles_file, ball_t))
    for (unsigned i = 1; i < ball_t.size(); ++i)
      for (unsigned j = 0; j < i; ++j)
        for (unsigned k = 0; k < j; ++k)
        {
          ballmesh.add_triangle(ballnodes[ball_t[i]], ballnodes[ball_t[j]], ballnodes[ball_t[k]]);
      }

  // Set masses of nodes equal to 1/N where N is the number of nodes
  // and the initial velocities to 0. Also, get the indexes of
  // the nodes at positions (1,0,0) and (0,0,0)
  for(auto it=ballmesh.node_begin(); it != ballmesh.node_end(); ++ it)
  {
      (*it).value().mass = mass_spring::total_mass/ballmesh.num_nodes();
      (*it).value().velocity = initialVel;
  }

  // Set spring constants for each node equal to spring_const variable
  // and set initial length values equal to lengths of edges prior to
  // running the symplectic Euler steps
  for(auto it=ballmesh.edge_begin(); it != ballmesh.edge_end(); ++ it)
  {
      (*it).value().spring_constant = mass_spring::spring_const;
      (*it).value().initial_length = (*it).length();
  }


  // Set the triangle direction values so that we can determine which
  // way to point normal vectors. This part assumes a convex shape
  Point center = mass_spring::get_center(ballmesh);
  for(auto it=ballmesh.tri_begin(); it != ballmesh.tri_end(); ++ it)
  {
  	//std::cout << (*it).index() << std::endl;
  	mass_spring::set_normal_direction((*it),center);
  }

  // Print out the stats
  std::cout << ballmesh.num_nodes() << " " << ballmesh.num_edges() << std::endl;

  std::cout << "Center: " << mass_spring::get_center(ballmesh) << std::endl;
  // Launch the SDLViewer
  auto ball_node_map = viewer.empty_node_map(ballmesh);

  viewer.add_nodes(ballmesh.node_begin(), ballmesh.node_end(), mass_spring::Node_Color(1),mass_spring::NodePosition(),ball_node_map);
  viewer.add_edges(ballmesh.edge_begin(), ballmesh.edge_end(), ball_node_map);
  int launchBall = 0;

        // add listener
  my_listener<mass_spring::MeshType>* l = new my_listener<mass_spring::MeshType>(viewer,ballmesh, launchBall, initialVel);
  viewer.add_listener(l);

  viewer.center_view();



  // Initialize constraints
  //mass_spring::PlaneConstraint c1;
  //SelfCollisionConstraint c2;
  //auto combined_constraints = make_combined_constraints(c1,c2);
   int reset_ball=0;
   mass_spring::PlaneConstraint c1( min_h.value().h ,  ballmeshOrig, Balltarget, initialVel, launchBall, reset_ball);

   viewer.center_view();


/**********************************************************
END OF working on the ball now
*************************************************************/



	// Begin the time stepping
	for (double t = t_start; t < t_end; t += dt) {
		// Step forward in time with forward Euler
		systime = t;
		hyperbolic_step_ext(mesh, f, t, dt, 1000.0, boat_loc, pressure, &launchBall, &reset_ball, mass_spring::get_center(ballmesh) );

		// Update node values with triangle-averaged values
		post_process(mesh);
		// Update the viewer with new node positions
		// HW4B: Need to define node_iterators before these can be used!
		viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
				Node_Color(max_h.value().h), NodePosition(), node_map);
		viewer.add_nodes(mesh_boat.node_begin(), mesh_boat.node_end(),
				CS207::DefaultColor(), BoatNodePosition(), boat_node_map);

        if(launchBall)
        {

        mass_spring::MassSpringForce ms_force;
        //PressureForce p_force = PressureForce(0.0);
        mass_spring::DampingForce d_force = mass_spring::DampingForce(ballmesh.num_nodes());
        mass_spring::GravityForce g_force;
        //WindForce w_force;

        (void) d_force; // prevents compiler from throwing error for unused variable

        //auto combined_forces = make_combined_force(ms_force, p_force, w_force, g_force);
        auto combined_forces = mass_spring ::make_combined_force(ms_force,  g_force);
        mass_spring::symp_euler_step(ballmesh, t, dt, combined_forces, c1);
        // Update viewer with nodes' new positions
        viewer.add_nodes(ballmesh.node_begin(), ballmesh.node_end(), ball_node_map);
        min_h = *std::min_element(mesh.node_begin(), mesh.node_end(), HeightComparator);

        }
        else if (reset_ball)
        {
            auto it=ballmesh.node_begin();
            auto cit = ballmeshOrig.node_begin();
            for(; it != ballmesh.node_end() && cit!= ballmeshOrig.node_end(); ++ it, ++cit)
            {
                (*it).set_position((*cit).position());
                (*it).value().velocity = initialVel;
            }
            viewer.add_nodes(ballmesh.node_begin(), ballmesh.node_end(), mass_spring::Node_Color(1),mass_spring::NodePosition(), ball_node_map);
            reset_ball=0;
        }

		viewer.set_label(t);

		boat_step(dt, boat_loc);

		// These lines slow down the animation for small meshes.
		// Feel free to remove them or tweak the constants.
		if (mesh.num_nodes() < 100)
			CS207::sleep(0.05);

	}

	return 0;
}
