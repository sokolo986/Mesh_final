/**
 * @file poisson.cpp
 * Test script for treating the Graph as a MTL Matrix
 * and solving a Poisson equation.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles.
 * Second file: Eges (one per line) defined by 2 indices into the point list
 *              of the first file.
 *
 * Launches an SDLViewer to visualize the solution.
 */

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"

#include "Graph.hpp"
#include <math.h>
#include <fstream>
#include "BoundingBox.hpp"
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include "CS207/Color.hpp"
using namespace std;

typedef Graph<double, double> GraphType;
typedef GraphType::node_type Node;
typedef GraphType::edge_type Edge;
typedef GraphType::size_type size_type;


bool Boundary(const Point& p, double& value){
/* define a boundary condition for poisson process
 * @param[in] a valid point p, double value
 * @pre a valid point p
 * @post return a boolean and assign the value to a double
		if infinity-norm of p ==1 , return true, value = 0
		if infinity-norm of p-(+-0.6,+-0.6,0)  < 0.2, return true, value = -0.2
		if x is in the bounding box of (-0.6,-0.2,-1) and (0.6,0.2,1) return true, value = 1
*/	
	const Point B=Point(0.6,0.2,1.0);
	const Point A=Point(-0.6,-0.2,-1.0);
	BoundingBox	BB(A,B);
		
		if (norm_inf(p) == 1){
			value = 0.0;
			return true;
		}
		else if ((norm_inf( p - Point(0.6,0.6,0)) < 0.2) ||
				(norm_inf( p - Point(-0.6,0.6,0)) < 0.2) ||
				(norm_inf( p - Point(0.6,-0.6,0)) < 0.2) ||
				(norm_inf( p - Point(-0.6,-0.6,0)) < 0.2) ){
			value = -0.2;
			return true;
		}
		else if  (BB.contains(p)){
			value = 1.0;
			return true;
		}
		else
			return false;
	
}




struct GraphSymmetricMatrix {
	/** Define a GraphSymmetricMatrix using Graph concept
	/** Helper function to perform multiplication . Allows for delayed
	* evaluation of results and various assignment operations such
	* as += , -= , and =.
	* @pre @a size (v) == size (w) */
	GraphSymmetricMatrix(GraphType* g) : g_(g) {}
		
	template <typename VectorIn , typename VectorOut , typename Assign >
	void mult ( const VectorIn & v, VectorOut & w, Assign ) const{
		assert(size(v) == size(w));
		size_type m = g_->size();
		for (size_type i = 0; i<m; i++){
			double val = 0.0;
			auto nodea = g_->node(i);
			auto position = nodea.position();
			double temp;
			if (Boundary(position,temp) )
				val += v[i];
			else{
				val += -(1.0)*nodea.degree() * v[i];
				for (auto it = nodea.edge_begin(); it != nodea.edge_end();++it){	
					// get adj node
					auto adj = nodea;
						if (nodea == (*it).node1())
							adj = (*it).node2();
						else
							adj = (*it).node1();
					// check if the adj node is in the boundary
					if(!Boundary(adj.position(), temp))
						val += v[size_type(adj.index())];
				}
			}
			
			Assign::apply(w[i],val);
		}
	}


	template <typename VectorIn >
	mtl::vector::mat_cvec_multiplier < GraphSymmetricMatrix , VectorIn >
	operator *( const VectorIn & v) const {
		return mtl :: vector :: mat_cvec_multiplier<GraphSymmetricMatrix , VectorIn >(* this , v);
	}

	GraphType* g_;
};





/** The number of elements in the matrix . */
inline std :: size_t size ( const GraphSymmetricMatrix & A){return A.g_->size() * A.g_->size();}

/** The number of rows in the matrix . */
inline std :: size_t num_rows ( const GraphSymmetricMatrix & A){return A.g_->size();}

/** The number of columns in the matrix . */
inline std :: size_t num_cols ( const GraphSymmetricMatrix & A){return A.g_->size();}

/*
Name space used for interface with MTL
*/
namespace mtl { 

    template <>
    struct Collection<GraphSymmetricMatrix >{
        typedef double value_type;
        typedef unsigned    size_type;
    };

	
    namespace ashape {
        template <> 
		struct ashape_aux<GraphSymmetricMatrix >{
			typedef nonscal type;    
		};
    }
}




template <typename VectorOut>
void constructBVectorFunc (const GraphType& g, VectorOut & w ){
/* Construct a valid B vector
* @param[in,out] a graph @g, a vector @w
* @pre length of w = g.size()
* @post return a valid B vector. B(i) = g(x_i) if node i is on the boundary
                                 B(i) = h*h * f(x_i) - \sum (g(xj)) for all the adjacent nodes of i if j is on boundary
								 function f = 5cos(||xi||_1)
								 function g is the boundary function
*/


	size_type m = g.size();
	
	for (size_type i = 0; i<m; i++){
		double val = 0;
		double adjval = 0;
		auto nodea = g.node(i);
		auto position = nodea.position();
		if (Boundary(position, val) ){
			w[i] = val;  // check if node i is at boundary
		}
		else{
			
			auto it = nodea.edge_begin();
			val = 5.0*cos(norm_1(position)) * double((*it).length()) * double((*it).length());  // calculate the part of 5cos(||x||_1) * h *h
			
				
			for (auto it = nodea.edge_begin(); it != nodea.edge_end();++it){
				// here we get the valid adj node
				auto adj = nodea;
				if (nodea == (*it).node1())
					adj = (*it).node2();
				else
					adj = (*it).node1();
				
				if ( Boundary(adj.position(),adjval))  // if node j is on boundary, substract the value from val
					val += -adjval;
			}
			w[i] = val;
		}
	}
}

template<typename VectorIn>
struct vector_Color{
/*
* @param[in] vec a valid vector with [] operator, a node object
* @ pre a valid vector. length of vector >= max(node.index())
* @ post returns a valid Color type using the information of a node and its correspoding information in the vector
*/
    VectorIn vec_;
	double max,min;
    vector_Color(VectorIn& vec): vec_(vec), max(*std::max_element(vec.begin(),vec.end())), min(*std::min_element(vec.begin(),vec.end())){};
    template <typename NODE>
    CS207::Color operator()(NODE& n)
    {
		auto sum = (vec_[n.index()] -min)/ (max-min+0.01);
		return CS207::Color::make_heat(sqrt(sum));
	}
};

template<typename VectorIn>
struct NodePositionFunc{
/*
* @param[in] vec a valid vector with [] operator, a node object
* @ pre a valid vector. length of vector >= max(node.index())
* @ post returns node position with x, y from the old node position, and z from vec value
*/
	VectorIn vec_;
	NodePositionFunc(VectorIn& vec): vec_(vec){};
	template<typename NODE>
	Point operator()(NODE& n)
	{
		return Point(n.position().x, n.position().y, vec_[n.index()]); 
	}

};

namespace itl {
	template <typename VectorIn, typename Nodemaptype,  class Real, class OStream = std::ostream>
	/*
	* inheritance of class cyclic_iteration where we can out put the visualization on the fly
	* @param[in] GraphType g, which is a graph. CS207::SDLViewer viewer, which is the SDL viewer.
				Nodemaptype nodemap, which is the nodemap used in SDL viewer
				VectorIn x, which is the result of pde x 
				the result of input is the same as super class
	  @pre  a valid graph, viewr, nodemap and x value
	  @post visualize the solution while printing everything else from super class
	*/
	class visual_iteration : public cyclic_iteration<Real, OStream> 
	{
		typedef cyclic_iteration<Real, OStream> super;
		typedef visual_iteration self;
		
		public:
	    template <class Vector>
		// Constructor
		visual_iteration(const Vector& r0, int max_iter_, Real tol_,
						 GraphType& g, CS207::SDLViewer& viewer,
						Nodemaptype& nodemap, VectorIn& x, 
						Real atol_ = Real(0), int cycle_ = 100, OStream& out = std::cout)
						: super(r0, max_iter_, tol_, atol_ , cycle_ , out),
						g_(&g), viewer_(&viewer), nodemap_(&nodemap),x(&x) {}
						

		bool finished() { return super::finished(); }
		template <typename T>
        bool finished(const T& r) 
        {
           bool ret= super::finished(r);
           visual();
           return ret;
        }
		void visual(){
			viewer_->add_nodes(g_->node_begin(),g_->node_end(),vector_Color<VectorIn>(*x),NodePositionFunc<VectorIn>(*x),*nodemap_);
		}
		protected:
			Nodemaptype*       nodemap_;
			CS207::SDLViewer*  viewer_;
			GraphType*         g_;
			VectorIn*          x;
	};

}// namespace itl



int main(int argc, char** argv)
{
	// Check arguments
	if (argc < 2) {
	std::cerr << "Usage: " << argv[0] << " NODES_FILE EDGES_FILE\n";
	exit(1);
	}

	GraphType graph;
	std::vector<Node> nodes;

	// Create a nodes_file from the first input argument
	std::ifstream nodes_file(argv[1]);
	// Interpret each line of the nodes_file as a 3D Point and add to the Graph
	Point p;
	while (CS207::getline_parsed(nodes_file, p))
	nodes.push_back(graph.add_node(p));

	// Create a tets_file from the second input argument
	std::ifstream tets_file(argv[2]);
	// Interpret each line of the tets_file as four ints which refer to nodes
	std::array<int,2> t;
	while (CS207::getline_parsed(tets_file, t))
		graph.add_edge(nodes[t[0]], nodes[t[1]]);
	// Print out the stats
	std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;
	// Launch the SDLViewer
	CS207::SDLViewer viewer;
	//typedef std::map<Node, unsigned> itl::nodemaptype;
 	viewer.launch();
	auto node_map=viewer.empty_node_map(graph);
	viewer.add_nodes(graph.node_begin(),graph.node_end(),node_map);
	viewer.add_edges(graph.edge_begin(),graph.edge_end(),node_map);
    viewer.center_view();	
	
    GraphSymmetricMatrix    A(&graph);  // construct A using graph
	typedef mtl::dense_vector<double> dvector;
	dvector b(graph.num_nodes()),x(graph.num_nodes(),0.0);    // construct vector b and x
	constructBVectorFunc< dvector >(graph, b);  // update vector b value using constructBVectorFunc function
    //itl::cyclic_iteration<double>             iter(b, 500, 1.e-10, 0.0, 50);
    itl::visual_iteration<dvector,std::map<Node, unsigned> ,double>  iter(b, 500, 1.e-10, graph, viewer, node_map, x, 0.0, 50); // construct visual_interation

					
	itl::cg(A, x, b, iter);  // call the cg function under itl namespace
	viewer.add_nodes(graph.node_begin(),graph.node_end(),vector_Color<dvector>(x),NodePositionFunc<dvector>(x),node_map); // final result

    return 0;
}
