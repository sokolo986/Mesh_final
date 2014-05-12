#include "CS207/Util.hpp"
#include "Point.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <mlpack/core.hpp>
#include <mlpack/methods/neighbor_search/neighbor_search.hpp>
//#include "omp.h"
using namespace std;

template <typename disType>
class Nearest{
    private:
	arma::Mat<size_t> n_;
	arma::mat d_;
	arma::mat new_arma_;
	unsigned num_neighbors_;

    public:
	typedef  mlpack::neighbor::NeighborSearch<NearestNeighborSort, mlpack::metric::EuclideanDistance> nn; 

	Nearest(const unsigned num_neighbors, std::vector<disType> x){
		num_neighbors_ = num_neighbors;
		ofstream myfile;
		myfile.open("nodeData.csv");
		for(unsigned i=0;i<x.size();++i)
			myfile << x[i] << "\n";
		myfile.close();
		mlpack::data::Load("nodeData.csv",new_arma_);
		nn a(new_arma_);
		a.Search(num_neighbors_,n_,d_);
	}

	std::vector<unsigned> n(const unsigned node_idx){
		assert((node_idx*num_neighbors_) < n_.size());
		std::vector<unsigned> v;
		for(unsigned j=0; j<num_neighbors_; ++j)
			v.push_back(n_[num_neighbors_*node_idx+j]);
		return v;
	}

	std::vector<double> d(const unsigned node_idx){
		assert((node_idx*num_neighbors_) < n_.size());
		std::vector<double> v;
		for(unsigned j=0; j<num_neighbors_; ++j)
			v.push_back(d_[num_neighbors_*node_idx+j]);
		return v;
	}

};


int main(){
	std::vector<double> x;
	x.push_back(23.2);
	x.push_back(23.83);
	x.push_back(63.3);
	x.push_back(35.48);
	x.push_back(56.3);
	x.push_back(73.2);
	x.push_back(13.63);
	x.push_back(12.33);
	x.push_back(3.2);
	x.push_back(33.3);
	x.push_back(63.3);

	int num_neighbors = 3;
	unsigned node_idx = 3;

	Nearest<double> my_n(num_neighbors,x);
	std::vector<unsigned> neighbor_idx = my_n.n(node_idx);
	std::vector<double> distances = my_n.d(node_idx);

	for (int i = 0; i < neighbor_idx.size(); ++i)
		cout << i << " neighbor is = " << neighbor_idx[i] << " and dis = " << distances[i] << endl;

	return 0;
}
