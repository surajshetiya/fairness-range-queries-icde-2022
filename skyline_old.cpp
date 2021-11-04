/*
Algorithm
===============
1. Data strcuture to create - Store points in r-tree for filtering during query. 
2. Lookup segment tree data structure for finding next point outside range. Or can use sorted per dimension data structure during pre processing for both removal and addition.
3. Preprocessing:
	a. Create R-tree data strcuture.
	b. Create either segment tree or sorted data structure.
4. Query processing:
	a. Create a heap.
	b. Add input range to heap.
	c. Compute neighbours and add each neighbour to heap along with blue and red counts.
		1. Look up dimension wise data structure/ segment tree data structure to find out neighbours.
		2. Get bounding box and lookup r-tree data structure to find points.
		3. Create a skyline out of the points.
		4. Add all neighbouring ranges to heap.
	d. Return from function once a fair range is found.
	
Lookup data structures:
=======================
1. R-tree - https://www.boost.org/doc/libs/1_75_0/libs/geometry/doc/html/geometry/reference/spatial_indexes/boost__geometry__index__rtree.html
2. Segment tree
*/

#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <random>

int main() {
	typedef bool value;
	const size_t dimension = 2;
	typedef boost::geometry::cs::cartesian coordinate_system_type;
	typedef boost::geometry::model::point<coordinate_type, dimension, coordinate_system_type> point;
	typedef boost::geometry::index::rtree<point, value> rtree;
	std::mt19937_64 random_engine;
	auto random_point = [&](){return point(random_distribution_coordinate(random_engine), random_distribution_coordinate(random_engine));};
	auto random_value = [&](){return random_distribution_value(random_engine);};

	rtree r;

	r.print();

	for(auto n = 0; n < 10; ++n)
		r.insert( random_point(), random_value() );
	
	r.print();

	std::cerr<<"---[find]---\n";
	for(const auto& v: r.find(box(point(0.25, 0.25), point(0.75, 0.75))))
		std::cerr<<v<<"\n";

	return 0;
}