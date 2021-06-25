#include "clipper_utils.hpp"

using namespace std;

ClipperLib::Path
Polygon_to_ClipperPath(const Polygon &poly)
{
	ClipperLib::Path retval;
	const std::vector<v3> &points = poly.getPoints();
	for (size_t i=0; i<points.size();++i)
		retval.push_back(ClipperLib::IntPoint((ClipperLib::cInt)points[i].x, (ClipperLib::cInt)points[i].y));
	return retval;
}


ClipperLib::Paths
Layer_to_ClipperPaths(const Layer &lay)
{
	ClipperLib::Paths retval;
	std::vector<Polygon> poly = lay.getPoly();
	for (size_t i=0; i<poly.size(); ++i)
		retval.push_back(Polygon_to_ClipperPath(poly[i]));
	return retval;
}