#ifndef clipper_utils_hpp
#define clipper_utils_hpp

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include "FIBR3Dapp.hpp"
#include "clipper.hpp"

ClipperLib::Path Polygon_to_ClipperPath(const Polygon &poly);


ClipperLib::Paths Layer_to_ClipperPaths(const Layer &lay);

#endif