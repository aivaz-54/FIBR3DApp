/** 
 * @file  FillingAlgorithms.hpp
 * @author Bruna
 * @date Created on 19 de Julho de 2019, 11:31
 * @brief Methods to create filling patterns and generate corresponding GCode
 */

#ifndef FILLINGALGORITHMS_HPP
#define FILLINGALGORITHMS_HPP

#include "FIBR3Dapp.hpp"
#include "PartsOrientation.hpp"
#include <algorithm>

/**
 * @brief Main method that starts the sliccing and filling process.
 * @details The method gets all possible printing sequences and, for each sequence, goes part by part to calculate the slices and the filling sequence of each slice.
 * @details The filling sequence is generated by \ref slice_filling method.
 * @note If there are more than one mesh to be printed simultaneously, these meshes are combined into one and then passed to \ref slice_filling
 * @param[in] Lf final possible printing sequences
 * @param[in] meshes vector of meshes of each part that compose the entire object
 * @param[in] modelFileName filename of the object
 */
void slice_filling_master(vector<vector<vector<pair<int,int>>>> *L_f, vector<Mesh> *meshes, char *modelFileName);

/**
 * @brief Create the slices for the mesh, then calculates the path to deposit the material by getting parallel segments
 * @param[in] m Mesh with information of the parts to be sliced and printed
 * @param[in] filename name of the file (incuding its directory) to store the generated GCode instructions
 */
void slice_filling(vector<Mesh> *m, vector<pair<int,int>> L, FILE *f, float z_offset);

/**
 * @brief Calculate the length of the segment considering the cartesian distance
 * @param[in] segment LineSegment with points to calculate the distance
 * @return cartesian distance
 */
double lengthSegment(LineSegment * segment);

/**
 * @brief Find the longest segment in a vector of segments
 * @param[in] segments vector of LineSegment to search for the longest
 * @return the longest segment
 */
LineSegment * longestSegment(vector<LineSegment*> segments);

/**
 * @brief Determine the intersection point of two segments
 * @param[in] s1 represents a line segment
 * @param[in] s2 represents a line (points are used to define a line)
 * @param[out] point intersection point
 * @return integer value determining if there is an intersection between the line and the segment
 * @retval 0 no intersection
 * @retval 1 intersection occurs at point
 * @retval 2 segment and line are colinear
 */
int segmentLineIntersect(LineSegment * s1, LineSegment * s2, v3 * point);

/**
 * @brief Create a parallel segment to a specific one at a given distance
 * @param[in] segment line segment from which the parallel is going to be calculated
 * @param[in] d distance of the parallel segment
 * @param[in] direction if set to true, parallels are calculated to the right, left otherwise
 * @return LineSegment representing a parallel segment
 */
LineSegment *findParallelSegment(LineSegment * segment, double d, bool direction);

/**
 * @brief Sorts the points in the vector in order to create a new vector where points have minimal betwee each other
 * @param[in] points vector of unordered v3 points
 * @return vector with ordered v3 points
 */
vector<v3*> fill_sort(vector<v3*> points);

///**
// * @brief Method that creates the filling segments to print the object
// * @param[in] S vector of vectors with LineSegment that contains the set of segments to print. Each position of the outter vector corresponds to a different poygon
// * @param[in] area vector of double with the same number of element as the outter vector of S with the area value of each polygon. Note that a negative area represent a hole
// * @param[in] d double setting the filament spaccing 
// * @param[in] filename the file to store the final result in GCODE
// * @param[in] rotation_y mesh rotation long y axis
// * @param[in] rotation_z mesh rotation long z axis
// */
//void fillingAlgorithm(vector<vector<LineSegment *>> S, vector<double> area, double d, char* filename, float rotation_y, float rotation_z);

/**
 * @brief Method that creates the filling segments to print the object
 * @param[in] SL vector of vector of vectors with LineSegment that contains the set of segments to print. The outter vector represent the layers; each position of the intermediate vector corresponds to a different poygon, which is a vector of LineSegment.
 * @param[in] area_layers vector of vector of double with the same number of element as the outter vector of S with the area value of each polygon. Note that a negative area represent a hole. The outter vector represents the layers whereas the inner vector is the area of eah polygon.
 * @param[in] d double setting the filament spaccing 
 * @param[in] filename the file to store the final result in GCODE
 * @param[in] rotation_y mesh rotation long y axis
 * @param[in] rotation_z mesh rotation long z axis
 */
void fillingAlgorithm(vector<vector<vector<LineSegment *>>> SL, vector<vector<double>> area_layers, double d, FILE *f);

/**
 * @brief Create a vector of LineSegment for a given polygon considering possible holes and the previously calculated points that intersect the horizontal plane
 * @param[in] polygon vector of LineSegment that define a polygon
 * @param[in] holes vector of LineSegment that define the hole(s) of the polygon
 * @param[in] points vector of v3 points that were previously calculated through the parallel segment intersection
 * @return vector of all LineSegment that contain, at least, one point
 */
vector<LineSegment *> createSegments(vector<LineSegment *> polygon, vector<vector<LineSegment*>> holes, vector<v3*> points);

/**
 * @brief Verify if a point is inside a polygon
 * @param[in] polygon vector of LineSegment that defines the polygon
 * @param[in] point v3 with coordinates to check if is inside the polygon
 * @return 1 if the point is inside, 0 otherwise
 */
int pointInsidePolygon(vector<LineSegment*> polygon, v3* point);

/**
 * @brief Create the middle point of the LineSegment
 * @param[in] segment LineSegment from which middle point will be generated
 * @return v3 with the coordinates of the middle point
 */
v3 getMiddlePoint(LineSegment* segment);

/**
 * @brief Calculate the perpendicular segment that crosses the original segment in its middle point
 * @param[in] segment LineSegment from where the perpendicular will be generated
 * @return perpendicular LineSegment
 */
LineSegment * findPerpendicularSegment(LineSegment * segment);

/**
 * @brief Adjust the z coodinate of the line used to generate the parallels
 * @param[in] segment longest segment calculated for the first layer
 * @param[in] h the z coordinate of current layer
 * @return LineSegment with adjusted z coordinate. Remaining coordinates stay unchanged.
 */
LineSegment * adjustH(LineSegment * segment, float h);

/**
 * @brief Verify if a v3 point is on a LineSegment
 * @param[in] p segment starting point
 * @param[in] q intermediate point to check if is inside the segment
 * @param[in] r segment ending point
 * @return true value if point q is on the segment, false otherwise
 */
bool onSegment(v3 * p, v3 * q, v3 * r);

/**
 * @brief Create the bounding box for the current layer
 * @param[in] S vector of polygons for the layer. The polygons are represents as a vector of LineSegment
 * @return vector with the LineSegment that wrap the current layer
 */
vector<LineSegment *> createBoundingLayer(vector<vector<LineSegment *>> S);

/**
 * @brief Verify if the LineSegment intersects the layer bounding box
 * @param[in] box vector of LineSegment that represent the bounding box of the layer
 * @param[in] parallel LineSegment that will be tested
 * @return true if the LineSegment intersects the bounding box, false otherwise
 */
bool intersectBoundingLayer(vector<LineSegment *> box, LineSegment * parallel);

/**
 * @brief Export G0 command to GCode filling file
 * @details Print a G0 GCode instruction stating the position to where the head shoud deposit material.
 * @verbatim
G00 F<speed> X<x> Y<y> Z<z> B<b> C<c> A0.00
@endverbatim
 * @param[in] f FILE* specifying the destination file
 * @param[in] speed float value setting the speed which the head will move
 * @param[in] point IPoint with all information regarding x,y,z
 * @param[in] line integer setting if the line number should be printed in the GCODE (-1 if false, 0 otherwise)
 * @param[in] rotation_y rotation degree value along y-axis
 * @param[in] rotation_z rotation degree value along z-axis
 */
void fill_exportG0command(FILE *f, float speed, v3 * point, int * line, float B, float C);

/**
 * @brief Export G1 command to GCode filling file
 * @details Print a G1 GCode instruction to move the deposition layer to the next starting deposition position.
 * @verbatim
G01 X<x> Y<y> Z<z> B0.00 C0.00 A0.00 
@endverbatim
 * @param[in] f FILE* specifying the destination file
 * @param[in] point IPoint with all information regarding x,y,z position and head rotation
 * @param[in] rate material deposition rate
 * @param[in] line integer setting if the line number should be printed in the GCODE (-1 if false, 0 otherwise)
 * @param[in] rotation_y rotation degree value along y-axis
 * @param[in] rotation_z rotation degree value along z-axis
 * @param[in] extrusion quantity of material that will be extruded during the current path
 */
void fill_exportG1command(FILE *f, v3 * point, float rate, int * line, float B, float C, float extrusion);

/**
 * @brief Generate the GCode file with the instructions to print the object
 * @param[in] path vector of LineSegment with the location to deposit the printing material
 * @param[in] filename file (with directory) specifying the file to create or append the layer GCode instructions
 * @param[in] line print the instruction line number. 1 - print, 0 otherwise
 * @param[in] rotation_y rotation degree value along y-axis
 * @param[in] rotation_z rotation degree value along z-axis
 * @return 0 if the instructions were correctly written to the file
 */
int fill_exportGCode5Axes_layer(vector<LineSegment*> path, FILE *f, int * line);

/**
 * @brief Order the segments that compose each polygon in order to deposit the material from one extremity of the layer to the other
 * @param[in] segments vector of LineSegment that define each polygon of the layer
 * @return ordered vector of LineSegment
 */
vector<LineSegment*> orderSegmentsClosest(vector<vector<LineSegment*>> *segments);

/**
 * @brief Create a string with the directory and filename to create the filling GCode sequence
 * @param[in] stl_file file with the original stl_file
 * @param[in] directory location of stl_file
 * @param[in] extension the extension to add to the new file (ex. .gcode)
 * @return char* with the file and directory to store the results
 */
char* createDirectoryFile_filling(char* stl_file, char* directory, char* extension);

/**
 * @brief Calcultes the cartesian distance between two v3 points
 * @param[in] from v3 representing the initial point
 * @param[in] to v3 representing the final point
 * @return cartesian distance
 */
double fill_distance(v3 from, v3 to);

/**
 * @brief Remove duplicated points in the vector
 * @param[in,out] points vector of v3 points
 */
void removeRepeatedPoints(vector<v3*> * points);

void printPolygon_Scad(char *filename, vector<vector<LineSegment *>> polygons);
#endif /* FILLINGALGORITHMS_HPP */
