/** 
 * @file Inspection.hpp
 * @author Bruna Ramos
 * @date Created on 31 de Outubro de 2018, 14:47
 * @brief Methods and procedures related to the inspection of the object
 */

#ifndef INSPECTION_HPP
#define INSPECTION_HPP

#include <vector>
#include "FIBR3Dapp.hpp"
#include <exception>
#include <map>
#include <algorithm>

/** @def PI
 *  @brief PI value.
 */
#define PI 3.141592653589793f

/** @def RAD_TO_DEG(x)
    \brief A macro that returns the radian value of \a x degrees.
 */
#define RAD_TO_DEG(x) (x* 180 / PI)

/**
 * @brief Structure to represent the inpection details, including details of the inspection camera
 */
struct inspection {
    float camLeftYSize = 190.05f; /**< @brief left Y size from the camera (0,0,0) point */
    float camRightYSize = 333.45f; /**< @brief right Y size from the camera (0,0,0) point */
    float camXSize = 280.55f; /**< @brief X size from the camera (0,0,0) point */
    float camUpperZSize = 80.50f; /**< @brief upper Z size from the camera (0,0,0) point */
    float camDownZSize = 80.50f; /**< @brief down Z size from the camera (0,0,0) point */
    float distance_to_part = 50.0f; /**< @brief distance between inspection camera and inspected object */
    float slicing = 50.0f; /**< @brief distance between slices */
    float sampling_distance = 50.0f; /**< @brief distance between inspection points */
    float speed = 200.0f; /**< @brief speed of the inspection */
    v3 min_axis = v3(0.0f, 0.0f, 0.0f); /**< @brief minimum v3 point of the inspection plate */
    v3 max_axis = v3(1015.0f, 1500.0f, 1000.0f); /**< @brief maximum v3 point of the inspection plate */

    float adjust_z = 285.5f; /**< @brief adjust over the z axe the difference between the inspection point and the GCode instructions */
    float adjust_x = -149.0f; /**< @brief adjust over the x axe the difference between the inspection point and the GCode instructions */
    int printGcodeLines = 0; /**< @brief print line number in GCode. Default is 0 (False), 1 otherwise*/
    int MIP = 0; /**< @brief create the inspection path using MIP approach. Default is 0 (False), 1 otherwise */
    int CombF = 0; /**< @brief create the inspection path using CombF approach. Default is 0 (False), 1 otherwise */
    int k_NNSH = 0; /**< @brief create the inspection path using k_NNSH approach. Default is 0 (False), k otherwise */
    int k_NNGH = 0; /**< @brief create the inspection path using k_NNGH approach. Default is 0 (False), k otherwise */
    int NNH = 0; /**< @brief create the inspection path using NNH approach. Default is 0 (False), 1 otherwise */
    bool action_print = false; /**< @brief print to console. Default is False*/
    int alpha = 2; /**< @brief set the maximum number of points to create sub-paths. @note must be greater than 1*/
    int graph_directed = 1; /**< @brief defines if the graph is to be considered oriented or not. Value 1 sets an oriented graph, 0 othewise.*/

    bool debug = true; /**< @brief used to set the debug level. When this value is true then the application prints more execution details */
    vector<SLayer> SLayers ;/**< @brief store the vector of SLayer created during the slicing methods */
};

/**
 * @class IPoint
 * @brief A class for the Inspection Point (IPoint)
 * @author Bruna Ramos
 */
class IPoint {
public:

    /**
     * @brief Create an empty IPoint
     */
    IPoint() {
    }

    /**
     * @brief Create an Inspection Point with values
     * @param[in] point v3 with the coordinates of the inspection point
     * @param[in] normal v3 with the normal values
     * @param[in] index integer with an id for the inspection point
     */
    IPoint(const v3 point, const v3 normal, int index) {
        this->point = point;
        this->normal = normal;
        this->collide = 0;
        this->MrotZ_inspection();
        this->MrotY_inspection();
        this->index = index;
    }

    /**
    @brief Calculate the Z rotation for inspection and updates the c_angle value
    @details The inspection machine has an -180 to 180-degree rotation along the Z axis. This function calculates the rotation angle in degrees considering the x and y coordinates of the normal vector associated to the inspection point.
    @image html MrotZ.PNG "Figure: Rotation along the Z axis" width=700px
    @image latex MrotZ.PNG "Rotation along the Z axis" width=10cm
    @note For the particular case where \a i or \a j is equal to \a 0, the return value is 0.
     */
    void MrotZ_inspection() {
        v3 inverted = this->normal * -1;
        if (inverted.x != 0 || inverted.y != 0) {
            c_angle = RAD_TO_DEG(atan2(inverted.y, inverted.x));
        } else {
            c_angle = 0;
        }
    }

    /**
    @brief Calculate the Y rotation for inspection and updates the b_angle value
    @details The inspection machine has an -30 to 30-degree rotation along the Y axis. This function calculates the rotation angle in degrees considering the z coordinate of the normal vector associated to the inspection point.
    @image html MrotY.PNG "Figure: Rotation along the Y axis" width=700px
    @image latex MrotY.PNG "Rotation along the Y axis" width=10cm
     */
    void MrotY_inspection() {
        v3 inverted = this->normal * -1;
        b_angle = RAD_TO_DEG(-asin(inverted.z));
        if (b_angle > 30) {
            b_angle = 30;
        } else if (b_angle < -30) {
            b_angle = -30;
        }
    }

    /**
     * @brief Calculate the euclidian distance between two inspection points
     * @param[in] ip IPoint to which distance will be measured
     * @return float with the euclidian distance
     */
    float evalTo(IPoint * ip) {
        if(ip->index == -10 || this->index == -10){
            return 0;
        }
        return sqrt(pow(ip->point.x - this->point.x, 2) + pow(ip->point.y - this->point.y, 2) + pow(ip->point.z - this->point.z, 2));
    }

    /**
     * @brief Prints the IPoint information in a human readable format
     * @details Example
     * @verbatim
      [collision][id:index] p(x,y,z) n(x,y,z) B(degrees) C(degrees)
      @endverbatim
     */
    void print() const {
        printf("[%d][id:%d] p(%f,%f,%f) n(%f,%f,%f) B(%f) C(%f)\n", collide, index, point.x, point.y, point.z, normal.x, normal.y, normal.z, b_angle, c_angle);
    }

    /**
     * @brief Print the IPoint in a SCad format
     * @details This method prints the point and a line that simulates the normal vector of the IPoint
     */
    void printSCad() {
        v3 temp2 = point + normal * 50;
        printf("point([%f,%f,%f],5);//inpection point\n", point.x, point.y, point.z);
        printf("line([%f,%f,%f],[%f,%f,%f],3);//inspection point segment\n", point.x, point.y, point.z, temp2.x, temp2.y, temp2.z);
    }
    
    /**
     * @brief Print the IPoint in a SCad format
     * @details This method prints the point and a line that simulates the normal vector of the IPoint. It may also include the xyz coordinates next to the point
     * @param i integer setting if the point is to be printed. 1 is true, 0 otherwise
     */
    void printSCad(int i) {
        v3 temp2 = point + normal * 50;
        printf("point([%f,%f,%f],5,%d);//inpection point\n", point.x, point.y, point.z,i);
        printf("line([%f,%f,%f],[%f,%f,%f],3);//inspection point segment\n", point.x, point.y, point.z, temp2.x, temp2.y, temp2.z);
    }


    /**
     * @brief v3 with the coordinates of the IPoint
     */
    v3 point;

    /** 
     * @brief v3 with the normal values of the IPoint 
     */
    v3 normal;

    /**
     * @brief Collision status of the IPoint.
     * @details Collide is set to zero (0) when no collision occurs, 1 if there is a collision.
     */
    int collide;

    /**
     * @brief Rotation angle over the Y-axe (B angle)
     */
    float b_angle;

    /**
     * @brief Rotation angle over the Z-axe (C angle)
     */
    float c_angle;

    /**
     * @brief Id of the IPoint
     */
    int index;

    /**
     * @brief Position of the IPoint in the layer polygon's vector
     */
    int polygon;

    /**
     * @brief Position of the IPoint in the ISubPath vector
     */
    int subPath;

    /**
     * @brief Position of the IPoint in the IPoint vector
     */
    int position;
    
    /**
     * @brief List of IPoint* setting the possible links between the current IPoint to the next
     */
    vector<IPoint*> adjacentList;
};

/**
 * @class ISubPath
 * @brief A class for the Inspection Sub Path (ISubPath)
 * @author Bruna Ramos
 */
class ISubPath {
public:

    /**
     * @brief Create an empty ISubPath
     */
    ISubPath() {
    }

    /**
     * @brief Get the number of IPoint in the subpath
     * @return size_t with the size of subpath
     */
    size_t size() const {
        return this->subPath.size();
    }

    /**
     * @brief Add a new IPoint to the subpath
     * @param[in] ip IPoint to add to the current subpath
     */
    void push_back(IPoint* ip) {
        this->subPath.push_back(ip);
    }

    /**
     * @brief Print the Subpath
     * @details Initially prints all IPoint in human readable form and then prints the point as segments in SCad format
     */
    void print() const {
        for (int i = 0; i < subPath.size(); i++) {
            printf("#%d ", i);
            subPath.at(i)->print();
        }

        for (int i = 0; i < subPath.size() - 1; i++) {
            IPoint* orig = subPath.at(i);
            IPoint* dest = subPath.at(i + 1);
            printf("line([%f,%f,%f],[%f,%f,%f],1);\n", orig->point.x, orig->point.y, orig->point.z, dest->point.x, dest->point.y, dest->point.z);
        }
    }

    /**
     * @brief Remove an IPoint from the subpath vector
     * @param[in] position integer with the position of the IPoint to be removed from the subpath
     */
    void removeIPoint(const int position) {
        this->subPath.erase(this->subPath.begin() + position);
    }

    /**
     * @brief vector of IPoint that represent an inspection sub path
     */
    vector<IPoint*> subPath;

};

/**
 * @class IPoly
 * @brief A class for the Inspection polygon (IPoly)
 * @author Bruna Ramos
 */
class IPoly {
public:

    /**
     * @brief Print the polygon as a set of lines in SCad format
     */
    void printSCad() const{
        printf("//polygon\n");
        for(auto subPath:this->subPaths){
            for(int i = 0;i<subPath->subPath.size()-1;i++){
                IPoint * orig = subPath->subPath.at(i);
                IPoint * dest = subPath->subPath.at(i+1);
                printf("line([%f,%f,%f],[%f,%f,%f],1);\n", orig->point.x, orig->point.y, orig->point.z, dest->point.x, dest->point.y, dest->point.z);
            }
        }
    }
  
    /**
     * @brief Print the current polygon
     * @details Prints, at an initial step, the IPoint in an human readable format and then all subpaths details
     */
    void print() const {

        printf("Polygon Inspection SubPaths (#%d)\n", (int) this->subPaths.size());
        for (int i = 0; i < subPaths.size(); i++) {
            printf("#subpath%d\n", i);
            subPaths.at(i)->print();
        }
    }

    /**
     * @brief Calculate the total number of valid IPoint that compose the IPoly
     * @details Sums the number of IPoint for each ISubPath
     * @return int value with the sum of valid IPoint
     */
    int getValidIPointSize(){
        int count = 0;
        for(int i=0;i<subPaths.size();i++){
            count += subPaths.at(i)->size();
        }
        return count;
    }

    /**
     * @brief Add an IPoint to the begining of the IPoint vector for the polygon
     * @param p IPoint to be added
     */
    void push_first(IPoint * p){
        this->inspectionPoints.insert(this->inspectionPoints.begin(),p);
    }
    
    /**
     * @brief Get the biggest index from all IPoint that made the polygon
     * @return the biggest index
     */
    int getBiggestIndex(){
        if(inspectionPoints.front()->index>inspectionPoints.back()->index){
            return inspectionPoints.front()->index;
        }
        return inspectionPoints.back()->index;
    }
    
    /**
     * @brief vector of all IPoint that constitute the inspection polygon
     * @note this includes points in which a collision may occur
     */
    vector<IPoint*> inspectionPoints;

    /**
     * @brief vector of all ISubPath
     * @details IPoint are sorted according to the inspection sub path and considering their possible collisions
     */
    vector<ISubPath*> subPaths;
};

/**
 * @class IPath
 * @author Bruna Ramos
 * @brief A class for the Inspection Path (IPath)
 */
class IPath {
public:

    /**
     * @brief Empty constructor
     */
    IPath() {
    }

    /**
     * @brief Clone an inspection path
     * @param[in] path IPath to be cloned
     */
    IPath(const IPath& path) {
        this->eval = path.eval;
        this->path = vector<IPoint*>(path.path);
    }

    /**
     * @brief update IPath eval value considering that the path is closed (return to the same point)
     */
    void calculateEvalClosedPath(){
        eval = 0.0f;
        for(int i=0;i<path.size()-1;i++){
            eval += path.at(i)->evalTo(path.at(i+1));
        }
        eval += path.back()->evalTo(path.front());
    }
    
    /**
     * @brief Update the eval value for the IPath considering the euclidian distance between IPoint
     */
    void calculateEval(){
        eval = 0.0f;
        for(int i=0;i<path.size()-1;i++){
            eval += path.at(i)->evalTo(path.at(i+1));
        }
    }

    /**
     * @brief Print the current path in SCAD format
     * @details This method prints all lines and, at the end, the evaluation i.e. the total time for the IPath
     */
    void print() const {
        for (int i = 0; i < path.size() - 1; i++) {
            IPoint *from = path.at(i);
            IPoint *to = path.at(i + 1);
            printf("line([%f,%f,%f],[%f,%f,%f],1);\n", from->point.x, from->point.y, from->point.z, to->point.x, to->point.y, to->point.z);
        }
        printf("//eval: %f\n", eval);
    }

    /**
     * @brief Create a file in SCAD format with the current path
     * @param[in] filename char* with the name of the file
     * @param[in] nextLayer integer used to define some printing details. Please refer to notes
     * @note Values for nextLayer are 
      <table style="width:25%">
        <tr><th>Value</th><th>Result</th></tr>
        <tr><th>1</th><th>Print the line '//up to the next layer'</th></tr>
        <tr><th>2</th><th>Print the line 'color(colorPath){'</th></tr>
        <tr><th>3</th><th>Print the line '}'</th></tr>
        <tr><th>4</th><th>Print both lines defines in point 2 and 3</th></tr>
      </table>
     */
    void print(char* filename, int nextLayer) const {
        FILE * f = NULL;
        fopen_s(&f, filename, "a");
        if (!f) {
            return;
        }

        if (nextLayer == 1) {
            fprintf(f, "//up to the next layer\n");
        }
        if (nextLayer == 2 || nextLayer == 4) {
            fprintf(f, "color(colorPath){\n");
        }

        for (int i = 0; i < path.size() - 1; i++) {
            IPoint *from = path.at(i);
            IPoint *to = path.at(i + 1);
            fprintf(f, "line([%f,%f,%f],[%f,%f,%f],N_path);\n", from->point.x, from->point.y, from->point.z, to->point.x, to->point.y, to->point.z);
        }
        fprintf(f, "//eval: %f\n", eval);
        if (nextLayer == 3 || nextLayer == 4) {
            fprintf(f, "}\n");
        }
        fclose(f);
    }

    /**
     * @brief Print the path using the IPoint index to define the connections
     */
    void print2() const {
        for (int i = 0; i < path.size(); i++) {
            printf("%d\t", path.at(i)->index);
        }
        printf("\n");
    }

    /**
     * @brief Print the path using all information in the IPoint
     */
    void print3() const {
        for (int i = 0; i < path.size(); i++) {
            path.at(i)->print();
        }
        printf("\n");
    }

    /**
     * @brief Store the current path. Note that the IPoint are considered to be in order
     */
    vector<IPoint*> path;

    /**
     * @brief Defines the cost (evaluation) of the path
     */
    float eval;
};

/**
 * @class ILayer
 * @author Bruna Ramos
 * @brief A class for the Inspection layer (ILayer)
 */
class ILayer {
public:

    /**
     * @brief Print the layer as a set of inspection polygons where each polygon is printed in an human readable format
     */
    void print() const {
        for (int i = 0; i < inspectionPolygons.size(); i++) {
            printf("Polygon %d\n", i);
            inspectionPolygons.at(i)->print();
        }
    }
    
    /**
     * @brief Print the ILayer in SCAD format.
     * @details print all polygons and the equivalent IPoint in SCad format
     */
    void printSCad() const {
        for(auto poly : inspectionPolygons){
            poly->printSCad();
            for(auto isp: poly->subPaths){
                for(auto p : isp->subPath){
                    p->printSCad();
                }
            }
        }
    }
    
    /**
     * @brief Print the ILayer in SCAD format.
     * @details print all polygons and the equivalent IPoint in SCad format. Equivalent to printSCad() but here one may set if the xyz coordinates need to be printed.
     * @param[in] i set if the xyz IPoint coordinates should be printed near the point
     */
    void printSCad(int i) const {
        for(auto poly : inspectionPolygons){
            poly->printSCad();
            for(auto isp: poly->subPaths){
                for(auto p : isp->subPath){
                    p->printSCad(i);
                }
            }
        }
    }

    /**
     * @brief get the total number of valid IPoint
     * @details counts the number of valid (without collisions) IPoint
     * @return total number of valid points
     */
    int getValidIPointSize(){
        int count = 0;
        for(int i=0;i<this->inspectionPolygons.size();i++){
            count += inspectionPolygons.at(i)->getValidIPointSize();
        }
        return count;
    }
    
    /**
     * @brief get the biggest index of all IPoint that made the ILayer
     * @return biggest index value (integer)
     */
    int getBiggestIndex(){
        int biggest = -1;
        for(IPoly* p : inspectionPolygons){
            int temp = p->getBiggestIndex();
            if(temp > biggest){
                biggest = temp;
            }
        }
        return biggest;
    }
    
    /**
     * @brief vector of inspection polygons (IPoly)
     */
    vector<IPoly*> inspectionPolygons;
    
    /**
     * @brief IPath used to store the calculated trajectory for the current ILayer
     */
    IPath* layerPath;
        
};

/**
 * @class ILayers
 * @author Bruna Ramos
 * @brief A class for the Inspection Layers (ILayers)
 */
class ILayers {
public:

    /**
     * @brief print ILayers information in SCAD format
     * @param[in] i set if the xyz coordinates should be printed near the point. 1 sets to true, 0 otherwise.
     */
    void printSCad(int i){
        int l=1;
        for(auto layer : this->inspectionLayers){
            printf("//Layer %d\n",l++);
            layer->printSCad(i);
        }
    }
    
    /**
     * @brief Print a layer information by invoking the print command for the selected ILayer
     * @param[in] layer integer setting the ILayer to be printed
     */
    void print(const int layer) const {
        printf("Layer %i\n", layer);
        inspectionLayers.at(layer)->print();
    }

    /**
     * @brief Print all layers by calling the print function for each ILayer
     */
    void print() const {
        for (int i = 0; i<this->inspectionLayers.size(); i++) {
            printf("====== Layer %d =======\n", i);
            inspectionLayers.at(i)->print();
        }
    }

    /**
     * @brief Print the complete path, i.e., each path for all ILayer and the path that connects them
     * @details To separate each layer there is a commented line with '// up to the next layer' where the following line is the connection between layers
     */
    void printCompletePath() const {
        for (int i = 0; i<this->inspectionLayers.size(); i++) {
            this->inspectionLayers.at(i)->layerPath->print();
            if (i != this->inspectionLayers.size() - 1) {
                printf("// up to the next layer\n");
                IPoint *from = inspectionLayers.at(i)->layerPath->path.back();
                IPoint *to = inspectionLayers.at(i + 1)->layerPath->path.front();
                printf("line([%f,%f,%f],[%f,%f,%f],1);\n", from->point.x, from->point.y, from->point.z, to->point.x, to->point.y, to->point.z);
            }
        }
    }

    /**
     * @brief calculate if the IPath is complete, i.e., its size is not 0 or the value is not null
     * @return value defining if the path is complete
     * @retval 1 path is complete
     * @retval 0 otherwise
     */
    int isPathComplete(){
        for(auto layer : this->inspectionLayers){
            if(layer->layerPath == NULL || layer->layerPath->path.size()==0){
                return 0;
            }
        }
        return 1;
    }
    
    /**
     * @brief vector to store all ILayer regaring the inspection layers created during slicing
     */
    vector<ILayer*> inspectionLayers;
    
    /**
     * @brief store the information of the inspected object.
     * @details this information is used in the functions that check for collisions
     */
    stl_file *stl_local;
    
    /**
     * @brief store the information of the inpscetion head/camera.
     * @details this information is used in the functions that check for collisions
     */
    stl_file *complete_head;
};

//------------------------------------------------------------------------------
//-------------- COMMON METHODS ------------------------------------------------
//------------------------------------------------------------------------------
/**
 * @brief Print approach statistics such total time, number of valid IPoint, number of created links
 * @param[in] layers created {@link ILayers} with all information regarding IPoint and generated IPath
 * @param[in] ms total time necessary to execute the method
 * @param[in] subPath integer defining if the method used sub-paths or not. 1 print sub-path information, 0 otherwise
 */
void printStatistics(ILayers * layers, int ms, int subPath);

/**
 * @brief Method that starts the inspection process
 * @details Slice for inspection a given mesh and the filename of the inspected object. Filename is also used to get the name of the object in order to create the filename of the export file
 * @param[in] m Mesh with the inspected object
 * @param[in] filename char* with the name (and directory) of the inspected object
 */
void Slice_inspection(Mesh * m, char *filename);

/**
 * @brief convert the initial created sub-paths to multiple considering the \a alpha value
 * @param layers ILayers with information regaring IPoints (location and normal vector)
 * @param alpha maximum number of points that made a sub-path
 */
void divideISubPath(ILayers* layers, int alpha);

/**
 * @brief Creates a OpenSCad file with the slices, layed out in 3D form. It also includes the inspection points and their corresponding normal vectors.
 * @image html exampleSCad.PNG "Figure: Example of an exported object to .scad format" width=50%
 * @image latex exampleSCad.PNG "Example of an exported object to .scad format" width=10cm
 * @param[in] stl_object stl_file with the object to be inspected
 * @param[in] filename char* specifying the filename of the scdad file to be created
 * @return int value identifying the state of the file creation
 * @retval 0 file sucessfully created
 * @retval 1 error creating/opening file
 */
int exportSLayersOpenSCadinspection(stl_file* stl_object, char *filename);

/**
 * @brief Create an output file enabling to set the storing directory and the extension of the new file
 * @param[in] stl_file char* with the name of the file containing the object in the stl format
 * @param[in] directory char* setting the location to store the new file
 * @param[in] extension char* with the extension of the file to be created
 * @return char* with the correctly formated name file (and directory) to be created 
 */
char* createDirectoryFile(char* stl_file, char* directory, char* extension);

/**
 * @brief Create the inspection layers (slicing process) from the object stl_file
 * @details Inspection layers are created and each IPoint is checked for collsions in order to set the Polygon and corresponding subpath(s)
 * @param[in] stl_object stl_file of the inspected object
 * @param[in] stl_inspectionHead stl_file of the inspection camera/head
 * @return ILayers object with all layers information
 * @note this method considers that the graph is not oriented
 */
ILayers createInspectionLayers(stl_file *stl_object, stl_file *stl_inspectionHead);

/**
 * @brief Create the inspection layers (slicing process) from the object stl_file
 * @details Inspection layers are created and each IPoint is checked for collsions in order to set the Polygon and corresponding subpath(s)
 * @param[in] stl_object stl_file of the inspected object
 * @param[in] stl_inspectionHead stl_file of the inspection camera/head
 * @return ILayers object with all layers information
 * @note this method considers that the graph is oriented
 */
ILayers createInspectionLayersDirected(stl_file *stl_object, stl_file *stl_inspectionHead);

/**
@brief Adjust the position of the inspection head considering the normal vector of the inspection point 
@details 
The .stl object represents the entire inspection head where the thermographic camera is placed at the point (0,0,0) to ensure that the rotations and translations are properly performed. 
According to the value of the normal vector that is associated to the inspection point, two different rotation angles are calculated. 
The B rotation angle is performed along the y-axis and the C rotation angle is performed along the z-axis. 
In addition to these rotations a translation is also executed according to the coordinate values of the inspection point. 
The inspection head is correctly placed according to the inspection point under analyzis.
@image html camera_adjust.PNG "Figure: Inspection head adjustment example" width=1000px
@image latex camera_adjust.PNG "Inspection head adjustment example" width=10cm
@param[in] stl_inspectionHead stl_file* with the information of the inspection head object 
@param[in] inspectionPoint IPoint with the inspection point information
@param[out] stl_adjusted_head stl_file* with the information of the adjusted inspection head object
@note The stl_file structure is included in the admesh library
\code{.h}
typedef struct {
        FILE          *fp;
        stl_facet     *facet_start;
        stl_edge      *edge_start;
        stl_hash_edge **heads;
        stl_hash_edge *tail;
        int           M;
        stl_neighbors *neighbors_start;
        v_indices_struct *v_indices;
        stl_vertex    *v_shared;
        stl_stats     stats;
        char          error;
} stl_file;
\endcode
@return integer value representing the status of the inspection head adjustment
@retval 0 inspection head correctly adjusted
@retval -1 error loading the facets of the stl object
 */
int adjust_inspection_head(stl_file *stl_inspectionHead, IPoint* inspectionPoint, stl_file *stl_adjusted_head);

/**
 * @brief create ISubPath considering the previously created vector of LineSegment
 * @details in this method the IPoint is checked for collision to verify if it is a valid point
 * @param[in,out] polygon IPoly reference with created information
 * @param[in] polygonIndex index of the current polygon in the ILayer vector
 * @param[in] segments vector<LineSegments> with the segments that made the polygon
 * @param[in] stl_inspectionHead stl_file with the inspection head
 * @param[in] stl_object stl_file with the object
 */
void createISubPaths(IPoly * polygon, int polygonIndex, vector<LineSegment> segments, stl_file * stl_inspectionHead, stl_file * stl_object);

/**
 * @brief Create the graph (adjacent list) with valid connection between IPoint
 * @param[in,out] layer to create the graph
 * @param[in] segments vector<LineSegments> with the segments that made the polygon
 * @param[in] stl_inspectionHead stl_file with the inspection head
 * @param[in] stl_object stl_file with the object
 * @note this method consider the directed graph
 */
void createAdjacentList(ILayer * layer, vector < LineSegment > segments, stl_file * stl_inspectionHead, stl_file * stl_object);

/**
 * @brief Create the graph (adjacent list) with valid connection between IPoint
 * @param[in,out] layer to create the graph
 * @param[in] segments vector<LineSegments> with the segments that made the polygon
 * @param[in] stl_inspectionHead stl_file with the inspection head
 * @param[in] stl_object stl_file with the object
 * @note this method consider the bidirectional graph
 */
void createAdjacentListDirected(ILayer * layer, vector < LineSegment > segments, stl_file * stl_inspectionHead, stl_file * stl_object);


//------------------------------------------------------------------------------
//-------------- COLLISION DETECTION -------------------------------------------
//------------------------------------------------------------------------------
/**
@brief Check if two stl objects collide. (check if a triangle of a stl intersect other triangle from a different stl)
@image html stlCollision.PNG "Figure: stl object collision example" width=1000px
@image latex stlCollision.PNG "stl object collision example" width=10cm
@details In a first phase this function determines if there is a collision between the bounding boxes of the two .stl objects.
If the bounding boxes do not collide then there is no possible collision. 
When there is a collision of bounding boxes, each triangle of each object is tested to see if there is a possible collision with another triangle of another object.
Once a collision is detected the function ends with collision status 1. This method is used and described by Moller in @cite Moller1997.
@image html TriangleIntersection.PNG "Figure: triangle intersection example adapted from [a]" width=600px
@image latex TriangleIntersection.PNG triangle intersection example adapted from [a]" width=10cm

@param[in] stl1 stl_file* with the stl information of object 1
@param[in] stl2 stl_file* with the stl information of object 2
@note The stl_file structure is included in the admesh library
\code{.h}
typedef struct {
        FILE          *fp;
        stl_facet     *facet_start;
        stl_edge      *edge_start;
        stl_hash_edge **heads;
        stl_hash_edge *tail;
        int           M;
        stl_neighbors *neighbors_start;
        v_indices_struct *v_indices;
        stl_vertex    *v_shared;
        stl_stats     stats;
        char          error;
} stl_file;
\endcode
@return integer value that indicates the collision status
@retval 1 collision detected
@retval 0 otherwise
 */
int checkSTLCollision(stl_file *stl1, stl_file *stl2);


int checkSTLCollision2(stl_file *stl1, stl_file *stl2);

/**
@brief Check if a line segment is overlapping another one
@details Verify in 1D axis if two line segments have any overlapping section
@image html overlapping1d.PNG "Figure: Segment overlapping example" width=700px
@image latex overlapping1d.PNG "Segment overlapping example" width=10cm
@param[in] seg1min float minimum value of the line segment 1
@param[in] seg1max float maximum value of the line segment 1
@param[in] seg2min float minimum value of the line segment 2
@param[in] seg2max float maximum value of the line segment 2
@return bool value indicating if there is an overlapping
@retval true the two line segments have an overlapping region
@retval false otherwise
 */
bool overlapping1D(float seg1min, float seg1max, float seg2min, float seg2max);

/**
@brief Check if a 3D box is overlapping another one
@details Verify in 3D axis if two boxes have any overlapping area (intersect each other)
@image html overlapping3d.PNG "Figure: 3D box overlapping example" width=900px
@image latex overlapping3d.PNG "3D box overlapping example" width=10cm
@param[in] box1min stl_vertex minimum value of the line segment 1
@param[in] box1max stl_vertex maximum value of the line segment 1
@param[in] box2min stl_vertex minimum value of the line segment 2
@param[in] box2max stl_vertex maximum value of the line segment 2
@note The stl_vertex structure is included in the admesh library 
\code{.h}
typedef struct {
        float x;
        float y;
        float z;
} stl_vertex;
\endcode
@return bool value indicating if there is an overlapping
@retval true the two line segments have an overlapping region
@retval false otherwise
 */
bool overlapping3D(stl_vertex box1min, stl_vertex box1max, stl_vertex box2min, stl_vertex box2max);

/**
 * @brief Given three colinear points p, q, r, the function checks if point q lies on line segment 'pr' 
 * @param[in] p v3 with information regarding point p
 * @param[in] q v3 with information regarding point q
 * @param[in] r v3 with information regarding point r
 * @return boolean value stating if the point q is on segment pr
 * @retval true point q is on segment pr
 * @retval false otherwise
 */
bool onSegment(v3 p, v3 q, v3 r);

/**
 * @brief Find orientation of ordered triplet (a, b, c).
 * @details  Orientation of an ordered triplet of points in the plane can be counterclockwise, clockwise, or colinear.
 * @details The following diagram shows different possible orientations of (a, b, c)
 * @image html orientation-of-3-order-points-1.PNG "Figure: Orientation of 3 ordered points" width=700px
 * @image latex orientation-of-3-order-points-1.PNG "Orientation of three ordered points" width=10cm
 * @param[in] p v3 with information regarding point p
 * @param[in] q v3 with information regarding point q
 * @param[in] r v3 with information regarding point r
 * @return int value stating the orientation of the triplet
 * @retval 0 p, q and r are colinear 
 * @retval 1 Clockwise
 * @retval 2 Counterclockwise
 * @see [More details at]: https://www.geeksforgeeks.org/orientation-3-ordered-points/
 */
int orientation(v3 p, v3 q, v3 r);

/**
 * @brief The main function that returns true if line segment (p1,q1) and (p2,q2) intersect. 
 * @details Two segments (p1,q1) and (p2,q2) intersect if and only if one of the following two conditions is verified
 * @details
 * General Case:
 * 
 *  – (p1, q1, p2) and (p1, q1, q2) have different orientations and
 * 
 *  – (p2, q2, p1) and (p2, q2, q1) have different orientations.
 * 
 * @image html linesegments1.PNG "Figure: Example of different types of intersections" width=700px
 * @image latex linesegments1.PNG "Example of different types of intersections" width=10cm
 * 
 * Special Case:
 * 
 *  – (p1, q1, p2), (p1, q1, q2), (p2, q2, p1), and (p2, q2, q1) are all collinear and
 * 
 *  – the x-projections of (p1, q1) and (p2, q2) intersect
 * 
 *  – the y-projections of (p1, q1) and (p2, q2) intersect.
 * 
 * @image html linesegments2.PNG "Figure: Example of special types of intersections" width=700px
 * @image latex linesegments2.PNG "Example of special types of intersections" width=10cm
 * 
 * @param p1 v3 with information regarding point p1
 * @param q1 v3 with information regarding point q1
 * @param p2 v3 with information regarding point p2
 * @param q2 v3 with information regarding point q2
 * @return boolean value stating if the two line segments do intersect
 * @retval true the two segments intersect
 * @retval false otherwise
 * 
 * @see [External source]:https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
 */
bool doIntersect(v3 p1, v3 q1, v3 p2, v3 q2);

/**
 * @brief check if there is a collision between two IPoint
 * @param[in] from initial point of the link
 * @param[in] to end point of the link
 * @param[in] segments vector<LineSegments> to check if the link crosses any segment
 * @param[in] stl_inspectionHead stl_file with the inspection head
 * @param[in] stl_object stl_file with the object
 * @return integer identifying if a collision occurs
 * @retval 1 there is a collision
 * @retval 0 otherwise
 * @note if no collision link is identified, the method automatically calls collisionPoint and return its value
 */
int collisionLink(IPoint* from, IPoint* to, vector<LineSegment> segments, stl_file * stl_inspectionHead, stl_file * stl_object);

/**
 * @brief verify if it is possible to place the inspection head at a given point without colliding with the object
 * @param[in] stl_inspectionHead stl_file with the inspection head
 * @param[in] point IPoint to verify if placing the inspection camera in that location will collide with the object
 * @param[in] stl_object stl_file with the object
 * @retval 1 there is a collision
 * @retval 0 otherwise
 */
int collisionPoint(stl_file * stl_inspectionHead, IPoint * point, stl_file * stl_object);

//------------------------------------------------------------------------------
//-------------- NNGH ----------------------------------------------------------
//------------------------------------------------------------------------------
/**
 * @brief Method to create path using the Nearest Neighbor Graph-orientation Heuristic
 * @details The Nearest Neighbor Graph-orientation Heuristic (k-NNGH ) is a graph-based
approach, since it uses the information of the inspection graph to avoid collisions
to the object. The k parameter represents the number of inspection points that
may be used to extend an inspection trajectory at each iteration. When k takes
the value 1, only the closest point in each node is considered in the trajectory.
For k > 1 the path is completed with k nearest inspection points, creating several
alternative trajectories in each step.
 * 
 * @image html NNGH-algorithm.JPG "Figure: Nearest Neighbor Heuristic with Graph orientation (k-NNGH)" width=500px
 * @image latex NNGH-algorithm.JPG "Nearest Neighbor Heuristic with Graph orientation (k-NNGH)" width=10cm
 * 
 * At each iteration each partial path is replicated k-1 times, creating k-1 new
partial trajectories. This procedure is repeated until all inspection points have
been visited. If a partial path visits all possible inspection points then it becomes 
a final inspection trajectory. The aim of this heuristic is to find a valid trajectory
of good quality to inspect an object. It should be noted that when k = N, where
N is the total number of inspection points, the algorithm finds all valid possible
solutions, becoming a complete combinatorial algorithm. Finally, the k-NNGH
determines the most attractive path of the list of possible paths, i.e., the path
that minimizes the total time of inspection.
 * 
 * @param[in,out] layers created layers with information regarind IPoint and the created graph. During this method the path is updated inside the corresponding class
 * @param[in] k integer defining how many new possible paths will be searched
 */
void k_NNGH(ILayers *layers, int k);

/**
 * @brief Method to create path using the Nearest Neighbor Graph-orientation Heuristic for a specific ILayer, see k_NNGH() for details
 * @param[in,out] layer ILayer to create the path
 * @param[in] startPoint IPoint setting the starting point of the path that will be generated
 * @param[in] numToInsert integer defining how many new possible paths will be searched
 * @return the generated IPath
 * @note when starting the NNGH path generation, the start point is set to NULL
 */
void k_NNGH_generatePath(ILayer * layer, IPoint * startPoint, int numToInsert);

//------------------------------------------------------------------------------
//-------------- NNSH ----------------------------------------------------------
//------------------------------------------------------------------------------
/**
 * @brief Method to create path using the Nearest Neighbor Sub-path Heuristic
 * 
 * The k parameter defines how many closest subpaths
will be taken into account when creating a final inspection path. When k
takes the value 1 the algorithm takes a local view of the problem and attempts
to complete the inspection path by looking for the nearest sub-path of a polygon
that has not yet been visited, as long as it does not represent a possible
collision. For k > 1 this approach has a more global view of the problem, not
only concerned with the closest path in each step, but also with the k possible
paths that has not been visited in each decision. 
 * 
 * @image html NNSH-algorithm.PNG "Figure: Nearest Neighbor Sub-path Heuristic (k-NNSH)" width=500px
 * @image latex NNSH-algorithm.PNG "Nearest Neighbor Sub-path Heuristic (k-NNSH )" width=10cm
 *  
 Although in algorithmic terms this difference in the value of the
parameter k may seem irrelevant, this may lead to different final inspections
paths, depending on the object to be inspected. In a first step each inspection
point is analyzed in order to verify whether this point will be a valid point. A
point is considered valid when there is no collision between the positioned inspection
head and the object to be inspected. When an invalid point is found
within a particular path, determined with the aid of a spline, then this point is
removed and the path is divided into smaller inspection paths. These sub-paths
are created individually respecting the polygon being inspected. The sub-paths
are then stored and recombined. A sub-path may be recombined with a sub-path
of the same polygon or different polygons as long as there is no collision in the
path that connects the sub-paths. While the k-NNSH with k = 1 finds the nearest
path to the next sub-path, the k-NNSH with k > 1 completes the current
path with k valid possibilities to other sub-paths, creating several alternative
paths in each iteration. These will be completed over future iterations. Finally,
the most promising path will be the one that inspects the object in the shortest
possible time, one of all alternatives that represent valid inspection paths. This
algorithm allows efficient recombination of sub-paths by avoiding node-to-node
linking to create a new inspection path. However, for particular cases, it may
 * 
 * @param[in,out] layers created layers with information regarind IPoint. During this method the path is updated inside the corresponding class
 * @param[in] numToInsert integer defining how many new possible paths will be searched
 */
void k_NNSH(ILayers *layers, int numToInsert);

/**
 * @brief Method to create path using the Nearest Neighbor Sub-path Heuristic for a specific ILayer, see k_NNSH() for details
 * @param[in,out] layer ILayer with information with IPoint, collisions and splines
 * @param[in] startPoint initial IPoint for the layer
 * @param[in] numToInsert numToInsert integer defining how many new possible paths will be searched
 * @return the generated IPath
 * @note when starting the NNGH path generation, the start point is set to NULL
 * @warning in this method the splines are considered to be bi-directional
 */
void k_NNSH_generatePath(ILayer * layer, IPoint * startPoint, int numToInsert);

/**
 * @brief Method to create path using the Nearest Neighbor Sub-path Heuristic for a specific ILayer, see k_NNSH() for details
 * @param[in,out] layer ILayer with information with IPoint, collisions and splines
 * @param[in] startPoint initial IPoint for the layer
 * @param[in] numToInsert numToInsert integer defining how many new possible paths will be searched
 * @return the generated IPath
 * @note when starting the NNGH path generation, the start point is set to NULL
 * @warning in this method the splines are considered to be directed
 */
void k_NNSH_generatePathDirected(ILayer * layer, IPoint* startPoint, int numToInsert);

//------------------------------------------------------------------------------
//-------------- NNH -----------------------------------------------------------
//------------------------------------------------------------------------------
/**
 * @brief Method to create the path using the Nearest Neighbor Heuristic (NNH)
 * @details The Nearest Neighbor Heuristic (NNH) is a greedy heuristic whose main objective
is to find a valid inspection path in a faster way. This heuristic is simple,
since it does not take into account a global view of the problem nor the existence
of possible collisions. At each iteration, the algorithm selects an initial
point to start the inspection path. After selecting a starting point, the algorithm
successively explores the closest inspection point until no further points may be
inspected. At the end of the procedure, the minimum cost path is select, since
the algorithm generates one path for each initial inspection point. The algorithm
that represents the NNH approach is presented in the next figure.
 * 
 * @image html NNH-algorithm.PNG "Figure: Nearest Neighbor Heuristic (NNH)" width=500px
 * @image latex NNH-algorithm.PNG "Nearest Neighbor Heuristic (NNH)" width=10cm
 * 
 * The NNH receives
a list of points to inspect and selects an initial point to start the procedure
that tries to find the closest point. The initial point is placed in the first position
of the list and then the nearest inspection point, in terms of time, is searched
throughout the list. When this point is found it changes its position with the
next point that was analyzed. When the list is finished and no further exchange
is possible, the closest path to that starting point has been found. The algorithm
continues to look for more paths, in this case, one for each starting point. This
algorithm ends when it is not possible to determine any new paths, selecting the
most attractive way to inspect the object, i.e., the fastest way to inspect the
object.
 * 
 * @param[in,out] layers created layers with information regarind IPoint. During this method the path is updated inside the corresponding class 
 */
void NNH(ILayers * layers);

/**
 * @brief Method to create path using the Nearest Neighbor Heuristic for a specific ILayer, see NNH() for details
 * @param[in,out] layer ILayer with information with IPoint, collisions and splines
 * @param[in] startPoint initial IPoint for the layer
 * @return the generated IPath
 * @note when starting the NNH path generation, the start point is set to NULL
 */
void NNH_generatePath(ILayer * layer, IPoint * startPoint);

//------------------------------------------------------------------------------
//-------------- HEURISTIC RELATED ---------------------------------------------
//------------------------------------------------------------------------------
/**
 * @brief Find the closest point at a different layer
 * @param[in] p IPoint with the origin point
 * @param[in] layer ILayer to find the best IPoint
 * @return closest IPoint in the current ILayer considering the origin IPoint
 */
IPoint * findNextLayerInitialPoint(IPoint * p, ILayer * layer);

//------------------------------------------------------------------------------
//-------------- MIP -----------------------------------------------------------
//------------------------------------------------------------------------------
/**
 * @brief Method to create path using the Mathematical Integer Programming Model
 * 
 * In order to guarantee the optimality of a solution it is necessary to use exact
models to find the best solution. One of the approaches that ensures accurate
searching for an optimal solution is a Mathematical Integer Programming (MIP)
model. The trajectory inspection problem is similar to the classic travelling salesman
problem concerning the mathematical modelling, since the main goal is to
visit all inspection points minimizing the total travel time. In this problem there
is a set N with n inspection points (ip), where \f$N = \{ip_1, ... , ip_n\}\f$, and a travel
graph \f$T_{ij}\f$ where \f$(i,j) \in \{1..n\}\f$, \f$T_{ij} = T_{ji}\f$ and \f$T_{ii} = 0\f$. The links contained in
the graph \f$T_{ij}\f$ do not allow any type of collision with the object being inspected,
either in the layer itself or in adjacent layers. The parameters and the decision
variables are presented below:
 * 
 * <b>Parameters: </b>
 * 
 * \f{align*}{
        x_{ij}  &=  \begin{cases}  1  \text{ if inspection head $h$ goes from inspection point $i$ to inspection point $j$}, \nonumber\\ \forall i\in N \;and\; \forall j\in N \backslash \{i\}, \nonumber\\ 
                         0  \text{ otherwise} \end{cases} \nonumber\\
    u_{i}  &= \text{support variable which allows the elimination of sub-tours}, \forall i\in N\nonumber
   \f}
 * 
The mathematical model presented in (0) - (5) is based on the Miller-Tucker-Zemlin \cite Bektas2014, \cite Velednitsky2017 formulation. This model takes into account all the points that must be inspected according to the information provided by the graph. The main objective is to minimize the total time of an inspection trajectory (0).
 * 
<b> Minimize </b>
\f{align}{
\sum_{i=1 \in N}\sum_{j\neq i,j = 1 \in N} t_{ij}x_{ij} & \;\;\; (0)
\f}

 * <b>Subject to:</b>
\f{align}{
& \sum_{i=1,i \neq j \in N} x_{ij}=1,\quad \forall i \in N, & (1)\\
& \sum_{j=1, j \neq i \in N} x_{ij}=1,\quad \forall j \in N,& (2)\\
& u_{i} - u_{j} + Nx_{ij} \leq N-1,\quad \forall i \in N \;\text{,}\; \forall j \in N \;\text{and}\;  2 < i\neq j \leq N ,& (3)\\
& 0 \leq u_{i} \leq N-1, \quad \forall i \in N \;\text{and}\;  2 < i \leq N & (4)\\
& 0 \leq x_{ij} \leq 1, \quad \forall i \in N \;\text{and}\; \forall j \in N & (5)
\f}

The set containing  constraints (1) ensures that each inspection point is visited only once from a previous inspection point, i.e., is preceded by exactly one inspection point.  The second set of constraints determines that from one inspection point there is only one departure to another inspection point (1). 
 * The sets of constraints (3) and (4) guarantee that there is only one inspection trajectory that starts and ends at the same inspection point and that visits all points exactly once. The last set of constraints (5) defines the bounds of variables \f$x_{ij}\f$.

It is important to note that for the first inspection layer there is no defined starting point, so a random starting point is used. The connection with the highest cost of the optimal inspection path is subsequently removed. This decision, which takes into consideration the result of the MIP, avoids running the model for each possible initial point of inspection, since the model returns a closed path that begins and ends at the same inspection point. 
The following layers already have a defined starting point since they start the trajectory at the closest point to the end of the previous path.

 * 
 * @param[in,out] layers created layers with information regarind IPoint. During this method the path is updated inside the corresponding class
 */
void MIP(ILayers * layers);

/**
 * @brief Method to create path using the Mathematical Integer Programming Model for a specific ILayer, see MIP() for details
 * @param[in,out] layer ILayer with information with IPoint, collisions and splines
 * @param[in] startPoint IPoint with the starting point
 * @note when starting the MIP path generation, the start point is set to NULL
 */
void MIP_generatePath(ILayer * layer, IPoint * startPoint);

//------------------------------------------------------------------------------
//-------------- COMBF ---------------------------------------------------------
//------------------------------------------------------------------------------
/**
 * @brief Method to create path using the combinatorial formulation
 * 
 * 
 * @details The combinatorial formulation (CombF) allows to generate all possible inspection paths in order to find the inspection path with minimum cost. Given a set \f$N\f$ with \f$n\f$ inspection points (ip), where \f$N = \{ip_{1}, \ldots, ip_{n}\}\f$, 
 and a travel time matrix \f$T_{ij}\f$ where \f$(i,j) \in \{1...n\}\f$, \f$T_{ij} = T_{ji}\f$ 
  and \f$T_{ii} = 0\f$, the main goal is to determine the permutation \f$\pi \in P_{n} = \{p:\{1,\ldots,n\} \rightarrow \{1,\ldots,n\}\}\f$ 
  where the objective function \f$f: P_{n} \rightarrow \mathbb{R}\f$ 
   
\f{align*}{
	f(\pi) &= \sum_{i=1}^{n-1} T_{\pi(i), \pi(i+1)} + T_{\pi(n), \pi(1)}
\f}

 * 
has its minimum value. The main drawback of this formulation is its dependency with the number of inspection points, since the solution space exponentially increases. Due to the particular features of the inspection problem it is possible to decrease the computational effort by reducing the number of possible connections between inspection points. This is made by introducing collision constraints that detect the collision between the inspection head and the polygons being inspected. Despite this reduction, the high computational complexity depends on the type of object being inspected and the inspection features used. 
 * 
 * @image html CombF-algorithm.PNG "Figure: Combinatorial Formulation Algorithm (CombF)" width=500px
 * @image latex CombF-algorithm.PNG "Combinatorial Formulation Algorithm (CombF)" width=10cm
 
Previous Algorithm is a combinatorial algorithm that allows the creation of all possible inspection paths according to the inspection network and to the points to be inspected. 
 * The graph used by the algorithm is generated according to valid inspection points and collision constraints. There are two possible collision types, the first is when the inspection head collides with a polygon that may be in the layer itself or in one of the adjacent layers. 
 * The second is when part of the inspection network traverses a polygon. The combinatorial algorithm uses a list of partial paths and tries to complete the path while there are inspection points to visit. At each iteration of the algorithm a partial path is replicated and completed as many times as the nodes that have not yet been visited and are reachable. In this way, at each inspection point all possible combinations are created for the next inspection node. %These partial paths are temporarily saved. 
When a partial path has already visited all possible inspection points, it becomes a final inspection path. This procedure is repeated while it is possible to find partial inspection paths.
 * 
 * 
 * @param[in,out] layers created layers with information regarind IPoint. During this method the path is updated inside the corresponding class
 */
void CombF(ILayers *layers);

/**
 * @brief Method to create path using the combinatorial formulation for a specific ILayer, see CombF() for details
 * @param[in] layer ILayer with information with IPoint, collisions and splines
 * @param[in] startPoint IPoint with information of the starting point for the path of the current layer
 * @note when starting the CombF path generation, the start point is set to NULL
 * @return vector with only one generated IPath
 */
vector<IPath*> CombF_generatePath(ILayer * layer, IPoint * startPoint);

//------------------------------------------------------------------------------
//-------------- GCODE ---------------------------------------------------------
//------------------------------------------------------------------------------
/**
 * @brief Create an ipoint in the middle of two defined ipoints 
 * @param[in] from IPoint defining the origin point
 * @param[in] to IPoint defining the ending point
 * @param[out] middle the intermidiate created IPoint
 */
void createMiddleIPoint(IPoint * from, IPoint * to, IPoint * middle);

/**
@brief Export preamble of inspection to GCode inspection file
@details Print a commented GCode instruction specifying the beginning of the inspection to the Gcode file.
@verbatim
; setup the inspection machine
@endverbatim
@param[in] f FILE* specifying the destination file
@param[in] line integer setting if the line number should be printed in the GCODE (-1 if false, 0 otherwise)
 */
void exportPreamble_inspection(FILE *f, int * line);

/**
@brief Export epilogue of inspection to GCode inspection file
@details Print a commented GCode instruction specifying the ending of the inspection to the Gcode file.
@verbatim
; end inspection machine
@endverbatim
@param[in] f FILE* specifying the destination file
@param[in] line integer setting if the line number should be printed in the GCODE (-1 if false, 0 otherwise)
 */
void exportEpilogue_inspection(FILE *f, int * line);

/**
@brief Export G0 command to GCode inspection file
@details Print a G0 GCode instruction stating the position to which the camera should move.
@verbatim
G00 F<speed> X<x> Y<y> Z<z> B<b> C<c>
@endverbatim
@param[in] f FILE* specifying the destination file
@param[in] speed float value setting the speed which the head will move
@param[in] point IPoint with all information regarding x,y,z position and head rotation
@param[in] line integer setting if the line number should be printed in the GCODE (-1 if false, 0 otherwise)
 */
void exportG0command_inspection(FILE *f, float speed, IPoint * point, int * line);

/**
@brief Export M42 activation command to GCode inspection file
@details Print a M42 activation GCode instruction to activate the inspection. This is performed after setting the position of the inspection head.
@verbatim
M42P6S255 ; activate inspection
@endverbatim
@param[in] f FILE* specifying the destination file
@param[in] line integer setting if the line number should be printed in the GCODE (-1 if false, 0 otherwise)
 */
void exportM42activateInspectionCommand_inspection(FILE *f, int * line);

/**
@brief Export M42 deactivation command to GCode inspection file
@details Print a M42 deactivation GCode instruction to deactivate the inspection. This is performed to end the current active inspection.
@verbatim
M42P6S0 ; deactivate inspection
@endverbatim
@param[in] f FILE* specifying the destination file
@param[in] line integer setting if the line number should be printed in the GCODE (-1 if false, 0 otherwise)
 */
void exportM42deactivateInspectionCommand_inspection(FILE *f, int * line);

/**
@brief Export M84 deactivation command to GCode inspection file
@details Print a M84 deactivation GCode instruction to deactivate the motors and wait for the end of the current instruction.
@verbatim
M84 ; deactivate motors and wait for the end of the current instruction
@endverbatim
@param[in] f FILE* specifying the destination file
@param[in] line integer setting if the line number should be printed in the GCODE (-1 if false, 0 otherwise)
 */
void exportM84deactivateMotorsCommand_inspection(FILE *f, int* line);

/**
@brief Export G1E1 command to GCode inspection file
@details Print a G1E1 GCode instruction with the angle along the z-axis to which the inspection head should rotate to.
@verbatim
G1E1 <C> ; MrotZ rotation in degrees
@endverbatim
@param[in] f FILE* specifying the destination file
@param[in] C float value with the z-axis rotation angle
@param[in] line integer setting if the line number should be printed in the GCODE (-1 if false, 0 otherwise)
 */
void exportG1E1_MrotZ_inspection(FILE *f, float C, int* line);

/**
@brief Export G1E2 command to GCode inspection file
@details Print a G1E2 GCode instruction with the angle along the y-axis to which the inspection head should rotate to.
@verbatim
G1E2 <B> ; MrotY rotation in degrees
@endverbatim
@param[in] f FILE* specifying the destination file
@param[in] B float value with the y-axis rotation angle
@param[in] line integer setting if the line number should be printed in the GCODE (-1 if false, 0 otherwise)
 */
void exportG1E2_MrotY_inspection(FILE *f, float B, int * line);

/**
@brief Create a GCODE file with inspection instructions.
@details The inspection file contains a set of instructions that enables the activation and deactivation of the motors, the activation and deactivation of the inspection, and the positioning of the inspection head.
The inspection head positioning is done considering its coordinates (translation to the inspection point) and the corresponding rotation in y-axis and z-axis.

Example of a set of GCODE instructions:
\code{.gcode}
G00 F200.00 X-114.58 Y55.44 Z260.67 B0.04 C-75.79
M84 ; deactivate motors and wait for the end of the current instruction
G1E1 -75.79 ; MrotZ rotation in degrees
M84 ; deactivate motors and wait for the end of the current instruction
G1E2 0.04 ; MrotY rotation in degrees
M84 ; deactivate motors and wait for the end of the current instruction
M42P6S255 ; activate inspection
M84 ; deactivate motors and wait for the end of the current instruction
M42P6S0 ; deactivate inspection
M84 ; deactivate motors and wait for the end of the current instruction
\endcode

@param[in] path IPath with a set of ordered IPoint and their corresponding information such as position and rotation angles
@param[in] filename char* setting the name of the destination file with ".gcode" extension
@param[in] line integer setting if the line number should be printed in the GCODE (-1 if false, 0 otherwise)
@return int value with the status of the file creation
@retval 0 the file was correctly created
@retval 1 otherwise
 */
int exportGCode5Axes_layer(IPath * path, char *filename, int * line);

/**
 * @brief Creates a scad file with the trajectory and object information and the GCODE file with all instructions
 * @details For the files creation this method calls the exportGCode5Axes_layer() for each individual layer
 * @param[in] layers ILayers with all information regarding IPath and IPoint
 * @param[in] filename char* with the name of the file to be created
 * @param[in] lineNumbering integer setting the printing of the line number in the generated file (0 print lines, -1 otherwise)
 * @return int value with the status of the file creation
 * @retval 0 the file was correctly created
 * @retval 1 otherwise
 */
int exportGCode5Axes(ILayers * layers, char * filename, int lineNumbering);

#endif /* INSPECTION_HPP */

