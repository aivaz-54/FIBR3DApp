/**
 * file Inspection.cpp
 * author Bruna Ramos
 * date Created on 31 de Outubro de 2018, 14:47
 * brief Methods and procedures related to the inspection of the object
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include "clipper_utils.hpp"
#include "triangleIntersection/NoDivTriTriIsect.c"
#include <lp_lib.h>
#include "Inspection.hpp"
#include <chrono>

inspection inspection_char; /**< @brief inspection structure to store inspection details */

extern int print_level;

void Slice_inspection(Mesh * m, char *filename) {
    setbuf(stdout, NULL);

    stl_file stl_local;

    // initialize local stl structure
    stl_initialize(&stl_local);

    // make a copy, so rotate can be undone
    stl_local.stats = m->stl.stats;
    // allocate memory for facets
    stl_local.facet_start = (stl_facet*) malloc(stl_local.stats.number_of_facets * sizeof (stl_facet));
    if (stl_local.facet_start == NULL)
        perror("Inspection Slice");
    // copy facets
    memcpy(stl_local.facet_start, m->stl.facet_start, m->stl.stats.number_of_facets * sizeof (stl_facet));

    // rotate along x and y axis
    stl_rotate_x(&stl_local, m->optimal_rotations[0].x);
    stl_rotate_y(&stl_local, m->optimal_rotations[0].y);

    //put the object into positive coordinates for the inspection machine
    float x_trans = -stl_local.stats.min.x + inspection_char.distance_to_part; // + 1 * inspection_char.camXSize;
    float y_trans = -stl_local.stats.min.y + inspection_char.distance_to_part; // + 1 * (inspection_char.camLeftYSize + inspection_char.camRightYSize);
    stl_translate_relative(&stl_local, x_trans, y_trans, -stl_local.stats.min.z);

    //printf("%f %f\n", x_trans, y_trans);

    //read the inspection head stl file
    Mesh complete_head;
    stl_open(&complete_head.stl, (char*) "PrinterParts\\inspection_head.stl");
    stl_exit_on_error(&complete_head.stl);

    clock_t timeLayers = clock();
    ILayers layers;
    if (inspection_char.graph_directed) {
        layers = createInspectionLayersDirected(&stl_local, &complete_head.stl);
    } else {
        layers = createInspectionLayers(&stl_local, &complete_head.stl);
    }
    timeLayers = clock() - timeLayers;
    int msLayers = double(timeLayers) / CLOCKS_PER_SEC * 1000;
    if(print_level)
        printf("Layers/Graph time: %d ms\n", msLayers);

    //    for (auto polygon : layers.inspectionLayers.at(1)->inspectionPolygons) {
    //        for (auto isp : polygon->subPaths) {
    //            for (auto po : isp->subPath) {
    //                for (auto pa : po->adjacentList) {
    //                    printf("from %d  to %d\n", po->index, pa->index);
    //                }
    //            }
    //        }
    //    }

    if (inspection_char.MIP == 1) {
        if(print_level){
            printf("//------------------------------------------------------------------------------\n");
            printf("//-------------- MIP -----------------------------------------------------------\n");
            printf("//------------------------------------------------------------------------------\n");
        }
        clock_t time = clock();
        MIP(&layers);
        time = clock() - time;
        int ms = double(time) / CLOCKS_PER_SEC * 1000;
        if(print_level)
            printf("MIP time: %d ms\n", ms);
        //printf("Total time: %d ms\n", ms + graph_time + layer_time);
        if (layers.isPathComplete()) {
            char * fname = createDirectoryFile((char*) filename, (char*) "Gcode_files/", (char*) "_MIP.gcode");
            exportGCode5Axes(&layers, fname, inspection_char.printGcodeLines);
            delete fname;
        }
        if(print_level)
            printf("//------------------------------------------------------------------------------\n\n");

//        if (layers.isPathComplete()) {
//            printStatistics(&layers, ms, 0);
//        }
    }

    if (inspection_char.CombF == 1) {
        if(print_level){
            printf("//------------------------------------------------------------------------------\n");
            printf("//-------------- COMBF ---------------------------------------------------------\n");
            printf("//------------------------------------------------------------------------------\n");
        }
        clock_t time = clock();
        CombF(&layers);
        time = clock() - time;
        int ms = double(time) / CLOCKS_PER_SEC * 1000;
        if(print_level)
            printf("CombF time: %d ms\n", ms);
        //printf("Total time: %d ms\n", ms + graph_time + layer_time);
        if (layers.isPathComplete()) {
            char * fname = createDirectoryFile((char*) filename, (char*) "Gcode_files/", (char*) "_CombF.gcode");
            exportGCode5Axes(&layers, fname, inspection_char.printGcodeLines);
            delete fname;
        }
        if(print_level)
            printf("//------------------------------------------------------------------------------\n\n");
//        if (layers.isPathComplete()) {
//            printStatistics(&layers, ms, 0);
//        }
    }

    if (inspection_char.k_NNGH >= 1) {
        if(print_level){
            printf("//------------------------------------------------------------------------------\n");
            printf("//-------------- %d-NNGH --------------------------------------------------------\n", inspection_char.k_NNGH);
            printf("//------------------------------------------------------------------------------\n");
        }
        clock_t time = clock();
        k_NNGH(&layers, inspection_char.k_NNGH);
        time = clock() - time;
        int ms = double(time) / CLOCKS_PER_SEC * 1000;
        if(print_level)
            printf("%d-NNGH time: %d ms\n", inspection_char.k_NNGH, ms);
        //printf("Total time: %d ms\n", ms + graph_time + layer_time);
        if (layers.isPathComplete()) {
            char * fname = createDirectoryFile((char*) filename, (char*) "Gcode_files/", (char*) "_NNGH.gcode");
            exportGCode5Axes(&layers, fname, inspection_char.printGcodeLines);
            delete fname;
        }
        if(print_level)
            printf("//------------------------------------------------------------------------------\n\n");
//        if (layers.isPathComplete()) {
//            printStatistics(&layers, ms, 0);
//        }
    }

    if (inspection_char.NNH == 1) {
        if(print_level){
            printf("//------------------------------------------------------------------------------\n");
            printf("//-------------- NNH -----------------------------------------------------------\n");
            printf("//------------------------------------------------------------------------------\n");
        }
        clock_t time = clock();
        NNH(&layers);
        time = clock() - time;
        int ms = double(time) / CLOCKS_PER_SEC * 1000;
        if(print_level)
            printf("NNH time: %d ms\n", ms);
        //printf("Total time: %d ms\n", ms + layer_time);
        if (layers.isPathComplete()) {
            char * fname = createDirectoryFile((char*) filename, (char*) "Gcode_files/", (char*) "_NNH.gcode");
            exportGCode5Axes(&layers, fname, inspection_char.printGcodeLines);
            delete fname;
        }
        if(print_level)
            printf("//------------------------------------------------------------------------------\n\n");
//        if (layers.isPathComplete()) {
//            printStatistics(&layers, ms, 0);
//        }
    }

    if (inspection_char.k_NNSH > 0) {
        divideISubPath(&layers, inspection_char.alpha);
        if(print_level){
            printf("//------------------------------------------------------------------------------\n");
            printf("//-------------- %d-NNSH --------------------------------------------------------\n", inspection_char.k_NNSH);
            printf("//------------------------------------------------------------------------------\n");
        }
        time_t time = clock();
        k_NNSH(&layers, inspection_char.k_NNSH);
        time = clock() - time;
        int ms = double(time) / CLOCKS_PER_SEC * 1000;
        if(print_level)
            printf("%d-NNSH time: %d ms\n", inspection_char.k_NNSH, ms);
        //printf("Total time: %d ms\n", ms + layer_time);
        if (layers.isPathComplete()) {
            char * fname = createDirectoryFile((char*) filename, (char*) "Gcode_files/", (char*) "_NNSH.gcode");
            exportGCode5Axes(&layers, fname, inspection_char.printGcodeLines);
            delete fname;
        }
        if(print_level)
            printf("//------------------------------------------------------------------------------\n\n");
//        if (layers.isPathComplete()) {
//            printStatistics(&layers, ms, 1);
//        }
    }

    if (layers.isPathComplete()) {

        //show two points for positioning
        printf("Rotation x:%f  y:%f\n", m->optimal_rotations[0].x, m->optimal_rotations[0].y);

        IPoint * one = layers.inspectionLayers.front()->layerPath->path.front();
        printf("Point #1\n");
        v3 onev3 = one->point - one->normal * inspection_char.distance_to_part;
        one->point = onev3;
        //one->printSCad();
        printf("inspection machine: [%f,%f,%f]\n", one->point.x, one->point.y, one->point.z);

        onev3.x = onev3.x - x_trans;
        onev3.y = onev3.y - y_trans;

        //un-rotate in x
        double radian_angle = -(m->optimal_rotations[0].x / 180.0) * M_PI;
        double r = sqrt((onev3.y * onev3.y) + (onev3.z * onev3.z));
        double theta = atan2(onev3.z, onev3.y);
        onev3.y = r * cos(theta + radian_angle);
        onev3.z = r * sin(theta + radian_angle);

        //un-rotate in y
        radian_angle = -(m->optimal_rotations[0].y / 180.0) * M_PI;
        r = sqrt((onev3.z * onev3.z) + (onev3.x * onev3.x));
        theta = atan2(onev3.x, onev3.z);
        onev3.z = r * cos(theta + radian_angle);
        onev3.x = r * sin(theta + radian_angle);

        onev3 = onev3 + m->displacement;
        one->point = onev3;
        //one.printSCad();
        printf("original stl file: [%f,%f,%f]\n", one->point.x, one->point.y, one->point.z);

        IPoint * two = layers.inspectionLayers.front()->layerPath->path.at(layers.inspectionLayers.front()->layerPath->path.size() / 2);
        printf("Point #2\n");
        v3 twov3 = two->point - two->normal * inspection_char.distance_to_part;
        two->point = (twov3);
        //two->printSCad();
        printf("inspection machine: [%f,%f,%f]\n", two->point.x, two->point.y, two->point.z);
        //    printf("Is point in stl file: \n");

        twov3.x = twov3.x - x_trans;
        twov3.y = twov3.y - y_trans;

        //un-rotate in x
        radian_angle = -(m->optimal_rotations[0].x / 180.0) * M_PI;
        r = sqrt((twov3.y * twov3.y) + (twov3.z * twov3.z));
        theta = atan2(twov3.z, twov3.y);
        twov3.y = r * cos(theta + radian_angle);
        twov3.z = r * sin(theta + radian_angle);

        //un-rotate in y
        radian_angle = -(m->optimal_rotations[0].y / 180.0) * M_PI;
        r = sqrt((twov3.z * twov3.z) + (twov3.x * twov3.x));
        theta = atan2(twov3.x, twov3.z);
        twov3.z = r * cos(theta + radian_angle);
        twov3.x = r * sin(theta + radian_angle);

        twov3 = twov3 + m->displacement;
        two->point = (twov3);
        //two.printSCad();
        printf("original stl file: [%f,%f,%f]\n", two->point.x, two->point.y, two->point.z);
    }
    //free memory
    for (auto layer : layers.inspectionLayers) {
        for (auto polygon : layer->inspectionPolygons) {
            for (auto isp : polygon->subPaths) {
                delete isp;
            }
            for (auto point : polygon->inspectionPoints) {
                delete point;
            }
            delete polygon;
        }
        delete layer;
    }
}

void printStatistics(ILayers * layers, int ms, int subPath) {
    //layer => #IP => #ligações grafo => #sub paths (apenas para o NNSH) => eval path layer => eval de subida => eval total => tempo de execução global
    int l = 0;
    float totalEval = 0.0;
    for (auto layer : layers->inspectionLayers) {
        printf("%d ", ++l);

        int totalIP = 0;
        int totalLinks = 0;
        int totalSubPaths = 0;
        for (auto polygon : layer->inspectionPolygons) {
            for (auto isp : polygon->subPaths) {
                totalSubPaths++;
                for (auto point : isp->subPath) {
                    totalIP++;
                    totalLinks += point->adjacentList.size();
                }
            }
        }

        printf("%d %d ", totalIP, totalLinks);

        if (subPath) {
            printf("%d ", totalSubPaths); //apenas imprime se nnsh
        } else {
            printf("- ");
        }


        printf("%.2f ", layer->layerPath->eval);
        totalEval += layer->layerPath->eval;

        if (l != layers->inspectionLayers.size()) {
            float evalUp = layers->inspectionLayers.at(l - 1)->layerPath->path.back()->evalTo(layers->inspectionLayers.at(l)->layerPath->path.front());
            totalEval += evalUp;
            printf("%.2f ", evalUp);
        } else {
            printf("- ");
        }

        if (l == layers->inspectionLayers.size()) {
            printf("%.2f %d\n", totalEval, ms);
        } else {
            printf("- -\n");
        }

    }
}

int exportSLayersOpenSCadinspection(stl_file* stl_object, char *filename) {

    FILE *f2 = NULL;
    int i, j, k, ls, ps, ss, m, pos, segment;
    vector<SPolygon> spoly;
    vector<Spline> splines;
    float x_end, t;

    fopen_s(&f2, filename, "w");
    if (!f2) return 1;

    fprintf(f2, "N_polygon = 1; \nN_path = 3; \nN_inspectionP = 2; \ncolorPath = \"red\"; \ncolorPoly = \"blue\";\n\n");

    // Get spline layers with polygons
    vector<SLayer> SLayers = inspection_char.SLayers;

    fprintf(f2, "color(colorPoly){\n");

    // go for all layers
    ls = (int) SLayers.size();
    for (i = 0; i < ls; ++i) {
        // for each Layer
        // get the spolygon
        spoly = SLayers[i].getSPolygon();
        // get the number of Splines
        ps = (int) spoly.size();
        for (j = 0; j < ps; ++j) {
            if (spoly[j].getArea() > 0) {
                //printf("Outside Polygon\n");
                splines = spoly[j].getSplines();

                ss = (int) splines.size();
                //printf("Number of splines: %i\n", ss);

                //for each spline of a polygon
                for (k = 0; k < ss; ++k) {
                    //exterior points 
                    vector<v3> points = splines[k].get_y_data();

                    pos = (int) points.size();
                    for (m = 0; m < pos - 1; ++m) {
                        fprintf(f2, "\nline([%f,%f,%f],[%f,%f,%f],N_polygon);", points[m].x, points[m].y, points[m].z, points[m + 1].x, points[m + 1].y, points[m + 1].z);
                    }
                }
            }
        }
    }

    fprintf(f2, "}\n");

    for (i = 0; i < ls; ++i) {
        // for each Layer
        // get the spolygon
        spoly = SLayers[i].getSPolygon();
        fprintf(f2, "\n\n\n // Layer %d \n", i);
        // get the number of Splines
        ps = (int) spoly.size();
        for (j = 0; j < ps; ++j) {
            if (spoly[j].getArea() > 0) {
                splines = spoly[j].getSplines();

                ss = (int) splines.size();
                //printf("Number of splines: %i\n", ss);

                x_end = splines[ss - 1].get_back();

                //for each distance of inspection
                for (t = 0.0f; t < x_end; t += inspection_char.sampling_distance) {
                    for (k = 0; k < ss; ++k) {

                        //find the value t
                        if (t >= splines[k].get_front() && t <= splines[k].get_back()) {
                            // found it
                            vector<v3> Points = splines[k].get_y_data();
                            vector<v3> Normals = splines[k].get_normals();
                            pos = (int) Points.size();

                            //evaluate the spline value at the inspection position t
                            v3 point = spline_eval(splines[k], t, &segment);

                            //calculate approximate normal
                            v3 normal = Normals[segment - 1];

                            //calculate new inspection point
                            point = point + normal * inspection_char.distance_to_part;

                            fprintf(f2, "\ncolor(colorPath){line([%f,%f,%f],[%f,%f,%f],N_inspectionP);}", point.x, point.y, point.z, point.x + 30 * normal.x, point.y + 30 * normal.y, point.z + 30 * normal.z);
                            fprintf(f2, "\npoint([%f,%f,%f],N_inspectionP);//inpection point", point.x, point.y, point.z);

                            break; //when the element was found the search is stopped 
                        }
                    }
                }
            }
        }
    }

    fprintf(f2, "\n\nmodule line (a,b,r) { \n \t  v=b-a; \n \t translate(a){\n \t \t    rotate(a=atan2(v[1],v[0])){\n \t \t \t   rotate(v=[0,1,0],a=atan2(sqrt(v[0]*v[0]+v[1]*v[1]),v[2])){\n \t \t \t \t       cylinder(h=sqrt(v*v),r=r, $fn=30);\n \t \t \t   }\n \t \t   } \n  \t  }\n}\n\n");
    fprintf(f2, "module point(point1,radius){\n\t color(colorPath)\n\t translate(point1) sphere(1.5*radius, $fn=5);\n}\n");

    fclose(f2);
    return 0;
}

char* createDirectoryFile(char* stl_file, char* directory, char* extension) {
    char * filename = (char*) malloc(strlen(stl_file) + 1);
    strcpy(filename, stl_file);
    char * stl_object_name_temp, * stl_object_name;
    stl_object_name_temp = strtok(filename, "\\/");

    while (stl_object_name_temp != NULL) {
        stl_object_name_temp = strtok(NULL, "\\/");
        if (stl_object_name_temp != NULL) {
            stl_object_name = stl_object_name_temp;
        }
    }

    char * newFile;
    newFile = (char*) malloc(strlen(directory) + 1 + strlen(stl_object_name) + strlen(extension) + 1);
    strcpy(newFile, directory);
    strcat(newFile, stl_object_name);
    strcat(newFile, extension);
    return newFile;
}

ILayers createInspectionLayers(stl_file *stl_object, stl_file *stl_inspectionHead) {
    ILayers layers;

    int i, j, k, ls, ps, ss, segment, collision;
    vector<SPolygon> spoly;
    vector<Spline> splines;

    float x_end, t;
    float sliceSize = inspection_char.slicing;

    // Generate Slices
    vector< vector < LineSegment > > slicesWithLineSegments;
    triMeshSlicer(stl_object, slicesWithLineSegments, sliceSize);

    // Build layers with polygons
    vector<Layer> Layers;
    buildLayers(slicesWithLineSegments, Layers);

    // Build spline layers with polygons
    vector<SLayer> SLayers;
    buildSplineLayers(Layers, SLayers);

    inspection_char.SLayers = SLayers;

    // go for all layers
    ls = (int) SLayers.size();
    int indexP;
    for (i = 0; i < ls; ++i) {
        ILayer * layer = new ILayer();
        indexP = 0;

        // for each Layer get the spolygon
        spoly = SLayers[i].getSPolygon();

        // get the number of Splines
        ps = (int) spoly.size();
        int polygonIndex = -1;
        for (j = 0; j < ps; ++j) {

            if (spoly[j].getArea() > 0) {//outside polygon
                polygonIndex++;
                IPoly * polygon = new IPoly();

                splines = spoly[j].getSplines();
                ss = (int) splines.size();

                x_end = splines[ss - 1].get_back();

                //for each distance of inspection
                for (t = 0.0f; t < x_end; t += inspection_char.sampling_distance) {
                    for (k = 0; k < ss; ++k) {

                        //find the value t
                        if (t >= splines[k].get_front() && t <= splines[k].get_back()) {
                            // found it
                            vector<v3> Normals = splines[k].get_normals();

                            //evaluate the spline value at the inspection position t
                            v3 point = spline_eval(splines[k], t, &segment);

                            //calculate approximate normal
                            v3 normal = Normals[segment - 1];

                            //calculate new inspection point
                            point = point + normal * inspection_char.distance_to_part;

                            IPoint * inspectionPoint = new IPoint(point, normal, indexP++);

                            //check if the inspection point is valid (does not have collisions)
                            inspectionPoint->collide = collisionPoint(stl_inspectionHead, inspectionPoint, stl_object);

                            polygon->inspectionPoints.push_back(inspectionPoint);
                            break; //when the element was found the search is stopped 
                        }
                    }
                }

                createISubPaths(polygon, polygonIndex, slicesWithLineSegments[i], stl_inspectionHead, stl_object);
                layer->inspectionPolygons.push_back(polygon);
            }
        }
        createAdjacentList(layer, slicesWithLineSegments[i], stl_inspectionHead, stl_object);
        layers.inspectionLayers.push_back(layer);
    }

    layers.complete_head = stl_inspectionHead;
    layers.stl_local = stl_object;

    return layers;
}

ILayers createInspectionLayersDirected(stl_file *stl_object, stl_file *stl_inspectionHead) {
    ILayers layers;

    int i, j, k, ls, ps, ss, segment, collision;
    vector<SPolygon> spoly;
    vector<Spline> splines;

    float x_end, t;
    float sliceSize = inspection_char.slicing;

    // Generate Slices
    vector< vector < LineSegment > > slicesWithLineSegments;
    triMeshSlicer(stl_object, slicesWithLineSegments, sliceSize);

    // Build layers with polygons
    vector<Layer> Layers;
    buildLayers(slicesWithLineSegments, Layers);

    // Build spline layers with polygons
    vector<SLayer> SLayers;
    buildSplineLayers(Layers, SLayers);

    inspection_char.SLayers = SLayers;

    // go for all layers
    ls = (int) SLayers.size();
    int indexP;
    for (i = 0; i < ls; ++i) {
        ILayer * layer = new ILayer();
        indexP = 0;

        // for each Layer get the spolygon
        spoly = SLayers[i].getSPolygon();

        // get the number of Splines
        ps = (int) spoly.size();
        int polygonIndex = -1;
        for (j = 0; j < ps; ++j) {

            if (spoly[j].getArea() > 0) {//outside polygon
                polygonIndex++;
                IPoly * polygon = new IPoly();

                splines = spoly[j].getSplines();
                ss = (int) splines.size();

                x_end = splines[ss - 1].get_back();

                //for each distance of inspection
                for (t = 0.0f; t < x_end; t += inspection_char.sampling_distance) {
                    for (k = 0; k < ss; ++k) {

                        //find the value t
                        if (t >= splines[k].get_front() && t <= splines[k].get_back()) {
                            // found it
                            vector<v3> Normals = splines[k].get_normals();

                            //evaluate the spline value at the inspection position t
                            v3 point = spline_eval(splines[k], t, &segment);

                            //calculate approximate normal
                            v3 normal = Normals[segment - 1];

                            //calculate new inspection point
                            point = point + normal * inspection_char.distance_to_part;

                            IPoint * inspectionPoint = new IPoint(point, normal, indexP++);

                            //check if the inspection point is valid (does not have collisions)
                            inspectionPoint->collide = collisionPoint(stl_inspectionHead, inspectionPoint, stl_object);

                            if ((i + 1) % 2 == 0) {
                                polygon->push_first(inspectionPoint);
                            } else {
                                polygon->inspectionPoints.push_back(inspectionPoint);
                            }
                            break; //when the element was found the search is stopped 
                        }
                    }
                }

                createISubPaths(polygon, polygonIndex, slicesWithLineSegments[i], stl_inspectionHead, stl_object);
                layer->inspectionPolygons.push_back(polygon);
            }
        }
        createAdjacentListDirected(layer, slicesWithLineSegments[i], stl_inspectionHead, stl_object);
        layers.inspectionLayers.push_back(layer);
    }

    layers.complete_head = stl_inspectionHead;
    layers.stl_local = stl_object;

    return layers;
}

void divideISubPath(ILayers* layers, int alpha) {
    for (auto layer : layers->inspectionLayers) {
        for (auto polygon : layer->inspectionPolygons) {
            vector<ISubPath*> currentSubPath = polygon->subPaths;

            polygon->subPaths.clear();
            ISubPath * newISP = new ISubPath();
            for (auto subpath : currentSubPath) {

                vector<IPoint*> currentSubPath = subpath->subPath;
                //vector<IPoint*> newSubPath;// = new vector<>();
                for (int i = 0; i < currentSubPath.size(); i++) {
                    IPoint * point = currentSubPath.at(i);
                    point->subPath = polygon->subPaths.size();
                    point->position = newISP->subPath.size();

                    newISP->push_back(point);
                    if (newISP->subPath.size() == alpha) {
                        polygon->subPaths.push_back(newISP);
                        newISP = new ISubPath();
                    }
                }
                if (newISP->subPath.size() != 0) {
                    polygon->subPaths.push_back(newISP);
                } else {
                    delete newISP;
                }

                delete subpath;
            }
        }
    }
}

int adjust_inspection_head(stl_file *stl_inspectionHead, IPoint * inspectionPoint, stl_file *stl_adjusted_head) {
    // initialize local stl structure
    stl_initialize(stl_adjusted_head);

    // make a copy of stats, so rotate can be undone
    stl_adjusted_head->stats = stl_inspectionHead->stats;

    // allocate memory for facets
    stl_adjusted_head->facet_start = (stl_facet*) malloc(stl_adjusted_head->stats.number_of_facets * sizeof (stl_facet));
    if (stl_adjusted_head->facet_start == NULL) {
        perror("Inspection Slice");
        return -1;
    }

    // copy facets
    memcpy(stl_adjusted_head->facet_start, stl_inspectionHead->facet_start, stl_inspectionHead->stats.number_of_facets * sizeof (stl_facet));



    stl_rotate_y(stl_adjusted_head, inspectionPoint->b_angle);
    stl_rotate_z(stl_adjusted_head, inspectionPoint->c_angle);
    stl_translate_relative(stl_adjusted_head, inspectionPoint->point.x, inspectionPoint->point.y, inspectionPoint->point.z);
    //printf("[%f,%f,%f]\t[0,%f,%f]\n",inspectionPoint->point.x, inspectionPoint->point.y, inspectionPoint->point.z, inspectionPoint->b_angle,inspectionPoint->c_angle);

    return 0;
}

int checkSTLCollision(stl_file *stl1, stl_file *stl2) {
    int i, j;
    float v0[3], v1[3], v2[3], u0[3], u1[3], u2[3];

    //check if the bounding boxes collide
    if (!overlapping3D(stl1->stats.min, stl1->stats.max, stl2->stats.min, stl2->stats.max)) {
        return 0; //bounding boxes do not collide
    }

    //check if there is a real collision (despite the bounding boxes collision)
    for (i = 0; i < stl1->stats.number_of_facets; i++) {

        //convert the facet points to the Moller format
        v0[0] = stl1->facet_start[i].vertex[0].x;
        v0[1] = stl1->facet_start[i].vertex[0].y;
        v0[2] = stl1->facet_start[i].vertex[0].z;

        v1[0] = stl1->facet_start[i].vertex[1].x;
        v1[1] = stl1->facet_start[i].vertex[1].y;
        v1[2] = stl1->facet_start[i].vertex[1].z;

        v2[0] = stl1->facet_start[i].vertex[2].x;
        v2[1] = stl1->facet_start[i].vertex[2].y;
        v2[2] = stl1->facet_start[i].vertex[2].z;

        for (j = 0; j < stl2->stats.number_of_facets; j++) {
            //convert the facet points to the Moller format
            u0[0] = stl2->facet_start[j].vertex[0].x;
            u0[1] = stl2->facet_start[j].vertex[0].y;
            u0[2] = stl2->facet_start[j].vertex[0].z;

            u1[0] = stl2->facet_start[j].vertex[1].x;
            u1[1] = stl2->facet_start[j].vertex[1].y;
            u1[2] = stl2->facet_start[j].vertex[1].z;

            u2[0] = stl2->facet_start[j].vertex[2].x;
            u2[1] = stl2->facet_start[j].vertex[2].y;
            u2[2] = stl2->facet_start[j].vertex[2].z;

            if (NoDivTriTriIsect(v0, v1, v2, u0, u1, u2)) {
//                printf("line([%f,%f], [%f,%f] ,[%f,%f]);\n",v0[0],v1[0],v0[1],v1[1],v0[2],v1[2]);
//                printf("line([%f,%f], [%f,%f] ,[%f,%f]);\n",v1[0],v2[0],v1[1],v2[1],v1[2],v2[2]);
//                printf("line([%f,%f], [%f,%f] ,[%f,%f]);\n",v2[0],v0[0],v2[1],v0[1],v2[2],v0[2]);
//                printf("line([%f,%f], [%f,%f] ,[%f,%f]);\n",u0[0],u1[0],u0[1],u1[1],u0[2],u1[2]);
//                printf("line([%f,%f], [%f,%f] ,[%f,%f]);\n",u1[0],u2[0],u1[1],u2[1],u1[2],u2[2]);
//                printf("line([%f,%f], [%f,%f] ,[%f,%f]);\n",u2[0],u0[0],u2[1],u0[1],u2[2],u0[2]);
                return 1;
            }
        }
    }

    return 0;
}

int checkSTLCollision2(stl_file *stl1, stl_file *stl2) {
    int i, j;
    float v0[3], v1[3], v2[3], u0[3], u1[3], u2[3];

    //check if the bounding boxes collide
    if (!overlapping3D(stl1->stats.min, stl1->stats.max, stl2->stats.min, stl2->stats.max)) {
        return 0; //bounding boxes do not collide
    }

    //check if there is a real collision (despite the bounding boxes collision)
    for (i = 0; i < stl1->stats.number_of_facets; i++) {

        //convert the facet points to the Moller format
        v0[0] = stl1->facet_start[i].vertex[0].x;
        v0[1] = stl1->facet_start[i].vertex[0].y;
        v0[2] = stl1->facet_start[i].vertex[0].z;

        v1[0] = stl1->facet_start[i].vertex[1].x;
        v1[1] = stl1->facet_start[i].vertex[1].y;
        v1[2] = stl1->facet_start[i].vertex[1].z;

        v2[0] = stl1->facet_start[i].vertex[2].x;
        v2[1] = stl1->facet_start[i].vertex[2].y;
        v2[2] = stl1->facet_start[i].vertex[2].z;

        for (j = 0; j < stl2->stats.number_of_facets; j++) {
            //convert the facet points to the Moller format
            u0[0] = stl2->facet_start[j].vertex[0].x;
            u0[1] = stl2->facet_start[j].vertex[0].y;
            u0[2] = stl2->facet_start[j].vertex[0].z;

            u1[0] = stl2->facet_start[j].vertex[1].x;
            u1[1] = stl2->facet_start[j].vertex[1].y;
            u1[2] = stl2->facet_start[j].vertex[1].z;

            u2[0] = stl2->facet_start[j].vertex[2].x;
            u2[1] = stl2->facet_start[j].vertex[2].y;
            u2[2] = stl2->facet_start[j].vertex[2].z;

            if (NoDivTriTriIsect2(v0, v1, v2, u0, u1, u2)) {
                // MATLAB
                printf("line([%f,%f], [%f,%f] ,[%f,%f]);\n",v0[0],v1[0],v0[1],v1[1],v0[2],v1[2]);
                printf("line([%f,%f], [%f,%f] ,[%f,%f]);\n",v1[0],v2[0],v1[1],v2[1],v1[2],v2[2]);
                printf("line([%f,%f], [%f,%f] ,[%f,%f]);\n",v2[0],v0[0],v2[1],v0[1],v2[2],v0[2]);
                printf("line([%f,%f], [%f,%f] ,[%f,%f]);\n",u0[0],u1[0],u0[1],u1[1],u0[2],u1[2]);
                printf("line([%f,%f], [%f,%f] ,[%f,%f]);\n",u1[0],u2[0],u1[1],u2[1],u1[2],u2[2]);
                printf("line([%f,%f], [%f,%f] ,[%f,%f]);\n",u2[0],u0[0],u2[1],u0[1],u2[2],u0[2]);
//                printf("line([%f,%f,%f],[%f,%f,%f],N_polygon);\n",v0[0],v0[1],v0[2],v1[0],v1[1],v1[2]);
//                printf("line([%f,%f,%f],[%f,%f,%f],N_polygon);\n",v1[0],v1[1],v1[2],v2[0],v2[1],v2[2]);
//                printf("line([%f,%f,%f],[%f,%f,%f],N_polygon);\n",v2[0],v2[1],v2[2],v0[0],v0[1],v0[2]);
//                printf("line([%f,%f,%f],[%f,%f,%f],N_polygon);\n",u0[0],u0[1],u0[2],u1[0],u1[1],u1[2]);
//                printf("line([%f,%f,%f],[%f,%f,%f],N_polygon);\n",u1[0],u1[1],u1[2],u2[0],u2[1],u2[2]);
//                printf("line([%f,%f,%f],[%f,%f,%f],N_polygon);\n",u2[0],u2[1],u2[2],u0[0],u0[1],u0[2]);
                return 1;
            }
        }
    }

    return 0;
}

bool overlapping1D(float seg1min, float seg1max, float seg2min, float seg2max) {
    return ((seg1max >= seg2min) && (seg2max >= seg1min));
}

bool overlapping3D(stl_vertex box1min, stl_vertex box1max, stl_vertex box2min, stl_vertex box2max) {
    return (overlapping1D(box1min.x, box1max.x, box2min.x, box2max.x) &&
            overlapping1D(box1min.y, box1max.y, box2min.y, box2max.y) &&
            overlapping1D(box1min.z, box1max.z, box2min.z, box2max.z));
}

void createISubPaths(IPoly * polygon, int polygonIndex, vector<LineSegment> segments, stl_file * stl_inspectionHead, stl_file * stl_object) {
    //find the first collision point
    int cp = -1;
    for (int s = 0; s < polygon->inspectionPoints.size(); s++) {
        if (polygon->inspectionPoints.at(s)->collide == 1) {
            cp = s;
            break;
        }
    }

    if (cp == -1) {
        int size = polygon->inspectionPoints.size();
        int pos = 0;
        ISubPath * newISP = new ISubPath();
        int j = 0;
        for (j = 0; j < size; j++) {
            if (newISP->subPath.size() != 0) {

                if (collisionLink(polygon->inspectionPoints.at(pos), newISP->subPath.back(), segments, stl_inspectionHead, stl_object) == 0) {

                    polygon->inspectionPoints.at(pos)->position = j;
                    polygon->inspectionPoints.at(pos)->subPath = polygon->subPaths.size();
                    polygon->inspectionPoints.at(pos)->polygon = polygonIndex;

                    newISP->subPath.push_back(polygon->inspectionPoints.at(pos++));
                }
            } else {
                polygon->inspectionPoints.at(pos)->position = j;
                polygon->inspectionPoints.at(pos)->subPath = polygon->subPaths.size();
                polygon->inspectionPoints.at(pos)->polygon = polygonIndex;

                newISP->subPath.push_back(polygon->inspectionPoints.at(pos++));
            }
        }

        polygon->subPaths.push_back(newISP);
    } else {
        ISubPath * newSubPath = new ISubPath();
        int temp = 0;
        for (int s = cp + 1; s < cp + polygon->inspectionPoints.size(); s++) {
            if (polygon->inspectionPoints.at(s % polygon->inspectionPoints.size())->collide == 0) {

                polygon->inspectionPoints.at(s % polygon->inspectionPoints.size())->polygon = polygonIndex;
                polygon->inspectionPoints.at(s % polygon->inspectionPoints.size())->subPath = temp;
                polygon->inspectionPoints.at(s % polygon->inspectionPoints.size())->position = newSubPath->subPath.size();

                newSubPath->subPath.push_back(polygon->inspectionPoints.at(s % polygon->inspectionPoints.size()));
            } else {
                if (newSubPath->size() != 0) {
                    polygon->subPaths.push_back(newSubPath);
                    newSubPath = new ISubPath();
                }
            }
        }
        if (newSubPath->size() != 0) {
            polygon->subPaths.push_back(newSubPath);
        } else {
            free(newSubPath);
        }
    }
}

void createAdjacentList(ILayer * layer, vector < LineSegment > segments, stl_file * stl_inspectionHead, stl_file * stl_object) {
    for (auto polygon : layer->inspectionPolygons) {
        for (auto isp : polygon->subPaths) {

            for (int p = 0; p < isp->subPath.size() - 1; p++) {
                IPoint * from = isp->subPath.at(p);
                IPoint * to = isp->subPath.at(p + 1);

                if (!collisionLink(from, to, segments, stl_inspectionHead, stl_object)) {
                    from->adjacentList.push_back(to);
                    to->adjacentList.push_back(from);
                }
            }

            if (isp->subPath.size() > 2) {
                IPoint * from = isp->subPath.front();
                IPoint * to = isp->subPath.back();

                if (!collisionLink(from, to, segments, stl_inspectionHead, stl_object)) {
                    from->adjacentList.push_back(to);
                    to->adjacentList.push_back(from);
                }
            }

            for (auto point : isp->subPath) {
                for (auto polygon2 : layer->inspectionPolygons) {
                    for (auto isp2 : polygon2->subPaths) {
                        if (isp2 != isp) {
                            for (auto point2 : isp2->subPath) {
                                if (point != point2) {
                                    if (!collisionLink(point, point2, segments, stl_inspectionHead, stl_object)) {
                                        point->adjacentList.push_back(point2);
                                    }
                                }
                            }
                        }
                    }
                }

            }
        }
    }

    //    for (int i = 0; i < layer->inspectionPolygons.size(); i++) {
    //        IPoly * polygon = layer->inspectionPolygons.at(i);
    //        for (int j = 0; j < polygon->subPaths.size(); j++) {
    //            ISubPath * isp = polygon->subPaths.at(j);
    //            for (int w = 0; w < isp->subPath.size(); w++) {
    //                IPoint * point = isp->subPath.at(w);
    //                if (isp->subPath.size() >= 2) {
    //                    if (w == 0) {
    //                        if (isp->subPath.size() > 2) {
    //                            point->adjacentList.push_back(isp->subPath.at(w + 1));
    //                        }
    //                        if (!collisionLink(point, isp->subPath.back(), segments, stl_inspectionHead, stl_object)) {
    //                            point->adjacentList.push_back(isp->subPath.back());
    //                            isp->subPath.back()->adjacentList.push_back(point);
    //                        }
    //                    } else if ((isp->subPath.size() > 2) && (w == isp->subPath.size() - 1)) {
    //                        point->adjacentList.push_back(isp->subPath.at(w - 1));
    //                    } else if ((isp->subPath.size() > 2)) {
    //                        point->adjacentList.push_back(isp->subPath.at(w + 1));
    //                        point->adjacentList.push_back(isp->subPath.at(w - 1));
    //                    }
    //                }
    //                // check links to points in other polygons/subtpahs
    //                for (int x = 0; x < layer->inspectionPolygons.size(); x++) {
    //                    for (int y = 0; y < layer->inspectionPolygons.at(x)->subPaths.size(); y++) {
    //                        if (!(i == x && j == y)) {
    //                            //we dont want to link to the same subpath
    //                            ISubPath * ispTo = layer->inspectionPolygons.at(x)->subPaths.at(y);
    //                            for (int ww = 0; ww < ispTo->subPath.size(); ww++) {
    //                                IPoint * pointTo = ispTo->subPath.at(ww);
    //                                if (find(point->adjacentList.begin(), point->adjacentList.end(), pointTo) == point->adjacentList.end()) {
    //                                    if (!collisionLink(point, pointTo, segments, stl_inspectionHead, stl_object)) {
    //                                        point->adjacentList.push_back(pointTo);
    //                                        pointTo->adjacentList.push_back(point);
    //                                    }
    //                                }
    //                            }
    //                        }
    //                    }
    //                }
    //            }
    //        }
    //    }
}

void createAdjacentListDirected(ILayer * layer, vector < LineSegment > segments, stl_file * stl_inspectionHead, stl_file * stl_object) {
    //    for (int i = 0; i < layer->inspectionPolygons.size(); i++) {
    //        IPoly * polygon = layer->inspectionPolygons.at(i);
    //        for (int j = 0; j < polygon->subPaths.size(); j++) {
    //            ISubPath * isp = polygon->subPaths.at(j);
    //            for (int w = 0; w < isp->subPath.size(); w++) {
    //                IPoint * point = isp->subPath.at(w);
    //                if (isp->subPath.size() >= 2) {
    //                    if (w == isp->subPath.size() - 1) {
    //                        if (!collisionLink(point, isp->subPath.back(), segments, stl_inspectionHead, stl_object)) {
    //                            point->adjacentList.push_back(isp->subPath.front());
    //                        }
    //                    } else {
    //                        point->adjacentList.push_back(isp->subPath.at(w + 1));
    //                    }
    //                }
    //                // check links to points in other polygons/subtpahs
    //                for (int x = 0; x < layer->inspectionPolygons.size(); x++) {
    //                    for (int y = 0; y < layer->inspectionPolygons.at(x)->subPaths.size(); y++) {
    //                        if (!(i == x && j == y)) {
    //                            //we dont want to link to the same subpath
    //                            ISubPath * ispTo = layer->inspectionPolygons.at(x)->subPaths.at(y);
    //                            for (int ww = 0; ww < ispTo->subPath.size(); ww++) {
    //                                IPoint * pointTo = ispTo->subPath.at(ww);
    //                                if (find(point->adjacentList.begin(), point->adjacentList.end(), pointTo) == point->adjacentList.end()) {
    //                                    if (!collisionLink(point, pointTo, segments, stl_inspectionHead, stl_object)) {
    //                                        point->adjacentList.push_back(pointTo);
    //                                    }
    //                                }
    //                            }
    //                        }
    //                    }
    //                }
    //            }
    //        }
    //    }

    for (auto polygon : layer->inspectionPolygons) {
        for (auto isp : polygon->subPaths) {

            for (int p = 0; p < isp->subPath.size() - 1; p++) {
                IPoint * from = isp->subPath.at(p);
                IPoint * to = isp->subPath.at(p + 1);

                if (!collisionLink(from, to, segments, stl_inspectionHead, stl_object)) {
                    from->adjacentList.push_back(to);
                    //to->adjacentList.push_back(from);
                }
            }

            if (isp->subPath.size() > 2) {
                IPoint * from = isp->subPath.front();
                IPoint * to = isp->subPath.back();

                if (!collisionLink(from, to, segments, stl_inspectionHead, stl_object)) {
                    //from->adjacentList.push_back(to);
                    to->adjacentList.push_back(from);
                }
            }

            for (auto point : isp->subPath) {
                for (auto polygon2 : layer->inspectionPolygons) {
                    for (auto isp2 : polygon2->subPaths) {
                        if (isp2 != isp) {
                            for (auto point2 : isp2->subPath) {
                                if (point != point2) {
                                    if (!collisionLink(point, point2, segments, stl_inspectionHead, stl_object)) {
                                        point->adjacentList.push_back(point2);
                                    }
                                }
                            }
                        }
                    }
                }

            }
        }
    }

}

bool onSegment(v3 p, v3 q, v3 r) {
    if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) &&
            q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))
        return true;

    return false;
}

int orientation(v3 p, v3 q, v3 r) {
    int val = (q.y - p.y) * (r.x - q.x) -
            (q.x - p.x) * (r.y - q.y);

    if (val == 0) return 0; // colinear 

    return (val > 0) ? 1 : 2; // clock or counterclock wise 
}

bool doIntersect(v3 p1, v3 q1, v3 p2, v3 q2) {
    // Find the four orientations needed for general and 
    // special cases 
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    // General case 
    if (o1 != o2 && o3 != o4)
        return true;

    // Special Cases 
    // p1, q1 and p2 are colinear and p2 lies on segment p1q1 
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;

    // p1, q1 and q2 are colinear and q2 lies on segment p1q1 
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;

    // p2, q2 and p1 are colinear and p1 lies on segment p2q2 
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;

    // p2, q2 and q1 are colinear and q1 lies on segment p2q2 
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;

    return false; // Doesn't fall in any of the above cases 
}

int collisionLink(IPoint* from, IPoint* to, vector<LineSegment> segments, stl_file * stl_inspectionHead, stl_file * stl_object) {
    for (int i = 0; i < segments.size(); i++) {
        if (doIntersect(from->point, to->point, segments.at(i).v[0], segments.at(i).v[1])) {
            return 1;
        }
    }

    IPoint middle;
    createMiddleIPoint(from, to, &middle);

    return collisionPoint(stl_inspectionHead, &middle, stl_object);
}

int collisionPoint(stl_file * stl_inspectionHead, IPoint * point, stl_file * stl_object) {
    stl_file stl_adjusted_head;
    adjust_inspection_head(stl_inspectionHead, point, &stl_adjusted_head);

    int collide = 0;
    int collision = checkSTLCollision(&stl_adjusted_head, stl_object);
    if (collision) {
        collide = 1;
    }
    free(stl_adjusted_head.facet_start);
    return collide;
}

void k_NNGH(ILayers *layers, int numToInsert) {
    IPoint * start = NULL;
    //    vector<IPath*> possibilities;
    for (int i = 0; i < layers->inspectionLayers.size(); i++) {
        if (i > 0) {
            start = findNextLayerInitialPoint(layers->inspectionLayers.at(i - 1)->layerPath->path.back(), layers->inspectionLayers.at(i));

            //            int best = -1;
            //            float ev = numeric_limits<float>::max();
            //            for(int j = 0;j<possibilities.size();j++){
            //                IPoint * temp = findNextLayerInitialPoint(layers->inspectionLayers.at(i - 1)->layerPath->path.back(), layers->inspectionLayers.at(i));
            //                float e = temp->evalTo(layers->inspectionLayers.at(i - 1)->layerPath->path.back());
            //                if(e < ev){
            //                    if(best != -1){
            //                        delete possibilities.at(best);
            //                    }
            //                    start = temp;
            //                    ev = e;
            //                    best = j;
            //                }else{
            //                    delete possibilities.at(j);
            //                }
            //            }
            //            

        }
        /*possibilities =*/ k_NNGH_generatePath(layers->inspectionLayers.at(i), start, numToInsert);

        if (layers->inspectionLayers.at(i)->layerPath == NULL || layers->inspectionLayers.at(i)->layerPath->path.size() == 0) {
            //if(possibilities.size() == 0){
            printf("Cannot find valid path. Try again with new parameters...\n");
            break;
        }
    }
}

void k_NNGH_generatePath(ILayer* layer, IPoint* startPoint, int numToInsert) {
    vector<IPath*> finalPaths;
    vector<IPath*> partialPaths;

    int numberValidPoints = layer->getValidIPointSize();

    //create initial partial path(s) according to the possible set inital point
    if (startPoint != NULL) {
        IPoint * start = startPoint;
        //printf("%d", start->adjacentList.size()); getchar();

        for (auto next : start->adjacentList) {
            IPath * temp = new IPath();
            temp->path.push_back(start);
            temp->path.push_back(next);
            temp->eval = start->evalTo(next);

            //printf("[%d %d]EVAL AKI: %f\n",start->index,next->index,temp->eval);
            //getchar();


            //add the best k partial paths
            int pos = 0;
            for (int i = 0; i < partialPaths.size(); i++, pos++) {
                if (temp->eval < partialPaths.at(i)->eval) {

                    //printf("[%f]<[%f]",temp->eval,partialPaths.at(i)->eval);

                    break;
                }
            }

            if (pos < numToInsert) {
                partialPaths.insert(partialPaths.begin() + pos, temp);
                if (partialPaths.size() > numToInsert) {
                    delete partialPaths.back();
                    partialPaths.pop_back();
                }
            } else {
                delete temp;
            }
        }
    } else {
        for (auto polygon : layer->inspectionPolygons) {
            for (auto isp : polygon->subPaths) {
                for (auto start : isp->subPath) {

                    vector<IPath*> partialTemp;


                    for (auto next : start->adjacentList) {
                        IPath * temp = new IPath();
                        temp->path.push_back(start);
                        temp->path.push_back(next);
                        //temp->eval = start->evalTo(next);
                        temp->calculateEval();

                        //add the best k partial paths
                        int pos = 0;
                        for (int i = 0; i < partialTemp.size(); i++, pos++) {
                            if (temp->eval < partialTemp.at(i)->eval) {
                                break;
                            }
                        }

                        if (pos < numToInsert) {
                            partialTemp.insert(partialTemp.begin() + pos, temp);
                            if (partialTemp.size() > numToInsert) {
                                delete partialTemp.back();
                                partialTemp.pop_back();
                            }
                        } else {
                            delete temp;
                        }
                    }

                    for (auto ip : partialTemp) {
                        partialPaths.push_back(ip);
                    }
                }
            }
        }
    }

    while (partialPaths.size() != 0) {
        IPath * current = partialPaths.back();
        partialPaths.pop_back();

        vector<IPoint*> bestKPoints;
        vector<float> bestKEval;

        //current->print2();
        for (auto p : current->path.back()->adjacentList) {
            //printf("from %d   to %d  (of %d)\n",current->path.back()->index,p->index,current->path.back()->adjacentList.size());
            //get the k ordered best points that haven't been visited
            int valid = 1;
            for (int i = 0; i < current->path.size(); i++) {
                if (p->index == current->path.at(i)->index) {
                    valid = 0;
                    break;
                }
            }

            if (valid) {
                float eval = current->path.back()->evalTo(p);
                int pos = 0;
                for (int i = 0; i < bestKPoints.size(); i++, pos++) {
                    if (eval < bestKEval.at(i)) {
                        break;
                    }
                }

                if (pos < numToInsert) {
                    bestKPoints.insert(bestKPoints.begin() + pos, p);
                    bestKEval.insert(bestKEval.begin() + pos, eval);
                    if (bestKPoints.size() > numToInsert) {
                        bestKPoints.pop_back();
                        bestKEval.pop_back();
                    }
                }
            }
        }

        for (int i = 0; i < bestKPoints.size(); i++) {
            IPath * newPath = new IPath();
            newPath->path.insert(newPath->path.begin(), current->path.begin(), current->path.end());
            newPath->path.push_back(bestKPoints.at(i));
            //newPath->eval = current->eval + bestKEval.at(i);
            newPath->calculateEval();

            if (newPath->path.size() == numberValidPoints) {
                //new final path has been found. insert it in the correct position of the final vector
                int pos = 0;
                for (int j = 0; j < finalPaths.size(); j++, pos++) {
                    if (newPath->eval < finalPaths.at(j)->eval) {
                        break;
                    }
                }

                finalPaths.insert(finalPaths.begin() + pos, newPath);
            } else {
                partialPaths.push_back(newPath);
            }
        }

    }


    if (finalPaths.size() != 0) {
        layer->layerPath = finalPaths.front();
        //delete remaining paths
        for (int i = 1; i < finalPaths.size(); i++) {
            delete finalPaths.at(i);
        }
    } else {
        layer->layerPath = NULL;
    }
}

IPoint * findNextLayerInitialPoint(IPoint * p, ILayer * layer) {
    float best = std::numeric_limits<float>::max();
    IPoint * toPoint;
    for (int i = 0; i < layer->inspectionPolygons.size(); i++) {
        IPoly * ip = layer->inspectionPolygons.at(i);
        for (int j = 0; j < ip->subPaths.size(); j++) {
            ISubPath * isp = ip->subPaths.at(j);
            for (int point = 0; point < isp->subPath.size(); point++) {
                float d = p->evalTo(isp->subPath.at(point));
                if (d < best) {
                    best = d;
                    toPoint = isp->subPath.at(point);
                }
            }
        }
    }
    return toPoint;
}

void k_NNSH(ILayers *layers, int k) {
    IPoint * start = NULL;
    for (int i = 0; i < layers->inspectionLayers.size(); i++) {
        if (i > 0) {
            start = findNextLayerInitialPoint(layers->inspectionLayers.at(i - 1)->layerPath->path.back(), layers->inspectionLayers.at(i));
            //start->print();
            //check if the new starting point is in the middle of a subpath
            if (start->position != 0) {
                //printf("break\n");
                //layers->inspectionLayers.at(i)->print();
                //split current subpath in order to create two new ones where the startpoint is the initial point of a subpath
                vector<IPoint*> currentSubPath = layers->inspectionLayers.at(i)->inspectionPolygons.at(start->polygon)->subPaths.at(start->subPath)->subPath;

                vector<IPoint*> temp1;
                temp1.insert(temp1.begin(), currentSubPath.begin(), currentSubPath.begin() + start->position);
                layers->inspectionLayers.at(i)->inspectionPolygons.at(start->polygon)->subPaths.at(start->subPath)->subPath = temp1;

                vector<IPoint*> temp2;
                temp2.insert(temp2.begin(), currentSubPath.begin() + start->position, currentSubPath.end());
                for (int p = 0; p < temp2.size(); p++) {
                    temp2.at(p)->subPath = layers->inspectionLayers.at(i)->inspectionPolygons.at(start->polygon)->subPaths.size();
                    temp2.at(p)->position = p;
                }
                ISubPath * t = new ISubPath();
                t->subPath = temp2;
                layers->inspectionLayers.at(i)->inspectionPolygons.at(start->polygon)->subPaths.push_back(t);

                //printf("DEPOIS\n");
                //layers->inspectionLayers.at(i)->print();

                //start = layers->inspectionLayers.at(i)->inspectionPolygons.at(start->polygon)->subPaths.at(start->subPath)->subPath.front();
            }

            //            for (auto polygon : layers->inspectionLayers.at(i)->inspectionPolygons) {
            //                for (auto isp : polygon->subPaths) {
            //                    for (auto po : isp->subPath) {
            //                        for (auto pa : po->adjacentList) {
            //                            printf("from %d  to %d (%d)\n", po->index, pa->index,pa->position);
            //                        }
            //                    }
            //                }
            //            }
            //            printf("bla\n");

        }
        if (inspection_char.graph_directed) {
            k_NNSH_generatePathDirected(layers->inspectionLayers.at(i), start, k);
        } else {
            k_NNSH_generatePath(layers->inspectionLayers.at(i), start, k);
        }
        if (layers->inspectionLayers.at(i)->layerPath == NULL || layers->inspectionLayers.at(i)->layerPath->path.size() == 0) {
            printf("Cannot find valid path. Try again with new parameters...\n");
            break;
        }
    }
}

void k_NNSH_generatePath(ILayer * layer, IPoint* startPoint, int numToInsert) {
    vector<IPath*> finalPaths;
    vector<IPath*> partialPaths;

    int numberValidPoints = layer->getValidIPointSize();

    //create initial partial path(s) according to the possible set inital point
    if (startPoint != NULL) {
        IPoint * start = startPoint;
        IPath * initial = new IPath();
        initial->path.insert(initial->path.begin(), layer->inspectionPolygons.at(start->polygon)->subPaths.at(start->subPath)->subPath.begin(), layer->inspectionPolygons.at(start->polygon)->subPaths.at(start->subPath)->subPath.end());

        partialPaths.push_back(initial);
    } else {
        for (auto polygon : layer->inspectionPolygons) {
            for (auto isp : polygon->subPaths) {
                for (int w = 0; w < 2; w++) {
                    IPath * initial = new IPath();
                    vector<IPoint*> tt;
                    if (w == 0) {
                        tt.insert(tt.begin(), isp->subPath.begin(), isp->subPath.end());
                    } else {
                        tt.insert(tt.begin(), isp->subPath.crbegin(), isp->subPath.crend());
                    }
                    initial->path = tt;
                    vector<IPoint *> bestK;
                    vector<float> bestKe;
                    for (auto next : initial->path.back()->adjacentList) {
                        if (next->position == 0 || next->position == layer->inspectionPolygons.at(next->polygon)->subPaths.at(next->subPath)->subPath.size() - 1) {
                            int valid = 1;
                            for (auto pp : initial->path) {
                                if (pp->index == next->index) {
                                    valid = 0;
                                    break;
                                }
                            }

                            if (!valid) {
                                continue;
                            }

                            float eval = initial->path.back()->evalTo(next);

                            //add the best k partial paths
                            int pos = 0;
                            for (int i = 0; i < bestKe.size(); i++, pos++) {
                                if (eval < bestKe.at(i)) {
                                    break;
                                }
                            }
                            bestK.insert(bestK.begin() + pos, next);
                            bestKe.insert(bestKe.begin() + pos, eval);
                            if (bestKe.size() > numToInsert) {
                                bestKe.pop_back();
                                bestK.pop_back();
                            }
                        }
                    }

                    for (auto next : bestK) {
                        IPath * temp = new IPath();
                        temp->path.insert(temp->path.begin(), initial->path.begin(), initial->path.end());
                        if (next->position == 0) {
                            temp->path.insert(temp->path.end(), layer->inspectionPolygons.at(next->polygon)->subPaths.at(next->subPath)->subPath.begin(), layer->inspectionPolygons.at(next->polygon)->subPaths.at(next->subPath)->subPath.end());
                        } else {
                            temp->path.insert(temp->path.end(), layer->inspectionPolygons.at(next->polygon)->subPaths.at(next->subPath)->subPath.crbegin(), layer->inspectionPolygons.at(next->polygon)->subPaths.at(next->subPath)->subPath.crend());
                        }
                        temp->calculateEval();
                        partialPaths.push_back(temp);
                    }

                    delete initial;
                }
            }
        }
    }


    while (partialPaths.size() != 0) {
        IPath * current = partialPaths.back();
        partialPaths.pop_back();

        vector<IPoint*> bestKPoints;
        vector<float> bestKEval;
        for (auto p : current->path.back()->adjacentList) {
            if (p->position == 0 || p->position == layer->inspectionPolygons.at(p->polygon)->subPaths.at(p->subPath)->subPath.size() - 1) {

                //get the k ordered best points that haven't been visited
                int valid = 1;
                for (int i = 0; i < current->path.size(); i++) {
                    if (p->index == current->path.at(i)->index) {
                        valid = 0;
                        break;
                    }
                }

                if (!valid) {
                    continue;
                }
                float eval = current->path.back()->evalTo(p);
                int pos = 0;
                for (int i = 0; i < bestKPoints.size(); i++, pos++) {
                    if (eval < bestKEval.at(i)) {
                        break;
                    }
                }

                if (pos < numToInsert) {
                    bestKPoints.insert(bestKPoints.begin() + pos, p);
                    bestKEval.insert(bestKEval.begin() + pos, eval);
                    if (bestKPoints.size() > numToInsert) {
                        bestKPoints.pop_back();
                        bestKEval.pop_back();
                    }
                }
            }
        }

        for (int i = 0; i < bestKPoints.size(); i++) {
            IPoint * next = bestKPoints.at(i);
            IPath * newPath = new IPath();
            newPath->path.insert(newPath->path.begin(), current->path.begin(), current->path.end());
            if (next->position == 0) {
                newPath->path.insert(newPath->path.end(), layer->inspectionPolygons.at(next->polygon)->subPaths.at(next->subPath)->subPath.begin(), layer->inspectionPolygons.at(next->polygon)->subPaths.at(next->subPath)->subPath.end());
            } else {
                newPath->path.insert(newPath->path.end(), layer->inspectionPolygons.at(next->polygon)->subPaths.at(next->subPath)->subPath.crbegin(), layer->inspectionPolygons.at(next->polygon)->subPaths.at(next->subPath)->subPath.crend());
            }

            newPath->calculateEval();

            if (newPath->path.size() == numberValidPoints) {
                //new final path has been found. insert it in the correct position of the final vector
                int pos = 0;
                for (int j = 0; j < finalPaths.size(); j++, pos++) {
                    if (newPath->eval < finalPaths.at(j)->eval) {
                        break;
                    }
                }

                finalPaths.insert(finalPaths.begin() + pos, newPath);
            } else {
                partialPaths.push_back(newPath);
            }
        }
    }
    if (finalPaths.size() != 0) {
        layer->layerPath = finalPaths.front();

        //delete remaining paths
        for (int i = 1; i < finalPaths.size(); i++) {
            delete finalPaths.at(i);
        }
    } else {
        layer->layerPath = NULL;
    }
}

void k_NNSH_generatePathDirected(ILayer * layer, IPoint* startPoint, int numToInsert) {
    vector<IPath*> finalPaths;
    vector<IPath*> partialPaths;

    int numberValidPoints = layer->getValidIPointSize();

    //create initial partial path(s) according to the possible set inital point
    if (startPoint != NULL) {
        IPoint * start = startPoint;
        IPath * initial = new IPath();
        initial->path.insert(initial->path.begin(), layer->inspectionPolygons.at(start->polygon)->subPaths.at(start->subPath)->subPath.begin(), layer->inspectionPolygons.at(start->polygon)->subPaths.at(start->subPath)->subPath.end());

        partialPaths.push_back(initial);
    } else {
        for (auto polygon : layer->inspectionPolygons) {
            for (auto isp : polygon->subPaths) {

                IPath * initial = new IPath();
                vector<IPoint*> tt;

                tt.insert(tt.begin(), isp->subPath.begin(), isp->subPath.end());

                initial->path = tt;
                vector<IPoint *> bestK;
                vector<float> bestKe;
                for (auto next : initial->path.back()->adjacentList) {
                    if (next->position == 0) {
                        int valid = 1;
                        for (auto pp : initial->path) {
                            if (pp->index == next->index) {
                                valid = 0;
                                break;
                            }
                        }

                        if (!valid) {
                            continue;
                        }

                        float eval = initial->path.back()->evalTo(next);

                        //add the best k partial paths
                        int pos = 0;
                        for (int i = 0; i < bestKe.size(); i++, pos++) {
                            if (eval < bestKe.at(i)) {
                                break;
                            }
                        }
                        bestK.insert(bestK.begin() + pos, next);
                        bestKe.insert(bestKe.begin() + pos, eval);
                        if (bestKe.size() > numToInsert) {
                            bestKe.pop_back();
                            bestK.pop_back();
                        }
                    }
                }

                for (auto next : bestK) {
                    IPath * temp = new IPath();
                    temp->path.insert(temp->path.begin(), initial->path.begin(), initial->path.end());

                    temp->path.insert(temp->path.end(), layer->inspectionPolygons.at(next->polygon)->subPaths.at(next->subPath)->subPath.begin(), layer->inspectionPolygons.at(next->polygon)->subPaths.at(next->subPath)->subPath.end());

                    temp->calculateEval();
                    partialPaths.push_back(temp);
                }

                delete initial;

            }
        }
    }


    while (partialPaths.size() != 0) {
        IPath * current = partialPaths.back();
        partialPaths.pop_back();

        // current->print2();

        vector<IPoint*> bestKPoints;
        vector<float> bestKEval;
        for (auto p : current->path.back()->adjacentList) {
            //printf("from %d   to %d\n", current->path.back()->index, p->index);
            if (p->position == 0) {

                //get the k ordered best points that haven't been visited
                int valid = 1;
                for (int i = 0; i < current->path.size(); i++) {
                    if (p->index == current->path.at(i)->index) {
                        valid = 0;
                        break;
                    }
                }

                if (!valid) {
                    continue;
                }
                float eval = current->path.back()->evalTo(p);
                int pos = 0;
                for (int i = 0; i < bestKPoints.size(); i++, pos++) {
                    if (eval < bestKEval.at(i)) {
                        break;
                    }
                }

                if (pos < numToInsert) {
                    bestKPoints.insert(bestKPoints.begin() + pos, p);
                    bestKEval.insert(bestKEval.begin() + pos, eval);
                    if (bestKPoints.size() > numToInsert) {
                        bestKPoints.pop_back();
                        bestKEval.pop_back();
                    }
                }
            }
        }

        for (int i = 0; i < bestKPoints.size(); i++) {
            IPoint * next = bestKPoints.at(i);
            IPath * newPath = new IPath();
            newPath->path.insert(newPath->path.begin(), current->path.begin(), current->path.end());

            newPath->path.insert(newPath->path.end(), layer->inspectionPolygons.at(next->polygon)->subPaths.at(next->subPath)->subPath.begin(), layer->inspectionPolygons.at(next->polygon)->subPaths.at(next->subPath)->subPath.end());

            newPath->calculateEval();

            if (newPath->path.size() == numberValidPoints) {
                //new final path has been found. insert it in the correct position of the final vector
                int pos = 0;
                for (int j = 0; j < finalPaths.size(); j++, pos++) {
                    if (newPath->eval < finalPaths.at(j)->eval) {
                        break;
                    }
                }

                finalPaths.insert(finalPaths.begin() + pos, newPath);
            } else {
                partialPaths.push_back(newPath);
            }
        }
    }
    if (finalPaths.size() != 0) {
        layer->layerPath = finalPaths.front();
        //layer->layerPath->print2();
        //delete remaining paths
        for (int i = 1; i < finalPaths.size(); i++) {
            delete finalPaths.at(i);
        }
    } else {
        layer->layerPath = NULL;
    }
}

void NNH(ILayers * layers) {
    IPoint * start = NULL;
    for (int i = 0; i < layers->inspectionLayers.size(); i++) {
        //printf("layer %d\n",i);
        if (i > 0) {
            start = findNextLayerInitialPoint(layers->inspectionLayers.at(i - 1)->layerPath->path.back(), layers->inspectionLayers.at(i));
        }
        NNH_generatePath(layers->inspectionLayers.at(i), start);


        if (layers->inspectionLayers.at(i)->layerPath == NULL || layers->inspectionLayers.at(i)->layerPath->path.size() == 0) {
            printf("Cannot find valid path. Try again with new parameters...\n");
            break;
        }
    }
}

void NNH_generatePath(ILayer * layer, IPoint * startPoint) {
    IPath * bestPath = NULL;

    vector<IPoint*> generatedPathOriginal;
    for (auto polygon : layer->inspectionPolygons) {
        for (auto sp : polygon->subPaths) {
            for (auto point : sp->subPath) {
                generatedPathOriginal.push_back(point);
            }
        }
    }

    int size = 1;
    if (startPoint == NULL) {
        size = generatedPathOriginal.size();
    } else {
        //check if the startPoint is in the first position. if not change it to be the first
        if (generatedPathOriginal.front()->index != startPoint->index) {
            for (int i = 1; i < generatedPathOriginal.size(); i++) {
                if (generatedPathOriginal.at(i)->index == startPoint->index) {
                    iter_swap(generatedPathOriginal.begin(), generatedPathOriginal.begin() + i);
                    break;
                }
            }
        }
    }

    for (int k = 0; k < size; k++) {
        vector<IPoint*> generatedPath; //(generatedPathOriginal);
        generatedPath.insert(generatedPath.begin(), generatedPathOriginal.begin(), generatedPathOriginal.end());

        if (k != 0) {
            //only enters here if the startPoint == NULL
            //only chage if the point isnt already in the correct position...
            iter_swap(generatedPath.begin(), generatedPath.begin() + k);
        }

        int valid = 1;
        for (int i = 0; i < generatedPath.size() - 1; i++) {
            IPoint * to = NULL;
            float eval = std::numeric_limits<float>::max();
            //using the adjacent list get the closest point to point i
            for (auto toPoint : generatedPath.at(i)->adjacentList) {
                int cont = 1;
                //check if toPoint hasn't already be used...
                for (int x = 0; x <= i; x++) {
                    if (generatedPath.at(x)->index == toPoint->index) {
                        cont = 0;
                        break;
                    }
                }

                if (cont) {
                    float temp = generatedPath.at(i)->evalTo(toPoint);
                    if (temp < eval) {
                        eval = temp;
                        to = toPoint;
                    }
                }
            }

            if (to == NULL) {
                //this path is invalid...
                valid = 0;
                break;
            }

            //get the index of the next best point
            int index = -1;
            for (int j = i + 1; j < generatedPath.size(); j++) {
                if (generatedPath.at(j)->index == to->index) {
                    index = j;
                    break;
                }
            }

            if (index != (i + 1)) {//confirmar este if no caso de já se encontrar ordenado ele vai trocar por ele proprio sendo desnecessario
                iter_swap(generatedPath.begin() + index, generatedPath.begin() + (i + 1));
            }
        }

        if (valid == 1) {
            IPath * tempPath = new IPath();
            tempPath->path = generatedPath;
            tempPath->calculateEval();

            if (bestPath != NULL) {
                if (tempPath->eval < bestPath->eval) {
                    delete bestPath;
                    bestPath = tempPath;
                } else {
                    delete tempPath;
                }
            } else {
                bestPath = tempPath;
            }
        }
    }

    layer->layerPath = bestPath;
}

void MIP(ILayers * layers) {
    IPoint * start = NULL;
    for (int i = 0; i < layers->inspectionLayers.size(); i++) {
        IPoint * dummy = new IPoint();
        dummy->index = -10;

        if (i != 0) {
            start = findNextLayerInitialPoint(layers->inspectionLayers.at(i - 1)->layerPath->path.back(), layers->inspectionLayers.at(i));
            dummy->adjacentList.push_back(start);
        }


        for (auto polygon : layers->inspectionLayers.at(i)->inspectionPolygons) {
            for (auto isp : polygon->subPaths) {
                for (auto point : isp->subPath) {
                    point->adjacentList.push_back(dummy);
                    if (start == NULL) {
                        dummy->adjacentList.push_back(point);
                    }
                }
            }
        }

        ISubPath * isp = new ISubPath();
        isp->subPath.push_back(dummy);
        layers->inspectionLayers.at(i)->inspectionPolygons.back()->subPaths.push_back(isp);

        MIP_generatePath(layers->inspectionLayers.at(i), start);

        //remove the dummy points and links
        for (auto polygon : layers->inspectionLayers.at(i)->inspectionPolygons) {
            for (auto isp : polygon->subPaths) {
                for (auto point : isp->subPath) {
                    point->adjacentList.pop_back();
                }
            }
        }

        delete layers->inspectionLayers.at(i)->inspectionPolygons.back()->subPaths.back()->subPath.front();
        delete layers->inspectionLayers.at(i)->inspectionPolygons.back()->subPaths.back();
        layers->inspectionLayers.at(i)->inspectionPolygons.back()->subPaths.pop_back();

        if (layers->inspectionLayers.at(i)->layerPath == NULL || layers->inspectionLayers.at(i)->layerPath->path.size() == 0) {
            printf("Cannot find valid path. Try again with new parameters...\n");
            break;
        }

    }
}

void MIP_generatePath(ILayer * layer, IPoint * startPoint) {
    vector<IPoint*> validPoints;

    for (int i = 0; i < layer->inspectionPolygons.size(); i++) {
        IPoly * poly = layer->inspectionPolygons.at(i);
        for (int j = 0; j < poly->subPaths.size(); j++) {
            ISubPath * subPath = poly->subPaths.at(j);
            for (int k = 0; k < subPath->subPath.size(); k++) {
                IPoint * temp = subPath->subPath.at(k);
                validPoints.push_back(temp);
            }
        }
    }

    double matrix[validPoints.size()][validPoints.size()];

    for (int i = 0; i < validPoints.size(); i++) {
        for (int j = 0; j < validPoints.size() && j <= i; j++) {
            matrix[i][j] = -1;
            matrix[j][i] = -1;
        }
    }

    int start = -1;
    for (int i = 0; i < validPoints.size(); i++) {
        if (startPoint != NULL && validPoints.at(i)->index == startPoint->index) {
            start = i;
        }
        vector<IPoint*> temp = validPoints.at(i)->adjacentList;

        for (int x = 0; x < temp.size(); x++) {
            int j = -1;

            for (int z = 0; z < validPoints.size(); z++) {
                if (validPoints.at(z)->index == temp.at(x)->index) {
                    j = z;
                    break;
                }
            }

            matrix[i][j] = validPoints.at(i)->evalTo(validPoints.at(j));
        }
    }

    int colIndex[validPoints.size()][validPoints.size()];
    int colIndex2[validPoints.size()];
    int ncols = 0; //xAB,xBA...
    for (int i = 0; i < validPoints.size(); i++) {
        for (int j = 0; j < validPoints.size(); j++) {
            if (i == j || matrix[i][j] == -1) {
                colIndex[i][j] = -1;
            } else {
                colIndex[i][j] = ++ncols;
            }
        }
    }

    for (int i = 0; i < validPoints.size(); i++) {
        if (start == -1 && i == 0) {
            colIndex2[i] = -1;
        } else if (i != start) {
            colIndex2[i] = ++ncols;
        } else {
            colIndex2[i] = -1;
        }
    }

    // create LP
    lprec *lp;
    int status = 0;
    lp = make_lp(0, ncols);
    if (lp == NULL) {
        puts("Error creating lp...");
        status = 1;
    }

    int *colno = NULL;
    REAL *row = NULL;

    if (status == 0) {
        //let us name our variables. Not required, but can be useful for debugging 
        for (int i = 0; i < validPoints.size(); i++) {
            for (int j = 0; j < validPoints.size(); j++) {
                if (matrix[i][j] != -1) {
                    char name[100];

                    string from = "dummy";
                    if (validPoints.at(i)->index != -10) {
                        from = to_string(validPoints.at(i)->index);
                    }
                    char fromC[from.size() + 1];
                    strcpy(fromC, from.c_str());

                    string to = "dummy";
                    if (validPoints.at(j)->index != -10) {
                        to = to_string(validPoints.at(j)->index);
                    }
                    char toC[to.size() + 1];
                    strcpy(toC, to.c_str());

                    snprintf(name, 100, "x%s_%s", fromC, toC);
                    set_col_name(lp, colIndex[i][j], name);
                    set_bounds(lp, colIndex[i][j], 0.0, 1.0);
                    set_int(lp, colIndex[i][j], TRUE);
                }
            }
        }

        for (int i = 0; i < validPoints.size(); i++) {
            if (!(start == -1 && i == 0) && i != start) {
                char name[25];

                string from;
                if (validPoints.at(i)->index != -10) {
                    from = std::to_string(validPoints.at(i)->index);
                } else {
                    from = "dummy";
                }

                char fromC[from.size() + 1];
                strcpy(fromC, from.c_str());
                snprintf(name, 25, "u%s", fromC);
                set_col_name(lp, colIndex2[i], name);
                set_bounds(lp, colIndex2[i], 0.0, validPoints.size() - 1); //N-1
                set_int(lp, colIndex2[i], TRUE);
            }
        }

        colno = (int *) malloc(ncols * sizeof (*colno));
        row = (REAL *) malloc(ncols * sizeof (*row));

        if ((colno == NULL) || (row == NULL)) {
            status = 2;
        }

        set_add_rowmode(lp, TRUE); // makes building the model faster if it is done rows by row
    }

    if (status == 0) {
        //------------------------CONSTRAINTS 1 (OUT)--------------------------------
        //create space large enough for one row
        for (int c = 0; c < validPoints.size(); c++) {
            int nz = 0;
            for (int j = 0; j < validPoints.size(); j++) {
                if (c != j) {
                    if (matrix[c][j] != -1) {
                        colno[nz] = colIndex[c][j];
                        row[nz] = 1;
                        nz++;
                    }
                }
            }
            status = !add_constraintex(lp, nz, row, colno, EQ, 1);
        }
    }

    if (status == 0) {
        //------------------------CONSTRAINTS 2 (IN)--------------------------------
        //create space large enough for one row
        for (int c = 0; c < validPoints.size(); c++) {
            int nz = 0;
            for (int j = 0; j < validPoints.size(); j++) {
                if (c != j) {
                    if (matrix[j][c] != -1) {
                        colno[nz] = colIndex[j][c];
                        row[nz] = 1;
                        nz++;
                    }
                }
            }
            status = !add_constraintex(lp, nz, row, colno, EQ, 1);
        }
    }

    if (status == 0) {
        //------------------------CONSTRAINTS 3 (dummy variables)--------------------------------
        for (int c = 0; c < validPoints.size(); c++) {
            if (!(start == -1 && c == 0) && c != start) {
                for (int j = 0; j < validPoints.size(); j++) {
                    int nz = 0;
                    if ((!(start == -1 && j == 0) && j != start) && colIndex[c][j] != -1) {
                        colno[nz] = colIndex2[c];
                        row[nz] = 1;
                        nz++;

                        colno[nz] = colIndex2[j];
                        row[nz] = -1;
                        nz++;

                        colno[nz] = colIndex[c][j];
                        row[nz] = validPoints.size(); //N
                        nz++;

                        status = !add_constraintex(lp, nz, row, colno, LE, validPoints.size() - 1); //N-1
                    }
                }
            }
        }
    }

    if (status == 0) {
        //----------------- Objective function -----------------------
        set_add_rowmode(lp, FALSE); // rowmode should be turned off again when done building the model

        int nz = 0;
        for (int c = 0; c < validPoints.size(); c++) {
            for (int j = 0; j < validPoints.size(); j++) {
                if (matrix[c][j] != -1) {
                    colno[nz] = colIndex[c][j];
                    row[nz] = matrix[c][j];
                    nz++;
                }
            }
        }

        //set the objective in lpsolve 
        if (!set_obj_fnex(lp, nz, row, colno)) {
            status = 4;
        }
    }

    if (status == 0) {
        //set the object direction to maximize
        set_minim(lp);

        //just out of curioucity, now show the model in lp format on screen 
        // this only works if this is a console application. If not, use write_lp and a filename 
        //write_LP(lp, (char*)"exemplo.lp");
//        write_lp(lp, (char*) "LP_files/model.lp"); //TODO set the name to the object
        //printf("any key to continue...");getchar();

        // I only want to see important messages on screen while solving 
        set_verbose(lp, IMPORTANT);

        //NEUTRAL (0) 	Only some specific debug messages in de debug print routines are reported.
        //CRITICAL (1) 	Only critical messages are reported. Hard errors like instability, out of memory, ...
        //SEVERE (2) 	Only severe messages are reported. Errors.
        //IMPORTANT (3) Only important messages are reported. Warnings and Errors.
        //NORMAL (4) 	Normal messages are reported. This is the default.
        //DETAILED (5) 	Detailed messages are reported. Like model size, continuing B&B improvements, ...
        //FULL (6) 	All messages are reported. Useful for debugging purposes and small models.

        set_solutionlimit(lp, 10);

        // Now let lpsolve calculate a solution 
        status = solve(lp);
//        printf("Solution limit: %d  SOLUTION COUNT: %d\n", get_solutionlimit(lp), get_solutioncount(lp));
        if (status == OPTIMAL || status == SUBOPTIMAL) {
            status = 0;
        }
        //else {
        //    printf("No solution found!\n");
        //}
    }

    if (status == 0) {
        // a solution is calculated, now lets get some results

        // objective value
        // printf("Objective value: %f\n", get_objective(lp));

        // variable values 
        get_variables(lp, row);
        //         for (int j = 0; j < ncols; j++) {
        //         printf("%s: %f   (id:%d)\n", get_col_name(lp, j + 1), row[j], j + 1);
        //         }
        //        getchar();
        // we are done now

        vector<IPoint*> path_temp;

        if (startPoint == NULL) {
            start = 0;
        }

        path_temp.push_back(validPoints.at(start));

        int i = start;
        for (int ip = start; ip < validPoints.size() + start; ip++) {
            for (int j = 0; j < validPoints.size() /*&& !found*/; j++) {
                if (row[colIndex[i][j] - 1] == 1) {
                    IPoint * point_j = validPoints.at(j);

                    path_temp.push_back(point_j);
                    i = j;
                    break;
                }
            }
        }

        //remove the last point since the path is closed
        path_temp.pop_back();

        if (path_temp.size() == validPoints.size()) {



            vector<IPoint*> finalPath;

            int pos = 0;
            for (int i = 0; i < path_temp.size(); i++, pos++) {
                if (path_temp.at(i)->index == -10) {
                    path_temp.erase(path_temp.begin() + i);
                    break;
                }
            }

            finalPath.insert(finalPath.end(), path_temp.begin() + pos, path_temp.end());
            finalPath.insert(finalPath.end(), path_temp.begin(), path_temp.begin() + pos);

            IPath * path = new IPath();
            path->path = finalPath;
            path->calculateEval();

            layer->layerPath = path;
        }
    } else {
        layer->layerPath = NULL;
    }

    // free allocated memory
    if (row != NULL)
        free(row);
    if (colno != NULL)
        free(colno);

    if (lp != NULL) {
        //clean up such that all used memory by lpsolve is freed 
        delete_lp(lp);
    }
}

void createMiddleIPoint(IPoint * from, IPoint * to, IPoint * middle) {
    v3 newNormal = (from->normal + to->normal) / 2;
    v3 newPoint = ((from->point + to->point) / 2); // + newNormal * 20;
    IPoint m(newPoint, newNormal, -2);
    *middle = m;
}

void exportPreamble_inspection(FILE *f, int * line) {
    //    fprintf(f, "N1 G53\n");
    //    fprintf(f, "N2 G74 X2 Z1\n");
    //    fprintf(f, "N3 G01 Z15 F500\n");
    //    fprintf(f, "N4 G74 Y1 B2 ; home all axes\n");
    //    fprintf(f, "N5 G54\n");
    //    fprintf(f, "N6 G01 X200 Y0 Z50 B0 C0 F400\n");
    //    fprintf(f, "N7 M61=2500 ; turn on heated extruder\n");
    //    fprintf(f, "N8 M51=11000 ; turn on heated bed\n");
    //    fprintf(f, "N9 M40 ;wait for temperatures to settle\n");
    //    fprintf(f, "N10 G92 AV.A.ACT_POS.A\n");
    //    fprintf(f, "N11 G01 A-2.00000 F2400.00000 ; retract\n");
    if ((*line) == -1) {
        fprintf(f, " ; setup the inspection machine\n");
    } else {
        fprintf(f, "N%d ; setup the inspection machine\n", (*line)++);
    }
}

void exportEpilogue_inspection(FILE *f, int * line) {
    //    fprintf(f, "N%d M41=0 ; turn off heated chamber\n", line_number++);
    //    fprintf(f, "N%d M51=0 ; turn off heated bed\n", line_number++);
    //    fprintf(f, "N%d M61=0 ; turn off heated extruder\n", line_number++);
    //    fprintf(f, "N%d G53 G90\n", line_number++);
    //    fprintf(f, "N%d $IF V.A.ACT_POS.Z<0\n", line_number++);
    //    fprintf(f, "N%d     G01 Z0 F1000\n", line_number++);
    //    fprintf(f, "N%d $ENDIF\n", line_number++);
    //    fprintf(f, "N%d X101.271 F1000\n", line_number++);
    //    fprintf(f, "N%d Y-173.3\n", line_number++);
    //    fprintf(f, "N%d Z-150\n", line_number++);
    //    fprintf(f, "N%d Z-207 B0 C0 F500\n", line_number++);
    //    fprintf(f, "N%d Z-211.133 F100\n", line_number++);
    //    fprintf(f, "N%d M30\n", line_number++);
    if ((*line) == -1) {
        fprintf(f, "; end inspection machine\n");
    } else {
        fprintf(f, "N%d; end inspection machine\n", (*line)++);
    }
}

void exportG0command_inspection(FILE *f, float speed, IPoint * point, int * line) {
    if ((*line) == -1) {
        fprintf(f, "G00 F%.2f X%.2f Y%.2f Z%.2f \n", speed, point->point.x + inspection_char.adjust_x, point->point.y, point->point.z + inspection_char.adjust_z);
    } else {
        fprintf(f, "N%d G00 F%.2f X%.2f Y%.2f Z%.2f \n", (*line)++, speed, point->point.x + inspection_char.adjust_x, point->point.y, point->point.z + inspection_char.adjust_z);
    }
}

void exportM42activateInspectionCommand_inspection(FILE *f, int * line) {
    if ((*line) == -1) {
        fprintf(f, "M42P6S255 ; activate inspection\n");
    } else {
        fprintf(f, "N%d M42P6S255 ; activate inspection\n", (*line)++);
    }
}

void exportM42deactivateInspectionCommand_inspection(FILE *f, int * line) {
    if ((*line) == -1) {
        fprintf(f, "M42P6S0 ; deactivate inspection\n");
    } else {
        fprintf(f, "N%d M42P6S0 ; deactivate inspection\n", (*line)++);
    }
}

void exportM84deactivateMotorsCommand_inspection(FILE *f, int* line) {
    if ((*line) == -1) {
        fprintf(f, "M84 ; deactivate motors and wait for the end of the current instruction \n");
    } else {
        fprintf(f, "N%d M84 ; deactivate motors and wait for the end of the current instruction \n", (*line)++);
    }
}

void exportG1E1_MrotZ_inspection(FILE *f, float C, int* line) {
    if ((*line) == -1) {
        fprintf(f, "G1E1 %.2f ; MrotZ rotation in degrees \n", C);
    } else {
        fprintf(f, "N%d G1E1 %.2f ; MrotZ rotation in degrees \n", (*line)++, C);
    }
}

void exportG1E2_MrotY_inspection(FILE *f, float B, int * line) {
    if ((*line) == -1) {
        fprintf(f, "G1E2 %.2f ; MrotY rotation in degrees\n", B);
    } else {
        fprintf(f, "N%d G1E2 %.2f ; MrotY rotation in degrees\n", (*line)++, B);
    }
}

int exportGCode5Axes_layer(IPath * path, char *filename, int * line) {
    FILE *f = NULL;
    fopen_s(&f, filename, "a");
    if (!f) return 1;

    exportPreamble_inspection(f, line);

    vector<IPoint*> points = path->path;
    for (int i = 0; i < points.size(); i++) {
        float B = points[i]->b_angle;
        float C = points[i]->c_angle;

        exportG0command_inspection(f, inspection_char.speed, points[i], line);
        exportM84deactivateMotorsCommand_inspection(f, line);
        exportG1E1_MrotZ_inspection(f, C, line);
        exportM84deactivateMotorsCommand_inspection(f, line);
        exportG1E2_MrotY_inspection(f, B, line);
        exportM84deactivateMotorsCommand_inspection(f, line);
        if (points[i]->index != -2) {
            exportM42activateInspectionCommand_inspection(f, line);
            exportM84deactivateMotorsCommand_inspection(f, line);
            exportM42deactivateInspectionCommand_inspection(f, line);
            exportM84deactivateMotorsCommand_inspection(f, line);
        }
    }

    exportEpilogue_inspection(f, line);

    fclose(f);
    return 0;
}

int exportGCode5Axes(ILayers * layers, char * filename, int lineNumbering) {
    char ff[1024];
    strcpy(ff, filename);
    char * fff;
    if (inspection_char.debug) {
        if (!strncmp(ff + strlen(ff) - 6, ".gcode", 6))
            ff[strlen(ff) - 6] = 0;

        fff = (char*) createDirectoryFile((char*) ff, (char*) "Results/", (char*) ".scad");
        remove(fff);
        if (exportSLayersOpenSCadinspection(layers->stl_local, fff) == 1) {
            return 1;
        }
    }

    //    remove(filename);
    float eval = 0.0;

    int line = 0;
    if (lineNumbering == 0) {
        line = -1;
    }

    for (int layerPos = 0; layerPos < layers->inspectionLayers.size(); layerPos++) {
        ILayer * layer = layers->inspectionLayers.at(layerPos);
        IPath * originalPath = layer->layerPath;

        if (originalPath->path.size() == 0) {
            break;
        }

        if (exportGCode5Axes_layer(originalPath, filename, &line) == 1) {
            return 1;
        }
        if (inspection_char.action_print) {
            originalPath->print();
        }
        if (inspection_char.debug && layers->inspectionLayers.size() == 1) {
            originalPath->print(fff, 4);
        } else if (inspection_char.debug && layerPos == 0) {
            originalPath->print(fff, 2);
        } else if (inspection_char.debug && layerPos == layers->inspectionLayers.size() - 1) {
            originalPath->print(fff, 3);
        } else if (inspection_char.debug) {
            originalPath->print(fff, 0);
        }
        eval += originalPath->eval;
        if (layerPos != layers->inspectionLayers.size() - 1 && layers->inspectionLayers.at(layerPos + 1)->layerPath->path.size() != 0) {
            IPath * linkLayers = new IPath;
            linkLayers->path.push_back(originalPath->path.back());
            linkLayers->path.push_back(layers->inspectionLayers.at(layerPos + 1)->layerPath->path.front());
            linkLayers->eval = linkLayers->path.front()->evalTo(linkLayers->path.back());

            if (inspection_char.action_print) {
                printf("//up to the next layer\n");
                linkLayers->print();
            }
            if (inspection_char.debug) {
                linkLayers->print(fff, 1);
            }
            eval += linkLayers->eval;
            exportGCode5Axes_layer(linkLayers, filename, &line);
            delete linkLayers;
        }
    }

    printf("Total eval: %f\n", eval);
    printf("End creating GCode file\n");
    return 0;
}

void CombF(ILayers * layers) {
    IPoint * start = NULL;
    vector<IPath*> possibilities;
    for (int i = 0; i < layers->inspectionLayers.size(); i++) {
        if (i > 0) {
            start = findNextLayerInitialPoint(possibilities.front()->path.back(), layers->inspectionLayers.at(i));
            layers->inspectionLayers.at(i - 1)->layerPath = possibilities.front();
            //            int best = -1;
            //            float ev = numeric_limits<float>::max();
            //            for(int j = 0;j<possibilities.size();j++){
            //                IPoint * temp = findNextLayerInitialPoint(possibilities.at(j)->path.back(), layers->inspectionLayers.at(i));
            //                float e = temp->evalTo(possibilities.at(j)->path.back());
            //                if(e < ev){
            //                    if(best != -1){
            //                        delete possibilities.at(best);
            //                    }
            //                    start = temp;
            //                    ev = e;
            //                    best = j;
            //                    layers->inspectionLayers.at(i-1)->layerPath = possibilities.at(j);
            //                }else{
            //                    delete possibilities.at(j);
            //                }
            //            }
        }
        possibilities = CombF_generatePath(layers->inspectionLayers.at(i), start);
        //printf("possibilities %d\n",possibilities.size());


        //if (layers->inspectionLayers.at(i)->layerPath->path.size() == 0) {
        if (possibilities.size() == 0) {
            printf("Cannot find valid path. Try again with new parameters...\n");
            break;
        } else if (i == layers->inspectionLayers.size() - 1) {
            layers->inspectionLayers.at(i)->layerPath = possibilities.front();
        }
    }
}

vector<IPath*> CombF_generatePath(ILayer * layer, IPoint * startPoint) {
    int totalIPoint = layer->getValidIPointSize();

    vector<IPath*> results;

    //IPath * bestPath = new IPath();
    //bestPath->eval = numeric_limits<float>::max();

    vector<vector < IPoint*>> partialPaths;
    vector<vector<int>> partialPathsVisited;

    if (startPoint != NULL) {
        IPoint * i = i = startPoint;
        for (auto next : i->adjacentList) {
            vector<IPoint*> p;
            IPoint * j = next;
            p.push_back(i);
            p.push_back(j);
            partialPaths.push_back(p);
        }
    } else {
        //create possible paths from every point
        for (auto polygon : layer->inspectionPolygons) {
            for (auto isp : polygon->subPaths) {
                for (auto i : isp->subPath) {
                    for (auto next : i->adjacentList) {
                        vector<IPoint*> p;
                        IPoint * j = next;
                        p.push_back(i);
                        p.push_back(j);
                        partialPaths.push_back(p);
                    }
                }
            }
        }
    }

    //initialization
    int size = layer->getBiggestIndex() + 1;
    for (auto points : partialPaths) {
        vector<int> temp;
        for (int y = 0; y < size; y++) {
            if (y == points.front()->index || y == points.back()->index) {
                temp.push_back(1);
            } else {
                temp.push_back(0);
            }
        }
        partialPathsVisited.push_back(temp);
    }

    while (partialPaths.size() != 0) {
        vector<IPoint*> currentPartialPath = partialPaths.back();
        partialPaths.pop_back();

        vector<int> currentPartialVisited = partialPathsVisited.back();
        partialPathsVisited.pop_back();

        for (auto next : currentPartialPath.back()->adjacentList) {
            if (currentPartialVisited.at(next->index) != 1) {
                vector<IPoint*> newPath(currentPartialPath);
                newPath.push_back(next);

                vector<int> newVisited(currentPartialVisited);
                newVisited.at(next->index) = 1;

                if (newPath.size() == totalIPoint) {
                    IPath * path = new IPath();
                    path->path = newPath;

                    path->calculateEval();


                    int pos = 0;
                    for (int i = 0; i < results.size(); i++, pos++) {
                        if (path->eval < results.at(i)->eval) {
                            break;
                        }
                    }

                    results.insert(results.begin() + pos, path);

                    //                    if (path->eval < bestPath->eval) {
                    //                        delete bestPath;
                    //                        bestPath = path;
                    //                    } else {
                    //                        delete path;
                    //                    }
                } else {
                    partialPaths.push_back(newPath);
                    partialPathsVisited.push_back(newVisited);
                }
            }
        }
    }
    if (results.size() != 0) {
        float best = results.front()->eval;
        int pos = 0;
        for (int i = 0; i < results.size(); i++, pos++) {
            if (best < results.at(i)->eval) {
                break;
            }
        }

        for (int w = results.size() - 1; w > pos; w--) {
            delete results.at(w);
            results.erase(results.begin() + w);
        }
    }
    //layer->layerPath = bestPath;
    return results;
}