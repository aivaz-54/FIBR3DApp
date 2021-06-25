#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <set>
#include "admesh/stl.h"
#include "FIBR3Dapp.hpp"

extern printer printer_char;
float last_extrusion_speed = -1.0f; // last speed exported
float last_extrusion_rate = -1.0f; // last rate exported
float total_extruded = 0.0f; // total extruded from begining
int line_number;
float lastCposition=0.0f; // last C printer position 

void ANGLES(float vx, float vy, float vz, float *B, float *C, float *lastb, float *lastc) {
//    float tmpC;

    //*B=0.75f;
    //*C=0.0f;
    //return;
#define TOLF 1e-8f

    *B = acos(vz);
    
//    if (*B < printer_char.min_B)
//        *B=printer_char.min_B;
//    else if(*B > printer_char.max_B){
//        *B=printer_char.max_B;
//    }

    if (fabs(*B) > DEG_TO_RAD(1)) {
        // B is not equal to zero
        if (fabs(vx) <= TOLF) {
            // vx is zero,
            if (fabs(vy) <= TOLF) {
                *C = *lastc;
                *B = 0;
            } else {
                if (vy > 0)
                    *C = PI / 2;
                else
                    *C = -PI / 2;
            }
        } else {
            // vx != 0
            if (vx > 0) {
                if (fabs(vy) <= TOLF) {
                    *C = PI;
                } else if (vy > 0) {
                    *C = PI - atan(vy / vx);
                } else {
                    *C = -PI - atan(vy / vx);
                }
            } else {
                if (fabs(vy) <= TOLF) {
                    *C = *lastc;
                    *B = 0;
                } else {
                    *C = -atan(vy / vx);
                }
            }
        }
    } else {
        // if B==0 then we have an infinite possibility for C, but the best choice is to keep it as lastC
        *C = *lastc;
    }

    return;
}

void ANGLES_Brim(float x, float y, float *C, float *lastc) {

    if (x == 0) {
        if (y == 0) {
            *C = *lastc;
        } else {
            if (y > 0)
                *C = PI / 2;
            else
                *C = -PI / 2;
        }
    } else {
        // x != 0
        if (x > 0) {
            if (y == 0) {
                *C = PI;
            } else if (y > 0) {
                *C = PI - atan(y / x);
            } else {
                *C = -PI - atan(y / x);
            }
        } else {
            if (y == 0)
                *C = 0;
            else
                *C = -atan(y / x);
        }
    };

}


///* adapts cartesian coordinates to spherical coordinates (B,C) */
void carttransform(float x, float y, float z, float B, float C, float *nx, float *ny, float *nz, float *lastC, float *lastB) {
    // do nothing
    //    *nx = x;
    //    *ny = y;
    //    *nz = z;
    //
    //    return;

    // do some sanity checkup
    if (B < printer_char.min_B || B > printer_char.max_B){
        //printf("B angle out of printer range!\n");
        B = *lastB; // last sucessful B
    }

    // add printer_char.plate_offset to z, and remove it later
    *nx = x * cos(B) * cos(C) - y * sin(C) * cos(B) + (z + printer_char.plate_offset)* sin(B);
    *ny = x * sin(C) + y * cos(C);
    *nz = -x * sin(B) * cos(C) + y * sin(B) * sin(C) + (z + printer_char.plate_offset) * cos(B) - printer_char.plate_offset;

    if (*nx < printer_char.min_axis.x || *nx > printer_char.max_axis.x || *ny < printer_char.min_axis.y ||
            *ny > printer_char.max_axis.y || *nz < printer_char.min_axis.z || *nz > printer_char.max_axis.z) {
        // check to see if we can fix it
        B = *lastB;
        *nx = x * cos(B) * cos(C) - y * sin(C) * cos(B) + (z + printer_char.plate_offset)* sin(B);
        *ny = x * sin(C) + y * cos(C);
        *nz = -x * sin(B) * cos(C) + y * sin(B) * sin(C) + (z + printer_char.plate_offset) * cos(B) - printer_char.plate_offset;
        if (*nx < printer_char.min_axis.x || *nx > printer_char.max_axis.x || *ny < printer_char.min_axis.y ||
                *ny > printer_char.max_axis.y || *nz < printer_char.min_axis.z || *nz > printer_char.max_axis.z) {
            C = *lastC;
            *nx = x * cos(B) * cos(C) - y * sin(C) * cos(B) + (z + printer_char.plate_offset)* sin(B);
            *ny = x * sin(C) + y * cos(C);
            *nz = -x * sin(B) * cos(C) + y * sin(B) * sin(C) + (z + printer_char.plate_offset) * cos(B) - printer_char.plate_offset;
        }
    }

    // update last sucessful position
    *lastB = B;
    *lastC = C;
    
//    if (*nx < printer_char.min_axis.x || *nx > printer_char.max_axis.x || *ny < printer_char.min_axis.y ||
//            *ny > printer_char.max_axis.y || *nz < printer_char.min_axis.z || *nz > printer_char.max_axis.z ||
//            *lastB < printer_char.min_B || *lastB > printer_char.max_B)
//        printf("Point is out of printer range!\n");

    

    // clockwise rotation ...
    /*
     *nx = x*cos(C)*cos(B) - y*sin(C)*cos(B) + z*sin(B);
     *ny = x*sin(C) + y*cos(C);
     *nz = -x*sin(B)*cos(C) + y*sin(C)*sin(B) + z*cos(B);
     */

    return;
}


bool check_outer_spline(SPolygon spoly, v3 point) {
    double minX, maxX, minY, maxY;
    bool outside;
    int ss, L;
    vector<v3> Ptmp;

    // check if point is out on every spline
    spoly.get_bounding_box(&minX, &maxX, &minY, &maxY);
    if (point.x > minX && point.x < maxX && point.y > minY && point.y < maxY) {
        // inside of the outer bounding box
        outside = true; // assume outside
        
        // join all points for the polygon
        //vector<Spline> splines_outer = spoly.getSplines();
        //ss=splines_outer.size();
        vector<v3> Points = spoly.get_spline_points(); //splines_outer[0].get_y_data();
        //for (int k = 1; k < ss; k++) {
        //    Ptmp = splines_outer[k].get_y_data();
        //    L = Ptmp.size();
        //    for (int i = 1; i < L; i++)
        //        Points.push_back(Ptmp[i]);
        //}

        // points are connected (first == last)
        L = Points.size();
        for (int i = 0; i < L - 1; i++) {
            if ((Points[i].y > point.y) != (Points[i + 1].y > point.y) &&
                    (point.x < (Points[i].x - Points[i + 1].x)*(point.y - Points[i + 1].y) / (Points[i].y - Points[i + 1].y) + Points[i + 1].x)) {
                outside = !outside;
            }
        }

//        if (outside)
//            printf("Outside!\n");

        return outside; // return true if is outside (oposite of inside)
    }

    // not inside of bounding box
    return true; // is out
}





// output a single MATLAB file with the slices, layed out in 3D form

int exportSingleMATLABFormat3D(
        const std::vector<std::vector<LineSegment>> &slicesWithLineSegments,
        char *filename) {

    FILE *f = NULL;
    v3 p0, p1;
    int i, j;
    const int nSlices = (int) slicesWithLineSegments.size();
    std::vector<LineSegment> lss;


    fopen_s(&f, filename, "w");
    if (!f) return 1;

    fprintf(f, "figure;\nhold on;\n");

    for (i = 0; i < nSlices; ++i) {
        fprintf(f, "%% Slice %d\n\n", i);
        lss = slicesWithLineSegments[i];
        for (j = 0; j < (int) lss.size(); ++j) {
            p0 = lss[j].v[0];
            p1 = lss[j].v[1];
            fprintf(f, "\nline([%f,%f],[%f,%f],[%f,%f]);\n", p0.x, p1.x, p0.y, p1.y, p0.z, p1.z);
        }
    }
    fclose(f);
    return 0;
}


// output a single MATLAB file with the slices, layed out in 3D form

int exportFeatureEdgesMATLABFormat3D(
        const std::vector<feature_edge> &feature_edges,
        char *filename_convex,
        char *filename_concave) {

    FILE *f_convex = NULL, *f_concave = NULL;
    int i;
    const int nEdges = (int) feature_edges.size();


    fopen_s(&f_convex, filename_convex, "w");
    fopen_s(&f_concave, filename_concave, "w");
    if (!f_convex || !f_concave) return 1;

    fprintf(f_convex, "figure;\nhold on;\n");
    fprintf(f_concave, "figure;\nhold on;\n");

    for (i = 0; i < nEdges; ++i) {
        if (feature_edges[i].convex) {
            fprintf(f_convex, "%% Convex edge %d\n\n", i);
            fprintf(f_convex, "\nline([%f,%f],[%f,%f],[%f,%f]);\n", feature_edges[i].p[0].x, feature_edges[i].p[1].x,
                    feature_edges[i].p[0].y, feature_edges[i].p[1].y,
                    feature_edges[i].p[0].z, feature_edges[i].p[1].z);
        } else {
            fprintf(f_concave, "%% Convex edge %d\n\n", i);
            fprintf(f_concave, "\nline([%f,%f],[%f,%f],[%f,%f]);\n", feature_edges[i].p[0].x, feature_edges[i].p[1].x,
                    feature_edges[i].p[0].y, feature_edges[i].p[1].y,
                    feature_edges[i].p[0].z, feature_edges[i].p[1].z);
        }
    }

    fclose(f_convex);
    fclose(f_concave);
    return 0;
}


// output a single MATLAB file with the slices, layed out in 3D form

int exportLayersMATLABFormat3D(
        std::vector<Layer> &Layers,
        char *filename) {

    FILE *f = NULL;
    int i, j, k, ls, ps, pos;
    vector<Polygon> poly;
    vector<v3> points, normals;


    fopen_s(&f, filename, "w");
    if (!f) return 1;

    ls = (int) Layers.size();
    for (i = 0; i < ls; ++i) {
        poly = Layers[i].getPoly();
        fprintf(f, "Layer%d=[", i);
        ps = (int) poly.size();
        for (j = 0; j < ps; ++j) {
            points = poly[j].getPoints();
            //fprintf(f, "Area:%f\n", poly[j].getArea());
            pos = (int) points.size();
            for (k = 0; k < pos; ++k)
                fprintf(f, "%f,%f,%f;", points[k].x, points[k].y, points[k].z);
            fprintf(f, "NaN,NaN,NaN;");
        }
        fprintf(f, "];\n\n");
        fprintf(f, "figure;mapshow(Layer%d(:,1),Layer%d(:,2),'DisplayType','polygon');\n\n", i, i);
        fprintf(f, "hold on;\n\n");
        ps = (int) poly.size();
        for (j = 0; j < ps; ++j) {
            normals = poly[j].getNormals();
            points = poly[j].getPoints();
            //fprintf(f, "Area:%f\n", poly[j].getArea());
            pos = (int) points.size();
            fprintf(f, "plot3(%f,%f,%f,'ok');\n", points[0].x, points[0].y, points[0].z);
            for (k = 0; k < pos; ++k)
                fprintf(f, "line([%f,%f],[%f,%f],[%f,%f]);", points[k].x, points[k].x + normals[k].x, points[k].y, points[k].y + normals[k].y, 0.0f, normals[k].z); // mapshow uses z=0
        }
        fprintf(f, "\n\n");

    }
    fclose(f);
    return 0;
}




// output a single MATLAB file with the slices, layed out in 3D form

int exportSLayersMATLABFormat3D(
        std::vector<SLayer> &SLayers,
        char *filename) {

    FILE *f = NULL;
    int i, j, k, ls, ps, ss, m, pos;
    vector<SPolygon> spoly;
    vector<Spline> splines;
    vector<v3> points, normals;
    v3 y, head_normal;


    fopen_s(&f, filename, "w");
    if (!f) return 1;

    ls = (int) SLayers.size();
    for (i = 0; i < ls; ++i) {
        // for each Layer
        // get the polygon
        spoly = SLayers[i].getSPolygon();
        fprintf(f, "Layer%d=[", i);
        // get the number of Splines
        ps = (int) spoly.size();
        for (j = 0; j < ps; ++j) { // ps em vez de 1
            // get splines
            splines = spoly[j].getSplines();
            ss = (int) splines.size();
            for (k = 0; k < ss; ++k) { // ss em vez de 1
                // for all splines
                //if (!splines[i].is_cubic()) {
                /*for (t = splines[k].get_front(); t <= splines[k].get_back(); t += 0.5f) {
                        y = spline_eval(splines[k], t, &segment);
                        fprintf(f, "%f,%f,%f,%f;", t, y.x, y.y, y.z);
                }*/
                points = splines[k].get_y_data();

                pos = (int) points.size();
                for (m = 0; m < pos - 1; ++m) // spline is connected to the next spline, so use pos-1
                    fprintf(f, "%f,%f,%f;", points[m].x, points[m].y, points[m].z);
                //}
            }
            // end of polygon
            fprintf(f, "NaN,NaN,NaN;");
        }
        fprintf(f, "];\n\n");
        fprintf(f, "figure;mapshow(Layer%d(:,1),Layer%d(:,2),'DisplayType','polygon');hold on;axis([-170 170 -170 170])\n\n", i, i);

        ps = (int) spoly.size();
        for (j = 0; j < ps; ++j) { // ps em vez de 1
            // get splines
            splines = spoly[j].getSplines();
            ss = (int) splines.size();
            points = splines[0].get_y_data();
            fprintf(f, "plot3(%f,%f,%f,'ok');\n", points[0].x, points[0].y, points[0].z);
            for (k = 0; k < ss; ++k) { // ss em vez de 1
                // for all splines
                //if (!splines[i].is_cubic()) {
                /*for (t = splines[k].get_front(); t <= splines[k].get_back(); t += 0.5f) {
                        y = spline_eval(splines[k], t, &segment);
                        fprintf(f, "%f,%f,%f,%f;", t, y.x, y.y, y.z);
                }*/
                //}
                normals = splines[k].get_normals();
                points = splines[k].get_y_data();

                pos = (int) points.size();
                for (m = 0; m < pos; ++m) {
                    fprintf(f, "line([%f,%f],[%f,%f],[%f,%f],'Color',[0 0 1]);\n", points[m].x, points[m].x + 2 * normals[m].x, points[m].y, points[m].y + 2 * normals[m].y, 0.0f, 2 * normals[m].z); // mapshow uses z=0
                    //fprintf(f, "X=[X;[%f,%f,%f]];\n", points[m].x, points[m].y, points[m].z); // mapshow uses z=0
                    //fprintf(f, "D=[D;[%f,%f,%f]];\n", normals[m].x, normals[m].y, normals[m].z); // mapshow uses z=0
                    if (m < pos - 1) {
                        head_normal = points[m + 1] - points[m];
                        head_normal = head_normal.crossproduct(normals[m]);
                        head_normal = head_normal / head_normal.norm();
                        if (head_normal.z < 0)
                            head_normal = v3(-head_normal.x, -head_normal.y, -head_normal.z);
                        fprintf(f, "line([%f,%f],[%f,%f],[%f,%f],'Color',[1 0 0]);\n", points[m].x, points[m].x + 2 * head_normal.x, points[m].y, points[m].y + 2 * head_normal.y, 0.0f, 2 * head_normal.z); // mapshow uses z=0
                        //fprintf(f, "X=[X;[%f,%f,%f]];\n", points[m].x, points[m].y, points[m].z); // mapshow uses z=0
                        //fprintf(f, "D=[D;[%f,%f,%f]];\n", head_normal.x, head_normal.y, head_normal.z); // mapshow uses z=0
                    }
                }
            }
        }

        fprintf(f, "\n\n");
    }
    fclose(f);
    return 0;
}


// output a single MATLAB file with the slices, layed out in 3D form

int exportSLayersGCode3Axes(
        std::vector<SLayer> &SLayers,
        char *filename) {

    FILE *f = NULL;
    int i, j, k, l, ls, ps, ss, pos, segment;
    vector<SPolygon> spoly;
    vector<Spline> splines;
    vector<v3> points;
    char command[3] = "G1";


    fopen_s(&f, filename, "w");
    if (!f) return 1;

    fprintf(f, "G0 F1000\n");

    ls = (int) SLayers.size();
    for (i = 0; i < ls; ++i) {
        // for each Layer
        // get the polygon
        spoly = SLayers[i].getSPolygon();
        // get the number of SPolygons
        ps = (int) spoly.size();
        for (j = 0; j < ps; ++j) {
            if (spoly[j].getArea() < 0) { // we have a inner loop, so start with it
                // get splines
                splines = spoly[j].getSplines();
                ss = (int) splines.size();

                for (float offset = 0.0f; offset < 21.0f; offset += 20.5f) {
                    for (k = 0; k < ss; ++k) {
                        // for all splines

                        if (splines[k].is_cubic()) {
                            // follow arcs
                            vector<v3> Points = splines[k].get_y_data();
                            vector<float> x_data = splines[k].get_x_data();
                            vector<v3> Normals = splines[k].get_normals();
                            pos = (int) Points.size();
                            for (l = 0; l < pos; l++)
                                Points[l] = Normals[l] * offset + Points[l];
                            splines[k].set_y_data(Points);
                            computeSplineM(splines[k]);
                            for (l = 0; l < pos - 1; l++) {
                                v3 tmp = cubic(splines[k], (x_data[l] + x_data[l + 1]) / 2, &segment);
                                float ma = (tmp.y - Points[l].y) / (tmp.x - Points[l].x);
                                float mb = (Points[l + 1].y - tmp.y) / (Points[l + 1].x - tmp.x);
                                float x = (ma * mb * (Points[l].y - Points[l + 1].y) + mb * (Points[l].x + tmp.x) - ma * (tmp.x + Points[l + 1].x)) / (2 * (mb - ma));
                                float y = (-1 / ma)*(x - (Points[l].x + tmp.y) / 2) + (Points[l].y + tmp.y) / 2;
                                float radius2 = sqrt(pow(x - Points[l].x, 2.0f) + pow(y - Points[l].y, 2.0f));

                                fprintf(f, "G2 X%.2f Y%.2f Z%.2f R%.2f\n", Points[l + 1].x, Points[l + 1].y, Points[l + 1].z, radius2);
                            }
                        } else {
                            vector<v3> Points = splines[k].get_y_data();
                            vector<v3> Normals = splines[k].get_normals();
                            pos = (int) Points.size();
                            for (l = 0; l < pos; l++)
                                Points[l] = Normals[l] * offset + Points[l];
                            for (l = 0; l < pos - 1; l++) {
                                fprintf(f, "%s X%.2f Y%.2f Z%.2f\n", command, Points[l + 1].x, Points[l + 1].y, Points[l + 1].z);
                            }
                        }
                    }
                    // closing point -- needs to be checked
                }
            }
        }
    }

    fclose(f);
    return 0;
}

bool check_collisions(float *x, float *y, float *z, float *B, float *C) {
    // are we inside bounds?
    if (printer_char.min_axis.x > *x)
        *x = printer_char.min_axis.x;
    else if (printer_char.max_axis.x < *x)
        *x = printer_char.max_axis.x;

    if (printer_char.min_axis.y > *y)
        *y = printer_char.min_axis.y;
    else if (printer_char.max_axis.y < *y)
        *y = printer_char.max_axis.y;

    if (printer_char.min_axis.z > *z)
        *z = printer_char.min_axis.z;
    else if (printer_char.max_axis.z < *z)
        *z = printer_char.max_axis.z;

    if (printer_char.min_B > *B)
        *B = printer_char.min_B;
    else if (printer_char.max_B < *B)
        *B = printer_char.max_B;


    if (printer_char.min_C > *C)
        *C = printer_char.min_C;
    else if (printer_char.max_C < *C)
        *C = printer_char.max_C;

    // return false to print command
    return false;
}

void exportG12command(FILE *f, float extrusion_speed, float x, float y, float z, float B, float C, float radius2, float total_extrusion) {
    char command[1024];
    char command2[1024];

    //if(check_collisions(&x, &y, &z, &B, &C, &radius2))return; // check for collisions and change if necessary

    radius2 = 0.0f; // force to be a linear spline....
    command2[0]=0; // no command2 until we find that we need it
    
    //if(line_number>=100)
    //    printf("here");

    if (radius2 > 0.0f)
        snprintf(command, 1024, "G02");
    else
        snprintf(command, 1024, "G01");

    if (extrusion_speed != last_extrusion_speed) {
        snprintf(command, 1024, "%s F%.2f", command, extrusion_speed);
        last_extrusion_speed = extrusion_speed;
    }

    //std::cout << "Last C=" << lastCposition << " C=" << C << endl;
    
    if(C-lastCposition>PI) // correct a bug in the simulator
        snprintf(command2, 1024, "%s X%.2f Y%.2f Z%.2f B%.2f C%.2f", command, x, y, z, RAD_TO_DEG(B), RAD_TO_DEG(C)-360.0f); // -180
    else if (C-lastCposition<-PI)
            snprintf(command2, 1024, "%s X%.2f Y%.2f Z%.2f B%.2f C%.2f", command, x, y, z, RAD_TO_DEG(B), RAD_TO_DEG(C)+360.0f); // +180
    
    snprintf(command, 1024, "%s X%.2f Y%.2f Z%.2f B%.2f C%.2f", command, x, y, z, RAD_TO_DEG(B), RAD_TO_DEG(C)); // +180

    lastCposition=C;
    
    if (radius2 > 0.0f) {
        snprintf(command, 1024, "%s R%.2f ", command, radius2);
        if(command2[0])
            snprintf(command2, 1024, "%s R%.2f ", command2, radius2);
    }


    //if (extrusion_rate != last_extrusion_rate) {
    snprintf(command, 1024, "%s A%.2f ", command, total_extrusion);
    if(command2[0])
        snprintf(command2, 1024, "%s A%.2f ", command2, total_extrusion);
    //		last_extrusion_rate = extrusion_rate;
    //}

    if(command2[0]) // command 2 first to avoid bug
        fprintf(f, "; G1 extra command to avoid simulator bug\nN%d %s\n", line_number++, command2);
    fprintf(f, "N%d %s\n", line_number++, command);
    //fprintf(f,"P=[P;[%.2f %.2f %.2f %.2f %.2f]];\n", x, y, z, RAD_TO_DEG(B), RAD_TO_DEG(C));
    
    fflush(f);

}

void exportG0command(FILE *f, float speed, float x, float y, float z, float B, float C) {
    
    //if(line_number>=24941)
    //    printf("here");
    
    if(C-lastCposition>PI){ // correct a bug in the simulator
        fprintf(f, "; G0 extra command to correct bug\nN%d G00 F%.2f X%.2f Y%.2f Z%.2f B%.2f C%.2f\n", line_number++, speed, x, y, z, RAD_TO_DEG(B), -RAD_TO_DEG(C)); // 180
        // bug in the simulator
        fprintf(f, "; G0 extra command to correct gap\nN%d G01 F%.2f X%.2f Y%.2f Z%.2f B%.2f C%.2f A%.2f\n", line_number++, speed, x, y, z, RAD_TO_DEG(B), -RAD_TO_DEG(C), total_extruded); // 180
    } else if (C-lastCposition<-PI){
        fprintf(f, "; G0 extra command to correct bug\nN%d G00 F%.2f X%.2f Y%.2f Z%.2f B%.2f C%.2f\n", line_number++, speed, x, y, z, RAD_TO_DEG(B), RAD_TO_DEG(C)+360.0f);
        // bug in the simulator
        fprintf(f, "; G0 extra command to correct gap\nN%d G01 F%.2f X%.2f Y%.2f Z%.2f B%.2f C%.2f A%.2f\n", line_number++, speed, x, y, z, RAD_TO_DEG(B), RAD_TO_DEG(C)+360.0f, total_extruded);
    }

    fprintf(f, "N%d G00 F%.2f X%.2f Y%.2f Z%.2f B%.2f C%.2f\n", line_number++, speed, x, y, z, RAD_TO_DEG(B), RAD_TO_DEG(C));
    // bug in the simulator
    fprintf(f, "; G0 extra command to correct bug\nN%d G01 F%.2f X%.2f Y%.2f Z%.2f B%.2f C%.2f A%.2f\n", line_number++, speed, x, y, z, RAD_TO_DEG(B), RAD_TO_DEG(C), total_extruded);

    lastCposition=C;
    
    fflush(f);
}

void exportPreamble(FILE *f) {

    fprintf(f, "N1 G53\n");
    fprintf(f, "N2 G74 X2 Z1\n");
    fprintf(f, "N3 G01 Z15 F500\n");
    fprintf(f, "N4 G74 Y1 B2 ; home all axes\n");
    fprintf(f, "N5 G54\n");
    fprintf(f, "N6 G01 X200 Y0 Z50 B0 C0 F400\n");
    fprintf(f, "N7 M61=2500 ; turn on heated extruder\n");
    fprintf(f, "N8 M51=11000 ; turn on heated bed\n");
    fprintf(f, "N9 M40 ;wait for temperatures to settle\n");
    fprintf(f, "N10 G92 AV.A.ACT_POS.A\n");
    fprintf(f, "N11 G01 A-2.00000 F2400.00000 ; retract\n");

    line_number = 12;
    total_extruded=0.0f;

}

void exportFullBrim(FILE *f, float part_size, float z_offset) {
    float t, max_t, A = 0.0f, B = 0.0f, C = 0.0f, lastC = 0.0f, lastB = 0.0f;
    float brim_step = printer_char.extrusion_width * 2.0f;
    float x, y, lastx = 0.0f, lasty = 0.0f, n_x, n_y, n_z;
    float distance;

    // make a spiral brim of the size of the part
    max_t = (part_size/2.0f + printer_char.brim_size)*2.0f * PI / printer_char.extrusion_width;


    // go to center
    
    exportG0command(f, printer_char.speed, printer_char.extrusion_width, 0.0f, z_offset, B, C);

    for (t = printer_char.extrusion_width; t < max_t; t += 2.0f*PI*brim_step/(t*printer_char.extrusion_width)) {
        x = printer_char.extrusion_width * t * sin(t) / (2.0f * PI);
        y = printer_char.extrusion_width * t * cos(t) / (2.0f * PI);
        distance = sqrt(pow(lastx - x, 2.0f) + pow(lasty - y, 2.0f));
        total_extruded += distance * printer_char.extrusion_multiplier * (printer_char.extrusion_height * (printer_char.extrusion_width - printer_char.extrusion_height) + PI * pow(printer_char.extrusion_height / 2.0f, 2.0f)) /
                (PI * pow(printer_char.filament_radius, 2.0f) / 4.0f);
        lastx = x;
        lasty = y;
        //ANGLES_Brim(x, y, &C, &lastC);
        carttransform(x, y, z_offset, B, C, &n_x, &n_y, &n_z, &lastC, &lastB);
        exportG12command(f, printer_char.extrusion_speed, n_x, n_y, n_z, B, C, 0.0, total_extruded);
    }
 
}

void exportEpilogue(FILE *f) {

    fprintf(f, "N%d M41=0 ; turn off heated chamber\n", line_number++);
    fprintf(f, "N%d M51=0 ; turn off heated bed\n", line_number++);
    fprintf(f, "N%d M61=0 ; turn off heated extruder\n", line_number++);
    fprintf(f, "N%d G53 G90\n", line_number++);
    fprintf(f, "N%d $IF V.A.ACT_POS.Z<0\n", line_number++);
    fprintf(f, "N%d     G01 Z0 F1000\n", line_number++);
    fprintf(f, "N%d $ENDIF\n", line_number++);
    fprintf(f, "N%d X101.271 F1000\n", line_number++);
    fprintf(f, "N%d Y-173.3\n", line_number++);
    fprintf(f, "N%d Z-150\n", line_number++);
    fprintf(f, "N%d Z-207 B0 C0 F500\n", line_number++);
    fprintf(f, "N%d Z-211.133 F100\n", line_number++);
    fprintf(f, "N%d M30\n", line_number++);

}

float expandSLayers(std::vector<SLayer> &SLayers, std::vector<SLayer> &ESLayers, float offset) {
    int i, j, k, l, ls, ps, ss, pos;
    vector<SPolygon> spoly;
    vector<Spline> splines;
    float x_end = 0.0f;

    // empty ESLayers
    ESLayers.clear();

    // go for all layers
    ls = (int) SLayers.size();
    for (i = 0; i < ls; ++i) {
        SLayer tmplay;
        // for each Layer
        // get the polygon
        spoly = SLayers[i].getSPolygon();
        // get the number of SPolygons
        ps = (int) spoly.size();
        // go for all polygons
        for (j = 0; j < ps; ++j) {
            if (spoly[j].getArea() < 0) { // we have a inner loop, so go for it
                SPolygon tmppoly;
                float arc_length = 0.0f;

                splines = spoly[j].getSplines(); // get splines
                ss = (int) splines.size();

                // for all splines
                for (k = 0; k < ss; ++k) {
                    vector<v3> Points_orig = splines[k].get_y_data();
                    vector<v3> Normals = splines[k].get_normals();
                    vector<float> x_data;
                    vector<v3> Points;
                    v3 tmpd, tmpp, tmppprev, normal;

                    // Adapt all points to the new offset
                    normal = Normals[0];
                    normal.z = 0;
                    normal = normal / normal.norm();
                    tmppprev = Points_orig[0] - normal * offset; // repeat first so arc_length=0
                    // Adapt all points to the new offset
                    pos = (int) Points_orig.size();
                    for (l = 0; l < pos; l++) {
                        normal = Normals[l];
                        normal.z = 0;
                        normal = normal / normal.norm();
                        tmpp = Points_orig[l] - normal * offset;
                        Points.push_back(tmpp);
                        tmpd = tmpp - tmppprev;
                        arc_length += tmpd.norm();
                        x_data.push_back(arc_length);
                        tmppprev = tmpp;
                    }

                    splines[k].set_y_data(Points);
                    splines[k].set_x_data(x_data);

                    if (x_end < arc_length)
                        x_end = arc_length;

                    if (splines[k].is_cubic())
                        computeSplineM(splines[k]);
                    
                    //splines[k].set_bounding_box();

                    tmppoly.push_back(splines[k]);
                }
                tmppoly.computeArea();
                tmppoly.set_bounding_box();
                tmplay.push_back(tmppoly);
                //espoly.push_back(tmppoly);
            } else {
                //SPolygon tmppoly;
                
                //splines = spoly[j].getSplines(); // get splines
                //ss = (int) splines.size();

                // for all splines
                //for (k = 0; k < ss; ++k){
                //    splines[k].set_bounding_box();
                //    tmppoly.push_back(splines[k]);
                //}
                
                //tmppoly.computeArea();
                spoly[j].set_bounding_box();
                tmplay.push_back(spoly[j]);
                //espoly.push_back(tmppoly); // outer spoly, just copy
            }
        }
        ESLayers.push_back(tmplay);
    }

    return x_end;
}

void expandSimpleSpline(std::vector<SimpleSpline> &SSpline, std::vector<SimpleSpline> &ESSpline, float offset) {
    int i, k, ps, ss;

    // empty ESLayers
    ESSpline.clear();

    // go for all Splines
    ss = (int) SSpline.size();
    for (i = 0; i < ss; ++i) {
        SimpleSpline S = SSpline[i], tmpS;

        vector<v3> Points_orig = S.get_y_data();
        vector<v3> Normals = S.get_normals();
        vector<SPolygon*> SPoly = S.get_spoly();

        tmpS.push_back(S.get_x_data()); // save x_data

        ps = Points_orig.size();
        // for all points
        for (k = 0; k < ps; ++k) {
            v3 normal = Normals[k];
            normal.z = 0;
            normal = normal / normal.norm();
            v3 tmpp=Points_orig[k];// - normal * offset;
            tmpS.push_back(tmpp, Normals[k], SPoly[k]);
        }

        ESSpline.push_back(tmpS);
    }

}


float buildSSpline(std::vector<SLayer> &SLayers, vector<SimpleSpline> &NSSpline, float radius_offset) {
    int i, j, k, ls, ps, ss, pos, segment;
    vector<SPolygon> spoly;
    vector<Spline> splines;
    SimpleSpline tmpSpline;
    float t, x_start=0.0f, x_end=0.0f;
    bool nonstop=true;

    // start at x_start (0.0f) and go to the longest spline
    for (t = x_start; nonstop; t += radius_offset) {
        nonstop=false; // we will stop if no point is found
        tmpSpline.clear();
        tmpSpline.push_back(t);
        // go for all layers, from top to bottom (since we have finished at top)
        ls = (int) SLayers.size();
        for (i = ls - 1; i >= 0; i--) {
            // for each Layer
            // get the polygon
            spoly = SLayers[i].getSPolygon();
            // get the number of SPolygons
            ps = (int) spoly.size();
            // for all polygons
            for (j = 0; j < ps; ++j) {
                if (spoly[j].getArea() < 0) {// we have the inner surface
                    splines = spoly[j].getSplines(); // splines to follow
                    ss = (int) splines.size();
                    // for all splines
                    for (k = 0; k < ss; ++k) {
                        if (t >= splines[k].get_front() && t <= splines[k].get_back()) {
                            // found it
                            nonstop=true;
                            //vector<v3> Points = splines[k].get_y_data();
                            vector<v3> Normals = splines[k].get_normals();
                            //pos = (int) Points.size();

                            v3 point = spline_eval(splines[k], t, &segment);
                            v3 normal = Normals[segment - 1];
                            
                            tmpSpline.push_back(point, normal, SLayers[i].getSPolygon((j+1)%ps));
                            break;
                        }
                    }
                }
            }
        }
        NSSpline.push_back(tmpSpline);
    }
    x_end=t-radius_offset; // last t used
    return(x_end);
}


void exportBrim(FILE *f, std::vector<SLayer> &SLayers, float offset, float z_offset){
    std::vector<SLayer> ESLayers;
    vector<SPolygon> spoly;
    vector<Spline> splines;
    int j, k, l, ps, ss, pos, segment;
    float x, y, z, lastx, lasty, B=0.0f, C=0.0f, lastB=0.0f, lastC=0.0f, ma, mb, radius2, distance;

    expandSLayers(SLayers, ESLayers, offset);
    // repeat first inner layer
    spoly = ESLayers[0].getSPolygon();

    ps = (int) spoly.size();
    // go for all polygons
    for (j = 0; j < ps; ++j) {
        if (spoly[j].getArea() < 0) {
            vector<v3> Points;
            vector<float> x_data;

            splines = spoly[j].getSplines(); // get splines
            ss = (int) splines.size();

            Points = splines[0].get_y_data();
            x = Points[0].x;
            y = Points[0].y;
            z = z_offset;

            exportG0command(f, printer_char.speed, x, y, z, lastB, lastC);

            lastx = Points[0].x;
            lasty = Points[0].y;


            for (k = 0; k < ss; ++k) {
                Points = splines[k].get_y_data();
                x_data = splines[k].get_x_data();

                pos = (int) Points.size();
                if (pos < 2) continue; // a segment with just one point, should not happen

                if (false && splines[k].is_cubic()) {
                    // only complanar arcs are used!
                    // we have a cubic spline and arcs should be used

                    for (l = 1; l < pos; l++) {
                        v3 tmp = cubic(splines[k], (x_data[l - 1] + x_data[l]) / 2.0f, &segment);
                        // compute the spline radius
                        if ((tmp.x - Points[l - 1].x) != 0 && (Points[l].x - tmp.x) != 0) {
                            ma = (tmp.y - Points[l - 1].y) / (tmp.x - Points[l - 1].x);
                            mb = (Points[l].y - tmp.y) / (Points[l].x - tmp.x);
                            x = (ma * mb * (Points[l - 1].y - Points[l].y) + mb * (Points[l - 1].x + tmp.x) - ma * (tmp.x + Points[l].x)) / (2 * (mb - ma));
                            y = (-1 / ma)*(x - (Points[l - 1].x + tmp.y) / 2) + (Points[l - 1].y + tmp.y) / 2;
                            radius2 = sqrt(pow(x - Points[l - 1].x, 2.0f) + pow(y - Points[l - 1].y, 2.0f));
                        } else {
                            radius2 = 0.0f;
                        }

                        x = Points[l].x;
                        y = Points[l].y;
                        // z already computed

                        distance = sqrt(pow(lastx - x, 2.0f) + pow(lasty - y, 2.0f));
                        total_extruded += distance * printer_char.extrusion_multiplier * (printer_char.extrusion_height * (printer_char.extrusion_width - printer_char.extrusion_height) + PI * pow(printer_char.extrusion_height / 2.0f, 2.0f)) /
                                (PI * pow(printer_char.filament_radius, 2.0f) / 4.0f);
                        lastx = Points[l].x;
                        lasty = Points[l].y;

                        exportG12command(f, printer_char.extrusion_speed, x, y, z, lastB, lastC, radius2, total_extruded);
                    }
                } else {
                    // we have a linear spline
                    pos = (int) Points.size();
                    // go for all points in spline
                    for (l = 1; l < pos; l++) {
                        x = Points[l].x;
                        y = Points[l].y;
                        // z already computed

                        distance = sqrt(pow(lastx - x, 2.0f) + pow(lasty - y, 2.0f));
                        total_extruded += distance * printer_char.extrusion_multiplier * (printer_char.extrusion_height * (printer_char.extrusion_width - printer_char.extrusion_height) + PI * pow(printer_char.extrusion_height / 2.0f, 2.0f)) /
                                (PI * pow(printer_char.filament_radius, 2.0f) / 4.0f);
                        lastx = Points[l].x;
                        lasty = Points[l].y;

                        exportG12command(f, printer_char.extrusion_speed, x, y, z, lastB, lastC, 0.0f, total_extruded);
                    }
                }
            }
        }
    }
    
}


// output G code for 5-axis shell type objects
int exportSLayersGCode5Axes_shell(
        std::vector<SLayer> &SLayers,
        float min_z, float part_size, bool GCode5D,
        FILE *f, float z_offset) {

    //FILE *f = NULL;
    int i, j, k, l, ls, ps, ss, pos, segment, lay, lays, counter;
    std::vector<SLayer> ESLayers, E2SLayers;
    vector<SPolygon> spoly;
    vector<Spline> splines; //, splines_outer;
    bool first_gcode = true, out_of_object;
    float x_start = 0.0f, x_end = -1.0f, offset, t; // z_offset;
    float ma, mb, x, y, z, lastx = 0.0f, lasty = 0.0f, lastz = 0.0f, radius2 = 0.0f;
    float lastC = 0.0f, lastB = 0.0f, B = 0.0f, C = 0.0f, distance, radius_offset;
    vector<GCode> GCodes;
    char FileName[1024];
    vector<SLayer> *Pointer;


    //fopen_s(&f, filename, "w");
    if (!f) return 1;

    // preamble was previously exported
    total_extruded=0.0f; // do it here, if not done in the preamble
    
    // print the full brim
    if (printer_char.export_full_brim) {
        // export first brim (full base of object at printer base level)
        z_offset += printer_char.extrusion_height;
        fprintf(f, ";Start full brim\n");
        exportFullBrim(f, part_size, z_offset);
        fprintf(f, ";End full brim\n");
        // next z_offset
        z_offset += printer_char.extrusion_height;
    }

    // print the brim
    if (printer_char.export_brim) {
        radius_offset = printer_char.extrusion_width - printer_char.extrusion_height * (1 - PI / 4);
        fprintf(f, ";Start brim\n");
        for (offset = printer_char.extrusion_width / 2.0f; offset < printer_char.brim_size; offset += radius_offset) {
            exportBrim(f,SLayers,offset,z_offset);
        }
        fprintf(f, ";End brim\n");
        z_offset+=printer_char.extrusion_height;
    }

    // print the support
    if (printer_char.export_support) {
        radius_offset = printer_char.extrusion_width - printer_char.extrusion_height * (1 - PI / 4);
        fprintf(f, ";Start support\n");
        for (i = 0; i < (int) (printer_char.support_height / printer_char.extrusion_height); i++) {
            offset = printer_char.extrusion_width / 2.0f;
            exportBrim(f, SLayers, offset, z_offset);
            offset += radius_offset;
            exportBrim(f, SLayers, offset, z_offset);
            z_offset += printer_char.extrusion_height;
        }
        fprintf(f, ";End support\n");
    }
    

     fprintf(f, ";Start core\n");
    // build a core for the object
    // we assume that a core is always possible, i.e., part is thick enough
    /* Export first two object layers for support */
     
    // start in the middle of the first extrusion and do as many loops as requested
    offset = printer_char.extrusion_width / 2.0f;
    expandSLayers(SLayers, ESLayers, offset);
    radius_offset = printer_char.extrusion_width - printer_char.extrusion_height * (1 - PI / 4);
    //radius_offset=0.5-0.5*(1-PI/4);
    offset += radius_offset;
    expandSLayers(SLayers, E2SLayers, offset);

    //snprintf(FileName, 1024, "layers.m");
    //exportSLayersMATLABFormat3D(SLayers, FileName);
    //snprintf(FileName, 1024, "layers_1.m");
    //exportSLayersMATLABFormat3D(ESLayers, FileName);
    //snprintf(FileName, 1024, "layers_2.m");
    //exportSLayersMATLABFormat3D(E2SLayers, FileName);

    // go for all layers
    ls = (int) ESLayers.size();
    lay = 0;
    for (i = 0; i < ls; i++) {
        // for each Layer
        // get the polygon
        for (lays = 0; lays < printer_char.extrusion_number_of_layers; lays++) {
            //printf("Lay=%d\n", lay);
            if (lay) spoly = E2SLayers[i].getSPolygon();
            else spoly = ESLayers[i].getSPolygon();

            //fprintf(f,"Spline Layer %d\n", lay);
            // get the number of SPolygons
            ps = (int) spoly.size();
            // go for all polygons
            for (j = 0; j < ps; ++j) {
                if (spoly[j].getArea() < 0) { // we have an inner loop, so start with it
                    vector<v3> Points;
                    vector<v3> Normals;
                    vector<float> x_data;
                    v3 head_normal;

                    splines = spoly[j].getSplines(); // get splines
                    ss = (int) splines.size();

                    // get nozzel to beginning of spline with a G00 command
                    Points = splines[0].get_y_data();
                    Normals = splines[0].get_normals();
                    x_data = splines[0].get_x_data();
                    head_normal = Points[1] - Points[0];
                    head_normal = head_normal.crossproduct(Normals[0]);
                    head_normal = head_normal / head_normal.norm();
                    if (head_normal.z < 0)
                        head_normal = v3(-head_normal.x, -head_normal.y, -head_normal.z);

                    if (GCode5D) {
                        // compute B and C table angles
                        ANGLES(head_normal.x, head_normal.y, head_normal.z, &B, &C, &lastB, &lastC);
                        // tranform x,y,z coordinates to new axis (rotation at origin)
                        carttransform(Points[0].x, Points[0].y, Points[0].z + z_offset, B, C, &x, &y, &z, &lastC, &lastB); // - min_z + printer_char.min_axis.z + printer_char.plate_offset
                        // and change plate offset
                        //z -= printer_char.plate_offset;
                    } else {
                        x = Points[0].x;
                        y = Points[0].y;
                        z = Points[0].z + z_offset; //- min_z + printer_char.min_axis.z 
                    }

                    // lastB and lastC were updated in carttransform to be the actual B and C
                    exportG0command(f, printer_char.speed, x, y, z, lastB, lastC);
                    //fprintf(f,"D=[D;[%.2f %.2f %.2f %.2f %.2f]];\n", head_normal.x, head_normal.y, head_normal.z, RAD_TO_DEG(B), RAD_TO_DEG(C));
                    //fprintf(f,"X=[X;[%.2f %.2f %.2f %.2f %.2f]];\n", Points[0].x, Points[0].y, Points[0].z - min_z + printer_char.min_axis.z + printer_char.plate_offset, RAD_TO_DEG(B), RAD_TO_DEG(C));

                    // last point before transform
                    lastx = Points[0].x;
                    lasty = Points[0].y;
                    lastz = Points[0].z;

                    // go for all spline segments
                    for (k = 0; k < ss; ++k) {
                        Points = splines[k].get_y_data();
                        Normals = splines[k].get_normals();
                        x_data = splines[k].get_x_data();

                        pos = (int) Points.size();
                        if (pos < 2) continue; // a segment with just one point, should not happen

                        if (false && splines[k].is_cubic()) {
                            // only complanar arcs are used!
                            // we have a cubic spline and arcs should be used

                            for (l = 1; l < pos; l++) {
                                v3 tmp = cubic(splines[k], (x_data[l - 1] + x_data[l]) / 2.0f, &segment);
                                // compute the spline radius
                                if ((tmp.x - Points[l - 1].x) != 0 && (Points[l].x - tmp.x) != 0) {
                                    ma = (tmp.y - Points[l - 1].y) / (tmp.x - Points[l - 1].x);
                                    mb = (Points[l].y - tmp.y) / (Points[l].x - tmp.x);
                                    x = (ma * mb * (Points[l - 1].y - Points[l].y) + mb * (Points[l - 1].x + tmp.x) - ma * (tmp.x + Points[l].x)) / (2 * (mb - ma));
                                    y = (-1 / ma)*(x - (Points[l - 1].x + tmp.y) / 2) + (Points[l - 1].y + tmp.y) / 2;
                                    radius2 = sqrt(pow(x - Points[l - 1].x, 2.0f) + pow(y - Points[l - 1].y, 2.0f));
                                } else {
                                    radius2 = 0.0f;
                                }
                                // compute the perpendicular to the normal
                                head_normal = Points[l] - Points[l - 1];
                                head_normal = head_normal.crossproduct(Normals[l]);
                                head_normal = head_normal / head_normal.norm();
                                if (head_normal.z < 0) // are we pointing right?
                                    head_normal = v3(-head_normal.x, -head_normal.y, -head_normal.z);

                                distance = sqrt(pow(lastx - Points[l].x, 2.0f) + pow(lasty - Points[l].y, 2.0f) + pow(lastz - Points[l].z, 2.0f));
                                total_extruded += distance * printer_char.extrusion_multiplier * (printer_char.extrusion_height * (printer_char.extrusion_width - printer_char.extrusion_height) + PI * pow(printer_char.extrusion_height / 2.0f, 2.0f)) /
                                        (PI * pow(printer_char.filament_radius, 2.0f) / 4.0f);
                                
                                // compute spherical if 5 axis
                                if (GCode5D) {
                                    ANGLES(head_normal.x, head_normal.y, head_normal.z, &B, &C, &lastB, &lastC);
                                    carttransform(Points[l].x, Points[l].y, Points[l].z + z_offset, B, C, &x, &y, &z, &lastC, &lastB); //- min_z + printer_char.min_axis.z + printer_char.plate_offset
                                    //z -= printer_char.plate_offset;
                                } else {
                                    // A and C are already set to zero
                                    x = Points[l].x;
                                    y = Points[l].y;
                                    z = Points[l].z + z_offset; // - min_z + printer_char.min_axis.z
                                }
                                lastx = Points[l].x;
                                lasty = Points[l].y;
                                lastz = Points[l].z;
                                exportG12command(f, printer_char.extrusion_speed, x, y, z, lastB, lastC, radius2, total_extruded);
                                //fprintf(f,"D=[D;[%.2f %.2f %.2f %.2f %.2f]];\n", head_normal.x, head_normal.y, head_normal.z, RAD_TO_DEG(B), RAD_TO_DEG(C));
                                //fprintf(f,"X=[X;[%.2f %.2f %.2f %.2f %.2f]];\n", Points[l].x, Points[l].y, Points[l].z - min_z + printer_char.min_axis.z + printer_char.plate_offset, RAD_TO_DEG(B), RAD_TO_DEG(C));

                            }
                        } else {
                            // we have a linear spline
                            pos = (int) Points.size();
                            // go for all points in spline
                            for (l = 1; l < pos; l++) {
                                // compute perpendicular to the normal
                                head_normal = Points[l] - Points[l - 1];
                                head_normal = head_normal.crossproduct(Normals[l - 1]);
                                head_normal = head_normal / head_normal.norm();
                                if (head_normal.z < 0)
                                    head_normal = v3(-head_normal.x, -head_normal.y, -head_normal.z);
                                
                                distance = sqrt(pow(lastx - Points[l].x, 2.0f) + pow(lasty - Points[l].y, 2.0f) + pow(lastz - Points[l].z, 2.0f));
                                total_extruded += distance * printer_char.extrusion_multiplier * (printer_char.extrusion_height * (printer_char.extrusion_width - printer_char.extrusion_height) + PI * pow(printer_char.extrusion_height / 2.0f, 2.0f)) /
                                        (PI * pow(printer_char.filament_radius, 2.0f) / 4.0f);

                                // compute spherical if 5 axis
                                if (GCode5D) {
                                    ANGLES(head_normal.x, head_normal.y, head_normal.z, &B, &C, &lastB, &lastC);
                                    carttransform(Points[l].x, Points[l].y, Points[l].z + z_offset, B, C, &x, &y, &z, &lastC, &lastB); // - min_z + printer_char.min_axis.z + printer_char.plate_offset
                                    //z -= printer_char.plate_offset;
                                } else {
                                    // A and C are already set to zero
                                    x = Points[l].x;
                                    y = Points[l].y;
                                    z = Points[l].z + z_offset; // - min_z + printer_char.min_axis.z 
                                }

                                lastx = Points[l].x;
                                lasty = Points[l].y;
                                lastz = Points[l].z;
                                exportG12command(f, printer_char.extrusion_speed, x, y, z, lastB, lastC, 0.0f, total_extruded);
                                //fprintf(f,"D=[D;[%.2f %.2f %.2f %.2f %.2f]];\n", head_normal.x, head_normal.y, head_normal.z, RAD_TO_DEG(B), RAD_TO_DEG(C));
                                //fprintf(f,"X=[X;[%.2f %.2f %.2f %.2f %.2f]];\n", Points[l].x, Points[l].y, Points[l].z- min_z + printer_char.min_axis.z + printer_char.plate_offset, RAD_TO_DEG(B), RAD_TO_DEG(C));
                            }
                        }
                    }
                }
            }
            lay++;
            lay %= printer_char.extrusion_number_of_layers;
        }
        lay++;
        lay %= printer_char.extrusion_number_of_layers;
    }

    fprintf(f, ";End Core\n");
    //fflush(f);

    counter=0;

    // Build CoreSpline
    vector<SimpleSpline> CoreSpline;
    x_end = buildSSpline(SLayers, CoreSpline, radius_offset);

    bool end_of_object = false;
    while (!end_of_object) {
        end_of_object=true; // we will finish next time, except if we found one G12 command
        
        snprintf(FileName,1024,"exterior_points%d.m",counter++);
        FILE *fexterior=fopen(FileName,"w");
        if(fexterior)
            fprintf(fexterior,"%%starting first layer\nPoints=[");
        
        GCodes.clear();
        
        vector<SimpleSpline> ESSpline, E2SSpline;
        
        // advance the offset and do as many loops as requested
        radius2 = 0.0f; // no G2 commands
        offset += radius_offset;
        expandSLayers(SLayers, ESLayers, offset);
        //expandSimpleSpline(CoreSpline, ESSpline, offset);
        offset += radius_offset;
        x_end = expandSLayers(SLayers, E2SLayers, offset);
        //expandSimpleSpline(CoreSpline, E2SSpline, offset);


//        //snprintf(FileName, 1024, "layers_3.m");
//        //exportSLayersMATLABFormat3D(ESLayers, FileName);
//        //snprintf(FileName, 1024, "layers_4.m");
//        //exportSLayersMATLABFormat3D(E2SLayers, FileName);
//
//        // we are in a new stage
//        lastB = 0;
//        lastC = 0;
//
//        fprintf(f, ";Starting first layer\n");
//
//        ss = ESSpline.size();
//        for (k = 0; k < ss; ++k) {
//            first_gcode = true;
//            vector<v3> Points = ESSpline[k].get_y_data();
//            vector<v3> Normals = ESSpline[k].get_normals();
//            vector<SPolygon*> SPoly = ESSpline[k].get_spoly();
//
//            pos = Points.size();
//            for (ps = 0; ps < pos; ps++) {
//                v3 point = Points[ps];
//                v3 normal = -Normals[ps];
//                SPolygon *thisSPoly = SPoly[ps];
//
//                if (GCode5D) {
//                    ANGLES(normal.x, normal.y, normal.z, &B, &C, &lastB, &lastC);
//                    carttransform(point.x, point.y, point.z + z_offset, B, C, &x, &y, &z, &lastC, &lastB); // - min_z + printer_char.min_axis.z + printer_char.plate_offset 
//                    //z -= printer_char.plate_offset;
//                } else {
//                    // A and C are already set to zero
//                    x = point.x;
//                    y = point.y;
//                    z = point.z + z_offset; // - min_z + printer_char.min_axis.z 
//                }
//
//
//                if (first_gcode) {
//                    first_gcode = false;
//                    lastx = point.x;
//                    lasty = point.y;
//                    lastz = point.z;
//                    //exportG0command(f, printer_char.speed, x, y, z, lastB, lastC);
//                    GCode tmp;
//                    tmp.type = 0;
//                    tmp.speed = printer_char.speed;
//                    tmp.x = x;
//                    tmp.y = y;
//                    tmp.z = z;
//                    tmp.B = lastB;
//                    tmp.C = lastC;
//                    tmp.radius2 = 0.0f;
//                    tmp.total_extruded = 0.0f;
//                    GCodes.push_back(tmp);
//                } else {
//                    distance = sqrt(pow(lastx - point.x, 2.0f) + pow(lasty - point.y, 2.0f) + pow(lastz - point.z, 2.0f));
//
//                    lastx = point.x;
//                    lasty = point.y;
//                    lastz = point.z;
//
//                    if (check_outer_spline(*thisSPoly, point)) {
//                        // go to point, but do not extrude
//                        //exportG0command(f, printer_char.speed, x, y, z, lastB, lastC);
//                        GCode tmp;
//                        tmp.type = 0;
//                        tmp.speed = printer_char.speed;
//                        tmp.x = x;
//                        tmp.y = y;
//                        tmp.z = z;
//                        tmp.B = lastB;
//                        tmp.C = lastC;
//                        tmp.radius2 = 0.0f;
//                        tmp.total_extruded = 0.0f;
//                        GCodes.push_back(tmp);
//                    } else {
//                        total_extruded += distance * printer_char.extrusion_multiplier * (printer_char.extrusion_height * (printer_char.extrusion_width - printer_char.extrusion_height) + PI * pow(printer_char.extrusion_height / 2.0f, 2.0f)) /
//                                (PI * pow(printer_char.filament_radius, 2.0f) / 4.0f);
//
//                        // go to point extruding
//                        //exportG12command(f, printer_char.extrusion_speed, x, y, z, lastB, lastC, radius2, total_extruded);
//                        GCode tmp;
//                        tmp.type = 1;
//                        tmp.speed = printer_char.extrusion_speed;
//                        tmp.x = x;
//                        tmp.y = y;
//                        tmp.z = z;
//                        tmp.B = lastB;
//                        tmp.C = lastC;
//                        tmp.radius2 = radius2;
//                        tmp.total_extruded = total_extruded;
//                        GCodes.push_back(tmp);
//
//                        if (fexterior){
//                            fprintf(fexterior, "%f,%f,%f;", point.x, point.y, point.z);
//                            fflush(fexterior);
//                        }
//
//                        end_of_object = false; // do not stop, since we have a inner point
//                    }
//                }
//                
//            
//            }
//            
//            
//            first_gcode = true;
//            Points = E2SSpline[k].get_y_data();
//            Normals = E2SSpline[k].get_normals();
//            SPoly = E2SSpline[k].get_spoly();
//
//            pos = Points.size();
//            for (ps = pos-1; ps > 0; ps--) {
//                v3 point = Points[ps];
//                v3 normal = -Normals[ps];
//                SPolygon *thisSPoly = SPoly[ps];
//
//                if (GCode5D) {
//                    ANGLES(normal.x, normal.y, normal.z, &B, &C, &lastB, &lastC);
//                    carttransform(point.x, point.y, point.z + z_offset, B, C, &x, &y, &z, &lastC, &lastB); // - min_z + printer_char.min_axis.z + printer_char.plate_offset 
//                    //z -= printer_char.plate_offset;
//                } else {
//                    // A and C are already set to zero
//                    x = point.x;
//                    y = point.y;
//                    z = point.z + z_offset; // - min_z + printer_char.min_axis.z 
//                }
//
//
//                if (first_gcode) {
//                    first_gcode = false;
//                    lastx = point.x;
//                    lasty = point.y;
//                    lastz = point.z;
//                    //exportG0command(f, printer_char.speed, x, y, z, lastB, lastC);
//                    GCode tmp;
//                    tmp.type = 0;
//                    tmp.speed = printer_char.speed;
//                    tmp.x = x;
//                    tmp.y = y;
//                    tmp.z = z;
//                    tmp.B = lastB;
//                    tmp.C = lastC;
//                    tmp.radius2 = 0.0f;
//                    tmp.total_extruded = 0.0f;
//                    GCodes.push_back(tmp);
//                } else {
//                    distance = sqrt(pow(lastx - point.x, 2.0f) + pow(lasty - point.y, 2.0f) + pow(lastz - point.z, 2.0f));
//
//                    lastx = point.x;
//                    lasty = point.y;
//                    lastz = point.z;
//
//                    if (check_outer_spline(*thisSPoly, point)) {
//                        // go to point, but do not extrude
//                        //exportG0command(f, printer_char.speed, x, y, z, lastB, lastC);
//                        GCode tmp;
//                        tmp.type = 0;
//                        tmp.speed = printer_char.speed;
//                        tmp.x = x;
//                        tmp.y = y;
//                        tmp.z = z;
//                        tmp.B = lastB;
//                        tmp.C = lastC;
//                        tmp.radius2 = 0.0f;
//                        tmp.total_extruded = 0.0f;
//                        GCodes.push_back(tmp);
//                    } else {
//                        total_extruded += distance * printer_char.extrusion_multiplier * (printer_char.extrusion_height * (printer_char.extrusion_width - printer_char.extrusion_height) + PI * pow(printer_char.extrusion_height / 2.0f, 2.0f)) /
//                                (PI * pow(printer_char.filament_radius, 2.0f) / 4.0f);
//
//                        // go to point extruding
//                        //exportG12command(f, printer_char.extrusion_speed, x, y, z, lastB, lastC, radius2, total_extruded);
//                        GCode tmp;
//                        tmp.type = 1;
//                        tmp.speed = printer_char.extrusion_speed;
//                        tmp.x = x;
//                        tmp.y = y;
//                        tmp.z = z;
//                        tmp.B = lastB;
//                        tmp.C = lastC;
//                        tmp.radius2 = radius2;
//                        tmp.total_extruded = total_extruded;
//                        GCodes.push_back(tmp);
//
//                        if (fexterior){
//                            fprintf(fexterior, "%f,%f,%f;", point.x, point.y, point.z);
//                            fflush(fexterior);
//                        }
//
//                        end_of_object = false; // do not stop, since we have a inner point
//                    }
//                }
//            }
//        }
//        
//        
        

        
        // Era aqui....
        
        
        // start at x_start (0.0f) and go to the longest spline
        for (t = x_start; t < x_end; t += radius_offset) {
            // go for all layers, from top to bottom (since we have finished at top)
            first_gcode = true;
            ls = (int) ESLayers.size();
            for (i = ls - 1; i >= 0; i--) {
                // for each Layer
                // get the polygon
                spoly = ESLayers[i].getSPolygon();
                // get the number of SPolygons
                ps = (int) spoly.size();
                // for all polygons
                for (j = 0; j < ps; ++j) {
                    if (spoly[j].getArea() < 0) {// we have the inner surface
                        splines = spoly[j].getSplines(); // splines to follow
                        
                        //splines_outer = spoly[(j+1)%ps].getSplines(); // external object splines
                        
                        ss = (int) splines.size();
                        // for all splines
                        for (k = 0; k < ss; ++k) {
                            if (t >= splines[k].get_front() && t <= splines[k].get_back()) {
                                // found it
                                vector<v3> Points = splines[k].get_y_data();
                                vector<v3> Normals = splines[k].get_normals();
                                pos = (int) Points.size();

                                v3 point = spline_eval(splines[k], t, &segment);
                                v3 normal = -Normals[segment - 1];


                                if (GCode5D) {
                                    ANGLES(normal.x, normal.y, normal.z, &B, &C, &lastB, &lastC);
                                    carttransform(point.x, point.y, point.z + z_offset, B, C, &x, &y, &z, &lastC, &lastB); // - min_z + printer_char.min_axis.z + printer_char.plate_offset 
                                    //z -= printer_char.plate_offset;
                                } else {
                                    // A and C are already set to zero
                                    x = point.x;
                                    y = point.y;
                                    z = point.z + z_offset; // - min_z + printer_char.min_axis.z 
                                }
                                

                                if (first_gcode) {
                                    first_gcode = false;
                                    lastx = point.x;
                                    lasty = point.y;
                                    lastz = point.z;
                                    //exportG0command(f, printer_char.speed, x, y, z, lastB, lastC);
                                    GCode tmp;
                                    tmp.type = 0;
                                    tmp.speed = printer_char.speed;
                                    tmp.x = x;
                                    tmp.y = y;
                                    tmp.z = z;
                                    tmp.B = lastB;
                                    tmp.C = lastC;
                                    tmp.radius2 = 0.0f;
                                    tmp.total_extruded = 0.0f;
                                    GCodes.push_back(tmp);
                                } else {
                                    distance = sqrt(pow(lastx - point.x, 2.0f) + pow(lasty - point.y, 2.0f) + pow(lastz - point.z, 2.0f));

                                    lastx = point.x;
                                    lasty = point.y;
                                    lastz = point.z;

                                    if(check_outer_spline(spoly[(j+1)%ps], point)){
                                        // go to point, but do not extrude
                                        //exportG0command(f, printer_char.speed, x, y, z, lastB, lastC);
                                        GCode tmp;
                                        tmp.type = 0;
                                        tmp.speed = printer_char.speed;
                                        tmp.x = x;
                                        tmp.y = y;
                                        tmp.z = z;
                                        tmp.B = lastB;
                                        tmp.C = lastC;
                                        tmp.radius2 = 0.0f;
                                        tmp.total_extruded = 0.0f;
                                        GCodes.push_back(tmp);
                                    } else {
                                        total_extruded += distance * printer_char.extrusion_multiplier * (printer_char.extrusion_height * (printer_char.extrusion_width - printer_char.extrusion_height) + PI * pow(printer_char.extrusion_height / 2.0f, 2.0f)) /
                                                (PI * pow(printer_char.filament_radius, 2.0f) / 4.0f);

                                        // go to point extruding
                                        //exportG12command(f, printer_char.extrusion_speed, x, y, z, lastB, lastC, radius2, total_extruded);
                                        GCode tmp;
                                        tmp.type = 1;
                                        tmp.speed = printer_char.extrusion_speed;
                                        tmp.x = x;
                                        tmp.y = y;
                                        tmp.z = z;
                                        tmp.B = lastB;
                                        tmp.C = lastC;
                                        tmp.radius2 = radius2;
                                        tmp.total_extruded = total_extruded;
                                        GCodes.push_back(tmp);
                                        
                                        if(fexterior){
                                            fprintf(fexterior,"%f,%f,%f;",point.x,point.y,point.z);
                                            fflush(fexterior);
                                        }
                                        
                                        end_of_object=false; // do not stop, since we have a inner point
                                    }
                                }

                                break;
                            }
                        }
                    }
                }
            }
            
            first_gcode = true;
            ls = (int) E2SLayers.size();
            for (i = 0; i < ls; i++) {
                // for each Layer
                // get the polygon
                spoly = E2SLayers[i].getSPolygon();
                // get the number of SPolygons
                ps = (int) spoly.size();
                // for all polygons
                for (j = 0; j < ps; ++j) {
                    if (spoly[j].getArea() < 0) {// we have the inner surface
                        splines = spoly[j].getSplines();
                        
                        //splines_outer = spoly[(j+1)%ps].getSplines(); // external object splines
                        
                        ss = (int) splines.size();
                        // for all splines
                        for (k = 0; k < ss; ++k) {
                            if (t >= splines[k].get_front() && t <= splines[k].get_back()) {
                                // found it
                                vector<v3> Points = splines[k].get_y_data();
                                vector<v3> Normals = splines[k].get_normals();
                                pos = (int) Points.size();

                                v3 point = spline_eval(splines[k], t, &segment);
                                v3 normal = -Normals[segment - 1];

                                if (GCode5D) {
                                    ANGLES(normal.x, normal.y, normal.z, &B, &C, &lastB, &lastC);
                                    carttransform(point.x, point.y, point.z + z_offset, B, C, &x, &y, &z, &lastC, &lastB); // - min_z + printer_char.min_axis.z + printer_char.plate_offset 
                                    //z -= printer_char.plate_offset;
                                } else {
                                    // A and C are already set to zero
                                    x = point.x;
                                    y = point.y;
                                    z = point.z + z_offset; // - min_z + printer_char.min_axis.z 
                                }
                                
                                if (first_gcode) {
                                    first_gcode = false;
                                    lastx = point.x;
                                    lasty = point.y;
                                    lastz = point.z;
                                    //exportG0command(f, printer_char.speed, x, y, z, lastB, lastC);
                                    GCode tmp;
                                    tmp.type = 0;
                                    tmp.speed = printer_char.speed;
                                    tmp.x = x;
                                    tmp.y = y;
                                    tmp.z = z;
                                    tmp.B = lastB;
                                    tmp.C = lastC;
                                    tmp.radius2 = radius2;
                                    tmp.total_extruded = total_extruded;
                                    GCodes.push_back(tmp);
                                    
                                } else {
                                    distance = sqrt(pow(lastx - point.x, 2.0f) + pow(lasty - point.y, 2.0f) + pow(lastz - point.z, 2.0f));
                                
                                    lastx = point.x;
                                    lasty = point.y;
                                    lastz = point.z;
                                    
                                    if(check_outer_spline(spoly[(j+1)%ps], point)){
                                        // go to point, but do not extrude
                                        //exportG0command(f, printer_char.speed, x, y, z, lastB, lastC);
                                        GCode tmp;
                                        tmp.type = 0;
                                        tmp.speed = printer_char.speed;
                                        tmp.x = x;
                                        tmp.y = y;
                                        tmp.z = z;
                                        tmp.B = lastB;
                                        tmp.C = lastC;
                                        tmp.radius2 = radius2;
                                        tmp.total_extruded = total_extruded;
                                        GCodes.push_back(tmp);
                                        
                                    } else {
                                        total_extruded += distance * printer_char.extrusion_multiplier * (printer_char.extrusion_height * (printer_char.extrusion_width - printer_char.extrusion_height) + PI * pow(printer_char.extrusion_height / 2.0f, 2.0f)) /
                                                (PI * pow(printer_char.filament_radius, 2.0f) / 4.0f);

                                        // go to point extruding
                                        //exportG12command(f, printer_char.extrusion_speed, x, y, z, lastB, lastC, radius2, total_extruded);
                                        GCode tmp;
                                        tmp.type = 1;
                                        tmp.speed = printer_char.extrusion_speed;
                                        tmp.x = x;
                                        tmp.y = y;
                                        tmp.z = z;
                                        tmp.B = lastB;
                                        tmp.C = lastC;
                                        tmp.radius2 = radius2;
                                        tmp.total_extruded = total_extruded;
                                        GCodes.push_back(tmp);
                                        
                                        if(fexterior){
                                            fprintf(fexterior,"%f,%f,%f;",point.x,point.y,point.z);
                                            fflush(fexterior);
                                        }
                                        
                                        end_of_object=false;
                                    }
                                }

                                break;
                            }
                        }
                    }
                }

            }
        }
        
        // output G codes
        {int k;
        for(k=0;k<GCodes.size()-1;k++){
            if(GCodes[k].type==0 && GCodes[k+1].type==0)
                continue;
            if(GCodes[k].type==0)
                exportG0command(f, GCodes[k].speed, GCodes[k].x, GCodes[k].y, GCodes[k].z, GCodes[k].B, GCodes[k].C);
            else
                exportG12command(f, GCodes[k].speed, GCodes[k].x, GCodes[k].y, GCodes[k].z, GCodes[k].B, GCodes[k].C, GCodes[k].radius2, GCodes[k].total_extruded);
        }
        if (GCodes[k].type != 0)
            exportG12command(f, GCodes[k].speed, GCodes[k].x, GCodes[k].y, GCodes[k].z, GCodes[k].B, GCodes[k].C, GCodes[k].radius2, GCodes[k].total_extruded);
        
        // no need to export last G0 comand
        // else exportG0command(f, GCodes[k].speed, GCodes[k].x, GCodes[k].y, GCodes[k].z, GCodes[k].B, GCodes[k].C);
        }
        
        fprintf(f, ";End first layer\n");
        
        if(fexterior){
            fprintf(fexterior,"];\nplot3(Points(:,1),Points(:,2),Points(:,3),'.k')\n");
            fclose(fexterior);
        }
        
        // don't go to another layer if object ended
        if(end_of_object)
            break;
        
        GCodes.clear();
        
        fprintf(f, ";Starting second layer\n");
        
        snprintf(FileName,1024,"exterior_points%d.m",counter++);
        fexterior=fopen(FileName,"w");
        if(fexterior)
            fprintf(fexterior,"%%starting second layer\nPoints=[");
        
        // build external layers using normal vectors
        // go for all layers
        ls = (int) ESLayers.size();
        for (i = 0; i < ls; i++) {
            // for each Layer
            // get the polygon
            lay = 0; // always begin with inner layer
            for (lays = 0; lays < printer_char.extrusion_number_of_layers; lays++) {
                //printf("Lay=%d\n", lay);
                if (lay) spoly = E2SLayers[i].getSPolygon(); // E2
                else spoly = ESLayers[i].getSPolygon(); // E

                // get the number of SPolygons
                ps = (int) spoly.size();
                // go for all polygons
                for (j = 0; j < ps; ++j) {
                    if (spoly[j].getArea() < 0) { // we have an inner loop, so start with it
                        vector<v3> Points;
                        vector<v3> Normals;
                        vector<float> x_data;
                        v3 head_normal;

                        splines = spoly[j].getSplines(); // get splines
                        
                        //splines_outer = spoly[(j+1)%ps].getSplines(); // external object splines
                                                
                        ss = (int) splines.size();

                        // get nozzel to beginning of spline with a G00 command
                        Points = splines[0].get_y_data();
                        Normals = splines[0].get_normals();
                        x_data = splines[0].get_x_data();
                        head_normal = -Normals[0];

                        if (GCode5D) {
                            // compute B and C table angles
                            ANGLES(head_normal.x, head_normal.y, head_normal.z, &B, &C, &lastB, &lastC);
                            // tranform x,y,z coordinates to new axis (rotation at origin)
                            carttransform(Points[0].x, Points[0].y, Points[0].z + z_offset, B, C, &x, &y, &z, &lastC, &lastB); // - min_z + printer_char.min_axis.z + printer_char.plate_offset 
                            // and change plate offset
                            //z -= printer_char.plate_offset;
                        } else {
                            x = Points[0].x;
                            y = Points[0].y;
                            z = Points[0].z + z_offset; // - min_z + printer_char.min_axis.z 
                        }

                        // lastB and lastC were updated in carttransform to be the actual B and C
                        //exportG0command(f, printer_char.speed, x, y, z, lastB, lastC);
                        GCode tmp;
                        tmp.type = 0;
                        tmp.speed = printer_char.speed;
                        tmp.x = x;
                        tmp.y = y;
                        tmp.z = z;
                        tmp.B = lastB;
                        tmp.C = lastC;
                        tmp.radius2 = radius2;
                        tmp.total_extruded = total_extruded;
                        GCodes.push_back(tmp);

                        //fprintf(f,"D=[D;[%.2f %.2f %.2f %.2f %.2f]];\n", head_normal.x, head_normal.y, head_normal.z, RAD_TO_DEG(B), RAD_TO_DEG(C));
                        //fprintf(f,"X=[X;[%.2f %.2f %.2f %.2f %.2f]];\n", Points[0].x, Points[0].y, Points[0].z - min_z + printer_char.min_axis.z + printer_char.plate_offset, RAD_TO_DEG(B), RAD_TO_DEG(C));

                        // last point before transform
                        lastx = Points[0].x;
                        lasty = Points[0].y;
                        lastz = Points[0].z;

                        // go for all spline segments
                        for (k = 0; k < ss; ++k) {
                            Points = splines[k].get_y_data();
                            Normals = splines[k].get_normals();
                            x_data = splines[k].get_x_data();

                            pos = (int) Points.size();
                            if (pos < 2) continue; // a segment with just one point, should not happen

                            // we have a linear spline
                            pos = (int) Points.size();
                            // go for all points in spline
                            for (l = 1; l < pos; l++) {

                                head_normal=-Normals[l];

                                // compute spherical if 5 axis
                                if (GCode5D) {
                                    ANGLES(head_normal.x, head_normal.y, head_normal.z, &B, &C, &lastB, &lastC);
                                    carttransform(Points[l].x, Points[l].y, Points[l].z + z_offset, B, C, &x, &y, &z, &lastC, &lastB); // - min_z + printer_char.min_axis.z + printer_char.plate_offset 
                                    //z -= printer_char.plate_offset;
                                } else {
                                    // A and C are already set to zero
                                    x = Points[l].x;
                                    y = Points[l].y;
                                    z = Points[l].z  + z_offset; //- min_z + printer_char.min_axis.z
                                }
                                
                                distance = sqrt(pow(lastx - Points[l].x, 2.0f) + pow(lasty - Points[l].y, 2.0f) + pow(lastz - Points[l].z, 2.0f));
                                
                                lastx = Points[l].x;
                                lasty = Points[l].y;
                                lastz = Points[l].z;
                                
                                if(check_outer_spline(spoly[(j+1)%ps], Points[l])){
                                    // go to point, but do not extrude
                                    //exportG0command(f, printer_char.speed, x, y, z, lastB, lastC);
                                    GCode tmp;
                                    tmp.type = 0;
                                    tmp.speed = printer_char.speed;
                                    tmp.x = x;
                                    tmp.y = y;
                                    tmp.z = z;
                                    tmp.B = lastB;
                                    tmp.C = lastC;
                                    tmp.radius2 = radius2;
                                    tmp.total_extruded = total_extruded;
                                    GCodes.push_back(tmp);
                                    
                                } else {
                                    total_extruded += distance * printer_char.extrusion_multiplier * (printer_char.extrusion_height * (printer_char.extrusion_width - printer_char.extrusion_height) + PI * pow(printer_char.extrusion_height / 2.0f, 2.0f)) /
                                        (PI * pow(printer_char.filament_radius, 2.0f) / 4.0f);
                                
                                    // go to point extruding
                                    //exportG12command(f, printer_char.extrusion_speed, x, y, z, lastB, lastC, 0.0f, total_extruded);
                                    GCode tmp;
                                    tmp.type = 1;
                                    tmp.speed = printer_char.extrusion_speed;
                                    tmp.x = x;
                                    tmp.y = y;
                                    tmp.z = z;
                                    tmp.B = lastB;
                                    tmp.C = lastC;
                                    tmp.radius2 = radius2;
                                    tmp.total_extruded = total_extruded;
                                    GCodes.push_back(tmp);
                                    
                                    if(fexterior)
                                        fprintf(fexterior,"%f,%f,%f;",Points[l].x,Points[l].y,Points[l].z);
                                    
                                    end_of_object=false;
                                }
                                //fprintf(f,"D=[D;[%.2f %.2f %.2f %.2f %.2f]];\n", head_normal.x, head_normal.y, head_normal.z, RAD_TO_DEG(B), RAD_TO_DEG(C));
                                //fprintf(f,"X=[X;[%.2f %.2f %.2f %.2f %.2f]];\n", Points[l].x, Points[l].y, Points[l].z- min_z + printer_char.min_axis.z + printer_char.plate_offset, RAD_TO_DEG(B), RAD_TO_DEG(C));
                            }
                        }
                    }
                }
                lay++;
                lay %= printer_char.extrusion_number_of_layers;
            }
        }

        // output G codes
        {int k;
        for(k=0;k<GCodes.size()-1;k++){
            if(GCodes[k].type==0 && GCodes[k+1].type==0)
                continue;
            if(GCodes[k].type==0)
                exportG0command(f, GCodes[k].speed, GCodes[k].x, GCodes[k].y, GCodes[k].z, GCodes[k].B, GCodes[k].C);
            else
                exportG12command(f, GCodes[k].speed, GCodes[k].x, GCodes[k].y, GCodes[k].z, GCodes[k].B, GCodes[k].C, GCodes[k].radius2, GCodes[k].total_extruded);
        }
        if (GCodes[k].type == 0)
            exportG12command(f, GCodes[k].speed, GCodes[k].x, GCodes[k].y, GCodes[k].z, GCodes[k].B, GCodes[k].C, GCodes[k].radius2, GCodes[k].total_extruded);
        
        // no need to export last G0 comand
        // else exportG0command(f, GCodes[k].speed, GCodes[k].x, GCodes[k].y, GCodes[k].z, GCodes[k].B, GCodes[k].C);
        }
        
        fprintf(f, ";End second layer\n");
        
        if(fexterior){
            fprintf(fexterior,"];\nplot3(Points(:,1),Points(:,2),Points(:,3),'.k')\n");
            fclose(fexterior);
        }
        
    }
    
    exportEpilogue(f);
    fclose(f);
    return 0;
}

// output a single MATLAB file with the slices, layed out in 3D form

int exportLayersGCode3Axes(
        std::vector<Layer> &Layers,
        char *filename) {

    FILE *f = NULL;
    int i, j, l, ls, ps, pos;
    vector<Polygon> poly;
    vector<Spline> splines;
    vector<v3> Points;
    char command[3] = "G1";


    fopen_s(&f, filename, "w");
    if (!f) return 1;

    fprintf(f, "G0 F1000\n");
    ls = (int) Layers.size();
    for (i = 0; i < ls; ++i) {
        // for each Layer
        // get the polygon
        poly = Layers[i].getPoly();
        // get the number of SPolygons
        ps = (int) poly.size();
        for (j = 0; j < ps; ++j) {
            if (poly[j].getArea() < 0) { // we have a inner loop, so start with it
                Points = poly[j].getPoints();
                pos = (int) Points.size();
                for (l = 0; l < pos; l++) {
                    fprintf(f, "%s X%.2f Y%.2f Z%.2f\n", command, Points[l].x, Points[l].y, Points[l].z);
                }
            }
        }
    }

    fclose(f);
    return 0;
}
