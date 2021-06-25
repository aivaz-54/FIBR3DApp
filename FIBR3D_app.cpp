#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>

#include "FIBR3Dapp.hpp"
#include "clipper_utils.hpp"
#include "pswarm.h"

#include "Inspection.hpp"
#include "PartsOrientation.hpp"
#include "FillingAlgorithms.hpp"

#include "MMDS.hpp"

extern "C" {
    int PSwarm(int n, void(*objf)(int n, int m, double *x, double *fx), double *llb, double *uub, int lincons, double *AA,
            double *b, double **sol, double *f, double *xx);
    extern struct Options opt;
}

// list of printing sequences
extern vector<vector<vector<pair<int, int>>>> L_f;

// a pointer to the STL being optimized
extern stl_file *stl_optimize;
extern Mesh complete_mesh;

// our printer characteristics
printer printer_char;
extern inspection inspection_char;
extern void printSTL_Scad_no(char *, stl_file *stl);
extern void ANGLES(float vx, float vy, float vz, float *B, float *C, float *lastb, float *lastc);
extern void my_rotate(float *x, float *y, float angle);
extern void my_stl_rotate(stl_file *stl, float B, float C);

// print level - 0 no printing - 1 printing
int print_level = 1;


// minimum of two values

const float mymin(const float &a, const float &b) {
    return (b > a) ? a : b;
}

// maximum of two values

const float mymax(const float &a, const float &b) {
    return (b < a) ? a : b;
}

// sum and subtraction of points

v3 operator-(const v3 &a, const v3 &b) {
    return v3(a.x - b.x, a.y - b.y, a.z - b.z);
}

v3 operator+(const v3 &a, const v3 &b) {
    return v3(a.x + b.x, a.y + b.y, a.z + b.z);
}


// intersect facet with horizontal plane
// given a facet return the line segment that intersect the z plane
// return 0 on success (intersection exists)

int intersectHorizontalPlane(stl_facet facet, const float &z, LineSegment &ls) {
    std::vector<v3> intersectPoints;
    std::vector<std::vector<v3>::size_type> points_on_layer;
    int start = 0, i;
    v3 a, b, bMinusa, c, d;
    float da, db, dc;

    // a triangle has 3 vertices that construct 3 line segments

    float minz = mymin(facet.vertex[0].z, mymin(facet.vertex[1].z, facet.vertex[2].z));
    float maxz = mymax(facet.vertex[0].z, mymax(facet.vertex[1].z, facet.vertex[2].z));

    if (minz == maxz) // we have a top or bottom facet
        return -1;

    if (maxz < z) // three points below plane, no intersection
        return -1;

    if (minz > z) // three points above plane, no intersection
        return 1;

    // we have an intersection
    if (facet.vertex[0].z == minz)
        start = 0;
    else if (facet.vertex[1].z == minz)
        start = 1;
    else if (facet.vertex[2].z == minz)
        start = 2;

    // used field is used to acknowledge segments on the plane
    ls.plane = false;

    for (i = start; (i - start) < 3; ++i) {
        a = v3(facet.vertex[i % 3].x, facet.vertex[i % 3].y, facet.vertex[i % 3].z);
        b = v3(facet.vertex[(i + 1) % 3].x, facet.vertex[(i + 1) % 3].y, facet.vertex[(i + 1) % 3].z);
        da = a.z - z;
        db = b.z - z;
        if (da == 0 && db == 0) { // edge is on plane
            c = v3(facet.vertex[(i + 2) % 3].x, facet.vertex[(i + 2) % 3].y, facet.vertex[(i + 2) % 3].z);
            dc = c.z - z;
            if (dc == 0) { //facet is on plane
                if (facet.normal.z < 0)
                    swap(a, b);

            } else {
                if (dc > 0)
                    swap(a, b);
            }

            intersectPoints.clear();
            points_on_layer.clear();
            intersectPoints.push_back(a);
            intersectPoints.push_back(b);
            // no point to continue
            // acknowledge linesegment on plane
            ls.plane = true;
            break;

        } else {
            if (da == 0) { // plane falls exactly on one of the three Triangle vertices
                //if (intersectPoints.size() < 2)
                intersectPoints.push_back(a);
                points_on_layer.push_back(intersectPoints.size() - 1);
            } else {
                if (db == 0) { // plane falls exactly on one of the three Triangle vertices
                    //if (intersectPoints.size() < 2)
                    intersectPoints.push_back(b);
                    points_on_layer.push_back(intersectPoints.size() - 1);
                } else {
                    if (da * db < 0) {
                        const float s = da / (da - db); // intersection factor (between 0 and 1)
                        bMinusa = b - a;
                        intersectPoints.push_back(a + bMinusa * s);
                    }
                }
            }
        }
    }

    if (!points_on_layer.empty()) {
        // remove duplicated vertices on layers
        intersectPoints.erase(intersectPoints.begin() + points_on_layer[1]);
    }


    if (intersectPoints.size() == 2) {
        // Output the intersecting line segment object
        ls.v[0] = intersectPoints[0];
        ls.v[1] = intersectPoints[1];
        ls.normal = v3(facet.normal.x, facet.normal.y, facet.normal.z);
        ls.rounding();
        ls.used = false;

        d = ls.v[1] - ls.v[0];
        // segment has size zero after rounding?
        if (d.norm() > 0) {
            return 0;
        }
    }

    // one or three point do not form a line segment!
    return -2;
}



// Take an input STL and fill the output parameter �slicesWithLineSegments� with line segments for
// each slice

void triMeshSlicer(stl_file *stl,
        std::vector<std::vector<LineSegment>> &slicesWithLineSegments, // the result slices
        const float sliceSize) { // slice size in 3D Model digital units
    int i, t;
    LineSegment ls, tmp;
    bool found;

    float z;

    z = stl->stats.min.z + printer_char.extrusion_height / 2.0f; // find the minimal z coordinate of the model (z0)

    for (; z <= stl->stats.max.z; z += sliceSize) { // position the plane according to slice index
        std::vector<LineSegment> linesegs; // the linesegs vector for each slice
        for (t = 0; t < stl->stats.number_of_facets; ++t) { // iterate all mesh triangles
            if (!intersectHorizontalPlane(stl->facet_start[t], z, ls)) {// the plane intersects the triangle?
                if (ls.plane) {
                    // we may only find coincident edges in the same layer
                    found = false;
                    for (i = 0; i < (int) linesegs.size(); i++) {
                        tmp = linesegs[i];
                        if (!tmp.plane)continue;
                        if (ls.equal(tmp)) { // found a duplicated one
                            found = true;
                            break;
                        }
                    }
                    if (!found)
                        linesegs.push_back(ls); // push a new Line Segment object to this slice
                } else {
                    linesegs.push_back(ls); // push a new Line Segment object to this slice
                }
            }
        }
        if (linesegs.size() > 0)
            slicesWithLineSegments.push_back(linesegs); // push this vector to the slices vector
    }
}


// build layers from a vector of line segments
// each layer is formed by a set of polygons

void buildLayers(
        std::vector<std::vector<LineSegment>> &slicesWithLineSegments, // given a vetor of slice line segments
        std::vector<Layer> &Layers) { // return layers
    int nSlices = (int) slicesWithLineSegments.size();
    int i, k, l, closest, lss_s;
    v3 tmp, last, first, lastNormal, firstNormal; // last position of the head
    float dist, closestdist;
    bool closestv0;

    last = v3(0.0f, 0.0f, 0.0f);

    for (i = 0; i < nSlices; ++i) {
        if (slicesWithLineSegments[i].size() <= 0) continue;

        std::vector<LineSegment> &lss = slicesWithLineSegments[i];
        Layer lay;
        Polygon poly;
        bool allnotused = true;
        lss_s = (int) lss.size();


        // first polygon line. find the closest point to start
        closest = -1;
        closestdist = 1.0e10f;
        for (l = 0; l < lss_s; ++l) {
            if (lss[l].used)continue;

            tmp = lss[l].v[0] - last;
            dist = tmp.norm();
            if (dist < closestdist) {
                closestdist = dist;
                closest = l;
            }
        }

        if (closest < 0) continue; // nothing to do in this layer?

        poly.push_back(lss[closest].v[0], lss[closest].normal);
        poly.push_back(lss[closest].v[1], lss[closest].normal);
        first = lss[closest].v[0];
        firstNormal = lss[closest].normal;
        last = lss[closest].v[1];
        lastNormal = lss[closest].normal;
        lss[closest].used = true;

        while (allnotused) {
            for (k = 0; k < lss_s; ++k) {
                if (lss[k].used)
                    continue;

                if (lss[k].v[0].x == last.x && lss[k].v[0].y == last.y) {// && lss[k].v[0].z == last.z) {
                    // add the other end of line
                    poly.correct_normal(lss[k].normal);
                    poly.push_back(lss[k].v[1], lss[k].normal);
                    last = lss[k].v[1];
                    lastNormal = lss[k].normal;
                    lss[k].used = true;
                    // we need to restart;
                    break;
                }
                if (lss[k].v[1].x == last.x && lss[k].v[1].y == last.y) {// && lss[k].v[1].z == last.z) {
                    // add the other end of line
                    poly.correct_normal(lss[k].normal);
                    poly.push_back(lss[k].v[0], lss[k].normal);
                    last = lss[k].v[0];
                    lastNormal = lss[k].normal;
                    lss[k].used = true;
                    // we need to restart;
                    break;
                }
            }
            // check for end of polygon
            if (k >= lss_s) {
                // we found no continuation to the polygon, so end this one or insert a continuation segment
                if (last.x == first.x && last.y == first.y) { // && last.z == first.z) { // we have a closed polygon
                    poly.correct_last_normal();
                    poly.computeArea();
                    if (fabs(poly.getArea()) > SLICE_EPS)
                        lay.push_back(poly);

                    // start a new polygon
                    closest = -1;
                    closestdist = 1.0e10f;
                    for (l = 0; l < lss_s; ++l) {
                        if (lss[l].used)
                            continue;

                        tmp = lss[l].v[0] - last;
                        dist = tmp.norm();
                        if (dist < closestdist) {
                            closestdist = dist;
                            closest = l;
                        }
                    }

                    if (closest >= 0) {
                        // start it
                        poly.clean();
                        poly.push_back(lss[closest].v[0], lss[closest].normal);
                        poly.push_back(lss[closest].v[1], lss[closest].normal);
                        first = lss[closest].v[0];
                        firstNormal = lss[closest].normal;
                        last = lss[closest].v[1];
                        lastNormal = lss[closest].normal;
                        lss[closest].used = true;
                    } else {
                        // no available point to start polygon
                        allnotused = false;
                    }
                } else {
                    // find the closest segment to the end point
                    closestv0 = true;
                    closestdist = 1.0e10f;
                    closest = -1;
                    for (l = 0; l < lss_s; ++l) {
                        if (lss[l].used)continue;

                        tmp = lss[l].v[0] - last;
                        dist = tmp.norm();
                        if (dist < closestdist) { // && lss[l].normal.dotproduct(lastNormal) > 0) {
                            closestdist = dist;
                            closest = l;
                            closestv0 = true;
                        }
                        tmp = lss[l].v[1] - last;
                        dist = tmp.norm();
                        if (dist < closestdist) { // && lss[l].normal.dotproduct(lastNormal) > 0) {
                            closestdist = dist;
                            closest = l;
                            closestv0 = false;
                        }
                    }
                    // the first was used, but may allow to close the polygon, if it is closer than the closest
                    tmp = first - last;
                    if (closest < 0 || tmp.norm() <= closestdist) { // && lastNormal.dotproduct(firstNormal)>0)) {
                        //printf("Open polygon!\n");
                        // close the polygon
                        poly.push_back(first, firstNormal);
                        // poly.correct_last_normal(); // no need to correct normal, since first is last
                        // and save it
                        poly.computeArea();
                        if (fabs(poly.getArea()) > SLICE_EPS)
                            lay.push_back(poly);

                        // start a new polygon
                        closest = -1;
                        closestdist = 1.0e10f;
                        // start a new polygon
                        for (l = 0; l < lss_s; ++l) {
                            if (lss[l].used)continue;

                            tmp = lss[l].v[0] - last;
                            dist = tmp.norm();
                            if (dist < closestdist) {
                                closestdist = dist;
                                closest = l;
                            }
                        }

                        if (closest >= 0) {
                            // start it
                            poly.clean();
                            poly.push_back(lss[closest].v[0], lss[closest].normal);
                            poly.push_back(lss[closest].v[1], lss[closest].normal);
                            first = lss[closest].v[0];
                            firstNormal = lss[closest].normal;
                            last = lss[closest].v[1];
                            lastNormal = lss[closest].normal;
                            lss[closest].used = true;
                        } else {
                            // no available point to start polygon
                            allnotused = false;
                        }
                    } else {
                        // add the closest point to the polygon
                        if (closestv0) {// add v0 and v1
                            poly.push_back(lss[closest].v[0], lss[closest].normal);
                            poly.push_back(lss[closest].v[1], lss[closest].normal);
                            last = lss[closest].v[1];
                            lastNormal = lss[closest].normal;
                            lss[closest].used = true;
                        } else {// add v1 and v0
                            poly.push_back(lss[closest].v[1], lss[closest].normal);
                            poly.push_back(lss[closest].v[0], lss[closest].normal);
                            last = lss[closest].v[0];
                            lastNormal = lss[closest].normal;
                            lss[closest].used = true;
                        }
                    }

                }
            }
        }
        if (lay.size() > 0)
            Layers.push_back(lay);
    }

}

// build splines layers from layers
void buildSplineLayers(std::vector<Layer> &Layers, std::vector<SLayer> &SLayers) {
    int i, j, ls, ps;
    vector<Polygon> poly;
    vector<v3> points, normals;
    
    ls = (int) Layers.size();
    for (i = 0; i < ls; ++i) {
        poly = Layers[i].getPoly();
        ps = (int) poly.size();
        SLayer slayer;
        for (j = 0; j < ps; ++j) {
            SPolygon spoly;
            points = poly[j].getPoints();
            normals = poly[j].getNormals();
            buildSpline(points, normals, spoly);
            computeSpline(spoly);
            slayer.push_back(spoly);
        }
        SLayers.push_back(slayer);
    }
}

// check if the part is a shell, i.e. formed by inner and outer polygons

bool is_shell(std::vector<SLayer> &SLayers) {
    int i, ls, ps;
    vector<SPolygon> spoly;

    ls = (int) SLayers.size();
    for (i = 0; i < ls; ++i) {
        spoly = SLayers[i].getSPolygon();
        ps = (int) spoly.size();
        if (ps != 2) // we do not have 2 polygons and therefore we do not have a shell part
            return false;
        if (spoly[0].getArea() * spoly[1].getArea() > 0) // both negative or both positive
            return false;
    }
    return true;
}


//// Slice a given STL using a rotation along the x and y
//// the given STL struct is not changed
//
//void Slice(stl_file *stl, float x, float y, char *filename) {
//    char FileName[1024];
//    float sliceSize = printer_char.extrusion_height;
//    stl_file stl_local;
//    bool shell;
//
//    // initialize local stl structure
//    stl_initialize(&stl_local);
//
//    // make a copy, so rotate can be undone
//    stl_local.stats = stl->stats;
//    // allocate memory for facets
//    stl_local.facet_start = (stl_facet*) malloc(stl_local.stats.number_of_facets * sizeof (stl_facet));
//    if (stl_local.facet_start == NULL)
//        perror("Slice");
//    // copy facets
//    memcpy(stl_local.facet_start, stl->facet_start, stl->stats.number_of_facets * sizeof (stl_facet));
//
//    // rotate along x and y axis
//    //stl_rotate_xy(&stl_local, x, y);
//    stl_rotate_x(&stl_local, x);
//    stl_rotate_y(&stl_local, y);
//
//    // Generate Slices
//    std::vector<std::vector < LineSegment>> slicesWithLineSegments;
//    triMeshSlicer(&stl_local, slicesWithLineSegments, sliceSize);
//
//    // export slicing to MATLAB
//    snprintf(FileName, 1024, "%s_out.m", filename);
//    exportSingleMATLABFormat3D(slicesWithLineSegments, FileName);
//
//    // Build layers with polygons
//    std::vector<Layer> Layers;
//    buildLayers(slicesWithLineSegments, Layers);
//
//    // export Layers to MATLAB
//    snprintf(FileName, 1024, "%s_layers.m", filename);
//    exportLayersMATLABFormat3D(Layers, FileName);
//
//    // Build spline layers with polygons
//    std::vector<SLayer> SLayers;
//    buildSplineLayers(Layers, SLayers);
//
//    // export SLayers to MATLAB
//    snprintf(FileName, 1024, "%s_slayers.m", filename);
//    //exportSLayersMATLABFormat3D(SLayers, FileName);
//
//
//    shell = is_shell(SLayers);
//
//
//    // export Spline Layers to MATLAB
//    // sprintf_s(FileName, 1024, "%s_splines.m", filename);
//    // exportSLayersMATLABFormat3D(SLayers, FileName);
//
//    //exportLayersGCode3Axes(Layers, "Poly.gcode");
//
//    //exportSLayersGCode3Axes(SLayers, "Splines.gcode");
//
//    //export G Code from Spline Layers
//    if (shell) {
//        snprintf(FileName, 1024, "%s_splines.gcode", filename);
//        exportSLayersGCode5Axes_shell(SLayers, stl_local.stats.min.z, max(stl_local.stats.size.x, stl_local.stats.size.y), true, FileName);
//    } else {
//        printf("Part %s is not a shell type part\n", filename);
//    }
//
//    free(stl_local.facet_start);
//}

// compute the are of a given facet (area of a triangle given by three points)

float facets_area(stl_facet facet) {

    return (float) 0.5 * sqrt(powf((facet.vertex[1].x * facet.vertex[0].y)
            - (facet.vertex[2].x * facet.vertex[0].y)
            - (facet.vertex[0].x * facet.vertex[1].y)
            + (facet.vertex[2].x * facet.vertex[1].y)
            + (facet.vertex[0].x * facet.vertex[2].y)
            - (facet.vertex[1].x * facet.vertex[2].y), 2.0f)
            + powf((facet.vertex[1].x * facet.vertex[0].z)
            - (facet.vertex[2].x * facet.vertex[0].z)
            - (facet.vertex[0].x * facet.vertex[1].z)
            + (facet.vertex[2].x * facet.vertex[1].z)
            + (facet.vertex[0].x * facet.vertex[2].z)
            - (facet.vertex[1].x * facet.vertex[2].z), 2.0f)
            + powf((facet.vertex[1].y * facet.vertex[0].z)
            - (facet.vertex[2].y * facet.vertex[0].z)
            - (facet.vertex[0].y * facet.vertex[1].z)
            + (facet.vertex[2].y * facet.vertex[1].z)
            + (facet.vertex[0].y * facet.vertex[2].z)
            - (facet.vertex[1].y * facet.vertex[2].z), 2.0f));

}

// check if facet in on floor (one vertex on floor)

bool floor_facet(stl_facet f1, float z) {
    int i;

    for (i = 0; i < 3; i++)
        if (f1.vertex[i].z < z + SLICE_EPS)
            return true;

    return false;
}

void print_facet(stl_facet f1) {
    int i;

    for (i = 0; i < 3; i++)
        printf("line([%f,%f],[%f,%f],[%f,%f]);\n", f1.vertex[i].x, f1.vertex[(i + 1) % 3].x,
            f1.vertex[i].y, f1.vertex[(i + 1) % 3].y, f1.vertex[i].z, f1.vertex[(i + 1) % 3].z);

}


// check if facet f1 is joint with facet f2

bool joint_facets(stl_facet f1, stl_facet f2) {
    int i;
    float d, total, norm, a1, a2, a3;
    float minx1, minx2, miny1, miny2, minz1, minz2;
    float maxx1, maxx2, maxy1, maxy2, maxz1, maxz2;
    v3 pa1, pa2, pa3;


    // Test bounding boxes, if they don't intersect then we're finished.
    minx1 = mymin(f1.vertex[0].x, mymin(f1.vertex[1].x, f1.vertex[2].x));
    maxx2 = mymax(f2.vertex[0].x, mymax(f2.vertex[1].x, f2.vertex[2].x));

    if (maxx2 < minx1 - SLICE_EPS) // They are apart in the x axes
        return false;

    maxx1 = mymax(f1.vertex[0].x, mymax(f1.vertex[1].x, f1.vertex[2].x));
    minx2 = mymin(f2.vertex[0].x, mymin(f2.vertex[1].x, f2.vertex[2].x));

    if (maxx1 < minx2 - SLICE_EPS) // They are apart in the x axes
        return false;

    miny1 = mymin(f1.vertex[0].y, mymin(f1.vertex[1].y, f1.vertex[2].y));
    maxy2 = mymax(f2.vertex[0].y, mymax(f2.vertex[1].y, f2.vertex[2].y));

    if (maxy2 < miny1 - SLICE_EPS) // They are apart in the y axes
        return false;

    maxy1 = mymax(f1.vertex[0].y, mymax(f1.vertex[1].y, f1.vertex[2].y));
    miny2 = mymin(f2.vertex[0].y, mymin(f2.vertex[1].y, f2.vertex[2].y));

    if (maxy1 < miny2 - SLICE_EPS) // They are apart in the y axes
        return false;

    minz1 = mymin(f1.vertex[0].z, mymin(f1.vertex[1].z, f1.vertex[2].z));
    maxz2 = mymax(f2.vertex[0].z, mymax(f2.vertex[1].z, f2.vertex[2].z));

    if (maxz2 < minz1 - SLICE_EPS) // They are apart in the z axes
        return false;

    maxz1 = mymax(f1.vertex[0].z, mymax(f1.vertex[1].z, f1.vertex[2].z));
    minz2 = mymin(f2.vertex[0].z, mymin(f2.vertex[1].z, f2.vertex[2].z));

    if (maxz1 < minz2 - SLICE_EPS) // They are apart in the z axes
        return false;

    //printf("Passei\n");
    // test whether all vertices of first polygon are on the same side as the plane described by the second, if they are then finish.

    // first facet plane is d=-(nx*x+ny*y+nz*z)
    d = -f1.normal.x * f1.vertex[0].x - f1.normal.y * f1.vertex[0].y - f1.normal.z * f1.vertex[0].z;

    // check if at least one f2 vertex is on the plane. We are looking for contiguous facets, so they cannot be intersection facets.
    if (fabs(f1.normal.x * f2.vertex[0].x + f1.normal.y * f2.vertex[0].y + f1.normal.z * f2.vertex[0].z + d) > SLICE_EPS &&
            fabs(f1.normal.x * f2.vertex[1].x + f1.normal.y * f2.vertex[1].y + f1.normal.z * f2.vertex[1].z + d) > SLICE_EPS &&
            fabs(f1.normal.x * f2.vertex[2].x + f1.normal.y * f2.vertex[2].y + f1.normal.z * f2.vertex[2].z + d) > SLICE_EPS)
        return false; // none is on the plane

    // at least one vertex is on plane, so check for the intersection
    // do the full line / polygon intersection test

    // f1 plane with f2 facet
    for (i = 0; i < 3; i++) {
        if (fabs(f1.normal.x * f2.vertex[i].x + f1.normal.y * f2.vertex[i].y + f1.normal.z * f2.vertex[i].z + d) > SLICE_EPS)
            continue; /* point is not on the plane */

        // point is on plane, check if it is inside the facet

        /* Determine whether or not the intersection point is bounded by pa,pb,pc */
        pa1.x = f1.vertex[i].x - f2.vertex[i].x;
        pa1.y = f1.vertex[i].y - f2.vertex[i].y;
        pa1.z = f1.vertex[i].z - f2.vertex[i].z;
        norm = pa1.norm();
        if (norm > 0.0f) {
            pa1.x /= norm;
            pa1.y /= norm;
            pa1.z /= norm;
        } else
            return true; // the point is coincident and we are neighbours facets
        pa2.x = f1.vertex[(i + 1) % 3].x - f2.vertex[i].x;
        pa2.y = f1.vertex[(i + 1) % 3].y - f2.vertex[i].y;
        pa2.z = f1.vertex[(i + 1) % 3].z - f2.vertex[i].z;
        norm = pa2.norm();
        if (norm > 0.0f) {
            pa2.x /= norm;
            pa2.y /= norm;
            pa2.z /= norm;
        } else
            return true;
        pa3.x = f1.vertex[(i + 2) % 3].x - f2.vertex[i].x;
        pa3.y = f1.vertex[(i + 2) % 3].y - f2.vertex[i].y;
        pa3.z = f1.vertex[(i + 2) % 3].z - f2.vertex[i].z;
        norm = pa3.norm();
        if (norm > 0.0f) {
            pa3.x /= norm;
            pa3.y /= norm;
            pa3.z /= norm;
        } else
            return true;

        a1 = pa1.x * pa2.x + pa1.y * pa2.y + pa1.z * pa2.z;
        a2 = pa2.x * pa3.x + pa2.y * pa3.y + pa2.z * pa3.z;
        a3 = pa3.x * pa1.x + pa3.y * pa1.y + pa3.z * pa1.z;
        total = RAD_TO_DEG((acos(a1) + acos(a2) + acos(a3)));
        if (fabs(total - 360) < SLICE_EPS)
            return true;
    }

    // second facet plane is d=-(nx*x+ny*y+nz*z)
    d = -f2.normal.x * f2.vertex[0].x - f2.normal.y * f2.vertex[0].y - f2.normal.z * f2.vertex[0].z;

    // f2 plane with f1 facet
    for (i = 0; i < 3; i++) {
        if (fabs(f2.normal.x * f1.vertex[i].x + f2.normal.y * f1.vertex[i].y + f2.normal.z * f1.vertex[i].z + d) > SLICE_EPS)
            continue; /* point is not on the plane */

        // point is on plane, check if it is inside the facet

        /* Determine whether or not the intersection point is bounded by pa, pb,pc */
        pa1.x = f2.vertex[i].x - f1.vertex[i].x;
        pa1.y = f2.vertex[i].y - f1.vertex[i].y;
        pa1.z = f2.vertex[i].z - f1.vertex[i].z;
        norm = pa1.norm();
        if (norm > 0.0f) {
            pa1.x /= norm;
            pa1.y /= norm;
            pa1.z /= norm;
        } else
            return true; // the point is coincident and we are neighbours facets
        pa2.x = f2.vertex[(i + 1) % 3].x - f1.vertex[i].x;
        pa2.y = f2.vertex[(i + 1) % 3].y - f1.vertex[i].y;
        pa2.z = f2.vertex[(i + 1) % 3].z - f1.vertex[i].z;
        norm = pa2.norm();
        if (norm > 0.0f) {
            pa2.x /= norm;
            pa2.y /= norm;
            pa2.z /= norm;
        } else
            return true;
        pa3.x = f2.vertex[(i + 2) % 3].x - f1.vertex[i].x;
        pa3.y = f2.vertex[(i + 2) % 3].y - f1.vertex[i].y;
        pa3.z = f2.vertex[(i + 2) % 3].z - f1.vertex[i].z;
        norm = pa3.norm();
        if (norm > 0.0f) {
            pa3.x /= norm;
            pa3.y /= norm;
            pa3.z /= norm;
        } else
            return true;


        a1 = pa1.x * pa2.x + pa1.y * pa2.y + pa1.z * pa2.z;
        a2 = pa2.x * pa3.x + pa2.y * pa3.y + pa2.z * pa3.z;
        a3 = pa3.x * pa1.x + pa3.y * pa1.y + pa3.z * pa1.z;
        total = RAD_TO_DEG((acos(a1) + acos(a2) + acos(a3)));
        if (fabs(total - 360) < SLICE_EPS)
            return true;
    }

    return false;
}

// given a STL object return a set of integers with the facets that belong to a object part

int split_parts(stl_file *stl, std::set<std::set<int>> &parts_idx, std::set<pair<int, int>> &connections, std::set<int> &parts_floor) {
    std::set<int> seen_facets;
    std::set<int>::iterator w, itp1, itp2;
    std::set<std::set<int>>::iterator it1, it2;
    bool connected;

    int n, neigh, i, j;

    // part number
    n = 0;
    for (i = 0; i < stl->stats.number_of_facets; ++i) {
        if (seen_facets.find(i) != seen_facets.end())continue;
        // we have an unseen facet, so start a part
        std::set<int> part;
        std::set<int> working;
        part.insert(i);
        seen_facets.insert(i);
        working.insert(i);
        while (!working.empty()) {
            // insert neighbors
            w = working.begin();
            for (j = 0; j < 3; ++j) {
                neigh = stl->neighbors_start[*w].neighbor[j];
                if (neigh != -1 && seen_facets.find(neigh) == seen_facets.end()) {
                    part.insert(neigh);
                    working.insert(neigh);
                    seen_facets.insert(neigh);
                }
            }
            working.erase(w);
        }
        parts_idx.insert(part);
        n++;
    }

    // we need to build a connection graph for the parts
    // for all parts
    for (i = 0, it1 = parts_idx.begin(); it1 != parts_idx.end(); i++, it1++) {
        // for all parts
        for (j = i+1, it2 = it1, it2++; it2 != parts_idx.end(); j++, it2++) {
            //if (it1 == it2)
            //    continue;
            //printf("%%Part\n");
            std::set<int> part1 = *it1;
            std::set<int> part2 = *it2;
            // for all facets check if there are two "equal"
            connected = false;
            for (itp1 = part1.begin(); itp1 != part1.end() && !connected; itp1++) {
                //printf("close all;figure;hold on;\n");
                //print_facet(stl->facet_start[*itp1]);

                for (itp2 = part2.begin(); itp2 != part2.end(); itp2++) {
                    //print_facet(stl->facet_start[*itp2]);
                    // compare facet indexed by itp1 to facet indexed by itp2
                    if (joint_facets(stl->facet_start[*itp1], stl->facet_start[*itp2])) {
                        connections.insert(std::make_pair(i, j));
                        printf("%d connected with %d\n", i, j);
                        connected = true;
                        break; // no need to proceed, since we already now the connection
                    }
                }
            }

        }
    }

    // check what parts are on floor
    for (i = 0, it1 = parts_idx.begin(); it1 != parts_idx.end(); i++, it1++) {
        // for all parts
        std::set<int> part1 = *it1;
        for (itp1 = part1.begin(); itp1 != part1.end(); itp1++) {
            //check if this facet is on the floor
            if (floor_facet(stl->facet_start[*itp1], stl->stats.min.z)) {
                parts_floor.insert(i);
                break;
            }

        }
    }


    return n;
}

void split_parts_layers(stl_file *stl, float x, float y) {
    float sliceSize = printer_char.extrusion_height, a;
    stl_file stl_local;
    int i, j, k, ls, ps, pos;
    vector<Polygon> poly;
    vector<v3> points;
    v3 d, d_previous;
    bool notfirst;
    std::vector<std::vector < LineSegment>> slicesWithLineSegments;
    std::vector<LineSegment> lines_segment;
    LineSegment line_segment;
    std::vector<Layer> Layers;
    FILE *f;


    stl_initialize(&stl_local);

    // make a copy, so rotate can be undone
    stl_local.stats = stl->stats;
    // allocate memory for facets
    stl_local.facet_start = (stl_facet*) malloc(stl_local.stats.number_of_facets * sizeof (stl_facet));
    if (stl_local.facet_start == NULL)
        perror("splite_parts_layers");
    // copy facets
    memcpy(stl_local.facet_start, stl->facet_start, stl->stats.number_of_facets * sizeof (stl_facet));

    // rotate along x and y axis
    //stl_rotate_xy(&stl_local, x, y);
    stl_rotate_x(&stl_local, x);
    stl_rotate_y(&stl_local, y);

    // Generate Slices
    triMeshSlicer(&stl_local, slicesWithLineSegments, sliceSize);

    // Build layers with polygons
    buildLayers(slicesWithLineSegments, Layers);

    fopen_s(&f, "Porra.m", "w");
    if (!f) return;

    ls = (int) Layers.size();
    ls = 1;
    for (i = 0; i < ls; ++i) {
        poly = Layers[i].getPoly();
        fprintf(f, "Layer%d=[", i);
        ps = (int) poly.size();
        for (j = 0; j < ps; ++j) {
            points = poly[j].getPoints();
            pos = (int) points.size();
            d_previous = v3(0, 0, 0);
            notfirst = false;
            for (k = 0; k < pos - 1; ++k) {
                fprintf(f, "%f,%f,%f;", points[k].x, points[k].y, points[k].z);
                d = points[k + 1] - points[k];
                a = (d.norm(), d_previous.norm()) * acos(abs(mymax(mymin(d.dotproduct(d_previous) / (d.norm(), d_previous.norm()), 1.0f), -1.0f)));
                if (a > DEG_TO_RAD(85)) {
                    if (notfirst) {
                        line_segment.v[1] = points[k];
                        notfirst = false;
                        lines_segment.push_back(line_segment);
                    } else {
                        line_segment.v[0] = points[k];
                        notfirst = true;
                    }
                }
                d_previous = d;
            }
            fprintf(f, "%f,%f,%f;", points[k].x, points[k].y, points[k].z);
            fprintf(f, "%f,%f,%f;", points[0].x, points[0].y, points[0].z);
            fprintf(f, "NaN,NaN,NaN;");
        }
        fprintf(f, "];\n\n");
        fprintf(f, "figure;mapshow(Layer%d(:,1),Layer%d(:,2),'DisplayType','polygon');hold on\n\n", i, i);
    }

    for (i = 0; i < (int) lines_segment.size(); i++) {
        fprintf(f, "line([%f,%f],[%f,%f],[%f,%f],'Color','r');\n", lines_segment[i].v[0].x, lines_segment[i].v[1].x, lines_segment[i].v[0].y, lines_segment[i].v[1].y, lines_segment[i].v[0].z, lines_segment[i].v[1].z);
    }

    free(stl_local.facet_start);
    fclose(f);
}

void usage(char *cmd) {
    printf("Usage:\n\n%s [-opt] stlfile.[stl,obj] \n\n", cmd);
    printf("  Where opt is:\n");
    printf("   -print          -> to produce the gcode for buinding\n");
    printf("   -inspection     -> to produce the gcode for inspection\n");
    printf("   -ex_height x    -> extruder height (current %f)\n", printer_char.extrusion_height);
    printf("   -ex_width x     -> extruder width (current %f)\n", printer_char.extrusion_width);
    printf("   -ex_mult x      -> extruder multiplier (current %f)\n", printer_char.extrusion_multiplier);
    printf("   -nozzle_dia x   -> nozzle diameter (current %f)\n", printer_char.nozzle_dia);
    printf("   -ins_slice x    -> extruder diameter (current %f)\n", inspection_char.slicing);
    printf("   -no_opt         -> disable optimal object orientation\n");
    printf("   -no_brim        -> disable brim export\n");
    printf("   -no_fullbrim    -> disable full brim export\n");
    printf("   -no_simulator   -> do not generate code for simulator (default)");
    printf("   -simulator      -> generate code for simulator (avoid simulator bugs and brim processing)");
    printf("   -ins_lift_off x -> distance between inspection head and object (current %f)\n", inspection_char.distance_to_part);
    printf("   -ins_sampling_distance x -> distance between sampling points (current %f)\n", inspection_char.sampling_distance);
    printf("   -ins_speed x    -> speed velocity for inspection (current %f)\n", inspection_char.speed);
    printf("   -MIP            -> use the MIP to create GCODE file\n");
    printf("   -CombF          -> use the CombF to create GCODE file\n");
    printf("   -NNH            -> use the NNH to create GCODE file\n");
    printf("   -NNSH x         -> use the NNSH to create GCODE file and specify the x number of possible path\n");
    printf("   -NNGH x         -> use the NNGH to create GCODE file and specify the x number of possible path\n");
    printf("   -ins_printLines -> print line number in GCODE file (current %d)\n", inspection_char.printGcodeLines);
    printf("   -ins_print      -> print SCAD info into console (current %d)\n", inspection_char.action_print);
}

int main(int argc, char **argv) {
    char modelFileName[1024];
    int modelFileNameLength = 0;
    int current_arg;
    bool action_print = false;
    bool action_inspection = false;
    bool action_optimization = true;
    double lb[2], ub[2], x0[2];
    int fixall_flag = 1; /* Default behavior is to fix all. */
    int exact_flag = 1; /* All checks turned off by default. */
    int tolerance_flag = 0; /* Is tolerance specified on cmdline */
    int nearby_flag = 0;
    int remove_unconnected_flag = 1;
    int fill_holes_flag = 1;
    int normal_directions_flag = 0;
    int normal_values_flag = 0;
    int reverse_all_flag = 0;
    float tolerance = 0;
    float increment = 0;
    int iterations = 2; /* Default number of iterations. */
    int increment_flag = 0;
    double *A = NULL, *b = NULL, *sol, f;
    int exit_code, i, part_number, first;
    //std::vector<feature_edge> feature_edges;


    if (argc < 2) {
        usage(argv[0]);
        return 0;
    }

    /* current argument to be dealed with */
    current_arg = 1;
    /* no default file name */
    modelFileName[0] = 0;
    while (argv[current_arg]) {
        if (!strncmp(argv[current_arg], "-print", 7)) {
            action_print = true;
            current_arg++;
        } else if (!strncmp(argv[current_arg], "-inspection", 12)) {
            action_inspection = true;
            current_arg++;
        } else if (!strncmp(argv[current_arg], "-no_opt", 8)) {
            action_optimization = false;
            current_arg++;
        } else if (!strncmp(argv[current_arg], "-no_simulator", 8)) {
            printer_char.simulator = false;
            current_arg++;
        } else if (!strncmp(argv[current_arg], "-simulator", 8)) {
            printer_char.simulator = true;
            current_arg++;
        } else if (!strncmp(argv[current_arg], "-no_brim", 9)) {
            printer_char.export_brim = false;
            current_arg++;
        } else if (!strncmp(argv[current_arg], "-no_fullbrim", 13)) {
            printer_char.export_full_brim = false;
            current_arg++;
        } else if (!strncmp(argv[current_arg], "-ex_height", 11)) {
            if (!argv[current_arg + 1]) {
                printf("\nMissing extrusion height\n\n");
                usage(argv[0]);
                return 1;
            }
            printer_char.extrusion_height = (float) atof(argv[current_arg + 1]);
            current_arg += 2;
        } else if (!strncmp(argv[current_arg], "-ex_width", 10)) {
            if (!argv[current_arg + 1]) {
                printf("\nMissing extrusion width\n\n");
                usage(argv[0]);
                return 1;
            }
            printer_char.extrusion_width = (float) atof(argv[current_arg + 1]);
            current_arg += 2;
        } else if (!strncmp(argv[current_arg], "-ex_mult", 9)) {
            if (!argv[current_arg + 1]) {
                printf("\nMissing extrusion multiplier\n\n");
                usage(argv[0]);
                return 1;
            }
            printer_char.extrusion_multiplier = (float) atof(argv[current_arg + 1]);
            current_arg += 2;
        } else if (!strncmp(argv[current_arg], "-nozzle_dia", 12)) {
            if (!argv[current_arg + 1]) {
                printf("\nMissing nozzle diameter\n\n");
                usage(argv[0]);
                return 1;
            }
            printer_char.nozzle_dia = (float) atof(argv[current_arg + 1]);
            current_arg += 2;
        } else if (!strncmp(argv[current_arg], "-ins_slice", 11)) {
            if (!argv[current_arg + 1]) {
                printf("\nMissing inspection slicing\n\n");
                usage(argv[0]);
                return 1;
            }
            inspection_char.slicing = (float) atof(argv[current_arg + 1]);
            current_arg += 2;
        } else if (!strncmp(argv[current_arg], "-ins_lift_off", 14)) {
            if (!argv[current_arg + 1]) {
                printf("\nMissing lift-off value\n\n");
                usage(argv[0]);
                return 1;
            }
            inspection_char.distance_to_part = (float) atof(argv[current_arg + 1]);
            current_arg += 2;
        } else if (!strncmp(argv[current_arg], "-ins_sampling_distance", 23)) {
            if (!argv[current_arg + 1]) {
                printf("\nMissing sampling distance value\n\n");
                usage(argv[0]);
                return 1;
            }
            inspection_char.sampling_distance = (float) atof(argv[current_arg + 1]);
            current_arg += 2;
        } else if (!strncmp(argv[current_arg], "-ins_speed", 11)) {
            if (!argv[current_arg + 1]) {
                printf("\nMissing speed value\n\n");
                usage(argv[0]);
                return 1;
            }
            inspection_char.speed = (float) atof(argv[current_arg + 1]);
            current_arg += 2;
        } else if (!strncmp(argv[current_arg], "-MIP", 5)) {
            inspection_char.MIP = 1;
            current_arg++;
        } else if (!strncmp(argv[current_arg], "-NNH", 5)) {
            inspection_char.NNH = 1;
            current_arg++;
        } else if (!strncmp(argv[current_arg], "-CombF", 7)) {
            inspection_char.CombF = 1;
            current_arg++;
        } else if (!strncmp(argv[current_arg], "-NNSH", 8)) {
            if (!argv[current_arg + 1]) {
                printf("\nMissing speed value\n\n");
                usage(argv[0]);
                return 1;
            }
            inspection_char.k_NNSH = (float) atof(argv[current_arg + 1]);
            if (inspection_char.k_NNSH <= 0) {
                usage(argv[0]);
                return 1;
            }
            current_arg += 2;
        } else if (!strncmp(argv[current_arg], "-NNGH", 8)) {
            if (!argv[current_arg + 1]) {
                printf("\nMissing speed value\n\n");
                usage(argv[0]);
                return 1;
            }
            inspection_char.k_NNGH = (float) atof(argv[current_arg + 1]);
            if (inspection_char.k_NNGH <= 0) {
                usage(argv[0]);
                return 1;
            }
            current_arg += 2;
        } else if (!strncmp(argv[current_arg], "-ins_printLines", 16)) {
            inspection_char.printGcodeLines = 1;
            current_arg++;
        } else if (!strncmp(argv[current_arg], "-ins_print", 11)) {
            inspection_char.action_print = true;
            current_arg++;
        } else if (argv[current_arg][0] == '-') {
            printf("\nInvalid option %s\n", argv[current_arg]);
            usage(argv[0]);
            return 1;
        } else if (argv[current_arg + 1]) {
            printf("\nSTL or OBJ filename must be the last argument.\n\n");
            usage(argv[0]);
            return 1;
        } else {
            strncpy(modelFileName, argv[current_arg], 1023);
            current_arg++;
        }

    }

    if (action_inspection &&
            inspection_char.MIP == 0 &&
            inspection_char.CombF == 0 &&
            inspection_char.k_NNSH == 0 &&
            inspection_char.k_NNGH == 0 &&
            inspection_char.NNH == 0) {
        printf("\nInspection selected, but no algorithm was selected...\nNothing to be done in inspection...\n\n");
        usage(argv[0]);
        return 0;
    }

    if (!action_print && !action_inspection) {
        printf("\nPrint or Inspection must be selected...\nNothing to be done...\n\n");
        usage(argv[0]);
        return 0;
    }

    if (!modelFileName[0]) {
        printf("\nSTL filename is mandatory.\n\n");
        usage(argv[0]);
        return 1;
    }

    modelFileNameLength = strlen(modelFileName);


    if (printer_char.extrusion_height <= 0) {
        printf("\nExtrusion height must be greater than zero.\n\n");
        usage(argv[0]);
        return 1;
    }

    if (printer_char.extrusion_width <= 0) {
        printf("\nExtrusion width must be greater than zero.\n\n");
        usage(argv[0]);
        return 1;
    }

    if (printer_char.extrusion_height > printer_char.extrusion_width) {
        printf("\nExtrusion width must be greater than or equal to the extruder height.\n\n");
        usage(argv[0]);
        return 1;
    }

    if (print_level)
        printf("\nUsing file %s with a Slice Step of %f.\n\n", modelFileName, printer_char.extrusion_height);

    // get extention to lower case
    if (modelFileNameLength > 4) {
        char *pointer = modelFileName + modelFileNameLength - 4;
        while (*pointer) {
            if (*pointer > 'A' && *pointer < 'Z')
                *pointer -= 'A' - 'a';
            pointer++;
        }
    }

    if (modelFileNameLength > 4 && !strncmp(modelFileName + modelFileNameLength - 4, ".stl", 4)) {
        // open STL file
        stl_open(&complete_mesh.stl, modelFileName);
        stl_exit_on_error(&complete_mesh.stl);
        if (complete_mesh.stl.error) {
            printf("\nError reading STL file.\n\n");
            return 1;
        }
    } else {
        if (modelFileNameLength > 4 && !strncmp(modelFileName + modelFileNameLength - 4, ".obj", 4)) {
            // open and load OBJ file
            stl_initialize(&complete_mesh.stl);

            if (complete_mesh.stl.error) {
                printf("\nError initializing stl struct.\n\n");
                return 1;
            }


            complete_mesh.stl.fp = fopen(modelFileName, "rb");

            if (complete_mesh.stl.fp == NULL) {
                complete_mesh.stl.error = 1;
                printf("\nError reading OBJ file.\n\n");
                return 1;
            } else {
                char lineHeader[1024];
                vector< unsigned int > vertexIndices, normalIndices;
                vector<v3> temp_vertices;
                vector<v3> temp_normals;

                while (fgets(lineHeader, 1024, complete_mesh.stl.fp) != NULL) {
                    //parse lineHeader
                    if (!strncmp(lineHeader, "v ", 2)) {
                        v3 vertex;
                        sscanf(lineHeader + 2, "%f %f %f", &vertex.x, &vertex.y, &vertex.z);
                        temp_vertices.push_back(vertex);
                    } else if (!strncmp(lineHeader, "vn ", 3)) {
                        v3 normal;
                        sscanf(lineHeader + 3, "%f %f %f", &normal.x, &normal.y, &normal.z);
                        temp_normals.push_back(normal);
                    } else if (!strncmp(lineHeader, "f ", 2)) {
                        unsigned int vertexIndex[3], uvIndex[3], normalIndex[3];
                        int matches = sscanf(lineHeader + 2, "%d/%d/%d %d/%d/%d %d/%d/%d", &vertexIndex[0], &uvIndex[0], &normalIndex[0], &vertexIndex[1], &uvIndex[1], &normalIndex[1], &vertexIndex[2], &uvIndex[2], &normalIndex[2]);
                        if (matches != 9) {
                            printf("File can't be read by our simple parser : (Try exporting with other options)\n");
                            usage(argv[0]);
                            fclose(complete_mesh.stl.fp);
                            return 1;
                        }
                        vertexIndices.push_back(vertexIndex[0]);
                        vertexIndices.push_back(vertexIndex[1]);
                        vertexIndices.push_back(vertexIndex[2]);
                        normalIndices.push_back(normalIndex[0]);
                        normalIndices.push_back(normalIndex[1]);
                        normalIndices.push_back(normalIndex[2]);
                    }
                }
                fclose(complete_mesh.stl.fp);

                complete_mesh.stl.stats.number_of_facets = vertexIndices.size() / 3;
                complete_mesh.stl.stats.original_num_facets = complete_mesh.stl.stats.number_of_facets;

                stl_allocate(&complete_mesh.stl);

                // convert vertices to STL
                int it, maxIt, facet_index, first = 1;

                maxIt = vertexIndices.size();
                for (it = 0, facet_index = 0; it < maxIt; it += 3, facet_index++) {
                    complete_mesh.stl.facet_start[facet_index].vertex[0].x = temp_vertices[vertexIndices[it] - 1].x;
                    complete_mesh.stl.facet_start[facet_index].vertex[0].y = temp_vertices[vertexIndices[it] - 1].y;
                    complete_mesh.stl.facet_start[facet_index].vertex[0].z = temp_vertices[vertexIndices[it] - 1].z;
                    complete_mesh.stl.facet_start[facet_index].vertex[1].x = temp_vertices[vertexIndices[it + 1] - 1].x;
                    complete_mesh.stl.facet_start[facet_index].vertex[1].y = temp_vertices[vertexIndices[it + 1] - 1].y;
                    complete_mesh.stl.facet_start[facet_index].vertex[1].z = temp_vertices[vertexIndices[it + 1] - 1].z;
                    complete_mesh.stl.facet_start[facet_index].vertex[2].x = temp_vertices[vertexIndices[it + 2] - 1].x;
                    complete_mesh.stl.facet_start[facet_index].vertex[2].y = temp_vertices[vertexIndices[it + 2] - 1].y;
                    complete_mesh.stl.facet_start[facet_index].vertex[2].z = temp_vertices[vertexIndices[it + 2] - 1].z;
                    complete_mesh.stl.facet_start[facet_index].normal.x = temp_normals[normalIndices[it] - 1].x;
                    complete_mesh.stl.facet_start[facet_index].normal.y = temp_normals[normalIndices[it] - 1].y;
                    complete_mesh.stl.facet_start[facet_index].normal.z = temp_normals[normalIndices[it] - 1].z;
                    stl_facet_stats(&complete_mesh.stl, complete_mesh.stl.facet_start[facet_index], first);
                    first = 0;
                }

                // CLEAN TMP VALUES
                //stl_stats_out(&complete_mesh.stl, stdout, modelFileName);
                //stl_write_ascii(&complete_mesh.stl, "Pois", modelFileName);
            }
        } else {
            printf("\nModel file must be a .stl or .obj file.\n\n");
            usage(argv[0]);
            return 1;
        }
    }
    

    // remove filename extension
    if (!strncmp(modelFileName + strlen(modelFileName) - 4, ".stl", 4) || !strncmp(modelFileName + strlen(modelFileName) - 4, ".obj", 4))
        modelFileName[strlen(modelFileName) - 4] = 0;

    if (print_level)
        fprintf(stderr, "\nObject has %d facets.\n\n", complete_mesh.stl.stats.number_of_facets);

    complete_mesh.normalize();

    stl_repair(&complete_mesh.stl,
            fixall_flag,
            exact_flag,
            tolerance_flag,
            tolerance,
            increment_flag,
            increment,
            nearby_flag,
            iterations,
            remove_unconnected_flag,
            fill_holes_flag,
            normal_directions_flag,
            normal_values_flag,
            reverse_all_flag,
            1);

    if (print_level)
        printf("\nObject has %d facets after repair.\n\n", complete_mesh.stl.stats.number_of_facets);

    if (action_optimization) {
        if (print_level)
            printf("\nGoing for object optimal orientation.\n\n");

        /* lower and upper bounds */
        lb[0] = lb[1] = 0.0;
        ub[0] = ub[1] = 180.0;
        /* initial guess - always consider no rotation as an initial guess*/
        x0[0] = 0.0;
        x0[1] = 0.0;

        stl_optimize = &complete_mesh.stl;
        opt.IPrint = 0;
        exit_code = PSwarm(2, &objfun, lb, ub, 0, A, b, &sol, &f, x0);

        if (print_level)
            printf("Object orientation optimal value of (%f,%f) degrees\n", sol[0], sol[1]);

        complete_mesh.optimal_rotations.push_back(v3((float) sol[0], (float) sol[1], 0.0f));

        // sol is allocated by pswarm
        free(sol);
    } else {
        printf("\nNo object optimal orientation.\n\n");
        complete_mesh.optimal_rotations.push_back(v3(0.0f, 0.0f, 0.0f));
    }

    
    if (action_inspection) {
        // output slice for inspection
        // we are inspecting the full object
        if (print_level)
            printf("\nStarting inspection... \n");

        Slice_inspection(&complete_mesh, modelFileName);

        if (print_level)
            printf("\nEnd of inspection.\n\n");
    }

    //split_parts_layers(&complete_mesh.stl, complete_mesh.optimal_rotation.x, complete_mesh.optimal_rotation.y);

    if (action_print) {
        // produce the GCODE for a 5-axis printer

        //check if optimizal rotation is different from the original stl
        if (complete_mesh.optimal_rotations[0].x ||
                complete_mesh.optimal_rotations[0].y ||
                complete_mesh.optimal_rotations[0].z) {
            printf("Rotating the object according to the optimal calculated (x:%f  y:%f  z:%f)\n",
                    complete_mesh.optimal_rotations[0].x,
                    complete_mesh.optimal_rotations[0].y,
                    complete_mesh.optimal_rotations[0].z);
            if (complete_mesh.optimal_rotations[0].x)
                stl_rotate_x(&complete_mesh.stl, complete_mesh.optimal_rotations[0].x);
            if (complete_mesh.optimal_rotations[0].y)
                stl_rotate_y(&complete_mesh.stl, complete_mesh.optimal_rotations[0].y);
            if (complete_mesh.optimal_rotations[0].z)
                stl_rotate_z(&complete_mesh.stl, complete_mesh.optimal_rotations[0].z);

            // we are done with the rotation
            complete_mesh.optimal_rotations[0].x = 0.0f;
            complete_mesh.optimal_rotations[0].y = 0.0f;
            complete_mesh.optimal_rotations[0].z = 0.0f;
            
            //complete_mesh.normalize();

            printf("rotation done!\n");

        } else {

            printf("Original position is the optimal one!\n");

        }

        // displacement to bring object to printing table (z=0), same as doing a translation
        complete_mesh.displacement.x = 0;
        complete_mesh.displacement.y = 0;
        complete_mesh.displacement.z = -complete_mesh.stl.stats.min.z;
        


        if (print_level) {
            printf("\nFinding object parts... ");
            // printf("\nComplete mesh bounding box: \n");
            // printf("min: %f, %f, %f\n", complete_mesh.stl.stats.min.x, complete_mesh.stl.stats.min.y, complete_mesh.stl.stats.min.z);
            // printf("min: %f, %f, %f\n", complete_mesh.stl.stats.max.x, complete_mesh.stl.stats.max.y, complete_mesh.stl.stats.max.z);
            // printf("displacement: %f,%f,%f\n", complete_mesh.displacement.x, complete_mesh.displacement.y, complete_mesh.displacement.z);
        }

        // start by spliting the object in several parts, if available
        std::set<std::set<int>> parts_idx;
        std::set<pair<int, int>> connections;
        std::set<int> parts_floor;

        // find parts, facets idx's, and connections
        split_parts(&complete_mesh.stl, parts_idx, connections, parts_floor);

        vector<Mesh> Meshes_part;
        //set<int> P_bar;
        vector<int> P_bar;

        if (print_level)
            printf(" found %d parts with %d connections, and %d part(s) is(are) on floor\n\n", (int) parts_idx.size(), (int) connections.size(), (int) parts_floor.size());
        
        //set the path for the table parts
        char * h1 = (char*) "PrinterParts\\head1.stl";
        char * h2 = (char*) "PrinterParts\\head2.stl";
        char * h3 = (char*) "PrinterParts\\table0.stl";

        //open head1
        stl_file head1;
        stl_open(&head1, h1);
        stl_exit_on_error(&head1);

        //open head2
        stl_file head2;
        stl_open(&head2, h2);
        stl_exit_on_error(&head2);

        //open table
        stl_file table;
        stl_open(&table, h3);
        stl_exit_on_error(&table);
        

        if (parts_idx.size() > 1) {
            // more than one part is available
            std::set<std::set<int>>::iterator it;

            // go for all parts
            for (part_number = 0, it = parts_idx.begin(); it != parts_idx.end(); it++, part_number++) {
                struct Mesh m_part;
                std::set<int> part = *it;
                std::set<int>::iterator itp;

                stl_initialize(&m_part.stl);
                m_part.stl.stats.number_of_facets = (int) part.size();
                stl_allocate(&m_part.stl);

                first = 1;
                for (i = 0, itp = part.begin(); itp != part.end(); itp++) {
                    m_part.stl.facet_start[i] = complete_mesh.stl.facet_start[*itp];
                    stl_facet_stats(&m_part.stl, m_part.stl.facet_start[i++], first);
                    first = 0;
                }

                // center part at origin
                m_part.normalize();

                if (print_level)
                    printf("\nPart %d has %d facets\n", part_number, m_part.stl.stats.number_of_facets);

                stl_repair(&m_part.stl,
                        fixall_flag,
                        exact_flag,
                        tolerance_flag,
                        tolerance,
                        increment_flag,
                        increment,
                        nearby_flag,
                        iterations,
                        remove_unconnected_flag,
                        fill_holes_flag,
                        normal_directions_flag,
                        normal_values_flag,
                        reverse_all_flag,
                        1);

                if (print_level)
                    printf("and %d facets after repair.\n\n", m_part.stl.stats.number_of_facets);

                // <editor-fold defaultstate="collapsed" desc="orientation for each part...">
                if (action_optimization) {
                    if (print_level)
                        printf("\nGoing for part optimal orientation.\n\n");

                    // MMDS call to find all the optimal rotations
                    // simple bounds and initial guess
                    arma::mat lbb(2, 1), ubb(2, 1), xx(2, 4); // xx(nParts * 2, 1);
                    // bounds
                    lbb.fill(0); ubb.fill(180);
                    // intitial guesses
                    xx.fill(0); xx(0,1)=180; xx(1,1)=180; xx(0,2)=0; xx(1,2)=90; xx(0,3)=90; xx(1,3)=0;
                    // the problem structure
                    Problem p(lbb, ubb, xx);
                    // use default options
                    MMDSOption o;
                    // statistics
                    Stat s;
                    // the part to optimize
                    stl_optimize = &m_part.stl;

                    // MMDS object
                    MMDS mm;
                    mm.objf = &objfun; // objective function
                    arma::mat fX; // columnwise solutions
                    arma::colvec ff, Alphas; // objective function values and final alpha values
                    mm.executeMMDS(p, o, fX, ff, Alphas, s);

                    if (print_level) {
                        std::cout << "part: " << part_number;
                        std::cout << ", found " << fX.n_cols << " optimal rotations" << endl;
                        //std::cout << "Optimal rotations "; 
                        //fX.print("fX: ");
                        //std::cout << "Optimal objective funtion values"; 
                        //ff.print("f: ");                      
                        //std::cout << part_number << " & " << s.func_eval << " & " << arma::min(ff) << " & " << fX.n_cols << " & "; // LaTex
                    }

                    // process the optimal rotations and save them in an ascending order
                    m_part.optimal_rotations.clear();
                    stl_translate_relative(&(m_part.stl), m_part.displacement.x, m_part.displacement.y, m_part.displacement.z + complete_mesh.displacement.z);
                    //stl_translate_relative(&(m_part.stl), 0.0f, 0.0f, -complete_mesh.stl.stats.min.z+0.01f);
                    arma::uvec indices = sort_index(ff);
                    vector<v3> testedRotations;
                    testedRotations.clear();
                    for (int j = 0; j < fX.n_cols; j++){
                        float B, C, lastb=0.0f, lastc=0.0f, vx=0.0f, vy=0.0f, vz=1.0f;;
    
                        // compute B and C from rotation of the bounding box vector, i.e., (0,0,1)
                        my_rotate(&vz, &vx, -fX(1, indices(j)));
                        my_rotate(&vy, &vz, -fX(0, indices(j)));
    
                        // This rotation must be in accordance with the optimization procedure.
                        // compute angles for printer   
                        ANGLES(vx, vy, vz, &B, &C, &lastb, &lastc);
    
                        // angles for STL manipolation are in degrees
                        v3 Rot=v3(B, C, 0.0f);

                        // verificar se B e C já existiam

                        std::cout << "B=" << B << " C=" << C << endl;

                        bool already = false;
                        for (vector<v3>::iterator itRot = testedRotations.begin();
                                itRot != testedRotations.end(); itRot++) {
                            if (B < (*itRot).x + DEG_TO_RAD(1) && B > (*itRot).x - DEG_TO_RAD(1) &&
                                    C < (*itRot).y + DEG_TO_RAD(1) && C > (*itRot).y - DEG_TO_RAD(1)) {
                                already = true;
                                break;
                            }
                        }

                        if (!already) {
                            testedRotations.push_back(Rot);
                            if (!checkStructureCollision(&table, &head1, &head2, &(m_part.stl), &Rot)) {
                                m_part.optimal_rotations.push_back(Rot);
                            } else {
                                std::cout << "Rotation not added" << endl;
                            }
                        } else {
                            std::cout << "Already tested" << endl;
                        }
                    }
                    
                    if(print_level){
                        std::cout << "Saved " << m_part.optimal_rotations.size() << " optima" << endl;
                        //std::cout << m_part.optimal_rotations.size() << "\\\\" << endl; // LaTeX
                    }
                } else {
                    //stl_translate_relative(&(m_part.stl), m_part.displacement.x, m_part.displacement.y, m_part.displacement.z + complete_mesh.displacement.z);
                    stl_translate_relative(&(m_part.stl), 0.0f, 0.0f, m_part.stl.stats.min.z + complete_mesh.stl.stats.min.z);
                    m_part.optimal_rotations.push_back(v3(0.0f, 0.0f, 0.0f));
                }

                // save it
                Meshes_part.push_back(m_part);
                //P_bar.insert(part_number);
                P_bar.push_back(part_number);
            }
        } else {
            struct Mesh m_part;
            
            stl_initialize(&m_part.stl);
            m_part.stl.stats.number_of_facets = complete_mesh.stl.stats.number_of_facets;
            stl_allocate(&m_part.stl);

            memcpy(m_part.stl.facet_start, complete_mesh.stl.facet_start, sizeof(stl_facet)*complete_mesh.stl.stats.number_of_facets);
            memcpy(m_part.stl.neighbors_start, complete_mesh.stl.neighbors_start, sizeof(stl_neighbors)*complete_mesh.stl.stats.number_of_facets);
            stl_get_size(&(m_part.stl));
            
            // center part at origin
            //m_part.normalize(); // it´s a unique part, so no point in normalize

            if (print_level)
                printf("\nUnique part has %d facets\n", m_part.stl.stats.number_of_facets);

            stl_repair(&m_part.stl,
                    fixall_flag,
                    exact_flag,
                    tolerance_flag,
                    tolerance,
                    increment_flag,
                    increment,
                    nearby_flag,
                    iterations,
                    remove_unconnected_flag,
                    fill_holes_flag,
                    normal_directions_flag,
                    normal_values_flag,
                    reverse_all_flag,
                    1);

            if (print_level)
                printf("and %d facets after repair.\n\n", m_part.stl.stats.number_of_facets);

            // part is the same as the complete_stl, but translated
            stl_translate_relative(&(m_part.stl), m_part.displacement.x, m_part.displacement.y, m_part.displacement.z + complete_mesh.displacement.z);
            m_part.optimal_rotations.push_back(v3(0.0f, 0.0f, 0.0f));
            
            Meshes_part.push_back(m_part);
            P_bar.push_back(0);
        }

        
        //now that we have all parts and corresponding optimal rotations,
        // try fo find an optimal printing sequence
        vector<vector<pair<int, int>>> Lt;
        vector<pair<int, int>> Lc;

        L_f.clear();
        level(&table, &head1, &head2, &Meshes_part, &parts_floor, &connections, Lt, Lc, P_bar);

        if (print_level) {
            printf("Found %d valid sequence(s)\n\n", (int) L_f.size());
            //print the result
            if(L_f.size())
                print_List_List(&Meshes_part);
        }

        if (print_level)
            printf("Going to generate GCode files for each sequence.\n");

        if(L_f.size())
            slice_filling_master(&L_f, &Meshes_part, modelFileName);

        if (print_level)
            printf("End creating filling GCode.\n");

        for (vector<Mesh>::iterator it = Meshes_part.begin(); it != Meshes_part.end(); it++) {
            stl_close(&((*it).stl));
        }
        stl_close(&complete_mesh.stl);
    }

    printf("press any key to exit...");
    getchar();
    return 0;
}