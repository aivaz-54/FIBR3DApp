/** 
 * file  FillingAlgorithms.cpp
 * author Bruna
 * date Created on 19 de Julho de 2019, 11:31
 * brief Methods to create filling patterns and generate corresponding GCode
 */


#include <limits>
#include "FillingAlgorithms.hpp"
#include "FIBR3Dapp.hpp"

extern printer printer_char;
extern void ANGLES(float vx, float vy, float vz, float *B, float *C, float *lastb, float *lastc);
extern void my_rotate(float *x, float *y, float angle);
extern void my_stl_rotate(stl_file *stl, float B, float C);
extern int line_number;
extern void exportPreamble(FILE *f);
extern int exportSLayersGCode5Axes_shell(std::vector<SLayer> &SLayers, float part_size, bool GCode5D, FILE *f, float z_offset);


float total_extruded_filling = 0.0f; // total extruded from begining

double fill_distance(v3 from, v3 to) {
    return sqrt(pow(from.x - to.x, 2) + pow(from.y - to.y, 2) + pow(from.z - to.z, 2));
}

double lengthSegment(LineSegment *segment) {
    v3 from = segment->v[0];
    v3 to = segment->v[1];
    return sqrt(pow(from.x - to.x, 2) + pow(from.y - to.y, 2) + pow(from.z - to.z, 2));
}

LineSegment *longestSegment(vector<LineSegment*> segments) {
    double longest = -1;
    double l;
    LineSegment *Smax;
    for (LineSegment *s : segments) {
        l = lengthSegment(s);
        if (longest < l) {
            Smax = s;
            longest = l;
        }
    }
    return Smax;
}

int segmentLineIntersect(LineSegment *s1, LineSegment *s2, v3 *point) {

    //line AB (s1) represented as a1x + b1y = c1
    double a1 = s1->v[1].y - s1->v[0].y;
    double b1 = s1->v[0].x - s1->v[1].x;
    double c1 = a1 * (s1->v[0].x) + b1 * (s1->v[0].y);

    //line CD (s2) represented as a2x + b2y = c2
    double a2 = s2->v[1].y - s2->v[0].y;
    double b2 = s2->v[0].x - s2->v[1].x;
    double c2 = a2 * (s2->v[0].x) + b2 * (s2->v[0].y);

    double determinant = a1 * b2 - a2*b1;

    if (determinant == 0) {
        if (a2 == 0 && s2->v[0].y != s1->v[0].y) {
            //parallel and horizontal
            return 0;
        } else if (a2 == 0 && s2->v[0].y == s1->v[0].y) {
            //collinear and horizontal
            return 2;
        } else if (b2 == 0 && s2->v[0].x != s1->v[0].x) {
            //parallel and vertical
            return 0;
        } else if (b2 == 0 && s2->v[0].x == s1->v[0].x) {
            //colinear and vertical
            return 2;
        } else {
            double tx = (s1->v[0].x - s2->v[0].x) / (-b1);
            double ty = (s1->v[0].y - s2->v[0].y) / (a1);

            if (tx == ty) {
                //colinear
                return 2;
            } else {
                //parallel
                return 0;
            }

        }
    } else {
        double x = (b2 * c1 - b1 * c2) / determinant;
        double y = (a1 * c2 - a2 * c1) / determinant;

        double minX = min(s1->v[0].x, s1->v[1].x);
        double maxX = max(s1->v[0].x, s1->v[1].x);
        double minY = min(s1->v[0].y, s1->v[1].y);
        double maxY = max(s1->v[0].y, s1->v[1].y);

        x = roundf(x * 1000000) / 1000000;
        y = roundf(y * 1000000) / 1000000;

        minX = roundf(minX * 1000000) / 1000000;
        maxX = roundf(maxX * 1000000) / 1000000;
        minY = roundf(minY * 1000000) / 1000000;
        maxY = roundf(maxY * 1000000) / 1000000;

        if (x >= minX &&
                x <= maxX &&
                y >= minY &&
                y <= maxY) {
            point->x = x;
            point->y = y;
            point->z = s1->v[0].z;
            return 1;
        } else {
            return 0;
        }
    }
}

LineSegment *findParallelSegment(LineSegment *segment, double d, bool direction) {
    double vABx = segment->v[1].y - segment->v[0].y;
    double vABy = segment->v[1].x - segment->v[0].x;

    if (direction) {
        vABy *= -1;
    } else {
        vABx *= -1;
    }

    double magnitude = sqrt(pow(vABx, 2) + pow(vABy, 2));
    vABx /= magnitude;
    vABy /= magnitude;

    double pAx = segment->v[0].x + d*vABx;
    double pAy = segment->v[0].y + d*vABy;
    double pBx = segment->v[1].x + d*vABx;
    double pBy = segment->v[1].y + d*vABy;

    LineSegment *s = new LineSegment(v3(pAx, pAy, segment->v[0].z), v3(pBx, pBy, segment->v[1].z), segment->normal); //v3(vABx, vABy, 0));
    return s;
}

LineSegment *findPerpendicularSegment(LineSegment *segment) {
    double vABy = -(segment->v[1].y - segment->v[0].y);
    double vABx = segment->v[1].x - segment->v[0].x;

    double middleX = (segment->v[1].y + segment->v[0].y) / 2;
    double middleY = (segment->v[1].x + segment->v[0].x) / 2;

    double d = 2;
    double pAx = middleY + d*vABy;
    double pAy = middleX + d*vABx;

    d = -2;
    double pBx = middleY + d*vABy;
    double pBy = middleX + d*vABx;


    LineSegment *s = new LineSegment(v3(pAx, pAy, segment->v[0].z), v3(pBx, pBy, segment->v[1].z), segment->normal); //v3(vABx, vABy, 0));
    return s;
}

vector<v3*> fill_sort(vector<v3*> points) {
    vector<v3*> sortedPoints;

    if (points.size() == 1) {
        sortedPoints.push_back(points.front());
    } else {
        v3 * start = points.at(0);
        v3 * maxP, * minP;
        int index, i;
        double max = 0;
        double min;
        double dist;

        for (i = 1; i < points.size(); i++) {
            dist = sqrt(pow(start->x - points.at(i)->x, 2) + pow(start->y - points.at(i)->y, 2));
            if (dist > max) {
                max = dist;
                maxP = points.at(i);
                index = i;
            }
        }

        sortedPoints.push_back(maxP);
        points.erase(points.begin() + index);

        while (!points.empty()) {
            min = std::numeric_limits<double>::max();
            for (i = 0; i < points.size(); i++) {
                dist = sqrt(pow(maxP->x - points.at(i)->x, 2) + pow(maxP->y - points.at(i)->y, 2));
                if (dist < min) {
                    min = dist;
                    minP = points.at(i);
                    index = i;
                }
            }
            sortedPoints.push_back(minP);
            points.erase(points.begin() + index);
        }
    }
    return sortedPoints;
}

bool onSegment(v3 *p, v3 *q, v3 *r) {

    double pr = sqrt(pow(p->x - r->x, 2) + pow(p->y - r->y, 2) + pow(p->z - r->z, 2));
    double pq = sqrt(pow(p->x - q->x, 2) + pow(p->y - q->y, 2) + pow(p->z - q->z, 2));
    double rq = sqrt(pow(r->x - q->x, 2) + pow(r->y - q->y, 2) + pow(r->z - q->z, 2));

    if (pr == pq + rq)
        return true;

    return false;
}

int pointInsidePolygon(vector<LineSegment*> polygon, v3 *point) {
    int c = 0;
    for (LineSegment *s : polygon) {
        if (((s->v[0].y > point->y) != (s->v[1].y > point->y)) &&
                (point->x < (s->v[1].x - s->v[0].x) * (point->y - s->v[0].y) / (s->v[1].y - s->v[0].y) + s->v[0].x)) {
            c = !c;
        }
    }
    return c;
}

v3 getMiddlePoint(LineSegment* segment) {
    double x = (segment->v[0].x + segment->v[1].x) / 2;
    double y = (segment->v[0].y + segment->v[1].y) / 2;
    double z = (segment->v[0].z + segment->v[1].z) / 2;
    return v3(x, y, z);
}

vector<LineSegment *> createSegments(vector<LineSegment *> polygon, vector<vector<LineSegment*>> holes, vector<v3*> points) {
    vector<LineSegment*> segments;
    v3 rotation;
    int i=0;

    // the same rotation for each segment in the polygon
    // so, get the first
    rotation=polygon.at(0)->normal;
    
    if (points.size() == 1) {
        LineSegment * s = new LineSegment(
                v3(points.at(i)->x, points.at(i)->y, points.at(i)->z),
                v3(points.at(i)->x, points.at(i)->y, points.at(i)->z),
                rotation);
        v3 middle = getMiddlePoint(s);
        if (pointInsidePolygon(polygon, &middle)) {
            bool insideHole = false;
            for (int h = 0; h < holes.size() && insideHole == false; h++) {
                insideHole = pointInsidePolygon(holes[h], &middle);
            }
            if (!insideHole) {
                segments.push_back(s);
            }
        }
    } else if (points.size() != 0) {


        for (i = 0; i < points.size() - 1; i++) {
            if (!(points.at(i)->x == points.at(i + 1)->x && points.at(i)->y == points.at(i + 1)->y && points.at(i)->z == points.at(i + 1)->z)) {
                LineSegment * s = new LineSegment(*points.at(i), *points.at(i + 1), rotation);
                v3 middle = getMiddlePoint(s);
                if (pointInsidePolygon(polygon, &middle)) {
                    bool insideHole = false;
                    for (int h = 0; h < holes.size() && insideHole == false; h++) {
                        insideHole = pointInsidePolygon(holes[h], &middle);
                    }
                    if (!insideHole) {
                        segments.push_back(s);
                    }
                }
            }
        }
    }
    return segments;

}

bool intersectBoundingLayer(vector<LineSegment *> box, LineSegment *parallel) {
    int j;
    v3 point;
    for (int i = 0; i < box.size(); i++) {
        j = segmentLineIntersect(box[i], parallel, &point);
        if (j == 1) {
            return true;
        }
    }
    return false;
}

void removeRepeatedPoints(vector<v3*> *points) {

    int s = points->size();
    int i;
    v3 *a, *b;
    for (i = s - 1; i > 0; i--) {
        a = points->at(i);
        b = points->at(i - 1);
        if ((a->x == b->x) && (a->y == b->y) && (a->z == b->z)) {
            points->erase(points->begin() + i);
        }
    }
}

LineSegment *adjustH(LineSegment *segment, float h) {
    LineSegment *s = new LineSegment();
    int i;
    for (i = 0; i < 2; i++) {
        s->v[i] = segment->v[i];
        s->v[i].z = h;
    }
    s->normal=segment->normal;
    return s;
}

void fillingAlgorithm(vector<vector<vector<LineSegment *>>> SL, vector<vector<double>> area_layers, double d, FILE *f) {
    // <editor-fold defaultstate="collapsed" desc="get the longest segment to create parallels">
    LineSegment * SmaxNormal = NULL;
    LineSegment * SmaxPerp = NULL;
    vector<vector<LineSegment *>> S;
    vector<double> area = area_layers.front();
    //vector<v3> rot = Rotations.front();
    double length = 0.0;

    if(SL.size()==0)
        return;
    S = SL.front();
    for (int i = 0; i < S.size(); i++) {
        //get the longest segment considering all polygons of the current layer
        //this excludes holes...
        if (area[i] > 0) {
            LineSegment *temp = longestSegment(S.at(i));
            double tempLength = lengthSegment(temp);
            if (tempLength > length) {
                SmaxNormal = temp;
                length = tempLength;
            }
        }
    }// </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="get the perpendicular of the longest segment">

    SmaxPerp = findPerpendicularSegment(SmaxNormal);
    // </editor-fold>

    int x, nLayers = SL.size();
    // <editor-fold defaultstate="collapsed" desc="calculate the lines for each layer">
    for (x = 0; x < nLayers; x++) {
        S = SL.at(x);
        area = area_layers.at(x);
        //rot=Rotations.at(x);

        LineSegment *Smax;
        if (x % 2 == 0) {
            Smax = adjustH(SmaxNormal, S.front().front()->v[0].z);
        } else {
            Smax = adjustH(SmaxPerp, S.front().front()->v[0].z);
        }

        vector<LineSegment*> box = createBoundingLayer(S);

        //create lines for each positive area
        for (int i = 0; i < S.size(); i++) {
            if (area[i] > 0) {
                // <editor-fold defaultstate="collapsed" desc="create vector with holes">
                //vector with holes...
                vector<vector<LineSegment*>> holes;
                for (int h = 0; h < area.size(); h++) {
                    if (h != i && area[h] < 0) {
                        if (pointInsidePolygon(S[i], &S[h][0]->v[0])) {
                            holes.push_back(S[h]);
                        }
                    }
                }// </editor-fold>

                // <editor-fold defaultstate="collapsed" desc="determine left/right parallels">
                LineSegment *SmaxL, *SmaxR;
                vector<vector<LineSegment *>> listOfParallelSegments;
                int countP = 0;

                // <editor-fold defaultstate="collapsed" desc="get left parallels">
                //find left parallels
                bool keep = true;
                do {
                    SmaxL = findParallelSegment(Smax, d*countP, false);
                    countP++;
                    vector<v3*> intersectionPoints;
                    vector<v3*> sortedIntersectionPoints;
                    int count = 0;
                    // <editor-fold defaultstate="collapsed" desc="parallels for the exterior polygon">
                    for (LineSegment *s : S[i]) {//parallels without holes

                        v3 *P = new v3();
                        int res = segmentLineIntersect(s, SmaxL, P);
                        if (res == 1) {
                            intersectionPoints.push_back(P);
                            count++;
                        } else if (res == 2) {
                            intersectionPoints.push_back(&(s->v[0]));
                            intersectionPoints.push_back(&(s->v[1]));
                            count += 2;
                        } else {
                            delete P;
                        }
                    }// </editor-fold>
                    // <editor-fold defaultstate="collapsed" desc="parallels for each hole">
                    for (int h = 0; h < holes.size(); h++) {//parallels for holes
                        for (LineSegment *hs : holes[h]) {

                            v3 *P = new v3();
                            int res = segmentLineIntersect(hs, SmaxL, P);
                            if (res == 1) {
                                intersectionPoints.push_back(P);
                                count++;
                            } else if (res == 2) {
                                intersectionPoints.push_back(&(hs->v[0]));
                                intersectionPoints.push_back(&(hs->v[1]));
                                count += 2;
                            } else {
                                delete P;
                            }
                        }
                    }// </editor-fold>


                    removeRepeatedPoints(&intersectionPoints);
                    if (count == 0) {
                        //check if the parallel is inside bounding box
                        if (!intersectBoundingLayer(box, SmaxL)) {
                            keep = false;
                        }
                    } else {
                        //only add vector if there are intersection points
                        //parallel may stand between polygon, so no intersection occurs
                        if (intersectionPoints.size() != 0) {
                            vector<LineSegment*> SegL;

                            sortedIntersectionPoints = fill_sort(intersectionPoints);

                            SegL = createSegments(S[i], holes, sortedIntersectionPoints);
                            if (SegL.size() != 0) {
                                listOfParallelSegments.push_back(SegL);
                            }
                        }
                    }
                } while (keep); // </editor-fold>

                //reverse the countP-1 elements
                reverse(listOfParallelSegments.begin(), listOfParallelSegments.end());
                // <editor-fold defaultstate="collapsed" desc="get right parallels">
                //find right parallels
                countP = 1;
                keep = true;
                do {
                    SmaxR = findParallelSegment(Smax, d*countP, true);
                    countP++;
                    vector<v3*> intersectionPoints;
                    vector<v3*> sortedIntersectionPoints;

                    int count = 0;
                    // <editor-fold defaultstate="collapsed" desc="intersection with the exterior polygon">
                    for (LineSegment *s : S[i]) {
                        v3 *P = new v3();
                        int res = segmentLineIntersect(s, SmaxR, P);
                        if (res == 1) {
                            intersectionPoints.push_back(P);
                            count++;
                        } else if (res == 2) {
                            intersectionPoints.push_back(&(s->v[0]));
                            intersectionPoints.push_back(&(s->v[1]));
                            count += 2;
                        } else {
                            delete P;
                        }
                    }// </editor-fold>

                    // <editor-fold defaultstate="collapsed" desc="intersection with holes">
                    for (int h = 0; h < holes.size(); h++) {//parallels for holes
                        for (LineSegment *hs : holes[h]) {
                            v3 *P = new v3();
                            int res = segmentLineIntersect(hs, SmaxR, P);
                            if (res == 1) {
                                intersectionPoints.push_back(P);
                                count++;
                            } else if (res == 2) {
                                intersectionPoints.push_back(&(hs->v[0]));
                                intersectionPoints.push_back(&(hs->v[1]));
                                count += 2;
                            } else {
                                delete P;
                            }
                        }
                    }// </editor-fold>
                    removeRepeatedPoints(&intersectionPoints);
                    if (count == 0) {
                        //check if the parallel is inside bounding box
                        if (!intersectBoundingLayer(box, SmaxR)) {
                            keep = false;
                        }
                    } else {
                        if (intersectionPoints.size() != 0) {
                            sortedIntersectionPoints = fill_sort(intersectionPoints);
                            vector<LineSegment*> SegR;
                            SegR = createSegments(S[i], holes, sortedIntersectionPoints);
                            if (SegR.size() != 0) {
                                listOfParallelSegments.push_back(SegR);
                            }
                        }
                    }

                } while (keep); // </editor-fold>

                // </editor-fold>

                vector<LineSegment*> ordered = orderSegmentsClosest(&listOfParallelSegments);

                fill_exportGCode5Axes_layer(ordered, f, &line_number);
            }
        }
        
        //free memory
        free(Smax);
    }
    // </editor-fold>

}


void fill_exportG0command(FILE *f, float speed, v3 *point, int *line, float B, float C) {
    // we should correct the printer bug (small to big angles rotates wrong)
    if ((*line) == -1) {
        fprintf(f, "G00 F%.2f X%.2f Y%.2f Z%.2f B%.2f C%.2f\n", speed, point->x, point->y, point->z, RAD_TO_DEG(B), RAD_TO_DEG(C)+180.0f);
    } else {
        fprintf(f, "N%d G00 F%.2f X%.2f Y%.2f Z%.2f B%.2f C%.2f\n", (*line)++, speed, point->x, point->y, point->z, RAD_TO_DEG(B), RAD_TO_DEG(C)+180.0f);
    }
}

void fill_exportG1command(FILE *f, v3 *point, float rate, int *line, float B, float C, float extrusion) {
    // we should correct the printer bug (small to big angles rotates wrong)
    if ((*line) == -1) {
        fprintf(f, "G01 F%.2f X%.2f Y%.2f Z%.2f B%.2f C%.2f A%.2f \n", printer_char.speed, point->x, point->y, point->z, RAD_TO_DEG(B), RAD_TO_DEG(C)+180.0f, extrusion);
    } else {
        fprintf(f, "N%d G01 F%.2f X%.2f Y%.2f Z%.2f B%.2f C%.2f A%.2f \n", (*line)++, printer_char.speed, point->x, point->y, point->z, RAD_TO_DEG(B), RAD_TO_DEG(C)+180.0f, extrusion);
    }
}

void fill_exportRcommand(FILE *f, float speed, float z) {
    
    fprintf(f, "N%d G00 F%.2f Z%.2f ;retract\n", line_number++, speed, z);
//    if(printer_char.simulator) // bug in the simulator
//        fprintf(f, "; G0 extra command to correct bug\nN%d G01 F%.2f Z%.2f ; retract\n", line_number++, speed, z);

    // fflush(f);
}

int fill_exportGCode5Axes_layer(vector<LineSegment*> path, FILE *f, int *line) {
   
    
    if (path.size() != 0) {
        float distance = sqrt(pow(path[0]->v[0].x - path[0]->v[1].x, 2) + pow(path[0]->v[0].y - path[0]->v[1].y, 2) + pow(path[0]->v[0].z - path[0]->v[1].z, 2));

        total_extruded_filling += distance * printer_char.extrusion_multiplier * (printer_char.extrusion_height * (printer_char.extrusion_width - printer_char.extrusion_height) + PI * pow(printer_char.extrusion_height / 2.0f, 2.0f)) /
                (PI * pow(printer_char.filament_radius, 2.0f) / 4.0f);

        fill_exportG0command(f, printer_char.speed, &path[0]->v[0], line, path[0]->normal.x, path[0]->normal.y);
        fill_exportG1command(f, &path[0]->v[1], printer_char.max_flow_rate, line, path[0]->normal.x, path[0]->normal.y, total_extruded_filling);



        for (int i = 1; i < path.size(); i++) {
            distance = sqrt(pow(path[i]->v[0].x - path[i]->v[1].x, 2) + pow(path[i]->v[0].y - path[i]->v[1].y, 2) + pow(path[i]->v[0].z - path[i]->v[1].z, 2));

            total_extruded_filling += distance * printer_char.extrusion_multiplier * (printer_char.extrusion_height * (printer_char.extrusion_width - printer_char.extrusion_height) + PI * pow(printer_char.extrusion_height / 2.0f, 2.0f)) /
                    (PI * pow(printer_char.filament_radius, 2.0f) / 4.0f);
            fill_exportG0command(f, printer_char.speed, &path[i]->v[0], line, path[i]->normal.x, path[i]->normal.y);
            fill_exportG1command(f, &path[i]->v[1], printer_char.max_flow_rate, line, path[i]->normal.x, path[i]->normal.y, total_extruded_filling);
        }
    }

    return 0;
}

vector<LineSegment*> orderSegmentsClosest(vector<vector<LineSegment*>> *segments) {
    vector<LineSegment*> list;

    v3 initial(0, 0, 0);

    for (int i = 0; i < segments->size(); i++) {
        vector<LineSegment*> temp = segments->at(i);

        if (temp.size() != 0) {

            //check who is closest. first or last?
            double distFront = fill_distance(initial, temp.front()->v[0]);
            double distBack = fill_distance(initial, temp.front()->v[1]);

            if (distFront > distBack) {
                //invert the list
                vector<LineSegment*> temp2;
                for (LineSegment * seg : temp) {
                    v3 pointTemp = seg->v[0];
                    seg->v[0] = seg->v[1];
                    seg->v[1] = pointTemp;
                    temp2.insert(temp2.begin(), seg);
                }
                for (LineSegment * s : temp2) {
                    list.push_back(s);
                }
            } else {
                for (LineSegment * s : temp) {
                    list.push_back(s);
                }
            }
        }
        initial = list.back()->v[1];
    }

    return list;
}

vector<LineSegment *> createBoundingLayer(vector<vector<LineSegment *>> S) {
    float minx = std::numeric_limits<float>::max();
    float miny = std::numeric_limits<float>::max();
    float maxx = std::numeric_limits<float>::min();
    float maxy = std::numeric_limits<float>::min();

    for (vector<LineSegment*> ls : S) {
        for (LineSegment* s : ls) {
            if (s->v[0].x < minx) {
                minx = s->v[0].x;
            }
            if (s->v[0].y < miny) {
                miny = s->v[0].y;
            }
            if (s->v[0].x > maxx) {
                maxx = s->v[0].x;
            }
            if (s->v[0].y > maxy) {
                maxy = s->v[0].y;
            }
        }
    }

    float z = S.front().front()->v[0].z;
    vector<LineSegment*> res;
    res.push_back(new LineSegment(v3(minx, miny, z), v3(minx, maxy, z), v3(0, 0, 0)));
    res.push_back(new LineSegment(v3(minx, maxy, z), v3(maxx, maxy, z), v3(0, 0, 0)));
    res.push_back(new LineSegment(v3(maxx, maxy, z), v3(maxx, miny, z), v3(0, 0, 0)));
    res.push_back(new LineSegment(v3(maxx, miny, z), v3(minx, miny, z), v3(0, 0, 0)));

    return res;
}

char* createDirectoryFile_filling(char* stl_file, char* directory, char* extension) {
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
    free(filename);
    return newFile;
}

// slice all L_f sequences, by slicing all L_t sequences
void slice_filling_master(vector<vector<vector<pair<int,int>>>> *L_f, vector<Mesh> *meshes, char *modelFileName) {
    //get the result
    int nBlock = 0, x;
    vector<vector<vector<pair<int,int>>>>::iterator k;
    FILE *file_desc;
    
    
    // go for all feasible sequencing that we found
    for (k = (*L_f).begin(), x = 0; k != (*L_f).end(); k++, x++) {

        char * integer_string = (char*) malloc(10 * sizeof (char));
        sprintf(integer_string, "%d", nBlock++);

        char* _modelFileName = (char*) malloc((strlen(modelFileName) + strlen(integer_string) + 1) * sizeof (char));
        strcpy(_modelFileName, modelFileName);
        strcat(_modelFileName, integer_string);

        char *_fname = createDirectoryFile_filling(_modelFileName, (char*) "Gcode_files/", (char*) "_filling.gcode");
        remove(_fname);

        free(integer_string);
        free(_modelFileName);
        
        file_desc = fopen(_fname, "w");
        free(_fname);
        if(!file_desc){
            std::cout << "Unable to create file " << _fname << "to export G-Code" << endl;
            return;
        }
        
        line_number=0;
        if(!printer_char.simulator)
            exportPreamble(file_desc); // preamble leads to problems in the simulator
        
        float z_offset=0.0f;
        // Brim and Support should be exported here
                
        vector<vector<pair<int,int>>> _Lt = (*k);
        int b = 0;
        vector<vector<pair<int,int>>>::iterator i;
        // for all printing levels
        for (b = 0, i = _Lt.begin(); i != _Lt.end(); i++, b++) {
            // slice each printing level (with possible many parts)
            slice_filling(meshes, (*i), file_desc, z_offset);
            // rectract between parts
            fill_exportRcommand(file_desc, printer_char.speed, printer_char.max_axis.z);
        }
    }
    fclose(file_desc);
}


void printPolygon_Scad(char *filename, vector<vector<vector<LineSegment *>>> polygons) {
    FILE *fp;

    //return;
    
    fp = fopen(filename, "w");
    if (fp){
        fprintf(fp,"module line(start, end, thickness = 1) {color(\"blue\") hull() {translate(start) sphere(thickness);translate(end) sphere(thickness);}}\n\n");
        for (vector<vector<vector<LineSegment*>>>::iterator p=polygons.begin(); p!=polygons.end(); p++) {
            for (vector<vector<LineSegment*>>::iterator i=(*p).begin(); i!=(*p).end(); i++) {
                for (vector<LineSegment*>::iterator j=(*i).begin(); j!=(*i).end(); j++) {
                    fprintf(fp, "line([%f,%f,%f],[%f,%f,%f],1);\n", (*j)->v[0].x, (*j)->v[0].y, (*j)->v[0].z,
                            (*j)->v[1].x, (*j)->v[1].y, (*j)->v[1].z);
                }
            }
        }
        fclose(fp);
    }
}


// slice parts in the same time slot
void slice_filling(vector<Mesh> *meshes, vector<pair<int,int>> L, FILE *f, float z_offset) {
    bool shell;
    float B, C;
    vector<vector<Layer>> LLayers;
    vector<v3> Rotations;
    v3 entryPoint = v3(0.0f,0.0f,0.0f);
    int maximumLayers=0;
    
    vector<pair<int,int>>::iterator j;
    for (j=L.begin(); j!=L.end();j++){
        Mesh *m=&(meshes->at((*j).first));
        stl_file *stl_local=&(m->stl);
        int Rot=(*j).second;

        B=m->optimal_rotations[Rot].x;          
        C=m->optimal_rotations[Rot].y;
    
        // z translate to printer table
        stl_translate_relative(stl_local, 0.0f, 0.0f, printer_char.plate_offset + z_offset); // printer_char.min_axis.z +
                
        // rotate along B and C printer angles
        my_stl_rotate(stl_local, B, C);
        
        // get to final position
        stl_translate_relative(stl_local, 0.0f, 0.0f, -printer_char.plate_offset);

        // Generate Slices
        vector<vector<LineSegment>> slicesWithLineSegments;
        triMeshSlicer(stl_local, slicesWithLineSegments, printer_char.extrusion_height);

        // Build layers with polygons
        vector<Layer> Layers;
        buildLayers(slicesWithLineSegments, Layers);
        correctStart(Layers, entryPoint);
        
        //exportLayersMATLABFormat3D(Layers,"Layers.m");
        
        LLayers.push_back(Layers);
        
        if(Layers.size()>maximumLayers)
            maximumLayers=Layers.size();
        Rotations.push_back(v3(B,C,0));
    }
    
    // Build spline layers with polygons
    std::vector<SLayer> SLayers;
    buildSplineLayers(LLayers.front(), SLayers);
    
    // export SLayers to MATLAB
    //exportSLayersMATLABFormat3D(SLayers, "SLayers.m");

    if(LLayers.size()==1 && (shell = is_shell(SLayers))){
        printf("IS SHELL!!!\n");
        stl_file *stl_local=&(meshes->at(L.begin()->first).stl);
                
        exportSLayersGCode5Axes_shell(SLayers,
                max(stl_local->stats.size.x, stl_local->stats.size.y), true, f, z_offset);

    } else {
        //create structures to store every layers
        vector<vector<vector<LineSegment *>>> polygons_layers;
        vector<vector<double>> area_layers;
        //vector<vector<v3>> LayersRot;
        
        // go for all layers
        for(int nLayer=0;nLayer<maximumLayers;nLayer++){

            vector<vector<Layer>>::iterator it;
            vector<v3>::iterator itRot;            
            vector<vector<LineSegment *>> polygons;
            vector<double> area;
            //vector<v3> rotations;

            // search all parts within this Layer
            for(itRot=Rotations.begin(), it=LLayers.begin(); itRot!=Rotations.end() && it!=LLayers.end(); itRot++, it++){
                vector<Layer> *Lay=&(*it);
                
                // do we have a part with a Layer at this altitude?!
                if(Lay->size()>nLayer){
                    vector<Polygon> vp = Lay->at(nLayer).getPoly();
                    int ps = vp.size();
                    for (int j = 0; j < ps; j++) {
                        vector<LineSegment *> segs;

                        Polygon p = vp.at(j);
                        int ks = p.getPoints().size();
                        for (int k = 0; k < ks - 1; k++) {
                            v3 from = p.getPoints().at(k);
                            v3 to = p.getPoints().at(k + 1);
                            LineSegment *lsegment = new LineSegment(from, to, (*itRot));
                            segs.push_back(lsegment);
                        }
                        polygons.push_back(segs);
                        area.push_back(p.getArea());
                        //rotations.push_back((*itRot));
                    }
                } // else, this part has exausted the Layers
            }

            polygons_layers.push_back(polygons);
            area_layers.push_back(area);
            //LayersRot.push_back(rotations);

        }
        
        //printPolygon_Scad("Polygons.scad", polygons_layers);

        fillingAlgorithm(polygons_layers, area_layers, printer_char.extrusion_width, f);

        //delete every lineSegment * 
        int a, b, c;
        for (a = polygons_layers.size() - 1; a >= 0; a--) {
            vector<vector<LineSegment*>> temp = polygons_layers.at(a);
            for (b = temp.size() - 1; b >= 0; b--) {
                vector<LineSegment*> temp2 = temp.at(b);
                for (c = temp2.size() - 1; c >= 0; c--) {
                    free(temp2.at(c));
                }
                temp2.clear();
            }
            temp.clear();
        }
        polygons_layers.clear();
        area_layers.clear();
    }
}
