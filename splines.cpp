#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <set>
#include "admesh/stl.h"
#include "FIBR3Dapp.hpp"

extern printer printer_char;


// compute the maximum distance from inner points to the line formed by outer points
// return the point index and the maximum distance

std::pair<int, float> findMaximumDistance(vector<v3>& Points) {
    int np = (int) Points.size(), i;
    v3 firstpoint = Points[0];
    v3 lastpoint = Points[np - 1];
    v3 pp;
    int index = 0; //index to be returned
    float Mdist = -1, p_norm, Dist; //the Maximum distance to be returned

    //distance calculation
    v3 p = lastpoint - firstpoint;
    // compute norm
    p_norm = p.norm();
    if (p_norm > 0.01f) { // we do not have a cycle...
        for (i = 1; i < np - 1; i++) { //traverse through second point to second last point
            pp = Points[i] - firstpoint;
            pp = p.crossproduct(pp);
            Dist = fabs(pp.norm()) / p_norm; //formula for point-to-line distance
            if (Dist > Mdist) {
                Mdist = Dist;
                index = i;
            }
        }
    }
    //else {
    //	printf("Too short\n");
    //}
    return std::make_pair(index, Mdist);
}


// Ramer-Douglas-Peucker-Algorithm
// remove points where distance to the line segments is lower than epsilon
// uses recursion

std::pair<vector<v3>, vector<v3>> simplifyWithRDP(vector<v3> Points, vector<v3> Normals, double epsilon) {

    if (Points.size() < 3) { //two points make a single line segment, so return with provided points
        return std::make_pair(Points, Normals);
    }

    // more than two points are provided
    // compute the maximum distance for all points
    std::pair<int, double> maxDistance = findMaximumDistance(Points);

    // where should we divide the set of points?
    int index = maxDistance.first;
    if (maxDistance.second < 0) {
        // by the middle, since we have a closed polygon
        index = (int) (Points.size() / 2);
    }
    // we have a point whose maximum distance is greater than epsilon, or we have a closed set of points
    if (maxDistance.second >= epsilon || maxDistance.second < 0) {
        // base  case 1, split at the index point and go recursively
        // the begin
        vector<v3>::iterator it = Points.begin();
        vector<v3>::iterator n_it = Normals.begin();
        // the two line segments
        vector<v3> path1(it, it + index + 1); //new path l1 from 0 to index
        vector<v3> path2(it + index, Points.end()); // new path l2 from index to last
        vector<v3> path1n(n_it, n_it + index + 1); //new path l1 from 0 to index
        vector<v3> path2n(n_it + index, Normals.end()); // new path l2 from index to last

        // recursively call
        std::pair<vector<v3>, vector < v3>> r1 = simplifyWithRDP(path1, path1n, epsilon);
        std::pair<vector<v3>, vector < v3>> r2 = simplifyWithRDP(path2, path2n, epsilon);

        //Concat simplified path1 and path2 together
        vector<v3> rs(r1.first);
        vector<v3> rn(r1.second);
        rn.pop_back();
        rs.pop_back();
        rs.insert(rs.end(), r2.first.begin(), r2.first.end());
        rn.insert(rn.end(), r2.second.begin(), r2.second.end());
        return std::make_pair(rs, rn);
    } else { //base case 2, all points between are to be removed.
        vector<v3> rs(1, Points[0]);
        rs.push_back(Points[Points.size() - 1]);
        vector<v3> rn(1, Normals[0]);
        rn.push_back(Normals[Normals.size() - 1]);
        return std::make_pair(rs, rn);
    }
}



/* solve a tridiagonal system */
// inf   - below diagonal
// diag  - diagonal
// sup   - above diagonal
// ind   - independent term
// sol   - tridiagonal linear system solution
// multi - multipliers list for warm start
// warm  - >0 if we have a warm start
// return 0 on success and >0 otherwise

int tridiagonal(int n, double *inf, double *diag, double *sup, double *ind,
        double *sol, double *multi, bool warm) {
    int i;
    double m;

    /* These are internal errors. It should not happen */
    if (!inf || !diag || !sup || !sol || n < 1)
        return 1;

    if (warm && !multi)
        return 1;

    /* Do inf elimination */
    for (i = 0; i < n - 1; i++) {
        if (warm) {
            m = multi[i];
        } else {
            m = -inf[i] / diag[i];
            if (multi)
                multi[i] = m;
            diag[i + 1] += m * sup[i];
        }
        ind[i + 1] += m * ind[i];
    }

    sol[n - 1] = ind[n - 1] / diag[n - 1];

    for (i = n - 2; i >= 0; i--)
        sol[i] = (ind[i] - sup[i] * sol[i + 1]) / diag[i];

    return 0;
}


// build splines (linear or cubic) from a set of points
// consider cubic splines for smooth angles and linear splines for sharp angles

void buildSpline(vector<v3> &points, vector<v3> &normals, SPolygon &spoly) {
    int i, j, np, nsp;
    v3 d, d_previous, previous_point;
    vector<v3> y_data, spline_normals;
    float d_norm, d_previous_norm, a, ignored;
    Spline spline, last_spline;
    bool cubic, start, last_cubic;

    np = (int) points.size();
    if (np < 2) // nothing to do with just one point
        return;

    // we start by assuming a linear spline
    cubic = false;
    last_cubic = true; // this will avoid consecutive linear splines
    start = true;
    ignored = 0.0f;

    //	d_previous = points[1] - points[0];
    //	d_previous_norm = d_previous.norm();
    //	arc_length = d_previous_norm;

    i = 0;
    do {
        d_previous = points[i + 1] - points[i];
        d_previous_norm = d_previous.norm();
        i++; // avoid leading equal points
    } while (d_previous_norm <= 0.00001);

    // starting spline
    spline.clear();
    // first point in spline
    spline.push_back(points[i - 1], normals[i - 1]);
    previous_point = points[i - 1];

    // go for all remaining points
    for (; i < np - 1; i++) {
        d = points[i + 1] - points[i];
        ignored += d.norm();
        if (ignored < 2 * printer_char.extrusion_height) // line segment too short
            continue;
        ignored = 0.0f;
        a = acos(abs(mymax(mymin(d.dotproduct(d_previous) / (d_previous_norm * d_norm), 1.0f), -1.0f)));
        if (a > DEG_TO_RAD(5)) {
            // we have a linear spline segment
            if (!start && cubic) {
                // we were saving a cubic spline and got a linear segment, so stop the cubic spline
                spline.push_back(points[i], normals[i]);

                spline.set_cubic(cubic);
                if (!spline.is_cubic() && !last_cubic) {// two consecutive
                    last_spline = spoly.back(); // get last spline
                    spoly.pop_back(); // remove it
                    y_data = spline.get_y_data();
                    spline_normals = spline.get_normals();
                    nsp = (int) y_data.size();
                    for (j = 1; j < nsp; j++) { // add current spline to last spline
                        last_spline.push_back(y_data[j], spline_normals[j]);
                    }
                    spoly.push_back(last_spline); // and insert it again
                    last_cubic = spline.is_cubic();
                } else {
                    spoly.push_back(spline);
                    last_cubic = spline.is_cubic();
                }

                spline.clear();
                spline.push_back(points[i], normals[i]);
                previous_point = points[i];
                cubic = false;
                start = true;
            } else {
                // we were saving a linear spline, so continue
                spline.push_back(points[i], normals[i]);
                previous_point = points[i];
                start = false;
                cubic = false;
            }
        } else {
            // we have a cubic spline segment
            if (!start && !cubic) {
                // we were saving a linear spline and got a cubic segment, so stop the linear spline
                spline.set_cubic(cubic);
                spline.push_back(points[i], normals[i]);
                previous_point = points[i];

                if (!spline.is_cubic() && !last_cubic) {// two consecutive
                    last_spline = spoly.back(); // get last spline
                    spoly.pop_back(); // remove it
                    y_data = spline.get_y_data();
                    spline_normals = spline.get_normals();
                    nsp = (int) y_data.size();
                    for (j = 1; j < nsp; j++) { // add current spline to last spline
                        last_spline.push_back(y_data[j], spline_normals[j]);
                    }
                    spoly.push_back(last_spline); // and insert it again
                    last_cubic = spline.is_cubic();
                } else {
                    spoly.push_back(spline);
                    last_cubic = spline.is_cubic();
                }

                spline.clear();
                spline.push_back(points[i], normals[i]);
                previous_point=points[i];
                cubic = true;
                start = true;
            } else {
                // we are saving a cubic spline, so continue
                spline.push_back(points[i], normals[i]);
                previous_point=points[i];
                cubic = true;
                start = false;
            }
        }
        d_previous = d;
        d_previous_norm = d.norm();
    }

    // insert last point to the current spline
    spline.set_cubic(cubic);
    spline.push_back(points[i], normals[i]);

    if (!spline.is_cubic() && !last_cubic) {// two consecutive
        last_spline = spoly.back(); // get last spline
        spoly.pop_back(); // remove it
        y_data = spline.get_y_data();
        spline_normals = spline.get_normals();
        nsp = (int) y_data.size();
        for (j = 1; j < nsp; j++) { // add current spline to last spline
            last_spline.push_back(y_data[j], spline_normals[j]);
        }
        spoly.push_back(last_spline); // and insert it again
    } else {
        spoly.push_back(spline);
    }
}

// build splines (linear or cubic) from a set of points
// consider cubic splines for smooth angles and linear splines for sharp angles
void buildSpline_new(vector<v3> &points, vector<v3> &normals, SPolygon &spoly) {
    int i, j, np, nsp;
    v3 d, d_previous;
    vector<v3> y_data, spline_normals;
    vector<float> x_data;
    float d_norm, d_previous_norm, arc_length, a; //, ignored;
    Spline spline, last_spline;
    bool cubic, start, last_cubic;

    np = (int) points.size();
    if (np < 2) // nothing to do with just one point
        return;

    // we start by assuming a linear spline
    cubic = false;
    last_cubic = true; // this will avoid consecutive linear splines
    start = true;
    //ignored = 0.0f;

    i=0;
    do {
        d_previous = points[i+1] - points[i];
        d_previous_norm = d_previous.norm();
        arc_length = d_previous_norm;
        i++; // avoid leading equal points
    } while(d_previous_norm==0);
    
    // starting spline
    spline.clear();
    // first point in spline
    spline.push_back(0.0f, points[i-1], normals[i-1]);

    // go for all remaining points
    for (; i < np - 1; i++) {
        d = points[i + 1] - points[i];
        d_norm = d.norm();
        if(d_norm==0.0f)
            continue;
        //ignored += d_norm;
        //if (ignored < 2 * printer_char.extrusion_height) // line segment too short
        //    continue;
        //ignored = 0.0f;
        a = acos(abs(mymax(mymin(d.dotproduct(d_previous) / (d_previous_norm * d_norm), 1.0f), -1.0f)));
        if (a > DEG_TO_RAD(5)) {
            // we have a linear spline segment
            if (!start && cubic) {
                // we were saving a cubic spline and got a linear segment, so stop the cubic spline
                spline.push_back(arc_length, points[i], normals[i]);
                /*if (spline.size() > 2) {// we need at least 3 point to form a cubic spline
                        spline.set_cubic(true);
                }
                else {
                        spline.set_cubic(false);
                }*/
                spline.set_cubic(cubic);
                if (!spline.is_cubic() && !last_cubic) {// two consecutive
                    last_spline = spoly.back(); // get last spline
                    spoly.pop_back(); // remove it
                    x_data = spline.get_x_data();
                    y_data = spline.get_y_data();
                    spline_normals = spline.get_normals();
                    nsp = (int) x_data.size();
                    for (j = 1; j < nsp; j++) { // add current spline to last spline
                        last_spline.push_back(x_data[j], y_data[j], spline_normals[j]);
                    }
                    spoly.push_back(last_spline); // and insert it again
                    last_cubic = spline.is_cubic();

                } else {
                    spoly.push_back(spline);
                    last_cubic = spline.is_cubic();
                }

                spline.clear();
                spline.push_back(arc_length, points[i], normals[i]);
                arc_length += d_norm;
                cubic = false;
                start = true;
            } else {
                // we were saving a linear spline, so continue
                spline.push_back(arc_length, points[i], normals[i]);
                arc_length += d_norm;
                start = false;
                cubic = false;
            }
        } else {
            // we have a cubic spline segment
            if (!start && !cubic) {
                // we were saving a linear spline and got a cubic segment, so stop the linear spline
                spline.set_cubic(cubic);
                spline.push_back(arc_length, points[i], normals[i]);

                if (!spline.is_cubic() && !last_cubic) {// two consecutive
                    last_spline = spoly.back(); // get last spline
                    spoly.pop_back(); // remove it
                    x_data = spline.get_x_data();
                    y_data = spline.get_y_data();
                    spline_normals = spline.get_normals();
                    nsp = (int) x_data.size();
                    for (j = 1; j < nsp; j++) { // add current spline to last spline
                        last_spline.push_back(x_data[j], y_data[j], spline_normals[j]);
                    }
                    spoly.push_back(last_spline); // and insert it again
                    last_cubic = spline.is_cubic();

                } else {
                    spoly.push_back(spline);
                    last_cubic = spline.is_cubic();
                }

                spline.clear();
                spline.push_back(arc_length, points[i], normals[i]);
                arc_length += d_norm;
                cubic = true;
                start = true;
            } else {
                // we are saving a cubic spline, so continue
                spline.push_back(arc_length, points[i], normals[i]);
                arc_length += d_norm;
                cubic = true;
                start = false;
            }
        }
        d_previous = d;
        d_previous_norm = d_norm;
    }

    // insert last point to the current spline
    spline.set_cubic(cubic);
    spline.push_back(arc_length, points[i], normals[i]);

    if (!spline.is_cubic() && !last_cubic) {// two consecutive
        last_spline = spoly.back(); // get last spline
        spoly.pop_back(); // remove it
        x_data = spline.get_x_data();
        y_data = spline.get_y_data();
        spline_normals = spline.get_normals();
        nsp = (int) x_data.size();
        for (j = 1; j < nsp; j++) { // add current spline to last spline
            last_spline.push_back(x_data[j], y_data[j], spline_normals[j]);
        }
        spoly.push_back(last_spline); // and insert it again
    } else {
        spoly.push_back(spline);
    }


}

void computeSplineM(Spline &spline) {
    vector<float> x_data = spline.get_x_data();
    vector<v3> y_data = spline.get_y_data();
    double *sup, *inf, *diag, *ind, *Mx, *My, *Mz, *Multi;
    v3 M0 = v3(0.0f, 0.0f, 0.0f), M;
    int np, j;

    // build the m_knots
    // allocate memory
    np = (int) x_data.size(); // number of point in the cubic spline
    sup = (double *) malloc((np - 3) * sizeof (double));
    inf = (double *) malloc((np - 3) * sizeof (double));
    diag = (double *) malloc((np - 2) * sizeof (double));
    ind = (double *) malloc((np - 2) * sizeof (double));
    Mx = (double *) malloc((np - 2) * sizeof (double));
    My = (double *) malloc((np - 2) * sizeof (double));
    Mz = (double *) malloc((np - 2) * sizeof (double));
    Multi = (double *) malloc((np - 2) * sizeof (double));


    // compute m_knots for the cubic splines

    /* line above and below diagonal */
    for (j = 0; j < np - 3; j++)
        sup[j] = inf[j] = (double) (x_data[j + 2] - x_data[j + 1]);

    /* diagonal */
    for (j = 0; j < np - 2; j++)
        diag[j] = 2.0 * (double) (x_data[j + 2] - x_data[j]);

    /* independent term */
    for (j = 0; j < np - 2; j++)
        ind[j] = 6.0 * (double) ((y_data[j + 2].x - y_data[j + 1].x) / (x_data[j + 2] - x_data[j + 1])
            - (y_data[j + 1].x - y_data[j].x) / (x_data[j + 1] - x_data[j]));

    tridiagonal(np - 2, inf, diag, sup, ind, Mx, Multi, false);

    /* independent term */
    for (j = 0; j < np - 2; j++)
        ind[j] = 6.0 * (double) ((y_data[j + 2].y - y_data[j + 1].y) / (x_data[j + 2] - x_data[j + 1])
            - (y_data[j + 1].y - y_data[j].y) / (x_data[j + 1] - x_data[j]));

    tridiagonal(np - 2, inf, diag, sup, ind, My, Multi, true);


    /* independent term */
    for (j = 0; j < np - 2; j++)
        ind[j] = 6.0 * (double) ((y_data[j + 2].z - y_data[j + 1].z) / (x_data[j + 2] - x_data[j + 1])
            - (y_data[j + 1].z - y_data[j].z) / (x_data[j + 1] - x_data[j]));

    tridiagonal(np - 2, inf, diag, sup, ind, Mz, Multi, true);

    // natural cubic spline
    spline.clear_m();
    spline.push_back_m(M0);
    for (j = 0; j < np - 2; j++) {
        M = v3((float) Mx[j], (float) My[j], (float) Mz[j]);
        spline.push_back_m(M);
    }
    spline.push_back_m(M0);

    free(sup);
    free(inf);
    free(diag);
    free(ind);
    free(Mx);
    free(My);
    free(Mz);
    free(Multi);

}

// compute the closest point to p in the line segment p1-p2
std::pair<float, v3> ClosestPoint(v3 p1, v3 p2, v3 p) {
    v3 pminusproj, projection;
    v3 p1minusp2 = p1 - p2;
    float t;
    
    t = p1minusp2.norm(); // ||p1-p2||
    if (t == 0.0f){
        pminusproj = p - p1;
        return make_pair(pminusproj.norm(), p1); // p1 == p2
    }

    t = std::max(0.0f, std::min(1.0f, p1minusp2.dotproduct(p - p2) / t));
    projection = p2 + p1minusp2 * t; // Projection falls on the segment
    pminusproj = p - projection;
    return make_pair(pminusproj.norm(), projection);
}


// compute the closest point to p in the line segment p1-p2
// not working!
std::pair<float, v3> PointPlane(v3 p1, v3 p2, v3 p, v3 v) {
    v3 pminusproj, pminusproj2, projection;
    v3 p1minusp2 = p1 - p2;
    v3 zVector=v3(0.0f,0.0f,1.0f);
    v3 normal = v.crossproduct(zVector);
    float d;
    
    
    d=normal.dotproduct(p);
    
    if (normal.dotproduct(p1minusp2) == 0) {
        // No intersection, the line is parallel to the plane
        pminusproj = p - p1;
        pminusproj2 = p -p2;
        if (pminusproj.norm()<pminusproj2.norm())
            return make_pair(pminusproj.norm(), p1);
        else
            return make_pair(pminusproj2.norm(), p2);
    }
    
    // Compute the X value for the directed line ray intersecting the plane
    float x = (d - normal.dotproduct(p1)) / normal.dotproduct(p1minusp2);

    projection= p1 + p1minusp2*(x/p1minusp2.norm()); //Make sure your ray vector is normalized
    pminusproj = p - projection;
    return make_pair(pminusproj.norm(), projection);
}

// compute the point in the line segment p1-p2 with the same p(2)/p(1)
// point in p1+alpha(p2-p1) with y/x equal do p.y/p.x
std::pair<float, v3> vyvxPoint(v3 p1, v3 p2, v3 p) {
    float num, den, alpha, retalpha;
    v3 projection, projminusp;
    
    den = p2.y*p.x-p2.x*p.y+p1.x*p.y-p1.y*p.x;
    
    if (fabs(den) <=1e-8) {
        // no den, so return +inf
        return make_pair(1e20f, p1);
    }
    
    // Compute the X value for the directed line ray intersecting the plane
    num = p2.x*p.y-p2.y*p.x;
    
    alpha=num/den;
    
//    retalpha=alpha*alpha;
//    if(alpha<0.0f)
//        retalpha+=1.0f;
    
    projection= p1 + (p2-p1)*std::max(0.0f, std::min(1.0f,alpha)); //Make sure your ray vector is normalized
    projminusp=projection-p;
    return make_pair(projminusp.norm(), projection);
}
    
    

void correctStart(vector<Layer> &Layers, v3 entrypoint){
    v3 closestpoint;
    float closestdist;
    std::vector<Polygon> poly;
    int nPoint, sl=Layers.size(), spo, sp, i, j, k;
    vector<v3> Points, NewPoints, Normals, NewNormals;
    std::pair<float, v3> tmpdist;

    //FILE *f = fopen("Closest.m", "w");
    //fprintf(f, "Points=[");
    
    for (i = 0; i < sl; i++) {
        poly = Layers[i].getPoly();
        spo = poly.size();

        for (j = 0; j < spo; j++) {
            if (poly[j].getArea() < 0) {
                closestdist = 1e20; // compute the minimum distance to all splines
                closestpoint = v3(0.0f, 0.0f, 0.0f);
                nPoint = 0;

                Points = poly[j].getPoints();
                Normals = poly[j].getNormals(); // get normals
                sp = Points.size();

                if (i == 0 && j == 0)
                    entrypoint = Points[0];
                else
                    entrypoint.z=Points[0].z;
                
                //if(Points[0].z>72 && Points[0].z<72.8)
                //    printf("Here\n");

//                FILE *ff = fopen("ClosestPoint.m", "w");
//                fprintf(ff, "Entry=[%f,%f,%f];plot3(Entry(1),Entry(2),Entry(3),'*k');hold on;\n", entrypoint.x, entrypoint.y, entrypoint.z);
//                fprintf(ff, "Points=[", entrypoint.x, entrypoint.y, entrypoint.z);

                for (k = 0; k < sp - 1; k++) { // last point is = to first
                    //tmpdist = ClosestPoint(Points[k], Points[k + 1], entrypoint);
                    tmpdist=vyvxPoint(Points[k], Points[k + 1], entrypoint);
                    //tmpdist=PointPlane(Points[k],Points[k+1],entrypoint,v3((float)sqrt(2),(float)sqrt(2),0.0f));
                    //fprintf(ff, "plot3(%f,%f,%f,'or');\n", tmpdist.second.x, tmpdist.second.y, tmpdist.second.z);
                    //fprintf(ff, "%f,%f,%f;", Points[k].x, Points[k].y, Points[k].z);
                    if (tmpdist.first < closestdist) {
                        closestdist = tmpdist.first;
                        closestpoint = tmpdist.second;
                        nPoint = k;
                    }
                }
                //fprintf(ff, "%f,%f,%f;", Points[k].x, Points[k].y, Points[k].z);

                //fclose(ff);

                //fprintf(f, "%f,%f,%f;", closestpoint.x, closestpoint.y, closestpoint.z);
                //fflush(f);


                //fprintf(ff,"];\n%% New spline\nNewPoints=[");
                NewPoints.clear();
                NewNormals.clear();
                NewPoints.push_back(closestpoint); // insert the closet point first
                NewNormals.push_back(Normals[nPoint]); // insert the normal of closet point
                //fprintf(ff, "%f,%f,%f;", closestpoint.x, closestpoint.y, closestpoint.z);
                for (k = (nPoint + 1) % (sp - 1); k != nPoint; k = (k + 1) % (sp - 1)) { // do not account for the last point, since it is the first one
                    NewPoints.push_back(Points[k]); // insert the next closest point
                    NewNormals.push_back(Normals[k]);
                    //fprintf(ff, "%f,%f,%f;", Points[k].x, Points[k].y, Points[k].z);
                    
                }
                // insert last segment (we have k==nPoint)
                NewPoints.push_back(Points[k]); // insert the next closet point
                NewNormals.push_back(Normals[k]);
                
                //fprintf(ff, "%f,%f,%f;", Points[k].x, Points[k].y, Points[k].z);

                NewPoints.push_back(closestpoint); // insert the closet point last to close poly
                NewNormals.push_back(Normals[k]);
                
                //fprintf(ff, "%f,%f,%f];\nhold on;\nplot3(Points(:,1),Points(:,2),Points(:,3))\nplot3(NewPoints(:,1),NewPoints(:,2),NewPoints(:,3))\n", closestpoint.x, closestpoint.y, closestpoint.z);

                //fclose(ff);
                poly[j].setPoints(NewPoints);
                poly[j].setNormals(NewNormals);

                entrypoint = closestpoint;
            }
        }
        Layers[i].setPoly(poly);
    }

    //fprintf(f,"];\nplot3(Points(:,1),Points(:,2),Points(:,3),'.r')");
    //fclose(f);

}


// simplify spline and compute M values for the cubic spline
void computeSpline(SPolygon &spoly) {
    int i, j, ss = (int) spoly.size(), sp; // , saved = 0, total = 0;
    std::vector<Spline> &splines = spoly.getSplines();
    std::pair<vector<v3>, vector <v3>> result;
    float x_end=0.0f;
    vector<float> x_data;
    v3 d;


    for (i = 0; i < ss; i++) {
        // simplify splines before building it
        //saved += (int)splines[i].size();
        //total += (int)splines[i].size();
        result = simplifyWithRDP(splines[i].get_y_data(), splines[i].get_normals(), printer_char.extrusion_height / 100);
        //result = std::make_pair(splines[i].get_x_data(), splines[i].get_y_data());
        //saved -= (int)result.first.size();
        splines[i].set_y_data(result.first);
        splines[i].set_normals(result.second);
        // Compute x_data
        sp=result.first.size();
        x_data.clear();
        x_data.push_back(x_end);
        for (j=0;j<sp-1;j++){
            d=result.first[j+1]-result.first[j];
            x_end+=d.norm();
            x_data.push_back(x_end);
        }
        splines[i].set_x_data(x_data);
        
        if (splines[i].is_cubic()) {
            if (result.first.size() < 3) {
                // we do not have enough points to consider a cubic spline
                // set it as a linear spline
                //printf("Converting cubic to linear spline\n");
                splines[i].set_cubic(false);
            } else {
                computeSplineM(splines[i]);
            }
        }
    }

    //printf("Saved %d points in %d segments with total %d points\n", saved, ss, total);
    spoly.computeArea();
}

v3 cubic(Spline spline, float t, int *segment) {
    int i, np;
    float h;
    vector<float> x_data = spline.get_x_data();
    vector<v3> y_data = spline.get_y_data();
    vector<v3> m_knots = spline.get_m_knots();

    //std::cout << "Cubic called!" << endl;

    np = (int) x_data.size();
    for (i = 1; t >= x_data[i] && i < np - 1; i++);

    h = x_data[i] - x_data[i - 1];
    *segment = i;

    return v3((float) (
            m_knots[i - 1].x / (6 * h) * pow(x_data[i] - t, 3.0) +
            m_knots[i].x / (6 * h) * pow(t - x_data[i - 1], 3.0) +
            (y_data[i - 1].x / h - (m_knots[i - 1].x * h) / 6)*(x_data[i] - t) +
            (y_data[i].x / h - (m_knots[i].x * h) / 6)*(t - x_data[i - 1])), (float) (
            m_knots[i - 1].y / (6 * h) * pow(x_data[i] - t, 3.0) +
            m_knots[i].y / (6 * h) * pow(t - x_data[i - 1], 3.0) +
            (y_data[i - 1].y / h - (m_knots[i - 1].y * h) / 6)*(x_data[i] - t) +
            (y_data[i].y / h - (m_knots[i].y * h) / 6)*(t - x_data[i - 1])), (float) (
            m_knots[i - 1].z / (6 * h) * pow(x_data[i] - t, 3.0) +
            m_knots[i].z / (6 * h) * pow(t - x_data[i - 1], 3.0) +
            (y_data[i - 1].z / h - (m_knots[i - 1].z * h) / 6)*(x_data[i] - t) +
            (y_data[i].z / h - (m_knots[i].z * h) / 6)*(t - x_data[i - 1]))
            );
}

v3 linear(Spline spline, float t, int *segment) {
    int i, np;
    float h;
    vector<float> x_data = spline.get_x_data();
    vector<v3> y_data = spline.get_y_data();


    np = (int) x_data.size();
    for (i = 1; t >= x_data[i] && i < np - 1; i++);

    if (t < x_data[i - 1] && t > x_data[i])
        printf("Something wrong\n");



    h = x_data[i] - x_data[i - 1];
    *segment = i;

    return v3(y_data[i - 1].x + ((y_data[i].x - y_data[i - 1].x)*(t - x_data[i - 1])) / h,
            y_data[i - 1].y + ((y_data[i].y - y_data[i - 1].y)*(t - x_data[i - 1])) / h,
            y_data[i - 1].z + ((y_data[i].z - y_data[i - 1].z)*(t - x_data[i - 1])) / h);
}

v3 spline_eval(Spline spline, float t, int *segment) {

    if (spline.is_cubic())
        return cubic(spline, t, segment);
    else
        return linear(spline, t, segment);

}
