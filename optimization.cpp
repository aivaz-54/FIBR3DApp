#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <set>
#include "admesh/stl.h"
#include "FIBR3Dapp.hpp"
#include "clipper_utils.hpp"

// our printer characteristics
extern printer printer_char;

stl_file *stl_optimize = NULL;
Mesh complete_mesh;

enum matlab_codes {
    OPEN_FILE = -1, NO_PART = -2, OK = 0
};

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197
#endif

/* setting up to be used on optimization -- used by MATLAB Mex */
enum matlab_codes optim_setup(char *modelFileName, float step, int n_part) {
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
    int modelFileNameLength = 0;
    std::set<std::set<int>> parts_idx;
    std::set<pair<int, int>> connections;
    std::set<int> parts_floor;

    printf("Using file %s with a Slice Step of %f\n", modelFileName, step);

    printer_char.extrusion_height = step;

    modelFileNameLength = (int)strlen(modelFileName);

    if (modelFileNameLength > 4 && !strncmp(modelFileName + modelFileNameLength - 4, ".stl", 4)) {
        // open STL file
        stl_open(&complete_mesh.stl, modelFileName);
        if (complete_mesh.stl.error) {
            fprintf(stderr,"\nError reading STL file.\n\n");
            return OPEN_FILE;
        }
    } else {
        if (modelFileNameLength > 4 && !strncmp(modelFileName + modelFileNameLength - 4, ".obj", 4)) {
            // open and load OBJ file
            stl_initialize(&complete_mesh.stl);

            if (complete_mesh.stl.error) {
                fprintf(stderr,"\nError initializing stl struct.\n\n");
                return OPEN_FILE;
            }


            complete_mesh.stl.fp = fopen(modelFileName, "rb");

            if (complete_mesh.stl.fp == NULL) {
                complete_mesh.stl.error = 1;
                printf("\nError reading OBJ file.\n\n");
                return OPEN_FILE;
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
                            fprintf(stderr, "File can't be read by our simple parser : (Try exporting with other options)\n");
                            fclose(complete_mesh.stl.fp);
                            return OPEN_FILE;
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

                complete_mesh.stl.stats.number_of_facets = (int)vertexIndices.size() / 3;

				if (complete_mesh.stl.stats.number_of_facets <= 0) {
					fprintf(stderr, "File is empty!\n");
					return OPEN_FILE;
				}

                complete_mesh.stl.stats.original_num_facets = complete_mesh.stl.stats.number_of_facets;

                stl_allocate(&complete_mesh.stl);

                // convert vertices to STL
                int it, maxIt, facet_index, first = 1;

                maxIt = (int)vertexIndices.size();
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
            fprintf(stderr, "\nModel file must be a .stl or .obj file.\n\n");
            return OPEN_FILE;
        }
    }

    // remove filename extension
    if (!strncmp(modelFileName + strlen(modelFileName) - 4, ".stl", 4))
        modelFileName[strlen(modelFileName) - 4] = 0;

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

    fprintf(stderr, "Object has %d facets after repair\n", complete_mesh.stl.stats.number_of_facets);

    if (n_part<0) {
        stl_optimize = &complete_mesh.stl;
        return OK;
    }

    split_parts(&complete_mesh.stl, parts_idx, connections, parts_floor);

    fprintf(stderr,"Found %d parts\n", (int) parts_idx.size());

    if (parts_idx.size() > 0 && n_part < parts_idx.size()) {
        std::set<std::set<int>>::iterator it;
        struct Mesh m_part;
        std::set<int>::iterator itp;
        std::set<int> part;
        int first, i;

        it = parts_idx.begin();
        for (i = 0; i < n_part; i++)
            it++;

        if (i > n_part)
            return NO_PART;
        else
            part = *it;

        stl_initialize(&m_part.stl);
        m_part.stl.stats.number_of_facets = (int) part.size();
        stl_allocate(&m_part.stl);

        first = 1;
        for (i = 0, itp = part.begin(); itp != part.end(); itp++) {
            m_part.stl.facet_start[i] = complete_mesh.stl.facet_start[*itp];
            stl_facet_stats(&m_part.stl, m_part.stl.facet_start[i++], first);
            first = 0;
        }

        m_part.normalize();

        fprintf(stderr, "Part %d has %d facets\n", n_part, m_part.stl.stats.number_of_facets);

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

        fprintf(stderr, "Part has %d facets after repair\n", m_part.stl.stats.number_of_facets);

        complete_mesh = m_part;
        stl_optimize = &complete_mesh.stl;
        return OK;

    }


    return NO_PART;
}

// Allow begin called by C

extern "C" int optim_setup_c(char *filename, float step, int part) {
    return optim_setup(filename, step, part);
}

void optim_clean() {
    if (stl_optimize)
        stl_close(stl_optimize);
    stl_optimize = NULL;
}

extern "C" void optim_clean_c() {
    optim_clean();
}

void Newobjfun(float x, float y, double *objvalue) {
    stl_file stl_local;
    int i, j, ls, APs, ss;

    stl_initialize(&stl_local);
    stl_local.stats = stl_optimize->stats;

    stl_local.facet_start = (stl_facet*) malloc(stl_local.stats.number_of_facets * sizeof (stl_facet));
    if (stl_local.facet_start == NULL)
        perror("Slice");
    // recover original facets
    memcpy(stl_local.facet_start, stl_optimize->facet_start, stl_local.stats.number_of_facets * sizeof (stl_facet));


    // rotate along the x and y axis
    //stl_rotate_xy(&stl_local, x, y);
    stl_rotate_x(&stl_local, x);
    stl_rotate_y(&stl_local, y);

    // Generate Slices
    std::vector<std::vector < LineSegment>> slicesWithLineSegments;
    triMeshSlicer(&stl_local, slicesWithLineSegments, printer_char.extrusion_height);

    // Build layers with polygons
    std::vector<Layer> Layers;
    buildLayers(slicesWithLineSegments, Layers);

    std::vector<ClipperLib::Paths> AllPaths;
    ls = (int) Layers.size();
    for (i = 0; i < ls; ++i) {
        AllPaths.push_back(Layer_to_ClipperPaths(Layers[i]));
    }

    *objvalue = 0.0;
    ClipperLib::Paths solution;
    ClipperLib::Clipper clpr;
    APs = (int) AllPaths.size();
    for (i = 0; i < APs - 1; ++i) {
        clpr.AddPaths(AllPaths[i + 1], ClipperLib::ptSubject, true);
        clpr.AddPaths(AllPaths[i], ClipperLib::ptClip, true);
        clpr.Execute(ClipperLib::ctDifference,
                solution, ClipperLib::pftEvenOdd, ClipperLib::pftEvenOdd);
        ss = (int) solution.size();
        for (j = 0; j < ss; ++j) {
            *objvalue += ClipperLib::Area(solution[j]);
        }
        clpr.Clear();
        solution.clear();
    }

}

void MOobjfun(double x, double y, double *SEobjvalue, double *SAobjvalue, double *PTobjvalue, double *se_grad, double *sa_grad, int *n_se, int *n_sa) {
    double dot, abs_dot, tmp_grad[2], nx, ny, nz, deg_to_rad, rx, ry, vx, vy, vz, v_bottom, v_top, v_pos, f_area;
    std::set<int>::iterator it;
    int t, i, at_bottom;


    if (SEobjvalue)*SEobjvalue = 0.0;
    if (SAobjvalue)*SAobjvalue = 0.0;
    // slicing along z
    if (PTobjvalue)*PTobjvalue = 0.0;

    // initialize gradients, if requested
    if (se_grad) {
        se_grad[0] = 0.0;
        se_grad[1] = 0.0;
    }
    if (sa_grad) {
        sa_grad[0] = 0.0;
        sa_grad[1] = 0.0;
    }

    // degrees to radian
    deg_to_rad = ((double) M_PI / (double) 180.0);

    // initialize near to nondifferentiable counters
    if (n_se)*n_se = 0;
    if (n_sa)*n_sa = 0;

    v_bottom = 1e20;
    v_top = -1e20;

    // convert to radians
    rx = deg_to_rad*x;
    ry = deg_to_rad*y;



    if (PTobjvalue || SAobjvalue) { // compute object size after rotation
        for (t = 0; t < stl_optimize->stats.number_of_facets; ++t) { // iterate over all mesh triangles
            for (i = 0; i < 3; ++i) {
                vx = (double) stl_optimize->facet_start[t].vertex[i].x;
                vy = (double) stl_optimize->facet_start[t].vertex[i].y;
                vz = (double) stl_optimize->facet_start[t].vertex[i].z;
                v_pos = -vx * sin(ry) + vy * sin(rx) * cos(ry) + vz * cos(rx) * cos(ry);
                if (v_pos < v_bottom) v_bottom = v_pos;
                if (v_pos > v_top) v_top = v_pos;
            }
        }
    }

    if (PTobjvalue) // compute object size if requested
        *PTobjvalue = v_top - v_bottom;

    if (!SEobjvalue && !SAobjvalue)
        return; // and return if we are done

    for (t = 0; t < stl_optimize->stats.number_of_facets; ++t) { // iterate over all mesh triangles
        // get normal
        nx = (double) stl_optimize->facet_start[t].normal.x;
        ny = (double) stl_optimize->facet_start[t].normal.y;
        nz = (double) stl_optimize->facet_start[t].normal.z;
        // compute triangle area
        f_area = (double) facets_area(stl_optimize->facet_start[t]);
        dot = -nx * sin(ry) + ny * sin(rx) * cos(ry) + nz * cos(rx) * cos(ry);
        if (SEobjvalue) {
            abs_dot = abs(dot);
            if (abs_dot < 1.0) { // do not account for fdot=1, since facet is horizontal
                *SEobjvalue += abs_dot*f_area; // facet is not "almost" horizontal, so staircase effect is present
                if (se_grad && !isnan(se_grad[0])) {// gradient requested
                    if (abs_dot == 0.0) {// non differentiable, return nan
                        se_grad[0] = NAN; // return nan
                        se_grad[1] = NAN;
                    } else {
                        if (n_se && abs_dot < 0.00001)(*n_se)++; // close to nondifferentiable
                        tmp_grad[0] = ny * deg_to_rad * cos(rx) * cos(ry) - nz * deg_to_rad * sin(rx) * cos(ry);
                        tmp_grad[1] = -nx * deg_to_rad * cos(ry) - ny * deg_to_rad * sin(rx) * sin(ry) - nz * deg_to_rad * cos(rx) * sin(ry);
                        tmp_grad[0] *= f_area;
                        tmp_grad[1] *= f_area;
                        if (dot > 0.0) {
                            se_grad[0] += tmp_grad[0];
                            se_grad[1] += tmp_grad[1];
                        } else {
                            se_grad[0] -= tmp_grad[0];
                            se_grad[1] -= tmp_grad[1];
                        }
                    }
                }
            }
        }
        if (SAobjvalue) { // compute object size after rotation
            if (dot <= 0.0) {
                at_bottom = 0;
                for (i = 0; i < 3; ++i) {
                    vx = (double) stl_optimize->facet_start[t].vertex[i].x;
                    vy = (double) stl_optimize->facet_start[t].vertex[i].y;
                    vz = (double) stl_optimize->facet_start[t].vertex[i].z;
                    v_pos = -vx * sin(ry) + vy * sin(rx) * cos(ry) + vz * cos(rx) * cos(ry);
                    if (v_pos <= v_bottom)at_bottom++;
                }
                if (at_bottom < 3) {
                    *SAobjvalue -= dot*f_area;
                    if (sa_grad && !isnan(sa_grad[0])) {// gradient requested	
                        if (dot == 0.0) {// non differentiable, return nan
                            sa_grad[0] = NAN; // return nan
                            sa_grad[1] = NAN;
                        } else {
                            if (n_sa && dot > -0.00001)(*n_sa)++;
                            sa_grad[0] -= (ny * deg_to_rad * cos(rx) * cos(ry) - nz * deg_to_rad * sin(rx) * cos(ry)) * f_area;
                            sa_grad[1] -= (-nx * deg_to_rad * cos(ry) - ny * deg_to_rad * sin(rx) * sin(ry) - nz * deg_to_rad * cos(rx) * sin(ry)) * f_area;
                        }
                    }
                }
            }
        }
    }


    if (SEobjvalue) {
        *SEobjvalue *= (double) (powf(printer_char.extrusion_height, 2.0f) / 2.0f);
        if (se_grad) {
            se_grad[0] *= (double) (powf(printer_char.extrusion_height, 2.0f) / 2.0f);
            se_grad[1] *= (double) (powf(printer_char.extrusion_height, 2.0f) / 2.0f);
        }
    }

    return;
}


extern "C" {

    void objfun(int n, int m, double *x, double *fx) {
        int j;

#pragma omp parallel for
        for (j = 0; j < m; j++) {
            //Newobjfun_c((float)(x[j*n]), (float)(x[j*n + 1]), &(fx[j]));
            MOobjfun((float) (x[j * n]), (float) (x[j * n + 1]), &(fx[j]), NULL, NULL, NULL, NULL, NULL, NULL);
        }


    }
}

