#ifndef FIBR3Dapp_hpp
#define FIBR3Dapp_hpp

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <set>
#include "admesh/stl.h"
#include <string.h>

#define DEG_TO_RAD(x) (x*0.0174532925199f)
#define PI 3.141592653589793f
#define RAD_TO_DEG(x) (x* 180 / PI)
#define SLICE_EPS 0.01f
#define SLICE_FILL 0.01f
#define SLICE_STEP 0.5f


#ifdef __unix
#define fopen_s(pFile,filename,mode) ((*(pFile))=fopen((filename),(mode)))==NULL
#endif

using namespace std;


const float mymin(const float &a, const float &b);

const float mymax(const float &a, const float &b);


// a point in 3D

struct v3 {
    // a vertice/vector is formed by three coordinates
    float x, y, z;
    // vertice/vector constructor

    v3(float _x = 0, float _y = 0, float _z = 0) : x(_x), y(_y), z(_z) {
    }
    // dotproduct between this vector and a given vector

    float dotproduct(const v3 &v) const {
        return x * v.x + y * v.y + z * v.z;
    }

    v3 crossproduct(v3&v) {
        return v3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }
    // operations with this vertice/vector

    v3& operator-=(const v3 &pt) {
        x -= pt.x;
        y -= pt.y;
        z -= pt.z;
        return *this;
    }
    
    v3 operator-(void) {
        return v3(-x, -y, -z);
    }

    v3 operator-(const v3 &pt) {
        return v3(x - pt.x, y - pt.y, z - pt.z);
    }

    v3 operator+(const v3 &pt) {
        return v3(x + pt.x, y + pt.y, z + pt.z);
    }

    v3 operator/(float a) {
        return v3(x / a, y / a, z / a);
    }

    v3 operator*(float a) {
        return v3(x*a, y*a, z * a);
    }

    float norm() {
        return sqrt(x * x + y * y + z * z);
    }

    void rounding() {
        x = round(x / SLICE_EPS) * SLICE_EPS;
        y = round(y / SLICE_EPS) * SLICE_EPS;
        z = round(z / SLICE_EPS) * SLICE_EPS;
    }
};

// operators between vectors
v3 operator-(const v3 &a, const v3 &b);
v3 operator+(const v3 &a, const v3 &b);

// a line segment between two points with a normal vector
struct LineSegment {
	LineSegment(v3 p0 = v3(), v3 p1 = v3(), v3 n = v3()) { v[0] = p0; v[1] = p1; normal = n; }
	void rounding() { v[0].rounding(); v[1].rounding(); }
	bool equal(LineSegment &ls) {
		if (ls.v[0].x == v[0].x && ls.v[0].y == v[0].y && ls.v[1].x == v[1].x && ls.v[1].y == v[1].y)
			return true;
		else return false;
	}

	// a line segment is formed by two points and a normal
	v3 v[2], normal;
	bool used;
	bool plane;
};


// in case we need to define a plane
//class Plane {
//public:
//	// initialize distance from plane to origin
//	Plane() : mDistance(0) {}
//	// return distance from plane to origin
//	float distance() const { return mDistance; }
//	// return distance from plane to a given point
//	// positive or negative?
//	float distanceToPoint(const v3 &vertex) const { return vertex.dotproduct(mNormal) - mDistance; }
//	// initialize plane
//	void setNormal(v3 normal) { mNormal = normal; }
//	void setDistance(float distance) { mDistance = distance; }
//protected:
//	// a plane is defined by the normal vector and the shortest distance to the origin
//	v3        mNormal;    // normalized Normal-Vector of the plane
//	float   mDistance;  // shortest distance from plane to Origin
//};


// a facet represented by a triangle
struct Triangle {
	// a triangle is formed by three vertices and a normalized normal vector pointing outside the object
	Triangle(v3 n, v3 v0, v3 v1, v3 v2) : normal(n) {
		v[0] = v0; v[1] = v1; v[2] = v2;
		area = 0.5*sqrt(powf((v1.x*v0.y) - (v2.x*v0.y) - (v0.x*v1.y) + (v2.x*v1.y) + (v0.x*v2.y) - (v1.x*v2.y), 2.0f)
			+ powf((v1.x*v0.z) - (v2.x*v0.z) - (v0.x*v1.z) + (v2.x*v1.z) + (v0.x*v2.z) - (v1.x*v2.z), 2.0f)
			+ powf((v1.y*v0.z) - (v2.y*v0.z) - (v0.y*v1.z) + (v2.y*v1.z) + (v0.y*v2.z) - (v1.y*v2.z), 2.0f));
	}
	// operator subtracting a vector to a triangle
	Triangle& operator-=(const v3 &pt) { v[0] -= pt; v[1] -= pt; v[2] -= pt; return *this; }
	float innerproductnorm() const {
		return (sqrt(powf(normal.y, 2.0f) + powf(normal.x, 2.0f)));
	}

	// data
	v3 v[3], normal;
	double area;
};


class TriangleMesh {
public:
	TriangleMesh() : bottomLeftVertex(999999, 999999, 999999), upperRightVertex(-999999, -999999, -999999) {}
	size_t size() const { return mesh.size(); }
	// move 3D Model coordinates to be center around COG(0,0,0)
	void normalize() {
		v3 halfBbox = (upperRightVertex - bottomLeftVertex) / 2.0f;
		v3 start = bottomLeftVertex + halfBbox;
		for (size_t i = 0; i < mesh.size(); ++i) {
			Triangle &triangle = mesh[i];
			triangle -= start;
		}
		bottomLeftVertex = halfBbox*-1.0f;
		upperRightVertex = halfBbox;
	}
	void clean() {
		mesh.clear();
		bottomLeftVertex = v3(999999, 999999, 999999);
		upperRightVertex = v3(-999999, -999999, -999999);
	}
	void push_back(const Triangle &t) {
		mesh.push_back(t);
		for (size_t i = 0; i < 3; ++i) {
			if (t.v[i].x < bottomLeftVertex.x)
				bottomLeftVertex.x = t.v[i].x;
			if (t.v[i].y < bottomLeftVertex.y)
				bottomLeftVertex.y = t.v[i].y;
			if (t.v[i].z < bottomLeftVertex.z)
				bottomLeftVertex.z = t.v[i].z;
			if (t.v[i].x > upperRightVertex.x)
				upperRightVertex.x = t.v[i].x;
			if (t.v[i].y > upperRightVertex.y)
				upperRightVertex.y = t.v[i].y;
			if (t.v[i].z > upperRightVertex.z)
				upperRightVertex.z = t.v[i].z;
		}
	}
	v3 meshBoxSize() const {
		return v3(upperRightVertex.x - bottomLeftVertex.x, upperRightVertex.y - bottomLeftVertex.y, upperRightVertex.z - bottomLeftVertex.z);
	}
	double bottom() const { return bottomLeftVertex.z; }
	std::vector<Triangle>& getMesh() { return mesh; }
	const std::vector<Triangle>& getMesh() const { return mesh; }
	v3 getBottomLeftVertex() const { return bottomLeftVertex; }
	v3 getUpperRightVertex() const { return upperRightVertex; }

protected:
	std::vector<Triangle> mesh;
	v3 bottomLeftVertex, upperRightVertex;
};


// a polygon defined by a list of connected line segments
class Polygon {
public:
	Polygon() {};
	void push_back(v3 const &p, v3 const &n) {
		points.push_back(p);
		normal.push_back(n);
	};
	void clean() { points.clear(); normal.clear(); }
	size_t size() const { return points.size(); }
	std::vector<v3>& getPoints() { return points; }
	const std::vector<v3>& getPoints() const { return points; }
	std::vector<v3>& getNormals() { return normal; }
        void setPoints(vector<v3> pp) { points=pp; };
        void setNormals(vector<v3> nn) { normal=nn; };
	const float getArea() { return area; }
	void correct_normal(v3 const &new_normal) {
		//v3 old_normal = normal.back();
		//normal.pop_back();
		//old_normal = (old_normal + new_normal) / 2.0f;
		//old_normal = old_normal / old_normal.norm();
		//normal.push_back(old_normal);
                normal[normal.size()-1]=new_normal;
	};
	void correct_last_normal(void){
		//v3 first = normal[0];
		//v3 last = normal[normal.size()-1];
		//v3 corrected = (first + last) / 2.0f;
		//corrected = corrected / corrected.norm();
		//normal[0] = corrected;
		//normal[normal.size() - 1] = corrected;
                normal[normal.size() - 1] = normal[0];
	}
	void computeArea() {
		if (points.size() < 3)
			area = 0.0f;

		float a = 0;
		int points_size = (int)points.size();
		for (int i = 0, j = points_size - 1; i < points_size; ++i)
		{
			a += (points[j].x + points[i].x) * (points[j].y - points[i].y);
			j = i;
		}
		area = a * 0.5f;
	}
protected:
	std::vector<v3> points, normal;
	float area;
};

class Spline {
public:

    Spline() {
    };

    void push_back(float x, v3 &y, v3 &n) {
        x_data.push_back(x);
        y_data.push_back(y);
        normals.push_back(n);
    }

    void push_back(v3 &y, v3 &n) {
        x_data.push_back(-1.0f); // undefined
        y_data.push_back(y);
        normals.push_back(n);
    }

    void push_back_m(v3 &M) {
        m_knots.push_back(M);
    }

    void clear() {
        x_data.clear();
        y_data.clear();
        normals.clear();
        m_knots.clear();
    }

    void clear_m() {
        m_knots.clear();
    }

    void set_cubic(bool c) {
        cubic = c;
    }

    bool is_cubic(void) {
        return cubic;
    }

    float get_back() {
        return x_data.back();
    }

    float get_front() {
        return x_data.front();
    }

    vector<float> get_x_data() {
        return x_data;
    }

    vector<v3> get_y_data() {
        return y_data;
    }

    vector<v3> get_normals() {
        return normals;
    }

    vector<v3> get_m_knots() {
        return m_knots;
    }

    void set_x_data(vector<float> x) {
        x_data = x;
    }

    void set_y_data(vector<v3> y) {
        y_data = y;
    }

    void set_normals(vector<v3> n) {
        normals = n;
    }

    size_t size() {
        return x_data.size();
    }

protected:
    bool cubic;
    std::vector<v3> y_data, m_knots, normals;
    std::vector<float> x_data;
};


// a spline polygon is defined by a list of splines

class SPolygon {
public:

    SPolygon() {
    };

    void push_back(Spline &p) {
        splines.push_back(p);
    };

    void pop_back() {
        splines.pop_back();
    };

    Spline& back() {
        return splines.back();
    }

    void clean() {
        splines.clear();
    }

    size_t size() const {
        return splines.size();
    }

    std::vector<Spline>& getSplines() {
        return splines;
    }

    const float getArea() {
        return area;
    }

    void computeArea() {
        // our polygon is defined by a set of Splines
        int i, j;
        int ns = (int) splines.size(), nt;
        vector<v3> points, tmp;


        points.clear();
        points=splines[0].get_y_data();

        for (i = 1; i < ns; i++) {
            tmp = splines[i].get_y_data();
            nt = (int) tmp.size();
            for (j = 1; j < nt; j++) {
                points.push_back(tmp[j]);
            }
        }
        
        //tmp = splines[i].get_y_data();
        //nt = (int) tmp.size();
        //for (j = 0; j < nt; j++) {
        //    points.push_back(tmp[j]);
        //}


        if (points.size() < 3)
            area = 0.0f;

        float a = 0.0f;
        int points_size = (int) points.size();
        for (i = 0, j = points_size - 1; i < points_size; ++i) {
            a += (points[j].x + points[i].x) * (points[j].y - points[i].y);
            j = i;
        }
        area = a * 0.5f;
    }

    void set_bounding_box() {
        // our polygon is defined by a set of Splines
        int i, j;
        int ns = (int) splines.size(), nt;
        vector<v3> tmp;
        
        minX = 1e20;
        maxX = -1e20;
        minY = 1e20;
        maxY = -1e20;
        spline_points.clear();
        
        for (i = 0; i < ns; i++) {
            tmp = splines[i].get_y_data();
            nt = (int) tmp.size();
            for (j = 0; j < nt-1; j++) { // do not copy last, since it is the same as the first
                minX = std::min((double) (tmp[j].x), minX);
                maxX = std::max((double) (tmp[j].x), maxX);
                minY = std::min((double) (tmp[j].y), minY);
                maxY = std::max((double) (tmp[j].y), maxY);
                spline_points.push_back(tmp[j]);
            }
        }
        spline_points.push_back(tmp[j]); // copy the last one to close the polygon
    }

    void get_bounding_box(double *min_X, double *max_X, double *min_Y, double *max_Y) {
        *min_X = minX;
        *max_X = maxX;
        *min_Y = minY;
        *max_Y = maxY;
    }
    
    void clear_points(){
        spline_points.clear();
    }
    
    std::vector<v3> get_spline_points(){
        return spline_points;
    }


protected:
    std::vector<Spline> splines;
    float area;
    double minX, maxX, minY, maxY;
    std::vector<v3> spline_points;
};



class SimpleSpline {
public:

    SimpleSpline() {
    };

    void push_back(v3 &y, v3 &n, SPolygon *p) {
        y_data.push_back(y);
        normals.push_back(n);
        SPoly.push_back(p);
        outer.push_back(false);
    }
    
    void push_back(float x) {
        x_data=x;
    }

    void clear() {
        x_data=0.0f;
        y_data.clear();
        normals.clear();
        outer.clear();
        SPoly.clear();
    }

    void set_outer(bool c, int j) {
        outer[j] = c;
    }

    vector<bool> get_outer(void) {
        return outer;
    }

    vector<v3> get_normals() {
        return normals;
    }
    
    vector<v3> get_y_data() {
        return y_data;
    }
    
    float get_x_data() {
        return x_data;
    }
    
    vector<SPolygon*> get_spoly() {
        return SPoly;
    }

protected:
    std::vector<SPolygon*> SPoly;
    std::vector<bool> outer;
    std::vector<v3> y_data, normals;
    float x_data;
};

class SimpleSpline2 {
public:

    SimpleSpline2() {
    };

    void push_back(v3 &y, v3 &n) {
        y_data.push_back(y);
        normals.push_back(n);
        outer.push_back(false);
    }

    void push_back(SPolygon *p) {
        SPoly=p;
    }

    void clear() {
        y_data.clear();
        normals.clear();
        outer.clear();
        SPoly=NULL;
    }

    void set_outer(bool c, int j) {
        outer[j] = c;
    }

    vector<bool> get_outer(void) {
        return outer;
    }

    vector<v3> get_normals() {
        return normals;
    }
    
    vector<v3> get_y_data() {
        return y_data;
    }
    
    SPolygon* get_spoly() {
        return SPoly;
    }

protected:
    SPolygon* SPoly;
    std::vector<bool> outer;
    std::vector<v3> y_data, normals;
};

// a layer is defined by a set of polygon
class Layer {
public:
	Layer() {};
	void push_back(Polygon &p) {
		poly.push_back(p);
	}
	size_t size() const {
		return poly.size();
	}
	std::vector<Polygon>& getPoly() { return poly; }
	const std::vector<Polygon>& getPoly() const { return poly; }
        void setPoly(std::vector<Polygon> polys) { poly = polys;};
private:
	std::vector<Polygon> poly;
};


// a spline layer is defined by a set of polygon described as splines
class SLayer {
public:
	SLayer() {};
	void push_back(SPolygon &p) {
                p.set_bounding_box();
		spolygon.push_back(p);
	}
	size_t size() const {
		return spolygon.size();
	}
	std::vector<SPolygon>& getSPolygon() { return spolygon; }
        SPolygon *getSPolygon(int j) {return &spolygon[j]; }
private:
	std::vector<SPolygon> spolygon;
};


struct Mesh {
	void normalize() {
		stl_vertex halfsize, start;
		int i, j;

		// normalize object
		halfsize.x = (stl.stats.max.x - stl.stats.min.x) / 2.0f;
		halfsize.y = (stl.stats.max.y - stl.stats.min.y) / 2.0f;
		halfsize.z = (stl.stats.max.z - stl.stats.min.z) / 2.0f;

		start.x = stl.stats.min.x + halfsize.x;
		start.y = stl.stats.min.y + halfsize.y;
		start.z = stl.stats.min.z + halfsize.z;

		for (i = 0; i < stl.stats.number_of_facets; ++i) {
			for (j = 0; j < 3; ++j) {
				stl.facet_start[i].vertex[j].x -= start.x;
				stl.facet_start[i].vertex[j].y -= start.y;
				stl.facet_start[i].vertex[j].z -= start.z;

				stl.facet_start[i].vertex[j].x = round(stl.facet_start[i].vertex[j].x / SLICE_EPS)*SLICE_EPS;
				stl.facet_start[i].vertex[j].y = round(stl.facet_start[i].vertex[j].y / SLICE_EPS)*SLICE_EPS;
				stl.facet_start[i].vertex[j].z = round(stl.facet_start[i].vertex[j].z / SLICE_EPS)*SLICE_EPS;
			}
		}
		stl.stats.max.x = halfsize.x;
		stl.stats.max.y = halfsize.y;
		stl.stats.max.z = halfsize.z;
		stl.stats.min.x = -halfsize.x;
		stl.stats.min.y = -halfsize.y;
		stl.stats.min.z = -halfsize.z;

		stl.stats.size.x = 2.0f*halfsize.x;
		stl.stats.size.y = 2.0f*halfsize.y;
		stl.stats.size.z = 2.0f*halfsize.z;

		stl.stats.bounding_diameter = sqrt(
			stl.stats.size.x * stl.stats.size.x +
			stl.stats.size.y * stl.stats.size.y +
			stl.stats.size.z * stl.stats.size.z
		);

		displacement = v3(start.x, start.y, start.z);
	}
	stl_file stl;
	v3 displacement;
	std::vector<v3> optimal_rotations;
};

struct feature_edge {
	int facet1_number;
	int facet2_number;
	int facet1_neighbor;
	stl_vertex p[2];
	float acos;
	double is_ccw;
	int convex;
};

struct facetstheta {
	float neigh[3];
	stl_vertex start[3];
	stl_vertex end[3];

};

enum GType {G0, G1, SHORT_RETRACT, RETRACT};

struct GCode {
    GType type;
    double speed, x, y, z, B, C, radius2, total_extruded;
};

/**
 * @brief structure that represents the printer
 */
struct printer {
    float extrusion_width = 0.48f/*0.4f*/; /**< @brief extrusion width */
    float extrusion_height = 0.2f; //2.0f/*0.3f*/; /**< @brief extrusion height */
    float filament_radius = 1.75f; /**< @brief filament radius */
    int extrusion_number_of_layers = 3; /**< @brief number of layers to build a surface */
    float extrusion_multiplier = 1.0f; /**< @brief extrusion multiplier */
    float extrusion_speed = 200.0f; /**< @brief speed when extrusion is taking place */
    float brim_size = 10.0f; /**< @brief brim size */
    float nozzle_dia = 0.5f; /**< @brief nozzle diameter */
    float speed = 200.0f; /**< @brief maximum printer speed without extrusion */
    float plate_offset = 67.98f;
    float min_rate;
    float max_rate;
    float min_flow_rate;
    float max_flow_rate;
    v3 min_axis = v3(-137.0f, -277.0f, -50.0f);
    v3 max_axis = v3(437.8f, 243.0f, 285.0f);
    float min_B = DEG_TO_RAD(-45.0f); /**< @brief minimum degree value for B rotation */
    float min_C = DEG_TO_RAD(-180.0f); /**< @brief minimum degree value for C rotation */
    float max_B = DEG_TO_RAD(105.0f); /**< @brief maximum degree value for B rotation*/
    float max_C = DEG_TO_RAD(180.0f); /**< @brief maximum degree value for C rotation*/
    bool export_full_brim = true; /**< @brief export the full brim */
    bool export_brim = true; /**< @brief export the brim */
    bool export_support = true; /**< @brief export the support */
    bool allow_simultaneous = true; /**< @brief allow Level algorithm to consider parts to be built simultaneous */
    float support_height = 16.0f; /**< @brief the support height */
    float max_change_factor = PI / 18.f; /**< @brief maximum change in angles (10 degree per 1 mm) */
    int max_Lf = 1; /**< @brief maximum of sequences allowed in L_f */
    bool simulator = false;
    int max_normals=5; /**< @brief maximum number of normal directions used to get average */
};


void my_stl_free_edges(stl_file *stl);
void my_stl_load_edge_exact(stl_file *stl, stl_hash_edge *edge,
        stl_vertex *a, stl_vertex *b);
void my_stl_initialize_facet_check_exact(stl_file *stl);
void my_insert_hash_edge(stl_file *stl, stl_hash_edge edge);
int my_stl_get_hash_for_edge(int M, stl_hash_edge *edge);
int my_stl_compare_function(stl_hash_edge *edge_a, stl_hash_edge *edge_b);

void find_edge_path(stl_file *stl, std::vector<feature_edge> &feature_edges);
void find_cutting_path(stl_file *stl, std::vector<feature_edge> &feature_edges);
void find_feature_edges(stl_file *stl, std::vector<feature_edge> &feature_edges);

int exportSingleMATLABFormat3D(
        const std::vector<std::vector<LineSegment>> &slicesWithLineSegments,
        char *filename);
int exportFeatureEdgesMATLABFormat3D(
        const std::vector<feature_edge> &feature_edges,
        char *filename_convex,
        char *filename_concave);
int exportLayersMATLABFormat3D(
        std::vector<Layer> &Layers,
        char *filename);


void Newobjfun(float x, float y, double *objvalue);
void MOobjfun(double x, double y, double *SEobjvalue, double *SAobjvalue, double *PTobjvalue, double *se_grad, double *sa_grad);

extern "C" {
    void objfun(int n, int m, double *x, double *fx);
}

void triMeshSlicer(stl_file *stl,
        std::vector<std::vector<LineSegment>> &slicesWithLineSegments, // the result slices
        const float sliceSize);
void buildLayers(// build Layers from a vetor of line segments
        std::vector<std::vector<LineSegment>> &slicesWithLineSegments, // given a vetor of slice line segments
        std::vector<Layer> &Layers);
int exportSLayersMATLABFormat3D(
        std::vector<SLayer> &SLayers,
        char *filename);
float facets_area(stl_facet facet);

void buildSpline(vector<v3> &points, vector<v3> &normals, SPolygon &spoly);
v3 spline_eval(Spline spline, float t, int *segment);
void computeSpline(SPolygon &spoly);
int exportSLayersGCode3Axes(std::vector<SLayer> &SLayers, char *filename);
int exportLayersGCode3Axes(std::vector<Layer> &Layers, char *filename);
v3 cubic(Spline spline, float t, int *segment);
void correctStart(vector<Layer> &Layers, v3 entrypoint);
void computeSplineM(Spline &spline);
int exportSLayersGCode5Axes_shell(std::vector<SLayer> &SLayers, float min_z, float part_size, bool GCode3D, char *filename);
void Slice_inspection(stl_file *stl, float x, float y, char *filename);
void buildSplineLayers(std::vector<Layer> &Layers, std::vector<SLayer> &SLayers);
int split_parts(stl_file *stl, std::set<std::set<int>> &parts_idx, std::set<pair<int, int>> &connections, std::set<int> &parts_floor);

/**
 * @brief check if the part is a shell, i.e. formed by inner and outer polygons
 * @param SLayers
 * @return true if it is considered a shell, false otherwise
 */
bool is_shell(std::vector<SLayer> &SLayers);
bool check_outer_spline(SPolygon spoly, v3 point);


#endif