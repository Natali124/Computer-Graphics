#include <algorithm>
#include <random>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>
#include <stdio.h>

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

# define PI 3.14
# define epsilon 0.1

class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) {
		data[0] = x;
		data[1] = y;
		data[2] = z;
	}
	double norm2() const {
		return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
	}
	double norm() const {
		return sqrt(norm2());
	}
	void normalize() {
		double n = norm();
		data[0] /= n;
		data[1] /= n;
		data[2] /= n;
	}
	double operator[](int i) const { return data[i]; };
	double& operator[](int i) { return data[i]; };
	double data[3];
};

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, const double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

class Ray {
public:
	// Data fields
	const Vector O;
	const Vector u;
	// Constructor
	Ray(const Vector &O, const Vector &u): O(O), u(u) {};
};

class Geometry {
public:
	bool mirror;
	double refindex;

	virtual double intersect(const Ray& r, Vector &P, Vector &N, Vector &albedo) = 0;

	Geometry() : mirror(false), refindex(0) {};
};

class BoundingBox {
public:
	Vector Bmin;
	Vector Bmax;

	BoundingBox () : Bmin(Vector(0,0,0)), Bmax(Vector(0,0,0)) {};

	void extendBox(Vector& point){
		Bmin = Vector(std::min(point[0], Bmin[0]), std::min(point[1], Bmin[1]), std::min(point[2], Bmin[2]));
		Bmax = Vector(std::max(point[0], Bmax[0]), std::max(point[1], Bmax[1]), std::max(point[2], Bmax[2]));
	}

	bool intersect(const Ray& r) {
		// x plane intersections
		double tx0 = (Bmin[0] - r.O[0]) / r.u[0];
		double tx1 = (Bmax[0] - r.O[0]) / r.u[0];
		double tx_entry = std::min(tx0, tx1);
    	double tx_exit = std::max(tx0, tx1);
		// y plane intersections
		double ty0 = (Bmin[1] - r.O[1]) / r.u[1];
		double ty1 = (Bmax[1] - r.O[1]) / r.u[1];
		double ty_entry = std::min(ty0, ty1);
    	double ty_exit = std::max(ty0, ty1);
		// z plane intersections
		double tz0 = (Bmin[2] - r.O[2]) / r.u[2];
		double tz1 = (Bmax[2] - r.O[2]) / r.u[2];
		double tz_entry = std::min(tz0, tz1);
    	double tz_exit = std::max(tz0, tz1);
		
		double t_total_entry = std::max(tx_entry, std::max(ty_entry, tz_entry));
		double t_total_exit = std::min(tx_exit, std::min(ty_exit, tz_exit));

		return (t_total_entry > 0 && t_total_entry < t_total_exit);
	}
};

class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;  // indices within the uv coordinates array
    int ni, nj, nk;  // indices within the normals array
    int group;       // face group
};

class TriangleMesh : public Geometry {
public:

	BoundingBox boundingBox;
  ~TriangleMesh() {}

    TriangleMesh() : boundingBox() {};

	void computeBoundingBox(BoundingBox& box) {
		box.Bmin = Vector(1e10, 1e10, 1e10);
		box.Bmax = Vector(-1e10, -1e10, -1e10);

		for (auto vertex: vertices) {
			box.extendBox(vertex);
		}
	}
    
    void readOBJ(const char* obj) {
 
        char matfile[255];
        char grp[255];
 
        FILE* f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f)) {
            char line[255];
            if (!fgets(line, 255, f)) break;
 
            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());
 
            if (line[0] == 'u' && line[1] == 's') {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }
 
            if (line[0] == 'v' && line[1] == ' ') {
                Vector vec;
 
                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));
 
                    vertices.push_back(vec);
                    vertexcolors.push_back(col);
 
                } else {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n') {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't') {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f') {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;
 
                char* consumedline = line + 1;
                int offset;
 
                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                    if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                    if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                    indices.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            indices.push_back(t);
                        } else {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                            if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                            if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }
 
                consumedline = consumedline + offset;
 
                while (true) {
                    if (consumedline[0] == '\n') break;
                    if (consumedline[0] == '\0') break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                        if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                        if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    } else {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                            if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                            if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                                if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                                if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;                             
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            } else {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                } else {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }
 
            }
 
        }
        fclose(f);

		computeBoundingBox(boundingBox);
 
    }
 
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
	
	void scale_and_translate(int factor, const Vector &coordinate){
		for (auto index: indices){
			vertices[index.vtxi] = coordinate + vertices[index.vtxi]*factor;
			vertices[index.vtxj] = coordinate + vertices[index.vtxj]*factor;
			vertices[index.vtxk] = coordinate + vertices[index.vtxk]*factor;
		}
	}

	// double intersectTriangle(const Ray& r, TriangleIndices& index, Vector &P, Vector &N, Vector &albedo){

	// 	Vector A = vertices[index.vtxi];
	// 	Vector B = vertices[index.vtxj];
	// 	Vector C = vertices[index.vtxk];

	// 	Vector e1 = B - A;
	// 	Vector e2 = C - A;
	// 	N = cross(e1, e2);

	// 	// std::cout << N[0] << N[1] << N[2] << std:: endl;

	// 	if (dot(r.u, N) == 0){
	// 		return -1;
	// 	}
	// 	double t = dot(A-r.O, N)/dot(r.u, N);

	// 	if (t < epsilon) return -1;

	// 	double beta = dot(e2, cross((A-r.O), r.u))/dot(r.u, N);
	// 	double gamma = -dot(e1, cross((A - r.O), r.u))/dot(r.u, N);
	// 	double alpha = 1-beta-gamma;

	// 	P = alpha*A + beta*B + gamma*C;
	// 	N.normalize();
	// 	if (dot(r.O-P, N) < 0){
	// 	 	N = Vector(0,0,0) - N;
	// 	}
	// 	albedo = Vector(255,255,255); // white by default
	// 	albedo.normalize();

	// 	return t;
	// }

	double intersectTriangle(const Ray& ray, TriangleIndices& index, Vector &P, Vector &N, Vector &albedo){
		Vector v0 = vertices[index.vtxi];
        Vector v1 = vertices[index.vtxj];
        Vector v2 = vertices[index.vtxk];

		Vector edge1 = v1 - v0;
        Vector edge2 = v2 - v0;
		
		Vector h = cross(ray.u, edge2);
        double a = dot(edge1, h);
        if (a > -1e-8 && a < 1e-8) return -1;

		double f = 1.0 / a;
        Vector s = ray.O - v0;
        double u = f * dot(s, h);
        if (u < 0.0 || u > 1.0) return -1;

		Vector q = cross(s, edge1);
        double v = f * dot(ray.u, q);
        if (v < 0.0 || u + v > 1.0) return -1;

		double t_temp = f * dot(edge2, q);
		if (t_temp < 1e-8) return -1;

		P = ray.O + t_temp * ray.u;
		N = cross(edge1, edge2);
		N.normalize();
		albedo = Vector(1, 1, 1);
		albedo.normalize();
		
		return t_temp;
	}

	virtual double intersect(const Ray& r, Vector &P, Vector &N, Vector &albedo){
		double bestt = 1e10;
		bool intersection = false;
		TriangleIndices best_index;

		double t = boundingBox.intersect(r);
		if (t == -1){
			return -1;
		}

		for (auto index: indices){
			double tmp = intersectTriangle(r, index, P, N, albedo);
			if (tmp == -1) continue;
			
			intersection = true;
			if (tmp <= bestt){
				bestt = tmp;
				best_index = index;
			}
		}
		if (intersection){
			double t = intersectTriangle(r, best_index, P, N, albedo);
			return t;
		} else {
			return -1;
		}
	}
};

class Sphere : public Geometry {
public:
	// Data fields
	const Vector C;
	double R;
	Vector albedo;
	// Constructor
	Sphere(const Vector &C, double R, const Vector &albedo): C(C), R(R), albedo(albedo) {};
	
	// Functions
	virtual double intersect(const Ray& r, Vector &P, Vector &N, Vector &albedo_out){
		
		double delta = dot(r.u, r.O - C)*dot(r.u, r.O - C) - ((r.O - C).norm2() - R*R);
		if (delta < 0) return -1;
		
		double t1 = dot(r.u, C-r.O) - sqrt(delta);
		double t2 = dot(r.u, C-r.O) + sqrt(delta);
	
		if (t2 < 0) {
			return -1; // No intersection
		}
		if (t1 < 0) {
			t1 = t2;
		}
		// Intersection happened
		P = r.O + t1*r.u;
		N = P - C;
		albedo_out = albedo;
		albedo_out.normalize();
		N.normalize();

		return t1;
	}
};

class Scene{
public:
	std::vector<Geometry*> objects;
	double I;
	Vector L;

	Scene(double intensity, Vector &light): objects(), I(intensity), L(light) {};
	
	void addSphere(Sphere &s){
		objects.push_back(&s);
	}

	void addTriangleMesh(TriangleMesh m){
		objects.push_back(&m);
	}
	
	int intersect(const Ray& r, Vector &P, Vector &N, Vector &albedo){
		// intersect function with same prototype, go through all the objects and check if there is an intersection with a ray.
		// make other intersect return t. find best intersection (smallest)
		// std::cout << "debuggg" << std::endl;
		double bestt = 1e10;
		int bestind = -1;

		for (int i = 0; i < objects.size(); i++){
			double tmp;
			// std::cout << (int)sizeof(objects) << (objects[i])->refindex << std::endl;
			tmp = objects[i]->intersect(r, P, N, albedo);
			if (tmp == -1){
				continue;
			}
			// std::cout << "debug" << std::endl;
			if (tmp <= bestt){
				bestt = tmp;
				bestind = i;
			}
		}
		if (bestind != -1){
			objects[bestind]->intersect(r, P, N, albedo);
		}
		return bestind;
	}

	bool check_shadow(Vector &P, Vector &N, Vector &albedo, Vector &L){
		Vector angle = L - P;
		angle.normalize();
		Ray r = Ray(P + N*epsilon, angle);
		Vector P_prime;
		Vector N_prime;

		double t = intersect(r, P_prime, N_prime, albedo);

		// std::cout << (P_prime - P).norm2() << " " << (L - P).norm2() << std::endl;
		// std::cout << "N is " << N[0] << N[1] << N[2] << std:: endl;

		if (t == -1){
			return false;
		}
		if ((P_prime - P).norm2() <= (L - P).norm2()){
			return true;
		}
		return false;
	}

	Vector randomDirection(Vector &N){
		std::random_device ran_dev;
		std::mt19937 gen(ran_dev());
		std::uniform_real_distribution<> unif_dis(0.0, 1.0);
		double r1 = unif_dis(gen);
		double r2 = unif_dis(gen);
		
		double x = cos(2*PI*r1)*sqrt(1-r2);
		double y = sin(2*PI*r1)*sqrt(1-r2);
		double z = sqrt(r2);

		Vector T1;

		if (N[0] <= N[1] && N[0] <= N[2]){
			T1 = Vector(0, -N[2], N[1]);
		} else if (N[1] <= N[0] && N[1] <= N[2]){
			T1 = Vector(-N[2], 0, N[0]);
		} else {
			T1 = Vector(-N[1], N[0], 0);
		}
		T1.normalize();
		Vector T2 = cross(N, T1);
		return x*T1 + y*T2 + z*N;
	}

	Vector getColor(const Ray& ray, int ray_depth){
		if (ray_depth < 0) return Vector(0, 0, 0);

		Vector P;
		Vector N;
		Vector albedo;
		
		int sphere_id = intersect(ray, P, N, albedo);
		
		if (sphere_id != -1){
			if (objects[sphere_id]->mirror) { // reflection
				Vector reflected_angle = ray.u - 2*dot(ray.u, N)*N;
				reflected_angle.normalize();
				Ray reflected_ray = Ray(P + N*epsilon, reflected_angle);
				return getColor(reflected_ray, ray_depth - 1);
			} 
			else if (objects[sphere_id]->refindex != 0) // refraction
			{
				double n2;
				Vector N_saved = N;

				n2 = objects[sphere_id]->refindex;
				
				if (dot(ray.u, N) >= 0){
					N = Vector(0,0,0) - N;
					n2 = 1/n2;
				}

				/* Fresnel */
				// Calculate how much refracted and how much reflected
				double k0 = ((1-n2)*(1-n2))/((1+n2)*(1+n2));
				double R = k0 + (1-k0)*pow((1-abs(dot(N, ray.u))), 5);

				std::random_device rd;
				std::mt19937 gen(rd());
				std::uniform_real_distribution<> dis(0.0, 1.0);
				double random_num = dis(gen);

				bool reflect = false;
				if (random_num < R) reflect = true;
				/*Done Frensnel*/

				if (reflect || 1-(1/(n2*n2))*(1-dot(ray.u, N)*dot(ray.u, N)) < 0){
					// same as reflection
					N = N_saved;
					Vector reflected_angle = ray.u - 2*dot(ray.u, N)*N;
					reflected_angle.normalize();
					Ray reflected_ray = Ray(P + N*epsilon, reflected_angle);

					return getColor(reflected_ray, ray_depth - 1);
				}

				double tN = -sqrt(1-(1/(n2*n2))*(1-dot(ray.u, N)*dot(ray.u, N)));
				Vector tTT = 1/n2*(ray.u-dot(ray.u, N)*N);
				
				Vector t = tN*N + tTT;
				
				Ray refracted_ray = Ray(P - N*epsilon, t);
				return getColor(refracted_ray, ray_depth);
			}
			else 
			{
				albedo.normalize();

				//std::cout << "N " << " " << N[0] <<" " << N[1] <<" " << N[2] << std::endl;
				//std::cout << "P " << P[0] <<" " << P[1] << " " <<P[2] << std::endl;
				//std::cout << "albedo " << albedo[0] << albedo[1] << albedo[2] << std::endl;

				Vector direct_light = I/(4*PI*(L-P).norm2())*(albedo/PI)*dot(N, (L-P)/(L-P).norm());

				if (check_shadow(P, N, albedo, L)){
					//std::cout << "shadow detected";
					direct_light = Vector(0,0,0);
				} else {
					// std::cout << "shadow not detected" << std::endl; 
				}
				Vector w = randomDirection(N);
				Ray indirect_ray = Ray(P + N*epsilon, w);
				Vector radiance = getColor(indirect_ray, ray_depth - 1);
				Vector indirect_light = Vector(albedo[0] * radiance[0], albedo[1] * radiance[1], albedo[2] * radiance[2]);

				return direct_light + indirect_light;
			}
		}
		return Vector(0, 0, 0);
	}
};

int main() {
	printf("main function\n");
	int W = 512; // width
	int H = 512; // height

	double fov=60*M_PI/180; // 60 degrees

	Vector C(0, 0, 55); // camera
	Vector L(-10,20,40); // light
	double I = 2e10; // intensity

	Scene s = Scene(I, L);


	TriangleMesh T;
	T.readOBJ("cat.obj");
	// T.scale_and_translate(0.8, Vector(0, -10, 0));
	std::cout<<"mesh loaded"<<std::endl;


	//s.objects.push_back(&T);
	std::cout<<"mesh added to the objects"<<std::endl;

	Sphere S(Vector(20, 0, 0), 10.0, Vector(1, 0, 1)); // sphere
	Sphere S_inside(Vector(20, 0, 0), 9.0, Vector(1, 1, 1)); // sphere
	Sphere S_more(Vector(0, 0, 0), 10.0, Vector(1, 1, 0)); // sphere
	Sphere S_3(Vector(-20, 0, 0), 10.0, Vector(1, 1, 1)); // sphere
	S.refindex = 1.5; // refraction index for glass
	S_inside.refindex = 1/1.5;
	S_3.mirror = true;
	// S_more.refindex = 1.5;

	Sphere S_up(Vector(0, 1000, 0), 940.0, Vector(1, 0, 0)); // sphere
	Sphere S_down(Vector(0, -1000, 0), 990.0, Vector(0, 0, 1)); // sphere
	Sphere S_left(Vector(0, 0, -1000), 940.0, Vector(5, 102, 8)); // sphere
	Sphere S_right(Vector(0, 0, 1000), 940.0, Vector(255, 20, 147)); // sphere
	Sphere S_side1(Vector(-1000, 0, 0), 940.0, Vector(255, 255, 0)); // sphere
	Sphere S_side2(Vector(1000, 0, 0), 940.0, Vector(0, 255, 255)); // sphere

	s.addSphere(S);
	s.addSphere(S_more);
	s.addSphere(S_inside);
	s.addSphere(S_3);
	s.addSphere(S_up);
	s.addSphere(S_down);
	s.addSphere(S_left);
	s.addSphere(S_right);
	s.addSphere(S_side1);
	s.addSphere(S_side2);

	s.objects.push_back(&T);

	std::vector<unsigned char> image(W * H * 3, 0);
#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			double z = -W/(2*tan(fov/2)); // field of view divided by 2
			
			Vector color(0,0,0);
			int K = 1;

			std::random_device ran_dev;
			std::mt19937 gen(ran_dev());
			std::uniform_real_distribution<> unif_dis(0.0, 1.0);

			for (int _i=0; _i<K; _i++){

				double r1 = unif_dis(gen);
				double r2 = unif_dis(gen);
				
				Vector u(j-W/2+0.5 + 0.4*sqrt(-2*log(r1))*cos(2*PI*r2), H/2-i-0.5 + 0.4*sqrt(-2*log(r1))*sin(2*PI*r2), z);
				u.normalize();
				Ray r(C, u);

				color = color + s.getColor(r, 5);
			}
			color = color / K;

			image[(i * W + j) * 3 + 0] = std::min(std::pow(color[0], 0.454545), 255.0); // gamma correction and capping
			image[(i * W + j) * 3 + 1] = std::min(std::pow(color[1], 0.454545), 255.0); // gamma correction and capping
			image[(i * W + j) * 3 + 2] = std::min(std::pow(color[2], 0.454545), 255.0); // gamma correction and capping
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);
	
	return 0;
}