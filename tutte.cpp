#include <algorithm>
#include <random>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>
#include <stdio.h>
#include <list>
#include <map>

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

# define PI 3.14159265359
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


class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;  // indices within the uv coordinates array
    int ni, nj, nk;  // indices within the normals array
    int group;       // face group
};

class Edge {
public:
    Edge(){
        vtxA = 0;
        vtxB = 0;
    }
    Edge(int a, int b){
        if (a <= b){
            vtxA = a;
            vtxB = b;
        } else {
            vtxA = b;
            vtxB = a;
        }
    }
    bool operator==(const Edge& other){
        return vtxA == other.vtxA && vtxB == other.vtxB;
    }
    
    int vtxA, vtxB;
};
bool operator< (const Edge& current, const Edge& other){
        if (current.vtxA < other.vtxA){
            return true;
        }
        if (current.vtxA > other.vtxA){
            return false;
        }
        if (current.vtxB < other.vtxB){
            return true;
        }
        else {
            return false;
        }
    }
class TriangleMesh {
public:

  ~TriangleMesh() {}

    TriangleMesh(){};
    
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
    }

    TriangleMesh Tutte() {
        // map edges to triangles

        std::map<Edge, std::vector<int>> edges_to_triangles;
        std::map<int, std::vector<Edge>> all_vtx_to_Edges;


        for (int i=0; i < indices.size(); i++){
            // for each edge there is 2 triangles to that edge
            edges_to_triangles[Edge(indices[i].vtxi, indices[i].vtxj)].push_back(i);
            edges_to_triangles[Edge(indices[i].vtxj, indices[i].vtxk)].push_back(i);
            edges_to_triangles[Edge(indices[i].vtxk, indices[i].vtxi)].push_back(i);

            // push edges opposite to a vertex (for triangle)
            all_vtx_to_Edges[indices[i].vtxi].push_back(Edge(indices[i].vtxj, indices[i].vtxk));
            all_vtx_to_Edges[indices[i].vtxj].push_back(Edge(indices[i].vtxk, indices[i].vtxi));
            all_vtx_to_Edges[indices[i].vtxk].push_back(Edge(indices[i].vtxi, indices[i].vtxj));
        }

        std::vector<Edge> boundaryEdges;
        std::map<int, std::vector<Edge>> vertex_to_edges;
        for (auto it = edges_to_triangles.begin(); it != edges_to_triangles.end(); it ++){
            if (it->second.size() == 1) {
                boundaryEdges.push_back(it->first);
                vertex_to_edges[it->first.vtxA].push_back(it->first);
                vertex_to_edges[it->first.vtxB].push_back(it->first);
            }
        }

        std::vector<Edge> orderedBoundary(boundaryEdges.size());
        orderedBoundary[0] = boundaryEdges[0];
        int veryfirstVertex = boundaryEdges[0].vtxA;
        int otherVertex = boundaryEdges[0].vtxB;
        int edgeIndex = 0;
        while(otherVertex != veryfirstVertex){
            edgeIndex++;
            Edge otherEdge;
            Edge other0 = vertex_to_edges[otherVertex][0];
            Edge other1 = vertex_to_edges[otherVertex][1];
            Edge previousEdge = orderedBoundary[edgeIndex-1]; // I changed boundaryEdges to orderedBoundary since it makes more sense

            if (previousEdge == other0){
                otherEdge = other1;
            } else {
                otherEdge = other0;
            }

            orderedBoundary[edgeIndex] = otherEdge;

            if (otherEdge.vtxA == otherVertex){
                otherVertex = otherEdge.vtxB;
            } else {
                otherVertex = otherEdge.vtxA;
            }
        }
        
        TriangleMesh result = *this;
        // set mapping of the boundary to a circle
        std::vector<bool> isOnBoundary(vertices.size(), false);
        int curVtx = orderedBoundary[0].vtxA;
        
        std::cout << orderedBoundary.size() << std::endl;
        
        for (int i = 0; i<orderedBoundary.size(); i++){
            Edge currentEdge = orderedBoundary[i];

            // std::cout << "A is " << currentEdge.vtxA << std::endl;
            // std::cout << "B is " << currentEdge.vtxB << std::endl;
            
            if (curVtx != currentEdge.vtxB){
                curVtx = currentEdge.vtxB;
            } else {
                curVtx = currentEdge.vtxA;
            }
            //std::cout << curVtx << std::endl;
            isOnBoundary[curVtx] = true;
            double theta = i*2*M_1_PI / (double)orderedBoundary.size();
            result.vertices[curVtx] = Vector(cos(theta), sin(theta), 0);
        }

        int Niter = 1000;
        
        std::vector<Vector> newPositions(vertices.size());
        
        // random initialization
        // for (int i = 0; i < vertices.size(); i++){
        //     if (!isOnBoundary[i]){
        //         result.vertices[i] = Vector(rand()/RAND_MAX, rand()/RAND_MAX, rand()/RAND_MAX);
        //     } else {
        //         std::cout << "is on boundary" << std::endl;
        //     }
        // }

        for (int iter = 0; iter<Niter; iter++){
            for (int i = 0; i<vertices.size(); i++){
                if (isOnBoundary[i]) {
                    newPositions[i] = result.vertices[i]; 
                    continue;
                }
                
                Vector averageVtx(0,0,0);
                for (int e=0; e< all_vtx_to_Edges[i].size(); e++){
                    const Edge &opposite_edge = all_vtx_to_Edges[i][e];
                    averageVtx = averageVtx + result.vertices[opposite_edge.vtxA];
                    averageVtx = averageVtx + result.vertices[opposite_edge.vtxB];
                }
                averageVtx = averageVtx / (2*all_vtx_to_Edges[i].size());
                newPositions[i] = averageVtx;
            }
            result.vertices = newPositions;
        }
        
        return result;
    }

    void writeOBJ(const char* obj) {

		FILE* f = fopen(obj, "w+");

		for (int i = 0; i < vertices.size(); i++) {
			fprintf(f, "v %f %f %f\n", vertices[i][0], vertices[i][1], vertices[i][2]);
		}

		for (int i = 0; i < indices.size(); i++) {
			fprintf(f, "f %u %u %u\n", indices[i].vtxi + 1, indices[i].vtxj + 1, indices[i].vtxk + 1);
		}

		fclose(f);

	}
	
	std::vector<unsigned char*> textures;
	std::vector<int> textW, textH;
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;

};

int main(){

    TriangleMesh T;
	T.readOBJ("goethe.obj");

    TriangleMesh result;
    result = T.Tutte();
    result.writeOBJ("resultTutte.obj");

    return 0;
}