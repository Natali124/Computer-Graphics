#include <algorithm>
#include <random>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>
#include <stdio.h>
#include <nanoflann.hpp>

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

// if the Polygon class name conflicts with a class in wingdi.h on Windows, use a namespace or change the name
class Polygon {  
public:
    std::vector<Vector> vertices;
};  
 
// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
    void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none") {
        FILE* f = fopen(filename.c_str(), "w+"); 
        fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
        for (int i=0; i<polygons.size(); i++) {
            fprintf(f, "<g>\n");
            fprintf(f, "<polygon points = \""); 
            for (int j = 0; j < polygons[i].vertices.size(); j++) {
                fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
            }
            fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
            fprintf(f, "</g>\n");
        }
        fprintf(f, "</svg>\n");
        fclose(f);
    }


class Voronoi{
public:
    Voronoi() {}

    // sutherland hodgman
    Polygon clip_by_bisector(const Polygon& cell, const Vector& P0, const Vector& Pi, double w0, double w1){
        Polygon result;
        int N = cell.vertices.size();
        Vector offset = (w0-w1)/(2. * (P0-Pi).norm2()) * (Pi-P0);
        Vector M_prime = 
        for (int i=0; i<cell.vertices.size(); i++){
            const Vector& A = cell.vertices[i==0 ? (N-1) : i - 1];
            const Vector& B = cell.vertices[i];

            if ((B-P0).norm2() <= (B - Pi).norm2()) { // B inside
                if ((A-P0).norm2() > (A-Pi).norm2()){ // A is outside
                    Vector M = (P0 + Pi)/2;
                    Vector M_prime = (w0-w1)/(2. * (P0-Pi).norm2()) * (Pi-P0);
                    double t = dot(M-A, Pi-P0)/dot(B-A, Pi-P0);
                    Vector P = A + t*(B - A);
                    result.vertices.push_back(P);
                }
                result.vertices.push_back(B);
            }
            else {
                if ((A-P0).norm2() <= (A-Pi).norm2()) { // A is inside
                    Vector M = (P0 + Pi)/2;
                    double t = dot(M-A, Pi-P0)/dot(B-A, Pi-P0);
                    Vector P = A + t*(B - A);
                    result.vertices.push_back(P);
                }
            }
        }
        return result;
    }


    void compute(){

        Polygon square;
        square.vertices.push_back(Vector(0,0,0));
        square.vertices.push_back(Vector(1,0,0));
        square.vertices.push_back(Vector(1,1,0));
        square.vertices.push_back(Vector(0,1,0));

        cells.resize(points.size());

        for (int i =0; i < points.size(); i++){
            // cell number i
            Polygon cell = square;
            #pragma omp parallel for schedule(dynamic, 1)
            for (int j=0; j<points.size(); j++){
                if (i == j){
                    continue;
                }
                cell = clip_by_bisector(cell, points[i], points[j]);
            }
            cells[i] = cell;
        }
    }

    std::vector<Vector> points;
    std::vector<Polygon> cells;
    std::vector<double> weights;
};
 

int main(){
    std::vector<Polygon> test;
    Polygon pol1, pol2;
    pol1.vertices.push_back(Vector(0,0,0));
    pol1.vertices.push_back(Vector(0,1,0));
    pol1.vertices.push_back(Vector(1,1,0));

    pol2.vertices.push_back(Vector(0.5,0.5,0));
    pol2.vertices.push_back(Vector(0.5,1,0));
    pol2.vertices.push_back(Vector(1,0.3,0));
    pol2.vertices.push_back(Vector(0.3,0.3,0));

    test.push_back(pol1);
    test.push_back(pol2);

    Voronoi vor;
    int N = 1000;
    vor.points.resize(N);
    vor.weights.resize(N);

    for (int i=0; i<N; i++){
        vor.points[i] = Vector(rand()/double(RAND_MAX), rand()/double(RAND_MAX));
        vor.weights[i] = rand()/double(RAND_MAX); // random initialization for testing
    }
    vor.compute();
    
    save_svg(vor.cells, "voronoi.svg");

    return 0;
}