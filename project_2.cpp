#include <algorithm>
#include <random>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>
#include <stdio.h>
#include <nanoflann.hpp>
#include "lbfgs.h"

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

    // Sutherland-Hodgman algorithm
    Polygon clip_by_bisector(const Polygon& cell, const Vector& P0, const Vector& Pi, const double &w0, const double &wi){
        Polygon result;
        int N = cell.vertices.size();
        Vector offset = (w0-wi)/(2. * (P0-Pi).norm2()) * (Pi-P0);
        Vector M = (P0 + Pi)/2;
        Vector M_prime = M + offset;
        for (int i=0; i<cell.vertices.size(); i++){
            const Vector& A = cell.vertices[i==0 ? (N-1) : i - 1];
            const Vector& B = cell.vertices[i];

            if ((B-P0).norm2() - w0 <= (B - Pi).norm2() - wi) { // B inside
                if ((A-P0).norm2() - w0 > (A-Pi).norm2() - wi){ // A is outside
                    double t = dot(M_prime-A, Pi-P0)/dot(B-A, Pi-P0);
                    Vector P = A + t*(B - A);
                    result.vertices.push_back(P);
                }
                result.vertices.push_back(B);
            }
            else {
                if ((A-P0).norm2() - w0 <= (A-Pi).norm2() - wi) { // A is inside
                    double t = dot(M_prime-A, Pi-P0)/dot(B-A, Pi-P0);
                    Vector P = A + t*(B - A);
                    result.vertices.push_back(P);
                }
            }
        }
        return result;
    }

    void compute(){

        cells.clear();

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
                cell = clip_by_bisector(cell, points[i], points[j], weights[i], weights[j]);
            }
            cells[i] = cell;
        }
    }

    std::vector<Vector> points;
    std::vector<Polygon> cells;
    std::vector<double> weights;
    std::vector<double> lambdas; // For TD7
};

class SemiDiscreteOT{
public: 

    Voronoi diagram;
    SemiDiscreteOT():diagram(){}

    // Integral calculation
    static double compute_integral_dist(const Polygon& cell, const Vector& P_i) {
        // Integral over the cell of y_i of distance from points to y_i
        // We use formula (4.12) in lecture notes
        double integral = 0.0;
        int N = cell.vertices.size();

        for (int k = 1; k <= N; k++){
            double x_k = cell.vertices[k % N][0];
            double x_k_1 = cell.vertices[(k-1)%N][0];
            double y_k = cell.vertices[k%N][1];
            double y_k_1 = cell.vertices[(k-1)%N][1];
            
            double first_term = x_k_1*y_k - x_k*y_k_1;
            double second_term = x_k_1*x_k_1 + x_k_1*x_k + x_k*x_k + y_k_1*y_k_1 + y_k_1*y_k + y_k*y_k;
            double third_term = -4*(P_i[0]*(x_k_1 + x_k) + P_i[1]*(y_k_1 + y_k)) + 6*P_i.norm2();
            integral += first_term * (second_term + third_term);
        }
        integral = integral / 12.0;
        return integral;
    }

    // Integral calculation
    static double compute_area(const Polygon& cell) {
        // Integral over the cell of y_i (we assume f = 1)
        // We use formula: A = 1/2 abs( sum x_i*y_{i+1} - x_{i+1}*y_i)
        double integral = 0.0;
        int N = cell.vertices.size();

        for (int i = 0; i < N; i++) {
            double x_i = cell.vertices[i][0];
            double x_i1 = cell.vertices[(i+1)%N][0];
            double y_i = cell.vertices[i][1];
            double y_i1 = cell.vertices[(i+1)%N][1];

            integral += x_i*y_i1 - x_i1*y_i;
        }
        integral = std::abs(integral)/2.0;
        return integral;
    }

    static lbfgsfloatval_t _evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
    ){
        SemiDiscreteOT* ot = reinterpret_cast<SemiDiscreteOT*>(instance);
        lbfgsfloatval_t fx = 0.0;

        // std::vector<double> weights_to_use = std::vector<double>(x, x + n);
        
        for (int i = 0; i < n; i ++){
            ot->diagram.weights[i] = x[i];
        }

        // ot->diagram.weights = weights_to_use;
        ot->diagram.compute();
        
        auto points_to_use = ot->diagram.points;
        auto lambdas_to_use = ot->diagram.lambdas;
        auto cells_to_use = ot->diagram.cells;

        for (int i = 0; i < points_to_use.size(); i++) {
            Vector y_i = points_to_use[i];
            Polygon cell_i = cells_to_use[i];
            double lambda_i = lambdas_to_use[i];

            double cell_area = compute_area(cell_i);

            fx += compute_integral_dist(cell_i, y_i) - x[i]*cell_area + lambda_i*x[i];
            g[i] = cell_area - lambda_i;
        }
        return -fx;
    }

    static int _progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,   
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
    )
    {
        printf("Iteration %d:\n", k);
        printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
        printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
        printf("\n");
        return 0;
    }

    void optimize(){
        int N = diagram.points.size();
        lbfgsfloatval_t *m_x = lbfgs_malloc(N);
        
        for (int i = 0; i < N; i++) {
            m_x[i] = (double)0.5; // initialize
        }

        for (int i = 0; i < N; i++) {
            if (diagram.points[i][0] < 0.1 || diagram.points[i][0] > 0.9 || diagram.points[i][1] < 0.1 || diagram.points[i][1] > 0.9 ) { // far from center
                m_x[i] = 0.99;
            } else {
                m_x[i] = 1; // more near center
            }
        }

        diagram.weights.resize(N);

        // for (int i = 0; i < N; i++){
        //     diagram.lambdas[i] = double(1)/double(N); // set diagram values to 1/N
        // }

        double objectivefct = -1;
        
        auto ret = lbfgs(N, m_x, &objectivefct, _evaluate, _progress, this, NULL);
        printf("L-BFGS optimization terminated with status code = %d\n", ret);
    }
};


/*
Gallouet Merigot scheme. Pass positions, velocity, and mass and updates positions and velocity
*/
// void GallouetMerigot(std::vector<Vector>&positions, std::vector<Vector>&velocity, const std::vector<double> &mass){
//     // optimize weights of W of cells of all particles
//     std::vector<Polygon>cells;
//     //
//     int N = positions.size();
//     for (int i=0; i<N; i++){
//         Vector F_spring = (centroid(cells[i]) - positions[i]);
//         F
//     }
// }

int main(){
    
    SemiDiscreteOT ot;

    Voronoi vor;
    int N = 100;
    vor.points.resize(N);
    vor.weights.resize(N);
    vor.lambdas.resize(N);

    for (int i=0; i<N; i++){
        vor.points[i] = Vector(rand()/double(RAND_MAX), rand()/double(RAND_MAX));
        vor.weights[i] = rand()/double(RAND_MAX);
        vor.lambdas[i] = 1.0/double(N);
    }

    double sum;
    Vector C (0.5, 0.5, 0); // center
    for (int j = 0; j < N; j++) {
        vor.lambdas[j] = std::exp(-(vor.points[j] - C).norm2()/0.02);
        sum += vor.lambdas[j];
    }
    for (int j = 0; j < N; j++) vor.lambdas[j] /= sum; // have to sum to 1

    // vor.compute();

    ot.diagram = vor;
    ot.optimize();
    
    save_svg(ot.diagram.cells, "voronoi2.svg");

    return 0;
}