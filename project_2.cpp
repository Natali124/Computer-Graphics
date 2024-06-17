#include <algorithm>
#include <random>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>
#include <stdio.h>
#include <nanoflann.hpp>
#include "lbfgs.h"
#include "stb_image_write.h"
#include "save_frame.h"


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

    Polygon clip_by_edge(Polygon &cell, const Vector &u, const Vector &v) {
        Polygon outPolygon;
        auto vertices = cell.vertices;
        size_t N = vertices.size();

        Vector Normal(v[1] - u[1], u[0] - v[0]);

        for (size_t i = 0; i < N; ++i) {
            Vector curVertex = vertices[i];
            Vector prevVertex = vertices[(i + N - 1) % N];

            bool current_inside = ((v[0] - u[0]) * (curVertex[1] - u[1]) > (v[1] - u[1]) * (curVertex[0] - u[0]));
            bool previous_inside = ((v[0] - u[0]) * (prevVertex[1] - u[1]) > (v[1] - u[1]) * (prevVertex[0] - u[0]));

            if (current_inside) {
                if (!previous_inside) {
                    double t = dot(Normal, (u - prevVertex) / dot(Normal, curVertex - prevVertex));
                    outPolygon.vertices.push_back(prevVertex + t*(curVertex - prevVertex));
                }
                outPolygon.vertices.push_back(curVertex);
            } else if (previous_inside) {
                double t = dot(Normal, (u - prevVertex) / dot(Normal, curVertex - prevVertex));
                outPolygon.vertices.push_back(prevVertex + t*(curVertex - prevVertex));
            }
        }

        return outPolygon;
    }

    Polygon clip_by_polygon(const Polygon &cell, const Polygon &polygon){
        Polygon clipped_polygon = cell;
        
        for (size_t i = 0; i < polygon.vertices.size(); i++)
            clipped_polygon = clip_by_edge(clipped_polygon, polygon.vertices[i], polygon.vertices[(i + 1) % polygon.vertices.size()]);

        return clipped_polygon;
    }

    Polygon clip_by_circle(Vector& center, Polygon &cell, double radius, int nvert = 100){
        Polygon circle;
        double angle_step = 2 * M_PI / nvert;

        for (int i = 0; i < nvert; i++) {
            double angle = i * angle_step;
            double x = center[0] + radius * std::cos(angle);
            double y = center[1] + radius * std::sin(angle);
            circle.vertices.push_back(Vector(x, y));
        }

        return clip_by_polygon(cell, circle);
    }

    void compute_fluid(){

        cells.clear();

        Polygon square;
        square.vertices.push_back(Vector(0,0,0));
        square.vertices.push_back(Vector(1,0,0));
        square.vertices.push_back(Vector(1,1,0));
        square.vertices.push_back(Vector(0,1,0));
        cells.resize(points.size());

        for (int i = 0; i < points.size(); i++){
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

        for (int i=0; i < points.size(); i++){
            cells[i] = clip_by_circle(points[i], cells[i], std::sqrt(weights[i] - weights[points.size()]));
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
    double fluid_volume;
    SemiDiscreteOT():diagram(), fluid_volume(0){}

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

    static lbfgsfloatval_t _evaluate_fluid(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
    ){
        SemiDiscreteOT* ot = reinterpret_cast<SemiDiscreteOT*>(instance);
        lbfgsfloatval_t fx = 0.0;
        size_t N = n-1;
        
        for (int i = 0; i < n; i++){
            ot->diagram.weights[i] = x[i]; // update weights by N weights for water particles + 1 weight for air
        }

        ot->diagram.compute_fluid(); // this gives N cells for water particles
        
        // Note: N points, N cells, N+1 weights, N+1 lambdas
        auto points_to_use = ot->diagram.points;
        auto cells_to_use = ot->diagram.cells;

        double estimated_fluid_volume = 0.0;

        for (int i = 0; i < points_to_use.size(); i++) {
            Vector y_i = points_to_use[i];
            Polygon cell_i = cells_to_use[i];

            double cell_area = compute_area(cell_i);
            estimated_fluid_volume += cell_area;

            fx += compute_integral_dist(cell_i, y_i) - x[i]*cell_area + ot->fluid_volume*x[i]/(double)N;
            g[i] = cell_area - ot->fluid_volume/double(N);
        }
        double desired_air_volume = 1-ot->fluid_volume;
        double estimated_air_volume = 1.0 - estimated_fluid_volume;
        fx += x[n-1]*(desired_air_volume - estimated_air_volume);
        g[n - 1] = estimated_air_volume - desired_air_volume / double(N);

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
        double objectivefct = -1;
        
        auto ret = lbfgs(N, m_x, &objectivefct, _evaluate, _progress, this, NULL);
        printf("L-BFGS optimization terminated with status code = %d\n", ret);

        lbfgs_free(m_x);
    }


    // Returns centroid of a passed polygon
    Vector find_centroid(Polygon& pol) {
        double C_x = 0.0;
        double C_y = 0.0;
        size_t N = pol.vertices.size();
        for (int i = 0; i < N; i++) {
            C_x += (pol.vertices[i][0] + pol.vertices[(i+1)%N][0]) * (pol.vertices[i][0]*pol.vertices[(i+1)%N][1] - pol.vertices[(i+1)%N][0]*pol.vertices[i][1]);
            C_y += (pol.vertices[i][1] + pol.vertices[(i+1)%N][1]) * (pol.vertices[i][0]*pol.vertices[(i+1)%N][1] - pol.vertices[(i+1)%N][0]*pol.vertices[i][1]);
        }
        double A = compute_area(pol);
        Vector centroid = Vector(C_x, C_y, 0.0) / (6.0*A);
        return centroid;
    }

    /*
    Gallouet Merigot scheme. Pass positions, velocity, and mass and updates positions and velocity
    */
   // void GallouetMerigot(std::vector<Vector>&positions, std::vector<Vector>&velocity, const std::vector<double> &mass){
    void GallouetMerigot(std::vector<Vector>&velocity, const std::vector<double> &mass){
        // optimize weights of W of cells of all particles
        // To do: optimize weights
        int N = diagram.points.size();
        lbfgsfloatval_t *m_x = lbfgs_malloc(N + 1);
        for (int i = 0; i < N + 1; i++) {
            m_x[i] = 1.0/double(5); // initialize
        }
        double objectivefct = -1;

        auto ret = lbfgs(N + 1, m_x, &objectivefct, _evaluate_fluid, _progress, this, NULL);
        printf("L-BFGS optimization terminated with status code = %d\n", ret);

        lbfgs_free(m_x);

        // // Apply forces
        double dt = 0.02;
        double eps = 0.004;
        Vector g(0.0, -20, 0.0); // downwards
        for (int i=0; i<N; i++){
            Vector F_spring = (find_centroid(diagram.cells[i]) - diagram.points[i])/(eps*eps);
            Vector F = F_spring + mass[i]*g;
            // printf("F is %f, %f, %f\n", F[0], F[1], F[2]);
            Vector v_new = velocity[i] + (dt/mass[i])*F;

            //printf("v_new is %f, %f, %f\n", v_new[0], v_new[1], v_new[2]);

            Vector x_new = diagram.points[i] + dt*velocity[i];

            bool bounce = false;
            if (x_new[0] <= 0.0 || x_new[0] >= 1){
                velocity[i] = Vector(0.0, v_new[1], 0);
                bounce = true;
            }
            if (x_new[1] <= 0.0 || x_new[1] >= 1){
                velocity[i] = Vector(v_new[0], 0.0, 0);
                bounce = true;
            }
            if (!bounce){
                diagram.points[i] = x_new;
                velocity[i] = v_new;
            }
        }
    }
};


int main(){
    
    SemiDiscreteOT ot;

    Voronoi vor;
    int N = 1000;
    vor.points.resize(N);
    vor.weights.resize(N);
    vor.lambdas.resize(N);

    for (int i=0; i<N; i++){
        vor.points[i] = Vector(rand()/double(RAND_MAX), rand()/double(RAND_MAX));
        vor.weights[i] = rand()/double(RAND_MAX);
        vor.lambdas[i] = 1.0/double(N);
    }

    // vor.compute();


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

    // Voronoi vor;
    // Polygon pol;

    // Vector init(0.1, 0.1);
    // Vector v1 = Vector(0, 0, 0) + init;
    // Vector v2 = Vector(0.5, 0, 0) + init;
    // Vector v3 = Vector(0.5, 0.5, 0) + init;
    // Vector v4 = Vector(0, 0.5, 0) + init;

    // pol.vertices.push_back(v1);
    // pol.vertices.push_back(v2);
    // pol.vertices.push_back(v3);
    // pol.vertices.push_back(v4);

    // Vector center = Vector(0.5, 0.5, 0) + init;

    // auto cell = vor.clip_by_circle(center, pol, 0.5);
    // std::vector<Polygon>cells;
    // cells.push_back(cell);

    // std::cout << cell.vertices.size();

    // save_svg(cells, "clipping_test.svg");

    /*
    
    FLUID TEST
    
    */
    // SemiDiscreteOT ot;
    // Voronoi vor;
    // int N = 10; // number of fluid cells
    // vor.points.resize(N);
    // vor.weights.resize(N + 1);

    // for (int i=0; i<N; i++){
    //     vor.points[i] = Vector(rand()/double(RAND_MAX), rand()/double(RAND_MAX));
    // }
    // ot.diagram = vor;
    // ot.fluid_volume = 0.2;

    // std::vector<Vector>velocity(N);
    // std::vector<double>mass(N);
    // for (int i=0; i<N; i++){
    //     velocity[i] = Vector(0.0, 0.0, 0.0);
    //     mass[i] = 200.0;
    // }

    // for (int iter = 0; iter < 100; iter++){
    //     ot.GallouetMerigot(velocity, mass);
    //     printf("Iter: %d\n", iter);
    //     save_frame(ot.diagram.cells, "gif_pngs/goodepstest", iter);
    // }
    return 0;
}