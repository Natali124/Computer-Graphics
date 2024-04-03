#include <algorithm>
#include <random>
#include <cstdlib>
#include <ctime>
#include <iostream>

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

class Sphere {
public:
	// Data fields
	const Vector C;
	double R;
	Vector albedo;
	bool mirror;
	double refindex;
	// Constructor
	Sphere(const Vector &C, double R, const Vector &albedo, bool mirror): C(C), R(R), albedo(albedo), mirror(mirror), refindex(0) {};
	
	// Functions
	double intersect(const Ray& r, Vector &P, Vector &N){
		
		double delta = dot(r.u, r.O - C)*dot(r.u, r.O - C) - ((r.O - C).norm2() - R*R);
		if (delta < 0) return -1;

		double t1 = dot(r.u, C-r.O) - sqrt(delta);
		double t2 = dot(r.u, C-r.O) + sqrt(delta);

		if (t2 < 0) return -1; // No intersection
		if (t1 < 0) {
			t1 = t2;
		}
		// Intersection happened
		P = r.O + t1*r.u;
		N = P - C;
		N.normalize();

		return t1;
	}
};

class Scene{
public:
	std::vector<Sphere> objects;
	double I;
	Vector L;

	Scene(double intensity, Vector &light): objects(), I(intensity), L(light) {};
	
	void addSphere(Sphere s){
		objects.push_back(s);
	}
	
	int intersect(const Ray& r, Vector &P, Vector &N){
		// intersect function with same prototype, go through all the objects and check if there is an intersection with a ray.
		// make other intersect return t. find best intersection (smallest)
		double bestt = 1e10;
		int bestind = -1;

		for (int i = 0; i < (int)sizeof(objects); i++){
			double tmp;
			tmp = objects[i].intersect(r, P, N);
			if (tmp == -1) continue;
			
			if (tmp <= bestt){
				bestt = tmp;
				bestind = i;
			}
		}
		objects[bestind].intersect(r, P, N);
		return bestind;
	}

	bool check_shadow(Vector &P, Vector &N, Vector &L){
		Vector angle = L - P;
		angle.normalize();
		Ray r = Ray(P + N*epsilon, angle);
		Vector P_prime;
		Vector N_prime;

		if (intersect(r, P_prime, N_prime) == -1){
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
		
		int sphere_id = intersect(ray, P, N);
		
		if (sphere_id != -1){
			if (objects[sphere_id].mirror) { // reflection
				Vector reflected_angle = ray.u - 2*dot(ray.u, N)*N;
				reflected_angle.normalize();
				Ray reflected_ray = Ray(P + N*epsilon, reflected_angle);
				return getColor(reflected_ray, ray_depth - 1);
			} 
			else if (objects[sphere_id].refindex != 0) // refraction
			{
				double n2;
				Vector N_saved = N;

				n2 = objects[sphere_id].refindex;
				
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
				Vector albedo = objects[sphere_id].albedo;
				albedo.normalize();
				Vector direct_light = I/(4*PI*(L-P).norm2())*(albedo/PI)*dot(N, (L-P)/(L-P).norm());
				if (check_shadow(P, N, L)){
					direct_light = Vector(0,0,0);
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

	Sphere S(Vector(20, 0, 0), 10.0, Vector(1, 0, 1), false); // sphere
	Sphere S_inside(Vector(20, 0, 0), 9.0, Vector(1, 1, 1), false); // sphere
	Sphere S_more(Vector(0, 0, 0), 10.0, Vector(1, 1, 0), false); // sphere
	Sphere S_3(Vector(-20, 0, 0), 10.0, Vector(1, 1, 1), true); // sphere
	S.refindex = 1.5; // refraction index for glass
	S_inside.refindex = 1/1.5;
	// S_more.refindex = 1.5;

	Sphere S_up(Vector(0, 1000, 0), 940.0, Vector(1, 0, 0), false); // sphere
	Sphere S_down(Vector(0, -1000, 0), 990.0, Vector(0, 0, 1), false); // sphere
	Sphere S_left(Vector(0, 0, -1000), 940.0, Vector(5, 102, 8), false); // sphere
	Sphere S_right(Vector(0, 0, 1000), 940.0, Vector(255, 20, 147), false); // sphere
	Sphere S_side1(Vector(-1000, 0, 0), 940.0, Vector(255, 255, 0), false); // sphere
	Sphere S_side2(Vector(1000, 0, 0), 940.0, Vector(0, 255, 255), false); // sphere

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

	std::vector<unsigned char> image(W * H * 3, 0);
#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			double z = -W/(2*tan(fov/2)); // field of view divided by 2
			
			Vector color(0,0,0);
			int K = 5;

			for (int _i=0; _i<K; _i++){
				
				Vector u(j-W/2+0.5, H/2-i-0.5, z);
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