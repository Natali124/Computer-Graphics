#include <algorithm>

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
	// Constructor
	Sphere(const Vector &C, double R, const Vector &albedo, bool mirror): C(C), R(R), albedo(albedo), mirror(mirror) {};
	
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

	// void getColor(int bestind, double I, Vector &L, Vector &P, Vector &N, Vector &Color){
	// 	Vector albedo = objects[bestind].albedo;
	// 	albedo.normalize();
	// 	Color = I/(4*PI*(L-P).norm2())*(albedo/PI)*dot(N, (L-P)/(L-P).norm());
	// }

	// void print(Vector &v){
	// 				printf("%f, %f, %f\n", v[0], v[1], v[2]);
	// 			}

	Vector getColor(const Ray& ray, int ray_depth){
		if (ray_depth < 0) return Vector(0, 0, 0);

		Vector P;
		Vector N;
		int sphere_id = intersect(ray, P, N);
		if (sphere_id != -1){
			if (objects[sphere_id].mirror) {
				Vector reflected_angle = ray.u - 2*dot(ray.u, N)*N;
				reflected_angle.normalize();
				Ray reflected_ray = Ray(P + N*epsilon, reflected_angle);
				return getColor(reflected_ray, ray_depth - 1);
			} else {
				Vector albedo = objects[sphere_id].albedo;
				albedo.normalize();
				return I/(4*PI*(L-P).norm2())*(albedo/PI)*dot(N, (L-P)/(L-P).norm());
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

	Sphere S(Vector(0, 0, 0), 10.0, Vector(1, 1, 1), true); // sphere
	Sphere S_more(Vector(-25, 0, 5), 10.0, Vector(1, 1, 1), false); // sphere

	Sphere S_up(Vector(0, 1000, 0), 940.0, Vector(1, 0, 0), false); // sphere
	Sphere S_down(Vector(0, -1000, 0), 990.0, Vector(0, 0, 1), false); // sphere
	Sphere S_left(Vector(0, 0, -1000), 940.0, Vector(5, 102, 8), false); // sphere
	Sphere S_right(Vector(0, 0, 1000), 940.0, Vector(255, 20, 147), false); // sphere

	s.addSphere(S);
	s.addSphere(S_more);
	s.addSphere(S_up);
	s.addSphere(S_down);
	s.addSphere(S_left);
	s.addSphere(S_right);

	std::vector<unsigned char> image(W * H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			double z = -W/(2*tan(fov/2)); // field of view divided by 2
			Vector u(j-W/2+0.5, H/2-i-0.5, z);
			u.normalize();
			Ray r(C, u);
			Vector P;
			Vector N;
			
			Vector color(0,0,0);
			int inter = s.intersect(r, P, N);
			if (inter != -1 && !s.check_shadow(P, N, L)) color = s.getColor(r, 3);

			image[(i * W + j) * 3 + 0] = std::min(std::pow(color[0], 0.454545), 255.0); // gamma correction and capping
			image[(i * W + j) * 3 + 1] = std::min(std::pow(color[1], 0.454545), 255.0); // gamma correction and capping
			image[(i * W + j) * 3 + 2] = std::min(std::pow(color[2], 0.454545), 255.0); // gamma correction and capping
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);
	
	return 0;
}