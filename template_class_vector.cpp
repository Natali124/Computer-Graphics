#include <algorithm>

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

# define PI 3.14

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
	// Constructor
	Sphere(const Vector &C, double R, const Vector &albedo): C(C), R(R), albedo(albedo) {};
	
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

	void getColor(double I, Vector &L, Vector &P, Vector &N, Vector &Color){
		Color = I/(4*PI*(L-P).norm2())*(albedo/PI)*dot(N, (L-P)/(L-P).norm());
	}
};

class Scene{
public:
	std::vector<Sphere> objects;
	void addSphere(Sphere s){
		objects.push_back(s);
	}
	// intersect function with same prototype, go through all the objects and check if there is an intersection with a ray.
	// make other intersect return t. find best intersection (smallest)
	bool intersect(const Ray& r, Vector &P, Vector &N){
		double bestt = 1e10;
		for (int i = 0; i < objects.size(); i++){
			objects[i].intersect(r);
		}
	}
};

int main() {
	printf("main function\n");
	int W = 512;
	int H = 512;
	double fov=60*M_PI/180; // 60 degrees
	Sphere S(Vector(0, 0, 0), 10.0, Vector(0, 1, 0)); // sphere
	Vector C(0, 0, 55); // camera
	Vector L(-10,20,40);
	double I = 2e10;

	std::vector<unsigned char> image(W * H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			
			double z = -W/(2*tan(fov/2)); // field of view divided by 2
			Vector u(j-W/2+0.5, H/2-i-0.5, z);
			u.normalize();
			Ray r(C, u);
			Vector P;
			Vector N;
			bool inter = S.intersect(r, P, N);
			Vector color(0,0,0);

			if (inter) {
				S.getColor(I, L, P, N, color);
			}

			image[(i * W + j) * 3 + 0] = std::min(std::pow(color[0], 0.454545), 255.0); // gamma correction and capping
			image[(i * W + j) * 3 + 1] = std::min(std::pow(color[1], 0.454545), 255.0); // gamma correction and capping
			image[(i * W + j) * 3 + 2] = std::min(std::pow(color[2], 0.454545), 255.0); // gamma correction and capping
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);
	
	return 0;
}