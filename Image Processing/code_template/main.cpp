#include <iostream>
#include <random>
#include <math.h>
#include <algorithm>

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


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


Vector random_direction(){
	std::random_device rd;
    std::mt19937 gen(rd());
    
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    double r1 = dis(gen);
    double r2 = dis(gen);

	double x = cos(2*M_PI*r1)*sqrt(r2*(1-r1));
	double y = sin(2*M_PI*r1)*sqrt(r2*(1-r2));
	double z = z = 1 - 2*r2;

	return Vector(x, y, z);
}



int main() {

	int W, H, C;
	
	// Load input image
	unsigned char *image_input_raw = stbi_load("8733654151_b9422bb2ec_k.jpg",
                                 &W,
                                 &H,
                                 &C,
                                 STBI_rgb);

	std::vector<double> input_image(W*H*3);
	for (int i=0; i<W*H*3; i++)
		input_image[i] = image_input_raw[i];

	// Load model image
	unsigned char *image_model_raw = stbi_load("redim.jpg",
                                 &W,
                                 &H,
                                 &C,
                                 STBI_rgb);

	std::vector<double> model_image(W*H*3);
	for (int i=0; i<W*H*3; i++)
		model_image[i] = image_model_raw[i];

	// Sliced optimal transport color transfer algorithm
	int K = 100;
	for (int k = 0; k < K; k++) {
		Vector rand_dir = random_direction();
		
		std::vector<std::pair<double, int>> projI;
		std::vector<std::pair<double, int>> projM;

		for (int i = 0; i < H; i++) {
			for (int j = 0; j < W; j++) {
				Vector color_input(input_image[(i*W + j) * 3 + 0], input_image[(i*W + j) * 3 + 1], input_image[(i*W + j) * 3 + 2]);
				Vector color_model(model_image[(i*W + j) * 3 + 0], model_image[(i*W + j) * 3 + 1], model_image[(i*W + j) * 3 + 2]);

				projI.emplace_back(dot(color_input, rand_dir), (i*W + j) * 3);
				projM.emplace_back(dot(color_model, rand_dir), (i*W + j) * 3);
			}
		}

		std::sort(projI.begin(), projI.end());
		std::sort(projM.begin(), projM.end());

		for (size_t i = 0; i < projI.size(); i++) {
			int idxI = projI[i].second;
			int idxM = projM[i].second;
			Vector color_shift = (projM[i].first - projI[i].first) * rand_dir;

			input_image[idxI + 0] += color_shift[0];
			input_image[idxI + 1] += color_shift[1];
			input_image[idxI + 2] += color_shift[2];
		}
	}

	// Result calculation
	
	std::vector<unsigned char> image_result(W*H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

			image_result[(i*W + j) * 3 + 0] = input_image[(i*W+j)*3+0];
			image_result[(i*W + j) * 3 + 1] = input_image[(i*W+j)*3+1];
			image_result[(i*W + j) * 3 + 2] = input_image[(i*W+j)*3+2];
		}
	}
	stbi_write_png("image.png", W, H, 3, &image_result[0], 0);

	return 0;
}
