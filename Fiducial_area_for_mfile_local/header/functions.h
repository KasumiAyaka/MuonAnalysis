#pragma once

#include <fstream>
#include <iostream>
#include <ios>
#include <iomanip> 
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <cassert>

namespace matrix_3D {
	class matrix_33 {
	public:
		double val[3][3];

		matrix_33();
		matrix_33(int axis, double rot_angle);
		void matrix_multiplication(matrix_33 other);
		void Print();
	};

	class vector_3D {
	public:
		double x, y, z;
		void matrix_multiplication(matrix_33 other);
	};

	vector_3D addition(vector_3D v0, vector_3D v1);
	double dot(vector_3D v0, vector_3D v1);
	vector_3D const_multiple(vector_3D v, double val);

	double minimum_distance(vector_3D pos0, vector_3D pos1, vector_3D dir0, vector_3D dir1);
	double minimum_distance(vector_3D pos0, vector_3D pos1, vector_3D dir0, vector_3D dir1, double point[3]);
	double minimum_distance(vector_3D pos0, vector_3D pos1, vector_3D dir0, vector_3D dir1, double& extra0, double& extra1);
	double minimum_distance(vector_3D pos0, vector_3D pos1, vector_3D dir0, vector_3D dir1, double z_range[2], double extra[2]);

	double distance(vector_3D v0, vector_3D v1);
	double opening_angle(vector_3D dir0, vector_3D dir1);
	double inpact_parameter(vector_3D pos0, vector_3D dir0, vector_3D pos1);

}
