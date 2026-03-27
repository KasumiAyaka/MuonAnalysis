// 2024/09/07
// kasumi
// based on "Check_upstream_base.cpp

#include <iomanip>
#include <cmath>
#include <tuple>
#include "header.h"
#define _USE_MATH_DEFINES
#include <math.h>

double matrix_3D::dot(vector_3D v0, vector_3D v1) {
	return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z;
}
double matrix_3D::distance(vector_3D v0, vector_3D v1) {
	return sqrt(pow(v0.x - v1.x, 2) + pow(v0.y - v1.y, 2) + pow(v0.z - v1.z, 2));

}
matrix_3D::vector_3D matrix_3D::const_multiple(vector_3D v, double val) {
	v.x = v.x * val;
	v.y = v.y * val;
	v.z = v.z * val;
	return v;
}
matrix_3D::vector_3D matrix_3D::addition(vector_3D v0, vector_3D v1) {
	vector_3D vec;
	vec.x = v0.x + v1.x;
	vec.y = v0.y + v1.y;
	vec.z = v0.z + v1.z;
	return vec;
}

double matrix_3D::minimum_distance(vector_3D pos0, vector_3D pos1, vector_3D dir0, vector_3D dir1) {

	double extra0_distance, extra1_distance, delta;
	vector_3D pos;
	pos.x = pos1.x - pos0.x;
	pos.y = pos1.y - pos0.y;
	pos.z = pos1.z - pos0.z;
	//ほぼ平行な場合
	if (opening_angle(dir0, dir1) < 0.0001) {
		extra0_distance = (pos1.z + pos0.z) / 2 - pos0.z;
		extra1_distance = (pos1.z + pos0.z) / 2 - pos1.z;
	}
	else {
		delta = dot(dir0, dir0) * dot(dir1, dir1) - pow(dot(dir0, dir1), 2.);
		extra0_distance = (+1 * dot(pos, dir0) * dot(dir1, dir1) - dot(dir0, dir1) * dot(pos, dir1)) / delta;
		extra1_distance = (-1 * dot(pos, dir1) * dot(dir0, dir0) + dot(dir0, dir1) * dot(pos, dir0)) / delta;
	}
	vector_3D extra0 = addition(pos0, const_multiple(dir0, extra0_distance));
	vector_3D extra1 = addition(pos1, const_multiple(dir1, extra1_distance));

	return distance(extra0, extra1);
}
double matrix_3D::minimum_distance(vector_3D pos0, vector_3D pos1, vector_3D dir0, vector_3D dir1, double point[3]) {
	double extra0_distance, extra1_distance, delta;
	vector_3D pos;
	pos.x = pos1.x - pos0.x;
	pos.y = pos1.y - pos0.y;
	pos.z = pos1.z - pos0.z;
	//ほぼ平行な場合
	if (opening_angle(dir0, dir1) < 0.0001) {
		extra0_distance = (pos1.z + pos0.z) / 2 - pos0.z;
		extra1_distance = (pos1.z + pos0.z) / 2 - pos1.z;
	}
	else {
		delta = dot(dir0, dir0) * dot(dir1, dir1) - pow(dot(dir0, dir1), 2.);
		extra0_distance = (+1 * dot(pos, dir0) * dot(dir1, dir1) - dot(dir0, dir1) * dot(pos, dir1)) / delta;
		extra1_distance = (-1 * dot(pos, dir1) * dot(dir0, dir0) + dot(dir0, dir1) * dot(pos, dir0)) / delta;
	}
	vector_3D extra0 = addition(pos0, const_multiple(dir0, extra0_distance));
	vector_3D extra1 = addition(pos1, const_multiple(dir1, extra1_distance));
	point[0] = (extra0.x + extra1.x) / 2;
	point[1] = (extra0.y + extra1.y) / 2;
	point[2] = (extra0.z + extra1.z) / 2;
	return distance(extra0, extra1);
}
double matrix_3D::minimum_distance(vector_3D pos0, vector_3D pos1, vector_3D dir0, vector_3D dir1, double& extra0_ret, double& extra1_ret) {

	double extra0_distance, extra1_distance, delta;
	vector_3D pos;
	pos.x = pos1.x - pos0.x;
	pos.y = pos1.y - pos0.y;
	pos.z = pos1.z - pos0.z;
	//ほぼ平行な場合
	if (opening_angle(dir0, dir1) < 0.0001) {
		extra0_distance = (pos1.z + pos0.z) / 2 - pos0.z;
		extra1_distance = (pos1.z + pos0.z) / 2 - pos1.z;
	}
	else {
		delta = dot(dir0, dir0) * dot(dir1, dir1) - pow(dot(dir0, dir1), 2.);
		extra0_distance = (+1 * dot(pos, dir0) * dot(dir1, dir1) - dot(dir0, dir1) * dot(pos, dir1)) / delta;
		extra1_distance = (-1 * dot(pos, dir1) * dot(dir0, dir0) + dot(dir0, dir1) * dot(pos, dir0)) / delta;
		if (fabs(extra0_distance) > fabs(pos0.z - pos1.z) || fabs(extra1_distance) > fabs(pos0.z - pos1.z)) {
			double distance0, distance1;
			extra0_distance = 0;
			extra1_distance = pos0.z - pos1.z;
			vector_3D extra0 = addition(pos0, const_multiple(dir0, extra0_distance));
			vector_3D extra1 = addition(pos1, const_multiple(dir1, extra1_distance));
			distance0 = distance(extra0, extra1);
			extra0_distance = pos1.z - pos0.z;
			extra1_distance = 0;
			extra0 = addition(pos0, const_multiple(dir0, extra0_distance));
			extra1 = addition(pos1, const_multiple(dir1, extra1_distance));
			distance1 = distance(extra0, extra1);
			if (distance0 < distance1) {
				extra0_distance = 0;
				extra1_distance = pos0.z - pos1.z;
			}
			else {
				extra0_distance = pos1.z - pos0.z;
				extra1_distance = 0;
			}

		}
	}
	extra0_ret = extra0_distance;
	extra1_ret = extra1_distance;

	vector_3D extra0 = addition(pos0, const_multiple(dir0, extra0_distance));
	vector_3D extra1 = addition(pos1, const_multiple(dir1, extra1_distance));

	return distance(extra0, extra1);
}
double matrix_3D::minimum_distance(vector_3D pos0, vector_3D pos1, vector_3D dir0, vector_3D dir1, double z_range[2], double extra[2]) {
	double extra0_distance, extra1_distance, delta;
	vector_3D pos;
	pos.x = pos1.x - pos0.x;
	pos.y = pos1.y - pos0.y;
	pos.z = pos1.z - pos0.z;
	//ほぼ平行な場合
	if (opening_angle(dir0, dir1) < 0.0001) {
		extra0_distance = (pos1.z + pos0.z) / 2 - pos0.z;
		extra1_distance = (pos1.z + pos0.z) / 2 - pos1.z;
	}
	else {
		delta = dot(dir0, dir0) * dot(dir1, dir1) - pow(dot(dir0, dir1), 2.);
		extra0_distance = (+1 * dot(pos, dir0) * dot(dir1, dir1) - dot(dir0, dir1) * dot(pos, dir1)) / delta;
		extra1_distance = (-1 * dot(pos, dir1) * dot(dir0, dir0) + dot(dir0, dir1) * dot(pos, dir0)) / delta;
	}
	//range[0]:小,range[1]:大
	if (z_range[0] > z_range[1]) {
		double tmp_d = z_range[0];
		z_range[0] = z_range[1];
		z_range[1] = tmp_d;
	}
	if (pos0.z + extra0_distance < z_range[0] || pos1.z + extra1_distance < z_range[0]) {
		extra0_distance = z_range[0] - pos0.z;
		extra1_distance = z_range[0] - pos1.z;
	}
	else if (pos0.z + extra0_distance > z_range[1] || pos1.z + extra1_distance > z_range[1]) {
		extra0_distance = z_range[1] - pos0.z;
		extra1_distance = z_range[1] - pos1.z;
	}

	extra[0] = extra0_distance;
	extra[1] = extra1_distance;
	vector_3D extra0 = addition(pos0, const_multiple(dir0, extra0_distance));
	vector_3D extra1 = addition(pos1, const_multiple(dir1, extra1_distance));

	return distance(extra0, extra1);

}
double matrix_3D::opening_angle(vector_3D dir0, vector_3D dir1) {
	//degreeで返す? radianか?
	double cos = dot(dir0, dir1) / (sqrt(dot(dir0, dir0)) * sqrt(dot(dir1, dir1)));
	if (cos > 1)return 0;
	if (cos < -1)return M_PI;
	return acos(cos);
}
double matrix_3D::inpact_parameter(vector_3D pos0, vector_3D dir0, vector_3D pos1) {
	//点と直線の距離
	double extra = matrix_3D::dot(matrix_3D::addition(pos1, matrix_3D::const_multiple(pos0, -1)), dir0) / matrix_3D::dot(dir0, dir0);
	matrix_3D::vector_3D extra_pos = matrix_3D::addition(pos0, matrix_3D::const_multiple(dir0, extra));

	return distance(pos1, extra_pos);

}
//double impact_parameter() {
//
//}



matrix_3D::matrix_33::matrix_33() {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			val[i][j] = 0;
		}
	}
}
matrix_3D::matrix_33::matrix_33(int axis, double rot_angle) {
	//axis=0:x軸回り
	//axis=1:y軸回り
	//axis=2:z軸回り
	val[axis][axis] = 1;

	val[(axis + 1) % 3][(axis + 1) % 3] = cos(rot_angle);
	val[(axis + 2) % 3][(axis + 2) % 3] = cos(rot_angle);

	val[(axis + 1) % 3][(axis + 2) % 3] = -1 * sin(rot_angle);
	val[(axis + 2) % 3][(axis + 1) % 3] = sin(rot_angle);

	val[axis][(axis + 1) % 3] = 0;
	val[axis][(axis + 2) % 3] = 0;
	val[(axis + 1) % 3][axis] = 0;
	val[(axis + 2) % 3][axis] = 0;

}
void matrix_3D::matrix_33::matrix_multiplication(matrix_3D::matrix_33 left) {
	//行列の積計算
	double calc[3][3];
	/* 行列の積（掛け算） */
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			calc[i][j] = 0;
			for (int k = 0; k < 3; k++) {
				calc[i][j] += left.val[i][k] * val[k][j];
			}
		}
	}
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			val[i][j] = calc[i][j];
		}
	}


}
void matrix_3D::vector_3D::matrix_multiplication(matrix_3D::matrix_33 other) {
	//行列の積計算
	double calc[3];
	/* 行列の積（掛け算） */
	for (int i = 0; i < 3; ++i) {
		calc[i] = other.val[i][0] * x + other.val[i][1] * y + other.val[i][2] * z;
	}
	x = calc[0];
	y = calc[1];
	z = calc[2];
}
void matrix_3D::matrix_33::Print() {
	printf(" %6.4lf %6.4lf %6.4lf\n", val[0][0], val[0][1], val[0][2]);
	printf(" %6.4lf %6.4lf %6.4lf\n", val[1][0], val[1][1], val[1][2]);
	printf(" %6.4lf %6.4lf %6.4lf\n", val[2][0], val[2][1], val[2][2]);
}

class output_format {
public:
	int groupid, chainid, pl, upl;
	int peke, count;
	double dal, dar, md, dz, dz2, oa;
};


struct Key {
	int pl, cid, gid;
};
bool operator<(const Key& lhs, const Key& rhs) {
	return std::tie(lhs.gid, lhs.cid, lhs.pl) < std::tie(rhs.gid, rhs.cid, rhs.pl);
}

struct Btrk0 {
	Key ky;
	int vph, rid;
	double x, y, z, ax, ay;
};
struct Btrk {
	int pl,vph, rid;
	double x, y, z, ax, ay;
};

void SetTracks(std::vector<Btrk>& trk);
void Calc_MinimumDistance(std::vector<Btrk>& trk);
void Calc_MinimumDistance_zrange(std::vector<Btrk>& trk);
void Calc_AngDiff(std::vector<Btrk>& trk);
void Calc_PosDiff(std::vector<Btrk>& trk);


int main(int argc, char** argv) {
	if (argc != 1) {
		fprintf(stderr, "usage:\n");
		exit(1);
	}
	std::vector<Btrk> trk;
	SetTracks(trk);

	std::cout << std::setfill(' ');


	std::cout << " * Calc md " << std::endl;
	Calc_MinimumDistance(trk);
	std::cout << " * Calc andle diff " << std::endl;
	Calc_AngDiff(trk);
	std::cout << " * Calc position diff " << std::endl;
	Calc_PosDiff(trk);
	std::cout << " * Calc md (with Z-Range) " << std::endl;
	Calc_MinimumDistance_zrange(trk);

	std::cout << " Finish process...!" << std::endl;
}


void SetTracks(std::vector<Btrk>&trk) {

	Btrk tmp;

	tmp.vph = 100;

	int i = 0;
	std::cout << "Input track1 info: pl ax ay rawid x y z" << std::endl;
	std::cin >> tmp.pl >> tmp.ax >> tmp.ay >> tmp.rid >>  tmp.x >> tmp.y >> tmp.z;
	trk.push_back(tmp);


	std::cout << "\nInput track2 info: pl ax ay rawid x y z" << std::endl;
	std::cin >> tmp.pl >> tmp.ax >> tmp.ay >> tmp.rid >> tmp.x >> tmp.y >> tmp.z;
	trk.push_back(tmp);

}

void Calc_MinimumDistance(std::vector<Btrk>& trk) {


	auto itr1 = trk.begin();
	auto itr2 = std::next(trk.begin());

	std::cout << " ( PL, rawid )" << std::endl;
	std::cout << std::fixed << std::right
		<< "trk1 :  (" << std::setw(3) << std::setprecision(0) << itr1->pl << ", "
		<< std::setw(12) << std::setprecision(0) << itr1->rid << ")\n"
		<< "trk2 :  (" << std::setw(3) << std::setprecision(0) << itr2->pl << ", "
		<< std::setw(12) << std::setprecision(0) << itr2->rid << ")\n"
		<< std::endl;

	matrix_3D::vector_3D pos0, pos1, dir0, dir1;

	pos0.x = itr1->x;
	pos0.y = itr1->y;
	pos0.z = itr1->z;
	dir0.x = itr1->ax;
	dir0.y = itr1->ay;
	dir0.z = 1;

	pos1.x = itr2->x;
	pos1.y = itr2->y;
	pos1.z = itr2->z;
	dir1.x = itr2->ax;
	dir1.y = itr2->ay;
	dir1.z = 1;


	double oa, md;
	double point[3];
	md = matrix_3D::minimum_distance(pos0, pos1, dir0, dir1, point);
	oa = matrix_3D::opening_angle(dir0, dir1);

	std::cout << std::fixed << std::right
		<< "   md    =  " << std::setw(10) << std::setprecision(1) << md << "\n"
		<< "   oa    =      " << std::setw(10) << std::setprecision(5) << oa << "\n"
		<< "   point = ("
		<< std::setw(10) << std::setprecision(1) << point[0] << ", "
		<< std::setw(10) << std::setprecision(1) << point[1] << ", "
		<< std::setw(10) << std::setprecision(1) << point[2] << ")\n"
		<< std::endl;

}
void Calc_MinimumDistance_zrange(std::vector<Btrk>& trk) {


	auto itr1 = trk.begin();
	auto itr2 = std::next(trk.begin());

	matrix_3D::vector_3D pos0, pos1, dir0, dir1;

	pos0.x = itr1->x;
	pos0.y = itr1->y;
	pos0.z = itr1->z;
	dir0.x = itr1->ax;
	dir0.y = itr1->ay;
	dir0.z = 1;

	pos1.x = itr2->x;
	pos1.y = itr2->y;
	pos1.z = itr2->z;
	dir1.x = itr2->ax;
	dir1.y = itr2->ay;
	dir1.z = 1;


	double oa, md, dz;
	double extra[2], z_range[2];
	z_range[0] = pos0.z;
	z_range[1] = pos1.z;
	oa = matrix_3D::opening_angle(dir0, dir1);
	md = matrix_3D::minimum_distance(pos0, pos1, dir0, dir1, z_range, extra);
	dz = (fabs(extra[0]) + fabs(extra[1])) / 2;

	std::cout << std::fixed << std::right
		<< "   md    = " << std::setw(10) << std::setprecision(1) << md << "\n"
		<< "   oa    =     " << std::setw(10) << std::setprecision(5) << oa << "\n"
		<< "   dz    = " << std::setw(10) << std::setprecision(1) << dz
		<< std::endl;

}

void Calc_AngDiff(std::vector<Btrk>& trk) {

	auto itr1 = trk.begin();
	auto itr2 = std::next(trk.begin());

	double ax, ay, dax, day, dlat, drad, angle;
	// dax,day
	dax = itr2->ax - itr1->ax;
	day = itr2->ay - itr1->ay;
	std::cout << std::fixed << std::right
		<< "   dax   = " << std::setw(10) << std::setprecision(5) << dax << "\n"
		<< "   day   = " << std::setw(10) << std::setprecision(5) << day
		<< std::endl;


	ax = itr2->ax + itr1->ax;
	ay = itr2->ay + itr1->ay;
	angle = sqrt(ax * ax + ay * ay);
	dlat = (-ay * dax + ax * day) / angle;
	drad = (ax * dax + ay * day) / angle;
	if (angle < 0.01) {
		dlat = dax;
		drad = day;
	}
	std::cout << std::fixed << std::right
		<< "   dal   = " << std::setw(10) << std::setprecision(5) << dlat << "\n"
		<< "   dar   = " << std::setw(10) << std::setprecision(5) << drad
		<< std::endl;

	// d(tan)
	double tan1, tan2, dtan;
	tan1 = sqrt(itr1->ax * itr1->ax + itr1->ay * itr1->ay);
	tan2 = sqrt(itr2->ax * itr2->ax + itr2->ay * itr2->ay);

	tan1 = atan(tan1);// theta [rad]
	tan2 = atan(tan2);// theta [rad]
	dtan = tan2 - tan1;// d(theta) [rad]



	std::cout << std::fixed << std::right
		<< "   dtan  = " << std::setw(10) << std::setprecision(5) << dtan << " [rad], ( tan = " << tan(dtan) << " )\n"
		<< std::endl;

}

void Calc_PosDiff(std::vector<Btrk>& trk) {

	auto itr1 = trk.begin();
	auto itr2 = std::next(trk.begin());
	double extra[2], z_range[2];
	matrix_3D::vector_3D pos0, pos1, point;

	pos0.x = -itr1->x;
	pos0.y = -itr1->y;
	pos0.z = -itr1->z;

	pos1.x = itr2->x;
	pos1.y = itr2->y;
	pos1.z = itr2->z;

	point = matrix_3D::addition(pos0, pos1);
	double d = matrix_3D::distance(pos0, pos1);

	std::cout << std::fixed << std::right
		<< "   distance   =  " << std::setw(10) << std::setprecision(1) << d << "\n"
		<< "   (dx,dy,dz) = ("
		<< std::setw(10) << std::setprecision(1) << point.x << ","
		<< std::setw(10) << std::setprecision(1) << point.y << ","
		<< std::setw(10) << std::setprecision(1) << point.z << ")\n"
		<< std::endl;

}

