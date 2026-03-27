// suzuki prg
// d(lateral)ÇÃåvéZåÎÇËÇÃÇΩÇﬂèCê≥
#define _CRT_SECURE_NO_WARNINGS
#pragma comment(lib,"FILE_structure.lib")
#pragma comment(lib,"functions.lib")
#include <FILE_structure.hpp>
#include <functions.hpp>


std::vector<netscan::linklet_t> lateral_cut_position(std::vector<netscan::linklet_t> link, double thr);
std::vector<netscan::linklet_t> lateral_cut_angle(std::vector<netscan::linklet_t> link, double thr);
std::vector<netscan::linklet_t> radial_cut_position(std::vector<netscan::linklet_t> link, double thr_intercept, double thr_slope = 0);
std::vector<netscan::linklet_t> radial_cut_angle(std::vector<netscan::linklet_t> link, double thr_intercept, double thr_slope = 0);
void position_difference(netscan::linklet_t link, double& dr, double& dl);
int main(int argc, char** argv) {
	if (argc != 9) {
		fprintf(stderr, "usage:prg in-linklet(bin) out-linklet(bin) lateral-cut-angle lateral-cut-position radial-cut-angle-intercept radial-cut-angle-slope radial-cut-position-intercept radial-cut-position-slope\n");
		exit(1);
	}
	std::string file_in_link = argv[1];
	std::string file_out_link = argv[2];
	double thr_angle_lat = std::stod(argv[3]);
	double thr_position_lat = std::stod(argv[4]);
	double thr_angle_rad_intercept = std::stod(argv[5]);
	double thr_angle_rad_slope = std::stod(argv[6]);
	double thr_position_rad_intercept = std::stod(argv[7]);
	double thr_position_rad_slope = std::stod(argv[8]);

	std::vector<netscan::linklet_t> link;
	netscan::read_linklet_bin(file_in_link, link);

	link = lateral_cut_angle(link, thr_angle_lat);
	link = lateral_cut_position(link, thr_position_lat);
	link = radial_cut_angle(link, thr_angle_rad_intercept, thr_angle_rad_slope);
	link = radial_cut_position(link, thr_position_rad_intercept, thr_position_rad_slope);

	netscan::write_linklet_bin(file_out_link, link);

}
std::vector<netscan::linklet_t> lateral_cut_angle(std::vector<netscan::linklet_t> link, double thr) {
	std::vector<netscan::linklet_t>ret;
	double d_lat;
	for (auto itr = link.begin(); itr != link.end(); itr++) {
		d_lat = ((itr->b[0].ax - itr->b[1].ax) * itr->b[0].ay - (itr->b[0].ay - itr->b[1].ay) * itr->b[0].ax) / sqrt(itr->b[0].ax * itr->b[0].ax + itr->b[0].ay * itr->b[0].ay);
		if (fabs(d_lat) > thr)continue;
		ret.push_back(*itr);
	}
	printf("lateral angle difference <= %5.4lf %d --> %d(%4.1lf%%)\n", thr, link.size(), ret.size(), ret.size() * 100. / link.size());
	return ret;
}
std::vector<netscan::linklet_t> lateral_cut_position(std::vector<netscan::linklet_t> link, double thr) {
	std::vector<netscan::linklet_t> ret;
	double dr, dl;
	for (auto itr = link.begin(); itr != link.end(); itr++) {
		position_difference(*itr, dr, dl);
		if (fabs(dl) > thr)continue;
		ret.push_back(*itr);
	}
	printf("lateral position difference <= %5.1lf %d --> %d(%4.1lf%%)\n", thr, link.size(), ret.size(), ret.size() * 100. / link.size());
	return ret;
}
std::vector<netscan::linklet_t> radial_cut_angle(std::vector<netscan::linklet_t> link, double thr_intercept, double thr_slope) {
	std::vector<netscan::linklet_t>ret;
	double d_rad, angle;
	for (auto itr = link.begin(); itr != link.end(); itr++) {
		angle = sqrt(itr->b[0].ax * itr->b[0].ax + itr->b[0].ay * itr->b[0].ay);
		d_rad = ((itr->b[0].ax - itr->b[1].ax) * itr->b[0].ax + (itr->b[0].ay - itr->b[1].ay) * itr->b[0].ay) / angle;
		if (fabs(d_rad) > thr_intercept + thr_slope * angle)continue;
		ret.push_back(*itr);
	}
	printf("radial angle difference <= %5.4lf+%5.4lf*angle %d --> %d(%4.1lf%%)\n", thr_intercept, thr_slope, link.size(), ret.size(), ret.size() * 100. / link.size());
	return ret;
}
std::vector<netscan::linklet_t> radial_cut_position(std::vector<netscan::linklet_t> link, double thr_intercept, double thr_slope) {
	std::vector<netscan::linklet_t> ret;
	double dr, dl, angle;
	for (auto itr = link.begin(); itr != link.end(); itr++) {
		angle = sqrt(itr->b[0].ax * itr->b[0].ax + itr->b[0].ay * itr->b[0].ay);
		position_difference(*itr, dr, dl);
		if (fabs(dr) > thr_intercept + thr_slope * angle)continue;
		ret.push_back(*itr);
	}
	printf("radial position difference <= %5.4lf+%5.4lf*angle %d --> %d(%4.1lf%%)\n", thr_intercept, thr_slope, link.size(), ret.size(), ret.size() * 100. / link.size());
	return ret;
}
void position_difference(netscan::linklet_t link, double& dr, double& dl) {
	using namespace matrix_3D;
	vector_3D pos0, pos1, dir0, dir1;
	pos0.x = link.b[0].x;
	pos0.y = link.b[0].y;
	pos0.z = link.b[0].z;
	dir0.x = link.b[0].ax;
	dir0.y = link.b[0].ay;
	dir0.z = 1;
	pos1.x = link.b[1].x;
	pos1.y = link.b[1].y;
	pos1.z = link.b[1].z;
	dir1.x = link.b[1].ax;
	dir1.y = link.b[1].ay;
	dir1.z = 1;

	vector_3D base_point, difference;
	//äOë}äÓèÄì_Ç1:1Ç…ì‡ï™ÇµÇΩì_Ç…ê›íË
	base_point = addition(const_multiple(pos0, 0.5), const_multiple(pos1, 0.5));
	difference = addition(const_multiple(pos0, -1), pos1);

	vector_3D extra0, extra1;
	double ratio0, ratio1;
	ratio0 = -1 * dot(addition(pos0, const_multiple(base_point, -1)), difference) / dot(dir0, difference);
	ratio1 = -1 * dot(addition(pos1, const_multiple(base_point, -1)), difference) / dot(dir1, difference);
	extra0 = addition(pos0, const_multiple(dir0, ratio0));
	extra1 = addition(pos1, const_multiple(dir1, ratio1));

	vector_3D unit_r, unit_l;
	unit_l.x = -1 * difference.y;
	unit_l.y = difference.x;
	unit_l.z = 0;
	unit_r.x = -1 * difference.x * difference.z;
	unit_r.y = -1 * difference.y * difference.z;
	unit_r.z = pow(difference.x, 2) + pow(difference.y, 2);

	double constant;
	constant = sqrt(pow(difference.x, 2) + pow(difference.y, 2));
	unit_l.x = unit_l.x / constant;
	unit_l.y = unit_l.y / constant;
	constant = sqrt((pow(difference.x, 2) + pow(difference.y, 2)) * (pow(difference.x, 2) + pow(difference.y, 2) + pow(difference.z, 2)));
	unit_r.x = unit_r.x / constant;
	unit_r.y = unit_r.y / constant;
	unit_r.z = unit_r.z / constant;

	dr = dot(addition(extra1, const_multiple(extra0, -1)), unit_r);
	dl = dot(addition(extra1, const_multiple(extra0, -1)), unit_l);

}

