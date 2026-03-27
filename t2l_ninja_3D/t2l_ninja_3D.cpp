#define _CRT_SECURE_NO_WARNINGS

#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>
#pragma comment(lib, "functions.lib")
#include <functions.hpp>
#pragma comment(lib, "VxxReader.lib")
#include "VxxReader.h"

#include <chrono>
#include <filesystem>
#include <set>
#include <omp.h>
#include <thread>
#include <picojson.h>


class t2l_param {
public:
	std::string file_in_ECC_path, file_in_align_path;
	std::string file_out_linklet_path, file_out_others_path;
	int pl0, pl1;
	double intercept_ax, intercept_ay, intercept_px, intercept_py;
	double slope_px, slope_py, slope_ax, slope_ay;
	double slope2_px, slope2_py, slope2_ax, slope2_ay;
	double intercept_ar, intercept_al, intercept_pr, intercept_pl;
	double slope_pr, slope_pl, slope_ar, slope_al;
	double slope2_pr, slope2_pl, slope2_ar, slope2_al;
	double position_hash, angle_hash;
	double angle_max;
	void Print_all();
};

class align_param {
public:
	int id, signal, ix, iy;
	//視野中心
	double x, y, z;
	//parameter(9);
	double dx, dy, dz, x_rot, y_rot, z_rot, x_shrink, y_shrink, z_shrink, yx_shear, zx_shear, zy_shear;

};
class align_param2 {
public:
	align_param* corr_p[3];
	//3点の視野中心の重心(回転中心)
	double x, y, z;
	//parameter(9);
	double dx, dy, dz, x_rot, y_rot, z_rot, x_shrink, y_shrink, z_shrink, yx_shear, zx_shear, zy_shear;
public:
	//3つのparameterから計算
	void Calc_9param();

};

class output_format_micro {
public:
	int pos, view, imager, zone, isg, ph, vph, px;
};
class output_format_base {
public:
	int pl, rawid;
	double ax, ay, x, y, z;
	output_format_micro m[2];
};

class output_format_link {
public:
	output_format_base b[2];
	double dax, day, dx, dy, dar, dal, dr, dl;
	void Calc_difference();

};
class microtrack_inf {
public:
	int zone, pos, col, row, isg, ph, pixelnum, hitnum;
};

double nominal_gap(std::string file_in_ECC, int pl[2]);
std::string Set_file_read_bvxx_path(std::string file_ECC_path, int pl, int area);
std::string Set_file_read_ali_path(std::string file_align_path, int pl[2], int area);
std::string Set_file_write_link_path(std::string file_link_path, int pl[2]);

std::vector<align_param> read_ali_param(std::string filename, bool output);
std::vector <align_param2 >DelaunayDivide(std::vector <align_param >& corr);
void GaussJorden(double in[3][3], double b[3], double c[3]);
std::vector <std::pair<vxx::base_track_t*, align_param2*>>track_affineparam_correspondence(std::vector<vxx::base_track_t>& base, std::vector <align_param2>& param);
std::vector<vxx::base_track_t*> connect_track(vxx::base_track_t& t, std::vector<vxx::base_track_t*>& connect_cand, t2l_param& param);
bool judge_connect_xy(vxx::base_track_t& t1, vxx::base_track_t& t2, t2l_param& param);
bool judge_connect_rl(vxx::base_track_t& t1, vxx::base_track_t& t2, t2l_param& param);
bool judge_connect_pb(vxx::base_track_t& t1, vxx::base_track_t& t2, t2l_param& param);
void Calc_position_difference(vxx::base_track_t& t1, vxx::base_track_t& t2, double& dr, double& dl);
void output_corrmap2(std::string filename, std::vector<align_param2>& corr);
void output_base_corrmap_pair(std::string filename, std::vector <std::pair<vxx::base_track_t*, align_param2*>>& pair_v);
void output_pair_txt(std::string filename, std::vector<output_format_link>& link);
void output_pair_bin(std::string filename, std::vector<output_format_link>& link);
output_format_link output_format(vxx::base_track_t& t1, vxx::base_track_t& t2);
t2l_param read_param_json(std::string filename);
double select_triangle_vale(align_param2* param, vxx::base_track_t& base);
align_param2* search_param(std::vector<align_param*>& param, vxx::base_track_t& base, std::multimap<int, align_param2*>& triangles);
std::vector<std::pair<vxx::base_track_t*, vxx::base_track_t*>> connect_track_all(std::vector<vxx::base_track_t>& base0, std::vector<vxx::base_track_t>& base1, t2l_param& param);
void trans_base_all(std::vector < std::pair<vxx::base_track_t*, align_param2*>>& track_pair);
void trans_base(std::vector<vxx::base_track_t*>& base, align_param2* param);
std::vector<output_format_link> basetrack_pair_to_linket(std::vector<std::pair<vxx::base_track_t*, vxx::base_track_t*>>& track_pair);
std::vector<microtrack_inf> read_microtrack_inf(std::string filename, bool output);
void Pick_up_pixel_count(std::vector<output_format_link>& link, std::string file_in_ECC);
int64_t file_size(std::string filename);
void zone_trans(std::map<std::tuple<int, int, int, int, int>, int>& pixel_inf);
void add_microtrack_inf_format0(std::string filename, std::multimap<std::tuple<int, int, int, int, int>, int*>& link_px, std::set<int>& zone_pos1, std::set<int>& zone_pos2);
void add_microtrack_inf_format1(std::string filename, std::multimap<std::tuple<int, int, int, int, int>, int*>& link_px);

int main(int argc, char** argv) {
	if (argc != 4) {
		fprintf(stderr, "prg param-json pl0 pl1\n");
		exit(1);
	}
	int area = 0;
	std::string file_in_json = argv[1];

	int pl0 = std::stoi(argv[2]);
	int pl1 = std::stoi(argv[3]);

	int pl[2];
	pl[0] = std::min(pl0, pl1);
	pl[1] = std::max(pl0, pl1);

	t2l_param param = read_param_json(file_in_json);
	std::cout << "ok!" << std::endl;
	param.pl0 = pl[0];
	param.pl1 = pl[1];
	param.Print_all();

	
	//double gap = nominal_gap(param.file_in_ECC_path, pl);
	//読み込みfile名設定
	std::string file_in_base[2];
	for (int i = 0; i < 2; i++) {
		file_in_base[i] = Set_file_read_bvxx_path(param.file_in_ECC_path, pl[i], area);
	}
	std::string file_in_align = Set_file_read_ali_path(param.file_in_align_path, pl, area);

	//fileの存在確認・読み込み
	vxx::BvxxReader br;
	std::vector<vxx::base_track_t> base[2];
	for (int i = 0; i < 2; i++) {
		if (!std::filesystem::exists(file_in_base[i])) {
			fprintf(stderr, "%s not exist\n", file_in_base[i].c_str());
			exit(1);
		}
		base[i] = br.ReadAll(file_in_base[i], pl[i], 0);
	}

	if (!std::filesystem::exists(file_in_align)) {
		fprintf(stderr, "%s not exist\n", file_in_align.c_str());
		exit(1);
	}
	std::vector <corrmap_3d::align_param > corr = corrmap_3d::read_ali_param(file_in_align, false);


	//delaunay3角形分割
	std::vector <corrmap_3d::align_param2 >corr2 = DelaunayDivide(corr);

	//trackとdelaunay3角形の対応
	std::vector < std::pair<vxx::base_track_t*, corrmap_3d::align_param2*>> track_param = track_affineparam_correspondence(base[1], corr2);
	//basetrackを変換
	corrmap_3d::trans_base_all(track_param);

	//飛跡接続
	//base0をhash(位置角度)
	std::vector<std::pair<vxx::base_track_t*, vxx::base_track_t*>> track_pair = connect_track_all(base[0], base[1], param);

	//format変換
	std::vector<output_format_link> link = basetrack_pair_to_linket(track_pair);

	//pixelcountを拾う
	Pick_up_pixel_count(link, param.file_in_ECC_path);

	//linklet出力
	output_pair_bin(Set_file_write_link_path(param.file_out_linklet_path, pl), link);
	//output_pair_txt(Set_file_write_link_path(param.file_out_linklet_path, pl), link);


	//std::string file_output_corr, file_output_base_corr;
	////corrmapの出力
	//output_corrmap2(file_output_corr, corr2);
	////basetrack-corrmapidの出力
	//output_base_corrmap_pair(file_output_base_corr, track_param);

}

//json fileでのparameter読み込み
t2l_param read_param_json(std::string filename) {

	std::ifstream ifs(filename, std::ios::in);
	if (ifs.fail()) {
		std::cerr << "failed to read test.json" << std::endl;
		std::cerr << "filename = " << filename << std::endl;
		exit(1);
	}
	const std::string json((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
	ifs.close();

	// JSONデータを解析する。
	picojson::value v;
	const std::string err = picojson::parse(v, json);
	if (err.empty() == false) {
		std::cerr << err << std::endl;
		exit(1);
	}
	t2l_param ret;
	picojson::object& all = v.get<picojson::object>();
	ret.file_in_ECC_path = all["file_in_ECC_path"].get<std::string>();
	ret.file_in_align_path = all["file_in_align_path"].get<std::string>();
	ret.file_out_linklet_path = all["file_out_linklet_path"].get<std::string>();
	ret.file_out_others_path = all["file_out_others"].get<std::string>();
	picojson::object& connect_param_angle = all["connect_param_angle"].get<picojson::object>();

	picojson::object& connect_param_position = all["connect_param_position"].get<picojson::object>();

	ret.intercept_ax = connect_param_angle["intercept_x"].get<double>();
	ret.intercept_ay = connect_param_angle["intercept_y"].get<double>();
	ret.intercept_ar = connect_param_angle["intercept_r"].get<double>();
	ret.intercept_al = connect_param_angle["intercept_l"].get<double>();
	ret.slope_ax = connect_param_angle["slope_x"].get<double>();
	ret.slope_ay = connect_param_angle["slope_y"].get<double>();
	ret.slope_ar = connect_param_angle["slope_r"].get<double>();
	ret.slope_al = connect_param_angle["slope_l"].get<double>();
	ret.slope2_ax = connect_param_angle["slope2_x"].get<double>();
	ret.slope2_ay = connect_param_angle["slope2_y"].get<double>();
	ret.slope2_ar = connect_param_angle["slope2_r"].get<double>();
	ret.slope2_al = connect_param_angle["slope2_l"].get<double>();


	ret.intercept_px = connect_param_position["intercept_x"].get<double>();
	ret.intercept_py = connect_param_position["intercept_y"].get<double>();
	ret.intercept_pr = connect_param_position["intercept_r"].get<double>();
	ret.intercept_pl = connect_param_position["intercept_l"].get<double>();
	ret.slope_px = connect_param_position["slope_x"].get<double>();
	ret.slope_py = connect_param_position["slope_y"].get<double>();
	ret.slope_pr = connect_param_position["slope_r"].get<double>();
	ret.slope_pl = connect_param_position["slope_l"].get<double>();
	ret.slope2_px = connect_param_position["slope2_x"].get<double>();
	ret.slope2_py = connect_param_position["slope2_y"].get<double>();
	ret.slope2_pr = connect_param_position["slope2_r"].get<double>();
	ret.slope2_pl = connect_param_position["slope2_l"].get<double>();

	ret.position_hash = all["position_hash"].get<double>();
	ret.angle_hash = all["angle_hash"].get<double>();
	ret.angle_max = all["angle_max"].get<double>();
	std::cout << "ok?" << std::endl;

	return ret;
}
void t2l_param::Print_all() {

	printf("file_in_ECC_path : %s\n", file_in_ECC_path.c_str());
	printf("file_in_align_path : %s\n", file_in_align_path.c_str());
	printf("file_out_linklet_path : %s\n", file_out_linklet_path.c_str());
	printf("file_out_others_path : %s\n", file_out_others_path.c_str());


	printf("PL0 : %03d \n", pl0);
	printf("PL1 : %03d \n", pl1);

	printf("position hash : %g[um]\n", position_hash);
	printf("angle hash : %\g\n", angle_hash);

	printf("dax: %.5lf * tan^2\u03B8x + %.5lf * tan\u03B8x + %.5lf\n", slope2_ax, slope_ax, intercept_ax);
	printf("day: %.5lf * tan^2\u03B8y + %.5lf * tan\u03B8y + %.5lf\n", slope2_ay, slope_ay, intercept_ay);
	printf("dar: %.5lf * tan^2\u03B8  + %.5lf * tan\u03B8  + %.5lf\n", slope2_ar, slope_ar, intercept_ar);
	printf("dal: %.5lf * tan^2\u03B8  + %.5lf * tan\u03B8  + %.5lf\n", slope2_al, slope_al, intercept_al);

	printf("dpx: %.1lf * tan^2\u03B8x + %.1lf * tan\u03B8x + %.1lf\n", slope2_px, slope_px, intercept_px);
	printf("dpy: %.1lf * tan^2\u03B8y + %.1lf * tan\u03B8y + %.1lf\n", slope2_py, slope_py, intercept_py);
	printf("dpr: %.1lf * tan^2\u03B8  + %.1lf * tan\u03B8  + %.1lf\n", slope2_pr, slope_pr, intercept_pr);
	printf("dpl: %.1lf * tan^2\u03B8  + %.1lf * tan\u03B8  + %.1lf\n", slope2_pl, slope_pl, intercept_pl);

}

std::string Set_file_read_bvxx_path(std::string file_ECC_path, int pl, int area) {
	std::stringstream file_in_base;
	file_in_base << file_ECC_path << "\\Area" << area << "\\PL" << std::setw(3) << std::setfill('0') << pl
		<< "\\b" << std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";
	return file_in_base.str();

}
std::string Set_file_read_ali_path(std::string file_align_path, int pl[2], int area) {
	std::stringstream file_in_ali;
	file_in_ali << file_align_path << "\\ali_"
		<< std::setw(3) << std::setfill('0') << pl[0] << "_"
		<< std::setw(3) << std::setfill('0') << pl[1] << "_interpolation.txt";
	return file_in_ali.str();
}
std::string Set_file_write_link_path(std::string file_link_path, int pl[2]) {
	std::stringstream file_out_link;
	file_out_link << file_link_path << "\\link_"
		<< std::setw(3) << std::setfill('0') << pl[0] << "_"
		<< std::setw(3) << std::setfill('0') << pl[1] << ".bin";
	return file_out_link.str();
}

/*
std::vector<align_param> read_ali_param(std::string filename, bool output) {

	std::vector<align_param> ret;
	align_param param_tmp;
	std::ifstream ifs(filename);

	while (ifs >> param_tmp.id >> param_tmp.ix >> param_tmp.iy >> param_tmp.signal
		>> param_tmp.x >> param_tmp.y >> param_tmp.z
		>> param_tmp.x_rot >> param_tmp.y_rot >> param_tmp.z_rot
		>> param_tmp.x_shrink >> param_tmp.y_shrink >> param_tmp.z_shrink
		>> param_tmp.yx_shear >> param_tmp.zx_shear >> param_tmp.zy_shear
		>> param_tmp.dx >> param_tmp.dy >> param_tmp.dz) {
		ret.push_back(param_tmp);
		//printf("ix %d iy%d\n", param_tmp.ix, param_tmp.iy);

	}
	if (output == 1) {
		fprintf(stderr, "%s input finish\n", filename.c_str());
	}
	if (ret.size() == 0) {
		fprintf(stderr, "%s alignment miss!\n", filename.c_str());
		exit(1);
	}
	return ret;

}

std::vector <align_param2 >DelaunayDivide(std::vector <align_param >&corr) {

	//delaunay分割
	std::vector<double> x, y;
	for (auto itr = corr.begin(); itr != corr.end(); itr++) {
		x.push_back(itr->x);
		y.push_back(itr->y);
	}

	delaunay::DelaunayTriangulation DT(x, y); // (std::vector<double> x, std::vector<double> y, uint32_t seed_)
	DT.execute(); // (double min_delta = 1e-6, double max_delta = 1e-5, int max_miss_count = 30)
	std::vector<delaunay::Edge> edge = DT.get_edges();

	std::multimap<int, int> edge_map;

	for (auto itr = edge.begin(); itr != edge.end(); itr++) {
		edge_map.insert(std::make_pair(std::min(itr->first, itr->second), std::max(itr->first, itr->second)));

	}
	std::set<std::tuple<int, int, int>>triangle;
	std::set<int> vertex;
	for (auto itr = edge_map.begin(); itr != edge_map.end(); itr++) {
		//itr->firstの点=aを通る三角形の探索
		vertex.clear();
		auto range = edge_map.equal_range(itr->first);
		//aを通りitr->secondの点=bに行く。bのsetを作成
		for (auto res = range.first; res != range.second; res++) {
			vertex.insert(res->second);
		}
		//bを通る線分の探索
		for (auto itr2 = vertex.begin(); itr2 != vertex.end(); itr2++) {
			if (edge_map.count(*itr2) == 0)continue;
			auto range2 = edge_map.equal_range(*itr2);
			//bを通る線分の中からaから始まる線分を探す
			for (auto res = range2.first; res != range2.second; res++) {
				if (vertex.count(res->second) == 1) {
					triangle.insert(std::make_tuple(itr->first, *itr2, res->second));
				}
			}

		}
	}

	std::vector <align_param2 > ret;
	for (auto itr = triangle.begin(); itr != triangle.end(); itr++) {
		//printf("delaunay triangle %d %d %d\n", std::get<0>(*itr), std::get<1>(*itr), std::get<2>(*itr));
		align_param2 param;
		param.corr_p[0] = &(corr[std::get<0>(*itr)]);
		param.corr_p[1] = &(corr[std::get<1>(*itr)]);
		param.corr_p[2] = &(corr[std::get<2>(*itr)]);
		param.x = 0;
		param.y = 0;
		param.z = 0;
		param.z_shrink = 0;
		for (int i = 0; i < 3; i++) {
			param.x += param.corr_p[i]->x;
			param.y += param.corr_p[i]->y;
			param.z += param.corr_p[i]->z;
			param.z_shrink += param.corr_p[i]->z_shrink;
		}
		param.x = param.x / 3;
		param.y = param.y / 3;
		param.z = param.z / 3;
		param.z_shrink = param.z_shrink / 3;
		param.zx_shear = 0;
		param.zy_shear = 0;

		param.Calc_9param();

		ret.push_back(param);
	}

	return ret;

}
void align_param2::Calc_9param() {


	double bp[3][3], ap[3][3], cos_z, sin_z;
	for (int i = 0; i < 3; i++) {
		bp[i][0] = corr_p[i]->x;
		bp[i][1] = corr_p[i]->y;
		bp[i][2] = corr_p[i]->z;

		cos_z = cos(corr_p[i]->z_rot);
		sin_z = sin(corr_p[i]->z_rot);


		ap[i][0] = corr_p[i]->x_shrink*cos_z*(corr_p[i]->x) - corr_p[i]->y_shrink*sin_z*(corr_p[i]->y) + corr_p[i]->dx;
		ap[i][1] = corr_p[i]->x_shrink*sin_z*(corr_p[i]->x) + corr_p[i]->y_shrink*cos_z*(corr_p[i]->y) + corr_p[i]->dy;
		ap[i][2] = corr_p[i]->z + corr_p[i]->dz;
		//printf("bp%d %8.1lf %8.1lf %8.1lf\n",i, bp[i][0], bp[i][1], bp[i][2]);
		//printf("ap%d %8.1lf %8.1lf %8.1lf\n", i,ap[i][0], ap[i][1], ap[i][2]);
	}
	//apの位置ずれvectorを定義
	double dp[2][3];
	for (int i = 0; i < 3; i++) {
		dp[0][i] = ap[1][i] - ap[0][i];
		dp[1][i] = ap[2][i] - ap[0][i];
	}
	//printf("0-->1 x,y,z : %.1lf %.1lf %.1lf\n", dp[0][0], dp[0][1], dp[0][2]);
	//printf("0-->2 x,y,z : %.1lf %.1lf %.1lf\n", dp[1][0], dp[1][1], dp[1][2]);
	//法線vector
	double n_v[3];
	n_v[0] = (dp[0][1] * dp[1][2] - dp[0][2] * dp[1][1]);
	n_v[1] = (dp[0][2] * dp[1][0] - dp[0][0] * dp[1][2]);
	n_v[2] = (dp[0][0] * dp[1][1] - dp[0][1] * dp[1][0]);

	//std::cout << "normal vector" << std::endl;
	//for (int i = 0; i < 3; i++) {
	//	std::cout << std::setw(14) << std::fixed << std::setprecision(10) << n_v[i] << std::endl;
	//}

	x_rot = atan(n_v[1] / n_v[2]);
	n_v[1] = cos(x_rot)*n_v[1] - sin(x_rot)*n_v[2];
	n_v[2] = sin(x_rot)*n_v[1] + cos(x_rot)*n_v[2];
	//std::cout << "normal vector" << std::endl;
	//for (int i = 0; i < 3; i++) {
	//	std::cout << std::setw(14) << std::fixed << std::setprecision(10) << n_v[i] << std::endl;
	//}
	y_rot = atan(-1 * n_v[0] / n_v[2]);
	n_v[0] = cos(y_rot)*n_v[0] + sin(y_rot)*n_v[2];
	n_v[2] = -1 * sin(y_rot)*n_v[0] + cos(y_rot)*n_v[2];

	//std::cout << "normal vector" << std::endl;
	//for (int i = 0; i < 3; i++) {
	//	std::cout << std::setw(14) << std::fixed << std::setprecision(10) << n_v[i] << std::endl;
	//}

	//printf("x rot:%.6lf\n", x_rot);
	//printf("y rot:%.6lf\n", y_rot);


	matrix_3D::matrix_33 x_rot_mat(0, x_rot), y_rot_mat(1, y_rot);
	matrix_3D::vector_3D ap_v[3];
	for (int i = 0; i < 3; i++) {
		ap_v[i].x = ap[i][0];
		ap_v[i].y = ap[i][1];
		ap_v[i].z = ap[i][2];
	}
	for (int i = 0; i < 3; i++) {
		ap_v[i].matrix_multiplication(x_rot_mat);
		ap_v[i].matrix_multiplication(y_rot_mat);
	}
	//for (int i = 0; i < 3; i++) {
	//	printf("point %d\n", i);
	//	printf("\t %.2lf %.2lf %.2lf\n", bp[i][0], bp[i][1], bp[i][2]);
	//	printf("\t %.2lf %.2lf %.2lf\n", ap_v[i].x, ap_v[i].y, ap_v[i].z);
	//}
	dz = (ap_v[0].z - bp[0][2] + ap_v[1].z - bp[1][2] + ap_v[2].z - bp[2][2]) / 3;
	//printf("dz=%.2lf\n", dz);
	//3元方程式を解く
	double a[2][3][3] = { { {bp[0][0],bp[0][1],1},{bp[1][0],bp[1][1],1},{bp[2][0],bp[2][1],1} },  { {bp[0][0],bp[0][1],1},{bp[1][0],bp[1][1],1},{bp[2][0],bp[2][1],1} } };
	double b[2][3] = { { ap_v[0].x,ap_v[1].x,ap_v[2].x },{ ap_v[0].y,ap_v[1].y,ap_v[2].y } };
	double c[2][3] = { {1, 1, 1},{1,1,1} };
	//gauss(a[0], b[0], c[0]);
	//gauss(a[1], b[1], c[1]);
	GaussJorden(a[0], b[0], c[0]);
	GaussJorden(a[1], b[1], c[1]);
	z_rot = atan(c[1][0] / c[0][0]);
	x_shrink = c[0][0] / cos(z_rot);
	y_shrink = (c[0][0] * c[1][1] - c[0][1] * c[1][0]) / (c[0][0] * cos(z_rot) + c[1][0] * sin(z_rot));
	yx_shear = (c[0][1] * cos(z_rot) + c[1][1] * sin(z_rot)) / (c[0][0] * cos(z_rot) + c[1][0] * sin(z_rot));

	dx = c[0][2];
	dy = c[1][2];
	matrix_3D::vector_3D dr;
	dr.x = c[0][2];
	dr.y = c[1][2];
	dr.z = dz;

	x_rot = x_rot * -1;
	y_rot = y_rot * -1;
	matrix_3D::matrix_33 x_rot_mat_inv(0, x_rot), y_rot_mat_inv(1, y_rot);


	dr.matrix_multiplication(y_rot_mat_inv);
	dr.matrix_multiplication(x_rot_mat_inv);

	dx = dr.x;
	dy = dr.y;
	dz = dr.z;

	//printf("x rot: %.6lf\n",x_rot);
	//printf("y rot: %.6lf\n",y_rot);
	//printf("z rot: %.6lf\n",z_rot);
	//printf("x shrink: %.6lf\n", x_shrink);
	//printf("y shrink: %.6lf\n", y_shrink);
	//printf("z shrink: %.6lf\n", z_shrink);
	//printf("x shift: %.5lf\n", dx);
	//printf("y shift: %.5lf\n", dy);
	//printf("z shift: %.5lf\n", dz);
	//printf("yx shear: %.6lf\n", yx_shear);
	//printf("zx shear: %.6lf\n", zx_shear);
	//printf("zy shear: %.6lf\n", zy_shear);

	//std::vector< matrix_3D::vector_3D >point,point_after;
	//for (int i = 0; i < 3; i++) {
	//	matrix_3D::vector_3D p;
	//	p.x = corr_p[i]->x;
	//	p.y = corr_p[i]->y;
	//	p.z = corr_p[i]->z;
	//	point.push_back(p);
	//	p.x = corr_p[i]->x + corr_p[i]->dx;
	//	p.y = corr_p[i]->y + corr_p[i]->dy;
	//	p.z = corr_p[i]->z + corr_p[i]->dz;
	//	point_after.push_back(p);
	//}
	//trans_9para(point, *this);
	//for (auto p : point_after) {
	//	printf("x:%10.1lf y:%10.1lf z:%10.1lf\n", p.x, p.y, p.z);
	//}
}

void GaussJorden(double in[3][3], double b[3], double c[3]) {


	double a[3][4];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 4; j++) {
			if (j < 3) {
				a[i][j] = in[i][j];
			}
			else {
				a[i][j] = b[i];
			}
		}
	}
	int N = 3;
	double p, d;         // ピボット係数、ピボット行ｘ係数
	double max, dummy;   // 最大絶対値、入れ替え時ダミー
	int s;

	//元の連立方程式をコンソール出力
   //for (int i = 0; i < N; i++) {
   //	for (int j = 0; j < N; j++)
   //		printf("%+fx%d ", a[i][j], j + 1);
   //	printf("= %+f\n", a[i][N]);
   //}

	for (int k = 0; k < N; k++) {
		// 行入れ替え
		max = 0; s = k;
		for (int j = k; j < N; j++) {
			if (fabs(a[j][k]) > max) {
				max = fabs(a[j][k]);
				s = j;
			}
		}
		if (max == 0) {
			printf("解けない！");
			exit(1);
		}
		for (int j = 0; j <= N; j++) {
			dummy = a[k][j];
			a[k][j] = a[s][j];
			a[s][j] = dummy;
		}

		// ピボット係数
		p = a[k][k];

		// ピボット行を p で除算
		for (int j = k; j < N + 1; j++)
			a[k][j] /= p;

		// ピボット列の掃き出し
		for (int i = 0; i < N; i++) {
			if (i != k) {
				d = a[i][k];
				for (int j = k; j < N + 1; j++)
					a[i][j] -= d * a[k][j];
			}
		}
	}

	// 結果出力
	for (int k = 0; k < N; k++) {
		c[k] = a[k][N];
		//printf("x%d = %f\n", k + 1, a[k][N]);
	}
}

//basetrack-alignment mapの対応
std::vector <std::pair<vxx::base_track_t*, align_param2*>>track_affineparam_correspondence(std::vector<vxx::base_track_t>&base, std::vector <align_param2> &param) {

	//local alignの視野中心を取り出して、位置でhash
	//local alignの視野中心の作るdelaunay三角形をmapで対応

	std::map<int, align_param*> view_center;
	std::multimap<int, align_param2*>triangles;
	double xmin = 999999, ymin = 999999, hash = 2000;
	for (auto itr = param.begin(); itr != param.end(); itr++) {
		for (int i = 0; i < 3; i++) {
			view_center.insert(std::make_pair(itr->corr_p[i]->id, (itr->corr_p[i])));
			triangles.insert(std::make_pair(itr->corr_p[i]->id, &(*itr)));
			xmin = std::min(itr->corr_p[i]->x, xmin);
			ymin = std::min(itr->corr_p[i]->y, ymin);
		}
	}
	std::multimap<std::pair<int, int>, align_param*> view_center_hash;
	std::pair<int, int>id;
	for (auto itr = view_center.begin(); itr != view_center.end(); itr++) {
		id.first = int((itr->second->x - xmin) / hash);
		id.second = int((itr->second->y - ymin) / hash);
		view_center_hash.insert(std::make_pair(id, itr->second));
	}

	std::vector < std::pair<vxx::base_track_t*, align_param2*>> ret;
	std::vector<align_param*> param_cand;
	int loop = 0, ix, iy, count = 0;
	for (auto itr = base.begin(); itr != base.end(); itr++) {
		if (count % 100000 == 0) {
			printf("\r search correspond triangles %d/%d(%4.1lf%%)", count, base.size(), count*100. / base.size());
		}
		count++;
		ix = (itr->x - xmin) / hash;
		iy = (itr->y - ymin) / hash;
		loop = 1;
		while (true) {
			param_cand.clear();
			for (int iix = ix - loop; iix <= ix + loop; iix++) {
				for (int iiy = iy - loop; iiy <= iy + loop; iiy++) {
					id.first = iix;
					id.second = iiy;
					if (view_center_hash.count(id) != 0) {
						auto range = view_center_hash.equal_range(id);
						for (auto res = range.first; res != range.second; res++) {
							param_cand.push_back(res->second);
						}
					}
				}
			}
			if (param_cand.size() > 2)break;
			loop++;
		}
		align_param2* param2 = search_param(param_cand, *itr, triangles);
		ret.push_back(std::make_pair(&(*itr), param2));
	}
	printf("\r search correspond triangles %d/%d(%4.1lf%%)\n", count, base.size(), count*100. / base.size());

	return ret;
}
align_param2* search_param(std::vector<align_param*> &param, vxx::base_track_t&base, std::multimap<int, align_param2*>&triangles) {
	//三角形内部
	//最近接三角形
	double dist = 0;
	std::map<double, align_param* > dist_map;
	//align_paramを近い順にsort
	for (auto itr = param.begin(); itr != param.end(); itr++) {
		dist = ((*itr)->x - base.x)*((*itr)->x - base.x) + ((*itr)->y - base.y)*((*itr)->y - base.y);
		dist_map.insert(std::make_pair(dist, (*itr)));
	}

	double sign[3];
	bool flg = false;
	int id;

	align_param2* ret = triangles.begin()->second;
	for (auto itr = dist_map.begin(); itr != dist_map.end(); itr++) {
		if (itr != dist_map.begin())continue;


		//corrmapのID
		id = itr->second->id;
		if (triangles.count(id) == 0) {
			fprintf(stderr, "alignment triangle ID=%d not found\n", id);
			exit(1);
		}
		//idの属する三角形を探索
		auto range = triangles.equal_range(id);
		for (auto itr2 = range.first; itr2 != range.second; itr2++) {
			sign[0] = (itr2->second->corr_p[1]->x - itr2->second->corr_p[0]->x)*(base.y - itr2->second->corr_p[1]->y) - (itr2->second->corr_p[1]->y - itr2->second->corr_p[0]->y)*(base.x - itr2->second->corr_p[1]->x);
			sign[1] = (itr2->second->corr_p[2]->x - itr2->second->corr_p[1]->x)*(base.y - itr2->second->corr_p[2]->y) - (itr2->second->corr_p[2]->y - itr2->second->corr_p[1]->y)*(base.x - itr2->second->corr_p[2]->x);
			sign[2] = (itr2->second->corr_p[0]->x - itr2->second->corr_p[2]->x)*(base.y - itr2->second->corr_p[0]->y) - (itr2->second->corr_p[0]->y - itr2->second->corr_p[2]->y)*(base.x - itr2->second->corr_p[0]->x);
			//printf("point %.lf,%.1lf\n", base.x, base.y);
			//printf("triangle %.1lf %.1lf\n", itr2->second->corr_p[0]->x, itr2->second->corr_p[0]->y);
			//printf("triangle %.1lf %.1lf\n", itr2->second->corr_p[1]->x, itr2->second->corr_p[1]->y);
			//printf("triangle %.1lf %.1lf\n", itr2->second->corr_p[2]->x, itr2->second->corr_p[2]->y);
			//printf("sign %.1lf %1.lf %.1lf\n", sign[0], sign[1], sign[2]);
			//printf("  signbit %d %d %d\n", std::signbit(sign[0]), std::signbit(sign[1]), std::signbit(sign[2]));
			//printf("n signbit %d %d %d\n", !std::signbit(sign[0]), !std::signbit(sign[1]), !std::signbit(sign[2]));
			//printf("judge %d\n", (std::signbit(sign[0]) && std::signbit(sign[1]) && std::signbit(sign[2])) || (!std::signbit(sign[0]) && !std::signbit(sign[1]) && !std::signbit(sign[2])));
			//printf("\n");

			//符号が3つとも一致でtrue
			if ((std::signbit(sign[0]) && std::signbit(sign[1]) && std::signbit(sign[2])) || (!std::signbit(sign[0]) && !std::signbit(sign[1]) && !std::signbit(sign[2]))) {
				ret = itr2->second;
				flg = true;
				break;
			}
		}
		if (flg)break;
	}
	if (flg) {
		//printf("point in trianlge\n");
		return ret;
	}

	//distが最小になるcorrmapをとってくる
	dist = -1;
	for (auto itr = dist_map.begin(); itr != dist_map.end(); itr++) {
		//corrmapのID
		id = itr->second->id;
		if (triangles.count(id) == 0) {
			fprintf(stderr, "alignment triangle ID=%d not found\n", id);
			exit(1);
		}
		//idの属する三角形を探索
		auto range = triangles.equal_range(id);
		for (auto itr2 = range.first; itr2 != range.second; itr2++) {
			if (dist<0 || dist>select_triangle_vale(itr2->second, base)) {
				dist = select_triangle_vale(itr2->second, base);
				ret = itr2->second;
			}
		}
	}
	//printf("point not in trianlge\n");
	return ret;
}
double select_triangle_vale(align_param2* param, vxx::base_track_t&base) {
	double x, y;
	double dist = 0;
	x = (param->corr_p[0]->x + param->corr_p[1]->x + param->corr_p[2]->x) / 3;
	y = (param->corr_p[0]->y + param->corr_p[1]->y + param->corr_p[2]->y) / 3;
	dist = (base.x - x)*(base.x - x) + (base.y - y)*(base.y - y);
	return dist;
}

//変換 zshrink補正-->9para変換
void trans_base_all(std::vector < std::pair<vxx::base_track_t*, align_param2*>>&track_pair) {
	std::map<std::tuple<int, int, int>, align_param2*> param_map;
	std::multimap<std::tuple<int, int, int>, vxx::base_track_t*>base_map;
	std::tuple<int, int, int>id;
	//三角形ごとにbasetrackをまとめる
	for (auto itr = track_pair.begin(); itr != track_pair.end(); itr++) {
		std::get<0>(id) = itr->second->corr_p[0]->id;
		std::get<1>(id) = itr->second->corr_p[1]->id;
		std::get<2>(id) = itr->second->corr_p[2]->id;
		param_map.insert(std::make_pair(id, itr->second));
		base_map.insert(std::make_pair(id, itr->first));
	}


	//ここで三角形ごとに変換
	int count = 0;
	std::vector<vxx::base_track_t*> t_base;
	for (auto itr = param_map.begin(); itr != param_map.end(); itr++) {
		if (count % 1000 == 0) {
			printf("\r basetrack trans num of triangles %d/%d(%4.1lf%%)", count, param_map.size(), count*100. / param_map.size());
		}
		count++;

		t_base.clear();

		if (base_map.count(itr->first) == 0)continue;
		auto range = base_map.equal_range(itr->first);
		for (auto res = range.first; res != range.second; res++) {
			t_base.push_back(res->second);
		}
		trans_base(t_base, itr->second);

	}
	printf("\r basetrack trans num of triangles %d/%d(%4.1lf%%)\n", count, param_map.size(), count*100. / param_map.size());

}
void trans_base(std::vector<vxx::base_track_t*>&base, align_param2 *param) {

	matrix_3D::matrix_33 x_rot_mat(0, param->x_rot), y_rot_mat(1, param->y_rot), z_rot_mat(2, param->z_rot), all_trans(0, 0), shear_mat(0, 0), shrink_mat(0, 0);

	shrink_mat.val[0][0] *= param->x_shrink;
	shrink_mat.val[1][1] *= param->y_shrink;
	//shrink_mat.val[2][2] *= param->z_shrink;
	shear_mat.val[0][1] = param->yx_shear;
	shear_mat.val[0][2] = param->zx_shear;
	shear_mat.val[1][2] = param->zy_shear;

	matrix_3D::vector_3D shift, center;
	center.x = param->x;
	center.y = param->y;
	center.z = param->z;
	shift.x = param->dx;
	shift.y = param->dy;
	shift.z = param->dz;

	all_trans.matrix_multiplication(shear_mat);
	all_trans.matrix_multiplication(shrink_mat);
	all_trans.matrix_multiplication(z_rot_mat);
	all_trans.matrix_multiplication(y_rot_mat);
	all_trans.matrix_multiplication(x_rot_mat);

	//all_trans.Print();
	matrix_3D::vector_3D base_p0, base_p1;
	for (auto itr = base.begin(); itr != base.end(); itr++) {
		base_p0.x = (*itr)->x;
		base_p0.y = (*itr)->y;
		base_p0.z = param->z;

		base_p1.x = (*itr)->x + (*itr)->ax*((*itr)->m[1].z - (*itr)->m[0].z);
		base_p1.y = (*itr)->y + (*itr)->ay*((*itr)->m[1].z - (*itr)->m[0].z);
		//角度shrink分はここでかける
		base_p1.z = param->z + ((*itr)->m[1].z - (*itr)->m[0].z) / param->z_shrink;

		//視野中心を原点に移動
		//base_p0 = matrix_3D::addition(base_p0, matrix_3D::const_multiple(center, -1));
		//base_p1 = matrix_3D::addition(base_p1, matrix_3D::const_multiple(center, -1));

		//変換の実行
		base_p0.matrix_multiplication(all_trans);
		base_p0 = matrix_3D::addition(base_p0, shift);
		base_p1.matrix_multiplication(all_trans);
		base_p1 = matrix_3D::addition(base_p1, shift);

		//原点をもとに戻す
		//base_p0 = matrix_3D::addition(base_p0, center);
		//base_p1 = matrix_3D::addition(base_p1, center);

		(*itr)->x = base_p0.x;
		(*itr)->y = base_p0.y;
		(*itr)->z = base_p0.z;

		//printf("ax:%.4lf --> %.4lf\n", (*itr)->ax, (base_p0.x - base_p1.x) / (base_p0.z - base_p1.z));
		//printf("ay:%.4lf --> %.4lf\n", (*itr)->ay, (base_p0.y - base_p1.y) / (base_p0.z - base_p1.z));

		(*itr)->ax = (base_p0.x - base_p1.x) / (base_p0.z - base_p1.z);
		(*itr)->ay = (base_p0.y - base_p1.y) / (base_p0.z - base_p1.z);

	}
}
*/

//接続関連
std::vector<std::pair<vxx::base_track_t*, vxx::base_track_t*>> connect_track_all(std::vector<vxx::base_track_t>& base0, std::vector<vxx::base_track_t>& base1, t2l_param& param) {

	std::vector<std::pair<vxx::base_track_t*, vxx::base_track_t*>> connected;
	const double hash_pos_size = param.position_hash;
	const double hash_ang_size = param.angle_hash;
	double min[4];
	for (auto itr = base0.begin(); itr != base0.end(); itr++) {
		if (itr == base0.begin()) {
			min[0] = itr->x;
			min[1] = itr->y;
			min[2] = itr->ax;
			min[3] = itr->ay;
		}
		min[0] = std::min(min[0], itr->x);
		min[1] = std::min(min[1], itr->y);
		min[2] = std::min(min[2], itr->ax);
		min[3] = std::min(min[3], itr->ay);
		itr->z = 0;
	}

	std::multimap< std::tuple<int, int, int, int>, vxx::base_track_t*>base_map;
	std::tuple<int, int, int, int> id;
	for (auto itr = base0.begin(); itr != base0.end(); itr++) {
		std::get<0>(id) = (itr->x - min[0]) / hash_pos_size;
		std::get<1>(id) = (itr->y - min[1]) / hash_pos_size;
		std::get<2>(id) = (itr->ax - min[2]) / hash_ang_size;
		std::get<3>(id) = (itr->ay - min[3]) / hash_ang_size;
		base_map.insert(std::make_pair(id, &(*itr)));
	}

	//base1
	double ex_x, ex_y;
	double x_min, x_max, y_min, y_max;
	double ax_min, ax_max, ay_min, ay_max;
	double angle_acc_x, angle_acc_y;
	double position_acc_x, position_acc_y;
	int i_ax_min, i_ax_max, i_ay_min, i_ay_max;
	int i_x_min, i_x_max, i_y_min, i_y_max;
	std::vector<vxx::base_track_t*> connect_cand;
	int count = 0;
	for (auto itr = base1.begin(); itr != base1.end(); itr++) {
		if (count % 10000 == 0) {
			printf("\r connect track %10d/%10d (%4.1lf%%)", count, base1.size(), count * 100. / base1.size());
		}
		count++;
		if (fabs(itr->ax) > param.angle_max)continue;
		if (fabs(itr->ay) > param.angle_max)continue;
		connect_cand.clear();
		//printf("%g %g %g %g\n", itr->ax, itr->ay, itr->x, itr->y);
		//allowance 
		angle_acc_x = param.intercept_ax + param.slope_ax * fabs(itr->ax) + param.slope2_ax * pow(itr->ax, 2);
		angle_acc_y = param.intercept_ay + param.slope_ay * fabs(itr->ay) + param.slope2_ay * pow(itr->ay, 2);
		position_acc_x = param.intercept_px + param.slope_px * fabs(itr->ax) + param.slope2_px * pow(itr->ax, 2);
		position_acc_y = param.intercept_py + param.slope_py * fabs(itr->ay) + param.slope2_py * pow(itr->ay, 2);
		//z=z1 --> (z0+z1)/2に外挿
		ex_x = itr->ax * (-1) * itr->z / 2 + itr->x;
		ex_y = itr->ay * (-1) * itr->z / 2 + itr->y;
		//中点に外挿-->位置ずれallowanceずらす-->角度ずれallowanceずらしてもう半分外挿
		x_min = ex_x - position_acc_x * 1.5 - (itr->ax - angle_acc_x * 1.5) * itr->z / 2;
		x_max = ex_x + position_acc_x * 1.5 - (itr->ax + angle_acc_x * 1.5) * itr->z / 2;
		y_min = ex_y - position_acc_y * 1.5 - (itr->ay - angle_acc_y * 1.5) * itr->z / 2;
		y_max = ex_y + position_acc_y * 1.5 - (itr->ay + angle_acc_y * 1.5) * itr->z / 2;
		//ang_center+-(angle_acc)*sigma
		ax_min = itr->ax - angle_acc_x * 1.5;
		ax_max = itr->ax + angle_acc_x * 1.5;
		ay_min = itr->ay - angle_acc_y * 1.5;
		ay_max = itr->ay + angle_acc_y * 1.5;

		//hash idに変換
		i_x_min = (x_min - min[0]) / hash_pos_size;
		i_x_max = (x_max - min[0]) / hash_pos_size;
		i_y_min = (y_min - min[1]) / hash_pos_size;
		i_y_max = (y_max - min[1]) / hash_pos_size;
		i_ax_min = (ax_min - min[2]) / hash_ang_size;
		i_ax_max = (ax_max - min[2]) / hash_ang_size;
		i_ay_min = (ay_min - min[3]) / hash_ang_size;
		i_ay_max = (ay_max - min[3]) / hash_ang_size;
		//printf("%d %d %d %d\n", i_x_min, i_x_max, i_y_min, i_y_max);
		//printf("%d %d %d %d\n", i_ax_min, i_ax_max, i_ay_min, i_ay_max);
		//hash mapの中から該当trackを探す
		for (int ix = i_x_min; ix <= i_x_max; ix++) {
			for (int iy = i_y_min; iy <= i_y_max; iy++) {
				for (int iax = i_ax_min; iax <= i_ax_max; iax++) {
					for (int iay = i_ay_min; iay <= i_ay_max; iay++) {
						std::get<0>(id) = ix;
						std::get<1>(id) = iy;
						std::get<2>(id) = iax;
						std::get<3>(id) = iay;
						if (base_map.count(id) == 0)continue;
						auto range = base_map.equal_range(id);
						for (auto res = range.first; res != range.second; res++) {
							connect_cand.push_back(res->second);
						}
					}
				}
			}
		}

		connect_cand = connect_track(*itr, connect_cand, param);

		if (connect_cand.size() == 0)continue;
		for (auto itr2 = connect_cand.begin(); itr2 != connect_cand.end(); itr2++) {
			connected.push_back(std::make_pair((*itr2), &(*itr)));
		}
	}
	printf("\r connect track %10d/%10d (%4.1lf%%)\n", count, base1.size(), count * 100. / base1.size());

	return connected;

}
std::vector<vxx::base_track_t*> connect_track(vxx::base_track_t& t, std::vector<vxx::base_track_t*>& connect_cand, t2l_param& param) {
	std::vector<vxx::base_track_t*> ret;
	for (auto itr = connect_cand.begin(); itr != connect_cand.end(); itr++) {
		//if (judge_connect_xy(t, *(*itr), param) && (judge_connect_rl(t, *(*itr), param)||judge_connect_pb(t, *(*itr), param))) {
		if (judge_connect_xy(*(*itr), t, param) && (judge_connect_rl(*(*itr), t, param))) {
			ret.push_back(*itr);
		}
	}
	return ret;
}
bool judge_connect_xy(vxx::base_track_t& t1, vxx::base_track_t& t2, t2l_param& param) {
	double angle, d_pos_x, d_pos_y, d_ang_x, d_ang_y;
	double all_pos_x, all_pos_y, all_ang_x, all_ang_y;

	all_ang_x = param.intercept_ax + param.slope_ax * fabs(t1.ax) + param.slope2_ax * pow(t1.ax, 2);
	all_ang_y = param.intercept_ay + param.slope_ay * fabs(t1.ay) + param.slope2_ay * pow(t1.ay, 2);

	all_pos_x = param.intercept_px + param.slope_px * fabs(t1.ax) + param.slope2_px * pow(t1.ax, 2);
	all_pos_y = param.intercept_py + param.slope_py * fabs(t1.ay) + param.slope2_py * pow(t1.ay, 2);

	d_pos_x = t2.x - t1.x - (t1.ax + t2.ax) / 2 * (t2.z - t1.z);
	d_pos_y = t2.y - t1.y - (t1.ay + t2.ay) / 2 * (t2.z - t1.z);

	d_ang_x = (t2.ax - t1.ax);
	d_ang_y = (t2.ay - t1.ay);

	if (fabs(d_ang_x) > fabs(all_ang_x))return false;
	if (fabs(d_ang_y) > fabs(all_ang_y))return false;

	if (fabs(d_pos_x) > fabs(all_pos_x))return false;
	if (fabs(d_pos_y) > fabs(all_pos_y))return false;

	return true;
}
bool judge_connect_rl(vxx::base_track_t& t1, vxx::base_track_t& t2, t2l_param& param) {
	double angle, d_pos_r, d_pos_l, d_ang_r, d_ang_l;
	double all_pos_r, all_pos_l, all_ang_r, all_ang_l;
	angle = sqrt(t1.ax * t1.ax + t1.ay * t1.ay);
	if (angle < 0.01)return true;

	all_ang_r = param.intercept_ar + param.slope_ar * angle + param.slope2_ar * angle * angle;
	all_ang_l = param.intercept_al + param.slope_al * angle + param.slope2_al * angle * angle;
	all_pos_r = param.intercept_pr + param.slope_pr * angle + param.slope2_pr * angle * angle;
	all_pos_l = param.intercept_pl + param.slope_pl * angle + param.slope2_pl * angle * angle;


	d_ang_r = ((t2.ax - t1.ax) * t1.ax + (t2.ay - t1.ay) * t1.ay);
	d_ang_l = ((t2.ax - t1.ax) * t1.ay - (t2.ay - t1.ay) * t1.ax);
	if (fabs(d_ang_r) > fabs(all_ang_r) * angle)return false;
	if (fabs(d_ang_l) > fabs(all_ang_l) * angle)return false;

	Calc_position_difference(t1, t2, d_pos_r, d_pos_l);
	if (fabs(d_pos_r) > fabs(all_pos_r))return false;
	if (fabs(d_pos_l) > fabs(all_pos_l))return false;

	//printf("angle %.4lf\n", angle);
	//printf("angle    radial : %.4lf %.4lf\n", d_ang_r / angle, all_ang_r);
	//printf("angle    lateral: %.4lf %.4lf\n", d_ang_l / angle, all_ang_l);
	//printf("position radial : %.4lf %.4lf\n", d_pos_r, all_pos_r);
	//printf("position lateral: %.4lf %.4lf\n", d_pos_l, all_pos_l);
	//printf("\n");

	return true;
}
bool judge_connect_pb(vxx::base_track_t& t1, vxx::base_track_t& t2, t2l_param& param) {
	//double angle, d_pos_r, d_pos_l, d_ang_r, d_ang_l;
	//double all_pos_r, all_pos_l, all_ang_r, all_ang_l;
	//angle = sqrt(t1.ax*t1.ax + t1.ay*t1.ay);
	//if (angle < 0.1)return true;

	//all_ang_r = param.intercept_ar + param.slope_ar*angle;
	//all_ang_l = param.intercept_al + param.slope_al*angle;
	//all_pos_r = param.intercept_pr + param.slope_pr*angle;
	//all_pos_l = param.intercept_pl + param.slope_pl*angle;


	//d_ang_r = ((t2.ax - t1.ax)*t1.ax + (t2.ay - t1.ay)*t1.ay);
	//d_ang_l = ((t2.ax - t1.ax)*t1.ay - (t2.ay - t1.ay)*t1.ax);
	//if (fabs(d_ang_r) > fabs(all_ang_r)*angle)return false;
	//if (fabs(d_ang_l) > fabs(all_ang_l)*angle)return false;

	//Calc_position_difference(t1, t2, d_pos_r, d_pos_l);
	//if (fabs(d_pos_r) > fabs(all_pos_r))return false;
	//if (fabs(d_pos_l) > fabs(all_pos_l))return false;

	//return true;
	return false;

}
void Calc_position_difference(vxx::base_track_t& t1, vxx::base_track_t& t2, double& dr, double& dl) {
	using namespace matrix_3D;
	vector_3D pos0, pos1, dir0, dir1;
	pos0.x = t1.x;
	pos0.y = t1.y;
	pos0.z = t1.z;
	dir0.x = t1.ax;
	dir0.y = t1.ay;
	dir0.z = 1;
	pos1.x = t2.x;
	pos1.y = t2.y;
	pos1.z = t2.z;
	dir1.x = t2.ax;
	dir1.y = t2.ay;
	dir1.z = 1;

	vector_3D base_point, difference;
	//外挿基準点を1:1に内分した点に設定
	base_point = addition(const_multiple(pos0, 0.5), const_multiple(pos1, 0.5));
	difference = addition(const_multiple(pos0, -1), pos1);

	vector_3D extra0, extra1;
	double ratio0, ratio1;
	ratio0 = -1 * dot(addition(pos0, const_multiple(base_point, -1)), difference) / dot(dir0, difference);
	ratio1 = -1 * dot(addition(pos1, const_multiple(base_point, -1)), difference) / dot(dir1, difference);
	extra0 = addition(pos0, const_multiple(dir0, ratio0));
	extra1 = addition(pos1, const_multiple(dir1, ratio1));

	vector_3D unit_r, unit_l;
	unit_l.x = -1 + difference.y;
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

void output_corrmap2(std::string filename, std::vector<align_param2>& corr) {

	std::ofstream ofs(filename);
	if (!ofs) {
		//file open 失敗
		fprintf(stderr, "File[%s] is not exist!!\n", filename.c_str());
		exit(1);
	}
	if (corr.size() == 0) {
		fprintf(stderr, "target linklet ... null\n");
		fprintf(stderr, "File[%s] has no text\n", filename.c_str());
	}
	else {
		int count = 0;
		std::cout << std::right << std::fixed;
		for (auto itr = corr.begin(); itr != corr.end(); itr++) {
			if (count % 10000 == 0) {
				fprintf(stderr, "\r Write corrmap2 ... %d/%d (%4.1lf%%)", count, int(corr.size()), count * 100. / corr.size());
			}
			count++;
			ofs << std::right << std::fixed
				<< std::setw(8) << std::setprecision(0) << itr->corr_p[0]->id << " "
				<< std::setw(8) << std::setprecision(0) << itr->corr_p[1]->id << " "
				<< std::setw(8) << std::setprecision(0) << itr->corr_p[2]->id << " "
				<< std::setw(8) << std::setprecision(1) << itr->x << " "
				<< std::setw(8) << std::setprecision(1) << itr->y << " "
				<< std::setw(8) << std::setprecision(1) << itr->z << " "
				<< std::setw(9) << std::setprecision(7) << itr->x_rot << " "
				<< std::setw(9) << std::setprecision(7) << itr->y_rot << " "
				<< std::setw(9) << std::setprecision(7) << itr->z_rot << " "
				<< std::setw(8) << std::setprecision(6) << itr->x_shrink << " "
				<< std::setw(8) << std::setprecision(6) << itr->y_shrink << " "
				<< std::setw(8) << std::setprecision(6) << itr->z_shrink << " "
				<< std::setw(8) << std::setprecision(6) << itr->yx_shear << " "
				<< std::setw(8) << std::setprecision(6) << itr->zx_shear << " "
				<< std::setw(8) << std::setprecision(6) << itr->zy_shear << " "
				<< std::setw(8) << std::setprecision(2) << itr->dx << " "
				<< std::setw(8) << std::setprecision(2) << itr->dy << " "
				<< std::setw(8) << std::setprecision(2) << itr->dz << std::endl;
		}
		fprintf(stderr, "\r Write corrmap2 ... %d/%d (%4.1lf%%)\n", count, int(corr.size()), count * 100. / corr.size());
	}



}

void output_base_corrmap_pair(std::string filename, std::vector <std::pair<vxx::base_track_t*, align_param2*>>& pair_v) {

	std::ofstream ofs(filename);
	if (!ofs) {
		//file open 失敗
		fprintf(stderr, "File[%s] is not exist!!\n", filename.c_str());
		exit(1);
	}
	if (pair_v.size() == 0) {
		fprintf(stderr, "target linklet ... null\n");
		fprintf(stderr, "File[%s] has no text\n", filename.c_str());
	}
	else {
		int count = 0;
		std::cout << std::right << std::fixed;
		for (auto itr = pair_v.begin(); itr != pair_v.end(); itr++) {
			if (count % 10000 == 0) {
				fprintf(stderr, "\r Write base corr pair ... %d/%d (%4.1lf%%)", count, int(pair_v.size()), count * 100. / pair_v.size());
			}
			count++;
			ofs << std::right << std::fixed
				<< std::setw(5) << std::setprecision(0) << itr->first->m[0].pos << " "
				<< std::setw(12) << std::setprecision(0) << itr->first->m[0].rawid << " "
				<< std::setw(5) << std::setprecision(0) << itr->first->m[1].pos << " "
				<< std::setw(12) << std::setprecision(0) << itr->first->m[1].rawid << " "
				<< std::setw(8) << std::setprecision(0) << itr->second->corr_p[0]->id << " "
				<< std::setw(8) << std::setprecision(0) << itr->second->corr_p[1]->id << " "
				<< std::setw(8) << std::setprecision(0) << itr->second->corr_p[2]->id << std::endl;
		}
		fprintf(stderr, "\r Write base corr pair ... %d/%d (%4.1lf%%)\n", count, int(pair_v.size()), count * 100. / pair_v.size());

	}
}

void output_pair_txt(std::string filename, std::vector<output_format_link>& link) {
	std::ofstream ofs(filename);
	int count = 0;
	for (auto l : link) {
		if (count % 10000 == 0) {
			printf("\r write linklet %10d/%10d(%4.1lf%%)", count, link.size(), count * 100. / link.size());
		}
		count++;

		ofs << std::right << std::fixed;
		for (int i = 0; i < 2; i++) {
			ofs << std::setw(4) << std::setprecision(0) << l.b[i].pl << " "
				<< std::setw(10) << std::setprecision(0) << l.b[i].rawid << " "
				<< std::setw(7) << std::setprecision(4) << l.b[i].ax << " "
				<< std::setw(7) << std::setprecision(4) << l.b[i].ay << " "
				<< std::setw(8) << std::setprecision(1) << l.b[i].x << " "
				<< std::setw(8) << std::setprecision(1) << l.b[i].y << " "
				<< std::setw(8) << std::setprecision(1) << l.b[i].z << " "
				<< std::setw(2) << std::setprecision(0) << l.b[i].m[0].ph << " "
				<< std::setw(2) << std::setprecision(0) << l.b[i].m[1].ph << " "
				<< std::setw(4) << std::setprecision(0) << l.b[i].m[0].vph << " "
				<< std::setw(4) << std::setprecision(0) << l.b[i].m[1].vph << " "
				<< std::setw(6) << std::setprecision(0) << l.b[i].m[0].px << " "
				<< std::setw(6) << std::setprecision(0) << l.b[i].m[1].px << " ";
		}
		ofs << std::setw(7) << std::setprecision(4) << l.dax << " "
			<< std::setw(7) << std::setprecision(4) << l.day << " "
			<< std::setw(6) << std::setprecision(1) << l.dx << " "
			<< std::setw(6) << std::setprecision(1) << l.dy << " "
			<< std::setw(7) << std::setprecision(4) << l.dar << " "
			<< std::setw(7) << std::setprecision(4) << l.dal << " "
			<< std::setw(6) << std::setprecision(1) << l.dr << " "
			<< std::setw(6) << std::setprecision(1) << l.dl << std::endl;
	}
	printf("\r write linklet %10d/%10d(%4.1lf%%)\n", count, link.size(), count * 100. / link.size());

}
void output_pair_bin(std::string filename, std::vector<output_format_link>& link) {

	std::ofstream ofs(filename, std::ios::binary);
	int count = 0;
	for (auto l : link) {
		if (count % 10000 == 0) {
			printf("\r write linklet %10d/%10d(%4.1lf%%)", count, link.size(), count * 100. / link.size());
		}
		count++;
		ofs.write((char*)&l, sizeof(output_format_link));
	}

	printf("\r write linklet %10d/%10d(%4.1lf%%)\n", count, link.size(), count * 100. / link.size());

}


output_format_link output_format(vxx::base_track_t& t1, vxx::base_track_t& t2) {
	int NumberOfImager = 72;
	uint32_t ShotID;
	output_format_link l;
	l.b[0].pl = t1.pl;
	l.b[0].rawid = t1.rawid;
	l.b[0].ax = t1.ax;
	l.b[0].ay = t1.ay;
	l.b[0].x = t1.x;
	l.b[0].y = t1.y;
	l.b[0].z = t1.z;

	l.b[0].m[0].pos = t1.m[0].pos;
	l.b[0].m[0].zone = t1.m[0].zone;
	l.b[0].m[0].isg = t1.m[0].isg;
	l.b[0].m[0].ph = t1.m[0].ph / 10000;
	l.b[0].m[0].vph = t1.m[0].ph % 10000;
	ShotID = ((uint32_t)(uint16_t)t1.m[0].row << 16) | ((uint32_t)(uint16_t)t1.m[0].col);
	l.b[0].m[0].view = ShotID / NumberOfImager;
	l.b[0].m[0].imager = ShotID % NumberOfImager;
	l.b[0].m[0].px = -1;

	l.b[0].m[1].pos = t1.m[1].pos;
	l.b[0].m[1].zone = t1.m[1].zone;
	l.b[0].m[1].isg = t1.m[1].isg;
	l.b[0].m[1].ph = t1.m[1].ph / 10000;
	l.b[0].m[1].vph = t1.m[1].ph % 10000;
	ShotID = ((uint32_t)(uint16_t)t1.m[1].row << 16) | ((uint32_t)(uint16_t)t1.m[1].col);
	l.b[0].m[1].view = ShotID / NumberOfImager;
	l.b[0].m[1].imager = ShotID % NumberOfImager;
	l.b[0].m[1].px = -1;

	l.b[1].pl = t2.pl;
	l.b[1].rawid = t2.rawid;
	l.b[1].ax = t2.ax;
	l.b[1].ay = t2.ay;
	l.b[1].x = t2.x;
	l.b[1].y = t2.y;
	l.b[1].z = t2.z;

	l.b[1].m[0].pos = t2.m[0].pos;
	l.b[1].m[0].zone = t2.m[0].zone;
	l.b[1].m[0].isg = t2.m[0].isg;
	l.b[1].m[0].ph = t2.m[0].ph / 10000;
	l.b[1].m[0].vph = t2.m[0].ph % 10000;
	ShotID = ((uint32_t)(uint16_t)t2.m[0].row << 16) | ((uint32_t)(uint16_t)t2.m[0].col);
	l.b[1].m[0].view = ShotID / NumberOfImager;
	l.b[1].m[0].imager = ShotID % NumberOfImager;
	l.b[1].m[0].px = -1;

	l.b[1].m[1].pos = t2.m[1].pos;
	l.b[1].m[1].zone = t2.m[1].zone;
	l.b[1].m[1].isg = t2.m[1].isg;
	l.b[1].m[1].ph = t2.m[1].ph / 10000;
	l.b[1].m[1].vph = t2.m[1].ph % 10000;
	ShotID = ((uint32_t)(uint16_t)t2.m[1].row << 16) | ((uint32_t)(uint16_t)t2.m[1].col);
	l.b[1].m[1].view = ShotID / NumberOfImager;
	l.b[1].m[1].imager = ShotID % NumberOfImager;
	l.b[1].m[1].px = -1;

	l.Calc_difference();

	return l;
}
void output_format_link::Calc_difference() {

	dax = b[1].ax - b[0].ax;
	day = b[1].ay - b[0].ay;
	dx = b[1].x - b[0].x - (b[0].ax + b[1].ax) / 2 * (b[1].z - b[0].z);
	dy = b[1].y - b[0].y - (b[0].ay + b[1].ay) / 2 * (b[1].z - b[0].z);
	dar = (dax * b[0].ax + day * b[0].ay) / sqrt(b[0].ax * b[0].ax + b[0].ay * b[0].ay);
	dal = (dax * b[0].ay - day * b[0].ax) / sqrt(b[0].ax * b[0].ax + b[0].ay * b[0].ay);


	//dr,dlの計算
	using namespace matrix_3D;
	vector_3D pos0, pos1, dir0, dir1;
	pos0.x = b[0].x;
	pos0.y = b[0].y;
	pos0.z = b[0].z;
	dir0.x = b[0].ax;
	dir0.y = b[0].ay;
	dir0.z = 1;
	pos1.x = b[1].x;
	pos1.y = b[1].y;
	pos1.z = b[1].z;
	dir1.x = b[1].ax;
	dir1.y = b[1].ay;
	dir1.z = 1;

	vector_3D base_point, difference;
	//外挿基準点を1:1に内分した点に設定
	base_point = addition(const_multiple(pos0, 0.5), const_multiple(pos1, 0.5));
	difference = addition(const_multiple(pos0, -1), pos1);

	vector_3D extra0, extra1;
	double ratio0, ratio1;
	ratio0 = -1 * dot(addition(pos0, const_multiple(base_point, -1)), difference) / dot(dir0, difference);
	ratio1 = -1 * dot(addition(pos1, const_multiple(base_point, -1)), difference) / dot(dir1, difference);
	extra0 = addition(pos0, const_multiple(dir0, ratio0));
	extra1 = addition(pos1, const_multiple(dir1, ratio1));

	vector_3D unit_r, unit_l;
	unit_l.x = -1 + difference.y;
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
std::vector<output_format_link> basetrack_pair_to_linket(std::vector<std::pair<vxx::base_track_t*, vxx::base_track_t*>>& track_pair) {
	std::vector<output_format_link> link;

	for (auto itr = track_pair.begin(); itr != track_pair.end(); itr++) {
		link.push_back(output_format(*(itr->first), *(itr->second)));
	}

	return link;
}

int64_t file_size(std::string filename) {
	std::ifstream ifs(filename, std::ios::binary);
	//filesize取得
	ifs.seekg(0, std::ios::end);
	int64_t eofpos = ifs.tellg();
	ifs.clear();
	ifs.seekg(0, std::ios::beg);
	int64_t begpos = ifs.tellg();
	int64_t nowpos = ifs.tellg();
	int64_t size2 = eofpos - begpos;
	return size2;
}

std::vector<microtrack_inf> read_microtrack_inf(std::string filename, bool output) {

	int64_t micro_num = file_size(filename) / sizeof(microtrack_inf);

	FILE* fp_in;
	if ((fp_in = fopen(filename.c_str(), "rb")) == NULL) {
		printf("%s file not open!\n", filename.c_str());
		exit(EXIT_FAILURE);
	}

	std::vector<microtrack_inf> ret;
	ret.reserve(micro_num);
	const int Read_Block = 10000;

	microtrack_inf m_buf[Read_Block];
	int64_t  now = 0;
	int read_num;
	bool flg = true;
	while (flg) {
		if (micro_num - now == Read_Block) {
			read_num = micro_num - now;
			flg = false;
		}
		else if (micro_num - now < Read_Block) {
			read_num = micro_num - now;
			flg = false;
		}
		else {
			read_num = Read_Block;
		}
		fread(&m_buf, sizeof(microtrack_inf), read_num, fp_in);
		for (int i = 0; i < read_num; i++) {
			ret.push_back(m_buf[i]);
		}
		now += read_num;
	}
	return ret;
}

void zone_trans(std::map<std::tuple<int, int, int, int, int>, int>& pixel_inf) {
	int pos, zone;
	std::tuple<int, int, int, int, int> newid;
	std::vector<std::pair<std::tuple<int, int, int, int, int>, int>> add_pixel_inf;
	add_pixel_inf.reserve(pixel_inf.size());

	for (auto itr = pixel_inf.begin(); itr != pixel_inf.end(); itr++) {
		pos = std::get<0>(itr->first);
		zone = std::get<1>(itr->first);
		newid = itr->first;
		if (pos % 10 == 1) {
			if ((zone - 1) / 6 == 0 || (zone - 1) / 6 == 2 || (zone - 1) / 6 == 4 || (zone - 1) / 6 == 6) {
				std::get<1>(newid) = zone + 6;
				add_pixel_inf.push_back(std::make_pair(newid, itr->second));
			}
			else {
				fprintf(stderr, "file format exception\n");
				fprintf(stderr, "pos %d zone %d\n", pos, zone);
				exit(1);
			}
		}
		else if (pos % 10 == 2) {
			if ((zone - 1) / 6 == 0 || (zone - 1) / 6 == 1 || (zone - 1) / 6 == 4 || (zone - 1) / 6 == 5) {
				std::get<1>(newid) = zone + 12;
				add_pixel_inf.push_back(std::make_pair(newid, itr->second));
			}
			else {
				fprintf(stderr, "file format exception\n");
				fprintf(stderr, "pos %d zone %d\n", pos, zone);
				exit(1);
			}
		}
	}
	int count = 0;
	for (auto itr = add_pixel_inf.begin(); itr != add_pixel_inf.end(); itr++) {
		if (count % 1000000 == 0) {
			fprintf(stderr, "\r pixel inf foramt change %d/%d (%4.1lf%%)", count, add_pixel_inf.size(), count * 100. / add_pixel_inf.size());
		}
		count++;

		pixel_inf.insert(*itr);
	}
	fprintf(stderr, "\r pixel inf foramt change %d/%d (%4.1lf%%)\n", count, add_pixel_inf.size(), count * 100. / add_pixel_inf.size());
}

//高速化の実装-->OK
void Pick_up_pixel_count(std::vector<output_format_link>& link, std::string file_in_ECC) {
	if (link.size() == 0)return;
	int pl[2] = { link.begin()->b[0].pl,link.begin()->b[1].pl };

	const int NumberOfImager = 72;
	uint32_t ShotID;
	int view, imager;

	std::tuple<int, int, int, int, int> id;
	//linkletのpxへのポインタをmapのvalにする
	std::multimap<std::tuple<int, int, int, int, int>, int*> link_px;
	for (auto itr = link.begin(); itr != link.end(); itr++) {
		for (int ib = 0; ib < 2; ib++) {
			for (int im = 0; im < 2; im++) {
				std::get<0>(id) = itr->b[ib].m[im].pos;
				std::get<1>(id) = itr->b[ib].m[im].zone;
				std::get<2>(id) = itr->b[ib].m[im].view;
				std::get<3>(id) = itr->b[ib].m[im].imager;
				std::get<4>(id) = itr->b[ib].m[im].isg;

				link_px.insert(std::make_pair(id, &(itr->b[ib].m[im].px)));
			}
		}
	}


	for (int ipl = 0; ipl <= 1; ipl++) {
		std::set<int>zone_pos1, zone_pos2;
		for (int Area = 1; Area <= 6; Area++) {
			for (int iscan = 0; iscan < 4; iscan++) {
				printf("PL%03d Area%d scan %d", pl[ipl], Area, iscan);
				//読み込むfile名
				std::stringstream infile_thick, infile_thin;
				infile_thick << file_in_ECC << "\\Area"
					<< std::setw(1) << Area << "\\PL"
					<< std::setw(3) << std::setfill('0') << pl[ipl] << "\\micro_inf_thick_"
					<< std::setw(1) << iscan;
				infile_thin << file_in_ECC << "\\Area"
					<< std::setw(1) << Area << "\\PL"
					<< std::setw(3) << std::setfill('0') << pl[ipl] << "\\micro_inf_thin_"
					<< std::setw(1) << iscan;
				//fileの存在確認
				if (!std::filesystem::exists(infile_thick.str())) {
					fprintf(stderr, "%s not exist\n", infile_thick.str().c_str());
					//exit(1);
					continue;
				}
				if (!std::filesystem::exists(infile_thin.str())) {
					fprintf(stderr, "%s not exist\n", infile_thin.str().c_str());
					//exit(1);
					continue;
				}
				printf(" thick %d , thin %d\n", file_size(infile_thick.str()) / sizeof(microtrack_inf), file_size(infile_thin.str()) / sizeof(microtrack_inf));
				//読み込み&pxの取得
				add_microtrack_inf_format0(infile_thick.str(), link_px, zone_pos1, zone_pos2);
				add_microtrack_inf_format0(infile_thin.str(), link_px, zone_pos1, zone_pos2);
			}
		}
		if (zone_pos1.size() == 48 && zone_pos2.size() == 48) {
			continue;
		}
		else if (zone_pos1.size() == 24 && zone_pos2.size() == 24) {
			for (int Area = 1; Area <= 6; Area++) {
				for (int iscan = 0; iscan < 4; iscan++) {
					printf("PL%03d Area%d scan %d", pl[ipl], Area, iscan);
					//読み込むfile名
					std::stringstream infile_thick, infile_thin;
					infile_thick << file_in_ECC << "\\Area"
						<< std::setw(1) << Area << "\\PL"
						<< std::setw(3) << std::setfill('0') << pl[ipl] << "\\micro_inf_thick_"
						<< std::setw(1) << iscan;
					infile_thin << file_in_ECC << "\\Area"
						<< std::setw(1) << Area << "\\PL"
						<< std::setw(3) << std::setfill('0') << pl[ipl] << "\\micro_inf_thin_"
						<< std::setw(1) << iscan;
					//fileの存在確認
					if (!std::filesystem::exists(infile_thick.str())) {
						fprintf(stderr, "%s not exist\n", infile_thick.str().c_str());
						exit(1);
					}
					if (!std::filesystem::exists(infile_thin.str())) {
						fprintf(stderr, "%s not exist\n", infile_thin.str().c_str());
						exit(1);
					}
					printf(" thick %d , thin %d\n", file_size(infile_thick.str()) / sizeof(microtrack_inf), file_size(infile_thin.str()) / sizeof(microtrack_inf));
					//読み込み&pxの取得
					add_microtrack_inf_format1(infile_thick.str(), link_px);
					add_microtrack_inf_format1(infile_thin.str(), link_px);
				}
			}

		}
		else {
			fprintf(stderr, "file format exception\n");
			fprintf(stderr, "zone1 size = %d zone2 size=%d\n", zone_pos1.size(), zone_pos2.size());
			//exit(1);
		}
	}

}

void add_microtrack_inf_format0(std::string filename, std::multimap<std::tuple<int, int, int, int, int>, int*>& link_px, std::set<int>& zone_pos1, std::set<int>& zone_pos2) {
	const int NumberOfImager = 72;
	uint32_t ShotID;
	int view, imager;

	int64_t micro_num = file_size(filename) / sizeof(microtrack_inf);

	FILE* fp_in;
	if ((fp_in = fopen(filename.c_str(), "rb")) == NULL) {
		printf("%s file not open!\n", filename.c_str());
		exit(EXIT_FAILURE);
	}

	std::vector<microtrack_inf> ret;
	ret.reserve(micro_num);
	const int Read_Block = 10000;
	std::tuple<int, int, int, int, int> id;

	microtrack_inf m_buf[Read_Block];
	int64_t  now = 0;
	int read_num;
	bool flg = true;
	while (flg) {
		if (micro_num - now == Read_Block) {
			read_num = micro_num - now;
			flg = false;
		}
		else if (micro_num - now < Read_Block) {
			read_num = micro_num - now;
			flg = false;
		}
		else {
			read_num = Read_Block;
		}
		fread(&m_buf, sizeof(microtrack_inf), read_num, fp_in);
		for (int i = 0; i < read_num; i++) {
			std::get<0>(id) = m_buf[i].pos;
			std::get<1>(id) = m_buf[i].zone;
			ShotID = ((uint32_t)(uint16_t)m_buf[i].row << 16) | ((uint32_t)(uint16_t)m_buf[i].col);
			view = ShotID / NumberOfImager;
			imager = ShotID % NumberOfImager;
			std::get<2>(id) = view;
			std::get<3>(id) = imager;
			std::get<4>(id) = m_buf[i].isg;
			if (m_buf[i].pos % 10 == 1)zone_pos1.insert(m_buf[i].zone);
			if (m_buf[i].pos % 10 == 2)zone_pos2.insert(m_buf[i].zone);
			if (link_px.count(id) == 0)continue;
			auto range = link_px.equal_range(id);
			for (auto res = range.first; res != range.second; res++) {
				(*res->second) = m_buf[i].hitnum;
			}
		}
		now += read_num;
	}
}

void add_microtrack_inf_format1(std::string filename, std::multimap<std::tuple<int, int, int, int, int>, int*>& link_px) {
	const int NumberOfImager = 72;
	uint32_t ShotID;
	int view, imager;

	int64_t micro_num = file_size(filename) / sizeof(microtrack_inf);

	FILE* fp_in;
	if ((fp_in = fopen(filename.c_str(), "rb")) == NULL) {
		printf("%s file not open!\n", filename.c_str());
		exit(EXIT_FAILURE);
	}

	std::vector<microtrack_inf> ret;
	ret.reserve(micro_num);
	const int Read_Block = 10000;
	std::tuple<int, int, int, int, int> id;

	microtrack_inf m_buf[Read_Block];
	int64_t  now = 0;
	int read_num;
	bool flg = true;
	while (flg) {
		if (micro_num - now == Read_Block) {
			read_num = micro_num - now;
			flg = false;
		}
		else if (micro_num - now < Read_Block) {
			read_num = micro_num - now;
			flg = false;
		}
		else {
			read_num = Read_Block;
		}
		fread(&m_buf, sizeof(microtrack_inf), read_num, fp_in);
		for (int i = 0; i < read_num; i++) {
			std::get<0>(id) = m_buf[i].pos;
			std::get<1>(id) = m_buf[i].zone;
			ShotID = ((uint32_t)(uint16_t)m_buf[i].row << 16) | ((uint32_t)(uint16_t)m_buf[i].col);
			view = ShotID / NumberOfImager;
			imager = ShotID % NumberOfImager;
			std::get<2>(id) = view;
			std::get<3>(id) = imager;
			std::get<4>(id) = m_buf[i].isg;
			if (link_px.count(id) != 0) {
				auto range = link_px.equal_range(id);
				for (auto res = range.first; res != range.second; res++) {
					(*res->second) = m_buf[i].hitnum;
				}
			}
			if (m_buf[i].pos % 10 == 1) {
				if ((m_buf[i].zone - 1) / 6 == 0 || (m_buf[i].zone - 1) / 6 == 2 || (m_buf[i].zone - 1) / 6 == 4 || (m_buf[i].zone - 1) / 6 == 6) {
					std::get<1>(id) = m_buf[i].zone + 6;
				}
				else {
					fprintf(stderr, "file format exception\n");
					fprintf(stderr, "pos %d zone %d\n", m_buf[i].pos, m_buf[i].zone);
					exit(1);
				}
			}
			else if (m_buf[i].pos % 10 == 2) {
				if ((m_buf[i].zone - 1) / 6 == 0 || (m_buf[i].zone - 1) / 6 == 1 || (m_buf[i].zone - 1) / 6 == 4 || (m_buf[i].zone - 1) / 6 == 5) {
					std::get<1>(id) = m_buf[i].zone + 12;
				}
				else {
					fprintf(stderr, "file format exception\n");
					fprintf(stderr, "pos %d zone %d\n", m_buf[i].pos, m_buf[i].zone);
					exit(1);
				}
			}
			if (link_px.count(id) != 0) {
				auto range = link_px.equal_range(id);
				for (auto res = range.first; res != range.second; res++) {
					(*res->second) = m_buf[i].hitnum;
				}
			}
		}
		now += read_num;
	}
}

