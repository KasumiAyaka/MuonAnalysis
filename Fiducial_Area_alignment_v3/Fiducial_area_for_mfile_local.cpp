#define _CRT_SECURE_NO_WARNINGS
#pragma comment(lib, "VxxReader.lib")
#include "VxxReader.h"
#pragma comment(lib,"FILE_structure.lib")
//#include <FILE_structure.h>
#include "File_structure.h"
#include <filesystem>
#include <set>
class Point {
public:
	double x, y, z;
};
class Fiducial_Area {
public:
	int pl;
	Point p[2];
	//double x0, y0, z0, x1, y1, z1;
};


std::map<int, std::vector<Fiducial_Area>> read_fiducial_Area(std::string filename);
void trans_mfile_cordinate(std::vector<corrmap0::Corrmap>& corr_abs, std::map<int, std::vector<Fiducial_Area>>& area, std::map<int, double>& z_map);
void output_Fiducial_Area_line(std::string fileanme, std::map<int, std::vector<Fiducial_Area>>& area);
void trans_base_all(std::vector < std::pair<Point*, corrmap_3d::align_param2*>>& track_pair);
std::vector <std::pair<Point*, corrmap_3d::align_param2*>>track_affineparam_correspondence(std::vector<Point*>& p, std::vector <corrmap_3d::align_param2>& param);
void trans_mfile_cordinate(std::vector<corrmap_3d::align_param2>& param, std::vector<Fiducial_Area>& area);

int main(int argc, char** argv) {
	if (argc != 4) {
		fprintf(stderr, "usage: prg fa.txt(txt) file-in-ECC output-path\n");
		exit(1);
	}

	std::string file_in_fa = argv[1];
	std::string file_in_ECC = argv[2];
	std::string file_out_fa = argv[3];

	//corrmap absの読み込み
	std::string file_in_corrmap = file_in_ECC + "\\Area0\\0\\align\\fine\\local\\corrmap-local-abs.lst";
	std::map<int, std::vector<corrmap_3d::align_param>>corrmap = corrmap_3d::read_ali_param_abs(file_in_corrmap, 1);
	std::map<int, std::vector<corrmap_3d::align_param2>>corrmap_dd = corrmap_3d::DelaunayDivide_map(corrmap);

	//std::vector<corrmap0::Corrmap> corr_abs;
	//corrmap0::read_cormap(file_in_corrabs, corr_abs);
	//fiducial areaの読み込み
	std::map<int, std::vector<Fiducial_Area>> area = read_fiducial_Area(file_in_fa);
	//gap nominal read
	//std::stringstream structure_path;
	//structure_path << file_in_ECC << "\\st\\st.dat";
	//chamber1::Chamber chamber;
	//chamber1::read_structure(structure_path.str(), chamber);
	//std::map<int, double> z_map = chamber1::base_z_convert(chamber);

	//corrmap absの適用
	for (auto itr = area.begin(); itr != area.end(); itr++) {
		if (corrmap_dd.count(itr->first) == 0) {
			fprintf(stderr, "corrmap local abs PL%03d not found\n", itr->first);
			exit(1);
		}
		std::vector<corrmap_3d::align_param2> param = corrmap_dd.at(itr->first);
		trans_mfile_cordinate(param, itr->second);
	}

	output_Fiducial_Area_line(file_out_fa, area);
}
std::map<int, std::vector<Fiducial_Area>> read_fiducial_Area(std::string filename) {

	std::ifstream ifs(filename);
	std::multimap<int, Fiducial_Area> fa_multi;
	std::map<int, std::vector<Fiducial_Area>> ret;
	Fiducial_Area fa;
	while (ifs >> fa.pl >> fa.p[0].x >> fa.p[0].y >> fa.p[0].z >> fa.p[1].x >> fa.p[1].y >> fa.p[1].z) {
		fa_multi.insert(std::make_pair(fa.pl, fa));
	}

	int count = 0;
	for (auto itr = fa_multi.begin(); itr != fa_multi.end(); itr++) {
		count = fa_multi.count(itr->first);
		auto range = fa_multi.equal_range(itr->first);
		std::vector<Fiducial_Area> fa_vec;
		for (auto itr2 = range.first; itr2 != range.second; itr2++) {
			fa_vec.push_back(itr2->second);
		}
		ret.insert(std::make_pair(itr->first, fa_vec));
	}

	return ret;

}


void trans_mfile_cordinate(std::vector<corrmap_3d::align_param2>& param, std::vector<Fiducial_Area>& area) {

	std::vector< Point*> p_trans;
	for (auto itr = area.begin(); itr != area.end(); itr++) {
		p_trans.push_back(&(itr->p[0]));
		p_trans.push_back(&(itr->p[1]));
	}
	std::vector <std::pair<Point*, corrmap_3d::align_param2*>> p_trans_map = track_affineparam_correspondence(p_trans, param);
	trans_base_all(p_trans_map);
}


void output_Fiducial_Area_line(std::string fileanme, std::map<int, std::vector<Fiducial_Area>>& area) {
	std::ofstream ofs(fileanme);
	if (area.size() == 0) {
		fprintf(stderr, "target corrmap ... null\n");
	}
	else {
		int count = 0;
		std::cout << std::right << std::fixed;
		for (auto itr = area.begin(); itr != area.end(); itr++) {
			for (auto itr2 = itr->second.begin(); itr2 != itr->second.end(); itr2++) {



				ofs << std::right << std::fixed
					<< std::setw(5) << std::setprecision(0) << itr2->pl << " "
					<< std::setw(8) << std::setprecision(1) << itr2->p[0].x << " "
					<< std::setw(8) << std::setprecision(1) << itr2->p[0].y << " "
					<< std::setw(8) << std::setprecision(1) << itr2->p[0].z << " "
					<< std::setw(8) << std::setprecision(1) << itr2->p[1].x << " "
					<< std::setw(8) << std::setprecision(1) << itr2->p[1].y << " "
					<< std::setw(8) << std::setprecision(1) << itr2->p[1].z << std::endl;
			}

		}
	}
}



//mfile0::M_Base
//basetrack-alignment mapの対応
double select_triangle_vale(corrmap_3d::align_param2* param, Point* p) {
	double x, y;
	double dist = 0;
	x = (param->corr_p[0]->x + param->corr_p[1]->x + param->corr_p[2]->x) / 3;
	y = (param->corr_p[0]->y + param->corr_p[1]->y + param->corr_p[2]->y) / 3;
	dist = (p->x - x) * (p->x - x) + (p->y - y) * (p->y - y);
	return dist;
}
corrmap_3d::align_param2* search_param(std::vector<corrmap_3d::align_param*>& param, Point* p, std::multimap<int, corrmap_3d::align_param2*>& triangles) {
	//三角形内部
	//最近接三角形
	double dist = 0;
	std::map<double, corrmap_3d::align_param* > dist_map;
	//align_paramを近い順にsort
	for (auto itr = param.begin(); itr != param.end(); itr++) {
		dist = ((*itr)->x - p->x) * ((*itr)->x - p->x) + ((*itr)->y - p->y) * ((*itr)->y - p->y);
		dist_map.insert(std::make_pair(dist, (*itr)));
	}

	double sign[3];
	bool flg = false;
	int id;

	corrmap_3d::align_param2* ret = triangles.begin()->second;
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
			sign[0] = (itr2->second->corr_p[1]->x - itr2->second->corr_p[0]->x) * (p->y - itr2->second->corr_p[1]->y) - (itr2->second->corr_p[1]->y - itr2->second->corr_p[0]->y) * (p->x - itr2->second->corr_p[1]->x);
			sign[1] = (itr2->second->corr_p[2]->x - itr2->second->corr_p[1]->x) * (p->y - itr2->second->corr_p[2]->y) - (itr2->second->corr_p[2]->y - itr2->second->corr_p[1]->y) * (p->x - itr2->second->corr_p[2]->x);
			sign[2] = (itr2->second->corr_p[0]->x - itr2->second->corr_p[2]->x) * (p->y - itr2->second->corr_p[0]->y) - (itr2->second->corr_p[0]->y - itr2->second->corr_p[2]->y) * (p->x - itr2->second->corr_p[0]->x);
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
			if (dist<0 || dist>select_triangle_vale(itr2->second, p)) {
				dist = select_triangle_vale(itr2->second, p);
				ret = itr2->second;
			}
		}
	}
	//printf("point not in trianlge\n");
	return ret;
}
std::vector <std::pair<Point*, corrmap_3d::align_param2*>>track_affineparam_correspondence(std::vector<Point*>& p, std::vector <corrmap_3d::align_param2>& param) {

	//local alignの視野中心を取り出して、位置でhash
	//local alignの視野中心の作るdelaunay三角形をmapで対応

	std::map<int, corrmap_3d::align_param*> view_center;
	std::multimap<int, corrmap_3d::align_param2*>triangles;
	double xmin = 999999, ymin = 999999, hash = 2000;
	for (auto itr = param.begin(); itr != param.end(); itr++) {
		for (int i = 0; i < 3; i++) {
			view_center.insert(std::make_pair(itr->corr_p[i]->id, (itr->corr_p[i])));
			triangles.insert(std::make_pair(itr->corr_p[i]->id, &(*itr)));
			xmin = std::min(itr->corr_p[i]->x, xmin);
			ymin = std::min(itr->corr_p[i]->y, ymin);
		}
	}
	std::multimap<std::pair<int, int>, corrmap_3d::align_param*> view_center_hash;
	std::pair<int, int>id;
	for (auto itr = view_center.begin(); itr != view_center.end(); itr++) {
		id.first = int((itr->second->x - xmin) / hash);
		id.second = int((itr->second->y - ymin) / hash);
		view_center_hash.insert(std::make_pair(id, itr->second));
	}

	std::vector < std::pair<Point*, corrmap_3d::align_param2*>> ret;
	std::vector<corrmap_3d::align_param*> param_cand;
	int loop = 0, ix, iy, count = 0;
	for (auto itr = p.begin(); itr != p.end(); itr++) {
		if (count % 100000 == 0) {
			printf("\r search correspond triangles %d/%d(%4.1lf%%)", count, p.size(), count * 100. / p.size());
		}
		count++;
		ix = ((*itr)->x - xmin) / hash;
		iy = ((*itr)->y - ymin) / hash;
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
			std::cout << "loop : " << loop << std::endl;
		}
		corrmap_3d::align_param2* param2 = search_param(param_cand, *itr, triangles);
		ret.push_back(std::make_pair((*itr), param2));
	}
	printf("\r search correspond triangles %d/%d(%4.1lf%%)\n", count, p.size(), count * 100. / p.size());

	return ret;
}
//変換 zshrink補正-->9para変換
void trans_base(std::vector<Point*>& p, corrmap_3d::align_param2* param) {

	matrix_3D::matrix_33 x_rot_mat(0, param->x_rot), y_rot_mat(1, param->y_rot), z_rot_mat(2, param->z_rot), all_trans(0, 0), shear_mat(0, 0), shrink_mat(0, 0);

	shrink_mat.val[0][0] *= param->x_shrink;
	shrink_mat.val[1][1] *= param->y_shrink;
	//shrink_mat.val[2][2] *= param->z_shrink;
	shear_mat.val[0][1] = param->yx_shear;
	//shear_mat.val[0][2] = param->zx_shear;
	//shear_mat.val[1][2] = param->zy_shear;

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
	matrix_3D::vector_3D base_p0;
	double base_thick = 210;
	for (auto itr = p.begin(); itr != p.end(); itr++) {
		base_p0.x = (*itr)->x;
		base_p0.y = (*itr)->y;
		base_p0.z = param->z;

		//base_p1.x = (*itr)->x + (*itr)->ax*base_thick;
		//base_p1.y = (*itr)->y + (*itr)->ay*base_thick;
		////角度shrink分はここでかける
		//base_p1.z = param->z + base_thick / param->z_shrink;

		//視野中心を原点に移動
		//base_p0 = matrix_3D::addition(base_p0, matrix_3D::const_multiple(center, -1));
		//base_p1 = matrix_3D::addition(base_p1, matrix_3D::const_multiple(center, -1));

		//変換の実行
		base_p0.matrix_multiplication(all_trans);
		base_p0 = matrix_3D::addition(base_p0, shift);
		//base_p1.matrix_multiplication(all_trans);
		//base_p1 = matrix_3D::addition(base_p1, shift);

		//原点をもとに戻す
		//base_p0 = matrix_3D::addition(base_p0, center);
		//base_p1 = matrix_3D::addition(base_p1, center);

		(*itr)->x = base_p0.x;
		(*itr)->y = base_p0.y;
		(*itr)->z = base_p0.z;

		//printf("ax:%.4lf --> %.4lf\n", (*itr)->ax, (base_p0.x - base_p1.x) / (base_p0.z - base_p1.z));
		//printf("ay:%.4lf --> %.4lf\n", (*itr)->ay, (base_p0.y - base_p1.y) / (base_p0.z - base_p1.z));

		//(*itr)->ax = (base_p0.x - base_p1.x) / (base_p0.z - base_p1.z) + param->zx_shear;
		//(*itr)->ay = (base_p0.y - base_p1.y) / (base_p0.z - base_p1.z) + param->zy_shear;

	}
}
void trans_base_all(std::vector < std::pair<Point*, corrmap_3d::align_param2*>>& track_pair) {
	std::map<std::tuple<int, int, int>, corrmap_3d::align_param2*> param_map;
	std::multimap<std::tuple<int, int, int>, Point*>base_map;
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
	std::vector<Point*> t_base;
	for (auto itr = param_map.begin(); itr != param_map.end(); itr++) {
		if (count % 1000 == 0) {
			printf("\r basetrack trans num of triangles %d/%d(%4.1lf%%)", count, param_map.size(), count * 100. / param_map.size());
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
	printf("\r basetrack trans num of triangles %d/%d(%4.1lf%%)\n", count, param_map.size(), count * 100. / param_map.size());

}

