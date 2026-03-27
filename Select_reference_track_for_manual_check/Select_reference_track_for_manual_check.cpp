//コンパイル失敗
#define _CRT_SECURE_NO_WARNINGS
#pragma comment(lib, "VxxReader.lib")
#include "VxxReader.h"
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

#include <filesystem>
#include <set>
//class Point {
//public:
//	double x, y, z;
//};
//class Fiducial_Area {
//public:
//	int pl;
//	Point p[2];
//	//double x0, y0, z0, x1, y1, z1;
//};
//


std::map<int, std::vector<Fiducial_Area::Fiducial_Area>> read_fiducial_Area(std::string filename);
void trans_mfile_cordinate(std::vector<corrmap0::Corrmap>& corr_abs, std::map<int, std::vector<Fiducial_Area::Fiducial_Area>>& area, std::map<int, double>& z_map);
std::vector<std::vector<mfile0::M_Chain>> group_divide(mfile0::Mfile& m);
std::vector<mfile0::M_Chain> divide_single_chain(std::vector<std::vector<mfile0::M_Chain>>& group);
std::map<std::pair<int, int>, mfile0::M_Base> divide_base(std::vector< mfile0::M_Chain>& chain);
bool unique_base_to_chain(std::map<std::pair<int, int>, mfile0::M_Base>& base_all, mfile0::M_Chain& c);
bool merge_chains(std::map<std::pair<int, int>, mfile0::M_Base>& base_all, mfile0::M_Chain& c);
mfile0::M_Chain select_best_chain(std::vector<mfile0::M_Chain>& group);
std::vector<mfile0::M_Chain> divide_penetrate(std::vector<mfile0::M_Chain>& all, std::vector<mfile0::M_Chain>& penetrate, int veto_pl);
std::vector<mfile0::M_Chain> divide_IronECC(std::vector<mfile0::M_Chain>& all, std::vector<mfile0::M_Chain>& iron, int veto_pl);
std::vector<mfile0::M_Chain> divide_edge_out(std::vector<mfile0::M_Chain>& all, std::vector<mfile0::M_Chain>& edge_out, std::map<int, std::vector<Fiducial_Area::Fiducial_Area>>& area, double edge_cut, int ex_pl_max);
bool judge_fiducial_area(std::vector<Fiducial_Area::Fiducial_Area>& area, mfile0::M_Base& b);
void output_mfile(std::string filename, std::vector<mfile0::M_Chain>& chain, mfile0::M_Header& header);
void output_upstream_track(std::ofstream& ofs, std::vector<mfile0::M_Chain>& chain);
void trans_base_all(std::vector < std::pair<Fiducial_Area::Point*, corrmap_3d::align_param2*>>& track_pair);
std::vector <std::pair<Fiducial_Area::Point*, corrmap_3d::align_param2*>>track_affineparam_correspondence(std::vector<Fiducial_Area::Point*>& p, std::vector <corrmap_3d::align_param2>& param);
void trans_mfile_cordinate(std::vector<corrmap_3d::align_param2>& param, std::vector<Fiducial_Area::Fiducial_Area>& area);

bool sort_M_base(const mfile0::M_Base& left, const mfile0::M_Base& right) {
	return left.pos < right.pos;
}

int main(int argc, char** argv) {
	if (argc != 5 && argc != 7) {
		fprintf(stderr, "usage1: prg file-in-mfile(txt) file-in-ECC fa.txt output-path\n");
		fprintf(stderr, "2025/03/10 I add\n");
		fprintf(stderr, "usage2: prg file-in-mfile(txt) file-in-ECC fa.txt output-path thr_pl extraporatePL \n");
		fprintf(stderr, "Default:thr_pl=131, extraporatePL=4\n");
		exit(1);
	}

	std::string file_in_mfile = argv[1];
	std::string file_in_ECC = argv[2];
	std::string file_in_area = argv[3];
	std::string file_out_path = argv[4];
	int thr_pl = 131;
	int extra_pl = 4;

	if (argc == 7) {
		thr_pl = std::stoi(argv[0]);
		extra_pl = std::stoi(argv[1]);
	}

	// set Area

	//corrmap absの読み込み
	std::string file_in_corrmap = file_in_ECC + "\\Area0\\0\\align\\fine\\local\\corrmap-local-abs.lst";
	std::map<int, std::vector<corrmap_3d::align_param>>corrmap = corrmap_3d::read_ali_param_abs(file_in_corrmap, 1);
	std::map<int, std::vector<corrmap_3d::align_param2>>corrmap_dd = corrmap_3d::DelaunayDivide_map(corrmap);

	//fiducial areaの読み込み
	std::map<int, std::vector<Fiducial_Area::Fiducial_Area>> area = read_fiducial_Area(file_in_area);
	for (auto itr = area.begin(); itr != area.end(); itr++) {
		if (corrmap_dd.count(itr->first) == 0) {
			fprintf(stderr, "corrmap local abs PL%03d not found\n", itr->first);
			exit(1);
		}
		std::vector<corrmap_3d::align_param2> param = corrmap_dd.at(itr->first);
		Fiducial_Area::trans_mfile_cordinate(param, itr->second);
	}

	mfile0::Mfile m;
	mfile0::read_mfile(file_in_mfile, m);

	printf("chain size=%d\n", m.chains.size());

	//penetrate check
	std::vector<mfile0::M_Chain> penetrate;
	//最上流がPL131以上
	m.chains = divide_penetrate(m.chains, penetrate, thr_pl);

	//edgeout check
	std::vector<mfile0::M_Chain> edge_out;
	//最上流から4PL外挿,edgeから5mm以内に入ったらedge out
	m.chains = divide_edge_out(m.chains, edge_out, area, 0, extra_pl);

	////Iron ECC event
	//std::vector<mfile0::M_Chain> iron_ecc;
	//single_chain_v = divide_IronECC(single_chain_v, iron_ecc, 15);


	printf("stop = %d\n", m.chains.size());

	std::string file_out_penetrate = file_out_path + "\\muon_penetrate.all";
	std::string file_out_edgeout = file_out_path + "\\muon_edgeout.all";
	std::string file_out_stop = file_out_path + "\\muon_stop.all";
	//std::string file_up_inf = file_out_path + "\\muon_stop.txt";

	std::ofstream ofs_log(file_out_path + "\\classification_log.txt");
	ofs_log << std::right << std::fixed
		<< "penetarte : " << std::setw(5) << std::setprecision(0) << penetrate.size() << " "
		<< "edge out  : " << std::setw(5) << std::setprecision(0) << edge_out.size() << " "
		<< "stop      : " << std::setw(5) << std::setprecision(0) << m.chains.size() << std::endl;
	ofs_log.close();

	output_mfile(file_out_penetrate, penetrate, m.header);
	output_mfile(file_out_edgeout, edge_out, m.header);
	output_mfile(file_out_stop, m.chains, m.header);

}


std::map<int, std::vector<Fiducial_Area::Fiducial_Area>> read_fiducial_Area(std::string filename) {

	std::ifstream ifs(filename);
	std::multimap<int, Fiducial_Area::Fiducial_Area> fa_multi;
	std::map<int, std::vector<Fiducial_Area::Fiducial_Area>> ret;
	Fiducial_Area::Fiducial_Area fa;
	while (ifs >> fa.pl >> fa.p[0].x >> fa.p[0].y >> fa.p[0].z >> fa.p[1].x >> fa.p[1].y >> fa.p[1].z) {
		fa_multi.insert(std::make_pair(fa.pl, fa));
	}

	int count = 0;
	for (auto itr = fa_multi.begin(); itr != fa_multi.end(); itr++) {
		count = fa_multi.count(itr->first);
		auto range = fa_multi.equal_range(itr->first);
		std::vector<Fiducial_Area::Fiducial_Area> fa_vec;
		for (auto itr2 = range.first; itr2 != range.second; itr2++) {
			fa_vec.push_back(itr2->second);
		}
		ret.insert(std::make_pair(itr->first, fa_vec));
	}

	return ret;

}

void trans_mfile_cordinate(std::vector<corrmap_3d::align_param2>& param, std::vector<Fiducial_Area::Fiducial_Area>& area) {

	std::vector< Fiducial_Area::Point*> p_trans;
	for (auto itr = area.begin(); itr != area.end(); itr++) {
		p_trans.push_back(&(itr->p[0]));
		p_trans.push_back(&(itr->p[1]));
	}
	std::vector <std::pair<Fiducial_Area::Point*, corrmap_3d::align_param2*>> p_trans_map = Fiducial_Area::track_affineparam_correspondence(p_trans, param);
	Fiducial_Area::trans_base_all(p_trans_map);
}

std::vector<mfile0::M_Chain> divide_penetrate(std::vector<mfile0::M_Chain>& all, std::vector<mfile0::M_Chain>& penetrate, int veto_pl) {
	std::vector<mfile0::M_Chain> ret;
	int up_pl;
	for (auto itr = all.begin(); itr != all.end(); itr++) {
		up_pl = itr->pos1 / 10;
		if (up_pl >= veto_pl) {
			penetrate.push_back(*itr);
		}
		else {
			ret.push_back(*itr);
		}
	}
	printf("penetrate = %d\n", penetrate.size());
	return ret;
}

std::vector<mfile0::M_Chain> divide_IronECC(std::vector<mfile0::M_Chain>& all, std::vector<mfile0::M_Chain>& iron, int veto_pl) {
	std::vector<mfile0::M_Chain> ret;
	int up_pl;
	for (auto itr = all.begin(); itr != all.end(); itr++) {
		up_pl = itr->pos1 / 10;
		if (up_pl <= veto_pl) {
			iron.push_back(*itr);
		}
		else {
			ret.push_back(*itr);
		}
	}
	printf("IronECC = %d\n", iron.size());
	return ret;
}

std::vector<mfile0::M_Chain> SelectArea(std::vector<mfile0::M_Chain>& all, std::vector<mfile0::M_Chain>& edge_out, std::map<int, std::vector<Fiducial_Area::Fiducial_Area>>& area, double edge_cut, int ex_pl_max) {


	double xmin, xmax, ymin, ymax;


	std::vector<mfile0::M_Chain> ret;
	int up_pl;
	double up_z, up_x, up_y, up_ax, up_ay, ex_z, ex_x, ex_y;
	bool flg = false;
	for (auto itr = all.begin(); itr != all.end(); itr++) {
		up_pl = itr->pos1 / 10;
		flg = false;
		for (int ex_pl = 0; ex_pl <= ex_pl_max; ex_pl++) {
			if (area.count(up_pl + ex_pl) == 0)continue;
			// *** for ecc6 ***
			int plnum = up_pl + ex_pl;
			//if (up_pl < 17 && up_pl + ex_pl >= 16) {
			//	plnum = plnum + 1;
			//}
			// *** ECC6 end ***
			if (!Fiducial_Area::judge_fiducial_area(area.at(plnum), *itr->basetracks.rbegin())) {
				flg = true;
			}
		}
		if (flg) {
			edge_out.push_back(*itr);
		}
		else {
			ret.push_back(*itr);
		}
	}
	printf("edge out = %d\n", edge_out.size());
	return ret;


}
bool judge_area(std::vector<Fiducial_Area::Fiducial_Area>& area, mfile0::M_Base& b) {

	std::map<double, Fiducial_Area::Point> point_map;
	double ex_x, ex_y, dist;
	for (auto itr = area.begin(); itr != area.end(); itr++) {
		ex_x = b.x + b.ax * (itr->p[0].z - b.z);
		ex_y = b.y + b.ay * (itr->p[0].z - b.z);
		dist = pow(ex_x - itr->p[0].x, 2) + pow(ex_y - itr->p[0].y, 2);
		point_map.insert(std::make_pair(dist, itr->p[0]));
	}
	//外挿先から距離の一番近い点のz座標を使用
	double z = point_map.begin()->second.z;
	double x = b.x + b.ax * (z - b.z);
	double y = b.y + b.ay * (z - b.z);


	//true でArea内　falseでarea外

	//点(x,y)からx軸性の方向に直線を引き、その直線と多角形の辺が何回交わるか。
	//下から上に交わったときwn+1
	//上から下に交わったときwn-1
	int wn = 0;
	double vt;
	for (auto itr = area.begin(); itr != area.end(); itr++) {
		// 上向きの辺、下向きの辺によって処理が分かれる。
	// 上向きの辺。点Pがy軸方向について、始点と終点の間にある。ただし、終点は含まない。(ルール1)
		if (itr->p[0].y <= y && itr->p[1].y > y) {
			// 辺は点pよりも右側にある。ただし、重ならない。(ルール4)
			// 辺が点pと同じ高さになる位置を特定し、その時のxの値と点pのxの値を比較する。
			vt = (y - itr->p[0].y) / (itr->p[1].y - itr->p[0].y);
			if (x < itr->p[0].x + vt * (itr->p[1].x - itr->p[0].x)) {
				++wn;  //ここが重要。上向きの辺と交差した場合は+1
			}
		}
		// 下向きの辺。点Pがy軸方向について、始点と終点の間にある。ただし、始点は含まない。(ルール2)
		else if (itr->p[0].y > y && itr->p[1].y <= y) {
			// 辺は点pよりも右側にある。ただし、重ならない。(ルール4)
			// 辺が点pと同じ高さになる位置を特定し、その時のxの値と点pのxの値を比較する。
			vt = (y - itr->p[0].y) / (itr->p[1].y - itr->p[0].y);
			if (x < itr->p[0].x + vt * (itr->p[1].x - itr->p[0].x)) {
				--wn;  //ここが重要。下向きの辺と交差した場合は-1
			}
		}
	}
	if (wn >= 1)return true;
	return false;
}


std::vector<mfile0::M_Chain> divide_edge_out(std::vector<mfile0::M_Chain>& all, std::vector<mfile0::M_Chain>& edge_out, std::map<int, std::vector<Fiducial_Area::Fiducial_Area>>& area, double edge_cut, int ex_pl_max) {


	double xmin, xmax, ymin, ymax;


	std::vector<mfile0::M_Chain> ret;
	int up_pl;
	double up_z, up_x, up_y, up_ax, up_ay, ex_z, ex_x, ex_y;
	bool flg = false;
	for (auto itr = all.begin(); itr != all.end(); itr++) {
		up_pl = itr->pos1 / 10;
		flg = false;
		for (int ex_pl = 0; ex_pl <= ex_pl_max; ex_pl++) {
			if (area.count(up_pl + ex_pl) == 0)continue;
			// *** for ecc6 ***
			int plnum = up_pl + ex_pl;
			//if (up_pl < 17 && up_pl + ex_pl >= 16) {
			//	plnum = plnum + 1;
			//}
			// *** ECC6 end ***
			if (!Fiducial_Area::judge_fiducial_area(area.at(plnum), *itr->basetracks.rbegin())) {
				flg = true;
			}
		}
		if (flg) {
			edge_out.push_back(*itr);
		}
		else {
			ret.push_back(*itr);
		}
	}
	printf("edge out = %d\n", edge_out.size());
	return ret;


}

bool judge_fiducial_area(std::vector<Fiducial_Area::Fiducial_Area>& area, mfile0::M_Base& b) {

	std::map<double, Fiducial_Area::Point> point_map;
	double ex_x, ex_y, dist;
	for (auto itr = area.begin(); itr != area.end(); itr++) {
		ex_x = b.x + b.ax * (itr->p[0].z - b.z);
		ex_y = b.y + b.ay * (itr->p[0].z - b.z);
		dist = pow(ex_x - itr->p[0].x, 2) + pow(ex_y - itr->p[0].y, 2);
		point_map.insert(std::make_pair(dist, itr->p[0]));
	}
	//外挿先から距離の一番近い点のz座標を使用
	double z = point_map.begin()->second.z;
	double x = b.x + b.ax * (z - b.z);
	double y = b.y + b.ay * (z - b.z);


	//true でArea内　falseでarea外

	//点(x,y)からx軸性の方向に直線を引き、その直線と多角形の辺が何回交わるか。
	//下から上に交わったときwn+1
	//上から下に交わったときwn-1
	int wn = 0;
	double vt;
	for (auto itr = area.begin(); itr != area.end(); itr++) {
		// 上向きの辺、下向きの辺によって処理が分かれる。
	// 上向きの辺。点Pがy軸方向について、始点と終点の間にある。ただし、終点は含まない。(ルール1)
		if (itr->p[0].y <= y && itr->p[1].y > y) {
			// 辺は点pよりも右側にある。ただし、重ならない。(ルール4)
			// 辺が点pと同じ高さになる位置を特定し、その時のxの値と点pのxの値を比較する。
			vt = (y - itr->p[0].y) / (itr->p[1].y - itr->p[0].y);
			if (x < itr->p[0].x + vt * (itr->p[1].x - itr->p[0].x)) {
				++wn;  //ここが重要。上向きの辺と交差した場合は+1
			}
		}
		// 下向きの辺。点Pがy軸方向について、始点と終点の間にある。ただし、始点は含まない。(ルール2)
		else if (itr->p[0].y > y && itr->p[1].y <= y) {
			// 辺は点pよりも右側にある。ただし、重ならない。(ルール4)
			// 辺が点pと同じ高さになる位置を特定し、その時のxの値と点pのxの値を比較する。
			vt = (y - itr->p[0].y) / (itr->p[1].y - itr->p[0].y);
			if (x < itr->p[0].x + vt * (itr->p[1].x - itr->p[0].x)) {
				--wn;  //ここが重要。下向きの辺と交差した場合は-1
			}
		}
	}
	if (wn >= 1)return true;
	return false;
}

void output_mfile(std::string filename, std::vector<mfile0::M_Chain>& chain, mfile0::M_Header& header) {

	mfile0::Mfile m;
	m.header = header;
	m.chains = chain;
	mfile0::write_mfile(filename, m);

}

void output_upstream_track(std::ofstream& ofs, std::vector<mfile0::M_Chain>& chain) {
	for (auto itr = chain.begin(); itr != chain.end(); itr++) {
		for (auto itr2 = itr->basetracks.begin(); itr2 != itr->basetracks.end(); itr2++) {
			if (itr2->flg_i[0] != 1)continue;
			ofs << std::right << std::fixed
				<< std::setw(10) << std::setprecision(0) << itr2->group_id << " "
				<< std::setw(5) << std::setprecision(0) << itr2->flg_i[1] << " "
				<< std::setw(5) << std::setprecision(0) << itr2->pos / 10 << " "
				<< std::setw(10) << std::setprecision(0) << itr2->rawid << std::endl;
		}
	}
}

std::vector<std::vector<mfile0::M_Base>> divide_mutli(std::vector<std::vector<mfile0::M_Base>>& all, std::vector<std::vector<mfile0::M_Base>>& single) {
	std::vector<std::vector<mfile0::M_Base>> ret;
	for (auto itr = all.begin(); itr != all.end(); itr++) {
		bool flg = false;
		std::multimap<int, mfile0::M_Base> base_map;
		for (auto itr2 = itr->begin(); itr2 != itr->end(); itr2++) {
			base_map.insert(std::make_pair(itr2->pos / 10, *itr2));
		}
		//最初はmuonなのでmultiになりえない
		mfile0::M_Base prev = base_map.begin()->second;
		int count;
		std::vector<mfile0::M_Base> sel;
		for (auto itr2 = base_map.begin(); itr2 != base_map.end(); itr2++) {
			count = base_map.count(itr2->first);
			if (count == 1) {
				prev = itr2->second;
				sel.push_back(prev);
			}
			//2track以上ある場合
			else {
				//距離を比較。離れていたらmutli-->flgをtrue
				std::vector<mfile0::M_Base> multi_base;
				auto range = base_map.equal_range(itr2->first);
				for (auto res = range.first; res != range.second; res++) {
					multi_base.push_back(res->second);
				}
				double dist;
				for (int i = 0; i < multi_base.size(); i++) {
					for (int j = i + 1; j < multi_base.size(); j++) {
						dist = sqrt(pow(multi_base[i].x - multi_base[j].x, 2) + pow(multi_base[i].x - multi_base[j].x, 2));
						if (dist > 500)flg = true;
					}
				}
				//近接なら良いほうをselにpush back
				if (!flg) {
					double ex_x, ex_y, gap;
					gap = multi_base.begin()->z - prev.z;
					ex_x = prev.x + prev.ax * (gap);
					ex_y = prev.y + prev.ay * (gap);
					for (int i = 0; i < multi_base.size(); i++) {
						if (i == 0) {
							dist = sqrt(pow(multi_base[i].x - ex_x, 2) + pow(multi_base[i].y - ex_y, 2));
							prev = multi_base[i];
						}
						if (dist > sqrt(pow(multi_base[i].x - ex_x, 2) + pow(multi_base[i].y - ex_y, 2))) {
							dist = sqrt(pow(multi_base[i].x - ex_x, 2) + pow(multi_base[i].y - ex_y, 2));
							prev = multi_base[i];
						}
					}
					sel.push_back(prev);
				}
			}

			itr2 = std::next(itr2, count - 1);
		}
		//近接にしかmultiがなかった場合
		if (!flg) {
			single.push_back(sel);
			//single.push_back(*itr);
		}
		else {
			ret.push_back(*itr);
		}
	}
	printf("mutli --> single %d\n", single.size());
	printf("mutli --> mutli  %d\n", ret.size());
	return ret;
}


//mfile0::M_Base
//basetrack-alignment mapの対応
double select_triangle_vale(corrmap_3d::align_param2* param, Fiducial_Area::Point* p) {
	double x, y;
	double dist = 0;
	x = (param->corr_p[0]->x + param->corr_p[1]->x + param->corr_p[2]->x) / 3;
	y = (param->corr_p[0]->y + param->corr_p[1]->y + param->corr_p[2]->y) / 3;
	dist = (p->x - x) * (p->x - x) + (p->y - y) * (p->y - y);
	return dist;
}
corrmap_3d::align_param2* search_param(std::vector<corrmap_3d::align_param*>& param, Fiducial_Area::Point* p, std::multimap<int, corrmap_3d::align_param2*>& triangles) {
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
			if (dist<0 || dist>Fiducial_Area::select_triangle_vale(itr2->second, p)) {
				dist = Fiducial_Area::select_triangle_vale(itr2->second, p);
				ret = itr2->second;
			}
		}
	}
	//printf("Fiducial_Area::Point not in trianlge\n");
	return ret;
}
std::vector <std::pair<Fiducial_Area::Point*, corrmap_3d::align_param2*>>track_affineparam_correspondence(std::vector<Fiducial_Area::Point*>& p, std::vector <corrmap_3d::align_param2>& param) {

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

	std::vector < std::pair<Fiducial_Area::Point*, corrmap_3d::align_param2*>> ret;
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
		}
		corrmap_3d::align_param2* param2 = Fiducial_Area::search_param(param_cand, *itr, triangles);
		ret.push_back(std::make_pair((*itr), param2));
	}
	printf("\r search correspond triangles %d/%d(%4.1lf%%)\n", count, p.size(), count * 100. / p.size());

	return ret;
}
//変換 zshrink補正-->9para変換
void trans_base(std::vector<Fiducial_Area::Point*>& p, corrmap_3d::align_param2* param) {

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
void trans_base_all(std::vector < std::pair<Fiducial_Area::Point*, corrmap_3d::align_param2*>>& track_pair) {
	std::map<std::tuple<int, int, int>, corrmap_3d::align_param2*> param_map;
	std::multimap<std::tuple<int, int, int>, Fiducial_Area::Point*>base_map;
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
	std::vector<Fiducial_Area::Point*> t_base;
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
		Fiducial_Area::trans_base(t_base, itr->second);

	}
	printf("\r basetrack trans num of triangles %d/%d(%4.1lf%%)\n", count, param_map.size(), count * 100. / param_map.size());

}

