
#define _CRT_SECURE_NO_WARNINGS
//#pragma comment(lib, "VxxReader.lib")
//#include "VxxReader.h"
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.h>
#include <filesystem>
#include <omp.h>

class basetrack_minimum {
public:
	int pl, rawid;
	float ax, ay, x, y;
	Momentum_recon::microtrack_minimum m[2];
};

bool sort_Event_information(Momentum_recon::Event_information& left, Momentum_recon::Event_information& right) {
	if (left.groupid != right.groupid)return left.groupid < right.groupid;
	return left.chains.begin()->chainid < right.chains.begin()->chainid;
}

void change_format(std::vector<mfile0::M_Chain>& chain, std::vector<Momentum_recon::Event_information>& m_c, std::map<int, std::multimap<int, Momentum_recon::Mom_basetrack*>>& read_list);
std::vector< basetrack_minimum> read_base_inf(std::string file_in_base_all, int pl, bool output_flg = false);
void add_basetrack_inf(std::multimap<int, Momentum_recon::Mom_basetrack*>& base_list, std::vector< basetrack_minimum>& base);
void SetPair(std::vector<Momentum_recon::Event_information>& m_c, std::set<std::pair<int, int>>& pl_pair, std::multimap<std::pair<int, int>, Momentum_recon::Mom_basetrack*>& pl_pair_to_base);
void SetBase(std::vector<Momentum_recon::Event_information>& m_c, std::multimap<int, Momentum_recon::Mom_basetrack*>& base_map);

std::vector <std::pair<Momentum_recon::Mom_basetrack*, corrmap_3d::align_param2*>>track_affineparam_correspondence(std::vector<Momentum_recon::Mom_basetrack*>& base, std::vector <corrmap_3d::align_param2>& param);
void trans_base_all(std::vector < std::pair<Momentum_recon::Mom_basetrack*, corrmap_3d::align_param2*>>& track_pair);

//void output_chains(std::string filename, std::vector<Momentum_recon::Mom_chain>&m_c);

std::vector<std::pair<std::pair<int, int>, std::string> >Get_alignment_filename(std::string file_in_ECC);
int use_thread(double ratio, bool output);
void add_direction(std::vector<Momentum_recon::Event_information>& momch);

int main(int argc, char** argv) {
	if (argc != 5) {
		fprintf(stderr, "usage:in-mfile in-base file-in-ECC-path out-mfile\n");
		exit(1);
	}

	std::string file_in_mfile = argv[1];
	std::string file_in_base = argv[2];
	std::string file_in_ECC = argv[3];
	std::string file_out_mfile = argv[4];

	mfile0::Mfile m;
	mfile1::read_mfile_extension(file_in_mfile, m);
	////debug
	//for (auto itr = m.chains.begin(); itr != m.chains.end(); itr++) {
	//	for (auto itr2 = itr->basetracks.begin(); itr2 != itr->basetracks.end(); itr2++) {
	//		if (std::next(itr2, 1) == itr->basetracks.end())continue;
	//		if (abs(itr2->pos / 10 - std::next(itr2, 1)->pos / 10) > 5) {
	//			printf("%d %d\n", itr2->group_id, itr->chain_id);
	//			printf("PL%03d rawid=%d\n", itr2->pos / 10, itr2->rawid);
	//			printf("PL%03d rawid=%d\n", std::next(itr2, 1)->pos / 10, std::next(itr2, 1)->rawid);
	//		}
	//	}
	//}
	//return 0;


	//std::string file_in_corrmap = "K:\\NINJA\\E71a\\work\\suzuki\\Global_ali\\ECC_align\\align\\fine\\corrmap-local-abs.lst";
	std::string file_in_corrmap = file_in_ECC + "\\Area0\\0\\align\\fine\\local\\corrmap-local-abs.lst";

	std::map<int, std::vector<corrmap_3d::align_param>>corrmap = corrmap_3d::read_ali_param_abs(file_in_corrmap, 1);
	std::map<int, std::vector<corrmap_3d::align_param2>>corrmap_dd = corrmap_3d::DelaunayDivide_map(corrmap);

	std::vector<Momentum_recon::Event_information>m_c;
	std::map<int, std::multimap<int, Momentum_recon::Mom_basetrack*>> read_list;

	change_format(m.chains, m_c, read_list);

	for (auto itr = read_list.begin(); itr != read_list.end(); itr++) {
		int pl = itr->first;
		std::vector< basetrack_minimum> base = read_base_inf(file_in_base, pl, true);

		add_basetrack_inf(itr->second, base);

		//for (auto itr2 = itr->second.begin(); itr2 != itr->second.end(); itr2++) {
		//	printf("PL%03d rawid=%d ax=%5.4lf ay=%5.4lf %5d %5d\n", itr->first, itr2->first, itr2->second->ax, itr2->second->ay, itr2->second->m[0].pixelnum, itr2->second->m[1].pixelnum);
		//}
	}
	std::set<std::pair<int, int>> pl_pair;
	std::multimap<std::pair<int, int>, Momentum_recon::Mom_basetrack*>pl_pair_to_base;
	SetPair(m_c, pl_pair, pl_pair_to_base);


	//alignment paramerter read
	//これはポインタで参照されるためとっておく
	std::vector < std::pair<std::pair<int, int>, std::vector <corrmap_3d::align_param >>>all_align_param;
	//std::vector<std::pair< std::pair<int, int>, std::vector <corrmap_3d::align_param2 >>>all_align_param2;
	std::map< std::pair<int, int>, std::vector <corrmap_3d::align_param2 >>all_align_param2;
	std::vector<std::pair<std::pair<int, int>, std::string> >  files_in_align = Get_alignment_filename(file_in_ECC);

	//通常alinment pl0-->pl1 読みこみ
	int count = 0, all = files_in_align.size();
#pragma omp parallel for num_threads(use_thread(0.4,false)) schedule(dynamic,1)
	for (int i = 0; i < files_in_align.size(); i++) {
		if (count % 10 == 0) {
#pragma omp critical
			printf("\r corrmap align read %d/%d", count, all);
		}
#pragma omp atomic
		count++;

		if (pl_pair.count(files_in_align[i].first) == 0)continue;

		std::vector <corrmap_3d::align_param > corr = corrmap_3d::read_ali_param(files_in_align[i].second, false);
#pragma omp critical
		all_align_param.push_back(std::make_pair(files_in_align[i].first, corr));
	}
	printf("\r corrmap align read %d/%d fin\n", count, all);


	//Delaunay3角形への分割
	count = 0, all = all_align_param.size();
#pragma omp parallel for num_threads(use_thread(0.4,false)) schedule(dynamic,1)
	for (int i = 0; i < all_align_param.size(); i++) {
		if (count % 10 == 0) {
#pragma omp critical
			printf("\r corrmap align calc delaunay %d/%d", count, all);
		}
#pragma omp atomic
		count++;

		std::vector <corrmap_3d::align_param2 >corr2 = DelaunayDivide(all_align_param[i].second);
#pragma omp critical
		all_align_param2.insert(std::make_pair(all_align_param[i].first, corr2));

	}
	printf("\r corrmap align calc delaunay %d/%d fin\n", count, all);

	//ここまでは終わっている
	all = pl_pair_to_base.size();
	int cnt = 0, nloop = 0;
	for (auto itr = pl_pair_to_base.begin(); itr != pl_pair_to_base.end(); itr++) {
		if (pl_pair_to_base.count(itr->first) == 0)continue;
		if (all_align_param2.count(itr->first) == 0) {
			fprintf(stderr, "alignment PL%03d PL%03d not found\n", itr->first.first, itr->first.second);
			exit(1);
		}
		printf("base pair trans PL%03d-PL%03d base:%d\n", itr->first.first, itr->first.second, pl_pair_to_base.count(itr->first));

		std::vector <corrmap_3d::align_param2 > ali_target = all_align_param2.at(itr->first);

		std::vector< Momentum_recon::Mom_basetrack*>trans_base;
		auto range = pl_pair_to_base.equal_range(itr->first);
		for (auto res = range.first; res != range.second; res++) {
			trans_base.push_back(res->second);
		}

		//trackとdelaunay3角形の対応
		std::vector <std::pair<Momentum_recon::Mom_basetrack*, corrmap_3d::align_param2*>>track_param = track_affineparam_correspondence(trans_base, ali_target);
		//basetrackを変換
		trans_base_all(track_param);

		itr = std::next(itr, pl_pair_to_base.count(itr->first) - 1);
	}

	std::multimap<int, Momentum_recon::Mom_basetrack*>base_map;
	SetBase(m_c, base_map);
	for (auto itr = base_map.begin(); itr != base_map.end(); itr++) {
		printf("base trans PL%03d\n", itr->first);
		int pl = itr->first;
		int count = base_map.count(pl);
		if (corrmap_dd.count(pl) == 0) {
			fprintf(stderr, "alignment abs local PL%03d not found\n", pl);
			exit(1);
		}
		std::vector <corrmap_3d::align_param2 > ali_target = corrmap_dd.at(pl);

		std::vector<Momentum_recon::Mom_basetrack*>trans_base;
		auto range = base_map.equal_range(pl);
		for (auto res = range.first; res != range.second; res++) {
			trans_base.push_back(res->second);
		}
		//trackとdelaunay3角形の対応
		std::vector <std::pair<Momentum_recon::Mom_basetrack*, corrmap_3d::align_param2*>>track_param = track_affineparam_correspondence(trans_base, ali_target);
		//basetrackを変換
		trans_base_all(track_param);

		itr = std::next(itr, count - 1);
	}
	printf("fin\n");


	//Momentum_recon::Write_mom_chain_txt(file_out_mfile, m_c);
	add_direction(m_c);

	Momentum_recon::Write_Event_information_extension(file_out_mfile, m_c);

	return 0;


}
void change_format(std::vector<mfile0::M_Chain>& chain, std::vector<Momentum_recon::Event_information>& m_c, std::map<int, std::multimap<int, Momentum_recon::Mom_basetrack*>>& read_list) {
	std::multimap<int, Momentum_recon::Event_information>ev_map;
	for (auto& c : chain) {
		Momentum_recon::Event_information ev;
		Momentum_recon::Mom_chain in_c;
		in_c.chainid = c.chain_id;
		ev.groupid = c.basetracks.begin()->group_id;

		for (auto& b : c.basetracks) {
			Momentum_recon::Mom_basetrack in_b;
			in_b.pl = b.pos / 10;
			in_b.rawid = b.rawid;
			in_c.base.push_back(in_b);
		}
		ev.chains.push_back(in_c);
		ev_map.insert(std::make_pair(ev.groupid, ev));
	}

	for (auto itr = ev_map.begin(); itr != ev_map.end(); itr++) {
		auto res = ev_map.find(itr->first);
		Momentum_recon::Event_information event_v = res->second;
		auto range = ev_map.equal_range(itr->first);
		for (auto itr = range.first; itr != range.second; itr++) {
			if (itr == res)continue;
			event_v.chains.push_back(itr->second.chains[0]);
		}
		m_c.push_back(event_v);
		itr = std::next(itr, ev_map.count(itr->first) - 1);
	}

	for (auto& ev : m_c) {
		for (auto& c : ev.chains) {
			for (int i = 0; i < c.base.size(); i++) {
				if (read_list.count(c.base[i].pl) == 0) {
					std::multimap<int, Momentum_recon::Mom_basetrack*>read_list_tmp;
					read_list.insert(std::make_pair(c.base[i].pl, read_list_tmp));
				}
				auto res = read_list.find(c.base[i].pl);
				res->second.insert(std::make_pair(c.base[i].rawid, &c.base[i]));
			}
		}
	}

}

void add_basetrack_inf(std::multimap<int, Momentum_recon::Mom_basetrack*>& base_list, std::vector< basetrack_minimum>& base) {
	int count = 0;
	if (base_list.size() == 0)return;

	for (auto itr = base.begin(); itr != base.end(); itr++) {
		if (base_list.count(itr->rawid) == 0)continue;
		auto range = base_list.equal_range(itr->rawid);
		for (auto res = range.first; res != range.second; res++) {
			res->second->ax = itr->ax;
			res->second->ay = itr->ay;
			res->second->x = itr->x;
			res->second->y = itr->y;
			res->second->z = 0;
			res->second->m[0] = itr->m[0];
			res->second->m[1] = itr->m[1];
			count++;
		}
	}
	printf("basetrack %d : add inf %d\n", base_list.size(), count);
	if (base_list.size() != count) {
		fprintf(stderr, "basetrack not found\n");
		exit(1);
	}
}


std::vector< basetrack_minimum> read_base_inf(std::string file_in_base_all, int pl, bool output_flg) {
	std::ifstream ifs(file_in_base_all, std::ios::binary);
	int read_pl;
	int64_t read_num;

	std::vector< basetrack_minimum> ret;
	while (1) {
		if (ifs.eof())break;
		if (!ifs.read((char*)&read_pl, sizeof(read_pl)))break;
		if (!ifs.read((char*)&read_num, sizeof(read_num)))break;


		if (pl == read_pl) {
			if (output_flg) {
				printf("Read PL%03d tracknum =%lld\n", read_pl, read_num);
			}
			ret.reserve(read_num);
			for (int64_t i = 0; i < read_num; i++) {
				basetrack_minimum b;
				ifs.read((char*)&b, sizeof(basetrack_minimum));
				ret.push_back(b);
			}
			break;
		}
		else {
			ifs.seekg(sizeof(basetrack_minimum) * read_num, std::ios_base::cur);
		}
	}
	ifs.close();

	return ret;
}

void SetPair(std::vector<Momentum_recon::Event_information>& m_c, std::set<std::pair<int, int>>& pl_pair, std::multimap<std::pair<int, int>, Momentum_recon::Mom_basetrack*>& pl_pair_to_base) {
	//for (auto &ev : m_c) {
	//	for (auto &c : ev.chains) {
	for (auto& ev : m_c) {
		for (auto& c : ev.chains) {
			for (int i = 0; i < c.base.size() - 1; i++) {
				c.base_pair.push_back(std::make_pair(c.base[i], c.base[i + 1]));
				pl_pair.insert(std::make_pair(c.base[i].pl, c.base[i + 1].pl));
			}
		}
	}
	for (auto& ev : m_c) {
		for (auto& c : ev.chains) {
			for (auto& p : c.base_pair) {
				pl_pair_to_base.insert(std::make_pair(std::make_pair(p.first.pl, p.second.pl), &(p.second)));
			}
		}
	}
}

void SetBase(std::vector<Momentum_recon::Event_information>& m_c, std::multimap<int, Momentum_recon::Mom_basetrack*>& base_map) {

	for (auto& ev : m_c) {
		for (auto& c : ev.chains) {
			for (auto itr = c.base.begin(); itr != c.base.end(); itr++) {
				base_map.insert(std::make_pair(itr->pl, &(*itr)));
			}
		}
	}
}


std::vector<std::pair<std::pair<int, int>, std::string> >Get_alignment_filename(std::string file_in_ECC) {
	std::vector<std::pair<std::pair<int, int>, std::string> > ret;
	std::map<std::pair<int, int>, std::string> ret_map;

	std::string file_in_align_path = file_in_ECC + "\\Area0\\0\\align\\fine";
	std::filesystem::directory_iterator iter(file_in_align_path), end;
	std::error_code err;

	std::vector<std::string > file_names;
	for (; iter != end && !err; iter.increment(err)) {
		const std::filesystem::directory_entry entry = *iter;

		file_names.push_back(entry.path().string());
		//printf("%s\n", file_names.back().c_str());
	}
	for (int i = 0; i < file_names.size(); i++) {
		int length = file_names[i].size();
		if (file_names[i].substr(length - 4, 4) != ".txt")continue;
		if (file_names[i].substr(length - 18, 14) != "_interpolation")continue;

		int pl0 = std::stoi(file_names[i].substr(length - 25, 3));
		int pl1 = std::stoi(file_names[i].substr(length - 21, 3));

		ret_map.insert(std::make_pair(std::make_pair(pl0, pl1), file_names[i]));
	}
	for (auto itr = ret_map.begin(); itr != ret_map.end(); itr++) {
		ret.push_back(*itr);
	}
	return ret;
}



//basetrack-alignment mapの対応
double select_triangle_vale(corrmap_3d::align_param2* param, Momentum_recon::Mom_basetrack* base) {
	double x, y;
	double dist = 0;
	x = (param->corr_p[0]->x + param->corr_p[1]->x + param->corr_p[2]->x) / 3;
	y = (param->corr_p[0]->y + param->corr_p[1]->y + param->corr_p[2]->y) / 3;
	dist = (base->x - x) * (base->x - x) + (base->y - y) * (base->y - y);
	return dist;
}
corrmap_3d::align_param2* search_param(std::vector<corrmap_3d::align_param*>& param, Momentum_recon::Mom_basetrack* base, std::multimap<int, corrmap_3d::align_param2*>& triangles) {
	//三角形内部
	//最近接三角形
	double dist = 0;
	std::map<double, corrmap_3d::align_param* > dist_map;
	//align_paramを近い順にsort
	for (auto itr = param.begin(); itr != param.end(); itr++) {
		dist = ((*itr)->x - base->x) * ((*itr)->x - base->x) + ((*itr)->y - base->y) * ((*itr)->y - base->y);
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
			sign[0] = (itr2->second->corr_p[1]->x - itr2->second->corr_p[0]->x) * (base->y - itr2->second->corr_p[1]->y) - (itr2->second->corr_p[1]->y - itr2->second->corr_p[0]->y) * (base->x - itr2->second->corr_p[1]->x);
			sign[1] = (itr2->second->corr_p[2]->x - itr2->second->corr_p[1]->x) * (base->y - itr2->second->corr_p[2]->y) - (itr2->second->corr_p[2]->y - itr2->second->corr_p[1]->y) * (base->x - itr2->second->corr_p[2]->x);
			sign[2] = (itr2->second->corr_p[0]->x - itr2->second->corr_p[2]->x) * (base->y - itr2->second->corr_p[0]->y) - (itr2->second->corr_p[0]->y - itr2->second->corr_p[2]->y) * (base->x - itr2->second->corr_p[0]->x);
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
std::vector <std::pair<Momentum_recon::Mom_basetrack*, corrmap_3d::align_param2*>>track_affineparam_correspondence(std::vector<Momentum_recon::Mom_basetrack*>& base, std::vector <corrmap_3d::align_param2>& param) {

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

	std::vector < std::pair<Momentum_recon::Mom_basetrack*, corrmap_3d::align_param2*>> ret;
	std::vector<corrmap_3d::align_param*> param_cand;
	int loop = 0, ix, iy, count = 0;
	for (auto itr = base.begin(); itr != base.end(); itr++) {
		//if (count % 100000 == 0) {
		//	printf("\r search correspond triangles %d/%d(%4.1lf%%)", count, base.size(), count*100. / base.size());
		//}
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
		corrmap_3d::align_param2* param2 = search_param(param_cand, *itr, triangles);
		ret.push_back(std::make_pair((*itr), param2));
	}
	//printf("\r search correspond triangles %d/%d(%4.1lf%%)\n", count, base.size(), count*100. / base.size());

	return ret;
}


//変換 zshrink補正-->9para変換
void trans_base(std::vector<Momentum_recon::Mom_basetrack*>& base, corrmap_3d::align_param2* param) {

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
	matrix_3D::vector_3D base_p0, base_p1;
	double base_thick = 210;
	for (auto itr = base.begin(); itr != base.end(); itr++) {
		base_p0.x = (*itr)->x;
		base_p0.y = (*itr)->y;
		base_p0.z = param->z;

		base_p1.x = (*itr)->x + (*itr)->ax * (base_thick);
		base_p1.y = (*itr)->y + (*itr)->ay * (base_thick);
		//角度shrink分はここでかける
		base_p1.z = param->z + (base_thick) / param->z_shrink;

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

		(*itr)->ax = (base_p0.x - base_p1.x) / (base_p0.z - base_p1.z) + param->zx_shear;
		(*itr)->ay = (base_p0.y - base_p1.y) / (base_p0.z - base_p1.z) + param->zy_shear;
	}
}
void trans_base_all(std::vector < std::pair<Momentum_recon::Mom_basetrack*, corrmap_3d::align_param2*>>& track_pair) {
	std::map<std::tuple<int, int, int>, corrmap_3d::align_param2*> param_map;
	std::multimap<std::tuple<int, int, int>, Momentum_recon::Mom_basetrack*>base_map;
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
	std::vector<Momentum_recon::Mom_basetrack*> t_base;
	for (auto itr = param_map.begin(); itr != param_map.end(); itr++) {
		//if (count % 1000 == 0) {
		//	printf("\r basetrack trans num of triangles %d/%d(%4.1lf%%)", count, param_map.size(), count*100. / param_map.size());
		//}
		count++;

		t_base.clear();

		if (base_map.count(itr->first) == 0)continue;
		auto range = base_map.equal_range(itr->first);
		for (auto res = range.first; res != range.second; res++) {
			t_base.push_back(res->second);
		}
		trans_base(t_base, itr->second);

	}
	//printf("\r basetrack trans num of triangles %d/%d(%4.1lf%%)\n", count, param_map.size(), count*100. / param_map.size());

}

int use_thread(double ratio, bool output) {
	if (ratio >= 1) {
		ratio = 1;
	}
	else if (ratio < 0.1) {
		ratio = 0.1;
	}
	int num_all_thread = omp_get_max_threads();
	if (output) {
		printf("max thread = %d\n", num_all_thread);
		printf("ratio      = %3.2lf\n", ratio);
		printf("using... %d thread\n", int(num_all_thread * ratio));
	}
	return (int)(num_all_thread * ratio);
}

void add_direction(std::vector<Momentum_recon::Event_information>& momch) {

	for (auto& ev : momch) {
		for (auto& c : ev.chains) {
			c.direction = 1;
		}
	}

}