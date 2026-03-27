#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>
#include <filesystem>
#include <omp.h>


class Difference {
public:
	int groupid, chainid, vertex_pl, target_pl[2], direction, add_nseg, edge_pl[2], edge_out_flg;
	double ax, ay, dax, day, dx, dy, dar, dal, dr, dl;
	double  s_dar, s_dal, s_dr, s_dl, s_pb_angle, s_pb_position;
	double mcs_pb[2], mcs_pb_error[2][2];
};
class Sigma_list {
public:
	int peke;
	double angle;
	double dal, dar, dr, dl;
};
class basetrack_minimum {
public:
	int pl, rawid;
	float ax, ay, x, y;
	Momentum_recon::microtrack_minimum m[2];
};


std::vector<Difference> Calc_difference(std::vector<Momentum_recon::Event_information>& momch_ev, std::map<int, std::map<int, Momentum_recon::Mom_chain>>& momch_mom_map, std::map<std::pair<int, int>, std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack>>& pair_map, std::map<int, std::vector<Fiducial_Area::Fiducial_Area>>& area);
void Calc_difference(Difference& diff, std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack>& pair);
void output_difference(std::string filename, std::vector<Difference>& diff);
std::map<int, std::map<double, Sigma_list>> Read_sigma_list(std::string filename);
void Calc_sigma(std::vector<Difference>& diff, std::map<int, std::map<double, Sigma_list>>& sigma_list_water, std::map<int, std::map<double, Sigma_list>>& sigma_list_iron);
std::set<std::pair<int, int>> Get_penetrate_list(std::vector<Difference>& diff);
bool judege_Black(Momentum_recon::Mom_chain& chain, std::map<double, std::pair<double, double>>& vph_mip, double thr_sigma);
std::vector<Momentum_recon::Event_information > cut_penetrate(std::vector<Momentum_recon::Event_information >& momch_ev, std::map<int, std::map<int, Momentum_recon::Mom_chain>>& momch_mom_map, std::set<std::pair<int, int>>& id_list, std::string in);
int judge_momch_edgeout(int direction, Momentum_recon::Mom_chain& chain, std::map<int, std::vector<Fiducial_Area::Fiducial_Area>>& area, double edge_cut, int ex_pl_max);
std::map<int, std::map<int, Momentum_recon::Mom_chain>> momch_format_change(std::vector<Momentum_recon::Event_information>& momch_mom);
void Calc_sigma_pb(Difference& diff, int index, std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack>& pair);
bool judge_fiducial_area(std::vector<Fiducial_Area::Fiducial_Area>& area, Momentum_recon::Mom_basetrack& b);
int count_pb_angle_diff(Momentum_recon::Mom_chain& momch);
std::map<std::pair<int, int>, std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack>> peke_pair_pickup(std::vector<Momentum_recon::Event_information>& momch_ev);
std::vector< basetrack_minimum> read_base_inf(std::string file_in_base_all, int pl, bool output_flg);
void add_basetrack_inf(std::multimap<int, Momentum_recon::Mom_basetrack*>& base_list, std::vector< basetrack_minimum>& base);
std::vector <std::pair<Momentum_recon::Mom_basetrack*, corrmap_3d::align_param2*>>track_affineparam_correspondence(std::vector<Momentum_recon::Mom_basetrack*>& base, std::vector <corrmap_3d::align_param2>& param);

void trans_base_all(std::vector < std::pair<Momentum_recon::Mom_basetrack*, corrmap_3d::align_param2*>>& track_pair);
std::vector<std::pair<std::pair<int, int>, std::string> >Get_alignment_filename(std::string file_in_ECC);

int use_thread(double ratio, bool output);
void Set_pair_information(std::map<std::pair<int, int>, std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack>>& pair_map, std::string file_in_ECC, std::string file_in_base);
std::map<double, std::pair<double, double>> vph_mip_distribution(std::string in);
////std::map<double, std::pair<double, double>> vph_mip_distribution();


int main(int argc, char** argv) {
	if (argc != 11) {
		printf("argc %d\n", argc);
		fprintf(stderr, "You need to be argc==11\nusage:file-in-momch-all file-in-momch-mom sigma_list_water sigma_list_iron file-in-ECC fa.txt file_in_base mip_vph_distribution_param file-out file_out_momch\n");
		exit(1);
	}
	std::string file_in_momch_all = argv[1];
	std::string file_in_momch_mom = argv[2];
	std::string file_in_sigma_water = argv[3];
	std::string file_in_sigma_iron = argv[4];
	std::string file_in_ECC = argv[5];
	std::string file_in_area = argv[6];
	std::string file_in_base = argv[7];

	std::string in = argv[8];
	std::string file_out_difference = argv[9];
	std::string file_out_momch = argv[10];

	//corrmap absの読み込み
	std::string file_in_corrmap = file_in_ECC + "\\Area0\\0\\align\\fine\\local\\corrmap-local-abs.lst";
	std::map<int, std::vector<corrmap_3d::align_param>>corrmap = corrmap_3d::read_ali_param_abs(file_in_corrmap, 1);
	std::map<int, std::vector<corrmap_3d::align_param2>>corrmap_dd = corrmap_3d::DelaunayDivide_map(corrmap);

	//fiducial areaの読み込み
	std::map<int, std::vector<Fiducial_Area::Fiducial_Area>> area = Fiducial_Area::read_fiducial_Area(file_in_area);
	for (auto itr = area.begin(); itr != area.end(); itr++) {
		if (corrmap_dd.count(itr->first) == 0) {
			fprintf(stderr, "corrmap local abs PL%03d not found\n", itr->first);
			exit(1);
		}
		std::vector<corrmap_3d::align_param2> param = corrmap_dd.at(itr->first);
		trans_mfile_cordinate(param, itr->second);
	}


	std::map<int, std::map<double, Sigma_list>> sigma_list_water = Read_sigma_list(file_in_sigma_water);
	std::map<int, std::map<double, Sigma_list>> sigma_list_iron = Read_sigma_list(file_in_sigma_iron);

	std::vector<Momentum_recon::Event_information> momch_all = Momentum_recon::Read_Event_information_extension(file_in_momch_all);
	std::vector<Momentum_recon::Event_information> momch_mom = Momentum_recon::Read_Event_information_extension(file_in_momch_mom);

	std::map<int, std::map<int, Momentum_recon::Mom_chain>> momch_mom_map = momch_format_change(momch_mom);
	std::map<std::pair<int, int>, std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack>> peke_pair = peke_pair_pickup(momch_all);
	Set_pair_information(peke_pair, file_in_ECC, file_in_base);

	//peke pairを使った実装
	std::vector<Difference>diff = Calc_difference(momch_all, momch_mom_map, peke_pair, area);
	Calc_sigma(diff, sigma_list_water, sigma_list_iron);
	output_difference(file_out_difference, diff);
	std::set<std::pair<int, int>> penetrate_id = Get_penetrate_list(diff);
	//peke pairを使った実装
	momch_all = cut_penetrate(momch_all, momch_mom_map, penetrate_id, in);
	Momentum_recon::Write_Event_information_extension(file_out_momch, momch_all);
}


std::map<int, std::map<int, Momentum_recon::Mom_chain>> momch_format_change(std::vector<Momentum_recon::Event_information>& momch_mom) {

	std::map<int, std::map<int, Momentum_recon::Mom_chain>> ret;
	for (auto& ev : momch_mom) {
		std::map<int, Momentum_recon::Mom_chain> chain_event;
		for (auto& ch : ev.chains) {
			chain_event.insert(std::make_pair(ch.chainid, ch));
		}
		ret.insert(std::make_pair(ev.groupid, chain_event));
	}
	return ret;


}


std::map<int, std::map<double, Sigma_list>> Read_sigma_list(std::string filename) {

	std::map<int, std::map<double, Sigma_list>> ret;
	std::ifstream ifs(filename);
	int peke, num;
	while (ifs >> peke >> num) {
		std::map<double, Sigma_list> sigma_map;
		Sigma_list sig;
		sig.peke = peke;
		double angle_min, angle_max;
		for (int i = 0; i < num; i++) {
			ifs >> sig.angle >> sig.dal >> sig.dar >> sig.dl >> sig.dr;
			sigma_map.insert(std::make_pair(sig.angle, sig));
		}
		ret.insert(std::make_pair(sig.peke, sigma_map));
	}
	return ret;

	//for (auto itr = ret.begin(); itr != ret.end(); itr++) {
	//	printf("peke:%d\n", itr->first);
	//	for (auto s : itr->second) {
	//		printf("%.1lf %.5lf %.5lf %.2lf %.2lf\n", s.second.angle, s.second.dal, s.second.dar, s.second.dl, s.second.dr);
	//	}
	//}

}
std::map<std::pair<int, int>, std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack>> peke_pair_pickup(std::vector<Momentum_recon::Event_information>& momch_ev) {
	std::map<std::pair<int, int>, std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack>> ret;
	std::pair<int, int> id;
	for (auto& ev : momch_ev) {
		id.first = ev.groupid;
		for (auto& c : ev.chains) {
			id.second = c.chainid;
			if (c.chainid == 0)continue;
			std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack> divide_pair;
			//forward
			int edge_pl = 0;
			int hit_count = 0;
			if (c.direction == 1) {
				//attachしているPLの決定
				for (auto& b : c.base) {
					if (b.pl <= ev.vertex_pl) {
						edge_pl = b.pl;
					}
				}
				for (auto& b : c.base) {
					if (b.pl == edge_pl) {
						break;
					}
					divide_pair.first = b;
					hit_count = 1;
				}
				for (auto& b : c.base) {
					if (b.pl > ev.vertex_pl) {
						divide_pair.second = b;
						hit_count++;
						break;
					}
				}
			}
			else if (c.direction == -1) {
				//attachしているPLの決定
				for (auto& b : c.base) {
					if (b.pl > ev.vertex_pl) {
						edge_pl = b.pl;
					}
				}
				for (auto& b : c.base) {
					if (b.pl <= ev.vertex_pl) {
						divide_pair.first = b;
						hit_count = 1;
					}
				}

				for (auto& b : c.base) {
					if (b.pl > edge_pl) {
						divide_pair.second = b;
						hit_count++;
						break;
					}
				}
			}

			if (hit_count == 2) {
				ret.insert(std::make_pair(id, divide_pair));
			}
		}
	}
	return ret;
}
void Set_pair_information(std::map<std::pair<int, int>, std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack>>& pair_map, std::string file_in_ECC, std::string file_in_base) {
	std::map<int, std::multimap<int, Momentum_recon::Mom_basetrack*>> read_list;
	std::set<std::pair<int, int>> pl_pair;
	std::multimap<std::pair<int, int>, Momentum_recon::Mom_basetrack*>pl_pair_to_base;

	for (auto itr = pair_map.begin(); itr != pair_map.end(); itr++) {
		pl_pair.insert(std::make_pair(itr->second.first.pl, itr->second.second.pl));
		pl_pair_to_base.insert(std::make_pair(std::make_pair(itr->second.first.pl, itr->second.second.pl), &itr->second.second));


		auto res = read_list.find(itr->second.first.pl);
		if (res == read_list.end()) {
			std::multimap<int, Momentum_recon::Mom_basetrack*> base_map_p;
			base_map_p.insert(std::make_pair(itr->second.first.rawid, &(itr->second.first)));
			read_list.insert(std::make_pair(itr->second.first.pl, base_map_p));
		}
		else {
			res->second.insert(std::make_pair(itr->second.first.rawid, &(itr->second.first)));
		}

		res = read_list.find(itr->second.second.pl);
		if (res == read_list.end()) {
			std::multimap<int, Momentum_recon::Mom_basetrack*> base_map_p;
			base_map_p.insert(std::make_pair(itr->second.second.rawid, &(itr->second.second)));
			read_list.insert(std::make_pair(itr->second.second.pl, base_map_p));
		}
		else {
			res->second.insert(std::make_pair(itr->second.second.rawid, &(itr->second.second)));
		}
	}
	for (auto itr = read_list.begin(); itr != read_list.end(); itr++) {
		int pl = itr->first;
		std::vector< basetrack_minimum> base = read_base_inf(file_in_base, pl, true);

		add_basetrack_inf(itr->second, base);

	}

	//for (auto &pl_p : pl_pair) {
	//	printf("PL%03d - PL%03d\n",pl_p.first,pl_p.second);
	//}

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

	std::set<std::pair<int, int>>del_pl_pair;

	for (auto itr = pl_pair_to_base.begin(); itr != pl_pair_to_base.end(); itr++) {
		if (pl_pair_to_base.count(itr->first) == 0)continue;
		if (all_align_param2.count(itr->first) == 0) {
			del_pl_pair.insert(itr->first);
			continue;
		}
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

	for (auto itr = pair_map.begin(); itr != pair_map.end();) {
		if (del_pl_pair.count(itr->first) == 1) {
			itr = pair_map.erase(itr);
		}
		else {
			itr++;
		}
	}

	for (auto& pair : del_pl_pair) {
		printf("del PL pair PL%03d PL%03d\n", pair.first, pair.second);
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
	//exit(1); //2026/01/20 commentout
	}
}



std::vector<Difference> Calc_difference(std::vector<Momentum_recon::Event_information>& momch_ev, std::map<int, std::map<int, Momentum_recon::Mom_chain>>& momch_mom_map, std::map<std::pair<int, int>, std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack>>& pair_map, std::map<int, std::vector<Fiducial_Area::Fiducial_Area>>& area) {
	std::vector<Difference> ret;

	for (auto& ev : momch_ev) {
		Difference diff;
		diff.groupid = ev.groupid;
		std::map<int, Momentum_recon::Mom_chain> mom_map, mom_map_inv;
		if (momch_mom_map.count(diff.groupid) == 0) {
			fprintf(stderr, "event %d partner not found\n", diff.groupid);
			continue;
		}
		if (momch_mom_map.count(diff.groupid * -1) == 0) {
			fprintf(stderr, "event %d partner inv not found\n", diff.groupid);
			continue;
		}
		mom_map = momch_mom_map.at(diff.groupid);
		mom_map_inv = momch_mom_map.at(-1 * diff.groupid);
		for (auto& c : ev.chains) {
			if (c.chainid == 0)continue;
			auto pair_itr = pair_map.find(std::make_pair(diff.groupid, c.chainid));
			if (pair_itr == pair_map.end())continue;
			diff.vertex_pl = ev.vertex_pl;
			diff.direction = c.direction;
			diff.chainid = c.chainid;
			diff.edge_pl[0] = c.base.begin()->pl;
			diff.edge_pl[1] = c.base.rbegin()->pl;
			diff.edge_out_flg = judge_momch_edgeout(diff.direction, c, area, 0, 4);
			diff.add_nseg = 0;
			diff.target_pl[0] = -1;
			diff.target_pl[1] = -1;
			diff.s_pb_angle = 0;
			diff.s_pb_position = 0;
			diff.mcs_pb[0] = -1;
			diff.mcs_pb[1] = -1;
			diff.mcs_pb_error[0][0] = -1;
			diff.mcs_pb_error[0][1] = -1;
			diff.mcs_pb_error[1][0] = -1;
			diff.mcs_pb_error[1][1] = -1;
			std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack> divide_pair;
			//forward
			if (diff.direction == 1) {
				diff.target_pl[0] = diff.vertex_pl;

				for (auto& b : c.base) {
					if (b.pl > diff.vertex_pl) {
						if (diff.target_pl[1] < 0) {
							diff.target_pl[1] = b.pl;
						}
						diff.add_nseg += 1;
					}
				}
				if (diff.target_pl[1] < 0)continue;
				for (auto& pair : c.base_pair) {
					if (pair.first.pl == diff.target_pl[0] && pair.second.pl == diff.target_pl[1]) {
						Calc_difference(diff, pair);
						divide_pair = pair;
						break;
					}
				}

			}
			//backward
			else if (diff.direction == -1) {
				for (auto b = c.base.rbegin(); b != c.base.rend(); b++) {
					if (b->pl <= diff.vertex_pl) {
						if (diff.target_pl[1] < 0) {
							diff.target_pl[0] = b->pl;
							diff.target_pl[1] = std::next(b, -1)->pl;
						}
						diff.add_nseg += 1;
					}
				}
				if (diff.target_pl[1] < 0)continue;
				for (auto& pair : c.base_pair) {
					if (pair.first.pl == diff.target_pl[0] && pair.second.pl == diff.target_pl[1]) {
						Calc_difference(diff, pair);
						divide_pair = pair;
						break;
					}
				}

			}
			diff.target_pl[0] = pair_itr->second.first.pl;
			diff.target_pl[1] = pair_itr->second.second.pl;
			Calc_difference(diff, pair_itr->second);

			double pb_error[2];
			Momentum_recon::Mom_chain mom, mom_inv;
			if (mom_map.count(c.chainid) == 0) {
				fprintf(stderr, "event %d chain %d partner not found\n", diff.groupid, c.chainid);
			}
			else {
				mom = mom_map.at(c.chainid);
				diff.mcs_pb[0] = mom.Get_muon_mcs_pb();
				mom.Get_muon_pb_mcs_error(pb_error);
				diff.mcs_pb_error[0][0] = pb_error[0];
				diff.mcs_pb_error[0][1] = pb_error[1];
			}

			if (mom_map_inv.count(c.chainid) == 0) {
				fprintf(stderr, "event %d chain %d partner inv not found\n", diff.groupid, c.chainid);
				continue;
			}
			else {
				mom_inv = mom_map_inv.at(c.chainid);
				diff.mcs_pb[1] = mom_inv.Get_muon_mcs_pb();
				mom_inv.Get_muon_pb_mcs_error(pb_error);
				diff.mcs_pb_error[1][0] = pb_error[0];
				diff.mcs_pb_error[1][1] = pb_error[1];
			}
			if (count_pb_angle_diff(mom) < 2 && count_pb_angle_diff(mom_inv) < 2) {
				Calc_sigma_pb(diff, 2, divide_pair);
			}
			else if (count_pb_angle_diff(mom) < count_pb_angle_diff(mom_inv)) {
				Calc_sigma_pb(diff, 1, divide_pair);
			}
			else {
				Calc_sigma_pb(diff, 0, pair_itr->second);
			}
			ret.push_back(diff);
		}
	}
	return ret;

}
int count_pb_angle_diff(Momentum_recon::Mom_chain& momch) {
	int ret = 0;
	for (auto& pair : momch.base_pair) {
		if (pair.second.pl - pair.first.pl != 1)continue;
		if (pair.first.pl <= 3)continue;
		else if (pair.first.pl <= 14)ret++;
		else if (pair.first.pl <= 15)continue;
		else if (pair.first.pl % 2 == 0)ret++;

	}
	return ret;
}
void Calc_sigma_pb(Difference& diff, int index, std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack>& pair) {
	double pb;
	if (index >= 2) {
		diff.s_pb_angle = 0;
		diff.s_pb_position = 0;
		return;
	}
	if (diff.mcs_pb[index] > 0) {
		pb = diff.mcs_pb[index];
	}
	else {
		diff.s_pb_angle = 0;
		diff.s_pb_position = 0;
		return;
	}
	double path = sqrt(pow(pair.first.x - pair.second.x, 2) + pow(pair.first.y - pair.second.y, 2) + pow(pair.first.z - pair.second.z, 2));
	double dz = fabs(pair.first.z - pair.second.z);
	int num_iron = 0, num_film = 0, num_water = 0;
	for (int pl = pair.first.pl; pl < pair.second.pl; pl++) {
		num_film++;
		if (pl == 3)continue;
		else if (pl <= 14)num_iron++;
		else if (pl == 15)continue;
		else if (pl % 2 == 0)num_iron++;
		else num_water++;
	}
	double thick_iron_ratio, thick_water_ratio, thick_base_ratio, thick_gel_ratio, thick_pack_ratio;
	thick_iron_ratio = 500 * num_iron / dz;
	thick_base_ratio = 210 * num_film / dz;
	thick_gel_ratio = 70 * num_film * 2 / dz;
	thick_pack_ratio = 100 * num_water * 2 / dz;
	thick_water_ratio = 1 - (thick_iron_ratio + thick_gel_ratio + thick_base_ratio + thick_pack_ratio);
	double radiation_length = 0;
	radiation_length += thick_iron_ratio * path / (17.57 * 1000);
	radiation_length += thick_gel_ratio * path / (30.3 * 1000);
	radiation_length += thick_base_ratio * path / (413.1 * 1000);
	radiation_length += thick_pack_ratio * path / (413.1 * 1000);
	radiation_length += thick_water_ratio * path / (360.8 * 1000);
	double mass = 105.65836668;
	double p2 = pb * pb / 2 * (1 + sqrt(1 + pow(2 * mass / pb, 2)));
	double beta = 1 / sqrt(1 + pow(mass / sqrt(p2), 2));

	diff.s_pb_angle = 13.6 / pb * sqrt(radiation_length) * (1 + 0.038 * log(radiation_length / beta));
	diff.s_pb_position = diff.s_pb_angle * path / (4 * sqrt(3));
}
int judge_pl_material(int pl, int peke) {
	//1:iron
	//2:water
	int ret = -1;
	if (pl <= 15 || pl % 2 == 0)ret = 1;
	else ret = 2;

	return ret;
}
void Calc_sigma(std::vector<Difference>& diff, std::map<int, std::map<double, Sigma_list>>& sigma_list_water, std::map<int, std::map<double, Sigma_list>>& sigma_list_iron) {

	int peke;
	double angle, sigma_dal, sigma_dar, sigma_dr, sigma_dl;
	Sigma_list sigma[2];
	for (auto itr = diff.begin(); itr != diff.end(); itr++) {
		peke = itr->target_pl[1] - itr->target_pl[0] - 1;
		angle = sqrt(pow(itr->ax, 2) + pow(itr->ay, 2));
		auto sigma_peke = sigma_list_iron.find(peke);

		if (judge_pl_material(itr->target_pl[0], peke) == 1) {
			auto sigma_peke = sigma_list_iron.find(peke);
			if (sigma_peke == sigma_list_iron.end()) {
				fprintf(stderr, "PL%03d - PL%03d\n\n", itr->target_pl[0], itr->target_pl[1]);
				fprintf(stderr, "peke=%d not found\n", peke);
				//itr->s_dal = 0;
				//itr->s_dar = 0;
				//itr->s_dl = 0;
				//itr->s_dr = 0;
				continue;
			}
		}
		else if (judge_pl_material(itr->target_pl[0], peke) == 2) {
			auto sigma_peke = sigma_list_water.find(peke);
			if (sigma_peke == sigma_list_water.end()) {
				fprintf(stderr, "PL%03d - PL%03d\n\n", itr->target_pl[0], itr->target_pl[1]);
				fprintf(stderr, "peke=%d not found\n", peke);
				//itr->s_dal = 0;
				//itr->s_dar = 0;
				//itr->s_dl = 0;
				//itr->s_dr = 0;
				continue;
			}
		}
		else continue;

		auto sigma_angle = sigma_peke->second.upper_bound(angle);
		if (sigma_angle == sigma_peke->second.begin()) {
			sigma[0] = sigma_angle->second;
			sigma[1] = sigma_angle->second;
			sigma_dal = sigma[0].dal;
			sigma_dar = sigma[0].dar;
			sigma_dl = sigma[0].dl;
			sigma_dr = sigma[0].dr;
		}
		else if (sigma_angle == sigma_peke->second.end()) {
			sigma[0] = sigma_peke->second.rbegin()->second;
			sigma[1] = sigma_peke->second.rbegin()->second;
			sigma_dal = sigma[0].dal;
			sigma_dar = sigma[0].dar;
			sigma_dl = sigma[0].dl;
			sigma_dr = sigma[0].dr;
		}
		else {
			sigma[0] = std::next(sigma_angle, -1)->second;
			sigma[1] = sigma_angle->second;
			sigma_dal = (sigma[1].dal - sigma[0].dal) / (sigma[1].angle - sigma[0].angle) * (angle - sigma[0].angle) + sigma[0].dal;
			sigma_dar = (sigma[1].dar - sigma[0].dar) / (sigma[1].angle - sigma[0].angle) * (angle - sigma[0].angle) + sigma[0].dar;
			sigma_dl = (sigma[1].dl - sigma[0].dl) / (sigma[1].angle - sigma[0].angle) * (angle - sigma[0].angle) + sigma[0].dl;
			sigma_dr = (sigma[1].dr - sigma[0].dr) / (sigma[1].angle - sigma[0].angle) * (angle - sigma[0].angle) + sigma[0].dr;

		}

		sigma_dal = sqrt(pow(sigma_dal, 2) + pow(itr->s_pb_angle, 2));
		sigma_dar = sqrt(pow(sigma_dar, 2) + pow(itr->s_pb_angle, 2));
		sigma_dl = sqrt(pow(sigma_dl, 2) + pow(itr->s_pb_position, 2));
		sigma_dr = sqrt(pow(sigma_dr, 2) + pow(itr->s_pb_position, 2));

		itr->s_dal = itr->dal / sigma_dal;
		itr->s_dar = itr->dar / sigma_dar;
		itr->s_dl = itr->dl / sigma_dl;
		itr->s_dr = itr->dr / sigma_dr;
	}


}
void Calc_difference(Difference& diff, std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack>& pair) {
	diff.ax = pair.first.ax;
	diff.ay = pair.first.ay;
	diff.dax = pair.second.ax - pair.first.ax;
	diff.day = pair.second.ay - pair.first.ay;

	diff.dx = pair.second.x - pair.first.x - (pair.first.ax + pair.second.ax) / 2 * (pair.second.z - pair.first.z);
	diff.dy = pair.second.y - pair.first.y - (pair.first.ay + pair.second.ay) / 2 * (pair.second.z - pair.first.z);
	diff.dar = (diff.dax * pair.first.ax + diff.day * pair.first.ay) / sqrt(pair.first.ax * pair.first.ax + pair.first.ay * pair.first.ay);
	diff.dal = (diff.dax * pair.first.ay - diff.day * pair.first.ax) / sqrt(pair.first.ax * pair.first.ax + pair.first.ay * pair.first.ay);

	double denominator, d_theta_r, d_theta_l;
	denominator = (pair.first.ax * pair.second.ax + pair.first.ay * pair.second.ay + 1) / (sqrt(pair.first.ax * pair.first.ax + pair.first.ay * pair.first.ay + 1));
	d_theta_r = (-1 * pair.first.ax * pair.second.ax - pair.first.ay * pair.second.ay + pair.first.ax * pair.first.ax + pair.first.ay * pair.first.ay)
		/
		(sqrt(pair.first.ax * pair.first.ax + pair.first.ay * pair.first.ay + pow(pair.first.ax * pair.first.ax + pair.first.ay * pair.first.ay, 2)));
	d_theta_l = (-1 * pair.first.ay * pair.second.ax + pair.first.ax * pair.second.ay) / sqrt(pair.first.ax * pair.first.ax + pair.first.ay * pair.first.ay);

	diff.dar = atan(d_theta_r / denominator);
	diff.dal = atan(d_theta_l / denominator);


	//dr,dlの計算
	using namespace matrix_3D;
	vector_3D pos0, pos1, dir0, dir1;
	pos0.x = pair.first.x;
	pos0.y = pair.first.y;
	pos0.z = pair.first.z;
	dir0.x = pair.first.ax;
	dir0.y = pair.first.ay;
	dir0.z = 1;
	pos1.x = pair.second.x;
	pos1.y = pair.second.y;
	pos1.z = pair.second.z;
	dir1.x = pair.second.ax;
	dir1.y = pair.second.ay;
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

	diff.dr = dot(addition(extra1, const_multiple(extra0, -1)), unit_r);
	diff.dl = dot(addition(extra1, const_multiple(extra0, -1)), unit_l);

	if (sqrt(pair.first.ax * pair.first.ax + pair.first.ay * pair.first.ay) < 0.001) {
		diff.dar = diff.day;
		diff.dal = diff.dax;
		diff.dr = diff.dy;
		diff.dl = diff.dx;
	}

}
void output_difference(std::string filename, std::vector<Difference>& diff) {

	std::ofstream ofs(filename);
	for (auto itr = diff.begin(); itr != diff.end(); itr++) {
		ofs << std::right << std::fixed
			<< std::setw(8) << std::setprecision(0) << itr->groupid << " "
			<< std::setw(8) << std::setprecision(0) << itr->chainid << " "
			<< std::setw(2) << std::setprecision(0) << itr->direction << " "
			<< std::setw(4) << std::setprecision(0) << itr->vertex_pl << " "
			<< std::setw(2) << std::setprecision(0) << itr->edge_out_flg << " "
			<< std::setw(4) << std::setprecision(0) << itr->target_pl[0] << " "
			<< std::setw(4) << std::setprecision(0) << itr->target_pl[1] << " "
			<< std::setw(4) << std::setprecision(0) << itr->edge_pl[0] << " "
			<< std::setw(4) << std::setprecision(0) << itr->edge_pl[1] << " "
			<< std::setw(4) << std::setprecision(0) << itr->add_nseg << " "
			<< std::setw(7) << std::setprecision(4) << itr->ax << " "
			<< std::setw(7) << std::setprecision(4) << itr->ay << " "
			<< std::setw(7) << std::setprecision(4) << itr->dax << " "
			<< std::setw(7) << std::setprecision(4) << itr->day << " "
			<< std::setw(8) << std::setprecision(1) << itr->dx << " "
			<< std::setw(8) << std::setprecision(1) << itr->dy << " "
			<< std::setw(7) << std::setprecision(4) << itr->dar << " "
			<< std::setw(7) << std::setprecision(4) << itr->dal << " "
			<< std::setw(8) << std::setprecision(1) << itr->dr << " "
			<< std::setw(8) << std::setprecision(1) << itr->dl << " "
			<< std::setw(7) << std::setprecision(4) << itr->s_dar << " "
			<< std::setw(7) << std::setprecision(4) << itr->s_dal << " "
			<< std::setw(7) << std::setprecision(4) << itr->s_dr << " "
			<< std::setw(7) << std::setprecision(4) << itr->s_dl << " "
			<< std::setw(7) << std::setprecision(4) << itr->s_pb_angle << " "
			<< std::setw(7) << std::setprecision(1) << itr->s_pb_position << " "
			<< std::setw(7) << std::setprecision(1) << itr->mcs_pb[0] << " "
			<< std::setw(7) << std::setprecision(1) << itr->mcs_pb_error[0][0] << " "
			<< std::setw(7) << std::setprecision(1) << itr->mcs_pb_error[0][1] << " "
			<< std::setw(7) << std::setprecision(1) << itr->mcs_pb[1] << " "
			<< std::setw(7) << std::setprecision(1) << itr->mcs_pb_error[1][0] << " "
			<< std::setw(7) << std::setprecision(1) << itr->mcs_pb_error[1][1] << std::endl;
	}
}

std::set<std::pair<int, int>> Get_penetrate_list(std::vector<Difference>& diff) {
	std::set<std::pair<int, int>> ret;
	for (auto itr = diff.begin(); itr != diff.end(); itr++) {
		double chi2 = pow(itr->s_dal, 2) + pow(itr->s_dar, 2) + pow(itr->s_dl, 2) + pow(itr->s_dr, 2);
		if (chi2 > 30)continue;

		ret.insert(std::make_pair(itr->groupid, itr->chainid));
	}
	return ret;
}

std::vector<Momentum_recon::Event_information > cut_penetrate(std::vector<Momentum_recon::Event_information >& momch_ev, std::map<int, std::map<int, Momentum_recon::Mom_chain>>& momch_mom_map, std::set<std::pair<int, int>>& id_list, std::string in) {
	std::vector<Momentum_recon::Event_information > ret;
	std::map<double, std::pair<double, double>> vph_mip = vph_mip_distribution(in);
	//int count[4] = {};

	//mfile0::Mfile m;
	//mfile1::read_mfile_extension("K:\\NINJA\\E71a\\work\\suzuki\\muon_analysis\\03_vertex_location\\re_connection\\group_penetrate_check\\m_out1.all", m);
	//std::set<std::pair<int, int>> target_chain;
	//for (auto &c : m.chains) {
	//	int gid = c.basetracks.begin()->group_id/100000;
	//	int cid = (c.basetracks.begin()->group_id%100000)/10;
	//	//150             13100451    268556    270098 - 0.0811   0.1893    61641.5    63264.6 - 12486  0  0  0  0 0.0000 0.0000

	//	target_chain.insert(std::make_pair(gid, cid));
	//}




	for (auto& ev : momch_ev) {
		Momentum_recon::Event_information chains = ev;
		std::map<int, Momentum_recon::Mom_chain> mom_map;
		if (momch_mom_map.count(ev.groupid) == 0) {
			fprintf(stderr, "event %d partner not found\n", ev.groupid);
		}
		mom_map = momch_mom_map.at(ev.groupid);

		//chains = ev;
		chains.chains.clear();
		int direction, vertex_pl;

		for (auto& c : ev.chains) {
			Momentum_recon::Mom_chain mom;
			if (c.chainid == 0) {
				if (mom_map.count(c.chainid) == 0) {
					fprintf(stderr, "event %d chain %d partner not found\n", ev.groupid, c.chainid);
				}
				else {
					mom = mom_map.at(c.chainid);
				}
				chains.chains.push_back(mom);
				continue;
			}
			vertex_pl = ev.vertex_pl;
			direction = c.direction;
			if (mom_map.count(c.chainid) == 0) {
				fprintf(stderr, "event %d chain %d partner not found\n", ev.groupid, c.chainid);
			}
			else {
				mom = mom_map.at(c.chainid);
			}

			Momentum_recon::Mom_chain chain = c;
			chain.base.clear();
			chain.base_pair.clear();
			//chainの分断
			//forward
			if (direction == 1) {
				for (auto& b : c.base) {
					if (b.pl <= vertex_pl) {
						chain.base.push_back(b);
					}
				}
				for (auto& pair : c.base_pair) {
					if (pair.second.pl <= vertex_pl) {
						chain.base_pair.push_back(pair);
					}
				}
			}
			else if (direction == -1) {
				for (auto& b : c.base) {
					if (b.pl > vertex_pl) {
						chain.base.push_back(b);
					}
				}
				for (auto& pair : c.base_pair) {
					if (pair.first.pl > vertex_pl) {
						chain.base_pair.push_back(pair);
					}
				}
			}

			//penetrate candidate
			if (id_list.count(std::make_pair(ev.groupid, chain.chainid)) == 1) {
				//blackの判定
				if (!judege_Black(chain, vph_mip, 5))continue;
			}

			//if (target_chain.count(std::make_pair(ev.groupid, chain.chainid)) == 1) {

			//	//if (c.base.begin()->pl <= vertex_pl && vertex_pl < c.base.rbegin()->pl) {
			//	if (id_list.count(std::make_pair(ev.groupid, chain.chainid)) == 1 && judege_Black(chain, vph_mip, 5))count[0]++;
			//	else if (id_list.count(std::make_pair(ev.groupid, chain.chainid)) != 1 && judege_Black(chain, vph_mip, 5))count[1]++;
			//	else if (id_list.count(std::make_pair(ev.groupid, chain.chainid)) == 1 && !judege_Black(chain, vph_mip, 5))count[2]++;
			//	else if (id_list.count(std::make_pair(ev.groupid, chain.chainid)) != 1 && !judege_Black(chain, vph_mip, 5))count[3]++;
			//	else {
			//		printf("exception\n");
			//	}
			//}

			//chains.chains.push_back(chain);
			chains.chains.push_back(mom);

		}
		ret.push_back(chains);

	}
	//printf(" penetrate  black %d\n", count[0]);
	//printf("!penetrate  black %d\n", count[1]);
	//printf(" penetrate !black %d\n", count[2]);
	//printf("!penetrate !black %d\n", count[3]);

	return ret;
}
bool judege_Black(Momentum_recon::Mom_chain& chain, std::map<double, std::pair<double, double>>& vph_mip, double thr_sigma) {
	double ax = 0, ay = 0, vph = 0;
	int count = 0, count2 = 0;
	for (auto itr = chain.base.begin(); itr != chain.base.end(); itr++) {
		ax += itr->ax;
		ay += itr->ay;
		count++;
		vph += itr->m[0].ph % 10000;
		count2++;
		vph += itr->m[1].ph % 10000;
		count2++;
	}
	ax /= count;
	ay /= count;
	vph /= count2;

	double angle = sqrt(ax * ax + ay * ay);

	double sigma = 0;
	auto vph_angle = vph_mip.upper_bound(angle);
	if (vph_angle == vph_mip.begin()) {
		sigma = (vph - vph_angle->second.first) / vph_angle->second.second;
	}
	else if (vph_angle == vph_mip.end()) {
		sigma = (vph - vph_mip.rbegin()->second.first) / vph_mip.rbegin()->second.second;
	}
	else {
		double vph_mean, vph_sigma;
		auto val0 = std::next(vph_angle, -1);
		auto val1 = vph_angle;
		vph_mean = (val1->second.first - val0->second.first) / (val1->first - val0->first) * (angle - val0->first) + val0->second.first;
		vph_sigma = (val1->second.second - val0->second.second) / (val1->first - val0->first) * (angle - val0->first) + val0->second.second;
		sigma = (vph - vph_mean) / vph_sigma;
	}
	if (sigma > thr_sigma)return true;
	return false;
}
//std::map<double, std::pair<double, double>> vph_mip_distribution() {
//	std::map<double, std::pair<double, double>> ret;
//	ret.insert(std::make_pair(0.05, std::make_pair(61.86, 6.8)));
//	ret.insert(std::make_pair(0.15, std::make_pair(49.68, 7.6)));
//	ret.insert(std::make_pair(0.25, std::make_pair(36.38, 5.0)));
//	ret.insert(std::make_pair(0.35, std::make_pair(30.85, 4.0)));
//	ret.insert(std::make_pair(0.45, std::make_pair(29.20, 3.8)));
//	ret.insert(std::make_pair(0.55, std::make_pair(27.17, 3.9)));
//	ret.insert(std::make_pair(0.65, std::make_pair(25.14, 3.8)));
//	ret.insert(std::make_pair(0.80, std::make_pair(22.26, 3.6)));
//	ret.insert(std::make_pair(1.00, std::make_pair(19.83, 3.3)));
//	ret.insert(std::make_pair(1.20, std::make_pair(18.39, 3.2)));
//	ret.insert(std::make_pair(1.40, std::make_pair(17.24, 3.1)));
//	ret.insert(std::make_pair(1.70, std::make_pair(16.02, 3.0)));
//	ret.insert(std::make_pair(2.10, std::make_pair(15.29, 3.0)));
//	ret.insert(std::make_pair(2.50, std::make_pair(14.19, 3.0)));
//	return ret;
//}
std::map<double, std::pair<double, double>> vph_mip_distribution(std::string in) {
	std::map<double, std::pair<double, double>> ret;

	std::ifstream ifs(in);
	if (!ifs) {
		std::cerr << "Failed to open " << in << std::endl;
		exit(0);
	}

	double p[3];
	while (ifs >> p[0] >> p[1] >> p[2]) {
		ret.insert(std::make_pair(p[0], std::make_pair(p[1], p[2])));
	}

	for (auto itr = ret.begin(); itr != ret.end(); itr++) {
		std::cout << itr->first << " " << itr->second.first << " " << itr->second.second << std::endl;
	}

	return ret;
}
int judge_momch_edgeout(int direction, Momentum_recon::Mom_chain& chain, std::map<int, std::vector<Fiducial_Area::Fiducial_Area>>& area, double edge_cut, int ex_pl_max) {

	int edge_pl;
	int return_flg = 0;
	if (direction == 1) {
		edge_pl = chain.base.begin()->pl;
		return_flg = 0;
		for (int ex_pl = 0; ex_pl <= ex_pl_max; ex_pl++) {
			if (area.count(edge_pl - ex_pl) == 0)continue;
			//ex_z = z_map.at(up_pl + ex_pl);
			//ex_x = up_x + up_ax * (ex_z - up_z);
			//ex_y = up_y + up_ay * (ex_z - up_z);

			if (!judge_fiducial_area(area.at(edge_pl - ex_pl), *chain.base.begin())) {
				return_flg = 2;
			}
		}
		if (edge_pl <= 4)return_flg = 1;
	}
	else if (direction == -1) {
		edge_pl = chain.base.rbegin()->pl;
		return_flg = 0;
		for (int ex_pl = 0; ex_pl <= ex_pl_max; ex_pl++) {
			if (area.count(edge_pl + ex_pl) == 0)continue;
			//ex_z = z_map.at(up_pl + ex_pl);
			//ex_x = up_x + up_ax * (ex_z - up_z);
			//ex_y = up_y + up_ay * (ex_z - up_z);

			if (!judge_fiducial_area(area.at(edge_pl + ex_pl), *chain.base.rbegin())) {
				return_flg = 2;
			}
		}
		if (edge_pl >= 132)return_flg = 1;

	}

	return return_flg;

}

bool judge_fiducial_area(std::vector<Fiducial_Area::Fiducial_Area>& area, Momentum_recon::Mom_basetrack& b) {

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
