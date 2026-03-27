// 2025/03/10
// bag修正版+FILEstructureのDelauneyDivideを修正したlibに差し替え

#define _CRT_SECURE_NO_WARNINGS
#pragma comment(lib, "VxxReader.lib")
#include "VxxReader.h"
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>
#pragma comment(lib,"functions.lib")
#include <functions.hpp>

#include <filesystem>
#include <set>
#include <omp.h>

class Track_file {
public:
	int eventid, trackid, pl, rawid;
};
bool operator<(const Track_file& left, const Track_file& right) {
	if (left.eventid != right.eventid)return left.eventid < right.eventid;
	else if (left.trackid != right.trackid)return left.trackid < right.trackid;
	else if (left.pl != right.pl)return left.pl < right.pl;
	else  return left.rawid < right.rawid;
}

class basetrack_minimum {
public:
	int pl, rawid, ph;
	float ax, ay, x, y;
};

class basetrack_minimum_z :public basetrack_minimum {
public:
	int trkid, black_flg;
	float z, ex_z0, ex_z1, md;
};

class output_track {
public:
	int eventid;
	basetrack_minimum_z muon;
	std::vector<basetrack_minimum_z> partner;
};

std::map<int, corrmap0::Corrmap> read_corrmap_abs(std::string file_in_ECC);
std::map<int, std::vector< basetrack_minimum_z>> read_base(std::string filename);
std::vector<std::pair<std::pair<int, int>, std::string> >Get_alignment_filename(std::string file_in_ECC);
std::vector<std::pair<std::pair<int, int>, std::string> >Get_alignment_inv_filename(std::string file_in_ECC);

vxx::base_track_t read_muon(std::string file_in_ECC, int pl, mfile0::M_Base mu_base);

std::vector <std::pair<basetrack_minimum_z*, corrmap_3d::align_param2*>>track_affineparam_correspondence(std::vector<basetrack_minimum_z>& base, std::vector <corrmap_3d::align_param2>& param);
void trans_base_all(std::vector < std::pair<basetrack_minimum_z*, corrmap_3d::align_param2*>>& track_pair);

void partner_search(vxx::base_track_t target, std::vector<basetrack_minimum_z>& all, double gap, double z_range[2], double accuracy, std::vector<basetrack_minimum_z>& connect);
void Set_search_z(double z_range[2], int muonPL);


std::pair<bool, vxx::base_track_t> basetrack_apply_local_ali(vxx::base_track_t t, std::string file_in_ECC, int pl0, int pl1);
vxx::base_track_t basetrack_tans(vxx::base_track_t t, std::vector<corrmap0::Corrmap> corr);
void basetrack_inverse_trans(vxx::base_track_t& t, corrmap0::Corrmap param);
void basetrack_affine_trans(vxx::base_track_t& t, corrmap0::Corrmap param);
void corrmap_area_inverse(std::vector<corrmap0::Corrmap>& corr);
void partner_search(vxx::base_track_t target, std::vector<basetrack_minimum_z>& all, double z_range[2], double angle_accuracy_intercept_mu, double angle_accuracy_slope_mu, double angle_accuracy_intercept, double angle_accuracy_slope, std::vector<basetrack_minimum_z>& connect, std::ofstream& ofs, int black_flg);
void pre_partner_search(vxx::base_track_t target, std::vector<basetrack_minimum_z>& all, double z_range[2], double angle_accuracy_intercept_mu, double angle_accuracy_slope_mu, double angle_accuracy_intercept, double angle_accuracy_slope, std::vector<std::pair<basetrack_minimum_z, basetrack_minimum_z>>& attach_candidate, int black_flg);
void partner_search_single(vxx::base_track_t target, basetrack_minimum_z& all, double z_range[2], double angle_accuracy_intercept_mu, double angle_accuracy_slope_mu, double angle_accuracy_intercept, double angle_accuracy_slope, std::vector<basetrack_minimum_z>& connect, std::ofstream& ofs, int black_flg);

std::vector<mfile0::M_Chain> Get_Mfile_chain(std::vector<vxx::base_track_t>& base, std::map<int, double> z_map, std::vector<corrmap0::Corrmap>& corr, mfile0::M_Chain mu);
void output_Track(std::string filename, std::vector<Track_file>& track);
void output_Track_inf(std::string filename, std::vector<output_track>& track);
std::vector<basetrack_minimum_z> target_base_selection(std::vector<basetrack_minimum_z>& base, vxx::base_track_t muon_base, std::map<int, corrmap0::Corrmap>& corrmap_abs, int stop_pl, int target_pl);
std::vector<basetrack_minimum_z> base_minimum_multidel(std::vector<basetrack_minimum_z>& b);
std::vector <std::pair< basetrack_minimum_z, double>> search_local_dz(std::vector<std::pair<basetrack_minimum_z, basetrack_minimum_z>>attach_candidate, std::vector <  basetrack_minimum_z>& local_dz);

int use_thread(double ratio, bool output);

int main(int argc, char** argv) {
	if (argc != 8) {
		fprintf(stderr, "usage:prg in-mu-mfile in-base_link in-base_black in-ECC-Area-path output-track-id output-track-file output_md_file\n");
		exit(1);
	}


	std::string file_in_mfile = argv[1];
	std::string file_in_base_link = argv[2];
	std::string file_in_base_black = argv[3];
	std::string file_in_ECC_path = argv[4];
	std::string file_out_track_id = argv[5];
	std::string file_out_track = argv[6];
	std::string file_out_md = argv[7];

	//gap nominal read
	std::stringstream structure_path;
	structure_path << file_in_ECC_path << "\\..\\st\\st.dat";
	chamber1::Chamber chamber;
	chamber1::read_structure(structure_path.str(), chamber);
	std::map<int, double> z_map = chamber1::base_z_convert(chamber);

	//corrmap abs read
	std::map<int, corrmap0::Corrmap> corrmap_abs = read_corrmap_abs(file_in_ECC_path);


	//alignment paramerter read
	//����̓|�C���^�ŎQ�Ƃ���邽�߂Ƃ��Ă���
	std::vector < std::pair<std::pair<int, int>, std::vector <corrmap_3d::align_param >>>all_align_param;
	//std::vector<std::pair< std::pair<int, int>, std::vector <corrmap_3d::align_param2 >>>all_align_param2;
	std::map< std::pair<int, int>, std::vector <corrmap_3d::align_param2 >>all_align_param2;
	std::vector<std::pair<std::pair<int, int>, std::string> >  files_in_align = Get_alignment_filename(file_in_ECC_path);
	std::vector<std::pair<std::pair<int, int>, std::string> >  files_in_align_inv = Get_alignment_inv_filename(file_in_ECC_path);

	//�ʏ�alinment pl0-->pl1 �ǂ݂���
	int count = 0, all = files_in_align.size();
#pragma omp parallel for num_threads(use_thread(0.4,false)) schedule(dynamic,1)
	for (int i = 0; i < files_in_align.size(); i++) {
		if (count % 10 == 0) {
#pragma omp critical
			printf("\r corrmap align read %d/%d", count, all);
		}
#pragma omp atomic
		count++;
		//�א�+1peke�����ǂ�
		if (files_in_align[i].first.second - files_in_align[i].first.first > 2)continue;

		std::vector <corrmap_3d::align_param > corr = corrmap_3d::read_ali_param(files_in_align[i].second, false);
#pragma omp critical
		all_align_param.push_back(std::make_pair(files_in_align[i].first, corr));
	}
	printf("\r corrmap align read %d/%d fin\n", count, all);

	//alinment inv pl1-->pl0 �ǂ݂���
	count = 0, all = files_in_align_inv.size();
#pragma omp parallel for num_threads(use_thread(0.4,false)) schedule(dynamic,1)
	for (int i = 0; i < files_in_align_inv.size(); i++) {
		if (count % 10 == 0) {
#pragma omp critical
			printf("\r corrmap align inv read %d/%d", count, all);
		}
#pragma omp atomic
		count++;
		//�אڂ����ǂ�
		if (files_in_align_inv[i].first.first - files_in_align_inv[i].first.second > 1)continue;

		std::vector <corrmap_3d::align_param > corr = corrmap_3d::read_ali_param(files_in_align_inv[i].second, false);
#pragma omp critical
		all_align_param.push_back(std::make_pair(files_in_align_inv[i].first, corr));
	}
	printf("\r corrmap align inv read %d/%d fin\n", count, all);

	//Delaunay3�p�`�ւ̕���
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

	//muon read 
	mfile0::Mfile muon;
	mfile0::read_mfile(file_in_mfile, muon);

	//basetrack read
	std::map<int, std::vector< basetrack_minimum_z>> base_link = read_base(file_in_base_link);
	std::map<int, std::vector< basetrack_minimum_z>> base_black = read_base(file_in_base_black);
	for (auto itr = base_black.begin(); itr != base_black.end(); itr++) {
		for (auto& b : itr->second) {
			b.black_flg = 1;
		}
	}
	/*
	std::map<std::pair<int, int>, std::vector< basetrack_minimum_z>>base_trans_all;
	count = 0, all = all_align_param2.size();
#pragma omp parallel for num_threads(use_thread(0.4,false)) schedule(dynamic,1)
	for (int i = 0; i < all_align_param2.size(); i++) {
		if (count % 10 == 0) {
#pragma omp critical
			printf("\r basetrack trans %d/%d", count, all);
		}
#pragma omp atomic
		count++;

		if (base_all.count(all_align_param2[i].first.second) == 0) {
			fprintf(stderr, "PL%03d basetrack not found\n", all_align_param2[i].first.second);
			exit(1);
		}
		std::vector<basetrack_minimum_z> base_target = base_all.at(all_align_param2[i].first.second);

		std::vector <corrmap_3d::align_param2 > ali_target = all_align_param2[i].second;
		//track��delaunay3�p�`�̑Ή�
		std::vector <std::pair<basetrack_minimum_z*, corrmap_3d::align_param2*>>track_param = track_affineparam_correspondence(base_target, ali_target);
		//basetrack��ϊ�
		trans_base_all(track_param);
#pragma omp critical
		base_trans_all.insert(std::make_pair(all_align_param2[i].first, base_target));
	}
	printf("\r basetrack trans %d/%d fin\n", count, all);
	*/

	//1�{���T�����B
	std::set<Track_file> partner_cand;
	std::vector<output_track> output_partner_cand;
	std::ofstream ofs(file_out_md);

	for (int i = 0; i < muon.chains.size(); i++) {
		printf("muon %4d/%4d start\n", i + 1, muon.chains.size());
		mfile0::M_Base up_muon = *(muon.chains[i].basetracks.rbegin());
		int stop_pl = up_muon.pos / 10;
		double mu_z = z_map.at(stop_pl);
		//basetrack�ǂݍ���
		//PL,rawid����muon file�̓ǂݍ���
		vxx::base_track_t muon_base = read_muon(file_in_ECC_path, stop_pl, up_muon);
		double z_range[2];
		//z_range�̐ݒ������
		Set_search_z(z_range, stop_pl);
		std::vector<basetrack_minimum_z> local_dz;
		{
			std::pair<int, int> pl_pair = std::make_pair(stop_pl, stop_pl + 1);
			if (all_align_param2.count(pl_pair) == 0) {
				fprintf(stderr, "corrmap align PL%03d - PL%03d not fnoud\n", pl_pair.first, pl_pair.second);
				continue;
			}
			std::vector <corrmap_3d::align_param2 > ali_target = all_align_param2.at(pl_pair);
			int ali_id = 0;
			for (auto& ali : ali_target) {
				basetrack_minimum_z point;
				point.x = ali.x;
				point.y = ali.y;
				point.z = 0;
				point.ax = 0;
				point.ay = 0;
				point.black_flg = 0;
				point.ex_z0 = 0;
				point.ex_z1 = 0;
				point.md = 0;
				point.ph = 0;
				point.rawid = ali_id;
				point.trkid = 0;
				point.pl = stop_pl + 1;
				ali_id++;
				local_dz.push_back(point);
			}
			//track��delaunay3�p�`�̑Ή�
			std::vector <std::pair<basetrack_minimum_z*, corrmap_3d::align_param2*>>track_param_dz = track_affineparam_correspondence(local_dz, ali_target);
			//basetrack��ϊ�
			trans_base_all(track_param_dz);

			//for (auto &dz : local_dz) {
			//	printf("PL%03d PL%03d dz=%.1lf\n", stop_pl, stop_pl + 1, dz.z);
			//}
		}

		//����PL�ł�attach�̒T��
		//nominal gap + corrmap�ŕϊ�
		//md�ŒT��
		std::vector<basetrack_minimum_z> connect;
		int search_plus = 1, search_minus = 0;
		for (int search_pl = stop_pl - search_minus; search_pl <= stop_pl + search_plus; search_pl++) {
			if (search_pl < 3 || search_pl>133)continue;
			printf("muon PL%03d search PL%03d\n", stop_pl, search_pl);

			std::vector<basetrack_minimum_z> target_base_link, target_base_black;
			if (base_link.count(search_pl) == 0) {
				fprintf(stderr, "PL%03d basetrack link not found\n", search_pl);
				exit(1);
			}
			if (base_black.count(search_pl) == 0) {
				fprintf(stderr, "PL%03d basetrack black not found\n", search_pl);
				exit(1);
			}
			target_base_link = base_link.at(search_pl);
			target_base_black = base_black.at(search_pl);

			//muon���ӗ̈�̐؂���
			target_base_link = target_base_selection(target_base_link, muon_base, corrmap_abs, stop_pl, search_pl);
			target_base_black = target_base_selection(target_base_black, muon_base, corrmap_abs, stop_pl, search_pl);

			//�ϊ�
			if (search_pl != stop_pl) {
				std::pair<int, int> pl_pair = std::make_pair(stop_pl, search_pl);
				if (all_align_param2.count(pl_pair) == 0) {
					fprintf(stderr, "corrmap align PL%03d - PL%03d not fnoud\n", pl_pair.first, pl_pair.second);
					continue;
				}
				std::vector <corrmap_3d::align_param2 > ali_target = all_align_param2.at(pl_pair);
				//track��delaunay3�p�`�̑Ή�
				std::vector <std::pair<basetrack_minimum_z*, corrmap_3d::align_param2*>>track_param_link = track_affineparam_correspondence(target_base_link, ali_target);
				std::vector <std::pair<basetrack_minimum_z*, corrmap_3d::align_param2*>>track_param_black = track_affineparam_correspondence(target_base_black, ali_target);
				//basetrack��ϊ�
				trans_base_all(track_param_link);
				trans_base_all(track_param_black);
			}

			//dz�͂��΂��΂�md search
			std::vector<std::pair<basetrack_minimum_z, basetrack_minimum_z>>attach_candidate_all, attach_candidate_black;
			pre_partner_search(muon_base, target_base_link, z_range, 0.04, 0.04, 0.04, 0.04, attach_candidate_all, 0);
			pre_partner_search(muon_base, target_base_black, z_range, 0.04, 0.04, 0.04, 0.04, attach_candidate_black, 1);

			//local_dz����md�ɋ߂�dz���擾-->�g�p
			//���̓�͑Ή�
			std::vector <std::pair< basetrack_minimum_z, double>>target_base_link_sel, target_base_black_sel;
			target_base_link_sel = search_local_dz(attach_candidate_all, local_dz);
			target_base_black_sel = search_local_dz(attach_candidate_black, local_dz);

			for (int j = 0; j < target_base_link_sel.size(); j++) {
				//dz�̐ݒ�(gap+emulsion+a)
				z_range[0] = target_base_link_sel[j].second - 70 - 30;
				partner_search_single(muon_base, target_base_link_sel[j].first, z_range, 0.04, 0.04, 0.04, 0.04, connect, ofs, 0);
			}
			for (int j = 0; j < target_base_black_sel.size(); j++) {
				//dz�̐ݒ�
				z_range[0] = target_base_black_sel[j].second - 70 - 30;
				partner_search_single(muon_base, target_base_black_sel[j].first, z_range, 0.04, 0.04, 0.04, 0.04, connect, ofs, 1);
			}
		}

		connect = base_minimum_multidel(connect);

		output_track out_tmp;
		out_tmp.eventid = muon.chains[i].basetracks.begin()->group_id;
		out_tmp.muon.ax = muon_base.ax;
		out_tmp.muon.ay = muon_base.ay;
		out_tmp.muon.x = muon_base.x;
		out_tmp.muon.y = muon_base.y;
		out_tmp.muon.z = 0;
		out_tmp.muon.pl = muon_base.pl;
		out_tmp.muon.rawid = muon_base.rawid;
		out_tmp.muon.ph = muon_base.m[0].ph + muon_base.m[1].ph;
		out_tmp.muon.trkid = 1;
		out_tmp.muon.ex_z0 = z_range[0];
		out_tmp.muon.ex_z1 = z_range[1];
		out_tmp.muon.md = -1;


		for (int j = 0; j < connect.size(); j++) {
			Track_file t;
			t.eventid = muon.chains[i].basetracks.begin()->group_id;
			t.trackid = j + 1;
			t.pl = connect[j].pl;
			t.rawid = connect[j].rawid;
			partner_cand.insert(t);
			out_tmp.partner.push_back(connect[j]);
		}
		output_partner_cand.push_back(out_tmp);

	}

	std::vector<Track_file>track_out(partner_cand.begin(), partner_cand.end());
	output_Track(file_out_track_id, track_out);
	output_Track_inf(file_out_track, output_partner_cand);

}
std::vector <std::pair< basetrack_minimum_z, double>> search_local_dz(std::vector<std::pair<basetrack_minimum_z, basetrack_minimum_z>>attach_candidate, std::vector <  basetrack_minimum_z>& local_dz) {
	std::vector <std::pair< basetrack_minimum_z, double>> ret;
	std::pair< basetrack_minimum_z, double> t_dz_pair;
	for (auto& t : attach_candidate) {

		double dist = DBL_MAX;
		double dist_tmp = -1;
		double dz_min;
		for (auto& dz : local_dz) {
			dist_tmp = pow(dz.x - t.second.x, 2) + pow(dz.y - t.second.y, 2);
			if (dist_tmp < dist) {
				dist = dist_tmp;
				dz_min = dz.z;
			}
		}
		//printf("distans %.1lf dz=%.1lf\n", sqrt(dist), dz_min);
		t_dz_pair.first = t.first;
		t_dz_pair.second = dz_min;
		ret.push_back(t_dz_pair);
	}
	return ret;
}

std::map<int, corrmap0::Corrmap> read_corrmap_abs(std::string file_in_ECC) {
	std::stringstream file_in_corr;
	file_in_corr << file_in_ECC << "\\0\\align\\corrmap-abs.lst";

	std::vector<corrmap0::Corrmap> corr;
	corrmap0::read_cormap(file_in_corr.str(), corr);

	std::map<int, corrmap0::Corrmap> ret;
	for (auto itr = corr.begin(); itr != corr.end(); itr++) {
		ret.insert(std::make_pair(itr->pos[0] / 10, *itr));
	}
	return ret;
}

int64_t Basetrack_num(std::string filename) {
	std::ifstream ifs(filename, std::ios::binary);
	//filesize�擾
	ifs.seekg(0, std::ios::end);
	int64_t eofpos = ifs.tellg();
	ifs.clear();
	ifs.seekg(0, std::ios::beg);
	int64_t begpos = ifs.tellg();
	int64_t nowpos = ifs.tellg();
	int64_t size2 = eofpos - begpos;
	return size2 / sizeof(basetrack_minimum);
}
std::map<int, std::vector< basetrack_minimum_z>> read_base(std::string filename) {
	const int PL_MAX = 133;
	std::vector< basetrack_minimum>  base_buf[PL_MAX];
	int64_t track_num = Basetrack_num(filename);
	for (int i = 0; i < PL_MAX; i++) {
		base_buf[i].reserve(track_num / PL_MAX * 2);
	}

	std::ifstream ifs(filename, std::ios::binary);
	//filesize�擾
	ifs.seekg(0, std::ios::end);
	int64_t eofpos = ifs.tellg();
	ifs.clear();
	ifs.seekg(0, std::ios::beg);
	int64_t begpos = ifs.tellg();
	int64_t nowpos = ifs.tellg();
	int64_t size2 = eofpos - begpos;
	int64_t GB = size2 / (1000 * 1000 * 1000);
	int64_t MB = (size2 - GB * 1000 * 1000 * 1000) / (1000 * 1000);
	int64_t KB = (size2 - GB * 1000 * 1000 * 1000 - MB * 1000 * 1000) / (1000);
	if (GB > 0) {
		std::cout << "FILE size :" << GB << "." << MB << " [GB]" << std::endl;
	}
	else {
		std::cout << "FILE size :" << MB << "." << KB << " [MB]" << std::endl;
	}
	int64_t count = 0;
	basetrack_minimum b;
	while (ifs.read((char*)&b, sizeof(basetrack_minimum))) {
		if (count % 10000 == 0) {
			nowpos = ifs.tellg();
			auto size1 = nowpos - begpos;
			std::cerr << std::right << std::fixed << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%";
		}
		count++;
		//����
		if (b.pl<1 || b.pl>PL_MAX) {
			fprintf(stderr, "PL exception\n");
			fprintf(stderr, "requrie 1<=PL<=P%d\n", PL_MAX);
			fprintf(stderr, "PL=%d exist\n", b.pl);
			exit(1);
		}
		base_buf[b.pl - 1].emplace_back(b);
	}
	auto size1 = eofpos - begpos;
	std::cerr << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%" << std::endl;;
	if (count == 0) {
		fprintf(stderr, "%s no basetrack!\n", filename.c_str());
		exit(1);
	}

	std::map<int, std::vector< basetrack_minimum_z>> ret;
	for (int i = 0; i < PL_MAX; i++) {
		int pl = i + 1;
		if (base_buf[i].size() == 0)continue;
		//base_buf[i].shrink_to_fit();
		std::vector< basetrack_minimum_z> base_z;
		base_z.reserve(base_buf[i].size());
		for (auto itr = base_buf[i].begin(); itr != base_buf[i].end(); itr++) {
			basetrack_minimum_z b_z;
			b_z.ax = itr->ax;
			b_z.ay = itr->ay;
			b_z.x = itr->x;
			b_z.y = itr->y;
			b_z.z = 0;
			b_z.ph = itr->ph;
			b_z.pl = itr->pl;
			b_z.rawid = itr->rawid;
			b_z.trkid = 0;
			b_z.black_flg = 0;
			b_z.ex_z0 = 0;
			b_z.ex_z1 = 0;

			base_z.push_back(b_z);
		}
		ret.insert(std::make_pair(pl, base_z));
	}
	return ret;

}

std::vector<std::pair<std::pair<int, int>, std::string> >Get_alignment_filename(std::string file_in_ECC) {
	std::vector<std::pair<std::pair<int, int>, std::string> > ret;
	std::map<std::pair<int, int>, std::string> ret_map;

	std::string file_in_align_path = file_in_ECC + "\\0\\align\\fine";
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

std::vector<std::pair<std::pair<int, int>, std::string> >Get_alignment_inv_filename(std::string file_in_ECC) {
	std::vector<std::pair<std::pair<int, int>, std::string> > ret;
	std::map<std::pair<int, int>, std::string> ret_map;

	std::string file_in_align_path = file_in_ECC + "\\0\\align_inv\\fine";
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

vxx::base_track_t read_muon(std::string file_in_ECC, int pl, mfile0::M_Base mu_base) {

	std::stringstream file_in_base;
	file_in_base << file_in_ECC << "\\PL" << std::setw(3) << std::setfill('0') << pl
		<< "\\b" << std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";

	vxx::BvxxReader br;
	std::array<int, 2> index = { mu_base.rawid,mu_base.rawid + 1 };//1234<=rawid<=5678�ł���悤�Ȃ��̂�����ǂށB

	std::vector<vxx::base_track_t >base = br.ReadAll(file_in_base.str(), pl, 0, vxx::opt::index = index);

	vxx::base_track_t ret;
	bool flg = false;
	for (auto itr = base.begin(); itr != base.end(); itr++) {
		if (itr->rawid == mu_base.rawid) {
			ret = *itr;
			flg = true;
			break;
		}
	}

	if (!flg) {
		printf("rawid=%d basetrack not found\n", mu_base.rawid);
		return ret;
	}
	return ret;
}



//basetrack-alignment map�̑Ή�
double select_triangle_vale(corrmap_3d::align_param2* param, basetrack_minimum_z& base) {
	double x, y;
	double dist = 0;
	x = (param->corr_p[0]->x + param->corr_p[1]->x + param->corr_p[2]->x) / 3;
	y = (param->corr_p[0]->y + param->corr_p[1]->y + param->corr_p[2]->y) / 3;
	dist = (base.x - x) * (base.x - x) + (base.y - y) * (base.y - y);
	return dist;
}
corrmap_3d::align_param2* search_param(std::vector<corrmap_3d::align_param*>& param, basetrack_minimum_z& base, std::multimap<int, corrmap_3d::align_param2*>& triangles) {
	//�O�p�`����
	//�ŋߐڎO�p�`
	double dist = 0;
	std::map<double, corrmap_3d::align_param* > dist_map;
	//align_param���߂�����sort
	for (auto itr = param.begin(); itr != param.end(); itr++) {
		dist = ((*itr)->x - base.x) * ((*itr)->x - base.x) + ((*itr)->y - base.y) * ((*itr)->y - base.y);
		dist_map.insert(std::make_pair(dist, (*itr)));
	}

	double sign[3];
	bool flg = false;
	int id;

	corrmap_3d::align_param2* ret = triangles.begin()->second;
	for (auto itr = dist_map.begin(); itr != dist_map.end(); itr++) {
		if (itr != dist_map.begin())continue;


		//corrmap��ID
		id = itr->second->id;
		if (triangles.count(id) == 0) {
			fprintf(stderr, "alignment triangle ID=%d not found\n", id);
			exit(1);
		}
		//id�̑�����O�p�`��T��
		auto range = triangles.equal_range(id);
		for (auto itr2 = range.first; itr2 != range.second; itr2++) {
			sign[0] = (itr2->second->corr_p[1]->x - itr2->second->corr_p[0]->x) * (base.y - itr2->second->corr_p[1]->y) - (itr2->second->corr_p[1]->y - itr2->second->corr_p[0]->y) * (base.x - itr2->second->corr_p[1]->x);
			sign[1] = (itr2->second->corr_p[2]->x - itr2->second->corr_p[1]->x) * (base.y - itr2->second->corr_p[2]->y) - (itr2->second->corr_p[2]->y - itr2->second->corr_p[1]->y) * (base.x - itr2->second->corr_p[2]->x);
			sign[2] = (itr2->second->corr_p[0]->x - itr2->second->corr_p[2]->x) * (base.y - itr2->second->corr_p[0]->y) - (itr2->second->corr_p[0]->y - itr2->second->corr_p[2]->y) * (base.x - itr2->second->corr_p[0]->x);
			//printf("point %.lf,%.1lf\n", base.x, base.y);
			//printf("triangle %.1lf %.1lf\n", itr2->second->corr_p[0]->x, itr2->second->corr_p[0]->y);
			//printf("triangle %.1lf %.1lf\n", itr2->second->corr_p[1]->x, itr2->second->corr_p[1]->y);
			//printf("triangle %.1lf %.1lf\n", itr2->second->corr_p[2]->x, itr2->second->corr_p[2]->y);
			//printf("sign %.1lf %1.lf %.1lf\n", sign[0], sign[1], sign[2]);
			//printf("  signbit %d %d %d\n", std::signbit(sign[0]), std::signbit(sign[1]), std::signbit(sign[2]));
			//printf("n signbit %d %d %d\n", !std::signbit(sign[0]), !std::signbit(sign[1]), !std::signbit(sign[2]));
			//printf("judge %d\n", (std::signbit(sign[0]) && std::signbit(sign[1]) && std::signbit(sign[2])) || (!std::signbit(sign[0]) && !std::signbit(sign[1]) && !std::signbit(sign[2])));
			//printf("\n");

			//������3�Ƃ���v��true
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

	//dist���ŏ��ɂȂ�corrmap���Ƃ��Ă���
	dist = -1;
	for (auto itr = dist_map.begin(); itr != dist_map.end(); itr++) {
		//corrmap��ID
		id = itr->second->id;
		if (triangles.count(id) == 0) {
			fprintf(stderr, "alignment triangle ID=%d not found\n", id);
			exit(1);
		}
		//id�̑�����O�p�`��T��
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
std::vector <std::pair<basetrack_minimum_z*, corrmap_3d::align_param2*>>track_affineparam_correspondence(std::vector<basetrack_minimum_z>& base, std::vector <corrmap_3d::align_param2>& param) {

	//local align�̎��쒆�S�����o���āA�ʒu��hash
	//local align�̎��쒆�S�̍��delaunay�O�p�`��map�őΉ�

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

	std::vector < std::pair<basetrack_minimum_z*, corrmap_3d::align_param2*>> ret;
	std::vector<corrmap_3d::align_param*> param_cand;
	int loop = 0, ix, iy, count = 0;
	for (auto itr = base.begin(); itr != base.end(); itr++) {
		//if (count % 100000 == 0) {
		//	printf("\r search correspond triangles %d/%d(%4.1lf%%)", count, base.size(), count*100. / base.size());
		//}
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
		corrmap_3d::align_param2* param2 = search_param(param_cand, *itr, triangles);
		ret.push_back(std::make_pair(&(*itr), param2));
	}
	//printf("\r search correspond triangles %d/%d(%4.1lf%%)\n", count, base.size(), count*100. / base.size());

	return ret;
}


//�ϊ� zshrink�␳-->9para�ϊ�
void trans_base(std::vector<basetrack_minimum_z*>& base, corrmap_3d::align_param2* param) {

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
		//�p�xshrink���͂����ł�����
		base_p1.z = param->z + (base_thick) / param->z_shrink;

		//���쒆�S�����_�Ɉړ�
		//base_p0 = matrix_3D::addition(base_p0, matrix_3D::const_multiple(center, -1));
		//base_p1 = matrix_3D::addition(base_p1, matrix_3D::const_multiple(center, -1));

		//�ϊ��̎��s
		base_p0.matrix_multiplication(all_trans);
		base_p0 = matrix_3D::addition(base_p0, shift);
		base_p1.matrix_multiplication(all_trans);
		base_p1 = matrix_3D::addition(base_p1, shift);

		//���_�����Ƃɖ߂�
		//base_p0 = matrix_3D::addition(base_p0, center);
		//base_p1 = matrix_3D::addition(base_p1, center);

		(*itr)->x = base_p0.x;
		(*itr)->y = base_p0.y;
		(*itr)->z = base_p0.z;

		//printf("ax:%.4lf --> %.4lf\n", (*itr)->ax, (base_p0.x - base_p1.x) / (base_p0.z - base_p1.z));
		//printf("ay:%.4lf --> %.4lf\n", (*itr)->ay, (base_p0.y - base_p1.y) / (base_p0.z - base_p1.z));

		(*itr)->ax = (base_p0.x - base_p1.x) / (base_p0.z - base_p1.z) + param->zx_shear;
		(*itr)->ay = (base_p0.y - base_p1.y) / (base_p0.z - base_p1.z) + param->zy_shear;
		(*itr)->trkid = -1;
	}
}
void trans_base_all(std::vector < std::pair<basetrack_minimum_z*, corrmap_3d::align_param2*>>& track_pair) {
	std::map<std::tuple<int, int, int>, corrmap_3d::align_param2*> param_map;
	std::multimap<std::tuple<int, int, int>, basetrack_minimum_z*>base_map;
	std::tuple<int, int, int>id;
	//�O�p�`���Ƃ�basetrack���܂Ƃ߂�
	for (auto itr = track_pair.begin(); itr != track_pair.end(); itr++) {
		std::get<0>(id) = itr->second->corr_p[0]->id;
		std::get<1>(id) = itr->second->corr_p[1]->id;
		std::get<2>(id) = itr->second->corr_p[2]->id;
		param_map.insert(std::make_pair(id, itr->second));
		base_map.insert(std::make_pair(id, itr->first));
	}


	//�����ŎO�p�`���Ƃɕϊ�
	int count = 0;
	std::vector<basetrack_minimum_z*> t_base;
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


std::vector<basetrack_minimum_z> target_base_selection(std::vector<basetrack_minimum_z>& base, vxx::base_track_t muon_base, std::map<int, corrmap0::Corrmap>& corrmap_abs, int stop_pl, int target_pl) {


	if (corrmap_abs.count(stop_pl) == 0) {
		fprintf(stderr, "corrmap abs PL%03d not found\n", stop_pl);
		exit(1);
	}
	if (corrmap_abs.count(target_pl) == 0) {
		fprintf(stderr, "corrmap abs PL%03d not found\n", target_pl);
		exit(1);
	}

	auto corr_stop = corrmap_abs.at(stop_pl);
	auto corr_target = corrmap_abs.at(target_pl);

	//stop_pl��muon��stop pl�Atarget_pl��partner��T������pl�ł���B

	//��،N�̃R�[�h�ł͐F�X�ƃo�O���������B
	//1. stop_pl - target_pl == 0�̂Ƃ��A�������ɓS�w�ƌ��߂��Ă���A�������̉\������������Ă����B
	//2. stop_pl % 2 == 1��S�����ƌ��Ȃ��Ă����B�t�B0�̂Ƃ����S�̂͂��B
	//3. ���A�S�����킹�ĕ����w�ʉ߂��Ă���ꍇ�ł��A���̌��݂�K�؂ɍl�����Ă��Ȃ��B�t�B���������������Ă���B
	//�̂ŁA�������C������B
	int flg_iron = 0;
	int flg_water = 0;
	if (stop_pl - target_pl == 0) {
		//forward
		if (stop_pl <= 16 || stop_pl % 2 == 0) {
			//interaction in iron
			flg_iron = 1;
		}
		else {
			//interaction in water
			flg_water = 1;
		}
	}
	else if (stop_pl - target_pl == 1) {
		//forward +1
		if (stop_pl <= 16) {
			//interaction in/near iron ECC
			//�S2���ʉ߂���ƍl����B
			flg_iron = 2;
		}
		else {
			//�������A�S�����ǂ���ł����Ă����ƓS��1�w���ʉ߂���B
			flg_iron = 1;
			flg_water = 1;
		}
	}
	else if (stop_pl - target_pl == -1) {
		//backward
		if (stop_pl <= 16 || stop_pl % 2 == 0) {
			//interaction in iron
			//�S�w���㗬�Ɍ������Ēʉ߂����ƍl����B
			flg_iron = 1;
		}
		else {
			//interaction in water
			//���w���㗬�Ɍ������Ēʉ߂����ƍl����B
			flg_water = 1;
		}
	}
	else if (stop_pl - target_pl == -2) {
		//backward+1
		if (target_pl <= 17) {
			//interaction in/near iron ECC
			flg_iron = 2;
		}
		else {
			flg_iron = 1;
			flg_water = 1;
		}
	}
	else {
		fprintf(stderr, "PL range error muon_PL%03d target_PL%03d\n", stop_pl, target_pl);
		exit(1);
	}

	double muon_x, muon_y;
	muon_x = muon_base.x * corr_stop.position[0] + muon_base.y * corr_stop.position[1] + corr_stop.position[4];
	muon_y = muon_base.x * corr_stop.position[2] + muon_base.y * corr_stop.position[3] + corr_stop.position[5];

	double target_x, target_y;
	std::vector<basetrack_minimum_z> ret;

	double range = flg_water * 42000. + flg_iron * 8000.;
	for (auto& b : base) {
		target_x = b.x * corr_target.position[0] + b.y * corr_target.position[1] + corr_target.position[4];
		target_y = b.x * corr_target.position[2] + b.y * corr_target.position[3] + corr_target.position[5];

		if (fabs(target_x - muon_x) > range) continue;
		if (fabs(target_y - muon_y) > range) continue;
		ret.push_back(b);
	}

	return ret;

}

std::vector<basetrack_minimum_z> base_minimum_multidel(std::vector<basetrack_minimum_z>& b) {
	std::vector<basetrack_minimum_z> ret;

	std::map<std::pair<int, int>, basetrack_minimum_z>b_map;
	for (auto itr = b.begin(); itr != b.end(); itr++) {
		int pl = itr->pl;
		int rawid = itr->rawid;
		auto res = b_map.insert(std::make_pair(std::make_pair(pl, rawid), *itr));
		if (!res.second) {
			if (res.first->second.black_flg == 1)continue;
			else {
				if (itr->black_flg == 1)res.first->second.black_flg = 1;
			}
		}
	}

	for (auto itr = b_map.begin(); itr != b_map.end(); itr++) {
		ret.push_back(itr->second);
	}
	return ret;

}


void partner_search(vxx::base_track_t target, std::vector<basetrack_minimum_z>& all, double z_range[2], double angle_accuracy_intercept_mu, double angle_accuracy_slope_mu, double angle_accuracy_intercept, double angle_accuracy_slope, std::vector<basetrack_minimum_z>& connect, std::ofstream& ofs, int black_flg) {

	matrix_3D::vector_3D pos, dir;
	pos.x = target.x;
	pos.y = target.y;
	pos.z = 0;

	dir.x = target.ax;
	dir.y = target.ay;
	dir.z = 1;

	double muon_angle_acc = angle_accuracy_intercept_mu + angle_accuracy_slope_mu * sqrt(dir.x * dir.x + dir.y * dir.y);
	double partner_angle_acc;


	double md, extra[2], allowance, angle;
	matrix_3D::vector_3D pos_all, dir_all;
	for (auto itr = all.begin(); itr != all.end(); itr++) {
		if (itr->pl == target.pl && itr->rawid == target.rawid)continue;
		pos_all.x = itr->x;
		pos_all.y = itr->y;
		pos_all.z = itr->z;
		dir_all.x = itr->ax;
		dir_all.y = itr->ay;
		dir_all.z = 1;
		if (black_flg == 0) {
			if (fabs(dir_all.x) > 4.5)continue;
			if (fabs(dir_all.y) > 4.5)continue;
		}
		md = matrix_3D::minimum_distance(pos, pos_all, dir, dir_all, z_range, extra);
		partner_angle_acc = angle_accuracy_intercept + angle_accuracy_slope * sqrt(dir_all.x * dir_all.x + dir_all.y * dir_all.y);

		allowance = pow(fabs(muon_angle_acc * extra[0]) + 5, 2) + pow(fabs(partner_angle_acc * extra[1]) + 5, 2);
		if (md * md < allowance) {
			itr->ex_z0 = pos_all.z;
			itr->ex_z1 = pos_all.z + extra[1];
			itr->md = md;
			connect.push_back(*itr);
			ofs << std::right << std::fixed
				<< std::setw(4) << std::setprecision(0) << target.pl << " "
				<< std::setw(10) << std::setprecision(0) << target.rawid << " "
				<< std::setw(7) << std::setprecision(4) << target.ax << " "
				<< std::setw(7) << std::setprecision(4) << target.ay << " "
				<< std::setw(8) << std::setprecision(1) << target.x << " "
				<< std::setw(8) << std::setprecision(1) << target.y << " "
				<< std::setw(8) << std::setprecision(1) << target.z << " "
				<< std::setw(8) << std::setprecision(1) << extra[0] << " "
				<< std::setw(4) << std::setprecision(0) << itr->pl << " "
				<< std::setw(10) << std::setprecision(0) << itr->rawid << " "
				<< std::setw(7) << std::setprecision(4) << itr->ax << " "
				<< std::setw(7) << std::setprecision(4) << itr->ay << " "
				<< std::setw(8) << std::setprecision(1) << itr->x << " "
				<< std::setw(8) << std::setprecision(1) << itr->y << " "
				<< std::setw(8) << std::setprecision(1) << itr->z << " "
				<< std::setw(7) << std::setprecision(1) << extra[1] << " "
				<< std::setw(2) << std::setprecision(0) << black_flg << " "
				<< std::setw(6) << std::setprecision(2) << md << " "
				<< std::setw(6) << std::setprecision(2) << sqrt(allowance) << std::endl;
		}
	}


}
void partner_search_single(vxx::base_track_t target, basetrack_minimum_z& all, double z_range[2], double angle_accuracy_intercept_mu, double angle_accuracy_slope_mu, double angle_accuracy_intercept, double angle_accuracy_slope, std::vector<basetrack_minimum_z>& connect, std::ofstream& ofs, int black_flg) {

	matrix_3D::vector_3D pos, dir;
	pos.x = target.x;
	pos.y = target.y;
	pos.z = 0;

	dir.x = target.ax;
	dir.y = target.ay;
	dir.z = 1;

	double muon_angle_acc = angle_accuracy_intercept_mu + angle_accuracy_slope_mu * sqrt(dir.x * dir.x + dir.y * dir.y);
	double partner_angle_acc;


	double md, extra[2], allowance, angle;
	matrix_3D::vector_3D pos_all, dir_all;
	pos_all.x = all.x;
	pos_all.y = all.y;
	pos_all.z = all.z;
	dir_all.x = all.ax;
	dir_all.y = all.ay;
	dir_all.z = 1;

	md = matrix_3D::minimum_distance(pos, pos_all, dir, dir_all, z_range, extra);
	partner_angle_acc = angle_accuracy_intercept + angle_accuracy_slope * sqrt(dir_all.x * dir_all.x + dir_all.y * dir_all.y);

	allowance = pow(fabs(muon_angle_acc * extra[0]) + 5, 2) + pow(fabs(partner_angle_acc * extra[1]) + 5, 2);
	if (md * md < allowance) {
		all.ex_z0 = pos_all.z;
		all.ex_z1 = pos_all.z + extra[1];
		all.md = md;
		connect.push_back(all);
		ofs << std::right << std::fixed
			<< std::setw(4) << std::setprecision(0) << target.pl << " "
			<< std::setw(10) << std::setprecision(0) << target.rawid << " "
			<< std::setw(7) << std::setprecision(4) << target.ax << " "
			<< std::setw(7) << std::setprecision(4) << target.ay << " "
			<< std::setw(8) << std::setprecision(1) << target.x << " "
			<< std::setw(8) << std::setprecision(1) << target.y << " "
			<< std::setw(8) << std::setprecision(1) << target.z << " "
			<< std::setw(8) << std::setprecision(1) << extra[0] << " "
			<< std::setw(4) << std::setprecision(0) << all.pl << " "
			<< std::setw(10) << std::setprecision(0) << all.rawid << " "
			<< std::setw(7) << std::setprecision(4) << all.ax << " "
			<< std::setw(7) << std::setprecision(4) << all.ay << " "
			<< std::setw(8) << std::setprecision(1) << all.x << " "
			<< std::setw(8) << std::setprecision(1) << all.y << " "
			<< std::setw(8) << std::setprecision(1) << all.z << " "
			<< std::setw(7) << std::setprecision(1) << extra[1] << " "
			<< std::setw(2) << std::setprecision(0) << black_flg << " "
			<< std::setw(6) << std::setprecision(2) << md << " "
			<< std::setw(6) << std::setprecision(2) << sqrt(allowance) << std::endl;
	}
}

void pre_partner_search(vxx::base_track_t target, std::vector<basetrack_minimum_z>& all, double z_range[2], double angle_accuracy_intercept_mu, double angle_accuracy_slope_mu, double angle_accuracy_intercept, double angle_accuracy_slope, std::vector<std::pair<basetrack_minimum_z, basetrack_minimum_z>>& attach_candidate, int black_flg) {

	matrix_3D::vector_3D pos, dir;
	pos.x = target.x;
	pos.y = target.y;
	pos.z = 0;

	dir.x = target.ax;
	dir.y = target.ay;
	dir.z = 1;

	double muon_angle_acc = angle_accuracy_intercept_mu + angle_accuracy_slope_mu * sqrt(dir.x * dir.x + dir.y * dir.y);
	double partner_angle_acc;
	basetrack_minimum_z md_point;

	double md, extra[2], allowance, angle;
	matrix_3D::vector_3D pos_all, dir_all;
	for (auto itr = all.begin(); itr != all.end(); itr++) {
		if (itr->pl == target.pl && itr->rawid == target.rawid)continue;
		pos_all.x = itr->x;
		pos_all.y = itr->y;
		pos_all.z = itr->z;
		dir_all.x = itr->ax;
		dir_all.y = itr->ay;
		dir_all.z = 1;
		if (black_flg == 0) {
			if (fabs(dir_all.x) > 4.5)continue;
			if (fabs(dir_all.y) > 4.5)continue;
		}
		md = matrix_3D::minimum_distance(pos, pos_all, dir, dir_all, z_range, extra);
		partner_angle_acc = angle_accuracy_intercept + angle_accuracy_slope * sqrt(dir_all.x * dir_all.x + dir_all.y * dir_all.y);

		allowance = pow(fabs(muon_angle_acc * extra[0]) + 5, 2) + pow(fabs(partner_angle_acc * extra[1]) + 5, 2);
		if (md * md < allowance) {
			md_point = *itr;
			//md��g��xy���W�̎擾
			matrix_3D::vector_3D extra0 = matrix_3D::addition(pos, const_multiple(dir, extra[0]));
			matrix_3D::vector_3D extra1 = matrix_3D::addition(pos_all, const_multiple(dir_all, extra[1]));

			md_point.x = (extra0.x + extra1.x) / 2;
			md_point.y = (extra0.y + extra1.y) / 2;
			md_point.z = 0;

			attach_candidate.push_back(std::make_pair(*itr, md_point));


		}
	}


}

void Set_search_z(double z_range[2], int muonPL) {
	bool iron = false;
	bool water = false;
	bool other = false;
	if (muonPL <= 15)iron = true;
	else if (16 <= muonPL && muonPL % 2 == 0)iron = true;
	else if (16 <= muonPL && muonPL % 2 == 1)water = true;
	else other = true;

	if (water) {
		//water-emulsion-iron-emulsion-water- (+a)
		z_range[0] = -2300 - 350 - 110 * 2 - 1000;
		//z_range[1] = +210 + 70 + 30;
		z_range[1] = 0;
	}
	else {
		//iron-emulsion-water- (+a)
		z_range[0] = -500 - 350 - 70 - 50;
		//z_range[1] = +210 + 70 + 30;
		z_range[1] = 0;
	}
}

void output_Track(std::string filename, std::vector<Track_file>& track) {
	std::ofstream ofs(filename);
	for (auto itr = track.begin(); itr != track.end(); itr++) {
		ofs << std::fixed << std::right
			<< std::setw(10) << std::setprecision(0) << itr->eventid << " "
			<< std::setw(5) << std::setprecision(0) << itr->trackid << " "
			<< std::setw(4) << std::setprecision(0) << itr->pl << " "
			<< std::setw(10) << std::setprecision(0) << itr->rawid << std::endl;
	}

}

void output_Track_inf(std::string filename, std::vector<output_track>& track) {
	std::ofstream ofs(filename);
	for (auto itr = track.begin(); itr != track.end(); itr++) {
		ofs << std::fixed << std::right
			<< std::setw(10) << std::setprecision(0) << itr->eventid << " "
			<< std::setw(5) << std::setprecision(0) << itr->partner.size() + 1 << std::endl;
		ofs << std::fixed << std::right
			<< std::setw(3) << std::setprecision(0) << itr->muon.trkid << " "
			<< std::setw(3) << std::setprecision(0) << itr->muon.pl << " "
			<< std::setw(10) << std::setprecision(0) << itr->muon.rawid << " "
			<< std::setw(7) << std::setprecision(0) << itr->muon.ph << " "
			<< std::setw(3) << std::setprecision(0) << itr->muon.black_flg << " "
			<< std::setw(8) << std::setprecision(4) << itr->muon.ax << " "
			<< std::setw(8) << std::setprecision(4) << itr->muon.ay << " "
			<< std::setw(10) << std::setprecision(1) << itr->muon.x << " "
			<< std::setw(10) << std::setprecision(1) << itr->muon.y << " "
			<< std::setw(10) << std::setprecision(1) << itr->muon.z << " "
			<< std::setw(10) << std::setprecision(1) << itr->muon.ex_z0 << " "
			<< std::setw(10) << std::setprecision(1) << itr->muon.ex_z1 << " "
			<< std::setw(10) << std::setprecision(1) << itr->muon.md << std::endl;
		for (int i = 0; i < itr->partner.size(); i++) {
			ofs << std::fixed << std::right
				<< std::setw(3) << std::setprecision(0) << itr->partner[i].trkid << " "
				<< std::setw(3) << std::setprecision(0) << itr->partner[i].pl << " "
				<< std::setw(10) << std::setprecision(0) << itr->partner[i].rawid << " "
				<< std::setw(7) << std::setprecision(0) << itr->partner[i].ph << " "
				<< std::setw(3) << std::setprecision(0) << itr->partner[i].black_flg << " "
				<< std::setw(8) << std::setprecision(4) << itr->partner[i].ax << " "
				<< std::setw(8) << std::setprecision(4) << itr->partner[i].ay << " "
				<< std::setw(10) << std::setprecision(1) << itr->partner[i].x << " "
				<< std::setw(10) << std::setprecision(1) << itr->partner[i].y << " "
				<< std::setw(10) << std::setprecision(1) << itr->partner[i].z << " "
				<< std::setw(10) << std::setprecision(1) << itr->partner[i].ex_z0 << " "
				<< std::setw(10) << std::setprecision(1) << itr->partner[i].ex_z1 << " "
				<< std::setw(10) << std::setprecision(1) << itr->partner[i].md << std::endl;
		}
	}

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