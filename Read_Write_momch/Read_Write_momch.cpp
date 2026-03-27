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
#include <map>

#include <list>
#include <set>
#include <unordered_map>
#include <utility>
#include <random>
#include <cmath>
#include <limits>
#include <chrono>
#include <cassert>

//Class定義
namespace Momentum_recon {

	class microtrack_minimum {
	public:
		int ph, zone, view, imager, pixelnum, hitnum;
		microtrack_minimum();
	};
	class Mom_basetrack {
	public:
		int pl, rawid;
		float ax, ay, x, y, z;
		microtrack_minimum m[2];
		Mom_basetrack();
	};
	class Mom_chain {
	public:
		//stop_flg=0 : penetarte / 1:bm stop / 2:ecc stop
		int chainid, stop_flg, particle_flg, direction, charge_sign;
		//[0]:muon [1]:proton
		double ecc_range_mom[2], ecc_mcs_mom[2], bm_range_mom, bm_curvature_mom;
		double ecc_range_mom_error[2][2], ecc_mcs_mom_error[2][2], bm_range_mom_error[2], bm_curvature_mom_error[2];
		//muon_likelihood   =(vph-  muon_mean)/  muon_sigma (符号付き)
		//proton_likelihood =(vph-proton_mean)/proton_sigma (符号付き)
		double muon_likelihood, proton_likelihood;
		//mfile座標系
		std::vector<Mom_basetrack> base;
		//localな座標系
		std::vector<std::pair<Mom_basetrack, Mom_basetrack>> base_pair;

		//コンストラクタ
		Mom_chain();
		double Get_muon_mcs_pb();
		void Get_muon_pb_mcs_error(double *pb_error);
		double Get_proton_mcs_pb();
		void Get_proton_pb_mcs_error(double *pb_error);
	};
	class Event_information {
	public:
		int groupid, unix_time, tracker_track_id, entry_in_daily_file, vertex_pl, ECC_id, vertex_material;
		double vertex_position[3], true_vertex_position[3], weight, nu_energy, nu_ax, nu_ay;
		std::vector<Mom_chain> chains;
		std::vector<Mom_chain> true_chains;

		//コンストラクタ
		Event_information();
		//0:water,2:iron -1:ini
	};

	void Write_Event_information_header(std::ofstream &ofs, Event_information&ev);
	bool Read_Event_information_header(std::ifstream &ifs, Event_information&ev, int &chian_num, int &true_chian_num);
	void Write_mom_chain_header(std::ofstream &ofs, Mom_chain&mom_chain);
	bool Read_mom_chain_header(std::ifstream &ifs, Mom_chain&mom_chain, int &base_num, int &base_pair_num);
	void Write_Event_information_bin(std::string filename, std::vector<Event_information>&ev_v);
	std::vector<Event_information> Read_Event_information_bin(std::string filename);
	void Write_Event_information_bin_block(std::ofstream &ofs, std::vector<Event_information>&ev_v);
	std::vector<Event_information> Read_Event_information_bin_block(std::ifstream &ifs, int max_event_num);

	void Write_Event_information_txt(std::string filename, std::vector<Event_information>&ev_v);
	std::vector<Event_information> Read_Event_information_txt(std::string filename);
	void input_basetrack_information(std::vector<Mom_basetrack> &base, std::vector<std::pair<Mom_basetrack, Mom_basetrack>> &base_pair);
	std::vector<Event_information> Read_Event_information_extension(std::string filename);
	void Write_Event_information_extension(std::string filename, std::vector<Event_information>&ev_v);
}

//実装
namespace Momentum_recon {
	void Write_Event_information_header(std::ofstream &ofs, Event_information&ev) {
		int chian_num = ev.chains.size();
		int true_chian_num = ev.true_chains.size();

		ofs.write((char*)& ev.groupid, sizeof(int));
		ofs.write((char*)& ev.unix_time, sizeof(int));
		ofs.write((char*)& ev.tracker_track_id, sizeof(int));
		ofs.write((char*)& ev.entry_in_daily_file, sizeof(int));
		ofs.write((char*)& ev.vertex_pl, sizeof(int));
		ofs.write((char*)& ev.ECC_id, sizeof(int));
		ofs.write((char*)& ev.vertex_material, sizeof(int));
		ofs.write((char*)& ev.vertex_position, sizeof(double) * 3);
		ofs.write((char*)& ev.true_vertex_position, sizeof(double) * 3);
		ofs.write((char*)& ev.weight, sizeof(double));
		ofs.write((char*)& ev.nu_energy, sizeof(double));
		ofs.write((char*)& ev.nu_ax, sizeof(double));
		ofs.write((char*)& ev.nu_ay, sizeof(double));
		ofs.write((char*)& chian_num, sizeof(int));
		ofs.write((char*)& true_chian_num, sizeof(int));

	}
	bool Read_Event_information_header(std::ifstream &ifs, Event_information&ev, int &chian_num, int &true_chian_num) {

		if (!ifs.read((char*)& ev.groupid, sizeof(int)))return false;
		if (!ifs.read((char*)& ev.unix_time, sizeof(int)))return false;
		if (!ifs.read((char*)& ev.tracker_track_id, sizeof(int)))return false;
		if (!ifs.read((char*)& ev.entry_in_daily_file, sizeof(int)))return false;
		if (!ifs.read((char*)& ev.vertex_pl, sizeof(int)))return false;
		if (!ifs.read((char*)& ev.ECC_id, sizeof(int)))return false;
		if (!ifs.read((char*)& ev.vertex_material, sizeof(int)))return false;
		if (!ifs.read((char*)& ev.vertex_position, sizeof(double) * 3))return false;
		if (!ifs.read((char*)& ev.true_vertex_position, sizeof(double) * 3))return false;
		if (!ifs.read((char*)& ev.weight, sizeof(double)))return false;
		if (!ifs.read((char*)& ev.nu_energy, sizeof(double)))return false;
		if (!ifs.read((char*)& ev.nu_ax, sizeof(double)))return false;
		if (!ifs.read((char*)& ev.nu_ay, sizeof(double)))return false;
		if (!ifs.read((char*)& chian_num, sizeof(int)))return false;
		if (!ifs.read((char*)& true_chian_num, sizeof(int)))return false;

		return true;
	}

	void Write_mom_chain_header(std::ofstream &ofs, Mom_chain&mom_chain) {
		int base_num = mom_chain.base.size();
		int base_pair_num = mom_chain.base_pair.size();

		ofs.write((char*)& mom_chain.chainid, sizeof(int));
		ofs.write((char*)& mom_chain.stop_flg, sizeof(int));
		ofs.write((char*)& mom_chain.particle_flg, sizeof(int));
		ofs.write((char*)& mom_chain.direction, sizeof(int));
		ofs.write((char*)& mom_chain.charge_sign, sizeof(int));
		ofs.write((char*)& mom_chain.ecc_range_mom, sizeof(double) * 2);
		ofs.write((char*)& mom_chain.ecc_mcs_mom, sizeof(double) * 2);
		ofs.write((char*)& mom_chain.bm_range_mom, sizeof(double));
		ofs.write((char*)& mom_chain.bm_curvature_mom, sizeof(double));
		ofs.write((char*)& mom_chain.ecc_range_mom_error, sizeof(double) * 2 * 2);
		ofs.write((char*)& mom_chain.ecc_mcs_mom_error, sizeof(double) * 2 * 2);
		ofs.write((char*)& mom_chain.bm_range_mom_error, sizeof(double) * 2);
		ofs.write((char*)& mom_chain.bm_curvature_mom_error, sizeof(double) * 2);
		ofs.write((char*)& mom_chain.muon_likelihood, sizeof(double));
		ofs.write((char*)& mom_chain.proton_likelihood, sizeof(double));
		ofs.write((char*)& base_num, sizeof(int));
		ofs.write((char*)& base_pair_num, sizeof(int));
	}
	bool Read_mom_chain_header(std::ifstream &ifs, Mom_chain&mom_chain, int &base_num, int &base_pair_num) {
		if (!ifs.read((char*)& mom_chain.chainid, sizeof(int)))return false;
		if (!ifs.read((char*)& mom_chain.stop_flg, sizeof(int)))return false;
		if (!ifs.read((char*)& mom_chain.particle_flg, sizeof(int)))return false;
		if (!ifs.read((char*)& mom_chain.direction, sizeof(int)))return false;
		if (!ifs.read((char*)& mom_chain.charge_sign, sizeof(int)))return false;
		if (!ifs.read((char*)& mom_chain.ecc_range_mom, sizeof(double) * 2))return false;
		if (!ifs.read((char*)& mom_chain.ecc_mcs_mom, sizeof(double) * 2))return false;
		if (!ifs.read((char*)& mom_chain.bm_range_mom, sizeof(double)))return false;
		if (!ifs.read((char*)& mom_chain.bm_curvature_mom, sizeof(double)))return false;
		if (!ifs.read((char*)& mom_chain.ecc_range_mom_error, sizeof(double) * 2 * 2))return false;
		if (!ifs.read((char*)& mom_chain.ecc_mcs_mom_error, sizeof(double) * 2 * 2))return false;
		if (!ifs.read((char*)& mom_chain.bm_range_mom_error, sizeof(double) * 2))return false;
		if (!ifs.read((char*)& mom_chain.bm_curvature_mom_error, sizeof(double) * 2))return false;
		if (!ifs.read((char*)& mom_chain.muon_likelihood, sizeof(double)))return false;
		if (!ifs.read((char*)& mom_chain.proton_likelihood, sizeof(double)))return false;
		if (!ifs.read((char*)& base_num, sizeof(int)))return false;
		if (!ifs.read((char*)& base_pair_num, sizeof(int)))return false;
		mom_chain.base.reserve(base_num);
		mom_chain.base_pair.reserve(base_pair_num);
		return true;
	}

	void Write_Event_information_bin(std::string filename, std::vector<Event_information>&ev_v) {
		std::ofstream ofs(filename, std::ios::binary);
		if (!ofs) {
			//file open 失敗
			fprintf(stderr, "File[%s] is not exist!!\n", filename.c_str());
			exit(1);
		}
		if (ev_v.size() == 0) {
			fprintf(stderr, "target data ... null\n");
			fprintf(stderr, "File[%s] has no text\n", filename.c_str());
		}
		else {
			int64_t count = 0;
			int64_t max = ev_v.size();
			int base_num, base_pair_num;
			int chain_num, true_chain_num;

			for (int i = 0; i < max; i++) {
				if (count % 100 == 0) {
					std::cerr << std::right << std::fixed << "\r now writing ..." << std::setw(4) << std::setprecision(1) << count * 100. / max << "%%";
				}
				count++;
				Write_Event_information_header(ofs, ev_v[i]);
				chain_num = ev_v[i].chains.size();
				true_chain_num = ev_v[i].true_chains.size();
				for (int j = 0; j < chain_num; j++) {
					base_num = ev_v[i].chains[j].base.size();
					base_pair_num = ev_v[i].chains[j].base_pair.size();
					Write_mom_chain_header(ofs, ev_v[i].chains[j]);

					for (int k = 0; k < base_num; k++) {
						ofs.write((char*)& ev_v[i].chains[j].base[k], sizeof(Mom_basetrack));
					}
					for (int k = 0; k < base_pair_num; k++) {
						ofs.write((char*)& ev_v[i].chains[j].base_pair[k].first, sizeof(Mom_basetrack));
						ofs.write((char*)& ev_v[i].chains[j].base_pair[k].second, sizeof(Mom_basetrack));
					}
				}
				for (int j = 0; j < true_chain_num; j++) {
					base_num = ev_v[i].true_chains[j].base.size();
					base_pair_num = ev_v[i].true_chains[j].base_pair.size();
					Write_mom_chain_header(ofs, ev_v[i].true_chains[j]);

					for (int k = 0; k < base_num; k++) {
						ofs.write((char*)& ev_v[i].true_chains[j].base[k], sizeof(Mom_basetrack));
					}
					for (int k = 0; k < base_pair_num; k++) {
						ofs.write((char*)& ev_v[i].true_chains[j].base_pair[k].first, sizeof(Mom_basetrack));
						ofs.write((char*)& ev_v[i].true_chains[j].base_pair[k].second, sizeof(Mom_basetrack));
					}
				}
			}
			std::cerr << std::right << std::fixed << "\r now writing ..." << std::setw(4) << std::setprecision(1) << count * 100. / max << "%%" << std::endl;
		}
	}
	std::vector<Event_information> Read_Event_information_bin(std::string filename) {
		std::vector<Event_information> ret;
		std::ifstream ifs(filename, std::ios::binary);
		//filesize取得
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
		Mom_chain t;
		Mom_basetrack m;
		std::pair<Mom_basetrack, Mom_basetrack> p;
		int base_num, base_pair_num;
		int chain_num, true_chain_num;
		int header[5];
		double mom_recon;
		Event_information ev;
		while (Read_Event_information_header(ifs, ev, chain_num, true_chain_num)) {

			if (count % 100 == 0) {
				nowpos = ifs.tellg();
				auto size1 = nowpos - begpos;
				std::cerr << std::right << std::fixed << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%";
			}
			count++;
			for (int i_ch = 0; i_ch < chain_num; i_ch++) {
				if (Read_mom_chain_header(ifs, t, base_num, base_pair_num)) {

					t.base.clear();
					t.base_pair.clear();
					t.base.reserve(base_num);
					t.base_pair.reserve(base_pair_num);
					for (int j = 0; j < base_num; j++) {
						ifs.read((char*)& m, sizeof(Mom_basetrack));
						t.base.push_back(m);
					}
					for (int j = 0; j < base_pair_num; j++) {
						ifs.read((char*)& p.first, sizeof(Mom_basetrack));
						ifs.read((char*)& p.second, sizeof(Mom_basetrack));
						t.base_pair.push_back(p);
					}
					ev.chains.push_back(t);
				}
			}
			for (int i_ch = 0; i_ch < true_chain_num; i_ch++) {
				if (Read_mom_chain_header(ifs, t, base_num, base_pair_num)) {
					t.base.clear();
					t.base_pair.clear();
					t.base.reserve(base_num);
					t.base_pair.reserve(base_pair_num);
					for (int j = 0; j < base_num; j++) {
						ifs.read((char*)& m, sizeof(Mom_basetrack));
						t.base.push_back(m);
					}
					for (int j = 0; j < base_pair_num; j++) {
						ifs.read((char*)& p.first, sizeof(Mom_basetrack));
						ifs.read((char*)& p.second, sizeof(Mom_basetrack));
						t.base_pair.push_back(p);
					}
					ev.true_chains.push_back(t);
				}
			}
			ret.push_back(ev);
			ev.chains.clear();
			ev.true_chains.clear();
		}

		auto size1 = eofpos - begpos;
		std::cerr << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%" << std::endl;;
		if (count == 0) {
			fprintf(stderr, "%s no data!\n", filename.c_str());
			exit(1);
		}
		return ret;

	}

	void Write_Event_information_txt(std::string filename, std::vector<Event_information>&ev_v) {

		std::ofstream ofs(filename);
		int count = 0, all = ev_v.size();
		for (auto &ev : ev_v) {
			if (count % 100 == 0) {
				fprintf(stderr, "\r write group %d/%d(%4.1lf%%)", count, all, count*100. / all);
			}
			count++;
			//event inforamtion header 書き出し
			ofs << std::right << std::fixed
				<< std::setw(10) << std::setprecision(0) << ev.groupid << " "
				<< std::setw(10) << std::setprecision(0) << ev.unix_time << " "
				<< std::setw(3) << std::setprecision(0) << ev.tracker_track_id << " "
				<< std::setw(5) << std::setprecision(0) << ev.entry_in_daily_file << " "
				<< std::setw(3) << std::setprecision(0) << ev.vertex_pl << " "
				<< std::setw(2) << std::setprecision(0) << ev.ECC_id << " "
				<< std::setw(2) << std::setprecision(0) << ev.vertex_material << " "
				<< std::setw(8) << std::setprecision(1) << ev.vertex_position[0] << " "
				<< std::setw(8) << std::setprecision(1) << ev.vertex_position[1] << " "
				<< std::setw(8) << std::setprecision(1) << ev.vertex_position[2] << " "
				<< std::setw(8) << std::setprecision(1) << ev.true_vertex_position[0] << " "
				<< std::setw(8) << std::setprecision(1) << ev.true_vertex_position[1] << " "
				<< std::setw(8) << std::setprecision(1) << ev.true_vertex_position[2] << " "
				<< std::setw(7) << std::setprecision(1) << ev.nu_energy << " "
				<< std::setw(7) << std::setprecision(4) << ev.nu_ax << " "
				<< std::setw(7) << std::setprecision(4) << ev.nu_ay << " "
				<< std::setw(5) << std::setprecision(0) << ev.chains.size() << " "
				<< std::setw(5) << std::setprecision(0) << ev.true_chains.size() << std::endl;
			for (int i = 0; i < ev.chains.size(); i++) {
				ofs << std::right << std::fixed
					<< std::setw(10) << std::setprecision(0) << ev.chains[i].chainid << " "
					<< std::setw(5) << std::setprecision(0) << ev.chains[i].stop_flg << " "
					<< std::setw(5) << std::setprecision(0) << ev.chains[i].particle_flg << " "
					<< std::setw(3) << std::setprecision(0) << ev.chains[i].direction << " "
					<< std::setw(2) << std::setprecision(0) << ev.chains[i].charge_sign << " "
					<< std::setw(7) << std::setprecision(1) << ev.chains[i].ecc_range_mom[0] << " "
					<< std::setw(7) << std::setprecision(1) << ev.chains[i].ecc_range_mom[1] << " "
					<< std::setw(7) << std::setprecision(1) << ev.chains[i].ecc_mcs_mom[0] << " "
					<< std::setw(7) << std::setprecision(1) << ev.chains[i].ecc_mcs_mom[1] << " "
					<< std::setw(7) << std::setprecision(1) << ev.chains[i].bm_range_mom << " "
					<< std::setw(7) << std::setprecision(1) << ev.chains[i].bm_curvature_mom << " "
					<< std::setw(7) << std::setprecision(1) << ev.chains[i].ecc_range_mom_error[0][0] << " "
					<< std::setw(7) << std::setprecision(1) << ev.chains[i].ecc_range_mom_error[0][1] << " "
					<< std::setw(7) << std::setprecision(1) << ev.chains[i].ecc_range_mom_error[1][0] << " "
					<< std::setw(7) << std::setprecision(1) << ev.chains[i].ecc_range_mom_error[1][1] << " "
					<< std::setw(7) << std::setprecision(1) << ev.chains[i].ecc_mcs_mom_error[0][0] << " "
					<< std::setw(7) << std::setprecision(1) << ev.chains[i].ecc_mcs_mom_error[0][1] << " "
					<< std::setw(7) << std::setprecision(1) << ev.chains[i].ecc_mcs_mom_error[1][0] << " "
					<< std::setw(7) << std::setprecision(1) << ev.chains[i].ecc_mcs_mom_error[1][1] << " "
					<< std::setw(7) << std::setprecision(1) << ev.chains[i].bm_range_mom_error[0] << " "
					<< std::setw(7) << std::setprecision(1) << ev.chains[i].bm_range_mom_error[1] << " "
					<< std::setw(7) << std::setprecision(1) << ev.chains[i].bm_curvature_mom_error[0] << " "
					<< std::setw(7) << std::setprecision(1) << ev.chains[i].bm_curvature_mom_error[1] << " "
					<< std::setw(7) << std::setprecision(4) << ev.chains[i].muon_likelihood << " "
					<< std::setw(7) << std::setprecision(4) << ev.chains[i].proton_likelihood << " "
					<< std::setw(3) << std::setprecision(0) << ev.chains[i].base.size() << " "
					<< std::setw(3) << std::setprecision(0) << ev.chains[i].base_pair.size() << std::endl;

				for (auto &b : ev.chains[i].base) {
					ofs << std::right << std::fixed
						<< std::setw(4) << std::setprecision(0) << b.pl << " "
						<< std::setw(10) << std::setprecision(0) << b.rawid << " "
						<< std::setw(7) << std::setprecision(4) << b.ax << " "
						<< std::setw(7) << std::setprecision(4) << b.ay << " "
						<< std::setw(8) << std::setprecision(1) << b.x << " "
						<< std::setw(8) << std::setprecision(1) << b.y << " "
						<< std::setw(8) << std::setprecision(1) << b.z << " "
						<< std::setw(3) << std::setprecision(0) << b.m[0].zone << " "
						<< std::setw(4) << std::setprecision(0) << b.m[0].view << " "
						<< std::setw(3) << std::setprecision(0) << b.m[0].imager << " "
						<< std::setw(7) << std::setprecision(0) << b.m[0].ph << " "
						<< std::setw(5) << std::setprecision(0) << b.m[0].pixelnum << " "
						<< std::setw(5) << std::setprecision(0) << b.m[0].hitnum << " "
						<< std::setw(3) << std::setprecision(0) << b.m[1].zone << " "
						<< std::setw(4) << std::setprecision(0) << b.m[1].view << " "
						<< std::setw(3) << std::setprecision(0) << b.m[1].imager << " "
						<< std::setw(7) << std::setprecision(0) << b.m[1].ph << " "
						<< std::setw(5) << std::setprecision(0) << b.m[1].pixelnum << " "
						<< std::setw(5) << std::setprecision(0) << b.m[1].hitnum << std::endl;
				}
				for (auto &p : ev.chains[i].base_pair) {
					ofs << std::right << std::fixed
						<< std::setw(4) << std::setprecision(0) << p.first.pl << " "
						<< std::setw(10) << std::setprecision(0) << p.first.rawid << " "
						<< std::setw(4) << std::setprecision(0) << p.second.pl << " "
						<< std::setw(10) << std::setprecision(0) << p.second.rawid << " "
						<< std::setw(7) << std::setprecision(4) << p.first.ax << " "
						<< std::setw(7) << std::setprecision(4) << p.first.ay << " "
						<< std::setw(8) << std::setprecision(1) << p.first.x << " "
						<< std::setw(8) << std::setprecision(1) << p.first.y << " "
						<< std::setw(8) << std::setprecision(1) << p.first.z << " "
						<< std::setw(7) << std::setprecision(4) << p.second.ax << " "
						<< std::setw(7) << std::setprecision(4) << p.second.ay << " "
						<< std::setw(8) << std::setprecision(1) << p.second.x << " "
						<< std::setw(8) << std::setprecision(1) << p.second.y << " "
						<< std::setw(8) << std::setprecision(1) << p.second.z << std::endl;
				}

			}

			for (int i = 0; i < ev.true_chains.size(); i++) {
				ofs << std::right << std::fixed
					<< std::setw(10) << std::setprecision(0) << ev.true_chains[i].chainid << " "
					<< std::setw(5) << std::setprecision(0) << ev.true_chains[i].stop_flg << " "
					<< std::setw(5) << std::setprecision(0) << ev.true_chains[i].particle_flg << " "
					<< std::setw(3) << std::setprecision(0) << ev.true_chains[i].direction << " "
					<< std::setw(2) << std::setprecision(0) << ev.true_chains[i].charge_sign << " "
					<< std::setw(7) << std::setprecision(1) << ev.true_chains[i].ecc_range_mom[0] << " "
					<< std::setw(7) << std::setprecision(1) << ev.true_chains[i].ecc_range_mom[1] << " "
					<< std::setw(7) << std::setprecision(1) << ev.true_chains[i].ecc_mcs_mom[0] << " "
					<< std::setw(7) << std::setprecision(1) << ev.true_chains[i].ecc_mcs_mom[1] << " "
					<< std::setw(7) << std::setprecision(1) << ev.true_chains[i].bm_range_mom << " "
					<< std::setw(7) << std::setprecision(1) << ev.true_chains[i].bm_curvature_mom << " "
					<< std::setw(7) << std::setprecision(1) << ev.true_chains[i].ecc_range_mom_error[0][0] << " "
					<< std::setw(7) << std::setprecision(1) << ev.true_chains[i].ecc_range_mom_error[0][1] << " "
					<< std::setw(7) << std::setprecision(1) << ev.true_chains[i].ecc_range_mom_error[1][0] << " "
					<< std::setw(7) << std::setprecision(1) << ev.true_chains[i].ecc_range_mom_error[1][1] << " "
					<< std::setw(7) << std::setprecision(1) << ev.true_chains[i].ecc_mcs_mom_error[0][0] << " "
					<< std::setw(7) << std::setprecision(1) << ev.true_chains[i].ecc_mcs_mom_error[0][1] << " "
					<< std::setw(7) << std::setprecision(1) << ev.true_chains[i].ecc_mcs_mom_error[1][0] << " "
					<< std::setw(7) << std::setprecision(1) << ev.true_chains[i].ecc_mcs_mom_error[1][1] << " "
					<< std::setw(7) << std::setprecision(1) << ev.true_chains[i].bm_range_mom_error[0] << " "
					<< std::setw(7) << std::setprecision(1) << ev.true_chains[i].bm_range_mom_error[1] << " "
					<< std::setw(7) << std::setprecision(1) << ev.true_chains[i].bm_curvature_mom_error[0] << " "
					<< std::setw(7) << std::setprecision(1) << ev.true_chains[i].bm_curvature_mom_error[1] << " "
					<< std::setw(7) << std::setprecision(4) << ev.true_chains[i].muon_likelihood << " "
<< std::setw(7) << std::setprecision(4) << ev.true_chains[i].proton_likelihood << " "
<< std::setw(3) << std::setprecision(0) << ev.true_chains[i].base.size() << " "
<< std::setw(3) << std::setprecision(0) << ev.true_chains[i].base_pair.size() << std::endl;

for (auto &b : ev.true_chains[i].base) {
	ofs << std::right << std::fixed
		<< std::setw(4) << std::setprecision(0) << b.pl << " "
		<< std::setw(10) << std::setprecision(0) << b.rawid << " "
		<< std::setw(7) << std::setprecision(4) << b.ax << " "
		<< std::setw(7) << std::setprecision(4) << b.ay << " "
		<< std::setw(8) << std::setprecision(1) << b.x << " "
		<< std::setw(8) << std::setprecision(1) << b.y << " "
		<< std::setw(8) << std::setprecision(1) << b.z << " "
		<< std::setw(3) << std::setprecision(0) << b.m[0].zone << " "
		<< std::setw(4) << std::setprecision(0) << b.m[0].view << " "
		<< std::setw(3) << std::setprecision(0) << b.m[0].imager << " "
		<< std::setw(7) << std::setprecision(0) << b.m[0].ph << " "
		<< std::setw(5) << std::setprecision(0) << b.m[0].pixelnum << " "
		<< std::setw(5) << std::setprecision(0) << b.m[0].hitnum << " "
		<< std::setw(3) << std::setprecision(0) << b.m[1].zone << " "
		<< std::setw(4) << std::setprecision(0) << b.m[1].view << " "
		<< std::setw(3) << std::setprecision(0) << b.m[1].imager << " "
		<< std::setw(7) << std::setprecision(0) << b.m[1].ph << " "
		<< std::setw(5) << std::setprecision(0) << b.m[1].pixelnum << " "
		<< std::setw(5) << std::setprecision(0) << b.m[1].hitnum << std::endl;
}
for (auto &p : ev.true_chains[i].base_pair) {
	ofs << std::right << std::fixed
		<< std::setw(4) << std::setprecision(0) << p.first.pl << " "
		<< std::setw(10) << std::setprecision(0) << p.first.rawid << " "
		<< std::setw(4) << std::setprecision(0) << p.second.pl << " "
		<< std::setw(10) << std::setprecision(0) << p.second.rawid << " "
		<< std::setw(7) << std::setprecision(4) << p.first.ax << " "
		<< std::setw(7) << std::setprecision(4) << p.first.ay << " "
		<< std::setw(8) << std::setprecision(1) << p.first.x << " "
		<< std::setw(8) << std::setprecision(1) << p.first.y << " "
		<< std::setw(8) << std::setprecision(1) << p.first.z << " "
		<< std::setw(7) << std::setprecision(4) << p.second.ax << " "
		<< std::setw(7) << std::setprecision(4) << p.second.ay << " "
		<< std::setw(8) << std::setprecision(1) << p.second.x << " "
		<< std::setw(8) << std::setprecision(1) << p.second.y << " "
		<< std::setw(8) << std::setprecision(1) << p.second.z << std::endl;
}

			}

		}
		fprintf(stderr, "\r write group %d/%d(%4.1lf%%)\n", count, all, count*100. / all);

	}
	std::vector<Event_information> Read_Event_information_txt(std::string filename) {
		std::vector<Event_information> ret;
		std::ifstream ifs(filename);
		//filesize取得
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
		Mom_chain t;
		Mom_basetrack m;
		Event_information ev;
		std::pair<Mom_basetrack, Mom_basetrack> p;
		int base_num, base_pair_num;
		int chain_num, true_chain_num;
		while (ifs >> ev.groupid >> ev.unix_time >> ev.tracker_track_id >> ev.entry_in_daily_file >> ev.vertex_pl >> ev.ECC_id >> ev.vertex_material
			>> ev.vertex_position[0] >> ev.vertex_position[1] >> ev.vertex_position[2]
			>> ev.true_vertex_position[0] >> ev.true_vertex_position[1] >> ev.true_vertex_position[2]
			>> ev.nu_energy >> ev.nu_ax >> ev.nu_ay >> chain_num >> true_chain_num) {
			if (count % 100 == 0) {
				nowpos = ifs.tellg();
				auto size1 = nowpos - begpos;
				std::cerr << std::right << std::fixed << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%";
			}
			count++;
			for (int i = 0; i < chain_num; i++) {
				t.base.clear();
				t.base_pair.clear();
				t.base.reserve(base_num);
				t.base_pair.reserve(base_pair_num);
				ifs >> t.chainid >> t.stop_flg >> t.particle_flg >> t.direction >> t.charge_sign
					>> t.ecc_range_mom[0] >> t.ecc_range_mom[1] >> t.ecc_mcs_mom[0] >> t.ecc_mcs_mom[1]
					>> t.bm_range_mom >> t.bm_curvature_mom
					>> t.ecc_range_mom_error[0][0] >> t.ecc_range_mom_error[0][1] >> t.ecc_range_mom_error[1][0] >> t.ecc_range_mom_error[1][1]
					>> t.ecc_mcs_mom_error[0][0] >> t.ecc_mcs_mom_error[0][1] >> t.ecc_mcs_mom_error[1][0] >> t.ecc_mcs_mom_error[1][1]
					>> t.bm_range_mom_error[0] >> t.bm_range_mom_error[1]
					>> t.bm_curvature_mom_error[0] >> t.bm_curvature_mom_error[1]
					>> t.muon_likelihood >> t.proton_likelihood
					>> base_num >> base_pair_num;
				for (int i = 0; i < base_num; i++) {
					ifs >> m.pl >> m.rawid >> m.ax >> m.ay >> m.x >> m.y >> m.z
						>> m.m[0].zone >> m.m[0].view >> m.m[0].imager >> m.m[0].ph >> m.m[0].pixelnum >> m.m[0].hitnum
						>> m.m[1].zone >> m.m[1].view >> m.m[1].imager >> m.m[1].ph >> m.m[1].pixelnum >> m.m[1].hitnum;
					t.base.push_back(m);
				}
				for (int i = 0; i < base_pair_num; i++) {
					ifs >> p.first.pl >> p.first.rawid >> p.second.pl >> p.second.rawid
						>> p.first.ax >> p.first.ay >> p.first.x >> p.first.y >> p.first.z
						>> p.second.ax >> p.second.ay >> p.second.x >> p.second.y >> p.second.z;
					t.base_pair.push_back(p);
				}
				input_basetrack_information(t.base, t.base_pair);
				ev.chains.push_back(t);
			}
			for (int i = 0; i < true_chain_num; i++) {
				t.base.clear();
				t.base_pair.clear();
				t.base.reserve(base_num);
				t.base_pair.reserve(base_pair_num);
				ifs >> t.chainid >> t.stop_flg >> t.particle_flg >> t.direction >> t.charge_sign
					>> t.ecc_range_mom[0] >> t.ecc_range_mom[1] >> t.ecc_mcs_mom[0] >> t.ecc_mcs_mom[1]
					>> t.bm_range_mom >> t.bm_curvature_mom
					>> t.ecc_range_mom_error[0][0] >> t.ecc_range_mom_error[0][1] >> t.ecc_range_mom_error[1][0] >> t.ecc_range_mom_error[1][1]
					>> t.ecc_mcs_mom_error[0][0] >> t.ecc_mcs_mom_error[0][1] >> t.ecc_mcs_mom_error[1][0] >> t.ecc_mcs_mom_error[1][1]
					>> t.bm_range_mom_error[0] >> t.bm_range_mom_error[1]
					>> t.bm_curvature_mom_error[0] >> t.bm_curvature_mom_error[1]
					>> t.muon_likelihood >> t.proton_likelihood
					>> base_num >> base_pair_num;
				for (int i = 0; i < base_num; i++) {
					ifs >> m.pl >> m.rawid >> m.ax >> m.ay >> m.x >> m.y >> m.z
						>> m.m[0].zone >> m.m[0].view >> m.m[0].imager >> m.m[0].ph >> m.m[0].pixelnum >> m.m[0].hitnum
						>> m.m[1].zone >> m.m[1].view >> m.m[1].imager >> m.m[1].ph >> m.m[1].pixelnum >> m.m[1].hitnum;
					t.base.push_back(m);
				}
				for (int i = 0; i < base_pair_num; i++) {
					ifs >> p.first.pl >> p.first.rawid >> p.second.pl >> p.second.rawid
						>> p.first.ax >> p.first.ay >> p.first.x >> p.first.y >> p.first.z
						>> p.second.ax >> p.second.ay >> p.second.x >> p.second.y >> p.second.z;
					t.base_pair.push_back(p);
				}
				input_basetrack_information(t.base, t.base_pair);
				ev.true_chains.push_back(t);
			}
			ret.push_back(ev);

			ev.chains.clear();
			ev.true_chains.clear();
		}
		auto size1 = eofpos - begpos;
		std::cerr << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%" << std::endl;;
		if (count == 0) {
			fprintf(stderr, "%s no data!\n", filename.c_str());
			exit(1);
		}
		return ret;
	}

	std::vector<Event_information> Read_Event_information_bin_block(std::ifstream &ifs,int max_event_num) {
		std::vector<Event_information> ret;
		int64_t count = 0;
		Mom_chain t;
		Mom_basetrack m;
		std::pair<Mom_basetrack, Mom_basetrack> p;
		int base_num, base_pair_num;
		int chain_num, true_chain_num;
		int header[5];
		double mom_recon;
		Event_information ev;
		while (Read_Event_information_header(ifs, ev, chain_num, true_chain_num)) {
			for (int i_ch = 0; i_ch < chain_num; i_ch++) {
				if (Read_mom_chain_header(ifs, t, base_num, base_pair_num)) {
					t.base.clear();
					t.base_pair.clear();
					t.base.reserve(base_num);
					t.base_pair.reserve(base_pair_num);
					for (int j = 0; j < base_num; j++) {
						ifs.read((char*)& m, sizeof(Mom_basetrack));
						t.base.push_back(m);
					}
					for (int j = 0; j < base_pair_num; j++) {
						ifs.read((char*)& p.first, sizeof(Mom_basetrack));
						ifs.read((char*)& p.second, sizeof(Mom_basetrack));
						t.base_pair.push_back(p);
					}
					ev.chains.push_back(t);
				}
			}
			for (int i_ch = 0; i_ch < true_chain_num; i_ch++) {
				if (Read_mom_chain_header(ifs, t, base_num, base_pair_num)) {
					t.base.clear();
					t.base_pair.clear();
					t.base.reserve(base_num);
					t.base_pair.reserve(base_pair_num);
					for (int j = 0; j < base_num; j++) {
						ifs.read((char*)& m, sizeof(Mom_basetrack));
						t.base.push_back(m);
					}
					for (int j = 0; j < base_pair_num; j++) {
						ifs.read((char*)& p.first, sizeof(Mom_basetrack));
						ifs.read((char*)& p.second, sizeof(Mom_basetrack));
						t.base_pair.push_back(p);
					}
					ev.true_chains.push_back(t);
				}
			}
			ret.push_back(ev);
			ev.chains.clear();
			ev.true_chains.clear();
			count++;
			if (max_event_num == count) {
				break;
			}
		}

		std::cerr << "\r read Event_information " << ret.size() << std::endl;
		return ret;

	}
	void Write_Event_information_bin_block(std::ofstream &ofs, std::vector<Event_information>&ev_v) {
		if (ev_v.size() == 0) {
			fprintf(stderr, "target data ... null\n");
		}
		else {
			int64_t count = 0;
			int64_t max = ev_v.size();
			int base_num, base_pair_num;
			int chain_num, true_chain_num;

			for (int i = 0; i < max; i++) {
				Write_Event_information_header(ofs, ev_v[i]);
				chain_num = ev_v[i].chains.size();
				true_chain_num = ev_v[i].true_chains.size();
				for (int j = 0; j < chain_num; j++) {
					base_num = ev_v[i].chains[j].base.size();
					base_pair_num = ev_v[i].chains[j].base_pair.size();
					Write_mom_chain_header(ofs, ev_v[i].chains[j]);

					for (int k = 0; k < base_num; k++) {
						ofs.write((char*)& ev_v[i].chains[j].base[k], sizeof(Mom_basetrack));
					}
					for (int k = 0; k < base_pair_num; k++) {
						ofs.write((char*)& ev_v[i].chains[j].base_pair[k].first, sizeof(Mom_basetrack));
						ofs.write((char*)& ev_v[i].chains[j].base_pair[k].second, sizeof(Mom_basetrack));
					}
				}
				for (int j = 0; j < true_chain_num; j++) {
					base_num = ev_v[i].true_chains[j].base.size();
					base_pair_num = ev_v[i].true_chains[j].base_pair.size();
					Write_mom_chain_header(ofs, ev_v[i].true_chains[j]);

					for (int k = 0; k < base_num; k++) {
						ofs.write((char*)& ev_v[i].true_chains[j].base[k], sizeof(Mom_basetrack));
					}
					for (int k = 0; k < base_pair_num; k++) {
						ofs.write((char*)& ev_v[i].true_chains[j].base_pair[k].first, sizeof(Mom_basetrack));
						ofs.write((char*)& ev_v[i].true_chains[j].base_pair[k].second, sizeof(Mom_basetrack));
					}
				}
				count++;
			}

			std::cerr << "\r write Event_information " << count << std::endl;
		}
	}

	void input_basetrack_information(std::vector<Mom_basetrack> &base,std::vector<std::pair<Mom_basetrack, Mom_basetrack>> &base_pair) {

		std::multimap<std::pair<int,int>, Mom_basetrack*> base_pair_p;
		for (auto itr = base_pair.begin(); itr != base_pair.end(); itr++) {
			base_pair_p.insert(std::make_pair(std::make_pair(itr->first.pl, itr->first.rawid), &(itr->first)));
			base_pair_p.insert(std::make_pair(std::make_pair(itr->second.pl, itr->second.rawid), &(itr->second)));
		}
		for (auto itr = base.begin(); itr != base.end(); itr++) {
			std::pair<int, int> id = std::make_pair(itr->pl, itr->rawid);
			if (base_pair_p.count(id) == 0)continue;
			auto range = base_pair_p.equal_range(id);
			for (auto res = range.first; res != range.second; res++) {
				res->second->m[0] = itr->m[0];
				res->second->m[1] = itr->m[1];
			}
		}


	}

	std::vector<Event_information> Read_Event_information_extension(std::string filename) {
		std::vector<Event_information> ret;
		std::string extension;
		extension = filename.substr(filename.size() - 5, 5);
		if (extension == "momch") {
			ret = Read_Event_information_bin(filename);
		}
		else {
			ret = Read_Event_information_txt(filename);
		}
		return ret;
	}
	void Write_Event_information_extension(std::string filename, std::vector<Event_information>&ev_v) {
		std::string extension;
		extension = filename.substr(filename.size() - 5, 5);
		if (extension == "momch") {
			Write_Event_information_bin(filename, ev_v);
		}
		else {
			Write_Event_information_txt(filename, ev_v);
		}
	}


	microtrack_minimum::microtrack_minimum() {
		ph = -1;
		zone = -1;
		imager = -1;
		pixelnum = -1;
		hitnum = -1;
	}
	Mom_basetrack::Mom_basetrack() {
		pl = -1;
		rawid = -1;
		ax = NAN;
		ay = NAN;
		x = NAN;
		y = NAN;
		z = NAN;
	}
	Mom_chain::Mom_chain() {
		chainid = -1;
		stop_flg = -1;
		particle_flg = -1;
		direction = 0;
		charge_sign = 0;
		ecc_range_mom[0] = -1;
		ecc_range_mom[1] = -1;
		ecc_mcs_mom[0] = -1;
		ecc_mcs_mom[1] = -1;
		bm_range_mom = -1;
		charge_sign = -1;
		ecc_range_mom_error[0][0] = -1;
		ecc_range_mom_error[0][1] = -1;
		ecc_range_mom_error[1][0] = -1;
		ecc_range_mom_error[1][1] = -1;
		ecc_mcs_mom_error[0][0] = -1;
		ecc_mcs_mom_error[0][1] = -1;
		ecc_mcs_mom_error[1][0] = -1;
		ecc_mcs_mom_error[1][1] = -1;
		bm_range_mom_error[0] = -1;
		bm_range_mom_error[1] = -1;
		bm_curvature_mom_error[0] = -1;
		bm_curvature_mom_error[1] = -1;
		muon_likelihood = NAN;
		proton_likelihood = NAN;
	}
	Event_information::Event_information() {
		groupid = -1;
		unix_time = -1;
		tracker_track_id = -1;
		entry_in_daily_file = -1;
		vertex_pl = -1;
		ECC_id = -1;
		vertex_material = -1;
		vertex_position[0] = NAN;
		vertex_position[1] = NAN;
		vertex_position[2] = NAN;
		true_vertex_position[0] = NAN;
		true_vertex_position[1] = NAN;
		true_vertex_position[2] = NAN;
		weight = -1;
		nu_energy = -1;
		nu_ax = NAN;
		nu_ay = NAN;
	}

	double Mom_chain::Get_muon_mcs_pb() {
		if (!isfinite(ecc_mcs_mom[0]))return -1;
		if (ecc_mcs_mom[0] < 0)return -1;
		else if (ecc_mcs_mom[0] < 0.0000001)return 0;

		double mass = 105.65836668;
		double beta = 1 / sqrt(1 + pow(mass / ecc_mcs_mom[0], 2));
		return ecc_mcs_mom[0] * beta;
	}
	void Mom_chain::Get_muon_pb_mcs_error(double *pb_error) {
		double mass = 105.65836668;
		double pb_edge,beta;

		if (!isfinite(ecc_mcs_mom[0])||ecc_mcs_mom[0] < 0) {
			pb_error[0] = -1;
			pb_error[1] = -1;
			return;
		}
		//下限の計算
		if (ecc_mcs_mom_error[0][0] < 0)pb_error[0] = -1;
		else {
			pb_edge = ecc_mcs_mom[0] - ecc_mcs_mom_error[0][0];
			if (pb_edge < 0)pb_error[0] = -1;
			else if (pb_edge < 0.0000001) {
				pb_error[0] = Get_muon_mcs_pb();
			}
			else {
				beta = 1 / sqrt(1 + pow(mass / pb_edge, 2));
				pb_error[0] = Get_muon_mcs_pb() - pb_edge * beta;
			}
		}
		//上限の計算
		if (ecc_mcs_mom_error[0][1] < 0)pb_error[1] = -1;
		else {
			pb_edge = ecc_mcs_mom[0] + ecc_mcs_mom_error[0][1];
			if (pb_edge < 0)pb_error[1] = -1;
			else if (pb_edge < 0.0000001) {
				pb_error[1] = Get_muon_mcs_pb();
			}
			else {
				beta = 1 / sqrt(1 + pow(mass / pb_edge, 2));
				pb_error[1] = pb_edge * beta- Get_muon_mcs_pb();
			}
		}
	}
	double Mom_chain::Get_proton_mcs_pb() {
		double mass = 938.2720813;
		if (!isfinite(ecc_mcs_mom[1]))return -1;
		if (ecc_mcs_mom[1] < 0)return -1;
		else if (ecc_mcs_mom[1] < 0.0000001)return 0;

		double beta = 1 / sqrt(1 + pow(mass / ecc_mcs_mom[1], 2));
		return ecc_mcs_mom[1] * beta;

	}

	void Mom_chain::Get_proton_pb_mcs_error(double *pb_error) {
		double mass = 938.2720813;
		double pb_edge, beta;

		if (!isfinite(ecc_mcs_mom[1])||ecc_mcs_mom[1] < 0) {
			pb_error[0] = -1;
			pb_error[1] = -1;
			return;
		}
		//下限の計算
		if (ecc_mcs_mom_error[1][0] < 0)pb_error[0] = -1;
		else {
			pb_edge = ecc_mcs_mom[1] - ecc_mcs_mom_error[1][0];
			if (pb_edge < 0)pb_error[0] = -1;
			else if (pb_edge < 0.0000001) {
				pb_error[0] = Get_proton_mcs_pb();
			}
			else {
				beta = 1 / sqrt(1 + pow(mass / pb_edge, 2));
				pb_error[0] = Get_proton_mcs_pb() - pb_edge * beta;
			}
		}
		//上限の計算
		if (ecc_mcs_mom_error[1][1] < 0)pb_error[1] = -1;
		else {
			pb_edge = ecc_mcs_mom[1] + ecc_mcs_mom_error[1][1];
			if (pb_edge < 0)pb_error[1] = -1;
			else if (pb_edge < 0.0000001) {
				pb_error[1] = Get_proton_mcs_pb();
			}
			else {
				beta = 1 / sqrt(1 + pow(mass / pb_edge, 2));
				pb_error[1] = pb_edge * beta - Get_proton_mcs_pb();
			}
		}
	}

}

int main(int argc,char**argv) {
	if (argc != 3) {
		fprintf(stderr, "usage\n");
		exit(1);
	}

	std::string file_in_momch = argv[1];
	std::string file_out_momch = argv[2];

	std::vector<Momentum_recon::Event_information> ev = Momentum_recon::Read_Event_information_extension(file_in_momch);

	Momentum_recon::Write_Event_information_extension(file_out_momch, ev);
}