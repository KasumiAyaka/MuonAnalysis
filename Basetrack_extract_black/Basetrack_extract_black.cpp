#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <set>
#include <map>

#pragma comment(lib, "VxxReader.lib")
#include "VxxReader.h"
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

//#include <boost/unordered_set.hpp>
//#include <boost/unordered_map.hpp>

class linklet_header {
public:
	int pos0, raw0, pos1, raw1;
};
class basetrack_minimum {
public:
	int pl, rawid, ph;
	float ax, ay, x, y;
};
void read_linklet_list(std::string filename, std::set<std::pair<int, int>>& base_id);
void write_base_min(std::string filename, std::vector<basetrack_minimum>& base_min);

void read_vph_cut_param(std::string filename, int pl_in, std::map<double, std::pair<int, int>>& ret);
std::vector<vxx::base_track_t> basetrack_vph_cut(std::vector<vxx::base_track_t>& base, std::map<double, std::pair<int, int>> thr);



int main(int argc, char** argv) {
	if (argc != 4) {
		fprintf(stderr, "usage:prg file_ECC file_cut_param file_out_base\n");
		exit(1);
	}
	std::string file_in_ECC = argv[1];
	std::string file_in_black_param = argv[2];
	std::string file_out_base = argv[3];


	std::vector<basetrack_minimum>base_min;



	vxx::BvxxReader br;
	for (int pl = 3; pl <= 133; pl++) {
		if (pl == 15) {
			printf("PL%03d was skiped.In case of ECC6", pl);
			continue;
		}
		printf("PL%03d basetrack read\n", pl);

		std::map<double, std::pair<int, int>>vph_cut;
		read_vph_cut_param(file_in_black_param, pl, vph_cut);

		std::stringstream file_in_base;
		file_in_base << file_in_ECC << "\\Area0\\PL"
			<< std::setw(3) << std::setfill('0') << pl << "\\b"
			<< std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";
		std::vector<vxx::base_track_t >base = br.ReadAll(file_in_base.str(), pl, 0);
		base = basetrack_vph_cut(base, vph_cut);

		for (auto itr_b = base.begin(); itr_b != base.end(); itr_b++) {
			basetrack_minimum b;
			b.ax = itr_b->ax;
			b.ay = itr_b->ay;
			b.x = itr_b->x;
			b.y = itr_b->y;
			b.pl = itr_b->pl;
			b.ph = itr_b->m[0].ph + itr_b->m[1].ph;
			b.rawid = itr_b->rawid;
			base_min.push_back(b);
		}
	}

	printf("all Blacck basetrack %lld\n", base_min.size());

	write_base_min(file_out_base, base_min);


}
void read_vph_cut_param(std::string filename, int pl_in, std::map<double, std::pair<int, int>>& ret) {
	std::ifstream ifs(filename.c_str());
	int pl;
	double angle_min, angle_max, vph_thr[2];

	//std::map<double, std::pair<int, int>> ret[2];

	while (ifs >> pl >> angle_min >> angle_max >> vph_thr[0] >> vph_thr[1]) {
		if (pl == pl_in) {
			ret.insert(std::make_pair(angle_max, std::make_pair(int(vph_thr[0]), int(vph_thr[1]))));
		}
	}
	if (ret.size() == 0) {
		fprintf(stderr, "PL%03d vph threshold not found\n", pl_in);
		exit(1);
	}
}
std::vector<vxx::base_track_t> basetrack_vph_cut(std::vector<vxx::base_track_t>& base, std::map<double, std::pair<int, int>> thr) {
	int64_t all = base.size();
	std::vector<vxx::base_track_t> ret;
	double angle;
	for (auto itr = base.begin(); itr != base.end(); itr++) {
		angle = sqrt(itr->ax * itr->ax + itr->ay * itr->ay);
		auto res = thr.upper_bound(angle);
		if (res == thr.end()) {
			res = std::prev(res, 1);
		}
		if (itr->m[0].ph % 10000 <= res->second.first)continue;
		if (itr->m[1].ph % 10000 <= res->second.second)continue;
		ret.push_back(*itr);

	}
	printf("%lld -->%lld(%4.1lf%%)\n", all, ret.size(), ret.size() * 100. / all);
	return ret;

}

void write_base_min(std::string filename, std::vector<basetrack_minimum>& base_min) {
	std::ofstream ofs(filename, std::ios::binary);
	if (!ofs) {
		//file open s
		fprintf(stderr, "File[%s] is not exist!!\n", filename.c_str());
		exit(1);
	}
	if (base_min.size() == 0) {
		fprintf(stderr, "target basetrack minimum ... null\n");
		fprintf(stderr, "File[%s] has no text\n", filename.c_str());
	}
	int64_t count = 0;
	int64_t max = base_min.size();
	for (int i = 0; i < base_min.size(); i++) {
		if (count % 10000 == 0) {
			std::cerr << std::right << std::fixed << "\r now writing ..." << std::setw(4) << std::setprecision(1) << count * 100. / max << "%";
		}
		count++;
		ofs.write((char*)&base_min[i], sizeof(basetrack_minimum));
	}
	std::cerr << std::right << std::fixed << "\r now writing ..." << std::setw(4) << std::setprecision(1) << count * 100. / max << "%" << std::endl;
	ofs.close();

}