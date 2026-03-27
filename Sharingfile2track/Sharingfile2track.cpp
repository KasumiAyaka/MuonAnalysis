#define _CRT_SECURE_NO_WARNINGS
#pragma comment(lib, "VxxReader.lib")
#include "VxxReader.h"
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

#include <filesystem>
#include <set>

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

//////////////////
void remove_nomuon_track(std::vector<Sharing_file::Sharing_file>& sf);
std::vector < std::pair<Sharing_file::Sharing_file, std::pair<int, int>>> pl4_track_multidelete(std::vector<Sharing_file::Sharing_file>& sf, std::string file_in_ECC, std::vector<output_format_link>& link);
void output_Track(std::string filename, std::vector < std::pair<Sharing_file::Sharing_file, std::pair<int, int>>>& sf_baseid);
void output_sharingfile(std::string filename, std::vector < std::pair<Sharing_file::Sharing_file, std::pair<int, int>>>& sf_baseid);
std::vector<output_format_link> read_link_bin(std::string filename);

int main(int argc, char** argv) {
	if (argc != 7) {
		fprintf(stderr, "usage:prg sharing-file file-mode ECC-Area-path output-track-file link-sel-flg output-sharing-file\n");
		exit(1);
	}
	std::string file_in_sharing = argv[1];
	int file_mode = std::stoi(argv[2]);
	std::string file_in_ECC = argv[3];
	std::string file_out_track = argv[4];
	int link_sel_mode = std::stoi(argv[5]);
	std::string file_out_gsf = argv[6];

	//sharing fileÇÃì«Ç›çûÇ›
	std::vector<Sharing_file::Sharing_file> sf;
	if (file_mode == 0)sf = Sharing_file::Read_sharing_file_bin(file_in_sharing);
	else if (file_mode == 1)sf = Sharing_file::Read_sharing_file_txt(file_in_sharing);
	else {
		fprintf(stderr, "filemode exception\n");
		fprintf(stderr, "mode=0 binary sharing file\n");
		fprintf(stderr, "mode=1 txt    sharing file\n");
		exit(1);
	}

	//linkletì«Ç›çûÇ›
	std::stringstream file_in_link;
	if (link_sel_mode == 1) {
		file_in_link << file_in_ECC << "\\0\\linklet_3d\\link_003_004.bin";
	}
	else if (link_sel_mode == 0) {
		file_in_link << file_in_ECC << "\\0\\linklet_3d\\link_003_004.sel.bin";
	}
	else {
		fprintf(stderr, "link sel flg exception\n");
		fprintf(stderr, "flg=0 use no selection linklet\n");
		fprintf(stderr, "flg=1 use    selected  linklet\n");
		exit(1);
	}
	std::vector<output_format_link> link = read_link_bin(file_in_link.str());

	remove_nomuon_track(sf);
	//PL003-PL004 linkletÇÃPL004ÇÃçÌèú+basetracak rawidïtó^
	std::vector < std::pair<Sharing_file::Sharing_file, std::pair<int, int>>> sf_baseid = pl4_track_multidelete(sf, file_in_ECC, link);

	output_Track(file_out_track, sf_baseid);
	output_sharingfile(file_out_gsf, sf_baseid);

}
void remove_nomuon_track(std::vector<Sharing_file::Sharing_file>& sf) {
	std::vector<Sharing_file::Sharing_file> ret;
	//eventidÇÃÇ¬Ç¢ÇƒÇ¢Ç»Ç¢trackÇÃçÌèú
	printf("Shifter connected track --> %d\n", sf.size());

	for (auto itr = sf.begin(); itr != sf.end(); itr++) {
		if (itr->eventid < 0) continue;
		ret.push_back(*itr);
	}
	std::swap(ret, sf);
	printf("Shifter connected muon --> %d\n", sf.size());

}

std::vector < std::pair<Sharing_file::Sharing_file, std::pair<int, int>>> pl4_track_multidelete(std::vector<Sharing_file::Sharing_file>& sf, std::string file_in_ECC, std::vector<output_format_link>& link) {

	std::vector<vxx::base_track_t> base_pl3, base_pl4;
	std::string file_in_base_pl3 = file_in_ECC + "\\PL003\\b003.sel.cor.vxx";
	std::string file_in_base_pl4 = file_in_ECC + "\\PL004\\b004.sel.cor.vxx";
	vxx::BvxxReader br;
	base_pl3 = br.ReadAll(file_in_base_pl3, 3, 0);
	base_pl4 = br.ReadAll(file_in_base_pl4, 4, 0);

	std::ofstream ofs("Sharingfile2track_log.txt");

	std::map<std::tuple<int, int, int, int>, std::pair<int, int>> base_pl3_id, base_pl4_id;
	for (auto& b : base_pl3) {
		base_pl3_id.insert(std::make_pair(std::make_tuple(b.m[0].zone, b.m[0].rawid, b.m[1].zone, b.m[1].rawid), std::make_pair(b.pl, b.rawid)));
	}
	for (auto& b : base_pl4) {
		base_pl4_id.insert(std::make_pair(std::make_tuple(b.m[0].zone, b.m[0].rawid, b.m[1].zone, b.m[1].rawid), std::make_pair(b.pl, b.rawid)));
	}


	std::multimap<std::pair<int, int>, std::pair<int, int>>linklet_id;
	for (auto& l : link) {
		linklet_id.insert(std::make_pair(std::make_pair(l.b[1].pl, l.b[1].rawid), std::make_pair(l.b[0].pl, l.b[0].rawid)));
	}

	std::vector < std::pair<Sharing_file::Sharing_file, std::pair<int, int>>>ret;
	for (auto itr = sf.begin(); itr != sf.end(); itr++) {
		if (itr->pl == 3) {
			auto res = base_pl3_id.find(std::make_tuple(itr->zone[0], itr->rawid[0], itr->zone[1], itr->rawid[1]));
			if (res == base_pl3_id.end()) {
				fprintf(stderr, "PL%03d %d %d %d %d not found\n", itr->pl, itr->zone[0], itr->rawid[0], itr->zone[1], itr->rawid[1]);
				ofs << "PL00" << itr->pl << " " << itr->zone[0] << " " << itr->rawid[0] << " " << itr->zone[1] << " " << itr->rawid[1] << std::endl;
				//exit(1);
			}
			else {
				ret.push_back(std::make_pair(*itr, res->second));
			}
		}
		else {
			auto res = base_pl4_id.find(std::make_tuple(itr->zone[0], itr->rawid[0], itr->zone[1], itr->rawid[1]));
			if (res == base_pl4_id.end()) {
				fprintf(stderr, "PL%03d %d %d %d %d not found\n", itr->pl, itr->zone[0], itr->rawid[0], itr->zone[1], itr->rawid[1]);
				ofs << "PL00" << itr->pl << " " << itr->zone[0] << " " << itr->rawid[0] << " " << itr->zone[1] << " " << itr->rawid[1] << std::endl;
				//exit(1);
			}
			if (res != base_pl4_id.end()&&linklet_id.count(res->second) == 0) {
				ret.push_back(std::make_pair(*itr, res->second));
				continue;
			}

			auto range = linklet_id.equal_range(res->second);
			bool flg = false;
			for (auto res2 = range.first; res2 != range.second; res2++) {
				Sharing_file::Sharing_file sf_tmp;
				bool sf_detect = false;
				for (auto itr2 = sf.begin(); itr2 != sf.end(); itr2++) {
					if (itr2->pl != 3)continue;
					auto res_pl3 = base_pl3_id.find(std::make_tuple(itr2->zone[0], itr2->rawid[0], itr2->zone[1], itr2->rawid[1]));
					if (res_pl3 == base_pl3_id.end()) {
						fprintf(stderr, "PL%03d %d %d %d %d not found\n", itr2->pl, itr2->zone[0], itr2->rawid[0], itr2->zone[1], itr->rawid[1]);
						//ofs << "PL00" << itr->pl << " " << itr->zone[0] << " " << itr->rawid[0] << " " << itr->zone[1] << " " << itr->rawid[1] << std::endl;
						//exit(1);
					}
					if (res_pl3 != base_pl3_id.end()&&res_pl3->second == res2->second) {
						sf_tmp = *itr2;
						sf_detect = true;
						break;
					}
				}
				if (sf_detect) {
					if (sf_tmp.unix_time != itr->unix_time)continue;
					if (sf_tmp.tracker_track_id != itr->tracker_track_id)continue;
					if (sf_tmp.fixedwall_id != itr->fixedwall_id)continue;
					if (sf_tmp.trackerwall_id != itr->trackerwall_id)continue;
					if (sf_tmp.oss_id != itr->oss_id)continue;
					flg = true;
					break;
				}
			}
			if (flg)continue;
			ret.push_back(std::make_pair(*itr, res->second));
		}
	}

	fprintf(stderr, "PL003-PL004 linklet delete %d-->%d\n", sf.size(), ret.size());
	return ret;


}
void output_Track(std::string filename, std::vector < std::pair<Sharing_file::Sharing_file, std::pair<int, int>>>& sf_baseid) {
	std::ofstream ofs(filename);
	int trackid = 0;
	for (auto itr = sf_baseid.begin(); itr != sf_baseid.end(); itr++) {
		ofs << std::fixed << std::right
			<< std::setw(10) << std::setprecision(0) << itr->first.eventid << " "
			<< std::setw(5) << std::setprecision(0) << trackid << " "
			<< std::setw(4) << std::setprecision(0) << itr->second.first << " "
			<< std::setw(10) << std::setprecision(0) << itr->second.second << std::endl;
	}

}

std::vector<output_format_link> read_link_bin(std::string filename) {
	std::vector<output_format_link> link;
	std::ifstream ifs(filename, std::ios::binary);
	//filesizeéÊìæ
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
	output_format_link l;
	while (ifs.read((char*)&l, sizeof(output_format_link))) {
		if (count % 10000 == 0) {
			nowpos = ifs.tellg();
			auto size1 = nowpos - begpos;
			std::cerr << std::right << std::fixed << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%";
		}
		count++;

		link.emplace_back(l);
	}
	auto size1 = eofpos - begpos;
	std::cerr << "\r now reading ..." << std::setw(4) << std::setprecision(1) << size1 * 100. / size2 << "%" << std::endl;;
	if (count == 0) {
		fprintf(stderr, "%s no linklet!\n", filename.c_str());
		exit(1);
	}
	return link;
}
void output_sharingfile(std::string filename, std::vector < std::pair<Sharing_file::Sharing_file, std::pair<int, int>>>& sf_baseid) {
	std::vector<Sharing_file::Sharing_file> sf;
	for (auto itr = sf_baseid.begin(); itr != sf_baseid.end(); itr++) {
		sf.push_back(itr->first);
	}

	Sharing_file::Write_sharing_file_txt(filename, sf);
}