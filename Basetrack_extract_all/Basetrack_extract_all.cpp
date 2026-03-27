#define _CRT_SECURE_NO_WARNINGS

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
#include <FILE_structure.h>

#include <filesystem>

class microtrack_minimum {
public:
	int ph, zone, view, imager, pixelnum, hitnum;
};
class basetrack_minimum {
public:
	int pl, rawid;
	float ax, ay, x, y;
	microtrack_minimum m[2];
};

class microtrack_inf {
public:
	int zone, pos, col, row, isg, ph, pixelnum, hitnum;
};

void Pick_up_pixel_count(std::vector<std::pair<vxx::base_track_t, std::pair<microtrack_minimum, microtrack_minimum>>>& base, std::string file_in_ECC);

void add_microtrack_inf_format0(std::string filename, std::multimap<std::tuple<int, int, int, int, int>, int*>& base_px_hit, std::multimap<std::tuple<int, int, int, int, int>, int*>& base_px_area, std::set<int>& zone_pos1, std::set<int>& zone_pos2);
void add_microtrack_inf_format1(std::string filename, std::multimap<std::tuple<int, int, int, int, int>, int*>& base_px_hit, std::multimap<std::tuple<int, int, int, int, int>, int*>& base_px_area);

std::vector< basetrack_minimum>change_format(std::vector<std::pair<vxx::base_track_t, std::pair<microtrack_minimum, microtrack_minimum>>>& base);
void out_base_inf(std::ofstream& ofs, std::vector< basetrack_minimum>& b_out);
std::vector<std::pair<vxx::base_track_t, std::pair<microtrack_minimum, microtrack_minimum>>> read_base(std::string filename, int pl);

int main(int argc, char** argv) {
	if (argc != 4) {
		fprintf(stderr, "usage:file-in-ECC file-in-microtrack file-out\n");
		exit(1);
	}

	std::string file_in_ECC = argv[1];//marged
	std::string file_basetrack = argv[2];//each area
	std::string file_out = argv[3];
	std::ofstream ofs(file_out, std::ios::binary);

	for (int pl = 1; pl <= 133; pl++) {
		printf("PL%03d start\n", pl);

		std::stringstream file_in_base;
		file_in_base << file_in_ECC << "\\Area0\\PL"
			<< std::setw(3) << std::setfill('0') << pl << "\\b"
			<< std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";
		printf("input file [%s]\n", file_in_base.str().c_str());
		if (!std::filesystem::exists(file_in_base.str()))continue;

		std::vector<std::pair<vxx::base_track_t, std::pair<microtrack_minimum, microtrack_minimum>>> base_w_micro = read_base(file_in_base.str(), pl);
		Pick_up_pixel_count(base_w_micro, file_basetrack);
		std::vector< basetrack_minimum> b_out = change_format(base_w_micro);
		out_base_inf(ofs, b_out);

	}
}
std::vector<std::pair<vxx::base_track_t, std::pair<microtrack_minimum, microtrack_minimum>>> read_base(std::string filename, int pl) {
	vxx::BvxxReader br;

	std::vector<vxx::base_track_t> base = br.ReadAll(filename, pl, 0);
	std::vector<std::pair<vxx::base_track_t, std::pair<microtrack_minimum, microtrack_minimum>>> ret;
	const int NumberOfImager = 72;
	uint32_t ShotID;
	int view, imager;
	for (auto& b : base) {
		std::pair<microtrack_minimum, microtrack_minimum> micro_pair;

		ShotID = ((uint32_t)(uint16_t)b.m[0].row << 16) | ((uint32_t)(uint16_t)b.m[0].col);
		view = ShotID / NumberOfImager;
		imager = ShotID % NumberOfImager;

		micro_pair.first.imager = imager;
		micro_pair.first.ph = b.m[0].ph;
		micro_pair.first.view = view;
		micro_pair.first.zone = b.m[0].zone;
		micro_pair.first.hitnum = -1;
		micro_pair.first.pixelnum = -1;

		ShotID = ((uint32_t)(uint16_t)b.m[1].row << 16) | ((uint32_t)(uint16_t)b.m[1].col);
		view = ShotID / NumberOfImager;
		imager = ShotID % NumberOfImager;

		micro_pair.second.imager = imager;
		micro_pair.second.ph = b.m[1].ph;
		micro_pair.second.view = view;
		micro_pair.second.zone = b.m[1].zone;
		micro_pair.second.hitnum = -1;
		micro_pair.second.pixelnum = -1;

		ret.push_back(std::make_pair(b, micro_pair));
	}
	return ret;
}

std::vector< basetrack_minimum>change_format(std::vector<std::pair<vxx::base_track_t, std::pair<microtrack_minimum, microtrack_minimum>>>& base) {
	std::vector< basetrack_minimum> ret;

	for (auto& b : base) {
		basetrack_minimum b_tmp;
		b_tmp.ax = b.first.ax;
		b_tmp.ay = b.first.ay;
		b_tmp.x = b.first.x;
		b_tmp.y = b.first.y;
		b_tmp.pl = b.first.pl;
		b_tmp.rawid = b.first.rawid;
		b_tmp.m[0] = b.second.first;
		b_tmp.m[1] = b.second.second;
		ret.push_back(b_tmp);
	}
	return ret;
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


//高速化の実装-->OK
//void Pick_up_pixel_count(std::vector<output_format_link>&link, std::string file_in_ECC) {
void Pick_up_pixel_count(std::vector<std::pair<vxx::base_track_t, std::pair<microtrack_minimum, microtrack_minimum>>>& base, std::string file_in_ECC) {
	if (base.size() == 0)return;
	int pl = base.begin()->first.pl;

	const int NumberOfImager = 72;
	uint32_t ShotID;
	int view, imager;

	std::tuple<int, int, int, int, int> id;
	//linkletのpxへのポインタをmapのvalにする
	std::multimap<std::tuple<int, int, int, int, int>, int*> base_px_hit;
	std::multimap<std::tuple<int, int, int, int, int>, int*> base_px_area;
	for (auto itr = base.begin(); itr != base.end(); itr++) {
		std::get<0>(id) = itr->first.m[0].pos;
		std::get<1>(id) = itr->second.first.zone;
		std::get<2>(id) = itr->second.first.view;
		std::get<3>(id) = itr->second.first.imager;
		std::get<4>(id) = itr->first.m[0].isg;

		base_px_hit.insert(std::make_pair(id, &(itr->second.first.hitnum)));
		base_px_area.insert(std::make_pair(id, &(itr->second.first.pixelnum)));

		std::get<0>(id) = itr->first.m[1].pos;
		std::get<1>(id) = itr->second.second.zone;
		std::get<2>(id) = itr->second.second.view;
		std::get<3>(id) = itr->second.second.imager;
		std::get<4>(id) = itr->first.m[1].isg;

		base_px_hit.insert(std::make_pair(id, &(itr->second.second.hitnum)));
		base_px_area.insert(std::make_pair(id, &(itr->second.second.pixelnum)));
	}


	std::set<int>zone_pos1, zone_pos2;
	for (int Area = 1; Area <= 6; Area++) {
		for (int iscan = 0; iscan < 4; iscan++) {
			printf("PL%03d Area%d scan %d", pl, Area, iscan);
			//読み込むfile名
			std::stringstream infile_thick, infile_thin;
			infile_thick << file_in_ECC << "\\Area"
				<< std::setw(1) << Area << "\\PL"
				<< std::setw(3) << std::setfill('0') << pl << "\\micro_inf_thick_"
				<< std::setw(1) << iscan;
			infile_thin << file_in_ECC << "\\Area"
				<< std::setw(1) << Area << "\\PL"
				<< std::setw(3) << std::setfill('0') << pl << "\\micro_inf_thin_"
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
			add_microtrack_inf_format0(infile_thick.str(), base_px_hit, base_px_area, zone_pos1, zone_pos2);
			add_microtrack_inf_format0(infile_thin.str(), base_px_hit, base_px_area, zone_pos1, zone_pos2);
		}
	}
	if (zone_pos1.size() == 48 && zone_pos2.size() == 48) {
		return;
	}
	else if (zone_pos1.size() == 24 && zone_pos2.size() == 24) {
		for (int Area = 1; Area <= 6; Area++) {
			for (int iscan = 0; iscan < 4; iscan++) {
				printf("PL%03d Area%d scan %d", pl, Area, iscan);
				//読み込むfile名
				std::stringstream infile_thick, infile_thin;
				infile_thick << file_in_ECC << "\\Area"
					<< std::setw(1) << Area << "\\PL"
					<< std::setw(3) << std::setfill('0') << pl << "\\micro_inf_thick_"
					<< std::setw(1) << iscan;
				infile_thin << file_in_ECC << "\\Area"
					<< std::setw(1) << Area << "\\PL"
					<< std::setw(3) << std::setfill('0') << pl << "\\micro_inf_thin_"
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
				add_microtrack_inf_format1(infile_thick.str(), base_px_hit, base_px_area);
				add_microtrack_inf_format1(infile_thin.str(), base_px_hit, base_px_area);
			}
		}

	}
	else {
		fprintf(stderr, "file format exception\n");
		fprintf(stderr, "zone1 size = %d zone2 size=%d\n", zone_pos1.size(), zone_pos2.size());
		//exit(1);
	}
}

void add_microtrack_inf_format0(std::string filename, std::multimap<std::tuple<int, int, int, int, int>, int*>& base_px_hit, std::multimap<std::tuple<int, int, int, int, int>, int*>& base_px_area, std::set<int>& zone_pos1, std::set<int>& zone_pos2) {
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
			if (base_px_hit.count(id) != 0) {
				auto range = base_px_hit.equal_range(id);
				for (auto res = range.first; res != range.second; res++) {
					(*res->second) = m_buf[i].hitnum;
				}
			}
			if (base_px_area.count(id) != 0) {
				auto range = base_px_area.equal_range(id);
				for (auto res = range.first; res != range.second; res++) {
					(*res->second) = m_buf[i].pixelnum;
				}
			}
		}
		now += read_num;
	}
	fclose(fp_in);
}

void add_microtrack_inf_format1(std::string filename, std::multimap<std::tuple<int, int, int, int, int>, int*>& base_px_hit, std::multimap<std::tuple<int, int, int, int, int>, int*>& base_px_area) {
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
			if (base_px_hit.count(id) != 0) {
				auto range = base_px_hit.equal_range(id);
				for (auto res = range.first; res != range.second; res++) {
					(*res->second) = m_buf[i].hitnum;
				}
			}
			if (base_px_area.count(id) != 0) {
				auto range = base_px_area.equal_range(id);
				for (auto res = range.first; res != range.second; res++) {
					(*res->second) = m_buf[i].pixelnum;
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
			if (base_px_hit.count(id) != 0) {
				auto range = base_px_hit.equal_range(id);
				for (auto res = range.first; res != range.second; res++) {
					(*res->second) = m_buf[i].hitnum;
				}
			}
			if (base_px_area.count(id) != 0) {
				auto range = base_px_area.equal_range(id);
				for (auto res = range.first; res != range.second; res++) {
					(*res->second) = m_buf[i].pixelnum;
				}
			}

		}
		now += read_num;
	}
	fclose(fp_in);

}


void out_base_inf(std::ofstream& ofs, std::vector< basetrack_minimum>& b_out) {
	if (b_out.size() == 0)return;
	int64_t num = b_out.size();
	int pl = b_out.begin()->pl;
	ofs.write((char*)&pl, sizeof(pl));
	ofs.write((char*)&num, sizeof(num));

	int64_t count = 0;
	for (int64_t i = 0; i < num; i++) {
		if (count % 100000 == 0) {
			fprintf(stderr, "\r write file %15lld/%15lld(%4.1lf%%)", count, num, count * 100. / num);
		}
		count++;

		ofs.write((char*)&b_out[i], sizeof(basetrack_minimum));

	}
	fprintf(stderr, "\r write file %15lld/%15lld(%4.1lf%%)\n", count, num, count * 100. / num);



}
