//muonのchainに対して、linkletから再度chainを再生成する
//1group 1chainが前提
#define _CRT_SECURE_NO_WARNINGS
//大体2^50
#define DEFAULT_CHAIN_UPPERLIM 1000000000000000

//#pragma comment(lib, "VxxReader.lib")
//#include "VxxReader.h"
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>
#pragma comment(lib, "LibL2c-x.lib")
#include <LibL2c-x.h>


#include <filesystem>
#include <set>
#include <netscan_data_types_ui.h>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <stack>
#include <unordered_set>
using namespace l2c;
// 自分で作った型
struct Segment {
	int16_t pos;
	int32_t rawid;
	bool operator==(const Segment& rhs) const
	{
		return pos == rhs.pos && rawid == rhs.rawid;
	}
	bool operator<(const Segment& rhs) const {
		if (pos == rhs.pos) {
			return rawid < rhs.rawid;
		}
		return pos < rhs.pos;
	}
};
class Line {
public:
	std::vector<int> points;
	int prev_edge;
	int next_edge;

};

// 自分で作った型をunorderdコンテナに入れたいときは、operator== の他に
// 型と同じ名前空間で hash_value 関数を定義

size_t hash_value(const Segment& d)
{
	// 複数の値のハッシュ値を組み合わせてハッシュ値を計算するには、
	// boost::hash_combine を使います。
	size_t h = 0;
	boost::hash_combine(h, d.pos);
	boost::hash_combine(h, d.rawid);
	return h;
}

class Chain_baselist {
public:
	int64_t groupid, trackid;
	std::pair<int32_t, int64_t> target_track;
	std::set<std::pair<int32_t, int64_t>> btset;
	std::set<std::tuple<int, int, int, int>>ltlist;
	std::vector<int32_t> usepos;
	int chain_num;
	std::vector<Linklet> make_link_list();
	void set_usepos();
};
class Chain_baselist_compress : public Chain_baselist
{
public:
	std::set<std::pair<int32_t, int64_t>> comp_btset;
	std::multimap<std::tuple<int, int, int, int>, std::vector<std::pair<int32_t, int64_t>>>comp_ltlist;
	std::set<std::tuple<int, int, int, int>>cut_ltlist;
	std::vector<Linklet> make_comp_link_list();
	std::vector<Linklet> make_comp_link_list(std::set<std::tuple<int, int, int, int>>& cut_list);
	void set_comp_usepos();

};

class output_format {
public:
	int groupid, trackid;
	int num_comfirmed_path, num_cut_path, num_select_path;
	std::vector<std::pair<int, std::tuple<int, int, int, int>>> comfirmed_path;
	std::vector<std::pair<int, std::tuple<int, int, int, int>>> cut_path;
	std::vector<std::pair<int, std::tuple<int, int, int, int>>> select_path;
};
class process_foramt {
	std::vector<output_format> out;
	std::vector<std::pair<int, std::tuple<int, int, int, int>>> cut_path;

};

bool sort_M_Base(const mfile0::M_Base& left, const mfile0::M_Base& right) {

	if (left.pos == right.pos) {
		return left.rawid < right.rawid;
	}
	return left.pos < right.pos;
}

std::vector < Chain_baselist > read_linklet_list2(std::string filename, int eventid, int trackid);
mfile0::Mfile SetMfile(std::vector < Chain_baselist >& c, std::multimap<int, mfile0::M_Base*>& base);


int main(int argc, char** argv) {
	if (argc != 6) {
		fprintf(stderr, "usage:prg file_in_link file-in-ECC eventid trackid file_out-mfile\n");

		exit(1);
	}

	std::string file_in_link_list = argv[1];
	std::string file_in_ECC = argv[2];
	int event_id = std::stoi(argv[3]);
	int track_id = std::stoi(argv[4]);
	std::string file_out_mfile = argv[5];

	std::vector<Chain_baselist> chain_list = read_linklet_list2(file_in_link_list, event_id, track_id);

	std::vector<corrmap0::Corrmap> corr;
	corrmap0::read_cormap(file_in_ECC + "\\Area0\\0\\align\\corrmap-abs.lst", corr);
	std::map<int, corrmap0::Corrmap> corr_map;
	for (auto itr = corr.begin(); itr != corr.end(); itr++) {
		corr_map.insert(std::make_pair(itr->pos[0] / 10, *itr));
	}
	//gap nominal read
	std::stringstream structure_path;
	structure_path << file_in_ECC << "\\st\\st.dat";
	chamber1::Chamber chamber;
	chamber1::read_structure(structure_path.str(), chamber);
	std::map<int, double> z_map = chamber1::base_z_convert(chamber);

	printf("z map ize=%d\n", z_map.size());

	std::multimap<int, mfile0::M_Base*>base;
	mfile0::Mfile m = SetMfile(chain_list, base);
	printf("base size =%d\n", base.size());
	vxx::BvxxReader br;
	for (auto itr = base.begin(); itr != base.end(); itr++) {
		if (corr_map.count(itr->first) == 0)continue;
		if (z_map.count(itr->first) == 0)continue;

		corrmap0::Corrmap param = corr_map.at(itr->first);
		double nominal_z = z_map.at(itr->first);


		std::stringstream file_in_base;
		file_in_base << file_in_ECC << "\\Area0\\PL"
			<< std::setw(3) << std::setfill('0') << itr->first << "\\b"
			<< std::setw(3) << std::setfill('0') << itr->first << ".sel.cor.vxx";

		std::array<int, 2> index = { itr->second->rawid,itr->second->rawid + 1 };//1234<=rawid<=5678であるようなものだけを読む。


		std::vector<vxx::base_track_t> base = br.ReadAll(file_in_base.str(), itr->first, 0, vxx::opt::index = index);

		itr->second->ax = base[0].ax * param.angle[0] + base[0].ay * param.angle[1] + param.angle[4];
		itr->second->ay = base[0].ax * param.angle[2] + base[0].ay * param.angle[3] + param.angle[5];
		itr->second->x = base[0].x * param.position[0] + base[0].y * param.position[1] + param.position[4];
		itr->second->y = base[0].x * param.position[2] + base[0].y * param.position[3] + param.position[5];
		itr->second->z = nominal_z + param.dz;
		itr->second->ph = base[0].m[0].ph + base[0].m[1].ph;
	}

	mfile0::write_mfile(file_out_mfile, m);

}


std::vector < Chain_baselist > read_linklet_list2(std::string filename, int eventid, int trackid) {
	std::vector < Chain_baselist > ret;
	std::ifstream ifs(filename);
	int gid, tid, muon_num, t_pl, t_rawid;
	int64_t link_num;
	int pl, raw;
	double weight;
	std::tuple<int, int, int, int> link;
	int count = 0;

	while (ifs >> gid >> tid >> t_pl >> t_rawid >> link_num) {
		printf("\r read group %d gid=%d link=%lld", count, gid, link_num);
		count++;
		Chain_baselist c;
		c.groupid = gid;
		c.trackid = tid;
		c.target_track.first = t_pl;
		c.target_track.second = t_rawid;
		for (int64_t i = 0; i < link_num; i++) {
			ifs >> std::get<0>(link) >> std::get<1>(link) >> std::get<2>(link) >> std::get<3>(link);
			c.btset.insert(std::make_pair(std::get<0>(link), std::get<2>(link)));
			c.btset.insert(std::make_pair(std::get<1>(link), std::get<3>(link)));
			c.ltlist.insert(link);
		}
		if (eventid != c.groupid) {
			continue;
		}
		if (trackid != c.trackid) {
			continue;
		}

		ret.push_back(c);
	}
	printf("\r read group %d gid=%d link=%lld\n", count, gid, link_num);

	return ret;

}

mfile0::Mfile SetMfile(std::vector < Chain_baselist >& c, std::multimap<int, mfile0::M_Base*>& base) {

	mfile0::Mfile ret;
	mfile0::set_header(3, 133, ret);
	int cid = 0;
	for (auto itr = c.begin(); itr != c.end(); itr++) {
		for (auto itr2 = itr->btset.begin(); itr2 != itr->btset.end(); itr2++) {
			mfile0::M_Chain c;
			mfile0::M_Base b;

			b.group_id = itr->groupid;
			b.pos = itr2->first + 1;
			b.rawid = itr2->second;
			b.flg_d[0] = 0;
			b.flg_d[1] = 0;
			b.flg_i[0] = 0;
			b.flg_i[1] = 0;
			b.flg_i[2] = 0;
			b.flg_i[3] = 0;

			c.basetracks.push_back(b);
			c.chain_id = cid;
			cid++;
			c.nseg = c.basetracks.size();
			c.pos0 = c.basetracks.begin()->pos;
			c.pos1 = c.basetracks.rbegin()->pos;
			ret.chains.push_back(c);
		}
	}

	for (auto itr = ret.chains.begin(); itr != ret.chains.end(); itr++) {
		for (auto itr2 = itr->basetracks.begin(); itr2 != itr->basetracks.end(); itr2++) {
			base.insert(std::make_pair(itr2->pos / 10, &(*itr2)));
		}
	}

	return ret;
}
