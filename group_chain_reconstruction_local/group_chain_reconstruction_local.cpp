//muonのchainに対して、linkletから再度chainを再生成する
//1group 1chainが前提
#define _CRT_SECURE_NO_WARNINGS
#define DEFAULT_CHAIN_UPPERLIM 1000000000000000

#pragma comment(lib, "VxxReader.lib")
#include "VxxReader.h"
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>
#pragma comment(lib, "Libl2c-x.lib")
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

class Track_file {
public:
	int eventid, trackid, pl, rawid;
};
class Group_file :public Track_file {
public:
	int link_num;
	std::vector<std::tuple<int, int, int, int>> linklet;
};
using namespace l2c;

l2c::Cdat l2c_x(std::vector<Linklet>& ltlist, std::vector<int32_t>& usepos, int64_t upperlim, bool output);
std::vector<mfile0::M_Chain> make_chain(Group_file& g);
mfile0::M_Chain make_chain_signle_base(Group_file& g);

std::map<int, corrmap0::Corrmap> read_corrmap_abs(std::string file_in_ECC);
std::vector <  Group_file > read_group_file(std::string filename);
mfile0::Mfile select_target_mfile(mfile0::Mfile& m);

void add_base_information(std::string file_in_ECC, std::multimap<int, mfile0::M_Base*>& base_map_single);
void trans_local(std::multimap<int, mfile0::M_Base*>& base_map_single, std::map<int, std::vector<corrmap_3d::align_param2>>& corr);

int main(int argc, char** argv) {
	if (argc != 5) {
		fprintf(stderr, "prg file-in-group file-in-ECC file-out-tragetchain file-out-all-chain\n");
		exit(1);
	}

	std::string file_in_group = argv[1];
	std::string file_in_ECC = argv[2];
	std::string file_out_mfile = argv[3];
	std::string file_out_mfile_all = argv[4];
	//ここ

	//	std::map<int, corrmap0::Corrmap> corrmap = read_corrmap_abs(file_in_ECC);
	//std::string file_in_corrmap = argv[3];

	std::string file_in_corrmap = file_in_ECC + "\\Area0\\0\\align\\fine\\local\\corrmap-local-abs.lst";
	std::map<int, std::vector<corrmap_3d::align_param>>corrmap = corrmap_3d::read_ali_param_abs(file_in_corrmap, 1);
	std::map<int, std::vector<corrmap_3d::align_param2>>corrmap_dd = corrmap_3d::DelaunayDivide_map(corrmap);


	std::vector <  Group_file >group = read_group_file(file_in_group);

	//std::multimap<std::pair<int, int>, mfile0::M_Base*>base_map_link;
	std::multimap<int, mfile0::M_Base*>base_map_pl;

	mfile0::Mfile m;
	mfile0::set_header(1, 133, m);
	for (int i = 0; i < group.size(); i++) {
		printf("group=%d %d/%d\r", group[i].eventid, i, group.size());
		if (group[i].linklet.size() == 0) {
			//single basetrackの場合
			mfile0::M_Chain chain = make_chain_signle_base(group[i]);
			m.chains.push_back(chain);
			mfile0::M_Chain* cp = &m.chains[m.chains.size() - 1];
			for (int k = 0; k < cp->basetracks.size(); k++) {
				base_map_pl.insert(std::make_pair(cp->basetracks[k].pos / 10, &cp->basetracks[k]));
			}

		}
		if (group[i].linklet.size() > 0) {
			std::vector<mfile0::M_Chain> chain = make_chain(group[i]);
			for (int j = 0; j < chain.size(); j++) {
				m.chains.push_back(chain[j]);
				mfile0::M_Chain* cp = &m.chains[m.chains.size() - 1];
				for (int k = 0; k < cp->basetracks.size(); k++) {
					base_map_pl.insert(std::make_pair(cp->basetracks[k].pos / 10, &cp->basetracks[k]));
				}
			}

		}
	}
	printf("\n");
	printf("chain size=%d\n", m.chains.size());
	add_base_information(file_in_ECC, base_map_pl);

	//ここ
	//trans_global(base_map_pl, corrmap, z_map);
	trans_local(base_map_pl, corrmap_dd);


	mfile0::Mfile m_sel = select_target_mfile(m);

	mfile0::write_mfile(file_out_mfile, m_sel);
	mfile0::write_mfile(file_out_mfile_all, m);
}

std::map<int, corrmap0::Corrmap> read_corrmap_abs(std::string file_in_ECC) {
	std::stringstream file_in_corr;
	file_in_corr << file_in_ECC << "\\Area0\\0\\align\\corrmap-abs.lst";

	std::vector<corrmap0::Corrmap> corr;
	corrmap0::read_cormap(file_in_corr.str(), corr);

	std::map<int, corrmap0::Corrmap> ret;
	for (auto itr = corr.begin(); itr != corr.end(); itr++) {
		ret.insert(std::make_pair(itr->pos[0] / 10, *itr));
	}
	return ret;
}
std::vector <  Group_file > read_group_file(std::string filename) {
	std::vector <  Group_file > ret;
	std::ifstream ifs(filename);
	int gid, tid, muon_num, link_num, t_pl, t_rawid;
	int pl, raw;
	double weight;
	std::tuple<int, int, int, int> link;
	int count = 0;

	while (ifs >> gid >> tid >> t_pl >> t_rawid >> link_num) {
		printf("\r read group %d", count);
		count++;
		Group_file g;
		g.eventid = gid;
		g.trackid = tid;
		g.pl = t_pl;
		g.rawid = t_rawid;
		for (int i = 0; i < link_num; i++) {
			ifs >> std::get<0>(link) >> std::get<1>(link) >> std::get<2>(link) >> std::get<3>(link);
			g.linklet.push_back(link);
		}
		g.link_num = g.linklet.size();
		ret.push_back(g);
	}
	printf("\r read group %d\n", count);

	return ret;

}

l2c::Cdat l2c_x(std::vector<Linklet>& ltlist, std::vector<int32_t>& usepos, int64_t upperlim, bool output) {
	try
	{
		size_t possize = usepos.size();

		//l2c関数本体を呼び出す。
		//出力されるCdatは基本的にl2c-xと同等のものになっているはずだが、若干の違いがある。
		//1. upperlimを超過したGroupは、Chainは持たず、全BaseTrackの情報だけを持った状態で出力される。
		//2. upperlim超過のGroupも含め、属す全BaseTrackの情報を持っている。
		//またplate枚数は最大64枚。変更するにはLibL2c-xのリビルドを要する。
		//デフォルトでoutput_isolated_linkletはfalse、upperlimは無効になっている。useposは与えないとエラーになる。
		if (output) {
			l2c::Cdat cdat = l2c::MakeCdat(ltlist, opt::usepos = usepos, opt::upperlim = upperlim, opt::output_isolated_linklet = true);
			return cdat;
		}
		l2c::Cdat cdat = l2c::MakeCdat_no_out(ltlist, opt::usepos = usepos, opt::upperlim = upperlim, opt::output_isolated_linklet = true);
		return cdat;

	}
	catch (const std::exception& e)
	{
		fprintf(stderr, "%s\n", e.what());
	}

}

std::vector<mfile0::M_Chain> make_chain(Group_file& g) {

	std::vector<Linklet> ltlist;
	std::set<int32_t> pos_list;


	ltlist.reserve(g.linklet.size());

	for (auto itr = g.linklet.begin(); itr != g.linklet.end(); itr++) {
		ltlist.emplace_back(std::get<0>(*itr), std::get<1>(*itr), std::get<2>(*itr), std::get<3>(*itr));
		pos_list.insert(std::get<0>(*itr));
		pos_list.insert(std::get<1>(*itr));
	}
	std::vector<int32_t> usepos(pos_list.begin(), pos_list.end());

	l2c::Cdat cdat = l2c_x(ltlist, usepos, DEFAULT_CHAIN_UPPERLIM, false);
	std::vector<mfile0::M_Chain> ret;
	size_t grsize = cdat.GetNumOfGroups();
	size_t possize = usepos.size();
	for (size_t grid = 0; grid < grsize; ++grid)
	{
		const Group& gr = cdat.GetGroup(grid);
		size_t chsize = gr.GetNumOfChains();
		if (grid != gr.GetID()) throw std::exception("grid != gr.GetID().");//gridはgr.GetID()の戻り値と基本的に等しいはずである。
		if (gr.IsOverUpperLim())
		{
			fprintf(stdout, "over upperlim.nchain:%-9lld\n", chsize);
			fprintf(stdout, "event:%d track%d\n", g.eventid, g.trackid);
			exit(1);
			//continue;
		}
		for (size_t ich = 0; ich < chsize; ++ich)
		{
			mfile0::M_Chain c;
			c.chain_id = ich;
			Chain ch = gr.GetChain(ich);//仕様上、Chainオブジェクトは一時変数として戻ってくる。
			int64_t chid = ch.GetID();
			int32_t nseg = ch.GetNSeg();
			int32_t spl = ch.GetStartPL();
			int32_t epl = ch.GetEndPL();

			for (size_t pl = 0; pl < possize; ++pl)
			{
				//ここで用いるplも上と同様に、useposの中で何番目のposかを意味する。
				//pl == 5だとすると、usepos中で（0から数えて）5番目、つまりpos==400のBaseTrack情報が返ってくる。
				BaseTrackID bt = ch.GetBaseTrack(pl);
				if (bt.IsEmpty()) continue;//IsEmptyがtrueのときは、このBaseTrackは歯抜け。
				int32_t btpl = bt.GetPL();
				int64_t btid = bt.GetRawID();
				//btplとplは厳密に一致しなければおかしい。
				if (btpl != pl) throw std::exception("BaseTrackID::PL != pl.");

				mfile0::M_Base b;
				b.pos = usepos[bt.GetPL()] + 1;
				b.rawid = btid;
				b.group_id = g.eventid;
				b.flg_i[0] = 0;
				b.flg_i[1] = g.trackid;
				b.flg_i[2] = 0;
				b.flg_i[3] = 0;
				b.flg_d[0] = 0;
				b.flg_d[1] = 0;
				if (g.pl == b.pos / 10 && g.rawid == b.rawid) {
					b.flg_i[0] = 1;
				}
				c.basetracks.push_back(b);
			}
			c.pos0 = c.basetracks.begin()->pos;
			c.pos1 = c.basetracks.rbegin()->pos;
			c.nseg = c.basetracks.size();
			ret.push_back(c);
		}
	}
	return ret;

}
mfile0::M_Chain make_chain_signle_base(Group_file& g) {
	mfile0::M_Chain c;
	mfile0::M_Base b;
	b.pos = g.pl * 10;
	b.rawid = g.rawid;
	b.group_id = g.eventid;
	b.flg_i[0] = 1;
	b.flg_i[1] = g.trackid;
	b.flg_i[2] = 0;
	b.flg_i[3] = 0;
	b.flg_d[0] = 0;
	b.flg_d[1] = 0;
	c.basetracks.push_back(b);
	c.chain_id = 0;
	c.pos0 = c.basetracks.begin()->pos;
	c.pos1 = c.basetracks.rbegin()->pos;
	c.nseg = c.basetracks.size();
	return c;
}

mfile0::Mfile select_target_mfile(mfile0::Mfile& m) {
	mfile0::Mfile ret;
	ret.header = m.header;
	for (auto& c : m.chains) {
		bool flg = false;
		for (auto& b : c.basetracks) {
			if (b.flg_i[0] == 1)flg = true;
		}
		if (flg) {
			ret.chains.push_back(c);
		}
	}
	return ret;
}


void add_base_information(std::string file_in_ECC, std::multimap<int, mfile0::M_Base*>& base_map_single) {
	int count = 0;
	int64_t pickup_base_num = 0;
	for (auto itr = base_map_single.begin(); itr != base_map_single.end(); itr++) {
		count = base_map_single.count(itr->first);
		int pl = itr->first;
		if (pl < 3 || pl>133) continue;//why PL000??
		printf("PL%03d basetrack read ", pl);

		std::multimap<int, mfile0::M_Base*>base_map;
		int raw_min = INT32_MAX, raw_max = 0;
		auto range = base_map_single.equal_range(itr->first);
		for (auto res = range.first; res != range.second; res++) {
			base_map.insert(std::make_pair(res->second->rawid, res->second));
			raw_min = std::min(raw_min, int(res->second->rawid));
			raw_max = std::max(raw_max, int(res->second->rawid));
		}

		std::stringstream file_in_base;
		file_in_base << file_in_ECC << "\\Area0\\PL"
			<< std::setw(3) << std::setfill('0') << pl << "\\b"
			<< std::setw(3) << std::setfill('0') << pl << ".sel.cor.vxx";
		std::vector<vxx::base_track_t >base;
		vxx::BvxxReader br;
		std::array<int, 2> index = { raw_min,raw_max + 1 };//1234<=rawid<=5678であるようなものだけを読む。
		printf("rawid %10d - %10d\n", raw_min, raw_max);
		base = br.ReadAll(file_in_base.str(), pl, 0, vxx::opt::index = index);

		for (auto itr = base.begin(); itr != base.end(); itr++) {
			auto res = base_map.find(itr->rawid);
			if (res == base_map.end())continue;
			auto range = base_map.equal_range(itr->rawid);
			for (auto res = range.first; res != range.second; res++) {
				res->second->ax = itr->ax;
				res->second->ay = itr->ay;
				res->second->x = itr->x;
				res->second->y = itr->y;
				res->second->z = 0;
				res->second->ph = itr->m[0].ph + itr->m[1].ph;
				pickup_base_num++;
			}
		}

		itr = std::next(itr, count - 1);
	}
}

void trans_local(std::multimap<int, mfile0::M_Base*>& base_map_single, std::map<int, std::vector<corrmap_3d::align_param2>>& corr) {
	int count = 0;
	for (auto itr = base_map_single.begin(); itr != base_map_single.end(); itr++) {
		count = base_map_single.count(itr->first);
		int pl = itr->first;
		printf("PL%03d basetrack tans\n", pl);
		if (pl < 3 || pl>133) continue;//why PL000??
		if (corr.count(pl) == 0) {
			fprintf(stderr, "PL%03d corrmap not found\n", pl);
		}
		std::vector<corrmap_3d::align_param2> param = corr.at(pl);

		std::vector< mfile0::M_Base*> base_trans;
		auto range = base_map_single.equal_range(itr->first);
		for (auto res = range.first; res != range.second; res++) {
			base_trans.push_back(res->second);
		}
		std::vector <std::pair<mfile0::M_Base*, corrmap_3d::align_param2*>> base_trans_map = corrmap_3d::track_affineparam_correspondence(base_trans, param);
		trans_base_all(base_trans_map);


		itr = std::next(itr, count - 1);
	}

}