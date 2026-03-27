#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

//muon‚Ìchain‚É‘Î‚µ‚ÄAlinklet‚©‚çÄ“xchain‚ğÄ¶¬‚·‚é
//1group 1chain‚ª‘O’ñ
#define _CRT_SECURE_NO_WARNINGS
#define DEFAULT_CHAIN_UPPERLIM 1000000000000000

//#pragma comment(lib, "VxxReader.lib")
//#include "VxxReader.h"
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>
#pragma comment(lib, "lib_l2cx.lib")
#include <LibL2c-x.h>

#include <filesystem>
#include <set>
#include <netscan_data_types_ui.h>
#include <algorithm>
#include <iostream>
#include <unordered_map>
//#include <boost/unordered_set.hpp>
//#include <boost/unordered_map.hpp>
#include <stack>
#include <unordered_set>

mfile0::M_Chain make_chain(Momentum_recon::Mom_chain& mom, int groupid);

std::map<int, corrmap0::Corrmap> read_corrmap_abs(std::string file_in_ECC);

void add_base_information(std::string file_in_ECC, std::multimap<int, mfile0::M_Base*>& base_map_single);
void trans_local(std::multimap<int, mfile0::M_Base*>& base_map_single, std::map<int, std::vector<corrmap_3d::align_param2>>& corr);

int main(int argc, char** argv) {
	if (argc != 4) {
		fprintf(stderr, "prg file-in-momch.momch file-out-mfile.all event-list.txt\n");
		fprintf(stderr, "event-list.txt:groupid chainid\n");
		exit(1);
	}

	std::string file_in_momch = argv[1];
	std::string file_out_mfile = argv[2];
	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch);

	std::string event_list = argv[3];

	std::set<int> eventid;
	int tmp[2];
	std::ifstream ifs(event_list);
	if (!ifs) {
		std::cerr << "failed to open event_list." << std::endl;
		exit(1);
	}
	while (ifs >> tmp[0] >> tmp[1]) {
		eventid.insert(tmp[0]);
	}

	//std::multimap<std::pair<int, int>, mfile0::M_Base*>base_map_link;
	std::multimap<int, mfile0::M_Base*>base_map_pl;

	mfile0::Mfile m;
	mfile0::set_header(1, 133, m);
	for (auto& ev_v : momch) {
		for (auto& c : ev_v.chains) {
			auto p = eventid.find(ev_v.groupid);
			if (p != eventid.end()) {
				//printf("%10d %10d %d\n", ev_v.groupid, c.chainid, c.base.size());
				mfile0::M_Chain chain = make_chain(c, ev_v.groupid);
				m.chains.push_back(chain);
			}
		}
	}
	mfile0::write_mfile(file_out_mfile, m);

}
std::set<int>list(std::string input){
	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "Fail reading event-list." << std::endl;
		exit(1);
	}

	int id[2];
	std::set<int> l;
	while (ifs >> id[0] >> id[1]) {
		l.insert(id[0]);
	}
	return l;
}

mfile0::M_Chain make_chain(Momentum_recon::Mom_chain& mom, int groupid) {
	mfile0::M_Chain c;
	for (auto itr = mom.base.begin(); itr != mom.base.end(); itr++) {
		mfile0::M_Base b;
		b.pos = itr->pl * 10;
		b.rawid = itr->rawid;
		b.group_id = groupid;
		b.ax = itr->ax;
		b.ay = itr->ay;
		b.x = itr->x;
		b.y = itr->y;
		b.z = itr->z;
		b.ph = itr->m[0].ph + itr->m[1].ph;
		b.flg_i[0] = 0;
		b.flg_i[1] = 0;
		b.flg_i[2] = 0;
		b.flg_i[3] = 0;
		b.flg_d[0] = 0;
		b.flg_d[1] = 0;
		c.basetracks.push_back(b);
	}
	c.chain_id = mom.chainid;
	c.pos0 = c.basetracks.begin()->pos;
	c.pos1 = c.basetracks.rbegin()->pos;
	c.nseg = c.basetracks.size();
	return c;
}

