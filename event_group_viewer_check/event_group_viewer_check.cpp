#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>


std::map<int, std::vector< mfile0::M_Chain>> group_divide(std::vector<mfile0::M_Chain>& chains);

std::multimap<std::pair<int, int>, int> base_to_flg1_map(std::vector< mfile0::M_Chain>& chains);
std::multimap< int, mfile0::M_Chain> flg1_chain_map(std::vector< mfile0::M_Chain>& chains);

std::pair<int, int> attach_base_pickup(mfile0::M_Chain& c, int pl);

int main(int argc, char** argv) {
	if (argc != 4) {
		fprintf(stderr, "usage:file-in-event-mfile file-in-gruop-mfile file-out-mfile\n");
		exit(1);
	}

	std::string file_in_event_mfile = argv[1];
	std::string file_in_group_mfile = argv[2];
	std::string file_out_mfile = argv[3];

	mfile0::Mfile ev, group;
	mfile1::read_mfile_extension(file_in_event_mfile, ev);
	mfile1::read_mfile_extension(file_in_group_mfile, group);
	mfile0::Mfile out;
	out.header = ev.header;

	std::map<int, std::vector< mfile0::M_Chain>> ev_group = group_divide(ev.chains);
	std::map<int, std::vector< mfile0::M_Chain>> group_group = group_divide(group.chains);

	for (auto itr = ev_group.begin(); itr != ev_group.end(); itr++) {
		//if (itr->first == 6302 || itr->first == 7155)continue;//for ECC1
		//if (itr->first == 4429)continue;//for ECC6
		if (group_group.count(itr->first) == 0) {
			fprintf(stderr, "event %d not found\n", itr->first);
			continue;
			exit(1);
		}
		std::vector< mfile0::M_Chain> g_g = group_group.at(itr->first);
		std::multimap<std::pair<int, int>, int> b2f = base_to_flg1_map(g_g);
		std::multimap< int, mfile0::M_Chain>f2c = flg1_chain_map(g_g);

		for (auto& partner : itr->second) {
			if (partner.chain_id == 0)continue;
			std::vector< mfile0::M_Chain> chains;
			mfile0::M_Chain muon = *itr->second.begin();
			chains.push_back(muon);
			chains.push_back(partner);
			int vertex_pl = muon.basetracks.rbegin()->pos / 10;
			//attach base‚Ì’Šo
			std::pair<int, int> b_id = attach_base_pickup(partner, vertex_pl);
			if (b2f.count(b_id) == 0) {
				fprintf(stderr, "event %5d chain %5d not found\n", itr->first, partner.chain_id);
				fprintf(stderr, "PL %5d rawid=%d\n", b_id.first, b_id.second);
				//exit(1);//2026/01/05 for ECC6
			}
			//attach base‚ğŠÜ‚Şgroup‚Ì’Šo
			std::set<int> flg1_set;
			auto range = b2f.equal_range(b_id);
			for (auto res = range.first; res != range.second; res++) {
				flg1_set.insert(res->second);
			}
			for (auto& f : flg1_set) {
				if (f2c.count(f) == 0) {
					fprintf(stderr, "event %5d chain %5d flg_i1 %d not found\n", itr->first, partner.chain_id, f);
					exit(1);
				}
				auto range2 = f2c.equal_range(f);
				for (auto res = range2.first; res != range2.second; res++) {
					chains.push_back(res->second);
				}
			}


			//groupid chainid‚ÌU‚è’¼‚µ
			int new_gid = itr->first * 100000 + partner.chain_id * 10;
			for (int i = 0; i < chains.size(); i++) {
				chains[i].chain_id = i;
				for (auto& b : chains[i].basetracks) {
					b.group_id = new_gid;
				}
			}
			//o—Ímfile‚É‘‚«‚İ
			for (int i = 0; i < 2; i++) {
				out.chains.push_back(chains[i]);
			}
			//groupid chainid‚ÌU‚è’¼‚µ
			new_gid = itr->first * 100000 + partner.chain_id * 10 + 1;
			for (int i = 0; i < chains.size(); i++) {
				chains[i].chain_id = i;
				for (auto& b : chains[i].basetracks) {
					b.group_id = new_gid;
				}
			}
			//o—Ímfile‚É‘‚«‚İ
			for (int i = 0; i < chains.size(); i++) {
				out.chains.push_back(chains[i]);
			}
		}
	}

	mfile1::write_mfile_extension(file_out_mfile, out);




}
std::map<int, std::vector< mfile0::M_Chain>> group_divide(std::vector<mfile0::M_Chain>& chains) {

	std::multimap<int, mfile0::M_Chain> group;
	for (auto& c : chains) {
		group.insert(std::make_pair(c.basetracks[0].group_id, c));
	}

	std::map<int, std::vector< mfile0::M_Chain>>  ret;

	for (auto itr = group.begin(); itr != group.end(); itr++) {
		std::vector< mfile0::M_Chain> c;
		auto range = group.equal_range(itr->first);
		for (auto res = range.first; res != range.second; res++) {
			c.push_back(res->second);
		}
		ret.insert(std::make_pair(itr->first, c));
		itr = std::next(itr, group.count(itr->first) - 1);
	}

	return ret;

}

std::multimap<std::pair<int, int>, int> base_to_flg1_map(std::vector< mfile0::M_Chain>& chains) {
	std::multimap<std::pair<int, int>, int> ret;
	for (auto& c : chains) {
		for (auto& b : c.basetracks) {
			ret.insert(std::make_pair(std::make_pair(b.pos / 10, b.rawid), b.flg_i[1]));
		}
	}
	return ret;
}
std::multimap< int, mfile0::M_Chain> flg1_chain_map(std::vector< mfile0::M_Chain>& chains) {
	std::multimap< int, mfile0::M_Chain> ret;
	for (auto& c : chains) {
		ret.insert(std::make_pair(c.basetracks[0].flg_i[1], c));
	}
	return ret;
}

std::pair<int, int> attach_base_pickup(mfile0::M_Chain& c, int pl) {
	std::pair<int, int> ret;
	ret.first = -1;
	ret.second = -1;
	for (auto itr = c.basetracks.begin(); itr != c.basetracks.end(); itr++) {
		if (pl == itr->pos / 10 || pl + 1 == itr->pos / 10) {
			ret.first = itr->pos / 10;
			ret.second = itr->rawid;
			return ret;
		}
	}
	fprintf(stderr, "attach base not found\n");
	return ret;
}