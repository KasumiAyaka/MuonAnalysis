// 2025/03/08
// kasumi
// edit "Event_reject_cc_1seg.cpp"

#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

class ChanceCoincidence {
public:
	int eventid, chainid, pl;
};
std::vector<ChanceCoincidence> read_chancecoincidence(std::string filename);
void reject_chancecoincidence(std::vector<mfile0::M_Chain>& chain, std::vector<ChanceCoincidence>& cc_list);
void reject_1seg_chain(std::vector<mfile0::M_Chain>& chain);

int main(int argc, char** argv) {
	if (argc != 4) {
		fprintf(stderr, "usage: in.all list.txt out.all\n");
		fprintf(stderr, "list.txt : eventid chainid pl\n");
		fprintf(stderr, "caution!!\n Erase basetracks eith a PL greater than the specified PL.(=erase if basetrack-pl > thr_pl)\n");
		exit(1);
	}

	std::string file_in_mfile = argv[1];
	std::string file_in_cclist = argv[2];
	std::string file_out_mfile = argv[3];

	std::vector<ChanceCoincidence> cc_list = read_chancecoincidence(file_in_cclist);
	mfile0::Mfile m;
	mfile1::read_mfile_extension(file_in_mfile, m);

	reject_chancecoincidence(m.chains, cc_list);
	//reject_1seg_chain(m.chains);

	mfile1::write_mfile_extension(file_out_mfile, m);
}
std::vector<ChanceCoincidence> read_chancecoincidence(std::string filename) {
	std::ifstream ifs(filename);
	std::vector<ChanceCoincidence> ret;
	ChanceCoincidence cc;
	while (ifs >> cc.eventid >> cc.chainid >> cc.pl) {
		ret.push_back(cc);
	}
	printf("reject track %d\n", ret.size());
	return ret;

}

void reject_chancecoincidence(std::vector<mfile0::M_Chain>& chain, std::vector<ChanceCoincidence>& cc_list) {
	int count = 0;
	for (auto& c : chain) {
		for (auto itr = c.basetracks.begin(); itr != c.basetracks.end(); ) {
			bool flg = false;

			for (auto& l : cc_list) {
				if (l.eventid != itr->group_id)continue;
				if (l.chainid != c.chain_id)continue;
				if (l.pl >= int(itr->pos / 10)) {
					//std::cout << itr->pos / 10 << " " << l.pl << std::endl;
					continue;
				}
				flg = true;
				count++;
			}
			if (flg) {
				itr = c.basetracks.erase(itr);
			}
			else {
				itr++;
			}
		}
		c.nseg = c.basetracks.size();
		c.pos0 = c.basetracks.begin()->pos;
		c.pos1 = c.basetracks.rbegin()->pos;
	}
	printf("remove basetrack %d\n", count);

}

void reject_1seg_chain(std::vector<mfile0::M_Chain>& chain) {
	for (auto itr = chain.begin(); itr != chain.end();) {
		if (itr->nseg == 1) {
			itr = chain.erase(itr);
		}
		else {
			itr++;
		}
	}
}
