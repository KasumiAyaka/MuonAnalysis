// 2024/07/28
// kasumi
// set most upsetream pl then divide muons whose most upstream muon is more than set pl.

// 2024/05/09
// kasumi
// based on "Check_upstream_base.cpp

#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

struct Lst {
	int gid, cid;
};
bool operator<(const Lst& lhs, const Lst& rhs) {
	return std::tie(lhs.gid, lhs.cid) < std::tie(rhs.gid, rhs.cid);
}


void add_momch(std::vector<Momentum_recon::Event_information>& momch, std::vector<Momentum_recon::Event_information>& momch_add, int add_event, int add_chain);
void read_list(std::string in, std::set<int>& list);
void erase_momch(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch, std::set<int>& set);
void Divide_momch(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch1, std::vector<Momentum_recon::Event_information>& momch2,int thr_penet);

int main(int argc, char** argv) {
	if (argc != 5) {
		fprintf(stderr, "usage:input.momch output.momch divided.momch penetrate-pl-thr(int)\n#pl>thr-->penetrate");

		exit(1);
	}
	std::string file_in_momch_ori = argv[1];//input momch
	std::string file_out_momch = argv[2];//erase list
	std::string file_penet_momch = argv[3];//output
	int thr = std::stoi(argv[4]);

	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch_ori);
	std::vector<Momentum_recon::Event_information> new_momch;
	std::vector<Momentum_recon::Event_information> penet_momch;
	std::cout << "ok1?" << std::endl;

	Divide_momch(momch, new_momch,penet_momch,thr);
	std::cout << "ok2?" << std::endl;

	Momentum_recon::Write_Event_information_extension(file_out_momch, new_momch);
	Momentum_recon::Write_Event_information_extension(file_penet_momch, penet_momch);
	std::cout << "fin." << std::endl;

}
void add_momch(std::vector<Momentum_recon::Event_information>& momch, std::vector<Momentum_recon::Event_information>& momch_add, int add_event, int add_chain) {
	Momentum_recon::Mom_chain chain;
	for (auto& ev : momch_add) {
		if (ev.groupid != add_event)continue;
		for (auto& c : ev.chains) {
			if (c.chainid != add_chain)continue;
			chain = c;
		}
	}

	for (auto& ev : momch) {
		if (ev.groupid != add_event)continue;
		ev.chains.push_back(chain);
	}

}
void read_list(std::string in, std::set<int>& list) {
	std::ifstream ifs(in);
	if (!ifs) {
		std::cerr << "reading file error" << std::endl;
	}

	int i;
	while (ifs >> i) {
		list.insert(i);
	}
	std::cout << "size:" << list.size() << std::endl;
}
void erase_momch(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch, std::set<int>& set) {
	Momentum_recon::Mom_chain chain;

	for (auto& ev : momch0) {
		for (auto itr = set.begin(); itr != set.end(); itr++) {
			if (*itr == ev.groupid) {
				std::cout << "erase " << ev.groupid << std::endl;
			}
			else {
				momch.push_back(ev);

			}
		}
		//auto itr = set.find(ev.groupid);
		//if (itr != set.end()) {
		//	momch.push_back(ev);
		//}
		//else {
		//	std::cout << ev.groupid << std::endl;
		//}
	}

}
void Divide_momch(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch1, std::vector<Momentum_recon::Event_information>& momch2, int thr_penet) {
	Momentum_recon::Mom_chain chain;
	int flg_penet = 0;
	int pl = 0;
	for (auto& ev : momch0) {
		for (auto& c : ev.chains) {

			for (auto itr = c.base.begin(); itr != c.base.end(); itr++) {
				if (itr->pl > thr_penet) {
					flg_penet++;
					pl = itr->pl;
				}
			}

			if (flg_penet == 0) {//original-penetrate
				momch1.push_back(ev);

			}
			else {// penetrate muon
				momch2.push_back(ev);
				std::cout << "Divide penetrate muon : groupid  " << ev.groupid << " PL" << pl << std::endl;
				flg_penet = 0;
			}
		}
	}


}