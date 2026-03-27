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
void erase_momch_by_unixtime(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch, std::set<int>& set);

int main(int argc, char** argv) {
	if (argc < 4||argc >5) {
		fprintf(stderr, "usage:input.momch elase-list.txt output.momch\neraselist:eventid\n");
		fprintf(stderr, "usage:input.momch elase-list.txt output.momch 0 \neraselist:unixtime\n");
		exit(1);
	}
	std::string file_in_momch_ori = argv[1];//input momch
	std::string file_in_momch_add = argv[2];//erase list
	//int add_event = std::stoi(argv[3]);
	//int add_chain = std::stoi(argv[4]);
	std::string file_out_momch = argv[3];//output

	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch_ori);
	std::vector<Momentum_recon::Event_information> new_momch;
	//std::vector<Momentum_recon::Event_information> momch_add = Momentum_recon::Read_Event_information_extension(file_in_momch_add);
	//add_momch(momch, momch_add, add_event, add_chain);
	std::set<int>list;
	read_list(file_in_momch_add, list);
	if (argc == 4) {
		erase_momch(momch, new_momch, list);
	}
	else if (argc == 5) {
		erase_momch_by_unixtime(momch, new_momch, list);
	}

	Momentum_recon::Write_Event_information_extension(file_out_momch, new_momch);
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
		std::cerr << "reading file error\n" << in << std::endl;
	}

	int i;
	while (ifs >> i) {
		list.insert(i);
	}
	std::cout << " * The # of erase events:" << list.size() << std::endl;
}
void erase_momch(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch, std::set<int>& set) {
	Momentum_recon::Mom_chain chain;

	std::cout << " * Before erase : " << momch0.size() << std::endl;
	for (auto& ev : momch0) {
		//for (auto itr = set.begin(); itr != set.end(); itr++) {
		//	if (*itr == ev.groupid) {
		//		std::cout << "erase "<<ev.groupid << std::endl;
		//	}
		//	else {
		//		momch.push_back(ev);
		//	}
		//}
		auto itr = set.find(ev.groupid);//erase
		if (itr != set.end()) {// found erase event
			std::cout << "    erase : " << ev.groupid << std::endl;
		}
		else {//remain
			momch.push_back(ev);
		}
	}
	std::cout << " * After  erase : " << momch.size() << std::endl;

}

void erase_momch_by_unixtime(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch, std::set<int>& set) {
	Momentum_recon::Mom_chain chain;

	int cnt = 0;
	for (auto& ev : momch0) {
		//for (auto itr = set.begin(); itr != set.end(); itr++) {
		//	if (*itr == ev.groupid) {
		//		std::cout << "erase "<<ev.groupid << std::endl;
		//	}
		//	else {
		//		momch.push_back(ev);
		//	}
		//}
		auto itr = set.find(ev.unix_time);//erase
		if (itr != set.end()) {// found erase event
			std::cout << "  erase : " << ev.groupid << std::endl;
			cnt++;
		}
		else {//remain
			momch.push_back(ev);
		}
	}
	std::cout << " * Before erase : " << momch0.size() << std::endl;
	std::cout << " * After  erase : " << momch.size() << std::endl;
	std::cout << " * Erased event num = " << cnt << std::endl;
}