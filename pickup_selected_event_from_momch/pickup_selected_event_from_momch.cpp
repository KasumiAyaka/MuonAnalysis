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
void pickup_momch(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch, std::set<int>& set);
void pickup_momch_by_unixtime(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch, std::set<int>& set);

int main(int argc, char** argv) {
	if (argc < 4||argc>5) {
		fprintf(stderr, "usage:input.momch pickup_event_list.txt output.momch \n");
		fprintf(stderr, "usage:input.momch pickup_event_list.txt output.momch 0\n");
		exit(1);
	}
	std::string file_in_momch_ori = argv[1];//input momch
	std::string file_in_momch_add = argv[2];//pickup event list
	std::string file_out_momch = argv[3];//output

	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch_ori);
	std::vector<Momentum_recon::Event_information> new_momch;

	std::set<int>list;
	read_list(file_in_momch_add, list);
	
	if (argc == 4) {
		pickup_momch(momch, new_momch, list);
	}
	else if (argc == 5) {
		pickup_momch_by_unixtime(momch, new_momch, list);
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
		std::cerr << "reading file error" << std::endl;
	}

	int i;
	while (ifs >> i) {
		list.insert(i);
	}
	std::cout << "list size:" << list.size() << std::endl;
}
void pickup_momch(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch, std::set<int>& set) {
	Momentum_recon::Mom_chain chain;

	int cnt = 0;
	for (auto& ev : momch0) {
		auto itr = set.find(ev.groupid);
		if (itr != set.end()) {
			momch.push_back(ev);
			cnt++;
			std::cout << "found : " << ev.groupid << std::endl;
		}
		else {
			//std::cout << "erase : " << ev.groupid << std::endl;
		}
	}
	std::cout << "pick up size : " << cnt << std::endl;

}
void pickup_momch_by_unixtime(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch, std::set<int>& set) {
	Momentum_recon::Mom_chain chain;

	int cnt = 0;
	for (auto& ev : momch0) {
		auto itr = set.find(ev.unix_time);
		if (itr != set.end()) {
			momch.push_back(ev);
			cnt++;
			std::cout << "found : " << ev.groupid << std::endl;
		}
		else {
			//std::cout << "erase : " << ev.groupid << std::endl;
		}
	}
	std::cout << "pick up size : " << cnt << std::endl;

}

