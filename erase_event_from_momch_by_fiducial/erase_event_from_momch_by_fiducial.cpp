//“r’†

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
void edgecut_momch(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch_in, std::vector<Momentum_recon::Event_information>& momch_out, double d);

int main(int argc, char** argv) {
	if (argc != 5) {
		fprintf(stderr, "usage:input.momch  output_inner_fiducial.momch output_outer_fiducial.momch  edgecut[cm]\n");
		exit(1);
	}
	std::string file_in_momch_ori = argv[1];//input momch
	std::string file_out_momch_inner = argv[2];//output
	std::string file_out_momch_outer = argv[3];//output
	double d = std::stod(argv[4]);

	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch_ori);
	std::vector<Momentum_recon::Event_information> new_momch_inner;
	std::vector<Momentum_recon::Event_information> new_momch_outer;


	edgecut_momch(momch, new_momch_inner, new_momch_outer,d);

	Momentum_recon::Write_Event_information_extension(file_out_momch_inner, new_momch_inner);
	Momentum_recon::Write_Event_information_extension(file_out_momch_outer, new_momch_outer);
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
		auto itr = set.find(ev.groupid);
		if (itr != set.end()) {
			momch.push_back(ev);
		}
		else {
			std::cout << "erase : " << ev.groupid << std::endl;
		}
	}

}
void edgecut_momch(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch_in, std::vector<Momentum_recon::Event_information>& momch_out, double d) {

	std::set<int>gid_list;
	double x, y;
	int pl,flg;
	for (auto& ev : momch0) {
		flg = 0;
		for (auto& c : ev.chains) {
			if (c.chainid == 0) {
				pl=c.base.rbegin()->pl;
				x=c.base.rbegin()->x;
				y=c.base.rbegin()->y;
				//std::cout << "check begin :" << c.base.begin()->pl << " " << c.base.begin()->x << ", " << c.base.begin()->y << ", " << std::endl;//pl003
				//std::cout << "check rbegin:" << pl << " " << x << ", " << y << ", " << std::endl;//upl
				if (x >= d*10000 && x <= 250000 - d * 10000 && y >= d * 10000 && y <= 25 * 10000 - d * 10000) {
					//std::cout << "check rbegin:" << pl << " " << x << ", " << y << ", " << std::endl;//upl

					// fkg in
					flg++;
				}
			}
		}
		if (flg > 0) {
			gid_list.insert(ev.groupid);
		}
	}

	for (auto& ev : momch0) {
		auto itr = gid_list.find(ev.groupid);//erase
		if (itr != gid_list.end()) {// found erase event
			momch_in.push_back(ev);
		}
		else {//remain
			momch_out.push_back(ev);
		}
	}

	std::cout << " * Central Area = " << momch_in.size() << std::endl;
	std::cout << " * Edge Area    = " << momch_out.size() << std::endl;

}