#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

struct Lst {
	int gid, cid;
};
bool operator<(const Lst& lhs, const Lst& rhs) {
	return std::tie(lhs.gid, lhs.cid) < std::tie(rhs.gid, rhs.cid);
}


void SelectEvents(std::vector<Momentum_recon::Event_information>& in_momch,	std::vector<Momentum_recon::Event_information>& out_momch, int trknum, int pnum, std::string log);

int main(int argc, char** argv) {
	if (argc != 6) {
		fprintf(stderr, "usage:input.momch out.momch #trk #proton select_event_list.txt\nselect condition: >#trk && >#proton");
		
		exit(1);
	}
	std::string file_in_momch = argv[1];//input momch
	std::string file_out_momch = argv[2];//output
	int trknum = std::stoi(argv[3]);
	int pnum = std::stoi(argv[4]);
	std::string file_out_txt = argv[5];//output

	std::vector<Momentum_recon::Event_information> in_momch = Momentum_recon::Read_Event_information_extension(file_in_momch);
	std::vector<Momentum_recon::Event_information> out_momch;

	SelectEvents(in_momch, out_momch, trknum, pnum,file_out_txt);
	Momentum_recon::Write_Event_information_extension(file_out_momch, out_momch);
}

void SelectEvents(std::vector<Momentum_recon::Event_information>& in_momch, 
	std::vector<Momentum_recon::Event_information>& out_momch, int trknum,int pnum, std::string log) {

	std::ofstream ofs(log);
	int cnt = 0;
	int flg;
	for (auto& ev : in_momch) {
		flg = 0;
		for (auto& c : ev.chains) {
			if (c.particle_flg == 13) {
				//muon
			}
			else if (c.particle_flg == 2212) {
				//proton 
				flg++;
			}
			else if (c.particle_flg == 211) {
				//
			}
			else {
				std::cout << "no pid?" << std::endl;
			}
		}

		if (ev.chains.size() > trknum && flg > pnum) {
			// #trk > 3, p > 1
			out_momch.push_back(ev);
			cnt++;
			ofs << ev.groupid << std::endl;
		}
	}

	std::cout << " * Selected events: " << cnt << std::endl;

}