// 2024/06/14
//kasumi
//I:\NINJA\E71a\work\kasumi\ECC\MuonAnalysis\x64\Release\erase_cahin_from_momch.exe  event_water_fin3.momch ..\temp_vtx_info\elaselist.txt event_water_fin3_ECC6.momch
//2024/07/12  list Çgid cidÇ…ïœçX
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

struct Lst {
	int gid, cid, pl, rawid;
};
bool operator<(const Lst& lhs, const Lst& rhs) {
	return std::tie(lhs.gid, lhs.cid, lhs.pl, lhs.rawid) < std::tie(rhs.gid, rhs.cid, rhs.pl, rhs.rawid);
}


void add_vertex_pl(std::vector<Momentum_recon::Event_information>& momch, std::vector<Momentum_recon::Event_information>& out_momch);

int main(int argc, char** argv) {
	if (argc != 3) {
		fprintf(stderr, "usage\n");
		exit(1);
	}
	std::string file_in_momch_ori = argv[1];//input momch
	std::string file_out_momch = argv[2];//erase list

	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch_ori);
	std::vector<Momentum_recon::Event_information> out_momch;
	add_vertex_pl(momch, out_momch);

	Momentum_recon::Write_Event_information_extension(file_out_momch, out_momch);
}
void add_vertex_pl(std::vector<Momentum_recon::Event_information>& momch, std::vector<Momentum_recon::Event_information>& out_momch) {
	
	Momentum_recon::Mom_chain chain;
	for (auto& ev : momch) {

		for (auto& c : ev.chains) {
			if (c.chainid != 0)continue;
			
			ev.vertex_pl = c.base.rbegin()->pl;
			
		}
	}

	out_momch = momch;
}