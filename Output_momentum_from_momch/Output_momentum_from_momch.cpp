#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

struct Lst {
	int gid, cid;
};
bool operator<(const Lst& lhs, const Lst& rhs) {
	return std::tie(lhs.gid, lhs.cid) < std::tie(rhs.gid, rhs.cid);
}


void Write(std::vector<Momentum_recon::Event_information>& in_momch, std::string output);

int main(int argc, char** argv) {
	if (argc != 3) {
		fprintf(stderr, "usage:input.momch out.txt");

		exit(1);
	}
	std::string file_in_momch = argv[1];//input momch
	std::string file_out_momch = argv[2];//output

	std::vector<Momentum_recon::Event_information> in_momch = Momentum_recon::Read_Event_information_extension(file_in_momch);

	Write(in_momch, file_out_momch);
}

void Write(std::vector<Momentum_recon::Event_information>& in_momch, std::string output) {

	std::ofstream ofs(output);
	int cnt = 0;
	int flg, mip;
	for (auto& ev : in_momch) {
		flg = 0;
		for (auto& c : ev.chains) {
			if (c.particle_flg == 13) {
				//muon
				ofs<< std::right << std::fixed
					<< std::setw(5) << std::setprecision(0) << ev.groupid << " "
					<< std::setw(5) << std::setprecision(0) << c.particle_flg << " "
					<< std::setw(10) << std::setprecision(1) << c.ecc_mcs_mom[0] << " "
					<< std::setw(10) << std::setprecision(1) << c.ecc_mcs_mom_error[0][0] << " "
					<< std::setw(10) << std::setprecision(1) << c.ecc_mcs_mom_error[0][1] << " "
					<< std::setw(10) << std::setprecision(1) << c.bm_range_mom << " "
					<< std::setw(10) << std::setprecision(1) << c.bm_range_mom_error[0] << " "
					<< std::setw(10) << std::setprecision(1) << c.bm_range_mom_error[1] << " "
					<< std::endl;
			}
			else if (c.particle_flg == 2212) {
				//proton 
				flg++;
			}
			else if (c.particle_flg == 211) {
				//
				mip++;
			}
			else {
				mip++;
			}
		}

	}


}