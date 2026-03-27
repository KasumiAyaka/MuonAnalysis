// 2024/06/14
//kasumi
//I:\NINJA\E71a\work\kasumi\ECC\MuonAnalysis\x64\Release\erase_cahin_from_momch.exe  event_water_fin3.momch ..\temp_vtx_info\elaselist.txt event_water_fin3_ECC6.momch
//2024/07/12  list をgid cidに変更
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

struct Lst {
	int gid, cid;
};
bool operator<(const Lst& lhs, const Lst& rhs) {
	return std::tie(lhs.gid, lhs.cid) < std::tie(rhs.gid, rhs.cid);
}


void erase_momch(std::vector<Momentum_recon::Event_information>& momch);

int main(int argc, char** argv) {
	if (argc != 3) {
		fprintf(stderr, "usage: in.momch out.momch\n");
		exit(1);
	}
	std::string file_in_momch_ori = argv[1];//input momch
	std::string file_out_momch = argv[2];//erase list

	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch_ori);
	erase_momch(momch);

	Momentum_recon::Write_Event_information_extension(file_out_momch, momch);
}

void erase_momch(std::vector<Momentum_recon::Event_information>& momch) {
	//std::vector<Momentum_recon::Event_information> mom;
	//Momentum_recon::Mom_chain chains;

	for (auto& ev : momch) {
		for (auto itr = ev.chains.begin(); itr != ev.chains.end();) {
			// 条件一致した要素を削除する
			if (itr->chainid!=0) {
				// 削除された要素の次を指すイテレータが返される。
				itr = ev.chains.erase(itr);
			}
			// 要素削除をしない場合に、イテレータを進める
			else {
				++itr;
			}
			//std::cout << ev.chains.size() << std::endl;
			//std::cout<<std::endl;
		}
	}
}