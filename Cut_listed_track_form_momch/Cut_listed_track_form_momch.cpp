// 2024/05/09
// kasumi
// based on "Check_upstream_base.cpp

#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>
#include <iomanip>

class output_format {
public:
	int groupid, chainid, pl, dpl;
	double sigma[2], chi2[4], md, dz;
};


void Cut_list(std::string input, std::multimap<int, int>& list);
void cut_upstream_base(Momentum_recon::Mom_chain& c, int pl);
void cut_upstream(std::vector<Momentum_recon::Event_information>& mom, std::multimap<int, int>& list);

int main(int argc, char** argv) {
	if (argc != 4) {
		fprintf(stderr, "usage:file-in-momch file-out-momch in_list\n");
		fprintf(stderr, "usage:file-in-momch.momch file-out-momch.txt file-out-momch.txt out_momch.momch\n");
		exit(1);
	}
	std::string file_in_momch = argv[1];
	std::string in_list = argv[3];
	std::string file_cut_momch = argv[2];

	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch);

	std::multimap<int, int> list;
	Cut_list(in_list, list);

	cut_upstream(momch, list);//chi2‚ªthr‚æ‚è‚‚¢‚à‚Ì‚Ííœ
	Momentum_recon::Write_Event_information_extension(file_cut_momch, momch);
}
void Cut_list(std::string input,std::multimap<int,int> &list){
	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "file open error" << std::endl;
		exit(1);
	}

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;

	int pl, gid;
	std::cout << input << std::endl;
	int count = 0;
	while (std::getline(ifs, str)) {

		str_v = StringSplit_with_tab(str);

		gid = std::stoi(str_v[0]);
		pl = std::stoi(str_v[1]);
		list.insert(std::make_pair(gid,pl));
		count++;
	}
	std::cout << "erase : " << count << std::endl;
	for (auto itr = list.begin(); itr != list.end(); itr++) {
		std::cout << itr->first << " " << itr->second << std::endl;
	}
}

void cut_upstream(std::vector<Momentum_recon::Event_information>& mom, std::multimap<int, int>& list) {
	std::vector<Momentum_recon::Mom_chain> ret;

	int pl=0;
	for (auto& ev : mom) {
		for (auto& c : ev.chains){

			if (list.find(ev.groupid) == list.end())continue;
			auto lst = list.equal_range(ev.groupid);

			for (auto itr = lst.first; itr != lst.second; itr++) {
				if (pl == 0) {
					pl = itr->second;
				}
				else if (pl > itr->second) {
					pl = itr->second;
				}
			}
			cut_upstream_base(c, pl);
			printf("eventid %d, cut PL %d\n", ev.groupid,pl);

		}
		pl = 0;
	}
}
void cut_upstream_base(Momentum_recon::Mom_chain& c, int pl) {
	for (auto itr = c.base.begin(); itr != c.base.end();) {
		if (itr->pl > pl) {
			itr = c.base.erase(itr);
		}
		else {
			itr++;
		}
	}
	for (auto itr = c.base_pair.begin(); itr != c.base_pair.end();) {
		if (itr->first.pl >= pl) {
			itr = c.base_pair.erase(itr);
		}
		else {
			itr++;
		}
	}
}