// Cut_Large_MD_partner_from_vtx
// usage1:md‚Ì‘å‚«‚¢partner‚ğœ‹
// usage2:md‚ª‘å‚«‚­‚Ä‚àVPH‚Ì‚‚¢track‚Íc‚·
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <utility>
#include <iomanip>
#include <set>
#include <sstream>


struct Key {
	int pl, rawid;
};
struct VTX_header {
	int gid, pl, trknum;
	double x, y, z;
};

struct Pair {
	Key k0, k1;
	double x, y, dz, oa, md;
};

struct Btrk {
	int pl0;
	Key k;
	int cid, gid, nseg, npl, ph;
	double ax, ay, x, y, z;
};

class Event {
	VTX_header v;
	std::vector<Pair>& p;
};

struct CutList {
	Key k0, k1;
	int gid;
};
bool operator<(const CutList& lhs, const CutList& rhs) {
	return std::tie(lhs.gid, lhs.k0.pl, lhs.k0.rawid, lhs.k1.pl, lhs.k1.rawid)
		< std::tie(rhs.gid, rhs.k0.pl, rhs.k0.rawid, rhs.k1.pl, rhs.k1.rawid);
}

struct List {
	Key k;
	int gid;
};
bool operator<(const List& lhs, const List& rhs) {
	return std::tie(lhs.gid, lhs.k.pl, lhs.k.rawid) < std::tie(rhs.gid, rhs.k.pl, rhs.k.rawid);
}


namespace {
	std::vector<std::string> StringSplit(std::string str) {
		std::stringstream ss{ str };
		std::vector<std::string> v;
		std::string buf;
		while (std::getline(ss, buf, ' ')) {
			if (buf != "") {
				v.push_back(buf);
			}
		}
		return v;
	}
	std::vector<std::string> StringSplit_with_tab(std::string str) {
		std::vector<std::string> v;

		std::vector<std::string> str_v = StringSplit(str);
		for (int i = 0; i < str_v.size(); i++) {
			std::stringstream ss{ str_v[i] };
			std::string buf;
			while (std::getline(ss, buf, '\t')) {
				if (buf != "") {
					v.push_back(buf);
				}
			}
		}
		return v;
	}
}


void set_list(std::string input, std::set<int>& eventlist);
void divide_event(std::string input, std::string output1, std::string output2, std::set<int>& eventlist);
int main(int argc, char** argv) {
	if (argc < 4) {
		fprintf(stderr, "usage1:input.txt event_list.txt output1.txt(listed_events) output2.txt");
		exit(1);
	}
	std::string input1 = argv[1];
	std::string input2 = argv[2];
	std::string output1 = argv[3];
	std::string output2 = argv[4];

	std::set<int> list;
	set_list(input2, list);
	divide_event(input1, output1, output2, list);
	std::cout << "Finished." << std::endl;

}
void set_list(std::string input, std::set<int>& eventlist) {

	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "File open error!" << std::endl;
		exit(0);
	}
	int a;
	while (ifs >> a) {
		eventlist.insert(a);
	}
}
void divide_event(std::string input,std::string output1, std::string output2, std::set<int>& eventlist) {
	std::ofstream ofs1(output1);
	std::ofstream ofs2(output2);
	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "File open error!" << std::endl;
		exit(0);
	}

	int gid = 0;
	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;
	int cnt0 = 0;
	int cnt1 = 0;

	while (std::getline(ifs, str)) {

		str_v = StringSplit_with_tab(str);
		gid = std::stoi(str_v[0]);
		//std::cout << str << std::endl;

		auto itr = eventlist.find(gid);
		
		if (itr != eventlist.end()) {
			//found
			ofs1 << str << std::endl;
			cnt0++;
		}
		else {
			ofs2 << str << std::endl;
			cnt1++;
		}

	}
	std::cout << cnt0 << "," << cnt1 << std::endl;

}
