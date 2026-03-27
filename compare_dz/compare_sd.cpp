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

struct Corrmap {
	//  3   7475.5   5018.0  -2980.0   9475.5   3018.0  -2980.0
	int pl;
	double x[2], y[2], z[2];
};

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


void set_list(std::string input, std::map<int, double>& list);
void Set_vertex(std::string input, std::string output, std::map<int, double>& list);

int main(int argc, char** argv) {
	if (argc < 4) {
		fprintf(stderr, "vtx.txt list.txt output.txt");
		exit(1);
	}
	std::string input = argv[1];
	std::string input2 = argv[2];
	std::string output = argv[3];

	std::map<int, double>list;
	set_list(input2, list);
	Set_vertex(input, output, list);
	std::cout << "Finished." << std::endl;

}
void set_list(std::string input, std::map<int,double>& list) {

	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "File open error!" << std::endl;
		exit(0);
	}
	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;

	int pl;
	double  z;
	while (std::getline(ifs, str)) {

		str_v = StringSplit_with_tab(str);
		//std::cout << str << std::endl;

		pl = std::stoi(str_v[0]);
		z = std::stod(str_v[1]);
		list.insert(std::make_pair(pl, z));
	}
}
void Set_vertex(std::string input, std::string output, std::map<int, double>& list) {
	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "File open error!" << std::endl;
		exit(0);
	}

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;

	std::ofstream ofs(output);
	VTX_header v;
	double dz;
	while (std::getline(ifs, str)) {

		str_v = StringSplit_with_tab(str);
		//std::cout << str << std::endl;
		if (str_v.size() == 6) {
			//header
			v.gid = std::stoi(str_v[0]);
			v.pl = std::stoi(str_v[1]);
			v.trknum = std::stoi(str_v[2]);
			v.x = std::stoi(str_v[3]);
			v.y = std::stoi(str_v[4]);
			v.z = std::stoi(str_v[5]);

			auto itr = list.find(v.pl);
			if (itr != list.end()) {

				dz = v.z - itr->second;
				std::cout << v.z << " " << itr->second << " " << dz << std::endl;

				ofs << std::right << std::fixed
					<< std::setw(5) << std::setprecision(0) << v.gid
					<< std::setw(5) << std::setprecision(0) << v.pl
					<< std::setw(5) << std::setprecision(0) << v.trknum
					<< std::setw(12) << std::setprecision(1) << v.x
					<< std::setw(12) << std::setprecision(1) << v.y
					<< std::setw(12) << std::setprecision(1) << dz
					<< std::endl;

			}
		}
		if (str_v.size() == 9) {
			//pair

		}
		if (str_v.size() == 13) {
			//track
		}

	}

}
