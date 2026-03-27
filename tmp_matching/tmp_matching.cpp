#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <set>
#include <map>
#include <math.h>
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>


struct Key {
	int unix_time, daily, bunch;
};
bool operator<(const Key& lhs, const Key& rhs) {
	return std::tie(lhs.unix_time, lhs.daily, lhs.bunch) <
		std::tie(rhs.unix_time, rhs.daily, rhs.bunch);
}

struct nt {
	Key k;
	int ecc;
};

struct ntbm {
	Key k;
	int eid;
	double ay, ax, y, x;
};
void set_ntbm(std::string input, std::map<Key, ntbm>& map);
void set_NT_BMpred(std::string input, std::map<Key, nt>& map);
void matching(std::map<Key, ntbm>& map, std::map<Key, nt>& map2);

int main(int argc, char** argv) {
	if (argc == 4) {
		fprintf(stderr, "usage : input_sharingfile list.txt out-sf.txt \n");
		exit(1);
	}
	std::string input = argv[1];//input momch
	std::string in2 = argv[2];//multi
	std::string output = argv[3];//uniqe

	std::map<Key, ntbm> map;
	std::map<Key, nt> map2;
	set_ntbm(input,map);
	set_NT_BMpred(in2,map2);
	matching(map,map2);

}
void set_ntbm(std::string input, std::map<Key, ntbm>& map) {

	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "file open error" << std::endl;
		exit(1);
	}
	ntbm nt;

	while (ifs >> nt.k.unix_time >> nt.k.daily >> nt.k.bunch >> nt.ay >> nt.ax >> nt.y >> nt.x >> nt.eid) {
		map.insert(std::make_pair(nt.k, nt));
	}

}
void set_NT_BMpred(std::string input, std::map<Key, nt>& map) {

	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "file open error" << std::endl;
		exit(1);
	}
	nt n;

	while (ifs >> n.k.unix_time >> n.k.daily >> n.k.bunch >> n.ecc) {
		map.insert(std::make_pair(n.k, n));
	}

}
void matching(std::map<Key, ntbm>& map, std::map<Key, nt>& map2) {
	
	for (auto itr = map2.begin(); itr != map2.end(); itr++) {
		Key k;
		k = itr->first;
		k.unix_time = itr->first.unix_time-1;
		auto p = map.find(itr->first);
		if (p != map.end()) {
			std::cout << std::fixed << std::right
				<< std::setw(12) << std::setprecision(0) << itr->first.unix_time << " "
				<< std::setw(12) << std::setprecision(0) << itr->first.daily << " "
				<< std::setw(12) << std::setprecision(0) << itr->first.bunch << " "
				<< std::setw(12) << std::setprecision(3) << p->second.ay << " "
				<< std::setw(12) << std::setprecision(3) << p->second.ax << " "
				<< std::setw(12) << std::setprecision(1) << p->second.x << " "
				<< std::setw(12) << std::setprecision(1) << p->second.y << " "
				<< std::setw(12) << std::setprecision(0) << p->second.eid << " "
				<< std::setw(12) << std::setprecision(0) << itr->second.ecc
				<< std::endl;
		}
		else if(map.find(k)!=map.end()) {
			std::cout << std::fixed << std::right
				<< std::setw(12) << std::setprecision(0) << itr->first.unix_time << " "
				<< std::setw(12) << std::setprecision(0) << itr->first.daily << " "
				<< std::setw(12) << std::setprecision(0) << itr->first.bunch << " "
				<< std::setw(12) << std::setprecision(3) << p->second.ay << " "
				<< std::setw(12) << std::setprecision(3) << p->second.ax << " "
				<< std::setw(12) << std::setprecision(1) << p->second.x << " "
				<< std::setw(12) << std::setprecision(1) << p->second.y << " "
				<< std::setw(12) << std::setprecision(0) << p->second.eid << " "
				<< std::setw(12) << std::setprecision(0) << itr->second.ecc
				<< std::endl;
		}
	}
}

