
//unixtime/I:bunch/I:ecc/I:ey/D:ex/D:NTgy/D:NTgx/D:BMgy/D:BMgx/D:NTy/D:NTx/D:NTay/D:NTax/D:BMy/D:BMx/D:daily_in_file/I");
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <sstream>
#include <istream>

struct Key {
	int utime, bunch;
};
bool operator<(const Key& lhs, const Key& rhs) {
	return std::tie(lhs.utime, lhs.bunch) < std::tie(rhs.utime, rhs.bunch);
}
struct MuonPrediction {
	Key k;
	int ecc;
};

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

void set_list(std::string input, std::map<Key, int>& key);
void set_MuonPrediction(std::string input, std::multimap<Key, int>& map);
void Matcing(std::map<Key,int>& set, std::multimap<Key, int>& map, std::string output);

int main(int argc, char** argv) {
		if (argc != 4) {
			fprintf(stderr, "usage:Timestamplist.txt MuonPredLint.txt output\n");
			exit(1);
		}
		std::string list = argv[1];
		std::string MuonPred = argv[2];
		std::string output = argv[3];
		std::map<Key, int> key;
		std::multimap<Key, int> map;

		set_list(list, key);
		set_MuonPrediction(MuonPred, map);
		Matcing(key, map, output);
}
void set_list(std::string input, std::map<Key, int>&key) {
	std::ifstream ifs(input);

	if (!ifs) {
		std::cerr << "File open error!" << std::endl;
		exit(0);
	}

	Key k;
	int flg;
	while (ifs >> k.utime >> k.bunch>>flg) {
		key.insert(std::make_pair(k,flg));
	}
	std::cout << "list size = " << key.size() << std::endl;
}

void set_MuonPrediction(std::string input, std::multimap<Key,int>& map) {
	std::ifstream ifs(input);

	if (!ifs) {
		std::cerr << "File open error!" << std::endl;
		exit(0);
	}

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;

	std::cout << input << std::endl;
	int count = 0;
	Key k;
	while (std::getline(ifs, str)) {

		str_v = StringSplit_with_tab(str);
		k.utime = std::stoi(str_v[0]);
		k.bunch = std::stoi(str_v[1]);

		map.insert(std::make_pair(k, std::stoi(str_v[2])));
	}
	std::cout << "MuonPrediction size = " << map.size() << std::endl;

}
void Matcing(std::map<Key, int>& set, std::multimap<Key, int>& map,std::string output) {

	std::ofstream ofs(output);

	int cnt = 0;
	for (auto itr = set.begin(); itr != set.end(); itr++) {//timestamp,bunch,single=1/multi=2
		auto p = map.equal_range(itr->first);
		if (p.second == map.end())continue;

		for (auto q = p.first; q != p.second; q++) {//timestamp,bunch,ecc
			ofs << q->first.utime << " " << q->first.bunch << " " << itr->second << " " << q->second << std::endl;
			cnt++;
		}
	}
	std::cout << "matched muon prediction num = " << cnt << std::endl;
}