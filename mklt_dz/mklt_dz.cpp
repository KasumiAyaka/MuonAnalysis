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
bool operator<(const Key& lhs, const Key& rhs) {
	return std::tie(lhs.pl, lhs.rawid) < std::tie(rhs.pl, rhs.rawid);
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


void set_list(std::string input, std::set<Key>& eventlist);
void divide_event(std::string input, std::string output,std::set<Key>& eventlist);
int main(int argc, char** argv) {
	if (argc < 4) {
		fprintf(stderr, "usage1:input.txt event_list.txt output1.txt(listed_events) output2.txt");
		exit(1);
	}
	std::string input1 = argv[1];
	std::string input2 = argv[2];
	std::string output1 = argv[3];

	std::set<Key> list;
	set_list(input2, list);
	divide_event(input1, output1, list);
	std::cout << "Finished." << std::endl;

}
void set_list(std::string input, std::set<Key>& eventlist) {

	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "File open error!" << std::endl;
		exit(0);
	}
	Key k;
	while (ifs >>k.pl>>k.rawid) {
		//#ecc #event
		eventlist.insert(k);

	}
}
void divide_event(std::string input, std::string output,std::set<Key>& eventlist) {
	std::ofstream ofs(output);

	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "File open error!" << std::endl;
		exit(0);
	}

	int cnt = 0;
	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;
	Key k;

	while (std::getline(ifs, str)) {

		str_v = StringSplit_with_tab(str);

		k.pl = std::stoi(str_v[0]);
		k.rawid = std::stoi(str_v[1]);

		auto itr = eventlist.find(k);
		if (itr != eventlist.end()) {
			ofs << str << std::endl;
		}
		cnt++;
	}
	std::cout << cnt << std::endl;

}
