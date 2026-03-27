#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <utility>
#include <iomanip>
#include <set>
#include <sstream>

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

struct mom {
	int ev, pid;
	double mcs, mceerr[2], bm, bmerr[2];
	int utime, bunch;
};

void read_list(std::string in, std::map<int, mom>& map);
void read_sf(std::string in, std::map<int, mom>& map, std::string out, int mode);

int main(int argc, char** argv) {
	if (argc < 4) {
		fprintf(stderr, "************************************************************");
		fprintf(stderr, "\tusage");
		fprintf(stderr, " 1  :in1.txt sf.txt out.txt");
		fprintf(stderr, " 2  :in1.txt sf.txt out.txt mode");
		fprintf(stderr, "************************************************************");
		fprintf(stderr, "\tmode");
		fprintf(stderr, " none or 0 : eid pid mcs -err +err BM -err +err");
		fprintf(stderr, " 1         : utime eid mcs -err +err BM -err +err");
		fprintf(stderr, " 2         : utime bunch eid pid mcs -err +err BM -err +err");
		fprintf(stderr, "************************************************************");
		exit(1);
	}
	std::string input = argv[1];
	std::string input2 = argv[2];
	std::string output = argv[3];
	int n = 0;
	if (argc != 4) {
		n = std::stoi(argv[4]);
	}

	std::map<int, mom> map;
	read_list(input, map);
	read_sf(input2, map, output,n);

	std::cout << "Finished." << std::endl;

}

void read_list(std::string in, std::map<int, mom>& map) {
	std::ifstream ifs(in);
	if (!ifs) {
		std::cerr << " ! " << std::endl;
		return ;
	}
	mom m;
	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;
	while (std::getline(ifs, str)) {

		str_v = StringSplit_with_tab(str);

		m.ev = std::stoi(str_v[0]);
		m.pid = std::stoi(str_v[1]);
		m.mcs = std::stod(str_v[2]);
		m.mceerr[0] = std::stod(str_v[3]);
		m.mceerr[1] = std::stod(str_v[4]);
		m.bm = std::stod(str_v[5]);
		m.bmerr[0] = std::stod(str_v[6]);
		m.bmerr[1] = std::stod(str_v[7]);

		m.bunch = 0;
		m.utime = 0;
		map.insert(std::make_pair(m.ev, m));
	}

}

void read_sf(std::string in,std::map<int,mom>& map,std::string out,int mode) {
	std::ifstream ifs(in);
	if (!ifs) {
		std::cerr << " ! " << std::endl;
		return;
	}
	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;
	int eid;
	double mom;
	while (std::getline(ifs, str)) {

		str_v = StringSplit_with_tab(str);
		auto itr = map.find(std::stoi(str_v[21]));

		if (itr != map.end()&&itr->second.bm==-1.0 &&itr->second.bmerr[0]==-1.0) {
			itr->second.bm = std::stoi(str_v[16]);
			itr->second.bunch = std::stoi(str_v[12]);
			itr->second.utime = std::stoi(str_v[10]);
		}
		
	}


	std::ofstream ofs(out);
	if(mode==0){
		for (auto itr = map.begin(); itr != map.end(); itr++) {
			ofs << std::right << std::fixed
				<< std::setw(5) << std::setprecision(0) << itr->second.ev << " "
				<< std::setw(5) << std::setprecision(0) << itr->second.pid << " "
				<< std::setw(10) << std::setprecision(1) << itr->second.mcs << " "
				<< std::setw(10) << std::setprecision(1) << itr->second.mceerr[0] << " "
				<< std::setw(10) << std::setprecision(1) << itr->second.mceerr[1] << " "
				<< std::setw(10) << std::setprecision(1) << itr->second.bm << " "
				<< std::setw(10) << std::setprecision(1) << itr->second.bmerr[0] << " "
				<< std::setw(10) << std::setprecision(1) << itr->second.bmerr[1] << " "
				<< std::endl;
		}

	}
	else if(mode==1){
		for (auto itr = map.begin(); itr != map.end(); itr++) {
			ofs << std::right << std::fixed
				<< std::setw(10) << std::setprecision(0) << itr->second.utime << " "
				<< std::setw(5) << std::setprecision(0) << itr->second.ev << " "
				<< std::setw(10) << std::setprecision(1) << itr->second.mcs << " "
				<< std::setw(10) << std::setprecision(1) << itr->second.mceerr[0] << " "
				<< std::setw(10) << std::setprecision(1) << itr->second.mceerr[1] << " "
				<< std::setw(10) << std::setprecision(1) << itr->second.bm << " "
				<< std::setw(10) << std::setprecision(1) << itr->second.bmerr[0] << " "
				<< std::setw(10) << std::setprecision(1) << itr->second.bmerr[1] << " "
				<< std::endl;


		}
	}
	else {
		for (auto itr = map.begin(); itr != map.end(); itr++) {
			ofs << std::right << std::fixed
				<< std::setw(10) << std::setprecision(0) << itr->second.utime << " "
				<< std::setw(1) << std::setprecision(0) << itr->second.bunch << " "
				<< std::setw(5) << std::setprecision(0) << itr->second.ev << " "
				<< std::setw(5) << std::setprecision(0) << itr->second.pid << " "
				<< std::setw(10) << std::setprecision(1) << itr->second.mcs << " "
				<< std::setw(10) << std::setprecision(1) << itr->second.mceerr[0] << " "
				<< std::setw(10) << std::setprecision(1) << itr->second.mceerr[1] << " "
				<< std::setw(10) << std::setprecision(1) << itr->second.bm << " "
				<< std::setw(10) << std::setprecision(1) << itr->second.bmerr[0] << " "
				<< std::setw(10) << std::setprecision(1) << itr->second.bmerr[1] << " "
				<< std::endl;
		}
	}
}