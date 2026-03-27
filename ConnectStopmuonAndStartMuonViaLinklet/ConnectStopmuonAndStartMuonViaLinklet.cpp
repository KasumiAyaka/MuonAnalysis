//2024/11/13 kasumi
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

struct Key {
	int pl, rawid;
};
bool operator<(const Key& lhs, const Key& rhs) {
	return std::tie(lhs.pl, lhs.rawid) < std::tie(rhs.pl, rhs.rawid);
}

struct Basetrack {
	Key k;
	double ax, ay, x, y;
};

struct Linklet {
	Basetrack b0, b1;
};
struct StartTrack {
	int gid, nseg, cid;
	Basetrack b0, b1;
	int ph0, ph1;
	double z0, z1;
};

struct Output {
	Key k;
	StartTrack st;
};
bool operator<(const Output& lhs, const Output& rhs) {
	return std::tie(lhs.k, lhs.st.gid, lhs.st.cid,lhs.st.nseg, lhs.st.b0.k, lhs.st.b1.k, lhs.st.ph0, 
		lhs.st.ph1, lhs.st.b0.ay, lhs.st.b0.ax, lhs.st.b0.x, lhs.st.b0.y, lhs.st.z0, lhs.st.z1)
		< std::tie(rhs.k, rhs.st.gid, rhs.st.cid, rhs.st.nseg, rhs.st.b0.k, rhs.st.b1.k, rhs.st.ph0,
			rhs.st.ph1, rhs.st.b0.ay, rhs.st.b0.ax, rhs.st.b0.x, rhs.st.b0.y, rhs.st.z0, rhs.st.z1);
}

void set_corrmap_Linklet(std::string input, std::multimap<Key, Linklet>& map);
void set_corrmap_StartTrack(std::string input, std::multimap<Key, StartTrack>& map);
void culc_difference(std::multimap<Key, Linklet>& map_new, std::multimap<Key, StartTrack>& map_old, std::string output);

int main(int argc, char** argv) {
	if (argc != 4) {
		fprintf(stderr, "usage : input-linklet.txt input-starttrack.txt output.txt \n");
		exit(1);
	}
	std::string in_old = argv[1];//input momch
	std::string in_new = argv[2];//multi
	std::string output = argv[3];//uniqe

	std::multimap<Key, Linklet> lnk;
	std::multimap<Key, StartTrack> stk;
	set_corrmap_Linklet(in_old, lnk);
	set_corrmap_StartTrack(in_new, stk);
	culc_difference(lnk, stk, output);


}

void set_corrmap_Linklet(std::string input, std::multimap<Key, Linklet> &map)
{
	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "file open error" << std::endl;
		exit(1);
	}

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;

	std::cout << input << std::endl;
	int count = 0;
	
	Linklet l;
	while (std::getline(ifs, str)) {

		str_v = StringSplit_with_tab(str);
		// downstream
		l.b0.k.pl = std::stoi(str_v[0]);
		l.b0.k.rawid = std::stoi(str_v[1]);
		l.b0.ax = std::stod(str_v[2]);
		l.b0.ay = std::stod(str_v[3]);
		l.b0.x = std::stod(str_v[4]);
		l.b0.y = std::stod(str_v[5]);
		// upstream
		l.b1.k.pl = std::stoi(str_v[6]);
		l.b1.k.rawid = std::stoi(str_v[7]);
		l.b1.ax = std::stod(str_v[8]);
		l.b1.ay = std::stod(str_v[9]);
		l.b1.x = std::stod(str_v[10]);
		l.b1.y = std::stod(str_v[11]);

		map.insert(std::make_pair(l.b1.k,l));
		count++;
	}
	std::cout << "\t# of Linklet : " << count << std::endl;

}

void set_corrmap_StartTrack(std::string input, std::multimap<Key, StartTrack>& map)
{
	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "file open error" << std::endl;
		exit(1);
	}

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;

	std::cout << input << std::endl;
	int count = 0;
	StartTrack s;
	while (std::getline(ifs, str)) {

		str_v = StringSplit_with_tab(str);
		s.cid = std::stoi(str_v[0]);
		s.nseg = std::stoi(str_v[1]);
		s.gid = std::stoi(str_v[2]);
		// downstream
		s.b0.k.pl = std::stoi(str_v[3]);
		s.b0.k.rawid = std::stoi(str_v[4]);
		s.ph0 = std::stod(str_v[5]);
		s.b0.ax = std::stod(str_v[6]);
		s.b0.ay = std::stod(str_v[7]);
		s.b0.x = std::stod(str_v[8]);
		s.b0.y = std::stod(str_v[9]);
		s.z0 = std::stod(str_v[10]);
		// upstream
		s.b1.k.pl = std::stoi(str_v[11]);
		s.b1.k.rawid = std::stoi(str_v[12]);
		s.ph1 = std::stoi(str_v[13]);
		s.b1.ax = std::stod(str_v[14]);
		s.b1.ay = std::stod(str_v[15]);
		s.b1.x = std::stod(str_v[16]);
		s.b1.y = std::stod(str_v[17]);
		s.z1 = std::stod(str_v[18]);

		map.insert(std::make_pair(s.b0.k, s));
		count++;
	}
	std::cout << "\t# of Starttrack : " << count << std::endl;

}

void culc_difference(std::multimap<Key, Linklet>& linklet, std::multimap<Key, StartTrack>& start, std::string output)
{
	std::ofstream ofs(output);

	Output o;
	std::set<Output> set;
	int count = 0;
	int cnd = 0;
	for (auto itr = linklet.begin(); itr != linklet.end(); itr++) {//Linklet

		auto p = start.equal_range(itr->first);
		if (p.first == start.end()) {
			//count++;
			std::cout << "not found..." << std::endl;
			++cnd;
		}
		else {
			count++;
			for (auto it = p.first; it != p.second; it++) {
				//std::cout << "found!" << std::endl;
				o.k = itr->second.b0.k;
				o.st = it->second;
				set.insert(o);
			}
		}
	}

	std::cout << "\nMatching\tLinklet total : " << linklet.size() << std::endl;//it should be count+cnd=linklet.size()
	std::cout << "\t# of Matched : " << count << ";\t# of not found : " << cnd << std::endl;
	std::cout<<"\tcheck : " << set.size() << std::endl;

	linklet.clear();
	start.clear();

	for (auto itr = set.begin(); itr != set.end(); itr++) {
		ofs << std::fixed << std::right
			<< std::setw(10) << std::setprecision(1) << itr->k.pl << " "//pl
			<< std::setw(10) << std::setprecision(1) << itr->k.rawid << " "//pos
			<< std::setw(10) << std::setprecision(1) << itr->st.gid << " "//start

			<< std::setw(10) << std::setprecision(1) << itr->st.nseg << " "
			<< std::setw(10) << std::setprecision(1) << itr->st.cid << " "
			<< std::setw(10) << std::setprecision(1) << itr->st.b0.k.pl << " "
			<< std::setw(10) << std::setprecision(1) << itr->st.b0.k.rawid << " "
			<< std::setw(10) << std::setprecision(1) << itr->st.ph0 << " "
			<< std::setw(10) << std::setprecision(1) << itr->st.b0.ax << " "
			<< std::setw(10) << std::setprecision(1) << itr->st.b0.ay << " "
			<< std::setw(10) << std::setprecision(1) << itr->st.b0.x << " "
			<< std::setw(10) << std::setprecision(1) << itr->st.b0.y << " "
			<< std::setw(10) << std::setprecision(1) << itr->st.z0 << " "

			<< std::setw(10) << std::setprecision(1) << itr->st.nseg << " "
			<< std::setw(10) << std::setprecision(1) << itr->st.cid << " "
			<< std::setw(10) << std::setprecision(1) << itr->st.b1.k.pl << " "
			<< std::setw(10) << std::setprecision(1) << itr->st.b1.k.rawid << " "
			<< std::setw(10) << std::setprecision(1) << itr->st.ph1 << " "
			<< std::setw(10) << std::setprecision(1) << itr->st.b1.ax << " "
			<< std::setw(10) << std::setprecision(1) << itr->st.b1.ay << " "
			<< std::setw(10) << std::setprecision(1) << itr->st.b1.x << " "
			<< std::setw(10) << std::setprecision(1) << itr->st.b1.y << " "
			<< std::setw(10) << std::setprecision(1) << itr->st.z1 << " "
			<< std::endl;
		count++;
	}
}