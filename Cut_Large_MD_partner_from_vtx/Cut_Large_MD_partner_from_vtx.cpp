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
	Key k0,k1;
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


void set_VTXfile(std::string input, std::string output, double md_thr, double VPH_thr);
void matching(std::multimap<int, Btrk>& trk, std::set<CutList>& cut_pair, std::map<int, Key>& mu, std::string output, double VPH_thr);

int main(int argc, char** argv) {
	if (argc < 4 || argc > 5) {
		fprintf(stderr, "usage1:in-vtx.txt out-list.txt MD-thresold\n");
		fprintf(stderr, "usage2:in-vtx.txt out-list.txt MD-thresold VPH-thr\n vph < VPH-thr => cut\n");
		exit(1);
	}
	std::string input = argv[1];
	std::string output = argv[2];
	double md_thr = std::stod(argv[3]);
	double VPH_thr = 2000.0001;
	if (argc == 5) {
		VPH_thr = std::stod(argv[4]);
	}
	set_VTXfile(input, output, md_thr, VPH_thr);

	std::cout << "Finished." << std::endl;

}

void set_VTXfile(std::string input, std::string output, double md_thr,double VPH_thr) {
	std::ifstream ifs(input);

	if (!ifs) {
		std::cerr << "File open error!" << std::endl;
		exit(0);
	}

	std::multimap<int, Btrk> trk;
	//std::multimap<int,Pair> cut_pair;
	std::set<CutList> cut_pair;
	std::map<int, Key> mu;

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;
	

	std::cout << input << " ; MD > " << md_thr << "‚ğcut."<< std::endl;
	int count = 0;
	CutList k;
	Btrk b;
	int gid=0;
	int vpl = 0;
	while (std::getline(ifs, str)) {

		str_v = StringSplit_with_tab(str);


		if (str_v.size() == 6) {
			// header
			gid = std::stoi(str_v[0]);
			vpl= std::stoi(str_v[1]);
		}
		else if (str_v.size() == 9) {
			// infomation of pair
			if (std::stod(str_v[8]) > md_thr) {
				k.k0.pl = std::stoi(str_v[0]);
				k.k0.rawid = std::stoi(str_v[1]);
				k.k1.pl = std::stoi(str_v[2]);
				k.k1.rawid = std::stoi(str_v[3]);
				k.gid = gid;
				cut_pair.insert(k);
				// Debag ECC1 event 321
				//if (k.k0.rawid == 1399727) {
				//	std::cout << k.gid << " " << std::stod(str_v[8]) << " : " << k.k1.pl << " " << k.k1.rawid << std::endl;
				//}else if(k.k1.rawid == 1399727) {
				//	std::cout << k.gid << " " << std::stod(str_v[8]) << " : " << k.k0.pl << " " << k.k0.rawid << std::endl;
				//}
			}
		}
		else if (str_v.size() == 13) {
			// information of basetrack
			if (std::stoi(str_v[1]) <= vpl) {
				b.pl0 = std::stoi(str_v[0]);
				b.k.pl = std::stoi(str_v[1]);

			}
			else {
				b.pl0 = std::stoi(str_v[1]);
				b.k.pl = std::stoi(str_v[0]);
			}

			b.k.rawid = std::stoi(str_v[2]);
			b.cid = std::stoi(str_v[3]);
			b.gid = std::stoi(str_v[4]);
			b.nseg = std::stoi(str_v[5]);
			b.npl = std::stoi(str_v[6]);
			b.ph = std::stoi(str_v[7]);
			b.ax = std::stod(str_v[8]);
			b.ay = std::stod(str_v[9]);
			b.y = std::stod(str_v[10]);
			b.x = std::stod(str_v[11]);
			b.z = std::stod(str_v[12]);
			trk.insert(std::make_pair(b.gid, b));
			if (b.cid == 0) {
				mu.insert(std::make_pair(b.gid,b.k));
			}
		}
		else {
			// unexpected error
			std::cout << str_v.size() << std::endl;
			exit(1);
		}
	}

	matching(trk, cut_pair, mu, output, VPH_thr);

}
void matching_old(std::multimap<int, Btrk>& trk, std::set<CutList>& cut_pair, std::map<int, Key>& mu,std::string output) {

	std::set<int> set;
	std::set<List> lst;//gid,pl,rid
	Pair p;
	List k;
	for (auto itr = cut_pair.begin(); itr != cut_pair.end(); itr++) {
		p.k0 = itr->k0;
		p.k1 = itr->k1;
		k.k = mu.find(itr->gid)->second;
		k.gid = itr->gid;

		if (p.k0.pl == k.k.pl && p.k0.rawid == k.k.rawid) {
			//muon
			set.insert(itr->gid);
			k.k = p.k1;
			lst.insert(k);

		}
		else if (p.k1.pl == k.k.pl && p.k1.rawid == k.k.rawid) {
			//muon
			set.insert(itr->gid);
			k.k = p.k0;
			lst.insert(k);
		}
		else {
			// pair of others

		}

	}
	cut_pair.clear();
	mu.clear();

	std::ofstream ofs(output);
	for (auto itr0 = set.begin(); itr0 != set.end(); itr0++) {
		// gid

		auto list = trk.equal_range(*itr0);
		for (auto itr = list.first; itr != list.second; itr++) {
			// loop by cid

			k.k = itr->second.k;
			k.gid = itr->first;

			if (lst.find(k) != lst.end()) {
				if (itr->second.cid != 0) {
					// except muon track
					//std::cout << k.gid << " " << k.k.pl << " " << k.k.rawid << std::endl;
						ofs << std::right << std::fixed
							<< std::setw(5) << std::setprecision(0) << k.gid << " "
							<< std::setw(3) << std::setprecision(0) << itr->second.cid << std::endl;
				}
			}
		}
	}
}

void matching(std::multimap<int, Btrk>& trk, std::set<CutList>& cut_pair, std::map<int, Key>& mu, std::string output,double VPH_thr) {

	//std::multimap<int, Btrk> & trk : gid, basetrack
	//std::set<CutList> & cut_pair : pair
	//std::map<int, Key> & mu : gid, muon-(pl,rawid)

	std::set<int> set;// gid
	std::set<List> lst;//gid,pl,rid
	Pair p;
	List k;
	for (auto itr = cut_pair.begin(); itr != cut_pair.end(); itr++) {
		p.k0 = itr->k0;
		p.k1 = itr->k1;
		
		k.gid = itr->gid;
		k.k = mu.find(itr->gid)->second;

		if (p.k0.pl == k.k.pl && p.k0.rawid == k.k.rawid) {
			// k0 -> muon
			set.insert(itr->gid);
			k.k = p.k1;
			lst.insert(k);

		}
		else if (p.k1.pl == k.k.pl && p.k1.rawid == k.k.rawid) {
			// k1 -> muon
			set.insert(itr->gid);
			k.k = p.k0;
			lst.insert(k);
		}
		else {
			// pair of others

		}

	}
	cut_pair.clear();
	mu.clear();

	std::cout << "VPH < " << VPH_thr << " ‚ğcut." << std::endl;
	std::ofstream ofs(output);
	for (auto itr0 = set.begin(); itr0 != set.end(); itr0++) {
		// gid

		auto list = trk.equal_range(*itr0);
		for (auto itr = list.first; itr != list.second; itr++) {
			// loop by cid

			k.gid = itr->first;
			k.k = itr->second.k;

			//std::cout << k.gid << " " << itr->second.cid << " " << k.k.pl << " " << k.k.rawid << std::endl;
			
			if (lst.find(k) != lst.end()) {
				if (itr->second.cid != 0) {
					// except muon track

					if (itr->second.ph % 10000 < VPH_thr) {// 2025/05/12 except proton
						std::cout << "cut : " << k.gid << " " << itr->second.cid << "\t" << k.k.pl << " " << k.k.rawid << " " << itr->second.ph % 10000 << std::endl;
						ofs << std::right << std::fixed
							<< std::setw(5) << std::setprecision(0) << k.gid << " "
							<< std::setw(3) << std::setprecision(0) << itr->second.cid << std::endl;
					}
				}
			}
		}
	}
}
