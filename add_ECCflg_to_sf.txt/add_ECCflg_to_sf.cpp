
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <sstream>
#include <istream>
#include <iomanip>

struct Key {
	int utime, bunch;
};
bool operator<(const Key& lhs, const Key& rhs) {
	return std::tie(lhs.utime, lhs.bunch) < std::tie(rhs.utime, rhs.bunch);
}

struct SharingFile {
	int ecc;
	int PL,
		ISS_rawid,
		OSS_rawid,
		fixed_rawid,
		tracker_rawid,
		shifter_spot,//spotA * 100 + spotB
		zone1,
		raw1,
		zone2,
		raw2,
		//unix_time,
		t_track_id,
		//BM_bunch,
		entry_in_daily_file,
		BM_nplane,
		charge;
	double	BM_mom,
		chi2_shifter[4];
	int	eventid,//(通し番号)
		PM_track_type,
		ECC_track_type;
	//momchとのマッチング
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

void set_list(std::string input, std::map<Key, int>& key, int ecc, int flg);
void set_sf(std::string input, std::multimap<Key, SharingFile>& map);
void Matcing(std::map<Key, int>& set, std::multimap<Key, SharingFile>& out);
void Write(std::multimap<Key, SharingFile>& map, std::string output);
void set_sf_normal_order(std::string input, std::multimap<Key, SharingFile>& map);
void Write_normal_order(std::multimap<Key, SharingFile>& map, std::string output);

int main(int argc, char** argv) {
	if (argc !=  7) {
		fprintf(stderr, "usage:sf.txt list-stop.txt list-pene.txt list-edge.txt #ECC output.txt\n");
		exit(1);
	}

	std::string list = argv[1];// eventid list
	std::string stop = argv[2]; //stop
	std::string pene = argv[3]; // pene
	std::string edge = argv[4];//edge
	int ecc = std::stoi(argv[5]);//ECC
	std::string output = argv[6];//edge

	std::map<Key, int> key;
	std::multimap<Key, SharingFile> out;

	set_list(stop, key, ecc, 1);
	set_list(pene, key, ecc, 2);
	set_list(edge, key, ecc, 3);
	set_sf_normal_order(list, out);
	Matcing(key, out);
	Write_normal_order(out, output);
}
void set_list(std::string input, std::map<Key, int>& key, int ecc, int flg) {
	std::ifstream ifs(input);

	if (!ifs) {
		std::cerr << "File open error!" << std::endl;
		exit(0);
	}

	Key k;
	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;

	std::cout << input << std::endl;
	int count = 0;
	int fg = flg;
	int j = 0;
	while (std::getline(ifs, str)) {
		str_v = StringSplit_with_tab(str);

		if (str_v.size() == 1) {
			std::cout << str << std::endl;
			k.utime = std::stoi(str_v[0]);//event
			k.bunch = ecc;
		}
		else if (str_v.size() == 3) {
			k.bunch = std::stoi(str_v[0]);//ecc
			k.utime = std::stoi(str_v[1]);//event
			fg = std::stoi(str_v[2]);//clasification:1=stop,2=edge,3=pene
			auto itr = key.find(k);
			if (key.find(k) != key.end()) {
				std::cout << k.bunch << " " << k.utime << " " << fg << std::endl;
				std::cout << itr->first.bunch << " " << itr->first.utime << " " << itr->second << std::endl;
				if (fg < itr->second) {
					itr->second = fg;
				}
				j++;
			}
		}
		else {
			std::cout << "unexpected error" << std::endl;
			std::cout << str_v.size() << std::endl;
		}

		key.insert(std::make_pair(k, fg));
	}

	std::cout << "\tSame ECC and EventID, different Classification = " << j << std::endl;
	std::cout << "\tlist size = " << key.size() << std::endl;
}
void set_sf(std::string input, std::multimap<Key, SharingFile>& map) {
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
	SharingFile s;
	while (std::getline(ifs, str)) {
		//std::cout << "ok?" << std::endl;
		str_v = StringSplit_with_tab(str);
		s.ecc = std::stoi(str_v[0]);
		k.utime = std::stoi(str_v[1]);
		k.bunch = std::stoi(str_v[2]);
		s.PL = std::stoi(str_v[3]);
		s.ISS_rawid = std::stoi(str_v[4]);
		s.OSS_rawid = std::stoi(str_v[5]);
		s.fixed_rawid = std::stoi(str_v[6]);
		s.tracker_rawid = std::stoi(str_v[7]);
		s.shifter_spot = std::stoi(str_v[8]);
		s.zone1 = std::stoi(str_v[9]);
		s.raw1 = std::stoi(str_v[10]);
		s.zone2 = std::stoi(str_v[11]);
		s.raw2 = std::stoi(str_v[12]);
		s.t_track_id = std::stoi(str_v[13]);
		s.entry_in_daily_file = std::stoi(str_v[14]);
		s.BM_nplane = std::stoi(str_v[15]);
		s.charge = std::stoi(str_v[16]);
		s.BM_mom = std::stod(str_v[17]);
		s.chi2_shifter[0] = std::stod(str_v[18]);
		s.chi2_shifter[1] = std::stod(str_v[19]);
		s.chi2_shifter[2] = std::stod(str_v[20]);
		s.chi2_shifter[3] = std::stod(str_v[21]);
		s.eventid = std::stoi(str_v[22]);
		s.PM_track_type = std::stoi(str_v[23]);
		s.ECC_track_type = std::stoi(str_v[24]);


		map.insert(std::make_pair(k, s));
	}
	std::cout << "\tMuonPrediction size = " << map.size() << std::endl;

}
void set_sf_normal_order(std::string input, std::multimap<Key, SharingFile>& map) {
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
	SharingFile s;
	while (std::getline(ifs, str)) {
		//std::cout << "ok?" << std::endl;
		str_v = StringSplit_with_tab(str);
		s.ecc = std::stoi(str_v[0]);
		s.PL = std::stoi(str_v[1]);
		s.ISS_rawid = std::stoi(str_v[2]);
		s.OSS_rawid = std::stoi(str_v[3]);
		s.fixed_rawid = std::stoi(str_v[4]);
		s.tracker_rawid = std::stoi(str_v[5]);
		s.shifter_spot = std::stoi(str_v[6]);
		s.zone1 = std::stoi(str_v[7]);
		s.raw1 = std::stoi(str_v[8]);
		s.zone2 = std::stoi(str_v[9]);
		s.raw2 = std::stoi(str_v[10]);
		k.utime = std::stoi(str_v[11]);
		s.t_track_id = std::stoi(str_v[12]);
		k.bunch = std::stoi(str_v[13]);
		s.entry_in_daily_file = std::stoi(str_v[14]);
		s.BM_nplane = std::stoi(str_v[15]);
		s.charge = std::stoi(str_v[16]);
		s.BM_mom = std::stod(str_v[17]);
		s.chi2_shifter[0] = std::stod(str_v[18]);
		s.chi2_shifter[1] = std::stod(str_v[19]);
		s.chi2_shifter[2] = std::stod(str_v[20]);
		s.chi2_shifter[3] = std::stod(str_v[21]);
		s.eventid = std::stoi(str_v[22]);
		s.PM_track_type = std::stoi(str_v[23]);
		s.ECC_track_type = std::stoi(str_v[24]);


		map.insert(std::make_pair(k, s));
	}
	std::cout << "\tMuonPrediction size = " << map.size() << std::endl;

}
void Matcing(std::map<Key, int>& set, std::multimap<Key, SharingFile>& out) {
	std::cout << "now matching..." << std::endl;
	int cnt = 0;
	Key k;
	for (auto itr = out.begin(); itr != out.end(); itr++) {
		k.bunch = itr->second.ecc;//ecc
		k.utime = itr->second.eventid;//eid
		//std::cout << k.bunch << " " << k.utime << std::endl;
		if (set.find(k) != set.end()) {
			int flg = set.find(k)->second;
			if (itr->second.ECC_track_type > 0) {
				std::cout << "!" << std::endl;
				itr->second.ECC_track_type = flg;
			}
			else {
				itr->second.ECC_track_type = flg;
			}
			cnt++;
		}
	}

	std::cout << "\tmatched muon prediction num = " << cnt << std::endl;
}
void Write(std::multimap<Key, SharingFile>& map, std::string output) {
	std::ofstream ofs(output);
	for (auto k = map.begin(); k != map.end(); k++) {
		ofs << std::fixed << std::right
			<< std::setw(2) << std::setprecision(0) << k->second.ecc << " "
			<< std::setw(10) << std::setprecision(0) << k->first.utime << " "
			<< std::setw(3) << std::setprecision(0) << k->first.bunch << " "
			<< std::setw(3) << std::setprecision(0) << k->second.PL << " "

			<< std::setw(10) << std::setprecision(0) << k->second.ISS_rawid << " "
			<< std::setw(10) << std::setprecision(0) << k->second.OSS_rawid << " "
			<< std::setw(10) << std::setprecision(0) << k->second.fixed_rawid << " "
			<< std::setw(10) << std::setprecision(0) << k->second.tracker_rawid << " "
			<< std::setw(5) << std::setprecision(0) << k->second.shifter_spot << " "

			<< std::setw(3) << std::setprecision(0) << k->second.zone1 << " "
			<< std::setw(10) << std::setprecision(0) << k->second.raw1 << " "
			<< std::setw(3) << std::setprecision(0) << k->second.zone2 << " "
			<< std::setw(10) << std::setprecision(0) << k->second.raw2 << " "

			<< std::setw(4) << std::setprecision(0) << k->second.t_track_id << " "
			<< std::setw(6) << std::setprecision(0) << k->second.entry_in_daily_file << " "
			<< std::setw(3) << std::setprecision(0) << k->second.BM_nplane << " "
			<< std::setw(3) << std::setprecision(0) << k->second.charge << " "
			<< std::setw(7) << std::setprecision(1) << k->second.BM_mom << " "

			<< std::setw(12) << std::setprecision(3) << k->second.chi2_shifter[0] << " "
			<< std::setw(12) << std::setprecision(3) << k->second.chi2_shifter[1] << " "
			<< std::setw(12) << std::setprecision(3) << k->second.chi2_shifter[2] << " "
			<< std::setw(12) << std::setprecision(3) << k->second.chi2_shifter[3] << " "
			<< std::setw(6) << std::setprecision(0) << k->second.eventid << " "
			<< std::setw(3) << std::setprecision(0) << k->second.PM_track_type << " "
			<< std::setw(3) << std::setprecision(0) << k->second.ECC_track_type
			<< std::endl;
	}

	std::cout << "finish writing " << output << std::endl;
}void Write_normal_order(std::multimap<Key, SharingFile>& map, std::string output) {
	std::ofstream ofs(output);
	for (auto k = map.begin(); k != map.end(); k++) {
		ofs << std::fixed << std::right
			<< std::setw(2) << std::setprecision(0) << k->second.ecc << " "
			<< std::setw(3) << std::setprecision(0) << k->second.PL << " "

			<< std::setw(10) << std::setprecision(0) << k->second.ISS_rawid << " "
			<< std::setw(10) << std::setprecision(0) << k->second.OSS_rawid << " "
			<< std::setw(10) << std::setprecision(0) << k->second.fixed_rawid << " "
			<< std::setw(10) << std::setprecision(0) << k->second.tracker_rawid << " "
			<< std::setw(5) << std::setprecision(0) << k->second.shifter_spot << " "

			<< std::setw(3) << std::setprecision(0) << k->second.zone1 << " "
			<< std::setw(10) << std::setprecision(0) << k->second.raw1 << " "
			<< std::setw(3) << std::setprecision(0) << k->second.zone2 << " "
			<< std::setw(10) << std::setprecision(0) << k->second.raw2 << " "
			<< std::setw(10) << std::setprecision(0) << k->first.utime << " "

			<< std::setw(4) << std::setprecision(0) << k->second.t_track_id << " "
			<< std::setw(3) << std::setprecision(0) << k->first.bunch << " "
			<< std::setw(6) << std::setprecision(0) << k->second.entry_in_daily_file << " "
			<< std::setw(3) << std::setprecision(0) << k->second.BM_nplane << " "
			<< std::setw(3) << std::setprecision(0) << k->second.charge << " "
			<< std::setw(7) << std::setprecision(1) << k->second.BM_mom << " "

			<< std::setw(12) << std::setprecision(3) << k->second.chi2_shifter[0] << " "
			<< std::setw(12) << std::setprecision(3) << k->second.chi2_shifter[1] << " "
			<< std::setw(12) << std::setprecision(3) << k->second.chi2_shifter[2] << " "
			<< std::setw(12) << std::setprecision(3) << k->second.chi2_shifter[3] << " "
			<< std::setw(6) << std::setprecision(0) << k->second.eventid << " "
			<< std::setw(3) << std::setprecision(0) << k->second.PM_track_type << " "
			<< std::setw(3) << std::setprecision(0) << k->second.ECC_track_type
			<< std::endl;
	}

	std::cout << "finish writing " << output << std::endl;
}