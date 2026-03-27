//2024/07/02 kasumi
//I:\NINJA\E71a\work\kasumi\ECC\MuonAnalysis\x64\Release\Output_duplicated_muon_in_SharingFile.exe R:\Analysis\NINJA\E71a\5\Shifter1\acryltrack\Timestamp\muon_ecc5_including_multi.txt ECC1_multi-sf.txt ECC1_unique-sf.txt

#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>


struct TrackerInfo {
	int unix_time, t_track_id, entry_in_daily_file;
	double	BM_mom;
};
bool operator<(const TrackerInfo& lhs, const TrackerInfo& rhs) {
	return std::tie(lhs.unix_time, lhs.t_track_id, lhs.entry_in_daily_file, lhs.BM_mom) <
		std::tie(rhs.unix_time, rhs.t_track_id, rhs.entry_in_daily_file, rhs.BM_mom);
}

struct SharingFile {
	TrackerInfo tk;
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
		//t_track_id,
		BM_bunch,
		//entry_in_daily_file,
		BM_nplane,
		charge;
	double	//BM_mom,
		chi2_shifter[4];
	int	eventid,//(通し番号)
		PM_track_type,
		ECC_track_type;
	//momchとのマッチング
	// unix_time,t_track_id,entry_in_daily_file
};
bool operator<(const SharingFile& lhs, const SharingFile& rhs) {
	return std::tie(lhs.tk, lhs.eventid, lhs.PL, lhs.ISS_rawid, lhs.OSS_rawid, lhs.fixed_rawid, lhs.tracker_rawid,
		lhs.shifter_spot, lhs.zone1, lhs.raw1, lhs.zone2, lhs.raw2, lhs.BM_bunch, lhs.BM_nplane, lhs.PM_track_type, lhs.ECC_track_type,
		lhs.chi2_shifter[0], lhs.chi2_shifter[1], lhs.chi2_shifter[2], lhs.chi2_shifter[3]) <
		std::tie(rhs.tk, rhs.eventid, rhs.PL, rhs.ISS_rawid, rhs.OSS_rawid, rhs.fixed_rawid, rhs.tracker_rawid,
			rhs.shifter_spot, rhs.zone1, rhs.raw1, rhs.zone2, rhs.raw2, rhs.BM_bunch, rhs.BM_nplane, rhs.PM_track_type, rhs.ECC_track_type,
			rhs.chi2_shifter[0], rhs.chi2_shifter[1], rhs.chi2_shifter[2], rhs.chi2_shifter[3]);
}

struct Lst {
	int gid, rawid3, rawid4;
};
bool operator<(const Lst& lhs, const Lst& rhs) {
	return std::tie(lhs.gid, lhs.rawid3, lhs.rawid4) < std::tie(rhs.gid, rhs.rawid3, rhs.rawid4);
}

void output_duplicate_muons(std::multimap<TrackerInfo, SharingFile >& map, std::string output, std::string output0);

void find_same_muon_prediction(std::string input, std::multimap<TrackerInfo, SharingFile >& map);
void set_stop_muon_rawid_list(std::string input, std::multimap<int, Lst>& list);
void matching_sf_stop_mu(std::multimap<TrackerInfo, SharingFile >& map, std::multimap<int, Lst>& list, std::string output);



void add_momch(std::vector<Momentum_recon::Event_information>& momch, std::vector<Momentum_recon::Event_information>& momch_add, int add_event, int add_chain);
void read_list(std::string in, std::set<int>& list);
void erase_momch(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch, std::set<int>& set);

int main(int argc, char** argv) {
	if (argc == 3) {
		fprintf(stderr, "usage : input_sharingfile multi-sf.txt unique-sf.txt \n");
		exit(1);
	}
	std::string input = argv[1];//input momch
	std::string output = argv[2];//multi
	std::string output2 = argv[3];//uniqe
	std::multimap<TrackerInfo, SharingFile > map;

	find_same_muon_prediction(input, map);
	output_duplicate_muons(map, output,output2);

}

void find_same_muon_prediction(std::string input, std::multimap<TrackerInfo, SharingFile >& map)
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
	SharingFile sf;
	while (std::getline(ifs, str)) {

		str_v = StringSplit_with_tab(str);

		sf.PL = std::stoi(str_v[0]);
		sf.ISS_rawid = std::stoi(str_v[1]);
		sf.OSS_rawid = std::stoi(str_v[2]);
		sf.fixed_rawid = std::stoi(str_v[3]);
		sf.tracker_rawid = std::stoi(str_v[4]);
		sf.shifter_spot = std::stoi(str_v[5]);
		sf.zone1 = std::stoi(str_v[6]);
		sf.raw1 = std::stoi(str_v[7]);
		sf.zone2 = std::stoi(str_v[8]);
		sf.raw2 = std::stoi(str_v[9]);
		sf.tk.unix_time = std::stoi(str_v[10]);
		sf.tk.t_track_id = std::stoi(str_v[11]);
		sf.BM_bunch = std::stoi(str_v[12]);
		sf.tk.entry_in_daily_file = std::stoi(str_v[13]);
		sf.BM_nplane = std::stoi(str_v[14]);
		sf.charge = std::stoi(str_v[15]);
		sf.tk.BM_mom = std::stoi(str_v[16]);
		sf.chi2_shifter[0] = std::stoi(str_v[17]);
		sf.chi2_shifter[1] = std::stoi(str_v[18]);
		sf.chi2_shifter[2] = std::stoi(str_v[19]);
		sf.chi2_shifter[3] = std::stoi(str_v[20]);
		sf.eventid = std::stoi(str_v[21]);
		sf.PM_track_type = std::stoi(str_v[22]);
		sf.ECC_track_type = std::stoi(str_v[23]);

		map.insert(std::make_pair(sf.tk, sf));
		count++;
	}
	std::cout << " sf has " << count << " lines" << std::endl;

}
void output_duplicate_muons(std::multimap<TrackerInfo, SharingFile >& map, std::string output, std::string output0)
{
	std::ofstream ofs(output);
	std::ofstream ofs0(output0);

	std::set<TrackerInfo>set;
	for (auto itr = map.begin(); itr != map.end(); itr++) {
		set.insert(itr->first);
	}

	int count = 0;
	int multi = 0;
	for (auto itr = set.begin(); itr != set.end(); itr++) {
		if (map.count(*itr) == 1) {
			auto k = map.find(*itr);
			ofs0 << std::fixed << std::right
				<< std::setw(3) << std::setprecision(0) << k->second.PL << " "
				<< std::setw(12) << std::setprecision(0) << k->second.ISS_rawid << " "
				<< std::setw(12) << std::setprecision(0) << k->second.OSS_rawid << " "
				<< std::setw(12) << std::setprecision(0) << k->second.fixed_rawid << " "
				<< std::setw(12) << std::setprecision(0) << k->second.tracker_rawid << " "
				<< std::setw(5) << std::setprecision(0) << k->second.shifter_spot << " "
				<< std::setw(10) << std::setprecision(0) << k->second.zone1 << " "
				<< std::setw(12) << std::setprecision(0) << k->second.raw1 << " "
				<< std::setw(10) << std::setprecision(0) << k->second.zone2 << " "
				<< std::setw(12) << std::setprecision(0) << k->second.raw2 << " "
				<< std::setw(20) << std::setprecision(0) << k->second.tk.unix_time << " "
				<< std::setw(4) << std::setprecision(0) << k->second.tk.t_track_id << " "
				<< std::setw(4) << std::setprecision(0) << k->second.BM_bunch << " "
				<< std::setw(10) << std::setprecision(0) << k->second.tk.entry_in_daily_file << " "
				<< std::setw(5) << std::setprecision(0) << k->second.BM_nplane << " "
				<< std::setw(3) << std::setprecision(0) << k->second.charge << " "
				<< std::setw(7) << std::setprecision(1) << k->second.tk.BM_mom << " "
				<< std::setw(5) << std::setprecision(1) << k->second.chi2_shifter[0] << " "
				<< std::setw(5) << std::setprecision(1) << k->second.chi2_shifter[1] << " "
				<< std::setw(5) << std::setprecision(1) << k->second.chi2_shifter[2] << " "
				<< std::setw(5) << std::setprecision(1) << k->second.chi2_shifter[3] << " "
				<< std::setw(10) << std::setprecision(0) << k->second.eventid << " "
				<< std::setw(3) << std::setprecision(0) << k->second.PM_track_type << " "
				<< std::setw(3) << std::setprecision(0) << k->second.ECC_track_type
				<< std::endl;
			count++;
		}
		else {
			auto p = map.equal_range(*itr);
			for (auto j = p.first; j != p.second; j++) {
				ofs << std::fixed << std::right
					<< std::setw(3) << std::setprecision(0) << j->second.PL << " "
					<< std::setw(12) << std::setprecision(0) << j->second.ISS_rawid << " "
					<< std::setw(12) << std::setprecision(0) << j->second.OSS_rawid << " "
					<< std::setw(12) << std::setprecision(0) << j->second.fixed_rawid << " "
					<< std::setw(12) << std::setprecision(0) << j->second.tracker_rawid << " "
					<< std::setw(5) << std::setprecision(0) << j->second.shifter_spot << " "
					<< std::setw(10) << std::setprecision(0) << j->second.zone1 << " "
					<< std::setw(12) << std::setprecision(0) << j->second.raw1 << " "
					<< std::setw(10) << std::setprecision(0) << j->second.zone2 << " "
					<< std::setw(12) << std::setprecision(0) << j->second.raw2 << " "
					<< std::setw(20) << std::setprecision(0) << j->second.tk.unix_time << " "
					<< std::setw(4) << std::setprecision(0) << j->second.tk.t_track_id << " "
					<< std::setw(4) << std::setprecision(0) << j->second.BM_bunch << " "
					<< std::setw(10) << std::setprecision(0) << j->second.tk.entry_in_daily_file << " "
					<< std::setw(5) << std::setprecision(0) << j->second.BM_nplane << " "
					<< std::setw(3) << std::setprecision(0) << j->second.charge << " "
					<< std::setw(7) << std::setprecision(1) << j->second.tk.BM_mom << " "
					<< std::setw(5) << std::setprecision(1) << j->second.chi2_shifter[0] << " "
					<< std::setw(5) << std::setprecision(1) << j->second.chi2_shifter[1] << " "
					<< std::setw(5) << std::setprecision(4) << j->second.chi2_shifter[2] << " "
					<< std::setw(5) << std::setprecision(1) << j->second.chi2_shifter[3] << " "
					<< std::setw(10) << std::setprecision(0) << j->second.eventid << " "
					<< std::setw(3) << std::setprecision(0) << j->second.PM_track_type << " "
					<< std::setw(3) << std::setprecision(0) << j->second.ECC_track_type
					<< std::endl;
			}
		}
	}
	std::cout << " unique / all = " << set.size() << " / " << map.size() << std::endl;
	std::cout << " only  unique = " << count << std::endl;

	//	std::cout << "total prediction                   : " << map.size()
	//		    << "\nuniqe muon prediction number       : " << set.size() <<
	//		       "\nnum of duplicated muon predictions : " << count <<" ("<<multi<<")"<<
	//		std::endl;
}
void set_stop_muon_rawid_list(std::string input, std::multimap<int, Lst>& list)
{
	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "file open error" << std::endl;
		exit(1);
	}

	Lst l;
	//while (ifs >> l.gid >> l.rawid3 >> l.rawid4) {
	//	list.insert(std::make_pair(l.gid, l));
	//}

	std::string str;					//1strein into
	std::vector<std::string> str_v;		//input 1 ward
	std::string buffer;

	while (std::getline(ifs, str)) {

		str_v = StringSplit_with_tab(str);

		l.gid = std::stoi(str_v[0]);
		l.rawid3 = std::stoi(str_v[1]);
		l.rawid4 = std::stoi(str_v[2]);
		list.insert(std::make_pair(l.gid, l));
		//std::cout << l.gid << std::endl;
	}
}
void matching_sf_stop_mu(std::multimap<TrackerInfo, SharingFile >& map, std::multimap<int, Lst>& list, std::string output)
{
	std::set<TrackerInfo>set;
	for (auto itr = map.begin(); itr != map.end(); itr++) {
		set.insert(itr->first);
	}

	std::ofstream ofs(output);

	int cnt = 0;
	for (auto itr = set.begin(); itr != set.end(); itr++) {//sharing file
		if (map.count(*itr) == 1)continue;
		auto p = map.equal_range(*itr);
		std::cout << "--- UnixTime : " << p.first->first.unix_time << " ---" << std::endl;

		for (auto j = p.first; j != p.second; j++) {//sfの重複
			//auto  r = list.find(j->second.eventid);
			//if (r != list.end())continue;
			if (list.count(j->second.eventid) == 1)continue;//dupliccatino in mfile
			for (auto r = list.begin(); r != list.end(); r++) {//mfile
				if (r->first == j->second.eventid) {
					std::cout << r->second.gid << "\t" << r->second.rawid3 << " " << r->second.rawid4 << std::endl;
					if (r->second.rawid3 > -1 && r->second.rawid4 > -1) {
						ofs << r->second.gid << "\t" << r->second.rawid3 << " " << r->second.rawid4 << std::endl;
					}
					cnt++;
				}
			}


			//std::cout << r->second.gid << "\t" << r->second.rawid3 << " " << r->second.rawid4 << std::endl;
		}
	}
	std::cout << cnt << std::endl;
}

void add_momch(std::vector<Momentum_recon::Event_information>& momch, std::vector<Momentum_recon::Event_information>& momch_add, int add_event, int add_chain) {
	Momentum_recon::Mom_chain chain;
	for (auto& ev : momch_add) {
		if (ev.groupid != add_event)continue;
		for (auto& c : ev.chains) {
			if (c.chainid != add_chain)continue;
			chain = c;
		}
	}

	for (auto& ev : momch) {
		if (ev.groupid != add_event)continue;
		ev.chains.push_back(chain);
	}

}
void read_list(std::string in, std::set<int>& list) {
	std::ifstream ifs(in);
	if (!ifs) {
		std::cerr << "reading file error" << std::endl;
	}

	int i;
	while (ifs >> i) {
		list.insert(i);
	}
	std::cout << "size:" << list.size() << std::endl;
}
void erase_momch(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch, std::set<int>& set) {
	Momentum_recon::Mom_chain chain;

	for (auto& ev : momch0) {
		auto itr = set.find(ev.groupid);
		if (itr != set.end()) {
			momch.push_back(ev);
		}
	}

}