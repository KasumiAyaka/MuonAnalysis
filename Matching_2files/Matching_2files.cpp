// 2024/08/26
// kasumi

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <utility>
#include <iomanip>

struct List1 {
	int eventid, pl, unixtime;
};

struct List2 {
	int bsd_spillnum;
	int unixtime;
};
bool operator<(const List2& lhs, const List2& rhs) {
	return std::tie(lhs.bsd_spillnum, lhs.unixtime) < std::tie(rhs.bsd_spillnum, rhs.unixtime);
}

std::vector<List1> ReadList1(std::string input);
std::multimap<int, List2> ReadList2(std::string input);

void Matching(std::vector<List1>& l1, std::multimap<int, List2>& l2, std::string output);
void Matching2(std::vector<List1>& l1, std::multimap<int, List2>& l2, std::string output);


int main(int argc, char** argv) {
	if (argc != 5) {
		fprintf(stderr, "usage:LocatedEvents EventViewerLists output1 output2\n");
		exit(1);
	}
	std::string input1 = argv[1];
	std::string input2 = argv[2];
	std::string output1 = argv[3];
	std::string output2 = argv[4];

	std::vector<List1>ls1;
	std::multimap<int, List2>ls2;

	ls1 = ReadList1(input1);//NETSCANBACK
	ls2 = ReadList2(input2);//EVENT DISPLAY

	Matching(ls1, ls2, output1);
	Matching2(ls1, ls2, output2);
	std::cout << "Finished." << std::endl;

}

std::vector<List1> ReadList1(std::string input) {
	std::vector<List1> ret;
	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "File open error : " << input << std::endl;
		exit(1);
	}

	List1 l = { 0 };
	std::cout << "now reading..." << std::endl;
	while (ifs>>l.eventid >> l.pl >> l.unixtime) {
		ret.push_back(l);
	}
	std::cout << ret.size() << std::endl;
	return ret;
}

std::multimap<int, List2> ReadList2(std::string input) {
	std::multimap<int, List2> ret;
	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "File open error : " << input << std::endl;
		exit(1);
	}

	std::cout << "now reading..." << std::endl;
	List2 l = { 0 };
	while (ifs >> l.bsd_spillnum >> l.unixtime) {
		ret.insert(std::make_pair(l.unixtime, l));
	}
	std::cout << ret.size() << std::endl;

	return ret;
}

void Matching(std::vector<List1>& l1, std::multimap<int, List2>& l2, std::string output) {

	std::ofstream ofs(output);
	std::cout << "now Matching..." << std::endl;

	for (auto itr = l1.begin(); itr != l1.end(); itr++) {//net scan back
		auto p = l2.find(itr->unixtime);//event display
		if (p != l2.end()) {
			ofs << std::right << std::fixed
				<< std::setw(5) << std::setprecision(0) << itr->eventid << " "
				<< std::setw(3) << std::setprecision(0) << itr->pl << " "
				<< std::setw(12) << std::setprecision(0) << itr->unixtime << " "
				<< std::setw(10) << std::setprecision(0) << p->second.bsd_spillnum << " "
				<< std::endl;
		}
		else {
			auto p = l2.find(itr->unixtime + 1);
			if (p != l2.end()) {
				ofs << std::right << std::fixed
					<< std::setw(5) << std::setprecision(0) << itr->eventid << " "
					<< std::setw(3) << std::setprecision(0) << itr->pl << " "
					<< std::setw(12) << std::setprecision(0) << itr->unixtime << " "
					<< std::setw(10) << std::setprecision(0) << p->second.bsd_spillnum << " "
					<< std::endl;
			}
			else {
				ofs << std::right << std::fixed
					<< std::setw(5) << std::setprecision(0) << itr->eventid << " "
					<< std::setw(3) << std::setprecision(0) << itr->pl << " "
					<< std::setw(12) << std::setprecision(0) << itr->unixtime << " "
					<< std::setw(10) << std::setprecision(0) << -1 << " "
					<< std::endl;

			}
		}
	}

}

void Matching2(std::vector<List1>& l1, std::multimap<int, List2>& l2, std::string output) {

	std::ofstream ofs(output);
	std::multimap<int, List1> lst;
	std::cout << "now Matching..." << std::endl;

	// located event list
	for (auto itr = l1.begin(); itr != l1.end(); itr++) {
		lst.insert(std::make_pair(itr->unixtime, *itr));//netscanback
	}
	
	// search from event_viewer-list
	for (auto itr = l2.begin(); itr != l2.end();itr++) {
		auto p = lst.find(itr->first);
		if (p != lst.end()) {
			ofs << std::right << std::fixed
				<< std::setw(12) << std::setprecision(0) << itr->first << " "
				<< std::setw(10) << std::setprecision(0) << itr->second.bsd_spillnum << " "
				<< std::setw(5) << std::setprecision(0) << p->second.eventid << " "
				<< std::setw(3) << std::setprecision(0) << p->second.pl << " "
				<< std::endl;
		}
		else {
			auto p = lst.find(itr->first - 1);
			if (p != lst.end()) {
				ofs << std::right << std::fixed
					<< std::setw(12) << std::setprecision(0) << itr->first << " "
					<< std::setw(10) << std::setprecision(0) << itr->second.bsd_spillnum << " "
					<< std::setw(5) << std::setprecision(0) << p->second.eventid << " "
					<< std::setw(3) << std::setprecision(0) << p->second.pl << " "
					<< std::endl;
			}
			else {
				ofs << std::right << std::fixed
					<< std::setw(12) << std::setprecision(0) << itr->first << " "
					<< std::setw(10) << std::setprecision(0) << itr->second.bsd_spillnum << " "
					<< std::setw(5) << std::setprecision(0) << -1 << " "
					<< std::setw(3) << std::setprecision(0) << -1 << " "
					<< std::endl;
			}
		}
	}



}
