// 2024/08/26
// kasumi

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <utility>
#include <iomanip>

struct Results {
	int utime, bunch, ecc, eid, pl, muflg, area;
};

struct List2 {
	int eid, a, pl;
};
bool operator<(const List2& lhs, const List2& rhs) {
	return std::tie(lhs.eid, lhs.pl) < std::tie(rhs.eid, rhs.pl);
}

std::vector<Results> ReadList1(std::string input);
std::multimap<int, List2> ReadList2(std::string input);
void Matching(std::vector<Results>& l1, std::multimap<int, List2>& l2);


int main(int argc, char** argv) {
	if (argc != 3) {
		fprintf(stderr, "usage: results list\n");
		exit(1);
	}
	std::string input1 = argv[1];//Results
	std::string input2 = argv[2];//list
	//std::string output1 = argv[3];

	std::vector<Results>ls1;
	std::multimap<int, List2>ls2;

	ls1 = ReadList1(input1);//Results
	ls2 = ReadList2(input2);//list

	Matching(ls1, ls2);
	std::cout << "Finished." << std::endl;

}

std::vector<Results> ReadList1(std::string input) {
	std::vector<Results> ret;
	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "File open error : " << input << std::endl;
		exit(1);
	}

	Results l = { 0 };
	std::cout << "now reading..." << std::endl;
	while (ifs >> l.utime>>l.bunch>>l.ecc>>l.eid>>l.pl>>l.muflg>>l.area) {
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
	while (ifs >> l.eid >> l.a>>l.pl) {
		ret.insert(std::make_pair(l.eid, l));
	}
	std::cout << ret.size() << std::endl;
	//for (auto itr = ret.begin(); itr != ret.end(); itr++) {
	//	std::cout << itr->second.eid << " " << itr->second.pl << std::endl;
	//}
	return ret;
}

void Matching(std::vector<Results>& l1, std::multimap<int, List2>& l2) {

	//std::ofstream ofs(output);
	std::cout << "now Matching..." << std::endl;

	int CentW = 0, CentF = 0, EdgeW = 0, EdgeF = 0;
	int pene = 0;
	for (auto itr = l1.begin(); itr != l1.end(); itr++) {
		auto p = l2.find(itr->eid);//event display
		if (p == l2.end())continue;
		if (itr->pl != p->second.pl) {
			std::cout << std::right << std::fixed
				<< std::setw(10) << std::setprecision(0) << itr->utime << " "
				<< std::setw(1) << std::setprecision(0) << itr->bunch << " "
				<< std::setw(1) << std::setprecision(0) << itr->ecc << " "
				<< std::setw(5) << std::setprecision(0) << itr->eid << " "
				<< std::setw(3) << std::setprecision(0) << itr->pl << " "
				<< std::setw(3) << std::setprecision(0) << itr->muflg << " "
				<< std::setw(3) << std::setprecision(0) << itr->area << " "
				<< std::setw(3) << std::setprecision(0) << p->second.pl << " "
				"different pl" << std::endl;
		}
		else {
			std::cout << std::right << std::fixed
				<< std::setw(10) << std::setprecision(0) << itr->utime << " "
				<< std::setw(1) << std::setprecision(0) << itr->bunch << " "
				<< std::setw(1) << std::setprecision(0) << itr->ecc << " "
				<< std::setw(5) << std::setprecision(0) << itr->eid << " "
				<< std::setw(3) << std::setprecision(0) << itr->pl << " "
				<< std::setw(3) << std::setprecision(0) << itr->muflg << " "
				<< std::setw(3) << std::setprecision(0) << itr->area << " "
				<< std::setw(3) << std::setprecision(0) << p->second.pl << std::endl;
		}


		if (itr->muflg < 0) {
			pene++;
		}
		else {
			if (itr->area > 0) {
				//Central
				if (itr->pl < 16 && itr->pl>130) {
					// penetrate
				}
				else {
					//stop
					if (itr->pl % 2 == 1) {
						//water
						CentW++;
					}
					else {
						CentF++;
					}
				}
			}
			else {
				//Edge
				if (itr->pl < 16 && itr->pl>130) {
					// penetrate
				}
				else {
					//stop
					if (itr->pl % 2 == 1) {
						//water
						EdgeW++;
					}
					else {
						EdgeF++;
					}
				}

			}
		}
	}
	std::cout <<" * Penetrate/Edgeout : "<<pene<<
		"\n * Central\n\tWater : " << CentW << ", Iron : " << CentF
		<< "\n* Edge\n\tWater : " << EdgeW << ", Iron : " << EdgeF << std::endl;
}

