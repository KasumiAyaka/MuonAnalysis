// 2025/03/05
// kasumi
// timestamp bunch‚Ì‘g‚É‘Î‚µ‚Äunique‚È•¨‚Ì”‚ğ”‚¦‚é

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <set>
#include <map>
#include <chrono>


struct Key {
	int bunch;
	int utime;
};
bool operator<(const Key& lhs, const Key& rhs) {
	return std::tie(lhs.utime, lhs.bunch) < std::tie(rhs.utime, rhs.bunch);
}

void Count(std::string input);

int main(int argc, char** argv) {

	if (argc != 2) {
		fprintf(stderr, "usage1\n file_in_mfile file_in_base\n");		exit(1);
	}


	std::string file_in_muon_mfile = argv[1];
	Count(file_in_muon_mfile);

}

void Count(std::string input) {

	std::ifstream ifs(input);
	if (!ifs) {
		std::cerr << "File open error." << std::endl;
		exit(0);
	}

	int cnt = 0;
	Key k;
	std::set<Key> set;
	while (ifs >> k.utime >> k.bunch) {
		set.insert(k);
		cnt++;
	}

	std::cout << "unique event = " << set.size() << " / " << cnt << std::endl;
}
