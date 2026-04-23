#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

std::vector<Momentum_recon::Event_information> merge_momch(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch1);

int main(int argc, char** argv) {
	if (argc < 4) {
		fprintf(stderr, "usage:file-in-momchain file-in-momchain file-out-momchain [overwrite mode]\nIf you set the overwrite mode flag to 1, \nany existing event will be overwritten with the values from the second file.\nIf not specified, the values from the first file will be output.\n");
		exit(1);
	}
	std::string file_in_momch0 = argv[1];
	std::string file_in_momch1 = argv[2];
	std::string file_out_momch = argv[3];
	int mode = -1;
	if (argc > 4) {
		mode = std::stoi(argv[4]);
	}

	std::vector<Momentum_recon::Event_information> momch0, momch1;
	momch0 = Momentum_recon::Read_Event_information_extension(file_in_momch0);
	momch1 = Momentum_recon::Read_Event_information_extension(file_in_momch1);

	std::vector<Momentum_recon::Event_information> m_momch = merge_momch(momch0, momch1);

	Momentum_recon::Write_Event_information_extension(file_out_momch, m_momch);

}

std::vector<Momentum_recon::Event_information> merge_momch(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch1) {
	std::map<int, Momentum_recon::Event_information> m_momch;


	for (auto& ev : momch0) {
		m_momch.insert(std::make_pair(ev.groupid, ev));
	}

	for (auto& ev : momch1) {
		auto res = m_momch.find(ev.groupid);
		if (res == m_momch.end()) {
			//fprintf(stderr, "same chain was input groupid=%d chainid=%d\n", ev.groupid, c.chainid);
			std::cout << "not existed event:" << ev.groupid << std::endl;
			m_momch.insert(std::make_pair(ev.groupid, ev));
		}

	}

	std::cout << "momch1 size: " << momch0.size() << std::endl;
	std::cout << "momch2 size: " << momch1.size() << std::endl;
	std::cout << "merged size: " << m_momch.size() << std::endl;


	std::vector<Momentum_recon::Event_information > ret;

	for (auto itr = m_momch.begin(); itr != m_momch.end(); itr++) {
		ret.push_back(itr->second);
	}
	return ret;

}

std::vector<Momentum_recon::Event_information> merge_momch_overwrite(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch1) {
	std::map<int, Momentum_recon::Event_information> m_momch;


	for (auto& ev : momch0) {
		m_momch.insert(std::make_pair(ev.groupid, ev));
	}

	int cnt = 0;
	for (auto& ev : momch1) {
		auto res = m_momch.find(ev.groupid);
		if (res == m_momch.end()) {
			std::cout << "not existed event:" << ev.groupid << std::endl;
			m_momch.insert(std::make_pair(ev.groupid, ev));
		}
		else {// duplicated
			std::cout << "overwrite :" << ev.groupid << std::endl;
			res->second = ev;
			cnt++;
		}

	}

	std::cout << "momch1 size: " << momch0.size() << std::endl;
	std::cout << "momch2 size: " << momch1.size() << std::endl;
	std::cout << "merged size: " << m_momch.size() << std::endl;
	std::cout << "overwrite size: " << cnt << std::endl;


	std::vector<Momentum_recon::Event_information > ret;

	for (auto itr = m_momch.begin(); itr != m_momch.end(); itr++) {
		ret.push_back(itr->second);
	}
	return ret;

}

std::vector<Momentum_recon::Event_information> merge_momch_old(std::vector<Momentum_recon::Event_information>& momch0, std::vector<Momentum_recon::Event_information>& momch1) {
	std::map<std::pair<int, int>, Momentum_recon::Event_information> m_momch;

	for (auto& ev : momch0) {
		for (auto& c : ev.chains) {
			m_momch.insert(std::make_pair(std::make_pair(ev.groupid, c.chainid), ev));
		}
	}

	for (auto& ev : momch1) {
		for (auto& c : ev.chains) {
			auto res = m_momch.insert(std::make_pair(std::make_pair(ev.groupid, c.chainid), ev));
			if (!res.second) {
				fprintf(stderr, "same chain input groupid=%d chainid=%d\n", ev.groupid, c.chainid);
			}

		}
	}

	std::vector<Momentum_recon::Event_information > ret;

	for (auto itr = m_momch.begin(); itr != m_momch.end(); itr++) {
		ret.push_back(itr->second);
	}
	return ret;

}