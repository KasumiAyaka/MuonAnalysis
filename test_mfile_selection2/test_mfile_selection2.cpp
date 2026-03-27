#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

#include <filesystem>

void output_chain_inf(std::string filename, std::vector<mfile0::M_Chain>& c);
inline bool judge_black_chain(mfile0::M_Chain& c);

int main(int argc, char** argv) {
	if (argc != 3) {
		fprintf(stderr, "usage:file-in-mfile output_mfile\n");
		exit(1);
	}
	std::string file_in_mfile = argv[1];
	std::string file_out_mfile = argv[2];

	mfile0::Mfile m;
	mfile1::read_mfile_extension(file_in_mfile, m);

	mfile0::Mfile m_all;
	m_all.header = m.header;
	int count = 0;
	for (auto itr = m.chains.begin(); itr != m.chains.end(); itr++) {
		if (itr->nseg <= 6)continue;
		m_all.chains.push_back(*itr);
	}

	mfile1::write_mfile_extension(file_out_mfile, m_all);

	//output_chain_inf(file_out_inf, m.chains);

}

inline bool judge_black_chain(mfile0::M_Chain& c) {

	double ax = mfile0::chain_ax(c);
	double ay = mfile0::chain_ay(c);
	double vph = mfile0::chain_vph(c);
	double angle = sqrt(ax * ax + ay * ay);
	if (angle < 0.4) {
		return vph > -200 * angle + 200;
	}
	else if (angle < 1.0) {
		return vph > (-100 * angle + 400) / 3;
	}
	return vph > 100;


}

void output_chain_inf(std::string filename, std::vector<mfile0::M_Chain>& c) {
	int count = 0, all = c.size();

	std::ofstream ofs(filename);
	double ax, ay, vph;
	for (auto itr = c.begin(); itr != c.end(); itr++) {
		if (count % 10000 == 0) {
			fprintf(stderr, "\r write chain inf %10d/%10d(%4.1lf%%)", count, all, count * 100. / all);
		}
		count++;


		ax = mfile0::chain_ax(*itr);
		ay = mfile0::chain_ay(*itr);
		vph = mfile0::chain_vph(*itr);
		ofs << std::right << std::fixed
			<< std::setw(12) << std::setprecision(0) << itr->chain_id << " "
			<< std::setw(3) << std::setprecision(0) << itr->nseg << " "
			<< std::setw(3) << std::setprecision(0) << itr->pos0 / 10 << " "
			<< std::setw(3) << std::setprecision(0) << itr->pos1 / 10 << " "
			<< std::setw(7) << std::setprecision(4) << ax << " "
			<< std::setw(7) << std::setprecision(4) << ay << " "
			<< std::setw(7) << std::setprecision(2) << vph << std::endl;

	}
	fprintf(stderr, "\r write chain inf %10d/%10d(%4.1lf%%)\n", count, all, count * 100. / all);

}