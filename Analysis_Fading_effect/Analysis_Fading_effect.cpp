#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <set>
#include <map>
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>
class output_format {
public:
	int groupid, chainid, nseg, npl, count, pl, face;
	double  ecc_mcs_mom, ax, ay, angle, average0, average1, error0, error1;
	double unixtime;
};

int judege_pl_id(int id);
std::vector<output_format> Calc_average_momch(std::vector<Momentum_recon::Event_information>& momch);
bool Calc_average_vph(Momentum_recon::Mom_chain& c, Momentum_recon::Mom_basetrack& base, int face, output_format& out);
void output(std::string file_out, std::vector<output_format>& ret);
void output_bin(std::string filename, std::vector<output_format>& ret);

int main(int argc, char** argv) {
	if (argc != 3) {
		fprintf(stderr, "usage:filename in.momch out.bin\n");
		exit(1);
	}
	std::string file_in_momch = argv[1];
	std::string file_out = argv[2];
	//std::string file_out_bin = argv[2];

	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch);

	//std::vector<Momentum_recon::Mom_chain> momch = Momentum_recon::Read_mom_chain_extension(file_in_momch);
	std::vector<output_format>ave = Calc_average_momch(momch);
	//output_bin(file_out_bin, ave);
	std::cout << ave.size() << std::endl;
	output(file_out, ave);

	exit(0);
}
std::vector<output_format> Calc_average_momch(std::vector<Momentum_recon::Event_information>& momch) {
	std::vector<output_format> ret;
	int all = momch.size(), now = 0;



	for (auto& ev : momch) {

		//std::cout << " EventID = " << ev.groupid << std::endl;
		for (auto& c : ev.chains) {
			if (now % 10000 == 0) {
				fprintf(stderr, "\r now average calc %d/%d(%4.1lf%%)", now, all, now * 100. / all);
			}
			now++;
			if (ev.unix_time <0)continue;

			output_format out;
			out.groupid = ev.groupid;
			out.chainid = c.chainid;
			out.ecc_mcs_mom = c.ecc_mcs_mom[0];
			out.nseg = c.base.size();
			out.npl = c.base.rbegin()->pl - c.base.begin()->pl + 1;
			out.ax = 0;
			out.ay = 0;
			out.unixtime = ev.unix_time;

			int count = 0;
			for (auto itr = c.base.begin(); itr != c.base.end(); itr++) {
				out.ax += itr->ax;
				out.ay += itr->ay;
				count++;
			}
			out.ax /= count;
			out.ay /= count;
			out.angle = sqrt(out.ax * out.ax + out.ay * out.ay);
			for (auto& b : c.base) {
				out.pl = b.pl;
				if (Calc_average_vph(c, b, 0, out)) {
					out.face = 0;
					ret.push_back(out);
				}
				if (Calc_average_vph(c, b, 1, out)) {
					out.face = 1;
					ret.push_back(out);
				}
			}
		}
	}

	fprintf(stderr, "\r now average calc %d/%d(%4.1lf%%)\n", now, all, now * 100. / all);

	return ret;

}
bool Calc_average_vph(Momentum_recon::Mom_chain& c, Momentum_recon::Mom_basetrack& base, int face, output_format& out) {
	bool detect_flg = false;
	double sum = 0, sum2 = 0, count = 0;
	int vph;
	for (auto& b : c.base) {
		for (int i = 0; i < 2; i++) {
			vph = b.m[i].ph % 10000;

			if (b.pl == base.pl && i == face) {

				out.average0 = vph;
				out.error0 = 0;
				detect_flg = true;
			}
			else {
				sum += vph;
				sum2 += pow(vph, 2);
				count += 1;
			}
		}
	}
	if (!detect_flg)return false;
	if (count < 3)return false;

	out.count = count;
	out.average1 = sum / count;
	double sig;
	sig = sum2 / count - out.average1 * out.average1;
	if (sig <= 0)return false;
	out.error1 = sqrt(sig) * sqrt(count) / (count - 1);
	out.error0 = sqrt(sig) * sqrt(count) / sqrt(count - 1);

	return true;
}
int judege_pl_id(int id) {
	if ((24 <= id && id <= 35) || id == 52)return 0;
	else return 1;
}

void output(std::string file_out, std::vector<output_format>& ret) {
	std::ofstream ofs(file_out);
	int count = 0, all = ret.size();
	for (auto itr = ret.begin(); itr != ret.end(); itr++) {
		if (count % 10000 == 0) {
			fprintf(stderr, "\r write file %10d/%10d(%4.1lf%%)", count, all, count * 100. / all);
		}
		count++;
		ofs << std::right << std::fixed
			<< std::setw(10) << std::setprecision(0) << itr->groupid << " "
			<< std::setw(10) << std::setprecision(0) << itr->chainid << " "
			<< std::setw(3) << std::setprecision(0) << itr->nseg << " "
			<< std::setw(3) << std::setprecision(0) << itr->npl << " "
			<< std::setw(7) << std::setprecision(4) << itr->ax << " "
			<< std::setw(7) << std::setprecision(4) << itr->ay << " "
			<< std::setw(7) << std::setprecision(4) << itr->angle << " "
			<< std::setw(10) << std::setprecision(3) << itr->ecc_mcs_mom << " "
			<< std::setw(5) << std::setprecision(0) << itr->count << " "
			<< std::setw(8) << std::setprecision(3) << itr->average0 << " "
			<< std::setw(8) << std::setprecision(3) << itr->error0 << " "
			<< std::setw(8) << std::setprecision(3) << itr->average1 << " "
			<< std::setw(8) << std::setprecision(3) << itr->error1 <<" "
			<< std::setw(12) << std::setprecision(3) << itr->unixtime
			<< std::endl;
	}
	fprintf(stderr, "\r write file %10d/%10d(%4.1lf%%)\n", count, all, count * 100. / all);

}
void output_bin(std::string filename, std::vector<output_format>& ret) {
	std::ofstream ofs(filename, std::ios::binary);
	if (!ofs) {
		//file open s
		fprintf(stderr, "File[%s] is not exist!!\n", filename.c_str());
		exit(1);
	}
	if (ret.size() == 0) {
		fprintf(stderr, "target linklet ... null\n");
		fprintf(stderr, "File[%s] has no text\n", filename.c_str());
	}
	int64_t count = 0;
	int64_t max = ret.size();
	for (int i = 0; i < ret.size(); i++) {
		if (count % 10000 == 0) {
			std::cerr << std::right << std::fixed << "\r now writing ..." << std::setw(4) << std::setprecision(1) << count * 100. / max << "%";
		}
		count++;
		ofs.write((char*)&ret[i], sizeof(output_format));
	}
	std::cerr << std::right << std::fixed << "\r now writing ..." << std::setw(4) << std::setprecision(1) << count * 100. / max << "%" << std::endl;
	ofs.close();

}
