#define _USE_MATH_DEFINES
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <set>
#include <map>
#include <math.h>
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

class pid_param {
public:
	int i_ang, i_mom;
	double e_param[3], e_error[3], pi_param[3], pi_error[3], p_param[3], p_error[3], p_thr;
	void Cacl_p_thr(double expected);
};
bool Calc_average(Momentum_recon::Mom_chain& c, int& count_vph, float& average_vph, int& count_pixel, float& average_pixel);
std::map<std::pair<int, int>, pid_param> read_param(std::string filename);
int ang_to_index(double angle);
int mom_to_index(double mom);
std::vector<Momentum_recon::Event_information> Divide_proton(std::vector<Momentum_recon::Event_information>& momch, std::map<std::pair<int, int>, pid_param>& param);
double sum_gaus(int min, int max, double param[3]);


int main(int argc, char** argv) {
	if (argc != 4) {
		fprintf(stderr, "filename");
		exit(1);
	}
	std::string file_in_momch = argv[1];
	std::string file_in_param = argv[2];
	//std::string file_out_momch_mip = argv[3];
	std::string file_out_momch_proton = argv[3];

	std::map<std::pair<int, int>, pid_param>param = read_param(file_in_param);
	
	std::cout << "\t* Now Reading Momch *" << std::endl;
	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch);
	std::vector<Momentum_recon::Event_information> proton = Divide_proton(momch, param);

	Momentum_recon::Write_Event_information_extension(file_out_momch_proton, proton);
	exit(0);
}
void pid_param::Cacl_p_thr(double sn) {
	if (e_param[1] > 0 && pi_param[1] > 0 && p_param[1] > 0) {
		for (int thr = p_param[1] - p_param[2] * 1.5; thr < 500; thr++) {
			double e_num = sum_gaus(thr, 500, e_param);
			double pi_num = sum_gaus(thr, 500, pi_param);
			double p_num = sum_gaus(thr, 500, p_param);
			if (p_num / (e_num + pi_num + p_num) > sn) {
				p_thr = thr;
				return;
			}
		}
		p_thr = 500;
	}
	else if (pi_param[1] > 0 && p_param[1] > 0) {
		for (int thr = p_param[1] - p_param[2] * 1.5; thr < 500; thr++) {
			double pi_num = sum_gaus(thr, 500, pi_param);
			double p_num = sum_gaus(thr, 500, p_param);
			if (p_num / (pi_num + p_num) > sn) {
				p_thr = thr;
				return;
			}
		}
		p_thr = 500;
	}
	else {
		p_thr = 500;
	}
}
double sum_gaus(int min, int max, double param[3]) {
	double ret = 0;
	for (int i = min; i < max; i++) {
		ret += param[0] / (sqrt(2 * M_PI) * param[2]) * exp(-1 * pow(i - param[1], 2) / (2 * param[2] * param[2]));
	}
	return ret;
}

std::map<std::pair<int, int>, pid_param> read_param(std::string filename) {
	std::ifstream ifs(filename);
	pid_param p;
	std::map<std::pair<int, int>, pid_param> ret;

	std::cout << "\t* Now Reading parameters *" << std::endl;
	while (ifs >> p.i_ang >> p.i_mom
		>> p.e_param[0] >> p.e_error[0] >> p.e_param[1] >> p.e_error[1] >> p.e_param[2] >> p.e_error[2]
		>> p.pi_param[0] >> p.pi_error[0] >> p.pi_param[1] >> p.pi_error[1] >> p.pi_param[2] >> p.pi_error[2]
		>> p.p_param[0] >> p.p_error[0] >> p.p_param[1] >> p.p_error[1] >> p.p_param[2] >> p.p_error[2]) {
		p.Cacl_p_thr(0.9999);
		ret.insert(std::make_pair(std::make_pair(p.i_ang, p.i_mom), p));
		//printf("%d %d %lf\n", p.i_ang, p.i_mom, p.p_thr);
	}
	std::cout << "\t  Fin!  " << std::endl;
	return ret;
}
int mom_to_index(double mom) {
	return int(mom / 100);
}
int ang_to_index(double angle) {
	int i_ang = 0;
	if (angle < 0.7) {
		i_ang = angle / 0.1;
	}
	else if (angle < 1.5) {
		i_ang = (angle - 0.7) / 0.2 + 7;
	}
	else if (angle < 3.1) {
		i_ang = (angle - 1.5) / 0.4 + 11;
	}
	else {
		i_ang = (angle - 3.1) / 0.6 + 14;
	}
	return i_ang;

}
bool Calc_average(Momentum_recon::Mom_chain& c, int& count_vph, float& average_vph, int& count_pixel, float& average_pixel) {
	//std::cout << "\t* Now Calculating averages *" << std::endl;
	count_vph = 0;
	count_pixel = 0;
	average_vph = 0;
	average_pixel = 0;
	for (auto itr = c.base.begin(); itr != c.base.end(); itr++) {
		for (int i = 0; i < 2; i++) {
			count_vph += 1;
			average_vph += itr->m[i].ph % 10000;
			if (itr->m[i].hitnum > 0) {
				count_pixel += 1;
				average_pixel += itr->m[i].hitnum;
			}
		}
	}
	if (count_pixel <= 0)return false;
	if (count_vph <= 0)return false;

	average_pixel /= count_pixel;
	average_vph /= count_vph;
	return true;

}
std::vector<Momentum_recon::Event_information> Divide_proton(std::vector<Momentum_recon::Event_information>& momch, std::map<std::pair<int, int>, pid_param>& param) {
	std::vector<Momentum_recon::Event_information>  ret;
	double ax, ay, angle;
	int count_vph, count_pixel;
	float average_vph, average_pixel;
	for (auto& ev : momch) {
		Momentum_recon::Event_information ev_inf;
		ev_inf = ev;
		ev_inf.chains.clear();
		for (auto& c : ev.chains) {
			int count = 0;
			for (auto itr = c.base.begin(); itr != c.base.end(); itr++) {
				ax += itr->ax;
				ay += itr->ay;
				count++;
			}
			ax /= count;
			ay /= count;
			angle = sqrt(ax * ax + ay * ay);
			if (Calc_average(c, count_vph, average_vph, count_pixel, average_pixel)) {
				int i_ang = ang_to_index(angle);
				int i_mom = mom_to_index(c.Get_muon_mcs_pb());
				auto p = param.find(std::make_pair(i_ang, i_mom));
				if (p == param.end())continue;
				if (p->second.p_thr <= average_vph) {
					ev_inf.chains.push_back(c);
				}
			}
		}
		if (ev_inf.chains.size() > 0) {
			ret.push_back(ev_inf);
		}
	}

	printf("proton num=%d\n", ret.size());
	return ret;

}
