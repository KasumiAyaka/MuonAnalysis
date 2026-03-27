// edit 2034/07/16

#pragma comment(lib, "FILE_structure.lib")
#include <FILE_structure.hpp>

#include <vector>
#include <map>
#include <fstream>
#include <ios>// std::left, std::right
#include <iomanip>
#include <math.h>
#include <algorithm>


#include "M:\data\NINJA\PID\PDG.h"
#include "M:\data\NINJA\PID\AbsorberMaterial_Emulsion_NINJA_Run6.h"


class DataPoint {
public:
	double input_ang_min, input_ang_max, input_mom_min, input_mom_max;
	double mean[3], mean_low[3], mean_hi[3], sigma[3], sigma_low[3], sigma_hi[3];
};
class Chain_inf {
public:
	int groupid, chainid, pid;
	double pb, vph, angle, proton_likelihood, pion_likelihood, likelihood_ratio;
};

class VPH_pb_func {
public:
	//keyはmomの真ん中
	std::map<double, DataPoint> datapoint;
	double vph_mean_param[4], vph_sigma_param[2], cross_point_pb, vph_sigma_thr;

	double Calc_VPH_mean_proton(double pb);
	double Calc_VPH_mean_pion(double pb);
	double Calc_VPH_sigma_proton(double pb);
	double Calc_VPH_sigma_pion(double pb);

	void SetCrossPoint_VPH();
};


std::vector<Chain_inf> Read_chainif(std::string filename);
void iang_to_angle_range(int i_ang, double& angle_min, double& angle_max);
int angle_to_iangle(double angle);

void Read_Theoretical_formula(std::string filename, double& angle_min, double& angle_max, double* vph_mean_param, double* vph_sigma_param);
void Read_Data_point(std::string filename, double& angle_min, double& angle_max, DataPoint* data, double& thr_sigma, int& DataNum);

int main(int argc, char** argv) {
	if (argc != 4&&argc!=6) {
		fprintf(stderr, "uasge:input-chain-inf output-pid　pion_thr\n");
		exit(1);
	}

	std::string file_in_chaininf = argv[1];//input
	std::string file_out_pidinf = argv[2];//output
	double pion_thr = std::stod(argv[3]);//pion threshold
	//old ver(for ECC5)
	std::string filename = "M:\\data\\NINJA\\PID\\average.bin.root_fitresult.txt";
	std::string fileame_fit = "M:\\data\\NINJA\\PID\\average.bin.root_fitresult.txt_fit_mom.txt";

	if (argc != 4)  {// other ECC(exept for ecc5)
		filename = argv[4];//inout
		fileame_fit = argv[5];//input
	}

	//std::vector<Chain_inf> c_inf;
	std::vector<Chain_inf> c_inf;
	c_inf = Read_chainif(file_in_chaininf);
	double angle_min, angle_max;
	//testdata生成
	/*
	int cid = 0;
	for (int i_ang = 0; i_ang < 10; i_ang++) {
		iang_to_angle_range(i_ang, angle_min, angle_max);
		for (double pb = 1; pb < 1500; pb += 5) {
			for (double vph = 1; vph < 500; vph += 2) {
				Chain_inf c;
				c.angle = (angle_max + angle_min) / 2;
				c.chainid = cid;
				c.groupid = cid;
				cid++;
				c.pb = pb;
				c.vph = vph;
				c.pid = 0;
				c.likelihood_ratio = 0;
				c.pion_likelihood = 0;
				c.proton_likelihood = 0;
				c_inf.push_back(c);
			}
		}
	}
	*/

	std::ofstream ofs(file_out_pidinf);

	double pth_muon = 0.10 * m_mu * (0.10 * m_mu / sqrt(m_mu * m_mu + 0.10 * 0.10 * m_mu * m_mu));
	double pth_proton = 0.10 * m_p * (0.10 * m_p / sqrt(m_p * m_p + 0.10 * 0.10 * m_p * m_p));
	std::map<int, VPH_pb_func>VPH_PID_param;
	for (int i_ang = 0; i_ang < 20; i_ang++) {
		iang_to_angle_range(i_ang, angle_min, angle_max);
		printf("%d angle %.1lf - %.1lf\n", i_ang, angle_min, angle_max);
		//理論式のパラメータを読み込む
		double vph_mean_param[4] = {}, vph_sigma_param[2] = {};
		Read_Theoretical_formula(fileame_fit, angle_min, angle_max, vph_mean_param, vph_sigma_param);
		//printf("%.1lf %.1lf %.1lf %.1lf\n", vph_mean_param[0], vph_mean_param[1], vph_mean_param[2], vph_mean_param[3]);
		double thr_sigma = -1;
		int DataNum = 0;
		DataPoint data[30] = {};
		//ここがバグ
		Read_Data_point(filename, angle_min, angle_max, data, thr_sigma, DataNum);
		printf("Data Num %d\n", DataNum);

		VPH_pb_func pid_param;
		for (int i = 0; i < DataNum; i++) {
			double pb_mean = (data[i].input_mom_min + data[i].input_mom_max) / 2;
			pid_param.datapoint.insert(std::make_pair(pb_mean, data[i]));
		}
		pid_param.vph_sigma_thr = thr_sigma;
		for (int i = 0; i < 4; i++) {
			pid_param.vph_mean_param[i] = vph_mean_param[i];
			if (i < 2) {
				pid_param.vph_sigma_param[i] = vph_sigma_param[i];
			}
		}

		pid_param.SetCrossPoint_VPH();

		VPH_PID_param.insert(std::make_pair(i_ang, pid_param));

	}
	int count_particle[3] = {};

	int all = c_inf.size(), cnt = 0;
	for (auto& c : c_inf) {
		if (cnt % 10000 == 0) {
			printf("\rnow calc likelihood %d/%d(%4.1lf%%)", cnt, all, cnt * 100. / all);
		}
		cnt++;
		int i_ang = angle_to_iangle(c.angle);

		VPH_pb_func vph_param;
		if (VPH_PID_param.count(i_ang) == 0) {
			vph_param = VPH_PID_param.rbegin()->second;
		}
		else {
			vph_param = VPH_PID_param.find(i_ang)->second;
		}

		c.proton_likelihood = (vph_param.Calc_VPH_mean_proton(c.pb) - c.vph) / vph_param.Calc_VPH_sigma_proton(c.pb);

		c.pion_likelihood = (c.vph - vph_param.Calc_VPH_mean_pion(c.pb)) / vph_param.Calc_VPH_sigma_pion(c.pb);

		c.likelihood_ratio = c.pion_likelihood / (c.pion_likelihood + c.proton_likelihood);

		if (c.pid == 0) {
			if (c.pion_likelihood < 0) {
				c.likelihood_ratio = 0;
				c.pid = 211;
				count_particle[1]++;
			}
			else if (c.proton_likelihood < 0) {
				c.likelihood_ratio = 1;
				c.pid = 2212;
				count_particle[2]++;
			}
			else if (c.likelihood_ratio < pion_thr) {
				c.pid = 211;
				count_particle[1]++;
			}
			else {
				c.pid = 2212;
				count_particle[2]++;
			}
		}
		else {
			count_particle[0]++;
		}
		ofs << std::setw(10) << c.groupid << " " << c.chainid << " " << c.angle << " " << c.pb << " " << c.vph << " "
			<< c.proton_likelihood << " " << c.pion_likelihood << " " << c.likelihood_ratio << " " << c.pid << std::endl;

	}
	printf("now calc likelihood %d/%d(%4.1lf%%)\n", cnt, all, cnt * 100. / all);

	printf("all:%d\n", c_inf.size());
	printf("\tmuon   :%d\n", count_particle[0]);
	printf("\tpion   :%d\n", count_particle[1]);
	printf("\tproton :%d\n", count_particle[2]);
}
std::vector<Chain_inf> Read_chainif(std::string filename) {
	std::vector<Chain_inf> ret;
	Chain_inf c_tmp;
	int num = 0;
	std::ifstream ifs(filename.c_str());
	while (ifs >> c_tmp.groupid >> c_tmp.chainid >> c_tmp.pid >>
		c_tmp.angle >> c_tmp.pb >> c_tmp.vph) {
		if (num % 100000 == 0) {
			printf("\r read chain %d", num);
		}
		num++;

		if (c_tmp.pb < 0)c_tmp.pb = 10;
		ret.push_back(c_tmp);
	}
	printf("\r read chain %d\n", num);

	ifs.close();
	return ret;
}

void iang_to_angle_range(int i_ang, double& angle_min, double& angle_max) {
	if (i_ang < 7) {
		angle_min = i_ang * 0.1;
		angle_max = (i_ang + 1) * 0.1;
	}
	else if (i_ang < 11) {
		angle_min = (i_ang - 7) * 0.2 + 0.7;
		angle_max = (i_ang - 7 + 1) * 0.2 + 0.7;
	}
	else if (i_ang < 15) {
		angle_min = (i_ang - 11) * 0.4 + 1.5;
		angle_max = (i_ang - 11 + 1) * 0.4 + 1.5;
	}
	else {
		angle_min = (i_ang - 15) * 0.6 + 3.1;
		angle_max = (i_ang - 15 + 1) * 0.6 + 3.1;
	}
}
int angle_to_iangle(double angle) {
	if (angle < 0.7) {
		return int(angle / 0.1 + 0.1);
	}
	else if (angle < 1.5) {
		return int((angle - 0.7) / 0.2 + 0.1) + 7;
	}
	else if (angle < 3.1) {
		return int((angle - 1.5) / 0.4 + 0.1) + 11;
	}
	else {
		return int((angle - 3.1) / 0.6 + 0.1) + 15;
	}
	return -1;
}
//理論式のパラメータを読み込む
void Read_Theoretical_formula(std::string filename, double& angle_min, double& angle_max, double* vph_mean_param, double* vph_sigma_param) {

	std::ifstream ifs(filename.c_str());
	double input_ang_min, input_ang_max, buf[10];
	//double vph_mean_param[4] = {}, vph_sigma_param[2] = {};
	//mass
	vph_mean_param[0] = m_p;
	//charge(abs)
	vph_mean_param[1] = 1;
	vph_mean_param[2] = 0;
	vph_mean_param[3] = 0;
	double all_angle_max = 0;
	while (ifs >> input_ang_min >> input_ang_max
		>> buf[0] >> buf[1] >> buf[2] >> buf[3] >> buf[4] >> buf[5] >> buf[6] >> buf[7]) {
		if (all_angle_max < (input_ang_min + input_ang_max) / 2)all_angle_max = (input_ang_min + input_ang_max) / 2;
		if ((input_ang_min + input_ang_max) / 2 < angle_min)continue;
		if ((input_ang_min + input_ang_max) / 2 > angle_max)continue;

		//if (all_angle_max < (input_ang_min + input_ang_max) / 2)all_angle_max = (input_ang_min + input_ang_max) / 2;
		//if ((input_ang_min + input_ang_max) / 2 < angle_min)continue;
		//if ((input_ang_min + input_ang_max) / 2 > angle_max)continue;
		vph_mean_param[2] = buf[0];
		vph_mean_param[3] = buf[2];
		vph_sigma_param[0] = buf[4];
		vph_sigma_param[1] = buf[6];
		break;
	}
	ifs.close();
	if (vph_mean_param[2] < 0.01) {
		ifs.open(filename.c_str());
		while (ifs >> input_ang_min >> input_ang_max
			>> buf[0] >> buf[1] >> buf[2] >> buf[3] >> buf[4] >> buf[5] >> buf[6] >> buf[7]) {
			if (all_angle_max > input_ang_max)continue;
			vph_mean_param[2] = buf[0];
			vph_mean_param[3] = buf[2];
			vph_sigma_param[0] = buf[4];
			vph_sigma_param[1] = buf[6];
			break;
		}
		ifs.close();
	}

}
//各点のデータを読み込む
void Read_Data_point(std::string filename, double& angle_min, double& angle_max, DataPoint* data, double& thr_sigma, int& DataNum) {
	std::ifstream ifs(filename.c_str());
	double input_mom_min, input_mom_max;
	double input_ang_min, input_ang_max;
	double mean[3], mean_low[3], mean_hi[3], sigma[3], sigma_low[3], sigma_hi[3];
	//Double_t thr_sigma = -1;
	thr_sigma = -1;
	//int DataNum = 0;
	DataNum = 0;
	double all_angle_max = 0;
	double input_ang_max_max = 0;
	while (ifs >> input_ang_min >> input_ang_max >> input_mom_min >> input_mom_max
		>> mean[0] >> mean_low[0] >> mean_hi[0] >> sigma[0] >> sigma_low[0] >> sigma_hi[0]
		>> mean[1] >> mean_low[1] >> mean_hi[1] >> sigma[1] >> sigma_low[1] >> sigma_hi[1]
		>> mean[2] >> mean_low[2] >> mean_hi[2] >> sigma[2] >> sigma_low[2] >> sigma_hi[2]) {
		//if (input_ang_max_max < input_ang_max)input_ang_max_max = input_ang_max;
		//if (input_ang_min > c.angle)continue;
		//if (input_ang_max <= c.angle)continue;

		if (all_angle_max < (input_ang_min + input_ang_max) / 2)all_angle_max = (input_ang_min + input_ang_max) / 2;
		if ((input_ang_min + input_ang_max) / 2 < angle_min)continue;
		if ((input_ang_min + input_ang_max) / 2 > angle_max)continue;
		data[DataNum].input_ang_min = input_ang_min;
		data[DataNum].input_ang_max = input_ang_max;
		data[DataNum].input_mom_min = input_mom_min;
		data[DataNum].input_mom_max = input_mom_max;
		data[DataNum].mean[0] = mean[0];
		data[DataNum].mean[1] = mean[1];
		data[DataNum].mean[2] = mean[2];
		data[DataNum].mean_low[0] = mean_low[0];
		data[DataNum].mean_low[1] = mean_low[1];
		data[DataNum].mean_low[2] = mean_low[2];
		data[DataNum].mean_hi[0] = mean_hi[0];
		data[DataNum].mean_hi[1] = mean_hi[1];
		data[DataNum].mean_hi[2] = mean_hi[2];
		data[DataNum].sigma[0] = sigma[0];
		data[DataNum].sigma[1] = sigma[1];
		data[DataNum].sigma[2] = sigma[2];
		data[DataNum].sigma_low[0] = sigma_low[0];
		data[DataNum].sigma_low[1] = sigma_low[1];
		data[DataNum].sigma_low[2] = sigma_low[2];
		data[DataNum].sigma_hi[0] = sigma_hi[0];
		data[DataNum].sigma_hi[1] = sigma_hi[1];
		data[DataNum].sigma_hi[2] = sigma_hi[2];
		DataNum++;
		if (1200 < (input_mom_min + input_mom_max) / 2 && (input_mom_min + input_mom_max) / 2 < 1300) {
			thr_sigma = sigma[1];
		}
	}
	ifs.close();
	if (DataNum == 0) {
		ifs.open(filename.c_str());
		while (ifs >> input_ang_min >> input_ang_max >> input_mom_min >> input_mom_max
			>> mean[0] >> mean_low[0] >> mean_hi[0] >> sigma[0] >> sigma_low[0] >> sigma_hi[0]
			>> mean[1] >> mean_low[1] >> mean_hi[1] >> sigma[1] >> sigma_low[1] >> sigma_hi[1]
			>> mean[2] >> mean_low[2] >> mean_hi[2] >> sigma[2] >> sigma_low[2] >> sigma_hi[2]) {
			if (all_angle_max > input_ang_max)continue;
			//if (c.angle < input_ang_max)continue;
			//if (input_ang_max_max - 0.1 > input_ang_max)continue;

			data[DataNum].input_ang_min = input_ang_min;
			data[DataNum].input_ang_max = input_ang_max;
			data[DataNum].input_mom_min = input_mom_min;
			data[DataNum].input_mom_max = input_mom_max;
			data[DataNum].mean[0] = mean[0];
			data[DataNum].mean[1] = mean[1];
			data[DataNum].mean[2] = mean[2];
			data[DataNum].mean_low[0] = mean_low[0];
			data[DataNum].mean_low[1] = mean_low[1];
			data[DataNum].mean_low[2] = mean_low[2];
			data[DataNum].mean_hi[0] = mean_hi[0];
			data[DataNum].mean_hi[1] = mean_hi[1];
			data[DataNum].mean_hi[2] = mean_hi[2];
			data[DataNum].sigma[0] = sigma[0];
			data[DataNum].sigma[1] = sigma[1];
			data[DataNum].sigma[2] = sigma[2];
			data[DataNum].sigma_low[0] = sigma_low[0];
			data[DataNum].sigma_low[1] = sigma_low[1];
			data[DataNum].sigma_low[2] = sigma_low[2];
			data[DataNum].sigma_hi[0] = sigma_hi[0];
			data[DataNum].sigma_hi[1] = sigma_hi[1];
			data[DataNum].sigma_hi[2] = sigma_hi[2];
			DataNum++;
			if (1200 < (input_mom_min + input_mom_max) / 2 && (input_mom_min + input_mom_max) / 2 < 1300) {
				thr_sigma = sigma[1];
			}
		}
		ifs.close();

	}

	if (thr_sigma < 0) {
		printf("sigma threshold not found\n");
		exit(1);
	}
}

double Calc_momentum_PHV_emulsion_Fit(double var, double* par) {

	double pbeta = var;

	double dEdx = 0.0;


	for (int i = 0; i < 12; i++) {

		switch (i) {
		case  0: Absorber = Ag; break;
		case  1: Absorber = Br; break;
		case  2: Absorber = I;  break;
		case  3: Absorber = C;  break;
		case  4: Absorber = N;  break;
		case  5: Absorber = O;  break;
		case  6: Absorber = H;  break;
		case  7: Absorber = S;  break;
		case  8: Absorber = Na; break;
		case  9: Absorber = Fe; break;
		case 10: Absorber = Au; break;
		case 11: Absorber = Cl; break;
		default: break;
		}


		double mass = par[0];
		double z = par[1];

		// information of constants //
		double me = m_e;                   // mass of electron [GeV]
		double mec2 = 0.510998928 / 1000.0;  // electron mass * c^2 [GeV]
		double K = 0.307075;               // 4πNA re^2 me c^2 [MeV mol^-1 cm^2]

		// information of absorber (Fe) //
		double Z = mtrl_Z[Absorber];              // atomic number of absorber           (PDG2017)
		double A = mtrl_A[Absorber];              // atomic mass of absorber [g mol^-1]  (PDG2017)
		double I = mtrl_I[Absorber] / 1000000000.0; // Mean excitation energy of Fe [GeV]  (PDG2017)
		double density = mtrl_Density[Absorber];  // density of Fe [g cm^-3]             (PDG2017)


		/* Calculation of dE/dx */

		// Energy, Beta, Lorentz factor //
		double E = (pbeta + sqrt(pbeta * pbeta + 4.0 * mass * mass)) / 2.0;   // Energy [GeV]
		double p = sqrt(pbeta * E);
		double beta = p / E;                         // Beta
		double gamma = 1.0 / sqrt(1.0 - beta * beta);  // Lorentz factor

		// Maximum Energy Transfer to an electron in a single collision //
		double Wmax = (2.0 * mec2 * beta * beta * gamma * gamma) / (1.0 + 2.0 * gamma * me / mass + (me / mass) * (me / mass)); // maximum energy transfer to an electron in a single collision

		// density effect //
		double  C_bar = mtrl_C_bar[Absorber];
		double      a = mtrl_a[Absorber];
		double      m = mtrl_m[Absorber];
		double     x1 = mtrl_x1[Absorber];
		double     x0 = mtrl_x0[Absorber];
		double delta0 = mtrl_delta0[Absorber];

		double x = log10(p / mass);
		double density_effect = 0.0;

		if (x >= x1)
			density_effect = (2.0 * log(10) * x - C_bar) / 2.0;
		else if (x0 <= x && x < x1)
			density_effect = (2.0 * log(10) * x - C_bar + a * pow(x1 - x, m)) / 2.0;
		else if (x < x0)
			density_effect = (delta0 * pow(10.0, 2.0 * (x - x0))) / 2.0; // conductors
			//density_effect = 0; // nonconductors // 180424 modified


		// weight ratio //
		double wratio = mtrl_wratio[Absorber];

		dEdx += wratio * (K * z * z * (Z / A) * (1.0 / (beta * beta)) * ((1.0 / 2.0) * log(2.0 * mec2 * beta * beta * gamma * gamma * Wmax / (I * I)) - beta * beta - density_effect));

	}


	double slope = par[2];     // 傾き
	double intercept = par[3]; // 切片

	return slope * dEdx + intercept;

}

double VPH_pb_func::Calc_VPH_mean_proton(double pb) {
	double ret_vph = 0;
	//各点データを使用
	if (pb < cross_point_pb) {
		DataPoint data[2];
		data[0] = datapoint.lower_bound(pb)->second;
		if (datapoint.lower_bound(pb) == datapoint.begin()) {
			data[1] = std::next(datapoint.lower_bound(pb), 1)->second;
		}
		else {
			data[1] = std::next(datapoint.lower_bound(pb), -1)->second;
		}
		double pb_mean[2];
		pb_mean[0] = (data[0].input_mom_min + data[0].input_mom_max) / 2;
		pb_mean[1] = (data[1].input_mom_min + data[1].input_mom_max) / 2;
		//data[0],[1]で直線を作って外挿
		ret_vph = (data[0].mean[2] - data[1].mean[2]) / (pb_mean[0] - pb_mean[1]) * (pb - pb_mean[0]) + data[0].mean[2];
	}
	else {
		double dE = Calc_momentum_PHV_emulsion_Fit(pb / 1000, vph_mean_param);
		ret_vph = dE;
	}

	return ret_vph;
}
double VPH_pb_func::Calc_VPH_sigma_proton(double pb) {
	double vph = Calc_momentum_PHV_emulsion_Fit(pb / 1000, vph_mean_param);
	double vph_thr = Calc_momentum_PHV_emulsion_Fit(cross_point_pb, vph_mean_param);
	double ret_sigma;

	if (pb < 200) {
		//VPH��mean thr���傫��-->��pb�̃f�[�^
		DataPoint data[2];
		if (datapoint.rbegin()->second.input_mom_min < pb) {
			data[0] = datapoint.rbegin()->second;
			data[1] = std::next(datapoint.rbegin(), 1)->second;
		}
		else {
			data[0] = datapoint.lower_bound(pb)->second;
			if (datapoint.lower_bound(pb) == datapoint.begin()) {
				data[1] = std::next(datapoint.lower_bound(pb), 1)->second;
			}
			else {
				data[1] = std::next(datapoint.lower_bound(pb), -1)->second;
			}
		}

		double pb_mean[2];
		pb_mean[0] = (data[0].input_mom_min + data[0].input_mom_max) / 2;
		pb_mean[1] = (data[1].input_mom_min + data[1].input_mom_max) / 2;
		//data[0],[1]で直線を作って外挿
		ret_sigma = (data[0].sigma[2] - data[1].sigma[2]) / (pb_mean[0] - pb_mean[1]) * (pb - pb_mean[0]) + data[0].sigma[2];
		return ret_sigma;

	}
	ret_sigma = vph_sigma_param[0] + vph_sigma_param[1] * sqrt(vph);
	return std::max(vph_sigma_thr, ret_sigma);

}

//ここを書く
void VPH_pb_func::SetCrossPoint_VPH() {
	cross_point_pb = -1;
	double pb_min = 1;
	double pb_max = 1200;
	double pich = 1;
	double vph_bete[2] = {};
	double vph_data[2] = {};
	int sign = 0;
	for (double pb = pb_min; pb <= pb_max; pb += pich) {
		vph_bete[0] = Calc_momentum_PHV_emulsion_Fit(pb / 1000, vph_mean_param);
		vph_bete[1] = Calc_momentum_PHV_emulsion_Fit((pb + pich) / 1000, vph_mean_param);
		DataPoint data[2];
		data[0] = datapoint.lower_bound(pb)->second;
		if (datapoint.lower_bound(pb) == datapoint.begin()) {
			data[1] = std::next(datapoint.lower_bound(pb), 1)->second;
		}
		else {
			data[1] = std::next(datapoint.lower_bound(pb), -1)->second;
		}
		double pb_mean[2];
		pb_mean[0] = (data[0].input_mom_min + data[0].input_mom_max) / 2;
		pb_mean[1] = (data[1].input_mom_min + data[1].input_mom_max) / 2;

		vph_data[0] = (data[0].mean[2] - data[1].mean[2]) / (pb_mean[0] - pb_mean[1]) * (pb - pb_mean[0]) + data[0].mean[2];
		vph_data[1] = (data[0].mean[2] - data[1].mean[2]) / (pb_mean[0] - pb_mean[1]) * (pb + pich - pb_mean[0]) + data[0].mean[2];

		//printf("%.1lf-%.1lf %.1lf %.1lf %.1lf %.1lf\n", pb, pb + pich, vph_bete[0], vph_bete[1], vph_data[0], vph_data[1]);
		//交点条件
		if (vph_bete[0] >= vph_data[0] && vph_bete[1] < vph_data[1]) {
			cross_point_pb = pb + std::min(pich, pich * (vph_bete[0] - vph_data[0]) / ((vph_data[1] - vph_bete[1]) + (vph_bete[0] - vph_data[0])));
			return;
		}
	}
	fprintf(stderr, "VPH corss point not found\n");
}

double VPH_pb_func::Calc_VPH_mean_pion(double pb) {
	double ret_vph = 0;

	DataPoint data[2];
	if (datapoint.rbegin()->second.input_mom_min < pb) {
		data[0] = datapoint.rbegin()->second;
		data[1] = std::next(datapoint.rbegin(), 1)->second;
	}
	else {
		data[0] = datapoint.lower_bound(pb)->second;
		if (datapoint.lower_bound(pb) == datapoint.begin()) {
			data[1] = std::next(datapoint.lower_bound(pb), 1)->second;
		}
		else {
			data[1] = std::next(datapoint.lower_bound(pb), -1)->second;
		}
	}
	double pb_mean[2];
	pb_mean[0] = (data[0].input_mom_min + data[0].input_mom_max) / 2;
	pb_mean[1] = (data[1].input_mom_min + data[1].input_mom_max) / 2;
	//data[0],[1]で直線を作って外挿
	ret_vph = (data[0].mean[1] - data[1].mean[1]) / (pb_mean[0] - pb_mean[1]) * (pb - pb_mean[0]) + data[0].mean[1];
	return ret_vph;
}
double VPH_pb_func::Calc_VPH_sigma_pion(double pb) {
	double ret_sigma = 0;

	DataPoint data[2];
	if (datapoint.rbegin()->second.input_mom_min < pb) {
		data[0] = datapoint.rbegin()->second;
		data[1] = std::next(datapoint.rbegin(), 1)->second;
	}
	else {
		data[0] = datapoint.lower_bound(pb)->second;
		if (datapoint.lower_bound(pb) == datapoint.begin()) {
			data[1] = std::next(datapoint.lower_bound(pb), 1)->second;
		}
		else {
			data[1] = std::next(datapoint.lower_bound(pb), -1)->second;
		}
	}
	double pb_mean[2];
	pb_mean[0] = (data[0].input_mom_min + data[0].input_mom_max) / 2;
	pb_mean[1] = (data[1].input_mom_min + data[1].input_mom_max) / 2;
	//data[0],[1]で直線を作って外挿
	ret_sigma = (data[0].sigma[1] - data[1].sigma[1]) / (pb_mean[0] - pb_mean[1]) * (pb - pb_mean[0]) + data[0].sigma[1];
	return ret_sigma;
}
