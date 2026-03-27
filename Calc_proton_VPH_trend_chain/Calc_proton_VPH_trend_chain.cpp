#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

#include "M:\data\NINJA\PID\PDG.h"
#include "M:\data\NINJA\PID\AbsorberMaterial_Emulsion_NINJA_Run6.h"
class DataPointParam {
public:
	double input_ang_min, input_ang_max, thr_sigma, input_mom_min, input_mom_max;
	double mean[3], mean_low[3], mean_hi[3], sigma[3], sigma_low[3], sigma_hi[3];
};
class DataPoint {
public:
	double input_ang_min, input_ang_max, thr_sigma;
	std::map<double, DataPointParam> data;
};
class VPH_dEdx_Param {
public:
	double input_ang_min, input_ang_max;
	double vph_mean_param[4] = {}, vph_sigma_param[2] = {};
};

bool sort_VPH_dEdx_Param(const VPH_dEdx_Param& left, const VPH_dEdx_Param& right) {
	return left.input_ang_max < right.input_ang_max;
}
bool sort_DataPoint(const DataPoint& left, const DataPoint& right) {
	return left.input_ang_max < right.input_ang_max;
}

class output_inf {
public:
	int groupid, chainid, num, pl;
	double angle, mom, vph;
};
class output_inf_chi2 {
public:
	int groupid, chainid, num, pl;
	double angle, mom, vph, chi2;
};
std::vector<VPH_dEdx_Param> Read_Theoretical_formula(std::string filename);
std::vector<DataPoint>Read_Data_point(std::string filename);

double Calc_angle(Momentum_recon::Mom_chain& c);

VPH_dEdx_Param Select_VPH_dEdx_Param(std::vector<VPH_dEdx_Param>& param, double angle);
DataPointParam Select_DataPoint(std::vector<DataPoint>& data, double angle, double mom);
DataPoint Select_DataPoint(std::vector<DataPoint>& data, double angle);

double search_vph_mean_cross_point(VPH_dEdx_Param& par, DataPoint& data);
double search_vph_sigma_cross_point(VPH_dEdx_Param& par, DataPoint& data, double vph_start);

bool judge_direction(Momentum_recon::Mom_chain& c);

void Calc_VPH_trend(std::vector<output_inf>& out, int groupid, Momentum_recon::Mom_chain& c, VPH_dEdx_Param& par);
void Calc_VPH_trend_inv(std::vector<output_inf>& out, int groupid, Momentum_recon::Mom_chain& c, VPH_dEdx_Param& par);

void Output_file(std::string filename, std::vector<output_inf>& out);
Momentum_recon::Event_information Momch_chain_sel(int eventid, int chainid, std::vector<Momentum_recon::Event_information>& momch_all, std::vector<Momentum_recon::Event_information>& momch_divide);

int main(int argc, char** argv) {
	if (argc < 6 || argc>7) {
		fprintf(stderr, "uasge:input-momch input-momch eventid chainid output-file [input_param_path = M:\\data\\NINJA\\PID ]\ninput_param is ...\n\taverage.bin.root_fitresult.txt\n\taverage.bin.root_fitresult.txt_fit_mom.txt");
		exit(1);
	}
	std::string file_in_momch_all = argv[1];
	std::string file_in_momch_divide = argv[2];
	int eventid = std::stoi(argv[3]);
	int chainid = std::stoi(argv[4]);
	std::string file_out = argv[5];
	std::string param_path = "M:\\data\\NINJA\\PID";
	if (argc == 7) {
		param_path = argv[6];
	}

	std::vector<Momentum_recon::Event_information> momch_all = Momentum_recon::Read_Event_information_extension(file_in_momch_all);
	std::vector<Momentum_recon::Event_information> momch_divide = Momentum_recon::Read_Event_information_extension(file_in_momch_divide);


	Momentum_recon::Event_information momch_sel = Momch_chain_sel(eventid, chainid, momch_all, momch_divide);

	std::string filename = param_path + "\\average.bin.root_fitresult.txt";
	std::string fileame_fit = param_path + "\\average.bin.root_fitresult.txt_fit_mom.txt";

	double pth_muon = 0.10 * m_mu * (0.10 * m_mu / sqrt(m_mu * m_mu + 0.10 * 0.10 * m_mu * m_mu));
	double pth_proton = 0.10 * m_p * (0.10 * m_p / sqrt(m_p * m_p + 0.10 * 0.10 * m_p * m_p));

	std::vector<VPH_dEdx_Param> Theoretical_formula = Read_Theoretical_formula(fileame_fit);
	std::vector<DataPoint>Data_point = Read_Data_point(filename);

	std::vector<output_inf> out;
	int64_t all = 0, cnt = 0;
	double angle = Calc_angle(momch_sel.chains[0]);

	for (auto& c : momch_sel.chains) {


		VPH_dEdx_Param sel_param = Select_VPH_dEdx_Param(Theoretical_formula, angle);
		DataPoint sel_data = Select_DataPoint(Data_point, angle);

		double mean_p_thr = search_vph_mean_cross_point(sel_param, sel_data);
		double sigma_p_thr = search_vph_sigma_cross_point(sel_param, sel_data, mean_p_thr);

		if (c.direction == -1) {
			Calc_VPH_trend(out, eventid, c, sel_param);
		}
		else if (c.direction == 1) {
			Calc_VPH_trend_inv(out, eventid, c, sel_param);
		}
	}


	Output_file(file_out, out);

}


Momentum_recon::Event_information Momch_chain_sel(int eventid, int chainid, std::vector<Momentum_recon::Event_information>& momch_all, std::vector<Momentum_recon::Event_information>& momch_divide) {

	Momentum_recon::Event_information ret;

	for (auto itr = momch_all.begin(); itr != momch_all.end(); itr++) {
		if (itr->groupid != eventid)continue;
		ret = *itr;
		ret.chains.clear();
		for (auto& c : itr->chains) {
			if (c.chainid / 10 != chainid)continue;
			ret.chains.push_back(c);
		}
	}
	for (auto itr = momch_divide.begin(); itr != momch_divide.end(); itr++) {
		if (itr->groupid != eventid)continue;
		for (auto& c : itr->chains) {
			if (c.chainid / 10 != chainid)continue;
			bool flg = false;
			for (auto& c2 : ret.chains) {
				if (c2.chainid == c.chainid) {
					flg = true;
				}
			}
			if (!flg) {
				ret.chains.push_back(c);
			}
		}


	}
	return ret;

}



//理論式のパラメータを読み込む
std::vector<VPH_dEdx_Param> Read_Theoretical_formula(std::string filename) {

	std::ifstream ifs(filename.c_str());
	std::vector<VPH_dEdx_Param> ret;
	VPH_dEdx_Param p;
	double buf[10];
	while (ifs >> p.input_ang_min >> p.input_ang_max
		>> buf[0] >> buf[1] >> buf[2] >> buf[3] >> buf[4] >> buf[5] >> buf[6] >> buf[7]) {
		//mass
		p.vph_mean_param[0] = m_p;
		//charge(abs)
		p.vph_mean_param[1] = 1;
		p.vph_mean_param[2] = 0;
		p.vph_mean_param[3] = 0;

		p.vph_mean_param[2] = buf[0];
		p.vph_mean_param[3] = buf[2];
		p.vph_sigma_param[0] = buf[4];
		p.vph_sigma_param[1] = buf[6];
		ret.push_back(p);
	}
	if (false) {
		for (auto& para : ret) {
			printf("%.2lf %.2lf %.1lf %.1lf %.1lf %.1lf\n"
				, para.input_ang_min, para.input_ang_max, para.vph_mean_param[2], para.vph_mean_param[3]
				, para.vph_sigma_param[0], para.vph_sigma_param[1]);
		}
	}
	std::sort(ret.begin(), ret.end(), sort_VPH_dEdx_Param);
	return ret;
}
//各点のデータを読み込む
std::vector<DataPoint>Read_Data_point(std::string filename) {
	std::ifstream ifs(filename.c_str());
	std::vector<DataPoint> ret;
	double input_ang_min, input_ang_max, angle;
	DataPointParam data;
	DataPoint* d_p = NULL;
	bool flg = false;
	while (ifs >> input_ang_min >> input_ang_max >> data.input_mom_min >> data.input_mom_max
		>> data.mean[0] >> data.mean_low[0] >> data.mean_hi[0] >> data.sigma[0] >> data.sigma_low[0] >> data.sigma_hi[0]
		>> data.mean[1] >> data.mean_low[1] >> data.mean_hi[1] >> data.sigma[1] >> data.sigma_low[1] >> data.sigma_hi[1]
		>> data.mean[2] >> data.mean_low[2] >> data.mean_hi[2] >> data.sigma[2] >> data.sigma_low[2] >> data.sigma_hi[2]) {
		flg = false;
		angle = (input_ang_min + input_ang_max) / 2;
		data.input_ang_min = input_ang_min;
		data.input_ang_max = input_ang_max;
		for (auto& d : ret) {
			if (angle < d.input_ang_min)continue;
			if (d.input_ang_max <= angle)continue;
			flg = true;
			d_p = &d;
		}
		if (!flg) {
			DataPoint data_tmp;
			data_tmp.data.insert(std::make_pair(data.input_mom_max, data));
			data_tmp.input_ang_min = input_ang_min;
			data_tmp.input_ang_max = input_ang_max;
			data_tmp.thr_sigma = -1;
			ret.push_back(data_tmp);
			for (auto& d : ret) {
				if (angle < d.input_ang_min)continue;
				if (d.input_ang_max <= angle)continue;
				flg = true;
				d_p = &d;
			}
			if (!flg) {
				fprintf(stderr, "exception\n");
				exit(1);
			}
		}
		else {
			d_p->data.insert(std::make_pair(data.input_mom_max, data));
		}

		if (1200 < (data.input_mom_min + data.input_mom_max) / 2 && (data.input_mom_min + data.input_mom_max) / 2 < 1300) {
			d_p->thr_sigma = data.sigma[1];
		}
	}

	if (false) {
		for (auto& para : ret) {
			printf("%.2lf %.2lf %.1lf %d\n"
				, para.input_ang_min, para.input_ang_max, para.thr_sigma, para.data.size());
			for (auto& d : para.data) {
				printf("\t%.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf %.1lf\n"
					, d.second.input_mom_min, d.second.input_mom_max, d.second.mean[0], d.second.mean[1], d.second.mean[2]
					, d.second.sigma[0], d.second.sigma[1], d.second.sigma[2]);
			}
		}
	}
	for (auto& p : ret) {
		double thr_sigma = p.thr_sigma;
		for (auto& d : p.data) {
			d.second.thr_sigma = thr_sigma;
		}
	}
	std::sort(ret.begin(), ret.end(), sort_DataPoint);
	return ret;
}

double Calc_angle(Momentum_recon::Mom_chain& c) {
	int num = 0;

	double ax = 0, ay = 0, anlge = 0;
	for (auto& b : c.base) {
		ax += b.ax;
		ay += b.ay;
		num++;
	}
	ax /= num;
	ay /= num;
	return sqrt(ax * ax + ay * ay);
}
VPH_dEdx_Param Select_VPH_dEdx_Param(std::vector<VPH_dEdx_Param>& param, double angle) {
	VPH_dEdx_Param ret;
	bool flg = false;
	for (auto& p : param) {
		if (angle < p.input_ang_min)continue;
		ret = p;
		flg = true;
	}
	if (!flg) {
		ret = *param.begin();
	}
	return ret;
}
DataPointParam Select_DataPoint(std::vector<DataPoint>& data, double angle, double mom) {
	DataPoint ret;
	bool flg = false;
	for (auto& p : data) {
		if (angle < p.input_ang_min)continue;
		ret = p;
		flg = true;
	}
	if (!flg) {
		ret = *data.begin();
	}

	DataPointParam ret2;
	flg = false;
	for (auto& p : ret.data) {
		if (mom < p.second.input_mom_min)continue;
		ret2 = p.second;
		flg = true;
	}
	if (!flg) {
		ret2 = ret.data.begin()->second;
	}
	return ret2;

}
DataPoint Select_DataPoint(std::vector<DataPoint>& data, double angle) {
	DataPoint ret;
	bool flg = false;
	for (auto& p : data) {
		if (angle < p.input_ang_min)continue;
		ret = p;
		flg = true;
	}
	if (!flg) {
		ret = *data.begin();
	}
	return ret;
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

double search_vph_mean_cross_point(VPH_dEdx_Param& par, DataPoint& data) {
	double mean_p_thr;
	double pb_min = 1;
	double pb_max = 300;
	double pich = 1;
	double vph = 0;
	int sign = 0;
	for (double pb = pb_min; pb <= pb_max; pb += pich) {
		vph = Calc_momentum_PHV_emulsion_Fit(pb / 1000, par.vph_mean_param);
		//ここで決めている
		if (pb >= 200) {
			mean_p_thr = vph;
			sign = -1;
			break;
		}
		else continue;


		double vph_data = 0;
		for (auto& d : data.data) {
			if (d.second.input_mom_max > pb) {
				vph_data = d.second.mean[2];
				break;
			}
		}

		//printf("pb:%lf %lf %lf\n",pb, vph_data, vph);
		if (vph_data < 1)continue;
		if (sign == 0) {
			if (vph - vph_data >= 0)sign = 1;
			else sign = -1;
		}
		else {
			if (sign * (vph - vph_data) < 0) {
				mean_p_thr = vph;
				break;
			}
		}
	}
	if (sign == 0) {
		printf("cross point not found\n");
		exit(1);
	}
	return mean_p_thr;
}
double search_vph_sigma_cross_point(VPH_dEdx_Param& par, DataPoint& data, double vph_start) {
	double sigma_p_thr = -1;
	double pich = 1;
	double vph = 0;
	double sigma = 0;
	for (double pb = vph_start; pb < 2000; pb += pich) {
		vph = Calc_momentum_PHV_emulsion_Fit(pb / 1000, par.vph_mean_param);
		sigma = par.vph_sigma_param[0] + par.vph_sigma_param[1] * sqrt(vph);
		if (sigma - data.thr_sigma < 1) {
			sigma_p_thr = vph;
			break;
		}
	}
	if (sigma_p_thr < 0) {
		printf("cross point not found\n");
		exit(1);
	}
	return sigma_p_thr;


}



void Calc_Energyloss_water(double& pbeta_mev, double angle, double* par, double thickness, double pbmin) {
	double density = 1;
	double pbeta = pbeta_mev / 1000;
	double mass = par[0];
	double z = par[1];
	double E_ini = (pbeta + sqrt(pbeta * pbeta + 4.0 * mass * mass)) / 2.0;   // Energy [GeV]

	double dEdx = 0.0;
	for (int i = 0; i < 2; i++) {

		switch (i) {
		case  0: Absorber = O;  break;
		case  1: Absorber = H;  break;
		default: break;
		}



		// information of constants //
		double me = m_e;                   // mass of electron [GeV]
		double mec2 = 0.510998928 / 1000.0;  // electron mass * c^2 [GeV]
		double K = 0.307075;               // 4πNA re^2 me c^2 [MeV mol^-1 cm^2]

		// information of absorber (Fe) //
		double Z = mtrl_Z[Absorber];              // atomic number of absorber           (PDG2017)
		double A = mtrl_A[Absorber];              // atomic mass of absorber [g mol^-1]  (PDG2017)
		double I = mtrl_I[Absorber] / 1000000000.0; // Mean excitation energy of Fe [GeV]  (PDG2017)
		//double density = mtrl_Density[Absorber];  // density of Fe [g cm^-3]             (PDG2017)


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
		double wratio = 0;
		if (Absorber == O)wratio = mtrl_A[O] / (mtrl_A[O] + 2 * mtrl_A[H]);
		else if (Absorber == H)wratio = 2 * mtrl_A[H] / (mtrl_A[O] + 2 * mtrl_A[H]);

		dEdx += wratio * (K * z * z * (Z / A) * (1.0 / (beta * beta)) * ((1.0 / 2.0) * log(2.0 * mec2 * beta * beta * gamma * gamma * Wmax / (I * I)) - beta * beta - density_effect));

	}
	// printf("water:%lf - %lf %lf\n", E_ini, dEdx * thickness * sqrt(angle*angle + 1), (E_ini*E_ini - mass * mass) / E_ini);
	 //printf("water:%.2lf - %.2lf --> ", E_ini*1000, dEdx* density * thickness * sqrt(angle*angle + 1));
	E_ini -= dEdx * density * thickness * sqrt(angle * angle + 1) / 1000;
	//printf("%.2lf \n", E_ini * 1000);
	pbeta = (E_ini * E_ini - mass * mass) / E_ini;
	pbeta_mev = pbeta * 1000;
	//printf("%lf\n", pbeta_mev);
	if (pbeta_mev < pbmin)	pbeta_mev = pbmin;

}
void Calc_Energyloss_iron(double& pbeta_mev, double angle, double* par, double thickness, double pbmin) {
	double density = 7.84;
	double pbeta = pbeta_mev / 1000;
	double mass = par[0];
	double z = par[1];
	double E_ini = (pbeta + sqrt(pbeta * pbeta + 4.0 * mass * mass)) / 2.0;   // Energy [GeV]

	double dEdx = 0.0;
	for (int i = 0; i < 1; i++) {

		switch (i) {
		case  0: Absorber = Fe;  break;
		default: break;
		}



		// information of constants //
		double me = m_e;                   // mass of electron [GeV]
		double mec2 = 0.510998928 / 1000.0;  // electron mass * c^2 [GeV]
		double K = 0.307075;               // 4πNA re^2 me c^2 [MeV mol^-1 cm^2]

		// information of absorber (Fe) //
		double Z = mtrl_Z[Absorber];              // atomic number of absorber           (PDG2017)
		double A = mtrl_A[Absorber];              // atomic mass of absorber [g mol^-1]  (PDG2017)
		double I = mtrl_I[Absorber] / 1000000000.0; // Mean excitation energy of Fe [GeV]  (PDG2017)
		//double density = mtrl_Density[Absorber];  // density of Fe [g cm^-3]             (PDG2017)


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
		double wratio = 1;

		dEdx += wratio * (K * z * z * (Z / A) * (1.0 / (beta * beta)) * ((1.0 / 2.0) * log(2.0 * mec2 * beta * beta * gamma * gamma * Wmax / (I * I)) - beta * beta - density_effect));

	}
	//printf("iron:%.2lf - %.2lf --> ", E_ini * 1000, dEdx * density* thickness * sqrt(angle*angle + 1));
	E_ini -= dEdx * density * thickness * sqrt(angle * angle + 1) / 1000;
	//printf("%.2lf \n", E_ini * 1000);
	// E_ini -= dEdx * thickness * sqrt(angle*angle + 1) / 1000;
	pbeta = (E_ini * E_ini - mass * mass) / E_ini;
	pbeta_mev = pbeta * 1000;
	if (pbeta_mev < pbmin)	pbeta_mev = pbmin;


}
void Calc_Energyloss_emulsion(double& pbeta_mev, double angle, double* par, double thickness, double pbmin) {
	double density = 2;
	double pbeta = pbeta_mev / 1000;
	double mass = par[0];
	double z = par[1];
	double E_ini = (pbeta + sqrt(pbeta * pbeta + 4.0 * mass * mass)) / 2.0;   // Energy [GeV]

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
		// information of constants //
		double me = m_e;                   // mass of electron [GeV]
		double mec2 = 0.510998928 / 1000.0;  // electron mass * c^2 [GeV]
		double K = 0.307075;               // 4πNA re^2 me c^2 [MeV mol^-1 cm^2]

		// information of absorber (Fe) //
		double Z = mtrl_Z[Absorber];              // atomic number of absorber           (PDG2017)
		double A = mtrl_A[Absorber];              // atomic mass of absorber [g mol^-1]  (PDG2017)
		double I = mtrl_I[Absorber] / 1000000000.0; // Mean excitation energy of Fe [GeV]  (PDG2017)
		//double density = mtrl_Density[Absorber];  // density of Fe [g cm^-3]             (PDG2017)


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
	//printf("emulsion:%.2lf - %.2lf --> ", E_ini * 1000, dEdx* density * thickness * sqrt(angle*angle + 1));
	E_ini -= dEdx * density * thickness * sqrt(angle * angle + 1) / 1000;
	//printf("%.2lf \n", E_ini * 1000);
	//E_ini -= dEdx * thickness * sqrt(angle*angle + 1) / 1000;
	pbeta = (E_ini * E_ini - mass * mass) / E_ini;
	pbeta_mev = pbeta * 1000;
	if (pbeta_mev < pbmin)	pbeta_mev = pbmin;


}

void Calc_VPH_trend(std::vector<output_inf>& out, int groupid, Momentum_recon::Mom_chain& c, VPH_dEdx_Param& par) {
	const double PB_MIN = 20;
	//proton
	bool flg = false;
	double pb_ini = c.Get_proton_mcs_pb(), angle;
	double evaluation_value = 0, vph_mean_theoretical, vph_sigma_theoretical;
	int pl, num = 0;
	output_inf out_tmp;
	out_tmp.groupid = groupid;
	for (auto& b : c.base) {
		//if (!flg && (c_inf.vph[pl][0] == 0 || c_inf.vph[pl][1] == 0))continue;
		if (!flg && (b.m[0].ph % 10000 < 20 && b.m[1].ph % 10000 < 20))continue;
		flg = true;
		pl = b.pl;
		angle = sqrt(b.ax * b.ax + b.ay * b.ay);
		out_tmp.angle = angle;
		out_tmp.chainid = c.chainid;
		out_tmp.mom = pb_ini;
		out_tmp.pl = pl;
		out_tmp.num = num;
		num++;
		out_tmp.vph = b.m[0].ph % 10000;
		out.push_back(out_tmp);
		//energy loss
	   // printf("::::::film::::::::::::\n");
		Calc_Energyloss_emulsion(pb_ini, angle, par.vph_mean_param, 0.007, PB_MIN);
		Calc_Energyloss_water(pb_ini, angle, par.vph_mean_param, 0.021, PB_MIN);
		Calc_Energyloss_emulsion(pb_ini, angle, par.vph_mean_param, 0.007, PB_MIN);
		//printf(":::::::::::::::::::\n");
		out_tmp.mom = pb_ini;
		out_tmp.num = num;
		num++;
		out_tmp.vph = b.m[1].ph % 10000;
		out.push_back(out_tmp);

		if (pl == 3 || pl == 15)continue;
		if (pl < 15 || pl % 2 == 0) {
			Calc_Energyloss_iron(pb_ini, angle, par.vph_mean_param, 0.05, PB_MIN);
		}
		else {
			Calc_Energyloss_water(pb_ini, angle, par.vph_mean_param, 0.25, PB_MIN);
		}

	}


}

void Calc_VPH_trend_inv(std::vector<output_inf>& out, int groupid, Momentum_recon::Mom_chain& c, VPH_dEdx_Param& par) {
	const double PB_MIN = 20;
	//proton
	bool flg = false;
	double pb_ini = c.Get_proton_mcs_pb(), angle;
	double evaluation_value = 0, vph_mean_theoretical, vph_sigma_theoretical;
	int pl, num = 0;
	output_inf out_tmp;
	out_tmp.groupid = groupid;
	for (auto itr = c.base.rbegin(); itr != c.base.rend(); itr++) {
		//if (!flg && (c_inf.vph[pl][0] == 0 || c_inf.vph[pl][1] == 0))continue;
		if (!flg && (itr->m[0].ph % 10000 < 20 && itr->m[1].ph % 10000 < 20))continue;
		flg = true;
		pl = itr->pl;
		angle = sqrt(itr->ax * itr->ax + itr->ay * itr->ay);
		out_tmp.angle = angle;
		out_tmp.chainid = c.chainid;
		out_tmp.mom = pb_ini;
		out_tmp.pl = pl;
		out_tmp.num = num;
		num++;
		out_tmp.vph = itr->m[1].ph % 10000;
		out.push_back(out_tmp);
		//energy loss
		//printf("::::::film::::::::::::\n");
		Calc_Energyloss_emulsion(pb_ini, angle, par.vph_mean_param, 0.007, PB_MIN);
		Calc_Energyloss_water(pb_ini, angle, par.vph_mean_param, 0.021, PB_MIN);
		Calc_Energyloss_emulsion(pb_ini, angle, par.vph_mean_param, 0.007, PB_MIN);
		//printf(":::::::::::::::::::::\n");
		out_tmp.mom = pb_ini;
		out_tmp.num = num;
		num++;
		out_tmp.vph = itr->m[0].ph % 10000;
		out.push_back(out_tmp);

		if (pl == 4 || pl == 16)continue;
		if (pl < 15 || pl % 2 == 1) {
			Calc_Energyloss_iron(pb_ini, angle, par.vph_mean_param, 0.05, PB_MIN);
		}
		else {
			Calc_Energyloss_water(pb_ini, angle, par.vph_mean_param, 0.25, PB_MIN);
		}

	}


}


bool judge_direction(Momentum_recon::Mom_chain& c) {
	double mean[2] = {};
	int count[2] = {};
	for (int i = 0; i < c.base.size() / 2; i++) {
		mean[0] += c.base[i].m[0].ph % 10000;
		count[0]++;
		mean[0] += c.base[i].m[1].ph % 10000;
		count[0]++;

		mean[1] += c.base[c.base.size() - i - 1].m[0].ph % 10000;
		count[1]++;
		mean[1] += c.base[c.base.size() - i - 1].m[1].ph % 10000;
		count[1]++;

	}
	mean[0] /= count[0];
	mean[1] /= count[1];
	return mean[0] < mean[1];

}
void Output_file(std::string filename, std::vector<output_inf>& out) {
	std::ofstream ofs(filename);
	int count = 0;
	for (auto o : out) {
		if (count % 10000 == 0) {
			fprintf(stderr, "\r Write VPH ... %d/%d (%4.1lf%%)", count, int(out.size()), count * 100. / out.size());
		}
		count++;

		ofs << std::right << std::fixed
			<< std::setw(10) << std::setprecision(0) << o.groupid << " "
			<< std::setw(10) << std::setprecision(0) << o.chainid << " "
			<< std::setw(4) << std::setprecision(0) << o.num << " "
			<< std::setw(4) << std::setprecision(0) << o.pl << " "
			<< std::setw(8) << std::setprecision(4) << o.angle << " "
			<< std::setw(8) << std::setprecision(2) << o.mom << " "
			<< std::setw(6) << std::setprecision(0) << o.vph << std::endl;
	}
	fprintf(stderr, "\r Write VPH ... %d/%d (%4.1lf%%)\n", count, int(out.size()), count * 100. / out.size());

}