#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

class output_format {
public:
	int groupid, chainid, pl;
	double sigma[2], chi2[4];
};

void Calc_chain_angle(Momentum_recon::Mom_chain& c, double& ax, double& ay);

void Calc_divide_vph(Momentum_recon::Mom_chain& c, int pl, double mean[2], double mean_err[2], double sd[2], int count[2]);
void Calc_divide_pixel(Momentum_recon::Mom_chain& c, int pl, double mean[2], double mean_err[2], double sd[2], int count[2]);
//void output_file(std::string filename, std::vector<Momentum_recon::Event_information>&mom);
void Calc_divide_angle_lateral(Momentum_recon::Mom_chain& c, int pl, double mean[2], double mean_err[2], double sd[2], int count[2], double ax, double ay);
void Calc_divide_angle_radial(Momentum_recon::Mom_chain& c, int pl, double mean[2], double mean_err[2], double sd[2], int count[2], double ax, double ay);
void output_file(std::string filename, std::vector<Momentum_recon::Event_information>& mom, std::vector<output_format>& out);
void output_file_best(std::string filename, std::vector<output_format>& out);
int main(int argc, char** argv) {
	if (argc != 4) {
		fprintf(stderr, "usage:file-in-momch file-out-momch.txt file-best_pl-momch.txt\n");
		exit(1);
	}
	std::string file_in_momch = argv[1];
	std::string file_out_momch = argv[2];
	std::string file_out_best_pl_momch = argv[3];

	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch);
	//std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_mom_chain_extension(file_in_momch);
	std::vector<output_format> out;
	output_file(file_out_momch, momch, out);
	output_file_best(file_out_best_pl_momch, out);


}


void output_file(std::string filename, std::vector<Momentum_recon::Event_information>& mom, std::vector<output_format>& out) {
	std::ofstream ofs(filename);
	output_format out_tmp;
	double mean[2], mean_err[2], sd[2], dal, rms_l, dar, rms_r, chi2[4], error, sigma;
	int count[2];
	double ax, ay;
	for (auto& ev: mom) {
		for (auto& c : ev.chains) {//add

			Calc_chain_angle(c, ax, ay);

			std::set<int> pl_set;
			for (auto itr = c.base.begin(); itr != c.base.end(); itr++) {
				pl_set.insert(itr->pl);
			}
			int pl_min = *pl_set.begin();
			int pl_max = *pl_set.rbegin();
			std::vector<int> pl_vec(pl_set.begin(), pl_set.end());
			//printf("PL range %d %d\n", pl_min, pl_max);
			for (int i = 0; i < pl_vec.size() - 1; i++) {
				int pl = pl_vec[i];
				out_tmp.groupid = ev.groupid;
				out_tmp.chainid = c.chainid;
				out_tmp.pl = pl;
				ofs << std::fixed << std::right
					<< std::setw(10) << std::setprecision(0) << ev.groupid << " "
					<< std::setw(10) << std::setprecision(0) << c.chainid << " "
					<< std::setw(7) << std::setprecision(4) << ax << " "
					<< std::setw(7) << std::setprecision(4) << ay << " "
					<< std::setw(3) << std::setprecision(0) << pl << " ";

				//lateralpx
				Calc_divide_angle_lateral(c, pl, mean, mean_err, sd, count, ax, ay);
				if (mean_err[0] > 0 && mean_err[1] > 0) {
					error = sqrt(pow(mean_err[1], 2) + pow(mean_err[0], 2));
					sigma = (mean[1] - mean[0]) / error;
					chi2[0] = sigma * sigma;
				}
				else {
					sigma = -1;
					chi2[0] = -1;
				}
				ofs << std::fixed << std::right
					<< std::setw(4) << std::setprecision(0) << count[0] << " "
					<< std::setw(8) << std::setprecision(4) << mean[0] << " "
					<< std::setw(6) << std::setprecision(5) << mean_err[0] << " "
					<< std::setw(4) << std::setprecision(0) << count[1] << " "
					<< std::setw(8) << std::setprecision(4) << mean[1] << " "
					<< std::setw(6) << std::setprecision(5) << mean_err[1] << " "
					<< std::setw(8) << std::setprecision(3) << sigma << " ";

				//radialpx
				Calc_divide_angle_radial(c, pl, mean, mean_err, sd, count, ax, ay);
				if (mean_err[0] > 0 && mean_err[1] > 0) {
					error = sqrt(pow(mean_err[1], 2) + pow(mean_err[0], 2));
					sigma = (mean[1] - mean[0]) / error;
					chi2[1] = sigma * sigma;
				}
				else {
					sigma = -1;
					chi2[1] = -1;
				}
				ofs << std::fixed << std::right
					<< std::setw(4) << std::setprecision(0) << count[0] << " "
					<< std::setw(8) << std::setprecision(4) << mean[0] << " "
					<< std::setw(6) << std::setprecision(5) << mean_err[0] << " "
					<< std::setw(4) << std::setprecision(0) << count[1] << " "
					<< std::setw(8) << std::setprecision(4) << mean[1] << " "
					<< std::setw(6) << std::setprecision(5) << mean_err[1] << " "
					<< std::setw(8) << std::setprecision(3) << sigma << " ";

				//vph
				Calc_divide_vph(c, pl, mean, mean_err, sd, count);
				if (mean_err[0] > 0 && mean_err[1] > 0) {
					//error = sqrt(pow(mean_err[1] / mean[0], 2) + pow(mean[1] / (mean[0] * mean[0])*mean_err[0], 2));
					//sigma = (mean[1] / mean[0] - 1) / error;
					//chi2[2] = sigma * sigma;
					error = sqrt(pow(mean_err[1], 2) + pow(mean_err[0], 2));
					sigma = (mean[1] - mean[0]) / error;
					chi2[2] = sigma * sigma;

				}
				else {
					sigma = -1;
					chi2[2] = -1;
				}
				ofs << std::fixed << std::right
					<< std::setw(4) << std::setprecision(0) << count[0] << " "
					<< std::setw(8) << std::setprecision(2) << mean[0] << " "
					<< std::setw(6) << std::setprecision(3) << mean_err[0] << " "
					<< std::setw(4) << std::setprecision(0) << count[1] << " "
					<< std::setw(8) << std::setprecision(2) << mean[1] << " "
					<< std::setw(6) << std::setprecision(3) << mean_err[1] << " "
					<< std::setw(8) << std::setprecision(3) << sigma << " ";
				out_tmp.sigma[0] = sigma;

				//pixel
				Calc_divide_pixel(c, pl, mean, mean_err, sd, count);
				if (mean_err[0] > 0 && mean_err[1] > 0) {
					error = sqrt(pow(mean_err[1] / mean[0], 2) + pow(mean[1] / (mean[0] * mean[0]) * mean_err[0], 2));
					sigma = (mean[1] / mean[0] - 1) / error;
					chi2[3] = sigma * sigma;
				}
				else {
					sigma = -1;
					chi2[3] = -1;
				}
				//ofs << std::fixed << std::right
				//	<< std::setw(4) << std::setprecision(0) << count[0] << " "
				//	<< std::setw(8) << std::setprecision(2) << mean[0] << " "
				//	<< std::setw(6) << std::setprecision(3) << mean_err[0] << " "
				//	<< std::setw(4) << std::setprecision(0) << count[1] << " "
				//	<< std::setw(8) << std::setprecision(2) << mean[1] << " "
				//	<< std::setw(6) << std::setprecision(3) << mean_err[1] << " "
				//	<< std::setw(8) << std::setprecision(3) << sigma << " ";
				out_tmp.sigma[1] = sigma;

				ofs << std::fixed << std::right
					<< std::setw(10) << std::setprecision(4) << chi2[0] << " "
					<< std::setw(10) << std::setprecision(4) << chi2[1] << " "
					<< std::setw(10) << std::setprecision(4) << chi2[2] << std::endl;
				//<< std::setw(10) << std::setprecision(4) << chi2[3] << std::endl;
				out_tmp.chi2[0] = chi2[0];
				out_tmp.chi2[1] = chi2[1];
				out_tmp.chi2[2] = chi2[2];
				out_tmp.chi2[3] = chi2[3];
				if (pl < pl_min + 5)continue;
				if (pl >= pl_max - 5)continue;
				out.push_back(out_tmp);
			}
		}

	}

	ofs.close();
}
void output_file_best(std::string filename, std::vector<output_format>& out) {
	std::multimap<int, output_format> out_map;
	for (auto itr = out.begin(); itr != out.end(); itr++) {
		//if (itr->sigma[0] < 0)continue;
		//if (itr->sigma[1] < 0)continue;
		out_map.insert(std::make_pair(itr->groupid, *itr));
	}
	std::ofstream ofs(filename);

	for (auto itr = out_map.begin(); itr != out_map.end(); itr = std::next(itr, out_map.count(itr->first))) {
		std::vector<output_format> out_v;
		auto range = out_map.equal_range(itr->first);
		for (auto res = range.first; res != range.second; res++) {
			out_v.push_back(res->second);
		}


		double max_chi2 = 0;
		output_format out_pl;
		for (auto itr2 = out_v.begin(); itr2 != out_v.end(); itr2++) {
			if (max_chi2 < itr2->chi2[0] + itr2->chi2[1] + itr2->chi2[2]) {
				max_chi2 = itr2->chi2[0] + itr2->chi2[1] + itr2->chi2[2];
				out_pl = *itr2;
			}
		}
		ofs << std::fixed << std::right
			<< std::setw(10) << std::setprecision(0) << out_pl.groupid << " "
			<< std::setw(10) << std::setprecision(0) << out_pl.chainid << " "
			<< std::setw(10) << std::setprecision(0) << out_pl.pl << " "
			<< std::setw(8) << std::setprecision(3) << out_pl.sigma[0] << " "
			<< std::setw(8) << std::setprecision(3) << out_pl.sigma[1] << " "
			<< std::setw(10) << std::setprecision(4) << out_pl.chi2[0] << " "
			<< std::setw(10) << std::setprecision(4) << out_pl.chi2[1] << " "
			<< std::setw(10) << std::setprecision(4) << out_pl.chi2[2] << " "
			<< std::setw(10) << std::setprecision(4) << out_pl.chi2[3] << std::endl;

	}

}
void Calc_chain_angle(Momentum_recon::Mom_chain& c, double& ax, double& ay) {
	ax = 0;
	ay = 0;
	int count = 0;
	for (auto itr = c.base.begin(); itr != c.base.end(); itr++) {
		ax += itr->ax;
		ay += itr->ay;
		count++;
	}
	ax = ax / count;
	ay = ay / count;


}
void Calc_divide_vph(Momentum_recon::Mom_chain& c, int pl, double mean[2], double mean_err[2], double sd[2], int count[2]) {
	count[0] = 0;
	count[1] = 0;
	double sum2[2] = {}, sum1[2] = {};

	for (auto itr = c.base.begin(); itr != c.base.end(); itr++) {
		if (itr->pl <= pl) {
			for (int i = 0; i < 2; i++) {
				count[0] += 1;
				sum1[0] += itr->m[i].ph % 10000;
				sum2[0] += pow(itr->m[i].ph % 10000, 2);
			}
		}
		else {
			for (int i = 0; i < 2; i++) {
				count[1] += 1;
				sum1[1] += itr->m[i].ph % 10000;
				sum2[1] += pow(itr->m[i].ph % 10000, 2);
			}
		}
	}

	for (int i = 0; i < 2; i++) {

		if (count[i] < 2) {
			mean[i] = -1;
			mean_err[i] = -1;
			sd[i] = -1;
		}
		else {
			mean[i] = sum1[i] / count[i];
			sd[i] = sum2[i] / count[i] - mean[i] * mean[i];
			sd[i] = std::max(sd[i], 1.);
			sd[i] = sqrt(sd[i]) * sqrt(count[i] - 1) / sqrt(count[i]);
			//vph mip rms
			sd[i] = std::max(sd[i], 10.);
			mean_err[i] = sd[i] / sqrt(count[i]);
		}
	}
}
void Calc_divide_pixel(Momentum_recon::Mom_chain& c, int pl, double mean[2], double mean_err[2], double sd[2], int count[2]) {
	count[0] = 0;
	count[1] = 0;
	double sum2[2] = {}, sum1[2] = {};

	for (auto itr = c.base.begin(); itr != c.base.end(); itr++) {
		if (itr->pl <= pl) {
			for (int i = 0; i < 2; i++) {
				if (itr->m[i].hitnum < 0)continue;
				count[0] += 1;
				sum1[0] += itr->m[i].hitnum;
				sum2[0] += pow(itr->m[i].hitnum, 2);
			}
		}
		else {
			for (int i = 0; i < 2; i++) {
				if (itr->m[i].hitnum < 0)continue;
				count[1] += 1;
				sum1[1] += itr->m[i].hitnum;
				sum2[1] += pow(itr->m[i].hitnum, 2);
			}
		}
	}

	for (int i = 0; i < 2; i++) {

		if (count[i] < 2) {
			mean[i] = -1;
			mean_err[i] = -1;
			sd[i] = -1;
		}
		else {
			mean[i] = sum1[i] / count[i];
			sd[i] = sum2[i] / count[i] - mean[i] * mean[i];
			if (sd[i] <= 0) {
				mean[i] = -1;
				mean_err[i] = -1;
				sd[i] = -1;
			}
			else {
				sd[i] = sqrt(sd[i]) * sqrt(count[i] - 1) / sqrt(count[i]);
				mean_err[i] = sd[i] / sqrt(count[i]);
			}
		}
	}
}
void Calc_divide_angle_lateral(Momentum_recon::Mom_chain& c, int pl, double mean[2], double mean_err[2], double sd[2], int count[2], double ax, double ay) {
	double angle = sqrt(ax * ax + ay * ay);
	count[0] = 0;
	count[1] = 0;
	double sum2[2] = {}, sum1[2] = {};

	double  dax[2], day[2], dang[2], diff;
	for (auto itr = c.base_pair.begin(); itr != c.base_pair.end(); itr++) {
		dax[0] = itr->first.ax - ax;
		day[0] = itr->first.ay - ay;
		dax[1] = itr->second.ax - ax;
		day[1] = itr->second.ay - ay;
		if (angle < 0.01) {
			dang[0] = dax[0];
			dang[1] = dax[1];
		}
		else {
			dang[0] = (dax[0] * ay - day[0] * ax) / angle;
			dang[1] = (dax[1] * ay - day[1] * ax) / angle;
		}
		diff = dang[1] - dang[0];
		if (itr->first.pl == pl)continue;
		else if (itr->first.pl < pl) {
			count[0] += 1;
			sum1[0] += diff;
			sum2[0] += pow(diff, 2);
		}
		else {
			count[1] += 1;
			sum1[1] += diff;
			sum2[1] += pow(diff, 2);
		}
	}

	for (int i = 0; i < 2; i++) {

		if (count[i] < 2) {
			mean[i] = -1;
			mean_err[i] = -1;
			sd[i] = -1;
		}
		else {
			mean[i] = sum1[i] / count[i];
			sd[i] = sum2[i] / count[i] - mean[i] * mean[i];
			sd[i] = std::max(sd[i], 0.001 * 0.001);
			sd[i] = sqrt(sd[i]) * sqrt(count[i] - 1) / sqrt(count[i]);
			//pxx
			sd[i] = std::max(sd[i], 0.002);
			mean_err[i] = sd[i] / sqrt(count[i]);
		}
	}
}
void Calc_divide_angle_radial(Momentum_recon::Mom_chain& c, int pl, double mean[2], double mean_err[2], double sd[2], int count[2], double ax, double ay) {
	double angle = sqrt(ax * ax + ay * ay);
	count[0] = 0;
	count[1] = 0;
	double sum2[2] = {}, sum1[2] = {};

	double  dax[2], day[2], dang[2], diff;
	for (auto itr = c.base_pair.begin(); itr != c.base_pair.end(); itr++) {
		dax[0] = itr->first.ax - ax;
		day[0] = itr->first.ay - ay;
		dax[1] = itr->second.ax - ax;
		day[1] = itr->second.ay - ay;
		if (angle < 0.01) {
			dang[0] = day[0];
			dang[1] = day[1];
		}
		else {
			dang[0] = (dax[0] * ax + day[0] * ay) / angle;
			dang[1] = (dax[1] * ax + day[1] * ay) / angle;
		}
		diff = dang[1] - dang[0];
		if (itr->first.pl <= pl) {
			count[0] += 1;
			sum1[0] += diff;
			sum2[0] += pow(diff, 2);
		}
		else {
			count[1] += 1;
			sum1[1] += diff;
			sum2[1] += pow(diff, 2);
		}
	}

	for (int i = 0; i < 2; i++) {

		if (count[i] < 2) {
			mean[i] = -1;
			mean_err[i] = -1;
			sd[i] = -1;
		}
		else {
			mean[i] = sum1[i] / count[i];
			sd[i] = sum2[i] / count[i] - mean[i] * mean[i];
			mean[i] = sum1[i] / count[i];
			sd[i] = sum2[i] / count[i] - mean[i] * mean[i];
			sd[i] = std::max(sd[i], 0.001 * 0.001);
			sd[i] = sqrt(sd[i]) * sqrt(count[i] - 1) / sqrt(count[i]);
			//pxx
			sd[i] = std::max(sd[i], 0.015 * angle + 0.002);
			mean_err[i] = sd[i] / sqrt(count[i]);
		}
	}
}
