#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

class output_format {
public:
	int groupid, chainid, pl, dpl;
	double sigma[2], chi2[4], md, dz;
};

void Calc_divide_vph(Momentum_recon::Mom_chain& c, int pl, double& mean, double& sd, int& count, double& vph);
void Calc_divide_pixel(Momentum_recon::Mom_chain& c, int pl, double& mean, double& sd, int& count, double& hitnum);
void output_file(std::string filename, std::vector<Momentum_recon::Event_information>& mom, std::vector<output_format>& out);
void Calc_divide_angle_lateral(Momentum_recon::Mom_chain& c, int pl, double& dal, double& rms);
void Calc_divide_angle_radial(Momentum_recon::Mom_chain& c, int pl, double& dar, double& rms);
void output_file_best(std::string filename, std::vector<output_format>& out);
double Calc_md(Momentum_recon::Mom_chain& c, int pl, double& dz);
void cut_upstream_base(Momentum_recon::Mom_chain& c, int pl);
void cut_upstream(std::vector<Momentum_recon::Event_information>& mom);

int main(int argc, char** argv) {
	if (argc != 5) {
		fprintf(stderr, "usage:file-in-momch file-out-momch out_momch\n");
		exit(1);
	}
	std::string file_in_momch = argv[1];
	std::string file_out_momch = argv[2];
	std::string file_out_best_pl_momch = argv[3];
	std::string file_cut_momch = argv[4];

	//std::vector<Momentum_recon::Mom_chain> momch = Momentum_recon::Read_mom_chain_extension(file_in_momch);
	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch);

	std::vector<output_format>out;
	output_file(file_out_momch, momch, out);
	output_file_best(file_out_best_pl_momch, out);

	cut_upstream(momch);//chi2Ç™thrÇÊÇËçÇÇ¢Ç‡ÇÃÇÕçÌèú
	//Momentum_recon::Write_mom_chain_extension(file_cut_momch, momch);
	Momentum_recon::Write_Event_information_extension(file_cut_momch, momch);
}
void output_file(std::string filename, std::vector<Momentum_recon::Event_information>& mom, std::vector<output_format>& out) {
	std::ofstream ofs(filename);

	double mean[2], mean_err[2], sd[2], dal, rms_l, dar, rms_r, chi2[4], error, sigma, vph, hitnum;
	int count[2];
	for (auto& ev : mom) {

		for (auto& c : ev.chains) {
			std::set<int> pl_set;
			for (auto itr = c.base.begin(); itr != c.base.end(); itr++) {
				pl_set.insert(itr->pl);
			}
			std::vector<int> pl_v(pl_set.begin(), pl_set.end());
			std::sort(pl_v.begin(), pl_v.end(), std::greater<int>());
			int pl;
			for (int i = 1; i < pl_v.size() && i <= 2; i++) {
				output_format out_tmp;
				out_tmp.groupid = ev.groupid;
				out_tmp.chainid = c.chainid;
				out_tmp.pl = pl_v[i];
				out_tmp.dpl = pl_v[i - 1] - pl_v[i];

				pl = pl_v[i];
				Calc_divide_angle_lateral(c, pl, dal, rms_l);
				Calc_divide_angle_radial(c, pl, dar, rms_r);
				if (!isfinite(dal) || !isfinite(dar) || !isfinite(rms_l) || !isfinite(rms_r))continue;
				chi2[0] = pow(dal / rms_l, 2);
				chi2[1] = pow(dar / rms_r, 2);

				ofs << std::fixed << std::right
					<< std::setw(10) << std::setprecision(0) << ev.groupid << " "
					<< std::setw(10) << std::setprecision(0) << c.chainid << " "
					<< std::setw(10) << std::setprecision(0) << pl << " ";

				ofs << std::fixed << std::right
					<< std::setw(8) << std::setprecision(5) << dal << " "
					<< std::setw(8) << std::setprecision(5) << rms_l << " "
					<< std::setw(8) << std::setprecision(5) << dar << " "
					<< std::setw(8) << std::setprecision(5) << rms_r << " ";


				Calc_divide_vph(c, pl, mean[0], sd[0], count[0], vph);
				Calc_divide_pixel(c, pl, mean[1], sd[1], count[1], hitnum);
				if (sd[0] > 0 && vph > 0) {
					sigma = (vph - mean[0]) / sd[0];
					chi2[2] = sigma * sigma;
				}
				else {
					sigma = -1;
					chi2[2] = -1;
				}
				out_tmp.sigma[0] = sigma;

				ofs << std::fixed << std::right
					<< std::setw(4) << std::setprecision(0) << count[0] << " "
					<< std::setw(8) << std::setprecision(2) << mean[0] << " "
					<< std::setw(6) << std::setprecision(3) << sd[0] << " "
					<< std::setw(8) << std::setprecision(2) << vph << " "
					<< std::setw(6) << std::setprecision(2) << sigma << " ";

				if (sd[1] > 0 && hitnum > 0) {
					sigma = (hitnum - mean[1]) / sd[1];
					chi2[3] = sigma * sigma;
				}
				else {
					sigma = -1;
					chi2[3] = -1;
				}
				out_tmp.sigma[1] = sigma;

				ofs << std::fixed << std::right
					<< std::setw(4) << std::setprecision(0) << count[1] << " "
					<< std::setw(8) << std::setprecision(2) << mean[1] << " "
					<< std::setw(6) << std::setprecision(3) << sd[1] << " "
					<< std::setw(8) << std::setprecision(2) << hitnum << " "
					<< std::setw(6) << std::setprecision(2) << sigma << " ";
				out_tmp.md = Calc_md(c, pl, out_tmp.dz);

				ofs << std::fixed << std::right
					<< std::setw(10) << std::setprecision(4) << chi2[0] << " "
					<< std::setw(10) << std::setprecision(4) << chi2[1] << " "
					<< std::setw(10) << std::setprecision(4) << chi2[2] << " "
					<< std::setw(10) << std::setprecision(4) << chi2[3] << " "
					<< std::setw(10) << std::setprecision(4) << out_tmp.md << " "
					<< std::setw(10) << std::setprecision(4) << out_tmp.dz << std::endl;

				out_tmp.chi2[0] = chi2[0];
				out_tmp.chi2[1] = chi2[1];
				out_tmp.chi2[2] = chi2[2];
				out_tmp.chi2[3] = chi2[3];
				//if (out_tmp.sigma[0] < 0)continue;
				//if (out_tmp.sigma[1] < 0)continue;
				out.push_back(out_tmp);

			}
		}

	}
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
			<< std::setw(10) << std::setprecision(0) << out_pl.dpl << " "
			<< std::setw(8) << std::setprecision(3) << out_pl.sigma[0] << " "
			<< std::setw(8) << std::setprecision(3) << out_pl.sigma[1] << " "
			<< std::setw(10) << std::setprecision(4) << out_pl.chi2[0] << " "
			<< std::setw(10) << std::setprecision(4) << out_pl.chi2[1] << " "
			<< std::setw(10) << std::setprecision(4) << out_pl.chi2[2] << " "
			<< std::setw(10) << std::setprecision(4) << out_pl.chi2[3] << " "
			<< std::setw(10) << std::setprecision(4) << out_pl.md << " "
			<< std::setw(10) << std::setprecision(4) << out_pl.dz << std::endl;

	}

}



void Calc_divide_vph(Momentum_recon::Mom_chain& c, int pl, double& mean, double& sd, int& count, double& vph) {
	vph = -1;
	count = 0;
	double sum2 = 0, sum1 = 0;
	bool flg = true;
	for (auto itr = c.base.begin(); itr != c.base.end(); itr++) {
		if (itr->pl <= pl) {
			for (int i = 0; i < 2; i++) {
				count += 1;
				sum1 += itr->m[i].ph % 10000;
				sum2 += pow(itr->m[i].ph % 10000, 2);
			}
		}
		else {
			if (flg) {
				vph = (itr->m[0].ph % 10000 + itr->m[1].ph % 10000) / 2;
				flg = false;

			}
		}
	}
	if (count < 2) {
		mean = -1;
		sd = -1;
	}
	else {
		mean = sum1 / count;
		sd = sum2 / count - mean * mean;
		if (sd <= 0) {
			mean = -1;
			sd = -1;
		}
		else {
			sd = sqrt(sd) * sqrt(count - 1) / sqrt(count);
		}

	}
}
void Calc_divide_pixel(Momentum_recon::Mom_chain& c, int pl, double& mean, double& sd, int& count, double& hitnum) {
	hitnum = -1;
	count = 0;
	double sum2 = 0., sum1 = 0;
	bool flg = true;
	for (auto itr = c.base.begin(); itr != c.base.end(); itr++) {
		if (itr->pl <= pl) {
			for (int i = 0; i < 2; i++) {
				if (itr->m[i].hitnum < 0)continue;
				count += 1;
				sum1 += itr->m[i].hitnum;
				sum2 += pow(itr->m[i].hitnum, 2);
			}
		}
		else {
			if (flg) {
				if (itr->m[0].hitnum > 0 && itr->m[1].hitnum > 0) {
					hitnum = (itr->m[0].hitnum + itr->m[1].hitnum) / 2;
				}
				else if (itr->m[0].hitnum > 0) {
					hitnum = itr->m[0].hitnum;
				}
				else if (itr->m[1].hitnum > 0) {
					hitnum = itr->m[1].hitnum;
				}
				else {
					hitnum = -1;
				}
				flg = false;
			}
		}
	}
	if (count < 2) {
		mean = -1;
		sd = -1;
	}
	else {
		mean = sum1 / count;
		sd = sum2 / count - mean * mean;
		if (sd <= 0) {
			mean = -1;
			sd = -1;
		}
		else {
			sd = sqrt(sd) * sqrt(count - 1) / sqrt(count);
		}

	}
}
void Calc_divide_angle_lateral(Momentum_recon::Mom_chain& c, int pl, double& dal, double& rms) {
	dal = NAN;

	double ax, ay, dax, day, dang, angle, sum2 = 0;
	int count = 0;
	for (auto itr = c.base_pair.begin(); itr != c.base_pair.end(); itr++) {
		ax = itr->first.ax;
		ay = itr->first.ay;
		angle = sqrt(ax * ax + ay * ay);
		dax = itr->second.ax - ax;
		day = itr->second.ay - ay;
		if (angle < 0.01) {
			dang = dax;
		}
		else {
			dang = (dax * ay - day * ax) / angle;
		}
		if (itr->first.pl == pl) {
			dal = dang;
		}
		else if (itr->first.pl < pl) {
			sum2 += dang * dang;
			count++;
		}
	}
	if (count <= 1) {
		rms = NAN;
	}
	else {
		rms = sqrt(sum2 / count);
	}


}
void Calc_divide_angle_radial(Momentum_recon::Mom_chain& c, int pl, double& dar, double& rms) {
	dar = NAN;
	double ax, ay, dax, day, dang, angle, sum2 = 0;
	int count = 0;
	for (auto itr = c.base_pair.begin(); itr != c.base_pair.end(); itr++) {
		ax = itr->first.ax;
		ay = itr->first.ay;
		angle = sqrt(ax * ax + ay * ay);
		dax = itr->second.ax - ax;
		day = itr->second.ay - ay;
		if (angle < 0.01) {
			dang = day;
		}
		else {
			dang = (dax * ax + day * ay) / angle;
		}
		if (itr->first.pl == pl) {
			dar = dang;
		}
		else if (itr->first.pl < pl) {
			sum2 += dang * dang;
			count++;
		}
	}
	if (count <= 1) {
		rms = NAN;
	}
	else {
		rms = sqrt(sum2 / count);
	}

}
double Calc_md(Momentum_recon::Mom_chain& c, int pl, double& dz) {
	double md = -1;
	for (auto itr = c.base_pair.begin(); itr != c.base_pair.end(); itr++) {
		if (pl == itr->first.pl) {
			matrix_3D::vector_3D pos0, pos1, dir0, dir1;
			pos0.x = itr->first.x;
			pos0.y = itr->first.y;
			pos0.z = itr->first.z;
			dir0.x = itr->first.ax;
			dir0.y = itr->first.ay;
			dir0.z = 1;
			pos1.x = itr->second.x;
			pos1.y = itr->second.y;
			pos1.z = itr->second.z;
			dir1.x = itr->second.ax;
			dir1.y = itr->second.ay;
			dir1.z = 1;
			double extra[2], z_range[2];
			z_range[0] = pos0.z;
			z_range[1] = pos1.z;
			md = matrix_3D::minimum_distance(pos0, pos1, dir0, dir1, z_range, extra);
			dz = (fabs(extra[0]) + fabs(extra[1])) / 2;
			return md;

		}
	}
	return md;
}

void cut_upstream(std::vector<Momentum_recon::Event_information>& mom) {
	std::vector<Momentum_recon::Mom_chain> ret;
	double mean[2], mean_err[2], sd[2], dal, rms_l, dar, rms_r, chi2[4], error, sigma, vph, hitnum, md;
	int count[2];
	for (auto& ev : mom) {
		for (auto& c : ev.chains) {
			std::set<int> pl_set;
			for (auto itr = c.base.begin(); itr != c.base.end(); itr++) {
				pl_set.insert(itr->pl);
			}
			std::vector<int> pl_v(pl_set.begin(), pl_set.end());
			std::sort(pl_v.begin(), pl_v.end(), std::greater<int>());
			int pl;
			for (int i = 1; i < pl_v.size() && i <= 2; i++) {
				output_format out_tmp;
				out_tmp.groupid = ev.groupid;
				out_tmp.chainid = c.chainid;
				out_tmp.pl = pl_v[i];
				out_tmp.dpl = pl_v[i - 1] - pl_v[i];

				pl = pl_v[i];
				Calc_divide_angle_lateral(c, pl, dal, rms_l);
				Calc_divide_angle_radial(c, pl, dar, rms_r);
				if (!isfinite(dal) || !isfinite(dar) || !isfinite(rms_l) || !isfinite(rms_r))continue;
				chi2[0] = pow(dal / rms_l, 2);
				chi2[1] = pow(dar / rms_r, 2);

				Calc_divide_vph(c, pl, mean[0], sd[0], count[0], vph);
				Calc_divide_pixel(c, pl, mean[1], sd[1], count[1], hitnum);
				if (sd[0] > 0 && vph > 0) {
					sigma = (vph - mean[0]) / sd[0];
					chi2[2] = sigma * sigma;
				}
				else {
					sigma = -1;
					chi2[2] = -1;
				}
				out_tmp.sigma[0] = sigma;

				if (sd[1] > 0 && hitnum > 0) {
					sigma = (hitnum - mean[1]) / sd[1];
					chi2[3] = sigma * sigma;
				}
				else {
					sigma = -1;
					chi2[3] = -1;
				}
				out_tmp.sigma[1] = sigma;
				double dz;
				md = Calc_md(c, pl, dz);

				out_tmp.chi2[0] = chi2[0];
				out_tmp.chi2[1] = chi2[1];
				out_tmp.chi2[2] = chi2[2];
				out_tmp.chi2[3] = chi2[3];
				if (chi2[0] + chi2[1] > 100 || chi2[2] > 25 || md > 100) {
					cut_upstream_base(c, pl);
					printf("eventid %d\n", ev.groupid);

				}
			}
		}
	}
}
void cut_upstream_base(Momentum_recon::Mom_chain& c, int pl) {
	for (auto itr = c.base.begin(); itr != c.base.end();) {
		if (itr->pl > pl) {
			itr = c.base.erase(itr);
		}
		else {
			itr++;
		}
	}
	for (auto itr = c.base_pair.begin(); itr != c.base_pair.end();) {
		if (itr->first.pl >= pl) {
			itr = c.base_pair.erase(itr);
		}
		else {
			itr++;
		}
	}
}