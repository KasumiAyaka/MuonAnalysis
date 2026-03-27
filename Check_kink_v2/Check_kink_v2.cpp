// 2024/05/09
// kasumi
// based on "Check_upstream_base.cpp

#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>
#include <iomanip>

class output_format {
public:
	int groupid, chainid, pl, dpl;
	double sigma[2], chi2[4], md, dz;
};

struct kink_cand {
	int pl, groupid;
};
bool operator<(const kink_cand& lhs, const kink_cand& rhs) {
	return std::tie(lhs.groupid, lhs.pl) < std::tie(rhs.groupid, rhs.pl);
}



void Calc_divide_vph(Momentum_recon::Mom_chain& c, int pl, double& mean, double& sd, int& count, double& vph);
void Calc_divide_pixel(Momentum_recon::Mom_chain& c, int pl, double& mean, double& sd, int& count, double& hitnum);
void Kink_search(std::vector<Momentum_recon::Event_information>& mom, std::vector<output_format>& out, std::multimap<int, int>& kink, double thr_ang, double thr_vph, double thr_md);
void Calc_divide_angle_lateral(Momentum_recon::Mom_chain& c, int pl, double& dal, double& rms);
void Calc_divide_angle_radial(Momentum_recon::Mom_chain& c, int pl, double& dar, double& rms);
void output_file_best(std::string filename, std::vector<output_format>& out);
double Calc_md(Momentum_recon::Mom_chain& c, int pl, double& dz);
void cut_upstream_base(Momentum_recon::Mom_chain& c, int pl);
void cut_upstream(std::vector<Momentum_recon::Event_information>& mom, std::multimap<int, int>& kink, std::string output, std::string output2, std::string output3);

void Calc_divide_vph_inv(Momentum_recon::Mom_chain& c, int pl, double& mean, double& sd, int& count, double& vph);
void Calc_divide_pixel_inv(Momentum_recon::Mom_chain& c, int pl, double& mean, double& sd, int& count, double& hitnum);
void Calc_divide_angle_lateral_inv(Momentum_recon::Mom_chain& c, int pl, double& dal, double& rms);
void Calc_divide_angle_radial_inv(Momentum_recon::Mom_chain& c, int pl, double& dar, double& rms);
void up_down_ave(std::vector<Momentum_recon::Event_information>& mom, std::multimap<int, int>& kink, std::string output);

int main(int argc, char** argv) {
	if (argc != 5&&argc!=8) {
		fprintf(stderr, "usage:file-in-momch file-out-momch file-out-momch2\n");
		fprintf(stderr, "usage:file-in-momch.momch diff.txt kink_cand.txt kink_point.txt\n");
		exit(1);
	}
	std::string file_in_momch = argv[1];
	//std::string file_out_momch = argv[2];
	std::string file_out_best_pl_momch = argv[2];
	std::string out2 = argv[3];//"kink_cand.txt"
	std::string out3 = argv[4];//"KinkPoint.txt"
	double thr_ChiAng = 80;
	double thr_ChiVph = 20;
	double thr_MD = 80;

	if (argc == 8) {
		thr_ChiAng = std::stod(argv[5]);
		thr_ChiVph = std::stod(argv[6]);
		thr_MD = std::stod(argv[7]);
	}

	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch);

	std::vector<output_format>out;
	std::multimap<int, int> kink;
	Kink_search( momch, out, kink,thr_ChiAng,thr_ChiVph,thr_MD);
	//output_file_best(file_out_best_pl_momch, out);

	cut_upstream(momch, kink, file_out_best_pl_momch,out2,out3);
	//Momentum_recon::Write_Event_information_extension(file_cut_momch, momch);
}
void Kink_search(std::vector<Momentum_recon::Event_information>& mom, std::vector<output_format>& out, std::multimap<int, int>& kink, double thr_ang, double thr_vph, double thr_md) {

	double mean[2], mean_err[2], sd[2], dal, rms_l, dar, rms_r, chi2[4], error, sigma, vph, hitnum;
	int count[2];
	kink_cand kc;
	for (auto& ev : mom) {

		for (auto& c : ev.chains) {
			std::set<int> pl_set;
			for (auto itr = c.base.begin(); itr != c.base.end(); itr++) {
				pl_set.insert(itr->pl);
			}
			std::vector<int> pl_v(pl_set.begin(), pl_set.end());
			std::sort(pl_v.begin(), pl_v.end(), std::greater<int>());
			int pl;
			for (int i = 1; i < pl_v.size(); i++) {// && i <= 2; i++) {
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

				out_tmp.md = Calc_md(c, pl, out_tmp.dz);

				out_tmp.chi2[0] = chi2[0];
				out_tmp.chi2[1] = chi2[1];
				out_tmp.chi2[2] = chi2[2];
				out_tmp.chi2[3] = chi2[3];
				//if (out_tmp.sigma[0] < 0)continue;
				//if (out_tmp.sigma[1] < 0)continue;
				if (chi2[0] + chi2[1] > thr_ang || chi2[2] > thr_vph || out_tmp.md > thr_md) {//–ÚŽ‹check‚Ì‚½‚ßŠÉ‚ß
				//	if (chi2[0] + chi2[1] > 100 || chi2[2] > 25 || out_tmp.md > 100) {
					kc.groupid = ev.groupid;
					kc.pl = pl;
					kink.insert(std::make_pair(kc.groupid, kc.pl));
				}

				out.push_back(out_tmp);

			}
		}

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

void cut_upstream(std::vector<Momentum_recon::Event_information>& mom, std::multimap<int, int>& kink, std::string output, std::string output2, std::string output3) {

	std::ofstream ofs(output);
	std::ofstream ofs1(output2);
	std::ofstream ofs2(output3);
	std::vector<Momentum_recon::Mom_chain> ret;
	int flg = 0;
	int cnt1 = 0;
	int cnt2 = 0;
	double ax, ay, dax, day, dang, angle;
	double dalat, darad, dtan;
	int cnt;

	for (auto& ev : mom) {
		for (auto& c : ev.chains) {
			std::set<int> pl_set;
			int pl0;
			auto itr0 = kink.equal_range(ev.groupid);

			for (auto itrk = itr0.first; itrk != itr0.second; itrk++) {
				cnt = 0;
				double avetan[3] = { 0 };
				double avedrad[3] = { 0 };
				double avedlat[3] = { 0 };
				std::cout << itrk->first << " " << itrk->second << std::endl;//gid,pl

				for (auto itr = c.base.begin(); itr != c.base.end(); itr++) {
					pl0 = itrk->second;

					if (cnt != 0) {
						dax = itr->ax - ax;
						day = itr->ay - ay;
						dtan = sqrt(itr->ax * itr->ax + itr->ay * itr->ay) - angle;

						if (angle < 0.01) {
							dalat = dax;
							darad = day;
						}
						else {
							darad = (dax * ax + day * ay) / angle;
							dalat = (dax * ay - day * ax) / angle;
						}

						ofs << std::fixed << std::right
							<< std::setw(10) << std::setprecision(0) << ev.groupid << " "
							<< std::setw(4) << std::setprecision(0) << itr->pl << " ";

						if (itr->pl > pl0) {
							ofs << std::setw(3) << std::setprecision(0) << 1 << " ";
							avetan[0] = avetan[0] + dtan;
							avedrad[0] = avedrad[0] + darad;
							avedlat[0] = avedlat[0] +dalat;
							cnt1++;

						}
						else if (itr->pl == pl0) {
							ofs << std::setw(3) << std::setprecision(0) << 0 << " ";
							avetan[2] = dtan;
							avedrad[2] =darad;
							avedlat[2] =dalat;

						}
						else {
							ofs << std::setw(3) << std::setprecision(0) << -1 << " ";
							cnt2++;
							avetan[1] = avetan[1] + dtan;
							avedrad[1] = avedrad[1] + darad;
							avedlat[1] = avedlat[1] + dalat;

						}

						ofs << std::fixed << std::right
							<< std::setw(8) << std::setprecision(5) << dax << " "
							<< std::setw(8) << std::setprecision(5) << day << " "
							<< std::setw(8) << std::setprecision(5) << dtan << " "
							<< std::setw(8) << std::setprecision(5) << dalat << " "
							<< std::setw(8) << std::setprecision(5) << darad << " "
							<< std::endl;

					}
					ax = itr->ax;
					ay = itr->ay;
					angle = sqrt(ax * ax + ay * ay);

					cnt++;
				}
				avetan[0] = avetan[0]/cnt1;
				avedrad[0] = avedrad[0] / cnt1;
				avedlat[0] = avedlat[0] / cnt1;
				avetan[1] = avetan[1] / cnt2;
				avedrad[1] = avedrad[1] / cnt2;
				avedlat[1] = avedlat[1] / cnt2;


				ofs2 << std::fixed << std::right 
					<< std::setw(10) << std::setprecision(0) << ev.groupid << " "
					<< std::setw(4) << std::setprecision(0) << pl0 << " "					
					<< std::setw(8) << std::setprecision(5) << avetan[0] << " "
					<< std::setw(8) << std::setprecision(5) << avetan[2] << " "
					<< std::setw(8) << std::setprecision(5) << avetan[1] << " "
					<< std::setw(8) << std::setprecision(5) << avedrad[0] << " "
					<< std::setw(8) << std::setprecision(5) << avedrad[2] << " "
					<< std::setw(8) << std::setprecision(5) << avedrad[1] << " "
					<< std::setw(8) << std::setprecision(5) << avedlat[0] << " "
					<< std::setw(8) << std::setprecision(5) << avedlat[2] << " "
					<< std::setw(8) << std::setprecision(5) << avedlat[1] << " "
					<< std::endl;

				if (avedlat[2] > 0.02) {
					ofs1 << std::setw(10) << std::setprecision(0) << ev.groupid << " "
						<< std::setw(4) << std::setprecision(0) << pl0 << " "
						<< std::endl;
				}
				cnt1 = 0;
				cnt2 = 0;
			}

		}
	}
}
void up_down_ave(std::vector<Momentum_recon::Event_information>& mom, std::multimap<int, int>& kink, std::string output) {


	std::ofstream ofs(output);
	std::ofstream ofs1("kink_cand.txt");
	std::ofstream ofs2("tmp.txt");
	std::vector<Momentum_recon::Mom_chain> ret;
	int flg = 0;
	int cnt1 = 0;
	int cnt2 = 0;
	double ax, ay, dax, day, dang, angle;
	double dalat, darad, dtan;
	int cnt;

	double sum[2], ave[2], sig[2];
	double al[2], dal[2];
	int upl,dpl,tan;

	for (auto& ev : mom) {
		for (auto& c : ev.chains) {
			std::set<int> pl_set;
			int pl0;
			auto itr0 = kink.equal_range(ev.groupid);

			for (auto itrk = itr0.first; itrk != itr0.second; itrk++) {//kink cand loop
				std::cout << itrk->first << " " << itrk->second << std::endl;//gid,pl

				pl0 = itrk->second;//kink angle
				cnt = 0;
				double avetan[3] = { 0 };
				double avedrad[3] = { 0 };
				double avedlat[3] = { 0 };

				
				for (auto itr = c.base.begin(); itr != c.base.end(); itr++) {//d-->u
					std::cout << "pl" << itr->pl << std::endl;//for check

					//downstream side
					if (cnt != 0) {
						dax = itr->ax - ax;
						day = itr->ay - ay;
						dtan = sqrt(itr->ax * itr->ax + itr->ay * itr->ay) - tan;
						if (angle < 0.01) {
							dalat = dax;
							darad = day;
						}
						else {
							darad = (dax * ax + day * ay) / tan;
							dalat = (dax * ay - day * ax) / tan;
						}

					}
					upl = itr->pl;
					dpl = itr->pl - 1;
					ax = itr->ax;
					ay = itr->ay;
					tan = sqrt(ax * ax + ay * ay);



					//upstream side
					if (cnt != 0) {
						dax = itr->ax - ax;
						day = itr->ay - ay;
						dtan = sqrt(itr->ax * itr->ax + itr->ay * itr->ay) - tan;
						if (angle < 0.01) {
							dalat = dax;
							darad = day;
						}
						else {
							darad = (dax * ax + day * ay) / tan;
							dalat = (dax * ay - day * ax) / tan;
						}

					}
					upl = itr->pl;
					dpl = itr->pl - 1;
					ax = itr->ax;
					ay = itr->ay;
					tan = sqrt(ax * ax + ay * ay);


					//downstream side

					if (cnt != 0) {
						dax = itr->ax - ax;
						day = itr->ay - ay;
						dtan = sqrt(itr->ax * itr->ax + itr->ay * itr->ay) - angle;

						if (angle < 0.01) {
							dalat = dax;
							darad = day;
						}
						else {
							darad = (dax * ax + day * ay) / angle;
							dalat = (dax * ay - day * ax) / angle;
						}

						ofs << std::fixed << std::right
							<< std::setw(10) << std::setprecision(0) << ev.groupid << " "
							<< std::setw(4) << std::setprecision(0) << itr->pl << " ";

						if (itr->pl > pl0) {
							ofs << std::setw(3) << std::setprecision(0) << 1 << " ";
							avetan[0] = avetan[0] + fabs(dtan);
							avedrad[0] = avedrad[0] + fabs(darad);
							avedlat[0] = avedlat[0] + fabs(dalat);
							cnt1++;

						}
						else if (itr->pl == pl0) {
							ofs << std::setw(3) << std::setprecision(0) << 0 << " ";
							avetan[2] = fabs(dtan);
							avedrad[2] = fabs(darad);
							avedlat[2] = fabs(dalat);

						}
						else {
							ofs << std::setw(3) << std::setprecision(0) << -1 << " ";
							cnt2++;
							avetan[1] = avetan[1] + fabs(dtan);
							avedrad[1] = avedrad[1] + fabs(darad);
							avedlat[1] = avedlat[1] + fabs(dalat);

						}

						ofs << std::fixed << std::right
							<< std::setw(8) << std::setprecision(5) << dax << " "
							<< std::setw(8) << std::setprecision(5) << day << " "
							<< std::setw(8) << std::setprecision(5) << dtan << " "
							<< std::setw(8) << std::setprecision(5) << dalat << " "
							<< std::setw(8) << std::setprecision(5) << darad << " "
							<< std::endl;

					}
					ax = itr->ax;
					ay = itr->ay;
					angle = sqrt(ax * ax + ay * ay);

					cnt++;
				}
				avetan[0] = avetan[0] / cnt1;
				avedrad[0] = avedrad[0] / cnt1;
				avedlat[0] = avedlat[0] / cnt1;
				avetan[1] = avetan[1] / cnt2;
				avedrad[1] = avedrad[1] / cnt2;
				avedlat[1] = avedlat[1] / cnt2;


				ofs2 << std::fixed << std::right
					<< std::setw(10) << std::setprecision(0) << ev.groupid << " "
					<< std::setw(4) << std::setprecision(0) << pl0 << " "
					<< std::setw(8) << std::setprecision(5) << avetan[0] << " "
					<< std::setw(8) << std::setprecision(5) << avetan[2] << " "
					<< std::setw(8) << std::setprecision(5) << avetan[1] << " "
					<< std::setw(8) << std::setprecision(5) << avedrad[0] << " "
					<< std::setw(8) << std::setprecision(5) << avedrad[2] << " "
					<< std::setw(8) << std::setprecision(5) << avedrad[1] << " "
					<< std::setw(8) << std::setprecision(5) << avedlat[0] << " "
					<< std::setw(8) << std::setprecision(5) << avedlat[2] << " "
					<< std::setw(8) << std::setprecision(5) << avedlat[1] << " "
					<< std::endl;

				if (avedlat[2] > 0.02) {
					ofs1 << std::setw(10) << std::setprecision(0) << ev.groupid << " "
						<< std::setw(4) << std::setprecision(0) << pl0 << " "
						<< std::endl;
				}
				cnt1 = 0;
				cnt2 = 0;
			}

		}
	}
}

void Calc_divide_vph_inv(Momentum_recon::Mom_chain& c, int pl, double& mean, double& sd, int& count, double& vph) {
	//upstream-->downstream
	vph = -1;
	count = 0;
	double sum2 = 0, sum1 = 0;
	bool flg = true;
	for (auto itr = c.base.rbegin(); itr != c.base.rend(); itr++) {
		if (itr->pl > pl) {
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
void Calc_divide_pixel_inv(Momentum_recon::Mom_chain& c, int pl, double& mean, double& sd, int& count, double& hitnum) {
	hitnum = -1;
	count = 0;
	double sum2 = 0., sum1 = 0;
	bool flg = true;
	for (auto itr = c.base.rbegin(); itr != c.base.rend(); itr++) {
		if (itr->pl > pl) {
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
void Calc_divide_angle_lateral_inv(Momentum_recon::Mom_chain& c, int pl, double& dal, double& rms) {
	dal = NAN;

	double ax, ay, dax, day, dang, angle, sum2 = 0;
	int count = 0;
	for (auto itr = c.base_pair.rbegin(); itr != c.base_pair.rend(); itr++) {
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
		else if (itr->first.pl > pl) {
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
void Calc_divide_angle_radial_inv(Momentum_recon::Mom_chain& c, int pl, double& dar, double& rms) {
	dar = NAN;
	double ax, ay, dax, day, dang, angle, sum2 = 0;
	int count = 0;
	for (auto itr = c.base_pair.rbegin(); itr != c.base_pair.rend(); itr++) {
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
		else if (itr->first.pl > pl) {
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
