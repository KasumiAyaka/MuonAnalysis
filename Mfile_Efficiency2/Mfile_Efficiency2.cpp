#define _CRT_SECURE_NO_WARNINGS
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>
#include <set>
void chain_nseg_selection(const mfile1::MFile_minimum& m, std::vector<uint64_t>& all, std::vector<uint64_t>& sel, const int nseg);
void reject_Fe_ECC(const mfile1::MFile_minimum& m, std::vector<uint64_t>& all, std::vector<uint64_t>& sel);
void chain_dlat_selection(mfile1::MFile_minimum& m, std::vector<uint64_t>& all, std::vector<uint64_t>& sel, double threshold);
void chain_angle_selection(mfile1::MFile_minimum& m, std::vector<uint64_t>& all, std::vector<uint64_t>& sel, double thr_ax, double thr_ay);

void CalcEfficiency(std::ofstream& ofs_eff, std::ofstream& ofs_pred, std::ofstream& ofs_id, mfile1::MFile_minimum& m, std::vector<uint64_t>& all, int PL, std::map<int, double> zmap);
int angle_divide(double angle);


struct Prediction {
	int PL, flg, rawid;
	double ax, ay, x, y;
};
int main(int argc, char** argv) {
	if (argc != 6) {
		fprintf(stderr, "usage:prg in-mfile out-file(eff) out-file(prediction) out_exist_base_id_list PLmax\n");
		printf("%s\n", argv[0]);
		printf("%s\n", argv[1]);
		printf("%s\n", argv[2]);
		printf("%s\n", argv[3]);
		printf("%s\n", argv[4]);
		printf("%s\n", argv[5]);
		exit(1);
	}

	//input value
	std::string file_in_mfile = argv[1];
	std::string file_out_eff = argv[2];
	std::string file_out_prediction = argv[3];
	std::string file_out_id_list = argv[4];
	int PLmax = std::stoi(argv[5]);

	mfile1::MFile_minimum m;
	mfile1::read_mfile(file_in_mfile, m, 15);
	//zmap
	std::map<int, double> zmap;
	for (auto c : m.all_basetracks) {
		for (auto b : c) {
			zmap.insert(std::make_pair(b.pos / 10, b.z));
		}
	}
	std::vector<uint64_t> all, sel;
	for (int i = 0; i < m.chains.size(); i++) {
		all.push_back(i);
	}

	reject_Fe_ECC(m, all, sel);
	all = sel;
	sel.clear();

	chain_angle_selection(m, all, sel, 4.0, 4.0);;
	all = sel;
	sel.clear();

	chain_dlat_selection(m, all, sel, 0.005);

	std::ofstream ofs_eff(file_out_eff);
	std::ofstream ofs_pred(file_out_prediction);
	std::ofstream ofs_id(file_out_id_list);
	for (int PL = 4; PL <= PLmax; PL++) {
		CalcEfficiency(ofs_eff, ofs_pred, ofs_id, m, sel, PL, zmap);
	}
}
void reject_Fe_ECC(const mfile1::MFile_minimum& m, std::vector<uint64_t>& all, std::vector<uint64_t>& sel) {
	for (auto i : all) {
		if (m.chains[i].pos0 / 10 > 15) {
			sel.push_back(i);
		}
		else {
			int flg = 0;
			for (auto b : m.all_basetracks[i]) {
				if (b.pos / 10 == 16)flg++;
				if (b.pos / 10 == 17)flg++;
				if (b.pos / 10 == 18)flg++;
				if (b.pos / 10 == 19)flg++;
				if (b.pos / 10 == 20)flg++;
				if (b.pos / 10 == 21)flg++;
				if (b.pos / 10 == 22)flg++;
				if (b.pos / 10 == 23)flg++;
			}
			if (flg >= 6 || m.chains[i].pos1 / 10 > 23) {
				sel.push_back(i);
			}
		}
	}
	printf("reject_Fe_ECC(PL1 >= 24): %lld --> %lld (%4.1lf%%)\n", all.size(), sel.size(), sel.size() * 100. / all.size());
	return;
}
void chain_angle_selection(mfile1::MFile_minimum& m, std::vector<uint64_t>& all, std::vector<uint64_t>& sel, double thr_ax, double thr_ay) {
	double ax, ay;
	for (auto i : all) {
		ax = mfile1::chain_ax(m.all_basetracks[i]);
		ay = mfile1::chain_ay(m.all_basetracks[i]);
		if (thr_ax - fabs(ax) < 0)continue;
		if (thr_ay - fabs(ay) < 0)continue;
		sel.push_back(i);

	}
	printf("chain average angle selection |ax|<=%3.1lf |ay|<=%3.1lf : %lld --> %lld (%4.1lf%%)\n", thr_ax, thr_ay, all.size(), sel.size(), sel.size() * 100. / all.size());
	return;
}
void chain_dlat_selection(mfile1::MFile_minimum& m, std::vector<uint64_t>& all, std::vector<uint64_t>& sel, double threshold) {
	for (auto i : all) {
		if (mfile1::angle_diff_dev_lat(m.all_basetracks[i]) > threshold)continue;
		sel.push_back(i);
	}
	printf("chain lateral selection <= %5.4lf : %lld --> %lld (%4.1lf%%)\n", threshold, all.size(), sel.size(), sel.size() * 100. / all.size());
	return;

}

void CalcEfficiency(std::ofstream& ofs_eff, std::ofstream& ofs_pred, std::ofstream& ofs_id, mfile1::MFile_minimum& m, std::vector<uint64_t>& all, int PL, std::map<int, double> zmap) {
	std::vector<uint64_t> sel;
	int flg = 0;
	double dz;
	int PL_min = 3;
	int PL_max = 133;
	for (auto c : all) {
		flg = 0;
		for (auto b : m.all_basetracks[c]) {
			if (PL_min == PL) {
				if (abs(b.pos / 10 - PL) != 0 && abs(b.pos / 10 - PL) <= 6) {
					flg++;
				}
			}
			else if (PL_min + 1 == PL) {
				if (abs(b.pos / 10 - PL) != 0 && abs(b.pos / 10 - PL) <= 5) {
					flg++;
				}
			}
			else if (PL_min + 2 == PL) {
				if (abs(b.pos / 10 - PL) != 0 && abs(b.pos / 10 - PL) <= 4) {
					flg++;
				}
			}
			else if (PL_min + 3 <= PL && PL <= PL_max - 3) {
				if (abs(b.pos / 10 - PL) != 0 && abs(b.pos / 10 - PL) <= 3) {
					flg++;
				}
			}
			else if (PL_max - 2 == PL) {
				if (abs(b.pos / 10 - PL) != 0 && abs(b.pos / 10 - PL) <= 4) {
					flg++;
				}
			}
			else if (PL_max - 1 == PL) {
				if (abs(b.pos / 10 - PL) != 0 && abs(b.pos / 10 - PL) <= 5) {
					flg++;
				}
			}
			else if (PL_max == PL) {
				if (abs(b.pos / 10 - PL) != 0 && abs(b.pos / 10 - PL) <= 6) {
					flg++;
				}
			}

		}
		if (flg == 6) {
			sel.push_back(c);
		}
	}

	printf("PL=%d prediction selection: %lld --> %lld (%4.1lf%%)\n", PL, all.size(), sel.size(), sel.size() * 100. / all.size());
	if (sel.size() == 0) {
		printf("PL%03d no prediction track\n", PL);
		return;
	}
	std::vector<Prediction> pred;
	for (auto i : sel) {
		Prediction pred_tmp;
		pred_tmp.PL = PL;
		pred_tmp.ax = mfile1::chain_ax(m.all_basetracks[i]);
		pred_tmp.ay = mfile1::chain_ay(m.all_basetracks[i]);
		pred_tmp.flg = 0;
		pred_tmp.rawid = 0;
		for (auto b : m.all_basetracks[i]) {
			if ((PL >= 5 && PL <= 15) || (PL >= 16 && PL % 2 == 1)) {
				if (b.pos / 10 == PL - 1) {
					pred_tmp.x = b.x + pred_tmp.ax * (zmap[PL] - zmap[PL - 1]);
					pred_tmp.y = b.y + pred_tmp.ay * (zmap[PL] - zmap[PL - 1]);
				}
			}
			else if (PL == 4 || (PL >= 16 && PL % 2 == 0)) {
				if (b.pos / 10 == PL + 1) {
					pred_tmp.x = b.x + pred_tmp.ax * (zmap[PL] - zmap[PL + 1]);
					pred_tmp.y = b.y + pred_tmp.ay * (zmap[PL] - zmap[PL + 1]);
				}
			}
			else {
				fprintf(stderr, "PL exception PL%03d\n", PL);
			}

			if (b.pos / 10 == PL) {
				pred_tmp.flg = 1;
				pred_tmp.rawid = b.rawid;
			}
		}
		pred.push_back(pred_tmp);
	}


	double xmin, xmax, ymin, ymax;
	for (auto itr = pred.begin(); itr != pred.end(); itr++) {
		if (itr == pred.begin()) {
			xmin = itr->x;
			xmax = itr->x;
			ymin = itr->y;
			ymax = itr->y;
		}
		xmin = std::min(itr->x, xmin);
		xmax = std::max(itr->x, xmax);
		ymin = std::min(itr->y, ymin);
		ymax = std::max(itr->y, ymax);
	}
	xmin = xmin + 5000;
	xmax = xmax - 5000;
	ymin = ymin + 3000;
	ymax = ymax - 3000;

	std::vector<Prediction> pred_sel;
	for (auto itr = pred.begin(); itr != pred.end(); itr++) {
		if (itr->x < xmin)continue;
		if (itr->x > xmax)continue;
		if (itr->y < ymin)continue;
		if (itr->y > ymax)continue;
		pred_sel.push_back(*itr);
	}
	if (pred_sel.size() == 0) {
		printf("PL%03d no prediction track\n", PL);
		return;
	}


	const int NUM = 16;
	int all_c[NUM] = {}, count[NUM] = {};
	double eff[NUM], eff_err[NUM];
	double angle;
	for (auto itr = pred_sel.begin(); itr != pred_sel.end(); itr++) {
		angle = sqrt(itr->ax * itr->ax + itr->ay * itr->ay);
		all_c[angle_divide(angle)]++;
		if (itr->flg == 1) {
			count[angle_divide(angle)]++;
		}
		ofs_pred << std::right << std::fixed << std::setfill(' ')
			<< std::setw(3) << std::setprecision(0) << itr->PL << " "
			<< std::setw(7) << std::setprecision(4) << itr->ax << " "
			<< std::setw(7) << std::setprecision(4) << itr->ay << " "
			<< std::setw(8) << std::setprecision(1) << itr->x << " "
			<< std::setw(8) << std::setprecision(1) << itr->y << " "
			<< std::setw(2) << std::setprecision(0) << itr->flg << std::endl;
		if (itr->flg == 1) {
			ofs_id << std::right << std::fixed << std::setfill(' ')
				<< std::setw(3) << std::setprecision(0) << itr->PL << " "
				<< std::setw(12) << std::setprecision(0) << itr->rawid << " "
				<< std::setw(7) << std::setprecision(4) << itr->ax << " "
				<< std::setw(7) << std::setprecision(4) << itr->ay << " "
				<< std::setw(8) << std::setprecision(1) << itr->x << " "
				<< std::setw(8) << std::setprecision(1) << itr->y << std::endl;

		}
	}

	for (int i = 0; i < 13; i++) {
		eff[i] = double(count[i]) / all_c[i];
		eff_err[i] = sqrt(count[i] * eff[i] * (1 - eff[i])) / all_c[i];
		ofs_eff << std::right << std::fixed << std::setfill(' ')
			<< std::setw(4) << std::setprecision(0) << PL << " ";

		if (i == 0)ofs_eff << "0.0 0.1 ";
		if (i == 1)ofs_eff << "0.1 0.3 ";
		if (i == 2)ofs_eff << "0.3 0.5 ";
		if (i == 3)ofs_eff << "0.5 0.7 ";
		if (i == 4)ofs_eff << "0.7 0.9 ";
		if (i == 5)ofs_eff << "0.9 1.1 ";
		if (i == 6)ofs_eff << "1.1 1.3 ";
		if (i == 7)ofs_eff << "1.3 1.5 ";
		if (i == 8)ofs_eff << "1.5 2.0 ";
		if (i == 9)ofs_eff << "2.0 2.5 ";
		if (i == 10)ofs_eff << "2.5 3.0 ";
		if (i == 11)ofs_eff << "3.0 3.5 ";
		if (i == 12)ofs_eff << "3.5 4.0 ";
		if (i == 13)ofs_eff << "4.0 4.5 ";
		if (i == 14)ofs_eff << "4.5 5.0 ";
		if (i == 15)ofs_eff << "5.0 10.0 ";

		ofs_eff << std::setw(8) << std::setprecision(0) << all_c[i] << " "
			<< std::setw(8) << std::setprecision(0) << count[i] << " "
			<< std::setw(5) << std::setprecision(4) << eff[i] << " "
			<< std::setw(5) << std::setprecision(4) << eff_err[i] << std::endl;

	}

}
int angle_divide(double angle) {
	if (angle < 0.1)return 0;
	if (0.1 <= angle && angle < 0.3)return 1;
	if (0.3 <= angle && angle < 0.5)return 2;
	if (0.5 <= angle && angle < 0.7)return 3;
	if (0.7 <= angle && angle < 0.9)return 4;
	if (0.9 <= angle && angle < 1.1)return 5;
	if (1.1 <= angle && angle < 1.3)return 6;
	if (1.3 <= angle && angle < 1.5)return 7;
	if (1.5 <= angle && angle < 2.0)return 8;
	if (2.0 <= angle && angle < 2.5)return 9;
	if (2.5 <= angle && angle < 3.0)return 10;
	if (3.0 <= angle && angle < 3.5)return 11;
	if (3.5 <= angle && angle < 4.0)return 12;
	if (4.0 <= angle && angle < 4.5)return 13;
	if (4.5 <= angle && angle < 5.0)return 14;
	if (5.0 <= angle)return 15;
	return -1;
}
