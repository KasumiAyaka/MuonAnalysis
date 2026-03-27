#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.h>
#pragma comment(lib, "functions.lib")
#include <functions.h>
#pragma comment(lib, "VxxReader.lib")
#include "VxxReader.h"

#include <chrono>
#include <filesystem>
#include <set>
#include "Cycle_enumerate.h"

class align_param {
public:
	int id;
	//視野中心
	double x, y, z;
	//parameter(9);
	double dx, dy, dz, x_rot, y_rot, z_rot, x_shrink, y_shrink, z_shrink, yx_shear, zx_shear, zy_shear;

};

std::string Set_file_read_ali_path(std::string file_ECC_path, int pl0, int pl1, int area);
std::string Set_file_read_ali_path_inv(std::string file_ECC_path, int pl0, int pl1, int area);

void signal_region(std::vector <corrmap_3d::align_param >& corr, double& center, double& sigma);
std::vector <corrmap_3d::align_param > select_high_signal(std::vector <corrmap_3d::align_param >& corr);
std::vector <corrmap_3d::align_param > select_edge(std::vector <corrmap_3d::align_param >& corr);
void remove_single_path(std::vector <corrmap_3d::align_param >& corr);
void remove_all_connect_path(std::vector <corrmap_3d::align_param >& corr);
void divide_connect_part(std::vector<std::set<int>>& connect, std::set<int>& point, std::multimap<int, int>& path_all);
void select_long_cycle(std::vector <corrmap_3d::align_param >& corr);
double Calc_area(std::vector <corrmap_3d::align_param >& corr, bool output);
void output_Fiducial_Area_line(std::ofstream& ofs, std::vector<corrmap_3d::align_param>& corr, int PL);

int main(int argc, char** argv) {
	if (argc != 5) {
		fprintf(stderr, "prg ECC_path plmin plmax output-file\n");
		exit(1);
	}
	int area = 0;
	int plmin = std::stoi(argv[2]);
	int plmax = std::stoi(argv[3]);

	std::string file_ECC_path = argv[1];
	std::string file_output = argv[4];


	std::map<int, std::vector<corrmap_3d::align_param>> all_PL_area;
	for (int pl = plmin; pl <= plmax; pl++) {

		//fileの存在確認・読み込み
		std::string file_in_ali;
		if (pl != plmax) {
			bool flg = false;
			for (int i_pl = 1; i_pl <= 3; i_pl++) {
				file_in_ali = Set_file_read_ali_path_inv(file_ECC_path, pl + i_pl, pl, area);
				if (!std::filesystem::exists(file_in_ali)) {
					fprintf(stderr, "%s not exist\n", file_in_ali.c_str());
					continue;
				}
				flg = true;
				break;
			}
			if (!flg)continue;
			printf("start PL%03d\n", pl);
			std::vector <corrmap_3d::align_param > corr = corrmap_3d::read_ali_param(file_in_ali, false);
			//signalでカット
			corr = select_high_signal(corr);
			remove_single_path(corr);
			remove_all_connect_path(corr);
			select_long_cycle(corr);

			all_PL_area.insert(std::make_pair(pl, corr));
		}
		else {

			bool flg = false;
			for (int i_pl = 1; i_pl <= 3; i_pl++) {
				file_in_ali = Set_file_read_ali_path(file_ECC_path, pl - i_pl, pl, area);
				if (!std::filesystem::exists(file_in_ali)) {
					fprintf(stderr, "%s not exist\n", file_in_ali.c_str());
					continue;
				}
				flg = true;
				break;
			}
			if (!flg)continue;

			printf("start PL%03d\n", pl);
			std::vector <corrmap_3d::align_param > corr = corrmap_3d::read_ali_param(file_in_ali, false);
			//signalでカット
			corr = select_high_signal(corr);
			remove_single_path(corr);
			remove_all_connect_path(corr);
			select_long_cycle(corr);

			all_PL_area.insert(std::make_pair(pl, corr));
		}
	}
	std::ofstream ofs(file_output);

	for (auto itr = all_PL_area.begin(); itr != all_PL_area.end(); itr++) {
		output_Fiducial_Area_line(ofs, itr->second, itr->first);
	}
	ofs.close();

}
std::string Set_file_read_ali_path(std::string file_ECC_path, int pl0, int pl1, int area) {
	std::stringstream file_in_ali;
	file_in_ali << file_ECC_path << "\\Area" << area << "\\0\\align\\fine\\ali_"
		<< std::setw(3) << std::setfill('0') << pl0 << "_"
		<< std::setw(3) << std::setfill('0') << pl1 << "_loop_multi.txt";
	return file_in_ali.str();
	//ali_003_004_loop_multi.txt
}

std::string Set_file_read_ali_path_inv(std::string file_ECC_path, int pl0, int pl1, int area) {
	std::stringstream file_in_ali;
	file_in_ali << file_ECC_path << "\\Area" << area << "\\0\\align_inv\\fine\\ali_"
		<< std::setw(3) << std::setfill('0') << pl0 << "_"
		<< std::setw(3) << std::setfill('0') << pl1 << "_loop_multi.txt";
	return file_in_ali.str();
	//ali_003_004_loop_multi.txt
}

void signal_region(std::vector <corrmap_3d::align_param >& corr, double& center, double& sigma) {

	int count = 0;
	double p_mean = 0, p_rms = DBL_MAX, mean = 0, rms = 0, rms_ratio = 0;
	int loop_num = 0;
	while (true) {

		count = 0;
		mean = 0;
		rms = 0;

		for (auto itr = corr.begin(); itr != corr.end(); itr++) {
			if (itr->signal < 1)continue;
			if (p_mean - 3 * p_rms < itr->signal && itr->signal < p_mean + 3 * p_rms) {
				count++;
				mean += itr->signal;
				rms += itr->signal * itr->signal;
			}

		}
		if (count < 2) {
			fprintf(stderr, "signal distributiuon exception(count = %d)\n", count);
			exit(1);
		}
		mean = mean / count;
		rms = sqrt(rms / count - mean * mean);
		rms_ratio = rms / p_rms;
		printf("count=%d %g +- %g (%g/%g = %g)\n", count, mean, rms, rms, p_rms, rms_ratio);

		if (fabs(1 - rms_ratio) < 0.1) {
			break;
		}
		p_rms = rms;
		p_mean = mean;
		loop_num++;
		if (loop_num > 100) {
			fprintf(stderr, "signal distributiuon exception\n");
			exit(1);
		}
	}
	center = mean;
	sigma = rms;
}

std::vector <corrmap_3d::align_param > select_high_signal(std::vector <corrmap_3d::align_param >& corr) {
	std::vector <corrmap_3d::align_param >ret;

	double signal_center, signal_rms, signal_min, signal_max;
	signal_region(corr, signal_center, signal_rms);
	printf("signal %.1lf %.1lf\n", signal_center, signal_rms);

	signal_min = signal_center - signal_rms * 3;
	signal_max = signal_center + signal_rms * 4;

	for (auto itr = corr.begin(); itr != corr.end(); itr++) {
		if (itr->signal < signal_min)continue;
		if (itr->signal > signal_max)continue;
		ret.push_back(*itr);
	}
	return ret;
}
std::vector <corrmap_3d::align_param > select_edge(std::vector <corrmap_3d::align_param >& corr) {
	std::multimap<int, corrmap_3d::align_param> map_x, map_y;
	std::map<std::pair<int, int>, corrmap_3d::align_param> edge_all;

	double x_ave = 0, y_ave = 0;
	for (auto itr = corr.begin(); itr != corr.end(); itr++) {
		map_x.insert(std::make_pair(itr->ix, *itr));
		map_y.insert(std::make_pair(itr->iy, *itr));
		x_ave += itr->x;
		y_ave += itr->y;
	}
	x_ave = x_ave / corr.size();
	y_ave = y_ave / corr.size();

	corrmap_3d::align_param edge_min, edge_max;
	//ixについてloop iyのmin/maxを抽出
	for (auto itr = map_x.begin(); itr != map_x.end(); itr++) {
		auto range = map_x.equal_range(itr->first);
		for (auto res = range.first; res != range.second; res++) {
			if (res == range.first) {
				edge_min = res->second;
				edge_max = res->second;
			}
			if (edge_min.iy > itr->second.iy)edge_min = res->second;
			if (edge_max.iy < res->second.iy)edge_max = res->second;
		}
		edge_all.insert(std::make_pair(std::make_pair(edge_min.ix, edge_min.iy), edge_min));
		edge_all.insert(std::make_pair(std::make_pair(edge_max.ix, edge_max.iy), edge_max));
	}
	//iyについてloop ixのmin/maxを抽出
	for (auto itr = map_y.begin(); itr != map_y.end(); itr++) {
		auto range = map_y.equal_range(itr->first);
		for (auto res = range.first; res != range.second; res++) {
			if (res == range.first) {
				edge_min = res->second;
				edge_max = res->second;
			}
			if (edge_min.ix > res->second.ix)edge_min = res->second;
			if (edge_max.ix < res->second.ix)edge_max = res->second;
		}
		edge_all.insert(std::make_pair(std::make_pair(edge_min.ix, edge_min.iy), edge_min));
		edge_all.insert(std::make_pair(std::make_pair(edge_max.ix, edge_max.iy), edge_max));
	}

	//edgeの点を中心から回して順番につなげる

	double base_x, base_y;
	base_x = edge_all.begin()->second.x;
	base_y = edge_all.begin()->second.y;
	double base_vec[2], dir_vec[2], base_length;
	base_vec[0] = base_x - x_ave;
	base_vec[1] = base_y - y_ave;
	base_length = sqrt(base_vec[0] * base_vec[0] + base_vec[1] * base_vec[1]);

	int count = 0;
	std::map<double, corrmap_3d::align_param> edge_rot;
	double dir_length, theta, c_product;
	for (auto itr = edge_all.begin(); itr != edge_all.end(); itr++) {
		dir_vec[0] = itr->second.x - x_ave;
		dir_vec[1] = itr->second.y - y_ave;
		dir_length = sqrt(dir_vec[0] * dir_vec[0] + dir_vec[1] * dir_vec[1]);
		theta = (base_vec[0] * dir_vec[0] + base_vec[1] * dir_vec[1]) / (base_length * dir_length);
		if (1 - fabs(theta) < 0.000000001) {
			//一直線上の場合
			//始点との距離を測る
			double dist = sqrt(pow(base_vec[0] - dir_vec[0], 2) + pow(base_vec[1] - dir_vec[1], 2));
			//中心からの距離<始点からの距離 -->0度
			if (dist < dir_length) {
				theta = 0;
			}
			else {
				theta = acos(-1);
			}
		}
		else {
			theta = acos(theta);
			c_product = base_vec[0] * dir_vec[1] - base_vec[1] * dir_vec[0];
			if (signbit(c_product)) {
				theta = 2 * acos(-1) - theta;
			}
		}
		count++;
		edge_rot.insert(std::make_pair(theta, itr->second));
	}




	std::vector<corrmap_3d::align_param> ret;
	for (auto itr = edge_rot.begin(); itr != edge_rot.end(); itr++) {
		//printf("%5.4lf\n", itr->first);
		ret.push_back(itr->second);
	}

	return ret;

}

void remove_point(std::vector <corrmap_3d::align_param >& corr) {
	std::multimap<std::pair<int, int>, std::pair<int, corrmap_3d::align_param>> corr_map;
	for (auto itr = corr.begin(); itr != corr.end(); itr++) {
		corr_map.insert(std::make_pair(std::make_pair(itr->ix, itr->iy), std::make_pair(0, *itr)));
	}
	std::pair<int, int> id;
	for (auto itr = corr.begin(); itr != corr.end(); itr++) {
		auto res0 = corr_map.find(std::make_pair(itr->ix - 1, itr->iy));
		auto res1 = corr_map.find(std::make_pair(itr->ix + 1, itr->iy));
		auto res2 = corr_map.find(std::make_pair(itr->ix, itr->iy - 1));
		auto res3 = corr_map.find(std::make_pair(itr->ix, itr->iy + 1));

		if (res0 != corr_map.end())res0->second.first++;
		if (res1 != corr_map.end())res1->second.first++;
		if (res2 != corr_map.end())res2->second.first++;
		if (res3 != corr_map.end())res3->second.first++;
	}
	std::vector <corrmap_3d::align_param > sel;
	for (auto itr = corr_map.begin(); itr != corr_map.end(); itr++) {
		if (itr->second.first <= 1)continue;
		sel.push_back(itr->second.second);
	}

	std::swap(corr, sel);
	//return corr.size();
}
void remove_single_path(std::vector <corrmap_3d::align_param >& corr) {
	int point_num = -1;
	while (corr.size() != point_num) {
		point_num = corr.size();
		remove_point(corr);
	}

}
void remove_all_connect_path(std::vector <corrmap_3d::align_param >& corr) {
	std::multimap<std::pair<int, int>, std::pair<int, corrmap_3d::align_param>> corr_map;
	for (auto itr = corr.begin(); itr != corr.end(); itr++) {
		corr_map.insert(std::make_pair(std::make_pair(itr->ix, itr->iy), std::make_pair(0, *itr)));
	}
	std::pair<int, int> id;
	for (auto itr = corr.begin(); itr != corr.end(); itr++) {
		auto res0 = corr_map.find(std::make_pair(itr->ix - 1, itr->iy - 1));
		auto res1 = corr_map.find(std::make_pair(itr->ix - 1, itr->iy));
		auto res2 = corr_map.find(std::make_pair(itr->ix - 1, itr->iy + 1));
		auto res3 = corr_map.find(std::make_pair(itr->ix, itr->iy - 1));
		auto res4 = corr_map.find(std::make_pair(itr->ix, itr->iy + 1));
		auto res5 = corr_map.find(std::make_pair(itr->ix + 1, itr->iy - 1));
		auto res6 = corr_map.find(std::make_pair(itr->ix + 1, itr->iy));
		auto res7 = corr_map.find(std::make_pair(itr->ix + 1, itr->iy + 1));

		if (res0 != corr_map.end())res0->second.first++;
		if (res1 != corr_map.end())res1->second.first++;
		if (res2 != corr_map.end())res2->second.first++;
		if (res3 != corr_map.end())res3->second.first++;
		if (res4 != corr_map.end())res4->second.first++;
		if (res5 != corr_map.end())res5->second.first++;
		if (res6 != corr_map.end())res6->second.first++;
		if (res7 != corr_map.end())res7->second.first++;
	}
	std::vector <corrmap_3d::align_param > sel;
	for (auto itr = corr_map.begin(); itr != corr_map.end(); itr++) {
		if (itr->second.first == 8)continue;
		sel.push_back(itr->second.second);
	}

	std::swap(corr, sel);
	//return corr.size();


}

void divide_connect_part(std::vector<std::set<int>>& connect, std::set<int>& point, std::multimap<int, int>& path_all) {
	std::set<int> finished;

	for (auto itr = point.begin(); itr != point.end(); itr++) {
		if (finished.count(*itr) == 1)continue;
		std::set<int> connect_part, add, now;
		add.insert(*itr);
		while (add.size() != 0) {
			now = add;
			add.clear();
			for (auto itr2 = now.begin(); itr2 != now.end(); itr2++) {
				if (path_all.count(*itr2) == 0)continue;
				auto range = path_all.equal_range(*itr2);
				for (auto res = range.first; res != range.second; res++) {
					if (connect_part.count(res->second) == 1)continue;
					if (add.count(res->second) == 1)continue;
					add.insert(res->second);
				}
			}
			for (auto itr2 = now.begin(); itr2 != now.end(); itr2++) {
				connect_part.insert(*itr2);
			}
		}

		for (auto itr = connect_part.begin(); itr != connect_part.end(); itr++) {
			finished.insert(*itr);
		}
		connect.push_back(connect_part);
	}
}

void select_long_cycle(std::vector <corrmap_3d::align_param >& corr) {
	std::map<int, corrmap_3d::align_param> corr_map;
	std::set<int> point;
	int num = 0;
	//corrmapに通し番号
	for (auto itr = corr.begin(); itr != corr.end(); itr++) {
		corr_map.insert(std::make_pair(num, *itr));
		point.insert(num);
		num++;
	}
	//隣接をつなぐ
	std::multimap<int, int> path_all, path_next, path_prev;
	for (auto itr = corr_map.begin(); itr != corr_map.end(); itr++) {
		for (auto itr2 = std::next(itr, 1); itr2 != corr_map.end(); itr2++) {
			int dix = abs(itr->second.ix - itr2->second.ix);
			int diy = abs(itr->second.iy - itr2->second.iy);
			if (dix + diy == 1) {
				path_all.insert(std::make_pair(itr->first, itr2->first));
				path_all.insert(std::make_pair(itr2->first, itr->first));
				path_next.insert(std::make_pair(itr->first, itr2->first));
				path_prev.insert(std::make_pair(itr2->first, itr->first));
			}
		}
	}
	//連結成分に分解
	std::vector<std::set<int>> connect;
	divide_connect_part(connect, point, path_all);

	//最大の連結部分を抽出
	std::set<int> connect_max;
	for (auto itr = connect.begin(); itr != connect.end(); itr++) {
		if (connect_max.size() < itr->size()) {
			connect_max = *itr;
		}
		//printf("connect part point %d\n", itr->size());
	}
	std::multimap<int, int> path_next_sel;
	for (auto itr = path_next.begin(); itr != path_next.end(); itr++) {
		if (connect_max.count(itr->first) == 0)continue;
		if (connect_max.count(itr->second) == 0)continue;
		path_next_sel.insert(*itr);
	}
	//閉路の列挙
	std::vector<std::vector<std::pair<int, int>>> all_cycle = cycle_enumerate(path_next_sel);
	std::vector<std::vector <corrmap_3d::align_param >>corr_cycle;
	for (int i = 0; i < all_cycle.size(); i++) {
		std::vector <corrmap_3d::align_param > cycle;
		for (int j = 0; j < all_cycle[i].size(); j++) {
			//printf("%d - %d\n", all_cycle[i][j].first, all_cycle[i][j].second);
			cycle.push_back(corr_map.at(all_cycle[i][j].first));
		}
		corr_cycle.push_back(cycle);
		//printf("cycle %d path %d\n", i, all_cycle[i].size());
	}
	std::vector <corrmap_3d::align_param > select_cycle;
	double area = 0;
	for (int i = 0; i < corr_cycle.size(); i++) {
		if (Calc_area(corr_cycle[i], false) > area) {
			select_cycle = corr_cycle[i];
			area = Calc_area(corr_cycle[i], false);
		}
	}
	corr = select_cycle;
}

double Calc_area(std::vector <corrmap_3d::align_param >& corr, bool output) {
	double area = 0;
	int n = corr.size();
	for (int i = 0; i <= n; i++) {
		area += corr[i % n].x * corr[(i + 1) % n].y - corr[i % n].y * corr[(i + 1) % n].x;
	}
	area = fabs(area / 2);
	area = area / (10 * 10 * 1000 * 1000);
	if (output) {
		printf("area %.4lf cm2\n", area);
		printf("ratio %g\n", area / (25 * 25));
	}
	return area;
}

void output_Fiducial_Area_line(std::ofstream& ofs, std::vector<corrmap_3d::align_param>& corr, int PL) {
	//std::ofstream ofs(filename);

	if (corr.size() == 0) {
		fprintf(stderr, "target corrmap ... null\n");
	}
	else {
		int count = 0;
		int n = corr.size();
		std::cout << std::right << std::fixed;
		for (int i = 0; i < n; i++) {
			if (count % 10000 == 0) {
				fprintf(stderr, "\r Write corrmap ... %d/%d (%4.1lf%%)", count, int(corr.size()), count * 100. / corr.size());
			}
			count++;
			ofs << std::right << std::fixed
				<< std::setw(5) << std::setprecision(0) << PL << " "
				<< std::setw(8) << std::setprecision(1) << corr[i % n].x << " "
				<< std::setw(8) << std::setprecision(1) << corr[i % n].y << " "
				<< std::setw(8) << std::setprecision(1) << corr[i % n].dz << " "
				<< std::setw(8) << std::setprecision(1) << corr[(i + 1) % n].x << " "
				<< std::setw(8) << std::setprecision(1) << corr[(i + 1) % n].y << " "
				<< std::setw(8) << std::setprecision(1) << corr[(i + 1) % n].dz << std::endl;
		}
		fprintf(stderr, "\r Write corrmap ... %d/%d (%4.1lf%%)\n", count, int(corr.size()), count * 100. / corr.size());

	}
}
