#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

//class Point {
//public:
//	double x, y, z;
//};
//class Fiducial_Area {
//public:
//	int pl;
//	Point p[2];
//	//double x0, y0, z0, x1, y1, z1;
//};

class Difference {
public:
	int groupid, chainid, vertex_pl, target_pl[2], direction, add_nseg, edge_pl[2], edge_out_flg;
	double ax, ay, dax, day, dx, dy, dar, dal, dr, dl;
	double  s_dar, s_dal, s_dr, s_dl;
};
class Sigma_list {
public:
	int peke;
	double angle;
	double dal, dar, dr, dl;
};
std::vector<std::pair<Momentum_recon::Mom_chain, std::vector<Momentum_recon::Mom_chain>>> divide_event(std::vector<Momentum_recon::Mom_chain>& momch);
std::vector<Momentum_recon::Mom_chain > format_change(std::vector<std::pair<Momentum_recon::Mom_chain, std::vector<Momentum_recon::Mom_chain>>>& mom_ev);
void Check_fill_factor_forward_water(std::string filename, std::vector<std::pair<Momentum_recon::Mom_chain, std::vector<Momentum_recon::Mom_chain>>>& momch_ev);
std::vector<Difference> Calc_difference(std::vector<Momentum_recon::Event_information>& momch_ev, std::map<int, std::vector<Fiducial_Area::Fiducial_Area>>& area);
void Calc_difference(Difference& diff, std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack>& pair);
void output_difference(std::string filename, std::vector<Difference>& diff);
std::map<int, std::map<double, Sigma_list>> Read_sigma_list(std::string filename);
void Calc_sigma(std::vector<Difference>& diff, std::map<int, std::map<double, Sigma_list>>& sigma_list);
std::set<std::pair<int, int>> Get_penetrate_list(std::vector<Difference>& diff);
std::map<double, std::pair<double, double>> vph_mip_distribution(std::string in);
bool judege_Black(Momentum_recon::Mom_chain& chain, std::map<double, std::pair<double, double>>& vph_mip, double thr_sigma);
std::vector<Momentum_recon::Event_information > cut_penetrate(std::vector<Momentum_recon::Event_information >& momch_ev, std::set<std::pair<int, int>>& id_list, std::string in);
void trans_base_all(std::vector < std::pair<Fiducial_Area::Point*, corrmap_3d::align_param2*>>& track_pair);
std::vector <std::pair<Fiducial_Area::Point*, corrmap_3d::align_param2*>>track_affineparam_correspondence(std::vector<Fiducial_Area::Point*>& p, std::vector <corrmap_3d::align_param2>& param);
void trans_mfile_cordinate(std::vector<corrmap_3d::align_param2>& param, std::vector<Fiducial_Area::Fiducial_Area>& area);
bool judge_fiducial_area(std::vector<Fiducial_Area::Fiducial_Area>& area, Momentum_recon::Mom_basetrack& b);
void trans_mfile_cordinate(std::vector<corrmap_3d::align_param2>& param, std::vector<Fiducial_Area::Fiducial_Area>& area);
std::map<int, std::vector<Fiducial_Area::Fiducial_Area>> read_fiducial_Area(std::string filename);
int judge_momch_edgeout(int direction, Momentum_recon::Mom_chain& chain, std::map<int, std::vector<Fiducial_Area::Fiducial_Area>>& area, double edge_cut, int ex_pl_max);

int main(int argc, char** argv) {
	if (argc != 7) {
		fprintf(stderr, "usage:file-in-momch sigma_list file-in-ECC fa.txt mip_distribution.txt  file-out file_out_momch\n");
		exit(1);
	}
	std::string file_in_momch = argv[1];
	std::string file_in_sigma = argv[2];
	std::string file_in_ECC = argv[3];
	std::string file_in_area = argv[4];

	std::string input = argv[5];
	std::string file_out_difference = argv[6];
	std::string file_out_momch = argv[7];

	//corrmap absの読み込み
	std::string file_in_corrmap = file_in_ECC + "\\Area0\\0\\align\\fine\\local\\corrmap-local-abs.lst";
	std::map<int, std::vector<corrmap_3d::align_param>>corrmap = corrmap_3d::read_ali_param_abs(file_in_corrmap, 1);
	std::map<int, std::vector<corrmap_3d::align_param2>>corrmap_dd = corrmap_3d::DelaunayDivide_map(corrmap);

	//fiducial areaの読み込み
	std::map<int, std::vector<Fiducial_Area::Fiducial_Area>> area = read_fiducial_Area(file_in_area);
	for (auto itr = area.begin(); itr != area.end(); itr++) {
		if (corrmap_dd.count(itr->first) == 0) {
			fprintf(stderr, "corrmap local abs PL%03d not found\n", itr->first);
			exit(1);
		}
		std::vector<corrmap_3d::align_param2> param = corrmap_dd.at(itr->first);
		Fiducial_Area::trans_mfile_cordinate(param, itr->second);
	}


	std::map<int, std::map<double, Sigma_list>> sigma_list = Read_sigma_list(file_in_sigma);

	std::vector<Momentum_recon::Event_information> momch_ev = Momentum_recon::Read_Event_information_extension(file_in_momch);

	std::vector<Difference>diff = Calc_difference(momch_ev, area);
	Calc_sigma(diff, sigma_list);
	output_difference(file_out_difference, diff);
	std::set<std::pair<int, int>> penetrate_id = Get_penetrate_list(diff);
	momch_ev = cut_penetrate(momch_ev, penetrate_id,input);
	Momentum_recon::Write_Event_information_extension(file_out_momch, momch_ev);
}

std::map<int, std::vector<Fiducial_Area::Fiducial_Area>> read_fiducial_Area(std::string filename) {

	std::ifstream ifs(filename);
	std::multimap<int, Fiducial_Area::Fiducial_Area> fa_multi;
	std::map<int, std::vector<Fiducial_Area::Fiducial_Area>> ret;
	Fiducial_Area::Fiducial_Area fa;
	while (ifs >> fa.pl >> fa.p[0].x >> fa.p[0].y >> fa.p[0].z >> fa.p[1].x >> fa.p[1].y >> fa.p[1].z) {
		fa_multi.insert(std::make_pair(fa.pl, fa));
	}

	int count = 0;
	for (auto itr = fa_multi.begin(); itr != fa_multi.end(); itr++) {
		count = fa_multi.count(itr->first);
		auto range = fa_multi.equal_range(itr->first);
		std::vector<Fiducial_Area::Fiducial_Area> fa_vec;
		for (auto itr2 = range.first; itr2 != range.second; itr2++) {
			fa_vec.push_back(itr2->second);
		}
		ret.insert(std::make_pair(itr->first, fa_vec));
	}

	return ret;

}

void trans_mfile_cordinate(std::vector<corrmap_3d::align_param2>& param, std::vector<Fiducial_Area::Fiducial_Area>& area) {

	std::vector< Fiducial_Area::Point*> p_trans;
	for (auto itr = area.begin(); itr != area.end(); itr++) {
		p_trans.push_back(&(itr->p[0]));
		p_trans.push_back(&(itr->p[1]));
	}
	std::vector <std::pair<Fiducial_Area::Point*, corrmap_3d::align_param2*>> p_trans_map = Fiducial_Area::track_affineparam_correspondence(p_trans, param);
	Fiducial_Area::trans_base_all(p_trans_map);
}



std::map<int, std::map<double, Sigma_list>> Read_sigma_list(std::string filename) {

	std::map<int, std::map<double, Sigma_list>> ret;
	std::ifstream ifs(filename);
	int peke, num;
	while (ifs >> peke >> num) {
		std::map<double, Sigma_list> sigma_map;
		Sigma_list sig;
		sig.peke = peke;
		double angle_min, angle_max;
		for (int i = 0; i < num; i++) {
			ifs >> sig.angle >> sig.dal >> sig.dar >> sig.dl >> sig.dr;
			sigma_map.insert(std::make_pair(sig.angle, sig));
		}
		ret.insert(std::make_pair(sig.peke, sigma_map));
	}
	return ret;

	//for (auto itr = ret.begin(); itr != ret.end(); itr++) {
	//	printf("peke:%d\n", itr->first);
	//	for (auto s : itr->second) {
	//		printf("%.1lf %.5lf %.5lf %.2lf %.2lf\n", s.second.angle, s.second.dal, s.second.dar, s.second.dl, s.second.dr);
	//	}
	//}

}

std::vector<Difference> Calc_difference(std::vector<Momentum_recon::Event_information>& momch_ev, std::map<int, std::vector<Fiducial_Area::Fiducial_Area>>& area) {
	std::vector<Difference> ret;

	for (auto& ev : momch_ev) {
		Difference diff;
		diff.groupid = ev.groupid;

		for (auto& c : ev.chains) {
			if (c.chainid == 0)continue;
			diff.vertex_pl = ev.vertex_pl;
			diff.direction = c.direction;
			diff.chainid = c.chainid;
			diff.edge_pl[0] = c.base.begin()->pl;
			diff.edge_pl[1] = c.base.rbegin()->pl;

			diff.add_nseg = 0;
			diff.target_pl[0] = -1;
			diff.target_pl[1] = -1;
			//forward
			if (diff.direction == 1) {
				diff.target_pl[0] = diff.vertex_pl;

				for (auto& b : c.base) {
					if (b.pl > diff.vertex_pl) {
						if (diff.target_pl[1] < 0) {
							diff.target_pl[1] = b.pl;
						}
						diff.add_nseg += 1;
					}
				}
				if (diff.target_pl[1] < 0)continue;
				for (auto& pair : c.base_pair) {
					if (pair.first.pl == diff.target_pl[0] && pair.second.pl == diff.target_pl[1]) {
						Calc_difference(diff, pair);
						ret.push_back(diff);
						break;
					}
				}

			}
			//backward
			else if (diff.direction == -1) {
				for (auto b = c.base.rbegin(); b != c.base.rend(); b++) {
					if (b->pl <= diff.vertex_pl) {
						if (diff.target_pl[1] < 0) {
							diff.target_pl[0] = b->pl;
							diff.target_pl[1] = std::next(b, -1)->pl;
						}
						diff.add_nseg += 1;
					}
				}
				if (diff.target_pl[1] < 0)continue;
				for (auto& pair : c.base_pair) {
					if (pair.first.pl == diff.target_pl[0] && pair.second.pl == diff.target_pl[1]) {
						Calc_difference(diff, pair);
						ret.push_back(diff);
						break;
					}
				}

			}

			diff.edge_out_flg = judge_momch_edgeout(diff.direction, c, area, 0, 4);

		}
	}
	return ret;

}
void Calc_sigma(std::vector<Difference>& diff, std::map<int, std::map<double, Sigma_list>>& sigma_list) {

	int peke;
	double angle, sigma_dal, sigma_dar, sigma_dr, sigma_dl;
	Sigma_list sigma[2];
	for (auto itr = diff.begin(); itr != diff.end(); itr++) {
		peke = itr->target_pl[1] - itr->target_pl[0] - 1;
		angle = sqrt(pow(itr->ax, 2) + pow(itr->ay, 2));

		auto sigma_peke = sigma_list.find(peke);
		if (sigma_peke == sigma_list.end()) {
			fprintf(stderr, "peke=%d not found\n", peke);
			exit(1);
		}
		auto sigma_angle = sigma_peke->second.upper_bound(angle);
		if (sigma_angle == sigma_peke->second.begin()) {
			sigma[0] = sigma_angle->second;
			sigma[1] = sigma_angle->second;
			sigma_dal = sigma[0].dal;
			sigma_dar = sigma[0].dar;
			sigma_dl = sigma[0].dl;
			sigma_dr = sigma[0].dr;
		}
		else if (sigma_angle == sigma_peke->second.end()) {
			sigma[0] = sigma_peke->second.rbegin()->second;
			sigma[1] = sigma_peke->second.rbegin()->second;
			sigma_dal = sigma[0].dal;
			sigma_dar = sigma[0].dar;
			sigma_dl = sigma[0].dl;
			sigma_dr = sigma[0].dr;
		}
		else {
			sigma[0] = std::next(sigma_angle, -1)->second;
			sigma[1] = sigma_angle->second;
			sigma_dal = (sigma[1].dal - sigma[0].dal) / (sigma[1].angle - sigma[0].angle) * (angle - sigma[0].angle) + sigma[0].dal;
			sigma_dar = (sigma[1].dar - sigma[0].dar) / (sigma[1].angle - sigma[0].angle) * (angle - sigma[0].angle) + sigma[0].dar;
			sigma_dl = (sigma[1].dl - sigma[0].dl) / (sigma[1].angle - sigma[0].angle) * (angle - sigma[0].angle) + sigma[0].dl;
			sigma_dr = (sigma[1].dr - sigma[0].dr) / (sigma[1].angle - sigma[0].angle) * (angle - sigma[0].angle) + sigma[0].dr;

		}


		itr->s_dal = itr->dal / sigma_dal;
		itr->s_dar = itr->dar / sigma_dar;
		itr->s_dl = itr->dl / sigma_dl;
		itr->s_dr = itr->dr / sigma_dr;
	}


}

void Calc_difference(Difference& diff, std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack>& pair) {
	diff.ax = pair.first.ax;
	diff.ay = pair.first.ay;
	diff.dax = pair.second.ax - pair.first.ax;
	diff.day = pair.second.ay - pair.first.ay;

	diff.dx = pair.second.x - pair.first.x - (pair.first.ax + pair.second.ax) / 2 * (pair.second.z - pair.first.z);
	diff.dy = pair.second.y - pair.first.y - (pair.first.ay + pair.second.ay) / 2 * (pair.second.z - pair.first.z);
	diff.dar = (diff.dax * pair.first.ax + diff.day * pair.first.ay) / sqrt(pair.first.ax * pair.first.ax + pair.first.ay * pair.first.ay);
	diff.dal = (diff.dax * pair.first.ay - diff.day * pair.first.ax) / sqrt(pair.first.ax * pair.first.ax + pair.first.ay * pair.first.ay);

	//dr,dlの計算
	using namespace matrix_3D;
	vector_3D pos0, pos1, dir0, dir1;
	pos0.x = pair.first.x;
	pos0.y = pair.first.y;
	pos0.z = pair.first.z;
	dir0.x = pair.first.ax;
	dir0.y = pair.first.ay;
	dir0.z = 1;
	pos1.x = pair.second.x;
	pos1.y = pair.second.y;
	pos1.z = pair.second.z;
	dir1.x = pair.second.ax;
	dir1.y = pair.second.ay;
	dir1.z = 1;

	vector_3D base_point, difference;
	//外挿基準点を1:1に内分した点に設定
	base_point = addition(const_multiple(pos0, 0.5), const_multiple(pos1, 0.5));
	difference = addition(const_multiple(pos0, -1), pos1);

	vector_3D extra0, extra1;
	double ratio0, ratio1;
	ratio0 = -1 * dot(addition(pos0, const_multiple(base_point, -1)), difference) / dot(dir0, difference);
	ratio1 = -1 * dot(addition(pos1, const_multiple(base_point, -1)), difference) / dot(dir1, difference);
	extra0 = addition(pos0, const_multiple(dir0, ratio0));
	extra1 = addition(pos1, const_multiple(dir1, ratio1));

	vector_3D unit_r, unit_l;
	unit_l.x = -1 + difference.y;
	unit_l.y = difference.x;
	unit_l.z = 0;
	unit_r.x = -1 * difference.x * difference.z;
	unit_r.y = -1 * difference.y * difference.z;
	unit_r.z = pow(difference.x, 2) + pow(difference.y, 2);

	double constant;
	constant = sqrt(pow(difference.x, 2) + pow(difference.y, 2));
	unit_l.x = unit_l.x / constant;
	unit_l.y = unit_l.y / constant;
	constant = sqrt((pow(difference.x, 2) + pow(difference.y, 2)) * (pow(difference.x, 2) + pow(difference.y, 2) + pow(difference.z, 2)));
	unit_r.x = unit_r.x / constant;
	unit_r.y = unit_r.y / constant;
	unit_r.z = unit_r.z / constant;

	diff.dr = dot(addition(extra1, const_multiple(extra0, -1)), unit_r);
	diff.dl = dot(addition(extra1, const_multiple(extra0, -1)), unit_l);

	if (sqrt(pair.first.ax * pair.first.ax + pair.first.ay * pair.first.ay) < 0.001) {
		diff.dar = diff.day;
		diff.dal = diff.dax;
		diff.dr = diff.dy;
		diff.dl = diff.dx;
	}

}
void output_difference(std::string filename, std::vector<Difference>& diff) {

	std::ofstream ofs(filename);
	for (auto itr = diff.begin(); itr != diff.end(); itr++) {
		ofs << std::right << std::fixed
			<< std::setw(8) << std::setprecision(0) << itr->groupid << " "
			<< std::setw(8) << std::setprecision(0) << itr->chainid << " "
			<< std::setw(2) << std::setprecision(0) << itr->direction << " "
			<< std::setw(4) << std::setprecision(0) << itr->vertex_pl << " "
			<< std::setw(2) << std::setprecision(0) << itr->edge_out_flg << " "
			<< std::setw(4) << std::setprecision(0) << itr->target_pl[0] << " "
			<< std::setw(4) << std::setprecision(0) << itr->target_pl[1] << " "
			<< std::setw(4) << std::setprecision(0) << itr->edge_pl[0] << " "
			<< std::setw(4) << std::setprecision(0) << itr->edge_pl[1] << " "
			<< std::setw(4) << std::setprecision(0) << itr->add_nseg << " "
			<< std::setw(7) << std::setprecision(4) << itr->ax << " "
			<< std::setw(7) << std::setprecision(4) << itr->ay << " "
			<< std::setw(7) << std::setprecision(4) << itr->dax << " "
			<< std::setw(7) << std::setprecision(4) << itr->day << " "
			<< std::setw(8) << std::setprecision(1) << itr->dx << " "
			<< std::setw(8) << std::setprecision(1) << itr->dy << " "
			<< std::setw(7) << std::setprecision(4) << itr->dar << " "
			<< std::setw(7) << std::setprecision(4) << itr->dal << " "
			<< std::setw(8) << std::setprecision(1) << itr->dr << " "
			<< std::setw(8) << std::setprecision(1) << itr->dl << " "
			<< std::setw(7) << std::setprecision(4) << itr->s_dar << " "
			<< std::setw(7) << std::setprecision(4) << itr->s_dal << " "
			<< std::setw(7) << std::setprecision(4) << itr->s_dr << " "
			<< std::setw(7) << std::setprecision(4) << itr->s_dl << std::endl;
	}
}

std::set<std::pair<int, int>> Get_penetrate_list(std::vector<Difference>& diff) {
	std::set<std::pair<int, int>> ret;
	for (auto itr = diff.begin(); itr != diff.end(); itr++) {
		double chi2 = pow(itr->s_dal, 2) + pow(itr->s_dar, 2) + pow(itr->s_dl, 2) + pow(itr->s_dr, 2);
		if (chi2 > 15)continue;

		ret.insert(std::make_pair(itr->groupid, itr->chainid));
	}
	return ret;
}

std::vector<Momentum_recon::Event_information > cut_penetrate(std::vector<Momentum_recon::Event_information >& momch_ev, std::set<std::pair<int, int>>& id_list,std::string in) {
	std::vector<Momentum_recon::Event_information > ret;
	std::map<double, std::pair<double, double>> vph_mip = vph_mip_distribution(in);

	for (auto& ev : momch_ev) {
		Momentum_recon::Event_information chains = ev;
		//chains = ev;
		chains.chains.clear();
		int direction, vertex_pl;

		for (auto& c : ev.chains) {
			if (c.chainid == 0) {
				chains.chains.push_back(c);
				continue;
			}
			vertex_pl = ev.vertex_pl;
			direction = c.direction;
			Momentum_recon::Mom_chain chain = c;
			chain.base.clear();
			chain.base_pair.clear();
			//chainの分断
			//forward
			if (direction == 1) {
				for (auto& b : c.base) {
					if (b.pl <= vertex_pl) {
						chain.base.push_back(b);
					}
				}
				for (auto& pair : c.base_pair) {
					if (pair.second.pl <= vertex_pl) {
						chain.base_pair.push_back(pair);
					}
				}
			}
			else if (direction == -1) {
				for (auto& b : c.base) {
					if (b.pl > vertex_pl) {
						chain.base.push_back(b);
					}
				}
				for (auto& pair : c.base_pair) {
					if (pair.first.pl > vertex_pl) {
						chain.base_pair.push_back(pair);
					}
				}
			}

			//penetrate candidate
			if (id_list.count(std::make_pair(ev.groupid, chain.chainid)) == 1) {
				//blackの判定
				if (!judege_Black(chain, vph_mip, 5))continue;
			}
			chains.chains.push_back(chain);
		}
		ret.push_back(chains);

	}
	return ret;
}
bool judege_Black(Momentum_recon::Mom_chain& chain, std::map<double, std::pair<double, double>>& vph_mip, double thr_sigma) {
	double ax = 0, ay = 0, vph = 0;
	int count = 0, count2 = 0;
	for (auto itr = chain.base.begin(); itr != chain.base.end(); itr++) {
		ax += itr->ax;
		ay += itr->ay;
		count++;
		vph += itr->m[0].ph % 10000;
		count2++;
		vph += itr->m[1].ph % 10000;
		count2++;
	}
	ax /= count;
	ay /= count;
	vph /= count2;

	double angle = sqrt(ax * ax + ay * ay);

	double sigma = 0;
	auto vph_angle = vph_mip.upper_bound(angle);
	if (vph_angle == vph_mip.begin()) {
		sigma = (vph - vph_angle->second.first) / vph_angle->second.second;
	}
	else if (vph_angle == vph_mip.end()) {
		sigma = (vph - vph_mip.rbegin()->second.first) / vph_mip.rbegin()->second.second;
	}
	else {
		double vph_mean, vph_sigma;
		auto val0 = std::next(vph_angle, -1);
		auto val1 = vph_angle;
		vph_mean = (val1->second.first - val0->second.first) / (val1->first - val0->first) * (angle - val0->first) + val0->second.first;
		vph_sigma = (val1->second.second - val0->second.second) / (val1->first - val0->first) * (angle - val0->first) + val0->second.second;
		sigma = (vph - vph_mean) / vph_sigma;
	}
	if (sigma > thr_sigma)return true;
	return false;
}
std::map<double, std::pair<double, double>> vph_mip_distribution(std::string in) {
	std::map<double, std::pair<double, double>> ret;

	std::ifstream ifs(in);
	if (!ifs) {
		std::cerr << "Failed to open " << in << std::endl;
		exit(0);
	}

	double p[3];
	while (ifs >> p[0] >> p[1] >> p[2]) {
		ret.insert(std::make_pair(p[0], std::make_pair(p[1], p[2])));
	}

	for (auto itr = ret.begin(); itr != ret.end(); itr++) {
		std::cout << itr->first << " " << itr->second.first << " " << itr->second.second << std::endl;
	}

	return ret;
}

int judge_momch_edgeout(int direction, Momentum_recon::Mom_chain& chain, std::map<int, std::vector<Fiducial_Area::Fiducial_Area>>& area, double edge_cut, int ex_pl_max) {

	int edge_pl;
	int return_flg = 0;
	if (direction == 1) {
		edge_pl = chain.base.begin()->pl;
		return_flg = 0;
		for (int ex_pl = 0; ex_pl <= ex_pl_max; ex_pl++) {
			if (area.count(edge_pl - ex_pl) == 0)continue;
			//ex_z = z_map.at(up_pl + ex_pl);
			//ex_x = up_x + up_ax * (ex_z - up_z);
			//ex_y = up_y + up_ay * (ex_z - up_z);

			if (!judge_fiducial_area(area.at(edge_pl - ex_pl), *chain.base.begin())) {
				return_flg = 2;
			}
		}
		if (edge_pl >= 132)return_flg = 1;
	}
	else if (direction == -1) {
		edge_pl = chain.base.rbegin()->pl;
		return_flg = 0;
		for (int ex_pl = 0; ex_pl <= ex_pl_max; ex_pl++) {
			if (area.count(edge_pl + ex_pl) == 0)continue;
			//ex_z = z_map.at(up_pl + ex_pl);
			//ex_x = up_x + up_ax * (ex_z - up_z);
			//ex_y = up_y + up_ay * (ex_z - up_z);

			if (!judge_fiducial_area(area.at(edge_pl + ex_pl), *chain.base.rbegin())) {
				return_flg = 2;
			}
		}
		if (edge_pl <= 4)return_flg = 1;

	}

	return return_flg;

}

bool judge_fiducial_area(std::vector<Fiducial_Area::Fiducial_Area>& area, Momentum_recon::Mom_basetrack& b) {

	std::map<double, Fiducial_Area::Point> point_map;
	double ex_x, ex_y, dist;
	for (auto itr = area.begin(); itr != area.end(); itr++) {
		ex_x = b.x + b.ax * (itr->p[0].z - b.z);
		ex_y = b.y + b.ay * (itr->p[0].z - b.z);
		dist = pow(ex_x - itr->p[0].x, 2) + pow(ex_y - itr->p[0].y, 2);
		point_map.insert(std::make_pair(dist, itr->p[0]));
	}
	//外挿先から距離の一番近い点のz座標を使用
	double z = point_map.begin()->second.z;
	double x = b.x + b.ax * (z - b.z);
	double y = b.y + b.ay * (z - b.z);


	//true でArea内　falseでarea外

	//点(x,y)からx軸性の方向に直線を引き、その直線と多角形の辺が何回交わるか。
	//下から上に交わったときwn+1
	//上から下に交わったときwn-1
	int wn = 0;
	double vt;
	for (auto itr = area.begin(); itr != area.end(); itr++) {
		// 上向きの辺、下向きの辺によって処理が分かれる。
	// 上向きの辺。点Pがy軸方向について、始点と終点の間にある。ただし、終点は含まない。(ルール1)
		if (itr->p[0].y <= y && itr->p[1].y > y) {
			// 辺は点pよりも右側にある。ただし、重ならない。(ルール4)
			// 辺が点pと同じ高さになる位置を特定し、その時のxの値と点pのxの値を比較する。
			vt = (y - itr->p[0].y) / (itr->p[1].y - itr->p[0].y);
			if (x < itr->p[0].x + vt * (itr->p[1].x - itr->p[0].x)) {
				++wn;  //ここが重要。上向きの辺と交差した場合は+1
			}
		}
		// 下向きの辺。点Pがy軸方向について、始点と終点の間にある。ただし、始点は含まない。(ルール2)
		else if (itr->p[0].y > y && itr->p[1].y <= y) {
			// 辺は点pよりも右側にある。ただし、重ならない。(ルール4)
			// 辺が点pと同じ高さになる位置を特定し、その時のxの値と点pのxの値を比較する。
			vt = (y - itr->p[0].y) / (itr->p[1].y - itr->p[0].y);
			if (x < itr->p[0].x + vt * (itr->p[1].x - itr->p[0].x)) {
				--wn;  //ここが重要。下向きの辺と交差した場合は-1
			}
		}
	}
	if (wn >= 1)return true;
	return false;
}

//mfile0::M_Base
//basetrack-alignment mapの対応
double select_triangle_vale(corrmap_3d::align_param2* param, Fiducial_Area::Point* p) {
	double x, y;
	double dist = 0;
	x = (param->corr_p[0]->x + param->corr_p[1]->x + param->corr_p[2]->x) / 3;
	y = (param->corr_p[0]->y + param->corr_p[1]->y + param->corr_p[2]->y) / 3;
	dist = (p->x - x) * (p->x - x) + (p->y - y) * (p->y - y);
	return dist;
}
corrmap_3d::align_param2* search_param(std::vector<corrmap_3d::align_param*>& param, Fiducial_Area::Point* p, std::multimap<int, corrmap_3d::align_param2*>& triangles) {
	//三角形内部
	//最近接三角形
	double dist = 0;
	std::map<double, corrmap_3d::align_param* > dist_map;
	//align_paramを近い順にsort
	for (auto itr = param.begin(); itr != param.end(); itr++) {
		dist = ((*itr)->x - p->x) * ((*itr)->x - p->x) + ((*itr)->y - p->y) * ((*itr)->y - p->y);
		dist_map.insert(std::make_pair(dist, (*itr)));
	}

	double sign[3];
	bool flg = false;
	int id;

	corrmap_3d::align_param2* ret = triangles.begin()->second;
	for (auto itr = dist_map.begin(); itr != dist_map.end(); itr++) {
		if (itr != dist_map.begin())continue;


		//corrmapのID
		id = itr->second->id;
		if (triangles.count(id) == 0) {
			fprintf(stderr, "alignment triangle ID=%d not found\n", id);
			exit(1);
		}
		//idの属する三角形を探索
		auto range = triangles.equal_range(id);
		for (auto itr2 = range.first; itr2 != range.second; itr2++) {
			sign[0] = (itr2->second->corr_p[1]->x - itr2->second->corr_p[0]->x) * (p->y - itr2->second->corr_p[1]->y) - (itr2->second->corr_p[1]->y - itr2->second->corr_p[0]->y) * (p->x - itr2->second->corr_p[1]->x);
			sign[1] = (itr2->second->corr_p[2]->x - itr2->second->corr_p[1]->x) * (p->y - itr2->second->corr_p[2]->y) - (itr2->second->corr_p[2]->y - itr2->second->corr_p[1]->y) * (p->x - itr2->second->corr_p[2]->x);
			sign[2] = (itr2->second->corr_p[0]->x - itr2->second->corr_p[2]->x) * (p->y - itr2->second->corr_p[0]->y) - (itr2->second->corr_p[0]->y - itr2->second->corr_p[2]->y) * (p->x - itr2->second->corr_p[0]->x);
			//printf("point %.lf,%.1lf\n", base.x, base.y);
			//printf("triangle %.1lf %.1lf\n", itr2->second->corr_p[0]->x, itr2->second->corr_p[0]->y);
			//printf("triangle %.1lf %.1lf\n", itr2->second->corr_p[1]->x, itr2->second->corr_p[1]->y);
			//printf("triangle %.1lf %.1lf\n", itr2->second->corr_p[2]->x, itr2->second->corr_p[2]->y);
			//printf("sign %.1lf %1.lf %.1lf\n", sign[0], sign[1], sign[2]);
			//printf("  signbit %d %d %d\n", std::signbit(sign[0]), std::signbit(sign[1]), std::signbit(sign[2]));
			//printf("n signbit %d %d %d\n", !std::signbit(sign[0]), !std::signbit(sign[1]), !std::signbit(sign[2]));
			//printf("judge %d\n", (std::signbit(sign[0]) && std::signbit(sign[1]) && std::signbit(sign[2])) || (!std::signbit(sign[0]) && !std::signbit(sign[1]) && !std::signbit(sign[2])));
			//printf("\n");

			//符号が3つとも一致でtrue
			if ((std::signbit(sign[0]) && std::signbit(sign[1]) && std::signbit(sign[2])) || (!std::signbit(sign[0]) && !std::signbit(sign[1]) && !std::signbit(sign[2]))) {
				ret = itr2->second;
				flg = true;
				break;
			}
		}
		if (flg)break;
	}
	if (flg) {
		//printf("point in trianlge\n");
		return ret;
	}

	//distが最小になるcorrmapをとってくる
	dist = -1;
	for (auto itr = dist_map.begin(); itr != dist_map.end(); itr++) {
		//corrmapのID
		id = itr->second->id;
		if (triangles.count(id) == 0) {
			fprintf(stderr, "alignment triangle ID=%d not found\n", id);
			exit(1);
		}
		//idの属する三角形を探索
		auto range = triangles.equal_range(id);
		for (auto itr2 = range.first; itr2 != range.second; itr2++) {
			if (dist<0 || dist>Fiducial_Area::select_triangle_vale(itr2->second, p)) {
				dist = Fiducial_Area::select_triangle_vale(itr2->second, p);
				ret = itr2->second;
			}
		}
	}
	//printf("point not in trianlge\n");
	return ret;
}
std::vector <std::pair<Fiducial_Area::Point*, corrmap_3d::align_param2*>>track_affineparam_correspondence(std::vector<Fiducial_Area::Point*>& p, std::vector <corrmap_3d::align_param2>& param) {

	//local alignの視野中心を取り出して、位置でhash
	//local alignの視野中心の作るdelaunay三角形をmapで対応

	std::map<int, corrmap_3d::align_param*> view_center;
	std::multimap<int, corrmap_3d::align_param2*>triangles;
	double xmin = 999999, ymin = 999999, hash = 2000;
	for (auto itr = param.begin(); itr != param.end(); itr++) {
		for (int i = 0; i < 3; i++) {
			view_center.insert(std::make_pair(itr->corr_p[i]->id, (itr->corr_p[i])));
			triangles.insert(std::make_pair(itr->corr_p[i]->id, &(*itr)));
			xmin = std::min(itr->corr_p[i]->x, xmin);
			ymin = std::min(itr->corr_p[i]->y, ymin);
		}
	}
	std::multimap<std::pair<int, int>, corrmap_3d::align_param*> view_center_hash;
	std::pair<int, int>id;
	for (auto itr = view_center.begin(); itr != view_center.end(); itr++) {
		id.first = int((itr->second->x - xmin) / hash);
		id.second = int((itr->second->y - ymin) / hash);
		view_center_hash.insert(std::make_pair(id, itr->second));
	}

	std::vector < std::pair<Fiducial_Area::Point*, corrmap_3d::align_param2*>> ret;
	std::vector<corrmap_3d::align_param*> param_cand;
	int loop = 0, ix, iy, count = 0;
	for (auto itr = p.begin(); itr != p.end(); itr++) {
		if (count % 100000 == 0) {
			printf("\r search correspond triangles %d/%d(%4.1lf%%)", count, p.size(), count * 100. / p.size());
		}
		count++;
		ix = ((*itr)->x - xmin) / hash;
		iy = ((*itr)->y - ymin) / hash;
		loop = 1;
		while (true) {
			param_cand.clear();
			for (int iix = ix - loop; iix <= ix + loop; iix++) {
				for (int iiy = iy - loop; iiy <= iy + loop; iiy++) {
					id.first = iix;
					id.second = iiy;
					if (view_center_hash.count(id) != 0) {
						auto range = view_center_hash.equal_range(id);
						for (auto res = range.first; res != range.second; res++) {
							param_cand.push_back(res->second);
						}
					}
				}
			}
			if (param_cand.size() > 2)break;
			loop++;
		}
		corrmap_3d::align_param2* param2 = Fiducial_Area::search_param(param_cand, *itr, triangles);
		ret.push_back(std::make_pair((*itr), param2));
	}
	printf("\r search correspond triangles %d/%d(%4.1lf%%)\n", count, p.size(), count * 100. / p.size());

	return ret;
}
//変換 zshrink補正-->9para変換
void trans_base(std::vector<Fiducial_Area::Point*>& p, corrmap_3d::align_param2* param) {

	matrix_3D::matrix_33 x_rot_mat(0, param->x_rot), y_rot_mat(1, param->y_rot), z_rot_mat(2, param->z_rot), all_trans(0, 0), shear_mat(0, 0), shrink_mat(0, 0);

	shrink_mat.val[0][0] *= param->x_shrink;
	shrink_mat.val[1][1] *= param->y_shrink;
	//shrink_mat.val[2][2] *= param->z_shrink;
	shear_mat.val[0][1] = param->yx_shear;
	//shear_mat.val[0][2] = param->zx_shear;
	//shear_mat.val[1][2] = param->zy_shear;

	matrix_3D::vector_3D shift, center;
	center.x = param->x;
	center.y = param->y;
	center.z = param->z;
	shift.x = param->dx;
	shift.y = param->dy;
	shift.z = param->dz;

	all_trans.matrix_multiplication(shear_mat);
	all_trans.matrix_multiplication(shrink_mat);
	all_trans.matrix_multiplication(z_rot_mat);
	all_trans.matrix_multiplication(y_rot_mat);
	all_trans.matrix_multiplication(x_rot_mat);

	//all_trans.Print();
	matrix_3D::vector_3D base_p0;
	double base_thick = 210;
	for (auto itr = p.begin(); itr != p.end(); itr++) {
		base_p0.x = (*itr)->x;
		base_p0.y = (*itr)->y;
		base_p0.z = param->z;

		//base_p1.x = (*itr)->x + (*itr)->ax*base_thick;
		//base_p1.y = (*itr)->y + (*itr)->ay*base_thick;
		////角度shrink分はここでかける
		//base_p1.z = param->z + base_thick / param->z_shrink;

		//視野中心を原点に移動
		//base_p0 = matrix_3D::addition(base_p0, matrix_3D::const_multiple(center, -1));
		//base_p1 = matrix_3D::addition(base_p1, matrix_3D::const_multiple(center, -1));

		//変換の実行
		base_p0.matrix_multiplication(all_trans);
		base_p0 = matrix_3D::addition(base_p0, shift);
		//base_p1.matrix_multiplication(all_trans);
		//base_p1 = matrix_3D::addition(base_p1, shift);

		//原点をもとに戻す
		//base_p0 = matrix_3D::addition(base_p0, center);
		//base_p1 = matrix_3D::addition(base_p1, center);

		(*itr)->x = base_p0.x;
		(*itr)->y = base_p0.y;
		(*itr)->z = base_p0.z;

		//printf("ax:%.4lf --> %.4lf\n", (*itr)->ax, (base_p0.x - base_p1.x) / (base_p0.z - base_p1.z));
		//printf("ay:%.4lf --> %.4lf\n", (*itr)->ay, (base_p0.y - base_p1.y) / (base_p0.z - base_p1.z));

		//(*itr)->ax = (base_p0.x - base_p1.x) / (base_p0.z - base_p1.z) + param->zx_shear;
		//(*itr)->ay = (base_p0.y - base_p1.y) / (base_p0.z - base_p1.z) + param->zy_shear;

	}
}
void trans_base_all(std::vector < std::pair<Fiducial_Area::Point*, corrmap_3d::align_param2*>>& track_pair) {
	std::map<std::tuple<int, int, int>, corrmap_3d::align_param2*> param_map;
	std::multimap<std::tuple<int, int, int>, Fiducial_Area::Point*>base_map;
	std::tuple<int, int, int>id;
	//三角形ごとにbasetrackをまとめる
	for (auto itr = track_pair.begin(); itr != track_pair.end(); itr++) {
		std::get<0>(id) = itr->second->corr_p[0]->id;
		std::get<1>(id) = itr->second->corr_p[1]->id;
		std::get<2>(id) = itr->second->corr_p[2]->id;
		param_map.insert(std::make_pair(id, itr->second));
		base_map.insert(std::make_pair(id, itr->first));
	}


	//ここで三角形ごとに変換
	int count = 0;
	std::vector<Fiducial_Area::Point*> t_base;
	for (auto itr = param_map.begin(); itr != param_map.end(); itr++) {
		if (count % 1000 == 0) {
			printf("\r basetrack trans num of triangles %d/%d(%4.1lf%%)", count, param_map.size(), count * 100. / param_map.size());
		}
		count++;

		t_base.clear();

		if (base_map.count(itr->first) == 0)continue;
		auto range = base_map.equal_range(itr->first);
		for (auto res = range.first; res != range.second; res++) {
			t_base.push_back(res->second);
		}
		Fiducial_Area::trans_base(t_base, itr->second);

	}
	printf("\r basetrack trans num of triangles %d/%d(%4.1lf%%)\n", count, param_map.size(), count * 100. / param_map.size());

}

