// 2025/03/11
// 2025/03/14 inutにmio_vph_distribusionを指定するよう変更
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>
//std::vector<Difference> Calc_difference

class Difference {
public:
	int groupid, chainid, vertex_pl, target_pl[2], direction, add_nseg, edge_pl[2], edge_out_flg;
	double ax, ay, dax, day, dx, dy, dar, dal, dr, dl;
	double  s_dar, s_dal, s_dr, s_dl, s_pb_angle, s_pb_position;
	double mcs_pb[2], mcs_pb_error[2][2];
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
std::vector<Difference> Calc_difference(std::vector<Momentum_recon::Event_information>& momch_ev, std::map<int, std::map<int, Momentum_recon::Mom_chain>>& momch_mom_map, std::map<int, std::vector<Fiducial_Area::Fiducial_Area>>& area);
void Calc_difference(Difference& diff, std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack>& pair);


void output_difference(std::string filename, std::vector<Difference>& diff);
std::map<int, std::map<double, Sigma_list>> Read_sigma_list(std::string filename);
void Calc_sigma(std::vector<Difference>& diff, std::map<int, std::map<double, Sigma_list>>& sigma_list);
std::set<std::pair<int, int>> Get_penetrate_list(std::vector<Difference>& diff);
std::map<double, std::pair<double, double>> vph_mip_distribution(std::string in);
bool judege_Black(Momentum_recon::Mom_chain& chain, std::map<double, std::pair<double, double>>& vph_mip, double thr_sigma);
std::vector<Momentum_recon::Event_information > cut_penetrate(std::vector<Momentum_recon::Event_information >& momch_ev, std::set<std::pair<int, int>>& id_list);
std::vector<Momentum_recon::Event_information > cut_penetrate(std::vector<Momentum_recon::Event_information >& momch_ev, std::map<int, std::map<int, Momentum_recon::Mom_chain>>& momch_mom_map, std::set<std::pair<int, int>>& id_list, std::string in);

int judge_momch_edgeout(int direction, Momentum_recon::Mom_chain& chain, std::map<int, std::vector<Fiducial_Area::Fiducial_Area>>& area, double edge_cut, int ex_pl_max);
std::map<int, std::map<int, Momentum_recon::Mom_chain>> momch_format_change(std::vector<Momentum_recon::Event_information>& momch_mom);
void Calc_sigma_pb(Difference& diff, int index, std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack>& pair);
bool judge_fiducial_area(std::vector<Fiducial_Area::Fiducial_Area>& area, Momentum_recon::Mom_basetrack& b);
int count_pb_angle_diff(Momentum_recon::Mom_chain& momch);

int main(int argc, char** argv) {
	if (argc != 9) {
		printf("argc %d\n", argc);
		fprintf(stderr, "usage:file-in-momch-all file-in-momch-mom sigma_list file-in-ECC fa.txt mip_vph_distribution.txt file-out file_out_momch\n");
		exit(1);
	}
	std::string file_in_momch_all = argv[1];
	std::string file_in_momch_mom = argv[2];
	std::string file_in_sigma = argv[3];
	std::string file_in_ECC = argv[4];
	std::string file_in_area = argv[5];
	std::string in = argv[6];

	std::string file_out_difference = argv[7];
	std::string file_out_momch = argv[8];

	//corrmap absの読み込み
	std::string file_in_corrmap = file_in_ECC + "\\Area0\\0\\align\\fine\\local\\corrmap-local-abs.lst";
	std::map<int, std::vector<corrmap_3d::align_param>>corrmap = corrmap_3d::read_ali_param_abs(file_in_corrmap, 1);
	std::map<int, std::vector<corrmap_3d::align_param2>>corrmap_dd = corrmap_3d::DelaunayDivide_map(corrmap);

	//fiducial areaの読み込み
	std::map<int, std::vector<Fiducial_Area::Fiducial_Area>> area = Fiducial_Area::read_fiducial_Area(file_in_area);
	for (auto itr = area.begin(); itr != area.end(); itr++) {
		if (corrmap_dd.count(itr->first) == 0) {
			fprintf(stderr, "corrmap local abs PL%03d not found\n", itr->first);
			exit(1);
		}
		std::vector<corrmap_3d::align_param2> param = corrmap_dd.at(itr->first);
		trans_mfile_cordinate(param, itr->second);
	}


	std::map<int, std::map<double, Sigma_list>> sigma_list = Read_sigma_list(file_in_sigma);

	std::vector<Momentum_recon::Event_information> momch_all = Momentum_recon::Read_Event_information_extension(file_in_momch_all);
	std::vector<Momentum_recon::Event_information> momch_mom = Momentum_recon::Read_Event_information_extension(file_in_momch_mom);

	std::map<int, std::map<int, Momentum_recon::Mom_chain>> momch_mom_map = momch_format_change(momch_mom);


	std::vector<Difference>diff = Calc_difference(momch_all, momch_mom_map, area);
	Calc_sigma(diff, sigma_list);
	output_difference(file_out_difference, diff);
	std::set<std::pair<int, int>> penetrate_id = Get_penetrate_list(diff);
	momch_all = cut_penetrate(momch_all, momch_mom_map, penetrate_id,in);
	Momentum_recon::Write_Event_information_extension(file_out_momch, momch_all);
}


std::map<int, std::map<int, Momentum_recon::Mom_chain>> momch_format_change(std::vector<Momentum_recon::Event_information>& momch_mom) {

	std::map<int, std::map<int, Momentum_recon::Mom_chain>> ret;
	for (auto& ev : momch_mom) {
		std::map<int, Momentum_recon::Mom_chain> chain_event;
		for (auto& ch : ev.chains) {
			chain_event.insert(std::make_pair(ch.chainid, ch));
		}
		ret.insert(std::make_pair(ev.groupid, chain_event));
	}
	return ret;


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

std::vector<Difference> Calc_difference(std::vector<Momentum_recon::Event_information>& momch_ev, std::map<int, std::map<int, Momentum_recon::Mom_chain>>& momch_mom_map, std::map<int, std::vector<Fiducial_Area::Fiducial_Area>>& area) {
	std::vector<Difference> ret;

	for (auto& ev : momch_ev) {
		Difference diff;
		diff.groupid = ev.groupid;
		//for ECC1 Water
		//if (ev.groupid == 6302 || ev.groupid == 7155)continue;
		//iron
		//if (ev.groupid == 6201)continue;
		//for ECC6 Water
		//if (ev.groupid == 4429)continue;

		std::map<int, Momentum_recon::Mom_chain> mom_map, mom_map_inv;
		if (momch_mom_map.count(diff.groupid) == 0) {
			fprintf(stderr, "event %d partner not found\n", diff.groupid);
			continue;//add 20251011
		}
		if (momch_mom_map.count(diff.groupid * -1) == 0) {
			fprintf(stderr, "event %d partner inv not found\n", diff.groupid);
			continue;//add 20251011 mapの処理で落ちるので
		}
		mom_map = momch_mom_map.at(diff.groupid);
		mom_map_inv = momch_mom_map.at(-1 * diff.groupid);
		for (auto& c : ev.chains) {
			if (c.chainid == 0)continue;
			diff.vertex_pl = ev.vertex_pl;
			diff.direction = c.direction;
			diff.chainid = c.chainid;
			diff.edge_pl[0] = c.base.begin()->pl;
			diff.edge_pl[1] = c.base.rbegin()->pl;
			diff.edge_out_flg = judge_momch_edgeout(diff.direction, c, area, 0, 4);
			diff.add_nseg = 0;
			diff.target_pl[0] = -1;
			diff.target_pl[1] = -1;
			diff.s_pb_angle = 0;
			diff.s_pb_position = 0;
			diff.mcs_pb[0] = -1;
			diff.mcs_pb[1] = -1;
			diff.mcs_pb_error[0][0] = -1;
			diff.mcs_pb_error[0][1] = -1;
			diff.mcs_pb_error[1][0] = -1;
			diff.mcs_pb_error[1][1] = -1;
			std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack> divide_pair;
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
						divide_pair = pair;
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
						divide_pair = pair;
						break;
					}
				}

			}

			double pb_error[2];
			Momentum_recon::Mom_chain mom, mom_inv;
			if (mom_map.count(c.chainid) == 0) {
				fprintf(stderr, "event %d chain %d partner not found\n", diff.groupid, c.chainid);
			}
			else {
				mom = mom_map.at(c.chainid);
				diff.mcs_pb[0] = mom.Get_muon_mcs_pb();
				mom.Get_muon_pb_mcs_error(pb_error);
				diff.mcs_pb_error[0][0] = pb_error[0];
				diff.mcs_pb_error[0][1] = pb_error[1];
			}

			if (mom_map_inv.count(c.chainid) == 0) {
				fprintf(stderr, "event %d chain %d partner inv not found\n", diff.groupid, c.chainid);
				continue;
			}
			else {
				mom_inv = mom_map_inv.at(c.chainid);
				diff.mcs_pb[1] = mom_inv.Get_muon_mcs_pb();
				mom_inv.Get_muon_pb_mcs_error(pb_error);
				diff.mcs_pb_error[1][0] = pb_error[0];
				diff.mcs_pb_error[1][1] = pb_error[1];
			}
			if (count_pb_angle_diff(mom) < 2 && count_pb_angle_diff(mom_inv) < 2) {
				Calc_sigma_pb(diff, 2, divide_pair);
			}
			else if (count_pb_angle_diff(mom) < count_pb_angle_diff(mom_inv)) {
				Calc_sigma_pb(diff, 1, divide_pair);
			}
			else {
				Calc_sigma_pb(diff, 0, divide_pair);
			}
			ret.push_back(diff);
		}
	}
	return ret;

}
int count_pb_angle_diff(Momentum_recon::Mom_chain& momch) {
	int ret = 0;
	for (auto& pair : momch.base_pair) {
		if (pair.second.pl - pair.first.pl != 1)continue;
		if (pair.first.pl <= 3)continue;
		else if (pair.first.pl <= 14)ret++;
		else if (pair.first.pl <= 15)continue;
		else if (pair.first.pl % 2 == 0)ret++;

	}
	return ret;
}
void Calc_sigma_pb(Difference& diff, int index, std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack>& pair) {
	double pb;
	if (index >= 2) {
		diff.s_pb_angle = 0;
		diff.s_pb_position = 0;
		return;
	}
	if (diff.mcs_pb[index] > 0) {
		pb = diff.mcs_pb[index];
	}
	else {
		diff.s_pb_angle = 0;
		diff.s_pb_position = 0;
		return;
	}
	double path = sqrt(pow(pair.first.x - pair.second.x, 2) + pow(pair.first.y - pair.second.y, 2) + pow(pair.first.z - pair.second.z, 2));
	double dz = fabs(pair.first.z - pair.second.z);
	int num_iron = 0, num_film = 0, num_water = 0;
	for (int pl = pair.first.pl; pl < pair.second.pl; pl++) {
		num_film++;
		if (pl == 3)continue;
		else if (pl <= 14)num_iron++;
		else if (pl == 15)continue;
		else if (pl % 2 == 0)num_iron++;
		else num_water++;
	}
	double thick_iron_ratio, thick_water_ratio, thick_base_ratio, thick_gel_ratio, thick_pack_ratio;
	thick_iron_ratio = 500 * num_iron / dz;
	thick_base_ratio = 210 * num_film / dz;
	thick_gel_ratio = 70 * num_film * 2 / dz;
	thick_pack_ratio = 100 * num_water * 2 / dz;
	thick_water_ratio = 1 - (thick_iron_ratio + thick_gel_ratio + thick_base_ratio + thick_pack_ratio);
	double radiation_length = 0;
	radiation_length += thick_iron_ratio * path / (17.57 * 1000);
	radiation_length += thick_gel_ratio * path / (30.3 * 1000);
	radiation_length += thick_base_ratio * path / (413.1 * 1000);
	radiation_length += thick_pack_ratio * path / (413.1 * 1000);
	radiation_length += thick_water_ratio * path / (360.8 * 1000);
	double mass = 105.65836668;
	double p2 = pb * pb / 2 * (1 + sqrt(1 + pow(2 * mass / pb, 2)));
	double beta = 1 / sqrt(1 + pow(mass / sqrt(p2), 2));

	diff.s_pb_angle = 13.6 / pb * sqrt(radiation_length) * (1 + 0.038 * log(radiation_length / beta));
	diff.s_pb_position = diff.s_pb_angle * path / (4 * sqrt(3));
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

		sigma_dal = sqrt(pow(sigma_dal, 2) + pow(itr->s_pb_angle, 2));
		sigma_dar = sqrt(pow(sigma_dar, 2) + pow(itr->s_pb_angle, 2));
		sigma_dl = sqrt(pow(sigma_dl, 2) + pow(itr->s_pb_position, 2));
		sigma_dr = sqrt(pow(sigma_dr, 2) + pow(itr->s_pb_position, 2));

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

	double denominator, d_theta_r, d_theta_l;
	denominator = (pair.first.ax * pair.second.ax + pair.first.ay * pair.second.ay + 1) / (sqrt(pair.first.ax * pair.first.ax + pair.first.ay * pair.first.ay + 1));
	d_theta_r = (-1 * pair.first.ax * pair.second.ax - pair.first.ay * pair.second.ay + pair.first.ax * pair.first.ax + pair.first.ay * pair.first.ay)
		/
		(sqrt(pair.first.ax * pair.first.ax + pair.first.ay * pair.first.ay + pow(pair.first.ax * pair.first.ax + pair.first.ay * pair.first.ay, 2)));
	d_theta_l = (-1 * pair.first.ay * pair.second.ax + pair.first.ax * pair.second.ay) / sqrt(pair.first.ax * pair.first.ax + pair.first.ay * pair.first.ay);

	diff.dar = atan(d_theta_r / denominator);
	diff.dal = atan(d_theta_l / denominator);


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
			<< std::setw(7) << std::setprecision(4) << itr->s_dl << " "
			<< std::setw(7) << std::setprecision(4) << itr->s_pb_angle << " "
			<< std::setw(7) << std::setprecision(1) << itr->s_pb_position << " "
			<< std::setw(7) << std::setprecision(1) << itr->mcs_pb[0] << " "
			<< std::setw(7) << std::setprecision(1) << itr->mcs_pb_error[0][0] << " "
			<< std::setw(7) << std::setprecision(1) << itr->mcs_pb_error[0][1] << " "
			<< std::setw(7) << std::setprecision(1) << itr->mcs_pb[1] << " "
			<< std::setw(7) << std::setprecision(1) << itr->mcs_pb_error[1][0] << " "
			<< std::setw(7) << std::setprecision(1) << itr->mcs_pb_error[1][1] << std::endl;
	}
}

std::set<std::pair<int, int>> Get_penetrate_list(std::vector<Difference>& diff) {
	std::set<std::pair<int, int>> ret;
	for (auto itr = diff.begin(); itr != diff.end(); itr++) {
		double chi2 = pow(itr->s_dal, 2) + pow(itr->s_dar, 2) + pow(itr->s_dl, 2) + pow(itr->s_dr, 2);
		if (chi2 > 30)continue;

		ret.insert(std::make_pair(itr->groupid, itr->chainid));
	}
	return ret;
}

std::vector<Momentum_recon::Event_information > cut_penetrate(std::vector<Momentum_recon::Event_information >& momch_ev, std::map<int, std::map<int, Momentum_recon::Mom_chain>>& momch_mom_map, std::set<std::pair<int, int>>& id_list, std::string in) {
	std::vector<Momentum_recon::Event_information > ret;
	std::map<double, std::pair<double, double>> vph_mip = vph_mip_distribution(in);

	for (auto& ev : momch_ev) {
		Momentum_recon::Event_information chains = ev;
		std::map<int, Momentum_recon::Mom_chain> mom_map;
		if (momch_mom_map.count(ev.groupid) == 0) {
			fprintf(stderr, "event %d partner not found\n", ev.groupid);
		}
		mom_map = momch_mom_map.at(ev.groupid);

		//chains = ev;
		chains.chains.clear();
		int direction, vertex_pl;

		for (auto& c : ev.chains) {
			Momentum_recon::Mom_chain mom;
			if (c.chainid == 0) {
				if (mom_map.count(c.chainid) == 0) {
					fprintf(stderr, "event %d chain %d partner not found\n", ev.groupid, c.chainid);
				}
				else {
					mom = mom_map.at(c.chainid);
				}
				chains.chains.push_back(mom);
				continue;
			}
			vertex_pl = ev.vertex_pl;
			direction = c.direction;
			if (mom_map.count(c.chainid) == 0) {
				fprintf(stderr, "event %d chain %d partner not found\n", ev.groupid, c.chainid);
			}
			else {
				mom = mom_map.at(c.chainid);
			}

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
				else {
					printf("%d %d\n", ev.groupid, chain.chainid);
				}
			}
			//attach側のみでcut
			//chains.chains.push_back(mom);
			//cut無し
			chains.chains.push_back(c);

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
		if (edge_pl <= 4)return_flg = 1;
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
		if (edge_pl >= 132)return_flg = 1;

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

