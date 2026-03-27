#pragma comment(lib,"FILE_structure.lib")
#pragma comment(lib,"functions.lib")
#include <FILE_structure.hpp>
#include <functions.hpp>
# include <iostream>
# include <random>

#include <list>

class stop_track {
public:
	int64_t chainid, groupid;
	int  nseg, npl, pl0, pl1, ph, rawid;
	double ax, ay, x, y, z;
	// ph-->pid
	int stoppl;
	int unixtime;

};
class track_pair {
public:
	int eventid;
	double x, y, z, md, oa;
	double dz;
	stop_track t[2];
};
class track_multi {
public:
	int eventid, pl;
	double x, y, z;
	std::vector< std::pair<double, stop_track>> trk;
	std::vector<track_pair>pair;
	int unixtime;
	double dz;
};
struct vtx_point {
	double x, y, z;
};

double minimum_distance_fixed(matrix_3D::vector_3D pos0, matrix_3D::vector_3D pos1, matrix_3D::vector_3D dir0, matrix_3D::vector_3D dir1, double z_range[2], double extra[2], double refz) {
	double extra0_distance, extra1_distance, delta;
	matrix_3D::vector_3D pos;
	pos.x = pos1.x - pos0.x;
	pos.y = pos1.y - pos0.y;
	pos.z = pos1.z - pos0.z;
	//ほぼ平行な場合
	if (opening_angle(dir0, dir1) < 0.0001) {
		extra0_distance = (pos1.z + pos0.z) / 2 - pos0.z;
		extra1_distance = (pos1.z + pos0.z) / 2 - pos1.z;
	}
	else {
		delta = dot(dir0, dir0) * dot(dir1, dir1) - pow(dot(dir0, dir1), 2.);
		extra0_distance = (+1 * dot(pos, dir0) * dot(dir1, dir1) - dot(dir0, dir1) * dot(pos, dir1)) / delta;
		extra1_distance = (-1 * dot(pos, dir1) * dot(dir0, dir0) + dot(dir0, dir1) * dot(pos, dir0)) / delta;
	}
	//range[0]:小,range[1]:大
	if (z_range[0] > z_range[1]) {
		double tmp_d = z_range[0];
		z_range[0] = z_range[1];
		z_range[1] = tmp_d;
	}

	matrix_3D::vector_3D extra0 = addition(pos0, const_multiple(dir0, extra0_distance));
	matrix_3D::vector_3D extra1 = addition(pos1, const_multiple(dir1, extra1_distance));

	if (extra0.z < refz + z_range[0] || extra1.z < refz + z_range[0]) {//2025/8/20 fixed
		extra0_distance = refz - pos0.z + z_range[0];
		extra1_distance = refz - pos1.z + z_range[0];
	}
	else if (extra0.z > refz + z_range[1] || extra1.z > refz + z_range[1]) {//2025/8/20 fixed
		extra0_distance = refz - pos0.z + z_range[1];
		extra1_distance = refz - pos1.z + z_range[1];
	}

	extra[0] = extra0_distance;
	extra[1] = extra1_distance;
	extra0 = addition(pos0, const_multiple(dir0, extra0_distance));
	extra1 = addition(pos1, const_multiple(dir1, extra1_distance));

	return distance(extra0, extra1);

}


void add_vtx_information(std::vector<Momentum_recon::Event_information>& momch, std::map<int, vtx_point>& vtx);
void read_stop_txt(std::vector<Momentum_recon::Event_information>& momch, std::multimap<int, stop_track>& tracks);
void clustering_2trk_vtx(std::multimap<int, stop_track>& tracks, int pl, std::map<int, vtx_point>& vtx);

int main(int argc, char** argv) {
	if (argc != 3) {
		fprintf(stderr, "usage:prg in-mfile.momch out-vtx.momch\n");
		exit(1);
	}
	std::string in_momch = argv[1];
	std::string out_momch = argv[2];

	//read momch
	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(in_momch);

	// reading stop.txt
	std::multimap<int, stop_track> stop;
	read_stop_txt(momch, stop);
	std::cout << "fin reading stoptrack." << std::endl;


	std::set<int> set;
	for (auto itr = stop.begin(); itr != stop.end(); itr++) {
		set.insert(itr->first);
	}
	std::cout << set.size() << std::endl;

	std::map<int, vtx_point> multi;
	for (auto ev = set.begin(); ev != set.end(); ev++) {
		auto tks = stop.equal_range(*ev);
		std::multimap<int, stop_track> rid;
		for (auto itr0 = tks.first; itr0 != tks.second; itr0++) {
			rid.insert(std::make_pair(itr0->second.rawid, itr0->second));
		}
		printf("event %5d, vtx search\n", *ev);
		clustering_2trk_vtx(rid, rid.begin()->second.stoppl,multi);

		rid.clear();
	}
		add_vtx_information(momch, multi);

	Momentum_recon::Write_Event_information_extension(out_momch, momch);

}
void read_stop_txt(std::vector<Momentum_recon::Event_information>& momch, std::multimap<int, stop_track>& tracks) {
	stop_track stop_tmp;

	int pos = 0;
	int cnt = 0;
	for (auto& ev : momch) {
		stop_tmp.stoppl = ev.vertex_pl;
		stop_tmp.groupid = ev.groupid;
		stop_tmp.unixtime = ev.unix_time;

		for (auto& c : ev.chains) {

			stop_tmp.chainid = c.chainid;
			stop_tmp.nseg = c.base.size();
			stop_tmp.pl0 = c.base.begin()->pl;//dounstream?
			stop_tmp.pl1 = c.base.rbegin()->pl;
			stop_tmp.npl = stop_tmp.pl1 - stop_tmp.pl0 + 1;
			stop_tmp.ph = c.particle_flg;
			if (stop_tmp.pl1 >= stop_tmp.stoppl) {
				//stop
				pos = stop_tmp.pl1;
				stop_tmp.rawid = c.base.rbegin()->rawid;
				stop_tmp.ax = c.base.rbegin()->ax;
				stop_tmp.ay = c.base.rbegin()->ay;
				stop_tmp.x = c.base.rbegin()->x;
				stop_tmp.y = c.base.rbegin()->y;
				stop_tmp.z = c.base.rbegin()->z;
			}
			if (stop_tmp.pl0 > stop_tmp.stoppl) {
				//start
				stop_tmp.rawid = c.base.begin()->rawid;
				stop_tmp.ax = c.base.begin()->ax;
				stop_tmp.ay = c.base.begin()->ay;
				stop_tmp.x = c.base.begin()->x;
				stop_tmp.y = c.base.begin()->y;
				stop_tmp.z = c.base.begin()->z;
			}
			tracks.insert(std::make_pair(stop_tmp.groupid, stop_tmp));

		}
		cnt++;
	}

	//printf("input fin %d track\n", cnt);

}
void clustering_2trk_vtx(std::multimap<int, stop_track>& tracks, int pl, std::map<int, vtx_point>& vtx) {
	double refz = 0; int utime;	for (auto itr1 = tracks.begin(); itr1 != tracks.end(); itr1++) {
		if (itr1->second.chainid == 0) {
			refz = itr1->second.z;
			utime = itr1->second.unixtime;
		}
	}

	//rawid,stop
	double zrange[2] = { 0,0 };
	if (pl <= 15 || (pl >= 16 && pl % 2 == 0)) {
		zrange[0] = -1000;
	}
	else if (pl % 2 == 1) {
		zrange[0] = -3200; //3500->3200
	}

	double extra[2];
	track_multi multi;
	multi.pl = pl;
	//全2trkのmd計算
	std::cout << tracks.size() << std::endl;
	for (auto itr1 = tracks.begin(); itr1 != tracks.end(); itr1++) {
		multi.eventid = itr1->second.groupid;
		multi.unixtime = utime;
		for (auto itr2 = std::next(itr1, 1); itr2 != tracks.end(); itr2++) {
			matrix_3D::vector_3D pos0, pos1, dir0, dir1;
			pos0.x = itr1->second.x;
			pos0.y = itr1->second.y;
			pos0.z = itr1->second.z;
			pos1.x = itr2->second.x;
			pos1.y = itr2->second.y;
			pos1.z = itr2->second.z;
			dir0.x = itr1->second.ax;
			dir0.y = itr1->second.ay;
			dir0.z = 1;
			dir1.x = itr2->second.ax;
			dir1.y = itr2->second.ay;
			dir1.z = 1;

			// pos0を基準にzrangeの範囲内で最近接距離をとる位置(extra)を探索
			double md = minimum_distance_fixed(pos0, pos1, dir0, dir1, zrange, extra, refz);
			track_pair pair_tmp;
			matrix_3D::vector_3D extra0 = addition(pos0, const_multiple(dir0, extra[0]));
			matrix_3D::vector_3D extra1 = addition(pos1, const_multiple(dir1, extra[1]));

			pair_tmp.x = (extra0.x + extra1.x) / 2;
			pair_tmp.y = (extra0.y + extra1.y) / 2;
			pair_tmp.z = (extra0.z + extra1.z) / 2;
			pair_tmp.dz = pair_tmp.z - refz;
			pair_tmp.eventid = multi.eventid;
			pair_tmp.md = md;
			pair_tmp.oa = matrix_3D::opening_angle(dir0, dir1);
			pair_tmp.t[0] = itr1->second;
			pair_tmp.t[1] = itr2->second;

			multi.pair.push_back(pair_tmp);
		}
	}


	matrix_3D::vector_3D p_vtx, pos, dir;
	vtx_point v;
	//加重平均でvtx pointの決定
	multi.x = 0;
	multi.y = 0;
	multi.z = 0;
	for (auto itr = multi.pair.begin(); itr != multi.pair.end(); itr++) {
		multi.x += itr->x;
		multi.y += itr->y;
		multi.z += itr->z;
	}
	v.x = multi.x / multi.pair.size();
	v.y = multi.y / multi.pair.size();
	v.z = multi.z / multi.pair.size();
	multi.dz = multi.z - refz;
	v.z = multi.dz;
	//各trkに対してIPの計算
	//for (auto itr = tracks.begin(); itr != tracks.end(); itr++) {
	//	matrix_3D::vector_3D pos0, pos1, dir0, dir1;
	//	pos0.x = itr->second.x;
	//	pos0.y = itr->second.y;
	//	pos0.z = itr->second.z;
	//	pos1.x = multi.x;
	//	pos1.y = multi.y;
	//	pos1.z = multi.z;
	//	dir0.x = itr->second.ax;
	//	dir0.y = itr->second.ay;
	//	dir0.z = 1;
	//	double ip = matrix_3D::inpact_parameter(pos0, dir0, pos1);
	//	multi.trk.push_back(std::make_pair(ip, itr->second));
	//}
	vtx.insert(std::make_pair(multi.eventid,v));


}
void add_vtx_information(std::vector<Momentum_recon::Event_information>& momch, std::map<int, vtx_point>& vtx) {
	stop_track stop_tmp;

	int pos = 0;
	int cnt = 0;
	for (auto& ev : momch) {
		auto itr = vtx.find(ev.groupid);
		ev.vertex_position[0] = itr->second.x;
		ev.vertex_position[1] = itr->second.y;
		ev.vertex_position[2] = itr->second.z;
		if (ev.chains.size() == 1) {
			ev.vertex_position[0] = -20000.0;
			ev.vertex_position[1] = -20000.0;
			ev.vertex_position[2] = -20000.0;

		}

	}
	//printf("input fin %d track\n", cnt);

}
