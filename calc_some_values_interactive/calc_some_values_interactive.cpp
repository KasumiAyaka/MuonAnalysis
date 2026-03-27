// 2024/09/07
// kasumi
// based on "Check_upstream_base.cpp

#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>
#include <iomanip>
#include <cmath>

class output_format {
public:
	int groupid, chainid, pl, upl;
	int peke, count;
	double dal, dar, md, dz, dz2, oa;
};


struct Key {
	int pl, cid, gid;
};
bool operator<(const Key& lhs, const Key& rhs) {
	return std::tie(lhs.gid, lhs.cid, lhs.pl) < std::tie(rhs.gid, rhs.cid, rhs.pl);
}

struct Btrk {
	Key ky;
	int vph, rid;
	double x, y, z, ax, ay;
};

void calc_md(mfile0::Mfile& m_all, std::map<int, int>& key);
void SetTracks(mfile0::Mfile& m_all, std::multimap<Key, Btrk>& trk);
void Calc_MinimumDistance(std::multimap<Key, Btrk>& trk, Key k[2]);
void Calc_MinimumDistance_zrange(std::multimap<Key, Btrk>& trk, Key k[2]);
void Calc_AngDiff(std::multimap<Key, Btrk>& trk, Key k[2]);
void Calc_PosDiff(std::multimap<Key, Btrk>& trk, Key k[2]);


int main(int argc, char** argv) {
	if (argc!=2) {
		fprintf(stderr, "usage:file-in.all\n");
		exit(1);
	}
	std::string file_in_mfile = argv[1];

	mfile0::Mfile all;
	mfile1::read_mfile_extension(file_in_mfile, all);

	std::multimap<Key, Btrk> trk;
	SetTracks(all, trk);

	int a = 0;
	int mode = 0;
	Key k[2];

	std::cout << std::setfill(' ');
	std::cout << "\n Input number : \n"
		<< "   finish    --> -1\n"
		<< "   continue  -->  1"
		<< std::endl;
	std::cin >> a;
	std::cout << std::endl;

	do {

		if (a > 0) {
			std::cout << " Select mode : \n"
				<< "  finish                         --> -1\n"
				<< " ---------------------------------------\n"
				<< "  Calculation\n"
				<< "    * Minimum Distane            --> 1\n"
				<< "    * Angle Distane              --> 2\n"
				<< "    * Position Distane           --> 3\n"
				<< "    * MD with z-range striction  --> 4\n"
				<< std::endl;

			std::cin >> mode;
			a = mode;
			std::cout << std::endl;

			if (mode > 0 && mode < 5) {
				std::cout << " Input basetrack1 :  [group-id] [chain-id] [pl]" << std::endl;
				std::cin >> k[0].gid >> k[0].cid >> k[0].pl;
				std::cout << "\n Input basetrack2 : [group-id] [chain-id] [pl]" << std::endl;
				std::cin >> k[1].gid >> k[1].cid >> k[1].pl;
				std::cout << std::endl;

				if (trk.find(k[0]) == trk.end() || trk.find(k[1]) == trk.end()) {
					std::cout << "\n   not exist that basetrack!\n" << std::endl;
					mode = 100;
				}

			}

			if (mode == 1) {
				std::cout << " * Calc md " << std::endl;
				Calc_MinimumDistance(trk, k);
			}
			else if (mode == 2) {
				std::cout << " * Calc andle diff " << std::endl;
				Calc_AngDiff(trk, k);
			}
			else if (mode == 3) {
				std::cout << " * Calc position diff " << std::endl;
				Calc_PosDiff(trk, k);
			}
			else if (mode == 4) {
				std::cout << " * Calc md (with Z-Range) " << std::endl;
				Calc_MinimumDistance_zrange(trk, k);
			}
			else if (mode == -1) {
				a == -1;
			}
			else {
				mode = 100;
			}

		}

	} while (a > 0);

	std::cout << " Finish process...!" << std::endl;
}

void SetTracks(mfile0::Mfile& m_all, std::multimap<Key, Btrk>& trk) {

	Key ky;
	Btrk tmp;

	int i = 0;
	for (auto itr = m_all.chains.begin(); itr != m_all.chains.end(); itr++) {
		ky.cid = itr->chain_id;
		ky.gid = itr->basetracks[0].group_id;

		for (auto itr_b = itr->basetracks.begin(); itr_b != itr->basetracks.end(); itr_b++) {
			ky.pl = itr_b->pos / 10;
			tmp.ax = itr_b->ax;
			tmp.ay = itr_b->ay;
			tmp.rid = itr_b->rawid;
			tmp.vph = itr_b->ph % 10000;
			tmp.x = itr_b->x;
			tmp.y = itr_b->y;
			tmp.z = itr_b->z;

			trk.insert(std::make_pair(ky, tmp));

		}
	}

}

void Calc_MinimumDistance(std::multimap<Key, Btrk>& trk, Key k[2]) {

	if (trk.count(k[0]) != 1) {
		std::cerr << "duplicated!" << std::endl;
		return;
	}

	auto itr1 = trk.find(k[0]);
	auto itr2 = trk.find(k[1]);

	std::cout << " ( PL, rawid )" << std::endl;
	std::cout << std::fixed << std::right
		<< "trk1 :  (" << std::setw(3) << std::setprecision(0) << itr1->first.pl << ", "
		<< std::setw(12) << std::setprecision(0) << itr1->second.rid << ")\n"
		<< "trk2 :  (" << std::setw(3) << std::setprecision(0) << itr2->first.pl << ", "
		<< std::setw(12) << std::setprecision(0) << itr2->second.rid << ")\n"
		<< std::endl;

	matrix_3D::vector_3D pos0, pos1, dir0, dir1;

	pos0.x = itr1->second.x;
	pos0.y = itr1->second.y;
	pos0.z = itr1->second.z;
	dir0.x = itr1->second.ax;
	dir0.y = itr1->second.ay;
	dir0.z = 1;

	pos1.x = itr2->second.x;
	pos1.y = itr2->second.y;
	pos1.z = itr2->second.z;
	dir1.x = itr2->second.ax;
	dir1.y = itr2->second.ay;
	dir1.z = 1;


	double oa, md;
	double point[3];
	md = matrix_3D::minimum_distance(pos0, pos1, dir0, dir1, point);
	oa = matrix_3D::opening_angle(dir0, dir1);

	std::cout << std::fixed << std::right
		<< "   md    =  " << std::setw(10) << std::setprecision(1) << md << "\n"
		<< "   oa    =      " << std::setw(10) << std::setprecision(5) << oa << "\n"
		<< "   point = ("
		<< std::setw(10) << std::setprecision(1) << point[0] << ", "
		<< std::setw(10) << std::setprecision(1) << point[1] << ", "
		<< std::setw(10) << std::setprecision(1) << point[2] << ")\n"
		<< std::endl;

}

void Calc_MinimumDistance_zrange(std::multimap<Key, Btrk>& trk, Key k[2]) {

	if (trk.count(k[0]) != 1) {
		std::cerr << "duplicated!" << std::endl;
		return;
	}

	auto itr1 = trk.find(k[0]);
	auto itr2 = trk.find(k[1]);
	std::cout << "PL :  " << itr1->first.pl << ", " << itr2->first.pl << std::endl;

	matrix_3D::vector_3D pos0, pos1, dir0, dir1;

	pos0.x = itr1->second.x;
	pos0.y = itr1->second.y;
	pos0.z = itr1->second.z;
	dir0.x = itr1->second.ax;
	dir0.y = itr1->second.ay;
	dir0.z = 1;

	pos1.x = itr2->second.x;
	pos1.y = itr2->second.y;
	pos1.z = itr2->second.z;
	dir1.x = itr2->second.ax;
	dir1.y = itr2->second.ay;
	dir1.z = 1;


	double oa, md,dz;
	double extra[2], z_range[2];
	z_range[0] = pos0.z;
	z_range[1] = pos1.z;
	oa = matrix_3D::opening_angle(dir0, dir1);
	md = matrix_3D::minimum_distance(pos0, pos1, dir0, dir1, z_range, extra);
	dz = (fabs(extra[0]) + fabs(extra[1])) / 2;

	std::cout << std::fixed << std::right
		<< "   md    = " << std::setw(10) << std::setprecision(1) << md << "\n"
		<< "   oa    =     " << std::setw(10) << std::setprecision(5) << oa << "\n"
		<< "   dz    = "<< std::setw(10) << std::setprecision(1) << dz 
		<< std::endl;

}

void Calc_AngDiff(std::multimap<Key, Btrk>& trk, Key k[2]) {

	if (trk.count(k[0]) != 1) {
		std::cerr << "duplicated!" << std::endl;
		return;
	}

	auto itr1 = trk.find(k[0]);
	auto itr2 = trk.find(k[1]);
	std::cout << "PL :  " << itr1->first.pl << ", " << itr2->first.pl << std::endl;

	double ax, ay, dax, day, dlat, drad,angle;
	// dax,day
	dax = itr2->second.ax - itr1->second.ax;
	day = itr2->second.ay - itr1->second.ay;
	std::cout << std::fixed << std::right
		<< "   dax   = " << std::setw(10) << std::setprecision(5) << dax << "\n"
		<< "   day   = " << std::setw(10) << std::setprecision(5) << day
		<< std::endl;


	// d(lat)
	/* 
	ax = itr1->second.ax;
	ay = itr1->second.ay;
	angle = sqrt(ax * ax + ay * ay);

	dax = itr2->second.ax - ax;
	day = itr2->second.ay - ay;
	if (angle < 0.01) {
		dlat = dax;
	}
	else {
		dlat = (dax * ay - day * ax) / angle;
	}

	// d(rad)
	if (angle < 0.01) {
		drad = day;
	}
	else {
		drad = (dax * ax + day * ay) / angle;
	}
	*/
	ax = itr2->second.ax + itr1->second.ax;
	ay = itr2->second.ay + itr1->second.ay;
	angle = sqrt(ax * ax + ay * ay);
	dlat = (-ay * dax + ax * day) / angle;
	drad = (ax * dax + ay * day) / angle;
	if (angle < 0.01) {
		dlat = dax;
		drad = day;
	}
	std::cout << std::fixed << std::right
		<< "   dal   = " << std::setw(10) << std::setprecision(5) << dlat << "\n"
		<< "   dar   = " << std::setw(10) << std::setprecision(5) << drad 
		<< std::endl;

	// d(tan)
	double tan1, tan2, dtan;
	tan1 = sqrt(itr1->second.ax * itr1->second.ax + itr1->second.ay * itr1->second.ay);
	tan2 = sqrt(itr2->second.ax * itr2->second.ax + itr2->second.ay * itr2->second.ay);

	tan1 = atan(tan1);// theta [rad]
	tan2 = atan(tan2);// theta [rad]
	dtan = tan2 - tan1;// d(theta) [rad]



	std::cout << std::fixed << std::right
		<< "   dtan  = " << std::setw(10) << std::setprecision(5) << dtan << " [rad], ( tan = " << tan(dtan) << " )\n"
		<< std::endl;

}

void Calc_PosDiff(std::multimap<Key, Btrk>& trk, Key k[2]) {

	if (trk.count(k[0]) != 1) {
		std::cerr << "duplicated!" << std::endl;
		return;
	}

	auto itr1 = trk.find(k[0]);
	auto itr2 = trk.find(k[1]);
	std::cout << "PL :  " << itr1->first.pl << ", " << itr2->first.pl << std::endl;

	double extra[2], z_range[2];
	matrix_3D::vector_3D pos0, pos1,point;

	pos0.x = -itr1->second.x;
	pos0.y = -itr1->second.y;
	pos0.z = -itr1->second.z;

	pos1.x = itr2->second.x;
	pos1.y = itr2->second.y;
	pos1.z = itr2->second.z;

	point=matrix_3D::addition(pos0, pos1);
	double d = matrix_3D::distance(pos0, pos1);

	std::cout << std::fixed << std::right
		<< "   distance   =  " << std::setw(10) << std::setprecision(1) << d << "\n"
		<< "   (dx,dy,dz) = ("
		<< std::setw(10) << std::setprecision(1) << point.x << ","
		<< std::setw(10) << std::setprecision(1) << point.y << ","
		<< std::setw(10) << std::setprecision(1) << point.z << ")\n"
		<< std::endl;

}

void calc_md(mfile0::Mfile& m_all, std::map<int, int>& key) {

	std::vector<mfile0::M_Chain> ret;
	matrix_3D::vector_3D pos0, pos1, dir0, dir1;
	int pl0, pl1;
	double z0, z1;

	Key ky;
	Btrk tmp, bt[2];
	std::multimap<Key, Btrk> btk;

	int i = 0;
	for (auto itr = m_all.chains.begin(); itr != m_all.chains.end(); itr++) {

		auto k = key.find(itr->chain_id);

		for (auto itr_b = itr->basetracks.begin(); itr_b != itr->basetracks.end(); itr_b++) {

			if (itr_b->pos / 10 == k->second) {
				tmp.ax = itr_b->ax;
				tmp.ay = itr_b->ay;
				tmp.rid = itr_b->rawid;
				tmp.vph = itr_b->ph % 10000;
				tmp.x = itr_b->x;
				tmp.y = itr_b->y;
				tmp.z = itr_b->z;

				bt[i] = tmp;
				i++;
			}
		}
		if (i == 2)break;
	}


	pos0.x = bt[0].x;
	pos0.y = bt[0].y;
	pos0.z = bt[0].z;
	dir0.x = bt[0].ax;
	dir0.y = bt[0].ay;
	dir0.z = 1;

	pos1.x = bt[1].x;
	pos1.y = bt[1].y;
	pos1.z = bt[1].z;
	dir1.x = bt[1].ax;
	dir1.y = bt[1].ay;
	dir1.z = 1;


	double oa, md;
	double extra[2], z_range[2];
	z_range[0] = pos0.z;
	z_range[1] = pos1.z;
	oa = matrix_3D::opening_angle(dir0, dir1);
	//md = matrix_3D::minimum_distance(pos0, pos1, dir0, dir1, z_range, extra);
	double point[3];
	md = matrix_3D::minimum_distance(pos0, pos1, dir0, dir1, point);

	std::cout << std::fixed << std::right
		<< "   md    = " << std::setw(10) << std::setprecision(1) << md << "\n"
		<< "   oa    = " << std::setw(10) << std::setprecision(5) << oa << "\n"
		<< "   point = "
		<< std::setw(10) << std::setprecision(1) << point[0] << ","
		<< std::setw(10) << std::setprecision(1) << point[1] << ","
		<< std::setw(10) << std::setprecision(1) << point[2] << "\n"
		<< std::endl;


}
