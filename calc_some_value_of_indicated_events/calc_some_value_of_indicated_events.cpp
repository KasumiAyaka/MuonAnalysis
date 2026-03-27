//calc_some_value_of_indicated_events
// 2025/01/22
// 2024/09/07
// kasumi
// based on "Check_upstream_base.cpp

#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>
#include <iomanip>

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

struct Basetrack {
	Key ky;
	int vph, rid;
	double x, y, z, ax, ay;
};
struct Btrk {
	Basetrack b0, b1;
};

//void calc_md(mfile0::Mfile& m_all, std::map<int, int>& key);
void SetTracks(mfile0::Mfile& m_all, std::multimap<int, Btrk>& trk);
//void Calc_MinimumDistance_zrange(std::multimap<Key, Btrk>& trk, Key k[2]);
//void Calc_AngDiff(std::multimap<Key, Btrk>& trk, Key k[2]);
//void Calc_PosDiff(std::multimap<Key, Btrk>& trk, Key k[2]);

void SetEventList(std::string input, std::set<int>& list, int eccnum);
void MatchTracks(std::set<int>& list, std::multimap<int, Btrk>& trk, std::string output, int eccnum);
void Calc_MinimumDistance(Btrk bk[2],  std::ofstream& ofs, int eccnum);

int main(int argc, char** argv) {
	if (argc != 5) {
		fprintf(stderr, "usage:file-in.all list.txt #ecc output.txt\n");
		exit(1);
	}
	std::string file_in_mfile = argv[1];
	std::string file_in_list = argv[2];
	int eccnum = std::stoi(argv[3]);
	std::string file_out = argv[4];

	mfile0::Mfile all;
	mfile1::read_mfile_extension(file_in_mfile, all);

	std::multimap<int, Btrk> trk;
	SetTracks(all, trk);

	int a = 0;
	int mode = 0;
	Key k[2];
	std::set<int> list;
	SetEventList(file_in_list, list,eccnum);
	MatchTracks(list, trk, file_out, eccnum);
}

void SetEventList(std::string input, std::set<int> &list,int eccnum) {

	std::ifstream ifs(input);
		if (!ifs) {
			std::cerr << "file open error" << std::endl;
			exit(1);
		}

		std::string str;					//1strein into
		std::vector<std::string> str_v;		//input 1 ward
		std::string buffer;

		std::cout << input << std::endl;
		int count = 0;
		while (std::getline(ifs, str)) {

			str_v = StringSplit_with_tab(str);
			if (std::stoi(str_v[0]) == eccnum) {
				list.insert(std::stoi(str_v[1]));
				//map.insert(std::make_pair(b.eid, b));
				count++;
				//debag
				//std::cout << std::stoi(str_v[1]) << std::endl;
			}
		}

}

void SetTracks(mfile0::Mfile& m_all, std::multimap<int, Btrk>& trk) {

	Btrk tmp;

	int i = 0;
	for (auto itr = m_all.chains.begin(); itr != m_all.chains.end(); itr++) {
		tmp.b0.ky.cid = itr->chain_id;
		tmp.b0.ky.gid = itr->basetracks[0].group_id;

		tmp.b0.ky.pl = itr->basetracks.begin()->pos / 10;
		tmp.b0.ax = itr->basetracks.begin()->ax;
		tmp.b0.ay = itr->basetracks.begin()->ay;
		tmp.b0.rid = itr->basetracks.begin()->rawid;
		tmp.b0.vph = itr->basetracks.begin()->ph % 10000;
		tmp.b0.x = itr->basetracks.begin()->x;
		tmp.b0.y = itr->basetracks.begin()->y;
		tmp.b0.z = itr->basetracks.begin()->z;

		auto itr_b = std::prev(itr->basetracks.end(), 1);
		tmp.b1.ky.pl = itr_b->pos / 10;
		tmp.b1.ky.cid = itr->chain_id;
		tmp.b1.ky.gid = itr->basetracks[0].group_id;
		tmp.b1.ax = itr_b->ax;
		tmp.b1.ay = itr_b->ay;
		tmp.b1.rid = itr_b->rawid;
		tmp.b1.vph = itr_b->ph % 10000;
		tmp.b1.x = itr_b->x;
		tmp.b1.y = itr_b->y;
		tmp.b1.z = itr_b->z;
		trk.insert(std::make_pair(tmp.b0.ky.gid, tmp));
		//std::cout << tmp.b0.ky.gid << std::endl;

	}

}

void MatchTracks(std::set<int>& list, std::multimap<int, Btrk>& trk, std::string output, int eccnum) {
	std::multimap<int, int> map;

	Btrk bk[2];
	int i;
	for (auto itr = list.begin(); itr != list.end(); itr++) {
		std::cout << *itr << std::endl;
		auto p = trk.equal_range(*itr);

		if (p.second == trk.end()) {
			std::cout << "not exist" << std::endl;
		}

		i = 0;
		for (auto itr0 = p.first; itr0 != p.second; itr0++) {
			std::cout<<itr0->second.b0.ky.gid <<std::endl;
			bk[i]= itr0->second;
			i++;
		}
		std::ofstream ofs(output,std::ios::app);

		Calc_MinimumDistance(bk,ofs,eccnum);

	}


}


void Calc_MinimumDistance(Btrk bk[2], std::ofstream &ofs, int eccnum) {

	matrix_3D::vector_3D pos0, pos1, dir0, dir1;

	if (bk[0].b0.ky.cid == 0) {//bk[0]-->muon
		//muon
		pos0.x = bk[0].b1.x;
		pos0.y = bk[0].b1.y;
		pos0.z = bk[0].b1.z;
		dir0.x = bk[0].b1.ax;
		dir0.y = bk[0].b1.ay;
		dir0.z = 1;

		//pion
		pos1.x = bk[1].b0.x;
		pos1.y = bk[1].b0.y;
		pos1.z = bk[1].b0.z;
		dir1.x = bk[1].b0.ax;
		dir1.y = bk[1].b0.ay;
		dir1.z = 1;
		std::cout << " ( Groupid, Chainid, PL, rawid )" << std::endl;
		std::cout << std::fixed << std::right << std::setfill(' ')
			<< "trk1 :  ("
			<< std::setw(12) << std::setprecision(0) << bk[0].b0.ky.gid << ", "
			<< std::setw(12) << std::setprecision(0) << bk[0].b1.ky.cid << ", "
			<< std::setw(3) << std::setprecision(0) << bk[0].b1.ky.pl << ", "
			<< std::setw(12) << std::setprecision(0) << bk[0].b1.rid << ")\n"
			<< "trk2 :  ("
			<< std::setw(12) << std::setprecision(0) << bk[1].b0.ky.gid << ", "
			<< std::setw(12) << std::setprecision(0) << bk[1].b0.ky.cid << ", "
			<< std::setw(3) << std::setprecision(0) << bk[1].b0.ky.pl << ", "
			<< std::setw(12) << std::setprecision(0) << bk[1].b0.rid << ")\n"
			<< std::endl;
	}
	else {//bk[1]-->muon
		//muon
		pos0.x = bk[1].b1.x;
		pos0.y = bk[1].b1.y;
		pos0.z = bk[1].b1.z;
		dir0.x = bk[1].b1.ax;
		dir0.y = bk[1].b1.ay;
		dir0.z = 1;

		//pion
		pos1.x = bk[0].b0.x;
		pos1.y = bk[0].b0.y;
		pos1.z = bk[0].b0.z;
		dir1.x = bk[0].b0.ax;
		dir1.y = bk[0].b0.ay;
		dir1.z = 1;
		std::cout << " ( Groupid, Chainid, PL, rawid )" << std::endl;
		std::cout << std::fixed << std::right
			<< "trk1 :  ("
			<< std::setw(12) << std::setprecision(0) << bk[1].b0.ky.gid << ", "
			<< std::setw(12) << std::setprecision(0) << bk[1].b0.ky.cid << ", "
			<< std::setw(3) << std::setprecision(0) << bk[1].b0.ky.pl << ", "
			<< std::setw(12) << std::setprecision(0) << bk[1].b0.rid << ")\n"
			<< "trk2 :  ("
			<< std::setw(12) << std::setprecision(0) << bk[0].b0.ky.gid << ", "
			<< std::setw(12) << std::setprecision(0) << bk[0].b1.ky.cid << ", "
			<< std::setw(3) << std::setprecision(0) << bk[0].b1.ky.pl << ", "
			<< std::setw(12) << std::setprecision(0) << bk[0].b1.rid << ")\n"
			<< std::endl;
	}




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

	ofs << std::fixed << std::right
		<< std::setw(2) << std::setprecision(0) << eccnum << " "
		<< std::setw(5) << std::setprecision(0) << bk[0].b0.ky.gid << " "
		<< std::setw(10) << std::setprecision(1) << md << " "
		<< std::setw(10) << std::setprecision(5) << oa << " "
		<< std::setw(10) << std::setprecision(1) << point[0] << " "
		<< std::setw(10) << std::setprecision(1) << point[1] << " "
		<< std::setw(10) << std::setprecision(1) << point[2]
		<< std::endl;

}

/*
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


	double oa, md, dz;
	double extra[2], z_range[2];
	z_range[0] = pos0.z;
	z_range[1] = pos1.z;
	oa = matrix_3D::opening_angle(dir0, dir1);
	md = matrix_3D::minimum_distance(pos0, pos1, dir0, dir1, z_range, extra);
	dz = (fabs(extra[0]) + fabs(extra[1])) / 2;

	std::cout << std::fixed << std::right
		<< "   md    = " << std::setw(10) << std::setprecision(1) << md << "\n"
		<< "   oa    =     " << std::setw(10) << std::setprecision(5) << oa << "\n"
		<< "   dz    = " << std::setw(10) << std::setprecision(1) << dz
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

	double ax, ay, dax, day, dlat, drad, angle;
	// dax,day
	dax = itr2->second.ax - itr1->second.ax;
	day = itr2->second.ay - itr1->second.ay;
	std::cout << std::fixed << std::right
		<< "   dax   = " << std::setw(10) << std::setprecision(5) << dax << "\n"
		<< "   day   = " << std::setw(10) << std::setprecision(5) << day
		<< std::endl;


	ax = itr2->second.ax + itr1->second.ax;
	ay = itr2->second.ay + itr1->second.ay;
	angle = sqrt(ax * ax + ay * ay);
	dlat = 0.5 * (-ay * dax + ax * day) / angle;
	drad = 0.5 * (ax * dax + ay * day) / angle;
	if (angle < 0.01) {
		dlat = dax;
		drad = day;
	}

	// d(tan)
	double tan1, tan2, dtan;
	tan1 = sqrt(itr1->second.ax * itr1->second.ax + itr1->second.ay * itr1->second.ay);
	tan2 = sqrt(itr2->second.ax * itr2->second.ax + itr2->second.ay * itr2->second.ay);
	dtan = tan2 - tan1;

	std::cout << std::fixed << std::right
		<< "   dal   = " << std::setw(10) << std::setprecision(5) << dlat << "\n"
		<< "   dar   = " << std::setw(10) << std::setprecision(5) << drad << "\n"
		<< "   dtan  = " << std::setw(10) << std::setprecision(5) << dtan << "\n"
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
	matrix_3D::vector_3D pos0, pos1, point;

	pos0.x = -itr1->second.x;
	pos0.y = -itr1->second.y;
	pos0.z = -itr1->second.z;

	pos1.x = itr2->second.x;
	pos1.y = itr2->second.y;
	pos1.z = itr2->second.z;


	double interporate;
	interporate = std::fabs(itr2->second.z - itr1->second.z);
	matrix_3D::vector_3D p0, p1;
	p0.x =


		std::cout << std::fixed << std::right
		<< "pos-difference"
		<< "  (dx,dy,dz) =  " << std::setw(10) << std::setprecision(1)
		<< std::setw(10) << std::setprecision(1) << pos1.x + pos0.x << ","
		<< std::setw(10) << std::setprecision(1) << pos1.y + pos0.y << ","
		<< std::setw(10) << std::setprecision(1) << pos1.z + pos0.z << ")\n"
		<< std::endl;


	point = matrix_3D::addition(pos0, pos1);
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
*/