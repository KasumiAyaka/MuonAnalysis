//2024/11/13 kasumi
#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>

struct Position {
	double x, y;
};
bool operator<(const Position& lhs, const Position& rhs) {
	return std::tie(lhs.x, lhs.y) <	std::tie(rhs.x, rhs.y);
}

struct ViewArea {
	Position min, max;
};
bool operator<(const ViewArea& lhs, const ViewArea& rhs) {
	return std::tie(lhs.min, lhs.max) < std::tie(rhs.min, rhs.max);
}

struct AffineParam {
	double a, b, c, d, p, q;
};
struct Corrmap {
		//01      ‹æ‰æ - id
		//02 - 03   pos0, pos1   pos0(base pos), pos1
		//04 - 07   xmin, xmax, ymin, ymax
		//08 - 13   a, b, c, d, p, q = > afp(affine parameter for position)
		//14 - 19   a, b, c, d, p, q = > aft(affine parameter for angle)
		//20      dz
		//21 - 23   signal, background, S / N
		//24 - 25   rms_x, rms_y
		//26 - 27   not used
		//28 - 29   rms_ax, rms_ay
		//30      not used
		//31 - 32   ix, iy   two dimensional index of the view
		//33 - 38   flags - int
		//39 - 41   flags - double

	int AreaId;
	int pos1, pos2;
	ViewArea a;
	AffineParam afp, aft;
	double dz;
	double sig, bkg, sn;
	Position rms_pos, rms_ang,iview;
};

struct TrackerInfo {
	int unix_time, t_track_id, entry_in_daily_file;
	double	BM_mom;
};
bool operator<(const TrackerInfo& lhs, const TrackerInfo& rhs) {
	return std::tie(lhs.unix_time, lhs.t_track_id, lhs.entry_in_daily_file, lhs.BM_mom) <
		std::tie(rhs.unix_time, rhs.t_track_id, rhs.entry_in_daily_file, rhs.BM_mom);
}

void set_corrmap_old(std::string input, std::map<Position, Corrmap>& map);
void set_corrmap_new(std::string input, std::map<ViewArea, Corrmap>& map);
void culc_difference(std::map<ViewArea, Corrmap>& map_new, std::map<Position, Corrmap>& map_old, std::string output);

int main(int argc, char** argv) {
	if (argc == 3) {
		fprintf(stderr, "usage : original-corrmap.txt new-corrmap.txt output.txt \n");
		exit(1);
	}
	std::string in_old = argv[1];//input momch
	std::string in_new = argv[2];//multi
	std::string output = argv[3];//uniqe

	std::map<Position, Corrmap> morg;
	std::map<ViewArea, Corrmap> mnew;
	set_corrmap_old(in_old, morg);
	set_corrmap_new(in_new, mnew);
	culc_difference(mnew, morg,output);


}

void set_corrmap_old(std::string input, std::map<Position, Corrmap>& map)
{
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
	Corrmap c;
	Position p;
	while (std::getline(ifs, str)) {

		str_v = StringSplit_with_tab(str);
		//01      ‹æ‰æ - id
		//02 - 03   pos0, pos1   pos0(base pos), pos1
		//04 - 07   xmin, xmax, ymin, ymax
		//08 - 13   a, b, c, d, p, q = > afp(affine parameter for position)
		//14 - 19   a, b, c, d, p, q = > aft(affine parameter for angle)
		//20      dz
		//21 - 23   signal, background, S / N
		//24 - 25   rms_x, rms_y
		//26 - 27   not used
		//28 - 29   rms_ax, rms_ay
		//30      not used
		//31 - 32   ix, iy   two dimensional index of the view
		//33 - 38   flags - int
		//39 - 41   flags - double

		c.AreaId = std::stoi(str_v[0]);
		c.pos1 = std::stoi(str_v[1]);
		c.pos2 = std::stoi(str_v[2]);

		c.a.min.x = std::stod(str_v[3]);
		c.a.max.x = std::stod(str_v[4]);
		c.a.min.y = std::stod(str_v[5]);
		c.a.max.y = std::stod(str_v[6]);

		c.afp.a = std::stod(str_v[7]);
		c.afp.b = std::stod(str_v[8]);
		c.afp.c = std::stod(str_v[9]);
		c.afp.d = std::stod(str_v[10]);
		c.afp.p = std::stod(str_v[11]);
		c.afp.q = std::stod(str_v[12]);

		c.aft.a = std::stod(str_v[13]);
		c.aft.b = std::stod(str_v[14]);
		c.aft.c = std::stod(str_v[15]);
		c.aft.d = std::stod(str_v[16]);
		c.aft.p = std::stod(str_v[17]);
		c.aft.q = std::stod(str_v[18]);

		c.dz = std::stod(str_v[19]);
		c.sig = std::stod(str_v[20]);
		c.bkg = std::stod(str_v[21]);
		c.sn = std::stod(str_v[22]);
		c.rms_pos.x = std::stod(str_v[23]);
		c.rms_pos.y = std::stod(str_v[24]);
		c.rms_ang.x = std::stod(str_v[27]);
		c.rms_ang.y = std::stod(str_v[28]);
		c.iview.x = std::stod(str_v[30]);
		c.iview.y = std::stod(str_v[31]);

		p.x = (c.a.max.x + c.a.min.x) * 0.5;
		p.y = (c.a.max.y + c.a.min.y) * 0.5;
		
		map.insert(std::make_pair(p, c));
		count++;
	}
	std::cout << "# of Original corrmap's View : " << count << std::endl;

}

void set_corrmap_new(std::string input,  std::map<ViewArea,Corrmap> &map)
{
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
	Corrmap c;
	while (std::getline(ifs, str)) {

		str_v = StringSplit_with_tab(str);
		//01      ‹æ‰æ - id
		//02 - 03   pos0, pos1   pos0(base pos), pos1
		//04 - 07   xmin, xmax, ymin, ymax
		//08 - 13   a, b, c, d, p, q = > afp(affine parameter for position)
		//14 - 19   a, b, c, d, p, q = > aft(affine parameter for angle)
		//20      dz
		//21 - 23   signal, background, S / N
		//24 - 25   rms_x, rms_y
		//26 - 27   not used
		//28 - 29   rms_ax, rms_ay
		//30      not used
		//31 - 32   ix, iy   two dimensional index of the view
		//33 - 38   flags - int
		//39 - 41   flags - double

		c.AreaId=std::stoi(str_v[0]);
		c.pos1 = std::stoi(str_v[1]);
		c.pos2 = std::stoi(str_v[2]);

		c.a.min.x = std::stod(str_v[3]);
		c.a.max.x = std::stod(str_v[4]);
		c.a.min.y = std::stod(str_v[5]);
		c.a.max.y = std::stod(str_v[6]);

		c.afp.a = std::stod(str_v[7]);
		c.afp.b = std::stod(str_v[8]);
		c.afp.c = std::stod(str_v[9]);
		c.afp.d = std::stod(str_v[10]);
		c.afp.p = std::stod(str_v[11]);
		c.afp.q = std::stod(str_v[12]);

		c.aft.a = std::stod(str_v[13]);
		c.aft.b = std::stod(str_v[14]);
		c.aft.c = std::stod(str_v[15]);
		c.aft.d = std::stod(str_v[16]);
		c.aft.p = std::stod(str_v[17]);
		c.aft.q = std::stod(str_v[18]);

		c.dz = std::stod(str_v[19]);
		c.sig = std::stod(str_v[20]);
		c.bkg = std::stod(str_v[21]);
		c.sn = std::stod(str_v[22]);
		c.rms_pos.x = std::stod(str_v[23]);
		c.rms_pos.y = std::stod(str_v[24]);
		c.rms_ang.x = std::stod(str_v[27]);
		c.rms_ang.y = std::stod(str_v[28]);
		c.iview.x = std::stod(str_v[30]);
		c.iview.y = std::stod(str_v[31]);

		map.insert(std::make_pair(c.a,c));
		count++;
	}
	std::cout << "# of Retaken corrmap's View : " << count << std::endl;

}

void culc_difference(std::map<ViewArea, Corrmap>& map_new, std::map<Position, Corrmap>& map_old, std::string output)
{
	std::ofstream ofs(output);

	int count = 0;
	for (auto itr = map_new.begin(); itr != map_new.end(); itr++) {
		auto p = map_old.lower_bound(itr->first.min);
		auto q = map_old.upper_bound(itr->first.max);

		for (auto it = p; it != q; it++) {
			if (it->first.x >= itr->first.min.x && it->first.x <= itr->first.max.x && it->first.y >= itr->first.min.y && it->first.y <= itr->first.max.y) {
				ofs << std::fixed << std::right
					<< std::setw(10) << std::setprecision(1) << it->first.x << " "
					<< std::setw(10) << std::setprecision(1) << it->first.y << " "
					<< std::setw(10) << std::setprecision(1) << it->second.a.min.x << " "
					<< std::setw(10) << std::setprecision(1) << it->second.a.max.x << " "
					<< std::setw(10) << std::setprecision(1) << it->second.a.min.y << " "
					<< std::setw(10) << std::setprecision(1) << it->second.a.max.y << " "
					<< std::setw(10) << std::setprecision(1) << (itr->first.min.x + itr->first.max.x) * 0.5 << " "
					<< std::setw(10) << std::setprecision(1) << (itr->first.min.y + itr->first.max.y) * 0.5 << " "
					<< std::setw(10) << std::setprecision(1) << itr->first.min.x << " "
					<< std::setw(10) << std::setprecision(1) << itr->first.max.x << " "
					<< std::setw(10) << std::setprecision(1) << itr->first.min.y << " "
					<< std::setw(10) << std::setprecision(1) << itr->first.max.y << " "
					<< std::setw(10) << std::setprecision(6) << it->second.afp.a - itr->second.afp.a << " "
					<< std::setw(10) << std::setprecision(6) << it->second.afp.b - itr->second.afp.b << " "
					<< std::setw(10) << std::setprecision(6) << it->second.afp.c - itr->second.afp.c << " "
					<< std::setw(10) << std::setprecision(6) << it->second.afp.d - itr->second.afp.d << " "
					<< std::setw(10) << std::setprecision(1) << it->second.afp.p - itr->second.afp.p << " "
					<< std::setw(10) << std::setprecision(1) << it->second.afp.q - itr->second.afp.q << " "
					<< std::setw(10) << std::setprecision(6) << it->second.aft.a - itr->second.aft.a << " "
					<< std::setw(10) << std::setprecision(6) << it->second.aft.b - itr->second.aft.b << " "
					<< std::setw(10) << std::setprecision(6) << it->second.aft.c - itr->second.aft.c << " "
					<< std::setw(10) << std::setprecision(6) << it->second.aft.d - itr->second.aft.d << " "
					<< std::setw(10) << std::setprecision(6) << it->second.aft.p - itr->second.aft.p << " "
					<< std::setw(10) << std::setprecision(6) << it->second.aft.q - itr->second.aft.q
					<< std::endl;
				count++;
			}

		}
		if (count == 0) {
			ofs << std::fixed << std::right
				<< std::setw(10) << std::setprecision(1) << -3000000 << " "
				<< std::setw(10) << std::setprecision(1) << -3000000 << " "
				<< std::setw(10) << std::setprecision(1) << -3000000 << " "
				<< std::setw(10) << std::setprecision(1) << -3000000 << " "
				<< std::setw(10) << std::setprecision(1) << -3000000 << " "
				<< std::setw(10) << std::setprecision(1) << -3000000 << " "
				<< std::setw(10) << std::setprecision(1) << (itr->first.min.x + itr->first.max.x) * 0.5 << " "
				<< std::setw(10) << std::setprecision(1) << (itr->first.min.y + itr->first.max.y) * 0.5 << " "
				<< std::setw(10) << std::setprecision(1) << itr->first.min.x << " "
				<< std::setw(10) << std::setprecision(1) << itr->first.max.x << " "
				<< std::setw(10) << std::setprecision(1) << itr->first.min.y << " "
				<< std::setw(10) << std::setprecision(1) << itr->first.max.y << " "
				<< std::setw(10) << std::setprecision(6) << -3000000 << " "
				<< std::setw(10) << std::setprecision(6) << -3000000 << " "
				<< std::setw(10) << std::setprecision(6) << -3000000 << " "
				<< std::setw(10) << std::setprecision(6) << -3000000 << " "
				<< std::setw(10) << std::setprecision(1) << -3000000 << " "
				<< std::setw(10) << std::setprecision(1) << -3000000 << " "
				<< std::setw(10) << std::setprecision(6) << -3000000 << " "
				<< std::setw(10) << std::setprecision(6) << -3000000 << " "
				<< std::setw(10) << std::setprecision(6) << -3000000 << " "
				<< std::setw(10) << std::setprecision(6) << -3000000 << " "
				<< std::setw(10) << std::setprecision(6) << -3000000 << " "
				<< std::setw(10) << std::setprecision(6) << -3000000
				<< std::endl;
			//std::cout << "non" << std::endl;
			//continue;
		}
		count = 0;
	}
}
