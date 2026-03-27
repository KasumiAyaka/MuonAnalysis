// suzuki prg
// d(lateral)ÇÃåvéZåÎÇËÇÃÇΩÇﬂèCê≥
#define _CRT_SECURE_NO_WARNINGS

#pragma comment(lib,"FILE_structure.lib")
#pragma comment(lib,"functions.lib")

#include <FILE_structure.hpp>
#include <functions.hpp>

#include <chrono>
#include <set>

class Linklet_rl {
public:
	netscan::linklet_t link;
	double dr, dl, md, oa;
};

std::vector<Linklet_rl> Radial_Laterai_information_addition(std::vector<netscan::linklet_t> link);
void position_difference(netscan::linklet_t link, double& dr, double& dl);
void Calc_MD_OA(netscan::linklet_t link, double& md, double& oa);

void write_linklet_txt(std::string filename, std::vector<Linklet_rl> l);
void write_linklet_bin(std::string filename, std::vector<Linklet_rl> l);

int main(int argc, char** argv) {
	if (argc != 3) {
		fprintf(stderr, "usage:prg input-linklet(bin) output-linklet(bin)\n");
		exit(1);
	}
	std::string file_in_link = argv[1];
	std::string file_out_link = argv[2];

	std::vector<netscan::linklet_t> link;
	netscan::read_linklet_bin(file_in_link, link);
	std::vector<Linklet_rl> link_rl = Radial_Laterai_information_addition(link);
	write_linklet_bin(file_out_link, link_rl);
	//write_linklet_txt(file_out_link, link_rl);

}
std::vector<Linklet_rl> Radial_Laterai_information_addition(std::vector<netscan::linklet_t> link) {
	std::vector<Linklet_rl> ret;
	int count = 0;
	for (auto itr = link.begin(); itr != link.end(); itr++) {
		if (count % 10000 == 0) {
			fprintf(stderr, "\r calculation ... %d/%d(%4.1lf%%)", count, link.size(), count * 100. / link.size());
		}
		count++;
		Linklet_rl tmp_link;
		tmp_link.link = *itr;
		position_difference(*itr, tmp_link.dr, tmp_link.dl);
		Calc_MD_OA(*itr, tmp_link.md, tmp_link.oa);
		ret.push_back(tmp_link);
	}
	fprintf(stderr, "\r calculation ... %d/%d(%4.1lf%%)\n", count, link.size(), count * 100. / link.size());

	return ret;

}
void position_difference(netscan::linklet_t link, double& dr, double& dl) {
	using namespace matrix_3D;
	vector_3D pos0, pos1, dir0, dir1;
	pos0.x = link.b[0].x;
	pos0.y = link.b[0].y;
	pos0.z = link.b[0].z;
	dir0.x = link.b[0].ax;
	dir0.y = link.b[0].ay;
	dir0.z = 1;
	pos1.x = link.b[1].x;
	pos1.y = link.b[1].y;
	pos1.z = link.b[1].z;
	dir1.x = link.b[1].ax;
	dir1.y = link.b[1].ay;
	dir1.z = 1;

	vector_3D base_point, difference;
	//äOë}äÓèÄì_Ç1:1Ç…ì‡ï™ÇµÇΩì_Ç…ê›íË
	base_point = addition(const_multiple(pos0, 0.5), const_multiple(pos1, 0.5));
	difference = addition(const_multiple(pos0, -1), pos1);

	vector_3D extra0, extra1;
	double ratio0, ratio1;
	ratio0 = -1 * dot(addition(pos0, const_multiple(base_point, -1)), difference) / dot(dir0, difference);
	ratio1 = -1 * dot(addition(pos1, const_multiple(base_point, -1)), difference) / dot(dir1, difference);
	extra0 = addition(pos0, const_multiple(dir0, ratio0));
	extra1 = addition(pos1, const_multiple(dir1, ratio1));

	vector_3D unit_r, unit_l;
	unit_l.x = -1 * difference.y;// + --> * Ç…èCê≥
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

	dr = dot(addition(extra1, const_multiple(extra0, -1)), unit_r);
	dl = dot(addition(extra1, const_multiple(extra0, -1)), unit_l);

}
void Calc_MD_OA(netscan::linklet_t link, double& md, double& oa) {
	using namespace matrix_3D;
	vector_3D pos0, pos1, dir0, dir1;
	pos0.x = link.b[0].x;
	pos0.y = link.b[0].y;
	pos0.z = link.b[0].z;
	dir0.x = link.b[0].ax;
	dir0.y = link.b[0].ay;
	dir0.z = 1;
	pos1.x = link.b[1].x;
	pos1.y = link.b[1].y;
	pos1.z = link.b[1].z;
	dir1.x = link.b[1].ax;
	dir1.y = link.b[1].ay;
	dir1.z = 1;

	oa = opening_angle(dir0, dir1);
	double z_range[2], extra[2];
	z_range[0] = pos0.z;
	z_range[1] = pos1.z;
	md = minimum_distance(pos0, pos1, dir0, dir1, z_range, extra);

}
void write_linklet_txt(std::string filename, std::vector<Linklet_rl> l) {
	std::ofstream ofs(filename);
	if (!ofs) {
		//file open é∏îs
		fprintf(stderr, "File[%s] is not exist!!\n", filename.c_str());
		exit(1);
	}
	if (l.size() == 0) {
		fprintf(stderr, "target linklet ... null\n");
		fprintf(stderr, "File[%s] has no text\n", filename.c_str());
	}
	else {
		int count = 0;
		std::cout << std::right << std::fixed;
		for (auto itr = l.begin(); itr != l.end(); itr++) {
			if (count % 10000 == 0) {
				fprintf(stderr, "\r Write Linklet ... %d/%d (%4.1lf%%)", count, int(l.size()), count * 100. / l.size());
			}
			count++;
			ofs << std::right << std::fixed
				<< std::setw(6) << itr->link.pos[0] << " "
				<< std::setw(15) << itr->link.b[0].rawid << " "
				<< std::setw(6) << itr->link.pos[1] << " "
				<< std::setw(15) << itr->link.b[1].rawid
				<< std::setw(8) << itr->link.b[0].m[0].ph + itr->link.b[0].m[1].ph << " "
				<< std::setw(8) << std::setprecision(4) << itr->link.b[0].ax << " "
				<< std::setw(8) << std::setprecision(4) << itr->link.b[0].ay << " "
				<< std::setw(10) << std::setprecision(1) << itr->link.b[0].x << " "
				<< std::setw(10) << std::setprecision(1) << itr->link.b[0].y << " "
				<< std::setw(8) << itr->link.b[1].m[0].ph + itr->link.b[1].m[1].ph << " "
				<< std::setw(8) << std::setprecision(4) << itr->link.b[1].ax << " "
				<< std::setw(8) << std::setprecision(4) << itr->link.b[1].ay << " "
				<< std::setw(10) << std::setprecision(1) << itr->link.b[1].x << " "
				<< std::setw(10) << std::setprecision(1) << itr->link.b[1].y << " "
				<< std::setw(11) << std::setprecision(1) << itr->link.b[0].z << " "
				<< std::setw(11) << std::setprecision(1) << itr->link.b[1].z << " "
				<< std::setw(11) << std::setprecision(1) << itr->link.zproj << " "
				<< std::setw(10) << std::setprecision(1) << itr->link.xc << " "
				<< std::setw(10) << std::setprecision(1) << itr->link.yc << " "
				<< std::setw(6) << itr->link.b[0].m[0].pos << " "
				<< std::setw(8) << itr->link.b[0].m[0].col << " "
				<< std::setw(10) << itr->link.b[0].m[0].row << " "
				<< std::setw(3) << itr->link.b[0].m[0].zone << " "
				<< std::setw(7) << itr->link.b[0].m[0].isg << " "
				<< std::setw(6) << itr->link.b[0].m[1].pos << " "
				<< std::setw(8) << itr->link.b[0].m[1].col << " "
				<< std::setw(10) << itr->link.b[0].m[1].row << " "
				<< std::setw(3) << itr->link.b[0].m[1].zone << " "
				<< std::setw(7) << itr->link.b[0].m[1].isg << " "
				<< std::setw(6) << itr->link.b[1].m[0].pos << " "
				<< std::setw(8) << itr->link.b[1].m[0].col << " "
				<< std::setw(10) << itr->link.b[1].m[0].row << " "
				<< std::setw(3) << itr->link.b[1].m[0].zone << " "
				<< std::setw(7) << itr->link.b[1].m[0].isg << " "
				<< std::setw(6) << itr->link.b[1].m[1].pos << " "
				<< std::setw(8) << itr->link.b[1].m[1].col << " "
				<< std::setw(10) << itr->link.b[1].m[1].row << " "
				<< std::setw(3) << itr->link.b[1].m[1].zone << " "
				<< std::setw(7) << itr->link.b[1].m[1].isg << " "
				<< "0" << " " << "0.0" << " "
				<< std::setw(15) << itr->link.b[0].m[0].rawid << " "
				<< std::setw(8) << itr->link.b[0].m[0].ph << " "
				<< std::setw(8) << std::setprecision(4) << itr->link.b[0].m[0].ax << " "
				<< std::setw(8) << std::setprecision(4) << itr->link.b[0].m[0].ay << " "
				<< std::setw(15) << itr->link.b[0].m[1].rawid << " "
				<< std::setw(8) << itr->link.b[0].m[1].ph << " "
				<< std::setw(8) << std::setprecision(4) << itr->link.b[0].m[1].ax << " "
				<< std::setw(8) << std::setprecision(4) << itr->link.b[0].m[1].ay << " "
				<< std::setw(15) << itr->link.b[1].m[0].rawid << " "
				<< std::setw(8) << itr->link.b[1].m[0].ph << " "
				<< std::setw(8) << std::setprecision(4) << itr->link.b[1].m[0].ax << " "
				<< std::setw(8) << std::setprecision(4) << itr->link.b[1].m[0].ay << " "
				<< std::setw(15) << itr->link.b[1].m[1].rawid << " "
				<< std::setw(8) << itr->link.b[1].m[1].ph << " "
				<< std::setw(8) << std::setprecision(4) << itr->link.b[1].m[1].ax << " "
				<< std::setw(8) << std::setprecision(4) << itr->link.b[1].m[1].ay << " "
				<< std::setw(10) << std::setprecision(1) << itr->link.dx << " "
				<< std::setw(10) << std::setprecision(1) << itr->link.dy << " "
				<< std::setw(10) << std::setprecision(1) << itr->dr << " "
				<< std::setw(10) << std::setprecision(1) << itr->dl << std::endl;
		}
		fprintf(stderr, "\r Write Linklet ... %d/%d (%4.1lf%%)\n", count, int(l.size()), count * 100. / l.size());
	}
	ofs.close();

}
void write_linklet_bin(std::string filename, std::vector<Linklet_rl> l) {
	std::ofstream ofs(filename, std::ios::binary);
	if (!ofs) {
		//file open é∏îs
		fprintf(stderr, "File[%s] is not exist!!\n", filename.c_str());
		exit(1);
	}
	if (l.size() == 0) {
		fprintf(stderr, "target linklet ... null\n");
		fprintf(stderr, "File[%s] has no text\n", filename.c_str());
	}
	int64_t count = 0;
	int64_t max = l.size();
	printf("Linklet Radial Lateral size = %d\n", sizeof(Linklet_rl));
	printf("calss description\n");
	printf(":::::::::::::::::::::::::::::::::::::::::::::::::\n");
	printf("class Linklet_rl {\n");
	printf(" public:\n");
	printf("  linklet_t link;\n");
	printf("  double dr, dl, md, oa;\n");
	printf("};\n");
	printf(":::::::::::::::::::::::::::::::::::::::::::::::::\n");

	for (int i = 0; i < l.size(); i++) {
		if (count % 10000 == 0) {
			std::cerr << std::right << std::fixed << "\r now writing ..." << std::setw(4) << std::setprecision(1) << count * 100. / max << "%";
		}
		count++;
		ofs.write((char*)&l[i], sizeof(Linklet_rl));
	}
	std::cerr << std::right << std::fixed << "\r now writing ..." << std::setw(4) << std::setprecision(1) << count * 100. / max << "%" << std::endl;
	ofs.close();
}
