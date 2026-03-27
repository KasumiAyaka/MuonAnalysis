// 2024/09/06
// kasumi
// based on "Check_upstream_base.cpp

#pragma comment(lib,"FILE_structure.lib")
#include <FILE_structure.hpp>
#include <iomanip>

class output_format {
public:
	int groupid, chainid, pl,upl;
	int peke, count;
	double dal,dar, md, dz,dz2,oa;
};

struct kink_cand {
	int pl, groupid;
};
bool operator<(const kink_cand& lhs, const kink_cand& rhs) {
	return std::tie(lhs.groupid, lhs.pl) < std::tie(rhs.groupid, rhs.pl);
}


void evaluate(std::vector<Momentum_recon::Event_information>& mom, std::vector<output_format>& out, int gid, int peke);
void evaluate_all(std::vector<Momentum_recon::Event_information>& mom, std::vector<output_format>& out, int peke);
void output(std::string output, std::vector<output_format>& out);


int main(int argc, char** argv) {
	if (argc < 4 || argc>5) {
		fprintf(stderr, "usage:file-in.momch file-out.txt groupid [ uplimit peke ]\n");
		exit(1);
	}
	std::string file_in_momch = argv[1];
	std::string file_out_txt = argv[2];
	int gid = std::stoi(argv[3]);
	int peke = 7;
	int mode = -1;
	if (argc == 5) {
		peke = std::stoi(argv[4]);
	}


	std::vector<Momentum_recon::Event_information> momch = Momentum_recon::Read_Event_information_extension(file_in_momch);
	std::vector<output_format>out;

	if (gid >0) {
		evaluate(momch, out, gid, peke);
	}
	else {
		evaluate_all(momch, out, peke);
	}
	output(file_out_txt, out);

	//Momentum_recon::Write_Event_information_extension(file_cut_momch, momch);
}
void evaluate(std::vector<Momentum_recon::Event_information>& mom, std::vector<output_format>& out, int gid, int peke) {
	std::cout << "Now calcurating ..." << std::endl;
	int count = 0;
	std::map<int, Momentum_recon::Mom_basetrack> map;
	output_format out_tmp;
	double ax, ay, dax, day, dlat, angle, drad;
	double extra[2], z_range[2];
	matrix_3D::vector_3D pos0, pos1, dir0, dir1;
	std::set<int> pl_set;

	for (auto& ev : mom) {
		if (ev.groupid == gid) {
			std::cout << " EventID = " << ev.groupid << std::endl;
			// only process indicated event
			for (auto& c : ev.chains) {
				out_tmp.groupid = gid;
				out_tmp.chainid = c.chainid;

				for (auto itr = c.base.begin(); itr != c.base.end(); itr++) {
					//std::cout << itr->pl << std::endl;

					pl_set.insert(itr->pl);
					map.insert(std::make_pair(itr->pl, *itr));
				}

				//std::vector<int> pl_v(pl_set.begin(), pl_set.end());
				//std::sort(pl_v.begin(), pl_v.end(), std::greater<int>());//became smaller

			}
		}

	}
	for (int pk = 0; pk < peke; pk++) {
		out_tmp.peke = pk;
		count = 0;

		//for (int i = 0; i < pl_v.size(); i++) {
		for (auto p = pl_set.begin(); p != pl_set.end();p++) {

			//auto itr1 = map.find(pl_v[i]);
			//auto itr2 = map.find(pl_v[i] + pk + 1);
			auto itr1 = map.find(*p);
			auto itr2 = map.find(*p + pk + 1);

			if (itr1 != map.end() || itr2 != map.end()) {

				out_tmp.pl = itr1->first;
				out_tmp.upl = itr2->first;

				if (out_tmp.upl - out_tmp.pl - 1 == pk) {
					std::cout << out_tmp.peke << " " << out_tmp.pl << " " << out_tmp.upl << "\n";

					// MD
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

					z_range[0] = pos0.z;
					z_range[1] = pos1.z;

					double point[3];
					out_tmp.md = matrix_3D::minimum_distance(pos0, pos1, dir0, dir1, point);
					out_tmp.oa = matrix_3D::opening_angle(dir0, dir1);
					// MD no tansaku hani wo Z-range de sitei
					//out_tmp.md = matrix_3D::minimum_distance(pos0, pos1, dir0, dir1, z_range, extra);
					//out_tmp.dz = (fabs(extra[0]) + fabs(extra[1])) / 2;
					out_tmp.dz = point[2];

					// d(lat)
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



					out_tmp.dal = dlat;
					out_tmp.dar = drad;
					out_tmp.count = count;

					out_tmp.dz2 = itr2->second.z - itr1->second.z;
					out.push_back(out_tmp);
					count++;
				//	std::cout << out_tmp.peke << " " << out_tmp.pl << " " << out_tmp.upl << std::endl;

				}
			}
		}
	}


}
void evaluate_all(std::vector<Momentum_recon::Event_information>& mom, std::vector<output_format>& out, int peke) {
	std::cout << "Now calcurating ..." << std::endl;
	int count = 0;
	output_format out_tmp;
	double ax, ay, dax, day, dlat, angle, drad;
	double extra[2], z_range[2];
	matrix_3D::vector_3D pos0, pos1, dir0, dir1;

	for (auto& ev : mom) {

		std::cout << " EventID = " << ev.groupid << std::endl;
		// only process indicated event
		for (auto& c : ev.chains) {
			out_tmp.groupid = ev.groupid;
			out_tmp.chainid = c.chainid;

			std::map<int, Momentum_recon::Mom_basetrack> map;
			std::set<int> pl_set;

			for (auto itr = c.base.begin(); itr != c.base.end(); itr++) {
				//std::cout << itr->pl << std::endl;

				pl_set.insert(itr->pl);
				map.insert(std::make_pair(itr->pl, *itr));
			}

			for (int pk = 0; pk < peke; pk++) {
				out_tmp.peke = pk;
				count = 0;

				//for (int i = 0; i < pl_v.size(); i++) {
				for (auto p = pl_set.begin(); p != pl_set.end(); p++) {

					//auto itr1 = map.find(pl_v[i]);
					//auto itr2 = map.find(pl_v[i] + pk + 1);
					auto itr1 = map.find(*p);
					auto itr2 = map.find(*p + pk + 1);

					if (itr1 != map.end() || itr2 != map.end()) {

						out_tmp.pl = itr1->first;
						out_tmp.upl = itr2->first;

						if (out_tmp.upl - out_tmp.pl - 1 == pk) {

							// MD
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

							z_range[0] = pos0.z;
							z_range[1] = pos1.z;

							out_tmp.oa = matrix_3D::opening_angle(dir0, dir1);
							// MD no tansaku hani wo Z-range de sitei
							//out_tmp.md = matrix_3D::minimum_distance(pos0, pos1, dir0, dir1, z_range, extra);
							//out_tmp.dz = (fabs(extra[0]) + fabs(extra[1])) / 2;
							double point[3];
							out_tmp.md = matrix_3D::minimum_distance(pos0, pos1, dir0, dir1, point);
							out_tmp.dz = point[2];

							// d(lat)
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



							out_tmp.dal = dlat;
							out_tmp.dar = drad;
							out_tmp.count = count;

							out_tmp.dz2 = itr2->second.z - itr1->second.z;
							out.push_back(out_tmp);
							count++;
							//	std::cout << out_tmp.peke << " " << out_tmp.pl << " " << out_tmp.upl << std::endl;

						}
					}
				}
			}

		}

	}



}


void output(std::string output, std::vector<output_format>& out) {

	std::ofstream ofs(output);
	if (!ofs) {
		std::cerr << "File open error : " << output << std::endl;
		exit(1);
	}
	int count = 0;

	for (auto itr = out.begin(); itr < out.end(); itr++) {
		if (itr->upl < 17) {
			count = 0;
			continue;

		}
		ofs<<std::fixed<<std::right
			<< std::setw(6) << std::setprecision(0) << itr->groupid << " "
			<< std::setw(3) << std::setprecision(0) << count << " "
			<< std::setw(3) << std::setprecision(0) <<  itr->peke<< " "
			<< std::setw(3) << std::setprecision(0) << itr->pl << " "
			<< std::setw(3) << std::setprecision(0) << itr->upl << " "
			<< std::setw(10) << std::setprecision(1) << itr->md << " "
			<< std::setw(10) << std::setprecision(5) << itr->oa << " "
			<< std::setw(10) << std::setprecision(1) << itr->dz << " "
			<< std::setw(10) << std::setprecision(5) <<  itr->dal<< " "
			<< std::setw(10) << std::setprecision(5) <<  itr->dar<< " "
			<< std::setw(10) << std::setprecision(1) << itr->dz2 << " "
			<< std::endl;
		count++;
	}

	std::cout << "writing file is finished!" << std::endl;
}